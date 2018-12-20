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
         mat(k,533) = -(rxt(k,356)*y(k,212))
         mat(k,1562) = -rxt(k,356)*y(k,1)
         mat(k,1942) = rxt(k,359)*y(k,189)
         mat(k,806) = rxt(k,359)*y(k,124)
         mat(k,522) = -(rxt(k,360)*y(k,212))
         mat(k,1561) = -rxt(k,360)*y(k,2)
         mat(k,805) = rxt(k,357)*y(k,201)
         mat(k,1344) = rxt(k,357)*y(k,189)
         mat(k,761) = -(rxt(k,439)*y(k,126) + rxt(k,440)*y(k,134) + rxt(k,441) &
                      *y(k,212))
         mat(k,1415) = -rxt(k,439)*y(k,6)
         mat(k,1800) = -rxt(k,440)*y(k,6)
         mat(k,1584) = -rxt(k,441)*y(k,6)
         mat(k,81) = -(rxt(k,398)*y(k,212))
         mat(k,1497) = -rxt(k,398)*y(k,7)
         mat(k,284) = -(rxt(k,401)*y(k,212))
         mat(k,1529) = -rxt(k,401)*y(k,8)
         mat(k,374) = rxt(k,399)*y(k,201)
         mat(k,1320) = rxt(k,399)*y(k,190)
         mat(k,82) = .120_r8*rxt(k,398)*y(k,212)
         mat(k,1498) = .120_r8*rxt(k,398)*y(k,7)
         mat(k,758) = .100_r8*rxt(k,440)*y(k,134)
         mat(k,784) = .100_r8*rxt(k,443)*y(k,134)
         mat(k,1789) = .100_r8*rxt(k,440)*y(k,6) + .100_r8*rxt(k,443)*y(k,110)
         mat(k,1930) = .500_r8*rxt(k,400)*y(k,190) + .200_r8*rxt(k,427)*y(k,218) &
                      + .060_r8*rxt(k,433)*y(k,220)
         mat(k,375) = .500_r8*rxt(k,400)*y(k,124)
         mat(k,590) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,606) = .060_r8*rxt(k,433)*y(k,124)
         mat(k,1923) = .200_r8*rxt(k,427)*y(k,218) + .200_r8*rxt(k,433)*y(k,220)
         mat(k,589) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,604) = .200_r8*rxt(k,433)*y(k,124)
         mat(k,1939) = .200_r8*rxt(k,427)*y(k,218) + .150_r8*rxt(k,433)*y(k,220)
         mat(k,592) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,607) = .150_r8*rxt(k,433)*y(k,124)
         mat(k,1925) = .210_r8*rxt(k,433)*y(k,220)
         mat(k,605) = .210_r8*rxt(k,433)*y(k,124)
         mat(k,158) = -(rxt(k,361)*y(k,212))
         mat(k,1510) = -rxt(k,361)*y(k,15)
         mat(k,757) = .050_r8*rxt(k,440)*y(k,134)
         mat(k,783) = .050_r8*rxt(k,443)*y(k,134)
         mat(k,1788) = .050_r8*rxt(k,440)*y(k,6) + .050_r8*rxt(k,443)*y(k,110)
         mat(k,258) = -(rxt(k,327)*y(k,126) + rxt(k,328)*y(k,212))
         mat(k,1409) = -rxt(k,327)*y(k,16)
         mat(k,1525) = -rxt(k,328)*y(k,16)
         mat(k,1250) = -(rxt(k,210)*y(k,42) + rxt(k,211)*y(k,201) + rxt(k,212) &
                      *y(k,134))
         mat(k,1902) = -rxt(k,210)*y(k,17)
         mat(k,1386) = -rxt(k,211)*y(k,17)
         mat(k,1825) = -rxt(k,212)*y(k,17)
         mat(k,1849) = 4.000_r8*rxt(k,213)*y(k,19) + (rxt(k,214)+rxt(k,215))*y(k,59) &
                      + rxt(k,218)*y(k,124) + rxt(k,221)*y(k,133) + rxt(k,466) &
                      *y(k,150) + rxt(k,222)*y(k,212)
         mat(k,1690) = (rxt(k,214)+rxt(k,215))*y(k,19)
         mat(k,696) = rxt(k,223)*y(k,133) + rxt(k,229)*y(k,211) + rxt(k,224)*y(k,212)
         mat(k,1980) = rxt(k,218)*y(k,19)
         mat(k,1879) = rxt(k,221)*y(k,19) + rxt(k,223)*y(k,81)
         mat(k,1084) = rxt(k,466)*y(k,19)
         mat(k,1467) = rxt(k,229)*y(k,81)
         mat(k,1614) = rxt(k,222)*y(k,19) + rxt(k,224)*y(k,81)
         mat(k,1843) = rxt(k,216)*y(k,59)
         mat(k,1684) = rxt(k,216)*y(k,19)
         mat(k,1289) = (rxt(k,516)+rxt(k,521))*y(k,91)
         mat(k,639) = (rxt(k,516)+rxt(k,521))*y(k,85)
         mat(k,1862) = -(4._r8*rxt(k,213)*y(k,19) + (rxt(k,214) + rxt(k,215) + rxt(k,216) &
                      ) * y(k,59) + rxt(k,217)*y(k,201) + rxt(k,218)*y(k,124) &
                      + rxt(k,219)*y(k,125) + rxt(k,221)*y(k,133) + rxt(k,222) &
                      *y(k,212) + rxt(k,466)*y(k,150))
         mat(k,1703) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,19)
         mat(k,1399) = -rxt(k,217)*y(k,19)
         mat(k,1993) = -rxt(k,218)*y(k,19)
         mat(k,1778) = -rxt(k,219)*y(k,19)
         mat(k,1892) = -rxt(k,221)*y(k,19)
         mat(k,1627) = -rxt(k,222)*y(k,19)
         mat(k,1092) = -rxt(k,466)*y(k,19)
         mat(k,1256) = rxt(k,212)*y(k,134)
         mat(k,449) = rxt(k,220)*y(k,133)
         mat(k,700) = rxt(k,230)*y(k,211)
         mat(k,645) = rxt(k,225)*y(k,133)
         mat(k,1892) = mat(k,1892) + rxt(k,220)*y(k,20) + rxt(k,225)*y(k,91)
         mat(k,1838) = rxt(k,212)*y(k,17)
         mat(k,1480) = rxt(k,230)*y(k,81)
         mat(k,443) = -(rxt(k,220)*y(k,133))
         mat(k,1869) = -rxt(k,220)*y(k,20)
         mat(k,1845) = rxt(k,219)*y(k,125)
         mat(k,1750) = rxt(k,219)*y(k,19)
         mat(k,164) = -(rxt(k,402)*y(k,212))
         mat(k,1511) = -rxt(k,402)*y(k,22)
         mat(k,1921) = rxt(k,405)*y(k,191)
         mat(k,332) = rxt(k,405)*y(k,124)
         mat(k,240) = -(rxt(k,404)*y(k,212))
         mat(k,1522) = -rxt(k,404)*y(k,23)
         mat(k,333) = rxt(k,403)*y(k,201)
         mat(k,1317) = rxt(k,403)*y(k,191)
         mat(k,199) = -(rxt(k,276)*y(k,56) + rxt(k,277)*y(k,212))
         mat(k,1709) = -rxt(k,276)*y(k,24)
         mat(k,1516) = -rxt(k,277)*y(k,24)
         mat(k,451) = -(rxt(k,278)*y(k,56) + rxt(k,279)*y(k,134) + rxt(k,304)*y(k,212))
         mat(k,1711) = -rxt(k,278)*y(k,25)
         mat(k,1793) = -rxt(k,279)*y(k,25)
         mat(k,1552) = -rxt(k,304)*y(k,25)
         mat(k,172) = -(rxt(k,284)*y(k,212))
         mat(k,1513) = -rxt(k,284)*y(k,26)
         mat(k,684) = .800_r8*rxt(k,280)*y(k,192) + .200_r8*rxt(k,281)*y(k,196)
         mat(k,1632) = .200_r8*rxt(k,281)*y(k,192)
         mat(k,245) = -(rxt(k,285)*y(k,212))
         mat(k,1523) = -rxt(k,285)*y(k,27)
         mat(k,685) = rxt(k,282)*y(k,201)
         mat(k,1318) = rxt(k,282)*y(k,192)
         mat(k,205) = -(rxt(k,286)*y(k,56) + rxt(k,287)*y(k,212))
         mat(k,1710) = -rxt(k,286)*y(k,28)
         mat(k,1517) = -rxt(k,287)*y(k,28)
         mat(k,841) = -(rxt(k,307)*y(k,126) + rxt(k,308)*y(k,134) + rxt(k,325) &
                      *y(k,212))
         mat(k,1419) = -rxt(k,307)*y(k,29)
         mat(k,1804) = -rxt(k,308)*y(k,29)
         mat(k,1589) = -rxt(k,325)*y(k,29)
         mat(k,710) = .130_r8*rxt(k,385)*y(k,134)
         mat(k,1804) = mat(k,1804) + .130_r8*rxt(k,385)*y(k,98)
         mat(k,314) = -(rxt(k,312)*y(k,212))
         mat(k,1533) = -rxt(k,312)*y(k,30)
         mat(k,665) = rxt(k,310)*y(k,201)
         mat(k,1324) = rxt(k,310)*y(k,193)
         mat(k,55) = -(rxt(k,313)*y(k,212))
         mat(k,1494) = -rxt(k,313)*y(k,31)
         mat(k,176) = -(rxt(k,408)*y(k,212))
         mat(k,1514) = -rxt(k,408)*y(k,32)
         mat(k,513) = rxt(k,406)*y(k,201)
         mat(k,1312) = rxt(k,406)*y(k,194)
         mat(k,1917) = -(rxt(k,174)*y(k,56) + rxt(k,210)*y(k,17) + rxt(k,254)*y(k,201) &
                      + rxt(k,255)*y(k,126) + rxt(k,256)*y(k,133) + rxt(k,257) &
                      *y(k,212))
         mat(k,1739) = -rxt(k,174)*y(k,42)
         mat(k,1258) = -rxt(k,210)*y(k,42)
         mat(k,1401) = -rxt(k,254)*y(k,42)
         mat(k,1458) = -rxt(k,255)*y(k,42)
         mat(k,1894) = -rxt(k,256)*y(k,42)
         mat(k,1629) = -rxt(k,257)*y(k,42)
         mat(k,542) = .400_r8*rxt(k,356)*y(k,212)
         mat(k,776) = .340_r8*rxt(k,440)*y(k,134)
         mat(k,265) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,458) = rxt(k,279)*y(k,134)
         mat(k,853) = .500_r8*rxt(k,308)*y(k,134)
         mat(k,395) = .500_r8*rxt(k,296)*y(k,212)
         mat(k,664) = rxt(k,262)*y(k,212)
         mat(k,324) = .300_r8*rxt(k,263)*y(k,212)
         mat(k,1705) = rxt(k,181)*y(k,196)
         mat(k,869) = .800_r8*rxt(k,301)*y(k,212)
         mat(k,723) = .910_r8*rxt(k,385)*y(k,134)
         mat(k,475) = .300_r8*rxt(k,376)*y(k,212)
         mat(k,1056) = .800_r8*rxt(k,380)*y(k,196)
         mat(k,1069) = .120_r8*rxt(k,338)*y(k,134)
         mat(k,422) = .500_r8*rxt(k,351)*y(k,212)
         mat(k,802) = .340_r8*rxt(k,443)*y(k,134)
         mat(k,1141) = .600_r8*rxt(k,352)*y(k,134)
         mat(k,1995) = .100_r8*rxt(k,358)*y(k,189) + rxt(k,261)*y(k,196) &
                      + .500_r8*rxt(k,329)*y(k,198) + .500_r8*rxt(k,298)*y(k,200) &
                      + .920_r8*rxt(k,368)*y(k,203) + .250_r8*rxt(k,336)*y(k,205) &
                      + rxt(k,345)*y(k,207) + rxt(k,319)*y(k,214) + rxt(k,323) &
                      *y(k,215) + .340_r8*rxt(k,452)*y(k,216) + .320_r8*rxt(k,457) &
                      *y(k,217) + .250_r8*rxt(k,393)*y(k,219)
         mat(k,1458) = mat(k,1458) + .500_r8*rxt(k,327)*y(k,16) + rxt(k,369)*y(k,203) &
                      + .250_r8*rxt(k,335)*y(k,205) + rxt(k,346)*y(k,207)
         mat(k,1840) = .340_r8*rxt(k,440)*y(k,6) + rxt(k,279)*y(k,25) &
                      + .500_r8*rxt(k,308)*y(k,29) + .910_r8*rxt(k,385)*y(k,98) &
                      + .120_r8*rxt(k,338)*y(k,105) + .340_r8*rxt(k,443)*y(k,110) &
                      + .600_r8*rxt(k,352)*y(k,111)
         mat(k,362) = rxt(k,303)*y(k,212)
         mat(k,893) = .680_r8*rxt(k,461)*y(k,212)
         mat(k,817) = .100_r8*rxt(k,358)*y(k,124)
         mat(k,693) = .700_r8*rxt(k,281)*y(k,196)
         mat(k,673) = rxt(k,309)*y(k,196)
         mat(k,1245) = rxt(k,292)*y(k,196) + rxt(k,365)*y(k,203) + .250_r8*rxt(k,332) &
                      *y(k,205) + rxt(k,341)*y(k,207) + .250_r8*rxt(k,390)*y(k,219)
         mat(k,1679) = rxt(k,181)*y(k,59) + .800_r8*rxt(k,380)*y(k,101) + rxt(k,261) &
                      *y(k,124) + .700_r8*rxt(k,281)*y(k,192) + rxt(k,309)*y(k,193) &
                      + rxt(k,292)*y(k,195) + (4.000_r8*rxt(k,258)+2.000_r8*rxt(k,259)) &
                      *y(k,196) + 1.500_r8*rxt(k,366)*y(k,203) + .750_r8*rxt(k,371) &
                      *y(k,204) + .880_r8*rxt(k,333)*y(k,205) + 2.000_r8*rxt(k,342) &
                      *y(k,207) + .750_r8*rxt(k,445)*y(k,210) + .800_r8*rxt(k,321) &
                      *y(k,215) + .930_r8*rxt(k,450)*y(k,216) + .950_r8*rxt(k,455) &
                      *y(k,217) + .800_r8*rxt(k,391)*y(k,219)
         mat(k,465) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,581) = .500_r8*rxt(k,298)*y(k,124)
         mat(k,1401) = mat(k,1401) + .450_r8*rxt(k,343)*y(k,207) + .150_r8*rxt(k,322) &
                      *y(k,215)
         mat(k,1121) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) + rxt(k,365) &
                      *y(k,195) + 1.500_r8*rxt(k,366)*y(k,196)
         mat(k,1197) = .750_r8*rxt(k,371)*y(k,196)
         mat(k,1163) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,195) + .880_r8*rxt(k,333)*y(k,196)
         mat(k,1215) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,341)*y(k,195) &
                      + 2.000_r8*rxt(k,342)*y(k,196) + .450_r8*rxt(k,343)*y(k,201) &
                      + 4.000_r8*rxt(k,344)*y(k,207)
         mat(k,979) = .750_r8*rxt(k,445)*y(k,196)
         mat(k,1629) = mat(k,1629) + .400_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,296) &
                      *y(k,51) + rxt(k,262)*y(k,52) + .300_r8*rxt(k,263)*y(k,53) &
                      + .800_r8*rxt(k,301)*y(k,74) + .300_r8*rxt(k,376)*y(k,99) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,139) &
                      + .680_r8*rxt(k,461)*y(k,178)
         mat(k,636) = rxt(k,319)*y(k,124)
         mat(k,992) = rxt(k,323)*y(k,124) + .800_r8*rxt(k,321)*y(k,196) &
                      + .150_r8*rxt(k,322)*y(k,201)
         mat(k,959) = .340_r8*rxt(k,452)*y(k,124) + .930_r8*rxt(k,450)*y(k,196)
         mat(k,940) = .320_r8*rxt(k,457)*y(k,124) + .950_r8*rxt(k,455)*y(k,196)
         mat(k,1033) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,390)*y(k,195) &
                      + .800_r8*rxt(k,391)*y(k,196)
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
         mat(k,996) = -(rxt(k,288)*y(k,126) + rxt(k,289)*y(k,212))
         mat(k,1431) = -rxt(k,288)*y(k,45)
         mat(k,1601) = -rxt(k,289)*y(k,45)
         mat(k,537) = .800_r8*rxt(k,356)*y(k,212)
         mat(k,261) = rxt(k,327)*y(k,126)
         mat(k,173) = rxt(k,284)*y(k,212)
         mat(k,247) = .500_r8*rxt(k,285)*y(k,212)
         mat(k,844) = .500_r8*rxt(k,308)*y(k,134)
         mat(k,1128) = .100_r8*rxt(k,352)*y(k,134)
         mat(k,1969) = .400_r8*rxt(k,358)*y(k,189) + rxt(k,283)*y(k,192) &
                      + .270_r8*rxt(k,311)*y(k,193) + rxt(k,329)*y(k,198) + rxt(k,348) &
                      *y(k,209) + rxt(k,319)*y(k,214)
         mat(k,1431) = mat(k,1431) + rxt(k,327)*y(k,16)
         mat(k,1814) = .500_r8*rxt(k,308)*y(k,29) + .100_r8*rxt(k,352)*y(k,111)
         mat(k,811) = .400_r8*rxt(k,358)*y(k,124)
         mat(k,688) = rxt(k,283)*y(k,124) + 3.200_r8*rxt(k,280)*y(k,192) &
                      + .800_r8*rxt(k,281)*y(k,196)
         mat(k,668) = .270_r8*rxt(k,311)*y(k,124)
         mat(k,1654) = .800_r8*rxt(k,281)*y(k,192)
         mat(k,462) = rxt(k,329)*y(k,124)
         mat(k,1374) = .200_r8*rxt(k,347)*y(k,209)
         mat(k,545) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,201)
         mat(k,1601) = mat(k,1601) + .800_r8*rxt(k,356)*y(k,1) + rxt(k,284)*y(k,26) &
                      + .500_r8*rxt(k,285)*y(k,27)
         mat(k,632) = rxt(k,319)*y(k,124)
         mat(k,52) = -(rxt(k,290)*y(k,212))
         mat(k,1493) = -rxt(k,290)*y(k,47)
         mat(k,819) = -(rxt(k,326)*y(k,212))
         mat(k,1587) = -rxt(k,326)*y(k,48)
         mat(k,536) = .800_r8*rxt(k,356)*y(k,212)
         mat(k,763) = .520_r8*rxt(k,440)*y(k,134)
         mat(k,260) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,789) = .520_r8*rxt(k,443)*y(k,134)
         mat(k,1957) = .250_r8*rxt(k,358)*y(k,189) + .820_r8*rxt(k,311)*y(k,193) &
                      + .500_r8*rxt(k,329)*y(k,198) + .270_r8*rxt(k,452)*y(k,216) &
                      + .040_r8*rxt(k,457)*y(k,217)
         mat(k,1418) = .500_r8*rxt(k,327)*y(k,16)
         mat(k,1803) = .520_r8*rxt(k,440)*y(k,6) + .520_r8*rxt(k,443)*y(k,110)
         mat(k,885) = .500_r8*rxt(k,461)*y(k,212)
         mat(k,810) = .250_r8*rxt(k,358)*y(k,124)
         mat(k,667) = .820_r8*rxt(k,311)*y(k,124) + .820_r8*rxt(k,309)*y(k,196)
         mat(k,1643) = .820_r8*rxt(k,309)*y(k,193) + .150_r8*rxt(k,450)*y(k,216) &
                      + .025_r8*rxt(k,455)*y(k,217)
         mat(k,460) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,1587) = mat(k,1587) + .800_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,461) &
                      *y(k,178)
         mat(k,946) = .270_r8*rxt(k,452)*y(k,124) + .150_r8*rxt(k,450)*y(k,196)
         mat(k,924) = .040_r8*rxt(k,457)*y(k,124) + .025_r8*rxt(k,455)*y(k,196)
         mat(k,1072) = -(rxt(k,314)*y(k,126) + rxt(k,315)*y(k,212))
         mat(k,1435) = -rxt(k,314)*y(k,49)
         mat(k,1606) = -rxt(k,315)*y(k,49)
         mat(k,916) = rxt(k,316)*y(k,212)
         mat(k,1061) = .880_r8*rxt(k,338)*y(k,134)
         mat(k,1129) = .500_r8*rxt(k,352)*y(k,134)
         mat(k,1973) = .170_r8*rxt(k,411)*y(k,197) + .050_r8*rxt(k,374)*y(k,204) &
                      + .250_r8*rxt(k,336)*y(k,205) + .170_r8*rxt(k,417)*y(k,208) &
                      + .400_r8*rxt(k,427)*y(k,218) + .250_r8*rxt(k,393)*y(k,219) &
                      + .540_r8*rxt(k,433)*y(k,220) + .510_r8*rxt(k,436)*y(k,221)
         mat(k,1435) = mat(k,1435) + .050_r8*rxt(k,375)*y(k,204) + .250_r8*rxt(k,335) &
                      *y(k,205) + .250_r8*rxt(k,394)*y(k,219)
         mat(k,733) = rxt(k,317)*y(k,212)
         mat(k,1817) = .880_r8*rxt(k,338)*y(k,105) + .500_r8*rxt(k,352)*y(k,111)
         mat(k,1230) = .250_r8*rxt(k,332)*y(k,205) + .250_r8*rxt(k,390)*y(k,219)
         mat(k,1658) = .240_r8*rxt(k,333)*y(k,205) + .500_r8*rxt(k,321)*y(k,215) &
                      + .100_r8*rxt(k,391)*y(k,219)
         mat(k,623) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,201)
         mat(k,1379) = .070_r8*rxt(k,410)*y(k,197) + .070_r8*rxt(k,416)*y(k,208)
         mat(k,1183) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1152) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,195) + .240_r8*rxt(k,333)*y(k,196)
         mat(k,746) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,201)
         mat(k,1606) = mat(k,1606) + rxt(k,316)*y(k,95) + rxt(k,317)*y(k,127)
         mat(k,986) = .500_r8*rxt(k,321)*y(k,196)
         mat(k,599) = .400_r8*rxt(k,427)*y(k,124)
         mat(k,1025) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,195) + .100_r8*rxt(k,391)*y(k,196)
         mat(k,615) = .540_r8*rxt(k,433)*y(k,124)
         mat(k,386) = .510_r8*rxt(k,436)*y(k,124)
         mat(k,431) = -(rxt(k,295)*y(k,212))
         mat(k,1550) = -rxt(k,295)*y(k,50)
         mat(k,837) = .120_r8*rxt(k,308)*y(k,134)
         mat(k,1792) = .120_r8*rxt(k,308)*y(k,29)
         mat(k,1221) = .100_r8*rxt(k,292)*y(k,196) + .150_r8*rxt(k,293)*y(k,201)
         mat(k,1636) = .100_r8*rxt(k,292)*y(k,195)
         mat(k,1337) = .150_r8*rxt(k,293)*y(k,195) + .150_r8*rxt(k,343)*y(k,207)
         mat(k,1202) = .150_r8*rxt(k,343)*y(k,201)
         mat(k,391) = -(rxt(k,296)*y(k,212))
         mat(k,1544) = -rxt(k,296)*y(k,51)
         mat(k,1220) = .400_r8*rxt(k,293)*y(k,201)
         mat(k,1335) = .400_r8*rxt(k,293)*y(k,195) + .400_r8*rxt(k,343)*y(k,207)
         mat(k,1200) = .400_r8*rxt(k,343)*y(k,201)
         mat(k,661) = -(rxt(k,262)*y(k,212))
         mat(k,1574) = -rxt(k,262)*y(k,52)
         mat(k,1038) = .200_r8*rxt(k,380)*y(k,196)
         mat(k,686) = .300_r8*rxt(k,281)*y(k,196)
         mat(k,1638) = .200_r8*rxt(k,380)*y(k,101) + .300_r8*rxt(k,281)*y(k,192) &
                      + 2.000_r8*rxt(k,259)*y(k,196) + .250_r8*rxt(k,366)*y(k,203) &
                      + .250_r8*rxt(k,371)*y(k,204) + .250_r8*rxt(k,333)*y(k,205) &
                      + .250_r8*rxt(k,445)*y(k,210) + .500_r8*rxt(k,321)*y(k,215) &
                      + .250_r8*rxt(k,450)*y(k,216) + .250_r8*rxt(k,455)*y(k,217) &
                      + .300_r8*rxt(k,391)*y(k,219)
         mat(k,1098) = .250_r8*rxt(k,366)*y(k,196)
         mat(k,1171) = .250_r8*rxt(k,371)*y(k,196)
         mat(k,1145) = .250_r8*rxt(k,333)*y(k,196)
         mat(k,964) = .250_r8*rxt(k,445)*y(k,196)
         mat(k,983) = .500_r8*rxt(k,321)*y(k,196)
         mat(k,945) = .250_r8*rxt(k,450)*y(k,196)
         mat(k,923) = .250_r8*rxt(k,455)*y(k,196)
         mat(k,1019) = .300_r8*rxt(k,391)*y(k,196)
         mat(k,320) = -(rxt(k,263)*y(k,212))
         mat(k,1534) = -rxt(k,263)*y(k,53)
         mat(k,1635) = rxt(k,260)*y(k,201)
         mat(k,1325) = rxt(k,260)*y(k,196)
         mat(k,1734) = -(rxt(k,174)*y(k,42) + rxt(k,176)*y(k,77) + rxt(k,177)*y(k,79) &
                      + (rxt(k,178) + rxt(k,179)) * y(k,201) + rxt(k,180)*y(k,134) &
                      + rxt(k,187)*y(k,60) + rxt(k,196)*y(k,92) + rxt(k,286)*y(k,28))
         mat(k,1912) = -rxt(k,174)*y(k,56)
         mat(k,1015) = -rxt(k,176)*y(k,56)
         mat(k,480) = -rxt(k,177)*y(k,56)
         mat(k,1396) = -(rxt(k,178) + rxt(k,179)) * y(k,56)
         mat(k,1835) = -rxt(k,180)*y(k,56)
         mat(k,833) = -rxt(k,187)*y(k,56)
         mat(k,681) = -rxt(k,196)*y(k,56)
         mat(k,209) = -rxt(k,286)*y(k,56)
         mat(k,1859) = rxt(k,215)*y(k,59)
         mat(k,1700) = rxt(k,215)*y(k,19) + (4.000_r8*rxt(k,182)+2.000_r8*rxt(k,184)) &
                      *y(k,59) + rxt(k,186)*y(k,124) + rxt(k,191)*y(k,133) &
                      + rxt(k,467)*y(k,150) + rxt(k,181)*y(k,196) + rxt(k,192) &
                      *y(k,212)
         mat(k,103) = rxt(k,236)*y(k,211)
         mat(k,1303) = rxt(k,194)*y(k,133) + rxt(k,206)*y(k,211) + rxt(k,195)*y(k,212)
         mat(k,1990) = rxt(k,186)*y(k,59)
         mat(k,1889) = rxt(k,191)*y(k,59) + rxt(k,194)*y(k,85)
         mat(k,1089) = rxt(k,467)*y(k,59)
         mat(k,1674) = rxt(k,181)*y(k,59)
         mat(k,1477) = rxt(k,236)*y(k,65) + rxt(k,206)*y(k,85)
         mat(k,1624) = rxt(k,192)*y(k,59) + rxt(k,195)*y(k,85)
         mat(k,1708) = rxt(k,187)*y(k,60)
         mat(k,1683) = 2.000_r8*rxt(k,183)*y(k,59)
         mat(k,825) = rxt(k,187)*y(k,56) + (rxt(k,514)+rxt(k,519)+rxt(k,524))*y(k,85)
         mat(k,1288) = (rxt(k,514)+rxt(k,519)+rxt(k,524))*y(k,60) + (rxt(k,509) &
                       +rxt(k,515)+rxt(k,520))*y(k,92)
         mat(k,676) = (rxt(k,509)+rxt(k,515)+rxt(k,520))*y(k,85)
         mat(k,1682) = 2.000_r8*rxt(k,208)*y(k,59)
         mat(k,1699) = -(rxt(k,181)*y(k,196) + (4._r8*rxt(k,182) + 4._r8*rxt(k,183) &
                      + 4._r8*rxt(k,184) + 4._r8*rxt(k,208)) * y(k,59) + rxt(k,185) &
                      *y(k,201) + rxt(k,186)*y(k,124) + rxt(k,188)*y(k,125) + rxt(k,191) &
                      *y(k,133) + (rxt(k,192) + rxt(k,193)) * y(k,212) + (rxt(k,214) &
                      + rxt(k,215) + rxt(k,216)) * y(k,19) + rxt(k,467)*y(k,150))
         mat(k,1673) = -rxt(k,181)*y(k,59)
         mat(k,1395) = -rxt(k,185)*y(k,59)
         mat(k,1989) = -rxt(k,186)*y(k,59)
         mat(k,1774) = -rxt(k,188)*y(k,59)
         mat(k,1888) = -rxt(k,191)*y(k,59)
         mat(k,1623) = -(rxt(k,192) + rxt(k,193)) * y(k,59)
         mat(k,1858) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,59)
         mat(k,1088) = -rxt(k,467)*y(k,59)
         mat(k,1733) = rxt(k,196)*y(k,92) + rxt(k,180)*y(k,134) + rxt(k,179)*y(k,201)
         mat(k,832) = rxt(k,189)*y(k,133)
         mat(k,1302) = rxt(k,207)*y(k,211)
         mat(k,680) = rxt(k,196)*y(k,56) + rxt(k,197)*y(k,133) + rxt(k,198)*y(k,212)
         mat(k,1888) = mat(k,1888) + rxt(k,189)*y(k,60) + rxt(k,197)*y(k,92)
         mat(k,1834) = rxt(k,180)*y(k,56)
         mat(k,232) = rxt(k,472)*y(k,150)
         mat(k,1088) = mat(k,1088) + rxt(k,472)*y(k,136)
         mat(k,1395) = mat(k,1395) + rxt(k,179)*y(k,56)
         mat(k,1476) = rxt(k,207)*y(k,85)
         mat(k,1623) = mat(k,1623) + rxt(k,198)*y(k,92)
         mat(k,827) = -(rxt(k,187)*y(k,56) + rxt(k,189)*y(k,133) + rxt(k,190)*y(k,212) &
                      + (rxt(k,514) + rxt(k,519) + rxt(k,524)) * y(k,85))
         mat(k,1718) = -rxt(k,187)*y(k,60)
         mat(k,1875) = -rxt(k,189)*y(k,60)
         mat(k,1588) = -rxt(k,190)*y(k,60)
         mat(k,1292) = -(rxt(k,514) + rxt(k,519) + rxt(k,524)) * y(k,60)
         mat(k,1688) = rxt(k,188)*y(k,125)
         mat(k,1758) = rxt(k,188)*y(k,59)
         mat(k,911) = -((rxt(k,265) + rxt(k,275)) * y(k,212))
         mat(k,1595) = -(rxt(k,265) + rxt(k,275)) * y(k,62)
         mat(k,766) = .230_r8*rxt(k,440)*y(k,134)
         mat(k,1249) = rxt(k,210)*y(k,42)
         mat(k,202) = .350_r8*rxt(k,277)*y(k,212)
         mat(k,454) = .630_r8*rxt(k,279)*y(k,134)
         mat(k,842) = .560_r8*rxt(k,308)*y(k,134)
         mat(k,1900) = rxt(k,210)*y(k,17) + rxt(k,174)*y(k,56) + rxt(k,255)*y(k,126) &
                      + rxt(k,256)*y(k,133) + rxt(k,257)*y(k,212)
         mat(k,1071) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,212)
         mat(k,1720) = rxt(k,174)*y(k,42)
         mat(k,740) = rxt(k,302)*y(k,212)
         mat(k,711) = .620_r8*rxt(k,385)*y(k,134)
         mat(k,1059) = .650_r8*rxt(k,338)*y(k,134)
         mat(k,792) = .230_r8*rxt(k,443)*y(k,134)
         mat(k,1126) = .560_r8*rxt(k,352)*y(k,134)
         mat(k,1963) = .170_r8*rxt(k,411)*y(k,197) + .220_r8*rxt(k,336)*y(k,205) &
                      + .400_r8*rxt(k,414)*y(k,206) + .350_r8*rxt(k,417)*y(k,208) &
                      + .225_r8*rxt(k,452)*y(k,216) + .250_r8*rxt(k,393)*y(k,219)
         mat(k,1425) = rxt(k,255)*y(k,42) + rxt(k,314)*y(k,49) + .220_r8*rxt(k,335) &
                      *y(k,205) + .500_r8*rxt(k,394)*y(k,219)
         mat(k,1876) = rxt(k,256)*y(k,42) + rxt(k,462)*y(k,137)
         mat(k,1808) = .230_r8*rxt(k,440)*y(k,6) + .630_r8*rxt(k,279)*y(k,25) &
                      + .560_r8*rxt(k,308)*y(k,29) + .620_r8*rxt(k,385)*y(k,98) &
                      + .650_r8*rxt(k,338)*y(k,105) + .230_r8*rxt(k,443)*y(k,110) &
                      + .560_r8*rxt(k,352)*y(k,111)
         mat(k,253) = rxt(k,462)*y(k,133) + rxt(k,463)*y(k,212)
         mat(k,887) = .700_r8*rxt(k,461)*y(k,212)
         mat(k,1225) = .220_r8*rxt(k,332)*y(k,205) + .250_r8*rxt(k,390)*y(k,219)
         mat(k,1648) = .110_r8*rxt(k,333)*y(k,205) + .125_r8*rxt(k,450)*y(k,216) &
                      + .200_r8*rxt(k,391)*y(k,219)
         mat(k,622) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,201)
         mat(k,1368) = .070_r8*rxt(k,410)*y(k,197) + .160_r8*rxt(k,413)*y(k,206) &
                      + .140_r8*rxt(k,416)*y(k,208)
         mat(k,1148) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,195) + .110_r8*rxt(k,333)*y(k,196)
         mat(k,585) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,201)
         mat(k,745) = .350_r8*rxt(k,417)*y(k,124) + .140_r8*rxt(k,416)*y(k,201)
         mat(k,1595) = mat(k,1595) + .350_r8*rxt(k,277)*y(k,24) + rxt(k,257)*y(k,42) &
                      + rxt(k,315)*y(k,49) + rxt(k,302)*y(k,75) + rxt(k,463)*y(k,137) &
                      + .700_r8*rxt(k,461)*y(k,178)
         mat(k,949) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,196)
         mat(k,1022) = .250_r8*rxt(k,393)*y(k,124) + .500_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,195) + .200_r8*rxt(k,391)*y(k,196)
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
         mat(k,59) = -(rxt(k,235)*y(k,211))
         mat(k,1461) = -rxt(k,235)*y(k,64)
         mat(k,100) = -(rxt(k,236)*y(k,211))
         mat(k,1463) = -rxt(k,236)*y(k,65)
         mat(k,116) = -(rxt(k,409)*y(k,212))
         mat(k,1503) = -rxt(k,409)*y(k,66)
         mat(k,110) = .180_r8*rxt(k,429)*y(k,212)
         mat(k,1503) = mat(k,1503) + .180_r8*rxt(k,429)*y(k,180)
         mat(k,193) = -(rxt(k,476)*y(k,126) + (rxt(k,477) + rxt(k,479)) * y(k,212))
         mat(k,1407) = -rxt(k,476)*y(k,67)
         mat(k,1515) = -(rxt(k,477) + rxt(k,479)) * y(k,67)
         mat(k,574) = rxt(k,297)*y(k,201)
         mat(k,1310) = rxt(k,297)*y(k,200)
         mat(k,649) = -(rxt(k,232)*y(k,77) + rxt(k,233)*y(k,222) + rxt(k,234)*y(k,89))
         mat(k,1006) = -rxt(k,232)*y(k,73)
         mat(k,2001) = -rxt(k,233)*y(k,73)
         mat(k,1261) = -rxt(k,234)*y(k,73)
         mat(k,60) = 2.000_r8*rxt(k,235)*y(k,211)
         mat(k,101) = rxt(k,236)*y(k,211)
         mat(k,1464) = 2.000_r8*rxt(k,235)*y(k,64) + rxt(k,236)*y(k,65)
         mat(k,865) = -(rxt(k,301)*y(k,212))
         mat(k,1591) = -rxt(k,301)*y(k,74)
         mat(k,468) = .700_r8*rxt(k,376)*y(k,212)
         mat(k,425) = .500_r8*rxt(k,377)*y(k,212)
         mat(k,280) = rxt(k,388)*y(k,212)
         mat(k,1959) = .050_r8*rxt(k,374)*y(k,204) + .530_r8*rxt(k,336)*y(k,205) &
                      + .225_r8*rxt(k,452)*y(k,216) + .250_r8*rxt(k,393)*y(k,219)
         mat(k,1421) = .050_r8*rxt(k,375)*y(k,204) + .530_r8*rxt(k,335)*y(k,205) &
                      + .250_r8*rxt(k,394)*y(k,219)
         mat(k,1223) = .530_r8*rxt(k,332)*y(k,205) + .250_r8*rxt(k,390)*y(k,219)
         mat(k,1645) = .260_r8*rxt(k,333)*y(k,205) + .125_r8*rxt(k,450)*y(k,216) &
                      + .100_r8*rxt(k,391)*y(k,219)
         mat(k,1175) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1146) = .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335)*y(k,126) &
                      + .530_r8*rxt(k,332)*y(k,195) + .260_r8*rxt(k,333)*y(k,196)
         mat(k,1591) = mat(k,1591) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + rxt(k,388)*y(k,115)
         mat(k,947) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,196)
         mat(k,1021) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,195) + .100_r8*rxt(k,391)*y(k,196)
         mat(k,739) = -(rxt(k,302)*y(k,212))
         mat(k,1582) = -rxt(k,302)*y(k,75)
         mat(k,201) = .650_r8*rxt(k,277)*y(k,212)
         mat(k,864) = .200_r8*rxt(k,301)*y(k,212)
         mat(k,872) = rxt(k,389)*y(k,212)
         mat(k,1954) = rxt(k,400)*y(k,190) + .050_r8*rxt(k,374)*y(k,204) &
                      + .400_r8*rxt(k,414)*y(k,206) + .170_r8*rxt(k,417)*y(k,208) &
                      + .700_r8*rxt(k,420)*y(k,213) + .600_r8*rxt(k,427)*y(k,218) &
                      + .250_r8*rxt(k,393)*y(k,219) + .340_r8*rxt(k,433)*y(k,220) &
                      + .170_r8*rxt(k,436)*y(k,221)
         mat(k,1414) = .050_r8*rxt(k,375)*y(k,204) + .250_r8*rxt(k,394)*y(k,219)
         mat(k,378) = rxt(k,400)*y(k,124)
         mat(k,1222) = .250_r8*rxt(k,390)*y(k,219)
         mat(k,1642) = .100_r8*rxt(k,391)*y(k,219)
         mat(k,1361) = .160_r8*rxt(k,413)*y(k,206) + .070_r8*rxt(k,416)*y(k,208)
         mat(k,1173) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,584) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,201)
         mat(k,743) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,201)
         mat(k,1582) = mat(k,1582) + .650_r8*rxt(k,277)*y(k,24) + .200_r8*rxt(k,301) &
                      *y(k,74) + rxt(k,389)*y(k,116)
         mat(k,348) = .700_r8*rxt(k,420)*y(k,124)
         mat(k,596) = .600_r8*rxt(k,427)*y(k,124)
         mat(k,1020) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,195) + .100_r8*rxt(k,391)*y(k,196)
         mat(k,612) = .340_r8*rxt(k,433)*y(k,124)
         mat(k,385) = .170_r8*rxt(k,436)*y(k,124)
         mat(k,1276) = -((rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,201) + rxt(k,140) &
                      *y(k,134))
         mat(k,1388) = -(rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,76)
         mat(k,1827) = -rxt(k,140)*y(k,76)
         mat(k,1904) = rxt(k,257)*y(k,212)
         mat(k,1726) = rxt(k,176)*y(k,77)
         mat(k,912) = rxt(k,275)*y(k,212)
         mat(k,652) = rxt(k,232)*y(k,77)
         mat(k,1009) = rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73) + rxt(k,132)*y(k,133) &
                      + rxt(k,124)*y(k,211) + rxt(k,141)*y(k,212)
         mat(k,697) = rxt(k,230)*y(k,211)
         mat(k,1295) = rxt(k,207)*y(k,211)
         mat(k,267) = rxt(k,162)*y(k,212)
         mat(k,1881) = rxt(k,132)*y(k,77) + rxt(k,144)*y(k,212)
         mat(k,255) = rxt(k,463)*y(k,212)
         mat(k,399) = rxt(k,468)*y(k,212)
         mat(k,1085) = rxt(k,473)*y(k,212)
         mat(k,1469) = rxt(k,124)*y(k,77) + rxt(k,230)*y(k,81) + rxt(k,207)*y(k,85)
         mat(k,1616) = rxt(k,257)*y(k,42) + rxt(k,275)*y(k,62) + rxt(k,141)*y(k,77) &
                      + rxt(k,162)*y(k,112) + rxt(k,144)*y(k,133) + rxt(k,463) &
                      *y(k,137) + rxt(k,468)*y(k,148) + rxt(k,473)*y(k,150)
         mat(k,1007) = -(rxt(k,124)*y(k,211) + rxt(k,132)*y(k,133) + rxt(k,141) &
                      *y(k,212) + rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73))
         mat(k,1466) = -rxt(k,124)*y(k,77)
         mat(k,1877) = -rxt(k,132)*y(k,77)
         mat(k,1602) = -rxt(k,141)*y(k,77)
         mat(k,1722) = -rxt(k,176)*y(k,77)
         mat(k,650) = -rxt(k,232)*y(k,77)
         mat(k,1274) = rxt(k,134)*y(k,201)
         mat(k,1375) = rxt(k,134)*y(k,76)
         mat(k,476) = -(rxt(k,133)*y(k,133) + rxt(k,142)*y(k,212) + rxt(k,177)*y(k,56))
         mat(k,1870) = -rxt(k,133)*y(k,79)
         mat(k,1555) = -rxt(k,142)*y(k,79)
         mat(k,1712) = -rxt(k,177)*y(k,79)
         mat(k,1339) = 2.000_r8*rxt(k,148)*y(k,201)
         mat(k,1555) = mat(k,1555) + 2.000_r8*rxt(k,147)*y(k,212)
         mat(k,167) = rxt(k,475)*y(k,222)
         mat(k,1998) = rxt(k,475)*y(k,152)
         mat(k,695) = -(rxt(k,223)*y(k,133) + rxt(k,224)*y(k,212) + (rxt(k,229) &
                      + rxt(k,230)) * y(k,211))
         mat(k,1873) = -rxt(k,223)*y(k,81)
         mat(k,1578) = -rxt(k,224)*y(k,81)
         mat(k,1465) = -(rxt(k,229) + rxt(k,230)) * y(k,81)
         mat(k,1248) = rxt(k,210)*y(k,42) + rxt(k,211)*y(k,201)
         mat(k,1899) = rxt(k,210)*y(k,17)
         mat(k,1359) = rxt(k,211)*y(k,17)
         mat(k,1296) = -(rxt(k,194)*y(k,133) + rxt(k,195)*y(k,212) + (rxt(k,206) &
                      + rxt(k,207)) * y(k,211) + (rxt(k,509) + rxt(k,515) + rxt(k,520) &
                      ) * y(k,92) + (rxt(k,514) + rxt(k,519) + rxt(k,524)) * y(k,60) &
                      + (rxt(k,516) + rxt(k,521)) * y(k,91))
         mat(k,1882) = -rxt(k,194)*y(k,85)
         mat(k,1617) = -rxt(k,195)*y(k,85)
         mat(k,1470) = -(rxt(k,206) + rxt(k,207)) * y(k,85)
         mat(k,678) = -(rxt(k,509) + rxt(k,515) + rxt(k,520)) * y(k,85)
         mat(k,829) = -(rxt(k,514) + rxt(k,519) + rxt(k,524)) * y(k,85)
         mat(k,642) = -(rxt(k,516) + rxt(k,521)) * y(k,85)
         mat(k,207) = rxt(k,286)*y(k,56)
         mat(k,1905) = rxt(k,174)*y(k,56)
         mat(k,1727) = rxt(k,286)*y(k,28) + rxt(k,174)*y(k,42) + rxt(k,176)*y(k,77) &
                      + rxt(k,177)*y(k,79) + rxt(k,196)*y(k,92) + rxt(k,178)*y(k,201)
         mat(k,1693) = rxt(k,193)*y(k,212)
         mat(k,1010) = rxt(k,176)*y(k,56)
         mat(k,477) = rxt(k,177)*y(k,56)
         mat(k,678) = mat(k,678) + rxt(k,196)*y(k,56)
         mat(k,1389) = rxt(k,178)*y(k,56)
         mat(k,1617) = mat(k,1617) + rxt(k,193)*y(k,59)
         mat(k,96) = -(rxt(k,266)*y(k,212) + rxt(k,274)*y(k,211))
         mat(k,1500) = -rxt(k,266)*y(k,86)
         mat(k,1462) = -rxt(k,274)*y(k,86)
         mat(k,657) = -(rxt(k,267)*y(k,212))
         mat(k,1573) = -rxt(k,267)*y(k,87)
         mat(k,759) = .050_r8*rxt(k,440)*y(k,134)
         mat(k,200) = .350_r8*rxt(k,277)*y(k,212)
         mat(k,453) = .370_r8*rxt(k,279)*y(k,134)
         mat(k,839) = .120_r8*rxt(k,308)*y(k,134)
         mat(k,708) = .110_r8*rxt(k,385)*y(k,134)
         mat(k,1058) = .330_r8*rxt(k,338)*y(k,134)
         mat(k,785) = .050_r8*rxt(k,443)*y(k,134)
         mat(k,1124) = .120_r8*rxt(k,352)*y(k,134)
         mat(k,1950) = rxt(k,270)*y(k,202)
         mat(k,1796) = .050_r8*rxt(k,440)*y(k,6) + .370_r8*rxt(k,279)*y(k,25) &
                      + .120_r8*rxt(k,308)*y(k,29) + .110_r8*rxt(k,385)*y(k,98) &
                      + .330_r8*rxt(k,338)*y(k,105) + .050_r8*rxt(k,443)*y(k,110) &
                      + .120_r8*rxt(k,352)*y(k,111)
         mat(k,1355) = rxt(k,268)*y(k,202)
         mat(k,341) = rxt(k,270)*y(k,124) + rxt(k,268)*y(k,201)
         mat(k,1573) = mat(k,1573) + .350_r8*rxt(k,277)*y(k,24)
         mat(k,648) = rxt(k,232)*y(k,77) + rxt(k,234)*y(k,89) + rxt(k,233)*y(k,222)
         mat(k,1005) = rxt(k,232)*y(k,73)
         mat(k,1260) = rxt(k,234)*y(k,73)
         mat(k,1999) = rxt(k,233)*y(k,73)
         mat(k,1263) = -(rxt(k,171)*y(k,212) + rxt(k,234)*y(k,73))
         mat(k,1615) = -rxt(k,171)*y(k,89)
         mat(k,651) = -rxt(k,234)*y(k,89)
         mat(k,1903) = rxt(k,255)*y(k,126)
         mat(k,998) = rxt(k,288)*y(k,126)
         mat(k,1074) = rxt(k,314)*y(k,126)
         mat(k,828) = (rxt(k,514)+rxt(k,519)+rxt(k,524))*y(k,85)
         mat(k,195) = rxt(k,476)*y(k,126)
         mat(k,1294) = (rxt(k,514)+rxt(k,519)+rxt(k,524))*y(k,60)
         mat(k,1766) = rxt(k,170)*y(k,212)
         mat(k,1444) = rxt(k,255)*y(k,42) + rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) &
                      + rxt(k,476)*y(k,67)
         mat(k,1615) = mat(k,1615) + rxt(k,170)*y(k,125)
         mat(k,272) = -(rxt(k,149)*y(k,212))
         mat(k,1527) = -rxt(k,149)*y(k,90)
         mat(k,1745) = rxt(k,168)*y(k,201)
         mat(k,1319) = rxt(k,168)*y(k,125)
         mat(k,640) = -(rxt(k,225)*y(k,133) + (rxt(k,516) + rxt(k,521)) * y(k,85))
         mat(k,1871) = -rxt(k,225)*y(k,91)
         mat(k,1290) = -(rxt(k,516) + rxt(k,521)) * y(k,91)
         mat(k,1846) = rxt(k,217)*y(k,201)
         mat(k,1354) = rxt(k,217)*y(k,19)
         mat(k,677) = -(rxt(k,196)*y(k,56) + rxt(k,197)*y(k,133) + rxt(k,198)*y(k,212) &
                      + (rxt(k,509) + rxt(k,515) + rxt(k,520)) * y(k,85))
         mat(k,1715) = -rxt(k,196)*y(k,92)
         mat(k,1872) = -rxt(k,197)*y(k,92)
         mat(k,1576) = -rxt(k,198)*y(k,92)
         mat(k,1291) = -(rxt(k,509) + rxt(k,515) + rxt(k,520)) * y(k,92)
         mat(k,1686) = rxt(k,185)*y(k,201)
         mat(k,826) = rxt(k,190)*y(k,212)
         mat(k,1357) = rxt(k,185)*y(k,59)
         mat(k,1576) = mat(k,1576) + rxt(k,190)*y(k,60)
         mat(k,898) = -(rxt(k,331)*y(k,212))
         mat(k,1594) = -rxt(k,331)*y(k,93)
         mat(k,469) = .300_r8*rxt(k,376)*y(k,212)
         mat(k,426) = .500_r8*rxt(k,377)*y(k,212)
         mat(k,1962) = rxt(k,330)*y(k,198) + rxt(k,337)*y(k,205)
         mat(k,461) = rxt(k,330)*y(k,124)
         mat(k,1147) = rxt(k,337)*y(k,124)
         mat(k,1594) = mat(k,1594) + .300_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100)
         mat(k,150) = -(rxt(k,362)*y(k,212))
         mat(k,1508) = -rxt(k,362)*y(k,94)
         mat(k,915) = -(rxt(k,316)*y(k,212))
         mat(k,1596) = -rxt(k,316)*y(k,95)
         mat(k,470) = .700_r8*rxt(k,376)*y(k,212)
         mat(k,427) = .500_r8*rxt(k,377)*y(k,212)
         mat(k,416) = .500_r8*rxt(k,351)*y(k,212)
         mat(k,1964) = .050_r8*rxt(k,374)*y(k,204) + .220_r8*rxt(k,336)*y(k,205) &
                      + .250_r8*rxt(k,393)*y(k,219)
         mat(k,1426) = .050_r8*rxt(k,375)*y(k,204) + .220_r8*rxt(k,335)*y(k,205) &
                      + .250_r8*rxt(k,394)*y(k,219)
         mat(k,437) = .500_r8*rxt(k,320)*y(k,212)
         mat(k,1226) = .220_r8*rxt(k,332)*y(k,205) + .250_r8*rxt(k,390)*y(k,219)
         mat(k,1649) = .230_r8*rxt(k,333)*y(k,205) + .200_r8*rxt(k,321)*y(k,215) &
                      + .100_r8*rxt(k,391)*y(k,219)
         mat(k,1178) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1149) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,195) + .230_r8*rxt(k,333)*y(k,196)
         mat(k,1596) = mat(k,1596) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + .500_r8*rxt(k,351)*y(k,109) + .500_r8*rxt(k,320) &
                      *y(k,146)
         mat(k,984) = .200_r8*rxt(k,321)*y(k,196)
         mat(k,1023) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,195) + .100_r8*rxt(k,391)*y(k,196)
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
         mat(k,211) = -(rxt(k,363)*y(k,212))
         mat(k,1518) = -rxt(k,363)*y(k,96)
         mat(k,1924) = .870_r8*rxt(k,374)*y(k,204)
         mat(k,1408) = .950_r8*rxt(k,375)*y(k,204)
         mat(k,1218) = rxt(k,370)*y(k,204)
         mat(k,1633) = .750_r8*rxt(k,371)*y(k,204)
         mat(k,1167) = .870_r8*rxt(k,374)*y(k,124) + .950_r8*rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,195) + .750_r8*rxt(k,371)*y(k,196)
         mat(k,68) = -(rxt(k,364)*y(k,212))
         mat(k,1496) = -rxt(k,364)*y(k,97)
         mat(k,551) = .600_r8*rxt(k,387)*y(k,212)
         mat(k,1496) = mat(k,1496) + .600_r8*rxt(k,387)*y(k,103)
         mat(k,709) = -(rxt(k,378)*y(k,126) + rxt(k,385)*y(k,134) + rxt(k,386) &
                      *y(k,212))
         mat(k,1411) = -rxt(k,378)*y(k,98)
         mat(k,1797) = -rxt(k,385)*y(k,98)
         mat(k,1579) = -rxt(k,386)*y(k,98)
         mat(k,467) = -(rxt(k,376)*y(k,212))
         mat(k,1554) = -rxt(k,376)*y(k,99)
         mat(k,1938) = .080_r8*rxt(k,368)*y(k,203)
         mat(k,1096) = .080_r8*rxt(k,368)*y(k,124)
         mat(k,423) = -(rxt(k,377)*y(k,212))
         mat(k,1549) = -rxt(k,377)*y(k,100)
         mat(k,1936) = .080_r8*rxt(k,374)*y(k,204)
         mat(k,1168) = .080_r8*rxt(k,374)*y(k,124)
         mat(k,1044) = -(rxt(k,379)*y(k,195) + rxt(k,380)*y(k,196) + rxt(k,381) &
                      *y(k,201) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126))
         mat(k,1228) = -rxt(k,379)*y(k,101)
         mat(k,1656) = -rxt(k,380)*y(k,101)
         mat(k,1377) = -rxt(k,381)*y(k,101)
         mat(k,1971) = -rxt(k,382)*y(k,101)
         mat(k,1433) = -rxt(k,383)*y(k,101)
         mat(k,712) = rxt(k,378)*y(k,126)
         mat(k,1433) = mat(k,1433) + rxt(k,378)*y(k,98)
         mat(k,308) = -(rxt(k,384)*y(k,212))
         mat(k,1532) = -rxt(k,384)*y(k,102)
         mat(k,1036) = rxt(k,381)*y(k,201)
         mat(k,1323) = rxt(k,381)*y(k,101)
         mat(k,552) = -(rxt(k,387)*y(k,212))
         mat(k,1564) = -rxt(k,387)*y(k,103)
         mat(k,1346) = rxt(k,367)*y(k,203) + rxt(k,372)*y(k,204)
         mat(k,1097) = rxt(k,367)*y(k,201)
         mat(k,1170) = rxt(k,372)*y(k,201)
         mat(k,39) = -(rxt(k,501)*y(k,212))
         mat(k,1490) = -rxt(k,501)*y(k,104)
         mat(k,1060) = -(rxt(k,338)*y(k,134) + rxt(k,339)*y(k,212))
         mat(k,1816) = -rxt(k,338)*y(k,105)
         mat(k,1605) = -rxt(k,339)*y(k,105)
         mat(k,713) = .300_r8*rxt(k,385)*y(k,134)
         mat(k,1972) = .360_r8*rxt(k,368)*y(k,203)
         mat(k,1434) = .400_r8*rxt(k,369)*y(k,203)
         mat(k,1816) = mat(k,1816) + .300_r8*rxt(k,385)*y(k,98)
         mat(k,1229) = .390_r8*rxt(k,365)*y(k,203)
         mat(k,1657) = .310_r8*rxt(k,366)*y(k,203)
         mat(k,1106) = .360_r8*rxt(k,368)*y(k,124) + .400_r8*rxt(k,369)*y(k,126) &
                      + .390_r8*rxt(k,365)*y(k,195) + .310_r8*rxt(k,366)*y(k,196)
         mat(k,214) = -(rxt(k,340)*y(k,212))
         mat(k,1519) = -rxt(k,340)*y(k,106)
         mat(k,1314) = rxt(k,334)*y(k,205)
         mat(k,1144) = rxt(k,334)*y(k,201)
         mat(k,403) = -(rxt(k,349)*y(k,212))
         mat(k,1546) = -rxt(k,349)*y(k,107)
         mat(k,1934) = .800_r8*rxt(k,358)*y(k,189)
         mat(k,804) = .800_r8*rxt(k,358)*y(k,124)
         mat(k,219) = -(rxt(k,350)*y(k,212))
         mat(k,1520) = -rxt(k,350)*y(k,108)
         mat(k,1315) = .800_r8*rxt(k,347)*y(k,209)
         mat(k,543) = .800_r8*rxt(k,347)*y(k,201)
         mat(k,415) = -(rxt(k,351)*y(k,212))
         mat(k,1548) = -rxt(k,351)*y(k,109)
         mat(k,1749) = rxt(k,354)*y(k,207)
         mat(k,1201) = rxt(k,354)*y(k,125)
         mat(k,787) = -(rxt(k,442)*y(k,126) + rxt(k,443)*y(k,134) + rxt(k,444) &
                      *y(k,212))
         mat(k,1416) = -rxt(k,442)*y(k,110)
         mat(k,1801) = -rxt(k,443)*y(k,110)
         mat(k,1585) = -rxt(k,444)*y(k,110)
         mat(k,1130) = -(rxt(k,352)*y(k,134) + rxt(k,353)*y(k,212))
         mat(k,1820) = -rxt(k,352)*y(k,111)
         mat(k,1609) = -rxt(k,353)*y(k,111)
         mat(k,715) = .200_r8*rxt(k,385)*y(k,134)
         mat(k,1975) = .560_r8*rxt(k,368)*y(k,203)
         mat(k,1438) = .600_r8*rxt(k,369)*y(k,203)
         mat(k,1820) = mat(k,1820) + .200_r8*rxt(k,385)*y(k,98)
         mat(k,1232) = .610_r8*rxt(k,365)*y(k,203)
         mat(k,1660) = .440_r8*rxt(k,366)*y(k,203)
         mat(k,1109) = .560_r8*rxt(k,368)*y(k,124) + .600_r8*rxt(k,369)*y(k,126) &
                      + .610_r8*rxt(k,365)*y(k,195) + .440_r8*rxt(k,366)*y(k,196)
         mat(k,266) = -(rxt(k,150)*y(k,124) + (rxt(k,151) + rxt(k,152) + rxt(k,153) &
                      ) * y(k,125) + rxt(k,162)*y(k,212))
         mat(k,1926) = -rxt(k,150)*y(k,112)
         mat(k,1744) = -(rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,112)
         mat(k,1526) = -rxt(k,162)*y(k,112)
         mat(k,1743) = rxt(k,169)*y(k,126)
         mat(k,1406) = rxt(k,169)*y(k,125)
         mat(k,278) = -(rxt(k,388)*y(k,212))
         mat(k,1528) = -rxt(k,388)*y(k,115)
         mat(k,1035) = .200_r8*rxt(k,380)*y(k,196)
         mat(k,1634) = .200_r8*rxt(k,380)*y(k,101)
         mat(k,874) = -(rxt(k,389)*y(k,212))
         mat(k,1592) = -rxt(k,389)*y(k,116)
         mat(k,1041) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) + rxt(k,379)*y(k,195) &
                      + .800_r8*rxt(k,380)*y(k,196)
         mat(k,1960) = rxt(k,382)*y(k,101)
         mat(k,1422) = rxt(k,383)*y(k,101)
         mat(k,1224) = rxt(k,379)*y(k,101)
         mat(k,1646) = .800_r8*rxt(k,380)*y(k,101)
         mat(k,49) = -(rxt(k,478)*y(k,212))
         mat(k,1492) = -rxt(k,478)*y(k,120)
         mat(k,1996) = -(rxt(k,150)*y(k,112) + rxt(k,159)*y(k,126) + rxt(k,163) &
                      *y(k,201) + rxt(k,164)*y(k,134) + rxt(k,165)*y(k,133) + rxt(k,186) &
                      *y(k,59) + rxt(k,218)*y(k,19) + rxt(k,261)*y(k,196) + rxt(k,270) &
                      *y(k,202) + rxt(k,283)*y(k,192) + rxt(k,294)*y(k,195) + rxt(k,298) &
                      *y(k,200) + rxt(k,311)*y(k,193) + rxt(k,319)*y(k,214) + rxt(k,323) &
                      *y(k,215) + (rxt(k,329) + rxt(k,330)) * y(k,198) + (rxt(k,336) &
                      + rxt(k,337)) * y(k,205) + rxt(k,345)*y(k,207) + rxt(k,348) &
                      *y(k,209) + (rxt(k,358) + rxt(k,359)) * y(k,189) + rxt(k,368) &
                      *y(k,203) + rxt(k,374)*y(k,204) + rxt(k,382)*y(k,101) + rxt(k,393) &
                      *y(k,219) + rxt(k,397)*y(k,188) + rxt(k,400)*y(k,190) + rxt(k,405) &
                      *y(k,191) + rxt(k,407)*y(k,194) + rxt(k,411)*y(k,197) + rxt(k,414) &
                      *y(k,206) + rxt(k,417)*y(k,208) + rxt(k,420)*y(k,213) + rxt(k,427) &
                      *y(k,218) + rxt(k,433)*y(k,220) + rxt(k,436)*y(k,221) + rxt(k,447) &
                      *y(k,210) + rxt(k,452)*y(k,216) + rxt(k,457)*y(k,217))
         mat(k,271) = -rxt(k,150)*y(k,124)
         mat(k,1459) = -rxt(k,159)*y(k,124)
         mat(k,1402) = -rxt(k,163)*y(k,124)
         mat(k,1841) = -rxt(k,164)*y(k,124)
         mat(k,1895) = -rxt(k,165)*y(k,124)
         mat(k,1706) = -rxt(k,186)*y(k,124)
         mat(k,1865) = -rxt(k,218)*y(k,124)
         mat(k,1680) = -rxt(k,261)*y(k,124)
         mat(k,345) = -rxt(k,270)*y(k,124)
         mat(k,694) = -rxt(k,283)*y(k,124)
         mat(k,1246) = -rxt(k,294)*y(k,124)
         mat(k,582) = -rxt(k,298)*y(k,124)
         mat(k,674) = -rxt(k,311)*y(k,124)
         mat(k,637) = -rxt(k,319)*y(k,124)
         mat(k,993) = -rxt(k,323)*y(k,124)
         mat(k,466) = -(rxt(k,329) + rxt(k,330)) * y(k,124)
         mat(k,1164) = -(rxt(k,336) + rxt(k,337)) * y(k,124)
         mat(k,1216) = -rxt(k,345)*y(k,124)
         mat(k,550) = -rxt(k,348)*y(k,124)
         mat(k,818) = -(rxt(k,358) + rxt(k,359)) * y(k,124)
         mat(k,1122) = -rxt(k,368)*y(k,124)
         mat(k,1198) = -rxt(k,374)*y(k,124)
         mat(k,1057) = -rxt(k,382)*y(k,124)
         mat(k,1034) = -rxt(k,393)*y(k,124)
         mat(k,414) = -rxt(k,397)*y(k,124)
         mat(k,382) = -rxt(k,400)*y(k,124)
         mat(k,339) = -rxt(k,405)*y(k,124)
         mat(k,520) = -rxt(k,407)*y(k,124)
         mat(k,628) = -rxt(k,411)*y(k,124)
         mat(k,588) = -rxt(k,414)*y(k,124)
         mat(k,751) = -rxt(k,417)*y(k,124)
         mat(k,352) = -rxt(k,420)*y(k,124)
         mat(k,603) = -rxt(k,427)*y(k,124)
         mat(k,620) = -rxt(k,433)*y(k,124)
         mat(k,390) = -rxt(k,436)*y(k,124)
         mat(k,980) = -rxt(k,447)*y(k,124)
         mat(k,960) = -rxt(k,452)*y(k,124)
         mat(k,941) = -rxt(k,457)*y(k,124)
         mat(k,271) = mat(k,271) + 2.000_r8*rxt(k,152)*y(k,125) + rxt(k,162)*y(k,212)
         mat(k,1781) = 2.000_r8*rxt(k,152)*y(k,112) + rxt(k,155)*y(k,133) + rxt(k,469) &
                      *y(k,150)
         mat(k,1895) = mat(k,1895) + rxt(k,155)*y(k,125)
         mat(k,1094) = rxt(k,469)*y(k,125)
         mat(k,1630) = rxt(k,162)*y(k,112)
         mat(k,1776) = -((rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,112) + (rxt(k,155) &
                      + rxt(k,157)) * y(k,133) + rxt(k,156)*y(k,134) + rxt(k,168) &
                      *y(k,201) + rxt(k,169)*y(k,126) + rxt(k,170)*y(k,212) + rxt(k,188) &
                      *y(k,59) + rxt(k,219)*y(k,19) + rxt(k,305)*y(k,195) + rxt(k,354) &
                      *y(k,207) + rxt(k,412)*y(k,197) + rxt(k,415)*y(k,206) + rxt(k,418) &
                      *y(k,208) + rxt(k,422)*y(k,141) + rxt(k,425)*y(k,188) + rxt(k,469) &
                      *y(k,150))
         mat(k,269) = -(rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,125)
         mat(k,1890) = -(rxt(k,155) + rxt(k,157)) * y(k,125)
         mat(k,1836) = -rxt(k,156)*y(k,125)
         mat(k,1397) = -rxt(k,168)*y(k,125)
         mat(k,1454) = -rxt(k,169)*y(k,125)
         mat(k,1625) = -rxt(k,170)*y(k,125)
         mat(k,1701) = -rxt(k,188)*y(k,125)
         mat(k,1860) = -rxt(k,219)*y(k,125)
         mat(k,1243) = -rxt(k,305)*y(k,125)
         mat(k,1213) = -rxt(k,354)*y(k,125)
         mat(k,627) = -rxt(k,412)*y(k,125)
         mat(k,587) = -rxt(k,415)*y(k,125)
         mat(k,750) = -rxt(k,418)*y(k,125)
         mat(k,365) = -rxt(k,422)*y(k,125)
         mat(k,413) = -rxt(k,425)*y(k,125)
         mat(k,1090) = -rxt(k,469)*y(k,125)
         mat(k,541) = rxt(k,356)*y(k,212)
         mat(k,264) = rxt(k,327)*y(k,126)
         mat(k,1860) = mat(k,1860) + rxt(k,218)*y(k,124)
         mat(k,1701) = mat(k,1701) + rxt(k,186)*y(k,124)
         mat(k,276) = rxt(k,149)*y(k,212)
         mat(k,474) = .700_r8*rxt(k,376)*y(k,212)
         mat(k,1055) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126)
         mat(k,1991) = rxt(k,218)*y(k,19) + rxt(k,186)*y(k,59) + rxt(k,382)*y(k,101) &
                      + 2.000_r8*rxt(k,159)*y(k,126) + rxt(k,165)*y(k,133) &
                      + rxt(k,164)*y(k,134) + rxt(k,397)*y(k,188) + rxt(k,358) &
                      *y(k,189) + rxt(k,400)*y(k,190) + rxt(k,405)*y(k,191) &
                      + rxt(k,283)*y(k,192) + rxt(k,311)*y(k,193) + rxt(k,407) &
                      *y(k,194) + rxt(k,294)*y(k,195) + rxt(k,261)*y(k,196) &
                      + rxt(k,411)*y(k,197) + rxt(k,329)*y(k,198) + rxt(k,298) &
                      *y(k,200) + rxt(k,163)*y(k,201) + rxt(k,270)*y(k,202) &
                      + .920_r8*rxt(k,368)*y(k,203) + .920_r8*rxt(k,374)*y(k,204) &
                      + rxt(k,336)*y(k,205) + rxt(k,414)*y(k,206) + rxt(k,345) &
                      *y(k,207) + rxt(k,417)*y(k,208) + rxt(k,348)*y(k,209) &
                      + 1.600_r8*rxt(k,447)*y(k,210) + rxt(k,420)*y(k,213) &
                      + rxt(k,319)*y(k,214) + rxt(k,323)*y(k,215) + .900_r8*rxt(k,452) &
                      *y(k,216) + .800_r8*rxt(k,457)*y(k,217) + rxt(k,427)*y(k,218) &
                      + rxt(k,393)*y(k,219) + rxt(k,433)*y(k,220) + rxt(k,436) &
                      *y(k,221)
         mat(k,1454) = mat(k,1454) + rxt(k,327)*y(k,16) + rxt(k,383)*y(k,101) &
                      + 2.000_r8*rxt(k,159)*y(k,124) + rxt(k,160)*y(k,133) &
                      + rxt(k,158)*y(k,201) + rxt(k,369)*y(k,203) + rxt(k,375) &
                      *y(k,204) + rxt(k,335)*y(k,205) + rxt(k,346)*y(k,207) &
                      + 2.000_r8*rxt(k,448)*y(k,210) + rxt(k,161)*y(k,212) &
                      + rxt(k,394)*y(k,219)
         mat(k,736) = rxt(k,317)*y(k,212)
         mat(k,1890) = mat(k,1890) + rxt(k,165)*y(k,124) + rxt(k,160)*y(k,126)
         mat(k,1836) = mat(k,1836) + rxt(k,164)*y(k,124)
         mat(k,512) = rxt(k,454)*y(k,212)
         mat(k,413) = mat(k,413) + rxt(k,397)*y(k,124)
         mat(k,816) = rxt(k,358)*y(k,124)
         mat(k,381) = rxt(k,400)*y(k,124)
         mat(k,338) = rxt(k,405)*y(k,124)
         mat(k,692) = rxt(k,283)*y(k,124)
         mat(k,672) = rxt(k,311)*y(k,124)
         mat(k,518) = rxt(k,407)*y(k,124)
         mat(k,1243) = mat(k,1243) + rxt(k,294)*y(k,124)
         mat(k,1675) = rxt(k,261)*y(k,124) + .500_r8*rxt(k,445)*y(k,210)
         mat(k,627) = mat(k,627) + rxt(k,411)*y(k,124)
         mat(k,464) = rxt(k,329)*y(k,124)
         mat(k,580) = rxt(k,298)*y(k,124)
         mat(k,1397) = mat(k,1397) + rxt(k,163)*y(k,124) + rxt(k,158)*y(k,126)
         mat(k,343) = rxt(k,270)*y(k,124)
         mat(k,1119) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126)
         mat(k,1195) = .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,1162) = rxt(k,336)*y(k,124) + rxt(k,335)*y(k,126)
         mat(k,587) = mat(k,587) + rxt(k,414)*y(k,124)
         mat(k,1213) = mat(k,1213) + rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126)
         mat(k,750) = mat(k,750) + rxt(k,417)*y(k,124)
         mat(k,549) = rxt(k,348)*y(k,124)
         mat(k,978) = 1.600_r8*rxt(k,447)*y(k,124) + 2.000_r8*rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,196)
         mat(k,1625) = mat(k,1625) + rxt(k,356)*y(k,1) + rxt(k,149)*y(k,90) &
                      + .700_r8*rxt(k,376)*y(k,99) + rxt(k,161)*y(k,126) + rxt(k,317) &
                      *y(k,127) + rxt(k,454)*y(k,175)
         mat(k,351) = rxt(k,420)*y(k,124)
         mat(k,635) = rxt(k,319)*y(k,124)
         mat(k,991) = rxt(k,323)*y(k,124)
         mat(k,958) = .900_r8*rxt(k,452)*y(k,124)
         mat(k,939) = .800_r8*rxt(k,457)*y(k,124)
         mat(k,602) = rxt(k,427)*y(k,124)
         mat(k,1032) = rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126)
         mat(k,619) = rxt(k,433)*y(k,124)
         mat(k,389) = rxt(k,436)*y(k,124)
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
         mat(k,1448) = -(rxt(k,158)*y(k,201) + rxt(k,159)*y(k,124) + rxt(k,160) &
                      *y(k,133) + rxt(k,161)*y(k,212) + rxt(k,169)*y(k,125) + rxt(k,255) &
                      *y(k,42) + rxt(k,288)*y(k,45) + rxt(k,307)*y(k,29) + rxt(k,314) &
                      *y(k,49) + rxt(k,327)*y(k,16) + rxt(k,335)*y(k,205) + rxt(k,346) &
                      *y(k,207) + rxt(k,369)*y(k,203) + rxt(k,375)*y(k,204) + rxt(k,378) &
                      *y(k,98) + rxt(k,383)*y(k,101) + rxt(k,394)*y(k,219) + rxt(k,439) &
                      *y(k,6) + rxt(k,442)*y(k,110) + rxt(k,448)*y(k,210) + rxt(k,459) &
                      *y(k,177) + rxt(k,476)*y(k,67))
         mat(k,1391) = -rxt(k,158)*y(k,126)
         mat(k,1985) = -rxt(k,159)*y(k,126)
         mat(k,1884) = -rxt(k,160)*y(k,126)
         mat(k,1619) = -rxt(k,161)*y(k,126)
         mat(k,1770) = -rxt(k,169)*y(k,126)
         mat(k,1907) = -rxt(k,255)*y(k,126)
         mat(k,1001) = -rxt(k,288)*y(k,126)
         mat(k,848) = -rxt(k,307)*y(k,126)
         mat(k,1077) = -rxt(k,314)*y(k,126)
         mat(k,262) = -rxt(k,327)*y(k,126)
         mat(k,1159) = -rxt(k,335)*y(k,126)
         mat(k,1210) = -rxt(k,346)*y(k,126)
         mat(k,1116) = -rxt(k,369)*y(k,126)
         mat(k,1192) = -rxt(k,375)*y(k,126)
         mat(k,719) = -rxt(k,378)*y(k,126)
         mat(k,1052) = -rxt(k,383)*y(k,126)
         mat(k,1029) = -rxt(k,394)*y(k,126)
         mat(k,773) = -rxt(k,439)*y(k,126)
         mat(k,799) = -rxt(k,442)*y(k,126)
         mat(k,975) = -rxt(k,448)*y(k,126)
         mat(k,862) = -rxt(k,459)*y(k,126)
         mat(k,197) = -rxt(k,476)*y(k,126)
         mat(k,447) = rxt(k,220)*y(k,133)
         mat(k,1729) = rxt(k,187)*y(k,60)
         mat(k,830) = rxt(k,187)*y(k,56) + rxt(k,189)*y(k,133) + rxt(k,190)*y(k,212)
         mat(k,653) = rxt(k,234)*y(k,89)
         mat(k,1266) = rxt(k,234)*y(k,73) + rxt(k,171)*y(k,212)
         mat(k,419) = .500_r8*rxt(k,351)*y(k,212)
         mat(k,1770) = mat(k,1770) + rxt(k,157)*y(k,133) + rxt(k,156)*y(k,134)
         mat(k,1884) = mat(k,1884) + rxt(k,220)*y(k,20) + rxt(k,189)*y(k,60) &
                      + rxt(k,157)*y(k,125)
         mat(k,1830) = rxt(k,156)*y(k,125)
         mat(k,358) = rxt(k,303)*y(k,212)
         mat(k,1619) = mat(k,1619) + rxt(k,190)*y(k,60) + rxt(k,171)*y(k,89) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,139)
         mat(k,732) = -(rxt(k,317)*y(k,212))
         mat(k,1581) = -rxt(k,317)*y(k,127)
         mat(k,840) = rxt(k,307)*y(k,126)
         mat(k,424) = .500_r8*rxt(k,377)*y(k,212)
         mat(k,310) = rxt(k,384)*y(k,212)
         mat(k,279) = rxt(k,388)*y(k,212)
         mat(k,871) = rxt(k,389)*y(k,212)
         mat(k,1413) = rxt(k,307)*y(k,29)
         mat(k,1581) = mat(k,1581) + .500_r8*rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) &
                      + rxt(k,388)*y(k,115) + rxt(k,389)*y(k,116)
         mat(k,296) = -(rxt(k,449)*y(k,212))
         mat(k,1530) = -rxt(k,449)*y(k,128)
         mat(k,1321) = rxt(k,446)*y(k,210)
         mat(k,962) = rxt(k,446)*y(k,201)
         mat(k,1893) = -(rxt(k,129)*y(k,134) + 4._r8*rxt(k,130)*y(k,133) + rxt(k,132) &
                      *y(k,77) + rxt(k,133)*y(k,79) + rxt(k,138)*y(k,201) + rxt(k,144) &
                      *y(k,212) + (rxt(k,155) + rxt(k,157)) * y(k,125) + rxt(k,160) &
                      *y(k,126) + rxt(k,165)*y(k,124) + rxt(k,189)*y(k,60) + rxt(k,191) &
                      *y(k,59) + rxt(k,194)*y(k,85) + rxt(k,197)*y(k,92) + rxt(k,220) &
                      *y(k,20) + rxt(k,221)*y(k,19) + rxt(k,223)*y(k,81) + rxt(k,225) &
                      *y(k,91) + rxt(k,256)*y(k,42) + rxt(k,462)*y(k,137))
         mat(k,1839) = -rxt(k,129)*y(k,133)
         mat(k,1016) = -rxt(k,132)*y(k,133)
         mat(k,481) = -rxt(k,133)*y(k,133)
         mat(k,1400) = -rxt(k,138)*y(k,133)
         mat(k,1628) = -rxt(k,144)*y(k,133)
         mat(k,1779) = -(rxt(k,155) + rxt(k,157)) * y(k,133)
         mat(k,1457) = -rxt(k,160)*y(k,133)
         mat(k,1994) = -rxt(k,165)*y(k,133)
         mat(k,835) = -rxt(k,189)*y(k,133)
         mat(k,1704) = -rxt(k,191)*y(k,133)
         mat(k,1307) = -rxt(k,194)*y(k,133)
         mat(k,682) = -rxt(k,197)*y(k,133)
         mat(k,450) = -rxt(k,220)*y(k,133)
         mat(k,1863) = -rxt(k,221)*y(k,133)
         mat(k,701) = -rxt(k,223)*y(k,133)
         mat(k,646) = -rxt(k,225)*y(k,133)
         mat(k,1916) = -rxt(k,256)*y(k,133)
         mat(k,257) = -rxt(k,462)*y(k,133)
         mat(k,1286) = rxt(k,136)*y(k,201)
         mat(k,270) = rxt(k,150)*y(k,124) + rxt(k,151)*y(k,125)
         mat(k,1994) = mat(k,1994) + rxt(k,150)*y(k,112)
         mat(k,1779) = mat(k,1779) + rxt(k,151)*y(k,112)
         mat(k,1400) = mat(k,1400) + rxt(k,136)*y(k,76)
         mat(k,1628) = mat(k,1628) + 2.000_r8*rxt(k,146)*y(k,212)
         mat(k,1837) = -(rxt(k,128)*y(k,211) + rxt(k,129)*y(k,133) + rxt(k,139) &
                      *y(k,201) + rxt(k,140)*y(k,76) + rxt(k,145)*y(k,212) + rxt(k,156) &
                      *y(k,125) + rxt(k,164)*y(k,124) + rxt(k,180)*y(k,56) + rxt(k,212) &
                      *y(k,17) + rxt(k,279)*y(k,25) + rxt(k,308)*y(k,29) + rxt(k,338) &
                      *y(k,105) + rxt(k,352)*y(k,111) + rxt(k,385)*y(k,98) + rxt(k,423) &
                      *y(k,141) + rxt(k,440)*y(k,6) + rxt(k,443)*y(k,110) + rxt(k,465) &
                      *y(k,148) + rxt(k,471)*y(k,150))
         mat(k,1479) = -rxt(k,128)*y(k,134)
         mat(k,1891) = -rxt(k,129)*y(k,134)
         mat(k,1398) = -rxt(k,139)*y(k,134)
         mat(k,1285) = -rxt(k,140)*y(k,134)
         mat(k,1626) = -rxt(k,145)*y(k,134)
         mat(k,1777) = -rxt(k,156)*y(k,134)
         mat(k,1992) = -rxt(k,164)*y(k,134)
         mat(k,1736) = -rxt(k,180)*y(k,134)
         mat(k,1255) = -rxt(k,212)*y(k,134)
         mat(k,457) = -rxt(k,279)*y(k,134)
         mat(k,852) = -rxt(k,308)*y(k,134)
         mat(k,1068) = -rxt(k,338)*y(k,134)
         mat(k,1140) = -rxt(k,352)*y(k,134)
         mat(k,722) = -rxt(k,385)*y(k,134)
         mat(k,366) = -rxt(k,423)*y(k,134)
         mat(k,775) = -rxt(k,440)*y(k,134)
         mat(k,801) = -rxt(k,443)*y(k,134)
         mat(k,401) = -rxt(k,465)*y(k,134)
         mat(k,1091) = -rxt(k,471)*y(k,134)
         mat(k,1244) = .150_r8*rxt(k,293)*y(k,201)
         mat(k,1398) = mat(k,1398) + .150_r8*rxt(k,293)*y(k,195) + .150_r8*rxt(k,343) &
                      *y(k,207)
         mat(k,1214) = .150_r8*rxt(k,343)*y(k,201)
         mat(k,229) = -(rxt(k,472)*y(k,150))
         mat(k,1080) = -rxt(k,472)*y(k,136)
         mat(k,1844) = rxt(k,214)*y(k,59)
         mat(k,1685) = rxt(k,214)*y(k,19) + 2.000_r8*rxt(k,184)*y(k,59)
         mat(k,250) = -(rxt(k,462)*y(k,133) + rxt(k,463)*y(k,212))
         mat(k,1867) = -rxt(k,462)*y(k,137)
         mat(k,1524) = -rxt(k,463)*y(k,137)
         mat(k,895) = rxt(k,331)*y(k,212)
         mat(k,1920) = .100_r8*rxt(k,452)*y(k,216)
         mat(k,1509) = rxt(k,331)*y(k,93)
         mat(k,943) = .100_r8*rxt(k,452)*y(k,124)
         mat(k,356) = -(rxt(k,303)*y(k,212))
         mat(k,1539) = -rxt(k,303)*y(k,139)
         mat(k,1746) = rxt(k,305)*y(k,195)
         mat(k,1219) = rxt(k,305)*y(k,125)
         mat(k,1742) = rxt(k,425)*y(k,188)
         mat(k,408) = rxt(k,425)*y(k,125)
         mat(k,363) = -(rxt(k,422)*y(k,125) + rxt(k,423)*y(k,134))
         mat(k,1747) = -rxt(k,422)*y(k,141)
         mat(k,1790) = -rxt(k,423)*y(k,141)
         mat(k,118) = .070_r8*rxt(k,409)*y(k,212)
         mat(k,1931) = rxt(k,407)*y(k,194)
         mat(k,93) = .060_r8*rxt(k,421)*y(k,212)
         mat(k,143) = .070_r8*rxt(k,437)*y(k,212)
         mat(k,514) = rxt(k,407)*y(k,124)
         mat(k,1540) = .070_r8*rxt(k,409)*y(k,66) + .060_r8*rxt(k,421)*y(k,142) &
                      + .070_r8*rxt(k,437)*y(k,184)
         mat(k,91) = -(rxt(k,421)*y(k,212))
         mat(k,1499) = -rxt(k,421)*y(k,142)
         mat(k,83) = .530_r8*rxt(k,398)*y(k,212)
         mat(k,1499) = mat(k,1499) + .530_r8*rxt(k,398)*y(k,7)
         mat(k,234) = -(rxt(k,424)*y(k,212))
         mat(k,1521) = -rxt(k,424)*y(k,143)
         mat(k,1316) = rxt(k,419)*y(k,213)
         mat(k,346) = rxt(k,419)*y(k,201)
         mat(k,435) = -(rxt(k,320)*y(k,212))
         mat(k,1551) = -rxt(k,320)*y(k,146)
         mat(k,1338) = rxt(k,318)*y(k,214)
         mat(k,629) = rxt(k,318)*y(k,201)
         mat(k,302) = -(rxt(k,324)*y(k,212))
         mat(k,1531) = -rxt(k,324)*y(k,147)
         mat(k,1322) = .850_r8*rxt(k,322)*y(k,215)
         mat(k,982) = .850_r8*rxt(k,322)*y(k,201)
         mat(k,397) = -(rxt(k,465)*y(k,134) + rxt(k,468)*y(k,212))
         mat(k,1791) = -rxt(k,465)*y(k,148)
         mat(k,1545) = -rxt(k,468)*y(k,148)
         mat(k,1083) = -(rxt(k,466)*y(k,19) + rxt(k,467)*y(k,59) + rxt(k,469)*y(k,125) &
                      + rxt(k,471)*y(k,134) + rxt(k,472)*y(k,136) + rxt(k,473) &
                      *y(k,212))
         mat(k,1848) = -rxt(k,466)*y(k,150)
         mat(k,1689) = -rxt(k,467)*y(k,150)
         mat(k,1762) = -rxt(k,469)*y(k,150)
         mat(k,1818) = -rxt(k,471)*y(k,150)
         mat(k,231) = -rxt(k,472)*y(k,150)
         mat(k,1607) = -rxt(k,473)*y(k,150)
         mat(k,1878) = rxt(k,462)*y(k,137)
         mat(k,1818) = mat(k,1818) + rxt(k,465)*y(k,148)
         mat(k,254) = rxt(k,462)*y(k,133)
         mat(k,398) = rxt(k,465)*y(k,134) + rxt(k,468)*y(k,212)
         mat(k,1607) = mat(k,1607) + rxt(k,468)*y(k,148)
         mat(k,726) = -(rxt(k,474)*y(k,212))
         mat(k,1580) = -rxt(k,474)*y(k,151)
         mat(k,1847) = rxt(k,466)*y(k,150)
         mat(k,1687) = rxt(k,467)*y(k,150)
         mat(k,194) = rxt(k,476)*y(k,126) + (rxt(k,477)+.500_r8*rxt(k,479))*y(k,212)
         mat(k,1755) = rxt(k,469)*y(k,150)
         mat(k,1412) = rxt(k,476)*y(k,67)
         mat(k,1798) = rxt(k,471)*y(k,150)
         mat(k,230) = rxt(k,472)*y(k,150)
         mat(k,252) = rxt(k,463)*y(k,212)
         mat(k,1082) = rxt(k,466)*y(k,19) + rxt(k,467)*y(k,59) + rxt(k,469)*y(k,125) &
                      + rxt(k,471)*y(k,134) + rxt(k,472)*y(k,136) + rxt(k,473) &
                      *y(k,212)
         mat(k,1580) = mat(k,1580) + (rxt(k,477)+.500_r8*rxt(k,479))*y(k,67) &
                      + rxt(k,463)*y(k,137) + rxt(k,473)*y(k,150)
         mat(k,168) = -(rxt(k,475)*y(k,222))
         mat(k,2000) = -rxt(k,475)*y(k,152)
         mat(k,725) = rxt(k,474)*y(k,212)
         mat(k,1512) = rxt(k,474)*y(k,151)
         mat(k,752) = .2202005_r8*rxt(k,495)*y(k,134) + .2202005_r8*rxt(k,496) &
                      *y(k,212)
         mat(k,76) = .0023005_r8*rxt(k,497)*y(k,212)
         mat(k,703) = .0031005_r8*rxt(k,500)*y(k,212)
         mat(k,34) = .2381005_r8*rxt(k,501)*y(k,212)
         mat(k,778) = .0508005_r8*rxt(k,503)*y(k,134) + .0508005_r8*rxt(k,504) &
                      *y(k,212)
         mat(k,1783) = .2202005_r8*rxt(k,495)*y(k,6) + .0508005_r8*rxt(k,503)*y(k,110)
         mat(k,40) = .5931005_r8*rxt(k,505)*y(k,212)
         mat(k,104) = .1364005_r8*rxt(k,506)*y(k,212)
         mat(k,128) = .1677005_r8*rxt(k,507)*y(k,212)
         mat(k,1485) = .2202005_r8*rxt(k,496)*y(k,6) + .0023005_r8*rxt(k,497)*y(k,7) &
                      + .0031005_r8*rxt(k,500)*y(k,98) + .2381005_r8*rxt(k,501) &
                      *y(k,104) + .0508005_r8*rxt(k,504)*y(k,110) &
                      + .5931005_r8*rxt(k,505)*y(k,172) + .1364005_r8*rxt(k,506) &
                      *y(k,180) + .1677005_r8*rxt(k,507)*y(k,182)
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
         mat(k,753) = .2067005_r8*rxt(k,495)*y(k,134) + .2067005_r8*rxt(k,496) &
                      *y(k,212)
         mat(k,77) = .0008005_r8*rxt(k,497)*y(k,212)
         mat(k,704) = .0035005_r8*rxt(k,500)*y(k,212)
         mat(k,35) = .1308005_r8*rxt(k,501)*y(k,212)
         mat(k,779) = .1149005_r8*rxt(k,503)*y(k,134) + .1149005_r8*rxt(k,504) &
                      *y(k,212)
         mat(k,1784) = .2067005_r8*rxt(k,495)*y(k,6) + .1149005_r8*rxt(k,503)*y(k,110)
         mat(k,41) = .1534005_r8*rxt(k,505)*y(k,212)
         mat(k,105) = .0101005_r8*rxt(k,506)*y(k,212)
         mat(k,129) = .0174005_r8*rxt(k,507)*y(k,212)
         mat(k,1486) = .2067005_r8*rxt(k,496)*y(k,6) + .0008005_r8*rxt(k,497)*y(k,7) &
                      + .0035005_r8*rxt(k,500)*y(k,98) + .1308005_r8*rxt(k,501) &
                      *y(k,104) + .1149005_r8*rxt(k,504)*y(k,110) &
                      + .1534005_r8*rxt(k,505)*y(k,172) + .0101005_r8*rxt(k,506) &
                      *y(k,180) + .0174005_r8*rxt(k,507)*y(k,182)
         mat(k,754) = .0653005_r8*rxt(k,495)*y(k,134) + .0653005_r8*rxt(k,496) &
                      *y(k,212)
         mat(k,78) = .0843005_r8*rxt(k,497)*y(k,212)
         mat(k,705) = .0003005_r8*rxt(k,500)*y(k,212)
         mat(k,36) = .0348005_r8*rxt(k,501)*y(k,212)
         mat(k,780) = .0348005_r8*rxt(k,503)*y(k,134) + .0348005_r8*rxt(k,504) &
                      *y(k,212)
         mat(k,1785) = .0653005_r8*rxt(k,495)*y(k,6) + .0348005_r8*rxt(k,503)*y(k,110)
         mat(k,42) = .0459005_r8*rxt(k,505)*y(k,212)
         mat(k,106) = .0763005_r8*rxt(k,506)*y(k,212)
         mat(k,130) = .086_r8*rxt(k,507)*y(k,212)
         mat(k,1487) = .0653005_r8*rxt(k,496)*y(k,6) + .0843005_r8*rxt(k,497)*y(k,7) &
                      + .0003005_r8*rxt(k,500)*y(k,98) + .0348005_r8*rxt(k,501) &
                      *y(k,104) + .0348005_r8*rxt(k,504)*y(k,110) &
                      + .0459005_r8*rxt(k,505)*y(k,172) + .0763005_r8*rxt(k,506) &
                      *y(k,180) + .086_r8*rxt(k,507)*y(k,182)
         mat(k,755) = .1749305_r8*rxt(k,494)*y(k,126) + .1284005_r8*rxt(k,495) &
                      *y(k,134) + .1284005_r8*rxt(k,496)*y(k,212)
         mat(k,79) = .0443005_r8*rxt(k,497)*y(k,212)
         mat(k,706) = .0590245_r8*rxt(k,498)*y(k,126) + .0033005_r8*rxt(k,499) &
                      *y(k,134) + .0271005_r8*rxt(k,500)*y(k,212)
         mat(k,37) = .0076005_r8*rxt(k,501)*y(k,212)
         mat(k,781) = .1749305_r8*rxt(k,502)*y(k,126) + .0554005_r8*rxt(k,503) &
                      *y(k,134) + .0554005_r8*rxt(k,504)*y(k,212)
         mat(k,1404) = .1749305_r8*rxt(k,494)*y(k,6) + .0590245_r8*rxt(k,498)*y(k,98) &
                      + .1749305_r8*rxt(k,502)*y(k,110)
         mat(k,1786) = .1284005_r8*rxt(k,495)*y(k,6) + .0033005_r8*rxt(k,499)*y(k,98) &
                      + .0554005_r8*rxt(k,503)*y(k,110)
         mat(k,43) = .0085005_r8*rxt(k,505)*y(k,212)
         mat(k,107) = .2157005_r8*rxt(k,506)*y(k,212)
         mat(k,131) = .0512005_r8*rxt(k,507)*y(k,212)
         mat(k,1488) = .1284005_r8*rxt(k,496)*y(k,6) + .0443005_r8*rxt(k,497)*y(k,7) &
                      + .0271005_r8*rxt(k,500)*y(k,98) + .0076005_r8*rxt(k,501) &
                      *y(k,104) + .0554005_r8*rxt(k,504)*y(k,110) &
                      + .0085005_r8*rxt(k,505)*y(k,172) + .2157005_r8*rxt(k,506) &
                      *y(k,180) + .0512005_r8*rxt(k,507)*y(k,182)
         mat(k,756) = .5901905_r8*rxt(k,494)*y(k,126) + .114_r8*rxt(k,495)*y(k,134) &
                      + .114_r8*rxt(k,496)*y(k,212)
         mat(k,80) = .1621005_r8*rxt(k,497)*y(k,212)
         mat(k,707) = .0250245_r8*rxt(k,498)*y(k,126) + .0474005_r8*rxt(k,500) &
                      *y(k,212)
         mat(k,38) = .0113005_r8*rxt(k,501)*y(k,212)
         mat(k,782) = .5901905_r8*rxt(k,502)*y(k,126) + .1278005_r8*rxt(k,503) &
                      *y(k,134) + .1278005_r8*rxt(k,504)*y(k,212)
         mat(k,1405) = .5901905_r8*rxt(k,494)*y(k,6) + .0250245_r8*rxt(k,498)*y(k,98) &
                      + .5901905_r8*rxt(k,502)*y(k,110)
         mat(k,1787) = .114_r8*rxt(k,495)*y(k,6) + .1278005_r8*rxt(k,503)*y(k,110)
         mat(k,44) = .0128005_r8*rxt(k,505)*y(k,212)
         mat(k,108) = .0232005_r8*rxt(k,506)*y(k,212)
         mat(k,132) = .1598005_r8*rxt(k,507)*y(k,212)
         mat(k,1489) = .114_r8*rxt(k,496)*y(k,6) + .1621005_r8*rxt(k,497)*y(k,7) &
                      + .0474005_r8*rxt(k,500)*y(k,98) + .0113005_r8*rxt(k,501) &
                      *y(k,104) + .1278005_r8*rxt(k,504)*y(k,110) &
                      + .0128005_r8*rxt(k,505)*y(k,172) + .0232005_r8*rxt(k,506) &
                      *y(k,180) + .1598005_r8*rxt(k,507)*y(k,182)
         mat(k,45) = -(rxt(k,505)*y(k,212))
         mat(k,1491) = -rxt(k,505)*y(k,172)
         mat(k,111) = .100_r8*rxt(k,429)*y(k,212)
         mat(k,133) = .230_r8*rxt(k,431)*y(k,212)
         mat(k,1504) = .100_r8*rxt(k,429)*y(k,180) + .230_r8*rxt(k,431)*y(k,182)
         mat(k,483) = -(rxt(k,453)*y(k,212))
         mat(k,1556) = -rxt(k,453)*y(k,174)
         mat(k,1340) = rxt(k,451)*y(k,216)
         mat(k,944) = rxt(k,451)*y(k,201)
         mat(k,507) = -(rxt(k,454)*y(k,212))
         mat(k,1559) = -rxt(k,454)*y(k,175)
         mat(k,1940) = .200_r8*rxt(k,447)*y(k,210) + .200_r8*rxt(k,457)*y(k,217)
         mat(k,1637) = .500_r8*rxt(k,445)*y(k,210)
         mat(k,963) = .200_r8*rxt(k,447)*y(k,124) + .500_r8*rxt(k,445)*y(k,196)
         mat(k,922) = .200_r8*rxt(k,457)*y(k,124)
         mat(k,367) = -(rxt(k,458)*y(k,212))
         mat(k,1541) = -rxt(k,458)*y(k,176)
         mat(k,1332) = rxt(k,456)*y(k,217)
         mat(k,921) = rxt(k,456)*y(k,201)
         mat(k,856) = -(rxt(k,459)*y(k,126) + rxt(k,460)*y(k,212))
         mat(k,1420) = -rxt(k,459)*y(k,177)
         mat(k,1590) = -rxt(k,460)*y(k,177)
         mat(k,764) = .330_r8*rxt(k,440)*y(k,134)
         mat(k,790) = .330_r8*rxt(k,443)*y(k,134)
         mat(k,1958) = .800_r8*rxt(k,447)*y(k,210) + .800_r8*rxt(k,457)*y(k,217)
         mat(k,1420) = mat(k,1420) + rxt(k,448)*y(k,210)
         mat(k,1805) = .330_r8*rxt(k,440)*y(k,6) + .330_r8*rxt(k,443)*y(k,110)
         mat(k,508) = rxt(k,454)*y(k,212)
         mat(k,1644) = .500_r8*rxt(k,445)*y(k,210) + rxt(k,455)*y(k,217)
         mat(k,965) = .800_r8*rxt(k,447)*y(k,124) + rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,196)
         mat(k,1590) = mat(k,1590) + rxt(k,454)*y(k,175)
         mat(k,925) = .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,196)
         mat(k,886) = -(rxt(k,461)*y(k,212))
         mat(k,1593) = -rxt(k,461)*y(k,178)
         mat(k,765) = .300_r8*rxt(k,440)*y(k,134)
         mat(k,791) = .300_r8*rxt(k,443)*y(k,134)
         mat(k,1961) = .900_r8*rxt(k,452)*y(k,216)
         mat(k,1807) = .300_r8*rxt(k,440)*y(k,6) + .300_r8*rxt(k,443)*y(k,110)
         mat(k,1647) = rxt(k,450)*y(k,216)
         mat(k,948) = .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,196)
         mat(k,494) = -(rxt(k,428)*y(k,212))
         mat(k,1557) = -rxt(k,428)*y(k,179)
         mat(k,1341) = rxt(k,426)*y(k,218)
         mat(k,591) = rxt(k,426)*y(k,201)
         mat(k,109) = -(rxt(k,429)*y(k,212))
         mat(k,1502) = -rxt(k,429)*y(k,180)
         mat(k,125) = -(rxt(k,395)*y(k,212))
         mat(k,1505) = -rxt(k,395)*y(k,181)
         mat(k,1311) = rxt(k,392)*y(k,219)
         mat(k,1018) = rxt(k,392)*y(k,201)
         mat(k,134) = -(rxt(k,431)*y(k,212))
         mat(k,1506) = -rxt(k,431)*y(k,182)
         mat(k,563) = -(rxt(k,434)*y(k,212))
         mat(k,1565) = -rxt(k,434)*y(k,183)
         mat(k,1347) = rxt(k,432)*y(k,220)
         mat(k,608) = rxt(k,432)*y(k,201)
         mat(k,142) = -(rxt(k,437)*y(k,212))
         mat(k,1507) = -rxt(k,437)*y(k,184)
         mat(k,135) = .150_r8*rxt(k,431)*y(k,212)
         mat(k,1507) = mat(k,1507) + .150_r8*rxt(k,431)*y(k,182)
         mat(k,326) = -(rxt(k,438)*y(k,212))
         mat(k,1535) = -rxt(k,438)*y(k,185)
         mat(k,1326) = rxt(k,435)*y(k,221)
         mat(k,383) = rxt(k,435)*y(k,201)
         mat(k,409) = -(rxt(k,396)*y(k,201) + rxt(k,397)*y(k,124) + rxt(k,425) &
                      *y(k,125))
         mat(k,1336) = -rxt(k,396)*y(k,188)
         mat(k,1935) = -rxt(k,397)*y(k,188)
         mat(k,1748) = -rxt(k,425)*y(k,188)
         mat(k,165) = rxt(k,402)*y(k,212)
         mat(k,1547) = rxt(k,402)*y(k,22)
         mat(k,809) = -(rxt(k,357)*y(k,201) + (rxt(k,358) + rxt(k,359)) * y(k,124))
         mat(k,1363) = -rxt(k,357)*y(k,189)
         mat(k,1956) = -(rxt(k,358) + rxt(k,359)) * y(k,189)
         mat(k,525) = rxt(k,360)*y(k,212)
         mat(k,159) = rxt(k,361)*y(k,212)
         mat(k,1586) = rxt(k,360)*y(k,2) + rxt(k,361)*y(k,15)
         mat(k,376) = -(rxt(k,399)*y(k,201) + rxt(k,400)*y(k,124))
         mat(k,1333) = -rxt(k,399)*y(k,190)
         mat(k,1932) = -rxt(k,400)*y(k,190)
         mat(k,84) = .350_r8*rxt(k,398)*y(k,212)
         mat(k,286) = rxt(k,401)*y(k,212)
         mat(k,1542) = .350_r8*rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8)
         mat(k,334) = -(rxt(k,403)*y(k,201) + rxt(k,405)*y(k,124))
         mat(k,1327) = -rxt(k,403)*y(k,191)
         mat(k,1927) = -rxt(k,405)*y(k,191)
         mat(k,241) = rxt(k,404)*y(k,212)
         mat(k,112) = .070_r8*rxt(k,429)*y(k,212)
         mat(k,136) = .060_r8*rxt(k,431)*y(k,212)
         mat(k,1536) = rxt(k,404)*y(k,23) + .070_r8*rxt(k,429)*y(k,180) &
                      + .060_r8*rxt(k,431)*y(k,182)
         mat(k,687) = -(4._r8*rxt(k,280)*y(k,192) + rxt(k,281)*y(k,196) + rxt(k,282) &
                      *y(k,201) + rxt(k,283)*y(k,124))
         mat(k,1640) = -rxt(k,281)*y(k,192)
         mat(k,1358) = -rxt(k,282)*y(k,192)
         mat(k,1952) = -rxt(k,283)*y(k,192)
         mat(k,246) = .500_r8*rxt(k,285)*y(k,212)
         mat(k,206) = rxt(k,286)*y(k,56) + rxt(k,287)*y(k,212)
         mat(k,1716) = rxt(k,286)*y(k,28)
         mat(k,1577) = .500_r8*rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28)
         mat(k,666) = -(rxt(k,309)*y(k,196) + rxt(k,310)*y(k,201) + rxt(k,311) &
                      *y(k,124))
         mat(k,1639) = -rxt(k,309)*y(k,193)
         mat(k,1356) = -rxt(k,310)*y(k,193)
         mat(k,1951) = -rxt(k,311)*y(k,193)
         mat(k,315) = rxt(k,312)*y(k,212)
         mat(k,56) = rxt(k,313)*y(k,212)
         mat(k,1575) = rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31)
         mat(k,515) = -(rxt(k,406)*y(k,201) + rxt(k,407)*y(k,124))
         mat(k,1343) = -rxt(k,406)*y(k,194)
         mat(k,1941) = -rxt(k,407)*y(k,194)
         mat(k,178) = rxt(k,408)*y(k,212)
         mat(k,1941) = mat(k,1941) + rxt(k,397)*y(k,188)
         mat(k,1794) = rxt(k,423)*y(k,141)
         mat(k,364) = rxt(k,423)*y(k,134)
         mat(k,410) = rxt(k,397)*y(k,124) + .400_r8*rxt(k,396)*y(k,201)
         mat(k,1343) = mat(k,1343) + .400_r8*rxt(k,396)*y(k,188)
         mat(k,1560) = rxt(k,408)*y(k,32)
         mat(k,1236) = -(4._r8*rxt(k,291)*y(k,195) + rxt(k,292)*y(k,196) + rxt(k,293) &
                      *y(k,201) + rxt(k,294)*y(k,124) + rxt(k,305)*y(k,125) + rxt(k,332) &
                      *y(k,205) + rxt(k,365)*y(k,203) + rxt(k,370)*y(k,204) + rxt(k,379) &
                      *y(k,101) + rxt(k,390)*y(k,219))
         mat(k,1664) = -rxt(k,292)*y(k,195)
         mat(k,1385) = -rxt(k,293)*y(k,195)
         mat(k,1979) = -rxt(k,294)*y(k,195)
         mat(k,1764) = -rxt(k,305)*y(k,195)
         mat(k,1155) = -rxt(k,332)*y(k,195)
         mat(k,1112) = -rxt(k,365)*y(k,195)
         mat(k,1188) = -rxt(k,370)*y(k,195)
         mat(k,1048) = -rxt(k,379)*y(k,195)
         mat(k,1026) = -rxt(k,390)*y(k,195)
         mat(k,771) = .060_r8*rxt(k,440)*y(k,134)
         mat(k,997) = rxt(k,288)*y(k,126) + rxt(k,289)*y(k,212)
         mat(k,1073) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,212)
         mat(k,392) = .500_r8*rxt(k,296)*y(k,212)
         mat(k,717) = .080_r8*rxt(k,385)*y(k,134)
         mat(k,1064) = .100_r8*rxt(k,338)*y(k,134)
         mat(k,797) = .060_r8*rxt(k,443)*y(k,134)
         mat(k,1132) = .280_r8*rxt(k,352)*y(k,134)
         mat(k,1979) = mat(k,1979) + .530_r8*rxt(k,336)*y(k,205) + rxt(k,345)*y(k,207) &
                      + rxt(k,348)*y(k,209) + rxt(k,323)*y(k,215)
         mat(k,1442) = rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) + .530_r8*rxt(k,335) &
                      *y(k,205) + rxt(k,346)*y(k,207)
         mat(k,1824) = .060_r8*rxt(k,440)*y(k,6) + .080_r8*rxt(k,385)*y(k,98) &
                      + .100_r8*rxt(k,338)*y(k,105) + .060_r8*rxt(k,443)*y(k,110) &
                      + .280_r8*rxt(k,352)*y(k,111)
         mat(k,889) = .650_r8*rxt(k,461)*y(k,212)
         mat(k,1236) = mat(k,1236) + .530_r8*rxt(k,332)*y(k,205)
         mat(k,1664) = mat(k,1664) + .260_r8*rxt(k,333)*y(k,205) + rxt(k,342)*y(k,207) &
                      + .300_r8*rxt(k,321)*y(k,215)
         mat(k,1385) = mat(k,1385) + .450_r8*rxt(k,343)*y(k,207) + .200_r8*rxt(k,347) &
                      *y(k,209) + .150_r8*rxt(k,322)*y(k,215)
         mat(k,1155) = mat(k,1155) + .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335) &
                      *y(k,126) + .530_r8*rxt(k,332)*y(k,195) + .260_r8*rxt(k,333) &
                      *y(k,196)
         mat(k,1206) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,342)*y(k,196) &
                      + .450_r8*rxt(k,343)*y(k,201) + 4.000_r8*rxt(k,344)*y(k,207)
         mat(k,546) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,201)
         mat(k,1613) = rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49) + .500_r8*rxt(k,296) &
                      *y(k,51) + .650_r8*rxt(k,461)*y(k,178)
         mat(k,987) = rxt(k,323)*y(k,124) + .300_r8*rxt(k,321)*y(k,196) &
                      + .150_r8*rxt(k,322)*y(k,201)
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
         mat(k,1672) = -(rxt(k,181)*y(k,59) + (4._r8*rxt(k,258) + 4._r8*rxt(k,259) &
                      ) * y(k,196) + rxt(k,260)*y(k,201) + rxt(k,261)*y(k,124) &
                      + rxt(k,281)*y(k,192) + rxt(k,292)*y(k,195) + rxt(k,309) &
                      *y(k,193) + rxt(k,321)*y(k,215) + rxt(k,333)*y(k,205) + rxt(k,342) &
                      *y(k,207) + rxt(k,366)*y(k,203) + rxt(k,371)*y(k,204) + rxt(k,380) &
                      *y(k,101) + rxt(k,391)*y(k,219) + rxt(k,445)*y(k,210) + rxt(k,450) &
                      *y(k,216) + rxt(k,455)*y(k,217))
         mat(k,1698) = -rxt(k,181)*y(k,196)
         mat(k,1394) = -rxt(k,260)*y(k,196)
         mat(k,1988) = -rxt(k,261)*y(k,196)
         mat(k,691) = -rxt(k,281)*y(k,196)
         mat(k,1242) = -rxt(k,292)*y(k,196)
         mat(k,671) = -rxt(k,309)*y(k,196)
         mat(k,990) = -rxt(k,321)*y(k,196)
         mat(k,1161) = -rxt(k,333)*y(k,196)
         mat(k,1212) = -rxt(k,342)*y(k,196)
         mat(k,1118) = -rxt(k,366)*y(k,196)
         mat(k,1194) = -rxt(k,371)*y(k,196)
         mat(k,1054) = -rxt(k,380)*y(k,196)
         mat(k,1031) = -rxt(k,391)*y(k,196)
         mat(k,977) = -rxt(k,445)*y(k,196)
         mat(k,957) = -rxt(k,450)*y(k,196)
         mat(k,938) = -rxt(k,455)*y(k,196)
         mat(k,850) = .280_r8*rxt(k,308)*y(k,134)
         mat(k,433) = rxt(k,295)*y(k,212)
         mat(k,323) = .700_r8*rxt(k,263)*y(k,212)
         mat(k,721) = .050_r8*rxt(k,385)*y(k,134)
         mat(k,1054) = mat(k,1054) + rxt(k,379)*y(k,195)
         mat(k,1988) = mat(k,1988) + rxt(k,294)*y(k,195) + .830_r8*rxt(k,411)*y(k,197) &
                      + .170_r8*rxt(k,417)*y(k,208)
         mat(k,1833) = .280_r8*rxt(k,308)*y(k,29) + .050_r8*rxt(k,385)*y(k,98)
         mat(k,1242) = mat(k,1242) + rxt(k,379)*y(k,101) + rxt(k,294)*y(k,124) &
                      + 4.000_r8*rxt(k,291)*y(k,195) + .900_r8*rxt(k,292)*y(k,196) &
                      + .450_r8*rxt(k,293)*y(k,201) + rxt(k,365)*y(k,203) + rxt(k,370) &
                      *y(k,204) + rxt(k,332)*y(k,205) + rxt(k,341)*y(k,207) &
                      + rxt(k,390)*y(k,219)
         mat(k,1672) = mat(k,1672) + .900_r8*rxt(k,292)*y(k,195)
         mat(k,626) = .830_r8*rxt(k,411)*y(k,124) + .330_r8*rxt(k,410)*y(k,201)
         mat(k,1394) = mat(k,1394) + .450_r8*rxt(k,293)*y(k,195) + .330_r8*rxt(k,410) &
                      *y(k,197) + .070_r8*rxt(k,416)*y(k,208)
         mat(k,1118) = mat(k,1118) + rxt(k,365)*y(k,195)
         mat(k,1194) = mat(k,1194) + rxt(k,370)*y(k,195)
         mat(k,1161) = mat(k,1161) + rxt(k,332)*y(k,195)
         mat(k,1212) = mat(k,1212) + rxt(k,341)*y(k,195)
         mat(k,749) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,201)
         mat(k,1622) = rxt(k,295)*y(k,50) + .700_r8*rxt(k,263)*y(k,53)
         mat(k,1031) = mat(k,1031) + rxt(k,390)*y(k,195)
         mat(k,621) = -(rxt(k,410)*y(k,201) + rxt(k,411)*y(k,124) + rxt(k,412) &
                      *y(k,125))
         mat(k,1352) = -rxt(k,410)*y(k,197)
         mat(k,1948) = -rxt(k,411)*y(k,197)
         mat(k,1753) = -rxt(k,412)*y(k,197)
         mat(k,459) = -((rxt(k,329) + rxt(k,330)) * y(k,124))
         mat(k,1937) = -(rxt(k,329) + rxt(k,330)) * y(k,198)
         mat(k,259) = rxt(k,328)*y(k,212)
         mat(k,1553) = rxt(k,328)*y(k,16)
         mat(k,1922) = .750_r8*rxt(k,298)*y(k,200)
         mat(k,575) = .750_r8*rxt(k,298)*y(k,124)
         mat(k,576) = -(rxt(k,297)*y(k,201) + rxt(k,298)*y(k,124))
         mat(k,1348) = -rxt(k,297)*y(k,200)
         mat(k,1944) = -rxt(k,298)*y(k,200)
         mat(k,452) = rxt(k,304)*y(k,212)
         mat(k,1566) = rxt(k,304)*y(k,25)
         mat(k,1390) = -((rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,76) + rxt(k,138) &
                      *y(k,133) + rxt(k,139)*y(k,134) + rxt(k,143)*y(k,212) &
                      + 4._r8*rxt(k,148)*y(k,201) + rxt(k,158)*y(k,126) + rxt(k,163) &
                      *y(k,124) + rxt(k,168)*y(k,125) + (rxt(k,178) + rxt(k,179) &
                      ) * y(k,56) + rxt(k,185)*y(k,59) + rxt(k,211)*y(k,17) + rxt(k,217) &
                      *y(k,19) + rxt(k,254)*y(k,42) + rxt(k,260)*y(k,196) + rxt(k,268) &
                      *y(k,202) + rxt(k,282)*y(k,192) + rxt(k,293)*y(k,195) + rxt(k,297) &
                      *y(k,200) + rxt(k,310)*y(k,193) + rxt(k,318)*y(k,214) + rxt(k,322) &
                      *y(k,215) + rxt(k,334)*y(k,205) + rxt(k,343)*y(k,207) + rxt(k,347) &
                      *y(k,209) + rxt(k,357)*y(k,189) + rxt(k,367)*y(k,203) + rxt(k,372) &
                      *y(k,204) + rxt(k,381)*y(k,101) + rxt(k,392)*y(k,219) + rxt(k,396) &
                      *y(k,188) + rxt(k,399)*y(k,190) + rxt(k,403)*y(k,191) + rxt(k,406) &
                      *y(k,194) + rxt(k,410)*y(k,197) + rxt(k,413)*y(k,206) + rxt(k,416) &
                      *y(k,208) + rxt(k,419)*y(k,213) + rxt(k,426)*y(k,218) + rxt(k,432) &
                      *y(k,220) + rxt(k,435)*y(k,221) + rxt(k,446)*y(k,210) + rxt(k,451) &
                      *y(k,216) + rxt(k,456)*y(k,217))
         mat(k,1278) = -(rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,201)
         mat(k,1883) = -rxt(k,138)*y(k,201)
         mat(k,1829) = -rxt(k,139)*y(k,201)
         mat(k,1618) = -rxt(k,143)*y(k,201)
         mat(k,1447) = -rxt(k,158)*y(k,201)
         mat(k,1984) = -rxt(k,163)*y(k,201)
         mat(k,1769) = -rxt(k,168)*y(k,201)
         mat(k,1728) = -(rxt(k,178) + rxt(k,179)) * y(k,201)
         mat(k,1694) = -rxt(k,185)*y(k,201)
         mat(k,1252) = -rxt(k,211)*y(k,201)
         mat(k,1853) = -rxt(k,217)*y(k,201)
         mat(k,1906) = -rxt(k,254)*y(k,201)
         mat(k,1668) = -rxt(k,260)*y(k,201)
         mat(k,342) = -rxt(k,268)*y(k,201)
         mat(k,689) = -rxt(k,282)*y(k,201)
         mat(k,1239) = -rxt(k,293)*y(k,201)
         mat(k,578) = -rxt(k,297)*y(k,201)
         mat(k,669) = -rxt(k,310)*y(k,201)
         mat(k,633) = -rxt(k,318)*y(k,201)
         mat(k,988) = -rxt(k,322)*y(k,201)
         mat(k,1158) = -rxt(k,334)*y(k,201)
         mat(k,1209) = -rxt(k,343)*y(k,201)
         mat(k,547) = -rxt(k,347)*y(k,201)
         mat(k,813) = -rxt(k,357)*y(k,201)
         mat(k,1115) = -rxt(k,367)*y(k,201)
         mat(k,1191) = -rxt(k,372)*y(k,201)
         mat(k,1051) = -rxt(k,381)*y(k,201)
         mat(k,1028) = -rxt(k,392)*y(k,201)
         mat(k,411) = -rxt(k,396)*y(k,201)
         mat(k,379) = -rxt(k,399)*y(k,201)
         mat(k,336) = -rxt(k,403)*y(k,201)
         mat(k,516) = -rxt(k,406)*y(k,201)
         mat(k,624) = -rxt(k,410)*y(k,201)
         mat(k,586) = -rxt(k,413)*y(k,201)
         mat(k,747) = -rxt(k,416)*y(k,201)
         mat(k,349) = -rxt(k,419)*y(k,201)
         mat(k,600) = -rxt(k,426)*y(k,201)
         mat(k,617) = -rxt(k,432)*y(k,201)
         mat(k,387) = -rxt(k,435)*y(k,201)
         mat(k,974) = -rxt(k,446)*y(k,201)
         mat(k,955) = -rxt(k,451)*y(k,201)
         mat(k,935) = -rxt(k,456)*y(k,201)
         mat(k,772) = .570_r8*rxt(k,440)*y(k,134)
         mat(k,85) = .650_r8*rxt(k,398)*y(k,212)
         mat(k,1252) = mat(k,1252) + rxt(k,210)*y(k,42)
         mat(k,1853) = mat(k,1853) + rxt(k,222)*y(k,212)
         mat(k,203) = .350_r8*rxt(k,277)*y(k,212)
         mat(k,455) = .130_r8*rxt(k,279)*y(k,134)
         mat(k,174) = rxt(k,284)*y(k,212)
         mat(k,847) = .280_r8*rxt(k,308)*y(k,134)
         mat(k,1906) = mat(k,1906) + rxt(k,210)*y(k,17) + rxt(k,174)*y(k,56) &
                      + rxt(k,255)*y(k,126) + rxt(k,256)*y(k,133)
         mat(k,53) = rxt(k,290)*y(k,212)
         mat(k,662) = rxt(k,262)*y(k,212)
         mat(k,1728) = mat(k,1728) + rxt(k,174)*y(k,42) + rxt(k,177)*y(k,79)
         mat(k,1694) = mat(k,1694) + rxt(k,181)*y(k,196) + rxt(k,192)*y(k,212)
         mat(k,913) = rxt(k,265)*y(k,212)
         mat(k,119) = .730_r8*rxt(k,409)*y(k,212)
         mat(k,196) = .500_r8*rxt(k,479)*y(k,212)
         mat(k,867) = rxt(k,301)*y(k,212)
         mat(k,741) = rxt(k,302)*y(k,212)
         mat(k,478) = rxt(k,177)*y(k,56) + rxt(k,133)*y(k,133) + rxt(k,142)*y(k,212)
         mat(k,97) = rxt(k,266)*y(k,212)
         mat(k,658) = rxt(k,267)*y(k,212)
         mat(k,905) = rxt(k,331)*y(k,212)
         mat(k,918) = rxt(k,316)*y(k,212)
         mat(k,718) = .370_r8*rxt(k,385)*y(k,134)
         mat(k,472) = .300_r8*rxt(k,376)*y(k,212)
         mat(k,429) = rxt(k,377)*y(k,212)
         mat(k,1051) = mat(k,1051) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) &
                      + rxt(k,379)*y(k,195) + 1.200_r8*rxt(k,380)*y(k,196)
         mat(k,311) = rxt(k,384)*y(k,212)
         mat(k,1066) = .140_r8*rxt(k,338)*y(k,134)
         mat(k,217) = .200_r8*rxt(k,340)*y(k,212)
         mat(k,418) = .500_r8*rxt(k,351)*y(k,212)
         mat(k,798) = .570_r8*rxt(k,443)*y(k,134)
         mat(k,1135) = .280_r8*rxt(k,352)*y(k,134)
         mat(k,282) = rxt(k,388)*y(k,212)
         mat(k,879) = rxt(k,389)*y(k,212)
         mat(k,1984) = mat(k,1984) + rxt(k,382)*y(k,101) + rxt(k,358)*y(k,189) &
                      + rxt(k,400)*y(k,190) + rxt(k,405)*y(k,191) + rxt(k,283) &
                      *y(k,192) + rxt(k,311)*y(k,193) + rxt(k,261)*y(k,196) &
                      + .170_r8*rxt(k,411)*y(k,197) + rxt(k,329)*y(k,198) &
                      + .250_r8*rxt(k,298)*y(k,200) + rxt(k,270)*y(k,202) &
                      + .920_r8*rxt(k,368)*y(k,203) + .920_r8*rxt(k,374)*y(k,204) &
                      + .470_r8*rxt(k,336)*y(k,205) + .400_r8*rxt(k,414)*y(k,206) &
                      + .830_r8*rxt(k,417)*y(k,208) + rxt(k,420)*y(k,213) + rxt(k,319) &
                      *y(k,214) + .900_r8*rxt(k,452)*y(k,216) + .800_r8*rxt(k,457) &
                      *y(k,217) + rxt(k,427)*y(k,218) + rxt(k,393)*y(k,219) &
                      + rxt(k,433)*y(k,220) + rxt(k,436)*y(k,221)
         mat(k,1447) = mat(k,1447) + rxt(k,255)*y(k,42) + rxt(k,383)*y(k,101) &
                      + rxt(k,369)*y(k,203) + rxt(k,375)*y(k,204) + .470_r8*rxt(k,335) &
                      *y(k,205) + rxt(k,161)*y(k,212) + rxt(k,394)*y(k,219)
         mat(k,1883) = mat(k,1883) + rxt(k,256)*y(k,42) + rxt(k,133)*y(k,79)
         mat(k,1829) = mat(k,1829) + .570_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .280_r8*rxt(k,308)*y(k,29) + .370_r8*rxt(k,385) &
                      *y(k,98) + .140_r8*rxt(k,338)*y(k,105) + .570_r8*rxt(k,443) &
                      *y(k,110) + .280_r8*rxt(k,352)*y(k,111) + rxt(k,145)*y(k,212)
         mat(k,94) = .800_r8*rxt(k,421)*y(k,212)
         mat(k,728) = rxt(k,474)*y(k,212)
         mat(k,890) = .200_r8*rxt(k,461)*y(k,212)
         mat(k,114) = .280_r8*rxt(k,429)*y(k,212)
         mat(k,140) = .380_r8*rxt(k,431)*y(k,212)
         mat(k,145) = .630_r8*rxt(k,437)*y(k,212)
         mat(k,813) = mat(k,813) + rxt(k,358)*y(k,124)
         mat(k,379) = mat(k,379) + rxt(k,400)*y(k,124)
         mat(k,336) = mat(k,336) + rxt(k,405)*y(k,124)
         mat(k,689) = mat(k,689) + rxt(k,283)*y(k,124) + 2.400_r8*rxt(k,280)*y(k,192) &
                      + rxt(k,281)*y(k,196)
         mat(k,669) = mat(k,669) + rxt(k,311)*y(k,124) + rxt(k,309)*y(k,196)
         mat(k,1239) = mat(k,1239) + rxt(k,379)*y(k,101) + .900_r8*rxt(k,292)*y(k,196) &
                      + rxt(k,365)*y(k,203) + rxt(k,370)*y(k,204) + .470_r8*rxt(k,332) &
                      *y(k,205) + rxt(k,390)*y(k,219)
         mat(k,1668) = mat(k,1668) + rxt(k,181)*y(k,59) + 1.200_r8*rxt(k,380)*y(k,101) &
                      + rxt(k,261)*y(k,124) + rxt(k,281)*y(k,192) + rxt(k,309) &
                      *y(k,193) + .900_r8*rxt(k,292)*y(k,195) + 4.000_r8*rxt(k,258) &
                      *y(k,196) + rxt(k,366)*y(k,203) + rxt(k,371)*y(k,204) &
                      + .730_r8*rxt(k,333)*y(k,205) + rxt(k,342)*y(k,207) &
                      + .500_r8*rxt(k,445)*y(k,210) + .300_r8*rxt(k,321)*y(k,215) &
                      + rxt(k,450)*y(k,216) + rxt(k,455)*y(k,217) + .800_r8*rxt(k,391) &
                      *y(k,219)
         mat(k,624) = mat(k,624) + .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410) &
                      *y(k,201)
         mat(k,463) = rxt(k,329)*y(k,124)
         mat(k,578) = mat(k,578) + .250_r8*rxt(k,298)*y(k,124)
         mat(k,1390) = mat(k,1390) + .070_r8*rxt(k,410)*y(k,197) + .160_r8*rxt(k,413) &
                      *y(k,206) + .330_r8*rxt(k,416)*y(k,208)
         mat(k,342) = mat(k,342) + rxt(k,270)*y(k,124)
         mat(k,1115) = mat(k,1115) + .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) &
                      + rxt(k,365)*y(k,195) + rxt(k,366)*y(k,196)
         mat(k,1191) = mat(k,1191) + .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,195) + rxt(k,371)*y(k,196)
         mat(k,1158) = mat(k,1158) + .470_r8*rxt(k,336)*y(k,124) + .470_r8*rxt(k,335) &
                      *y(k,126) + .470_r8*rxt(k,332)*y(k,195) + .730_r8*rxt(k,333) &
                      *y(k,196)
         mat(k,586) = mat(k,586) + .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413) &
                      *y(k,201)
         mat(k,1209) = mat(k,1209) + rxt(k,342)*y(k,196)
         mat(k,747) = mat(k,747) + .830_r8*rxt(k,417)*y(k,124) + .330_r8*rxt(k,416) &
                      *y(k,201)
         mat(k,974) = mat(k,974) + .500_r8*rxt(k,445)*y(k,196)
         mat(k,1618) = mat(k,1618) + .650_r8*rxt(k,398)*y(k,7) + rxt(k,222)*y(k,19) &
                      + .350_r8*rxt(k,277)*y(k,24) + rxt(k,284)*y(k,26) + rxt(k,290) &
                      *y(k,47) + rxt(k,262)*y(k,52) + rxt(k,192)*y(k,59) + rxt(k,265) &
                      *y(k,62) + .730_r8*rxt(k,409)*y(k,66) + .500_r8*rxt(k,479) &
                      *y(k,67) + rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) + rxt(k,142) &
                      *y(k,79) + rxt(k,266)*y(k,86) + rxt(k,267)*y(k,87) + rxt(k,331) &
                      *y(k,93) + rxt(k,316)*y(k,95) + .300_r8*rxt(k,376)*y(k,99) &
                      + rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) + .200_r8*rxt(k,340) &
                      *y(k,106) + .500_r8*rxt(k,351)*y(k,109) + rxt(k,388)*y(k,115) &
                      + rxt(k,389)*y(k,116) + rxt(k,161)*y(k,126) + rxt(k,145) &
                      *y(k,134) + .800_r8*rxt(k,421)*y(k,142) + rxt(k,474)*y(k,151) &
                      + .200_r8*rxt(k,461)*y(k,178) + .280_r8*rxt(k,429)*y(k,180) &
                      + .380_r8*rxt(k,431)*y(k,182) + .630_r8*rxt(k,437)*y(k,184)
         mat(k,349) = mat(k,349) + rxt(k,420)*y(k,124)
         mat(k,633) = mat(k,633) + rxt(k,319)*y(k,124)
         mat(k,988) = mat(k,988) + .300_r8*rxt(k,321)*y(k,196)
         mat(k,955) = mat(k,955) + .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,196)
         mat(k,935) = mat(k,935) + .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,196)
         mat(k,600) = mat(k,600) + rxt(k,427)*y(k,124)
         mat(k,1028) = mat(k,1028) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126) &
                      + rxt(k,390)*y(k,195) + .800_r8*rxt(k,391)*y(k,196)
         mat(k,617) = mat(k,617) + rxt(k,433)*y(k,124)
         mat(k,387) = mat(k,387) + rxt(k,436)*y(k,124)
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
         mat(k,340) = -(rxt(k,268)*y(k,201) + rxt(k,270)*y(k,124))
         mat(k,1328) = -rxt(k,268)*y(k,202)
         mat(k,1928) = -rxt(k,270)*y(k,202)
         mat(k,1897) = rxt(k,254)*y(k,201)
         mat(k,1328) = mat(k,1328) + rxt(k,254)*y(k,42)
         mat(k,1108) = -(rxt(k,365)*y(k,195) + rxt(k,366)*y(k,196) + rxt(k,367) &
                      *y(k,201) + rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126))
         mat(k,1231) = -rxt(k,365)*y(k,203)
         mat(k,1659) = -rxt(k,366)*y(k,203)
         mat(k,1380) = -rxt(k,367)*y(k,203)
         mat(k,1974) = -rxt(k,368)*y(k,203)
         mat(k,1437) = -rxt(k,369)*y(k,203)
         mat(k,714) = .600_r8*rxt(k,386)*y(k,212)
         mat(k,1608) = .600_r8*rxt(k,386)*y(k,98)
         mat(k,1186) = -(rxt(k,370)*y(k,195) + rxt(k,371)*y(k,196) + rxt(k,372) &
                      *y(k,201) + rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126))
         mat(k,1234) = -rxt(k,370)*y(k,204)
         mat(k,1662) = -rxt(k,371)*y(k,204)
         mat(k,1383) = -rxt(k,372)*y(k,204)
         mat(k,1977) = -rxt(k,374)*y(k,204)
         mat(k,1440) = -rxt(k,375)*y(k,204)
         mat(k,716) = .400_r8*rxt(k,386)*y(k,212)
         mat(k,1611) = .400_r8*rxt(k,386)*y(k,98)
         mat(k,1153) = -(rxt(k,332)*y(k,195) + rxt(k,333)*y(k,196) + rxt(k,334) &
                      *y(k,201) + rxt(k,335)*y(k,126) + (rxt(k,336) + rxt(k,337) &
                      ) * y(k,124))
         mat(k,1233) = -rxt(k,332)*y(k,205)
         mat(k,1661) = -rxt(k,333)*y(k,205)
         mat(k,1382) = -rxt(k,334)*y(k,205)
         mat(k,1439) = -rxt(k,335)*y(k,205)
         mat(k,1976) = -(rxt(k,336) + rxt(k,337)) * y(k,205)
         mat(k,1062) = .500_r8*rxt(k,339)*y(k,212)
         mat(k,215) = .200_r8*rxt(k,340)*y(k,212)
         mat(k,1131) = rxt(k,353)*y(k,212)
         mat(k,1610) = .500_r8*rxt(k,339)*y(k,105) + .200_r8*rxt(k,340)*y(k,106) &
                      + rxt(k,353)*y(k,111)
         mat(k,583) = -(rxt(k,413)*y(k,201) + rxt(k,414)*y(k,124) + rxt(k,415) &
                      *y(k,125))
         mat(k,1349) = -rxt(k,413)*y(k,206)
         mat(k,1945) = -rxt(k,414)*y(k,206)
         mat(k,1752) = -rxt(k,415)*y(k,206)
         mat(k,1205) = -(rxt(k,341)*y(k,195) + rxt(k,342)*y(k,196) + rxt(k,343) &
                      *y(k,201) + 4._r8*rxt(k,344)*y(k,207) + rxt(k,345)*y(k,124) &
                      + rxt(k,346)*y(k,126) + rxt(k,354)*y(k,125))
         mat(k,1235) = -rxt(k,341)*y(k,207)
         mat(k,1663) = -rxt(k,342)*y(k,207)
         mat(k,1384) = -rxt(k,343)*y(k,207)
         mat(k,1978) = -rxt(k,345)*y(k,207)
         mat(k,1441) = -rxt(k,346)*y(k,207)
         mat(k,1763) = -rxt(k,354)*y(k,207)
         mat(k,1063) = .500_r8*rxt(k,339)*y(k,212)
         mat(k,216) = .500_r8*rxt(k,340)*y(k,212)
         mat(k,1612) = .500_r8*rxt(k,339)*y(k,105) + .500_r8*rxt(k,340)*y(k,106)
         mat(k,744) = -(rxt(k,416)*y(k,201) + rxt(k,417)*y(k,124) + rxt(k,418) &
                      *y(k,125))
         mat(k,1362) = -rxt(k,416)*y(k,208)
         mat(k,1955) = -rxt(k,417)*y(k,208)
         mat(k,1757) = -rxt(k,418)*y(k,208)
         mat(k,544) = -(rxt(k,347)*y(k,201) + rxt(k,348)*y(k,124))
         mat(k,1345) = -rxt(k,347)*y(k,209)
         mat(k,1943) = -rxt(k,348)*y(k,209)
         mat(k,404) = rxt(k,349)*y(k,212)
         mat(k,220) = rxt(k,350)*y(k,212)
         mat(k,1563) = rxt(k,349)*y(k,107) + rxt(k,350)*y(k,108)
         mat(k,969) = -(rxt(k,445)*y(k,196) + rxt(k,446)*y(k,201) + rxt(k,447) &
                      *y(k,124) + rxt(k,448)*y(k,126))
         mat(k,1652) = -rxt(k,445)*y(k,210)
         mat(k,1372) = -rxt(k,446)*y(k,210)
         mat(k,1967) = -rxt(k,447)*y(k,210)
         mat(k,1429) = -rxt(k,448)*y(k,210)
         mat(k,768) = rxt(k,439)*y(k,126)
         mat(k,794) = rxt(k,442)*y(k,126)
         mat(k,1429) = mat(k,1429) + rxt(k,439)*y(k,6) + rxt(k,442)*y(k,110) &
                      + .500_r8*rxt(k,459)*y(k,177)
         mat(k,298) = rxt(k,449)*y(k,212)
         mat(k,860) = .500_r8*rxt(k,459)*y(k,126)
         mat(k,1599) = rxt(k,449)*y(k,128)
         mat(k,1473) = -(rxt(k,124)*y(k,77) + rxt(k,125)*y(k,222) + rxt(k,128) &
                      *y(k,134) + (rxt(k,206) + rxt(k,207)) * y(k,85) + (rxt(k,229) &
                      + rxt(k,230)) * y(k,81) + rxt(k,235)*y(k,64) + rxt(k,236) &
                      *y(k,65) + rxt(k,274)*y(k,86))
         mat(k,1012) = -rxt(k,124)*y(k,211)
         mat(k,2011) = -rxt(k,125)*y(k,211)
         mat(k,1831) = -rxt(k,128)*y(k,211)
         mat(k,1299) = -(rxt(k,206) + rxt(k,207)) * y(k,211)
         mat(k,698) = -(rxt(k,229) + rxt(k,230)) * y(k,211)
         mat(k,61) = -rxt(k,235)*y(k,211)
         mat(k,102) = -rxt(k,236)*y(k,211)
         mat(k,98) = -rxt(k,274)*y(k,211)
         mat(k,1621) = -(rxt(k,141)*y(k,77) + rxt(k,142)*y(k,79) + rxt(k,143)*y(k,201) &
                      + rxt(k,144)*y(k,133) + rxt(k,145)*y(k,134) + (4._r8*rxt(k,146) &
                      + 4._r8*rxt(k,147)) * y(k,212) + rxt(k,149)*y(k,90) + rxt(k,161) &
                      *y(k,126) + rxt(k,162)*y(k,112) + rxt(k,170)*y(k,125) + rxt(k,171) &
                      *y(k,89) + rxt(k,190)*y(k,60) + (rxt(k,192) + rxt(k,193) &
                      ) * y(k,59) + rxt(k,195)*y(k,85) + rxt(k,198)*y(k,92) + rxt(k,222) &
                      *y(k,19) + rxt(k,224)*y(k,81) + rxt(k,257)*y(k,42) + rxt(k,262) &
                      *y(k,52) + rxt(k,263)*y(k,53) + (rxt(k,265) + rxt(k,275) &
                      ) * y(k,62) + rxt(k,266)*y(k,86) + rxt(k,267)*y(k,87) + rxt(k,277) &
                      *y(k,24) + rxt(k,284)*y(k,26) + rxt(k,285)*y(k,27) + rxt(k,287) &
                      *y(k,28) + rxt(k,289)*y(k,45) + rxt(k,290)*y(k,47) + rxt(k,295) &
                      *y(k,50) + rxt(k,296)*y(k,51) + rxt(k,301)*y(k,74) + rxt(k,302) &
                      *y(k,75) + rxt(k,303)*y(k,139) + rxt(k,304)*y(k,25) + rxt(k,312) &
                      *y(k,30) + rxt(k,313)*y(k,31) + rxt(k,315)*y(k,49) + rxt(k,316) &
                      *y(k,95) + rxt(k,317)*y(k,127) + rxt(k,320)*y(k,146) + rxt(k,324) &
                      *y(k,147) + rxt(k,325)*y(k,29) + rxt(k,326)*y(k,48) + rxt(k,328) &
                      *y(k,16) + rxt(k,331)*y(k,93) + rxt(k,339)*y(k,105) + rxt(k,340) &
                      *y(k,106) + rxt(k,349)*y(k,107) + rxt(k,350)*y(k,108) + rxt(k,351) &
                      *y(k,109) + rxt(k,353)*y(k,111) + rxt(k,356)*y(k,1) + rxt(k,360) &
                      *y(k,2) + rxt(k,361)*y(k,15) + rxt(k,362)*y(k,94) + rxt(k,363) &
                      *y(k,96) + rxt(k,364)*y(k,97) + rxt(k,376)*y(k,99) + rxt(k,377) &
                      *y(k,100) + rxt(k,384)*y(k,102) + rxt(k,386)*y(k,98) + rxt(k,387) &
                      *y(k,103) + rxt(k,388)*y(k,115) + rxt(k,389)*y(k,116) + rxt(k,395) &
                      *y(k,181) + rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8) + rxt(k,402) &
                      *y(k,22) + rxt(k,404)*y(k,23) + rxt(k,408)*y(k,32) + rxt(k,409) &
                      *y(k,66) + rxt(k,421)*y(k,142) + rxt(k,424)*y(k,143) + rxt(k,428) &
                      *y(k,179) + rxt(k,429)*y(k,180) + rxt(k,431)*y(k,182) + rxt(k,434) &
                      *y(k,183) + rxt(k,437)*y(k,184) + rxt(k,438)*y(k,185) + rxt(k,441) &
                      *y(k,6) + rxt(k,444)*y(k,110) + rxt(k,449)*y(k,128) + rxt(k,453) &
                      *y(k,174) + rxt(k,454)*y(k,175) + rxt(k,458)*y(k,176) + rxt(k,460) &
                      *y(k,177) + rxt(k,461)*y(k,178) + rxt(k,463)*y(k,137) + rxt(k,468) &
                      *y(k,148) + rxt(k,473)*y(k,150) + rxt(k,474)*y(k,151) + (rxt(k,477) &
                      + rxt(k,479)) * y(k,67) + rxt(k,478)*y(k,120))
         mat(k,1013) = -rxt(k,141)*y(k,212)
         mat(k,479) = -rxt(k,142)*y(k,212)
         mat(k,1393) = -rxt(k,143)*y(k,212)
         mat(k,1886) = -rxt(k,144)*y(k,212)
         mat(k,1832) = -rxt(k,145)*y(k,212)
         mat(k,275) = -rxt(k,149)*y(k,212)
         mat(k,1450) = -rxt(k,161)*y(k,212)
         mat(k,268) = -rxt(k,162)*y(k,212)
         mat(k,1772) = -rxt(k,170)*y(k,212)
         mat(k,1268) = -rxt(k,171)*y(k,212)
         mat(k,831) = -rxt(k,190)*y(k,212)
         mat(k,1697) = -(rxt(k,192) + rxt(k,193)) * y(k,212)
         mat(k,1300) = -rxt(k,195)*y(k,212)
         mat(k,679) = -rxt(k,198)*y(k,212)
         mat(k,1856) = -rxt(k,222)*y(k,212)
         mat(k,699) = -rxt(k,224)*y(k,212)
         mat(k,1909) = -rxt(k,257)*y(k,212)
         mat(k,663) = -rxt(k,262)*y(k,212)
         mat(k,322) = -rxt(k,263)*y(k,212)
         mat(k,914) = -(rxt(k,265) + rxt(k,275)) * y(k,212)
         mat(k,99) = -rxt(k,266)*y(k,212)
         mat(k,659) = -rxt(k,267)*y(k,212)
         mat(k,204) = -rxt(k,277)*y(k,212)
         mat(k,175) = -rxt(k,284)*y(k,212)
         mat(k,249) = -rxt(k,285)*y(k,212)
         mat(k,208) = -rxt(k,287)*y(k,212)
         mat(k,1002) = -rxt(k,289)*y(k,212)
         mat(k,54) = -rxt(k,290)*y(k,212)
         mat(k,432) = -rxt(k,295)*y(k,212)
         mat(k,393) = -rxt(k,296)*y(k,212)
         mat(k,868) = -rxt(k,301)*y(k,212)
         mat(k,742) = -rxt(k,302)*y(k,212)
         mat(k,359) = -rxt(k,303)*y(k,212)
         mat(k,456) = -rxt(k,304)*y(k,212)
         mat(k,318) = -rxt(k,312)*y(k,212)
         mat(k,57) = -rxt(k,313)*y(k,212)
         mat(k,1078) = -rxt(k,315)*y(k,212)
         mat(k,919) = -rxt(k,316)*y(k,212)
         mat(k,735) = -rxt(k,317)*y(k,212)
         mat(k,440) = -rxt(k,320)*y(k,212)
         mat(k,305) = -rxt(k,324)*y(k,212)
         mat(k,849) = -rxt(k,325)*y(k,212)
         mat(k,822) = -rxt(k,326)*y(k,212)
         mat(k,263) = -rxt(k,328)*y(k,212)
         mat(k,906) = -rxt(k,331)*y(k,212)
         mat(k,1067) = -rxt(k,339)*y(k,212)
         mat(k,218) = -rxt(k,340)*y(k,212)
         mat(k,407) = -rxt(k,349)*y(k,212)
         mat(k,223) = -rxt(k,350)*y(k,212)
         mat(k,420) = -rxt(k,351)*y(k,212)
         mat(k,1137) = -rxt(k,353)*y(k,212)
         mat(k,540) = -rxt(k,356)*y(k,212)
         mat(k,530) = -rxt(k,360)*y(k,212)
         mat(k,160) = -rxt(k,361)*y(k,212)
         mat(k,154) = -rxt(k,362)*y(k,212)
         mat(k,213) = -rxt(k,363)*y(k,212)
         mat(k,70) = -rxt(k,364)*y(k,212)
         mat(k,473) = -rxt(k,376)*y(k,212)
         mat(k,430) = -rxt(k,377)*y(k,212)
         mat(k,312) = -rxt(k,384)*y(k,212)
         mat(k,720) = -rxt(k,386)*y(k,212)
         mat(k,557) = -rxt(k,387)*y(k,212)
         mat(k,283) = -rxt(k,388)*y(k,212)
         mat(k,880) = -rxt(k,389)*y(k,212)
         mat(k,127) = -rxt(k,395)*y(k,212)
         mat(k,86) = -rxt(k,398)*y(k,212)
         mat(k,289) = -rxt(k,401)*y(k,212)
         mat(k,166) = -rxt(k,402)*y(k,212)
         mat(k,244) = -rxt(k,404)*y(k,212)
         mat(k,179) = -rxt(k,408)*y(k,212)
         mat(k,120) = -rxt(k,409)*y(k,212)
         mat(k,95) = -rxt(k,421)*y(k,212)
         mat(k,238) = -rxt(k,424)*y(k,212)
         mat(k,502) = -rxt(k,428)*y(k,212)
         mat(k,115) = -rxt(k,429)*y(k,212)
         mat(k,141) = -rxt(k,431)*y(k,212)
         mat(k,573) = -rxt(k,434)*y(k,212)
         mat(k,146) = -rxt(k,437)*y(k,212)
         mat(k,331) = -rxt(k,438)*y(k,212)
         mat(k,774) = -rxt(k,441)*y(k,212)
         mat(k,800) = -rxt(k,444)*y(k,212)
         mat(k,300) = -rxt(k,449)*y(k,212)
         mat(k,490) = -rxt(k,453)*y(k,212)
         mat(k,511) = -rxt(k,454)*y(k,212)
         mat(k,372) = -rxt(k,458)*y(k,212)
         mat(k,863) = -rxt(k,460)*y(k,212)
         mat(k,891) = -rxt(k,461)*y(k,212)
         mat(k,256) = -rxt(k,463)*y(k,212)
         mat(k,400) = -rxt(k,468)*y(k,212)
         mat(k,1087) = -rxt(k,473)*y(k,212)
         mat(k,729) = -rxt(k,474)*y(k,212)
         mat(k,198) = -(rxt(k,477) + rxt(k,479)) * y(k,212)
         mat(k,50) = -rxt(k,478)*y(k,212)
         mat(k,774) = mat(k,774) + .630_r8*rxt(k,440)*y(k,134)
         mat(k,204) = mat(k,204) + .650_r8*rxt(k,277)*y(k,212)
         mat(k,456) = mat(k,456) + .130_r8*rxt(k,279)*y(k,134)
         mat(k,249) = mat(k,249) + .500_r8*rxt(k,285)*y(k,212)
         mat(k,849) = mat(k,849) + .360_r8*rxt(k,308)*y(k,134)
         mat(k,1909) = mat(k,1909) + rxt(k,256)*y(k,133)
         mat(k,322) = mat(k,322) + .300_r8*rxt(k,263)*y(k,212)
         mat(k,1731) = rxt(k,179)*y(k,201)
         mat(k,654) = rxt(k,233)*y(k,222)
         mat(k,1281) = rxt(k,140)*y(k,134) + 2.000_r8*rxt(k,135)*y(k,201)
         mat(k,1013) = mat(k,1013) + rxt(k,132)*y(k,133) + rxt(k,124)*y(k,211)
         mat(k,479) = mat(k,479) + rxt(k,133)*y(k,133)
         mat(k,699) = mat(k,699) + rxt(k,223)*y(k,133) + rxt(k,229)*y(k,211)
         mat(k,1300) = mat(k,1300) + rxt(k,194)*y(k,133) + rxt(k,206)*y(k,211)
         mat(k,99) = mat(k,99) + rxt(k,274)*y(k,211)
         mat(k,643) = rxt(k,225)*y(k,133)
         mat(k,679) = mat(k,679) + rxt(k,197)*y(k,133)
         mat(k,720) = mat(k,720) + .320_r8*rxt(k,385)*y(k,134)
         mat(k,557) = mat(k,557) + .600_r8*rxt(k,387)*y(k,212)
         mat(k,1067) = mat(k,1067) + .240_r8*rxt(k,338)*y(k,134)
         mat(k,218) = mat(k,218) + .100_r8*rxt(k,340)*y(k,212)
         mat(k,800) = mat(k,800) + .630_r8*rxt(k,443)*y(k,134)
         mat(k,1137) = mat(k,1137) + .360_r8*rxt(k,352)*y(k,134)
         mat(k,1987) = rxt(k,163)*y(k,201)
         mat(k,1450) = mat(k,1450) + rxt(k,158)*y(k,201)
         mat(k,1886) = mat(k,1886) + rxt(k,256)*y(k,42) + rxt(k,132)*y(k,77) &
                      + rxt(k,133)*y(k,79) + rxt(k,223)*y(k,81) + rxt(k,194)*y(k,85) &
                      + rxt(k,225)*y(k,91) + rxt(k,197)*y(k,92) + rxt(k,138)*y(k,201)
         mat(k,1832) = mat(k,1832) + .630_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .360_r8*rxt(k,308)*y(k,29) + rxt(k,140)*y(k,76) &
                      + .320_r8*rxt(k,385)*y(k,98) + .240_r8*rxt(k,338)*y(k,105) &
                      + .630_r8*rxt(k,443)*y(k,110) + .360_r8*rxt(k,352)*y(k,111) &
                      + rxt(k,139)*y(k,201)
         mat(k,440) = mat(k,440) + .500_r8*rxt(k,320)*y(k,212)
         mat(k,127) = mat(k,127) + .500_r8*rxt(k,395)*y(k,212)
         mat(k,412) = .400_r8*rxt(k,396)*y(k,201)
         mat(k,1241) = .450_r8*rxt(k,293)*y(k,201)
         mat(k,625) = .400_r8*rxt(k,410)*y(k,201)
         mat(k,1393) = mat(k,1393) + rxt(k,179)*y(k,56) + 2.000_r8*rxt(k,135)*y(k,76) &
                      + rxt(k,163)*y(k,124) + rxt(k,158)*y(k,126) + rxt(k,138) &
                      *y(k,133) + rxt(k,139)*y(k,134) + .400_r8*rxt(k,396)*y(k,188) &
                      + .450_r8*rxt(k,293)*y(k,195) + .400_r8*rxt(k,410)*y(k,197) &
                      + .450_r8*rxt(k,343)*y(k,207) + .400_r8*rxt(k,416)*y(k,208) &
                      + .200_r8*rxt(k,347)*y(k,209) + .150_r8*rxt(k,322)*y(k,215)
         mat(k,1211) = .450_r8*rxt(k,343)*y(k,201)
         mat(k,748) = .400_r8*rxt(k,416)*y(k,201)
         mat(k,548) = .200_r8*rxt(k,347)*y(k,201)
         mat(k,1474) = rxt(k,124)*y(k,77) + rxt(k,229)*y(k,81) + rxt(k,206)*y(k,85) &
                      + rxt(k,274)*y(k,86) + 2.000_r8*rxt(k,125)*y(k,222)
         mat(k,1621) = mat(k,1621) + .650_r8*rxt(k,277)*y(k,24) + .500_r8*rxt(k,285) &
                      *y(k,27) + .300_r8*rxt(k,263)*y(k,53) + .600_r8*rxt(k,387) &
                      *y(k,103) + .100_r8*rxt(k,340)*y(k,106) + .500_r8*rxt(k,320) &
                      *y(k,146) + .500_r8*rxt(k,395)*y(k,181)
         mat(k,989) = .150_r8*rxt(k,322)*y(k,201)
         mat(k,2012) = rxt(k,233)*y(k,73) + 2.000_r8*rxt(k,125)*y(k,211)
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
         mat(k,347) = -(rxt(k,419)*y(k,201) + rxt(k,420)*y(k,124))
         mat(k,1329) = -rxt(k,419)*y(k,213)
         mat(k,1929) = -rxt(k,420)*y(k,213)
         mat(k,117) = .200_r8*rxt(k,409)*y(k,212)
         mat(k,92) = .140_r8*rxt(k,421)*y(k,212)
         mat(k,235) = rxt(k,424)*y(k,212)
         mat(k,1537) = .200_r8*rxt(k,409)*y(k,66) + .140_r8*rxt(k,421)*y(k,142) &
                      + rxt(k,424)*y(k,143)
         mat(k,630) = -(rxt(k,318)*y(k,201) + rxt(k,319)*y(k,124))
         mat(k,1353) = -rxt(k,318)*y(k,214)
         mat(k,1949) = -rxt(k,319)*y(k,214)
         mat(k,838) = rxt(k,325)*y(k,212)
         mat(k,436) = .500_r8*rxt(k,320)*y(k,212)
         mat(k,1571) = rxt(k,325)*y(k,29) + .500_r8*rxt(k,320)*y(k,146)
         mat(k,985) = -(rxt(k,321)*y(k,196) + rxt(k,322)*y(k,201) + rxt(k,323) &
                      *y(k,124))
         mat(k,1653) = -rxt(k,321)*y(k,215)
         mat(k,1373) = -rxt(k,322)*y(k,215)
         mat(k,1968) = -rxt(k,323)*y(k,215)
         mat(k,769) = .060_r8*rxt(k,440)*y(k,134)
         mat(k,820) = rxt(k,326)*y(k,212)
         mat(k,795) = .060_r8*rxt(k,443)*y(k,134)
         mat(k,1813) = .060_r8*rxt(k,440)*y(k,6) + .060_r8*rxt(k,443)*y(k,110)
         mat(k,303) = rxt(k,324)*y(k,212)
         mat(k,888) = .150_r8*rxt(k,461)*y(k,212)
         mat(k,1600) = rxt(k,326)*y(k,48) + rxt(k,324)*y(k,147) + .150_r8*rxt(k,461) &
                      *y(k,178)
         mat(k,950) = -(rxt(k,450)*y(k,196) + rxt(k,451)*y(k,201) + rxt(k,452) &
                      *y(k,124))
         mat(k,1651) = -rxt(k,450)*y(k,216)
         mat(k,1371) = -rxt(k,451)*y(k,216)
         mat(k,1966) = -rxt(k,452)*y(k,216)
         mat(k,1428) = .500_r8*rxt(k,459)*y(k,177)
         mat(k,488) = rxt(k,453)*y(k,212)
         mat(k,859) = .500_r8*rxt(k,459)*y(k,126) + rxt(k,460)*y(k,212)
         mat(k,1598) = rxt(k,453)*y(k,174) + rxt(k,460)*y(k,177)
         mat(k,928) = -(rxt(k,455)*y(k,196) + rxt(k,456)*y(k,201) + rxt(k,457) &
                      *y(k,124))
         mat(k,1650) = -rxt(k,455)*y(k,217)
         mat(k,1370) = -rxt(k,456)*y(k,217)
         mat(k,1965) = -rxt(k,457)*y(k,217)
         mat(k,767) = rxt(k,441)*y(k,212)
         mat(k,793) = rxt(k,444)*y(k,212)
         mat(k,370) = rxt(k,458)*y(k,212)
         mat(k,1597) = rxt(k,441)*y(k,6) + rxt(k,444)*y(k,110) + rxt(k,458)*y(k,176)
         mat(k,594) = -(rxt(k,426)*y(k,201) + rxt(k,427)*y(k,124))
         mat(k,1350) = -rxt(k,426)*y(k,218)
         mat(k,1946) = -rxt(k,427)*y(k,218)
         mat(k,497) = rxt(k,428)*y(k,212)
         mat(k,113) = .650_r8*rxt(k,429)*y(k,212)
         mat(k,1568) = rxt(k,428)*y(k,179) + .650_r8*rxt(k,429)*y(k,180)
         mat(k,1024) = -(rxt(k,390)*y(k,195) + rxt(k,391)*y(k,196) + rxt(k,392) &
                      *y(k,201) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126))
         mat(k,1227) = -rxt(k,390)*y(k,219)
         mat(k,1655) = -rxt(k,391)*y(k,219)
         mat(k,1376) = -rxt(k,392)*y(k,219)
         mat(k,1970) = -rxt(k,393)*y(k,219)
         mat(k,1432) = -rxt(k,394)*y(k,219)
         mat(k,152) = rxt(k,362)*y(k,212)
         mat(k,212) = rxt(k,363)*y(k,212)
         mat(k,69) = rxt(k,364)*y(k,212)
         mat(k,553) = .400_r8*rxt(k,387)*y(k,212)
         mat(k,126) = .500_r8*rxt(k,395)*y(k,212)
         mat(k,1603) = rxt(k,362)*y(k,94) + rxt(k,363)*y(k,96) + rxt(k,364)*y(k,97) &
                      + .400_r8*rxt(k,387)*y(k,103) + .500_r8*rxt(k,395)*y(k,181)
         mat(k,610) = -(rxt(k,432)*y(k,201) + rxt(k,433)*y(k,124))
         mat(k,1351) = -rxt(k,432)*y(k,220)
         mat(k,1947) = -rxt(k,433)*y(k,220)
         mat(k,137) = .560_r8*rxt(k,431)*y(k,212)
         mat(k,565) = rxt(k,434)*y(k,212)
         mat(k,1569) = .560_r8*rxt(k,431)*y(k,182) + rxt(k,434)*y(k,183)
         mat(k,384) = -(rxt(k,435)*y(k,201) + rxt(k,436)*y(k,124))
         mat(k,1334) = -rxt(k,435)*y(k,221)
         mat(k,1933) = -rxt(k,436)*y(k,221)
         mat(k,144) = .300_r8*rxt(k,437)*y(k,212)
         mat(k,327) = rxt(k,438)*y(k,212)
         mat(k,1543) = .300_r8*rxt(k,437)*y(k,184) + rxt(k,438)*y(k,185)
         mat(k,2022) = -(rxt(k,125)*y(k,211) + rxt(k,233)*y(k,73) + rxt(k,475) &
                      *y(k,152))
         mat(k,1484) = -rxt(k,125)*y(k,222)
         mat(k,656) = -rxt(k,233)*y(k,222)
         mat(k,171) = -rxt(k,475)*y(k,222)
         mat(k,210) = rxt(k,287)*y(k,212)
         mat(k,319) = rxt(k,312)*y(k,212)
         mat(k,58) = rxt(k,313)*y(k,212)
         mat(k,1919) = rxt(k,257)*y(k,212)
         mat(k,1004) = rxt(k,289)*y(k,212)
         mat(k,824) = rxt(k,326)*y(k,212)
         mat(k,1079) = rxt(k,315)*y(k,212)
         mat(k,434) = rxt(k,295)*y(k,212)
         mat(k,396) = rxt(k,296)*y(k,212)
         mat(k,325) = rxt(k,263)*y(k,212)
         mat(k,1287) = rxt(k,136)*y(k,201)
         mat(k,1017) = rxt(k,141)*y(k,212)
         mat(k,482) = rxt(k,142)*y(k,212)
         mat(k,702) = rxt(k,224)*y(k,212)
         mat(k,1309) = (rxt(k,516)+rxt(k,521))*y(k,91) + (rxt(k,509)+rxt(k,515) &
                       +rxt(k,520))*y(k,92) + rxt(k,195)*y(k,212)
         mat(k,660) = rxt(k,267)*y(k,212)
         mat(k,1273) = rxt(k,171)*y(k,212)
         mat(k,277) = rxt(k,149)*y(k,212)
         mat(k,647) = (rxt(k,516)+rxt(k,521))*y(k,85)
         mat(k,683) = (rxt(k,509)+rxt(k,515)+rxt(k,520))*y(k,85) + rxt(k,198)*y(k,212)
         mat(k,1070) = .500_r8*rxt(k,339)*y(k,212)
         mat(k,51) = rxt(k,478)*y(k,212)
         mat(k,442) = rxt(k,320)*y(k,212)
         mat(k,307) = rxt(k,324)*y(k,212)
         mat(k,1403) = rxt(k,136)*y(k,76) + rxt(k,143)*y(k,212)
         mat(k,1631) = rxt(k,287)*y(k,28) + rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31) &
                      + rxt(k,257)*y(k,42) + rxt(k,289)*y(k,45) + rxt(k,326)*y(k,48) &
                      + rxt(k,315)*y(k,49) + rxt(k,295)*y(k,50) + rxt(k,296)*y(k,51) &
                      + rxt(k,263)*y(k,53) + rxt(k,141)*y(k,77) + rxt(k,142)*y(k,79) &
                      + rxt(k,224)*y(k,81) + rxt(k,195)*y(k,85) + rxt(k,267)*y(k,87) &
                      + rxt(k,171)*y(k,89) + rxt(k,149)*y(k,90) + rxt(k,198)*y(k,92) &
                      + .500_r8*rxt(k,339)*y(k,105) + rxt(k,478)*y(k,120) + rxt(k,320) &
                      *y(k,146) + rxt(k,324)*y(k,147) + rxt(k,143)*y(k,201) &
                      + 2.000_r8*rxt(k,146)*y(k,212)
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
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 87) = lmat(k, 87)
         mat(k, 88) = lmat(k, 88)
         mat(k, 89) = lmat(k, 89)
         mat(k, 90) = lmat(k, 90)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 121) = lmat(k, 121)
         mat(k, 122) = lmat(k, 122)
         mat(k, 123) = lmat(k, 123)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 147) = lmat(k, 147)
         mat(k, 148) = lmat(k, 148)
         mat(k, 149) = lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 151) = lmat(k, 151)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 155) = lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = lmat(k, 162)
         mat(k, 163) = lmat(k, 163)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 169) = lmat(k, 169)
         mat(k, 170) = lmat(k, 170)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 181) = lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = lmat(k, 185)
         mat(k, 186) = lmat(k, 186)
         mat(k, 187) = lmat(k, 187)
         mat(k, 188) = lmat(k, 188)
         mat(k, 189) = lmat(k, 189)
         mat(k, 190) = lmat(k, 190)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
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
         mat(k, 237) = lmat(k, 237)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 239) = lmat(k, 239)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 248) = lmat(k, 248)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 251) = lmat(k, 251)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 278) = mat(k, 278) + lmat(k, 278)
         mat(k, 281) = lmat(k, 281)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 285) = lmat(k, 285)
         mat(k, 287) = lmat(k, 287)
         mat(k, 288) = lmat(k, 288)
         mat(k, 289) = mat(k, 289) + lmat(k, 289)
         mat(k, 290) = lmat(k, 290)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = lmat(k, 292)
         mat(k, 293) = lmat(k, 293)
         mat(k, 294) = lmat(k, 294)
         mat(k, 295) = lmat(k, 295)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 297) = lmat(k, 297)
         mat(k, 299) = lmat(k, 299)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 301) = lmat(k, 301)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 304) = lmat(k, 304)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 313) = lmat(k, 313)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 316) = lmat(k, 316)
         mat(k, 317) = lmat(k, 317)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = mat(k, 322) + lmat(k, 322)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 328) = lmat(k, 328)
         mat(k, 329) = lmat(k, 329)
         mat(k, 330) = lmat(k, 330)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 344) = lmat(k, 344)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 353) = lmat(k, 353)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = lmat(k, 355)
         mat(k, 356) = mat(k, 356) + lmat(k, 356)
         mat(k, 357) = lmat(k, 357)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 360) = lmat(k, 360)
         mat(k, 361) = lmat(k, 361)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 368) = lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 371) = lmat(k, 371)
         mat(k, 372) = mat(k, 372) + lmat(k, 372)
         mat(k, 373) = lmat(k, 373)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 384) = mat(k, 384) + lmat(k, 384)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = lmat(k, 394)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 402) = lmat(k, 402)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 405) = lmat(k, 405)
         mat(k, 406) = lmat(k, 406)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 417) = lmat(k, 417)
         mat(k, 421) = lmat(k, 421)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 428) = lmat(k, 428)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 438) = lmat(k, 438)
         mat(k, 439) = lmat(k, 439)
         mat(k, 440) = mat(k, 440) + lmat(k, 440)
         mat(k, 441) = lmat(k, 441)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 444) = lmat(k, 444)
         mat(k, 445) = lmat(k, 445)
         mat(k, 446) = lmat(k, 446)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = lmat(k, 448)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 471) = lmat(k, 471)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 483) = mat(k, 483) + lmat(k, 483)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = lmat(k, 486)
         mat(k, 487) = lmat(k, 487)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = mat(k, 490) + lmat(k, 490)
         mat(k, 491) = lmat(k, 491)
         mat(k, 492) = lmat(k, 492)
         mat(k, 493) = lmat(k, 493)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 495) = lmat(k, 495)
         mat(k, 499) = lmat(k, 499)
         mat(k, 500) = lmat(k, 500)
         mat(k, 501) = lmat(k, 501)
         mat(k, 502) = mat(k, 502) + lmat(k, 502)
         mat(k, 503) = lmat(k, 503)
         mat(k, 504) = lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 509) = lmat(k, 509)
         mat(k, 510) = lmat(k, 510)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 521) = lmat(k, 521)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 526) = lmat(k, 526)
         mat(k, 527) = lmat(k, 527)
         mat(k, 529) = lmat(k, 529)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = lmat(k, 532)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 539) = lmat(k, 539)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 554) = lmat(k, 554)
         mat(k, 555) = lmat(k, 555)
         mat(k, 556) = lmat(k, 556)
         mat(k, 558) = lmat(k, 558)
         mat(k, 559) = lmat(k, 559)
         mat(k, 560) = lmat(k, 560)
         mat(k, 561) = lmat(k, 561)
         mat(k, 562) = lmat(k, 562)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 567) = lmat(k, 567)
         mat(k, 570) = lmat(k, 570)
         mat(k, 572) = lmat(k, 572)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 583) = mat(k, 583) + lmat(k, 583)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 621) = mat(k, 621) + lmat(k, 621)
         mat(k, 630) = mat(k, 630) + lmat(k, 630)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 641) = lmat(k, 641)
         mat(k, 643) = mat(k, 643) + lmat(k, 643)
         mat(k, 648) = mat(k, 648) + lmat(k, 648)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 655) = lmat(k, 655)
         mat(k, 657) = mat(k, 657) + lmat(k, 657)
         mat(k, 661) = mat(k, 661) + lmat(k, 661)
         mat(k, 666) = mat(k, 666) + lmat(k, 666)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 679) = mat(k, 679) + lmat(k, 679)
         mat(k, 681) = mat(k, 681) + lmat(k, 681)
         mat(k, 687) = mat(k, 687) + lmat(k, 687)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 696) = mat(k, 696) + lmat(k, 696)
         mat(k, 697) = mat(k, 697) + lmat(k, 697)
         mat(k, 709) = mat(k, 709) + lmat(k, 709)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 727) = lmat(k, 727)
         mat(k, 730) = lmat(k, 730)
         mat(k, 732) = mat(k, 732) + lmat(k, 732)
         mat(k, 734) = lmat(k, 734)
         mat(k, 736) = mat(k, 736) + lmat(k, 736)
         mat(k, 737) = lmat(k, 737)
         mat(k, 738) = lmat(k, 738)
         mat(k, 739) = mat(k, 739) + lmat(k, 739)
         mat(k, 740) = mat(k, 740) + lmat(k, 740)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 761) = mat(k, 761) + lmat(k, 761)
         mat(k, 787) = mat(k, 787) + lmat(k, 787)
         mat(k, 809) = mat(k, 809) + lmat(k, 809)
         mat(k, 819) = mat(k, 819) + lmat(k, 819)
         mat(k, 821) = lmat(k, 821)
         mat(k, 823) = lmat(k, 823)
         mat(k, 826) = mat(k, 826) + lmat(k, 826)
         mat(k, 827) = mat(k, 827) + lmat(k, 827)
         mat(k, 828) = mat(k, 828) + lmat(k, 828)
         mat(k, 830) = mat(k, 830) + lmat(k, 830)
         mat(k, 832) = mat(k, 832) + lmat(k, 832)
         mat(k, 833) = mat(k, 833) + lmat(k, 833)
         mat(k, 834) = lmat(k, 834)
         mat(k, 841) = mat(k, 841) + lmat(k, 841)
         mat(k, 856) = mat(k, 856) + lmat(k, 856)
         mat(k, 857) = lmat(k, 857)
         mat(k, 858) = lmat(k, 858)
         mat(k, 861) = lmat(k, 861)
         mat(k, 865) = mat(k, 865) + lmat(k, 865)
         mat(k, 866) = lmat(k, 866)
         mat(k, 867) = mat(k, 867) + lmat(k, 867)
         mat(k, 869) = mat(k, 869) + lmat(k, 869)
         mat(k, 870) = lmat(k, 870)
         mat(k, 874) = mat(k, 874) + lmat(k, 874)
         mat(k, 878) = lmat(k, 878)
         mat(k, 879) = mat(k, 879) + lmat(k, 879)
         mat(k, 882) = lmat(k, 882)
         mat(k, 885) = mat(k, 885) + lmat(k, 885)
         mat(k, 886) = mat(k, 886) + lmat(k, 886)
         mat(k, 887) = mat(k, 887) + lmat(k, 887)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 890) = mat(k, 890) + lmat(k, 890)
         mat(k, 893) = mat(k, 893) + lmat(k, 893)
         mat(k, 896) = lmat(k, 896)
         mat(k, 897) = lmat(k, 897)
         mat(k, 898) = mat(k, 898) + lmat(k, 898)
         mat(k, 899) = lmat(k, 899)
         mat(k, 900) = lmat(k, 900)
         mat(k, 902) = lmat(k, 902)
         mat(k, 903) = lmat(k, 903)
         mat(k, 904) = lmat(k, 904)
         mat(k, 905) = mat(k, 905) + lmat(k, 905)
         mat(k, 908) = lmat(k, 908)
         mat(k, 909) = lmat(k, 909)
         mat(k, 911) = mat(k, 911) + lmat(k, 911)
         mat(k, 915) = mat(k, 915) + lmat(k, 915)
         mat(k, 917) = lmat(k, 917)
         mat(k, 918) = mat(k, 918) + lmat(k, 918)
         mat(k, 920) = lmat(k, 920)
         mat(k, 928) = mat(k, 928) + lmat(k, 928)
         mat(k, 950) = mat(k, 950) + lmat(k, 950)
         mat(k, 969) = mat(k, 969) + lmat(k, 969)
         mat(k, 985) = mat(k, 985) + lmat(k, 985)
         mat(k, 995) = lmat(k, 995)
         mat(k, 996) = mat(k, 996) + lmat(k, 996)
         mat(k,1000) = lmat(k,1000)
         mat(k,1003) = lmat(k,1003)
         mat(k,1007) = mat(k,1007) + lmat(k,1007)
         mat(k,1024) = mat(k,1024) + lmat(k,1024)
         mat(k,1044) = mat(k,1044) + lmat(k,1044)
         mat(k,1059) = mat(k,1059) + lmat(k,1059)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1064) = mat(k,1064) + lmat(k,1064)
         mat(k,1066) = mat(k,1066) + lmat(k,1066)
         mat(k,1069) = mat(k,1069) + lmat(k,1069)
         mat(k,1071) = mat(k,1071) + lmat(k,1071)
         mat(k,1072) = mat(k,1072) + lmat(k,1072)
         mat(k,1073) = mat(k,1073) + lmat(k,1073)
         mat(k,1076) = lmat(k,1076)
         mat(k,1081) = lmat(k,1081)
         mat(k,1082) = mat(k,1082) + lmat(k,1082)
         mat(k,1083) = mat(k,1083) + lmat(k,1083)
         mat(k,1093) = lmat(k,1093)
         mat(k,1108) = mat(k,1108) + lmat(k,1108)
         mat(k,1125) = lmat(k,1125)
         mat(k,1126) = mat(k,1126) + lmat(k,1126)
         mat(k,1130) = mat(k,1130) + lmat(k,1130)
         mat(k,1132) = mat(k,1132) + lmat(k,1132)
         mat(k,1138) = lmat(k,1138)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1166) = lmat(k,1166)
         mat(k,1186) = mat(k,1186) + lmat(k,1186)
         mat(k,1191) = mat(k,1191) + lmat(k,1191)
         mat(k,1205) = mat(k,1205) + lmat(k,1205)
         mat(k,1236) = mat(k,1236) + lmat(k,1236)
         mat(k,1250) = mat(k,1250) + lmat(k,1250)
         mat(k,1263) = mat(k,1263) + lmat(k,1263)
         mat(k,1268) = mat(k,1268) + lmat(k,1268)
         mat(k,1271) = lmat(k,1271)
         mat(k,1276) = mat(k,1276) + lmat(k,1276)
         mat(k,1278) = mat(k,1278) + lmat(k,1278)
         mat(k,1295) = mat(k,1295) + lmat(k,1295)
         mat(k,1296) = mat(k,1296) + lmat(k,1296)
         mat(k,1303) = mat(k,1303) + lmat(k,1303)
         mat(k,1339) = mat(k,1339) + lmat(k,1339)
         mat(k,1390) = mat(k,1390) + lmat(k,1390)
         mat(k,1444) = mat(k,1444) + lmat(k,1444)
         mat(k,1448) = mat(k,1448) + lmat(k,1448)
         mat(k,1454) = mat(k,1454) + lmat(k,1454)
         mat(k,1457) = mat(k,1457) + lmat(k,1457)
         mat(k,1459) = mat(k,1459) + lmat(k,1459)
         mat(k,1461) = mat(k,1461) + lmat(k,1461)
         mat(k,1463) = mat(k,1463) + lmat(k,1463)
         mat(k,1464) = mat(k,1464) + lmat(k,1464)
         mat(k,1466) = mat(k,1466) + lmat(k,1466)
         mat(k,1467) = mat(k,1467) + lmat(k,1467)
         mat(k,1469) = mat(k,1469) + lmat(k,1469)
         mat(k,1471) = lmat(k,1471)
         mat(k,1473) = mat(k,1473) + lmat(k,1473)
         mat(k,1474) = mat(k,1474) + lmat(k,1474)
         mat(k,1475) = lmat(k,1475)
         mat(k,1477) = mat(k,1477) + lmat(k,1477)
         mat(k,1481) = lmat(k,1481)
         mat(k,1482) = lmat(k,1482)
         mat(k,1483) = lmat(k,1483)
         mat(k,1495) = lmat(k,1495)
         mat(k,1501) = lmat(k,1501)
         mat(k,1614) = mat(k,1614) + lmat(k,1614)
         mat(k,1618) = mat(k,1618) + lmat(k,1618)
         mat(k,1621) = mat(k,1621) + lmat(k,1621)
         mat(k,1622) = mat(k,1622) + lmat(k,1622)
         mat(k,1624) = mat(k,1624) + lmat(k,1624)
         mat(k,1631) = mat(k,1631) + lmat(k,1631)
         mat(k,1672) = mat(k,1672) + lmat(k,1672)
         mat(k,1699) = mat(k,1699) + lmat(k,1699)
         mat(k,1700) = mat(k,1700) + lmat(k,1700)
         mat(k,1704) = mat(k,1704) + lmat(k,1704)
         mat(k,1720) = mat(k,1720) + lmat(k,1720)
         mat(k,1724) = lmat(k,1724)
         mat(k,1727) = mat(k,1727) + lmat(k,1727)
         mat(k,1728) = mat(k,1728) + lmat(k,1728)
         mat(k,1732) = lmat(k,1732)
         mat(k,1734) = mat(k,1734) + lmat(k,1734)
         mat(k,1766) = mat(k,1766) + lmat(k,1766)
         mat(k,1772) = mat(k,1772) + lmat(k,1772)
         mat(k,1776) = mat(k,1776) + lmat(k,1776)
         mat(k,1779) = mat(k,1779) + lmat(k,1779)
         mat(k,1781) = mat(k,1781) + lmat(k,1781)
         mat(k,1831) = mat(k,1831) + lmat(k,1831)
         mat(k,1837) = mat(k,1837) + lmat(k,1837)
         mat(k,1839) = mat(k,1839) + lmat(k,1839)
         mat(k,1849) = mat(k,1849) + lmat(k,1849)
         mat(k,1862) = mat(k,1862) + lmat(k,1862)
         mat(k,1863) = mat(k,1863) + lmat(k,1863)
         mat(k,1891) = mat(k,1891) + lmat(k,1891)
         mat(k,1893) = mat(k,1893) + lmat(k,1893)
         mat(k,1900) = mat(k,1900) + lmat(k,1900)
         mat(k,1901) = lmat(k,1901)
         mat(k,1904) = mat(k,1904) + lmat(k,1904)
         mat(k,1917) = mat(k,1917) + lmat(k,1917)
         mat(k,1926) = mat(k,1926) + lmat(k,1926)
         mat(k,1994) = mat(k,1994) + lmat(k,1994)
         mat(k,1996) = mat(k,1996) + lmat(k,1996)
         mat(k,2003) = lmat(k,2003)
         mat(k,2007) = lmat(k,2007)
         mat(k,2011) = mat(k,2011) + lmat(k,2011)
         mat(k,2012) = mat(k,2012) + lmat(k,2012)
         mat(k,2019) = lmat(k,2019)
         mat(k,2022) = mat(k,2022) + lmat(k,2022)
         mat(k, 138) = 0._r8
         mat(k, 139) = 0._r8
         mat(k, 242) = 0._r8
         mat(k, 335) = 0._r8
         mat(k, 337) = 0._r8
         mat(k, 350) = 0._r8
         mat(k, 377) = 0._r8
         mat(k, 380) = 0._r8
         mat(k, 388) = 0._r8
         mat(k, 496) = 0._r8
         mat(k, 498) = 0._r8
         mat(k, 517) = 0._r8
         mat(k, 519) = 0._r8
         mat(k, 523) = 0._r8
         mat(k, 524) = 0._r8
         mat(k, 528) = 0._r8
         mat(k, 534) = 0._r8
         mat(k, 535) = 0._r8
         mat(k, 538) = 0._r8
         mat(k, 564) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 568) = 0._r8
         mat(k, 569) = 0._r8
         mat(k, 571) = 0._r8
         mat(k, 577) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 595) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 598) = 0._r8
         mat(k, 601) = 0._r8
         mat(k, 609) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 614) = 0._r8
         mat(k, 616) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 631) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 638) = 0._r8
         mat(k, 644) = 0._r8
         mat(k, 670) = 0._r8
         mat(k, 675) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 760) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 777) = 0._r8
         mat(k, 786) = 0._r8
         mat(k, 788) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 803) = 0._r8
         mat(k, 807) = 0._r8
         mat(k, 808) = 0._r8
         mat(k, 812) = 0._r8
         mat(k, 814) = 0._r8
         mat(k, 815) = 0._r8
         mat(k, 836) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 851) = 0._r8
         mat(k, 854) = 0._r8
         mat(k, 855) = 0._r8
         mat(k, 873) = 0._r8
         mat(k, 875) = 0._r8
         mat(k, 876) = 0._r8
         mat(k, 877) = 0._r8
         mat(k, 881) = 0._r8
         mat(k, 883) = 0._r8
         mat(k, 884) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 926) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 929) = 0._r8
         mat(k, 930) = 0._r8
         mat(k, 931) = 0._r8
         mat(k, 932) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 934) = 0._r8
         mat(k, 936) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 952) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 956) = 0._r8
         mat(k, 961) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 970) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 972) = 0._r8
         mat(k, 973) = 0._r8
         mat(k, 976) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1011) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1030) = 0._r8
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
         mat(k,1053) = 0._r8
         mat(k,1065) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1099) = 0._r8
         mat(k,1100) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1110) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1120) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1172) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1240) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1262) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1267) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1277) = 0._r8
         mat(k,1279) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1331) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1360) = 0._r8
         mat(k,1364) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1367) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1410) = 0._r8
         mat(k,1417) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1430) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1460) = 0._r8
         mat(k,1468) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1558) = 0._r8
         mat(k,1567) = 0._r8
         mat(k,1570) = 0._r8
         mat(k,1572) = 0._r8
         mat(k,1583) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1641) = 0._r8
         mat(k,1665) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1667) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1670) = 0._r8
         mat(k,1671) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1678) = 0._r8
         mat(k,1681) = 0._r8
         mat(k,1691) = 0._r8
         mat(k,1692) = 0._r8
         mat(k,1695) = 0._r8
         mat(k,1696) = 0._r8
         mat(k,1702) = 0._r8
         mat(k,1707) = 0._r8
         mat(k,1713) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1719) = 0._r8
         mat(k,1721) = 0._r8
         mat(k,1723) = 0._r8
         mat(k,1725) = 0._r8
         mat(k,1730) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1738) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1741) = 0._r8
         mat(k,1751) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1759) = 0._r8
         mat(k,1760) = 0._r8
         mat(k,1761) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1773) = 0._r8
         mat(k,1775) = 0._r8
         mat(k,1780) = 0._r8
         mat(k,1782) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1799) = 0._r8
         mat(k,1802) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1809) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1811) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1815) = 0._r8
         mat(k,1819) = 0._r8
         mat(k,1821) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1823) = 0._r8
         mat(k,1826) = 0._r8
         mat(k,1828) = 0._r8
         mat(k,1842) = 0._r8
         mat(k,1850) = 0._r8
         mat(k,1851) = 0._r8
         mat(k,1852) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1855) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1861) = 0._r8
         mat(k,1864) = 0._r8
         mat(k,1866) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1874) = 0._r8
         mat(k,1880) = 0._r8
         mat(k,1885) = 0._r8
         mat(k,1887) = 0._r8
         mat(k,1896) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1908) = 0._r8
         mat(k,1910) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1913) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1915) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1953) = 0._r8
         mat(k,1981) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1983) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1997) = 0._r8
         mat(k,2002) = 0._r8
         mat(k,2004) = 0._r8
         mat(k,2005) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2009) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2013) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2015) = 0._r8
         mat(k,2016) = 0._r8
         mat(k,2017) = 0._r8
         mat(k,2018) = 0._r8
         mat(k,2020) = 0._r8
         mat(k,2021) = 0._r8
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
         mat(k, 81) = mat(k, 81) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 91) = mat(k, 91) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 109) = mat(k, 109) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 121) = mat(k, 121) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 150) = mat(k, 150) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 158) = mat(k, 158) - dti(k)
         mat(k, 161) = mat(k, 161) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 172) = mat(k, 172) - dti(k)
         mat(k, 176) = mat(k, 176) - dti(k)
         mat(k, 180) = mat(k, 180) - dti(k)
         mat(k, 184) = mat(k, 184) - dti(k)
         mat(k, 190) = mat(k, 190) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 205) = mat(k, 205) - dti(k)
         mat(k, 211) = mat(k, 211) - dti(k)
         mat(k, 214) = mat(k, 214) - dti(k)
         mat(k, 219) = mat(k, 219) - dti(k)
         mat(k, 224) = mat(k, 224) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 234) = mat(k, 234) - dti(k)
         mat(k, 240) = mat(k, 240) - dti(k)
         mat(k, 245) = mat(k, 245) - dti(k)
         mat(k, 250) = mat(k, 250) - dti(k)
         mat(k, 258) = mat(k, 258) - dti(k)
         mat(k, 266) = mat(k, 266) - dti(k)
         mat(k, 272) = mat(k, 272) - dti(k)
         mat(k, 278) = mat(k, 278) - dti(k)
         mat(k, 284) = mat(k, 284) - dti(k)
         mat(k, 290) = mat(k, 290) - dti(k)
         mat(k, 296) = mat(k, 296) - dti(k)
         mat(k, 302) = mat(k, 302) - dti(k)
         mat(k, 308) = mat(k, 308) - dti(k)
         mat(k, 314) = mat(k, 314) - dti(k)
         mat(k, 320) = mat(k, 320) - dti(k)
         mat(k, 326) = mat(k, 326) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 353) = mat(k, 353) - dti(k)
         mat(k, 356) = mat(k, 356) - dti(k)
         mat(k, 363) = mat(k, 363) - dti(k)
         mat(k, 367) = mat(k, 367) - dti(k)
         mat(k, 376) = mat(k, 376) - dti(k)
         mat(k, 384) = mat(k, 384) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 397) = mat(k, 397) - dti(k)
         mat(k, 403) = mat(k, 403) - dti(k)
         mat(k, 409) = mat(k, 409) - dti(k)
         mat(k, 415) = mat(k, 415) - dti(k)
         mat(k, 423) = mat(k, 423) - dti(k)
         mat(k, 431) = mat(k, 431) - dti(k)
         mat(k, 435) = mat(k, 435) - dti(k)
         mat(k, 443) = mat(k, 443) - dti(k)
         mat(k, 451) = mat(k, 451) - dti(k)
         mat(k, 459) = mat(k, 459) - dti(k)
         mat(k, 467) = mat(k, 467) - dti(k)
         mat(k, 476) = mat(k, 476) - dti(k)
         mat(k, 483) = mat(k, 483) - dti(k)
         mat(k, 494) = mat(k, 494) - dti(k)
         mat(k, 503) = mat(k, 503) - dti(k)
         mat(k, 507) = mat(k, 507) - dti(k)
         mat(k, 515) = mat(k, 515) - dti(k)
         mat(k, 522) = mat(k, 522) - dti(k)
         mat(k, 533) = mat(k, 533) - dti(k)
         mat(k, 544) = mat(k, 544) - dti(k)
         mat(k, 552) = mat(k, 552) - dti(k)
         mat(k, 563) = mat(k, 563) - dti(k)
         mat(k, 576) = mat(k, 576) - dti(k)
         mat(k, 583) = mat(k, 583) - dti(k)
         mat(k, 594) = mat(k, 594) - dti(k)
         mat(k, 610) = mat(k, 610) - dti(k)
         mat(k, 621) = mat(k, 621) - dti(k)
         mat(k, 630) = mat(k, 630) - dti(k)
         mat(k, 640) = mat(k, 640) - dti(k)
         mat(k, 649) = mat(k, 649) - dti(k)
         mat(k, 657) = mat(k, 657) - dti(k)
         mat(k, 661) = mat(k, 661) - dti(k)
         mat(k, 666) = mat(k, 666) - dti(k)
         mat(k, 677) = mat(k, 677) - dti(k)
         mat(k, 687) = mat(k, 687) - dti(k)
         mat(k, 695) = mat(k, 695) - dti(k)
         mat(k, 709) = mat(k, 709) - dti(k)
         mat(k, 726) = mat(k, 726) - dti(k)
         mat(k, 732) = mat(k, 732) - dti(k)
         mat(k, 739) = mat(k, 739) - dti(k)
         mat(k, 744) = mat(k, 744) - dti(k)
         mat(k, 761) = mat(k, 761) - dti(k)
         mat(k, 787) = mat(k, 787) - dti(k)
         mat(k, 809) = mat(k, 809) - dti(k)
         mat(k, 819) = mat(k, 819) - dti(k)
         mat(k, 827) = mat(k, 827) - dti(k)
         mat(k, 841) = mat(k, 841) - dti(k)
         mat(k, 856) = mat(k, 856) - dti(k)
         mat(k, 865) = mat(k, 865) - dti(k)
         mat(k, 874) = mat(k, 874) - dti(k)
         mat(k, 886) = mat(k, 886) - dti(k)
         mat(k, 898) = mat(k, 898) - dti(k)
         mat(k, 911) = mat(k, 911) - dti(k)
         mat(k, 915) = mat(k, 915) - dti(k)
         mat(k, 928) = mat(k, 928) - dti(k)
         mat(k, 950) = mat(k, 950) - dti(k)
         mat(k, 969) = mat(k, 969) - dti(k)
         mat(k, 985) = mat(k, 985) - dti(k)
         mat(k, 996) = mat(k, 996) - dti(k)
         mat(k,1007) = mat(k,1007) - dti(k)
         mat(k,1024) = mat(k,1024) - dti(k)
         mat(k,1044) = mat(k,1044) - dti(k)
         mat(k,1060) = mat(k,1060) - dti(k)
         mat(k,1072) = mat(k,1072) - dti(k)
         mat(k,1083) = mat(k,1083) - dti(k)
         mat(k,1108) = mat(k,1108) - dti(k)
         mat(k,1130) = mat(k,1130) - dti(k)
         mat(k,1153) = mat(k,1153) - dti(k)
         mat(k,1186) = mat(k,1186) - dti(k)
         mat(k,1205) = mat(k,1205) - dti(k)
         mat(k,1236) = mat(k,1236) - dti(k)
         mat(k,1250) = mat(k,1250) - dti(k)
         mat(k,1263) = mat(k,1263) - dti(k)
         mat(k,1276) = mat(k,1276) - dti(k)
         mat(k,1296) = mat(k,1296) - dti(k)
         mat(k,1390) = mat(k,1390) - dti(k)
         mat(k,1448) = mat(k,1448) - dti(k)
         mat(k,1473) = mat(k,1473) - dti(k)
         mat(k,1621) = mat(k,1621) - dti(k)
         mat(k,1672) = mat(k,1672) - dti(k)
         mat(k,1699) = mat(k,1699) - dti(k)
         mat(k,1734) = mat(k,1734) - dti(k)
         mat(k,1776) = mat(k,1776) - dti(k)
         mat(k,1837) = mat(k,1837) - dti(k)
         mat(k,1862) = mat(k,1862) - dti(k)
         mat(k,1893) = mat(k,1893) - dti(k)
         mat(k,1917) = mat(k,1917) - dti(k)
         mat(k,1996) = mat(k,1996) - dti(k)
         mat(k,2022) = mat(k,2022) - dti(k)
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
