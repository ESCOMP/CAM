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
         mat(k,666) = -(rxt(k,356)*y(k,217))
         mat(k,1631) = -rxt(k,356)*y(k,1)
         mat(k,1795) = rxt(k,359)*y(k,189)
         mat(k,896) = rxt(k,359)*y(k,124)
         mat(k,632) = -(rxt(k,360)*y(k,217))
         mat(k,1628) = -rxt(k,360)*y(k,2)
         mat(k,895) = rxt(k,357)*y(k,203)
         mat(k,1895) = rxt(k,357)*y(k,189)
         mat(k,999) = -(rxt(k,439)*y(k,126) + rxt(k,440)*y(k,134) + rxt(k,441) &
                      *y(k,217))
         mat(k,1719) = -rxt(k,439)*y(k,6)
         mat(k,2099) = -rxt(k,440)*y(k,6)
         mat(k,1661) = -rxt(k,441)*y(k,6)
         mat(k,164) = -(rxt(k,398)*y(k,217))
         mat(k,1559) = -rxt(k,398)*y(k,7)
         mat(k,421) = -(rxt(k,401)*y(k,217))
         mat(k,1600) = -rxt(k,401)*y(k,8)
         mat(k,481) = rxt(k,399)*y(k,203)
         mat(k,1880) = rxt(k,399)*y(k,191)
         mat(k,165) = .120_r8*rxt(k,398)*y(k,217)
         mat(k,1560) = .120_r8*rxt(k,398)*y(k,7)
         mat(k,991) = .100_r8*rxt(k,440)*y(k,134)
         mat(k,942) = .100_r8*rxt(k,443)*y(k,134)
         mat(k,2083) = .100_r8*rxt(k,440)*y(k,6) + .100_r8*rxt(k,443)*y(k,110)
         mat(k,1782) = .500_r8*rxt(k,400)*y(k,191) + .200_r8*rxt(k,427)*y(k,223) &
                      + .060_r8*rxt(k,433)*y(k,226)
         mat(k,482) = .500_r8*rxt(k,400)*y(k,124)
         mat(k,728) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,744) = .060_r8*rxt(k,433)*y(k,124)
         mat(k,1776) = .200_r8*rxt(k,427)*y(k,223) + .200_r8*rxt(k,433)*y(k,226)
         mat(k,727) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,742) = .200_r8*rxt(k,433)*y(k,124)
         mat(k,1792) = .200_r8*rxt(k,427)*y(k,223) + .150_r8*rxt(k,433)*y(k,226)
         mat(k,729) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,745) = .150_r8*rxt(k,433)*y(k,124)
         mat(k,1778) = .210_r8*rxt(k,433)*y(k,226)
         mat(k,743) = .210_r8*rxt(k,433)*y(k,124)
         mat(k,238) = -(rxt(k,361)*y(k,217))
         mat(k,1573) = -rxt(k,361)*y(k,15)
         mat(k,990) = .050_r8*rxt(k,440)*y(k,134)
         mat(k,941) = .050_r8*rxt(k,443)*y(k,134)
         mat(k,2082) = .050_r8*rxt(k,440)*y(k,6) + .050_r8*rxt(k,443)*y(k,110)
         mat(k,349) = -(rxt(k,327)*y(k,126) + rxt(k,328)*y(k,217))
         mat(k,1709) = -rxt(k,327)*y(k,16)
         mat(k,1590) = -rxt(k,328)*y(k,16)
         mat(k,1415) = -(rxt(k,211)*y(k,42) + rxt(k,212)*y(k,203) + rxt(k,213) &
                      *y(k,134))
         mat(k,1481) = -rxt(k,211)*y(k,17)
         mat(k,1941) = -rxt(k,212)*y(k,17)
         mat(k,2119) = -rxt(k,213)*y(k,17)
         mat(k,2211) = 4.000_r8*rxt(k,214)*y(k,19) + (rxt(k,215)+rxt(k,216))*y(k,59) &
                      + rxt(k,219)*y(k,124) + rxt(k,222)*y(k,133) + rxt(k,469) &
                      *y(k,150) + rxt(k,223)*y(k,217)
         mat(k,143) = rxt(k,201)*y(k,216)
         mat(k,149) = rxt(k,227)*y(k,216)
         mat(k,468) = 2.000_r8*rxt(k,238)*y(k,56) + 2.000_r8*rxt(k,250)*y(k,216) &
                      + 2.000_r8*rxt(k,239)*y(k,217)
         mat(k,595) = rxt(k,240)*y(k,56) + rxt(k,251)*y(k,216) + rxt(k,241)*y(k,217)
         mat(k,448) = 3.000_r8*rxt(k,245)*y(k,56) + 3.000_r8*rxt(k,228)*y(k,216) &
                      + 3.000_r8*rxt(k,246)*y(k,217)
         mat(k,2006) = 2.000_r8*rxt(k,238)*y(k,41) + rxt(k,240)*y(k,43) &
                      + 3.000_r8*rxt(k,245)*y(k,55)
         mat(k,1968) = (rxt(k,215)+rxt(k,216))*y(k,19)
         mat(k,107) = 2.000_r8*rxt(k,229)*y(k,216)
         mat(k,805) = rxt(k,224)*y(k,133) + rxt(k,230)*y(k,216) + rxt(k,225)*y(k,217)
         mat(k,1834) = rxt(k,219)*y(k,19)
         mat(k,2241) = rxt(k,222)*y(k,19) + rxt(k,224)*y(k,81)
         mat(k,1233) = rxt(k,469)*y(k,19)
         mat(k,1521) = rxt(k,201)*y(k,34) + rxt(k,227)*y(k,35) + 2.000_r8*rxt(k,250) &
                      *y(k,41) + rxt(k,251)*y(k,43) + 3.000_r8*rxt(k,228)*y(k,55) &
                      + 2.000_r8*rxt(k,229)*y(k,78) + rxt(k,230)*y(k,81)
         mat(k,1685) = rxt(k,223)*y(k,19) + 2.000_r8*rxt(k,239)*y(k,41) + rxt(k,241) &
                      *y(k,43) + 3.000_r8*rxt(k,246)*y(k,55) + rxt(k,225)*y(k,81)
         mat(k,2205) = rxt(k,217)*y(k,59)
         mat(k,1962) = rxt(k,217)*y(k,19)
         mat(k,2139) = (rxt(k,530)+rxt(k,535))*y(k,91)
         mat(k,777) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,2226) = -(4._r8*rxt(k,214)*y(k,19) + (rxt(k,215) + rxt(k,216) + rxt(k,217) &
                      ) * y(k,59) + rxt(k,218)*y(k,203) + rxt(k,219)*y(k,124) &
                      + rxt(k,220)*y(k,125) + rxt(k,222)*y(k,133) + rxt(k,223) &
                      *y(k,217) + rxt(k,469)*y(k,150))
         mat(k,1983) = -(rxt(k,215) + rxt(k,216) + rxt(k,217)) * y(k,19)
         mat(k,1957) = -rxt(k,218)*y(k,19)
         mat(k,1850) = -rxt(k,219)*y(k,19)
         mat(k,2202) = -rxt(k,220)*y(k,19)
         mat(k,2257) = -rxt(k,222)*y(k,19)
         mat(k,1701) = -rxt(k,223)*y(k,19)
         mat(k,1242) = -rxt(k,469)*y(k,19)
         mat(k,1422) = rxt(k,213)*y(k,134)
         mat(k,546) = rxt(k,221)*y(k,133)
         mat(k,809) = rxt(k,231)*y(k,216)
         mat(k,783) = rxt(k,226)*y(k,133)
         mat(k,2257) = mat(k,2257) + rxt(k,221)*y(k,20) + rxt(k,226)*y(k,91)
         mat(k,2135) = rxt(k,213)*y(k,17)
         mat(k,1537) = rxt(k,231)*y(k,81)
         mat(k,540) = -(rxt(k,221)*y(k,133))
         mat(k,2231) = -rxt(k,221)*y(k,20)
         mat(k,2207) = rxt(k,220)*y(k,125)
         mat(k,2169) = rxt(k,220)*y(k,19)
         mat(k,253) = -(rxt(k,402)*y(k,217))
         mat(k,1576) = -rxt(k,402)*y(k,22)
         mat(k,1774) = rxt(k,405)*y(k,193)
         mat(k,433) = rxt(k,405)*y(k,124)
         mat(k,339) = -(rxt(k,404)*y(k,217))
         mat(k,1588) = -rxt(k,404)*y(k,23)
         mat(k,434) = rxt(k,403)*y(k,203)
         mat(k,1872) = rxt(k,403)*y(k,193)
         mat(k,285) = -(rxt(k,276)*y(k,56) + rxt(k,277)*y(k,217))
         mat(k,1987) = -rxt(k,276)*y(k,24)
         mat(k,1580) = -rxt(k,277)*y(k,24)
         mat(k,548) = -(rxt(k,278)*y(k,56) + rxt(k,279)*y(k,134) + rxt(k,304)*y(k,217))
         mat(k,1992) = -rxt(k,278)*y(k,25)
         mat(k,2086) = -rxt(k,279)*y(k,25)
         mat(k,1617) = -rxt(k,304)*y(k,25)
         mat(k,261) = -(rxt(k,284)*y(k,217))
         mat(k,1578) = -rxt(k,284)*y(k,26)
         mat(k,812) = .800_r8*rxt(k,280)*y(k,194) + .200_r8*rxt(k,281)*y(k,198)
         mat(k,2025) = .200_r8*rxt(k,281)*y(k,194)
         mat(k,344) = -(rxt(k,285)*y(k,217))
         mat(k,1589) = -rxt(k,285)*y(k,27)
         mat(k,813) = rxt(k,282)*y(k,203)
         mat(k,1873) = rxt(k,282)*y(k,194)
         mat(k,291) = -(rxt(k,286)*y(k,56) + rxt(k,287)*y(k,217))
         mat(k,1988) = -rxt(k,286)*y(k,28)
         mat(k,1581) = -rxt(k,287)*y(k,28)
         mat(k,1024) = -(rxt(k,307)*y(k,126) + rxt(k,308)*y(k,134) + rxt(k,325) &
                      *y(k,217))
         mat(k,1720) = -rxt(k,307)*y(k,29)
         mat(k,2100) = -rxt(k,308)*y(k,29)
         mat(k,1662) = -rxt(k,325)*y(k,29)
         mat(k,843) = .130_r8*rxt(k,385)*y(k,134)
         mat(k,2100) = mat(k,2100) + .130_r8*rxt(k,385)*y(k,98)
         mat(k,415) = -(rxt(k,312)*y(k,217))
         mat(k,1599) = -rxt(k,312)*y(k,30)
         mat(k,790) = rxt(k,310)*y(k,203)
         mat(k,1879) = rxt(k,310)*y(k,195)
         mat(k,109) = -(rxt(k,313)*y(k,217))
         mat(k,1556) = -rxt(k,313)*y(k,31)
         mat(k,265) = -(rxt(k,408)*y(k,217))
         mat(k,1579) = -rxt(k,408)*y(k,32)
         mat(k,623) = rxt(k,406)*y(k,203)
         mat(k,1867) = rxt(k,406)*y(k,196)
         mat(k,99) = -(rxt(k,200)*y(k,216))
         mat(k,1499) = -rxt(k,200)*y(k,33)
         mat(k,141) = -(rxt(k,201)*y(k,216))
         mat(k,1504) = -rxt(k,201)*y(k,34)
         mat(k,146) = -(rxt(k,227)*y(k,216))
         mat(k,1505) = -rxt(k,227)*y(k,35)
         mat(k,113) = -(rxt(k,202)*y(k,216))
         mat(k,1501) = -rxt(k,202)*y(k,36)
         mat(k,151) = -(rxt(k,203)*y(k,216))
         mat(k,1506) = -rxt(k,203)*y(k,37)
         mat(k,117) = -(rxt(k,204)*y(k,216))
         mat(k,1502) = -rxt(k,204)*y(k,38)
         mat(k,156) = -(rxt(k,205)*y(k,216))
         mat(k,1507) = -rxt(k,205)*y(k,39)
         mat(k,121) = -(rxt(k,206)*y(k,216))
         mat(k,1503) = -rxt(k,206)*y(k,40)
         mat(k,467) = -(rxt(k,238)*y(k,56) + rxt(k,239)*y(k,217) + rxt(k,250)*y(k,216))
         mat(k,1991) = -rxt(k,238)*y(k,41)
         mat(k,1607) = -rxt(k,239)*y(k,41)
         mat(k,1516) = -rxt(k,250)*y(k,41)
         mat(k,1485) = -(rxt(k,175)*y(k,56) + rxt(k,211)*y(k,17) + rxt(k,255)*y(k,203) &
                      + rxt(k,256)*y(k,126) + rxt(k,257)*y(k,133) + rxt(k,258) &
                      *y(k,217))
         mat(k,2010) = -rxt(k,175)*y(k,42)
         mat(k,1417) = -rxt(k,211)*y(k,42)
         mat(k,1945) = -rxt(k,255)*y(k,42)
         mat(k,1746) = -rxt(k,256)*y(k,42)
         mat(k,2245) = -rxt(k,257)*y(k,42)
         mat(k,1689) = -rxt(k,258)*y(k,42)
         mat(k,672) = .400_r8*rxt(k,356)*y(k,217)
         mat(k,1009) = .340_r8*rxt(k,440)*y(k,134)
         mat(k,353) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,552) = rxt(k,279)*y(k,134)
         mat(k,1031) = .500_r8*rxt(k,308)*y(k,134)
         mat(k,613) = .500_r8*rxt(k,296)*y(k,217)
         mat(k,787) = rxt(k,263)*y(k,217)
         mat(k,387) = .300_r8*rxt(k,264)*y(k,217)
         mat(k,1433) = (rxt(k,272)+rxt(k,273))*y(k,216)
         mat(k,1971) = rxt(k,182)*y(k,198)
         mat(k,1098) = .800_r8*rxt(k,301)*y(k,217)
         mat(k,851) = .910_r8*rxt(k,385)*y(k,134)
         mat(k,586) = .300_r8*rxt(k,376)*y(k,217)
         mat(k,1199) = .800_r8*rxt(k,380)*y(k,198)
         mat(k,1214) = .120_r8*rxt(k,338)*y(k,134)
         mat(k,576) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,959) = .340_r8*rxt(k,443)*y(k,134)
         mat(k,1337) = .600_r8*rxt(k,352)*y(k,134)
         mat(k,1838) = .100_r8*rxt(k,358)*y(k,189) + rxt(k,262)*y(k,198) &
                      + .500_r8*rxt(k,329)*y(k,200) + .500_r8*rxt(k,298)*y(k,202) &
                      + .920_r8*rxt(k,368)*y(k,205) + .250_r8*rxt(k,336)*y(k,209) &
                      + rxt(k,345)*y(k,211) + rxt(k,319)*y(k,219) + rxt(k,323) &
                      *y(k,220) + .340_r8*rxt(k,452)*y(k,221) + .320_r8*rxt(k,457) &
                      *y(k,222) + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1746) = mat(k,1746) + .500_r8*rxt(k,327)*y(k,16) + rxt(k,369)*y(k,205) &
                      + .250_r8*rxt(k,335)*y(k,209) + rxt(k,346)*y(k,211)
         mat(k,2123) = .340_r8*rxt(k,440)*y(k,6) + rxt(k,279)*y(k,25) &
                      + .500_r8*rxt(k,308)*y(k,29) + .910_r8*rxt(k,385)*y(k,98) &
                      + .120_r8*rxt(k,338)*y(k,105) + .340_r8*rxt(k,443)*y(k,110) &
                      + .600_r8*rxt(k,352)*y(k,111)
         mat(k,527) = rxt(k,303)*y(k,217)
         mat(k,1063) = .680_r8*rxt(k,461)*y(k,217)
         mat(k,903) = .100_r8*rxt(k,358)*y(k,124)
         mat(k,817) = .700_r8*rxt(k,281)*y(k,198)
         mat(k,794) = rxt(k,309)*y(k,198)
         mat(k,1389) = rxt(k,292)*y(k,198) + rxt(k,365)*y(k,205) + .250_r8*rxt(k,332) &
                      *y(k,209) + rxt(k,341)*y(k,211) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,2062) = rxt(k,182)*y(k,59) + .800_r8*rxt(k,380)*y(k,101) + rxt(k,262) &
                      *y(k,124) + .700_r8*rxt(k,281)*y(k,194) + rxt(k,309)*y(k,195) &
                      + rxt(k,292)*y(k,197) + (4.000_r8*rxt(k,259)+2.000_r8*rxt(k,260)) &
                      *y(k,198) + 1.500_r8*rxt(k,366)*y(k,205) + .750_r8*rxt(k,371) &
                      *y(k,206) + .880_r8*rxt(k,333)*y(k,209) + 2.000_r8*rxt(k,342) &
                      *y(k,211) + .750_r8*rxt(k,445)*y(k,215) + .800_r8*rxt(k,321) &
                      *y(k,220) + .930_r8*rxt(k,450)*y(k,221) + .950_r8*rxt(k,455) &
                      *y(k,222) + .800_r8*rxt(k,391)*y(k,225)
         mat(k,568) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,716) = .500_r8*rxt(k,298)*y(k,124)
         mat(k,1945) = mat(k,1945) + .450_r8*rxt(k,343)*y(k,211) + .150_r8*rxt(k,322) &
                      *y(k,220)
         mat(k,1262) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) + rxt(k,365) &
                      *y(k,197) + 1.500_r8*rxt(k,366)*y(k,198)
         mat(k,1294) = .750_r8*rxt(k,371)*y(k,198)
         mat(k,1315) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,197) + .880_r8*rxt(k,333)*y(k,198)
         mat(k,1357) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,341)*y(k,197) &
                      + 2.000_r8*rxt(k,342)*y(k,198) + .450_r8*rxt(k,343)*y(k,203) &
                      + 4.000_r8*rxt(k,344)*y(k,211)
         mat(k,1050) = .750_r8*rxt(k,445)*y(k,198)
         mat(k,1525) = (rxt(k,272)+rxt(k,273))*y(k,54)
         mat(k,1689) = mat(k,1689) + .400_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,296) &
                      *y(k,51) + rxt(k,263)*y(k,52) + .300_r8*rxt(k,264)*y(k,53) &
                      + .800_r8*rxt(k,301)*y(k,74) + .300_r8*rxt(k,376)*y(k,99) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,139) &
                      + .680_r8*rxt(k,461)*y(k,178)
         mat(k,771) = rxt(k,319)*y(k,124)
         mat(k,1133) = rxt(k,323)*y(k,124) + .800_r8*rxt(k,321)*y(k,198) &
                      + .150_r8*rxt(k,322)*y(k,203)
         mat(k,1119) = .340_r8*rxt(k,452)*y(k,124) + .930_r8*rxt(k,450)*y(k,198)
         mat(k,916) = .320_r8*rxt(k,457)*y(k,124) + .950_r8*rxt(k,455)*y(k,198)
         mat(k,1176) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,390)*y(k,197) &
                      + .800_r8*rxt(k,391)*y(k,198)
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
         mat(k,594) = -(rxt(k,240)*y(k,56) + rxt(k,241)*y(k,217) + rxt(k,251)*y(k,216))
         mat(k,1993) = -rxt(k,240)*y(k,43)
         mat(k,1623) = -rxt(k,241)*y(k,43)
         mat(k,1517) = -rxt(k,251)*y(k,43)
         mat(k,125) = -(rxt(k,242)*y(k,217))
         mat(k,1557) = -rxt(k,242)*y(k,44)
         mat(k,1085) = -(rxt(k,288)*y(k,126) + rxt(k,289)*y(k,217))
         mat(k,1724) = -rxt(k,288)*y(k,45)
         mat(k,1666) = -rxt(k,289)*y(k,45)
         mat(k,670) = .800_r8*rxt(k,356)*y(k,217)
         mat(k,352) = rxt(k,327)*y(k,126)
         mat(k,262) = rxt(k,284)*y(k,217)
         mat(k,346) = .500_r8*rxt(k,285)*y(k,217)
         mat(k,1025) = .500_r8*rxt(k,308)*y(k,134)
         mat(k,1327) = .100_r8*rxt(k,352)*y(k,134)
         mat(k,1817) = .400_r8*rxt(k,358)*y(k,189) + rxt(k,283)*y(k,194) &
                      + .270_r8*rxt(k,311)*y(k,195) + rxt(k,329)*y(k,200) + rxt(k,348) &
                      *y(k,213) + rxt(k,319)*y(k,219)
         mat(k,1724) = mat(k,1724) + rxt(k,327)*y(k,16)
         mat(k,2103) = .500_r8*rxt(k,308)*y(k,29) + .100_r8*rxt(k,352)*y(k,111)
         mat(k,901) = .400_r8*rxt(k,358)*y(k,124)
         mat(k,816) = rxt(k,283)*y(k,124) + 3.200_r8*rxt(k,280)*y(k,194) &
                      + .800_r8*rxt(k,281)*y(k,198)
         mat(k,793) = .270_r8*rxt(k,311)*y(k,124)
         mat(k,2043) = .800_r8*rxt(k,281)*y(k,194)
         mat(k,566) = rxt(k,329)*y(k,124)
         mat(k,1924) = .200_r8*rxt(k,347)*y(k,213)
         mat(k,678) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,203)
         mat(k,1666) = mat(k,1666) + .800_r8*rxt(k,356)*y(k,1) + rxt(k,284)*y(k,26) &
                      + .500_r8*rxt(k,285)*y(k,27)
         mat(k,769) = rxt(k,319)*y(k,124)
         mat(k,365) = -(rxt(k,243)*y(k,56) + rxt(k,244)*y(k,217))
         mat(k,1989) = -rxt(k,243)*y(k,46)
         mat(k,1592) = -rxt(k,244)*y(k,46)
         mat(k,102) = -(rxt(k,290)*y(k,217))
         mat(k,1555) = -rxt(k,290)*y(k,47)
         mat(k,922) = -(rxt(k,326)*y(k,217))
         mat(k,1656) = -rxt(k,326)*y(k,48)
         mat(k,669) = .800_r8*rxt(k,356)*y(k,217)
         mat(k,995) = .520_r8*rxt(k,440)*y(k,134)
         mat(k,351) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,946) = .520_r8*rxt(k,443)*y(k,134)
         mat(k,1810) = .250_r8*rxt(k,358)*y(k,189) + .820_r8*rxt(k,311)*y(k,195) &
                      + .500_r8*rxt(k,329)*y(k,200) + .270_r8*rxt(k,452)*y(k,221) &
                      + .040_r8*rxt(k,457)*y(k,222)
         mat(k,1714) = .500_r8*rxt(k,327)*y(k,16)
         mat(k,2094) = .520_r8*rxt(k,440)*y(k,6) + .520_r8*rxt(k,443)*y(k,110)
         mat(k,1058) = .500_r8*rxt(k,461)*y(k,217)
         mat(k,900) = .250_r8*rxt(k,358)*y(k,124)
         mat(k,792) = .820_r8*rxt(k,311)*y(k,124) + .820_r8*rxt(k,309)*y(k,198)
         mat(k,2037) = .820_r8*rxt(k,309)*y(k,195) + .150_r8*rxt(k,450)*y(k,221) &
                      + .025_r8*rxt(k,455)*y(k,222)
         mat(k,565) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,1656) = mat(k,1656) + .800_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,461) &
                      *y(k,178)
         mat(k,1111) = .270_r8*rxt(k,452)*y(k,124) + .150_r8*rxt(k,450)*y(k,198)
         mat(k,913) = .040_r8*rxt(k,457)*y(k,124) + .025_r8*rxt(k,455)*y(k,198)
         mat(k,1221) = -(rxt(k,314)*y(k,126) + rxt(k,315)*y(k,217))
         mat(k,1734) = -rxt(k,314)*y(k,49)
         mat(k,1676) = -rxt(k,315)*y(k,49)
         mat(k,1141) = rxt(k,316)*y(k,217)
         mat(k,1210) = .880_r8*rxt(k,338)*y(k,134)
         mat(k,1330) = .500_r8*rxt(k,352)*y(k,134)
         mat(k,1827) = .170_r8*rxt(k,411)*y(k,199) + .050_r8*rxt(k,374)*y(k,206) &
                      + .250_r8*rxt(k,336)*y(k,209) + .170_r8*rxt(k,417)*y(k,212) &
                      + .400_r8*rxt(k,427)*y(k,223) + .250_r8*rxt(k,393)*y(k,225) &
                      + .540_r8*rxt(k,433)*y(k,226) + .510_r8*rxt(k,436)*y(k,228)
         mat(k,1734) = mat(k,1734) + .050_r8*rxt(k,375)*y(k,206) + .250_r8*rxt(k,335) &
                      *y(k,209) + .250_r8*rxt(k,394)*y(k,225)
         mat(k,858) = rxt(k,317)*y(k,217)
         mat(k,2111) = .880_r8*rxt(k,338)*y(k,105) + .500_r8*rxt(k,352)*y(k,111)
         mat(k,1380) = .250_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,2052) = .240_r8*rxt(k,333)*y(k,209) + .500_r8*rxt(k,321)*y(k,220) &
                      + .100_r8*rxt(k,391)*y(k,225)
         mat(k,761) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,203)
         mat(k,1933) = .070_r8*rxt(k,410)*y(k,199) + .070_r8*rxt(k,416)*y(k,212)
         mat(k,1287) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1310) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,197) + .240_r8*rxt(k,333)*y(k,198)
         mat(k,876) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,203)
         mat(k,1676) = mat(k,1676) + rxt(k,316)*y(k,95) + rxt(k,317)*y(k,127)
         mat(k,1131) = .500_r8*rxt(k,321)*y(k,198)
         mat(k,737) = .400_r8*rxt(k,427)*y(k,124)
         mat(k,1174) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
         mat(k,753) = .540_r8*rxt(k,433)*y(k,124)
         mat(k,501) = .510_r8*rxt(k,436)*y(k,124)
         mat(k,684) = -(rxt(k,295)*y(k,217))
         mat(k,1633) = -rxt(k,295)*y(k,50)
         mat(k,1019) = .120_r8*rxt(k,308)*y(k,134)
         mat(k,2088) = .120_r8*rxt(k,308)*y(k,29)
         mat(k,1370) = .100_r8*rxt(k,292)*y(k,198) + .150_r8*rxt(k,293)*y(k,203)
         mat(k,2030) = .100_r8*rxt(k,292)*y(k,197)
         mat(k,1899) = .150_r8*rxt(k,293)*y(k,197) + .150_r8*rxt(k,343)*y(k,211)
         mat(k,1349) = .150_r8*rxt(k,343)*y(k,203)
         mat(k,610) = -(rxt(k,296)*y(k,217))
         mat(k,1625) = -rxt(k,296)*y(k,51)
         mat(k,1369) = .400_r8*rxt(k,293)*y(k,203)
         mat(k,1893) = .400_r8*rxt(k,293)*y(k,197) + .400_r8*rxt(k,343)*y(k,211)
         mat(k,1348) = .400_r8*rxt(k,343)*y(k,203)
         mat(k,786) = -(rxt(k,263)*y(k,217))
         mat(k,1642) = -rxt(k,263)*y(k,52)
         mat(k,1187) = .200_r8*rxt(k,380)*y(k,198)
         mat(k,814) = .300_r8*rxt(k,281)*y(k,198)
         mat(k,2031) = .200_r8*rxt(k,380)*y(k,101) + .300_r8*rxt(k,281)*y(k,194) &
                      + 2.000_r8*rxt(k,260)*y(k,198) + .250_r8*rxt(k,366)*y(k,205) &
                      + .250_r8*rxt(k,371)*y(k,206) + .250_r8*rxt(k,333)*y(k,209) &
                      + .250_r8*rxt(k,445)*y(k,215) + .500_r8*rxt(k,321)*y(k,220) &
                      + .250_r8*rxt(k,450)*y(k,221) + .250_r8*rxt(k,455)*y(k,222) &
                      + .300_r8*rxt(k,391)*y(k,225)
         mat(k,1247) = .250_r8*rxt(k,366)*y(k,198)
         mat(k,1277) = .250_r8*rxt(k,371)*y(k,198)
         mat(k,1305) = .250_r8*rxt(k,333)*y(k,198)
         mat(k,1043) = .250_r8*rxt(k,445)*y(k,198)
         mat(k,1128) = .500_r8*rxt(k,321)*y(k,198)
         mat(k,1109) = .250_r8*rxt(k,450)*y(k,198)
         mat(k,911) = .250_r8*rxt(k,455)*y(k,198)
         mat(k,1167) = .300_r8*rxt(k,391)*y(k,198)
         mat(k,385) = -(rxt(k,264)*y(k,217))
         mat(k,1594) = -rxt(k,264)*y(k,53)
         mat(k,2028) = rxt(k,261)*y(k,203)
         mat(k,1874) = rxt(k,261)*y(k,198)
         mat(k,1430) = -(rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73) + rxt(k,265)*y(k,217) &
                      + (rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,216))
         mat(k,2007) = -rxt(k,176)*y(k,54)
         mat(k,866) = -rxt(k,232)*y(k,54)
         mat(k,1686) = -rxt(k,265)*y(k,54)
         mat(k,1522) = -(rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,54)
         mat(k,1030) = .100_r8*rxt(k,308)*y(k,134)
         mat(k,2120) = .100_r8*rxt(k,308)*y(k,29)
         mat(k,447) = -(rxt(k,228)*y(k,216) + rxt(k,245)*y(k,56) + rxt(k,246)*y(k,217))
         mat(k,1515) = -rxt(k,228)*y(k,55)
         mat(k,1990) = -rxt(k,245)*y(k,55)
         mat(k,1603) = -rxt(k,246)*y(k,55)
         mat(k,2017) = -(rxt(k,175)*y(k,42) + rxt(k,176)*y(k,54) + rxt(k,177)*y(k,77) &
                      + rxt(k,178)*y(k,79) + (rxt(k,179) + rxt(k,180)) * y(k,203) &
                      + rxt(k,181)*y(k,134) + rxt(k,188)*y(k,60) + rxt(k,197)*y(k,92) &
                      + rxt(k,238)*y(k,41) + rxt(k,240)*y(k,43) + rxt(k,243)*y(k,46) &
                      + rxt(k,245)*y(k,55) + rxt(k,286)*y(k,28))
         mat(k,1491) = -rxt(k,175)*y(k,56)
         mat(k,1438) = -rxt(k,176)*y(k,56)
         mat(k,1408) = -rxt(k,177)*y(k,56)
         mat(k,606) = -rxt(k,178)*y(k,56)
         mat(k,1952) = -(rxt(k,179) + rxt(k,180)) * y(k,56)
         mat(k,2130) = -rxt(k,181)*y(k,56)
         mat(k,889) = -rxt(k,188)*y(k,56)
         mat(k,827) = -rxt(k,197)*y(k,56)
         mat(k,471) = -rxt(k,238)*y(k,56)
         mat(k,599) = -rxt(k,240)*y(k,56)
         mat(k,369) = -rxt(k,243)*y(k,56)
         mat(k,451) = -rxt(k,245)*y(k,56)
         mat(k,294) = -rxt(k,286)*y(k,56)
         mat(k,2221) = rxt(k,216)*y(k,59)
         mat(k,101) = 4.000_r8*rxt(k,200)*y(k,216)
         mat(k,145) = rxt(k,201)*y(k,216)
         mat(k,116) = 2.000_r8*rxt(k,202)*y(k,216)
         mat(k,155) = 2.000_r8*rxt(k,203)*y(k,216)
         mat(k,120) = 2.000_r8*rxt(k,204)*y(k,216)
         mat(k,160) = rxt(k,205)*y(k,216)
         mat(k,124) = 2.000_r8*rxt(k,206)*y(k,216)
         mat(k,127) = 3.000_r8*rxt(k,242)*y(k,217)
         mat(k,369) = mat(k,369) + rxt(k,244)*y(k,217)
         mat(k,1978) = rxt(k,216)*y(k,19) + (4.000_r8*rxt(k,183)+2.000_r8*rxt(k,185)) &
                      *y(k,59) + rxt(k,187)*y(k,124) + rxt(k,192)*y(k,133) &
                      + rxt(k,470)*y(k,150) + rxt(k,182)*y(k,198) + rxt(k,193) &
                      *y(k,217)
         mat(k,229) = rxt(k,237)*y(k,216)
         mat(k,225) = rxt(k,252)*y(k,216) + rxt(k,247)*y(k,217)
         mat(k,252) = rxt(k,253)*y(k,216) + rxt(k,248)*y(k,217)
         mat(k,308) = rxt(k,254)*y(k,216) + rxt(k,249)*y(k,217)
         mat(k,2153) = rxt(k,195)*y(k,133) + rxt(k,207)*y(k,216) + rxt(k,196)*y(k,217)
         mat(k,1845) = rxt(k,187)*y(k,59)
         mat(k,2252) = rxt(k,192)*y(k,59) + rxt(k,195)*y(k,85)
         mat(k,1239) = rxt(k,470)*y(k,59)
         mat(k,2069) = rxt(k,182)*y(k,59)
         mat(k,1532) = 4.000_r8*rxt(k,200)*y(k,33) + rxt(k,201)*y(k,34) &
                      + 2.000_r8*rxt(k,202)*y(k,36) + 2.000_r8*rxt(k,203)*y(k,37) &
                      + 2.000_r8*rxt(k,204)*y(k,38) + rxt(k,205)*y(k,39) &
                      + 2.000_r8*rxt(k,206)*y(k,40) + rxt(k,237)*y(k,65) + rxt(k,252) &
                      *y(k,82) + rxt(k,253)*y(k,83) + rxt(k,254)*y(k,84) + rxt(k,207) &
                      *y(k,85)
         mat(k,1696) = 3.000_r8*rxt(k,242)*y(k,44) + rxt(k,244)*y(k,46) + rxt(k,193) &
                      *y(k,59) + rxt(k,247)*y(k,82) + rxt(k,248)*y(k,83) + rxt(k,249) &
                      *y(k,84) + rxt(k,196)*y(k,85)
         mat(k,1986) = rxt(k,188)*y(k,60)
         mat(k,1961) = 2.000_r8*rxt(k,184)*y(k,59)
         mat(k,882) = rxt(k,188)*y(k,56) + (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,2138) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60) + (rxt(k,523) &
                       +rxt(k,529)+rxt(k,534))*y(k,92)
         mat(k,823) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85)
         mat(k,1960) = 2.000_r8*rxt(k,209)*y(k,59)
         mat(k,1977) = -(rxt(k,182)*y(k,198) + (4._r8*rxt(k,183) + 4._r8*rxt(k,184) &
                      + 4._r8*rxt(k,185) + 4._r8*rxt(k,209)) * y(k,59) + rxt(k,186) &
                      *y(k,203) + rxt(k,187)*y(k,124) + rxt(k,189)*y(k,125) + rxt(k,192) &
                      *y(k,133) + (rxt(k,193) + rxt(k,194)) * y(k,217) + (rxt(k,215) &
                      + rxt(k,216) + rxt(k,217)) * y(k,19) + rxt(k,470)*y(k,150))
         mat(k,2068) = -rxt(k,182)*y(k,59)
         mat(k,1951) = -rxt(k,186)*y(k,59)
         mat(k,1844) = -rxt(k,187)*y(k,59)
         mat(k,2196) = -rxt(k,189)*y(k,59)
         mat(k,2251) = -rxt(k,192)*y(k,59)
         mat(k,1695) = -(rxt(k,193) + rxt(k,194)) * y(k,59)
         mat(k,2220) = -(rxt(k,215) + rxt(k,216) + rxt(k,217)) * y(k,59)
         mat(k,1238) = -rxt(k,470)*y(k,59)
         mat(k,2016) = rxt(k,197)*y(k,92) + rxt(k,181)*y(k,134) + rxt(k,180)*y(k,203)
         mat(k,888) = rxt(k,190)*y(k,133)
         mat(k,2152) = rxt(k,208)*y(k,216)
         mat(k,826) = rxt(k,197)*y(k,56) + rxt(k,198)*y(k,133) + rxt(k,199)*y(k,217)
         mat(k,2251) = mat(k,2251) + rxt(k,190)*y(k,60) + rxt(k,198)*y(k,92)
         mat(k,2129) = rxt(k,181)*y(k,56)
         mat(k,331) = rxt(k,475)*y(k,150)
         mat(k,1238) = mat(k,1238) + rxt(k,475)*y(k,136)
         mat(k,1951) = mat(k,1951) + rxt(k,180)*y(k,56)
         mat(k,1531) = rxt(k,208)*y(k,85)
         mat(k,1695) = mat(k,1695) + rxt(k,199)*y(k,92)
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
         mat(k,884) = -(rxt(k,188)*y(k,56) + rxt(k,190)*y(k,133) + rxt(k,191)*y(k,217) &
                      + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85))
         mat(k,1998) = -rxt(k,188)*y(k,60)
         mat(k,2237) = -rxt(k,190)*y(k,60)
         mat(k,1653) = -rxt(k,191)*y(k,60)
         mat(k,2142) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60)
         mat(k,1966) = rxt(k,189)*y(k,125)
         mat(k,2178) = rxt(k,189)*y(k,59)
         mat(k,1103) = -(rxt(k,275)*y(k,217))
         mat(k,1668) = -rxt(k,275)*y(k,62)
         mat(k,1003) = .230_r8*rxt(k,440)*y(k,134)
         mat(k,1414) = rxt(k,211)*y(k,42)
         mat(k,288) = .350_r8*rxt(k,277)*y(k,217)
         mat(k,551) = .630_r8*rxt(k,279)*y(k,134)
         mat(k,1026) = .560_r8*rxt(k,308)*y(k,134)
         mat(k,1479) = rxt(k,211)*y(k,17) + rxt(k,175)*y(k,56) + rxt(k,256)*y(k,126) &
                      + rxt(k,257)*y(k,133) + rxt(k,258)*y(k,217)
         mat(k,366) = rxt(k,243)*y(k,56)
         mat(k,1220) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,217)
         mat(k,2003) = rxt(k,175)*y(k,42) + rxt(k,243)*y(k,46)
         mat(k,980) = rxt(k,302)*y(k,217)
         mat(k,844) = .620_r8*rxt(k,385)*y(k,134)
         mat(k,1208) = .650_r8*rxt(k,338)*y(k,134)
         mat(k,954) = .230_r8*rxt(k,443)*y(k,134)
         mat(k,1328) = .560_r8*rxt(k,352)*y(k,134)
         mat(k,1819) = .170_r8*rxt(k,411)*y(k,199) + .220_r8*rxt(k,336)*y(k,209) &
                      + .400_r8*rxt(k,414)*y(k,210) + .350_r8*rxt(k,417)*y(k,212) &
                      + .225_r8*rxt(k,452)*y(k,221) + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1726) = rxt(k,256)*y(k,42) + rxt(k,314)*y(k,49) + .220_r8*rxt(k,335) &
                      *y(k,209) + .500_r8*rxt(k,394)*y(k,225)
         mat(k,2238) = rxt(k,257)*y(k,42) + rxt(k,464)*y(k,137)
         mat(k,2105) = .230_r8*rxt(k,440)*y(k,6) + .630_r8*rxt(k,279)*y(k,25) &
                      + .560_r8*rxt(k,308)*y(k,29) + .620_r8*rxt(k,385)*y(k,98) &
                      + .650_r8*rxt(k,338)*y(k,105) + .230_r8*rxt(k,443)*y(k,110) &
                      + .560_r8*rxt(k,352)*y(k,111)
         mat(k,360) = rxt(k,464)*y(k,133) + rxt(k,465)*y(k,217)
         mat(k,1060) = .700_r8*rxt(k,461)*y(k,217)
         mat(k,1375) = .220_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,2045) = .110_r8*rxt(k,333)*y(k,209) + .125_r8*rxt(k,450)*y(k,221) &
                      + .200_r8*rxt(k,391)*y(k,225)
         mat(k,760) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,203)
         mat(k,1926) = .070_r8*rxt(k,410)*y(k,199) + .160_r8*rxt(k,413)*y(k,210) &
                      + .140_r8*rxt(k,416)*y(k,212)
         mat(k,1307) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,197) + .110_r8*rxt(k,333)*y(k,198)
         mat(k,723) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,203)
         mat(k,875) = .350_r8*rxt(k,417)*y(k,124) + .140_r8*rxt(k,416)*y(k,203)
         mat(k,1668) = mat(k,1668) + .350_r8*rxt(k,277)*y(k,24) + rxt(k,258)*y(k,42) &
                      + rxt(k,315)*y(k,49) + rxt(k,302)*y(k,75) + rxt(k,465)*y(k,137) &
                      + .700_r8*rxt(k,461)*y(k,178)
         mat(k,1114) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,198)
         mat(k,1171) = .250_r8*rxt(k,393)*y(k,124) + .500_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .200_r8*rxt(k,391)*y(k,198)
         mat(k,992) = .270_r8*rxt(k,440)*y(k,134)
         mat(k,1021) = .200_r8*rxt(k,308)*y(k,134)
         mat(k,685) = rxt(k,295)*y(k,217)
         mat(k,611) = .500_r8*rxt(k,296)*y(k,217)
         mat(k,1102) = rxt(k,275)*y(k,217)
         mat(k,1094) = .800_r8*rxt(k,301)*y(k,217)
         mat(k,978) = rxt(k,302)*y(k,217)
         mat(k,928) = rxt(k,267)*y(k,217)
         mat(k,573) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,943) = .270_r8*rxt(k,443)*y(k,134)
         mat(k,1324) = .100_r8*rxt(k,352)*y(k,134)
         mat(k,1804) = rxt(k,294)*y(k,197) + .900_r8*rxt(k,452)*y(k,221)
         mat(k,2090) = .270_r8*rxt(k,440)*y(k,6) + .200_r8*rxt(k,308)*y(k,29) &
                      + .270_r8*rxt(k,443)*y(k,110) + .100_r8*rxt(k,352)*y(k,111)
         mat(k,1057) = 1.800_r8*rxt(k,461)*y(k,217)
         mat(k,1371) = rxt(k,294)*y(k,124) + 4.000_r8*rxt(k,291)*y(k,197) &
                      + .900_r8*rxt(k,292)*y(k,198) + rxt(k,365)*y(k,205) &
                      + 2.000_r8*rxt(k,341)*y(k,211) + rxt(k,390)*y(k,225)
         mat(k,2033) = .900_r8*rxt(k,292)*y(k,197) + rxt(k,342)*y(k,211) &
                      + .500_r8*rxt(k,450)*y(k,221)
         mat(k,1910) = .450_r8*rxt(k,343)*y(k,211)
         mat(k,1248) = rxt(k,365)*y(k,197)
         mat(k,1350) = 2.000_r8*rxt(k,341)*y(k,197) + rxt(k,342)*y(k,198) &
                      + .450_r8*rxt(k,343)*y(k,203) + 4.000_r8*rxt(k,344)*y(k,211)
         mat(k,1644) = rxt(k,295)*y(k,50) + .500_r8*rxt(k,296)*y(k,51) + rxt(k,275) &
                      *y(k,62) + .800_r8*rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) &
                      + rxt(k,267)*y(k,87) + .500_r8*rxt(k,351)*y(k,109) &
                      + 1.800_r8*rxt(k,461)*y(k,178)
         mat(k,1110) = .900_r8*rxt(k,452)*y(k,124) + .500_r8*rxt(k,450)*y(k,198)
         mat(k,1168) = rxt(k,390)*y(k,197)
         mat(k,235) = -(rxt(k,236)*y(k,216))
         mat(k,1512) = -rxt(k,236)*y(k,64)
         mat(k,142) = rxt(k,201)*y(k,216)
         mat(k,147) = rxt(k,227)*y(k,216)
         mat(k,153) = rxt(k,203)*y(k,216)
         mat(k,118) = 2.000_r8*rxt(k,204)*y(k,216)
         mat(k,157) = 2.000_r8*rxt(k,205)*y(k,216)
         mat(k,122) = rxt(k,206)*y(k,216)
         mat(k,106) = 2.000_r8*rxt(k,229)*y(k,216)
         mat(k,247) = rxt(k,253)*y(k,216) + rxt(k,248)*y(k,217)
         mat(k,303) = rxt(k,254)*y(k,216) + rxt(k,249)*y(k,217)
         mat(k,1512) = mat(k,1512) + rxt(k,201)*y(k,34) + rxt(k,227)*y(k,35) &
                      + rxt(k,203)*y(k,37) + 2.000_r8*rxt(k,204)*y(k,38) &
                      + 2.000_r8*rxt(k,205)*y(k,39) + rxt(k,206)*y(k,40) &
                      + 2.000_r8*rxt(k,229)*y(k,78) + rxt(k,253)*y(k,83) + rxt(k,254) &
                      *y(k,84)
         mat(k,1572) = rxt(k,248)*y(k,83) + rxt(k,249)*y(k,84)
         mat(k,226) = -(rxt(k,237)*y(k,216))
         mat(k,1511) = -rxt(k,237)*y(k,65)
         mat(k,114) = rxt(k,202)*y(k,216)
         mat(k,152) = rxt(k,203)*y(k,216)
         mat(k,222) = rxt(k,252)*y(k,216) + rxt(k,247)*y(k,217)
         mat(k,1511) = mat(k,1511) + rxt(k,202)*y(k,36) + rxt(k,203)*y(k,37) &
                      + rxt(k,252)*y(k,82)
         mat(k,1570) = rxt(k,247)*y(k,82)
         mat(k,194) = -(rxt(k,409)*y(k,217))
         mat(k,1564) = -rxt(k,409)*y(k,66)
         mat(k,188) = .180_r8*rxt(k,429)*y(k,217)
         mat(k,1564) = mat(k,1564) + .180_r8*rxt(k,429)*y(k,180)
         mat(k,297) = -(rxt(k,462)*y(k,126) + (rxt(k,463) + rxt(k,477)) * y(k,217))
         mat(k,1707) = -rxt(k,462)*y(k,67)
         mat(k,1582) = -(rxt(k,463) + rxt(k,477)) * y(k,67)
         mat(k,712) = rxt(k,297)*y(k,203)
         mat(k,1865) = rxt(k,297)*y(k,202)
         mat(k,864) = -(rxt(k,232)*y(k,54) + rxt(k,233)*y(k,77) + rxt(k,234)*y(k,229) &
                      + rxt(k,235)*y(k,89))
         mat(k,1427) = -rxt(k,232)*y(k,73)
         mat(k,1400) = -rxt(k,233)*y(k,73)
         mat(k,2264) = -rxt(k,234)*y(k,73)
         mat(k,1444) = -rxt(k,235)*y(k,73)
         mat(k,148) = rxt(k,227)*y(k,216)
         mat(k,158) = rxt(k,205)*y(k,216)
         mat(k,236) = 2.000_r8*rxt(k,236)*y(k,216)
         mat(k,227) = rxt(k,237)*y(k,216)
         mat(k,1519) = rxt(k,227)*y(k,35) + rxt(k,205)*y(k,39) + 2.000_r8*rxt(k,236) &
                      *y(k,64) + rxt(k,237)*y(k,65)
         mat(k,1096) = -(rxt(k,301)*y(k,217))
         mat(k,1667) = -rxt(k,301)*y(k,74)
         mat(k,582) = .700_r8*rxt(k,376)*y(k,217)
         mat(k,558) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,375) = rxt(k,388)*y(k,217)
         mat(k,1818) = .050_r8*rxt(k,374)*y(k,206) + .530_r8*rxt(k,336)*y(k,209) &
                      + .225_r8*rxt(k,452)*y(k,221) + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1725) = .050_r8*rxt(k,375)*y(k,206) + .530_r8*rxt(k,335)*y(k,209) &
                      + .250_r8*rxt(k,394)*y(k,225)
         mat(k,1374) = .530_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,2044) = .260_r8*rxt(k,333)*y(k,209) + .125_r8*rxt(k,450)*y(k,221) &
                      + .100_r8*rxt(k,391)*y(k,225)
         mat(k,1281) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1306) = .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335)*y(k,126) &
                      + .530_r8*rxt(k,332)*y(k,197) + .260_r8*rxt(k,333)*y(k,198)
         mat(k,1667) = mat(k,1667) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + rxt(k,388)*y(k,115)
         mat(k,1113) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,198)
         mat(k,1170) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
         mat(k,979) = -(rxt(k,302)*y(k,217))
         mat(k,1660) = -rxt(k,302)*y(k,75)
         mat(k,287) = .650_r8*rxt(k,277)*y(k,217)
         mat(k,1095) = .200_r8*rxt(k,301)*y(k,217)
         mat(k,1072) = rxt(k,389)*y(k,217)
         mat(k,1813) = rxt(k,400)*y(k,191) + .050_r8*rxt(k,374)*y(k,206) &
                      + .400_r8*rxt(k,414)*y(k,210) + .170_r8*rxt(k,417)*y(k,212) &
                      + .700_r8*rxt(k,420)*y(k,218) + .600_r8*rxt(k,427)*y(k,223) &
                      + .250_r8*rxt(k,393)*y(k,225) + .340_r8*rxt(k,433)*y(k,226) &
                      + .170_r8*rxt(k,436)*y(k,228)
         mat(k,1718) = .050_r8*rxt(k,375)*y(k,206) + .250_r8*rxt(k,394)*y(k,225)
         mat(k,485) = rxt(k,400)*y(k,124)
         mat(k,1372) = .250_r8*rxt(k,390)*y(k,225)
         mat(k,2039) = .100_r8*rxt(k,391)*y(k,225)
         mat(k,1921) = .160_r8*rxt(k,413)*y(k,210) + .070_r8*rxt(k,416)*y(k,212)
         mat(k,1280) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,722) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,203)
         mat(k,874) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,203)
         mat(k,1660) = mat(k,1660) + .650_r8*rxt(k,277)*y(k,24) + .200_r8*rxt(k,301) &
                      *y(k,74) + rxt(k,389)*y(k,116)
         mat(k,455) = .700_r8*rxt(k,420)*y(k,124)
         mat(k,735) = .600_r8*rxt(k,427)*y(k,124)
         mat(k,1169) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
         mat(k,751) = .340_r8*rxt(k,433)*y(k,124)
         mat(k,500) = .170_r8*rxt(k,436)*y(k,124)
         mat(k,1463) = -((rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,203) + rxt(k,141) &
                      *y(k,134))
         mat(k,1944) = -(rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,76)
         mat(k,2122) = -rxt(k,141)*y(k,76)
         mat(k,1484) = rxt(k,258)*y(k,217)
         mat(k,1432) = rxt(k,272)*y(k,216)
         mat(k,2009) = rxt(k,177)*y(k,77)
         mat(k,868) = rxt(k,233)*y(k,77)
         mat(k,1404) = rxt(k,177)*y(k,56) + rxt(k,233)*y(k,73) + rxt(k,133)*y(k,133) &
                      + rxt(k,125)*y(k,216) + rxt(k,142)*y(k,217)
         mat(k,806) = rxt(k,231)*y(k,216)
         mat(k,2145) = rxt(k,208)*y(k,216)
         mat(k,492) = rxt(k,163)*y(k,217)
         mat(k,2244) = rxt(k,133)*y(k,77) + rxt(k,145)*y(k,217)
         mat(k,362) = rxt(k,465)*y(k,217)
         mat(k,513) = rxt(k,471)*y(k,217)
         mat(k,1234) = rxt(k,476)*y(k,217)
         mat(k,1524) = rxt(k,272)*y(k,54) + rxt(k,125)*y(k,77) + rxt(k,231)*y(k,81) &
                      + rxt(k,208)*y(k,85)
         mat(k,1688) = rxt(k,258)*y(k,42) + rxt(k,142)*y(k,77) + rxt(k,163)*y(k,112) &
                      + rxt(k,145)*y(k,133) + rxt(k,465)*y(k,137) + rxt(k,471) &
                      *y(k,148) + rxt(k,476)*y(k,150)
         mat(k,1401) = -(rxt(k,125)*y(k,216) + rxt(k,133)*y(k,133) + rxt(k,142) &
                      *y(k,217) + rxt(k,177)*y(k,56) + rxt(k,233)*y(k,73))
         mat(k,1520) = -rxt(k,125)*y(k,77)
         mat(k,2240) = -rxt(k,133)*y(k,77)
         mat(k,1684) = -rxt(k,142)*y(k,77)
         mat(k,2005) = -rxt(k,177)*y(k,77)
         mat(k,865) = -rxt(k,233)*y(k,77)
         mat(k,1429) = rxt(k,273)*y(k,216)
         mat(k,1460) = rxt(k,135)*y(k,203)
         mat(k,1940) = rxt(k,135)*y(k,76)
         mat(k,1520) = mat(k,1520) + rxt(k,273)*y(k,54)
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
         mat(k,105) = -(rxt(k,229)*y(k,216))
         mat(k,1500) = -rxt(k,229)*y(k,78)
         mat(k,603) = -(rxt(k,134)*y(k,133) + rxt(k,143)*y(k,217) + rxt(k,178)*y(k,56))
         mat(k,2232) = -rxt(k,134)*y(k,79)
         mat(k,1624) = -rxt(k,143)*y(k,79)
         mat(k,1994) = -rxt(k,178)*y(k,79)
         mat(k,1892) = 2.000_r8*rxt(k,149)*y(k,203)
         mat(k,1624) = mat(k,1624) + 2.000_r8*rxt(k,148)*y(k,217)
         mat(k,256) = rxt(k,478)*y(k,229)
         mat(k,2260) = rxt(k,478)*y(k,152)
         mat(k,804) = -(rxt(k,224)*y(k,133) + rxt(k,225)*y(k,217) + (rxt(k,230) &
                      + rxt(k,231)) * y(k,216))
         mat(k,2234) = -rxt(k,224)*y(k,81)
         mat(k,1645) = -rxt(k,225)*y(k,81)
         mat(k,1518) = -(rxt(k,230) + rxt(k,231)) * y(k,81)
         mat(k,1413) = rxt(k,211)*y(k,42) + rxt(k,212)*y(k,203)
         mat(k,1477) = rxt(k,211)*y(k,17)
         mat(k,1911) = rxt(k,212)*y(k,17)
         mat(k,221) = -(rxt(k,247)*y(k,217) + rxt(k,252)*y(k,216))
         mat(k,1569) = -rxt(k,247)*y(k,82)
         mat(k,1510) = -rxt(k,252)*y(k,82)
         mat(k,248) = -(rxt(k,248)*y(k,217) + rxt(k,253)*y(k,216))
         mat(k,1575) = -rxt(k,248)*y(k,83)
         mat(k,1513) = -rxt(k,253)*y(k,83)
         mat(k,304) = -(rxt(k,249)*y(k,217) + rxt(k,254)*y(k,216))
         mat(k,1583) = -rxt(k,249)*y(k,84)
         mat(k,1514) = -rxt(k,254)*y(k,84)
         mat(k,2156) = -(rxt(k,195)*y(k,133) + rxt(k,196)*y(k,217) + (rxt(k,207) &
                      + rxt(k,208)) * y(k,216) + (rxt(k,523) + rxt(k,529) + rxt(k,534) &
                      ) * y(k,92) + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60) &
                      + (rxt(k,530) + rxt(k,535)) * y(k,91))
         mat(k,2255) = -rxt(k,195)*y(k,85)
         mat(k,1699) = -rxt(k,196)*y(k,85)
         mat(k,1535) = -(rxt(k,207) + rxt(k,208)) * y(k,85)
         mat(k,828) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85)
         mat(k,890) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85)
         mat(k,782) = -(rxt(k,530) + rxt(k,535)) * y(k,85)
         mat(k,295) = rxt(k,286)*y(k,56)
         mat(k,472) = rxt(k,238)*y(k,56)
         mat(k,1494) = rxt(k,175)*y(k,56)
         mat(k,601) = rxt(k,240)*y(k,56)
         mat(k,371) = 2.000_r8*rxt(k,243)*y(k,56)
         mat(k,1440) = rxt(k,176)*y(k,56)
         mat(k,452) = rxt(k,245)*y(k,56)
         mat(k,2020) = rxt(k,286)*y(k,28) + rxt(k,238)*y(k,41) + rxt(k,175)*y(k,42) &
                      + rxt(k,240)*y(k,43) + 2.000_r8*rxt(k,243)*y(k,46) + rxt(k,176) &
                      *y(k,54) + rxt(k,245)*y(k,55) + rxt(k,177)*y(k,77) + rxt(k,178) &
                      *y(k,79) + rxt(k,197)*y(k,92) + rxt(k,179)*y(k,203)
         mat(k,1981) = rxt(k,194)*y(k,217)
         mat(k,1410) = rxt(k,177)*y(k,56)
         mat(k,607) = rxt(k,178)*y(k,56)
         mat(k,828) = mat(k,828) + rxt(k,197)*y(k,56)
         mat(k,1955) = rxt(k,179)*y(k,56)
         mat(k,1699) = mat(k,1699) + rxt(k,194)*y(k,59)
         mat(k,179) = -(rxt(k,266)*y(k,217) + rxt(k,274)*y(k,216))
         mat(k,1562) = -rxt(k,266)*y(k,86)
         mat(k,1508) = -rxt(k,274)*y(k,86)
         mat(k,929) = -(rxt(k,267)*y(k,217))
         mat(k,1657) = -rxt(k,267)*y(k,87)
         mat(k,996) = .050_r8*rxt(k,440)*y(k,134)
         mat(k,286) = .350_r8*rxt(k,277)*y(k,217)
         mat(k,550) = .370_r8*rxt(k,279)*y(k,134)
         mat(k,1023) = .120_r8*rxt(k,308)*y(k,134)
         mat(k,842) = .110_r8*rxt(k,385)*y(k,134)
         mat(k,1207) = .330_r8*rxt(k,338)*y(k,134)
         mat(k,947) = .050_r8*rxt(k,443)*y(k,134)
         mat(k,1325) = .120_r8*rxt(k,352)*y(k,134)
         mat(k,1811) = rxt(k,270)*y(k,204)
         mat(k,2095) = .050_r8*rxt(k,440)*y(k,6) + .370_r8*rxt(k,279)*y(k,25) &
                      + .120_r8*rxt(k,308)*y(k,29) + .110_r8*rxt(k,385)*y(k,98) &
                      + .330_r8*rxt(k,338)*y(k,105) + .050_r8*rxt(k,443)*y(k,110) &
                      + .120_r8*rxt(k,352)*y(k,111)
         mat(k,1919) = rxt(k,268)*y(k,204)
         mat(k,442) = rxt(k,270)*y(k,124) + rxt(k,268)*y(k,203)
         mat(k,1657) = mat(k,1657) + .350_r8*rxt(k,277)*y(k,24)
         mat(k,1425) = rxt(k,232)*y(k,73)
         mat(k,863) = rxt(k,232)*y(k,54) + rxt(k,233)*y(k,77) + rxt(k,235)*y(k,89) &
                      + rxt(k,234)*y(k,229)
         mat(k,1399) = rxt(k,233)*y(k,73)
         mat(k,1443) = rxt(k,235)*y(k,73)
         mat(k,2262) = rxt(k,234)*y(k,73)
         mat(k,1447) = -(rxt(k,172)*y(k,217) + rxt(k,235)*y(k,73))
         mat(k,1687) = -rxt(k,172)*y(k,89)
         mat(k,867) = -rxt(k,235)*y(k,89)
         mat(k,1483) = rxt(k,256)*y(k,126)
         mat(k,1088) = rxt(k,288)*y(k,126)
         mat(k,1223) = rxt(k,314)*y(k,126)
         mat(k,885) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,299) = rxt(k,462)*y(k,126)
         mat(k,2144) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60)
         mat(k,2188) = rxt(k,171)*y(k,217)
         mat(k,1744) = rxt(k,256)*y(k,42) + rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) &
                      + rxt(k,462)*y(k,67)
         mat(k,1687) = mat(k,1687) + rxt(k,171)*y(k,125)
         mat(k,403) = -(rxt(k,150)*y(k,217))
         mat(k,1597) = -rxt(k,150)*y(k,90)
         mat(k,2164) = rxt(k,169)*y(k,203)
         mat(k,1877) = rxt(k,169)*y(k,125)
         mat(k,778) = -(rxt(k,226)*y(k,133) + (rxt(k,530) + rxt(k,535)) * y(k,85))
         mat(k,2233) = -rxt(k,226)*y(k,91)
         mat(k,2140) = -(rxt(k,530) + rxt(k,535)) * y(k,91)
         mat(k,2208) = rxt(k,218)*y(k,203)
         mat(k,1908) = rxt(k,218)*y(k,19)
         mat(k,824) = -(rxt(k,197)*y(k,56) + rxt(k,198)*y(k,133) + rxt(k,199)*y(k,217) &
                      + (rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85))
         mat(k,1997) = -rxt(k,197)*y(k,92)
         mat(k,2235) = -rxt(k,198)*y(k,92)
         mat(k,1647) = -rxt(k,199)*y(k,92)
         mat(k,2141) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,92)
         mat(k,1964) = rxt(k,186)*y(k,203)
         mat(k,883) = rxt(k,191)*y(k,217)
         mat(k,1913) = rxt(k,186)*y(k,59)
         mat(k,1647) = mat(k,1647) + rxt(k,191)*y(k,60)
         mat(k,1153) = -(rxt(k,331)*y(k,217))
         mat(k,1672) = -rxt(k,331)*y(k,93)
         mat(k,584) = .300_r8*rxt(k,376)*y(k,217)
         mat(k,560) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,1823) = rxt(k,330)*y(k,200) + rxt(k,337)*y(k,209)
         mat(k,567) = rxt(k,330)*y(k,124)
         mat(k,1309) = rxt(k,337)*y(k,124)
         mat(k,1672) = mat(k,1672) + .300_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100)
         mat(k,230) = -(rxt(k,362)*y(k,217))
         mat(k,1571) = -rxt(k,362)*y(k,94)
         mat(k,1140) = -(rxt(k,316)*y(k,217))
         mat(k,1671) = -rxt(k,316)*y(k,95)
         mat(k,583) = .700_r8*rxt(k,376)*y(k,217)
         mat(k,559) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,574) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,1822) = .050_r8*rxt(k,374)*y(k,206) + .220_r8*rxt(k,336)*y(k,209) &
                      + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1729) = .050_r8*rxt(k,375)*y(k,206) + .220_r8*rxt(k,335)*y(k,209) &
                      + .250_r8*rxt(k,394)*y(k,225)
         mat(k,535) = .500_r8*rxt(k,320)*y(k,217)
         mat(k,1376) = .220_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,2048) = .230_r8*rxt(k,333)*y(k,209) + .200_r8*rxt(k,321)*y(k,220) &
                      + .100_r8*rxt(k,391)*y(k,225)
         mat(k,1283) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1308) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,197) + .230_r8*rxt(k,333)*y(k,198)
         mat(k,1671) = mat(k,1671) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + .500_r8*rxt(k,351)*y(k,109) + .500_r8*rxt(k,320) &
                      *y(k,146)
         mat(k,1130) = .200_r8*rxt(k,321)*y(k,198)
         mat(k,1172) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
         mat(k,320) = -(rxt(k,363)*y(k,217))
         mat(k,1586) = -rxt(k,363)*y(k,96)
         mat(k,1777) = .870_r8*rxt(k,374)*y(k,206)
         mat(k,1708) = .950_r8*rxt(k,375)*y(k,206)
         mat(k,1367) = rxt(k,370)*y(k,206)
         mat(k,2026) = .750_r8*rxt(k,371)*y(k,206)
         mat(k,1273) = .870_r8*rxt(k,374)*y(k,124) + .950_r8*rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,197) + .750_r8*rxt(k,371)*y(k,198)
         mat(k,135) = -(rxt(k,364)*y(k,217))
         mat(k,1558) = -rxt(k,364)*y(k,97)
         mat(k,689) = .600_r8*rxt(k,387)*y(k,217)
         mat(k,1558) = mat(k,1558) + .600_r8*rxt(k,387)*y(k,103)
         mat(k,841) = -(rxt(k,378)*y(k,126) + rxt(k,385)*y(k,134) + rxt(k,386) &
                      *y(k,217))
         mat(k,1712) = -rxt(k,378)*y(k,98)
         mat(k,2092) = -rxt(k,385)*y(k,98)
         mat(k,1649) = -rxt(k,386)*y(k,98)
         mat(k,581) = -(rxt(k,376)*y(k,217))
         mat(k,1621) = -rxt(k,376)*y(k,99)
         mat(k,1791) = .080_r8*rxt(k,368)*y(k,205)
         mat(k,1245) = .080_r8*rxt(k,368)*y(k,124)
         mat(k,556) = -(rxt(k,377)*y(k,217))
         mat(k,1618) = -rxt(k,377)*y(k,100)
         mat(k,1789) = .080_r8*rxt(k,374)*y(k,206)
         mat(k,1274) = .080_r8*rxt(k,374)*y(k,124)
         mat(k,1193) = -(rxt(k,379)*y(k,197) + rxt(k,380)*y(k,198) + rxt(k,381) &
                      *y(k,203) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126))
         mat(k,1378) = -rxt(k,379)*y(k,101)
         mat(k,2050) = -rxt(k,380)*y(k,101)
         mat(k,1931) = -rxt(k,381)*y(k,101)
         mat(k,1825) = -rxt(k,382)*y(k,101)
         mat(k,1732) = -rxt(k,383)*y(k,101)
         mat(k,845) = rxt(k,378)*y(k,126)
         mat(k,1732) = mat(k,1732) + rxt(k,378)*y(k,98)
         mat(k,397) = -(rxt(k,384)*y(k,217))
         mat(k,1596) = -rxt(k,384)*y(k,102)
         mat(k,1185) = rxt(k,381)*y(k,203)
         mat(k,1876) = rxt(k,381)*y(k,101)
         mat(k,690) = -(rxt(k,387)*y(k,217))
         mat(k,1634) = -rxt(k,387)*y(k,103)
         mat(k,1900) = rxt(k,367)*y(k,205) + rxt(k,372)*y(k,206)
         mat(k,1246) = rxt(k,367)*y(k,203)
         mat(k,1276) = rxt(k,372)*y(k,203)
         mat(k,74) = -(rxt(k,509)*y(k,217))
         mat(k,1550) = -rxt(k,509)*y(k,104)
         mat(k,1209) = -(rxt(k,338)*y(k,134) + rxt(k,339)*y(k,217))
         mat(k,2110) = -rxt(k,338)*y(k,105)
         mat(k,1675) = -rxt(k,339)*y(k,105)
         mat(k,846) = .300_r8*rxt(k,385)*y(k,134)
         mat(k,1826) = .360_r8*rxt(k,368)*y(k,205)
         mat(k,1733) = .400_r8*rxt(k,369)*y(k,205)
         mat(k,2110) = mat(k,2110) + .300_r8*rxt(k,385)*y(k,98)
         mat(k,1379) = .390_r8*rxt(k,365)*y(k,205)
         mat(k,2051) = .310_r8*rxt(k,366)*y(k,205)
         mat(k,1254) = .360_r8*rxt(k,368)*y(k,124) + .400_r8*rxt(k,369)*y(k,126) &
                      + .390_r8*rxt(k,365)*y(k,197) + .310_r8*rxt(k,366)*y(k,198)
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
         mat(k,310) = -(rxt(k,340)*y(k,217))
         mat(k,1584) = -rxt(k,340)*y(k,106)
         mat(k,1869) = rxt(k,334)*y(k,209)
         mat(k,1304) = rxt(k,334)*y(k,203)
         mat(k,506) = -(rxt(k,349)*y(k,217))
         mat(k,1612) = -rxt(k,349)*y(k,107)
         mat(k,1787) = .800_r8*rxt(k,358)*y(k,189)
         mat(k,894) = .800_r8*rxt(k,358)*y(k,124)
         mat(k,315) = -(rxt(k,350)*y(k,217))
         mat(k,1585) = -rxt(k,350)*y(k,108)
         mat(k,1870) = .800_r8*rxt(k,347)*y(k,213)
         mat(k,676) = .800_r8*rxt(k,347)*y(k,203)
         mat(k,572) = -(rxt(k,351)*y(k,217))
         mat(k,1620) = -rxt(k,351)*y(k,109)
         mat(k,2170) = rxt(k,354)*y(k,211)
         mat(k,1347) = rxt(k,354)*y(k,125)
         mat(k,948) = -(rxt(k,442)*y(k,126) + rxt(k,443)*y(k,134) + rxt(k,444) &
                      *y(k,217))
         mat(k,1716) = -rxt(k,442)*y(k,110)
         mat(k,2096) = -rxt(k,443)*y(k,110)
         mat(k,1658) = -rxt(k,444)*y(k,110)
         mat(k,1332) = -(rxt(k,352)*y(k,134) + rxt(k,353)*y(k,217))
         mat(k,2116) = -rxt(k,352)*y(k,111)
         mat(k,1681) = -rxt(k,353)*y(k,111)
         mat(k,849) = .200_r8*rxt(k,385)*y(k,134)
         mat(k,1831) = .560_r8*rxt(k,368)*y(k,205)
         mat(k,1739) = .600_r8*rxt(k,369)*y(k,205)
         mat(k,2116) = mat(k,2116) + .200_r8*rxt(k,385)*y(k,98)
         mat(k,1384) = .610_r8*rxt(k,365)*y(k,205)
         mat(k,2056) = .440_r8*rxt(k,366)*y(k,205)
         mat(k,1258) = .560_r8*rxt(k,368)*y(k,124) + .600_r8*rxt(k,369)*y(k,126) &
                      + .610_r8*rxt(k,365)*y(k,197) + .440_r8*rxt(k,366)*y(k,198)
         mat(k,491) = -(rxt(k,151)*y(k,124) + (rxt(k,152) + rxt(k,153) + rxt(k,154) &
                      ) * y(k,125) + rxt(k,163)*y(k,217))
         mat(k,1785) = -rxt(k,151)*y(k,112)
         mat(k,2166) = -(rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,112)
         mat(k,1610) = -rxt(k,163)*y(k,112)
         mat(k,183) = -((rxt(k,167) + rxt(k,168)) * y(k,216))
         mat(k,1509) = -(rxt(k,167) + rxt(k,168)) * y(k,113)
         mat(k,490) = rxt(k,152)*y(k,125)
         mat(k,2162) = rxt(k,152)*y(k,112)
         mat(k,2163) = rxt(k,170)*y(k,126)
         mat(k,1706) = rxt(k,170)*y(k,125)
         mat(k,373) = -(rxt(k,388)*y(k,217))
         mat(k,1593) = -rxt(k,388)*y(k,115)
         mat(k,1184) = .200_r8*rxt(k,380)*y(k,198)
         mat(k,2027) = .200_r8*rxt(k,380)*y(k,101)
         mat(k,1073) = -(rxt(k,389)*y(k,217))
         mat(k,1665) = -rxt(k,389)*y(k,116)
         mat(k,1189) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) + rxt(k,379)*y(k,197) &
                      + .800_r8*rxt(k,380)*y(k,198)
         mat(k,1816) = rxt(k,382)*y(k,101)
         mat(k,1723) = rxt(k,383)*y(k,101)
         mat(k,1373) = rxt(k,379)*y(k,101)
         mat(k,2042) = .800_r8*rxt(k,380)*y(k,101)
         mat(k,96) = -(rxt(k,479)*y(k,217))
         mat(k,1554) = -rxt(k,479)*y(k,120)
         mat(k,1842) = -(rxt(k,151)*y(k,112) + rxt(k,160)*y(k,126) + rxt(k,164) &
                      *y(k,203) + rxt(k,165)*y(k,134) + rxt(k,166)*y(k,133) + rxt(k,187) &
                      *y(k,59) + rxt(k,219)*y(k,19) + rxt(k,262)*y(k,198) + rxt(k,270) &
                      *y(k,204) + rxt(k,283)*y(k,194) + rxt(k,294)*y(k,197) + rxt(k,298) &
                      *y(k,202) + rxt(k,311)*y(k,195) + rxt(k,319)*y(k,219) + rxt(k,323) &
                      *y(k,220) + (rxt(k,329) + rxt(k,330)) * y(k,200) + (rxt(k,336) &
                      + rxt(k,337)) * y(k,209) + rxt(k,345)*y(k,211) + rxt(k,348) &
                      *y(k,213) + (rxt(k,358) + rxt(k,359)) * y(k,189) + rxt(k,368) &
                      *y(k,205) + rxt(k,374)*y(k,206) + rxt(k,382)*y(k,101) + rxt(k,393) &
                      *y(k,225) + rxt(k,397)*y(k,188) + rxt(k,400)*y(k,191) + rxt(k,405) &
                      *y(k,193) + rxt(k,407)*y(k,196) + rxt(k,411)*y(k,199) + rxt(k,414) &
                      *y(k,210) + rxt(k,417)*y(k,212) + rxt(k,420)*y(k,218) + rxt(k,427) &
                      *y(k,223) + rxt(k,433)*y(k,226) + rxt(k,436)*y(k,228) + rxt(k,447) &
                      *y(k,215) + rxt(k,452)*y(k,221) + rxt(k,457)*y(k,222))
         mat(k,495) = -rxt(k,151)*y(k,124)
         mat(k,1750) = -rxt(k,160)*y(k,124)
         mat(k,1949) = -rxt(k,164)*y(k,124)
         mat(k,2127) = -rxt(k,165)*y(k,124)
         mat(k,2249) = -rxt(k,166)*y(k,124)
         mat(k,1975) = -rxt(k,187)*y(k,124)
         mat(k,2218) = -rxt(k,219)*y(k,124)
         mat(k,2066) = -rxt(k,262)*y(k,124)
         mat(k,444) = -rxt(k,270)*y(k,124)
         mat(k,819) = -rxt(k,283)*y(k,124)
         mat(k,1392) = -rxt(k,294)*y(k,124)
         mat(k,718) = -rxt(k,298)*y(k,124)
         mat(k,796) = -rxt(k,311)*y(k,124)
         mat(k,773) = -rxt(k,319)*y(k,124)
         mat(k,1135) = -rxt(k,323)*y(k,124)
         mat(k,569) = -(rxt(k,329) + rxt(k,330)) * y(k,124)
         mat(k,1318) = -(rxt(k,336) + rxt(k,337)) * y(k,124)
         mat(k,1360) = -rxt(k,345)*y(k,124)
         mat(k,681) = -rxt(k,348)*y(k,124)
         mat(k,905) = -(rxt(k,358) + rxt(k,359)) * y(k,124)
         mat(k,1265) = -rxt(k,368)*y(k,124)
         mat(k,1297) = -rxt(k,374)*y(k,124)
         mat(k,1202) = -rxt(k,382)*y(k,124)
         mat(k,1179) = -rxt(k,393)*y(k,124)
         mat(k,521) = -rxt(k,397)*y(k,124)
         mat(k,487) = -rxt(k,400)*y(k,124)
         mat(k,438) = -rxt(k,405)*y(k,124)
         mat(k,627) = -rxt(k,407)*y(k,124)
         mat(k,763) = -rxt(k,411)*y(k,124)
         mat(k,724) = -rxt(k,414)*y(k,124)
         mat(k,878) = -rxt(k,417)*y(k,124)
         mat(k,457) = -rxt(k,420)*y(k,124)
         mat(k,739) = -rxt(k,427)*y(k,124)
         mat(k,756) = -rxt(k,433)*y(k,124)
         mat(k,503) = -rxt(k,436)*y(k,124)
         mat(k,1053) = -rxt(k,447)*y(k,124)
         mat(k,1121) = -rxt(k,452)*y(k,124)
         mat(k,918) = -rxt(k,457)*y(k,124)
         mat(k,495) = mat(k,495) + 2.000_r8*rxt(k,153)*y(k,125) + rxt(k,163)*y(k,217)
         mat(k,185) = 2.000_r8*rxt(k,167)*y(k,216)
         mat(k,2194) = 2.000_r8*rxt(k,153)*y(k,112) + rxt(k,156)*y(k,133) + rxt(k,472) &
                      *y(k,150)
         mat(k,2249) = mat(k,2249) + rxt(k,156)*y(k,125)
         mat(k,1236) = rxt(k,472)*y(k,125)
         mat(k,1529) = 2.000_r8*rxt(k,167)*y(k,113)
         mat(k,1693) = rxt(k,163)*y(k,112)
         mat(k,2201) = -((rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,112) + (rxt(k,156) &
                      + rxt(k,158)) * y(k,133) + rxt(k,157)*y(k,134) + rxt(k,169) &
                      *y(k,203) + rxt(k,170)*y(k,126) + rxt(k,171)*y(k,217) + rxt(k,189) &
                      *y(k,59) + rxt(k,220)*y(k,19) + rxt(k,305)*y(k,197) + rxt(k,354) &
                      *y(k,211) + rxt(k,412)*y(k,199) + rxt(k,415)*y(k,210) + rxt(k,418) &
                      *y(k,212) + rxt(k,422)*y(k,141) + rxt(k,425)*y(k,188) + rxt(k,472) &
                      *y(k,150))
         mat(k,496) = -(rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,125)
         mat(k,2256) = -(rxt(k,156) + rxt(k,158)) * y(k,125)
         mat(k,2134) = -rxt(k,157)*y(k,125)
         mat(k,1956) = -rxt(k,169)*y(k,125)
         mat(k,1757) = -rxt(k,170)*y(k,125)
         mat(k,1700) = -rxt(k,171)*y(k,125)
         mat(k,1982) = -rxt(k,189)*y(k,125)
         mat(k,2225) = -rxt(k,220)*y(k,125)
         mat(k,1396) = -rxt(k,305)*y(k,125)
         mat(k,1364) = -rxt(k,354)*y(k,125)
         mat(k,766) = -rxt(k,412)*y(k,125)
         mat(k,726) = -rxt(k,415)*y(k,125)
         mat(k,881) = -rxt(k,418)*y(k,125)
         mat(k,466) = -rxt(k,422)*y(k,125)
         mat(k,523) = -rxt(k,425)*y(k,125)
         mat(k,1241) = -rxt(k,472)*y(k,125)
         mat(k,675) = rxt(k,356)*y(k,217)
         mat(k,356) = rxt(k,327)*y(k,126)
         mat(k,2225) = mat(k,2225) + rxt(k,219)*y(k,124)
         mat(k,1982) = mat(k,1982) + rxt(k,187)*y(k,124)
         mat(k,407) = rxt(k,150)*y(k,217)
         mat(k,589) = .700_r8*rxt(k,376)*y(k,217)
         mat(k,1205) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126)
         mat(k,1849) = rxt(k,219)*y(k,19) + rxt(k,187)*y(k,59) + rxt(k,382)*y(k,101) &
                      + 2.000_r8*rxt(k,160)*y(k,126) + rxt(k,166)*y(k,133) &
                      + rxt(k,165)*y(k,134) + rxt(k,397)*y(k,188) + rxt(k,358) &
                      *y(k,189) + rxt(k,400)*y(k,191) + rxt(k,405)*y(k,193) &
                      + rxt(k,283)*y(k,194) + rxt(k,311)*y(k,195) + rxt(k,407) &
                      *y(k,196) + rxt(k,294)*y(k,197) + rxt(k,262)*y(k,198) &
                      + rxt(k,411)*y(k,199) + rxt(k,329)*y(k,200) + rxt(k,298) &
                      *y(k,202) + rxt(k,164)*y(k,203) + rxt(k,270)*y(k,204) &
                      + .920_r8*rxt(k,368)*y(k,205) + .920_r8*rxt(k,374)*y(k,206) &
                      + rxt(k,336)*y(k,209) + rxt(k,414)*y(k,210) + rxt(k,345) &
                      *y(k,211) + rxt(k,417)*y(k,212) + rxt(k,348)*y(k,213) &
                      + 1.600_r8*rxt(k,447)*y(k,215) + rxt(k,420)*y(k,218) &
                      + rxt(k,319)*y(k,219) + rxt(k,323)*y(k,220) + .900_r8*rxt(k,452) &
                      *y(k,221) + .800_r8*rxt(k,457)*y(k,222) + rxt(k,427)*y(k,223) &
                      + rxt(k,393)*y(k,225) + rxt(k,433)*y(k,226) + rxt(k,436) &
                      *y(k,228)
         mat(k,1757) = mat(k,1757) + rxt(k,327)*y(k,16) + rxt(k,383)*y(k,101) &
                      + 2.000_r8*rxt(k,160)*y(k,124) + rxt(k,161)*y(k,133) &
                      + rxt(k,159)*y(k,203) + rxt(k,369)*y(k,205) + rxt(k,375) &
                      *y(k,206) + rxt(k,335)*y(k,209) + rxt(k,346)*y(k,211) &
                      + 2.000_r8*rxt(k,448)*y(k,215) + rxt(k,162)*y(k,217) &
                      + rxt(k,394)*y(k,225)
         mat(k,862) = rxt(k,317)*y(k,217)
         mat(k,2256) = mat(k,2256) + rxt(k,166)*y(k,124) + rxt(k,161)*y(k,126)
         mat(k,2134) = mat(k,2134) + rxt(k,165)*y(k,124)
         mat(k,622) = rxt(k,454)*y(k,217)
         mat(k,523) = mat(k,523) + rxt(k,397)*y(k,124)
         mat(k,908) = rxt(k,358)*y(k,124)
         mat(k,489) = rxt(k,400)*y(k,124)
         mat(k,440) = rxt(k,405)*y(k,124)
         mat(k,822) = rxt(k,283)*y(k,124)
         mat(k,799) = rxt(k,311)*y(k,124)
         mat(k,630) = rxt(k,407)*y(k,124)
         mat(k,1396) = mat(k,1396) + rxt(k,294)*y(k,124)
         mat(k,2073) = rxt(k,262)*y(k,124) + .500_r8*rxt(k,445)*y(k,215)
         mat(k,766) = mat(k,766) + rxt(k,411)*y(k,124)
         mat(k,571) = rxt(k,329)*y(k,124)
         mat(k,720) = rxt(k,298)*y(k,124)
         mat(k,1956) = mat(k,1956) + rxt(k,164)*y(k,124) + rxt(k,159)*y(k,126)
         mat(k,446) = rxt(k,270)*y(k,124)
         mat(k,1269) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126)
         mat(k,1301) = .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,1321) = rxt(k,336)*y(k,124) + rxt(k,335)*y(k,126)
         mat(k,726) = mat(k,726) + rxt(k,414)*y(k,124)
         mat(k,1364) = mat(k,1364) + rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126)
         mat(k,881) = mat(k,881) + rxt(k,417)*y(k,124)
         mat(k,683) = rxt(k,348)*y(k,124)
         mat(k,1056) = 1.600_r8*rxt(k,447)*y(k,124) + 2.000_r8*rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,198)
         mat(k,1700) = mat(k,1700) + rxt(k,356)*y(k,1) + rxt(k,150)*y(k,90) &
                      + .700_r8*rxt(k,376)*y(k,99) + rxt(k,162)*y(k,126) + rxt(k,317) &
                      *y(k,127) + rxt(k,454)*y(k,175)
         mat(k,459) = rxt(k,420)*y(k,124)
         mat(k,775) = rxt(k,319)*y(k,124)
         mat(k,1138) = rxt(k,323)*y(k,124)
         mat(k,1124) = .900_r8*rxt(k,452)*y(k,124)
         mat(k,921) = .800_r8*rxt(k,457)*y(k,124)
         mat(k,741) = rxt(k,427)*y(k,124)
         mat(k,1182) = rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126)
         mat(k,758) = rxt(k,433)*y(k,124)
         mat(k,505) = rxt(k,436)*y(k,124)
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
         mat(k,1749) = -(rxt(k,159)*y(k,203) + rxt(k,160)*y(k,124) + rxt(k,161) &
                      *y(k,133) + rxt(k,162)*y(k,217) + rxt(k,170)*y(k,125) + rxt(k,256) &
                      *y(k,42) + rxt(k,288)*y(k,45) + rxt(k,307)*y(k,29) + rxt(k,314) &
                      *y(k,49) + rxt(k,327)*y(k,16) + rxt(k,335)*y(k,209) + rxt(k,346) &
                      *y(k,211) + rxt(k,369)*y(k,205) + rxt(k,375)*y(k,206) + rxt(k,378) &
                      *y(k,98) + rxt(k,383)*y(k,101) + rxt(k,394)*y(k,225) + rxt(k,439) &
                      *y(k,6) + rxt(k,442)*y(k,110) + rxt(k,448)*y(k,215) + rxt(k,459) &
                      *y(k,177) + rxt(k,462)*y(k,67))
         mat(k,1948) = -rxt(k,159)*y(k,126)
         mat(k,1841) = -rxt(k,160)*y(k,126)
         mat(k,2248) = -rxt(k,161)*y(k,126)
         mat(k,1692) = -rxt(k,162)*y(k,126)
         mat(k,2193) = -rxt(k,170)*y(k,126)
         mat(k,1488) = -rxt(k,256)*y(k,126)
         mat(k,1090) = -rxt(k,288)*y(k,126)
         mat(k,1033) = -rxt(k,307)*y(k,126)
         mat(k,1225) = -rxt(k,314)*y(k,126)
         mat(k,355) = -rxt(k,327)*y(k,126)
         mat(k,1317) = -rxt(k,335)*y(k,126)
         mat(k,1359) = -rxt(k,346)*y(k,126)
         mat(k,1264) = -rxt(k,369)*y(k,126)
         mat(k,1296) = -rxt(k,375)*y(k,126)
         mat(k,853) = -rxt(k,378)*y(k,126)
         mat(k,1201) = -rxt(k,383)*y(k,126)
         mat(k,1178) = -rxt(k,394)*y(k,126)
         mat(k,1011) = -rxt(k,439)*y(k,126)
         mat(k,961) = -rxt(k,442)*y(k,126)
         mat(k,1052) = -rxt(k,448)*y(k,126)
         mat(k,975) = -rxt(k,459)*y(k,126)
         mat(k,301) = -rxt(k,462)*y(k,126)
         mat(k,544) = rxt(k,221)*y(k,133)
         mat(k,2013) = rxt(k,188)*y(k,60)
         mat(k,887) = rxt(k,188)*y(k,56) + rxt(k,190)*y(k,133) + rxt(k,191)*y(k,217)
         mat(k,870) = rxt(k,235)*y(k,89)
         mat(k,1452) = rxt(k,235)*y(k,73) + rxt(k,172)*y(k,217)
         mat(k,578) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,2193) = mat(k,2193) + rxt(k,158)*y(k,133) + rxt(k,157)*y(k,134)
         mat(k,2248) = mat(k,2248) + rxt(k,221)*y(k,20) + rxt(k,190)*y(k,60) &
                      + rxt(k,158)*y(k,125)
         mat(k,2126) = rxt(k,157)*y(k,125)
         mat(k,529) = rxt(k,303)*y(k,217)
         mat(k,1692) = mat(k,1692) + rxt(k,191)*y(k,60) + rxt(k,172)*y(k,89) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,139)
         mat(k,857) = -(rxt(k,317)*y(k,217))
         mat(k,1650) = -rxt(k,317)*y(k,127)
         mat(k,1022) = rxt(k,307)*y(k,126)
         mat(k,557) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,399) = rxt(k,384)*y(k,217)
         mat(k,374) = rxt(k,388)*y(k,217)
         mat(k,1070) = rxt(k,389)*y(k,217)
         mat(k,1713) = rxt(k,307)*y(k,29)
         mat(k,1650) = mat(k,1650) + .500_r8*rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) &
                      + rxt(k,388)*y(k,115) + rxt(k,389)*y(k,116)
         mat(k,391) = -(rxt(k,449)*y(k,217))
         mat(k,1595) = -rxt(k,449)*y(k,128)
         mat(k,1875) = rxt(k,446)*y(k,215)
         mat(k,1041) = rxt(k,446)*y(k,203)
         mat(k,2258) = -(rxt(k,130)*y(k,134) + 4._r8*rxt(k,131)*y(k,133) + rxt(k,133) &
                      *y(k,77) + rxt(k,134)*y(k,79) + rxt(k,139)*y(k,203) + rxt(k,145) &
                      *y(k,217) + (rxt(k,156) + rxt(k,158)) * y(k,125) + rxt(k,161) &
                      *y(k,126) + rxt(k,166)*y(k,124) + rxt(k,190)*y(k,60) + rxt(k,192) &
                      *y(k,59) + rxt(k,195)*y(k,85) + rxt(k,198)*y(k,92) + rxt(k,221) &
                      *y(k,20) + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,81) + rxt(k,226) &
                      *y(k,91) + rxt(k,257)*y(k,42) + rxt(k,464)*y(k,137))
         mat(k,2136) = -rxt(k,130)*y(k,133)
         mat(k,1411) = -rxt(k,133)*y(k,133)
         mat(k,608) = -rxt(k,134)*y(k,133)
         mat(k,1958) = -rxt(k,139)*y(k,133)
         mat(k,1702) = -rxt(k,145)*y(k,133)
         mat(k,2203) = -(rxt(k,156) + rxt(k,158)) * y(k,133)
         mat(k,1759) = -rxt(k,161)*y(k,133)
         mat(k,1851) = -rxt(k,166)*y(k,133)
         mat(k,892) = -rxt(k,190)*y(k,133)
         mat(k,1984) = -rxt(k,192)*y(k,133)
         mat(k,2159) = -rxt(k,195)*y(k,133)
         mat(k,829) = -rxt(k,198)*y(k,133)
         mat(k,547) = -rxt(k,221)*y(k,133)
         mat(k,2227) = -rxt(k,222)*y(k,133)
         mat(k,810) = -rxt(k,224)*y(k,133)
         mat(k,784) = -rxt(k,226)*y(k,133)
         mat(k,1497) = -rxt(k,257)*y(k,133)
         mat(k,364) = -rxt(k,464)*y(k,133)
         mat(k,1474) = rxt(k,137)*y(k,203)
         mat(k,497) = rxt(k,151)*y(k,124) + rxt(k,152)*y(k,125)
         mat(k,1851) = mat(k,1851) + rxt(k,151)*y(k,112)
         mat(k,2203) = mat(k,2203) + rxt(k,152)*y(k,112)
         mat(k,1958) = mat(k,1958) + rxt(k,137)*y(k,76)
         mat(k,1702) = mat(k,1702) + 2.000_r8*rxt(k,147)*y(k,217)
         mat(k,2132) = -(rxt(k,129)*y(k,216) + rxt(k,130)*y(k,133) + rxt(k,140) &
                      *y(k,203) + rxt(k,141)*y(k,76) + rxt(k,146)*y(k,217) + rxt(k,157) &
                      *y(k,125) + rxt(k,165)*y(k,124) + rxt(k,181)*y(k,56) + rxt(k,213) &
                      *y(k,17) + rxt(k,279)*y(k,25) + rxt(k,308)*y(k,29) + rxt(k,338) &
                      *y(k,105) + rxt(k,352)*y(k,111) + rxt(k,385)*y(k,98) + rxt(k,423) &
                      *y(k,141) + rxt(k,440)*y(k,6) + rxt(k,443)*y(k,110) + rxt(k,468) &
                      *y(k,148) + rxt(k,474)*y(k,150))
         mat(k,1534) = -rxt(k,129)*y(k,134)
         mat(k,2254) = -rxt(k,130)*y(k,134)
         mat(k,1954) = -rxt(k,140)*y(k,134)
         mat(k,1471) = -rxt(k,141)*y(k,134)
         mat(k,1698) = -rxt(k,146)*y(k,134)
         mat(k,2199) = -rxt(k,157)*y(k,134)
         mat(k,1847) = -rxt(k,165)*y(k,134)
         mat(k,2019) = -rxt(k,181)*y(k,134)
         mat(k,1421) = -rxt(k,213)*y(k,134)
         mat(k,555) = -rxt(k,279)*y(k,134)
         mat(k,1037) = -rxt(k,308)*y(k,134)
         mat(k,1217) = -rxt(k,338)*y(k,134)
         mat(k,1343) = -rxt(k,352)*y(k,134)
         mat(k,856) = -rxt(k,385)*y(k,134)
         mat(k,465) = -rxt(k,423)*y(k,134)
         mat(k,1015) = -rxt(k,440)*y(k,134)
         mat(k,965) = -rxt(k,443)*y(k,134)
         mat(k,515) = -rxt(k,468)*y(k,134)
         mat(k,1240) = -rxt(k,474)*y(k,134)
         mat(k,1395) = .150_r8*rxt(k,293)*y(k,203)
         mat(k,1954) = mat(k,1954) + .150_r8*rxt(k,293)*y(k,197) + .150_r8*rxt(k,343) &
                      *y(k,211)
         mat(k,1363) = .150_r8*rxt(k,343)*y(k,203)
         mat(k,328) = -(rxt(k,475)*y(k,150))
         mat(k,1229) = -rxt(k,475)*y(k,136)
         mat(k,2206) = rxt(k,215)*y(k,59)
         mat(k,1963) = rxt(k,215)*y(k,19) + 2.000_r8*rxt(k,185)*y(k,59)
         mat(k,357) = -(rxt(k,464)*y(k,133) + rxt(k,465)*y(k,217))
         mat(k,2229) = -rxt(k,464)*y(k,137)
         mat(k,1591) = -rxt(k,465)*y(k,137)
         mat(k,1146) = rxt(k,331)*y(k,217)
         mat(k,1773) = .100_r8*rxt(k,452)*y(k,221)
         mat(k,1574) = rxt(k,331)*y(k,93)
         mat(k,1107) = .100_r8*rxt(k,452)*y(k,124)
         mat(k,524) = -(rxt(k,303)*y(k,217))
         mat(k,1615) = -rxt(k,303)*y(k,139)
         mat(k,2168) = rxt(k,305)*y(k,197)
         mat(k,1368) = rxt(k,305)*y(k,125)
         mat(k,2161) = rxt(k,425)*y(k,188)
         mat(k,517) = rxt(k,425)*y(k,125)
         mat(k,463) = -(rxt(k,422)*y(k,125) + rxt(k,423)*y(k,134))
         mat(k,2165) = -rxt(k,422)*y(k,141)
         mat(k,2084) = -rxt(k,423)*y(k,141)
         mat(k,196) = .070_r8*rxt(k,409)*y(k,217)
         mat(k,1783) = rxt(k,407)*y(k,196)
         mat(k,176) = .060_r8*rxt(k,421)*y(k,217)
         mat(k,217) = .070_r8*rxt(k,437)*y(k,217)
         mat(k,624) = rxt(k,407)*y(k,124)
         mat(k,1606) = .070_r8*rxt(k,409)*y(k,66) + .060_r8*rxt(k,421)*y(k,142) &
                      + .070_r8*rxt(k,437)*y(k,184)
         mat(k,174) = -(rxt(k,421)*y(k,217))
         mat(k,1561) = -rxt(k,421)*y(k,142)
         mat(k,166) = .530_r8*rxt(k,398)*y(k,217)
         mat(k,1561) = mat(k,1561) + .530_r8*rxt(k,398)*y(k,7)
         mat(k,333) = -(rxt(k,424)*y(k,217))
         mat(k,1587) = -rxt(k,424)*y(k,143)
         mat(k,1871) = rxt(k,419)*y(k,218)
         mat(k,453) = rxt(k,419)*y(k,203)
         mat(k,532) = -(rxt(k,320)*y(k,217))
         mat(k,1616) = -rxt(k,320)*y(k,146)
         mat(k,1891) = rxt(k,318)*y(k,219)
         mat(k,767) = rxt(k,318)*y(k,203)
         mat(k,409) = -(rxt(k,324)*y(k,217))
         mat(k,1598) = -rxt(k,324)*y(k,147)
         mat(k,1878) = .850_r8*rxt(k,322)*y(k,220)
         mat(k,1127) = .850_r8*rxt(k,322)*y(k,203)
         mat(k,511) = -(rxt(k,468)*y(k,134) + rxt(k,471)*y(k,217))
         mat(k,2085) = -rxt(k,468)*y(k,148)
         mat(k,1613) = -rxt(k,471)*y(k,148)
         mat(k,1232) = -(rxt(k,469)*y(k,19) + rxt(k,470)*y(k,59) + rxt(k,472)*y(k,125) &
                      + rxt(k,474)*y(k,134) + rxt(k,475)*y(k,136) + rxt(k,476) &
                      *y(k,217))
         mat(k,2210) = -rxt(k,469)*y(k,150)
         mat(k,1967) = -rxt(k,470)*y(k,150)
         mat(k,2183) = -rxt(k,472)*y(k,150)
         mat(k,2112) = -rxt(k,474)*y(k,150)
         mat(k,330) = -rxt(k,475)*y(k,150)
         mat(k,1677) = -rxt(k,476)*y(k,150)
         mat(k,2239) = rxt(k,464)*y(k,137)
         mat(k,2112) = mat(k,2112) + rxt(k,468)*y(k,148)
         mat(k,361) = rxt(k,464)*y(k,133)
         mat(k,512) = rxt(k,468)*y(k,134) + rxt(k,471)*y(k,217)
         mat(k,1677) = mat(k,1677) + rxt(k,471)*y(k,148)
         mat(k,832) = -(rxt(k,467)*y(k,217))
         mat(k,1648) = -rxt(k,467)*y(k,151)
         mat(k,2209) = rxt(k,469)*y(k,150)
         mat(k,1965) = rxt(k,470)*y(k,150)
         mat(k,298) = rxt(k,462)*y(k,126) + (rxt(k,463)+.500_r8*rxt(k,477))*y(k,217)
         mat(k,2176) = rxt(k,472)*y(k,150)
         mat(k,1711) = rxt(k,462)*y(k,67)
         mat(k,2091) = rxt(k,474)*y(k,150)
         mat(k,329) = rxt(k,475)*y(k,150)
         mat(k,359) = rxt(k,465)*y(k,217)
         mat(k,1231) = rxt(k,469)*y(k,19) + rxt(k,470)*y(k,59) + rxt(k,472)*y(k,125) &
                      + rxt(k,474)*y(k,134) + rxt(k,475)*y(k,136) + rxt(k,476) &
                      *y(k,217)
         mat(k,1648) = mat(k,1648) + (rxt(k,463)+.500_r8*rxt(k,477))*y(k,67) &
                      + rxt(k,465)*y(k,137) + rxt(k,476)*y(k,150)
         mat(k,257) = -(rxt(k,478)*y(k,229))
         mat(k,2261) = -rxt(k,478)*y(k,152)
         mat(k,831) = rxt(k,467)*y(k,217)
         mat(k,1577) = rxt(k,467)*y(k,151)
         mat(k,984) = .2202005_r8*rxt(k,497)*y(k,134)
         mat(k,935) = .0508005_r8*rxt(k,513)*y(k,134)
         mat(k,1761) = .1279005_r8*rxt(k,496)*y(k,190) + .0097005_r8*rxt(k,501) &
                      *y(k,192) + .0003005_r8*rxt(k,504)*y(k,207) &
                      + .1056005_r8*rxt(k,508)*y(k,208) + .0245005_r8*rxt(k,512) &
                      *y(k,214) + .0154005_r8*rxt(k,518)*y(k,224) &
                      + .0063005_r8*rxt(k,521)*y(k,227)
         mat(k,2077) = .2202005_r8*rxt(k,497)*y(k,6) + .0508005_r8*rxt(k,513)*y(k,110)
         mat(k,43) = .5931005_r8*rxt(k,515)*y(k,217)
         mat(k,49) = .1279005_r8*rxt(k,496)*y(k,124) + .2202005_r8*rxt(k,495)*y(k,203)
         mat(k,55) = .0097005_r8*rxt(k,501)*y(k,124) + .0023005_r8*rxt(k,500)*y(k,203)
         mat(k,1853) = .2202005_r8*rxt(k,495)*y(k,190) + .0023005_r8*rxt(k,500) &
                      *y(k,192) + .0031005_r8*rxt(k,503)*y(k,207) &
                      + .2381005_r8*rxt(k,507)*y(k,208) + .0508005_r8*rxt(k,511) &
                      *y(k,214) + .1364005_r8*rxt(k,517)*y(k,224) &
                      + .1677005_r8*rxt(k,520)*y(k,227)
         mat(k,61) = .0003005_r8*rxt(k,504)*y(k,124) + .0031005_r8*rxt(k,503)*y(k,203)
         mat(k,67) = .1056005_r8*rxt(k,508)*y(k,124) + .2381005_r8*rxt(k,507)*y(k,203)
         mat(k,75) = .0245005_r8*rxt(k,512)*y(k,124) + .0508005_r8*rxt(k,511)*y(k,203)
         mat(k,1540) = .5931005_r8*rxt(k,515)*y(k,172)
         mat(k,81) = .0154005_r8*rxt(k,518)*y(k,124) + .1364005_r8*rxt(k,517)*y(k,203)
         mat(k,87) = .0063005_r8*rxt(k,521)*y(k,124) + .1677005_r8*rxt(k,520)*y(k,203)
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
         mat(k,985) = .2067005_r8*rxt(k,497)*y(k,134)
         mat(k,936) = .1149005_r8*rxt(k,513)*y(k,134)
         mat(k,1762) = .1792005_r8*rxt(k,496)*y(k,190) + .0034005_r8*rxt(k,501) &
                      *y(k,192) + .0003005_r8*rxt(k,504)*y(k,207) &
                      + .1026005_r8*rxt(k,508)*y(k,208) + .0082005_r8*rxt(k,512) &
                      *y(k,214) + .0452005_r8*rxt(k,518)*y(k,224) &
                      + .0237005_r8*rxt(k,521)*y(k,227)
         mat(k,2078) = .2067005_r8*rxt(k,497)*y(k,6) + .1149005_r8*rxt(k,513)*y(k,110)
         mat(k,44) = .1534005_r8*rxt(k,515)*y(k,217)
         mat(k,50) = .1792005_r8*rxt(k,496)*y(k,124) + .2067005_r8*rxt(k,495)*y(k,203)
         mat(k,56) = .0034005_r8*rxt(k,501)*y(k,124) + .0008005_r8*rxt(k,500)*y(k,203)
         mat(k,1854) = .2067005_r8*rxt(k,495)*y(k,190) + .0008005_r8*rxt(k,500) &
                      *y(k,192) + .0035005_r8*rxt(k,503)*y(k,207) &
                      + .1308005_r8*rxt(k,507)*y(k,208) + .1149005_r8*rxt(k,511) &
                      *y(k,214) + .0101005_r8*rxt(k,517)*y(k,224) &
                      + .0174005_r8*rxt(k,520)*y(k,227)
         mat(k,62) = .0003005_r8*rxt(k,504)*y(k,124) + .0035005_r8*rxt(k,503)*y(k,203)
         mat(k,68) = .1026005_r8*rxt(k,508)*y(k,124) + .1308005_r8*rxt(k,507)*y(k,203)
         mat(k,76) = .0082005_r8*rxt(k,512)*y(k,124) + .1149005_r8*rxt(k,511)*y(k,203)
         mat(k,1541) = .1534005_r8*rxt(k,515)*y(k,172)
         mat(k,82) = .0452005_r8*rxt(k,518)*y(k,124) + .0101005_r8*rxt(k,517)*y(k,203)
         mat(k,88) = .0237005_r8*rxt(k,521)*y(k,124) + .0174005_r8*rxt(k,520)*y(k,203)
         mat(k,986) = .0653005_r8*rxt(k,497)*y(k,134)
         mat(k,937) = .0348005_r8*rxt(k,513)*y(k,134)
         mat(k,1763) = .0676005_r8*rxt(k,496)*y(k,190) + .1579005_r8*rxt(k,501) &
                      *y(k,192) + .0073005_r8*rxt(k,504)*y(k,207) &
                      + .0521005_r8*rxt(k,508)*y(k,208) + .0772005_r8*rxt(k,512) &
                      *y(k,214) + .0966005_r8*rxt(k,518)*y(k,224) &
                      + .0025005_r8*rxt(k,521)*y(k,227)
         mat(k,2079) = .0653005_r8*rxt(k,497)*y(k,6) + .0348005_r8*rxt(k,513)*y(k,110)
         mat(k,45) = .0459005_r8*rxt(k,515)*y(k,217)
         mat(k,51) = .0676005_r8*rxt(k,496)*y(k,124) + .0653005_r8*rxt(k,495)*y(k,203)
         mat(k,57) = .1579005_r8*rxt(k,501)*y(k,124) + .0843005_r8*rxt(k,500)*y(k,203)
         mat(k,1855) = .0653005_r8*rxt(k,495)*y(k,190) + .0843005_r8*rxt(k,500) &
                      *y(k,192) + .0003005_r8*rxt(k,503)*y(k,207) &
                      + .0348005_r8*rxt(k,507)*y(k,208) + .0348005_r8*rxt(k,511) &
                      *y(k,214) + .0763005_r8*rxt(k,517)*y(k,224) + .086_r8*rxt(k,520) &
                      *y(k,227)
         mat(k,63) = .0073005_r8*rxt(k,504)*y(k,124) + .0003005_r8*rxt(k,503)*y(k,203)
         mat(k,69) = .0521005_r8*rxt(k,508)*y(k,124) + .0348005_r8*rxt(k,507)*y(k,203)
         mat(k,77) = .0772005_r8*rxt(k,512)*y(k,124) + .0348005_r8*rxt(k,511)*y(k,203)
         mat(k,1542) = .0459005_r8*rxt(k,515)*y(k,172)
         mat(k,83) = .0966005_r8*rxt(k,518)*y(k,124) + .0763005_r8*rxt(k,517)*y(k,203)
         mat(k,89) = .0025005_r8*rxt(k,521)*y(k,124) + .086_r8*rxt(k,520)*y(k,203)
         mat(k,987) = .1749305_r8*rxt(k,494)*y(k,126) + .1284005_r8*rxt(k,497) &
                      *y(k,134)
         mat(k,838) = .0590245_r8*rxt(k,502)*y(k,126) + .0033005_r8*rxt(k,505) &
                      *y(k,134)
         mat(k,938) = .1749305_r8*rxt(k,510)*y(k,126) + .0554005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1764) = .079_r8*rxt(k,496)*y(k,190) + .0059005_r8*rxt(k,501)*y(k,192) &
                      + .0057005_r8*rxt(k,504)*y(k,207) + .0143005_r8*rxt(k,508) &
                      *y(k,208) + .0332005_r8*rxt(k,512)*y(k,214) &
                      + .0073005_r8*rxt(k,518)*y(k,224) + .011_r8*rxt(k,521)*y(k,227)
         mat(k,1704) = .1749305_r8*rxt(k,494)*y(k,6) + .0590245_r8*rxt(k,502)*y(k,98) &
                      + .1749305_r8*rxt(k,510)*y(k,110)
         mat(k,2080) = .1284005_r8*rxt(k,497)*y(k,6) + .0033005_r8*rxt(k,505)*y(k,98) &
                      + .0554005_r8*rxt(k,513)*y(k,110)
         mat(k,46) = .0085005_r8*rxt(k,515)*y(k,217)
         mat(k,52) = .079_r8*rxt(k,496)*y(k,124) + .1284005_r8*rxt(k,495)*y(k,203)
         mat(k,58) = .0059005_r8*rxt(k,501)*y(k,124) + .0443005_r8*rxt(k,500)*y(k,203)
         mat(k,1856) = .1284005_r8*rxt(k,495)*y(k,190) + .0443005_r8*rxt(k,500) &
                      *y(k,192) + .0271005_r8*rxt(k,503)*y(k,207) &
                      + .0076005_r8*rxt(k,507)*y(k,208) + .0554005_r8*rxt(k,511) &
                      *y(k,214) + .2157005_r8*rxt(k,517)*y(k,224) &
                      + .0512005_r8*rxt(k,520)*y(k,227)
         mat(k,64) = .0057005_r8*rxt(k,504)*y(k,124) + .0271005_r8*rxt(k,503)*y(k,203)
         mat(k,70) = .0143005_r8*rxt(k,508)*y(k,124) + .0076005_r8*rxt(k,507)*y(k,203)
         mat(k,78) = .0332005_r8*rxt(k,512)*y(k,124) + .0554005_r8*rxt(k,511)*y(k,203)
         mat(k,1543) = .0085005_r8*rxt(k,515)*y(k,172)
         mat(k,84) = .0073005_r8*rxt(k,518)*y(k,124) + .2157005_r8*rxt(k,517)*y(k,203)
         mat(k,90) = .011_r8*rxt(k,521)*y(k,124) + .0512005_r8*rxt(k,520)*y(k,203)
         mat(k,988) = .5901905_r8*rxt(k,494)*y(k,126) + .114_r8*rxt(k,497)*y(k,134)
         mat(k,839) = .0250245_r8*rxt(k,502)*y(k,126)
         mat(k,939) = .5901905_r8*rxt(k,510)*y(k,126) + .1278005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1765) = .1254005_r8*rxt(k,496)*y(k,190) + .0536005_r8*rxt(k,501) &
                      *y(k,192) + .0623005_r8*rxt(k,504)*y(k,207) &
                      + .0166005_r8*rxt(k,508)*y(k,208) + .130_r8*rxt(k,512)*y(k,214) &
                      + .238_r8*rxt(k,518)*y(k,224) + .1185005_r8*rxt(k,521)*y(k,227)
         mat(k,1705) = .5901905_r8*rxt(k,494)*y(k,6) + .0250245_r8*rxt(k,502)*y(k,98) &
                      + .5901905_r8*rxt(k,510)*y(k,110)
         mat(k,2081) = .114_r8*rxt(k,497)*y(k,6) + .1278005_r8*rxt(k,513)*y(k,110)
         mat(k,47) = .0128005_r8*rxt(k,515)*y(k,217)
         mat(k,53) = .1254005_r8*rxt(k,496)*y(k,124) + .114_r8*rxt(k,495)*y(k,203)
         mat(k,59) = .0536005_r8*rxt(k,501)*y(k,124) + .1621005_r8*rxt(k,500)*y(k,203)
         mat(k,1857) = .114_r8*rxt(k,495)*y(k,190) + .1621005_r8*rxt(k,500)*y(k,192) &
                      + .0474005_r8*rxt(k,503)*y(k,207) + .0113005_r8*rxt(k,507) &
                      *y(k,208) + .1278005_r8*rxt(k,511)*y(k,214) &
                      + .0738005_r8*rxt(k,517)*y(k,224) + .1598005_r8*rxt(k,520) &
                      *y(k,227)
         mat(k,65) = .0623005_r8*rxt(k,504)*y(k,124) + .0474005_r8*rxt(k,503)*y(k,203)
         mat(k,71) = .0166005_r8*rxt(k,508)*y(k,124) + .0113005_r8*rxt(k,507)*y(k,203)
         mat(k,79) = .130_r8*rxt(k,512)*y(k,124) + .1278005_r8*rxt(k,511)*y(k,203)
         mat(k,1544) = .0128005_r8*rxt(k,515)*y(k,172)
         mat(k,85) = .238_r8*rxt(k,518)*y(k,124) + .0738005_r8*rxt(k,517)*y(k,203)
         mat(k,91) = .1185005_r8*rxt(k,521)*y(k,124) + .1598005_r8*rxt(k,520)*y(k,203)
         mat(k,48) = -(rxt(k,515)*y(k,217))
         mat(k,1545) = -rxt(k,515)*y(k,172)
         mat(k,189) = .100_r8*rxt(k,429)*y(k,217)
         mat(k,207) = .230_r8*rxt(k,431)*y(k,217)
         mat(k,1565) = .100_r8*rxt(k,429)*y(k,180) + .230_r8*rxt(k,431)*y(k,182)
         mat(k,642) = -(rxt(k,453)*y(k,217))
         mat(k,1629) = -rxt(k,453)*y(k,174)
         mat(k,1896) = rxt(k,451)*y(k,221)
         mat(k,1108) = rxt(k,451)*y(k,203)
         mat(k,617) = -(rxt(k,454)*y(k,217))
         mat(k,1626) = -rxt(k,454)*y(k,175)
         mat(k,1793) = .200_r8*rxt(k,447)*y(k,215) + .200_r8*rxt(k,457)*y(k,222)
         mat(k,2029) = .500_r8*rxt(k,445)*y(k,215)
         mat(k,1042) = .200_r8*rxt(k,447)*y(k,124) + .500_r8*rxt(k,445)*y(k,198)
         mat(k,910) = .200_r8*rxt(k,457)*y(k,124)
         mat(k,474) = -(rxt(k,458)*y(k,217))
         mat(k,1608) = -rxt(k,458)*y(k,176)
         mat(k,1887) = rxt(k,456)*y(k,222)
         mat(k,909) = rxt(k,456)*y(k,203)
         mat(k,969) = -(rxt(k,459)*y(k,126) + rxt(k,460)*y(k,217))
         mat(k,1717) = -rxt(k,459)*y(k,177)
         mat(k,1659) = -rxt(k,460)*y(k,177)
         mat(k,997) = .330_r8*rxt(k,440)*y(k,134)
         mat(k,949) = .330_r8*rxt(k,443)*y(k,134)
         mat(k,1812) = .800_r8*rxt(k,447)*y(k,215) + .800_r8*rxt(k,457)*y(k,222)
         mat(k,1717) = mat(k,1717) + rxt(k,448)*y(k,215)
         mat(k,2097) = .330_r8*rxt(k,440)*y(k,6) + .330_r8*rxt(k,443)*y(k,110)
         mat(k,618) = rxt(k,454)*y(k,217)
         mat(k,2038) = .500_r8*rxt(k,445)*y(k,215) + rxt(k,455)*y(k,222)
         mat(k,1044) = .800_r8*rxt(k,447)*y(k,124) + rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,198)
         mat(k,1659) = mat(k,1659) + rxt(k,454)*y(k,175)
         mat(k,914) = .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,198)
         mat(k,1059) = -(rxt(k,461)*y(k,217))
         mat(k,1664) = -rxt(k,461)*y(k,178)
         mat(k,1001) = .300_r8*rxt(k,440)*y(k,134)
         mat(k,952) = .300_r8*rxt(k,443)*y(k,134)
         mat(k,1815) = .900_r8*rxt(k,452)*y(k,221)
         mat(k,2102) = .300_r8*rxt(k,440)*y(k,6) + .300_r8*rxt(k,443)*y(k,110)
         mat(k,2041) = rxt(k,450)*y(k,221)
         mat(k,1112) = .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,198)
         mat(k,655) = -(rxt(k,428)*y(k,217))
         mat(k,1630) = -rxt(k,428)*y(k,179)
         mat(k,1897) = rxt(k,426)*y(k,223)
         mat(k,730) = rxt(k,426)*y(k,203)
         mat(k,187) = -(rxt(k,429)*y(k,217))
         mat(k,1563) = -rxt(k,429)*y(k,180)
         mat(k,203) = -(rxt(k,395)*y(k,217))
         mat(k,1566) = -rxt(k,395)*y(k,181)
         mat(k,1866) = rxt(k,392)*y(k,225)
         mat(k,1166) = rxt(k,392)*y(k,203)
         mat(k,208) = -(rxt(k,431)*y(k,217))
         mat(k,1567) = -rxt(k,431)*y(k,182)
         mat(k,701) = -(rxt(k,434)*y(k,217))
         mat(k,1635) = -rxt(k,434)*y(k,183)
         mat(k,1901) = rxt(k,432)*y(k,226)
         mat(k,746) = rxt(k,432)*y(k,203)
         mat(k,216) = -(rxt(k,437)*y(k,217))
         mat(k,1568) = -rxt(k,437)*y(k,184)
         mat(k,209) = .150_r8*rxt(k,431)*y(k,217)
         mat(k,1568) = mat(k,1568) + .150_r8*rxt(k,431)*y(k,182)
         mat(k,427) = -(rxt(k,438)*y(k,217))
         mat(k,1601) = -rxt(k,438)*y(k,185)
         mat(k,1881) = rxt(k,435)*y(k,228)
         mat(k,498) = rxt(k,435)*y(k,203)
         mat(k,518) = -(rxt(k,396)*y(k,203) + rxt(k,397)*y(k,124) + rxt(k,425) &
                      *y(k,125))
         mat(k,1890) = -rxt(k,396)*y(k,188)
         mat(k,1788) = -rxt(k,397)*y(k,188)
         mat(k,2167) = -rxt(k,425)*y(k,188)
         mat(k,254) = rxt(k,402)*y(k,217)
         mat(k,1614) = rxt(k,402)*y(k,22)
         mat(k,899) = -(rxt(k,357)*y(k,203) + (rxt(k,358) + rxt(k,359)) * y(k,124))
         mat(k,1916) = -rxt(k,357)*y(k,189)
         mat(k,1808) = -(rxt(k,358) + rxt(k,359)) * y(k,189)
         mat(k,635) = rxt(k,360)*y(k,217)
         mat(k,239) = rxt(k,361)*y(k,217)
         mat(k,1654) = rxt(k,360)*y(k,2) + rxt(k,361)*y(k,15)
         mat(k,54) = -(rxt(k,495)*y(k,203) + rxt(k,496)*y(k,124))
         mat(k,1858) = -rxt(k,495)*y(k,190)
         mat(k,1766) = -rxt(k,496)*y(k,190)
         mat(k,989) = rxt(k,498)*y(k,217)
         mat(k,1546) = rxt(k,498)*y(k,6)
         mat(k,483) = -(rxt(k,399)*y(k,203) + rxt(k,400)*y(k,124))
         mat(k,1888) = -rxt(k,399)*y(k,191)
         mat(k,1784) = -rxt(k,400)*y(k,191)
         mat(k,167) = .350_r8*rxt(k,398)*y(k,217)
         mat(k,423) = rxt(k,401)*y(k,217)
         mat(k,1609) = .350_r8*rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8)
         mat(k,60) = -(rxt(k,500)*y(k,203) + rxt(k,501)*y(k,124))
         mat(k,1859) = -rxt(k,500)*y(k,192)
         mat(k,1767) = -rxt(k,501)*y(k,192)
         mat(k,163) = rxt(k,499)*y(k,217)
         mat(k,1547) = rxt(k,499)*y(k,7)
         mat(k,435) = -(rxt(k,403)*y(k,203) + rxt(k,405)*y(k,124))
         mat(k,1882) = -rxt(k,403)*y(k,193)
         mat(k,1779) = -rxt(k,405)*y(k,193)
         mat(k,340) = rxt(k,404)*y(k,217)
         mat(k,190) = .070_r8*rxt(k,429)*y(k,217)
         mat(k,210) = .060_r8*rxt(k,431)*y(k,217)
         mat(k,1602) = rxt(k,404)*y(k,23) + .070_r8*rxt(k,429)*y(k,180) &
                      + .060_r8*rxt(k,431)*y(k,182)
         mat(k,815) = -(4._r8*rxt(k,280)*y(k,194) + rxt(k,281)*y(k,198) + rxt(k,282) &
                      *y(k,203) + rxt(k,283)*y(k,124))
         mat(k,2034) = -rxt(k,281)*y(k,194)
         mat(k,1912) = -rxt(k,282)*y(k,194)
         mat(k,1805) = -rxt(k,283)*y(k,194)
         mat(k,345) = .500_r8*rxt(k,285)*y(k,217)
         mat(k,292) = rxt(k,286)*y(k,56) + rxt(k,287)*y(k,217)
         mat(k,1996) = rxt(k,286)*y(k,28)
         mat(k,1646) = .500_r8*rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28)
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
         mat(k,791) = -(rxt(k,309)*y(k,198) + rxt(k,310)*y(k,203) + rxt(k,311) &
                      *y(k,124))
         mat(k,2032) = -rxt(k,309)*y(k,195)
         mat(k,1909) = -rxt(k,310)*y(k,195)
         mat(k,1803) = -rxt(k,311)*y(k,195)
         mat(k,416) = rxt(k,312)*y(k,217)
         mat(k,110) = rxt(k,313)*y(k,217)
         mat(k,1643) = rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31)
         mat(k,625) = -(rxt(k,406)*y(k,203) + rxt(k,407)*y(k,124))
         mat(k,1894) = -rxt(k,406)*y(k,196)
         mat(k,1794) = -rxt(k,407)*y(k,196)
         mat(k,267) = rxt(k,408)*y(k,217)
         mat(k,1794) = mat(k,1794) + rxt(k,397)*y(k,188)
         mat(k,2087) = rxt(k,423)*y(k,141)
         mat(k,464) = rxt(k,423)*y(k,134)
         mat(k,519) = rxt(k,397)*y(k,124) + .400_r8*rxt(k,396)*y(k,203)
         mat(k,1894) = mat(k,1894) + .400_r8*rxt(k,396)*y(k,188)
         mat(k,1627) = rxt(k,408)*y(k,32)
         mat(k,1386) = -(4._r8*rxt(k,291)*y(k,197) + rxt(k,292)*y(k,198) + rxt(k,293) &
                      *y(k,203) + rxt(k,294)*y(k,124) + rxt(k,305)*y(k,125) + rxt(k,332) &
                      *y(k,209) + rxt(k,365)*y(k,205) + rxt(k,370)*y(k,206) + rxt(k,379) &
                      *y(k,101) + rxt(k,390)*y(k,225))
         mat(k,2058) = -rxt(k,292)*y(k,197)
         mat(k,1939) = -rxt(k,293)*y(k,197)
         mat(k,1833) = -rxt(k,294)*y(k,197)
         mat(k,2185) = -rxt(k,305)*y(k,197)
         mat(k,1313) = -rxt(k,332)*y(k,197)
         mat(k,1260) = -rxt(k,365)*y(k,197)
         mat(k,1292) = -rxt(k,370)*y(k,197)
         mat(k,1197) = -rxt(k,379)*y(k,197)
         mat(k,1175) = -rxt(k,390)*y(k,197)
         mat(k,1007) = .060_r8*rxt(k,440)*y(k,134)
         mat(k,1087) = rxt(k,288)*y(k,126) + rxt(k,289)*y(k,217)
         mat(k,1222) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,217)
         mat(k,612) = .500_r8*rxt(k,296)*y(k,217)
         mat(k,850) = .080_r8*rxt(k,385)*y(k,134)
         mat(k,1213) = .100_r8*rxt(k,338)*y(k,134)
         mat(k,957) = .060_r8*rxt(k,443)*y(k,134)
         mat(k,1334) = .280_r8*rxt(k,352)*y(k,134)
         mat(k,1833) = mat(k,1833) + .530_r8*rxt(k,336)*y(k,209) + rxt(k,345)*y(k,211) &
                      + rxt(k,348)*y(k,213) + rxt(k,323)*y(k,220)
         mat(k,1741) = rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) + .530_r8*rxt(k,335) &
                      *y(k,209) + rxt(k,346)*y(k,211)
         mat(k,2118) = .060_r8*rxt(k,440)*y(k,6) + .080_r8*rxt(k,385)*y(k,98) &
                      + .100_r8*rxt(k,338)*y(k,105) + .060_r8*rxt(k,443)*y(k,110) &
                      + .280_r8*rxt(k,352)*y(k,111)
         mat(k,1062) = .650_r8*rxt(k,461)*y(k,217)
         mat(k,1386) = mat(k,1386) + .530_r8*rxt(k,332)*y(k,209)
         mat(k,2058) = mat(k,2058) + .260_r8*rxt(k,333)*y(k,209) + rxt(k,342)*y(k,211) &
                      + .300_r8*rxt(k,321)*y(k,220)
         mat(k,1939) = mat(k,1939) + .450_r8*rxt(k,343)*y(k,211) + .200_r8*rxt(k,347) &
                      *y(k,213) + .150_r8*rxt(k,322)*y(k,220)
         mat(k,1313) = mat(k,1313) + .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335) &
                      *y(k,126) + .530_r8*rxt(k,332)*y(k,197) + .260_r8*rxt(k,333) &
                      *y(k,198)
         mat(k,1355) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,342)*y(k,198) &
                      + .450_r8*rxt(k,343)*y(k,203) + 4.000_r8*rxt(k,344)*y(k,211)
         mat(k,679) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,203)
         mat(k,1683) = rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49) + .500_r8*rxt(k,296) &
                      *y(k,51) + .650_r8*rxt(k,461)*y(k,178)
         mat(k,1132) = rxt(k,323)*y(k,124) + .300_r8*rxt(k,321)*y(k,198) &
                      + .150_r8*rxt(k,322)*y(k,203)
         mat(k,2070) = -(rxt(k,182)*y(k,59) + (4._r8*rxt(k,259) + 4._r8*rxt(k,260) &
                      ) * y(k,198) + rxt(k,261)*y(k,203) + rxt(k,262)*y(k,124) &
                      + rxt(k,281)*y(k,194) + rxt(k,292)*y(k,197) + rxt(k,309) &
                      *y(k,195) + rxt(k,321)*y(k,220) + rxt(k,333)*y(k,209) + rxt(k,342) &
                      *y(k,211) + rxt(k,366)*y(k,205) + rxt(k,371)*y(k,206) + rxt(k,380) &
                      *y(k,101) + rxt(k,391)*y(k,225) + rxt(k,445)*y(k,215) + rxt(k,450) &
                      *y(k,221) + rxt(k,455)*y(k,222))
         mat(k,1979) = -rxt(k,182)*y(k,198)
         mat(k,1953) = -rxt(k,261)*y(k,198)
         mat(k,1846) = -rxt(k,262)*y(k,198)
         mat(k,821) = -rxt(k,281)*y(k,198)
         mat(k,1394) = -rxt(k,292)*y(k,198)
         mat(k,798) = -rxt(k,309)*y(k,198)
         mat(k,1137) = -rxt(k,321)*y(k,198)
         mat(k,1320) = -rxt(k,333)*y(k,198)
         mat(k,1362) = -rxt(k,342)*y(k,198)
         mat(k,1267) = -rxt(k,366)*y(k,198)
         mat(k,1299) = -rxt(k,371)*y(k,198)
         mat(k,1204) = -rxt(k,380)*y(k,198)
         mat(k,1181) = -rxt(k,391)*y(k,198)
         mat(k,1055) = -rxt(k,445)*y(k,198)
         mat(k,1123) = -rxt(k,450)*y(k,198)
         mat(k,920) = -rxt(k,455)*y(k,198)
         mat(k,1036) = .280_r8*rxt(k,308)*y(k,134)
         mat(k,687) = rxt(k,295)*y(k,217)
         mat(k,389) = .700_r8*rxt(k,264)*y(k,217)
         mat(k,1439) = rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73) + rxt(k,271)*y(k,216) &
                      + rxt(k,265)*y(k,217)
         mat(k,2018) = rxt(k,176)*y(k,54)
         mat(k,871) = rxt(k,232)*y(k,54)
         mat(k,855) = .050_r8*rxt(k,385)*y(k,134)
         mat(k,1204) = mat(k,1204) + rxt(k,379)*y(k,197)
         mat(k,1846) = mat(k,1846) + rxt(k,294)*y(k,197) + .830_r8*rxt(k,411)*y(k,199) &
                      + .170_r8*rxt(k,417)*y(k,212)
         mat(k,2131) = .280_r8*rxt(k,308)*y(k,29) + .050_r8*rxt(k,385)*y(k,98)
         mat(k,1394) = mat(k,1394) + rxt(k,379)*y(k,101) + rxt(k,294)*y(k,124) &
                      + 4.000_r8*rxt(k,291)*y(k,197) + .900_r8*rxt(k,292)*y(k,198) &
                      + .450_r8*rxt(k,293)*y(k,203) + rxt(k,365)*y(k,205) + rxt(k,370) &
                      *y(k,206) + rxt(k,332)*y(k,209) + rxt(k,341)*y(k,211) &
                      + rxt(k,390)*y(k,225)
         mat(k,2070) = mat(k,2070) + .900_r8*rxt(k,292)*y(k,197)
         mat(k,765) = .830_r8*rxt(k,411)*y(k,124) + .330_r8*rxt(k,410)*y(k,203)
         mat(k,1953) = mat(k,1953) + .450_r8*rxt(k,293)*y(k,197) + .330_r8*rxt(k,410) &
                      *y(k,199) + .070_r8*rxt(k,416)*y(k,212)
         mat(k,1267) = mat(k,1267) + rxt(k,365)*y(k,197)
         mat(k,1299) = mat(k,1299) + rxt(k,370)*y(k,197)
         mat(k,1320) = mat(k,1320) + rxt(k,332)*y(k,197)
         mat(k,1362) = mat(k,1362) + rxt(k,341)*y(k,197)
         mat(k,880) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,203)
         mat(k,1533) = rxt(k,271)*y(k,54)
         mat(k,1697) = rxt(k,295)*y(k,50) + .700_r8*rxt(k,264)*y(k,53) + rxt(k,265) &
                      *y(k,54)
         mat(k,1181) = mat(k,1181) + rxt(k,390)*y(k,197)
         mat(k,759) = -(rxt(k,410)*y(k,203) + rxt(k,411)*y(k,124) + rxt(k,412) &
                      *y(k,125))
         mat(k,1906) = -rxt(k,410)*y(k,199)
         mat(k,1801) = -rxt(k,411)*y(k,199)
         mat(k,2173) = -rxt(k,412)*y(k,199)
         mat(k,564) = -((rxt(k,329) + rxt(k,330)) * y(k,124))
         mat(k,1790) = -(rxt(k,329) + rxt(k,330)) * y(k,200)
         mat(k,350) = rxt(k,328)*y(k,217)
         mat(k,1619) = rxt(k,328)*y(k,16)
         mat(k,1775) = .750_r8*rxt(k,298)*y(k,202)
         mat(k,713) = .750_r8*rxt(k,298)*y(k,124)
         mat(k,714) = -(rxt(k,297)*y(k,203) + rxt(k,298)*y(k,124))
         mat(k,1902) = -rxt(k,297)*y(k,202)
         mat(k,1797) = -rxt(k,298)*y(k,202)
         mat(k,549) = rxt(k,304)*y(k,217)
         mat(k,1636) = rxt(k,304)*y(k,25)
         mat(k,1950) = -((rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,76) + rxt(k,139) &
                      *y(k,133) + rxt(k,140)*y(k,134) + rxt(k,144)*y(k,217) &
                      + 4._r8*rxt(k,149)*y(k,203) + rxt(k,159)*y(k,126) + rxt(k,164) &
                      *y(k,124) + rxt(k,169)*y(k,125) + (rxt(k,179) + rxt(k,180) &
                      ) * y(k,56) + rxt(k,186)*y(k,59) + rxt(k,212)*y(k,17) + rxt(k,218) &
                      *y(k,19) + rxt(k,255)*y(k,42) + rxt(k,261)*y(k,198) + rxt(k,268) &
                      *y(k,204) + rxt(k,282)*y(k,194) + rxt(k,293)*y(k,197) + rxt(k,297) &
                      *y(k,202) + rxt(k,310)*y(k,195) + rxt(k,318)*y(k,219) + rxt(k,322) &
                      *y(k,220) + rxt(k,334)*y(k,209) + rxt(k,343)*y(k,211) + rxt(k,347) &
                      *y(k,213) + rxt(k,357)*y(k,189) + rxt(k,367)*y(k,205) + rxt(k,372) &
                      *y(k,206) + rxt(k,381)*y(k,101) + rxt(k,392)*y(k,225) + rxt(k,396) &
                      *y(k,188) + rxt(k,399)*y(k,191) + rxt(k,403)*y(k,193) + rxt(k,406) &
                      *y(k,196) + rxt(k,410)*y(k,199) + rxt(k,413)*y(k,210) + rxt(k,416) &
                      *y(k,212) + rxt(k,419)*y(k,218) + rxt(k,426)*y(k,223) + rxt(k,432) &
                      *y(k,226) + rxt(k,435)*y(k,228) + rxt(k,446)*y(k,215) + rxt(k,451) &
                      *y(k,221) + rxt(k,456)*y(k,222))
         mat(k,1468) = -(rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,203)
         mat(k,2250) = -rxt(k,139)*y(k,203)
         mat(k,2128) = -rxt(k,140)*y(k,203)
         mat(k,1694) = -rxt(k,144)*y(k,203)
         mat(k,1751) = -rxt(k,159)*y(k,203)
         mat(k,1843) = -rxt(k,164)*y(k,203)
         mat(k,2195) = -rxt(k,169)*y(k,203)
         mat(k,2015) = -(rxt(k,179) + rxt(k,180)) * y(k,203)
         mat(k,1976) = -rxt(k,186)*y(k,203)
         mat(k,1420) = -rxt(k,212)*y(k,203)
         mat(k,2219) = -rxt(k,218)*y(k,203)
         mat(k,1490) = -rxt(k,255)*y(k,203)
         mat(k,2067) = -rxt(k,261)*y(k,203)
         mat(k,445) = -rxt(k,268)*y(k,203)
         mat(k,820) = -rxt(k,282)*y(k,203)
         mat(k,1393) = -rxt(k,293)*y(k,203)
         mat(k,719) = -rxt(k,297)*y(k,203)
         mat(k,797) = -rxt(k,310)*y(k,203)
         mat(k,774) = -rxt(k,318)*y(k,203)
         mat(k,1136) = -rxt(k,322)*y(k,203)
         mat(k,1319) = -rxt(k,334)*y(k,203)
         mat(k,1361) = -rxt(k,343)*y(k,203)
         mat(k,682) = -rxt(k,347)*y(k,203)
         mat(k,906) = -rxt(k,357)*y(k,203)
         mat(k,1266) = -rxt(k,367)*y(k,203)
         mat(k,1298) = -rxt(k,372)*y(k,203)
         mat(k,1203) = -rxt(k,381)*y(k,203)
         mat(k,1180) = -rxt(k,392)*y(k,203)
         mat(k,522) = -rxt(k,396)*y(k,203)
         mat(k,488) = -rxt(k,399)*y(k,203)
         mat(k,439) = -rxt(k,403)*y(k,203)
         mat(k,628) = -rxt(k,406)*y(k,203)
         mat(k,764) = -rxt(k,410)*y(k,203)
         mat(k,725) = -rxt(k,413)*y(k,203)
         mat(k,879) = -rxt(k,416)*y(k,203)
         mat(k,458) = -rxt(k,419)*y(k,203)
         mat(k,740) = -rxt(k,426)*y(k,203)
         mat(k,757) = -rxt(k,432)*y(k,203)
         mat(k,504) = -rxt(k,435)*y(k,203)
         mat(k,1054) = -rxt(k,446)*y(k,203)
         mat(k,1122) = -rxt(k,451)*y(k,203)
         mat(k,919) = -rxt(k,456)*y(k,203)
         mat(k,1013) = .570_r8*rxt(k,440)*y(k,134)
         mat(k,169) = .650_r8*rxt(k,398)*y(k,217)
         mat(k,1420) = mat(k,1420) + rxt(k,211)*y(k,42)
         mat(k,2219) = mat(k,2219) + rxt(k,223)*y(k,217)
         mat(k,290) = .350_r8*rxt(k,277)*y(k,217)
         mat(k,554) = .130_r8*rxt(k,279)*y(k,134)
         mat(k,264) = rxt(k,284)*y(k,217)
         mat(k,1035) = .280_r8*rxt(k,308)*y(k,134)
         mat(k,1490) = mat(k,1490) + rxt(k,211)*y(k,17) + rxt(k,175)*y(k,56) &
                      + rxt(k,256)*y(k,126) + rxt(k,257)*y(k,133)
         mat(k,598) = rxt(k,240)*y(k,56) + rxt(k,241)*y(k,217)
         mat(k,368) = rxt(k,243)*y(k,56) + rxt(k,244)*y(k,217)
         mat(k,104) = rxt(k,290)*y(k,217)
         mat(k,789) = rxt(k,263)*y(k,217)
         mat(k,1437) = rxt(k,272)*y(k,216)
         mat(k,2015) = mat(k,2015) + rxt(k,175)*y(k,42) + rxt(k,240)*y(k,43) &
                      + rxt(k,243)*y(k,46) + rxt(k,178)*y(k,79)
         mat(k,1976) = mat(k,1976) + rxt(k,182)*y(k,198) + rxt(k,193)*y(k,217)
         mat(k,1105) = rxt(k,275)*y(k,217)
         mat(k,198) = .730_r8*rxt(k,409)*y(k,217)
         mat(k,302) = .500_r8*rxt(k,477)*y(k,217)
         mat(k,1100) = rxt(k,301)*y(k,217)
         mat(k,982) = rxt(k,302)*y(k,217)
         mat(k,605) = rxt(k,178)*y(k,56) + rxt(k,134)*y(k,133) + rxt(k,143)*y(k,217)
         mat(k,182) = rxt(k,266)*y(k,217)
         mat(k,932) = rxt(k,267)*y(k,217)
         mat(k,1161) = rxt(k,331)*y(k,217)
         mat(k,1145) = rxt(k,316)*y(k,217)
         mat(k,854) = .370_r8*rxt(k,385)*y(k,134)
         mat(k,588) = .300_r8*rxt(k,376)*y(k,217)
         mat(k,563) = rxt(k,377)*y(k,217)
         mat(k,1203) = mat(k,1203) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) &
                      + rxt(k,379)*y(k,197) + 1.200_r8*rxt(k,380)*y(k,198)
         mat(k,401) = rxt(k,384)*y(k,217)
         mat(k,1216) = .140_r8*rxt(k,338)*y(k,134)
         mat(k,314) = .200_r8*rxt(k,340)*y(k,217)
         mat(k,579) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,963) = .570_r8*rxt(k,443)*y(k,134)
         mat(k,1341) = .280_r8*rxt(k,352)*y(k,134)
         mat(k,378) = rxt(k,388)*y(k,217)
         mat(k,1081) = rxt(k,389)*y(k,217)
         mat(k,1843) = mat(k,1843) + rxt(k,382)*y(k,101) + rxt(k,358)*y(k,189) &
                      + rxt(k,400)*y(k,191) + rxt(k,405)*y(k,193) + rxt(k,283) &
                      *y(k,194) + rxt(k,311)*y(k,195) + rxt(k,262)*y(k,198) &
                      + .170_r8*rxt(k,411)*y(k,199) + rxt(k,329)*y(k,200) &
                      + .250_r8*rxt(k,298)*y(k,202) + rxt(k,270)*y(k,204) &
                      + .920_r8*rxt(k,368)*y(k,205) + .920_r8*rxt(k,374)*y(k,206) &
                      + .470_r8*rxt(k,336)*y(k,209) + .400_r8*rxt(k,414)*y(k,210) &
                      + .830_r8*rxt(k,417)*y(k,212) + rxt(k,420)*y(k,218) + rxt(k,319) &
                      *y(k,219) + .900_r8*rxt(k,452)*y(k,221) + .800_r8*rxt(k,457) &
                      *y(k,222) + rxt(k,427)*y(k,223) + rxt(k,393)*y(k,225) &
                      + rxt(k,433)*y(k,226) + rxt(k,436)*y(k,228)
         mat(k,1751) = mat(k,1751) + rxt(k,256)*y(k,42) + rxt(k,383)*y(k,101) &
                      + rxt(k,369)*y(k,205) + rxt(k,375)*y(k,206) + .470_r8*rxt(k,335) &
                      *y(k,209) + rxt(k,162)*y(k,217) + rxt(k,394)*y(k,225)
         mat(k,2250) = mat(k,2250) + rxt(k,257)*y(k,42) + rxt(k,134)*y(k,79)
         mat(k,2128) = mat(k,2128) + .570_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .280_r8*rxt(k,308)*y(k,29) + .370_r8*rxt(k,385) &
                      *y(k,98) + .140_r8*rxt(k,338)*y(k,105) + .570_r8*rxt(k,443) &
                      *y(k,110) + .280_r8*rxt(k,352)*y(k,111) + rxt(k,146)*y(k,217)
         mat(k,178) = .800_r8*rxt(k,421)*y(k,217)
         mat(k,835) = rxt(k,467)*y(k,217)
         mat(k,1065) = .200_r8*rxt(k,461)*y(k,217)
         mat(k,193) = .280_r8*rxt(k,429)*y(k,217)
         mat(k,215) = .380_r8*rxt(k,431)*y(k,217)
         mat(k,220) = .630_r8*rxt(k,437)*y(k,217)
         mat(k,906) = mat(k,906) + rxt(k,358)*y(k,124)
         mat(k,488) = mat(k,488) + rxt(k,400)*y(k,124)
         mat(k,439) = mat(k,439) + rxt(k,405)*y(k,124)
         mat(k,820) = mat(k,820) + rxt(k,283)*y(k,124) + 2.400_r8*rxt(k,280)*y(k,194) &
                      + rxt(k,281)*y(k,198)
         mat(k,797) = mat(k,797) + rxt(k,311)*y(k,124) + rxt(k,309)*y(k,198)
         mat(k,1393) = mat(k,1393) + rxt(k,379)*y(k,101) + .900_r8*rxt(k,292)*y(k,198) &
                      + rxt(k,365)*y(k,205) + rxt(k,370)*y(k,206) + .470_r8*rxt(k,332) &
                      *y(k,209) + rxt(k,390)*y(k,225)
         mat(k,2067) = mat(k,2067) + rxt(k,182)*y(k,59) + 1.200_r8*rxt(k,380)*y(k,101) &
                      + rxt(k,262)*y(k,124) + rxt(k,281)*y(k,194) + rxt(k,309) &
                      *y(k,195) + .900_r8*rxt(k,292)*y(k,197) + 4.000_r8*rxt(k,259) &
                      *y(k,198) + rxt(k,366)*y(k,205) + rxt(k,371)*y(k,206) &
                      + .730_r8*rxt(k,333)*y(k,209) + rxt(k,342)*y(k,211) &
                      + .500_r8*rxt(k,445)*y(k,215) + .300_r8*rxt(k,321)*y(k,220) &
                      + rxt(k,450)*y(k,221) + rxt(k,455)*y(k,222) + .800_r8*rxt(k,391) &
                      *y(k,225)
         mat(k,764) = mat(k,764) + .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410) &
                      *y(k,203)
         mat(k,570) = rxt(k,329)*y(k,124)
         mat(k,719) = mat(k,719) + .250_r8*rxt(k,298)*y(k,124)
         mat(k,1950) = mat(k,1950) + .070_r8*rxt(k,410)*y(k,199) + .160_r8*rxt(k,413) &
                      *y(k,210) + .330_r8*rxt(k,416)*y(k,212)
         mat(k,445) = mat(k,445) + rxt(k,270)*y(k,124)
         mat(k,1266) = mat(k,1266) + .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) &
                      + rxt(k,365)*y(k,197) + rxt(k,366)*y(k,198)
         mat(k,1298) = mat(k,1298) + .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,197) + rxt(k,371)*y(k,198)
         mat(k,1319) = mat(k,1319) + .470_r8*rxt(k,336)*y(k,124) + .470_r8*rxt(k,335) &
                      *y(k,126) + .470_r8*rxt(k,332)*y(k,197) + .730_r8*rxt(k,333) &
                      *y(k,198)
         mat(k,725) = mat(k,725) + .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413) &
                      *y(k,203)
         mat(k,1361) = mat(k,1361) + rxt(k,342)*y(k,198)
         mat(k,879) = mat(k,879) + .830_r8*rxt(k,417)*y(k,124) + .330_r8*rxt(k,416) &
                      *y(k,203)
         mat(k,1054) = mat(k,1054) + .500_r8*rxt(k,445)*y(k,198)
         mat(k,1530) = rxt(k,272)*y(k,54)
         mat(k,1694) = mat(k,1694) + .650_r8*rxt(k,398)*y(k,7) + rxt(k,223)*y(k,19) &
                      + .350_r8*rxt(k,277)*y(k,24) + rxt(k,284)*y(k,26) + rxt(k,241) &
                      *y(k,43) + rxt(k,244)*y(k,46) + rxt(k,290)*y(k,47) + rxt(k,263) &
                      *y(k,52) + rxt(k,193)*y(k,59) + rxt(k,275)*y(k,62) &
                      + .730_r8*rxt(k,409)*y(k,66) + .500_r8*rxt(k,477)*y(k,67) &
                      + rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) + rxt(k,143)*y(k,79) &
                      + rxt(k,266)*y(k,86) + rxt(k,267)*y(k,87) + rxt(k,331)*y(k,93) &
                      + rxt(k,316)*y(k,95) + .300_r8*rxt(k,376)*y(k,99) + rxt(k,377) &
                      *y(k,100) + rxt(k,384)*y(k,102) + .200_r8*rxt(k,340)*y(k,106) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,388)*y(k,115) + rxt(k,389) &
                      *y(k,116) + rxt(k,162)*y(k,126) + rxt(k,146)*y(k,134) &
                      + .800_r8*rxt(k,421)*y(k,142) + rxt(k,467)*y(k,151) &
                      + .200_r8*rxt(k,461)*y(k,178) + .280_r8*rxt(k,429)*y(k,180) &
                      + .380_r8*rxt(k,431)*y(k,182) + .630_r8*rxt(k,437)*y(k,184)
         mat(k,458) = mat(k,458) + rxt(k,420)*y(k,124)
         mat(k,774) = mat(k,774) + rxt(k,319)*y(k,124)
         mat(k,1136) = mat(k,1136) + .300_r8*rxt(k,321)*y(k,198)
         mat(k,1122) = mat(k,1122) + .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,198)
         mat(k,919) = mat(k,919) + .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,198)
         mat(k,740) = mat(k,740) + rxt(k,427)*y(k,124)
         mat(k,1180) = mat(k,1180) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126) &
                      + rxt(k,390)*y(k,197) + .800_r8*rxt(k,391)*y(k,198)
         mat(k,757) = mat(k,757) + rxt(k,433)*y(k,124)
         mat(k,504) = mat(k,504) + rxt(k,436)*y(k,124)
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
         mat(k,441) = -(rxt(k,268)*y(k,203) + rxt(k,270)*y(k,124))
         mat(k,1883) = -rxt(k,268)*y(k,204)
         mat(k,1780) = -rxt(k,270)*y(k,204)
         mat(k,1476) = rxt(k,255)*y(k,203)
         mat(k,1883) = mat(k,1883) + rxt(k,255)*y(k,42)
         mat(k,1256) = -(rxt(k,365)*y(k,197) + rxt(k,366)*y(k,198) + rxt(k,367) &
                      *y(k,203) + rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126))
         mat(k,1381) = -rxt(k,365)*y(k,205)
         mat(k,2053) = -rxt(k,366)*y(k,205)
         mat(k,1934) = -rxt(k,367)*y(k,205)
         mat(k,1828) = -rxt(k,368)*y(k,205)
         mat(k,1736) = -rxt(k,369)*y(k,205)
         mat(k,847) = .600_r8*rxt(k,386)*y(k,217)
         mat(k,1678) = .600_r8*rxt(k,386)*y(k,98)
         mat(k,1288) = -(rxt(k,370)*y(k,197) + rxt(k,371)*y(k,198) + rxt(k,372) &
                      *y(k,203) + rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126))
         mat(k,1382) = -rxt(k,370)*y(k,206)
         mat(k,2054) = -rxt(k,371)*y(k,206)
         mat(k,1935) = -rxt(k,372)*y(k,206)
         mat(k,1829) = -rxt(k,374)*y(k,206)
         mat(k,1737) = -rxt(k,375)*y(k,206)
         mat(k,848) = .400_r8*rxt(k,386)*y(k,217)
         mat(k,1679) = .400_r8*rxt(k,386)*y(k,98)
         mat(k,66) = -(rxt(k,503)*y(k,203) + rxt(k,504)*y(k,124))
         mat(k,1860) = -rxt(k,503)*y(k,207)
         mat(k,1768) = -rxt(k,504)*y(k,207)
         mat(k,840) = rxt(k,506)*y(k,217)
         mat(k,1548) = rxt(k,506)*y(k,98)
         mat(k,72) = -(rxt(k,507)*y(k,203) + rxt(k,508)*y(k,124))
         mat(k,1861) = -rxt(k,507)*y(k,208)
         mat(k,1769) = -rxt(k,508)*y(k,208)
         mat(k,73) = rxt(k,509)*y(k,217)
         mat(k,1549) = rxt(k,509)*y(k,104)
         mat(k,1311) = -(rxt(k,332)*y(k,197) + rxt(k,333)*y(k,198) + rxt(k,334) &
                      *y(k,203) + rxt(k,335)*y(k,126) + (rxt(k,336) + rxt(k,337) &
                      ) * y(k,124))
         mat(k,1383) = -rxt(k,332)*y(k,209)
         mat(k,2055) = -rxt(k,333)*y(k,209)
         mat(k,1936) = -rxt(k,334)*y(k,209)
         mat(k,1738) = -rxt(k,335)*y(k,209)
         mat(k,1830) = -(rxt(k,336) + rxt(k,337)) * y(k,209)
         mat(k,1211) = .500_r8*rxt(k,339)*y(k,217)
         mat(k,311) = .200_r8*rxt(k,340)*y(k,217)
         mat(k,1331) = rxt(k,353)*y(k,217)
         mat(k,1680) = .500_r8*rxt(k,339)*y(k,105) + .200_r8*rxt(k,340)*y(k,106) &
                      + rxt(k,353)*y(k,111)
         mat(k,721) = -(rxt(k,413)*y(k,203) + rxt(k,414)*y(k,124) + rxt(k,415) &
                      *y(k,125))
         mat(k,1903) = -rxt(k,413)*y(k,210)
         mat(k,1798) = -rxt(k,414)*y(k,210)
         mat(k,2172) = -rxt(k,415)*y(k,210)
         mat(k,1354) = -(rxt(k,341)*y(k,197) + rxt(k,342)*y(k,198) + rxt(k,343) &
                      *y(k,203) + 4._r8*rxt(k,344)*y(k,211) + rxt(k,345)*y(k,124) &
                      + rxt(k,346)*y(k,126) + rxt(k,354)*y(k,125))
         mat(k,1385) = -rxt(k,341)*y(k,211)
         mat(k,2057) = -rxt(k,342)*y(k,211)
         mat(k,1938) = -rxt(k,343)*y(k,211)
         mat(k,1832) = -rxt(k,345)*y(k,211)
         mat(k,1740) = -rxt(k,346)*y(k,211)
         mat(k,2184) = -rxt(k,354)*y(k,211)
         mat(k,1212) = .500_r8*rxt(k,339)*y(k,217)
         mat(k,312) = .500_r8*rxt(k,340)*y(k,217)
         mat(k,1682) = .500_r8*rxt(k,339)*y(k,105) + .500_r8*rxt(k,340)*y(k,106)
         mat(k,873) = -(rxt(k,416)*y(k,203) + rxt(k,417)*y(k,124) + rxt(k,418) &
                      *y(k,125))
         mat(k,1915) = -rxt(k,416)*y(k,212)
         mat(k,1807) = -rxt(k,417)*y(k,212)
         mat(k,2177) = -rxt(k,418)*y(k,212)
         mat(k,677) = -(rxt(k,347)*y(k,203) + rxt(k,348)*y(k,124))
         mat(k,1898) = -rxt(k,347)*y(k,213)
         mat(k,1796) = -rxt(k,348)*y(k,213)
         mat(k,507) = rxt(k,349)*y(k,217)
         mat(k,316) = rxt(k,350)*y(k,217)
         mat(k,1632) = rxt(k,349)*y(k,107) + rxt(k,350)*y(k,108)
         mat(k,80) = -(rxt(k,511)*y(k,203) + rxt(k,512)*y(k,124))
         mat(k,1862) = -rxt(k,511)*y(k,214)
         mat(k,1770) = -rxt(k,512)*y(k,214)
         mat(k,940) = rxt(k,514)*y(k,217)
         mat(k,1551) = rxt(k,514)*y(k,110)
         mat(k,1045) = -(rxt(k,445)*y(k,198) + rxt(k,446)*y(k,203) + rxt(k,447) &
                      *y(k,124) + rxt(k,448)*y(k,126))
         mat(k,2040) = -rxt(k,445)*y(k,215)
         mat(k,1922) = -rxt(k,446)*y(k,215)
         mat(k,1814) = -rxt(k,447)*y(k,215)
         mat(k,1721) = -rxt(k,448)*y(k,215)
         mat(k,1000) = rxt(k,439)*y(k,126)
         mat(k,951) = rxt(k,442)*y(k,126)
         mat(k,1721) = mat(k,1721) + rxt(k,439)*y(k,6) + rxt(k,442)*y(k,110) &
                      + .500_r8*rxt(k,459)*y(k,177)
         mat(k,393) = rxt(k,449)*y(k,217)
         mat(k,970) = .500_r8*rxt(k,459)*y(k,126)
         mat(k,1663) = rxt(k,449)*y(k,128)
         mat(k,1526) = -(rxt(k,125)*y(k,77) + rxt(k,126)*y(k,229) + rxt(k,129) &
                      *y(k,134) + (rxt(k,167) + rxt(k,168)) * y(k,113) + rxt(k,200) &
                      *y(k,33) + rxt(k,201)*y(k,34) + rxt(k,202)*y(k,36) + rxt(k,203) &
                      *y(k,37) + rxt(k,204)*y(k,38) + rxt(k,205)*y(k,39) + rxt(k,206) &
                      *y(k,40) + (rxt(k,207) + rxt(k,208)) * y(k,85) + rxt(k,227) &
                      *y(k,35) + rxt(k,228)*y(k,55) + rxt(k,229)*y(k,78) + (rxt(k,230) &
                      + rxt(k,231)) * y(k,81) + rxt(k,236)*y(k,64) + rxt(k,237) &
                      *y(k,65) + rxt(k,250)*y(k,41) + rxt(k,251)*y(k,43) + rxt(k,252) &
                      *y(k,82) + rxt(k,253)*y(k,83) + rxt(k,254)*y(k,84) + (rxt(k,271) &
                      + rxt(k,272) + rxt(k,273)) * y(k,54) + rxt(k,274)*y(k,86))
         mat(k,1405) = -rxt(k,125)*y(k,216)
         mat(k,2272) = -rxt(k,126)*y(k,216)
         mat(k,2124) = -rxt(k,129)*y(k,216)
         mat(k,184) = -(rxt(k,167) + rxt(k,168)) * y(k,216)
         mat(k,100) = -rxt(k,200)*y(k,216)
         mat(k,144) = -rxt(k,201)*y(k,216)
         mat(k,115) = -rxt(k,202)*y(k,216)
         mat(k,154) = -rxt(k,203)*y(k,216)
         mat(k,119) = -rxt(k,204)*y(k,216)
         mat(k,159) = -rxt(k,205)*y(k,216)
         mat(k,123) = -rxt(k,206)*y(k,216)
         mat(k,2147) = -(rxt(k,207) + rxt(k,208)) * y(k,216)
         mat(k,150) = -rxt(k,227)*y(k,216)
         mat(k,449) = -rxt(k,228)*y(k,216)
         mat(k,108) = -rxt(k,229)*y(k,216)
         mat(k,807) = -(rxt(k,230) + rxt(k,231)) * y(k,216)
         mat(k,237) = -rxt(k,236)*y(k,216)
         mat(k,228) = -rxt(k,237)*y(k,216)
         mat(k,469) = -rxt(k,250)*y(k,216)
         mat(k,596) = -rxt(k,251)*y(k,216)
         mat(k,223) = -rxt(k,252)*y(k,216)
         mat(k,250) = -rxt(k,253)*y(k,216)
         mat(k,306) = -rxt(k,254)*y(k,216)
         mat(k,1434) = -(rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,216)
         mat(k,180) = -rxt(k,274)*y(k,216)
         mat(k,1691) = -(rxt(k,142)*y(k,77) + rxt(k,143)*y(k,79) + rxt(k,144)*y(k,203) &
                      + rxt(k,145)*y(k,133) + rxt(k,146)*y(k,134) + (4._r8*rxt(k,147) &
                      + 4._r8*rxt(k,148)) * y(k,217) + rxt(k,150)*y(k,90) + rxt(k,162) &
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
                      *y(k,139) + rxt(k,304)*y(k,25) + rxt(k,312)*y(k,30) + rxt(k,313) &
                      *y(k,31) + rxt(k,315)*y(k,49) + rxt(k,316)*y(k,95) + rxt(k,317) &
                      *y(k,127) + rxt(k,320)*y(k,146) + rxt(k,324)*y(k,147) + rxt(k,325) &
                      *y(k,29) + rxt(k,326)*y(k,48) + rxt(k,328)*y(k,16) + rxt(k,331) &
                      *y(k,93) + rxt(k,339)*y(k,105) + rxt(k,340)*y(k,106) + rxt(k,349) &
                      *y(k,107) + rxt(k,350)*y(k,108) + rxt(k,351)*y(k,109) + rxt(k,353) &
                      *y(k,111) + rxt(k,356)*y(k,1) + rxt(k,360)*y(k,2) + rxt(k,361) &
                      *y(k,15) + rxt(k,362)*y(k,94) + rxt(k,363)*y(k,96) + rxt(k,364) &
                      *y(k,97) + rxt(k,376)*y(k,99) + rxt(k,377)*y(k,100) + rxt(k,384) &
                      *y(k,102) + rxt(k,386)*y(k,98) + rxt(k,387)*y(k,103) + rxt(k,388) &
                      *y(k,115) + rxt(k,389)*y(k,116) + rxt(k,395)*y(k,181) + rxt(k,398) &
                      *y(k,7) + rxt(k,401)*y(k,8) + rxt(k,402)*y(k,22) + rxt(k,404) &
                      *y(k,23) + rxt(k,408)*y(k,32) + rxt(k,409)*y(k,66) + rxt(k,421) &
                      *y(k,142) + rxt(k,424)*y(k,143) + rxt(k,428)*y(k,179) + rxt(k,429) &
                      *y(k,180) + rxt(k,431)*y(k,182) + rxt(k,434)*y(k,183) + rxt(k,437) &
                      *y(k,184) + rxt(k,438)*y(k,185) + rxt(k,441)*y(k,6) + rxt(k,444) &
                      *y(k,110) + rxt(k,449)*y(k,128) + rxt(k,453)*y(k,174) + rxt(k,454) &
                      *y(k,175) + rxt(k,458)*y(k,176) + rxt(k,460)*y(k,177) + rxt(k,461) &
                      *y(k,178) + (rxt(k,463) + rxt(k,477)) * y(k,67) + rxt(k,465) &
                      *y(k,137) + rxt(k,467)*y(k,151) + rxt(k,471)*y(k,148) + rxt(k,476) &
                      *y(k,150) + rxt(k,479)*y(k,120))
         mat(k,1406) = -rxt(k,142)*y(k,217)
         mat(k,604) = -rxt(k,143)*y(k,217)
         mat(k,1947) = -rxt(k,144)*y(k,217)
         mat(k,2247) = -rxt(k,145)*y(k,217)
         mat(k,2125) = -rxt(k,146)*y(k,217)
         mat(k,404) = -rxt(k,150)*y(k,217)
         mat(k,1748) = -rxt(k,162)*y(k,217)
         mat(k,494) = -rxt(k,163)*y(k,217)
         mat(k,2192) = -rxt(k,171)*y(k,217)
         mat(k,1451) = -rxt(k,172)*y(k,217)
         mat(k,886) = -rxt(k,191)*y(k,217)
         mat(k,1973) = -(rxt(k,193) + rxt(k,194)) * y(k,217)
         mat(k,2148) = -rxt(k,196)*y(k,217)
         mat(k,825) = -rxt(k,199)*y(k,217)
         mat(k,2216) = -rxt(k,223)*y(k,217)
         mat(k,808) = -rxt(k,225)*y(k,217)
         mat(k,470) = -rxt(k,239)*y(k,217)
         mat(k,597) = -rxt(k,241)*y(k,217)
         mat(k,126) = -rxt(k,242)*y(k,217)
         mat(k,367) = -rxt(k,244)*y(k,217)
         mat(k,450) = -rxt(k,246)*y(k,217)
         mat(k,224) = -rxt(k,247)*y(k,217)
         mat(k,251) = -rxt(k,248)*y(k,217)
         mat(k,307) = -rxt(k,249)*y(k,217)
         mat(k,1487) = -rxt(k,258)*y(k,217)
         mat(k,788) = -rxt(k,263)*y(k,217)
         mat(k,388) = -rxt(k,264)*y(k,217)
         mat(k,1435) = -rxt(k,265)*y(k,217)
         mat(k,181) = -rxt(k,266)*y(k,217)
         mat(k,931) = -rxt(k,267)*y(k,217)
         mat(k,1104) = -rxt(k,275)*y(k,217)
         mat(k,289) = -rxt(k,277)*y(k,217)
         mat(k,263) = -rxt(k,284)*y(k,217)
         mat(k,347) = -rxt(k,285)*y(k,217)
         mat(k,293) = -rxt(k,287)*y(k,217)
         mat(k,1089) = -rxt(k,289)*y(k,217)
         mat(k,103) = -rxt(k,290)*y(k,217)
         mat(k,686) = -rxt(k,295)*y(k,217)
         mat(k,614) = -rxt(k,296)*y(k,217)
         mat(k,1099) = -rxt(k,301)*y(k,217)
         mat(k,981) = -rxt(k,302)*y(k,217)
         mat(k,528) = -rxt(k,303)*y(k,217)
         mat(k,553) = -rxt(k,304)*y(k,217)
         mat(k,418) = -rxt(k,312)*y(k,217)
         mat(k,111) = -rxt(k,313)*y(k,217)
         mat(k,1224) = -rxt(k,315)*y(k,217)
         mat(k,1144) = -rxt(k,316)*y(k,217)
         mat(k,861) = -rxt(k,317)*y(k,217)
         mat(k,537) = -rxt(k,320)*y(k,217)
         mat(k,413) = -rxt(k,324)*y(k,217)
         mat(k,1032) = -rxt(k,325)*y(k,217)
         mat(k,925) = -rxt(k,326)*y(k,217)
         mat(k,354) = -rxt(k,328)*y(k,217)
         mat(k,1158) = -rxt(k,331)*y(k,217)
         mat(k,1215) = -rxt(k,339)*y(k,217)
         mat(k,313) = -rxt(k,340)*y(k,217)
         mat(k,510) = -rxt(k,349)*y(k,217)
         mat(k,319) = -rxt(k,350)*y(k,217)
         mat(k,577) = -rxt(k,351)*y(k,217)
         mat(k,1338) = -rxt(k,353)*y(k,217)
         mat(k,673) = -rxt(k,356)*y(k,217)
         mat(k,640) = -rxt(k,360)*y(k,217)
         mat(k,240) = -rxt(k,361)*y(k,217)
         mat(k,233) = -rxt(k,362)*y(k,217)
         mat(k,322) = -rxt(k,363)*y(k,217)
         mat(k,137) = -rxt(k,364)*y(k,217)
         mat(k,587) = -rxt(k,376)*y(k,217)
         mat(k,562) = -rxt(k,377)*y(k,217)
         mat(k,400) = -rxt(k,384)*y(k,217)
         mat(k,852) = -rxt(k,386)*y(k,217)
         mat(k,695) = -rxt(k,387)*y(k,217)
         mat(k,377) = -rxt(k,388)*y(k,217)
         mat(k,1079) = -rxt(k,389)*y(k,217)
         mat(k,205) = -rxt(k,395)*y(k,217)
         mat(k,168) = -rxt(k,398)*y(k,217)
         mat(k,425) = -rxt(k,401)*y(k,217)
         mat(k,255) = -rxt(k,402)*y(k,217)
         mat(k,342) = -rxt(k,404)*y(k,217)
         mat(k,268) = -rxt(k,408)*y(k,217)
         mat(k,197) = -rxt(k,409)*y(k,217)
         mat(k,177) = -rxt(k,421)*y(k,217)
         mat(k,336) = -rxt(k,424)*y(k,217)
         mat(k,663) = -rxt(k,428)*y(k,217)
         mat(k,192) = -rxt(k,429)*y(k,217)
         mat(k,214) = -rxt(k,431)*y(k,217)
         mat(k,710) = -rxt(k,434)*y(k,217)
         mat(k,219) = -rxt(k,437)*y(k,217)
         mat(k,431) = -rxt(k,438)*y(k,217)
         mat(k,1010) = -rxt(k,441)*y(k,217)
         mat(k,960) = -rxt(k,444)*y(k,217)
         mat(k,395) = -rxt(k,449)*y(k,217)
         mat(k,650) = -rxt(k,453)*y(k,217)
         mat(k,620) = -rxt(k,454)*y(k,217)
         mat(k,479) = -rxt(k,458)*y(k,217)
         mat(k,974) = -rxt(k,460)*y(k,217)
         mat(k,1064) = -rxt(k,461)*y(k,217)
         mat(k,300) = -(rxt(k,463) + rxt(k,477)) * y(k,217)
         mat(k,363) = -rxt(k,465)*y(k,217)
         mat(k,834) = -rxt(k,467)*y(k,217)
         mat(k,514) = -rxt(k,471)*y(k,217)
         mat(k,1235) = -rxt(k,476)*y(k,217)
         mat(k,97) = -rxt(k,479)*y(k,217)
         mat(k,1010) = mat(k,1010) + .630_r8*rxt(k,440)*y(k,134)
         mat(k,289) = mat(k,289) + .650_r8*rxt(k,277)*y(k,217)
         mat(k,553) = mat(k,553) + .130_r8*rxt(k,279)*y(k,134)
         mat(k,347) = mat(k,347) + .500_r8*rxt(k,285)*y(k,217)
         mat(k,1032) = mat(k,1032) + .360_r8*rxt(k,308)*y(k,134)
         mat(k,1487) = mat(k,1487) + rxt(k,257)*y(k,133)
         mat(k,388) = mat(k,388) + .300_r8*rxt(k,264)*y(k,217)
         mat(k,1435) = mat(k,1435) + rxt(k,271)*y(k,216)
         mat(k,2012) = rxt(k,180)*y(k,203)
         mat(k,869) = rxt(k,234)*y(k,229)
         mat(k,1466) = rxt(k,141)*y(k,134) + 2.000_r8*rxt(k,136)*y(k,203)
         mat(k,1406) = mat(k,1406) + rxt(k,133)*y(k,133) + rxt(k,125)*y(k,216)
         mat(k,604) = mat(k,604) + rxt(k,134)*y(k,133)
         mat(k,808) = mat(k,808) + rxt(k,224)*y(k,133) + rxt(k,230)*y(k,216)
         mat(k,2148) = mat(k,2148) + rxt(k,195)*y(k,133) + rxt(k,207)*y(k,216)
         mat(k,181) = mat(k,181) + rxt(k,274)*y(k,216)
         mat(k,780) = rxt(k,226)*y(k,133)
         mat(k,825) = mat(k,825) + rxt(k,198)*y(k,133)
         mat(k,852) = mat(k,852) + .320_r8*rxt(k,385)*y(k,134)
         mat(k,695) = mat(k,695) + .600_r8*rxt(k,387)*y(k,217)
         mat(k,1215) = mat(k,1215) + .240_r8*rxt(k,338)*y(k,134)
         mat(k,313) = mat(k,313) + .100_r8*rxt(k,340)*y(k,217)
         mat(k,960) = mat(k,960) + .630_r8*rxt(k,443)*y(k,134)
         mat(k,1338) = mat(k,1338) + .360_r8*rxt(k,352)*y(k,134)
         mat(k,1840) = rxt(k,164)*y(k,203)
         mat(k,1748) = mat(k,1748) + rxt(k,159)*y(k,203)
         mat(k,2247) = mat(k,2247) + rxt(k,257)*y(k,42) + rxt(k,133)*y(k,77) &
                      + rxt(k,134)*y(k,79) + rxt(k,224)*y(k,81) + rxt(k,195)*y(k,85) &
                      + rxt(k,226)*y(k,91) + rxt(k,198)*y(k,92) + rxt(k,139)*y(k,203)
         mat(k,2125) = mat(k,2125) + .630_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .360_r8*rxt(k,308)*y(k,29) + rxt(k,141)*y(k,76) &
                      + .320_r8*rxt(k,385)*y(k,98) + .240_r8*rxt(k,338)*y(k,105) &
                      + .630_r8*rxt(k,443)*y(k,110) + .360_r8*rxt(k,352)*y(k,111) &
                      + rxt(k,140)*y(k,203)
         mat(k,537) = mat(k,537) + .500_r8*rxt(k,320)*y(k,217)
         mat(k,205) = mat(k,205) + .500_r8*rxt(k,395)*y(k,217)
         mat(k,520) = .400_r8*rxt(k,396)*y(k,203)
         mat(k,1390) = .450_r8*rxt(k,293)*y(k,203)
         mat(k,762) = .400_r8*rxt(k,410)*y(k,203)
         mat(k,1947) = mat(k,1947) + rxt(k,180)*y(k,56) + 2.000_r8*rxt(k,136)*y(k,76) &
                      + rxt(k,164)*y(k,124) + rxt(k,159)*y(k,126) + rxt(k,139) &
                      *y(k,133) + rxt(k,140)*y(k,134) + .400_r8*rxt(k,396)*y(k,188) &
                      + .450_r8*rxt(k,293)*y(k,197) + .400_r8*rxt(k,410)*y(k,199) &
                      + .450_r8*rxt(k,343)*y(k,211) + .400_r8*rxt(k,416)*y(k,212) &
                      + .200_r8*rxt(k,347)*y(k,213) + .150_r8*rxt(k,322)*y(k,220)
         mat(k,1358) = .450_r8*rxt(k,343)*y(k,203)
         mat(k,877) = .400_r8*rxt(k,416)*y(k,203)
         mat(k,680) = .200_r8*rxt(k,347)*y(k,203)
         mat(k,1527) = rxt(k,271)*y(k,54) + rxt(k,125)*y(k,77) + rxt(k,230)*y(k,81) &
                      + rxt(k,207)*y(k,85) + rxt(k,274)*y(k,86) + 2.000_r8*rxt(k,126) &
                      *y(k,229)
         mat(k,1691) = mat(k,1691) + .650_r8*rxt(k,277)*y(k,24) + .500_r8*rxt(k,285) &
                      *y(k,27) + .300_r8*rxt(k,264)*y(k,53) + .600_r8*rxt(k,387) &
                      *y(k,103) + .100_r8*rxt(k,340)*y(k,106) + .500_r8*rxt(k,320) &
                      *y(k,146) + .500_r8*rxt(k,395)*y(k,181)
         mat(k,1134) = .150_r8*rxt(k,322)*y(k,203)
         mat(k,2273) = rxt(k,234)*y(k,73) + 2.000_r8*rxt(k,126)*y(k,216)
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
         mat(k,454) = -(rxt(k,419)*y(k,203) + rxt(k,420)*y(k,124))
         mat(k,1884) = -rxt(k,419)*y(k,218)
         mat(k,1781) = -rxt(k,420)*y(k,218)
         mat(k,195) = .200_r8*rxt(k,409)*y(k,217)
         mat(k,175) = .140_r8*rxt(k,421)*y(k,217)
         mat(k,334) = rxt(k,424)*y(k,217)
         mat(k,1604) = .200_r8*rxt(k,409)*y(k,66) + .140_r8*rxt(k,421)*y(k,142) &
                      + rxt(k,424)*y(k,143)
         mat(k,768) = -(rxt(k,318)*y(k,203) + rxt(k,319)*y(k,124))
         mat(k,1907) = -rxt(k,318)*y(k,219)
         mat(k,1802) = -rxt(k,319)*y(k,219)
         mat(k,1020) = rxt(k,325)*y(k,217)
         mat(k,533) = .500_r8*rxt(k,320)*y(k,217)
         mat(k,1641) = rxt(k,325)*y(k,29) + .500_r8*rxt(k,320)*y(k,146)
         mat(k,1129) = -(rxt(k,321)*y(k,198) + rxt(k,322)*y(k,203) + rxt(k,323) &
                      *y(k,124))
         mat(k,2047) = -rxt(k,321)*y(k,220)
         mat(k,1928) = -rxt(k,322)*y(k,220)
         mat(k,1821) = -rxt(k,323)*y(k,220)
         mat(k,1005) = .060_r8*rxt(k,440)*y(k,134)
         mat(k,923) = rxt(k,326)*y(k,217)
         mat(k,955) = .060_r8*rxt(k,443)*y(k,134)
         mat(k,2107) = .060_r8*rxt(k,440)*y(k,6) + .060_r8*rxt(k,443)*y(k,110)
         mat(k,410) = rxt(k,324)*y(k,217)
         mat(k,1061) = .150_r8*rxt(k,461)*y(k,217)
         mat(k,1670) = rxt(k,326)*y(k,48) + rxt(k,324)*y(k,147) + .150_r8*rxt(k,461) &
                      *y(k,178)
         mat(k,1115) = -(rxt(k,450)*y(k,198) + rxt(k,451)*y(k,203) + rxt(k,452) &
                      *y(k,124))
         mat(k,2046) = -rxt(k,450)*y(k,221)
         mat(k,1927) = -rxt(k,451)*y(k,221)
         mat(k,1820) = -rxt(k,452)*y(k,221)
         mat(k,1727) = .500_r8*rxt(k,459)*y(k,177)
         mat(k,648) = rxt(k,453)*y(k,217)
         mat(k,973) = .500_r8*rxt(k,459)*y(k,126) + rxt(k,460)*y(k,217)
         mat(k,1669) = rxt(k,453)*y(k,174) + rxt(k,460)*y(k,177)
         mat(k,912) = -(rxt(k,455)*y(k,198) + rxt(k,456)*y(k,203) + rxt(k,457) &
                      *y(k,124))
         mat(k,2036) = -rxt(k,455)*y(k,222)
         mat(k,1917) = -rxt(k,456)*y(k,222)
         mat(k,1809) = -rxt(k,457)*y(k,222)
         mat(k,994) = rxt(k,441)*y(k,217)
         mat(k,945) = rxt(k,444)*y(k,217)
         mat(k,475) = rxt(k,458)*y(k,217)
         mat(k,1655) = rxt(k,441)*y(k,6) + rxt(k,444)*y(k,110) + rxt(k,458)*y(k,176)
         mat(k,732) = -(rxt(k,426)*y(k,203) + rxt(k,427)*y(k,124))
         mat(k,1904) = -rxt(k,426)*y(k,223)
         mat(k,1799) = -rxt(k,427)*y(k,223)
         mat(k,657) = rxt(k,428)*y(k,217)
         mat(k,191) = .650_r8*rxt(k,429)*y(k,217)
         mat(k,1638) = rxt(k,428)*y(k,179) + .650_r8*rxt(k,429)*y(k,180)
         mat(k,86) = -(rxt(k,517)*y(k,203) + rxt(k,518)*y(k,124))
         mat(k,1863) = -rxt(k,517)*y(k,224)
         mat(k,1771) = -rxt(k,518)*y(k,224)
         mat(k,186) = rxt(k,516)*y(k,217)
         mat(k,1552) = rxt(k,516)*y(k,180)
         mat(k,1173) = -(rxt(k,390)*y(k,197) + rxt(k,391)*y(k,198) + rxt(k,392) &
                      *y(k,203) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126))
         mat(k,1377) = -rxt(k,390)*y(k,225)
         mat(k,2049) = -rxt(k,391)*y(k,225)
         mat(k,1930) = -rxt(k,392)*y(k,225)
         mat(k,1824) = -rxt(k,393)*y(k,225)
         mat(k,1731) = -rxt(k,394)*y(k,225)
         mat(k,232) = rxt(k,362)*y(k,217)
         mat(k,321) = rxt(k,363)*y(k,217)
         mat(k,136) = rxt(k,364)*y(k,217)
         mat(k,691) = .400_r8*rxt(k,387)*y(k,217)
         mat(k,204) = .500_r8*rxt(k,395)*y(k,217)
         mat(k,1673) = rxt(k,362)*y(k,94) + rxt(k,363)*y(k,96) + rxt(k,364)*y(k,97) &
                      + .400_r8*rxt(k,387)*y(k,103) + .500_r8*rxt(k,395)*y(k,181)
         mat(k,748) = -(rxt(k,432)*y(k,203) + rxt(k,433)*y(k,124))
         mat(k,1905) = -rxt(k,432)*y(k,226)
         mat(k,1800) = -rxt(k,433)*y(k,226)
         mat(k,211) = .560_r8*rxt(k,431)*y(k,217)
         mat(k,703) = rxt(k,434)*y(k,217)
         mat(k,1639) = .560_r8*rxt(k,431)*y(k,182) + rxt(k,434)*y(k,183)
         mat(k,92) = -(rxt(k,520)*y(k,203) + rxt(k,521)*y(k,124))
         mat(k,1864) = -rxt(k,520)*y(k,227)
         mat(k,1772) = -rxt(k,521)*y(k,227)
         mat(k,206) = rxt(k,519)*y(k,217)
         mat(k,1553) = rxt(k,519)*y(k,182)
         mat(k,499) = -(rxt(k,435)*y(k,203) + rxt(k,436)*y(k,124))
         mat(k,1889) = -rxt(k,435)*y(k,228)
         mat(k,1786) = -rxt(k,436)*y(k,228)
         mat(k,218) = .300_r8*rxt(k,437)*y(k,217)
         mat(k,428) = rxt(k,438)*y(k,217)
         mat(k,1611) = .300_r8*rxt(k,437)*y(k,184) + rxt(k,438)*y(k,185)
         mat(k,2285) = -(rxt(k,126)*y(k,216) + rxt(k,234)*y(k,73) + rxt(k,478) &
                      *y(k,152))
         mat(k,1539) = -rxt(k,126)*y(k,229)
         mat(k,872) = -rxt(k,234)*y(k,229)
         mat(k,260) = -rxt(k,478)*y(k,229)
         mat(k,296) = rxt(k,287)*y(k,217)
         mat(k,420) = rxt(k,312)*y(k,217)
         mat(k,112) = rxt(k,313)*y(k,217)
         mat(k,473) = rxt(k,239)*y(k,217)
         mat(k,1498) = rxt(k,258)*y(k,217)
         mat(k,602) = rxt(k,241)*y(k,217)
         mat(k,128) = rxt(k,242)*y(k,217)
         mat(k,1093) = rxt(k,289)*y(k,217)
         mat(k,372) = rxt(k,244)*y(k,217)
         mat(k,927) = rxt(k,326)*y(k,217)
         mat(k,1228) = rxt(k,315)*y(k,217)
         mat(k,688) = rxt(k,295)*y(k,217)
         mat(k,616) = rxt(k,296)*y(k,217)
         mat(k,390) = rxt(k,264)*y(k,217)
         mat(k,1442) = rxt(k,265)*y(k,217)
         mat(k,1475) = rxt(k,137)*y(k,203)
         mat(k,1412) = rxt(k,142)*y(k,217)
         mat(k,609) = rxt(k,143)*y(k,217)
         mat(k,811) = rxt(k,225)*y(k,217)
         mat(k,309) = rxt(k,249)*y(k,217)
         mat(k,2160) = (rxt(k,530)+rxt(k,535))*y(k,91) + (rxt(k,523)+rxt(k,529) &
                       +rxt(k,534))*y(k,92) + rxt(k,196)*y(k,217)
         mat(k,934) = rxt(k,267)*y(k,217)
         mat(k,1459) = rxt(k,172)*y(k,217)
         mat(k,408) = rxt(k,150)*y(k,217)
         mat(k,785) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,830) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85) + rxt(k,199)*y(k,217)
         mat(k,1219) = .500_r8*rxt(k,339)*y(k,217)
         mat(k,98) = rxt(k,479)*y(k,217)
         mat(k,539) = rxt(k,320)*y(k,217)
         mat(k,414) = rxt(k,324)*y(k,217)
         mat(k,1959) = rxt(k,137)*y(k,76) + rxt(k,144)*y(k,217)
         mat(k,1703) = rxt(k,287)*y(k,28) + rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31) &
                      + rxt(k,239)*y(k,41) + rxt(k,258)*y(k,42) + rxt(k,241)*y(k,43) &
                      + rxt(k,242)*y(k,44) + rxt(k,289)*y(k,45) + rxt(k,244)*y(k,46) &
                      + rxt(k,326)*y(k,48) + rxt(k,315)*y(k,49) + rxt(k,295)*y(k,50) &
                      + rxt(k,296)*y(k,51) + rxt(k,264)*y(k,53) + rxt(k,265)*y(k,54) &
                      + rxt(k,142)*y(k,77) + rxt(k,143)*y(k,79) + rxt(k,225)*y(k,81) &
                      + rxt(k,249)*y(k,84) + rxt(k,196)*y(k,85) + rxt(k,267)*y(k,87) &
                      + rxt(k,172)*y(k,89) + rxt(k,150)*y(k,90) + rxt(k,199)*y(k,92) &
                      + .500_r8*rxt(k,339)*y(k,105) + rxt(k,479)*y(k,120) + rxt(k,320) &
                      *y(k,146) + rxt(k,324)*y(k,147) + rxt(k,144)*y(k,203) &
                      + 2.000_r8*rxt(k,147)*y(k,217)
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
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 138) = lmat(k, 138)
         mat(k, 139) = lmat(k, 139)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = lmat(k, 162)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 170) = lmat(k, 170)
         mat(k, 171) = lmat(k, 171)
         mat(k, 172) = lmat(k, 172)
         mat(k, 173) = lmat(k, 173)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
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
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 230) = mat(k, 230) + lmat(k, 230)
         mat(k, 231) = lmat(k, 231)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 241) = lmat(k, 241)
         mat(k, 242) = lmat(k, 242)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 258) = lmat(k, 258)
         mat(k, 259) = lmat(k, 259)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 266) = lmat(k, 266)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
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
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 317) = lmat(k, 317)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = lmat(k, 325)
         mat(k, 326) = lmat(k, 326)
         mat(k, 327) = lmat(k, 327)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 335) = lmat(k, 335)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 357) = mat(k, 357) + lmat(k, 357)
         mat(k, 358) = lmat(k, 358)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 370) = lmat(k, 370)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 376) = lmat(k, 376)
         mat(k, 379) = lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 394) = lmat(k, 394)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 402) = lmat(k, 402)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 405) = lmat(k, 405)
         mat(k, 406) = lmat(k, 406)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 411) = lmat(k, 411)
         mat(k, 412) = lmat(k, 412)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 417) = lmat(k, 417)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 419) = lmat(k, 419)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 422) = lmat(k, 422)
         mat(k, 424) = lmat(k, 424)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 426) = lmat(k, 426)
         mat(k, 427) = mat(k, 427) + lmat(k, 427)
         mat(k, 429) = lmat(k, 429)
         mat(k, 430) = lmat(k, 430)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 432) = lmat(k, 432)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 443) = lmat(k, 443)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 460) = lmat(k, 460)
         mat(k, 461) = lmat(k, 461)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 468) = mat(k, 468) + lmat(k, 468)
         mat(k, 474) = mat(k, 474) + lmat(k, 474)
         mat(k, 476) = lmat(k, 476)
         mat(k, 477) = lmat(k, 477)
         mat(k, 478) = lmat(k, 478)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 480) = lmat(k, 480)
         mat(k, 483) = mat(k, 483) + lmat(k, 483)
         mat(k, 491) = mat(k, 491) + lmat(k, 491)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 499) = mat(k, 499) + lmat(k, 499)
         mat(k, 506) = mat(k, 506) + lmat(k, 506)
         mat(k, 508) = lmat(k, 508)
         mat(k, 509) = lmat(k, 509)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 516) = lmat(k, 516)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 525) = lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 530) = lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 534) = lmat(k, 534)
         mat(k, 536) = lmat(k, 536)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 538) = lmat(k, 538)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 541) = lmat(k, 541)
         mat(k, 542) = lmat(k, 542)
         mat(k, 543) = lmat(k, 543)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 561) = lmat(k, 561)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 572) = mat(k, 572) + lmat(k, 572)
         mat(k, 575) = lmat(k, 575)
         mat(k, 580) = lmat(k, 580)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 585) = lmat(k, 585)
         mat(k, 590) = lmat(k, 590)
         mat(k, 591) = lmat(k, 591)
         mat(k, 592) = lmat(k, 592)
         mat(k, 593) = lmat(k, 593)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 595) = mat(k, 595) + lmat(k, 595)
         mat(k, 600) = lmat(k, 600)
         mat(k, 603) = mat(k, 603) + lmat(k, 603)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 611) = mat(k, 611) + lmat(k, 611)
         mat(k, 614) = mat(k, 614) + lmat(k, 614)
         mat(k, 615) = lmat(k, 615)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 619) = lmat(k, 619)
         mat(k, 621) = lmat(k, 621)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 631) = lmat(k, 631)
         mat(k, 632) = mat(k, 632) + lmat(k, 632)
         mat(k, 636) = lmat(k, 636)
         mat(k, 637) = lmat(k, 637)
         mat(k, 639) = lmat(k, 639)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 641) = lmat(k, 641)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 643) = lmat(k, 643)
         mat(k, 644) = lmat(k, 644)
         mat(k, 645) = lmat(k, 645)
         mat(k, 646) = lmat(k, 646)
         mat(k, 647) = lmat(k, 647)
         mat(k, 649) = lmat(k, 649)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 660) = lmat(k, 660)
         mat(k, 662) = lmat(k, 662)
         mat(k, 663) = mat(k, 663) + lmat(k, 663)
         mat(k, 664) = lmat(k, 664)
         mat(k, 665) = lmat(k, 665)
         mat(k, 666) = mat(k, 666) + lmat(k, 666)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 674) = lmat(k, 674)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 684) = mat(k, 684) + lmat(k, 684)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 692) = lmat(k, 692)
         mat(k, 693) = lmat(k, 693)
         mat(k, 694) = lmat(k, 694)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 696) = lmat(k, 696)
         mat(k, 697) = lmat(k, 697)
         mat(k, 698) = lmat(k, 698)
         mat(k, 699) = lmat(k, 699)
         mat(k, 700) = lmat(k, 700)
         mat(k, 701) = mat(k, 701) + lmat(k, 701)
         mat(k, 706) = lmat(k, 706)
         mat(k, 708) = lmat(k, 708)
         mat(k, 710) = mat(k, 710) + lmat(k, 710)
         mat(k, 711) = lmat(k, 711)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 732) = mat(k, 732) + lmat(k, 732)
         mat(k, 748) = mat(k, 748) + lmat(k, 748)
         mat(k, 759) = mat(k, 759) + lmat(k, 759)
         mat(k, 768) = mat(k, 768) + lmat(k, 768)
         mat(k, 778) = mat(k, 778) + lmat(k, 778)
         mat(k, 779) = lmat(k, 779)
         mat(k, 780) = mat(k, 780) + lmat(k, 780)
         mat(k, 786) = mat(k, 786) + lmat(k, 786)
         mat(k, 791) = mat(k, 791) + lmat(k, 791)
         mat(k, 801) = lmat(k, 801)
         mat(k, 802) = lmat(k, 802)
         mat(k, 803) = lmat(k, 803)
         mat(k, 804) = mat(k, 804) + lmat(k, 804)
         mat(k, 805) = mat(k, 805) + lmat(k, 805)
         mat(k, 806) = mat(k, 806) + lmat(k, 806)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 824) = mat(k, 824) + lmat(k, 824)
         mat(k, 825) = mat(k, 825) + lmat(k, 825)
         mat(k, 827) = mat(k, 827) + lmat(k, 827)
         mat(k, 832) = mat(k, 832) + lmat(k, 832)
         mat(k, 833) = lmat(k, 833)
         mat(k, 836) = lmat(k, 836)
         mat(k, 841) = mat(k, 841) + lmat(k, 841)
         mat(k, 857) = mat(k, 857) + lmat(k, 857)
         mat(k, 859) = lmat(k, 859)
         mat(k, 860) = lmat(k, 860)
         mat(k, 862) = mat(k, 862) + lmat(k, 862)
         mat(k, 864) = mat(k, 864) + lmat(k, 864)
         mat(k, 873) = mat(k, 873) + lmat(k, 873)
         mat(k, 883) = mat(k, 883) + lmat(k, 883)
         mat(k, 884) = mat(k, 884) + lmat(k, 884)
         mat(k, 885) = mat(k, 885) + lmat(k, 885)
         mat(k, 887) = mat(k, 887) + lmat(k, 887)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 891) = lmat(k, 891)
         mat(k, 899) = mat(k, 899) + lmat(k, 899)
         mat(k, 912) = mat(k, 912) + lmat(k, 912)
         mat(k, 922) = mat(k, 922) + lmat(k, 922)
         mat(k, 924) = lmat(k, 924)
         mat(k, 926) = lmat(k, 926)
         mat(k, 929) = mat(k, 929) + lmat(k, 929)
         mat(k, 948) = mat(k, 948) + lmat(k, 948)
         mat(k, 969) = mat(k, 969) + lmat(k, 969)
         mat(k, 971) = lmat(k, 971)
         mat(k, 972) = lmat(k, 972)
         mat(k, 976) = lmat(k, 976)
         mat(k, 977) = lmat(k, 977)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 980) = mat(k, 980) + lmat(k, 980)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 999) = mat(k, 999) + lmat(k, 999)
         mat(k,1024) = mat(k,1024) + lmat(k,1024)
         mat(k,1045) = mat(k,1045) + lmat(k,1045)
         mat(k,1057) = mat(k,1057) + lmat(k,1057)
         mat(k,1058) = mat(k,1058) + lmat(k,1058)
         mat(k,1059) = mat(k,1059) + lmat(k,1059)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1061) = mat(k,1061) + lmat(k,1061)
         mat(k,1062) = mat(k,1062) + lmat(k,1062)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1065) = mat(k,1065) + lmat(k,1065)
         mat(k,1069) = lmat(k,1069)
         mat(k,1073) = mat(k,1073) + lmat(k,1073)
         mat(k,1077) = lmat(k,1077)
         mat(k,1081) = mat(k,1081) + lmat(k,1081)
         mat(k,1083) = lmat(k,1083)
         mat(k,1085) = mat(k,1085) + lmat(k,1085)
         mat(k,1086) = lmat(k,1086)
         mat(k,1091) = lmat(k,1091)
         mat(k,1092) = lmat(k,1092)
         mat(k,1096) = mat(k,1096) + lmat(k,1096)
         mat(k,1097) = lmat(k,1097)
         mat(k,1098) = mat(k,1098) + lmat(k,1098)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1129) = mat(k,1129) + lmat(k,1129)
         mat(k,1140) = mat(k,1140) + lmat(k,1140)
         mat(k,1142) = lmat(k,1142)
         mat(k,1143) = lmat(k,1143)
         mat(k,1145) = mat(k,1145) + lmat(k,1145)
         mat(k,1147) = lmat(k,1147)
         mat(k,1148) = lmat(k,1148)
         mat(k,1149) = lmat(k,1149)
         mat(k,1150) = lmat(k,1150)
         mat(k,1152) = lmat(k,1152)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1155) = lmat(k,1155)
         mat(k,1156) = lmat(k,1156)
         mat(k,1157) = lmat(k,1157)
         mat(k,1161) = mat(k,1161) + lmat(k,1161)
         mat(k,1163) = lmat(k,1163)
         mat(k,1173) = mat(k,1173) + lmat(k,1173)
         mat(k,1193) = mat(k,1193) + lmat(k,1193)
         mat(k,1208) = mat(k,1208) + lmat(k,1208)
         mat(k,1209) = mat(k,1209) + lmat(k,1209)
         mat(k,1212) = mat(k,1212) + lmat(k,1212)
         mat(k,1213) = mat(k,1213) + lmat(k,1213)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1216) = mat(k,1216) + lmat(k,1216)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1221) = mat(k,1221) + lmat(k,1221)
         mat(k,1222) = mat(k,1222) + lmat(k,1222)
         mat(k,1226) = lmat(k,1226)
         mat(k,1230) = lmat(k,1230)
         mat(k,1231) = mat(k,1231) + lmat(k,1231)
         mat(k,1232) = mat(k,1232) + lmat(k,1232)
         mat(k,1243) = lmat(k,1243)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1272) = lmat(k,1272)
         mat(k,1288) = mat(k,1288) + lmat(k,1288)
         mat(k,1298) = mat(k,1298) + lmat(k,1298)
         mat(k,1311) = mat(k,1311) + lmat(k,1311)
         mat(k,1326) = lmat(k,1326)
         mat(k,1328) = mat(k,1328) + lmat(k,1328)
         mat(k,1332) = mat(k,1332) + lmat(k,1332)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1342) = lmat(k,1342)
         mat(k,1354) = mat(k,1354) + lmat(k,1354)
         mat(k,1386) = mat(k,1386) + lmat(k,1386)
         mat(k,1401) = mat(k,1401) + lmat(k,1401)
         mat(k,1415) = mat(k,1415) + lmat(k,1415)
         mat(k,1426) = lmat(k,1426)
         mat(k,1428) = lmat(k,1428)
         mat(k,1429) = mat(k,1429) + lmat(k,1429)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1432) = mat(k,1432) + lmat(k,1432)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1435) = mat(k,1435) + lmat(k,1435)
         mat(k,1439) = mat(k,1439) + lmat(k,1439)
         mat(k,1441) = lmat(k,1441)
         mat(k,1442) = mat(k,1442) + lmat(k,1442)
         mat(k,1447) = mat(k,1447) + lmat(k,1447)
         mat(k,1451) = mat(k,1451) + lmat(k,1451)
         mat(k,1457) = lmat(k,1457)
         mat(k,1463) = mat(k,1463) + lmat(k,1463)
         mat(k,1468) = mat(k,1468) + lmat(k,1468)
         mat(k,1479) = mat(k,1479) + lmat(k,1479)
         mat(k,1480) = lmat(k,1480)
         mat(k,1484) = mat(k,1484) + lmat(k,1484)
         mat(k,1485) = mat(k,1485) + lmat(k,1485)
         mat(k,1526) = mat(k,1526) + lmat(k,1526)
         mat(k,1538) = lmat(k,1538)
         mat(k,1691) = mat(k,1691) + lmat(k,1691)
         mat(k,1744) = mat(k,1744) + lmat(k,1744)
         mat(k,1749) = mat(k,1749) + lmat(k,1749)
         mat(k,1750) = mat(k,1750) + lmat(k,1750)
         mat(k,1757) = mat(k,1757) + lmat(k,1757)
         mat(k,1759) = mat(k,1759) + lmat(k,1759)
         mat(k,1785) = mat(k,1785) + lmat(k,1785)
         mat(k,1842) = mat(k,1842) + lmat(k,1842)
         mat(k,1851) = mat(k,1851) + lmat(k,1851)
         mat(k,1950) = mat(k,1950) + lmat(k,1950)
         mat(k,1959) = mat(k,1959) + lmat(k,1959)
         mat(k,1977) = mat(k,1977) + lmat(k,1977)
         mat(k,1978) = mat(k,1978) + lmat(k,1978)
         mat(k,1984) = mat(k,1984) + lmat(k,1984)
         mat(k,2017) = mat(k,2017) + lmat(k,2017)
         mat(k,2070) = mat(k,2070) + lmat(k,2070)
         mat(k,2124) = mat(k,2124) + lmat(k,2124)
         mat(k,2132) = mat(k,2132) + lmat(k,2132)
         mat(k,2136) = mat(k,2136) + lmat(k,2136)
         mat(k,2145) = mat(k,2145) + lmat(k,2145)
         mat(k,2153) = mat(k,2153) + lmat(k,2153)
         mat(k,2156) = mat(k,2156) + lmat(k,2156)
         mat(k,2188) = mat(k,2188) + lmat(k,2188)
         mat(k,2192) = mat(k,2192) + lmat(k,2192)
         mat(k,2194) = mat(k,2194) + lmat(k,2194)
         mat(k,2201) = mat(k,2201) + lmat(k,2201)
         mat(k,2203) = mat(k,2203) + lmat(k,2203)
         mat(k,2211) = mat(k,2211) + lmat(k,2211)
         mat(k,2226) = mat(k,2226) + lmat(k,2226)
         mat(k,2227) = mat(k,2227) + lmat(k,2227)
         mat(k,2254) = mat(k,2254) + lmat(k,2254)
         mat(k,2258) = mat(k,2258) + lmat(k,2258)
         mat(k,2266) = lmat(k,2266)
         mat(k,2270) = lmat(k,2270)
         mat(k,2272) = mat(k,2272) + lmat(k,2272)
         mat(k,2273) = mat(k,2273) + lmat(k,2273)
         mat(k,2284) = lmat(k,2284)
         mat(k,2285) = mat(k,2285) + lmat(k,2285)
         mat(k, 212) = 0._r8
         mat(k, 213) = 0._r8
         mat(k, 249) = 0._r8
         mat(k, 305) = 0._r8
         mat(k, 341) = 0._r8
         mat(k, 436) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 456) = 0._r8
         mat(k, 484) = 0._r8
         mat(k, 486) = 0._r8
         mat(k, 493) = 0._r8
         mat(k, 502) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 633) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 638) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 658) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 667) = 0._r8
         mat(k, 668) = 0._r8
         mat(k, 671) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 704) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 707) = 0._r8
         mat(k, 709) = 0._r8
         mat(k, 715) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 747) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 776) = 0._r8
         mat(k, 781) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 930) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 956) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 983) = 0._r8
         mat(k, 993) = 0._r8
         mat(k, 998) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1016) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1029) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1084) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1118) = 0._r8
         mat(k,1120) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1192) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1227) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1278) = 0._r8
         mat(k,1279) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1329) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1351) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1397) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1402) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1407) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1416) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1461) = 0._r8
         mat(k,1462) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1467) = 0._r8
         mat(k,1469) = 0._r8
         mat(k,1470) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1486) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1493) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1536) = 0._r8
         mat(k,1605) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1640) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1652) = 0._r8
         mat(k,1674) = 0._r8
         mat(k,1690) = 0._r8
         mat(k,1710) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1722) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1730) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1752) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1758) = 0._r8
         mat(k,1760) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1836) = 0._r8
         mat(k,1837) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1848) = 0._r8
         mat(k,1852) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1885) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1925) = 0._r8
         mat(k,1929) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1937) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1943) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2000) = 0._r8
         mat(k,2001) = 0._r8
         mat(k,2002) = 0._r8
         mat(k,2004) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2011) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2021) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2023) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2059) = 0._r8
         mat(k,2060) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2064) = 0._r8
         mat(k,2065) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2075) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2089) = 0._r8
         mat(k,2093) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2101) = 0._r8
         mat(k,2104) = 0._r8
         mat(k,2106) = 0._r8
         mat(k,2108) = 0._r8
         mat(k,2109) = 0._r8
         mat(k,2113) = 0._r8
         mat(k,2114) = 0._r8
         mat(k,2115) = 0._r8
         mat(k,2117) = 0._r8
         mat(k,2121) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2143) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2149) = 0._r8
         mat(k,2150) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2154) = 0._r8
         mat(k,2155) = 0._r8
         mat(k,2157) = 0._r8
         mat(k,2158) = 0._r8
         mat(k,2171) = 0._r8
         mat(k,2174) = 0._r8
         mat(k,2175) = 0._r8
         mat(k,2179) = 0._r8
         mat(k,2180) = 0._r8
         mat(k,2181) = 0._r8
         mat(k,2182) = 0._r8
         mat(k,2186) = 0._r8
         mat(k,2187) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2190) = 0._r8
         mat(k,2191) = 0._r8
         mat(k,2197) = 0._r8
         mat(k,2198) = 0._r8
         mat(k,2200) = 0._r8
         mat(k,2204) = 0._r8
         mat(k,2212) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2214) = 0._r8
         mat(k,2215) = 0._r8
         mat(k,2217) = 0._r8
         mat(k,2222) = 0._r8
         mat(k,2223) = 0._r8
         mat(k,2224) = 0._r8
         mat(k,2228) = 0._r8
         mat(k,2230) = 0._r8
         mat(k,2236) = 0._r8
         mat(k,2242) = 0._r8
         mat(k,2243) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2253) = 0._r8
         mat(k,2259) = 0._r8
         mat(k,2263) = 0._r8
         mat(k,2265) = 0._r8
         mat(k,2267) = 0._r8
         mat(k,2268) = 0._r8
         mat(k,2269) = 0._r8
         mat(k,2271) = 0._r8
         mat(k,2274) = 0._r8
         mat(k,2275) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2277) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2279) = 0._r8
         mat(k,2280) = 0._r8
         mat(k,2281) = 0._r8
         mat(k,2282) = 0._r8
         mat(k,2283) = 0._r8
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
         mat(k, 129) = mat(k, 129) - dti(k)
         mat(k, 132) = mat(k, 132) - dti(k)
         mat(k, 135) = mat(k, 135) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 141) = mat(k, 141) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 151) = mat(k, 151) - dti(k)
         mat(k, 156) = mat(k, 156) - dti(k)
         mat(k, 161) = mat(k, 161) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 187) = mat(k, 187) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 203) = mat(k, 203) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 226) = mat(k, 226) - dti(k)
         mat(k, 230) = mat(k, 230) - dti(k)
         mat(k, 235) = mat(k, 235) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 241) = mat(k, 241) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 261) = mat(k, 261) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 297) = mat(k, 297) - dti(k)
         mat(k, 304) = mat(k, 304) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 315) = mat(k, 315) - dti(k)
         mat(k, 320) = mat(k, 320) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 344) = mat(k, 344) - dti(k)
         mat(k, 349) = mat(k, 349) - dti(k)
         mat(k, 357) = mat(k, 357) - dti(k)
         mat(k, 365) = mat(k, 365) - dti(k)
         mat(k, 373) = mat(k, 373) - dti(k)
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
         mat(k, 447) = mat(k, 447) - dti(k)
         mat(k, 454) = mat(k, 454) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 463) = mat(k, 463) - dti(k)
         mat(k, 467) = mat(k, 467) - dti(k)
         mat(k, 474) = mat(k, 474) - dti(k)
         mat(k, 483) = mat(k, 483) - dti(k)
         mat(k, 491) = mat(k, 491) - dti(k)
         mat(k, 499) = mat(k, 499) - dti(k)
         mat(k, 506) = mat(k, 506) - dti(k)
         mat(k, 511) = mat(k, 511) - dti(k)
         mat(k, 518) = mat(k, 518) - dti(k)
         mat(k, 524) = mat(k, 524) - dti(k)
         mat(k, 532) = mat(k, 532) - dti(k)
         mat(k, 540) = mat(k, 540) - dti(k)
         mat(k, 548) = mat(k, 548) - dti(k)
         mat(k, 556) = mat(k, 556) - dti(k)
         mat(k, 564) = mat(k, 564) - dti(k)
         mat(k, 572) = mat(k, 572) - dti(k)
         mat(k, 581) = mat(k, 581) - dti(k)
         mat(k, 590) = mat(k, 590) - dti(k)
         mat(k, 594) = mat(k, 594) - dti(k)
         mat(k, 603) = mat(k, 603) - dti(k)
         mat(k, 610) = mat(k, 610) - dti(k)
         mat(k, 617) = mat(k, 617) - dti(k)
         mat(k, 625) = mat(k, 625) - dti(k)
         mat(k, 632) = mat(k, 632) - dti(k)
         mat(k, 642) = mat(k, 642) - dti(k)
         mat(k, 655) = mat(k, 655) - dti(k)
         mat(k, 666) = mat(k, 666) - dti(k)
         mat(k, 677) = mat(k, 677) - dti(k)
         mat(k, 684) = mat(k, 684) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 701) = mat(k, 701) - dti(k)
         mat(k, 714) = mat(k, 714) - dti(k)
         mat(k, 721) = mat(k, 721) - dti(k)
         mat(k, 732) = mat(k, 732) - dti(k)
         mat(k, 748) = mat(k, 748) - dti(k)
         mat(k, 759) = mat(k, 759) - dti(k)
         mat(k, 768) = mat(k, 768) - dti(k)
         mat(k, 778) = mat(k, 778) - dti(k)
         mat(k, 786) = mat(k, 786) - dti(k)
         mat(k, 791) = mat(k, 791) - dti(k)
         mat(k, 801) = mat(k, 801) - dti(k)
         mat(k, 804) = mat(k, 804) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 824) = mat(k, 824) - dti(k)
         mat(k, 832) = mat(k, 832) - dti(k)
         mat(k, 841) = mat(k, 841) - dti(k)
         mat(k, 857) = mat(k, 857) - dti(k)
         mat(k, 864) = mat(k, 864) - dti(k)
         mat(k, 873) = mat(k, 873) - dti(k)
         mat(k, 884) = mat(k, 884) - dti(k)
         mat(k, 899) = mat(k, 899) - dti(k)
         mat(k, 912) = mat(k, 912) - dti(k)
         mat(k, 922) = mat(k, 922) - dti(k)
         mat(k, 929) = mat(k, 929) - dti(k)
         mat(k, 948) = mat(k, 948) - dti(k)
         mat(k, 969) = mat(k, 969) - dti(k)
         mat(k, 979) = mat(k, 979) - dti(k)
         mat(k, 999) = mat(k, 999) - dti(k)
         mat(k,1024) = mat(k,1024) - dti(k)
         mat(k,1045) = mat(k,1045) - dti(k)
         mat(k,1059) = mat(k,1059) - dti(k)
         mat(k,1073) = mat(k,1073) - dti(k)
         mat(k,1085) = mat(k,1085) - dti(k)
         mat(k,1096) = mat(k,1096) - dti(k)
         mat(k,1103) = mat(k,1103) - dti(k)
         mat(k,1115) = mat(k,1115) - dti(k)
         mat(k,1129) = mat(k,1129) - dti(k)
         mat(k,1140) = mat(k,1140) - dti(k)
         mat(k,1153) = mat(k,1153) - dti(k)
         mat(k,1173) = mat(k,1173) - dti(k)
         mat(k,1193) = mat(k,1193) - dti(k)
         mat(k,1209) = mat(k,1209) - dti(k)
         mat(k,1221) = mat(k,1221) - dti(k)
         mat(k,1232) = mat(k,1232) - dti(k)
         mat(k,1256) = mat(k,1256) - dti(k)
         mat(k,1288) = mat(k,1288) - dti(k)
         mat(k,1311) = mat(k,1311) - dti(k)
         mat(k,1332) = mat(k,1332) - dti(k)
         mat(k,1354) = mat(k,1354) - dti(k)
         mat(k,1386) = mat(k,1386) - dti(k)
         mat(k,1401) = mat(k,1401) - dti(k)
         mat(k,1415) = mat(k,1415) - dti(k)
         mat(k,1430) = mat(k,1430) - dti(k)
         mat(k,1447) = mat(k,1447) - dti(k)
         mat(k,1463) = mat(k,1463) - dti(k)
         mat(k,1485) = mat(k,1485) - dti(k)
         mat(k,1526) = mat(k,1526) - dti(k)
         mat(k,1691) = mat(k,1691) - dti(k)
         mat(k,1749) = mat(k,1749) - dti(k)
         mat(k,1842) = mat(k,1842) - dti(k)
         mat(k,1950) = mat(k,1950) - dti(k)
         mat(k,1977) = mat(k,1977) - dti(k)
         mat(k,2017) = mat(k,2017) - dti(k)
         mat(k,2070) = mat(k,2070) - dti(k)
         mat(k,2132) = mat(k,2132) - dti(k)
         mat(k,2156) = mat(k,2156) - dti(k)
         mat(k,2201) = mat(k,2201) - dti(k)
         mat(k,2226) = mat(k,2226) - dti(k)
         mat(k,2258) = mat(k,2258) - dti(k)
         mat(k,2285) = mat(k,2285) - dti(k)
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
