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
         mat(k,561) = -(rxt(k,356)*y(k,217))
         mat(k,1450) = -rxt(k,356)*y(k,1)
         mat(k,1595) = rxt(k,359)*y(k,189)
         mat(k,834) = rxt(k,359)*y(k,124)
         mat(k,537) = -(rxt(k,360)*y(k,217))
         mat(k,1448) = -rxt(k,360)*y(k,2)
         mat(k,833) = rxt(k,357)*y(k,203)
         mat(k,1791) = rxt(k,357)*y(k,189)
         mat(k,788) = -(rxt(k,439)*y(k,126) + rxt(k,440)*y(k,134) + rxt(k,441) &
                      *y(k,217))
         mat(k,1950) = -rxt(k,439)*y(k,6)
         mat(k,2013) = -rxt(k,440)*y(k,6)
         mat(k,1472) = -rxt(k,441)*y(k,6)
         mat(k,113) = -(rxt(k,398)*y(k,217))
         mat(k,1385) = -rxt(k,398)*y(k,7)
         mat(k,304) = -(rxt(k,401)*y(k,217))
         mat(k,1416) = -rxt(k,401)*y(k,8)
         mat(k,400) = rxt(k,399)*y(k,203)
         mat(k,1768) = rxt(k,399)*y(k,191)
         mat(k,114) = .120_r8*rxt(k,398)*y(k,217)
         mat(k,1386) = .120_r8*rxt(k,398)*y(k,7)
         mat(k,785) = .100_r8*rxt(k,440)*y(k,134)
         mat(k,812) = .100_r8*rxt(k,443)*y(k,134)
         mat(k,2002) = .100_r8*rxt(k,440)*y(k,6) + .100_r8*rxt(k,443)*y(k,110)
         mat(k,1583) = .500_r8*rxt(k,400)*y(k,191) + .200_r8*rxt(k,427)*y(k,223) &
                      + .060_r8*rxt(k,433)*y(k,226)
         mat(k,401) = .500_r8*rxt(k,400)*y(k,124)
         mat(k,618) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,634) = .060_r8*rxt(k,433)*y(k,124)
         mat(k,1576) = .200_r8*rxt(k,427)*y(k,223) + .200_r8*rxt(k,433)*y(k,226)
         mat(k,617) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,632) = .200_r8*rxt(k,433)*y(k,124)
         mat(k,1591) = .200_r8*rxt(k,427)*y(k,223) + .150_r8*rxt(k,433)*y(k,226)
         mat(k,619) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,635) = .150_r8*rxt(k,433)*y(k,124)
         mat(k,1577) = .210_r8*rxt(k,433)*y(k,226)
         mat(k,633) = .210_r8*rxt(k,433)*y(k,124)
         mat(k,184) = -(rxt(k,361)*y(k,217))
         mat(k,1398) = -rxt(k,361)*y(k,15)
         mat(k,784) = .050_r8*rxt(k,440)*y(k,134)
         mat(k,811) = .050_r8*rxt(k,443)*y(k,134)
         mat(k,2001) = .050_r8*rxt(k,440)*y(k,6) + .050_r8*rxt(k,443)*y(k,110)
         mat(k,284) = -(rxt(k,327)*y(k,126) + rxt(k,328)*y(k,217))
         mat(k,1944) = -rxt(k,327)*y(k,16)
         mat(k,1413) = -rxt(k,328)*y(k,16)
         mat(k,1278) = -(rxt(k,210)*y(k,42) + rxt(k,211)*y(k,203) + rxt(k,212) &
                      *y(k,134))
         mat(k,1730) = -rxt(k,210)*y(k,17)
         mat(k,1834) = -rxt(k,211)*y(k,17)
         mat(k,2038) = -rxt(k,212)*y(k,17)
         mat(k,1707) = 4.000_r8*rxt(k,213)*y(k,19) + (rxt(k,214)+rxt(k,215))*y(k,59) &
                      + rxt(k,218)*y(k,124) + rxt(k,221)*y(k,133) + rxt(k,468) &
                      *y(k,150) + rxt(k,222)*y(k,217)
         mat(k,1324) = (rxt(k,214)+rxt(k,215))*y(k,19)
         mat(k,705) = rxt(k,223)*y(k,133) + rxt(k,229)*y(k,216) + rxt(k,224)*y(k,217)
         mat(k,1633) = rxt(k,218)*y(k,19)
         mat(k,1864) = rxt(k,221)*y(k,19) + rxt(k,223)*y(k,81)
         mat(k,1112) = rxt(k,468)*y(k,19)
         mat(k,1348) = rxt(k,229)*y(k,81)
         mat(k,1502) = rxt(k,222)*y(k,19) + rxt(k,224)*y(k,81)
         mat(k,1701) = rxt(k,216)*y(k,59)
         mat(k,1318) = rxt(k,216)*y(k,19)
         mat(k,1883) = (rxt(k,530)+rxt(k,535))*y(k,91)
         mat(k,667) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,1716) = -(4._r8*rxt(k,213)*y(k,19) + (rxt(k,214) + rxt(k,215) + rxt(k,216) &
                      ) * y(k,59) + rxt(k,217)*y(k,203) + rxt(k,218)*y(k,124) &
                      + rxt(k,219)*y(k,125) + rxt(k,221)*y(k,133) + rxt(k,222) &
                      *y(k,217) + rxt(k,468)*y(k,150))
         mat(k,1333) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,19)
         mat(k,1843) = -rxt(k,217)*y(k,19)
         mat(k,1642) = -rxt(k,218)*y(k,19)
         mat(k,1552) = -rxt(k,219)*y(k,19)
         mat(k,1873) = -rxt(k,221)*y(k,19)
         mat(k,1511) = -rxt(k,222)*y(k,19)
         mat(k,1118) = -rxt(k,468)*y(k,19)
         mat(k,1282) = rxt(k,212)*y(k,134)
         mat(k,454) = rxt(k,220)*y(k,133)
         mat(k,709) = rxt(k,230)*y(k,216)
         mat(k,671) = rxt(k,225)*y(k,133)
         mat(k,1873) = mat(k,1873) + rxt(k,220)*y(k,20) + rxt(k,225)*y(k,91)
         mat(k,2047) = rxt(k,212)*y(k,17)
         mat(k,1357) = rxt(k,230)*y(k,81)
         mat(k,449) = -(rxt(k,220)*y(k,133))
         mat(k,1854) = -rxt(k,220)*y(k,20)
         mat(k,1703) = rxt(k,219)*y(k,125)
         mat(k,1527) = rxt(k,219)*y(k,19)
         mat(k,190) = -(rxt(k,402)*y(k,217))
         mat(k,1399) = -rxt(k,402)*y(k,22)
         mat(k,1574) = rxt(k,405)*y(k,193)
         mat(k,358) = rxt(k,405)*y(k,124)
         mat(k,263) = -(rxt(k,404)*y(k,217))
         mat(k,1409) = -rxt(k,404)*y(k,23)
         mat(k,359) = rxt(k,403)*y(k,203)
         mat(k,1766) = rxt(k,403)*y(k,193)
         mat(k,225) = -(rxt(k,276)*y(k,56) + rxt(k,277)*y(k,217))
         mat(k,1906) = -rxt(k,276)*y(k,24)
         mat(k,1404) = -rxt(k,277)*y(k,24)
         mat(k,465) = -(rxt(k,278)*y(k,56) + rxt(k,279)*y(k,134) + rxt(k,304)*y(k,217))
         mat(k,1908) = -rxt(k,278)*y(k,25)
         mat(k,2005) = -rxt(k,279)*y(k,25)
         mat(k,1438) = -rxt(k,304)*y(k,25)
         mat(k,198) = -(rxt(k,284)*y(k,217))
         mat(k,1401) = -rxt(k,284)*y(k,26)
         mat(k,712) = .800_r8*rxt(k,280)*y(k,194) + .200_r8*rxt(k,281)*y(k,198)
         mat(k,1651) = .200_r8*rxt(k,281)*y(k,194)
         mat(k,268) = -(rxt(k,285)*y(k,217))
         mat(k,1410) = -rxt(k,285)*y(k,27)
         mat(k,713) = rxt(k,282)*y(k,203)
         mat(k,1767) = rxt(k,282)*y(k,194)
         mat(k,231) = -(rxt(k,286)*y(k,56) + rxt(k,287)*y(k,217))
         mat(k,1907) = -rxt(k,286)*y(k,28)
         mat(k,1405) = -rxt(k,287)*y(k,28)
         mat(k,869) = -(rxt(k,307)*y(k,126) + rxt(k,308)*y(k,134) + rxt(k,325) &
                      *y(k,217))
         mat(k,1954) = -rxt(k,307)*y(k,29)
         mat(k,2017) = -rxt(k,308)*y(k,29)
         mat(k,1477) = -rxt(k,325)*y(k,29)
         mat(k,736) = .130_r8*rxt(k,385)*y(k,134)
         mat(k,2017) = mat(k,2017) + .130_r8*rxt(k,385)*y(k,98)
         mat(k,346) = -(rxt(k,312)*y(k,217))
         mat(k,1422) = -rxt(k,312)*y(k,30)
         mat(k,685) = rxt(k,310)*y(k,203)
         mat(k,1774) = rxt(k,310)*y(k,195)
         mat(k,93) = -(rxt(k,313)*y(k,217))
         mat(k,1382) = -rxt(k,313)*y(k,31)
         mat(k,202) = -(rxt(k,408)*y(k,217))
         mat(k,1402) = -rxt(k,408)*y(k,32)
         mat(k,528) = rxt(k,406)*y(k,203)
         mat(k,1761) = rxt(k,406)*y(k,196)
         mat(k,1739) = -(rxt(k,174)*y(k,56) + rxt(k,210)*y(k,17) + rxt(k,254)*y(k,203) &
                      + rxt(k,255)*y(k,126) + rxt(k,256)*y(k,133) + rxt(k,257) &
                      *y(k,217))
         mat(k,1931) = -rxt(k,174)*y(k,42)
         mat(k,1283) = -rxt(k,210)*y(k,42)
         mat(k,1844) = -rxt(k,254)*y(k,42)
         mat(k,1988) = -rxt(k,255)*y(k,42)
         mat(k,1874) = -rxt(k,256)*y(k,42)
         mat(k,1512) = -rxt(k,257)*y(k,42)
         mat(k,569) = .400_r8*rxt(k,356)*y(k,217)
         mat(k,800) = .340_r8*rxt(k,440)*y(k,134)
         mat(k,290) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,470) = rxt(k,279)*y(k,134)
         mat(k,879) = .500_r8*rxt(k,308)*y(k,134)
         mat(k,439) = .500_r8*rxt(k,296)*y(k,217)
         mat(k,698) = rxt(k,262)*y(k,217)
         mat(k,326) = .300_r8*rxt(k,263)*y(k,217)
         mat(k,1334) = rxt(k,181)*y(k,198)
         mat(k,911) = .800_r8*rxt(k,301)*y(k,217)
         mat(k,746) = .910_r8*rxt(k,385)*y(k,134)
         mat(k,520) = .300_r8*rxt(k,376)*y(k,217)
         mat(k,1083) = .800_r8*rxt(k,380)*y(k,198)
         mat(k,1095) = .120_r8*rxt(k,338)*y(k,134)
         mat(k,490) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,827) = .340_r8*rxt(k,443)*y(k,134)
         mat(k,1173) = .600_r8*rxt(k,352)*y(k,134)
         mat(k,1643) = .100_r8*rxt(k,358)*y(k,189) + rxt(k,261)*y(k,198) &
                      + .500_r8*rxt(k,329)*y(k,200) + .500_r8*rxt(k,298)*y(k,202) &
                      + .920_r8*rxt(k,368)*y(k,205) + .250_r8*rxt(k,336)*y(k,209) &
                      + rxt(k,345)*y(k,211) + rxt(k,319)*y(k,219) + rxt(k,323) &
                      *y(k,220) + .340_r8*rxt(k,452)*y(k,221) + .320_r8*rxt(k,457) &
                      *y(k,222) + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1988) = mat(k,1988) + .500_r8*rxt(k,327)*y(k,16) + rxt(k,369)*y(k,205) &
                      + .250_r8*rxt(k,335)*y(k,209) + rxt(k,346)*y(k,211)
         mat(k,2048) = .340_r8*rxt(k,440)*y(k,6) + rxt(k,279)*y(k,25) &
                      + .500_r8*rxt(k,308)*y(k,29) + .910_r8*rxt(k,385)*y(k,98) &
                      + .120_r8*rxt(k,338)*y(k,105) + .340_r8*rxt(k,443)*y(k,110) &
                      + .600_r8*rxt(k,352)*y(k,111)
         mat(k,387) = rxt(k,303)*y(k,217)
         mat(k,936) = .680_r8*rxt(k,461)*y(k,217)
         mat(k,845) = .100_r8*rxt(k,358)*y(k,124)
         mat(k,721) = .700_r8*rxt(k,281)*y(k,198)
         mat(k,693) = rxt(k,309)*y(k,198)
         mat(k,1271) = rxt(k,292)*y(k,198) + rxt(k,365)*y(k,205) + .250_r8*rxt(k,332) &
                      *y(k,209) + rxt(k,341)*y(k,211) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,1693) = rxt(k,181)*y(k,59) + .800_r8*rxt(k,380)*y(k,101) + rxt(k,261) &
                      *y(k,124) + .700_r8*rxt(k,281)*y(k,194) + rxt(k,309)*y(k,195) &
                      + rxt(k,292)*y(k,197) + (4.000_r8*rxt(k,258)+2.000_r8*rxt(k,259)) &
                      *y(k,198) + 1.500_r8*rxt(k,366)*y(k,205) + .750_r8*rxt(k,371) &
                      *y(k,206) + .880_r8*rxt(k,333)*y(k,209) + 2.000_r8*rxt(k,342) &
                      *y(k,211) + .750_r8*rxt(k,445)*y(k,215) + .800_r8*rxt(k,321) &
                      *y(k,220) + .930_r8*rxt(k,450)*y(k,221) + .950_r8*rxt(k,455) &
                      *y(k,222) + .800_r8*rxt(k,391)*y(k,225)
         mat(k,483) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,609) = .500_r8*rxt(k,298)*y(k,124)
         mat(k,1844) = mat(k,1844) + .450_r8*rxt(k,343)*y(k,211) + .150_r8*rxt(k,322) &
                      *y(k,220)
         mat(k,1223) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) + rxt(k,365) &
                      *y(k,197) + 1.500_r8*rxt(k,366)*y(k,198)
         mat(k,1153) = .750_r8*rxt(k,371)*y(k,198)
         mat(k,1196) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,197) + .880_r8*rxt(k,333)*y(k,198)
         mat(k,1241) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,341)*y(k,197) &
                      + 2.000_r8*rxt(k,342)*y(k,198) + .450_r8*rxt(k,343)*y(k,203) &
                      + 4.000_r8*rxt(k,344)*y(k,211)
         mat(k,1006) = .750_r8*rxt(k,445)*y(k,198)
         mat(k,1512) = mat(k,1512) + .400_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,296) &
                      *y(k,51) + rxt(k,262)*y(k,52) + .300_r8*rxt(k,263)*y(k,53) &
                      + .800_r8*rxt(k,301)*y(k,74) + .300_r8*rxt(k,376)*y(k,99) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,139) &
                      + .680_r8*rxt(k,461)*y(k,178)
         mat(k,664) = rxt(k,319)*y(k,124)
         mat(k,1020) = rxt(k,323)*y(k,124) + .800_r8*rxt(k,321)*y(k,198) &
                      + .150_r8*rxt(k,322)*y(k,203)
         mat(k,987) = .340_r8*rxt(k,452)*y(k,124) + .930_r8*rxt(k,450)*y(k,198)
         mat(k,967) = .320_r8*rxt(k,457)*y(k,124) + .950_r8*rxt(k,455)*y(k,198)
         mat(k,1060) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,390)*y(k,197) &
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
         mat(k,1024) = -(rxt(k,288)*y(k,126) + rxt(k,289)*y(k,217))
         mat(k,1966) = -rxt(k,288)*y(k,45)
         mat(k,1489) = -rxt(k,289)*y(k,45)
         mat(k,565) = .800_r8*rxt(k,356)*y(k,217)
         mat(k,287) = rxt(k,327)*y(k,126)
         mat(k,199) = rxt(k,284)*y(k,217)
         mat(k,270) = .500_r8*rxt(k,285)*y(k,217)
         mat(k,872) = .500_r8*rxt(k,308)*y(k,134)
         mat(k,1162) = .100_r8*rxt(k,352)*y(k,134)
         mat(k,1622) = .400_r8*rxt(k,358)*y(k,189) + rxt(k,283)*y(k,194) &
                      + .270_r8*rxt(k,311)*y(k,195) + rxt(k,329)*y(k,200) + rxt(k,348) &
                      *y(k,213) + rxt(k,319)*y(k,219)
         mat(k,1966) = mat(k,1966) + rxt(k,327)*y(k,16)
         mat(k,2027) = .500_r8*rxt(k,308)*y(k,29) + .100_r8*rxt(k,352)*y(k,111)
         mat(k,839) = .400_r8*rxt(k,358)*y(k,124)
         mat(k,716) = rxt(k,283)*y(k,124) + 3.200_r8*rxt(k,280)*y(k,194) &
                      + .800_r8*rxt(k,281)*y(k,198)
         mat(k,688) = .270_r8*rxt(k,311)*y(k,124)
         mat(k,1673) = .800_r8*rxt(k,281)*y(k,194)
         mat(k,480) = rxt(k,329)*y(k,124)
         mat(k,1822) = .200_r8*rxt(k,347)*y(k,213)
         mat(k,573) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,203)
         mat(k,1489) = mat(k,1489) + .800_r8*rxt(k,356)*y(k,1) + rxt(k,284)*y(k,26) &
                      + .500_r8*rxt(k,285)*y(k,27)
         mat(k,660) = rxt(k,319)*y(k,124)
         mat(k,87) = -(rxt(k,290)*y(k,217))
         mat(k,1381) = -rxt(k,290)*y(k,47)
         mat(k,847) = -(rxt(k,326)*y(k,217))
         mat(k,1475) = -rxt(k,326)*y(k,48)
         mat(k,564) = .800_r8*rxt(k,356)*y(k,217)
         mat(k,790) = .520_r8*rxt(k,440)*y(k,134)
         mat(k,286) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,817) = .520_r8*rxt(k,443)*y(k,134)
         mat(k,1610) = .250_r8*rxt(k,358)*y(k,189) + .820_r8*rxt(k,311)*y(k,195) &
                      + .500_r8*rxt(k,329)*y(k,200) + .270_r8*rxt(k,452)*y(k,221) &
                      + .040_r8*rxt(k,457)*y(k,222)
         mat(k,1953) = .500_r8*rxt(k,327)*y(k,16)
         mat(k,2016) = .520_r8*rxt(k,440)*y(k,6) + .520_r8*rxt(k,443)*y(k,110)
         mat(k,929) = .500_r8*rxt(k,461)*y(k,217)
         mat(k,838) = .250_r8*rxt(k,358)*y(k,124)
         mat(k,687) = .820_r8*rxt(k,311)*y(k,124) + .820_r8*rxt(k,309)*y(k,198)
         mat(k,1662) = .820_r8*rxt(k,309)*y(k,195) + .150_r8*rxt(k,450)*y(k,221) &
                      + .025_r8*rxt(k,455)*y(k,222)
         mat(k,478) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,1475) = mat(k,1475) + .800_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,461) &
                      *y(k,178)
         mat(k,974) = .270_r8*rxt(k,452)*y(k,124) + .150_r8*rxt(k,450)*y(k,198)
         mat(k,952) = .040_r8*rxt(k,457)*y(k,124) + .025_r8*rxt(k,455)*y(k,198)
         mat(k,1100) = -(rxt(k,314)*y(k,126) + rxt(k,315)*y(k,217))
         mat(k,1970) = -rxt(k,314)*y(k,49)
         mat(k,1494) = -rxt(k,315)*y(k,49)
         mat(k,944) = rxt(k,316)*y(k,217)
         mat(k,1089) = .880_r8*rxt(k,338)*y(k,134)
         mat(k,1163) = .500_r8*rxt(k,352)*y(k,134)
         mat(k,1626) = .170_r8*rxt(k,411)*y(k,199) + .050_r8*rxt(k,374)*y(k,206) &
                      + .250_r8*rxt(k,336)*y(k,209) + .170_r8*rxt(k,417)*y(k,212) &
                      + .400_r8*rxt(k,427)*y(k,223) + .250_r8*rxt(k,393)*y(k,225) &
                      + .540_r8*rxt(k,433)*y(k,226) + .510_r8*rxt(k,436)*y(k,228)
         mat(k,1970) = mat(k,1970) + .050_r8*rxt(k,375)*y(k,206) + .250_r8*rxt(k,335) &
                      *y(k,209) + .250_r8*rxt(k,394)*y(k,225)
         mat(k,752) = rxt(k,317)*y(k,217)
         mat(k,2030) = .880_r8*rxt(k,338)*y(k,105) + .500_r8*rxt(k,352)*y(k,111)
         mat(k,1258) = .250_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,1677) = .240_r8*rxt(k,333)*y(k,209) + .500_r8*rxt(k,321)*y(k,220) &
                      + .100_r8*rxt(k,391)*y(k,225)
         mat(k,651) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,203)
         mat(k,1827) = .070_r8*rxt(k,410)*y(k,199) + .070_r8*rxt(k,416)*y(k,212)
         mat(k,1141) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1186) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,197) + .240_r8*rxt(k,333)*y(k,198)
         mat(k,772) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,203)
         mat(k,1494) = mat(k,1494) + rxt(k,316)*y(k,95) + rxt(k,317)*y(k,127)
         mat(k,1014) = .500_r8*rxt(k,321)*y(k,198)
         mat(k,627) = .400_r8*rxt(k,427)*y(k,124)
         mat(k,1053) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
         mat(k,643) = .540_r8*rxt(k,433)*y(k,124)
         mat(k,412) = .510_r8*rxt(k,436)*y(k,124)
         mat(k,473) = -(rxt(k,295)*y(k,217))
         mat(k,1439) = -rxt(k,295)*y(k,50)
         mat(k,865) = .120_r8*rxt(k,308)*y(k,134)
         mat(k,2006) = .120_r8*rxt(k,308)*y(k,29)
         mat(k,1249) = .100_r8*rxt(k,292)*y(k,198) + .150_r8*rxt(k,293)*y(k,203)
         mat(k,1655) = .100_r8*rxt(k,292)*y(k,197)
         mat(k,1787) = .150_r8*rxt(k,293)*y(k,197) + .150_r8*rxt(k,343)*y(k,211)
         mat(k,1229) = .150_r8*rxt(k,343)*y(k,203)
         mat(k,435) = -(rxt(k,296)*y(k,217))
         mat(k,1435) = -rxt(k,296)*y(k,51)
         mat(k,1248) = .400_r8*rxt(k,293)*y(k,203)
         mat(k,1785) = .400_r8*rxt(k,293)*y(k,197) + .400_r8*rxt(k,343)*y(k,211)
         mat(k,1228) = .400_r8*rxt(k,343)*y(k,203)
         mat(k,696) = -(rxt(k,262)*y(k,217))
         mat(k,1462) = -rxt(k,262)*y(k,52)
         mat(k,1066) = .200_r8*rxt(k,380)*y(k,198)
         mat(k,714) = .300_r8*rxt(k,281)*y(k,198)
         mat(k,1658) = .200_r8*rxt(k,380)*y(k,101) + .300_r8*rxt(k,281)*y(k,194) &
                      + 2.000_r8*rxt(k,259)*y(k,198) + .250_r8*rxt(k,366)*y(k,205) &
                      + .250_r8*rxt(k,371)*y(k,206) + .250_r8*rxt(k,333)*y(k,209) &
                      + .250_r8*rxt(k,445)*y(k,215) + .500_r8*rxt(k,321)*y(k,220) &
                      + .250_r8*rxt(k,450)*y(k,221) + .250_r8*rxt(k,455)*y(k,222) &
                      + .300_r8*rxt(k,391)*y(k,225)
         mat(k,1202) = .250_r8*rxt(k,366)*y(k,198)
         mat(k,1129) = .250_r8*rxt(k,371)*y(k,198)
         mat(k,1179) = .250_r8*rxt(k,333)*y(k,198)
         mat(k,992) = .250_r8*rxt(k,445)*y(k,198)
         mat(k,1011) = .500_r8*rxt(k,321)*y(k,198)
         mat(k,973) = .250_r8*rxt(k,450)*y(k,198)
         mat(k,951) = .250_r8*rxt(k,455)*y(k,198)
         mat(k,1047) = .300_r8*rxt(k,391)*y(k,198)
         mat(k,322) = -(rxt(k,263)*y(k,217))
         mat(k,1418) = -rxt(k,263)*y(k,53)
         mat(k,1654) = rxt(k,260)*y(k,203)
         mat(k,1770) = rxt(k,260)*y(k,198)
         mat(k,1935) = -(rxt(k,174)*y(k,42) + rxt(k,176)*y(k,77) + rxt(k,177)*y(k,79) &
                      + (rxt(k,178) + rxt(k,179)) * y(k,203) + rxt(k,180)*y(k,134) &
                      + rxt(k,187)*y(k,60) + rxt(k,196)*y(k,92) + rxt(k,286)*y(k,28))
         mat(k,1743) = -rxt(k,174)*y(k,56)
         mat(k,1043) = -rxt(k,176)*y(k,56)
         mat(k,502) = -rxt(k,177)*y(k,56)
         mat(k,1848) = -(rxt(k,178) + rxt(k,179)) * y(k,56)
         mat(k,2052) = -rxt(k,180)*y(k,56)
         mat(k,862) = -rxt(k,187)*y(k,56)
         mat(k,729) = -rxt(k,196)*y(k,56)
         mat(k,235) = -rxt(k,286)*y(k,56)
         mat(k,1721) = rxt(k,215)*y(k,59)
         mat(k,1338) = rxt(k,215)*y(k,19) + (4.000_r8*rxt(k,182)+2.000_r8*rxt(k,184)) &
                      *y(k,59) + rxt(k,186)*y(k,124) + rxt(k,191)*y(k,133) &
                      + rxt(k,469)*y(k,150) + rxt(k,181)*y(k,198) + rxt(k,192) &
                      *y(k,217)
         mat(k,137) = rxt(k,236)*y(k,216)
         mat(k,1901) = rxt(k,194)*y(k,133) + rxt(k,206)*y(k,216) + rxt(k,195)*y(k,217)
         mat(k,1647) = rxt(k,186)*y(k,59)
         mat(k,1878) = rxt(k,191)*y(k,59) + rxt(k,194)*y(k,85)
         mat(k,1121) = rxt(k,469)*y(k,59)
         mat(k,1697) = rxt(k,181)*y(k,59)
         mat(k,1362) = rxt(k,236)*y(k,65) + rxt(k,206)*y(k,85)
         mat(k,1516) = rxt(k,192)*y(k,59) + rxt(k,195)*y(k,85)
         mat(k,1905) = rxt(k,187)*y(k,60)
         mat(k,1317) = 2.000_r8*rxt(k,183)*y(k,59)
         mat(k,853) = rxt(k,187)*y(k,56) + (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,1882) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60) + (rxt(k,523) &
                       +rxt(k,529)+rxt(k,534))*y(k,92)
         mat(k,723) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85)
         mat(k,1316) = 2.000_r8*rxt(k,208)*y(k,59)
         mat(k,1327) = -(rxt(k,181)*y(k,198) + (4._r8*rxt(k,182) + 4._r8*rxt(k,183) &
                      + 4._r8*rxt(k,184) + 4._r8*rxt(k,208)) * y(k,59) + rxt(k,185) &
                      *y(k,203) + rxt(k,186)*y(k,124) + rxt(k,188)*y(k,125) + rxt(k,191) &
                      *y(k,133) + (rxt(k,192) + rxt(k,193)) * y(k,217) + (rxt(k,214) &
                      + rxt(k,215) + rxt(k,216)) * y(k,19) + rxt(k,469)*y(k,150))
         mat(k,1686) = -rxt(k,181)*y(k,59)
         mat(k,1837) = -rxt(k,185)*y(k,59)
         mat(k,1636) = -rxt(k,186)*y(k,59)
         mat(k,1546) = -rxt(k,188)*y(k,59)
         mat(k,1867) = -rxt(k,191)*y(k,59)
         mat(k,1505) = -(rxt(k,192) + rxt(k,193)) * y(k,59)
         mat(k,1710) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,59)
         mat(k,1114) = -rxt(k,469)*y(k,59)
         mat(k,1924) = rxt(k,196)*y(k,92) + rxt(k,180)*y(k,134) + rxt(k,179)*y(k,203)
         mat(k,857) = rxt(k,189)*y(k,133)
         mat(k,1890) = rxt(k,207)*y(k,216)
         mat(k,725) = rxt(k,196)*y(k,56) + rxt(k,197)*y(k,133) + rxt(k,198)*y(k,217)
         mat(k,1867) = mat(k,1867) + rxt(k,189)*y(k,60) + rxt(k,197)*y(k,92)
         mat(k,2041) = rxt(k,180)*y(k,56)
         mat(k,255) = rxt(k,474)*y(k,150)
         mat(k,1114) = mat(k,1114) + rxt(k,474)*y(k,136)
         mat(k,1837) = mat(k,1837) + rxt(k,179)*y(k,56)
         mat(k,1351) = rxt(k,207)*y(k,85)
         mat(k,1505) = mat(k,1505) + rxt(k,198)*y(k,92)
         mat(k,855) = -(rxt(k,187)*y(k,56) + rxt(k,189)*y(k,133) + rxt(k,190)*y(k,217) &
                      + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85))
         mat(k,1915) = -rxt(k,187)*y(k,60)
         mat(k,1860) = -rxt(k,189)*y(k,60)
         mat(k,1476) = -rxt(k,190)*y(k,60)
         mat(k,1886) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60)
         mat(k,1322) = rxt(k,188)*y(k,125)
         mat(k,1536) = rxt(k,188)*y(k,59)
         mat(k,939) = -((rxt(k,265) + rxt(k,275)) * y(k,217))
         mat(k,1483) = -(rxt(k,265) + rxt(k,275)) * y(k,62)
         mat(k,793) = .230_r8*rxt(k,440)*y(k,134)
         mat(k,1277) = rxt(k,210)*y(k,42)
         mat(k,228) = .350_r8*rxt(k,277)*y(k,217)
         mat(k,468) = .630_r8*rxt(k,279)*y(k,134)
         mat(k,870) = .560_r8*rxt(k,308)*y(k,134)
         mat(k,1728) = rxt(k,210)*y(k,17) + rxt(k,174)*y(k,56) + rxt(k,255)*y(k,126) &
                      + rxt(k,256)*y(k,133) + rxt(k,257)*y(k,217)
         mat(k,1099) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,217)
         mat(k,1917) = rxt(k,174)*y(k,42)
         mat(k,766) = rxt(k,302)*y(k,217)
         mat(k,737) = .620_r8*rxt(k,385)*y(k,134)
         mat(k,1087) = .650_r8*rxt(k,338)*y(k,134)
         mat(k,820) = .230_r8*rxt(k,443)*y(k,134)
         mat(k,1160) = .560_r8*rxt(k,352)*y(k,134)
         mat(k,1616) = .170_r8*rxt(k,411)*y(k,199) + .220_r8*rxt(k,336)*y(k,209) &
                      + .400_r8*rxt(k,414)*y(k,210) + .350_r8*rxt(k,417)*y(k,212) &
                      + .225_r8*rxt(k,452)*y(k,221) + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1960) = rxt(k,255)*y(k,42) + rxt(k,314)*y(k,49) + .220_r8*rxt(k,335) &
                      *y(k,209) + .500_r8*rxt(k,394)*y(k,225)
         mat(k,1861) = rxt(k,256)*y(k,42) + rxt(k,464)*y(k,137)
         mat(k,2021) = .230_r8*rxt(k,440)*y(k,6) + .630_r8*rxt(k,279)*y(k,25) &
                      + .560_r8*rxt(k,308)*y(k,29) + .620_r8*rxt(k,385)*y(k,98) &
                      + .650_r8*rxt(k,338)*y(k,105) + .230_r8*rxt(k,443)*y(k,110) &
                      + .560_r8*rxt(k,352)*y(k,111)
         mat(k,279) = rxt(k,464)*y(k,133) + rxt(k,465)*y(k,217)
         mat(k,931) = .700_r8*rxt(k,461)*y(k,217)
         mat(k,1253) = .220_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,1667) = .110_r8*rxt(k,333)*y(k,209) + .125_r8*rxt(k,450)*y(k,221) &
                      + .200_r8*rxt(k,391)*y(k,225)
         mat(k,650) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,203)
         mat(k,1816) = .070_r8*rxt(k,410)*y(k,199) + .160_r8*rxt(k,413)*y(k,210) &
                      + .140_r8*rxt(k,416)*y(k,212)
         mat(k,1182) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,197) + .110_r8*rxt(k,333)*y(k,198)
         mat(k,613) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,203)
         mat(k,771) = .350_r8*rxt(k,417)*y(k,124) + .140_r8*rxt(k,416)*y(k,203)
         mat(k,1483) = mat(k,1483) + .350_r8*rxt(k,277)*y(k,24) + rxt(k,257)*y(k,42) &
                      + rxt(k,315)*y(k,49) + rxt(k,302)*y(k,75) + rxt(k,465)*y(k,137) &
                      + .700_r8*rxt(k,461)*y(k,178)
         mat(k,977) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,198)
         mat(k,1050) = .250_r8*rxt(k,393)*y(k,124) + .500_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .200_r8*rxt(k,391)*y(k,198)
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
         mat(k,109) = -(rxt(k,235)*y(k,216))
         mat(k,1342) = -rxt(k,235)*y(k,64)
         mat(k,134) = -(rxt(k,236)*y(k,216))
         mat(k,1344) = -rxt(k,236)*y(k,65)
         mat(k,146) = -(rxt(k,409)*y(k,217))
         mat(k,1391) = -rxt(k,409)*y(k,66)
         mat(k,140) = .180_r8*rxt(k,429)*y(k,217)
         mat(k,1391) = mat(k,1391) + .180_r8*rxt(k,429)*y(k,180)
         mat(k,210) = -(rxt(k,462)*y(k,126) + (rxt(k,463) + rxt(k,476)) * y(k,217))
         mat(k,1941) = -rxt(k,462)*y(k,67)
         mat(k,1403) = -(rxt(k,463) + rxt(k,476)) * y(k,67)
         mat(k,602) = rxt(k,297)*y(k,203)
         mat(k,1759) = rxt(k,297)*y(k,202)
         mat(k,677) = -(rxt(k,232)*y(k,77) + rxt(k,233)*y(k,229) + rxt(k,234)*y(k,89))
         mat(k,1034) = -rxt(k,232)*y(k,73)
         mat(k,2059) = -rxt(k,233)*y(k,73)
         mat(k,1289) = -rxt(k,234)*y(k,73)
         mat(k,110) = 2.000_r8*rxt(k,235)*y(k,216)
         mat(k,135) = rxt(k,236)*y(k,216)
         mat(k,1345) = 2.000_r8*rxt(k,235)*y(k,64) + rxt(k,236)*y(k,65)
         mat(k,908) = -(rxt(k,301)*y(k,217))
         mat(k,1480) = -rxt(k,301)*y(k,74)
         mat(k,514) = .700_r8*rxt(k,376)*y(k,217)
         mat(k,459) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,300) = rxt(k,388)*y(k,217)
         mat(k,1613) = .050_r8*rxt(k,374)*y(k,206) + .530_r8*rxt(k,336)*y(k,209) &
                      + .225_r8*rxt(k,452)*y(k,221) + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1957) = .050_r8*rxt(k,375)*y(k,206) + .530_r8*rxt(k,335)*y(k,209) &
                      + .250_r8*rxt(k,394)*y(k,225)
         mat(k,1252) = .530_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,1665) = .260_r8*rxt(k,333)*y(k,209) + .125_r8*rxt(k,450)*y(k,221) &
                      + .100_r8*rxt(k,391)*y(k,225)
         mat(k,1133) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1180) = .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335)*y(k,126) &
                      + .530_r8*rxt(k,332)*y(k,197) + .260_r8*rxt(k,333)*y(k,198)
         mat(k,1480) = mat(k,1480) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + rxt(k,388)*y(k,115)
         mat(k,975) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,198)
         mat(k,1049) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
         mat(k,765) = -(rxt(k,302)*y(k,217))
         mat(k,1470) = -rxt(k,302)*y(k,75)
         mat(k,227) = .650_r8*rxt(k,277)*y(k,217)
         mat(k,907) = .200_r8*rxt(k,301)*y(k,217)
         mat(k,894) = rxt(k,389)*y(k,217)
         mat(k,1607) = rxt(k,400)*y(k,191) + .050_r8*rxt(k,374)*y(k,206) &
                      + .400_r8*rxt(k,414)*y(k,210) + .170_r8*rxt(k,417)*y(k,212) &
                      + .700_r8*rxt(k,420)*y(k,218) + .600_r8*rxt(k,427)*y(k,223) &
                      + .250_r8*rxt(k,393)*y(k,225) + .340_r8*rxt(k,433)*y(k,226) &
                      + .170_r8*rxt(k,436)*y(k,228)
         mat(k,1949) = .050_r8*rxt(k,375)*y(k,206) + .250_r8*rxt(k,394)*y(k,225)
         mat(k,404) = rxt(k,400)*y(k,124)
         mat(k,1250) = .250_r8*rxt(k,390)*y(k,225)
         mat(k,1661) = .100_r8*rxt(k,391)*y(k,225)
         mat(k,1809) = .160_r8*rxt(k,413)*y(k,210) + .070_r8*rxt(k,416)*y(k,212)
         mat(k,1131) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,612) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,203)
         mat(k,769) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,203)
         mat(k,1470) = mat(k,1470) + .650_r8*rxt(k,277)*y(k,24) + .200_r8*rxt(k,301) &
                      *y(k,74) + rxt(k,389)*y(k,116)
         mat(k,374) = .700_r8*rxt(k,420)*y(k,124)
         mat(k,624) = .600_r8*rxt(k,427)*y(k,124)
         mat(k,1048) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
         mat(k,640) = .340_r8*rxt(k,433)*y(k,124)
         mat(k,411) = .170_r8*rxt(k,436)*y(k,124)
         mat(k,1304) = -((rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,203) + rxt(k,140) &
                      *y(k,134))
         mat(k,1836) = -(rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,76)
         mat(k,2040) = -rxt(k,140)*y(k,76)
         mat(k,1732) = rxt(k,257)*y(k,217)
         mat(k,1923) = rxt(k,176)*y(k,77)
         mat(k,940) = rxt(k,275)*y(k,217)
         mat(k,680) = rxt(k,232)*y(k,77)
         mat(k,1037) = rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73) + rxt(k,132)*y(k,133) &
                      + rxt(k,124)*y(k,216) + rxt(k,141)*y(k,217)
         mat(k,706) = rxt(k,230)*y(k,216)
         mat(k,1889) = rxt(k,207)*y(k,216)
         mat(k,293) = rxt(k,162)*y(k,217)
         mat(k,1866) = rxt(k,132)*y(k,77) + rxt(k,144)*y(k,217)
         mat(k,281) = rxt(k,465)*y(k,217)
         mat(k,419) = rxt(k,470)*y(k,217)
         mat(k,1113) = rxt(k,475)*y(k,217)
         mat(k,1350) = rxt(k,124)*y(k,77) + rxt(k,230)*y(k,81) + rxt(k,207)*y(k,85)
         mat(k,1504) = rxt(k,257)*y(k,42) + rxt(k,275)*y(k,62) + rxt(k,141)*y(k,77) &
                      + rxt(k,162)*y(k,112) + rxt(k,144)*y(k,133) + rxt(k,465) &
                      *y(k,137) + rxt(k,470)*y(k,148) + rxt(k,475)*y(k,150)
         mat(k,1035) = -(rxt(k,124)*y(k,216) + rxt(k,132)*y(k,133) + rxt(k,141) &
                      *y(k,217) + rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73))
         mat(k,1347) = -rxt(k,124)*y(k,77)
         mat(k,1862) = -rxt(k,132)*y(k,77)
         mat(k,1490) = -rxt(k,141)*y(k,77)
         mat(k,1919) = -rxt(k,176)*y(k,77)
         mat(k,678) = -rxt(k,232)*y(k,77)
         mat(k,1302) = rxt(k,134)*y(k,203)
         mat(k,1823) = rxt(k,134)*y(k,76)
         mat(k,497) = -(rxt(k,133)*y(k,133) + rxt(k,142)*y(k,217) + rxt(k,177)*y(k,56))
         mat(k,1855) = -rxt(k,133)*y(k,79)
         mat(k,1443) = -rxt(k,142)*y(k,79)
         mat(k,1909) = -rxt(k,177)*y(k,79)
         mat(k,1788) = 2.000_r8*rxt(k,148)*y(k,203)
         mat(k,1443) = mat(k,1443) + 2.000_r8*rxt(k,147)*y(k,217)
         mat(k,193) = rxt(k,478)*y(k,229)
         mat(k,2056) = rxt(k,478)*y(k,152)
         mat(k,704) = -(rxt(k,223)*y(k,133) + rxt(k,224)*y(k,217) + (rxt(k,229) &
                      + rxt(k,230)) * y(k,216))
         mat(k,1857) = -rxt(k,223)*y(k,81)
         mat(k,1464) = -rxt(k,224)*y(k,81)
         mat(k,1346) = -(rxt(k,229) + rxt(k,230)) * y(k,81)
         mat(k,1276) = rxt(k,210)*y(k,42) + rxt(k,211)*y(k,203)
         mat(k,1727) = rxt(k,210)*y(k,17)
         mat(k,1805) = rxt(k,211)*y(k,17)
         mat(k,1900) = -(rxt(k,194)*y(k,133) + rxt(k,195)*y(k,217) + (rxt(k,206) &
                      + rxt(k,207)) * y(k,216) + (rxt(k,523) + rxt(k,529) + rxt(k,534) &
                      ) * y(k,92) + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60) &
                      + (rxt(k,530) + rxt(k,535)) * y(k,91))
         mat(k,1877) = -rxt(k,194)*y(k,85)
         mat(k,1515) = -rxt(k,195)*y(k,85)
         mat(k,1361) = -(rxt(k,206) + rxt(k,207)) * y(k,85)
         mat(k,728) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85)
         mat(k,861) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85)
         mat(k,673) = -(rxt(k,530) + rxt(k,535)) * y(k,85)
         mat(k,234) = rxt(k,286)*y(k,56)
         mat(k,1742) = rxt(k,174)*y(k,56)
         mat(k,1934) = rxt(k,286)*y(k,28) + rxt(k,174)*y(k,42) + rxt(k,176)*y(k,77) &
                      + rxt(k,177)*y(k,79) + rxt(k,196)*y(k,92) + rxt(k,178)*y(k,203)
         mat(k,1337) = rxt(k,193)*y(k,217)
         mat(k,1042) = rxt(k,176)*y(k,56)
         mat(k,501) = rxt(k,177)*y(k,56)
         mat(k,728) = mat(k,728) + rxt(k,196)*y(k,56)
         mat(k,1847) = rxt(k,178)*y(k,56)
         mat(k,1515) = mat(k,1515) + rxt(k,193)*y(k,59)
         mat(k,130) = -(rxt(k,266)*y(k,217) + rxt(k,274)*y(k,216))
         mat(k,1388) = -rxt(k,266)*y(k,86)
         mat(k,1343) = -rxt(k,274)*y(k,86)
         mat(k,700) = -(rxt(k,267)*y(k,217))
         mat(k,1463) = -rxt(k,267)*y(k,87)
         mat(k,786) = .050_r8*rxt(k,440)*y(k,134)
         mat(k,226) = .350_r8*rxt(k,277)*y(k,217)
         mat(k,467) = .370_r8*rxt(k,279)*y(k,134)
         mat(k,867) = .120_r8*rxt(k,308)*y(k,134)
         mat(k,734) = .110_r8*rxt(k,385)*y(k,134)
         mat(k,1086) = .330_r8*rxt(k,338)*y(k,134)
         mat(k,813) = .050_r8*rxt(k,443)*y(k,134)
         mat(k,1158) = .120_r8*rxt(k,352)*y(k,134)
         mat(k,1604) = rxt(k,270)*y(k,204)
         mat(k,2009) = .050_r8*rxt(k,440)*y(k,6) + .370_r8*rxt(k,279)*y(k,25) &
                      + .120_r8*rxt(k,308)*y(k,29) + .110_r8*rxt(k,385)*y(k,98) &
                      + .330_r8*rxt(k,338)*y(k,105) + .050_r8*rxt(k,443)*y(k,110) &
                      + .120_r8*rxt(k,352)*y(k,111)
         mat(k,1804) = rxt(k,268)*y(k,204)
         mat(k,367) = rxt(k,270)*y(k,124) + rxt(k,268)*y(k,203)
         mat(k,1463) = mat(k,1463) + .350_r8*rxt(k,277)*y(k,24)
         mat(k,676) = rxt(k,232)*y(k,77) + rxt(k,234)*y(k,89) + rxt(k,233)*y(k,229)
         mat(k,1033) = rxt(k,232)*y(k,73)
         mat(k,1288) = rxt(k,234)*y(k,73)
         mat(k,2057) = rxt(k,233)*y(k,73)
         mat(k,1291) = -(rxt(k,171)*y(k,217) + rxt(k,234)*y(k,73))
         mat(k,1503) = -rxt(k,171)*y(k,89)
         mat(k,679) = -rxt(k,234)*y(k,89)
         mat(k,1731) = rxt(k,255)*y(k,126)
         mat(k,1026) = rxt(k,288)*y(k,126)
         mat(k,1102) = rxt(k,314)*y(k,126)
         mat(k,856) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,212) = rxt(k,462)*y(k,126)
         mat(k,1888) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60)
         mat(k,1544) = rxt(k,170)*y(k,217)
         mat(k,1979) = rxt(k,255)*y(k,42) + rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) &
                      + rxt(k,462)*y(k,67)
         mat(k,1503) = mat(k,1503) + rxt(k,170)*y(k,125)
         mat(k,334) = -(rxt(k,149)*y(k,217))
         mat(k,1420) = -rxt(k,149)*y(k,90)
         mat(k,1523) = rxt(k,168)*y(k,203)
         mat(k,1772) = rxt(k,168)*y(k,125)
         mat(k,668) = -(rxt(k,225)*y(k,133) + (rxt(k,530) + rxt(k,535)) * y(k,85))
         mat(k,1856) = -rxt(k,225)*y(k,91)
         mat(k,1884) = -(rxt(k,530) + rxt(k,535)) * y(k,91)
         mat(k,1704) = rxt(k,217)*y(k,203)
         mat(k,1802) = rxt(k,217)*y(k,19)
         mat(k,724) = -(rxt(k,196)*y(k,56) + rxt(k,197)*y(k,133) + rxt(k,198)*y(k,217) &
                      + (rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85))
         mat(k,1913) = -rxt(k,196)*y(k,92)
         mat(k,1858) = -rxt(k,197)*y(k,92)
         mat(k,1466) = -rxt(k,198)*y(k,92)
         mat(k,1885) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,92)
         mat(k,1320) = rxt(k,185)*y(k,203)
         mat(k,854) = rxt(k,190)*y(k,217)
         mat(k,1807) = rxt(k,185)*y(k,59)
         mat(k,1466) = mat(k,1466) + rxt(k,190)*y(k,60)
         mat(k,916) = -(rxt(k,331)*y(k,217))
         mat(k,1481) = -rxt(k,331)*y(k,93)
         mat(k,515) = .300_r8*rxt(k,376)*y(k,217)
         mat(k,460) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,1614) = rxt(k,330)*y(k,200) + rxt(k,337)*y(k,209)
         mat(k,479) = rxt(k,330)*y(k,124)
         mat(k,1181) = rxt(k,337)*y(k,124)
         mat(k,1481) = mat(k,1481) + .300_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100)
         mat(k,176) = -(rxt(k,362)*y(k,217))
         mat(k,1396) = -rxt(k,362)*y(k,94)
         mat(k,943) = -(rxt(k,316)*y(k,217))
         mat(k,1484) = -rxt(k,316)*y(k,95)
         mat(k,516) = .700_r8*rxt(k,376)*y(k,217)
         mat(k,461) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,486) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,1617) = .050_r8*rxt(k,374)*y(k,206) + .220_r8*rxt(k,336)*y(k,209) &
                      + .250_r8*rxt(k,393)*y(k,225)
         mat(k,1961) = .050_r8*rxt(k,375)*y(k,206) + .220_r8*rxt(k,335)*y(k,209) &
                      + .250_r8*rxt(k,394)*y(k,225)
         mat(k,443) = .500_r8*rxt(k,320)*y(k,217)
         mat(k,1254) = .220_r8*rxt(k,332)*y(k,209) + .250_r8*rxt(k,390)*y(k,225)
         mat(k,1668) = .230_r8*rxt(k,333)*y(k,209) + .200_r8*rxt(k,321)*y(k,220) &
                      + .100_r8*rxt(k,391)*y(k,225)
         mat(k,1136) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1183) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,197) + .230_r8*rxt(k,333)*y(k,198)
         mat(k,1484) = mat(k,1484) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + .500_r8*rxt(k,351)*y(k,109) + .500_r8*rxt(k,320) &
                      *y(k,146)
         mat(k,1012) = .200_r8*rxt(k,321)*y(k,198)
         mat(k,1051) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,197) + .100_r8*rxt(k,391)*y(k,198)
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
         mat(k,273) = -(rxt(k,363)*y(k,217))
         mat(k,1411) = -rxt(k,363)*y(k,96)
         mat(k,1578) = .870_r8*rxt(k,374)*y(k,206)
         mat(k,1943) = .950_r8*rxt(k,375)*y(k,206)
         mat(k,1246) = rxt(k,370)*y(k,206)
         mat(k,1652) = .750_r8*rxt(k,371)*y(k,206)
         mat(k,1125) = .870_r8*rxt(k,374)*y(k,124) + .950_r8*rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,197) + .750_r8*rxt(k,371)*y(k,198)
         mat(k,103) = -(rxt(k,364)*y(k,217))
         mat(k,1383) = -rxt(k,364)*y(k,97)
         mat(k,579) = .600_r8*rxt(k,387)*y(k,217)
         mat(k,1383) = mat(k,1383) + .600_r8*rxt(k,387)*y(k,103)
         mat(k,735) = -(rxt(k,378)*y(k,126) + rxt(k,385)*y(k,134) + rxt(k,386) &
                      *y(k,217))
         mat(k,1946) = -rxt(k,378)*y(k,98)
         mat(k,2010) = -rxt(k,385)*y(k,98)
         mat(k,1467) = -rxt(k,386)*y(k,98)
         mat(k,513) = -(rxt(k,376)*y(k,217))
         mat(k,1445) = -rxt(k,376)*y(k,99)
         mat(k,1592) = .080_r8*rxt(k,368)*y(k,205)
         mat(k,1200) = .080_r8*rxt(k,368)*y(k,124)
         mat(k,457) = -(rxt(k,377)*y(k,217))
         mat(k,1437) = -rxt(k,377)*y(k,100)
         mat(k,1589) = .080_r8*rxt(k,374)*y(k,206)
         mat(k,1126) = .080_r8*rxt(k,374)*y(k,124)
         mat(k,1072) = -(rxt(k,379)*y(k,197) + rxt(k,380)*y(k,198) + rxt(k,381) &
                      *y(k,203) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126))
         mat(k,1256) = -rxt(k,379)*y(k,101)
         mat(k,1675) = -rxt(k,380)*y(k,101)
         mat(k,1825) = -rxt(k,381)*y(k,101)
         mat(k,1624) = -rxt(k,382)*y(k,101)
         mat(k,1968) = -rxt(k,383)*y(k,101)
         mat(k,738) = rxt(k,378)*y(k,126)
         mat(k,1968) = mat(k,1968) + rxt(k,378)*y(k,98)
         mat(k,328) = -(rxt(k,384)*y(k,217))
         mat(k,1419) = -rxt(k,384)*y(k,102)
         mat(k,1064) = rxt(k,381)*y(k,203)
         mat(k,1771) = rxt(k,381)*y(k,101)
         mat(k,580) = -(rxt(k,387)*y(k,217))
         mat(k,1452) = -rxt(k,387)*y(k,103)
         mat(k,1794) = rxt(k,367)*y(k,205) + rxt(k,372)*y(k,206)
         mat(k,1201) = rxt(k,367)*y(k,203)
         mat(k,1128) = rxt(k,372)*y(k,203)
         mat(k,65) = -(rxt(k,509)*y(k,217))
         mat(k,1376) = -rxt(k,509)*y(k,104)
         mat(k,1088) = -(rxt(k,338)*y(k,134) + rxt(k,339)*y(k,217))
         mat(k,2029) = -rxt(k,338)*y(k,105)
         mat(k,1493) = -rxt(k,339)*y(k,105)
         mat(k,739) = .300_r8*rxt(k,385)*y(k,134)
         mat(k,1625) = .360_r8*rxt(k,368)*y(k,205)
         mat(k,1969) = .400_r8*rxt(k,369)*y(k,205)
         mat(k,2029) = mat(k,2029) + .300_r8*rxt(k,385)*y(k,98)
         mat(k,1257) = .390_r8*rxt(k,365)*y(k,205)
         mat(k,1676) = .310_r8*rxt(k,366)*y(k,205)
         mat(k,1210) = .360_r8*rxt(k,368)*y(k,124) + .400_r8*rxt(k,369)*y(k,126) &
                      + .390_r8*rxt(k,365)*y(k,197) + .310_r8*rxt(k,366)*y(k,198)
         mat(k,237) = -(rxt(k,340)*y(k,217))
         mat(k,1406) = -rxt(k,340)*y(k,106)
         mat(k,1763) = rxt(k,334)*y(k,209)
         mat(k,1178) = rxt(k,334)*y(k,203)
         mat(k,423) = -(rxt(k,349)*y(k,217))
         mat(k,1433) = -rxt(k,349)*y(k,107)
         mat(k,1587) = .800_r8*rxt(k,358)*y(k,189)
         mat(k,832) = .800_r8*rxt(k,358)*y(k,124)
         mat(k,242) = -(rxt(k,350)*y(k,217))
         mat(k,1407) = -rxt(k,350)*y(k,108)
         mat(k,1764) = .800_r8*rxt(k,347)*y(k,213)
         mat(k,571) = .800_r8*rxt(k,347)*y(k,203)
         mat(k,485) = -(rxt(k,351)*y(k,217))
         mat(k,1441) = -rxt(k,351)*y(k,109)
         mat(k,1528) = rxt(k,354)*y(k,211)
         mat(k,1230) = rxt(k,354)*y(k,125)
         mat(k,815) = -(rxt(k,442)*y(k,126) + rxt(k,443)*y(k,134) + rxt(k,444) &
                      *y(k,217))
         mat(k,1951) = -rxt(k,442)*y(k,110)
         mat(k,2014) = -rxt(k,443)*y(k,110)
         mat(k,1473) = -rxt(k,444)*y(k,110)
         mat(k,1164) = -(rxt(k,352)*y(k,134) + rxt(k,353)*y(k,217))
         mat(k,2033) = -rxt(k,352)*y(k,111)
         mat(k,1497) = -rxt(k,353)*y(k,111)
         mat(k,741) = .200_r8*rxt(k,385)*y(k,134)
         mat(k,1628) = .560_r8*rxt(k,368)*y(k,205)
         mat(k,1973) = .600_r8*rxt(k,369)*y(k,205)
         mat(k,2033) = mat(k,2033) + .200_r8*rxt(k,385)*y(k,98)
         mat(k,1260) = .610_r8*rxt(k,365)*y(k,205)
         mat(k,1679) = .440_r8*rxt(k,366)*y(k,205)
         mat(k,1212) = .560_r8*rxt(k,368)*y(k,124) + .600_r8*rxt(k,369)*y(k,126) &
                      + .610_r8*rxt(k,365)*y(k,197) + .440_r8*rxt(k,366)*y(k,198)
         mat(k,292) = -(rxt(k,150)*y(k,124) + (rxt(k,151) + rxt(k,152) + rxt(k,153) &
                      ) * y(k,125) + rxt(k,162)*y(k,217))
         mat(k,1579) = -rxt(k,150)*y(k,112)
         mat(k,1522) = -(rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,112)
         mat(k,1414) = -rxt(k,162)*y(k,112)
         mat(k,1521) = rxt(k,169)*y(k,126)
         mat(k,1942) = rxt(k,169)*y(k,125)
         mat(k,298) = -(rxt(k,388)*y(k,217))
         mat(k,1415) = -rxt(k,388)*y(k,115)
         mat(k,1063) = .200_r8*rxt(k,380)*y(k,198)
         mat(k,1653) = .200_r8*rxt(k,380)*y(k,101)
         mat(k,896) = -(rxt(k,389)*y(k,217))
         mat(k,1479) = -rxt(k,389)*y(k,116)
         mat(k,1068) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) + rxt(k,379)*y(k,197) &
                      + .800_r8*rxt(k,380)*y(k,198)
         mat(k,1612) = rxt(k,382)*y(k,101)
         mat(k,1956) = rxt(k,383)*y(k,101)
         mat(k,1251) = rxt(k,379)*y(k,101)
         mat(k,1664) = .800_r8*rxt(k,380)*y(k,101)
         mat(k,84) = -(rxt(k,479)*y(k,217))
         mat(k,1380) = -rxt(k,479)*y(k,120)
         mat(k,1640) = -(rxt(k,150)*y(k,112) + rxt(k,159)*y(k,126) + rxt(k,163) &
                      *y(k,203) + rxt(k,164)*y(k,134) + rxt(k,165)*y(k,133) + rxt(k,186) &
                      *y(k,59) + rxt(k,218)*y(k,19) + rxt(k,261)*y(k,198) + rxt(k,270) &
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
         mat(k,296) = -rxt(k,150)*y(k,124)
         mat(k,1985) = -rxt(k,159)*y(k,124)
         mat(k,1841) = -rxt(k,163)*y(k,124)
         mat(k,2045) = -rxt(k,164)*y(k,124)
         mat(k,1871) = -rxt(k,165)*y(k,124)
         mat(k,1331) = -rxt(k,186)*y(k,124)
         mat(k,1714) = -rxt(k,218)*y(k,124)
         mat(k,1690) = -rxt(k,261)*y(k,124)
         mat(k,369) = -rxt(k,270)*y(k,124)
         mat(k,719) = -rxt(k,283)*y(k,124)
         mat(k,1269) = -rxt(k,294)*y(k,124)
         mat(k,608) = -rxt(k,298)*y(k,124)
         mat(k,691) = -rxt(k,311)*y(k,124)
         mat(k,663) = -rxt(k,319)*y(k,124)
         mat(k,1018) = -rxt(k,323)*y(k,124)
         mat(k,482) = -(rxt(k,329) + rxt(k,330)) * y(k,124)
         mat(k,1194) = -(rxt(k,336) + rxt(k,337)) * y(k,124)
         mat(k,1239) = -rxt(k,345)*y(k,124)
         mat(k,577) = -rxt(k,348)*y(k,124)
         mat(k,843) = -(rxt(k,358) + rxt(k,359)) * y(k,124)
         mat(k,1221) = -rxt(k,368)*y(k,124)
         mat(k,1151) = -rxt(k,374)*y(k,124)
         mat(k,1081) = -rxt(k,382)*y(k,124)
         mat(k,1058) = -rxt(k,393)*y(k,124)
         mat(k,433) = -rxt(k,397)*y(k,124)
         mat(k,407) = -rxt(k,400)*y(k,124)
         mat(k,364) = -rxt(k,405)*y(k,124)
         mat(k,533) = -rxt(k,407)*y(k,124)
         mat(k,654) = -rxt(k,411)*y(k,124)
         mat(k,615) = -rxt(k,414)*y(k,124)
         mat(k,775) = -rxt(k,417)*y(k,124)
         mat(k,377) = -rxt(k,420)*y(k,124)
         mat(k,630) = -rxt(k,427)*y(k,124)
         mat(k,647) = -rxt(k,433)*y(k,124)
         mat(k,415) = -rxt(k,436)*y(k,124)
         mat(k,1004) = -rxt(k,447)*y(k,124)
         mat(k,985) = -rxt(k,452)*y(k,124)
         mat(k,965) = -rxt(k,457)*y(k,124)
         mat(k,296) = mat(k,296) + 2.000_r8*rxt(k,152)*y(k,125) + rxt(k,162)*y(k,217)
         mat(k,1550) = 2.000_r8*rxt(k,152)*y(k,112) + rxt(k,155)*y(k,133) + rxt(k,471) &
                      *y(k,150)
         mat(k,1871) = mat(k,1871) + rxt(k,155)*y(k,125)
         mat(k,1117) = rxt(k,471)*y(k,125)
         mat(k,1509) = rxt(k,162)*y(k,112)
         mat(k,1549) = -((rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,112) + (rxt(k,155) &
                      + rxt(k,157)) * y(k,133) + rxt(k,156)*y(k,134) + rxt(k,168) &
                      *y(k,203) + rxt(k,169)*y(k,126) + rxt(k,170)*y(k,217) + rxt(k,188) &
                      *y(k,59) + rxt(k,219)*y(k,19) + rxt(k,305)*y(k,197) + rxt(k,354) &
                      *y(k,211) + rxt(k,412)*y(k,199) + rxt(k,415)*y(k,210) + rxt(k,418) &
                      *y(k,212) + rxt(k,422)*y(k,141) + rxt(k,425)*y(k,188) + rxt(k,471) &
                      *y(k,150))
         mat(k,295) = -(rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,125)
         mat(k,1870) = -(rxt(k,155) + rxt(k,157)) * y(k,125)
         mat(k,2044) = -rxt(k,156)*y(k,125)
         mat(k,1840) = -rxt(k,168)*y(k,125)
         mat(k,1984) = -rxt(k,169)*y(k,125)
         mat(k,1508) = -rxt(k,170)*y(k,125)
         mat(k,1330) = -rxt(k,188)*y(k,125)
         mat(k,1713) = -rxt(k,219)*y(k,125)
         mat(k,1268) = -rxt(k,305)*y(k,125)
         mat(k,1238) = -rxt(k,354)*y(k,125)
         mat(k,653) = -rxt(k,412)*y(k,125)
         mat(k,614) = -rxt(k,415)*y(k,125)
         mat(k,774) = -rxt(k,418)*y(k,125)
         mat(k,391) = -rxt(k,422)*y(k,125)
         mat(k,432) = -rxt(k,425)*y(k,125)
         mat(k,1116) = -rxt(k,471)*y(k,125)
         mat(k,568) = rxt(k,356)*y(k,217)
         mat(k,289) = rxt(k,327)*y(k,126)
         mat(k,1713) = mat(k,1713) + rxt(k,218)*y(k,124)
         mat(k,1330) = mat(k,1330) + rxt(k,186)*y(k,124)
         mat(k,336) = rxt(k,149)*y(k,217)
         mat(k,519) = .700_r8*rxt(k,376)*y(k,217)
         mat(k,1080) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126)
         mat(k,1639) = rxt(k,218)*y(k,19) + rxt(k,186)*y(k,59) + rxt(k,382)*y(k,101) &
                      + 2.000_r8*rxt(k,159)*y(k,126) + rxt(k,165)*y(k,133) &
                      + rxt(k,164)*y(k,134) + rxt(k,397)*y(k,188) + rxt(k,358) &
                      *y(k,189) + rxt(k,400)*y(k,191) + rxt(k,405)*y(k,193) &
                      + rxt(k,283)*y(k,194) + rxt(k,311)*y(k,195) + rxt(k,407) &
                      *y(k,196) + rxt(k,294)*y(k,197) + rxt(k,261)*y(k,198) &
                      + rxt(k,411)*y(k,199) + rxt(k,329)*y(k,200) + rxt(k,298) &
                      *y(k,202) + rxt(k,163)*y(k,203) + rxt(k,270)*y(k,204) &
                      + .920_r8*rxt(k,368)*y(k,205) + .920_r8*rxt(k,374)*y(k,206) &
                      + rxt(k,336)*y(k,209) + rxt(k,414)*y(k,210) + rxt(k,345) &
                      *y(k,211) + rxt(k,417)*y(k,212) + rxt(k,348)*y(k,213) &
                      + 1.600_r8*rxt(k,447)*y(k,215) + rxt(k,420)*y(k,218) &
                      + rxt(k,319)*y(k,219) + rxt(k,323)*y(k,220) + .900_r8*rxt(k,452) &
                      *y(k,221) + .800_r8*rxt(k,457)*y(k,222) + rxt(k,427)*y(k,223) &
                      + rxt(k,393)*y(k,225) + rxt(k,433)*y(k,226) + rxt(k,436) &
                      *y(k,228)
         mat(k,1984) = mat(k,1984) + rxt(k,327)*y(k,16) + rxt(k,383)*y(k,101) &
                      + 2.000_r8*rxt(k,159)*y(k,124) + rxt(k,160)*y(k,133) &
                      + rxt(k,158)*y(k,203) + rxt(k,369)*y(k,205) + rxt(k,375) &
                      *y(k,206) + rxt(k,335)*y(k,209) + rxt(k,346)*y(k,211) &
                      + 2.000_r8*rxt(k,448)*y(k,215) + rxt(k,161)*y(k,217) &
                      + rxt(k,394)*y(k,225)
         mat(k,755) = rxt(k,317)*y(k,217)
         mat(k,1870) = mat(k,1870) + rxt(k,165)*y(k,124) + rxt(k,160)*y(k,126)
         mat(k,2044) = mat(k,2044) + rxt(k,164)*y(k,124)
         mat(k,526) = rxt(k,454)*y(k,217)
         mat(k,432) = mat(k,432) + rxt(k,397)*y(k,124)
         mat(k,842) = rxt(k,358)*y(k,124)
         mat(k,406) = rxt(k,400)*y(k,124)
         mat(k,363) = rxt(k,405)*y(k,124)
         mat(k,718) = rxt(k,283)*y(k,124)
         mat(k,690) = rxt(k,311)*y(k,124)
         mat(k,532) = rxt(k,407)*y(k,124)
         mat(k,1268) = mat(k,1268) + rxt(k,294)*y(k,124)
         mat(k,1689) = rxt(k,261)*y(k,124) + .500_r8*rxt(k,445)*y(k,215)
         mat(k,653) = mat(k,653) + rxt(k,411)*y(k,124)
         mat(k,481) = rxt(k,329)*y(k,124)
         mat(k,607) = rxt(k,298)*y(k,124)
         mat(k,1840) = mat(k,1840) + rxt(k,163)*y(k,124) + rxt(k,158)*y(k,126)
         mat(k,368) = rxt(k,270)*y(k,124)
         mat(k,1220) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126)
         mat(k,1150) = .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,1193) = rxt(k,336)*y(k,124) + rxt(k,335)*y(k,126)
         mat(k,614) = mat(k,614) + rxt(k,414)*y(k,124)
         mat(k,1238) = mat(k,1238) + rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126)
         mat(k,774) = mat(k,774) + rxt(k,417)*y(k,124)
         mat(k,576) = rxt(k,348)*y(k,124)
         mat(k,1003) = 1.600_r8*rxt(k,447)*y(k,124) + 2.000_r8*rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,198)
         mat(k,1508) = mat(k,1508) + rxt(k,356)*y(k,1) + rxt(k,149)*y(k,90) &
                      + .700_r8*rxt(k,376)*y(k,99) + rxt(k,161)*y(k,126) + rxt(k,317) &
                      *y(k,127) + rxt(k,454)*y(k,175)
         mat(k,376) = rxt(k,420)*y(k,124)
         mat(k,662) = rxt(k,319)*y(k,124)
         mat(k,1017) = rxt(k,323)*y(k,124)
         mat(k,984) = .900_r8*rxt(k,452)*y(k,124)
         mat(k,964) = .800_r8*rxt(k,457)*y(k,124)
         mat(k,629) = rxt(k,427)*y(k,124)
         mat(k,1057) = rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126)
         mat(k,646) = rxt(k,433)*y(k,124)
         mat(k,414) = rxt(k,436)*y(k,124)
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
         mat(k,1993) = -(rxt(k,158)*y(k,203) + rxt(k,159)*y(k,124) + rxt(k,160) &
                      *y(k,133) + rxt(k,161)*y(k,217) + rxt(k,169)*y(k,125) + rxt(k,255) &
                      *y(k,42) + rxt(k,288)*y(k,45) + rxt(k,307)*y(k,29) + rxt(k,314) &
                      *y(k,49) + rxt(k,327)*y(k,16) + rxt(k,335)*y(k,209) + rxt(k,346) &
                      *y(k,211) + rxt(k,369)*y(k,205) + rxt(k,375)*y(k,206) + rxt(k,378) &
                      *y(k,98) + rxt(k,383)*y(k,101) + rxt(k,394)*y(k,225) + rxt(k,439) &
                      *y(k,6) + rxt(k,442)*y(k,110) + rxt(k,448)*y(k,215) + rxt(k,459) &
                      *y(k,177) + rxt(k,462)*y(k,67))
         mat(k,1849) = -rxt(k,158)*y(k,126)
         mat(k,1648) = -rxt(k,159)*y(k,126)
         mat(k,1879) = -rxt(k,160)*y(k,126)
         mat(k,1517) = -rxt(k,161)*y(k,126)
         mat(k,1558) = -rxt(k,169)*y(k,126)
         mat(k,1744) = -rxt(k,255)*y(k,126)
         mat(k,1031) = -rxt(k,288)*y(k,126)
         mat(k,881) = -rxt(k,307)*y(k,126)
         mat(k,1106) = -rxt(k,314)*y(k,126)
         mat(k,291) = -rxt(k,327)*y(k,126)
         mat(k,1198) = -rxt(k,335)*y(k,126)
         mat(k,1243) = -rxt(k,346)*y(k,126)
         mat(k,1225) = -rxt(k,369)*y(k,126)
         mat(k,1155) = -rxt(k,375)*y(k,126)
         mat(k,748) = -rxt(k,378)*y(k,126)
         mat(k,1085) = -rxt(k,383)*y(k,126)
         mat(k,1062) = -rxt(k,394)*y(k,126)
         mat(k,802) = -rxt(k,439)*y(k,126)
         mat(k,829) = -rxt(k,442)*y(k,126)
         mat(k,1008) = -rxt(k,448)*y(k,126)
         mat(k,891) = -rxt(k,459)*y(k,126)
         mat(k,215) = -rxt(k,462)*y(k,126)
         mat(k,456) = rxt(k,220)*y(k,133)
         mat(k,1936) = rxt(k,187)*y(k,60)
         mat(k,863) = rxt(k,187)*y(k,56) + rxt(k,189)*y(k,133) + rxt(k,190)*y(k,217)
         mat(k,683) = rxt(k,234)*y(k,89)
         mat(k,1300) = rxt(k,234)*y(k,73) + rxt(k,171)*y(k,217)
         mat(k,492) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,1558) = mat(k,1558) + rxt(k,157)*y(k,133) + rxt(k,156)*y(k,134)
         mat(k,1879) = mat(k,1879) + rxt(k,220)*y(k,20) + rxt(k,189)*y(k,60) &
                      + rxt(k,157)*y(k,125)
         mat(k,2053) = rxt(k,156)*y(k,125)
         mat(k,388) = rxt(k,303)*y(k,217)
         mat(k,1517) = mat(k,1517) + rxt(k,190)*y(k,60) + rxt(k,171)*y(k,89) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,139)
         mat(k,751) = -(rxt(k,317)*y(k,217))
         mat(k,1468) = -rxt(k,317)*y(k,127)
         mat(k,868) = rxt(k,307)*y(k,126)
         mat(k,458) = .500_r8*rxt(k,377)*y(k,217)
         mat(k,330) = rxt(k,384)*y(k,217)
         mat(k,299) = rxt(k,388)*y(k,217)
         mat(k,893) = rxt(k,389)*y(k,217)
         mat(k,1947) = rxt(k,307)*y(k,29)
         mat(k,1468) = mat(k,1468) + .500_r8*rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) &
                      + rxt(k,388)*y(k,115) + rxt(k,389)*y(k,116)
         mat(k,316) = -(rxt(k,449)*y(k,217))
         mat(k,1417) = -rxt(k,449)*y(k,128)
         mat(k,1769) = rxt(k,446)*y(k,215)
         mat(k,990) = rxt(k,446)*y(k,203)
         mat(k,1876) = -(rxt(k,129)*y(k,134) + 4._r8*rxt(k,130)*y(k,133) + rxt(k,132) &
                      *y(k,77) + rxt(k,133)*y(k,79) + rxt(k,138)*y(k,203) + rxt(k,144) &
                      *y(k,217) + (rxt(k,155) + rxt(k,157)) * y(k,125) + rxt(k,160) &
                      *y(k,126) + rxt(k,165)*y(k,124) + rxt(k,189)*y(k,60) + rxt(k,191) &
                      *y(k,59) + rxt(k,194)*y(k,85) + rxt(k,197)*y(k,92) + rxt(k,220) &
                      *y(k,20) + rxt(k,221)*y(k,19) + rxt(k,223)*y(k,81) + rxt(k,225) &
                      *y(k,91) + rxt(k,256)*y(k,42) + rxt(k,464)*y(k,137))
         mat(k,2050) = -rxt(k,129)*y(k,133)
         mat(k,1041) = -rxt(k,132)*y(k,133)
         mat(k,500) = -rxt(k,133)*y(k,133)
         mat(k,1846) = -rxt(k,138)*y(k,133)
         mat(k,1514) = -rxt(k,144)*y(k,133)
         mat(k,1555) = -(rxt(k,155) + rxt(k,157)) * y(k,133)
         mat(k,1990) = -rxt(k,160)*y(k,133)
         mat(k,1645) = -rxt(k,165)*y(k,133)
         mat(k,860) = -rxt(k,189)*y(k,133)
         mat(k,1336) = -rxt(k,191)*y(k,133)
         mat(k,1899) = -rxt(k,194)*y(k,133)
         mat(k,727) = -rxt(k,197)*y(k,133)
         mat(k,455) = -rxt(k,220)*y(k,133)
         mat(k,1719) = -rxt(k,221)*y(k,133)
         mat(k,710) = -rxt(k,223)*y(k,133)
         mat(k,672) = -rxt(k,225)*y(k,133)
         mat(k,1741) = -rxt(k,256)*y(k,133)
         mat(k,283) = -rxt(k,464)*y(k,133)
         mat(k,1310) = rxt(k,136)*y(k,203)
         mat(k,297) = rxt(k,150)*y(k,124) + rxt(k,151)*y(k,125)
         mat(k,1645) = mat(k,1645) + rxt(k,150)*y(k,112)
         mat(k,1555) = mat(k,1555) + rxt(k,151)*y(k,112)
         mat(k,1846) = mat(k,1846) + rxt(k,136)*y(k,76)
         mat(k,1514) = mat(k,1514) + 2.000_r8*rxt(k,146)*y(k,217)
         mat(k,2054) = -(rxt(k,128)*y(k,216) + rxt(k,129)*y(k,133) + rxt(k,139) &
                      *y(k,203) + rxt(k,140)*y(k,76) + rxt(k,145)*y(k,217) + rxt(k,156) &
                      *y(k,125) + rxt(k,164)*y(k,124) + rxt(k,180)*y(k,56) + rxt(k,212) &
                      *y(k,17) + rxt(k,279)*y(k,25) + rxt(k,308)*y(k,29) + rxt(k,338) &
                      *y(k,105) + rxt(k,352)*y(k,111) + rxt(k,385)*y(k,98) + rxt(k,423) &
                      *y(k,141) + rxt(k,440)*y(k,6) + rxt(k,443)*y(k,110) + rxt(k,467) &
                      *y(k,148) + rxt(k,473)*y(k,150))
         mat(k,1364) = -rxt(k,128)*y(k,134)
         mat(k,1880) = -rxt(k,129)*y(k,134)
         mat(k,1850) = -rxt(k,139)*y(k,134)
         mat(k,1314) = -rxt(k,140)*y(k,134)
         mat(k,1518) = -rxt(k,145)*y(k,134)
         mat(k,1559) = -rxt(k,156)*y(k,134)
         mat(k,1649) = -rxt(k,164)*y(k,134)
         mat(k,1937) = -rxt(k,180)*y(k,134)
         mat(k,1286) = -rxt(k,212)*y(k,134)
         mat(k,472) = -rxt(k,279)*y(k,134)
         mat(k,882) = -rxt(k,308)*y(k,134)
         mat(k,1097) = -rxt(k,338)*y(k,134)
         mat(k,1176) = -rxt(k,352)*y(k,134)
         mat(k,749) = -rxt(k,385)*y(k,134)
         mat(k,392) = -rxt(k,423)*y(k,134)
         mat(k,803) = -rxt(k,440)*y(k,134)
         mat(k,830) = -rxt(k,443)*y(k,134)
         mat(k,422) = -rxt(k,467)*y(k,134)
         mat(k,1122) = -rxt(k,473)*y(k,134)
         mat(k,1274) = .150_r8*rxt(k,293)*y(k,203)
         mat(k,1850) = mat(k,1850) + .150_r8*rxt(k,293)*y(k,197) + .150_r8*rxt(k,343) &
                      *y(k,211)
         mat(k,1244) = .150_r8*rxt(k,343)*y(k,203)
         mat(k,252) = -(rxt(k,474)*y(k,150))
         mat(k,1108) = -rxt(k,474)*y(k,136)
         mat(k,1702) = rxt(k,214)*y(k,59)
         mat(k,1319) = rxt(k,214)*y(k,19) + 2.000_r8*rxt(k,184)*y(k,59)
         mat(k,276) = -(rxt(k,464)*y(k,133) + rxt(k,465)*y(k,217))
         mat(k,1852) = -rxt(k,464)*y(k,137)
         mat(k,1412) = -rxt(k,465)*y(k,137)
         mat(k,913) = rxt(k,331)*y(k,217)
         mat(k,1573) = .100_r8*rxt(k,452)*y(k,221)
         mat(k,1397) = rxt(k,331)*y(k,93)
         mat(k,971) = .100_r8*rxt(k,452)*y(k,124)
         mat(k,382) = -(rxt(k,303)*y(k,217))
         mat(k,1427) = -rxt(k,303)*y(k,139)
         mat(k,1524) = rxt(k,305)*y(k,197)
         mat(k,1247) = rxt(k,305)*y(k,125)
         mat(k,1520) = rxt(k,425)*y(k,188)
         mat(k,428) = rxt(k,425)*y(k,125)
         mat(k,389) = -(rxt(k,422)*y(k,125) + rxt(k,423)*y(k,134))
         mat(k,1525) = -rxt(k,422)*y(k,141)
         mat(k,2003) = -rxt(k,423)*y(k,141)
         mat(k,148) = .070_r8*rxt(k,409)*y(k,217)
         mat(k,1584) = rxt(k,407)*y(k,196)
         mat(k,125) = .060_r8*rxt(k,421)*y(k,217)
         mat(k,169) = .070_r8*rxt(k,437)*y(k,217)
         mat(k,529) = rxt(k,407)*y(k,124)
         mat(k,1428) = .070_r8*rxt(k,409)*y(k,66) + .060_r8*rxt(k,421)*y(k,142) &
                      + .070_r8*rxt(k,437)*y(k,184)
         mat(k,123) = -(rxt(k,421)*y(k,217))
         mat(k,1387) = -rxt(k,421)*y(k,142)
         mat(k,115) = .530_r8*rxt(k,398)*y(k,217)
         mat(k,1387) = mat(k,1387) + .530_r8*rxt(k,398)*y(k,7)
         mat(k,257) = -(rxt(k,424)*y(k,217))
         mat(k,1408) = -rxt(k,424)*y(k,143)
         mat(k,1765) = rxt(k,419)*y(k,218)
         mat(k,372) = rxt(k,419)*y(k,203)
         mat(k,441) = -(rxt(k,320)*y(k,217))
         mat(k,1436) = -rxt(k,320)*y(k,146)
         mat(k,1786) = rxt(k,318)*y(k,219)
         mat(k,657) = rxt(k,318)*y(k,203)
         mat(k,340) = -(rxt(k,324)*y(k,217))
         mat(k,1421) = -rxt(k,324)*y(k,147)
         mat(k,1773) = .850_r8*rxt(k,322)*y(k,220)
         mat(k,1010) = .850_r8*rxt(k,322)*y(k,203)
         mat(k,417) = -(rxt(k,467)*y(k,134) + rxt(k,470)*y(k,217))
         mat(k,2004) = -rxt(k,467)*y(k,148)
         mat(k,1432) = -rxt(k,470)*y(k,148)
         mat(k,1111) = -(rxt(k,468)*y(k,19) + rxt(k,469)*y(k,59) + rxt(k,471)*y(k,125) &
                      + rxt(k,473)*y(k,134) + rxt(k,474)*y(k,136) + rxt(k,475) &
                      *y(k,217))
         mat(k,1706) = -rxt(k,468)*y(k,150)
         mat(k,1323) = -rxt(k,469)*y(k,150)
         mat(k,1540) = -rxt(k,471)*y(k,150)
         mat(k,2031) = -rxt(k,473)*y(k,150)
         mat(k,254) = -rxt(k,474)*y(k,150)
         mat(k,1495) = -rxt(k,475)*y(k,150)
         mat(k,1863) = rxt(k,464)*y(k,137)
         mat(k,2031) = mat(k,2031) + rxt(k,467)*y(k,148)
         mat(k,280) = rxt(k,464)*y(k,133)
         mat(k,418) = rxt(k,467)*y(k,134) + rxt(k,470)*y(k,217)
         mat(k,1495) = mat(k,1495) + rxt(k,470)*y(k,148)
         mat(k,758) = -(rxt(k,477)*y(k,217))
         mat(k,1469) = -rxt(k,477)*y(k,151)
         mat(k,1705) = rxt(k,468)*y(k,150)
         mat(k,1321) = rxt(k,469)*y(k,150)
         mat(k,211) = rxt(k,462)*y(k,126) + (rxt(k,463)+.500_r8*rxt(k,476))*y(k,217)
         mat(k,1533) = rxt(k,471)*y(k,150)
         mat(k,1948) = rxt(k,462)*y(k,67)
         mat(k,2011) = rxt(k,473)*y(k,150)
         mat(k,253) = rxt(k,474)*y(k,150)
         mat(k,278) = rxt(k,465)*y(k,217)
         mat(k,1110) = rxt(k,468)*y(k,19) + rxt(k,469)*y(k,59) + rxt(k,471)*y(k,125) &
                      + rxt(k,473)*y(k,134) + rxt(k,474)*y(k,136) + rxt(k,475) &
                      *y(k,217)
         mat(k,1469) = mat(k,1469) + (rxt(k,463)+.500_r8*rxt(k,476))*y(k,67) &
                      + rxt(k,465)*y(k,137) + rxt(k,475)*y(k,150)
         mat(k,194) = -(rxt(k,478)*y(k,229))
         mat(k,2058) = -rxt(k,478)*y(k,152)
         mat(k,757) = rxt(k,477)*y(k,217)
         mat(k,1400) = rxt(k,477)*y(k,151)
         mat(k,778) = .2202005_r8*rxt(k,497)*y(k,134)
         mat(k,805) = .0508005_r8*rxt(k,513)*y(k,134)
         mat(k,1561) = .1279005_r8*rxt(k,496)*y(k,190) + .0097005_r8*rxt(k,501) &
                      *y(k,192) + .0003005_r8*rxt(k,504)*y(k,207) &
                      + .1056005_r8*rxt(k,508)*y(k,208) + .0245005_r8*rxt(k,512) &
                      *y(k,214) + .0154005_r8*rxt(k,518)*y(k,224) &
                      + .0063005_r8*rxt(k,521)*y(k,227)
         mat(k,1996) = .2202005_r8*rxt(k,497)*y(k,6) + .0508005_r8*rxt(k,513)*y(k,110)
         mat(k,34) = .5931005_r8*rxt(k,515)*y(k,217)
         mat(k,40) = .1279005_r8*rxt(k,496)*y(k,124) + .2202005_r8*rxt(k,495)*y(k,203)
         mat(k,46) = .0097005_r8*rxt(k,501)*y(k,124) + .0023005_r8*rxt(k,500)*y(k,203)
         mat(k,1747) = .2202005_r8*rxt(k,495)*y(k,190) + .0023005_r8*rxt(k,500) &
                      *y(k,192) + .0031005_r8*rxt(k,503)*y(k,207) &
                      + .2381005_r8*rxt(k,507)*y(k,208) + .0508005_r8*rxt(k,511) &
                      *y(k,214) + .1364005_r8*rxt(k,517)*y(k,224) &
                      + .1677005_r8*rxt(k,520)*y(k,227)
         mat(k,52) = .0003005_r8*rxt(k,504)*y(k,124) + .0031005_r8*rxt(k,503)*y(k,203)
         mat(k,58) = .1056005_r8*rxt(k,508)*y(k,124) + .2381005_r8*rxt(k,507)*y(k,203)
         mat(k,66) = .0245005_r8*rxt(k,512)*y(k,124) + .0508005_r8*rxt(k,511)*y(k,203)
         mat(k,1366) = .5931005_r8*rxt(k,515)*y(k,172)
         mat(k,72) = .0154005_r8*rxt(k,518)*y(k,124) + .1364005_r8*rxt(k,517)*y(k,203)
         mat(k,78) = .0063005_r8*rxt(k,521)*y(k,124) + .1677005_r8*rxt(k,520)*y(k,203)
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
         mat(k,779) = .2067005_r8*rxt(k,497)*y(k,134)
         mat(k,806) = .1149005_r8*rxt(k,513)*y(k,134)
         mat(k,1562) = .1792005_r8*rxt(k,496)*y(k,190) + .0034005_r8*rxt(k,501) &
                      *y(k,192) + .0003005_r8*rxt(k,504)*y(k,207) &
                      + .1026005_r8*rxt(k,508)*y(k,208) + .0082005_r8*rxt(k,512) &
                      *y(k,214) + .0452005_r8*rxt(k,518)*y(k,224) &
                      + .0237005_r8*rxt(k,521)*y(k,227)
         mat(k,1997) = .2067005_r8*rxt(k,497)*y(k,6) + .1149005_r8*rxt(k,513)*y(k,110)
         mat(k,35) = .1534005_r8*rxt(k,515)*y(k,217)
         mat(k,41) = .1792005_r8*rxt(k,496)*y(k,124) + .2067005_r8*rxt(k,495)*y(k,203)
         mat(k,47) = .0034005_r8*rxt(k,501)*y(k,124) + .0008005_r8*rxt(k,500)*y(k,203)
         mat(k,1748) = .2067005_r8*rxt(k,495)*y(k,190) + .0008005_r8*rxt(k,500) &
                      *y(k,192) + .0035005_r8*rxt(k,503)*y(k,207) &
                      + .1308005_r8*rxt(k,507)*y(k,208) + .1149005_r8*rxt(k,511) &
                      *y(k,214) + .0101005_r8*rxt(k,517)*y(k,224) &
                      + .0174005_r8*rxt(k,520)*y(k,227)
         mat(k,53) = .0003005_r8*rxt(k,504)*y(k,124) + .0035005_r8*rxt(k,503)*y(k,203)
         mat(k,59) = .1026005_r8*rxt(k,508)*y(k,124) + .1308005_r8*rxt(k,507)*y(k,203)
         mat(k,67) = .0082005_r8*rxt(k,512)*y(k,124) + .1149005_r8*rxt(k,511)*y(k,203)
         mat(k,1367) = .1534005_r8*rxt(k,515)*y(k,172)
         mat(k,73) = .0452005_r8*rxt(k,518)*y(k,124) + .0101005_r8*rxt(k,517)*y(k,203)
         mat(k,79) = .0237005_r8*rxt(k,521)*y(k,124) + .0174005_r8*rxt(k,520)*y(k,203)
         mat(k,780) = .0653005_r8*rxt(k,497)*y(k,134)
         mat(k,807) = .0348005_r8*rxt(k,513)*y(k,134)
         mat(k,1563) = .0676005_r8*rxt(k,496)*y(k,190) + .1579005_r8*rxt(k,501) &
                      *y(k,192) + .0073005_r8*rxt(k,504)*y(k,207) &
                      + .0521005_r8*rxt(k,508)*y(k,208) + .0772005_r8*rxt(k,512) &
                      *y(k,214) + .0966005_r8*rxt(k,518)*y(k,224) &
                      + .0025005_r8*rxt(k,521)*y(k,227)
         mat(k,1998) = .0653005_r8*rxt(k,497)*y(k,6) + .0348005_r8*rxt(k,513)*y(k,110)
         mat(k,36) = .0459005_r8*rxt(k,515)*y(k,217)
         mat(k,42) = .0676005_r8*rxt(k,496)*y(k,124) + .0653005_r8*rxt(k,495)*y(k,203)
         mat(k,48) = .1579005_r8*rxt(k,501)*y(k,124) + .0843005_r8*rxt(k,500)*y(k,203)
         mat(k,1749) = .0653005_r8*rxt(k,495)*y(k,190) + .0843005_r8*rxt(k,500) &
                      *y(k,192) + .0003005_r8*rxt(k,503)*y(k,207) &
                      + .0348005_r8*rxt(k,507)*y(k,208) + .0348005_r8*rxt(k,511) &
                      *y(k,214) + .0763005_r8*rxt(k,517)*y(k,224) + .086_r8*rxt(k,520) &
                      *y(k,227)
         mat(k,54) = .0073005_r8*rxt(k,504)*y(k,124) + .0003005_r8*rxt(k,503)*y(k,203)
         mat(k,60) = .0521005_r8*rxt(k,508)*y(k,124) + .0348005_r8*rxt(k,507)*y(k,203)
         mat(k,68) = .0772005_r8*rxt(k,512)*y(k,124) + .0348005_r8*rxt(k,511)*y(k,203)
         mat(k,1368) = .0459005_r8*rxt(k,515)*y(k,172)
         mat(k,74) = .0966005_r8*rxt(k,518)*y(k,124) + .0763005_r8*rxt(k,517)*y(k,203)
         mat(k,80) = .0025005_r8*rxt(k,521)*y(k,124) + .086_r8*rxt(k,520)*y(k,203)
         mat(k,781) = .1749305_r8*rxt(k,494)*y(k,126) + .1284005_r8*rxt(k,497) &
                      *y(k,134)
         mat(k,731) = .0590245_r8*rxt(k,502)*y(k,126) + .0033005_r8*rxt(k,505) &
                      *y(k,134)
         mat(k,808) = .1749305_r8*rxt(k,510)*y(k,126) + .0554005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1564) = .079_r8*rxt(k,496)*y(k,190) + .0059005_r8*rxt(k,501)*y(k,192) &
                      + .0057005_r8*rxt(k,504)*y(k,207) + .0143005_r8*rxt(k,508) &
                      *y(k,208) + .0332005_r8*rxt(k,512)*y(k,214) &
                      + .0073005_r8*rxt(k,518)*y(k,224) + .011_r8*rxt(k,521)*y(k,227)
         mat(k,1939) = .1749305_r8*rxt(k,494)*y(k,6) + .0590245_r8*rxt(k,502)*y(k,98) &
                      + .1749305_r8*rxt(k,510)*y(k,110)
         mat(k,1999) = .1284005_r8*rxt(k,497)*y(k,6) + .0033005_r8*rxt(k,505)*y(k,98) &
                      + .0554005_r8*rxt(k,513)*y(k,110)
         mat(k,37) = .0085005_r8*rxt(k,515)*y(k,217)
         mat(k,43) = .079_r8*rxt(k,496)*y(k,124) + .1284005_r8*rxt(k,495)*y(k,203)
         mat(k,49) = .0059005_r8*rxt(k,501)*y(k,124) + .0443005_r8*rxt(k,500)*y(k,203)
         mat(k,1750) = .1284005_r8*rxt(k,495)*y(k,190) + .0443005_r8*rxt(k,500) &
                      *y(k,192) + .0271005_r8*rxt(k,503)*y(k,207) &
                      + .0076005_r8*rxt(k,507)*y(k,208) + .0554005_r8*rxt(k,511) &
                      *y(k,214) + .2157005_r8*rxt(k,517)*y(k,224) &
                      + .0512005_r8*rxt(k,520)*y(k,227)
         mat(k,55) = .0057005_r8*rxt(k,504)*y(k,124) + .0271005_r8*rxt(k,503)*y(k,203)
         mat(k,61) = .0143005_r8*rxt(k,508)*y(k,124) + .0076005_r8*rxt(k,507)*y(k,203)
         mat(k,69) = .0332005_r8*rxt(k,512)*y(k,124) + .0554005_r8*rxt(k,511)*y(k,203)
         mat(k,1369) = .0085005_r8*rxt(k,515)*y(k,172)
         mat(k,75) = .0073005_r8*rxt(k,518)*y(k,124) + .2157005_r8*rxt(k,517)*y(k,203)
         mat(k,81) = .011_r8*rxt(k,521)*y(k,124) + .0512005_r8*rxt(k,520)*y(k,203)
         mat(k,782) = .5901905_r8*rxt(k,494)*y(k,126) + .114_r8*rxt(k,497)*y(k,134)
         mat(k,732) = .0250245_r8*rxt(k,502)*y(k,126)
         mat(k,809) = .5901905_r8*rxt(k,510)*y(k,126) + .1278005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1565) = .1254005_r8*rxt(k,496)*y(k,190) + .0536005_r8*rxt(k,501) &
                      *y(k,192) + .0623005_r8*rxt(k,504)*y(k,207) &
                      + .0166005_r8*rxt(k,508)*y(k,208) + .130_r8*rxt(k,512)*y(k,214) &
                      + .238_r8*rxt(k,518)*y(k,224) + .1185005_r8*rxt(k,521)*y(k,227)
         mat(k,1940) = .5901905_r8*rxt(k,494)*y(k,6) + .0250245_r8*rxt(k,502)*y(k,98) &
                      + .5901905_r8*rxt(k,510)*y(k,110)
         mat(k,2000) = .114_r8*rxt(k,497)*y(k,6) + .1278005_r8*rxt(k,513)*y(k,110)
         mat(k,38) = .0128005_r8*rxt(k,515)*y(k,217)
         mat(k,44) = .1254005_r8*rxt(k,496)*y(k,124) + .114_r8*rxt(k,495)*y(k,203)
         mat(k,50) = .0536005_r8*rxt(k,501)*y(k,124) + .1621005_r8*rxt(k,500)*y(k,203)
         mat(k,1751) = .114_r8*rxt(k,495)*y(k,190) + .1621005_r8*rxt(k,500)*y(k,192) &
                      + .0474005_r8*rxt(k,503)*y(k,207) + .0113005_r8*rxt(k,507) &
                      *y(k,208) + .1278005_r8*rxt(k,511)*y(k,214) &
                      + .0738005_r8*rxt(k,517)*y(k,224) + .1598005_r8*rxt(k,520) &
                      *y(k,227)
         mat(k,56) = .0623005_r8*rxt(k,504)*y(k,124) + .0474005_r8*rxt(k,503)*y(k,203)
         mat(k,62) = .0166005_r8*rxt(k,508)*y(k,124) + .0113005_r8*rxt(k,507)*y(k,203)
         mat(k,70) = .130_r8*rxt(k,512)*y(k,124) + .1278005_r8*rxt(k,511)*y(k,203)
         mat(k,1370) = .0128005_r8*rxt(k,515)*y(k,172)
         mat(k,76) = .238_r8*rxt(k,518)*y(k,124) + .0738005_r8*rxt(k,517)*y(k,203)
         mat(k,82) = .1185005_r8*rxt(k,521)*y(k,124) + .1598005_r8*rxt(k,520)*y(k,203)
         mat(k,39) = -(rxt(k,515)*y(k,217))
         mat(k,1371) = -rxt(k,515)*y(k,172)
         mat(k,141) = .100_r8*rxt(k,429)*y(k,217)
         mat(k,159) = .230_r8*rxt(k,431)*y(k,217)
         mat(k,1392) = .100_r8*rxt(k,429)*y(k,180) + .230_r8*rxt(k,431)*y(k,182)
         mat(k,504) = -(rxt(k,453)*y(k,217))
         mat(k,1444) = -rxt(k,453)*y(k,174)
         mat(k,1789) = rxt(k,451)*y(k,221)
         mat(k,972) = rxt(k,451)*y(k,203)
         mat(k,522) = -(rxt(k,454)*y(k,217))
         mat(k,1446) = -rxt(k,454)*y(k,175)
         mat(k,1593) = .200_r8*rxt(k,447)*y(k,215) + .200_r8*rxt(k,457)*y(k,222)
         mat(k,1656) = .500_r8*rxt(k,445)*y(k,215)
         mat(k,991) = .200_r8*rxt(k,447)*y(k,124) + .500_r8*rxt(k,445)*y(k,198)
         mat(k,950) = .200_r8*rxt(k,457)*y(k,124)
         mat(k,393) = -(rxt(k,458)*y(k,217))
         mat(k,1429) = -rxt(k,458)*y(k,176)
         mat(k,1781) = rxt(k,456)*y(k,222)
         mat(k,949) = rxt(k,456)*y(k,203)
         mat(k,884) = -(rxt(k,459)*y(k,126) + rxt(k,460)*y(k,217))
         mat(k,1955) = -rxt(k,459)*y(k,177)
         mat(k,1478) = -rxt(k,460)*y(k,177)
         mat(k,791) = .330_r8*rxt(k,440)*y(k,134)
         mat(k,818) = .330_r8*rxt(k,443)*y(k,134)
         mat(k,1611) = .800_r8*rxt(k,447)*y(k,215) + .800_r8*rxt(k,457)*y(k,222)
         mat(k,1955) = mat(k,1955) + rxt(k,448)*y(k,215)
         mat(k,2018) = .330_r8*rxt(k,440)*y(k,6) + .330_r8*rxt(k,443)*y(k,110)
         mat(k,523) = rxt(k,454)*y(k,217)
         mat(k,1663) = .500_r8*rxt(k,445)*y(k,215) + rxt(k,455)*y(k,222)
         mat(k,993) = .800_r8*rxt(k,447)*y(k,124) + rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,198)
         mat(k,1478) = mat(k,1478) + rxt(k,454)*y(k,175)
         mat(k,953) = .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,198)
         mat(k,930) = -(rxt(k,461)*y(k,217))
         mat(k,1482) = -rxt(k,461)*y(k,178)
         mat(k,792) = .300_r8*rxt(k,440)*y(k,134)
         mat(k,819) = .300_r8*rxt(k,443)*y(k,134)
         mat(k,1615) = .900_r8*rxt(k,452)*y(k,221)
         mat(k,2020) = .300_r8*rxt(k,440)*y(k,6) + .300_r8*rxt(k,443)*y(k,110)
         mat(k,1666) = rxt(k,450)*y(k,221)
         mat(k,976) = .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,198)
         mat(k,550) = -(rxt(k,428)*y(k,217))
         mat(k,1449) = -rxt(k,428)*y(k,179)
         mat(k,1792) = rxt(k,426)*y(k,223)
         mat(k,620) = rxt(k,426)*y(k,203)
         mat(k,139) = -(rxt(k,429)*y(k,217))
         mat(k,1390) = -rxt(k,429)*y(k,180)
         mat(k,155) = -(rxt(k,395)*y(k,217))
         mat(k,1393) = -rxt(k,395)*y(k,181)
         mat(k,1760) = rxt(k,392)*y(k,225)
         mat(k,1046) = rxt(k,392)*y(k,203)
         mat(k,160) = -(rxt(k,431)*y(k,217))
         mat(k,1394) = -rxt(k,431)*y(k,182)
         mat(k,591) = -(rxt(k,434)*y(k,217))
         mat(k,1453) = -rxt(k,434)*y(k,183)
         mat(k,1795) = rxt(k,432)*y(k,226)
         mat(k,636) = rxt(k,432)*y(k,203)
         mat(k,168) = -(rxt(k,437)*y(k,217))
         mat(k,1395) = -rxt(k,437)*y(k,184)
         mat(k,161) = .150_r8*rxt(k,431)*y(k,217)
         mat(k,1395) = mat(k,1395) + .150_r8*rxt(k,431)*y(k,182)
         mat(k,352) = -(rxt(k,438)*y(k,217))
         mat(k,1423) = -rxt(k,438)*y(k,185)
         mat(k,1775) = rxt(k,435)*y(k,228)
         mat(k,409) = rxt(k,435)*y(k,203)
         mat(k,429) = -(rxt(k,396)*y(k,203) + rxt(k,397)*y(k,124) + rxt(k,425) &
                      *y(k,125))
         mat(k,1784) = -rxt(k,396)*y(k,188)
         mat(k,1588) = -rxt(k,397)*y(k,188)
         mat(k,1526) = -rxt(k,425)*y(k,188)
         mat(k,191) = rxt(k,402)*y(k,217)
         mat(k,1434) = rxt(k,402)*y(k,22)
         mat(k,837) = -(rxt(k,357)*y(k,203) + (rxt(k,358) + rxt(k,359)) * y(k,124))
         mat(k,1811) = -rxt(k,357)*y(k,189)
         mat(k,1609) = -(rxt(k,358) + rxt(k,359)) * y(k,189)
         mat(k,540) = rxt(k,360)*y(k,217)
         mat(k,185) = rxt(k,361)*y(k,217)
         mat(k,1474) = rxt(k,360)*y(k,2) + rxt(k,361)*y(k,15)
         mat(k,45) = -(rxt(k,495)*y(k,203) + rxt(k,496)*y(k,124))
         mat(k,1752) = -rxt(k,495)*y(k,190)
         mat(k,1566) = -rxt(k,496)*y(k,190)
         mat(k,783) = rxt(k,498)*y(k,217)
         mat(k,1372) = rxt(k,498)*y(k,6)
         mat(k,402) = -(rxt(k,399)*y(k,203) + rxt(k,400)*y(k,124))
         mat(k,1782) = -rxt(k,399)*y(k,191)
         mat(k,1585) = -rxt(k,400)*y(k,191)
         mat(k,116) = .350_r8*rxt(k,398)*y(k,217)
         mat(k,306) = rxt(k,401)*y(k,217)
         mat(k,1430) = .350_r8*rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8)
         mat(k,51) = -(rxt(k,500)*y(k,203) + rxt(k,501)*y(k,124))
         mat(k,1753) = -rxt(k,500)*y(k,192)
         mat(k,1567) = -rxt(k,501)*y(k,192)
         mat(k,112) = rxt(k,499)*y(k,217)
         mat(k,1373) = rxt(k,499)*y(k,7)
         mat(k,360) = -(rxt(k,403)*y(k,203) + rxt(k,405)*y(k,124))
         mat(k,1776) = -rxt(k,403)*y(k,193)
         mat(k,1580) = -rxt(k,405)*y(k,193)
         mat(k,264) = rxt(k,404)*y(k,217)
         mat(k,142) = .070_r8*rxt(k,429)*y(k,217)
         mat(k,162) = .060_r8*rxt(k,431)*y(k,217)
         mat(k,1424) = rxt(k,404)*y(k,23) + .070_r8*rxt(k,429)*y(k,180) &
                      + .060_r8*rxt(k,431)*y(k,182)
         mat(k,715) = -(4._r8*rxt(k,280)*y(k,194) + rxt(k,281)*y(k,198) + rxt(k,282) &
                      *y(k,203) + rxt(k,283)*y(k,124))
         mat(k,1659) = -rxt(k,281)*y(k,194)
         mat(k,1806) = -rxt(k,282)*y(k,194)
         mat(k,1605) = -rxt(k,283)*y(k,194)
         mat(k,269) = .500_r8*rxt(k,285)*y(k,217)
         mat(k,232) = rxt(k,286)*y(k,56) + rxt(k,287)*y(k,217)
         mat(k,1912) = rxt(k,286)*y(k,28)
         mat(k,1465) = .500_r8*rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28)
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
         mat(k,686) = -(rxt(k,309)*y(k,198) + rxt(k,310)*y(k,203) + rxt(k,311) &
                      *y(k,124))
         mat(k,1657) = -rxt(k,309)*y(k,195)
         mat(k,1803) = -rxt(k,310)*y(k,195)
         mat(k,1603) = -rxt(k,311)*y(k,195)
         mat(k,347) = rxt(k,312)*y(k,217)
         mat(k,94) = rxt(k,313)*y(k,217)
         mat(k,1461) = rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31)
         mat(k,530) = -(rxt(k,406)*y(k,203) + rxt(k,407)*y(k,124))
         mat(k,1790) = -rxt(k,406)*y(k,196)
         mat(k,1594) = -rxt(k,407)*y(k,196)
         mat(k,204) = rxt(k,408)*y(k,217)
         mat(k,1594) = mat(k,1594) + rxt(k,397)*y(k,188)
         mat(k,2007) = rxt(k,423)*y(k,141)
         mat(k,390) = rxt(k,423)*y(k,134)
         mat(k,430) = rxt(k,397)*y(k,124) + .400_r8*rxt(k,396)*y(k,203)
         mat(k,1790) = mat(k,1790) + .400_r8*rxt(k,396)*y(k,188)
         mat(k,1447) = rxt(k,408)*y(k,32)
         mat(k,1264) = -(4._r8*rxt(k,291)*y(k,197) + rxt(k,292)*y(k,198) + rxt(k,293) &
                      *y(k,203) + rxt(k,294)*y(k,124) + rxt(k,305)*y(k,125) + rxt(k,332) &
                      *y(k,209) + rxt(k,365)*y(k,205) + rxt(k,370)*y(k,206) + rxt(k,379) &
                      *y(k,101) + rxt(k,390)*y(k,225))
         mat(k,1683) = -rxt(k,292)*y(k,197)
         mat(k,1833) = -rxt(k,293)*y(k,197)
         mat(k,1632) = -rxt(k,294)*y(k,197)
         mat(k,1542) = -rxt(k,305)*y(k,197)
         mat(k,1189) = -rxt(k,332)*y(k,197)
         mat(k,1216) = -rxt(k,365)*y(k,197)
         mat(k,1146) = -rxt(k,370)*y(k,197)
         mat(k,1076) = -rxt(k,379)*y(k,197)
         mat(k,1054) = -rxt(k,390)*y(k,197)
         mat(k,798) = .060_r8*rxt(k,440)*y(k,134)
         mat(k,1025) = rxt(k,288)*y(k,126) + rxt(k,289)*y(k,217)
         mat(k,1101) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,217)
         mat(k,436) = .500_r8*rxt(k,296)*y(k,217)
         mat(k,743) = .080_r8*rxt(k,385)*y(k,134)
         mat(k,1092) = .100_r8*rxt(k,338)*y(k,134)
         mat(k,825) = .060_r8*rxt(k,443)*y(k,134)
         mat(k,1166) = .280_r8*rxt(k,352)*y(k,134)
         mat(k,1632) = mat(k,1632) + .530_r8*rxt(k,336)*y(k,209) + rxt(k,345)*y(k,211) &
                      + rxt(k,348)*y(k,213) + rxt(k,323)*y(k,220)
         mat(k,1977) = rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) + .530_r8*rxt(k,335) &
                      *y(k,209) + rxt(k,346)*y(k,211)
         mat(k,2037) = .060_r8*rxt(k,440)*y(k,6) + .080_r8*rxt(k,385)*y(k,98) &
                      + .100_r8*rxt(k,338)*y(k,105) + .060_r8*rxt(k,443)*y(k,110) &
                      + .280_r8*rxt(k,352)*y(k,111)
         mat(k,933) = .650_r8*rxt(k,461)*y(k,217)
         mat(k,1264) = mat(k,1264) + .530_r8*rxt(k,332)*y(k,209)
         mat(k,1683) = mat(k,1683) + .260_r8*rxt(k,333)*y(k,209) + rxt(k,342)*y(k,211) &
                      + .300_r8*rxt(k,321)*y(k,220)
         mat(k,1833) = mat(k,1833) + .450_r8*rxt(k,343)*y(k,211) + .200_r8*rxt(k,347) &
                      *y(k,213) + .150_r8*rxt(k,322)*y(k,220)
         mat(k,1189) = mat(k,1189) + .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335) &
                      *y(k,126) + .530_r8*rxt(k,332)*y(k,197) + .260_r8*rxt(k,333) &
                      *y(k,198)
         mat(k,1234) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,342)*y(k,198) &
                      + .450_r8*rxt(k,343)*y(k,203) + 4.000_r8*rxt(k,344)*y(k,211)
         mat(k,574) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,203)
         mat(k,1501) = rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49) + .500_r8*rxt(k,296) &
                      *y(k,51) + .650_r8*rxt(k,461)*y(k,178)
         mat(k,1015) = rxt(k,323)*y(k,124) + .300_r8*rxt(k,321)*y(k,198) &
                      + .150_r8*rxt(k,322)*y(k,203)
         mat(k,1691) = -(rxt(k,181)*y(k,59) + (4._r8*rxt(k,258) + 4._r8*rxt(k,259) &
                      ) * y(k,198) + rxt(k,260)*y(k,203) + rxt(k,261)*y(k,124) &
                      + rxt(k,281)*y(k,194) + rxt(k,292)*y(k,197) + rxt(k,309) &
                      *y(k,195) + rxt(k,321)*y(k,220) + rxt(k,333)*y(k,209) + rxt(k,342) &
                      *y(k,211) + rxt(k,366)*y(k,205) + rxt(k,371)*y(k,206) + rxt(k,380) &
                      *y(k,101) + rxt(k,391)*y(k,225) + rxt(k,445)*y(k,215) + rxt(k,450) &
                      *y(k,221) + rxt(k,455)*y(k,222))
         mat(k,1332) = -rxt(k,181)*y(k,198)
         mat(k,1842) = -rxt(k,260)*y(k,198)
         mat(k,1641) = -rxt(k,261)*y(k,198)
         mat(k,720) = -rxt(k,281)*y(k,198)
         mat(k,1270) = -rxt(k,292)*y(k,198)
         mat(k,692) = -rxt(k,309)*y(k,198)
         mat(k,1019) = -rxt(k,321)*y(k,198)
         mat(k,1195) = -rxt(k,333)*y(k,198)
         mat(k,1240) = -rxt(k,342)*y(k,198)
         mat(k,1222) = -rxt(k,366)*y(k,198)
         mat(k,1152) = -rxt(k,371)*y(k,198)
         mat(k,1082) = -rxt(k,380)*y(k,198)
         mat(k,1059) = -rxt(k,391)*y(k,198)
         mat(k,1005) = -rxt(k,445)*y(k,198)
         mat(k,986) = -rxt(k,450)*y(k,198)
         mat(k,966) = -rxt(k,455)*y(k,198)
         mat(k,878) = .280_r8*rxt(k,308)*y(k,134)
         mat(k,475) = rxt(k,295)*y(k,217)
         mat(k,325) = .700_r8*rxt(k,263)*y(k,217)
         mat(k,745) = .050_r8*rxt(k,385)*y(k,134)
         mat(k,1082) = mat(k,1082) + rxt(k,379)*y(k,197)
         mat(k,1641) = mat(k,1641) + rxt(k,294)*y(k,197) + .830_r8*rxt(k,411)*y(k,199) &
                      + .170_r8*rxt(k,417)*y(k,212)
         mat(k,2046) = .280_r8*rxt(k,308)*y(k,29) + .050_r8*rxt(k,385)*y(k,98)
         mat(k,1270) = mat(k,1270) + rxt(k,379)*y(k,101) + rxt(k,294)*y(k,124) &
                      + 4.000_r8*rxt(k,291)*y(k,197) + .900_r8*rxt(k,292)*y(k,198) &
                      + .450_r8*rxt(k,293)*y(k,203) + rxt(k,365)*y(k,205) + rxt(k,370) &
                      *y(k,206) + rxt(k,332)*y(k,209) + rxt(k,341)*y(k,211) &
                      + rxt(k,390)*y(k,225)
         mat(k,1691) = mat(k,1691) + .900_r8*rxt(k,292)*y(k,197)
         mat(k,655) = .830_r8*rxt(k,411)*y(k,124) + .330_r8*rxt(k,410)*y(k,203)
         mat(k,1842) = mat(k,1842) + .450_r8*rxt(k,293)*y(k,197) + .330_r8*rxt(k,410) &
                      *y(k,199) + .070_r8*rxt(k,416)*y(k,212)
         mat(k,1222) = mat(k,1222) + rxt(k,365)*y(k,197)
         mat(k,1152) = mat(k,1152) + rxt(k,370)*y(k,197)
         mat(k,1195) = mat(k,1195) + rxt(k,332)*y(k,197)
         mat(k,1240) = mat(k,1240) + rxt(k,341)*y(k,197)
         mat(k,776) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,203)
         mat(k,1510) = rxt(k,295)*y(k,50) + .700_r8*rxt(k,263)*y(k,53)
         mat(k,1059) = mat(k,1059) + rxt(k,390)*y(k,197)
         mat(k,649) = -(rxt(k,410)*y(k,203) + rxt(k,411)*y(k,124) + rxt(k,412) &
                      *y(k,125))
         mat(k,1800) = -rxt(k,410)*y(k,199)
         mat(k,1601) = -rxt(k,411)*y(k,199)
         mat(k,1531) = -rxt(k,412)*y(k,199)
         mat(k,477) = -((rxt(k,329) + rxt(k,330)) * y(k,124))
         mat(k,1590) = -(rxt(k,329) + rxt(k,330)) * y(k,200)
         mat(k,285) = rxt(k,328)*y(k,217)
         mat(k,1440) = rxt(k,328)*y(k,16)
         mat(k,1575) = .750_r8*rxt(k,298)*y(k,202)
         mat(k,603) = .750_r8*rxt(k,298)*y(k,124)
         mat(k,604) = -(rxt(k,297)*y(k,203) + rxt(k,298)*y(k,124))
         mat(k,1796) = -rxt(k,297)*y(k,202)
         mat(k,1597) = -rxt(k,298)*y(k,202)
         mat(k,466) = rxt(k,304)*y(k,217)
         mat(k,1454) = rxt(k,304)*y(k,25)
         mat(k,1845) = -((rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,76) + rxt(k,138) &
                      *y(k,133) + rxt(k,139)*y(k,134) + rxt(k,143)*y(k,217) &
                      + 4._r8*rxt(k,148)*y(k,203) + rxt(k,158)*y(k,126) + rxt(k,163) &
                      *y(k,124) + rxt(k,168)*y(k,125) + (rxt(k,178) + rxt(k,179) &
                      ) * y(k,56) + rxt(k,185)*y(k,59) + rxt(k,211)*y(k,17) + rxt(k,217) &
                      *y(k,19) + rxt(k,254)*y(k,42) + rxt(k,260)*y(k,198) + rxt(k,268) &
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
         mat(k,1309) = -(rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,203)
         mat(k,1875) = -rxt(k,138)*y(k,203)
         mat(k,2049) = -rxt(k,139)*y(k,203)
         mat(k,1513) = -rxt(k,143)*y(k,203)
         mat(k,1989) = -rxt(k,158)*y(k,203)
         mat(k,1644) = -rxt(k,163)*y(k,203)
         mat(k,1554) = -rxt(k,168)*y(k,203)
         mat(k,1932) = -(rxt(k,178) + rxt(k,179)) * y(k,203)
         mat(k,1335) = -rxt(k,185)*y(k,203)
         mat(k,1284) = -rxt(k,211)*y(k,203)
         mat(k,1718) = -rxt(k,217)*y(k,203)
         mat(k,1740) = -rxt(k,254)*y(k,203)
         mat(k,1694) = -rxt(k,260)*y(k,203)
         mat(k,371) = -rxt(k,268)*y(k,203)
         mat(k,722) = -rxt(k,282)*y(k,203)
         mat(k,1272) = -rxt(k,293)*y(k,203)
         mat(k,610) = -rxt(k,297)*y(k,203)
         mat(k,694) = -rxt(k,310)*y(k,203)
         mat(k,665) = -rxt(k,318)*y(k,203)
         mat(k,1021) = -rxt(k,322)*y(k,203)
         mat(k,1197) = -rxt(k,334)*y(k,203)
         mat(k,1242) = -rxt(k,343)*y(k,203)
         mat(k,578) = -rxt(k,347)*y(k,203)
         mat(k,846) = -rxt(k,357)*y(k,203)
         mat(k,1224) = -rxt(k,367)*y(k,203)
         mat(k,1154) = -rxt(k,372)*y(k,203)
         mat(k,1084) = -rxt(k,381)*y(k,203)
         mat(k,1061) = -rxt(k,392)*y(k,203)
         mat(k,434) = -rxt(k,396)*y(k,203)
         mat(k,408) = -rxt(k,399)*y(k,203)
         mat(k,365) = -rxt(k,403)*y(k,203)
         mat(k,534) = -rxt(k,406)*y(k,203)
         mat(k,656) = -rxt(k,410)*y(k,203)
         mat(k,616) = -rxt(k,413)*y(k,203)
         mat(k,777) = -rxt(k,416)*y(k,203)
         mat(k,378) = -rxt(k,419)*y(k,203)
         mat(k,631) = -rxt(k,426)*y(k,203)
         mat(k,648) = -rxt(k,432)*y(k,203)
         mat(k,416) = -rxt(k,435)*y(k,203)
         mat(k,1007) = -rxt(k,446)*y(k,203)
         mat(k,988) = -rxt(k,451)*y(k,203)
         mat(k,968) = -rxt(k,456)*y(k,203)
         mat(k,801) = .570_r8*rxt(k,440)*y(k,134)
         mat(k,118) = .650_r8*rxt(k,398)*y(k,217)
         mat(k,1284) = mat(k,1284) + rxt(k,210)*y(k,42)
         mat(k,1718) = mat(k,1718) + rxt(k,222)*y(k,217)
         mat(k,230) = .350_r8*rxt(k,277)*y(k,217)
         mat(k,471) = .130_r8*rxt(k,279)*y(k,134)
         mat(k,201) = rxt(k,284)*y(k,217)
         mat(k,880) = .280_r8*rxt(k,308)*y(k,134)
         mat(k,1740) = mat(k,1740) + rxt(k,210)*y(k,17) + rxt(k,174)*y(k,56) &
                      + rxt(k,255)*y(k,126) + rxt(k,256)*y(k,133)
         mat(k,89) = rxt(k,290)*y(k,217)
         mat(k,699) = rxt(k,262)*y(k,217)
         mat(k,1932) = mat(k,1932) + rxt(k,174)*y(k,42) + rxt(k,177)*y(k,79)
         mat(k,1335) = mat(k,1335) + rxt(k,181)*y(k,198) + rxt(k,192)*y(k,217)
         mat(k,942) = rxt(k,265)*y(k,217)
         mat(k,150) = .730_r8*rxt(k,409)*y(k,217)
         mat(k,214) = .500_r8*rxt(k,476)*y(k,217)
         mat(k,912) = rxt(k,301)*y(k,217)
         mat(k,768) = rxt(k,302)*y(k,217)
         mat(k,499) = rxt(k,177)*y(k,56) + rxt(k,133)*y(k,133) + rxt(k,142)*y(k,217)
         mat(k,133) = rxt(k,266)*y(k,217)
         mat(k,702) = rxt(k,267)*y(k,217)
         mat(k,927) = rxt(k,331)*y(k,217)
         mat(k,948) = rxt(k,316)*y(k,217)
         mat(k,747) = .370_r8*rxt(k,385)*y(k,134)
         mat(k,521) = .300_r8*rxt(k,376)*y(k,217)
         mat(k,464) = rxt(k,377)*y(k,217)
         mat(k,1084) = mat(k,1084) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) &
                      + rxt(k,379)*y(k,197) + 1.200_r8*rxt(k,380)*y(k,198)
         mat(k,333) = rxt(k,384)*y(k,217)
         mat(k,1096) = .140_r8*rxt(k,338)*y(k,134)
         mat(k,241) = .200_r8*rxt(k,340)*y(k,217)
         mat(k,491) = .500_r8*rxt(k,351)*y(k,217)
         mat(k,828) = .570_r8*rxt(k,443)*y(k,134)
         mat(k,1174) = .280_r8*rxt(k,352)*y(k,134)
         mat(k,303) = rxt(k,388)*y(k,217)
         mat(k,906) = rxt(k,389)*y(k,217)
         mat(k,1644) = mat(k,1644) + rxt(k,382)*y(k,101) + rxt(k,358)*y(k,189) &
                      + rxt(k,400)*y(k,191) + rxt(k,405)*y(k,193) + rxt(k,283) &
                      *y(k,194) + rxt(k,311)*y(k,195) + rxt(k,261)*y(k,198) &
                      + .170_r8*rxt(k,411)*y(k,199) + rxt(k,329)*y(k,200) &
                      + .250_r8*rxt(k,298)*y(k,202) + rxt(k,270)*y(k,204) &
                      + .920_r8*rxt(k,368)*y(k,205) + .920_r8*rxt(k,374)*y(k,206) &
                      + .470_r8*rxt(k,336)*y(k,209) + .400_r8*rxt(k,414)*y(k,210) &
                      + .830_r8*rxt(k,417)*y(k,212) + rxt(k,420)*y(k,218) + rxt(k,319) &
                      *y(k,219) + .900_r8*rxt(k,452)*y(k,221) + .800_r8*rxt(k,457) &
                      *y(k,222) + rxt(k,427)*y(k,223) + rxt(k,393)*y(k,225) &
                      + rxt(k,433)*y(k,226) + rxt(k,436)*y(k,228)
         mat(k,1989) = mat(k,1989) + rxt(k,255)*y(k,42) + rxt(k,383)*y(k,101) &
                      + rxt(k,369)*y(k,205) + rxt(k,375)*y(k,206) + .470_r8*rxt(k,335) &
                      *y(k,209) + rxt(k,161)*y(k,217) + rxt(k,394)*y(k,225)
         mat(k,1875) = mat(k,1875) + rxt(k,256)*y(k,42) + rxt(k,133)*y(k,79)
         mat(k,2049) = mat(k,2049) + .570_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .280_r8*rxt(k,308)*y(k,29) + .370_r8*rxt(k,385) &
                      *y(k,98) + .140_r8*rxt(k,338)*y(k,105) + .570_r8*rxt(k,443) &
                      *y(k,110) + .280_r8*rxt(k,352)*y(k,111) + rxt(k,145)*y(k,217)
         mat(k,127) = .800_r8*rxt(k,421)*y(k,217)
         mat(k,761) = rxt(k,477)*y(k,217)
         mat(k,937) = .200_r8*rxt(k,461)*y(k,217)
         mat(k,145) = .280_r8*rxt(k,429)*y(k,217)
         mat(k,167) = .380_r8*rxt(k,431)*y(k,217)
         mat(k,172) = .630_r8*rxt(k,437)*y(k,217)
         mat(k,846) = mat(k,846) + rxt(k,358)*y(k,124)
         mat(k,408) = mat(k,408) + rxt(k,400)*y(k,124)
         mat(k,365) = mat(k,365) + rxt(k,405)*y(k,124)
         mat(k,722) = mat(k,722) + rxt(k,283)*y(k,124) + 2.400_r8*rxt(k,280)*y(k,194) &
                      + rxt(k,281)*y(k,198)
         mat(k,694) = mat(k,694) + rxt(k,311)*y(k,124) + rxt(k,309)*y(k,198)
         mat(k,1272) = mat(k,1272) + rxt(k,379)*y(k,101) + .900_r8*rxt(k,292)*y(k,198) &
                      + rxt(k,365)*y(k,205) + rxt(k,370)*y(k,206) + .470_r8*rxt(k,332) &
                      *y(k,209) + rxt(k,390)*y(k,225)
         mat(k,1694) = mat(k,1694) + rxt(k,181)*y(k,59) + 1.200_r8*rxt(k,380)*y(k,101) &
                      + rxt(k,261)*y(k,124) + rxt(k,281)*y(k,194) + rxt(k,309) &
                      *y(k,195) + .900_r8*rxt(k,292)*y(k,197) + 4.000_r8*rxt(k,258) &
                      *y(k,198) + rxt(k,366)*y(k,205) + rxt(k,371)*y(k,206) &
                      + .730_r8*rxt(k,333)*y(k,209) + rxt(k,342)*y(k,211) &
                      + .500_r8*rxt(k,445)*y(k,215) + .300_r8*rxt(k,321)*y(k,220) &
                      + rxt(k,450)*y(k,221) + rxt(k,455)*y(k,222) + .800_r8*rxt(k,391) &
                      *y(k,225)
         mat(k,656) = mat(k,656) + .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410) &
                      *y(k,203)
         mat(k,484) = rxt(k,329)*y(k,124)
         mat(k,610) = mat(k,610) + .250_r8*rxt(k,298)*y(k,124)
         mat(k,1845) = mat(k,1845) + .070_r8*rxt(k,410)*y(k,199) + .160_r8*rxt(k,413) &
                      *y(k,210) + .330_r8*rxt(k,416)*y(k,212)
         mat(k,371) = mat(k,371) + rxt(k,270)*y(k,124)
         mat(k,1224) = mat(k,1224) + .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) &
                      + rxt(k,365)*y(k,197) + rxt(k,366)*y(k,198)
         mat(k,1154) = mat(k,1154) + .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,197) + rxt(k,371)*y(k,198)
         mat(k,1197) = mat(k,1197) + .470_r8*rxt(k,336)*y(k,124) + .470_r8*rxt(k,335) &
                      *y(k,126) + .470_r8*rxt(k,332)*y(k,197) + .730_r8*rxt(k,333) &
                      *y(k,198)
         mat(k,616) = mat(k,616) + .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413) &
                      *y(k,203)
         mat(k,1242) = mat(k,1242) + rxt(k,342)*y(k,198)
         mat(k,777) = mat(k,777) + .830_r8*rxt(k,417)*y(k,124) + .330_r8*rxt(k,416) &
                      *y(k,203)
         mat(k,1007) = mat(k,1007) + .500_r8*rxt(k,445)*y(k,198)
         mat(k,1513) = mat(k,1513) + .650_r8*rxt(k,398)*y(k,7) + rxt(k,222)*y(k,19) &
                      + .350_r8*rxt(k,277)*y(k,24) + rxt(k,284)*y(k,26) + rxt(k,290) &
                      *y(k,47) + rxt(k,262)*y(k,52) + rxt(k,192)*y(k,59) + rxt(k,265) &
                      *y(k,62) + .730_r8*rxt(k,409)*y(k,66) + .500_r8*rxt(k,476) &
                      *y(k,67) + rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) + rxt(k,142) &
                      *y(k,79) + rxt(k,266)*y(k,86) + rxt(k,267)*y(k,87) + rxt(k,331) &
                      *y(k,93) + rxt(k,316)*y(k,95) + .300_r8*rxt(k,376)*y(k,99) &
                      + rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) + .200_r8*rxt(k,340) &
                      *y(k,106) + .500_r8*rxt(k,351)*y(k,109) + rxt(k,388)*y(k,115) &
                      + rxt(k,389)*y(k,116) + rxt(k,161)*y(k,126) + rxt(k,145) &
                      *y(k,134) + .800_r8*rxt(k,421)*y(k,142) + rxt(k,477)*y(k,151) &
                      + .200_r8*rxt(k,461)*y(k,178) + .280_r8*rxt(k,429)*y(k,180) &
                      + .380_r8*rxt(k,431)*y(k,182) + .630_r8*rxt(k,437)*y(k,184)
         mat(k,378) = mat(k,378) + rxt(k,420)*y(k,124)
         mat(k,665) = mat(k,665) + rxt(k,319)*y(k,124)
         mat(k,1021) = mat(k,1021) + .300_r8*rxt(k,321)*y(k,198)
         mat(k,988) = mat(k,988) + .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,198)
         mat(k,968) = mat(k,968) + .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,198)
         mat(k,631) = mat(k,631) + rxt(k,427)*y(k,124)
         mat(k,1061) = mat(k,1061) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126) &
                      + rxt(k,390)*y(k,197) + .800_r8*rxt(k,391)*y(k,198)
         mat(k,648) = mat(k,648) + rxt(k,433)*y(k,124)
         mat(k,416) = mat(k,416) + rxt(k,436)*y(k,124)
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
         mat(k,366) = -(rxt(k,268)*y(k,203) + rxt(k,270)*y(k,124))
         mat(k,1777) = -rxt(k,268)*y(k,204)
         mat(k,1581) = -rxt(k,270)*y(k,204)
         mat(k,1725) = rxt(k,254)*y(k,203)
         mat(k,1777) = mat(k,1777) + rxt(k,254)*y(k,42)
         mat(k,1214) = -(rxt(k,365)*y(k,197) + rxt(k,366)*y(k,198) + rxt(k,367) &
                      *y(k,203) + rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126))
         mat(k,1262) = -rxt(k,365)*y(k,205)
         mat(k,1681) = -rxt(k,366)*y(k,205)
         mat(k,1831) = -rxt(k,367)*y(k,205)
         mat(k,1630) = -rxt(k,368)*y(k,205)
         mat(k,1975) = -rxt(k,369)*y(k,205)
         mat(k,742) = .600_r8*rxt(k,386)*y(k,217)
         mat(k,1499) = .600_r8*rxt(k,386)*y(k,98)
         mat(k,1142) = -(rxt(k,370)*y(k,197) + rxt(k,371)*y(k,198) + rxt(k,372) &
                      *y(k,203) + rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126))
         mat(k,1259) = -rxt(k,370)*y(k,206)
         mat(k,1678) = -rxt(k,371)*y(k,206)
         mat(k,1828) = -rxt(k,372)*y(k,206)
         mat(k,1627) = -rxt(k,374)*y(k,206)
         mat(k,1972) = -rxt(k,375)*y(k,206)
         mat(k,740) = .400_r8*rxt(k,386)*y(k,217)
         mat(k,1496) = .400_r8*rxt(k,386)*y(k,98)
         mat(k,57) = -(rxt(k,503)*y(k,203) + rxt(k,504)*y(k,124))
         mat(k,1754) = -rxt(k,503)*y(k,207)
         mat(k,1568) = -rxt(k,504)*y(k,207)
         mat(k,733) = rxt(k,506)*y(k,217)
         mat(k,1374) = rxt(k,506)*y(k,98)
         mat(k,63) = -(rxt(k,507)*y(k,203) + rxt(k,508)*y(k,124))
         mat(k,1755) = -rxt(k,507)*y(k,208)
         mat(k,1569) = -rxt(k,508)*y(k,208)
         mat(k,64) = rxt(k,509)*y(k,217)
         mat(k,1375) = rxt(k,509)*y(k,104)
         mat(k,1187) = -(rxt(k,332)*y(k,197) + rxt(k,333)*y(k,198) + rxt(k,334) &
                      *y(k,203) + rxt(k,335)*y(k,126) + (rxt(k,336) + rxt(k,337) &
                      ) * y(k,124))
         mat(k,1261) = -rxt(k,332)*y(k,209)
         mat(k,1680) = -rxt(k,333)*y(k,209)
         mat(k,1830) = -rxt(k,334)*y(k,209)
         mat(k,1974) = -rxt(k,335)*y(k,209)
         mat(k,1629) = -(rxt(k,336) + rxt(k,337)) * y(k,209)
         mat(k,1090) = .500_r8*rxt(k,339)*y(k,217)
         mat(k,238) = .200_r8*rxt(k,340)*y(k,217)
         mat(k,1165) = rxt(k,353)*y(k,217)
         mat(k,1498) = .500_r8*rxt(k,339)*y(k,105) + .200_r8*rxt(k,340)*y(k,106) &
                      + rxt(k,353)*y(k,111)
         mat(k,611) = -(rxt(k,413)*y(k,203) + rxt(k,414)*y(k,124) + rxt(k,415) &
                      *y(k,125))
         mat(k,1797) = -rxt(k,413)*y(k,210)
         mat(k,1598) = -rxt(k,414)*y(k,210)
         mat(k,1530) = -rxt(k,415)*y(k,210)
         mat(k,1233) = -(rxt(k,341)*y(k,197) + rxt(k,342)*y(k,198) + rxt(k,343) &
                      *y(k,203) + 4._r8*rxt(k,344)*y(k,211) + rxt(k,345)*y(k,124) &
                      + rxt(k,346)*y(k,126) + rxt(k,354)*y(k,125))
         mat(k,1263) = -rxt(k,341)*y(k,211)
         mat(k,1682) = -rxt(k,342)*y(k,211)
         mat(k,1832) = -rxt(k,343)*y(k,211)
         mat(k,1631) = -rxt(k,345)*y(k,211)
         mat(k,1976) = -rxt(k,346)*y(k,211)
         mat(k,1541) = -rxt(k,354)*y(k,211)
         mat(k,1091) = .500_r8*rxt(k,339)*y(k,217)
         mat(k,239) = .500_r8*rxt(k,340)*y(k,217)
         mat(k,1500) = .500_r8*rxt(k,339)*y(k,105) + .500_r8*rxt(k,340)*y(k,106)
         mat(k,770) = -(rxt(k,416)*y(k,203) + rxt(k,417)*y(k,124) + rxt(k,418) &
                      *y(k,125))
         mat(k,1810) = -rxt(k,416)*y(k,212)
         mat(k,1608) = -rxt(k,417)*y(k,212)
         mat(k,1535) = -rxt(k,418)*y(k,212)
         mat(k,572) = -(rxt(k,347)*y(k,203) + rxt(k,348)*y(k,124))
         mat(k,1793) = -rxt(k,347)*y(k,213)
         mat(k,1596) = -rxt(k,348)*y(k,213)
         mat(k,424) = rxt(k,349)*y(k,217)
         mat(k,243) = rxt(k,350)*y(k,217)
         mat(k,1451) = rxt(k,349)*y(k,107) + rxt(k,350)*y(k,108)
         mat(k,71) = -(rxt(k,511)*y(k,203) + rxt(k,512)*y(k,124))
         mat(k,1756) = -rxt(k,511)*y(k,214)
         mat(k,1570) = -rxt(k,512)*y(k,214)
         mat(k,810) = rxt(k,514)*y(k,217)
         mat(k,1377) = rxt(k,514)*y(k,110)
         mat(k,997) = -(rxt(k,445)*y(k,198) + rxt(k,446)*y(k,203) + rxt(k,447) &
                      *y(k,124) + rxt(k,448)*y(k,126))
         mat(k,1671) = -rxt(k,445)*y(k,215)
         mat(k,1820) = -rxt(k,446)*y(k,215)
         mat(k,1620) = -rxt(k,447)*y(k,215)
         mat(k,1964) = -rxt(k,448)*y(k,215)
         mat(k,795) = rxt(k,439)*y(k,126)
         mat(k,822) = rxt(k,442)*y(k,126)
         mat(k,1964) = mat(k,1964) + rxt(k,439)*y(k,6) + rxt(k,442)*y(k,110) &
                      + .500_r8*rxt(k,459)*y(k,177)
         mat(k,318) = rxt(k,449)*y(k,217)
         mat(k,888) = .500_r8*rxt(k,459)*y(k,126)
         mat(k,1487) = rxt(k,449)*y(k,128)
         mat(k,1352) = -(rxt(k,124)*y(k,77) + rxt(k,125)*y(k,229) + rxt(k,128) &
                      *y(k,134) + (rxt(k,206) + rxt(k,207)) * y(k,85) + (rxt(k,229) &
                      + rxt(k,230)) * y(k,81) + rxt(k,235)*y(k,64) + rxt(k,236) &
                      *y(k,65) + rxt(k,274)*y(k,86))
         mat(k,1038) = -rxt(k,124)*y(k,216)
         mat(k,2067) = -rxt(k,125)*y(k,216)
         mat(k,2042) = -rxt(k,128)*y(k,216)
         mat(k,1891) = -(rxt(k,206) + rxt(k,207)) * y(k,216)
         mat(k,707) = -(rxt(k,229) + rxt(k,230)) * y(k,216)
         mat(k,111) = -rxt(k,235)*y(k,216)
         mat(k,136) = -rxt(k,236)*y(k,216)
         mat(k,131) = -rxt(k,274)*y(k,216)
         mat(k,1507) = -(rxt(k,141)*y(k,77) + rxt(k,142)*y(k,79) + rxt(k,143)*y(k,203) &
                      + rxt(k,144)*y(k,133) + rxt(k,145)*y(k,134) + (4._r8*rxt(k,146) &
                      + 4._r8*rxt(k,147)) * y(k,217) + rxt(k,149)*y(k,90) + rxt(k,161) &
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
                      *y(k,177) + rxt(k,461)*y(k,178) + (rxt(k,463) + rxt(k,476) &
                      ) * y(k,67) + rxt(k,465)*y(k,137) + rxt(k,470)*y(k,148) &
                      + rxt(k,475)*y(k,150) + rxt(k,477)*y(k,151) + rxt(k,479) &
                      *y(k,120))
         mat(k,1039) = -rxt(k,141)*y(k,217)
         mat(k,498) = -rxt(k,142)*y(k,217)
         mat(k,1839) = -rxt(k,143)*y(k,217)
         mat(k,1869) = -rxt(k,144)*y(k,217)
         mat(k,2043) = -rxt(k,145)*y(k,217)
         mat(k,335) = -rxt(k,149)*y(k,217)
         mat(k,1983) = -rxt(k,161)*y(k,217)
         mat(k,294) = -rxt(k,162)*y(k,217)
         mat(k,1548) = -rxt(k,170)*y(k,217)
         mat(k,1294) = -rxt(k,171)*y(k,217)
         mat(k,858) = -rxt(k,190)*y(k,217)
         mat(k,1329) = -(rxt(k,192) + rxt(k,193)) * y(k,217)
         mat(k,1892) = -rxt(k,195)*y(k,217)
         mat(k,726) = -rxt(k,198)*y(k,217)
         mat(k,1712) = -rxt(k,222)*y(k,217)
         mat(k,708) = -rxt(k,224)*y(k,217)
         mat(k,1734) = -rxt(k,257)*y(k,217)
         mat(k,697) = -rxt(k,262)*y(k,217)
         mat(k,324) = -rxt(k,263)*y(k,217)
         mat(k,941) = -(rxt(k,265) + rxt(k,275)) * y(k,217)
         mat(k,132) = -rxt(k,266)*y(k,217)
         mat(k,701) = -rxt(k,267)*y(k,217)
         mat(k,229) = -rxt(k,277)*y(k,217)
         mat(k,200) = -rxt(k,284)*y(k,217)
         mat(k,271) = -rxt(k,285)*y(k,217)
         mat(k,233) = -rxt(k,287)*y(k,217)
         mat(k,1028) = -rxt(k,289)*y(k,217)
         mat(k,88) = -rxt(k,290)*y(k,217)
         mat(k,474) = -rxt(k,295)*y(k,217)
         mat(k,437) = -rxt(k,296)*y(k,217)
         mat(k,910) = -rxt(k,301)*y(k,217)
         mat(k,767) = -rxt(k,302)*y(k,217)
         mat(k,384) = -rxt(k,303)*y(k,217)
         mat(k,469) = -rxt(k,304)*y(k,217)
         mat(k,349) = -rxt(k,312)*y(k,217)
         mat(k,95) = -rxt(k,313)*y(k,217)
         mat(k,1104) = -rxt(k,315)*y(k,217)
         mat(k,946) = -rxt(k,316)*y(k,217)
         mat(k,754) = -rxt(k,317)*y(k,217)
         mat(k,445) = -rxt(k,320)*y(k,217)
         mat(k,343) = -rxt(k,324)*y(k,217)
         mat(k,875) = -rxt(k,325)*y(k,217)
         mat(k,850) = -rxt(k,326)*y(k,217)
         mat(k,288) = -rxt(k,328)*y(k,217)
         mat(k,923) = -rxt(k,331)*y(k,217)
         mat(k,1094) = -rxt(k,339)*y(k,217)
         mat(k,240) = -rxt(k,340)*y(k,217)
         mat(k,427) = -rxt(k,349)*y(k,217)
         mat(k,246) = -rxt(k,350)*y(k,217)
         mat(k,488) = -rxt(k,351)*y(k,217)
         mat(k,1169) = -rxt(k,353)*y(k,217)
         mat(k,567) = -rxt(k,356)*y(k,217)
         mat(k,544) = -rxt(k,360)*y(k,217)
         mat(k,186) = -rxt(k,361)*y(k,217)
         mat(k,179) = -rxt(k,362)*y(k,217)
         mat(k,275) = -rxt(k,363)*y(k,217)
         mat(k,105) = -rxt(k,364)*y(k,217)
         mat(k,518) = -rxt(k,376)*y(k,217)
         mat(k,463) = -rxt(k,377)*y(k,217)
         mat(k,331) = -rxt(k,384)*y(k,217)
         mat(k,744) = -rxt(k,386)*y(k,217)
         mat(k,584) = -rxt(k,387)*y(k,217)
         mat(k,302) = -rxt(k,388)*y(k,217)
         mat(k,901) = -rxt(k,389)*y(k,217)
         mat(k,157) = -rxt(k,395)*y(k,217)
         mat(k,117) = -rxt(k,398)*y(k,217)
         mat(k,308) = -rxt(k,401)*y(k,217)
         mat(k,192) = -rxt(k,402)*y(k,217)
         mat(k,266) = -rxt(k,404)*y(k,217)
         mat(k,205) = -rxt(k,408)*y(k,217)
         mat(k,149) = -rxt(k,409)*y(k,217)
         mat(k,126) = -rxt(k,421)*y(k,217)
         mat(k,260) = -rxt(k,424)*y(k,217)
         mat(k,558) = -rxt(k,428)*y(k,217)
         mat(k,144) = -rxt(k,429)*y(k,217)
         mat(k,166) = -rxt(k,431)*y(k,217)
         mat(k,600) = -rxt(k,434)*y(k,217)
         mat(k,171) = -rxt(k,437)*y(k,217)
         mat(k,356) = -rxt(k,438)*y(k,217)
         mat(k,799) = -rxt(k,441)*y(k,217)
         mat(k,826) = -rxt(k,444)*y(k,217)
         mat(k,320) = -rxt(k,449)*y(k,217)
         mat(k,510) = -rxt(k,453)*y(k,217)
         mat(k,525) = -rxt(k,454)*y(k,217)
         mat(k,397) = -rxt(k,458)*y(k,217)
         mat(k,889) = -rxt(k,460)*y(k,217)
         mat(k,934) = -rxt(k,461)*y(k,217)
         mat(k,213) = -(rxt(k,463) + rxt(k,476)) * y(k,217)
         mat(k,282) = -rxt(k,465)*y(k,217)
         mat(k,420) = -rxt(k,470)*y(k,217)
         mat(k,1115) = -rxt(k,475)*y(k,217)
         mat(k,760) = -rxt(k,477)*y(k,217)
         mat(k,85) = -rxt(k,479)*y(k,217)
         mat(k,799) = mat(k,799) + .630_r8*rxt(k,440)*y(k,134)
         mat(k,229) = mat(k,229) + .650_r8*rxt(k,277)*y(k,217)
         mat(k,469) = mat(k,469) + .130_r8*rxt(k,279)*y(k,134)
         mat(k,271) = mat(k,271) + .500_r8*rxt(k,285)*y(k,217)
         mat(k,875) = mat(k,875) + .360_r8*rxt(k,308)*y(k,134)
         mat(k,1734) = mat(k,1734) + rxt(k,256)*y(k,133)
         mat(k,324) = mat(k,324) + .300_r8*rxt(k,263)*y(k,217)
         mat(k,1926) = rxt(k,179)*y(k,203)
         mat(k,681) = rxt(k,233)*y(k,229)
         mat(k,1306) = rxt(k,140)*y(k,134) + 2.000_r8*rxt(k,135)*y(k,203)
         mat(k,1039) = mat(k,1039) + rxt(k,132)*y(k,133) + rxt(k,124)*y(k,216)
         mat(k,498) = mat(k,498) + rxt(k,133)*y(k,133)
         mat(k,708) = mat(k,708) + rxt(k,223)*y(k,133) + rxt(k,229)*y(k,216)
         mat(k,1892) = mat(k,1892) + rxt(k,194)*y(k,133) + rxt(k,206)*y(k,216)
         mat(k,132) = mat(k,132) + rxt(k,274)*y(k,216)
         mat(k,670) = rxt(k,225)*y(k,133)
         mat(k,726) = mat(k,726) + rxt(k,197)*y(k,133)
         mat(k,744) = mat(k,744) + .320_r8*rxt(k,385)*y(k,134)
         mat(k,584) = mat(k,584) + .600_r8*rxt(k,387)*y(k,217)
         mat(k,1094) = mat(k,1094) + .240_r8*rxt(k,338)*y(k,134)
         mat(k,240) = mat(k,240) + .100_r8*rxt(k,340)*y(k,217)
         mat(k,826) = mat(k,826) + .630_r8*rxt(k,443)*y(k,134)
         mat(k,1169) = mat(k,1169) + .360_r8*rxt(k,352)*y(k,134)
         mat(k,1638) = rxt(k,163)*y(k,203)
         mat(k,1983) = mat(k,1983) + rxt(k,158)*y(k,203)
         mat(k,1869) = mat(k,1869) + rxt(k,256)*y(k,42) + rxt(k,132)*y(k,77) &
                      + rxt(k,133)*y(k,79) + rxt(k,223)*y(k,81) + rxt(k,194)*y(k,85) &
                      + rxt(k,225)*y(k,91) + rxt(k,197)*y(k,92) + rxt(k,138)*y(k,203)
         mat(k,2043) = mat(k,2043) + .630_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .360_r8*rxt(k,308)*y(k,29) + rxt(k,140)*y(k,76) &
                      + .320_r8*rxt(k,385)*y(k,98) + .240_r8*rxt(k,338)*y(k,105) &
                      + .630_r8*rxt(k,443)*y(k,110) + .360_r8*rxt(k,352)*y(k,111) &
                      + rxt(k,139)*y(k,203)
         mat(k,445) = mat(k,445) + .500_r8*rxt(k,320)*y(k,217)
         mat(k,157) = mat(k,157) + .500_r8*rxt(k,395)*y(k,217)
         mat(k,431) = .400_r8*rxt(k,396)*y(k,203)
         mat(k,1267) = .450_r8*rxt(k,293)*y(k,203)
         mat(k,652) = .400_r8*rxt(k,410)*y(k,203)
         mat(k,1839) = mat(k,1839) + rxt(k,179)*y(k,56) + 2.000_r8*rxt(k,135)*y(k,76) &
                      + rxt(k,163)*y(k,124) + rxt(k,158)*y(k,126) + rxt(k,138) &
                      *y(k,133) + rxt(k,139)*y(k,134) + .400_r8*rxt(k,396)*y(k,188) &
                      + .450_r8*rxt(k,293)*y(k,197) + .400_r8*rxt(k,410)*y(k,199) &
                      + .450_r8*rxt(k,343)*y(k,211) + .400_r8*rxt(k,416)*y(k,212) &
                      + .200_r8*rxt(k,347)*y(k,213) + .150_r8*rxt(k,322)*y(k,220)
         mat(k,1237) = .450_r8*rxt(k,343)*y(k,203)
         mat(k,773) = .400_r8*rxt(k,416)*y(k,203)
         mat(k,575) = .200_r8*rxt(k,347)*y(k,203)
         mat(k,1353) = rxt(k,124)*y(k,77) + rxt(k,229)*y(k,81) + rxt(k,206)*y(k,85) &
                      + rxt(k,274)*y(k,86) + 2.000_r8*rxt(k,125)*y(k,229)
         mat(k,1507) = mat(k,1507) + .650_r8*rxt(k,277)*y(k,24) + .500_r8*rxt(k,285) &
                      *y(k,27) + .300_r8*rxt(k,263)*y(k,53) + .600_r8*rxt(k,387) &
                      *y(k,103) + .100_r8*rxt(k,340)*y(k,106) + .500_r8*rxt(k,320) &
                      *y(k,146) + .500_r8*rxt(k,395)*y(k,181)
         mat(k,1016) = .150_r8*rxt(k,322)*y(k,203)
         mat(k,2068) = rxt(k,233)*y(k,73) + 2.000_r8*rxt(k,125)*y(k,216)
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
         mat(k,373) = -(rxt(k,419)*y(k,203) + rxt(k,420)*y(k,124))
         mat(k,1778) = -rxt(k,419)*y(k,218)
         mat(k,1582) = -rxt(k,420)*y(k,218)
         mat(k,147) = .200_r8*rxt(k,409)*y(k,217)
         mat(k,124) = .140_r8*rxt(k,421)*y(k,217)
         mat(k,258) = rxt(k,424)*y(k,217)
         mat(k,1425) = .200_r8*rxt(k,409)*y(k,66) + .140_r8*rxt(k,421)*y(k,142) &
                      + rxt(k,424)*y(k,143)
         mat(k,658) = -(rxt(k,318)*y(k,203) + rxt(k,319)*y(k,124))
         mat(k,1801) = -rxt(k,318)*y(k,219)
         mat(k,1602) = -rxt(k,319)*y(k,219)
         mat(k,866) = rxt(k,325)*y(k,217)
         mat(k,442) = .500_r8*rxt(k,320)*y(k,217)
         mat(k,1459) = rxt(k,325)*y(k,29) + .500_r8*rxt(k,320)*y(k,146)
         mat(k,1013) = -(rxt(k,321)*y(k,198) + rxt(k,322)*y(k,203) + rxt(k,323) &
                      *y(k,124))
         mat(k,1672) = -rxt(k,321)*y(k,220)
         mat(k,1821) = -rxt(k,322)*y(k,220)
         mat(k,1621) = -rxt(k,323)*y(k,220)
         mat(k,796) = .060_r8*rxt(k,440)*y(k,134)
         mat(k,848) = rxt(k,326)*y(k,217)
         mat(k,823) = .060_r8*rxt(k,443)*y(k,134)
         mat(k,2026) = .060_r8*rxt(k,440)*y(k,6) + .060_r8*rxt(k,443)*y(k,110)
         mat(k,341) = rxt(k,324)*y(k,217)
         mat(k,932) = .150_r8*rxt(k,461)*y(k,217)
         mat(k,1488) = rxt(k,326)*y(k,48) + rxt(k,324)*y(k,147) + .150_r8*rxt(k,461) &
                      *y(k,178)
         mat(k,978) = -(rxt(k,450)*y(k,198) + rxt(k,451)*y(k,203) + rxt(k,452) &
                      *y(k,124))
         mat(k,1670) = -rxt(k,450)*y(k,221)
         mat(k,1819) = -rxt(k,451)*y(k,221)
         mat(k,1619) = -rxt(k,452)*y(k,221)
         mat(k,1963) = .500_r8*rxt(k,459)*y(k,177)
         mat(k,509) = rxt(k,453)*y(k,217)
         mat(k,887) = .500_r8*rxt(k,459)*y(k,126) + rxt(k,460)*y(k,217)
         mat(k,1486) = rxt(k,453)*y(k,174) + rxt(k,460)*y(k,177)
         mat(k,956) = -(rxt(k,455)*y(k,198) + rxt(k,456)*y(k,203) + rxt(k,457) &
                      *y(k,124))
         mat(k,1669) = -rxt(k,455)*y(k,222)
         mat(k,1818) = -rxt(k,456)*y(k,222)
         mat(k,1618) = -rxt(k,457)*y(k,222)
         mat(k,794) = rxt(k,441)*y(k,217)
         mat(k,821) = rxt(k,444)*y(k,217)
         mat(k,396) = rxt(k,458)*y(k,217)
         mat(k,1485) = rxt(k,441)*y(k,6) + rxt(k,444)*y(k,110) + rxt(k,458)*y(k,176)
         mat(k,622) = -(rxt(k,426)*y(k,203) + rxt(k,427)*y(k,124))
         mat(k,1798) = -rxt(k,426)*y(k,223)
         mat(k,1599) = -rxt(k,427)*y(k,223)
         mat(k,552) = rxt(k,428)*y(k,217)
         mat(k,143) = .650_r8*rxt(k,429)*y(k,217)
         mat(k,1456) = rxt(k,428)*y(k,179) + .650_r8*rxt(k,429)*y(k,180)
         mat(k,77) = -(rxt(k,517)*y(k,203) + rxt(k,518)*y(k,124))
         mat(k,1757) = -rxt(k,517)*y(k,224)
         mat(k,1571) = -rxt(k,518)*y(k,224)
         mat(k,138) = rxt(k,516)*y(k,217)
         mat(k,1378) = rxt(k,516)*y(k,180)
         mat(k,1052) = -(rxt(k,390)*y(k,197) + rxt(k,391)*y(k,198) + rxt(k,392) &
                      *y(k,203) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126))
         mat(k,1255) = -rxt(k,390)*y(k,225)
         mat(k,1674) = -rxt(k,391)*y(k,225)
         mat(k,1824) = -rxt(k,392)*y(k,225)
         mat(k,1623) = -rxt(k,393)*y(k,225)
         mat(k,1967) = -rxt(k,394)*y(k,225)
         mat(k,178) = rxt(k,362)*y(k,217)
         mat(k,274) = rxt(k,363)*y(k,217)
         mat(k,104) = rxt(k,364)*y(k,217)
         mat(k,581) = .400_r8*rxt(k,387)*y(k,217)
         mat(k,156) = .500_r8*rxt(k,395)*y(k,217)
         mat(k,1491) = rxt(k,362)*y(k,94) + rxt(k,363)*y(k,96) + rxt(k,364)*y(k,97) &
                      + .400_r8*rxt(k,387)*y(k,103) + .500_r8*rxt(k,395)*y(k,181)
         mat(k,638) = -(rxt(k,432)*y(k,203) + rxt(k,433)*y(k,124))
         mat(k,1799) = -rxt(k,432)*y(k,226)
         mat(k,1600) = -rxt(k,433)*y(k,226)
         mat(k,163) = .560_r8*rxt(k,431)*y(k,217)
         mat(k,593) = rxt(k,434)*y(k,217)
         mat(k,1457) = .560_r8*rxt(k,431)*y(k,182) + rxt(k,434)*y(k,183)
         mat(k,83) = -(rxt(k,520)*y(k,203) + rxt(k,521)*y(k,124))
         mat(k,1758) = -rxt(k,520)*y(k,227)
         mat(k,1572) = -rxt(k,521)*y(k,227)
         mat(k,158) = rxt(k,519)*y(k,217)
         mat(k,1379) = rxt(k,519)*y(k,182)
         mat(k,410) = -(rxt(k,435)*y(k,203) + rxt(k,436)*y(k,124))
         mat(k,1783) = -rxt(k,435)*y(k,228)
         mat(k,1586) = -rxt(k,436)*y(k,228)
         mat(k,170) = .300_r8*rxt(k,437)*y(k,217)
         mat(k,353) = rxt(k,438)*y(k,217)
         mat(k,1431) = .300_r8*rxt(k,437)*y(k,184) + rxt(k,438)*y(k,185)
         mat(k,2080) = -(rxt(k,125)*y(k,216) + rxt(k,233)*y(k,73) + rxt(k,478) &
                      *y(k,152))
         mat(k,1365) = -rxt(k,125)*y(k,229)
         mat(k,684) = -rxt(k,233)*y(k,229)
         mat(k,197) = -rxt(k,478)*y(k,229)
         mat(k,236) = rxt(k,287)*y(k,217)
         mat(k,351) = rxt(k,312)*y(k,217)
         mat(k,96) = rxt(k,313)*y(k,217)
         mat(k,1746) = rxt(k,257)*y(k,217)
         mat(k,1032) = rxt(k,289)*y(k,217)
         mat(k,852) = rxt(k,326)*y(k,217)
         mat(k,1107) = rxt(k,315)*y(k,217)
         mat(k,476) = rxt(k,295)*y(k,217)
         mat(k,440) = rxt(k,296)*y(k,217)
         mat(k,327) = rxt(k,263)*y(k,217)
         mat(k,1315) = rxt(k,136)*y(k,203)
         mat(k,1045) = rxt(k,141)*y(k,217)
         mat(k,503) = rxt(k,142)*y(k,217)
         mat(k,711) = rxt(k,224)*y(k,217)
         mat(k,1904) = (rxt(k,530)+rxt(k,535))*y(k,91) + (rxt(k,523)+rxt(k,529) &
                       +rxt(k,534))*y(k,92) + rxt(k,195)*y(k,217)
         mat(k,703) = rxt(k,267)*y(k,217)
         mat(k,1301) = rxt(k,171)*y(k,217)
         mat(k,339) = rxt(k,149)*y(k,217)
         mat(k,675) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,730) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85) + rxt(k,198)*y(k,217)
         mat(k,1098) = .500_r8*rxt(k,339)*y(k,217)
         mat(k,86) = rxt(k,479)*y(k,217)
         mat(k,448) = rxt(k,320)*y(k,217)
         mat(k,345) = rxt(k,324)*y(k,217)
         mat(k,1851) = rxt(k,136)*y(k,76) + rxt(k,143)*y(k,217)
         mat(k,1519) = rxt(k,287)*y(k,28) + rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31) &
                      + rxt(k,257)*y(k,42) + rxt(k,289)*y(k,45) + rxt(k,326)*y(k,48) &
                      + rxt(k,315)*y(k,49) + rxt(k,295)*y(k,50) + rxt(k,296)*y(k,51) &
                      + rxt(k,263)*y(k,53) + rxt(k,141)*y(k,77) + rxt(k,142)*y(k,79) &
                      + rxt(k,224)*y(k,81) + rxt(k,195)*y(k,85) + rxt(k,267)*y(k,87) &
                      + rxt(k,171)*y(k,89) + rxt(k,149)*y(k,90) + rxt(k,198)*y(k,92) &
                      + .500_r8*rxt(k,339)*y(k,105) + rxt(k,479)*y(k,120) + rxt(k,320) &
                      *y(k,146) + rxt(k,324)*y(k,147) + rxt(k,143)*y(k,203) &
                      + 2.000_r8*rxt(k,146)*y(k,217)
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
         mat(k, 51) = mat(k, 51) + lmat(k, 51)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 90) = lmat(k, 90)
         mat(k, 91) = lmat(k, 91)
         mat(k, 92) = lmat(k, 92)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = lmat(k, 99)
         mat(k, 100) = lmat(k, 100)
         mat(k, 101) = lmat(k, 101)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = lmat(k, 107)
         mat(k, 108) = lmat(k, 108)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 119) = lmat(k, 119)
         mat(k, 120) = lmat(k, 120)
         mat(k, 121) = lmat(k, 121)
         mat(k, 122) = lmat(k, 122)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 128) = lmat(k, 128)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 151) = lmat(k, 151)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 173) = lmat(k, 173)
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = lmat(k, 175)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 181) = lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 187) = lmat(k, 187)
         mat(k, 188) = lmat(k, 188)
         mat(k, 189) = lmat(k, 189)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = lmat(k, 195)
         mat(k, 196) = lmat(k, 196)
         mat(k, 198) = mat(k, 198) + lmat(k, 198)
         mat(k, 202) = mat(k, 202) + lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = lmat(k, 209)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 216) = lmat(k, 216)
         mat(k, 217) = lmat(k, 217)
         mat(k, 218) = lmat(k, 218)
         mat(k, 219) = lmat(k, 219)
         mat(k, 220) = lmat(k, 220)
         mat(k, 221) = lmat(k, 221)
         mat(k, 222) = lmat(k, 222)
         mat(k, 223) = lmat(k, 223)
         mat(k, 224) = lmat(k, 224)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 247) = lmat(k, 247)
         mat(k, 248) = lmat(k, 248)
         mat(k, 249) = lmat(k, 249)
         mat(k, 250) = lmat(k, 250)
         mat(k, 251) = lmat(k, 251)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 256) = lmat(k, 256)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 259) = lmat(k, 259)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 261) = lmat(k, 261)
         mat(k, 262) = lmat(k, 262)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 267) = lmat(k, 267)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 273) = mat(k, 273) + lmat(k, 273)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 301) = lmat(k, 301)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 305) = lmat(k, 305)
         mat(k, 307) = lmat(k, 307)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 310) = lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 312) = lmat(k, 312)
         mat(k, 313) = lmat(k, 313)
         mat(k, 314) = lmat(k, 314)
         mat(k, 315) = lmat(k, 315)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 317) = lmat(k, 317)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = mat(k, 322) + lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 329) = lmat(k, 329)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = mat(k, 343) + lmat(k, 343)
         mat(k, 344) = lmat(k, 344)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 350) = lmat(k, 350)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = lmat(k, 355)
         mat(k, 356) = mat(k, 356) + lmat(k, 356)
         mat(k, 357) = lmat(k, 357)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 370) = lmat(k, 370)
         mat(k, 371) = mat(k, 371) + lmat(k, 371)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 379) = lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = mat(k, 382) + lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = lmat(k, 394)
         mat(k, 395) = lmat(k, 395)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 399) = lmat(k, 399)
         mat(k, 402) = mat(k, 402) + lmat(k, 402)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 421) = lmat(k, 421)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 425) = lmat(k, 425)
         mat(k, 426) = lmat(k, 426)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 437) = mat(k, 437) + lmat(k, 437)
         mat(k, 438) = lmat(k, 438)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 444) = lmat(k, 444)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 446) = lmat(k, 446)
         mat(k, 447) = lmat(k, 447)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 450) = lmat(k, 450)
         mat(k, 451) = lmat(k, 451)
         mat(k, 452) = lmat(k, 452)
         mat(k, 453) = lmat(k, 453)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 462) = lmat(k, 462)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 485) = mat(k, 485) + lmat(k, 485)
         mat(k, 487) = lmat(k, 487)
         mat(k, 489) = lmat(k, 489)
         mat(k, 493) = lmat(k, 493)
         mat(k, 494) = lmat(k, 494)
         mat(k, 495) = lmat(k, 495)
         mat(k, 496) = lmat(k, 496)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 498) = mat(k, 498) + lmat(k, 498)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = lmat(k, 508)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 511) = lmat(k, 511)
         mat(k, 512) = lmat(k, 512)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 517) = lmat(k, 517)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 524) = lmat(k, 524)
         mat(k, 526) = mat(k, 526) + lmat(k, 526)
         mat(k, 527) = lmat(k, 527)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 536) = lmat(k, 536)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 541) = lmat(k, 541)
         mat(k, 542) = lmat(k, 542)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = lmat(k, 546)
         mat(k, 547) = lmat(k, 547)
         mat(k, 548) = lmat(k, 548)
         mat(k, 549) = lmat(k, 549)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 554) = lmat(k, 554)
         mat(k, 557) = lmat(k, 557)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 559) = lmat(k, 559)
         mat(k, 560) = lmat(k, 560)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 565) = mat(k, 565) + lmat(k, 565)
         mat(k, 568) = mat(k, 568) + lmat(k, 568)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
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
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 696) = mat(k, 696) + lmat(k, 696)
         mat(k, 700) = mat(k, 700) + lmat(k, 700)
         mat(k, 704) = mat(k, 704) + lmat(k, 704)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 706) = mat(k, 706) + lmat(k, 706)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 724) = mat(k, 724) + lmat(k, 724)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 729) = mat(k, 729) + lmat(k, 729)
         mat(k, 735) = mat(k, 735) + lmat(k, 735)
         mat(k, 751) = mat(k, 751) + lmat(k, 751)
         mat(k, 753) = lmat(k, 753)
         mat(k, 755) = mat(k, 755) + lmat(k, 755)
         mat(k, 756) = lmat(k, 756)
         mat(k, 758) = mat(k, 758) + lmat(k, 758)
         mat(k, 759) = lmat(k, 759)
         mat(k, 762) = lmat(k, 762)
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
         mat(k, 857) = mat(k, 857) + lmat(k, 857)
         mat(k, 859) = lmat(k, 859)
         mat(k, 862) = mat(k, 862) + lmat(k, 862)
         mat(k, 863) = mat(k, 863) + lmat(k, 863)
         mat(k, 869) = mat(k, 869) + lmat(k, 869)
         mat(k, 884) = mat(k, 884) + lmat(k, 884)
         mat(k, 885) = lmat(k, 885)
         mat(k, 886) = lmat(k, 886)
         mat(k, 890) = lmat(k, 890)
         mat(k, 892) = lmat(k, 892)
         mat(k, 896) = mat(k, 896) + lmat(k, 896)
         mat(k, 900) = lmat(k, 900)
         mat(k, 902) = lmat(k, 902)
         mat(k, 906) = mat(k, 906) + lmat(k, 906)
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
         mat(k, 924) = lmat(k, 924)
         mat(k, 926) = lmat(k, 926)
         mat(k, 927) = mat(k, 927) + lmat(k, 927)
         mat(k, 929) = mat(k, 929) + lmat(k, 929)
         mat(k, 930) = mat(k, 930) + lmat(k, 930)
         mat(k, 931) = mat(k, 931) + lmat(k, 931)
         mat(k, 932) = mat(k, 932) + lmat(k, 932)
         mat(k, 933) = mat(k, 933) + lmat(k, 933)
         mat(k, 936) = mat(k, 936) + lmat(k, 936)
         mat(k, 937) = mat(k, 937) + lmat(k, 937)
         mat(k, 939) = mat(k, 939) + lmat(k, 939)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 945) = lmat(k, 945)
         mat(k, 947) = lmat(k, 947)
         mat(k, 948) = mat(k, 948) + lmat(k, 948)
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
         mat(k,1095) = mat(k,1095) + lmat(k,1095)
         mat(k,1096) = mat(k,1096) + lmat(k,1096)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1101) = mat(k,1101) + lmat(k,1101)
         mat(k,1105) = lmat(k,1105)
         mat(k,1109) = lmat(k,1109)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1111) = mat(k,1111) + lmat(k,1111)
         mat(k,1120) = lmat(k,1120)
         mat(k,1124) = lmat(k,1124)
         mat(k,1142) = mat(k,1142) + lmat(k,1142)
         mat(k,1154) = mat(k,1154) + lmat(k,1154)
         mat(k,1159) = lmat(k,1159)
         mat(k,1160) = mat(k,1160) + lmat(k,1160)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1166) = mat(k,1166) + lmat(k,1166)
         mat(k,1172) = lmat(k,1172)
         mat(k,1187) = mat(k,1187) + lmat(k,1187)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1264) = mat(k,1264) + lmat(k,1264)
         mat(k,1278) = mat(k,1278) + lmat(k,1278)
         mat(k,1291) = mat(k,1291) + lmat(k,1291)
         mat(k,1294) = mat(k,1294) + lmat(k,1294)
         mat(k,1295) = lmat(k,1295)
         mat(k,1304) = mat(k,1304) + lmat(k,1304)
         mat(k,1309) = mat(k,1309) + lmat(k,1309)
         mat(k,1327) = mat(k,1327) + lmat(k,1327)
         mat(k,1336) = mat(k,1336) + lmat(k,1336)
         mat(k,1338) = mat(k,1338) + lmat(k,1338)
         mat(k,1342) = mat(k,1342) + lmat(k,1342)
         mat(k,1344) = mat(k,1344) + lmat(k,1344)
         mat(k,1345) = mat(k,1345) + lmat(k,1345)
         mat(k,1347) = mat(k,1347) + lmat(k,1347)
         mat(k,1348) = mat(k,1348) + lmat(k,1348)
         mat(k,1350) = mat(k,1350) + lmat(k,1350)
         mat(k,1352) = mat(k,1352) + lmat(k,1352)
         mat(k,1353) = mat(k,1353) + lmat(k,1353)
         mat(k,1355) = lmat(k,1355)
         mat(k,1356) = lmat(k,1356)
         mat(k,1358) = lmat(k,1358)
         mat(k,1359) = lmat(k,1359)
         mat(k,1360) = lmat(k,1360)
         mat(k,1362) = mat(k,1362) + lmat(k,1362)
         mat(k,1384) = lmat(k,1384)
         mat(k,1389) = lmat(k,1389)
         mat(k,1502) = mat(k,1502) + lmat(k,1502)
         mat(k,1507) = mat(k,1507) + lmat(k,1507)
         mat(k,1510) = mat(k,1510) + lmat(k,1510)
         mat(k,1513) = mat(k,1513) + lmat(k,1513)
         mat(k,1516) = mat(k,1516) + lmat(k,1516)
         mat(k,1519) = mat(k,1519) + lmat(k,1519)
         mat(k,1544) = mat(k,1544) + lmat(k,1544)
         mat(k,1548) = mat(k,1548) + lmat(k,1548)
         mat(k,1549) = mat(k,1549) + lmat(k,1549)
         mat(k,1550) = mat(k,1550) + lmat(k,1550)
         mat(k,1555) = mat(k,1555) + lmat(k,1555)
         mat(k,1579) = mat(k,1579) + lmat(k,1579)
         mat(k,1640) = mat(k,1640) + lmat(k,1640)
         mat(k,1645) = mat(k,1645) + lmat(k,1645)
         mat(k,1691) = mat(k,1691) + lmat(k,1691)
         mat(k,1707) = mat(k,1707) + lmat(k,1707)
         mat(k,1716) = mat(k,1716) + lmat(k,1716)
         mat(k,1719) = mat(k,1719) + lmat(k,1719)
         mat(k,1728) = mat(k,1728) + lmat(k,1728)
         mat(k,1729) = lmat(k,1729)
         mat(k,1732) = mat(k,1732) + lmat(k,1732)
         mat(k,1739) = mat(k,1739) + lmat(k,1739)
         mat(k,1788) = mat(k,1788) + lmat(k,1788)
         mat(k,1845) = mat(k,1845) + lmat(k,1845)
         mat(k,1876) = mat(k,1876) + lmat(k,1876)
         mat(k,1880) = mat(k,1880) + lmat(k,1880)
         mat(k,1889) = mat(k,1889) + lmat(k,1889)
         mat(k,1900) = mat(k,1900) + lmat(k,1900)
         mat(k,1901) = mat(k,1901) + lmat(k,1901)
         mat(k,1917) = mat(k,1917) + lmat(k,1917)
         mat(k,1921) = lmat(k,1921)
         mat(k,1929) = lmat(k,1929)
         mat(k,1932) = mat(k,1932) + lmat(k,1932)
         mat(k,1934) = mat(k,1934) + lmat(k,1934)
         mat(k,1935) = mat(k,1935) + lmat(k,1935)
         mat(k,1979) = mat(k,1979) + lmat(k,1979)
         mat(k,1984) = mat(k,1984) + lmat(k,1984)
         mat(k,1985) = mat(k,1985) + lmat(k,1985)
         mat(k,1990) = mat(k,1990) + lmat(k,1990)
         mat(k,1993) = mat(k,1993) + lmat(k,1993)
         mat(k,2042) = mat(k,2042) + lmat(k,2042)
         mat(k,2050) = mat(k,2050) + lmat(k,2050)
         mat(k,2054) = mat(k,2054) + lmat(k,2054)
         mat(k,2061) = lmat(k,2061)
         mat(k,2065) = lmat(k,2065)
         mat(k,2067) = mat(k,2067) + lmat(k,2067)
         mat(k,2068) = mat(k,2068) + lmat(k,2068)
         mat(k,2075) = lmat(k,2075)
         mat(k,2080) = mat(k,2080) + lmat(k,2080)
         mat(k, 164) = 0._r8
         mat(k, 165) = 0._r8
         mat(k, 265) = 0._r8
         mat(k, 361) = 0._r8
         mat(k, 362) = 0._r8
         mat(k, 375) = 0._r8
         mat(k, 403) = 0._r8
         mat(k, 405) = 0._r8
         mat(k, 413) = 0._r8
         mat(k, 531) = 0._r8
         mat(k, 535) = 0._r8
         mat(k, 538) = 0._r8
         mat(k, 539) = 0._r8
         mat(k, 543) = 0._r8
         mat(k, 551) = 0._r8
         mat(k, 553) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 556) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 563) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 594) = 0._r8
         mat(k, 596) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 599) = 0._r8
         mat(k, 605) = 0._r8
         mat(k, 606) = 0._r8
         mat(k, 621) = 0._r8
         mat(k, 623) = 0._r8
         mat(k, 625) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 628) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 641) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 644) = 0._r8
         mat(k, 645) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 674) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 695) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 750) = 0._r8
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
         mat(k, 841) = 0._r8
         mat(k, 844) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 871) = 0._r8
         mat(k, 873) = 0._r8
         mat(k, 874) = 0._r8
         mat(k, 876) = 0._r8
         mat(k, 877) = 0._r8
         mat(k, 883) = 0._r8
         mat(k, 895) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 903) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 925) = 0._r8
         mat(k, 928) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 955) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 960) = 0._r8
         mat(k, 961) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 969) = 0._r8
         mat(k, 970) = 0._r8
         mat(k, 979) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 982) = 0._r8
         mat(k, 983) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 995) = 0._r8
         mat(k, 996) = 0._r8
         mat(k, 998) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1000) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1022) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1044) = 0._r8
         mat(k,1055) = 0._r8
         mat(k,1056) = 0._r8
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
         mat(k,1079) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1119) = 0._r8
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
         mat(k,1149) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1171) = 0._r8
         mat(k,1175) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1192) = 0._r8
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
         mat(k,1219) = 0._r8
         mat(k,1226) = 0._r8
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
         mat(k,1296) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1328) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1341) = 0._r8
         mat(k,1349) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1363) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1442) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1460) = 0._r8
         mat(k,1471) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1532) = 0._r8
         mat(k,1534) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1545) = 0._r8
         mat(k,1547) = 0._r8
         mat(k,1551) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1556) = 0._r8
         mat(k,1557) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1634) = 0._r8
         mat(k,1635) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1646) = 0._r8
         mat(k,1650) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1684) = 0._r8
         mat(k,1685) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1692) = 0._r8
         mat(k,1695) = 0._r8
         mat(k,1696) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1699) = 0._r8
         mat(k,1700) = 0._r8
         mat(k,1708) = 0._r8
         mat(k,1709) = 0._r8
         mat(k,1711) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1720) = 0._r8
         mat(k,1722) = 0._r8
         mat(k,1723) = 0._r8
         mat(k,1724) = 0._r8
         mat(k,1726) = 0._r8
         mat(k,1733) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1736) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1738) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1779) = 0._r8
         mat(k,1780) = 0._r8
         mat(k,1808) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1813) = 0._r8
         mat(k,1814) = 0._r8
         mat(k,1815) = 0._r8
         mat(k,1817) = 0._r8
         mat(k,1826) = 0._r8
         mat(k,1829) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1853) = 0._r8
         mat(k,1859) = 0._r8
         mat(k,1865) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1872) = 0._r8
         mat(k,1881) = 0._r8
         mat(k,1887) = 0._r8
         mat(k,1893) = 0._r8
         mat(k,1894) = 0._r8
         mat(k,1895) = 0._r8
         mat(k,1896) = 0._r8
         mat(k,1897) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1902) = 0._r8
         mat(k,1903) = 0._r8
         mat(k,1910) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1922) = 0._r8
         mat(k,1925) = 0._r8
         mat(k,1927) = 0._r8
         mat(k,1928) = 0._r8
         mat(k,1930) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1938) = 0._r8
         mat(k,1945) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1965) = 0._r8
         mat(k,1971) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1981) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1987) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1992) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2012) = 0._r8
         mat(k,2015) = 0._r8
         mat(k,2019) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2023) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2025) = 0._r8
         mat(k,2028) = 0._r8
         mat(k,2032) = 0._r8
         mat(k,2034) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2036) = 0._r8
         mat(k,2039) = 0._r8
         mat(k,2051) = 0._r8
         mat(k,2055) = 0._r8
         mat(k,2060) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2064) = 0._r8
         mat(k,2066) = 0._r8
         mat(k,2069) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2073) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2079) = 0._r8
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
         mat(k, 51) = mat(k, 51) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 63) = mat(k, 63) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 84) = mat(k, 84) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 106) = mat(k, 106) - dti(k)
         mat(k, 109) = mat(k, 109) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 123) = mat(k, 123) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 151) = mat(k, 151) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 173) = mat(k, 173) - dti(k)
         mat(k, 176) = mat(k, 176) - dti(k)
         mat(k, 181) = mat(k, 181) - dti(k)
         mat(k, 184) = mat(k, 184) - dti(k)
         mat(k, 187) = mat(k, 187) - dti(k)
         mat(k, 190) = mat(k, 190) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 198) = mat(k, 198) - dti(k)
         mat(k, 202) = mat(k, 202) - dti(k)
         mat(k, 206) = mat(k, 206) - dti(k)
         mat(k, 210) = mat(k, 210) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 222) = mat(k, 222) - dti(k)
         mat(k, 225) = mat(k, 225) - dti(k)
         mat(k, 231) = mat(k, 231) - dti(k)
         mat(k, 237) = mat(k, 237) - dti(k)
         mat(k, 242) = mat(k, 242) - dti(k)
         mat(k, 247) = mat(k, 247) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 263) = mat(k, 263) - dti(k)
         mat(k, 268) = mat(k, 268) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 284) = mat(k, 284) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 298) = mat(k, 298) - dti(k)
         mat(k, 304) = mat(k, 304) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 322) = mat(k, 322) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 346) = mat(k, 346) - dti(k)
         mat(k, 352) = mat(k, 352) - dti(k)
         mat(k, 360) = mat(k, 360) - dti(k)
         mat(k, 366) = mat(k, 366) - dti(k)
         mat(k, 373) = mat(k, 373) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 382) = mat(k, 382) - dti(k)
         mat(k, 389) = mat(k, 389) - dti(k)
         mat(k, 393) = mat(k, 393) - dti(k)
         mat(k, 402) = mat(k, 402) - dti(k)
         mat(k, 410) = mat(k, 410) - dti(k)
         mat(k, 417) = mat(k, 417) - dti(k)
         mat(k, 423) = mat(k, 423) - dti(k)
         mat(k, 429) = mat(k, 429) - dti(k)
         mat(k, 435) = mat(k, 435) - dti(k)
         mat(k, 441) = mat(k, 441) - dti(k)
         mat(k, 449) = mat(k, 449) - dti(k)
         mat(k, 457) = mat(k, 457) - dti(k)
         mat(k, 465) = mat(k, 465) - dti(k)
         mat(k, 473) = mat(k, 473) - dti(k)
         mat(k, 477) = mat(k, 477) - dti(k)
         mat(k, 485) = mat(k, 485) - dti(k)
         mat(k, 493) = mat(k, 493) - dti(k)
         mat(k, 497) = mat(k, 497) - dti(k)
         mat(k, 504) = mat(k, 504) - dti(k)
         mat(k, 513) = mat(k, 513) - dti(k)
         mat(k, 522) = mat(k, 522) - dti(k)
         mat(k, 530) = mat(k, 530) - dti(k)
         mat(k, 537) = mat(k, 537) - dti(k)
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
         mat(k, 686) = mat(k, 686) - dti(k)
         mat(k, 696) = mat(k, 696) - dti(k)
         mat(k, 700) = mat(k, 700) - dti(k)
         mat(k, 704) = mat(k, 704) - dti(k)
         mat(k, 715) = mat(k, 715) - dti(k)
         mat(k, 724) = mat(k, 724) - dti(k)
         mat(k, 735) = mat(k, 735) - dti(k)
         mat(k, 751) = mat(k, 751) - dti(k)
         mat(k, 758) = mat(k, 758) - dti(k)
         mat(k, 765) = mat(k, 765) - dti(k)
         mat(k, 770) = mat(k, 770) - dti(k)
         mat(k, 788) = mat(k, 788) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 837) = mat(k, 837) - dti(k)
         mat(k, 847) = mat(k, 847) - dti(k)
         mat(k, 855) = mat(k, 855) - dti(k)
         mat(k, 869) = mat(k, 869) - dti(k)
         mat(k, 884) = mat(k, 884) - dti(k)
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
         mat(k,1327) = mat(k,1327) - dti(k)
         mat(k,1352) = mat(k,1352) - dti(k)
         mat(k,1507) = mat(k,1507) - dti(k)
         mat(k,1549) = mat(k,1549) - dti(k)
         mat(k,1640) = mat(k,1640) - dti(k)
         mat(k,1691) = mat(k,1691) - dti(k)
         mat(k,1716) = mat(k,1716) - dti(k)
         mat(k,1739) = mat(k,1739) - dti(k)
         mat(k,1845) = mat(k,1845) - dti(k)
         mat(k,1876) = mat(k,1876) - dti(k)
         mat(k,1900) = mat(k,1900) - dti(k)
         mat(k,1935) = mat(k,1935) - dti(k)
         mat(k,1993) = mat(k,1993) - dti(k)
         mat(k,2054) = mat(k,2054) - dti(k)
         mat(k,2080) = mat(k,2080) - dti(k)
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
