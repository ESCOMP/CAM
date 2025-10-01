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
         mat(k,732) = -(rxt(k,419)*y(k,265))
         mat(k,1792) = -rxt(k,419)*y(k,1)
         mat(k,1913) = rxt(k,422)*y(k,228)
         mat(k,1084) = rxt(k,422)*y(k,125)
         mat(k,711) = -(rxt(k,423)*y(k,265))
         mat(k,1790) = -rxt(k,423)*y(k,2)
         mat(k,1083) = rxt(k,420)*y(k,243)
         mat(k,2322) = rxt(k,420)*y(k,228)
         mat(k,1063) = -(rxt(k,502)*y(k,127) + rxt(k,503)*y(k,137) + rxt(k,504) &
                      *y(k,265))
         mat(k,2399) = -rxt(k,502)*y(k,6)
         mat(k,2083) = -rxt(k,503)*y(k,6)
         mat(k,1818) = -rxt(k,504)*y(k,6)
         mat(k,210) = -(rxt(k,461)*y(k,265))
         mat(k,1721) = -rxt(k,461)*y(k,7)
         mat(k,457) = -(rxt(k,464)*y(k,265))
         mat(k,1761) = -rxt(k,464)*y(k,8)
         mat(k,530) = rxt(k,462)*y(k,243)
         mat(k,2302) = rxt(k,462)*y(k,230)
         mat(k,211) = .120_r8*rxt(k,461)*y(k,265)
         mat(k,1722) = .120_r8*rxt(k,461)*y(k,7)
         mat(k,1061) = .100_r8*rxt(k,503)*y(k,137)
         mat(k,1005) = .100_r8*rxt(k,506)*y(k,137)
         mat(k,2072) = .100_r8*rxt(k,503)*y(k,6) + .100_r8*rxt(k,506)*y(k,111)
         mat(k,1899) = .500_r8*rxt(k,463)*y(k,230) + .200_r8*rxt(k,490)*y(k,272) &
                      + .060_r8*rxt(k,496)*y(k,275)
         mat(k,531) = .500_r8*rxt(k,463)*y(k,125)
         mat(k,784) = .200_r8*rxt(k,490)*y(k,125)
         mat(k,808) = .060_r8*rxt(k,496)*y(k,125)
         mat(k,1893) = .200_r8*rxt(k,490)*y(k,272) + .200_r8*rxt(k,496)*y(k,275)
         mat(k,783) = .200_r8*rxt(k,490)*y(k,125)
         mat(k,806) = .200_r8*rxt(k,496)*y(k,125)
         mat(k,1909) = .200_r8*rxt(k,490)*y(k,272) + .150_r8*rxt(k,496)*y(k,275)
         mat(k,786) = .200_r8*rxt(k,490)*y(k,125)
         mat(k,809) = .150_r8*rxt(k,496)*y(k,125)
         mat(k,1895) = .210_r8*rxt(k,496)*y(k,275)
         mat(k,807) = .210_r8*rxt(k,496)*y(k,125)
         mat(k,269) = -(rxt(k,424)*y(k,265))
         mat(k,1731) = -rxt(k,424)*y(k,15)
         mat(k,1060) = .050_r8*rxt(k,503)*y(k,137)
         mat(k,1004) = .050_r8*rxt(k,506)*y(k,137)
         mat(k,2071) = .050_r8*rxt(k,503)*y(k,6) + .050_r8*rxt(k,506)*y(k,111)
         mat(k,403) = -(rxt(k,390)*y(k,127) + rxt(k,391)*y(k,265))
         mat(k,2392) = -rxt(k,390)*y(k,16)
         mat(k,1753) = -rxt(k,391)*y(k,16)
         mat(k,1556) = -(rxt(k,273)*y(k,42) + rxt(k,274)*y(k,243) + rxt(k,275) &
                      *y(k,137))
         mat(k,2454) = -rxt(k,273)*y(k,17)
         mat(k,2368) = -rxt(k,274)*y(k,17)
         mat(k,2109) = -rxt(k,275)*y(k,17)
         mat(k,1608) = 4.000_r8*rxt(k,276)*y(k,19) + (rxt(k,277)+rxt(k,278))*y(k,59) &
                      + rxt(k,281)*y(k,125) + rxt(k,284)*y(k,135) + rxt(k,532) &
                      *y(k,157) + rxt(k,285)*y(k,265)
         mat(k,180) = rxt(k,263)*y(k,261)
         mat(k,186) = rxt(k,289)*y(k,261)
         mat(k,518) = 2.000_r8*rxt(k,300)*y(k,56) + 2.000_r8*rxt(k,312)*y(k,261) &
                      + 2.000_r8*rxt(k,301)*y(k,265)
         mat(k,683) = rxt(k,302)*y(k,56) + rxt(k,313)*y(k,261) + rxt(k,303)*y(k,265)
         mat(k,453) = 3.000_r8*rxt(k,307)*y(k,56) + 3.000_r8*rxt(k,290)*y(k,261) &
                      + 3.000_r8*rxt(k,308)*y(k,265)
         mat(k,2002) = 2.000_r8*rxt(k,300)*y(k,41) + rxt(k,302)*y(k,43) &
                      + 3.000_r8*rxt(k,307)*y(k,55)
         mat(k,2482) = (rxt(k,277)+rxt(k,278))*y(k,19)
         mat(k,173) = 2.000_r8*rxt(k,291)*y(k,261)
         mat(k,874) = rxt(k,286)*y(k,135) + rxt(k,292)*y(k,261) + rxt(k,287)*y(k,265)
         mat(k,1956) = rxt(k,281)*y(k,19)
         mat(k,2045) = rxt(k,284)*y(k,19) + rxt(k,286)*y(k,81)
         mat(k,1522) = rxt(k,532)*y(k,19)
         mat(k,2225) = rxt(k,263)*y(k,34) + rxt(k,289)*y(k,35) + 2.000_r8*rxt(k,312) &
                      *y(k,41) + rxt(k,313)*y(k,43) + 3.000_r8*rxt(k,290)*y(k,55) &
                      + 2.000_r8*rxt(k,291)*y(k,78) + rxt(k,292)*y(k,81)
         mat(k,1849) = rxt(k,285)*y(k,19) + 2.000_r8*rxt(k,301)*y(k,41) + rxt(k,303) &
                      *y(k,43) + 3.000_r8*rxt(k,308)*y(k,55) + rxt(k,287)*y(k,81)
         mat(k,1601) = rxt(k,279)*y(k,59)
         mat(k,2475) = rxt(k,279)*y(k,19)
         mat(k,1536) = (rxt(k,597)+rxt(k,602))*y(k,91)
         mat(k,831) = (rxt(k,597)+rxt(k,602))*y(k,85)
         mat(k,1610) = -(4._r8*rxt(k,276)*y(k,19) + (rxt(k,277) + rxt(k,278) + rxt(k,279) &
                      ) * y(k,59) + rxt(k,280)*y(k,243) + rxt(k,281)*y(k,125) &
                      + rxt(k,282)*y(k,126) + rxt(k,284)*y(k,135) + rxt(k,285) &
                      *y(k,265) + rxt(k,532)*y(k,157))
         mat(k,2484) = -(rxt(k,277) + rxt(k,278) + rxt(k,279)) * y(k,19)
         mat(k,2370) = -rxt(k,280)*y(k,19)
         mat(k,1958) = -rxt(k,281)*y(k,19)
         mat(k,1678) = -rxt(k,282)*y(k,19)
         mat(k,2047) = -rxt(k,284)*y(k,19)
         mat(k,1851) = -rxt(k,285)*y(k,19)
         mat(k,1524) = -rxt(k,532)*y(k,19)
         mat(k,1558) = rxt(k,275)*y(k,137)
         mat(k,620) = rxt(k,283)*y(k,135)
         mat(k,875) = rxt(k,293)*y(k,261)
         mat(k,835) = rxt(k,288)*y(k,135)
         mat(k,2047) = mat(k,2047) + rxt(k,283)*y(k,20) + rxt(k,288)*y(k,91)
         mat(k,2111) = rxt(k,275)*y(k,17)
         mat(k,2227) = rxt(k,293)*y(k,81)
         mat(k,617) = -(rxt(k,283)*y(k,135))
         mat(k,2025) = -rxt(k,283)*y(k,20)
         mat(k,1603) = rxt(k,282)*y(k,126)
         mat(k,1657) = rxt(k,282)*y(k,19)
         mat(k,278) = -(rxt(k,465)*y(k,265))
         mat(k,1733) = -rxt(k,465)*y(k,22)
         mat(k,1892) = rxt(k,468)*y(k,232)
         mat(k,475) = rxt(k,468)*y(k,125)
         mat(k,370) = -(rxt(k,467)*y(k,265))
         mat(k,1748) = -rxt(k,467)*y(k,23)
         mat(k,476) = rxt(k,466)*y(k,243)
         mat(k,2294) = rxt(k,466)*y(k,232)
         mat(k,321) = -(rxt(k,338)*y(k,56) + rxt(k,339)*y(k,265))
         mat(k,1976) = -rxt(k,338)*y(k,24)
         mat(k,1742) = -rxt(k,339)*y(k,24)
         mat(k,577) = -(rxt(k,340)*y(k,56) + rxt(k,341)*y(k,137) + rxt(k,366)*y(k,265))
         mat(k,1982) = -rxt(k,340)*y(k,25)
         mat(k,2074) = -rxt(k,341)*y(k,25)
         mat(k,1776) = -rxt(k,366)*y(k,25)
         mat(k,304) = -(rxt(k,346)*y(k,265))
         mat(k,1739) = -rxt(k,346)*y(k,26)
         mat(k,943) = .800_r8*rxt(k,342)*y(k,233) + .200_r8*rxt(k,343)*y(k,237)
         mat(k,2149) = .200_r8*rxt(k,343)*y(k,233)
         mat(k,375) = -(rxt(k,347)*y(k,265))
         mat(k,1749) = -rxt(k,347)*y(k,27)
         mat(k,944) = rxt(k,344)*y(k,243)
         mat(k,2295) = rxt(k,344)*y(k,233)
         mat(k,327) = -(rxt(k,348)*y(k,56) + rxt(k,349)*y(k,265))
         mat(k,1977) = -rxt(k,348)*y(k,28)
         mat(k,1743) = -rxt(k,349)*y(k,28)
         mat(k,1178) = -(rxt(k,369)*y(k,127) + rxt(k,370)*y(k,137) + rxt(k,388) &
                      *y(k,265))
         mat(k,2408) = -rxt(k,369)*y(k,29)
         mat(k,2091) = -rxt(k,370)*y(k,29)
         mat(k,1827) = -rxt(k,388)*y(k,29)
         mat(k,929) = .130_r8*rxt(k,448)*y(k,137)
         mat(k,2091) = mat(k,2091) + .130_r8*rxt(k,448)*y(k,98)
         mat(k,433) = -(rxt(k,374)*y(k,265))
         mat(k,1757) = -rxt(k,374)*y(k,30)
         mat(k,979) = rxt(k,372)*y(k,243)
         mat(k,2299) = rxt(k,372)*y(k,234)
         mat(k,333) = -(rxt(k,375)*y(k,265) + rxt(k,378)*y(k,56))
         mat(k,1744) = -rxt(k,375)*y(k,31)
         mat(k,1978) = -rxt(k,378)*y(k,31)
         mat(k,313) = -(rxt(k,471)*y(k,265))
         mat(k,1741) = -rxt(k,471)*y(k,32)
         mat(k,696) = rxt(k,469)*y(k,243)
         mat(k,2292) = rxt(k,469)*y(k,235)
         mat(k,146) = -(rxt(k,262)*y(k,261))
         mat(k,2201) = -rxt(k,262)*y(k,33)
         mat(k,178) = -(rxt(k,263)*y(k,261))
         mat(k,2206) = -rxt(k,263)*y(k,34)
         mat(k,183) = -(rxt(k,289)*y(k,261))
         mat(k,2207) = -rxt(k,289)*y(k,35)
         mat(k,155) = -(rxt(k,264)*y(k,261))
         mat(k,2202) = -rxt(k,264)*y(k,36)
         mat(k,188) = -(rxt(k,265)*y(k,261))
         mat(k,2208) = -rxt(k,265)*y(k,37)
         mat(k,159) = -(rxt(k,266)*y(k,261))
         mat(k,2203) = -rxt(k,266)*y(k,38)
         mat(k,193) = -(rxt(k,267)*y(k,261))
         mat(k,2209) = -rxt(k,267)*y(k,39)
         mat(k,163) = -(rxt(k,268)*y(k,261))
         mat(k,2204) = -rxt(k,268)*y(k,40)
         mat(k,516) = -(rxt(k,300)*y(k,56) + rxt(k,301)*y(k,265) + rxt(k,312)*y(k,261))
         mat(k,1981) = -rxt(k,300)*y(k,41)
         mat(k,1768) = -rxt(k,301)*y(k,41)
         mat(k,2219) = -rxt(k,312)*y(k,41)
         mat(k,2470) = -(rxt(k,237)*y(k,56) + rxt(k,273)*y(k,17) + rxt(k,317)*y(k,243) &
                      + rxt(k,318)*y(k,127) + rxt(k,319)*y(k,135) + rxt(k,320) &
                      *y(k,265))
         mat(k,2018) = -rxt(k,237)*y(k,42)
         mat(k,1567) = -rxt(k,273)*y(k,42)
         mat(k,2384) = -rxt(k,317)*y(k,42)
         mat(k,2444) = -rxt(k,318)*y(k,42)
         mat(k,2061) = -rxt(k,319)*y(k,42)
         mat(k,1865) = -rxt(k,320)*y(k,42)
         mat(k,741) = .400_r8*rxt(k,419)*y(k,265)
         mat(k,1080) = .340_r8*rxt(k,503)*y(k,137)
         mat(k,410) = .500_r8*rxt(k,390)*y(k,127)
         mat(k,584) = rxt(k,341)*y(k,137)
         mat(k,1194) = .500_r8*rxt(k,370)*y(k,137)
         mat(k,670) = .500_r8*rxt(k,358)*y(k,265)
         mat(k,872) = rxt(k,325)*y(k,265)
         mat(k,449) = .300_r8*rxt(k,326)*y(k,265)
         mat(k,1647) = (rxt(k,334)+rxt(k,335))*y(k,261)
         mat(k,2497) = rxt(k,244)*y(k,237)
         mat(k,1215) = .800_r8*rxt(k,363)*y(k,265)
         mat(k,942) = .910_r8*rxt(k,448)*y(k,137)
         mat(k,642) = .300_r8*rxt(k,439)*y(k,265)
         mat(k,1311) = .800_r8*rxt(k,443)*y(k,237)
         mat(k,1323) = .120_r8*rxt(k,401)*y(k,137)
         mat(k,633) = .500_r8*rxt(k,414)*y(k,265)
         mat(k,1024) = .340_r8*rxt(k,506)*y(k,137)
         mat(k,1435) = .600_r8*rxt(k,415)*y(k,137)
         mat(k,1972) = .100_r8*rxt(k,421)*y(k,228) + rxt(k,324)*y(k,237) &
                      + .500_r8*rxt(k,392)*y(k,240) + .500_r8*rxt(k,360)*y(k,242) &
                      + .920_r8*rxt(k,431)*y(k,245) + .250_r8*rxt(k,399)*y(k,250) &
                      + rxt(k,408)*y(k,252) + rxt(k,382)*y(k,268) + rxt(k,386) &
                      *y(k,269) + .340_r8*rxt(k,515)*y(k,270) + .320_r8*rxt(k,520) &
                      *y(k,271) + .250_r8*rxt(k,456)*y(k,274)
         mat(k,2444) = mat(k,2444) + .500_r8*rxt(k,390)*y(k,16) + rxt(k,432)*y(k,245) &
                      + .250_r8*rxt(k,398)*y(k,250) + rxt(k,409)*y(k,252)
         mat(k,2125) = .340_r8*rxt(k,503)*y(k,6) + rxt(k,341)*y(k,25) &
                      + .500_r8*rxt(k,370)*y(k,29) + .910_r8*rxt(k,448)*y(k,98) &
                      + .120_r8*rxt(k,401)*y(k,106) + .340_r8*rxt(k,506)*y(k,111) &
                      + .600_r8*rxt(k,415)*y(k,112)
         mat(k,592) = rxt(k,365)*y(k,265)
         mat(k,1172) = .680_r8*rxt(k,524)*y(k,265)
         mat(k,1097) = .100_r8*rxt(k,421)*y(k,125)
         mat(k,954) = .700_r8*rxt(k,343)*y(k,237)
         mat(k,989) = rxt(k,371)*y(k,237)
         mat(k,1485) = rxt(k,354)*y(k,237) + rxt(k,428)*y(k,245) + .250_r8*rxt(k,395) &
                      *y(k,250) + rxt(k,404)*y(k,252) + .250_r8*rxt(k,453)*y(k,274)
         mat(k,2198) = rxt(k,244)*y(k,59) + .800_r8*rxt(k,443)*y(k,101) + rxt(k,324) &
                      *y(k,125) + .700_r8*rxt(k,343)*y(k,233) + rxt(k,371)*y(k,234) &
                      + rxt(k,354)*y(k,236) + (4.000_r8*rxt(k,321)+2.000_r8*rxt(k,322)) &
                      *y(k,237) + 1.500_r8*rxt(k,429)*y(k,245) + .750_r8*rxt(k,434) &
                      *y(k,246) + .880_r8*rxt(k,396)*y(k,250) + 2.000_r8*rxt(k,405) &
                      *y(k,252) + .750_r8*rxt(k,508)*y(k,260) + .800_r8*rxt(k,384) &
                      *y(k,269) + .930_r8*rxt(k,513)*y(k,270) + .950_r8*rxt(k,518) &
                      *y(k,271) + .800_r8*rxt(k,454)*y(k,274)
         mat(k,616) = .500_r8*rxt(k,392)*y(k,125)
         mat(k,849) = .500_r8*rxt(k,360)*y(k,125)
         mat(k,2384) = mat(k,2384) + .450_r8*rxt(k,406)*y(k,252) + .150_r8*rxt(k,385) &
                      *y(k,269)
         mat(k,1358) = .920_r8*rxt(k,431)*y(k,125) + rxt(k,432)*y(k,127) + rxt(k,428) &
                      *y(k,236) + 1.500_r8*rxt(k,429)*y(k,237)
         mat(k,1391) = .750_r8*rxt(k,434)*y(k,237)
         mat(k,1412) = .250_r8*rxt(k,399)*y(k,125) + .250_r8*rxt(k,398)*y(k,127) &
                      + .250_r8*rxt(k,395)*y(k,236) + .880_r8*rxt(k,396)*y(k,237)
         mat(k,1453) = rxt(k,408)*y(k,125) + rxt(k,409)*y(k,127) + rxt(k,404)*y(k,236) &
                      + 2.000_r8*rxt(k,405)*y(k,237) + .450_r8*rxt(k,406)*y(k,243) &
                      + 4.000_r8*rxt(k,407)*y(k,252)
         mat(k,1162) = .750_r8*rxt(k,508)*y(k,237)
         mat(k,2241) = (rxt(k,334)+rxt(k,335))*y(k,54)
         mat(k,1865) = mat(k,1865) + .400_r8*rxt(k,419)*y(k,1) + .500_r8*rxt(k,358) &
                      *y(k,51) + rxt(k,325)*y(k,52) + .300_r8*rxt(k,326)*y(k,53) &
                      + .800_r8*rxt(k,363)*y(k,74) + .300_r8*rxt(k,439)*y(k,99) &
                      + .500_r8*rxt(k,414)*y(k,110) + rxt(k,365)*y(k,142) &
                      + .680_r8*rxt(k,524)*y(k,217)
         mat(k,867) = rxt(k,382)*y(k,125)
         mat(k,1270) = rxt(k,386)*y(k,125) + .800_r8*rxt(k,384)*y(k,237) &
                      + .150_r8*rxt(k,385)*y(k,243)
         mat(k,1233) = .340_r8*rxt(k,515)*y(k,125) + .930_r8*rxt(k,513)*y(k,237)
         mat(k,1116) = .320_r8*rxt(k,520)*y(k,125) + .950_r8*rxt(k,518)*y(k,237)
         mat(k,1288) = .250_r8*rxt(k,456)*y(k,125) + .250_r8*rxt(k,453)*y(k,236) &
                      + .800_r8*rxt(k,454)*y(k,237)
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
         mat(k,681) = -(rxt(k,302)*y(k,56) + rxt(k,303)*y(k,265) + rxt(k,313)*y(k,261))
         mat(k,1984) = -rxt(k,302)*y(k,43)
         mat(k,1787) = -rxt(k,303)*y(k,43)
         mat(k,2220) = -rxt(k,313)*y(k,43)
         mat(k,167) = -(rxt(k,304)*y(k,265))
         mat(k,1719) = -rxt(k,304)*y(k,44)
         mat(k,1196) = -(rxt(k,350)*y(k,127) + rxt(k,351)*y(k,265))
         mat(k,2409) = -rxt(k,350)*y(k,45)
         mat(k,1828) = -rxt(k,351)*y(k,45)
         mat(k,736) = .800_r8*rxt(k,419)*y(k,265)
         mat(k,406) = rxt(k,390)*y(k,127)
         mat(k,305) = rxt(k,346)*y(k,265)
         mat(k,377) = .500_r8*rxt(k,347)*y(k,265)
         mat(k,1179) = .500_r8*rxt(k,370)*y(k,137)
         mat(k,1416) = .100_r8*rxt(k,415)*y(k,137)
         mat(k,1938) = .400_r8*rxt(k,421)*y(k,228) + rxt(k,345)*y(k,233) &
                      + .270_r8*rxt(k,373)*y(k,234) + rxt(k,392)*y(k,240) + rxt(k,411) &
                      *y(k,254) + rxt(k,382)*y(k,268)
         mat(k,2409) = mat(k,2409) + rxt(k,390)*y(k,16)
         mat(k,2092) = .500_r8*rxt(k,370)*y(k,29) + .100_r8*rxt(k,415)*y(k,112)
         mat(k,1089) = .400_r8*rxt(k,421)*y(k,125)
         mat(k,947) = rxt(k,345)*y(k,125) + 3.200_r8*rxt(k,342)*y(k,233) &
                      + .800_r8*rxt(k,343)*y(k,237)
         mat(k,982) = .270_r8*rxt(k,373)*y(k,125)
         mat(k,2166) = .800_r8*rxt(k,343)*y(k,233)
         mat(k,611) = rxt(k,392)*y(k,125)
         mat(k,2349) = .200_r8*rxt(k,410)*y(k,254)
         mat(k,744) = rxt(k,411)*y(k,125) + .200_r8*rxt(k,410)*y(k,243)
         mat(k,1828) = mat(k,1828) + .800_r8*rxt(k,419)*y(k,1) + rxt(k,346)*y(k,26) &
                      + .500_r8*rxt(k,347)*y(k,27)
         mat(k,860) = rxt(k,382)*y(k,125)
         mat(k,411) = -(rxt(k,305)*y(k,56) + rxt(k,306)*y(k,265))
         mat(k,1979) = -rxt(k,305)*y(k,46)
         mat(k,1754) = -rxt(k,306)*y(k,46)
         mat(k,149) = -(rxt(k,352)*y(k,265))
         mat(k,1718) = -rxt(k,352)*y(k,47)
         mat(k,1125) = -(rxt(k,389)*y(k,265))
         mat(k,1823) = -rxt(k,389)*y(k,48)
         mat(k,735) = .800_r8*rxt(k,419)*y(k,265)
         mat(k,1068) = .520_r8*rxt(k,503)*y(k,137)
         mat(k,405) = .500_r8*rxt(k,390)*y(k,127)
         mat(k,1012) = .520_r8*rxt(k,506)*y(k,137)
         mat(k,1934) = .250_r8*rxt(k,421)*y(k,228) + .820_r8*rxt(k,373)*y(k,234) &
                      + .500_r8*rxt(k,392)*y(k,240) + .270_r8*rxt(k,515)*y(k,270) &
                      + .040_r8*rxt(k,520)*y(k,271)
         mat(k,2404) = .500_r8*rxt(k,390)*y(k,16)
         mat(k,2088) = .520_r8*rxt(k,503)*y(k,6) + .520_r8*rxt(k,506)*y(k,111)
         mat(k,1163) = .500_r8*rxt(k,524)*y(k,265)
         mat(k,1088) = .250_r8*rxt(k,421)*y(k,125)
         mat(k,981) = .820_r8*rxt(k,373)*y(k,125) + .820_r8*rxt(k,371)*y(k,237)
         mat(k,2162) = .820_r8*rxt(k,371)*y(k,234) + .150_r8*rxt(k,513)*y(k,270) &
                      + .025_r8*rxt(k,518)*y(k,271)
         mat(k,610) = .500_r8*rxt(k,392)*y(k,125)
         mat(k,1823) = mat(k,1823) + .800_r8*rxt(k,419)*y(k,1) + .500_r8*rxt(k,524) &
                      *y(k,217)
         mat(k,1219) = .270_r8*rxt(k,515)*y(k,125) + .150_r8*rxt(k,513)*y(k,237)
         mat(k,1109) = .040_r8*rxt(k,520)*y(k,125) + .025_r8*rxt(k,518)*y(k,237)
         mat(k,1326) = -(rxt(k,376)*y(k,127) + rxt(k,377)*y(k,265))
         mat(k,2419) = -rxt(k,376)*y(k,49)
         mat(k,1838) = -rxt(k,377)*y(k,49)
         mat(k,1254) = rxt(k,379)*y(k,265)
         mat(k,1315) = .880_r8*rxt(k,401)*y(k,137)
         mat(k,1419) = .500_r8*rxt(k,415)*y(k,137)
         mat(k,1948) = .170_r8*rxt(k,474)*y(k,238) + .050_r8*rxt(k,437)*y(k,246) &
                      + .250_r8*rxt(k,399)*y(k,250) + .170_r8*rxt(k,480)*y(k,253) &
                      + .400_r8*rxt(k,490)*y(k,272) + .250_r8*rxt(k,456)*y(k,274) &
                      + .540_r8*rxt(k,496)*y(k,275) + .510_r8*rxt(k,499)*y(k,277)
         mat(k,2419) = mat(k,2419) + .050_r8*rxt(k,438)*y(k,246) + .250_r8*rxt(k,398) &
                      *y(k,250) + .250_r8*rxt(k,457)*y(k,274)
         mat(k,919) = rxt(k,380)*y(k,265)
         mat(k,2100) = .880_r8*rxt(k,401)*y(k,106) + .500_r8*rxt(k,415)*y(k,112)
         mat(k,1467) = .250_r8*rxt(k,395)*y(k,250) + .250_r8*rxt(k,453)*y(k,274)
         mat(k,2175) = .240_r8*rxt(k,396)*y(k,250) + .500_r8*rxt(k,384)*y(k,269) &
                      + .100_r8*rxt(k,454)*y(k,274)
         mat(k,825) = .170_r8*rxt(k,474)*y(k,125) + .070_r8*rxt(k,473)*y(k,243)
         mat(k,2358) = .070_r8*rxt(k,473)*y(k,238) + .070_r8*rxt(k,479)*y(k,253)
         mat(k,1376) = .050_r8*rxt(k,437)*y(k,125) + .050_r8*rxt(k,438)*y(k,127)
         mat(k,1400) = .250_r8*rxt(k,399)*y(k,125) + .250_r8*rxt(k,398)*y(k,127) &
                      + .250_r8*rxt(k,395)*y(k,236) + .240_r8*rxt(k,396)*y(k,237)
         mat(k,963) = .170_r8*rxt(k,480)*y(k,125) + .070_r8*rxt(k,479)*y(k,243)
         mat(k,1838) = mat(k,1838) + rxt(k,379)*y(k,95) + rxt(k,380)*y(k,128)
         mat(k,1263) = .500_r8*rxt(k,384)*y(k,237)
         mat(k,793) = .400_r8*rxt(k,490)*y(k,125)
         mat(k,1279) = .250_r8*rxt(k,456)*y(k,125) + .250_r8*rxt(k,457)*y(k,127) &
                      + .250_r8*rxt(k,453)*y(k,236) + .100_r8*rxt(k,454)*y(k,237)
         mat(k,817) = .540_r8*rxt(k,496)*y(k,125)
         mat(k,560) = .510_r8*rxt(k,499)*y(k,125)
         mat(k,750) = -(rxt(k,357)*y(k,265))
         mat(k,1794) = -rxt(k,357)*y(k,50)
         mat(k,1174) = .120_r8*rxt(k,370)*y(k,137)
         mat(k,2076) = .120_r8*rxt(k,370)*y(k,29)
         mat(k,1458) = .100_r8*rxt(k,354)*y(k,237) + .150_r8*rxt(k,355)*y(k,243)
         mat(k,2154) = .100_r8*rxt(k,354)*y(k,236)
         mat(k,2325) = .150_r8*rxt(k,355)*y(k,236) + .150_r8*rxt(k,406)*y(k,252)
         mat(k,1439) = .150_r8*rxt(k,406)*y(k,243)
         mat(k,665) = -(rxt(k,358)*y(k,265))
         mat(k,1786) = -rxt(k,358)*y(k,51)
         mat(k,1457) = .400_r8*rxt(k,355)*y(k,243)
         mat(k,2320) = .400_r8*rxt(k,355)*y(k,236) + .400_r8*rxt(k,406)*y(k,252)
         mat(k,1438) = .400_r8*rxt(k,406)*y(k,243)
         mat(k,869) = -(rxt(k,325)*y(k,265))
         mat(k,1804) = -rxt(k,325)*y(k,52)
         mat(k,1292) = .200_r8*rxt(k,443)*y(k,237)
         mat(k,945) = .300_r8*rxt(k,343)*y(k,237)
         mat(k,2155) = .200_r8*rxt(k,443)*y(k,101) + .300_r8*rxt(k,343)*y(k,233) &
                      + 2.000_r8*rxt(k,322)*y(k,237) + .250_r8*rxt(k,429)*y(k,245) &
                      + .250_r8*rxt(k,434)*y(k,246) + .250_r8*rxt(k,396)*y(k,250) &
                      + .250_r8*rxt(k,508)*y(k,260) + .500_r8*rxt(k,384)*y(k,269) &
                      + .250_r8*rxt(k,513)*y(k,270) + .250_r8*rxt(k,518)*y(k,271) &
                      + .300_r8*rxt(k,454)*y(k,274)
         mat(k,1336) = .250_r8*rxt(k,429)*y(k,237)
         mat(k,1365) = .250_r8*rxt(k,434)*y(k,237)
         mat(k,1394) = .250_r8*rxt(k,396)*y(k,237)
         mat(k,1149) = .250_r8*rxt(k,508)*y(k,237)
         mat(k,1260) = .500_r8*rxt(k,384)*y(k,237)
         mat(k,1218) = .250_r8*rxt(k,513)*y(k,237)
         mat(k,1106) = .250_r8*rxt(k,518)*y(k,237)
         mat(k,1273) = .300_r8*rxt(k,454)*y(k,237)
         mat(k,445) = -(rxt(k,326)*y(k,265))
         mat(k,1759) = -rxt(k,326)*y(k,53)
         mat(k,2151) = rxt(k,323)*y(k,243)
         mat(k,2301) = rxt(k,323)*y(k,237)
         mat(k,1634) = -(rxt(k,238)*y(k,56) + rxt(k,294)*y(k,73) + rxt(k,327)*y(k,265) &
                      + (rxt(k,333) + rxt(k,334) + rxt(k,335)) * y(k,261))
         mat(k,2005) = -rxt(k,238)*y(k,54)
         mat(k,972) = -rxt(k,294)*y(k,54)
         mat(k,1852) = -rxt(k,327)*y(k,54)
         mat(k,2228) = -(rxt(k,333) + rxt(k,334) + rxt(k,335)) * y(k,54)
         mat(k,1186) = .100_r8*rxt(k,370)*y(k,137)
         mat(k,2112) = .100_r8*rxt(k,370)*y(k,29)
         mat(k,451) = -(rxt(k,290)*y(k,261) + rxt(k,307)*y(k,56) + rxt(k,308)*y(k,265))
         mat(k,2218) = -rxt(k,290)*y(k,55)
         mat(k,1980) = -rxt(k,307)*y(k,55)
         mat(k,1760) = -rxt(k,308)*y(k,55)
         mat(k,2009) = -(rxt(k,237)*y(k,42) + rxt(k,238)*y(k,54) + rxt(k,239)*y(k,77) &
                      + rxt(k,240)*y(k,79) + (rxt(k,241) + rxt(k,242)) * y(k,243) &
                      + rxt(k,243)*y(k,137) + rxt(k,250)*y(k,60) + rxt(k,259)*y(k,92) &
                      + rxt(k,300)*y(k,41) + rxt(k,302)*y(k,43) + rxt(k,305)*y(k,46) &
                      + rxt(k,307)*y(k,55) + rxt(k,348)*y(k,28) + rxt(k,378)*y(k,31))
         mat(k,2461) = -rxt(k,237)*y(k,56)
         mat(k,1638) = -rxt(k,238)*y(k,56)
         mat(k,1510) = -rxt(k,239)*y(k,56)
         mat(k,646) = -rxt(k,240)*y(k,56)
         mat(k,2375) = -(rxt(k,241) + rxt(k,242)) * y(k,56)
         mat(k,2116) = -rxt(k,243)*y(k,56)
         mat(k,1048) = -rxt(k,250)*y(k,56)
         mat(k,885) = -rxt(k,259)*y(k,56)
         mat(k,520) = -rxt(k,300)*y(k,56)
         mat(k,685) = -rxt(k,302)*y(k,56)
         mat(k,415) = -rxt(k,305)*y(k,56)
         mat(k,455) = -rxt(k,307)*y(k,56)
         mat(k,331) = -rxt(k,348)*y(k,56)
         mat(k,337) = -rxt(k,378)*y(k,56)
         mat(k,1614) = rxt(k,278)*y(k,59)
         mat(k,147) = 4.000_r8*rxt(k,262)*y(k,261)
         mat(k,181) = rxt(k,263)*y(k,261)
         mat(k,157) = 2.000_r8*rxt(k,264)*y(k,261)
         mat(k,191) = 2.000_r8*rxt(k,265)*y(k,261)
         mat(k,161) = 2.000_r8*rxt(k,266)*y(k,261)
         mat(k,196) = rxt(k,267)*y(k,261)
         mat(k,165) = 2.000_r8*rxt(k,268)*y(k,261)
         mat(k,169) = 3.000_r8*rxt(k,304)*y(k,265)
         mat(k,415) = mat(k,415) + rxt(k,306)*y(k,265)
         mat(k,2488) = rxt(k,278)*y(k,19) + (4.000_r8*rxt(k,245)+2.000_r8*rxt(k,247)) &
                      *y(k,59) + rxt(k,249)*y(k,125) + rxt(k,254)*y(k,135) &
                      + rxt(k,533)*y(k,157) + rxt(k,244)*y(k,237) + rxt(k,255) &
                      *y(k,265)
         mat(k,296) = rxt(k,299)*y(k,261)
         mat(k,292) = rxt(k,314)*y(k,261) + rxt(k,309)*y(k,265)
         mat(k,302) = rxt(k,315)*y(k,261) + rxt(k,310)*y(k,265)
         mat(k,358) = rxt(k,316)*y(k,261) + rxt(k,311)*y(k,265)
         mat(k,1545) = rxt(k,257)*y(k,135) + rxt(k,269)*y(k,261) + rxt(k,258)*y(k,265)
         mat(k,1963) = rxt(k,249)*y(k,59)
         mat(k,2052) = rxt(k,254)*y(k,59) + rxt(k,257)*y(k,85)
         mat(k,1528) = rxt(k,533)*y(k,59)
         mat(k,2189) = rxt(k,244)*y(k,59)
         mat(k,2232) = 4.000_r8*rxt(k,262)*y(k,33) + rxt(k,263)*y(k,34) &
                      + 2.000_r8*rxt(k,264)*y(k,36) + 2.000_r8*rxt(k,265)*y(k,37) &
                      + 2.000_r8*rxt(k,266)*y(k,38) + rxt(k,267)*y(k,39) &
                      + 2.000_r8*rxt(k,268)*y(k,40) + rxt(k,299)*y(k,65) + rxt(k,314) &
                      *y(k,82) + rxt(k,315)*y(k,83) + rxt(k,316)*y(k,84) + rxt(k,269) &
                      *y(k,85)
         mat(k,1856) = 3.000_r8*rxt(k,304)*y(k,44) + rxt(k,306)*y(k,46) + rxt(k,255) &
                      *y(k,59) + rxt(k,309)*y(k,82) + rxt(k,310)*y(k,83) + rxt(k,311) &
                      *y(k,84) + rxt(k,258)*y(k,85)
         mat(k,1975) = rxt(k,250)*y(k,60)
         mat(k,2474) = 2.000_r8*rxt(k,246)*y(k,59)
         mat(k,1042) = rxt(k,250)*y(k,56) + (rxt(k,595)+rxt(k,600)+rxt(k,605))*y(k,85)
         mat(k,1535) = (rxt(k,595)+rxt(k,600)+rxt(k,605))*y(k,60) + (rxt(k,590) &
                       +rxt(k,596)+rxt(k,601))*y(k,92)
         mat(k,881) = (rxt(k,590)+rxt(k,596)+rxt(k,601))*y(k,85)
         mat(k,2473) = 2.000_r8*rxt(k,271)*y(k,59)
         mat(k,2498) = -(rxt(k,244)*y(k,237) + (4._r8*rxt(k,245) + 4._r8*rxt(k,246) &
                      + 4._r8*rxt(k,247) + 4._r8*rxt(k,271)) * y(k,59) + rxt(k,248) &
                      *y(k,243) + rxt(k,249)*y(k,125) + rxt(k,251)*y(k,126) + rxt(k,254) &
                      *y(k,135) + (rxt(k,255) + rxt(k,256)) * y(k,265) + (rxt(k,277) &
                      + rxt(k,278) + rxt(k,279)) * y(k,19) + rxt(k,533)*y(k,157))
         mat(k,2199) = -rxt(k,244)*y(k,59)
         mat(k,2385) = -rxt(k,248)*y(k,59)
         mat(k,1973) = -rxt(k,249)*y(k,59)
         mat(k,1693) = -rxt(k,251)*y(k,59)
         mat(k,2062) = -rxt(k,254)*y(k,59)
         mat(k,1866) = -(rxt(k,255) + rxt(k,256)) * y(k,59)
         mat(k,1623) = -(rxt(k,277) + rxt(k,278) + rxt(k,279)) * y(k,59)
         mat(k,1533) = -rxt(k,533)*y(k,59)
         mat(k,2019) = rxt(k,259)*y(k,92) + rxt(k,243)*y(k,137) + rxt(k,242)*y(k,243)
         mat(k,1052) = rxt(k,252)*y(k,135)
         mat(k,1551) = rxt(k,270)*y(k,261)
         mat(k,887) = rxt(k,259)*y(k,56) + rxt(k,260)*y(k,135) + rxt(k,261)*y(k,265)
         mat(k,2062) = mat(k,2062) + rxt(k,252)*y(k,60) + rxt(k,260)*y(k,92)
         mat(k,2126) = rxt(k,243)*y(k,56)
         mat(k,397) = rxt(k,538)*y(k,157)
         mat(k,1533) = mat(k,1533) + rxt(k,538)*y(k,139)
         mat(k,2385) = mat(k,2385) + rxt(k,242)*y(k,56)
         mat(k,2242) = rxt(k,270)*y(k,85)
         mat(k,1866) = mat(k,1866) + rxt(k,261)*y(k,92)
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
         mat(k,1044) = -(rxt(k,250)*y(k,56) + rxt(k,252)*y(k,135) + rxt(k,253) &
                      *y(k,265) + (rxt(k,595) + rxt(k,600) + rxt(k,605)) * y(k,85))
         mat(k,1990) = -rxt(k,250)*y(k,60)
         mat(k,2039) = -rxt(k,252)*y(k,60)
         mat(k,1817) = -rxt(k,253)*y(k,60)
         mat(k,1539) = -(rxt(k,595) + rxt(k,600) + rxt(k,605)) * y(k,60)
         mat(k,2479) = rxt(k,251)*y(k,126)
         mat(k,1666) = rxt(k,251)*y(k,59)
         mat(k,1205) = -(rxt(k,337)*y(k,265))
         mat(k,1829) = -rxt(k,337)*y(k,62)
         mat(k,1071) = .230_r8*rxt(k,503)*y(k,137)
         mat(k,1554) = rxt(k,273)*y(k,42)
         mat(k,324) = .350_r8*rxt(k,339)*y(k,265)
         mat(k,580) = .630_r8*rxt(k,341)*y(k,137)
         mat(k,1180) = .560_r8*rxt(k,370)*y(k,137)
         mat(k,2450) = rxt(k,273)*y(k,17) + rxt(k,237)*y(k,56) + rxt(k,318)*y(k,127) &
                      + rxt(k,319)*y(k,135) + rxt(k,320)*y(k,265)
         mat(k,412) = rxt(k,305)*y(k,56)
         mat(k,1325) = rxt(k,376)*y(k,127) + rxt(k,377)*y(k,265)
         mat(k,1994) = rxt(k,237)*y(k,42) + rxt(k,305)*y(k,46)
         mat(k,1494) = rxt(k,618)*y(k,266)
         mat(k,1100) = rxt(k,364)*y(k,265)
         mat(k,930) = .620_r8*rxt(k,448)*y(k,137)
         mat(k,1313) = .650_r8*rxt(k,401)*y(k,137)
         mat(k,1015) = .230_r8*rxt(k,506)*y(k,137)
         mat(k,1417) = .560_r8*rxt(k,415)*y(k,137)
         mat(k,1939) = .170_r8*rxt(k,474)*y(k,238) + .220_r8*rxt(k,399)*y(k,250) &
                      + .400_r8*rxt(k,477)*y(k,251) + .350_r8*rxt(k,480)*y(k,253) &
                      + .225_r8*rxt(k,515)*y(k,270) + .250_r8*rxt(k,456)*y(k,274)
         mat(k,2410) = rxt(k,318)*y(k,42) + rxt(k,376)*y(k,49) + .220_r8*rxt(k,398) &
                      *y(k,250) + .500_r8*rxt(k,457)*y(k,274)
         mat(k,2040) = rxt(k,319)*y(k,42) + rxt(k,527)*y(k,140)
         mat(k,2093) = .230_r8*rxt(k,503)*y(k,6) + .630_r8*rxt(k,341)*y(k,25) &
                      + .560_r8*rxt(k,370)*y(k,29) + .620_r8*rxt(k,448)*y(k,98) &
                      + .650_r8*rxt(k,401)*y(k,106) + .230_r8*rxt(k,506)*y(k,111) &
                      + .560_r8*rxt(k,415)*y(k,112)
         mat(k,422) = rxt(k,527)*y(k,135) + rxt(k,528)*y(k,265)
         mat(k,1165) = .700_r8*rxt(k,524)*y(k,265)
         mat(k,1461) = .220_r8*rxt(k,395)*y(k,250) + .250_r8*rxt(k,453)*y(k,274)
         mat(k,2167) = .110_r8*rxt(k,396)*y(k,250) + .125_r8*rxt(k,513)*y(k,270) &
                      + .200_r8*rxt(k,454)*y(k,274)
         mat(k,824) = .170_r8*rxt(k,474)*y(k,125) + .070_r8*rxt(k,473)*y(k,243)
         mat(k,2350) = .070_r8*rxt(k,473)*y(k,238) + .160_r8*rxt(k,476)*y(k,251) &
                      + .140_r8*rxt(k,479)*y(k,253)
         mat(k,1395) = .220_r8*rxt(k,399)*y(k,125) + .220_r8*rxt(k,398)*y(k,127) &
                      + .220_r8*rxt(k,395)*y(k,236) + .110_r8*rxt(k,396)*y(k,237)
         mat(k,779) = .400_r8*rxt(k,477)*y(k,125) + .160_r8*rxt(k,476)*y(k,243)
         mat(k,962) = .350_r8*rxt(k,480)*y(k,125) + .140_r8*rxt(k,479)*y(k,243)
         mat(k,1829) = mat(k,1829) + .350_r8*rxt(k,339)*y(k,24) + rxt(k,320)*y(k,42) &
                      + rxt(k,377)*y(k,49) + rxt(k,364)*y(k,75) + rxt(k,528)*y(k,140) &
                      + .700_r8*rxt(k,524)*y(k,217)
         mat(k,854) = rxt(k,618)*y(k,63)
         mat(k,1221) = .225_r8*rxt(k,515)*y(k,125) + .125_r8*rxt(k,513)*y(k,237)
         mat(k,1275) = .250_r8*rxt(k,456)*y(k,125) + .500_r8*rxt(k,457)*y(k,127) &
                      + .250_r8*rxt(k,453)*y(k,236) + .200_r8*rxt(k,454)*y(k,237)
         mat(k,1495) = -(rxt(k,618)*y(k,266))
         mat(k,855) = -rxt(k,618)*y(k,63)
         mat(k,1075) = .270_r8*rxt(k,503)*y(k,137)
         mat(k,1184) = .200_r8*rxt(k,370)*y(k,137)
         mat(k,751) = rxt(k,357)*y(k,265)
         mat(k,667) = .500_r8*rxt(k,358)*y(k,265)
         mat(k,1206) = rxt(k,337)*y(k,265)
         mat(k,1212) = .800_r8*rxt(k,363)*y(k,265)
         mat(k,1101) = rxt(k,364)*y(k,265)
         mat(k,956) = rxt(k,329)*y(k,265)
         mat(k,628) = .500_r8*rxt(k,414)*y(k,265)
         mat(k,1019) = .270_r8*rxt(k,506)*y(k,137)
         mat(k,1424) = .100_r8*rxt(k,415)*y(k,137)
         mat(k,1955) = rxt(k,356)*y(k,236) + .900_r8*rxt(k,515)*y(k,270)
         mat(k,2107) = .270_r8*rxt(k,503)*y(k,6) + .200_r8*rxt(k,370)*y(k,29) &
                      + .270_r8*rxt(k,506)*y(k,111) + .100_r8*rxt(k,415)*y(k,112)
         mat(k,1168) = 1.800_r8*rxt(k,524)*y(k,265)
         mat(k,1474) = rxt(k,356)*y(k,125) + 4.000_r8*rxt(k,353)*y(k,236) &
                      + .900_r8*rxt(k,354)*y(k,237) + rxt(k,428)*y(k,245) &
                      + 2.000_r8*rxt(k,404)*y(k,252) + rxt(k,453)*y(k,274)
         mat(k,2182) = .900_r8*rxt(k,354)*y(k,236) + rxt(k,405)*y(k,252) &
                      + .500_r8*rxt(k,513)*y(k,270)
         mat(k,2365) = .450_r8*rxt(k,406)*y(k,252)
         mat(k,1349) = rxt(k,428)*y(k,236)
         mat(k,1444) = 2.000_r8*rxt(k,404)*y(k,236) + rxt(k,405)*y(k,237) &
                      + .450_r8*rxt(k,406)*y(k,243) + 4.000_r8*rxt(k,407)*y(k,252)
         mat(k,1845) = rxt(k,357)*y(k,50) + .500_r8*rxt(k,358)*y(k,51) + rxt(k,337) &
                      *y(k,62) + .800_r8*rxt(k,363)*y(k,74) + rxt(k,364)*y(k,75) &
                      + rxt(k,329)*y(k,87) + .500_r8*rxt(k,414)*y(k,110) &
                      + 1.800_r8*rxt(k,524)*y(k,217)
         mat(k,1226) = .900_r8*rxt(k,515)*y(k,125) + .500_r8*rxt(k,513)*y(k,237)
         mat(k,1281) = rxt(k,453)*y(k,236)
         mat(k,286) = -(rxt(k,298)*y(k,261))
         mat(k,2212) = -rxt(k,298)*y(k,64)
         mat(k,179) = rxt(k,263)*y(k,261)
         mat(k,184) = rxt(k,289)*y(k,261)
         mat(k,189) = rxt(k,265)*y(k,261)
         mat(k,160) = 2.000_r8*rxt(k,266)*y(k,261)
         mat(k,194) = 2.000_r8*rxt(k,267)*y(k,261)
         mat(k,164) = rxt(k,268)*y(k,261)
         mat(k,172) = 2.000_r8*rxt(k,291)*y(k,261)
         mat(k,298) = rxt(k,315)*y(k,261) + rxt(k,310)*y(k,265)
         mat(k,354) = rxt(k,316)*y(k,261) + rxt(k,311)*y(k,265)
         mat(k,2212) = mat(k,2212) + rxt(k,263)*y(k,34) + rxt(k,289)*y(k,35) &
                      + rxt(k,265)*y(k,37) + 2.000_r8*rxt(k,266)*y(k,38) &
                      + 2.000_r8*rxt(k,267)*y(k,39) + rxt(k,268)*y(k,40) &
                      + 2.000_r8*rxt(k,291)*y(k,78) + rxt(k,315)*y(k,83) + rxt(k,316) &
                      *y(k,84)
         mat(k,1735) = rxt(k,310)*y(k,83) + rxt(k,311)*y(k,84)
         mat(k,294) = -(rxt(k,299)*y(k,261))
         mat(k,2214) = -rxt(k,299)*y(k,65)
         mat(k,156) = rxt(k,264)*y(k,261)
         mat(k,190) = rxt(k,265)*y(k,261)
         mat(k,290) = rxt(k,314)*y(k,261) + rxt(k,309)*y(k,265)
         mat(k,2214) = mat(k,2214) + rxt(k,264)*y(k,36) + rxt(k,265)*y(k,37) &
                      + rxt(k,314)*y(k,82)
         mat(k,1737) = rxt(k,309)*y(k,82)
         mat(k,238) = -(rxt(k,472)*y(k,265))
         mat(k,1725) = -rxt(k,472)*y(k,66)
         mat(k,232) = .180_r8*rxt(k,492)*y(k,265)
         mat(k,1725) = mat(k,1725) + .180_r8*rxt(k,492)*y(k,219)
         mat(k,345) = -(rxt(k,525)*y(k,127) + (rxt(k,526) + rxt(k,540)) * y(k,265))
         mat(k,2390) = -rxt(k,525)*y(k,67)
         mat(k,1745) = -(rxt(k,526) + rxt(k,540)) * y(k,67)
         mat(k,840) = rxt(k,359)*y(k,243)
         mat(k,2290) = rxt(k,359)*y(k,242)
         mat(k,970) = -(rxt(k,294)*y(k,54) + rxt(k,295)*y(k,77) + rxt(k,296)*y(k,278) &
                      + rxt(k,297)*y(k,89))
         mat(k,1626) = -rxt(k,294)*y(k,73)
         mat(k,1505) = -rxt(k,295)*y(k,73)
         mat(k,2503) = -rxt(k,296)*y(k,73)
         mat(k,2245) = -rxt(k,297)*y(k,73)
         mat(k,185) = rxt(k,289)*y(k,261)
         mat(k,195) = rxt(k,267)*y(k,261)
         mat(k,287) = 2.000_r8*rxt(k,298)*y(k,261)
         mat(k,295) = rxt(k,299)*y(k,261)
         mat(k,2222) = rxt(k,289)*y(k,35) + rxt(k,267)*y(k,39) + 2.000_r8*rxt(k,298) &
                      *y(k,64) + rxt(k,299)*y(k,65)
         mat(k,1211) = -(rxt(k,363)*y(k,265))
         mat(k,1830) = -rxt(k,363)*y(k,74)
         mat(k,635) = .700_r8*rxt(k,439)*y(k,265)
         mat(k,595) = .500_r8*rxt(k,440)*y(k,265)
         mat(k,465) = rxt(k,451)*y(k,265)
         mat(k,1940) = .050_r8*rxt(k,437)*y(k,246) + .530_r8*rxt(k,399)*y(k,250) &
                      + .225_r8*rxt(k,515)*y(k,270) + .250_r8*rxt(k,456)*y(k,274)
         mat(k,2411) = .050_r8*rxt(k,438)*y(k,246) + .530_r8*rxt(k,398)*y(k,250) &
                      + .250_r8*rxt(k,457)*y(k,274)
         mat(k,1583) = rxt(k,362)*y(k,241)
         mat(k,1462) = .530_r8*rxt(k,395)*y(k,250) + .250_r8*rxt(k,453)*y(k,274)
         mat(k,2168) = .260_r8*rxt(k,396)*y(k,250) + .125_r8*rxt(k,513)*y(k,270) &
                      + .100_r8*rxt(k,454)*y(k,274)
         mat(k,512) = rxt(k,362)*y(k,136)
         mat(k,1370) = .050_r8*rxt(k,437)*y(k,125) + .050_r8*rxt(k,438)*y(k,127)
         mat(k,1396) = .530_r8*rxt(k,399)*y(k,125) + .530_r8*rxt(k,398)*y(k,127) &
                      + .530_r8*rxt(k,395)*y(k,236) + .260_r8*rxt(k,396)*y(k,237)
         mat(k,1830) = mat(k,1830) + .700_r8*rxt(k,439)*y(k,99) + .500_r8*rxt(k,440) &
                      *y(k,100) + rxt(k,451)*y(k,116)
         mat(k,1222) = .225_r8*rxt(k,515)*y(k,125) + .125_r8*rxt(k,513)*y(k,237)
         mat(k,1276) = .250_r8*rxt(k,456)*y(k,125) + .250_r8*rxt(k,457)*y(k,127) &
                      + .250_r8*rxt(k,453)*y(k,236) + .100_r8*rxt(k,454)*y(k,237)
         mat(k,1099) = -(rxt(k,364)*y(k,265))
         mat(k,1820) = -rxt(k,364)*y(k,75)
         mat(k,323) = .650_r8*rxt(k,339)*y(k,265)
         mat(k,1209) = .200_r8*rxt(k,363)*y(k,265)
         mat(k,1134) = rxt(k,452)*y(k,265)
         mat(k,1931) = rxt(k,463)*y(k,230) + .050_r8*rxt(k,437)*y(k,246) &
                      + .400_r8*rxt(k,477)*y(k,251) + .170_r8*rxt(k,480)*y(k,253) &
                      + .700_r8*rxt(k,483)*y(k,267) + .600_r8*rxt(k,490)*y(k,272) &
                      + .250_r8*rxt(k,456)*y(k,274) + .340_r8*rxt(k,496)*y(k,275) &
                      + .170_r8*rxt(k,499)*y(k,277)
         mat(k,2401) = .050_r8*rxt(k,438)*y(k,246) + .250_r8*rxt(k,457)*y(k,274)
         mat(k,534) = rxt(k,463)*y(k,125)
         mat(k,1459) = .250_r8*rxt(k,453)*y(k,274)
         mat(k,2159) = .100_r8*rxt(k,454)*y(k,274)
         mat(k,2343) = .160_r8*rxt(k,476)*y(k,251) + .070_r8*rxt(k,479)*y(k,253)
         mat(k,1368) = .050_r8*rxt(k,437)*y(k,125) + .050_r8*rxt(k,438)*y(k,127)
         mat(k,778) = .400_r8*rxt(k,477)*y(k,125) + .160_r8*rxt(k,476)*y(k,243)
         mat(k,961) = .170_r8*rxt(k,480)*y(k,125) + .070_r8*rxt(k,479)*y(k,243)
         mat(k,1820) = mat(k,1820) + .650_r8*rxt(k,339)*y(k,24) + .200_r8*rxt(k,363) &
                      *y(k,74) + rxt(k,452)*y(k,117)
         mat(k,503) = .700_r8*rxt(k,483)*y(k,125)
         mat(k,791) = .600_r8*rxt(k,490)*y(k,125)
         mat(k,1274) = .250_r8*rxt(k,456)*y(k,125) + .250_r8*rxt(k,457)*y(k,127) &
                      + .250_r8*rxt(k,453)*y(k,236) + .100_r8*rxt(k,454)*y(k,237)
         mat(k,815) = .340_r8*rxt(k,496)*y(k,125)
         mat(k,559) = .170_r8*rxt(k,499)*y(k,125)
         mat(k,2140) = -((rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,243) + rxt(k,195) &
                      *y(k,136) + rxt(k,198)*y(k,137))
         mat(k,2378) = -(rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,76)
         mat(k,1595) = -rxt(k,195)*y(k,76)
         mat(k,2119) = -rxt(k,198)*y(k,76)
         mat(k,2464) = rxt(k,320)*y(k,265)
         mat(k,1641) = rxt(k,334)*y(k,261)
         mat(k,2012) = rxt(k,239)*y(k,77)
         mat(k,974) = rxt(k,295)*y(k,77)
         mat(k,1512) = rxt(k,239)*y(k,56) + rxt(k,295)*y(k,73) + rxt(k,190)*y(k,135) &
                      + rxt(k,173)*y(k,261) + rxt(k,199)*y(k,265)
         mat(k,878) = rxt(k,293)*y(k,261)
         mat(k,1547) = rxt(k,270)*y(k,261)
         mat(k,1040) = rxt(k,222)*y(k,265)
         mat(k,2055) = rxt(k,190)*y(k,77) + rxt(k,202)*y(k,265)
         mat(k,426) = rxt(k,528)*y(k,265)
         mat(k,761) = rxt(k,534)*y(k,265)
         mat(k,1531) = rxt(k,539)*y(k,265)
         mat(k,2235) = rxt(k,334)*y(k,54) + rxt(k,173)*y(k,77) + rxt(k,293)*y(k,81) &
                      + rxt(k,270)*y(k,85)
         mat(k,1859) = rxt(k,320)*y(k,42) + rxt(k,199)*y(k,77) + rxt(k,222)*y(k,113) &
                      + rxt(k,202)*y(k,135) + rxt(k,528)*y(k,140) + rxt(k,534) &
                      *y(k,155) + rxt(k,539)*y(k,157)
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
         mat(k,1506) = -(rxt(k,173)*y(k,261) + rxt(k,190)*y(k,135) + rxt(k,199) &
                      *y(k,265) + rxt(k,239)*y(k,56) + rxt(k,295)*y(k,73))
         mat(k,2223) = -rxt(k,173)*y(k,77)
         mat(k,2042) = -rxt(k,190)*y(k,77)
         mat(k,1846) = -rxt(k,199)*y(k,77)
         mat(k,2000) = -rxt(k,239)*y(k,77)
         mat(k,971) = -rxt(k,295)*y(k,77)
         mat(k,1629) = rxt(k,335)*y(k,261)
         mat(k,2128) = rxt(k,192)*y(k,243)
         mat(k,2366) = rxt(k,192)*y(k,76)
         mat(k,2223) = mat(k,2223) + rxt(k,335)*y(k,54)
         mat(k,171) = -(rxt(k,291)*y(k,261))
         mat(k,2205) = -rxt(k,291)*y(k,78)
         mat(k,643) = -(rxt(k,191)*y(k,135) + rxt(k,200)*y(k,265) + rxt(k,240)*y(k,56))
         mat(k,2026) = -rxt(k,191)*y(k,79)
         mat(k,1783) = -rxt(k,200)*y(k,79)
         mat(k,1983) = -rxt(k,240)*y(k,79)
         mat(k,2317) = 2.000_r8*rxt(k,206)*y(k,243)
         mat(k,1783) = mat(k,1783) + 2.000_r8*rxt(k,205)*y(k,265)
         mat(k,308) = rxt(k,541)*y(k,278)
         mat(k,2500) = rxt(k,541)*y(k,159)
         mat(k,873) = -(rxt(k,286)*y(k,135) + rxt(k,287)*y(k,265) + (rxt(k,292) &
                      + rxt(k,293)) * y(k,261))
         mat(k,2032) = -rxt(k,286)*y(k,81)
         mat(k,1805) = -rxt(k,287)*y(k,81)
         mat(k,2221) = -(rxt(k,292) + rxt(k,293)) * y(k,81)
         mat(k,1553) = rxt(k,273)*y(k,42) + rxt(k,274)*y(k,243)
         mat(k,2448) = rxt(k,273)*y(k,17)
         mat(k,2335) = rxt(k,274)*y(k,17)
         mat(k,289) = -(rxt(k,309)*y(k,265) + rxt(k,314)*y(k,261))
         mat(k,1736) = -rxt(k,309)*y(k,82)
         mat(k,2213) = -rxt(k,314)*y(k,82)
         mat(k,299) = -(rxt(k,310)*y(k,265) + rxt(k,315)*y(k,261))
         mat(k,1738) = -rxt(k,310)*y(k,83)
         mat(k,2215) = -rxt(k,315)*y(k,83)
         mat(k,355) = -(rxt(k,311)*y(k,265) + rxt(k,316)*y(k,261))
         mat(k,1746) = -rxt(k,311)*y(k,84)
         mat(k,2217) = -rxt(k,316)*y(k,84)
         mat(k,1540) = -(rxt(k,257)*y(k,135) + rxt(k,258)*y(k,265) + (rxt(k,269) &
                      + rxt(k,270)) * y(k,261) + (rxt(k,590) + rxt(k,596) + rxt(k,601) &
                      ) * y(k,92) + (rxt(k,595) + rxt(k,600) + rxt(k,605)) * y(k,60) &
                      + (rxt(k,597) + rxt(k,602)) * y(k,91))
         mat(k,2044) = -rxt(k,257)*y(k,85)
         mat(k,1848) = -rxt(k,258)*y(k,85)
         mat(k,2224) = -(rxt(k,269) + rxt(k,270)) * y(k,85)
         mat(k,883) = -(rxt(k,590) + rxt(k,596) + rxt(k,601)) * y(k,85)
         mat(k,1045) = -(rxt(k,595) + rxt(k,600) + rxt(k,605)) * y(k,85)
         mat(k,833) = -(rxt(k,597) + rxt(k,602)) * y(k,85)
         mat(k,329) = rxt(k,348)*y(k,56)
         mat(k,335) = rxt(k,378)*y(k,56)
         mat(k,517) = rxt(k,300)*y(k,56)
         mat(k,2453) = rxt(k,237)*y(k,56)
         mat(k,682) = rxt(k,302)*y(k,56)
         mat(k,413) = 2.000_r8*rxt(k,305)*y(k,56)
         mat(k,1630) = rxt(k,238)*y(k,56)
         mat(k,452) = rxt(k,307)*y(k,56)
         mat(k,2001) = rxt(k,348)*y(k,28) + rxt(k,378)*y(k,31) + rxt(k,300)*y(k,41) &
                      + rxt(k,237)*y(k,42) + rxt(k,302)*y(k,43) + 2.000_r8*rxt(k,305) &
                      *y(k,46) + rxt(k,238)*y(k,54) + rxt(k,307)*y(k,55) + rxt(k,239) &
                      *y(k,77) + rxt(k,240)*y(k,79) + rxt(k,259)*y(k,92) + rxt(k,241) &
                      *y(k,243)
         mat(k,2481) = rxt(k,256)*y(k,265)
         mat(k,1507) = rxt(k,239)*y(k,56)
         mat(k,644) = rxt(k,240)*y(k,56)
         mat(k,883) = mat(k,883) + rxt(k,259)*y(k,56)
         mat(k,2367) = rxt(k,241)*y(k,56)
         mat(k,1848) = mat(k,1848) + rxt(k,256)*y(k,59)
         mat(k,250) = -(rxt(k,328)*y(k,265) + rxt(k,336)*y(k,261))
         mat(k,1728) = -rxt(k,328)*y(k,86)
         mat(k,2211) = -rxt(k,336)*y(k,86)
         mat(k,955) = -(rxt(k,329)*y(k,265))
         mat(k,1810) = -rxt(k,329)*y(k,87)
         mat(k,1062) = .050_r8*rxt(k,503)*y(k,137)
         mat(k,322) = .350_r8*rxt(k,339)*y(k,265)
         mat(k,579) = .370_r8*rxt(k,341)*y(k,137)
         mat(k,1177) = .120_r8*rxt(k,370)*y(k,137)
         mat(k,928) = .110_r8*rxt(k,448)*y(k,137)
         mat(k,1312) = .330_r8*rxt(k,401)*y(k,137)
         mat(k,1006) = .050_r8*rxt(k,506)*y(k,137)
         mat(k,1414) = .120_r8*rxt(k,415)*y(k,137)
         mat(k,1926) = rxt(k,332)*y(k,244)
         mat(k,2080) = .050_r8*rxt(k,503)*y(k,6) + .370_r8*rxt(k,341)*y(k,25) &
                      + .120_r8*rxt(k,370)*y(k,29) + .110_r8*rxt(k,448)*y(k,98) &
                      + .330_r8*rxt(k,401)*y(k,106) + .050_r8*rxt(k,506)*y(k,111) &
                      + .120_r8*rxt(k,415)*y(k,112)
         mat(k,2339) = rxt(k,330)*y(k,244)
         mat(k,484) = rxt(k,332)*y(k,125) + rxt(k,330)*y(k,243)
         mat(k,1810) = mat(k,1810) + .350_r8*rxt(k,339)*y(k,24)
         mat(k,1625) = rxt(k,294)*y(k,73)
         mat(k,969) = rxt(k,294)*y(k,54) + rxt(k,295)*y(k,77) + rxt(k,297)*y(k,89) &
                      + rxt(k,296)*y(k,278)
         mat(k,1504) = rxt(k,295)*y(k,73)
         mat(k,2244) = rxt(k,297)*y(k,73)
         mat(k,2502) = rxt(k,296)*y(k,73)
         mat(k,2261) = -(rxt(k,231)*y(k,265) + rxt(k,297)*y(k,73))
         mat(k,1862) = -rxt(k,231)*y(k,89)
         mat(k,976) = -rxt(k,297)*y(k,89)
         mat(k,2467) = rxt(k,318)*y(k,127)
         mat(k,1201) = rxt(k,350)*y(k,127)
         mat(k,1330) = rxt(k,376)*y(k,127)
         mat(k,1050) = (rxt(k,595)+rxt(k,600)+rxt(k,605))*y(k,85)
         mat(k,348) = rxt(k,525)*y(k,127)
         mat(k,1549) = (rxt(k,595)+rxt(k,600)+rxt(k,605))*y(k,60)
         mat(k,1689) = rxt(k,230)*y(k,265)
         mat(k,2441) = rxt(k,318)*y(k,42) + rxt(k,350)*y(k,45) + rxt(k,376)*y(k,49) &
                      + rxt(k,525)*y(k,67)
         mat(k,1862) = mat(k,1862) + rxt(k,230)*y(k,126)
         mat(k,550) = -(rxt(k,208)*y(k,265))
         mat(k,1772) = -rxt(k,208)*y(k,90)
         mat(k,1654) = rxt(k,228)*y(k,243)
         mat(k,2313) = rxt(k,228)*y(k,126)
         mat(k,832) = -(rxt(k,288)*y(k,135) + (rxt(k,597) + rxt(k,602)) * y(k,85))
         mat(k,2030) = -rxt(k,288)*y(k,91)
         mat(k,1537) = -(rxt(k,597) + rxt(k,602)) * y(k,91)
         mat(k,1604) = rxt(k,280)*y(k,243)
         mat(k,2332) = rxt(k,280)*y(k,19)
         mat(k,882) = -(rxt(k,259)*y(k,56) + rxt(k,260)*y(k,135) + rxt(k,261)*y(k,265) &
                      + (rxt(k,590) + rxt(k,596) + rxt(k,601)) * y(k,85))
         mat(k,1986) = -rxt(k,259)*y(k,92)
         mat(k,2033) = -rxt(k,260)*y(k,92)
         mat(k,1806) = -rxt(k,261)*y(k,92)
         mat(k,1538) = -(rxt(k,590) + rxt(k,596) + rxt(k,601)) * y(k,92)
         mat(k,2477) = rxt(k,248)*y(k,243)
         mat(k,1043) = rxt(k,253)*y(k,265)
         mat(k,2336) = rxt(k,248)*y(k,59)
         mat(k,1806) = mat(k,1806) + rxt(k,253)*y(k,60)
         mat(k,1240) = -(rxt(k,394)*y(k,265))
         mat(k,1832) = -rxt(k,394)*y(k,93)
         mat(k,636) = .300_r8*rxt(k,439)*y(k,265)
         mat(k,596) = .500_r8*rxt(k,440)*y(k,265)
         mat(k,1942) = rxt(k,393)*y(k,240) + rxt(k,400)*y(k,250)
         mat(k,612) = rxt(k,393)*y(k,125)
         mat(k,1397) = rxt(k,400)*y(k,125)
         mat(k,1832) = mat(k,1832) + .300_r8*rxt(k,439)*y(k,99) + .500_r8*rxt(k,440) &
                      *y(k,100)
         mat(k,281) = -(rxt(k,425)*y(k,265))
         mat(k,1734) = -rxt(k,425)*y(k,94)
         mat(k,1253) = -(rxt(k,379)*y(k,265))
         mat(k,1833) = -rxt(k,379)*y(k,95)
         mat(k,637) = .700_r8*rxt(k,439)*y(k,265)
         mat(k,597) = .500_r8*rxt(k,440)*y(k,265)
         mat(k,626) = .500_r8*rxt(k,414)*y(k,265)
         mat(k,1943) = .050_r8*rxt(k,437)*y(k,246) + .220_r8*rxt(k,399)*y(k,250) &
                      + .250_r8*rxt(k,456)*y(k,274)
         mat(k,2414) = .050_r8*rxt(k,438)*y(k,246) + .220_r8*rxt(k,398)*y(k,250) &
                      + .250_r8*rxt(k,457)*y(k,274)
         mat(k,604) = .500_r8*rxt(k,383)*y(k,265)
         mat(k,1463) = .220_r8*rxt(k,395)*y(k,250) + .250_r8*rxt(k,453)*y(k,274)
         mat(k,2170) = .230_r8*rxt(k,396)*y(k,250) + .200_r8*rxt(k,384)*y(k,269) &
                      + .100_r8*rxt(k,454)*y(k,274)
         mat(k,1372) = .050_r8*rxt(k,437)*y(k,125) + .050_r8*rxt(k,438)*y(k,127)
         mat(k,1398) = .220_r8*rxt(k,399)*y(k,125) + .220_r8*rxt(k,398)*y(k,127) &
                      + .220_r8*rxt(k,395)*y(k,236) + .230_r8*rxt(k,396)*y(k,237)
         mat(k,1833) = mat(k,1833) + .700_r8*rxt(k,439)*y(k,99) + .500_r8*rxt(k,440) &
                      *y(k,100) + .500_r8*rxt(k,414)*y(k,110) + .500_r8*rxt(k,383) &
                      *y(k,153)
         mat(k,1261) = .200_r8*rxt(k,384)*y(k,237)
         mat(k,1277) = .250_r8*rxt(k,456)*y(k,125) + .250_r8*rxt(k,457)*y(k,127) &
                      + .250_r8*rxt(k,453)*y(k,236) + .100_r8*rxt(k,454)*y(k,237)
         mat(k,380) = -(rxt(k,426)*y(k,265))
         mat(k,1750) = -rxt(k,426)*y(k,96)
         mat(k,1894) = .870_r8*rxt(k,437)*y(k,246)
         mat(k,2391) = .950_r8*rxt(k,438)*y(k,246)
         mat(k,1455) = rxt(k,433)*y(k,246)
         mat(k,2150) = .750_r8*rxt(k,434)*y(k,246)
         mat(k,1361) = .870_r8*rxt(k,437)*y(k,125) + .950_r8*rxt(k,438)*y(k,127) &
                      + rxt(k,433)*y(k,236) + .750_r8*rxt(k,434)*y(k,237)
         mat(k,198) = -(rxt(k,427)*y(k,265))
         mat(k,1720) = -rxt(k,427)*y(k,97)
         mat(k,798) = .600_r8*rxt(k,450)*y(k,265)
         mat(k,1720) = mat(k,1720) + .600_r8*rxt(k,450)*y(k,103)
         mat(k,927) = -(rxt(k,441)*y(k,127) + rxt(k,448)*y(k,137) + rxt(k,449) &
                      *y(k,265))
         mat(k,2395) = -rxt(k,441)*y(k,98)
         mat(k,2079) = -rxt(k,448)*y(k,98)
         mat(k,1808) = -rxt(k,449)*y(k,98)
         mat(k,634) = -(rxt(k,439)*y(k,265))
         mat(k,1782) = -rxt(k,439)*y(k,99)
         mat(k,1908) = .080_r8*rxt(k,431)*y(k,245)
         mat(k,1334) = .080_r8*rxt(k,431)*y(k,125)
         mat(k,593) = -(rxt(k,440)*y(k,265))
         mat(k,1778) = -rxt(k,440)*y(k,100)
         mat(k,1906) = .080_r8*rxt(k,437)*y(k,246)
         mat(k,1362) = .080_r8*rxt(k,437)*y(k,125)
         mat(k,1298) = -(rxt(k,442)*y(k,236) + rxt(k,443)*y(k,237) + rxt(k,444) &
                      *y(k,243) + rxt(k,445)*y(k,125) + rxt(k,446)*y(k,127))
         mat(k,1465) = -rxt(k,442)*y(k,101)
         mat(k,2173) = -rxt(k,443)*y(k,101)
         mat(k,2356) = -rxt(k,444)*y(k,101)
         mat(k,1946) = -rxt(k,445)*y(k,101)
         mat(k,2417) = -rxt(k,446)*y(k,101)
         mat(k,931) = rxt(k,441)*y(k,127)
         mat(k,2417) = mat(k,2417) + rxt(k,441)*y(k,98)
         mat(k,439) = -(rxt(k,447)*y(k,265))
         mat(k,1758) = -rxt(k,447)*y(k,102)
         mat(k,1289) = rxt(k,444)*y(k,243)
         mat(k,2300) = rxt(k,444)*y(k,101)
         mat(k,799) = -(rxt(k,450)*y(k,265))
         mat(k,1799) = -rxt(k,450)*y(k,103)
         mat(k,2329) = rxt(k,430)*y(k,245) + rxt(k,435)*y(k,246)
         mat(k,1335) = rxt(k,430)*y(k,243)
         mat(k,1364) = rxt(k,435)*y(k,243)
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
         mat(k,116) = -(rxt(k,571)*y(k,265))
         mat(k,1711) = -rxt(k,571)*y(k,104)
         mat(k,124) = -(rxt(k,574)*y(k,265))
         mat(k,1713) = -rxt(k,574)*y(k,105)
         mat(k,1314) = -(rxt(k,401)*y(k,137) + rxt(k,402)*y(k,265))
         mat(k,2099) = -rxt(k,401)*y(k,106)
         mat(k,1837) = -rxt(k,402)*y(k,106)
         mat(k,932) = .300_r8*rxt(k,448)*y(k,137)
         mat(k,1947) = .360_r8*rxt(k,431)*y(k,245)
         mat(k,2418) = .400_r8*rxt(k,432)*y(k,245)
         mat(k,2099) = mat(k,2099) + .300_r8*rxt(k,448)*y(k,98)
         mat(k,1466) = .390_r8*rxt(k,428)*y(k,245)
         mat(k,2174) = .310_r8*rxt(k,429)*y(k,245)
         mat(k,1342) = .360_r8*rxt(k,431)*y(k,125) + .400_r8*rxt(k,432)*y(k,127) &
                      + .390_r8*rxt(k,428)*y(k,236) + .310_r8*rxt(k,429)*y(k,237)
         mat(k,388) = -(rxt(k,403)*y(k,265))
         mat(k,1752) = -rxt(k,403)*y(k,107)
         mat(k,2297) = rxt(k,397)*y(k,250)
         mat(k,1393) = rxt(k,397)*y(k,243)
         mat(k,565) = -(rxt(k,412)*y(k,265))
         mat(k,1774) = -rxt(k,412)*y(k,108)
         mat(k,1904) = .800_r8*rxt(k,421)*y(k,228)
         mat(k,1082) = .800_r8*rxt(k,421)*y(k,125)
         mat(k,383) = -(rxt(k,413)*y(k,265))
         mat(k,1751) = -rxt(k,413)*y(k,109)
         mat(k,2296) = .800_r8*rxt(k,410)*y(k,254)
         mat(k,742) = .800_r8*rxt(k,410)*y(k,243)
         mat(k,625) = -(rxt(k,414)*y(k,265))
         mat(k,1781) = -rxt(k,414)*y(k,110)
         mat(k,1658) = rxt(k,417)*y(k,252)
         mat(k,1437) = rxt(k,417)*y(k,126)
         mat(k,1007) = -(rxt(k,505)*y(k,127) + rxt(k,506)*y(k,137) + rxt(k,507) &
                      *y(k,265))
         mat(k,2398) = -rxt(k,505)*y(k,111)
         mat(k,2082) = -rxt(k,506)*y(k,111)
         mat(k,1815) = -rxt(k,507)*y(k,111)
         mat(k,1421) = -(rxt(k,415)*y(k,137) + rxt(k,416)*y(k,265))
         mat(k,2104) = -rxt(k,415)*y(k,112)
         mat(k,1842) = -rxt(k,416)*y(k,112)
         mat(k,935) = .200_r8*rxt(k,448)*y(k,137)
         mat(k,1952) = .560_r8*rxt(k,431)*y(k,245)
         mat(k,2423) = .600_r8*rxt(k,432)*y(k,245)
         mat(k,2104) = mat(k,2104) + .200_r8*rxt(k,448)*y(k,98)
         mat(k,1471) = .610_r8*rxt(k,428)*y(k,245)
         mat(k,2179) = .440_r8*rxt(k,429)*y(k,245)
         mat(k,1346) = .560_r8*rxt(k,431)*y(k,125) + .600_r8*rxt(k,432)*y(k,127) &
                      + .610_r8*rxt(k,428)*y(k,236) + .440_r8*rxt(k,429)*y(k,237)
         mat(k,1032) = -(rxt(k,211)*y(k,125) + (rxt(k,212) + rxt(k,213) + rxt(k,214) &
                      ) * y(k,126) + rxt(k,222)*y(k,265) + rxt(k,236)*y(k,136) &
                      + rxt(k,615)*y(k,264))
         mat(k,1929) = -rxt(k,211)*y(k,113)
         mat(k,1665) = -(rxt(k,212) + rxt(k,213) + rxt(k,214)) * y(k,113)
         mat(k,1816) = -rxt(k,222)*y(k,113)
         mat(k,1581) = -rxt(k,236)*y(k,113)
         mat(k,900) = -rxt(k,615)*y(k,113)
         mat(k,2038) = rxt(k,210)*y(k,256) + rxt(k,612)*y(k,259)
         mat(k,1581) = mat(k,1581) + rxt(k,613)*y(k,259)
         mat(k,911) = rxt(k,233)*y(k,256) + 1.100_r8*rxt(k,608)*y(k,257) &
                      + .200_r8*rxt(k,606)*y(k,258)
         mat(k,705) = rxt(k,210)*y(k,135) + rxt(k,233)*y(k,239)
         mat(k,678) = 1.100_r8*rxt(k,608)*y(k,239)
         mat(k,892) = .200_r8*rxt(k,606)*y(k,239)
         mat(k,547) = rxt(k,612)*y(k,135) + rxt(k,613)*y(k,136)
         mat(k,317) = -((rxt(k,226) + rxt(k,227)) * y(k,261))
         mat(k,2216) = -(rxt(k,226) + rxt(k,227)) * y(k,114)
         mat(k,1026) = rxt(k,212)*y(k,126)
         mat(k,1651) = rxt(k,212)*y(k,113)
         mat(k,1652) = rxt(k,229)*y(k,127)
         mat(k,2389) = rxt(k,229)*y(k,126)
         mat(k,463) = -(rxt(k,451)*y(k,265))
         mat(k,1762) = -rxt(k,451)*y(k,116)
         mat(k,1290) = .200_r8*rxt(k,443)*y(k,237)
         mat(k,2152) = .200_r8*rxt(k,443)*y(k,101)
         mat(k,1135) = -(rxt(k,452)*y(k,265))
         mat(k,1824) = -rxt(k,452)*y(k,117)
         mat(k,1294) = rxt(k,445)*y(k,125) + rxt(k,446)*y(k,127) + rxt(k,442)*y(k,236) &
                      + .800_r8*rxt(k,443)*y(k,237)
         mat(k,1935) = rxt(k,445)*y(k,101)
         mat(k,2405) = rxt(k,446)*y(k,101)
         mat(k,1460) = rxt(k,442)*y(k,101)
         mat(k,2163) = .800_r8*rxt(k,443)*y(k,101)
         mat(k,143) = -(rxt(k,542)*y(k,265))
         mat(k,1717) = -rxt(k,542)*y(k,121)
         mat(k,1962) = -(rxt(k,209)*y(k,256) + rxt(k,211)*y(k,113) + rxt(k,219) &
                      *y(k,127) + rxt(k,223)*y(k,243) + rxt(k,224)*y(k,137) + rxt(k,225) &
                      *y(k,135) + rxt(k,249)*y(k,59) + rxt(k,281)*y(k,19) + rxt(k,324) &
                      *y(k,237) + rxt(k,332)*y(k,244) + rxt(k,345)*y(k,233) + rxt(k,356) &
                      *y(k,236) + rxt(k,360)*y(k,242) + rxt(k,373)*y(k,234) + rxt(k,382) &
                      *y(k,268) + rxt(k,386)*y(k,269) + (rxt(k,392) + rxt(k,393) &
                      ) * y(k,240) + (rxt(k,399) + rxt(k,400)) * y(k,250) + rxt(k,408) &
                      *y(k,252) + rxt(k,411)*y(k,254) + (rxt(k,421) + rxt(k,422) &
                      ) * y(k,228) + rxt(k,431)*y(k,245) + rxt(k,437)*y(k,246) &
                      + rxt(k,445)*y(k,101) + rxt(k,456)*y(k,274) + rxt(k,460) &
                      *y(k,227) + rxt(k,463)*y(k,230) + rxt(k,468)*y(k,232) + rxt(k,470) &
                      *y(k,235) + rxt(k,474)*y(k,238) + rxt(k,477)*y(k,251) + rxt(k,480) &
                      *y(k,253) + rxt(k,483)*y(k,267) + rxt(k,490)*y(k,272) + rxt(k,496) &
                      *y(k,275) + rxt(k,499)*y(k,277) + rxt(k,510)*y(k,260) + rxt(k,515) &
                      *y(k,270) + rxt(k,520)*y(k,271) + rxt(k,617)*y(k,264))
         mat(k,707) = -rxt(k,209)*y(k,125)
         mat(k,1038) = -rxt(k,211)*y(k,125)
         mat(k,2434) = -rxt(k,219)*y(k,125)
         mat(k,2374) = -rxt(k,223)*y(k,125)
         mat(k,2115) = -rxt(k,224)*y(k,125)
         mat(k,2051) = -rxt(k,225)*y(k,125)
         mat(k,2487) = -rxt(k,249)*y(k,125)
         mat(k,1613) = -rxt(k,281)*y(k,125)
         mat(k,2188) = -rxt(k,324)*y(k,125)
         mat(k,486) = -rxt(k,332)*y(k,125)
         mat(k,951) = -rxt(k,345)*y(k,125)
         mat(k,1479) = -rxt(k,356)*y(k,125)
         mat(k,847) = -rxt(k,360)*y(k,125)
         mat(k,986) = -rxt(k,373)*y(k,125)
         mat(k,865) = -rxt(k,382)*y(k,125)
         mat(k,1267) = -rxt(k,386)*y(k,125)
         mat(k,614) = -(rxt(k,392) + rxt(k,393)) * y(k,125)
         mat(k,1407) = -(rxt(k,399) + rxt(k,400)) * y(k,125)
         mat(k,1447) = -rxt(k,408)*y(k,125)
         mat(k,748) = -rxt(k,411)*y(k,125)
         mat(k,1094) = -(rxt(k,421) + rxt(k,422)) * y(k,125)
         mat(k,1352) = -rxt(k,431)*y(k,125)
         mat(k,1385) = -rxt(k,437)*y(k,125)
         mat(k,1306) = -rxt(k,445)*y(k,125)
         mat(k,1284) = -rxt(k,456)*y(k,125)
         mat(k,575) = -rxt(k,460)*y(k,125)
         mat(k,537) = -rxt(k,463)*y(k,125)
         mat(k,481) = -rxt(k,468)*y(k,125)
         mat(k,701) = -rxt(k,470)*y(k,125)
         mat(k,828) = -rxt(k,474)*y(k,125)
         mat(k,781) = -rxt(k,477)*y(k,125)
         mat(k,966) = -rxt(k,480)*y(k,125)
         mat(k,506) = -rxt(k,483)*y(k,125)
         mat(k,796) = -rxt(k,490)*y(k,125)
         mat(k,821) = -rxt(k,496)*y(k,125)
         mat(k,563) = -rxt(k,499)*y(k,125)
         mat(k,1157) = -rxt(k,510)*y(k,125)
         mat(k,1229) = -rxt(k,515)*y(k,125)
         mat(k,1112) = -rxt(k,520)*y(k,125)
         mat(k,902) = -rxt(k,617)*y(k,125)
         mat(k,1038) = mat(k,1038) + 2.000_r8*rxt(k,213)*y(k,126) + rxt(k,236) &
                      *y(k,136) + rxt(k,222)*y(k,265)
         mat(k,319) = 2.000_r8*rxt(k,226)*y(k,261)
         mat(k,1682) = 2.000_r8*rxt(k,213)*y(k,113) + rxt(k,215)*y(k,135) + rxt(k,535) &
                      *y(k,157)
         mat(k,2051) = mat(k,2051) + rxt(k,215)*y(k,126)
         mat(k,1591) = rxt(k,236)*y(k,113) + rxt(k,234)*y(k,256)
         mat(k,1527) = rxt(k,535)*y(k,126)
         mat(k,707) = mat(k,707) + rxt(k,234)*y(k,136)
         mat(k,2231) = 2.000_r8*rxt(k,226)*y(k,114)
         mat(k,1855) = rxt(k,222)*y(k,113)
         mat(k,1680) = -((rxt(k,212) + rxt(k,213) + rxt(k,214)) * y(k,113) + (rxt(k,215) &
                      + rxt(k,217)) * y(k,135) + rxt(k,216)*y(k,137) + rxt(k,228) &
                      *y(k,243) + rxt(k,229)*y(k,127) + rxt(k,230)*y(k,265) + rxt(k,251) &
                      *y(k,59) + rxt(k,282)*y(k,19) + rxt(k,367)*y(k,236) + rxt(k,417) &
                      *y(k,252) + rxt(k,475)*y(k,238) + rxt(k,478)*y(k,251) + rxt(k,481) &
                      *y(k,253) + rxt(k,485)*y(k,144) + rxt(k,488)*y(k,227) + rxt(k,535) &
                      *y(k,157))
         mat(k,1036) = -(rxt(k,212) + rxt(k,213) + rxt(k,214)) * y(k,126)
         mat(k,2049) = -(rxt(k,215) + rxt(k,217)) * y(k,126)
         mat(k,2113) = -rxt(k,216)*y(k,126)
         mat(k,2372) = -rxt(k,228)*y(k,126)
         mat(k,2432) = -rxt(k,229)*y(k,126)
         mat(k,1853) = -rxt(k,230)*y(k,126)
         mat(k,2485) = -rxt(k,251)*y(k,126)
         mat(k,1611) = -rxt(k,282)*y(k,126)
         mat(k,1477) = -rxt(k,367)*y(k,126)
         mat(k,1445) = -rxt(k,417)*y(k,126)
         mat(k,826) = -rxt(k,475)*y(k,126)
         mat(k,780) = -rxt(k,478)*y(k,126)
         mat(k,964) = -rxt(k,481)*y(k,126)
         mat(k,541) = -rxt(k,485)*y(k,126)
         mat(k,573) = -rxt(k,488)*y(k,126)
         mat(k,1525) = -rxt(k,535)*y(k,126)
         mat(k,738) = rxt(k,419)*y(k,265)
         mat(k,407) = rxt(k,390)*y(k,127)
         mat(k,1611) = mat(k,1611) + rxt(k,281)*y(k,125)
         mat(k,2485) = mat(k,2485) + rxt(k,249)*y(k,125)
         mat(k,552) = rxt(k,208)*y(k,265)
         mat(k,638) = .700_r8*rxt(k,439)*y(k,265)
         mat(k,1304) = rxt(k,445)*y(k,125) + rxt(k,446)*y(k,127)
         mat(k,1960) = rxt(k,281)*y(k,19) + rxt(k,249)*y(k,59) + rxt(k,445)*y(k,101) &
                      + 2.000_r8*rxt(k,219)*y(k,127) + rxt(k,225)*y(k,135) &
                      + rxt(k,224)*y(k,137) + rxt(k,460)*y(k,227) + rxt(k,421) &
                      *y(k,228) + rxt(k,463)*y(k,230) + rxt(k,468)*y(k,232) &
                      + rxt(k,345)*y(k,233) + rxt(k,373)*y(k,234) + rxt(k,470) &
                      *y(k,235) + rxt(k,356)*y(k,236) + rxt(k,324)*y(k,237) &
                      + rxt(k,474)*y(k,238) + rxt(k,392)*y(k,240) + rxt(k,360) &
                      *y(k,242) + rxt(k,223)*y(k,243) + rxt(k,332)*y(k,244) &
                      + .920_r8*rxt(k,431)*y(k,245) + .920_r8*rxt(k,437)*y(k,246) &
                      + rxt(k,399)*y(k,250) + rxt(k,477)*y(k,251) + rxt(k,408) &
                      *y(k,252) + rxt(k,480)*y(k,253) + rxt(k,411)*y(k,254) &
                      + 1.600_r8*rxt(k,510)*y(k,260) + rxt(k,483)*y(k,267) &
                      + rxt(k,382)*y(k,268) + rxt(k,386)*y(k,269) + .900_r8*rxt(k,515) &
                      *y(k,270) + .800_r8*rxt(k,520)*y(k,271) + rxt(k,490)*y(k,272) &
                      + rxt(k,456)*y(k,274) + rxt(k,496)*y(k,275) + rxt(k,499) &
                      *y(k,277)
         mat(k,2432) = mat(k,2432) + rxt(k,390)*y(k,16) + rxt(k,446)*y(k,101) &
                      + 2.000_r8*rxt(k,219)*y(k,125) + rxt(k,220)*y(k,135) &
                      + rxt(k,218)*y(k,243) + rxt(k,432)*y(k,245) + rxt(k,438) &
                      *y(k,246) + rxt(k,398)*y(k,250) + rxt(k,409)*y(k,252) &
                      + 2.000_r8*rxt(k,511)*y(k,260) + rxt(k,221)*y(k,265) &
                      + rxt(k,457)*y(k,274)
         mat(k,921) = rxt(k,380)*y(k,265)
         mat(k,2049) = mat(k,2049) + rxt(k,225)*y(k,125) + rxt(k,220)*y(k,127)
         mat(k,2113) = mat(k,2113) + rxt(k,224)*y(k,125)
         mat(k,692) = rxt(k,517)*y(k,265)
         mat(k,573) = mat(k,573) + rxt(k,460)*y(k,125)
         mat(k,1092) = rxt(k,421)*y(k,125)
         mat(k,535) = rxt(k,463)*y(k,125)
         mat(k,479) = rxt(k,468)*y(k,125)
         mat(k,949) = rxt(k,345)*y(k,125)
         mat(k,984) = rxt(k,373)*y(k,125)
         mat(k,699) = rxt(k,470)*y(k,125)
         mat(k,1477) = mat(k,1477) + rxt(k,356)*y(k,125)
         mat(k,2186) = rxt(k,324)*y(k,125) + .500_r8*rxt(k,508)*y(k,260)
         mat(k,826) = mat(k,826) + rxt(k,474)*y(k,125)
         mat(k,613) = rxt(k,392)*y(k,125)
         mat(k,845) = rxt(k,360)*y(k,125)
         mat(k,2372) = mat(k,2372) + rxt(k,223)*y(k,125) + rxt(k,218)*y(k,127)
         mat(k,485) = rxt(k,332)*y(k,125)
         mat(k,1350) = .920_r8*rxt(k,431)*y(k,125) + rxt(k,432)*y(k,127)
         mat(k,1383) = .920_r8*rxt(k,437)*y(k,125) + rxt(k,438)*y(k,127)
         mat(k,1405) = rxt(k,399)*y(k,125) + rxt(k,398)*y(k,127)
         mat(k,780) = mat(k,780) + rxt(k,477)*y(k,125)
         mat(k,1445) = mat(k,1445) + rxt(k,408)*y(k,125) + rxt(k,409)*y(k,127)
         mat(k,964) = mat(k,964) + rxt(k,480)*y(k,125)
         mat(k,746) = rxt(k,411)*y(k,125)
         mat(k,1155) = 1.600_r8*rxt(k,510)*y(k,125) + 2.000_r8*rxt(k,511)*y(k,127) &
                      + .500_r8*rxt(k,508)*y(k,237)
         mat(k,1853) = mat(k,1853) + rxt(k,419)*y(k,1) + rxt(k,208)*y(k,90) &
                      + .700_r8*rxt(k,439)*y(k,99) + rxt(k,221)*y(k,127) + rxt(k,380) &
                      *y(k,128) + rxt(k,517)*y(k,214)
         mat(k,504) = rxt(k,483)*y(k,125)
         mat(k,863) = rxt(k,382)*y(k,125)
         mat(k,1265) = rxt(k,386)*y(k,125)
         mat(k,1227) = .900_r8*rxt(k,515)*y(k,125)
         mat(k,1110) = .800_r8*rxt(k,520)*y(k,125)
         mat(k,794) = rxt(k,490)*y(k,125)
         mat(k,1282) = rxt(k,456)*y(k,125) + rxt(k,457)*y(k,127)
         mat(k,819) = rxt(k,496)*y(k,125)
         mat(k,561) = rxt(k,499)*y(k,125)
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
         mat(k,2443) = -(rxt(k,218)*y(k,243) + rxt(k,219)*y(k,125) + rxt(k,220) &
                      *y(k,135) + rxt(k,221)*y(k,265) + rxt(k,229)*y(k,126) + rxt(k,318) &
                      *y(k,42) + rxt(k,350)*y(k,45) + rxt(k,369)*y(k,29) + rxt(k,376) &
                      *y(k,49) + rxt(k,390)*y(k,16) + rxt(k,398)*y(k,250) + rxt(k,409) &
                      *y(k,252) + rxt(k,432)*y(k,245) + rxt(k,438)*y(k,246) + rxt(k,441) &
                      *y(k,98) + rxt(k,446)*y(k,101) + rxt(k,457)*y(k,274) + rxt(k,502) &
                      *y(k,6) + rxt(k,505)*y(k,111) + rxt(k,511)*y(k,260) + rxt(k,522) &
                      *y(k,216) + rxt(k,525)*y(k,67))
         mat(k,2383) = -rxt(k,218)*y(k,127)
         mat(k,1971) = -rxt(k,219)*y(k,127)
         mat(k,2060) = -rxt(k,220)*y(k,127)
         mat(k,1864) = -rxt(k,221)*y(k,127)
         mat(k,1691) = -rxt(k,229)*y(k,127)
         mat(k,2469) = -rxt(k,318)*y(k,127)
         mat(k,1203) = -rxt(k,350)*y(k,127)
         mat(k,1193) = -rxt(k,369)*y(k,127)
         mat(k,1332) = -rxt(k,376)*y(k,127)
         mat(k,409) = -rxt(k,390)*y(k,127)
         mat(k,1411) = -rxt(k,398)*y(k,127)
         mat(k,1452) = -rxt(k,409)*y(k,127)
         mat(k,1357) = -rxt(k,432)*y(k,127)
         mat(k,1390) = -rxt(k,438)*y(k,127)
         mat(k,941) = -rxt(k,441)*y(k,127)
         mat(k,1310) = -rxt(k,446)*y(k,127)
         mat(k,1287) = -rxt(k,457)*y(k,127)
         mat(k,1079) = -rxt(k,502)*y(k,127)
         mat(k,1023) = -rxt(k,505)*y(k,127)
         mat(k,1161) = -rxt(k,511)*y(k,127)
         mat(k,1124) = -rxt(k,522)*y(k,127)
         mat(k,350) = -rxt(k,525)*y(k,127)
         mat(k,624) = rxt(k,283)*y(k,135)
         mat(k,2017) = rxt(k,250)*y(k,60)
         mat(k,1051) = rxt(k,250)*y(k,56) + rxt(k,252)*y(k,135) + rxt(k,253)*y(k,265)
         mat(k,977) = rxt(k,297)*y(k,89)
         mat(k,2263) = rxt(k,297)*y(k,73) + rxt(k,231)*y(k,265)
         mat(k,632) = .500_r8*rxt(k,414)*y(k,265)
         mat(k,1691) = mat(k,1691) + rxt(k,217)*y(k,135) + rxt(k,216)*y(k,137)
         mat(k,2060) = mat(k,2060) + rxt(k,283)*y(k,20) + rxt(k,252)*y(k,60) &
                      + rxt(k,217)*y(k,126)
         mat(k,2124) = rxt(k,216)*y(k,126)
         mat(k,591) = rxt(k,365)*y(k,265)
         mat(k,1864) = mat(k,1864) + rxt(k,253)*y(k,60) + rxt(k,231)*y(k,89) &
                      + .500_r8*rxt(k,414)*y(k,110) + rxt(k,365)*y(k,142)
         mat(k,918) = -(rxt(k,380)*y(k,265))
         mat(k,1807) = -rxt(k,380)*y(k,128)
         mat(k,1176) = rxt(k,369)*y(k,127)
         mat(k,594) = .500_r8*rxt(k,440)*y(k,265)
         mat(k,441) = rxt(k,447)*y(k,265)
         mat(k,464) = rxt(k,451)*y(k,265)
         mat(k,1132) = rxt(k,452)*y(k,265)
         mat(k,2394) = rxt(k,369)*y(k,29)
         mat(k,1807) = mat(k,1807) + .500_r8*rxt(k,440)*y(k,100) + rxt(k,447)*y(k,102) &
                      + rxt(k,451)*y(k,116) + rxt(k,452)*y(k,117)
         mat(k,489) = -(rxt(k,512)*y(k,265))
         mat(k,1765) = -rxt(k,512)*y(k,129)
         mat(k,2306) = rxt(k,509)*y(k,260)
         mat(k,1147) = rxt(k,509)*y(k,243)
         mat(k,2053) = -(rxt(k,187)*y(k,137) + 4._r8*rxt(k,188)*y(k,135) + rxt(k,189) &
                      *y(k,136) + rxt(k,190)*y(k,77) + rxt(k,191)*y(k,79) + rxt(k,196) &
                      *y(k,243) + rxt(k,202)*y(k,265) + (rxt(k,215) + rxt(k,217) &
                      ) * y(k,126) + rxt(k,220)*y(k,127) + rxt(k,225)*y(k,125) &
                      + rxt(k,252)*y(k,60) + rxt(k,254)*y(k,59) + rxt(k,257)*y(k,85) &
                      + rxt(k,260)*y(k,92) + rxt(k,283)*y(k,20) + rxt(k,284)*y(k,19) &
                      + rxt(k,286)*y(k,81) + rxt(k,288)*y(k,91) + rxt(k,319)*y(k,42) &
                      + rxt(k,527)*y(k,140) + (rxt(k,610) + rxt(k,611)) * y(k,257) &
                      + rxt(k,612)*y(k,259))
         mat(k,2117) = -rxt(k,187)*y(k,135)
         mat(k,1593) = -rxt(k,189)*y(k,135)
         mat(k,1511) = -rxt(k,190)*y(k,135)
         mat(k,647) = -rxt(k,191)*y(k,135)
         mat(k,2376) = -rxt(k,196)*y(k,135)
         mat(k,1857) = -rxt(k,202)*y(k,135)
         mat(k,1684) = -(rxt(k,215) + rxt(k,217)) * y(k,135)
         mat(k,2436) = -rxt(k,220)*y(k,135)
         mat(k,1964) = -rxt(k,225)*y(k,135)
         mat(k,1049) = -rxt(k,252)*y(k,135)
         mat(k,2489) = -rxt(k,254)*y(k,135)
         mat(k,1546) = -rxt(k,257)*y(k,135)
         mat(k,886) = -rxt(k,260)*y(k,135)
         mat(k,622) = -rxt(k,283)*y(k,135)
         mat(k,1615) = -rxt(k,284)*y(k,135)
         mat(k,877) = -rxt(k,286)*y(k,135)
         mat(k,838) = -rxt(k,288)*y(k,135)
         mat(k,2462) = -rxt(k,319)*y(k,135)
         mat(k,425) = -rxt(k,527)*y(k,135)
         mat(k,680) = -(rxt(k,610) + rxt(k,611)) * y(k,135)
         mat(k,549) = -rxt(k,612)*y(k,135)
         mat(k,2138) = rxt(k,194)*y(k,243)
         mat(k,1039) = rxt(k,211)*y(k,125) + rxt(k,212)*y(k,126) + rxt(k,236)*y(k,136) &
                      + rxt(k,615)*y(k,264)
         mat(k,1964) = mat(k,1964) + rxt(k,211)*y(k,113) + rxt(k,209)*y(k,256)
         mat(k,1684) = mat(k,1684) + rxt(k,212)*y(k,113)
         mat(k,1593) = mat(k,1593) + rxt(k,236)*y(k,113) + rxt(k,529)*y(k,155) &
                      + rxt(k,536)*y(k,157) + rxt(k,614)*y(k,259) + (rxt(k,176) &
                       +rxt(k,177))*y(k,261) + rxt(k,620)*y(k,266)
         mat(k,759) = rxt(k,529)*y(k,136)
         mat(k,1529) = rxt(k,536)*y(k,136)
         mat(k,916) = rxt(k,606)*y(k,258) + 1.150_r8*rxt(k,607)*y(k,264)
         mat(k,2376) = mat(k,2376) + rxt(k,194)*y(k,76)
         mat(k,708) = rxt(k,209)*y(k,125)
         mat(k,895) = rxt(k,606)*y(k,239)
         mat(k,549) = mat(k,549) + rxt(k,614)*y(k,136)
         mat(k,2233) = (rxt(k,176)+rxt(k,177))*y(k,136)
         mat(k,903) = rxt(k,615)*y(k,113) + 1.150_r8*rxt(k,607)*y(k,239)
         mat(k,1857) = mat(k,1857) + 2.000_r8*rxt(k,204)*y(k,265)
         mat(k,857) = rxt(k,620)*y(k,136)
         mat(k,1587) = -(rxt(k,176)*y(k,261) + rxt(k,181)*y(k,262) + rxt(k,189) &
                      *y(k,135) + rxt(k,195)*y(k,76) + rxt(k,234)*y(k,256) + rxt(k,236) &
                      *y(k,113) + rxt(k,362)*y(k,241) + rxt(k,529)*y(k,155) + rxt(k,536) &
                      *y(k,157) + rxt(k,609)*y(k,257) + (rxt(k,613) + rxt(k,614) &
                      ) * y(k,259) + rxt(k,620)*y(k,266))
         mat(k,2226) = -rxt(k,176)*y(k,136)
         mat(k,226) = -rxt(k,181)*y(k,136)
         mat(k,2046) = -rxt(k,189)*y(k,136)
         mat(k,2131) = -rxt(k,195)*y(k,136)
         mat(k,706) = -rxt(k,234)*y(k,136)
         mat(k,1035) = -rxt(k,236)*y(k,136)
         mat(k,513) = -rxt(k,362)*y(k,136)
         mat(k,757) = -rxt(k,529)*y(k,136)
         mat(k,1523) = -rxt(k,536)*y(k,136)
         mat(k,679) = -rxt(k,609)*y(k,136)
         mat(k,548) = -(rxt(k,613) + rxt(k,614)) * y(k,136)
         mat(k,856) = -rxt(k,620)*y(k,136)
         mat(k,1557) = rxt(k,275)*y(k,137) + rxt(k,274)*y(k,243)
         mat(k,1609) = 2.000_r8*rxt(k,276)*y(k,19) + (rxt(k,278)+rxt(k,279))*y(k,59) &
                      + rxt(k,284)*y(k,135) + rxt(k,280)*y(k,243)
         mat(k,2003) = rxt(k,243)*y(k,137) + rxt(k,241)*y(k,243)
         mat(k,2483) = (rxt(k,278)+rxt(k,279))*y(k,19) + (2.000_r8*rxt(k,245) &
                       +2.000_r8*rxt(k,246))*y(k,59) + rxt(k,254)*y(k,135) &
                      + rxt(k,248)*y(k,243) + rxt(k,256)*y(k,265)
         mat(k,2131) = mat(k,2131) + rxt(k,198)*y(k,137) + rxt(k,192)*y(k,243)
         mat(k,551) = rxt(k,208)*y(k,265)
         mat(k,1035) = mat(k,1035) + rxt(k,214)*y(k,126)
         mat(k,318) = rxt(k,227)*y(k,261)
         mat(k,1957) = rxt(k,224)*y(k,137) + rxt(k,617)*y(k,264)
         mat(k,1677) = rxt(k,214)*y(k,113) + rxt(k,215)*y(k,135) + rxt(k,216)*y(k,137)
         mat(k,2429) = rxt(k,220)*y(k,135) + rxt(k,218)*y(k,243)
         mat(k,2046) = mat(k,2046) + rxt(k,284)*y(k,19) + rxt(k,254)*y(k,59) &
                      + rxt(k,215)*y(k,126) + rxt(k,220)*y(k,127) &
                      + 2.000_r8*rxt(k,188)*y(k,135) + 2.000_r8*rxt(k,187)*y(k,137) &
                      + rxt(k,196)*y(k,243) + rxt(k,180)*y(k,262) + rxt(k,202) &
                      *y(k,265)
         mat(k,1587) = mat(k,1587) + 2.000_r8*rxt(k,181)*y(k,262)
         mat(k,2110) = rxt(k,275)*y(k,17) + rxt(k,243)*y(k,56) + rxt(k,198)*y(k,76) &
                      + rxt(k,224)*y(k,125) + rxt(k,216)*y(k,126) &
                      + 2.000_r8*rxt(k,187)*y(k,135) + rxt(k,531)*y(k,155) &
                      + rxt(k,537)*y(k,157) + 2.000_r8*rxt(k,197)*y(k,243) &
                      + 2.000_r8*rxt(k,178)*y(k,261) + rxt(k,203)*y(k,265)
         mat(k,757) = mat(k,757) + rxt(k,531)*y(k,137)
         mat(k,1523) = mat(k,1523) + rxt(k,537)*y(k,137)
         mat(k,948) = rxt(k,344)*y(k,243)
         mat(k,983) = rxt(k,372)*y(k,243)
         mat(k,2183) = rxt(k,323)*y(k,243)
         mat(k,2369) = rxt(k,274)*y(k,17) + rxt(k,280)*y(k,19) + rxt(k,241)*y(k,56) &
                      + rxt(k,248)*y(k,59) + rxt(k,192)*y(k,76) + rxt(k,218)*y(k,127) &
                      + rxt(k,196)*y(k,135) + 2.000_r8*rxt(k,197)*y(k,137) &
                      + rxt(k,344)*y(k,233) + rxt(k,372)*y(k,234) + rxt(k,323) &
                      *y(k,237) + 2.000_r8*rxt(k,206)*y(k,243) + rxt(k,201)*y(k,265) &
                      + rxt(k,381)*y(k,268)
         mat(k,2226) = mat(k,2226) + rxt(k,227)*y(k,114) + 2.000_r8*rxt(k,178) &
                      *y(k,137)
         mat(k,226) = mat(k,226) + rxt(k,180)*y(k,135) + 2.000_r8*rxt(k,181)*y(k,136)
         mat(k,901) = rxt(k,617)*y(k,125)
         mat(k,1850) = rxt(k,256)*y(k,59) + rxt(k,208)*y(k,90) + rxt(k,202)*y(k,135) &
                      + rxt(k,203)*y(k,137) + rxt(k,201)*y(k,243)
         mat(k,862) = rxt(k,381)*y(k,243)
         mat(k,2118) = -(rxt(k,178)*y(k,261) + rxt(k,187)*y(k,135) + rxt(k,197) &
                      *y(k,243) + rxt(k,198)*y(k,76) + rxt(k,203)*y(k,265) + rxt(k,216) &
                      *y(k,126) + rxt(k,224)*y(k,125) + rxt(k,243)*y(k,56) + rxt(k,275) &
                      *y(k,17) + rxt(k,341)*y(k,25) + rxt(k,370)*y(k,29) + rxt(k,401) &
                      *y(k,106) + rxt(k,415)*y(k,112) + rxt(k,448)*y(k,98) + rxt(k,486) &
                      *y(k,144) + rxt(k,503)*y(k,6) + rxt(k,506)*y(k,111) + rxt(k,531) &
                      *y(k,155) + rxt(k,537)*y(k,157))
         mat(k,2234) = -rxt(k,178)*y(k,137)
         mat(k,2054) = -rxt(k,187)*y(k,137)
         mat(k,2377) = -rxt(k,197)*y(k,137)
         mat(k,2139) = -rxt(k,198)*y(k,137)
         mat(k,1858) = -rxt(k,203)*y(k,137)
         mat(k,1685) = -rxt(k,216)*y(k,137)
         mat(k,1965) = -rxt(k,224)*y(k,137)
         mat(k,2011) = -rxt(k,243)*y(k,137)
         mat(k,1563) = -rxt(k,275)*y(k,137)
         mat(k,582) = -rxt(k,341)*y(k,137)
         mat(k,1190) = -rxt(k,370)*y(k,137)
         mat(k,1321) = -rxt(k,401)*y(k,137)
         mat(k,1430) = -rxt(k,415)*y(k,137)
         mat(k,938) = -rxt(k,448)*y(k,137)
         mat(k,542) = -rxt(k,486)*y(k,137)
         mat(k,1077) = -rxt(k,503)*y(k,137)
         mat(k,1021) = -rxt(k,506)*y(k,137)
         mat(k,760) = -rxt(k,531)*y(k,137)
         mat(k,1530) = -rxt(k,537)*y(k,137)
         mat(k,2054) = mat(k,2054) + rxt(k,189)*y(k,136)
         mat(k,1594) = rxt(k,189)*y(k,135)
         mat(k,1480) = .150_r8*rxt(k,355)*y(k,243)
         mat(k,2377) = mat(k,2377) + .150_r8*rxt(k,355)*y(k,236) + .150_r8*rxt(k,406) &
                      *y(k,252)
         mat(k,1448) = .150_r8*rxt(k,406)*y(k,243)
         mat(k,393) = -(rxt(k,538)*y(k,157))
         mat(k,1518) = -rxt(k,538)*y(k,139)
         mat(k,1602) = rxt(k,277)*y(k,59)
         mat(k,2476) = rxt(k,277)*y(k,19) + 2.000_r8*rxt(k,247)*y(k,59)
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
         mat(k,419) = -(rxt(k,527)*y(k,135) + rxt(k,528)*y(k,265))
         mat(k,2023) = -rxt(k,527)*y(k,140)
         mat(k,1755) = -rxt(k,528)*y(k,140)
         mat(k,1235) = rxt(k,394)*y(k,265)
         mat(k,1891) = .100_r8*rxt(k,515)*y(k,270)
         mat(k,1732) = rxt(k,394)*y(k,93)
         mat(k,1216) = .100_r8*rxt(k,515)*y(k,125)
         mat(k,585) = -(rxt(k,365)*y(k,265))
         mat(k,1777) = -rxt(k,365)*y(k,142)
         mat(k,1656) = rxt(k,367)*y(k,236)
         mat(k,1456) = rxt(k,367)*y(k,126)
         mat(k,1650) = rxt(k,488)*y(k,227)
         mat(k,570) = rxt(k,488)*y(k,126)
         mat(k,539) = -(rxt(k,485)*y(k,126) + rxt(k,486)*y(k,137))
         mat(k,1653) = -rxt(k,485)*y(k,144)
         mat(k,2073) = -rxt(k,486)*y(k,144)
         mat(k,240) = .070_r8*rxt(k,472)*y(k,265)
         mat(k,1902) = rxt(k,470)*y(k,235)
         mat(k,222) = .060_r8*rxt(k,484)*y(k,265)
         mat(k,265) = .070_r8*rxt(k,500)*y(k,265)
         mat(k,697) = rxt(k,470)*y(k,125)
         mat(k,1771) = .070_r8*rxt(k,472)*y(k,66) + .060_r8*rxt(k,484)*y(k,145) &
                      + .070_r8*rxt(k,500)*y(k,223)
         mat(k,220) = -(rxt(k,484)*y(k,265))
         mat(k,1723) = -rxt(k,484)*y(k,145)
         mat(k,212) = .530_r8*rxt(k,461)*y(k,265)
         mat(k,1723) = mat(k,1723) + .530_r8*rxt(k,461)*y(k,7)
         mat(k,364) = -(rxt(k,487)*y(k,265))
         mat(k,1747) = -rxt(k,487)*y(k,146)
         mat(k,2293) = rxt(k,482)*y(k,267)
         mat(k,501) = rxt(k,482)*y(k,243)
         mat(k,601) = -(rxt(k,383)*y(k,265))
         mat(k,1779) = -rxt(k,383)*y(k,153)
         mat(k,2316) = rxt(k,381)*y(k,268)
         mat(k,858) = rxt(k,381)*y(k,243)
         mat(k,427) = -(rxt(k,387)*y(k,265))
         mat(k,1756) = -rxt(k,387)*y(k,154)
         mat(k,2298) = .850_r8*rxt(k,385)*y(k,269)
         mat(k,1259) = .850_r8*rxt(k,385)*y(k,243)
         mat(k,755) = -(rxt(k,529)*y(k,136) + rxt(k,531)*y(k,137) + rxt(k,534) &
                      *y(k,265))
         mat(k,1575) = -rxt(k,529)*y(k,155)
         mat(k,2077) = -rxt(k,531)*y(k,155)
         mat(k,1795) = -rxt(k,534)*y(k,155)
         mat(k,1521) = -(rxt(k,532)*y(k,19) + rxt(k,533)*y(k,59) + rxt(k,535)*y(k,126) &
                      + rxt(k,536)*y(k,136) + rxt(k,537)*y(k,137) + rxt(k,538) &
                      *y(k,139) + rxt(k,539)*y(k,265))
         mat(k,1606) = -rxt(k,532)*y(k,157)
         mat(k,2480) = -rxt(k,533)*y(k,157)
         mat(k,1674) = -rxt(k,535)*y(k,157)
         mat(k,1585) = -rxt(k,536)*y(k,157)
         mat(k,2108) = -rxt(k,537)*y(k,157)
         mat(k,395) = -rxt(k,538)*y(k,157)
         mat(k,1847) = -rxt(k,539)*y(k,157)
         mat(k,2043) = rxt(k,527)*y(k,140)
         mat(k,1585) = mat(k,1585) + rxt(k,529)*y(k,155)
         mat(k,2108) = mat(k,2108) + rxt(k,531)*y(k,155)
         mat(k,423) = rxt(k,527)*y(k,135)
         mat(k,756) = rxt(k,529)*y(k,136) + rxt(k,531)*y(k,137) + rxt(k,534)*y(k,265)
         mat(k,1847) = mat(k,1847) + rxt(k,534)*y(k,155)
         mat(k,992) = -(rxt(k,530)*y(k,265))
         mat(k,1814) = -rxt(k,530)*y(k,158)
         mat(k,1605) = rxt(k,532)*y(k,157)
         mat(k,2478) = rxt(k,533)*y(k,157)
         mat(k,346) = rxt(k,525)*y(k,127) + (rxt(k,526)+.500_r8*rxt(k,540))*y(k,265)
         mat(k,1664) = rxt(k,535)*y(k,157)
         mat(k,2397) = rxt(k,525)*y(k,67)
         mat(k,1580) = rxt(k,536)*y(k,157)
         mat(k,2081) = rxt(k,537)*y(k,157)
         mat(k,394) = rxt(k,538)*y(k,157)
         mat(k,421) = rxt(k,528)*y(k,265)
         mat(k,1520) = rxt(k,532)*y(k,19) + rxt(k,533)*y(k,59) + rxt(k,535)*y(k,126) &
                      + rxt(k,536)*y(k,136) + rxt(k,537)*y(k,137) + rxt(k,538) &
                      *y(k,139) + rxt(k,539)*y(k,265)
         mat(k,1814) = mat(k,1814) + (rxt(k,526)+.500_r8*rxt(k,540))*y(k,67) &
                      + rxt(k,528)*y(k,140) + rxt(k,539)*y(k,157)
         mat(k,309) = -(rxt(k,541)*y(k,278))
         mat(k,2501) = -rxt(k,541)*y(k,159)
         mat(k,991) = rxt(k,530)*y(k,265)
         mat(k,1740) = rxt(k,530)*y(k,158)
         mat(k,1868) = .1056005_r8*rxt(k,570)*y(k,248)
         mat(k,79) = .5931005_r8*rxt(k,580)*y(k,265)
         mat(k,2267) = .2381005_r8*rxt(k,569)*y(k,248)
         mat(k,109) = .1056005_r8*rxt(k,570)*y(k,125) + .2381005_r8*rxt(k,569) &
                      *y(k,243)
         mat(k,1695) = .5931005_r8*rxt(k,580)*y(k,210)
         mat(k,1869) = .1026005_r8*rxt(k,570)*y(k,248)
         mat(k,80) = .1534005_r8*rxt(k,580)*y(k,265)
         mat(k,2268) = .1308005_r8*rxt(k,569)*y(k,248)
         mat(k,110) = .1026005_r8*rxt(k,570)*y(k,125) + .1308005_r8*rxt(k,569) &
                      *y(k,243)
         mat(k,1696) = .1534005_r8*rxt(k,580)*y(k,210)
         mat(k,1870) = .0521005_r8*rxt(k,570)*y(k,248)
         mat(k,81) = .0459005_r8*rxt(k,580)*y(k,265)
         mat(k,2269) = .0348005_r8*rxt(k,569)*y(k,248)
         mat(k,111) = .0521005_r8*rxt(k,570)*y(k,125) + .0348005_r8*rxt(k,569) &
                      *y(k,243)
         mat(k,1697) = .0459005_r8*rxt(k,580)*y(k,210)
         mat(k,1871) = .0143005_r8*rxt(k,570)*y(k,248)
         mat(k,82) = .0085005_r8*rxt(k,580)*y(k,265)
         mat(k,2270) = .0076005_r8*rxt(k,569)*y(k,248)
         mat(k,112) = .0143005_r8*rxt(k,570)*y(k,125) + .0076005_r8*rxt(k,569) &
                      *y(k,243)
         mat(k,1698) = .0085005_r8*rxt(k,580)*y(k,210)
         mat(k,1872) = .0166005_r8*rxt(k,570)*y(k,248)
         mat(k,83) = .0128005_r8*rxt(k,580)*y(k,265)
         mat(k,2271) = .0113005_r8*rxt(k,569)*y(k,248)
         mat(k,113) = .0166005_r8*rxt(k,570)*y(k,125) + .0113005_r8*rxt(k,569) &
                      *y(k,243)
         mat(k,1699) = .0128005_r8*rxt(k,580)*y(k,210)
         mat(k,1054) = .2202005_r8*rxt(k,559)*y(k,137)
         mat(k,998) = .0508005_r8*rxt(k,578)*y(k,137)
         mat(k,1873) = .1279005_r8*rxt(k,558)*y(k,229) + .0003005_r8*rxt(k,566) &
                      *y(k,247) + .0245005_r8*rxt(k,577)*y(k,255)
         mat(k,2064) = .2202005_r8*rxt(k,559)*y(k,6) + .0508005_r8*rxt(k,578)*y(k,111)
         mat(k,91) = .1279005_r8*rxt(k,558)*y(k,125) + .2202005_r8*rxt(k,557)*y(k,243)
         mat(k,2272) = .2202005_r8*rxt(k,557)*y(k,229) + .0031005_r8*rxt(k,565) &
                      *y(k,247) + .0508005_r8*rxt(k,576)*y(k,255)
         mat(k,103) = .0003005_r8*rxt(k,566)*y(k,125) + .0031005_r8*rxt(k,565) &
                      *y(k,243)
         mat(k,125) = .0245005_r8*rxt(k,577)*y(k,125) + .0508005_r8*rxt(k,576) &
                      *y(k,243)
         mat(k,1055) = .2067005_r8*rxt(k,559)*y(k,137)
         mat(k,999) = .1149005_r8*rxt(k,578)*y(k,137)
         mat(k,1874) = .1792005_r8*rxt(k,558)*y(k,229) + .0003005_r8*rxt(k,566) &
                      *y(k,247) + .0082005_r8*rxt(k,577)*y(k,255)
         mat(k,2065) = .2067005_r8*rxt(k,559)*y(k,6) + .1149005_r8*rxt(k,578)*y(k,111)
         mat(k,92) = .1792005_r8*rxt(k,558)*y(k,125) + .2067005_r8*rxt(k,557)*y(k,243)
         mat(k,2273) = .2067005_r8*rxt(k,557)*y(k,229) + .0035005_r8*rxt(k,565) &
                      *y(k,247) + .1149005_r8*rxt(k,576)*y(k,255)
         mat(k,104) = .0003005_r8*rxt(k,566)*y(k,125) + .0035005_r8*rxt(k,565) &
                      *y(k,243)
         mat(k,126) = .0082005_r8*rxt(k,577)*y(k,125) + .1149005_r8*rxt(k,576) &
                      *y(k,243)
         mat(k,1056) = .0653005_r8*rxt(k,559)*y(k,137)
         mat(k,1000) = .0348005_r8*rxt(k,578)*y(k,137)
         mat(k,1875) = .0676005_r8*rxt(k,558)*y(k,229) + .0073005_r8*rxt(k,566) &
                      *y(k,247) + .0772005_r8*rxt(k,577)*y(k,255)
         mat(k,2066) = .0653005_r8*rxt(k,559)*y(k,6) + .0348005_r8*rxt(k,578)*y(k,111)
         mat(k,93) = .0676005_r8*rxt(k,558)*y(k,125) + .0653005_r8*rxt(k,557)*y(k,243)
         mat(k,2274) = .0653005_r8*rxt(k,557)*y(k,229) + .0003005_r8*rxt(k,565) &
                      *y(k,247) + .0348005_r8*rxt(k,576)*y(k,255)
         mat(k,105) = .0073005_r8*rxt(k,566)*y(k,125) + .0003005_r8*rxt(k,565) &
                      *y(k,243)
         mat(k,127) = .0772005_r8*rxt(k,577)*y(k,125) + .0348005_r8*rxt(k,576) &
                      *y(k,243)
         mat(k,1057) = .1749305_r8*rxt(k,556)*y(k,127) + .1284005_r8*rxt(k,559) &
                      *y(k,137)
         mat(k,924) = .0590245_r8*rxt(k,564)*y(k,127) + .0033005_r8*rxt(k,567) &
                      *y(k,137)
         mat(k,1001) = .1749305_r8*rxt(k,575)*y(k,127) + .0554005_r8*rxt(k,578) &
                      *y(k,137)
         mat(k,1876) = .079_r8*rxt(k,558)*y(k,229) + .0057005_r8*rxt(k,566)*y(k,247) &
                      + .0332005_r8*rxt(k,577)*y(k,255)
         mat(k,2387) = .1749305_r8*rxt(k,556)*y(k,6) + .0590245_r8*rxt(k,564)*y(k,98) &
                      + .1749305_r8*rxt(k,575)*y(k,111)
         mat(k,2067) = .1284005_r8*rxt(k,559)*y(k,6) + .0033005_r8*rxt(k,567)*y(k,98) &
                      + .0554005_r8*rxt(k,578)*y(k,111)
         mat(k,94) = .079_r8*rxt(k,558)*y(k,125) + .1284005_r8*rxt(k,557)*y(k,243)
         mat(k,2275) = .1284005_r8*rxt(k,557)*y(k,229) + .0271005_r8*rxt(k,565) &
                      *y(k,247) + .0554005_r8*rxt(k,576)*y(k,255)
         mat(k,106) = .0057005_r8*rxt(k,566)*y(k,125) + .0271005_r8*rxt(k,565) &
                      *y(k,243)
         mat(k,128) = .0332005_r8*rxt(k,577)*y(k,125) + .0554005_r8*rxt(k,576) &
                      *y(k,243)
         mat(k,1058) = .5901905_r8*rxt(k,556)*y(k,127) + .114_r8*rxt(k,559)*y(k,137)
         mat(k,925) = .0250245_r8*rxt(k,564)*y(k,127)
         mat(k,1002) = .5901905_r8*rxt(k,575)*y(k,127) + .1278005_r8*rxt(k,578) &
                      *y(k,137)
         mat(k,1877) = .1254005_r8*rxt(k,558)*y(k,229) + .0623005_r8*rxt(k,566) &
                      *y(k,247) + .130_r8*rxt(k,577)*y(k,255)
         mat(k,2388) = .5901905_r8*rxt(k,556)*y(k,6) + .0250245_r8*rxt(k,564)*y(k,98) &
                      + .5901905_r8*rxt(k,575)*y(k,111)
         mat(k,2068) = .114_r8*rxt(k,559)*y(k,6) + .1278005_r8*rxt(k,578)*y(k,111)
         mat(k,95) = .1254005_r8*rxt(k,558)*y(k,125) + .114_r8*rxt(k,557)*y(k,243)
         mat(k,2276) = .114_r8*rxt(k,557)*y(k,229) + .0474005_r8*rxt(k,565)*y(k,247) &
                      + .1278005_r8*rxt(k,576)*y(k,255)
         mat(k,107) = .0623005_r8*rxt(k,566)*y(k,125) + .0474005_r8*rxt(k,565) &
                      *y(k,243)
         mat(k,129) = .130_r8*rxt(k,577)*y(k,125) + .1278005_r8*rxt(k,576)*y(k,243)
         mat(k,1878) = .0097005_r8*rxt(k,563)*y(k,231) + .1056005_r8*rxt(k,573) &
                      *y(k,249) + .0154005_r8*rxt(k,584)*y(k,273) &
                      + .0063005_r8*rxt(k,588)*y(k,276)
         mat(k,85) = .5931005_r8*rxt(k,581)*y(k,265)
         mat(k,97) = .0097005_r8*rxt(k,563)*y(k,125) + .0023005_r8*rxt(k,562)*y(k,243)
         mat(k,2277) = .0023005_r8*rxt(k,562)*y(k,231) + .2381005_r8*rxt(k,572) &
                      *y(k,249) + .1364005_r8*rxt(k,583)*y(k,273) &
                      + .1677005_r8*rxt(k,587)*y(k,276)
         mat(k,117) = .1056005_r8*rxt(k,573)*y(k,125) + .2381005_r8*rxt(k,572) &
                      *y(k,243)
         mat(k,1700) = .5931005_r8*rxt(k,581)*y(k,211)
         mat(k,131) = .0154005_r8*rxt(k,584)*y(k,125) + .1364005_r8*rxt(k,583) &
                      *y(k,243)
         mat(k,137) = .0063005_r8*rxt(k,588)*y(k,125) + .1677005_r8*rxt(k,587) &
                      *y(k,243)
         mat(k,1879) = .0034005_r8*rxt(k,563)*y(k,231) + .1026005_r8*rxt(k,573) &
                      *y(k,249) + .0452005_r8*rxt(k,584)*y(k,273) &
                      + .0237005_r8*rxt(k,588)*y(k,276)
         mat(k,86) = .1534005_r8*rxt(k,581)*y(k,265)
         mat(k,98) = .0034005_r8*rxt(k,563)*y(k,125) + .0008005_r8*rxt(k,562)*y(k,243)
         mat(k,2278) = .0008005_r8*rxt(k,562)*y(k,231) + .1308005_r8*rxt(k,572) &
                      *y(k,249) + .0101005_r8*rxt(k,583)*y(k,273) &
                      + .0174005_r8*rxt(k,587)*y(k,276)
         mat(k,118) = .1026005_r8*rxt(k,573)*y(k,125) + .1308005_r8*rxt(k,572) &
                      *y(k,243)
         mat(k,1701) = .1534005_r8*rxt(k,581)*y(k,211)
         mat(k,132) = .0452005_r8*rxt(k,584)*y(k,125) + .0101005_r8*rxt(k,583) &
                      *y(k,243)
         mat(k,138) = .0237005_r8*rxt(k,588)*y(k,125) + .0174005_r8*rxt(k,587) &
                      *y(k,243)
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
         mat(k,1880) = .1579005_r8*rxt(k,563)*y(k,231) + .0521005_r8*rxt(k,573) &
                      *y(k,249) + .0966005_r8*rxt(k,584)*y(k,273) &
                      + .0025005_r8*rxt(k,588)*y(k,276)
         mat(k,87) = .0459005_r8*rxt(k,581)*y(k,265)
         mat(k,99) = .1579005_r8*rxt(k,563)*y(k,125) + .0843005_r8*rxt(k,562)*y(k,243)
         mat(k,2279) = .0843005_r8*rxt(k,562)*y(k,231) + .0348005_r8*rxt(k,572) &
                      *y(k,249) + .0763005_r8*rxt(k,583)*y(k,273) + .086_r8*rxt(k,587) &
                      *y(k,276)
         mat(k,119) = .0521005_r8*rxt(k,573)*y(k,125) + .0348005_r8*rxt(k,572) &
                      *y(k,243)
         mat(k,1702) = .0459005_r8*rxt(k,581)*y(k,211)
         mat(k,133) = .0966005_r8*rxt(k,584)*y(k,125) + .0763005_r8*rxt(k,583) &
                      *y(k,243)
         mat(k,139) = .0025005_r8*rxt(k,588)*y(k,125) + .086_r8*rxt(k,587)*y(k,243)
         mat(k,1881) = .0059005_r8*rxt(k,563)*y(k,231) + .0143005_r8*rxt(k,573) &
                      *y(k,249) + .0073005_r8*rxt(k,584)*y(k,273) + .011_r8*rxt(k,588) &
                      *y(k,276)
         mat(k,88) = .0085005_r8*rxt(k,581)*y(k,265)
         mat(k,100) = .0059005_r8*rxt(k,563)*y(k,125) + .0443005_r8*rxt(k,562) &
                      *y(k,243)
         mat(k,2280) = .0443005_r8*rxt(k,562)*y(k,231) + .0076005_r8*rxt(k,572) &
                      *y(k,249) + .2157005_r8*rxt(k,583)*y(k,273) &
                      + .0512005_r8*rxt(k,587)*y(k,276)
         mat(k,120) = .0143005_r8*rxt(k,573)*y(k,125) + .0076005_r8*rxt(k,572) &
                      *y(k,243)
         mat(k,1703) = .0085005_r8*rxt(k,581)*y(k,211)
         mat(k,134) = .0073005_r8*rxt(k,584)*y(k,125) + .2157005_r8*rxt(k,583) &
                      *y(k,243)
         mat(k,140) = .011_r8*rxt(k,588)*y(k,125) + .0512005_r8*rxt(k,587)*y(k,243)
         mat(k,1882) = .0536005_r8*rxt(k,563)*y(k,231) + .0166005_r8*rxt(k,573) &
                      *y(k,249) + .238_r8*rxt(k,584)*y(k,273) + .1185005_r8*rxt(k,588) &
                      *y(k,276)
         mat(k,89) = .0128005_r8*rxt(k,581)*y(k,265)
         mat(k,101) = .0536005_r8*rxt(k,563)*y(k,125) + .1621005_r8*rxt(k,562) &
                      *y(k,243)
         mat(k,2281) = .1621005_r8*rxt(k,562)*y(k,231) + .0113005_r8*rxt(k,572) &
                      *y(k,249) + .0738005_r8*rxt(k,583)*y(k,273) &
                      + .1598005_r8*rxt(k,587)*y(k,276)
         mat(k,121) = .0166005_r8*rxt(k,573)*y(k,125) + .0113005_r8*rxt(k,572) &
                      *y(k,243)
         mat(k,1704) = .0128005_r8*rxt(k,581)*y(k,211)
         mat(k,135) = .238_r8*rxt(k,584)*y(k,125) + .0738005_r8*rxt(k,583)*y(k,243)
         mat(k,141) = .1185005_r8*rxt(k,588)*y(k,125) + .1598005_r8*rxt(k,587) &
                      *y(k,243)
         mat(k,84) = -(rxt(k,580)*y(k,265))
         mat(k,1705) = -rxt(k,580)*y(k,210)
         mat(k,90) = -(rxt(k,581)*y(k,265))
         mat(k,1706) = -rxt(k,581)*y(k,211)
         mat(k,233) = .100_r8*rxt(k,492)*y(k,265)
         mat(k,255) = .230_r8*rxt(k,494)*y(k,265)
         mat(k,1726) = .100_r8*rxt(k,492)*y(k,219) + .230_r8*rxt(k,494)*y(k,221)
         mat(k,721) = -(rxt(k,516)*y(k,265))
         mat(k,1791) = -rxt(k,516)*y(k,213)
         mat(k,2323) = rxt(k,514)*y(k,270)
         mat(k,1217) = rxt(k,514)*y(k,243)
         mat(k,690) = -(rxt(k,517)*y(k,265))
         mat(k,1788) = -rxt(k,517)*y(k,214)
         mat(k,1910) = .200_r8*rxt(k,510)*y(k,260) + .200_r8*rxt(k,520)*y(k,271)
         mat(k,2153) = .500_r8*rxt(k,508)*y(k,260)
         mat(k,1148) = .200_r8*rxt(k,510)*y(k,125) + .500_r8*rxt(k,508)*y(k,237)
         mat(k,1105) = .200_r8*rxt(k,520)*y(k,125)
         mat(k,523) = -(rxt(k,521)*y(k,265))
         mat(k,1769) = -rxt(k,521)*y(k,215)
         mat(k,2310) = rxt(k,519)*y(k,271)
         mat(k,1104) = rxt(k,519)*y(k,243)
         mat(k,1117) = -(rxt(k,522)*y(k,127) + rxt(k,523)*y(k,265))
         mat(k,2403) = -rxt(k,522)*y(k,216)
         mat(k,1822) = -rxt(k,523)*y(k,216)
         mat(k,1067) = .330_r8*rxt(k,503)*y(k,137)
         mat(k,1011) = .330_r8*rxt(k,506)*y(k,137)
         mat(k,1933) = .800_r8*rxt(k,510)*y(k,260) + .800_r8*rxt(k,520)*y(k,271)
         mat(k,2403) = mat(k,2403) + rxt(k,511)*y(k,260)
         mat(k,2087) = .330_r8*rxt(k,503)*y(k,6) + .330_r8*rxt(k,506)*y(k,111)
         mat(k,691) = rxt(k,517)*y(k,265)
         mat(k,2161) = .500_r8*rxt(k,508)*y(k,260) + rxt(k,518)*y(k,271)
         mat(k,1150) = .800_r8*rxt(k,510)*y(k,125) + rxt(k,511)*y(k,127) &
                      + .500_r8*rxt(k,508)*y(k,237)
         mat(k,1822) = mat(k,1822) + rxt(k,517)*y(k,214)
         mat(k,1108) = .800_r8*rxt(k,520)*y(k,125) + rxt(k,518)*y(k,237)
         mat(k,1164) = -(rxt(k,524)*y(k,265))
         mat(k,1826) = -rxt(k,524)*y(k,217)
         mat(k,1070) = .300_r8*rxt(k,503)*y(k,137)
         mat(k,1014) = .300_r8*rxt(k,506)*y(k,137)
         mat(k,1937) = .900_r8*rxt(k,515)*y(k,270)
         mat(k,2090) = .300_r8*rxt(k,503)*y(k,6) + .300_r8*rxt(k,506)*y(k,111)
         mat(k,2165) = rxt(k,513)*y(k,270)
         mat(k,1220) = .900_r8*rxt(k,515)*y(k,125) + rxt(k,513)*y(k,237)
         mat(k,652) = -(rxt(k,491)*y(k,265))
         mat(k,1784) = -rxt(k,491)*y(k,218)
         mat(k,2318) = rxt(k,489)*y(k,272)
         mat(k,785) = rxt(k,489)*y(k,243)
         mat(k,231) = -((rxt(k,492) + rxt(k,582)) * y(k,265))
         mat(k,1724) = -(rxt(k,492) + rxt(k,582)) * y(k,219)
         mat(k,247) = -(rxt(k,458)*y(k,265))
         mat(k,1727) = -rxt(k,458)*y(k,220)
         mat(k,2291) = rxt(k,455)*y(k,274)
         mat(k,1272) = rxt(k,455)*y(k,243)
         mat(k,256) = -(rxt(k,494)*y(k,265))
         mat(k,1729) = -rxt(k,494)*y(k,221)
         mat(k,766) = -(rxt(k,497)*y(k,265))
         mat(k,1796) = -rxt(k,497)*y(k,222)
         mat(k,2326) = rxt(k,495)*y(k,275)
         mat(k,810) = rxt(k,495)*y(k,243)
         mat(k,264) = -(rxt(k,500)*y(k,265))
         mat(k,1730) = -rxt(k,500)*y(k,223)
         mat(k,257) = .150_r8*rxt(k,494)*y(k,265)
         mat(k,1730) = mat(k,1730) + .150_r8*rxt(k,494)*y(k,221)
         mat(k,469) = -(rxt(k,501)*y(k,265))
         mat(k,1763) = -rxt(k,501)*y(k,224)
         mat(k,2303) = rxt(k,498)*y(k,277)
         mat(k,557) = rxt(k,498)*y(k,243)
         mat(k,571) = -(rxt(k,459)*y(k,243) + rxt(k,460)*y(k,125) + rxt(k,488) &
                      *y(k,126))
         mat(k,2315) = -rxt(k,459)*y(k,227)
         mat(k,1905) = -rxt(k,460)*y(k,227)
         mat(k,1655) = -rxt(k,488)*y(k,227)
         mat(k,279) = rxt(k,465)*y(k,265)
         mat(k,1775) = rxt(k,465)*y(k,22)
         mat(k,1087) = -(rxt(k,420)*y(k,243) + (rxt(k,421) + rxt(k,422)) * y(k,125))
         mat(k,2342) = -rxt(k,420)*y(k,228)
         mat(k,1930) = -(rxt(k,421) + rxt(k,422)) * y(k,228)
         mat(k,714) = rxt(k,423)*y(k,265)
         mat(k,270) = rxt(k,424)*y(k,265)
         mat(k,1819) = rxt(k,423)*y(k,2) + rxt(k,424)*y(k,15)
         mat(k,96) = -(rxt(k,557)*y(k,243) + rxt(k,558)*y(k,125))
         mat(k,2282) = -rxt(k,557)*y(k,229)
         mat(k,1883) = -rxt(k,558)*y(k,229)
         mat(k,1059) = rxt(k,560)*y(k,265)
         mat(k,1707) = rxt(k,560)*y(k,6)
         mat(k,532) = -(rxt(k,462)*y(k,243) + rxt(k,463)*y(k,125))
         mat(k,2311) = -rxt(k,462)*y(k,230)
         mat(k,1901) = -rxt(k,463)*y(k,230)
         mat(k,213) = .350_r8*rxt(k,461)*y(k,265)
         mat(k,459) = rxt(k,464)*y(k,265)
         mat(k,1770) = .350_r8*rxt(k,461)*y(k,7) + rxt(k,464)*y(k,8)
         mat(k,102) = -(rxt(k,562)*y(k,243) + rxt(k,563)*y(k,125))
         mat(k,2283) = -rxt(k,562)*y(k,231)
         mat(k,1884) = -rxt(k,563)*y(k,231)
         mat(k,209) = rxt(k,561)*y(k,265)
         mat(k,1708) = rxt(k,561)*y(k,7)
         mat(k,477) = -(rxt(k,466)*y(k,243) + rxt(k,468)*y(k,125))
         mat(k,2304) = -rxt(k,466)*y(k,232)
         mat(k,1896) = -rxt(k,468)*y(k,232)
         mat(k,371) = rxt(k,467)*y(k,265)
         mat(k,234) = .070_r8*rxt(k,492)*y(k,265)
         mat(k,258) = .060_r8*rxt(k,494)*y(k,265)
         mat(k,1764) = rxt(k,467)*y(k,23) + .070_r8*rxt(k,492)*y(k,219) &
                      + .060_r8*rxt(k,494)*y(k,221)
         mat(k,946) = -(4._r8*rxt(k,342)*y(k,233) + rxt(k,343)*y(k,237) + rxt(k,344) &
                      *y(k,243) + rxt(k,345)*y(k,125))
         mat(k,2157) = -rxt(k,343)*y(k,233)
         mat(k,2338) = -rxt(k,344)*y(k,233)
         mat(k,1925) = -rxt(k,345)*y(k,233)
         mat(k,376) = .500_r8*rxt(k,347)*y(k,265)
         mat(k,328) = rxt(k,348)*y(k,56) + rxt(k,349)*y(k,265)
         mat(k,1987) = rxt(k,348)*y(k,28)
         mat(k,1809) = .500_r8*rxt(k,347)*y(k,27) + rxt(k,349)*y(k,28)
         mat(k,980) = -(rxt(k,371)*y(k,237) + rxt(k,372)*y(k,243) + rxt(k,373) &
                      *y(k,125))
         mat(k,2158) = -rxt(k,371)*y(k,234)
         mat(k,2341) = -rxt(k,372)*y(k,234)
         mat(k,1928) = -rxt(k,373)*y(k,234)
         mat(k,434) = rxt(k,374)*y(k,265)
         mat(k,334) = rxt(k,378)*y(k,56) + rxt(k,375)*y(k,265)
         mat(k,1989) = rxt(k,378)*y(k,31)
         mat(k,1813) = rxt(k,374)*y(k,30) + rxt(k,375)*y(k,31)
         mat(k,698) = -(rxt(k,469)*y(k,243) + rxt(k,470)*y(k,125))
         mat(k,2321) = -rxt(k,469)*y(k,235)
         mat(k,1911) = -rxt(k,470)*y(k,235)
         mat(k,315) = rxt(k,471)*y(k,265)
         mat(k,1911) = mat(k,1911) + rxt(k,460)*y(k,227)
         mat(k,2075) = rxt(k,486)*y(k,144)
         mat(k,540) = rxt(k,486)*y(k,137)
         mat(k,572) = rxt(k,460)*y(k,125) + .400_r8*rxt(k,459)*y(k,243)
         mat(k,2321) = mat(k,2321) + .400_r8*rxt(k,459)*y(k,227)
         mat(k,1789) = rxt(k,471)*y(k,32)
         mat(k,1473) = -(4._r8*rxt(k,353)*y(k,236) + rxt(k,354)*y(k,237) + rxt(k,355) &
                      *y(k,243) + rxt(k,356)*y(k,125) + rxt(k,367)*y(k,126) + rxt(k,395) &
                      *y(k,250) + rxt(k,428)*y(k,245) + rxt(k,433)*y(k,246) + rxt(k,442) &
                      *y(k,101) + rxt(k,453)*y(k,274))
         mat(k,2181) = -rxt(k,354)*y(k,236)
         mat(k,2364) = -rxt(k,355)*y(k,236)
         mat(k,1954) = -rxt(k,356)*y(k,236)
         mat(k,1672) = -rxt(k,367)*y(k,236)
         mat(k,1403) = -rxt(k,395)*y(k,236)
         mat(k,1348) = -rxt(k,428)*y(k,236)
         mat(k,1381) = -rxt(k,433)*y(k,236)
         mat(k,1302) = -rxt(k,442)*y(k,236)
         mat(k,1280) = -rxt(k,453)*y(k,236)
         mat(k,1074) = .060_r8*rxt(k,503)*y(k,137)
         mat(k,1198) = rxt(k,350)*y(k,127) + rxt(k,351)*y(k,265)
         mat(k,1327) = rxt(k,376)*y(k,127) + rxt(k,377)*y(k,265)
         mat(k,666) = .500_r8*rxt(k,358)*y(k,265)
         mat(k,936) = .080_r8*rxt(k,448)*y(k,137)
         mat(k,1318) = .100_r8*rxt(k,401)*y(k,137)
         mat(k,1018) = .060_r8*rxt(k,506)*y(k,137)
         mat(k,1423) = .280_r8*rxt(k,415)*y(k,137)
         mat(k,1954) = mat(k,1954) + .530_r8*rxt(k,399)*y(k,250) + rxt(k,408)*y(k,252) &
                      + rxt(k,411)*y(k,254) + rxt(k,386)*y(k,269)
         mat(k,2425) = rxt(k,350)*y(k,45) + rxt(k,376)*y(k,49) + .530_r8*rxt(k,398) &
                      *y(k,250) + rxt(k,409)*y(k,252)
         mat(k,2106) = .060_r8*rxt(k,503)*y(k,6) + .080_r8*rxt(k,448)*y(k,98) &
                      + .100_r8*rxt(k,401)*y(k,106) + .060_r8*rxt(k,506)*y(k,111) &
                      + .280_r8*rxt(k,415)*y(k,112)
         mat(k,1167) = .650_r8*rxt(k,524)*y(k,265)
         mat(k,1473) = mat(k,1473) + .530_r8*rxt(k,395)*y(k,250)
         mat(k,2181) = mat(k,2181) + .260_r8*rxt(k,396)*y(k,250) + rxt(k,405)*y(k,252) &
                      + .300_r8*rxt(k,384)*y(k,269)
         mat(k,2364) = mat(k,2364) + .450_r8*rxt(k,406)*y(k,252) + .200_r8*rxt(k,410) &
                      *y(k,254) + .150_r8*rxt(k,385)*y(k,269)
         mat(k,1403) = mat(k,1403) + .530_r8*rxt(k,399)*y(k,125) + .530_r8*rxt(k,398) &
                      *y(k,127) + .530_r8*rxt(k,395)*y(k,236) + .260_r8*rxt(k,396) &
                      *y(k,237)
         mat(k,1443) = rxt(k,408)*y(k,125) + rxt(k,409)*y(k,127) + rxt(k,405)*y(k,237) &
                      + .450_r8*rxt(k,406)*y(k,243) + 4.000_r8*rxt(k,407)*y(k,252)
         mat(k,745) = rxt(k,411)*y(k,125) + .200_r8*rxt(k,410)*y(k,243)
         mat(k,1844) = rxt(k,351)*y(k,45) + rxt(k,377)*y(k,49) + .500_r8*rxt(k,358) &
                      *y(k,51) + .650_r8*rxt(k,524)*y(k,217)
         mat(k,1264) = rxt(k,386)*y(k,125) + .300_r8*rxt(k,384)*y(k,237) &
                      + .150_r8*rxt(k,385)*y(k,243)
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
         mat(k,2193) = -(rxt(k,244)*y(k,59) + (4._r8*rxt(k,321) + 4._r8*rxt(k,322) &
                      ) * y(k,237) + rxt(k,323)*y(k,243) + rxt(k,324)*y(k,125) &
                      + rxt(k,343)*y(k,233) + rxt(k,354)*y(k,236) + rxt(k,371) &
                      *y(k,234) + rxt(k,384)*y(k,269) + rxt(k,396)*y(k,250) + rxt(k,405) &
                      *y(k,252) + rxt(k,429)*y(k,245) + rxt(k,434)*y(k,246) + rxt(k,443) &
                      *y(k,101) + rxt(k,454)*y(k,274) + rxt(k,508)*y(k,260) + rxt(k,513) &
                      *y(k,270) + rxt(k,518)*y(k,271))
         mat(k,2492) = -rxt(k,244)*y(k,237)
         mat(k,2379) = -rxt(k,323)*y(k,237)
         mat(k,1967) = -rxt(k,324)*y(k,237)
         mat(k,952) = -rxt(k,343)*y(k,237)
         mat(k,1481) = -rxt(k,354)*y(k,237)
         mat(k,987) = -rxt(k,371)*y(k,237)
         mat(k,1268) = -rxt(k,384)*y(k,237)
         mat(k,1408) = -rxt(k,396)*y(k,237)
         mat(k,1449) = -rxt(k,405)*y(k,237)
         mat(k,1354) = -rxt(k,429)*y(k,237)
         mat(k,1387) = -rxt(k,434)*y(k,237)
         mat(k,1307) = -rxt(k,443)*y(k,237)
         mat(k,1285) = -rxt(k,454)*y(k,237)
         mat(k,1158) = -rxt(k,508)*y(k,237)
         mat(k,1230) = -rxt(k,513)*y(k,237)
         mat(k,1113) = -rxt(k,518)*y(k,237)
         mat(k,1191) = .280_r8*rxt(k,370)*y(k,137)
         mat(k,753) = rxt(k,357)*y(k,265)
         mat(k,448) = .700_r8*rxt(k,326)*y(k,265)
         mat(k,1642) = rxt(k,238)*y(k,56) + rxt(k,294)*y(k,73) + rxt(k,333)*y(k,261) &
                      + rxt(k,327)*y(k,265)
         mat(k,2013) = rxt(k,238)*y(k,54)
         mat(k,975) = rxt(k,294)*y(k,54)
         mat(k,939) = .050_r8*rxt(k,448)*y(k,137)
         mat(k,1307) = mat(k,1307) + rxt(k,442)*y(k,236)
         mat(k,1967) = mat(k,1967) + rxt(k,356)*y(k,236) + .830_r8*rxt(k,474)*y(k,238) &
                      + .170_r8*rxt(k,480)*y(k,253)
         mat(k,2120) = .280_r8*rxt(k,370)*y(k,29) + .050_r8*rxt(k,448)*y(k,98)
         mat(k,1481) = mat(k,1481) + rxt(k,442)*y(k,101) + rxt(k,356)*y(k,125) &
                      + 4.000_r8*rxt(k,353)*y(k,236) + .900_r8*rxt(k,354)*y(k,237) &
                      + .450_r8*rxt(k,355)*y(k,243) + rxt(k,428)*y(k,245) + rxt(k,433) &
                      *y(k,246) + rxt(k,395)*y(k,250) + rxt(k,404)*y(k,252) &
                      + rxt(k,453)*y(k,274)
         mat(k,2193) = mat(k,2193) + .900_r8*rxt(k,354)*y(k,236)
         mat(k,829) = .830_r8*rxt(k,474)*y(k,125) + .330_r8*rxt(k,473)*y(k,243)
         mat(k,2379) = mat(k,2379) + .450_r8*rxt(k,355)*y(k,236) + .330_r8*rxt(k,473) &
                      *y(k,238) + .070_r8*rxt(k,479)*y(k,253)
         mat(k,1354) = mat(k,1354) + rxt(k,428)*y(k,236)
         mat(k,1387) = mat(k,1387) + rxt(k,433)*y(k,236)
         mat(k,1408) = mat(k,1408) + rxt(k,395)*y(k,236)
         mat(k,1449) = mat(k,1449) + rxt(k,404)*y(k,236)
         mat(k,967) = .170_r8*rxt(k,480)*y(k,125) + .070_r8*rxt(k,479)*y(k,243)
         mat(k,2236) = rxt(k,333)*y(k,54)
         mat(k,1860) = rxt(k,357)*y(k,50) + .700_r8*rxt(k,326)*y(k,53) + rxt(k,327) &
                      *y(k,54)
         mat(k,1285) = mat(k,1285) + rxt(k,453)*y(k,236)
         mat(k,823) = -(rxt(k,473)*y(k,243) + rxt(k,474)*y(k,125) + rxt(k,475) &
                      *y(k,126))
         mat(k,2331) = -rxt(k,473)*y(k,238)
         mat(k,1918) = -rxt(k,474)*y(k,238)
         mat(k,1661) = -rxt(k,475)*y(k,238)
         mat(k,910) = -(rxt(k,606)*y(k,258) + rxt(k,607)*y(k,264) + rxt(k,608) &
                      *y(k,257))
         mat(k,891) = -rxt(k,606)*y(k,239)
         mat(k,899) = -rxt(k,607)*y(k,239)
         mat(k,677) = -rxt(k,608)*y(k,239)
         mat(k,609) = -((rxt(k,392) + rxt(k,393)) * y(k,125))
         mat(k,1907) = -(rxt(k,392) + rxt(k,393)) * y(k,240)
         mat(k,404) = rxt(k,391)*y(k,265)
         mat(k,1780) = rxt(k,391)*y(k,16)
         mat(k,511) = -(rxt(k,362)*y(k,136))
         mat(k,1571) = -rxt(k,362)*y(k,241)
         mat(k,1900) = .750_r8*rxt(k,360)*y(k,242)
         mat(k,841) = .750_r8*rxt(k,360)*y(k,125)
         mat(k,842) = -(rxt(k,359)*y(k,243) + rxt(k,360)*y(k,125))
         mat(k,2333) = -rxt(k,359)*y(k,242)
         mat(k,1919) = -rxt(k,360)*y(k,242)
         mat(k,578) = rxt(k,366)*y(k,265)
         mat(k,1802) = rxt(k,366)*y(k,25)
         mat(k,2382) = -((rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,76) + rxt(k,196) &
                      *y(k,135) + rxt(k,197)*y(k,137) + rxt(k,201)*y(k,265) &
                      + 4._r8*rxt(k,206)*y(k,243) + rxt(k,218)*y(k,127) + rxt(k,223) &
                      *y(k,125) + rxt(k,228)*y(k,126) + (rxt(k,241) + rxt(k,242) &
                      ) * y(k,56) + rxt(k,248)*y(k,59) + rxt(k,274)*y(k,17) + rxt(k,280) &
                      *y(k,19) + rxt(k,317)*y(k,42) + rxt(k,323)*y(k,237) + rxt(k,330) &
                      *y(k,244) + rxt(k,344)*y(k,233) + rxt(k,355)*y(k,236) + rxt(k,359) &
                      *y(k,242) + rxt(k,372)*y(k,234) + rxt(k,381)*y(k,268) + rxt(k,385) &
                      *y(k,269) + rxt(k,397)*y(k,250) + rxt(k,406)*y(k,252) + rxt(k,410) &
                      *y(k,254) + rxt(k,420)*y(k,228) + rxt(k,430)*y(k,245) + rxt(k,435) &
                      *y(k,246) + rxt(k,444)*y(k,101) + rxt(k,455)*y(k,274) + rxt(k,459) &
                      *y(k,227) + rxt(k,462)*y(k,230) + rxt(k,466)*y(k,232) + rxt(k,469) &
                      *y(k,235) + rxt(k,473)*y(k,238) + rxt(k,476)*y(k,251) + rxt(k,479) &
                      *y(k,253) + rxt(k,482)*y(k,267) + rxt(k,489)*y(k,272) + rxt(k,495) &
                      *y(k,275) + rxt(k,498)*y(k,277) + rxt(k,509)*y(k,260) + rxt(k,514) &
                      *y(k,270) + rxt(k,519)*y(k,271))
         mat(k,2144) = -(rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,243)
         mat(k,2059) = -rxt(k,196)*y(k,243)
         mat(k,2123) = -rxt(k,197)*y(k,243)
         mat(k,1863) = -rxt(k,201)*y(k,243)
         mat(k,2442) = -rxt(k,218)*y(k,243)
         mat(k,1970) = -rxt(k,223)*y(k,243)
         mat(k,1690) = -rxt(k,228)*y(k,243)
         mat(k,2016) = -(rxt(k,241) + rxt(k,242)) * y(k,243)
         mat(k,2495) = -rxt(k,248)*y(k,243)
         mat(k,1566) = -rxt(k,274)*y(k,243)
         mat(k,1620) = -rxt(k,280)*y(k,243)
         mat(k,2468) = -rxt(k,317)*y(k,243)
         mat(k,2196) = -rxt(k,323)*y(k,243)
         mat(k,487) = -rxt(k,330)*y(k,243)
         mat(k,953) = -rxt(k,344)*y(k,243)
         mat(k,1483) = -rxt(k,355)*y(k,243)
         mat(k,848) = -rxt(k,359)*y(k,243)
         mat(k,988) = -rxt(k,372)*y(k,243)
         mat(k,866) = -rxt(k,381)*y(k,243)
         mat(k,1269) = -rxt(k,385)*y(k,243)
         mat(k,1410) = -rxt(k,397)*y(k,243)
         mat(k,1451) = -rxt(k,406)*y(k,243)
         mat(k,749) = -rxt(k,410)*y(k,243)
         mat(k,1096) = -rxt(k,420)*y(k,243)
         mat(k,1356) = -rxt(k,430)*y(k,243)
         mat(k,1389) = -rxt(k,435)*y(k,243)
         mat(k,1309) = -rxt(k,444)*y(k,243)
         mat(k,1286) = -rxt(k,455)*y(k,243)
         mat(k,576) = -rxt(k,459)*y(k,243)
         mat(k,538) = -rxt(k,462)*y(k,243)
         mat(k,482) = -rxt(k,466)*y(k,243)
         mat(k,703) = -rxt(k,469)*y(k,243)
         mat(k,830) = -rxt(k,473)*y(k,243)
         mat(k,782) = -rxt(k,476)*y(k,243)
         mat(k,968) = -rxt(k,479)*y(k,243)
         mat(k,507) = -rxt(k,482)*y(k,243)
         mat(k,797) = -rxt(k,489)*y(k,243)
         mat(k,822) = -rxt(k,495)*y(k,243)
         mat(k,564) = -rxt(k,498)*y(k,243)
         mat(k,1160) = -rxt(k,509)*y(k,243)
         mat(k,1232) = -rxt(k,514)*y(k,243)
         mat(k,1115) = -rxt(k,519)*y(k,243)
         mat(k,1078) = .570_r8*rxt(k,503)*y(k,137)
         mat(k,215) = .650_r8*rxt(k,461)*y(k,265)
         mat(k,1566) = mat(k,1566) + rxt(k,273)*y(k,42)
         mat(k,1620) = mat(k,1620) + rxt(k,285)*y(k,265)
         mat(k,326) = .350_r8*rxt(k,339)*y(k,265)
         mat(k,583) = .130_r8*rxt(k,341)*y(k,137)
         mat(k,307) = rxt(k,346)*y(k,265)
         mat(k,1192) = .280_r8*rxt(k,370)*y(k,137)
         mat(k,2468) = mat(k,2468) + rxt(k,273)*y(k,17) + rxt(k,237)*y(k,56) &
                      + rxt(k,318)*y(k,127) + rxt(k,319)*y(k,135)
         mat(k,688) = rxt(k,302)*y(k,56) + rxt(k,303)*y(k,265)
         mat(k,417) = rxt(k,305)*y(k,56) + rxt(k,306)*y(k,265)
         mat(k,151) = rxt(k,352)*y(k,265)
         mat(k,871) = rxt(k,325)*y(k,265)
         mat(k,1645) = rxt(k,334)*y(k,261)
         mat(k,2016) = mat(k,2016) + rxt(k,237)*y(k,42) + rxt(k,302)*y(k,43) &
                      + rxt(k,305)*y(k,46) + rxt(k,240)*y(k,79)
         mat(k,2495) = mat(k,2495) + rxt(k,244)*y(k,237) + rxt(k,255)*y(k,265)
         mat(k,1208) = rxt(k,337)*y(k,265)
         mat(k,242) = .730_r8*rxt(k,472)*y(k,265)
         mat(k,349) = .500_r8*rxt(k,540)*y(k,265)
         mat(k,1214) = rxt(k,363)*y(k,265)
         mat(k,1103) = rxt(k,364)*y(k,265)
         mat(k,2144) = mat(k,2144) + rxt(k,195)*y(k,136)
         mat(k,648) = rxt(k,240)*y(k,56) + rxt(k,191)*y(k,135) + rxt(k,200)*y(k,265)
         mat(k,253) = rxt(k,328)*y(k,265)
         mat(k,958) = rxt(k,329)*y(k,265)
         mat(k,1249) = rxt(k,394)*y(k,265)
         mat(k,1257) = rxt(k,379)*y(k,265)
         mat(k,940) = .370_r8*rxt(k,448)*y(k,137)
         mat(k,641) = .300_r8*rxt(k,439)*y(k,265)
         mat(k,600) = rxt(k,440)*y(k,265)
         mat(k,1309) = mat(k,1309) + rxt(k,445)*y(k,125) + rxt(k,446)*y(k,127) &
                      + rxt(k,442)*y(k,236) + 1.200_r8*rxt(k,443)*y(k,237)
         mat(k,444) = rxt(k,447)*y(k,265)
         mat(k,1322) = .140_r8*rxt(k,401)*y(k,137)
         mat(k,392) = .200_r8*rxt(k,403)*y(k,265)
         mat(k,631) = .500_r8*rxt(k,414)*y(k,265)
         mat(k,1022) = .570_r8*rxt(k,506)*y(k,137)
         mat(k,1433) = .280_r8*rxt(k,415)*y(k,137)
         mat(k,468) = rxt(k,451)*y(k,265)
         mat(k,1145) = rxt(k,452)*y(k,265)
         mat(k,1970) = mat(k,1970) + rxt(k,445)*y(k,101) + rxt(k,421)*y(k,228) &
                      + rxt(k,463)*y(k,230) + rxt(k,468)*y(k,232) + rxt(k,345) &
                      *y(k,233) + rxt(k,373)*y(k,234) + rxt(k,324)*y(k,237) &
                      + .170_r8*rxt(k,474)*y(k,238) + rxt(k,392)*y(k,240) &
                      + .250_r8*rxt(k,360)*y(k,242) + rxt(k,332)*y(k,244) &
                      + .920_r8*rxt(k,431)*y(k,245) + .920_r8*rxt(k,437)*y(k,246) &
                      + .470_r8*rxt(k,399)*y(k,250) + .400_r8*rxt(k,477)*y(k,251) &
                      + .830_r8*rxt(k,480)*y(k,253) + rxt(k,483)*y(k,267) + rxt(k,382) &
                      *y(k,268) + .900_r8*rxt(k,515)*y(k,270) + .800_r8*rxt(k,520) &
                      *y(k,271) + rxt(k,490)*y(k,272) + rxt(k,456)*y(k,274) &
                      + rxt(k,496)*y(k,275) + rxt(k,499)*y(k,277)
         mat(k,2442) = mat(k,2442) + rxt(k,318)*y(k,42) + rxt(k,446)*y(k,101) &
                      + rxt(k,432)*y(k,245) + rxt(k,438)*y(k,246) + .470_r8*rxt(k,398) &
                      *y(k,250) + rxt(k,221)*y(k,265) + rxt(k,457)*y(k,274)
         mat(k,2059) = mat(k,2059) + rxt(k,319)*y(k,42) + rxt(k,191)*y(k,79)
         mat(k,1597) = rxt(k,195)*y(k,76) + rxt(k,362)*y(k,241)
         mat(k,2123) = mat(k,2123) + .570_r8*rxt(k,503)*y(k,6) + .130_r8*rxt(k,341) &
                      *y(k,25) + .280_r8*rxt(k,370)*y(k,29) + .370_r8*rxt(k,448) &
                      *y(k,98) + .140_r8*rxt(k,401)*y(k,106) + .570_r8*rxt(k,506) &
                      *y(k,111) + .280_r8*rxt(k,415)*y(k,112) + rxt(k,203)*y(k,265)
         mat(k,224) = .800_r8*rxt(k,484)*y(k,265)
         mat(k,996) = rxt(k,530)*y(k,265)
         mat(k,1171) = .200_r8*rxt(k,524)*y(k,265)
         mat(k,237) = .280_r8*rxt(k,492)*y(k,265)
         mat(k,263) = .380_r8*rxt(k,494)*y(k,265)
         mat(k,268) = .630_r8*rxt(k,500)*y(k,265)
         mat(k,1096) = mat(k,1096) + rxt(k,421)*y(k,125)
         mat(k,538) = mat(k,538) + rxt(k,463)*y(k,125)
         mat(k,482) = mat(k,482) + rxt(k,468)*y(k,125)
         mat(k,953) = mat(k,953) + rxt(k,345)*y(k,125) + 2.400_r8*rxt(k,342)*y(k,233) &
                      + rxt(k,343)*y(k,237)
         mat(k,988) = mat(k,988) + rxt(k,373)*y(k,125) + rxt(k,371)*y(k,237)
         mat(k,1483) = mat(k,1483) + rxt(k,442)*y(k,101) + .900_r8*rxt(k,354)*y(k,237) &
                      + rxt(k,428)*y(k,245) + rxt(k,433)*y(k,246) + .470_r8*rxt(k,395) &
                      *y(k,250) + rxt(k,453)*y(k,274)
         mat(k,2196) = mat(k,2196) + rxt(k,244)*y(k,59) + 1.200_r8*rxt(k,443)*y(k,101) &
                      + rxt(k,324)*y(k,125) + rxt(k,343)*y(k,233) + rxt(k,371) &
                      *y(k,234) + .900_r8*rxt(k,354)*y(k,236) + 4.000_r8*rxt(k,321) &
                      *y(k,237) + rxt(k,429)*y(k,245) + rxt(k,434)*y(k,246) &
                      + .730_r8*rxt(k,396)*y(k,250) + rxt(k,405)*y(k,252) &
                      + .500_r8*rxt(k,508)*y(k,260) + .300_r8*rxt(k,384)*y(k,269) &
                      + rxt(k,513)*y(k,270) + rxt(k,518)*y(k,271) + .800_r8*rxt(k,454) &
                      *y(k,274)
         mat(k,830) = mat(k,830) + .170_r8*rxt(k,474)*y(k,125) + .070_r8*rxt(k,473) &
                      *y(k,243)
         mat(k,615) = rxt(k,392)*y(k,125)
         mat(k,514) = rxt(k,362)*y(k,136)
         mat(k,848) = mat(k,848) + .250_r8*rxt(k,360)*y(k,125)
         mat(k,2382) = mat(k,2382) + .070_r8*rxt(k,473)*y(k,238) + .160_r8*rxt(k,476) &
                      *y(k,251) + .330_r8*rxt(k,479)*y(k,253)
         mat(k,487) = mat(k,487) + rxt(k,332)*y(k,125)
         mat(k,1356) = mat(k,1356) + .920_r8*rxt(k,431)*y(k,125) + rxt(k,432)*y(k,127) &
                      + rxt(k,428)*y(k,236) + rxt(k,429)*y(k,237)
         mat(k,1389) = mat(k,1389) + .920_r8*rxt(k,437)*y(k,125) + rxt(k,438)*y(k,127) &
                      + rxt(k,433)*y(k,236) + rxt(k,434)*y(k,237)
         mat(k,1410) = mat(k,1410) + .470_r8*rxt(k,399)*y(k,125) + .470_r8*rxt(k,398) &
                      *y(k,127) + .470_r8*rxt(k,395)*y(k,236) + .730_r8*rxt(k,396) &
                      *y(k,237)
         mat(k,782) = mat(k,782) + .400_r8*rxt(k,477)*y(k,125) + .160_r8*rxt(k,476) &
                      *y(k,243)
         mat(k,1451) = mat(k,1451) + rxt(k,405)*y(k,237)
         mat(k,968) = mat(k,968) + .830_r8*rxt(k,480)*y(k,125) + .330_r8*rxt(k,479) &
                      *y(k,243)
         mat(k,1160) = mat(k,1160) + .500_r8*rxt(k,508)*y(k,237)
         mat(k,2239) = rxt(k,334)*y(k,54)
         mat(k,1863) = mat(k,1863) + .650_r8*rxt(k,461)*y(k,7) + rxt(k,285)*y(k,19) &
                      + .350_r8*rxt(k,339)*y(k,24) + rxt(k,346)*y(k,26) + rxt(k,303) &
                      *y(k,43) + rxt(k,306)*y(k,46) + rxt(k,352)*y(k,47) + rxt(k,325) &
                      *y(k,52) + rxt(k,255)*y(k,59) + rxt(k,337)*y(k,62) &
                      + .730_r8*rxt(k,472)*y(k,66) + .500_r8*rxt(k,540)*y(k,67) &
                      + rxt(k,363)*y(k,74) + rxt(k,364)*y(k,75) + rxt(k,200)*y(k,79) &
                      + rxt(k,328)*y(k,86) + rxt(k,329)*y(k,87) + rxt(k,394)*y(k,93) &
                      + rxt(k,379)*y(k,95) + .300_r8*rxt(k,439)*y(k,99) + rxt(k,440) &
                      *y(k,100) + rxt(k,447)*y(k,102) + .200_r8*rxt(k,403)*y(k,107) &
                      + .500_r8*rxt(k,414)*y(k,110) + rxt(k,451)*y(k,116) + rxt(k,452) &
                      *y(k,117) + rxt(k,221)*y(k,127) + rxt(k,203)*y(k,137) &
                      + .800_r8*rxt(k,484)*y(k,145) + rxt(k,530)*y(k,158) &
                      + .200_r8*rxt(k,524)*y(k,217) + .280_r8*rxt(k,492)*y(k,219) &
                      + .380_r8*rxt(k,494)*y(k,221) + .630_r8*rxt(k,500)*y(k,223)
         mat(k,507) = mat(k,507) + rxt(k,483)*y(k,125)
         mat(k,866) = mat(k,866) + rxt(k,382)*y(k,125)
         mat(k,1269) = mat(k,1269) + .300_r8*rxt(k,384)*y(k,237)
         mat(k,1232) = mat(k,1232) + .900_r8*rxt(k,515)*y(k,125) + rxt(k,513)*y(k,237)
         mat(k,1115) = mat(k,1115) + .800_r8*rxt(k,520)*y(k,125) + rxt(k,518)*y(k,237)
         mat(k,797) = mat(k,797) + rxt(k,490)*y(k,125)
         mat(k,1286) = mat(k,1286) + rxt(k,456)*y(k,125) + rxt(k,457)*y(k,127) &
                      + rxt(k,453)*y(k,236) + .800_r8*rxt(k,454)*y(k,237)
         mat(k,822) = mat(k,822) + rxt(k,496)*y(k,125)
         mat(k,564) = mat(k,564) + rxt(k,499)*y(k,125)
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
         mat(k,483) = -(rxt(k,330)*y(k,243) + rxt(k,332)*y(k,125))
         mat(k,2305) = -rxt(k,330)*y(k,244)
         mat(k,1897) = -rxt(k,332)*y(k,244)
         mat(k,2447) = rxt(k,317)*y(k,243)
         mat(k,2305) = mat(k,2305) + rxt(k,317)*y(k,42)
         mat(k,1344) = -(rxt(k,428)*y(k,236) + rxt(k,429)*y(k,237) + rxt(k,430) &
                      *y(k,243) + rxt(k,431)*y(k,125) + rxt(k,432)*y(k,127))
         mat(k,1468) = -rxt(k,428)*y(k,245)
         mat(k,2176) = -rxt(k,429)*y(k,245)
         mat(k,2359) = -rxt(k,430)*y(k,245)
         mat(k,1949) = -rxt(k,431)*y(k,245)
         mat(k,2420) = -rxt(k,432)*y(k,245)
         mat(k,933) = .600_r8*rxt(k,449)*y(k,265)
         mat(k,1839) = .600_r8*rxt(k,449)*y(k,98)
         mat(k,1377) = -(rxt(k,433)*y(k,236) + rxt(k,434)*y(k,237) + rxt(k,435) &
                      *y(k,243) + rxt(k,437)*y(k,125) + rxt(k,438)*y(k,127))
         mat(k,1469) = -rxt(k,433)*y(k,246)
         mat(k,2177) = -rxt(k,434)*y(k,246)
         mat(k,2360) = -rxt(k,435)*y(k,246)
         mat(k,1950) = -rxt(k,437)*y(k,246)
         mat(k,2421) = -rxt(k,438)*y(k,246)
         mat(k,934) = .400_r8*rxt(k,449)*y(k,265)
         mat(k,1840) = .400_r8*rxt(k,449)*y(k,98)
         mat(k,108) = -(rxt(k,565)*y(k,243) + rxt(k,566)*y(k,125))
         mat(k,2284) = -rxt(k,565)*y(k,247)
         mat(k,1885) = -rxt(k,566)*y(k,247)
         mat(k,926) = rxt(k,568)*y(k,265)
         mat(k,1709) = rxt(k,568)*y(k,98)
         mat(k,114) = -(rxt(k,569)*y(k,243) + rxt(k,570)*y(k,125))
         mat(k,2285) = -rxt(k,569)*y(k,248)
         mat(k,1886) = -rxt(k,570)*y(k,248)
         mat(k,115) = rxt(k,571)*y(k,265)
         mat(k,1710) = rxt(k,571)*y(k,104)
         mat(k,122) = -(rxt(k,572)*y(k,243) + rxt(k,573)*y(k,125))
         mat(k,2286) = -rxt(k,572)*y(k,249)
         mat(k,1887) = -rxt(k,573)*y(k,249)
         mat(k,123) = rxt(k,574)*y(k,265)
         mat(k,1712) = rxt(k,574)*y(k,105)
         mat(k,1401) = -(rxt(k,395)*y(k,236) + rxt(k,396)*y(k,237) + rxt(k,397) &
                      *y(k,243) + rxt(k,398)*y(k,127) + (rxt(k,399) + rxt(k,400) &
                      ) * y(k,125))
         mat(k,1470) = -rxt(k,395)*y(k,250)
         mat(k,2178) = -rxt(k,396)*y(k,250)
         mat(k,2361) = -rxt(k,397)*y(k,250)
         mat(k,2422) = -rxt(k,398)*y(k,250)
         mat(k,1951) = -(rxt(k,399) + rxt(k,400)) * y(k,250)
         mat(k,1316) = .500_r8*rxt(k,402)*y(k,265)
         mat(k,389) = .200_r8*rxt(k,403)*y(k,265)
         mat(k,1420) = rxt(k,416)*y(k,265)
         mat(k,1841) = .500_r8*rxt(k,402)*y(k,106) + .200_r8*rxt(k,403)*y(k,107) &
                      + rxt(k,416)*y(k,112)
         mat(k,777) = -(rxt(k,476)*y(k,243) + rxt(k,477)*y(k,125) + rxt(k,478) &
                      *y(k,126))
         mat(k,2327) = -rxt(k,476)*y(k,251)
         mat(k,1915) = -rxt(k,477)*y(k,251)
         mat(k,1660) = -rxt(k,478)*y(k,251)
         mat(k,1442) = -(rxt(k,404)*y(k,236) + rxt(k,405)*y(k,237) + rxt(k,406) &
                      *y(k,243) + 4._r8*rxt(k,407)*y(k,252) + rxt(k,408)*y(k,125) &
                      + rxt(k,409)*y(k,127) + rxt(k,417)*y(k,126))
         mat(k,1472) = -rxt(k,404)*y(k,252)
         mat(k,2180) = -rxt(k,405)*y(k,252)
         mat(k,2363) = -rxt(k,406)*y(k,252)
         mat(k,1953) = -rxt(k,408)*y(k,252)
         mat(k,2424) = -rxt(k,409)*y(k,252)
         mat(k,1671) = -rxt(k,417)*y(k,252)
         mat(k,1317) = .500_r8*rxt(k,402)*y(k,265)
         mat(k,390) = .500_r8*rxt(k,403)*y(k,265)
         mat(k,1843) = .500_r8*rxt(k,402)*y(k,106) + .500_r8*rxt(k,403)*y(k,107)
         mat(k,960) = -(rxt(k,479)*y(k,243) + rxt(k,480)*y(k,125) + rxt(k,481) &
                      *y(k,126))
         mat(k,2340) = -rxt(k,479)*y(k,253)
         mat(k,1927) = -rxt(k,480)*y(k,253)
         mat(k,1663) = -rxt(k,481)*y(k,253)
         mat(k,743) = -(rxt(k,410)*y(k,243) + rxt(k,411)*y(k,125))
         mat(k,2324) = -rxt(k,410)*y(k,254)
         mat(k,1914) = -rxt(k,411)*y(k,254)
         mat(k,566) = rxt(k,412)*y(k,265)
         mat(k,384) = rxt(k,413)*y(k,265)
         mat(k,1793) = rxt(k,412)*y(k,108) + rxt(k,413)*y(k,109)
         mat(k,130) = -(rxt(k,576)*y(k,243) + rxt(k,577)*y(k,125))
         mat(k,2287) = -rxt(k,576)*y(k,255)
         mat(k,1888) = -rxt(k,577)*y(k,255)
         mat(k,1003) = rxt(k,579)*y(k,265)
         mat(k,1714) = rxt(k,579)*y(k,111)
         mat(k,704) = -(rxt(k,209)*y(k,125) + rxt(k,210)*y(k,135) + rxt(k,233) &
                      *y(k,239) + rxt(k,234)*y(k,136))
         mat(k,1912) = -rxt(k,209)*y(k,256)
         mat(k,2028) = -rxt(k,210)*y(k,256)
         mat(k,906) = -rxt(k,233)*y(k,256)
         mat(k,1574) = -rxt(k,234)*y(k,256)
         mat(k,2028) = mat(k,2028) + rxt(k,610)*y(k,257)
         mat(k,906) = mat(k,906) + .900_r8*rxt(k,608)*y(k,257) + .800_r8*rxt(k,606) &
                      *y(k,258)
         mat(k,673) = rxt(k,610)*y(k,135) + .900_r8*rxt(k,608)*y(k,239)
         mat(k,889) = .800_r8*rxt(k,606)*y(k,239)
         mat(k,672) = -(rxt(k,608)*y(k,239) + rxt(k,609)*y(k,136) + (rxt(k,610) &
                      + rxt(k,611)) * y(k,135))
         mat(k,905) = -rxt(k,608)*y(k,257)
         mat(k,1573) = -rxt(k,609)*y(k,257)
         mat(k,2027) = -(rxt(k,610) + rxt(k,611)) * y(k,257)
         mat(k,890) = -(rxt(k,606)*y(k,239))
         mat(k,908) = -rxt(k,606)*y(k,258)
         mat(k,1029) = rxt(k,615)*y(k,264)
         mat(k,1921) = rxt(k,617)*y(k,264)
         mat(k,2034) = rxt(k,610)*y(k,257)
         mat(k,1577) = rxt(k,614)*y(k,259)
         mat(k,675) = rxt(k,610)*y(k,135)
         mat(k,545) = rxt(k,614)*y(k,136)
         mat(k,897) = rxt(k,615)*y(k,113) + rxt(k,617)*y(k,125)
         mat(k,543) = -(rxt(k,612)*y(k,135) + (rxt(k,613) + rxt(k,614)) * y(k,136))
         mat(k,2024) = -rxt(k,612)*y(k,259)
         mat(k,1572) = -(rxt(k,613) + rxt(k,614)) * y(k,259)
         mat(k,1151) = -(rxt(k,508)*y(k,237) + rxt(k,509)*y(k,243) + rxt(k,510) &
                      *y(k,125) + rxt(k,511)*y(k,127))
         mat(k,2164) = -rxt(k,508)*y(k,260)
         mat(k,2347) = -rxt(k,509)*y(k,260)
         mat(k,1936) = -rxt(k,510)*y(k,260)
         mat(k,2406) = -rxt(k,511)*y(k,260)
         mat(k,1069) = rxt(k,502)*y(k,127)
         mat(k,1013) = rxt(k,505)*y(k,127)
         mat(k,2406) = mat(k,2406) + rxt(k,502)*y(k,6) + rxt(k,505)*y(k,111) &
                      + .500_r8*rxt(k,522)*y(k,216)
         mat(k,491) = rxt(k,512)*y(k,265)
         mat(k,1118) = .500_r8*rxt(k,522)*y(k,127)
         mat(k,1825) = rxt(k,512)*y(k,129)
         mat(k,2237) = -(rxt(k,173)*y(k,77) + rxt(k,174)*y(k,278) + (rxt(k,176) &
                      + rxt(k,177)) * y(k,136) + rxt(k,178)*y(k,137) + (rxt(k,226) &
                      + rxt(k,227)) * y(k,114) + rxt(k,262)*y(k,33) + rxt(k,263) &
                      *y(k,34) + rxt(k,264)*y(k,36) + rxt(k,265)*y(k,37) + rxt(k,266) &
                      *y(k,38) + rxt(k,267)*y(k,39) + rxt(k,268)*y(k,40) + (rxt(k,269) &
                      + rxt(k,270)) * y(k,85) + rxt(k,289)*y(k,35) + rxt(k,290) &
                      *y(k,55) + rxt(k,291)*y(k,78) + (rxt(k,292) + rxt(k,293) &
                      ) * y(k,81) + rxt(k,298)*y(k,64) + rxt(k,299)*y(k,65) + rxt(k,312) &
                      *y(k,41) + rxt(k,313)*y(k,43) + rxt(k,314)*y(k,82) + rxt(k,315) &
                      *y(k,83) + rxt(k,316)*y(k,84) + (rxt(k,333) + rxt(k,334) &
                      + rxt(k,335)) * y(k,54) + rxt(k,336)*y(k,86))
         mat(k,1514) = -rxt(k,173)*y(k,261)
         mat(k,2520) = -rxt(k,174)*y(k,261)
         mat(k,1596) = -(rxt(k,176) + rxt(k,177)) * y(k,261)
         mat(k,2121) = -rxt(k,178)*y(k,261)
         mat(k,320) = -(rxt(k,226) + rxt(k,227)) * y(k,261)
         mat(k,148) = -rxt(k,262)*y(k,261)
         mat(k,182) = -rxt(k,263)*y(k,261)
         mat(k,158) = -rxt(k,264)*y(k,261)
         mat(k,192) = -rxt(k,265)*y(k,261)
         mat(k,162) = -rxt(k,266)*y(k,261)
         mat(k,197) = -rxt(k,267)*y(k,261)
         mat(k,166) = -rxt(k,268)*y(k,261)
         mat(k,1548) = -(rxt(k,269) + rxt(k,270)) * y(k,261)
         mat(k,187) = -rxt(k,289)*y(k,261)
         mat(k,456) = -rxt(k,290)*y(k,261)
         mat(k,174) = -rxt(k,291)*y(k,261)
         mat(k,879) = -(rxt(k,292) + rxt(k,293)) * y(k,261)
         mat(k,288) = -rxt(k,298)*y(k,261)
         mat(k,297) = -rxt(k,299)*y(k,261)
         mat(k,521) = -rxt(k,312)*y(k,261)
         mat(k,687) = -rxt(k,313)*y(k,261)
         mat(k,293) = -rxt(k,314)*y(k,261)
         mat(k,303) = -rxt(k,315)*y(k,261)
         mat(k,359) = -rxt(k,316)*y(k,261)
         mat(k,1643) = -(rxt(k,333) + rxt(k,334) + rxt(k,335)) * y(k,261)
         mat(k,252) = -rxt(k,336)*y(k,261)
         mat(k,1596) = mat(k,1596) + rxt(k,234)*y(k,256)
         mat(k,917) = .850_r8*rxt(k,607)*y(k,264)
         mat(k,709) = rxt(k,234)*y(k,136)
         mat(k,904) = .850_r8*rxt(k,607)*y(k,239)
         mat(k,225) = -(rxt(k,180)*y(k,135) + rxt(k,181)*y(k,136))
         mat(k,2021) = -rxt(k,180)*y(k,262)
         mat(k,1569) = -rxt(k,181)*y(k,262)
         mat(k,1487) = rxt(k,182)*y(k,263)
         mat(k,2021) = mat(k,2021) + rxt(k,184)*y(k,263)
         mat(k,1569) = mat(k,1569) + rxt(k,185)*y(k,263)
         mat(k,2069) = rxt(k,186)*y(k,263)
         mat(k,227) = rxt(k,182)*y(k,63) + rxt(k,184)*y(k,135) + rxt(k,185)*y(k,136) &
                      + rxt(k,186)*y(k,137)
         mat(k,228) = -(rxt(k,182)*y(k,63) + rxt(k,184)*y(k,135) + rxt(k,185)*y(k,136) &
                      + rxt(k,186)*y(k,137))
         mat(k,1488) = -rxt(k,182)*y(k,263)
         mat(k,2022) = -rxt(k,184)*y(k,263)
         mat(k,1570) = -rxt(k,185)*y(k,263)
         mat(k,2070) = -rxt(k,186)*y(k,263)
         mat(k,1570) = mat(k,1570) + rxt(k,176)*y(k,261)
         mat(k,2210) = rxt(k,176)*y(k,136)
         mat(k,898) = -(rxt(k,607)*y(k,239) + rxt(k,615)*y(k,113) + rxt(k,617) &
                      *y(k,125))
         mat(k,909) = -rxt(k,607)*y(k,264)
         mat(k,1030) = -rxt(k,615)*y(k,264)
         mat(k,1922) = -rxt(k,617)*y(k,264)
         mat(k,1491) = rxt(k,618)*y(k,266)
         mat(k,1578) = rxt(k,609)*y(k,257) + rxt(k,613)*y(k,259) + rxt(k,620)*y(k,266)
         mat(k,676) = rxt(k,609)*y(k,136)
         mat(k,546) = rxt(k,613)*y(k,136)
         mat(k,852) = rxt(k,618)*y(k,63) + rxt(k,620)*y(k,136)
         mat(k,1854) = -(rxt(k,199)*y(k,77) + rxt(k,200)*y(k,79) + rxt(k,201)*y(k,243) &
                      + rxt(k,202)*y(k,135) + rxt(k,203)*y(k,137) + (4._r8*rxt(k,204) &
                      + 4._r8*rxt(k,205)) * y(k,265) + rxt(k,208)*y(k,90) + rxt(k,221) &
                      *y(k,127) + rxt(k,222)*y(k,113) + rxt(k,230)*y(k,126) + rxt(k,231) &
                      *y(k,89) + rxt(k,253)*y(k,60) + (rxt(k,255) + rxt(k,256) &
                      ) * y(k,59) + rxt(k,258)*y(k,85) + rxt(k,261)*y(k,92) + rxt(k,285) &
                      *y(k,19) + rxt(k,287)*y(k,81) + rxt(k,301)*y(k,41) + rxt(k,303) &
                      *y(k,43) + rxt(k,304)*y(k,44) + rxt(k,306)*y(k,46) + rxt(k,308) &
                      *y(k,55) + rxt(k,309)*y(k,82) + rxt(k,310)*y(k,83) + rxt(k,311) &
                      *y(k,84) + rxt(k,320)*y(k,42) + rxt(k,325)*y(k,52) + rxt(k,326) &
                      *y(k,53) + rxt(k,327)*y(k,54) + rxt(k,328)*y(k,86) + rxt(k,329) &
                      *y(k,87) + rxt(k,337)*y(k,62) + rxt(k,339)*y(k,24) + rxt(k,346) &
                      *y(k,26) + rxt(k,347)*y(k,27) + rxt(k,349)*y(k,28) + rxt(k,351) &
                      *y(k,45) + rxt(k,352)*y(k,47) + rxt(k,357)*y(k,50) + rxt(k,358) &
                      *y(k,51) + rxt(k,363)*y(k,74) + rxt(k,364)*y(k,75) + rxt(k,365) &
                      *y(k,142) + rxt(k,366)*y(k,25) + rxt(k,374)*y(k,30) + rxt(k,375) &
                      *y(k,31) + rxt(k,377)*y(k,49) + rxt(k,379)*y(k,95) + rxt(k,380) &
                      *y(k,128) + rxt(k,383)*y(k,153) + rxt(k,387)*y(k,154) + rxt(k,388) &
                      *y(k,29) + rxt(k,389)*y(k,48) + rxt(k,391)*y(k,16) + rxt(k,394) &
                      *y(k,93) + rxt(k,402)*y(k,106) + rxt(k,403)*y(k,107) + rxt(k,412) &
                      *y(k,108) + rxt(k,413)*y(k,109) + rxt(k,414)*y(k,110) + rxt(k,416) &
                      *y(k,112) + rxt(k,419)*y(k,1) + rxt(k,423)*y(k,2) + rxt(k,424) &
                      *y(k,15) + rxt(k,425)*y(k,94) + rxt(k,426)*y(k,96) + rxt(k,427) &
                      *y(k,97) + rxt(k,439)*y(k,99) + rxt(k,440)*y(k,100) + rxt(k,447) &
                      *y(k,102) + rxt(k,449)*y(k,98) + rxt(k,450)*y(k,103) + rxt(k,451) &
                      *y(k,116) + rxt(k,452)*y(k,117) + rxt(k,458)*y(k,220) + rxt(k,461) &
                      *y(k,7) + rxt(k,464)*y(k,8) + rxt(k,465)*y(k,22) + rxt(k,467) &
                      *y(k,23) + rxt(k,471)*y(k,32) + rxt(k,472)*y(k,66) + rxt(k,484) &
                      *y(k,145) + rxt(k,487)*y(k,146) + rxt(k,491)*y(k,218) + (rxt(k,492) &
                      + rxt(k,582)) * y(k,219) + rxt(k,494)*y(k,221) + rxt(k,497) &
                      *y(k,222) + rxt(k,500)*y(k,223) + rxt(k,501)*y(k,224) + rxt(k,504) &
                      *y(k,6) + rxt(k,507)*y(k,111) + rxt(k,512)*y(k,129) + rxt(k,516) &
                      *y(k,213) + rxt(k,517)*y(k,214) + rxt(k,521)*y(k,215) + rxt(k,523) &
                      *y(k,216) + rxt(k,524)*y(k,217) + (rxt(k,526) + rxt(k,540) &
                      ) * y(k,67) + rxt(k,528)*y(k,140) + rxt(k,530)*y(k,158) &
                      + rxt(k,534)*y(k,155) + rxt(k,539)*y(k,157) + rxt(k,542) &
                      *y(k,121))
         mat(k,1509) = -rxt(k,199)*y(k,265)
         mat(k,645) = -rxt(k,200)*y(k,265)
         mat(k,2373) = -rxt(k,201)*y(k,265)
         mat(k,2050) = -rxt(k,202)*y(k,265)
         mat(k,2114) = -rxt(k,203)*y(k,265)
         mat(k,553) = -rxt(k,208)*y(k,265)
         mat(k,2433) = -rxt(k,221)*y(k,265)
         mat(k,1037) = -rxt(k,222)*y(k,265)
         mat(k,1681) = -rxt(k,230)*y(k,265)
         mat(k,2253) = -rxt(k,231)*y(k,265)
         mat(k,1047) = -rxt(k,253)*y(k,265)
         mat(k,2486) = -(rxt(k,255) + rxt(k,256)) * y(k,265)
         mat(k,1544) = -rxt(k,258)*y(k,265)
         mat(k,884) = -rxt(k,261)*y(k,265)
         mat(k,1612) = -rxt(k,285)*y(k,265)
         mat(k,876) = -rxt(k,287)*y(k,265)
         mat(k,519) = -rxt(k,301)*y(k,265)
         mat(k,684) = -rxt(k,303)*y(k,265)
         mat(k,168) = -rxt(k,304)*y(k,265)
         mat(k,414) = -rxt(k,306)*y(k,265)
         mat(k,454) = -rxt(k,308)*y(k,265)
         mat(k,291) = -rxt(k,309)*y(k,265)
         mat(k,301) = -rxt(k,310)*y(k,265)
         mat(k,357) = -rxt(k,311)*y(k,265)
         mat(k,2459) = -rxt(k,320)*y(k,265)
         mat(k,870) = -rxt(k,325)*y(k,265)
         mat(k,446) = -rxt(k,326)*y(k,265)
         mat(k,1636) = -rxt(k,327)*y(k,265)
         mat(k,251) = -rxt(k,328)*y(k,265)
         mat(k,957) = -rxt(k,329)*y(k,265)
         mat(k,1207) = -rxt(k,337)*y(k,265)
         mat(k,325) = -rxt(k,339)*y(k,265)
         mat(k,306) = -rxt(k,346)*y(k,265)
         mat(k,378) = -rxt(k,347)*y(k,265)
         mat(k,330) = -rxt(k,349)*y(k,265)
         mat(k,1199) = -rxt(k,351)*y(k,265)
         mat(k,150) = -rxt(k,352)*y(k,265)
         mat(k,752) = -rxt(k,357)*y(k,265)
         mat(k,668) = -rxt(k,358)*y(k,265)
         mat(k,1213) = -rxt(k,363)*y(k,265)
         mat(k,1102) = -rxt(k,364)*y(k,265)
         mat(k,589) = -rxt(k,365)*y(k,265)
         mat(k,581) = -rxt(k,366)*y(k,265)
         mat(k,436) = -rxt(k,374)*y(k,265)
         mat(k,336) = -rxt(k,375)*y(k,265)
         mat(k,1329) = -rxt(k,377)*y(k,265)
         mat(k,1256) = -rxt(k,379)*y(k,265)
         mat(k,922) = -rxt(k,380)*y(k,265)
         mat(k,605) = -rxt(k,383)*y(k,265)
         mat(k,430) = -rxt(k,387)*y(k,265)
         mat(k,1188) = -rxt(k,388)*y(k,265)
         mat(k,1128) = -rxt(k,389)*y(k,265)
         mat(k,408) = -rxt(k,391)*y(k,265)
         mat(k,1246) = -rxt(k,394)*y(k,265)
         mat(k,1320) = -rxt(k,402)*y(k,265)
         mat(k,391) = -rxt(k,403)*y(k,265)
         mat(k,569) = -rxt(k,412)*y(k,265)
         mat(k,387) = -rxt(k,413)*y(k,265)
         mat(k,630) = -rxt(k,414)*y(k,265)
         mat(k,1428) = -rxt(k,416)*y(k,265)
         mat(k,739) = -rxt(k,419)*y(k,265)
         mat(k,718) = -rxt(k,423)*y(k,265)
         mat(k,271) = -rxt(k,424)*y(k,265)
         mat(k,284) = -rxt(k,425)*y(k,265)
         mat(k,382) = -rxt(k,426)*y(k,265)
         mat(k,200) = -rxt(k,427)*y(k,265)
         mat(k,639) = -rxt(k,439)*y(k,265)
         mat(k,598) = -rxt(k,440)*y(k,265)
         mat(k,443) = -rxt(k,447)*y(k,265)
         mat(k,937) = -rxt(k,449)*y(k,265)
         mat(k,803) = -rxt(k,450)*y(k,265)
         mat(k,466) = -rxt(k,451)*y(k,265)
         mat(k,1141) = -rxt(k,452)*y(k,265)
         mat(k,249) = -rxt(k,458)*y(k,265)
         mat(k,214) = -rxt(k,461)*y(k,265)
         mat(k,461) = -rxt(k,464)*y(k,265)
         mat(k,280) = -rxt(k,465)*y(k,265)
         mat(k,373) = -rxt(k,467)*y(k,265)
         mat(k,316) = -rxt(k,471)*y(k,265)
         mat(k,241) = -rxt(k,472)*y(k,265)
         mat(k,223) = -rxt(k,484)*y(k,265)
         mat(k,367) = -rxt(k,487)*y(k,265)
         mat(k,659) = -rxt(k,491)*y(k,265)
         mat(k,236) = -(rxt(k,492) + rxt(k,582)) * y(k,265)
         mat(k,262) = -rxt(k,494)*y(k,265)
         mat(k,775) = -rxt(k,497)*y(k,265)
         mat(k,267) = -rxt(k,500)*y(k,265)
         mat(k,473) = -rxt(k,501)*y(k,265)
         mat(k,1076) = -rxt(k,504)*y(k,265)
         mat(k,1020) = -rxt(k,507)*y(k,265)
         mat(k,493) = -rxt(k,512)*y(k,265)
         mat(k,728) = -rxt(k,516)*y(k,265)
         mat(k,693) = -rxt(k,517)*y(k,265)
         mat(k,527) = -rxt(k,521)*y(k,265)
         mat(k,1122) = -rxt(k,523)*y(k,265)
         mat(k,1169) = -rxt(k,524)*y(k,265)
         mat(k,347) = -(rxt(k,526) + rxt(k,540)) * y(k,265)
         mat(k,424) = -rxt(k,528)*y(k,265)
         mat(k,994) = -rxt(k,530)*y(k,265)
         mat(k,758) = -rxt(k,534)*y(k,265)
         mat(k,1526) = -rxt(k,539)*y(k,265)
         mat(k,144) = -rxt(k,542)*y(k,265)
         mat(k,1076) = mat(k,1076) + .630_r8*rxt(k,503)*y(k,137)
         mat(k,325) = mat(k,325) + .650_r8*rxt(k,339)*y(k,265)
         mat(k,581) = mat(k,581) + .130_r8*rxt(k,341)*y(k,137)
         mat(k,378) = mat(k,378) + .500_r8*rxt(k,347)*y(k,265)
         mat(k,1188) = mat(k,1188) + .360_r8*rxt(k,370)*y(k,137)
         mat(k,2459) = mat(k,2459) + rxt(k,319)*y(k,135)
         mat(k,446) = mat(k,446) + .300_r8*rxt(k,326)*y(k,265)
         mat(k,1636) = mat(k,1636) + rxt(k,333)*y(k,261)
         mat(k,2007) = rxt(k,242)*y(k,243)
         mat(k,973) = rxt(k,296)*y(k,278)
         mat(k,2135) = rxt(k,198)*y(k,137) + 2.000_r8*rxt(k,193)*y(k,243)
         mat(k,1509) = mat(k,1509) + rxt(k,190)*y(k,135) + rxt(k,173)*y(k,261)
         mat(k,645) = mat(k,645) + rxt(k,191)*y(k,135)
         mat(k,876) = mat(k,876) + rxt(k,286)*y(k,135) + rxt(k,292)*y(k,261)
         mat(k,1544) = mat(k,1544) + rxt(k,257)*y(k,135) + rxt(k,269)*y(k,261)
         mat(k,251) = mat(k,251) + rxt(k,336)*y(k,261)
         mat(k,836) = rxt(k,288)*y(k,135)
         mat(k,884) = mat(k,884) + rxt(k,260)*y(k,135)
         mat(k,937) = mat(k,937) + .320_r8*rxt(k,448)*y(k,137)
         mat(k,803) = mat(k,803) + .600_r8*rxt(k,450)*y(k,265)
         mat(k,1320) = mat(k,1320) + .240_r8*rxt(k,401)*y(k,137)
         mat(k,391) = mat(k,391) + .100_r8*rxt(k,403)*y(k,265)
         mat(k,1020) = mat(k,1020) + .630_r8*rxt(k,506)*y(k,137)
         mat(k,1428) = mat(k,1428) + .360_r8*rxt(k,415)*y(k,137)
         mat(k,1961) = rxt(k,223)*y(k,243)
         mat(k,2433) = mat(k,2433) + rxt(k,218)*y(k,243)
         mat(k,2050) = mat(k,2050) + rxt(k,319)*y(k,42) + rxt(k,190)*y(k,77) &
                      + rxt(k,191)*y(k,79) + rxt(k,286)*y(k,81) + rxt(k,257)*y(k,85) &
                      + rxt(k,288)*y(k,91) + rxt(k,260)*y(k,92) + rxt(k,196)*y(k,243)
         mat(k,2114) = mat(k,2114) + .630_r8*rxt(k,503)*y(k,6) + .130_r8*rxt(k,341) &
                      *y(k,25) + .360_r8*rxt(k,370)*y(k,29) + rxt(k,198)*y(k,76) &
                      + .320_r8*rxt(k,448)*y(k,98) + .240_r8*rxt(k,401)*y(k,106) &
                      + .630_r8*rxt(k,506)*y(k,111) + .360_r8*rxt(k,415)*y(k,112) &
                      + rxt(k,197)*y(k,243)
         mat(k,605) = mat(k,605) + .500_r8*rxt(k,383)*y(k,265)
         mat(k,249) = mat(k,249) + .500_r8*rxt(k,458)*y(k,265)
         mat(k,574) = .400_r8*rxt(k,459)*y(k,243)
         mat(k,1478) = .450_r8*rxt(k,355)*y(k,243)
         mat(k,827) = .400_r8*rxt(k,473)*y(k,243)
         mat(k,2373) = mat(k,2373) + rxt(k,242)*y(k,56) + 2.000_r8*rxt(k,193)*y(k,76) &
                      + rxt(k,223)*y(k,125) + rxt(k,218)*y(k,127) + rxt(k,196) &
                      *y(k,135) + rxt(k,197)*y(k,137) + .400_r8*rxt(k,459)*y(k,227) &
                      + .450_r8*rxt(k,355)*y(k,236) + .400_r8*rxt(k,473)*y(k,238) &
                      + .450_r8*rxt(k,406)*y(k,252) + .400_r8*rxt(k,479)*y(k,253) &
                      + .200_r8*rxt(k,410)*y(k,254) + .150_r8*rxt(k,385)*y(k,269)
         mat(k,1446) = .450_r8*rxt(k,406)*y(k,243)
         mat(k,965) = .400_r8*rxt(k,479)*y(k,243)
         mat(k,747) = .200_r8*rxt(k,410)*y(k,243)
         mat(k,2230) = rxt(k,333)*y(k,54) + rxt(k,173)*y(k,77) + rxt(k,292)*y(k,81) &
                      + rxt(k,269)*y(k,85) + rxt(k,336)*y(k,86) + 2.000_r8*rxt(k,174) &
                      *y(k,278)
         mat(k,1854) = mat(k,1854) + .650_r8*rxt(k,339)*y(k,24) + .500_r8*rxt(k,347) &
                      *y(k,27) + .300_r8*rxt(k,326)*y(k,53) + .600_r8*rxt(k,450) &
                      *y(k,103) + .100_r8*rxt(k,403)*y(k,107) + .500_r8*rxt(k,383) &
                      *y(k,153) + .500_r8*rxt(k,458)*y(k,220)
         mat(k,1266) = .150_r8*rxt(k,385)*y(k,243)
         mat(k,2513) = rxt(k,296)*y(k,73) + 2.000_r8*rxt(k,174)*y(k,261)
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
         mat(k,850) = -(rxt(k,618)*y(k,63) + rxt(k,620)*y(k,136))
         mat(k,1489) = -rxt(k,618)*y(k,266)
         mat(k,1576) = -rxt(k,620)*y(k,266)
         mat(k,2031) = rxt(k,611)*y(k,257) + rxt(k,612)*y(k,259)
         mat(k,674) = rxt(k,611)*y(k,135)
         mat(k,544) = rxt(k,612)*y(k,135)
         mat(k,502) = -(rxt(k,482)*y(k,243) + rxt(k,483)*y(k,125))
         mat(k,2307) = -rxt(k,482)*y(k,267)
         mat(k,1898) = -rxt(k,483)*y(k,267)
         mat(k,239) = .200_r8*rxt(k,472)*y(k,265)
         mat(k,221) = .140_r8*rxt(k,484)*y(k,265)
         mat(k,365) = rxt(k,487)*y(k,265)
         mat(k,1766) = .200_r8*rxt(k,472)*y(k,66) + .140_r8*rxt(k,484)*y(k,145) &
                      + rxt(k,487)*y(k,146)
         mat(k,859) = -(rxt(k,381)*y(k,243) + rxt(k,382)*y(k,125))
         mat(k,2334) = -rxt(k,381)*y(k,268)
         mat(k,1920) = -rxt(k,382)*y(k,268)
         mat(k,1175) = rxt(k,388)*y(k,265)
         mat(k,602) = .500_r8*rxt(k,383)*y(k,265)
         mat(k,1803) = rxt(k,388)*y(k,29) + .500_r8*rxt(k,383)*y(k,153)
         mat(k,1262) = -(rxt(k,384)*y(k,237) + rxt(k,385)*y(k,243) + rxt(k,386) &
                      *y(k,125))
         mat(k,2171) = -rxt(k,384)*y(k,269)
         mat(k,2354) = -rxt(k,385)*y(k,269)
         mat(k,1944) = -rxt(k,386)*y(k,269)
         mat(k,1072) = .060_r8*rxt(k,503)*y(k,137)
         mat(k,1126) = rxt(k,389)*y(k,265)
         mat(k,1016) = .060_r8*rxt(k,506)*y(k,137)
         mat(k,2097) = .060_r8*rxt(k,503)*y(k,6) + .060_r8*rxt(k,506)*y(k,111)
         mat(k,428) = rxt(k,387)*y(k,265)
         mat(k,1166) = .150_r8*rxt(k,524)*y(k,265)
         mat(k,1834) = rxt(k,389)*y(k,48) + rxt(k,387)*y(k,154) + .150_r8*rxt(k,524) &
                      *y(k,217)
         mat(k,1223) = -(rxt(k,513)*y(k,237) + rxt(k,514)*y(k,243) + rxt(k,515) &
                      *y(k,125))
         mat(k,2169) = -rxt(k,513)*y(k,270)
         mat(k,2352) = -rxt(k,514)*y(k,270)
         mat(k,1941) = -rxt(k,515)*y(k,270)
         mat(k,2412) = .500_r8*rxt(k,522)*y(k,216)
         mat(k,726) = rxt(k,516)*y(k,265)
         mat(k,1121) = .500_r8*rxt(k,522)*y(k,127) + rxt(k,523)*y(k,265)
         mat(k,1831) = rxt(k,516)*y(k,213) + rxt(k,523)*y(k,216)
         mat(k,1107) = -(rxt(k,518)*y(k,237) + rxt(k,519)*y(k,243) + rxt(k,520) &
                      *y(k,125))
         mat(k,2160) = -rxt(k,518)*y(k,271)
         mat(k,2344) = -rxt(k,519)*y(k,271)
         mat(k,1932) = -rxt(k,520)*y(k,271)
         mat(k,1066) = rxt(k,504)*y(k,265)
         mat(k,1010) = rxt(k,507)*y(k,265)
         mat(k,524) = rxt(k,521)*y(k,265)
         mat(k,1821) = rxt(k,504)*y(k,6) + rxt(k,507)*y(k,111) + rxt(k,521)*y(k,215)
         mat(k,788) = -(rxt(k,489)*y(k,243) + rxt(k,490)*y(k,125))
         mat(k,2328) = -rxt(k,489)*y(k,272)
         mat(k,1916) = -rxt(k,490)*y(k,272)
         mat(k,655) = rxt(k,491)*y(k,265)
         mat(k,235) = (.650_r8*rxt(k,492)+rxt(k,582))*y(k,265)
         mat(k,1798) = rxt(k,491)*y(k,218) + (.650_r8*rxt(k,492)+rxt(k,582))*y(k,219)
         mat(k,136) = -(rxt(k,583)*y(k,243) + rxt(k,584)*y(k,125))
         mat(k,2288) = -rxt(k,583)*y(k,273)
         mat(k,1889) = -rxt(k,584)*y(k,273)
         mat(k,230) = rxt(k,582)*y(k,265)
         mat(k,1715) = rxt(k,582)*y(k,219)
         mat(k,1278) = -(rxt(k,453)*y(k,236) + rxt(k,454)*y(k,237) + rxt(k,455) &
                      *y(k,243) + rxt(k,456)*y(k,125) + rxt(k,457)*y(k,127))
         mat(k,1464) = -rxt(k,453)*y(k,274)
         mat(k,2172) = -rxt(k,454)*y(k,274)
         mat(k,2355) = -rxt(k,455)*y(k,274)
         mat(k,1945) = -rxt(k,456)*y(k,274)
         mat(k,2416) = -rxt(k,457)*y(k,274)
         mat(k,283) = rxt(k,425)*y(k,265)
         mat(k,381) = rxt(k,426)*y(k,265)
         mat(k,199) = rxt(k,427)*y(k,265)
         mat(k,800) = .400_r8*rxt(k,450)*y(k,265)
         mat(k,248) = .500_r8*rxt(k,458)*y(k,265)
         mat(k,1835) = rxt(k,425)*y(k,94) + rxt(k,426)*y(k,96) + rxt(k,427)*y(k,97) &
                      + .400_r8*rxt(k,450)*y(k,103) + .500_r8*rxt(k,458)*y(k,220)
         mat(k,812) = -(rxt(k,495)*y(k,243) + rxt(k,496)*y(k,125))
         mat(k,2330) = -rxt(k,495)*y(k,275)
         mat(k,1917) = -rxt(k,496)*y(k,275)
         mat(k,259) = .560_r8*rxt(k,494)*y(k,265)
         mat(k,768) = rxt(k,497)*y(k,265)
         mat(k,1800) = .560_r8*rxt(k,494)*y(k,221) + rxt(k,497)*y(k,222)
         mat(k,142) = -(rxt(k,587)*y(k,243) + rxt(k,588)*y(k,125))
         mat(k,2289) = -rxt(k,587)*y(k,276)
         mat(k,1890) = -rxt(k,588)*y(k,276)
         mat(k,254) = rxt(k,586)*y(k,265)
         mat(k,1716) = rxt(k,586)*y(k,221)
         mat(k,558) = -(rxt(k,498)*y(k,243) + rxt(k,499)*y(k,125))
         mat(k,2314) = -rxt(k,498)*y(k,277)
         mat(k,1903) = -rxt(k,499)*y(k,277)
         mat(k,266) = .300_r8*rxt(k,500)*y(k,265)
         mat(k,470) = rxt(k,501)*y(k,265)
         mat(k,1773) = .300_r8*rxt(k,500)*y(k,223) + rxt(k,501)*y(k,224)
         mat(k,2526) = -(rxt(k,174)*y(k,261) + rxt(k,296)*y(k,73) + rxt(k,541) &
                      *y(k,159))
         mat(k,2243) = -rxt(k,174)*y(k,278)
         mat(k,978) = -rxt(k,296)*y(k,278)
         mat(k,312) = -rxt(k,541)*y(k,278)
         mat(k,332) = rxt(k,349)*y(k,265)
         mat(k,438) = rxt(k,374)*y(k,265)
         mat(k,338) = rxt(k,375)*y(k,265)
         mat(k,522) = rxt(k,301)*y(k,265)
         mat(k,2472) = rxt(k,320)*y(k,265)
         mat(k,689) = rxt(k,303)*y(k,265)
         mat(k,170) = rxt(k,304)*y(k,265)
         mat(k,1204) = rxt(k,351)*y(k,265)
         mat(k,418) = rxt(k,306)*y(k,265)
         mat(k,1130) = rxt(k,389)*y(k,265)
         mat(k,1333) = rxt(k,377)*y(k,265)
         mat(k,754) = rxt(k,357)*y(k,265)
         mat(k,671) = rxt(k,358)*y(k,265)
         mat(k,450) = rxt(k,326)*y(k,265)
         mat(k,1649) = rxt(k,327)*y(k,265)
         mat(k,2148) = rxt(k,194)*y(k,243)
         mat(k,1517) = rxt(k,199)*y(k,265)
         mat(k,649) = rxt(k,200)*y(k,265)
         mat(k,880) = rxt(k,287)*y(k,265)
         mat(k,360) = rxt(k,311)*y(k,265)
         mat(k,1552) = (rxt(k,597)+rxt(k,602))*y(k,91) + (rxt(k,590)+rxt(k,596) &
                       +rxt(k,601))*y(k,92) + rxt(k,258)*y(k,265)
         mat(k,959) = rxt(k,329)*y(k,265)
         mat(k,2266) = rxt(k,231)*y(k,265)
         mat(k,556) = rxt(k,208)*y(k,265)
         mat(k,839) = (rxt(k,597)+rxt(k,602))*y(k,85)
         mat(k,888) = (rxt(k,590)+rxt(k,596)+rxt(k,601))*y(k,85) + rxt(k,261)*y(k,265)
         mat(k,1324) = .500_r8*rxt(k,402)*y(k,265)
         mat(k,145) = rxt(k,542)*y(k,265)
         mat(k,608) = rxt(k,383)*y(k,265)
         mat(k,432) = rxt(k,387)*y(k,265)
         mat(k,2386) = rxt(k,194)*y(k,76) + rxt(k,201)*y(k,265)
         mat(k,1867) = rxt(k,349)*y(k,28) + rxt(k,374)*y(k,30) + rxt(k,375)*y(k,31) &
                      + rxt(k,301)*y(k,41) + rxt(k,320)*y(k,42) + rxt(k,303)*y(k,43) &
                      + rxt(k,304)*y(k,44) + rxt(k,351)*y(k,45) + rxt(k,306)*y(k,46) &
                      + rxt(k,389)*y(k,48) + rxt(k,377)*y(k,49) + rxt(k,357)*y(k,50) &
                      + rxt(k,358)*y(k,51) + rxt(k,326)*y(k,53) + rxt(k,327)*y(k,54) &
                      + rxt(k,199)*y(k,77) + rxt(k,200)*y(k,79) + rxt(k,287)*y(k,81) &
                      + rxt(k,311)*y(k,84) + rxt(k,258)*y(k,85) + rxt(k,329)*y(k,87) &
                      + rxt(k,231)*y(k,89) + rxt(k,208)*y(k,90) + rxt(k,261)*y(k,92) &
                      + .500_r8*rxt(k,402)*y(k,106) + rxt(k,542)*y(k,121) + rxt(k,383) &
                      *y(k,153) + rxt(k,387)*y(k,154) + rxt(k,201)*y(k,243) &
                      + 2.000_r8*rxt(k,204)*y(k,265)
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
         mat(k, 45) = lmat(k, 45)
         mat(k, 46) = lmat(k, 46)
         mat(k, 47) = lmat(k, 47)
         mat(k, 48) = lmat(k, 48)
         mat(k, 49) = lmat(k, 49)
         mat(k, 50) = lmat(k, 50)
         mat(k, 51) = lmat(k, 51)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = lmat(k, 53)
         mat(k, 54) = lmat(k, 54)
         mat(k, 55) = lmat(k, 55)
         mat(k, 56) = lmat(k, 56)
         mat(k, 57) = lmat(k, 57)
         mat(k, 58) = lmat(k, 58)
         mat(k, 59) = lmat(k, 59)
         mat(k, 60) = lmat(k, 60)
         mat(k, 61) = lmat(k, 61)
         mat(k, 62) = lmat(k, 62)
         mat(k, 63) = lmat(k, 63)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 69) = lmat(k, 69)
         mat(k, 70) = lmat(k, 70)
         mat(k, 71) = lmat(k, 71)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = lmat(k, 78)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 165) = mat(k, 165) + lmat(k, 165)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 175) = lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 198) = mat(k, 198) + lmat(k, 198)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 216) = lmat(k, 216)
         mat(k, 217) = lmat(k, 217)
         mat(k, 218) = lmat(k, 218)
         mat(k, 219) = lmat(k, 219)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 229) = lmat(k, 229)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 272) = lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = mat(k, 278) + lmat(k, 278)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 285) = lmat(k, 285)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 289) = mat(k, 289) + lmat(k, 289)
         mat(k, 290) = mat(k, 290) + lmat(k, 290)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 295) = mat(k, 295) + lmat(k, 295)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 309) = mat(k, 309) + lmat(k, 309)
         mat(k, 310) = lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 314) = lmat(k, 314)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 339) = lmat(k, 339)
         mat(k, 340) = lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 351) = lmat(k, 351)
         mat(k, 352) = lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 354) = mat(k, 354) + lmat(k, 354)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 361) = lmat(k, 361)
         mat(k, 362) = lmat(k, 362)
         mat(k, 363) = lmat(k, 363)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 366) = lmat(k, 366)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 368) = lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 374) = lmat(k, 374)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 378) = mat(k, 378) + lmat(k, 378)
         mat(k, 379) = lmat(k, 379)
         mat(k, 380) = mat(k, 380) + lmat(k, 380)
         mat(k, 383) = mat(k, 383) + lmat(k, 383)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 399) = lmat(k, 399)
         mat(k, 400) = lmat(k, 400)
         mat(k, 401) = lmat(k, 401)
         mat(k, 402) = lmat(k, 402)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 416) = lmat(k, 416)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 420) = lmat(k, 420)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 427) = mat(k, 427) + lmat(k, 427)
         mat(k, 429) = lmat(k, 429)
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 431) = lmat(k, 431)
         mat(k, 433) = mat(k, 433) + lmat(k, 433)
         mat(k, 435) = lmat(k, 435)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 437) = lmat(k, 437)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 440) = lmat(k, 440)
         mat(k, 442) = lmat(k, 442)
         mat(k, 444) = mat(k, 444) + lmat(k, 444)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 447) = lmat(k, 447)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 458) = lmat(k, 458)
         mat(k, 460) = lmat(k, 460)
         mat(k, 461) = mat(k, 461) + lmat(k, 461)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 467) = lmat(k, 467)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 471) = lmat(k, 471)
         mat(k, 472) = lmat(k, 472)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 474) = lmat(k, 474)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 483) = mat(k, 483) + lmat(k, 483)
         mat(k, 487) = mat(k, 487) + lmat(k, 487)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = mat(k, 489) + lmat(k, 489)
         mat(k, 490) = lmat(k, 490)
         mat(k, 492) = lmat(k, 492)
         mat(k, 493) = mat(k, 493) + lmat(k, 493)
         mat(k, 494) = lmat(k, 494)
         mat(k, 495) = lmat(k, 495)
         mat(k, 496) = lmat(k, 496)
         mat(k, 497) = lmat(k, 497)
         mat(k, 498) = lmat(k, 498)
         mat(k, 499) = lmat(k, 499)
         mat(k, 500) = lmat(k, 500)
         mat(k, 502) = mat(k, 502) + lmat(k, 502)
         mat(k, 508) = lmat(k, 508)
         mat(k, 509) = lmat(k, 509)
         mat(k, 510) = lmat(k, 510)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 515) = lmat(k, 515)
         mat(k, 516) = mat(k, 516) + lmat(k, 516)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 525) = lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 528) = lmat(k, 528)
         mat(k, 529) = lmat(k, 529)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 553) = mat(k, 553) + lmat(k, 553)
         mat(k, 554) = lmat(k, 554)
         mat(k, 555) = lmat(k, 555)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 565) = mat(k, 565) + lmat(k, 565)
         mat(k, 567) = lmat(k, 567)
         mat(k, 568) = lmat(k, 568)
         mat(k, 571) = mat(k, 571) + lmat(k, 571)
         mat(k, 577) = mat(k, 577) + lmat(k, 577)
         mat(k, 585) = mat(k, 585) + lmat(k, 585)
         mat(k, 586) = lmat(k, 586)
         mat(k, 587) = lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 590) = lmat(k, 590)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 593) = mat(k, 593) + lmat(k, 593)
         mat(k, 599) = lmat(k, 599)
         mat(k, 601) = mat(k, 601) + lmat(k, 601)
         mat(k, 603) = lmat(k, 603)
         mat(k, 605) = mat(k, 605) + lmat(k, 605)
         mat(k, 606) = lmat(k, 606)
         mat(k, 607) = lmat(k, 607)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 618) = lmat(k, 618)
         mat(k, 619) = lmat(k, 619)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 621) = lmat(k, 621)
         mat(k, 623) = lmat(k, 623)
         mat(k, 624) = mat(k, 624) + lmat(k, 624)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 627) = lmat(k, 627)
         mat(k, 629) = lmat(k, 629)
         mat(k, 634) = mat(k, 634) + lmat(k, 634)
         mat(k, 640) = lmat(k, 640)
         mat(k, 643) = mat(k, 643) + lmat(k, 643)
         mat(k, 645) = mat(k, 645) + lmat(k, 645)
         mat(k, 650) = lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = mat(k, 652) + lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 657) = lmat(k, 657)
         mat(k, 658) = lmat(k, 658)
         mat(k, 659) = mat(k, 659) + lmat(k, 659)
         mat(k, 660) = lmat(k, 660)
         mat(k, 661) = lmat(k, 661)
         mat(k, 662) = lmat(k, 662)
         mat(k, 663) = lmat(k, 663)
         mat(k, 664) = lmat(k, 664)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 667) = mat(k, 667) + lmat(k, 667)
         mat(k, 668) = mat(k, 668) + lmat(k, 668)
         mat(k, 669) = lmat(k, 669)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 681) = mat(k, 681) + lmat(k, 681)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 686) = lmat(k, 686)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 691) = mat(k, 691) + lmat(k, 691)
         mat(k, 692) = mat(k, 692) + lmat(k, 692)
         mat(k, 694) = lmat(k, 694)
         mat(k, 695) = lmat(k, 695)
         mat(k, 698) = mat(k, 698) + lmat(k, 698)
         mat(k, 704) = mat(k, 704) + lmat(k, 704)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 710) = lmat(k, 710)
         mat(k, 711) = mat(k, 711) + lmat(k, 711)
         mat(k, 715) = lmat(k, 715)
         mat(k, 716) = lmat(k, 716)
         mat(k, 718) = mat(k, 718) + lmat(k, 718)
         mat(k, 719) = lmat(k, 719)
         mat(k, 720) = lmat(k, 720)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 722) = lmat(k, 722)
         mat(k, 723) = lmat(k, 723)
         mat(k, 724) = lmat(k, 724)
         mat(k, 725) = lmat(k, 725)
         mat(k, 727) = lmat(k, 727)
         mat(k, 728) = mat(k, 728) + lmat(k, 728)
         mat(k, 729) = lmat(k, 729)
         mat(k, 730) = lmat(k, 730)
         mat(k, 731) = lmat(k, 731)
         mat(k, 732) = mat(k, 732) + lmat(k, 732)
         mat(k, 735) = mat(k, 735) + lmat(k, 735)
         mat(k, 736) = mat(k, 736) + lmat(k, 736)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 740) = lmat(k, 740)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 743) = mat(k, 743) + lmat(k, 743)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 755) = mat(k, 755) + lmat(k, 755)
         mat(k, 762) = lmat(k, 762)
         mat(k, 763) = lmat(k, 763)
         mat(k, 764) = lmat(k, 764)
         mat(k, 765) = lmat(k, 765)
         mat(k, 766) = mat(k, 766) + lmat(k, 766)
         mat(k, 771) = lmat(k, 771)
         mat(k, 773) = lmat(k, 773)
         mat(k, 775) = mat(k, 775) + lmat(k, 775)
         mat(k, 776) = lmat(k, 776)
         mat(k, 777) = mat(k, 777) + lmat(k, 777)
         mat(k, 788) = mat(k, 788) + lmat(k, 788)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 801) = lmat(k, 801)
         mat(k, 802) = lmat(k, 802)
         mat(k, 803) = mat(k, 803) + lmat(k, 803)
         mat(k, 804) = lmat(k, 804)
         mat(k, 805) = lmat(k, 805)
         mat(k, 812) = mat(k, 812) + lmat(k, 812)
         mat(k, 823) = mat(k, 823) + lmat(k, 823)
         mat(k, 832) = mat(k, 832) + lmat(k, 832)
         mat(k, 834) = lmat(k, 834)
         mat(k, 836) = mat(k, 836) + lmat(k, 836)
         mat(k, 842) = mat(k, 842) + lmat(k, 842)
         mat(k, 850) = mat(k, 850) + lmat(k, 850)
         mat(k, 851) = lmat(k, 851)
         mat(k, 853) = lmat(k, 853)
         mat(k, 859) = mat(k, 859) + lmat(k, 859)
         mat(k, 869) = mat(k, 869) + lmat(k, 869)
         mat(k, 873) = mat(k, 873) + lmat(k, 873)
         mat(k, 874) = mat(k, 874) + lmat(k, 874)
         mat(k, 878) = mat(k, 878) + lmat(k, 878)
         mat(k, 882) = mat(k, 882) + lmat(k, 882)
         mat(k, 884) = mat(k, 884) + lmat(k, 884)
         mat(k, 885) = mat(k, 885) + lmat(k, 885)
         mat(k, 890) = mat(k, 890) + lmat(k, 890)
         mat(k, 897) = mat(k, 897) + lmat(k, 897)
         mat(k, 898) = mat(k, 898) + lmat(k, 898)
         mat(k, 902) = mat(k, 902) + lmat(k, 902)
         mat(k, 910) = mat(k, 910) + lmat(k, 910)
         mat(k, 918) = mat(k, 918) + lmat(k, 918)
         mat(k, 920) = lmat(k, 920)
         mat(k, 921) = mat(k, 921) + lmat(k, 921)
         mat(k, 923) = lmat(k, 923)
         mat(k, 927) = mat(k, 927) + lmat(k, 927)
         mat(k, 946) = mat(k, 946) + lmat(k, 946)
         mat(k, 955) = mat(k, 955) + lmat(k, 955)
         mat(k, 960) = mat(k, 960) + lmat(k, 960)
         mat(k, 970) = mat(k, 970) + lmat(k, 970)
         mat(k, 980) = mat(k, 980) + lmat(k, 980)
         mat(k, 992) = mat(k, 992) + lmat(k, 992)
         mat(k, 993) = lmat(k, 993)
         mat(k, 995) = lmat(k, 995)
         mat(k,1007) = mat(k,1007) + lmat(k,1007)
         mat(k,1027) = lmat(k,1027)
         mat(k,1031) = lmat(k,1031)
         mat(k,1032) = mat(k,1032) + lmat(k,1032)
         mat(k,1043) = mat(k,1043) + lmat(k,1043)
         mat(k,1044) = mat(k,1044) + lmat(k,1044)
         mat(k,1046) = lmat(k,1046)
         mat(k,1048) = mat(k,1048) + lmat(k,1048)
         mat(k,1050) = mat(k,1050) + lmat(k,1050)
         mat(k,1051) = mat(k,1051) + lmat(k,1051)
         mat(k,1052) = mat(k,1052) + lmat(k,1052)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1087) = mat(k,1087) + lmat(k,1087)
         mat(k,1098) = lmat(k,1098)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1107) = mat(k,1107) + lmat(k,1107)
         mat(k,1117) = mat(k,1117) + lmat(k,1117)
         mat(k,1119) = lmat(k,1119)
         mat(k,1120) = lmat(k,1120)
         mat(k,1123) = lmat(k,1123)
         mat(k,1125) = mat(k,1125) + lmat(k,1125)
         mat(k,1127) = lmat(k,1127)
         mat(k,1129) = lmat(k,1129)
         mat(k,1131) = lmat(k,1131)
         mat(k,1135) = mat(k,1135) + lmat(k,1135)
         mat(k,1140) = lmat(k,1140)
         mat(k,1144) = lmat(k,1144)
         mat(k,1145) = mat(k,1145) + lmat(k,1145)
         mat(k,1151) = mat(k,1151) + lmat(k,1151)
         mat(k,1163) = mat(k,1163) + lmat(k,1163)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1165) = mat(k,1165) + lmat(k,1165)
         mat(k,1166) = mat(k,1166) + lmat(k,1166)
         mat(k,1167) = mat(k,1167) + lmat(k,1167)
         mat(k,1168) = mat(k,1168) + lmat(k,1168)
         mat(k,1171) = mat(k,1171) + lmat(k,1171)
         mat(k,1172) = mat(k,1172) + lmat(k,1172)
         mat(k,1178) = mat(k,1178) + lmat(k,1178)
         mat(k,1196) = mat(k,1196) + lmat(k,1196)
         mat(k,1197) = lmat(k,1197)
         mat(k,1200) = lmat(k,1200)
         mat(k,1202) = lmat(k,1202)
         mat(k,1205) = mat(k,1205) + lmat(k,1205)
         mat(k,1210) = lmat(k,1210)
         mat(k,1211) = mat(k,1211) + lmat(k,1211)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1215) = mat(k,1215) + lmat(k,1215)
         mat(k,1223) = mat(k,1223) + lmat(k,1223)
         mat(k,1236) = lmat(k,1236)
         mat(k,1237) = lmat(k,1237)
         mat(k,1238) = lmat(k,1238)
         mat(k,1239) = lmat(k,1239)
         mat(k,1240) = mat(k,1240) + lmat(k,1240)
         mat(k,1241) = lmat(k,1241)
         mat(k,1243) = lmat(k,1243)
         mat(k,1245) = lmat(k,1245)
         mat(k,1248) = lmat(k,1248)
         mat(k,1249) = mat(k,1249) + lmat(k,1249)
         mat(k,1251) = lmat(k,1251)
         mat(k,1253) = mat(k,1253) + lmat(k,1253)
         mat(k,1255) = lmat(k,1255)
         mat(k,1257) = mat(k,1257) + lmat(k,1257)
         mat(k,1258) = lmat(k,1258)
         mat(k,1262) = mat(k,1262) + lmat(k,1262)
         mat(k,1278) = mat(k,1278) + lmat(k,1278)
         mat(k,1298) = mat(k,1298) + lmat(k,1298)
         mat(k,1313) = mat(k,1313) + lmat(k,1313)
         mat(k,1314) = mat(k,1314) + lmat(k,1314)
         mat(k,1317) = mat(k,1317) + lmat(k,1317)
         mat(k,1318) = mat(k,1318) + lmat(k,1318)
         mat(k,1322) = mat(k,1322) + lmat(k,1322)
         mat(k,1323) = mat(k,1323) + lmat(k,1323)
         mat(k,1325) = mat(k,1325) + lmat(k,1325)
         mat(k,1326) = mat(k,1326) + lmat(k,1326)
         mat(k,1327) = mat(k,1327) + lmat(k,1327)
         mat(k,1331) = lmat(k,1331)
         mat(k,1344) = mat(k,1344) + lmat(k,1344)
         mat(k,1360) = lmat(k,1360)
         mat(k,1377) = mat(k,1377) + lmat(k,1377)
         mat(k,1389) = mat(k,1389) + lmat(k,1389)
         mat(k,1401) = mat(k,1401) + lmat(k,1401)
         mat(k,1415) = lmat(k,1415)
         mat(k,1417) = mat(k,1417) + lmat(k,1417)
         mat(k,1421) = mat(k,1421) + lmat(k,1421)
         mat(k,1423) = mat(k,1423) + lmat(k,1423)
         mat(k,1431) = lmat(k,1431)
         mat(k,1442) = mat(k,1442) + lmat(k,1442)
         mat(k,1473) = mat(k,1473) + lmat(k,1473)
         mat(k,1494) = mat(k,1494) + lmat(k,1494)
         mat(k,1495) = mat(k,1495) + lmat(k,1495)
         mat(k,1500) = lmat(k,1500)
         mat(k,1506) = mat(k,1506) + lmat(k,1506)
         mat(k,1519) = lmat(k,1519)
         mat(k,1521) = mat(k,1521) + lmat(k,1521)
         mat(k,1529) = mat(k,1529) + lmat(k,1529)
         mat(k,1540) = mat(k,1540) + lmat(k,1540)
         mat(k,1545) = mat(k,1545) + lmat(k,1545)
         mat(k,1547) = mat(k,1547) + lmat(k,1547)
         mat(k,1556) = mat(k,1556) + lmat(k,1556)
         mat(k,1576) = mat(k,1576) + lmat(k,1576)
         mat(k,1578) = mat(k,1578) + lmat(k,1578)
         mat(k,1579) = lmat(k,1579)
         mat(k,1587) = mat(k,1587) + lmat(k,1587)
         mat(k,1593) = mat(k,1593) + lmat(k,1593)
         mat(k,1596) = mat(k,1596) + lmat(k,1596)
         mat(k,1608) = mat(k,1608) + lmat(k,1608)
         mat(k,1610) = mat(k,1610) + lmat(k,1610)
         mat(k,1615) = mat(k,1615) + lmat(k,1615)
         mat(k,1627) = lmat(k,1627)
         mat(k,1628) = lmat(k,1628)
         mat(k,1629) = mat(k,1629) + lmat(k,1629)
         mat(k,1634) = mat(k,1634) + lmat(k,1634)
         mat(k,1636) = mat(k,1636) + lmat(k,1636)
         mat(k,1639) = lmat(k,1639)
         mat(k,1641) = mat(k,1641) + lmat(k,1641)
         mat(k,1642) = mat(k,1642) + lmat(k,1642)
         mat(k,1647) = mat(k,1647) + lmat(k,1647)
         mat(k,1649) = mat(k,1649) + lmat(k,1649)
         mat(k,1680) = mat(k,1680) + lmat(k,1680)
         mat(k,1681) = mat(k,1681) + lmat(k,1681)
         mat(k,1682) = mat(k,1682) + lmat(k,1682)
         mat(k,1684) = mat(k,1684) + lmat(k,1684)
         mat(k,1689) = mat(k,1689) + lmat(k,1689)
         mat(k,1854) = mat(k,1854) + lmat(k,1854)
         mat(k,1921) = mat(k,1921) + lmat(k,1921)
         mat(k,1923) = lmat(k,1923)
         mat(k,1929) = mat(k,1929) + lmat(k,1929)
         mat(k,1962) = mat(k,1962) + lmat(k,1962)
         mat(k,1964) = mat(k,1964) + lmat(k,1964)
         mat(k,2009) = mat(k,2009) + lmat(k,2009)
         mat(k,2031) = mat(k,2031) + lmat(k,2031)
         mat(k,2036) = lmat(k,2036)
         mat(k,2053) = mat(k,2053) + lmat(k,2053)
         mat(k,2069) = mat(k,2069) + lmat(k,2069)
         mat(k,2110) = mat(k,2110) + lmat(k,2110)
         mat(k,2117) = mat(k,2117) + lmat(k,2117)
         mat(k,2118) = mat(k,2118) + lmat(k,2118)
         mat(k,2121) = mat(k,2121) + lmat(k,2121)
         mat(k,2140) = mat(k,2140) + lmat(k,2140)
         mat(k,2193) = mat(k,2193) + lmat(k,2193)
         mat(k,2233) = mat(k,2233) + lmat(k,2233)
         mat(k,2237) = mat(k,2237) + lmat(k,2237)
         mat(k,2252) = lmat(k,2252)
         mat(k,2253) = mat(k,2253) + lmat(k,2253)
         mat(k,2261) = mat(k,2261) + lmat(k,2261)
         mat(k,2382) = mat(k,2382) + lmat(k,2382)
         mat(k,2386) = mat(k,2386) + lmat(k,2386)
         mat(k,2429) = mat(k,2429) + lmat(k,2429)
         mat(k,2432) = mat(k,2432) + lmat(k,2432)
         mat(k,2434) = mat(k,2434) + lmat(k,2434)
         mat(k,2436) = mat(k,2436) + lmat(k,2436)
         mat(k,2441) = mat(k,2441) + lmat(k,2441)
         mat(k,2443) = mat(k,2443) + lmat(k,2443)
         mat(k,2450) = mat(k,2450) + lmat(k,2450)
         mat(k,2452) = lmat(k,2452)
         mat(k,2464) = mat(k,2464) + lmat(k,2464)
         mat(k,2470) = mat(k,2470) + lmat(k,2470)
         mat(k,2488) = mat(k,2488) + lmat(k,2488)
         mat(k,2489) = mat(k,2489) + lmat(k,2489)
         mat(k,2498) = mat(k,2498) + lmat(k,2498)
         mat(k,2505) = lmat(k,2505)
         mat(k,2513) = mat(k,2513) + lmat(k,2513)
         mat(k,2516) = lmat(k,2516)
         mat(k,2518) = lmat(k,2518)
         mat(k,2520) = mat(k,2520) + lmat(k,2520)
         mat(k,2526) = mat(k,2526) + lmat(k,2526)
         mat(k, 260) = 0._r8
         mat(k, 261) = 0._r8
         mat(k, 300) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 372) = 0._r8
         mat(k, 478) = 0._r8
         mat(k, 480) = 0._r8
         mat(k, 505) = 0._r8
         mat(k, 533) = 0._r8
         mat(k, 536) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 654) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 712) = 0._r8
         mat(k, 713) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 767) = 0._r8
         mat(k, 769) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 774) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 814) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 820) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 844) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 868) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 896) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 912) = 0._r8
         mat(k, 913) = 0._r8
         mat(k, 914) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 990) = 0._r8
         mat(k, 997) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1025) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1033) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1041) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1064) = 0._r8
         mat(k,1065) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1081) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1224) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1328) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1341) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1351) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1363) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1367) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1374) = 0._r8
         mat(k,1375) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1380) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1384) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1402) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1413) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1429) = 0._r8
         mat(k,1432) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1441) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1476) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1486) = 0._r8
         mat(k,1490) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1493) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1502) = 0._r8
         mat(k,1503) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1513) = 0._r8
         mat(k,1515) = 0._r8
         mat(k,1516) = 0._r8
         mat(k,1532) = 0._r8
         mat(k,1534) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1550) = 0._r8
         mat(k,1555) = 0._r8
         mat(k,1559) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1562) = 0._r8
         mat(k,1564) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1568) = 0._r8
         mat(k,1582) = 0._r8
         mat(k,1584) = 0._r8
         mat(k,1586) = 0._r8
         mat(k,1588) = 0._r8
         mat(k,1589) = 0._r8
         mat(k,1590) = 0._r8
         mat(k,1592) = 0._r8
         mat(k,1598) = 0._r8
         mat(k,1599) = 0._r8
         mat(k,1600) = 0._r8
         mat(k,1607) = 0._r8
         mat(k,1616) = 0._r8
         mat(k,1617) = 0._r8
         mat(k,1618) = 0._r8
         mat(k,1619) = 0._r8
         mat(k,1621) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1624) = 0._r8
         mat(k,1631) = 0._r8
         mat(k,1632) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1635) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1640) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1646) = 0._r8
         mat(k,1648) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1662) = 0._r8
         mat(k,1667) = 0._r8
         mat(k,1668) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1670) = 0._r8
         mat(k,1673) = 0._r8
         mat(k,1675) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1679) = 0._r8
         mat(k,1683) = 0._r8
         mat(k,1686) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1692) = 0._r8
         mat(k,1694) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1797) = 0._r8
         mat(k,1801) = 0._r8
         mat(k,1811) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1836) = 0._r8
         mat(k,1861) = 0._r8
         mat(k,1924) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1966) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1988) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1992) = 0._r8
         mat(k,1993) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1996) = 0._r8
         mat(k,1997) = 0._r8
         mat(k,1998) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2004) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2015) = 0._r8
         mat(k,2020) = 0._r8
         mat(k,2029) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2037) = 0._r8
         mat(k,2041) = 0._r8
         mat(k,2048) = 0._r8
         mat(k,2056) = 0._r8
         mat(k,2057) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2084) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2089) = 0._r8
         mat(k,2094) = 0._r8
         mat(k,2095) = 0._r8
         mat(k,2096) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2101) = 0._r8
         mat(k,2102) = 0._r8
         mat(k,2103) = 0._r8
         mat(k,2105) = 0._r8
         mat(k,2122) = 0._r8
         mat(k,2127) = 0._r8
         mat(k,2129) = 0._r8
         mat(k,2130) = 0._r8
         mat(k,2132) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2134) = 0._r8
         mat(k,2136) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2141) = 0._r8
         mat(k,2142) = 0._r8
         mat(k,2143) = 0._r8
         mat(k,2145) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2147) = 0._r8
         mat(k,2156) = 0._r8
         mat(k,2184) = 0._r8
         mat(k,2185) = 0._r8
         mat(k,2187) = 0._r8
         mat(k,2190) = 0._r8
         mat(k,2191) = 0._r8
         mat(k,2192) = 0._r8
         mat(k,2194) = 0._r8
         mat(k,2195) = 0._r8
         mat(k,2197) = 0._r8
         mat(k,2200) = 0._r8
         mat(k,2229) = 0._r8
         mat(k,2238) = 0._r8
         mat(k,2240) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2247) = 0._r8
         mat(k,2248) = 0._r8
         mat(k,2249) = 0._r8
         mat(k,2250) = 0._r8
         mat(k,2251) = 0._r8
         mat(k,2254) = 0._r8
         mat(k,2255) = 0._r8
         mat(k,2256) = 0._r8
         mat(k,2257) = 0._r8
         mat(k,2258) = 0._r8
         mat(k,2259) = 0._r8
         mat(k,2260) = 0._r8
         mat(k,2262) = 0._r8
         mat(k,2264) = 0._r8
         mat(k,2265) = 0._r8
         mat(k,2308) = 0._r8
         mat(k,2309) = 0._r8
         mat(k,2312) = 0._r8
         mat(k,2319) = 0._r8
         mat(k,2337) = 0._r8
         mat(k,2345) = 0._r8
         mat(k,2346) = 0._r8
         mat(k,2348) = 0._r8
         mat(k,2351) = 0._r8
         mat(k,2353) = 0._r8
         mat(k,2357) = 0._r8
         mat(k,2362) = 0._r8
         mat(k,2371) = 0._r8
         mat(k,2380) = 0._r8
         mat(k,2381) = 0._r8
         mat(k,2393) = 0._r8
         mat(k,2396) = 0._r8
         mat(k,2400) = 0._r8
         mat(k,2402) = 0._r8
         mat(k,2407) = 0._r8
         mat(k,2413) = 0._r8
         mat(k,2415) = 0._r8
         mat(k,2426) = 0._r8
         mat(k,2427) = 0._r8
         mat(k,2428) = 0._r8
         mat(k,2430) = 0._r8
         mat(k,2431) = 0._r8
         mat(k,2435) = 0._r8
         mat(k,2437) = 0._r8
         mat(k,2438) = 0._r8
         mat(k,2439) = 0._r8
         mat(k,2440) = 0._r8
         mat(k,2445) = 0._r8
         mat(k,2446) = 0._r8
         mat(k,2449) = 0._r8
         mat(k,2451) = 0._r8
         mat(k,2455) = 0._r8
         mat(k,2456) = 0._r8
         mat(k,2457) = 0._r8
         mat(k,2458) = 0._r8
         mat(k,2460) = 0._r8
         mat(k,2463) = 0._r8
         mat(k,2465) = 0._r8
         mat(k,2466) = 0._r8
         mat(k,2471) = 0._r8
         mat(k,2490) = 0._r8
         mat(k,2491) = 0._r8
         mat(k,2493) = 0._r8
         mat(k,2494) = 0._r8
         mat(k,2496) = 0._r8
         mat(k,2499) = 0._r8
         mat(k,2504) = 0._r8
         mat(k,2506) = 0._r8
         mat(k,2507) = 0._r8
         mat(k,2508) = 0._r8
         mat(k,2509) = 0._r8
         mat(k,2510) = 0._r8
         mat(k,2511) = 0._r8
         mat(k,2512) = 0._r8
         mat(k,2514) = 0._r8
         mat(k,2515) = 0._r8
         mat(k,2517) = 0._r8
         mat(k,2519) = 0._r8
         mat(k,2521) = 0._r8
         mat(k,2522) = 0._r8
         mat(k,2523) = 0._r8
         mat(k,2524) = 0._r8
         mat(k,2525) = 0._r8
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
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 46) = mat(k, 46) - dti(k)
         mat(k, 47) = mat(k, 47) - dti(k)
         mat(k, 48) = mat(k, 48) - dti(k)
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 50) = mat(k, 50) - dti(k)
         mat(k, 51) = mat(k, 51) - dti(k)
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 54) = mat(k, 54) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 56) = mat(k, 56) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 58) = mat(k, 58) - dti(k)
         mat(k, 59) = mat(k, 59) - dti(k)
         mat(k, 60) = mat(k, 60) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 62) = mat(k, 62) - dti(k)
         mat(k, 63) = mat(k, 63) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 66) = mat(k, 66) - dti(k)
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 68) = mat(k, 68) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 70) = mat(k, 70) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 72) = mat(k, 72) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 78) = mat(k, 78) - dti(k)
         mat(k, 84) = mat(k, 84) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 114) = mat(k, 114) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 124) = mat(k, 124) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 143) = mat(k, 143) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 149) = mat(k, 149) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 159) = mat(k, 159) - dti(k)
         mat(k, 163) = mat(k, 163) - dti(k)
         mat(k, 167) = mat(k, 167) - dti(k)
         mat(k, 171) = mat(k, 171) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 188) = mat(k, 188) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 198) = mat(k, 198) - dti(k)
         mat(k, 201) = mat(k, 201) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 210) = mat(k, 210) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 220) = mat(k, 220) - dti(k)
         mat(k, 225) = mat(k, 225) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 231) = mat(k, 231) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 243) = mat(k, 243) - dti(k)
         mat(k, 247) = mat(k, 247) - dti(k)
         mat(k, 250) = mat(k, 250) - dti(k)
         mat(k, 256) = mat(k, 256) - dti(k)
         mat(k, 264) = mat(k, 264) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 272) = mat(k, 272) - dti(k)
         mat(k, 275) = mat(k, 275) - dti(k)
         mat(k, 278) = mat(k, 278) - dti(k)
         mat(k, 281) = mat(k, 281) - dti(k)
         mat(k, 286) = mat(k, 286) - dti(k)
         mat(k, 289) = mat(k, 289) - dti(k)
         mat(k, 294) = mat(k, 294) - dti(k)
         mat(k, 299) = mat(k, 299) - dti(k)
         mat(k, 304) = mat(k, 304) - dti(k)
         mat(k, 309) = mat(k, 309) - dti(k)
         mat(k, 313) = mat(k, 313) - dti(k)
         mat(k, 317) = mat(k, 317) - dti(k)
         mat(k, 321) = mat(k, 321) - dti(k)
         mat(k, 327) = mat(k, 327) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 345) = mat(k, 345) - dti(k)
         mat(k, 351) = mat(k, 351) - dti(k)
         mat(k, 355) = mat(k, 355) - dti(k)
         mat(k, 361) = mat(k, 361) - dti(k)
         mat(k, 364) = mat(k, 364) - dti(k)
         mat(k, 370) = mat(k, 370) - dti(k)
         mat(k, 375) = mat(k, 375) - dti(k)
         mat(k, 380) = mat(k, 380) - dti(k)
         mat(k, 383) = mat(k, 383) - dti(k)
         mat(k, 388) = mat(k, 388) - dti(k)
         mat(k, 393) = mat(k, 393) - dti(k)
         mat(k, 398) = mat(k, 398) - dti(k)
         mat(k, 403) = mat(k, 403) - dti(k)
         mat(k, 411) = mat(k, 411) - dti(k)
         mat(k, 419) = mat(k, 419) - dti(k)
         mat(k, 427) = mat(k, 427) - dti(k)
         mat(k, 433) = mat(k, 433) - dti(k)
         mat(k, 439) = mat(k, 439) - dti(k)
         mat(k, 445) = mat(k, 445) - dti(k)
         mat(k, 451) = mat(k, 451) - dti(k)
         mat(k, 457) = mat(k, 457) - dti(k)
         mat(k, 463) = mat(k, 463) - dti(k)
         mat(k, 469) = mat(k, 469) - dti(k)
         mat(k, 477) = mat(k, 477) - dti(k)
         mat(k, 483) = mat(k, 483) - dti(k)
         mat(k, 489) = mat(k, 489) - dti(k)
         mat(k, 495) = mat(k, 495) - dti(k)
         mat(k, 502) = mat(k, 502) - dti(k)
         mat(k, 508) = mat(k, 508) - dti(k)
         mat(k, 511) = mat(k, 511) - dti(k)
         mat(k, 516) = mat(k, 516) - dti(k)
         mat(k, 523) = mat(k, 523) - dti(k)
         mat(k, 532) = mat(k, 532) - dti(k)
         mat(k, 539) = mat(k, 539) - dti(k)
         mat(k, 543) = mat(k, 543) - dti(k)
         mat(k, 550) = mat(k, 550) - dti(k)
         mat(k, 558) = mat(k, 558) - dti(k)
         mat(k, 565) = mat(k, 565) - dti(k)
         mat(k, 571) = mat(k, 571) - dti(k)
         mat(k, 577) = mat(k, 577) - dti(k)
         mat(k, 585) = mat(k, 585) - dti(k)
         mat(k, 593) = mat(k, 593) - dti(k)
         mat(k, 601) = mat(k, 601) - dti(k)
         mat(k, 609) = mat(k, 609) - dti(k)
         mat(k, 617) = mat(k, 617) - dti(k)
         mat(k, 625) = mat(k, 625) - dti(k)
         mat(k, 634) = mat(k, 634) - dti(k)
         mat(k, 643) = mat(k, 643) - dti(k)
         mat(k, 652) = mat(k, 652) - dti(k)
         mat(k, 661) = mat(k, 661) - dti(k)
         mat(k, 665) = mat(k, 665) - dti(k)
         mat(k, 672) = mat(k, 672) - dti(k)
         mat(k, 681) = mat(k, 681) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 698) = mat(k, 698) - dti(k)
         mat(k, 704) = mat(k, 704) - dti(k)
         mat(k, 711) = mat(k, 711) - dti(k)
         mat(k, 721) = mat(k, 721) - dti(k)
         mat(k, 732) = mat(k, 732) - dti(k)
         mat(k, 743) = mat(k, 743) - dti(k)
         mat(k, 750) = mat(k, 750) - dti(k)
         mat(k, 755) = mat(k, 755) - dti(k)
         mat(k, 766) = mat(k, 766) - dti(k)
         mat(k, 777) = mat(k, 777) - dti(k)
         mat(k, 788) = mat(k, 788) - dti(k)
         mat(k, 799) = mat(k, 799) - dti(k)
         mat(k, 812) = mat(k, 812) - dti(k)
         mat(k, 823) = mat(k, 823) - dti(k)
         mat(k, 832) = mat(k, 832) - dti(k)
         mat(k, 842) = mat(k, 842) - dti(k)
         mat(k, 850) = mat(k, 850) - dti(k)
         mat(k, 859) = mat(k, 859) - dti(k)
         mat(k, 869) = mat(k, 869) - dti(k)
         mat(k, 873) = mat(k, 873) - dti(k)
         mat(k, 882) = mat(k, 882) - dti(k)
         mat(k, 890) = mat(k, 890) - dti(k)
         mat(k, 898) = mat(k, 898) - dti(k)
         mat(k, 910) = mat(k, 910) - dti(k)
         mat(k, 918) = mat(k, 918) - dti(k)
         mat(k, 927) = mat(k, 927) - dti(k)
         mat(k, 946) = mat(k, 946) - dti(k)
         mat(k, 955) = mat(k, 955) - dti(k)
         mat(k, 960) = mat(k, 960) - dti(k)
         mat(k, 970) = mat(k, 970) - dti(k)
         mat(k, 980) = mat(k, 980) - dti(k)
         mat(k, 992) = mat(k, 992) - dti(k)
         mat(k,1007) = mat(k,1007) - dti(k)
         mat(k,1032) = mat(k,1032) - dti(k)
         mat(k,1044) = mat(k,1044) - dti(k)
         mat(k,1063) = mat(k,1063) - dti(k)
         mat(k,1087) = mat(k,1087) - dti(k)
         mat(k,1099) = mat(k,1099) - dti(k)
         mat(k,1107) = mat(k,1107) - dti(k)
         mat(k,1117) = mat(k,1117) - dti(k)
         mat(k,1125) = mat(k,1125) - dti(k)
         mat(k,1135) = mat(k,1135) - dti(k)
         mat(k,1151) = mat(k,1151) - dti(k)
         mat(k,1164) = mat(k,1164) - dti(k)
         mat(k,1178) = mat(k,1178) - dti(k)
         mat(k,1196) = mat(k,1196) - dti(k)
         mat(k,1205) = mat(k,1205) - dti(k)
         mat(k,1211) = mat(k,1211) - dti(k)
         mat(k,1223) = mat(k,1223) - dti(k)
         mat(k,1240) = mat(k,1240) - dti(k)
         mat(k,1253) = mat(k,1253) - dti(k)
         mat(k,1262) = mat(k,1262) - dti(k)
         mat(k,1278) = mat(k,1278) - dti(k)
         mat(k,1298) = mat(k,1298) - dti(k)
         mat(k,1314) = mat(k,1314) - dti(k)
         mat(k,1326) = mat(k,1326) - dti(k)
         mat(k,1344) = mat(k,1344) - dti(k)
         mat(k,1377) = mat(k,1377) - dti(k)
         mat(k,1401) = mat(k,1401) - dti(k)
         mat(k,1421) = mat(k,1421) - dti(k)
         mat(k,1442) = mat(k,1442) - dti(k)
         mat(k,1473) = mat(k,1473) - dti(k)
         mat(k,1495) = mat(k,1495) - dti(k)
         mat(k,1506) = mat(k,1506) - dti(k)
         mat(k,1521) = mat(k,1521) - dti(k)
         mat(k,1540) = mat(k,1540) - dti(k)
         mat(k,1556) = mat(k,1556) - dti(k)
         mat(k,1587) = mat(k,1587) - dti(k)
         mat(k,1610) = mat(k,1610) - dti(k)
         mat(k,1634) = mat(k,1634) - dti(k)
         mat(k,1680) = mat(k,1680) - dti(k)
         mat(k,1854) = mat(k,1854) - dti(k)
         mat(k,1962) = mat(k,1962) - dti(k)
         mat(k,2009) = mat(k,2009) - dti(k)
         mat(k,2053) = mat(k,2053) - dti(k)
         mat(k,2118) = mat(k,2118) - dti(k)
         mat(k,2140) = mat(k,2140) - dti(k)
         mat(k,2193) = mat(k,2193) - dti(k)
         mat(k,2237) = mat(k,2237) - dti(k)
         mat(k,2261) = mat(k,2261) - dti(k)
         mat(k,2382) = mat(k,2382) - dti(k)
         mat(k,2443) = mat(k,2443) - dti(k)
         mat(k,2470) = mat(k,2470) - dti(k)
         mat(k,2498) = mat(k,2498) - dti(k)
         mat(k,2526) = mat(k,2526) - dti(k)
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
