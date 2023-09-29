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
         mat(k,699) = -(rxt(k,416)*y(k,265))
         mat(k,1940) = -rxt(k,416)*y(k,1)
         mat(k,1740) = rxt(k,419)*y(k,234)
         mat(k,1082) = rxt(k,419)*y(k,131)
         mat(k,710) = -(rxt(k,420)*y(k,265))
         mat(k,1941) = -rxt(k,420)*y(k,2)
         mat(k,1083) = rxt(k,417)*y(k,247)
         mat(k,2158) = rxt(k,417)*y(k,234)
         mat(k,1062) = -(rxt(k,499)*y(k,133) + rxt(k,500)*y(k,143) + rxt(k,501) &
                      *y(k,265))
         mat(k,2030) = -rxt(k,499)*y(k,6)
         mat(k,2384) = -rxt(k,500)*y(k,6)
         mat(k,1968) = -rxt(k,501)*y(k,6)
         mat(k,96) = -(rxt(k,554)*y(k,247) + rxt(k,555)*y(k,131))
         mat(k,2116) = -rxt(k,554)*y(k,7)
         mat(k,1709) = -rxt(k,555)*y(k,7)
         mat(k,1058) = rxt(k,557)*y(k,265)
         mat(k,1854) = rxt(k,557)*y(k,6)
         mat(k,208) = -(rxt(k,458)*y(k,265))
         mat(k,1871) = -rxt(k,458)*y(k,8)
         mat(k,107) = -(rxt(k,559)*y(k,247) + rxt(k,560)*y(k,131))
         mat(k,2122) = -rxt(k,559)*y(k,9)
         mat(k,1715) = -rxt(k,560)*y(k,9)
         mat(k,207) = rxt(k,558)*y(k,265)
         mat(k,1860) = rxt(k,558)*y(k,8)
         mat(k,449) = -(rxt(k,461)*y(k,265))
         mat(k,1910) = -rxt(k,461)*y(k,10)
         mat(k,539) = rxt(k,459)*y(k,247)
         mat(k,2138) = rxt(k,459)*y(k,235)
         mat(k,209) = .120_r8*rxt(k,458)*y(k,265)
         mat(k,1872) = .120_r8*rxt(k,458)*y(k,8)
         mat(k,1060) = .100_r8*rxt(k,500)*y(k,143)
         mat(k,1004) = .100_r8*rxt(k,503)*y(k,143)
         mat(k,2373) = .100_r8*rxt(k,500)*y(k,6) + .100_r8*rxt(k,503)*y(k,116)
         mat(k,1727) = .500_r8*rxt(k,460)*y(k,235) + .200_r8*rxt(k,487)*y(k,272) &
                      + .060_r8*rxt(k,493)*y(k,274)
         mat(k,540) = .500_r8*rxt(k,460)*y(k,131)
         mat(k,791) = .200_r8*rxt(k,487)*y(k,131)
         mat(k,807) = .060_r8*rxt(k,493)*y(k,131)
         mat(k,1721) = .200_r8*rxt(k,487)*y(k,272) + .200_r8*rxt(k,493)*y(k,274)
         mat(k,790) = .200_r8*rxt(k,487)*y(k,131)
         mat(k,805) = .200_r8*rxt(k,493)*y(k,131)
         mat(k,1737) = .200_r8*rxt(k,487)*y(k,272) + .150_r8*rxt(k,493)*y(k,274)
         mat(k,793) = .200_r8*rxt(k,487)*y(k,131)
         mat(k,808) = .150_r8*rxt(k,493)*y(k,131)
         mat(k,1722) = .210_r8*rxt(k,493)*y(k,274)
         mat(k,806) = .210_r8*rxt(k,493)*y(k,131)
         mat(k,270) = -(rxt(k,421)*y(k,265))
         mat(k,1882) = -rxt(k,421)*y(k,17)
         mat(k,1059) = .050_r8*rxt(k,500)*y(k,143)
         mat(k,1003) = .050_r8*rxt(k,503)*y(k,143)
         mat(k,2372) = .050_r8*rxt(k,500)*y(k,6) + .050_r8*rxt(k,503)*y(k,116)
         mat(k,409) = -(rxt(k,387)*y(k,133) + rxt(k,388)*y(k,265))
         mat(k,2023) = -rxt(k,387)*y(k,18)
         mat(k,1904) = -rxt(k,388)*y(k,18)
         mat(k,1555) = -(rxt(k,270)*y(k,44) + rxt(k,271)*y(k,247) + rxt(k,272) &
                      *y(k,143))
         mat(k,2479) = -rxt(k,270)*y(k,19)
         mat(k,2204) = -rxt(k,271)*y(k,19)
         mat(k,2410) = -rxt(k,272)*y(k,19)
         mat(k,1607) = 4.000_r8*rxt(k,273)*y(k,21) + (rxt(k,274)+rxt(k,275))*y(k,61) &
                      + rxt(k,278)*y(k,131) + rxt(k,281)*y(k,141) + rxt(k,529) &
                      *y(k,161) + rxt(k,282)*y(k,265)
         mat(k,187) = rxt(k,260)*y(k,261)
         mat(k,193) = rxt(k,286)*y(k,261)
         mat(k,516) = 2.000_r8*rxt(k,297)*y(k,58) + 2.000_r8*rxt(k,309)*y(k,261) &
                      + 2.000_r8*rxt(k,298)*y(k,265)
         mat(k,630) = rxt(k,299)*y(k,58) + rxt(k,310)*y(k,261) + rxt(k,300)*y(k,265)
         mat(k,463) = 3.000_r8*rxt(k,304)*y(k,58) + 3.000_r8*rxt(k,287)*y(k,261) &
                      + 3.000_r8*rxt(k,305)*y(k,265)
         mat(k,2250) = 2.000_r8*rxt(k,297)*y(k,43) + rxt(k,299)*y(k,45) &
                      + 3.000_r8*rxt(k,304)*y(k,57)
         mat(k,1633) = (rxt(k,274)+rxt(k,275))*y(k,21)
         mat(k,171) = 2.000_r8*rxt(k,288)*y(k,261)
         mat(k,881) = rxt(k,283)*y(k,141) + rxt(k,289)*y(k,261) + rxt(k,284)*y(k,265)
         mat(k,1783) = rxt(k,278)*y(k,21)
         mat(k,2453) = rxt(k,281)*y(k,21) + rxt(k,283)*y(k,83)
         mat(k,1521) = rxt(k,529)*y(k,21)
         mat(k,1826) = rxt(k,260)*y(k,36) + rxt(k,286)*y(k,37) + 2.000_r8*rxt(k,309) &
                      *y(k,43) + rxt(k,310)*y(k,45) + 3.000_r8*rxt(k,287)*y(k,57) &
                      + 2.000_r8*rxt(k,288)*y(k,80) + rxt(k,289)*y(k,83)
         mat(k,1999) = rxt(k,282)*y(k,21) + 2.000_r8*rxt(k,298)*y(k,43) + rxt(k,300) &
                      *y(k,45) + 3.000_r8*rxt(k,305)*y(k,57) + rxt(k,284)*y(k,83)
         mat(k,1600) = rxt(k,276)*y(k,61)
         mat(k,1626) = rxt(k,276)*y(k,21)
         mat(k,1535) = (rxt(k,594)+rxt(k,599))*y(k,93)
         mat(k,840) = (rxt(k,594)+rxt(k,599))*y(k,87)
         mat(k,1609) = -(4._r8*rxt(k,273)*y(k,21) + (rxt(k,274) + rxt(k,275) + rxt(k,276) &
                      ) * y(k,61) + rxt(k,277)*y(k,247) + rxt(k,278)*y(k,131) &
                      + rxt(k,279)*y(k,132) + rxt(k,281)*y(k,141) + rxt(k,282) &
                      *y(k,265) + rxt(k,529)*y(k,161))
         mat(k,1635) = -(rxt(k,274) + rxt(k,275) + rxt(k,276)) * y(k,21)
         mat(k,2206) = -rxt(k,277)*y(k,21)
         mat(k,1785) = -rxt(k,278)*y(k,21)
         mat(k,1679) = -rxt(k,279)*y(k,21)
         mat(k,2455) = -rxt(k,281)*y(k,21)
         mat(k,2001) = -rxt(k,282)*y(k,21)
         mat(k,1523) = -rxt(k,529)*y(k,21)
         mat(k,1557) = rxt(k,272)*y(k,143)
         mat(k,607) = rxt(k,280)*y(k,141)
         mat(k,882) = rxt(k,290)*y(k,261)
         mat(k,844) = rxt(k,285)*y(k,141)
         mat(k,2455) = mat(k,2455) + rxt(k,280)*y(k,22) + rxt(k,285)*y(k,93)
         mat(k,2412) = rxt(k,272)*y(k,19)
         mat(k,1828) = rxt(k,290)*y(k,83)
         mat(k,604) = -(rxt(k,280)*y(k,141))
         mat(k,2434) = -rxt(k,280)*y(k,22)
         mat(k,1602) = rxt(k,279)*y(k,132)
         mat(k,1658) = rxt(k,279)*y(k,21)
         mat(k,279) = -(rxt(k,462)*y(k,265))
         mat(k,1884) = -rxt(k,462)*y(k,24)
         mat(k,1720) = rxt(k,465)*y(k,236)
         mat(k,479) = rxt(k,465)*y(k,131)
         mat(k,373) = -(rxt(k,464)*y(k,265))
         mat(k,1898) = -rxt(k,464)*y(k,25)
         mat(k,480) = rxt(k,463)*y(k,247)
         mat(k,2130) = rxt(k,463)*y(k,236)
         mat(k,325) = -(rxt(k,335)*y(k,58) + rxt(k,336)*y(k,265))
         mat(k,2224) = -rxt(k,335)*y(k,26)
         mat(k,1892) = -rxt(k,336)*y(k,26)
         mat(k,620) = -(rxt(k,337)*y(k,58) + rxt(k,338)*y(k,143) + rxt(k,363)*y(k,265))
         mat(k,2230) = -rxt(k,337)*y(k,27)
         mat(k,2375) = -rxt(k,338)*y(k,27)
         mat(k,1930) = -rxt(k,363)*y(k,27)
         mat(k,311) = -(rxt(k,343)*y(k,265))
         mat(k,1890) = -rxt(k,343)*y(k,28)
         mat(k,947) = .800_r8*rxt(k,339)*y(k,237) + .200_r8*rxt(k,340)*y(k,241)
         mat(k,2313) = .200_r8*rxt(k,340)*y(k,237)
         mat(k,378) = -(rxt(k,344)*y(k,265))
         mat(k,1899) = -rxt(k,344)*y(k,29)
         mat(k,948) = rxt(k,341)*y(k,247)
         mat(k,2131) = rxt(k,341)*y(k,237)
         mat(k,334) = -(rxt(k,345)*y(k,58) + rxt(k,346)*y(k,265))
         mat(k,2225) = -rxt(k,345)*y(k,30)
         mat(k,1893) = -rxt(k,346)*y(k,30)
         mat(k,1177) = -(rxt(k,366)*y(k,133) + rxt(k,367)*y(k,143) + rxt(k,385) &
                      *y(k,265))
         mat(k,2039) = -rxt(k,366)*y(k,31)
         mat(k,2392) = -rxt(k,367)*y(k,31)
         mat(k,1977) = -rxt(k,385)*y(k,31)
         mat(k,922) = .130_r8*rxt(k,445)*y(k,143)
         mat(k,2392) = mat(k,2392) + .130_r8*rxt(k,445)*y(k,100)
         mat(k,437) = -(rxt(k,371)*y(k,265))
         mat(k,1908) = -rxt(k,371)*y(k,32)
         mat(k,978) = rxt(k,369)*y(k,247)
         mat(k,2136) = rxt(k,369)*y(k,238)
         mat(k,340) = -(rxt(k,372)*y(k,265) + rxt(k,375)*y(k,58))
         mat(k,1894) = -rxt(k,372)*y(k,33)
         mat(k,2226) = -rxt(k,375)*y(k,33)
         mat(k,315) = -(rxt(k,468)*y(k,265))
         mat(k,1891) = -rxt(k,468)*y(k,34)
         mat(k,690) = rxt(k,466)*y(k,247)
         mat(k,2128) = rxt(k,466)*y(k,239)
         mat(k,144) = -(rxt(k,259)*y(k,261))
         mat(k,1802) = -rxt(k,259)*y(k,35)
         mat(k,185) = -(rxt(k,260)*y(k,261))
         mat(k,1807) = -rxt(k,260)*y(k,36)
         mat(k,190) = -(rxt(k,286)*y(k,261))
         mat(k,1808) = -rxt(k,286)*y(k,37)
         mat(k,153) = -(rxt(k,261)*y(k,261))
         mat(k,1803) = -rxt(k,261)*y(k,38)
         mat(k,195) = -(rxt(k,262)*y(k,261))
         mat(k,1809) = -rxt(k,262)*y(k,39)
         mat(k,157) = -(rxt(k,263)*y(k,261))
         mat(k,1804) = -rxt(k,263)*y(k,40)
         mat(k,200) = -(rxt(k,264)*y(k,261))
         mat(k,1810) = -rxt(k,264)*y(k,41)
         mat(k,161) = -(rxt(k,265)*y(k,261))
         mat(k,1805) = -rxt(k,265)*y(k,42)
         mat(k,514) = -(rxt(k,297)*y(k,58) + rxt(k,298)*y(k,265) + rxt(k,309)*y(k,261))
         mat(k,2229) = -rxt(k,297)*y(k,43)
         mat(k,1918) = -rxt(k,298)*y(k,43)
         mat(k,1820) = -rxt(k,309)*y(k,43)
         mat(k,2496) = -(rxt(k,234)*y(k,58) + rxt(k,270)*y(k,19) + rxt(k,314)*y(k,247) &
                      + rxt(k,315)*y(k,133) + rxt(k,316)*y(k,141) + rxt(k,317) &
                      *y(k,265))
         mat(k,2267) = -rxt(k,234)*y(k,44)
         mat(k,1566) = -rxt(k,270)*y(k,44)
         mat(k,2221) = -rxt(k,314)*y(k,44)
         mat(k,2076) = -rxt(k,315)*y(k,44)
         mat(k,2470) = -rxt(k,316)*y(k,44)
         mat(k,2016) = -rxt(k,317)*y(k,44)
         mat(k,708) = .400_r8*rxt(k,416)*y(k,265)
         mat(k,1079) = .340_r8*rxt(k,500)*y(k,143)
         mat(k,416) = .500_r8*rxt(k,387)*y(k,133)
         mat(k,627) = rxt(k,338)*y(k,143)
         mat(k,1193) = .500_r8*rxt(k,367)*y(k,143)
         mat(k,642) = .500_r8*rxt(k,355)*y(k,265)
         mat(k,860) = rxt(k,322)*y(k,265)
         mat(k,459) = .300_r8*rxt(k,323)*y(k,265)
         mat(k,2101) = (rxt(k,331)+rxt(k,332))*y(k,261)
         mat(k,1649) = rxt(k,241)*y(k,241)
         mat(k,1214) = .800_r8*rxt(k,360)*y(k,265)
         mat(k,935) = .910_r8*rxt(k,445)*y(k,143)
         mat(k,659) = .300_r8*rxt(k,436)*y(k,265)
         mat(k,1310) = .800_r8*rxt(k,440)*y(k,241)
         mat(k,1322) = .120_r8*rxt(k,398)*y(k,143)
         mat(k,679) = .500_r8*rxt(k,411)*y(k,265)
         mat(k,1023) = .340_r8*rxt(k,503)*y(k,143)
         mat(k,1434) = .600_r8*rxt(k,412)*y(k,143)
         mat(k,1800) = .100_r8*rxt(k,418)*y(k,234) + rxt(k,321)*y(k,241) &
                      + .500_r8*rxt(k,389)*y(k,244) + .500_r8*rxt(k,357)*y(k,246) &
                      + .920_r8*rxt(k,428)*y(k,249) + .250_r8*rxt(k,396)*y(k,251) &
                      + rxt(k,405)*y(k,253) + rxt(k,379)*y(k,268) + rxt(k,383) &
                      *y(k,269) + .340_r8*rxt(k,512)*y(k,270) + .320_r8*rxt(k,517) &
                      *y(k,271) + .250_r8*rxt(k,453)*y(k,273)
         mat(k,2076) = mat(k,2076) + .500_r8*rxt(k,387)*y(k,18) + rxt(k,429)*y(k,249) &
                      + .250_r8*rxt(k,395)*y(k,251) + rxt(k,406)*y(k,253)
         mat(k,2427) = .340_r8*rxt(k,500)*y(k,6) + rxt(k,338)*y(k,27) &
                      + .500_r8*rxt(k,367)*y(k,31) + .910_r8*rxt(k,445)*y(k,100) &
                      + .120_r8*rxt(k,398)*y(k,111) + .340_r8*rxt(k,503)*y(k,116) &
                      + .600_r8*rxt(k,412)*y(k,118)
         mat(k,587) = rxt(k,362)*y(k,265)
         mat(k,1171) = .680_r8*rxt(k,521)*y(k,265)
         mat(k,1096) = .100_r8*rxt(k,418)*y(k,131)
         mat(k,958) = .700_r8*rxt(k,340)*y(k,241)
         mat(k,988) = rxt(k,368)*y(k,241)
         mat(k,1484) = rxt(k,351)*y(k,241) + rxt(k,425)*y(k,249) + .250_r8*rxt(k,392) &
                      *y(k,251) + rxt(k,401)*y(k,253) + .250_r8*rxt(k,450)*y(k,273)
         mat(k,2363) = rxt(k,241)*y(k,61) + .800_r8*rxt(k,440)*y(k,103) + rxt(k,321) &
                      *y(k,131) + .700_r8*rxt(k,340)*y(k,237) + rxt(k,368)*y(k,238) &
                      + rxt(k,351)*y(k,240) + (4.000_r8*rxt(k,318)+2.000_r8*rxt(k,319)) &
                      *y(k,241) + 1.500_r8*rxt(k,426)*y(k,249) + .750_r8*rxt(k,431) &
                      *y(k,250) + .880_r8*rxt(k,393)*y(k,251) + 2.000_r8*rxt(k,402) &
                      *y(k,253) + .750_r8*rxt(k,505)*y(k,260) + .800_r8*rxt(k,381) &
                      *y(k,269) + .930_r8*rxt(k,510)*y(k,270) + .950_r8*rxt(k,515) &
                      *y(k,271) + .800_r8*rxt(k,451)*y(k,273)
         mat(k,619) = .500_r8*rxt(k,389)*y(k,131)
         mat(k,839) = .500_r8*rxt(k,357)*y(k,131)
         mat(k,2221) = mat(k,2221) + .450_r8*rxt(k,403)*y(k,253) + .150_r8*rxt(k,382) &
                      *y(k,269)
         mat(k,1357) = .920_r8*rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133) + rxt(k,425) &
                      *y(k,240) + 1.500_r8*rxt(k,426)*y(k,241)
         mat(k,1390) = .750_r8*rxt(k,431)*y(k,241)
         mat(k,1411) = .250_r8*rxt(k,396)*y(k,131) + .250_r8*rxt(k,395)*y(k,133) &
                      + .250_r8*rxt(k,392)*y(k,240) + .880_r8*rxt(k,393)*y(k,241)
         mat(k,1452) = rxt(k,405)*y(k,131) + rxt(k,406)*y(k,133) + rxt(k,401)*y(k,240) &
                      + 2.000_r8*rxt(k,402)*y(k,241) + .450_r8*rxt(k,403)*y(k,247) &
                      + 4.000_r8*rxt(k,404)*y(k,253)
         mat(k,1161) = .750_r8*rxt(k,505)*y(k,241)
         mat(k,1843) = (rxt(k,331)+rxt(k,332))*y(k,56)
         mat(k,2016) = mat(k,2016) + .400_r8*rxt(k,416)*y(k,1) + .500_r8*rxt(k,355) &
                      *y(k,53) + rxt(k,322)*y(k,54) + .300_r8*rxt(k,323)*y(k,55) &
                      + .800_r8*rxt(k,360)*y(k,76) + .300_r8*rxt(k,436)*y(k,101) &
                      + .500_r8*rxt(k,411)*y(k,115) + rxt(k,362)*y(k,148) &
                      + .680_r8*rxt(k,521)*y(k,221)
         mat(k,870) = rxt(k,379)*y(k,131)
         mat(k,1269) = rxt(k,383)*y(k,131) + .800_r8*rxt(k,381)*y(k,241) &
                      + .150_r8*rxt(k,382)*y(k,247)
         mat(k,1232) = .340_r8*rxt(k,512)*y(k,131) + .930_r8*rxt(k,510)*y(k,241)
         mat(k,1115) = .320_r8*rxt(k,517)*y(k,131) + .950_r8*rxt(k,515)*y(k,241)
         mat(k,1287) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,450)*y(k,240) &
                      + .800_r8*rxt(k,451)*y(k,241)
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
         mat(k,628) = -(rxt(k,299)*y(k,58) + rxt(k,300)*y(k,265) + rxt(k,310)*y(k,261))
         mat(k,2231) = -rxt(k,299)*y(k,45)
         mat(k,1931) = -rxt(k,300)*y(k,45)
         mat(k,1821) = -rxt(k,310)*y(k,45)
         mat(k,165) = -(rxt(k,301)*y(k,265))
         mat(k,1869) = -rxt(k,301)*y(k,46)
         mat(k,1195) = -(rxt(k,347)*y(k,133) + rxt(k,348)*y(k,265))
         mat(k,2040) = -rxt(k,347)*y(k,47)
         mat(k,1978) = -rxt(k,348)*y(k,47)
         mat(k,703) = .800_r8*rxt(k,416)*y(k,265)
         mat(k,412) = rxt(k,387)*y(k,133)
         mat(k,312) = rxt(k,343)*y(k,265)
         mat(k,380) = .500_r8*rxt(k,344)*y(k,265)
         mat(k,1178) = .500_r8*rxt(k,367)*y(k,143)
         mat(k,1415) = .100_r8*rxt(k,412)*y(k,143)
         mat(k,1765) = .400_r8*rxt(k,418)*y(k,234) + rxt(k,342)*y(k,237) &
                      + .270_r8*rxt(k,370)*y(k,238) + rxt(k,389)*y(k,244) + rxt(k,408) &
                      *y(k,255) + rxt(k,379)*y(k,268)
         mat(k,2040) = mat(k,2040) + rxt(k,387)*y(k,18)
         mat(k,2393) = .500_r8*rxt(k,367)*y(k,31) + .100_r8*rxt(k,412)*y(k,118)
         mat(k,1088) = .400_r8*rxt(k,418)*y(k,131)
         mat(k,951) = rxt(k,342)*y(k,131) + 3.200_r8*rxt(k,339)*y(k,237) &
                      + .800_r8*rxt(k,340)*y(k,241)
         mat(k,981) = .270_r8*rxt(k,370)*y(k,131)
         mat(k,2330) = .800_r8*rxt(k,340)*y(k,237)
         mat(k,614) = rxt(k,389)*y(k,131)
         mat(k,2185) = .200_r8*rxt(k,407)*y(k,255)
         mat(k,722) = rxt(k,408)*y(k,131) + .200_r8*rxt(k,407)*y(k,247)
         mat(k,1978) = mat(k,1978) + .800_r8*rxt(k,416)*y(k,1) + rxt(k,343)*y(k,28) &
                      + .500_r8*rxt(k,344)*y(k,29)
         mat(k,863) = rxt(k,379)*y(k,131)
         mat(k,401) = -(rxt(k,302)*y(k,58) + rxt(k,303)*y(k,265))
         mat(k,2227) = -rxt(k,302)*y(k,48)
         mat(k,1903) = -rxt(k,303)*y(k,48)
         mat(k,147) = -(rxt(k,349)*y(k,265))
         mat(k,1868) = -rxt(k,349)*y(k,49)
         mat(k,1124) = -(rxt(k,386)*y(k,265))
         mat(k,1973) = -rxt(k,386)*y(k,50)
         mat(k,702) = .800_r8*rxt(k,416)*y(k,265)
         mat(k,1067) = .520_r8*rxt(k,500)*y(k,143)
         mat(k,411) = .500_r8*rxt(k,387)*y(k,133)
         mat(k,1011) = .520_r8*rxt(k,503)*y(k,143)
         mat(k,1761) = .250_r8*rxt(k,418)*y(k,234) + .820_r8*rxt(k,370)*y(k,238) &
                      + .500_r8*rxt(k,389)*y(k,244) + .270_r8*rxt(k,512)*y(k,270) &
                      + .040_r8*rxt(k,517)*y(k,271)
         mat(k,2035) = .500_r8*rxt(k,387)*y(k,18)
         mat(k,2389) = .520_r8*rxt(k,500)*y(k,6) + .520_r8*rxt(k,503)*y(k,116)
         mat(k,1162) = .500_r8*rxt(k,521)*y(k,265)
         mat(k,1087) = .250_r8*rxt(k,418)*y(k,131)
         mat(k,980) = .820_r8*rxt(k,370)*y(k,131) + .820_r8*rxt(k,368)*y(k,241)
         mat(k,2326) = .820_r8*rxt(k,368)*y(k,238) + .150_r8*rxt(k,510)*y(k,270) &
                      + .025_r8*rxt(k,515)*y(k,271)
         mat(k,613) = .500_r8*rxt(k,389)*y(k,131)
         mat(k,1973) = mat(k,1973) + .800_r8*rxt(k,416)*y(k,1) + .500_r8*rxt(k,521) &
                      *y(k,221)
         mat(k,1218) = .270_r8*rxt(k,512)*y(k,131) + .150_r8*rxt(k,510)*y(k,241)
         mat(k,1108) = .040_r8*rxt(k,517)*y(k,131) + .025_r8*rxt(k,515)*y(k,241)
         mat(k,1325) = -(rxt(k,373)*y(k,133) + rxt(k,374)*y(k,265))
         mat(k,2050) = -rxt(k,373)*y(k,51)
         mat(k,1988) = -rxt(k,374)*y(k,51)
         mat(k,1253) = rxt(k,376)*y(k,265)
         mat(k,1314) = .880_r8*rxt(k,398)*y(k,143)
         mat(k,1418) = .500_r8*rxt(k,412)*y(k,143)
         mat(k,1775) = .170_r8*rxt(k,471)*y(k,242) + .050_r8*rxt(k,434)*y(k,250) &
                      + .250_r8*rxt(k,396)*y(k,251) + .170_r8*rxt(k,477)*y(k,254) &
                      + .400_r8*rxt(k,487)*y(k,272) + .250_r8*rxt(k,453)*y(k,273) &
                      + .540_r8*rxt(k,493)*y(k,274) + .510_r8*rxt(k,496)*y(k,275)
         mat(k,2050) = mat(k,2050) + .050_r8*rxt(k,435)*y(k,250) + .250_r8*rxt(k,395) &
                      *y(k,251) + .250_r8*rxt(k,454)*y(k,273)
         mat(k,937) = rxt(k,377)*y(k,265)
         mat(k,2401) = .880_r8*rxt(k,398)*y(k,111) + .500_r8*rxt(k,412)*y(k,118)
         mat(k,1466) = .250_r8*rxt(k,392)*y(k,251) + .250_r8*rxt(k,450)*y(k,273)
         mat(k,2339) = .240_r8*rxt(k,393)*y(k,251) + .500_r8*rxt(k,381)*y(k,269) &
                      + .100_r8*rxt(k,451)*y(k,273)
         mat(k,824) = .170_r8*rxt(k,471)*y(k,131) + .070_r8*rxt(k,470)*y(k,247)
         mat(k,2194) = .070_r8*rxt(k,470)*y(k,242) + .070_r8*rxt(k,476)*y(k,254)
         mat(k,1375) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,1399) = .250_r8*rxt(k,396)*y(k,131) + .250_r8*rxt(k,395)*y(k,133) &
                      + .250_r8*rxt(k,392)*y(k,240) + .240_r8*rxt(k,393)*y(k,241)
         mat(k,962) = .170_r8*rxt(k,477)*y(k,131) + .070_r8*rxt(k,476)*y(k,247)
         mat(k,1988) = mat(k,1988) + rxt(k,376)*y(k,97) + rxt(k,377)*y(k,134)
         mat(k,1262) = .500_r8*rxt(k,381)*y(k,241)
         mat(k,800) = .400_r8*rxt(k,487)*y(k,131)
         mat(k,1278) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,240) + .100_r8*rxt(k,451)*y(k,241)
         mat(k,816) = .540_r8*rxt(k,493)*y(k,131)
         mat(k,558) = .510_r8*rxt(k,496)*y(k,131)
         mat(k,749) = -(rxt(k,354)*y(k,265))
         mat(k,1944) = -rxt(k,354)*y(k,52)
         mat(k,1173) = .120_r8*rxt(k,367)*y(k,143)
         mat(k,2377) = .120_r8*rxt(k,367)*y(k,31)
         mat(k,1457) = .100_r8*rxt(k,351)*y(k,241) + .150_r8*rxt(k,352)*y(k,247)
         mat(k,2318) = .100_r8*rxt(k,351)*y(k,240)
         mat(k,2161) = .150_r8*rxt(k,352)*y(k,240) + .150_r8*rxt(k,403)*y(k,253)
         mat(k,1438) = .150_r8*rxt(k,403)*y(k,247)
         mat(k,637) = -(rxt(k,355)*y(k,265))
         mat(k,1932) = -rxt(k,355)*y(k,53)
         mat(k,1456) = .400_r8*rxt(k,352)*y(k,247)
         mat(k,2153) = .400_r8*rxt(k,352)*y(k,240) + .400_r8*rxt(k,403)*y(k,253)
         mat(k,1436) = .400_r8*rxt(k,403)*y(k,247)
         mat(k,857) = -(rxt(k,322)*y(k,265))
         mat(k,1953) = -rxt(k,322)*y(k,54)
         mat(k,1291) = .200_r8*rxt(k,440)*y(k,241)
         mat(k,949) = .300_r8*rxt(k,340)*y(k,241)
         mat(k,2319) = .200_r8*rxt(k,440)*y(k,103) + .300_r8*rxt(k,340)*y(k,237) &
                      + 2.000_r8*rxt(k,319)*y(k,241) + .250_r8*rxt(k,426)*y(k,249) &
                      + .250_r8*rxt(k,431)*y(k,250) + .250_r8*rxt(k,393)*y(k,251) &
                      + .250_r8*rxt(k,505)*y(k,260) + .500_r8*rxt(k,381)*y(k,269) &
                      + .250_r8*rxt(k,510)*y(k,270) + .250_r8*rxt(k,515)*y(k,271) &
                      + .300_r8*rxt(k,451)*y(k,273)
         mat(k,1335) = .250_r8*rxt(k,426)*y(k,241)
         mat(k,1364) = .250_r8*rxt(k,431)*y(k,241)
         mat(k,1393) = .250_r8*rxt(k,393)*y(k,241)
         mat(k,1148) = .250_r8*rxt(k,505)*y(k,241)
         mat(k,1259) = .500_r8*rxt(k,381)*y(k,241)
         mat(k,1217) = .250_r8*rxt(k,510)*y(k,241)
         mat(k,1105) = .250_r8*rxt(k,515)*y(k,241)
         mat(k,1272) = .300_r8*rxt(k,451)*y(k,241)
         mat(k,455) = -(rxt(k,323)*y(k,265))
         mat(k,1911) = -rxt(k,323)*y(k,55)
         mat(k,2315) = rxt(k,320)*y(k,247)
         mat(k,2139) = rxt(k,320)*y(k,241)
         mat(k,2093) = -(rxt(k,235)*y(k,58) + rxt(k,291)*y(k,75) + rxt(k,324)*y(k,265) &
                      + (rxt(k,330) + rxt(k,331) + rxt(k,332)) * y(k,261))
         mat(k,2259) = -rxt(k,235)*y(k,56)
         mat(k,973) = -rxt(k,291)*y(k,56)
         mat(k,2008) = -rxt(k,324)*y(k,56)
         mat(k,1835) = -(rxt(k,330) + rxt(k,331) + rxt(k,332)) * y(k,56)
         mat(k,1189) = .100_r8*rxt(k,367)*y(k,143)
         mat(k,2419) = .100_r8*rxt(k,367)*y(k,31)
         mat(k,461) = -(rxt(k,287)*y(k,261) + rxt(k,304)*y(k,58) + rxt(k,305)*y(k,265))
         mat(k,1819) = -rxt(k,287)*y(k,57)
         mat(k,2228) = -rxt(k,304)*y(k,57)
         mat(k,1912) = -rxt(k,305)*y(k,57)
         mat(k,2261) = -(rxt(k,234)*y(k,44) + rxt(k,235)*y(k,56) + rxt(k,236)*y(k,79) &
                      + rxt(k,237)*y(k,81) + (rxt(k,238) + rxt(k,239)) * y(k,247) &
                      + rxt(k,240)*y(k,143) + rxt(k,247)*y(k,62) + rxt(k,256)*y(k,94) &
                      + rxt(k,297)*y(k,43) + rxt(k,299)*y(k,45) + rxt(k,302)*y(k,48) &
                      + rxt(k,304)*y(k,57) + rxt(k,345)*y(k,30) + rxt(k,375)*y(k,33))
         mat(k,2490) = -rxt(k,234)*y(k,58)
         mat(k,2095) = -rxt(k,235)*y(k,58)
         mat(k,1511) = -rxt(k,236)*y(k,58)
         mat(k,648) = -rxt(k,237)*y(k,58)
         mat(k,2215) = -(rxt(k,238) + rxt(k,239)) * y(k,58)
         mat(k,2421) = -rxt(k,240)*y(k,58)
         mat(k,1033) = -rxt(k,247)*y(k,58)
         mat(k,877) = -rxt(k,256)*y(k,58)
         mat(k,519) = -rxt(k,297)*y(k,58)
         mat(k,634) = -rxt(k,299)*y(k,58)
         mat(k,406) = -rxt(k,302)*y(k,58)
         mat(k,466) = -rxt(k,304)*y(k,58)
         mat(k,338) = -rxt(k,345)*y(k,58)
         mat(k,344) = -rxt(k,375)*y(k,58)
         mat(k,1617) = rxt(k,275)*y(k,61)
         mat(k,146) = 4.000_r8*rxt(k,259)*y(k,261)
         mat(k,189) = rxt(k,260)*y(k,261)
         mat(k,156) = 2.000_r8*rxt(k,261)*y(k,261)
         mat(k,199) = 2.000_r8*rxt(k,262)*y(k,261)
         mat(k,160) = 2.000_r8*rxt(k,263)*y(k,261)
         mat(k,204) = rxt(k,264)*y(k,261)
         mat(k,164) = 2.000_r8*rxt(k,265)*y(k,261)
         mat(k,167) = 3.000_r8*rxt(k,301)*y(k,265)
         mat(k,406) = mat(k,406) + rxt(k,303)*y(k,265)
         mat(k,1643) = rxt(k,275)*y(k,21) + (4.000_r8*rxt(k,242)+2.000_r8*rxt(k,244)) &
                      *y(k,61) + rxt(k,246)*y(k,131) + rxt(k,251)*y(k,141) &
                      + rxt(k,530)*y(k,161) + rxt(k,241)*y(k,241) + rxt(k,252) &
                      *y(k,265)
         mat(k,295) = rxt(k,296)*y(k,261)
         mat(k,291) = rxt(k,311)*y(k,261) + rxt(k,306)*y(k,265)
         mat(k,301) = rxt(k,312)*y(k,261) + rxt(k,307)*y(k,265)
         mat(k,360) = rxt(k,313)*y(k,261) + rxt(k,308)*y(k,265)
         mat(k,1547) = rxt(k,254)*y(k,141) + rxt(k,266)*y(k,261) + rxt(k,255)*y(k,265)
         mat(k,1794) = rxt(k,246)*y(k,61)
         mat(k,2464) = rxt(k,251)*y(k,61) + rxt(k,254)*y(k,87)
         mat(k,1529) = rxt(k,530)*y(k,61)
         mat(k,2357) = rxt(k,241)*y(k,61)
         mat(k,1837) = 4.000_r8*rxt(k,259)*y(k,35) + rxt(k,260)*y(k,36) &
                      + 2.000_r8*rxt(k,261)*y(k,38) + 2.000_r8*rxt(k,262)*y(k,39) &
                      + 2.000_r8*rxt(k,263)*y(k,40) + rxt(k,264)*y(k,41) &
                      + 2.000_r8*rxt(k,265)*y(k,42) + rxt(k,296)*y(k,67) + rxt(k,311) &
                      *y(k,84) + rxt(k,312)*y(k,85) + rxt(k,313)*y(k,86) + rxt(k,266) &
                      *y(k,87)
         mat(k,2010) = 3.000_r8*rxt(k,301)*y(k,46) + rxt(k,303)*y(k,48) + rxt(k,252) &
                      *y(k,61) + rxt(k,306)*y(k,84) + rxt(k,307)*y(k,85) + rxt(k,308) &
                      *y(k,86) + rxt(k,255)*y(k,87)
         mat(k,2223) = rxt(k,247)*y(k,62)
         mat(k,1625) = 2.000_r8*rxt(k,243)*y(k,61)
         mat(k,1025) = rxt(k,247)*y(k,58) + (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,87)
         mat(k,1534) = (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,62) + (rxt(k,587) &
                       +rxt(k,593)+rxt(k,598))*y(k,94)
         mat(k,872) = (rxt(k,587)+rxt(k,593)+rxt(k,598))*y(k,87)
         mat(k,1624) = 2.000_r8*rxt(k,268)*y(k,61)
         mat(k,1636) = -(rxt(k,241)*y(k,241) + (4._r8*rxt(k,242) + 4._r8*rxt(k,243) &
                      + 4._r8*rxt(k,244) + 4._r8*rxt(k,268)) * y(k,61) + rxt(k,245) &
                      *y(k,247) + rxt(k,246)*y(k,131) + rxt(k,248)*y(k,132) + rxt(k,251) &
                      *y(k,141) + (rxt(k,252) + rxt(k,253)) * y(k,265) + (rxt(k,274) &
                      + rxt(k,275) + rxt(k,276)) * y(k,21) + rxt(k,530)*y(k,161))
         mat(k,2349) = -rxt(k,241)*y(k,61)
         mat(k,2207) = -rxt(k,245)*y(k,61)
         mat(k,1786) = -rxt(k,246)*y(k,61)
         mat(k,1680) = -rxt(k,248)*y(k,61)
         mat(k,2456) = -rxt(k,251)*y(k,61)
         mat(k,2002) = -(rxt(k,252) + rxt(k,253)) * y(k,61)
         mat(k,1610) = -(rxt(k,274) + rxt(k,275) + rxt(k,276)) * y(k,61)
         mat(k,1524) = -rxt(k,530)*y(k,61)
         mat(k,2253) = rxt(k,256)*y(k,94) + rxt(k,240)*y(k,143) + rxt(k,239)*y(k,247)
         mat(k,1029) = rxt(k,249)*y(k,141)
         mat(k,1542) = rxt(k,267)*y(k,261)
         mat(k,875) = rxt(k,256)*y(k,58) + rxt(k,257)*y(k,141) + rxt(k,258)*y(k,265)
         mat(k,2456) = mat(k,2456) + rxt(k,249)*y(k,62) + rxt(k,257)*y(k,94)
         mat(k,2413) = rxt(k,240)*y(k,58)
         mat(k,365) = rxt(k,535)*y(k,161)
         mat(k,1524) = mat(k,1524) + rxt(k,535)*y(k,145)
         mat(k,2207) = mat(k,2207) + rxt(k,239)*y(k,58)
         mat(k,1829) = rxt(k,267)*y(k,87)
         mat(k,2002) = mat(k,2002) + rxt(k,258)*y(k,94)
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
         mat(k,1027) = -(rxt(k,247)*y(k,58) + rxt(k,249)*y(k,141) + rxt(k,250) &
                      *y(k,265) + (rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,87))
         mat(k,2238) = -rxt(k,247)*y(k,62)
         mat(k,2446) = -rxt(k,249)*y(k,62)
         mat(k,1966) = -rxt(k,250)*y(k,62)
         mat(k,1538) = -(rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,62)
         mat(k,1630) = rxt(k,248)*y(k,132)
         mat(k,1666) = rxt(k,248)*y(k,61)
         mat(k,1204) = -(rxt(k,334)*y(k,265))
         mat(k,1979) = -rxt(k,334)*y(k,64)
         mat(k,1070) = .230_r8*rxt(k,500)*y(k,143)
         mat(k,1553) = rxt(k,270)*y(k,44)
         mat(k,328) = .350_r8*rxt(k,336)*y(k,265)
         mat(k,623) = .630_r8*rxt(k,338)*y(k,143)
         mat(k,1179) = .560_r8*rxt(k,367)*y(k,143)
         mat(k,2475) = rxt(k,270)*y(k,19) + rxt(k,234)*y(k,58) + rxt(k,315)*y(k,133) &
                      + rxt(k,316)*y(k,141) + rxt(k,317)*y(k,265)
         mat(k,402) = rxt(k,302)*y(k,58)
         mat(k,1324) = rxt(k,373)*y(k,133) + rxt(k,374)*y(k,265)
         mat(k,2242) = rxt(k,234)*y(k,44) + rxt(k,302)*y(k,48)
         mat(k,1493) = rxt(k,615)*y(k,266)
         mat(k,1099) = rxt(k,361)*y(k,265)
         mat(k,923) = .620_r8*rxt(k,445)*y(k,143)
         mat(k,1312) = .650_r8*rxt(k,398)*y(k,143)
         mat(k,1014) = .230_r8*rxt(k,503)*y(k,143)
         mat(k,1416) = .560_r8*rxt(k,412)*y(k,143)
         mat(k,1766) = .170_r8*rxt(k,471)*y(k,242) + .220_r8*rxt(k,396)*y(k,251) &
                      + .400_r8*rxt(k,474)*y(k,252) + .350_r8*rxt(k,477)*y(k,254) &
                      + .225_r8*rxt(k,512)*y(k,270) + .250_r8*rxt(k,453)*y(k,273)
         mat(k,2041) = rxt(k,315)*y(k,44) + rxt(k,373)*y(k,51) + .220_r8*rxt(k,395) &
                      *y(k,251) + .500_r8*rxt(k,454)*y(k,273)
         mat(k,2448) = rxt(k,316)*y(k,44) + rxt(k,524)*y(k,146)
         mat(k,2394) = .230_r8*rxt(k,500)*y(k,6) + .630_r8*rxt(k,338)*y(k,27) &
                      + .560_r8*rxt(k,367)*y(k,31) + .620_r8*rxt(k,445)*y(k,100) &
                      + .650_r8*rxt(k,398)*y(k,111) + .230_r8*rxt(k,503)*y(k,116) &
                      + .560_r8*rxt(k,412)*y(k,118)
         mat(k,420) = rxt(k,524)*y(k,141) + rxt(k,525)*y(k,265)
         mat(k,1164) = .700_r8*rxt(k,521)*y(k,265)
         mat(k,1460) = .220_r8*rxt(k,392)*y(k,251) + .250_r8*rxt(k,450)*y(k,273)
         mat(k,2331) = .110_r8*rxt(k,393)*y(k,251) + .125_r8*rxt(k,510)*y(k,270) &
                      + .200_r8*rxt(k,451)*y(k,273)
         mat(k,823) = .170_r8*rxt(k,471)*y(k,131) + .070_r8*rxt(k,470)*y(k,247)
         mat(k,2186) = .070_r8*rxt(k,470)*y(k,242) + .160_r8*rxt(k,473)*y(k,252) &
                      + .140_r8*rxt(k,476)*y(k,254)
         mat(k,1394) = .220_r8*rxt(k,396)*y(k,131) + .220_r8*rxt(k,395)*y(k,133) &
                      + .220_r8*rxt(k,392)*y(k,240) + .110_r8*rxt(k,393)*y(k,241)
         mat(k,786) = .400_r8*rxt(k,474)*y(k,131) + .160_r8*rxt(k,473)*y(k,247)
         mat(k,961) = .350_r8*rxt(k,477)*y(k,131) + .140_r8*rxt(k,476)*y(k,247)
         mat(k,1979) = mat(k,1979) + .350_r8*rxt(k,336)*y(k,26) + rxt(k,317)*y(k,44) &
                      + rxt(k,374)*y(k,51) + rxt(k,361)*y(k,77) + rxt(k,525)*y(k,146) &
                      + .700_r8*rxt(k,521)*y(k,221)
         mat(k,853) = rxt(k,615)*y(k,65)
         mat(k,1220) = .225_r8*rxt(k,512)*y(k,131) + .125_r8*rxt(k,510)*y(k,241)
         mat(k,1274) = .250_r8*rxt(k,453)*y(k,131) + .500_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,240) + .200_r8*rxt(k,451)*y(k,241)
         mat(k,1494) = -(rxt(k,615)*y(k,266))
         mat(k,854) = -rxt(k,615)*y(k,65)
         mat(k,1074) = .270_r8*rxt(k,500)*y(k,143)
         mat(k,1183) = .200_r8*rxt(k,367)*y(k,143)
         mat(k,750) = rxt(k,354)*y(k,265)
         mat(k,639) = .500_r8*rxt(k,355)*y(k,265)
         mat(k,1205) = rxt(k,334)*y(k,265)
         mat(k,1211) = .800_r8*rxt(k,360)*y(k,265)
         mat(k,1100) = rxt(k,361)*y(k,265)
         mat(k,943) = rxt(k,326)*y(k,265)
         mat(k,674) = .500_r8*rxt(k,411)*y(k,265)
         mat(k,1018) = .270_r8*rxt(k,503)*y(k,143)
         mat(k,1423) = .100_r8*rxt(k,412)*y(k,143)
         mat(k,1782) = rxt(k,353)*y(k,240) + .900_r8*rxt(k,512)*y(k,270)
         mat(k,2408) = .270_r8*rxt(k,500)*y(k,6) + .200_r8*rxt(k,367)*y(k,31) &
                      + .270_r8*rxt(k,503)*y(k,116) + .100_r8*rxt(k,412)*y(k,118)
         mat(k,1167) = 1.800_r8*rxt(k,521)*y(k,265)
         mat(k,1473) = rxt(k,353)*y(k,131) + 4.000_r8*rxt(k,350)*y(k,240) &
                      + .900_r8*rxt(k,351)*y(k,241) + rxt(k,425)*y(k,249) &
                      + 2.000_r8*rxt(k,401)*y(k,253) + rxt(k,450)*y(k,273)
         mat(k,2346) = .900_r8*rxt(k,351)*y(k,240) + rxt(k,402)*y(k,253) &
                      + .500_r8*rxt(k,510)*y(k,270)
         mat(k,2201) = .450_r8*rxt(k,403)*y(k,253)
         mat(k,1348) = rxt(k,425)*y(k,240)
         mat(k,1443) = 2.000_r8*rxt(k,401)*y(k,240) + rxt(k,402)*y(k,241) &
                      + .450_r8*rxt(k,403)*y(k,247) + 4.000_r8*rxt(k,404)*y(k,253)
         mat(k,1995) = rxt(k,354)*y(k,52) + .500_r8*rxt(k,355)*y(k,53) + rxt(k,334) &
                      *y(k,64) + .800_r8*rxt(k,360)*y(k,76) + rxt(k,361)*y(k,77) &
                      + rxt(k,326)*y(k,89) + .500_r8*rxt(k,411)*y(k,115) &
                      + 1.800_r8*rxt(k,521)*y(k,221)
         mat(k,1225) = .900_r8*rxt(k,512)*y(k,131) + .500_r8*rxt(k,510)*y(k,241)
         mat(k,1280) = rxt(k,450)*y(k,240)
         mat(k,267) = -(rxt(k,295)*y(k,261))
         mat(k,1813) = -rxt(k,295)*y(k,66)
         mat(k,186) = rxt(k,260)*y(k,261)
         mat(k,191) = rxt(k,286)*y(k,261)
         mat(k,196) = rxt(k,262)*y(k,261)
         mat(k,158) = 2.000_r8*rxt(k,263)*y(k,261)
         mat(k,201) = 2.000_r8*rxt(k,264)*y(k,261)
         mat(k,162) = rxt(k,265)*y(k,261)
         mat(k,170) = 2.000_r8*rxt(k,288)*y(k,261)
         mat(k,296) = rxt(k,312)*y(k,261) + rxt(k,307)*y(k,265)
         mat(k,355) = rxt(k,313)*y(k,261) + rxt(k,308)*y(k,265)
         mat(k,1813) = mat(k,1813) + rxt(k,260)*y(k,36) + rxt(k,286)*y(k,37) &
                      + rxt(k,262)*y(k,39) + 2.000_r8*rxt(k,263)*y(k,40) &
                      + 2.000_r8*rxt(k,264)*y(k,41) + rxt(k,265)*y(k,42) &
                      + 2.000_r8*rxt(k,288)*y(k,80) + rxt(k,312)*y(k,85) + rxt(k,313) &
                      *y(k,86)
         mat(k,1881) = rxt(k,307)*y(k,85) + rxt(k,308)*y(k,86)
         mat(k,292) = -(rxt(k,296)*y(k,261))
         mat(k,1815) = -rxt(k,296)*y(k,67)
         mat(k,154) = rxt(k,261)*y(k,261)
         mat(k,197) = rxt(k,262)*y(k,261)
         mat(k,288) = rxt(k,311)*y(k,261) + rxt(k,306)*y(k,265)
         mat(k,1815) = mat(k,1815) + rxt(k,261)*y(k,38) + rxt(k,262)*y(k,39) &
                      + rxt(k,311)*y(k,84)
         mat(k,1887) = rxt(k,306)*y(k,84)
         mat(k,236) = -(rxt(k,469)*y(k,265))
         mat(k,1875) = -rxt(k,469)*y(k,68)
         mat(k,230) = .180_r8*rxt(k,489)*y(k,265)
         mat(k,1875) = mat(k,1875) + .180_r8*rxt(k,489)*y(k,223)
         mat(k,349) = -(rxt(k,522)*y(k,133) + (rxt(k,523) + rxt(k,537)) * y(k,265))
         mat(k,2021) = -rxt(k,522)*y(k,69)
         mat(k,1895) = -(rxt(k,523) + rxt(k,537)) * y(k,69)
         mat(k,830) = rxt(k,356)*y(k,247)
         mat(k,2126) = rxt(k,356)*y(k,246)
         mat(k,969) = -(rxt(k,291)*y(k,56) + rxt(k,292)*y(k,79) + rxt(k,293)*y(k,276) &
                      + rxt(k,294)*y(k,91))
         mat(k,2079) = -rxt(k,291)*y(k,75)
         mat(k,1504) = -rxt(k,292)*y(k,75)
         mat(k,2501) = -rxt(k,293)*y(k,75)
         mat(k,2291) = -rxt(k,294)*y(k,75)
         mat(k,192) = rxt(k,286)*y(k,261)
         mat(k,202) = rxt(k,264)*y(k,261)
         mat(k,268) = 2.000_r8*rxt(k,295)*y(k,261)
         mat(k,293) = rxt(k,296)*y(k,261)
         mat(k,1823) = rxt(k,286)*y(k,37) + rxt(k,264)*y(k,41) + 2.000_r8*rxt(k,295) &
                      *y(k,66) + rxt(k,296)*y(k,67)
         mat(k,1210) = -(rxt(k,360)*y(k,265))
         mat(k,1980) = -rxt(k,360)*y(k,76)
         mat(k,652) = .700_r8*rxt(k,436)*y(k,265)
         mat(k,598) = .500_r8*rxt(k,437)*y(k,265)
         mat(k,495) = rxt(k,448)*y(k,265)
         mat(k,1767) = .050_r8*rxt(k,434)*y(k,250) + .530_r8*rxt(k,396)*y(k,251) &
                      + .225_r8*rxt(k,512)*y(k,270) + .250_r8*rxt(k,453)*y(k,273)
         mat(k,2042) = .050_r8*rxt(k,435)*y(k,250) + .530_r8*rxt(k,395)*y(k,251) &
                      + .250_r8*rxt(k,454)*y(k,273)
         mat(k,1582) = rxt(k,359)*y(k,245)
         mat(k,1461) = .530_r8*rxt(k,392)*y(k,251) + .250_r8*rxt(k,450)*y(k,273)
         mat(k,2332) = .260_r8*rxt(k,393)*y(k,251) + .125_r8*rxt(k,510)*y(k,270) &
                      + .100_r8*rxt(k,451)*y(k,273)
         mat(k,510) = rxt(k,359)*y(k,142)
         mat(k,1369) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,1395) = .530_r8*rxt(k,396)*y(k,131) + .530_r8*rxt(k,395)*y(k,133) &
                      + .530_r8*rxt(k,392)*y(k,240) + .260_r8*rxt(k,393)*y(k,241)
         mat(k,1980) = mat(k,1980) + .700_r8*rxt(k,436)*y(k,101) + .500_r8*rxt(k,437) &
                      *y(k,102) + rxt(k,448)*y(k,122)
         mat(k,1221) = .225_r8*rxt(k,512)*y(k,131) + .125_r8*rxt(k,510)*y(k,241)
         mat(k,1275) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,240) + .100_r8*rxt(k,451)*y(k,241)
         mat(k,1098) = -(rxt(k,361)*y(k,265))
         mat(k,1970) = -rxt(k,361)*y(k,77)
         mat(k,327) = .650_r8*rxt(k,336)*y(k,265)
         mat(k,1208) = .200_r8*rxt(k,360)*y(k,265)
         mat(k,1133) = rxt(k,449)*y(k,265)
         mat(k,1758) = rxt(k,460)*y(k,235) + .050_r8*rxt(k,434)*y(k,250) &
                      + .400_r8*rxt(k,474)*y(k,252) + .170_r8*rxt(k,477)*y(k,254) &
                      + .700_r8*rxt(k,480)*y(k,267) + .600_r8*rxt(k,487)*y(k,272) &
                      + .250_r8*rxt(k,453)*y(k,273) + .340_r8*rxt(k,493)*y(k,274) &
                      + .170_r8*rxt(k,496)*y(k,275)
         mat(k,2032) = .050_r8*rxt(k,435)*y(k,250) + .250_r8*rxt(k,454)*y(k,273)
         mat(k,543) = rxt(k,460)*y(k,131)
         mat(k,1458) = .250_r8*rxt(k,450)*y(k,273)
         mat(k,2323) = .100_r8*rxt(k,451)*y(k,273)
         mat(k,2179) = .160_r8*rxt(k,473)*y(k,252) + .070_r8*rxt(k,476)*y(k,254)
         mat(k,1367) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,785) = .400_r8*rxt(k,474)*y(k,131) + .160_r8*rxt(k,473)*y(k,247)
         mat(k,960) = .170_r8*rxt(k,477)*y(k,131) + .070_r8*rxt(k,476)*y(k,247)
         mat(k,1970) = mat(k,1970) + .650_r8*rxt(k,336)*y(k,26) + .200_r8*rxt(k,360) &
                      *y(k,76) + rxt(k,449)*y(k,123)
         mat(k,501) = .700_r8*rxt(k,480)*y(k,131)
         mat(k,798) = .600_r8*rxt(k,487)*y(k,131)
         mat(k,1273) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,240) + .100_r8*rxt(k,451)*y(k,241)
         mat(k,814) = .340_r8*rxt(k,493)*y(k,131)
         mat(k,557) = .170_r8*rxt(k,496)*y(k,131)
         mat(k,2283) = -((rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,247) + rxt(k,195) &
                      *y(k,142) + rxt(k,198)*y(k,143))
         mat(k,2216) = -(rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,78)
         mat(k,1595) = -rxt(k,195)*y(k,78)
         mat(k,2422) = -rxt(k,198)*y(k,78)
         mat(k,2491) = rxt(k,317)*y(k,265)
         mat(k,2096) = rxt(k,331)*y(k,261)
         mat(k,2262) = rxt(k,236)*y(k,79)
         mat(k,974) = rxt(k,292)*y(k,79)
         mat(k,1512) = rxt(k,236)*y(k,58) + rxt(k,292)*y(k,75) + rxt(k,190)*y(k,141) &
                      + rxt(k,173)*y(k,261) + rxt(k,199)*y(k,265)
         mat(k,885) = rxt(k,290)*y(k,261)
         mat(k,1548) = rxt(k,267)*y(k,261)
         mat(k,1051) = rxt(k,222)*y(k,265)
         mat(k,2465) = rxt(k,190)*y(k,79) + rxt(k,202)*y(k,265)
         mat(k,423) = rxt(k,525)*y(k,265)
         mat(k,781) = rxt(k,531)*y(k,265)
         mat(k,1530) = rxt(k,536)*y(k,265)
         mat(k,1838) = rxt(k,331)*y(k,56) + rxt(k,173)*y(k,79) + rxt(k,290)*y(k,83) &
                      + rxt(k,267)*y(k,87)
         mat(k,2011) = rxt(k,317)*y(k,44) + rxt(k,199)*y(k,79) + rxt(k,222)*y(k,119) &
                      + rxt(k,202)*y(k,141) + rxt(k,525)*y(k,146) + rxt(k,531) &
                      *y(k,159) + rxt(k,536)*y(k,161)
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
         mat(k,1505) = -(rxt(k,173)*y(k,261) + rxt(k,190)*y(k,141) + rxt(k,199) &
                      *y(k,265) + rxt(k,236)*y(k,58) + rxt(k,292)*y(k,75))
         mat(k,1824) = -rxt(k,173)*y(k,79)
         mat(k,2450) = -rxt(k,190)*y(k,79)
         mat(k,1996) = -rxt(k,199)*y(k,79)
         mat(k,2248) = -rxt(k,236)*y(k,79)
         mat(k,970) = -rxt(k,292)*y(k,79)
         mat(k,2082) = rxt(k,332)*y(k,261)
         mat(k,2269) = rxt(k,192)*y(k,247)
         mat(k,2202) = rxt(k,192)*y(k,78)
         mat(k,1824) = mat(k,1824) + rxt(k,332)*y(k,56)
         mat(k,169) = -(rxt(k,288)*y(k,261))
         mat(k,1806) = -rxt(k,288)*y(k,80)
         mat(k,644) = -(rxt(k,191)*y(k,141) + rxt(k,200)*y(k,265) + rxt(k,237)*y(k,58))
         mat(k,2435) = -rxt(k,191)*y(k,81)
         mat(k,1933) = -rxt(k,200)*y(k,81)
         mat(k,2232) = -rxt(k,237)*y(k,81)
         mat(k,2154) = 2.000_r8*rxt(k,206)*y(k,247)
         mat(k,1933) = mat(k,1933) + 2.000_r8*rxt(k,205)*y(k,265)
         mat(k,306) = rxt(k,538)*y(k,276)
         mat(k,2498) = rxt(k,538)*y(k,163)
         mat(k,880) = -(rxt(k,283)*y(k,141) + rxt(k,284)*y(k,265) + (rxt(k,289) &
                      + rxt(k,290)) * y(k,261))
         mat(k,2441) = -rxt(k,283)*y(k,83)
         mat(k,1956) = -rxt(k,284)*y(k,83)
         mat(k,1822) = -(rxt(k,289) + rxt(k,290)) * y(k,83)
         mat(k,1552) = rxt(k,270)*y(k,44) + rxt(k,271)*y(k,247)
         mat(k,2473) = rxt(k,270)*y(k,19)
         mat(k,2172) = rxt(k,271)*y(k,19)
         mat(k,287) = -(rxt(k,306)*y(k,265) + rxt(k,311)*y(k,261))
         mat(k,1886) = -rxt(k,306)*y(k,84)
         mat(k,1814) = -rxt(k,311)*y(k,84)
         mat(k,297) = -(rxt(k,307)*y(k,265) + rxt(k,312)*y(k,261))
         mat(k,1888) = -rxt(k,307)*y(k,85)
         mat(k,1816) = -rxt(k,312)*y(k,85)
         mat(k,356) = -(rxt(k,308)*y(k,265) + rxt(k,313)*y(k,261))
         mat(k,1896) = -rxt(k,308)*y(k,86)
         mat(k,1818) = -rxt(k,313)*y(k,86)
         mat(k,1539) = -(rxt(k,254)*y(k,141) + rxt(k,255)*y(k,265) + (rxt(k,266) &
                      + rxt(k,267)) * y(k,261) + (rxt(k,587) + rxt(k,593) + rxt(k,598) &
                      ) * y(k,94) + (rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,62) &
                      + (rxt(k,594) + rxt(k,599)) * y(k,93))
         mat(k,2452) = -rxt(k,254)*y(k,87)
         mat(k,1998) = -rxt(k,255)*y(k,87)
         mat(k,1825) = -(rxt(k,266) + rxt(k,267)) * y(k,87)
         mat(k,874) = -(rxt(k,587) + rxt(k,593) + rxt(k,598)) * y(k,87)
         mat(k,1028) = -(rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,87)
         mat(k,842) = -(rxt(k,594) + rxt(k,599)) * y(k,87)
         mat(k,336) = rxt(k,345)*y(k,58)
         mat(k,342) = rxt(k,375)*y(k,58)
         mat(k,515) = rxt(k,297)*y(k,58)
         mat(k,2478) = rxt(k,234)*y(k,58)
         mat(k,629) = rxt(k,299)*y(k,58)
         mat(k,403) = 2.000_r8*rxt(k,302)*y(k,58)
         mat(k,2083) = rxt(k,235)*y(k,58)
         mat(k,462) = rxt(k,304)*y(k,58)
         mat(k,2249) = rxt(k,345)*y(k,30) + rxt(k,375)*y(k,33) + rxt(k,297)*y(k,43) &
                      + rxt(k,234)*y(k,44) + rxt(k,299)*y(k,45) + 2.000_r8*rxt(k,302) &
                      *y(k,48) + rxt(k,235)*y(k,56) + rxt(k,304)*y(k,57) + rxt(k,236) &
                      *y(k,79) + rxt(k,237)*y(k,81) + rxt(k,256)*y(k,94) + rxt(k,238) &
                      *y(k,247)
         mat(k,1632) = rxt(k,253)*y(k,265)
         mat(k,1506) = rxt(k,236)*y(k,58)
         mat(k,645) = rxt(k,237)*y(k,58)
         mat(k,874) = mat(k,874) + rxt(k,256)*y(k,58)
         mat(k,2203) = rxt(k,238)*y(k,58)
         mat(k,1998) = mat(k,1998) + rxt(k,253)*y(k,61)
         mat(k,248) = -(rxt(k,325)*y(k,265) + rxt(k,333)*y(k,261))
         mat(k,1878) = -rxt(k,325)*y(k,88)
         mat(k,1812) = -rxt(k,333)*y(k,88)
         mat(k,942) = -(rxt(k,326)*y(k,265))
         mat(k,1959) = -rxt(k,326)*y(k,89)
         mat(k,1061) = .050_r8*rxt(k,500)*y(k,143)
         mat(k,326) = .350_r8*rxt(k,336)*y(k,265)
         mat(k,622) = .370_r8*rxt(k,338)*y(k,143)
         mat(k,1176) = .120_r8*rxt(k,367)*y(k,143)
         mat(k,921) = .110_r8*rxt(k,445)*y(k,143)
         mat(k,1311) = .330_r8*rxt(k,398)*y(k,143)
         mat(k,1005) = .050_r8*rxt(k,503)*y(k,143)
         mat(k,1413) = .120_r8*rxt(k,412)*y(k,143)
         mat(k,1752) = rxt(k,329)*y(k,248)
         mat(k,2381) = .050_r8*rxt(k,500)*y(k,6) + .370_r8*rxt(k,338)*y(k,27) &
                      + .120_r8*rxt(k,367)*y(k,31) + .110_r8*rxt(k,445)*y(k,100) &
                      + .330_r8*rxt(k,398)*y(k,111) + .050_r8*rxt(k,503)*y(k,116) &
                      + .120_r8*rxt(k,412)*y(k,118)
         mat(k,2174) = rxt(k,327)*y(k,248)
         mat(k,488) = rxt(k,329)*y(k,131) + rxt(k,327)*y(k,247)
         mat(k,1959) = mat(k,1959) + .350_r8*rxt(k,336)*y(k,26)
         mat(k,2078) = rxt(k,291)*y(k,75)
         mat(k,968) = rxt(k,291)*y(k,56) + rxt(k,292)*y(k,79) + rxt(k,294)*y(k,91) &
                      + rxt(k,293)*y(k,276)
         mat(k,1503) = rxt(k,292)*y(k,75)
         mat(k,2290) = rxt(k,294)*y(k,75)
         mat(k,2500) = rxt(k,293)*y(k,75)
         mat(k,2307) = -(rxt(k,231)*y(k,265) + rxt(k,294)*y(k,75))
         mat(k,2012) = -rxt(k,231)*y(k,91)
         mat(k,975) = -rxt(k,294)*y(k,91)
         mat(k,2492) = rxt(k,315)*y(k,133)
         mat(k,1201) = rxt(k,347)*y(k,133)
         mat(k,1331) = rxt(k,373)*y(k,133)
         mat(k,1034) = (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,87)
         mat(k,354) = rxt(k,522)*y(k,133)
         mat(k,1549) = (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,62)
         mat(k,1690) = rxt(k,230)*y(k,265)
         mat(k,2072) = rxt(k,315)*y(k,44) + rxt(k,347)*y(k,47) + rxt(k,373)*y(k,51) &
                      + rxt(k,522)*y(k,69)
         mat(k,2012) = mat(k,2012) + rxt(k,230)*y(k,132)
         mat(k,521) = -(rxt(k,207)*y(k,265))
         mat(k,1919) = -rxt(k,207)*y(k,92)
         mat(k,1654) = rxt(k,228)*y(k,247)
         mat(k,2146) = rxt(k,228)*y(k,132)
         mat(k,841) = -(rxt(k,285)*y(k,141) + (rxt(k,594) + rxt(k,599)) * y(k,87))
         mat(k,2438) = -rxt(k,285)*y(k,93)
         mat(k,1536) = -(rxt(k,594) + rxt(k,599)) * y(k,93)
         mat(k,1603) = rxt(k,277)*y(k,247)
         mat(k,2169) = rxt(k,277)*y(k,21)
         mat(k,873) = -(rxt(k,256)*y(k,58) + rxt(k,257)*y(k,141) + rxt(k,258)*y(k,265) &
                      + (rxt(k,587) + rxt(k,593) + rxt(k,598)) * y(k,87))
         mat(k,2234) = -rxt(k,256)*y(k,94)
         mat(k,2440) = -rxt(k,257)*y(k,94)
         mat(k,1955) = -rxt(k,258)*y(k,94)
         mat(k,1537) = -(rxt(k,587) + rxt(k,593) + rxt(k,598)) * y(k,94)
         mat(k,1628) = rxt(k,245)*y(k,247)
         mat(k,1026) = rxt(k,250)*y(k,265)
         mat(k,2171) = rxt(k,245)*y(k,61)
         mat(k,1955) = mat(k,1955) + rxt(k,250)*y(k,62)
         mat(k,1239) = -(rxt(k,391)*y(k,265))
         mat(k,1982) = -rxt(k,391)*y(k,95)
         mat(k,653) = .300_r8*rxt(k,436)*y(k,265)
         mat(k,599) = .500_r8*rxt(k,437)*y(k,265)
         mat(k,1769) = rxt(k,390)*y(k,244) + rxt(k,397)*y(k,251)
         mat(k,615) = rxt(k,390)*y(k,131)
         mat(k,1396) = rxt(k,397)*y(k,131)
         mat(k,1982) = mat(k,1982) + .300_r8*rxt(k,436)*y(k,101) + .500_r8*rxt(k,437) &
                      *y(k,102)
         mat(k,282) = -(rxt(k,422)*y(k,265))
         mat(k,1885) = -rxt(k,422)*y(k,96)
         mat(k,1252) = -(rxt(k,376)*y(k,265))
         mat(k,1983) = -rxt(k,376)*y(k,97)
         mat(k,654) = .700_r8*rxt(k,436)*y(k,265)
         mat(k,600) = .500_r8*rxt(k,437)*y(k,265)
         mat(k,672) = .500_r8*rxt(k,411)*y(k,265)
         mat(k,1770) = .050_r8*rxt(k,434)*y(k,250) + .220_r8*rxt(k,396)*y(k,251) &
                      + .250_r8*rxt(k,453)*y(k,273)
         mat(k,2045) = .050_r8*rxt(k,435)*y(k,250) + .220_r8*rxt(k,395)*y(k,251) &
                      + .250_r8*rxt(k,454)*y(k,273)
         mat(k,591) = .500_r8*rxt(k,380)*y(k,265)
         mat(k,1462) = .220_r8*rxt(k,392)*y(k,251) + .250_r8*rxt(k,450)*y(k,273)
         mat(k,2334) = .230_r8*rxt(k,393)*y(k,251) + .200_r8*rxt(k,381)*y(k,269) &
                      + .100_r8*rxt(k,451)*y(k,273)
         mat(k,1371) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,1397) = .220_r8*rxt(k,396)*y(k,131) + .220_r8*rxt(k,395)*y(k,133) &
                      + .220_r8*rxt(k,392)*y(k,240) + .230_r8*rxt(k,393)*y(k,241)
         mat(k,1983) = mat(k,1983) + .700_r8*rxt(k,436)*y(k,101) + .500_r8*rxt(k,437) &
                      *y(k,102) + .500_r8*rxt(k,411)*y(k,115) + .500_r8*rxt(k,380) &
                      *y(k,157)
         mat(k,1260) = .200_r8*rxt(k,381)*y(k,241)
         mat(k,1276) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,240) + .100_r8*rxt(k,451)*y(k,241)
         mat(k,398) = -(rxt(k,423)*y(k,265))
         mat(k,1902) = -rxt(k,423)*y(k,98)
         mat(k,1723) = .870_r8*rxt(k,434)*y(k,250)
         mat(k,2022) = .950_r8*rxt(k,435)*y(k,250)
         mat(k,1454) = rxt(k,430)*y(k,250)
         mat(k,2314) = .750_r8*rxt(k,431)*y(k,250)
         mat(k,1360) = .870_r8*rxt(k,434)*y(k,131) + .950_r8*rxt(k,435)*y(k,133) &
                      + rxt(k,430)*y(k,240) + .750_r8*rxt(k,431)*y(k,241)
         mat(k,173) = -(rxt(k,424)*y(k,265))
         mat(k,1870) = -rxt(k,424)*y(k,99)
         mat(k,754) = .600_r8*rxt(k,447)*y(k,265)
         mat(k,1870) = mat(k,1870) + .600_r8*rxt(k,447)*y(k,106)
         mat(k,920) = -(rxt(k,438)*y(k,133) + rxt(k,445)*y(k,143) + rxt(k,446) &
                      *y(k,265))
         mat(k,2025) = -rxt(k,438)*y(k,100)
         mat(k,2380) = -rxt(k,445)*y(k,100)
         mat(k,1957) = -rxt(k,446)*y(k,100)
         mat(k,651) = -(rxt(k,436)*y(k,265))
         mat(k,1934) = -rxt(k,436)*y(k,101)
         mat(k,1736) = .080_r8*rxt(k,428)*y(k,249)
         mat(k,1333) = .080_r8*rxt(k,428)*y(k,131)
         mat(k,596) = -(rxt(k,437)*y(k,265))
         mat(k,1928) = -rxt(k,437)*y(k,102)
         mat(k,1734) = .080_r8*rxt(k,434)*y(k,250)
         mat(k,1361) = .080_r8*rxt(k,434)*y(k,131)
         mat(k,1297) = -(rxt(k,439)*y(k,240) + rxt(k,440)*y(k,241) + rxt(k,441) &
                      *y(k,247) + rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133))
         mat(k,1464) = -rxt(k,439)*y(k,103)
         mat(k,2337) = -rxt(k,440)*y(k,103)
         mat(k,2192) = -rxt(k,441)*y(k,103)
         mat(k,1773) = -rxt(k,442)*y(k,103)
         mat(k,2048) = -rxt(k,443)*y(k,103)
         mat(k,924) = rxt(k,438)*y(k,133)
         mat(k,2048) = mat(k,2048) + rxt(k,438)*y(k,100)
         mat(k,443) = -(rxt(k,444)*y(k,265))
         mat(k,1909) = -rxt(k,444)*y(k,104)
         mat(k,1288) = rxt(k,441)*y(k,247)
         mat(k,2137) = rxt(k,441)*y(k,103)
         mat(k,84) = -(rxt(k,562)*y(k,247) + rxt(k,563)*y(k,131))
         mat(k,2114) = -rxt(k,562)*y(k,105)
         mat(k,1707) = -rxt(k,563)*y(k,105)
         mat(k,919) = rxt(k,565)*y(k,265)
         mat(k,1852) = rxt(k,565)*y(k,100)
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
         mat(k,755) = -(rxt(k,447)*y(k,265))
         mat(k,1945) = -rxt(k,447)*y(k,106)
         mat(k,2162) = rxt(k,427)*y(k,249) + rxt(k,432)*y(k,250)
         mat(k,1334) = rxt(k,427)*y(k,247)
         mat(k,1363) = rxt(k,432)*y(k,247)
         mat(k,73) = -(rxt(k,568)*y(k,265))
         mat(k,1851) = -rxt(k,568)*y(k,107)
         mat(k,71) = -(rxt(k,566)*y(k,247) + rxt(k,567)*y(k,131))
         mat(k,2108) = -rxt(k,566)*y(k,108)
         mat(k,1701) = -rxt(k,567)*y(k,108)
         mat(k,72) = rxt(k,568)*y(k,265)
         mat(k,1850) = rxt(k,568)*y(k,107)
         mat(k,115) = -(rxt(k,571)*y(k,265))
         mat(k,1862) = -rxt(k,571)*y(k,109)
         mat(k,113) = -(rxt(k,569)*y(k,247) + rxt(k,570)*y(k,131))
         mat(k,2123) = -rxt(k,569)*y(k,110)
         mat(k,1716) = -rxt(k,570)*y(k,110)
         mat(k,114) = rxt(k,571)*y(k,265)
         mat(k,1861) = rxt(k,571)*y(k,109)
         mat(k,1313) = -(rxt(k,398)*y(k,143) + rxt(k,399)*y(k,265))
         mat(k,2400) = -rxt(k,398)*y(k,111)
         mat(k,1987) = -rxt(k,399)*y(k,111)
         mat(k,925) = .300_r8*rxt(k,445)*y(k,143)
         mat(k,1774) = .360_r8*rxt(k,428)*y(k,249)
         mat(k,2049) = .400_r8*rxt(k,429)*y(k,249)
         mat(k,2400) = mat(k,2400) + .300_r8*rxt(k,445)*y(k,100)
         mat(k,1465) = .390_r8*rxt(k,425)*y(k,249)
         mat(k,2338) = .310_r8*rxt(k,426)*y(k,249)
         mat(k,1341) = .360_r8*rxt(k,428)*y(k,131) + .400_r8*rxt(k,429)*y(k,133) &
                      + .390_r8*rxt(k,425)*y(k,240) + .310_r8*rxt(k,426)*y(k,241)
         mat(k,393) = -(rxt(k,400)*y(k,265))
         mat(k,1901) = -rxt(k,400)*y(k,112)
         mat(k,2133) = rxt(k,394)*y(k,251)
         mat(k,1392) = rxt(k,394)*y(k,247)
         mat(k,563) = -(rxt(k,409)*y(k,265))
         mat(k,1924) = -rxt(k,409)*y(k,113)
         mat(k,1732) = .800_r8*rxt(k,418)*y(k,234)
         mat(k,1081) = .800_r8*rxt(k,418)*y(k,131)
         mat(k,388) = -(rxt(k,410)*y(k,265))
         mat(k,1900) = -rxt(k,410)*y(k,114)
         mat(k,2132) = .800_r8*rxt(k,407)*y(k,255)
         mat(k,720) = .800_r8*rxt(k,407)*y(k,247)
         mat(k,671) = -(rxt(k,411)*y(k,265))
         mat(k,1936) = -rxt(k,411)*y(k,115)
         mat(k,1659) = rxt(k,414)*y(k,253)
         mat(k,1437) = rxt(k,414)*y(k,132)
         mat(k,1006) = -(rxt(k,502)*y(k,133) + rxt(k,503)*y(k,143) + rxt(k,504) &
                      *y(k,265))
         mat(k,2029) = -rxt(k,502)*y(k,116)
         mat(k,2383) = -rxt(k,503)*y(k,116)
         mat(k,1965) = -rxt(k,504)*y(k,116)
         mat(k,90) = -(rxt(k,573)*y(k,247) + rxt(k,574)*y(k,131))
         mat(k,2115) = -rxt(k,573)*y(k,117)
         mat(k,1708) = -rxt(k,574)*y(k,117)
         mat(k,1002) = rxt(k,576)*y(k,265)
         mat(k,1853) = rxt(k,576)*y(k,116)
         mat(k,1420) = -(rxt(k,412)*y(k,143) + rxt(k,413)*y(k,265))
         mat(k,2405) = -rxt(k,412)*y(k,118)
         mat(k,1992) = -rxt(k,413)*y(k,118)
         mat(k,928) = .200_r8*rxt(k,445)*y(k,143)
         mat(k,1779) = .560_r8*rxt(k,428)*y(k,249)
         mat(k,2054) = .600_r8*rxt(k,429)*y(k,249)
         mat(k,2405) = mat(k,2405) + .200_r8*rxt(k,445)*y(k,100)
         mat(k,1470) = .610_r8*rxt(k,425)*y(k,249)
         mat(k,2343) = .440_r8*rxt(k,426)*y(k,249)
         mat(k,1345) = .560_r8*rxt(k,428)*y(k,131) + .600_r8*rxt(k,429)*y(k,133) &
                      + .610_r8*rxt(k,425)*y(k,240) + .440_r8*rxt(k,426)*y(k,241)
         mat(k,1043) = -(rxt(k,210)*y(k,131) + (rxt(k,211) + rxt(k,212) + rxt(k,213) &
                      ) * y(k,132) + rxt(k,214)*y(k,142) + rxt(k,222)*y(k,265) &
                      + rxt(k,612)*y(k,264))
         mat(k,1756) = -rxt(k,210)*y(k,119)
         mat(k,1667) = -(rxt(k,211) + rxt(k,212) + rxt(k,213)) * y(k,119)
         mat(k,1580) = -rxt(k,214)*y(k,119)
         mat(k,1967) = -rxt(k,222)*y(k,119)
         mat(k,899) = -rxt(k,612)*y(k,119)
         mat(k,2447) = rxt(k,208)*y(k,256) + rxt(k,609)*y(k,259)
         mat(k,1580) = mat(k,1580) + rxt(k,610)*y(k,259)
         mat(k,910) = 1.100_r8*rxt(k,605)*y(k,257) + .200_r8*rxt(k,603)*y(k,258)
         mat(k,576) = rxt(k,208)*y(k,141)
         mat(k,744) = 1.100_r8*rxt(k,605)*y(k,243)
         mat(k,891) = .200_r8*rxt(k,603)*y(k,243)
         mat(k,552) = rxt(k,609)*y(k,141) + rxt(k,610)*y(k,142)
         mat(k,302) = -((rxt(k,226) + rxt(k,227)) * y(k,261))
         mat(k,1817) = -(rxt(k,226) + rxt(k,227)) * y(k,120)
         mat(k,1037) = rxt(k,211)*y(k,132)
         mat(k,1652) = rxt(k,211)*y(k,119)
         mat(k,1653) = rxt(k,229)*y(k,133)
         mat(k,2020) = rxt(k,229)*y(k,132)
         mat(k,493) = -(rxt(k,448)*y(k,265))
         mat(k,1915) = -rxt(k,448)*y(k,122)
         mat(k,1289) = .200_r8*rxt(k,440)*y(k,241)
         mat(k,2316) = .200_r8*rxt(k,440)*y(k,103)
         mat(k,1134) = -(rxt(k,449)*y(k,265))
         mat(k,1974) = -rxt(k,449)*y(k,123)
         mat(k,1293) = rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133) + rxt(k,439)*y(k,240) &
                      + .800_r8*rxt(k,440)*y(k,241)
         mat(k,1762) = rxt(k,442)*y(k,103)
         mat(k,2036) = rxt(k,443)*y(k,103)
         mat(k,1459) = rxt(k,439)*y(k,103)
         mat(k,2327) = .800_r8*rxt(k,440)*y(k,103)
         mat(k,141) = -(rxt(k,539)*y(k,265))
         mat(k,1867) = -rxt(k,539)*y(k,127)
         mat(k,1788) = -(rxt(k,210)*y(k,119) + rxt(k,219)*y(k,133) + rxt(k,223) &
                      *y(k,247) + rxt(k,224)*y(k,143) + rxt(k,225)*y(k,141) + rxt(k,246) &
                      *y(k,61) + rxt(k,278)*y(k,21) + rxt(k,321)*y(k,241) + rxt(k,329) &
                      *y(k,248) + rxt(k,342)*y(k,237) + rxt(k,353)*y(k,240) + rxt(k,357) &
                      *y(k,246) + rxt(k,370)*y(k,238) + rxt(k,379)*y(k,268) + rxt(k,383) &
                      *y(k,269) + (rxt(k,389) + rxt(k,390)) * y(k,244) + (rxt(k,396) &
                      + rxt(k,397)) * y(k,251) + rxt(k,405)*y(k,253) + rxt(k,408) &
                      *y(k,255) + (rxt(k,418) + rxt(k,419)) * y(k,234) + rxt(k,428) &
                      *y(k,249) + rxt(k,434)*y(k,250) + rxt(k,442)*y(k,103) + rxt(k,453) &
                      *y(k,273) + rxt(k,457)*y(k,233) + rxt(k,460)*y(k,235) + rxt(k,465) &
                      *y(k,236) + rxt(k,467)*y(k,239) + rxt(k,471)*y(k,242) + rxt(k,474) &
                      *y(k,252) + rxt(k,477)*y(k,254) + rxt(k,480)*y(k,267) + rxt(k,487) &
                      *y(k,272) + rxt(k,493)*y(k,274) + rxt(k,496)*y(k,275) + rxt(k,507) &
                      *y(k,260) + rxt(k,512)*y(k,270) + rxt(k,517)*y(k,271) + rxt(k,614) &
                      *y(k,264))
         mat(k,1048) = -rxt(k,210)*y(k,131)
         mat(k,2064) = -rxt(k,219)*y(k,131)
         mat(k,2209) = -rxt(k,223)*y(k,131)
         mat(k,2415) = -rxt(k,224)*y(k,131)
         mat(k,2458) = -rxt(k,225)*y(k,131)
         mat(k,1638) = -rxt(k,246)*y(k,131)
         mat(k,1612) = -rxt(k,278)*y(k,131)
         mat(k,2351) = -rxt(k,321)*y(k,131)
         mat(k,490) = -rxt(k,329)*y(k,131)
         mat(k,954) = -rxt(k,342)*y(k,131)
         mat(k,1476) = -rxt(k,353)*y(k,131)
         mat(k,836) = -rxt(k,357)*y(k,131)
         mat(k,984) = -rxt(k,370)*y(k,131)
         mat(k,867) = -rxt(k,379)*y(k,131)
         mat(k,1265) = -rxt(k,383)*y(k,131)
         mat(k,617) = -(rxt(k,389) + rxt(k,390)) * y(k,131)
         mat(k,1405) = -(rxt(k,396) + rxt(k,397)) * y(k,131)
         mat(k,1445) = -rxt(k,405)*y(k,131)
         mat(k,725) = -rxt(k,408)*y(k,131)
         mat(k,1092) = -(rxt(k,418) + rxt(k,419)) * y(k,131)
         mat(k,1350) = -rxt(k,428)*y(k,131)
         mat(k,1383) = -rxt(k,434)*y(k,131)
         mat(k,1304) = -rxt(k,442)*y(k,131)
         mat(k,1282) = -rxt(k,453)*y(k,131)
         mat(k,572) = -rxt(k,457)*y(k,131)
         mat(k,545) = -rxt(k,460)*y(k,131)
         mat(k,484) = -rxt(k,465)*y(k,131)
         mat(k,694) = -rxt(k,467)*y(k,131)
         mat(k,826) = -rxt(k,471)*y(k,131)
         mat(k,788) = -rxt(k,474)*y(k,131)
         mat(k,964) = -rxt(k,477)*y(k,131)
         mat(k,503) = -rxt(k,480)*y(k,131)
         mat(k,802) = -rxt(k,487)*y(k,131)
         mat(k,819) = -rxt(k,493)*y(k,131)
         mat(k,560) = -rxt(k,496)*y(k,131)
         mat(k,1155) = -rxt(k,507)*y(k,131)
         mat(k,1227) = -rxt(k,512)*y(k,131)
         mat(k,1110) = -rxt(k,517)*y(k,131)
         mat(k,901) = -rxt(k,614)*y(k,131)
         mat(k,1048) = mat(k,1048) + 2.000_r8*rxt(k,212)*y(k,132) + rxt(k,214) &
                      *y(k,142) + rxt(k,222)*y(k,265)
         mat(k,304) = 2.000_r8*rxt(k,226)*y(k,261)
         mat(k,1682) = 2.000_r8*rxt(k,212)*y(k,119) + rxt(k,215)*y(k,141) + rxt(k,532) &
                      *y(k,161)
         mat(k,2458) = mat(k,2458) + rxt(k,215)*y(k,132)
         mat(k,1590) = rxt(k,214)*y(k,119) + rxt(k,209)*y(k,256)
         mat(k,1526) = rxt(k,532)*y(k,132)
         mat(k,578) = rxt(k,209)*y(k,142)
         mat(k,1831) = 2.000_r8*rxt(k,226)*y(k,120)
         mat(k,2004) = rxt(k,222)*y(k,119)
         mat(k,1681) = -((rxt(k,211) + rxt(k,212) + rxt(k,213)) * y(k,119) + (rxt(k,215) &
                      + rxt(k,217)) * y(k,141) + rxt(k,216)*y(k,143) + rxt(k,228) &
                      *y(k,247) + rxt(k,229)*y(k,133) + rxt(k,230)*y(k,265) + rxt(k,248) &
                      *y(k,61) + rxt(k,279)*y(k,21) + rxt(k,364)*y(k,240) + rxt(k,414) &
                      *y(k,253) + rxt(k,472)*y(k,242) + rxt(k,475)*y(k,252) + rxt(k,478) &
                      *y(k,254) + rxt(k,482)*y(k,150) + rxt(k,485)*y(k,233) + rxt(k,532) &
                      *y(k,161))
         mat(k,1047) = -(rxt(k,211) + rxt(k,212) + rxt(k,213)) * y(k,132)
         mat(k,2457) = -(rxt(k,215) + rxt(k,217)) * y(k,132)
         mat(k,2414) = -rxt(k,216)*y(k,132)
         mat(k,2208) = -rxt(k,228)*y(k,132)
         mat(k,2063) = -rxt(k,229)*y(k,132)
         mat(k,2003) = -rxt(k,230)*y(k,132)
         mat(k,1637) = -rxt(k,248)*y(k,132)
         mat(k,1611) = -rxt(k,279)*y(k,132)
         mat(k,1475) = -rxt(k,364)*y(k,132)
         mat(k,1444) = -rxt(k,414)*y(k,132)
         mat(k,825) = -rxt(k,472)*y(k,132)
         mat(k,787) = -rxt(k,475)*y(k,132)
         mat(k,963) = -rxt(k,478)*y(k,132)
         mat(k,537) = -rxt(k,482)*y(k,132)
         mat(k,571) = -rxt(k,485)*y(k,132)
         mat(k,1525) = -rxt(k,532)*y(k,132)
         mat(k,705) = rxt(k,416)*y(k,265)
         mat(k,413) = rxt(k,387)*y(k,133)
         mat(k,1611) = mat(k,1611) + rxt(k,278)*y(k,131)
         mat(k,1637) = mat(k,1637) + rxt(k,246)*y(k,131)
         mat(k,523) = rxt(k,207)*y(k,265)
         mat(k,655) = .700_r8*rxt(k,436)*y(k,265)
         mat(k,1303) = rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133)
         mat(k,1787) = rxt(k,278)*y(k,21) + rxt(k,246)*y(k,61) + rxt(k,442)*y(k,103) &
                      + 2.000_r8*rxt(k,219)*y(k,133) + rxt(k,225)*y(k,141) &
                      + rxt(k,224)*y(k,143) + rxt(k,457)*y(k,233) + rxt(k,418) &
                      *y(k,234) + rxt(k,460)*y(k,235) + rxt(k,465)*y(k,236) &
                      + rxt(k,342)*y(k,237) + rxt(k,370)*y(k,238) + rxt(k,467) &
                      *y(k,239) + rxt(k,353)*y(k,240) + rxt(k,321)*y(k,241) &
                      + rxt(k,471)*y(k,242) + rxt(k,389)*y(k,244) + rxt(k,357) &
                      *y(k,246) + rxt(k,223)*y(k,247) + rxt(k,329)*y(k,248) &
                      + .920_r8*rxt(k,428)*y(k,249) + .920_r8*rxt(k,434)*y(k,250) &
                      + rxt(k,396)*y(k,251) + rxt(k,474)*y(k,252) + rxt(k,405) &
                      *y(k,253) + rxt(k,477)*y(k,254) + rxt(k,408)*y(k,255) &
                      + 1.600_r8*rxt(k,507)*y(k,260) + rxt(k,480)*y(k,267) &
                      + rxt(k,379)*y(k,268) + rxt(k,383)*y(k,269) + .900_r8*rxt(k,512) &
                      *y(k,270) + .800_r8*rxt(k,517)*y(k,271) + rxt(k,487)*y(k,272) &
                      + rxt(k,453)*y(k,273) + rxt(k,493)*y(k,274) + rxt(k,496) &
                      *y(k,275)
         mat(k,2063) = mat(k,2063) + rxt(k,387)*y(k,18) + rxt(k,443)*y(k,103) &
                      + 2.000_r8*rxt(k,219)*y(k,131) + rxt(k,220)*y(k,141) &
                      + rxt(k,218)*y(k,247) + rxt(k,429)*y(k,249) + rxt(k,435) &
                      *y(k,250) + rxt(k,395)*y(k,251) + rxt(k,406)*y(k,253) &
                      + 2.000_r8*rxt(k,508)*y(k,260) + rxt(k,221)*y(k,265) &
                      + rxt(k,454)*y(k,273)
         mat(k,939) = rxt(k,377)*y(k,265)
         mat(k,2457) = mat(k,2457) + rxt(k,225)*y(k,131) + rxt(k,220)*y(k,133)
         mat(k,2414) = mat(k,2414) + rxt(k,224)*y(k,131)
         mat(k,686) = rxt(k,514)*y(k,265)
         mat(k,571) = mat(k,571) + rxt(k,457)*y(k,131)
         mat(k,1091) = rxt(k,418)*y(k,131)
         mat(k,544) = rxt(k,460)*y(k,131)
         mat(k,483) = rxt(k,465)*y(k,131)
         mat(k,953) = rxt(k,342)*y(k,131)
         mat(k,983) = rxt(k,370)*y(k,131)
         mat(k,693) = rxt(k,467)*y(k,131)
         mat(k,1475) = mat(k,1475) + rxt(k,353)*y(k,131)
         mat(k,2350) = rxt(k,321)*y(k,131) + .500_r8*rxt(k,505)*y(k,260)
         mat(k,825) = mat(k,825) + rxt(k,471)*y(k,131)
         mat(k,616) = rxt(k,389)*y(k,131)
         mat(k,835) = rxt(k,357)*y(k,131)
         mat(k,2208) = mat(k,2208) + rxt(k,223)*y(k,131) + rxt(k,218)*y(k,133)
         mat(k,489) = rxt(k,329)*y(k,131)
         mat(k,1349) = .920_r8*rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133)
         mat(k,1382) = .920_r8*rxt(k,434)*y(k,131) + rxt(k,435)*y(k,133)
         mat(k,1404) = rxt(k,396)*y(k,131) + rxt(k,395)*y(k,133)
         mat(k,787) = mat(k,787) + rxt(k,474)*y(k,131)
         mat(k,1444) = mat(k,1444) + rxt(k,405)*y(k,131) + rxt(k,406)*y(k,133)
         mat(k,963) = mat(k,963) + rxt(k,477)*y(k,131)
         mat(k,724) = rxt(k,408)*y(k,131)
         mat(k,1154) = 1.600_r8*rxt(k,507)*y(k,131) + 2.000_r8*rxt(k,508)*y(k,133) &
                      + .500_r8*rxt(k,505)*y(k,241)
         mat(k,2003) = mat(k,2003) + rxt(k,416)*y(k,1) + rxt(k,207)*y(k,92) &
                      + .700_r8*rxt(k,436)*y(k,101) + rxt(k,221)*y(k,133) + rxt(k,377) &
                      *y(k,134) + rxt(k,514)*y(k,218)
         mat(k,502) = rxt(k,480)*y(k,131)
         mat(k,866) = rxt(k,379)*y(k,131)
         mat(k,1264) = rxt(k,383)*y(k,131)
         mat(k,1226) = .900_r8*rxt(k,512)*y(k,131)
         mat(k,1109) = .800_r8*rxt(k,517)*y(k,131)
         mat(k,801) = rxt(k,487)*y(k,131)
         mat(k,1281) = rxt(k,453)*y(k,131) + rxt(k,454)*y(k,133)
         mat(k,818) = rxt(k,493)*y(k,131)
         mat(k,559) = rxt(k,496)*y(k,131)
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
         mat(k,2067) = -(rxt(k,218)*y(k,247) + rxt(k,219)*y(k,131) + rxt(k,220) &
                      *y(k,141) + rxt(k,221)*y(k,265) + rxt(k,229)*y(k,132) + rxt(k,315) &
                      *y(k,44) + rxt(k,347)*y(k,47) + rxt(k,366)*y(k,31) + rxt(k,373) &
                      *y(k,51) + rxt(k,387)*y(k,18) + rxt(k,395)*y(k,251) + rxt(k,406) &
                      *y(k,253) + rxt(k,429)*y(k,249) + rxt(k,435)*y(k,250) + rxt(k,438) &
                      *y(k,100) + rxt(k,443)*y(k,103) + rxt(k,454)*y(k,273) + rxt(k,499) &
                      *y(k,6) + rxt(k,502)*y(k,116) + rxt(k,508)*y(k,260) + rxt(k,519) &
                      *y(k,220) + rxt(k,522)*y(k,69))
         mat(k,2212) = -rxt(k,218)*y(k,133)
         mat(k,1791) = -rxt(k,219)*y(k,133)
         mat(k,2461) = -rxt(k,220)*y(k,133)
         mat(k,2007) = -rxt(k,221)*y(k,133)
         mat(k,1685) = -rxt(k,229)*y(k,133)
         mat(k,2487) = -rxt(k,315)*y(k,133)
         mat(k,1199) = -rxt(k,347)*y(k,133)
         mat(k,1188) = -rxt(k,366)*y(k,133)
         mat(k,1329) = -rxt(k,373)*y(k,133)
         mat(k,415) = -rxt(k,387)*y(k,133)
         mat(k,1407) = -rxt(k,395)*y(k,133)
         mat(k,1447) = -rxt(k,406)*y(k,133)
         mat(k,1352) = -rxt(k,429)*y(k,133)
         mat(k,1385) = -rxt(k,435)*y(k,133)
         mat(k,931) = -rxt(k,438)*y(k,133)
         mat(k,1306) = -rxt(k,443)*y(k,133)
         mat(k,1284) = -rxt(k,454)*y(k,133)
         mat(k,1076) = -rxt(k,499)*y(k,133)
         mat(k,1020) = -rxt(k,502)*y(k,133)
         mat(k,1157) = -rxt(k,508)*y(k,133)
         mat(k,1122) = -rxt(k,519)*y(k,133)
         mat(k,352) = -rxt(k,522)*y(k,133)
         mat(k,609) = rxt(k,280)*y(k,141)
         mat(k,2258) = rxt(k,247)*y(k,62)
         mat(k,1032) = rxt(k,247)*y(k,58) + rxt(k,249)*y(k,141) + rxt(k,250)*y(k,265)
         mat(k,972) = rxt(k,294)*y(k,91)
         mat(k,2302) = rxt(k,294)*y(k,75) + rxt(k,231)*y(k,265)
         mat(k,677) = .500_r8*rxt(k,411)*y(k,265)
         mat(k,1685) = mat(k,1685) + rxt(k,217)*y(k,141) + rxt(k,216)*y(k,143)
         mat(k,2461) = mat(k,2461) + rxt(k,280)*y(k,22) + rxt(k,249)*y(k,62) &
                      + rxt(k,217)*y(k,132)
         mat(k,2418) = rxt(k,216)*y(k,132)
         mat(k,585) = rxt(k,362)*y(k,265)
         mat(k,2007) = mat(k,2007) + rxt(k,250)*y(k,62) + rxt(k,231)*y(k,91) &
                      + .500_r8*rxt(k,411)*y(k,115) + rxt(k,362)*y(k,148)
         mat(k,936) = -(rxt(k,377)*y(k,265))
         mat(k,1958) = -rxt(k,377)*y(k,134)
         mat(k,1175) = rxt(k,366)*y(k,133)
         mat(k,597) = .500_r8*rxt(k,437)*y(k,265)
         mat(k,445) = rxt(k,444)*y(k,265)
         mat(k,494) = rxt(k,448)*y(k,265)
         mat(k,1131) = rxt(k,449)*y(k,265)
         mat(k,2026) = rxt(k,366)*y(k,31)
         mat(k,1958) = mat(k,1958) + .500_r8*rxt(k,437)*y(k,102) + rxt(k,444)*y(k,104) &
                      + rxt(k,448)*y(k,122) + rxt(k,449)*y(k,123)
         mat(k,431) = -(rxt(k,509)*y(k,265))
         mat(k,1907) = -rxt(k,509)*y(k,135)
         mat(k,2135) = rxt(k,506)*y(k,260)
         mat(k,1146) = rxt(k,506)*y(k,247)
         mat(k,2469) = -(rxt(k,187)*y(k,143) + 4._r8*rxt(k,188)*y(k,141) + rxt(k,189) &
                      *y(k,142) + rxt(k,190)*y(k,79) + rxt(k,191)*y(k,81) + rxt(k,196) &
                      *y(k,247) + rxt(k,202)*y(k,265) + (rxt(k,215) + rxt(k,217) &
                      ) * y(k,132) + rxt(k,220)*y(k,133) + rxt(k,225)*y(k,131) &
                      + rxt(k,249)*y(k,62) + rxt(k,251)*y(k,61) + rxt(k,254)*y(k,87) &
                      + rxt(k,257)*y(k,94) + rxt(k,280)*y(k,22) + rxt(k,281)*y(k,21) &
                      + rxt(k,283)*y(k,83) + rxt(k,285)*y(k,93) + rxt(k,316)*y(k,44) &
                      + rxt(k,524)*y(k,146) + (rxt(k,607) + rxt(k,608)) * y(k,257) &
                      + rxt(k,609)*y(k,259))
         mat(k,2426) = -rxt(k,187)*y(k,141)
         mat(k,1597) = -rxt(k,189)*y(k,141)
         mat(k,1515) = -rxt(k,190)*y(k,141)
         mat(k,649) = -rxt(k,191)*y(k,141)
         mat(k,2220) = -rxt(k,196)*y(k,141)
         mat(k,2015) = -rxt(k,202)*y(k,141)
         mat(k,1693) = -(rxt(k,215) + rxt(k,217)) * y(k,141)
         mat(k,2075) = -rxt(k,220)*y(k,141)
         mat(k,1799) = -rxt(k,225)*y(k,141)
         mat(k,1035) = -rxt(k,249)*y(k,141)
         mat(k,1648) = -rxt(k,251)*y(k,141)
         mat(k,1550) = -rxt(k,254)*y(k,141)
         mat(k,878) = -rxt(k,257)*y(k,141)
         mat(k,611) = -rxt(k,280)*y(k,141)
         mat(k,1621) = -rxt(k,281)*y(k,141)
         mat(k,886) = -rxt(k,283)*y(k,141)
         mat(k,847) = -rxt(k,285)*y(k,141)
         mat(k,2495) = -rxt(k,316)*y(k,141)
         mat(k,424) = -rxt(k,524)*y(k,141)
         mat(k,748) = -(rxt(k,607) + rxt(k,608)) * y(k,141)
         mat(k,554) = -rxt(k,609)*y(k,141)
         mat(k,2287) = rxt(k,194)*y(k,247)
         mat(k,1052) = rxt(k,210)*y(k,131) + rxt(k,211)*y(k,132) + rxt(k,214)*y(k,142) &
                      + rxt(k,612)*y(k,264)
         mat(k,1799) = mat(k,1799) + rxt(k,210)*y(k,119)
         mat(k,1693) = mat(k,1693) + rxt(k,211)*y(k,119)
         mat(k,1597) = mat(k,1597) + rxt(k,214)*y(k,119) + rxt(k,526)*y(k,159) &
                      + rxt(k,533)*y(k,161) + rxt(k,611)*y(k,259) + (rxt(k,176) &
                       +rxt(k,177))*y(k,261) + rxt(k,617)*y(k,266)
         mat(k,783) = rxt(k,526)*y(k,142)
         mat(k,1532) = rxt(k,533)*y(k,142)
         mat(k,916) = rxt(k,603)*y(k,258) + 1.150_r8*rxt(k,604)*y(k,264)
         mat(k,2220) = mat(k,2220) + rxt(k,194)*y(k,78)
         mat(k,895) = rxt(k,603)*y(k,243)
         mat(k,554) = mat(k,554) + rxt(k,611)*y(k,142)
         mat(k,1842) = (rxt(k,176)+rxt(k,177))*y(k,142)
         mat(k,903) = rxt(k,612)*y(k,119) + 1.150_r8*rxt(k,604)*y(k,243)
         mat(k,2015) = mat(k,2015) + 2.000_r8*rxt(k,204)*y(k,265)
         mat(k,856) = rxt(k,617)*y(k,142)
         mat(k,1586) = -(rxt(k,176)*y(k,261) + rxt(k,181)*y(k,262) + rxt(k,189) &
                      *y(k,141) + rxt(k,195)*y(k,78) + rxt(k,209)*y(k,256) + rxt(k,214) &
                      *y(k,119) + rxt(k,359)*y(k,245) + rxt(k,526)*y(k,159) + rxt(k,533) &
                      *y(k,161) + rxt(k,606)*y(k,257) + (rxt(k,610) + rxt(k,611) &
                      ) * y(k,259) + rxt(k,617)*y(k,266))
         mat(k,1827) = -rxt(k,176)*y(k,142)
         mat(k,224) = -rxt(k,181)*y(k,142)
         mat(k,2454) = -rxt(k,189)*y(k,142)
         mat(k,2272) = -rxt(k,195)*y(k,142)
         mat(k,577) = -rxt(k,209)*y(k,142)
         mat(k,1046) = -rxt(k,214)*y(k,142)
         mat(k,511) = -rxt(k,359)*y(k,142)
         mat(k,779) = -rxt(k,526)*y(k,142)
         mat(k,1522) = -rxt(k,533)*y(k,142)
         mat(k,745) = -rxt(k,606)*y(k,142)
         mat(k,553) = -(rxt(k,610) + rxt(k,611)) * y(k,142)
         mat(k,855) = -rxt(k,617)*y(k,142)
         mat(k,1556) = rxt(k,272)*y(k,143) + rxt(k,271)*y(k,247)
         mat(k,1608) = 2.000_r8*rxt(k,273)*y(k,21) + (rxt(k,275)+rxt(k,276))*y(k,61) &
                      + rxt(k,281)*y(k,141) + rxt(k,277)*y(k,247)
         mat(k,2251) = rxt(k,240)*y(k,143) + rxt(k,238)*y(k,247)
         mat(k,1634) = (rxt(k,275)+rxt(k,276))*y(k,21) + (2.000_r8*rxt(k,242) &
                       +2.000_r8*rxt(k,243))*y(k,61) + rxt(k,251)*y(k,141) &
                      + rxt(k,245)*y(k,247) + rxt(k,253)*y(k,265)
         mat(k,2272) = mat(k,2272) + rxt(k,198)*y(k,143) + rxt(k,192)*y(k,247)
         mat(k,522) = rxt(k,207)*y(k,265)
         mat(k,1046) = mat(k,1046) + rxt(k,213)*y(k,132)
         mat(k,303) = rxt(k,227)*y(k,261)
         mat(k,1784) = rxt(k,224)*y(k,143) + rxt(k,614)*y(k,264)
         mat(k,1678) = rxt(k,213)*y(k,119) + rxt(k,215)*y(k,141) + rxt(k,216)*y(k,143)
         mat(k,2060) = rxt(k,220)*y(k,141) + rxt(k,218)*y(k,247)
         mat(k,2454) = mat(k,2454) + rxt(k,281)*y(k,21) + rxt(k,251)*y(k,61) &
                      + rxt(k,215)*y(k,132) + rxt(k,220)*y(k,133) &
                      + 2.000_r8*rxt(k,188)*y(k,141) + 2.000_r8*rxt(k,187)*y(k,143) &
                      + rxt(k,196)*y(k,247) + rxt(k,180)*y(k,262) + rxt(k,202) &
                      *y(k,265)
         mat(k,1586) = mat(k,1586) + 2.000_r8*rxt(k,181)*y(k,262)
         mat(k,2411) = rxt(k,272)*y(k,19) + rxt(k,240)*y(k,58) + rxt(k,198)*y(k,78) &
                      + rxt(k,224)*y(k,131) + rxt(k,216)*y(k,132) &
                      + 2.000_r8*rxt(k,187)*y(k,141) + rxt(k,528)*y(k,159) &
                      + rxt(k,534)*y(k,161) + 2.000_r8*rxt(k,197)*y(k,247) &
                      + 2.000_r8*rxt(k,178)*y(k,261) + rxt(k,203)*y(k,265)
         mat(k,779) = mat(k,779) + rxt(k,528)*y(k,143)
         mat(k,1522) = mat(k,1522) + rxt(k,534)*y(k,143)
         mat(k,952) = rxt(k,341)*y(k,247)
         mat(k,982) = rxt(k,369)*y(k,247)
         mat(k,2347) = rxt(k,320)*y(k,247)
         mat(k,2205) = rxt(k,271)*y(k,19) + rxt(k,277)*y(k,21) + rxt(k,238)*y(k,58) &
                      + rxt(k,245)*y(k,61) + rxt(k,192)*y(k,78) + rxt(k,218)*y(k,133) &
                      + rxt(k,196)*y(k,141) + 2.000_r8*rxt(k,197)*y(k,143) &
                      + rxt(k,341)*y(k,237) + rxt(k,369)*y(k,238) + rxt(k,320) &
                      *y(k,241) + 2.000_r8*rxt(k,206)*y(k,247) + rxt(k,201)*y(k,265) &
                      + rxt(k,378)*y(k,268)
         mat(k,1827) = mat(k,1827) + rxt(k,227)*y(k,120) + 2.000_r8*rxt(k,178) &
                      *y(k,143)
         mat(k,224) = mat(k,224) + rxt(k,180)*y(k,141) + 2.000_r8*rxt(k,181)*y(k,142)
         mat(k,900) = rxt(k,614)*y(k,131)
         mat(k,2000) = rxt(k,253)*y(k,61) + rxt(k,207)*y(k,92) + rxt(k,202)*y(k,141) &
                      + rxt(k,203)*y(k,143) + rxt(k,201)*y(k,247)
         mat(k,865) = rxt(k,378)*y(k,247)
         mat(k,2425) = -(rxt(k,178)*y(k,261) + rxt(k,187)*y(k,141) + rxt(k,197) &
                      *y(k,247) + rxt(k,198)*y(k,78) + rxt(k,203)*y(k,265) + rxt(k,216) &
                      *y(k,132) + rxt(k,224)*y(k,131) + rxt(k,240)*y(k,58) + rxt(k,272) &
                      *y(k,19) + rxt(k,338)*y(k,27) + rxt(k,367)*y(k,31) + rxt(k,398) &
                      *y(k,111) + rxt(k,412)*y(k,118) + rxt(k,445)*y(k,100) + rxt(k,483) &
                      *y(k,150) + rxt(k,500)*y(k,6) + rxt(k,503)*y(k,116) + rxt(k,528) &
                      *y(k,159) + rxt(k,534)*y(k,161))
         mat(k,1841) = -rxt(k,178)*y(k,143)
         mat(k,2468) = -rxt(k,187)*y(k,143)
         mat(k,2219) = -rxt(k,197)*y(k,143)
         mat(k,2286) = -rxt(k,198)*y(k,143)
         mat(k,2014) = -rxt(k,203)*y(k,143)
         mat(k,1692) = -rxt(k,216)*y(k,143)
         mat(k,1798) = -rxt(k,224)*y(k,143)
         mat(k,2265) = -rxt(k,240)*y(k,143)
         mat(k,1564) = -rxt(k,272)*y(k,143)
         mat(k,626) = -rxt(k,338)*y(k,143)
         mat(k,1192) = -rxt(k,367)*y(k,143)
         mat(k,1321) = -rxt(k,398)*y(k,143)
         mat(k,1433) = -rxt(k,412)*y(k,143)
         mat(k,934) = -rxt(k,445)*y(k,143)
         mat(k,538) = -rxt(k,483)*y(k,143)
         mat(k,1078) = -rxt(k,500)*y(k,143)
         mat(k,1022) = -rxt(k,503)*y(k,143)
         mat(k,782) = -rxt(k,528)*y(k,143)
         mat(k,1531) = -rxt(k,534)*y(k,143)
         mat(k,2468) = mat(k,2468) + rxt(k,189)*y(k,142)
         mat(k,1596) = rxt(k,189)*y(k,141)
         mat(k,1483) = .150_r8*rxt(k,352)*y(k,247)
         mat(k,2219) = mat(k,2219) + .150_r8*rxt(k,352)*y(k,240) + .150_r8*rxt(k,403) &
                      *y(k,253)
         mat(k,1451) = .150_r8*rxt(k,403)*y(k,247)
         mat(k,362) = -(rxt(k,535)*y(k,161))
         mat(k,1517) = -rxt(k,535)*y(k,145)
         mat(k,1601) = rxt(k,274)*y(k,61)
         mat(k,1627) = rxt(k,274)*y(k,21) + 2.000_r8*rxt(k,244)*y(k,61)
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
         mat(k,417) = -(rxt(k,524)*y(k,141) + rxt(k,525)*y(k,265))
         mat(k,2431) = -rxt(k,524)*y(k,146)
         mat(k,1905) = -rxt(k,525)*y(k,146)
         mat(k,1234) = rxt(k,391)*y(k,265)
         mat(k,1719) = .100_r8*rxt(k,512)*y(k,270)
         mat(k,1883) = rxt(k,391)*y(k,95)
         mat(k,1215) = .100_r8*rxt(k,512)*y(k,131)
         mat(k,580) = -(rxt(k,362)*y(k,265))
         mat(k,1926) = -rxt(k,362)*y(k,148)
         mat(k,1657) = rxt(k,364)*y(k,240)
         mat(k,1455) = rxt(k,364)*y(k,132)
         mat(k,1651) = rxt(k,485)*y(k,233)
         mat(k,568) = rxt(k,485)*y(k,132)
         mat(k,535) = -(rxt(k,482)*y(k,132) + rxt(k,483)*y(k,143))
         mat(k,1655) = -rxt(k,482)*y(k,150)
         mat(k,2374) = -rxt(k,483)*y(k,150)
         mat(k,238) = .070_r8*rxt(k,469)*y(k,265)
         mat(k,1729) = rxt(k,467)*y(k,239)
         mat(k,220) = .060_r8*rxt(k,481)*y(k,265)
         mat(k,263) = .070_r8*rxt(k,497)*y(k,265)
         mat(k,691) = rxt(k,467)*y(k,131)
         mat(k,1921) = .070_r8*rxt(k,469)*y(k,68) + .060_r8*rxt(k,481)*y(k,151) &
                      + .070_r8*rxt(k,497)*y(k,229)
         mat(k,218) = -(rxt(k,481)*y(k,265))
         mat(k,1873) = -rxt(k,481)*y(k,151)
         mat(k,210) = .530_r8*rxt(k,458)*y(k,265)
         mat(k,1873) = mat(k,1873) + .530_r8*rxt(k,458)*y(k,8)
         mat(k,367) = -(rxt(k,484)*y(k,265))
         mat(k,1897) = -rxt(k,484)*y(k,152)
         mat(k,2129) = rxt(k,479)*y(k,267)
         mat(k,499) = rxt(k,479)*y(k,247)
         mat(k,588) = -(rxt(k,380)*y(k,265))
         mat(k,1927) = -rxt(k,380)*y(k,157)
         mat(k,2152) = rxt(k,378)*y(k,268)
         mat(k,861) = rxt(k,378)*y(k,247)
         mat(k,425) = -(rxt(k,384)*y(k,265))
         mat(k,1906) = -rxt(k,384)*y(k,158)
         mat(k,2134) = .850_r8*rxt(k,382)*y(k,269)
         mat(k,1258) = .850_r8*rxt(k,382)*y(k,247)
         mat(k,777) = -(rxt(k,526)*y(k,142) + rxt(k,528)*y(k,143) + rxt(k,531) &
                      *y(k,265))
         mat(k,1574) = -rxt(k,526)*y(k,159)
         mat(k,2378) = -rxt(k,528)*y(k,159)
         mat(k,1947) = -rxt(k,531)*y(k,159)
         mat(k,1520) = -(rxt(k,529)*y(k,21) + rxt(k,530)*y(k,61) + rxt(k,532)*y(k,132) &
                      + rxt(k,533)*y(k,142) + rxt(k,534)*y(k,143) + rxt(k,535) &
                      *y(k,145) + rxt(k,536)*y(k,265))
         mat(k,1605) = -rxt(k,529)*y(k,161)
         mat(k,1631) = -rxt(k,530)*y(k,161)
         mat(k,1675) = -rxt(k,532)*y(k,161)
         mat(k,1584) = -rxt(k,533)*y(k,161)
         mat(k,2409) = -rxt(k,534)*y(k,161)
         mat(k,364) = -rxt(k,535)*y(k,161)
         mat(k,1997) = -rxt(k,536)*y(k,161)
         mat(k,2451) = rxt(k,524)*y(k,146)
         mat(k,1584) = mat(k,1584) + rxt(k,526)*y(k,159)
         mat(k,2409) = mat(k,2409) + rxt(k,528)*y(k,159)
         mat(k,421) = rxt(k,524)*y(k,141)
         mat(k,778) = rxt(k,526)*y(k,142) + rxt(k,528)*y(k,143) + rxt(k,531)*y(k,265)
         mat(k,1997) = mat(k,1997) + rxt(k,531)*y(k,159)
         mat(k,991) = -(rxt(k,527)*y(k,265))
         mat(k,1964) = -rxt(k,527)*y(k,162)
         mat(k,1604) = rxt(k,529)*y(k,161)
         mat(k,1629) = rxt(k,530)*y(k,161)
         mat(k,350) = rxt(k,522)*y(k,133) + (rxt(k,523)+.500_r8*rxt(k,537))*y(k,265)
         mat(k,1665) = rxt(k,532)*y(k,161)
         mat(k,2028) = rxt(k,522)*y(k,69)
         mat(k,1579) = rxt(k,533)*y(k,161)
         mat(k,2382) = rxt(k,534)*y(k,161)
         mat(k,363) = rxt(k,535)*y(k,161)
         mat(k,419) = rxt(k,525)*y(k,265)
         mat(k,1519) = rxt(k,529)*y(k,21) + rxt(k,530)*y(k,61) + rxt(k,532)*y(k,132) &
                      + rxt(k,533)*y(k,142) + rxt(k,534)*y(k,143) + rxt(k,535) &
                      *y(k,145) + rxt(k,536)*y(k,265)
         mat(k,1964) = mat(k,1964) + (rxt(k,523)+.500_r8*rxt(k,537))*y(k,69) &
                      + rxt(k,525)*y(k,146) + rxt(k,536)*y(k,161)
         mat(k,307) = -(rxt(k,538)*y(k,276))
         mat(k,2499) = -rxt(k,538)*y(k,163)
         mat(k,990) = rxt(k,527)*y(k,265)
         mat(k,1889) = rxt(k,527)*y(k,162)
         mat(k,66) = .1056005_r8*rxt(k,567)*y(k,131) + .2381005_r8*rxt(k,566)*y(k,247)
         mat(k,1696) = .1056005_r8*rxt(k,567)*y(k,108)
         mat(k,117) = .5931005_r8*rxt(k,577)*y(k,265)
         mat(k,2103) = .2381005_r8*rxt(k,566)*y(k,108)
         mat(k,1845) = .5931005_r8*rxt(k,577)*y(k,214)
         mat(k,67) = .1026005_r8*rxt(k,567)*y(k,131) + .1308005_r8*rxt(k,566)*y(k,247)
         mat(k,1697) = .1026005_r8*rxt(k,567)*y(k,108)
         mat(k,118) = .1534005_r8*rxt(k,577)*y(k,265)
         mat(k,2104) = .1308005_r8*rxt(k,566)*y(k,108)
         mat(k,1846) = .1534005_r8*rxt(k,577)*y(k,214)
         mat(k,68) = .0521005_r8*rxt(k,567)*y(k,131) + .0348005_r8*rxt(k,566)*y(k,247)
         mat(k,1698) = .0521005_r8*rxt(k,567)*y(k,108)
         mat(k,119) = .0459005_r8*rxt(k,577)*y(k,265)
         mat(k,2105) = .0348005_r8*rxt(k,566)*y(k,108)
         mat(k,1847) = .0459005_r8*rxt(k,577)*y(k,214)
         mat(k,69) = .0143005_r8*rxt(k,567)*y(k,131) + .0076005_r8*rxt(k,566)*y(k,247)
         mat(k,1699) = .0143005_r8*rxt(k,567)*y(k,108)
         mat(k,120) = .0085005_r8*rxt(k,577)*y(k,265)
         mat(k,2106) = .0076005_r8*rxt(k,566)*y(k,108)
         mat(k,1848) = .0085005_r8*rxt(k,577)*y(k,214)
         mat(k,70) = .0166005_r8*rxt(k,567)*y(k,131) + .0113005_r8*rxt(k,566)*y(k,247)
         mat(k,1700) = .0166005_r8*rxt(k,567)*y(k,108)
         mat(k,121) = .0128005_r8*rxt(k,577)*y(k,265)
         mat(k,2107) = .0113005_r8*rxt(k,566)*y(k,108)
         mat(k,1849) = .0128005_r8*rxt(k,577)*y(k,214)
         mat(k,1053) = .2202005_r8*rxt(k,556)*y(k,143)
         mat(k,91) = .1279005_r8*rxt(k,555)*y(k,131) + .2202005_r8*rxt(k,554)*y(k,247)
         mat(k,79) = .0003005_r8*rxt(k,563)*y(k,131) + .0031005_r8*rxt(k,562)*y(k,247)
         mat(k,997) = .0508005_r8*rxt(k,575)*y(k,143)
         mat(k,85) = .0245005_r8*rxt(k,574)*y(k,131) + .0508005_r8*rxt(k,573)*y(k,247)
         mat(k,1702) = .1279005_r8*rxt(k,555)*y(k,7) + .0003005_r8*rxt(k,563)*y(k,105) &
                      + .0245005_r8*rxt(k,574)*y(k,117)
         mat(k,2365) = .2202005_r8*rxt(k,556)*y(k,6) + .0508005_r8*rxt(k,575)*y(k,116)
         mat(k,2109) = .2202005_r8*rxt(k,554)*y(k,7) + .0031005_r8*rxt(k,562)*y(k,105) &
                      + .0508005_r8*rxt(k,573)*y(k,117)
         mat(k,1054) = .2067005_r8*rxt(k,556)*y(k,143)
         mat(k,92) = .1792005_r8*rxt(k,555)*y(k,131) + .2067005_r8*rxt(k,554)*y(k,247)
         mat(k,80) = .0003005_r8*rxt(k,563)*y(k,131) + .0035005_r8*rxt(k,562)*y(k,247)
         mat(k,998) = .1149005_r8*rxt(k,575)*y(k,143)
         mat(k,86) = .0082005_r8*rxt(k,574)*y(k,131) + .1149005_r8*rxt(k,573)*y(k,247)
         mat(k,1703) = .1792005_r8*rxt(k,555)*y(k,7) + .0003005_r8*rxt(k,563)*y(k,105) &
                      + .0082005_r8*rxt(k,574)*y(k,117)
         mat(k,2366) = .2067005_r8*rxt(k,556)*y(k,6) + .1149005_r8*rxt(k,575)*y(k,116)
         mat(k,2110) = .2067005_r8*rxt(k,554)*y(k,7) + .0035005_r8*rxt(k,562)*y(k,105) &
                      + .1149005_r8*rxt(k,573)*y(k,117)
         mat(k,1055) = .0653005_r8*rxt(k,556)*y(k,143)
         mat(k,93) = .0676005_r8*rxt(k,555)*y(k,131) + .0653005_r8*rxt(k,554)*y(k,247)
         mat(k,81) = .0073005_r8*rxt(k,563)*y(k,131) + .0003005_r8*rxt(k,562)*y(k,247)
         mat(k,999) = .0348005_r8*rxt(k,575)*y(k,143)
         mat(k,87) = .0772005_r8*rxt(k,574)*y(k,131) + .0348005_r8*rxt(k,573)*y(k,247)
         mat(k,1704) = .0676005_r8*rxt(k,555)*y(k,7) + .0073005_r8*rxt(k,563)*y(k,105) &
                      + .0772005_r8*rxt(k,574)*y(k,117)
         mat(k,2367) = .0653005_r8*rxt(k,556)*y(k,6) + .0348005_r8*rxt(k,575)*y(k,116)
         mat(k,2111) = .0653005_r8*rxt(k,554)*y(k,7) + .0003005_r8*rxt(k,562)*y(k,105) &
                      + .0348005_r8*rxt(k,573)*y(k,117)
         mat(k,1056) = .1749305_r8*rxt(k,553)*y(k,133) + .1284005_r8*rxt(k,556) &
                      *y(k,143)
         mat(k,94) = .079_r8*rxt(k,555)*y(k,131) + .1284005_r8*rxt(k,554)*y(k,247)
         mat(k,917) = .0590245_r8*rxt(k,561)*y(k,133) + .0033005_r8*rxt(k,564) &
                      *y(k,143)
         mat(k,82) = .0057005_r8*rxt(k,563)*y(k,131) + .0271005_r8*rxt(k,562)*y(k,247)
         mat(k,1000) = .1749305_r8*rxt(k,572)*y(k,133) + .0554005_r8*rxt(k,575) &
                      *y(k,143)
         mat(k,88) = .0332005_r8*rxt(k,574)*y(k,131) + .0554005_r8*rxt(k,573)*y(k,247)
         mat(k,1705) = .079_r8*rxt(k,555)*y(k,7) + .0057005_r8*rxt(k,563)*y(k,105) &
                      + .0332005_r8*rxt(k,574)*y(k,117)
         mat(k,2018) = .1749305_r8*rxt(k,553)*y(k,6) + .0590245_r8*rxt(k,561)*y(k,100) &
                      + .1749305_r8*rxt(k,572)*y(k,116)
         mat(k,2368) = .1284005_r8*rxt(k,556)*y(k,6) + .0033005_r8*rxt(k,564)*y(k,100) &
                      + .0554005_r8*rxt(k,575)*y(k,116)
         mat(k,2112) = .1284005_r8*rxt(k,554)*y(k,7) + .0271005_r8*rxt(k,562)*y(k,105) &
                      + .0554005_r8*rxt(k,573)*y(k,117)
         mat(k,1057) = .5901905_r8*rxt(k,553)*y(k,133) + .114_r8*rxt(k,556)*y(k,143)
         mat(k,95) = .1254005_r8*rxt(k,555)*y(k,131) + .114_r8*rxt(k,554)*y(k,247)
         mat(k,918) = .0250245_r8*rxt(k,561)*y(k,133)
         mat(k,83) = .0623005_r8*rxt(k,563)*y(k,131) + .0474005_r8*rxt(k,562)*y(k,247)
         mat(k,1001) = .5901905_r8*rxt(k,572)*y(k,133) + .1278005_r8*rxt(k,575) &
                      *y(k,143)
         mat(k,89) = .130_r8*rxt(k,574)*y(k,131) + .1278005_r8*rxt(k,573)*y(k,247)
         mat(k,1706) = .1254005_r8*rxt(k,555)*y(k,7) + .0623005_r8*rxt(k,563)*y(k,105) &
                      + .130_r8*rxt(k,574)*y(k,117)
         mat(k,2019) = .5901905_r8*rxt(k,553)*y(k,6) + .0250245_r8*rxt(k,561)*y(k,100) &
                      + .5901905_r8*rxt(k,572)*y(k,116)
         mat(k,2369) = .114_r8*rxt(k,556)*y(k,6) + .1278005_r8*rxt(k,575)*y(k,116)
         mat(k,2113) = .114_r8*rxt(k,554)*y(k,7) + .0474005_r8*rxt(k,562)*y(k,105) &
                      + .1278005_r8*rxt(k,573)*y(k,117)
         mat(k,102) = .0097005_r8*rxt(k,560)*y(k,131) + .0023005_r8*rxt(k,559) &
                      *y(k,247)
         mat(k,108) = .1056005_r8*rxt(k,570)*y(k,131) + .2381005_r8*rxt(k,569) &
                      *y(k,247)
         mat(k,1710) = .0097005_r8*rxt(k,560)*y(k,9) + .1056005_r8*rxt(k,570)*y(k,110) &
                      + .0154005_r8*rxt(k,581)*y(k,224) + .0063005_r8*rxt(k,585) &
                      *y(k,228)
         mat(k,123) = .5931005_r8*rxt(k,578)*y(k,265)
         mat(k,129) = .0154005_r8*rxt(k,581)*y(k,131) + .1364005_r8*rxt(k,580) &
                      *y(k,247)
         mat(k,135) = .0063005_r8*rxt(k,585)*y(k,131) + .1677005_r8*rxt(k,584) &
                      *y(k,247)
         mat(k,2117) = .0023005_r8*rxt(k,559)*y(k,9) + .2381005_r8*rxt(k,569)*y(k,110) &
                      + .1364005_r8*rxt(k,580)*y(k,224) + .1677005_r8*rxt(k,584) &
                      *y(k,228)
         mat(k,1855) = .5931005_r8*rxt(k,578)*y(k,215)
         mat(k,103) = .0034005_r8*rxt(k,560)*y(k,131) + .0008005_r8*rxt(k,559) &
                      *y(k,247)
         mat(k,109) = .1026005_r8*rxt(k,570)*y(k,131) + .1308005_r8*rxt(k,569) &
                      *y(k,247)
         mat(k,1711) = .0034005_r8*rxt(k,560)*y(k,9) + .1026005_r8*rxt(k,570)*y(k,110) &
                      + .0452005_r8*rxt(k,581)*y(k,224) + .0237005_r8*rxt(k,585) &
                      *y(k,228)
         mat(k,124) = .1534005_r8*rxt(k,578)*y(k,265)
         mat(k,130) = .0452005_r8*rxt(k,581)*y(k,131) + .0101005_r8*rxt(k,580) &
                      *y(k,247)
         mat(k,136) = .0237005_r8*rxt(k,585)*y(k,131) + .0174005_r8*rxt(k,584) &
                      *y(k,247)
         mat(k,2118) = .0008005_r8*rxt(k,559)*y(k,9) + .1308005_r8*rxt(k,569)*y(k,110) &
                      + .0101005_r8*rxt(k,580)*y(k,224) + .0174005_r8*rxt(k,584) &
                      *y(k,228)
         mat(k,1856) = .1534005_r8*rxt(k,578)*y(k,215)
         mat(k,104) = .1579005_r8*rxt(k,560)*y(k,131) + .0843005_r8*rxt(k,559) &
                      *y(k,247)
         mat(k,110) = .0521005_r8*rxt(k,570)*y(k,131) + .0348005_r8*rxt(k,569) &
                      *y(k,247)
         mat(k,1712) = .1579005_r8*rxt(k,560)*y(k,9) + .0521005_r8*rxt(k,570)*y(k,110) &
                      + .0966005_r8*rxt(k,581)*y(k,224) + .0025005_r8*rxt(k,585) &
                      *y(k,228)
         mat(k,125) = .0459005_r8*rxt(k,578)*y(k,265)
         mat(k,131) = .0966005_r8*rxt(k,581)*y(k,131) + .0763005_r8*rxt(k,580) &
                      *y(k,247)
         mat(k,137) = .0025005_r8*rxt(k,585)*y(k,131) + .086_r8*rxt(k,584)*y(k,247)
         mat(k,2119) = .0843005_r8*rxt(k,559)*y(k,9) + .0348005_r8*rxt(k,569)*y(k,110) &
                      + .0763005_r8*rxt(k,580)*y(k,224) + .086_r8*rxt(k,584)*y(k,228)
         mat(k,1857) = .0459005_r8*rxt(k,578)*y(k,215)
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
         mat(k,105) = .0059005_r8*rxt(k,560)*y(k,131) + .0443005_r8*rxt(k,559) &
                      *y(k,247)
         mat(k,111) = .0143005_r8*rxt(k,570)*y(k,131) + .0076005_r8*rxt(k,569) &
                      *y(k,247)
         mat(k,1713) = .0059005_r8*rxt(k,560)*y(k,9) + .0143005_r8*rxt(k,570)*y(k,110) &
                      + .0073005_r8*rxt(k,581)*y(k,224) + .011_r8*rxt(k,585)*y(k,228)
         mat(k,126) = .0085005_r8*rxt(k,578)*y(k,265)
         mat(k,132) = .0073005_r8*rxt(k,581)*y(k,131) + .2157005_r8*rxt(k,580) &
                      *y(k,247)
         mat(k,138) = .011_r8*rxt(k,585)*y(k,131) + .0512005_r8*rxt(k,584)*y(k,247)
         mat(k,2120) = .0443005_r8*rxt(k,559)*y(k,9) + .0076005_r8*rxt(k,569)*y(k,110) &
                      + .2157005_r8*rxt(k,580)*y(k,224) + .0512005_r8*rxt(k,584) &
                      *y(k,228)
         mat(k,1858) = .0085005_r8*rxt(k,578)*y(k,215)
         mat(k,106) = .0536005_r8*rxt(k,560)*y(k,131) + .1621005_r8*rxt(k,559) &
                      *y(k,247)
         mat(k,112) = .0166005_r8*rxt(k,570)*y(k,131) + .0113005_r8*rxt(k,569) &
                      *y(k,247)
         mat(k,1714) = .0536005_r8*rxt(k,560)*y(k,9) + .0166005_r8*rxt(k,570)*y(k,110) &
                      + .238_r8*rxt(k,581)*y(k,224) + .1185005_r8*rxt(k,585)*y(k,228)
         mat(k,127) = .0128005_r8*rxt(k,578)*y(k,265)
         mat(k,133) = .238_r8*rxt(k,581)*y(k,131) + .0738005_r8*rxt(k,580)*y(k,247)
         mat(k,139) = .1185005_r8*rxt(k,585)*y(k,131) + .1598005_r8*rxt(k,584) &
                      *y(k,247)
         mat(k,2121) = .1621005_r8*rxt(k,559)*y(k,9) + .0113005_r8*rxt(k,569)*y(k,110) &
                      + .0738005_r8*rxt(k,580)*y(k,224) + .1598005_r8*rxt(k,584) &
                      *y(k,228)
         mat(k,1859) = .0128005_r8*rxt(k,578)*y(k,215)
         mat(k,122) = -(rxt(k,577)*y(k,265))
         mat(k,1863) = -rxt(k,577)*y(k,214)
         mat(k,128) = -(rxt(k,578)*y(k,265))
         mat(k,1864) = -rxt(k,578)*y(k,215)
         mat(k,231) = .100_r8*rxt(k,489)*y(k,265)
         mat(k,253) = .230_r8*rxt(k,491)*y(k,265)
         mat(k,1876) = .100_r8*rxt(k,489)*y(k,223) + .230_r8*rxt(k,491)*y(k,226)
         mat(k,728) = -(rxt(k,513)*y(k,265))
         mat(k,1943) = -rxt(k,513)*y(k,217)
         mat(k,2160) = rxt(k,511)*y(k,270)
         mat(k,1216) = rxt(k,511)*y(k,247)
         mat(k,684) = -(rxt(k,514)*y(k,265))
         mat(k,1938) = -rxt(k,514)*y(k,218)
         mat(k,1738) = .200_r8*rxt(k,507)*y(k,260) + .200_r8*rxt(k,517)*y(k,271)
         mat(k,2317) = .500_r8*rxt(k,505)*y(k,260)
         mat(k,1147) = .200_r8*rxt(k,507)*y(k,131) + .500_r8*rxt(k,505)*y(k,241)
         mat(k,1104) = .200_r8*rxt(k,517)*y(k,131)
         mat(k,528) = -(rxt(k,518)*y(k,265))
         mat(k,1920) = -rxt(k,518)*y(k,219)
         mat(k,2147) = rxt(k,516)*y(k,271)
         mat(k,1103) = rxt(k,516)*y(k,247)
         mat(k,1116) = -(rxt(k,519)*y(k,133) + rxt(k,520)*y(k,265))
         mat(k,2034) = -rxt(k,519)*y(k,220)
         mat(k,1972) = -rxt(k,520)*y(k,220)
         mat(k,1066) = .330_r8*rxt(k,500)*y(k,143)
         mat(k,1010) = .330_r8*rxt(k,503)*y(k,143)
         mat(k,1760) = .800_r8*rxt(k,507)*y(k,260) + .800_r8*rxt(k,517)*y(k,271)
         mat(k,2034) = mat(k,2034) + rxt(k,508)*y(k,260)
         mat(k,2388) = .330_r8*rxt(k,500)*y(k,6) + .330_r8*rxt(k,503)*y(k,116)
         mat(k,685) = rxt(k,514)*y(k,265)
         mat(k,2325) = .500_r8*rxt(k,505)*y(k,260) + rxt(k,515)*y(k,271)
         mat(k,1149) = .800_r8*rxt(k,507)*y(k,131) + rxt(k,508)*y(k,133) &
                      + .500_r8*rxt(k,505)*y(k,241)
         mat(k,1972) = mat(k,1972) + rxt(k,514)*y(k,218)
         mat(k,1107) = .800_r8*rxt(k,517)*y(k,131) + rxt(k,515)*y(k,241)
         mat(k,1163) = -(rxt(k,521)*y(k,265))
         mat(k,1976) = -rxt(k,521)*y(k,221)
         mat(k,1069) = .300_r8*rxt(k,500)*y(k,143)
         mat(k,1013) = .300_r8*rxt(k,503)*y(k,143)
         mat(k,1764) = .900_r8*rxt(k,512)*y(k,270)
         mat(k,2391) = .300_r8*rxt(k,500)*y(k,6) + .300_r8*rxt(k,503)*y(k,116)
         mat(k,2329) = rxt(k,510)*y(k,270)
         mat(k,1219) = .900_r8*rxt(k,512)*y(k,131) + rxt(k,510)*y(k,241)
         mat(k,662) = -(rxt(k,488)*y(k,265))
         mat(k,1935) = -rxt(k,488)*y(k,222)
         mat(k,2155) = rxt(k,486)*y(k,272)
         mat(k,792) = rxt(k,486)*y(k,247)
         mat(k,229) = -((rxt(k,489) + rxt(k,579)) * y(k,265))
         mat(k,1874) = -(rxt(k,489) + rxt(k,579)) * y(k,223)
         mat(k,134) = -(rxt(k,580)*y(k,247) + rxt(k,581)*y(k,131))
         mat(k,2124) = -rxt(k,580)*y(k,224)
         mat(k,1717) = -rxt(k,581)*y(k,224)
         mat(k,228) = rxt(k,579)*y(k,265)
         mat(k,1865) = rxt(k,579)*y(k,223)
         mat(k,245) = -(rxt(k,455)*y(k,265))
         mat(k,1877) = -rxt(k,455)*y(k,225)
         mat(k,2127) = rxt(k,452)*y(k,273)
         mat(k,1271) = rxt(k,452)*y(k,247)
         mat(k,254) = -(rxt(k,491)*y(k,265))
         mat(k,1879) = -rxt(k,491)*y(k,226)
         mat(k,766) = -(rxt(k,494)*y(k,265))
         mat(k,1946) = -rxt(k,494)*y(k,227)
         mat(k,2163) = rxt(k,492)*y(k,274)
         mat(k,809) = rxt(k,492)*y(k,247)
         mat(k,140) = -(rxt(k,584)*y(k,247) + rxt(k,585)*y(k,131))
         mat(k,2125) = -rxt(k,584)*y(k,228)
         mat(k,1718) = -rxt(k,585)*y(k,228)
         mat(k,252) = rxt(k,583)*y(k,265)
         mat(k,1866) = rxt(k,583)*y(k,226)
         mat(k,262) = -(rxt(k,497)*y(k,265))
         mat(k,1880) = -rxt(k,497)*y(k,229)
         mat(k,255) = .150_r8*rxt(k,491)*y(k,265)
         mat(k,1880) = mat(k,1880) + .150_r8*rxt(k,491)*y(k,226)
         mat(k,473) = -(rxt(k,498)*y(k,265))
         mat(k,1913) = -rxt(k,498)*y(k,230)
         mat(k,2140) = rxt(k,495)*y(k,275)
         mat(k,555) = rxt(k,495)*y(k,247)
         mat(k,569) = -(rxt(k,456)*y(k,247) + rxt(k,457)*y(k,131) + rxt(k,485) &
                      *y(k,132))
         mat(k,2151) = -rxt(k,456)*y(k,233)
         mat(k,1733) = -rxt(k,457)*y(k,233)
         mat(k,1656) = -rxt(k,485)*y(k,233)
         mat(k,280) = rxt(k,462)*y(k,265)
         mat(k,1925) = rxt(k,462)*y(k,24)
         mat(k,1086) = -(rxt(k,417)*y(k,247) + (rxt(k,418) + rxt(k,419)) * y(k,131))
         mat(k,2178) = -rxt(k,417)*y(k,234)
         mat(k,1757) = -(rxt(k,418) + rxt(k,419)) * y(k,234)
         mat(k,713) = rxt(k,420)*y(k,265)
         mat(k,271) = rxt(k,421)*y(k,265)
         mat(k,1969) = rxt(k,420)*y(k,2) + rxt(k,421)*y(k,17)
         mat(k,541) = -(rxt(k,459)*y(k,247) + rxt(k,460)*y(k,131))
         mat(k,2149) = -rxt(k,459)*y(k,235)
         mat(k,1730) = -rxt(k,460)*y(k,235)
         mat(k,211) = .350_r8*rxt(k,458)*y(k,265)
         mat(k,451) = rxt(k,461)*y(k,265)
         mat(k,1922) = .350_r8*rxt(k,458)*y(k,8) + rxt(k,461)*y(k,10)
         mat(k,481) = -(rxt(k,463)*y(k,247) + rxt(k,465)*y(k,131))
         mat(k,2141) = -rxt(k,463)*y(k,236)
         mat(k,1724) = -rxt(k,465)*y(k,236)
         mat(k,374) = rxt(k,464)*y(k,265)
         mat(k,232) = .070_r8*rxt(k,489)*y(k,265)
         mat(k,256) = .060_r8*rxt(k,491)*y(k,265)
         mat(k,1914) = rxt(k,464)*y(k,25) + .070_r8*rxt(k,489)*y(k,223) &
                      + .060_r8*rxt(k,491)*y(k,226)
         mat(k,950) = -(4._r8*rxt(k,339)*y(k,237) + rxt(k,340)*y(k,241) + rxt(k,341) &
                      *y(k,247) + rxt(k,342)*y(k,131))
         mat(k,2321) = -rxt(k,340)*y(k,237)
         mat(k,2175) = -rxt(k,341)*y(k,237)
         mat(k,1753) = -rxt(k,342)*y(k,237)
         mat(k,379) = .500_r8*rxt(k,344)*y(k,265)
         mat(k,335) = rxt(k,345)*y(k,58) + rxt(k,346)*y(k,265)
         mat(k,2236) = rxt(k,345)*y(k,30)
         mat(k,1960) = .500_r8*rxt(k,344)*y(k,29) + rxt(k,346)*y(k,30)
         mat(k,979) = -(rxt(k,368)*y(k,241) + rxt(k,369)*y(k,247) + rxt(k,370) &
                      *y(k,131))
         mat(k,2322) = -rxt(k,368)*y(k,238)
         mat(k,2177) = -rxt(k,369)*y(k,238)
         mat(k,1755) = -rxt(k,370)*y(k,238)
         mat(k,438) = rxt(k,371)*y(k,265)
         mat(k,341) = rxt(k,375)*y(k,58) + rxt(k,372)*y(k,265)
         mat(k,2237) = rxt(k,375)*y(k,33)
         mat(k,1963) = rxt(k,371)*y(k,32) + rxt(k,372)*y(k,33)
         mat(k,692) = -(rxt(k,466)*y(k,247) + rxt(k,467)*y(k,131))
         mat(k,2157) = -rxt(k,466)*y(k,239)
         mat(k,1739) = -rxt(k,467)*y(k,239)
         mat(k,317) = rxt(k,468)*y(k,265)
         mat(k,1739) = mat(k,1739) + rxt(k,457)*y(k,233)
         mat(k,2376) = rxt(k,483)*y(k,150)
         mat(k,536) = rxt(k,483)*y(k,143)
         mat(k,570) = rxt(k,457)*y(k,131) + .400_r8*rxt(k,456)*y(k,247)
         mat(k,2157) = mat(k,2157) + .400_r8*rxt(k,456)*y(k,233)
         mat(k,1939) = rxt(k,468)*y(k,34)
         mat(k,1472) = -(4._r8*rxt(k,350)*y(k,240) + rxt(k,351)*y(k,241) + rxt(k,352) &
                      *y(k,247) + rxt(k,353)*y(k,131) + rxt(k,364)*y(k,132) + rxt(k,392) &
                      *y(k,251) + rxt(k,425)*y(k,249) + rxt(k,430)*y(k,250) + rxt(k,439) &
                      *y(k,103) + rxt(k,450)*y(k,273))
         mat(k,2345) = -rxt(k,351)*y(k,240)
         mat(k,2200) = -rxt(k,352)*y(k,240)
         mat(k,1781) = -rxt(k,353)*y(k,240)
         mat(k,1673) = -rxt(k,364)*y(k,240)
         mat(k,1402) = -rxt(k,392)*y(k,240)
         mat(k,1347) = -rxt(k,425)*y(k,240)
         mat(k,1380) = -rxt(k,430)*y(k,240)
         mat(k,1301) = -rxt(k,439)*y(k,240)
         mat(k,1279) = -rxt(k,450)*y(k,240)
         mat(k,1073) = .060_r8*rxt(k,500)*y(k,143)
         mat(k,1197) = rxt(k,347)*y(k,133) + rxt(k,348)*y(k,265)
         mat(k,1326) = rxt(k,373)*y(k,133) + rxt(k,374)*y(k,265)
         mat(k,638) = .500_r8*rxt(k,355)*y(k,265)
         mat(k,929) = .080_r8*rxt(k,445)*y(k,143)
         mat(k,1317) = .100_r8*rxt(k,398)*y(k,143)
         mat(k,1017) = .060_r8*rxt(k,503)*y(k,143)
         mat(k,1422) = .280_r8*rxt(k,412)*y(k,143)
         mat(k,1781) = mat(k,1781) + .530_r8*rxt(k,396)*y(k,251) + rxt(k,405)*y(k,253) &
                      + rxt(k,408)*y(k,255) + rxt(k,383)*y(k,269)
         mat(k,2056) = rxt(k,347)*y(k,47) + rxt(k,373)*y(k,51) + .530_r8*rxt(k,395) &
                      *y(k,251) + rxt(k,406)*y(k,253)
         mat(k,2407) = .060_r8*rxt(k,500)*y(k,6) + .080_r8*rxt(k,445)*y(k,100) &
                      + .100_r8*rxt(k,398)*y(k,111) + .060_r8*rxt(k,503)*y(k,116) &
                      + .280_r8*rxt(k,412)*y(k,118)
         mat(k,1166) = .650_r8*rxt(k,521)*y(k,265)
         mat(k,1472) = mat(k,1472) + .530_r8*rxt(k,392)*y(k,251)
         mat(k,2345) = mat(k,2345) + .260_r8*rxt(k,393)*y(k,251) + rxt(k,402)*y(k,253) &
                      + .300_r8*rxt(k,381)*y(k,269)
         mat(k,2200) = mat(k,2200) + .450_r8*rxt(k,403)*y(k,253) + .200_r8*rxt(k,407) &
                      *y(k,255) + .150_r8*rxt(k,382)*y(k,269)
         mat(k,1402) = mat(k,1402) + .530_r8*rxt(k,396)*y(k,131) + .530_r8*rxt(k,395) &
                      *y(k,133) + .530_r8*rxt(k,392)*y(k,240) + .260_r8*rxt(k,393) &
                      *y(k,241)
         mat(k,1442) = rxt(k,405)*y(k,131) + rxt(k,406)*y(k,133) + rxt(k,402)*y(k,241) &
                      + .450_r8*rxt(k,403)*y(k,247) + 4.000_r8*rxt(k,404)*y(k,253)
         mat(k,723) = rxt(k,408)*y(k,131) + .200_r8*rxt(k,407)*y(k,247)
         mat(k,1994) = rxt(k,348)*y(k,47) + rxt(k,374)*y(k,51) + .500_r8*rxt(k,355) &
                      *y(k,53) + .650_r8*rxt(k,521)*y(k,221)
         mat(k,1263) = rxt(k,383)*y(k,131) + .300_r8*rxt(k,381)*y(k,241) &
                      + .150_r8*rxt(k,382)*y(k,247)
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
         mat(k,2360) = -(rxt(k,241)*y(k,61) + (4._r8*rxt(k,318) + 4._r8*rxt(k,319) &
                      ) * y(k,241) + rxt(k,320)*y(k,247) + rxt(k,321)*y(k,131) &
                      + rxt(k,340)*y(k,237) + rxt(k,351)*y(k,240) + rxt(k,368) &
                      *y(k,238) + rxt(k,381)*y(k,269) + rxt(k,393)*y(k,251) + rxt(k,402) &
                      *y(k,253) + rxt(k,426)*y(k,249) + rxt(k,431)*y(k,250) + rxt(k,440) &
                      *y(k,103) + rxt(k,451)*y(k,273) + rxt(k,505)*y(k,260) + rxt(k,510) &
                      *y(k,270) + rxt(k,515)*y(k,271))
         mat(k,1646) = -rxt(k,241)*y(k,241)
         mat(k,2218) = -rxt(k,320)*y(k,241)
         mat(k,1797) = -rxt(k,321)*y(k,241)
         mat(k,957) = -rxt(k,340)*y(k,241)
         mat(k,1482) = -rxt(k,351)*y(k,241)
         mat(k,987) = -rxt(k,368)*y(k,241)
         mat(k,1268) = -rxt(k,381)*y(k,241)
         mat(k,1410) = -rxt(k,393)*y(k,241)
         mat(k,1450) = -rxt(k,402)*y(k,241)
         mat(k,1355) = -rxt(k,426)*y(k,241)
         mat(k,1388) = -rxt(k,431)*y(k,241)
         mat(k,1309) = -rxt(k,440)*y(k,241)
         mat(k,1286) = -rxt(k,451)*y(k,241)
         mat(k,1160) = -rxt(k,505)*y(k,241)
         mat(k,1231) = -rxt(k,510)*y(k,241)
         mat(k,1114) = -rxt(k,515)*y(k,241)
         mat(k,1191) = .280_r8*rxt(k,367)*y(k,143)
         mat(k,752) = rxt(k,354)*y(k,265)
         mat(k,458) = .700_r8*rxt(k,323)*y(k,265)
         mat(k,2098) = rxt(k,235)*y(k,58) + rxt(k,291)*y(k,75) + rxt(k,330)*y(k,261) &
                      + rxt(k,324)*y(k,265)
         mat(k,2264) = rxt(k,235)*y(k,56)
         mat(k,976) = rxt(k,291)*y(k,56)
         mat(k,933) = .050_r8*rxt(k,445)*y(k,143)
         mat(k,1309) = mat(k,1309) + rxt(k,439)*y(k,240)
         mat(k,1797) = mat(k,1797) + rxt(k,353)*y(k,240) + .830_r8*rxt(k,471)*y(k,242) &
                      + .170_r8*rxt(k,477)*y(k,254)
         mat(k,2424) = .280_r8*rxt(k,367)*y(k,31) + .050_r8*rxt(k,445)*y(k,100)
         mat(k,1482) = mat(k,1482) + rxt(k,439)*y(k,103) + rxt(k,353)*y(k,131) &
                      + 4.000_r8*rxt(k,350)*y(k,240) + .900_r8*rxt(k,351)*y(k,241) &
                      + .450_r8*rxt(k,352)*y(k,247) + rxt(k,425)*y(k,249) + rxt(k,430) &
                      *y(k,250) + rxt(k,392)*y(k,251) + rxt(k,401)*y(k,253) &
                      + rxt(k,450)*y(k,273)
         mat(k,2360) = mat(k,2360) + .900_r8*rxt(k,351)*y(k,240)
         mat(k,829) = .830_r8*rxt(k,471)*y(k,131) + .330_r8*rxt(k,470)*y(k,247)
         mat(k,2218) = mat(k,2218) + .450_r8*rxt(k,352)*y(k,240) + .330_r8*rxt(k,470) &
                      *y(k,242) + .070_r8*rxt(k,476)*y(k,254)
         mat(k,1355) = mat(k,1355) + rxt(k,425)*y(k,240)
         mat(k,1388) = mat(k,1388) + rxt(k,430)*y(k,240)
         mat(k,1410) = mat(k,1410) + rxt(k,392)*y(k,240)
         mat(k,1450) = mat(k,1450) + rxt(k,401)*y(k,240)
         mat(k,967) = .170_r8*rxt(k,477)*y(k,131) + .070_r8*rxt(k,476)*y(k,247)
         mat(k,1840) = rxt(k,330)*y(k,56)
         mat(k,2013) = rxt(k,354)*y(k,52) + .700_r8*rxt(k,323)*y(k,55) + rxt(k,324) &
                      *y(k,56)
         mat(k,1286) = mat(k,1286) + rxt(k,450)*y(k,240)
         mat(k,822) = -(rxt(k,470)*y(k,247) + rxt(k,471)*y(k,131) + rxt(k,472) &
                      *y(k,132))
         mat(k,2167) = -rxt(k,470)*y(k,242)
         mat(k,1745) = -rxt(k,471)*y(k,242)
         mat(k,1662) = -rxt(k,472)*y(k,242)
         mat(k,909) = -(rxt(k,603)*y(k,258) + rxt(k,604)*y(k,264) + rxt(k,605) &
                      *y(k,257))
         mat(k,890) = -rxt(k,603)*y(k,243)
         mat(k,898) = -rxt(k,604)*y(k,243)
         mat(k,743) = -rxt(k,605)*y(k,243)
         mat(k,612) = -((rxt(k,389) + rxt(k,390)) * y(k,131))
         mat(k,1735) = -(rxt(k,389) + rxt(k,390)) * y(k,244)
         mat(k,410) = rxt(k,388)*y(k,265)
         mat(k,1929) = rxt(k,388)*y(k,18)
         mat(k,509) = -(rxt(k,359)*y(k,142))
         mat(k,1570) = -rxt(k,359)*y(k,245)
         mat(k,1728) = .750_r8*rxt(k,357)*y(k,246)
         mat(k,831) = .750_r8*rxt(k,357)*y(k,131)
         mat(k,832) = -(rxt(k,356)*y(k,247) + rxt(k,357)*y(k,131))
         mat(k,2168) = -rxt(k,356)*y(k,246)
         mat(k,1746) = -rxt(k,357)*y(k,246)
         mat(k,621) = rxt(k,363)*y(k,265)
         mat(k,1952) = rxt(k,363)*y(k,27)
         mat(k,2214) = -((rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,78) + rxt(k,196) &
                      *y(k,141) + rxt(k,197)*y(k,143) + rxt(k,201)*y(k,265) &
                      + 4._r8*rxt(k,206)*y(k,247) + rxt(k,218)*y(k,133) + rxt(k,223) &
                      *y(k,131) + rxt(k,228)*y(k,132) + (rxt(k,238) + rxt(k,239) &
                      ) * y(k,58) + rxt(k,245)*y(k,61) + rxt(k,271)*y(k,19) + rxt(k,277) &
                      *y(k,21) + rxt(k,314)*y(k,44) + rxt(k,320)*y(k,241) + rxt(k,327) &
                      *y(k,248) + rxt(k,341)*y(k,237) + rxt(k,352)*y(k,240) + rxt(k,356) &
                      *y(k,246) + rxt(k,369)*y(k,238) + rxt(k,378)*y(k,268) + rxt(k,382) &
                      *y(k,269) + rxt(k,394)*y(k,251) + rxt(k,403)*y(k,253) + rxt(k,407) &
                      *y(k,255) + rxt(k,417)*y(k,234) + rxt(k,427)*y(k,249) + rxt(k,432) &
                      *y(k,250) + rxt(k,441)*y(k,103) + rxt(k,452)*y(k,273) + rxt(k,456) &
                      *y(k,233) + rxt(k,459)*y(k,235) + rxt(k,463)*y(k,236) + rxt(k,466) &
                      *y(k,239) + rxt(k,470)*y(k,242) + rxt(k,473)*y(k,252) + rxt(k,476) &
                      *y(k,254) + rxt(k,479)*y(k,267) + rxt(k,486)*y(k,272) + rxt(k,492) &
                      *y(k,274) + rxt(k,495)*y(k,275) + rxt(k,506)*y(k,260) + rxt(k,511) &
                      *y(k,270) + rxt(k,516)*y(k,271))
         mat(k,2281) = -(rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,247)
         mat(k,2463) = -rxt(k,196)*y(k,247)
         mat(k,2420) = -rxt(k,197)*y(k,247)
         mat(k,2009) = -rxt(k,201)*y(k,247)
         mat(k,2069) = -rxt(k,218)*y(k,247)
         mat(k,1793) = -rxt(k,223)*y(k,247)
         mat(k,1687) = -rxt(k,228)*y(k,247)
         mat(k,2260) = -(rxt(k,238) + rxt(k,239)) * y(k,247)
         mat(k,1642) = -rxt(k,245)*y(k,247)
         mat(k,1562) = -rxt(k,271)*y(k,247)
         mat(k,1616) = -rxt(k,277)*y(k,247)
         mat(k,2489) = -rxt(k,314)*y(k,247)
         mat(k,2356) = -rxt(k,320)*y(k,247)
         mat(k,491) = -rxt(k,327)*y(k,247)
         mat(k,956) = -rxt(k,341)*y(k,247)
         mat(k,1480) = -rxt(k,352)*y(k,247)
         mat(k,838) = -rxt(k,356)*y(k,247)
         mat(k,986) = -rxt(k,369)*y(k,247)
         mat(k,869) = -rxt(k,378)*y(k,247)
         mat(k,1267) = -rxt(k,382)*y(k,247)
         mat(k,1408) = -rxt(k,394)*y(k,247)
         mat(k,1448) = -rxt(k,403)*y(k,247)
         mat(k,727) = -rxt(k,407)*y(k,247)
         mat(k,1094) = -rxt(k,417)*y(k,247)
         mat(k,1353) = -rxt(k,427)*y(k,247)
         mat(k,1386) = -rxt(k,432)*y(k,247)
         mat(k,1307) = -rxt(k,441)*y(k,247)
         mat(k,1285) = -rxt(k,452)*y(k,247)
         mat(k,574) = -rxt(k,456)*y(k,247)
         mat(k,547) = -rxt(k,459)*y(k,247)
         mat(k,486) = -rxt(k,463)*y(k,247)
         mat(k,696) = -rxt(k,466)*y(k,247)
         mat(k,828) = -rxt(k,470)*y(k,247)
         mat(k,789) = -rxt(k,473)*y(k,247)
         mat(k,966) = -rxt(k,476)*y(k,247)
         mat(k,505) = -rxt(k,479)*y(k,247)
         mat(k,804) = -rxt(k,486)*y(k,247)
         mat(k,821) = -rxt(k,492)*y(k,247)
         mat(k,562) = -rxt(k,495)*y(k,247)
         mat(k,1158) = -rxt(k,506)*y(k,247)
         mat(k,1229) = -rxt(k,511)*y(k,247)
         mat(k,1112) = -rxt(k,516)*y(k,247)
         mat(k,1077) = .570_r8*rxt(k,500)*y(k,143)
         mat(k,213) = .650_r8*rxt(k,458)*y(k,265)
         mat(k,1562) = mat(k,1562) + rxt(k,270)*y(k,44)
         mat(k,1616) = mat(k,1616) + rxt(k,282)*y(k,265)
         mat(k,330) = .350_r8*rxt(k,336)*y(k,265)
         mat(k,625) = .130_r8*rxt(k,338)*y(k,143)
         mat(k,314) = rxt(k,343)*y(k,265)
         mat(k,1190) = .280_r8*rxt(k,367)*y(k,143)
         mat(k,2489) = mat(k,2489) + rxt(k,270)*y(k,19) + rxt(k,234)*y(k,58) &
                      + rxt(k,315)*y(k,133) + rxt(k,316)*y(k,141)
         mat(k,633) = rxt(k,299)*y(k,58) + rxt(k,300)*y(k,265)
         mat(k,405) = rxt(k,302)*y(k,58) + rxt(k,303)*y(k,265)
         mat(k,149) = rxt(k,349)*y(k,265)
         mat(k,859) = rxt(k,322)*y(k,265)
         mat(k,2094) = rxt(k,331)*y(k,261)
         mat(k,2260) = mat(k,2260) + rxt(k,234)*y(k,44) + rxt(k,299)*y(k,45) &
                      + rxt(k,302)*y(k,48) + rxt(k,237)*y(k,81)
         mat(k,1642) = mat(k,1642) + rxt(k,241)*y(k,241) + rxt(k,252)*y(k,265)
         mat(k,1207) = rxt(k,334)*y(k,265)
         mat(k,240) = .730_r8*rxt(k,469)*y(k,265)
         mat(k,353) = .500_r8*rxt(k,537)*y(k,265)
         mat(k,1213) = rxt(k,360)*y(k,265)
         mat(k,1102) = rxt(k,361)*y(k,265)
         mat(k,2281) = mat(k,2281) + rxt(k,195)*y(k,142)
         mat(k,647) = rxt(k,237)*y(k,58) + rxt(k,191)*y(k,141) + rxt(k,200)*y(k,265)
         mat(k,251) = rxt(k,325)*y(k,265)
         mat(k,945) = rxt(k,326)*y(k,265)
         mat(k,1247) = rxt(k,391)*y(k,265)
         mat(k,1256) = rxt(k,376)*y(k,265)
         mat(k,932) = .370_r8*rxt(k,445)*y(k,143)
         mat(k,657) = .300_r8*rxt(k,436)*y(k,265)
         mat(k,602) = rxt(k,437)*y(k,265)
         mat(k,1307) = mat(k,1307) + rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133) &
                      + rxt(k,439)*y(k,240) + 1.200_r8*rxt(k,440)*y(k,241)
         mat(k,448) = rxt(k,444)*y(k,265)
         mat(k,1320) = .140_r8*rxt(k,398)*y(k,143)
         mat(k,397) = .200_r8*rxt(k,400)*y(k,265)
         mat(k,678) = .500_r8*rxt(k,411)*y(k,265)
         mat(k,1021) = .570_r8*rxt(k,503)*y(k,143)
         mat(k,1430) = .280_r8*rxt(k,412)*y(k,143)
         mat(k,497) = rxt(k,448)*y(k,265)
         mat(k,1142) = rxt(k,449)*y(k,265)
         mat(k,1793) = mat(k,1793) + rxt(k,442)*y(k,103) + rxt(k,418)*y(k,234) &
                      + rxt(k,460)*y(k,235) + rxt(k,465)*y(k,236) + rxt(k,342) &
                      *y(k,237) + rxt(k,370)*y(k,238) + rxt(k,321)*y(k,241) &
                      + .170_r8*rxt(k,471)*y(k,242) + rxt(k,389)*y(k,244) &
                      + .250_r8*rxt(k,357)*y(k,246) + rxt(k,329)*y(k,248) &
                      + .920_r8*rxt(k,428)*y(k,249) + .920_r8*rxt(k,434)*y(k,250) &
                      + .470_r8*rxt(k,396)*y(k,251) + .400_r8*rxt(k,474)*y(k,252) &
                      + .830_r8*rxt(k,477)*y(k,254) + rxt(k,480)*y(k,267) + rxt(k,379) &
                      *y(k,268) + .900_r8*rxt(k,512)*y(k,270) + .800_r8*rxt(k,517) &
                      *y(k,271) + rxt(k,487)*y(k,272) + rxt(k,453)*y(k,273) &
                      + rxt(k,493)*y(k,274) + rxt(k,496)*y(k,275)
         mat(k,2069) = mat(k,2069) + rxt(k,315)*y(k,44) + rxt(k,443)*y(k,103) &
                      + rxt(k,429)*y(k,249) + rxt(k,435)*y(k,250) + .470_r8*rxt(k,395) &
                      *y(k,251) + rxt(k,221)*y(k,265) + rxt(k,454)*y(k,273)
         mat(k,2463) = mat(k,2463) + rxt(k,316)*y(k,44) + rxt(k,191)*y(k,81)
         mat(k,1593) = rxt(k,195)*y(k,78) + rxt(k,359)*y(k,245)
         mat(k,2420) = mat(k,2420) + .570_r8*rxt(k,500)*y(k,6) + .130_r8*rxt(k,338) &
                      *y(k,27) + .280_r8*rxt(k,367)*y(k,31) + .370_r8*rxt(k,445) &
                      *y(k,100) + .140_r8*rxt(k,398)*y(k,111) + .570_r8*rxt(k,503) &
                      *y(k,116) + .280_r8*rxt(k,412)*y(k,118) + rxt(k,203)*y(k,265)
         mat(k,222) = .800_r8*rxt(k,481)*y(k,265)
         mat(k,994) = rxt(k,527)*y(k,265)
         mat(k,1169) = .200_r8*rxt(k,521)*y(k,265)
         mat(k,235) = .280_r8*rxt(k,489)*y(k,265)
         mat(k,261) = .380_r8*rxt(k,491)*y(k,265)
         mat(k,266) = .630_r8*rxt(k,497)*y(k,265)
         mat(k,1094) = mat(k,1094) + rxt(k,418)*y(k,131)
         mat(k,547) = mat(k,547) + rxt(k,460)*y(k,131)
         mat(k,486) = mat(k,486) + rxt(k,465)*y(k,131)
         mat(k,956) = mat(k,956) + rxt(k,342)*y(k,131) + 2.400_r8*rxt(k,339)*y(k,237) &
                      + rxt(k,340)*y(k,241)
         mat(k,986) = mat(k,986) + rxt(k,370)*y(k,131) + rxt(k,368)*y(k,241)
         mat(k,1480) = mat(k,1480) + rxt(k,439)*y(k,103) + .900_r8*rxt(k,351)*y(k,241) &
                      + rxt(k,425)*y(k,249) + rxt(k,430)*y(k,250) + .470_r8*rxt(k,392) &
                      *y(k,251) + rxt(k,450)*y(k,273)
         mat(k,2356) = mat(k,2356) + rxt(k,241)*y(k,61) + 1.200_r8*rxt(k,440)*y(k,103) &
                      + rxt(k,321)*y(k,131) + rxt(k,340)*y(k,237) + rxt(k,368) &
                      *y(k,238) + .900_r8*rxt(k,351)*y(k,240) + 4.000_r8*rxt(k,318) &
                      *y(k,241) + rxt(k,426)*y(k,249) + rxt(k,431)*y(k,250) &
                      + .730_r8*rxt(k,393)*y(k,251) + rxt(k,402)*y(k,253) &
                      + .500_r8*rxt(k,505)*y(k,260) + .300_r8*rxt(k,381)*y(k,269) &
                      + rxt(k,510)*y(k,270) + rxt(k,515)*y(k,271) + .800_r8*rxt(k,451) &
                      *y(k,273)
         mat(k,828) = mat(k,828) + .170_r8*rxt(k,471)*y(k,131) + .070_r8*rxt(k,470) &
                      *y(k,247)
         mat(k,618) = rxt(k,389)*y(k,131)
         mat(k,512) = rxt(k,359)*y(k,142)
         mat(k,838) = mat(k,838) + .250_r8*rxt(k,357)*y(k,131)
         mat(k,2214) = mat(k,2214) + .070_r8*rxt(k,470)*y(k,242) + .160_r8*rxt(k,473) &
                      *y(k,252) + .330_r8*rxt(k,476)*y(k,254)
         mat(k,491) = mat(k,491) + rxt(k,329)*y(k,131)
         mat(k,1353) = mat(k,1353) + .920_r8*rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133) &
                      + rxt(k,425)*y(k,240) + rxt(k,426)*y(k,241)
         mat(k,1386) = mat(k,1386) + .920_r8*rxt(k,434)*y(k,131) + rxt(k,435)*y(k,133) &
                      + rxt(k,430)*y(k,240) + rxt(k,431)*y(k,241)
         mat(k,1408) = mat(k,1408) + .470_r8*rxt(k,396)*y(k,131) + .470_r8*rxt(k,395) &
                      *y(k,133) + .470_r8*rxt(k,392)*y(k,240) + .730_r8*rxt(k,393) &
                      *y(k,241)
         mat(k,789) = mat(k,789) + .400_r8*rxt(k,474)*y(k,131) + .160_r8*rxt(k,473) &
                      *y(k,247)
         mat(k,1448) = mat(k,1448) + rxt(k,402)*y(k,241)
         mat(k,966) = mat(k,966) + .830_r8*rxt(k,477)*y(k,131) + .330_r8*rxt(k,476) &
                      *y(k,247)
         mat(k,1158) = mat(k,1158) + .500_r8*rxt(k,505)*y(k,241)
         mat(k,1836) = rxt(k,331)*y(k,56)
         mat(k,2009) = mat(k,2009) + .650_r8*rxt(k,458)*y(k,8) + rxt(k,282)*y(k,21) &
                      + .350_r8*rxt(k,336)*y(k,26) + rxt(k,343)*y(k,28) + rxt(k,300) &
                      *y(k,45) + rxt(k,303)*y(k,48) + rxt(k,349)*y(k,49) + rxt(k,322) &
                      *y(k,54) + rxt(k,252)*y(k,61) + rxt(k,334)*y(k,64) &
                      + .730_r8*rxt(k,469)*y(k,68) + .500_r8*rxt(k,537)*y(k,69) &
                      + rxt(k,360)*y(k,76) + rxt(k,361)*y(k,77) + rxt(k,200)*y(k,81) &
                      + rxt(k,325)*y(k,88) + rxt(k,326)*y(k,89) + rxt(k,391)*y(k,95) &
                      + rxt(k,376)*y(k,97) + .300_r8*rxt(k,436)*y(k,101) + rxt(k,437) &
                      *y(k,102) + rxt(k,444)*y(k,104) + .200_r8*rxt(k,400)*y(k,112) &
                      + .500_r8*rxt(k,411)*y(k,115) + rxt(k,448)*y(k,122) + rxt(k,449) &
                      *y(k,123) + rxt(k,221)*y(k,133) + rxt(k,203)*y(k,143) &
                      + .800_r8*rxt(k,481)*y(k,151) + rxt(k,527)*y(k,162) &
                      + .200_r8*rxt(k,521)*y(k,221) + .280_r8*rxt(k,489)*y(k,223) &
                      + .380_r8*rxt(k,491)*y(k,226) + .630_r8*rxt(k,497)*y(k,229)
         mat(k,505) = mat(k,505) + rxt(k,480)*y(k,131)
         mat(k,869) = mat(k,869) + rxt(k,379)*y(k,131)
         mat(k,1267) = mat(k,1267) + .300_r8*rxt(k,381)*y(k,241)
         mat(k,1229) = mat(k,1229) + .900_r8*rxt(k,512)*y(k,131) + rxt(k,510)*y(k,241)
         mat(k,1112) = mat(k,1112) + .800_r8*rxt(k,517)*y(k,131) + rxt(k,515)*y(k,241)
         mat(k,804) = mat(k,804) + rxt(k,487)*y(k,131)
         mat(k,1285) = mat(k,1285) + rxt(k,453)*y(k,131) + rxt(k,454)*y(k,133) &
                      + rxt(k,450)*y(k,240) + .800_r8*rxt(k,451)*y(k,241)
         mat(k,821) = mat(k,821) + rxt(k,493)*y(k,131)
         mat(k,562) = mat(k,562) + rxt(k,496)*y(k,131)
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
         mat(k,487) = -(rxt(k,327)*y(k,247) + rxt(k,329)*y(k,131))
         mat(k,2142) = -rxt(k,327)*y(k,248)
         mat(k,1725) = -rxt(k,329)*y(k,248)
         mat(k,2472) = rxt(k,314)*y(k,247)
         mat(k,2142) = mat(k,2142) + rxt(k,314)*y(k,44)
         mat(k,1343) = -(rxt(k,425)*y(k,240) + rxt(k,426)*y(k,241) + rxt(k,427) &
                      *y(k,247) + rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133))
         mat(k,1467) = -rxt(k,425)*y(k,249)
         mat(k,2340) = -rxt(k,426)*y(k,249)
         mat(k,2195) = -rxt(k,427)*y(k,249)
         mat(k,1776) = -rxt(k,428)*y(k,249)
         mat(k,2051) = -rxt(k,429)*y(k,249)
         mat(k,926) = .600_r8*rxt(k,446)*y(k,265)
         mat(k,1989) = .600_r8*rxt(k,446)*y(k,100)
         mat(k,1376) = -(rxt(k,430)*y(k,240) + rxt(k,431)*y(k,241) + rxt(k,432) &
                      *y(k,247) + rxt(k,434)*y(k,131) + rxt(k,435)*y(k,133))
         mat(k,1468) = -rxt(k,430)*y(k,250)
         mat(k,2341) = -rxt(k,431)*y(k,250)
         mat(k,2196) = -rxt(k,432)*y(k,250)
         mat(k,1777) = -rxt(k,434)*y(k,250)
         mat(k,2052) = -rxt(k,435)*y(k,250)
         mat(k,927) = .400_r8*rxt(k,446)*y(k,265)
         mat(k,1990) = .400_r8*rxt(k,446)*y(k,100)
         mat(k,1400) = -(rxt(k,392)*y(k,240) + rxt(k,393)*y(k,241) + rxt(k,394) &
                      *y(k,247) + rxt(k,395)*y(k,133) + (rxt(k,396) + rxt(k,397) &
                      ) * y(k,131))
         mat(k,1469) = -rxt(k,392)*y(k,251)
         mat(k,2342) = -rxt(k,393)*y(k,251)
         mat(k,2197) = -rxt(k,394)*y(k,251)
         mat(k,2053) = -rxt(k,395)*y(k,251)
         mat(k,1778) = -(rxt(k,396) + rxt(k,397)) * y(k,251)
         mat(k,1315) = .500_r8*rxt(k,399)*y(k,265)
         mat(k,394) = .200_r8*rxt(k,400)*y(k,265)
         mat(k,1419) = rxt(k,413)*y(k,265)
         mat(k,1991) = .500_r8*rxt(k,399)*y(k,111) + .200_r8*rxt(k,400)*y(k,112) &
                      + rxt(k,413)*y(k,118)
         mat(k,784) = -(rxt(k,473)*y(k,247) + rxt(k,474)*y(k,131) + rxt(k,475) &
                      *y(k,132))
         mat(k,2164) = -rxt(k,473)*y(k,252)
         mat(k,1742) = -rxt(k,474)*y(k,252)
         mat(k,1661) = -rxt(k,475)*y(k,252)
         mat(k,1441) = -(rxt(k,401)*y(k,240) + rxt(k,402)*y(k,241) + rxt(k,403) &
                      *y(k,247) + 4._r8*rxt(k,404)*y(k,253) + rxt(k,405)*y(k,131) &
                      + rxt(k,406)*y(k,133) + rxt(k,414)*y(k,132))
         mat(k,1471) = -rxt(k,401)*y(k,253)
         mat(k,2344) = -rxt(k,402)*y(k,253)
         mat(k,2199) = -rxt(k,403)*y(k,253)
         mat(k,1780) = -rxt(k,405)*y(k,253)
         mat(k,2055) = -rxt(k,406)*y(k,253)
         mat(k,1672) = -rxt(k,414)*y(k,253)
         mat(k,1316) = .500_r8*rxt(k,399)*y(k,265)
         mat(k,395) = .500_r8*rxt(k,400)*y(k,265)
         mat(k,1993) = .500_r8*rxt(k,399)*y(k,111) + .500_r8*rxt(k,400)*y(k,112)
         mat(k,959) = -(rxt(k,476)*y(k,247) + rxt(k,477)*y(k,131) + rxt(k,478) &
                      *y(k,132))
         mat(k,2176) = -rxt(k,476)*y(k,254)
         mat(k,1754) = -rxt(k,477)*y(k,254)
         mat(k,1664) = -rxt(k,478)*y(k,254)
         mat(k,721) = -(rxt(k,407)*y(k,247) + rxt(k,408)*y(k,131))
         mat(k,2159) = -rxt(k,407)*y(k,255)
         mat(k,1741) = -rxt(k,408)*y(k,255)
         mat(k,564) = rxt(k,409)*y(k,265)
         mat(k,389) = rxt(k,410)*y(k,265)
         mat(k,1942) = rxt(k,409)*y(k,113) + rxt(k,410)*y(k,114)
         mat(k,575) = -(rxt(k,208)*y(k,141) + rxt(k,209)*y(k,142))
         mat(k,2433) = -rxt(k,208)*y(k,256)
         mat(k,1572) = -rxt(k,209)*y(k,256)
         mat(k,2433) = mat(k,2433) + rxt(k,607)*y(k,257)
         mat(k,904) = .900_r8*rxt(k,605)*y(k,257) + .800_r8*rxt(k,603)*y(k,258)
         mat(k,738) = rxt(k,607)*y(k,141) + .900_r8*rxt(k,605)*y(k,243)
         mat(k,888) = .800_r8*rxt(k,603)*y(k,243)
         mat(k,739) = -(rxt(k,605)*y(k,243) + rxt(k,606)*y(k,142) + (rxt(k,607) &
                      + rxt(k,608)) * y(k,141))
         mat(k,905) = -rxt(k,605)*y(k,257)
         mat(k,1573) = -rxt(k,606)*y(k,257)
         mat(k,2436) = -(rxt(k,607) + rxt(k,608)) * y(k,257)
         mat(k,889) = -(rxt(k,603)*y(k,243))
         mat(k,907) = -rxt(k,603)*y(k,258)
         mat(k,1040) = rxt(k,612)*y(k,264)
         mat(k,1748) = rxt(k,614)*y(k,264)
         mat(k,2442) = rxt(k,607)*y(k,257)
         mat(k,1576) = rxt(k,611)*y(k,259)
         mat(k,741) = rxt(k,607)*y(k,141)
         mat(k,550) = rxt(k,611)*y(k,142)
         mat(k,896) = rxt(k,612)*y(k,119) + rxt(k,614)*y(k,131)
         mat(k,548) = -(rxt(k,609)*y(k,141) + (rxt(k,610) + rxt(k,611)) * y(k,142))
         mat(k,2432) = -rxt(k,609)*y(k,259)
         mat(k,1571) = -(rxt(k,610) + rxt(k,611)) * y(k,259)
         mat(k,1150) = -(rxt(k,505)*y(k,241) + rxt(k,506)*y(k,247) + rxt(k,507) &
                      *y(k,131) + rxt(k,508)*y(k,133))
         mat(k,2328) = -rxt(k,505)*y(k,260)
         mat(k,2183) = -rxt(k,506)*y(k,260)
         mat(k,1763) = -rxt(k,507)*y(k,260)
         mat(k,2037) = -rxt(k,508)*y(k,260)
         mat(k,1068) = rxt(k,499)*y(k,133)
         mat(k,1012) = rxt(k,502)*y(k,133)
         mat(k,2037) = mat(k,2037) + rxt(k,499)*y(k,6) + rxt(k,502)*y(k,116) &
                      + .500_r8*rxt(k,519)*y(k,220)
         mat(k,433) = rxt(k,509)*y(k,265)
         mat(k,1117) = .500_r8*rxt(k,519)*y(k,133)
         mat(k,1975) = rxt(k,509)*y(k,135)
         mat(k,1832) = -(rxt(k,173)*y(k,79) + rxt(k,174)*y(k,276) + (rxt(k,176) &
                      + rxt(k,177)) * y(k,142) + rxt(k,178)*y(k,143) + (rxt(k,226) &
                      + rxt(k,227)) * y(k,120) + rxt(k,259)*y(k,35) + rxt(k,260) &
                      *y(k,36) + rxt(k,261)*y(k,38) + rxt(k,262)*y(k,39) + rxt(k,263) &
                      *y(k,40) + rxt(k,264)*y(k,41) + rxt(k,265)*y(k,42) + (rxt(k,266) &
                      + rxt(k,267)) * y(k,87) + rxt(k,286)*y(k,37) + rxt(k,287) &
                      *y(k,57) + rxt(k,288)*y(k,80) + (rxt(k,289) + rxt(k,290) &
                      ) * y(k,83) + rxt(k,295)*y(k,66) + rxt(k,296)*y(k,67) + rxt(k,309) &
                      *y(k,43) + rxt(k,310)*y(k,45) + rxt(k,311)*y(k,84) + rxt(k,312) &
                      *y(k,85) + rxt(k,313)*y(k,86) + (rxt(k,330) + rxt(k,331) &
                      + rxt(k,332)) * y(k,56) + rxt(k,333)*y(k,88))
         mat(k,1507) = -rxt(k,173)*y(k,261)
         mat(k,2512) = -rxt(k,174)*y(k,261)
         mat(k,1591) = -(rxt(k,176) + rxt(k,177)) * y(k,261)
         mat(k,2416) = -rxt(k,178)*y(k,261)
         mat(k,305) = -(rxt(k,226) + rxt(k,227)) * y(k,261)
         mat(k,145) = -rxt(k,259)*y(k,261)
         mat(k,188) = -rxt(k,260)*y(k,261)
         mat(k,155) = -rxt(k,261)*y(k,261)
         mat(k,198) = -rxt(k,262)*y(k,261)
         mat(k,159) = -rxt(k,263)*y(k,261)
         mat(k,203) = -rxt(k,264)*y(k,261)
         mat(k,163) = -rxt(k,265)*y(k,261)
         mat(k,1544) = -(rxt(k,266) + rxt(k,267)) * y(k,261)
         mat(k,194) = -rxt(k,286)*y(k,261)
         mat(k,464) = -rxt(k,287)*y(k,261)
         mat(k,172) = -rxt(k,288)*y(k,261)
         mat(k,883) = -(rxt(k,289) + rxt(k,290)) * y(k,261)
         mat(k,269) = -rxt(k,295)*y(k,261)
         mat(k,294) = -rxt(k,296)*y(k,261)
         mat(k,517) = -rxt(k,309)*y(k,261)
         mat(k,631) = -rxt(k,310)*y(k,261)
         mat(k,289) = -rxt(k,311)*y(k,261)
         mat(k,299) = -rxt(k,312)*y(k,261)
         mat(k,358) = -rxt(k,313)*y(k,261)
         mat(k,2090) = -(rxt(k,330) + rxt(k,331) + rxt(k,332)) * y(k,261)
         mat(k,249) = -rxt(k,333)*y(k,261)
         mat(k,1591) = mat(k,1591) + rxt(k,209)*y(k,256)
         mat(k,915) = .850_r8*rxt(k,604)*y(k,264)
         mat(k,579) = rxt(k,209)*y(k,142)
         mat(k,902) = .850_r8*rxt(k,604)*y(k,243)
         mat(k,223) = -(rxt(k,180)*y(k,141) + rxt(k,181)*y(k,142))
         mat(k,2429) = -rxt(k,180)*y(k,262)
         mat(k,1568) = -rxt(k,181)*y(k,262)
         mat(k,1486) = rxt(k,182)*y(k,263)
         mat(k,2429) = mat(k,2429) + rxt(k,184)*y(k,263)
         mat(k,1568) = mat(k,1568) + rxt(k,185)*y(k,263)
         mat(k,2370) = rxt(k,186)*y(k,263)
         mat(k,225) = rxt(k,182)*y(k,65) + rxt(k,184)*y(k,141) + rxt(k,185)*y(k,142) &
                      + rxt(k,186)*y(k,143)
         mat(k,226) = -(rxt(k,182)*y(k,65) + rxt(k,184)*y(k,141) + rxt(k,185)*y(k,142) &
                      + rxt(k,186)*y(k,143))
         mat(k,1487) = -rxt(k,182)*y(k,263)
         mat(k,2430) = -rxt(k,184)*y(k,263)
         mat(k,1569) = -rxt(k,185)*y(k,263)
         mat(k,2371) = -rxt(k,186)*y(k,263)
         mat(k,1569) = mat(k,1569) + rxt(k,176)*y(k,261)
         mat(k,1811) = rxt(k,176)*y(k,142)
         mat(k,897) = -(rxt(k,604)*y(k,243) + rxt(k,612)*y(k,119) + rxt(k,614) &
                      *y(k,131))
         mat(k,908) = -rxt(k,604)*y(k,264)
         mat(k,1041) = -rxt(k,612)*y(k,264)
         mat(k,1749) = -rxt(k,614)*y(k,264)
         mat(k,1490) = rxt(k,615)*y(k,266)
         mat(k,1577) = rxt(k,606)*y(k,257) + rxt(k,610)*y(k,259) + rxt(k,617)*y(k,266)
         mat(k,742) = rxt(k,606)*y(k,142)
         mat(k,551) = rxt(k,610)*y(k,142)
         mat(k,851) = rxt(k,615)*y(k,65) + rxt(k,617)*y(k,142)
         mat(k,2006) = -(rxt(k,199)*y(k,79) + rxt(k,200)*y(k,81) + rxt(k,201)*y(k,247) &
                      + rxt(k,202)*y(k,141) + rxt(k,203)*y(k,143) + (4._r8*rxt(k,204) &
                      + 4._r8*rxt(k,205)) * y(k,265) + rxt(k,207)*y(k,92) + rxt(k,221) &
                      *y(k,133) + rxt(k,222)*y(k,119) + rxt(k,230)*y(k,132) + rxt(k,231) &
                      *y(k,91) + rxt(k,250)*y(k,62) + (rxt(k,252) + rxt(k,253) &
                      ) * y(k,61) + rxt(k,255)*y(k,87) + rxt(k,258)*y(k,94) + rxt(k,282) &
                      *y(k,21) + rxt(k,284)*y(k,83) + rxt(k,298)*y(k,43) + rxt(k,300) &
                      *y(k,45) + rxt(k,301)*y(k,46) + rxt(k,303)*y(k,48) + rxt(k,305) &
                      *y(k,57) + rxt(k,306)*y(k,84) + rxt(k,307)*y(k,85) + rxt(k,308) &
                      *y(k,86) + rxt(k,317)*y(k,44) + rxt(k,322)*y(k,54) + rxt(k,323) &
                      *y(k,55) + rxt(k,324)*y(k,56) + rxt(k,325)*y(k,88) + rxt(k,326) &
                      *y(k,89) + rxt(k,334)*y(k,64) + rxt(k,336)*y(k,26) + rxt(k,343) &
                      *y(k,28) + rxt(k,344)*y(k,29) + rxt(k,346)*y(k,30) + rxt(k,348) &
                      *y(k,47) + rxt(k,349)*y(k,49) + rxt(k,354)*y(k,52) + rxt(k,355) &
                      *y(k,53) + rxt(k,360)*y(k,76) + rxt(k,361)*y(k,77) + rxt(k,362) &
                      *y(k,148) + rxt(k,363)*y(k,27) + rxt(k,371)*y(k,32) + rxt(k,372) &
                      *y(k,33) + rxt(k,374)*y(k,51) + rxt(k,376)*y(k,97) + rxt(k,377) &
                      *y(k,134) + rxt(k,380)*y(k,157) + rxt(k,384)*y(k,158) + rxt(k,385) &
                      *y(k,31) + rxt(k,386)*y(k,50) + rxt(k,388)*y(k,18) + rxt(k,391) &
                      *y(k,95) + rxt(k,399)*y(k,111) + rxt(k,400)*y(k,112) + rxt(k,409) &
                      *y(k,113) + rxt(k,410)*y(k,114) + rxt(k,411)*y(k,115) + rxt(k,413) &
                      *y(k,118) + rxt(k,416)*y(k,1) + rxt(k,420)*y(k,2) + rxt(k,421) &
                      *y(k,17) + rxt(k,422)*y(k,96) + rxt(k,423)*y(k,98) + rxt(k,424) &
                      *y(k,99) + rxt(k,436)*y(k,101) + rxt(k,437)*y(k,102) + rxt(k,444) &
                      *y(k,104) + rxt(k,446)*y(k,100) + rxt(k,447)*y(k,106) + rxt(k,448) &
                      *y(k,122) + rxt(k,449)*y(k,123) + rxt(k,455)*y(k,225) + rxt(k,458) &
                      *y(k,8) + rxt(k,461)*y(k,10) + rxt(k,462)*y(k,24) + rxt(k,464) &
                      *y(k,25) + rxt(k,468)*y(k,34) + rxt(k,469)*y(k,68) + rxt(k,481) &
                      *y(k,151) + rxt(k,484)*y(k,152) + rxt(k,488)*y(k,222) + (rxt(k,489) &
                      + rxt(k,579)) * y(k,223) + rxt(k,491)*y(k,226) + rxt(k,494) &
                      *y(k,227) + rxt(k,497)*y(k,229) + rxt(k,498)*y(k,230) + rxt(k,501) &
                      *y(k,6) + rxt(k,504)*y(k,116) + rxt(k,509)*y(k,135) + rxt(k,513) &
                      *y(k,217) + rxt(k,514)*y(k,218) + rxt(k,518)*y(k,219) + rxt(k,520) &
                      *y(k,220) + rxt(k,521)*y(k,221) + (rxt(k,523) + rxt(k,537) &
                      ) * y(k,69) + rxt(k,525)*y(k,146) + rxt(k,527)*y(k,162) &
                      + rxt(k,531)*y(k,159) + rxt(k,536)*y(k,161) + rxt(k,539) &
                      *y(k,127))
         mat(k,1508) = -rxt(k,199)*y(k,265)
         mat(k,646) = -rxt(k,200)*y(k,265)
         mat(k,2211) = -rxt(k,201)*y(k,265)
         mat(k,2460) = -rxt(k,202)*y(k,265)
         mat(k,2417) = -rxt(k,203)*y(k,265)
         mat(k,524) = -rxt(k,207)*y(k,265)
         mat(k,2066) = -rxt(k,221)*y(k,265)
         mat(k,1050) = -rxt(k,222)*y(k,265)
         mat(k,1684) = -rxt(k,230)*y(k,265)
         mat(k,2301) = -rxt(k,231)*y(k,265)
         mat(k,1031) = -rxt(k,250)*y(k,265)
         mat(k,1640) = -(rxt(k,252) + rxt(k,253)) * y(k,265)
         mat(k,1545) = -rxt(k,255)*y(k,265)
         mat(k,876) = -rxt(k,258)*y(k,265)
         mat(k,1614) = -rxt(k,282)*y(k,265)
         mat(k,884) = -rxt(k,284)*y(k,265)
         mat(k,518) = -rxt(k,298)*y(k,265)
         mat(k,632) = -rxt(k,300)*y(k,265)
         mat(k,166) = -rxt(k,301)*y(k,265)
         mat(k,404) = -rxt(k,303)*y(k,265)
         mat(k,465) = -rxt(k,305)*y(k,265)
         mat(k,290) = -rxt(k,306)*y(k,265)
         mat(k,300) = -rxt(k,307)*y(k,265)
         mat(k,359) = -rxt(k,308)*y(k,265)
         mat(k,2486) = -rxt(k,317)*y(k,265)
         mat(k,858) = -rxt(k,322)*y(k,265)
         mat(k,456) = -rxt(k,323)*y(k,265)
         mat(k,2091) = -rxt(k,324)*y(k,265)
         mat(k,250) = -rxt(k,325)*y(k,265)
         mat(k,944) = -rxt(k,326)*y(k,265)
         mat(k,1206) = -rxt(k,334)*y(k,265)
         mat(k,329) = -rxt(k,336)*y(k,265)
         mat(k,313) = -rxt(k,343)*y(k,265)
         mat(k,381) = -rxt(k,344)*y(k,265)
         mat(k,337) = -rxt(k,346)*y(k,265)
         mat(k,1198) = -rxt(k,348)*y(k,265)
         mat(k,148) = -rxt(k,349)*y(k,265)
         mat(k,751) = -rxt(k,354)*y(k,265)
         mat(k,640) = -rxt(k,355)*y(k,265)
         mat(k,1212) = -rxt(k,360)*y(k,265)
         mat(k,1101) = -rxt(k,361)*y(k,265)
         mat(k,584) = -rxt(k,362)*y(k,265)
         mat(k,624) = -rxt(k,363)*y(k,265)
         mat(k,440) = -rxt(k,371)*y(k,265)
         mat(k,343) = -rxt(k,372)*y(k,265)
         mat(k,1328) = -rxt(k,374)*y(k,265)
         mat(k,1255) = -rxt(k,376)*y(k,265)
         mat(k,940) = -rxt(k,377)*y(k,265)
         mat(k,592) = -rxt(k,380)*y(k,265)
         mat(k,428) = -rxt(k,384)*y(k,265)
         mat(k,1187) = -rxt(k,385)*y(k,265)
         mat(k,1127) = -rxt(k,386)*y(k,265)
         mat(k,414) = -rxt(k,388)*y(k,265)
         mat(k,1245) = -rxt(k,391)*y(k,265)
         mat(k,1319) = -rxt(k,399)*y(k,265)
         mat(k,396) = -rxt(k,400)*y(k,265)
         mat(k,567) = -rxt(k,409)*y(k,265)
         mat(k,392) = -rxt(k,410)*y(k,265)
         mat(k,676) = -rxt(k,411)*y(k,265)
         mat(k,1427) = -rxt(k,413)*y(k,265)
         mat(k,706) = -rxt(k,416)*y(k,265)
         mat(k,717) = -rxt(k,420)*y(k,265)
         mat(k,272) = -rxt(k,421)*y(k,265)
         mat(k,285) = -rxt(k,422)*y(k,265)
         mat(k,400) = -rxt(k,423)*y(k,265)
         mat(k,175) = -rxt(k,424)*y(k,265)
         mat(k,656) = -rxt(k,436)*y(k,265)
         mat(k,601) = -rxt(k,437)*y(k,265)
         mat(k,447) = -rxt(k,444)*y(k,265)
         mat(k,930) = -rxt(k,446)*y(k,265)
         mat(k,759) = -rxt(k,447)*y(k,265)
         mat(k,496) = -rxt(k,448)*y(k,265)
         mat(k,1141) = -rxt(k,449)*y(k,265)
         mat(k,247) = -rxt(k,455)*y(k,265)
         mat(k,212) = -rxt(k,458)*y(k,265)
         mat(k,453) = -rxt(k,461)*y(k,265)
         mat(k,281) = -rxt(k,462)*y(k,265)
         mat(k,376) = -rxt(k,464)*y(k,265)
         mat(k,318) = -rxt(k,468)*y(k,265)
         mat(k,239) = -rxt(k,469)*y(k,265)
         mat(k,221) = -rxt(k,481)*y(k,265)
         mat(k,370) = -rxt(k,484)*y(k,265)
         mat(k,669) = -rxt(k,488)*y(k,265)
         mat(k,234) = -(rxt(k,489) + rxt(k,579)) * y(k,265)
         mat(k,260) = -rxt(k,491)*y(k,265)
         mat(k,775) = -rxt(k,494)*y(k,265)
         mat(k,265) = -rxt(k,497)*y(k,265)
         mat(k,477) = -rxt(k,498)*y(k,265)
         mat(k,1075) = -rxt(k,501)*y(k,265)
         mat(k,1019) = -rxt(k,504)*y(k,265)
         mat(k,435) = -rxt(k,509)*y(k,265)
         mat(k,735) = -rxt(k,513)*y(k,265)
         mat(k,687) = -rxt(k,514)*y(k,265)
         mat(k,532) = -rxt(k,518)*y(k,265)
         mat(k,1121) = -rxt(k,520)*y(k,265)
         mat(k,1168) = -rxt(k,521)*y(k,265)
         mat(k,351) = -(rxt(k,523) + rxt(k,537)) * y(k,265)
         mat(k,422) = -rxt(k,525)*y(k,265)
         mat(k,993) = -rxt(k,527)*y(k,265)
         mat(k,780) = -rxt(k,531)*y(k,265)
         mat(k,1527) = -rxt(k,536)*y(k,265)
         mat(k,142) = -rxt(k,539)*y(k,265)
         mat(k,1075) = mat(k,1075) + .630_r8*rxt(k,500)*y(k,143)
         mat(k,329) = mat(k,329) + .650_r8*rxt(k,336)*y(k,265)
         mat(k,624) = mat(k,624) + .130_r8*rxt(k,338)*y(k,143)
         mat(k,381) = mat(k,381) + .500_r8*rxt(k,344)*y(k,265)
         mat(k,1187) = mat(k,1187) + .360_r8*rxt(k,367)*y(k,143)
         mat(k,2486) = mat(k,2486) + rxt(k,316)*y(k,141)
         mat(k,456) = mat(k,456) + .300_r8*rxt(k,323)*y(k,265)
         mat(k,2091) = mat(k,2091) + rxt(k,330)*y(k,261)
         mat(k,2257) = rxt(k,239)*y(k,247)
         mat(k,971) = rxt(k,293)*y(k,276)
         mat(k,2278) = rxt(k,198)*y(k,143) + 2.000_r8*rxt(k,193)*y(k,247)
         mat(k,1508) = mat(k,1508) + rxt(k,190)*y(k,141) + rxt(k,173)*y(k,261)
         mat(k,646) = mat(k,646) + rxt(k,191)*y(k,141)
         mat(k,884) = mat(k,884) + rxt(k,283)*y(k,141) + rxt(k,289)*y(k,261)
         mat(k,1545) = mat(k,1545) + rxt(k,254)*y(k,141) + rxt(k,266)*y(k,261)
         mat(k,250) = mat(k,250) + rxt(k,333)*y(k,261)
         mat(k,845) = rxt(k,285)*y(k,141)
         mat(k,876) = mat(k,876) + rxt(k,257)*y(k,141)
         mat(k,930) = mat(k,930) + .320_r8*rxt(k,445)*y(k,143)
         mat(k,759) = mat(k,759) + .600_r8*rxt(k,447)*y(k,265)
         mat(k,1319) = mat(k,1319) + .240_r8*rxt(k,398)*y(k,143)
         mat(k,396) = mat(k,396) + .100_r8*rxt(k,400)*y(k,265)
         mat(k,1019) = mat(k,1019) + .630_r8*rxt(k,503)*y(k,143)
         mat(k,1427) = mat(k,1427) + .360_r8*rxt(k,412)*y(k,143)
         mat(k,1790) = rxt(k,223)*y(k,247)
         mat(k,2066) = mat(k,2066) + rxt(k,218)*y(k,247)
         mat(k,2460) = mat(k,2460) + rxt(k,316)*y(k,44) + rxt(k,190)*y(k,79) &
                      + rxt(k,191)*y(k,81) + rxt(k,283)*y(k,83) + rxt(k,254)*y(k,87) &
                      + rxt(k,285)*y(k,93) + rxt(k,257)*y(k,94) + rxt(k,196)*y(k,247)
         mat(k,2417) = mat(k,2417) + .630_r8*rxt(k,500)*y(k,6) + .130_r8*rxt(k,338) &
                      *y(k,27) + .360_r8*rxt(k,367)*y(k,31) + rxt(k,198)*y(k,78) &
                      + .320_r8*rxt(k,445)*y(k,100) + .240_r8*rxt(k,398)*y(k,111) &
                      + .630_r8*rxt(k,503)*y(k,116) + .360_r8*rxt(k,412)*y(k,118) &
                      + rxt(k,197)*y(k,247)
         mat(k,592) = mat(k,592) + .500_r8*rxt(k,380)*y(k,265)
         mat(k,247) = mat(k,247) + .500_r8*rxt(k,455)*y(k,265)
         mat(k,573) = .400_r8*rxt(k,456)*y(k,247)
         mat(k,1477) = .450_r8*rxt(k,352)*y(k,247)
         mat(k,827) = .400_r8*rxt(k,470)*y(k,247)
         mat(k,2211) = mat(k,2211) + rxt(k,239)*y(k,58) + 2.000_r8*rxt(k,193)*y(k,78) &
                      + rxt(k,223)*y(k,131) + rxt(k,218)*y(k,133) + rxt(k,196) &
                      *y(k,141) + rxt(k,197)*y(k,143) + .400_r8*rxt(k,456)*y(k,233) &
                      + .450_r8*rxt(k,352)*y(k,240) + .400_r8*rxt(k,470)*y(k,242) &
                      + .450_r8*rxt(k,403)*y(k,253) + .400_r8*rxt(k,476)*y(k,254) &
                      + .200_r8*rxt(k,407)*y(k,255) + .150_r8*rxt(k,382)*y(k,269)
         mat(k,1446) = .450_r8*rxt(k,403)*y(k,247)
         mat(k,965) = .400_r8*rxt(k,476)*y(k,247)
         mat(k,726) = .200_r8*rxt(k,407)*y(k,247)
         mat(k,1833) = rxt(k,330)*y(k,56) + rxt(k,173)*y(k,79) + rxt(k,289)*y(k,83) &
                      + rxt(k,266)*y(k,87) + rxt(k,333)*y(k,88) + 2.000_r8*rxt(k,174) &
                      *y(k,276)
         mat(k,2006) = mat(k,2006) + .650_r8*rxt(k,336)*y(k,26) + .500_r8*rxt(k,344) &
                      *y(k,29) + .300_r8*rxt(k,323)*y(k,55) + .600_r8*rxt(k,447) &
                      *y(k,106) + .100_r8*rxt(k,400)*y(k,112) + .500_r8*rxt(k,380) &
                      *y(k,157) + .500_r8*rxt(k,455)*y(k,225)
         mat(k,1266) = .150_r8*rxt(k,382)*y(k,247)
         mat(k,2513) = rxt(k,293)*y(k,75) + 2.000_r8*rxt(k,174)*y(k,261)
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
         mat(k,849) = -(rxt(k,615)*y(k,65) + rxt(k,617)*y(k,142))
         mat(k,1488) = -rxt(k,615)*y(k,266)
         mat(k,1575) = -rxt(k,617)*y(k,266)
         mat(k,2439) = rxt(k,608)*y(k,257) + rxt(k,609)*y(k,259)
         mat(k,740) = rxt(k,608)*y(k,141)
         mat(k,549) = rxt(k,609)*y(k,141)
         mat(k,500) = -(rxt(k,479)*y(k,247) + rxt(k,480)*y(k,131))
         mat(k,2143) = -rxt(k,479)*y(k,267)
         mat(k,1726) = -rxt(k,480)*y(k,267)
         mat(k,237) = .200_r8*rxt(k,469)*y(k,265)
         mat(k,219) = .140_r8*rxt(k,481)*y(k,265)
         mat(k,368) = rxt(k,484)*y(k,265)
         mat(k,1916) = .200_r8*rxt(k,469)*y(k,68) + .140_r8*rxt(k,481)*y(k,151) &
                      + rxt(k,484)*y(k,152)
         mat(k,862) = -(rxt(k,378)*y(k,247) + rxt(k,379)*y(k,131))
         mat(k,2170) = -rxt(k,378)*y(k,268)
         mat(k,1747) = -rxt(k,379)*y(k,268)
         mat(k,1174) = rxt(k,385)*y(k,265)
         mat(k,589) = .500_r8*rxt(k,380)*y(k,265)
         mat(k,1954) = rxt(k,385)*y(k,31) + .500_r8*rxt(k,380)*y(k,157)
         mat(k,1261) = -(rxt(k,381)*y(k,241) + rxt(k,382)*y(k,247) + rxt(k,383) &
                      *y(k,131))
         mat(k,2335) = -rxt(k,381)*y(k,269)
         mat(k,2190) = -rxt(k,382)*y(k,269)
         mat(k,1771) = -rxt(k,383)*y(k,269)
         mat(k,1071) = .060_r8*rxt(k,500)*y(k,143)
         mat(k,1125) = rxt(k,386)*y(k,265)
         mat(k,1015) = .060_r8*rxt(k,503)*y(k,143)
         mat(k,2398) = .060_r8*rxt(k,500)*y(k,6) + .060_r8*rxt(k,503)*y(k,116)
         mat(k,426) = rxt(k,384)*y(k,265)
         mat(k,1165) = .150_r8*rxt(k,521)*y(k,265)
         mat(k,1984) = rxt(k,386)*y(k,50) + rxt(k,384)*y(k,158) + .150_r8*rxt(k,521) &
                      *y(k,221)
         mat(k,1222) = -(rxt(k,510)*y(k,241) + rxt(k,511)*y(k,247) + rxt(k,512) &
                      *y(k,131))
         mat(k,2333) = -rxt(k,510)*y(k,270)
         mat(k,2188) = -rxt(k,511)*y(k,270)
         mat(k,1768) = -rxt(k,512)*y(k,270)
         mat(k,2043) = .500_r8*rxt(k,519)*y(k,220)
         mat(k,733) = rxt(k,513)*y(k,265)
         mat(k,1120) = .500_r8*rxt(k,519)*y(k,133) + rxt(k,520)*y(k,265)
         mat(k,1981) = rxt(k,513)*y(k,217) + rxt(k,520)*y(k,220)
         mat(k,1106) = -(rxt(k,515)*y(k,241) + rxt(k,516)*y(k,247) + rxt(k,517) &
                      *y(k,131))
         mat(k,2324) = -rxt(k,515)*y(k,271)
         mat(k,2180) = -rxt(k,516)*y(k,271)
         mat(k,1759) = -rxt(k,517)*y(k,271)
         mat(k,1065) = rxt(k,501)*y(k,265)
         mat(k,1009) = rxt(k,504)*y(k,265)
         mat(k,529) = rxt(k,518)*y(k,265)
         mat(k,1971) = rxt(k,501)*y(k,6) + rxt(k,504)*y(k,116) + rxt(k,518)*y(k,219)
         mat(k,795) = -(rxt(k,486)*y(k,247) + rxt(k,487)*y(k,131))
         mat(k,2165) = -rxt(k,486)*y(k,272)
         mat(k,1743) = -rxt(k,487)*y(k,272)
         mat(k,665) = rxt(k,488)*y(k,265)
         mat(k,233) = (.650_r8*rxt(k,489)+rxt(k,579))*y(k,265)
         mat(k,1949) = rxt(k,488)*y(k,222) + (.650_r8*rxt(k,489)+rxt(k,579))*y(k,223)
         mat(k,1277) = -(rxt(k,450)*y(k,240) + rxt(k,451)*y(k,241) + rxt(k,452) &
                      *y(k,247) + rxt(k,453)*y(k,131) + rxt(k,454)*y(k,133))
         mat(k,1463) = -rxt(k,450)*y(k,273)
         mat(k,2336) = -rxt(k,451)*y(k,273)
         mat(k,2191) = -rxt(k,452)*y(k,273)
         mat(k,1772) = -rxt(k,453)*y(k,273)
         mat(k,2047) = -rxt(k,454)*y(k,273)
         mat(k,284) = rxt(k,422)*y(k,265)
         mat(k,399) = rxt(k,423)*y(k,265)
         mat(k,174) = rxt(k,424)*y(k,265)
         mat(k,756) = .400_r8*rxt(k,447)*y(k,265)
         mat(k,246) = .500_r8*rxt(k,455)*y(k,265)
         mat(k,1985) = rxt(k,422)*y(k,96) + rxt(k,423)*y(k,98) + rxt(k,424)*y(k,99) &
                      + .400_r8*rxt(k,447)*y(k,106) + .500_r8*rxt(k,455)*y(k,225)
         mat(k,811) = -(rxt(k,492)*y(k,247) + rxt(k,493)*y(k,131))
         mat(k,2166) = -rxt(k,492)*y(k,274)
         mat(k,1744) = -rxt(k,493)*y(k,274)
         mat(k,257) = .560_r8*rxt(k,491)*y(k,265)
         mat(k,768) = rxt(k,494)*y(k,265)
         mat(k,1950) = .560_r8*rxt(k,491)*y(k,226) + rxt(k,494)*y(k,227)
         mat(k,556) = -(rxt(k,495)*y(k,247) + rxt(k,496)*y(k,131))
         mat(k,2150) = -rxt(k,495)*y(k,275)
         mat(k,1731) = -rxt(k,496)*y(k,275)
         mat(k,264) = .300_r8*rxt(k,497)*y(k,265)
         mat(k,474) = rxt(k,498)*y(k,265)
         mat(k,1923) = .300_r8*rxt(k,497)*y(k,229) + rxt(k,498)*y(k,230)
         mat(k,2524) = -(rxt(k,174)*y(k,261) + rxt(k,293)*y(k,75) + rxt(k,538) &
                      *y(k,163))
         mat(k,1844) = -rxt(k,174)*y(k,276)
         mat(k,977) = -rxt(k,293)*y(k,276)
         mat(k,310) = -rxt(k,538)*y(k,276)
         mat(k,339) = rxt(k,346)*y(k,265)
         mat(k,442) = rxt(k,371)*y(k,265)
         mat(k,345) = rxt(k,372)*y(k,265)
         mat(k,520) = rxt(k,298)*y(k,265)
         mat(k,2497) = rxt(k,317)*y(k,265)
         mat(k,636) = rxt(k,300)*y(k,265)
         mat(k,168) = rxt(k,301)*y(k,265)
         mat(k,1203) = rxt(k,348)*y(k,265)
         mat(k,408) = rxt(k,303)*y(k,265)
         mat(k,1129) = rxt(k,386)*y(k,265)
         mat(k,1332) = rxt(k,374)*y(k,265)
         mat(k,753) = rxt(k,354)*y(k,265)
         mat(k,643) = rxt(k,355)*y(k,265)
         mat(k,460) = rxt(k,323)*y(k,265)
         mat(k,2102) = rxt(k,324)*y(k,265)
         mat(k,2289) = rxt(k,194)*y(k,247)
         mat(k,1516) = rxt(k,199)*y(k,265)
         mat(k,650) = rxt(k,200)*y(k,265)
         mat(k,887) = rxt(k,284)*y(k,265)
         mat(k,361) = rxt(k,308)*y(k,265)
         mat(k,1551) = (rxt(k,594)+rxt(k,599))*y(k,93) + (rxt(k,587)+rxt(k,593) &
                       +rxt(k,598))*y(k,94) + rxt(k,255)*y(k,265)
         mat(k,946) = rxt(k,326)*y(k,265)
         mat(k,2312) = rxt(k,231)*y(k,265)
         mat(k,527) = rxt(k,207)*y(k,265)
         mat(k,848) = (rxt(k,594)+rxt(k,599))*y(k,87)
         mat(k,879) = (rxt(k,587)+rxt(k,593)+rxt(k,598))*y(k,87) + rxt(k,258)*y(k,265)
         mat(k,1323) = .500_r8*rxt(k,399)*y(k,265)
         mat(k,143) = rxt(k,539)*y(k,265)
         mat(k,595) = rxt(k,380)*y(k,265)
         mat(k,430) = rxt(k,384)*y(k,265)
         mat(k,2222) = rxt(k,194)*y(k,78) + rxt(k,201)*y(k,265)
         mat(k,2017) = rxt(k,346)*y(k,30) + rxt(k,371)*y(k,32) + rxt(k,372)*y(k,33) &
                      + rxt(k,298)*y(k,43) + rxt(k,317)*y(k,44) + rxt(k,300)*y(k,45) &
                      + rxt(k,301)*y(k,46) + rxt(k,348)*y(k,47) + rxt(k,303)*y(k,48) &
                      + rxt(k,386)*y(k,50) + rxt(k,374)*y(k,51) + rxt(k,354)*y(k,52) &
                      + rxt(k,355)*y(k,53) + rxt(k,323)*y(k,55) + rxt(k,324)*y(k,56) &
                      + rxt(k,199)*y(k,79) + rxt(k,200)*y(k,81) + rxt(k,284)*y(k,83) &
                      + rxt(k,308)*y(k,86) + rxt(k,255)*y(k,87) + rxt(k,326)*y(k,89) &
                      + rxt(k,231)*y(k,91) + rxt(k,207)*y(k,92) + rxt(k,258)*y(k,94) &
                      + .500_r8*rxt(k,399)*y(k,111) + rxt(k,539)*y(k,127) + rxt(k,380) &
                      *y(k,157) + rxt(k,384)*y(k,158) + rxt(k,201)*y(k,247) &
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
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 73) = mat(k, 73) + lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = lmat(k, 78)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = lmat(k, 99)
         mat(k, 100) = lmat(k, 100)
         mat(k, 101) = lmat(k, 101)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 116) = lmat(k, 116)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 150) = lmat(k, 150)
         mat(k, 151) = lmat(k, 151)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 165) = mat(k, 165) + lmat(k, 165)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = lmat(k, 178)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 181) = lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 201) = mat(k, 201) + lmat(k, 201)
         mat(k, 202) = mat(k, 202) + lmat(k, 202)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = lmat(k, 206)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 214) = lmat(k, 214)
         mat(k, 215) = lmat(k, 215)
         mat(k, 216) = lmat(k, 216)
         mat(k, 217) = lmat(k, 217)
         mat(k, 218) = mat(k, 218) + lmat(k, 218)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 227) = lmat(k, 227)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 241) = lmat(k, 241)
         mat(k, 242) = lmat(k, 242)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 283) = lmat(k, 283)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 286) = lmat(k, 286)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 293) = mat(k, 293) + lmat(k, 293)
         mat(k, 295) = mat(k, 295) + lmat(k, 295)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 301) = mat(k, 301) + lmat(k, 301)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 308) = lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 316) = lmat(k, 316)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 331) = lmat(k, 331)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = lmat(k, 333)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 346) = lmat(k, 346)
         mat(k, 347) = lmat(k, 347)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 356) = mat(k, 356) + lmat(k, 356)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 366) = lmat(k, 366)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 371) = lmat(k, 371)
         mat(k, 372) = lmat(k, 372)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 377) = lmat(k, 377)
         mat(k, 378) = mat(k, 378) + lmat(k, 378)
         mat(k, 380) = mat(k, 380) + lmat(k, 380)
         mat(k, 381) = mat(k, 381) + lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 387) = lmat(k, 387)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 390) = lmat(k, 390)
         mat(k, 391) = lmat(k, 391)
         mat(k, 392) = mat(k, 392) + lmat(k, 392)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 407) = lmat(k, 407)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 418) = lmat(k, 418)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 427) = lmat(k, 427)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 429) = lmat(k, 429)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 432) = lmat(k, 432)
         mat(k, 434) = lmat(k, 434)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 436) = lmat(k, 436)
         mat(k, 437) = mat(k, 437) + lmat(k, 437)
         mat(k, 439) = lmat(k, 439)
         mat(k, 440) = mat(k, 440) + lmat(k, 440)
         mat(k, 441) = lmat(k, 441)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 444) = lmat(k, 444)
         mat(k, 446) = lmat(k, 446)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 450) = lmat(k, 450)
         mat(k, 452) = lmat(k, 452)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 454) = lmat(k, 454)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 457) = lmat(k, 457)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 461) = mat(k, 461) + lmat(k, 461)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 467) = lmat(k, 467)
         mat(k, 468) = lmat(k, 468)
         mat(k, 469) = lmat(k, 469)
         mat(k, 470) = lmat(k, 470)
         mat(k, 471) = lmat(k, 471)
         mat(k, 472) = lmat(k, 472)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 475) = lmat(k, 475)
         mat(k, 476) = lmat(k, 476)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 478) = lmat(k, 478)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 487) = mat(k, 487) + lmat(k, 487)
         mat(k, 491) = mat(k, 491) + lmat(k, 491)
         mat(k, 492) = lmat(k, 492)
         mat(k, 493) = mat(k, 493) + lmat(k, 493)
         mat(k, 498) = lmat(k, 498)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = lmat(k, 508)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 513) = lmat(k, 513)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 516) = mat(k, 516) + lmat(k, 516)
         mat(k, 521) = mat(k, 521) + lmat(k, 521)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 525) = lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 528) = mat(k, 528) + lmat(k, 528)
         mat(k, 530) = lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 533) = lmat(k, 533)
         mat(k, 534) = lmat(k, 534)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 565) = lmat(k, 565)
         mat(k, 566) = lmat(k, 566)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 580) = mat(k, 580) + lmat(k, 580)
         mat(k, 581) = lmat(k, 581)
         mat(k, 582) = lmat(k, 582)
         mat(k, 583) = lmat(k, 583)
         mat(k, 585) = mat(k, 585) + lmat(k, 585)
         mat(k, 586) = lmat(k, 586)
         mat(k, 588) = mat(k, 588) + lmat(k, 588)
         mat(k, 590) = lmat(k, 590)
         mat(k, 592) = mat(k, 592) + lmat(k, 592)
         mat(k, 593) = lmat(k, 593)
         mat(k, 594) = lmat(k, 594)
         mat(k, 596) = mat(k, 596) + lmat(k, 596)
         mat(k, 603) = lmat(k, 603)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 605) = lmat(k, 605)
         mat(k, 606) = lmat(k, 606)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 608) = lmat(k, 608)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 610) = lmat(k, 610)
         mat(k, 612) = mat(k, 612) + lmat(k, 612)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 630) = mat(k, 630) + lmat(k, 630)
         mat(k, 635) = lmat(k, 635)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 641) = lmat(k, 641)
         mat(k, 644) = mat(k, 644) + lmat(k, 644)
         mat(k, 646) = mat(k, 646) + lmat(k, 646)
         mat(k, 651) = mat(k, 651) + lmat(k, 651)
         mat(k, 658) = lmat(k, 658)
         mat(k, 660) = lmat(k, 660)
         mat(k, 661) = lmat(k, 661)
         mat(k, 662) = mat(k, 662) + lmat(k, 662)
         mat(k, 663) = lmat(k, 663)
         mat(k, 667) = lmat(k, 667)
         mat(k, 668) = lmat(k, 668)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 670) = lmat(k, 670)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 673) = lmat(k, 673)
         mat(k, 675) = lmat(k, 675)
         mat(k, 680) = lmat(k, 680)
         mat(k, 681) = lmat(k, 681)
         mat(k, 682) = lmat(k, 682)
         mat(k, 683) = lmat(k, 683)
         mat(k, 684) = mat(k, 684) + lmat(k, 684)
         mat(k, 685) = mat(k, 685) + lmat(k, 685)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 688) = lmat(k, 688)
         mat(k, 689) = lmat(k, 689)
         mat(k, 692) = mat(k, 692) + lmat(k, 692)
         mat(k, 698) = lmat(k, 698)
         mat(k, 699) = mat(k, 699) + lmat(k, 699)
         mat(k, 702) = mat(k, 702) + lmat(k, 702)
         mat(k, 703) = mat(k, 703) + lmat(k, 703)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 707) = lmat(k, 707)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 709) = lmat(k, 709)
         mat(k, 710) = mat(k, 710) + lmat(k, 710)
         mat(k, 714) = lmat(k, 714)
         mat(k, 715) = lmat(k, 715)
         mat(k, 717) = mat(k, 717) + lmat(k, 717)
         mat(k, 718) = lmat(k, 718)
         mat(k, 719) = lmat(k, 719)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 728) = mat(k, 728) + lmat(k, 728)
         mat(k, 729) = lmat(k, 729)
         mat(k, 730) = lmat(k, 730)
         mat(k, 731) = lmat(k, 731)
         mat(k, 732) = lmat(k, 732)
         mat(k, 734) = lmat(k, 734)
         mat(k, 735) = mat(k, 735) + lmat(k, 735)
         mat(k, 736) = lmat(k, 736)
         mat(k, 737) = lmat(k, 737)
         mat(k, 739) = mat(k, 739) + lmat(k, 739)
         mat(k, 749) = mat(k, 749) + lmat(k, 749)
         mat(k, 755) = mat(k, 755) + lmat(k, 755)
         mat(k, 757) = lmat(k, 757)
         mat(k, 758) = lmat(k, 758)
         mat(k, 759) = mat(k, 759) + lmat(k, 759)
         mat(k, 760) = lmat(k, 760)
         mat(k, 761) = lmat(k, 761)
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
         mat(k, 784) = mat(k, 784) + lmat(k, 784)
         mat(k, 795) = mat(k, 795) + lmat(k, 795)
         mat(k, 811) = mat(k, 811) + lmat(k, 811)
         mat(k, 822) = mat(k, 822) + lmat(k, 822)
         mat(k, 832) = mat(k, 832) + lmat(k, 832)
         mat(k, 841) = mat(k, 841) + lmat(k, 841)
         mat(k, 843) = lmat(k, 843)
         mat(k, 845) = mat(k, 845) + lmat(k, 845)
         mat(k, 849) = mat(k, 849) + lmat(k, 849)
         mat(k, 850) = lmat(k, 850)
         mat(k, 852) = lmat(k, 852)
         mat(k, 857) = mat(k, 857) + lmat(k, 857)
         mat(k, 862) = mat(k, 862) + lmat(k, 862)
         mat(k, 873) = mat(k, 873) + lmat(k, 873)
         mat(k, 876) = mat(k, 876) + lmat(k, 876)
         mat(k, 877) = mat(k, 877) + lmat(k, 877)
         mat(k, 880) = mat(k, 880) + lmat(k, 880)
         mat(k, 881) = mat(k, 881) + lmat(k, 881)
         mat(k, 885) = mat(k, 885) + lmat(k, 885)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 896) = mat(k, 896) + lmat(k, 896)
         mat(k, 897) = mat(k, 897) + lmat(k, 897)
         mat(k, 901) = mat(k, 901) + lmat(k, 901)
         mat(k, 909) = mat(k, 909) + lmat(k, 909)
         mat(k, 920) = mat(k, 920) + lmat(k, 920)
         mat(k, 936) = mat(k, 936) + lmat(k, 936)
         mat(k, 938) = lmat(k, 938)
         mat(k, 939) = mat(k, 939) + lmat(k, 939)
         mat(k, 941) = lmat(k, 941)
         mat(k, 942) = mat(k, 942) + lmat(k, 942)
         mat(k, 950) = mat(k, 950) + lmat(k, 950)
         mat(k, 959) = mat(k, 959) + lmat(k, 959)
         mat(k, 969) = mat(k, 969) + lmat(k, 969)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 991) = mat(k, 991) + lmat(k, 991)
         mat(k, 992) = lmat(k, 992)
         mat(k, 995) = lmat(k, 995)
         mat(k,1006) = mat(k,1006) + lmat(k,1006)
         mat(k,1026) = mat(k,1026) + lmat(k,1026)
         mat(k,1027) = mat(k,1027) + lmat(k,1027)
         mat(k,1029) = mat(k,1029) + lmat(k,1029)
         mat(k,1030) = lmat(k,1030)
         mat(k,1032) = mat(k,1032) + lmat(k,1032)
         mat(k,1033) = mat(k,1033) + lmat(k,1033)
         mat(k,1034) = mat(k,1034) + lmat(k,1034)
         mat(k,1038) = lmat(k,1038)
         mat(k,1042) = lmat(k,1042)
         mat(k,1043) = mat(k,1043) + lmat(k,1043)
         mat(k,1062) = mat(k,1062) + lmat(k,1062)
         mat(k,1086) = mat(k,1086) + lmat(k,1086)
         mat(k,1097) = lmat(k,1097)
         mat(k,1098) = mat(k,1098) + lmat(k,1098)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1102) = mat(k,1102) + lmat(k,1102)
         mat(k,1106) = mat(k,1106) + lmat(k,1106)
         mat(k,1116) = mat(k,1116) + lmat(k,1116)
         mat(k,1118) = lmat(k,1118)
         mat(k,1119) = lmat(k,1119)
         mat(k,1123) = lmat(k,1123)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1126) = lmat(k,1126)
         mat(k,1128) = lmat(k,1128)
         mat(k,1130) = lmat(k,1130)
         mat(k,1134) = mat(k,1134) + lmat(k,1134)
         mat(k,1139) = lmat(k,1139)
         mat(k,1142) = mat(k,1142) + lmat(k,1142)
         mat(k,1143) = lmat(k,1143)
         mat(k,1150) = mat(k,1150) + lmat(k,1150)
         mat(k,1162) = mat(k,1162) + lmat(k,1162)
         mat(k,1163) = mat(k,1163) + lmat(k,1163)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1165) = mat(k,1165) + lmat(k,1165)
         mat(k,1166) = mat(k,1166) + lmat(k,1166)
         mat(k,1167) = mat(k,1167) + lmat(k,1167)
         mat(k,1169) = mat(k,1169) + lmat(k,1169)
         mat(k,1171) = mat(k,1171) + lmat(k,1171)
         mat(k,1177) = mat(k,1177) + lmat(k,1177)
         mat(k,1195) = mat(k,1195) + lmat(k,1195)
         mat(k,1196) = lmat(k,1196)
         mat(k,1200) = lmat(k,1200)
         mat(k,1202) = lmat(k,1202)
         mat(k,1204) = mat(k,1204) + lmat(k,1204)
         mat(k,1209) = lmat(k,1209)
         mat(k,1210) = mat(k,1210) + lmat(k,1210)
         mat(k,1213) = mat(k,1213) + lmat(k,1213)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1222) = mat(k,1222) + lmat(k,1222)
         mat(k,1235) = lmat(k,1235)
         mat(k,1236) = lmat(k,1236)
         mat(k,1237) = lmat(k,1237)
         mat(k,1238) = lmat(k,1238)
         mat(k,1239) = mat(k,1239) + lmat(k,1239)
         mat(k,1240) = lmat(k,1240)
         mat(k,1242) = lmat(k,1242)
         mat(k,1244) = lmat(k,1244)
         mat(k,1247) = mat(k,1247) + lmat(k,1247)
         mat(k,1248) = lmat(k,1248)
         mat(k,1250) = lmat(k,1250)
         mat(k,1252) = mat(k,1252) + lmat(k,1252)
         mat(k,1254) = lmat(k,1254)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1257) = lmat(k,1257)
         mat(k,1261) = mat(k,1261) + lmat(k,1261)
         mat(k,1277) = mat(k,1277) + lmat(k,1277)
         mat(k,1297) = mat(k,1297) + lmat(k,1297)
         mat(k,1312) = mat(k,1312) + lmat(k,1312)
         mat(k,1313) = mat(k,1313) + lmat(k,1313)
         mat(k,1316) = mat(k,1316) + lmat(k,1316)
         mat(k,1317) = mat(k,1317) + lmat(k,1317)
         mat(k,1320) = mat(k,1320) + lmat(k,1320)
         mat(k,1322) = mat(k,1322) + lmat(k,1322)
         mat(k,1324) = mat(k,1324) + lmat(k,1324)
         mat(k,1325) = mat(k,1325) + lmat(k,1325)
         mat(k,1326) = mat(k,1326) + lmat(k,1326)
         mat(k,1330) = lmat(k,1330)
         mat(k,1343) = mat(k,1343) + lmat(k,1343)
         mat(k,1359) = lmat(k,1359)
         mat(k,1376) = mat(k,1376) + lmat(k,1376)
         mat(k,1386) = mat(k,1386) + lmat(k,1386)
         mat(k,1400) = mat(k,1400) + lmat(k,1400)
         mat(k,1414) = lmat(k,1414)
         mat(k,1416) = mat(k,1416) + lmat(k,1416)
         mat(k,1420) = mat(k,1420) + lmat(k,1420)
         mat(k,1422) = mat(k,1422) + lmat(k,1422)
         mat(k,1432) = lmat(k,1432)
         mat(k,1441) = mat(k,1441) + lmat(k,1441)
         mat(k,1472) = mat(k,1472) + lmat(k,1472)
         mat(k,1493) = mat(k,1493) + lmat(k,1493)
         mat(k,1494) = mat(k,1494) + lmat(k,1494)
         mat(k,1502) = lmat(k,1502)
         mat(k,1505) = mat(k,1505) + lmat(k,1505)
         mat(k,1518) = lmat(k,1518)
         mat(k,1520) = mat(k,1520) + lmat(k,1520)
         mat(k,1532) = mat(k,1532) + lmat(k,1532)
         mat(k,1539) = mat(k,1539) + lmat(k,1539)
         mat(k,1547) = mat(k,1547) + lmat(k,1547)
         mat(k,1548) = mat(k,1548) + lmat(k,1548)
         mat(k,1555) = mat(k,1555) + lmat(k,1555)
         mat(k,1575) = mat(k,1575) + lmat(k,1575)
         mat(k,1577) = mat(k,1577) + lmat(k,1577)
         mat(k,1578) = lmat(k,1578)
         mat(k,1586) = mat(k,1586) + lmat(k,1586)
         mat(k,1591) = mat(k,1591) + lmat(k,1591)
         mat(k,1597) = mat(k,1597) + lmat(k,1597)
         mat(k,1607) = mat(k,1607) + lmat(k,1607)
         mat(k,1609) = mat(k,1609) + lmat(k,1609)
         mat(k,1621) = mat(k,1621) + lmat(k,1621)
         mat(k,1636) = mat(k,1636) + lmat(k,1636)
         mat(k,1643) = mat(k,1643) + lmat(k,1643)
         mat(k,1648) = mat(k,1648) + lmat(k,1648)
         mat(k,1681) = mat(k,1681) + lmat(k,1681)
         mat(k,1682) = mat(k,1682) + lmat(k,1682)
         mat(k,1684) = mat(k,1684) + lmat(k,1684)
         mat(k,1690) = mat(k,1690) + lmat(k,1690)
         mat(k,1693) = mat(k,1693) + lmat(k,1693)
         mat(k,1748) = mat(k,1748) + lmat(k,1748)
         mat(k,1750) = lmat(k,1750)
         mat(k,1756) = mat(k,1756) + lmat(k,1756)
         mat(k,1788) = mat(k,1788) + lmat(k,1788)
         mat(k,1799) = mat(k,1799) + lmat(k,1799)
         mat(k,1832) = mat(k,1832) + lmat(k,1832)
         mat(k,1842) = mat(k,1842) + lmat(k,1842)
         mat(k,2006) = mat(k,2006) + lmat(k,2006)
         mat(k,2060) = mat(k,2060) + lmat(k,2060)
         mat(k,2063) = mat(k,2063) + lmat(k,2063)
         mat(k,2064) = mat(k,2064) + lmat(k,2064)
         mat(k,2067) = mat(k,2067) + lmat(k,2067)
         mat(k,2072) = mat(k,2072) + lmat(k,2072)
         mat(k,2075) = mat(k,2075) + lmat(k,2075)
         mat(k,2080) = lmat(k,2080)
         mat(k,2081) = lmat(k,2081)
         mat(k,2082) = mat(k,2082) + lmat(k,2082)
         mat(k,2091) = mat(k,2091) + lmat(k,2091)
         mat(k,2093) = mat(k,2093) + lmat(k,2093)
         mat(k,2096) = mat(k,2096) + lmat(k,2096)
         mat(k,2098) = mat(k,2098) + lmat(k,2098)
         mat(k,2100) = lmat(k,2100)
         mat(k,2101) = mat(k,2101) + lmat(k,2101)
         mat(k,2102) = mat(k,2102) + lmat(k,2102)
         mat(k,2214) = mat(k,2214) + lmat(k,2214)
         mat(k,2222) = mat(k,2222) + lmat(k,2222)
         mat(k,2261) = mat(k,2261) + lmat(k,2261)
         mat(k,2283) = mat(k,2283) + lmat(k,2283)
         mat(k,2298) = lmat(k,2298)
         mat(k,2301) = mat(k,2301) + lmat(k,2301)
         mat(k,2307) = mat(k,2307) + lmat(k,2307)
         mat(k,2360) = mat(k,2360) + lmat(k,2360)
         mat(k,2370) = mat(k,2370) + lmat(k,2370)
         mat(k,2411) = mat(k,2411) + lmat(k,2411)
         mat(k,2416) = mat(k,2416) + lmat(k,2416)
         mat(k,2425) = mat(k,2425) + lmat(k,2425)
         mat(k,2426) = mat(k,2426) + lmat(k,2426)
         mat(k,2439) = mat(k,2439) + lmat(k,2439)
         mat(k,2444) = lmat(k,2444)
         mat(k,2469) = mat(k,2469) + lmat(k,2469)
         mat(k,2475) = mat(k,2475) + lmat(k,2475)
         mat(k,2477) = lmat(k,2477)
         mat(k,2491) = mat(k,2491) + lmat(k,2491)
         mat(k,2496) = mat(k,2496) + lmat(k,2496)
         mat(k,2503) = lmat(k,2503)
         mat(k,2512) = mat(k,2512) + lmat(k,2512)
         mat(k,2513) = mat(k,2513) + lmat(k,2513)
         mat(k,2518) = lmat(k,2518)
         mat(k,2522) = lmat(k,2522)
         mat(k,2524) = mat(k,2524) + lmat(k,2524)
         mat(k, 258) = 0._r8
         mat(k, 259) = 0._r8
         mat(k, 298) = 0._r8
         mat(k, 357) = 0._r8
         mat(k, 375) = 0._r8
         mat(k, 482) = 0._r8
         mat(k, 485) = 0._r8
         mat(k, 504) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 561) = 0._r8
         mat(k, 664) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 695) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 701) = 0._r8
         mat(k, 704) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 712) = 0._r8
         mat(k, 716) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 747) = 0._r8
         mat(k, 767) = 0._r8
         mat(k, 769) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 774) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 799) = 0._r8
         mat(k, 803) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 812) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 815) = 0._r8
         mat(k, 817) = 0._r8
         mat(k, 820) = 0._r8
         mat(k, 833) = 0._r8
         mat(k, 834) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 868) = 0._r8
         mat(k, 871) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 912) = 0._r8
         mat(k, 913) = 0._r8
         mat(k, 914) = 0._r8
         mat(k, 955) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 996) = 0._r8
         mat(k,1007) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1016) = 0._r8
         mat(k,1024) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1044) = 0._r8
         mat(k,1045) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1063) = 0._r8
         mat(k,1064) = 0._r8
         mat(k,1072) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1084) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1132) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1172) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1224) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1230) = 0._r8
         mat(k,1233) = 0._r8
         mat(k,1241) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1318) = 0._r8
         mat(k,1327) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1351) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1362) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1370) = 0._r8
         mat(k,1372) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1374) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1384) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1412) = 0._r8
         mat(k,1417) = 0._r8
         mat(k,1421) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1429) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1439) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1474) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1481) = 0._r8
         mat(k,1485) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1510) = 0._r8
         mat(k,1513) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1540) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1546) = 0._r8
         mat(k,1554) = 0._r8
         mat(k,1558) = 0._r8
         mat(k,1559) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1563) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1567) = 0._r8
         mat(k,1581) = 0._r8
         mat(k,1583) = 0._r8
         mat(k,1585) = 0._r8
         mat(k,1587) = 0._r8
         mat(k,1588) = 0._r8
         mat(k,1589) = 0._r8
         mat(k,1592) = 0._r8
         mat(k,1594) = 0._r8
         mat(k,1598) = 0._r8
         mat(k,1599) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1615) = 0._r8
         mat(k,1618) = 0._r8
         mat(k,1619) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1623) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1641) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1650) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1663) = 0._r8
         mat(k,1668) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1670) = 0._r8
         mat(k,1671) = 0._r8
         mat(k,1674) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1683) = 0._r8
         mat(k,1686) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1691) = 0._r8
         mat(k,1694) = 0._r8
         mat(k,1695) = 0._r8
         mat(k,1751) = 0._r8
         mat(k,1789) = 0._r8
         mat(k,1792) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1796) = 0._r8
         mat(k,1801) = 0._r8
         mat(k,1830) = 0._r8
         mat(k,1834) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1937) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1951) = 0._r8
         mat(k,1961) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,2005) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2027) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2038) = 0._r8
         mat(k,2044) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2057) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2059) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2065) = 0._r8
         mat(k,2068) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2073) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2084) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2088) = 0._r8
         mat(k,2089) = 0._r8
         mat(k,2092) = 0._r8
         mat(k,2097) = 0._r8
         mat(k,2099) = 0._r8
         mat(k,2144) = 0._r8
         mat(k,2145) = 0._r8
         mat(k,2148) = 0._r8
         mat(k,2156) = 0._r8
         mat(k,2173) = 0._r8
         mat(k,2181) = 0._r8
         mat(k,2182) = 0._r8
         mat(k,2184) = 0._r8
         mat(k,2187) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2193) = 0._r8
         mat(k,2198) = 0._r8
         mat(k,2210) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2217) = 0._r8
         mat(k,2233) = 0._r8
         mat(k,2235) = 0._r8
         mat(k,2239) = 0._r8
         mat(k,2240) = 0._r8
         mat(k,2241) = 0._r8
         mat(k,2243) = 0._r8
         mat(k,2244) = 0._r8
         mat(k,2245) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2247) = 0._r8
         mat(k,2252) = 0._r8
         mat(k,2254) = 0._r8
         mat(k,2255) = 0._r8
         mat(k,2256) = 0._r8
         mat(k,2263) = 0._r8
         mat(k,2266) = 0._r8
         mat(k,2268) = 0._r8
         mat(k,2270) = 0._r8
         mat(k,2271) = 0._r8
         mat(k,2273) = 0._r8
         mat(k,2274) = 0._r8
         mat(k,2275) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2277) = 0._r8
         mat(k,2279) = 0._r8
         mat(k,2280) = 0._r8
         mat(k,2282) = 0._r8
         mat(k,2284) = 0._r8
         mat(k,2285) = 0._r8
         mat(k,2288) = 0._r8
         mat(k,2292) = 0._r8
         mat(k,2293) = 0._r8
         mat(k,2294) = 0._r8
         mat(k,2295) = 0._r8
         mat(k,2296) = 0._r8
         mat(k,2297) = 0._r8
         mat(k,2299) = 0._r8
         mat(k,2300) = 0._r8
         mat(k,2303) = 0._r8
         mat(k,2304) = 0._r8
         mat(k,2305) = 0._r8
         mat(k,2306) = 0._r8
         mat(k,2308) = 0._r8
         mat(k,2309) = 0._r8
         mat(k,2310) = 0._r8
         mat(k,2311) = 0._r8
         mat(k,2320) = 0._r8
         mat(k,2348) = 0._r8
         mat(k,2352) = 0._r8
         mat(k,2353) = 0._r8
         mat(k,2354) = 0._r8
         mat(k,2355) = 0._r8
         mat(k,2358) = 0._r8
         mat(k,2359) = 0._r8
         mat(k,2361) = 0._r8
         mat(k,2362) = 0._r8
         mat(k,2364) = 0._r8
         mat(k,2379) = 0._r8
         mat(k,2385) = 0._r8
         mat(k,2386) = 0._r8
         mat(k,2387) = 0._r8
         mat(k,2390) = 0._r8
         mat(k,2395) = 0._r8
         mat(k,2396) = 0._r8
         mat(k,2397) = 0._r8
         mat(k,2399) = 0._r8
         mat(k,2402) = 0._r8
         mat(k,2403) = 0._r8
         mat(k,2404) = 0._r8
         mat(k,2406) = 0._r8
         mat(k,2423) = 0._r8
         mat(k,2428) = 0._r8
         mat(k,2437) = 0._r8
         mat(k,2443) = 0._r8
         mat(k,2445) = 0._r8
         mat(k,2449) = 0._r8
         mat(k,2459) = 0._r8
         mat(k,2462) = 0._r8
         mat(k,2466) = 0._r8
         mat(k,2467) = 0._r8
         mat(k,2471) = 0._r8
         mat(k,2474) = 0._r8
         mat(k,2476) = 0._r8
         mat(k,2480) = 0._r8
         mat(k,2481) = 0._r8
         mat(k,2482) = 0._r8
         mat(k,2483) = 0._r8
         mat(k,2484) = 0._r8
         mat(k,2485) = 0._r8
         mat(k,2488) = 0._r8
         mat(k,2493) = 0._r8
         mat(k,2494) = 0._r8
         mat(k,2502) = 0._r8
         mat(k,2504) = 0._r8
         mat(k,2505) = 0._r8
         mat(k,2506) = 0._r8
         mat(k,2507) = 0._r8
         mat(k,2508) = 0._r8
         mat(k,2509) = 0._r8
         mat(k,2510) = 0._r8
         mat(k,2511) = 0._r8
         mat(k,2514) = 0._r8
         mat(k,2515) = 0._r8
         mat(k,2516) = 0._r8
         mat(k,2517) = 0._r8
         mat(k,2519) = 0._r8
         mat(k,2520) = 0._r8
         mat(k,2521) = 0._r8
         mat(k,2523) = 0._r8
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
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 78) = mat(k, 78) - dti(k)
         mat(k, 84) = mat(k, 84) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 98) = mat(k, 98) - dti(k)
         mat(k, 99) = mat(k, 99) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 107) = mat(k, 107) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 140) = mat(k, 140) - dti(k)
         mat(k, 141) = mat(k, 141) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 150) = mat(k, 150) - dti(k)
         mat(k, 153) = mat(k, 153) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 161) = mat(k, 161) - dti(k)
         mat(k, 165) = mat(k, 165) - dti(k)
         mat(k, 169) = mat(k, 169) - dti(k)
         mat(k, 173) = mat(k, 173) - dti(k)
         mat(k, 176) = mat(k, 176) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 182) = mat(k, 182) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 190) = mat(k, 190) - dti(k)
         mat(k, 195) = mat(k, 195) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 205) = mat(k, 205) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 214) = mat(k, 214) - dti(k)
         mat(k, 218) = mat(k, 218) - dti(k)
         mat(k, 223) = mat(k, 223) - dti(k)
         mat(k, 226) = mat(k, 226) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 236) = mat(k, 236) - dti(k)
         mat(k, 241) = mat(k, 241) - dti(k)
         mat(k, 245) = mat(k, 245) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 254) = mat(k, 254) - dti(k)
         mat(k, 262) = mat(k, 262) - dti(k)
         mat(k, 267) = mat(k, 267) - dti(k)
         mat(k, 270) = mat(k, 270) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 297) = mat(k, 297) - dti(k)
         mat(k, 302) = mat(k, 302) - dti(k)
         mat(k, 307) = mat(k, 307) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 315) = mat(k, 315) - dti(k)
         mat(k, 319) = mat(k, 319) - dti(k)
         mat(k, 325) = mat(k, 325) - dti(k)
         mat(k, 331) = mat(k, 331) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 346) = mat(k, 346) - dti(k)
         mat(k, 349) = mat(k, 349) - dti(k)
         mat(k, 356) = mat(k, 356) - dti(k)
         mat(k, 362) = mat(k, 362) - dti(k)
         mat(k, 367) = mat(k, 367) - dti(k)
         mat(k, 373) = mat(k, 373) - dti(k)
         mat(k, 378) = mat(k, 378) - dti(k)
         mat(k, 383) = mat(k, 383) - dti(k)
         mat(k, 388) = mat(k, 388) - dti(k)
         mat(k, 393) = mat(k, 393) - dti(k)
         mat(k, 398) = mat(k, 398) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 409) = mat(k, 409) - dti(k)
         mat(k, 417) = mat(k, 417) - dti(k)
         mat(k, 425) = mat(k, 425) - dti(k)
         mat(k, 431) = mat(k, 431) - dti(k)
         mat(k, 437) = mat(k, 437) - dti(k)
         mat(k, 443) = mat(k, 443) - dti(k)
         mat(k, 449) = mat(k, 449) - dti(k)
         mat(k, 455) = mat(k, 455) - dti(k)
         mat(k, 461) = mat(k, 461) - dti(k)
         mat(k, 467) = mat(k, 467) - dti(k)
         mat(k, 473) = mat(k, 473) - dti(k)
         mat(k, 481) = mat(k, 481) - dti(k)
         mat(k, 487) = mat(k, 487) - dti(k)
         mat(k, 493) = mat(k, 493) - dti(k)
         mat(k, 500) = mat(k, 500) - dti(k)
         mat(k, 506) = mat(k, 506) - dti(k)
         mat(k, 509) = mat(k, 509) - dti(k)
         mat(k, 514) = mat(k, 514) - dti(k)
         mat(k, 521) = mat(k, 521) - dti(k)
         mat(k, 528) = mat(k, 528) - dti(k)
         mat(k, 535) = mat(k, 535) - dti(k)
         mat(k, 541) = mat(k, 541) - dti(k)
         mat(k, 548) = mat(k, 548) - dti(k)
         mat(k, 556) = mat(k, 556) - dti(k)
         mat(k, 563) = mat(k, 563) - dti(k)
         mat(k, 569) = mat(k, 569) - dti(k)
         mat(k, 575) = mat(k, 575) - dti(k)
         mat(k, 580) = mat(k, 580) - dti(k)
         mat(k, 588) = mat(k, 588) - dti(k)
         mat(k, 596) = mat(k, 596) - dti(k)
         mat(k, 604) = mat(k, 604) - dti(k)
         mat(k, 612) = mat(k, 612) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 628) = mat(k, 628) - dti(k)
         mat(k, 637) = mat(k, 637) - dti(k)
         mat(k, 644) = mat(k, 644) - dti(k)
         mat(k, 651) = mat(k, 651) - dti(k)
         mat(k, 662) = mat(k, 662) - dti(k)
         mat(k, 671) = mat(k, 671) - dti(k)
         mat(k, 680) = mat(k, 680) - dti(k)
         mat(k, 684) = mat(k, 684) - dti(k)
         mat(k, 692) = mat(k, 692) - dti(k)
         mat(k, 699) = mat(k, 699) - dti(k)
         mat(k, 710) = mat(k, 710) - dti(k)
         mat(k, 721) = mat(k, 721) - dti(k)
         mat(k, 728) = mat(k, 728) - dti(k)
         mat(k, 739) = mat(k, 739) - dti(k)
         mat(k, 749) = mat(k, 749) - dti(k)
         mat(k, 755) = mat(k, 755) - dti(k)
         mat(k, 766) = mat(k, 766) - dti(k)
         mat(k, 777) = mat(k, 777) - dti(k)
         mat(k, 784) = mat(k, 784) - dti(k)
         mat(k, 795) = mat(k, 795) - dti(k)
         mat(k, 811) = mat(k, 811) - dti(k)
         mat(k, 822) = mat(k, 822) - dti(k)
         mat(k, 832) = mat(k, 832) - dti(k)
         mat(k, 841) = mat(k, 841) - dti(k)
         mat(k, 849) = mat(k, 849) - dti(k)
         mat(k, 857) = mat(k, 857) - dti(k)
         mat(k, 862) = mat(k, 862) - dti(k)
         mat(k, 873) = mat(k, 873) - dti(k)
         mat(k, 880) = mat(k, 880) - dti(k)
         mat(k, 889) = mat(k, 889) - dti(k)
         mat(k, 897) = mat(k, 897) - dti(k)
         mat(k, 909) = mat(k, 909) - dti(k)
         mat(k, 920) = mat(k, 920) - dti(k)
         mat(k, 936) = mat(k, 936) - dti(k)
         mat(k, 942) = mat(k, 942) - dti(k)
         mat(k, 950) = mat(k, 950) - dti(k)
         mat(k, 959) = mat(k, 959) - dti(k)
         mat(k, 969) = mat(k, 969) - dti(k)
         mat(k, 979) = mat(k, 979) - dti(k)
         mat(k, 991) = mat(k, 991) - dti(k)
         mat(k,1006) = mat(k,1006) - dti(k)
         mat(k,1027) = mat(k,1027) - dti(k)
         mat(k,1043) = mat(k,1043) - dti(k)
         mat(k,1062) = mat(k,1062) - dti(k)
         mat(k,1086) = mat(k,1086) - dti(k)
         mat(k,1098) = mat(k,1098) - dti(k)
         mat(k,1106) = mat(k,1106) - dti(k)
         mat(k,1116) = mat(k,1116) - dti(k)
         mat(k,1124) = mat(k,1124) - dti(k)
         mat(k,1134) = mat(k,1134) - dti(k)
         mat(k,1150) = mat(k,1150) - dti(k)
         mat(k,1163) = mat(k,1163) - dti(k)
         mat(k,1177) = mat(k,1177) - dti(k)
         mat(k,1195) = mat(k,1195) - dti(k)
         mat(k,1204) = mat(k,1204) - dti(k)
         mat(k,1210) = mat(k,1210) - dti(k)
         mat(k,1222) = mat(k,1222) - dti(k)
         mat(k,1239) = mat(k,1239) - dti(k)
         mat(k,1252) = mat(k,1252) - dti(k)
         mat(k,1261) = mat(k,1261) - dti(k)
         mat(k,1277) = mat(k,1277) - dti(k)
         mat(k,1297) = mat(k,1297) - dti(k)
         mat(k,1313) = mat(k,1313) - dti(k)
         mat(k,1325) = mat(k,1325) - dti(k)
         mat(k,1343) = mat(k,1343) - dti(k)
         mat(k,1376) = mat(k,1376) - dti(k)
         mat(k,1400) = mat(k,1400) - dti(k)
         mat(k,1420) = mat(k,1420) - dti(k)
         mat(k,1441) = mat(k,1441) - dti(k)
         mat(k,1472) = mat(k,1472) - dti(k)
         mat(k,1494) = mat(k,1494) - dti(k)
         mat(k,1505) = mat(k,1505) - dti(k)
         mat(k,1520) = mat(k,1520) - dti(k)
         mat(k,1539) = mat(k,1539) - dti(k)
         mat(k,1555) = mat(k,1555) - dti(k)
         mat(k,1586) = mat(k,1586) - dti(k)
         mat(k,1609) = mat(k,1609) - dti(k)
         mat(k,1636) = mat(k,1636) - dti(k)
         mat(k,1681) = mat(k,1681) - dti(k)
         mat(k,1788) = mat(k,1788) - dti(k)
         mat(k,1832) = mat(k,1832) - dti(k)
         mat(k,2006) = mat(k,2006) - dti(k)
         mat(k,2067) = mat(k,2067) - dti(k)
         mat(k,2093) = mat(k,2093) - dti(k)
         mat(k,2214) = mat(k,2214) - dti(k)
         mat(k,2261) = mat(k,2261) - dti(k)
         mat(k,2283) = mat(k,2283) - dti(k)
         mat(k,2307) = mat(k,2307) - dti(k)
         mat(k,2360) = mat(k,2360) - dti(k)
         mat(k,2425) = mat(k,2425) - dti(k)
         mat(k,2469) = mat(k,2469) - dti(k)
         mat(k,2496) = mat(k,2496) - dti(k)
         mat(k,2524) = mat(k,2524) - dti(k)
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
