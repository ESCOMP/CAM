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
         mat(k,757) = -(rxt(k,450)*y(k,248))
         mat(k,2098) = -rxt(k,450)*y(k,1)
         mat(k,2613) = rxt(k,453)*y(k,219)
         mat(k,995) = rxt(k,453)*y(k,153)
         mat(k,736) = -(rxt(k,454)*y(k,248))
         mat(k,2096) = -rxt(k,454)*y(k,2)
         mat(k,994) = rxt(k,451)*y(k,233)
         mat(k,2280) = rxt(k,451)*y(k,219)
         mat(k,1064) = -(rxt(k,533)*y(k,155) + rxt(k,534)*y(k,163) + rxt(k,535) &
                      *y(k,248))
         mat(k,1772) = -rxt(k,533)*y(k,6)
         mat(k,1906) = -rxt(k,534)*y(k,6)
         mat(k,2124) = -rxt(k,535)*y(k,6)
         mat(k,188) = -(rxt(k,492)*y(k,248))
         mat(k,2018) = -rxt(k,492)*y(k,7)
         mat(k,494) = -(rxt(k,495)*y(k,248))
         mat(k,2065) = -rxt(k,495)*y(k,8)
         mat(k,566) = rxt(k,493)*y(k,233)
         mat(k,2263) = rxt(k,493)*y(k,221)
         mat(k,189) = .120_r8*rxt(k,492)*y(k,248)
         mat(k,2019) = .120_r8*rxt(k,492)*y(k,7)
         mat(k,1059) = .100_r8*rxt(k,534)*y(k,163)
         mat(k,972) = .100_r8*rxt(k,537)*y(k,163)
         mat(k,1891) = .100_r8*rxt(k,534)*y(k,6) + .100_r8*rxt(k,537)*y(k,139)
         mat(k,2599) = .500_r8*rxt(k,494)*y(k,221) + .200_r8*rxt(k,521)*y(k,254) &
                      + .060_r8*rxt(k,527)*y(k,257)
         mat(k,567) = .500_r8*rxt(k,494)*y(k,153)
         mat(k,819) = .200_r8*rxt(k,521)*y(k,153)
         mat(k,835) = .060_r8*rxt(k,527)*y(k,153)
         mat(k,2593) = .200_r8*rxt(k,521)*y(k,254) + .200_r8*rxt(k,527)*y(k,257)
         mat(k,818) = .200_r8*rxt(k,521)*y(k,153)
         mat(k,833) = .200_r8*rxt(k,527)*y(k,153)
         mat(k,2610) = .200_r8*rxt(k,521)*y(k,254) + .150_r8*rxt(k,527)*y(k,257)
         mat(k,821) = .200_r8*rxt(k,521)*y(k,153)
         mat(k,836) = .150_r8*rxt(k,527)*y(k,153)
         mat(k,2595) = .210_r8*rxt(k,527)*y(k,257)
         mat(k,834) = .210_r8*rxt(k,527)*y(k,153)
         mat(k,277) = -(rxt(k,455)*y(k,248))
         mat(k,2034) = -rxt(k,455)*y(k,15)
         mat(k,1058) = .050_r8*rxt(k,534)*y(k,163)
         mat(k,971) = .050_r8*rxt(k,537)*y(k,163)
         mat(k,1890) = .050_r8*rxt(k,534)*y(k,6) + .050_r8*rxt(k,537)*y(k,139)
         mat(k,411) = -(rxt(k,421)*y(k,155) + rxt(k,422)*y(k,248))
         mat(k,1759) = -rxt(k,421)*y(k,16)
         mat(k,2054) = -rxt(k,422)*y(k,16)
         mat(k,2492) = -(rxt(k,243)*y(k,51) + rxt(k,244)*y(k,233) + rxt(k,245) &
                      *y(k,154) + rxt(k,246)*y(k,163) + rxt(k,253)*y(k,22) + rxt(k,282) &
                      *y(k,125))
         mat(k,1703) = -rxt(k,243)*y(k,17)
         mat(k,2345) = -rxt(k,244)*y(k,17)
         mat(k,2572) = -rxt(k,245)*y(k,17)
         mat(k,1946) = -rxt(k,246)*y(k,17)
         mat(k,914) = -rxt(k,253)*y(k,17)
         mat(k,2230) = -rxt(k,282)*y(k,17)
         mat(k,544) = rxt(k,242)*y(k,248)
         mat(k,2464) = 4.000_r8*rxt(k,247)*y(k,21) + (rxt(k,248)+rxt(k,249))*y(k,74) &
                      + rxt(k,556)*y(k,83) + rxt(k,272)*y(k,115) + (rxt(k,283) &
                       +rxt(k,284))*y(k,125) + rxt(k,252)*y(k,153) + rxt(k,257) &
                      *y(k,162) + rxt(k,567)*y(k,180) + rxt(k,258)*y(k,248)
         mat(k,171) = rxt(k,232)*y(k,247)
         mat(k,176) = rxt(k,262)*y(k,247)
         mat(k,557) = 2.000_r8*rxt(k,316)*y(k,70) + 2.000_r8*rxt(k,343)*y(k,247) &
                      + 2.000_r8*rxt(k,317)*y(k,248)
         mat(k,147) = rxt(k,318)*y(k,248)
         mat(k,672) = rxt(k,321)*y(k,70) + rxt(k,344)*y(k,247) + rxt(k,322)*y(k,248)
         mat(k,119) = 2.000_r8*rxt(k,328)*y(k,248)
         mat(k,505) = 3.000_r8*rxt(k,329)*y(k,70) + 3.000_r8*rxt(k,263)*y(k,247) &
                      + 3.000_r8*rxt(k,330)*y(k,248)
         mat(k,123) = rxt(k,331)*y(k,248)
         mat(k,1880) = 2.000_r8*rxt(k,316)*y(k,45) + rxt(k,321)*y(k,52) &
                      + 3.000_r8*rxt(k,329)*y(k,66)
         mat(k,2374) = (rxt(k,248)+rxt(k,249))*y(k,21)
         mat(k,1037) = rxt(k,556)*y(k,21)
         mat(k,127) = 2.000_r8*rxt(k,264)*y(k,247)
         mat(k,1517) = rxt(k,259)*y(k,162) + rxt(k,265)*y(k,247) + rxt(k,260)*y(k,248)
         mat(k,2402) = rxt(k,272)*y(k,21)
         mat(k,2230) = mat(k,2230) + (rxt(k,283)+rxt(k,284))*y(k,21)
         mat(k,2673) = rxt(k,252)*y(k,21)
         mat(k,2437) = rxt(k,257)*y(k,21) + rxt(k,259)*y(k,97)
         mat(k,1557) = rxt(k,567)*y(k,21)
         mat(k,1992) = rxt(k,232)*y(k,38) + rxt(k,262)*y(k,39) + 2.000_r8*rxt(k,343) &
                      *y(k,45) + rxt(k,344)*y(k,52) + 3.000_r8*rxt(k,263)*y(k,66) &
                      + 2.000_r8*rxt(k,264)*y(k,94) + rxt(k,265)*y(k,97)
         mat(k,2174) = rxt(k,242)*y(k,18) + rxt(k,258)*y(k,21) + 2.000_r8*rxt(k,317) &
                      *y(k,45) + rxt(k,318)*y(k,46) + rxt(k,322)*y(k,52) &
                      + 2.000_r8*rxt(k,328)*y(k,65) + 3.000_r8*rxt(k,330)*y(k,66) &
                      + rxt(k,331)*y(k,67) + rxt(k,260)*y(k,97)
         mat(k,541) = -(rxt(k,242)*y(k,248))
         mat(k,2071) = -rxt(k,242)*y(k,18)
         mat(k,2470) = rxt(k,253)*y(k,22)
         mat(k,904) = rxt(k,253)*y(k,17)
         mat(k,1505) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,108)
         mat(k,1657) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,97)
         mat(k,2444) = rxt(k,250)*y(k,74)
         mat(k,905) = rxt(k,254)*y(k,70)
         mat(k,1838) = rxt(k,254)*y(k,22)
         mat(k,2354) = rxt(k,250)*y(k,21)
         mat(k,1506) = (rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,109)
         mat(k,1734) = (rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,108)
         mat(k,1658) = (rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,101)
         mat(k,1710) = (rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,97)
         mat(k,2469) = rxt(k,245)*y(k,154)
         mat(k,2523) = rxt(k,245)*y(k,17)
         mat(k,2463) = -(4._r8*rxt(k,247)*y(k,21) + (rxt(k,248) + rxt(k,249) + rxt(k,250) &
                      ) * y(k,74) + rxt(k,251)*y(k,233) + rxt(k,252)*y(k,153) &
                      + rxt(k,255)*y(k,154) + rxt(k,257)*y(k,162) + rxt(k,258) &
                      *y(k,248) + rxt(k,272)*y(k,115) + (rxt(k,283) + rxt(k,284) &
                      ) * y(k,125) + rxt(k,556)*y(k,83) + rxt(k,567)*y(k,180))
         mat(k,2373) = -(rxt(k,248) + rxt(k,249) + rxt(k,250)) * y(k,21)
         mat(k,2344) = -rxt(k,251)*y(k,21)
         mat(k,2672) = -rxt(k,252)*y(k,21)
         mat(k,2571) = -rxt(k,255)*y(k,21)
         mat(k,2436) = -rxt(k,257)*y(k,21)
         mat(k,2173) = -rxt(k,258)*y(k,21)
         mat(k,2401) = -rxt(k,272)*y(k,21)
         mat(k,2229) = -(rxt(k,283) + rxt(k,284)) * y(k,21)
         mat(k,1036) = -rxt(k,556)*y(k,21)
         mat(k,1556) = -rxt(k,567)*y(k,21)
         mat(k,2491) = rxt(k,282)*y(k,125) + rxt(k,246)*y(k,163)
         mat(k,913) = rxt(k,256)*y(k,162)
         mat(k,1516) = rxt(k,266)*y(k,247)
         mat(k,1674) = rxt(k,261)*y(k,162)
         mat(k,2229) = mat(k,2229) + rxt(k,282)*y(k,17)
         mat(k,2436) = mat(k,2436) + rxt(k,256)*y(k,22) + rxt(k,261)*y(k,108)
         mat(k,1945) = rxt(k,246)*y(k,17)
         mat(k,1991) = rxt(k,266)*y(k,97)
         mat(k,906) = -(rxt(k,253)*y(k,17) + rxt(k,254)*y(k,70) + rxt(k,256)*y(k,162))
         mat(k,2472) = -rxt(k,253)*y(k,22)
         mat(k,1845) = -rxt(k,254)*y(k,22)
         mat(k,2410) = -rxt(k,256)*y(k,22)
         mat(k,2445) = rxt(k,255)*y(k,154)
         mat(k,2541) = rxt(k,255)*y(k,21)
         mat(k,280) = -(rxt(k,496)*y(k,248))
         mat(k,2035) = -rxt(k,496)*y(k,24)
         mat(k,2590) = rxt(k,499)*y(k,223)
         mat(k,512) = rxt(k,499)*y(k,153)
         mat(k,382) = -(rxt(k,498)*y(k,248))
         mat(k,2050) = -rxt(k,498)*y(k,25)
         mat(k,513) = rxt(k,497)*y(k,233)
         mat(k,2254) = rxt(k,497)*y(k,223)
         mat(k,216) = -(rxt(k,312)*y(k,70) + rxt(k,313)*y(k,248))
         mat(k,1826) = -rxt(k,312)*y(k,26)
         mat(k,2022) = -rxt(k,313)*y(k,26)
         mat(k,321) = -(rxt(k,369)*y(k,70) + rxt(k,370)*y(k,248))
         mat(k,1828) = -rxt(k,369)*y(k,27)
         mat(k,2042) = -rxt(k,370)*y(k,27)
         mat(k,641) = -(rxt(k,371)*y(k,70) + rxt(k,372)*y(k,163) + rxt(k,397)*y(k,248))
         mat(k,1840) = -rxt(k,371)*y(k,28)
         mat(k,1895) = -rxt(k,372)*y(k,28)
         mat(k,2084) = -rxt(k,397)*y(k,28)
         mat(k,286) = -(rxt(k,314)*y(k,70) + rxt(k,315)*y(k,248))
         mat(k,1827) = -rxt(k,314)*y(k,29)
         mat(k,2037) = -rxt(k,315)*y(k,29)
         mat(k,296) = -(rxt(k,377)*y(k,248))
         mat(k,2039) = -rxt(k,377)*y(k,30)
         mat(k,875) = .800_r8*rxt(k,373)*y(k,224) + .200_r8*rxt(k,374)*y(k,228)
         mat(k,1602) = .200_r8*rxt(k,374)*y(k,224)
         mat(k,387) = -(rxt(k,378)*y(k,248))
         mat(k,2051) = -rxt(k,378)*y(k,31)
         mat(k,876) = rxt(k,375)*y(k,233)
         mat(k,2255) = rxt(k,375)*y(k,224)
         mat(k,327) = -(rxt(k,379)*y(k,70) + rxt(k,380)*y(k,248))
         mat(k,1829) = -rxt(k,379)*y(k,32)
         mat(k,2043) = -rxt(k,380)*y(k,32)
         mat(k,1129) = -(rxt(k,400)*y(k,155) + rxt(k,401)*y(k,163) + rxt(k,419) &
                      *y(k,248))
         mat(k,1777) = -rxt(k,400)*y(k,33)
         mat(k,1910) = -rxt(k,401)*y(k,33)
         mat(k,2129) = -rxt(k,419)*y(k,33)
         mat(k,921) = .130_r8*rxt(k,479)*y(k,163)
         mat(k,1910) = mat(k,1910) + .130_r8*rxt(k,479)*y(k,127)
         mat(k,488) = -(rxt(k,405)*y(k,248))
         mat(k,2064) = -rxt(k,405)*y(k,34)
         mat(k,935) = rxt(k,403)*y(k,233)
         mat(k,2262) = rxt(k,403)*y(k,225)
         mat(k,333) = -(rxt(k,406)*y(k,248) + rxt(k,409)*y(k,70))
         mat(k,2044) = -rxt(k,406)*y(k,35)
         mat(k,1830) = -rxt(k,409)*y(k,35)
         mat(k,300) = -(rxt(k,502)*y(k,248))
         mat(k,2040) = -rxt(k,502)*y(k,36)
         mat(k,727) = rxt(k,500)*y(k,233)
         mat(k,2249) = rxt(k,500)*y(k,226)
         mat(k,113) = -(rxt(k,231)*y(k,247))
         mat(k,1951) = -rxt(k,231)*y(k,37)
         mat(k,167) = -(rxt(k,232)*y(k,247))
         mat(k,1956) = -rxt(k,232)*y(k,38)
         mat(k,172) = -(rxt(k,262)*y(k,247))
         mat(k,1957) = -rxt(k,262)*y(k,39)
         mat(k,132) = -(rxt(k,233)*y(k,247))
         mat(k,1953) = -rxt(k,233)*y(k,40)
         mat(k,177) = -(rxt(k,234)*y(k,247))
         mat(k,1958) = -rxt(k,234)*y(k,41)
         mat(k,136) = -(rxt(k,235)*y(k,247))
         mat(k,1954) = -rxt(k,235)*y(k,42)
         mat(k,182) = -(rxt(k,236)*y(k,247))
         mat(k,1959) = -rxt(k,236)*y(k,43)
         mat(k,140) = -(rxt(k,237)*y(k,247))
         mat(k,1955) = -rxt(k,237)*y(k,44)
         mat(k,552) = -(rxt(k,316)*y(k,70) + rxt(k,317)*y(k,248) + rxt(k,343)*y(k,247))
         mat(k,1837) = -rxt(k,316)*y(k,45)
         mat(k,2073) = -rxt(k,317)*y(k,45)
         mat(k,1968) = -rxt(k,343)*y(k,45)
         mat(k,144) = -(rxt(k,318)*y(k,248))
         mat(k,2015) = -rxt(k,318)*y(k,46)
         mat(k,342) = -(rxt(k,319)*y(k,70) + rxt(k,320)*y(k,248))
         mat(k,1831) = -rxt(k,319)*y(k,47)
         mat(k,2045) = -rxt(k,320)*y(k,47)
         mat(k,1688) = -(rxt(k,204)*y(k,70) + rxt(k,243)*y(k,17) + rxt(k,348)*y(k,233) &
                      + rxt(k,349)*y(k,155) + rxt(k,350)*y(k,162) + rxt(k,351) &
                      *y(k,248))
         mat(k,1865) = -rxt(k,204)*y(k,51)
         mat(k,2477) = -rxt(k,243)*y(k,51)
         mat(k,2330) = -rxt(k,348)*y(k,51)
         mat(k,1805) = -rxt(k,349)*y(k,51)
         mat(k,2422) = -rxt(k,350)*y(k,51)
         mat(k,2159) = -rxt(k,351)*y(k,51)
         mat(k,763) = .400_r8*rxt(k,450)*y(k,248)
         mat(k,1076) = .340_r8*rxt(k,534)*y(k,163)
         mat(k,415) = .500_r8*rxt(k,421)*y(k,155)
         mat(k,645) = rxt(k,372)*y(k,163)
         mat(k,1137) = .500_r8*rxt(k,401)*y(k,163)
         mat(k,678) = .500_r8*rxt(k,389)*y(k,248)
         mat(k,869) = rxt(k,356)*y(k,248)
         mat(k,478) = .300_r8*rxt(k,357)*y(k,248)
         mat(k,1569) = (rxt(k,365)+rxt(k,366))*y(k,247)
         mat(k,1116) = rxt(k,332)*y(k,228)
         mat(k,2359) = rxt(k,213)*y(k,228)
         mat(k,1208) = .800_r8*rxt(k,394)*y(k,248)
         mat(k,930) = .910_r8*rxt(k,479)*y(k,163)
         mat(k,694) = .300_r8*rxt(k,470)*y(k,248)
         mat(k,1332) = .120_r8*rxt(k,432)*y(k,163)
         mat(k,685) = .500_r8*rxt(k,445)*y(k,248)
         mat(k,987) = .340_r8*rxt(k,537)*y(k,163)
         mat(k,1442) = .600_r8*rxt(k,446)*y(k,163)
         mat(k,2658) = .100_r8*rxt(k,452)*y(k,219) + rxt(k,355)*y(k,228) &
                      + .500_r8*rxt(k,423)*y(k,230) + .500_r8*rxt(k,391)*y(k,232) &
                      + .920_r8*rxt(k,462)*y(k,235) + .250_r8*rxt(k,430)*y(k,240) &
                      + rxt(k,439)*y(k,242) + rxt(k,413)*y(k,250) + rxt(k,417) &
                      *y(k,251) + .340_r8*rxt(k,546)*y(k,252) + .320_r8*rxt(k,551) &
                      *y(k,253) + .250_r8*rxt(k,487)*y(k,256)
         mat(k,1805) = mat(k,1805) + .500_r8*rxt(k,421)*y(k,16) + rxt(k,463)*y(k,235) &
                      + .250_r8*rxt(k,429)*y(k,240) + rxt(k,440)*y(k,242)
         mat(k,1931) = .340_r8*rxt(k,534)*y(k,6) + rxt(k,372)*y(k,28) &
                      + .500_r8*rxt(k,401)*y(k,33) + .910_r8*rxt(k,479)*y(k,127) &
                      + .120_r8*rxt(k,432)*y(k,134) + .340_r8*rxt(k,537)*y(k,139) &
                      + .600_r8*rxt(k,446)*y(k,140)
         mat(k,629) = rxt(k,396)*y(k,248)
         mat(k,1185) = .680_r8*rxt(k,555)*y(k,248)
         mat(k,1003) = .100_r8*rxt(k,452)*y(k,153)
         mat(k,881) = .700_r8*rxt(k,374)*y(k,228)
         mat(k,940) = rxt(k,402)*y(k,228)
         mat(k,1494) = rxt(k,385)*y(k,228) + rxt(k,459)*y(k,235) + .250_r8*rxt(k,426) &
                      *y(k,240) + rxt(k,435)*y(k,242) + .250_r8*rxt(k,484)*y(k,256)
         mat(k,1640) = rxt(k,332)*y(k,68) + rxt(k,213)*y(k,74) + rxt(k,355)*y(k,153) &
                      + .700_r8*rxt(k,374)*y(k,224) + rxt(k,402)*y(k,225) + rxt(k,385) &
                      *y(k,227) + (4.000_r8*rxt(k,352)+2.000_r8*rxt(k,353))*y(k,228) &
                      + 1.500_r8*rxt(k,460)*y(k,235) + .750_r8*rxt(k,465)*y(k,236) &
                      + .800_r8*rxt(k,474)*y(k,237) + .880_r8*rxt(k,427)*y(k,240) &
                      + 2.000_r8*rxt(k,436)*y(k,242) + .750_r8*rxt(k,539)*y(k,246) &
                      + .800_r8*rxt(k,415)*y(k,251) + .930_r8*rxt(k,544)*y(k,252) &
                      + .950_r8*rxt(k,549)*y(k,253) + .800_r8*rxt(k,485)*y(k,256)
         mat(k,653) = .500_r8*rxt(k,423)*y(k,153)
         mat(k,799) = .500_r8*rxt(k,391)*y(k,153)
         mat(k,2330) = mat(k,2330) + .450_r8*rxt(k,437)*y(k,242) + .150_r8*rxt(k,416) &
                      *y(k,251)
         mat(k,1365) = .920_r8*rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155) + rxt(k,459) &
                      *y(k,227) + 1.500_r8*rxt(k,460)*y(k,228)
         mat(k,1398) = .750_r8*rxt(k,465)*y(k,228)
         mat(k,1317) = .800_r8*rxt(k,474)*y(k,228)
         mat(k,1420) = .250_r8*rxt(k,430)*y(k,153) + .250_r8*rxt(k,429)*y(k,155) &
                      + .250_r8*rxt(k,426)*y(k,227) + .880_r8*rxt(k,427)*y(k,228)
         mat(k,1462) = rxt(k,439)*y(k,153) + rxt(k,440)*y(k,155) + rxt(k,435)*y(k,227) &
                      + 2.000_r8*rxt(k,436)*y(k,228) + .450_r8*rxt(k,437)*y(k,233) &
                      + 4.000_r8*rxt(k,438)*y(k,242)
         mat(k,1171) = .750_r8*rxt(k,539)*y(k,228)
         mat(k,1977) = (rxt(k,365)+rxt(k,366))*y(k,64)
         mat(k,2159) = mat(k,2159) + .400_r8*rxt(k,450)*y(k,1) + .500_r8*rxt(k,389) &
                      *y(k,60) + rxt(k,356)*y(k,62) + .300_r8*rxt(k,357)*y(k,63) &
                      + .800_r8*rxt(k,394)*y(k,90) + .300_r8*rxt(k,470)*y(k,128) &
                      + .500_r8*rxt(k,445)*y(k,138) + rxt(k,396)*y(k,169) &
                      + .680_r8*rxt(k,555)*y(k,208)
         mat(k,862) = rxt(k,413)*y(k,153)
         mat(k,1278) = rxt(k,417)*y(k,153) + .800_r8*rxt(k,415)*y(k,228) &
                      + .150_r8*rxt(k,416)*y(k,233)
         mat(k,1224) = .340_r8*rxt(k,546)*y(k,153) + .930_r8*rxt(k,544)*y(k,228)
         mat(k,1092) = .320_r8*rxt(k,551)*y(k,153) + .950_r8*rxt(k,549)*y(k,228)
         mat(k,1295) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,484)*y(k,227) &
                      + .800_r8*rxt(k,485)*y(k,228)
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
         mat(k,665) = -(rxt(k,321)*y(k,70) + rxt(k,322)*y(k,248) + rxt(k,344)*y(k,247))
         mat(k,1841) = -rxt(k,321)*y(k,52)
         mat(k,2087) = -rxt(k,322)*y(k,52)
         mat(k,1969) = -rxt(k,344)*y(k,52)
         mat(k,148) = -(rxt(k,323)*y(k,248))
         mat(k,2016) = -rxt(k,323)*y(k,53)
         mat(k,1190) = -(rxt(k,381)*y(k,155) + rxt(k,382)*y(k,248))
         mat(k,1781) = -rxt(k,381)*y(k,54)
         mat(k,2133) = -rxt(k,382)*y(k,54)
         mat(k,761) = .800_r8*rxt(k,450)*y(k,248)
         mat(k,414) = rxt(k,421)*y(k,155)
         mat(k,297) = rxt(k,377)*y(k,248)
         mat(k,389) = .500_r8*rxt(k,378)*y(k,248)
         mat(k,1130) = .500_r8*rxt(k,401)*y(k,163)
         mat(k,1432) = .100_r8*rxt(k,446)*y(k,163)
         mat(k,2636) = .400_r8*rxt(k,452)*y(k,219) + rxt(k,376)*y(k,224) &
                      + .270_r8*rxt(k,404)*y(k,225) + rxt(k,423)*y(k,230) + rxt(k,442) &
                      *y(k,244) + rxt(k,413)*y(k,250)
         mat(k,1781) = mat(k,1781) + rxt(k,421)*y(k,16)
         mat(k,1913) = .500_r8*rxt(k,401)*y(k,33) + .100_r8*rxt(k,446)*y(k,140)
         mat(k,1000) = .400_r8*rxt(k,452)*y(k,153)
         mat(k,879) = rxt(k,376)*y(k,153) + 3.200_r8*rxt(k,373)*y(k,224) &
                      + .800_r8*rxt(k,374)*y(k,228)
         mat(k,938) = .270_r8*rxt(k,404)*y(k,153)
         mat(k,1622) = .800_r8*rxt(k,374)*y(k,224)
         mat(k,651) = rxt(k,423)*y(k,153)
         mat(k,2306) = .200_r8*rxt(k,441)*y(k,244)
         mat(k,769) = rxt(k,442)*y(k,153) + .200_r8*rxt(k,441)*y(k,233)
         mat(k,2133) = mat(k,2133) + .800_r8*rxt(k,450)*y(k,1) + rxt(k,377)*y(k,30) &
                      + .500_r8*rxt(k,378)*y(k,31)
         mat(k,860) = rxt(k,413)*y(k,153)
         mat(k,438) = -(rxt(k,324)*y(k,70) + rxt(k,325)*y(k,248))
         mat(k,1835) = -rxt(k,324)*y(k,55)
         mat(k,2057) = -rxt(k,325)*y(k,55)
         mat(k,107) = -(rxt(k,383)*y(k,248))
         mat(k,2011) = -rxt(k,383)*y(k,56)
         mat(k,1106) = -(rxt(k,420)*y(k,248))
         mat(k,2127) = -rxt(k,420)*y(k,57)
         mat(k,760) = .800_r8*rxt(k,450)*y(k,248)
         mat(k,1067) = .520_r8*rxt(k,534)*y(k,163)
         mat(k,413) = .500_r8*rxt(k,421)*y(k,155)
         mat(k,980) = .520_r8*rxt(k,537)*y(k,163)
         mat(k,2631) = .250_r8*rxt(k,452)*y(k,219) + .820_r8*rxt(k,404)*y(k,225) &
                      + .500_r8*rxt(k,423)*y(k,230) + .270_r8*rxt(k,546)*y(k,252) &
                      + .040_r8*rxt(k,551)*y(k,253)
         mat(k,1775) = .500_r8*rxt(k,421)*y(k,16)
         mat(k,1909) = .520_r8*rxt(k,534)*y(k,6) + .520_r8*rxt(k,537)*y(k,139)
         mat(k,1179) = .500_r8*rxt(k,555)*y(k,248)
         mat(k,999) = .250_r8*rxt(k,452)*y(k,153)
         mat(k,937) = .820_r8*rxt(k,404)*y(k,153) + .820_r8*rxt(k,402)*y(k,228)
         mat(k,1617) = .820_r8*rxt(k,402)*y(k,225) + .150_r8*rxt(k,544)*y(k,252) &
                      + .025_r8*rxt(k,549)*y(k,253)
         mat(k,650) = .500_r8*rxt(k,423)*y(k,153)
         mat(k,2127) = mat(k,2127) + .800_r8*rxt(k,450)*y(k,1) + .500_r8*rxt(k,555) &
                      *y(k,208)
         mat(k,1216) = .270_r8*rxt(k,546)*y(k,153) + .150_r8*rxt(k,544)*y(k,228)
         mat(k,1090) = .040_r8*rxt(k,551)*y(k,153) + .025_r8*rxt(k,549)*y(k,228)
         mat(k,1339) = -(rxt(k,407)*y(k,155) + rxt(k,408)*y(k,248))
         mat(k,1792) = -rxt(k,407)*y(k,58)
         mat(k,2144) = -rxt(k,408)*y(k,58)
         mat(k,1251) = rxt(k,410)*y(k,248)
         mat(k,1328) = .880_r8*rxt(k,432)*y(k,163)
         mat(k,1435) = .500_r8*rxt(k,446)*y(k,163)
         mat(k,2646) = .170_r8*rxt(k,505)*y(k,229) + .050_r8*rxt(k,468)*y(k,236) &
                      + .250_r8*rxt(k,430)*y(k,240) + .170_r8*rxt(k,511)*y(k,243) &
                      + .400_r8*rxt(k,521)*y(k,254) + .250_r8*rxt(k,487)*y(k,256) &
                      + .540_r8*rxt(k,527)*y(k,257) + .510_r8*rxt(k,530)*y(k,259)
         mat(k,1792) = mat(k,1792) + .050_r8*rxt(k,469)*y(k,236) + .250_r8*rxt(k,429) &
                      *y(k,240) + .250_r8*rxt(k,488)*y(k,256)
         mat(k,899) = rxt(k,411)*y(k,248)
         mat(k,1921) = .880_r8*rxt(k,432)*y(k,134) + .500_r8*rxt(k,446)*y(k,140)
         mat(k,1485) = .250_r8*rxt(k,426)*y(k,240) + .250_r8*rxt(k,484)*y(k,256)
         mat(k,1631) = .240_r8*rxt(k,427)*y(k,240) + .500_r8*rxt(k,415)*y(k,251) &
                      + .100_r8*rxt(k,485)*y(k,256)
         mat(k,852) = .170_r8*rxt(k,505)*y(k,153) + .070_r8*rxt(k,504)*y(k,233)
         mat(k,2316) = .070_r8*rxt(k,504)*y(k,229) + .070_r8*rxt(k,510)*y(k,243)
         mat(k,1391) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,1415) = .250_r8*rxt(k,430)*y(k,153) + .250_r8*rxt(k,429)*y(k,155) &
                      + .250_r8*rxt(k,426)*y(k,227) + .240_r8*rxt(k,427)*y(k,228)
         mat(k,959) = .170_r8*rxt(k,511)*y(k,153) + .070_r8*rxt(k,510)*y(k,233)
         mat(k,2144) = mat(k,2144) + rxt(k,410)*y(k,113) + rxt(k,411)*y(k,156)
         mat(k,1275) = .500_r8*rxt(k,415)*y(k,228)
         mat(k,828) = .400_r8*rxt(k,521)*y(k,153)
         mat(k,1292) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,227) + .100_r8*rxt(k,485)*y(k,228)
         mat(k,844) = .540_r8*rxt(k,527)*y(k,153)
         mat(k,586) = .510_r8*rxt(k,530)*y(k,153)
         mat(k,775) = -(rxt(k,388)*y(k,248))
         mat(k,2100) = -rxt(k,388)*y(k,59)
         mat(k,1124) = .120_r8*rxt(k,401)*y(k,163)
         mat(k,1897) = .120_r8*rxt(k,401)*y(k,33)
         mat(k,1475) = .100_r8*rxt(k,385)*y(k,228) + .150_r8*rxt(k,386)*y(k,233)
         mat(k,1608) = .100_r8*rxt(k,385)*y(k,227)
         mat(k,2283) = .150_r8*rxt(k,386)*y(k,227) + .150_r8*rxt(k,437)*y(k,242)
         mat(k,1454) = .150_r8*rxt(k,437)*y(k,233)
         mat(k,674) = -(rxt(k,389)*y(k,248))
         mat(k,2088) = -rxt(k,389)*y(k,60)
         mat(k,1474) = .400_r8*rxt(k,386)*y(k,233)
         mat(k,2275) = .400_r8*rxt(k,386)*y(k,227) + .400_r8*rxt(k,437)*y(k,242)
         mat(k,1452) = .400_r8*rxt(k,437)*y(k,233)
         mat(k,430) = -(rxt(k,326)*y(k,70) + rxt(k,327)*y(k,248))
         mat(k,1834) = -rxt(k,326)*y(k,61)
         mat(k,2056) = -rxt(k,327)*y(k,61)
         mat(k,868) = -(rxt(k,356)*y(k,248))
         mat(k,2109) = -rxt(k,356)*y(k,62)
         mat(k,877) = .300_r8*rxt(k,374)*y(k,228)
         mat(k,1609) = .300_r8*rxt(k,374)*y(k,224) + 2.000_r8*rxt(k,353)*y(k,228) &
                      + .250_r8*rxt(k,460)*y(k,235) + .250_r8*rxt(k,465)*y(k,236) &
                      + .200_r8*rxt(k,474)*y(k,237) + .250_r8*rxt(k,427)*y(k,240) &
                      + .250_r8*rxt(k,539)*y(k,246) + .500_r8*rxt(k,415)*y(k,251) &
                      + .250_r8*rxt(k,544)*y(k,252) + .250_r8*rxt(k,549)*y(k,253) &
                      + .300_r8*rxt(k,485)*y(k,256)
         mat(k,1349) = .250_r8*rxt(k,460)*y(k,228)
         mat(k,1380) = .250_r8*rxt(k,465)*y(k,228)
         mat(k,1305) = .200_r8*rxt(k,474)*y(k,228)
         mat(k,1409) = .250_r8*rxt(k,427)*y(k,228)
         mat(k,1164) = .250_r8*rxt(k,539)*y(k,228)
         mat(k,1272) = .500_r8*rxt(k,415)*y(k,228)
         mat(k,1214) = .250_r8*rxt(k,544)*y(k,228)
         mat(k,1087) = .250_r8*rxt(k,549)*y(k,228)
         mat(k,1285) = .300_r8*rxt(k,485)*y(k,228)
         mat(k,476) = -(rxt(k,357)*y(k,248))
         mat(k,2062) = -rxt(k,357)*y(k,63)
         mat(k,1606) = rxt(k,354)*y(k,233)
         mat(k,2260) = rxt(k,354)*y(k,228)
         mat(k,1567) = -(rxt(k,205)*y(k,70) + rxt(k,306)*y(k,89) + rxt(k,358)*y(k,248) &
                      + (rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,247))
         mat(k,1862) = -rxt(k,205)*y(k,64)
         mat(k,949) = -rxt(k,306)*y(k,64)
         mat(k,2155) = -rxt(k,358)*y(k,64)
         mat(k,1973) = -(rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,64)
         mat(k,1135) = .100_r8*rxt(k,401)*y(k,163)
         mat(k,1929) = .100_r8*rxt(k,401)*y(k,33)
         mat(k,116) = -(rxt(k,328)*y(k,248))
         mat(k,2013) = -rxt(k,328)*y(k,65)
         mat(k,500) = -(rxt(k,263)*y(k,247) + rxt(k,329)*y(k,70) + rxt(k,330)*y(k,248))
         mat(k,1967) = -rxt(k,263)*y(k,66)
         mat(k,1836) = -rxt(k,329)*y(k,66)
         mat(k,2066) = -rxt(k,330)*y(k,66)
         mat(k,120) = -(rxt(k,331)*y(k,248))
         mat(k,2014) = -rxt(k,331)*y(k,67)
         mat(k,1113) = -((rxt(k,332) + rxt(k,333)) * y(k,228) + (rxt(k,334) + rxt(k,335) &
                      ) * y(k,233) + rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155))
         mat(k,1618) = -(rxt(k,332) + rxt(k,333)) * y(k,68)
         mat(k,2303) = -(rxt(k,334) + rxt(k,335)) * y(k,68)
         mat(k,2632) = -rxt(k,336)*y(k,68)
         mat(k,1776) = -rxt(k,337)*y(k,68)
         mat(k,343) = rxt(k,319)*y(k,70) + rxt(k,320)*y(k,248)
         mat(k,1853) = rxt(k,319)*y(k,47)
         mat(k,2128) = rxt(k,320)*y(k,47)
         mat(k,396) = -(rxt(k,338)*y(k,70) + rxt(k,339)*y(k,248))
         mat(k,1833) = -rxt(k,338)*y(k,69)
         mat(k,2053) = -rxt(k,339)*y(k,69)
         mat(k,1869) = -(rxt(k,204)*y(k,51) + rxt(k,205)*y(k,64) + rxt(k,206)*y(k,93) &
                      + rxt(k,207)*y(k,95) + (rxt(k,208) + rxt(k,209)) * y(k,233) &
                      + rxt(k,210)*y(k,154) + rxt(k,212)*y(k,163) + rxt(k,219)*y(k,75) &
                      + rxt(k,228)*y(k,109) + rxt(k,254)*y(k,22) + rxt(k,312)*y(k,26) &
                      + rxt(k,314)*y(k,29) + rxt(k,316)*y(k,45) + rxt(k,319)*y(k,47) &
                      + rxt(k,321)*y(k,52) + rxt(k,324)*y(k,55) + rxt(k,326)*y(k,61) &
                      + rxt(k,329)*y(k,66) + rxt(k,379)*y(k,32) + rxt(k,409)*y(k,35) &
                      + (rxt(k,557) + rxt(k,558)) * y(k,83))
         mat(k,1692) = -rxt(k,204)*y(k,70)
         mat(k,1572) = -rxt(k,205)*y(k,70)
         mat(k,1539) = -rxt(k,206)*y(k,70)
         mat(k,716) = -rxt(k,207)*y(k,70)
         mat(k,2334) = -(rxt(k,208) + rxt(k,209)) * y(k,70)
         mat(k,2561) = -rxt(k,210)*y(k,70)
         mat(k,1935) = -rxt(k,212)*y(k,70)
         mat(k,1013) = -rxt(k,219)*y(k,70)
         mat(k,1719) = -rxt(k,228)*y(k,70)
         mat(k,909) = -rxt(k,254)*y(k,70)
         mat(k,218) = -rxt(k,312)*y(k,70)
         mat(k,288) = -rxt(k,314)*y(k,70)
         mat(k,554) = -rxt(k,316)*y(k,70)
         mat(k,345) = -rxt(k,319)*y(k,70)
         mat(k,668) = -rxt(k,321)*y(k,70)
         mat(k,442) = -rxt(k,324)*y(k,70)
         mat(k,433) = -rxt(k,326)*y(k,70)
         mat(k,502) = -rxt(k,329)*y(k,70)
         mat(k,330) = -rxt(k,379)*y(k,70)
         mat(k,336) = -rxt(k,409)*y(k,70)
         mat(k,1030) = -(rxt(k,557) + rxt(k,558)) * y(k,70)
         mat(k,2453) = rxt(k,249)*y(k,74)
         mat(k,218) = mat(k,218) + 5.000_r8*rxt(k,312)*y(k,70) + 3.060_r8*rxt(k,313) &
                      *y(k,248)
         mat(k,288) = mat(k,288) + 2.000_r8*rxt(k,314)*y(k,70) + 2.000_r8*rxt(k,315) &
                      *y(k,248)
         mat(k,114) = 4.000_r8*rxt(k,231)*y(k,247)
         mat(k,169) = rxt(k,232)*y(k,247)
         mat(k,134) = 2.000_r8*rxt(k,233)*y(k,247)
         mat(k,180) = 2.000_r8*rxt(k,234)*y(k,247)
         mat(k,138) = 2.000_r8*rxt(k,235)*y(k,247)
         mat(k,185) = rxt(k,236)*y(k,247)
         mat(k,142) = 2.000_r8*rxt(k,237)*y(k,247)
         mat(k,145) = rxt(k,318)*y(k,248)
         mat(k,149) = 3.000_r8*rxt(k,323)*y(k,248)
         mat(k,442) = mat(k,442) + rxt(k,325)*y(k,248)
         mat(k,117) = rxt(k,328)*y(k,248)
         mat(k,121) = 2.000_r8*rxt(k,331)*y(k,248)
         mat(k,1119) = 2.000_r8*rxt(k,336)*y(k,153) + 2.000_r8*rxt(k,337)*y(k,155) &
                      + 2.000_r8*rxt(k,332)*y(k,228) + rxt(k,335)*y(k,233)
         mat(k,400) = rxt(k,339)*y(k,248)
         mat(k,1869) = mat(k,1869) + 5.000_r8*rxt(k,312)*y(k,26) + 2.000_r8*rxt(k,314) &
                      *y(k,29)
         mat(k,2363) = rxt(k,249)*y(k,21) + (4.000_r8*rxt(k,214)+2.000_r8*rxt(k,216)) &
                      *y(k,74) + rxt(k,286)*y(k,125) + rxt(k,218)*y(k,153) &
                      + rxt(k,223)*y(k,162) + rxt(k,568)*y(k,180) + rxt(k,213) &
                      *y(k,228) + rxt(k,224)*y(k,248)
         mat(k,262) = rxt(k,311)*y(k,247)
         mat(k,257) = rxt(k,345)*y(k,247) + rxt(k,340)*y(k,248)
         mat(k,266) = rxt(k,346)*y(k,247) + rxt(k,341)*y(k,248)
         mat(k,317) = rxt(k,347)*y(k,247) + rxt(k,342)*y(k,248)
         mat(k,1742) = rxt(k,226)*y(k,162) + rxt(k,238)*y(k,247) + rxt(k,227)*y(k,248)
         mat(k,2219) = rxt(k,286)*y(k,74)
         mat(k,2662) = 2.000_r8*rxt(k,336)*y(k,68) + rxt(k,218)*y(k,74)
         mat(k,1809) = 2.000_r8*rxt(k,337)*y(k,68)
         mat(k,2426) = rxt(k,223)*y(k,74) + rxt(k,226)*y(k,101)
         mat(k,1550) = rxt(k,568)*y(k,74)
         mat(k,1644) = 2.000_r8*rxt(k,332)*y(k,68) + rxt(k,213)*y(k,74)
         mat(k,2334) = mat(k,2334) + rxt(k,335)*y(k,68)
         mat(k,1981) = 4.000_r8*rxt(k,231)*y(k,37) + rxt(k,232)*y(k,38) &
                      + 2.000_r8*rxt(k,233)*y(k,40) + 2.000_r8*rxt(k,234)*y(k,41) &
                      + 2.000_r8*rxt(k,235)*y(k,42) + rxt(k,236)*y(k,43) &
                      + 2.000_r8*rxt(k,237)*y(k,44) + rxt(k,311)*y(k,81) + rxt(k,345) &
                      *y(k,98) + rxt(k,346)*y(k,99) + rxt(k,347)*y(k,100) + rxt(k,238) &
                      *y(k,101)
         mat(k,2163) = 3.060_r8*rxt(k,313)*y(k,26) + 2.000_r8*rxt(k,315)*y(k,29) &
                      + rxt(k,318)*y(k,46) + 3.000_r8*rxt(k,323)*y(k,53) + rxt(k,325) &
                      *y(k,55) + rxt(k,328)*y(k,65) + 2.000_r8*rxt(k,331)*y(k,67) &
                      + rxt(k,339)*y(k,69) + rxt(k,224)*y(k,74) + rxt(k,340)*y(k,98) &
                      + rxt(k,341)*y(k,99) + rxt(k,342)*y(k,100) + rxt(k,227)*y(k,101)
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
         mat(k,1825) = rxt(k,219)*y(k,75)
         mat(k,2351) = 2.000_r8*rxt(k,215)*y(k,74)
         mat(k,1008) = rxt(k,219)*y(k,70) + (rxt(k,666)+rxt(k,675)+rxt(k,684)) &
                      *y(k,101)
         mat(k,1732) = (rxt(k,666)+rxt(k,675)+rxt(k,684))*y(k,75) + (rxt(k,585) &
                       +rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,109)
         mat(k,1708) = (rxt(k,585)+rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,101)
         mat(k,2350) = 2.000_r8*rxt(k,240)*y(k,74)
         mat(k,620) = -(rxt(k,211)*y(k,248))
         mat(k,2081) = -rxt(k,211)*y(k,73)
         mat(k,1839) = rxt(k,210)*y(k,154)
         mat(k,1735) = rxt(k,601)*y(k,143)
         mat(k,404) = rxt(k,601)*y(k,101)
         mat(k,2533) = rxt(k,210)*y(k,70)
         mat(k,2370) = -(rxt(k,213)*y(k,228) + (4._r8*rxt(k,214) + 4._r8*rxt(k,215) &
                      + 4._r8*rxt(k,216) + 4._r8*rxt(k,240)) * y(k,74) + rxt(k,217) &
                      *y(k,233) + rxt(k,218)*y(k,153) + rxt(k,220)*y(k,154) + rxt(k,223) &
                      *y(k,162) + (rxt(k,224) + rxt(k,225)) * y(k,248) + (rxt(k,248) &
                      + rxt(k,249) + rxt(k,250)) * y(k,21) + (rxt(k,285) + rxt(k,286) &
                      + rxt(k,287)) * y(k,125) + rxt(k,568)*y(k,180))
         mat(k,1650) = -rxt(k,213)*y(k,74)
         mat(k,2341) = -rxt(k,217)*y(k,74)
         mat(k,2669) = -rxt(k,218)*y(k,74)
         mat(k,2568) = -rxt(k,220)*y(k,74)
         mat(k,2433) = -rxt(k,223)*y(k,74)
         mat(k,2170) = -(rxt(k,224) + rxt(k,225)) * y(k,74)
         mat(k,2460) = -(rxt(k,248) + rxt(k,249) + rxt(k,250)) * y(k,74)
         mat(k,2226) = -(rxt(k,285) + rxt(k,286) + rxt(k,287)) * y(k,74)
         mat(k,1554) = -rxt(k,568)*y(k,74)
         mat(k,1876) = rxt(k,228)*y(k,109) + rxt(k,212)*y(k,163) + rxt(k,209)*y(k,233)
         mat(k,1016) = rxt(k,221)*y(k,162)
         mat(k,1747) = rxt(k,239)*y(k,247)
         mat(k,1724) = rxt(k,228)*y(k,70) + rxt(k,229)*y(k,162) + rxt(k,230)*y(k,248)
         mat(k,2433) = mat(k,2433) + rxt(k,221)*y(k,75) + rxt(k,229)*y(k,109)
         mat(k,1942) = rxt(k,212)*y(k,70)
         mat(k,536) = rxt(k,573)*y(k,180)
         mat(k,1554) = mat(k,1554) + rxt(k,573)*y(k,165)
         mat(k,2341) = mat(k,2341) + rxt(k,209)*y(k,70)
         mat(k,1988) = rxt(k,239)*y(k,101)
         mat(k,2170) = mat(k,2170) + rxt(k,230)*y(k,109)
         mat(k,1009) = -(rxt(k,219)*y(k,70) + rxt(k,221)*y(k,162) + rxt(k,222) &
                      *y(k,248) + (rxt(k,666) + rxt(k,675) + rxt(k,684)) * y(k,101))
         mat(k,1847) = -rxt(k,219)*y(k,75)
         mat(k,2411) = -rxt(k,221)*y(k,75)
         mat(k,2119) = -rxt(k,222)*y(k,75)
         mat(k,1736) = -(rxt(k,666) + rxt(k,675) + rxt(k,684)) * y(k,75)
         mat(k,2355) = rxt(k,220)*y(k,154)
         mat(k,2543) = rxt(k,220)*y(k,74)
         mat(k,1200) = -(rxt(k,368)*y(k,248))
         mat(k,2134) = -rxt(k,368)*y(k,77)
         mat(k,1071) = .230_r8*rxt(k,534)*y(k,163)
         mat(k,2473) = rxt(k,243)*y(k,51)
         mat(k,324) = .350_r8*rxt(k,370)*y(k,248)
         mat(k,644) = .630_r8*rxt(k,372)*y(k,163)
         mat(k,1131) = .560_r8*rxt(k,401)*y(k,163)
         mat(k,1681) = rxt(k,243)*y(k,17) + rxt(k,204)*y(k,70) + rxt(k,349)*y(k,155) &
                      + rxt(k,350)*y(k,162) + rxt(k,351)*y(k,248)
         mat(k,439) = rxt(k,324)*y(k,70)
         mat(k,1338) = rxt(k,407)*y(k,155) + rxt(k,408)*y(k,248)
         mat(k,1114) = rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155) + (rxt(k,332) &
                       +rxt(k,333))*y(k,228) + rxt(k,335)*y(k,233)
         mat(k,1855) = rxt(k,204)*y(k,51) + rxt(k,324)*y(k,55)
         mat(k,1048) = rxt(k,395)*y(k,248)
         mat(k,922) = .620_r8*rxt(k,479)*y(k,163)
         mat(k,1326) = .650_r8*rxt(k,432)*y(k,163)
         mat(k,983) = .230_r8*rxt(k,537)*y(k,163)
         mat(k,1433) = .560_r8*rxt(k,446)*y(k,163)
         mat(k,2637) = rxt(k,336)*y(k,68) + .170_r8*rxt(k,505)*y(k,229) &
                      + .220_r8*rxt(k,430)*y(k,240) + .400_r8*rxt(k,508)*y(k,241) &
                      + .350_r8*rxt(k,511)*y(k,243) + .225_r8*rxt(k,546)*y(k,252) &
                      + .250_r8*rxt(k,487)*y(k,256)
         mat(k,1782) = rxt(k,349)*y(k,51) + rxt(k,407)*y(k,58) + rxt(k,337)*y(k,68) &
                      + .220_r8*rxt(k,429)*y(k,240) + .500_r8*rxt(k,488)*y(k,256)
         mat(k,2413) = rxt(k,350)*y(k,51) + rxt(k,562)*y(k,166)
         mat(k,1914) = .230_r8*rxt(k,534)*y(k,6) + .630_r8*rxt(k,372)*y(k,28) &
                      + .560_r8*rxt(k,401)*y(k,33) + .620_r8*rxt(k,479)*y(k,127) &
                      + .650_r8*rxt(k,432)*y(k,134) + .230_r8*rxt(k,537)*y(k,139) &
                      + .560_r8*rxt(k,446)*y(k,140)
         mat(k,422) = rxt(k,562)*y(k,162) + rxt(k,563)*y(k,248)
         mat(k,1181) = .700_r8*rxt(k,555)*y(k,248)
         mat(k,1479) = .220_r8*rxt(k,426)*y(k,240) + .250_r8*rxt(k,484)*y(k,256)
         mat(k,1623) = (rxt(k,332)+rxt(k,333))*y(k,68) + .110_r8*rxt(k,427)*y(k,240) &
                      + .125_r8*rxt(k,544)*y(k,252) + .200_r8*rxt(k,485)*y(k,256)
         mat(k,851) = .170_r8*rxt(k,505)*y(k,153) + .070_r8*rxt(k,504)*y(k,233)
         mat(k,2307) = rxt(k,335)*y(k,68) + .070_r8*rxt(k,504)*y(k,229) &
                      + .160_r8*rxt(k,507)*y(k,241) + .140_r8*rxt(k,510)*y(k,243)
         mat(k,1410) = .220_r8*rxt(k,430)*y(k,153) + .220_r8*rxt(k,429)*y(k,155) &
                      + .220_r8*rxt(k,426)*y(k,227) + .110_r8*rxt(k,427)*y(k,228)
         mat(k,814) = .400_r8*rxt(k,508)*y(k,153) + .160_r8*rxt(k,507)*y(k,233)
         mat(k,958) = .350_r8*rxt(k,511)*y(k,153) + .140_r8*rxt(k,510)*y(k,233)
         mat(k,2134) = mat(k,2134) + .350_r8*rxt(k,370)*y(k,27) + rxt(k,351)*y(k,51) &
                      + rxt(k,408)*y(k,58) + rxt(k,395)*y(k,91) + rxt(k,563)*y(k,166) &
                      + .700_r8*rxt(k,555)*y(k,208)
         mat(k,1218) = .225_r8*rxt(k,546)*y(k,153) + .125_r8*rxt(k,544)*y(k,228)
         mat(k,1288) = .250_r8*rxt(k,487)*y(k,153) + .500_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,227) + .200_r8*rxt(k,485)*y(k,228)
         mat(k,1060) = .270_r8*rxt(k,534)*y(k,163)
         mat(k,1126) = .200_r8*rxt(k,401)*y(k,163)
         mat(k,776) = rxt(k,388)*y(k,248)
         mat(k,675) = .500_r8*rxt(k,389)*y(k,248)
         mat(k,1199) = rxt(k,368)*y(k,248)
         mat(k,1204) = .800_r8*rxt(k,394)*y(k,248)
         mat(k,1046) = rxt(k,395)*y(k,248)
         mat(k,1019) = rxt(k,360)*y(k,248)
         mat(k,682) = .500_r8*rxt(k,445)*y(k,248)
         mat(k,973) = .270_r8*rxt(k,537)*y(k,163)
         mat(k,1429) = .100_r8*rxt(k,446)*y(k,163)
         mat(k,2621) = rxt(k,387)*y(k,227) + .900_r8*rxt(k,546)*y(k,252)
         mat(k,1899) = .270_r8*rxt(k,534)*y(k,6) + .200_r8*rxt(k,401)*y(k,33) &
                      + .270_r8*rxt(k,537)*y(k,139) + .100_r8*rxt(k,446)*y(k,140)
         mat(k,1178) = 1.800_r8*rxt(k,555)*y(k,248)
         mat(k,1476) = rxt(k,387)*y(k,153) + 4.000_r8*rxt(k,384)*y(k,227) &
                      + .900_r8*rxt(k,385)*y(k,228) + rxt(k,459)*y(k,235) &
                      + 2.000_r8*rxt(k,435)*y(k,242) + rxt(k,484)*y(k,256)
         mat(k,1610) = .900_r8*rxt(k,385)*y(k,227) + rxt(k,436)*y(k,242) &
                      + .500_r8*rxt(k,544)*y(k,252)
         mat(k,2292) = .450_r8*rxt(k,437)*y(k,242)
         mat(k,1350) = rxt(k,459)*y(k,227)
         mat(k,1455) = 2.000_r8*rxt(k,435)*y(k,227) + rxt(k,436)*y(k,228) &
                      + .450_r8*rxt(k,437)*y(k,233) + 4.000_r8*rxt(k,438)*y(k,242)
         mat(k,2110) = rxt(k,388)*y(k,59) + .500_r8*rxt(k,389)*y(k,60) + rxt(k,368) &
                      *y(k,77) + .800_r8*rxt(k,394)*y(k,90) + rxt(k,395)*y(k,91) &
                      + rxt(k,360)*y(k,103) + .500_r8*rxt(k,445)*y(k,138) &
                      + 1.800_r8*rxt(k,555)*y(k,208)
         mat(k,1215) = .900_r8*rxt(k,546)*y(k,153) + .500_r8*rxt(k,544)*y(k,228)
         mat(k,1286) = rxt(k,484)*y(k,227)
         mat(k,217) = .470_r8*rxt(k,313)*y(k,248)
         mat(k,1112) = rxt(k,333)*y(k,228) + rxt(k,334)*y(k,233)
         mat(k,395) = rxt(k,338)*y(k,70) + rxt(k,339)*y(k,248)
         mat(k,1832) = rxt(k,338)*y(k,69)
         mat(k,1604) = rxt(k,333)*y(k,68)
         mat(k,2256) = rxt(k,334)*y(k,68)
         mat(k,2052) = .470_r8*rxt(k,313)*y(k,26) + rxt(k,339)*y(k,69)
         mat(k,269) = -(rxt(k,310)*y(k,247))
         mat(k,1965) = -rxt(k,310)*y(k,80)
         mat(k,168) = rxt(k,232)*y(k,247)
         mat(k,173) = rxt(k,262)*y(k,247)
         mat(k,179) = rxt(k,234)*y(k,247)
         mat(k,137) = 2.000_r8*rxt(k,235)*y(k,247)
         mat(k,183) = 2.000_r8*rxt(k,236)*y(k,247)
         mat(k,141) = rxt(k,237)*y(k,247)
         mat(k,125) = 2.000_r8*rxt(k,264)*y(k,247)
         mat(k,265) = rxt(k,346)*y(k,247) + rxt(k,341)*y(k,248)
         mat(k,314) = rxt(k,347)*y(k,247) + rxt(k,342)*y(k,248)
         mat(k,1965) = mat(k,1965) + rxt(k,232)*y(k,38) + rxt(k,262)*y(k,39) &
                      + rxt(k,234)*y(k,41) + 2.000_r8*rxt(k,235)*y(k,42) &
                      + 2.000_r8*rxt(k,236)*y(k,43) + rxt(k,237)*y(k,44) &
                      + 2.000_r8*rxt(k,264)*y(k,94) + rxt(k,346)*y(k,99) + rxt(k,347) &
                      *y(k,100)
         mat(k,2032) = rxt(k,341)*y(k,99) + rxt(k,342)*y(k,100)
         mat(k,260) = -(rxt(k,311)*y(k,247))
         mat(k,1963) = -rxt(k,311)*y(k,81)
         mat(k,133) = rxt(k,233)*y(k,247)
         mat(k,178) = rxt(k,234)*y(k,247)
         mat(k,256) = rxt(k,345)*y(k,247) + rxt(k,340)*y(k,248)
         mat(k,1963) = mat(k,1963) + rxt(k,233)*y(k,40) + rxt(k,234)*y(k,41) &
                      + rxt(k,345)*y(k,98)
         mat(k,2030) = rxt(k,340)*y(k,98)
         mat(k,228) = -(rxt(k,503)*y(k,248))
         mat(k,2024) = -rxt(k,503)*y(k,82)
         mat(k,222) = .180_r8*rxt(k,523)*y(k,248)
         mat(k,2024) = mat(k,2024) + .180_r8*rxt(k,523)*y(k,210)
         mat(k,1026) = -(rxt(k,556)*y(k,21) + (rxt(k,557) + rxt(k,558)) * y(k,70) &
                      + rxt(k,559)*y(k,125) + rxt(k,560)*y(k,155) + (rxt(k,561) &
                      + rxt(k,575)) * y(k,248))
         mat(k,2446) = -rxt(k,556)*y(k,83)
         mat(k,1849) = -(rxt(k,557) + rxt(k,558)) * y(k,83)
         mat(k,2210) = -rxt(k,559)*y(k,83)
         mat(k,1769) = -rxt(k,560)*y(k,83)
         mat(k,2121) = -(rxt(k,561) + rxt(k,575)) * y(k,83)
         mat(k,795) = rxt(k,390)*y(k,233)
         mat(k,2247) = rxt(k,390)*y(k,232)
         mat(k,947) = -(rxt(k,306)*y(k,64) + rxt(k,307)*y(k,93) + rxt(k,308)*y(k,260) &
                      + rxt(k,309)*y(k,106))
         mat(k,1564) = -rxt(k,306)*y(k,89)
         mat(k,1533) = -rxt(k,307)*y(k,89)
         mat(k,2681) = -rxt(k,308)*y(k,89)
         mat(k,2180) = -rxt(k,309)*y(k,89)
         mat(k,174) = rxt(k,262)*y(k,247)
         mat(k,184) = rxt(k,236)*y(k,247)
         mat(k,270) = 2.000_r8*rxt(k,310)*y(k,247)
         mat(k,261) = rxt(k,311)*y(k,247)
         mat(k,1970) = rxt(k,262)*y(k,39) + rxt(k,236)*y(k,43) + 2.000_r8*rxt(k,310) &
                      *y(k,80) + rxt(k,311)*y(k,81)
         mat(k,1207) = -(rxt(k,394)*y(k,248))
         mat(k,2135) = -rxt(k,394)*y(k,90)
         mat(k,691) = .700_r8*rxt(k,470)*y(k,248)
         mat(k,659) = .500_r8*rxt(k,471)*y(k,248)
         mat(k,460) = rxt(k,482)*y(k,248)
         mat(k,2638) = .050_r8*rxt(k,468)*y(k,236) + .530_r8*rxt(k,430)*y(k,240) &
                      + .225_r8*rxt(k,546)*y(k,252) + .250_r8*rxt(k,487)*y(k,256)
         mat(k,1783) = .050_r8*rxt(k,469)*y(k,236) + .530_r8*rxt(k,429)*y(k,240) &
                      + .250_r8*rxt(k,488)*y(k,256)
         mat(k,1480) = .530_r8*rxt(k,426)*y(k,240) + .250_r8*rxt(k,484)*y(k,256)
         mat(k,1624) = .260_r8*rxt(k,427)*y(k,240) + .125_r8*rxt(k,544)*y(k,252) &
                      + .100_r8*rxt(k,485)*y(k,256)
         mat(k,1385) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,1411) = .530_r8*rxt(k,430)*y(k,153) + .530_r8*rxt(k,429)*y(k,155) &
                      + .530_r8*rxt(k,426)*y(k,227) + .260_r8*rxt(k,427)*y(k,228)
         mat(k,2135) = mat(k,2135) + .700_r8*rxt(k,470)*y(k,128) + .500_r8*rxt(k,471) &
                      *y(k,129) + rxt(k,482)*y(k,144)
         mat(k,1219) = .225_r8*rxt(k,546)*y(k,153) + .125_r8*rxt(k,544)*y(k,228)
         mat(k,1289) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,227) + .100_r8*rxt(k,485)*y(k,228)
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
         mat(k,1047) = -(rxt(k,395)*y(k,248))
         mat(k,2123) = -rxt(k,395)*y(k,91)
         mat(k,323) = .650_r8*rxt(k,370)*y(k,248)
         mat(k,1205) = .200_r8*rxt(k,394)*y(k,248)
         mat(k,1149) = rxt(k,483)*y(k,248)
         mat(k,2628) = rxt(k,494)*y(k,221) + .050_r8*rxt(k,468)*y(k,236) &
                      + .400_r8*rxt(k,508)*y(k,241) + .170_r8*rxt(k,511)*y(k,243) &
                      + .700_r8*rxt(k,514)*y(k,249) + .600_r8*rxt(k,521)*y(k,254) &
                      + .250_r8*rxt(k,487)*y(k,256) + .340_r8*rxt(k,527)*y(k,257) &
                      + .170_r8*rxt(k,530)*y(k,259)
         mat(k,1771) = .050_r8*rxt(k,469)*y(k,236) + .250_r8*rxt(k,488)*y(k,256)
         mat(k,570) = rxt(k,494)*y(k,153)
         mat(k,1477) = .250_r8*rxt(k,484)*y(k,256)
         mat(k,1614) = .100_r8*rxt(k,485)*y(k,256)
         mat(k,2299) = .160_r8*rxt(k,507)*y(k,241) + .070_r8*rxt(k,510)*y(k,243)
         mat(k,1383) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,813) = .400_r8*rxt(k,508)*y(k,153) + .160_r8*rxt(k,507)*y(k,233)
         mat(k,957) = .170_r8*rxt(k,511)*y(k,153) + .070_r8*rxt(k,510)*y(k,233)
         mat(k,2123) = mat(k,2123) + .650_r8*rxt(k,370)*y(k,27) + .200_r8*rxt(k,394) &
                      *y(k,90) + rxt(k,483)*y(k,145)
         mat(k,528) = .700_r8*rxt(k,514)*y(k,153)
         mat(k,826) = .600_r8*rxt(k,521)*y(k,153)
         mat(k,1287) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,227) + .100_r8*rxt(k,485)*y(k,228)
         mat(k,842) = .340_r8*rxt(k,527)*y(k,153)
         mat(k,585) = .170_r8*rxt(k,530)*y(k,153)
         mat(k,2516) = -((rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,233) + rxt(k,170) &
                      *y(k,163))
         mat(k,2346) = -(rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,92)
         mat(k,1947) = -rxt(k,170)*y(k,92)
         mat(k,1704) = rxt(k,351)*y(k,248)
         mat(k,1578) = rxt(k,365)*y(k,247)
         mat(k,1881) = rxt(k,206)*y(k,93)
         mat(k,954) = rxt(k,307)*y(k,93)
         mat(k,1544) = rxt(k,206)*y(k,70) + rxt(k,307)*y(k,89) + rxt(k,162)*y(k,162) &
                      + rxt(k,153)*y(k,247) + rxt(k,171)*y(k,248)
         mat(k,1518) = rxt(k,266)*y(k,247)
         mat(k,1752) = rxt(k,239)*y(k,247)
         mat(k,580) = rxt(k,192)*y(k,248)
         mat(k,2438) = rxt(k,162)*y(k,93) + rxt(k,174)*y(k,248)
         mat(k,426) = rxt(k,563)*y(k,248)
         mat(k,608) = rxt(k,569)*y(k,248)
         mat(k,1558) = rxt(k,574)*y(k,248)
         mat(k,1993) = rxt(k,365)*y(k,64) + rxt(k,153)*y(k,93) + rxt(k,266)*y(k,97) &
                      + rxt(k,239)*y(k,101)
         mat(k,2175) = rxt(k,351)*y(k,51) + rxt(k,171)*y(k,93) + rxt(k,192)*y(k,141) &
                      + rxt(k,174)*y(k,162) + rxt(k,563)*y(k,166) + rxt(k,569) &
                      *y(k,178) + rxt(k,574)*y(k,180)
         mat(k,1534) = -(rxt(k,153)*y(k,247) + rxt(k,162)*y(k,162) + rxt(k,171) &
                      *y(k,248) + rxt(k,206)*y(k,70) + rxt(k,307)*y(k,89))
         mat(k,1972) = -rxt(k,153)*y(k,93)
         mat(k,2416) = -rxt(k,162)*y(k,93)
         mat(k,2153) = -rxt(k,171)*y(k,93)
         mat(k,1860) = -rxt(k,206)*y(k,93)
         mat(k,948) = -rxt(k,307)*y(k,93)
         mat(k,1566) = rxt(k,366)*y(k,247)
         mat(k,2497) = rxt(k,164)*y(k,233)
         mat(k,2325) = rxt(k,164)*y(k,92)
         mat(k,1972) = mat(k,1972) + rxt(k,366)*y(k,64)
         mat(k,124) = -(rxt(k,264)*y(k,247))
         mat(k,1952) = -rxt(k,264)*y(k,94)
         mat(k,714) = -(rxt(k,163)*y(k,162) + rxt(k,172)*y(k,248) + rxt(k,207)*y(k,70))
         mat(k,2409) = -rxt(k,163)*y(k,95)
         mat(k,2093) = -rxt(k,172)*y(k,95)
         mat(k,1842) = -rxt(k,207)*y(k,95)
         mat(k,2278) = 2.000_r8*rxt(k,178)*y(k,233)
         mat(k,2093) = mat(k,2093) + 2.000_r8*rxt(k,177)*y(k,248)
         mat(k,291) = rxt(k,576)*y(k,260)
         mat(k,2678) = rxt(k,576)*y(k,182)
         mat(k,1507) = -(rxt(k,259)*y(k,162) + rxt(k,260)*y(k,248) + (rxt(k,265) &
                      + rxt(k,266)) * y(k,247) + (rxt(k,583) + rxt(k,660) + rxt(k,668) &
                      + rxt(k,677)) * y(k,109) + (rxt(k,584) + rxt(k,658) + rxt(k,671) &
                      + rxt(k,680)) * y(k,108) + (rxt(k,591) + rxt(k,687) + rxt(k,691) &
                      + rxt(k,695)) * y(k,110))
         mat(k,2414) = -rxt(k,259)*y(k,97)
         mat(k,2151) = -rxt(k,260)*y(k,97)
         mat(k,1971) = -(rxt(k,265) + rxt(k,266)) * y(k,97)
         mat(k,1712) = -(rxt(k,583) + rxt(k,660) + rxt(k,668) + rxt(k,677)) * y(k,97)
         mat(k,1660) = -(rxt(k,584) + rxt(k,658) + rxt(k,671) + rxt(k,680)) * y(k,97)
         mat(k,1583) = -(rxt(k,591) + rxt(k,687) + rxt(k,691) + rxt(k,695)) * y(k,97)
         mat(k,2474) = rxt(k,243)*y(k,51) + rxt(k,244)*y(k,233)
         mat(k,1682) = rxt(k,243)*y(k,17)
         mat(k,2323) = rxt(k,244)*y(k,17)
         mat(k,255) = -(rxt(k,340)*y(k,248) + rxt(k,345)*y(k,247))
         mat(k,2029) = -rxt(k,340)*y(k,98)
         mat(k,1962) = -rxt(k,345)*y(k,98)
         mat(k,264) = -(rxt(k,341)*y(k,248) + rxt(k,346)*y(k,247))
         mat(k,2031) = -rxt(k,341)*y(k,99)
         mat(k,1964) = -rxt(k,346)*y(k,99)
         mat(k,315) = -(rxt(k,342)*y(k,248) + rxt(k,347)*y(k,247))
         mat(k,2041) = -rxt(k,342)*y(k,100)
         mat(k,1966) = -rxt(k,347)*y(k,100)
         mat(k,1740) = -(rxt(k,226)*y(k,162) + rxt(k,227)*y(k,248) + (rxt(k,238) &
                      + rxt(k,239)) * y(k,247) + (rxt(k,585) + rxt(k,656) + rxt(k,667) &
                      + rxt(k,676)) * y(k,109) + (rxt(k,586) + rxt(k,657) + rxt(k,670) &
                      + rxt(k,679)) * y(k,108) + (rxt(k,590) + rxt(k,686) + rxt(k,690) &
                      + rxt(k,694)) * y(k,110) + rxt(k,601)*y(k,143) + (rxt(k,666) &
                      + rxt(k,675) + rxt(k,684)) * y(k,75))
         mat(k,2424) = -rxt(k,226)*y(k,101)
         mat(k,2161) = -rxt(k,227)*y(k,101)
         mat(k,1979) = -(rxt(k,238) + rxt(k,239)) * y(k,101)
         mat(k,1717) = -(rxt(k,585) + rxt(k,656) + rxt(k,667) + rxt(k,676)) * y(k,101)
         mat(k,1665) = -(rxt(k,586) + rxt(k,657) + rxt(k,670) + rxt(k,679)) * y(k,101)
         mat(k,1588) = -(rxt(k,590) + rxt(k,686) + rxt(k,690) + rxt(k,694)) * y(k,101)
         mat(k,405) = -rxt(k,601)*y(k,101)
         mat(k,1011) = -(rxt(k,666) + rxt(k,675) + rxt(k,684)) * y(k,101)
         mat(k,287) = rxt(k,314)*y(k,70)
         mat(k,329) = rxt(k,379)*y(k,70)
         mat(k,335) = rxt(k,409)*y(k,70)
         mat(k,553) = rxt(k,316)*y(k,70)
         mat(k,344) = rxt(k,319)*y(k,70)
         mat(k,1690) = rxt(k,204)*y(k,70)
         mat(k,667) = rxt(k,321)*y(k,70)
         mat(k,441) = 2.000_r8*rxt(k,324)*y(k,70)
         mat(k,432) = rxt(k,326)*y(k,70)
         mat(k,1570) = rxt(k,205)*y(k,70)
         mat(k,501) = rxt(k,329)*y(k,70)
         mat(k,399) = rxt(k,338)*y(k,70)
         mat(k,1867) = rxt(k,314)*y(k,29) + rxt(k,379)*y(k,32) + rxt(k,409)*y(k,35) &
                      + rxt(k,316)*y(k,45) + rxt(k,319)*y(k,47) + rxt(k,204)*y(k,51) &
                      + rxt(k,321)*y(k,52) + 2.000_r8*rxt(k,324)*y(k,55) + rxt(k,326) &
                      *y(k,61) + rxt(k,205)*y(k,64) + rxt(k,329)*y(k,66) + rxt(k,338) &
                      *y(k,69) + rxt(k,558)*y(k,83) + rxt(k,206)*y(k,93) + rxt(k,207) &
                      *y(k,95) + rxt(k,228)*y(k,109) + rxt(k,208)*y(k,233)
         mat(k,2361) = rxt(k,225)*y(k,248)
         mat(k,1028) = rxt(k,558)*y(k,70)
         mat(k,1537) = rxt(k,206)*y(k,70)
         mat(k,715) = rxt(k,207)*y(k,70)
         mat(k,1717) = mat(k,1717) + rxt(k,228)*y(k,70)
         mat(k,2332) = rxt(k,208)*y(k,70)
         mat(k,2161) = mat(k,2161) + rxt(k,225)*y(k,74)
         mat(k,205) = -(rxt(k,359)*y(k,248) + rxt(k,367)*y(k,247))
         mat(k,2021) = -rxt(k,359)*y(k,102)
         mat(k,1960) = -rxt(k,367)*y(k,102)
         mat(k,1020) = -(rxt(k,360)*y(k,248))
         mat(k,2120) = -rxt(k,360)*y(k,103)
         mat(k,1062) = .050_r8*rxt(k,534)*y(k,163)
         mat(k,322) = .350_r8*rxt(k,370)*y(k,248)
         mat(k,643) = .370_r8*rxt(k,372)*y(k,163)
         mat(k,1128) = .120_r8*rxt(k,401)*y(k,163)
         mat(k,920) = .110_r8*rxt(k,479)*y(k,163)
         mat(k,1325) = .330_r8*rxt(k,432)*y(k,163)
         mat(k,976) = .050_r8*rxt(k,537)*y(k,163)
         mat(k,1430) = .120_r8*rxt(k,446)*y(k,163)
         mat(k,2627) = rxt(k,363)*y(k,234)
         mat(k,1903) = .050_r8*rxt(k,534)*y(k,6) + .370_r8*rxt(k,372)*y(k,28) &
                      + .120_r8*rxt(k,401)*y(k,33) + .110_r8*rxt(k,479)*y(k,127) &
                      + .330_r8*rxt(k,432)*y(k,134) + .050_r8*rxt(k,537)*y(k,139) &
                      + .120_r8*rxt(k,446)*y(k,140)
         mat(k,2298) = rxt(k,361)*y(k,234)
         mat(k,521) = rxt(k,363)*y(k,153) + rxt(k,361)*y(k,233)
         mat(k,2120) = mat(k,2120) + .350_r8*rxt(k,370)*y(k,27)
         mat(k,1562) = rxt(k,306)*y(k,89)
         mat(k,946) = rxt(k,306)*y(k,64) + rxt(k,307)*y(k,93) + rxt(k,309)*y(k,106) &
                      + rxt(k,308)*y(k,260)
         mat(k,1532) = rxt(k,307)*y(k,89)
         mat(k,2179) = rxt(k,309)*y(k,89)
         mat(k,2680) = rxt(k,308)*y(k,89)
         mat(k,1258) = -(rxt(k,267)*y(k,155) + rxt(k,295)*y(k,248) + (rxt(k,587) &
                      + rxt(k,661) + rxt(k,669) + rxt(k,678)) * y(k,109) + (rxt(k,588) &
                      + rxt(k,659) + rxt(k,672) + rxt(k,681)) * y(k,108) + (rxt(k,592) &
                      + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,110))
         mat(k,1787) = -rxt(k,267)*y(k,105)
         mat(k,2139) = -rxt(k,295)*y(k,105)
         mat(k,1711) = -(rxt(k,587) + rxt(k,661) + rxt(k,669) + rxt(k,678)) * y(k,105)
         mat(k,1659) = -(rxt(k,588) + rxt(k,659) + rxt(k,672) + rxt(k,681)) * y(k,105)
         mat(k,1582) = -(rxt(k,592) + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,105)
         mat(k,2384) = rxt(k,273)*y(k,233)
         mat(k,2311) = rxt(k,273)*y(k,115)
         mat(k,2192) = -(rxt(k,201)*y(k,248) + rxt(k,309)*y(k,89))
         mat(k,2167) = -rxt(k,201)*y(k,106)
         mat(k,953) = -rxt(k,309)*y(k,106)
         mat(k,1696) = rxt(k,349)*y(k,155)
         mat(k,1196) = rxt(k,381)*y(k,155)
         mat(k,1343) = rxt(k,407)*y(k,155)
         mat(k,1015) = (rxt(k,666)+rxt(k,675)+rxt(k,684))*y(k,101)
         mat(k,1032) = rxt(k,560)*y(k,155)
         mat(k,1745) = (rxt(k,666)+rxt(k,675)+rxt(k,684))*y(k,75) + rxt(k,601) &
                      *y(k,143)
         mat(k,1266) = rxt(k,267)*y(k,155)
         mat(k,1593) = rxt(k,297)*y(k,155)
         mat(k,407) = rxt(k,601)*y(k,101)
         mat(k,2565) = rxt(k,200)*y(k,248)
         mat(k,1813) = rxt(k,349)*y(k,51) + rxt(k,381)*y(k,54) + rxt(k,407)*y(k,58) &
                      + rxt(k,560)*y(k,83) + rxt(k,267)*y(k,105) + rxt(k,297)*y(k,110)
         mat(k,2167) = mat(k,2167) + rxt(k,200)*y(k,154)
         mat(k,470) = -(rxt(k,179)*y(k,248))
         mat(k,2061) = -rxt(k,179)*y(k,107)
         mat(k,2526) = rxt(k,198)*y(k,233)
         mat(k,2259) = rxt(k,198)*y(k,154)
         mat(k,1663) = -(rxt(k,261)*y(k,162) + (rxt(k,584) + rxt(k,658) + rxt(k,671) &
                      + rxt(k,680)) * y(k,97) + (rxt(k,586) + rxt(k,657) + rxt(k,670) &
                      + rxt(k,679)) * y(k,101) + (rxt(k,588) + rxt(k,659) + rxt(k,672) &
                      + rxt(k,681)) * y(k,105))
         mat(k,2421) = -rxt(k,261)*y(k,108)
         mat(k,1509) = -(rxt(k,584) + rxt(k,658) + rxt(k,671) + rxt(k,680)) * y(k,108)
         mat(k,1738) = -(rxt(k,586) + rxt(k,657) + rxt(k,670) + rxt(k,679)) * y(k,108)
         mat(k,1261) = -(rxt(k,588) + rxt(k,659) + rxt(k,672) + rxt(k,681)) * y(k,108)
         mat(k,542) = rxt(k,242)*y(k,248)
         mat(k,2449) = rxt(k,251)*y(k,233)
         mat(k,2329) = rxt(k,251)*y(k,21)
         mat(k,2158) = rxt(k,242)*y(k,18)
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
         mat(k,1716) = -(rxt(k,228)*y(k,70) + rxt(k,229)*y(k,162) + rxt(k,230) &
                      *y(k,248) + (rxt(k,583) + rxt(k,660) + rxt(k,668) + rxt(k,677) &
                      ) * y(k,97) + (rxt(k,585) + rxt(k,656) + rxt(k,667) + rxt(k,676) &
                      ) * y(k,101) + (rxt(k,587) + rxt(k,661) + rxt(k,669) + rxt(k,678) &
                      ) * y(k,105))
         mat(k,1866) = -rxt(k,228)*y(k,109)
         mat(k,2423) = -rxt(k,229)*y(k,109)
         mat(k,2160) = -rxt(k,230)*y(k,109)
         mat(k,1510) = -(rxt(k,583) + rxt(k,660) + rxt(k,668) + rxt(k,677)) * y(k,109)
         mat(k,1739) = -(rxt(k,585) + rxt(k,656) + rxt(k,667) + rxt(k,676)) * y(k,109)
         mat(k,1262) = -(rxt(k,587) + rxt(k,661) + rxt(k,669) + rxt(k,678)) * y(k,109)
         mat(k,1117) = rxt(k,335)*y(k,233)
         mat(k,621) = rxt(k,211)*y(k,248)
         mat(k,2360) = rxt(k,217)*y(k,233)
         mat(k,1010) = rxt(k,222)*y(k,248)
         mat(k,2331) = rxt(k,335)*y(k,68) + rxt(k,217)*y(k,74)
         mat(k,2160) = mat(k,2160) + rxt(k,211)*y(k,73) + rxt(k,222)*y(k,75)
         mat(k,1585) = -(rxt(k,268)*y(k,248) + rxt(k,297)*y(k,155) + (rxt(k,590) &
                      + rxt(k,686) + rxt(k,690) + rxt(k,694)) * y(k,101) + (rxt(k,591) &
                      + rxt(k,687) + rxt(k,691) + rxt(k,695)) * y(k,97) + (rxt(k,592) &
                      + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,105))
         mat(k,2156) = -rxt(k,268)*y(k,110)
         mat(k,1802) = -rxt(k,297)*y(k,110)
         mat(k,1737) = -(rxt(k,590) + rxt(k,686) + rxt(k,690) + rxt(k,694)) * y(k,110)
         mat(k,1508) = -(rxt(k,591) + rxt(k,687) + rxt(k,691) + rxt(k,695)) * y(k,110)
         mat(k,1260) = -(rxt(k,592) + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,110)
         mat(k,1522) = rxt(k,271)*y(k,248)
         mat(k,2214) = rxt(k,288)*y(k,233)
         mat(k,2327) = rxt(k,288)*y(k,125)
         mat(k,2156) = mat(k,2156) + rxt(k,271)*y(k,116)
         mat(k,1237) = -(rxt(k,425)*y(k,248))
         mat(k,2137) = -rxt(k,425)*y(k,111)
         mat(k,692) = .300_r8*rxt(k,470)*y(k,248)
         mat(k,660) = .500_r8*rxt(k,471)*y(k,248)
         mat(k,2640) = rxt(k,424)*y(k,230) + rxt(k,431)*y(k,240)
         mat(k,652) = rxt(k,424)*y(k,153)
         mat(k,1412) = rxt(k,431)*y(k,153)
         mat(k,2137) = mat(k,2137) + .300_r8*rxt(k,470)*y(k,128) + .500_r8*rxt(k,471) &
                      *y(k,129)
         mat(k,272) = -(rxt(k,456)*y(k,248))
         mat(k,2033) = -rxt(k,456)*y(k,112)
         mat(k,1250) = -(rxt(k,410)*y(k,248))
         mat(k,2138) = -rxt(k,410)*y(k,113)
         mat(k,693) = .700_r8*rxt(k,470)*y(k,248)
         mat(k,661) = .500_r8*rxt(k,471)*y(k,248)
         mat(k,683) = .500_r8*rxt(k,445)*y(k,248)
         mat(k,2641) = .050_r8*rxt(k,468)*y(k,236) + .220_r8*rxt(k,430)*y(k,240) &
                      + .250_r8*rxt(k,487)*y(k,256)
         mat(k,1786) = .050_r8*rxt(k,469)*y(k,236) + .220_r8*rxt(k,429)*y(k,240) &
                      + .250_r8*rxt(k,488)*y(k,256)
         mat(k,636) = .500_r8*rxt(k,414)*y(k,248)
         mat(k,1481) = .220_r8*rxt(k,426)*y(k,240) + .250_r8*rxt(k,484)*y(k,256)
         mat(k,1626) = .230_r8*rxt(k,427)*y(k,240) + .200_r8*rxt(k,415)*y(k,251) &
                      + .100_r8*rxt(k,485)*y(k,256)
         mat(k,1387) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,1413) = .220_r8*rxt(k,430)*y(k,153) + .220_r8*rxt(k,429)*y(k,155) &
                      + .220_r8*rxt(k,426)*y(k,227) + .230_r8*rxt(k,427)*y(k,228)
         mat(k,2138) = mat(k,2138) + .700_r8*rxt(k,470)*y(k,128) + .500_r8*rxt(k,471) &
                      *y(k,129) + .500_r8*rxt(k,445)*y(k,138) + .500_r8*rxt(k,414) &
                      *y(k,176)
         mat(k,1273) = .200_r8*rxt(k,415)*y(k,228)
         mat(k,1290) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,227) + .100_r8*rxt(k,485)*y(k,228)
         mat(k,351) = -(rxt(k,457)*y(k,248))
         mat(k,2046) = -rxt(k,457)*y(k,114)
         mat(k,2594) = .870_r8*rxt(k,468)*y(k,236)
         mat(k,1757) = .950_r8*rxt(k,469)*y(k,236)
         mat(k,1472) = rxt(k,464)*y(k,236)
         mat(k,1603) = .750_r8*rxt(k,465)*y(k,236)
         mat(k,1376) = .870_r8*rxt(k,468)*y(k,153) + .950_r8*rxt(k,469)*y(k,155) &
                      + rxt(k,464)*y(k,227) + .750_r8*rxt(k,465)*y(k,228)
         mat(k,2399) = -(rxt(k,272)*y(k,21) + rxt(k,273)*y(k,233) + rxt(k,274) &
                      *y(k,126) + rxt(k,276)*y(k,154) + rxt(k,278)*y(k,155) + rxt(k,280) &
                      *y(k,153) + rxt(k,281)*y(k,163))
         mat(k,2461) = -rxt(k,272)*y(k,115)
         mat(k,2342) = -rxt(k,273)*y(k,115)
         mat(k,895) = -rxt(k,274)*y(k,115)
         mat(k,2569) = -rxt(k,276)*y(k,115)
         mat(k,1817) = -rxt(k,278)*y(k,115)
         mat(k,2670) = -rxt(k,280)*y(k,115)
         mat(k,1943) = -rxt(k,281)*y(k,115)
         mat(k,2489) = rxt(k,282)*y(k,125)
         mat(k,2461) = mat(k,2461) + rxt(k,283)*y(k,125)
         mat(k,436) = rxt(k,326)*y(k,70) + rxt(k,327)*y(k,248)
         mat(k,1877) = rxt(k,326)*y(k,61)
         mat(k,2371) = (rxt(k,285)+rxt(k,286))*y(k,125)
         mat(k,1035) = rxt(k,559)*y(k,125)
         mat(k,1267) = rxt(k,267)*y(k,155) + rxt(k,295)*y(k,248)
         mat(k,1528) = rxt(k,269)*y(k,155) + rxt(k,270)*y(k,162) + rxt(k,271)*y(k,248)
         mat(k,2227) = rxt(k,282)*y(k,17) + rxt(k,283)*y(k,21) + (rxt(k,285) &
                       +rxt(k,286))*y(k,74) + rxt(k,559)*y(k,83) + 2.000_r8*rxt(k,301) &
                      *y(k,125) + rxt(k,289)*y(k,153) + rxt(k,292)*y(k,162) &
                      + rxt(k,294)*y(k,248)
         mat(k,2670) = mat(k,2670) + rxt(k,289)*y(k,125)
         mat(k,1817) = mat(k,1817) + rxt(k,267)*y(k,105) + rxt(k,269)*y(k,116)
         mat(k,2434) = rxt(k,270)*y(k,116) + rxt(k,292)*y(k,125)
         mat(k,2171) = rxt(k,327)*y(k,61) + rxt(k,295)*y(k,105) + rxt(k,271)*y(k,116) &
                      + rxt(k,294)*y(k,125)
         mat(k,1521) = -(rxt(k,269)*y(k,155) + rxt(k,270)*y(k,162) + rxt(k,271) &
                      *y(k,248))
         mat(k,1799) = -rxt(k,269)*y(k,116)
         mat(k,2415) = -rxt(k,270)*y(k,116)
         mat(k,2152) = -rxt(k,271)*y(k,116)
         mat(k,1259) = (rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,110)
         mat(k,1584) = (rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,105)
         mat(k,2385) = rxt(k,274)*y(k,126)
         mat(k,210) = 2.000_r8*rxt(k,279)*y(k,123)
         mat(k,311) = 2.000_r8*rxt(k,275)*y(k,124)
         mat(k,889) = rxt(k,274)*y(k,115)
         mat(k,2204) = 2.000_r8*rxt(k,302)*y(k,125)
         mat(k,2205) = rxt(k,304)*y(k,167)
         mat(k,592) = rxt(k,304)*y(k,125)
         mat(k,591) = 2.000_r8*rxt(k,305)*y(k,167)
         mat(k,1504) = (rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,110)
         mat(k,1256) = (rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,108)
         mat(k,1656) = (rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,105)
         mat(k,1580) = (rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,97)
         mat(k,2353) = rxt(k,287)*y(k,125)
         mat(k,1733) = (rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,110)
         mat(k,1257) = (rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,109)
         mat(k,1709) = (rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,105)
         mat(k,1581) = (rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,101)
         mat(k,2207) = rxt(k,287)*y(k,74)
         mat(k,161) = -(rxt(k,458)*y(k,248))
         mat(k,2017) = -rxt(k,458)*y(k,122)
         mat(k,804) = .600_r8*rxt(k,481)*y(k,248)
         mat(k,2017) = mat(k,2017) + .600_r8*rxt(k,481)*y(k,131)
         mat(k,209) = -(4._r8*rxt(k,279)*y(k,123))
         mat(k,2379) = rxt(k,280)*y(k,153)
         mat(k,2589) = rxt(k,280)*y(k,115)
         mat(k,308) = -(4._r8*rxt(k,275)*y(k,124))
         mat(k,2380) = rxt(k,276)*y(k,154)
         mat(k,2522) = rxt(k,276)*y(k,115)
         mat(k,2224) = -(rxt(k,282)*y(k,17) + (rxt(k,283) + rxt(k,284)) * y(k,21) &
                      + (rxt(k,285) + rxt(k,286) + rxt(k,287)) * y(k,74) + rxt(k,288) &
                      *y(k,233) + rxt(k,289)*y(k,153) + rxt(k,290)*y(k,154) + rxt(k,291) &
                      *y(k,155) + rxt(k,292)*y(k,162) + rxt(k,293)*y(k,163) + rxt(k,294) &
                      *y(k,248) + (4._r8*rxt(k,301) + 4._r8*rxt(k,302)) * y(k,125) &
                      + rxt(k,304)*y(k,167) + rxt(k,559)*y(k,83))
         mat(k,2486) = -rxt(k,282)*y(k,125)
         mat(k,2458) = -(rxt(k,283) + rxt(k,284)) * y(k,125)
         mat(k,2368) = -(rxt(k,285) + rxt(k,286) + rxt(k,287)) * y(k,125)
         mat(k,2339) = -rxt(k,288)*y(k,125)
         mat(k,2667) = -rxt(k,289)*y(k,125)
         mat(k,2566) = -rxt(k,290)*y(k,125)
         mat(k,1814) = -rxt(k,291)*y(k,125)
         mat(k,2431) = -rxt(k,292)*y(k,125)
         mat(k,1940) = -rxt(k,293)*y(k,125)
         mat(k,2168) = -rxt(k,294)*y(k,125)
         mat(k,594) = -rxt(k,304)*y(k,125)
         mat(k,1033) = -rxt(k,559)*y(k,125)
         mat(k,2458) = mat(k,2458) + rxt(k,272)*y(k,115)
         mat(k,1594) = rxt(k,297)*y(k,155) + rxt(k,268)*y(k,248)
         mat(k,2396) = rxt(k,272)*y(k,21) + rxt(k,278)*y(k,155) + rxt(k,281)*y(k,163)
         mat(k,1527) = rxt(k,270)*y(k,162)
         mat(k,2667) = mat(k,2667) + rxt(k,296)*y(k,167)
         mat(k,1814) = mat(k,1814) + rxt(k,297)*y(k,110) + rxt(k,278)*y(k,115)
         mat(k,2431) = mat(k,2431) + rxt(k,270)*y(k,116)
         mat(k,1940) = mat(k,1940) + rxt(k,281)*y(k,115)
         mat(k,594) = mat(k,594) + rxt(k,296)*y(k,153)
         mat(k,2168) = mat(k,2168) + rxt(k,268)*y(k,110)
         mat(k,888) = -(rxt(k,274)*y(k,115))
         mat(k,2383) = -rxt(k,274)*y(k,126)
         mat(k,1520) = rxt(k,269)*y(k,155)
         mat(k,2209) = rxt(k,290)*y(k,154)
         mat(k,2540) = rxt(k,290)*y(k,125)
         mat(k,1763) = rxt(k,269)*y(k,116)
         mat(k,919) = -(rxt(k,472)*y(k,155) + rxt(k,479)*y(k,163) + rxt(k,480) &
                      *y(k,248))
         mat(k,1765) = -rxt(k,472)*y(k,127)
         mat(k,1900) = -rxt(k,479)*y(k,127)
         mat(k,2113) = -rxt(k,480)*y(k,127)
         mat(k,690) = -(rxt(k,470)*y(k,248))
         mat(k,2090) = -rxt(k,470)*y(k,128)
         mat(k,2609) = .080_r8*rxt(k,462)*y(k,235)
         mat(k,1347) = .080_r8*rxt(k,462)*y(k,153)
         mat(k,657) = -(rxt(k,471)*y(k,248))
         mat(k,2086) = -rxt(k,471)*y(k,129)
         mat(k,2608) = .080_r8*rxt(k,468)*y(k,236)
         mat(k,1377) = .080_r8*rxt(k,468)*y(k,153)
         mat(k,446) = -(rxt(k,478)*y(k,248))
         mat(k,2058) = -rxt(k,478)*y(k,130)
         mat(k,2257) = rxt(k,475)*y(k,237)
         mat(k,1302) = rxt(k,475)*y(k,233)
         mat(k,805) = -(rxt(k,481)*y(k,248))
         mat(k,2103) = -rxt(k,481)*y(k,131)
         mat(k,2286) = rxt(k,461)*y(k,235) + rxt(k,466)*y(k,236)
         mat(k,1348) = rxt(k,461)*y(k,233)
         mat(k,1379) = rxt(k,466)*y(k,233)
         mat(k,83) = -(rxt(k,641)*y(k,248))
         mat(k,2007) = -rxt(k,641)*y(k,132)
         mat(k,1327) = -(rxt(k,432)*y(k,163) + rxt(k,433)*y(k,248))
         mat(k,1920) = -rxt(k,432)*y(k,134)
         mat(k,2143) = -rxt(k,433)*y(k,134)
         mat(k,924) = .300_r8*rxt(k,479)*y(k,163)
         mat(k,2645) = .360_r8*rxt(k,462)*y(k,235)
         mat(k,1791) = .400_r8*rxt(k,463)*y(k,235)
         mat(k,1920) = mat(k,1920) + .300_r8*rxt(k,479)*y(k,127)
         mat(k,1484) = .390_r8*rxt(k,459)*y(k,235)
         mat(k,1630) = .310_r8*rxt(k,460)*y(k,235)
         mat(k,1357) = .360_r8*rxt(k,462)*y(k,153) + .400_r8*rxt(k,463)*y(k,155) &
                      + .390_r8*rxt(k,459)*y(k,227) + .310_r8*rxt(k,460)*y(k,228)
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
         mat(k,354) = -(rxt(k,434)*y(k,248))
         mat(k,2047) = -rxt(k,434)*y(k,135)
         mat(k,2251) = rxt(k,428)*y(k,240)
         mat(k,1408) = rxt(k,428)*y(k,233)
         mat(k,598) = -(rxt(k,443)*y(k,248))
         mat(k,2078) = -rxt(k,443)*y(k,136)
         mat(k,2605) = .800_r8*rxt(k,452)*y(k,219)
         mat(k,993) = .800_r8*rxt(k,452)*y(k,153)
         mat(k,359) = -(rxt(k,444)*y(k,248))
         mat(k,2048) = -rxt(k,444)*y(k,137)
         mat(k,2252) = .800_r8*rxt(k,441)*y(k,244)
         mat(k,767) = .800_r8*rxt(k,441)*y(k,233)
         mat(k,681) = -(rxt(k,445)*y(k,248))
         mat(k,2089) = -rxt(k,445)*y(k,138)
         mat(k,2535) = rxt(k,448)*y(k,242)
         mat(k,1453) = rxt(k,448)*y(k,154)
         mat(k,974) = -(rxt(k,536)*y(k,155) + rxt(k,537)*y(k,163) + rxt(k,538) &
                      *y(k,248))
         mat(k,1766) = -rxt(k,536)*y(k,139)
         mat(k,1901) = -rxt(k,537)*y(k,139)
         mat(k,2117) = -rxt(k,538)*y(k,139)
         mat(k,1437) = -(rxt(k,446)*y(k,163) + rxt(k,447)*y(k,248))
         mat(k,1925) = -rxt(k,446)*y(k,140)
         mat(k,2148) = -rxt(k,447)*y(k,140)
         mat(k,927) = .200_r8*rxt(k,479)*y(k,163)
         mat(k,2650) = .560_r8*rxt(k,462)*y(k,235)
         mat(k,1796) = .600_r8*rxt(k,463)*y(k,235)
         mat(k,1925) = mat(k,1925) + .200_r8*rxt(k,479)*y(k,127)
         mat(k,1489) = .610_r8*rxt(k,459)*y(k,235)
         mat(k,1635) = .440_r8*rxt(k,460)*y(k,235)
         mat(k,1361) = .560_r8*rxt(k,462)*y(k,153) + .600_r8*rxt(k,463)*y(k,155) &
                      + .610_r8*rxt(k,459)*y(k,227) + .440_r8*rxt(k,460)*y(k,228)
         mat(k,576) = -(rxt(k,180)*y(k,153) + (rxt(k,181) + rxt(k,182) + rxt(k,183) &
                      ) * y(k,154) + rxt(k,192)*y(k,248))
         mat(k,2602) = -rxt(k,180)*y(k,141)
         mat(k,2530) = -(rxt(k,181) + rxt(k,182) + rxt(k,183)) * y(k,141)
         mat(k,2076) = -rxt(k,192)*y(k,141)
         mat(k,213) = -((rxt(k,196) + rxt(k,197)) * y(k,247))
         mat(k,1961) = -(rxt(k,196) + rxt(k,197)) * y(k,142)
         mat(k,575) = rxt(k,181)*y(k,154)
         mat(k,2521) = rxt(k,181)*y(k,141)
         mat(k,2524) = rxt(k,199)*y(k,155)
         mat(k,1758) = rxt(k,199)*y(k,154)
         mat(k,458) = -(rxt(k,482)*y(k,248))
         mat(k,2059) = -rxt(k,482)*y(k,144)
         mat(k,1605) = .200_r8*rxt(k,474)*y(k,237)
         mat(k,1303) = .200_r8*rxt(k,474)*y(k,228)
         mat(k,1150) = -(rxt(k,483)*y(k,248))
         mat(k,2130) = -rxt(k,483)*y(k,145)
         mat(k,2633) = rxt(k,476)*y(k,237)
         mat(k,1778) = rxt(k,477)*y(k,237)
         mat(k,1478) = rxt(k,473)*y(k,237)
         mat(k,1619) = .800_r8*rxt(k,474)*y(k,237)
         mat(k,1307) = rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155) + rxt(k,473)*y(k,227) &
                      + .800_r8*rxt(k,474)*y(k,228)
         mat(k,110) = -(rxt(k,594)*y(k,248))
         mat(k,2012) = -rxt(k,594)*y(k,149)
         mat(k,2676) = -(rxt(k,180)*y(k,141) + rxt(k,189)*y(k,155) + rxt(k,193) &
                      *y(k,233) + rxt(k,194)*y(k,163) + rxt(k,195)*y(k,162) + rxt(k,218) &
                      *y(k,74) + rxt(k,252)*y(k,21) + rxt(k,280)*y(k,115) + rxt(k,289) &
                      *y(k,125) + rxt(k,296)*y(k,167) + rxt(k,336)*y(k,68) + rxt(k,355) &
                      *y(k,228) + rxt(k,363)*y(k,234) + rxt(k,376)*y(k,224) + rxt(k,387) &
                      *y(k,227) + rxt(k,391)*y(k,232) + rxt(k,404)*y(k,225) + rxt(k,413) &
                      *y(k,250) + rxt(k,417)*y(k,251) + (rxt(k,423) + rxt(k,424) &
                      ) * y(k,230) + (rxt(k,430) + rxt(k,431)) * y(k,240) + rxt(k,439) &
                      *y(k,242) + rxt(k,442)*y(k,244) + (rxt(k,452) + rxt(k,453) &
                      ) * y(k,219) + rxt(k,462)*y(k,235) + rxt(k,468)*y(k,236) &
                      + rxt(k,476)*y(k,237) + rxt(k,487)*y(k,256) + rxt(k,491) &
                      *y(k,218) + rxt(k,494)*y(k,221) + rxt(k,499)*y(k,223) + rxt(k,501) &
                      *y(k,226) + rxt(k,505)*y(k,229) + rxt(k,508)*y(k,241) + rxt(k,511) &
                      *y(k,243) + rxt(k,514)*y(k,249) + rxt(k,521)*y(k,254) + rxt(k,527) &
                      *y(k,257) + rxt(k,530)*y(k,259) + rxt(k,541)*y(k,246) + rxt(k,546) &
                      *y(k,252) + rxt(k,551)*y(k,253))
         mat(k,582) = -rxt(k,180)*y(k,153)
         mat(k,1823) = -rxt(k,189)*y(k,153)
         mat(k,2348) = -rxt(k,193)*y(k,153)
         mat(k,1949) = -rxt(k,194)*y(k,153)
         mat(k,2440) = -rxt(k,195)*y(k,153)
         mat(k,2377) = -rxt(k,218)*y(k,153)
         mat(k,2467) = -rxt(k,252)*y(k,153)
         mat(k,2405) = -rxt(k,280)*y(k,153)
         mat(k,2233) = -rxt(k,289)*y(k,153)
         mat(k,597) = -rxt(k,296)*y(k,153)
         mat(k,1122) = -rxt(k,336)*y(k,153)
         mat(k,1654) = -rxt(k,355)*y(k,153)
         mat(k,525) = -rxt(k,363)*y(k,153)
         mat(k,885) = -rxt(k,376)*y(k,153)
         mat(k,1502) = -rxt(k,387)*y(k,153)
         mat(k,803) = -rxt(k,391)*y(k,153)
         mat(k,944) = -rxt(k,404)*y(k,153)
         mat(k,866) = -rxt(k,413)*y(k,153)
         mat(k,1282) = -rxt(k,417)*y(k,153)
         mat(k,656) = -(rxt(k,423) + rxt(k,424)) * y(k,153)
         mat(k,1427) = -(rxt(k,430) + rxt(k,431)) * y(k,153)
         mat(k,1470) = -rxt(k,439)*y(k,153)
         mat(k,774) = -rxt(k,442)*y(k,153)
         mat(k,1007) = -(rxt(k,452) + rxt(k,453)) * y(k,153)
         mat(k,1373) = -rxt(k,462)*y(k,153)
         mat(k,1406) = -rxt(k,468)*y(k,153)
         mat(k,1324) = -rxt(k,476)*y(k,153)
         mat(k,1301) = -rxt(k,487)*y(k,153)
         mat(k,615) = -rxt(k,491)*y(k,153)
         mat(k,574) = -rxt(k,494)*y(k,153)
         mat(k,519) = -rxt(k,499)*y(k,153)
         mat(k,734) = -rxt(k,501)*y(k,153)
         mat(k,857) = -rxt(k,505)*y(k,153)
         mat(k,817) = -rxt(k,508)*y(k,153)
         mat(k,964) = -rxt(k,511)*y(k,153)
         mat(k,532) = -rxt(k,514)*y(k,153)
         mat(k,832) = -rxt(k,521)*y(k,153)
         mat(k,849) = -rxt(k,527)*y(k,153)
         mat(k,590) = -rxt(k,530)*y(k,153)
         mat(k,1177) = -rxt(k,541)*y(k,153)
         mat(k,1230) = -rxt(k,546)*y(k,153)
         mat(k,1097) = -rxt(k,551)*y(k,153)
         mat(k,212) = 4.000_r8*rxt(k,279)*y(k,123)
         mat(k,582) = mat(k,582) + 2.000_r8*rxt(k,182)*y(k,154) + rxt(k,192)*y(k,248)
         mat(k,215) = 2.000_r8*rxt(k,196)*y(k,247)
         mat(k,2575) = 2.000_r8*rxt(k,182)*y(k,141) + rxt(k,185)*y(k,162) + rxt(k,570) &
                      *y(k,180)
         mat(k,2440) = mat(k,2440) + rxt(k,185)*y(k,154)
         mat(k,1560) = rxt(k,570)*y(k,154)
         mat(k,1995) = 2.000_r8*rxt(k,196)*y(k,142)
         mat(k,2177) = rxt(k,192)*y(k,141)
         mat(k,2574) = -((rxt(k,181) + rxt(k,182) + rxt(k,183)) * y(k,141) + (rxt(k,185) &
                      + rxt(k,187)) * y(k,162) + rxt(k,186)*y(k,163) + rxt(k,198) &
                      *y(k,233) + rxt(k,199)*y(k,155) + rxt(k,200)*y(k,248) + rxt(k,210) &
                      *y(k,70) + rxt(k,220)*y(k,74) + rxt(k,245)*y(k,17) + rxt(k,255) &
                      *y(k,21) + rxt(k,276)*y(k,115) + rxt(k,290)*y(k,125) + rxt(k,398) &
                      *y(k,227) + rxt(k,448)*y(k,242) + rxt(k,506)*y(k,229) + rxt(k,509) &
                      *y(k,241) + rxt(k,512)*y(k,243) + rxt(k,516)*y(k,171) + rxt(k,519) &
                      *y(k,218) + rxt(k,570)*y(k,180))
         mat(k,581) = -(rxt(k,181) + rxt(k,182) + rxt(k,183)) * y(k,154)
         mat(k,2439) = -(rxt(k,185) + rxt(k,187)) * y(k,154)
         mat(k,1948) = -rxt(k,186)*y(k,154)
         mat(k,2347) = -rxt(k,198)*y(k,154)
         mat(k,1822) = -rxt(k,199)*y(k,154)
         mat(k,2176) = -rxt(k,200)*y(k,154)
         mat(k,1882) = -rxt(k,210)*y(k,154)
         mat(k,2376) = -rxt(k,220)*y(k,154)
         mat(k,2494) = -rxt(k,245)*y(k,154)
         mat(k,2466) = -rxt(k,255)*y(k,154)
         mat(k,2404) = -rxt(k,276)*y(k,154)
         mat(k,2232) = -rxt(k,290)*y(k,154)
         mat(k,1501) = -rxt(k,398)*y(k,154)
         mat(k,1469) = -rxt(k,448)*y(k,154)
         mat(k,856) = -rxt(k,506)*y(k,154)
         mat(k,816) = -rxt(k,509)*y(k,154)
         mat(k,963) = -rxt(k,512)*y(k,154)
         mat(k,548) = -rxt(k,516)*y(k,154)
         mat(k,614) = -rxt(k,519)*y(k,154)
         mat(k,1559) = -rxt(k,570)*y(k,154)
         mat(k,766) = rxt(k,450)*y(k,248)
         mat(k,418) = rxt(k,421)*y(k,155)
         mat(k,2466) = mat(k,2466) + rxt(k,252)*y(k,153)
         mat(k,1121) = rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155)
         mat(k,624) = rxt(k,211)*y(k,248)
         mat(k,2376) = mat(k,2376) + rxt(k,218)*y(k,153)
         mat(k,474) = rxt(k,179)*y(k,248)
         mat(k,2404) = mat(k,2404) + rxt(k,278)*y(k,155)
         mat(k,313) = 4.000_r8*rxt(k,275)*y(k,124)
         mat(k,2232) = mat(k,2232) + rxt(k,289)*y(k,153) + rxt(k,291)*y(k,155)
         mat(k,698) = .700_r8*rxt(k,470)*y(k,248)
         mat(k,2675) = rxt(k,252)*y(k,21) + rxt(k,336)*y(k,68) + rxt(k,218)*y(k,74) &
                      + rxt(k,289)*y(k,125) + 2.000_r8*rxt(k,189)*y(k,155) &
                      + rxt(k,195)*y(k,162) + rxt(k,194)*y(k,163) + rxt(k,296) &
                      *y(k,167) + rxt(k,491)*y(k,218) + rxt(k,452)*y(k,219) &
                      + rxt(k,494)*y(k,221) + rxt(k,499)*y(k,223) + rxt(k,376) &
                      *y(k,224) + rxt(k,404)*y(k,225) + rxt(k,501)*y(k,226) &
                      + rxt(k,387)*y(k,227) + rxt(k,355)*y(k,228) + rxt(k,505) &
                      *y(k,229) + rxt(k,423)*y(k,230) + rxt(k,391)*y(k,232) &
                      + rxt(k,193)*y(k,233) + rxt(k,363)*y(k,234) + .920_r8*rxt(k,462) &
                      *y(k,235) + .920_r8*rxt(k,468)*y(k,236) + rxt(k,476)*y(k,237) &
                      + rxt(k,430)*y(k,240) + rxt(k,508)*y(k,241) + rxt(k,439) &
                      *y(k,242) + rxt(k,511)*y(k,243) + rxt(k,442)*y(k,244) &
                      + 1.600_r8*rxt(k,541)*y(k,246) + rxt(k,514)*y(k,249) &
                      + rxt(k,413)*y(k,250) + rxt(k,417)*y(k,251) + .900_r8*rxt(k,546) &
                      *y(k,252) + .800_r8*rxt(k,551)*y(k,253) + rxt(k,521)*y(k,254) &
                      + rxt(k,487)*y(k,256) + rxt(k,527)*y(k,257) + rxt(k,530) &
                      *y(k,259)
         mat(k,1822) = mat(k,1822) + rxt(k,421)*y(k,16) + rxt(k,337)*y(k,68) &
                      + rxt(k,278)*y(k,115) + rxt(k,291)*y(k,125) &
                      + 2.000_r8*rxt(k,189)*y(k,153) + rxt(k,190)*y(k,162) &
                      + rxt(k,188)*y(k,233) + rxt(k,463)*y(k,235) + rxt(k,469) &
                      *y(k,236) + rxt(k,477)*y(k,237) + rxt(k,429)*y(k,240) &
                      + rxt(k,440)*y(k,242) + 2.000_r8*rxt(k,542)*y(k,246) &
                      + rxt(k,191)*y(k,248) + rxt(k,488)*y(k,256)
         mat(k,903) = rxt(k,411)*y(k,248)
         mat(k,2439) = mat(k,2439) + rxt(k,195)*y(k,153) + rxt(k,190)*y(k,155)
         mat(k,1948) = mat(k,1948) + rxt(k,194)*y(k,153)
         mat(k,596) = rxt(k,296)*y(k,153)
         mat(k,726) = rxt(k,548)*y(k,248)
         mat(k,614) = mat(k,614) + rxt(k,491)*y(k,153)
         mat(k,1006) = rxt(k,452)*y(k,153)
         mat(k,573) = rxt(k,494)*y(k,153)
         mat(k,518) = rxt(k,499)*y(k,153)
         mat(k,884) = rxt(k,376)*y(k,153)
         mat(k,943) = rxt(k,404)*y(k,153)
         mat(k,733) = rxt(k,501)*y(k,153)
         mat(k,1501) = mat(k,1501) + rxt(k,387)*y(k,153)
         mat(k,1653) = rxt(k,355)*y(k,153) + .500_r8*rxt(k,539)*y(k,246)
         mat(k,856) = mat(k,856) + rxt(k,505)*y(k,153)
         mat(k,655) = rxt(k,423)*y(k,153)
         mat(k,802) = rxt(k,391)*y(k,153)
         mat(k,2347) = mat(k,2347) + rxt(k,193)*y(k,153) + rxt(k,188)*y(k,155)
         mat(k,524) = rxt(k,363)*y(k,153)
         mat(k,1372) = .920_r8*rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155)
         mat(k,1405) = .920_r8*rxt(k,468)*y(k,153) + rxt(k,469)*y(k,155)
         mat(k,1323) = rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155)
         mat(k,1426) = rxt(k,430)*y(k,153) + rxt(k,429)*y(k,155)
         mat(k,816) = mat(k,816) + rxt(k,508)*y(k,153)
         mat(k,1469) = mat(k,1469) + rxt(k,439)*y(k,153) + rxt(k,440)*y(k,155)
         mat(k,963) = mat(k,963) + rxt(k,511)*y(k,153)
         mat(k,773) = rxt(k,442)*y(k,153)
         mat(k,1176) = 1.600_r8*rxt(k,541)*y(k,153) + 2.000_r8*rxt(k,542)*y(k,155) &
                      + .500_r8*rxt(k,539)*y(k,228)
         mat(k,2176) = mat(k,2176) + rxt(k,450)*y(k,1) + rxt(k,211)*y(k,73) &
                      + rxt(k,179)*y(k,107) + .700_r8*rxt(k,470)*y(k,128) + rxt(k,191) &
                      *y(k,155) + rxt(k,411)*y(k,156) + rxt(k,548)*y(k,205)
         mat(k,531) = rxt(k,514)*y(k,153)
         mat(k,865) = rxt(k,413)*y(k,153)
         mat(k,1281) = rxt(k,417)*y(k,153)
         mat(k,1229) = .900_r8*rxt(k,546)*y(k,153)
         mat(k,1096) = .800_r8*rxt(k,551)*y(k,153)
         mat(k,831) = rxt(k,521)*y(k,153)
         mat(k,1300) = rxt(k,487)*y(k,153) + rxt(k,488)*y(k,155)
         mat(k,848) = rxt(k,527)*y(k,153)
         mat(k,589) = rxt(k,530)*y(k,153)
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
         mat(k,1808) = -(rxt(k,188)*y(k,233) + rxt(k,189)*y(k,153) + rxt(k,190) &
                      *y(k,162) + rxt(k,191)*y(k,248) + rxt(k,199)*y(k,154) + rxt(k,267) &
                      *y(k,105) + rxt(k,269)*y(k,116) + rxt(k,278)*y(k,115) + rxt(k,291) &
                      *y(k,125) + rxt(k,297)*y(k,110) + rxt(k,337)*y(k,68) + rxt(k,349) &
                      *y(k,51) + rxt(k,381)*y(k,54) + rxt(k,400)*y(k,33) + rxt(k,407) &
                      *y(k,58) + rxt(k,421)*y(k,16) + rxt(k,429)*y(k,240) + rxt(k,440) &
                      *y(k,242) + rxt(k,463)*y(k,235) + rxt(k,469)*y(k,236) + rxt(k,472) &
                      *y(k,127) + rxt(k,477)*y(k,237) + rxt(k,488)*y(k,256) + rxt(k,533) &
                      *y(k,6) + rxt(k,536)*y(k,139) + rxt(k,542)*y(k,246) + rxt(k,553) &
                      *y(k,207) + rxt(k,560)*y(k,83))
         mat(k,2333) = -rxt(k,188)*y(k,155)
         mat(k,2661) = -rxt(k,189)*y(k,155)
         mat(k,2425) = -rxt(k,190)*y(k,155)
         mat(k,2162) = -rxt(k,191)*y(k,155)
         mat(k,2560) = -rxt(k,199)*y(k,155)
         mat(k,1263) = -rxt(k,267)*y(k,155)
         mat(k,1523) = -rxt(k,269)*y(k,155)
         mat(k,2390) = -rxt(k,278)*y(k,155)
         mat(k,2218) = -rxt(k,291)*y(k,155)
         mat(k,1589) = -rxt(k,297)*y(k,155)
         mat(k,1118) = -rxt(k,337)*y(k,155)
         mat(k,1691) = -rxt(k,349)*y(k,155)
         mat(k,1194) = -rxt(k,381)*y(k,155)
         mat(k,1138) = -rxt(k,400)*y(k,155)
         mat(k,1341) = -rxt(k,407)*y(k,155)
         mat(k,416) = -rxt(k,421)*y(k,155)
         mat(k,1421) = -rxt(k,429)*y(k,155)
         mat(k,1463) = -rxt(k,440)*y(k,155)
         mat(k,1366) = -rxt(k,463)*y(k,155)
         mat(k,1399) = -rxt(k,469)*y(k,155)
         mat(k,931) = -rxt(k,472)*y(k,155)
         mat(k,1318) = -rxt(k,477)*y(k,155)
         mat(k,1296) = -rxt(k,488)*y(k,155)
         mat(k,1077) = -rxt(k,533)*y(k,155)
         mat(k,988) = -rxt(k,536)*y(k,155)
         mat(k,1172) = -rxt(k,542)*y(k,155)
         mat(k,1103) = -rxt(k,553)*y(k,155)
         mat(k,1029) = -rxt(k,560)*y(k,155)
         mat(k,2480) = rxt(k,253)*y(k,22)
         mat(k,908) = rxt(k,253)*y(k,17) + rxt(k,254)*y(k,70) + rxt(k,256)*y(k,162)
         mat(k,1868) = rxt(k,254)*y(k,22) + rxt(k,219)*y(k,75)
         mat(k,1012) = rxt(k,219)*y(k,70) + rxt(k,221)*y(k,162) + rxt(k,222)*y(k,248)
         mat(k,951) = rxt(k,309)*y(k,106)
         mat(k,2187) = rxt(k,309)*y(k,89) + rxt(k,201)*y(k,248)
         mat(k,2390) = mat(k,2390) + rxt(k,274)*y(k,126)
         mat(k,891) = rxt(k,274)*y(k,115)
         mat(k,686) = .500_r8*rxt(k,445)*y(k,248)
         mat(k,2560) = mat(k,2560) + rxt(k,187)*y(k,162) + rxt(k,186)*y(k,163)
         mat(k,2425) = mat(k,2425) + rxt(k,256)*y(k,22) + rxt(k,221)*y(k,75) &
                      + rxt(k,187)*y(k,154)
         mat(k,1934) = rxt(k,186)*y(k,154)
         mat(k,630) = rxt(k,396)*y(k,248)
         mat(k,2162) = mat(k,2162) + rxt(k,222)*y(k,75) + rxt(k,201)*y(k,106) &
                      + .500_r8*rxt(k,445)*y(k,138) + rxt(k,396)*y(k,169)
         mat(k,898) = -(rxt(k,411)*y(k,248))
         mat(k,2112) = -rxt(k,411)*y(k,156)
         mat(k,1127) = rxt(k,400)*y(k,155)
         mat(k,658) = .500_r8*rxt(k,471)*y(k,248)
         mat(k,448) = rxt(k,478)*y(k,248)
         mat(k,459) = rxt(k,482)*y(k,248)
         mat(k,1147) = rxt(k,483)*y(k,248)
         mat(k,1764) = rxt(k,400)*y(k,33)
         mat(k,2112) = mat(k,2112) + .500_r8*rxt(k,471)*y(k,129) + rxt(k,478)*y(k,130) &
                      + rxt(k,482)*y(k,144) + rxt(k,483)*y(k,145)
         mat(k,464) = -(rxt(k,543)*y(k,248))
         mat(k,2060) = -rxt(k,543)*y(k,157)
         mat(k,2258) = rxt(k,540)*y(k,246)
         mat(k,1162) = rxt(k,540)*y(k,233)
         mat(k,2435) = -(rxt(k,159)*y(k,163) + 4._r8*rxt(k,160)*y(k,162) + rxt(k,162) &
                      *y(k,93) + rxt(k,163)*y(k,95) + rxt(k,168)*y(k,233) + rxt(k,174) &
                      *y(k,248) + (rxt(k,185) + rxt(k,187)) * y(k,154) + rxt(k,190) &
                      *y(k,155) + rxt(k,195)*y(k,153) + rxt(k,221)*y(k,75) + rxt(k,223) &
                      *y(k,74) + rxt(k,226)*y(k,101) + rxt(k,229)*y(k,109) + rxt(k,256) &
                      *y(k,22) + rxt(k,257)*y(k,21) + rxt(k,259)*y(k,97) + rxt(k,261) &
                      *y(k,108) + rxt(k,270)*y(k,116) + rxt(k,292)*y(k,125) + rxt(k,350) &
                      *y(k,51) + rxt(k,562)*y(k,166))
         mat(k,1944) = -rxt(k,159)*y(k,162)
         mat(k,1543) = -rxt(k,162)*y(k,162)
         mat(k,719) = -rxt(k,163)*y(k,162)
         mat(k,2343) = -rxt(k,168)*y(k,162)
         mat(k,2172) = -rxt(k,174)*y(k,162)
         mat(k,2570) = -(rxt(k,185) + rxt(k,187)) * y(k,162)
         mat(k,1818) = -rxt(k,190)*y(k,162)
         mat(k,2671) = -rxt(k,195)*y(k,162)
         mat(k,1017) = -rxt(k,221)*y(k,162)
         mat(k,2372) = -rxt(k,223)*y(k,162)
         mat(k,1749) = -rxt(k,226)*y(k,162)
         mat(k,1726) = -rxt(k,229)*y(k,162)
         mat(k,912) = -rxt(k,256)*y(k,162)
         mat(k,2462) = -rxt(k,257)*y(k,162)
         mat(k,1515) = -rxt(k,259)*y(k,162)
         mat(k,1673) = -rxt(k,261)*y(k,162)
         mat(k,1529) = -rxt(k,270)*y(k,162)
         mat(k,2228) = -rxt(k,292)*y(k,162)
         mat(k,1701) = -rxt(k,350)*y(k,162)
         mat(k,425) = -rxt(k,562)*y(k,162)
         mat(k,2513) = rxt(k,166)*y(k,233)
         mat(k,579) = rxt(k,180)*y(k,153) + rxt(k,181)*y(k,154)
         mat(k,2671) = mat(k,2671) + rxt(k,180)*y(k,141)
         mat(k,2570) = mat(k,2570) + rxt(k,181)*y(k,141)
         mat(k,1944) = mat(k,1944) + 2.000_r8*rxt(k,158)*y(k,247)
         mat(k,2343) = mat(k,2343) + rxt(k,166)*y(k,92)
         mat(k,1990) = 2.000_r8*rxt(k,158)*y(k,163)
         mat(k,2172) = mat(k,2172) + 2.000_r8*rxt(k,176)*y(k,248)
         mat(k,1936) = -((rxt(k,157) + rxt(k,158)) * y(k,247) + rxt(k,159)*y(k,162) &
                      + rxt(k,169)*y(k,233) + rxt(k,170)*y(k,92) + rxt(k,175)*y(k,248) &
                      + rxt(k,186)*y(k,154) + rxt(k,194)*y(k,153) + rxt(k,212)*y(k,70) &
                      + rxt(k,246)*y(k,17) + rxt(k,281)*y(k,115) + rxt(k,293)*y(k,125) &
                      + rxt(k,372)*y(k,28) + rxt(k,401)*y(k,33) + rxt(k,432)*y(k,134) &
                      + rxt(k,446)*y(k,140) + rxt(k,479)*y(k,127) + rxt(k,517) &
                      *y(k,171) + rxt(k,534)*y(k,6) + rxt(k,537)*y(k,139) + rxt(k,566) &
                      *y(k,178) + rxt(k,572)*y(k,180))
         mat(k,1982) = -(rxt(k,157) + rxt(k,158)) * y(k,163)
         mat(k,2427) = -rxt(k,159)*y(k,163)
         mat(k,2335) = -rxt(k,169)*y(k,163)
         mat(k,2505) = -rxt(k,170)*y(k,163)
         mat(k,2164) = -rxt(k,175)*y(k,163)
         mat(k,2562) = -rxt(k,186)*y(k,163)
         mat(k,2663) = -rxt(k,194)*y(k,163)
         mat(k,1870) = -rxt(k,212)*y(k,163)
         mat(k,2482) = -rxt(k,246)*y(k,163)
         mat(k,2392) = -rxt(k,281)*y(k,163)
         mat(k,2220) = -rxt(k,293)*y(k,163)
         mat(k,646) = -rxt(k,372)*y(k,163)
         mat(k,1139) = -rxt(k,401)*y(k,163)
         mat(k,1333) = -rxt(k,432)*y(k,163)
         mat(k,1444) = -rxt(k,446)*y(k,163)
         mat(k,932) = -rxt(k,479)*y(k,163)
         mat(k,547) = -rxt(k,517)*y(k,163)
         mat(k,1078) = -rxt(k,534)*y(k,163)
         mat(k,989) = -rxt(k,537)*y(k,163)
         mat(k,605) = -rxt(k,566)*y(k,163)
         mat(k,1551) = -rxt(k,572)*y(k,163)
         mat(k,1496) = .150_r8*rxt(k,386)*y(k,233)
         mat(k,2335) = mat(k,2335) + .150_r8*rxt(k,386)*y(k,227) + .150_r8*rxt(k,437) &
                      *y(k,242)
         mat(k,1464) = .150_r8*rxt(k,437)*y(k,233)
         mat(k,533) = -(rxt(k,573)*y(k,180))
         mat(k,1546) = -rxt(k,573)*y(k,165)
         mat(k,2442) = rxt(k,248)*y(k,74)
         mat(k,2352) = rxt(k,248)*y(k,21) + 2.000_r8*rxt(k,216)*y(k,74) + rxt(k,285) &
                      *y(k,125)
         mat(k,2206) = rxt(k,285)*y(k,74)
         mat(k,419) = -(rxt(k,562)*y(k,162) + rxt(k,563)*y(k,248))
         mat(k,2407) = -rxt(k,562)*y(k,166)
         mat(k,2055) = -rxt(k,563)*y(k,166)
         mat(k,593) = -(rxt(k,296)*y(k,153) + rxt(k,304)*y(k,125) + 4._r8*rxt(k,305) &
                      *y(k,167))
         mat(k,2604) = -rxt(k,296)*y(k,167)
         mat(k,2208) = -rxt(k,304)*y(k,167)
         mat(k,2443) = rxt(k,284)*y(k,125)
         mat(k,2208) = mat(k,2208) + rxt(k,284)*y(k,21) + 2.000_r8*rxt(k,301)*y(k,125) &
                      + rxt(k,291)*y(k,155) + rxt(k,293)*y(k,163)
         mat(k,1760) = rxt(k,291)*y(k,125)
         mat(k,1893) = rxt(k,293)*y(k,125)
         mat(k,1232) = rxt(k,425)*y(k,248)
         mat(k,2591) = .100_r8*rxt(k,546)*y(k,252)
         mat(k,2036) = rxt(k,425)*y(k,111)
         mat(k,1212) = .100_r8*rxt(k,546)*y(k,153)
         mat(k,625) = -(rxt(k,396)*y(k,248))
         mat(k,2082) = -rxt(k,396)*y(k,169)
         mat(k,2534) = rxt(k,398)*y(k,227)
         mat(k,1473) = rxt(k,398)*y(k,154)
         mat(k,2520) = rxt(k,519)*y(k,218)
         mat(k,609) = rxt(k,519)*y(k,154)
         mat(k,545) = -(rxt(k,516)*y(k,154) + rxt(k,517)*y(k,163))
         mat(k,2528) = -rxt(k,516)*y(k,171)
         mat(k,1892) = -rxt(k,517)*y(k,171)
         mat(k,230) = .070_r8*rxt(k,503)*y(k,248)
         mat(k,2600) = rxt(k,501)*y(k,226)
         mat(k,200) = .060_r8*rxt(k,515)*y(k,248)
         mat(k,251) = .070_r8*rxt(k,531)*y(k,248)
         mat(k,728) = rxt(k,501)*y(k,153)
         mat(k,2072) = .070_r8*rxt(k,503)*y(k,82) + .060_r8*rxt(k,515)*y(k,172) &
                      + .070_r8*rxt(k,531)*y(k,214)
         mat(k,198) = -(rxt(k,515)*y(k,248))
         mat(k,2020) = -rxt(k,515)*y(k,172)
         mat(k,190) = .530_r8*rxt(k,492)*y(k,248)
         mat(k,2020) = mat(k,2020) + .530_r8*rxt(k,492)*y(k,7)
         mat(k,376) = -(rxt(k,518)*y(k,248))
         mat(k,2049) = -rxt(k,518)*y(k,173)
         mat(k,2253) = rxt(k,513)*y(k,249)
         mat(k,526) = rxt(k,513)*y(k,233)
         mat(k,633) = -(rxt(k,414)*y(k,248))
         mat(k,2083) = -rxt(k,414)*y(k,176)
         mat(k,2274) = rxt(k,412)*y(k,250)
         mat(k,858) = rxt(k,412)*y(k,233)
         mat(k,482) = -(rxt(k,418)*y(k,248))
         mat(k,2063) = -rxt(k,418)*y(k,177)
         mat(k,2261) = .850_r8*rxt(k,416)*y(k,251)
         mat(k,1271) = .850_r8*rxt(k,416)*y(k,233)
         mat(k,603) = -(rxt(k,566)*y(k,163) + rxt(k,569)*y(k,248))
         mat(k,1894) = -rxt(k,566)*y(k,178)
         mat(k,2079) = -rxt(k,569)*y(k,178)
         mat(k,1549) = -(rxt(k,567)*y(k,21) + rxt(k,568)*y(k,74) + rxt(k,570)*y(k,154) &
                      + rxt(k,572)*y(k,163) + rxt(k,573)*y(k,165) + rxt(k,574) &
                      *y(k,248))
         mat(k,2448) = -rxt(k,567)*y(k,180)
         mat(k,2357) = -rxt(k,568)*y(k,180)
         mat(k,2552) = -rxt(k,570)*y(k,180)
         mat(k,1928) = -rxt(k,572)*y(k,180)
         mat(k,535) = -rxt(k,573)*y(k,180)
         mat(k,2154) = -rxt(k,574)*y(k,180)
         mat(k,2417) = rxt(k,562)*y(k,166)
         mat(k,1928) = mat(k,1928) + rxt(k,566)*y(k,178)
         mat(k,423) = rxt(k,562)*y(k,162)
         mat(k,604) = rxt(k,566)*y(k,163) + rxt(k,569)*y(k,248)
         mat(k,2154) = mat(k,2154) + rxt(k,569)*y(k,178)
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
         mat(k,1039) = -(rxt(k,565)*y(k,248))
         mat(k,2122) = -rxt(k,565)*y(k,181)
         mat(k,2447) = rxt(k,556)*y(k,83) + rxt(k,567)*y(k,180)
         mat(k,1850) = rxt(k,558)*y(k,83)
         mat(k,2356) = rxt(k,568)*y(k,180)
         mat(k,1027) = rxt(k,556)*y(k,21) + rxt(k,558)*y(k,70) + rxt(k,559)*y(k,125) &
                      + rxt(k,560)*y(k,155) + (rxt(k,561)+.500_r8*rxt(k,575))*y(k,248)
         mat(k,2211) = rxt(k,559)*y(k,83)
         mat(k,2544) = rxt(k,570)*y(k,180)
         mat(k,1770) = rxt(k,560)*y(k,83)
         mat(k,1904) = rxt(k,572)*y(k,180)
         mat(k,534) = rxt(k,573)*y(k,180)
         mat(k,421) = rxt(k,563)*y(k,248)
         mat(k,1548) = rxt(k,567)*y(k,21) + rxt(k,568)*y(k,74) + rxt(k,570)*y(k,154) &
                      + rxt(k,572)*y(k,163) + rxt(k,573)*y(k,165) + rxt(k,574) &
                      *y(k,248)
         mat(k,2122) = mat(k,2122) + (rxt(k,561)+.500_r8*rxt(k,575))*y(k,83) &
                      + rxt(k,563)*y(k,166) + rxt(k,574)*y(k,180)
         mat(k,292) = -(rxt(k,576)*y(k,260))
         mat(k,2679) = -rxt(k,576)*y(k,182)
         mat(k,1038) = rxt(k,565)*y(k,248)
         mat(k,2038) = rxt(k,565)*y(k,181)
         mat(k,1052) = .2202005_r8*rxt(k,629)*y(k,163)
         mat(k,965) = .0508005_r8*rxt(k,645)*y(k,163)
         mat(k,2577) = .1279005_r8*rxt(k,628)*y(k,220) + .0097005_r8*rxt(k,633) &
                      *y(k,222) + .0003005_r8*rxt(k,636)*y(k,238) &
                      + .1056005_r8*rxt(k,640)*y(k,239) + .0245005_r8*rxt(k,644) &
                      *y(k,245) + .0154005_r8*rxt(k,650)*y(k,255) &
                      + .0063005_r8*rxt(k,654)*y(k,258)
         mat(k,1885) = .2202005_r8*rxt(k,629)*y(k,6) + .0508005_r8*rxt(k,645)*y(k,139)
         mat(k,52) = .5931005_r8*rxt(k,647)*y(k,248)
         mat(k,58) = .1279005_r8*rxt(k,628)*y(k,153) + .2202005_r8*rxt(k,627)*y(k,233)
         mat(k,64) = .0097005_r8*rxt(k,633)*y(k,153) + .0023005_r8*rxt(k,632)*y(k,233)
         mat(k,2235) = .2202005_r8*rxt(k,627)*y(k,220) + .0023005_r8*rxt(k,632) &
                      *y(k,222) + .0031005_r8*rxt(k,635)*y(k,238) &
                      + .2381005_r8*rxt(k,639)*y(k,239) + .0508005_r8*rxt(k,643) &
                      *y(k,245) + .1364005_r8*rxt(k,649)*y(k,255) &
                      + .1677005_r8*rxt(k,653)*y(k,258)
         mat(k,70) = .0003005_r8*rxt(k,636)*y(k,153) + .0031005_r8*rxt(k,635)*y(k,233)
         mat(k,76) = .1056005_r8*rxt(k,640)*y(k,153) + .2381005_r8*rxt(k,639)*y(k,233)
         mat(k,84) = .0245005_r8*rxt(k,644)*y(k,153) + .0508005_r8*rxt(k,643)*y(k,233)
         mat(k,1997) = .5931005_r8*rxt(k,647)*y(k,202)
         mat(k,90) = .0154005_r8*rxt(k,650)*y(k,153) + .1364005_r8*rxt(k,649)*y(k,233)
         mat(k,96) = .0063005_r8*rxt(k,654)*y(k,153) + .1677005_r8*rxt(k,653)*y(k,233)
         mat(k,1053) = .2067005_r8*rxt(k,629)*y(k,163)
         mat(k,966) = .1149005_r8*rxt(k,645)*y(k,163)
         mat(k,2578) = .1792005_r8*rxt(k,628)*y(k,220) + .0034005_r8*rxt(k,633) &
                      *y(k,222) + .0003005_r8*rxt(k,636)*y(k,238) &
                      + .1026005_r8*rxt(k,640)*y(k,239) + .0082005_r8*rxt(k,644) &
                      *y(k,245) + .0452005_r8*rxt(k,650)*y(k,255) &
                      + .0237005_r8*rxt(k,654)*y(k,258)
         mat(k,1886) = .2067005_r8*rxt(k,629)*y(k,6) + .1149005_r8*rxt(k,645)*y(k,139)
         mat(k,53) = .1534005_r8*rxt(k,647)*y(k,248)
         mat(k,59) = .1792005_r8*rxt(k,628)*y(k,153) + .2067005_r8*rxt(k,627)*y(k,233)
         mat(k,65) = .0034005_r8*rxt(k,633)*y(k,153) + .0008005_r8*rxt(k,632)*y(k,233)
         mat(k,2236) = .2067005_r8*rxt(k,627)*y(k,220) + .0008005_r8*rxt(k,632) &
                      *y(k,222) + .0035005_r8*rxt(k,635)*y(k,238) &
                      + .1308005_r8*rxt(k,639)*y(k,239) + .1149005_r8*rxt(k,643) &
                      *y(k,245) + .0101005_r8*rxt(k,649)*y(k,255) &
                      + .0174005_r8*rxt(k,653)*y(k,258)
         mat(k,71) = .0003005_r8*rxt(k,636)*y(k,153) + .0035005_r8*rxt(k,635)*y(k,233)
         mat(k,77) = .1026005_r8*rxt(k,640)*y(k,153) + .1308005_r8*rxt(k,639)*y(k,233)
         mat(k,85) = .0082005_r8*rxt(k,644)*y(k,153) + .1149005_r8*rxt(k,643)*y(k,233)
         mat(k,1998) = .1534005_r8*rxt(k,647)*y(k,202)
         mat(k,91) = .0452005_r8*rxt(k,650)*y(k,153) + .0101005_r8*rxt(k,649)*y(k,233)
         mat(k,97) = .0237005_r8*rxt(k,654)*y(k,153) + .0174005_r8*rxt(k,653)*y(k,233)
         mat(k,1054) = .0653005_r8*rxt(k,629)*y(k,163)
         mat(k,967) = .0348005_r8*rxt(k,645)*y(k,163)
         mat(k,2579) = .0676005_r8*rxt(k,628)*y(k,220) + .1579005_r8*rxt(k,633) &
                      *y(k,222) + .0073005_r8*rxt(k,636)*y(k,238) &
                      + .0521005_r8*rxt(k,640)*y(k,239) + .0772005_r8*rxt(k,644) &
                      *y(k,245) + .0966005_r8*rxt(k,650)*y(k,255) &
                      + .0025005_r8*rxt(k,654)*y(k,258)
         mat(k,1887) = .0653005_r8*rxt(k,629)*y(k,6) + .0348005_r8*rxt(k,645)*y(k,139)
         mat(k,54) = .0459005_r8*rxt(k,647)*y(k,248)
         mat(k,60) = .0676005_r8*rxt(k,628)*y(k,153) + .0653005_r8*rxt(k,627)*y(k,233)
         mat(k,66) = .1579005_r8*rxt(k,633)*y(k,153) + .0843005_r8*rxt(k,632)*y(k,233)
         mat(k,2237) = .0653005_r8*rxt(k,627)*y(k,220) + .0843005_r8*rxt(k,632) &
                      *y(k,222) + .0003005_r8*rxt(k,635)*y(k,238) &
                      + .0348005_r8*rxt(k,639)*y(k,239) + .0348005_r8*rxt(k,643) &
                      *y(k,245) + .0763005_r8*rxt(k,649)*y(k,255) + .086_r8*rxt(k,653) &
                      *y(k,258)
         mat(k,72) = .0073005_r8*rxt(k,636)*y(k,153) + .0003005_r8*rxt(k,635)*y(k,233)
         mat(k,78) = .0521005_r8*rxt(k,640)*y(k,153) + .0348005_r8*rxt(k,639)*y(k,233)
         mat(k,86) = .0772005_r8*rxt(k,644)*y(k,153) + .0348005_r8*rxt(k,643)*y(k,233)
         mat(k,1999) = .0459005_r8*rxt(k,647)*y(k,202)
         mat(k,92) = .0966005_r8*rxt(k,650)*y(k,153) + .0763005_r8*rxt(k,649)*y(k,233)
         mat(k,98) = .0025005_r8*rxt(k,654)*y(k,153) + .086_r8*rxt(k,653)*y(k,233)
         mat(k,1055) = .1749305_r8*rxt(k,626)*y(k,155) + .1284005_r8*rxt(k,629) &
                      *y(k,163)
         mat(k,916) = .0590245_r8*rxt(k,634)*y(k,155) + .0033005_r8*rxt(k,637) &
                      *y(k,163)
         mat(k,968) = .1749305_r8*rxt(k,642)*y(k,155) + .0554005_r8*rxt(k,645) &
                      *y(k,163)
         mat(k,2580) = .079_r8*rxt(k,628)*y(k,220) + .0059005_r8*rxt(k,633)*y(k,222) &
                      + .0057005_r8*rxt(k,636)*y(k,238) + .0143005_r8*rxt(k,640) &
                      *y(k,239) + .0332005_r8*rxt(k,644)*y(k,245) &
                      + .0073005_r8*rxt(k,650)*y(k,255) + .011_r8*rxt(k,654)*y(k,258)
         mat(k,1755) = .1749305_r8*rxt(k,626)*y(k,6) + .0590245_r8*rxt(k,634)*y(k,127) &
                      + .1749305_r8*rxt(k,642)*y(k,139)
         mat(k,1888) = .1284005_r8*rxt(k,629)*y(k,6) + .0033005_r8*rxt(k,637)*y(k,127) &
                      + .0554005_r8*rxt(k,645)*y(k,139)
         mat(k,55) = .0085005_r8*rxt(k,647)*y(k,248)
         mat(k,61) = .079_r8*rxt(k,628)*y(k,153) + .1284005_r8*rxt(k,627)*y(k,233)
         mat(k,67) = .0059005_r8*rxt(k,633)*y(k,153) + .0443005_r8*rxt(k,632)*y(k,233)
         mat(k,2238) = .1284005_r8*rxt(k,627)*y(k,220) + .0443005_r8*rxt(k,632) &
                      *y(k,222) + .0271005_r8*rxt(k,635)*y(k,238) &
                      + .0076005_r8*rxt(k,639)*y(k,239) + .0554005_r8*rxt(k,643) &
                      *y(k,245) + .2157005_r8*rxt(k,649)*y(k,255) &
                      + .0512005_r8*rxt(k,653)*y(k,258)
         mat(k,73) = .0057005_r8*rxt(k,636)*y(k,153) + .0271005_r8*rxt(k,635)*y(k,233)
         mat(k,79) = .0143005_r8*rxt(k,640)*y(k,153) + .0076005_r8*rxt(k,639)*y(k,233)
         mat(k,87) = .0332005_r8*rxt(k,644)*y(k,153) + .0554005_r8*rxt(k,643)*y(k,233)
         mat(k,2000) = .0085005_r8*rxt(k,647)*y(k,202)
         mat(k,93) = .0073005_r8*rxt(k,650)*y(k,153) + .2157005_r8*rxt(k,649)*y(k,233)
         mat(k,99) = .011_r8*rxt(k,654)*y(k,153) + .0512005_r8*rxt(k,653)*y(k,233)
         mat(k,1056) = .5901905_r8*rxt(k,626)*y(k,155) + .114_r8*rxt(k,629)*y(k,163)
         mat(k,917) = .0250245_r8*rxt(k,634)*y(k,155)
         mat(k,969) = .5901905_r8*rxt(k,642)*y(k,155) + .1278005_r8*rxt(k,645) &
                      *y(k,163)
         mat(k,2581) = .1254005_r8*rxt(k,628)*y(k,220) + .0536005_r8*rxt(k,633) &
                      *y(k,222) + .0623005_r8*rxt(k,636)*y(k,238) &
                      + .0166005_r8*rxt(k,640)*y(k,239) + .130_r8*rxt(k,644)*y(k,245) &
                      + .238_r8*rxt(k,650)*y(k,255) + .1185005_r8*rxt(k,654)*y(k,258)
         mat(k,1756) = .5901905_r8*rxt(k,626)*y(k,6) + .0250245_r8*rxt(k,634)*y(k,127) &
                      + .5901905_r8*rxt(k,642)*y(k,139)
         mat(k,1889) = .114_r8*rxt(k,629)*y(k,6) + .1278005_r8*rxt(k,645)*y(k,139)
         mat(k,56) = .0128005_r8*rxt(k,647)*y(k,248)
         mat(k,62) = .1254005_r8*rxt(k,628)*y(k,153) + .114_r8*rxt(k,627)*y(k,233)
         mat(k,68) = .0536005_r8*rxt(k,633)*y(k,153) + .1621005_r8*rxt(k,632)*y(k,233)
         mat(k,2239) = .114_r8*rxt(k,627)*y(k,220) + .1621005_r8*rxt(k,632)*y(k,222) &
                      + .0474005_r8*rxt(k,635)*y(k,238) + .0113005_r8*rxt(k,639) &
                      *y(k,239) + .1278005_r8*rxt(k,643)*y(k,245) &
                      + .0738005_r8*rxt(k,649)*y(k,255) + .1598005_r8*rxt(k,653) &
                      *y(k,258)
         mat(k,74) = .0623005_r8*rxt(k,636)*y(k,153) + .0474005_r8*rxt(k,635)*y(k,233)
         mat(k,80) = .0166005_r8*rxt(k,640)*y(k,153) + .0113005_r8*rxt(k,639)*y(k,233)
         mat(k,88) = .130_r8*rxt(k,644)*y(k,153) + .1278005_r8*rxt(k,643)*y(k,233)
         mat(k,2001) = .0128005_r8*rxt(k,647)*y(k,202)
         mat(k,94) = .238_r8*rxt(k,650)*y(k,153) + .0738005_r8*rxt(k,649)*y(k,233)
         mat(k,100) = .1185005_r8*rxt(k,654)*y(k,153) + .1598005_r8*rxt(k,653) &
                      *y(k,233)
         mat(k,57) = -(rxt(k,647)*y(k,248))
         mat(k,2002) = -rxt(k,647)*y(k,202)
         mat(k,223) = .100_r8*rxt(k,523)*y(k,248)
         mat(k,241) = .230_r8*rxt(k,525)*y(k,248)
         mat(k,2025) = .100_r8*rxt(k,523)*y(k,210) + .230_r8*rxt(k,525)*y(k,212)
         mat(k,746) = -(rxt(k,547)*y(k,248))
         mat(k,2097) = -rxt(k,547)*y(k,204)
         mat(k,2281) = rxt(k,545)*y(k,252)
         mat(k,1213) = rxt(k,545)*y(k,233)
         mat(k,721) = -(rxt(k,548)*y(k,248))
         mat(k,2094) = -rxt(k,548)*y(k,205)
         mat(k,2611) = .200_r8*rxt(k,541)*y(k,246) + .200_r8*rxt(k,551)*y(k,253)
         mat(k,1607) = .500_r8*rxt(k,539)*y(k,246)
         mat(k,1163) = .200_r8*rxt(k,541)*y(k,153) + .500_r8*rxt(k,539)*y(k,228)
         mat(k,1086) = .200_r8*rxt(k,551)*y(k,153)
         mat(k,559) = -(rxt(k,552)*y(k,248))
         mat(k,2074) = -rxt(k,552)*y(k,206)
         mat(k,2270) = rxt(k,550)*y(k,253)
         mat(k,1085) = rxt(k,550)*y(k,233)
         mat(k,1098) = -(rxt(k,553)*y(k,155) + rxt(k,554)*y(k,248))
         mat(k,1774) = -rxt(k,553)*y(k,207)
         mat(k,2126) = -rxt(k,554)*y(k,207)
         mat(k,1066) = .330_r8*rxt(k,534)*y(k,163)
         mat(k,979) = .330_r8*rxt(k,537)*y(k,163)
         mat(k,2630) = .800_r8*rxt(k,541)*y(k,246) + .800_r8*rxt(k,551)*y(k,253)
         mat(k,1774) = mat(k,1774) + rxt(k,542)*y(k,246)
         mat(k,1908) = .330_r8*rxt(k,534)*y(k,6) + .330_r8*rxt(k,537)*y(k,139)
         mat(k,722) = rxt(k,548)*y(k,248)
         mat(k,1616) = .500_r8*rxt(k,539)*y(k,246) + rxt(k,549)*y(k,253)
         mat(k,1165) = .800_r8*rxt(k,541)*y(k,153) + rxt(k,542)*y(k,155) &
                      + .500_r8*rxt(k,539)*y(k,228)
         mat(k,2126) = mat(k,2126) + rxt(k,548)*y(k,205)
         mat(k,1089) = .800_r8*rxt(k,551)*y(k,153) + rxt(k,549)*y(k,228)
         mat(k,1180) = -(rxt(k,555)*y(k,248))
         mat(k,2132) = -rxt(k,555)*y(k,208)
         mat(k,1069) = .300_r8*rxt(k,534)*y(k,163)
         mat(k,982) = .300_r8*rxt(k,537)*y(k,163)
         mat(k,2635) = .900_r8*rxt(k,546)*y(k,252)
         mat(k,1912) = .300_r8*rxt(k,534)*y(k,6) + .300_r8*rxt(k,537)*y(k,139)
         mat(k,1621) = rxt(k,544)*y(k,252)
         mat(k,1217) = .900_r8*rxt(k,546)*y(k,153) + rxt(k,544)*y(k,228)
         mat(k,701) = -(rxt(k,522)*y(k,248))
         mat(k,2091) = -rxt(k,522)*y(k,209)
         mat(k,2276) = rxt(k,520)*y(k,254)
         mat(k,820) = rxt(k,520)*y(k,233)
         mat(k,221) = -(rxt(k,523)*y(k,248))
         mat(k,2023) = -rxt(k,523)*y(k,210)
         mat(k,237) = -(rxt(k,489)*y(k,248))
         mat(k,2026) = -rxt(k,489)*y(k,211)
         mat(k,2248) = rxt(k,486)*y(k,256)
         mat(k,1284) = rxt(k,486)*y(k,233)
         mat(k,242) = -(rxt(k,525)*y(k,248))
         mat(k,2027) = -rxt(k,525)*y(k,212)
         mat(k,784) = -(rxt(k,528)*y(k,248))
         mat(k,2101) = -rxt(k,528)*y(k,213)
         mat(k,2284) = rxt(k,526)*y(k,257)
         mat(k,837) = rxt(k,526)*y(k,233)
         mat(k,250) = -(rxt(k,531)*y(k,248))
         mat(k,2028) = -rxt(k,531)*y(k,214)
         mat(k,243) = .150_r8*rxt(k,525)*y(k,248)
         mat(k,2028) = mat(k,2028) + .150_r8*rxt(k,525)*y(k,212)
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
         mat(k,506) = -(rxt(k,532)*y(k,248))
         mat(k,2067) = -rxt(k,532)*y(k,215)
         mat(k,2264) = rxt(k,529)*y(k,259)
         mat(k,583) = rxt(k,529)*y(k,233)
         mat(k,610) = -(rxt(k,490)*y(k,233) + rxt(k,491)*y(k,153) + rxt(k,519) &
                      *y(k,154))
         mat(k,2273) = -rxt(k,490)*y(k,218)
         mat(k,2606) = -rxt(k,491)*y(k,218)
         mat(k,2531) = -rxt(k,519)*y(k,218)
         mat(k,281) = rxt(k,496)*y(k,248)
         mat(k,2080) = rxt(k,496)*y(k,24)
         mat(k,998) = -(rxt(k,451)*y(k,233) + (rxt(k,452) + rxt(k,453)) * y(k,153))
         mat(k,2297) = -rxt(k,451)*y(k,219)
         mat(k,2626) = -(rxt(k,452) + rxt(k,453)) * y(k,219)
         mat(k,739) = rxt(k,454)*y(k,248)
         mat(k,278) = rxt(k,455)*y(k,248)
         mat(k,2118) = rxt(k,454)*y(k,2) + rxt(k,455)*y(k,15)
         mat(k,63) = -(rxt(k,627)*y(k,233) + rxt(k,628)*y(k,153))
         mat(k,2240) = -rxt(k,627)*y(k,220)
         mat(k,2582) = -rxt(k,628)*y(k,220)
         mat(k,1057) = rxt(k,630)*y(k,248)
         mat(k,2003) = rxt(k,630)*y(k,6)
         mat(k,568) = -(rxt(k,493)*y(k,233) + rxt(k,494)*y(k,153))
         mat(k,2271) = -rxt(k,493)*y(k,221)
         mat(k,2601) = -rxt(k,494)*y(k,221)
         mat(k,191) = .350_r8*rxt(k,492)*y(k,248)
         mat(k,496) = rxt(k,495)*y(k,248)
         mat(k,2075) = .350_r8*rxt(k,492)*y(k,7) + rxt(k,495)*y(k,8)
         mat(k,69) = -(rxt(k,632)*y(k,233) + rxt(k,633)*y(k,153))
         mat(k,2241) = -rxt(k,632)*y(k,222)
         mat(k,2583) = -rxt(k,633)*y(k,222)
         mat(k,187) = rxt(k,631)*y(k,248)
         mat(k,2004) = rxt(k,631)*y(k,7)
         mat(k,514) = -(rxt(k,497)*y(k,233) + rxt(k,499)*y(k,153))
         mat(k,2265) = -rxt(k,497)*y(k,223)
         mat(k,2596) = -rxt(k,499)*y(k,223)
         mat(k,383) = rxt(k,498)*y(k,248)
         mat(k,224) = .070_r8*rxt(k,523)*y(k,248)
         mat(k,244) = .060_r8*rxt(k,525)*y(k,248)
         mat(k,2068) = rxt(k,498)*y(k,25) + .070_r8*rxt(k,523)*y(k,210) &
                      + .060_r8*rxt(k,525)*y(k,212)
         mat(k,878) = -(4._r8*rxt(k,373)*y(k,224) + rxt(k,374)*y(k,228) + rxt(k,375) &
                      *y(k,233) + rxt(k,376)*y(k,153))
         mat(k,1611) = -rxt(k,374)*y(k,224)
         mat(k,2293) = -rxt(k,375)*y(k,224)
         mat(k,2622) = -rxt(k,376)*y(k,224)
         mat(k,388) = .500_r8*rxt(k,378)*y(k,248)
         mat(k,328) = rxt(k,379)*y(k,70) + rxt(k,380)*y(k,248)
         mat(k,1844) = rxt(k,379)*y(k,32)
         mat(k,2111) = .500_r8*rxt(k,378)*y(k,31) + rxt(k,380)*y(k,32)
         mat(k,936) = -(rxt(k,402)*y(k,228) + rxt(k,403)*y(k,233) + rxt(k,404) &
                      *y(k,153))
         mat(k,1613) = -rxt(k,402)*y(k,225)
         mat(k,2295) = -rxt(k,403)*y(k,225)
         mat(k,2624) = -rxt(k,404)*y(k,225)
         mat(k,489) = rxt(k,405)*y(k,248)
         mat(k,334) = rxt(k,409)*y(k,70) + rxt(k,406)*y(k,248)
         mat(k,1846) = rxt(k,409)*y(k,35)
         mat(k,2114) = rxt(k,405)*y(k,34) + rxt(k,406)*y(k,35)
         mat(k,729) = -(rxt(k,500)*y(k,233) + rxt(k,501)*y(k,153))
         mat(k,2279) = -rxt(k,500)*y(k,226)
         mat(k,2612) = -rxt(k,501)*y(k,226)
         mat(k,302) = rxt(k,502)*y(k,248)
         mat(k,2612) = mat(k,2612) + rxt(k,491)*y(k,218)
         mat(k,1896) = rxt(k,517)*y(k,171)
         mat(k,546) = rxt(k,517)*y(k,163)
         mat(k,611) = rxt(k,491)*y(k,153) + .400_r8*rxt(k,490)*y(k,233)
         mat(k,2279) = mat(k,2279) + .400_r8*rxt(k,490)*y(k,218)
         mat(k,2095) = rxt(k,502)*y(k,36)
         mat(k,1491) = -(4._r8*rxt(k,384)*y(k,227) + rxt(k,385)*y(k,228) + rxt(k,386) &
                      *y(k,233) + rxt(k,387)*y(k,153) + rxt(k,398)*y(k,154) + rxt(k,426) &
                      *y(k,240) + rxt(k,459)*y(k,235) + rxt(k,464)*y(k,236) + rxt(k,473) &
                      *y(k,237) + rxt(k,484)*y(k,256))
         mat(k,1637) = -rxt(k,385)*y(k,227)
         mat(k,2322) = -rxt(k,386)*y(k,227)
         mat(k,2652) = -rxt(k,387)*y(k,227)
         mat(k,2550) = -rxt(k,398)*y(k,227)
         mat(k,1418) = -rxt(k,426)*y(k,227)
         mat(k,1363) = -rxt(k,459)*y(k,227)
         mat(k,1396) = -rxt(k,464)*y(k,227)
         mat(k,1315) = -rxt(k,473)*y(k,227)
         mat(k,1293) = -rxt(k,484)*y(k,227)
         mat(k,1074) = .060_r8*rxt(k,534)*y(k,163)
         mat(k,1192) = rxt(k,381)*y(k,155) + rxt(k,382)*y(k,248)
         mat(k,1340) = rxt(k,407)*y(k,155) + rxt(k,408)*y(k,248)
         mat(k,676) = .500_r8*rxt(k,389)*y(k,248)
         mat(k,928) = .080_r8*rxt(k,479)*y(k,163)
         mat(k,1331) = .100_r8*rxt(k,432)*y(k,163)
         mat(k,986) = .060_r8*rxt(k,537)*y(k,163)
         mat(k,1439) = .280_r8*rxt(k,446)*y(k,163)
         mat(k,2652) = mat(k,2652) + .530_r8*rxt(k,430)*y(k,240) + rxt(k,439)*y(k,242) &
                      + rxt(k,442)*y(k,244) + rxt(k,417)*y(k,251)
         mat(k,1798) = rxt(k,381)*y(k,54) + rxt(k,407)*y(k,58) + .530_r8*rxt(k,429) &
                      *y(k,240) + rxt(k,440)*y(k,242)
         mat(k,1927) = .060_r8*rxt(k,534)*y(k,6) + .080_r8*rxt(k,479)*y(k,127) &
                      + .100_r8*rxt(k,432)*y(k,134) + .060_r8*rxt(k,537)*y(k,139) &
                      + .280_r8*rxt(k,446)*y(k,140)
         mat(k,1183) = .650_r8*rxt(k,555)*y(k,248)
         mat(k,1491) = mat(k,1491) + .530_r8*rxt(k,426)*y(k,240)
         mat(k,1637) = mat(k,1637) + .260_r8*rxt(k,427)*y(k,240) + rxt(k,436)*y(k,242) &
                      + .300_r8*rxt(k,415)*y(k,251)
         mat(k,2322) = mat(k,2322) + .450_r8*rxt(k,437)*y(k,242) + .200_r8*rxt(k,441) &
                      *y(k,244) + .150_r8*rxt(k,416)*y(k,251)
         mat(k,1418) = mat(k,1418) + .530_r8*rxt(k,430)*y(k,153) + .530_r8*rxt(k,429) &
                      *y(k,155) + .530_r8*rxt(k,426)*y(k,227) + .260_r8*rxt(k,427) &
                      *y(k,228)
         mat(k,1460) = rxt(k,439)*y(k,153) + rxt(k,440)*y(k,155) + rxt(k,436)*y(k,228) &
                      + .450_r8*rxt(k,437)*y(k,233) + 4.000_r8*rxt(k,438)*y(k,242)
         mat(k,770) = rxt(k,442)*y(k,153) + .200_r8*rxt(k,441)*y(k,233)
         mat(k,2150) = rxt(k,382)*y(k,54) + rxt(k,408)*y(k,58) + .500_r8*rxt(k,389) &
                      *y(k,60) + .650_r8*rxt(k,555)*y(k,208)
         mat(k,1276) = rxt(k,417)*y(k,153) + .300_r8*rxt(k,415)*y(k,228) &
                      + .150_r8*rxt(k,416)*y(k,233)
         mat(k,1639) = -(rxt(k,213)*y(k,74) + (rxt(k,332) + rxt(k,333)) * y(k,68) &
                      + (4._r8*rxt(k,352) + 4._r8*rxt(k,353)) * y(k,228) + rxt(k,354) &
                      *y(k,233) + rxt(k,355)*y(k,153) + rxt(k,374)*y(k,224) + rxt(k,385) &
                      *y(k,227) + rxt(k,402)*y(k,225) + rxt(k,415)*y(k,251) + rxt(k,427) &
                      *y(k,240) + rxt(k,436)*y(k,242) + rxt(k,460)*y(k,235) + rxt(k,465) &
                      *y(k,236) + rxt(k,474)*y(k,237) + rxt(k,485)*y(k,256) + rxt(k,539) &
                      *y(k,246) + rxt(k,544)*y(k,252) + rxt(k,549)*y(k,253))
         mat(k,2358) = -rxt(k,213)*y(k,228)
         mat(k,1115) = -(rxt(k,332) + rxt(k,333)) * y(k,228)
         mat(k,2328) = -rxt(k,354)*y(k,228)
         mat(k,2656) = -rxt(k,355)*y(k,228)
         mat(k,880) = -rxt(k,374)*y(k,228)
         mat(k,1493) = -rxt(k,385)*y(k,228)
         mat(k,939) = -rxt(k,402)*y(k,228)
         mat(k,1277) = -rxt(k,415)*y(k,228)
         mat(k,1419) = -rxt(k,427)*y(k,228)
         mat(k,1461) = -rxt(k,436)*y(k,228)
         mat(k,1364) = -rxt(k,460)*y(k,228)
         mat(k,1397) = -rxt(k,465)*y(k,228)
         mat(k,1316) = -rxt(k,474)*y(k,228)
         mat(k,1294) = -rxt(k,485)*y(k,228)
         mat(k,1170) = -rxt(k,539)*y(k,228)
         mat(k,1223) = -rxt(k,544)*y(k,228)
         mat(k,1091) = -rxt(k,549)*y(k,228)
         mat(k,1136) = .280_r8*rxt(k,401)*y(k,163)
         mat(k,777) = rxt(k,388)*y(k,248)
         mat(k,477) = .700_r8*rxt(k,357)*y(k,248)
         mat(k,1568) = rxt(k,205)*y(k,70) + rxt(k,306)*y(k,89) + rxt(k,364)*y(k,247) &
                      + rxt(k,358)*y(k,248)
         mat(k,1863) = rxt(k,205)*y(k,64)
         mat(k,950) = rxt(k,306)*y(k,64)
         mat(k,929) = .050_r8*rxt(k,479)*y(k,163)
         mat(k,2656) = mat(k,2656) + rxt(k,387)*y(k,227) + .830_r8*rxt(k,505)*y(k,229) &
                      + .170_r8*rxt(k,511)*y(k,243)
         mat(k,1930) = .280_r8*rxt(k,401)*y(k,33) + .050_r8*rxt(k,479)*y(k,127)
         mat(k,1493) = mat(k,1493) + rxt(k,387)*y(k,153) + 4.000_r8*rxt(k,384) &
                      *y(k,227) + .900_r8*rxt(k,385)*y(k,228) + .450_r8*rxt(k,386) &
                      *y(k,233) + rxt(k,459)*y(k,235) + rxt(k,464)*y(k,236) &
                      + rxt(k,473)*y(k,237) + rxt(k,426)*y(k,240) + rxt(k,435) &
                      *y(k,242) + rxt(k,484)*y(k,256)
         mat(k,1639) = mat(k,1639) + .900_r8*rxt(k,385)*y(k,227)
         mat(k,853) = .830_r8*rxt(k,505)*y(k,153) + .330_r8*rxt(k,504)*y(k,233)
         mat(k,2328) = mat(k,2328) + .450_r8*rxt(k,386)*y(k,227) + .330_r8*rxt(k,504) &
                      *y(k,229) + .070_r8*rxt(k,510)*y(k,243)
         mat(k,1364) = mat(k,1364) + rxt(k,459)*y(k,227)
         mat(k,1397) = mat(k,1397) + rxt(k,464)*y(k,227)
         mat(k,1316) = mat(k,1316) + rxt(k,473)*y(k,227)
         mat(k,1419) = mat(k,1419) + rxt(k,426)*y(k,227)
         mat(k,1461) = mat(k,1461) + rxt(k,435)*y(k,227)
         mat(k,960) = .170_r8*rxt(k,511)*y(k,153) + .070_r8*rxt(k,510)*y(k,233)
         mat(k,1975) = rxt(k,364)*y(k,64)
         mat(k,2157) = rxt(k,388)*y(k,59) + .700_r8*rxt(k,357)*y(k,63) + rxt(k,358) &
                      *y(k,64)
         mat(k,1294) = mat(k,1294) + rxt(k,484)*y(k,227)
         mat(k,850) = -(rxt(k,504)*y(k,233) + rxt(k,505)*y(k,153) + rxt(k,506) &
                      *y(k,154))
         mat(k,2290) = -rxt(k,504)*y(k,229)
         mat(k,2619) = -rxt(k,505)*y(k,229)
         mat(k,2538) = -rxt(k,506)*y(k,229)
         mat(k,649) = -((rxt(k,423) + rxt(k,424)) * y(k,153))
         mat(k,2607) = -(rxt(k,423) + rxt(k,424)) * y(k,230)
         mat(k,412) = rxt(k,422)*y(k,248)
         mat(k,2085) = rxt(k,422)*y(k,16)
         mat(k,2592) = .750_r8*rxt(k,391)*y(k,232)
         mat(k,796) = .750_r8*rxt(k,391)*y(k,153)
         mat(k,797) = -(rxt(k,390)*y(k,233) + rxt(k,391)*y(k,153))
         mat(k,2285) = -rxt(k,390)*y(k,232)
         mat(k,2615) = -rxt(k,391)*y(k,232)
         mat(k,642) = rxt(k,397)*y(k,248)
         mat(k,2102) = rxt(k,397)*y(k,28)
         mat(k,2340) = -((rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,92) + rxt(k,168) &
                      *y(k,162) + rxt(k,169)*y(k,163) + rxt(k,173)*y(k,248) &
                      + 4._r8*rxt(k,178)*y(k,233) + rxt(k,188)*y(k,155) + rxt(k,193) &
                      *y(k,153) + rxt(k,198)*y(k,154) + (rxt(k,208) + rxt(k,209) &
                      ) * y(k,70) + rxt(k,217)*y(k,74) + rxt(k,244)*y(k,17) + rxt(k,251) &
                      *y(k,21) + rxt(k,273)*y(k,115) + rxt(k,288)*y(k,125) + rxt(k,334) &
                      *y(k,68) + rxt(k,348)*y(k,51) + rxt(k,354)*y(k,228) + rxt(k,361) &
                      *y(k,234) + rxt(k,375)*y(k,224) + rxt(k,386)*y(k,227) + rxt(k,390) &
                      *y(k,232) + rxt(k,403)*y(k,225) + rxt(k,412)*y(k,250) + rxt(k,416) &
                      *y(k,251) + rxt(k,428)*y(k,240) + rxt(k,437)*y(k,242) + rxt(k,441) &
                      *y(k,244) + rxt(k,451)*y(k,219) + rxt(k,461)*y(k,235) + rxt(k,466) &
                      *y(k,236) + rxt(k,475)*y(k,237) + rxt(k,486)*y(k,256) + rxt(k,490) &
                      *y(k,218) + rxt(k,493)*y(k,221) + rxt(k,497)*y(k,223) + rxt(k,500) &
                      *y(k,226) + rxt(k,504)*y(k,229) + rxt(k,507)*y(k,241) + rxt(k,510) &
                      *y(k,243) + rxt(k,513)*y(k,249) + rxt(k,520)*y(k,254) + rxt(k,526) &
                      *y(k,257) + rxt(k,529)*y(k,259) + rxt(k,540)*y(k,246) + rxt(k,545) &
                      *y(k,252) + rxt(k,550)*y(k,253))
         mat(k,2510) = -(rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,233)
         mat(k,2432) = -rxt(k,168)*y(k,233)
         mat(k,1941) = -rxt(k,169)*y(k,233)
         mat(k,2169) = -rxt(k,173)*y(k,233)
         mat(k,1815) = -rxt(k,188)*y(k,233)
         mat(k,2668) = -rxt(k,193)*y(k,233)
         mat(k,2567) = -rxt(k,198)*y(k,233)
         mat(k,1875) = -(rxt(k,208) + rxt(k,209)) * y(k,233)
         mat(k,2369) = -rxt(k,217)*y(k,233)
         mat(k,2487) = -rxt(k,244)*y(k,233)
         mat(k,2459) = -rxt(k,251)*y(k,233)
         mat(k,2397) = -rxt(k,273)*y(k,233)
         mat(k,2225) = -rxt(k,288)*y(k,233)
         mat(k,1120) = -rxt(k,334)*y(k,233)
         mat(k,1698) = -rxt(k,348)*y(k,233)
         mat(k,1649) = -rxt(k,354)*y(k,233)
         mat(k,523) = -rxt(k,361)*y(k,233)
         mat(k,883) = -rxt(k,375)*y(k,233)
         mat(k,1499) = -rxt(k,386)*y(k,233)
         mat(k,801) = -rxt(k,390)*y(k,233)
         mat(k,942) = -rxt(k,403)*y(k,233)
         mat(k,864) = -rxt(k,412)*y(k,233)
         mat(k,1280) = -rxt(k,416)*y(k,233)
         mat(k,1424) = -rxt(k,428)*y(k,233)
         mat(k,1467) = -rxt(k,437)*y(k,233)
         mat(k,772) = -rxt(k,441)*y(k,233)
         mat(k,1005) = -rxt(k,451)*y(k,233)
         mat(k,1370) = -rxt(k,461)*y(k,233)
         mat(k,1403) = -rxt(k,466)*y(k,233)
         mat(k,1321) = -rxt(k,475)*y(k,233)
         mat(k,1298) = -rxt(k,486)*y(k,233)
         mat(k,613) = -rxt(k,490)*y(k,233)
         mat(k,572) = -rxt(k,493)*y(k,233)
         mat(k,517) = -rxt(k,497)*y(k,233)
         mat(k,732) = -rxt(k,500)*y(k,233)
         mat(k,855) = -rxt(k,504)*y(k,233)
         mat(k,815) = -rxt(k,507)*y(k,233)
         mat(k,962) = -rxt(k,510)*y(k,233)
         mat(k,530) = -rxt(k,513)*y(k,233)
         mat(k,830) = -rxt(k,520)*y(k,233)
         mat(k,847) = -rxt(k,526)*y(k,233)
         mat(k,588) = -rxt(k,529)*y(k,233)
         mat(k,1175) = -rxt(k,540)*y(k,233)
         mat(k,1227) = -rxt(k,545)*y(k,233)
         mat(k,1095) = -rxt(k,550)*y(k,233)
         mat(k,1080) = .570_r8*rxt(k,534)*y(k,163)
         mat(k,193) = .650_r8*rxt(k,492)*y(k,248)
         mat(k,2487) = mat(k,2487) + rxt(k,243)*y(k,51)
         mat(k,2459) = mat(k,2459) + rxt(k,258)*y(k,248)
         mat(k,326) = .350_r8*rxt(k,370)*y(k,248)
         mat(k,648) = .130_r8*rxt(k,372)*y(k,163)
         mat(k,299) = rxt(k,377)*y(k,248)
         mat(k,1141) = .280_r8*rxt(k,401)*y(k,163)
         mat(k,1698) = mat(k,1698) + rxt(k,243)*y(k,17) + rxt(k,204)*y(k,70) &
                      + rxt(k,349)*y(k,155) + rxt(k,350)*y(k,162)
         mat(k,671) = rxt(k,321)*y(k,70) + rxt(k,322)*y(k,248)
         mat(k,444) = rxt(k,324)*y(k,70) + rxt(k,325)*y(k,248)
         mat(k,109) = rxt(k,383)*y(k,248)
         mat(k,435) = rxt(k,326)*y(k,70) + rxt(k,327)*y(k,248)
         mat(k,871) = rxt(k,356)*y(k,248)
         mat(k,1576) = rxt(k,365)*y(k,247)
         mat(k,1120) = mat(k,1120) + rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155) + ( &
                      + 2.000_r8*rxt(k,332)+rxt(k,333))*y(k,228)
         mat(k,1875) = mat(k,1875) + rxt(k,204)*y(k,51) + rxt(k,321)*y(k,52) &
                      + rxt(k,324)*y(k,55) + rxt(k,326)*y(k,61) + rxt(k,207)*y(k,95)
         mat(k,2369) = mat(k,2369) + rxt(k,213)*y(k,228) + rxt(k,224)*y(k,248)
         mat(k,1202) = rxt(k,368)*y(k,248)
         mat(k,232) = .730_r8*rxt(k,503)*y(k,248)
         mat(k,1034) = .500_r8*rxt(k,575)*y(k,248)
         mat(k,1210) = rxt(k,394)*y(k,248)
         mat(k,1050) = rxt(k,395)*y(k,248)
         mat(k,718) = rxt(k,207)*y(k,70) + rxt(k,163)*y(k,162) + rxt(k,172)*y(k,248)
         mat(k,208) = rxt(k,359)*y(k,248)
         mat(k,1023) = rxt(k,360)*y(k,248)
         mat(k,1246) = rxt(k,425)*y(k,248)
         mat(k,1255) = rxt(k,410)*y(k,248)
         mat(k,2225) = mat(k,2225) + rxt(k,294)*y(k,248)
         mat(k,934) = .370_r8*rxt(k,479)*y(k,163)
         mat(k,697) = .300_r8*rxt(k,470)*y(k,248)
         mat(k,664) = rxt(k,471)*y(k,248)
         mat(k,450) = rxt(k,478)*y(k,248)
         mat(k,1335) = .140_r8*rxt(k,432)*y(k,163)
         mat(k,358) = .200_r8*rxt(k,434)*y(k,248)
         mat(k,688) = .500_r8*rxt(k,445)*y(k,248)
         mat(k,991) = .570_r8*rxt(k,537)*y(k,163)
         mat(k,1447) = .280_r8*rxt(k,446)*y(k,163)
         mat(k,463) = rxt(k,482)*y(k,248)
         mat(k,1158) = rxt(k,483)*y(k,248)
         mat(k,2668) = mat(k,2668) + rxt(k,336)*y(k,68) + rxt(k,452)*y(k,219) &
                      + rxt(k,494)*y(k,221) + rxt(k,499)*y(k,223) + rxt(k,376) &
                      *y(k,224) + rxt(k,404)*y(k,225) + rxt(k,355)*y(k,228) &
                      + .170_r8*rxt(k,505)*y(k,229) + rxt(k,423)*y(k,230) &
                      + .250_r8*rxt(k,391)*y(k,232) + rxt(k,363)*y(k,234) &
                      + .920_r8*rxt(k,462)*y(k,235) + .920_r8*rxt(k,468)*y(k,236) &
                      + rxt(k,476)*y(k,237) + .470_r8*rxt(k,430)*y(k,240) &
                      + .400_r8*rxt(k,508)*y(k,241) + .830_r8*rxt(k,511)*y(k,243) &
                      + rxt(k,514)*y(k,249) + rxt(k,413)*y(k,250) + .900_r8*rxt(k,546) &
                      *y(k,252) + .800_r8*rxt(k,551)*y(k,253) + rxt(k,521)*y(k,254) &
                      + rxt(k,487)*y(k,256) + rxt(k,527)*y(k,257) + rxt(k,530) &
                      *y(k,259)
         mat(k,1815) = mat(k,1815) + rxt(k,349)*y(k,51) + rxt(k,337)*y(k,68) &
                      + rxt(k,463)*y(k,235) + rxt(k,469)*y(k,236) + rxt(k,477) &
                      *y(k,237) + .470_r8*rxt(k,429)*y(k,240) + rxt(k,191)*y(k,248) &
                      + rxt(k,488)*y(k,256)
         mat(k,2432) = mat(k,2432) + rxt(k,350)*y(k,51) + rxt(k,163)*y(k,95)
         mat(k,1941) = mat(k,1941) + .570_r8*rxt(k,534)*y(k,6) + .130_r8*rxt(k,372) &
                      *y(k,28) + .280_r8*rxt(k,401)*y(k,33) + .370_r8*rxt(k,479) &
                      *y(k,127) + .140_r8*rxt(k,432)*y(k,134) + .570_r8*rxt(k,537) &
                      *y(k,139) + .280_r8*rxt(k,446)*y(k,140) + rxt(k,175)*y(k,248)
         mat(k,202) = .800_r8*rxt(k,515)*y(k,248)
         mat(k,1042) = rxt(k,565)*y(k,248)
         mat(k,1187) = .200_r8*rxt(k,555)*y(k,248)
         mat(k,227) = .280_r8*rxt(k,523)*y(k,248)
         mat(k,249) = .380_r8*rxt(k,525)*y(k,248)
         mat(k,254) = .630_r8*rxt(k,531)*y(k,248)
         mat(k,1005) = mat(k,1005) + rxt(k,452)*y(k,153)
         mat(k,572) = mat(k,572) + rxt(k,494)*y(k,153)
         mat(k,517) = mat(k,517) + rxt(k,499)*y(k,153)
         mat(k,883) = mat(k,883) + rxt(k,376)*y(k,153) + 2.400_r8*rxt(k,373)*y(k,224) &
                      + rxt(k,374)*y(k,228)
         mat(k,942) = mat(k,942) + rxt(k,404)*y(k,153) + rxt(k,402)*y(k,228)
         mat(k,1499) = mat(k,1499) + .900_r8*rxt(k,385)*y(k,228) + rxt(k,459)*y(k,235) &
                      + rxt(k,464)*y(k,236) + rxt(k,473)*y(k,237) + .470_r8*rxt(k,426) &
                      *y(k,240) + rxt(k,484)*y(k,256)
         mat(k,1649) = mat(k,1649) + (2.000_r8*rxt(k,332)+rxt(k,333))*y(k,68) &
                      + rxt(k,213)*y(k,74) + rxt(k,355)*y(k,153) + rxt(k,374)*y(k,224) &
                      + rxt(k,402)*y(k,225) + .900_r8*rxt(k,385)*y(k,227) &
                      + 4.000_r8*rxt(k,352)*y(k,228) + rxt(k,460)*y(k,235) &
                      + rxt(k,465)*y(k,236) + 1.200_r8*rxt(k,474)*y(k,237) &
                      + .730_r8*rxt(k,427)*y(k,240) + rxt(k,436)*y(k,242) &
                      + .500_r8*rxt(k,539)*y(k,246) + .300_r8*rxt(k,415)*y(k,251) &
                      + rxt(k,544)*y(k,252) + rxt(k,549)*y(k,253) + .800_r8*rxt(k,485) &
                      *y(k,256)
         mat(k,855) = mat(k,855) + .170_r8*rxt(k,505)*y(k,153) + .070_r8*rxt(k,504) &
                      *y(k,233)
         mat(k,654) = rxt(k,423)*y(k,153)
         mat(k,801) = mat(k,801) + .250_r8*rxt(k,391)*y(k,153)
         mat(k,2340) = mat(k,2340) + .070_r8*rxt(k,504)*y(k,229) + .160_r8*rxt(k,507) &
                      *y(k,241) + .330_r8*rxt(k,510)*y(k,243)
         mat(k,523) = mat(k,523) + rxt(k,363)*y(k,153)
         mat(k,1370) = mat(k,1370) + .920_r8*rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155) &
                      + rxt(k,459)*y(k,227) + rxt(k,460)*y(k,228)
         mat(k,1403) = mat(k,1403) + .920_r8*rxt(k,468)*y(k,153) + rxt(k,469)*y(k,155) &
                      + rxt(k,464)*y(k,227) + rxt(k,465)*y(k,228)
         mat(k,1321) = mat(k,1321) + rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155) &
                      + rxt(k,473)*y(k,227) + 1.200_r8*rxt(k,474)*y(k,228)
         mat(k,1424) = mat(k,1424) + .470_r8*rxt(k,430)*y(k,153) + .470_r8*rxt(k,429) &
                      *y(k,155) + .470_r8*rxt(k,426)*y(k,227) + .730_r8*rxt(k,427) &
                      *y(k,228)
         mat(k,815) = mat(k,815) + .400_r8*rxt(k,508)*y(k,153) + .160_r8*rxt(k,507) &
                      *y(k,233)
         mat(k,1467) = mat(k,1467) + rxt(k,436)*y(k,228)
         mat(k,962) = mat(k,962) + .830_r8*rxt(k,511)*y(k,153) + .330_r8*rxt(k,510) &
                      *y(k,233)
         mat(k,1175) = mat(k,1175) + .500_r8*rxt(k,539)*y(k,228)
         mat(k,1987) = rxt(k,365)*y(k,64)
         mat(k,2169) = mat(k,2169) + .650_r8*rxt(k,492)*y(k,7) + rxt(k,258)*y(k,21) &
                      + .350_r8*rxt(k,370)*y(k,27) + rxt(k,377)*y(k,30) + rxt(k,322) &
                      *y(k,52) + rxt(k,325)*y(k,55) + rxt(k,383)*y(k,56) + rxt(k,327) &
                      *y(k,61) + rxt(k,356)*y(k,62) + rxt(k,224)*y(k,74) + rxt(k,368) &
                      *y(k,77) + .730_r8*rxt(k,503)*y(k,82) + .500_r8*rxt(k,575) &
                      *y(k,83) + rxt(k,394)*y(k,90) + rxt(k,395)*y(k,91) + rxt(k,172) &
                      *y(k,95) + rxt(k,359)*y(k,102) + rxt(k,360)*y(k,103) &
                      + rxt(k,425)*y(k,111) + rxt(k,410)*y(k,113) + rxt(k,294) &
                      *y(k,125) + .300_r8*rxt(k,470)*y(k,128) + rxt(k,471)*y(k,129) &
                      + rxt(k,478)*y(k,130) + .200_r8*rxt(k,434)*y(k,135) &
                      + .500_r8*rxt(k,445)*y(k,138) + rxt(k,482)*y(k,144) + rxt(k,483) &
                      *y(k,145) + rxt(k,191)*y(k,155) + rxt(k,175)*y(k,163) &
                      + .800_r8*rxt(k,515)*y(k,172) + rxt(k,565)*y(k,181) &
                      + .200_r8*rxt(k,555)*y(k,208) + .280_r8*rxt(k,523)*y(k,210) &
                      + .380_r8*rxt(k,525)*y(k,212) + .630_r8*rxt(k,531)*y(k,214)
         mat(k,530) = mat(k,530) + rxt(k,514)*y(k,153)
         mat(k,864) = mat(k,864) + rxt(k,413)*y(k,153)
         mat(k,1280) = mat(k,1280) + .300_r8*rxt(k,415)*y(k,228)
         mat(k,1227) = mat(k,1227) + .900_r8*rxt(k,546)*y(k,153) + rxt(k,544)*y(k,228)
         mat(k,1095) = mat(k,1095) + .800_r8*rxt(k,551)*y(k,153) + rxt(k,549)*y(k,228)
         mat(k,830) = mat(k,830) + rxt(k,521)*y(k,153)
         mat(k,1298) = mat(k,1298) + rxt(k,487)*y(k,153) + rxt(k,488)*y(k,155) &
                      + rxt(k,484)*y(k,227) + .800_r8*rxt(k,485)*y(k,228)
         mat(k,847) = mat(k,847) + rxt(k,527)*y(k,153)
         mat(k,588) = mat(k,588) + rxt(k,530)*y(k,153)
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
         mat(k,520) = -(rxt(k,361)*y(k,233) + rxt(k,363)*y(k,153))
         mat(k,2266) = -rxt(k,361)*y(k,234)
         mat(k,2597) = -rxt(k,363)*y(k,234)
         mat(k,1679) = rxt(k,348)*y(k,233)
         mat(k,2266) = mat(k,2266) + rxt(k,348)*y(k,51)
         mat(k,1359) = -(rxt(k,459)*y(k,227) + rxt(k,460)*y(k,228) + rxt(k,461) &
                      *y(k,233) + rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155))
         mat(k,1486) = -rxt(k,459)*y(k,235)
         mat(k,1632) = -rxt(k,460)*y(k,235)
         mat(k,2317) = -rxt(k,461)*y(k,235)
         mat(k,2647) = -rxt(k,462)*y(k,235)
         mat(k,1793) = -rxt(k,463)*y(k,235)
         mat(k,925) = .600_r8*rxt(k,480)*y(k,248)
         mat(k,2145) = .600_r8*rxt(k,480)*y(k,127)
         mat(k,1392) = -(rxt(k,464)*y(k,227) + rxt(k,465)*y(k,228) + rxt(k,466) &
                      *y(k,233) + rxt(k,468)*y(k,153) + rxt(k,469)*y(k,155))
         mat(k,1487) = -rxt(k,464)*y(k,236)
         mat(k,1633) = -rxt(k,465)*y(k,236)
         mat(k,2318) = -rxt(k,466)*y(k,236)
         mat(k,2648) = -rxt(k,468)*y(k,236)
         mat(k,1794) = -rxt(k,469)*y(k,236)
         mat(k,926) = .400_r8*rxt(k,480)*y(k,248)
         mat(k,2146) = .400_r8*rxt(k,480)*y(k,127)
         mat(k,1311) = -(rxt(k,473)*y(k,227) + rxt(k,474)*y(k,228) + rxt(k,475) &
                      *y(k,233) + rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155))
         mat(k,1483) = -rxt(k,473)*y(k,237)
         mat(k,1629) = -rxt(k,474)*y(k,237)
         mat(k,2314) = -rxt(k,475)*y(k,237)
         mat(k,2644) = -rxt(k,476)*y(k,237)
         mat(k,1790) = -rxt(k,477)*y(k,237)
         mat(k,923) = rxt(k,472)*y(k,155)
         mat(k,1790) = mat(k,1790) + rxt(k,472)*y(k,127)
         mat(k,75) = -(rxt(k,635)*y(k,233) + rxt(k,636)*y(k,153))
         mat(k,2242) = -rxt(k,635)*y(k,238)
         mat(k,2584) = -rxt(k,636)*y(k,238)
         mat(k,918) = rxt(k,638)*y(k,248)
         mat(k,2005) = rxt(k,638)*y(k,127)
         mat(k,81) = -(rxt(k,639)*y(k,233) + rxt(k,640)*y(k,153))
         mat(k,2243) = -rxt(k,639)*y(k,239)
         mat(k,2585) = -rxt(k,640)*y(k,239)
         mat(k,82) = rxt(k,641)*y(k,248)
         mat(k,2006) = rxt(k,641)*y(k,132)
         mat(k,1416) = -(rxt(k,426)*y(k,227) + rxt(k,427)*y(k,228) + rxt(k,428) &
                      *y(k,233) + rxt(k,429)*y(k,155) + (rxt(k,430) + rxt(k,431) &
                      ) * y(k,153))
         mat(k,1488) = -rxt(k,426)*y(k,240)
         mat(k,1634) = -rxt(k,427)*y(k,240)
         mat(k,2319) = -rxt(k,428)*y(k,240)
         mat(k,1795) = -rxt(k,429)*y(k,240)
         mat(k,2649) = -(rxt(k,430) + rxt(k,431)) * y(k,240)
         mat(k,1329) = .500_r8*rxt(k,433)*y(k,248)
         mat(k,355) = .200_r8*rxt(k,434)*y(k,248)
         mat(k,1436) = rxt(k,447)*y(k,248)
         mat(k,2147) = .500_r8*rxt(k,433)*y(k,134) + .200_r8*rxt(k,434)*y(k,135) &
                      + rxt(k,447)*y(k,140)
         mat(k,812) = -(rxt(k,507)*y(k,233) + rxt(k,508)*y(k,153) + rxt(k,509) &
                      *y(k,154))
         mat(k,2287) = -rxt(k,507)*y(k,241)
         mat(k,2616) = -rxt(k,508)*y(k,241)
         mat(k,2537) = -rxt(k,509)*y(k,241)
         mat(k,1459) = -(rxt(k,435)*y(k,227) + rxt(k,436)*y(k,228) + rxt(k,437) &
                      *y(k,233) + 4._r8*rxt(k,438)*y(k,242) + rxt(k,439)*y(k,153) &
                      + rxt(k,440)*y(k,155) + rxt(k,448)*y(k,154))
         mat(k,1490) = -rxt(k,435)*y(k,242)
         mat(k,1636) = -rxt(k,436)*y(k,242)
         mat(k,2321) = -rxt(k,437)*y(k,242)
         mat(k,2651) = -rxt(k,439)*y(k,242)
         mat(k,1797) = -rxt(k,440)*y(k,242)
         mat(k,2549) = -rxt(k,448)*y(k,242)
         mat(k,1330) = .500_r8*rxt(k,433)*y(k,248)
         mat(k,356) = .500_r8*rxt(k,434)*y(k,248)
         mat(k,2149) = .500_r8*rxt(k,433)*y(k,134) + .500_r8*rxt(k,434)*y(k,135)
         mat(k,956) = -(rxt(k,510)*y(k,233) + rxt(k,511)*y(k,153) + rxt(k,512) &
                      *y(k,154))
         mat(k,2296) = -rxt(k,510)*y(k,243)
         mat(k,2625) = -rxt(k,511)*y(k,243)
         mat(k,2542) = -rxt(k,512)*y(k,243)
         mat(k,768) = -(rxt(k,441)*y(k,233) + rxt(k,442)*y(k,153))
         mat(k,2282) = -rxt(k,441)*y(k,244)
         mat(k,2614) = -rxt(k,442)*y(k,244)
         mat(k,599) = rxt(k,443)*y(k,248)
         mat(k,360) = rxt(k,444)*y(k,248)
         mat(k,2099) = rxt(k,443)*y(k,136) + rxt(k,444)*y(k,137)
         mat(k,89) = -(rxt(k,643)*y(k,233) + rxt(k,644)*y(k,153))
         mat(k,2244) = -rxt(k,643)*y(k,245)
         mat(k,2586) = -rxt(k,644)*y(k,245)
         mat(k,970) = rxt(k,646)*y(k,248)
         mat(k,2008) = rxt(k,646)*y(k,139)
         mat(k,1166) = -(rxt(k,539)*y(k,228) + rxt(k,540)*y(k,233) + rxt(k,541) &
                      *y(k,153) + rxt(k,542)*y(k,155))
         mat(k,1620) = -rxt(k,539)*y(k,246)
         mat(k,2304) = -rxt(k,540)*y(k,246)
         mat(k,2634) = -rxt(k,541)*y(k,246)
         mat(k,1779) = -rxt(k,542)*y(k,246)
         mat(k,1068) = rxt(k,533)*y(k,155)
         mat(k,981) = rxt(k,536)*y(k,155)
         mat(k,1779) = mat(k,1779) + rxt(k,533)*y(k,6) + rxt(k,536)*y(k,139) &
                      + .500_r8*rxt(k,553)*y(k,207)
         mat(k,466) = rxt(k,543)*y(k,248)
         mat(k,1099) = .500_r8*rxt(k,553)*y(k,155)
         mat(k,2131) = rxt(k,543)*y(k,157)
         mat(k,1983) = -(rxt(k,153)*y(k,93) + rxt(k,154)*y(k,260) + (rxt(k,157) &
                      + rxt(k,158)) * y(k,163) + (rxt(k,196) + rxt(k,197)) * y(k,142) &
                      + rxt(k,231)*y(k,37) + rxt(k,232)*y(k,38) + rxt(k,233)*y(k,40) &
                      + rxt(k,234)*y(k,41) + rxt(k,235)*y(k,42) + rxt(k,236)*y(k,43) &
                      + rxt(k,237)*y(k,44) + (rxt(k,238) + rxt(k,239)) * y(k,101) &
                      + rxt(k,262)*y(k,39) + rxt(k,263)*y(k,66) + rxt(k,264)*y(k,94) &
                      + (rxt(k,265) + rxt(k,266)) * y(k,97) + rxt(k,310)*y(k,80) &
                      + rxt(k,311)*y(k,81) + rxt(k,343)*y(k,45) + rxt(k,344)*y(k,52) &
                      + rxt(k,345)*y(k,98) + rxt(k,346)*y(k,99) + rxt(k,347)*y(k,100) &
                      + (rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,64) + rxt(k,367) &
                      *y(k,102))
         mat(k,1540) = -rxt(k,153)*y(k,247)
         mat(k,2693) = -rxt(k,154)*y(k,247)
         mat(k,1937) = -(rxt(k,157) + rxt(k,158)) * y(k,247)
         mat(k,214) = -(rxt(k,196) + rxt(k,197)) * y(k,247)
         mat(k,115) = -rxt(k,231)*y(k,247)
         mat(k,170) = -rxt(k,232)*y(k,247)
         mat(k,135) = -rxt(k,233)*y(k,247)
         mat(k,181) = -rxt(k,234)*y(k,247)
         mat(k,139) = -rxt(k,235)*y(k,247)
         mat(k,186) = -rxt(k,236)*y(k,247)
         mat(k,143) = -rxt(k,237)*y(k,247)
         mat(k,1743) = -(rxt(k,238) + rxt(k,239)) * y(k,247)
         mat(k,175) = -rxt(k,262)*y(k,247)
         mat(k,503) = -rxt(k,263)*y(k,247)
         mat(k,126) = -rxt(k,264)*y(k,247)
         mat(k,1512) = -(rxt(k,265) + rxt(k,266)) * y(k,247)
         mat(k,271) = -rxt(k,310)*y(k,247)
         mat(k,263) = -rxt(k,311)*y(k,247)
         mat(k,555) = -rxt(k,343)*y(k,247)
         mat(k,669) = -rxt(k,344)*y(k,247)
         mat(k,258) = -rxt(k,345)*y(k,247)
         mat(k,267) = -rxt(k,346)*y(k,247)
         mat(k,318) = -rxt(k,347)*y(k,247)
         mat(k,1573) = -(rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,247)
         mat(k,206) = -rxt(k,367)*y(k,247)
         mat(k,2166) = -(rxt(k,171)*y(k,93) + rxt(k,172)*y(k,95) + rxt(k,173)*y(k,233) &
                      + rxt(k,174)*y(k,162) + rxt(k,175)*y(k,163) + (4._r8*rxt(k,176) &
                      + 4._r8*rxt(k,177)) * y(k,248) + rxt(k,179)*y(k,107) + rxt(k,191) &
                      *y(k,155) + rxt(k,192)*y(k,141) + rxt(k,200)*y(k,154) + rxt(k,201) &
                      *y(k,106) + rxt(k,211)*y(k,73) + rxt(k,222)*y(k,75) + (rxt(k,224) &
                      + rxt(k,225)) * y(k,74) + rxt(k,227)*y(k,101) + rxt(k,230) &
                      *y(k,109) + rxt(k,242)*y(k,18) + rxt(k,258)*y(k,21) + rxt(k,260) &
                      *y(k,97) + rxt(k,268)*y(k,110) + rxt(k,271)*y(k,116) + rxt(k,294) &
                      *y(k,125) + rxt(k,295)*y(k,105) + rxt(k,313)*y(k,26) + rxt(k,315) &
                      *y(k,29) + rxt(k,317)*y(k,45) + rxt(k,318)*y(k,46) + rxt(k,320) &
                      *y(k,47) + rxt(k,322)*y(k,52) + rxt(k,323)*y(k,53) + rxt(k,325) &
                      *y(k,55) + rxt(k,327)*y(k,61) + rxt(k,328)*y(k,65) + rxt(k,330) &
                      *y(k,66) + rxt(k,331)*y(k,67) + rxt(k,339)*y(k,69) + rxt(k,340) &
                      *y(k,98) + rxt(k,341)*y(k,99) + rxt(k,342)*y(k,100) + rxt(k,351) &
                      *y(k,51) + rxt(k,356)*y(k,62) + rxt(k,357)*y(k,63) + rxt(k,358) &
                      *y(k,64) + rxt(k,359)*y(k,102) + rxt(k,360)*y(k,103) + rxt(k,368) &
                      *y(k,77) + rxt(k,370)*y(k,27) + rxt(k,377)*y(k,30) + rxt(k,378) &
                      *y(k,31) + rxt(k,380)*y(k,32) + rxt(k,382)*y(k,54) + rxt(k,383) &
                      *y(k,56) + rxt(k,388)*y(k,59) + rxt(k,389)*y(k,60) + rxt(k,394) &
                      *y(k,90) + rxt(k,395)*y(k,91) + rxt(k,396)*y(k,169) + rxt(k,397) &
                      *y(k,28) + rxt(k,405)*y(k,34) + rxt(k,406)*y(k,35) + rxt(k,408) &
                      *y(k,58) + rxt(k,410)*y(k,113) + rxt(k,411)*y(k,156) + rxt(k,414) &
                      *y(k,176) + rxt(k,418)*y(k,177) + rxt(k,419)*y(k,33) + rxt(k,420) &
                      *y(k,57) + rxt(k,422)*y(k,16) + rxt(k,425)*y(k,111) + rxt(k,433) &
                      *y(k,134) + rxt(k,434)*y(k,135) + rxt(k,443)*y(k,136) + rxt(k,444) &
                      *y(k,137) + rxt(k,445)*y(k,138) + rxt(k,447)*y(k,140) + rxt(k,450) &
                      *y(k,1) + rxt(k,454)*y(k,2) + rxt(k,455)*y(k,15) + rxt(k,456) &
                      *y(k,112) + rxt(k,457)*y(k,114) + rxt(k,458)*y(k,122) + rxt(k,470) &
                      *y(k,128) + rxt(k,471)*y(k,129) + rxt(k,478)*y(k,130) + rxt(k,480) &
                      *y(k,127) + rxt(k,481)*y(k,131) + rxt(k,482)*y(k,144) + rxt(k,483) &
                      *y(k,145) + rxt(k,489)*y(k,211) + rxt(k,492)*y(k,7) + rxt(k,495) &
                      *y(k,8) + rxt(k,496)*y(k,24) + rxt(k,498)*y(k,25) + rxt(k,502) &
                      *y(k,36) + rxt(k,503)*y(k,82) + rxt(k,515)*y(k,172) + rxt(k,518) &
                      *y(k,173) + rxt(k,522)*y(k,209) + rxt(k,523)*y(k,210) + rxt(k,525) &
                      *y(k,212) + rxt(k,528)*y(k,213) + rxt(k,531)*y(k,214) + rxt(k,532) &
                      *y(k,215) + rxt(k,535)*y(k,6) + rxt(k,538)*y(k,139) + rxt(k,543) &
                      *y(k,157) + rxt(k,547)*y(k,204) + rxt(k,548)*y(k,205) + rxt(k,552) &
                      *y(k,206) + rxt(k,554)*y(k,207) + rxt(k,555)*y(k,208) + (rxt(k,561) &
                      + rxt(k,575)) * y(k,83) + rxt(k,563)*y(k,166) + rxt(k,565) &
                      *y(k,181) + rxt(k,569)*y(k,178) + rxt(k,574)*y(k,180) + rxt(k,594) &
                      *y(k,149))
         mat(k,1541) = -rxt(k,171)*y(k,248)
         mat(k,717) = -rxt(k,172)*y(k,248)
         mat(k,2337) = -rxt(k,173)*y(k,248)
         mat(k,2429) = -rxt(k,174)*y(k,248)
         mat(k,1938) = -rxt(k,175)*y(k,248)
         mat(k,472) = -rxt(k,179)*y(k,248)
         mat(k,1812) = -rxt(k,191)*y(k,248)
         mat(k,578) = -rxt(k,192)*y(k,248)
         mat(k,2564) = -rxt(k,200)*y(k,248)
         mat(k,2191) = -rxt(k,201)*y(k,248)
         mat(k,623) = -rxt(k,211)*y(k,248)
         mat(k,1014) = -rxt(k,222)*y(k,248)
         mat(k,2366) = -(rxt(k,224) + rxt(k,225)) * y(k,248)
         mat(k,1744) = -rxt(k,227)*y(k,248)
         mat(k,1721) = -rxt(k,230)*y(k,248)
         mat(k,543) = -rxt(k,242)*y(k,248)
         mat(k,2456) = -rxt(k,258)*y(k,248)
         mat(k,1513) = -rxt(k,260)*y(k,248)
         mat(k,1592) = -rxt(k,268)*y(k,248)
         mat(k,1525) = -rxt(k,271)*y(k,248)
         mat(k,2222) = -rxt(k,294)*y(k,248)
         mat(k,1265) = -rxt(k,295)*y(k,248)
         mat(k,219) = -rxt(k,313)*y(k,248)
         mat(k,289) = -rxt(k,315)*y(k,248)
         mat(k,556) = -rxt(k,317)*y(k,248)
         mat(k,146) = -rxt(k,318)*y(k,248)
         mat(k,346) = -rxt(k,320)*y(k,248)
         mat(k,670) = -rxt(k,322)*y(k,248)
         mat(k,150) = -rxt(k,323)*y(k,248)
         mat(k,443) = -rxt(k,325)*y(k,248)
         mat(k,434) = -rxt(k,327)*y(k,248)
         mat(k,118) = -rxt(k,328)*y(k,248)
         mat(k,504) = -rxt(k,330)*y(k,248)
         mat(k,122) = -rxt(k,331)*y(k,248)
         mat(k,401) = -rxt(k,339)*y(k,248)
         mat(k,259) = -rxt(k,340)*y(k,248)
         mat(k,268) = -rxt(k,341)*y(k,248)
         mat(k,319) = -rxt(k,342)*y(k,248)
         mat(k,1695) = -rxt(k,351)*y(k,248)
         mat(k,870) = -rxt(k,356)*y(k,248)
         mat(k,479) = -rxt(k,357)*y(k,248)
         mat(k,1574) = -rxt(k,358)*y(k,248)
         mat(k,207) = -rxt(k,359)*y(k,248)
         mat(k,1022) = -rxt(k,360)*y(k,248)
         mat(k,1201) = -rxt(k,368)*y(k,248)
         mat(k,325) = -rxt(k,370)*y(k,248)
         mat(k,298) = -rxt(k,377)*y(k,248)
         mat(k,390) = -rxt(k,378)*y(k,248)
         mat(k,331) = -rxt(k,380)*y(k,248)
         mat(k,1195) = -rxt(k,382)*y(k,248)
         mat(k,108) = -rxt(k,383)*y(k,248)
         mat(k,778) = -rxt(k,388)*y(k,248)
         mat(k,679) = -rxt(k,389)*y(k,248)
         mat(k,1209) = -rxt(k,394)*y(k,248)
         mat(k,1049) = -rxt(k,395)*y(k,248)
         mat(k,631) = -rxt(k,396)*y(k,248)
         mat(k,647) = -rxt(k,397)*y(k,248)
         mat(k,491) = -rxt(k,405)*y(k,248)
         mat(k,337) = -rxt(k,406)*y(k,248)
         mat(k,1342) = -rxt(k,408)*y(k,248)
         mat(k,1254) = -rxt(k,410)*y(k,248)
         mat(k,902) = -rxt(k,411)*y(k,248)
         mat(k,638) = -rxt(k,414)*y(k,248)
         mat(k,486) = -rxt(k,418)*y(k,248)
         mat(k,1140) = -rxt(k,419)*y(k,248)
         mat(k,1110) = -rxt(k,420)*y(k,248)
         mat(k,417) = -rxt(k,422)*y(k,248)
         mat(k,1244) = -rxt(k,425)*y(k,248)
         mat(k,1334) = -rxt(k,433)*y(k,248)
         mat(k,357) = -rxt(k,434)*y(k,248)
         mat(k,602) = -rxt(k,443)*y(k,248)
         mat(k,363) = -rxt(k,444)*y(k,248)
         mat(k,687) = -rxt(k,445)*y(k,248)
         mat(k,1445) = -rxt(k,447)*y(k,248)
         mat(k,764) = -rxt(k,450)*y(k,248)
         mat(k,744) = -rxt(k,454)*y(k,248)
         mat(k,279) = -rxt(k,455)*y(k,248)
         mat(k,275) = -rxt(k,456)*y(k,248)
         mat(k,353) = -rxt(k,457)*y(k,248)
         mat(k,163) = -rxt(k,458)*y(k,248)
         mat(k,695) = -rxt(k,470)*y(k,248)
         mat(k,662) = -rxt(k,471)*y(k,248)
         mat(k,449) = -rxt(k,478)*y(k,248)
         mat(k,933) = -rxt(k,480)*y(k,248)
         mat(k,810) = -rxt(k,481)*y(k,248)
         mat(k,461) = -rxt(k,482)*y(k,248)
         mat(k,1156) = -rxt(k,483)*y(k,248)
         mat(k,239) = -rxt(k,489)*y(k,248)
         mat(k,192) = -rxt(k,492)*y(k,248)
         mat(k,498) = -rxt(k,495)*y(k,248)
         mat(k,282) = -rxt(k,496)*y(k,248)
         mat(k,385) = -rxt(k,498)*y(k,248)
         mat(k,303) = -rxt(k,502)*y(k,248)
         mat(k,231) = -rxt(k,503)*y(k,248)
         mat(k,201) = -rxt(k,515)*y(k,248)
         mat(k,379) = -rxt(k,518)*y(k,248)
         mat(k,708) = -rxt(k,522)*y(k,248)
         mat(k,226) = -rxt(k,523)*y(k,248)
         mat(k,248) = -rxt(k,525)*y(k,248)
         mat(k,793) = -rxt(k,528)*y(k,248)
         mat(k,253) = -rxt(k,531)*y(k,248)
         mat(k,510) = -rxt(k,532)*y(k,248)
         mat(k,1079) = -rxt(k,535)*y(k,248)
         mat(k,990) = -rxt(k,538)*y(k,248)
         mat(k,467) = -rxt(k,543)*y(k,248)
         mat(k,754) = -rxt(k,547)*y(k,248)
         mat(k,723) = -rxt(k,548)*y(k,248)
         mat(k,564) = -rxt(k,552)*y(k,248)
         mat(k,1104) = -rxt(k,554)*y(k,248)
         mat(k,1186) = -rxt(k,555)*y(k,248)
         mat(k,1031) = -(rxt(k,561) + rxt(k,575)) * y(k,248)
         mat(k,424) = -rxt(k,563)*y(k,248)
         mat(k,1041) = -rxt(k,565)*y(k,248)
         mat(k,606) = -rxt(k,569)*y(k,248)
         mat(k,1552) = -rxt(k,574)*y(k,248)
         mat(k,111) = -rxt(k,594)*y(k,248)
         mat(k,1079) = mat(k,1079) + .630_r8*rxt(k,534)*y(k,163)
         mat(k,325) = mat(k,325) + .650_r8*rxt(k,370)*y(k,248)
         mat(k,647) = mat(k,647) + .130_r8*rxt(k,372)*y(k,163)
         mat(k,390) = mat(k,390) + .500_r8*rxt(k,378)*y(k,248)
         mat(k,1140) = mat(k,1140) + .360_r8*rxt(k,401)*y(k,163)
         mat(k,1695) = mat(k,1695) + rxt(k,350)*y(k,162)
         mat(k,479) = mat(k,479) + .300_r8*rxt(k,357)*y(k,248)
         mat(k,1574) = mat(k,1574) + rxt(k,364)*y(k,247)
         mat(k,1872) = rxt(k,209)*y(k,233)
         mat(k,952) = rxt(k,308)*y(k,260)
         mat(k,2507) = rxt(k,170)*y(k,163) + 2.000_r8*rxt(k,165)*y(k,233)
         mat(k,1541) = mat(k,1541) + rxt(k,162)*y(k,162) + rxt(k,153)*y(k,247)
         mat(k,717) = mat(k,717) + rxt(k,163)*y(k,162)
         mat(k,1513) = mat(k,1513) + rxt(k,259)*y(k,162) + rxt(k,265)*y(k,247)
         mat(k,1744) = mat(k,1744) + rxt(k,226)*y(k,162) + rxt(k,238)*y(k,247)
         mat(k,207) = mat(k,207) + rxt(k,367)*y(k,247)
         mat(k,1669) = rxt(k,261)*y(k,162)
         mat(k,1721) = mat(k,1721) + rxt(k,229)*y(k,162)
         mat(k,933) = mat(k,933) + .320_r8*rxt(k,479)*y(k,163)
         mat(k,810) = mat(k,810) + .600_r8*rxt(k,481)*y(k,248)
         mat(k,1334) = mat(k,1334) + .240_r8*rxt(k,432)*y(k,163)
         mat(k,357) = mat(k,357) + .100_r8*rxt(k,434)*y(k,248)
         mat(k,990) = mat(k,990) + .630_r8*rxt(k,537)*y(k,163)
         mat(k,1445) = mat(k,1445) + .360_r8*rxt(k,446)*y(k,163)
         mat(k,2665) = rxt(k,193)*y(k,233)
         mat(k,1812) = mat(k,1812) + rxt(k,188)*y(k,233)
         mat(k,2429) = mat(k,2429) + rxt(k,350)*y(k,51) + rxt(k,162)*y(k,93) &
                      + rxt(k,163)*y(k,95) + rxt(k,259)*y(k,97) + rxt(k,226)*y(k,101) &
                      + rxt(k,261)*y(k,108) + rxt(k,229)*y(k,109) + rxt(k,168) &
                      *y(k,233)
         mat(k,1938) = mat(k,1938) + .630_r8*rxt(k,534)*y(k,6) + .130_r8*rxt(k,372) &
                      *y(k,28) + .360_r8*rxt(k,401)*y(k,33) + rxt(k,170)*y(k,92) &
                      + .320_r8*rxt(k,479)*y(k,127) + .240_r8*rxt(k,432)*y(k,134) &
                      + .630_r8*rxt(k,537)*y(k,139) + .360_r8*rxt(k,446)*y(k,140) &
                      + rxt(k,169)*y(k,233)
         mat(k,638) = mat(k,638) + .500_r8*rxt(k,414)*y(k,248)
         mat(k,239) = mat(k,239) + .500_r8*rxt(k,489)*y(k,248)
         mat(k,612) = .400_r8*rxt(k,490)*y(k,233)
         mat(k,1497) = .450_r8*rxt(k,386)*y(k,233)
         mat(k,854) = .400_r8*rxt(k,504)*y(k,233)
         mat(k,2337) = mat(k,2337) + rxt(k,209)*y(k,70) + 2.000_r8*rxt(k,165)*y(k,92) &
                      + rxt(k,193)*y(k,153) + rxt(k,188)*y(k,155) + rxt(k,168) &
                      *y(k,162) + rxt(k,169)*y(k,163) + .400_r8*rxt(k,490)*y(k,218) &
                      + .450_r8*rxt(k,386)*y(k,227) + .400_r8*rxt(k,504)*y(k,229) &
                      + .450_r8*rxt(k,437)*y(k,242) + .400_r8*rxt(k,510)*y(k,243) &
                      + .200_r8*rxt(k,441)*y(k,244) + .150_r8*rxt(k,416)*y(k,251)
         mat(k,1465) = .450_r8*rxt(k,437)*y(k,233)
         mat(k,961) = .400_r8*rxt(k,510)*y(k,233)
         mat(k,771) = .200_r8*rxt(k,441)*y(k,233)
         mat(k,1984) = rxt(k,364)*y(k,64) + rxt(k,153)*y(k,93) + rxt(k,265)*y(k,97) &
                      + rxt(k,238)*y(k,101) + rxt(k,367)*y(k,102) &
                      + 2.000_r8*rxt(k,154)*y(k,260)
         mat(k,2166) = mat(k,2166) + .650_r8*rxt(k,370)*y(k,27) + .500_r8*rxt(k,378) &
                      *y(k,31) + .300_r8*rxt(k,357)*y(k,63) + .600_r8*rxt(k,481) &
                      *y(k,131) + .100_r8*rxt(k,434)*y(k,135) + .500_r8*rxt(k,414) &
                      *y(k,176) + .500_r8*rxt(k,489)*y(k,211)
         mat(k,1279) = .150_r8*rxt(k,416)*y(k,233)
         mat(k,2694) = rxt(k,308)*y(k,89) + 2.000_r8*rxt(k,154)*y(k,247)
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
         mat(k,527) = -(rxt(k,513)*y(k,233) + rxt(k,514)*y(k,153))
         mat(k,2267) = -rxt(k,513)*y(k,249)
         mat(k,2598) = -rxt(k,514)*y(k,249)
         mat(k,229) = .200_r8*rxt(k,503)*y(k,248)
         mat(k,199) = .140_r8*rxt(k,515)*y(k,248)
         mat(k,377) = rxt(k,518)*y(k,248)
         mat(k,2069) = .200_r8*rxt(k,503)*y(k,82) + .140_r8*rxt(k,515)*y(k,172) &
                      + rxt(k,518)*y(k,173)
         mat(k,859) = -(rxt(k,412)*y(k,233) + rxt(k,413)*y(k,153))
         mat(k,2291) = -rxt(k,412)*y(k,250)
         mat(k,2620) = -rxt(k,413)*y(k,250)
         mat(k,1125) = rxt(k,419)*y(k,248)
         mat(k,634) = .500_r8*rxt(k,414)*y(k,248)
         mat(k,2108) = rxt(k,419)*y(k,33) + .500_r8*rxt(k,414)*y(k,176)
         mat(k,1274) = -(rxt(k,415)*y(k,228) + rxt(k,416)*y(k,233) + rxt(k,417) &
                      *y(k,153))
         mat(k,1627) = -rxt(k,415)*y(k,251)
         mat(k,2312) = -rxt(k,416)*y(k,251)
         mat(k,2642) = -rxt(k,417)*y(k,251)
         mat(k,1072) = .060_r8*rxt(k,534)*y(k,163)
         mat(k,1107) = rxt(k,420)*y(k,248)
         mat(k,984) = .060_r8*rxt(k,537)*y(k,163)
         mat(k,1918) = .060_r8*rxt(k,534)*y(k,6) + .060_r8*rxt(k,537)*y(k,139)
         mat(k,483) = rxt(k,418)*y(k,248)
         mat(k,1182) = .150_r8*rxt(k,555)*y(k,248)
         mat(k,2140) = rxt(k,420)*y(k,57) + rxt(k,418)*y(k,177) + .150_r8*rxt(k,555) &
                      *y(k,208)
         mat(k,1220) = -(rxt(k,544)*y(k,228) + rxt(k,545)*y(k,233) + rxt(k,546) &
                      *y(k,153))
         mat(k,1625) = -rxt(k,544)*y(k,252)
         mat(k,2309) = -rxt(k,545)*y(k,252)
         mat(k,2639) = -rxt(k,546)*y(k,252)
         mat(k,1784) = .500_r8*rxt(k,553)*y(k,207)
         mat(k,752) = rxt(k,547)*y(k,248)
         mat(k,1102) = .500_r8*rxt(k,553)*y(k,155) + rxt(k,554)*y(k,248)
         mat(k,2136) = rxt(k,547)*y(k,204) + rxt(k,554)*y(k,207)
         mat(k,1088) = -(rxt(k,549)*y(k,228) + rxt(k,550)*y(k,233) + rxt(k,551) &
                      *y(k,153))
         mat(k,1615) = -rxt(k,549)*y(k,253)
         mat(k,2300) = -rxt(k,550)*y(k,253)
         mat(k,2629) = -rxt(k,551)*y(k,253)
         mat(k,1065) = rxt(k,535)*y(k,248)
         mat(k,978) = rxt(k,538)*y(k,248)
         mat(k,560) = rxt(k,552)*y(k,248)
         mat(k,2125) = rxt(k,535)*y(k,6) + rxt(k,538)*y(k,139) + rxt(k,552)*y(k,206)
         mat(k,823) = -(rxt(k,520)*y(k,233) + rxt(k,521)*y(k,153))
         mat(k,2288) = -rxt(k,520)*y(k,254)
         mat(k,2617) = -rxt(k,521)*y(k,254)
         mat(k,704) = rxt(k,522)*y(k,248)
         mat(k,225) = .650_r8*rxt(k,523)*y(k,248)
         mat(k,2105) = rxt(k,522)*y(k,209) + .650_r8*rxt(k,523)*y(k,210)
         mat(k,95) = -(rxt(k,649)*y(k,233) + rxt(k,650)*y(k,153))
         mat(k,2245) = -rxt(k,649)*y(k,255)
         mat(k,2587) = -rxt(k,650)*y(k,255)
         mat(k,220) = rxt(k,648)*y(k,248)
         mat(k,2009) = rxt(k,648)*y(k,210)
         mat(k,1291) = -(rxt(k,484)*y(k,227) + rxt(k,485)*y(k,228) + rxt(k,486) &
                      *y(k,233) + rxt(k,487)*y(k,153) + rxt(k,488)*y(k,155))
         mat(k,1482) = -rxt(k,484)*y(k,256)
         mat(k,1628) = -rxt(k,485)*y(k,256)
         mat(k,2313) = -rxt(k,486)*y(k,256)
         mat(k,2643) = -rxt(k,487)*y(k,256)
         mat(k,1789) = -rxt(k,488)*y(k,256)
         mat(k,274) = rxt(k,456)*y(k,248)
         mat(k,352) = rxt(k,457)*y(k,248)
         mat(k,162) = rxt(k,458)*y(k,248)
         mat(k,806) = .400_r8*rxt(k,481)*y(k,248)
         mat(k,238) = .500_r8*rxt(k,489)*y(k,248)
         mat(k,2141) = rxt(k,456)*y(k,112) + rxt(k,457)*y(k,114) + rxt(k,458)*y(k,122) &
                      + .400_r8*rxt(k,481)*y(k,131) + .500_r8*rxt(k,489)*y(k,211)
         mat(k,839) = -(rxt(k,526)*y(k,233) + rxt(k,527)*y(k,153))
         mat(k,2289) = -rxt(k,526)*y(k,257)
         mat(k,2618) = -rxt(k,527)*y(k,257)
         mat(k,245) = .560_r8*rxt(k,525)*y(k,248)
         mat(k,786) = rxt(k,528)*y(k,248)
         mat(k,2106) = .560_r8*rxt(k,525)*y(k,212) + rxt(k,528)*y(k,213)
         mat(k,101) = -(rxt(k,653)*y(k,233) + rxt(k,654)*y(k,153))
         mat(k,2246) = -rxt(k,653)*y(k,258)
         mat(k,2588) = -rxt(k,654)*y(k,258)
         mat(k,240) = rxt(k,652)*y(k,248)
         mat(k,2010) = rxt(k,652)*y(k,212)
         mat(k,584) = -(rxt(k,529)*y(k,233) + rxt(k,530)*y(k,153))
         mat(k,2272) = -rxt(k,529)*y(k,259)
         mat(k,2603) = -rxt(k,530)*y(k,259)
         mat(k,252) = .300_r8*rxt(k,531)*y(k,248)
         mat(k,507) = rxt(k,532)*y(k,248)
         mat(k,2077) = .300_r8*rxt(k,531)*y(k,214) + rxt(k,532)*y(k,215)
         mat(k,2706) = -(rxt(k,154)*y(k,247) + rxt(k,308)*y(k,89) + rxt(k,576) &
                      *y(k,182))
         mat(k,1996) = -rxt(k,154)*y(k,260)
         mat(k,955) = -rxt(k,308)*y(k,260)
         mat(k,295) = -rxt(k,576)*y(k,260)
         mat(k,290) = rxt(k,315)*y(k,248)
         mat(k,332) = rxt(k,380)*y(k,248)
         mat(k,493) = rxt(k,405)*y(k,248)
         mat(k,338) = rxt(k,406)*y(k,248)
         mat(k,558) = rxt(k,317)*y(k,248)
         mat(k,347) = rxt(k,320)*y(k,248)
         mat(k,1707) = rxt(k,351)*y(k,248)
         mat(k,673) = rxt(k,322)*y(k,248)
         mat(k,151) = rxt(k,323)*y(k,248)
         mat(k,1198) = rxt(k,382)*y(k,248)
         mat(k,445) = rxt(k,325)*y(k,248)
         mat(k,1111) = rxt(k,420)*y(k,248)
         mat(k,1346) = rxt(k,408)*y(k,248)
         mat(k,779) = rxt(k,388)*y(k,248)
         mat(k,680) = rxt(k,389)*y(k,248)
         mat(k,437) = rxt(k,327)*y(k,248)
         mat(k,481) = rxt(k,357)*y(k,248)
         mat(k,1579) = rxt(k,358)*y(k,248)
         mat(k,1123) = rxt(k,334)*y(k,233)
         mat(k,402) = rxt(k,339)*y(k,248)
         mat(k,2519) = rxt(k,166)*y(k,233)
         mat(k,1545) = rxt(k,171)*y(k,248)
         mat(k,720) = rxt(k,172)*y(k,248)
         mat(k,1519) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,108) + ( &
                      + rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,109) + ( &
                      + rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,110) &
                      + rxt(k,260)*y(k,248)
         mat(k,320) = rxt(k,342)*y(k,248)
         mat(k,1754) = (rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,108) + ( &
                      + rxt(k,585)+rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,109) + ( &
                      + rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,110) &
                      + rxt(k,227)*y(k,248)
         mat(k,1025) = rxt(k,360)*y(k,248)
         mat(k,1270) = (rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,108) + ( &
                      + rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,109) + ( &
                      + rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,110) &
                      + rxt(k,295)*y(k,248)
         mat(k,2203) = rxt(k,201)*y(k,248)
         mat(k,475) = rxt(k,179)*y(k,248)
         mat(k,1678) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,97) + ( &
                      + rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,101) + ( &
                      + rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,105)
         mat(k,1731) = (rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,97) + ( &
                      + rxt(k,585)+rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,101) + ( &
                      + rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,105) &
                      + rxt(k,230)*y(k,248)
         mat(k,1601) = (rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,97) + ( &
                      + rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,101) + ( &
                      + rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,105) &
                      + rxt(k,268)*y(k,248)
         mat(k,1337) = .500_r8*rxt(k,433)*y(k,248)
         mat(k,112) = rxt(k,594)*y(k,248)
         mat(k,640) = rxt(k,414)*y(k,248)
         mat(k,487) = rxt(k,418)*y(k,248)
         mat(k,2349) = rxt(k,334)*y(k,68) + rxt(k,166)*y(k,92) + rxt(k,173)*y(k,248)
         mat(k,2178) = rxt(k,315)*y(k,29) + rxt(k,380)*y(k,32) + rxt(k,405)*y(k,34) &
                      + rxt(k,406)*y(k,35) + rxt(k,317)*y(k,45) + rxt(k,320)*y(k,47) &
                      + rxt(k,351)*y(k,51) + rxt(k,322)*y(k,52) + rxt(k,323)*y(k,53) &
                      + rxt(k,382)*y(k,54) + rxt(k,325)*y(k,55) + rxt(k,420)*y(k,57) &
                      + rxt(k,408)*y(k,58) + rxt(k,388)*y(k,59) + rxt(k,389)*y(k,60) &
                      + rxt(k,327)*y(k,61) + rxt(k,357)*y(k,63) + rxt(k,358)*y(k,64) &
                      + rxt(k,339)*y(k,69) + rxt(k,171)*y(k,93) + rxt(k,172)*y(k,95) &
                      + rxt(k,260)*y(k,97) + rxt(k,342)*y(k,100) + rxt(k,227)*y(k,101) &
                      + rxt(k,360)*y(k,103) + rxt(k,295)*y(k,105) + rxt(k,201) &
                      *y(k,106) + rxt(k,179)*y(k,107) + rxt(k,230)*y(k,109) &
                      + rxt(k,268)*y(k,110) + .500_r8*rxt(k,433)*y(k,134) + rxt(k,594) &
                      *y(k,149) + rxt(k,414)*y(k,176) + rxt(k,418)*y(k,177) &
                      + rxt(k,173)*y(k,233) + 2.000_r8*rxt(k,176)*y(k,248)
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
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = lmat(k, 103)
         mat(k, 104) = lmat(k, 104)
         mat(k, 105) = lmat(k, 105)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 128) = lmat(k, 128)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = lmat(k, 157)
         mat(k, 158) = lmat(k, 158)
         mat(k, 159) = lmat(k, 159)
         mat(k, 160) = lmat(k, 160)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 164) = lmat(k, 164)
         mat(k, 165) = lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 194) = lmat(k, 194)
         mat(k, 195) = lmat(k, 195)
         mat(k, 196) = lmat(k, 196)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = mat(k, 198) + lmat(k, 198)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 211) = lmat(k, 211)
         mat(k, 212) = mat(k, 212) + lmat(k, 212)
         mat(k, 213) = mat(k, 213) + lmat(k, 213)
         mat(k, 214) = mat(k, 214) + lmat(k, 214)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 218) = mat(k, 218) + lmat(k, 218)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = lmat(k, 285)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 293) = lmat(k, 293)
         mat(k, 294) = lmat(k, 294)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 301) = lmat(k, 301)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 304) = lmat(k, 304)
         mat(k, 305) = lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 307) = lmat(k, 307)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 310) = lmat(k, 310)
         mat(k, 312) = lmat(k, 312)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 339) = lmat(k, 339)
         mat(k, 340) = lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = lmat(k, 349)
         mat(k, 350) = lmat(k, 350)
         mat(k, 351) = mat(k, 351) + lmat(k, 351)
         mat(k, 354) = mat(k, 354) + lmat(k, 354)
         mat(k, 359) = mat(k, 359) + lmat(k, 359)
         mat(k, 361) = lmat(k, 361)
         mat(k, 362) = lmat(k, 362)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 364) = lmat(k, 364)
         mat(k, 365) = lmat(k, 365)
         mat(k, 366) = lmat(k, 366)
         mat(k, 367) = lmat(k, 367)
         mat(k, 368) = lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = lmat(k, 370)
         mat(k, 371) = lmat(k, 371)
         mat(k, 372) = lmat(k, 372)
         mat(k, 373) = lmat(k, 373)
         mat(k, 374) = lmat(k, 374)
         mat(k, 375) = lmat(k, 375)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 378) = lmat(k, 378)
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = mat(k, 382) + lmat(k, 382)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 390) = mat(k, 390) + lmat(k, 390)
         mat(k, 391) = lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 393) = lmat(k, 393)
         mat(k, 394) = lmat(k, 394)
         mat(k, 396) = mat(k, 396) + lmat(k, 396)
         mat(k, 397) = lmat(k, 397)
         mat(k, 400) = mat(k, 400) + lmat(k, 400)
         mat(k, 403) = lmat(k, 403)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 406) = lmat(k, 406)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 408) = lmat(k, 408)
         mat(k, 409) = lmat(k, 409)
         mat(k, 410) = lmat(k, 410)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 420) = lmat(k, 420)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 427) = lmat(k, 427)
         mat(k, 428) = lmat(k, 428)
         mat(k, 429) = lmat(k, 429)
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 431) = lmat(k, 431)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 440) = lmat(k, 440)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 447) = lmat(k, 447)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 451) = lmat(k, 451)
         mat(k, 452) = lmat(k, 452)
         mat(k, 453) = lmat(k, 453)
         mat(k, 454) = lmat(k, 454)
         mat(k, 455) = lmat(k, 455)
         mat(k, 456) = lmat(k, 456)
         mat(k, 457) = lmat(k, 457)
         mat(k, 458) = mat(k, 458) + lmat(k, 458)
         mat(k, 462) = lmat(k, 462)
         mat(k, 464) = mat(k, 464) + lmat(k, 464)
         mat(k, 465) = lmat(k, 465)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 468) = lmat(k, 468)
         mat(k, 469) = lmat(k, 469)
         mat(k, 470) = mat(k, 470) + lmat(k, 470)
         mat(k, 471) = lmat(k, 471)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 473) = lmat(k, 473)
         mat(k, 474) = mat(k, 474) + lmat(k, 474)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 478) = mat(k, 478) + lmat(k, 478)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 480) = lmat(k, 480)
         mat(k, 482) = mat(k, 482) + lmat(k, 482)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = mat(k, 486) + lmat(k, 486)
         mat(k, 488) = mat(k, 488) + lmat(k, 488)
         mat(k, 490) = lmat(k, 490)
         mat(k, 491) = mat(k, 491) + lmat(k, 491)
         mat(k, 492) = lmat(k, 492)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 495) = lmat(k, 495)
         mat(k, 497) = lmat(k, 497)
         mat(k, 498) = mat(k, 498) + lmat(k, 498)
         mat(k, 499) = lmat(k, 499)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 506) = mat(k, 506) + lmat(k, 506)
         mat(k, 508) = lmat(k, 508)
         mat(k, 509) = lmat(k, 509)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 511) = lmat(k, 511)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 520) = mat(k, 520) + lmat(k, 520)
         mat(k, 522) = lmat(k, 522)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 537) = lmat(k, 537)
         mat(k, 538) = lmat(k, 538)
         mat(k, 539) = lmat(k, 539)
         mat(k, 540) = lmat(k, 540)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 545) = mat(k, 545) + lmat(k, 545)
         mat(k, 549) = lmat(k, 549)
         mat(k, 550) = lmat(k, 550)
         mat(k, 551) = lmat(k, 551)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 557) = mat(k, 557) + lmat(k, 557)
         mat(k, 559) = mat(k, 559) + lmat(k, 559)
         mat(k, 561) = lmat(k, 561)
         mat(k, 562) = lmat(k, 562)
         mat(k, 563) = lmat(k, 563)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 565) = lmat(k, 565)
         mat(k, 568) = mat(k, 568) + lmat(k, 568)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 579) = mat(k, 579) + lmat(k, 579)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 593) = mat(k, 593) + lmat(k, 593)
         mat(k, 595) = lmat(k, 595)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 600) = lmat(k, 600)
         mat(k, 601) = lmat(k, 601)
         mat(k, 603) = mat(k, 603) + lmat(k, 603)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 607) = lmat(k, 607)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 616) = lmat(k, 616)
         mat(k, 617) = lmat(k, 617)
         mat(k, 618) = lmat(k, 618)
         mat(k, 619) = lmat(k, 619)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 622) = lmat(k, 622)
         mat(k, 624) = mat(k, 624) + lmat(k, 624)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 626) = lmat(k, 626)
         mat(k, 627) = lmat(k, 627)
         mat(k, 628) = lmat(k, 628)
         mat(k, 630) = mat(k, 630) + lmat(k, 630)
         mat(k, 632) = lmat(k, 632)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 635) = lmat(k, 635)
         mat(k, 637) = lmat(k, 637)
         mat(k, 638) = mat(k, 638) + lmat(k, 638)
         mat(k, 639) = lmat(k, 639)
         mat(k, 641) = mat(k, 641) + lmat(k, 641)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 657) = mat(k, 657) + lmat(k, 657)
         mat(k, 663) = lmat(k, 663)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 666) = lmat(k, 666)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 674) = mat(k, 674) + lmat(k, 674)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 677) = lmat(k, 677)
         mat(k, 679) = mat(k, 679) + lmat(k, 679)
         mat(k, 681) = mat(k, 681) + lmat(k, 681)
         mat(k, 684) = lmat(k, 684)
         mat(k, 689) = lmat(k, 689)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 696) = lmat(k, 696)
         mat(k, 699) = lmat(k, 699)
         mat(k, 700) = lmat(k, 700)
         mat(k, 701) = mat(k, 701) + lmat(k, 701)
         mat(k, 702) = lmat(k, 702)
         mat(k, 706) = lmat(k, 706)
         mat(k, 707) = lmat(k, 707)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 709) = lmat(k, 709)
         mat(k, 710) = lmat(k, 710)
         mat(k, 711) = lmat(k, 711)
         mat(k, 712) = lmat(k, 712)
         mat(k, 713) = lmat(k, 713)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 717) = mat(k, 717) + lmat(k, 717)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 722) = mat(k, 722) + lmat(k, 722)
         mat(k, 724) = lmat(k, 724)
         mat(k, 725) = lmat(k, 725)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 729) = mat(k, 729) + lmat(k, 729)
         mat(k, 735) = lmat(k, 735)
         mat(k, 736) = mat(k, 736) + lmat(k, 736)
         mat(k, 740) = lmat(k, 740)
         mat(k, 741) = lmat(k, 741)
         mat(k, 743) = lmat(k, 743)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 745) = lmat(k, 745)
         mat(k, 746) = mat(k, 746) + lmat(k, 746)
         mat(k, 747) = lmat(k, 747)
         mat(k, 748) = lmat(k, 748)
         mat(k, 749) = lmat(k, 749)
         mat(k, 750) = lmat(k, 750)
         mat(k, 751) = lmat(k, 751)
         mat(k, 753) = lmat(k, 753)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 755) = lmat(k, 755)
         mat(k, 756) = lmat(k, 756)
         mat(k, 757) = mat(k, 757) + lmat(k, 757)
         mat(k, 760) = mat(k, 760) + lmat(k, 760)
         mat(k, 761) = mat(k, 761) + lmat(k, 761)
         mat(k, 763) = mat(k, 763) + lmat(k, 763)
         mat(k, 765) = lmat(k, 765)
         mat(k, 766) = mat(k, 766) + lmat(k, 766)
         mat(k, 768) = mat(k, 768) + lmat(k, 768)
         mat(k, 775) = mat(k, 775) + lmat(k, 775)
         mat(k, 780) = lmat(k, 780)
         mat(k, 781) = lmat(k, 781)
         mat(k, 782) = lmat(k, 782)
         mat(k, 783) = lmat(k, 783)
         mat(k, 784) = mat(k, 784) + lmat(k, 784)
         mat(k, 789) = lmat(k, 789)
         mat(k, 791) = lmat(k, 791)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 794) = lmat(k, 794)
         mat(k, 797) = mat(k, 797) + lmat(k, 797)
         mat(k, 805) = mat(k, 805) + lmat(k, 805)
         mat(k, 807) = lmat(k, 807)
         mat(k, 808) = lmat(k, 808)
         mat(k, 809) = lmat(k, 809)
         mat(k, 810) = mat(k, 810) + lmat(k, 810)
         mat(k, 811) = lmat(k, 811)
         mat(k, 812) = mat(k, 812) + lmat(k, 812)
         mat(k, 823) = mat(k, 823) + lmat(k, 823)
         mat(k, 839) = mat(k, 839) + lmat(k, 839)
         mat(k, 850) = mat(k, 850) + lmat(k, 850)
         mat(k, 859) = mat(k, 859) + lmat(k, 859)
         mat(k, 868) = mat(k, 868) + lmat(k, 868)
         mat(k, 872) = lmat(k, 872)
         mat(k, 873) = lmat(k, 873)
         mat(k, 874) = lmat(k, 874)
         mat(k, 878) = mat(k, 878) + lmat(k, 878)
         mat(k, 886) = lmat(k, 886)
         mat(k, 887) = lmat(k, 887)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 890) = lmat(k, 890)
         mat(k, 891) = mat(k, 891) + lmat(k, 891)
         mat(k, 893) = lmat(k, 893)
         mat(k, 894) = lmat(k, 894)
         mat(k, 895) = mat(k, 895) + lmat(k, 895)
         mat(k, 897) = lmat(k, 897)
         mat(k, 898) = mat(k, 898) + lmat(k, 898)
         mat(k, 900) = lmat(k, 900)
         mat(k, 901) = lmat(k, 901)
         mat(k, 903) = mat(k, 903) + lmat(k, 903)
         mat(k, 904) = mat(k, 904) + lmat(k, 904)
         mat(k, 905) = mat(k, 905) + lmat(k, 905)
         mat(k, 906) = mat(k, 906) + lmat(k, 906)
         mat(k, 907) = lmat(k, 907)
         mat(k, 908) = mat(k, 908) + lmat(k, 908)
         mat(k, 911) = lmat(k, 911)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 915) = lmat(k, 915)
         mat(k, 919) = mat(k, 919) + lmat(k, 919)
         mat(k, 936) = mat(k, 936) + lmat(k, 936)
         mat(k, 947) = mat(k, 947) + lmat(k, 947)
         mat(k, 956) = mat(k, 956) + lmat(k, 956)
         mat(k, 974) = mat(k, 974) + lmat(k, 974)
         mat(k, 998) = mat(k, 998) + lmat(k, 998)
         mat(k,1008) = mat(k,1008) + lmat(k,1008)
         mat(k,1009) = mat(k,1009) + lmat(k,1009)
         mat(k,1010) = mat(k,1010) + lmat(k,1010)
         mat(k,1012) = mat(k,1012) + lmat(k,1012)
         mat(k,1013) = mat(k,1013) + lmat(k,1013)
         mat(k,1015) = mat(k,1015) + lmat(k,1015)
         mat(k,1016) = mat(k,1016) + lmat(k,1016)
         mat(k,1018) = lmat(k,1018)
         mat(k,1020) = mat(k,1020) + lmat(k,1020)
         mat(k,1026) = mat(k,1026) + lmat(k,1026)
         mat(k,1039) = mat(k,1039) + lmat(k,1039)
         mat(k,1040) = lmat(k,1040)
         mat(k,1043) = lmat(k,1043)
         mat(k,1045) = lmat(k,1045)
         mat(k,1047) = mat(k,1047) + lmat(k,1047)
         mat(k,1048) = mat(k,1048) + lmat(k,1048)
         mat(k,1050) = mat(k,1050) + lmat(k,1050)
         mat(k,1064) = mat(k,1064) + lmat(k,1064)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1098) = mat(k,1098) + lmat(k,1098)
         mat(k,1100) = lmat(k,1100)
         mat(k,1101) = lmat(k,1101)
         mat(k,1105) = lmat(k,1105)
         mat(k,1106) = mat(k,1106) + lmat(k,1106)
         mat(k,1108) = lmat(k,1108)
         mat(k,1109) = lmat(k,1109)
         mat(k,1113) = mat(k,1113) + lmat(k,1113)
         mat(k,1129) = mat(k,1129) + lmat(k,1129)
         mat(k,1146) = lmat(k,1146)
         mat(k,1150) = mat(k,1150) + lmat(k,1150)
         mat(k,1157) = lmat(k,1157)
         mat(k,1158) = mat(k,1158) + lmat(k,1158)
         mat(k,1160) = lmat(k,1160)
         mat(k,1166) = mat(k,1166) + lmat(k,1166)
         mat(k,1178) = mat(k,1178) + lmat(k,1178)
         mat(k,1179) = mat(k,1179) + lmat(k,1179)
         mat(k,1180) = mat(k,1180) + lmat(k,1180)
         mat(k,1181) = mat(k,1181) + lmat(k,1181)
         mat(k,1182) = mat(k,1182) + lmat(k,1182)
         mat(k,1183) = mat(k,1183) + lmat(k,1183)
         mat(k,1185) = mat(k,1185) + lmat(k,1185)
         mat(k,1187) = mat(k,1187) + lmat(k,1187)
         mat(k,1190) = mat(k,1190) + lmat(k,1190)
         mat(k,1191) = lmat(k,1191)
         mat(k,1193) = lmat(k,1193)
         mat(k,1197) = lmat(k,1197)
         mat(k,1200) = mat(k,1200) + lmat(k,1200)
         mat(k,1206) = lmat(k,1206)
         mat(k,1207) = mat(k,1207) + lmat(k,1207)
         mat(k,1208) = mat(k,1208) + lmat(k,1208)
         mat(k,1210) = mat(k,1210) + lmat(k,1210)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1233) = lmat(k,1233)
         mat(k,1234) = lmat(k,1234)
         mat(k,1235) = lmat(k,1235)
         mat(k,1236) = lmat(k,1236)
         mat(k,1237) = mat(k,1237) + lmat(k,1237)
         mat(k,1238) = lmat(k,1238)
         mat(k,1240) = lmat(k,1240)
         mat(k,1242) = lmat(k,1242)
         mat(k,1245) = lmat(k,1245)
         mat(k,1246) = mat(k,1246) + lmat(k,1246)
         mat(k,1248) = lmat(k,1248)
         mat(k,1250) = mat(k,1250) + lmat(k,1250)
         mat(k,1252) = lmat(k,1252)
         mat(k,1253) = lmat(k,1253)
         mat(k,1255) = mat(k,1255) + lmat(k,1255)
         mat(k,1258) = mat(k,1258) + lmat(k,1258)
         mat(k,1267) = mat(k,1267) + lmat(k,1267)
         mat(k,1269) = lmat(k,1269)
         mat(k,1274) = mat(k,1274) + lmat(k,1274)
         mat(k,1291) = mat(k,1291) + lmat(k,1291)
         mat(k,1311) = mat(k,1311) + lmat(k,1311)
         mat(k,1326) = mat(k,1326) + lmat(k,1326)
         mat(k,1327) = mat(k,1327) + lmat(k,1327)
         mat(k,1330) = mat(k,1330) + lmat(k,1330)
         mat(k,1331) = mat(k,1331) + lmat(k,1331)
         mat(k,1332) = mat(k,1332) + lmat(k,1332)
         mat(k,1335) = mat(k,1335) + lmat(k,1335)
         mat(k,1338) = mat(k,1338) + lmat(k,1338)
         mat(k,1339) = mat(k,1339) + lmat(k,1339)
         mat(k,1340) = mat(k,1340) + lmat(k,1340)
         mat(k,1344) = lmat(k,1344)
         mat(k,1359) = mat(k,1359) + lmat(k,1359)
         mat(k,1375) = lmat(k,1375)
         mat(k,1392) = mat(k,1392) + lmat(k,1392)
         mat(k,1403) = mat(k,1403) + lmat(k,1403)
         mat(k,1416) = mat(k,1416) + lmat(k,1416)
         mat(k,1431) = lmat(k,1431)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1437) = mat(k,1437) + lmat(k,1437)
         mat(k,1439) = mat(k,1439) + lmat(k,1439)
         mat(k,1441) = lmat(k,1441)
         mat(k,1459) = mat(k,1459) + lmat(k,1459)
         mat(k,1491) = mat(k,1491) + lmat(k,1491)
         mat(k,1507) = mat(k,1507) + lmat(k,1507)
         mat(k,1517) = mat(k,1517) + lmat(k,1517)
         mat(k,1518) = mat(k,1518) + lmat(k,1518)
         mat(k,1521) = mat(k,1521) + lmat(k,1521)
         mat(k,1528) = mat(k,1528) + lmat(k,1528)
         mat(k,1534) = mat(k,1534) + lmat(k,1534)
         mat(k,1547) = lmat(k,1547)
         mat(k,1548) = mat(k,1548) + lmat(k,1548)
         mat(k,1549) = mat(k,1549) + lmat(k,1549)
         mat(k,1555) = lmat(k,1555)
         mat(k,1563) = lmat(k,1563)
         mat(k,1565) = lmat(k,1565)
         mat(k,1566) = mat(k,1566) + lmat(k,1566)
         mat(k,1567) = mat(k,1567) + lmat(k,1567)
         mat(k,1568) = mat(k,1568) + lmat(k,1568)
         mat(k,1569) = mat(k,1569) + lmat(k,1569)
         mat(k,1574) = mat(k,1574) + lmat(k,1574)
         mat(k,1577) = lmat(k,1577)
         mat(k,1578) = mat(k,1578) + lmat(k,1578)
         mat(k,1579) = mat(k,1579) + lmat(k,1579)
         mat(k,1580) = mat(k,1580) + lmat(k,1580)
         mat(k,1581) = mat(k,1581) + lmat(k,1581)
         mat(k,1585) = mat(k,1585) + lmat(k,1585)
         mat(k,1592) = mat(k,1592) + lmat(k,1592)
         mat(k,1595) = lmat(k,1595)
         mat(k,1639) = mat(k,1639) + lmat(k,1639)
         mat(k,1657) = mat(k,1657) + lmat(k,1657)
         mat(k,1658) = mat(k,1658) + lmat(k,1658)
         mat(k,1663) = mat(k,1663) + lmat(k,1663)
         mat(k,1669) = mat(k,1669) + lmat(k,1669)
         mat(k,1675) = lmat(k,1675)
         mat(k,1681) = mat(k,1681) + lmat(k,1681)
         mat(k,1683) = lmat(k,1683)
         mat(k,1688) = mat(k,1688) + lmat(k,1688)
         mat(k,1704) = mat(k,1704) + lmat(k,1704)
         mat(k,1708) = mat(k,1708) + lmat(k,1708)
         mat(k,1716) = mat(k,1716) + lmat(k,1716)
         mat(k,1719) = mat(k,1719) + lmat(k,1719)
         mat(k,1721) = mat(k,1721) + lmat(k,1721)
         mat(k,1740) = mat(k,1740) + lmat(k,1740)
         mat(k,1742) = mat(k,1742) + lmat(k,1742)
         mat(k,1752) = mat(k,1752) + lmat(k,1752)
         mat(k,1808) = mat(k,1808) + lmat(k,1808)
         mat(k,1813) = mat(k,1813) + lmat(k,1813)
         mat(k,1818) = mat(k,1818) + lmat(k,1818)
         mat(k,1822) = mat(k,1822) + lmat(k,1822)
         mat(k,1823) = mat(k,1823) + lmat(k,1823)
         mat(k,1869) = mat(k,1869) + lmat(k,1869)
         mat(k,1936) = mat(k,1936) + lmat(k,1936)
         mat(k,1937) = mat(k,1937) + lmat(k,1937)
         mat(k,1944) = mat(k,1944) + lmat(k,1944)
         mat(k,1983) = mat(k,1983) + lmat(k,1983)
         mat(k,1990) = mat(k,1990) + lmat(k,1990)
         mat(k,2166) = mat(k,2166) + lmat(k,2166)
         mat(k,2186) = lmat(k,2186)
         mat(k,2191) = mat(k,2191) + lmat(k,2191)
         mat(k,2192) = mat(k,2192) + lmat(k,2192)
         mat(k,2201) = lmat(k,2201)
         mat(k,2224) = mat(k,2224) + lmat(k,2224)
         mat(k,2227) = mat(k,2227) + lmat(k,2227)
         mat(k,2228) = mat(k,2228) + lmat(k,2228)
         mat(k,2340) = mat(k,2340) + lmat(k,2340)
         mat(k,2349) = mat(k,2349) + lmat(k,2349)
         mat(k,2363) = mat(k,2363) + lmat(k,2363)
         mat(k,2370) = mat(k,2370) + lmat(k,2370)
         mat(k,2372) = mat(k,2372) + lmat(k,2372)
         mat(k,2399) = mat(k,2399) + lmat(k,2399)
         mat(k,2427) = mat(k,2427) + lmat(k,2427)
         mat(k,2435) = mat(k,2435) + lmat(k,2435)
         mat(k,2462) = mat(k,2462) + lmat(k,2462)
         mat(k,2463) = mat(k,2463) + lmat(k,2463)
         mat(k,2464) = mat(k,2464) + lmat(k,2464)
         mat(k,2492) = mat(k,2492) + lmat(k,2492)
         mat(k,2510) = mat(k,2510) + lmat(k,2510)
         mat(k,2516) = mat(k,2516) + lmat(k,2516)
         mat(k,2564) = mat(k,2564) + lmat(k,2564)
         mat(k,2565) = mat(k,2565) + lmat(k,2565)
         mat(k,2570) = mat(k,2570) + lmat(k,2570)
         mat(k,2574) = mat(k,2574) + lmat(k,2574)
         mat(k,2575) = mat(k,2575) + lmat(k,2575)
         mat(k,2602) = mat(k,2602) + lmat(k,2602)
         mat(k,2671) = mat(k,2671) + lmat(k,2671)
         mat(k,2676) = mat(k,2676) + lmat(k,2676)
         mat(k,2683) = lmat(k,2683)
         mat(k,2693) = mat(k,2693) + lmat(k,2693)
         mat(k,2694) = mat(k,2694) + lmat(k,2694)
         mat(k,2700) = lmat(k,2700)
         mat(k,2703) = lmat(k,2703)
         mat(k,2706) = mat(k,2706) + lmat(k,2706)
         mat(k, 246) = 0._r8
         mat(k, 247) = 0._r8
         mat(k, 316) = 0._r8
         mat(k, 384) = 0._r8
         mat(k, 398) = 0._r8
         mat(k, 515) = 0._r8
         mat(k, 516) = 0._r8
         mat(k, 529) = 0._r8
         mat(k, 569) = 0._r8
         mat(k, 571) = 0._r8
         mat(k, 577) = 0._r8
         mat(k, 587) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 730) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 759) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 785) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 788) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 798) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 822) = 0._r8
         mat(k, 824) = 0._r8
         mat(k, 825) = 0._r8
         mat(k, 827) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 838) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 841) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 863) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 882) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 896) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 975) = 0._r8
         mat(k, 977) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 992) = 0._r8
         mat(k, 996) = 0._r8
         mat(k, 997) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1021) = 0._r8
         mat(k,1024) = 0._r8
         mat(k,1044) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1061) = 0._r8
         mat(k,1063) = 0._r8
         mat(k,1070) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1081) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1083) = 0._r8
         mat(k,1084) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1132) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1211) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1222) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1241) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1320) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1351) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1360) = 0._r8
         mat(k,1362) = 0._r8
         mat(k,1367) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1374) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1384) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1390) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1394) = 0._r8
         mat(k,1395) = 0._r8
         mat(k,1400) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1402) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1407) = 0._r8
         mat(k,1414) = 0._r8
         mat(k,1417) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1466) = 0._r8
         mat(k,1468) = 0._r8
         mat(k,1471) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1503) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1524) = 0._r8
         mat(k,1526) = 0._r8
         mat(k,1530) = 0._r8
         mat(k,1531) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1536) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1571) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1586) = 0._r8
         mat(k,1587) = 0._r8
         mat(k,1590) = 0._r8
         mat(k,1591) = 0._r8
         mat(k,1596) = 0._r8
         mat(k,1597) = 0._r8
         mat(k,1598) = 0._r8
         mat(k,1599) = 0._r8
         mat(k,1600) = 0._r8
         mat(k,1612) = 0._r8
         mat(k,1638) = 0._r8
         mat(k,1641) = 0._r8
         mat(k,1642) = 0._r8
         mat(k,1643) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1646) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1648) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1652) = 0._r8
         mat(k,1655) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1662) = 0._r8
         mat(k,1664) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1667) = 0._r8
         mat(k,1668) = 0._r8
         mat(k,1670) = 0._r8
         mat(k,1671) = 0._r8
         mat(k,1672) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1680) = 0._r8
         mat(k,1684) = 0._r8
         mat(k,1685) = 0._r8
         mat(k,1686) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1693) = 0._r8
         mat(k,1694) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1699) = 0._r8
         mat(k,1700) = 0._r8
         mat(k,1702) = 0._r8
         mat(k,1705) = 0._r8
         mat(k,1706) = 0._r8
         mat(k,1713) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1718) = 0._r8
         mat(k,1720) = 0._r8
         mat(k,1722) = 0._r8
         mat(k,1723) = 0._r8
         mat(k,1725) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1729) = 0._r8
         mat(k,1730) = 0._r8
         mat(k,1741) = 0._r8
         mat(k,1746) = 0._r8
         mat(k,1748) = 0._r8
         mat(k,1750) = 0._r8
         mat(k,1751) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1761) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1773) = 0._r8
         mat(k,1780) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1788) = 0._r8
         mat(k,1800) = 0._r8
         mat(k,1801) = 0._r8
         mat(k,1803) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1807) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1811) = 0._r8
         mat(k,1816) = 0._r8
         mat(k,1819) = 0._r8
         mat(k,1820) = 0._r8
         mat(k,1821) = 0._r8
         mat(k,1824) = 0._r8
         mat(k,1843) = 0._r8
         mat(k,1848) = 0._r8
         mat(k,1851) = 0._r8
         mat(k,1852) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1856) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1858) = 0._r8
         mat(k,1859) = 0._r8
         mat(k,1861) = 0._r8
         mat(k,1864) = 0._r8
         mat(k,1871) = 0._r8
         mat(k,1873) = 0._r8
         mat(k,1874) = 0._r8
         mat(k,1878) = 0._r8
         mat(k,1879) = 0._r8
         mat(k,1883) = 0._r8
         mat(k,1884) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1902) = 0._r8
         mat(k,1905) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1915) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1919) = 0._r8
         mat(k,1922) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1924) = 0._r8
         mat(k,1926) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1939) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1989) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2092) = 0._r8
         mat(k,2104) = 0._r8
         mat(k,2107) = 0._r8
         mat(k,2115) = 0._r8
         mat(k,2116) = 0._r8
         mat(k,2142) = 0._r8
         mat(k,2165) = 0._r8
         mat(k,2181) = 0._r8
         mat(k,2182) = 0._r8
         mat(k,2183) = 0._r8
         mat(k,2184) = 0._r8
         mat(k,2185) = 0._r8
         mat(k,2188) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2190) = 0._r8
         mat(k,2193) = 0._r8
         mat(k,2194) = 0._r8
         mat(k,2195) = 0._r8
         mat(k,2196) = 0._r8
         mat(k,2197) = 0._r8
         mat(k,2198) = 0._r8
         mat(k,2199) = 0._r8
         mat(k,2200) = 0._r8
         mat(k,2202) = 0._r8
         mat(k,2212) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2215) = 0._r8
         mat(k,2216) = 0._r8
         mat(k,2217) = 0._r8
         mat(k,2221) = 0._r8
         mat(k,2223) = 0._r8
         mat(k,2231) = 0._r8
         mat(k,2234) = 0._r8
         mat(k,2250) = 0._r8
         mat(k,2268) = 0._r8
         mat(k,2269) = 0._r8
         mat(k,2277) = 0._r8
         mat(k,2294) = 0._r8
         mat(k,2301) = 0._r8
         mat(k,2302) = 0._r8
         mat(k,2305) = 0._r8
         mat(k,2308) = 0._r8
         mat(k,2310) = 0._r8
         mat(k,2315) = 0._r8
         mat(k,2320) = 0._r8
         mat(k,2324) = 0._r8
         mat(k,2326) = 0._r8
         mat(k,2336) = 0._r8
         mat(k,2338) = 0._r8
         mat(k,2362) = 0._r8
         mat(k,2364) = 0._r8
         mat(k,2365) = 0._r8
         mat(k,2367) = 0._r8
         mat(k,2375) = 0._r8
         mat(k,2378) = 0._r8
         mat(k,2381) = 0._r8
         mat(k,2382) = 0._r8
         mat(k,2386) = 0._r8
         mat(k,2387) = 0._r8
         mat(k,2388) = 0._r8
         mat(k,2389) = 0._r8
         mat(k,2391) = 0._r8
         mat(k,2393) = 0._r8
         mat(k,2394) = 0._r8
         mat(k,2395) = 0._r8
         mat(k,2398) = 0._r8
         mat(k,2400) = 0._r8
         mat(k,2403) = 0._r8
         mat(k,2406) = 0._r8
         mat(k,2408) = 0._r8
         mat(k,2412) = 0._r8
         mat(k,2418) = 0._r8
         mat(k,2419) = 0._r8
         mat(k,2420) = 0._r8
         mat(k,2428) = 0._r8
         mat(k,2430) = 0._r8
         mat(k,2441) = 0._r8
         mat(k,2450) = 0._r8
         mat(k,2451) = 0._r8
         mat(k,2452) = 0._r8
         mat(k,2454) = 0._r8
         mat(k,2455) = 0._r8
         mat(k,2457) = 0._r8
         mat(k,2465) = 0._r8
         mat(k,2468) = 0._r8
         mat(k,2471) = 0._r8
         mat(k,2475) = 0._r8
         mat(k,2476) = 0._r8
         mat(k,2478) = 0._r8
         mat(k,2479) = 0._r8
         mat(k,2481) = 0._r8
         mat(k,2483) = 0._r8
         mat(k,2484) = 0._r8
         mat(k,2485) = 0._r8
         mat(k,2488) = 0._r8
         mat(k,2490) = 0._r8
         mat(k,2493) = 0._r8
         mat(k,2495) = 0._r8
         mat(k,2496) = 0._r8
         mat(k,2498) = 0._r8
         mat(k,2499) = 0._r8
         mat(k,2500) = 0._r8
         mat(k,2501) = 0._r8
         mat(k,2502) = 0._r8
         mat(k,2503) = 0._r8
         mat(k,2504) = 0._r8
         mat(k,2506) = 0._r8
         mat(k,2508) = 0._r8
         mat(k,2509) = 0._r8
         mat(k,2511) = 0._r8
         mat(k,2512) = 0._r8
         mat(k,2514) = 0._r8
         mat(k,2515) = 0._r8
         mat(k,2517) = 0._r8
         mat(k,2518) = 0._r8
         mat(k,2525) = 0._r8
         mat(k,2527) = 0._r8
         mat(k,2529) = 0._r8
         mat(k,2532) = 0._r8
         mat(k,2536) = 0._r8
         mat(k,2539) = 0._r8
         mat(k,2545) = 0._r8
         mat(k,2546) = 0._r8
         mat(k,2547) = 0._r8
         mat(k,2548) = 0._r8
         mat(k,2551) = 0._r8
         mat(k,2553) = 0._r8
         mat(k,2554) = 0._r8
         mat(k,2555) = 0._r8
         mat(k,2556) = 0._r8
         mat(k,2557) = 0._r8
         mat(k,2558) = 0._r8
         mat(k,2559) = 0._r8
         mat(k,2563) = 0._r8
         mat(k,2573) = 0._r8
         mat(k,2576) = 0._r8
         mat(k,2623) = 0._r8
         mat(k,2653) = 0._r8
         mat(k,2654) = 0._r8
         mat(k,2655) = 0._r8
         mat(k,2657) = 0._r8
         mat(k,2659) = 0._r8
         mat(k,2660) = 0._r8
         mat(k,2664) = 0._r8
         mat(k,2666) = 0._r8
         mat(k,2674) = 0._r8
         mat(k,2677) = 0._r8
         mat(k,2682) = 0._r8
         mat(k,2684) = 0._r8
         mat(k,2685) = 0._r8
         mat(k,2686) = 0._r8
         mat(k,2687) = 0._r8
         mat(k,2688) = 0._r8
         mat(k,2689) = 0._r8
         mat(k,2690) = 0._r8
         mat(k,2691) = 0._r8
         mat(k,2692) = 0._r8
         mat(k,2695) = 0._r8
         mat(k,2696) = 0._r8
         mat(k,2697) = 0._r8
         mat(k,2698) = 0._r8
         mat(k,2699) = 0._r8
         mat(k,2701) = 0._r8
         mat(k,2702) = 0._r8
         mat(k,2704) = 0._r8
         mat(k,2705) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 2) = mat(k, 2) - dti(k)
         mat(k, 3) = mat(k, 3) - dti(k)
         mat(k, 4) = mat(k, 4) - dti(k)
         mat(k, 5) = mat(k, 5) - dti(k)
         mat(k, 7) = mat(k, 7) - dti(k)
         mat(k, 10) = mat(k, 10) - dti(k)
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
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 63) = mat(k, 63) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 81) = mat(k, 81) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 89) = mat(k, 89) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 104) = mat(k, 104) - dti(k)
         mat(k, 107) = mat(k, 107) - dti(k)
         mat(k, 110) = mat(k, 110) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 120) = mat(k, 120) - dti(k)
         mat(k, 124) = mat(k, 124) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 132) = mat(k, 132) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 140) = mat(k, 140) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 148) = mat(k, 148) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 158) = mat(k, 158) - dti(k)
         mat(k, 161) = mat(k, 161) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 167) = mat(k, 167) - dti(k)
         mat(k, 172) = mat(k, 172) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 182) = mat(k, 182) - dti(k)
         mat(k, 188) = mat(k, 188) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 198) = mat(k, 198) - dti(k)
         mat(k, 203) = mat(k, 203) - dti(k)
         mat(k, 205) = mat(k, 205) - dti(k)
         mat(k, 209) = mat(k, 209) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 233) = mat(k, 233) - dti(k)
         mat(k, 237) = mat(k, 237) - dti(k)
         mat(k, 242) = mat(k, 242) - dti(k)
         mat(k, 250) = mat(k, 250) - dti(k)
         mat(k, 255) = mat(k, 255) - dti(k)
         mat(k, 260) = mat(k, 260) - dti(k)
         mat(k, 264) = mat(k, 264) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 272) = mat(k, 272) - dti(k)
         mat(k, 277) = mat(k, 277) - dti(k)
         mat(k, 280) = mat(k, 280) - dti(k)
         mat(k, 283) = mat(k, 283) - dti(k)
         mat(k, 286) = mat(k, 286) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 296) = mat(k, 296) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 304) = mat(k, 304) - dti(k)
         mat(k, 308) = mat(k, 308) - dti(k)
         mat(k, 315) = mat(k, 315) - dti(k)
         mat(k, 321) = mat(k, 321) - dti(k)
         mat(k, 327) = mat(k, 327) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 342) = mat(k, 342) - dti(k)
         mat(k, 348) = mat(k, 348) - dti(k)
         mat(k, 351) = mat(k, 351) - dti(k)
         mat(k, 354) = mat(k, 354) - dti(k)
         mat(k, 359) = mat(k, 359) - dti(k)
         mat(k, 364) = mat(k, 364) - dti(k)
         mat(k, 369) = mat(k, 369) - dti(k)
         mat(k, 376) = mat(k, 376) - dti(k)
         mat(k, 382) = mat(k, 382) - dti(k)
         mat(k, 387) = mat(k, 387) - dti(k)
         mat(k, 392) = mat(k, 392) - dti(k)
         mat(k, 396) = mat(k, 396) - dti(k)
         mat(k, 403) = mat(k, 403) - dti(k)
         mat(k, 411) = mat(k, 411) - dti(k)
         mat(k, 419) = mat(k, 419) - dti(k)
         mat(k, 427) = mat(k, 427) - dti(k)
         mat(k, 430) = mat(k, 430) - dti(k)
         mat(k, 438) = mat(k, 438) - dti(k)
         mat(k, 446) = mat(k, 446) - dti(k)
         mat(k, 452) = mat(k, 452) - dti(k)
         mat(k, 458) = mat(k, 458) - dti(k)
         mat(k, 464) = mat(k, 464) - dti(k)
         mat(k, 470) = mat(k, 470) - dti(k)
         mat(k, 476) = mat(k, 476) - dti(k)
         mat(k, 482) = mat(k, 482) - dti(k)
         mat(k, 488) = mat(k, 488) - dti(k)
         mat(k, 494) = mat(k, 494) - dti(k)
         mat(k, 500) = mat(k, 500) - dti(k)
         mat(k, 506) = mat(k, 506) - dti(k)
         mat(k, 514) = mat(k, 514) - dti(k)
         mat(k, 520) = mat(k, 520) - dti(k)
         mat(k, 527) = mat(k, 527) - dti(k)
         mat(k, 533) = mat(k, 533) - dti(k)
         mat(k, 538) = mat(k, 538) - dti(k)
         mat(k, 541) = mat(k, 541) - dti(k)
         mat(k, 545) = mat(k, 545) - dti(k)
         mat(k, 549) = mat(k, 549) - dti(k)
         mat(k, 552) = mat(k, 552) - dti(k)
         mat(k, 559) = mat(k, 559) - dti(k)
         mat(k, 568) = mat(k, 568) - dti(k)
         mat(k, 576) = mat(k, 576) - dti(k)
         mat(k, 584) = mat(k, 584) - dti(k)
         mat(k, 593) = mat(k, 593) - dti(k)
         mat(k, 598) = mat(k, 598) - dti(k)
         mat(k, 603) = mat(k, 603) - dti(k)
         mat(k, 610) = mat(k, 610) - dti(k)
         mat(k, 616) = mat(k, 616) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 625) = mat(k, 625) - dti(k)
         mat(k, 633) = mat(k, 633) - dti(k)
         mat(k, 641) = mat(k, 641) - dti(k)
         mat(k, 649) = mat(k, 649) - dti(k)
         mat(k, 657) = mat(k, 657) - dti(k)
         mat(k, 665) = mat(k, 665) - dti(k)
         mat(k, 674) = mat(k, 674) - dti(k)
         mat(k, 681) = mat(k, 681) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 701) = mat(k, 701) - dti(k)
         mat(k, 710) = mat(k, 710) - dti(k)
         mat(k, 714) = mat(k, 714) - dti(k)
         mat(k, 721) = mat(k, 721) - dti(k)
         mat(k, 729) = mat(k, 729) - dti(k)
         mat(k, 736) = mat(k, 736) - dti(k)
         mat(k, 746) = mat(k, 746) - dti(k)
         mat(k, 757) = mat(k, 757) - dti(k)
         mat(k, 768) = mat(k, 768) - dti(k)
         mat(k, 775) = mat(k, 775) - dti(k)
         mat(k, 784) = mat(k, 784) - dti(k)
         mat(k, 797) = mat(k, 797) - dti(k)
         mat(k, 805) = mat(k, 805) - dti(k)
         mat(k, 812) = mat(k, 812) - dti(k)
         mat(k, 823) = mat(k, 823) - dti(k)
         mat(k, 839) = mat(k, 839) - dti(k)
         mat(k, 850) = mat(k, 850) - dti(k)
         mat(k, 859) = mat(k, 859) - dti(k)
         mat(k, 868) = mat(k, 868) - dti(k)
         mat(k, 872) = mat(k, 872) - dti(k)
         mat(k, 878) = mat(k, 878) - dti(k)
         mat(k, 888) = mat(k, 888) - dti(k)
         mat(k, 898) = mat(k, 898) - dti(k)
         mat(k, 906) = mat(k, 906) - dti(k)
         mat(k, 919) = mat(k, 919) - dti(k)
         mat(k, 936) = mat(k, 936) - dti(k)
         mat(k, 947) = mat(k, 947) - dti(k)
         mat(k, 956) = mat(k, 956) - dti(k)
         mat(k, 974) = mat(k, 974) - dti(k)
         mat(k, 998) = mat(k, 998) - dti(k)
         mat(k,1009) = mat(k,1009) - dti(k)
         mat(k,1020) = mat(k,1020) - dti(k)
         mat(k,1026) = mat(k,1026) - dti(k)
         mat(k,1039) = mat(k,1039) - dti(k)
         mat(k,1047) = mat(k,1047) - dti(k)
         mat(k,1064) = mat(k,1064) - dti(k)
         mat(k,1088) = mat(k,1088) - dti(k)
         mat(k,1098) = mat(k,1098) - dti(k)
         mat(k,1106) = mat(k,1106) - dti(k)
         mat(k,1113) = mat(k,1113) - dti(k)
         mat(k,1129) = mat(k,1129) - dti(k)
         mat(k,1150) = mat(k,1150) - dti(k)
         mat(k,1166) = mat(k,1166) - dti(k)
         mat(k,1180) = mat(k,1180) - dti(k)
         mat(k,1190) = mat(k,1190) - dti(k)
         mat(k,1200) = mat(k,1200) - dti(k)
         mat(k,1207) = mat(k,1207) - dti(k)
         mat(k,1220) = mat(k,1220) - dti(k)
         mat(k,1237) = mat(k,1237) - dti(k)
         mat(k,1250) = mat(k,1250) - dti(k)
         mat(k,1258) = mat(k,1258) - dti(k)
         mat(k,1274) = mat(k,1274) - dti(k)
         mat(k,1291) = mat(k,1291) - dti(k)
         mat(k,1311) = mat(k,1311) - dti(k)
         mat(k,1327) = mat(k,1327) - dti(k)
         mat(k,1339) = mat(k,1339) - dti(k)
         mat(k,1359) = mat(k,1359) - dti(k)
         mat(k,1392) = mat(k,1392) - dti(k)
         mat(k,1416) = mat(k,1416) - dti(k)
         mat(k,1437) = mat(k,1437) - dti(k)
         mat(k,1459) = mat(k,1459) - dti(k)
         mat(k,1491) = mat(k,1491) - dti(k)
         mat(k,1507) = mat(k,1507) - dti(k)
         mat(k,1521) = mat(k,1521) - dti(k)
         mat(k,1534) = mat(k,1534) - dti(k)
         mat(k,1549) = mat(k,1549) - dti(k)
         mat(k,1567) = mat(k,1567) - dti(k)
         mat(k,1585) = mat(k,1585) - dti(k)
         mat(k,1639) = mat(k,1639) - dti(k)
         mat(k,1663) = mat(k,1663) - dti(k)
         mat(k,1688) = mat(k,1688) - dti(k)
         mat(k,1716) = mat(k,1716) - dti(k)
         mat(k,1740) = mat(k,1740) - dti(k)
         mat(k,1808) = mat(k,1808) - dti(k)
         mat(k,1869) = mat(k,1869) - dti(k)
         mat(k,1936) = mat(k,1936) - dti(k)
         mat(k,1983) = mat(k,1983) - dti(k)
         mat(k,2166) = mat(k,2166) - dti(k)
         mat(k,2192) = mat(k,2192) - dti(k)
         mat(k,2224) = mat(k,2224) - dti(k)
         mat(k,2340) = mat(k,2340) - dti(k)
         mat(k,2370) = mat(k,2370) - dti(k)
         mat(k,2399) = mat(k,2399) - dti(k)
         mat(k,2435) = mat(k,2435) - dti(k)
         mat(k,2463) = mat(k,2463) - dti(k)
         mat(k,2492) = mat(k,2492) - dti(k)
         mat(k,2516) = mat(k,2516) - dti(k)
         mat(k,2574) = mat(k,2574) - dti(k)
         mat(k,2676) = mat(k,2676) - dti(k)
         mat(k,2706) = mat(k,2706) - dti(k)
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
