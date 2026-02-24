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
         mat(k,738) = -(rxt(k,450)*y(k,250))
         mat(k,2079) = -rxt(k,450)*y(k,1)
         mat(k,1870) = rxt(k,453)*y(k,221)
         mat(k,1024) = rxt(k,453)*y(k,153)
         mat(k,759) = -(rxt(k,454)*y(k,250))
         mat(k,2081) = -rxt(k,454)*y(k,2)
         mat(k,1025) = rxt(k,451)*y(k,235)
         mat(k,2265) = rxt(k,451)*y(k,221)
         mat(k,1004) = -(rxt(k,533)*y(k,155) + rxt(k,534)*y(k,164) + rxt(k,535) &
                      *y(k,250))
         mat(k,2617) = -rxt(k,533)*y(k,6)
         mat(k,2556) = -rxt(k,534)*y(k,6)
         mat(k,2101) = -rxt(k,535)*y(k,6)
         mat(k,190) = -(rxt(k,492)*y(k,250))
         mat(k,2001) = -rxt(k,492)*y(k,7)
         mat(k,466) = -(rxt(k,495)*y(k,250))
         mat(k,2044) = -rxt(k,495)*y(k,8)
         mat(k,576) = rxt(k,493)*y(k,235)
         mat(k,2242) = rxt(k,493)*y(k,223)
         mat(k,191) = .120_r8*rxt(k,492)*y(k,250)
         mat(k,2002) = .120_r8*rxt(k,492)*y(k,7)
         mat(k,1002) = .100_r8*rxt(k,534)*y(k,164)
         mat(k,974) = .100_r8*rxt(k,537)*y(k,164)
         mat(k,2545) = .100_r8*rxt(k,534)*y(k,6) + .100_r8*rxt(k,537)*y(k,139)
         mat(k,1856) = .500_r8*rxt(k,494)*y(k,223) + .200_r8*rxt(k,521)*y(k,256) &
                      + .060_r8*rxt(k,527)*y(k,259)
         mat(k,577) = .500_r8*rxt(k,494)*y(k,153)
         mat(k,813) = .200_r8*rxt(k,521)*y(k,153)
         mat(k,837) = .060_r8*rxt(k,527)*y(k,153)
         mat(k,1850) = .200_r8*rxt(k,521)*y(k,256) + .200_r8*rxt(k,527)*y(k,259)
         mat(k,812) = .200_r8*rxt(k,521)*y(k,153)
         mat(k,835) = .200_r8*rxt(k,527)*y(k,153)
         mat(k,1867) = .200_r8*rxt(k,521)*y(k,256) + .150_r8*rxt(k,527)*y(k,259)
         mat(k,815) = .200_r8*rxt(k,521)*y(k,153)
         mat(k,838) = .150_r8*rxt(k,527)*y(k,153)
         mat(k,1851) = .210_r8*rxt(k,527)*y(k,259)
         mat(k,836) = .210_r8*rxt(k,527)*y(k,153)
         mat(k,279) = -(rxt(k,455)*y(k,250))
         mat(k,2017) = -rxt(k,455)*y(k,15)
         mat(k,1001) = .050_r8*rxt(k,534)*y(k,164)
         mat(k,973) = .050_r8*rxt(k,537)*y(k,164)
         mat(k,2544) = .050_r8*rxt(k,534)*y(k,6) + .050_r8*rxt(k,537)*y(k,139)
         mat(k,413) = -(rxt(k,421)*y(k,155) + rxt(k,422)*y(k,250))
         mat(k,2609) = -rxt(k,421)*y(k,16)
         mat(k,2037) = -rxt(k,422)*y(k,16)
         mat(k,2535) = -(rxt(k,243)*y(k,51) + rxt(k,244)*y(k,235) + rxt(k,245) &
                      *y(k,154) + rxt(k,246)*y(k,164) + rxt(k,253)*y(k,22) + rxt(k,282) &
                      *y(k,125))
         mat(k,1701) = -rxt(k,243)*y(k,17)
         mat(k,2330) = -rxt(k,244)*y(k,17)
         mat(k,1830) = -rxt(k,245)*y(k,17)
         mat(k,2601) = -rxt(k,246)*y(k,17)
         mat(k,923) = -rxt(k,253)*y(k,17)
         mat(k,2421) = -rxt(k,282)*y(k,17)
         mat(k,546) = rxt(k,242)*y(k,250)
         mat(k,2186) = 4.000_r8*rxt(k,247)*y(k,21) + (rxt(k,248)+rxt(k,249))*y(k,74) &
                      + rxt(k,556)*y(k,83) + rxt(k,272)*y(k,115) + (rxt(k,283) &
                       +rxt(k,284))*y(k,125) + rxt(k,252)*y(k,153) + rxt(k,257) &
                      *y(k,163) + rxt(k,567)*y(k,181) + rxt(k,258)*y(k,250)
         mat(k,173) = rxt(k,232)*y(k,249)
         mat(k,178) = rxt(k,262)*y(k,249)
         mat(k,559) = 2.000_r8*rxt(k,316)*y(k,70) + 2.000_r8*rxt(k,343)*y(k,249) &
                      + 2.000_r8*rxt(k,317)*y(k,250)
         mat(k,149) = rxt(k,318)*y(k,250)
         mat(k,721) = rxt(k,321)*y(k,70) + rxt(k,344)*y(k,249) + rxt(k,322)*y(k,250)
         mat(k,121) = 2.000_r8*rxt(k,328)*y(k,250)
         mat(k,453) = 3.000_r8*rxt(k,329)*y(k,70) + 3.000_r8*rxt(k,263)*y(k,249) &
                      + 3.000_r8*rxt(k,330)*y(k,250)
         mat(k,125) = rxt(k,331)*y(k,250)
         mat(k,2391) = 2.000_r8*rxt(k,316)*y(k,45) + rxt(k,321)*y(k,52) &
                      + 3.000_r8*rxt(k,329)*y(k,66)
         mat(k,2215) = (rxt(k,248)+rxt(k,249))*y(k,21)
         mat(k,1100) = rxt(k,556)*y(k,21)
         mat(k,129) = 2.000_r8*rxt(k,264)*y(k,249)
         mat(k,1515) = rxt(k,259)*y(k,163) + rxt(k,265)*y(k,249) + rxt(k,260)*y(k,250)
         mat(k,2449) = rxt(k,272)*y(k,21)
         mat(k,2421) = mat(k,2421) + (rxt(k,283)+rxt(k,284))*y(k,21)
         mat(k,1930) = rxt(k,252)*y(k,21)
         mat(k,2507) = rxt(k,257)*y(k,21) + rxt(k,259)*y(k,97)
         mat(k,1556) = rxt(k,567)*y(k,21)
         mat(k,1976) = rxt(k,232)*y(k,38) + rxt(k,262)*y(k,39) + 2.000_r8*rxt(k,343) &
                      *y(k,45) + rxt(k,344)*y(k,52) + 3.000_r8*rxt(k,263)*y(k,66) &
                      + 2.000_r8*rxt(k,264)*y(k,94) + rxt(k,265)*y(k,97)
         mat(k,2158) = rxt(k,242)*y(k,18) + rxt(k,258)*y(k,21) + 2.000_r8*rxt(k,317) &
                      *y(k,45) + rxt(k,318)*y(k,46) + rxt(k,322)*y(k,52) &
                      + 2.000_r8*rxt(k,328)*y(k,65) + 3.000_r8*rxt(k,330)*y(k,66) &
                      + rxt(k,331)*y(k,67) + rxt(k,260)*y(k,97)
         mat(k,543) = -(rxt(k,242)*y(k,250))
         mat(k,2054) = -rxt(k,242)*y(k,18)
         mat(k,2512) = rxt(k,253)*y(k,22)
         mat(k,913) = rxt(k,253)*y(k,17)
         mat(k,1502) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,108)
         mat(k,1578) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,97)
         mat(k,2164) = rxt(k,250)*y(k,74)
         mat(k,914) = rxt(k,254)*y(k,70)
         mat(k,2347) = rxt(k,254)*y(k,22)
         mat(k,2194) = rxt(k,250)*y(k,21)
         mat(k,1503) = (rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,109)
         mat(k,1707) = (rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,108)
         mat(k,1579) = (rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,101)
         mat(k,1730) = (rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,97)
         mat(k,2511) = rxt(k,245)*y(k,154)
         mat(k,1780) = rxt(k,245)*y(k,17)
         mat(k,2178) = -(4._r8*rxt(k,247)*y(k,21) + (rxt(k,248) + rxt(k,249) + rxt(k,250) &
                      ) * y(k,74) + rxt(k,251)*y(k,235) + rxt(k,252)*y(k,153) &
                      + rxt(k,255)*y(k,154) + rxt(k,257)*y(k,163) + rxt(k,258) &
                      *y(k,250) + rxt(k,272)*y(k,115) + (rxt(k,283) + rxt(k,284) &
                      ) * y(k,125) + rxt(k,556)*y(k,83) + rxt(k,567)*y(k,181))
         mat(k,2207) = -(rxt(k,248) + rxt(k,249) + rxt(k,250)) * y(k,21)
         mat(k,2322) = -rxt(k,251)*y(k,21)
         mat(k,1922) = -rxt(k,252)*y(k,21)
         mat(k,1822) = -rxt(k,255)*y(k,21)
         mat(k,2499) = -rxt(k,257)*y(k,21)
         mat(k,2150) = -rxt(k,258)*y(k,21)
         mat(k,2441) = -rxt(k,272)*y(k,21)
         mat(k,2413) = -(rxt(k,283) + rxt(k,284)) * y(k,21)
         mat(k,1095) = -rxt(k,556)*y(k,21)
         mat(k,1550) = -rxt(k,567)*y(k,21)
         mat(k,2527) = rxt(k,282)*y(k,125) + rxt(k,246)*y(k,164)
         mat(k,920) = rxt(k,256)*y(k,163)
         mat(k,1510) = rxt(k,266)*y(k,249)
         mat(k,1591) = rxt(k,261)*y(k,163)
         mat(k,2413) = mat(k,2413) + rxt(k,282)*y(k,17)
         mat(k,2499) = mat(k,2499) + rxt(k,256)*y(k,22) + rxt(k,261)*y(k,108)
         mat(k,2593) = rxt(k,246)*y(k,17)
         mat(k,1968) = rxt(k,266)*y(k,97)
         mat(k,915) = -(rxt(k,253)*y(k,17) + rxt(k,254)*y(k,70) + rxt(k,256)*y(k,163))
         mat(k,2514) = -rxt(k,253)*y(k,22)
         mat(k,2354) = -rxt(k,254)*y(k,22)
         mat(k,2479) = -rxt(k,256)*y(k,22)
         mat(k,2165) = rxt(k,255)*y(k,154)
         mat(k,1797) = rxt(k,255)*y(k,21)
         mat(k,282) = -(rxt(k,496)*y(k,250))
         mat(k,2018) = -rxt(k,496)*y(k,24)
         mat(k,1847) = rxt(k,499)*y(k,225)
         mat(k,508) = rxt(k,499)*y(k,153)
         mat(k,376) = -(rxt(k,498)*y(k,250))
         mat(k,2031) = -rxt(k,498)*y(k,25)
         mat(k,509) = rxt(k,497)*y(k,235)
         mat(k,2237) = rxt(k,497)*y(k,225)
         mat(k,216) = -(rxt(k,312)*y(k,70) + rxt(k,313)*y(k,250))
         mat(k,2334) = -rxt(k,312)*y(k,26)
         mat(k,2005) = -rxt(k,313)*y(k,26)
         mat(k,326) = -(rxt(k,369)*y(k,70) + rxt(k,370)*y(k,250))
         mat(k,2337) = -rxt(k,369)*y(k,27)
         mat(k,2025) = -rxt(k,370)*y(k,27)
         mat(k,651) = -(rxt(k,371)*y(k,70) + rxt(k,372)*y(k,164) + rxt(k,397)*y(k,250))
         mat(k,2349) = -rxt(k,371)*y(k,28)
         mat(k,2549) = -rxt(k,372)*y(k,28)
         mat(k,2068) = -rxt(k,397)*y(k,28)
         mat(k,288) = -(rxt(k,314)*y(k,70) + rxt(k,315)*y(k,250))
         mat(k,2336) = -rxt(k,314)*y(k,29)
         mat(k,2020) = -rxt(k,315)*y(k,29)
         mat(k,298) = -(rxt(k,377)*y(k,250))
         mat(k,2022) = -rxt(k,377)*y(k,30)
         mat(k,877) = .800_r8*rxt(k,373)*y(k,226) + .200_r8*rxt(k,374)*y(k,230)
         mat(k,1600) = .200_r8*rxt(k,374)*y(k,226)
         mat(k,386) = -(rxt(k,378)*y(k,250))
         mat(k,2033) = -rxt(k,378)*y(k,31)
         mat(k,878) = rxt(k,375)*y(k,235)
         mat(k,2239) = rxt(k,375)*y(k,226)
         mat(k,335) = -(rxt(k,379)*y(k,70) + rxt(k,380)*y(k,250))
         mat(k,2338) = -rxt(k,379)*y(k,32)
         mat(k,2026) = -rxt(k,380)*y(k,32)
         mat(k,1126) = -(rxt(k,400)*y(k,155) + rxt(k,401)*y(k,164) + rxt(k,419) &
                      *y(k,250))
         mat(k,2627) = -rxt(k,400)*y(k,33)
         mat(k,2564) = -rxt(k,401)*y(k,33)
         mat(k,2112) = -rxt(k,419)*y(k,33)
         mat(k,893) = .130_r8*rxt(k,479)*y(k,164)
         mat(k,2564) = mat(k,2564) + .130_r8*rxt(k,479)*y(k,127)
         mat(k,478) = -(rxt(k,405)*y(k,250))
         mat(k,2046) = -rxt(k,405)*y(k,34)
         mat(k,925) = rxt(k,403)*y(k,235)
         mat(k,2244) = rxt(k,403)*y(k,227)
         mat(k,341) = -(rxt(k,406)*y(k,250) + rxt(k,409)*y(k,70))
         mat(k,2027) = -rxt(k,406)*y(k,35)
         mat(k,2339) = -rxt(k,409)*y(k,35)
         mat(k,302) = -(rxt(k,502)*y(k,250))
         mat(k,2023) = -rxt(k,502)*y(k,36)
         mat(k,729) = rxt(k,500)*y(k,235)
         mat(k,2233) = rxt(k,500)*y(k,228)
         mat(k,112) = -(rxt(k,231)*y(k,249))
         mat(k,1934) = -rxt(k,231)*y(k,37)
         mat(k,169) = -(rxt(k,232)*y(k,249))
         mat(k,1939) = -rxt(k,232)*y(k,38)
         mat(k,174) = -(rxt(k,262)*y(k,249))
         mat(k,1940) = -rxt(k,262)*y(k,39)
         mat(k,134) = -(rxt(k,233)*y(k,249))
         mat(k,1936) = -rxt(k,233)*y(k,40)
         mat(k,179) = -(rxt(k,234)*y(k,249))
         mat(k,1941) = -rxt(k,234)*y(k,41)
         mat(k,138) = -(rxt(k,235)*y(k,249))
         mat(k,1937) = -rxt(k,235)*y(k,42)
         mat(k,184) = -(rxt(k,236)*y(k,249))
         mat(k,1942) = -rxt(k,236)*y(k,43)
         mat(k,142) = -(rxt(k,237)*y(k,249))
         mat(k,1938) = -rxt(k,237)*y(k,44)
         mat(k,554) = -(rxt(k,316)*y(k,70) + rxt(k,317)*y(k,250) + rxt(k,343)*y(k,249))
         mat(k,2346) = -rxt(k,316)*y(k,45)
         mat(k,2056) = -rxt(k,317)*y(k,45)
         mat(k,1951) = -rxt(k,343)*y(k,45)
         mat(k,146) = -(rxt(k,318)*y(k,250))
         mat(k,1998) = -rxt(k,318)*y(k,46)
         mat(k,347) = -(rxt(k,319)*y(k,70) + rxt(k,320)*y(k,250))
         mat(k,2340) = -rxt(k,319)*y(k,47)
         mat(k,2028) = -rxt(k,320)*y(k,47)
         mat(k,1685) = -(rxt(k,204)*y(k,70) + rxt(k,243)*y(k,17) + rxt(k,348)*y(k,235) &
                      + rxt(k,349)*y(k,155) + rxt(k,350)*y(k,163) + rxt(k,351) &
                      *y(k,250))
         mat(k,2375) = -rxt(k,204)*y(k,51)
         mat(k,2519) = -rxt(k,243)*y(k,51)
         mat(k,2314) = -rxt(k,348)*y(k,51)
         mat(k,2655) = -rxt(k,349)*y(k,51)
         mat(k,2491) = -rxt(k,350)*y(k,51)
         mat(k,2142) = -rxt(k,351)*y(k,51)
         mat(k,744) = .400_r8*rxt(k,450)*y(k,250)
         mat(k,1017) = .340_r8*rxt(k,534)*y(k,164)
         mat(k,417) = .500_r8*rxt(k,421)*y(k,155)
         mat(k,655) = rxt(k,372)*y(k,164)
         mat(k,1134) = .500_r8*rxt(k,401)*y(k,164)
         mat(k,671) = .500_r8*rxt(k,389)*y(k,250)
         mat(k,871) = rxt(k,356)*y(k,250)
         mat(k,498) = .300_r8*rxt(k,357)*y(k,250)
         mat(k,1566) = (rxt(k,365)+rxt(k,366))*y(k,249)
         mat(k,1113) = rxt(k,332)*y(k,230)
         mat(k,2199) = rxt(k,213)*y(k,230)
         mat(k,1200) = .800_r8*rxt(k,394)*y(k,250)
         mat(k,902) = .910_r8*rxt(k,479)*y(k,164)
         mat(k,678) = .300_r8*rxt(k,470)*y(k,250)
         mat(k,1329) = .120_r8*rxt(k,432)*y(k,164)
         mat(k,687) = .500_r8*rxt(k,445)*y(k,250)
         mat(k,989) = .340_r8*rxt(k,537)*y(k,164)
         mat(k,1439) = .600_r8*rxt(k,446)*y(k,164)
         mat(k,1914) = .100_r8*rxt(k,452)*y(k,221) + rxt(k,355)*y(k,230) &
                      + .500_r8*rxt(k,423)*y(k,232) + .500_r8*rxt(k,391)*y(k,234) &
                      + .920_r8*rxt(k,462)*y(k,237) + .250_r8*rxt(k,430)*y(k,242) &
                      + rxt(k,439)*y(k,244) + rxt(k,413)*y(k,252) + rxt(k,417) &
                      *y(k,253) + .340_r8*rxt(k,546)*y(k,254) + .320_r8*rxt(k,551) &
                      *y(k,255) + .250_r8*rxt(k,487)*y(k,258)
         mat(k,2655) = mat(k,2655) + .500_r8*rxt(k,421)*y(k,16) + rxt(k,463)*y(k,237) &
                      + .250_r8*rxt(k,429)*y(k,242) + rxt(k,440)*y(k,244)
         mat(k,2585) = .340_r8*rxt(k,534)*y(k,6) + rxt(k,372)*y(k,28) &
                      + .500_r8*rxt(k,401)*y(k,33) + .910_r8*rxt(k,479)*y(k,127) &
                      + .120_r8*rxt(k,432)*y(k,134) + .340_r8*rxt(k,537)*y(k,139) &
                      + .600_r8*rxt(k,446)*y(k,140)
         mat(k,639) = rxt(k,396)*y(k,250)
         mat(k,1166) = .680_r8*rxt(k,555)*y(k,250)
         mat(k,1033) = .100_r8*rxt(k,452)*y(k,153)
         mat(k,883) = .700_r8*rxt(k,374)*y(k,230)
         mat(k,930) = rxt(k,402)*y(k,230)
         mat(k,1491) = rxt(k,385)*y(k,230) + rxt(k,459)*y(k,237) + .250_r8*rxt(k,426) &
                      *y(k,242) + rxt(k,435)*y(k,244) + .250_r8*rxt(k,484)*y(k,258)
         mat(k,1638) = rxt(k,332)*y(k,68) + rxt(k,213)*y(k,74) + rxt(k,355)*y(k,153) &
                      + .700_r8*rxt(k,374)*y(k,226) + rxt(k,402)*y(k,227) + rxt(k,385) &
                      *y(k,229) + (4.000_r8*rxt(k,352)+2.000_r8*rxt(k,353))*y(k,230) &
                      + 1.500_r8*rxt(k,460)*y(k,237) + .750_r8*rxt(k,465)*y(k,238) &
                      + .800_r8*rxt(k,474)*y(k,239) + .880_r8*rxt(k,427)*y(k,242) &
                      + 2.000_r8*rxt(k,436)*y(k,244) + .750_r8*rxt(k,539)*y(k,248) &
                      + .800_r8*rxt(k,415)*y(k,253) + .930_r8*rxt(k,544)*y(k,254) &
                      + .950_r8*rxt(k,549)*y(k,255) + .800_r8*rxt(k,485)*y(k,258)
         mat(k,663) = .500_r8*rxt(k,423)*y(k,153)
         mat(k,801) = .500_r8*rxt(k,391)*y(k,153)
         mat(k,2314) = mat(k,2314) + .450_r8*rxt(k,437)*y(k,244) + .150_r8*rxt(k,416) &
                      *y(k,253)
         mat(k,1362) = .920_r8*rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155) + rxt(k,459) &
                      *y(k,229) + 1.500_r8*rxt(k,460)*y(k,230)
         mat(k,1395) = .750_r8*rxt(k,465)*y(k,230)
         mat(k,1314) = .800_r8*rxt(k,474)*y(k,230)
         mat(k,1417) = .250_r8*rxt(k,430)*y(k,153) + .250_r8*rxt(k,429)*y(k,155) &
                      + .250_r8*rxt(k,426)*y(k,229) + .880_r8*rxt(k,427)*y(k,230)
         mat(k,1459) = rxt(k,439)*y(k,153) + rxt(k,440)*y(k,155) + rxt(k,435)*y(k,229) &
                      + 2.000_r8*rxt(k,436)*y(k,230) + .450_r8*rxt(k,437)*y(k,235) &
                      + 4.000_r8*rxt(k,438)*y(k,244)
         mat(k,1152) = .750_r8*rxt(k,539)*y(k,230)
         mat(k,1960) = (rxt(k,365)+rxt(k,366))*y(k,64)
         mat(k,2142) = mat(k,2142) + .400_r8*rxt(k,450)*y(k,1) + .500_r8*rxt(k,389) &
                      *y(k,60) + rxt(k,356)*y(k,62) + .300_r8*rxt(k,357)*y(k,63) &
                      + .800_r8*rxt(k,394)*y(k,90) + .300_r8*rxt(k,470)*y(k,128) &
                      + .500_r8*rxt(k,445)*y(k,138) + rxt(k,396)*y(k,170) &
                      + .680_r8*rxt(k,555)*y(k,210)
         mat(k,864) = rxt(k,413)*y(k,153)
         mat(k,1275) = rxt(k,417)*y(k,153) + .800_r8*rxt(k,415)*y(k,230) &
                      + .150_r8*rxt(k,416)*y(k,235)
         mat(k,1221) = .340_r8*rxt(k,546)*y(k,153) + .930_r8*rxt(k,544)*y(k,230)
         mat(k,1059) = .320_r8*rxt(k,551)*y(k,153) + .950_r8*rxt(k,549)*y(k,230)
         mat(k,1292) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,484)*y(k,229) &
                      + .800_r8*rxt(k,485)*y(k,230)
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
         mat(k,714) = -(rxt(k,321)*y(k,70) + rxt(k,322)*y(k,250) + rxt(k,344)*y(k,249))
         mat(k,2351) = -rxt(k,321)*y(k,52)
         mat(k,2076) = -rxt(k,322)*y(k,52)
         mat(k,1952) = -rxt(k,344)*y(k,52)
         mat(k,150) = -(rxt(k,323)*y(k,250))
         mat(k,1999) = -rxt(k,323)*y(k,53)
         mat(k,1187) = -(rxt(k,381)*y(k,155) + rxt(k,382)*y(k,250))
         mat(k,2631) = -rxt(k,381)*y(k,54)
         mat(k,2116) = -rxt(k,382)*y(k,54)
         mat(k,742) = .800_r8*rxt(k,450)*y(k,250)
         mat(k,416) = rxt(k,421)*y(k,155)
         mat(k,299) = rxt(k,377)*y(k,250)
         mat(k,388) = .500_r8*rxt(k,378)*y(k,250)
         mat(k,1127) = .500_r8*rxt(k,401)*y(k,164)
         mat(k,1429) = .100_r8*rxt(k,446)*y(k,164)
         mat(k,1893) = .400_r8*rxt(k,452)*y(k,221) + rxt(k,376)*y(k,226) &
                      + .270_r8*rxt(k,404)*y(k,227) + rxt(k,423)*y(k,232) + rxt(k,442) &
                      *y(k,246) + rxt(k,413)*y(k,252)
         mat(k,2631) = mat(k,2631) + rxt(k,421)*y(k,16)
         mat(k,2567) = .500_r8*rxt(k,401)*y(k,33) + .100_r8*rxt(k,446)*y(k,140)
         mat(k,1030) = .400_r8*rxt(k,452)*y(k,153)
         mat(k,881) = rxt(k,376)*y(k,153) + 3.200_r8*rxt(k,373)*y(k,226) &
                      + .800_r8*rxt(k,374)*y(k,230)
         mat(k,928) = .270_r8*rxt(k,404)*y(k,153)
         mat(k,1620) = .800_r8*rxt(k,374)*y(k,226)
         mat(k,661) = rxt(k,423)*y(k,153)
         mat(k,2290) = .200_r8*rxt(k,441)*y(k,246)
         mat(k,771) = rxt(k,442)*y(k,153) + .200_r8*rxt(k,441)*y(k,235)
         mat(k,2116) = mat(k,2116) + .800_r8*rxt(k,450)*y(k,1) + rxt(k,377)*y(k,30) &
                      + .500_r8*rxt(k,378)*y(k,31)
         mat(k,862) = rxt(k,413)*y(k,153)
         mat(k,440) = -(rxt(k,324)*y(k,70) + rxt(k,325)*y(k,250))
         mat(k,2344) = -rxt(k,324)*y(k,55)
         mat(k,2040) = -rxt(k,325)*y(k,55)
         mat(k,115) = -(rxt(k,383)*y(k,250))
         mat(k,1995) = -rxt(k,383)*y(k,56)
         mat(k,1073) = -(rxt(k,420)*y(k,250))
         mat(k,2107) = -rxt(k,420)*y(k,57)
         mat(k,741) = .800_r8*rxt(k,450)*y(k,250)
         mat(k,1010) = .520_r8*rxt(k,534)*y(k,164)
         mat(k,415) = .500_r8*rxt(k,421)*y(k,155)
         mat(k,982) = .520_r8*rxt(k,537)*y(k,164)
         mat(k,1888) = .250_r8*rxt(k,452)*y(k,221) + .820_r8*rxt(k,404)*y(k,227) &
                      + .500_r8*rxt(k,423)*y(k,232) + .270_r8*rxt(k,546)*y(k,254) &
                      + .040_r8*rxt(k,551)*y(k,255)
         mat(k,2623) = .500_r8*rxt(k,421)*y(k,16)
         mat(k,2562) = .520_r8*rxt(k,534)*y(k,6) + .520_r8*rxt(k,537)*y(k,139)
         mat(k,1160) = .500_r8*rxt(k,555)*y(k,250)
         mat(k,1029) = .250_r8*rxt(k,452)*y(k,153)
         mat(k,927) = .820_r8*rxt(k,404)*y(k,153) + .820_r8*rxt(k,402)*y(k,230)
         mat(k,1615) = .820_r8*rxt(k,402)*y(k,227) + .150_r8*rxt(k,544)*y(k,254) &
                      + .025_r8*rxt(k,549)*y(k,255)
         mat(k,660) = .500_r8*rxt(k,423)*y(k,153)
         mat(k,2107) = mat(k,2107) + .800_r8*rxt(k,450)*y(k,1) + .500_r8*rxt(k,555) &
                      *y(k,210)
         mat(k,1213) = .270_r8*rxt(k,546)*y(k,153) + .150_r8*rxt(k,544)*y(k,230)
         mat(k,1057) = .040_r8*rxt(k,551)*y(k,153) + .025_r8*rxt(k,549)*y(k,230)
         mat(k,1336) = -(rxt(k,407)*y(k,155) + rxt(k,408)*y(k,250))
         mat(k,2642) = -rxt(k,407)*y(k,58)
         mat(k,2127) = -rxt(k,408)*y(k,58)
         mat(k,1248) = rxt(k,410)*y(k,250)
         mat(k,1325) = .880_r8*rxt(k,432)*y(k,164)
         mat(k,1432) = .500_r8*rxt(k,446)*y(k,164)
         mat(k,1903) = .170_r8*rxt(k,505)*y(k,231) + .050_r8*rxt(k,468)*y(k,238) &
                      + .250_r8*rxt(k,430)*y(k,242) + .170_r8*rxt(k,511)*y(k,245) &
                      + .400_r8*rxt(k,521)*y(k,256) + .250_r8*rxt(k,487)*y(k,258) &
                      + .540_r8*rxt(k,527)*y(k,259) + .510_r8*rxt(k,530)*y(k,261)
         mat(k,2642) = mat(k,2642) + .050_r8*rxt(k,469)*y(k,238) + .250_r8*rxt(k,429) &
                      *y(k,242) + .250_r8*rxt(k,488)*y(k,258)
         mat(k,908) = rxt(k,411)*y(k,250)
         mat(k,2575) = .880_r8*rxt(k,432)*y(k,134) + .500_r8*rxt(k,446)*y(k,140)
         mat(k,1482) = .250_r8*rxt(k,426)*y(k,242) + .250_r8*rxt(k,484)*y(k,258)
         mat(k,1629) = .240_r8*rxt(k,427)*y(k,242) + .500_r8*rxt(k,415)*y(k,253) &
                      + .100_r8*rxt(k,485)*y(k,258)
         mat(k,854) = .170_r8*rxt(k,505)*y(k,153) + .070_r8*rxt(k,504)*y(k,235)
         mat(k,2300) = .070_r8*rxt(k,504)*y(k,231) + .070_r8*rxt(k,510)*y(k,245)
         mat(k,1388) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,1412) = .250_r8*rxt(k,430)*y(k,153) + .250_r8*rxt(k,429)*y(k,155) &
                      + .250_r8*rxt(k,426)*y(k,229) + .240_r8*rxt(k,427)*y(k,230)
         mat(k,961) = .170_r8*rxt(k,511)*y(k,153) + .070_r8*rxt(k,510)*y(k,235)
         mat(k,2127) = mat(k,2127) + rxt(k,410)*y(k,113) + rxt(k,411)*y(k,156)
         mat(k,1272) = .500_r8*rxt(k,415)*y(k,230)
         mat(k,822) = .400_r8*rxt(k,521)*y(k,153)
         mat(k,1289) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,229) + .100_r8*rxt(k,485)*y(k,230)
         mat(k,846) = .540_r8*rxt(k,527)*y(k,153)
         mat(k,588) = .510_r8*rxt(k,530)*y(k,153)
         mat(k,777) = -(rxt(k,388)*y(k,250))
         mat(k,2083) = -rxt(k,388)*y(k,59)
         mat(k,1121) = .120_r8*rxt(k,401)*y(k,164)
         mat(k,2551) = .120_r8*rxt(k,401)*y(k,33)
         mat(k,1472) = .100_r8*rxt(k,385)*y(k,230) + .150_r8*rxt(k,386)*y(k,235)
         mat(k,1606) = .100_r8*rxt(k,385)*y(k,229)
         mat(k,2267) = .150_r8*rxt(k,386)*y(k,229) + .150_r8*rxt(k,437)*y(k,244)
         mat(k,1451) = .150_r8*rxt(k,437)*y(k,235)
         mat(k,667) = -(rxt(k,389)*y(k,250))
         mat(k,2070) = -rxt(k,389)*y(k,60)
         mat(k,1471) = .400_r8*rxt(k,386)*y(k,235)
         mat(k,2259) = .400_r8*rxt(k,386)*y(k,229) + .400_r8*rxt(k,437)*y(k,244)
         mat(k,1449) = .400_r8*rxt(k,437)*y(k,235)
         mat(k,432) = -(rxt(k,326)*y(k,70) + rxt(k,327)*y(k,250))
         mat(k,2343) = -rxt(k,326)*y(k,61)
         mat(k,2039) = -rxt(k,327)*y(k,61)
         mat(k,870) = -(rxt(k,356)*y(k,250))
         mat(k,2092) = -rxt(k,356)*y(k,62)
         mat(k,879) = .300_r8*rxt(k,374)*y(k,230)
         mat(k,1607) = .300_r8*rxt(k,374)*y(k,226) + 2.000_r8*rxt(k,353)*y(k,230) &
                      + .250_r8*rxt(k,460)*y(k,237) + .250_r8*rxt(k,465)*y(k,238) &
                      + .200_r8*rxt(k,474)*y(k,239) + .250_r8*rxt(k,427)*y(k,242) &
                      + .250_r8*rxt(k,539)*y(k,248) + .500_r8*rxt(k,415)*y(k,253) &
                      + .250_r8*rxt(k,544)*y(k,254) + .250_r8*rxt(k,549)*y(k,255) &
                      + .300_r8*rxt(k,485)*y(k,258)
         mat(k,1346) = .250_r8*rxt(k,460)*y(k,230)
         mat(k,1377) = .250_r8*rxt(k,465)*y(k,230)
         mat(k,1302) = .200_r8*rxt(k,474)*y(k,230)
         mat(k,1406) = .250_r8*rxt(k,427)*y(k,230)
         mat(k,1145) = .250_r8*rxt(k,539)*y(k,230)
         mat(k,1269) = .500_r8*rxt(k,415)*y(k,230)
         mat(k,1211) = .250_r8*rxt(k,544)*y(k,230)
         mat(k,1054) = .250_r8*rxt(k,549)*y(k,230)
         mat(k,1282) = .300_r8*rxt(k,485)*y(k,230)
         mat(k,496) = -(rxt(k,357)*y(k,250))
         mat(k,2048) = -rxt(k,357)*y(k,63)
         mat(k,1604) = rxt(k,354)*y(k,235)
         mat(k,2246) = rxt(k,354)*y(k,230)
         mat(k,1564) = -(rxt(k,205)*y(k,70) + rxt(k,306)*y(k,89) + rxt(k,358)*y(k,250) &
                      + (rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,249))
         mat(k,2371) = -rxt(k,205)*y(k,64)
         mat(k,951) = -rxt(k,306)*y(k,64)
         mat(k,2138) = -rxt(k,358)*y(k,64)
         mat(k,1956) = -(rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,64)
         mat(k,1132) = .100_r8*rxt(k,401)*y(k,164)
         mat(k,2583) = .100_r8*rxt(k,401)*y(k,33)
         mat(k,118) = -(rxt(k,328)*y(k,250))
         mat(k,1996) = -rxt(k,328)*y(k,65)
         mat(k,448) = -(rxt(k,263)*y(k,249) + rxt(k,329)*y(k,70) + rxt(k,330)*y(k,250))
         mat(k,1950) = -rxt(k,263)*y(k,66)
         mat(k,2345) = -rxt(k,329)*y(k,66)
         mat(k,2041) = -rxt(k,330)*y(k,66)
         mat(k,122) = -(rxt(k,331)*y(k,250))
         mat(k,1997) = -rxt(k,331)*y(k,67)
         mat(k,1110) = -((rxt(k,332) + rxt(k,333)) * y(k,230) + (rxt(k,334) + rxt(k,335) &
                      ) * y(k,235) + rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155))
         mat(k,1616) = -(rxt(k,332) + rxt(k,333)) * y(k,68)
         mat(k,2287) = -(rxt(k,334) + rxt(k,335)) * y(k,68)
         mat(k,1889) = -rxt(k,336)*y(k,68)
         mat(k,2626) = -rxt(k,337)*y(k,68)
         mat(k,348) = rxt(k,319)*y(k,70) + rxt(k,320)*y(k,250)
         mat(k,2362) = rxt(k,319)*y(k,47)
         mat(k,2111) = rxt(k,320)*y(k,47)
         mat(k,398) = -(rxt(k,338)*y(k,70) + rxt(k,339)*y(k,250))
         mat(k,2342) = -rxt(k,338)*y(k,69)
         mat(k,2036) = -rxt(k,339)*y(k,69)
         mat(k,2386) = -(rxt(k,204)*y(k,51) + rxt(k,205)*y(k,64) + rxt(k,206)*y(k,93) &
                      + rxt(k,207)*y(k,95) + (rxt(k,208) + rxt(k,209)) * y(k,235) &
                      + rxt(k,210)*y(k,154) + rxt(k,212)*y(k,164) + rxt(k,219)*y(k,75) &
                      + rxt(k,228)*y(k,109) + rxt(k,254)*y(k,22) + rxt(k,312)*y(k,26) &
                      + rxt(k,314)*y(k,29) + rxt(k,316)*y(k,45) + rxt(k,319)*y(k,47) &
                      + rxt(k,321)*y(k,52) + rxt(k,324)*y(k,55) + rxt(k,326)*y(k,61) &
                      + rxt(k,329)*y(k,66) + rxt(k,379)*y(k,32) + rxt(k,409)*y(k,35) &
                      + (rxt(k,557) + rxt(k,558)) * y(k,83))
         mat(k,1696) = -rxt(k,204)*y(k,70)
         mat(k,1572) = -rxt(k,205)*y(k,70)
         mat(k,1526) = -rxt(k,206)*y(k,70)
         mat(k,696) = -rxt(k,207)*y(k,70)
         mat(k,2325) = -(rxt(k,208) + rxt(k,209)) * y(k,70)
         mat(k,1825) = -rxt(k,210)*y(k,70)
         mat(k,2596) = -rxt(k,212)*y(k,70)
         mat(k,1087) = -rxt(k,219)*y(k,70)
         mat(k,1744) = -rxt(k,228)*y(k,70)
         mat(k,921) = -rxt(k,254)*y(k,70)
         mat(k,219) = -rxt(k,312)*y(k,70)
         mat(k,291) = -rxt(k,314)*y(k,70)
         mat(k,558) = -rxt(k,316)*y(k,70)
         mat(k,351) = -rxt(k,319)*y(k,70)
         mat(k,720) = -rxt(k,321)*y(k,70)
         mat(k,446) = -rxt(k,324)*y(k,70)
         mat(k,437) = -rxt(k,326)*y(k,70)
         mat(k,452) = -rxt(k,329)*y(k,70)
         mat(k,339) = -rxt(k,379)*y(k,70)
         mat(k,345) = -rxt(k,409)*y(k,70)
         mat(k,1097) = -(rxt(k,557) + rxt(k,558)) * y(k,70)
         mat(k,2181) = rxt(k,249)*y(k,74)
         mat(k,219) = mat(k,219) + 5.000_r8*rxt(k,312)*y(k,70) + 3.060_r8*rxt(k,313) &
                      *y(k,250)
         mat(k,291) = mat(k,291) + 2.000_r8*rxt(k,314)*y(k,70) + 2.000_r8*rxt(k,315) &
                      *y(k,250)
         mat(k,114) = 4.000_r8*rxt(k,231)*y(k,249)
         mat(k,172) = rxt(k,232)*y(k,249)
         mat(k,137) = 2.000_r8*rxt(k,233)*y(k,249)
         mat(k,183) = 2.000_r8*rxt(k,234)*y(k,249)
         mat(k,141) = 2.000_r8*rxt(k,235)*y(k,249)
         mat(k,188) = rxt(k,236)*y(k,249)
         mat(k,145) = 2.000_r8*rxt(k,237)*y(k,249)
         mat(k,148) = rxt(k,318)*y(k,250)
         mat(k,152) = 3.000_r8*rxt(k,323)*y(k,250)
         mat(k,446) = mat(k,446) + rxt(k,325)*y(k,250)
         mat(k,120) = rxt(k,328)*y(k,250)
         mat(k,124) = 2.000_r8*rxt(k,331)*y(k,250)
         mat(k,1118) = 2.000_r8*rxt(k,336)*y(k,153) + 2.000_r8*rxt(k,337)*y(k,155) &
                      + 2.000_r8*rxt(k,332)*y(k,230) + rxt(k,335)*y(k,235)
         mat(k,403) = rxt(k,339)*y(k,250)
         mat(k,2386) = mat(k,2386) + 5.000_r8*rxt(k,312)*y(k,26) + 2.000_r8*rxt(k,314) &
                      *y(k,29)
         mat(k,2210) = rxt(k,249)*y(k,21) + (4.000_r8*rxt(k,214)+2.000_r8*rxt(k,216)) &
                      *y(k,74) + rxt(k,286)*y(k,125) + rxt(k,218)*y(k,153) &
                      + rxt(k,223)*y(k,163) + rxt(k,568)*y(k,181) + rxt(k,213) &
                      *y(k,230) + rxt(k,224)*y(k,250)
         mat(k,265) = rxt(k,311)*y(k,249)
         mat(k,261) = rxt(k,345)*y(k,249) + rxt(k,340)*y(k,250)
         mat(k,270) = rxt(k,346)*y(k,249) + rxt(k,341)*y(k,250)
         mat(k,321) = rxt(k,347)*y(k,249) + rxt(k,342)*y(k,250)
         mat(k,1720) = rxt(k,226)*y(k,163) + rxt(k,238)*y(k,249) + rxt(k,227)*y(k,250)
         mat(k,2416) = rxt(k,286)*y(k,74)
         mat(k,1925) = 2.000_r8*rxt(k,336)*y(k,68) + rxt(k,218)*y(k,74)
         mat(k,2666) = 2.000_r8*rxt(k,337)*y(k,68)
         mat(k,2502) = rxt(k,223)*y(k,74) + rxt(k,226)*y(k,101)
         mat(k,1553) = rxt(k,568)*y(k,74)
         mat(k,1648) = 2.000_r8*rxt(k,332)*y(k,68) + rxt(k,213)*y(k,74)
         mat(k,2325) = mat(k,2325) + rxt(k,335)*y(k,68)
         mat(k,1971) = 4.000_r8*rxt(k,231)*y(k,37) + rxt(k,232)*y(k,38) &
                      + 2.000_r8*rxt(k,233)*y(k,40) + 2.000_r8*rxt(k,234)*y(k,41) &
                      + 2.000_r8*rxt(k,235)*y(k,42) + rxt(k,236)*y(k,43) &
                      + 2.000_r8*rxt(k,237)*y(k,44) + rxt(k,311)*y(k,81) + rxt(k,345) &
                      *y(k,98) + rxt(k,346)*y(k,99) + rxt(k,347)*y(k,100) + rxt(k,238) &
                      *y(k,101)
         mat(k,2153) = 3.060_r8*rxt(k,313)*y(k,26) + 2.000_r8*rxt(k,315)*y(k,29) &
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
         mat(k,2335) = rxt(k,219)*y(k,75)
         mat(k,2191) = 2.000_r8*rxt(k,215)*y(k,74)
         mat(k,1079) = rxt(k,219)*y(k,70) + (rxt(k,666)+rxt(k,675)+rxt(k,684)) &
                      *y(k,101)
         mat(k,1705) = (rxt(k,666)+rxt(k,675)+rxt(k,684))*y(k,75) + (rxt(k,585) &
                       +rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,109)
         mat(k,1728) = (rxt(k,585)+rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,101)
         mat(k,2190) = 2.000_r8*rxt(k,240)*y(k,74)
         mat(k,622) = -(rxt(k,211)*y(k,250))
         mat(k,2064) = -rxt(k,211)*y(k,73)
         mat(k,2348) = rxt(k,210)*y(k,154)
         mat(k,1708) = rxt(k,601)*y(k,143)
         mat(k,406) = rxt(k,601)*y(k,101)
         mat(k,1790) = rxt(k,210)*y(k,70)
         mat(k,2208) = -(rxt(k,213)*y(k,230) + (4._r8*rxt(k,214) + 4._r8*rxt(k,215) &
                      + 4._r8*rxt(k,216) + 4._r8*rxt(k,240)) * y(k,74) + rxt(k,217) &
                      *y(k,235) + rxt(k,218)*y(k,153) + rxt(k,220)*y(k,154) + rxt(k,223) &
                      *y(k,163) + (rxt(k,224) + rxt(k,225)) * y(k,250) + (rxt(k,248) &
                      + rxt(k,249) + rxt(k,250)) * y(k,21) + (rxt(k,285) + rxt(k,286) &
                      + rxt(k,287)) * y(k,125) + rxt(k,568)*y(k,181))
         mat(k,1646) = -rxt(k,213)*y(k,74)
         mat(k,2323) = -rxt(k,217)*y(k,74)
         mat(k,1923) = -rxt(k,218)*y(k,74)
         mat(k,1823) = -rxt(k,220)*y(k,74)
         mat(k,2500) = -rxt(k,223)*y(k,74)
         mat(k,2151) = -(rxt(k,224) + rxt(k,225)) * y(k,74)
         mat(k,2179) = -(rxt(k,248) + rxt(k,249) + rxt(k,250)) * y(k,74)
         mat(k,2414) = -(rxt(k,285) + rxt(k,286) + rxt(k,287)) * y(k,74)
         mat(k,1551) = -rxt(k,568)*y(k,74)
         mat(k,2384) = rxt(k,228)*y(k,109) + rxt(k,212)*y(k,164) + rxt(k,209)*y(k,235)
         mat(k,1086) = rxt(k,221)*y(k,163)
         mat(k,1719) = rxt(k,239)*y(k,249)
         mat(k,1743) = rxt(k,228)*y(k,70) + rxt(k,229)*y(k,163) + rxt(k,230)*y(k,250)
         mat(k,2500) = mat(k,2500) + rxt(k,221)*y(k,75) + rxt(k,229)*y(k,109)
         mat(k,2594) = rxt(k,212)*y(k,70)
         mat(k,538) = rxt(k,573)*y(k,181)
         mat(k,1551) = mat(k,1551) + rxt(k,573)*y(k,166)
         mat(k,2323) = mat(k,2323) + rxt(k,209)*y(k,70)
         mat(k,1969) = rxt(k,239)*y(k,101)
         mat(k,2151) = mat(k,2151) + rxt(k,230)*y(k,109)
         mat(k,1080) = -(rxt(k,219)*y(k,70) + rxt(k,221)*y(k,163) + rxt(k,222) &
                      *y(k,250) + (rxt(k,666) + rxt(k,675) + rxt(k,684)) * y(k,101))
         mat(k,2359) = -rxt(k,219)*y(k,75)
         mat(k,2480) = -rxt(k,221)*y(k,75)
         mat(k,2108) = -rxt(k,222)*y(k,75)
         mat(k,1709) = -(rxt(k,666) + rxt(k,675) + rxt(k,684)) * y(k,75)
         mat(k,2195) = rxt(k,220)*y(k,154)
         mat(k,1801) = rxt(k,220)*y(k,74)
         mat(k,1205) = -(rxt(k,368)*y(k,250))
         mat(k,2118) = -rxt(k,368)*y(k,77)
         mat(k,1013) = .230_r8*rxt(k,534)*y(k,164)
         mat(k,2515) = rxt(k,243)*y(k,51)
         mat(k,329) = .350_r8*rxt(k,370)*y(k,250)
         mat(k,654) = .630_r8*rxt(k,372)*y(k,164)
         mat(k,1128) = .560_r8*rxt(k,401)*y(k,164)
         mat(k,1678) = rxt(k,243)*y(k,17) + rxt(k,204)*y(k,70) + rxt(k,349)*y(k,155) &
                      + rxt(k,350)*y(k,163) + rxt(k,351)*y(k,250)
         mat(k,441) = rxt(k,324)*y(k,70)
         mat(k,1335) = rxt(k,407)*y(k,155) + rxt(k,408)*y(k,250)
         mat(k,1111) = rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155) + (rxt(k,332) &
                       +rxt(k,333))*y(k,230) + rxt(k,335)*y(k,235)
         mat(k,2365) = rxt(k,204)*y(k,51) + rxt(k,324)*y(k,55)
         mat(k,1048) = rxt(k,395)*y(k,250)
         mat(k,894) = .620_r8*rxt(k,479)*y(k,164)
         mat(k,1323) = .650_r8*rxt(k,432)*y(k,164)
         mat(k,985) = .230_r8*rxt(k,537)*y(k,164)
         mat(k,1430) = .560_r8*rxt(k,446)*y(k,164)
         mat(k,1895) = rxt(k,336)*y(k,68) + .170_r8*rxt(k,505)*y(k,231) &
                      + .220_r8*rxt(k,430)*y(k,242) + .400_r8*rxt(k,508)*y(k,243) &
                      + .350_r8*rxt(k,511)*y(k,245) + .225_r8*rxt(k,546)*y(k,254) &
                      + .250_r8*rxt(k,487)*y(k,258)
         mat(k,2633) = rxt(k,349)*y(k,51) + rxt(k,407)*y(k,58) + rxt(k,337)*y(k,68) &
                      + .220_r8*rxt(k,429)*y(k,242) + .500_r8*rxt(k,488)*y(k,258)
         mat(k,2482) = rxt(k,350)*y(k,51) + rxt(k,562)*y(k,167)
         mat(k,2569) = .230_r8*rxt(k,534)*y(k,6) + .630_r8*rxt(k,372)*y(k,28) &
                      + .560_r8*rxt(k,401)*y(k,33) + .620_r8*rxt(k,479)*y(k,127) &
                      + .650_r8*rxt(k,432)*y(k,134) + .230_r8*rxt(k,537)*y(k,139) &
                      + .560_r8*rxt(k,446)*y(k,140)
         mat(k,427) = rxt(k,562)*y(k,163) + rxt(k,563)*y(k,250)
         mat(k,1162) = .700_r8*rxt(k,555)*y(k,250)
         mat(k,1477) = .220_r8*rxt(k,426)*y(k,242) + .250_r8*rxt(k,484)*y(k,258)
         mat(k,1622) = (rxt(k,332)+rxt(k,333))*y(k,68) + .110_r8*rxt(k,427)*y(k,242) &
                      + .125_r8*rxt(k,544)*y(k,254) + .200_r8*rxt(k,485)*y(k,258)
         mat(k,853) = .170_r8*rxt(k,505)*y(k,153) + .070_r8*rxt(k,504)*y(k,235)
         mat(k,2292) = rxt(k,335)*y(k,68) + .070_r8*rxt(k,504)*y(k,231) &
                      + .160_r8*rxt(k,507)*y(k,243) + .140_r8*rxt(k,510)*y(k,245)
         mat(k,1408) = .220_r8*rxt(k,430)*y(k,153) + .220_r8*rxt(k,429)*y(k,155) &
                      + .220_r8*rxt(k,426)*y(k,229) + .110_r8*rxt(k,427)*y(k,230)
         mat(k,808) = .400_r8*rxt(k,508)*y(k,153) + .160_r8*rxt(k,507)*y(k,235)
         mat(k,960) = .350_r8*rxt(k,511)*y(k,153) + .140_r8*rxt(k,510)*y(k,235)
         mat(k,2118) = mat(k,2118) + .350_r8*rxt(k,370)*y(k,27) + rxt(k,351)*y(k,51) &
                      + rxt(k,408)*y(k,58) + rxt(k,395)*y(k,91) + rxt(k,563)*y(k,167) &
                      + .700_r8*rxt(k,555)*y(k,210)
         mat(k,1216) = .225_r8*rxt(k,546)*y(k,153) + .125_r8*rxt(k,544)*y(k,230)
         mat(k,1286) = .250_r8*rxt(k,487)*y(k,153) + .500_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,229) + .200_r8*rxt(k,485)*y(k,230)
         mat(k,1003) = .270_r8*rxt(k,534)*y(k,164)
         mat(k,1123) = .200_r8*rxt(k,401)*y(k,164)
         mat(k,778) = rxt(k,388)*y(k,250)
         mat(k,668) = .500_r8*rxt(k,389)*y(k,250)
         mat(k,1204) = rxt(k,368)*y(k,250)
         mat(k,1196) = .800_r8*rxt(k,394)*y(k,250)
         mat(k,1046) = rxt(k,395)*y(k,250)
         mat(k,1038) = rxt(k,360)*y(k,250)
         mat(k,684) = .500_r8*rxt(k,445)*y(k,250)
         mat(k,975) = .270_r8*rxt(k,537)*y(k,164)
         mat(k,1426) = .100_r8*rxt(k,446)*y(k,164)
         mat(k,1878) = rxt(k,387)*y(k,229) + .900_r8*rxt(k,546)*y(k,254)
         mat(k,2553) = .270_r8*rxt(k,534)*y(k,6) + .200_r8*rxt(k,401)*y(k,33) &
                      + .270_r8*rxt(k,537)*y(k,139) + .100_r8*rxt(k,446)*y(k,140)
         mat(k,1159) = 1.800_r8*rxt(k,555)*y(k,250)
         mat(k,1473) = rxt(k,387)*y(k,153) + 4.000_r8*rxt(k,384)*y(k,229) &
                      + .900_r8*rxt(k,385)*y(k,230) + rxt(k,459)*y(k,237) &
                      + 2.000_r8*rxt(k,435)*y(k,244) + rxt(k,484)*y(k,258)
         mat(k,1608) = .900_r8*rxt(k,385)*y(k,229) + rxt(k,436)*y(k,244) &
                      + .500_r8*rxt(k,544)*y(k,254)
         mat(k,2276) = .450_r8*rxt(k,437)*y(k,244)
         mat(k,1347) = rxt(k,459)*y(k,229)
         mat(k,1452) = 2.000_r8*rxt(k,435)*y(k,229) + rxt(k,436)*y(k,230) &
                      + .450_r8*rxt(k,437)*y(k,235) + 4.000_r8*rxt(k,438)*y(k,244)
         mat(k,2093) = rxt(k,388)*y(k,59) + .500_r8*rxt(k,389)*y(k,60) + rxt(k,368) &
                      *y(k,77) + .800_r8*rxt(k,394)*y(k,90) + rxt(k,395)*y(k,91) &
                      + rxt(k,360)*y(k,103) + .500_r8*rxt(k,445)*y(k,138) &
                      + 1.800_r8*rxt(k,555)*y(k,210)
         mat(k,1212) = .900_r8*rxt(k,546)*y(k,153) + .500_r8*rxt(k,544)*y(k,230)
         mat(k,1283) = rxt(k,484)*y(k,229)
         mat(k,217) = .470_r8*rxt(k,313)*y(k,250)
         mat(k,1109) = rxt(k,333)*y(k,230) + rxt(k,334)*y(k,235)
         mat(k,397) = rxt(k,338)*y(k,70) + rxt(k,339)*y(k,250)
         mat(k,2341) = rxt(k,338)*y(k,69)
         mat(k,1602) = rxt(k,333)*y(k,68)
         mat(k,2240) = rxt(k,334)*y(k,68)
         mat(k,2035) = .470_r8*rxt(k,313)*y(k,26) + rxt(k,339)*y(k,69)
         mat(k,271) = -(rxt(k,310)*y(k,249))
         mat(k,1948) = -rxt(k,310)*y(k,80)
         mat(k,170) = rxt(k,232)*y(k,249)
         mat(k,175) = rxt(k,262)*y(k,249)
         mat(k,181) = rxt(k,234)*y(k,249)
         mat(k,139) = 2.000_r8*rxt(k,235)*y(k,249)
         mat(k,185) = 2.000_r8*rxt(k,236)*y(k,249)
         mat(k,143) = rxt(k,237)*y(k,249)
         mat(k,127) = 2.000_r8*rxt(k,264)*y(k,249)
         mat(k,267) = rxt(k,346)*y(k,249) + rxt(k,341)*y(k,250)
         mat(k,316) = rxt(k,347)*y(k,249) + rxt(k,342)*y(k,250)
         mat(k,1948) = mat(k,1948) + rxt(k,232)*y(k,38) + rxt(k,262)*y(k,39) &
                      + rxt(k,234)*y(k,41) + 2.000_r8*rxt(k,235)*y(k,42) &
                      + 2.000_r8*rxt(k,236)*y(k,43) + rxt(k,237)*y(k,44) &
                      + 2.000_r8*rxt(k,264)*y(k,94) + rxt(k,346)*y(k,99) + rxt(k,347) &
                      *y(k,100)
         mat(k,2015) = rxt(k,341)*y(k,99) + rxt(k,342)*y(k,100)
         mat(k,262) = -(rxt(k,311)*y(k,249))
         mat(k,1946) = -rxt(k,311)*y(k,81)
         mat(k,135) = rxt(k,233)*y(k,249)
         mat(k,180) = rxt(k,234)*y(k,249)
         mat(k,258) = rxt(k,345)*y(k,249) + rxt(k,340)*y(k,250)
         mat(k,1946) = mat(k,1946) + rxt(k,233)*y(k,40) + rxt(k,234)*y(k,41) &
                      + rxt(k,345)*y(k,98)
         mat(k,2013) = rxt(k,340)*y(k,98)
         mat(k,228) = -(rxt(k,503)*y(k,250))
         mat(k,2007) = -rxt(k,503)*y(k,82)
         mat(k,222) = .180_r8*rxt(k,523)*y(k,250)
         mat(k,2007) = mat(k,2007) + .180_r8*rxt(k,523)*y(k,212)
         mat(k,1090) = -(rxt(k,556)*y(k,21) + (rxt(k,557) + rxt(k,558)) * y(k,70) &
                      + rxt(k,559)*y(k,125) + rxt(k,560)*y(k,155) + (rxt(k,561) &
                      + rxt(k,575)) * y(k,250))
         mat(k,2166) = -rxt(k,556)*y(k,83)
         mat(k,2360) = -(rxt(k,557) + rxt(k,558)) * y(k,83)
         mat(k,2401) = -rxt(k,559)*y(k,83)
         mat(k,2624) = -rxt(k,560)*y(k,83)
         mat(k,2109) = -(rxt(k,561) + rxt(k,575)) * y(k,83)
         mat(k,797) = rxt(k,390)*y(k,235)
         mat(k,2231) = rxt(k,390)*y(k,234)
         mat(k,949) = -(rxt(k,306)*y(k,64) + rxt(k,307)*y(k,93) + rxt(k,308)*y(k,262) &
                      + rxt(k,309)*y(k,106))
         mat(k,1561) = -rxt(k,306)*y(k,89)
         mat(k,1518) = -rxt(k,307)*y(k,89)
         mat(k,2678) = -rxt(k,308)*y(k,89)
         mat(k,1753) = -rxt(k,309)*y(k,89)
         mat(k,176) = rxt(k,262)*y(k,249)
         mat(k,186) = rxt(k,236)*y(k,249)
         mat(k,272) = 2.000_r8*rxt(k,310)*y(k,249)
         mat(k,263) = rxt(k,311)*y(k,249)
         mat(k,1953) = rxt(k,262)*y(k,39) + rxt(k,236)*y(k,43) + 2.000_r8*rxt(k,310) &
                      *y(k,80) + rxt(k,311)*y(k,81)
         mat(k,1198) = -(rxt(k,394)*y(k,250))
         mat(k,2117) = -rxt(k,394)*y(k,90)
         mat(k,675) = .700_r8*rxt(k,470)*y(k,250)
         mat(k,629) = .500_r8*rxt(k,471)*y(k,250)
         mat(k,456) = rxt(k,482)*y(k,250)
         mat(k,1894) = .050_r8*rxt(k,468)*y(k,238) + .530_r8*rxt(k,430)*y(k,242) &
                      + .225_r8*rxt(k,546)*y(k,254) + .250_r8*rxt(k,487)*y(k,258)
         mat(k,2632) = .050_r8*rxt(k,469)*y(k,238) + .530_r8*rxt(k,429)*y(k,242) &
                      + .250_r8*rxt(k,488)*y(k,258)
         mat(k,1476) = .530_r8*rxt(k,426)*y(k,242) + .250_r8*rxt(k,484)*y(k,258)
         mat(k,1621) = .260_r8*rxt(k,427)*y(k,242) + .125_r8*rxt(k,544)*y(k,254) &
                      + .100_r8*rxt(k,485)*y(k,258)
         mat(k,1381) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,1407) = .530_r8*rxt(k,430)*y(k,153) + .530_r8*rxt(k,429)*y(k,155) &
                      + .530_r8*rxt(k,426)*y(k,229) + .260_r8*rxt(k,427)*y(k,230)
         mat(k,2117) = mat(k,2117) + .700_r8*rxt(k,470)*y(k,128) + .500_r8*rxt(k,471) &
                      *y(k,129) + rxt(k,482)*y(k,144)
         mat(k,1215) = .225_r8*rxt(k,546)*y(k,153) + .125_r8*rxt(k,544)*y(k,230)
         mat(k,1285) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,229) + .100_r8*rxt(k,485)*y(k,230)
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
         mat(k,1047) = -(rxt(k,395)*y(k,250))
         mat(k,2104) = -rxt(k,395)*y(k,91)
         mat(k,328) = .650_r8*rxt(k,370)*y(k,250)
         mat(k,1197) = .200_r8*rxt(k,394)*y(k,250)
         mat(k,1174) = rxt(k,483)*y(k,250)
         mat(k,1885) = rxt(k,494)*y(k,223) + .050_r8*rxt(k,468)*y(k,238) &
                      + .400_r8*rxt(k,508)*y(k,243) + .170_r8*rxt(k,511)*y(k,245) &
                      + .700_r8*rxt(k,514)*y(k,251) + .600_r8*rxt(k,521)*y(k,256) &
                      + .250_r8*rxt(k,487)*y(k,258) + .340_r8*rxt(k,527)*y(k,259) &
                      + .170_r8*rxt(k,530)*y(k,261)
         mat(k,2620) = .050_r8*rxt(k,469)*y(k,238) + .250_r8*rxt(k,488)*y(k,258)
         mat(k,580) = rxt(k,494)*y(k,153)
         mat(k,1474) = .250_r8*rxt(k,484)*y(k,258)
         mat(k,1612) = .100_r8*rxt(k,485)*y(k,258)
         mat(k,2283) = .160_r8*rxt(k,507)*y(k,243) + .070_r8*rxt(k,510)*y(k,245)
         mat(k,1380) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,807) = .400_r8*rxt(k,508)*y(k,153) + .160_r8*rxt(k,507)*y(k,235)
         mat(k,959) = .170_r8*rxt(k,511)*y(k,153) + .070_r8*rxt(k,510)*y(k,235)
         mat(k,2104) = mat(k,2104) + .650_r8*rxt(k,370)*y(k,27) + .200_r8*rxt(k,394) &
                      *y(k,90) + rxt(k,483)*y(k,145)
         mat(k,530) = .700_r8*rxt(k,514)*y(k,153)
         mat(k,820) = .600_r8*rxt(k,521)*y(k,153)
         mat(k,1284) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,229) + .100_r8*rxt(k,485)*y(k,230)
         mat(k,844) = .340_r8*rxt(k,527)*y(k,153)
         mat(k,587) = .170_r8*rxt(k,530)*y(k,153)
         mat(k,2470) = -((rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,235) + rxt(k,170) &
                      *y(k,164))
         mat(k,2328) = -(rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,92)
         mat(k,2599) = -rxt(k,170)*y(k,92)
         mat(k,1699) = rxt(k,351)*y(k,250)
         mat(k,1573) = rxt(k,365)*y(k,249)
         mat(k,2389) = rxt(k,206)*y(k,93)
         mat(k,955) = rxt(k,307)*y(k,93)
         mat(k,1527) = rxt(k,206)*y(k,70) + rxt(k,307)*y(k,89) + rxt(k,162)*y(k,163) &
                      + rxt(k,153)*y(k,249) + rxt(k,171)*y(k,250)
         mat(k,1513) = rxt(k,266)*y(k,249)
         mat(k,1723) = rxt(k,239)*y(k,249)
         mat(k,574) = rxt(k,192)*y(k,250)
         mat(k,2505) = rxt(k,162)*y(k,93) + rxt(k,174)*y(k,250)
         mat(k,430) = rxt(k,563)*y(k,250)
         mat(k,608) = rxt(k,569)*y(k,250)
         mat(k,1554) = rxt(k,574)*y(k,250)
         mat(k,1974) = rxt(k,365)*y(k,64) + rxt(k,153)*y(k,93) + rxt(k,266)*y(k,97) &
                      + rxt(k,239)*y(k,101)
         mat(k,2156) = rxt(k,351)*y(k,51) + rxt(k,171)*y(k,93) + rxt(k,192)*y(k,141) &
                      + rxt(k,174)*y(k,163) + rxt(k,563)*y(k,167) + rxt(k,569) &
                      *y(k,179) + rxt(k,574)*y(k,181)
         mat(k,1519) = -(rxt(k,153)*y(k,249) + rxt(k,162)*y(k,163) + rxt(k,171) &
                      *y(k,250) + rxt(k,206)*y(k,70) + rxt(k,307)*y(k,89))
         mat(k,1955) = -rxt(k,153)*y(k,93)
         mat(k,2484) = -rxt(k,162)*y(k,93)
         mat(k,2135) = -rxt(k,171)*y(k,93)
         mat(k,2369) = -rxt(k,206)*y(k,93)
         mat(k,950) = -rxt(k,307)*y(k,93)
         mat(k,1563) = rxt(k,366)*y(k,249)
         mat(k,2453) = rxt(k,164)*y(k,235)
         mat(k,2308) = rxt(k,164)*y(k,92)
         mat(k,1955) = mat(k,1955) + rxt(k,366)*y(k,64)
         mat(k,126) = -(rxt(k,264)*y(k,249))
         mat(k,1935) = -rxt(k,264)*y(k,94)
         mat(k,692) = -(rxt(k,163)*y(k,163) + rxt(k,172)*y(k,250) + rxt(k,207)*y(k,70))
         mat(k,2478) = -rxt(k,163)*y(k,95)
         mat(k,2073) = -rxt(k,172)*y(k,95)
         mat(k,2350) = -rxt(k,207)*y(k,95)
         mat(k,2260) = 2.000_r8*rxt(k,178)*y(k,235)
         mat(k,2073) = mat(k,2073) + 2.000_r8*rxt(k,177)*y(k,250)
         mat(k,293) = rxt(k,576)*y(k,262)
         mat(k,2675) = rxt(k,576)*y(k,183)
         mat(k,1504) = -(rxt(k,259)*y(k,163) + rxt(k,260)*y(k,250) + (rxt(k,265) &
                      + rxt(k,266)) * y(k,249) + (rxt(k,583) + rxt(k,660) + rxt(k,668) &
                      + rxt(k,677)) * y(k,109) + (rxt(k,584) + rxt(k,658) + rxt(k,671) &
                      + rxt(k,680)) * y(k,108) + (rxt(k,591) + rxt(k,687) + rxt(k,691) &
                      + rxt(k,695)) * y(k,110))
         mat(k,2483) = -rxt(k,259)*y(k,97)
         mat(k,2134) = -rxt(k,260)*y(k,97)
         mat(k,1954) = -(rxt(k,265) + rxt(k,266)) * y(k,97)
         mat(k,1732) = -(rxt(k,583) + rxt(k,660) + rxt(k,668) + rxt(k,677)) * y(k,97)
         mat(k,1581) = -(rxt(k,584) + rxt(k,658) + rxt(k,671) + rxt(k,680)) * y(k,97)
         mat(k,1657) = -(rxt(k,591) + rxt(k,687) + rxt(k,691) + rxt(k,695)) * y(k,97)
         mat(k,2516) = rxt(k,243)*y(k,51) + rxt(k,244)*y(k,235)
         mat(k,1679) = rxt(k,243)*y(k,17)
         mat(k,2307) = rxt(k,244)*y(k,17)
         mat(k,257) = -(rxt(k,340)*y(k,250) + rxt(k,345)*y(k,249))
         mat(k,2012) = -rxt(k,340)*y(k,98)
         mat(k,1945) = -rxt(k,345)*y(k,98)
         mat(k,266) = -(rxt(k,341)*y(k,250) + rxt(k,346)*y(k,249))
         mat(k,2014) = -rxt(k,341)*y(k,99)
         mat(k,1947) = -rxt(k,346)*y(k,99)
         mat(k,317) = -(rxt(k,342)*y(k,250) + rxt(k,347)*y(k,249))
         mat(k,2024) = -rxt(k,342)*y(k,100)
         mat(k,1949) = -rxt(k,347)*y(k,100)
         mat(k,1712) = -(rxt(k,226)*y(k,163) + rxt(k,227)*y(k,250) + (rxt(k,238) &
                      + rxt(k,239)) * y(k,249) + (rxt(k,585) + rxt(k,656) + rxt(k,667) &
                      + rxt(k,676)) * y(k,109) + (rxt(k,586) + rxt(k,657) + rxt(k,670) &
                      + rxt(k,679)) * y(k,108) + (rxt(k,590) + rxt(k,686) + rxt(k,690) &
                      + rxt(k,694)) * y(k,110) + rxt(k,601)*y(k,143) + (rxt(k,666) &
                      + rxt(k,675) + rxt(k,684)) * y(k,75))
         mat(k,2492) = -rxt(k,226)*y(k,101)
         mat(k,2143) = -rxt(k,227)*y(k,101)
         mat(k,1961) = -(rxt(k,238) + rxt(k,239)) * y(k,101)
         mat(k,1736) = -(rxt(k,585) + rxt(k,656) + rxt(k,667) + rxt(k,676)) * y(k,101)
         mat(k,1585) = -(rxt(k,586) + rxt(k,657) + rxt(k,670) + rxt(k,679)) * y(k,101)
         mat(k,1661) = -(rxt(k,590) + rxt(k,686) + rxt(k,690) + rxt(k,694)) * y(k,101)
         mat(k,407) = -rxt(k,601)*y(k,101)
         mat(k,1081) = -(rxt(k,666) + rxt(k,675) + rxt(k,684)) * y(k,101)
         mat(k,289) = rxt(k,314)*y(k,70)
         mat(k,337) = rxt(k,379)*y(k,70)
         mat(k,343) = rxt(k,409)*y(k,70)
         mat(k,555) = rxt(k,316)*y(k,70)
         mat(k,349) = rxt(k,319)*y(k,70)
         mat(k,1686) = rxt(k,204)*y(k,70)
         mat(k,716) = rxt(k,321)*y(k,70)
         mat(k,443) = 2.000_r8*rxt(k,324)*y(k,70)
         mat(k,434) = rxt(k,326)*y(k,70)
         mat(k,1567) = rxt(k,205)*y(k,70)
         mat(k,449) = rxt(k,329)*y(k,70)
         mat(k,401) = rxt(k,338)*y(k,70)
         mat(k,2376) = rxt(k,314)*y(k,29) + rxt(k,379)*y(k,32) + rxt(k,409)*y(k,35) &
                      + rxt(k,316)*y(k,45) + rxt(k,319)*y(k,47) + rxt(k,204)*y(k,51) &
                      + rxt(k,321)*y(k,52) + 2.000_r8*rxt(k,324)*y(k,55) + rxt(k,326) &
                      *y(k,61) + rxt(k,205)*y(k,64) + rxt(k,329)*y(k,66) + rxt(k,338) &
                      *y(k,69) + rxt(k,558)*y(k,83) + rxt(k,206)*y(k,93) + rxt(k,207) &
                      *y(k,95) + rxt(k,228)*y(k,109) + rxt(k,208)*y(k,235)
         mat(k,2200) = rxt(k,225)*y(k,250)
         mat(k,1092) = rxt(k,558)*y(k,70)
         mat(k,1522) = rxt(k,206)*y(k,70)
         mat(k,693) = rxt(k,207)*y(k,70)
         mat(k,1736) = mat(k,1736) + rxt(k,228)*y(k,70)
         mat(k,2315) = rxt(k,208)*y(k,70)
         mat(k,2143) = mat(k,2143) + rxt(k,225)*y(k,74)
         mat(k,205) = -(rxt(k,359)*y(k,250) + rxt(k,367)*y(k,249))
         mat(k,2004) = -rxt(k,359)*y(k,102)
         mat(k,1943) = -rxt(k,367)*y(k,102)
         mat(k,1039) = -(rxt(k,360)*y(k,250))
         mat(k,2103) = -rxt(k,360)*y(k,103)
         mat(k,1006) = .050_r8*rxt(k,534)*y(k,164)
         mat(k,327) = .350_r8*rxt(k,370)*y(k,250)
         mat(k,653) = .370_r8*rxt(k,372)*y(k,164)
         mat(k,1125) = .120_r8*rxt(k,401)*y(k,164)
         mat(k,892) = .110_r8*rxt(k,479)*y(k,164)
         mat(k,1322) = .330_r8*rxt(k,432)*y(k,164)
         mat(k,978) = .050_r8*rxt(k,537)*y(k,164)
         mat(k,1427) = .120_r8*rxt(k,446)*y(k,164)
         mat(k,1884) = rxt(k,363)*y(k,236)
         mat(k,2558) = .050_r8*rxt(k,534)*y(k,6) + .370_r8*rxt(k,372)*y(k,28) &
                      + .120_r8*rxt(k,401)*y(k,33) + .110_r8*rxt(k,479)*y(k,127) &
                      + .330_r8*rxt(k,432)*y(k,134) + .050_r8*rxt(k,537)*y(k,139) &
                      + .120_r8*rxt(k,446)*y(k,140)
         mat(k,2282) = rxt(k,361)*y(k,236)
         mat(k,523) = rxt(k,363)*y(k,153) + rxt(k,361)*y(k,235)
         mat(k,2103) = mat(k,2103) + .350_r8*rxt(k,370)*y(k,27)
         mat(k,1559) = rxt(k,306)*y(k,89)
         mat(k,948) = rxt(k,306)*y(k,64) + rxt(k,307)*y(k,93) + rxt(k,309)*y(k,106) &
                      + rxt(k,308)*y(k,262)
         mat(k,1517) = rxt(k,307)*y(k,89)
         mat(k,1752) = rxt(k,309)*y(k,89)
         mat(k,2677) = rxt(k,308)*y(k,89)
         mat(k,1255) = -(rxt(k,267)*y(k,155) + rxt(k,295)*y(k,250) + (rxt(k,587) &
                      + rxt(k,661) + rxt(k,669) + rxt(k,678)) * y(k,109) + (rxt(k,588) &
                      + rxt(k,659) + rxt(k,672) + rxt(k,681)) * y(k,108) + (rxt(k,592) &
                      + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,110))
         mat(k,2637) = -rxt(k,267)*y(k,105)
         mat(k,2122) = -rxt(k,295)*y(k,105)
         mat(k,1731) = -(rxt(k,587) + rxt(k,661) + rxt(k,669) + rxt(k,678)) * y(k,105)
         mat(k,1580) = -(rxt(k,588) + rxt(k,659) + rxt(k,672) + rxt(k,681)) * y(k,105)
         mat(k,1656) = -(rxt(k,592) + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,105)
         mat(k,2430) = rxt(k,273)*y(k,235)
         mat(k,2295) = rxt(k,273)*y(k,115)
         mat(k,1760) = -(rxt(k,201)*y(k,250) + rxt(k,309)*y(k,89))
         mat(k,2145) = -rxt(k,201)*y(k,106)
         mat(k,953) = -rxt(k,309)*y(k,106)
         mat(k,1688) = rxt(k,349)*y(k,155)
         mat(k,1191) = rxt(k,381)*y(k,155)
         mat(k,1338) = rxt(k,407)*y(k,155)
         mat(k,1083) = (rxt(k,666)+rxt(k,675)+rxt(k,684))*y(k,101)
         mat(k,1093) = rxt(k,560)*y(k,155)
         mat(k,1714) = (rxt(k,666)+rxt(k,675)+rxt(k,684))*y(k,75) + rxt(k,601) &
                      *y(k,143)
         mat(k,1260) = rxt(k,267)*y(k,155)
         mat(k,1663) = rxt(k,297)*y(k,155)
         mat(k,408) = rxt(k,601)*y(k,101)
         mat(k,1817) = rxt(k,200)*y(k,250)
         mat(k,2658) = rxt(k,349)*y(k,51) + rxt(k,381)*y(k,54) + rxt(k,407)*y(k,58) &
                      + rxt(k,560)*y(k,83) + rxt(k,267)*y(k,105) + rxt(k,297)*y(k,110)
         mat(k,2145) = mat(k,2145) + rxt(k,200)*y(k,154)
         mat(k,516) = -(rxt(k,179)*y(k,250))
         mat(k,2051) = -rxt(k,179)*y(k,107)
         mat(k,1783) = rxt(k,198)*y(k,235)
         mat(k,2249) = rxt(k,198)*y(k,154)
         mat(k,1583) = -(rxt(k,261)*y(k,163) + (rxt(k,584) + rxt(k,658) + rxt(k,671) &
                      + rxt(k,680)) * y(k,97) + (rxt(k,586) + rxt(k,657) + rxt(k,670) &
                      + rxt(k,679)) * y(k,101) + (rxt(k,588) + rxt(k,659) + rxt(k,672) &
                      + rxt(k,681)) * y(k,105))
         mat(k,2488) = -rxt(k,261)*y(k,108)
         mat(k,1505) = -(rxt(k,584) + rxt(k,658) + rxt(k,671) + rxt(k,680)) * y(k,108)
         mat(k,1710) = -(rxt(k,586) + rxt(k,657) + rxt(k,670) + rxt(k,679)) * y(k,108)
         mat(k,1257) = -(rxt(k,588) + rxt(k,659) + rxt(k,672) + rxt(k,681)) * y(k,108)
         mat(k,544) = rxt(k,242)*y(k,250)
         mat(k,2169) = rxt(k,251)*y(k,235)
         mat(k,2311) = rxt(k,251)*y(k,21)
         mat(k,2139) = rxt(k,242)*y(k,18)
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
         mat(k,1737) = -(rxt(k,228)*y(k,70) + rxt(k,229)*y(k,163) + rxt(k,230) &
                      *y(k,250) + (rxt(k,583) + rxt(k,660) + rxt(k,668) + rxt(k,677) &
                      ) * y(k,97) + (rxt(k,585) + rxt(k,656) + rxt(k,667) + rxt(k,676) &
                      ) * y(k,101) + (rxt(k,587) + rxt(k,661) + rxt(k,669) + rxt(k,678) &
                      ) * y(k,105))
         mat(k,2377) = -rxt(k,228)*y(k,109)
         mat(k,2493) = -rxt(k,229)*y(k,109)
         mat(k,2144) = -rxt(k,230)*y(k,109)
         mat(k,1507) = -(rxt(k,583) + rxt(k,660) + rxt(k,668) + rxt(k,677)) * y(k,109)
         mat(k,1713) = -(rxt(k,585) + rxt(k,656) + rxt(k,667) + rxt(k,676)) * y(k,109)
         mat(k,1259) = -(rxt(k,587) + rxt(k,661) + rxt(k,669) + rxt(k,678)) * y(k,109)
         mat(k,1114) = rxt(k,335)*y(k,235)
         mat(k,623) = rxt(k,211)*y(k,250)
         mat(k,2201) = rxt(k,217)*y(k,235)
         mat(k,1082) = rxt(k,222)*y(k,250)
         mat(k,2316) = rxt(k,335)*y(k,68) + rxt(k,217)*y(k,74)
         mat(k,2144) = mat(k,2144) + rxt(k,211)*y(k,73) + rxt(k,222)*y(k,75)
         mat(k,1660) = -(rxt(k,268)*y(k,250) + rxt(k,297)*y(k,155) + (rxt(k,590) &
                      + rxt(k,686) + rxt(k,690) + rxt(k,694)) * y(k,101) + (rxt(k,591) &
                      + rxt(k,687) + rxt(k,691) + rxt(k,695)) * y(k,97) + (rxt(k,592) &
                      + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,105))
         mat(k,2141) = -rxt(k,268)*y(k,110)
         mat(k,2654) = -rxt(k,297)*y(k,110)
         mat(k,1711) = -(rxt(k,590) + rxt(k,686) + rxt(k,690) + rxt(k,694)) * y(k,110)
         mat(k,1506) = -(rxt(k,591) + rxt(k,687) + rxt(k,691) + rxt(k,695)) * y(k,110)
         mat(k,1258) = -(rxt(k,592) + rxt(k,688) + rxt(k,692) + rxt(k,696)) * y(k,110)
         mat(k,1533) = rxt(k,271)*y(k,250)
         mat(k,2405) = rxt(k,288)*y(k,235)
         mat(k,2313) = rxt(k,288)*y(k,125)
         mat(k,2141) = mat(k,2141) + rxt(k,271)*y(k,116)
         mat(k,1234) = -(rxt(k,425)*y(k,250))
         mat(k,2120) = -rxt(k,425)*y(k,111)
         mat(k,676) = .300_r8*rxt(k,470)*y(k,250)
         mat(k,630) = .500_r8*rxt(k,471)*y(k,250)
         mat(k,1897) = rxt(k,424)*y(k,232) + rxt(k,431)*y(k,242)
         mat(k,662) = rxt(k,424)*y(k,153)
         mat(k,1409) = rxt(k,431)*y(k,153)
         mat(k,2120) = mat(k,2120) + .300_r8*rxt(k,470)*y(k,128) + .500_r8*rxt(k,471) &
                      *y(k,129)
         mat(k,274) = -(rxt(k,456)*y(k,250))
         mat(k,2016) = -rxt(k,456)*y(k,112)
         mat(k,1247) = -(rxt(k,410)*y(k,250))
         mat(k,2121) = -rxt(k,410)*y(k,113)
         mat(k,677) = .700_r8*rxt(k,470)*y(k,250)
         mat(k,631) = .500_r8*rxt(k,471)*y(k,250)
         mat(k,685) = .500_r8*rxt(k,445)*y(k,250)
         mat(k,1898) = .050_r8*rxt(k,468)*y(k,238) + .220_r8*rxt(k,430)*y(k,242) &
                      + .250_r8*rxt(k,487)*y(k,258)
         mat(k,2636) = .050_r8*rxt(k,469)*y(k,238) + .220_r8*rxt(k,429)*y(k,242) &
                      + .250_r8*rxt(k,488)*y(k,258)
         mat(k,646) = .500_r8*rxt(k,414)*y(k,250)
         mat(k,1478) = .220_r8*rxt(k,426)*y(k,242) + .250_r8*rxt(k,484)*y(k,258)
         mat(k,1624) = .230_r8*rxt(k,427)*y(k,242) + .200_r8*rxt(k,415)*y(k,253) &
                      + .100_r8*rxt(k,485)*y(k,258)
         mat(k,1384) = .050_r8*rxt(k,468)*y(k,153) + .050_r8*rxt(k,469)*y(k,155)
         mat(k,1410) = .220_r8*rxt(k,430)*y(k,153) + .220_r8*rxt(k,429)*y(k,155) &
                      + .220_r8*rxt(k,426)*y(k,229) + .230_r8*rxt(k,427)*y(k,230)
         mat(k,2121) = mat(k,2121) + .700_r8*rxt(k,470)*y(k,128) + .500_r8*rxt(k,471) &
                      *y(k,129) + .500_r8*rxt(k,445)*y(k,138) + .500_r8*rxt(k,414) &
                      *y(k,177)
         mat(k,1270) = .200_r8*rxt(k,415)*y(k,230)
         mat(k,1287) = .250_r8*rxt(k,487)*y(k,153) + .250_r8*rxt(k,488)*y(k,155) &
                      + .250_r8*rxt(k,484)*y(k,229) + .100_r8*rxt(k,485)*y(k,230)
         mat(k,391) = -(rxt(k,457)*y(k,250))
         mat(k,2034) = -rxt(k,457)*y(k,114)
         mat(k,1852) = .870_r8*rxt(k,468)*y(k,238)
         mat(k,2607) = .950_r8*rxt(k,469)*y(k,238)
         mat(k,1469) = rxt(k,464)*y(k,238)
         mat(k,1601) = .750_r8*rxt(k,465)*y(k,238)
         mat(k,1373) = .870_r8*rxt(k,468)*y(k,153) + .950_r8*rxt(k,469)*y(k,155) &
                      + rxt(k,464)*y(k,229) + .750_r8*rxt(k,465)*y(k,230)
         mat(k,2446) = -(rxt(k,272)*y(k,21) + rxt(k,273)*y(k,235) + rxt(k,274) &
                      *y(k,126) + rxt(k,276)*y(k,154) + rxt(k,278)*y(k,155) + rxt(k,280) &
                      *y(k,153) + rxt(k,281)*y(k,164))
         mat(k,2183) = -rxt(k,272)*y(k,115)
         mat(k,2327) = -rxt(k,273)*y(k,115)
         mat(k,945) = -rxt(k,274)*y(k,115)
         mat(k,1827) = -rxt(k,276)*y(k,115)
         mat(k,2668) = -rxt(k,278)*y(k,115)
         mat(k,1927) = -rxt(k,280)*y(k,115)
         mat(k,2598) = -rxt(k,281)*y(k,115)
         mat(k,2532) = rxt(k,282)*y(k,125)
         mat(k,2183) = mat(k,2183) + rxt(k,283)*y(k,125)
         mat(k,438) = rxt(k,326)*y(k,70) + rxt(k,327)*y(k,250)
         mat(k,2388) = rxt(k,326)*y(k,61)
         mat(k,2212) = (rxt(k,285)+rxt(k,286))*y(k,125)
         mat(k,1099) = rxt(k,559)*y(k,125)
         mat(k,1263) = rxt(k,267)*y(k,155) + rxt(k,295)*y(k,250)
         mat(k,1539) = rxt(k,269)*y(k,155) + rxt(k,270)*y(k,163) + rxt(k,271)*y(k,250)
         mat(k,2418) = rxt(k,282)*y(k,17) + rxt(k,283)*y(k,21) + (rxt(k,285) &
                       +rxt(k,286))*y(k,74) + rxt(k,559)*y(k,83) + 2.000_r8*rxt(k,301) &
                      *y(k,125) + rxt(k,289)*y(k,153) + rxt(k,292)*y(k,163) &
                      + rxt(k,294)*y(k,250)
         mat(k,1927) = mat(k,1927) + rxt(k,289)*y(k,125)
         mat(k,2668) = mat(k,2668) + rxt(k,267)*y(k,105) + rxt(k,269)*y(k,116)
         mat(k,2504) = rxt(k,270)*y(k,116) + rxt(k,292)*y(k,125)
         mat(k,2155) = rxt(k,327)*y(k,61) + rxt(k,295)*y(k,105) + rxt(k,271)*y(k,116) &
                      + rxt(k,294)*y(k,125)
         mat(k,1532) = -(rxt(k,269)*y(k,155) + rxt(k,270)*y(k,163) + rxt(k,271) &
                      *y(k,250))
         mat(k,2649) = -rxt(k,269)*y(k,116)
         mat(k,2485) = -rxt(k,270)*y(k,116)
         mat(k,2136) = -rxt(k,271)*y(k,116)
         mat(k,1256) = (rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,110)
         mat(k,1658) = (rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,105)
         mat(k,2431) = rxt(k,274)*y(k,126)
         mat(k,210) = 2.000_r8*rxt(k,279)*y(k,123)
         mat(k,313) = 2.000_r8*rxt(k,275)*y(k,124)
         mat(k,939) = rxt(k,274)*y(k,115)
         mat(k,2395) = 2.000_r8*rxt(k,302)*y(k,125)
         mat(k,2396) = rxt(k,304)*y(k,168)
         mat(k,594) = rxt(k,304)*y(k,125)
         mat(k,593) = 2.000_r8*rxt(k,305)*y(k,168)
         mat(k,1501) = (rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,110)
         mat(k,1253) = (rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,108)
         mat(k,1577) = (rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,105)
         mat(k,1654) = (rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,97)
         mat(k,2193) = rxt(k,287)*y(k,125)
         mat(k,1706) = (rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,110)
         mat(k,1254) = (rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,109)
         mat(k,1729) = (rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,105)
         mat(k,1655) = (rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,101)
         mat(k,2398) = rxt(k,287)*y(k,74)
         mat(k,163) = -(rxt(k,458)*y(k,250))
         mat(k,2000) = -rxt(k,458)*y(k,122)
         mat(k,827) = .600_r8*rxt(k,481)*y(k,250)
         mat(k,2000) = mat(k,2000) + .600_r8*rxt(k,481)*y(k,131)
         mat(k,209) = -(4._r8*rxt(k,279)*y(k,123))
         mat(k,2425) = rxt(k,280)*y(k,153)
         mat(k,1846) = rxt(k,280)*y(k,115)
         mat(k,310) = -(4._r8*rxt(k,275)*y(k,124))
         mat(k,2426) = rxt(k,276)*y(k,154)
         mat(k,1779) = rxt(k,276)*y(k,115)
         mat(k,2417) = -(rxt(k,282)*y(k,17) + (rxt(k,283) + rxt(k,284)) * y(k,21) &
                      + (rxt(k,285) + rxt(k,286) + rxt(k,287)) * y(k,74) + rxt(k,288) &
                      *y(k,235) + rxt(k,289)*y(k,153) + rxt(k,290)*y(k,154) + rxt(k,291) &
                      *y(k,155) + rxt(k,292)*y(k,163) + rxt(k,293)*y(k,164) + rxt(k,294) &
                      *y(k,250) + (4._r8*rxt(k,301) + 4._r8*rxt(k,302)) * y(k,125) &
                      + rxt(k,304)*y(k,168) + rxt(k,559)*y(k,83))
         mat(k,2531) = -rxt(k,282)*y(k,125)
         mat(k,2182) = -(rxt(k,283) + rxt(k,284)) * y(k,125)
         mat(k,2211) = -(rxt(k,285) + rxt(k,286) + rxt(k,287)) * y(k,125)
         mat(k,2326) = -rxt(k,288)*y(k,125)
         mat(k,1926) = -rxt(k,289)*y(k,125)
         mat(k,1826) = -rxt(k,290)*y(k,125)
         mat(k,2667) = -rxt(k,291)*y(k,125)
         mat(k,2503) = -rxt(k,292)*y(k,125)
         mat(k,2597) = -rxt(k,293)*y(k,125)
         mat(k,2154) = -rxt(k,294)*y(k,125)
         mat(k,598) = -rxt(k,304)*y(k,125)
         mat(k,1098) = -rxt(k,559)*y(k,125)
         mat(k,2182) = mat(k,2182) + rxt(k,272)*y(k,115)
         mat(k,1669) = rxt(k,297)*y(k,155) + rxt(k,268)*y(k,250)
         mat(k,2445) = rxt(k,272)*y(k,21) + rxt(k,278)*y(k,155) + rxt(k,281)*y(k,164)
         mat(k,1538) = rxt(k,270)*y(k,163)
         mat(k,1926) = mat(k,1926) + rxt(k,296)*y(k,168)
         mat(k,2667) = mat(k,2667) + rxt(k,297)*y(k,110) + rxt(k,278)*y(k,115)
         mat(k,2503) = mat(k,2503) + rxt(k,270)*y(k,116)
         mat(k,2597) = mat(k,2597) + rxt(k,281)*y(k,115)
         mat(k,598) = mat(k,598) + rxt(k,296)*y(k,153)
         mat(k,2154) = mat(k,2154) + rxt(k,268)*y(k,110)
         mat(k,938) = -(rxt(k,274)*y(k,115))
         mat(k,2429) = -rxt(k,274)*y(k,126)
         mat(k,1531) = rxt(k,269)*y(k,155)
         mat(k,2400) = rxt(k,290)*y(k,154)
         mat(k,1798) = rxt(k,290)*y(k,125)
         mat(k,2615) = rxt(k,269)*y(k,116)
         mat(k,891) = -(rxt(k,472)*y(k,155) + rxt(k,479)*y(k,164) + rxt(k,480) &
                      *y(k,250))
         mat(k,2613) = -rxt(k,472)*y(k,127)
         mat(k,2554) = -rxt(k,479)*y(k,127)
         mat(k,2095) = -rxt(k,480)*y(k,127)
         mat(k,674) = -(rxt(k,470)*y(k,250))
         mat(k,2071) = -rxt(k,470)*y(k,128)
         mat(k,1866) = .080_r8*rxt(k,462)*y(k,237)
         mat(k,1344) = .080_r8*rxt(k,462)*y(k,153)
         mat(k,627) = -(rxt(k,471)*y(k,250))
         mat(k,2065) = -rxt(k,471)*y(k,129)
         mat(k,1864) = .080_r8*rxt(k,468)*y(k,238)
         mat(k,1374) = .080_r8*rxt(k,468)*y(k,153)
         mat(k,490) = -(rxt(k,478)*y(k,250))
         mat(k,2047) = -rxt(k,478)*y(k,130)
         mat(k,2245) = rxt(k,475)*y(k,239)
         mat(k,1300) = rxt(k,475)*y(k,235)
         mat(k,828) = -(rxt(k,481)*y(k,250))
         mat(k,2088) = -rxt(k,481)*y(k,131)
         mat(k,2272) = rxt(k,461)*y(k,237) + rxt(k,466)*y(k,238)
         mat(k,1345) = rxt(k,461)*y(k,235)
         mat(k,1376) = rxt(k,466)*y(k,235)
         mat(k,85) = -(rxt(k,641)*y(k,250))
         mat(k,1990) = -rxt(k,641)*y(k,132)
         mat(k,1324) = -(rxt(k,432)*y(k,164) + rxt(k,433)*y(k,250))
         mat(k,2574) = -rxt(k,432)*y(k,134)
         mat(k,2126) = -rxt(k,433)*y(k,134)
         mat(k,896) = .300_r8*rxt(k,479)*y(k,164)
         mat(k,1902) = .360_r8*rxt(k,462)*y(k,237)
         mat(k,2641) = .400_r8*rxt(k,463)*y(k,237)
         mat(k,2574) = mat(k,2574) + .300_r8*rxt(k,479)*y(k,127)
         mat(k,1481) = .390_r8*rxt(k,459)*y(k,237)
         mat(k,1628) = .310_r8*rxt(k,460)*y(k,237)
         mat(k,1354) = .360_r8*rxt(k,462)*y(k,153) + .400_r8*rxt(k,463)*y(k,155) &
                      + .390_r8*rxt(k,459)*y(k,229) + .310_r8*rxt(k,460)*y(k,230)
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
         mat(k,353) = -(rxt(k,434)*y(k,250))
         mat(k,2029) = -rxt(k,434)*y(k,135)
         mat(k,2235) = rxt(k,428)*y(k,242)
         mat(k,1405) = rxt(k,428)*y(k,235)
         mat(k,600) = -(rxt(k,443)*y(k,250))
         mat(k,2061) = -rxt(k,443)*y(k,136)
         mat(k,1862) = .800_r8*rxt(k,452)*y(k,221)
         mat(k,1023) = .800_r8*rxt(k,452)*y(k,153)
         mat(k,358) = -(rxt(k,444)*y(k,250))
         mat(k,2030) = -rxt(k,444)*y(k,137)
         mat(k,2236) = .800_r8*rxt(k,441)*y(k,246)
         mat(k,769) = .800_r8*rxt(k,441)*y(k,235)
         mat(k,683) = -(rxt(k,445)*y(k,250))
         mat(k,2072) = -rxt(k,445)*y(k,138)
         mat(k,1792) = rxt(k,448)*y(k,244)
         mat(k,1450) = rxt(k,448)*y(k,154)
         mat(k,976) = -(rxt(k,536)*y(k,155) + rxt(k,537)*y(k,164) + rxt(k,538) &
                      *y(k,250))
         mat(k,2616) = -rxt(k,536)*y(k,139)
         mat(k,2555) = -rxt(k,537)*y(k,139)
         mat(k,2100) = -rxt(k,538)*y(k,139)
         mat(k,1434) = -(rxt(k,446)*y(k,164) + rxt(k,447)*y(k,250))
         mat(k,2579) = -rxt(k,446)*y(k,140)
         mat(k,2131) = -rxt(k,447)*y(k,140)
         mat(k,899) = .200_r8*rxt(k,479)*y(k,164)
         mat(k,1907) = .560_r8*rxt(k,462)*y(k,237)
         mat(k,2646) = .600_r8*rxt(k,463)*y(k,237)
         mat(k,2579) = mat(k,2579) + .200_r8*rxt(k,479)*y(k,127)
         mat(k,1486) = .610_r8*rxt(k,459)*y(k,237)
         mat(k,1633) = .440_r8*rxt(k,460)*y(k,237)
         mat(k,1358) = .560_r8*rxt(k,462)*y(k,153) + .600_r8*rxt(k,463)*y(k,155) &
                      + .610_r8*rxt(k,459)*y(k,229) + .440_r8*rxt(k,460)*y(k,230)
         mat(k,569) = -(rxt(k,180)*y(k,153) + (rxt(k,181) + rxt(k,182) + rxt(k,183) &
                      ) * y(k,154) + rxt(k,192)*y(k,250))
         mat(k,1858) = -rxt(k,180)*y(k,141)
         mat(k,1787) = -(rxt(k,181) + rxt(k,182) + rxt(k,183)) * y(k,141)
         mat(k,2058) = -rxt(k,192)*y(k,141)
         mat(k,213) = -((rxt(k,196) + rxt(k,197)) * y(k,249))
         mat(k,1944) = -(rxt(k,196) + rxt(k,197)) * y(k,142)
         mat(k,568) = rxt(k,181)*y(k,154)
         mat(k,1778) = rxt(k,181)*y(k,141)
         mat(k,1781) = rxt(k,199)*y(k,155)
         mat(k,2608) = rxt(k,199)*y(k,154)
         mat(k,454) = -(rxt(k,482)*y(k,250))
         mat(k,2042) = -rxt(k,482)*y(k,144)
         mat(k,1603) = .200_r8*rxt(k,474)*y(k,239)
         mat(k,1299) = .200_r8*rxt(k,474)*y(k,230)
         mat(k,1175) = -(rxt(k,483)*y(k,250))
         mat(k,2115) = -rxt(k,483)*y(k,145)
         mat(k,1892) = rxt(k,476)*y(k,239)
         mat(k,2630) = rxt(k,477)*y(k,239)
         mat(k,1475) = rxt(k,473)*y(k,239)
         mat(k,1619) = .800_r8*rxt(k,474)*y(k,239)
         mat(k,1304) = rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155) + rxt(k,473)*y(k,229) &
                      + .800_r8*rxt(k,474)*y(k,230)
         mat(k,109) = -(rxt(k,594)*y(k,250))
         mat(k,1994) = -rxt(k,594)*y(k,149)
         mat(k,1919) = -(rxt(k,180)*y(k,141) + rxt(k,189)*y(k,155) + rxt(k,193) &
                      *y(k,235) + rxt(k,194)*y(k,164) + rxt(k,195)*y(k,163) + rxt(k,218) &
                      *y(k,74) + rxt(k,252)*y(k,21) + rxt(k,280)*y(k,115) + rxt(k,289) &
                      *y(k,125) + rxt(k,296)*y(k,168) + rxt(k,336)*y(k,68) + rxt(k,355) &
                      *y(k,230) + rxt(k,363)*y(k,236) + rxt(k,376)*y(k,226) + rxt(k,387) &
                      *y(k,229) + rxt(k,391)*y(k,234) + rxt(k,404)*y(k,227) + rxt(k,413) &
                      *y(k,252) + rxt(k,417)*y(k,253) + (rxt(k,423) + rxt(k,424) &
                      ) * y(k,232) + (rxt(k,430) + rxt(k,431)) * y(k,242) + rxt(k,439) &
                      *y(k,244) + rxt(k,442)*y(k,246) + (rxt(k,452) + rxt(k,453) &
                      ) * y(k,221) + rxt(k,462)*y(k,237) + rxt(k,468)*y(k,238) &
                      + rxt(k,476)*y(k,239) + rxt(k,487)*y(k,258) + rxt(k,491) &
                      *y(k,220) + rxt(k,494)*y(k,223) + rxt(k,499)*y(k,225) + rxt(k,501) &
                      *y(k,228) + rxt(k,505)*y(k,231) + rxt(k,508)*y(k,243) + rxt(k,511) &
                      *y(k,245) + rxt(k,514)*y(k,251) + rxt(k,521)*y(k,256) + rxt(k,527) &
                      *y(k,259) + rxt(k,530)*y(k,261) + rxt(k,541)*y(k,248) + rxt(k,546) &
                      *y(k,254) + rxt(k,551)*y(k,255))
         mat(k,571) = -rxt(k,180)*y(k,153)
         mat(k,2660) = -rxt(k,189)*y(k,153)
         mat(k,2319) = -rxt(k,193)*y(k,153)
         mat(k,2590) = -rxt(k,194)*y(k,153)
         mat(k,2496) = -rxt(k,195)*y(k,153)
         mat(k,2204) = -rxt(k,218)*y(k,153)
         mat(k,2175) = -rxt(k,252)*y(k,153)
         mat(k,2438) = -rxt(k,280)*y(k,153)
         mat(k,2410) = -rxt(k,289)*y(k,153)
         mat(k,597) = -rxt(k,296)*y(k,153)
         mat(k,1116) = -rxt(k,336)*y(k,153)
         mat(k,1643) = -rxt(k,355)*y(k,153)
         mat(k,526) = -rxt(k,363)*y(k,153)
         mat(k,885) = -rxt(k,376)*y(k,153)
         mat(k,1494) = -rxt(k,387)*y(k,153)
         mat(k,803) = -rxt(k,391)*y(k,153)
         mat(k,932) = -rxt(k,404)*y(k,153)
         mat(k,866) = -rxt(k,413)*y(k,153)
         mat(k,1277) = -rxt(k,417)*y(k,153)
         mat(k,665) = -(rxt(k,423) + rxt(k,424)) * y(k,153)
         mat(k,1420) = -(rxt(k,430) + rxt(k,431)) * y(k,153)
         mat(k,1462) = -rxt(k,439)*y(k,153)
         mat(k,774) = -rxt(k,442)*y(k,153)
         mat(k,1035) = -(rxt(k,452) + rxt(k,453)) * y(k,153)
         mat(k,1365) = -rxt(k,462)*y(k,153)
         mat(k,1398) = -rxt(k,468)*y(k,153)
         mat(k,1317) = -rxt(k,476)*y(k,153)
         mat(k,1294) = -rxt(k,487)*y(k,153)
         mat(k,618) = -rxt(k,491)*y(k,153)
         mat(k,582) = -rxt(k,494)*y(k,153)
         mat(k,513) = -rxt(k,499)*y(k,153)
         mat(k,733) = -rxt(k,501)*y(k,153)
         mat(k,857) = -rxt(k,505)*y(k,153)
         mat(k,810) = -rxt(k,508)*y(k,153)
         mat(k,964) = -rxt(k,511)*y(k,153)
         mat(k,532) = -rxt(k,514)*y(k,153)
         mat(k,824) = -rxt(k,521)*y(k,153)
         mat(k,849) = -rxt(k,527)*y(k,153)
         mat(k,590) = -rxt(k,530)*y(k,153)
         mat(k,1155) = -rxt(k,541)*y(k,153)
         mat(k,1224) = -rxt(k,546)*y(k,153)
         mat(k,1062) = -rxt(k,551)*y(k,153)
         mat(k,211) = 4.000_r8*rxt(k,279)*y(k,123)
         mat(k,571) = mat(k,571) + 2.000_r8*rxt(k,182)*y(k,154) + rxt(k,192)*y(k,250)
         mat(k,214) = 2.000_r8*rxt(k,196)*y(k,249)
         mat(k,1819) = 2.000_r8*rxt(k,182)*y(k,141) + rxt(k,185)*y(k,163) + rxt(k,570) &
                      *y(k,181)
         mat(k,2496) = mat(k,2496) + rxt(k,185)*y(k,154)
         mat(k,1548) = rxt(k,570)*y(k,154)
         mat(k,1965) = 2.000_r8*rxt(k,196)*y(k,142)
         mat(k,2147) = rxt(k,192)*y(k,141)
         mat(k,1818) = -((rxt(k,181) + rxt(k,182) + rxt(k,183)) * y(k,141) + (rxt(k,185) &
                      + rxt(k,187)) * y(k,163) + rxt(k,186)*y(k,164) + rxt(k,198) &
                      *y(k,235) + rxt(k,199)*y(k,155) + rxt(k,200)*y(k,250) + rxt(k,210) &
                      *y(k,70) + rxt(k,220)*y(k,74) + rxt(k,245)*y(k,17) + rxt(k,255) &
                      *y(k,21) + rxt(k,276)*y(k,115) + rxt(k,290)*y(k,125) + rxt(k,398) &
                      *y(k,229) + rxt(k,448)*y(k,244) + rxt(k,506)*y(k,231) + rxt(k,509) &
                      *y(k,243) + rxt(k,512)*y(k,245) + rxt(k,516)*y(k,172) + rxt(k,519) &
                      *y(k,220) + rxt(k,570)*y(k,181))
         mat(k,570) = -(rxt(k,181) + rxt(k,182) + rxt(k,183)) * y(k,154)
         mat(k,2495) = -(rxt(k,185) + rxt(k,187)) * y(k,154)
         mat(k,2589) = -rxt(k,186)*y(k,154)
         mat(k,2318) = -rxt(k,198)*y(k,154)
         mat(k,2659) = -rxt(k,199)*y(k,154)
         mat(k,2146) = -rxt(k,200)*y(k,154)
         mat(k,2379) = -rxt(k,210)*y(k,154)
         mat(k,2203) = -rxt(k,220)*y(k,154)
         mat(k,2523) = -rxt(k,245)*y(k,154)
         mat(k,2174) = -rxt(k,255)*y(k,154)
         mat(k,2437) = -rxt(k,276)*y(k,154)
         mat(k,2409) = -rxt(k,290)*y(k,154)
         mat(k,1493) = -rxt(k,398)*y(k,154)
         mat(k,1461) = -rxt(k,448)*y(k,154)
         mat(k,856) = -rxt(k,506)*y(k,154)
         mat(k,809) = -rxt(k,509)*y(k,154)
         mat(k,963) = -rxt(k,512)*y(k,154)
         mat(k,549) = -rxt(k,516)*y(k,154)
         mat(k,617) = -rxt(k,519)*y(k,154)
         mat(k,1547) = -rxt(k,570)*y(k,154)
         mat(k,745) = rxt(k,450)*y(k,250)
         mat(k,418) = rxt(k,421)*y(k,155)
         mat(k,2174) = mat(k,2174) + rxt(k,252)*y(k,153)
         mat(k,1115) = rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155)
         mat(k,624) = rxt(k,211)*y(k,250)
         mat(k,2203) = mat(k,2203) + rxt(k,218)*y(k,153)
         mat(k,517) = rxt(k,179)*y(k,250)
         mat(k,2437) = mat(k,2437) + rxt(k,278)*y(k,155)
         mat(k,314) = 4.000_r8*rxt(k,275)*y(k,124)
         mat(k,2409) = mat(k,2409) + rxt(k,289)*y(k,153) + rxt(k,291)*y(k,155)
         mat(k,680) = .700_r8*rxt(k,470)*y(k,250)
         mat(k,1918) = rxt(k,252)*y(k,21) + rxt(k,336)*y(k,68) + rxt(k,218)*y(k,74) &
                      + rxt(k,289)*y(k,125) + 2.000_r8*rxt(k,189)*y(k,155) &
                      + rxt(k,195)*y(k,163) + rxt(k,194)*y(k,164) + rxt(k,296) &
                      *y(k,168) + rxt(k,491)*y(k,220) + rxt(k,452)*y(k,221) &
                      + rxt(k,494)*y(k,223) + rxt(k,499)*y(k,225) + rxt(k,376) &
                      *y(k,226) + rxt(k,404)*y(k,227) + rxt(k,501)*y(k,228) &
                      + rxt(k,387)*y(k,229) + rxt(k,355)*y(k,230) + rxt(k,505) &
                      *y(k,231) + rxt(k,423)*y(k,232) + rxt(k,391)*y(k,234) &
                      + rxt(k,193)*y(k,235) + rxt(k,363)*y(k,236) + .920_r8*rxt(k,462) &
                      *y(k,237) + .920_r8*rxt(k,468)*y(k,238) + rxt(k,476)*y(k,239) &
                      + rxt(k,430)*y(k,242) + rxt(k,508)*y(k,243) + rxt(k,439) &
                      *y(k,244) + rxt(k,511)*y(k,245) + rxt(k,442)*y(k,246) &
                      + 1.600_r8*rxt(k,541)*y(k,248) + rxt(k,514)*y(k,251) &
                      + rxt(k,413)*y(k,252) + rxt(k,417)*y(k,253) + .900_r8*rxt(k,546) &
                      *y(k,254) + .800_r8*rxt(k,551)*y(k,255) + rxt(k,521)*y(k,256) &
                      + rxt(k,487)*y(k,258) + rxt(k,527)*y(k,259) + rxt(k,530) &
                      *y(k,261)
         mat(k,2659) = mat(k,2659) + rxt(k,421)*y(k,16) + rxt(k,337)*y(k,68) &
                      + rxt(k,278)*y(k,115) + rxt(k,291)*y(k,125) &
                      + 2.000_r8*rxt(k,189)*y(k,153) + rxt(k,190)*y(k,163) &
                      + rxt(k,188)*y(k,235) + rxt(k,463)*y(k,237) + rxt(k,469) &
                      *y(k,238) + rxt(k,477)*y(k,239) + rxt(k,429)*y(k,242) &
                      + rxt(k,440)*y(k,244) + 2.000_r8*rxt(k,542)*y(k,248) &
                      + rxt(k,191)*y(k,250) + rxt(k,488)*y(k,258)
         mat(k,911) = rxt(k,411)*y(k,250)
         mat(k,2495) = mat(k,2495) + rxt(k,195)*y(k,153) + rxt(k,190)*y(k,155)
         mat(k,2589) = mat(k,2589) + rxt(k,194)*y(k,153)
         mat(k,596) = rxt(k,296)*y(k,153)
         mat(k,726) = rxt(k,548)*y(k,250)
         mat(k,617) = mat(k,617) + rxt(k,491)*y(k,153)
         mat(k,1034) = rxt(k,452)*y(k,153)
         mat(k,581) = rxt(k,494)*y(k,153)
         mat(k,512) = rxt(k,499)*y(k,153)
         mat(k,884) = rxt(k,376)*y(k,153)
         mat(k,931) = rxt(k,404)*y(k,153)
         mat(k,732) = rxt(k,501)*y(k,153)
         mat(k,1493) = mat(k,1493) + rxt(k,387)*y(k,153)
         mat(k,1642) = rxt(k,355)*y(k,153) + .500_r8*rxt(k,539)*y(k,248)
         mat(k,856) = mat(k,856) + rxt(k,505)*y(k,153)
         mat(k,664) = rxt(k,423)*y(k,153)
         mat(k,802) = rxt(k,391)*y(k,153)
         mat(k,2318) = mat(k,2318) + rxt(k,193)*y(k,153) + rxt(k,188)*y(k,155)
         mat(k,525) = rxt(k,363)*y(k,153)
         mat(k,1364) = .920_r8*rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155)
         mat(k,1397) = .920_r8*rxt(k,468)*y(k,153) + rxt(k,469)*y(k,155)
         mat(k,1316) = rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155)
         mat(k,1419) = rxt(k,430)*y(k,153) + rxt(k,429)*y(k,155)
         mat(k,809) = mat(k,809) + rxt(k,508)*y(k,153)
         mat(k,1461) = mat(k,1461) + rxt(k,439)*y(k,153) + rxt(k,440)*y(k,155)
         mat(k,963) = mat(k,963) + rxt(k,511)*y(k,153)
         mat(k,773) = rxt(k,442)*y(k,153)
         mat(k,1154) = 1.600_r8*rxt(k,541)*y(k,153) + 2.000_r8*rxt(k,542)*y(k,155) &
                      + .500_r8*rxt(k,539)*y(k,230)
         mat(k,2146) = mat(k,2146) + rxt(k,450)*y(k,1) + rxt(k,211)*y(k,73) &
                      + rxt(k,179)*y(k,107) + .700_r8*rxt(k,470)*y(k,128) + rxt(k,191) &
                      *y(k,155) + rxt(k,411)*y(k,156) + rxt(k,548)*y(k,207)
         mat(k,531) = rxt(k,514)*y(k,153)
         mat(k,865) = rxt(k,413)*y(k,153)
         mat(k,1276) = rxt(k,417)*y(k,153)
         mat(k,1223) = .900_r8*rxt(k,546)*y(k,153)
         mat(k,1061) = .800_r8*rxt(k,551)*y(k,153)
         mat(k,823) = rxt(k,521)*y(k,153)
         mat(k,1293) = rxt(k,487)*y(k,153) + rxt(k,488)*y(k,155)
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
         mat(k,2673) = -(rxt(k,188)*y(k,235) + rxt(k,189)*y(k,153) + rxt(k,190) &
                      *y(k,163) + rxt(k,191)*y(k,250) + rxt(k,199)*y(k,154) + rxt(k,267) &
                      *y(k,105) + rxt(k,269)*y(k,116) + rxt(k,278)*y(k,115) + rxt(k,291) &
                      *y(k,125) + rxt(k,297)*y(k,110) + rxt(k,337)*y(k,68) + rxt(k,349) &
                      *y(k,51) + rxt(k,381)*y(k,54) + rxt(k,400)*y(k,33) + rxt(k,407) &
                      *y(k,58) + rxt(k,421)*y(k,16) + rxt(k,429)*y(k,242) + rxt(k,440) &
                      *y(k,244) + rxt(k,463)*y(k,237) + rxt(k,469)*y(k,238) + rxt(k,472) &
                      *y(k,127) + rxt(k,477)*y(k,239) + rxt(k,488)*y(k,258) + rxt(k,533) &
                      *y(k,6) + rxt(k,536)*y(k,139) + rxt(k,542)*y(k,248) + rxt(k,553) &
                      *y(k,209) + rxt(k,560)*y(k,83))
         mat(k,2332) = -rxt(k,188)*y(k,155)
         mat(k,1932) = -rxt(k,189)*y(k,155)
         mat(k,2509) = -rxt(k,190)*y(k,155)
         mat(k,2160) = -rxt(k,191)*y(k,155)
         mat(k,1832) = -rxt(k,199)*y(k,155)
         mat(k,1266) = -rxt(k,267)*y(k,155)
         mat(k,1542) = -rxt(k,269)*y(k,155)
         mat(k,2451) = -rxt(k,278)*y(k,155)
         mat(k,2423) = -rxt(k,291)*y(k,155)
         mat(k,1674) = -rxt(k,297)*y(k,155)
         mat(k,1119) = -rxt(k,337)*y(k,155)
         mat(k,1703) = -rxt(k,349)*y(k,155)
         mat(k,1194) = -rxt(k,381)*y(k,155)
         mat(k,1141) = -rxt(k,400)*y(k,155)
         mat(k,1342) = -rxt(k,407)*y(k,155)
         mat(k,420) = -rxt(k,421)*y(k,155)
         mat(k,1424) = -rxt(k,429)*y(k,155)
         mat(k,1467) = -rxt(k,440)*y(k,155)
         mat(k,1370) = -rxt(k,463)*y(k,155)
         mat(k,1403) = -rxt(k,469)*y(k,155)
         mat(k,906) = -rxt(k,472)*y(k,155)
         mat(k,1321) = -rxt(k,477)*y(k,155)
         mat(k,1298) = -rxt(k,488)*y(k,155)
         mat(k,1022) = -rxt(k,533)*y(k,155)
         mat(k,994) = -rxt(k,536)*y(k,155)
         mat(k,1158) = -rxt(k,542)*y(k,155)
         mat(k,1072) = -rxt(k,553)*y(k,155)
         mat(k,1101) = -rxt(k,560)*y(k,155)
         mat(k,2537) = rxt(k,253)*y(k,22)
         mat(k,924) = rxt(k,253)*y(k,17) + rxt(k,254)*y(k,70) + rxt(k,256)*y(k,163)
         mat(k,2393) = rxt(k,254)*y(k,22) + rxt(k,219)*y(k,75)
         mat(k,1089) = rxt(k,219)*y(k,70) + rxt(k,221)*y(k,163) + rxt(k,222)*y(k,250)
         mat(k,956) = rxt(k,309)*y(k,106)
         mat(k,1775) = rxt(k,309)*y(k,89) + rxt(k,201)*y(k,250)
         mat(k,2451) = mat(k,2451) + rxt(k,274)*y(k,126)
         mat(k,947) = rxt(k,274)*y(k,115)
         mat(k,691) = .500_r8*rxt(k,445)*y(k,250)
         mat(k,1832) = mat(k,1832) + rxt(k,187)*y(k,163) + rxt(k,186)*y(k,164)
         mat(k,2509) = mat(k,2509) + rxt(k,256)*y(k,22) + rxt(k,221)*y(k,75) &
                      + rxt(k,187)*y(k,154)
         mat(k,2603) = rxt(k,186)*y(k,154)
         mat(k,642) = rxt(k,396)*y(k,250)
         mat(k,2160) = mat(k,2160) + rxt(k,222)*y(k,75) + rxt(k,201)*y(k,106) &
                      + .500_r8*rxt(k,445)*y(k,138) + rxt(k,396)*y(k,170)
         mat(k,907) = -(rxt(k,411)*y(k,250))
         mat(k,2096) = -rxt(k,411)*y(k,156)
         mat(k,1124) = rxt(k,400)*y(k,155)
         mat(k,628) = .500_r8*rxt(k,471)*y(k,250)
         mat(k,492) = rxt(k,478)*y(k,250)
         mat(k,455) = rxt(k,482)*y(k,250)
         mat(k,1172) = rxt(k,483)*y(k,250)
         mat(k,2614) = rxt(k,400)*y(k,33)
         mat(k,2096) = mat(k,2096) + .500_r8*rxt(k,471)*y(k,129) + rxt(k,478)*y(k,130) &
                      + rxt(k,482)*y(k,144) + rxt(k,483)*y(k,145)
         mat(k,460) = -(rxt(k,543)*y(k,250))
         mat(k,2043) = -rxt(k,543)*y(k,157)
         mat(k,2241) = rxt(k,540)*y(k,248)
         mat(k,1143) = rxt(k,540)*y(k,235)
         mat(k,2506) = -(rxt(k,159)*y(k,164) + 4._r8*rxt(k,160)*y(k,163) + rxt(k,162) &
                      *y(k,93) + rxt(k,163)*y(k,95) + rxt(k,168)*y(k,235) + rxt(k,174) &
                      *y(k,250) + (rxt(k,185) + rxt(k,187)) * y(k,154) + rxt(k,190) &
                      *y(k,155) + rxt(k,195)*y(k,153) + rxt(k,221)*y(k,75) + rxt(k,223) &
                      *y(k,74) + rxt(k,226)*y(k,101) + rxt(k,229)*y(k,109) + rxt(k,256) &
                      *y(k,22) + rxt(k,257)*y(k,21) + rxt(k,259)*y(k,97) + rxt(k,261) &
                      *y(k,108) + rxt(k,270)*y(k,116) + rxt(k,292)*y(k,125) + rxt(k,350) &
                      *y(k,51) + rxt(k,562)*y(k,167))
         mat(k,2600) = -rxt(k,159)*y(k,163)
         mat(k,1528) = -rxt(k,162)*y(k,163)
         mat(k,697) = -rxt(k,163)*y(k,163)
         mat(k,2329) = -rxt(k,168)*y(k,163)
         mat(k,2157) = -rxt(k,174)*y(k,163)
         mat(k,1829) = -(rxt(k,185) + rxt(k,187)) * y(k,163)
         mat(k,2670) = -rxt(k,190)*y(k,163)
         mat(k,1929) = -rxt(k,195)*y(k,163)
         mat(k,1088) = -rxt(k,221)*y(k,163)
         mat(k,2214) = -rxt(k,223)*y(k,163)
         mat(k,1724) = -rxt(k,226)*y(k,163)
         mat(k,1748) = -rxt(k,229)*y(k,163)
         mat(k,922) = -rxt(k,256)*y(k,163)
         mat(k,2185) = -rxt(k,257)*y(k,163)
         mat(k,1514) = -rxt(k,259)*y(k,163)
         mat(k,1596) = -rxt(k,261)*y(k,163)
         mat(k,1540) = -rxt(k,270)*y(k,163)
         mat(k,2420) = -rxt(k,292)*y(k,163)
         mat(k,1700) = -rxt(k,350)*y(k,163)
         mat(k,431) = -rxt(k,562)*y(k,163)
         mat(k,2471) = rxt(k,166)*y(k,235)
         mat(k,575) = rxt(k,180)*y(k,153) + rxt(k,181)*y(k,154)
         mat(k,1929) = mat(k,1929) + rxt(k,180)*y(k,141)
         mat(k,1829) = mat(k,1829) + rxt(k,181)*y(k,141)
         mat(k,2600) = mat(k,2600) + 2.000_r8*rxt(k,158)*y(k,249)
         mat(k,2329) = mat(k,2329) + rxt(k,166)*y(k,92)
         mat(k,1975) = 2.000_r8*rxt(k,158)*y(k,164)
         mat(k,2157) = mat(k,2157) + 2.000_r8*rxt(k,176)*y(k,250)
         mat(k,2602) = -((rxt(k,157) + rxt(k,158)) * y(k,249) + rxt(k,159)*y(k,163) &
                      + rxt(k,169)*y(k,235) + rxt(k,170)*y(k,92) + rxt(k,175)*y(k,250) &
                      + rxt(k,186)*y(k,154) + rxt(k,194)*y(k,153) + rxt(k,212)*y(k,70) &
                      + rxt(k,246)*y(k,17) + rxt(k,281)*y(k,115) + rxt(k,293)*y(k,125) &
                      + rxt(k,372)*y(k,28) + rxt(k,401)*y(k,33) + rxt(k,432)*y(k,134) &
                      + rxt(k,446)*y(k,140) + rxt(k,479)*y(k,127) + rxt(k,517) &
                      *y(k,172) + rxt(k,534)*y(k,6) + rxt(k,537)*y(k,139) + rxt(k,566) &
                      *y(k,179) + rxt(k,572)*y(k,181))
         mat(k,1977) = -(rxt(k,157) + rxt(k,158)) * y(k,164)
         mat(k,2508) = -rxt(k,159)*y(k,164)
         mat(k,2331) = -rxt(k,169)*y(k,164)
         mat(k,2473) = -rxt(k,170)*y(k,164)
         mat(k,2159) = -rxt(k,175)*y(k,164)
         mat(k,1831) = -rxt(k,186)*y(k,164)
         mat(k,1931) = -rxt(k,194)*y(k,164)
         mat(k,2392) = -rxt(k,212)*y(k,164)
         mat(k,2536) = -rxt(k,246)*y(k,164)
         mat(k,2450) = -rxt(k,281)*y(k,164)
         mat(k,2422) = -rxt(k,293)*y(k,164)
         mat(k,658) = -rxt(k,372)*y(k,164)
         mat(k,1140) = -rxt(k,401)*y(k,164)
         mat(k,1333) = -rxt(k,432)*y(k,164)
         mat(k,1446) = -rxt(k,446)*y(k,164)
         mat(k,905) = -rxt(k,479)*y(k,164)
         mat(k,550) = -rxt(k,517)*y(k,164)
         mat(k,1021) = -rxt(k,534)*y(k,164)
         mat(k,993) = -rxt(k,537)*y(k,164)
         mat(k,610) = -rxt(k,566)*y(k,164)
         mat(k,1557) = -rxt(k,572)*y(k,164)
         mat(k,1498) = .150_r8*rxt(k,386)*y(k,235)
         mat(k,2331) = mat(k,2331) + .150_r8*rxt(k,386)*y(k,229) + .150_r8*rxt(k,437) &
                      *y(k,244)
         mat(k,1466) = .150_r8*rxt(k,437)*y(k,235)
         mat(k,535) = -(rxt(k,573)*y(k,181))
         mat(k,1543) = -rxt(k,573)*y(k,166)
         mat(k,2162) = rxt(k,248)*y(k,74)
         mat(k,2192) = rxt(k,248)*y(k,21) + 2.000_r8*rxt(k,216)*y(k,74) + rxt(k,285) &
                      *y(k,125)
         mat(k,2397) = rxt(k,285)*y(k,74)
         mat(k,424) = -(rxt(k,562)*y(k,163) + rxt(k,563)*y(k,250))
         mat(k,2476) = -rxt(k,562)*y(k,167)
         mat(k,2038) = -rxt(k,563)*y(k,167)
         mat(k,595) = -(rxt(k,296)*y(k,153) + rxt(k,304)*y(k,125) + 4._r8*rxt(k,305) &
                      *y(k,168))
         mat(k,1861) = -rxt(k,296)*y(k,168)
         mat(k,2399) = -rxt(k,304)*y(k,168)
         mat(k,2163) = rxt(k,284)*y(k,125)
         mat(k,2399) = mat(k,2399) + rxt(k,284)*y(k,21) + 2.000_r8*rxt(k,301)*y(k,125) &
                      + rxt(k,291)*y(k,155) + rxt(k,293)*y(k,164)
         mat(k,2610) = rxt(k,291)*y(k,125)
         mat(k,2547) = rxt(k,293)*y(k,125)
         mat(k,1229) = rxt(k,425)*y(k,250)
         mat(k,1848) = .100_r8*rxt(k,546)*y(k,254)
         mat(k,2019) = rxt(k,425)*y(k,111)
         mat(k,1209) = .100_r8*rxt(k,546)*y(k,153)
         mat(k,635) = -(rxt(k,396)*y(k,250))
         mat(k,2066) = -rxt(k,396)*y(k,170)
         mat(k,1791) = rxt(k,398)*y(k,229)
         mat(k,1470) = rxt(k,398)*y(k,154)
         mat(k,1777) = rxt(k,519)*y(k,220)
         mat(k,614) = rxt(k,519)*y(k,154)
         mat(k,547) = -(rxt(k,516)*y(k,154) + rxt(k,517)*y(k,164))
         mat(k,1785) = -rxt(k,516)*y(k,172)
         mat(k,2546) = -rxt(k,517)*y(k,172)
         mat(k,230) = .070_r8*rxt(k,503)*y(k,250)
         mat(k,1857) = rxt(k,501)*y(k,228)
         mat(k,202) = .060_r8*rxt(k,515)*y(k,250)
         mat(k,253) = .070_r8*rxt(k,531)*y(k,250)
         mat(k,730) = rxt(k,501)*y(k,153)
         mat(k,2055) = .070_r8*rxt(k,503)*y(k,82) + .060_r8*rxt(k,515)*y(k,173) &
                      + .070_r8*rxt(k,531)*y(k,216)
         mat(k,200) = -(rxt(k,515)*y(k,250))
         mat(k,2003) = -rxt(k,515)*y(k,173)
         mat(k,192) = .530_r8*rxt(k,492)*y(k,250)
         mat(k,2003) = mat(k,2003) + .530_r8*rxt(k,492)*y(k,7)
         mat(k,381) = -(rxt(k,518)*y(k,250))
         mat(k,2032) = -rxt(k,518)*y(k,174)
         mat(k,2238) = rxt(k,513)*y(k,251)
         mat(k,528) = rxt(k,513)*y(k,235)
         mat(k,643) = -(rxt(k,414)*y(k,250))
         mat(k,2067) = -rxt(k,414)*y(k,177)
         mat(k,2258) = rxt(k,412)*y(k,252)
         mat(k,860) = rxt(k,412)*y(k,235)
         mat(k,472) = -(rxt(k,418)*y(k,250))
         mat(k,2045) = -rxt(k,418)*y(k,178)
         mat(k,2243) = .850_r8*rxt(k,416)*y(k,253)
         mat(k,1268) = .850_r8*rxt(k,416)*y(k,235)
         mat(k,605) = -(rxt(k,566)*y(k,164) + rxt(k,569)*y(k,250))
         mat(k,2548) = -rxt(k,566)*y(k,179)
         mat(k,2062) = -rxt(k,569)*y(k,179)
         mat(k,1546) = -(rxt(k,567)*y(k,21) + rxt(k,568)*y(k,74) + rxt(k,570)*y(k,154) &
                      + rxt(k,572)*y(k,164) + rxt(k,573)*y(k,166) + rxt(k,574) &
                      *y(k,250))
         mat(k,2168) = -rxt(k,567)*y(k,181)
         mat(k,2197) = -rxt(k,568)*y(k,181)
         mat(k,1809) = -rxt(k,570)*y(k,181)
         mat(k,2582) = -rxt(k,572)*y(k,181)
         mat(k,537) = -rxt(k,573)*y(k,181)
         mat(k,2137) = -rxt(k,574)*y(k,181)
         mat(k,2486) = rxt(k,562)*y(k,167)
         mat(k,2582) = mat(k,2582) + rxt(k,566)*y(k,179)
         mat(k,428) = rxt(k,562)*y(k,163)
         mat(k,606) = rxt(k,566)*y(k,164) + rxt(k,569)*y(k,250)
         mat(k,2137) = mat(k,2137) + rxt(k,569)*y(k,179)
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
         mat(k,1103) = -(rxt(k,565)*y(k,250))
         mat(k,2110) = -rxt(k,565)*y(k,182)
         mat(k,2167) = rxt(k,556)*y(k,83) + rxt(k,567)*y(k,181)
         mat(k,2361) = rxt(k,558)*y(k,83)
         mat(k,2196) = rxt(k,568)*y(k,181)
         mat(k,1091) = rxt(k,556)*y(k,21) + rxt(k,558)*y(k,70) + rxt(k,559)*y(k,125) &
                      + rxt(k,560)*y(k,155) + (rxt(k,561)+.500_r8*rxt(k,575))*y(k,250)
         mat(k,2402) = rxt(k,559)*y(k,83)
         mat(k,1802) = rxt(k,570)*y(k,181)
         mat(k,2625) = rxt(k,560)*y(k,83)
         mat(k,2563) = rxt(k,572)*y(k,181)
         mat(k,536) = rxt(k,573)*y(k,181)
         mat(k,426) = rxt(k,563)*y(k,250)
         mat(k,1545) = rxt(k,567)*y(k,21) + rxt(k,568)*y(k,74) + rxt(k,570)*y(k,154) &
                      + rxt(k,572)*y(k,164) + rxt(k,573)*y(k,166) + rxt(k,574) &
                      *y(k,250)
         mat(k,2110) = mat(k,2110) + (rxt(k,561)+.500_r8*rxt(k,575))*y(k,83) &
                      + rxt(k,563)*y(k,167) + rxt(k,574)*y(k,181)
         mat(k,294) = -(rxt(k,576)*y(k,262))
         mat(k,2676) = -rxt(k,576)*y(k,183)
         mat(k,1102) = rxt(k,565)*y(k,250)
         mat(k,2021) = rxt(k,565)*y(k,182)
         mat(k,995) = .2202005_r8*rxt(k,629)*y(k,164)
         mat(k,967) = .0508005_r8*rxt(k,645)*y(k,164)
         mat(k,1834) = .1279005_r8*rxt(k,628)*y(k,222) + .0097005_r8*rxt(k,633) &
                      *y(k,224) + .0003005_r8*rxt(k,636)*y(k,240) &
                      + .1056005_r8*rxt(k,640)*y(k,241) + .0245005_r8*rxt(k,644) &
                      *y(k,247) + .0154005_r8*rxt(k,650)*y(k,257) &
                      + .0063005_r8*rxt(k,654)*y(k,260)
         mat(k,2539) = .2202005_r8*rxt(k,629)*y(k,6) + .0508005_r8*rxt(k,645)*y(k,139)
         mat(k,54) = .5931005_r8*rxt(k,647)*y(k,250)
         mat(k,60) = .1279005_r8*rxt(k,628)*y(k,153) + .2202005_r8*rxt(k,627)*y(k,235)
         mat(k,66) = .0097005_r8*rxt(k,633)*y(k,153) + .0023005_r8*rxt(k,632)*y(k,235)
         mat(k,2219) = .2202005_r8*rxt(k,627)*y(k,222) + .0023005_r8*rxt(k,632) &
                      *y(k,224) + .0031005_r8*rxt(k,635)*y(k,240) &
                      + .2381005_r8*rxt(k,639)*y(k,241) + .0508005_r8*rxt(k,643) &
                      *y(k,247) + .1364005_r8*rxt(k,649)*y(k,257) &
                      + .1677005_r8*rxt(k,653)*y(k,260)
         mat(k,72) = .0003005_r8*rxt(k,636)*y(k,153) + .0031005_r8*rxt(k,635)*y(k,235)
         mat(k,78) = .1056005_r8*rxt(k,640)*y(k,153) + .2381005_r8*rxt(k,639)*y(k,235)
         mat(k,86) = .0245005_r8*rxt(k,644)*y(k,153) + .0508005_r8*rxt(k,643)*y(k,235)
         mat(k,1980) = .5931005_r8*rxt(k,647)*y(k,204)
         mat(k,92) = .0154005_r8*rxt(k,650)*y(k,153) + .1364005_r8*rxt(k,649)*y(k,235)
         mat(k,98) = .0063005_r8*rxt(k,654)*y(k,153) + .1677005_r8*rxt(k,653)*y(k,235)
         mat(k,996) = .2067005_r8*rxt(k,629)*y(k,164)
         mat(k,968) = .1149005_r8*rxt(k,645)*y(k,164)
         mat(k,1835) = .1792005_r8*rxt(k,628)*y(k,222) + .0034005_r8*rxt(k,633) &
                      *y(k,224) + .0003005_r8*rxt(k,636)*y(k,240) &
                      + .1026005_r8*rxt(k,640)*y(k,241) + .0082005_r8*rxt(k,644) &
                      *y(k,247) + .0452005_r8*rxt(k,650)*y(k,257) &
                      + .0237005_r8*rxt(k,654)*y(k,260)
         mat(k,2540) = .2067005_r8*rxt(k,629)*y(k,6) + .1149005_r8*rxt(k,645)*y(k,139)
         mat(k,55) = .1534005_r8*rxt(k,647)*y(k,250)
         mat(k,61) = .1792005_r8*rxt(k,628)*y(k,153) + .2067005_r8*rxt(k,627)*y(k,235)
         mat(k,67) = .0034005_r8*rxt(k,633)*y(k,153) + .0008005_r8*rxt(k,632)*y(k,235)
         mat(k,2220) = .2067005_r8*rxt(k,627)*y(k,222) + .0008005_r8*rxt(k,632) &
                      *y(k,224) + .0035005_r8*rxt(k,635)*y(k,240) &
                      + .1308005_r8*rxt(k,639)*y(k,241) + .1149005_r8*rxt(k,643) &
                      *y(k,247) + .0101005_r8*rxt(k,649)*y(k,257) &
                      + .0174005_r8*rxt(k,653)*y(k,260)
         mat(k,73) = .0003005_r8*rxt(k,636)*y(k,153) + .0035005_r8*rxt(k,635)*y(k,235)
         mat(k,79) = .1026005_r8*rxt(k,640)*y(k,153) + .1308005_r8*rxt(k,639)*y(k,235)
         mat(k,87) = .0082005_r8*rxt(k,644)*y(k,153) + .1149005_r8*rxt(k,643)*y(k,235)
         mat(k,1981) = .1534005_r8*rxt(k,647)*y(k,204)
         mat(k,93) = .0452005_r8*rxt(k,650)*y(k,153) + .0101005_r8*rxt(k,649)*y(k,235)
         mat(k,99) = .0237005_r8*rxt(k,654)*y(k,153) + .0174005_r8*rxt(k,653)*y(k,235)
         mat(k,997) = .0653005_r8*rxt(k,629)*y(k,164)
         mat(k,969) = .0348005_r8*rxt(k,645)*y(k,164)
         mat(k,1836) = .0676005_r8*rxt(k,628)*y(k,222) + .1579005_r8*rxt(k,633) &
                      *y(k,224) + .0073005_r8*rxt(k,636)*y(k,240) &
                      + .0521005_r8*rxt(k,640)*y(k,241) + .0772005_r8*rxt(k,644) &
                      *y(k,247) + .0966005_r8*rxt(k,650)*y(k,257) &
                      + .0025005_r8*rxt(k,654)*y(k,260)
         mat(k,2541) = .0653005_r8*rxt(k,629)*y(k,6) + .0348005_r8*rxt(k,645)*y(k,139)
         mat(k,56) = .0459005_r8*rxt(k,647)*y(k,250)
         mat(k,62) = .0676005_r8*rxt(k,628)*y(k,153) + .0653005_r8*rxt(k,627)*y(k,235)
         mat(k,68) = .1579005_r8*rxt(k,633)*y(k,153) + .0843005_r8*rxt(k,632)*y(k,235)
         mat(k,2221) = .0653005_r8*rxt(k,627)*y(k,222) + .0843005_r8*rxt(k,632) &
                      *y(k,224) + .0003005_r8*rxt(k,635)*y(k,240) &
                      + .0348005_r8*rxt(k,639)*y(k,241) + .0348005_r8*rxt(k,643) &
                      *y(k,247) + .0763005_r8*rxt(k,649)*y(k,257) + .086_r8*rxt(k,653) &
                      *y(k,260)
         mat(k,74) = .0073005_r8*rxt(k,636)*y(k,153) + .0003005_r8*rxt(k,635)*y(k,235)
         mat(k,80) = .0521005_r8*rxt(k,640)*y(k,153) + .0348005_r8*rxt(k,639)*y(k,235)
         mat(k,88) = .0772005_r8*rxt(k,644)*y(k,153) + .0348005_r8*rxt(k,643)*y(k,235)
         mat(k,1982) = .0459005_r8*rxt(k,647)*y(k,204)
         mat(k,94) = .0966005_r8*rxt(k,650)*y(k,153) + .0763005_r8*rxt(k,649)*y(k,235)
         mat(k,100) = .0025005_r8*rxt(k,654)*y(k,153) + .086_r8*rxt(k,653)*y(k,235)
         mat(k,998) = .1749305_r8*rxt(k,626)*y(k,155) + .1284005_r8*rxt(k,629) &
                      *y(k,164)
         mat(k,888) = .0590245_r8*rxt(k,634)*y(k,155) + .0033005_r8*rxt(k,637) &
                      *y(k,164)
         mat(k,970) = .1749305_r8*rxt(k,642)*y(k,155) + .0554005_r8*rxt(k,645) &
                      *y(k,164)
         mat(k,1837) = .079_r8*rxt(k,628)*y(k,222) + .0059005_r8*rxt(k,633)*y(k,224) &
                      + .0057005_r8*rxt(k,636)*y(k,240) + .0143005_r8*rxt(k,640) &
                      *y(k,241) + .0332005_r8*rxt(k,644)*y(k,247) &
                      + .0073005_r8*rxt(k,650)*y(k,257) + .011_r8*rxt(k,654)*y(k,260)
         mat(k,2605) = .1749305_r8*rxt(k,626)*y(k,6) + .0590245_r8*rxt(k,634)*y(k,127) &
                      + .1749305_r8*rxt(k,642)*y(k,139)
         mat(k,2542) = .1284005_r8*rxt(k,629)*y(k,6) + .0033005_r8*rxt(k,637)*y(k,127) &
                      + .0554005_r8*rxt(k,645)*y(k,139)
         mat(k,57) = .0085005_r8*rxt(k,647)*y(k,250)
         mat(k,63) = .079_r8*rxt(k,628)*y(k,153) + .1284005_r8*rxt(k,627)*y(k,235)
         mat(k,69) = .0059005_r8*rxt(k,633)*y(k,153) + .0443005_r8*rxt(k,632)*y(k,235)
         mat(k,2222) = .1284005_r8*rxt(k,627)*y(k,222) + .0443005_r8*rxt(k,632) &
                      *y(k,224) + .0271005_r8*rxt(k,635)*y(k,240) &
                      + .0076005_r8*rxt(k,639)*y(k,241) + .0554005_r8*rxt(k,643) &
                      *y(k,247) + .2157005_r8*rxt(k,649)*y(k,257) &
                      + .0512005_r8*rxt(k,653)*y(k,260)
         mat(k,75) = .0057005_r8*rxt(k,636)*y(k,153) + .0271005_r8*rxt(k,635)*y(k,235)
         mat(k,81) = .0143005_r8*rxt(k,640)*y(k,153) + .0076005_r8*rxt(k,639)*y(k,235)
         mat(k,89) = .0332005_r8*rxt(k,644)*y(k,153) + .0554005_r8*rxt(k,643)*y(k,235)
         mat(k,1983) = .0085005_r8*rxt(k,647)*y(k,204)
         mat(k,95) = .0073005_r8*rxt(k,650)*y(k,153) + .2157005_r8*rxt(k,649)*y(k,235)
         mat(k,101) = .011_r8*rxt(k,654)*y(k,153) + .0512005_r8*rxt(k,653)*y(k,235)
         mat(k,999) = .5901905_r8*rxt(k,626)*y(k,155) + .114_r8*rxt(k,629)*y(k,164)
         mat(k,889) = .0250245_r8*rxt(k,634)*y(k,155)
         mat(k,971) = .5901905_r8*rxt(k,642)*y(k,155) + .1278005_r8*rxt(k,645) &
                      *y(k,164)
         mat(k,1838) = .1254005_r8*rxt(k,628)*y(k,222) + .0536005_r8*rxt(k,633) &
                      *y(k,224) + .0623005_r8*rxt(k,636)*y(k,240) &
                      + .0166005_r8*rxt(k,640)*y(k,241) + .130_r8*rxt(k,644)*y(k,247) &
                      + .238_r8*rxt(k,650)*y(k,257) + .1185005_r8*rxt(k,654)*y(k,260)
         mat(k,2606) = .5901905_r8*rxt(k,626)*y(k,6) + .0250245_r8*rxt(k,634)*y(k,127) &
                      + .5901905_r8*rxt(k,642)*y(k,139)
         mat(k,2543) = .114_r8*rxt(k,629)*y(k,6) + .1278005_r8*rxt(k,645)*y(k,139)
         mat(k,58) = .0128005_r8*rxt(k,647)*y(k,250)
         mat(k,64) = .1254005_r8*rxt(k,628)*y(k,153) + .114_r8*rxt(k,627)*y(k,235)
         mat(k,70) = .0536005_r8*rxt(k,633)*y(k,153) + .1621005_r8*rxt(k,632)*y(k,235)
         mat(k,2223) = .114_r8*rxt(k,627)*y(k,222) + .1621005_r8*rxt(k,632)*y(k,224) &
                      + .0474005_r8*rxt(k,635)*y(k,240) + .0113005_r8*rxt(k,639) &
                      *y(k,241) + .1278005_r8*rxt(k,643)*y(k,247) &
                      + .0738005_r8*rxt(k,649)*y(k,257) + .1598005_r8*rxt(k,653) &
                      *y(k,260)
         mat(k,76) = .0623005_r8*rxt(k,636)*y(k,153) + .0474005_r8*rxt(k,635)*y(k,235)
         mat(k,82) = .0166005_r8*rxt(k,640)*y(k,153) + .0113005_r8*rxt(k,639)*y(k,235)
         mat(k,90) = .130_r8*rxt(k,644)*y(k,153) + .1278005_r8*rxt(k,643)*y(k,235)
         mat(k,1984) = .0128005_r8*rxt(k,647)*y(k,204)
         mat(k,96) = .238_r8*rxt(k,650)*y(k,153) + .0738005_r8*rxt(k,649)*y(k,235)
         mat(k,102) = .1185005_r8*rxt(k,654)*y(k,153) + .1598005_r8*rxt(k,653) &
                      *y(k,235)
         mat(k,59) = -(rxt(k,647)*y(k,250))
         mat(k,1985) = -rxt(k,647)*y(k,204)
         mat(k,223) = .100_r8*rxt(k,523)*y(k,250)
         mat(k,243) = .230_r8*rxt(k,525)*y(k,250)
         mat(k,2008) = .100_r8*rxt(k,523)*y(k,212) + .230_r8*rxt(k,525)*y(k,214)
         mat(k,748) = -(rxt(k,547)*y(k,250))
         mat(k,2080) = -rxt(k,547)*y(k,206)
         mat(k,2264) = rxt(k,545)*y(k,254)
         mat(k,1210) = rxt(k,545)*y(k,235)
         mat(k,723) = -(rxt(k,548)*y(k,250))
         mat(k,2077) = -rxt(k,548)*y(k,207)
         mat(k,1868) = .200_r8*rxt(k,541)*y(k,248) + .200_r8*rxt(k,551)*y(k,255)
         mat(k,1605) = .500_r8*rxt(k,539)*y(k,248)
         mat(k,1144) = .200_r8*rxt(k,541)*y(k,153) + .500_r8*rxt(k,539)*y(k,230)
         mat(k,1053) = .200_r8*rxt(k,551)*y(k,153)
         mat(k,561) = -(rxt(k,552)*y(k,250))
         mat(k,2057) = -rxt(k,552)*y(k,208)
         mat(k,2254) = rxt(k,550)*y(k,255)
         mat(k,1052) = rxt(k,550)*y(k,235)
         mat(k,1065) = -(rxt(k,553)*y(k,155) + rxt(k,554)*y(k,250))
         mat(k,2622) = -rxt(k,553)*y(k,209)
         mat(k,2106) = -rxt(k,554)*y(k,209)
         mat(k,1009) = .330_r8*rxt(k,534)*y(k,164)
         mat(k,981) = .330_r8*rxt(k,537)*y(k,164)
         mat(k,1887) = .800_r8*rxt(k,541)*y(k,248) + .800_r8*rxt(k,551)*y(k,255)
         mat(k,2622) = mat(k,2622) + rxt(k,542)*y(k,248)
         mat(k,2561) = .330_r8*rxt(k,534)*y(k,6) + .330_r8*rxt(k,537)*y(k,139)
         mat(k,724) = rxt(k,548)*y(k,250)
         mat(k,1614) = .500_r8*rxt(k,539)*y(k,248) + rxt(k,549)*y(k,255)
         mat(k,1146) = .800_r8*rxt(k,541)*y(k,153) + rxt(k,542)*y(k,155) &
                      + .500_r8*rxt(k,539)*y(k,230)
         mat(k,2106) = mat(k,2106) + rxt(k,548)*y(k,207)
         mat(k,1056) = .800_r8*rxt(k,551)*y(k,153) + rxt(k,549)*y(k,230)
         mat(k,1161) = -(rxt(k,555)*y(k,250))
         mat(k,2114) = -rxt(k,555)*y(k,210)
         mat(k,1012) = .300_r8*rxt(k,534)*y(k,164)
         mat(k,984) = .300_r8*rxt(k,537)*y(k,164)
         mat(k,1891) = .900_r8*rxt(k,546)*y(k,254)
         mat(k,2566) = .300_r8*rxt(k,534)*y(k,6) + .300_r8*rxt(k,537)*y(k,139)
         mat(k,1618) = rxt(k,544)*y(k,254)
         mat(k,1214) = .900_r8*rxt(k,546)*y(k,153) + rxt(k,544)*y(k,230)
         mat(k,701) = -(rxt(k,522)*y(k,250))
         mat(k,2074) = -rxt(k,522)*y(k,211)
         mat(k,2261) = rxt(k,520)*y(k,256)
         mat(k,814) = rxt(k,520)*y(k,235)
         mat(k,221) = -(rxt(k,523)*y(k,250))
         mat(k,2006) = -rxt(k,523)*y(k,212)
         mat(k,237) = -(rxt(k,489)*y(k,250))
         mat(k,2009) = -rxt(k,489)*y(k,213)
         mat(k,2232) = rxt(k,486)*y(k,258)
         mat(k,1281) = rxt(k,486)*y(k,235)
         mat(k,244) = -(rxt(k,525)*y(k,250))
         mat(k,2010) = -rxt(k,525)*y(k,214)
         mat(k,786) = -(rxt(k,528)*y(k,250))
         mat(k,2084) = -rxt(k,528)*y(k,215)
         mat(k,2268) = rxt(k,526)*y(k,259)
         mat(k,839) = rxt(k,526)*y(k,235)
         mat(k,252) = -(rxt(k,531)*y(k,250))
         mat(k,2011) = -rxt(k,531)*y(k,216)
         mat(k,245) = .150_r8*rxt(k,525)*y(k,250)
         mat(k,2011) = mat(k,2011) + .150_r8*rxt(k,525)*y(k,214)
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
         mat(k,502) = -(rxt(k,532)*y(k,250))
         mat(k,2049) = -rxt(k,532)*y(k,217)
         mat(k,2247) = rxt(k,529)*y(k,261)
         mat(k,585) = rxt(k,529)*y(k,235)
         mat(k,615) = -(rxt(k,490)*y(k,235) + rxt(k,491)*y(k,153) + rxt(k,519) &
                      *y(k,154))
         mat(k,2257) = -rxt(k,490)*y(k,220)
         mat(k,1863) = -rxt(k,491)*y(k,220)
         mat(k,1789) = -rxt(k,519)*y(k,220)
         mat(k,283) = rxt(k,496)*y(k,250)
         mat(k,2063) = rxt(k,496)*y(k,24)
         mat(k,1028) = -(rxt(k,451)*y(k,235) + (rxt(k,452) + rxt(k,453)) * y(k,153))
         mat(k,2281) = -rxt(k,451)*y(k,221)
         mat(k,1883) = -(rxt(k,452) + rxt(k,453)) * y(k,221)
         mat(k,762) = rxt(k,454)*y(k,250)
         mat(k,280) = rxt(k,455)*y(k,250)
         mat(k,2102) = rxt(k,454)*y(k,2) + rxt(k,455)*y(k,15)
         mat(k,65) = -(rxt(k,627)*y(k,235) + rxt(k,628)*y(k,153))
         mat(k,2224) = -rxt(k,627)*y(k,222)
         mat(k,1839) = -rxt(k,628)*y(k,222)
         mat(k,1000) = rxt(k,630)*y(k,250)
         mat(k,1986) = rxt(k,630)*y(k,6)
         mat(k,578) = -(rxt(k,493)*y(k,235) + rxt(k,494)*y(k,153))
         mat(k,2255) = -rxt(k,493)*y(k,223)
         mat(k,1859) = -rxt(k,494)*y(k,223)
         mat(k,193) = .350_r8*rxt(k,492)*y(k,250)
         mat(k,468) = rxt(k,495)*y(k,250)
         mat(k,2059) = .350_r8*rxt(k,492)*y(k,7) + rxt(k,495)*y(k,8)
         mat(k,71) = -(rxt(k,632)*y(k,235) + rxt(k,633)*y(k,153))
         mat(k,2225) = -rxt(k,632)*y(k,224)
         mat(k,1840) = -rxt(k,633)*y(k,224)
         mat(k,189) = rxt(k,631)*y(k,250)
         mat(k,1987) = rxt(k,631)*y(k,7)
         mat(k,510) = -(rxt(k,497)*y(k,235) + rxt(k,499)*y(k,153))
         mat(k,2248) = -rxt(k,497)*y(k,225)
         mat(k,1853) = -rxt(k,499)*y(k,225)
         mat(k,377) = rxt(k,498)*y(k,250)
         mat(k,224) = .070_r8*rxt(k,523)*y(k,250)
         mat(k,246) = .060_r8*rxt(k,525)*y(k,250)
         mat(k,2050) = rxt(k,498)*y(k,25) + .070_r8*rxt(k,523)*y(k,212) &
                      + .060_r8*rxt(k,525)*y(k,214)
         mat(k,880) = -(4._r8*rxt(k,373)*y(k,226) + rxt(k,374)*y(k,230) + rxt(k,375) &
                      *y(k,235) + rxt(k,376)*y(k,153))
         mat(k,1609) = -rxt(k,374)*y(k,226)
         mat(k,2277) = -rxt(k,375)*y(k,226)
         mat(k,1879) = -rxt(k,376)*y(k,226)
         mat(k,387) = .500_r8*rxt(k,378)*y(k,250)
         mat(k,336) = rxt(k,379)*y(k,70) + rxt(k,380)*y(k,250)
         mat(k,2353) = rxt(k,379)*y(k,32)
         mat(k,2094) = .500_r8*rxt(k,378)*y(k,31) + rxt(k,380)*y(k,32)
         mat(k,926) = -(rxt(k,402)*y(k,230) + rxt(k,403)*y(k,235) + rxt(k,404) &
                      *y(k,153))
         mat(k,1611) = -rxt(k,402)*y(k,227)
         mat(k,2279) = -rxt(k,403)*y(k,227)
         mat(k,1881) = -rxt(k,404)*y(k,227)
         mat(k,479) = rxt(k,405)*y(k,250)
         mat(k,342) = rxt(k,409)*y(k,70) + rxt(k,406)*y(k,250)
         mat(k,2355) = rxt(k,409)*y(k,35)
         mat(k,2097) = rxt(k,405)*y(k,34) + rxt(k,406)*y(k,35)
         mat(k,731) = -(rxt(k,500)*y(k,235) + rxt(k,501)*y(k,153))
         mat(k,2263) = -rxt(k,500)*y(k,228)
         mat(k,1869) = -rxt(k,501)*y(k,228)
         mat(k,304) = rxt(k,502)*y(k,250)
         mat(k,1869) = mat(k,1869) + rxt(k,491)*y(k,220)
         mat(k,2550) = rxt(k,517)*y(k,172)
         mat(k,548) = rxt(k,517)*y(k,164)
         mat(k,616) = rxt(k,491)*y(k,153) + .400_r8*rxt(k,490)*y(k,235)
         mat(k,2263) = mat(k,2263) + .400_r8*rxt(k,490)*y(k,220)
         mat(k,2078) = rxt(k,502)*y(k,36)
         mat(k,1488) = -(4._r8*rxt(k,384)*y(k,229) + rxt(k,385)*y(k,230) + rxt(k,386) &
                      *y(k,235) + rxt(k,387)*y(k,153) + rxt(k,398)*y(k,154) + rxt(k,426) &
                      *y(k,242) + rxt(k,459)*y(k,237) + rxt(k,464)*y(k,238) + rxt(k,473) &
                      *y(k,239) + rxt(k,484)*y(k,258))
         mat(k,1635) = -rxt(k,385)*y(k,229)
         mat(k,2306) = -rxt(k,386)*y(k,229)
         mat(k,1909) = -rxt(k,387)*y(k,229)
         mat(k,1807) = -rxt(k,398)*y(k,229)
         mat(k,1415) = -rxt(k,426)*y(k,229)
         mat(k,1360) = -rxt(k,459)*y(k,229)
         mat(k,1393) = -rxt(k,464)*y(k,229)
         mat(k,1312) = -rxt(k,473)*y(k,229)
         mat(k,1290) = -rxt(k,484)*y(k,229)
         mat(k,1016) = .060_r8*rxt(k,534)*y(k,164)
         mat(k,1189) = rxt(k,381)*y(k,155) + rxt(k,382)*y(k,250)
         mat(k,1337) = rxt(k,407)*y(k,155) + rxt(k,408)*y(k,250)
         mat(k,669) = .500_r8*rxt(k,389)*y(k,250)
         mat(k,900) = .080_r8*rxt(k,479)*y(k,164)
         mat(k,1328) = .100_r8*rxt(k,432)*y(k,164)
         mat(k,988) = .060_r8*rxt(k,537)*y(k,164)
         mat(k,1436) = .280_r8*rxt(k,446)*y(k,164)
         mat(k,1909) = mat(k,1909) + .530_r8*rxt(k,430)*y(k,242) + rxt(k,439)*y(k,244) &
                      + rxt(k,442)*y(k,246) + rxt(k,417)*y(k,253)
         mat(k,2648) = rxt(k,381)*y(k,54) + rxt(k,407)*y(k,58) + .530_r8*rxt(k,429) &
                      *y(k,242) + rxt(k,440)*y(k,244)
         mat(k,2581) = .060_r8*rxt(k,534)*y(k,6) + .080_r8*rxt(k,479)*y(k,127) &
                      + .100_r8*rxt(k,432)*y(k,134) + .060_r8*rxt(k,537)*y(k,139) &
                      + .280_r8*rxt(k,446)*y(k,140)
         mat(k,1164) = .650_r8*rxt(k,555)*y(k,250)
         mat(k,1488) = mat(k,1488) + .530_r8*rxt(k,426)*y(k,242)
         mat(k,1635) = mat(k,1635) + .260_r8*rxt(k,427)*y(k,242) + rxt(k,436)*y(k,244) &
                      + .300_r8*rxt(k,415)*y(k,253)
         mat(k,2306) = mat(k,2306) + .450_r8*rxt(k,437)*y(k,244) + .200_r8*rxt(k,441) &
                      *y(k,246) + .150_r8*rxt(k,416)*y(k,253)
         mat(k,1415) = mat(k,1415) + .530_r8*rxt(k,430)*y(k,153) + .530_r8*rxt(k,429) &
                      *y(k,155) + .530_r8*rxt(k,426)*y(k,229) + .260_r8*rxt(k,427) &
                      *y(k,230)
         mat(k,1457) = rxt(k,439)*y(k,153) + rxt(k,440)*y(k,155) + rxt(k,436)*y(k,230) &
                      + .450_r8*rxt(k,437)*y(k,235) + 4.000_r8*rxt(k,438)*y(k,244)
         mat(k,772) = rxt(k,442)*y(k,153) + .200_r8*rxt(k,441)*y(k,235)
         mat(k,2133) = rxt(k,382)*y(k,54) + rxt(k,408)*y(k,58) + .500_r8*rxt(k,389) &
                      *y(k,60) + .650_r8*rxt(k,555)*y(k,210)
         mat(k,1273) = rxt(k,417)*y(k,153) + .300_r8*rxt(k,415)*y(k,230) &
                      + .150_r8*rxt(k,416)*y(k,235)
         mat(k,1637) = -(rxt(k,213)*y(k,74) + (rxt(k,332) + rxt(k,333)) * y(k,68) &
                      + (4._r8*rxt(k,352) + 4._r8*rxt(k,353)) * y(k,230) + rxt(k,354) &
                      *y(k,235) + rxt(k,355)*y(k,153) + rxt(k,374)*y(k,226) + rxt(k,385) &
                      *y(k,229) + rxt(k,402)*y(k,227) + rxt(k,415)*y(k,253) + rxt(k,427) &
                      *y(k,242) + rxt(k,436)*y(k,244) + rxt(k,460)*y(k,237) + rxt(k,465) &
                      *y(k,238) + rxt(k,474)*y(k,239) + rxt(k,485)*y(k,258) + rxt(k,539) &
                      *y(k,248) + rxt(k,544)*y(k,254) + rxt(k,549)*y(k,255))
         mat(k,2198) = -rxt(k,213)*y(k,230)
         mat(k,1112) = -(rxt(k,332) + rxt(k,333)) * y(k,230)
         mat(k,2312) = -rxt(k,354)*y(k,230)
         mat(k,1912) = -rxt(k,355)*y(k,230)
         mat(k,882) = -rxt(k,374)*y(k,230)
         mat(k,1490) = -rxt(k,385)*y(k,230)
         mat(k,929) = -rxt(k,402)*y(k,230)
         mat(k,1274) = -rxt(k,415)*y(k,230)
         mat(k,1416) = -rxt(k,427)*y(k,230)
         mat(k,1458) = -rxt(k,436)*y(k,230)
         mat(k,1361) = -rxt(k,460)*y(k,230)
         mat(k,1394) = -rxt(k,465)*y(k,230)
         mat(k,1313) = -rxt(k,474)*y(k,230)
         mat(k,1291) = -rxt(k,485)*y(k,230)
         mat(k,1151) = -rxt(k,539)*y(k,230)
         mat(k,1220) = -rxt(k,544)*y(k,230)
         mat(k,1058) = -rxt(k,549)*y(k,230)
         mat(k,1133) = .280_r8*rxt(k,401)*y(k,164)
         mat(k,779) = rxt(k,388)*y(k,250)
         mat(k,497) = .700_r8*rxt(k,357)*y(k,250)
         mat(k,1565) = rxt(k,205)*y(k,70) + rxt(k,306)*y(k,89) + rxt(k,364)*y(k,249) &
                      + rxt(k,358)*y(k,250)
         mat(k,2373) = rxt(k,205)*y(k,64)
         mat(k,952) = rxt(k,306)*y(k,64)
         mat(k,901) = .050_r8*rxt(k,479)*y(k,164)
         mat(k,1912) = mat(k,1912) + rxt(k,387)*y(k,229) + .830_r8*rxt(k,505)*y(k,231) &
                      + .170_r8*rxt(k,511)*y(k,245)
         mat(k,2584) = .280_r8*rxt(k,401)*y(k,33) + .050_r8*rxt(k,479)*y(k,127)
         mat(k,1490) = mat(k,1490) + rxt(k,387)*y(k,153) + 4.000_r8*rxt(k,384) &
                      *y(k,229) + .900_r8*rxt(k,385)*y(k,230) + .450_r8*rxt(k,386) &
                      *y(k,235) + rxt(k,459)*y(k,237) + rxt(k,464)*y(k,238) &
                      + rxt(k,473)*y(k,239) + rxt(k,426)*y(k,242) + rxt(k,435) &
                      *y(k,244) + rxt(k,484)*y(k,258)
         mat(k,1637) = mat(k,1637) + .900_r8*rxt(k,385)*y(k,229)
         mat(k,855) = .830_r8*rxt(k,505)*y(k,153) + .330_r8*rxt(k,504)*y(k,235)
         mat(k,2312) = mat(k,2312) + .450_r8*rxt(k,386)*y(k,229) + .330_r8*rxt(k,504) &
                      *y(k,231) + .070_r8*rxt(k,510)*y(k,245)
         mat(k,1361) = mat(k,1361) + rxt(k,459)*y(k,229)
         mat(k,1394) = mat(k,1394) + rxt(k,464)*y(k,229)
         mat(k,1313) = mat(k,1313) + rxt(k,473)*y(k,229)
         mat(k,1416) = mat(k,1416) + rxt(k,426)*y(k,229)
         mat(k,1458) = mat(k,1458) + rxt(k,435)*y(k,229)
         mat(k,962) = .170_r8*rxt(k,511)*y(k,153) + .070_r8*rxt(k,510)*y(k,235)
         mat(k,1958) = rxt(k,364)*y(k,64)
         mat(k,2140) = rxt(k,388)*y(k,59) + .700_r8*rxt(k,357)*y(k,63) + rxt(k,358) &
                      *y(k,64)
         mat(k,1291) = mat(k,1291) + rxt(k,484)*y(k,229)
         mat(k,852) = -(rxt(k,504)*y(k,235) + rxt(k,505)*y(k,153) + rxt(k,506) &
                      *y(k,154))
         mat(k,2274) = -rxt(k,504)*y(k,231)
         mat(k,1876) = -rxt(k,505)*y(k,231)
         mat(k,1795) = -rxt(k,506)*y(k,231)
         mat(k,659) = -((rxt(k,423) + rxt(k,424)) * y(k,153))
         mat(k,1865) = -(rxt(k,423) + rxt(k,424)) * y(k,232)
         mat(k,414) = rxt(k,422)*y(k,250)
         mat(k,2069) = rxt(k,422)*y(k,16)
         mat(k,1849) = .750_r8*rxt(k,391)*y(k,234)
         mat(k,798) = .750_r8*rxt(k,391)*y(k,153)
         mat(k,799) = -(rxt(k,390)*y(k,235) + rxt(k,391)*y(k,153))
         mat(k,2269) = -rxt(k,390)*y(k,234)
         mat(k,1872) = -rxt(k,391)*y(k,234)
         mat(k,652) = rxt(k,397)*y(k,250)
         mat(k,2085) = rxt(k,397)*y(k,28)
         mat(k,2324) = -((rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,92) + rxt(k,168) &
                      *y(k,163) + rxt(k,169)*y(k,164) + rxt(k,173)*y(k,250) &
                      + 4._r8*rxt(k,178)*y(k,235) + rxt(k,188)*y(k,155) + rxt(k,193) &
                      *y(k,153) + rxt(k,198)*y(k,154) + (rxt(k,208) + rxt(k,209) &
                      ) * y(k,70) + rxt(k,217)*y(k,74) + rxt(k,244)*y(k,17) + rxt(k,251) &
                      *y(k,21) + rxt(k,273)*y(k,115) + rxt(k,288)*y(k,125) + rxt(k,334) &
                      *y(k,68) + rxt(k,348)*y(k,51) + rxt(k,354)*y(k,230) + rxt(k,361) &
                      *y(k,236) + rxt(k,375)*y(k,226) + rxt(k,386)*y(k,229) + rxt(k,390) &
                      *y(k,234) + rxt(k,403)*y(k,227) + rxt(k,412)*y(k,252) + rxt(k,416) &
                      *y(k,253) + rxt(k,428)*y(k,242) + rxt(k,437)*y(k,244) + rxt(k,441) &
                      *y(k,246) + rxt(k,451)*y(k,221) + rxt(k,461)*y(k,237) + rxt(k,466) &
                      *y(k,238) + rxt(k,475)*y(k,239) + rxt(k,486)*y(k,258) + rxt(k,490) &
                      *y(k,220) + rxt(k,493)*y(k,223) + rxt(k,497)*y(k,225) + rxt(k,500) &
                      *y(k,228) + rxt(k,504)*y(k,231) + rxt(k,507)*y(k,243) + rxt(k,510) &
                      *y(k,245) + rxt(k,513)*y(k,251) + rxt(k,520)*y(k,256) + rxt(k,526) &
                      *y(k,259) + rxt(k,529)*y(k,261) + rxt(k,540)*y(k,248) + rxt(k,545) &
                      *y(k,254) + rxt(k,550)*y(k,255))
         mat(k,2466) = -(rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,235)
         mat(k,2501) = -rxt(k,168)*y(k,235)
         mat(k,2595) = -rxt(k,169)*y(k,235)
         mat(k,2152) = -rxt(k,173)*y(k,235)
         mat(k,2665) = -rxt(k,188)*y(k,235)
         mat(k,1924) = -rxt(k,193)*y(k,235)
         mat(k,1824) = -rxt(k,198)*y(k,235)
         mat(k,2385) = -(rxt(k,208) + rxt(k,209)) * y(k,235)
         mat(k,2209) = -rxt(k,217)*y(k,235)
         mat(k,2529) = -rxt(k,244)*y(k,235)
         mat(k,2180) = -rxt(k,251)*y(k,235)
         mat(k,2443) = -rxt(k,273)*y(k,235)
         mat(k,2415) = -rxt(k,288)*y(k,235)
         mat(k,1117) = -rxt(k,334)*y(k,235)
         mat(k,1695) = -rxt(k,348)*y(k,235)
         mat(k,1647) = -rxt(k,354)*y(k,235)
         mat(k,527) = -rxt(k,361)*y(k,235)
         mat(k,887) = -rxt(k,375)*y(k,235)
         mat(k,1496) = -rxt(k,386)*y(k,235)
         mat(k,805) = -rxt(k,390)*y(k,235)
         mat(k,934) = -rxt(k,403)*y(k,235)
         mat(k,868) = -rxt(k,412)*y(k,235)
         mat(k,1279) = -rxt(k,416)*y(k,235)
         mat(k,1422) = -rxt(k,428)*y(k,235)
         mat(k,1464) = -rxt(k,437)*y(k,235)
         mat(k,776) = -rxt(k,441)*y(k,235)
         mat(k,1037) = -rxt(k,451)*y(k,235)
         mat(k,1367) = -rxt(k,461)*y(k,235)
         mat(k,1400) = -rxt(k,466)*y(k,235)
         mat(k,1319) = -rxt(k,475)*y(k,235)
         mat(k,1296) = -rxt(k,486)*y(k,235)
         mat(k,620) = -rxt(k,490)*y(k,235)
         mat(k,584) = -rxt(k,493)*y(k,235)
         mat(k,515) = -rxt(k,497)*y(k,235)
         mat(k,735) = -rxt(k,500)*y(k,235)
         mat(k,859) = -rxt(k,504)*y(k,235)
         mat(k,811) = -rxt(k,507)*y(k,235)
         mat(k,966) = -rxt(k,510)*y(k,235)
         mat(k,534) = -rxt(k,513)*y(k,235)
         mat(k,826) = -rxt(k,520)*y(k,235)
         mat(k,851) = -rxt(k,526)*y(k,235)
         mat(k,592) = -rxt(k,529)*y(k,235)
         mat(k,1157) = -rxt(k,540)*y(k,235)
         mat(k,1226) = -rxt(k,545)*y(k,235)
         mat(k,1064) = -rxt(k,550)*y(k,235)
         mat(k,1019) = .570_r8*rxt(k,534)*y(k,164)
         mat(k,195) = .650_r8*rxt(k,492)*y(k,250)
         mat(k,2529) = mat(k,2529) + rxt(k,243)*y(k,51)
         mat(k,2180) = mat(k,2180) + rxt(k,258)*y(k,250)
         mat(k,331) = .350_r8*rxt(k,370)*y(k,250)
         mat(k,657) = .130_r8*rxt(k,372)*y(k,164)
         mat(k,301) = rxt(k,377)*y(k,250)
         mat(k,1138) = .280_r8*rxt(k,401)*y(k,164)
         mat(k,1695) = mat(k,1695) + rxt(k,243)*y(k,17) + rxt(k,204)*y(k,70) &
                      + rxt(k,349)*y(k,155) + rxt(k,350)*y(k,163)
         mat(k,719) = rxt(k,321)*y(k,70) + rxt(k,322)*y(k,250)
         mat(k,445) = rxt(k,324)*y(k,70) + rxt(k,325)*y(k,250)
         mat(k,117) = rxt(k,383)*y(k,250)
         mat(k,436) = rxt(k,326)*y(k,70) + rxt(k,327)*y(k,250)
         mat(k,873) = rxt(k,356)*y(k,250)
         mat(k,1571) = rxt(k,365)*y(k,249)
         mat(k,1117) = mat(k,1117) + rxt(k,336)*y(k,153) + rxt(k,337)*y(k,155) + ( &
                      + 2.000_r8*rxt(k,332)+rxt(k,333))*y(k,230)
         mat(k,2385) = mat(k,2385) + rxt(k,204)*y(k,51) + rxt(k,321)*y(k,52) &
                      + rxt(k,324)*y(k,55) + rxt(k,326)*y(k,61) + rxt(k,207)*y(k,95)
         mat(k,2209) = mat(k,2209) + rxt(k,213)*y(k,230) + rxt(k,224)*y(k,250)
         mat(k,1207) = rxt(k,368)*y(k,250)
         mat(k,232) = .730_r8*rxt(k,503)*y(k,250)
         mat(k,1096) = .500_r8*rxt(k,575)*y(k,250)
         mat(k,1202) = rxt(k,394)*y(k,250)
         mat(k,1050) = rxt(k,395)*y(k,250)
         mat(k,695) = rxt(k,207)*y(k,70) + rxt(k,163)*y(k,163) + rxt(k,172)*y(k,250)
         mat(k,208) = rxt(k,359)*y(k,250)
         mat(k,1042) = rxt(k,360)*y(k,250)
         mat(k,1243) = rxt(k,425)*y(k,250)
         mat(k,1252) = rxt(k,410)*y(k,250)
         mat(k,2415) = mat(k,2415) + rxt(k,294)*y(k,250)
         mat(k,904) = .370_r8*rxt(k,479)*y(k,164)
         mat(k,682) = .300_r8*rxt(k,470)*y(k,250)
         mat(k,634) = rxt(k,471)*y(k,250)
         mat(k,495) = rxt(k,478)*y(k,250)
         mat(k,1331) = .140_r8*rxt(k,432)*y(k,164)
         mat(k,357) = .200_r8*rxt(k,434)*y(k,250)
         mat(k,690) = .500_r8*rxt(k,445)*y(k,250)
         mat(k,991) = .570_r8*rxt(k,537)*y(k,164)
         mat(k,1444) = .280_r8*rxt(k,446)*y(k,164)
         mat(k,459) = rxt(k,482)*y(k,250)
         mat(k,1185) = rxt(k,483)*y(k,250)
         mat(k,1924) = mat(k,1924) + rxt(k,336)*y(k,68) + rxt(k,452)*y(k,221) &
                      + rxt(k,494)*y(k,223) + rxt(k,499)*y(k,225) + rxt(k,376) &
                      *y(k,226) + rxt(k,404)*y(k,227) + rxt(k,355)*y(k,230) &
                      + .170_r8*rxt(k,505)*y(k,231) + rxt(k,423)*y(k,232) &
                      + .250_r8*rxt(k,391)*y(k,234) + rxt(k,363)*y(k,236) &
                      + .920_r8*rxt(k,462)*y(k,237) + .920_r8*rxt(k,468)*y(k,238) &
                      + rxt(k,476)*y(k,239) + .470_r8*rxt(k,430)*y(k,242) &
                      + .400_r8*rxt(k,508)*y(k,243) + .830_r8*rxt(k,511)*y(k,245) &
                      + rxt(k,514)*y(k,251) + rxt(k,413)*y(k,252) + .900_r8*rxt(k,546) &
                      *y(k,254) + .800_r8*rxt(k,551)*y(k,255) + rxt(k,521)*y(k,256) &
                      + rxt(k,487)*y(k,258) + rxt(k,527)*y(k,259) + rxt(k,530) &
                      *y(k,261)
         mat(k,2665) = mat(k,2665) + rxt(k,349)*y(k,51) + rxt(k,337)*y(k,68) &
                      + rxt(k,463)*y(k,237) + rxt(k,469)*y(k,238) + rxt(k,477) &
                      *y(k,239) + .470_r8*rxt(k,429)*y(k,242) + rxt(k,191)*y(k,250) &
                      + rxt(k,488)*y(k,258)
         mat(k,2501) = mat(k,2501) + rxt(k,350)*y(k,51) + rxt(k,163)*y(k,95)
         mat(k,2595) = mat(k,2595) + .570_r8*rxt(k,534)*y(k,6) + .130_r8*rxt(k,372) &
                      *y(k,28) + .280_r8*rxt(k,401)*y(k,33) + .370_r8*rxt(k,479) &
                      *y(k,127) + .140_r8*rxt(k,432)*y(k,134) + .570_r8*rxt(k,537) &
                      *y(k,139) + .280_r8*rxt(k,446)*y(k,140) + rxt(k,175)*y(k,250)
         mat(k,204) = .800_r8*rxt(k,515)*y(k,250)
         mat(k,1106) = rxt(k,565)*y(k,250)
         mat(k,1168) = .200_r8*rxt(k,555)*y(k,250)
         mat(k,227) = .280_r8*rxt(k,523)*y(k,250)
         mat(k,251) = .380_r8*rxt(k,525)*y(k,250)
         mat(k,256) = .630_r8*rxt(k,531)*y(k,250)
         mat(k,1037) = mat(k,1037) + rxt(k,452)*y(k,153)
         mat(k,584) = mat(k,584) + rxt(k,494)*y(k,153)
         mat(k,515) = mat(k,515) + rxt(k,499)*y(k,153)
         mat(k,887) = mat(k,887) + rxt(k,376)*y(k,153) + 2.400_r8*rxt(k,373)*y(k,226) &
                      + rxt(k,374)*y(k,230)
         mat(k,934) = mat(k,934) + rxt(k,404)*y(k,153) + rxt(k,402)*y(k,230)
         mat(k,1496) = mat(k,1496) + .900_r8*rxt(k,385)*y(k,230) + rxt(k,459)*y(k,237) &
                      + rxt(k,464)*y(k,238) + rxt(k,473)*y(k,239) + .470_r8*rxt(k,426) &
                      *y(k,242) + rxt(k,484)*y(k,258)
         mat(k,1647) = mat(k,1647) + (2.000_r8*rxt(k,332)+rxt(k,333))*y(k,68) &
                      + rxt(k,213)*y(k,74) + rxt(k,355)*y(k,153) + rxt(k,374)*y(k,226) &
                      + rxt(k,402)*y(k,227) + .900_r8*rxt(k,385)*y(k,229) &
                      + 4.000_r8*rxt(k,352)*y(k,230) + rxt(k,460)*y(k,237) &
                      + rxt(k,465)*y(k,238) + 1.200_r8*rxt(k,474)*y(k,239) &
                      + .730_r8*rxt(k,427)*y(k,242) + rxt(k,436)*y(k,244) &
                      + .500_r8*rxt(k,539)*y(k,248) + .300_r8*rxt(k,415)*y(k,253) &
                      + rxt(k,544)*y(k,254) + rxt(k,549)*y(k,255) + .800_r8*rxt(k,485) &
                      *y(k,258)
         mat(k,859) = mat(k,859) + .170_r8*rxt(k,505)*y(k,153) + .070_r8*rxt(k,504) &
                      *y(k,235)
         mat(k,666) = rxt(k,423)*y(k,153)
         mat(k,805) = mat(k,805) + .250_r8*rxt(k,391)*y(k,153)
         mat(k,2324) = mat(k,2324) + .070_r8*rxt(k,504)*y(k,231) + .160_r8*rxt(k,507) &
                      *y(k,243) + .330_r8*rxt(k,510)*y(k,245)
         mat(k,527) = mat(k,527) + rxt(k,363)*y(k,153)
         mat(k,1367) = mat(k,1367) + .920_r8*rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155) &
                      + rxt(k,459)*y(k,229) + rxt(k,460)*y(k,230)
         mat(k,1400) = mat(k,1400) + .920_r8*rxt(k,468)*y(k,153) + rxt(k,469)*y(k,155) &
                      + rxt(k,464)*y(k,229) + rxt(k,465)*y(k,230)
         mat(k,1319) = mat(k,1319) + rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155) &
                      + rxt(k,473)*y(k,229) + 1.200_r8*rxt(k,474)*y(k,230)
         mat(k,1422) = mat(k,1422) + .470_r8*rxt(k,430)*y(k,153) + .470_r8*rxt(k,429) &
                      *y(k,155) + .470_r8*rxt(k,426)*y(k,229) + .730_r8*rxt(k,427) &
                      *y(k,230)
         mat(k,811) = mat(k,811) + .400_r8*rxt(k,508)*y(k,153) + .160_r8*rxt(k,507) &
                      *y(k,235)
         mat(k,1464) = mat(k,1464) + rxt(k,436)*y(k,230)
         mat(k,966) = mat(k,966) + .830_r8*rxt(k,511)*y(k,153) + .330_r8*rxt(k,510) &
                      *y(k,235)
         mat(k,1157) = mat(k,1157) + .500_r8*rxt(k,539)*y(k,230)
         mat(k,1970) = rxt(k,365)*y(k,64)
         mat(k,2152) = mat(k,2152) + .650_r8*rxt(k,492)*y(k,7) + rxt(k,258)*y(k,21) &
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
                      *y(k,145) + rxt(k,191)*y(k,155) + rxt(k,175)*y(k,164) &
                      + .800_r8*rxt(k,515)*y(k,173) + rxt(k,565)*y(k,182) &
                      + .200_r8*rxt(k,555)*y(k,210) + .280_r8*rxt(k,523)*y(k,212) &
                      + .380_r8*rxt(k,525)*y(k,214) + .630_r8*rxt(k,531)*y(k,216)
         mat(k,534) = mat(k,534) + rxt(k,514)*y(k,153)
         mat(k,868) = mat(k,868) + rxt(k,413)*y(k,153)
         mat(k,1279) = mat(k,1279) + .300_r8*rxt(k,415)*y(k,230)
         mat(k,1226) = mat(k,1226) + .900_r8*rxt(k,546)*y(k,153) + rxt(k,544)*y(k,230)
         mat(k,1064) = mat(k,1064) + .800_r8*rxt(k,551)*y(k,153) + rxt(k,549)*y(k,230)
         mat(k,826) = mat(k,826) + rxt(k,521)*y(k,153)
         mat(k,1296) = mat(k,1296) + rxt(k,487)*y(k,153) + rxt(k,488)*y(k,155) &
                      + rxt(k,484)*y(k,229) + .800_r8*rxt(k,485)*y(k,230)
         mat(k,851) = mat(k,851) + rxt(k,527)*y(k,153)
         mat(k,592) = mat(k,592) + rxt(k,530)*y(k,153)
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
         mat(k,522) = -(rxt(k,361)*y(k,235) + rxt(k,363)*y(k,153))
         mat(k,2250) = -rxt(k,361)*y(k,236)
         mat(k,1854) = -rxt(k,363)*y(k,236)
         mat(k,1676) = rxt(k,348)*y(k,235)
         mat(k,2250) = mat(k,2250) + rxt(k,348)*y(k,51)
         mat(k,1356) = -(rxt(k,459)*y(k,229) + rxt(k,460)*y(k,230) + rxt(k,461) &
                      *y(k,235) + rxt(k,462)*y(k,153) + rxt(k,463)*y(k,155))
         mat(k,1483) = -rxt(k,459)*y(k,237)
         mat(k,1630) = -rxt(k,460)*y(k,237)
         mat(k,2301) = -rxt(k,461)*y(k,237)
         mat(k,1904) = -rxt(k,462)*y(k,237)
         mat(k,2643) = -rxt(k,463)*y(k,237)
         mat(k,897) = .600_r8*rxt(k,480)*y(k,250)
         mat(k,2128) = .600_r8*rxt(k,480)*y(k,127)
         mat(k,1389) = -(rxt(k,464)*y(k,229) + rxt(k,465)*y(k,230) + rxt(k,466) &
                      *y(k,235) + rxt(k,468)*y(k,153) + rxt(k,469)*y(k,155))
         mat(k,1484) = -rxt(k,464)*y(k,238)
         mat(k,1631) = -rxt(k,465)*y(k,238)
         mat(k,2302) = -rxt(k,466)*y(k,238)
         mat(k,1905) = -rxt(k,468)*y(k,238)
         mat(k,2644) = -rxt(k,469)*y(k,238)
         mat(k,898) = .400_r8*rxt(k,480)*y(k,250)
         mat(k,2129) = .400_r8*rxt(k,480)*y(k,127)
         mat(k,1308) = -(rxt(k,473)*y(k,229) + rxt(k,474)*y(k,230) + rxt(k,475) &
                      *y(k,235) + rxt(k,476)*y(k,153) + rxt(k,477)*y(k,155))
         mat(k,1480) = -rxt(k,473)*y(k,239)
         mat(k,1627) = -rxt(k,474)*y(k,239)
         mat(k,2298) = -rxt(k,475)*y(k,239)
         mat(k,1901) = -rxt(k,476)*y(k,239)
         mat(k,2640) = -rxt(k,477)*y(k,239)
         mat(k,895) = rxt(k,472)*y(k,155)
         mat(k,2640) = mat(k,2640) + rxt(k,472)*y(k,127)
         mat(k,77) = -(rxt(k,635)*y(k,235) + rxt(k,636)*y(k,153))
         mat(k,2226) = -rxt(k,635)*y(k,240)
         mat(k,1841) = -rxt(k,636)*y(k,240)
         mat(k,890) = rxt(k,638)*y(k,250)
         mat(k,1988) = rxt(k,638)*y(k,127)
         mat(k,83) = -(rxt(k,639)*y(k,235) + rxt(k,640)*y(k,153))
         mat(k,2227) = -rxt(k,639)*y(k,241)
         mat(k,1842) = -rxt(k,640)*y(k,241)
         mat(k,84) = rxt(k,641)*y(k,250)
         mat(k,1989) = rxt(k,641)*y(k,132)
         mat(k,1413) = -(rxt(k,426)*y(k,229) + rxt(k,427)*y(k,230) + rxt(k,428) &
                      *y(k,235) + rxt(k,429)*y(k,155) + (rxt(k,430) + rxt(k,431) &
                      ) * y(k,153))
         mat(k,1485) = -rxt(k,426)*y(k,242)
         mat(k,1632) = -rxt(k,427)*y(k,242)
         mat(k,2303) = -rxt(k,428)*y(k,242)
         mat(k,2645) = -rxt(k,429)*y(k,242)
         mat(k,1906) = -(rxt(k,430) + rxt(k,431)) * y(k,242)
         mat(k,1326) = .500_r8*rxt(k,433)*y(k,250)
         mat(k,354) = .200_r8*rxt(k,434)*y(k,250)
         mat(k,1433) = rxt(k,447)*y(k,250)
         mat(k,2130) = .500_r8*rxt(k,433)*y(k,134) + .200_r8*rxt(k,434)*y(k,135) &
                      + rxt(k,447)*y(k,140)
         mat(k,806) = -(rxt(k,507)*y(k,235) + rxt(k,508)*y(k,153) + rxt(k,509) &
                      *y(k,154))
         mat(k,2270) = -rxt(k,507)*y(k,243)
         mat(k,1873) = -rxt(k,508)*y(k,243)
         mat(k,1794) = -rxt(k,509)*y(k,243)
         mat(k,1456) = -(rxt(k,435)*y(k,229) + rxt(k,436)*y(k,230) + rxt(k,437) &
                      *y(k,235) + 4._r8*rxt(k,438)*y(k,244) + rxt(k,439)*y(k,153) &
                      + rxt(k,440)*y(k,155) + rxt(k,448)*y(k,154))
         mat(k,1487) = -rxt(k,435)*y(k,244)
         mat(k,1634) = -rxt(k,436)*y(k,244)
         mat(k,2305) = -rxt(k,437)*y(k,244)
         mat(k,1908) = -rxt(k,439)*y(k,244)
         mat(k,2647) = -rxt(k,440)*y(k,244)
         mat(k,1806) = -rxt(k,448)*y(k,244)
         mat(k,1327) = .500_r8*rxt(k,433)*y(k,250)
         mat(k,355) = .500_r8*rxt(k,434)*y(k,250)
         mat(k,2132) = .500_r8*rxt(k,433)*y(k,134) + .500_r8*rxt(k,434)*y(k,135)
         mat(k,958) = -(rxt(k,510)*y(k,235) + rxt(k,511)*y(k,153) + rxt(k,512) &
                      *y(k,154))
         mat(k,2280) = -rxt(k,510)*y(k,245)
         mat(k,1882) = -rxt(k,511)*y(k,245)
         mat(k,1799) = -rxt(k,512)*y(k,245)
         mat(k,770) = -(rxt(k,441)*y(k,235) + rxt(k,442)*y(k,153))
         mat(k,2266) = -rxt(k,441)*y(k,246)
         mat(k,1871) = -rxt(k,442)*y(k,246)
         mat(k,601) = rxt(k,443)*y(k,250)
         mat(k,359) = rxt(k,444)*y(k,250)
         mat(k,2082) = rxt(k,443)*y(k,136) + rxt(k,444)*y(k,137)
         mat(k,91) = -(rxt(k,643)*y(k,235) + rxt(k,644)*y(k,153))
         mat(k,2228) = -rxt(k,643)*y(k,247)
         mat(k,1843) = -rxt(k,644)*y(k,247)
         mat(k,972) = rxt(k,646)*y(k,250)
         mat(k,1991) = rxt(k,646)*y(k,139)
         mat(k,1147) = -(rxt(k,539)*y(k,230) + rxt(k,540)*y(k,235) + rxt(k,541) &
                      *y(k,153) + rxt(k,542)*y(k,155))
         mat(k,1617) = -rxt(k,539)*y(k,248)
         mat(k,2288) = -rxt(k,540)*y(k,248)
         mat(k,1890) = -rxt(k,541)*y(k,248)
         mat(k,2628) = -rxt(k,542)*y(k,248)
         mat(k,1011) = rxt(k,533)*y(k,155)
         mat(k,983) = rxt(k,536)*y(k,155)
         mat(k,2628) = mat(k,2628) + rxt(k,533)*y(k,6) + rxt(k,536)*y(k,139) &
                      + .500_r8*rxt(k,553)*y(k,209)
         mat(k,462) = rxt(k,543)*y(k,250)
         mat(k,1066) = .500_r8*rxt(k,553)*y(k,155)
         mat(k,2113) = rxt(k,543)*y(k,157)
         mat(k,1966) = -(rxt(k,153)*y(k,93) + rxt(k,154)*y(k,262) + (rxt(k,157) &
                      + rxt(k,158)) * y(k,164) + (rxt(k,196) + rxt(k,197)) * y(k,142) &
                      + rxt(k,231)*y(k,37) + rxt(k,232)*y(k,38) + rxt(k,233)*y(k,40) &
                      + rxt(k,234)*y(k,41) + rxt(k,235)*y(k,42) + rxt(k,236)*y(k,43) &
                      + rxt(k,237)*y(k,44) + (rxt(k,238) + rxt(k,239)) * y(k,101) &
                      + rxt(k,262)*y(k,39) + rxt(k,263)*y(k,66) + rxt(k,264)*y(k,94) &
                      + (rxt(k,265) + rxt(k,266)) * y(k,97) + rxt(k,310)*y(k,80) &
                      + rxt(k,311)*y(k,81) + rxt(k,343)*y(k,45) + rxt(k,344)*y(k,52) &
                      + rxt(k,345)*y(k,98) + rxt(k,346)*y(k,99) + rxt(k,347)*y(k,100) &
                      + (rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,64) + rxt(k,367) &
                      *y(k,102))
         mat(k,1524) = -rxt(k,153)*y(k,249)
         mat(k,2690) = -rxt(k,154)*y(k,249)
         mat(k,2591) = -(rxt(k,157) + rxt(k,158)) * y(k,249)
         mat(k,215) = -(rxt(k,196) + rxt(k,197)) * y(k,249)
         mat(k,113) = -rxt(k,231)*y(k,249)
         mat(k,171) = -rxt(k,232)*y(k,249)
         mat(k,136) = -rxt(k,233)*y(k,249)
         mat(k,182) = -rxt(k,234)*y(k,249)
         mat(k,140) = -rxt(k,235)*y(k,249)
         mat(k,187) = -rxt(k,236)*y(k,249)
         mat(k,144) = -rxt(k,237)*y(k,249)
         mat(k,1716) = -(rxt(k,238) + rxt(k,239)) * y(k,249)
         mat(k,177) = -rxt(k,262)*y(k,249)
         mat(k,450) = -rxt(k,263)*y(k,249)
         mat(k,128) = -rxt(k,264)*y(k,249)
         mat(k,1508) = -(rxt(k,265) + rxt(k,266)) * y(k,249)
         mat(k,273) = -rxt(k,310)*y(k,249)
         mat(k,264) = -rxt(k,311)*y(k,249)
         mat(k,556) = -rxt(k,343)*y(k,249)
         mat(k,717) = -rxt(k,344)*y(k,249)
         mat(k,259) = -rxt(k,345)*y(k,249)
         mat(k,268) = -rxt(k,346)*y(k,249)
         mat(k,319) = -rxt(k,347)*y(k,249)
         mat(k,1569) = -(rxt(k,364) + rxt(k,365) + rxt(k,366)) * y(k,249)
         mat(k,206) = -rxt(k,367)*y(k,249)
         mat(k,2149) = -(rxt(k,171)*y(k,93) + rxt(k,172)*y(k,95) + rxt(k,173)*y(k,235) &
                      + rxt(k,174)*y(k,163) + rxt(k,175)*y(k,164) + (4._r8*rxt(k,176) &
                      + 4._r8*rxt(k,177)) * y(k,250) + rxt(k,179)*y(k,107) + rxt(k,191) &
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
                      *y(k,90) + rxt(k,395)*y(k,91) + rxt(k,396)*y(k,170) + rxt(k,397) &
                      *y(k,28) + rxt(k,405)*y(k,34) + rxt(k,406)*y(k,35) + rxt(k,408) &
                      *y(k,58) + rxt(k,410)*y(k,113) + rxt(k,411)*y(k,156) + rxt(k,414) &
                      *y(k,177) + rxt(k,418)*y(k,178) + rxt(k,419)*y(k,33) + rxt(k,420) &
                      *y(k,57) + rxt(k,422)*y(k,16) + rxt(k,425)*y(k,111) + rxt(k,433) &
                      *y(k,134) + rxt(k,434)*y(k,135) + rxt(k,443)*y(k,136) + rxt(k,444) &
                      *y(k,137) + rxt(k,445)*y(k,138) + rxt(k,447)*y(k,140) + rxt(k,450) &
                      *y(k,1) + rxt(k,454)*y(k,2) + rxt(k,455)*y(k,15) + rxt(k,456) &
                      *y(k,112) + rxt(k,457)*y(k,114) + rxt(k,458)*y(k,122) + rxt(k,470) &
                      *y(k,128) + rxt(k,471)*y(k,129) + rxt(k,478)*y(k,130) + rxt(k,480) &
                      *y(k,127) + rxt(k,481)*y(k,131) + rxt(k,482)*y(k,144) + rxt(k,483) &
                      *y(k,145) + rxt(k,489)*y(k,213) + rxt(k,492)*y(k,7) + rxt(k,495) &
                      *y(k,8) + rxt(k,496)*y(k,24) + rxt(k,498)*y(k,25) + rxt(k,502) &
                      *y(k,36) + rxt(k,503)*y(k,82) + rxt(k,515)*y(k,173) + rxt(k,518) &
                      *y(k,174) + rxt(k,522)*y(k,211) + rxt(k,523)*y(k,212) + rxt(k,525) &
                      *y(k,214) + rxt(k,528)*y(k,215) + rxt(k,531)*y(k,216) + rxt(k,532) &
                      *y(k,217) + rxt(k,535)*y(k,6) + rxt(k,538)*y(k,139) + rxt(k,543) &
                      *y(k,157) + rxt(k,547)*y(k,206) + rxt(k,548)*y(k,207) + rxt(k,552) &
                      *y(k,208) + rxt(k,554)*y(k,209) + rxt(k,555)*y(k,210) + (rxt(k,561) &
                      + rxt(k,575)) * y(k,83) + rxt(k,563)*y(k,167) + rxt(k,565) &
                      *y(k,182) + rxt(k,569)*y(k,179) + rxt(k,574)*y(k,181) + rxt(k,594) &
                      *y(k,149))
         mat(k,1525) = -rxt(k,171)*y(k,250)
         mat(k,694) = -rxt(k,172)*y(k,250)
         mat(k,2321) = -rxt(k,173)*y(k,250)
         mat(k,2498) = -rxt(k,174)*y(k,250)
         mat(k,2592) = -rxt(k,175)*y(k,250)
         mat(k,518) = -rxt(k,179)*y(k,250)
         mat(k,2662) = -rxt(k,191)*y(k,250)
         mat(k,573) = -rxt(k,192)*y(k,250)
         mat(k,1821) = -rxt(k,200)*y(k,250)
         mat(k,1764) = -rxt(k,201)*y(k,250)
         mat(k,625) = -rxt(k,211)*y(k,250)
         mat(k,1085) = -rxt(k,222)*y(k,250)
         mat(k,2206) = -(rxt(k,224) + rxt(k,225)) * y(k,250)
         mat(k,1717) = -rxt(k,227)*y(k,250)
         mat(k,1741) = -rxt(k,230)*y(k,250)
         mat(k,545) = -rxt(k,242)*y(k,250)
         mat(k,2177) = -rxt(k,258)*y(k,250)
         mat(k,1509) = -rxt(k,260)*y(k,250)
         mat(k,1666) = -rxt(k,268)*y(k,250)
         mat(k,1536) = -rxt(k,271)*y(k,250)
         mat(k,2412) = -rxt(k,294)*y(k,250)
         mat(k,1261) = -rxt(k,295)*y(k,250)
         mat(k,218) = -rxt(k,313)*y(k,250)
         mat(k,290) = -rxt(k,315)*y(k,250)
         mat(k,557) = -rxt(k,317)*y(k,250)
         mat(k,147) = -rxt(k,318)*y(k,250)
         mat(k,350) = -rxt(k,320)*y(k,250)
         mat(k,718) = -rxt(k,322)*y(k,250)
         mat(k,151) = -rxt(k,323)*y(k,250)
         mat(k,444) = -rxt(k,325)*y(k,250)
         mat(k,435) = -rxt(k,327)*y(k,250)
         mat(k,119) = -rxt(k,328)*y(k,250)
         mat(k,451) = -rxt(k,330)*y(k,250)
         mat(k,123) = -rxt(k,331)*y(k,250)
         mat(k,402) = -rxt(k,339)*y(k,250)
         mat(k,260) = -rxt(k,340)*y(k,250)
         mat(k,269) = -rxt(k,341)*y(k,250)
         mat(k,320) = -rxt(k,342)*y(k,250)
         mat(k,1692) = -rxt(k,351)*y(k,250)
         mat(k,872) = -rxt(k,356)*y(k,250)
         mat(k,499) = -rxt(k,357)*y(k,250)
         mat(k,1570) = -rxt(k,358)*y(k,250)
         mat(k,207) = -rxt(k,359)*y(k,250)
         mat(k,1041) = -rxt(k,360)*y(k,250)
         mat(k,1206) = -rxt(k,368)*y(k,250)
         mat(k,330) = -rxt(k,370)*y(k,250)
         mat(k,300) = -rxt(k,377)*y(k,250)
         mat(k,389) = -rxt(k,378)*y(k,250)
         mat(k,338) = -rxt(k,380)*y(k,250)
         mat(k,1192) = -rxt(k,382)*y(k,250)
         mat(k,116) = -rxt(k,383)*y(k,250)
         mat(k,780) = -rxt(k,388)*y(k,250)
         mat(k,672) = -rxt(k,389)*y(k,250)
         mat(k,1201) = -rxt(k,394)*y(k,250)
         mat(k,1049) = -rxt(k,395)*y(k,250)
         mat(k,641) = -rxt(k,396)*y(k,250)
         mat(k,656) = -rxt(k,397)*y(k,250)
         mat(k,481) = -rxt(k,405)*y(k,250)
         mat(k,344) = -rxt(k,406)*y(k,250)
         mat(k,1339) = -rxt(k,408)*y(k,250)
         mat(k,1251) = -rxt(k,410)*y(k,250)
         mat(k,912) = -rxt(k,411)*y(k,250)
         mat(k,648) = -rxt(k,414)*y(k,250)
         mat(k,476) = -rxt(k,418)*y(k,250)
         mat(k,1137) = -rxt(k,419)*y(k,250)
         mat(k,1077) = -rxt(k,420)*y(k,250)
         mat(k,419) = -rxt(k,422)*y(k,250)
         mat(k,1242) = -rxt(k,425)*y(k,250)
         mat(k,1330) = -rxt(k,433)*y(k,250)
         mat(k,356) = -rxt(k,434)*y(k,250)
         mat(k,604) = -rxt(k,443)*y(k,250)
         mat(k,362) = -rxt(k,444)*y(k,250)
         mat(k,689) = -rxt(k,445)*y(k,250)
         mat(k,1443) = -rxt(k,447)*y(k,250)
         mat(k,746) = -rxt(k,450)*y(k,250)
         mat(k,767) = -rxt(k,454)*y(k,250)
         mat(k,281) = -rxt(k,455)*y(k,250)
         mat(k,277) = -rxt(k,456)*y(k,250)
         mat(k,393) = -rxt(k,457)*y(k,250)
         mat(k,165) = -rxt(k,458)*y(k,250)
         mat(k,681) = -rxt(k,470)*y(k,250)
         mat(k,633) = -rxt(k,471)*y(k,250)
         mat(k,494) = -rxt(k,478)*y(k,250)
         mat(k,903) = -rxt(k,480)*y(k,250)
         mat(k,833) = -rxt(k,481)*y(k,250)
         mat(k,458) = -rxt(k,482)*y(k,250)
         mat(k,1184) = -rxt(k,483)*y(k,250)
         mat(k,239) = -rxt(k,489)*y(k,250)
         mat(k,194) = -rxt(k,492)*y(k,250)
         mat(k,470) = -rxt(k,495)*y(k,250)
         mat(k,284) = -rxt(k,496)*y(k,250)
         mat(k,379) = -rxt(k,498)*y(k,250)
         mat(k,305) = -rxt(k,502)*y(k,250)
         mat(k,231) = -rxt(k,503)*y(k,250)
         mat(k,203) = -rxt(k,515)*y(k,250)
         mat(k,384) = -rxt(k,518)*y(k,250)
         mat(k,708) = -rxt(k,522)*y(k,250)
         mat(k,226) = -rxt(k,523)*y(k,250)
         mat(k,250) = -rxt(k,525)*y(k,250)
         mat(k,795) = -rxt(k,528)*y(k,250)
         mat(k,255) = -rxt(k,531)*y(k,250)
         mat(k,506) = -rxt(k,532)*y(k,250)
         mat(k,1018) = -rxt(k,535)*y(k,250)
         mat(k,990) = -rxt(k,538)*y(k,250)
         mat(k,465) = -rxt(k,543)*y(k,250)
         mat(k,756) = -rxt(k,547)*y(k,250)
         mat(k,727) = -rxt(k,548)*y(k,250)
         mat(k,566) = -rxt(k,552)*y(k,250)
         mat(k,1070) = -rxt(k,554)*y(k,250)
         mat(k,1167) = -rxt(k,555)*y(k,250)
         mat(k,1094) = -(rxt(k,561) + rxt(k,575)) * y(k,250)
         mat(k,429) = -rxt(k,563)*y(k,250)
         mat(k,1105) = -rxt(k,565)*y(k,250)
         mat(k,607) = -rxt(k,569)*y(k,250)
         mat(k,1549) = -rxt(k,574)*y(k,250)
         mat(k,110) = -rxt(k,594)*y(k,250)
         mat(k,1018) = mat(k,1018) + .630_r8*rxt(k,534)*y(k,164)
         mat(k,330) = mat(k,330) + .650_r8*rxt(k,370)*y(k,250)
         mat(k,656) = mat(k,656) + .130_r8*rxt(k,372)*y(k,164)
         mat(k,389) = mat(k,389) + .500_r8*rxt(k,378)*y(k,250)
         mat(k,1137) = mat(k,1137) + .360_r8*rxt(k,401)*y(k,164)
         mat(k,1692) = mat(k,1692) + rxt(k,350)*y(k,163)
         mat(k,499) = mat(k,499) + .300_r8*rxt(k,357)*y(k,250)
         mat(k,1570) = mat(k,1570) + rxt(k,364)*y(k,249)
         mat(k,2382) = rxt(k,209)*y(k,235)
         mat(k,954) = rxt(k,308)*y(k,262)
         mat(k,2463) = rxt(k,170)*y(k,164) + 2.000_r8*rxt(k,165)*y(k,235)
         mat(k,1525) = mat(k,1525) + rxt(k,162)*y(k,163) + rxt(k,153)*y(k,249)
         mat(k,694) = mat(k,694) + rxt(k,163)*y(k,163)
         mat(k,1509) = mat(k,1509) + rxt(k,259)*y(k,163) + rxt(k,265)*y(k,249)
         mat(k,1717) = mat(k,1717) + rxt(k,226)*y(k,163) + rxt(k,238)*y(k,249)
         mat(k,207) = mat(k,207) + rxt(k,367)*y(k,249)
         mat(k,1590) = rxt(k,261)*y(k,163)
         mat(k,1741) = mat(k,1741) + rxt(k,229)*y(k,163)
         mat(k,903) = mat(k,903) + .320_r8*rxt(k,479)*y(k,164)
         mat(k,833) = mat(k,833) + .600_r8*rxt(k,481)*y(k,250)
         mat(k,1330) = mat(k,1330) + .240_r8*rxt(k,432)*y(k,164)
         mat(k,356) = mat(k,356) + .100_r8*rxt(k,434)*y(k,250)
         mat(k,990) = mat(k,990) + .630_r8*rxt(k,537)*y(k,164)
         mat(k,1443) = mat(k,1443) + .360_r8*rxt(k,446)*y(k,164)
         mat(k,1921) = rxt(k,193)*y(k,235)
         mat(k,2662) = mat(k,2662) + rxt(k,188)*y(k,235)
         mat(k,2498) = mat(k,2498) + rxt(k,350)*y(k,51) + rxt(k,162)*y(k,93) &
                      + rxt(k,163)*y(k,95) + rxt(k,259)*y(k,97) + rxt(k,226)*y(k,101) &
                      + rxt(k,261)*y(k,108) + rxt(k,229)*y(k,109) + rxt(k,168) &
                      *y(k,235)
         mat(k,2592) = mat(k,2592) + .630_r8*rxt(k,534)*y(k,6) + .130_r8*rxt(k,372) &
                      *y(k,28) + .360_r8*rxt(k,401)*y(k,33) + rxt(k,170)*y(k,92) &
                      + .320_r8*rxt(k,479)*y(k,127) + .240_r8*rxt(k,432)*y(k,134) &
                      + .630_r8*rxt(k,537)*y(k,139) + .360_r8*rxt(k,446)*y(k,140) &
                      + rxt(k,169)*y(k,235)
         mat(k,648) = mat(k,648) + .500_r8*rxt(k,414)*y(k,250)
         mat(k,239) = mat(k,239) + .500_r8*rxt(k,489)*y(k,250)
         mat(k,619) = .400_r8*rxt(k,490)*y(k,235)
         mat(k,1495) = .450_r8*rxt(k,386)*y(k,235)
         mat(k,858) = .400_r8*rxt(k,504)*y(k,235)
         mat(k,2321) = mat(k,2321) + rxt(k,209)*y(k,70) + 2.000_r8*rxt(k,165)*y(k,92) &
                      + rxt(k,193)*y(k,153) + rxt(k,188)*y(k,155) + rxt(k,168) &
                      *y(k,163) + rxt(k,169)*y(k,164) + .400_r8*rxt(k,490)*y(k,220) &
                      + .450_r8*rxt(k,386)*y(k,229) + .400_r8*rxt(k,504)*y(k,231) &
                      + .450_r8*rxt(k,437)*y(k,244) + .400_r8*rxt(k,510)*y(k,245) &
                      + .200_r8*rxt(k,441)*y(k,246) + .150_r8*rxt(k,416)*y(k,253)
         mat(k,1463) = .450_r8*rxt(k,437)*y(k,235)
         mat(k,965) = .400_r8*rxt(k,510)*y(k,235)
         mat(k,775) = .200_r8*rxt(k,441)*y(k,235)
         mat(k,1967) = rxt(k,364)*y(k,64) + rxt(k,153)*y(k,93) + rxt(k,265)*y(k,97) &
                      + rxt(k,238)*y(k,101) + rxt(k,367)*y(k,102) &
                      + 2.000_r8*rxt(k,154)*y(k,262)
         mat(k,2149) = mat(k,2149) + .650_r8*rxt(k,370)*y(k,27) + .500_r8*rxt(k,378) &
                      *y(k,31) + .300_r8*rxt(k,357)*y(k,63) + .600_r8*rxt(k,481) &
                      *y(k,131) + .100_r8*rxt(k,434)*y(k,135) + .500_r8*rxt(k,414) &
                      *y(k,177) + .500_r8*rxt(k,489)*y(k,213)
         mat(k,1278) = .150_r8*rxt(k,416)*y(k,235)
         mat(k,2691) = rxt(k,308)*y(k,89) + 2.000_r8*rxt(k,154)*y(k,249)
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
         mat(k,529) = -(rxt(k,513)*y(k,235) + rxt(k,514)*y(k,153))
         mat(k,2251) = -rxt(k,513)*y(k,251)
         mat(k,1855) = -rxt(k,514)*y(k,251)
         mat(k,229) = .200_r8*rxt(k,503)*y(k,250)
         mat(k,201) = .140_r8*rxt(k,515)*y(k,250)
         mat(k,382) = rxt(k,518)*y(k,250)
         mat(k,2052) = .200_r8*rxt(k,503)*y(k,82) + .140_r8*rxt(k,515)*y(k,173) &
                      + rxt(k,518)*y(k,174)
         mat(k,861) = -(rxt(k,412)*y(k,235) + rxt(k,413)*y(k,153))
         mat(k,2275) = -rxt(k,412)*y(k,252)
         mat(k,1877) = -rxt(k,413)*y(k,252)
         mat(k,1122) = rxt(k,419)*y(k,250)
         mat(k,644) = .500_r8*rxt(k,414)*y(k,250)
         mat(k,2091) = rxt(k,419)*y(k,33) + .500_r8*rxt(k,414)*y(k,177)
         mat(k,1271) = -(rxt(k,415)*y(k,230) + rxt(k,416)*y(k,235) + rxt(k,417) &
                      *y(k,153))
         mat(k,1625) = -rxt(k,415)*y(k,253)
         mat(k,2296) = -rxt(k,416)*y(k,253)
         mat(k,1899) = -rxt(k,417)*y(k,253)
         mat(k,1014) = .060_r8*rxt(k,534)*y(k,164)
         mat(k,1074) = rxt(k,420)*y(k,250)
         mat(k,986) = .060_r8*rxt(k,537)*y(k,164)
         mat(k,2572) = .060_r8*rxt(k,534)*y(k,6) + .060_r8*rxt(k,537)*y(k,139)
         mat(k,473) = rxt(k,418)*y(k,250)
         mat(k,1163) = .150_r8*rxt(k,555)*y(k,250)
         mat(k,2123) = rxt(k,420)*y(k,57) + rxt(k,418)*y(k,178) + .150_r8*rxt(k,555) &
                      *y(k,210)
         mat(k,1217) = -(rxt(k,544)*y(k,230) + rxt(k,545)*y(k,235) + rxt(k,546) &
                      *y(k,153))
         mat(k,1623) = -rxt(k,544)*y(k,254)
         mat(k,2293) = -rxt(k,545)*y(k,254)
         mat(k,1896) = -rxt(k,546)*y(k,254)
         mat(k,2634) = .500_r8*rxt(k,553)*y(k,209)
         mat(k,754) = rxt(k,547)*y(k,250)
         mat(k,1069) = .500_r8*rxt(k,553)*y(k,155) + rxt(k,554)*y(k,250)
         mat(k,2119) = rxt(k,547)*y(k,206) + rxt(k,554)*y(k,209)
         mat(k,1055) = -(rxt(k,549)*y(k,230) + rxt(k,550)*y(k,235) + rxt(k,551) &
                      *y(k,153))
         mat(k,1613) = -rxt(k,549)*y(k,255)
         mat(k,2284) = -rxt(k,550)*y(k,255)
         mat(k,1886) = -rxt(k,551)*y(k,255)
         mat(k,1008) = rxt(k,535)*y(k,250)
         mat(k,980) = rxt(k,538)*y(k,250)
         mat(k,562) = rxt(k,552)*y(k,250)
         mat(k,2105) = rxt(k,535)*y(k,6) + rxt(k,538)*y(k,139) + rxt(k,552)*y(k,208)
         mat(k,817) = -(rxt(k,520)*y(k,235) + rxt(k,521)*y(k,153))
         mat(k,2271) = -rxt(k,520)*y(k,256)
         mat(k,1874) = -rxt(k,521)*y(k,256)
         mat(k,704) = rxt(k,522)*y(k,250)
         mat(k,225) = .650_r8*rxt(k,523)*y(k,250)
         mat(k,2087) = rxt(k,522)*y(k,211) + .650_r8*rxt(k,523)*y(k,212)
         mat(k,97) = -(rxt(k,649)*y(k,235) + rxt(k,650)*y(k,153))
         mat(k,2229) = -rxt(k,649)*y(k,257)
         mat(k,1844) = -rxt(k,650)*y(k,257)
         mat(k,220) = rxt(k,648)*y(k,250)
         mat(k,1992) = rxt(k,648)*y(k,212)
         mat(k,1288) = -(rxt(k,484)*y(k,229) + rxt(k,485)*y(k,230) + rxt(k,486) &
                      *y(k,235) + rxt(k,487)*y(k,153) + rxt(k,488)*y(k,155))
         mat(k,1479) = -rxt(k,484)*y(k,258)
         mat(k,1626) = -rxt(k,485)*y(k,258)
         mat(k,2297) = -rxt(k,486)*y(k,258)
         mat(k,1900) = -rxt(k,487)*y(k,258)
         mat(k,2639) = -rxt(k,488)*y(k,258)
         mat(k,276) = rxt(k,456)*y(k,250)
         mat(k,392) = rxt(k,457)*y(k,250)
         mat(k,164) = rxt(k,458)*y(k,250)
         mat(k,829) = .400_r8*rxt(k,481)*y(k,250)
         mat(k,238) = .500_r8*rxt(k,489)*y(k,250)
         mat(k,2124) = rxt(k,456)*y(k,112) + rxt(k,457)*y(k,114) + rxt(k,458)*y(k,122) &
                      + .400_r8*rxt(k,481)*y(k,131) + .500_r8*rxt(k,489)*y(k,213)
         mat(k,841) = -(rxt(k,526)*y(k,235) + rxt(k,527)*y(k,153))
         mat(k,2273) = -rxt(k,526)*y(k,259)
         mat(k,1875) = -rxt(k,527)*y(k,259)
         mat(k,247) = .560_r8*rxt(k,525)*y(k,250)
         mat(k,788) = rxt(k,528)*y(k,250)
         mat(k,2089) = .560_r8*rxt(k,525)*y(k,214) + rxt(k,528)*y(k,215)
         mat(k,103) = -(rxt(k,653)*y(k,235) + rxt(k,654)*y(k,153))
         mat(k,2230) = -rxt(k,653)*y(k,260)
         mat(k,1845) = -rxt(k,654)*y(k,260)
         mat(k,242) = rxt(k,652)*y(k,250)
         mat(k,1993) = rxt(k,652)*y(k,214)
         mat(k,586) = -(rxt(k,529)*y(k,235) + rxt(k,530)*y(k,153))
         mat(k,2256) = -rxt(k,529)*y(k,261)
         mat(k,1860) = -rxt(k,530)*y(k,261)
         mat(k,254) = .300_r8*rxt(k,531)*y(k,250)
         mat(k,503) = rxt(k,532)*y(k,250)
         mat(k,2060) = .300_r8*rxt(k,531)*y(k,216) + rxt(k,532)*y(k,217)
         mat(k,2703) = -(rxt(k,154)*y(k,249) + rxt(k,308)*y(k,89) + rxt(k,576) &
                      *y(k,183))
         mat(k,1979) = -rxt(k,154)*y(k,262)
         mat(k,957) = -rxt(k,308)*y(k,262)
         mat(k,297) = -rxt(k,576)*y(k,262)
         mat(k,292) = rxt(k,315)*y(k,250)
         mat(k,340) = rxt(k,380)*y(k,250)
         mat(k,483) = rxt(k,405)*y(k,250)
         mat(k,346) = rxt(k,406)*y(k,250)
         mat(k,560) = rxt(k,317)*y(k,250)
         mat(k,352) = rxt(k,320)*y(k,250)
         mat(k,1704) = rxt(k,351)*y(k,250)
         mat(k,722) = rxt(k,322)*y(k,250)
         mat(k,153) = rxt(k,323)*y(k,250)
         mat(k,1195) = rxt(k,382)*y(k,250)
         mat(k,447) = rxt(k,325)*y(k,250)
         mat(k,1078) = rxt(k,420)*y(k,250)
         mat(k,1343) = rxt(k,408)*y(k,250)
         mat(k,781) = rxt(k,388)*y(k,250)
         mat(k,673) = rxt(k,389)*y(k,250)
         mat(k,439) = rxt(k,327)*y(k,250)
         mat(k,501) = rxt(k,357)*y(k,250)
         mat(k,1576) = rxt(k,358)*y(k,250)
         mat(k,1120) = rxt(k,334)*y(k,235)
         mat(k,404) = rxt(k,339)*y(k,250)
         mat(k,2475) = rxt(k,166)*y(k,235)
         mat(k,1530) = rxt(k,171)*y(k,250)
         mat(k,698) = rxt(k,172)*y(k,250)
         mat(k,1516) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,108) + ( &
                      + rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,109) + ( &
                      + rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,110) &
                      + rxt(k,260)*y(k,250)
         mat(k,322) = rxt(k,342)*y(k,250)
         mat(k,1727) = (rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,108) + ( &
                      + rxt(k,585)+rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,109) + ( &
                      + rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,110) &
                      + rxt(k,227)*y(k,250)
         mat(k,1044) = rxt(k,360)*y(k,250)
         mat(k,1267) = (rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,108) + ( &
                      + rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,109) + ( &
                      + rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,110) &
                      + rxt(k,295)*y(k,250)
         mat(k,1776) = rxt(k,201)*y(k,250)
         mat(k,521) = rxt(k,179)*y(k,250)
         mat(k,1599) = (rxt(k,584)+rxt(k,658)+rxt(k,671)+rxt(k,680))*y(k,97) + ( &
                      + rxt(k,586)+rxt(k,657)+rxt(k,670)+rxt(k,679))*y(k,101) + ( &
                      + rxt(k,588)+rxt(k,659)+rxt(k,672)+rxt(k,681))*y(k,105)
         mat(k,1751) = (rxt(k,583)+rxt(k,660)+rxt(k,668)+rxt(k,677))*y(k,97) + ( &
                      + rxt(k,585)+rxt(k,656)+rxt(k,667)+rxt(k,676))*y(k,101) + ( &
                      + rxt(k,587)+rxt(k,661)+rxt(k,669)+rxt(k,678))*y(k,105) &
                      + rxt(k,230)*y(k,250)
         mat(k,1675) = (rxt(k,591)+rxt(k,687)+rxt(k,691)+rxt(k,695))*y(k,97) + ( &
                      + rxt(k,590)+rxt(k,686)+rxt(k,690)+rxt(k,694))*y(k,101) + ( &
                      + rxt(k,592)+rxt(k,688)+rxt(k,692)+rxt(k,696))*y(k,105) &
                      + rxt(k,268)*y(k,250)
         mat(k,1334) = .500_r8*rxt(k,433)*y(k,250)
         mat(k,111) = rxt(k,594)*y(k,250)
         mat(k,650) = rxt(k,414)*y(k,250)
         mat(k,477) = rxt(k,418)*y(k,250)
         mat(k,2333) = rxt(k,334)*y(k,68) + rxt(k,166)*y(k,92) + rxt(k,173)*y(k,250)
         mat(k,2161) = rxt(k,315)*y(k,29) + rxt(k,380)*y(k,32) + rxt(k,405)*y(k,34) &
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
                      *y(k,149) + rxt(k,414)*y(k,177) + rxt(k,418)*y(k,178) &
                      + rxt(k,173)*y(k,235) + 2.000_r8*rxt(k,176)*y(k,250)
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
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 85) = mat(k, 85) + lmat(k, 85)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = lmat(k, 104)
         mat(k, 105) = lmat(k, 105)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = lmat(k, 107)
         mat(k, 108) = lmat(k, 108)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 129) = mat(k, 129) + lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = lmat(k, 157)
         mat(k, 158) = lmat(k, 158)
         mat(k, 159) = lmat(k, 159)
         mat(k, 160) = lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = lmat(k, 162)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 196) = lmat(k, 196)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 212) = lmat(k, 212)
         mat(k, 213) = mat(k, 213) + lmat(k, 213)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 240) = lmat(k, 240)
         mat(k, 241) = lmat(k, 241)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 274) = mat(k, 274) + lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 285) = lmat(k, 285)
         mat(k, 286) = lmat(k, 286)
         mat(k, 287) = lmat(k, 287)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 295) = lmat(k, 295)
         mat(k, 296) = lmat(k, 296)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 303) = lmat(k, 303)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 307) = lmat(k, 307)
         mat(k, 308) = lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 312) = lmat(k, 312)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 315) = lmat(k, 315)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = lmat(k, 325)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = lmat(k, 333)
         mat(k, 334) = lmat(k, 334)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 341) = mat(k, 341) + lmat(k, 341)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 351) = mat(k, 351) + lmat(k, 351)
         mat(k, 353) = mat(k, 353) + lmat(k, 353)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 360) = lmat(k, 360)
         mat(k, 361) = lmat(k, 361)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 363) = lmat(k, 363)
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
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = mat(k, 381) + lmat(k, 381)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = mat(k, 384) + lmat(k, 384)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = mat(k, 386) + lmat(k, 386)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 390) = lmat(k, 390)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 394) = lmat(k, 394)
         mat(k, 395) = lmat(k, 395)
         mat(k, 396) = lmat(k, 396)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 399) = lmat(k, 399)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 405) = lmat(k, 405)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 408) = mat(k, 408) + lmat(k, 408)
         mat(k, 409) = lmat(k, 409)
         mat(k, 410) = lmat(k, 410)
         mat(k, 411) = lmat(k, 411)
         mat(k, 412) = lmat(k, 412)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 421) = lmat(k, 421)
         mat(k, 422) = lmat(k, 422)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 425) = lmat(k, 425)
         mat(k, 427) = mat(k, 427) + lmat(k, 427)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 433) = lmat(k, 433)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 440) = mat(k, 440) + lmat(k, 440)
         mat(k, 442) = lmat(k, 442)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 457) = lmat(k, 457)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 461) = lmat(k, 461)
         mat(k, 463) = lmat(k, 463)
         mat(k, 464) = lmat(k, 464)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 466) = mat(k, 466) + lmat(k, 466)
         mat(k, 467) = lmat(k, 467)
         mat(k, 469) = lmat(k, 469)
         mat(k, 470) = mat(k, 470) + lmat(k, 470)
         mat(k, 471) = lmat(k, 471)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 474) = lmat(k, 474)
         mat(k, 475) = lmat(k, 475)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 478) = mat(k, 478) + lmat(k, 478)
         mat(k, 480) = lmat(k, 480)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 482) = lmat(k, 482)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = lmat(k, 486)
         mat(k, 487) = lmat(k, 487)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = mat(k, 490) + lmat(k, 490)
         mat(k, 491) = lmat(k, 491)
         mat(k, 493) = lmat(k, 493)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 498) = mat(k, 498) + lmat(k, 498)
         mat(k, 499) = mat(k, 499) + lmat(k, 499)
         mat(k, 500) = lmat(k, 500)
         mat(k, 502) = mat(k, 502) + lmat(k, 502)
         mat(k, 504) = lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = mat(k, 506) + lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 516) = mat(k, 516) + lmat(k, 516)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 519) = lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 524) = lmat(k, 524)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 538) = mat(k, 538) + lmat(k, 538)
         mat(k, 539) = lmat(k, 539)
         mat(k, 540) = lmat(k, 540)
         mat(k, 541) = lmat(k, 541)
         mat(k, 542) = lmat(k, 542)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 547) = mat(k, 547) + lmat(k, 547)
         mat(k, 551) = lmat(k, 551)
         mat(k, 552) = lmat(k, 552)
         mat(k, 553) = lmat(k, 553)
         mat(k, 554) = mat(k, 554) + lmat(k, 554)
         mat(k, 559) = mat(k, 559) + lmat(k, 559)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 563) = lmat(k, 563)
         mat(k, 564) = lmat(k, 564)
         mat(k, 565) = lmat(k, 565)
         mat(k, 566) = mat(k, 566) + lmat(k, 566)
         mat(k, 567) = lmat(k, 567)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 571) = mat(k, 571) + lmat(k, 571)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 586) = mat(k, 586) + lmat(k, 586)
         mat(k, 595) = mat(k, 595) + lmat(k, 595)
         mat(k, 599) = lmat(k, 599)
         mat(k, 600) = mat(k, 600) + lmat(k, 600)
         mat(k, 602) = lmat(k, 602)
         mat(k, 603) = lmat(k, 603)
         mat(k, 605) = mat(k, 605) + lmat(k, 605)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 609) = lmat(k, 609)
         mat(k, 611) = lmat(k, 611)
         mat(k, 612) = lmat(k, 612)
         mat(k, 613) = lmat(k, 613)
         mat(k, 615) = mat(k, 615) + lmat(k, 615)
         mat(k, 621) = lmat(k, 621)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 624) = mat(k, 624) + lmat(k, 624)
         mat(k, 626) = lmat(k, 626)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 632) = lmat(k, 632)
         mat(k, 635) = mat(k, 635) + lmat(k, 635)
         mat(k, 636) = lmat(k, 636)
         mat(k, 637) = lmat(k, 637)
         mat(k, 638) = lmat(k, 638)
         mat(k, 640) = lmat(k, 640)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 643) = mat(k, 643) + lmat(k, 643)
         mat(k, 645) = lmat(k, 645)
         mat(k, 647) = lmat(k, 647)
         mat(k, 648) = mat(k, 648) + lmat(k, 648)
         mat(k, 649) = lmat(k, 649)
         mat(k, 651) = mat(k, 651) + lmat(k, 651)
         mat(k, 659) = mat(k, 659) + lmat(k, 659)
         mat(k, 667) = mat(k, 667) + lmat(k, 667)
         mat(k, 668) = mat(k, 668) + lmat(k, 668)
         mat(k, 670) = lmat(k, 670)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 674) = mat(k, 674) + lmat(k, 674)
         mat(k, 679) = lmat(k, 679)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 686) = lmat(k, 686)
         mat(k, 688) = lmat(k, 688)
         mat(k, 692) = mat(k, 692) + lmat(k, 692)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
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
         mat(k, 715) = lmat(k, 715)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 723) = mat(k, 723) + lmat(k, 723)
         mat(k, 724) = mat(k, 724) + lmat(k, 724)
         mat(k, 725) = lmat(k, 725)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 728) = lmat(k, 728)
         mat(k, 731) = mat(k, 731) + lmat(k, 731)
         mat(k, 737) = lmat(k, 737)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 742) = mat(k, 742) + lmat(k, 742)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 745) = mat(k, 745) + lmat(k, 745)
         mat(k, 747) = lmat(k, 747)
         mat(k, 748) = mat(k, 748) + lmat(k, 748)
         mat(k, 749) = lmat(k, 749)
         mat(k, 750) = lmat(k, 750)
         mat(k, 751) = lmat(k, 751)
         mat(k, 752) = lmat(k, 752)
         mat(k, 753) = lmat(k, 753)
         mat(k, 755) = lmat(k, 755)
         mat(k, 756) = mat(k, 756) + lmat(k, 756)
         mat(k, 757) = lmat(k, 757)
         mat(k, 758) = lmat(k, 758)
         mat(k, 759) = mat(k, 759) + lmat(k, 759)
         mat(k, 763) = lmat(k, 763)
         mat(k, 764) = lmat(k, 764)
         mat(k, 766) = lmat(k, 766)
         mat(k, 767) = mat(k, 767) + lmat(k, 767)
         mat(k, 768) = lmat(k, 768)
         mat(k, 770) = mat(k, 770) + lmat(k, 770)
         mat(k, 777) = mat(k, 777) + lmat(k, 777)
         mat(k, 782) = lmat(k, 782)
         mat(k, 783) = lmat(k, 783)
         mat(k, 784) = lmat(k, 784)
         mat(k, 785) = lmat(k, 785)
         mat(k, 786) = mat(k, 786) + lmat(k, 786)
         mat(k, 791) = lmat(k, 791)
         mat(k, 793) = lmat(k, 793)
         mat(k, 795) = mat(k, 795) + lmat(k, 795)
         mat(k, 796) = lmat(k, 796)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 806) = mat(k, 806) + lmat(k, 806)
         mat(k, 817) = mat(k, 817) + lmat(k, 817)
         mat(k, 828) = mat(k, 828) + lmat(k, 828)
         mat(k, 830) = lmat(k, 830)
         mat(k, 831) = lmat(k, 831)
         mat(k, 832) = lmat(k, 832)
         mat(k, 833) = mat(k, 833) + lmat(k, 833)
         mat(k, 834) = lmat(k, 834)
         mat(k, 841) = mat(k, 841) + lmat(k, 841)
         mat(k, 852) = mat(k, 852) + lmat(k, 852)
         mat(k, 861) = mat(k, 861) + lmat(k, 861)
         mat(k, 870) = mat(k, 870) + lmat(k, 870)
         mat(k, 874) = lmat(k, 874)
         mat(k, 875) = lmat(k, 875)
         mat(k, 876) = lmat(k, 876)
         mat(k, 880) = mat(k, 880) + lmat(k, 880)
         mat(k, 891) = mat(k, 891) + lmat(k, 891)
         mat(k, 907) = mat(k, 907) + lmat(k, 907)
         mat(k, 909) = lmat(k, 909)
         mat(k, 910) = lmat(k, 910)
         mat(k, 911) = mat(k, 911) + lmat(k, 911)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 915) = mat(k, 915) + lmat(k, 915)
         mat(k, 916) = lmat(k, 916)
         mat(k, 917) = lmat(k, 917)
         mat(k, 918) = lmat(k, 918)
         mat(k, 920) = mat(k, 920) + lmat(k, 920)
         mat(k, 923) = mat(k, 923) + lmat(k, 923)
         mat(k, 924) = mat(k, 924) + lmat(k, 924)
         mat(k, 926) = mat(k, 926) + lmat(k, 926)
         mat(k, 936) = lmat(k, 936)
         mat(k, 937) = lmat(k, 937)
         mat(k, 938) = mat(k, 938) + lmat(k, 938)
         mat(k, 940) = lmat(k, 940)
         mat(k, 941) = lmat(k, 941)
         mat(k, 942) = lmat(k, 942)
         mat(k, 944) = lmat(k, 944)
         mat(k, 945) = mat(k, 945) + lmat(k, 945)
         mat(k, 947) = mat(k, 947) + lmat(k, 947)
         mat(k, 949) = mat(k, 949) + lmat(k, 949)
         mat(k, 958) = mat(k, 958) + lmat(k, 958)
         mat(k, 976) = mat(k, 976) + lmat(k, 976)
         mat(k,1004) = mat(k,1004) + lmat(k,1004)
         mat(k,1028) = mat(k,1028) + lmat(k,1028)
         mat(k,1039) = mat(k,1039) + lmat(k,1039)
         mat(k,1045) = lmat(k,1045)
         mat(k,1047) = mat(k,1047) + lmat(k,1047)
         mat(k,1048) = mat(k,1048) + lmat(k,1048)
         mat(k,1050) = mat(k,1050) + lmat(k,1050)
         mat(k,1055) = mat(k,1055) + lmat(k,1055)
         mat(k,1065) = mat(k,1065) + lmat(k,1065)
         mat(k,1067) = lmat(k,1067)
         mat(k,1068) = lmat(k,1068)
         mat(k,1071) = lmat(k,1071)
         mat(k,1073) = mat(k,1073) + lmat(k,1073)
         mat(k,1075) = lmat(k,1075)
         mat(k,1076) = lmat(k,1076)
         mat(k,1079) = mat(k,1079) + lmat(k,1079)
         mat(k,1080) = mat(k,1080) + lmat(k,1080)
         mat(k,1082) = mat(k,1082) + lmat(k,1082)
         mat(k,1083) = mat(k,1083) + lmat(k,1083)
         mat(k,1084) = lmat(k,1084)
         mat(k,1086) = mat(k,1086) + lmat(k,1086)
         mat(k,1087) = mat(k,1087) + lmat(k,1087)
         mat(k,1089) = mat(k,1089) + lmat(k,1089)
         mat(k,1090) = mat(k,1090) + lmat(k,1090)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1104) = lmat(k,1104)
         mat(k,1107) = lmat(k,1107)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1126) = mat(k,1126) + lmat(k,1126)
         mat(k,1147) = mat(k,1147) + lmat(k,1147)
         mat(k,1159) = mat(k,1159) + lmat(k,1159)
         mat(k,1160) = mat(k,1160) + lmat(k,1160)
         mat(k,1161) = mat(k,1161) + lmat(k,1161)
         mat(k,1162) = mat(k,1162) + lmat(k,1162)
         mat(k,1163) = mat(k,1163) + lmat(k,1163)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1166) = mat(k,1166) + lmat(k,1166)
         mat(k,1168) = mat(k,1168) + lmat(k,1168)
         mat(k,1171) = lmat(k,1171)
         mat(k,1175) = mat(k,1175) + lmat(k,1175)
         mat(k,1181) = lmat(k,1181)
         mat(k,1182) = lmat(k,1182)
         mat(k,1185) = mat(k,1185) + lmat(k,1185)
         mat(k,1187) = mat(k,1187) + lmat(k,1187)
         mat(k,1188) = lmat(k,1188)
         mat(k,1190) = lmat(k,1190)
         mat(k,1193) = lmat(k,1193)
         mat(k,1198) = mat(k,1198) + lmat(k,1198)
         mat(k,1199) = lmat(k,1199)
         mat(k,1200) = mat(k,1200) + lmat(k,1200)
         mat(k,1202) = mat(k,1202) + lmat(k,1202)
         mat(k,1205) = mat(k,1205) + lmat(k,1205)
         mat(k,1217) = mat(k,1217) + lmat(k,1217)
         mat(k,1230) = lmat(k,1230)
         mat(k,1231) = lmat(k,1231)
         mat(k,1232) = lmat(k,1232)
         mat(k,1233) = lmat(k,1233)
         mat(k,1234) = mat(k,1234) + lmat(k,1234)
         mat(k,1235) = lmat(k,1235)
         mat(k,1237) = lmat(k,1237)
         mat(k,1239) = lmat(k,1239)
         mat(k,1240) = lmat(k,1240)
         mat(k,1241) = lmat(k,1241)
         mat(k,1243) = mat(k,1243) + lmat(k,1243)
         mat(k,1247) = mat(k,1247) + lmat(k,1247)
         mat(k,1249) = lmat(k,1249)
         mat(k,1250) = lmat(k,1250)
         mat(k,1252) = mat(k,1252) + lmat(k,1252)
         mat(k,1255) = mat(k,1255) + lmat(k,1255)
         mat(k,1263) = mat(k,1263) + lmat(k,1263)
         mat(k,1264) = lmat(k,1264)
         mat(k,1271) = mat(k,1271) + lmat(k,1271)
         mat(k,1288) = mat(k,1288) + lmat(k,1288)
         mat(k,1308) = mat(k,1308) + lmat(k,1308)
         mat(k,1323) = mat(k,1323) + lmat(k,1323)
         mat(k,1324) = mat(k,1324) + lmat(k,1324)
         mat(k,1327) = mat(k,1327) + lmat(k,1327)
         mat(k,1328) = mat(k,1328) + lmat(k,1328)
         mat(k,1329) = mat(k,1329) + lmat(k,1329)
         mat(k,1331) = mat(k,1331) + lmat(k,1331)
         mat(k,1335) = mat(k,1335) + lmat(k,1335)
         mat(k,1336) = mat(k,1336) + lmat(k,1336)
         mat(k,1337) = mat(k,1337) + lmat(k,1337)
         mat(k,1340) = lmat(k,1340)
         mat(k,1356) = mat(k,1356) + lmat(k,1356)
         mat(k,1372) = lmat(k,1372)
         mat(k,1389) = mat(k,1389) + lmat(k,1389)
         mat(k,1400) = mat(k,1400) + lmat(k,1400)
         mat(k,1413) = mat(k,1413) + lmat(k,1413)
         mat(k,1428) = lmat(k,1428)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1434) = mat(k,1434) + lmat(k,1434)
         mat(k,1436) = mat(k,1436) + lmat(k,1436)
         mat(k,1438) = lmat(k,1438)
         mat(k,1456) = mat(k,1456) + lmat(k,1456)
         mat(k,1488) = mat(k,1488) + lmat(k,1488)
         mat(k,1504) = mat(k,1504) + lmat(k,1504)
         mat(k,1513) = mat(k,1513) + lmat(k,1513)
         mat(k,1515) = mat(k,1515) + lmat(k,1515)
         mat(k,1519) = mat(k,1519) + lmat(k,1519)
         mat(k,1532) = mat(k,1532) + lmat(k,1532)
         mat(k,1539) = mat(k,1539) + lmat(k,1539)
         mat(k,1544) = lmat(k,1544)
         mat(k,1545) = mat(k,1545) + lmat(k,1545)
         mat(k,1546) = mat(k,1546) + lmat(k,1546)
         mat(k,1555) = lmat(k,1555)
         mat(k,1560) = lmat(k,1560)
         mat(k,1562) = lmat(k,1562)
         mat(k,1563) = mat(k,1563) + lmat(k,1563)
         mat(k,1564) = mat(k,1564) + lmat(k,1564)
         mat(k,1565) = mat(k,1565) + lmat(k,1565)
         mat(k,1566) = mat(k,1566) + lmat(k,1566)
         mat(k,1570) = mat(k,1570) + lmat(k,1570)
         mat(k,1573) = mat(k,1573) + lmat(k,1573)
         mat(k,1574) = lmat(k,1574)
         mat(k,1576) = mat(k,1576) + lmat(k,1576)
         mat(k,1578) = mat(k,1578) + lmat(k,1578)
         mat(k,1579) = mat(k,1579) + lmat(k,1579)
         mat(k,1583) = mat(k,1583) + lmat(k,1583)
         mat(k,1590) = mat(k,1590) + lmat(k,1590)
         mat(k,1597) = lmat(k,1597)
         mat(k,1637) = mat(k,1637) + lmat(k,1637)
         mat(k,1654) = mat(k,1654) + lmat(k,1654)
         mat(k,1655) = mat(k,1655) + lmat(k,1655)
         mat(k,1660) = mat(k,1660) + lmat(k,1660)
         mat(k,1666) = mat(k,1666) + lmat(k,1666)
         mat(k,1670) = lmat(k,1670)
         mat(k,1678) = mat(k,1678) + lmat(k,1678)
         mat(k,1680) = lmat(k,1680)
         mat(k,1685) = mat(k,1685) + lmat(k,1685)
         mat(k,1699) = mat(k,1699) + lmat(k,1699)
         mat(k,1712) = mat(k,1712) + lmat(k,1712)
         mat(k,1720) = mat(k,1720) + lmat(k,1720)
         mat(k,1723) = mat(k,1723) + lmat(k,1723)
         mat(k,1728) = mat(k,1728) + lmat(k,1728)
         mat(k,1737) = mat(k,1737) + lmat(k,1737)
         mat(k,1741) = mat(k,1741) + lmat(k,1741)
         mat(k,1744) = mat(k,1744) + lmat(k,1744)
         mat(k,1758) = lmat(k,1758)
         mat(k,1760) = mat(k,1760) + lmat(k,1760)
         mat(k,1761) = lmat(k,1761)
         mat(k,1764) = mat(k,1764) + lmat(k,1764)
         mat(k,1817) = mat(k,1817) + lmat(k,1817)
         mat(k,1818) = mat(k,1818) + lmat(k,1818)
         mat(k,1819) = mat(k,1819) + lmat(k,1819)
         mat(k,1821) = mat(k,1821) + lmat(k,1821)
         mat(k,1829) = mat(k,1829) + lmat(k,1829)
         mat(k,1858) = mat(k,1858) + lmat(k,1858)
         mat(k,1919) = mat(k,1919) + lmat(k,1919)
         mat(k,1929) = mat(k,1929) + lmat(k,1929)
         mat(k,1966) = mat(k,1966) + lmat(k,1966)
         mat(k,1975) = mat(k,1975) + lmat(k,1975)
         mat(k,2149) = mat(k,2149) + lmat(k,2149)
         mat(k,2178) = mat(k,2178) + lmat(k,2178)
         mat(k,2185) = mat(k,2185) + lmat(k,2185)
         mat(k,2186) = mat(k,2186) + lmat(k,2186)
         mat(k,2208) = mat(k,2208) + lmat(k,2208)
         mat(k,2210) = mat(k,2210) + lmat(k,2210)
         mat(k,2214) = mat(k,2214) + lmat(k,2214)
         mat(k,2324) = mat(k,2324) + lmat(k,2324)
         mat(k,2333) = mat(k,2333) + lmat(k,2333)
         mat(k,2386) = mat(k,2386) + lmat(k,2386)
         mat(k,2417) = mat(k,2417) + lmat(k,2417)
         mat(k,2418) = mat(k,2418) + lmat(k,2418)
         mat(k,2420) = mat(k,2420) + lmat(k,2420)
         mat(k,2446) = mat(k,2446) + lmat(k,2446)
         mat(k,2466) = mat(k,2466) + lmat(k,2466)
         mat(k,2470) = mat(k,2470) + lmat(k,2470)
         mat(k,2506) = mat(k,2506) + lmat(k,2506)
         mat(k,2508) = mat(k,2508) + lmat(k,2508)
         mat(k,2535) = mat(k,2535) + lmat(k,2535)
         mat(k,2591) = mat(k,2591) + lmat(k,2591)
         mat(k,2600) = mat(k,2600) + lmat(k,2600)
         mat(k,2602) = mat(k,2602) + lmat(k,2602)
         mat(k,2658) = mat(k,2658) + lmat(k,2658)
         mat(k,2659) = mat(k,2659) + lmat(k,2659)
         mat(k,2660) = mat(k,2660) + lmat(k,2660)
         mat(k,2670) = mat(k,2670) + lmat(k,2670)
         mat(k,2673) = mat(k,2673) + lmat(k,2673)
         mat(k,2680) = lmat(k,2680)
         mat(k,2690) = mat(k,2690) + lmat(k,2690)
         mat(k,2691) = mat(k,2691) + lmat(k,2691)
         mat(k,2698) = lmat(k,2698)
         mat(k,2699) = lmat(k,2699)
         mat(k,2703) = mat(k,2703) + lmat(k,2703)
         mat(k, 248) = 0._r8
         mat(k, 249) = 0._r8
         mat(k, 318) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 400) = 0._r8
         mat(k, 511) = 0._r8
         mat(k, 514) = 0._r8
         mat(k, 533) = 0._r8
         mat(k, 572) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 583) = 0._r8
         mat(k, 591) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 739) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 760) = 0._r8
         mat(k, 761) = 0._r8
         mat(k, 765) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 804) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 819) = 0._r8
         mat(k, 821) = 0._r8
         mat(k, 825) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 842) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 847) = 0._r8
         mat(k, 850) = 0._r8
         mat(k, 863) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 886) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 943) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 977) = 0._r8
         mat(k, 979) = 0._r8
         mat(k, 987) = 0._r8
         mat(k, 992) = 0._r8
         mat(k,1005) = 0._r8
         mat(k,1007) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1020) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1032) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1043) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1060) = 0._r8
         mat(k,1063) = 0._r8
         mat(k,1108) = 0._r8
         mat(k,1129) = 0._r8
         mat(k,1130) = 0._r8
         mat(k,1131) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1219) = 0._r8
         mat(k,1222) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1227) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1262) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1318) = 0._r8
         mat(k,1320) = 0._r8
         mat(k,1332) = 0._r8
         mat(k,1341) = 0._r8
         mat(k,1348) = 0._r8
         mat(k,1349) = 0._r8
         mat(k,1350) = 0._r8
         mat(k,1351) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1363) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1375) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1390) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1402) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1411) = 0._r8
         mat(k,1414) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1421) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1441) = 0._r8
         mat(k,1442) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1447) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1460) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1468) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1512) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1534) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1552) = 0._r8
         mat(k,1558) = 0._r8
         mat(k,1568) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1582) = 0._r8
         mat(k,1584) = 0._r8
         mat(k,1586) = 0._r8
         mat(k,1587) = 0._r8
         mat(k,1588) = 0._r8
         mat(k,1589) = 0._r8
         mat(k,1592) = 0._r8
         mat(k,1593) = 0._r8
         mat(k,1594) = 0._r8
         mat(k,1595) = 0._r8
         mat(k,1598) = 0._r8
         mat(k,1610) = 0._r8
         mat(k,1636) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1640) = 0._r8
         mat(k,1641) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1649) = 0._r8
         mat(k,1650) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1652) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1662) = 0._r8
         mat(k,1664) = 0._r8
         mat(k,1665) = 0._r8
         mat(k,1667) = 0._r8
         mat(k,1668) = 0._r8
         mat(k,1671) = 0._r8
         mat(k,1672) = 0._r8
         mat(k,1673) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1681) = 0._r8
         mat(k,1682) = 0._r8
         mat(k,1683) = 0._r8
         mat(k,1684) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1690) = 0._r8
         mat(k,1691) = 0._r8
         mat(k,1693) = 0._r8
         mat(k,1694) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1702) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1718) = 0._r8
         mat(k,1721) = 0._r8
         mat(k,1722) = 0._r8
         mat(k,1725) = 0._r8
         mat(k,1726) = 0._r8
         mat(k,1733) = 0._r8
         mat(k,1734) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1738) = 0._r8
         mat(k,1739) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1746) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1749) = 0._r8
         mat(k,1750) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1759) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1763) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1769) = 0._r8
         mat(k,1770) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1772) = 0._r8
         mat(k,1773) = 0._r8
         mat(k,1774) = 0._r8
         mat(k,1782) = 0._r8
         mat(k,1784) = 0._r8
         mat(k,1786) = 0._r8
         mat(k,1788) = 0._r8
         mat(k,1793) = 0._r8
         mat(k,1796) = 0._r8
         mat(k,1800) = 0._r8
         mat(k,1803) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1805) = 0._r8
         mat(k,1808) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1811) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1813) = 0._r8
         mat(k,1814) = 0._r8
         mat(k,1815) = 0._r8
         mat(k,1816) = 0._r8
         mat(k,1820) = 0._r8
         mat(k,1828) = 0._r8
         mat(k,1833) = 0._r8
         mat(k,1880) = 0._r8
         mat(k,1910) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1913) = 0._r8
         mat(k,1915) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1928) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1963) = 0._r8
         mat(k,1964) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1973) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,2053) = 0._r8
         mat(k,2075) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2090) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2099) = 0._r8
         mat(k,2125) = 0._r8
         mat(k,2148) = 0._r8
         mat(k,2170) = 0._r8
         mat(k,2171) = 0._r8
         mat(k,2172) = 0._r8
         mat(k,2173) = 0._r8
         mat(k,2176) = 0._r8
         mat(k,2184) = 0._r8
         mat(k,2187) = 0._r8
         mat(k,2188) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2202) = 0._r8
         mat(k,2205) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2216) = 0._r8
         mat(k,2217) = 0._r8
         mat(k,2218) = 0._r8
         mat(k,2234) = 0._r8
         mat(k,2252) = 0._r8
         mat(k,2253) = 0._r8
         mat(k,2262) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2285) = 0._r8
         mat(k,2286) = 0._r8
         mat(k,2289) = 0._r8
         mat(k,2291) = 0._r8
         mat(k,2294) = 0._r8
         mat(k,2299) = 0._r8
         mat(k,2304) = 0._r8
         mat(k,2309) = 0._r8
         mat(k,2310) = 0._r8
         mat(k,2317) = 0._r8
         mat(k,2320) = 0._r8
         mat(k,2352) = 0._r8
         mat(k,2356) = 0._r8
         mat(k,2357) = 0._r8
         mat(k,2358) = 0._r8
         mat(k,2363) = 0._r8
         mat(k,2364) = 0._r8
         mat(k,2366) = 0._r8
         mat(k,2367) = 0._r8
         mat(k,2368) = 0._r8
         mat(k,2370) = 0._r8
         mat(k,2372) = 0._r8
         mat(k,2374) = 0._r8
         mat(k,2378) = 0._r8
         mat(k,2380) = 0._r8
         mat(k,2381) = 0._r8
         mat(k,2383) = 0._r8
         mat(k,2387) = 0._r8
         mat(k,2390) = 0._r8
         mat(k,2394) = 0._r8
         mat(k,2403) = 0._r8
         mat(k,2404) = 0._r8
         mat(k,2406) = 0._r8
         mat(k,2407) = 0._r8
         mat(k,2408) = 0._r8
         mat(k,2411) = 0._r8
         mat(k,2419) = 0._r8
         mat(k,2424) = 0._r8
         mat(k,2427) = 0._r8
         mat(k,2428) = 0._r8
         mat(k,2432) = 0._r8
         mat(k,2433) = 0._r8
         mat(k,2434) = 0._r8
         mat(k,2435) = 0._r8
         mat(k,2436) = 0._r8
         mat(k,2439) = 0._r8
         mat(k,2440) = 0._r8
         mat(k,2442) = 0._r8
         mat(k,2444) = 0._r8
         mat(k,2447) = 0._r8
         mat(k,2448) = 0._r8
         mat(k,2452) = 0._r8
         mat(k,2454) = 0._r8
         mat(k,2455) = 0._r8
         mat(k,2456) = 0._r8
         mat(k,2457) = 0._r8
         mat(k,2458) = 0._r8
         mat(k,2459) = 0._r8
         mat(k,2460) = 0._r8
         mat(k,2461) = 0._r8
         mat(k,2462) = 0._r8
         mat(k,2464) = 0._r8
         mat(k,2465) = 0._r8
         mat(k,2467) = 0._r8
         mat(k,2468) = 0._r8
         mat(k,2469) = 0._r8
         mat(k,2472) = 0._r8
         mat(k,2474) = 0._r8
         mat(k,2477) = 0._r8
         mat(k,2481) = 0._r8
         mat(k,2487) = 0._r8
         mat(k,2489) = 0._r8
         mat(k,2490) = 0._r8
         mat(k,2494) = 0._r8
         mat(k,2497) = 0._r8
         mat(k,2510) = 0._r8
         mat(k,2513) = 0._r8
         mat(k,2517) = 0._r8
         mat(k,2518) = 0._r8
         mat(k,2520) = 0._r8
         mat(k,2521) = 0._r8
         mat(k,2522) = 0._r8
         mat(k,2524) = 0._r8
         mat(k,2525) = 0._r8
         mat(k,2526) = 0._r8
         mat(k,2528) = 0._r8
         mat(k,2530) = 0._r8
         mat(k,2533) = 0._r8
         mat(k,2534) = 0._r8
         mat(k,2538) = 0._r8
         mat(k,2552) = 0._r8
         mat(k,2557) = 0._r8
         mat(k,2559) = 0._r8
         mat(k,2560) = 0._r8
         mat(k,2565) = 0._r8
         mat(k,2568) = 0._r8
         mat(k,2570) = 0._r8
         mat(k,2571) = 0._r8
         mat(k,2573) = 0._r8
         mat(k,2576) = 0._r8
         mat(k,2577) = 0._r8
         mat(k,2578) = 0._r8
         mat(k,2580) = 0._r8
         mat(k,2586) = 0._r8
         mat(k,2587) = 0._r8
         mat(k,2588) = 0._r8
         mat(k,2604) = 0._r8
         mat(k,2611) = 0._r8
         mat(k,2612) = 0._r8
         mat(k,2618) = 0._r8
         mat(k,2619) = 0._r8
         mat(k,2621) = 0._r8
         mat(k,2629) = 0._r8
         mat(k,2635) = 0._r8
         mat(k,2638) = 0._r8
         mat(k,2650) = 0._r8
         mat(k,2651) = 0._r8
         mat(k,2652) = 0._r8
         mat(k,2653) = 0._r8
         mat(k,2656) = 0._r8
         mat(k,2657) = 0._r8
         mat(k,2661) = 0._r8
         mat(k,2663) = 0._r8
         mat(k,2664) = 0._r8
         mat(k,2669) = 0._r8
         mat(k,2671) = 0._r8
         mat(k,2672) = 0._r8
         mat(k,2674) = 0._r8
         mat(k,2679) = 0._r8
         mat(k,2681) = 0._r8
         mat(k,2682) = 0._r8
         mat(k,2683) = 0._r8
         mat(k,2684) = 0._r8
         mat(k,2685) = 0._r8
         mat(k,2686) = 0._r8
         mat(k,2687) = 0._r8
         mat(k,2688) = 0._r8
         mat(k,2689) = 0._r8
         mat(k,2692) = 0._r8
         mat(k,2693) = 0._r8
         mat(k,2694) = 0._r8
         mat(k,2695) = 0._r8
         mat(k,2696) = 0._r8
         mat(k,2697) = 0._r8
         mat(k,2700) = 0._r8
         mat(k,2701) = 0._r8
         mat(k,2702) = 0._r8
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
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 59) = mat(k, 59) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 85) = mat(k, 85) - dti(k)
         mat(k, 91) = mat(k, 91) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 104) = mat(k, 104) - dti(k)
         mat(k, 106) = mat(k, 106) - dti(k)
         mat(k, 109) = mat(k, 109) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
         mat(k, 118) = mat(k, 118) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 126) = mat(k, 126) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 150) = mat(k, 150) - dti(k)
         mat(k, 154) = mat(k, 154) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 163) = mat(k, 163) - dti(k)
         mat(k, 166) = mat(k, 166) - dti(k)
         mat(k, 169) = mat(k, 169) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 184) = mat(k, 184) - dti(k)
         mat(k, 190) = mat(k, 190) - dti(k)
         mat(k, 196) = mat(k, 196) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 205) = mat(k, 205) - dti(k)
         mat(k, 209) = mat(k, 209) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 233) = mat(k, 233) - dti(k)
         mat(k, 237) = mat(k, 237) - dti(k)
         mat(k, 240) = mat(k, 240) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 262) = mat(k, 262) - dti(k)
         mat(k, 266) = mat(k, 266) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 274) = mat(k, 274) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 288) = mat(k, 288) - dti(k)
         mat(k, 294) = mat(k, 294) - dti(k)
         mat(k, 298) = mat(k, 298) - dti(k)
         mat(k, 302) = mat(k, 302) - dti(k)
         mat(k, 306) = mat(k, 306) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 317) = mat(k, 317) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 326) = mat(k, 326) - dti(k)
         mat(k, 332) = mat(k, 332) - dti(k)
         mat(k, 335) = mat(k, 335) - dti(k)
         mat(k, 341) = mat(k, 341) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 353) = mat(k, 353) - dti(k)
         mat(k, 358) = mat(k, 358) - dti(k)
         mat(k, 363) = mat(k, 363) - dti(k)
         mat(k, 368) = mat(k, 368) - dti(k)
         mat(k, 376) = mat(k, 376) - dti(k)
         mat(k, 381) = mat(k, 381) - dti(k)
         mat(k, 386) = mat(k, 386) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 394) = mat(k, 394) - dti(k)
         mat(k, 398) = mat(k, 398) - dti(k)
         mat(k, 405) = mat(k, 405) - dti(k)
         mat(k, 413) = mat(k, 413) - dti(k)
         mat(k, 421) = mat(k, 421) - dti(k)
         mat(k, 424) = mat(k, 424) - dti(k)
         mat(k, 432) = mat(k, 432) - dti(k)
         mat(k, 440) = mat(k, 440) - dti(k)
         mat(k, 448) = mat(k, 448) - dti(k)
         mat(k, 454) = mat(k, 454) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 466) = mat(k, 466) - dti(k)
         mat(k, 472) = mat(k, 472) - dti(k)
         mat(k, 478) = mat(k, 478) - dti(k)
         mat(k, 484) = mat(k, 484) - dti(k)
         mat(k, 490) = mat(k, 490) - dti(k)
         mat(k, 496) = mat(k, 496) - dti(k)
         mat(k, 502) = mat(k, 502) - dti(k)
         mat(k, 510) = mat(k, 510) - dti(k)
         mat(k, 516) = mat(k, 516) - dti(k)
         mat(k, 522) = mat(k, 522) - dti(k)
         mat(k, 529) = mat(k, 529) - dti(k)
         mat(k, 535) = mat(k, 535) - dti(k)
         mat(k, 540) = mat(k, 540) - dti(k)
         mat(k, 543) = mat(k, 543) - dti(k)
         mat(k, 547) = mat(k, 547) - dti(k)
         mat(k, 551) = mat(k, 551) - dti(k)
         mat(k, 554) = mat(k, 554) - dti(k)
         mat(k, 561) = mat(k, 561) - dti(k)
         mat(k, 569) = mat(k, 569) - dti(k)
         mat(k, 578) = mat(k, 578) - dti(k)
         mat(k, 586) = mat(k, 586) - dti(k)
         mat(k, 595) = mat(k, 595) - dti(k)
         mat(k, 600) = mat(k, 600) - dti(k)
         mat(k, 605) = mat(k, 605) - dti(k)
         mat(k, 611) = mat(k, 611) - dti(k)
         mat(k, 615) = mat(k, 615) - dti(k)
         mat(k, 622) = mat(k, 622) - dti(k)
         mat(k, 627) = mat(k, 627) - dti(k)
         mat(k, 635) = mat(k, 635) - dti(k)
         mat(k, 643) = mat(k, 643) - dti(k)
         mat(k, 651) = mat(k, 651) - dti(k)
         mat(k, 659) = mat(k, 659) - dti(k)
         mat(k, 667) = mat(k, 667) - dti(k)
         mat(k, 674) = mat(k, 674) - dti(k)
         mat(k, 683) = mat(k, 683) - dti(k)
         mat(k, 692) = mat(k, 692) - dti(k)
         mat(k, 701) = mat(k, 701) - dti(k)
         mat(k, 710) = mat(k, 710) - dti(k)
         mat(k, 714) = mat(k, 714) - dti(k)
         mat(k, 723) = mat(k, 723) - dti(k)
         mat(k, 731) = mat(k, 731) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 748) = mat(k, 748) - dti(k)
         mat(k, 759) = mat(k, 759) - dti(k)
         mat(k, 770) = mat(k, 770) - dti(k)
         mat(k, 777) = mat(k, 777) - dti(k)
         mat(k, 786) = mat(k, 786) - dti(k)
         mat(k, 799) = mat(k, 799) - dti(k)
         mat(k, 806) = mat(k, 806) - dti(k)
         mat(k, 817) = mat(k, 817) - dti(k)
         mat(k, 828) = mat(k, 828) - dti(k)
         mat(k, 841) = mat(k, 841) - dti(k)
         mat(k, 852) = mat(k, 852) - dti(k)
         mat(k, 861) = mat(k, 861) - dti(k)
         mat(k, 870) = mat(k, 870) - dti(k)
         mat(k, 874) = mat(k, 874) - dti(k)
         mat(k, 880) = mat(k, 880) - dti(k)
         mat(k, 891) = mat(k, 891) - dti(k)
         mat(k, 907) = mat(k, 907) - dti(k)
         mat(k, 915) = mat(k, 915) - dti(k)
         mat(k, 926) = mat(k, 926) - dti(k)
         mat(k, 938) = mat(k, 938) - dti(k)
         mat(k, 949) = mat(k, 949) - dti(k)
         mat(k, 958) = mat(k, 958) - dti(k)
         mat(k, 976) = mat(k, 976) - dti(k)
         mat(k,1004) = mat(k,1004) - dti(k)
         mat(k,1028) = mat(k,1028) - dti(k)
         mat(k,1039) = mat(k,1039) - dti(k)
         mat(k,1047) = mat(k,1047) - dti(k)
         mat(k,1055) = mat(k,1055) - dti(k)
         mat(k,1065) = mat(k,1065) - dti(k)
         mat(k,1073) = mat(k,1073) - dti(k)
         mat(k,1080) = mat(k,1080) - dti(k)
         mat(k,1090) = mat(k,1090) - dti(k)
         mat(k,1103) = mat(k,1103) - dti(k)
         mat(k,1110) = mat(k,1110) - dti(k)
         mat(k,1126) = mat(k,1126) - dti(k)
         mat(k,1147) = mat(k,1147) - dti(k)
         mat(k,1161) = mat(k,1161) - dti(k)
         mat(k,1175) = mat(k,1175) - dti(k)
         mat(k,1187) = mat(k,1187) - dti(k)
         mat(k,1198) = mat(k,1198) - dti(k)
         mat(k,1205) = mat(k,1205) - dti(k)
         mat(k,1217) = mat(k,1217) - dti(k)
         mat(k,1234) = mat(k,1234) - dti(k)
         mat(k,1247) = mat(k,1247) - dti(k)
         mat(k,1255) = mat(k,1255) - dti(k)
         mat(k,1271) = mat(k,1271) - dti(k)
         mat(k,1288) = mat(k,1288) - dti(k)
         mat(k,1308) = mat(k,1308) - dti(k)
         mat(k,1324) = mat(k,1324) - dti(k)
         mat(k,1336) = mat(k,1336) - dti(k)
         mat(k,1356) = mat(k,1356) - dti(k)
         mat(k,1389) = mat(k,1389) - dti(k)
         mat(k,1413) = mat(k,1413) - dti(k)
         mat(k,1434) = mat(k,1434) - dti(k)
         mat(k,1456) = mat(k,1456) - dti(k)
         mat(k,1488) = mat(k,1488) - dti(k)
         mat(k,1504) = mat(k,1504) - dti(k)
         mat(k,1519) = mat(k,1519) - dti(k)
         mat(k,1532) = mat(k,1532) - dti(k)
         mat(k,1546) = mat(k,1546) - dti(k)
         mat(k,1564) = mat(k,1564) - dti(k)
         mat(k,1583) = mat(k,1583) - dti(k)
         mat(k,1637) = mat(k,1637) - dti(k)
         mat(k,1660) = mat(k,1660) - dti(k)
         mat(k,1685) = mat(k,1685) - dti(k)
         mat(k,1712) = mat(k,1712) - dti(k)
         mat(k,1737) = mat(k,1737) - dti(k)
         mat(k,1760) = mat(k,1760) - dti(k)
         mat(k,1818) = mat(k,1818) - dti(k)
         mat(k,1919) = mat(k,1919) - dti(k)
         mat(k,1966) = mat(k,1966) - dti(k)
         mat(k,2149) = mat(k,2149) - dti(k)
         mat(k,2178) = mat(k,2178) - dti(k)
         mat(k,2208) = mat(k,2208) - dti(k)
         mat(k,2324) = mat(k,2324) - dti(k)
         mat(k,2386) = mat(k,2386) - dti(k)
         mat(k,2417) = mat(k,2417) - dti(k)
         mat(k,2446) = mat(k,2446) - dti(k)
         mat(k,2470) = mat(k,2470) - dti(k)
         mat(k,2506) = mat(k,2506) - dti(k)
         mat(k,2535) = mat(k,2535) - dti(k)
         mat(k,2602) = mat(k,2602) - dti(k)
         mat(k,2673) = mat(k,2673) - dti(k)
         mat(k,2703) = mat(k,2703) - dti(k)
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
