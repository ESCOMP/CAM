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
         mat(k,708) = -(rxt(k,416)*y(k,263))
         mat(k,1839) = -rxt(k,416)*y(k,1)
         mat(k,2183) = rxt(k,419)*y(k,232)
         mat(k,1069) = rxt(k,419)*y(k,131)
         mat(k,697) = -(rxt(k,420)*y(k,263))
         mat(k,1838) = -rxt(k,420)*y(k,2)
         mat(k,1068) = rxt(k,417)*y(k,245)
         mat(k,2360) = rxt(k,417)*y(k,232)
         mat(k,1048) = -(rxt(k,499)*y(k,133) + rxt(k,500)*y(k,142) + rxt(k,501) &
                      *y(k,263))
         mat(k,2257) = -rxt(k,499)*y(k,6)
         mat(k,2067) = -rxt(k,500)*y(k,6)
         mat(k,1865) = -rxt(k,501)*y(k,6)
         mat(k,82) = -(rxt(k,554)*y(k,245) + rxt(k,555)*y(k,131))
         mat(k,2316) = -rxt(k,554)*y(k,7)
         mat(k,2150) = -rxt(k,555)*y(k,7)
         mat(k,1044) = rxt(k,557)*y(k,263)
         mat(k,1750) = rxt(k,557)*y(k,6)
         mat(k,204) = -(rxt(k,458)*y(k,263))
         mat(k,1769) = -rxt(k,458)*y(k,8)
         mat(k,113) = -(rxt(k,559)*y(k,245) + rxt(k,560)*y(k,131))
         mat(k,2325) = -rxt(k,559)*y(k,9)
         mat(k,2159) = -rxt(k,560)*y(k,9)
         mat(k,203) = rxt(k,558)*y(k,263)
         mat(k,1760) = rxt(k,558)*y(k,8)
         mat(k,465) = -(rxt(k,461)*y(k,263))
         mat(k,1810) = -rxt(k,461)*y(k,10)
         mat(k,537) = rxt(k,459)*y(k,245)
         mat(k,2340) = rxt(k,459)*y(k,233)
         mat(k,205) = .120_r8*rxt(k,458)*y(k,263)
         mat(k,1770) = .120_r8*rxt(k,458)*y(k,8)
         mat(k,1046) = .100_r8*rxt(k,500)*y(k,142)
         mat(k,1018) = .100_r8*rxt(k,503)*y(k,142)
         mat(k,2056) = .100_r8*rxt(k,500)*y(k,6) + .100_r8*rxt(k,503)*y(k,116)
         mat(k,2170) = .500_r8*rxt(k,460)*y(k,233) + .200_r8*rxt(k,487)*y(k,270) &
                      + .060_r8*rxt(k,493)*y(k,272)
         mat(k,538) = .500_r8*rxt(k,460)*y(k,131)
         mat(k,789) = .200_r8*rxt(k,487)*y(k,131)
         mat(k,805) = .060_r8*rxt(k,493)*y(k,131)
         mat(k,2164) = .200_r8*rxt(k,487)*y(k,270) + .200_r8*rxt(k,493)*y(k,272)
         mat(k,788) = .200_r8*rxt(k,487)*y(k,131)
         mat(k,803) = .200_r8*rxt(k,493)*y(k,131)
         mat(k,2180) = .200_r8*rxt(k,487)*y(k,270) + .150_r8*rxt(k,493)*y(k,272)
         mat(k,791) = .200_r8*rxt(k,487)*y(k,131)
         mat(k,806) = .150_r8*rxt(k,493)*y(k,131)
         mat(k,2166) = .210_r8*rxt(k,493)*y(k,272)
         mat(k,804) = .210_r8*rxt(k,493)*y(k,131)
         mat(k,265) = -(rxt(k,421)*y(k,263))
         mat(k,1779) = -rxt(k,421)*y(k,17)
         mat(k,1045) = .050_r8*rxt(k,500)*y(k,142)
         mat(k,1017) = .050_r8*rxt(k,503)*y(k,142)
         mat(k,2055) = .050_r8*rxt(k,500)*y(k,6) + .050_r8*rxt(k,503)*y(k,116)
         mat(k,407) = -(rxt(k,387)*y(k,133) + rxt(k,388)*y(k,263))
         mat(k,2250) = -rxt(k,387)*y(k,18)
         mat(k,1802) = -rxt(k,388)*y(k,18)
         mat(k,1553) = -(rxt(k,270)*y(k,44) + rxt(k,271)*y(k,245) + rxt(k,272) &
                      *y(k,142))
         mat(k,2432) = -rxt(k,270)*y(k,19)
         mat(k,2406) = -rxt(k,271)*y(k,19)
         mat(k,2093) = -rxt(k,272)*y(k,19)
         mat(k,1605) = 4.000_r8*rxt(k,273)*y(k,21) + (rxt(k,274)+rxt(k,275))*y(k,61) &
                      + rxt(k,278)*y(k,131) + rxt(k,281)*y(k,140) + rxt(k,529) &
                      *y(k,160) + rxt(k,282)*y(k,263)
         mat(k,179) = rxt(k,260)*y(k,259)
         mat(k,185) = rxt(k,286)*y(k,259)
         mat(k,521) = 2.000_r8*rxt(k,297)*y(k,58) + 2.000_r8*rxt(k,309)*y(k,259) &
                      + 2.000_r8*rxt(k,298)*y(k,263)
         mat(k,628) = rxt(k,299)*y(k,58) + rxt(k,310)*y(k,259) + rxt(k,300)*y(k,263)
         mat(k,449) = 3.000_r8*rxt(k,304)*y(k,58) + 3.000_r8*rxt(k,287)*y(k,259) &
                      + 3.000_r8*rxt(k,305)*y(k,263)
         mat(k,1986) = 2.000_r8*rxt(k,297)*y(k,43) + rxt(k,299)*y(k,45) &
                      + 3.000_r8*rxt(k,304)*y(k,57)
         mat(k,2121) = (rxt(k,274)+rxt(k,275))*y(k,21)
         mat(k,153) = 2.000_r8*rxt(k,288)*y(k,259)
         mat(k,871) = rxt(k,283)*y(k,140) + rxt(k,289)*y(k,259) + rxt(k,284)*y(k,263)
         mat(k,2226) = rxt(k,278)*y(k,21)
         mat(k,1940) = rxt(k,281)*y(k,21) + rxt(k,283)*y(k,83)
         mat(k,1519) = rxt(k,529)*y(k,21)
         mat(k,2029) = rxt(k,260)*y(k,36) + rxt(k,286)*y(k,37) + 2.000_r8*rxt(k,309) &
                      *y(k,43) + rxt(k,310)*y(k,45) + 3.000_r8*rxt(k,287)*y(k,57) &
                      + 2.000_r8*rxt(k,288)*y(k,80) + rxt(k,289)*y(k,83)
         mat(k,1897) = rxt(k,282)*y(k,21) + 2.000_r8*rxt(k,298)*y(k,43) + rxt(k,300) &
                      *y(k,45) + 3.000_r8*rxt(k,305)*y(k,57) + rxt(k,284)*y(k,83)
         mat(k,1598) = rxt(k,276)*y(k,61)
         mat(k,2114) = rxt(k,276)*y(k,21)
         mat(k,1533) = (rxt(k,594)+rxt(k,599))*y(k,93)
         mat(k,828) = (rxt(k,594)+rxt(k,599))*y(k,87)
         mat(k,1607) = -(4._r8*rxt(k,273)*y(k,21) + (rxt(k,274) + rxt(k,275) + rxt(k,276) &
                      ) * y(k,61) + rxt(k,277)*y(k,245) + rxt(k,278)*y(k,131) &
                      + rxt(k,279)*y(k,132) + rxt(k,281)*y(k,140) + rxt(k,282) &
                      *y(k,263) + rxt(k,529)*y(k,160))
         mat(k,2123) = -(rxt(k,274) + rxt(k,275) + rxt(k,276)) * y(k,21)
         mat(k,2408) = -rxt(k,277)*y(k,21)
         mat(k,2228) = -rxt(k,278)*y(k,21)
         mat(k,2479) = -rxt(k,279)*y(k,21)
         mat(k,1942) = -rxt(k,281)*y(k,21)
         mat(k,1899) = -rxt(k,282)*y(k,21)
         mat(k,1521) = -rxt(k,529)*y(k,21)
         mat(k,1555) = rxt(k,272)*y(k,142)
         mat(k,589) = rxt(k,280)*y(k,140)
         mat(k,872) = rxt(k,290)*y(k,259)
         mat(k,832) = rxt(k,285)*y(k,140)
         mat(k,1942) = mat(k,1942) + rxt(k,280)*y(k,22) + rxt(k,285)*y(k,93)
         mat(k,2095) = rxt(k,272)*y(k,19)
         mat(k,2031) = rxt(k,290)*y(k,83)
         mat(k,586) = -(rxt(k,280)*y(k,140))
         mat(k,1921) = -rxt(k,280)*y(k,22)
         mat(k,1600) = rxt(k,279)*y(k,132)
         mat(k,2457) = rxt(k,279)*y(k,21)
         mat(k,274) = -(rxt(k,462)*y(k,263))
         mat(k,1781) = -rxt(k,462)*y(k,24)
         mat(k,2163) = rxt(k,465)*y(k,234)
         mat(k,483) = rxt(k,465)*y(k,131)
         mat(k,371) = -(rxt(k,464)*y(k,263))
         mat(k,1796) = -rxt(k,464)*y(k,25)
         mat(k,484) = rxt(k,463)*y(k,245)
         mat(k,2332) = rxt(k,463)*y(k,234)
         mat(k,317) = -(rxt(k,335)*y(k,58) + rxt(k,336)*y(k,263))
         mat(k,1960) = -rxt(k,335)*y(k,26)
         mat(k,1790) = -rxt(k,336)*y(k,26)
         mat(k,578) = -(rxt(k,337)*y(k,58) + rxt(k,338)*y(k,142) + rxt(k,363)*y(k,263))
         mat(k,1966) = -rxt(k,337)*y(k,27)
         mat(k,2058) = -rxt(k,338)*y(k,27)
         mat(k,1824) = -rxt(k,363)*y(k,27)
         mat(k,309) = -(rxt(k,343)*y(k,263))
         mat(k,1788) = -rxt(k,343)*y(k,28)
         mat(k,950) = .800_r8*rxt(k,339)*y(k,235) + .200_r8*rxt(k,340)*y(k,239)
         mat(k,1647) = .200_r8*rxt(k,340)*y(k,235)
         mat(k,376) = -(rxt(k,344)*y(k,263))
         mat(k,1797) = -rxt(k,344)*y(k,29)
         mat(k,951) = rxt(k,341)*y(k,245)
         mat(k,2333) = rxt(k,341)*y(k,235)
         mat(k,323) = -(rxt(k,345)*y(k,58) + rxt(k,346)*y(k,263))
         mat(k,1961) = -rxt(k,345)*y(k,30)
         mat(k,1791) = -rxt(k,346)*y(k,30)
         mat(k,1175) = -(rxt(k,366)*y(k,133) + rxt(k,367)*y(k,142) + rxt(k,385) &
                      *y(k,263))
         mat(k,2266) = -rxt(k,366)*y(k,31)
         mat(k,2075) = -rxt(k,367)*y(k,31)
         mat(k,1875) = -rxt(k,385)*y(k,31)
         mat(k,926) = .130_r8*rxt(k,445)*y(k,142)
         mat(k,2075) = mat(k,2075) + .130_r8*rxt(k,445)*y(k,100)
         mat(k,435) = -(rxt(k,371)*y(k,263))
         mat(k,1805) = -rxt(k,371)*y(k,32)
         mat(k,976) = rxt(k,369)*y(k,245)
         mat(k,2337) = rxt(k,369)*y(k,236)
         mat(k,329) = -(rxt(k,372)*y(k,263) + rxt(k,375)*y(k,58))
         mat(k,1792) = -rxt(k,372)*y(k,33)
         mat(k,1962) = -rxt(k,375)*y(k,33)
         mat(k,313) = -(rxt(k,468)*y(k,263))
         mat(k,1789) = -rxt(k,468)*y(k,34)
         mat(k,688) = rxt(k,466)*y(k,245)
         mat(k,2330) = rxt(k,466)*y(k,237)
         mat(k,142) = -(rxt(k,259)*y(k,259))
         mat(k,2005) = -rxt(k,259)*y(k,35)
         mat(k,177) = -(rxt(k,260)*y(k,259))
         mat(k,2010) = -rxt(k,260)*y(k,36)
         mat(k,182) = -(rxt(k,286)*y(k,259))
         mat(k,2011) = -rxt(k,286)*y(k,37)
         mat(k,155) = -(rxt(k,261)*y(k,259))
         mat(k,2007) = -rxt(k,261)*y(k,38)
         mat(k,187) = -(rxt(k,262)*y(k,259))
         mat(k,2012) = -rxt(k,262)*y(k,39)
         mat(k,159) = -(rxt(k,263)*y(k,259))
         mat(k,2008) = -rxt(k,263)*y(k,40)
         mat(k,192) = -(rxt(k,264)*y(k,259))
         mat(k,2013) = -rxt(k,264)*y(k,41)
         mat(k,163) = -(rxt(k,265)*y(k,259))
         mat(k,2009) = -rxt(k,265)*y(k,42)
         mat(k,519) = -(rxt(k,297)*y(k,58) + rxt(k,298)*y(k,263) + rxt(k,309)*y(k,259))
         mat(k,1965) = -rxt(k,297)*y(k,43)
         mat(k,1817) = -rxt(k,298)*y(k,43)
         mat(k,2023) = -rxt(k,309)*y(k,43)
         mat(k,2448) = -(rxt(k,234)*y(k,58) + rxt(k,270)*y(k,19) + rxt(k,314)*y(k,245) &
                      + rxt(k,315)*y(k,133) + rxt(k,316)*y(k,140) + rxt(k,317) &
                      *y(k,263))
         mat(k,2002) = -rxt(k,234)*y(k,44)
         mat(k,1563) = -rxt(k,270)*y(k,44)
         mat(k,2422) = -rxt(k,314)*y(k,44)
         mat(k,2302) = -rxt(k,315)*y(k,44)
         mat(k,1956) = -rxt(k,316)*y(k,44)
         mat(k,1913) = -rxt(k,317)*y(k,44)
         mat(k,716) = .400_r8*rxt(k,416)*y(k,263)
         mat(k,1065) = .340_r8*rxt(k,500)*y(k,142)
         mat(k,413) = .500_r8*rxt(k,387)*y(k,133)
         mat(k,585) = rxt(k,338)*y(k,142)
         mat(k,1190) = .500_r8*rxt(k,367)*y(k,142)
         mat(k,640) = .500_r8*rxt(k,355)*y(k,263)
         mat(k,858) = rxt(k,322)*y(k,263)
         mat(k,445) = .300_r8*rxt(k,323)*y(k,263)
         mat(k,1644) = (rxt(k,331)+rxt(k,332))*y(k,259)
         mat(k,2136) = rxt(k,241)*y(k,239)
         mat(k,1212) = .800_r8*rxt(k,360)*y(k,263)
         mat(k,939) = .910_r8*rxt(k,445)*y(k,142)
         mat(k,649) = .300_r8*rxt(k,436)*y(k,263)
         mat(k,1307) = .800_r8*rxt(k,440)*y(k,239)
         mat(k,1320) = .120_r8*rxt(k,398)*y(k,142)
         mat(k,669) = .500_r8*rxt(k,411)*y(k,263)
         mat(k,1037) = .340_r8*rxt(k,503)*y(k,142)
         mat(k,1431) = .600_r8*rxt(k,412)*y(k,142)
         mat(k,2242) = .100_r8*rxt(k,418)*y(k,232) + rxt(k,321)*y(k,239) &
                      + .500_r8*rxt(k,389)*y(k,242) + .500_r8*rxt(k,357)*y(k,244) &
                      + .920_r8*rxt(k,428)*y(k,247) + .250_r8*rxt(k,396)*y(k,249) &
                      + rxt(k,405)*y(k,251) + rxt(k,379)*y(k,266) + rxt(k,383) &
                      *y(k,267) + .340_r8*rxt(k,512)*y(k,268) + .320_r8*rxt(k,517) &
                      *y(k,269) + .250_r8*rxt(k,453)*y(k,271)
         mat(k,2302) = mat(k,2302) + .500_r8*rxt(k,387)*y(k,18) + rxt(k,429)*y(k,247) &
                      + .250_r8*rxt(k,395)*y(k,249) + rxt(k,406)*y(k,251)
         mat(k,2109) = .340_r8*rxt(k,500)*y(k,6) + rxt(k,338)*y(k,27) &
                      + .500_r8*rxt(k,367)*y(k,31) + .910_r8*rxt(k,445)*y(k,100) &
                      + .120_r8*rxt(k,398)*y(k,111) + .340_r8*rxt(k,503)*y(k,116) &
                      + .600_r8*rxt(k,412)*y(k,118)
         mat(k,608) = rxt(k,362)*y(k,263)
         mat(k,1153) = .680_r8*rxt(k,521)*y(k,263)
         mat(k,1081) = .100_r8*rxt(k,418)*y(k,131)
         mat(k,960) = .700_r8*rxt(k,340)*y(k,239)
         mat(k,985) = rxt(k,368)*y(k,239)
         mat(k,1481) = rxt(k,351)*y(k,239) + rxt(k,425)*y(k,247) + .250_r8*rxt(k,392) &
                      *y(k,249) + rxt(k,401)*y(k,251) + .250_r8*rxt(k,450)*y(k,271)
         mat(k,1696) = rxt(k,241)*y(k,61) + .800_r8*rxt(k,440)*y(k,103) + rxt(k,321) &
                      *y(k,131) + .700_r8*rxt(k,340)*y(k,235) + rxt(k,368)*y(k,236) &
                      + rxt(k,351)*y(k,238) + (4.000_r8*rxt(k,318)+2.000_r8*rxt(k,319)) &
                      *y(k,239) + 1.500_r8*rxt(k,426)*y(k,247) + .750_r8*rxt(k,431) &
                      *y(k,248) + .880_r8*rxt(k,393)*y(k,249) + 2.000_r8*rxt(k,402) &
                      *y(k,251) + .750_r8*rxt(k,505)*y(k,258) + .800_r8*rxt(k,381) &
                      *y(k,267) + .930_r8*rxt(k,510)*y(k,268) + .950_r8*rxt(k,515) &
                      *y(k,269) + .800_r8*rxt(k,451)*y(k,271)
         mat(k,624) = .500_r8*rxt(k,389)*y(k,131)
         mat(k,845) = .500_r8*rxt(k,357)*y(k,131)
         mat(k,2422) = mat(k,2422) + .450_r8*rxt(k,403)*y(k,251) + .150_r8*rxt(k,382) &
                      *y(k,267)
         mat(k,1354) = .920_r8*rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133) + rxt(k,425) &
                      *y(k,238) + 1.500_r8*rxt(k,426)*y(k,239)
         mat(k,1387) = .750_r8*rxt(k,431)*y(k,239)
         mat(k,1408) = .250_r8*rxt(k,396)*y(k,131) + .250_r8*rxt(k,395)*y(k,133) &
                      + .250_r8*rxt(k,392)*y(k,238) + .880_r8*rxt(k,393)*y(k,239)
         mat(k,1449) = rxt(k,405)*y(k,131) + rxt(k,406)*y(k,133) + rxt(k,401)*y(k,238) &
                      + 2.000_r8*rxt(k,402)*y(k,239) + .450_r8*rxt(k,403)*y(k,245) &
                      + 4.000_r8*rxt(k,404)*y(k,251)
         mat(k,1142) = .750_r8*rxt(k,505)*y(k,239)
         mat(k,2045) = (rxt(k,331)+rxt(k,332))*y(k,56)
         mat(k,1913) = mat(k,1913) + .400_r8*rxt(k,416)*y(k,1) + .500_r8*rxt(k,355) &
                      *y(k,53) + rxt(k,322)*y(k,54) + .300_r8*rxt(k,323)*y(k,55) &
                      + .800_r8*rxt(k,360)*y(k,76) + .300_r8*rxt(k,436)*y(k,101) &
                      + .500_r8*rxt(k,411)*y(k,115) + rxt(k,362)*y(k,147) &
                      + .680_r8*rxt(k,521)*y(k,219)
         mat(k,867) = rxt(k,379)*y(k,131)
         mat(k,1266) = rxt(k,383)*y(k,131) + .800_r8*rxt(k,381)*y(k,239) &
                      + .150_r8*rxt(k,382)*y(k,245)
         mat(k,1229) = .340_r8*rxt(k,512)*y(k,131) + .930_r8*rxt(k,510)*y(k,239)
         mat(k,1100) = .320_r8*rxt(k,517)*y(k,131) + .950_r8*rxt(k,515)*y(k,239)
         mat(k,1284) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,450)*y(k,238) &
                      + .800_r8*rxt(k,451)*y(k,239)
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
         mat(k,626) = -(rxt(k,299)*y(k,58) + rxt(k,300)*y(k,263) + rxt(k,310)*y(k,259))
         mat(k,1967) = -rxt(k,299)*y(k,45)
         mat(k,1829) = -rxt(k,300)*y(k,45)
         mat(k,2024) = -rxt(k,310)*y(k,45)
         mat(k,167) = -(rxt(k,301)*y(k,263))
         mat(k,1767) = -rxt(k,301)*y(k,46)
         mat(k,1193) = -(rxt(k,347)*y(k,133) + rxt(k,348)*y(k,263))
         mat(k,2267) = -rxt(k,347)*y(k,47)
         mat(k,1876) = -rxt(k,348)*y(k,47)
         mat(k,712) = .800_r8*rxt(k,416)*y(k,263)
         mat(k,410) = rxt(k,387)*y(k,133)
         mat(k,310) = rxt(k,343)*y(k,263)
         mat(k,378) = .500_r8*rxt(k,344)*y(k,263)
         mat(k,1176) = .500_r8*rxt(k,367)*y(k,142)
         mat(k,1413) = .100_r8*rxt(k,412)*y(k,142)
         mat(k,2208) = .400_r8*rxt(k,418)*y(k,232) + rxt(k,342)*y(k,235) &
                      + .270_r8*rxt(k,370)*y(k,236) + rxt(k,389)*y(k,242) + rxt(k,408) &
                      *y(k,253) + rxt(k,379)*y(k,266)
         mat(k,2267) = mat(k,2267) + rxt(k,387)*y(k,18)
         mat(k,2076) = .500_r8*rxt(k,367)*y(k,31) + .100_r8*rxt(k,412)*y(k,118)
         mat(k,1074) = .400_r8*rxt(k,418)*y(k,131)
         mat(k,954) = rxt(k,342)*y(k,131) + 3.200_r8*rxt(k,339)*y(k,235) &
                      + .800_r8*rxt(k,340)*y(k,239)
         mat(k,979) = .270_r8*rxt(k,370)*y(k,131)
         mat(k,1664) = .800_r8*rxt(k,340)*y(k,235)
         mat(k,620) = rxt(k,389)*y(k,131)
         mat(k,2387) = .200_r8*rxt(k,407)*y(k,253)
         mat(k,720) = rxt(k,408)*y(k,131) + .200_r8*rxt(k,407)*y(k,245)
         mat(k,1876) = mat(k,1876) + .800_r8*rxt(k,416)*y(k,1) + rxt(k,343)*y(k,28) &
                      + .500_r8*rxt(k,344)*y(k,29)
         mat(k,861) = rxt(k,379)*y(k,131)
         mat(k,399) = -(rxt(k,302)*y(k,58) + rxt(k,303)*y(k,263))
         mat(k,1963) = -rxt(k,302)*y(k,48)
         mat(k,1801) = -rxt(k,303)*y(k,48)
         mat(k,145) = -(rxt(k,349)*y(k,263))
         mat(k,1766) = -rxt(k,349)*y(k,49)
         mat(k,1110) = -(rxt(k,386)*y(k,263))
         mat(k,1870) = -rxt(k,386)*y(k,50)
         mat(k,711) = .800_r8*rxt(k,416)*y(k,263)
         mat(k,1053) = .520_r8*rxt(k,500)*y(k,142)
         mat(k,409) = .500_r8*rxt(k,387)*y(k,133)
         mat(k,1025) = .520_r8*rxt(k,503)*y(k,142)
         mat(k,2204) = .250_r8*rxt(k,418)*y(k,232) + .820_r8*rxt(k,370)*y(k,236) &
                      + .500_r8*rxt(k,389)*y(k,242) + .270_r8*rxt(k,512)*y(k,268) &
                      + .040_r8*rxt(k,517)*y(k,269)
         mat(k,2262) = .500_r8*rxt(k,387)*y(k,18)
         mat(k,2072) = .520_r8*rxt(k,500)*y(k,6) + .520_r8*rxt(k,503)*y(k,116)
         mat(k,1144) = .500_r8*rxt(k,521)*y(k,263)
         mat(k,1073) = .250_r8*rxt(k,418)*y(k,131)
         mat(k,978) = .820_r8*rxt(k,370)*y(k,131) + .820_r8*rxt(k,368)*y(k,239)
         mat(k,1660) = .820_r8*rxt(k,368)*y(k,236) + .150_r8*rxt(k,510)*y(k,268) &
                      + .025_r8*rxt(k,515)*y(k,269)
         mat(k,619) = .500_r8*rxt(k,389)*y(k,131)
         mat(k,1870) = mat(k,1870) + .800_r8*rxt(k,416)*y(k,1) + .500_r8*rxt(k,521) &
                      *y(k,219)
         mat(k,1216) = .270_r8*rxt(k,512)*y(k,131) + .150_r8*rxt(k,510)*y(k,239)
         mat(k,1094) = .040_r8*rxt(k,517)*y(k,131) + .025_r8*rxt(k,515)*y(k,239)
         mat(k,1323) = -(rxt(k,373)*y(k,133) + rxt(k,374)*y(k,263))
         mat(k,2277) = -rxt(k,373)*y(k,51)
         mat(k,1886) = -rxt(k,374)*y(k,51)
         mat(k,1251) = rxt(k,376)*y(k,263)
         mat(k,1312) = .880_r8*rxt(k,398)*y(k,142)
         mat(k,1416) = .500_r8*rxt(k,412)*y(k,142)
         mat(k,2218) = .170_r8*rxt(k,471)*y(k,240) + .050_r8*rxt(k,434)*y(k,248) &
                      + .250_r8*rxt(k,396)*y(k,249) + .170_r8*rxt(k,477)*y(k,252) &
                      + .400_r8*rxt(k,487)*y(k,270) + .250_r8*rxt(k,453)*y(k,271) &
                      + .540_r8*rxt(k,493)*y(k,272) + .510_r8*rxt(k,496)*y(k,273)
         mat(k,2277) = mat(k,2277) + .050_r8*rxt(k,435)*y(k,248) + .250_r8*rxt(k,395) &
                      *y(k,249) + .250_r8*rxt(k,454)*y(k,271)
         mat(k,916) = rxt(k,377)*y(k,263)
         mat(k,2084) = .880_r8*rxt(k,398)*y(k,111) + .500_r8*rxt(k,412)*y(k,118)
         mat(k,1464) = .250_r8*rxt(k,392)*y(k,249) + .250_r8*rxt(k,450)*y(k,271)
         mat(k,1673) = .240_r8*rxt(k,393)*y(k,249) + .500_r8*rxt(k,381)*y(k,267) &
                      + .100_r8*rxt(k,451)*y(k,271)
         mat(k,822) = .170_r8*rxt(k,471)*y(k,131) + .070_r8*rxt(k,470)*y(k,245)
         mat(k,2396) = .070_r8*rxt(k,470)*y(k,240) + .070_r8*rxt(k,476)*y(k,252)
         mat(k,1373) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,1397) = .250_r8*rxt(k,396)*y(k,131) + .250_r8*rxt(k,395)*y(k,133) &
                      + .250_r8*rxt(k,392)*y(k,238) + .240_r8*rxt(k,393)*y(k,239)
         mat(k,970) = .170_r8*rxt(k,477)*y(k,131) + .070_r8*rxt(k,476)*y(k,245)
         mat(k,1886) = mat(k,1886) + rxt(k,376)*y(k,97) + rxt(k,377)*y(k,134)
         mat(k,1260) = .500_r8*rxt(k,381)*y(k,239)
         mat(k,798) = .400_r8*rxt(k,487)*y(k,131)
         mat(k,1276) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,238) + .100_r8*rxt(k,451)*y(k,239)
         mat(k,814) = .540_r8*rxt(k,493)*y(k,131)
         mat(k,556) = .510_r8*rxt(k,496)*y(k,131)
         mat(k,747) = -(rxt(k,354)*y(k,263))
         mat(k,1842) = -rxt(k,354)*y(k,52)
         mat(k,1171) = .120_r8*rxt(k,367)*y(k,142)
         mat(k,2060) = .120_r8*rxt(k,367)*y(k,31)
         mat(k,1455) = .100_r8*rxt(k,351)*y(k,239) + .150_r8*rxt(k,352)*y(k,245)
         mat(k,1652) = .100_r8*rxt(k,351)*y(k,238)
         mat(k,2363) = .150_r8*rxt(k,352)*y(k,238) + .150_r8*rxt(k,403)*y(k,251)
         mat(k,1436) = .150_r8*rxt(k,403)*y(k,245)
         mat(k,635) = -(rxt(k,355)*y(k,263))
         mat(k,1830) = -rxt(k,355)*y(k,53)
         mat(k,1454) = .400_r8*rxt(k,352)*y(k,245)
         mat(k,2355) = .400_r8*rxt(k,352)*y(k,238) + .400_r8*rxt(k,403)*y(k,251)
         mat(k,1434) = .400_r8*rxt(k,403)*y(k,245)
         mat(k,855) = -(rxt(k,322)*y(k,263))
         mat(k,1851) = -rxt(k,322)*y(k,54)
         mat(k,1289) = .200_r8*rxt(k,440)*y(k,239)
         mat(k,952) = .300_r8*rxt(k,340)*y(k,239)
         mat(k,1653) = .200_r8*rxt(k,440)*y(k,103) + .300_r8*rxt(k,340)*y(k,235) &
                      + 2.000_r8*rxt(k,319)*y(k,239) + .250_r8*rxt(k,426)*y(k,247) &
                      + .250_r8*rxt(k,431)*y(k,248) + .250_r8*rxt(k,393)*y(k,249) &
                      + .250_r8*rxt(k,505)*y(k,258) + .500_r8*rxt(k,381)*y(k,267) &
                      + .250_r8*rxt(k,510)*y(k,268) + .250_r8*rxt(k,515)*y(k,269) &
                      + .300_r8*rxt(k,451)*y(k,271)
         mat(k,1333) = .250_r8*rxt(k,426)*y(k,239)
         mat(k,1362) = .250_r8*rxt(k,431)*y(k,239)
         mat(k,1391) = .250_r8*rxt(k,393)*y(k,239)
         mat(k,1130) = .250_r8*rxt(k,505)*y(k,239)
         mat(k,1257) = .500_r8*rxt(k,381)*y(k,239)
         mat(k,1215) = .250_r8*rxt(k,510)*y(k,239)
         mat(k,1091) = .250_r8*rxt(k,515)*y(k,239)
         mat(k,1270) = .300_r8*rxt(k,451)*y(k,239)
         mat(k,441) = -(rxt(k,323)*y(k,263))
         mat(k,1806) = -rxt(k,323)*y(k,55)
         mat(k,1649) = rxt(k,320)*y(k,245)
         mat(k,2338) = rxt(k,320)*y(k,239)
         mat(k,1631) = -(rxt(k,235)*y(k,58) + rxt(k,291)*y(k,75) + rxt(k,324)*y(k,263) &
                      + (rxt(k,330) + rxt(k,331) + rxt(k,332)) * y(k,259))
         mat(k,1989) = -rxt(k,235)*y(k,56)
         mat(k,943) = -rxt(k,291)*y(k,56)
         mat(k,1900) = -rxt(k,324)*y(k,56)
         mat(k,2032) = -(rxt(k,330) + rxt(k,331) + rxt(k,332)) * y(k,56)
         mat(k,1183) = .100_r8*rxt(k,367)*y(k,142)
         mat(k,2096) = .100_r8*rxt(k,367)*y(k,31)
         mat(k,447) = -(rxt(k,287)*y(k,259) + rxt(k,304)*y(k,58) + rxt(k,305)*y(k,263))
         mat(k,2022) = -rxt(k,287)*y(k,57)
         mat(k,1964) = -rxt(k,304)*y(k,57)
         mat(k,1807) = -rxt(k,305)*y(k,57)
         mat(k,1995) = -(rxt(k,234)*y(k,44) + rxt(k,235)*y(k,56) + rxt(k,236)*y(k,79) &
                      + rxt(k,237)*y(k,81) + (rxt(k,238) + rxt(k,239)) * y(k,245) &
                      + rxt(k,240)*y(k,142) + rxt(k,247)*y(k,62) + rxt(k,256)*y(k,94) &
                      + rxt(k,297)*y(k,43) + rxt(k,299)*y(k,45) + rxt(k,302)*y(k,48) &
                      + rxt(k,304)*y(k,57) + rxt(k,345)*y(k,30) + rxt(k,375)*y(k,33))
         mat(k,2441) = -rxt(k,234)*y(k,58)
         mat(k,1637) = -rxt(k,235)*y(k,58)
         mat(k,1511) = -rxt(k,236)*y(k,58)
         mat(k,675) = -rxt(k,237)*y(k,58)
         mat(k,2415) = -(rxt(k,238) + rxt(k,239)) * y(k,58)
         mat(k,2102) = -rxt(k,240)*y(k,58)
         mat(k,1123) = -rxt(k,247)*y(k,58)
         mat(k,883) = -rxt(k,256)*y(k,58)
         mat(k,523) = -rxt(k,297)*y(k,58)
         mat(k,631) = -rxt(k,299)*y(k,58)
         mat(k,404) = -rxt(k,302)*y(k,58)
         mat(k,451) = -rxt(k,304)*y(k,58)
         mat(k,327) = -rxt(k,345)*y(k,58)
         mat(k,333) = -rxt(k,375)*y(k,58)
         mat(k,1612) = rxt(k,275)*y(k,61)
         mat(k,143) = 4.000_r8*rxt(k,259)*y(k,259)
         mat(k,180) = rxt(k,260)*y(k,259)
         mat(k,157) = 2.000_r8*rxt(k,261)*y(k,259)
         mat(k,190) = 2.000_r8*rxt(k,262)*y(k,259)
         mat(k,161) = 2.000_r8*rxt(k,263)*y(k,259)
         mat(k,195) = rxt(k,264)*y(k,259)
         mat(k,165) = 2.000_r8*rxt(k,265)*y(k,259)
         mat(k,169) = 3.000_r8*rxt(k,301)*y(k,263)
         mat(k,404) = mat(k,404) + rxt(k,303)*y(k,263)
         mat(k,2129) = rxt(k,275)*y(k,21) + (4.000_r8*rxt(k,242)+2.000_r8*rxt(k,244)) &
                      *y(k,61) + rxt(k,246)*y(k,131) + rxt(k,251)*y(k,140) &
                      + rxt(k,530)*y(k,160) + rxt(k,241)*y(k,239) + rxt(k,252) &
                      *y(k,263)
         mat(k,292) = rxt(k,296)*y(k,259)
         mat(k,288) = rxt(k,311)*y(k,259) + rxt(k,306)*y(k,263)
         mat(k,298) = rxt(k,312)*y(k,259) + rxt(k,307)*y(k,263)
         mat(k,357) = rxt(k,313)*y(k,259) + rxt(k,308)*y(k,263)
         mat(k,1544) = rxt(k,254)*y(k,140) + rxt(k,266)*y(k,259) + rxt(k,255)*y(k,263)
         mat(k,2235) = rxt(k,246)*y(k,61)
         mat(k,1949) = rxt(k,251)*y(k,61) + rxt(k,254)*y(k,87)
         mat(k,1525) = rxt(k,530)*y(k,61)
         mat(k,1689) = rxt(k,241)*y(k,61)
         mat(k,2038) = 4.000_r8*rxt(k,259)*y(k,35) + rxt(k,260)*y(k,36) &
                      + 2.000_r8*rxt(k,261)*y(k,38) + 2.000_r8*rxt(k,262)*y(k,39) &
                      + 2.000_r8*rxt(k,263)*y(k,40) + rxt(k,264)*y(k,41) &
                      + 2.000_r8*rxt(k,265)*y(k,42) + rxt(k,296)*y(k,67) + rxt(k,311) &
                      *y(k,84) + rxt(k,312)*y(k,85) + rxt(k,313)*y(k,86) + rxt(k,266) &
                      *y(k,87)
         mat(k,1906) = 3.000_r8*rxt(k,301)*y(k,46) + rxt(k,303)*y(k,48) + rxt(k,252) &
                      *y(k,61) + rxt(k,306)*y(k,84) + rxt(k,307)*y(k,85) + rxt(k,308) &
                      *y(k,86) + rxt(k,255)*y(k,87)
         mat(k,1959) = rxt(k,247)*y(k,62)
         mat(k,2113) = 2.000_r8*rxt(k,243)*y(k,61)
         mat(k,1116) = rxt(k,247)*y(k,58) + (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,87)
         mat(k,1532) = (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,62) + (rxt(k,587) &
                       +rxt(k,593)+rxt(k,598))*y(k,94)
         mat(k,878) = (rxt(k,587)+rxt(k,593)+rxt(k,598))*y(k,87)
         mat(k,2112) = 2.000_r8*rxt(k,268)*y(k,61)
         mat(k,2132) = -(rxt(k,241)*y(k,239) + (4._r8*rxt(k,242) + 4._r8*rxt(k,243) &
                      + 4._r8*rxt(k,244) + 4._r8*rxt(k,268)) * y(k,61) + rxt(k,245) &
                      *y(k,245) + rxt(k,246)*y(k,131) + rxt(k,248)*y(k,132) + rxt(k,251) &
                      *y(k,140) + (rxt(k,252) + rxt(k,253)) * y(k,263) + (rxt(k,274) &
                      + rxt(k,275) + rxt(k,276)) * y(k,21) + rxt(k,530)*y(k,160))
         mat(k,1692) = -rxt(k,241)*y(k,61)
         mat(k,2418) = -rxt(k,245)*y(k,61)
         mat(k,2238) = -rxt(k,246)*y(k,61)
         mat(k,2489) = -rxt(k,248)*y(k,61)
         mat(k,1952) = -rxt(k,251)*y(k,61)
         mat(k,1909) = -(rxt(k,252) + rxt(k,253)) * y(k,61)
         mat(k,1615) = -(rxt(k,274) + rxt(k,275) + rxt(k,276)) * y(k,61)
         mat(k,1527) = -rxt(k,530)*y(k,61)
         mat(k,1998) = rxt(k,256)*y(k,94) + rxt(k,240)*y(k,142) + rxt(k,239)*y(k,245)
         mat(k,1124) = rxt(k,249)*y(k,140)
         mat(k,1546) = rxt(k,267)*y(k,259)
         mat(k,884) = rxt(k,256)*y(k,58) + rxt(k,257)*y(k,140) + rxt(k,258)*y(k,263)
         mat(k,1952) = mat(k,1952) + rxt(k,249)*y(k,62) + rxt(k,257)*y(k,94)
         mat(k,2105) = rxt(k,240)*y(k,58)
         mat(k,364) = rxt(k,535)*y(k,160)
         mat(k,1527) = mat(k,1527) + rxt(k,535)*y(k,144)
         mat(k,2418) = mat(k,2418) + rxt(k,239)*y(k,58)
         mat(k,2041) = rxt(k,267)*y(k,87)
         mat(k,1909) = mat(k,1909) + rxt(k,258)*y(k,94)
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
         mat(k,1118) = -(rxt(k,247)*y(k,58) + rxt(k,249)*y(k,140) + rxt(k,250) &
                      *y(k,263) + (rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,87))
         mat(k,1976) = -rxt(k,247)*y(k,62)
         mat(k,1934) = -rxt(k,249)*y(k,62)
         mat(k,1871) = -rxt(k,250)*y(k,62)
         mat(k,1536) = -(rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,62)
         mat(k,2118) = rxt(k,248)*y(k,132)
         mat(k,2468) = rxt(k,248)*y(k,61)
         mat(k,1202) = -(rxt(k,334)*y(k,263))
         mat(k,1877) = -rxt(k,334)*y(k,64)
         mat(k,1056) = .230_r8*rxt(k,500)*y(k,142)
         mat(k,1551) = rxt(k,270)*y(k,44)
         mat(k,320) = .350_r8*rxt(k,336)*y(k,263)
         mat(k,581) = .630_r8*rxt(k,338)*y(k,142)
         mat(k,1177) = .560_r8*rxt(k,367)*y(k,142)
         mat(k,2428) = rxt(k,270)*y(k,19) + rxt(k,234)*y(k,58) + rxt(k,315)*y(k,133) &
                      + rxt(k,316)*y(k,140) + rxt(k,317)*y(k,263)
         mat(k,400) = rxt(k,302)*y(k,58)
         mat(k,1322) = rxt(k,373)*y(k,133) + rxt(k,374)*y(k,263)
         mat(k,1978) = rxt(k,234)*y(k,44) + rxt(k,302)*y(k,48)
         mat(k,1491) = rxt(k,615)*y(k,264)
         mat(k,1085) = rxt(k,361)*y(k,263)
         mat(k,927) = .620_r8*rxt(k,445)*y(k,142)
         mat(k,1310) = .650_r8*rxt(k,398)*y(k,142)
         mat(k,1028) = .230_r8*rxt(k,503)*y(k,142)
         mat(k,1414) = .560_r8*rxt(k,412)*y(k,142)
         mat(k,2209) = .170_r8*rxt(k,471)*y(k,240) + .220_r8*rxt(k,396)*y(k,249) &
                      + .400_r8*rxt(k,474)*y(k,250) + .350_r8*rxt(k,477)*y(k,252) &
                      + .225_r8*rxt(k,512)*y(k,268) + .250_r8*rxt(k,453)*y(k,271)
         mat(k,2268) = rxt(k,315)*y(k,44) + rxt(k,373)*y(k,51) + .220_r8*rxt(k,395) &
                      *y(k,249) + .500_r8*rxt(k,454)*y(k,271)
         mat(k,1935) = rxt(k,316)*y(k,44) + rxt(k,524)*y(k,145)
         mat(k,2077) = .230_r8*rxt(k,500)*y(k,6) + .630_r8*rxt(k,338)*y(k,27) &
                      + .560_r8*rxt(k,367)*y(k,31) + .620_r8*rxt(k,445)*y(k,100) &
                      + .650_r8*rxt(k,398)*y(k,111) + .230_r8*rxt(k,503)*y(k,116) &
                      + .560_r8*rxt(k,412)*y(k,118)
         mat(k,418) = rxt(k,524)*y(k,140) + rxt(k,525)*y(k,263)
         mat(k,1146) = .700_r8*rxt(k,521)*y(k,263)
         mat(k,1458) = .220_r8*rxt(k,392)*y(k,249) + .250_r8*rxt(k,450)*y(k,271)
         mat(k,1665) = .110_r8*rxt(k,393)*y(k,249) + .125_r8*rxt(k,510)*y(k,268) &
                      + .200_r8*rxt(k,451)*y(k,271)
         mat(k,821) = .170_r8*rxt(k,471)*y(k,131) + .070_r8*rxt(k,470)*y(k,245)
         mat(k,2388) = .070_r8*rxt(k,470)*y(k,240) + .160_r8*rxt(k,473)*y(k,250) &
                      + .140_r8*rxt(k,476)*y(k,252)
         mat(k,1392) = .220_r8*rxt(k,396)*y(k,131) + .220_r8*rxt(k,395)*y(k,133) &
                      + .220_r8*rxt(k,392)*y(k,238) + .110_r8*rxt(k,393)*y(k,239)
         mat(k,784) = .400_r8*rxt(k,474)*y(k,131) + .160_r8*rxt(k,473)*y(k,245)
         mat(k,969) = .350_r8*rxt(k,477)*y(k,131) + .140_r8*rxt(k,476)*y(k,245)
         mat(k,1877) = mat(k,1877) + .350_r8*rxt(k,336)*y(k,26) + rxt(k,317)*y(k,44) &
                      + rxt(k,374)*y(k,51) + rxt(k,361)*y(k,77) + rxt(k,525)*y(k,145) &
                      + .700_r8*rxt(k,521)*y(k,219)
         mat(k,851) = rxt(k,615)*y(k,65)
         mat(k,1218) = .225_r8*rxt(k,512)*y(k,131) + .125_r8*rxt(k,510)*y(k,239)
         mat(k,1272) = .250_r8*rxt(k,453)*y(k,131) + .500_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,238) + .200_r8*rxt(k,451)*y(k,239)
         mat(k,1492) = -(rxt(k,615)*y(k,264))
         mat(k,852) = -rxt(k,615)*y(k,65)
         mat(k,1060) = .270_r8*rxt(k,500)*y(k,142)
         mat(k,1181) = .200_r8*rxt(k,367)*y(k,142)
         mat(k,748) = rxt(k,354)*y(k,263)
         mat(k,637) = .500_r8*rxt(k,355)*y(k,263)
         mat(k,1203) = rxt(k,334)*y(k,263)
         mat(k,1209) = .800_r8*rxt(k,360)*y(k,263)
         mat(k,1086) = rxt(k,361)*y(k,263)
         mat(k,963) = rxt(k,326)*y(k,263)
         mat(k,665) = .500_r8*rxt(k,411)*y(k,263)
         mat(k,1032) = .270_r8*rxt(k,503)*y(k,142)
         mat(k,1421) = .100_r8*rxt(k,412)*y(k,142)
         mat(k,2225) = rxt(k,353)*y(k,238) + .900_r8*rxt(k,512)*y(k,268)
         mat(k,2091) = .270_r8*rxt(k,500)*y(k,6) + .200_r8*rxt(k,367)*y(k,31) &
                      + .270_r8*rxt(k,503)*y(k,116) + .100_r8*rxt(k,412)*y(k,118)
         mat(k,1149) = 1.800_r8*rxt(k,521)*y(k,263)
         mat(k,1471) = rxt(k,353)*y(k,131) + 4.000_r8*rxt(k,350)*y(k,238) &
                      + .900_r8*rxt(k,351)*y(k,239) + rxt(k,425)*y(k,247) &
                      + 2.000_r8*rxt(k,401)*y(k,251) + rxt(k,450)*y(k,271)
         mat(k,1680) = .900_r8*rxt(k,351)*y(k,238) + rxt(k,402)*y(k,251) &
                      + .500_r8*rxt(k,510)*y(k,268)
         mat(k,2403) = .450_r8*rxt(k,403)*y(k,251)
         mat(k,1346) = rxt(k,425)*y(k,238)
         mat(k,1441) = 2.000_r8*rxt(k,401)*y(k,238) + rxt(k,402)*y(k,239) &
                      + .450_r8*rxt(k,403)*y(k,245) + 4.000_r8*rxt(k,404)*y(k,251)
         mat(k,1893) = rxt(k,354)*y(k,52) + .500_r8*rxt(k,355)*y(k,53) + rxt(k,334) &
                      *y(k,64) + .800_r8*rxt(k,360)*y(k,76) + rxt(k,361)*y(k,77) &
                      + rxt(k,326)*y(k,89) + .500_r8*rxt(k,411)*y(k,115) &
                      + 1.800_r8*rxt(k,521)*y(k,219)
         mat(k,1223) = .900_r8*rxt(k,512)*y(k,131) + .500_r8*rxt(k,510)*y(k,239)
         mat(k,1278) = rxt(k,450)*y(k,238)
         mat(k,282) = -(rxt(k,295)*y(k,259))
         mat(k,2016) = -rxt(k,295)*y(k,66)
         mat(k,178) = rxt(k,260)*y(k,259)
         mat(k,183) = rxt(k,286)*y(k,259)
         mat(k,188) = rxt(k,262)*y(k,259)
         mat(k,160) = 2.000_r8*rxt(k,263)*y(k,259)
         mat(k,193) = 2.000_r8*rxt(k,264)*y(k,259)
         mat(k,164) = rxt(k,265)*y(k,259)
         mat(k,152) = 2.000_r8*rxt(k,288)*y(k,259)
         mat(k,294) = rxt(k,312)*y(k,259) + rxt(k,307)*y(k,263)
         mat(k,353) = rxt(k,313)*y(k,259) + rxt(k,308)*y(k,263)
         mat(k,2016) = mat(k,2016) + rxt(k,260)*y(k,36) + rxt(k,286)*y(k,37) &
                      + rxt(k,262)*y(k,39) + 2.000_r8*rxt(k,263)*y(k,40) &
                      + 2.000_r8*rxt(k,264)*y(k,41) + rxt(k,265)*y(k,42) &
                      + 2.000_r8*rxt(k,288)*y(k,80) + rxt(k,312)*y(k,85) + rxt(k,313) &
                      *y(k,86)
         mat(k,1783) = rxt(k,307)*y(k,85) + rxt(k,308)*y(k,86)
         mat(k,290) = -(rxt(k,296)*y(k,259))
         mat(k,2018) = -rxt(k,296)*y(k,67)
         mat(k,156) = rxt(k,261)*y(k,259)
         mat(k,189) = rxt(k,262)*y(k,259)
         mat(k,286) = rxt(k,311)*y(k,259) + rxt(k,306)*y(k,263)
         mat(k,2018) = mat(k,2018) + rxt(k,261)*y(k,38) + rxt(k,262)*y(k,39) &
                      + rxt(k,311)*y(k,84)
         mat(k,1785) = rxt(k,306)*y(k,84)
         mat(k,238) = -(rxt(k,469)*y(k,263))
         mat(k,1774) = -rxt(k,469)*y(k,68)
         mat(k,232) = .180_r8*rxt(k,489)*y(k,263)
         mat(k,1774) = mat(k,1774) + .180_r8*rxt(k,489)*y(k,221)
         mat(k,347) = -(rxt(k,522)*y(k,133) + (rxt(k,523) + rxt(k,537)) * y(k,263))
         mat(k,2248) = -rxt(k,522)*y(k,69)
         mat(k,1793) = -(rxt(k,523) + rxt(k,537)) * y(k,69)
         mat(k,837) = rxt(k,356)*y(k,245)
         mat(k,2328) = rxt(k,356)*y(k,244)
         mat(k,941) = -(rxt(k,291)*y(k,56) + rxt(k,292)*y(k,79) + rxt(k,293)*y(k,274) &
                      + rxt(k,294)*y(k,91))
         mat(k,1623) = -rxt(k,291)*y(k,75)
         mat(k,1502) = -rxt(k,292)*y(k,75)
         mat(k,2499) = -rxt(k,293)*y(k,75)
         mat(k,1700) = -rxt(k,294)*y(k,75)
         mat(k,184) = rxt(k,286)*y(k,259)
         mat(k,194) = rxt(k,264)*y(k,259)
         mat(k,283) = 2.000_r8*rxt(k,295)*y(k,259)
         mat(k,291) = rxt(k,296)*y(k,259)
         mat(k,2026) = rxt(k,286)*y(k,37) + rxt(k,264)*y(k,41) + 2.000_r8*rxt(k,295) &
                      *y(k,66) + rxt(k,296)*y(k,67)
         mat(k,1208) = -(rxt(k,360)*y(k,263))
         mat(k,1878) = -rxt(k,360)*y(k,76)
         mat(k,643) = .700_r8*rxt(k,436)*y(k,263)
         mat(k,612) = .500_r8*rxt(k,437)*y(k,263)
         mat(k,461) = rxt(k,448)*y(k,263)
         mat(k,2210) = .050_r8*rxt(k,434)*y(k,248) + .530_r8*rxt(k,396)*y(k,249) &
                      + .225_r8*rxt(k,512)*y(k,268) + .250_r8*rxt(k,453)*y(k,271)
         mat(k,2269) = .050_r8*rxt(k,435)*y(k,248) + .530_r8*rxt(k,395)*y(k,249) &
                      + .250_r8*rxt(k,454)*y(k,271)
         mat(k,1580) = rxt(k,359)*y(k,243)
         mat(k,1459) = .530_r8*rxt(k,392)*y(k,249) + .250_r8*rxt(k,450)*y(k,271)
         mat(k,1666) = .260_r8*rxt(k,393)*y(k,249) + .125_r8*rxt(k,510)*y(k,268) &
                      + .100_r8*rxt(k,451)*y(k,271)
         mat(k,508) = rxt(k,359)*y(k,141)
         mat(k,1367) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,1393) = .530_r8*rxt(k,396)*y(k,131) + .530_r8*rxt(k,395)*y(k,133) &
                      + .530_r8*rxt(k,392)*y(k,238) + .260_r8*rxt(k,393)*y(k,239)
         mat(k,1878) = mat(k,1878) + .700_r8*rxt(k,436)*y(k,101) + .500_r8*rxt(k,437) &
                      *y(k,102) + rxt(k,448)*y(k,122)
         mat(k,1219) = .225_r8*rxt(k,512)*y(k,131) + .125_r8*rxt(k,510)*y(k,239)
         mat(k,1273) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,238) + .100_r8*rxt(k,451)*y(k,239)
         mat(k,1084) = -(rxt(k,361)*y(k,263))
         mat(k,1867) = -rxt(k,361)*y(k,77)
         mat(k,319) = .650_r8*rxt(k,336)*y(k,263)
         mat(k,1206) = .200_r8*rxt(k,360)*y(k,263)
         mat(k,1158) = rxt(k,449)*y(k,263)
         mat(k,2201) = rxt(k,460)*y(k,233) + .050_r8*rxt(k,434)*y(k,248) &
                      + .400_r8*rxt(k,474)*y(k,250) + .170_r8*rxt(k,477)*y(k,252) &
                      + .700_r8*rxt(k,480)*y(k,265) + .600_r8*rxt(k,487)*y(k,270) &
                      + .250_r8*rxt(k,453)*y(k,271) + .340_r8*rxt(k,493)*y(k,272) &
                      + .170_r8*rxt(k,496)*y(k,273)
         mat(k,2259) = .050_r8*rxt(k,435)*y(k,248) + .250_r8*rxt(k,454)*y(k,271)
         mat(k,541) = rxt(k,460)*y(k,131)
         mat(k,1456) = .250_r8*rxt(k,450)*y(k,271)
         mat(k,1657) = .100_r8*rxt(k,451)*y(k,271)
         mat(k,2381) = .160_r8*rxt(k,473)*y(k,250) + .070_r8*rxt(k,476)*y(k,252)
         mat(k,1365) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,783) = .400_r8*rxt(k,474)*y(k,131) + .160_r8*rxt(k,473)*y(k,245)
         mat(k,968) = .170_r8*rxt(k,477)*y(k,131) + .070_r8*rxt(k,476)*y(k,245)
         mat(k,1867) = mat(k,1867) + .650_r8*rxt(k,336)*y(k,26) + .200_r8*rxt(k,360) &
                      *y(k,76) + rxt(k,449)*y(k,123)
         mat(k,499) = .700_r8*rxt(k,480)*y(k,131)
         mat(k,796) = .600_r8*rxt(k,487)*y(k,131)
         mat(k,1271) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,238) + .100_r8*rxt(k,451)*y(k,239)
         mat(k,812) = .340_r8*rxt(k,493)*y(k,131)
         mat(k,555) = .170_r8*rxt(k,496)*y(k,131)
         mat(k,1730) = -((rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,245) + rxt(k,195) &
                      *y(k,141) + rxt(k,198)*y(k,142))
         mat(k,2412) = -(rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,78)
         mat(k,1586) = -rxt(k,195)*y(k,78)
         mat(k,2099) = -rxt(k,198)*y(k,78)
         mat(k,2438) = rxt(k,317)*y(k,263)
         mat(k,1634) = rxt(k,331)*y(k,259)
         mat(k,1992) = rxt(k,236)*y(k,79)
         mat(k,946) = rxt(k,292)*y(k,79)
         mat(k,1508) = rxt(k,236)*y(k,58) + rxt(k,292)*y(k,75) + rxt(k,190)*y(k,140) &
                      + rxt(k,173)*y(k,259) + rxt(k,199)*y(k,263)
         mat(k,873) = rxt(k,290)*y(k,259)
         mat(k,1541) = rxt(k,267)*y(k,259)
         mat(k,1005) = rxt(k,222)*y(k,263)
         mat(k,1946) = rxt(k,190)*y(k,79) + rxt(k,202)*y(k,263)
         mat(k,420) = rxt(k,525)*y(k,263)
         mat(k,778) = rxt(k,531)*y(k,263)
         mat(k,1522) = rxt(k,536)*y(k,263)
         mat(k,2035) = rxt(k,331)*y(k,56) + rxt(k,173)*y(k,79) + rxt(k,290)*y(k,83) &
                      + rxt(k,267)*y(k,87)
         mat(k,1903) = rxt(k,317)*y(k,44) + rxt(k,199)*y(k,79) + rxt(k,222)*y(k,119) &
                      + rxt(k,202)*y(k,140) + rxt(k,525)*y(k,145) + rxt(k,531) &
                      *y(k,158) + rxt(k,536)*y(k,160)
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
         mat(k,1503) = -(rxt(k,173)*y(k,259) + rxt(k,190)*y(k,140) + rxt(k,199) &
                      *y(k,263) + rxt(k,236)*y(k,58) + rxt(k,292)*y(k,75))
         mat(k,2027) = -rxt(k,173)*y(k,79)
         mat(k,1937) = -rxt(k,190)*y(k,79)
         mat(k,1894) = -rxt(k,199)*y(k,79)
         mat(k,1984) = -rxt(k,236)*y(k,79)
         mat(k,942) = -rxt(k,292)*y(k,79)
         mat(k,1626) = rxt(k,332)*y(k,259)
         mat(k,1722) = rxt(k,192)*y(k,245)
         mat(k,2404) = rxt(k,192)*y(k,78)
         mat(k,2027) = mat(k,2027) + rxt(k,332)*y(k,56)
         mat(k,151) = -(rxt(k,288)*y(k,259))
         mat(k,2006) = -rxt(k,288)*y(k,80)
         mat(k,671) = -(rxt(k,191)*y(k,140) + rxt(k,200)*y(k,263) + rxt(k,237)*y(k,58))
         mat(k,1922) = -rxt(k,191)*y(k,81)
         mat(k,1834) = -rxt(k,200)*y(k,81)
         mat(k,1968) = -rxt(k,237)*y(k,81)
         mat(k,2357) = 2.000_r8*rxt(k,206)*y(k,245)
         mat(k,1834) = mat(k,1834) + 2.000_r8*rxt(k,205)*y(k,263)
         mat(k,304) = rxt(k,538)*y(k,274)
         mat(k,2496) = rxt(k,538)*y(k,162)
         mat(k,870) = -(rxt(k,283)*y(k,140) + rxt(k,284)*y(k,263) + (rxt(k,289) &
                      + rxt(k,290)) * y(k,259))
         mat(k,1927) = -rxt(k,283)*y(k,83)
         mat(k,1853) = -rxt(k,284)*y(k,83)
         mat(k,2025) = -(rxt(k,289) + rxt(k,290)) * y(k,83)
         mat(k,1550) = rxt(k,270)*y(k,44) + rxt(k,271)*y(k,245)
         mat(k,2426) = rxt(k,270)*y(k,19)
         mat(k,2373) = rxt(k,271)*y(k,19)
         mat(k,285) = -(rxt(k,306)*y(k,263) + rxt(k,311)*y(k,259))
         mat(k,1784) = -rxt(k,306)*y(k,84)
         mat(k,2017) = -rxt(k,311)*y(k,84)
         mat(k,295) = -(rxt(k,307)*y(k,263) + rxt(k,312)*y(k,259))
         mat(k,1786) = -rxt(k,307)*y(k,85)
         mat(k,2019) = -rxt(k,312)*y(k,85)
         mat(k,354) = -(rxt(k,308)*y(k,263) + rxt(k,313)*y(k,259))
         mat(k,1794) = -rxt(k,308)*y(k,86)
         mat(k,2021) = -rxt(k,313)*y(k,86)
         mat(k,1537) = -(rxt(k,254)*y(k,140) + rxt(k,255)*y(k,263) + (rxt(k,266) &
                      + rxt(k,267)) * y(k,259) + (rxt(k,587) + rxt(k,593) + rxt(k,598) &
                      ) * y(k,94) + (rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,62) &
                      + (rxt(k,594) + rxt(k,599)) * y(k,93))
         mat(k,1939) = -rxt(k,254)*y(k,87)
         mat(k,1896) = -rxt(k,255)*y(k,87)
         mat(k,2028) = -(rxt(k,266) + rxt(k,267)) * y(k,87)
         mat(k,880) = -(rxt(k,587) + rxt(k,593) + rxt(k,598)) * y(k,87)
         mat(k,1119) = -(rxt(k,592) + rxt(k,597) + rxt(k,602)) * y(k,87)
         mat(k,830) = -(rxt(k,594) + rxt(k,599)) * y(k,87)
         mat(k,325) = rxt(k,345)*y(k,58)
         mat(k,331) = rxt(k,375)*y(k,58)
         mat(k,520) = rxt(k,297)*y(k,58)
         mat(k,2431) = rxt(k,234)*y(k,58)
         mat(k,627) = rxt(k,299)*y(k,58)
         mat(k,401) = 2.000_r8*rxt(k,302)*y(k,58)
         mat(k,1627) = rxt(k,235)*y(k,58)
         mat(k,448) = rxt(k,304)*y(k,58)
         mat(k,1985) = rxt(k,345)*y(k,30) + rxt(k,375)*y(k,33) + rxt(k,297)*y(k,43) &
                      + rxt(k,234)*y(k,44) + rxt(k,299)*y(k,45) + 2.000_r8*rxt(k,302) &
                      *y(k,48) + rxt(k,235)*y(k,56) + rxt(k,304)*y(k,57) + rxt(k,236) &
                      *y(k,79) + rxt(k,237)*y(k,81) + rxt(k,256)*y(k,94) + rxt(k,238) &
                      *y(k,245)
         mat(k,2120) = rxt(k,253)*y(k,263)
         mat(k,1504) = rxt(k,236)*y(k,58)
         mat(k,672) = rxt(k,237)*y(k,58)
         mat(k,880) = mat(k,880) + rxt(k,256)*y(k,58)
         mat(k,2405) = rxt(k,238)*y(k,58)
         mat(k,1896) = mat(k,1896) + rxt(k,253)*y(k,61)
         mat(k,226) = -(rxt(k,325)*y(k,263) + rxt(k,333)*y(k,259))
         mat(k,1772) = -rxt(k,325)*y(k,88)
         mat(k,2015) = -rxt(k,333)*y(k,88)
         mat(k,962) = -(rxt(k,326)*y(k,263))
         mat(k,1859) = -rxt(k,326)*y(k,89)
         mat(k,1047) = .050_r8*rxt(k,500)*y(k,142)
         mat(k,318) = .350_r8*rxt(k,336)*y(k,263)
         mat(k,580) = .370_r8*rxt(k,338)*y(k,142)
         mat(k,1174) = .120_r8*rxt(k,367)*y(k,142)
         mat(k,925) = .110_r8*rxt(k,445)*y(k,142)
         mat(k,1309) = .330_r8*rxt(k,398)*y(k,142)
         mat(k,1019) = .050_r8*rxt(k,503)*y(k,142)
         mat(k,1411) = .120_r8*rxt(k,412)*y(k,142)
         mat(k,2196) = rxt(k,329)*y(k,246)
         mat(k,2064) = .050_r8*rxt(k,500)*y(k,6) + .370_r8*rxt(k,338)*y(k,27) &
                      + .120_r8*rxt(k,367)*y(k,31) + .110_r8*rxt(k,445)*y(k,100) &
                      + .330_r8*rxt(k,398)*y(k,111) + .050_r8*rxt(k,503)*y(k,116) &
                      + .120_r8*rxt(k,412)*y(k,118)
         mat(k,2377) = rxt(k,327)*y(k,246)
         mat(k,492) = rxt(k,329)*y(k,131) + rxt(k,327)*y(k,245)
         mat(k,1859) = mat(k,1859) + .350_r8*rxt(k,336)*y(k,26)
         mat(k,1622) = rxt(k,291)*y(k,75)
         mat(k,940) = rxt(k,291)*y(k,56) + rxt(k,292)*y(k,79) + rxt(k,294)*y(k,91) &
                      + rxt(k,293)*y(k,274)
         mat(k,1501) = rxt(k,292)*y(k,75)
         mat(k,1699) = rxt(k,294)*y(k,75)
         mat(k,2498) = rxt(k,293)*y(k,75)
         mat(k,1708) = -(rxt(k,231)*y(k,263) + rxt(k,294)*y(k,75))
         mat(k,1902) = -rxt(k,231)*y(k,91)
         mat(k,945) = -rxt(k,294)*y(k,91)
         mat(k,2437) = rxt(k,315)*y(k,133)
         mat(k,1197) = rxt(k,347)*y(k,133)
         mat(k,1326) = rxt(k,373)*y(k,133)
         mat(k,1120) = (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,87)
         mat(k,349) = rxt(k,522)*y(k,133)
         mat(k,1540) = (rxt(k,592)+rxt(k,597)+rxt(k,602))*y(k,62)
         mat(k,2482) = rxt(k,230)*y(k,263)
         mat(k,2291) = rxt(k,315)*y(k,44) + rxt(k,347)*y(k,47) + rxt(k,373)*y(k,51) &
                      + rxt(k,522)*y(k,69)
         mat(k,1902) = mat(k,1902) + rxt(k,230)*y(k,132)
         mat(k,512) = -(rxt(k,207)*y(k,263))
         mat(k,1816) = -rxt(k,207)*y(k,92)
         mat(k,2454) = rxt(k,228)*y(k,245)
         mat(k,2348) = rxt(k,228)*y(k,132)
         mat(k,829) = -(rxt(k,285)*y(k,140) + (rxt(k,594) + rxt(k,599)) * y(k,87))
         mat(k,1925) = -rxt(k,285)*y(k,93)
         mat(k,1534) = -(rxt(k,594) + rxt(k,599)) * y(k,93)
         mat(k,1601) = rxt(k,277)*y(k,245)
         mat(k,2370) = rxt(k,277)*y(k,21)
         mat(k,879) = -(rxt(k,256)*y(k,58) + rxt(k,257)*y(k,140) + rxt(k,258)*y(k,263) &
                      + (rxt(k,587) + rxt(k,593) + rxt(k,598)) * y(k,87))
         mat(k,1970) = -rxt(k,256)*y(k,94)
         mat(k,1928) = -rxt(k,257)*y(k,94)
         mat(k,1854) = -rxt(k,258)*y(k,94)
         mat(k,1535) = -(rxt(k,587) + rxt(k,593) + rxt(k,598)) * y(k,94)
         mat(k,2116) = rxt(k,245)*y(k,245)
         mat(k,1117) = rxt(k,250)*y(k,263)
         mat(k,2374) = rxt(k,245)*y(k,61)
         mat(k,1854) = mat(k,1854) + rxt(k,250)*y(k,62)
         mat(k,1237) = -(rxt(k,391)*y(k,263))
         mat(k,1880) = -rxt(k,391)*y(k,95)
         mat(k,644) = .300_r8*rxt(k,436)*y(k,263)
         mat(k,613) = .500_r8*rxt(k,437)*y(k,263)
         mat(k,2212) = rxt(k,390)*y(k,242) + rxt(k,397)*y(k,249)
         mat(k,621) = rxt(k,390)*y(k,131)
         mat(k,1394) = rxt(k,397)*y(k,131)
         mat(k,1880) = mat(k,1880) + .300_r8*rxt(k,436)*y(k,101) + .500_r8*rxt(k,437) &
                      *y(k,102)
         mat(k,277) = -(rxt(k,422)*y(k,263))
         mat(k,1782) = -rxt(k,422)*y(k,96)
         mat(k,1250) = -(rxt(k,376)*y(k,263))
         mat(k,1881) = -rxt(k,376)*y(k,97)
         mat(k,645) = .700_r8*rxt(k,436)*y(k,263)
         mat(k,614) = .500_r8*rxt(k,437)*y(k,263)
         mat(k,663) = .500_r8*rxt(k,411)*y(k,263)
         mat(k,2213) = .050_r8*rxt(k,434)*y(k,248) + .220_r8*rxt(k,396)*y(k,249) &
                      + .250_r8*rxt(k,453)*y(k,271)
         mat(k,2272) = .050_r8*rxt(k,435)*y(k,248) + .220_r8*rxt(k,395)*y(k,249) &
                      + .250_r8*rxt(k,454)*y(k,271)
         mat(k,597) = .500_r8*rxt(k,380)*y(k,263)
         mat(k,1460) = .220_r8*rxt(k,392)*y(k,249) + .250_r8*rxt(k,450)*y(k,271)
         mat(k,1668) = .230_r8*rxt(k,393)*y(k,249) + .200_r8*rxt(k,381)*y(k,267) &
                      + .100_r8*rxt(k,451)*y(k,271)
         mat(k,1369) = .050_r8*rxt(k,434)*y(k,131) + .050_r8*rxt(k,435)*y(k,133)
         mat(k,1395) = .220_r8*rxt(k,396)*y(k,131) + .220_r8*rxt(k,395)*y(k,133) &
                      + .220_r8*rxt(k,392)*y(k,238) + .230_r8*rxt(k,393)*y(k,239)
         mat(k,1881) = mat(k,1881) + .700_r8*rxt(k,436)*y(k,101) + .500_r8*rxt(k,437) &
                      *y(k,102) + .500_r8*rxt(k,411)*y(k,115) + .500_r8*rxt(k,380) &
                      *y(k,156)
         mat(k,1258) = .200_r8*rxt(k,381)*y(k,239)
         mat(k,1274) = .250_r8*rxt(k,453)*y(k,131) + .250_r8*rxt(k,454)*y(k,133) &
                      + .250_r8*rxt(k,450)*y(k,238) + .100_r8*rxt(k,451)*y(k,239)
         mat(k,381) = -(rxt(k,423)*y(k,263))
         mat(k,1798) = -rxt(k,423)*y(k,98)
         mat(k,2165) = .870_r8*rxt(k,434)*y(k,248)
         mat(k,2249) = .950_r8*rxt(k,435)*y(k,248)
         mat(k,1452) = rxt(k,430)*y(k,248)
         mat(k,1648) = .750_r8*rxt(k,431)*y(k,248)
         mat(k,1358) = .870_r8*rxt(k,434)*y(k,131) + .950_r8*rxt(k,435)*y(k,133) &
                      + rxt(k,430)*y(k,238) + .750_r8*rxt(k,431)*y(k,239)
         mat(k,171) = -(rxt(k,424)*y(k,263))
         mat(k,1768) = -rxt(k,424)*y(k,99)
         mat(k,752) = .600_r8*rxt(k,447)*y(k,263)
         mat(k,1768) = mat(k,1768) + .600_r8*rxt(k,447)*y(k,106)
         mat(k,924) = -(rxt(k,438)*y(k,133) + rxt(k,445)*y(k,142) + rxt(k,446) &
                      *y(k,263))
         mat(k,2253) = -rxt(k,438)*y(k,100)
         mat(k,2063) = -rxt(k,445)*y(k,100)
         mat(k,1856) = -rxt(k,446)*y(k,100)
         mat(k,642) = -(rxt(k,436)*y(k,263))
         mat(k,1831) = -rxt(k,436)*y(k,101)
         mat(k,2179) = .080_r8*rxt(k,428)*y(k,247)
         mat(k,1331) = .080_r8*rxt(k,428)*y(k,131)
         mat(k,610) = -(rxt(k,437)*y(k,263))
         mat(k,1827) = -rxt(k,437)*y(k,102)
         mat(k,2177) = .080_r8*rxt(k,434)*y(k,248)
         mat(k,1359) = .080_r8*rxt(k,434)*y(k,131)
         mat(k,1295) = -(rxt(k,439)*y(k,238) + rxt(k,440)*y(k,239) + rxt(k,441) &
                      *y(k,245) + rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133))
         mat(k,1462) = -rxt(k,439)*y(k,103)
         mat(k,1671) = -rxt(k,440)*y(k,103)
         mat(k,2394) = -rxt(k,441)*y(k,103)
         mat(k,2216) = -rxt(k,442)*y(k,103)
         mat(k,2275) = -rxt(k,443)*y(k,103)
         mat(k,928) = rxt(k,438)*y(k,133)
         mat(k,2275) = mat(k,2275) + rxt(k,438)*y(k,100)
         mat(k,471) = -(rxt(k,444)*y(k,263))
         mat(k,1811) = -rxt(k,444)*y(k,104)
         mat(k,1287) = rxt(k,441)*y(k,245)
         mat(k,2341) = rxt(k,441)*y(k,103)
         mat(k,88) = -(rxt(k,562)*y(k,245) + rxt(k,563)*y(k,131))
         mat(k,2317) = -rxt(k,562)*y(k,105)
         mat(k,2151) = -rxt(k,563)*y(k,105)
         mat(k,923) = rxt(k,565)*y(k,263)
         mat(k,1751) = rxt(k,565)*y(k,100)
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
         mat(k,753) = -(rxt(k,447)*y(k,263))
         mat(k,1843) = -rxt(k,447)*y(k,106)
         mat(k,2364) = rxt(k,427)*y(k,247) + rxt(k,432)*y(k,248)
         mat(k,1332) = rxt(k,427)*y(k,245)
         mat(k,1361) = rxt(k,432)*y(k,245)
         mat(k,71) = -(rxt(k,568)*y(k,263))
         mat(k,1749) = -rxt(k,568)*y(k,107)
         mat(k,69) = -(rxt(k,566)*y(k,245) + rxt(k,567)*y(k,131))
         mat(k,2310) = -rxt(k,566)*y(k,108)
         mat(k,2144) = -rxt(k,567)*y(k,108)
         mat(k,70) = rxt(k,568)*y(k,263)
         mat(k,1748) = rxt(k,568)*y(k,107)
         mat(k,107) = -(rxt(k,571)*y(k,263))
         mat(k,1759) = -rxt(k,571)*y(k,109)
         mat(k,105) = -(rxt(k,569)*y(k,245) + rxt(k,570)*y(k,131))
         mat(k,2324) = -rxt(k,569)*y(k,110)
         mat(k,2158) = -rxt(k,570)*y(k,110)
         mat(k,106) = rxt(k,571)*y(k,263)
         mat(k,1758) = rxt(k,571)*y(k,109)
         mat(k,1311) = -(rxt(k,398)*y(k,142) + rxt(k,399)*y(k,263))
         mat(k,2083) = -rxt(k,398)*y(k,111)
         mat(k,1885) = -rxt(k,399)*y(k,111)
         mat(k,929) = .300_r8*rxt(k,445)*y(k,142)
         mat(k,2217) = .360_r8*rxt(k,428)*y(k,247)
         mat(k,2276) = .400_r8*rxt(k,429)*y(k,247)
         mat(k,2083) = mat(k,2083) + .300_r8*rxt(k,445)*y(k,100)
         mat(k,1463) = .390_r8*rxt(k,425)*y(k,247)
         mat(k,1672) = .310_r8*rxt(k,426)*y(k,247)
         mat(k,1339) = .360_r8*rxt(k,428)*y(k,131) + .400_r8*rxt(k,429)*y(k,133) &
                      + .390_r8*rxt(k,425)*y(k,238) + .310_r8*rxt(k,426)*y(k,239)
         mat(k,384) = -(rxt(k,400)*y(k,263))
         mat(k,1799) = -rxt(k,400)*y(k,112)
         mat(k,2334) = rxt(k,394)*y(k,249)
         mat(k,1390) = rxt(k,394)*y(k,245)
         mat(k,561) = -(rxt(k,409)*y(k,263))
         mat(k,1822) = -rxt(k,409)*y(k,113)
         mat(k,2175) = .800_r8*rxt(k,418)*y(k,232)
         mat(k,1067) = .800_r8*rxt(k,418)*y(k,131)
         mat(k,394) = -(rxt(k,410)*y(k,263))
         mat(k,1800) = -rxt(k,410)*y(k,114)
         mat(k,2335) = .800_r8*rxt(k,407)*y(k,253)
         mat(k,718) = .800_r8*rxt(k,407)*y(k,245)
         mat(k,662) = -(rxt(k,411)*y(k,263))
         mat(k,1833) = -rxt(k,411)*y(k,115)
         mat(k,2459) = rxt(k,414)*y(k,251)
         mat(k,1435) = rxt(k,414)*y(k,132)
         mat(k,1020) = -(rxt(k,502)*y(k,133) + rxt(k,503)*y(k,142) + rxt(k,504) &
                      *y(k,263))
         mat(k,2256) = -rxt(k,502)*y(k,116)
         mat(k,2066) = -rxt(k,503)*y(k,116)
         mat(k,1864) = -rxt(k,504)*y(k,116)
         mat(k,94) = -(rxt(k,573)*y(k,245) + rxt(k,574)*y(k,131))
         mat(k,2318) = -rxt(k,573)*y(k,117)
         mat(k,2152) = -rxt(k,574)*y(k,117)
         mat(k,1016) = rxt(k,576)*y(k,263)
         mat(k,1752) = rxt(k,576)*y(k,116)
         mat(k,1418) = -(rxt(k,412)*y(k,142) + rxt(k,413)*y(k,263))
         mat(k,2088) = -rxt(k,412)*y(k,118)
         mat(k,1890) = -rxt(k,413)*y(k,118)
         mat(k,932) = .200_r8*rxt(k,445)*y(k,142)
         mat(k,2222) = .560_r8*rxt(k,428)*y(k,247)
         mat(k,2281) = .600_r8*rxt(k,429)*y(k,247)
         mat(k,2088) = mat(k,2088) + .200_r8*rxt(k,445)*y(k,100)
         mat(k,1468) = .610_r8*rxt(k,425)*y(k,247)
         mat(k,1677) = .440_r8*rxt(k,426)*y(k,247)
         mat(k,1343) = .560_r8*rxt(k,428)*y(k,131) + .600_r8*rxt(k,429)*y(k,133) &
                      + .610_r8*rxt(k,425)*y(k,238) + .440_r8*rxt(k,426)*y(k,239)
         mat(k,1001) = -(rxt(k,210)*y(k,131) + (rxt(k,211) + rxt(k,212) + rxt(k,213) &
                      ) * y(k,132) + rxt(k,214)*y(k,141) + rxt(k,222)*y(k,263) &
                      + rxt(k,612)*y(k,262))
         mat(k,2199) = -rxt(k,210)*y(k,119)
         mat(k,2466) = -(rxt(k,211) + rxt(k,212) + rxt(k,213)) * y(k,119)
         mat(k,1578) = -rxt(k,214)*y(k,119)
         mat(k,1863) = -rxt(k,222)*y(k,119)
         mat(k,897) = -rxt(k,612)*y(k,119)
         mat(k,1933) = rxt(k,208)*y(k,254) + rxt(k,609)*y(k,257)
         mat(k,1578) = mat(k,1578) + rxt(k,610)*y(k,257)
         mat(k,908) = 1.100_r8*rxt(k,605)*y(k,255) + .200_r8*rxt(k,603)*y(k,256)
         mat(k,574) = rxt(k,208)*y(k,140)
         mat(k,742) = 1.100_r8*rxt(k,605)*y(k,241)
         mat(k,889) = .200_r8*rxt(k,603)*y(k,241)
         mat(k,550) = rxt(k,609)*y(k,140) + rxt(k,610)*y(k,141)
         mat(k,300) = -((rxt(k,226) + rxt(k,227)) * y(k,259))
         mat(k,2020) = -(rxt(k,226) + rxt(k,227)) * y(k,120)
         mat(k,995) = rxt(k,211)*y(k,132)
         mat(k,2452) = rxt(k,211)*y(k,119)
         mat(k,2453) = rxt(k,229)*y(k,133)
         mat(k,2247) = rxt(k,229)*y(k,132)
         mat(k,459) = -(rxt(k,448)*y(k,263))
         mat(k,1809) = -rxt(k,448)*y(k,122)
         mat(k,1286) = .200_r8*rxt(k,440)*y(k,239)
         mat(k,1650) = .200_r8*rxt(k,440)*y(k,103)
         mat(k,1159) = -(rxt(k,449)*y(k,263))
         mat(k,1874) = -rxt(k,449)*y(k,123)
         mat(k,1291) = rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133) + rxt(k,439)*y(k,238) &
                      + .800_r8*rxt(k,440)*y(k,239)
         mat(k,2207) = rxt(k,442)*y(k,103)
         mat(k,2265) = rxt(k,443)*y(k,103)
         mat(k,1457) = rxt(k,439)*y(k,103)
         mat(k,1663) = .800_r8*rxt(k,440)*y(k,103)
         mat(k,139) = -(rxt(k,539)*y(k,263))
         mat(k,1765) = -rxt(k,539)*y(k,127)
         mat(k,2239) = -(rxt(k,210)*y(k,119) + rxt(k,219)*y(k,133) + rxt(k,223) &
                      *y(k,245) + rxt(k,224)*y(k,142) + rxt(k,225)*y(k,140) + rxt(k,246) &
                      *y(k,61) + rxt(k,278)*y(k,21) + rxt(k,321)*y(k,239) + rxt(k,329) &
                      *y(k,246) + rxt(k,342)*y(k,235) + rxt(k,353)*y(k,238) + rxt(k,357) &
                      *y(k,244) + rxt(k,370)*y(k,236) + rxt(k,379)*y(k,266) + rxt(k,383) &
                      *y(k,267) + (rxt(k,389) + rxt(k,390)) * y(k,242) + (rxt(k,396) &
                      + rxt(k,397)) * y(k,249) + rxt(k,405)*y(k,251) + rxt(k,408) &
                      *y(k,253) + (rxt(k,418) + rxt(k,419)) * y(k,232) + rxt(k,428) &
                      *y(k,247) + rxt(k,434)*y(k,248) + rxt(k,442)*y(k,103) + rxt(k,453) &
                      *y(k,271) + rxt(k,457)*y(k,231) + rxt(k,460)*y(k,233) + rxt(k,465) &
                      *y(k,234) + rxt(k,467)*y(k,237) + rxt(k,471)*y(k,240) + rxt(k,474) &
                      *y(k,250) + rxt(k,477)*y(k,252) + rxt(k,480)*y(k,265) + rxt(k,487) &
                      *y(k,270) + rxt(k,493)*y(k,272) + rxt(k,496)*y(k,273) + rxt(k,507) &
                      *y(k,258) + rxt(k,512)*y(k,268) + rxt(k,517)*y(k,269) + rxt(k,614) &
                      *y(k,262))
         mat(k,1009) = -rxt(k,210)*y(k,131)
         mat(k,2299) = -rxt(k,219)*y(k,131)
         mat(k,2419) = -rxt(k,223)*y(k,131)
         mat(k,2106) = -rxt(k,224)*y(k,131)
         mat(k,1953) = -rxt(k,225)*y(k,131)
         mat(k,2133) = -rxt(k,246)*y(k,131)
         mat(k,1616) = -rxt(k,278)*y(k,131)
         mat(k,1693) = -rxt(k,321)*y(k,131)
         mat(k,493) = -rxt(k,329)*y(k,131)
         mat(k,958) = -rxt(k,342)*y(k,131)
         mat(k,1478) = -rxt(k,353)*y(k,131)
         mat(k,843) = -rxt(k,357)*y(k,131)
         mat(k,983) = -rxt(k,370)*y(k,131)
         mat(k,865) = -rxt(k,379)*y(k,131)
         mat(k,1264) = -rxt(k,383)*y(k,131)
         mat(k,622) = -(rxt(k,389) + rxt(k,390)) * y(k,131)
         mat(k,1405) = -(rxt(k,396) + rxt(k,397)) * y(k,131)
         mat(k,1446) = -rxt(k,405)*y(k,131)
         mat(k,723) = -rxt(k,408)*y(k,131)
         mat(k,1079) = -(rxt(k,418) + rxt(k,419)) * y(k,131)
         mat(k,1351) = -rxt(k,428)*y(k,131)
         mat(k,1384) = -rxt(k,434)*y(k,131)
         mat(k,1304) = -rxt(k,442)*y(k,131)
         mat(k,1281) = -rxt(k,453)*y(k,131)
         mat(k,570) = -rxt(k,457)*y(k,131)
         mat(k,543) = -rxt(k,460)*y(k,131)
         mat(k,488) = -rxt(k,465)*y(k,131)
         mat(k,693) = -rxt(k,467)*y(k,131)
         mat(k,825) = -rxt(k,471)*y(k,131)
         mat(k,785) = -rxt(k,474)*y(k,131)
         mat(k,973) = -rxt(k,477)*y(k,131)
         mat(k,501) = -rxt(k,480)*y(k,131)
         mat(k,800) = -rxt(k,487)*y(k,131)
         mat(k,817) = -rxt(k,493)*y(k,131)
         mat(k,558) = -rxt(k,496)*y(k,131)
         mat(k,1139) = -rxt(k,507)*y(k,131)
         mat(k,1227) = -rxt(k,512)*y(k,131)
         mat(k,1098) = -rxt(k,517)*y(k,131)
         mat(k,901) = -rxt(k,614)*y(k,131)
         mat(k,1009) = mat(k,1009) + 2.000_r8*rxt(k,212)*y(k,132) + rxt(k,214) &
                      *y(k,141) + rxt(k,222)*y(k,263)
         mat(k,303) = 2.000_r8*rxt(k,226)*y(k,259)
         mat(k,2490) = 2.000_r8*rxt(k,212)*y(k,119) + rxt(k,215)*y(k,140) + rxt(k,532) &
                      *y(k,160)
         mat(k,1953) = mat(k,1953) + rxt(k,215)*y(k,132)
         mat(k,1593) = rxt(k,214)*y(k,119) + rxt(k,209)*y(k,254)
         mat(k,1528) = rxt(k,532)*y(k,132)
         mat(k,577) = rxt(k,209)*y(k,141)
         mat(k,2042) = 2.000_r8*rxt(k,226)*y(k,120)
         mat(k,1910) = rxt(k,222)*y(k,119)
         mat(k,2494) = -((rxt(k,211) + rxt(k,212) + rxt(k,213)) * y(k,119) + (rxt(k,215) &
                      + rxt(k,217)) * y(k,140) + rxt(k,216)*y(k,142) + rxt(k,228) &
                      *y(k,245) + rxt(k,229)*y(k,133) + rxt(k,230)*y(k,263) + rxt(k,248) &
                      *y(k,61) + rxt(k,279)*y(k,21) + rxt(k,364)*y(k,238) + rxt(k,414) &
                      *y(k,251) + rxt(k,472)*y(k,240) + rxt(k,475)*y(k,250) + rxt(k,478) &
                      *y(k,252) + rxt(k,482)*y(k,149) + rxt(k,485)*y(k,231) + rxt(k,532) &
                      *y(k,160))
         mat(k,1010) = -(rxt(k,211) + rxt(k,212) + rxt(k,213)) * y(k,132)
         mat(k,1957) = -(rxt(k,215) + rxt(k,217)) * y(k,132)
         mat(k,2110) = -rxt(k,216)*y(k,132)
         mat(k,2423) = -rxt(k,228)*y(k,132)
         mat(k,2303) = -rxt(k,229)*y(k,132)
         mat(k,1914) = -rxt(k,230)*y(k,132)
         mat(k,2137) = -rxt(k,248)*y(k,132)
         mat(k,1620) = -rxt(k,279)*y(k,132)
         mat(k,1482) = -rxt(k,364)*y(k,132)
         mat(k,1450) = -rxt(k,414)*y(k,132)
         mat(k,827) = -rxt(k,472)*y(k,132)
         mat(k,787) = -rxt(k,475)*y(k,132)
         mat(k,975) = -rxt(k,478)*y(k,132)
         mat(k,529) = -rxt(k,482)*y(k,132)
         mat(k,572) = -rxt(k,485)*y(k,132)
         mat(k,1530) = -rxt(k,532)*y(k,132)
         mat(k,717) = rxt(k,416)*y(k,263)
         mat(k,414) = rxt(k,387)*y(k,133)
         mat(k,1620) = mat(k,1620) + rxt(k,278)*y(k,131)
         mat(k,2137) = mat(k,2137) + rxt(k,246)*y(k,131)
         mat(k,517) = rxt(k,207)*y(k,263)
         mat(k,650) = .700_r8*rxt(k,436)*y(k,263)
         mat(k,1308) = rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133)
         mat(k,2243) = rxt(k,278)*y(k,21) + rxt(k,246)*y(k,61) + rxt(k,442)*y(k,103) &
                      + 2.000_r8*rxt(k,219)*y(k,133) + rxt(k,225)*y(k,140) &
                      + rxt(k,224)*y(k,142) + rxt(k,457)*y(k,231) + rxt(k,418) &
                      *y(k,232) + rxt(k,460)*y(k,233) + rxt(k,465)*y(k,234) &
                      + rxt(k,342)*y(k,235) + rxt(k,370)*y(k,236) + rxt(k,467) &
                      *y(k,237) + rxt(k,353)*y(k,238) + rxt(k,321)*y(k,239) &
                      + rxt(k,471)*y(k,240) + rxt(k,389)*y(k,242) + rxt(k,357) &
                      *y(k,244) + rxt(k,223)*y(k,245) + rxt(k,329)*y(k,246) &
                      + .920_r8*rxt(k,428)*y(k,247) + .920_r8*rxt(k,434)*y(k,248) &
                      + rxt(k,396)*y(k,249) + rxt(k,474)*y(k,250) + rxt(k,405) &
                      *y(k,251) + rxt(k,477)*y(k,252) + rxt(k,408)*y(k,253) &
                      + 1.600_r8*rxt(k,507)*y(k,258) + rxt(k,480)*y(k,265) &
                      + rxt(k,379)*y(k,266) + rxt(k,383)*y(k,267) + .900_r8*rxt(k,512) &
                      *y(k,268) + .800_r8*rxt(k,517)*y(k,269) + rxt(k,487)*y(k,270) &
                      + rxt(k,453)*y(k,271) + rxt(k,493)*y(k,272) + rxt(k,496) &
                      *y(k,273)
         mat(k,2303) = mat(k,2303) + rxt(k,387)*y(k,18) + rxt(k,443)*y(k,103) &
                      + 2.000_r8*rxt(k,219)*y(k,131) + rxt(k,220)*y(k,140) &
                      + rxt(k,218)*y(k,245) + rxt(k,429)*y(k,247) + rxt(k,435) &
                      *y(k,248) + rxt(k,395)*y(k,249) + rxt(k,406)*y(k,251) &
                      + 2.000_r8*rxt(k,508)*y(k,258) + rxt(k,221)*y(k,263) &
                      + rxt(k,454)*y(k,271)
         mat(k,920) = rxt(k,377)*y(k,263)
         mat(k,1957) = mat(k,1957) + rxt(k,225)*y(k,131) + rxt(k,220)*y(k,133)
         mat(k,2110) = mat(k,2110) + rxt(k,224)*y(k,131)
         mat(k,687) = rxt(k,514)*y(k,263)
         mat(k,572) = mat(k,572) + rxt(k,457)*y(k,131)
         mat(k,1082) = rxt(k,418)*y(k,131)
         mat(k,545) = rxt(k,460)*y(k,131)
         mat(k,490) = rxt(k,465)*y(k,131)
         mat(k,961) = rxt(k,342)*y(k,131)
         mat(k,986) = rxt(k,370)*y(k,131)
         mat(k,695) = rxt(k,467)*y(k,131)
         mat(k,1482) = mat(k,1482) + rxt(k,353)*y(k,131)
         mat(k,1697) = rxt(k,321)*y(k,131) + .500_r8*rxt(k,505)*y(k,258)
         mat(k,827) = mat(k,827) + rxt(k,471)*y(k,131)
         mat(k,625) = rxt(k,389)*y(k,131)
         mat(k,846) = rxt(k,357)*y(k,131)
         mat(k,2423) = mat(k,2423) + rxt(k,223)*y(k,131) + rxt(k,218)*y(k,133)
         mat(k,496) = rxt(k,329)*y(k,131)
         mat(k,1355) = .920_r8*rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133)
         mat(k,1388) = .920_r8*rxt(k,434)*y(k,131) + rxt(k,435)*y(k,133)
         mat(k,1409) = rxt(k,396)*y(k,131) + rxt(k,395)*y(k,133)
         mat(k,787) = mat(k,787) + rxt(k,474)*y(k,131)
         mat(k,1450) = mat(k,1450) + rxt(k,405)*y(k,131) + rxt(k,406)*y(k,133)
         mat(k,975) = mat(k,975) + rxt(k,477)*y(k,131)
         mat(k,725) = rxt(k,408)*y(k,131)
         mat(k,1143) = 1.600_r8*rxt(k,507)*y(k,131) + 2.000_r8*rxt(k,508)*y(k,133) &
                      + .500_r8*rxt(k,505)*y(k,239)
         mat(k,1914) = mat(k,1914) + rxt(k,416)*y(k,1) + rxt(k,207)*y(k,92) &
                      + .700_r8*rxt(k,436)*y(k,101) + rxt(k,221)*y(k,133) + rxt(k,377) &
                      *y(k,134) + rxt(k,514)*y(k,216)
         mat(k,503) = rxt(k,480)*y(k,131)
         mat(k,868) = rxt(k,379)*y(k,131)
         mat(k,1267) = rxt(k,383)*y(k,131)
         mat(k,1230) = .900_r8*rxt(k,512)*y(k,131)
         mat(k,1101) = .800_r8*rxt(k,517)*y(k,131)
         mat(k,802) = rxt(k,487)*y(k,131)
         mat(k,1285) = rxt(k,453)*y(k,131) + rxt(k,454)*y(k,133)
         mat(k,819) = rxt(k,493)*y(k,131)
         mat(k,560) = rxt(k,496)*y(k,131)
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
         mat(k,2300) = -(rxt(k,218)*y(k,245) + rxt(k,219)*y(k,131) + rxt(k,220) &
                      *y(k,140) + rxt(k,221)*y(k,263) + rxt(k,229)*y(k,132) + rxt(k,315) &
                      *y(k,44) + rxt(k,347)*y(k,47) + rxt(k,366)*y(k,31) + rxt(k,373) &
                      *y(k,51) + rxt(k,387)*y(k,18) + rxt(k,395)*y(k,249) + rxt(k,406) &
                      *y(k,251) + rxt(k,429)*y(k,247) + rxt(k,435)*y(k,248) + rxt(k,438) &
                      *y(k,100) + rxt(k,443)*y(k,103) + rxt(k,454)*y(k,271) + rxt(k,499) &
                      *y(k,6) + rxt(k,502)*y(k,116) + rxt(k,508)*y(k,258) + rxt(k,519) &
                      *y(k,218) + rxt(k,522)*y(k,69))
         mat(k,2420) = -rxt(k,218)*y(k,133)
         mat(k,2240) = -rxt(k,219)*y(k,133)
         mat(k,1954) = -rxt(k,220)*y(k,133)
         mat(k,1911) = -rxt(k,221)*y(k,133)
         mat(k,2491) = -rxt(k,229)*y(k,133)
         mat(k,2446) = -rxt(k,315)*y(k,133)
         mat(k,1199) = -rxt(k,347)*y(k,133)
         mat(k,1188) = -rxt(k,366)*y(k,133)
         mat(k,1328) = -rxt(k,373)*y(k,133)
         mat(k,412) = -rxt(k,387)*y(k,133)
         mat(k,1406) = -rxt(k,395)*y(k,133)
         mat(k,1447) = -rxt(k,406)*y(k,133)
         mat(k,1352) = -rxt(k,429)*y(k,133)
         mat(k,1385) = -rxt(k,435)*y(k,133)
         mat(k,937) = -rxt(k,438)*y(k,133)
         mat(k,1305) = -rxt(k,443)*y(k,133)
         mat(k,1282) = -rxt(k,454)*y(k,133)
         mat(k,1063) = -rxt(k,499)*y(k,133)
         mat(k,1035) = -rxt(k,502)*y(k,133)
         mat(k,1140) = -rxt(k,508)*y(k,133)
         mat(k,1108) = -rxt(k,519)*y(k,133)
         mat(k,351) = -rxt(k,522)*y(k,133)
         mat(k,592) = rxt(k,280)*y(k,140)
         mat(k,2000) = rxt(k,247)*y(k,62)
         mat(k,1125) = rxt(k,247)*y(k,58) + rxt(k,249)*y(k,140) + rxt(k,250)*y(k,263)
         mat(k,948) = rxt(k,294)*y(k,91)
         mat(k,1717) = rxt(k,294)*y(k,75) + rxt(k,231)*y(k,263)
         mat(k,667) = .500_r8*rxt(k,411)*y(k,263)
         mat(k,2491) = mat(k,2491) + rxt(k,217)*y(k,140) + rxt(k,216)*y(k,142)
         mat(k,1954) = mat(k,1954) + rxt(k,280)*y(k,22) + rxt(k,249)*y(k,62) &
                      + rxt(k,217)*y(k,132)
         mat(k,2107) = rxt(k,216)*y(k,132)
         mat(k,607) = rxt(k,362)*y(k,263)
         mat(k,1911) = mat(k,1911) + rxt(k,250)*y(k,62) + rxt(k,231)*y(k,91) &
                      + .500_r8*rxt(k,411)*y(k,115) + rxt(k,362)*y(k,147)
         mat(k,915) = -(rxt(k,377)*y(k,263))
         mat(k,1855) = -rxt(k,377)*y(k,134)
         mat(k,1173) = rxt(k,366)*y(k,133)
         mat(k,611) = .500_r8*rxt(k,437)*y(k,263)
         mat(k,473) = rxt(k,444)*y(k,263)
         mat(k,460) = rxt(k,448)*y(k,263)
         mat(k,1156) = rxt(k,449)*y(k,263)
         mat(k,2252) = rxt(k,366)*y(k,31)
         mat(k,1855) = mat(k,1855) + .500_r8*rxt(k,437)*y(k,102) + rxt(k,444)*y(k,104) &
                      + rxt(k,448)*y(k,122) + rxt(k,449)*y(k,123)
         mat(k,453) = -(rxt(k,509)*y(k,263))
         mat(k,1808) = -rxt(k,509)*y(k,135)
         mat(k,2339) = rxt(k,506)*y(k,258)
         mat(k,1128) = rxt(k,506)*y(k,245)
         mat(k,1948) = -(rxt(k,187)*y(k,142) + 4._r8*rxt(k,188)*y(k,140) + rxt(k,189) &
                      *y(k,141) + rxt(k,190)*y(k,79) + rxt(k,191)*y(k,81) + rxt(k,196) &
                      *y(k,245) + rxt(k,202)*y(k,263) + (rxt(k,215) + rxt(k,217) &
                      ) * y(k,132) + rxt(k,220)*y(k,133) + rxt(k,225)*y(k,131) &
                      + rxt(k,249)*y(k,62) + rxt(k,251)*y(k,61) + rxt(k,254)*y(k,87) &
                      + rxt(k,257)*y(k,94) + rxt(k,280)*y(k,22) + rxt(k,281)*y(k,21) &
                      + rxt(k,283)*y(k,83) + rxt(k,285)*y(k,93) + rxt(k,316)*y(k,44) &
                      + rxt(k,524)*y(k,145) + (rxt(k,607) + rxt(k,608)) * y(k,255) &
                      + rxt(k,609)*y(k,257))
         mat(k,2101) = -rxt(k,187)*y(k,140)
         mat(k,1588) = -rxt(k,189)*y(k,140)
         mat(k,1510) = -rxt(k,190)*y(k,140)
         mat(k,674) = -rxt(k,191)*y(k,140)
         mat(k,2414) = -rxt(k,196)*y(k,140)
         mat(k,1905) = -rxt(k,202)*y(k,140)
         mat(k,2485) = -(rxt(k,215) + rxt(k,217)) * y(k,140)
         mat(k,2294) = -rxt(k,220)*y(k,140)
         mat(k,2234) = -rxt(k,225)*y(k,140)
         mat(k,1122) = -rxt(k,249)*y(k,140)
         mat(k,2128) = -rxt(k,251)*y(k,140)
         mat(k,1543) = -rxt(k,254)*y(k,140)
         mat(k,882) = -rxt(k,257)*y(k,140)
         mat(k,591) = -rxt(k,280)*y(k,140)
         mat(k,1611) = -rxt(k,281)*y(k,140)
         mat(k,875) = -rxt(k,283)*y(k,140)
         mat(k,834) = -rxt(k,285)*y(k,140)
         mat(k,2440) = -rxt(k,316)*y(k,140)
         mat(k,422) = -rxt(k,524)*y(k,140)
         mat(k,744) = -(rxt(k,607) + rxt(k,608)) * y(k,140)
         mat(k,552) = -rxt(k,609)*y(k,140)
         mat(k,1732) = rxt(k,194)*y(k,245)
         mat(k,1007) = rxt(k,210)*y(k,131) + rxt(k,211)*y(k,132) + rxt(k,214)*y(k,141) &
                      + rxt(k,612)*y(k,262)
         mat(k,2234) = mat(k,2234) + rxt(k,210)*y(k,119)
         mat(k,2485) = mat(k,2485) + rxt(k,211)*y(k,119)
         mat(k,1588) = mat(k,1588) + rxt(k,214)*y(k,119) + rxt(k,526)*y(k,158) &
                      + rxt(k,533)*y(k,160) + rxt(k,611)*y(k,257) + (rxt(k,176) &
                       +rxt(k,177))*y(k,259) + rxt(k,617)*y(k,264)
         mat(k,780) = rxt(k,526)*y(k,141)
         mat(k,1524) = rxt(k,533)*y(k,141)
         mat(k,912) = rxt(k,603)*y(k,256) + 1.150_r8*rxt(k,604)*y(k,262)
         mat(k,2414) = mat(k,2414) + rxt(k,194)*y(k,78)
         mat(k,891) = rxt(k,603)*y(k,241)
         mat(k,552) = mat(k,552) + rxt(k,611)*y(k,141)
         mat(k,2037) = (rxt(k,176)+rxt(k,177))*y(k,141)
         mat(k,899) = rxt(k,612)*y(k,119) + 1.150_r8*rxt(k,604)*y(k,241)
         mat(k,1905) = mat(k,1905) + 2.000_r8*rxt(k,204)*y(k,263)
         mat(k,854) = rxt(k,617)*y(k,141)
         mat(k,1584) = -(rxt(k,176)*y(k,259) + rxt(k,181)*y(k,260) + rxt(k,189) &
                      *y(k,140) + rxt(k,195)*y(k,78) + rxt(k,209)*y(k,254) + rxt(k,214) &
                      *y(k,119) + rxt(k,359)*y(k,243) + rxt(k,526)*y(k,158) + rxt(k,533) &
                      *y(k,160) + rxt(k,606)*y(k,255) + (rxt(k,610) + rxt(k,611) &
                      ) * y(k,257) + rxt(k,617)*y(k,264))
         mat(k,2030) = -rxt(k,176)*y(k,141)
         mat(k,222) = -rxt(k,181)*y(k,141)
         mat(k,1941) = -rxt(k,189)*y(k,141)
         mat(k,1725) = -rxt(k,195)*y(k,141)
         mat(k,575) = -rxt(k,209)*y(k,141)
         mat(k,1004) = -rxt(k,214)*y(k,141)
         mat(k,509) = -rxt(k,359)*y(k,141)
         mat(k,777) = -rxt(k,526)*y(k,141)
         mat(k,1520) = -rxt(k,533)*y(k,141)
         mat(k,743) = -rxt(k,606)*y(k,141)
         mat(k,551) = -(rxt(k,610) + rxt(k,611)) * y(k,141)
         mat(k,853) = -rxt(k,617)*y(k,141)
         mat(k,1554) = rxt(k,272)*y(k,142) + rxt(k,271)*y(k,245)
         mat(k,1606) = 2.000_r8*rxt(k,273)*y(k,21) + (rxt(k,275)+rxt(k,276))*y(k,61) &
                      + rxt(k,281)*y(k,140) + rxt(k,277)*y(k,245)
         mat(k,1987) = rxt(k,240)*y(k,142) + rxt(k,238)*y(k,245)
         mat(k,2122) = (rxt(k,275)+rxt(k,276))*y(k,21) + (2.000_r8*rxt(k,242) &
                       +2.000_r8*rxt(k,243))*y(k,61) + rxt(k,251)*y(k,140) &
                      + rxt(k,245)*y(k,245) + rxt(k,253)*y(k,263)
         mat(k,1725) = mat(k,1725) + rxt(k,198)*y(k,142) + rxt(k,192)*y(k,245)
         mat(k,513) = rxt(k,207)*y(k,263)
         mat(k,1004) = mat(k,1004) + rxt(k,213)*y(k,132)
         mat(k,301) = rxt(k,227)*y(k,259)
         mat(k,2227) = rxt(k,224)*y(k,142) + rxt(k,614)*y(k,262)
         mat(k,2478) = rxt(k,213)*y(k,119) + rxt(k,215)*y(k,140) + rxt(k,216)*y(k,142)
         mat(k,2287) = rxt(k,220)*y(k,140) + rxt(k,218)*y(k,245)
         mat(k,1941) = mat(k,1941) + rxt(k,281)*y(k,21) + rxt(k,251)*y(k,61) &
                      + rxt(k,215)*y(k,132) + rxt(k,220)*y(k,133) &
                      + 2.000_r8*rxt(k,188)*y(k,140) + 2.000_r8*rxt(k,187)*y(k,142) &
                      + rxt(k,196)*y(k,245) + rxt(k,180)*y(k,260) + rxt(k,202) &
                      *y(k,263)
         mat(k,1584) = mat(k,1584) + 2.000_r8*rxt(k,181)*y(k,260)
         mat(k,2094) = rxt(k,272)*y(k,19) + rxt(k,240)*y(k,58) + rxt(k,198)*y(k,78) &
                      + rxt(k,224)*y(k,131) + rxt(k,216)*y(k,132) &
                      + 2.000_r8*rxt(k,187)*y(k,140) + rxt(k,528)*y(k,158) &
                      + rxt(k,534)*y(k,160) + 2.000_r8*rxt(k,197)*y(k,245) &
                      + 2.000_r8*rxt(k,178)*y(k,259) + rxt(k,203)*y(k,263)
         mat(k,777) = mat(k,777) + rxt(k,528)*y(k,142)
         mat(k,1520) = mat(k,1520) + rxt(k,534)*y(k,142)
         mat(k,955) = rxt(k,341)*y(k,245)
         mat(k,980) = rxt(k,369)*y(k,245)
         mat(k,1681) = rxt(k,320)*y(k,245)
         mat(k,2407) = rxt(k,271)*y(k,19) + rxt(k,277)*y(k,21) + rxt(k,238)*y(k,58) &
                      + rxt(k,245)*y(k,61) + rxt(k,192)*y(k,78) + rxt(k,218)*y(k,133) &
                      + rxt(k,196)*y(k,140) + 2.000_r8*rxt(k,197)*y(k,142) &
                      + rxt(k,341)*y(k,235) + rxt(k,369)*y(k,236) + rxt(k,320) &
                      *y(k,239) + 2.000_r8*rxt(k,206)*y(k,245) + rxt(k,201)*y(k,263) &
                      + rxt(k,378)*y(k,266)
         mat(k,2030) = mat(k,2030) + rxt(k,227)*y(k,120) + 2.000_r8*rxt(k,178) &
                      *y(k,142)
         mat(k,222) = mat(k,222) + rxt(k,180)*y(k,140) + 2.000_r8*rxt(k,181)*y(k,141)
         mat(k,898) = rxt(k,614)*y(k,131)
         mat(k,1898) = rxt(k,253)*y(k,61) + rxt(k,207)*y(k,92) + rxt(k,202)*y(k,140) &
                      + rxt(k,203)*y(k,142) + rxt(k,201)*y(k,245)
         mat(k,863) = rxt(k,378)*y(k,245)
         mat(k,2104) = -(rxt(k,178)*y(k,259) + rxt(k,187)*y(k,140) + rxt(k,197) &
                      *y(k,245) + rxt(k,198)*y(k,78) + rxt(k,203)*y(k,263) + rxt(k,216) &
                      *y(k,132) + rxt(k,224)*y(k,131) + rxt(k,240)*y(k,58) + rxt(k,272) &
                      *y(k,19) + rxt(k,338)*y(k,27) + rxt(k,367)*y(k,31) + rxt(k,398) &
                      *y(k,111) + rxt(k,412)*y(k,118) + rxt(k,445)*y(k,100) + rxt(k,483) &
                      *y(k,149) + rxt(k,500)*y(k,6) + rxt(k,503)*y(k,116) + rxt(k,528) &
                      *y(k,158) + rxt(k,534)*y(k,160))
         mat(k,2040) = -rxt(k,178)*y(k,142)
         mat(k,1951) = -rxt(k,187)*y(k,142)
         mat(k,2417) = -rxt(k,197)*y(k,142)
         mat(k,1735) = -rxt(k,198)*y(k,142)
         mat(k,1908) = -rxt(k,203)*y(k,142)
         mat(k,2488) = -rxt(k,216)*y(k,142)
         mat(k,2237) = -rxt(k,224)*y(k,142)
         mat(k,1997) = -rxt(k,240)*y(k,142)
         mat(k,1560) = -rxt(k,272)*y(k,142)
         mat(k,583) = -rxt(k,338)*y(k,142)
         mat(k,1186) = -rxt(k,367)*y(k,142)
         mat(k,1318) = -rxt(k,398)*y(k,142)
         mat(k,1427) = -rxt(k,412)*y(k,142)
         mat(k,936) = -rxt(k,445)*y(k,142)
         mat(k,528) = -rxt(k,483)*y(k,142)
         mat(k,1062) = -rxt(k,500)*y(k,142)
         mat(k,1034) = -rxt(k,503)*y(k,142)
         mat(k,781) = -rxt(k,528)*y(k,142)
         mat(k,1526) = -rxt(k,534)*y(k,142)
         mat(k,1951) = mat(k,1951) + rxt(k,189)*y(k,141)
         mat(k,1591) = rxt(k,189)*y(k,140)
         mat(k,1477) = .150_r8*rxt(k,352)*y(k,245)
         mat(k,2417) = mat(k,2417) + .150_r8*rxt(k,352)*y(k,238) + .150_r8*rxt(k,403) &
                      *y(k,251)
         mat(k,1445) = .150_r8*rxt(k,403)*y(k,245)
         mat(k,360) = -(rxt(k,535)*y(k,160))
         mat(k,1515) = -rxt(k,535)*y(k,144)
         mat(k,1599) = rxt(k,274)*y(k,61)
         mat(k,2115) = rxt(k,274)*y(k,21) + 2.000_r8*rxt(k,244)*y(k,61)
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
         mat(k,415) = -(rxt(k,524)*y(k,140) + rxt(k,525)*y(k,263))
         mat(k,1918) = -rxt(k,524)*y(k,145)
         mat(k,1803) = -rxt(k,525)*y(k,145)
         mat(k,1232) = rxt(k,391)*y(k,263)
         mat(k,2162) = .100_r8*rxt(k,512)*y(k,268)
         mat(k,1780) = rxt(k,391)*y(k,95)
         mat(k,1213) = .100_r8*rxt(k,512)*y(k,131)
         mat(k,602) = -(rxt(k,362)*y(k,263))
         mat(k,1826) = -rxt(k,362)*y(k,147)
         mat(k,2458) = rxt(k,364)*y(k,238)
         mat(k,1453) = rxt(k,364)*y(k,132)
         mat(k,2451) = rxt(k,485)*y(k,231)
         mat(k,566) = rxt(k,485)*y(k,132)
         mat(k,526) = -(rxt(k,482)*y(k,132) + rxt(k,483)*y(k,142))
         mat(k,2455) = -rxt(k,482)*y(k,149)
         mat(k,2057) = -rxt(k,483)*y(k,149)
         mat(k,240) = .070_r8*rxt(k,469)*y(k,263)
         mat(k,2172) = rxt(k,467)*y(k,237)
         mat(k,216) = .060_r8*rxt(k,481)*y(k,263)
         mat(k,261) = .070_r8*rxt(k,497)*y(k,263)
         mat(k,689) = rxt(k,467)*y(k,131)
         mat(k,1818) = .070_r8*rxt(k,469)*y(k,68) + .060_r8*rxt(k,481)*y(k,150) &
                      + .070_r8*rxt(k,497)*y(k,227)
         mat(k,214) = -(rxt(k,481)*y(k,263))
         mat(k,1771) = -rxt(k,481)*y(k,150)
         mat(k,206) = .530_r8*rxt(k,458)*y(k,263)
         mat(k,1771) = mat(k,1771) + .530_r8*rxt(k,458)*y(k,8)
         mat(k,365) = -(rxt(k,484)*y(k,263))
         mat(k,1795) = -rxt(k,484)*y(k,151)
         mat(k,2331) = rxt(k,479)*y(k,265)
         mat(k,497) = rxt(k,479)*y(k,245)
         mat(k,594) = -(rxt(k,380)*y(k,263))
         mat(k,1825) = -rxt(k,380)*y(k,156)
         mat(k,2354) = rxt(k,378)*y(k,266)
         mat(k,859) = rxt(k,378)*y(k,245)
         mat(k,429) = -(rxt(k,384)*y(k,263))
         mat(k,1804) = -rxt(k,384)*y(k,157)
         mat(k,2336) = .850_r8*rxt(k,382)*y(k,267)
         mat(k,1256) = .850_r8*rxt(k,382)*y(k,245)
         mat(k,775) = -(rxt(k,526)*y(k,141) + rxt(k,528)*y(k,142) + rxt(k,531) &
                      *y(k,263))
         mat(k,1572) = -rxt(k,526)*y(k,158)
         mat(k,2061) = -rxt(k,528)*y(k,158)
         mat(k,1845) = -rxt(k,531)*y(k,158)
         mat(k,1518) = -(rxt(k,529)*y(k,21) + rxt(k,530)*y(k,61) + rxt(k,532)*y(k,132) &
                      + rxt(k,533)*y(k,141) + rxt(k,534)*y(k,142) + rxt(k,535) &
                      *y(k,144) + rxt(k,536)*y(k,263))
         mat(k,1603) = -rxt(k,529)*y(k,160)
         mat(k,2119) = -rxt(k,530)*y(k,160)
         mat(k,2475) = -rxt(k,532)*y(k,160)
         mat(k,1582) = -rxt(k,533)*y(k,160)
         mat(k,2092) = -rxt(k,534)*y(k,160)
         mat(k,362) = -rxt(k,535)*y(k,160)
         mat(k,1895) = -rxt(k,536)*y(k,160)
         mat(k,1938) = rxt(k,524)*y(k,145)
         mat(k,1582) = mat(k,1582) + rxt(k,526)*y(k,158)
         mat(k,2092) = mat(k,2092) + rxt(k,528)*y(k,158)
         mat(k,419) = rxt(k,524)*y(k,140)
         mat(k,776) = rxt(k,526)*y(k,141) + rxt(k,528)*y(k,142) + rxt(k,531)*y(k,263)
         mat(k,1895) = mat(k,1895) + rxt(k,531)*y(k,158)
         mat(k,989) = -(rxt(k,527)*y(k,263))
         mat(k,1862) = -rxt(k,527)*y(k,161)
         mat(k,1602) = rxt(k,529)*y(k,160)
         mat(k,2117) = rxt(k,530)*y(k,160)
         mat(k,348) = rxt(k,522)*y(k,133) + (rxt(k,523)+.500_r8*rxt(k,537))*y(k,263)
         mat(k,2465) = rxt(k,532)*y(k,160)
         mat(k,2255) = rxt(k,522)*y(k,69)
         mat(k,1577) = rxt(k,533)*y(k,160)
         mat(k,2065) = rxt(k,534)*y(k,160)
         mat(k,361) = rxt(k,535)*y(k,160)
         mat(k,417) = rxt(k,525)*y(k,263)
         mat(k,1517) = rxt(k,529)*y(k,21) + rxt(k,530)*y(k,61) + rxt(k,532)*y(k,132) &
                      + rxt(k,533)*y(k,141) + rxt(k,534)*y(k,142) + rxt(k,535) &
                      *y(k,144) + rxt(k,536)*y(k,263)
         mat(k,1862) = mat(k,1862) + (rxt(k,523)+.500_r8*rxt(k,537))*y(k,69) &
                      + rxt(k,525)*y(k,145) + rxt(k,536)*y(k,160)
         mat(k,305) = -(rxt(k,538)*y(k,274))
         mat(k,2497) = -rxt(k,538)*y(k,162)
         mat(k,988) = rxt(k,527)*y(k,263)
         mat(k,1787) = rxt(k,527)*y(k,161)
         mat(k,64) = .1056005_r8*rxt(k,567)*y(k,131) + .2381005_r8*rxt(k,566)*y(k,245)
         mat(k,2139) = .1056005_r8*rxt(k,567)*y(k,108)
         mat(k,115) = .5931005_r8*rxt(k,577)*y(k,263)
         mat(k,2305) = .2381005_r8*rxt(k,566)*y(k,108)
         mat(k,1743) = .5931005_r8*rxt(k,577)*y(k,212)
         mat(k,65) = .1026005_r8*rxt(k,567)*y(k,131) + .1308005_r8*rxt(k,566)*y(k,245)
         mat(k,2140) = .1026005_r8*rxt(k,567)*y(k,108)
         mat(k,116) = .1534005_r8*rxt(k,577)*y(k,263)
         mat(k,2306) = .1308005_r8*rxt(k,566)*y(k,108)
         mat(k,1744) = .1534005_r8*rxt(k,577)*y(k,212)
         mat(k,66) = .0521005_r8*rxt(k,567)*y(k,131) + .0348005_r8*rxt(k,566)*y(k,245)
         mat(k,2141) = .0521005_r8*rxt(k,567)*y(k,108)
         mat(k,117) = .0459005_r8*rxt(k,577)*y(k,263)
         mat(k,2307) = .0348005_r8*rxt(k,566)*y(k,108)
         mat(k,1745) = .0459005_r8*rxt(k,577)*y(k,212)
         mat(k,67) = .0143005_r8*rxt(k,567)*y(k,131) + .0076005_r8*rxt(k,566)*y(k,245)
         mat(k,2142) = .0143005_r8*rxt(k,567)*y(k,108)
         mat(k,118) = .0085005_r8*rxt(k,577)*y(k,263)
         mat(k,2308) = .0076005_r8*rxt(k,566)*y(k,108)
         mat(k,1746) = .0085005_r8*rxt(k,577)*y(k,212)
         mat(k,68) = .0166005_r8*rxt(k,567)*y(k,131) + .0113005_r8*rxt(k,566)*y(k,245)
         mat(k,2143) = .0166005_r8*rxt(k,567)*y(k,108)
         mat(k,119) = .0128005_r8*rxt(k,577)*y(k,263)
         mat(k,2309) = .0113005_r8*rxt(k,566)*y(k,108)
         mat(k,1747) = .0128005_r8*rxt(k,577)*y(k,212)
         mat(k,1039) = .2202005_r8*rxt(k,556)*y(k,142)
         mat(k,77) = .1279005_r8*rxt(k,555)*y(k,131) + .2202005_r8*rxt(k,554)*y(k,245)
         mat(k,83) = .0003005_r8*rxt(k,563)*y(k,131) + .0031005_r8*rxt(k,562)*y(k,245)
         mat(k,1011) = .0508005_r8*rxt(k,575)*y(k,142)
         mat(k,89) = .0245005_r8*rxt(k,574)*y(k,131) + .0508005_r8*rxt(k,573)*y(k,245)
         mat(k,2145) = .1279005_r8*rxt(k,555)*y(k,7) + .0003005_r8*rxt(k,563)*y(k,105) &
                      + .0245005_r8*rxt(k,574)*y(k,117)
         mat(k,2048) = .2202005_r8*rxt(k,556)*y(k,6) + .0508005_r8*rxt(k,575)*y(k,116)
         mat(k,2311) = .2202005_r8*rxt(k,554)*y(k,7) + .0031005_r8*rxt(k,562)*y(k,105) &
                      + .0508005_r8*rxt(k,573)*y(k,117)
         mat(k,1040) = .2067005_r8*rxt(k,556)*y(k,142)
         mat(k,78) = .1792005_r8*rxt(k,555)*y(k,131) + .2067005_r8*rxt(k,554)*y(k,245)
         mat(k,84) = .0003005_r8*rxt(k,563)*y(k,131) + .0035005_r8*rxt(k,562)*y(k,245)
         mat(k,1012) = .1149005_r8*rxt(k,575)*y(k,142)
         mat(k,90) = .0082005_r8*rxt(k,574)*y(k,131) + .1149005_r8*rxt(k,573)*y(k,245)
         mat(k,2146) = .1792005_r8*rxt(k,555)*y(k,7) + .0003005_r8*rxt(k,563)*y(k,105) &
                      + .0082005_r8*rxt(k,574)*y(k,117)
         mat(k,2049) = .2067005_r8*rxt(k,556)*y(k,6) + .1149005_r8*rxt(k,575)*y(k,116)
         mat(k,2312) = .2067005_r8*rxt(k,554)*y(k,7) + .0035005_r8*rxt(k,562)*y(k,105) &
                      + .1149005_r8*rxt(k,573)*y(k,117)
         mat(k,1041) = .0653005_r8*rxt(k,556)*y(k,142)
         mat(k,79) = .0676005_r8*rxt(k,555)*y(k,131) + .0653005_r8*rxt(k,554)*y(k,245)
         mat(k,85) = .0073005_r8*rxt(k,563)*y(k,131) + .0003005_r8*rxt(k,562)*y(k,245)
         mat(k,1013) = .0348005_r8*rxt(k,575)*y(k,142)
         mat(k,91) = .0772005_r8*rxt(k,574)*y(k,131) + .0348005_r8*rxt(k,573)*y(k,245)
         mat(k,2147) = .0676005_r8*rxt(k,555)*y(k,7) + .0073005_r8*rxt(k,563)*y(k,105) &
                      + .0772005_r8*rxt(k,574)*y(k,117)
         mat(k,2050) = .0653005_r8*rxt(k,556)*y(k,6) + .0348005_r8*rxt(k,575)*y(k,116)
         mat(k,2313) = .0653005_r8*rxt(k,554)*y(k,7) + .0003005_r8*rxt(k,562)*y(k,105) &
                      + .0348005_r8*rxt(k,573)*y(k,117)
         mat(k,1042) = .1749305_r8*rxt(k,553)*y(k,133) + .1284005_r8*rxt(k,556) &
                      *y(k,142)
         mat(k,80) = .079_r8*rxt(k,555)*y(k,131) + .1284005_r8*rxt(k,554)*y(k,245)
         mat(k,921) = .0590245_r8*rxt(k,561)*y(k,133) + .0033005_r8*rxt(k,564) &
                      *y(k,142)
         mat(k,86) = .0057005_r8*rxt(k,563)*y(k,131) + .0271005_r8*rxt(k,562)*y(k,245)
         mat(k,1014) = .1749305_r8*rxt(k,572)*y(k,133) + .0554005_r8*rxt(k,575) &
                      *y(k,142)
         mat(k,92) = .0332005_r8*rxt(k,574)*y(k,131) + .0554005_r8*rxt(k,573)*y(k,245)
         mat(k,2148) = .079_r8*rxt(k,555)*y(k,7) + .0057005_r8*rxt(k,563)*y(k,105) &
                      + .0332005_r8*rxt(k,574)*y(k,117)
         mat(k,2245) = .1749305_r8*rxt(k,553)*y(k,6) + .0590245_r8*rxt(k,561)*y(k,100) &
                      + .1749305_r8*rxt(k,572)*y(k,116)
         mat(k,2051) = .1284005_r8*rxt(k,556)*y(k,6) + .0033005_r8*rxt(k,564)*y(k,100) &
                      + .0554005_r8*rxt(k,575)*y(k,116)
         mat(k,2314) = .1284005_r8*rxt(k,554)*y(k,7) + .0271005_r8*rxt(k,562)*y(k,105) &
                      + .0554005_r8*rxt(k,573)*y(k,117)
         mat(k,1043) = .5901905_r8*rxt(k,553)*y(k,133) + .114_r8*rxt(k,556)*y(k,142)
         mat(k,81) = .1254005_r8*rxt(k,555)*y(k,131) + .114_r8*rxt(k,554)*y(k,245)
         mat(k,922) = .0250245_r8*rxt(k,561)*y(k,133)
         mat(k,87) = .0623005_r8*rxt(k,563)*y(k,131) + .0474005_r8*rxt(k,562)*y(k,245)
         mat(k,1015) = .5901905_r8*rxt(k,572)*y(k,133) + .1278005_r8*rxt(k,575) &
                      *y(k,142)
         mat(k,93) = .130_r8*rxt(k,574)*y(k,131) + .1278005_r8*rxt(k,573)*y(k,245)
         mat(k,2149) = .1254005_r8*rxt(k,555)*y(k,7) + .0623005_r8*rxt(k,563)*y(k,105) &
                      + .130_r8*rxt(k,574)*y(k,117)
         mat(k,2246) = .5901905_r8*rxt(k,553)*y(k,6) + .0250245_r8*rxt(k,561)*y(k,100) &
                      + .5901905_r8*rxt(k,572)*y(k,116)
         mat(k,2052) = .114_r8*rxt(k,556)*y(k,6) + .1278005_r8*rxt(k,575)*y(k,116)
         mat(k,2315) = .114_r8*rxt(k,554)*y(k,7) + .0474005_r8*rxt(k,562)*y(k,105) &
                      + .1278005_r8*rxt(k,573)*y(k,117)
         mat(k,108) = .0097005_r8*rxt(k,560)*y(k,131) + .0023005_r8*rxt(k,559) &
                      *y(k,245)
         mat(k,100) = .1056005_r8*rxt(k,570)*y(k,131) + .2381005_r8*rxt(k,569) &
                      *y(k,245)
         mat(k,2153) = .0097005_r8*rxt(k,560)*y(k,9) + .1056005_r8*rxt(k,570)*y(k,110) &
                      + .0154005_r8*rxt(k,581)*y(k,222) + .0063005_r8*rxt(k,585) &
                      *y(k,226)
         mat(k,121) = .5931005_r8*rxt(k,578)*y(k,263)
         mat(k,127) = .0154005_r8*rxt(k,581)*y(k,131) + .1364005_r8*rxt(k,580) &
                      *y(k,245)
         mat(k,133) = .0063005_r8*rxt(k,585)*y(k,131) + .1677005_r8*rxt(k,584) &
                      *y(k,245)
         mat(k,2319) = .0023005_r8*rxt(k,559)*y(k,9) + .2381005_r8*rxt(k,569)*y(k,110) &
                      + .1364005_r8*rxt(k,580)*y(k,222) + .1677005_r8*rxt(k,584) &
                      *y(k,226)
         mat(k,1753) = .5931005_r8*rxt(k,578)*y(k,213)
         mat(k,109) = .0034005_r8*rxt(k,560)*y(k,131) + .0008005_r8*rxt(k,559) &
                      *y(k,245)
         mat(k,101) = .1026005_r8*rxt(k,570)*y(k,131) + .1308005_r8*rxt(k,569) &
                      *y(k,245)
         mat(k,2154) = .0034005_r8*rxt(k,560)*y(k,9) + .1026005_r8*rxt(k,570)*y(k,110) &
                      + .0452005_r8*rxt(k,581)*y(k,222) + .0237005_r8*rxt(k,585) &
                      *y(k,226)
         mat(k,122) = .1534005_r8*rxt(k,578)*y(k,263)
         mat(k,128) = .0452005_r8*rxt(k,581)*y(k,131) + .0101005_r8*rxt(k,580) &
                      *y(k,245)
         mat(k,134) = .0237005_r8*rxt(k,585)*y(k,131) + .0174005_r8*rxt(k,584) &
                      *y(k,245)
         mat(k,2320) = .0008005_r8*rxt(k,559)*y(k,9) + .1308005_r8*rxt(k,569)*y(k,110) &
                      + .0101005_r8*rxt(k,580)*y(k,222) + .0174005_r8*rxt(k,584) &
                      *y(k,226)
         mat(k,1754) = .1534005_r8*rxt(k,578)*y(k,213)
         mat(k,110) = .1579005_r8*rxt(k,560)*y(k,131) + .0843005_r8*rxt(k,559) &
                      *y(k,245)
         mat(k,102) = .0521005_r8*rxt(k,570)*y(k,131) + .0348005_r8*rxt(k,569) &
                      *y(k,245)
         mat(k,2155) = .1579005_r8*rxt(k,560)*y(k,9) + .0521005_r8*rxt(k,570)*y(k,110) &
                      + .0966005_r8*rxt(k,581)*y(k,222) + .0025005_r8*rxt(k,585) &
                      *y(k,226)
         mat(k,123) = .0459005_r8*rxt(k,578)*y(k,263)
         mat(k,129) = .0966005_r8*rxt(k,581)*y(k,131) + .0763005_r8*rxt(k,580) &
                      *y(k,245)
         mat(k,135) = .0025005_r8*rxt(k,585)*y(k,131) + .086_r8*rxt(k,584)*y(k,245)
         mat(k,2321) = .0843005_r8*rxt(k,559)*y(k,9) + .0348005_r8*rxt(k,569)*y(k,110) &
                      + .0763005_r8*rxt(k,580)*y(k,222) + .086_r8*rxt(k,584)*y(k,226)
         mat(k,1755) = .0459005_r8*rxt(k,578)*y(k,213)
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
         mat(k,111) = .0059005_r8*rxt(k,560)*y(k,131) + .0443005_r8*rxt(k,559) &
                      *y(k,245)
         mat(k,103) = .0143005_r8*rxt(k,570)*y(k,131) + .0076005_r8*rxt(k,569) &
                      *y(k,245)
         mat(k,2156) = .0059005_r8*rxt(k,560)*y(k,9) + .0143005_r8*rxt(k,570)*y(k,110) &
                      + .0073005_r8*rxt(k,581)*y(k,222) + .011_r8*rxt(k,585)*y(k,226)
         mat(k,124) = .0085005_r8*rxt(k,578)*y(k,263)
         mat(k,130) = .0073005_r8*rxt(k,581)*y(k,131) + .2157005_r8*rxt(k,580) &
                      *y(k,245)
         mat(k,136) = .011_r8*rxt(k,585)*y(k,131) + .0512005_r8*rxt(k,584)*y(k,245)
         mat(k,2322) = .0443005_r8*rxt(k,559)*y(k,9) + .0076005_r8*rxt(k,569)*y(k,110) &
                      + .2157005_r8*rxt(k,580)*y(k,222) + .0512005_r8*rxt(k,584) &
                      *y(k,226)
         mat(k,1756) = .0085005_r8*rxt(k,578)*y(k,213)
         mat(k,112) = .0536005_r8*rxt(k,560)*y(k,131) + .1621005_r8*rxt(k,559) &
                      *y(k,245)
         mat(k,104) = .0166005_r8*rxt(k,570)*y(k,131) + .0113005_r8*rxt(k,569) &
                      *y(k,245)
         mat(k,2157) = .0536005_r8*rxt(k,560)*y(k,9) + .0166005_r8*rxt(k,570)*y(k,110) &
                      + .238_r8*rxt(k,581)*y(k,222) + .1185005_r8*rxt(k,585)*y(k,226)
         mat(k,125) = .0128005_r8*rxt(k,578)*y(k,263)
         mat(k,131) = .238_r8*rxt(k,581)*y(k,131) + .0738005_r8*rxt(k,580)*y(k,245)
         mat(k,137) = .1185005_r8*rxt(k,585)*y(k,131) + .1598005_r8*rxt(k,584) &
                      *y(k,245)
         mat(k,2323) = .1621005_r8*rxt(k,559)*y(k,9) + .0113005_r8*rxt(k,569)*y(k,110) &
                      + .0738005_r8*rxt(k,580)*y(k,222) + .1598005_r8*rxt(k,584) &
                      *y(k,226)
         mat(k,1757) = .0128005_r8*rxt(k,578)*y(k,213)
         mat(k,120) = -(rxt(k,577)*y(k,263))
         mat(k,1761) = -rxt(k,577)*y(k,212)
         mat(k,126) = -(rxt(k,578)*y(k,263))
         mat(k,1762) = -rxt(k,578)*y(k,213)
         mat(k,233) = .100_r8*rxt(k,489)*y(k,263)
         mat(k,251) = .230_r8*rxt(k,491)*y(k,263)
         mat(k,1775) = .100_r8*rxt(k,489)*y(k,221) + .230_r8*rxt(k,491)*y(k,224)
         mat(k,726) = -(rxt(k,513)*y(k,263))
         mat(k,1841) = -rxt(k,513)*y(k,215)
         mat(k,2362) = rxt(k,511)*y(k,268)
         mat(k,1214) = rxt(k,511)*y(k,245)
         mat(k,682) = -(rxt(k,514)*y(k,263))
         mat(k,1836) = -rxt(k,514)*y(k,216)
         mat(k,2181) = .200_r8*rxt(k,507)*y(k,258) + .200_r8*rxt(k,517)*y(k,269)
         mat(k,1651) = .500_r8*rxt(k,505)*y(k,258)
         mat(k,1129) = .200_r8*rxt(k,507)*y(k,131) + .500_r8*rxt(k,505)*y(k,239)
         mat(k,1090) = .200_r8*rxt(k,517)*y(k,131)
         mat(k,530) = -(rxt(k,518)*y(k,263))
         mat(k,1819) = -rxt(k,518)*y(k,217)
         mat(k,2350) = rxt(k,516)*y(k,269)
         mat(k,1089) = rxt(k,516)*y(k,245)
         mat(k,1102) = -(rxt(k,519)*y(k,133) + rxt(k,520)*y(k,263))
         mat(k,2261) = -rxt(k,519)*y(k,218)
         mat(k,1869) = -rxt(k,520)*y(k,218)
         mat(k,1052) = .330_r8*rxt(k,500)*y(k,142)
         mat(k,1024) = .330_r8*rxt(k,503)*y(k,142)
         mat(k,2203) = .800_r8*rxt(k,507)*y(k,258) + .800_r8*rxt(k,517)*y(k,269)
         mat(k,2261) = mat(k,2261) + rxt(k,508)*y(k,258)
         mat(k,2071) = .330_r8*rxt(k,500)*y(k,6) + .330_r8*rxt(k,503)*y(k,116)
         mat(k,683) = rxt(k,514)*y(k,263)
         mat(k,1659) = .500_r8*rxt(k,505)*y(k,258) + rxt(k,515)*y(k,269)
         mat(k,1131) = .800_r8*rxt(k,507)*y(k,131) + rxt(k,508)*y(k,133) &
                      + .500_r8*rxt(k,505)*y(k,239)
         mat(k,1869) = mat(k,1869) + rxt(k,514)*y(k,216)
         mat(k,1093) = .800_r8*rxt(k,517)*y(k,131) + rxt(k,515)*y(k,239)
         mat(k,1145) = -(rxt(k,521)*y(k,263))
         mat(k,1873) = -rxt(k,521)*y(k,219)
         mat(k,1055) = .300_r8*rxt(k,500)*y(k,142)
         mat(k,1027) = .300_r8*rxt(k,503)*y(k,142)
         mat(k,2206) = .900_r8*rxt(k,512)*y(k,268)
         mat(k,2074) = .300_r8*rxt(k,500)*y(k,6) + .300_r8*rxt(k,503)*y(k,116)
         mat(k,1662) = rxt(k,510)*y(k,268)
         mat(k,1217) = .900_r8*rxt(k,512)*y(k,131) + rxt(k,510)*y(k,239)
         mat(k,653) = -(rxt(k,488)*y(k,263))
         mat(k,1832) = -rxt(k,488)*y(k,220)
         mat(k,2356) = rxt(k,486)*y(k,270)
         mat(k,790) = rxt(k,486)*y(k,245)
         mat(k,231) = -((rxt(k,489) + rxt(k,579)) * y(k,263))
         mat(k,1773) = -(rxt(k,489) + rxt(k,579)) * y(k,221)
         mat(k,132) = -(rxt(k,580)*y(k,245) + rxt(k,581)*y(k,131))
         mat(k,2326) = -rxt(k,580)*y(k,222)
         mat(k,2160) = -rxt(k,581)*y(k,222)
         mat(k,230) = rxt(k,579)*y(k,263)
         mat(k,1763) = rxt(k,579)*y(k,221)
         mat(k,247) = -(rxt(k,455)*y(k,263))
         mat(k,1776) = -rxt(k,455)*y(k,223)
         mat(k,2329) = rxt(k,452)*y(k,271)
         mat(k,1269) = rxt(k,452)*y(k,245)
         mat(k,252) = -(rxt(k,491)*y(k,263))
         mat(k,1777) = -rxt(k,491)*y(k,224)
         mat(k,764) = -(rxt(k,494)*y(k,263))
         mat(k,1844) = -rxt(k,494)*y(k,225)
         mat(k,2365) = rxt(k,492)*y(k,272)
         mat(k,807) = rxt(k,492)*y(k,245)
         mat(k,138) = -(rxt(k,584)*y(k,245) + rxt(k,585)*y(k,131))
         mat(k,2327) = -rxt(k,584)*y(k,226)
         mat(k,2161) = -rxt(k,585)*y(k,226)
         mat(k,250) = rxt(k,583)*y(k,263)
         mat(k,1764) = rxt(k,583)*y(k,224)
         mat(k,260) = -(rxt(k,497)*y(k,263))
         mat(k,1778) = -rxt(k,497)*y(k,227)
         mat(k,253) = .150_r8*rxt(k,491)*y(k,263)
         mat(k,1778) = mat(k,1778) + .150_r8*rxt(k,491)*y(k,224)
         mat(k,477) = -(rxt(k,498)*y(k,263))
         mat(k,1812) = -rxt(k,498)*y(k,228)
         mat(k,2342) = rxt(k,495)*y(k,273)
         mat(k,553) = rxt(k,495)*y(k,245)
         mat(k,567) = -(rxt(k,456)*y(k,245) + rxt(k,457)*y(k,131) + rxt(k,485) &
                      *y(k,132))
         mat(k,2353) = -rxt(k,456)*y(k,231)
         mat(k,2176) = -rxt(k,457)*y(k,231)
         mat(k,2456) = -rxt(k,485)*y(k,231)
         mat(k,275) = rxt(k,462)*y(k,263)
         mat(k,1823) = rxt(k,462)*y(k,24)
         mat(k,1072) = -(rxt(k,417)*y(k,245) + (rxt(k,418) + rxt(k,419)) * y(k,131))
         mat(k,2380) = -rxt(k,417)*y(k,232)
         mat(k,2200) = -(rxt(k,418) + rxt(k,419)) * y(k,232)
         mat(k,700) = rxt(k,420)*y(k,263)
         mat(k,266) = rxt(k,421)*y(k,263)
         mat(k,1866) = rxt(k,420)*y(k,2) + rxt(k,421)*y(k,17)
         mat(k,539) = -(rxt(k,459)*y(k,245) + rxt(k,460)*y(k,131))
         mat(k,2351) = -rxt(k,459)*y(k,233)
         mat(k,2173) = -rxt(k,460)*y(k,233)
         mat(k,207) = .350_r8*rxt(k,458)*y(k,263)
         mat(k,467) = rxt(k,461)*y(k,263)
         mat(k,1820) = .350_r8*rxt(k,458)*y(k,8) + rxt(k,461)*y(k,10)
         mat(k,485) = -(rxt(k,463)*y(k,245) + rxt(k,465)*y(k,131))
         mat(k,2343) = -rxt(k,463)*y(k,234)
         mat(k,2167) = -rxt(k,465)*y(k,234)
         mat(k,372) = rxt(k,464)*y(k,263)
         mat(k,234) = .070_r8*rxt(k,489)*y(k,263)
         mat(k,254) = .060_r8*rxt(k,491)*y(k,263)
         mat(k,1813) = rxt(k,464)*y(k,25) + .070_r8*rxt(k,489)*y(k,221) &
                      + .060_r8*rxt(k,491)*y(k,224)
         mat(k,953) = -(4._r8*rxt(k,339)*y(k,235) + rxt(k,340)*y(k,239) + rxt(k,341) &
                      *y(k,245) + rxt(k,342)*y(k,131))
         mat(k,1655) = -rxt(k,340)*y(k,235)
         mat(k,2376) = -rxt(k,341)*y(k,235)
         mat(k,2195) = -rxt(k,342)*y(k,235)
         mat(k,377) = .500_r8*rxt(k,344)*y(k,263)
         mat(k,324) = rxt(k,345)*y(k,58) + rxt(k,346)*y(k,263)
         mat(k,1971) = rxt(k,345)*y(k,30)
         mat(k,1858) = .500_r8*rxt(k,344)*y(k,29) + rxt(k,346)*y(k,30)
         mat(k,977) = -(rxt(k,368)*y(k,239) + rxt(k,369)*y(k,245) + rxt(k,370) &
                      *y(k,131))
         mat(k,1656) = -rxt(k,368)*y(k,236)
         mat(k,2379) = -rxt(k,369)*y(k,236)
         mat(k,2198) = -rxt(k,370)*y(k,236)
         mat(k,436) = rxt(k,371)*y(k,263)
         mat(k,330) = rxt(k,375)*y(k,58) + rxt(k,372)*y(k,263)
         mat(k,1973) = rxt(k,375)*y(k,33)
         mat(k,1861) = rxt(k,371)*y(k,32) + rxt(k,372)*y(k,33)
         mat(k,690) = -(rxt(k,466)*y(k,245) + rxt(k,467)*y(k,131))
         mat(k,2359) = -rxt(k,466)*y(k,237)
         mat(k,2182) = -rxt(k,467)*y(k,237)
         mat(k,315) = rxt(k,468)*y(k,263)
         mat(k,2182) = mat(k,2182) + rxt(k,457)*y(k,231)
         mat(k,2059) = rxt(k,483)*y(k,149)
         mat(k,527) = rxt(k,483)*y(k,142)
         mat(k,568) = rxt(k,457)*y(k,131) + .400_r8*rxt(k,456)*y(k,245)
         mat(k,2359) = mat(k,2359) + .400_r8*rxt(k,456)*y(k,231)
         mat(k,1837) = rxt(k,468)*y(k,34)
         mat(k,1470) = -(4._r8*rxt(k,350)*y(k,238) + rxt(k,351)*y(k,239) + rxt(k,352) &
                      *y(k,245) + rxt(k,353)*y(k,131) + rxt(k,364)*y(k,132) + rxt(k,392) &
                      *y(k,249) + rxt(k,425)*y(k,247) + rxt(k,430)*y(k,248) + rxt(k,439) &
                      *y(k,103) + rxt(k,450)*y(k,271))
         mat(k,1679) = -rxt(k,351)*y(k,238)
         mat(k,2402) = -rxt(k,352)*y(k,238)
         mat(k,2224) = -rxt(k,353)*y(k,238)
         mat(k,2473) = -rxt(k,364)*y(k,238)
         mat(k,1400) = -rxt(k,392)*y(k,238)
         mat(k,1345) = -rxt(k,425)*y(k,238)
         mat(k,1378) = -rxt(k,430)*y(k,238)
         mat(k,1299) = -rxt(k,439)*y(k,238)
         mat(k,1277) = -rxt(k,450)*y(k,238)
         mat(k,1059) = .060_r8*rxt(k,500)*y(k,142)
         mat(k,1195) = rxt(k,347)*y(k,133) + rxt(k,348)*y(k,263)
         mat(k,1324) = rxt(k,373)*y(k,133) + rxt(k,374)*y(k,263)
         mat(k,636) = .500_r8*rxt(k,355)*y(k,263)
         mat(k,933) = .080_r8*rxt(k,445)*y(k,142)
         mat(k,1315) = .100_r8*rxt(k,398)*y(k,142)
         mat(k,1031) = .060_r8*rxt(k,503)*y(k,142)
         mat(k,1420) = .280_r8*rxt(k,412)*y(k,142)
         mat(k,2224) = mat(k,2224) + .530_r8*rxt(k,396)*y(k,249) + rxt(k,405)*y(k,251) &
                      + rxt(k,408)*y(k,253) + rxt(k,383)*y(k,267)
         mat(k,2283) = rxt(k,347)*y(k,47) + rxt(k,373)*y(k,51) + .530_r8*rxt(k,395) &
                      *y(k,249) + rxt(k,406)*y(k,251)
         mat(k,2090) = .060_r8*rxt(k,500)*y(k,6) + .080_r8*rxt(k,445)*y(k,100) &
                      + .100_r8*rxt(k,398)*y(k,111) + .060_r8*rxt(k,503)*y(k,116) &
                      + .280_r8*rxt(k,412)*y(k,118)
         mat(k,1148) = .650_r8*rxt(k,521)*y(k,263)
         mat(k,1470) = mat(k,1470) + .530_r8*rxt(k,392)*y(k,249)
         mat(k,1679) = mat(k,1679) + .260_r8*rxt(k,393)*y(k,249) + rxt(k,402)*y(k,251) &
                      + .300_r8*rxt(k,381)*y(k,267)
         mat(k,2402) = mat(k,2402) + .450_r8*rxt(k,403)*y(k,251) + .200_r8*rxt(k,407) &
                      *y(k,253) + .150_r8*rxt(k,382)*y(k,267)
         mat(k,1400) = mat(k,1400) + .530_r8*rxt(k,396)*y(k,131) + .530_r8*rxt(k,395) &
                      *y(k,133) + .530_r8*rxt(k,392)*y(k,238) + .260_r8*rxt(k,393) &
                      *y(k,239)
         mat(k,1440) = rxt(k,405)*y(k,131) + rxt(k,406)*y(k,133) + rxt(k,402)*y(k,239) &
                      + .450_r8*rxt(k,403)*y(k,245) + 4.000_r8*rxt(k,404)*y(k,251)
         mat(k,721) = rxt(k,408)*y(k,131) + .200_r8*rxt(k,407)*y(k,245)
         mat(k,1892) = rxt(k,348)*y(k,47) + rxt(k,374)*y(k,51) + .500_r8*rxt(k,355) &
                      *y(k,53) + .650_r8*rxt(k,521)*y(k,219)
         mat(k,1261) = rxt(k,383)*y(k,131) + .300_r8*rxt(k,381)*y(k,239) &
                      + .150_r8*rxt(k,382)*y(k,245)
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
         mat(k,1684) = -(rxt(k,241)*y(k,61) + (4._r8*rxt(k,318) + 4._r8*rxt(k,319) &
                      ) * y(k,239) + rxt(k,320)*y(k,245) + rxt(k,321)*y(k,131) &
                      + rxt(k,340)*y(k,235) + rxt(k,351)*y(k,238) + rxt(k,368) &
                      *y(k,236) + rxt(k,381)*y(k,267) + rxt(k,393)*y(k,249) + rxt(k,402) &
                      *y(k,251) + rxt(k,426)*y(k,247) + rxt(k,431)*y(k,248) + rxt(k,440) &
                      *y(k,103) + rxt(k,451)*y(k,271) + rxt(k,505)*y(k,258) + rxt(k,510) &
                      *y(k,268) + rxt(k,515)*y(k,269))
         mat(k,2124) = -rxt(k,241)*y(k,239)
         mat(k,2410) = -rxt(k,320)*y(k,239)
         mat(k,2230) = -rxt(k,321)*y(k,239)
         mat(k,956) = -rxt(k,340)*y(k,239)
         mat(k,1474) = -rxt(k,351)*y(k,239)
         mat(k,981) = -rxt(k,368)*y(k,239)
         mat(k,1262) = -rxt(k,381)*y(k,239)
         mat(k,1402) = -rxt(k,393)*y(k,239)
         mat(k,1442) = -rxt(k,402)*y(k,239)
         mat(k,1347) = -rxt(k,426)*y(k,239)
         mat(k,1380) = -rxt(k,431)*y(k,239)
         mat(k,1301) = -rxt(k,440)*y(k,239)
         mat(k,1279) = -rxt(k,451)*y(k,239)
         mat(k,1136) = -rxt(k,505)*y(k,239)
         mat(k,1224) = -rxt(k,510)*y(k,239)
         mat(k,1095) = -rxt(k,515)*y(k,239)
         mat(k,1184) = .280_r8*rxt(k,367)*y(k,142)
         mat(k,749) = rxt(k,354)*y(k,263)
         mat(k,442) = .700_r8*rxt(k,323)*y(k,263)
         mat(k,1632) = rxt(k,235)*y(k,58) + rxt(k,291)*y(k,75) + rxt(k,330)*y(k,259) &
                      + rxt(k,324)*y(k,263)
         mat(k,1990) = rxt(k,235)*y(k,56)
         mat(k,944) = rxt(k,291)*y(k,56)
         mat(k,934) = .050_r8*rxt(k,445)*y(k,142)
         mat(k,1301) = mat(k,1301) + rxt(k,439)*y(k,238)
         mat(k,2230) = mat(k,2230) + rxt(k,353)*y(k,238) + .830_r8*rxt(k,471)*y(k,240) &
                      + .170_r8*rxt(k,477)*y(k,252)
         mat(k,2097) = .280_r8*rxt(k,367)*y(k,31) + .050_r8*rxt(k,445)*y(k,100)
         mat(k,1474) = mat(k,1474) + rxt(k,439)*y(k,103) + rxt(k,353)*y(k,131) &
                      + 4.000_r8*rxt(k,350)*y(k,238) + .900_r8*rxt(k,351)*y(k,239) &
                      + .450_r8*rxt(k,352)*y(k,245) + rxt(k,425)*y(k,247) + rxt(k,430) &
                      *y(k,248) + rxt(k,392)*y(k,249) + rxt(k,401)*y(k,251) &
                      + rxt(k,450)*y(k,271)
         mat(k,1684) = mat(k,1684) + .900_r8*rxt(k,351)*y(k,238)
         mat(k,823) = .830_r8*rxt(k,471)*y(k,131) + .330_r8*rxt(k,470)*y(k,245)
         mat(k,2410) = mat(k,2410) + .450_r8*rxt(k,352)*y(k,238) + .330_r8*rxt(k,470) &
                      *y(k,240) + .070_r8*rxt(k,476)*y(k,252)
         mat(k,1347) = mat(k,1347) + rxt(k,425)*y(k,238)
         mat(k,1380) = mat(k,1380) + rxt(k,430)*y(k,238)
         mat(k,1402) = mat(k,1402) + rxt(k,392)*y(k,238)
         mat(k,1442) = mat(k,1442) + rxt(k,401)*y(k,238)
         mat(k,971) = .170_r8*rxt(k,477)*y(k,131) + .070_r8*rxt(k,476)*y(k,245)
         mat(k,2033) = rxt(k,330)*y(k,56)
         mat(k,1901) = rxt(k,354)*y(k,52) + .700_r8*rxt(k,323)*y(k,55) + rxt(k,324) &
                      *y(k,56)
         mat(k,1279) = mat(k,1279) + rxt(k,450)*y(k,238)
         mat(k,820) = -(rxt(k,470)*y(k,245) + rxt(k,471)*y(k,131) + rxt(k,472) &
                      *y(k,132))
         mat(k,2369) = -rxt(k,470)*y(k,240)
         mat(k,2188) = -rxt(k,471)*y(k,240)
         mat(k,2462) = -rxt(k,472)*y(k,240)
         mat(k,907) = -(rxt(k,603)*y(k,256) + rxt(k,604)*y(k,262) + rxt(k,605) &
                      *y(k,255))
         mat(k,888) = -rxt(k,603)*y(k,241)
         mat(k,896) = -rxt(k,604)*y(k,241)
         mat(k,741) = -rxt(k,605)*y(k,241)
         mat(k,618) = -((rxt(k,389) + rxt(k,390)) * y(k,131))
         mat(k,2178) = -(rxt(k,389) + rxt(k,390)) * y(k,242)
         mat(k,408) = rxt(k,388)*y(k,263)
         mat(k,1828) = rxt(k,388)*y(k,18)
         mat(k,507) = -(rxt(k,359)*y(k,141))
         mat(k,1568) = -rxt(k,359)*y(k,243)
         mat(k,2171) = .750_r8*rxt(k,357)*y(k,244)
         mat(k,838) = .750_r8*rxt(k,357)*y(k,131)
         mat(k,839) = -(rxt(k,356)*y(k,245) + rxt(k,357)*y(k,131))
         mat(k,2371) = -rxt(k,356)*y(k,244)
         mat(k,2189) = -rxt(k,357)*y(k,244)
         mat(k,579) = rxt(k,363)*y(k,263)
         mat(k,1850) = rxt(k,363)*y(k,27)
         mat(k,2421) = -((rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,78) + rxt(k,196) &
                      *y(k,140) + rxt(k,197)*y(k,142) + rxt(k,201)*y(k,263) &
                      + 4._r8*rxt(k,206)*y(k,245) + rxt(k,218)*y(k,133) + rxt(k,223) &
                      *y(k,131) + rxt(k,228)*y(k,132) + (rxt(k,238) + rxt(k,239) &
                      ) * y(k,58) + rxt(k,245)*y(k,61) + rxt(k,271)*y(k,19) + rxt(k,277) &
                      *y(k,21) + rxt(k,314)*y(k,44) + rxt(k,320)*y(k,239) + rxt(k,327) &
                      *y(k,246) + rxt(k,341)*y(k,235) + rxt(k,352)*y(k,238) + rxt(k,356) &
                      *y(k,244) + rxt(k,369)*y(k,236) + rxt(k,378)*y(k,266) + rxt(k,382) &
                      *y(k,267) + rxt(k,394)*y(k,249) + rxt(k,403)*y(k,251) + rxt(k,407) &
                      *y(k,253) + rxt(k,417)*y(k,232) + rxt(k,427)*y(k,247) + rxt(k,432) &
                      *y(k,248) + rxt(k,441)*y(k,103) + rxt(k,452)*y(k,271) + rxt(k,456) &
                      *y(k,231) + rxt(k,459)*y(k,233) + rxt(k,463)*y(k,234) + rxt(k,466) &
                      *y(k,237) + rxt(k,470)*y(k,240) + rxt(k,473)*y(k,250) + rxt(k,476) &
                      *y(k,252) + rxt(k,479)*y(k,265) + rxt(k,486)*y(k,270) + rxt(k,492) &
                      *y(k,272) + rxt(k,495)*y(k,273) + rxt(k,506)*y(k,258) + rxt(k,511) &
                      *y(k,268) + rxt(k,516)*y(k,269))
         mat(k,1739) = -(rxt(k,192) + rxt(k,193) + rxt(k,194)) * y(k,245)
         mat(k,1955) = -rxt(k,196)*y(k,245)
         mat(k,2108) = -rxt(k,197)*y(k,245)
         mat(k,1912) = -rxt(k,201)*y(k,245)
         mat(k,2301) = -rxt(k,218)*y(k,245)
         mat(k,2241) = -rxt(k,223)*y(k,245)
         mat(k,2492) = -rxt(k,228)*y(k,245)
         mat(k,2001) = -(rxt(k,238) + rxt(k,239)) * y(k,245)
         mat(k,2135) = -rxt(k,245)*y(k,245)
         mat(k,1562) = -rxt(k,271)*y(k,245)
         mat(k,1618) = -rxt(k,277)*y(k,245)
         mat(k,2447) = -rxt(k,314)*y(k,245)
         mat(k,1695) = -rxt(k,320)*y(k,245)
         mat(k,494) = -rxt(k,327)*y(k,245)
         mat(k,959) = -rxt(k,341)*y(k,245)
         mat(k,1480) = -rxt(k,352)*y(k,245)
         mat(k,844) = -rxt(k,356)*y(k,245)
         mat(k,984) = -rxt(k,369)*y(k,245)
         mat(k,866) = -rxt(k,378)*y(k,245)
         mat(k,1265) = -rxt(k,382)*y(k,245)
         mat(k,1407) = -rxt(k,394)*y(k,245)
         mat(k,1448) = -rxt(k,403)*y(k,245)
         mat(k,724) = -rxt(k,407)*y(k,245)
         mat(k,1080) = -rxt(k,417)*y(k,245)
         mat(k,1353) = -rxt(k,427)*y(k,245)
         mat(k,1386) = -rxt(k,432)*y(k,245)
         mat(k,1306) = -rxt(k,441)*y(k,245)
         mat(k,1283) = -rxt(k,452)*y(k,245)
         mat(k,571) = -rxt(k,456)*y(k,245)
         mat(k,544) = -rxt(k,459)*y(k,245)
         mat(k,489) = -rxt(k,463)*y(k,245)
         mat(k,694) = -rxt(k,466)*y(k,245)
         mat(k,826) = -rxt(k,470)*y(k,245)
         mat(k,786) = -rxt(k,473)*y(k,245)
         mat(k,974) = -rxt(k,476)*y(k,245)
         mat(k,502) = -rxt(k,479)*y(k,245)
         mat(k,801) = -rxt(k,486)*y(k,245)
         mat(k,818) = -rxt(k,492)*y(k,245)
         mat(k,559) = -rxt(k,495)*y(k,245)
         mat(k,1141) = -rxt(k,506)*y(k,245)
         mat(k,1228) = -rxt(k,511)*y(k,245)
         mat(k,1099) = -rxt(k,516)*y(k,245)
         mat(k,1064) = .570_r8*rxt(k,500)*y(k,142)
         mat(k,209) = .650_r8*rxt(k,458)*y(k,263)
         mat(k,1562) = mat(k,1562) + rxt(k,270)*y(k,44)
         mat(k,1618) = mat(k,1618) + rxt(k,282)*y(k,263)
         mat(k,322) = .350_r8*rxt(k,336)*y(k,263)
         mat(k,584) = .130_r8*rxt(k,338)*y(k,142)
         mat(k,312) = rxt(k,343)*y(k,263)
         mat(k,1189) = .280_r8*rxt(k,367)*y(k,142)
         mat(k,2447) = mat(k,2447) + rxt(k,270)*y(k,19) + rxt(k,234)*y(k,58) &
                      + rxt(k,315)*y(k,133) + rxt(k,316)*y(k,140)
         mat(k,633) = rxt(k,299)*y(k,58) + rxt(k,300)*y(k,263)
         mat(k,405) = rxt(k,302)*y(k,58) + rxt(k,303)*y(k,263)
         mat(k,147) = rxt(k,349)*y(k,263)
         mat(k,857) = rxt(k,322)*y(k,263)
         mat(k,1643) = rxt(k,331)*y(k,259)
         mat(k,2001) = mat(k,2001) + rxt(k,234)*y(k,44) + rxt(k,299)*y(k,45) &
                      + rxt(k,302)*y(k,48) + rxt(k,237)*y(k,81)
         mat(k,2135) = mat(k,2135) + rxt(k,241)*y(k,239) + rxt(k,252)*y(k,263)
         mat(k,1205) = rxt(k,334)*y(k,263)
         mat(k,242) = .730_r8*rxt(k,469)*y(k,263)
         mat(k,352) = .500_r8*rxt(k,537)*y(k,263)
         mat(k,1211) = rxt(k,360)*y(k,263)
         mat(k,1088) = rxt(k,361)*y(k,263)
         mat(k,1739) = mat(k,1739) + rxt(k,195)*y(k,141)
         mat(k,676) = rxt(k,237)*y(k,58) + rxt(k,191)*y(k,140) + rxt(k,200)*y(k,263)
         mat(k,229) = rxt(k,325)*y(k,263)
         mat(k,965) = rxt(k,326)*y(k,263)
         mat(k,1246) = rxt(k,391)*y(k,263)
         mat(k,1254) = rxt(k,376)*y(k,263)
         mat(k,938) = .370_r8*rxt(k,445)*y(k,142)
         mat(k,648) = .300_r8*rxt(k,436)*y(k,263)
         mat(k,617) = rxt(k,437)*y(k,263)
         mat(k,1306) = mat(k,1306) + rxt(k,442)*y(k,131) + rxt(k,443)*y(k,133) &
                      + rxt(k,439)*y(k,238) + 1.200_r8*rxt(k,440)*y(k,239)
         mat(k,475) = rxt(k,444)*y(k,263)
         mat(k,1319) = .140_r8*rxt(k,398)*y(k,142)
         mat(k,388) = .200_r8*rxt(k,400)*y(k,263)
         mat(k,668) = .500_r8*rxt(k,411)*y(k,263)
         mat(k,1036) = .570_r8*rxt(k,503)*y(k,142)
         mat(k,1430) = .280_r8*rxt(k,412)*y(k,142)
         mat(k,464) = rxt(k,448)*y(k,263)
         mat(k,1168) = rxt(k,449)*y(k,263)
         mat(k,2241) = mat(k,2241) + rxt(k,442)*y(k,103) + rxt(k,418)*y(k,232) &
                      + rxt(k,460)*y(k,233) + rxt(k,465)*y(k,234) + rxt(k,342) &
                      *y(k,235) + rxt(k,370)*y(k,236) + rxt(k,321)*y(k,239) &
                      + .170_r8*rxt(k,471)*y(k,240) + rxt(k,389)*y(k,242) &
                      + .250_r8*rxt(k,357)*y(k,244) + rxt(k,329)*y(k,246) &
                      + .920_r8*rxt(k,428)*y(k,247) + .920_r8*rxt(k,434)*y(k,248) &
                      + .470_r8*rxt(k,396)*y(k,249) + .400_r8*rxt(k,474)*y(k,250) &
                      + .830_r8*rxt(k,477)*y(k,252) + rxt(k,480)*y(k,265) + rxt(k,379) &
                      *y(k,266) + .900_r8*rxt(k,512)*y(k,268) + .800_r8*rxt(k,517) &
                      *y(k,269) + rxt(k,487)*y(k,270) + rxt(k,453)*y(k,271) &
                      + rxt(k,493)*y(k,272) + rxt(k,496)*y(k,273)
         mat(k,2301) = mat(k,2301) + rxt(k,315)*y(k,44) + rxt(k,443)*y(k,103) &
                      + rxt(k,429)*y(k,247) + rxt(k,435)*y(k,248) + .470_r8*rxt(k,395) &
                      *y(k,249) + rxt(k,221)*y(k,263) + rxt(k,454)*y(k,271)
         mat(k,1955) = mat(k,1955) + rxt(k,316)*y(k,44) + rxt(k,191)*y(k,81)
         mat(k,1594) = rxt(k,195)*y(k,78) + rxt(k,359)*y(k,243)
         mat(k,2108) = mat(k,2108) + .570_r8*rxt(k,500)*y(k,6) + .130_r8*rxt(k,338) &
                      *y(k,27) + .280_r8*rxt(k,367)*y(k,31) + .370_r8*rxt(k,445) &
                      *y(k,100) + .140_r8*rxt(k,398)*y(k,111) + .570_r8*rxt(k,503) &
                      *y(k,116) + .280_r8*rxt(k,412)*y(k,118) + rxt(k,203)*y(k,263)
         mat(k,218) = .800_r8*rxt(k,481)*y(k,263)
         mat(k,993) = rxt(k,527)*y(k,263)
         mat(k,1152) = .200_r8*rxt(k,521)*y(k,263)
         mat(k,237) = .280_r8*rxt(k,489)*y(k,263)
         mat(k,259) = .380_r8*rxt(k,491)*y(k,263)
         mat(k,264) = .630_r8*rxt(k,497)*y(k,263)
         mat(k,1080) = mat(k,1080) + rxt(k,418)*y(k,131)
         mat(k,544) = mat(k,544) + rxt(k,460)*y(k,131)
         mat(k,489) = mat(k,489) + rxt(k,465)*y(k,131)
         mat(k,959) = mat(k,959) + rxt(k,342)*y(k,131) + 2.400_r8*rxt(k,339)*y(k,235) &
                      + rxt(k,340)*y(k,239)
         mat(k,984) = mat(k,984) + rxt(k,370)*y(k,131) + rxt(k,368)*y(k,239)
         mat(k,1480) = mat(k,1480) + rxt(k,439)*y(k,103) + .900_r8*rxt(k,351)*y(k,239) &
                      + rxt(k,425)*y(k,247) + rxt(k,430)*y(k,248) + .470_r8*rxt(k,392) &
                      *y(k,249) + rxt(k,450)*y(k,271)
         mat(k,1695) = mat(k,1695) + rxt(k,241)*y(k,61) + 1.200_r8*rxt(k,440)*y(k,103) &
                      + rxt(k,321)*y(k,131) + rxt(k,340)*y(k,235) + rxt(k,368) &
                      *y(k,236) + .900_r8*rxt(k,351)*y(k,238) + 4.000_r8*rxt(k,318) &
                      *y(k,239) + rxt(k,426)*y(k,247) + rxt(k,431)*y(k,248) &
                      + .730_r8*rxt(k,393)*y(k,249) + rxt(k,402)*y(k,251) &
                      + .500_r8*rxt(k,505)*y(k,258) + .300_r8*rxt(k,381)*y(k,267) &
                      + rxt(k,510)*y(k,268) + rxt(k,515)*y(k,269) + .800_r8*rxt(k,451) &
                      *y(k,271)
         mat(k,826) = mat(k,826) + .170_r8*rxt(k,471)*y(k,131) + .070_r8*rxt(k,470) &
                      *y(k,245)
         mat(k,623) = rxt(k,389)*y(k,131)
         mat(k,510) = rxt(k,359)*y(k,141)
         mat(k,844) = mat(k,844) + .250_r8*rxt(k,357)*y(k,131)
         mat(k,2421) = mat(k,2421) + .070_r8*rxt(k,470)*y(k,240) + .160_r8*rxt(k,473) &
                      *y(k,250) + .330_r8*rxt(k,476)*y(k,252)
         mat(k,494) = mat(k,494) + rxt(k,329)*y(k,131)
         mat(k,1353) = mat(k,1353) + .920_r8*rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133) &
                      + rxt(k,425)*y(k,238) + rxt(k,426)*y(k,239)
         mat(k,1386) = mat(k,1386) + .920_r8*rxt(k,434)*y(k,131) + rxt(k,435)*y(k,133) &
                      + rxt(k,430)*y(k,238) + rxt(k,431)*y(k,239)
         mat(k,1407) = mat(k,1407) + .470_r8*rxt(k,396)*y(k,131) + .470_r8*rxt(k,395) &
                      *y(k,133) + .470_r8*rxt(k,392)*y(k,238) + .730_r8*rxt(k,393) &
                      *y(k,239)
         mat(k,786) = mat(k,786) + .400_r8*rxt(k,474)*y(k,131) + .160_r8*rxt(k,473) &
                      *y(k,245)
         mat(k,1448) = mat(k,1448) + rxt(k,402)*y(k,239)
         mat(k,974) = mat(k,974) + .830_r8*rxt(k,477)*y(k,131) + .330_r8*rxt(k,476) &
                      *y(k,245)
         mat(k,1141) = mat(k,1141) + .500_r8*rxt(k,505)*y(k,239)
         mat(k,2044) = rxt(k,331)*y(k,56)
         mat(k,1912) = mat(k,1912) + .650_r8*rxt(k,458)*y(k,8) + rxt(k,282)*y(k,21) &
                      + .350_r8*rxt(k,336)*y(k,26) + rxt(k,343)*y(k,28) + rxt(k,300) &
                      *y(k,45) + rxt(k,303)*y(k,48) + rxt(k,349)*y(k,49) + rxt(k,322) &
                      *y(k,54) + rxt(k,252)*y(k,61) + rxt(k,334)*y(k,64) &
                      + .730_r8*rxt(k,469)*y(k,68) + .500_r8*rxt(k,537)*y(k,69) &
                      + rxt(k,360)*y(k,76) + rxt(k,361)*y(k,77) + rxt(k,200)*y(k,81) &
                      + rxt(k,325)*y(k,88) + rxt(k,326)*y(k,89) + rxt(k,391)*y(k,95) &
                      + rxt(k,376)*y(k,97) + .300_r8*rxt(k,436)*y(k,101) + rxt(k,437) &
                      *y(k,102) + rxt(k,444)*y(k,104) + .200_r8*rxt(k,400)*y(k,112) &
                      + .500_r8*rxt(k,411)*y(k,115) + rxt(k,448)*y(k,122) + rxt(k,449) &
                      *y(k,123) + rxt(k,221)*y(k,133) + rxt(k,203)*y(k,142) &
                      + .800_r8*rxt(k,481)*y(k,150) + rxt(k,527)*y(k,161) &
                      + .200_r8*rxt(k,521)*y(k,219) + .280_r8*rxt(k,489)*y(k,221) &
                      + .380_r8*rxt(k,491)*y(k,224) + .630_r8*rxt(k,497)*y(k,227)
         mat(k,502) = mat(k,502) + rxt(k,480)*y(k,131)
         mat(k,866) = mat(k,866) + rxt(k,379)*y(k,131)
         mat(k,1265) = mat(k,1265) + .300_r8*rxt(k,381)*y(k,239)
         mat(k,1228) = mat(k,1228) + .900_r8*rxt(k,512)*y(k,131) + rxt(k,510)*y(k,239)
         mat(k,1099) = mat(k,1099) + .800_r8*rxt(k,517)*y(k,131) + rxt(k,515)*y(k,239)
         mat(k,801) = mat(k,801) + rxt(k,487)*y(k,131)
         mat(k,1283) = mat(k,1283) + rxt(k,453)*y(k,131) + rxt(k,454)*y(k,133) &
                      + rxt(k,450)*y(k,238) + .800_r8*rxt(k,451)*y(k,239)
         mat(k,818) = mat(k,818) + rxt(k,493)*y(k,131)
         mat(k,559) = mat(k,559) + rxt(k,496)*y(k,131)
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
         mat(k,491) = -(rxt(k,327)*y(k,245) + rxt(k,329)*y(k,131))
         mat(k,2344) = -rxt(k,327)*y(k,246)
         mat(k,2168) = -rxt(k,329)*y(k,246)
         mat(k,2425) = rxt(k,314)*y(k,245)
         mat(k,2344) = mat(k,2344) + rxt(k,314)*y(k,44)
         mat(k,1341) = -(rxt(k,425)*y(k,238) + rxt(k,426)*y(k,239) + rxt(k,427) &
                      *y(k,245) + rxt(k,428)*y(k,131) + rxt(k,429)*y(k,133))
         mat(k,1465) = -rxt(k,425)*y(k,247)
         mat(k,1674) = -rxt(k,426)*y(k,247)
         mat(k,2397) = -rxt(k,427)*y(k,247)
         mat(k,2219) = -rxt(k,428)*y(k,247)
         mat(k,2278) = -rxt(k,429)*y(k,247)
         mat(k,930) = .600_r8*rxt(k,446)*y(k,263)
         mat(k,1887) = .600_r8*rxt(k,446)*y(k,100)
         mat(k,1374) = -(rxt(k,430)*y(k,238) + rxt(k,431)*y(k,239) + rxt(k,432) &
                      *y(k,245) + rxt(k,434)*y(k,131) + rxt(k,435)*y(k,133))
         mat(k,1466) = -rxt(k,430)*y(k,248)
         mat(k,1675) = -rxt(k,431)*y(k,248)
         mat(k,2398) = -rxt(k,432)*y(k,248)
         mat(k,2220) = -rxt(k,434)*y(k,248)
         mat(k,2279) = -rxt(k,435)*y(k,248)
         mat(k,931) = .400_r8*rxt(k,446)*y(k,263)
         mat(k,1888) = .400_r8*rxt(k,446)*y(k,100)
         mat(k,1398) = -(rxt(k,392)*y(k,238) + rxt(k,393)*y(k,239) + rxt(k,394) &
                      *y(k,245) + rxt(k,395)*y(k,133) + (rxt(k,396) + rxt(k,397) &
                      ) * y(k,131))
         mat(k,1467) = -rxt(k,392)*y(k,249)
         mat(k,1676) = -rxt(k,393)*y(k,249)
         mat(k,2399) = -rxt(k,394)*y(k,249)
         mat(k,2280) = -rxt(k,395)*y(k,249)
         mat(k,2221) = -(rxt(k,396) + rxt(k,397)) * y(k,249)
         mat(k,1313) = .500_r8*rxt(k,399)*y(k,263)
         mat(k,385) = .200_r8*rxt(k,400)*y(k,263)
         mat(k,1417) = rxt(k,413)*y(k,263)
         mat(k,1889) = .500_r8*rxt(k,399)*y(k,111) + .200_r8*rxt(k,400)*y(k,112) &
                      + rxt(k,413)*y(k,118)
         mat(k,782) = -(rxt(k,473)*y(k,245) + rxt(k,474)*y(k,131) + rxt(k,475) &
                      *y(k,132))
         mat(k,2366) = -rxt(k,473)*y(k,250)
         mat(k,2185) = -rxt(k,474)*y(k,250)
         mat(k,2461) = -rxt(k,475)*y(k,250)
         mat(k,1439) = -(rxt(k,401)*y(k,238) + rxt(k,402)*y(k,239) + rxt(k,403) &
                      *y(k,245) + 4._r8*rxt(k,404)*y(k,251) + rxt(k,405)*y(k,131) &
                      + rxt(k,406)*y(k,133) + rxt(k,414)*y(k,132))
         mat(k,1469) = -rxt(k,401)*y(k,251)
         mat(k,1678) = -rxt(k,402)*y(k,251)
         mat(k,2401) = -rxt(k,403)*y(k,251)
         mat(k,2223) = -rxt(k,405)*y(k,251)
         mat(k,2282) = -rxt(k,406)*y(k,251)
         mat(k,2472) = -rxt(k,414)*y(k,251)
         mat(k,1314) = .500_r8*rxt(k,399)*y(k,263)
         mat(k,386) = .500_r8*rxt(k,400)*y(k,263)
         mat(k,1891) = .500_r8*rxt(k,399)*y(k,111) + .500_r8*rxt(k,400)*y(k,112)
         mat(k,967) = -(rxt(k,476)*y(k,245) + rxt(k,477)*y(k,131) + rxt(k,478) &
                      *y(k,132))
         mat(k,2378) = -rxt(k,476)*y(k,252)
         mat(k,2197) = -rxt(k,477)*y(k,252)
         mat(k,2464) = -rxt(k,478)*y(k,252)
         mat(k,719) = -(rxt(k,407)*y(k,245) + rxt(k,408)*y(k,131))
         mat(k,2361) = -rxt(k,407)*y(k,253)
         mat(k,2184) = -rxt(k,408)*y(k,253)
         mat(k,562) = rxt(k,409)*y(k,263)
         mat(k,395) = rxt(k,410)*y(k,263)
         mat(k,1840) = rxt(k,409)*y(k,113) + rxt(k,410)*y(k,114)
         mat(k,573) = -(rxt(k,208)*y(k,140) + rxt(k,209)*y(k,141))
         mat(k,1920) = -rxt(k,208)*y(k,254)
         mat(k,1570) = -rxt(k,209)*y(k,254)
         mat(k,1920) = mat(k,1920) + rxt(k,607)*y(k,255)
         mat(k,902) = .900_r8*rxt(k,605)*y(k,255) + .800_r8*rxt(k,603)*y(k,256)
         mat(k,736) = rxt(k,607)*y(k,140) + .900_r8*rxt(k,605)*y(k,241)
         mat(k,886) = .800_r8*rxt(k,603)*y(k,241)
         mat(k,737) = -(rxt(k,605)*y(k,241) + rxt(k,606)*y(k,141) + (rxt(k,607) &
                      + rxt(k,608)) * y(k,140))
         mat(k,903) = -rxt(k,605)*y(k,255)
         mat(k,1571) = -rxt(k,606)*y(k,255)
         mat(k,1923) = -(rxt(k,607) + rxt(k,608)) * y(k,255)
         mat(k,887) = -(rxt(k,603)*y(k,241))
         mat(k,905) = -rxt(k,603)*y(k,256)
         mat(k,998) = rxt(k,612)*y(k,262)
         mat(k,2191) = rxt(k,614)*y(k,262)
         mat(k,1929) = rxt(k,607)*y(k,255)
         mat(k,1574) = rxt(k,611)*y(k,257)
         mat(k,739) = rxt(k,607)*y(k,140)
         mat(k,548) = rxt(k,611)*y(k,141)
         mat(k,894) = rxt(k,612)*y(k,119) + rxt(k,614)*y(k,131)
         mat(k,546) = -(rxt(k,609)*y(k,140) + (rxt(k,610) + rxt(k,611)) * y(k,141))
         mat(k,1919) = -rxt(k,609)*y(k,257)
         mat(k,1569) = -(rxt(k,610) + rxt(k,611)) * y(k,257)
         mat(k,1132) = -(rxt(k,505)*y(k,239) + rxt(k,506)*y(k,245) + rxt(k,507) &
                      *y(k,131) + rxt(k,508)*y(k,133))
         mat(k,1661) = -rxt(k,505)*y(k,258)
         mat(k,2385) = -rxt(k,506)*y(k,258)
         mat(k,2205) = -rxt(k,507)*y(k,258)
         mat(k,2263) = -rxt(k,508)*y(k,258)
         mat(k,1054) = rxt(k,499)*y(k,133)
         mat(k,1026) = rxt(k,502)*y(k,133)
         mat(k,2263) = mat(k,2263) + rxt(k,499)*y(k,6) + rxt(k,502)*y(k,116) &
                      + .500_r8*rxt(k,519)*y(k,218)
         mat(k,455) = rxt(k,509)*y(k,263)
         mat(k,1103) = .500_r8*rxt(k,519)*y(k,133)
         mat(k,1872) = rxt(k,509)*y(k,135)
         mat(k,2039) = -(rxt(k,173)*y(k,79) + rxt(k,174)*y(k,274) + (rxt(k,176) &
                      + rxt(k,177)) * y(k,141) + rxt(k,178)*y(k,142) + (rxt(k,226) &
                      + rxt(k,227)) * y(k,120) + rxt(k,259)*y(k,35) + rxt(k,260) &
                      *y(k,36) + rxt(k,261)*y(k,38) + rxt(k,262)*y(k,39) + rxt(k,263) &
                      *y(k,40) + rxt(k,264)*y(k,41) + rxt(k,265)*y(k,42) + (rxt(k,266) &
                      + rxt(k,267)) * y(k,87) + rxt(k,286)*y(k,37) + rxt(k,287) &
                      *y(k,57) + rxt(k,288)*y(k,80) + (rxt(k,289) + rxt(k,290) &
                      ) * y(k,83) + rxt(k,295)*y(k,66) + rxt(k,296)*y(k,67) + rxt(k,309) &
                      *y(k,43) + rxt(k,310)*y(k,45) + rxt(k,311)*y(k,84) + rxt(k,312) &
                      *y(k,85) + rxt(k,313)*y(k,86) + (rxt(k,330) + rxt(k,331) &
                      + rxt(k,332)) * y(k,56) + rxt(k,333)*y(k,88))
         mat(k,1512) = -rxt(k,173)*y(k,259)
         mat(k,2514) = -rxt(k,174)*y(k,259)
         mat(k,1590) = -(rxt(k,176) + rxt(k,177)) * y(k,259)
         mat(k,2103) = -rxt(k,178)*y(k,259)
         mat(k,302) = -(rxt(k,226) + rxt(k,227)) * y(k,259)
         mat(k,144) = -rxt(k,259)*y(k,259)
         mat(k,181) = -rxt(k,260)*y(k,259)
         mat(k,158) = -rxt(k,261)*y(k,259)
         mat(k,191) = -rxt(k,262)*y(k,259)
         mat(k,162) = -rxt(k,263)*y(k,259)
         mat(k,196) = -rxt(k,264)*y(k,259)
         mat(k,166) = -rxt(k,265)*y(k,259)
         mat(k,1545) = -(rxt(k,266) + rxt(k,267)) * y(k,259)
         mat(k,186) = -rxt(k,286)*y(k,259)
         mat(k,452) = -rxt(k,287)*y(k,259)
         mat(k,154) = -rxt(k,288)*y(k,259)
         mat(k,876) = -(rxt(k,289) + rxt(k,290)) * y(k,259)
         mat(k,284) = -rxt(k,295)*y(k,259)
         mat(k,293) = -rxt(k,296)*y(k,259)
         mat(k,524) = -rxt(k,309)*y(k,259)
         mat(k,632) = -rxt(k,310)*y(k,259)
         mat(k,289) = -rxt(k,311)*y(k,259)
         mat(k,299) = -rxt(k,312)*y(k,259)
         mat(k,358) = -rxt(k,313)*y(k,259)
         mat(k,1638) = -(rxt(k,330) + rxt(k,331) + rxt(k,332)) * y(k,259)
         mat(k,228) = -rxt(k,333)*y(k,259)
         mat(k,1590) = mat(k,1590) + rxt(k,209)*y(k,254)
         mat(k,913) = .850_r8*rxt(k,604)*y(k,262)
         mat(k,576) = rxt(k,209)*y(k,141)
         mat(k,900) = .850_r8*rxt(k,604)*y(k,241)
         mat(k,221) = -(rxt(k,180)*y(k,140) + rxt(k,181)*y(k,141))
         mat(k,1916) = -rxt(k,180)*y(k,260)
         mat(k,1566) = -rxt(k,181)*y(k,260)
         mat(k,1484) = rxt(k,182)*y(k,261)
         mat(k,1916) = mat(k,1916) + rxt(k,184)*y(k,261)
         mat(k,1566) = mat(k,1566) + rxt(k,185)*y(k,261)
         mat(k,2053) = rxt(k,186)*y(k,261)
         mat(k,223) = rxt(k,182)*y(k,65) + rxt(k,184)*y(k,140) + rxt(k,185)*y(k,141) &
                      + rxt(k,186)*y(k,142)
         mat(k,224) = -(rxt(k,182)*y(k,65) + rxt(k,184)*y(k,140) + rxt(k,185)*y(k,141) &
                      + rxt(k,186)*y(k,142))
         mat(k,1485) = -rxt(k,182)*y(k,261)
         mat(k,1917) = -rxt(k,184)*y(k,261)
         mat(k,1567) = -rxt(k,185)*y(k,261)
         mat(k,2054) = -rxt(k,186)*y(k,261)
         mat(k,1567) = mat(k,1567) + rxt(k,176)*y(k,259)
         mat(k,2014) = rxt(k,176)*y(k,141)
         mat(k,895) = -(rxt(k,604)*y(k,241) + rxt(k,612)*y(k,119) + rxt(k,614) &
                      *y(k,131))
         mat(k,906) = -rxt(k,604)*y(k,262)
         mat(k,999) = -rxt(k,612)*y(k,262)
         mat(k,2192) = -rxt(k,614)*y(k,262)
         mat(k,1488) = rxt(k,615)*y(k,264)
         mat(k,1575) = rxt(k,606)*y(k,255) + rxt(k,610)*y(k,257) + rxt(k,617)*y(k,264)
         mat(k,740) = rxt(k,606)*y(k,141)
         mat(k,549) = rxt(k,610)*y(k,141)
         mat(k,849) = rxt(k,615)*y(k,65) + rxt(k,617)*y(k,141)
         mat(k,1904) = -(rxt(k,199)*y(k,79) + rxt(k,200)*y(k,81) + rxt(k,201)*y(k,245) &
                      + rxt(k,202)*y(k,140) + rxt(k,203)*y(k,142) + (4._r8*rxt(k,204) &
                      + 4._r8*rxt(k,205)) * y(k,263) + rxt(k,207)*y(k,92) + rxt(k,221) &
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
                      *y(k,147) + rxt(k,363)*y(k,27) + rxt(k,371)*y(k,32) + rxt(k,372) &
                      *y(k,33) + rxt(k,374)*y(k,51) + rxt(k,376)*y(k,97) + rxt(k,377) &
                      *y(k,134) + rxt(k,380)*y(k,156) + rxt(k,384)*y(k,157) + rxt(k,385) &
                      *y(k,31) + rxt(k,386)*y(k,50) + rxt(k,388)*y(k,18) + rxt(k,391) &
                      *y(k,95) + rxt(k,399)*y(k,111) + rxt(k,400)*y(k,112) + rxt(k,409) &
                      *y(k,113) + rxt(k,410)*y(k,114) + rxt(k,411)*y(k,115) + rxt(k,413) &
                      *y(k,118) + rxt(k,416)*y(k,1) + rxt(k,420)*y(k,2) + rxt(k,421) &
                      *y(k,17) + rxt(k,422)*y(k,96) + rxt(k,423)*y(k,98) + rxt(k,424) &
                      *y(k,99) + rxt(k,436)*y(k,101) + rxt(k,437)*y(k,102) + rxt(k,444) &
                      *y(k,104) + rxt(k,446)*y(k,100) + rxt(k,447)*y(k,106) + rxt(k,448) &
                      *y(k,122) + rxt(k,449)*y(k,123) + rxt(k,455)*y(k,223) + rxt(k,458) &
                      *y(k,8) + rxt(k,461)*y(k,10) + rxt(k,462)*y(k,24) + rxt(k,464) &
                      *y(k,25) + rxt(k,468)*y(k,34) + rxt(k,469)*y(k,68) + rxt(k,481) &
                      *y(k,150) + rxt(k,484)*y(k,151) + rxt(k,488)*y(k,220) + (rxt(k,489) &
                      + rxt(k,579)) * y(k,221) + rxt(k,491)*y(k,224) + rxt(k,494) &
                      *y(k,225) + rxt(k,497)*y(k,227) + rxt(k,498)*y(k,228) + rxt(k,501) &
                      *y(k,6) + rxt(k,504)*y(k,116) + rxt(k,509)*y(k,135) + rxt(k,513) &
                      *y(k,215) + rxt(k,514)*y(k,216) + rxt(k,518)*y(k,217) + rxt(k,520) &
                      *y(k,218) + rxt(k,521)*y(k,219) + (rxt(k,523) + rxt(k,537) &
                      ) * y(k,69) + rxt(k,525)*y(k,145) + rxt(k,527)*y(k,161) &
                      + rxt(k,531)*y(k,158) + rxt(k,536)*y(k,160) + rxt(k,539) &
                      *y(k,127))
         mat(k,1509) = -rxt(k,199)*y(k,263)
         mat(k,673) = -rxt(k,200)*y(k,263)
         mat(k,2413) = -rxt(k,201)*y(k,263)
         mat(k,1947) = -rxt(k,202)*y(k,263)
         mat(k,2100) = -rxt(k,203)*y(k,263)
         mat(k,514) = -rxt(k,207)*y(k,263)
         mat(k,2293) = -rxt(k,221)*y(k,263)
         mat(k,1006) = -rxt(k,222)*y(k,263)
         mat(k,2484) = -rxt(k,230)*y(k,263)
         mat(k,1710) = -rxt(k,231)*y(k,263)
         mat(k,1121) = -rxt(k,250)*y(k,263)
         mat(k,2127) = -(rxt(k,252) + rxt(k,253)) * y(k,263)
         mat(k,1542) = -rxt(k,255)*y(k,263)
         mat(k,881) = -rxt(k,258)*y(k,263)
         mat(k,1610) = -rxt(k,282)*y(k,263)
         mat(k,874) = -rxt(k,284)*y(k,263)
         mat(k,522) = -rxt(k,298)*y(k,263)
         mat(k,630) = -rxt(k,300)*y(k,263)
         mat(k,168) = -rxt(k,301)*y(k,263)
         mat(k,403) = -rxt(k,303)*y(k,263)
         mat(k,450) = -rxt(k,305)*y(k,263)
         mat(k,287) = -rxt(k,306)*y(k,263)
         mat(k,297) = -rxt(k,307)*y(k,263)
         mat(k,356) = -rxt(k,308)*y(k,263)
         mat(k,2439) = -rxt(k,317)*y(k,263)
         mat(k,856) = -rxt(k,322)*y(k,263)
         mat(k,444) = -rxt(k,323)*y(k,263)
         mat(k,1635) = -rxt(k,324)*y(k,263)
         mat(k,227) = -rxt(k,325)*y(k,263)
         mat(k,964) = -rxt(k,326)*y(k,263)
         mat(k,1204) = -rxt(k,334)*y(k,263)
         mat(k,321) = -rxt(k,336)*y(k,263)
         mat(k,311) = -rxt(k,343)*y(k,263)
         mat(k,379) = -rxt(k,344)*y(k,263)
         mat(k,326) = -rxt(k,346)*y(k,263)
         mat(k,1198) = -rxt(k,348)*y(k,263)
         mat(k,146) = -rxt(k,349)*y(k,263)
         mat(k,750) = -rxt(k,354)*y(k,263)
         mat(k,639) = -rxt(k,355)*y(k,263)
         mat(k,1210) = -rxt(k,360)*y(k,263)
         mat(k,1087) = -rxt(k,361)*y(k,263)
         mat(k,606) = -rxt(k,362)*y(k,263)
         mat(k,582) = -rxt(k,363)*y(k,263)
         mat(k,438) = -rxt(k,371)*y(k,263)
         mat(k,332) = -rxt(k,372)*y(k,263)
         mat(k,1327) = -rxt(k,374)*y(k,263)
         mat(k,1253) = -rxt(k,376)*y(k,263)
         mat(k,918) = -rxt(k,377)*y(k,263)
         mat(k,598) = -rxt(k,380)*y(k,263)
         mat(k,432) = -rxt(k,384)*y(k,263)
         mat(k,1185) = -rxt(k,385)*y(k,263)
         mat(k,1114) = -rxt(k,386)*y(k,263)
         mat(k,411) = -rxt(k,388)*y(k,263)
         mat(k,1244) = -rxt(k,391)*y(k,263)
         mat(k,1317) = -rxt(k,399)*y(k,263)
         mat(k,387) = -rxt(k,400)*y(k,263)
         mat(k,565) = -rxt(k,409)*y(k,263)
         mat(k,398) = -rxt(k,410)*y(k,263)
         mat(k,666) = -rxt(k,411)*y(k,263)
         mat(k,1426) = -rxt(k,413)*y(k,263)
         mat(k,714) = -rxt(k,416)*y(k,263)
         mat(k,704) = -rxt(k,420)*y(k,263)
         mat(k,267) = -rxt(k,421)*y(k,263)
         mat(k,280) = -rxt(k,422)*y(k,263)
         mat(k,383) = -rxt(k,423)*y(k,263)
         mat(k,173) = -rxt(k,424)*y(k,263)
         mat(k,647) = -rxt(k,436)*y(k,263)
         mat(k,616) = -rxt(k,437)*y(k,263)
         mat(k,474) = -rxt(k,444)*y(k,263)
         mat(k,935) = -rxt(k,446)*y(k,263)
         mat(k,757) = -rxt(k,447)*y(k,263)
         mat(k,463) = -rxt(k,448)*y(k,263)
         mat(k,1166) = -rxt(k,449)*y(k,263)
         mat(k,249) = -rxt(k,455)*y(k,263)
         mat(k,208) = -rxt(k,458)*y(k,263)
         mat(k,469) = -rxt(k,461)*y(k,263)
         mat(k,276) = -rxt(k,462)*y(k,263)
         mat(k,374) = -rxt(k,464)*y(k,263)
         mat(k,316) = -rxt(k,468)*y(k,263)
         mat(k,241) = -rxt(k,469)*y(k,263)
         mat(k,217) = -rxt(k,481)*y(k,263)
         mat(k,368) = -rxt(k,484)*y(k,263)
         mat(k,660) = -rxt(k,488)*y(k,263)
         mat(k,236) = -(rxt(k,489) + rxt(k,579)) * y(k,263)
         mat(k,258) = -rxt(k,491)*y(k,263)
         mat(k,773) = -rxt(k,494)*y(k,263)
         mat(k,263) = -rxt(k,497)*y(k,263)
         mat(k,481) = -rxt(k,498)*y(k,263)
         mat(k,1061) = -rxt(k,501)*y(k,263)
         mat(k,1033) = -rxt(k,504)*y(k,263)
         mat(k,457) = -rxt(k,509)*y(k,263)
         mat(k,733) = -rxt(k,513)*y(k,263)
         mat(k,685) = -rxt(k,514)*y(k,263)
         mat(k,534) = -rxt(k,518)*y(k,263)
         mat(k,1107) = -rxt(k,520)*y(k,263)
         mat(k,1151) = -rxt(k,521)*y(k,263)
         mat(k,350) = -(rxt(k,523) + rxt(k,537)) * y(k,263)
         mat(k,421) = -rxt(k,525)*y(k,263)
         mat(k,991) = -rxt(k,527)*y(k,263)
         mat(k,779) = -rxt(k,531)*y(k,263)
         mat(k,1523) = -rxt(k,536)*y(k,263)
         mat(k,140) = -rxt(k,539)*y(k,263)
         mat(k,1061) = mat(k,1061) + .630_r8*rxt(k,500)*y(k,142)
         mat(k,321) = mat(k,321) + .650_r8*rxt(k,336)*y(k,263)
         mat(k,582) = mat(k,582) + .130_r8*rxt(k,338)*y(k,142)
         mat(k,379) = mat(k,379) + .500_r8*rxt(k,344)*y(k,263)
         mat(k,1185) = mat(k,1185) + .360_r8*rxt(k,367)*y(k,142)
         mat(k,2439) = mat(k,2439) + rxt(k,316)*y(k,140)
         mat(k,444) = mat(k,444) + .300_r8*rxt(k,323)*y(k,263)
         mat(k,1635) = mat(k,1635) + rxt(k,330)*y(k,259)
         mat(k,1993) = rxt(k,239)*y(k,245)
         mat(k,947) = rxt(k,293)*y(k,274)
         mat(k,1731) = rxt(k,198)*y(k,142) + 2.000_r8*rxt(k,193)*y(k,245)
         mat(k,1509) = mat(k,1509) + rxt(k,190)*y(k,140) + rxt(k,173)*y(k,259)
         mat(k,673) = mat(k,673) + rxt(k,191)*y(k,140)
         mat(k,874) = mat(k,874) + rxt(k,283)*y(k,140) + rxt(k,289)*y(k,259)
         mat(k,1542) = mat(k,1542) + rxt(k,254)*y(k,140) + rxt(k,266)*y(k,259)
         mat(k,227) = mat(k,227) + rxt(k,333)*y(k,259)
         mat(k,833) = rxt(k,285)*y(k,140)
         mat(k,881) = mat(k,881) + rxt(k,257)*y(k,140)
         mat(k,935) = mat(k,935) + .320_r8*rxt(k,445)*y(k,142)
         mat(k,757) = mat(k,757) + .600_r8*rxt(k,447)*y(k,263)
         mat(k,1317) = mat(k,1317) + .240_r8*rxt(k,398)*y(k,142)
         mat(k,387) = mat(k,387) + .100_r8*rxt(k,400)*y(k,263)
         mat(k,1033) = mat(k,1033) + .630_r8*rxt(k,503)*y(k,142)
         mat(k,1426) = mat(k,1426) + .360_r8*rxt(k,412)*y(k,142)
         mat(k,2233) = rxt(k,223)*y(k,245)
         mat(k,2293) = mat(k,2293) + rxt(k,218)*y(k,245)
         mat(k,1947) = mat(k,1947) + rxt(k,316)*y(k,44) + rxt(k,190)*y(k,79) &
                      + rxt(k,191)*y(k,81) + rxt(k,283)*y(k,83) + rxt(k,254)*y(k,87) &
                      + rxt(k,285)*y(k,93) + rxt(k,257)*y(k,94) + rxt(k,196)*y(k,245)
         mat(k,2100) = mat(k,2100) + .630_r8*rxt(k,500)*y(k,6) + .130_r8*rxt(k,338) &
                      *y(k,27) + .360_r8*rxt(k,367)*y(k,31) + rxt(k,198)*y(k,78) &
                      + .320_r8*rxt(k,445)*y(k,100) + .240_r8*rxt(k,398)*y(k,111) &
                      + .630_r8*rxt(k,503)*y(k,116) + .360_r8*rxt(k,412)*y(k,118) &
                      + rxt(k,197)*y(k,245)
         mat(k,598) = mat(k,598) + .500_r8*rxt(k,380)*y(k,263)
         mat(k,249) = mat(k,249) + .500_r8*rxt(k,455)*y(k,263)
         mat(k,569) = .400_r8*rxt(k,456)*y(k,245)
         mat(k,1476) = .450_r8*rxt(k,352)*y(k,245)
         mat(k,824) = .400_r8*rxt(k,470)*y(k,245)
         mat(k,2413) = mat(k,2413) + rxt(k,239)*y(k,58) + 2.000_r8*rxt(k,193)*y(k,78) &
                      + rxt(k,223)*y(k,131) + rxt(k,218)*y(k,133) + rxt(k,196) &
                      *y(k,140) + rxt(k,197)*y(k,142) + .400_r8*rxt(k,456)*y(k,231) &
                      + .450_r8*rxt(k,352)*y(k,238) + .400_r8*rxt(k,470)*y(k,240) &
                      + .450_r8*rxt(k,403)*y(k,251) + .400_r8*rxt(k,476)*y(k,252) &
                      + .200_r8*rxt(k,407)*y(k,253) + .150_r8*rxt(k,382)*y(k,267)
         mat(k,1444) = .450_r8*rxt(k,403)*y(k,245)
         mat(k,972) = .400_r8*rxt(k,476)*y(k,245)
         mat(k,722) = .200_r8*rxt(k,407)*y(k,245)
         mat(k,2036) = rxt(k,330)*y(k,56) + rxt(k,173)*y(k,79) + rxt(k,289)*y(k,83) &
                      + rxt(k,266)*y(k,87) + rxt(k,333)*y(k,88) + 2.000_r8*rxt(k,174) &
                      *y(k,274)
         mat(k,1904) = mat(k,1904) + .650_r8*rxt(k,336)*y(k,26) + .500_r8*rxt(k,344) &
                      *y(k,29) + .300_r8*rxt(k,323)*y(k,55) + .600_r8*rxt(k,447) &
                      *y(k,106) + .100_r8*rxt(k,400)*y(k,112) + .500_r8*rxt(k,380) &
                      *y(k,156) + .500_r8*rxt(k,455)*y(k,223)
         mat(k,1263) = .150_r8*rxt(k,382)*y(k,245)
         mat(k,2511) = rxt(k,293)*y(k,75) + 2.000_r8*rxt(k,174)*y(k,259)
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
         mat(k,847) = -(rxt(k,615)*y(k,65) + rxt(k,617)*y(k,141))
         mat(k,1486) = -rxt(k,615)*y(k,264)
         mat(k,1573) = -rxt(k,617)*y(k,264)
         mat(k,1926) = rxt(k,608)*y(k,255) + rxt(k,609)*y(k,257)
         mat(k,738) = rxt(k,608)*y(k,140)
         mat(k,547) = rxt(k,609)*y(k,140)
         mat(k,498) = -(rxt(k,479)*y(k,245) + rxt(k,480)*y(k,131))
         mat(k,2345) = -rxt(k,479)*y(k,265)
         mat(k,2169) = -rxt(k,480)*y(k,265)
         mat(k,239) = .200_r8*rxt(k,469)*y(k,263)
         mat(k,215) = .140_r8*rxt(k,481)*y(k,263)
         mat(k,366) = rxt(k,484)*y(k,263)
         mat(k,1814) = .200_r8*rxt(k,469)*y(k,68) + .140_r8*rxt(k,481)*y(k,150) &
                      + rxt(k,484)*y(k,151)
         mat(k,860) = -(rxt(k,378)*y(k,245) + rxt(k,379)*y(k,131))
         mat(k,2372) = -rxt(k,378)*y(k,266)
         mat(k,2190) = -rxt(k,379)*y(k,266)
         mat(k,1172) = rxt(k,385)*y(k,263)
         mat(k,595) = .500_r8*rxt(k,380)*y(k,263)
         mat(k,1852) = rxt(k,385)*y(k,31) + .500_r8*rxt(k,380)*y(k,156)
         mat(k,1259) = -(rxt(k,381)*y(k,239) + rxt(k,382)*y(k,245) + rxt(k,383) &
                      *y(k,131))
         mat(k,1669) = -rxt(k,381)*y(k,267)
         mat(k,2392) = -rxt(k,382)*y(k,267)
         mat(k,2214) = -rxt(k,383)*y(k,267)
         mat(k,1057) = .060_r8*rxt(k,500)*y(k,142)
         mat(k,1111) = rxt(k,386)*y(k,263)
         mat(k,1029) = .060_r8*rxt(k,503)*y(k,142)
         mat(k,2081) = .060_r8*rxt(k,500)*y(k,6) + .060_r8*rxt(k,503)*y(k,116)
         mat(k,430) = rxt(k,384)*y(k,263)
         mat(k,1147) = .150_r8*rxt(k,521)*y(k,263)
         mat(k,1882) = rxt(k,386)*y(k,50) + rxt(k,384)*y(k,157) + .150_r8*rxt(k,521) &
                      *y(k,219)
         mat(k,1220) = -(rxt(k,510)*y(k,239) + rxt(k,511)*y(k,245) + rxt(k,512) &
                      *y(k,131))
         mat(k,1667) = -rxt(k,510)*y(k,268)
         mat(k,2390) = -rxt(k,511)*y(k,268)
         mat(k,2211) = -rxt(k,512)*y(k,268)
         mat(k,2270) = .500_r8*rxt(k,519)*y(k,218)
         mat(k,731) = rxt(k,513)*y(k,263)
         mat(k,1106) = .500_r8*rxt(k,519)*y(k,133) + rxt(k,520)*y(k,263)
         mat(k,1879) = rxt(k,513)*y(k,215) + rxt(k,520)*y(k,218)
         mat(k,1092) = -(rxt(k,515)*y(k,239) + rxt(k,516)*y(k,245) + rxt(k,517) &
                      *y(k,131))
         mat(k,1658) = -rxt(k,515)*y(k,269)
         mat(k,2382) = -rxt(k,516)*y(k,269)
         mat(k,2202) = -rxt(k,517)*y(k,269)
         mat(k,1051) = rxt(k,501)*y(k,263)
         mat(k,1023) = rxt(k,504)*y(k,263)
         mat(k,531) = rxt(k,518)*y(k,263)
         mat(k,1868) = rxt(k,501)*y(k,6) + rxt(k,504)*y(k,116) + rxt(k,518)*y(k,217)
         mat(k,793) = -(rxt(k,486)*y(k,245) + rxt(k,487)*y(k,131))
         mat(k,2367) = -rxt(k,486)*y(k,270)
         mat(k,2186) = -rxt(k,487)*y(k,270)
         mat(k,656) = rxt(k,488)*y(k,263)
         mat(k,235) = (.650_r8*rxt(k,489)+rxt(k,579))*y(k,263)
         mat(k,1847) = rxt(k,488)*y(k,220) + (.650_r8*rxt(k,489)+rxt(k,579))*y(k,221)
         mat(k,1275) = -(rxt(k,450)*y(k,238) + rxt(k,451)*y(k,239) + rxt(k,452) &
                      *y(k,245) + rxt(k,453)*y(k,131) + rxt(k,454)*y(k,133))
         mat(k,1461) = -rxt(k,450)*y(k,271)
         mat(k,1670) = -rxt(k,451)*y(k,271)
         mat(k,2393) = -rxt(k,452)*y(k,271)
         mat(k,2215) = -rxt(k,453)*y(k,271)
         mat(k,2274) = -rxt(k,454)*y(k,271)
         mat(k,279) = rxt(k,422)*y(k,263)
         mat(k,382) = rxt(k,423)*y(k,263)
         mat(k,172) = rxt(k,424)*y(k,263)
         mat(k,754) = .400_r8*rxt(k,447)*y(k,263)
         mat(k,248) = .500_r8*rxt(k,455)*y(k,263)
         mat(k,1883) = rxt(k,422)*y(k,96) + rxt(k,423)*y(k,98) + rxt(k,424)*y(k,99) &
                      + .400_r8*rxt(k,447)*y(k,106) + .500_r8*rxt(k,455)*y(k,223)
         mat(k,809) = -(rxt(k,492)*y(k,245) + rxt(k,493)*y(k,131))
         mat(k,2368) = -rxt(k,492)*y(k,272)
         mat(k,2187) = -rxt(k,493)*y(k,272)
         mat(k,255) = .560_r8*rxt(k,491)*y(k,263)
         mat(k,766) = rxt(k,494)*y(k,263)
         mat(k,1848) = .560_r8*rxt(k,491)*y(k,224) + rxt(k,494)*y(k,225)
         mat(k,554) = -(rxt(k,495)*y(k,245) + rxt(k,496)*y(k,131))
         mat(k,2352) = -rxt(k,495)*y(k,273)
         mat(k,2174) = -rxt(k,496)*y(k,273)
         mat(k,262) = .300_r8*rxt(k,497)*y(k,263)
         mat(k,478) = rxt(k,498)*y(k,263)
         mat(k,1821) = .300_r8*rxt(k,497)*y(k,227) + rxt(k,498)*y(k,228)
         mat(k,2522) = -(rxt(k,174)*y(k,259) + rxt(k,293)*y(k,75) + rxt(k,538) &
                      *y(k,162))
         mat(k,2047) = -rxt(k,174)*y(k,274)
         mat(k,949) = -rxt(k,293)*y(k,274)
         mat(k,308) = -rxt(k,538)*y(k,274)
         mat(k,328) = rxt(k,346)*y(k,263)
         mat(k,440) = rxt(k,371)*y(k,263)
         mat(k,334) = rxt(k,372)*y(k,263)
         mat(k,525) = rxt(k,298)*y(k,263)
         mat(k,2450) = rxt(k,317)*y(k,263)
         mat(k,634) = rxt(k,300)*y(k,263)
         mat(k,170) = rxt(k,301)*y(k,263)
         mat(k,1201) = rxt(k,348)*y(k,263)
         mat(k,406) = rxt(k,303)*y(k,263)
         mat(k,1115) = rxt(k,386)*y(k,263)
         mat(k,1330) = rxt(k,374)*y(k,263)
         mat(k,751) = rxt(k,354)*y(k,263)
         mat(k,641) = rxt(k,355)*y(k,263)
         mat(k,446) = rxt(k,323)*y(k,263)
         mat(k,1646) = rxt(k,324)*y(k,263)
         mat(k,1742) = rxt(k,194)*y(k,245)
         mat(k,1514) = rxt(k,199)*y(k,263)
         mat(k,677) = rxt(k,200)*y(k,263)
         mat(k,877) = rxt(k,284)*y(k,263)
         mat(k,359) = rxt(k,308)*y(k,263)
         mat(k,1549) = (rxt(k,594)+rxt(k,599))*y(k,93) + (rxt(k,587)+rxt(k,593) &
                       +rxt(k,598))*y(k,94) + rxt(k,255)*y(k,263)
         mat(k,966) = rxt(k,326)*y(k,263)
         mat(k,1721) = rxt(k,231)*y(k,263)
         mat(k,518) = rxt(k,207)*y(k,263)
         mat(k,836) = (rxt(k,594)+rxt(k,599))*y(k,87)
         mat(k,885) = (rxt(k,587)+rxt(k,593)+rxt(k,598))*y(k,87) + rxt(k,258)*y(k,263)
         mat(k,1321) = .500_r8*rxt(k,399)*y(k,263)
         mat(k,141) = rxt(k,539)*y(k,263)
         mat(k,601) = rxt(k,380)*y(k,263)
         mat(k,434) = rxt(k,384)*y(k,263)
         mat(k,2424) = rxt(k,194)*y(k,78) + rxt(k,201)*y(k,263)
         mat(k,1915) = rxt(k,346)*y(k,30) + rxt(k,371)*y(k,32) + rxt(k,372)*y(k,33) &
                      + rxt(k,298)*y(k,43) + rxt(k,317)*y(k,44) + rxt(k,300)*y(k,45) &
                      + rxt(k,301)*y(k,46) + rxt(k,348)*y(k,47) + rxt(k,303)*y(k,48) &
                      + rxt(k,386)*y(k,50) + rxt(k,374)*y(k,51) + rxt(k,354)*y(k,52) &
                      + rxt(k,355)*y(k,53) + rxt(k,323)*y(k,55) + rxt(k,324)*y(k,56) &
                      + rxt(k,199)*y(k,79) + rxt(k,200)*y(k,81) + rxt(k,284)*y(k,83) &
                      + rxt(k,308)*y(k,86) + rxt(k,255)*y(k,87) + rxt(k,326)*y(k,89) &
                      + rxt(k,231)*y(k,91) + rxt(k,207)*y(k,92) + rxt(k,258)*y(k,94) &
                      + .500_r8*rxt(k,399)*y(k,111) + rxt(k,539)*y(k,127) + rxt(k,380) &
                      *y(k,156) + rxt(k,384)*y(k,157) + rxt(k,201)*y(k,245) &
                      + 2.000_r8*rxt(k,204)*y(k,263)
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
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 82) = mat(k, 82) + lmat(k, 82)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = lmat(k, 99)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 114) = lmat(k, 114)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 148) = lmat(k, 148)
         mat(k, 149) = lmat(k, 149)
         mat(k, 150) = lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
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
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = lmat(k, 202)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 210) = lmat(k, 210)
         mat(k, 211) = lmat(k, 211)
         mat(k, 212) = lmat(k, 212)
         mat(k, 213) = lmat(k, 213)
         mat(k, 214) = mat(k, 214) + lmat(k, 214)
         mat(k, 219) = lmat(k, 219)
         mat(k, 220) = lmat(k, 220)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = lmat(k, 225)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 268) = lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = mat(k, 274) + lmat(k, 274)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 290) = mat(k, 290) + lmat(k, 290)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 295) = mat(k, 295) + lmat(k, 295)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 307) = lmat(k, 307)
         mat(k, 309) = mat(k, 309) + lmat(k, 309)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 314) = lmat(k, 314)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 323) = mat(k, 323) + lmat(k, 323)
         mat(k, 329) = mat(k, 329) + lmat(k, 329)
         mat(k, 335) = lmat(k, 335)
         mat(k, 336) = lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = lmat(k, 339)
         mat(k, 340) = lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = lmat(k, 345)
         mat(k, 346) = lmat(k, 346)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 353) = mat(k, 353) + lmat(k, 353)
         mat(k, 354) = mat(k, 354) + lmat(k, 354)
         mat(k, 357) = mat(k, 357) + lmat(k, 357)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 363) = lmat(k, 363)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 367) = lmat(k, 367)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = lmat(k, 370)
         mat(k, 371) = mat(k, 371) + lmat(k, 371)
         mat(k, 374) = mat(k, 374) + lmat(k, 374)
         mat(k, 375) = lmat(k, 375)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 378) = mat(k, 378) + lmat(k, 378)
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = mat(k, 381) + lmat(k, 381)
         mat(k, 384) = mat(k, 384) + lmat(k, 384)
         mat(k, 389) = lmat(k, 389)
         mat(k, 390) = lmat(k, 390)
         mat(k, 391) = lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 393) = lmat(k, 393)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = lmat(k, 397)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 402) = lmat(k, 402)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 416) = lmat(k, 416)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = lmat(k, 424)
         mat(k, 425) = lmat(k, 425)
         mat(k, 426) = lmat(k, 426)
         mat(k, 427) = lmat(k, 427)
         mat(k, 428) = lmat(k, 428)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 431) = lmat(k, 431)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 433) = lmat(k, 433)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 437) = lmat(k, 437)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 439) = lmat(k, 439)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 443) = lmat(k, 443)
         mat(k, 444) = mat(k, 444) + lmat(k, 444)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 454) = lmat(k, 454)
         mat(k, 456) = lmat(k, 456)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 458) = lmat(k, 458)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 462) = lmat(k, 462)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 466) = lmat(k, 466)
         mat(k, 468) = lmat(k, 468)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 470) = lmat(k, 470)
         mat(k, 471) = mat(k, 471) + lmat(k, 471)
         mat(k, 472) = lmat(k, 472)
         mat(k, 475) = mat(k, 475) + lmat(k, 475)
         mat(k, 476) = lmat(k, 476)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 479) = lmat(k, 479)
         mat(k, 480) = lmat(k, 480)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 482) = lmat(k, 482)
         mat(k, 485) = mat(k, 485) + lmat(k, 485)
         mat(k, 491) = mat(k, 491) + lmat(k, 491)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 495) = lmat(k, 495)
         mat(k, 498) = mat(k, 498) + lmat(k, 498)
         mat(k, 504) = lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 511) = lmat(k, 511)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 515) = lmat(k, 515)
         mat(k, 516) = lmat(k, 516)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 521) = mat(k, 521) + lmat(k, 521)
         mat(k, 526) = mat(k, 526) + lmat(k, 526)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 532) = lmat(k, 532)
         mat(k, 533) = lmat(k, 533)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 535) = lmat(k, 535)
         mat(k, 536) = lmat(k, 536)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 554) = mat(k, 554) + lmat(k, 554)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 563) = lmat(k, 563)
         mat(k, 564) = lmat(k, 564)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 586) = mat(k, 586) + lmat(k, 586)
         mat(k, 587) = lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 590) = lmat(k, 590)
         mat(k, 592) = mat(k, 592) + lmat(k, 592)
         mat(k, 593) = lmat(k, 593)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 596) = lmat(k, 596)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 599) = lmat(k, 599)
         mat(k, 600) = lmat(k, 600)
         mat(k, 602) = mat(k, 602) + lmat(k, 602)
         mat(k, 603) = lmat(k, 603)
         mat(k, 604) = lmat(k, 604)
         mat(k, 605) = lmat(k, 605)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 609) = lmat(k, 609)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 615) = lmat(k, 615)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 626) = mat(k, 626) + lmat(k, 626)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 629) = lmat(k, 629)
         mat(k, 635) = mat(k, 635) + lmat(k, 635)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 638) = lmat(k, 638)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 646) = lmat(k, 646)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = lmat(k, 652)
         mat(k, 653) = mat(k, 653) + lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 658) = lmat(k, 658)
         mat(k, 659) = lmat(k, 659)
         mat(k, 660) = mat(k, 660) + lmat(k, 660)
         mat(k, 661) = lmat(k, 661)
         mat(k, 662) = mat(k, 662) + lmat(k, 662)
         mat(k, 664) = lmat(k, 664)
         mat(k, 670) = lmat(k, 670)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 673) = mat(k, 673) + lmat(k, 673)
         mat(k, 678) = lmat(k, 678)
         mat(k, 679) = lmat(k, 679)
         mat(k, 680) = lmat(k, 680)
         mat(k, 681) = lmat(k, 681)
         mat(k, 682) = mat(k, 682) + lmat(k, 682)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 684) = lmat(k, 684)
         mat(k, 686) = lmat(k, 686)
         mat(k, 687) = mat(k, 687) + lmat(k, 687)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 696) = lmat(k, 696)
         mat(k, 697) = mat(k, 697) + lmat(k, 697)
         mat(k, 701) = lmat(k, 701)
         mat(k, 702) = lmat(k, 702)
         mat(k, 704) = mat(k, 704) + lmat(k, 704)
         mat(k, 705) = lmat(k, 705)
         mat(k, 706) = lmat(k, 706)
         mat(k, 707) = lmat(k, 707)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 711) = mat(k, 711) + lmat(k, 711)
         mat(k, 712) = mat(k, 712) + lmat(k, 712)
         mat(k, 715) = lmat(k, 715)
         mat(k, 716) = mat(k, 716) + lmat(k, 716)
         mat(k, 717) = mat(k, 717) + lmat(k, 717)
         mat(k, 719) = mat(k, 719) + lmat(k, 719)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 727) = lmat(k, 727)
         mat(k, 728) = lmat(k, 728)
         mat(k, 729) = lmat(k, 729)
         mat(k, 730) = lmat(k, 730)
         mat(k, 732) = lmat(k, 732)
         mat(k, 733) = mat(k, 733) + lmat(k, 733)
         mat(k, 734) = lmat(k, 734)
         mat(k, 735) = lmat(k, 735)
         mat(k, 737) = mat(k, 737) + lmat(k, 737)
         mat(k, 747) = mat(k, 747) + lmat(k, 747)
         mat(k, 753) = mat(k, 753) + lmat(k, 753)
         mat(k, 755) = lmat(k, 755)
         mat(k, 756) = lmat(k, 756)
         mat(k, 757) = mat(k, 757) + lmat(k, 757)
         mat(k, 758) = lmat(k, 758)
         mat(k, 759) = lmat(k, 759)
         mat(k, 760) = lmat(k, 760)
         mat(k, 761) = lmat(k, 761)
         mat(k, 762) = lmat(k, 762)
         mat(k, 763) = lmat(k, 763)
         mat(k, 764) = mat(k, 764) + lmat(k, 764)
         mat(k, 769) = lmat(k, 769)
         mat(k, 771) = lmat(k, 771)
         mat(k, 773) = mat(k, 773) + lmat(k, 773)
         mat(k, 774) = lmat(k, 774)
         mat(k, 775) = mat(k, 775) + lmat(k, 775)
         mat(k, 782) = mat(k, 782) + lmat(k, 782)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 809) = mat(k, 809) + lmat(k, 809)
         mat(k, 820) = mat(k, 820) + lmat(k, 820)
         mat(k, 829) = mat(k, 829) + lmat(k, 829)
         mat(k, 831) = lmat(k, 831)
         mat(k, 833) = mat(k, 833) + lmat(k, 833)
         mat(k, 839) = mat(k, 839) + lmat(k, 839)
         mat(k, 847) = mat(k, 847) + lmat(k, 847)
         mat(k, 848) = lmat(k, 848)
         mat(k, 850) = lmat(k, 850)
         mat(k, 855) = mat(k, 855) + lmat(k, 855)
         mat(k, 860) = mat(k, 860) + lmat(k, 860)
         mat(k, 870) = mat(k, 870) + lmat(k, 870)
         mat(k, 871) = mat(k, 871) + lmat(k, 871)
         mat(k, 873) = mat(k, 873) + lmat(k, 873)
         mat(k, 879) = mat(k, 879) + lmat(k, 879)
         mat(k, 881) = mat(k, 881) + lmat(k, 881)
         mat(k, 883) = mat(k, 883) + lmat(k, 883)
         mat(k, 887) = mat(k, 887) + lmat(k, 887)
         mat(k, 894) = mat(k, 894) + lmat(k, 894)
         mat(k, 895) = mat(k, 895) + lmat(k, 895)
         mat(k, 901) = mat(k, 901) + lmat(k, 901)
         mat(k, 907) = mat(k, 907) + lmat(k, 907)
         mat(k, 915) = mat(k, 915) + lmat(k, 915)
         mat(k, 917) = lmat(k, 917)
         mat(k, 919) = lmat(k, 919)
         mat(k, 920) = mat(k, 920) + lmat(k, 920)
         mat(k, 924) = mat(k, 924) + lmat(k, 924)
         mat(k, 941) = mat(k, 941) + lmat(k, 941)
         mat(k, 953) = mat(k, 953) + lmat(k, 953)
         mat(k, 962) = mat(k, 962) + lmat(k, 962)
         mat(k, 967) = mat(k, 967) + lmat(k, 967)
         mat(k, 977) = mat(k, 977) + lmat(k, 977)
         mat(k, 989) = mat(k, 989) + lmat(k, 989)
         mat(k, 990) = lmat(k, 990)
         mat(k, 992) = lmat(k, 992)
         mat(k, 996) = lmat(k, 996)
         mat(k,1000) = lmat(k,1000)
         mat(k,1001) = mat(k,1001) + lmat(k,1001)
         mat(k,1020) = mat(k,1020) + lmat(k,1020)
         mat(k,1048) = mat(k,1048) + lmat(k,1048)
         mat(k,1072) = mat(k,1072) + lmat(k,1072)
         mat(k,1083) = lmat(k,1083)
         mat(k,1084) = mat(k,1084) + lmat(k,1084)
         mat(k,1085) = mat(k,1085) + lmat(k,1085)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1092) = mat(k,1092) + lmat(k,1092)
         mat(k,1102) = mat(k,1102) + lmat(k,1102)
         mat(k,1104) = lmat(k,1104)
         mat(k,1105) = lmat(k,1105)
         mat(k,1109) = lmat(k,1109)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1112) = lmat(k,1112)
         mat(k,1113) = lmat(k,1113)
         mat(k,1117) = mat(k,1117) + lmat(k,1117)
         mat(k,1118) = mat(k,1118) + lmat(k,1118)
         mat(k,1120) = mat(k,1120) + lmat(k,1120)
         mat(k,1123) = mat(k,1123) + lmat(k,1123)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1125) = mat(k,1125) + lmat(k,1125)
         mat(k,1126) = lmat(k,1126)
         mat(k,1132) = mat(k,1132) + lmat(k,1132)
         mat(k,1144) = mat(k,1144) + lmat(k,1144)
         mat(k,1145) = mat(k,1145) + lmat(k,1145)
         mat(k,1146) = mat(k,1146) + lmat(k,1146)
         mat(k,1147) = mat(k,1147) + lmat(k,1147)
         mat(k,1148) = mat(k,1148) + lmat(k,1148)
         mat(k,1149) = mat(k,1149) + lmat(k,1149)
         mat(k,1152) = mat(k,1152) + lmat(k,1152)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1155) = lmat(k,1155)
         mat(k,1159) = mat(k,1159) + lmat(k,1159)
         mat(k,1165) = lmat(k,1165)
         mat(k,1168) = mat(k,1168) + lmat(k,1168)
         mat(k,1170) = lmat(k,1170)
         mat(k,1175) = mat(k,1175) + lmat(k,1175)
         mat(k,1193) = mat(k,1193) + lmat(k,1193)
         mat(k,1194) = lmat(k,1194)
         mat(k,1196) = lmat(k,1196)
         mat(k,1200) = lmat(k,1200)
         mat(k,1202) = mat(k,1202) + lmat(k,1202)
         mat(k,1207) = lmat(k,1207)
         mat(k,1208) = mat(k,1208) + lmat(k,1208)
         mat(k,1211) = mat(k,1211) + lmat(k,1211)
         mat(k,1212) = mat(k,1212) + lmat(k,1212)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1233) = lmat(k,1233)
         mat(k,1234) = lmat(k,1234)
         mat(k,1235) = lmat(k,1235)
         mat(k,1236) = lmat(k,1236)
         mat(k,1237) = mat(k,1237) + lmat(k,1237)
         mat(k,1238) = lmat(k,1238)
         mat(k,1240) = lmat(k,1240)
         mat(k,1243) = lmat(k,1243)
         mat(k,1246) = mat(k,1246) + lmat(k,1246)
         mat(k,1247) = lmat(k,1247)
         mat(k,1248) = lmat(k,1248)
         mat(k,1250) = mat(k,1250) + lmat(k,1250)
         mat(k,1252) = lmat(k,1252)
         mat(k,1254) = mat(k,1254) + lmat(k,1254)
         mat(k,1255) = lmat(k,1255)
         mat(k,1259) = mat(k,1259) + lmat(k,1259)
         mat(k,1275) = mat(k,1275) + lmat(k,1275)
         mat(k,1295) = mat(k,1295) + lmat(k,1295)
         mat(k,1310) = mat(k,1310) + lmat(k,1310)
         mat(k,1311) = mat(k,1311) + lmat(k,1311)
         mat(k,1314) = mat(k,1314) + lmat(k,1314)
         mat(k,1315) = mat(k,1315) + lmat(k,1315)
         mat(k,1319) = mat(k,1319) + lmat(k,1319)
         mat(k,1320) = mat(k,1320) + lmat(k,1320)
         mat(k,1322) = mat(k,1322) + lmat(k,1322)
         mat(k,1323) = mat(k,1323) + lmat(k,1323)
         mat(k,1324) = mat(k,1324) + lmat(k,1324)
         mat(k,1329) = lmat(k,1329)
         mat(k,1341) = mat(k,1341) + lmat(k,1341)
         mat(k,1357) = lmat(k,1357)
         mat(k,1374) = mat(k,1374) + lmat(k,1374)
         mat(k,1386) = mat(k,1386) + lmat(k,1386)
         mat(k,1398) = mat(k,1398) + lmat(k,1398)
         mat(k,1412) = lmat(k,1412)
         mat(k,1414) = mat(k,1414) + lmat(k,1414)
         mat(k,1418) = mat(k,1418) + lmat(k,1418)
         mat(k,1420) = mat(k,1420) + lmat(k,1420)
         mat(k,1424) = lmat(k,1424)
         mat(k,1439) = mat(k,1439) + lmat(k,1439)
         mat(k,1470) = mat(k,1470) + lmat(k,1470)
         mat(k,1491) = mat(k,1491) + lmat(k,1491)
         mat(k,1492) = mat(k,1492) + lmat(k,1492)
         mat(k,1496) = lmat(k,1496)
         mat(k,1503) = mat(k,1503) + lmat(k,1503)
         mat(k,1516) = lmat(k,1516)
         mat(k,1518) = mat(k,1518) + lmat(k,1518)
         mat(k,1524) = mat(k,1524) + lmat(k,1524)
         mat(k,1537) = mat(k,1537) + lmat(k,1537)
         mat(k,1541) = mat(k,1541) + lmat(k,1541)
         mat(k,1544) = mat(k,1544) + lmat(k,1544)
         mat(k,1553) = mat(k,1553) + lmat(k,1553)
         mat(k,1573) = mat(k,1573) + lmat(k,1573)
         mat(k,1575) = mat(k,1575) + lmat(k,1575)
         mat(k,1576) = lmat(k,1576)
         mat(k,1584) = mat(k,1584) + lmat(k,1584)
         mat(k,1588) = mat(k,1588) + lmat(k,1588)
         mat(k,1590) = mat(k,1590) + lmat(k,1590)
         mat(k,1605) = mat(k,1605) + lmat(k,1605)
         mat(k,1607) = mat(k,1607) + lmat(k,1607)
         mat(k,1611) = mat(k,1611) + lmat(k,1611)
         mat(k,1624) = lmat(k,1624)
         mat(k,1625) = lmat(k,1625)
         mat(k,1626) = mat(k,1626) + lmat(k,1626)
         mat(k,1631) = mat(k,1631) + lmat(k,1631)
         mat(k,1632) = mat(k,1632) + lmat(k,1632)
         mat(k,1634) = mat(k,1634) + lmat(k,1634)
         mat(k,1635) = mat(k,1635) + lmat(k,1635)
         mat(k,1636) = lmat(k,1636)
         mat(k,1644) = mat(k,1644) + lmat(k,1644)
         mat(k,1646) = mat(k,1646) + lmat(k,1646)
         mat(k,1684) = mat(k,1684) + lmat(k,1684)
         mat(k,1708) = mat(k,1708) + lmat(k,1708)
         mat(k,1710) = mat(k,1710) + lmat(k,1710)
         mat(k,1720) = lmat(k,1720)
         mat(k,1730) = mat(k,1730) + lmat(k,1730)
         mat(k,1904) = mat(k,1904) + lmat(k,1904)
         mat(k,1926) = mat(k,1926) + lmat(k,1926)
         mat(k,1931) = lmat(k,1931)
         mat(k,1948) = mat(k,1948) + lmat(k,1948)
         mat(k,1995) = mat(k,1995) + lmat(k,1995)
         mat(k,2037) = mat(k,2037) + lmat(k,2037)
         mat(k,2039) = mat(k,2039) + lmat(k,2039)
         mat(k,2053) = mat(k,2053) + lmat(k,2053)
         mat(k,2094) = mat(k,2094) + lmat(k,2094)
         mat(k,2101) = mat(k,2101) + lmat(k,2101)
         mat(k,2103) = mat(k,2103) + lmat(k,2103)
         mat(k,2104) = mat(k,2104) + lmat(k,2104)
         mat(k,2128) = mat(k,2128) + lmat(k,2128)
         mat(k,2129) = mat(k,2129) + lmat(k,2129)
         mat(k,2132) = mat(k,2132) + lmat(k,2132)
         mat(k,2191) = mat(k,2191) + lmat(k,2191)
         mat(k,2193) = lmat(k,2193)
         mat(k,2199) = mat(k,2199) + lmat(k,2199)
         mat(k,2234) = mat(k,2234) + lmat(k,2234)
         mat(k,2239) = mat(k,2239) + lmat(k,2239)
         mat(k,2287) = mat(k,2287) + lmat(k,2287)
         mat(k,2291) = mat(k,2291) + lmat(k,2291)
         mat(k,2294) = mat(k,2294) + lmat(k,2294)
         mat(k,2299) = mat(k,2299) + lmat(k,2299)
         mat(k,2300) = mat(k,2300) + lmat(k,2300)
         mat(k,2303) = mat(k,2303) + lmat(k,2303)
         mat(k,2421) = mat(k,2421) + lmat(k,2421)
         mat(k,2424) = mat(k,2424) + lmat(k,2424)
         mat(k,2428) = mat(k,2428) + lmat(k,2428)
         mat(k,2430) = lmat(k,2430)
         mat(k,2438) = mat(k,2438) + lmat(k,2438)
         mat(k,2448) = mat(k,2448) + lmat(k,2448)
         mat(k,2482) = mat(k,2482) + lmat(k,2482)
         mat(k,2484) = mat(k,2484) + lmat(k,2484)
         mat(k,2485) = mat(k,2485) + lmat(k,2485)
         mat(k,2490) = mat(k,2490) + lmat(k,2490)
         mat(k,2494) = mat(k,2494) + lmat(k,2494)
         mat(k,2501) = lmat(k,2501)
         mat(k,2510) = lmat(k,2510)
         mat(k,2511) = mat(k,2511) + lmat(k,2511)
         mat(k,2512) = lmat(k,2512)
         mat(k,2514) = mat(k,2514) + lmat(k,2514)
         mat(k,2522) = mat(k,2522) + lmat(k,2522)
         mat(k, 256) = 0._r8
         mat(k, 257) = 0._r8
         mat(k, 296) = 0._r8
         mat(k, 355) = 0._r8
         mat(k, 373) = 0._r8
         mat(k, 486) = 0._r8
         mat(k, 487) = 0._r8
         mat(k, 500) = 0._r8
         mat(k, 540) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 557) = 0._r8
         mat(k, 655) = 0._r8
         mat(k, 657) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 692) = 0._r8
         mat(k, 698) = 0._r8
         mat(k, 699) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 709) = 0._r8
         mat(k, 710) = 0._r8
         mat(k, 713) = 0._r8
         mat(k, 745) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 765) = 0._r8
         mat(k, 767) = 0._r8
         mat(k, 768) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 799) = 0._r8
         mat(k, 808) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 815) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 841) = 0._r8
         mat(k, 842) = 0._r8
         mat(k, 862) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 890) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 909) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 914) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 982) = 0._r8
         mat(k, 987) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 997) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1003) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1021) = 0._r8
         mat(k,1022) = 0._r8
         mat(k,1030) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1070) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1077) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1096) = 0._r8
         mat(k,1097) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1163) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1192) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1222) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1241) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1288) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1334) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1348) = 0._r8
         mat(k,1349) = 0._r8
         mat(k,1350) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1360) = 0._r8
         mat(k,1363) = 0._r8
         mat(k,1364) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1370) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1372) = 0._r8
         mat(k,1375) = 0._r8
         mat(k,1376) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1410) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1429) = 0._r8
         mat(k,1432) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1483) = 0._r8
         mat(k,1487) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1490) = 0._r8
         mat(k,1493) = 0._r8
         mat(k,1494) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1507) = 0._r8
         mat(k,1513) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1531) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1547) = 0._r8
         mat(k,1548) = 0._r8
         mat(k,1552) = 0._r8
         mat(k,1556) = 0._r8
         mat(k,1557) = 0._r8
         mat(k,1558) = 0._r8
         mat(k,1559) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1564) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1579) = 0._r8
         mat(k,1581) = 0._r8
         mat(k,1583) = 0._r8
         mat(k,1585) = 0._r8
         mat(k,1587) = 0._r8
         mat(k,1589) = 0._r8
         mat(k,1592) = 0._r8
         mat(k,1595) = 0._r8
         mat(k,1596) = 0._r8
         mat(k,1597) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1608) = 0._r8
         mat(k,1609) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1614) = 0._r8
         mat(k,1617) = 0._r8
         mat(k,1619) = 0._r8
         mat(k,1621) = 0._r8
         mat(k,1628) = 0._r8
         mat(k,1629) = 0._r8
         mat(k,1630) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1640) = 0._r8
         mat(k,1641) = 0._r8
         mat(k,1642) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1654) = 0._r8
         mat(k,1682) = 0._r8
         mat(k,1683) = 0._r8
         mat(k,1685) = 0._r8
         mat(k,1686) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1690) = 0._r8
         mat(k,1691) = 0._r8
         mat(k,1694) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1701) = 0._r8
         mat(k,1702) = 0._r8
         mat(k,1703) = 0._r8
         mat(k,1704) = 0._r8
         mat(k,1705) = 0._r8
         mat(k,1706) = 0._r8
         mat(k,1707) = 0._r8
         mat(k,1709) = 0._r8
         mat(k,1711) = 0._r8
         mat(k,1712) = 0._r8
         mat(k,1713) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1716) = 0._r8
         mat(k,1718) = 0._r8
         mat(k,1719) = 0._r8
         mat(k,1723) = 0._r8
         mat(k,1724) = 0._r8
         mat(k,1726) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1729) = 0._r8
         mat(k,1733) = 0._r8
         mat(k,1734) = 0._r8
         mat(k,1736) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1738) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1741) = 0._r8
         mat(k,1815) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1846) = 0._r8
         mat(k,1849) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1860) = 0._r8
         mat(k,1884) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1924) = 0._r8
         mat(k,1930) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1936) = 0._r8
         mat(k,1943) = 0._r8
         mat(k,1944) = 0._r8
         mat(k,1945) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1975) = 0._r8
         mat(k,1977) = 0._r8
         mat(k,1979) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1981) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1983) = 0._r8
         mat(k,1988) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1996) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2003) = 0._r8
         mat(k,2004) = 0._r8
         mat(k,2034) = 0._r8
         mat(k,2043) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2068) = 0._r8
         mat(k,2069) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2073) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2079) = 0._r8
         mat(k,2080) = 0._r8
         mat(k,2082) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2089) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2111) = 0._r8
         mat(k,2125) = 0._r8
         mat(k,2126) = 0._r8
         mat(k,2130) = 0._r8
         mat(k,2131) = 0._r8
         mat(k,2134) = 0._r8
         mat(k,2138) = 0._r8
         mat(k,2194) = 0._r8
         mat(k,2229) = 0._r8
         mat(k,2231) = 0._r8
         mat(k,2232) = 0._r8
         mat(k,2236) = 0._r8
         mat(k,2244) = 0._r8
         mat(k,2251) = 0._r8
         mat(k,2254) = 0._r8
         mat(k,2258) = 0._r8
         mat(k,2260) = 0._r8
         mat(k,2264) = 0._r8
         mat(k,2271) = 0._r8
         mat(k,2273) = 0._r8
         mat(k,2284) = 0._r8
         mat(k,2285) = 0._r8
         mat(k,2286) = 0._r8
         mat(k,2288) = 0._r8
         mat(k,2289) = 0._r8
         mat(k,2290) = 0._r8
         mat(k,2292) = 0._r8
         mat(k,2295) = 0._r8
         mat(k,2296) = 0._r8
         mat(k,2297) = 0._r8
         mat(k,2298) = 0._r8
         mat(k,2304) = 0._r8
         mat(k,2346) = 0._r8
         mat(k,2347) = 0._r8
         mat(k,2349) = 0._r8
         mat(k,2358) = 0._r8
         mat(k,2375) = 0._r8
         mat(k,2383) = 0._r8
         mat(k,2384) = 0._r8
         mat(k,2386) = 0._r8
         mat(k,2389) = 0._r8
         mat(k,2391) = 0._r8
         mat(k,2395) = 0._r8
         mat(k,2400) = 0._r8
         mat(k,2409) = 0._r8
         mat(k,2411) = 0._r8
         mat(k,2416) = 0._r8
         mat(k,2427) = 0._r8
         mat(k,2429) = 0._r8
         mat(k,2433) = 0._r8
         mat(k,2434) = 0._r8
         mat(k,2435) = 0._r8
         mat(k,2436) = 0._r8
         mat(k,2442) = 0._r8
         mat(k,2443) = 0._r8
         mat(k,2444) = 0._r8
         mat(k,2445) = 0._r8
         mat(k,2449) = 0._r8
         mat(k,2460) = 0._r8
         mat(k,2463) = 0._r8
         mat(k,2467) = 0._r8
         mat(k,2469) = 0._r8
         mat(k,2470) = 0._r8
         mat(k,2471) = 0._r8
         mat(k,2474) = 0._r8
         mat(k,2476) = 0._r8
         mat(k,2477) = 0._r8
         mat(k,2480) = 0._r8
         mat(k,2481) = 0._r8
         mat(k,2483) = 0._r8
         mat(k,2486) = 0._r8
         mat(k,2487) = 0._r8
         mat(k,2493) = 0._r8
         mat(k,2495) = 0._r8
         mat(k,2500) = 0._r8
         mat(k,2502) = 0._r8
         mat(k,2503) = 0._r8
         mat(k,2504) = 0._r8
         mat(k,2505) = 0._r8
         mat(k,2506) = 0._r8
         mat(k,2507) = 0._r8
         mat(k,2508) = 0._r8
         mat(k,2509) = 0._r8
         mat(k,2513) = 0._r8
         mat(k,2515) = 0._r8
         mat(k,2516) = 0._r8
         mat(k,2517) = 0._r8
         mat(k,2518) = 0._r8
         mat(k,2519) = 0._r8
         mat(k,2520) = 0._r8
         mat(k,2521) = 0._r8
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
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 72) = mat(k, 72) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 82) = mat(k, 82) - dti(k)
         mat(k, 88) = mat(k, 88) - dti(k)
         mat(k, 94) = mat(k, 94) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 98) = mat(k, 98) - dti(k)
         mat(k, 99) = mat(k, 99) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 107) = mat(k, 107) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 114) = mat(k, 114) - dti(k)
         mat(k, 120) = mat(k, 120) - dti(k)
         mat(k, 126) = mat(k, 126) - dti(k)
         mat(k, 132) = mat(k, 132) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 145) = mat(k, 145) - dti(k)
         mat(k, 148) = mat(k, 148) - dti(k)
         mat(k, 151) = mat(k, 151) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 159) = mat(k, 159) - dti(k)
         mat(k, 163) = mat(k, 163) - dti(k)
         mat(k, 167) = mat(k, 167) - dti(k)
         mat(k, 171) = mat(k, 171) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 182) = mat(k, 182) - dti(k)
         mat(k, 187) = mat(k, 187) - dti(k)
         mat(k, 192) = mat(k, 192) - dti(k)
         mat(k, 197) = mat(k, 197) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 210) = mat(k, 210) - dti(k)
         mat(k, 214) = mat(k, 214) - dti(k)
         mat(k, 219) = mat(k, 219) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 224) = mat(k, 224) - dti(k)
         mat(k, 226) = mat(k, 226) - dti(k)
         mat(k, 231) = mat(k, 231) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 243) = mat(k, 243) - dti(k)
         mat(k, 247) = mat(k, 247) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 260) = mat(k, 260) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 268) = mat(k, 268) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 274) = mat(k, 274) - dti(k)
         mat(k, 277) = mat(k, 277) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 290) = mat(k, 290) - dti(k)
         mat(k, 295) = mat(k, 295) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 305) = mat(k, 305) - dti(k)
         mat(k, 309) = mat(k, 309) - dti(k)
         mat(k, 313) = mat(k, 313) - dti(k)
         mat(k, 317) = mat(k, 317) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 329) = mat(k, 329) - dti(k)
         mat(k, 335) = mat(k, 335) - dti(k)
         mat(k, 338) = mat(k, 338) - dti(k)
         mat(k, 344) = mat(k, 344) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 354) = mat(k, 354) - dti(k)
         mat(k, 360) = mat(k, 360) - dti(k)
         mat(k, 365) = mat(k, 365) - dti(k)
         mat(k, 371) = mat(k, 371) - dti(k)
         mat(k, 376) = mat(k, 376) - dti(k)
         mat(k, 381) = mat(k, 381) - dti(k)
         mat(k, 384) = mat(k, 384) - dti(k)
         mat(k, 389) = mat(k, 389) - dti(k)
         mat(k, 394) = mat(k, 394) - dti(k)
         mat(k, 399) = mat(k, 399) - dti(k)
         mat(k, 407) = mat(k, 407) - dti(k)
         mat(k, 415) = mat(k, 415) - dti(k)
         mat(k, 423) = mat(k, 423) - dti(k)
         mat(k, 429) = mat(k, 429) - dti(k)
         mat(k, 435) = mat(k, 435) - dti(k)
         mat(k, 441) = mat(k, 441) - dti(k)
         mat(k, 447) = mat(k, 447) - dti(k)
         mat(k, 453) = mat(k, 453) - dti(k)
         mat(k, 459) = mat(k, 459) - dti(k)
         mat(k, 465) = mat(k, 465) - dti(k)
         mat(k, 471) = mat(k, 471) - dti(k)
         mat(k, 477) = mat(k, 477) - dti(k)
         mat(k, 485) = mat(k, 485) - dti(k)
         mat(k, 491) = mat(k, 491) - dti(k)
         mat(k, 498) = mat(k, 498) - dti(k)
         mat(k, 504) = mat(k, 504) - dti(k)
         mat(k, 507) = mat(k, 507) - dti(k)
         mat(k, 512) = mat(k, 512) - dti(k)
         mat(k, 519) = mat(k, 519) - dti(k)
         mat(k, 526) = mat(k, 526) - dti(k)
         mat(k, 530) = mat(k, 530) - dti(k)
         mat(k, 539) = mat(k, 539) - dti(k)
         mat(k, 546) = mat(k, 546) - dti(k)
         mat(k, 554) = mat(k, 554) - dti(k)
         mat(k, 561) = mat(k, 561) - dti(k)
         mat(k, 567) = mat(k, 567) - dti(k)
         mat(k, 573) = mat(k, 573) - dti(k)
         mat(k, 578) = mat(k, 578) - dti(k)
         mat(k, 586) = mat(k, 586) - dti(k)
         mat(k, 594) = mat(k, 594) - dti(k)
         mat(k, 602) = mat(k, 602) - dti(k)
         mat(k, 610) = mat(k, 610) - dti(k)
         mat(k, 618) = mat(k, 618) - dti(k)
         mat(k, 626) = mat(k, 626) - dti(k)
         mat(k, 635) = mat(k, 635) - dti(k)
         mat(k, 642) = mat(k, 642) - dti(k)
         mat(k, 653) = mat(k, 653) - dti(k)
         mat(k, 662) = mat(k, 662) - dti(k)
         mat(k, 671) = mat(k, 671) - dti(k)
         mat(k, 678) = mat(k, 678) - dti(k)
         mat(k, 682) = mat(k, 682) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 697) = mat(k, 697) - dti(k)
         mat(k, 708) = mat(k, 708) - dti(k)
         mat(k, 719) = mat(k, 719) - dti(k)
         mat(k, 726) = mat(k, 726) - dti(k)
         mat(k, 737) = mat(k, 737) - dti(k)
         mat(k, 747) = mat(k, 747) - dti(k)
         mat(k, 753) = mat(k, 753) - dti(k)
         mat(k, 764) = mat(k, 764) - dti(k)
         mat(k, 775) = mat(k, 775) - dti(k)
         mat(k, 782) = mat(k, 782) - dti(k)
         mat(k, 793) = mat(k, 793) - dti(k)
         mat(k, 809) = mat(k, 809) - dti(k)
         mat(k, 820) = mat(k, 820) - dti(k)
         mat(k, 829) = mat(k, 829) - dti(k)
         mat(k, 839) = mat(k, 839) - dti(k)
         mat(k, 847) = mat(k, 847) - dti(k)
         mat(k, 855) = mat(k, 855) - dti(k)
         mat(k, 860) = mat(k, 860) - dti(k)
         mat(k, 870) = mat(k, 870) - dti(k)
         mat(k, 879) = mat(k, 879) - dti(k)
         mat(k, 887) = mat(k, 887) - dti(k)
         mat(k, 895) = mat(k, 895) - dti(k)
         mat(k, 907) = mat(k, 907) - dti(k)
         mat(k, 915) = mat(k, 915) - dti(k)
         mat(k, 924) = mat(k, 924) - dti(k)
         mat(k, 941) = mat(k, 941) - dti(k)
         mat(k, 953) = mat(k, 953) - dti(k)
         mat(k, 962) = mat(k, 962) - dti(k)
         mat(k, 967) = mat(k, 967) - dti(k)
         mat(k, 977) = mat(k, 977) - dti(k)
         mat(k, 989) = mat(k, 989) - dti(k)
         mat(k,1001) = mat(k,1001) - dti(k)
         mat(k,1020) = mat(k,1020) - dti(k)
         mat(k,1048) = mat(k,1048) - dti(k)
         mat(k,1072) = mat(k,1072) - dti(k)
         mat(k,1084) = mat(k,1084) - dti(k)
         mat(k,1092) = mat(k,1092) - dti(k)
         mat(k,1102) = mat(k,1102) - dti(k)
         mat(k,1110) = mat(k,1110) - dti(k)
         mat(k,1118) = mat(k,1118) - dti(k)
         mat(k,1132) = mat(k,1132) - dti(k)
         mat(k,1145) = mat(k,1145) - dti(k)
         mat(k,1159) = mat(k,1159) - dti(k)
         mat(k,1175) = mat(k,1175) - dti(k)
         mat(k,1193) = mat(k,1193) - dti(k)
         mat(k,1202) = mat(k,1202) - dti(k)
         mat(k,1208) = mat(k,1208) - dti(k)
         mat(k,1220) = mat(k,1220) - dti(k)
         mat(k,1237) = mat(k,1237) - dti(k)
         mat(k,1250) = mat(k,1250) - dti(k)
         mat(k,1259) = mat(k,1259) - dti(k)
         mat(k,1275) = mat(k,1275) - dti(k)
         mat(k,1295) = mat(k,1295) - dti(k)
         mat(k,1311) = mat(k,1311) - dti(k)
         mat(k,1323) = mat(k,1323) - dti(k)
         mat(k,1341) = mat(k,1341) - dti(k)
         mat(k,1374) = mat(k,1374) - dti(k)
         mat(k,1398) = mat(k,1398) - dti(k)
         mat(k,1418) = mat(k,1418) - dti(k)
         mat(k,1439) = mat(k,1439) - dti(k)
         mat(k,1470) = mat(k,1470) - dti(k)
         mat(k,1492) = mat(k,1492) - dti(k)
         mat(k,1503) = mat(k,1503) - dti(k)
         mat(k,1518) = mat(k,1518) - dti(k)
         mat(k,1537) = mat(k,1537) - dti(k)
         mat(k,1553) = mat(k,1553) - dti(k)
         mat(k,1584) = mat(k,1584) - dti(k)
         mat(k,1607) = mat(k,1607) - dti(k)
         mat(k,1631) = mat(k,1631) - dti(k)
         mat(k,1684) = mat(k,1684) - dti(k)
         mat(k,1708) = mat(k,1708) - dti(k)
         mat(k,1730) = mat(k,1730) - dti(k)
         mat(k,1904) = mat(k,1904) - dti(k)
         mat(k,1948) = mat(k,1948) - dti(k)
         mat(k,1995) = mat(k,1995) - dti(k)
         mat(k,2039) = mat(k,2039) - dti(k)
         mat(k,2104) = mat(k,2104) - dti(k)
         mat(k,2132) = mat(k,2132) - dti(k)
         mat(k,2239) = mat(k,2239) - dti(k)
         mat(k,2300) = mat(k,2300) - dti(k)
         mat(k,2421) = mat(k,2421) - dti(k)
         mat(k,2448) = mat(k,2448) - dti(k)
         mat(k,2494) = mat(k,2494) - dti(k)
         mat(k,2522) = mat(k,2522) - dti(k)
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
