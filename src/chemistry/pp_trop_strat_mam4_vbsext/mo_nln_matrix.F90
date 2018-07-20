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
         mat(k,566) = -(rxt(k,376)*y(k,241))
         mat(k,1638) = -rxt(k,376)*y(k,1)
         mat(k,1354) = rxt(k,379)*y(k,218)
         mat(k,849) = rxt(k,379)*y(k,122)
         mat(k,577) = -(rxt(k,380)*y(k,241))
         mat(k,1639) = -rxt(k,380)*y(k,2)
         mat(k,850) = rxt(k,377)*y(k,230)
         mat(k,1466) = rxt(k,377)*y(k,218)
         mat(k,831) = -(rxt(k,459)*y(k,124) + rxt(k,460)*y(k,132) + rxt(k,461) &
                      *y(k,241))
         mat(k,1785) = -rxt(k,459)*y(k,5)
         mat(k,1898) = -rxt(k,460)*y(k,5)
         mat(k,1662) = -rxt(k,461)*y(k,5)
         mat(k,123) = -(rxt(k,418)*y(k,241))
         mat(k,1574) = -rxt(k,418)*y(k,6)
         mat(k,346) = -(rxt(k,421)*y(k,241))
         mat(k,1608) = -rxt(k,421)*y(k,7)
         mat(k,418) = rxt(k,419)*y(k,230)
         mat(k,1445) = rxt(k,419)*y(k,219)
         mat(k,124) = .120_r8*rxt(k,418)*y(k,241)
         mat(k,1575) = .120_r8*rxt(k,418)*y(k,6)
         mat(k,828) = .100_r8*rxt(k,460)*y(k,132)
         mat(k,802) = .100_r8*rxt(k,463)*y(k,132)
         mat(k,1886) = .100_r8*rxt(k,460)*y(k,5) + .100_r8*rxt(k,463)*y(k,110)
         mat(k,1342) = .500_r8*rxt(k,420)*y(k,219) + .200_r8*rxt(k,447)*y(k,247) &
                      + .060_r8*rxt(k,453)*y(k,249)
         mat(k,419) = .500_r8*rxt(k,420)*y(k,122)
         mat(k,634) = .200_r8*rxt(k,447)*y(k,122)
         mat(k,650) = .060_r8*rxt(k,453)*y(k,122)
         mat(k,1335) = .200_r8*rxt(k,447)*y(k,247) + .200_r8*rxt(k,453)*y(k,249)
         mat(k,633) = .200_r8*rxt(k,447)*y(k,122)
         mat(k,648) = .200_r8*rxt(k,453)*y(k,122)
         mat(k,1351) = .200_r8*rxt(k,447)*y(k,247) + .150_r8*rxt(k,453)*y(k,249)
         mat(k,636) = .200_r8*rxt(k,447)*y(k,122)
         mat(k,651) = .150_r8*rxt(k,453)*y(k,122)
         mat(k,1336) = .210_r8*rxt(k,453)*y(k,249)
         mat(k,649) = .210_r8*rxt(k,453)*y(k,122)
         mat(k,194) = -(rxt(k,381)*y(k,241))
         mat(k,1586) = -rxt(k,381)*y(k,14)
         mat(k,827) = .050_r8*rxt(k,460)*y(k,132)
         mat(k,801) = .050_r8*rxt(k,463)*y(k,132)
         mat(k,1885) = .050_r8*rxt(k,460)*y(k,5) + .050_r8*rxt(k,463)*y(k,110)
         mat(k,294) = -(rxt(k,347)*y(k,124) + rxt(k,348)*y(k,241))
         mat(k,1778) = -rxt(k,347)*y(k,15)
         mat(k,1601) = -rxt(k,348)*y(k,15)
         mat(k,1294) = -(rxt(k,230)*y(k,41) + rxt(k,231)*y(k,230) + rxt(k,232) &
                      *y(k,132))
         mat(k,1714) = -rxt(k,230)*y(k,16)
         mat(k,1508) = -rxt(k,231)*y(k,16)
         mat(k,1922) = -rxt(k,232)*y(k,16)
         mat(k,1946) = 4.000_r8*rxt(k,233)*y(k,18) + (rxt(k,234)+rxt(k,235))*y(k,58) &
                      + rxt(k,238)*y(k,122) + rxt(k,241)*y(k,131) + rxt(k,486) &
                      *y(k,149) + rxt(k,242)*y(k,241)
         mat(k,1972) = (rxt(k,234)+rxt(k,235))*y(k,18)
         mat(k,721) = rxt(k,243)*y(k,131) + rxt(k,249)*y(k,240) + rxt(k,244)*y(k,241)
         mat(k,1392) = rxt(k,238)*y(k,18)
         mat(k,2002) = rxt(k,241)*y(k,18) + rxt(k,243)*y(k,80)
         mat(k,1128) = rxt(k,486)*y(k,18)
         mat(k,1532) = rxt(k,249)*y(k,80)
         mat(k,1691) = rxt(k,242)*y(k,18) + rxt(k,244)*y(k,80)
         mat(k,1940) = rxt(k,236)*y(k,58)
         mat(k,1966) = rxt(k,236)*y(k,18)
         mat(k,1411) = (rxt(k,538)+rxt(k,543))*y(k,90)
         mat(k,683) = (rxt(k,538)+rxt(k,543))*y(k,84)
         mat(k,1959) = -(4._r8*rxt(k,233)*y(k,18) + (rxt(k,234) + rxt(k,235) + rxt(k,236) &
                      ) * y(k,58) + rxt(k,237)*y(k,230) + rxt(k,238)*y(k,122) &
                      + rxt(k,239)*y(k,123) + rxt(k,241)*y(k,131) + rxt(k,242) &
                      *y(k,241) + rxt(k,486)*y(k,149))
         mat(k,1985) = -(rxt(k,234) + rxt(k,235) + rxt(k,236)) * y(k,18)
         mat(k,1521) = -rxt(k,237)*y(k,18)
         mat(k,1405) = -rxt(k,238)*y(k,18)
         mat(k,1768) = -rxt(k,239)*y(k,18)
         mat(k,2015) = -rxt(k,241)*y(k,18)
         mat(k,1704) = -rxt(k,242)*y(k,18)
         mat(k,1135) = -rxt(k,486)*y(k,18)
         mat(k,1301) = rxt(k,232)*y(k,132)
         mat(k,509) = rxt(k,240)*y(k,131)
         mat(k,725) = rxt(k,250)*y(k,240)
         mat(k,688) = rxt(k,245)*y(k,131)
         mat(k,2015) = mat(k,2015) + rxt(k,240)*y(k,19) + rxt(k,245)*y(k,90)
         mat(k,1935) = rxt(k,232)*y(k,16)
         mat(k,1545) = rxt(k,250)*y(k,80)
         mat(k,503) = -(rxt(k,240)*y(k,131))
         mat(k,1992) = -rxt(k,240)*y(k,19)
         mat(k,1942) = rxt(k,239)*y(k,123)
         mat(k,1740) = rxt(k,239)*y(k,18)
         mat(k,200) = -(rxt(k,422)*y(k,241))
         mat(k,1587) = -rxt(k,422)*y(k,21)
         mat(k,1333) = rxt(k,425)*y(k,220)
         mat(k,376) = rxt(k,425)*y(k,122)
         mat(k,271) = -(rxt(k,424)*y(k,241))
         mat(k,1596) = -rxt(k,424)*y(k,22)
         mat(k,377) = rxt(k,423)*y(k,230)
         mat(k,1437) = rxt(k,423)*y(k,220)
         mat(k,237) = -(rxt(k,296)*y(k,55) + rxt(k,297)*y(k,241))
         mat(k,2021) = -rxt(k,296)*y(k,23)
         mat(k,1592) = -rxt(k,297)*y(k,23)
         mat(k,487) = -(rxt(k,298)*y(k,55) + rxt(k,299)*y(k,132) + rxt(k,324)*y(k,241))
         mat(k,2023) = -rxt(k,298)*y(k,24)
         mat(k,1890) = -rxt(k,299)*y(k,24)
         mat(k,1629) = -rxt(k,324)*y(k,24)
         mat(k,216) = -(rxt(k,304)*y(k,241))
         mat(k,1590) = -rxt(k,304)*y(k,25)
         mat(k,736) = .800_r8*rxt(k,300)*y(k,221) + .200_r8*rxt(k,301)*y(k,225)
         mat(k,1830) = .200_r8*rxt(k,301)*y(k,221)
         mat(k,276) = -(rxt(k,305)*y(k,241))
         mat(k,1597) = -rxt(k,305)*y(k,26)
         mat(k,737) = rxt(k,302)*y(k,230)
         mat(k,1438) = rxt(k,302)*y(k,221)
         mat(k,243) = -(rxt(k,306)*y(k,55) + rxt(k,307)*y(k,241))
         mat(k,2022) = -rxt(k,306)*y(k,27)
         mat(k,1593) = -rxt(k,307)*y(k,27)
         mat(k,893) = -(rxt(k,327)*y(k,124) + rxt(k,328)*y(k,132) + rxt(k,345) &
                      *y(k,241))
         mat(k,1789) = -rxt(k,327)*y(k,28)
         mat(k,1902) = -rxt(k,328)*y(k,28)
         mat(k,1667) = -rxt(k,345)*y(k,28)
         mat(k,767) = .130_r8*rxt(k,405)*y(k,132)
         mat(k,1902) = mat(k,1902) + .130_r8*rxt(k,405)*y(k,97)
         mat(k,340) = -(rxt(k,332)*y(k,241))
         mat(k,1607) = -rxt(k,332)*y(k,29)
         mat(k,709) = rxt(k,330)*y(k,230)
         mat(k,1444) = rxt(k,330)*y(k,222)
         mat(k,99) = -(rxt(k,333)*y(k,241))
         mat(k,1571) = -rxt(k,333)*y(k,30)
         mat(k,220) = -(rxt(k,428)*y(k,241))
         mat(k,1591) = -rxt(k,428)*y(k,31)
         mat(k,557) = rxt(k,426)*y(k,230)
         mat(k,1434) = rxt(k,426)*y(k,223)
         mat(k,1722) = -(rxt(k,194)*y(k,55) + rxt(k,230)*y(k,16) + rxt(k,274)*y(k,230) &
                      + rxt(k,275)*y(k,124) + rxt(k,276)*y(k,131) + rxt(k,277) &
                      *y(k,241))
         mat(k,2044) = -rxt(k,194)*y(k,41)
         mat(k,1299) = -rxt(k,230)*y(k,41)
         mat(k,1516) = -rxt(k,274)*y(k,41)
         mat(k,1820) = -rxt(k,275)*y(k,41)
         mat(k,2010) = -rxt(k,276)*y(k,41)
         mat(k,1699) = -rxt(k,277)*y(k,41)
         mat(k,574) = .400_r8*rxt(k,376)*y(k,241)
         mat(k,844) = .340_r8*rxt(k,460)*y(k,132)
         mat(k,299) = .500_r8*rxt(k,347)*y(k,124)
         mat(k,493) = rxt(k,299)*y(k,132)
         mat(k,902) = .500_r8*rxt(k,328)*y(k,132)
         mat(k,438) = .500_r8*rxt(k,316)*y(k,241)
         mat(k,704) = rxt(k,282)*y(k,241)
         mat(k,355) = .300_r8*rxt(k,283)*y(k,241)
         mat(k,1980) = rxt(k,201)*y(k,225)
         mat(k,928) = .800_r8*rxt(k,321)*y(k,241)
         mat(k,777) = .910_r8*rxt(k,405)*y(k,132)
         mat(k,518) = .300_r8*rxt(k,396)*y(k,241)
         mat(k,1098) = .800_r8*rxt(k,400)*y(k,225)
         mat(k,1112) = .120_r8*rxt(k,358)*y(k,132)
         mat(k,484) = .500_r8*rxt(k,371)*y(k,241)
         mat(k,818) = .340_r8*rxt(k,463)*y(k,132)
         mat(k,1182) = .600_r8*rxt(k,372)*y(k,132)
         mat(k,1400) = .100_r8*rxt(k,378)*y(k,218) + rxt(k,281)*y(k,225) &
                      + .500_r8*rxt(k,349)*y(k,227) + .500_r8*rxt(k,318)*y(k,229) &
                      + .920_r8*rxt(k,388)*y(k,232) + .250_r8*rxt(k,356)*y(k,234) &
                      + rxt(k,365)*y(k,236) + rxt(k,339)*y(k,243) + rxt(k,343) &
                      *y(k,244) + .340_r8*rxt(k,472)*y(k,245) + .320_r8*rxt(k,477) &
                      *y(k,246) + .250_r8*rxt(k,413)*y(k,248)
         mat(k,1820) = mat(k,1820) + .500_r8*rxt(k,347)*y(k,15) + rxt(k,389)*y(k,232) &
                      + .250_r8*rxt(k,355)*y(k,234) + rxt(k,366)*y(k,236)
         mat(k,1930) = .340_r8*rxt(k,460)*y(k,5) + rxt(k,299)*y(k,24) &
                      + .500_r8*rxt(k,328)*y(k,28) + .910_r8*rxt(k,405)*y(k,97) &
                      + .120_r8*rxt(k,358)*y(k,105) + .340_r8*rxt(k,463)*y(k,110) &
                      + .600_r8*rxt(k,372)*y(k,111)
         mat(k,407) = rxt(k,323)*y(k,241)
         mat(k,952) = .680_r8*rxt(k,481)*y(k,241)
         mat(k,860) = .100_r8*rxt(k,378)*y(k,122)
         mat(k,744) = .700_r8*rxt(k,301)*y(k,225)
         mat(k,716) = rxt(k,329)*y(k,225)
         mat(k,1286) = rxt(k,312)*y(k,225) + rxt(k,385)*y(k,232) + .250_r8*rxt(k,352) &
                      *y(k,234) + rxt(k,361)*y(k,236) + .250_r8*rxt(k,410)*y(k,248)
         mat(k,1870) = rxt(k,201)*y(k,58) + .800_r8*rxt(k,400)*y(k,100) + rxt(k,281) &
                      *y(k,122) + .700_r8*rxt(k,301)*y(k,221) + rxt(k,329)*y(k,222) &
                      + rxt(k,312)*y(k,224) + (4.000_r8*rxt(k,278)+2.000_r8*rxt(k,279)) &
                      *y(k,225) + 1.500_r8*rxt(k,386)*y(k,232) + .750_r8*rxt(k,391) &
                      *y(k,233) + .880_r8*rxt(k,353)*y(k,234) + 2.000_r8*rxt(k,362) &
                      *y(k,236) + .750_r8*rxt(k,465)*y(k,239) + .800_r8*rxt(k,341) &
                      *y(k,244) + .930_r8*rxt(k,470)*y(k,245) + .950_r8*rxt(k,475) &
                      *y(k,246) + .800_r8*rxt(k,411)*y(k,248)
         mat(k,501) = .500_r8*rxt(k,349)*y(k,122)
         mat(k,625) = .500_r8*rxt(k,318)*y(k,122)
         mat(k,1516) = mat(k,1516) + .450_r8*rxt(k,363)*y(k,236) + .150_r8*rxt(k,342) &
                      *y(k,244)
         mat(k,1162) = .920_r8*rxt(k,388)*y(k,122) + rxt(k,389)*y(k,124) + rxt(k,385) &
                      *y(k,224) + 1.500_r8*rxt(k,386)*y(k,225)
         mat(k,1238) = .750_r8*rxt(k,391)*y(k,225)
         mat(k,1205) = .250_r8*rxt(k,356)*y(k,122) + .250_r8*rxt(k,355)*y(k,124) &
                      + .250_r8*rxt(k,352)*y(k,224) + .880_r8*rxt(k,353)*y(k,225)
         mat(k,1256) = rxt(k,365)*y(k,122) + rxt(k,366)*y(k,124) + rxt(k,361)*y(k,224) &
                      + 2.000_r8*rxt(k,362)*y(k,225) + .450_r8*rxt(k,363)*y(k,230) &
                      + 4.000_r8*rxt(k,364)*y(k,236)
         mat(k,1021) = .750_r8*rxt(k,465)*y(k,225)
         mat(k,1699) = mat(k,1699) + .400_r8*rxt(k,376)*y(k,1) + .500_r8*rxt(k,316) &
                      *y(k,50) + rxt(k,282)*y(k,51) + .300_r8*rxt(k,283)*y(k,52) &
                      + .800_r8*rxt(k,321)*y(k,73) + .300_r8*rxt(k,396)*y(k,98) &
                      + .500_r8*rxt(k,371)*y(k,109) + rxt(k,323)*y(k,136) &
                      + .680_r8*rxt(k,481)*y(k,207)
         mat(k,680) = rxt(k,339)*y(k,122)
         mat(k,1035) = rxt(k,343)*y(k,122) + .800_r8*rxt(k,341)*y(k,225) &
                      + .150_r8*rxt(k,342)*y(k,230)
         mat(k,1002) = .340_r8*rxt(k,472)*y(k,122) + .930_r8*rxt(k,470)*y(k,225)
         mat(k,982) = .320_r8*rxt(k,477)*y(k,122) + .950_r8*rxt(k,475)*y(k,225)
         mat(k,1075) = .250_r8*rxt(k,413)*y(k,122) + .250_r8*rxt(k,410)*y(k,224) &
                      + .800_r8*rxt(k,411)*y(k,225)
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
         mat(k,1040) = -(rxt(k,308)*y(k,124) + rxt(k,309)*y(k,241))
         mat(k,1800) = -rxt(k,308)*y(k,44)
         mat(k,1678) = -rxt(k,309)*y(k,44)
         mat(k,570) = .800_r8*rxt(k,376)*y(k,241)
         mat(k,297) = rxt(k,347)*y(k,124)
         mat(k,217) = rxt(k,304)*y(k,241)
         mat(k,278) = .500_r8*rxt(k,305)*y(k,241)
         mat(k,896) = .500_r8*rxt(k,328)*y(k,132)
         mat(k,1172) = .100_r8*rxt(k,372)*y(k,132)
         mat(k,1381) = .400_r8*rxt(k,378)*y(k,218) + rxt(k,303)*y(k,221) &
                      + .270_r8*rxt(k,331)*y(k,222) + rxt(k,349)*y(k,227) + rxt(k,368) &
                      *y(k,238) + rxt(k,339)*y(k,243)
         mat(k,1800) = mat(k,1800) + rxt(k,347)*y(k,15)
         mat(k,1911) = .500_r8*rxt(k,328)*y(k,28) + .100_r8*rxt(k,372)*y(k,111)
         mat(k,855) = .400_r8*rxt(k,378)*y(k,122)
         mat(k,740) = rxt(k,303)*y(k,122) + 3.200_r8*rxt(k,300)*y(k,221) &
                      + .800_r8*rxt(k,301)*y(k,225)
         mat(k,712) = .270_r8*rxt(k,331)*y(k,122)
         mat(k,1852) = .800_r8*rxt(k,301)*y(k,221)
         mat(k,498) = rxt(k,349)*y(k,122)
         mat(k,1496) = .200_r8*rxt(k,367)*y(k,238)
         mat(k,589) = rxt(k,368)*y(k,122) + .200_r8*rxt(k,367)*y(k,230)
         mat(k,1678) = mat(k,1678) + .800_r8*rxt(k,376)*y(k,1) + rxt(k,304)*y(k,25) &
                      + .500_r8*rxt(k,305)*y(k,26)
         mat(k,676) = rxt(k,339)*y(k,122)
         mat(k,93) = -(rxt(k,310)*y(k,241))
         mat(k,1570) = -rxt(k,310)*y(k,46)
         mat(k,863) = -(rxt(k,346)*y(k,241))
         mat(k,1664) = -rxt(k,346)*y(k,47)
         mat(k,569) = .800_r8*rxt(k,376)*y(k,241)
         mat(k,833) = .520_r8*rxt(k,460)*y(k,132)
         mat(k,296) = .500_r8*rxt(k,347)*y(k,124)
         mat(k,807) = .520_r8*rxt(k,463)*y(k,132)
         mat(k,1369) = .250_r8*rxt(k,378)*y(k,218) + .820_r8*rxt(k,331)*y(k,222) &
                      + .500_r8*rxt(k,349)*y(k,227) + .270_r8*rxt(k,472)*y(k,245) &
                      + .040_r8*rxt(k,477)*y(k,246)
         mat(k,1787) = .500_r8*rxt(k,347)*y(k,15)
         mat(k,1900) = .520_r8*rxt(k,460)*y(k,5) + .520_r8*rxt(k,463)*y(k,110)
         mat(k,945) = .500_r8*rxt(k,481)*y(k,241)
         mat(k,854) = .250_r8*rxt(k,378)*y(k,122)
         mat(k,711) = .820_r8*rxt(k,331)*y(k,122) + .820_r8*rxt(k,329)*y(k,225)
         mat(k,1841) = .820_r8*rxt(k,329)*y(k,222) + .150_r8*rxt(k,470)*y(k,245) &
                      + .025_r8*rxt(k,475)*y(k,246)
         mat(k,496) = .500_r8*rxt(k,349)*y(k,122)
         mat(k,1664) = mat(k,1664) + .800_r8*rxt(k,376)*y(k,1) + .500_r8*rxt(k,481) &
                      *y(k,207)
         mat(k,990) = .270_r8*rxt(k,472)*y(k,122) + .150_r8*rxt(k,470)*y(k,225)
         mat(k,968) = .040_r8*rxt(k,477)*y(k,122) + .025_r8*rxt(k,475)*y(k,225)
         mat(k,1116) = -(rxt(k,334)*y(k,124) + rxt(k,335)*y(k,241))
         mat(k,1804) = -rxt(k,334)*y(k,48)
         mat(k,1683) = -rxt(k,335)*y(k,48)
         mat(k,960) = rxt(k,336)*y(k,241)
         mat(k,1105) = .880_r8*rxt(k,358)*y(k,132)
         mat(k,1173) = .500_r8*rxt(k,372)*y(k,132)
         mat(k,1385) = .170_r8*rxt(k,431)*y(k,226) + .050_r8*rxt(k,394)*y(k,233) &
                      + .250_r8*rxt(k,356)*y(k,234) + .170_r8*rxt(k,437)*y(k,237) &
                      + .400_r8*rxt(k,447)*y(k,247) + .250_r8*rxt(k,413)*y(k,248) &
                      + .540_r8*rxt(k,453)*y(k,249) + .510_r8*rxt(k,456)*y(k,250)
         mat(k,1804) = mat(k,1804) + .050_r8*rxt(k,395)*y(k,233) + .250_r8*rxt(k,355) &
                      *y(k,234) + .250_r8*rxt(k,414)*y(k,248)
         mat(k,748) = rxt(k,337)*y(k,241)
         mat(k,1914) = .880_r8*rxt(k,358)*y(k,105) + .500_r8*rxt(k,372)*y(k,111)
         mat(k,1274) = .250_r8*rxt(k,352)*y(k,234) + .250_r8*rxt(k,410)*y(k,248)
         mat(k,1856) = .240_r8*rxt(k,353)*y(k,234) + .500_r8*rxt(k,341)*y(k,244) &
                      + .100_r8*rxt(k,411)*y(k,248)
         mat(k,667) = .170_r8*rxt(k,431)*y(k,122) + .070_r8*rxt(k,430)*y(k,230)
         mat(k,1501) = .070_r8*rxt(k,430)*y(k,226) + .070_r8*rxt(k,436)*y(k,237)
         mat(k,1227) = .050_r8*rxt(k,394)*y(k,122) + .050_r8*rxt(k,395)*y(k,124)
         mat(k,1196) = .250_r8*rxt(k,356)*y(k,122) + .250_r8*rxt(k,355)*y(k,124) &
                      + .250_r8*rxt(k,352)*y(k,224) + .240_r8*rxt(k,353)*y(k,225)
         mat(k,785) = .170_r8*rxt(k,437)*y(k,122) + .070_r8*rxt(k,436)*y(k,230)
         mat(k,1683) = mat(k,1683) + rxt(k,336)*y(k,94) + rxt(k,337)*y(k,125)
         mat(k,1030) = .500_r8*rxt(k,341)*y(k,225)
         mat(k,643) = .400_r8*rxt(k,447)*y(k,122)
         mat(k,1069) = .250_r8*rxt(k,413)*y(k,122) + .250_r8*rxt(k,414)*y(k,124) &
                      + .250_r8*rxt(k,410)*y(k,224) + .100_r8*rxt(k,411)*y(k,225)
         mat(k,659) = .540_r8*rxt(k,453)*y(k,122)
         mat(k,430) = .510_r8*rxt(k,456)*y(k,122)
         mat(k,467) = -(rxt(k,315)*y(k,241))
         mat(k,1626) = -rxt(k,315)*y(k,49)
         mat(k,889) = .120_r8*rxt(k,328)*y(k,132)
         mat(k,1889) = .120_r8*rxt(k,328)*y(k,28)
         mat(k,1265) = .100_r8*rxt(k,312)*y(k,225) + .150_r8*rxt(k,313)*y(k,230)
         mat(k,1834) = .100_r8*rxt(k,312)*y(k,224)
         mat(k,1460) = .150_r8*rxt(k,313)*y(k,224) + .150_r8*rxt(k,363)*y(k,236)
         mat(k,1245) = .150_r8*rxt(k,363)*y(k,230)
         mat(k,435) = -(rxt(k,316)*y(k,241))
         mat(k,1621) = -rxt(k,316)*y(k,50)
         mat(k,1264) = .400_r8*rxt(k,313)*y(k,230)
         mat(k,1457) = .400_r8*rxt(k,313)*y(k,224) + .400_r8*rxt(k,363)*y(k,236)
         mat(k,1244) = .400_r8*rxt(k,363)*y(k,230)
         mat(k,701) = -(rxt(k,282)*y(k,241))
         mat(k,1650) = -rxt(k,282)*y(k,51)
         mat(k,1082) = .200_r8*rxt(k,400)*y(k,225)
         mat(k,738) = .300_r8*rxt(k,301)*y(k,225)
         mat(k,1836) = .200_r8*rxt(k,400)*y(k,100) + .300_r8*rxt(k,301)*y(k,221) &
                      + 2.000_r8*rxt(k,279)*y(k,225) + .250_r8*rxt(k,386)*y(k,232) &
                      + .250_r8*rxt(k,391)*y(k,233) + .250_r8*rxt(k,353)*y(k,234) &
                      + .250_r8*rxt(k,465)*y(k,239) + .500_r8*rxt(k,341)*y(k,244) &
                      + .250_r8*rxt(k,470)*y(k,245) + .250_r8*rxt(k,475)*y(k,246) &
                      + .300_r8*rxt(k,411)*y(k,248)
         mat(k,1142) = .250_r8*rxt(k,386)*y(k,225)
         mat(k,1215) = .250_r8*rxt(k,391)*y(k,225)
         mat(k,1189) = .250_r8*rxt(k,353)*y(k,225)
         mat(k,1008) = .250_r8*rxt(k,465)*y(k,225)
         mat(k,1027) = .500_r8*rxt(k,341)*y(k,225)
         mat(k,989) = .250_r8*rxt(k,470)*y(k,225)
         mat(k,967) = .250_r8*rxt(k,475)*y(k,225)
         mat(k,1063) = .300_r8*rxt(k,411)*y(k,225)
         mat(k,352) = -(rxt(k,283)*y(k,241))
         mat(k,1609) = -rxt(k,283)*y(k,52)
         mat(k,1833) = rxt(k,280)*y(k,230)
         mat(k,1446) = rxt(k,280)*y(k,225)
         mat(k,2052) = -(rxt(k,194)*y(k,41) + rxt(k,196)*y(k,76) + rxt(k,197)*y(k,78) &
                      + (rxt(k,198) + rxt(k,199)) * y(k,230) + rxt(k,200)*y(k,132) &
                      + rxt(k,207)*y(k,59) + rxt(k,216)*y(k,91) + rxt(k,306)*y(k,27))
         mat(k,1730) = -rxt(k,194)*y(k,55)
         mat(k,1060) = -rxt(k,196)*y(k,55)
         mat(k,525) = -rxt(k,197)*y(k,55)
         mat(k,1524) = -(rxt(k,198) + rxt(k,199)) * y(k,55)
         mat(k,1938) = -rxt(k,200)*y(k,55)
         mat(k,879) = -rxt(k,207)*y(k,55)
         mat(k,734) = -rxt(k,216)*y(k,55)
         mat(k,247) = -rxt(k,306)*y(k,55)
         mat(k,1962) = rxt(k,235)*y(k,58)
         mat(k,1988) = rxt(k,235)*y(k,18) + (4.000_r8*rxt(k,202)+2.000_r8*rxt(k,204)) &
                      *y(k,58) + rxt(k,206)*y(k,122) + rxt(k,211)*y(k,131) &
                      + rxt(k,487)*y(k,149) + rxt(k,201)*y(k,225) + rxt(k,212) &
                      *y(k,241)
         mat(k,143) = rxt(k,256)*y(k,240)
         mat(k,1430) = rxt(k,214)*y(k,131) + rxt(k,226)*y(k,240) + rxt(k,215)*y(k,241)
         mat(k,1408) = rxt(k,206)*y(k,58)
         mat(k,2018) = rxt(k,211)*y(k,58) + rxt(k,214)*y(k,84)
         mat(k,1138) = rxt(k,487)*y(k,58)
         mat(k,1878) = rxt(k,201)*y(k,58)
         mat(k,1548) = rxt(k,256)*y(k,64) + rxt(k,226)*y(k,84)
         mat(k,1707) = rxt(k,212)*y(k,58) + rxt(k,215)*y(k,84)
         mat(k,2020) = rxt(k,207)*y(k,59)
         mat(k,1965) = 2.000_r8*rxt(k,203)*y(k,58)
         mat(k,869) = rxt(k,207)*y(k,55) + (rxt(k,536)+rxt(k,541)+rxt(k,546))*y(k,84)
         mat(k,1410) = (rxt(k,536)+rxt(k,541)+rxt(k,546))*y(k,59) + (rxt(k,531) &
                       +rxt(k,537)+rxt(k,542))*y(k,91)
         mat(k,728) = (rxt(k,531)+rxt(k,537)+rxt(k,542))*y(k,84)
         mat(k,1964) = 2.000_r8*rxt(k,228)*y(k,58)
         mat(k,1986) = -(rxt(k,201)*y(k,225) + (4._r8*rxt(k,202) + 4._r8*rxt(k,203) &
                      + 4._r8*rxt(k,204) + 4._r8*rxt(k,228)) * y(k,58) + rxt(k,205) &
                      *y(k,230) + rxt(k,206)*y(k,122) + rxt(k,208)*y(k,123) + rxt(k,211) &
                      *y(k,131) + (rxt(k,212) + rxt(k,213)) * y(k,241) + (rxt(k,234) &
                      + rxt(k,235) + rxt(k,236)) * y(k,18) + rxt(k,487)*y(k,149))
         mat(k,1876) = -rxt(k,201)*y(k,58)
         mat(k,1522) = -rxt(k,205)*y(k,58)
         mat(k,1406) = -rxt(k,206)*y(k,58)
         mat(k,1769) = -rxt(k,208)*y(k,58)
         mat(k,2016) = -rxt(k,211)*y(k,58)
         mat(k,1705) = -(rxt(k,212) + rxt(k,213)) * y(k,58)
         mat(k,1960) = -(rxt(k,234) + rxt(k,235) + rxt(k,236)) * y(k,58)
         mat(k,1136) = -rxt(k,487)*y(k,58)
         mat(k,2050) = rxt(k,216)*y(k,91) + rxt(k,200)*y(k,132) + rxt(k,199)*y(k,230)
         mat(k,877) = rxt(k,209)*y(k,131)
         mat(k,1428) = rxt(k,227)*y(k,240)
         mat(k,732) = rxt(k,216)*y(k,55) + rxt(k,217)*y(k,131) + rxt(k,218)*y(k,241)
         mat(k,2016) = mat(k,2016) + rxt(k,209)*y(k,59) + rxt(k,217)*y(k,91)
         mat(k,1936) = rxt(k,200)*y(k,55)
         mat(k,263) = rxt(k,492)*y(k,149)
         mat(k,1136) = mat(k,1136) + rxt(k,492)*y(k,133)
         mat(k,1522) = mat(k,1522) + rxt(k,199)*y(k,55)
         mat(k,1546) = rxt(k,227)*y(k,84)
         mat(k,1705) = mat(k,1705) + rxt(k,218)*y(k,91)
         mat(k,871) = -(rxt(k,207)*y(k,55) + rxt(k,209)*y(k,131) + rxt(k,210)*y(k,241) &
                      + (rxt(k,536) + rxt(k,541) + rxt(k,546)) * y(k,84))
         mat(k,2030) = -rxt(k,207)*y(k,59)
         mat(k,1998) = -rxt(k,209)*y(k,59)
         mat(k,1665) = -rxt(k,210)*y(k,59)
         mat(k,1414) = -(rxt(k,536) + rxt(k,541) + rxt(k,546)) * y(k,59)
         mat(k,1970) = rxt(k,208)*y(k,123)
         mat(k,1748) = rxt(k,208)*y(k,58)
         mat(k,955) = -((rxt(k,285) + rxt(k,295)) * y(k,241))
         mat(k,1672) = -(rxt(k,285) + rxt(k,295)) * y(k,61)
         mat(k,836) = .230_r8*rxt(k,460)*y(k,132)
         mat(k,1293) = rxt(k,230)*y(k,41)
         mat(k,240) = .350_r8*rxt(k,297)*y(k,241)
         mat(k,490) = .630_r8*rxt(k,299)*y(k,132)
         mat(k,894) = .560_r8*rxt(k,328)*y(k,132)
         mat(k,1712) = rxt(k,230)*y(k,16) + rxt(k,194)*y(k,55) + rxt(k,275)*y(k,124) &
                      + rxt(k,276)*y(k,131) + rxt(k,277)*y(k,241)
         mat(k,1115) = rxt(k,334)*y(k,124) + rxt(k,335)*y(k,241)
         mat(k,2032) = rxt(k,194)*y(k,41)
         mat(k,793) = rxt(k,322)*y(k,241)
         mat(k,768) = .620_r8*rxt(k,405)*y(k,132)
         mat(k,1103) = .650_r8*rxt(k,358)*y(k,132)
         mat(k,810) = .230_r8*rxt(k,463)*y(k,132)
         mat(k,1170) = .560_r8*rxt(k,372)*y(k,132)
         mat(k,1375) = .170_r8*rxt(k,431)*y(k,226) + .220_r8*rxt(k,356)*y(k,234) &
                      + .400_r8*rxt(k,434)*y(k,235) + .350_r8*rxt(k,437)*y(k,237) &
                      + .225_r8*rxt(k,472)*y(k,245) + .250_r8*rxt(k,413)*y(k,248)
         mat(k,1794) = rxt(k,275)*y(k,41) + rxt(k,334)*y(k,48) + .220_r8*rxt(k,355) &
                      *y(k,234) + .500_r8*rxt(k,414)*y(k,248)
         mat(k,1999) = rxt(k,276)*y(k,41) + rxt(k,482)*y(k,134)
         mat(k,1905) = .230_r8*rxt(k,460)*y(k,5) + .630_r8*rxt(k,299)*y(k,24) &
                      + .560_r8*rxt(k,328)*y(k,28) + .620_r8*rxt(k,405)*y(k,97) &
                      + .650_r8*rxt(k,358)*y(k,105) + .230_r8*rxt(k,463)*y(k,110) &
                      + .560_r8*rxt(k,372)*y(k,111)
         mat(k,305) = rxt(k,482)*y(k,131) + rxt(k,483)*y(k,241)
         mat(k,947) = .700_r8*rxt(k,481)*y(k,241)
         mat(k,1269) = .220_r8*rxt(k,352)*y(k,234) + .250_r8*rxt(k,410)*y(k,248)
         mat(k,1846) = .110_r8*rxt(k,353)*y(k,234) + .125_r8*rxt(k,470)*y(k,245) &
                      + .200_r8*rxt(k,411)*y(k,248)
         mat(k,666) = .170_r8*rxt(k,431)*y(k,122) + .070_r8*rxt(k,430)*y(k,230)
         mat(k,1490) = .070_r8*rxt(k,430)*y(k,226) + .160_r8*rxt(k,433)*y(k,235) &
                      + .140_r8*rxt(k,436)*y(k,237)
         mat(k,1192) = .220_r8*rxt(k,356)*y(k,122) + .220_r8*rxt(k,355)*y(k,124) &
                      + .220_r8*rxt(k,352)*y(k,224) + .110_r8*rxt(k,353)*y(k,225)
         mat(k,629) = .400_r8*rxt(k,434)*y(k,122) + .160_r8*rxt(k,433)*y(k,230)
         mat(k,784) = .350_r8*rxt(k,437)*y(k,122) + .140_r8*rxt(k,436)*y(k,230)
         mat(k,1672) = mat(k,1672) + .350_r8*rxt(k,297)*y(k,23) + rxt(k,277)*y(k,41) &
                      + rxt(k,335)*y(k,48) + rxt(k,322)*y(k,74) + rxt(k,483)*y(k,134) &
                      + .700_r8*rxt(k,481)*y(k,207)
         mat(k,993) = .225_r8*rxt(k,472)*y(k,122) + .125_r8*rxt(k,470)*y(k,225)
         mat(k,1066) = .250_r8*rxt(k,413)*y(k,122) + .500_r8*rxt(k,414)*y(k,124) &
                      + .250_r8*rxt(k,410)*y(k,224) + .200_r8*rxt(k,411)*y(k,225)
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
         mat(k,109) = -(rxt(k,255)*y(k,240))
         mat(k,1526) = -rxt(k,255)*y(k,63)
         mat(k,140) = -(rxt(k,256)*y(k,240))
         mat(k,1527) = -rxt(k,256)*y(k,64)
         mat(k,160) = -(rxt(k,429)*y(k,241))
         mat(k,1580) = -rxt(k,429)*y(k,65)
         mat(k,154) = .180_r8*rxt(k,449)*y(k,241)
         mat(k,1580) = mat(k,1580) + .180_r8*rxt(k,449)*y(k,209)
         mat(k,249) = -(rxt(k,496)*y(k,124) + (rxt(k,497) + rxt(k,499)) * y(k,241))
         mat(k,1776) = -rxt(k,496)*y(k,66)
         mat(k,1594) = -(rxt(k,497) + rxt(k,499)) * y(k,66)
         mat(k,618) = rxt(k,317)*y(k,230)
         mat(k,1432) = rxt(k,317)*y(k,229)
         mat(k,693) = -(rxt(k,252)*y(k,76) + rxt(k,253)*y(k,251) + rxt(k,254)*y(k,88))
         mat(k,1050) = -rxt(k,252)*y(k,72)
         mat(k,2057) = -rxt(k,253)*y(k,72)
         mat(k,1305) = -rxt(k,254)*y(k,72)
         mat(k,110) = 2.000_r8*rxt(k,255)*y(k,240)
         mat(k,141) = rxt(k,256)*y(k,240)
         mat(k,1529) = 2.000_r8*rxt(k,255)*y(k,63) + rxt(k,256)*y(k,64)
         mat(k,924) = -(rxt(k,321)*y(k,241))
         mat(k,1669) = -rxt(k,321)*y(k,73)
         mat(k,512) = .700_r8*rxt(k,396)*y(k,241)
         mat(k,473) = .500_r8*rxt(k,397)*y(k,241)
         mat(k,312) = rxt(k,408)*y(k,241)
         mat(k,1372) = .050_r8*rxt(k,394)*y(k,233) + .530_r8*rxt(k,356)*y(k,234) &
                      + .225_r8*rxt(k,472)*y(k,245) + .250_r8*rxt(k,413)*y(k,248)
         mat(k,1791) = .050_r8*rxt(k,395)*y(k,233) + .530_r8*rxt(k,355)*y(k,234) &
                      + .250_r8*rxt(k,414)*y(k,248)
         mat(k,1268) = .530_r8*rxt(k,352)*y(k,234) + .250_r8*rxt(k,410)*y(k,248)
         mat(k,1844) = .260_r8*rxt(k,353)*y(k,234) + .125_r8*rxt(k,470)*y(k,245) &
                      + .100_r8*rxt(k,411)*y(k,248)
         mat(k,1219) = .050_r8*rxt(k,394)*y(k,122) + .050_r8*rxt(k,395)*y(k,124)
         mat(k,1190) = .530_r8*rxt(k,356)*y(k,122) + .530_r8*rxt(k,355)*y(k,124) &
                      + .530_r8*rxt(k,352)*y(k,224) + .260_r8*rxt(k,353)*y(k,225)
         mat(k,1669) = mat(k,1669) + .700_r8*rxt(k,396)*y(k,98) + .500_r8*rxt(k,397) &
                      *y(k,99) + rxt(k,408)*y(k,115)
         mat(k,991) = .225_r8*rxt(k,472)*y(k,122) + .125_r8*rxt(k,470)*y(k,225)
         mat(k,1065) = .250_r8*rxt(k,413)*y(k,122) + .250_r8*rxt(k,414)*y(k,124) &
                      + .250_r8*rxt(k,410)*y(k,224) + .100_r8*rxt(k,411)*y(k,225)
         mat(k,792) = -(rxt(k,322)*y(k,241))
         mat(k,1660) = -rxt(k,322)*y(k,74)
         mat(k,239) = .650_r8*rxt(k,297)*y(k,241)
         mat(k,923) = .200_r8*rxt(k,321)*y(k,241)
         mat(k,911) = rxt(k,409)*y(k,241)
         mat(k,1367) = rxt(k,420)*y(k,219) + .050_r8*rxt(k,394)*y(k,233) &
                      + .400_r8*rxt(k,434)*y(k,235) + .170_r8*rxt(k,437)*y(k,237) &
                      + .700_r8*rxt(k,440)*y(k,242) + .600_r8*rxt(k,447)*y(k,247) &
                      + .250_r8*rxt(k,413)*y(k,248) + .340_r8*rxt(k,453)*y(k,249) &
                      + .170_r8*rxt(k,456)*y(k,250)
         mat(k,1783) = .050_r8*rxt(k,395)*y(k,233) + .250_r8*rxt(k,414)*y(k,248)
         mat(k,422) = rxt(k,420)*y(k,122)
         mat(k,1266) = .250_r8*rxt(k,410)*y(k,248)
         mat(k,1840) = .100_r8*rxt(k,411)*y(k,248)
         mat(k,1484) = .160_r8*rxt(k,433)*y(k,235) + .070_r8*rxt(k,436)*y(k,237)
         mat(k,1218) = .050_r8*rxt(k,394)*y(k,122) + .050_r8*rxt(k,395)*y(k,124)
         mat(k,628) = .400_r8*rxt(k,434)*y(k,122) + .160_r8*rxt(k,433)*y(k,230)
         mat(k,783) = .170_r8*rxt(k,437)*y(k,122) + .070_r8*rxt(k,436)*y(k,230)
         mat(k,1660) = mat(k,1660) + .650_r8*rxt(k,297)*y(k,23) + .200_r8*rxt(k,321) &
                      *y(k,73) + rxt(k,409)*y(k,116)
         mat(k,392) = .700_r8*rxt(k,440)*y(k,122)
         mat(k,641) = .600_r8*rxt(k,447)*y(k,122)
         mat(k,1064) = .250_r8*rxt(k,413)*y(k,122) + .250_r8*rxt(k,414)*y(k,124) &
                      + .250_r8*rxt(k,410)*y(k,224) + .100_r8*rxt(k,411)*y(k,225)
         mat(k,657) = .340_r8*rxt(k,453)*y(k,122)
         mat(k,429) = .170_r8*rxt(k,456)*y(k,122)
         mat(k,1320) = -((rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,230) + rxt(k,160) &
                      *y(k,132))
         mat(k,1510) = -(rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,75)
         mat(k,1924) = -rxt(k,160)*y(k,75)
         mat(k,1716) = rxt(k,277)*y(k,241)
         mat(k,2038) = rxt(k,196)*y(k,76)
         mat(k,956) = rxt(k,295)*y(k,241)
         mat(k,696) = rxt(k,252)*y(k,76)
         mat(k,1053) = rxt(k,196)*y(k,55) + rxt(k,252)*y(k,72) + rxt(k,152)*y(k,131) &
                      + rxt(k,144)*y(k,240) + rxt(k,161)*y(k,241)
         mat(k,722) = rxt(k,250)*y(k,240)
         mat(k,1417) = rxt(k,227)*y(k,240)
         mat(k,365) = rxt(k,182)*y(k,241)
         mat(k,2004) = rxt(k,152)*y(k,76) + rxt(k,164)*y(k,241)
         mat(k,307) = rxt(k,483)*y(k,241)
         mat(k,443) = rxt(k,488)*y(k,241)
         mat(k,1129) = rxt(k,493)*y(k,241)
         mat(k,1534) = rxt(k,144)*y(k,76) + rxt(k,250)*y(k,80) + rxt(k,227)*y(k,84)
         mat(k,1693) = rxt(k,277)*y(k,41) + rxt(k,295)*y(k,61) + rxt(k,161)*y(k,76) &
                      + rxt(k,182)*y(k,112) + rxt(k,164)*y(k,131) + rxt(k,483) &
                      *y(k,134) + rxt(k,488)*y(k,147) + rxt(k,493)*y(k,149)
         mat(k,1051) = -(rxt(k,144)*y(k,240) + rxt(k,152)*y(k,131) + rxt(k,161) &
                      *y(k,241) + rxt(k,196)*y(k,55) + rxt(k,252)*y(k,72))
         mat(k,1531) = -rxt(k,144)*y(k,76)
         mat(k,2000) = -rxt(k,152)*y(k,76)
         mat(k,1679) = -rxt(k,161)*y(k,76)
         mat(k,2034) = -rxt(k,196)*y(k,76)
         mat(k,694) = -rxt(k,252)*y(k,76)
         mat(k,1318) = rxt(k,154)*y(k,230)
         mat(k,1497) = rxt(k,154)*y(k,75)
         mat(k,520) = -(rxt(k,153)*y(k,131) + rxt(k,162)*y(k,241) + rxt(k,197)*y(k,55))
         mat(k,1993) = -rxt(k,153)*y(k,78)
         mat(k,1632) = -rxt(k,162)*y(k,78)
         mat(k,2024) = -rxt(k,197)*y(k,78)
         mat(k,1461) = 2.000_r8*rxt(k,168)*y(k,230)
         mat(k,1632) = mat(k,1632) + 2.000_r8*rxt(k,167)*y(k,241)
         mat(k,211) = rxt(k,495)*y(k,251)
         mat(k,2054) = rxt(k,495)*y(k,151)
         mat(k,720) = -(rxt(k,243)*y(k,131) + rxt(k,244)*y(k,241) + (rxt(k,249) &
                      + rxt(k,250)) * y(k,240))
         mat(k,1995) = -rxt(k,243)*y(k,80)
         mat(k,1653) = -rxt(k,244)*y(k,80)
         mat(k,1530) = -(rxt(k,249) + rxt(k,250)) * y(k,80)
         mat(k,1292) = rxt(k,230)*y(k,41) + rxt(k,231)*y(k,230)
         mat(k,1711) = rxt(k,230)*y(k,16)
         mat(k,1479) = rxt(k,231)*y(k,16)
         mat(k,1418) = -(rxt(k,214)*y(k,131) + rxt(k,215)*y(k,241) + (rxt(k,226) &
                      + rxt(k,227)) * y(k,240) + (rxt(k,531) + rxt(k,537) + rxt(k,542) &
                      ) * y(k,91) + (rxt(k,536) + rxt(k,541) + rxt(k,546)) * y(k,59) &
                      + (rxt(k,538) + rxt(k,543)) * y(k,90))
         mat(k,2006) = -rxt(k,214)*y(k,84)
         mat(k,1695) = -rxt(k,215)*y(k,84)
         mat(k,1536) = -(rxt(k,226) + rxt(k,227)) * y(k,84)
         mat(k,730) = -(rxt(k,531) + rxt(k,537) + rxt(k,542)) * y(k,84)
         mat(k,873) = -(rxt(k,536) + rxt(k,541) + rxt(k,546)) * y(k,84)
         mat(k,686) = -(rxt(k,538) + rxt(k,543)) * y(k,84)
         mat(k,245) = rxt(k,306)*y(k,55)
         mat(k,1718) = rxt(k,194)*y(k,55)
         mat(k,2040) = rxt(k,306)*y(k,27) + rxt(k,194)*y(k,41) + rxt(k,196)*y(k,76) &
                      + rxt(k,197)*y(k,78) + rxt(k,216)*y(k,91) + rxt(k,198)*y(k,230)
         mat(k,1976) = rxt(k,213)*y(k,241)
         mat(k,1054) = rxt(k,196)*y(k,55)
         mat(k,521) = rxt(k,197)*y(k,55)
         mat(k,730) = mat(k,730) + rxt(k,216)*y(k,55)
         mat(k,1512) = rxt(k,198)*y(k,55)
         mat(k,1695) = mat(k,1695) + rxt(k,213)*y(k,58)
         mat(k,144) = -(rxt(k,286)*y(k,241) + rxt(k,294)*y(k,240))
         mat(k,1578) = -rxt(k,286)*y(k,85)
         mat(k,1528) = -rxt(k,294)*y(k,85)
         mat(k,705) = -(rxt(k,287)*y(k,241))
         mat(k,1651) = -rxt(k,287)*y(k,86)
         mat(k,829) = .050_r8*rxt(k,460)*y(k,132)
         mat(k,238) = .350_r8*rxt(k,297)*y(k,241)
         mat(k,489) = .370_r8*rxt(k,299)*y(k,132)
         mat(k,891) = .120_r8*rxt(k,328)*y(k,132)
         mat(k,765) = .110_r8*rxt(k,405)*y(k,132)
         mat(k,1102) = .330_r8*rxt(k,358)*y(k,132)
         mat(k,803) = .050_r8*rxt(k,463)*y(k,132)
         mat(k,1168) = .120_r8*rxt(k,372)*y(k,132)
         mat(k,1362) = rxt(k,290)*y(k,231)
         mat(k,1893) = .050_r8*rxt(k,460)*y(k,5) + .370_r8*rxt(k,299)*y(k,24) &
                      + .120_r8*rxt(k,328)*y(k,28) + .110_r8*rxt(k,405)*y(k,97) &
                      + .330_r8*rxt(k,358)*y(k,105) + .050_r8*rxt(k,463)*y(k,110) &
                      + .120_r8*rxt(k,372)*y(k,111)
         mat(k,1477) = rxt(k,288)*y(k,231)
         mat(k,385) = rxt(k,290)*y(k,122) + rxt(k,288)*y(k,230)
         mat(k,1651) = mat(k,1651) + .350_r8*rxt(k,297)*y(k,23)
         mat(k,692) = rxt(k,252)*y(k,76) + rxt(k,254)*y(k,88) + rxt(k,253)*y(k,251)
         mat(k,1049) = rxt(k,252)*y(k,72)
         mat(k,1304) = rxt(k,254)*y(k,72)
         mat(k,2055) = rxt(k,253)*y(k,72)
         mat(k,1307) = -(rxt(k,191)*y(k,241) + rxt(k,254)*y(k,72))
         mat(k,1692) = -rxt(k,191)*y(k,88)
         mat(k,695) = -rxt(k,254)*y(k,88)
         mat(k,1715) = rxt(k,275)*y(k,124)
         mat(k,1042) = rxt(k,308)*y(k,124)
         mat(k,1118) = rxt(k,334)*y(k,124)
         mat(k,872) = (rxt(k,536)+rxt(k,541)+rxt(k,546))*y(k,84)
         mat(k,251) = rxt(k,496)*y(k,124)
         mat(k,1416) = (rxt(k,536)+rxt(k,541)+rxt(k,546))*y(k,59)
         mat(k,1756) = rxt(k,190)*y(k,241)
         mat(k,1813) = rxt(k,275)*y(k,41) + rxt(k,308)*y(k,44) + rxt(k,334)*y(k,48) &
                      + rxt(k,496)*y(k,66)
         mat(k,1692) = mat(k,1692) + rxt(k,190)*y(k,123)
         mat(k,358) = -(rxt(k,169)*y(k,241))
         mat(k,1610) = -rxt(k,169)*y(k,89)
         mat(k,1734) = rxt(k,188)*y(k,230)
         mat(k,1447) = rxt(k,188)*y(k,123)
         mat(k,684) = -(rxt(k,245)*y(k,131) + (rxt(k,538) + rxt(k,543)) * y(k,84))
         mat(k,1994) = -rxt(k,245)*y(k,90)
         mat(k,1412) = -(rxt(k,538) + rxt(k,543)) * y(k,90)
         mat(k,1943) = rxt(k,237)*y(k,230)
         mat(k,1476) = rxt(k,237)*y(k,18)
         mat(k,729) = -(rxt(k,216)*y(k,55) + rxt(k,217)*y(k,131) + rxt(k,218)*y(k,241) &
                      + (rxt(k,531) + rxt(k,537) + rxt(k,542)) * y(k,84))
         mat(k,2027) = -rxt(k,216)*y(k,91)
         mat(k,1996) = -rxt(k,217)*y(k,91)
         mat(k,1654) = -rxt(k,218)*y(k,91)
         mat(k,1413) = -(rxt(k,531) + rxt(k,537) + rxt(k,542)) * y(k,91)
         mat(k,1968) = rxt(k,205)*y(k,230)
         mat(k,870) = rxt(k,210)*y(k,241)
         mat(k,1480) = rxt(k,205)*y(k,58)
         mat(k,1654) = mat(k,1654) + rxt(k,210)*y(k,59)
         mat(k,932) = -(rxt(k,351)*y(k,241))
         mat(k,1670) = -rxt(k,351)*y(k,92)
         mat(k,513) = .300_r8*rxt(k,396)*y(k,241)
         mat(k,474) = .500_r8*rxt(k,397)*y(k,241)
         mat(k,1373) = rxt(k,350)*y(k,227) + rxt(k,357)*y(k,234)
         mat(k,497) = rxt(k,350)*y(k,122)
         mat(k,1191) = rxt(k,357)*y(k,122)
         mat(k,1670) = mat(k,1670) + .300_r8*rxt(k,396)*y(k,98) + .500_r8*rxt(k,397) &
                      *y(k,99)
         mat(k,203) = -(rxt(k,382)*y(k,241))
         mat(k,1588) = -rxt(k,382)*y(k,93)
         mat(k,959) = -(rxt(k,336)*y(k,241))
         mat(k,1673) = -rxt(k,336)*y(k,94)
         mat(k,514) = .700_r8*rxt(k,396)*y(k,241)
         mat(k,475) = .500_r8*rxt(k,397)*y(k,241)
         mat(k,480) = .500_r8*rxt(k,371)*y(k,241)
         mat(k,1376) = .050_r8*rxt(k,394)*y(k,233) + .220_r8*rxt(k,356)*y(k,234) &
                      + .250_r8*rxt(k,413)*y(k,248)
         mat(k,1795) = .050_r8*rxt(k,395)*y(k,233) + .220_r8*rxt(k,355)*y(k,234) &
                      + .250_r8*rxt(k,414)*y(k,248)
         mat(k,461) = .500_r8*rxt(k,340)*y(k,241)
         mat(k,1270) = .220_r8*rxt(k,352)*y(k,234) + .250_r8*rxt(k,410)*y(k,248)
         mat(k,1847) = .230_r8*rxt(k,353)*y(k,234) + .200_r8*rxt(k,341)*y(k,244) &
                      + .100_r8*rxt(k,411)*y(k,248)
         mat(k,1222) = .050_r8*rxt(k,394)*y(k,122) + .050_r8*rxt(k,395)*y(k,124)
         mat(k,1193) = .220_r8*rxt(k,356)*y(k,122) + .220_r8*rxt(k,355)*y(k,124) &
                      + .220_r8*rxt(k,352)*y(k,224) + .230_r8*rxt(k,353)*y(k,225)
         mat(k,1673) = mat(k,1673) + .700_r8*rxt(k,396)*y(k,98) + .500_r8*rxt(k,397) &
                      *y(k,99) + .500_r8*rxt(k,371)*y(k,109) + .500_r8*rxt(k,340) &
                      *y(k,145)
         mat(k,1028) = .200_r8*rxt(k,341)*y(k,225)
         mat(k,1067) = .250_r8*rxt(k,413)*y(k,122) + .250_r8*rxt(k,414)*y(k,124) &
                      + .250_r8*rxt(k,410)*y(k,224) + .100_r8*rxt(k,411)*y(k,225)
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
         mat(k,281) = -(rxt(k,383)*y(k,241))
         mat(k,1598) = -rxt(k,383)*y(k,95)
         mat(k,1337) = .870_r8*rxt(k,394)*y(k,233)
         mat(k,1777) = .950_r8*rxt(k,395)*y(k,233)
         mat(k,1262) = rxt(k,390)*y(k,233)
         mat(k,1831) = .750_r8*rxt(k,391)*y(k,233)
         mat(k,1211) = .870_r8*rxt(k,394)*y(k,122) + .950_r8*rxt(k,395)*y(k,124) &
                      + rxt(k,390)*y(k,224) + .750_r8*rxt(k,391)*y(k,225)
         mat(k,103) = -(rxt(k,384)*y(k,241))
         mat(k,1572) = -rxt(k,384)*y(k,96)
         mat(k,595) = .600_r8*rxt(k,407)*y(k,241)
         mat(k,1572) = mat(k,1572) + .600_r8*rxt(k,407)*y(k,102)
         mat(k,766) = -(rxt(k,398)*y(k,124) + rxt(k,405)*y(k,132) + rxt(k,406) &
                      *y(k,241))
         mat(k,1782) = -rxt(k,398)*y(k,97)
         mat(k,1895) = -rxt(k,405)*y(k,97)
         mat(k,1658) = -rxt(k,406)*y(k,97)
         mat(k,511) = -(rxt(k,396)*y(k,241))
         mat(k,1631) = -rxt(k,396)*y(k,98)
         mat(k,1350) = .080_r8*rxt(k,388)*y(k,232)
         mat(k,1140) = .080_r8*rxt(k,388)*y(k,122)
         mat(k,471) = -(rxt(k,397)*y(k,241))
         mat(k,1627) = -rxt(k,397)*y(k,99)
         mat(k,1348) = .080_r8*rxt(k,394)*y(k,233)
         mat(k,1212) = .080_r8*rxt(k,394)*y(k,122)
         mat(k,1088) = -(rxt(k,399)*y(k,224) + rxt(k,400)*y(k,225) + rxt(k,401) &
                      *y(k,230) + rxt(k,402)*y(k,122) + rxt(k,403)*y(k,124))
         mat(k,1272) = -rxt(k,399)*y(k,100)
         mat(k,1854) = -rxt(k,400)*y(k,100)
         mat(k,1499) = -rxt(k,401)*y(k,100)
         mat(k,1383) = -rxt(k,402)*y(k,100)
         mat(k,1802) = -rxt(k,403)*y(k,100)
         mat(k,769) = rxt(k,398)*y(k,124)
         mat(k,1802) = mat(k,1802) + rxt(k,398)*y(k,97)
         mat(k,322) = -(rxt(k,404)*y(k,241))
         mat(k,1605) = -rxt(k,404)*y(k,101)
         mat(k,1080) = rxt(k,401)*y(k,230)
         mat(k,1442) = rxt(k,401)*y(k,100)
         mat(k,596) = -(rxt(k,407)*y(k,241))
         mat(k,1641) = -rxt(k,407)*y(k,102)
         mat(k,1468) = rxt(k,387)*y(k,232) + rxt(k,392)*y(k,233)
         mat(k,1141) = rxt(k,387)*y(k,230)
         mat(k,1214) = rxt(k,392)*y(k,230)
         mat(k,61) = -(rxt(k,521)*y(k,241))
         mat(k,1555) = -rxt(k,521)*y(k,103)
         mat(k,77) = -(rxt(k,522)*y(k,241))
         mat(k,1566) = -rxt(k,522)*y(k,104)
         mat(k,1104) = -(rxt(k,358)*y(k,132) + rxt(k,359)*y(k,241))
         mat(k,1913) = -rxt(k,358)*y(k,105)
         mat(k,1682) = -rxt(k,359)*y(k,105)
         mat(k,770) = .300_r8*rxt(k,405)*y(k,132)
         mat(k,1384) = .360_r8*rxt(k,388)*y(k,232)
         mat(k,1803) = .400_r8*rxt(k,389)*y(k,232)
         mat(k,1913) = mat(k,1913) + .300_r8*rxt(k,405)*y(k,97)
         mat(k,1273) = .390_r8*rxt(k,385)*y(k,232)
         mat(k,1855) = .310_r8*rxt(k,386)*y(k,232)
         mat(k,1150) = .360_r8*rxt(k,388)*y(k,122) + .400_r8*rxt(k,389)*y(k,124) &
                      + .390_r8*rxt(k,385)*y(k,224) + .310_r8*rxt(k,386)*y(k,225)
         mat(k,289) = -(rxt(k,360)*y(k,241))
         mat(k,1600) = -rxt(k,360)*y(k,106)
         mat(k,1440) = rxt(k,354)*y(k,234)
         mat(k,1188) = rxt(k,354)*y(k,230)
         mat(k,447) = -(rxt(k,369)*y(k,241))
         mat(k,1623) = -rxt(k,369)*y(k,107)
         mat(k,1346) = .800_r8*rxt(k,378)*y(k,218)
         mat(k,848) = .800_r8*rxt(k,378)*y(k,122)
         mat(k,284) = -(rxt(k,370)*y(k,241))
         mat(k,1599) = -rxt(k,370)*y(k,108)
         mat(k,1439) = .800_r8*rxt(k,367)*y(k,238)
         mat(k,587) = .800_r8*rxt(k,367)*y(k,230)
         mat(k,479) = -(rxt(k,371)*y(k,241))
         mat(k,1628) = -rxt(k,371)*y(k,109)
         mat(k,1739) = rxt(k,374)*y(k,236)
         mat(k,1246) = rxt(k,374)*y(k,123)
         mat(k,805) = -(rxt(k,462)*y(k,124) + rxt(k,463)*y(k,132) + rxt(k,464) &
                      *y(k,241))
         mat(k,1784) = -rxt(k,462)*y(k,110)
         mat(k,1897) = -rxt(k,463)*y(k,110)
         mat(k,1661) = -rxt(k,464)*y(k,110)
         mat(k,1174) = -(rxt(k,372)*y(k,132) + rxt(k,373)*y(k,241))
         mat(k,1917) = -rxt(k,372)*y(k,111)
         mat(k,1686) = -rxt(k,373)*y(k,111)
         mat(k,772) = .200_r8*rxt(k,405)*y(k,132)
         mat(k,1387) = .560_r8*rxt(k,388)*y(k,232)
         mat(k,1807) = .600_r8*rxt(k,389)*y(k,232)
         mat(k,1917) = mat(k,1917) + .200_r8*rxt(k,405)*y(k,97)
         mat(k,1276) = .610_r8*rxt(k,385)*y(k,232)
         mat(k,1858) = .440_r8*rxt(k,386)*y(k,232)
         mat(k,1153) = .560_r8*rxt(k,388)*y(k,122) + .600_r8*rxt(k,389)*y(k,124) &
                      + .610_r8*rxt(k,385)*y(k,224) + .440_r8*rxt(k,386)*y(k,225)
         mat(k,364) = -(rxt(k,170)*y(k,122) + (rxt(k,171) + rxt(k,172) + rxt(k,173) &
                      ) * y(k,123) + rxt(k,182)*y(k,241))
         mat(k,1338) = -rxt(k,170)*y(k,112)
         mat(k,1735) = -(rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,112)
         mat(k,1611) = -rxt(k,182)*y(k,112)
         mat(k,1733) = rxt(k,189)*y(k,124)
         mat(k,1775) = rxt(k,189)*y(k,123)
         mat(k,310) = -(rxt(k,408)*y(k,241))
         mat(k,1603) = -rxt(k,408)*y(k,115)
         mat(k,1079) = .200_r8*rxt(k,400)*y(k,225)
         mat(k,1832) = .200_r8*rxt(k,400)*y(k,100)
         mat(k,912) = -(rxt(k,409)*y(k,241))
         mat(k,1668) = -rxt(k,409)*y(k,116)
         mat(k,1084) = rxt(k,402)*y(k,122) + rxt(k,403)*y(k,124) + rxt(k,399)*y(k,224) &
                      + .800_r8*rxt(k,400)*y(k,225)
         mat(k,1371) = rxt(k,402)*y(k,100)
         mat(k,1790) = rxt(k,403)*y(k,100)
         mat(k,1267) = rxt(k,399)*y(k,100)
         mat(k,1843) = .800_r8*rxt(k,400)*y(k,100)
         mat(k,90) = -(rxt(k,498)*y(k,241))
         mat(k,1569) = -rxt(k,498)*y(k,120)
         mat(k,1395) = -(rxt(k,170)*y(k,112) + rxt(k,179)*y(k,124) + rxt(k,183) &
                      *y(k,230) + rxt(k,184)*y(k,132) + rxt(k,185)*y(k,131) + rxt(k,206) &
                      *y(k,58) + rxt(k,238)*y(k,18) + rxt(k,281)*y(k,225) + rxt(k,290) &
                      *y(k,231) + rxt(k,303)*y(k,221) + rxt(k,314)*y(k,224) + rxt(k,318) &
                      *y(k,229) + rxt(k,331)*y(k,222) + rxt(k,339)*y(k,243) + rxt(k,343) &
                      *y(k,244) + (rxt(k,349) + rxt(k,350)) * y(k,227) + (rxt(k,356) &
                      + rxt(k,357)) * y(k,234) + rxt(k,365)*y(k,236) + rxt(k,368) &
                      *y(k,238) + (rxt(k,378) + rxt(k,379)) * y(k,218) + rxt(k,388) &
                      *y(k,232) + rxt(k,394)*y(k,233) + rxt(k,402)*y(k,100) + rxt(k,413) &
                      *y(k,248) + rxt(k,417)*y(k,217) + rxt(k,420)*y(k,219) + rxt(k,425) &
                      *y(k,220) + rxt(k,427)*y(k,223) + rxt(k,431)*y(k,226) + rxt(k,434) &
                      *y(k,235) + rxt(k,437)*y(k,237) + rxt(k,440)*y(k,242) + rxt(k,447) &
                      *y(k,247) + rxt(k,453)*y(k,249) + rxt(k,456)*y(k,250) + rxt(k,467) &
                      *y(k,239) + rxt(k,472)*y(k,245) + rxt(k,477)*y(k,246))
         mat(k,366) = -rxt(k,170)*y(k,122)
         mat(k,1815) = -rxt(k,179)*y(k,122)
         mat(k,1511) = -rxt(k,183)*y(k,122)
         mat(k,1925) = -rxt(k,184)*y(k,122)
         mat(k,2005) = -rxt(k,185)*y(k,122)
         mat(k,1975) = -rxt(k,206)*y(k,122)
         mat(k,1949) = -rxt(k,238)*y(k,122)
         mat(k,1865) = -rxt(k,281)*y(k,122)
         mat(k,386) = -rxt(k,290)*y(k,122)
         mat(k,741) = -rxt(k,303)*y(k,122)
         mat(k,1283) = -rxt(k,314)*y(k,122)
         mat(k,622) = -rxt(k,318)*y(k,122)
         mat(k,713) = -rxt(k,331)*y(k,122)
         mat(k,677) = -rxt(k,339)*y(k,122)
         mat(k,1032) = -rxt(k,343)*y(k,122)
         mat(k,499) = -(rxt(k,349) + rxt(k,350)) * y(k,122)
         mat(k,1202) = -(rxt(k,356) + rxt(k,357)) * y(k,122)
         mat(k,1253) = -rxt(k,365)*y(k,122)
         mat(k,591) = -rxt(k,368)*y(k,122)
         mat(k,857) = -(rxt(k,378) + rxt(k,379)) * y(k,122)
         mat(k,1159) = -rxt(k,388)*y(k,122)
         mat(k,1235) = -rxt(k,394)*y(k,122)
         mat(k,1095) = -rxt(k,402)*y(k,122)
         mat(k,1072) = -rxt(k,413)*y(k,122)
         mat(k,455) = -rxt(k,417)*y(k,122)
         mat(k,423) = -rxt(k,420)*y(k,122)
         mat(k,380) = -rxt(k,425)*y(k,122)
         mat(k,560) = -rxt(k,427)*y(k,122)
         mat(k,668) = -rxt(k,431)*y(k,122)
         mat(k,630) = -rxt(k,434)*y(k,122)
         mat(k,786) = -rxt(k,437)*y(k,122)
         mat(k,393) = -rxt(k,440)*y(k,122)
         mat(k,644) = -rxt(k,447)*y(k,122)
         mat(k,661) = -rxt(k,453)*y(k,122)
         mat(k,431) = -rxt(k,456)*y(k,122)
         mat(k,1018) = -rxt(k,467)*y(k,122)
         mat(k,999) = -rxt(k,472)*y(k,122)
         mat(k,979) = -rxt(k,477)*y(k,122)
         mat(k,366) = mat(k,366) + 2.000_r8*rxt(k,172)*y(k,123) + rxt(k,182)*y(k,241)
         mat(k,1758) = 2.000_r8*rxt(k,172)*y(k,112) + rxt(k,175)*y(k,131) + rxt(k,489) &
                      *y(k,149)
         mat(k,2005) = mat(k,2005) + rxt(k,175)*y(k,123)
         mat(k,1130) = rxt(k,489)*y(k,123)
         mat(k,1694) = rxt(k,182)*y(k,112)
         mat(k,1764) = -((rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,112) + (rxt(k,175) &
                      + rxt(k,177)) * y(k,131) + rxt(k,176)*y(k,132) + rxt(k,188) &
                      *y(k,230) + rxt(k,189)*y(k,124) + rxt(k,190)*y(k,241) + rxt(k,208) &
                      *y(k,58) + rxt(k,239)*y(k,18) + rxt(k,325)*y(k,224) + rxt(k,374) &
                      *y(k,236) + rxt(k,432)*y(k,226) + rxt(k,435)*y(k,235) + rxt(k,438) &
                      *y(k,237) + rxt(k,442)*y(k,138) + rxt(k,445)*y(k,217) + rxt(k,489) &
                      *y(k,149))
         mat(k,368) = -(rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,123)
         mat(k,2011) = -(rxt(k,175) + rxt(k,177)) * y(k,123)
         mat(k,1931) = -rxt(k,176)*y(k,123)
         mat(k,1517) = -rxt(k,188)*y(k,123)
         mat(k,1821) = -rxt(k,189)*y(k,123)
         mat(k,1700) = -rxt(k,190)*y(k,123)
         mat(k,1981) = -rxt(k,208)*y(k,123)
         mat(k,1955) = -rxt(k,239)*y(k,123)
         mat(k,1287) = -rxt(k,325)*y(k,123)
         mat(k,1257) = -rxt(k,374)*y(k,123)
         mat(k,671) = -rxt(k,432)*y(k,123)
         mat(k,632) = -rxt(k,435)*y(k,123)
         mat(k,789) = -rxt(k,438)*y(k,123)
         mat(k,402) = -rxt(k,442)*y(k,123)
         mat(k,458) = -rxt(k,445)*y(k,123)
         mat(k,1133) = -rxt(k,489)*y(k,123)
         mat(k,575) = rxt(k,376)*y(k,241)
         mat(k,300) = rxt(k,347)*y(k,124)
         mat(k,1955) = mat(k,1955) + rxt(k,238)*y(k,122)
         mat(k,1981) = mat(k,1981) + rxt(k,206)*y(k,122)
         mat(k,361) = rxt(k,169)*y(k,241)
         mat(k,519) = .700_r8*rxt(k,396)*y(k,241)
         mat(k,1099) = rxt(k,402)*y(k,122) + rxt(k,403)*y(k,124)
         mat(k,1401) = rxt(k,238)*y(k,18) + rxt(k,206)*y(k,58) + rxt(k,402)*y(k,100) &
                      + 2.000_r8*rxt(k,179)*y(k,124) + rxt(k,185)*y(k,131) &
                      + rxt(k,184)*y(k,132) + rxt(k,417)*y(k,217) + rxt(k,378) &
                      *y(k,218) + rxt(k,420)*y(k,219) + rxt(k,425)*y(k,220) &
                      + rxt(k,303)*y(k,221) + rxt(k,331)*y(k,222) + rxt(k,427) &
                      *y(k,223) + rxt(k,314)*y(k,224) + rxt(k,281)*y(k,225) &
                      + rxt(k,431)*y(k,226) + rxt(k,349)*y(k,227) + rxt(k,318) &
                      *y(k,229) + rxt(k,183)*y(k,230) + rxt(k,290)*y(k,231) &
                      + .920_r8*rxt(k,388)*y(k,232) + .920_r8*rxt(k,394)*y(k,233) &
                      + rxt(k,356)*y(k,234) + rxt(k,434)*y(k,235) + rxt(k,365) &
                      *y(k,236) + rxt(k,437)*y(k,237) + rxt(k,368)*y(k,238) &
                      + 1.600_r8*rxt(k,467)*y(k,239) + rxt(k,440)*y(k,242) &
                      + rxt(k,339)*y(k,243) + rxt(k,343)*y(k,244) + .900_r8*rxt(k,472) &
                      *y(k,245) + .800_r8*rxt(k,477)*y(k,246) + rxt(k,447)*y(k,247) &
                      + rxt(k,413)*y(k,248) + rxt(k,453)*y(k,249) + rxt(k,456) &
                      *y(k,250)
         mat(k,1821) = mat(k,1821) + rxt(k,347)*y(k,15) + rxt(k,403)*y(k,100) &
                      + 2.000_r8*rxt(k,179)*y(k,122) + rxt(k,180)*y(k,131) &
                      + rxt(k,178)*y(k,230) + rxt(k,389)*y(k,232) + rxt(k,395) &
                      *y(k,233) + rxt(k,355)*y(k,234) + rxt(k,366)*y(k,236) &
                      + 2.000_r8*rxt(k,468)*y(k,239) + rxt(k,181)*y(k,241) &
                      + rxt(k,414)*y(k,248)
         mat(k,752) = rxt(k,337)*y(k,241)
         mat(k,2011) = mat(k,2011) + rxt(k,185)*y(k,122) + rxt(k,180)*y(k,124)
         mat(k,1931) = mat(k,1931) + rxt(k,184)*y(k,122)
         mat(k,556) = rxt(k,474)*y(k,241)
         mat(k,458) = mat(k,458) + rxt(k,417)*y(k,122)
         mat(k,861) = rxt(k,378)*y(k,122)
         mat(k,426) = rxt(k,420)*y(k,122)
         mat(k,383) = rxt(k,425)*y(k,122)
         mat(k,745) = rxt(k,303)*y(k,122)
         mat(k,717) = rxt(k,331)*y(k,122)
         mat(k,563) = rxt(k,427)*y(k,122)
         mat(k,1287) = mat(k,1287) + rxt(k,314)*y(k,122)
         mat(k,1871) = rxt(k,281)*y(k,122) + .500_r8*rxt(k,465)*y(k,239)
         mat(k,671) = mat(k,671) + rxt(k,431)*y(k,122)
         mat(k,502) = rxt(k,349)*y(k,122)
         mat(k,626) = rxt(k,318)*y(k,122)
         mat(k,1517) = mat(k,1517) + rxt(k,183)*y(k,122) + rxt(k,178)*y(k,124)
         mat(k,389) = rxt(k,290)*y(k,122)
         mat(k,1163) = .920_r8*rxt(k,388)*y(k,122) + rxt(k,389)*y(k,124)
         mat(k,1239) = .920_r8*rxt(k,394)*y(k,122) + rxt(k,395)*y(k,124)
         mat(k,1206) = rxt(k,356)*y(k,122) + rxt(k,355)*y(k,124)
         mat(k,632) = mat(k,632) + rxt(k,434)*y(k,122)
         mat(k,1257) = mat(k,1257) + rxt(k,365)*y(k,122) + rxt(k,366)*y(k,124)
         mat(k,789) = mat(k,789) + rxt(k,437)*y(k,122)
         mat(k,594) = rxt(k,368)*y(k,122)
         mat(k,1022) = 1.600_r8*rxt(k,467)*y(k,122) + 2.000_r8*rxt(k,468)*y(k,124) &
                      + .500_r8*rxt(k,465)*y(k,225)
         mat(k,1700) = mat(k,1700) + rxt(k,376)*y(k,1) + rxt(k,169)*y(k,89) &
                      + .700_r8*rxt(k,396)*y(k,98) + rxt(k,181)*y(k,124) + rxt(k,337) &
                      *y(k,125) + rxt(k,474)*y(k,204)
         mat(k,396) = rxt(k,440)*y(k,122)
         mat(k,681) = rxt(k,339)*y(k,122)
         mat(k,1036) = rxt(k,343)*y(k,122)
         mat(k,1003) = .900_r8*rxt(k,472)*y(k,122)
         mat(k,983) = .800_r8*rxt(k,477)*y(k,122)
         mat(k,647) = rxt(k,447)*y(k,122)
         mat(k,1076) = rxt(k,413)*y(k,122) + rxt(k,414)*y(k,124)
         mat(k,664) = rxt(k,453)*y(k,122)
         mat(k,434) = rxt(k,456)*y(k,122)
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
         mat(k,1822) = -(rxt(k,178)*y(k,230) + rxt(k,179)*y(k,122) + rxt(k,180) &
                      *y(k,131) + rxt(k,181)*y(k,241) + rxt(k,189)*y(k,123) + rxt(k,275) &
                      *y(k,41) + rxt(k,308)*y(k,44) + rxt(k,327)*y(k,28) + rxt(k,334) &
                      *y(k,48) + rxt(k,347)*y(k,15) + rxt(k,355)*y(k,234) + rxt(k,366) &
                      *y(k,236) + rxt(k,389)*y(k,232) + rxt(k,395)*y(k,233) + rxt(k,398) &
                      *y(k,97) + rxt(k,403)*y(k,100) + rxt(k,414)*y(k,248) + rxt(k,459) &
                      *y(k,5) + rxt(k,462)*y(k,110) + rxt(k,468)*y(k,239) + rxt(k,479) &
                      *y(k,206) + rxt(k,496)*y(k,66))
         mat(k,1518) = -rxt(k,178)*y(k,124)
         mat(k,1402) = -rxt(k,179)*y(k,124)
         mat(k,2012) = -rxt(k,180)*y(k,124)
         mat(k,1701) = -rxt(k,181)*y(k,124)
         mat(k,1765) = -rxt(k,189)*y(k,124)
         mat(k,1724) = -rxt(k,275)*y(k,124)
         mat(k,1046) = -rxt(k,308)*y(k,124)
         mat(k,904) = -rxt(k,327)*y(k,124)
         mat(k,1122) = -rxt(k,334)*y(k,124)
         mat(k,301) = -rxt(k,347)*y(k,124)
         mat(k,1207) = -rxt(k,355)*y(k,124)
         mat(k,1258) = -rxt(k,366)*y(k,124)
         mat(k,1164) = -rxt(k,389)*y(k,124)
         mat(k,1240) = -rxt(k,395)*y(k,124)
         mat(k,778) = -rxt(k,398)*y(k,124)
         mat(k,1100) = -rxt(k,403)*y(k,124)
         mat(k,1077) = -rxt(k,414)*y(k,124)
         mat(k,845) = -rxt(k,459)*y(k,124)
         mat(k,819) = -rxt(k,462)*y(k,124)
         mat(k,1023) = -rxt(k,468)*y(k,124)
         mat(k,888) = -rxt(k,479)*y(k,124)
         mat(k,254) = -rxt(k,496)*y(k,124)
         mat(k,508) = rxt(k,240)*y(k,131)
         mat(k,2046) = rxt(k,207)*y(k,59)
         mat(k,876) = rxt(k,207)*y(k,55) + rxt(k,209)*y(k,131) + rxt(k,210)*y(k,241)
         mat(k,698) = rxt(k,254)*y(k,88)
         mat(k,1313) = rxt(k,254)*y(k,72) + rxt(k,191)*y(k,241)
         mat(k,486) = .500_r8*rxt(k,371)*y(k,241)
         mat(k,1765) = mat(k,1765) + rxt(k,177)*y(k,131) + rxt(k,176)*y(k,132)
         mat(k,2012) = mat(k,2012) + rxt(k,240)*y(k,19) + rxt(k,209)*y(k,59) &
                      + rxt(k,177)*y(k,123)
         mat(k,1932) = rxt(k,176)*y(k,123)
         mat(k,409) = rxt(k,323)*y(k,241)
         mat(k,1701) = mat(k,1701) + rxt(k,210)*y(k,59) + rxt(k,191)*y(k,88) &
                      + .500_r8*rxt(k,371)*y(k,109) + rxt(k,323)*y(k,136)
         mat(k,747) = -(rxt(k,337)*y(k,241))
         mat(k,1656) = -rxt(k,337)*y(k,125)
         mat(k,892) = rxt(k,327)*y(k,124)
         mat(k,472) = .500_r8*rxt(k,397)*y(k,241)
         mat(k,324) = rxt(k,404)*y(k,241)
         mat(k,311) = rxt(k,408)*y(k,241)
         mat(k,909) = rxt(k,409)*y(k,241)
         mat(k,1780) = rxt(k,327)*y(k,28)
         mat(k,1656) = mat(k,1656) + .500_r8*rxt(k,397)*y(k,99) + rxt(k,404)*y(k,101) &
                      + rxt(k,408)*y(k,115) + rxt(k,409)*y(k,116)
         mat(k,328) = -(rxt(k,469)*y(k,241))
         mat(k,1606) = -rxt(k,469)*y(k,126)
         mat(k,1443) = rxt(k,466)*y(k,239)
         mat(k,1006) = rxt(k,466)*y(k,230)
         mat(k,2017) = -(rxt(k,149)*y(k,132) + 4._r8*rxt(k,150)*y(k,131) + rxt(k,152) &
                      *y(k,76) + rxt(k,153)*y(k,78) + rxt(k,158)*y(k,230) + rxt(k,164) &
                      *y(k,241) + (rxt(k,175) + rxt(k,177)) * y(k,123) + rxt(k,180) &
                      *y(k,124) + rxt(k,185)*y(k,122) + rxt(k,209)*y(k,59) + rxt(k,211) &
                      *y(k,58) + rxt(k,214)*y(k,84) + rxt(k,217)*y(k,91) + rxt(k,240) &
                      *y(k,19) + rxt(k,241)*y(k,18) + rxt(k,243)*y(k,80) + rxt(k,245) &
                      *y(k,90) + rxt(k,276)*y(k,41) + rxt(k,482)*y(k,134))
         mat(k,1937) = -rxt(k,149)*y(k,131)
         mat(k,1059) = -rxt(k,152)*y(k,131)
         mat(k,524) = -rxt(k,153)*y(k,131)
         mat(k,1523) = -rxt(k,158)*y(k,131)
         mat(k,1706) = -rxt(k,164)*y(k,131)
         mat(k,1770) = -(rxt(k,175) + rxt(k,177)) * y(k,131)
         mat(k,1827) = -rxt(k,180)*y(k,131)
         mat(k,1407) = -rxt(k,185)*y(k,131)
         mat(k,878) = -rxt(k,209)*y(k,131)
         mat(k,1987) = -rxt(k,211)*y(k,131)
         mat(k,1429) = -rxt(k,214)*y(k,131)
         mat(k,733) = -rxt(k,217)*y(k,131)
         mat(k,510) = -rxt(k,240)*y(k,131)
         mat(k,1961) = -rxt(k,241)*y(k,131)
         mat(k,726) = -rxt(k,243)*y(k,131)
         mat(k,689) = -rxt(k,245)*y(k,131)
         mat(k,1729) = -rxt(k,276)*y(k,131)
         mat(k,309) = -rxt(k,482)*y(k,131)
         mat(k,1329) = rxt(k,156)*y(k,230)
         mat(k,369) = rxt(k,170)*y(k,122) + rxt(k,171)*y(k,123)
         mat(k,1407) = mat(k,1407) + rxt(k,170)*y(k,112)
         mat(k,1770) = mat(k,1770) + rxt(k,171)*y(k,112)
         mat(k,1523) = mat(k,1523) + rxt(k,156)*y(k,75)
         mat(k,1706) = mat(k,1706) + 2.000_r8*rxt(k,166)*y(k,241)
         mat(k,1934) = -(rxt(k,148)*y(k,240) + rxt(k,149)*y(k,131) + rxt(k,159) &
                      *y(k,230) + rxt(k,160)*y(k,75) + rxt(k,165)*y(k,241) + rxt(k,176) &
                      *y(k,123) + rxt(k,184)*y(k,122) + rxt(k,200)*y(k,55) + rxt(k,232) &
                      *y(k,16) + rxt(k,299)*y(k,24) + rxt(k,328)*y(k,28) + rxt(k,358) &
                      *y(k,105) + rxt(k,372)*y(k,111) + rxt(k,405)*y(k,97) + rxt(k,443) &
                      *y(k,138) + rxt(k,460)*y(k,5) + rxt(k,463)*y(k,110) + rxt(k,485) &
                      *y(k,147) + rxt(k,491)*y(k,149))
         mat(k,1544) = -rxt(k,148)*y(k,132)
         mat(k,2014) = -rxt(k,149)*y(k,132)
         mat(k,1520) = -rxt(k,159)*y(k,132)
         mat(k,1328) = -rxt(k,160)*y(k,132)
         mat(k,1703) = -rxt(k,165)*y(k,132)
         mat(k,1767) = -rxt(k,176)*y(k,132)
         mat(k,1404) = -rxt(k,184)*y(k,132)
         mat(k,2048) = -rxt(k,200)*y(k,132)
         mat(k,1300) = -rxt(k,232)*y(k,132)
         mat(k,494) = -rxt(k,299)*y(k,132)
         mat(k,906) = -rxt(k,328)*y(k,132)
         mat(k,1113) = -rxt(k,358)*y(k,132)
         mat(k,1186) = -rxt(k,372)*y(k,132)
         mat(k,780) = -rxt(k,405)*y(k,132)
         mat(k,403) = -rxt(k,443)*y(k,132)
         mat(k,846) = -rxt(k,460)*y(k,132)
         mat(k,820) = -rxt(k,463)*y(k,132)
         mat(k,445) = -rxt(k,485)*y(k,132)
         mat(k,1134) = -rxt(k,491)*y(k,132)
         mat(k,1290) = .150_r8*rxt(k,313)*y(k,230)
         mat(k,1520) = mat(k,1520) + .150_r8*rxt(k,313)*y(k,224) + .150_r8*rxt(k,363) &
                      *y(k,236)
         mat(k,1260) = .150_r8*rxt(k,363)*y(k,230)
         mat(k,260) = -(rxt(k,492)*y(k,149))
         mat(k,1124) = -rxt(k,492)*y(k,133)
         mat(k,1941) = rxt(k,234)*y(k,58)
         mat(k,1967) = rxt(k,234)*y(k,18) + 2.000_r8*rxt(k,204)*y(k,58)
         mat(k,302) = -(rxt(k,482)*y(k,131) + rxt(k,483)*y(k,241))
         mat(k,1990) = -rxt(k,482)*y(k,134)
         mat(k,1602) = -rxt(k,483)*y(k,134)
         mat(k,929) = rxt(k,351)*y(k,241)
         mat(k,1332) = .100_r8*rxt(k,472)*y(k,245)
         mat(k,1585) = rxt(k,351)*y(k,92)
         mat(k,987) = .100_r8*rxt(k,472)*y(k,122)
         mat(k,404) = -(rxt(k,323)*y(k,241))
         mat(k,1617) = -rxt(k,323)*y(k,136)
         mat(k,1737) = rxt(k,325)*y(k,224)
         mat(k,1263) = rxt(k,325)*y(k,123)
         mat(k,1732) = rxt(k,445)*y(k,217)
         mat(k,452) = rxt(k,445)*y(k,123)
         mat(k,400) = -(rxt(k,442)*y(k,123) + rxt(k,443)*y(k,132))
         mat(k,1736) = -rxt(k,442)*y(k,138)
         mat(k,1887) = -rxt(k,443)*y(k,138)
         mat(k,162) = .070_r8*rxt(k,429)*y(k,241)
         mat(k,1343) = rxt(k,427)*y(k,223)
         mat(k,135) = .060_r8*rxt(k,441)*y(k,241)
         mat(k,187) = .070_r8*rxt(k,457)*y(k,241)
         mat(k,558) = rxt(k,427)*y(k,122)
         mat(k,1616) = .070_r8*rxt(k,429)*y(k,65) + .060_r8*rxt(k,441)*y(k,139) &
                      + .070_r8*rxt(k,457)*y(k,213)
         mat(k,133) = -(rxt(k,441)*y(k,241))
         mat(k,1576) = -rxt(k,441)*y(k,139)
         mat(k,125) = .530_r8*rxt(k,418)*y(k,241)
         mat(k,1576) = mat(k,1576) + .530_r8*rxt(k,418)*y(k,6)
         mat(k,265) = -(rxt(k,444)*y(k,241))
         mat(k,1595) = -rxt(k,444)*y(k,140)
         mat(k,1436) = rxt(k,439)*y(k,242)
         mat(k,390) = rxt(k,439)*y(k,230)
         mat(k,459) = -(rxt(k,340)*y(k,241))
         mat(k,1625) = -rxt(k,340)*y(k,145)
         mat(k,1459) = rxt(k,338)*y(k,243)
         mat(k,673) = rxt(k,338)*y(k,230)
         mat(k,316) = -(rxt(k,344)*y(k,241))
         mat(k,1604) = -rxt(k,344)*y(k,146)
         mat(k,1441) = .850_r8*rxt(k,342)*y(k,244)
         mat(k,1026) = .850_r8*rxt(k,342)*y(k,230)
         mat(k,441) = -(rxt(k,485)*y(k,132) + rxt(k,488)*y(k,241))
         mat(k,1888) = -rxt(k,485)*y(k,147)
         mat(k,1622) = -rxt(k,488)*y(k,147)
         mat(k,1127) = -(rxt(k,486)*y(k,18) + rxt(k,487)*y(k,58) + rxt(k,489)*y(k,123) &
                      + rxt(k,491)*y(k,132) + rxt(k,492)*y(k,133) + rxt(k,493) &
                      *y(k,241))
         mat(k,1945) = -rxt(k,486)*y(k,149)
         mat(k,1971) = -rxt(k,487)*y(k,149)
         mat(k,1752) = -rxt(k,489)*y(k,149)
         mat(k,1915) = -rxt(k,491)*y(k,149)
         mat(k,262) = -rxt(k,492)*y(k,149)
         mat(k,1684) = -rxt(k,493)*y(k,149)
         mat(k,2001) = rxt(k,482)*y(k,134)
         mat(k,1915) = mat(k,1915) + rxt(k,485)*y(k,147)
         mat(k,306) = rxt(k,482)*y(k,131)
         mat(k,442) = rxt(k,485)*y(k,132) + rxt(k,488)*y(k,241)
         mat(k,1684) = mat(k,1684) + rxt(k,488)*y(k,147)
         mat(k,754) = -(rxt(k,494)*y(k,241))
         mat(k,1657) = -rxt(k,494)*y(k,150)
         mat(k,1944) = rxt(k,486)*y(k,149)
         mat(k,1969) = rxt(k,487)*y(k,149)
         mat(k,250) = rxt(k,496)*y(k,124) + (rxt(k,497)+.500_r8*rxt(k,499))*y(k,241)
         mat(k,1745) = rxt(k,489)*y(k,149)
         mat(k,1781) = rxt(k,496)*y(k,66)
         mat(k,1894) = rxt(k,491)*y(k,149)
         mat(k,261) = rxt(k,492)*y(k,149)
         mat(k,304) = rxt(k,483)*y(k,241)
         mat(k,1126) = rxt(k,486)*y(k,18) + rxt(k,487)*y(k,58) + rxt(k,489)*y(k,123) &
                      + rxt(k,491)*y(k,132) + rxt(k,492)*y(k,133) + rxt(k,493) &
                      *y(k,241)
         mat(k,1657) = mat(k,1657) + (rxt(k,497)+.500_r8*rxt(k,499))*y(k,66) &
                      + rxt(k,483)*y(k,134) + rxt(k,493)*y(k,149)
         mat(k,212) = -(rxt(k,495)*y(k,251))
         mat(k,2056) = -rxt(k,495)*y(k,151)
         mat(k,753) = rxt(k,494)*y(k,241)
         mat(k,1589) = rxt(k,494)*y(k,150)
         mat(k,56) = .2381005_r8*rxt(k,521)*y(k,241)
         mat(k,78) = .5931005_r8*rxt(k,526)*y(k,241)
         mat(k,1550) = .2381005_r8*rxt(k,521)*y(k,103) + .5931005_r8*rxt(k,526) &
                      *y(k,200)
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
         mat(k,57) = .1308005_r8*rxt(k,521)*y(k,241)
         mat(k,79) = .1534005_r8*rxt(k,526)*y(k,241)
         mat(k,1551) = .1308005_r8*rxt(k,521)*y(k,103) + .1534005_r8*rxt(k,526) &
                      *y(k,200)
         mat(k,58) = .0348005_r8*rxt(k,521)*y(k,241)
         mat(k,80) = .0459005_r8*rxt(k,526)*y(k,241)
         mat(k,1552) = .0348005_r8*rxt(k,521)*y(k,103) + .0459005_r8*rxt(k,526) &
                      *y(k,200)
         mat(k,59) = .0076005_r8*rxt(k,521)*y(k,241)
         mat(k,81) = .0085005_r8*rxt(k,526)*y(k,241)
         mat(k,1553) = .0076005_r8*rxt(k,521)*y(k,103) + .0085005_r8*rxt(k,526) &
                      *y(k,200)
         mat(k,60) = .0113005_r8*rxt(k,521)*y(k,241)
         mat(k,82) = .0128005_r8*rxt(k,526)*y(k,241)
         mat(k,1554) = .0113005_r8*rxt(k,521)*y(k,103) + .0128005_r8*rxt(k,526) &
                      *y(k,200)
         mat(k,822) = .2202005_r8*rxt(k,515)*y(k,132) + .2202005_r8*rxt(k,516) &
                      *y(k,241)
         mat(k,760) = .0031005_r8*rxt(k,520)*y(k,241)
         mat(k,796) = .0508005_r8*rxt(k,524)*y(k,132) + .0508005_r8*rxt(k,525) &
                      *y(k,241)
         mat(k,1880) = .2202005_r8*rxt(k,515)*y(k,5) + .0508005_r8*rxt(k,524)*y(k,110)
         mat(k,1556) = .2202005_r8*rxt(k,516)*y(k,5) + .0031005_r8*rxt(k,520)*y(k,97) &
                      + .0508005_r8*rxt(k,525)*y(k,110)
         mat(k,823) = .2067005_r8*rxt(k,515)*y(k,132) + .2067005_r8*rxt(k,516) &
                      *y(k,241)
         mat(k,761) = .0035005_r8*rxt(k,520)*y(k,241)
         mat(k,797) = .1149005_r8*rxt(k,524)*y(k,132) + .1149005_r8*rxt(k,525) &
                      *y(k,241)
         mat(k,1881) = .2067005_r8*rxt(k,515)*y(k,5) + .1149005_r8*rxt(k,524)*y(k,110)
         mat(k,1557) = .2067005_r8*rxt(k,516)*y(k,5) + .0035005_r8*rxt(k,520)*y(k,97) &
                      + .1149005_r8*rxt(k,525)*y(k,110)
         mat(k,824) = .0653005_r8*rxt(k,515)*y(k,132) + .0653005_r8*rxt(k,516) &
                      *y(k,241)
         mat(k,762) = .0003005_r8*rxt(k,520)*y(k,241)
         mat(k,798) = .0348005_r8*rxt(k,524)*y(k,132) + .0348005_r8*rxt(k,525) &
                      *y(k,241)
         mat(k,1882) = .0653005_r8*rxt(k,515)*y(k,5) + .0348005_r8*rxt(k,524)*y(k,110)
         mat(k,1558) = .0653005_r8*rxt(k,516)*y(k,5) + .0003005_r8*rxt(k,520)*y(k,97) &
                      + .0348005_r8*rxt(k,525)*y(k,110)
         mat(k,825) = .1749305_r8*rxt(k,514)*y(k,124) + .1284005_r8*rxt(k,515) &
                      *y(k,132) + .1284005_r8*rxt(k,516)*y(k,241)
         mat(k,763) = .0590245_r8*rxt(k,518)*y(k,124) + .0033005_r8*rxt(k,519) &
                      *y(k,132) + .0271005_r8*rxt(k,520)*y(k,241)
         mat(k,799) = .1749305_r8*rxt(k,523)*y(k,124) + .0554005_r8*rxt(k,524) &
                      *y(k,132) + .0554005_r8*rxt(k,525)*y(k,241)
         mat(k,1773) = .1749305_r8*rxt(k,514)*y(k,5) + .0590245_r8*rxt(k,518)*y(k,97) &
                      + .1749305_r8*rxt(k,523)*y(k,110)
         mat(k,1883) = .1284005_r8*rxt(k,515)*y(k,5) + .0033005_r8*rxt(k,519)*y(k,97) &
                      + .0554005_r8*rxt(k,524)*y(k,110)
         mat(k,1559) = .1284005_r8*rxt(k,516)*y(k,5) + .0271005_r8*rxt(k,520)*y(k,97) &
                      + .0554005_r8*rxt(k,525)*y(k,110)
         mat(k,826) = .5901905_r8*rxt(k,514)*y(k,124) + .114_r8*rxt(k,515)*y(k,132) &
                      + .114_r8*rxt(k,516)*y(k,241)
         mat(k,764) = .0250245_r8*rxt(k,518)*y(k,124) + .0474005_r8*rxt(k,520) &
                      *y(k,241)
         mat(k,800) = .5901905_r8*rxt(k,523)*y(k,124) + .1278005_r8*rxt(k,524) &
                      *y(k,132) + .1278005_r8*rxt(k,525)*y(k,241)
         mat(k,1774) = .5901905_r8*rxt(k,514)*y(k,5) + .0250245_r8*rxt(k,518)*y(k,97) &
                      + .5901905_r8*rxt(k,523)*y(k,110)
         mat(k,1884) = .114_r8*rxt(k,515)*y(k,5) + .1278005_r8*rxt(k,524)*y(k,110)
         mat(k,1560) = .114_r8*rxt(k,516)*y(k,5) + .0474005_r8*rxt(k,520)*y(k,97) &
                      + .1278005_r8*rxt(k,525)*y(k,110)
         mat(k,118) = .0023005_r8*rxt(k,517)*y(k,241)
         mat(k,72) = .2381005_r8*rxt(k,522)*y(k,241)
         mat(k,84) = .5931005_r8*rxt(k,527)*y(k,241)
         mat(k,148) = .1364005_r8*rxt(k,528)*y(k,241)
         mat(k,172) = .1677005_r8*rxt(k,529)*y(k,241)
         mat(k,1561) = .0023005_r8*rxt(k,517)*y(k,6) + .2381005_r8*rxt(k,522)*y(k,104) &
                      + .5931005_r8*rxt(k,527)*y(k,201) + .1364005_r8*rxt(k,528) &
                      *y(k,209) + .1677005_r8*rxt(k,529)*y(k,211)
         mat(k,119) = .0008005_r8*rxt(k,517)*y(k,241)
         mat(k,73) = .1308005_r8*rxt(k,522)*y(k,241)
         mat(k,85) = .1534005_r8*rxt(k,527)*y(k,241)
         mat(k,149) = .0101005_r8*rxt(k,528)*y(k,241)
         mat(k,173) = .0174005_r8*rxt(k,529)*y(k,241)
         mat(k,1562) = .0008005_r8*rxt(k,517)*y(k,6) + .1308005_r8*rxt(k,522)*y(k,104) &
                      + .1534005_r8*rxt(k,527)*y(k,201) + .0101005_r8*rxt(k,528) &
                      *y(k,209) + .0174005_r8*rxt(k,529)*y(k,211)
         mat(k,120) = .0843005_r8*rxt(k,517)*y(k,241)
         mat(k,74) = .0348005_r8*rxt(k,522)*y(k,241)
         mat(k,86) = .0459005_r8*rxt(k,527)*y(k,241)
         mat(k,150) = .0763005_r8*rxt(k,528)*y(k,241)
         mat(k,174) = .086_r8*rxt(k,529)*y(k,241)
         mat(k,1563) = .0843005_r8*rxt(k,517)*y(k,6) + .0348005_r8*rxt(k,522)*y(k,104) &
                      + .0459005_r8*rxt(k,527)*y(k,201) + .0763005_r8*rxt(k,528) &
                      *y(k,209) + .086_r8*rxt(k,529)*y(k,211)
         mat(k,121) = .0443005_r8*rxt(k,517)*y(k,241)
         mat(k,75) = .0076005_r8*rxt(k,522)*y(k,241)
         mat(k,87) = .0085005_r8*rxt(k,527)*y(k,241)
         mat(k,151) = .2157005_r8*rxt(k,528)*y(k,241)
         mat(k,175) = .0512005_r8*rxt(k,529)*y(k,241)
         mat(k,1564) = .0443005_r8*rxt(k,517)*y(k,6) + .0076005_r8*rxt(k,522)*y(k,104) &
                      + .0085005_r8*rxt(k,527)*y(k,201) + .2157005_r8*rxt(k,528) &
                      *y(k,209) + .0512005_r8*rxt(k,529)*y(k,211)
         mat(k,122) = .1621005_r8*rxt(k,517)*y(k,241)
         mat(k,76) = .0113005_r8*rxt(k,522)*y(k,241)
         mat(k,88) = .0128005_r8*rxt(k,527)*y(k,241)
         mat(k,152) = .0232005_r8*rxt(k,528)*y(k,241)
         mat(k,176) = .1598005_r8*rxt(k,529)*y(k,241)
         mat(k,1565) = .1621005_r8*rxt(k,517)*y(k,6) + .0113005_r8*rxt(k,522)*y(k,104) &
                      + .0128005_r8*rxt(k,527)*y(k,201) + .0232005_r8*rxt(k,528) &
                      *y(k,209) + .1598005_r8*rxt(k,529)*y(k,211)
         mat(k,83) = -(rxt(k,526)*y(k,241))
         mat(k,1567) = -rxt(k,526)*y(k,200)
         mat(k,89) = -(rxt(k,527)*y(k,241))
         mat(k,1568) = -rxt(k,527)*y(k,201)
         mat(k,155) = .100_r8*rxt(k,449)*y(k,241)
         mat(k,177) = .230_r8*rxt(k,451)*y(k,241)
         mat(k,1581) = .100_r8*rxt(k,449)*y(k,209) + .230_r8*rxt(k,451)*y(k,211)
         mat(k,527) = -(rxt(k,473)*y(k,241))
         mat(k,1633) = -rxt(k,473)*y(k,203)
         mat(k,1462) = rxt(k,471)*y(k,245)
         mat(k,988) = rxt(k,471)*y(k,230)
         mat(k,551) = -(rxt(k,474)*y(k,241))
         mat(k,1636) = -rxt(k,474)*y(k,204)
         mat(k,1352) = .200_r8*rxt(k,467)*y(k,239) + .200_r8*rxt(k,477)*y(k,246)
         mat(k,1835) = .500_r8*rxt(k,465)*y(k,239)
         mat(k,1007) = .200_r8*rxt(k,467)*y(k,122) + .500_r8*rxt(k,465)*y(k,225)
         mat(k,966) = .200_r8*rxt(k,477)*y(k,122)
         mat(k,411) = -(rxt(k,478)*y(k,241))
         mat(k,1618) = -rxt(k,478)*y(k,205)
         mat(k,1454) = rxt(k,476)*y(k,246)
         mat(k,965) = rxt(k,476)*y(k,230)
         mat(k,881) = -(rxt(k,479)*y(k,124) + rxt(k,480)*y(k,241))
         mat(k,1788) = -rxt(k,479)*y(k,206)
         mat(k,1666) = -rxt(k,480)*y(k,206)
         mat(k,834) = .330_r8*rxt(k,460)*y(k,132)
         mat(k,808) = .330_r8*rxt(k,463)*y(k,132)
         mat(k,1370) = .800_r8*rxt(k,467)*y(k,239) + .800_r8*rxt(k,477)*y(k,246)
         mat(k,1788) = mat(k,1788) + rxt(k,468)*y(k,239)
         mat(k,1901) = .330_r8*rxt(k,460)*y(k,5) + .330_r8*rxt(k,463)*y(k,110)
         mat(k,552) = rxt(k,474)*y(k,241)
         mat(k,1842) = .500_r8*rxt(k,465)*y(k,239) + rxt(k,475)*y(k,246)
         mat(k,1009) = .800_r8*rxt(k,467)*y(k,122) + rxt(k,468)*y(k,124) &
                      + .500_r8*rxt(k,465)*y(k,225)
         mat(k,1666) = mat(k,1666) + rxt(k,474)*y(k,204)
         mat(k,969) = .800_r8*rxt(k,477)*y(k,122) + rxt(k,475)*y(k,225)
         mat(k,946) = -(rxt(k,481)*y(k,241))
         mat(k,1671) = -rxt(k,481)*y(k,207)
         mat(k,835) = .300_r8*rxt(k,460)*y(k,132)
         mat(k,809) = .300_r8*rxt(k,463)*y(k,132)
         mat(k,1374) = .900_r8*rxt(k,472)*y(k,245)
         mat(k,1904) = .300_r8*rxt(k,460)*y(k,5) + .300_r8*rxt(k,463)*y(k,110)
         mat(k,1845) = rxt(k,470)*y(k,245)
         mat(k,992) = .900_r8*rxt(k,472)*y(k,122) + rxt(k,470)*y(k,225)
         mat(k,538) = -(rxt(k,448)*y(k,241))
         mat(k,1634) = -rxt(k,448)*y(k,208)
         mat(k,1463) = rxt(k,446)*y(k,247)
         mat(k,635) = rxt(k,446)*y(k,230)
         mat(k,153) = -(rxt(k,449)*y(k,241))
         mat(k,1579) = -rxt(k,449)*y(k,209)
         mat(k,169) = -(rxt(k,415)*y(k,241))
         mat(k,1582) = -rxt(k,415)*y(k,210)
         mat(k,1433) = rxt(k,412)*y(k,248)
         mat(k,1062) = rxt(k,412)*y(k,230)
         mat(k,178) = -(rxt(k,451)*y(k,241))
         mat(k,1583) = -rxt(k,451)*y(k,211)
         mat(k,607) = -(rxt(k,454)*y(k,241))
         mat(k,1642) = -rxt(k,454)*y(k,212)
         mat(k,1469) = rxt(k,452)*y(k,249)
         mat(k,652) = rxt(k,452)*y(k,230)
         mat(k,186) = -(rxt(k,457)*y(k,241))
         mat(k,1584) = -rxt(k,457)*y(k,213)
         mat(k,179) = .150_r8*rxt(k,451)*y(k,241)
         mat(k,1584) = mat(k,1584) + .150_r8*rxt(k,451)*y(k,211)
         mat(k,370) = -(rxt(k,458)*y(k,241))
         mat(k,1612) = -rxt(k,458)*y(k,214)
         mat(k,1448) = rxt(k,455)*y(k,250)
         mat(k,427) = rxt(k,455)*y(k,230)
         mat(k,453) = -(rxt(k,416)*y(k,230) + rxt(k,417)*y(k,122) + rxt(k,445) &
                      *y(k,123))
         mat(k,1458) = -rxt(k,416)*y(k,217)
         mat(k,1347) = -rxt(k,417)*y(k,217)
         mat(k,1738) = -rxt(k,445)*y(k,217)
         mat(k,201) = rxt(k,422)*y(k,241)
         mat(k,1624) = rxt(k,422)*y(k,21)
         mat(k,853) = -(rxt(k,377)*y(k,230) + (rxt(k,378) + rxt(k,379)) * y(k,122))
         mat(k,1485) = -rxt(k,377)*y(k,218)
         mat(k,1368) = -(rxt(k,378) + rxt(k,379)) * y(k,218)
         mat(k,580) = rxt(k,380)*y(k,241)
         mat(k,195) = rxt(k,381)*y(k,241)
         mat(k,1663) = rxt(k,380)*y(k,2) + rxt(k,381)*y(k,14)
         mat(k,420) = -(rxt(k,419)*y(k,230) + rxt(k,420)*y(k,122))
         mat(k,1455) = -rxt(k,419)*y(k,219)
         mat(k,1344) = -rxt(k,420)*y(k,219)
         mat(k,126) = .350_r8*rxt(k,418)*y(k,241)
         mat(k,348) = rxt(k,421)*y(k,241)
         mat(k,1619) = .350_r8*rxt(k,418)*y(k,6) + rxt(k,421)*y(k,7)
         mat(k,378) = -(rxt(k,423)*y(k,230) + rxt(k,425)*y(k,122))
         mat(k,1449) = -rxt(k,423)*y(k,220)
         mat(k,1339) = -rxt(k,425)*y(k,220)
         mat(k,272) = rxt(k,424)*y(k,241)
         mat(k,156) = .070_r8*rxt(k,449)*y(k,241)
         mat(k,180) = .060_r8*rxt(k,451)*y(k,241)
         mat(k,1613) = rxt(k,424)*y(k,22) + .070_r8*rxt(k,449)*y(k,209) &
                      + .060_r8*rxt(k,451)*y(k,211)
         mat(k,739) = -(4._r8*rxt(k,300)*y(k,221) + rxt(k,301)*y(k,225) + rxt(k,302) &
                      *y(k,230) + rxt(k,303)*y(k,122))
         mat(k,1838) = -rxt(k,301)*y(k,221)
         mat(k,1481) = -rxt(k,302)*y(k,221)
         mat(k,1364) = -rxt(k,303)*y(k,221)
         mat(k,277) = .500_r8*rxt(k,305)*y(k,241)
         mat(k,244) = rxt(k,306)*y(k,55) + rxt(k,307)*y(k,241)
         mat(k,2028) = rxt(k,306)*y(k,27)
         mat(k,1655) = .500_r8*rxt(k,305)*y(k,26) + rxt(k,307)*y(k,27)
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
         mat(k,710) = -(rxt(k,329)*y(k,225) + rxt(k,330)*y(k,230) + rxt(k,331) &
                      *y(k,122))
         mat(k,1837) = -rxt(k,329)*y(k,222)
         mat(k,1478) = -rxt(k,330)*y(k,222)
         mat(k,1363) = -rxt(k,331)*y(k,222)
         mat(k,341) = rxt(k,332)*y(k,241)
         mat(k,100) = rxt(k,333)*y(k,241)
         mat(k,1652) = rxt(k,332)*y(k,29) + rxt(k,333)*y(k,30)
         mat(k,559) = -(rxt(k,426)*y(k,230) + rxt(k,427)*y(k,122))
         mat(k,1465) = -rxt(k,426)*y(k,223)
         mat(k,1353) = -rxt(k,427)*y(k,223)
         mat(k,222) = rxt(k,428)*y(k,241)
         mat(k,1353) = mat(k,1353) + rxt(k,417)*y(k,217)
         mat(k,1891) = rxt(k,443)*y(k,138)
         mat(k,401) = rxt(k,443)*y(k,132)
         mat(k,454) = rxt(k,417)*y(k,122) + .400_r8*rxt(k,416)*y(k,230)
         mat(k,1465) = mat(k,1465) + .400_r8*rxt(k,416)*y(k,217)
         mat(k,1637) = rxt(k,428)*y(k,31)
         mat(k,1280) = -(4._r8*rxt(k,311)*y(k,224) + rxt(k,312)*y(k,225) + rxt(k,313) &
                      *y(k,230) + rxt(k,314)*y(k,122) + rxt(k,325)*y(k,123) + rxt(k,352) &
                      *y(k,234) + rxt(k,385)*y(k,232) + rxt(k,390)*y(k,233) + rxt(k,399) &
                      *y(k,100) + rxt(k,410)*y(k,248))
         mat(k,1862) = -rxt(k,312)*y(k,224)
         mat(k,1507) = -rxt(k,313)*y(k,224)
         mat(k,1391) = -rxt(k,314)*y(k,224)
         mat(k,1754) = -rxt(k,325)*y(k,224)
         mat(k,1199) = -rxt(k,352)*y(k,224)
         mat(k,1156) = -rxt(k,385)*y(k,224)
         mat(k,1232) = -rxt(k,390)*y(k,224)
         mat(k,1092) = -rxt(k,399)*y(k,224)
         mat(k,1070) = -rxt(k,410)*y(k,224)
         mat(k,841) = .060_r8*rxt(k,460)*y(k,132)
         mat(k,1041) = rxt(k,308)*y(k,124) + rxt(k,309)*y(k,241)
         mat(k,1117) = rxt(k,334)*y(k,124) + rxt(k,335)*y(k,241)
         mat(k,436) = .500_r8*rxt(k,316)*y(k,241)
         mat(k,774) = .080_r8*rxt(k,405)*y(k,132)
         mat(k,1108) = .100_r8*rxt(k,358)*y(k,132)
         mat(k,815) = .060_r8*rxt(k,463)*y(k,132)
         mat(k,1176) = .280_r8*rxt(k,372)*y(k,132)
         mat(k,1391) = mat(k,1391) + .530_r8*rxt(k,356)*y(k,234) + rxt(k,365)*y(k,236) &
                      + rxt(k,368)*y(k,238) + rxt(k,343)*y(k,244)
         mat(k,1811) = rxt(k,308)*y(k,44) + rxt(k,334)*y(k,48) + .530_r8*rxt(k,355) &
                      *y(k,234) + rxt(k,366)*y(k,236)
         mat(k,1921) = .060_r8*rxt(k,460)*y(k,5) + .080_r8*rxt(k,405)*y(k,97) &
                      + .100_r8*rxt(k,358)*y(k,105) + .060_r8*rxt(k,463)*y(k,110) &
                      + .280_r8*rxt(k,372)*y(k,111)
         mat(k,949) = .650_r8*rxt(k,481)*y(k,241)
         mat(k,1280) = mat(k,1280) + .530_r8*rxt(k,352)*y(k,234)
         mat(k,1862) = mat(k,1862) + .260_r8*rxt(k,353)*y(k,234) + rxt(k,362)*y(k,236) &
                      + .300_r8*rxt(k,341)*y(k,244)
         mat(k,1507) = mat(k,1507) + .450_r8*rxt(k,363)*y(k,236) + .200_r8*rxt(k,367) &
                      *y(k,238) + .150_r8*rxt(k,342)*y(k,244)
         mat(k,1199) = mat(k,1199) + .530_r8*rxt(k,356)*y(k,122) + .530_r8*rxt(k,355) &
                      *y(k,124) + .530_r8*rxt(k,352)*y(k,224) + .260_r8*rxt(k,353) &
                      *y(k,225)
         mat(k,1250) = rxt(k,365)*y(k,122) + rxt(k,366)*y(k,124) + rxt(k,362)*y(k,225) &
                      + .450_r8*rxt(k,363)*y(k,230) + 4.000_r8*rxt(k,364)*y(k,236)
         mat(k,590) = rxt(k,368)*y(k,122) + .200_r8*rxt(k,367)*y(k,230)
         mat(k,1690) = rxt(k,309)*y(k,44) + rxt(k,335)*y(k,48) + .500_r8*rxt(k,316) &
                      *y(k,50) + .650_r8*rxt(k,481)*y(k,207)
         mat(k,1031) = rxt(k,343)*y(k,122) + .300_r8*rxt(k,341)*y(k,225) &
                      + .150_r8*rxt(k,342)*y(k,230)
         mat(k,1873) = -(rxt(k,201)*y(k,58) + (4._r8*rxt(k,278) + 4._r8*rxt(k,279) &
                      ) * y(k,225) + rxt(k,280)*y(k,230) + rxt(k,281)*y(k,122) &
                      + rxt(k,301)*y(k,221) + rxt(k,312)*y(k,224) + rxt(k,329) &
                      *y(k,222) + rxt(k,341)*y(k,244) + rxt(k,353)*y(k,234) + rxt(k,362) &
                      *y(k,236) + rxt(k,386)*y(k,232) + rxt(k,391)*y(k,233) + rxt(k,400) &
                      *y(k,100) + rxt(k,411)*y(k,248) + rxt(k,465)*y(k,239) + rxt(k,470) &
                      *y(k,245) + rxt(k,475)*y(k,246))
         mat(k,1983) = -rxt(k,201)*y(k,225)
         mat(k,1519) = -rxt(k,280)*y(k,225)
         mat(k,1403) = -rxt(k,281)*y(k,225)
         mat(k,746) = -rxt(k,301)*y(k,225)
         mat(k,1289) = -rxt(k,312)*y(k,225)
         mat(k,718) = -rxt(k,329)*y(k,225)
         mat(k,1037) = -rxt(k,341)*y(k,225)
         mat(k,1208) = -rxt(k,353)*y(k,225)
         mat(k,1259) = -rxt(k,362)*y(k,225)
         mat(k,1165) = -rxt(k,386)*y(k,225)
         mat(k,1241) = -rxt(k,391)*y(k,225)
         mat(k,1101) = -rxt(k,400)*y(k,225)
         mat(k,1078) = -rxt(k,411)*y(k,225)
         mat(k,1024) = -rxt(k,465)*y(k,225)
         mat(k,1004) = -rxt(k,470)*y(k,225)
         mat(k,985) = -rxt(k,475)*y(k,225)
         mat(k,905) = .280_r8*rxt(k,328)*y(k,132)
         mat(k,469) = rxt(k,315)*y(k,241)
         mat(k,356) = .700_r8*rxt(k,283)*y(k,241)
         mat(k,779) = .050_r8*rxt(k,405)*y(k,132)
         mat(k,1101) = mat(k,1101) + rxt(k,399)*y(k,224)
         mat(k,1403) = mat(k,1403) + rxt(k,314)*y(k,224) + .830_r8*rxt(k,431)*y(k,226) &
                      + .170_r8*rxt(k,437)*y(k,237)
         mat(k,1933) = .280_r8*rxt(k,328)*y(k,28) + .050_r8*rxt(k,405)*y(k,97)
         mat(k,1289) = mat(k,1289) + rxt(k,399)*y(k,100) + rxt(k,314)*y(k,122) &
                      + 4.000_r8*rxt(k,311)*y(k,224) + .900_r8*rxt(k,312)*y(k,225) &
                      + .450_r8*rxt(k,313)*y(k,230) + rxt(k,385)*y(k,232) + rxt(k,390) &
                      *y(k,233) + rxt(k,352)*y(k,234) + rxt(k,361)*y(k,236) &
                      + rxt(k,410)*y(k,248)
         mat(k,1873) = mat(k,1873) + .900_r8*rxt(k,312)*y(k,224)
         mat(k,672) = .830_r8*rxt(k,431)*y(k,122) + .330_r8*rxt(k,430)*y(k,230)
         mat(k,1519) = mat(k,1519) + .450_r8*rxt(k,313)*y(k,224) + .330_r8*rxt(k,430) &
                      *y(k,226) + .070_r8*rxt(k,436)*y(k,237)
         mat(k,1165) = mat(k,1165) + rxt(k,385)*y(k,224)
         mat(k,1241) = mat(k,1241) + rxt(k,390)*y(k,224)
         mat(k,1208) = mat(k,1208) + rxt(k,352)*y(k,224)
         mat(k,1259) = mat(k,1259) + rxt(k,361)*y(k,224)
         mat(k,790) = .170_r8*rxt(k,437)*y(k,122) + .070_r8*rxt(k,436)*y(k,230)
         mat(k,1702) = rxt(k,315)*y(k,49) + .700_r8*rxt(k,283)*y(k,52)
         mat(k,1078) = mat(k,1078) + rxt(k,410)*y(k,224)
         mat(k,665) = -(rxt(k,430)*y(k,230) + rxt(k,431)*y(k,122) + rxt(k,432) &
                      *y(k,123))
         mat(k,1474) = -rxt(k,430)*y(k,226)
         mat(k,1360) = -rxt(k,431)*y(k,226)
         mat(k,1743) = -rxt(k,432)*y(k,226)
         mat(k,495) = -((rxt(k,349) + rxt(k,350)) * y(k,122))
         mat(k,1349) = -(rxt(k,349) + rxt(k,350)) * y(k,227)
         mat(k,295) = rxt(k,348)*y(k,241)
         mat(k,1630) = rxt(k,348)*y(k,15)
         mat(k,1334) = .750_r8*rxt(k,318)*y(k,229)
         mat(k,619) = .750_r8*rxt(k,318)*y(k,122)
         mat(k,620) = -(rxt(k,317)*y(k,230) + rxt(k,318)*y(k,122))
         mat(k,1470) = -rxt(k,317)*y(k,229)
         mat(k,1356) = -rxt(k,318)*y(k,229)
         mat(k,488) = rxt(k,324)*y(k,241)
         mat(k,1643) = rxt(k,324)*y(k,24)
         mat(k,1513) = -((rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,75) + rxt(k,158) &
                      *y(k,131) + rxt(k,159)*y(k,132) + rxt(k,163)*y(k,241) &
                      + 4._r8*rxt(k,168)*y(k,230) + rxt(k,178)*y(k,124) + rxt(k,183) &
                      *y(k,122) + rxt(k,188)*y(k,123) + (rxt(k,198) + rxt(k,199) &
                      ) * y(k,55) + rxt(k,205)*y(k,58) + rxt(k,231)*y(k,16) + rxt(k,237) &
                      *y(k,18) + rxt(k,274)*y(k,41) + rxt(k,280)*y(k,225) + rxt(k,288) &
                      *y(k,231) + rxt(k,302)*y(k,221) + rxt(k,313)*y(k,224) + rxt(k,317) &
                      *y(k,229) + rxt(k,330)*y(k,222) + rxt(k,338)*y(k,243) + rxt(k,342) &
                      *y(k,244) + rxt(k,354)*y(k,234) + rxt(k,363)*y(k,236) + rxt(k,367) &
                      *y(k,238) + rxt(k,377)*y(k,218) + rxt(k,387)*y(k,232) + rxt(k,392) &
                      *y(k,233) + rxt(k,401)*y(k,100) + rxt(k,412)*y(k,248) + rxt(k,416) &
                      *y(k,217) + rxt(k,419)*y(k,219) + rxt(k,423)*y(k,220) + rxt(k,426) &
                      *y(k,223) + rxt(k,430)*y(k,226) + rxt(k,433)*y(k,235) + rxt(k,436) &
                      *y(k,237) + rxt(k,439)*y(k,242) + rxt(k,446)*y(k,247) + rxt(k,452) &
                      *y(k,249) + rxt(k,455)*y(k,250) + rxt(k,466)*y(k,239) + rxt(k,471) &
                      *y(k,245) + rxt(k,476)*y(k,246))
         mat(k,1322) = -(rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,230)
         mat(k,2007) = -rxt(k,158)*y(k,230)
         mat(k,1927) = -rxt(k,159)*y(k,230)
         mat(k,1696) = -rxt(k,163)*y(k,230)
         mat(k,1817) = -rxt(k,178)*y(k,230)
         mat(k,1397) = -rxt(k,183)*y(k,230)
         mat(k,1760) = -rxt(k,188)*y(k,230)
         mat(k,2041) = -(rxt(k,198) + rxt(k,199)) * y(k,230)
         mat(k,1977) = -rxt(k,205)*y(k,230)
         mat(k,1296) = -rxt(k,231)*y(k,230)
         mat(k,1951) = -rxt(k,237)*y(k,230)
         mat(k,1719) = -rxt(k,274)*y(k,230)
         mat(k,1867) = -rxt(k,280)*y(k,230)
         mat(k,387) = -rxt(k,288)*y(k,230)
         mat(k,742) = -rxt(k,302)*y(k,230)
         mat(k,1284) = -rxt(k,313)*y(k,230)
         mat(k,623) = -rxt(k,317)*y(k,230)
         mat(k,714) = -rxt(k,330)*y(k,230)
         mat(k,678) = -rxt(k,338)*y(k,230)
         mat(k,1033) = -rxt(k,342)*y(k,230)
         mat(k,1203) = -rxt(k,354)*y(k,230)
         mat(k,1254) = -rxt(k,363)*y(k,230)
         mat(k,592) = -rxt(k,367)*y(k,230)
         mat(k,858) = -rxt(k,377)*y(k,230)
         mat(k,1160) = -rxt(k,387)*y(k,230)
         mat(k,1236) = -rxt(k,392)*y(k,230)
         mat(k,1096) = -rxt(k,401)*y(k,230)
         mat(k,1073) = -rxt(k,412)*y(k,230)
         mat(k,456) = -rxt(k,416)*y(k,230)
         mat(k,424) = -rxt(k,419)*y(k,230)
         mat(k,381) = -rxt(k,423)*y(k,230)
         mat(k,561) = -rxt(k,426)*y(k,230)
         mat(k,669) = -rxt(k,430)*y(k,230)
         mat(k,631) = -rxt(k,433)*y(k,230)
         mat(k,787) = -rxt(k,436)*y(k,230)
         mat(k,394) = -rxt(k,439)*y(k,230)
         mat(k,645) = -rxt(k,446)*y(k,230)
         mat(k,662) = -rxt(k,452)*y(k,230)
         mat(k,432) = -rxt(k,455)*y(k,230)
         mat(k,1019) = -rxt(k,466)*y(k,230)
         mat(k,1000) = -rxt(k,471)*y(k,230)
         mat(k,980) = -rxt(k,476)*y(k,230)
         mat(k,842) = .570_r8*rxt(k,460)*y(k,132)
         mat(k,127) = .650_r8*rxt(k,418)*y(k,241)
         mat(k,1296) = mat(k,1296) + rxt(k,230)*y(k,41)
         mat(k,1951) = mat(k,1951) + rxt(k,242)*y(k,241)
         mat(k,241) = .350_r8*rxt(k,297)*y(k,241)
         mat(k,491) = .130_r8*rxt(k,299)*y(k,132)
         mat(k,218) = rxt(k,304)*y(k,241)
         mat(k,900) = .280_r8*rxt(k,328)*y(k,132)
         mat(k,1719) = mat(k,1719) + rxt(k,230)*y(k,16) + rxt(k,194)*y(k,55) &
                      + rxt(k,275)*y(k,124) + rxt(k,276)*y(k,131)
         mat(k,94) = rxt(k,310)*y(k,241)
         mat(k,702) = rxt(k,282)*y(k,241)
         mat(k,2041) = mat(k,2041) + rxt(k,194)*y(k,41) + rxt(k,197)*y(k,78)
         mat(k,1977) = mat(k,1977) + rxt(k,201)*y(k,225) + rxt(k,212)*y(k,241)
         mat(k,957) = rxt(k,285)*y(k,241)
         mat(k,163) = .730_r8*rxt(k,429)*y(k,241)
         mat(k,252) = .500_r8*rxt(k,499)*y(k,241)
         mat(k,926) = rxt(k,321)*y(k,241)
         mat(k,794) = rxt(k,322)*y(k,241)
         mat(k,522) = rxt(k,197)*y(k,55) + rxt(k,153)*y(k,131) + rxt(k,162)*y(k,241)
         mat(k,145) = rxt(k,286)*y(k,241)
         mat(k,706) = rxt(k,287)*y(k,241)
         mat(k,939) = rxt(k,351)*y(k,241)
         mat(k,962) = rxt(k,336)*y(k,241)
         mat(k,775) = .370_r8*rxt(k,405)*y(k,132)
         mat(k,516) = .300_r8*rxt(k,396)*y(k,241)
         mat(k,477) = rxt(k,397)*y(k,241)
         mat(k,1096) = mat(k,1096) + rxt(k,402)*y(k,122) + rxt(k,403)*y(k,124) &
                      + rxt(k,399)*y(k,224) + 1.200_r8*rxt(k,400)*y(k,225)
         mat(k,325) = rxt(k,404)*y(k,241)
         mat(k,1110) = .140_r8*rxt(k,358)*y(k,132)
         mat(k,292) = .200_r8*rxt(k,360)*y(k,241)
         mat(k,482) = .500_r8*rxt(k,371)*y(k,241)
         mat(k,816) = .570_r8*rxt(k,463)*y(k,132)
         mat(k,1180) = .280_r8*rxt(k,372)*y(k,132)
         mat(k,314) = rxt(k,408)*y(k,241)
         mat(k,918) = rxt(k,409)*y(k,241)
         mat(k,1397) = mat(k,1397) + rxt(k,402)*y(k,100) + rxt(k,378)*y(k,218) &
                      + rxt(k,420)*y(k,219) + rxt(k,425)*y(k,220) + rxt(k,303) &
                      *y(k,221) + rxt(k,331)*y(k,222) + rxt(k,281)*y(k,225) &
                      + .170_r8*rxt(k,431)*y(k,226) + rxt(k,349)*y(k,227) &
                      + .250_r8*rxt(k,318)*y(k,229) + rxt(k,290)*y(k,231) &
                      + .920_r8*rxt(k,388)*y(k,232) + .920_r8*rxt(k,394)*y(k,233) &
                      + .470_r8*rxt(k,356)*y(k,234) + .400_r8*rxt(k,434)*y(k,235) &
                      + .830_r8*rxt(k,437)*y(k,237) + rxt(k,440)*y(k,242) + rxt(k,339) &
                      *y(k,243) + .900_r8*rxt(k,472)*y(k,245) + .800_r8*rxt(k,477) &
                      *y(k,246) + rxt(k,447)*y(k,247) + rxt(k,413)*y(k,248) &
                      + rxt(k,453)*y(k,249) + rxt(k,456)*y(k,250)
         mat(k,1817) = mat(k,1817) + rxt(k,275)*y(k,41) + rxt(k,403)*y(k,100) &
                      + rxt(k,389)*y(k,232) + rxt(k,395)*y(k,233) + .470_r8*rxt(k,355) &
                      *y(k,234) + rxt(k,181)*y(k,241) + rxt(k,414)*y(k,248)
         mat(k,2007) = mat(k,2007) + rxt(k,276)*y(k,41) + rxt(k,153)*y(k,78)
         mat(k,1927) = mat(k,1927) + .570_r8*rxt(k,460)*y(k,5) + .130_r8*rxt(k,299) &
                      *y(k,24) + .280_r8*rxt(k,328)*y(k,28) + .370_r8*rxt(k,405) &
                      *y(k,97) + .140_r8*rxt(k,358)*y(k,105) + .570_r8*rxt(k,463) &
                      *y(k,110) + .280_r8*rxt(k,372)*y(k,111) + rxt(k,165)*y(k,241)
         mat(k,136) = .800_r8*rxt(k,441)*y(k,241)
         mat(k,756) = rxt(k,494)*y(k,241)
         mat(k,950) = .200_r8*rxt(k,481)*y(k,241)
         mat(k,158) = .280_r8*rxt(k,449)*y(k,241)
         mat(k,184) = .380_r8*rxt(k,451)*y(k,241)
         mat(k,189) = .630_r8*rxt(k,457)*y(k,241)
         mat(k,858) = mat(k,858) + rxt(k,378)*y(k,122)
         mat(k,424) = mat(k,424) + rxt(k,420)*y(k,122)
         mat(k,381) = mat(k,381) + rxt(k,425)*y(k,122)
         mat(k,742) = mat(k,742) + rxt(k,303)*y(k,122) + 2.400_r8*rxt(k,300)*y(k,221) &
                      + rxt(k,301)*y(k,225)
         mat(k,714) = mat(k,714) + rxt(k,331)*y(k,122) + rxt(k,329)*y(k,225)
         mat(k,1284) = mat(k,1284) + rxt(k,399)*y(k,100) + .900_r8*rxt(k,312)*y(k,225) &
                      + rxt(k,385)*y(k,232) + rxt(k,390)*y(k,233) + .470_r8*rxt(k,352) &
                      *y(k,234) + rxt(k,410)*y(k,248)
         mat(k,1867) = mat(k,1867) + rxt(k,201)*y(k,58) + 1.200_r8*rxt(k,400)*y(k,100) &
                      + rxt(k,281)*y(k,122) + rxt(k,301)*y(k,221) + rxt(k,329) &
                      *y(k,222) + .900_r8*rxt(k,312)*y(k,224) + 4.000_r8*rxt(k,278) &
                      *y(k,225) + rxt(k,386)*y(k,232) + rxt(k,391)*y(k,233) &
                      + .730_r8*rxt(k,353)*y(k,234) + rxt(k,362)*y(k,236) &
                      + .500_r8*rxt(k,465)*y(k,239) + .300_r8*rxt(k,341)*y(k,244) &
                      + rxt(k,470)*y(k,245) + rxt(k,475)*y(k,246) + .800_r8*rxt(k,411) &
                      *y(k,248)
         mat(k,669) = mat(k,669) + .170_r8*rxt(k,431)*y(k,122) + .070_r8*rxt(k,430) &
                      *y(k,230)
         mat(k,500) = rxt(k,349)*y(k,122)
         mat(k,623) = mat(k,623) + .250_r8*rxt(k,318)*y(k,122)
         mat(k,1513) = mat(k,1513) + .070_r8*rxt(k,430)*y(k,226) + .160_r8*rxt(k,433) &
                      *y(k,235) + .330_r8*rxt(k,436)*y(k,237)
         mat(k,387) = mat(k,387) + rxt(k,290)*y(k,122)
         mat(k,1160) = mat(k,1160) + .920_r8*rxt(k,388)*y(k,122) + rxt(k,389)*y(k,124) &
                      + rxt(k,385)*y(k,224) + rxt(k,386)*y(k,225)
         mat(k,1236) = mat(k,1236) + .920_r8*rxt(k,394)*y(k,122) + rxt(k,395)*y(k,124) &
                      + rxt(k,390)*y(k,224) + rxt(k,391)*y(k,225)
         mat(k,1203) = mat(k,1203) + .470_r8*rxt(k,356)*y(k,122) + .470_r8*rxt(k,355) &
                      *y(k,124) + .470_r8*rxt(k,352)*y(k,224) + .730_r8*rxt(k,353) &
                      *y(k,225)
         mat(k,631) = mat(k,631) + .400_r8*rxt(k,434)*y(k,122) + .160_r8*rxt(k,433) &
                      *y(k,230)
         mat(k,1254) = mat(k,1254) + rxt(k,362)*y(k,225)
         mat(k,787) = mat(k,787) + .830_r8*rxt(k,437)*y(k,122) + .330_r8*rxt(k,436) &
                      *y(k,230)
         mat(k,1019) = mat(k,1019) + .500_r8*rxt(k,465)*y(k,225)
         mat(k,1696) = mat(k,1696) + .650_r8*rxt(k,418)*y(k,6) + rxt(k,242)*y(k,18) &
                      + .350_r8*rxt(k,297)*y(k,23) + rxt(k,304)*y(k,25) + rxt(k,310) &
                      *y(k,46) + rxt(k,282)*y(k,51) + rxt(k,212)*y(k,58) + rxt(k,285) &
                      *y(k,61) + .730_r8*rxt(k,429)*y(k,65) + .500_r8*rxt(k,499) &
                      *y(k,66) + rxt(k,321)*y(k,73) + rxt(k,322)*y(k,74) + rxt(k,162) &
                      *y(k,78) + rxt(k,286)*y(k,85) + rxt(k,287)*y(k,86) + rxt(k,351) &
                      *y(k,92) + rxt(k,336)*y(k,94) + .300_r8*rxt(k,396)*y(k,98) &
                      + rxt(k,397)*y(k,99) + rxt(k,404)*y(k,101) + .200_r8*rxt(k,360) &
                      *y(k,106) + .500_r8*rxt(k,371)*y(k,109) + rxt(k,408)*y(k,115) &
                      + rxt(k,409)*y(k,116) + rxt(k,181)*y(k,124) + rxt(k,165) &
                      *y(k,132) + .800_r8*rxt(k,441)*y(k,139) + rxt(k,494)*y(k,150) &
                      + .200_r8*rxt(k,481)*y(k,207) + .280_r8*rxt(k,449)*y(k,209) &
                      + .380_r8*rxt(k,451)*y(k,211) + .630_r8*rxt(k,457)*y(k,213)
         mat(k,394) = mat(k,394) + rxt(k,440)*y(k,122)
         mat(k,678) = mat(k,678) + rxt(k,339)*y(k,122)
         mat(k,1033) = mat(k,1033) + .300_r8*rxt(k,341)*y(k,225)
         mat(k,1000) = mat(k,1000) + .900_r8*rxt(k,472)*y(k,122) + rxt(k,470)*y(k,225)
         mat(k,980) = mat(k,980) + .800_r8*rxt(k,477)*y(k,122) + rxt(k,475)*y(k,225)
         mat(k,645) = mat(k,645) + rxt(k,447)*y(k,122)
         mat(k,1073) = mat(k,1073) + rxt(k,413)*y(k,122) + rxt(k,414)*y(k,124) &
                      + rxt(k,410)*y(k,224) + .800_r8*rxt(k,411)*y(k,225)
         mat(k,662) = mat(k,662) + rxt(k,453)*y(k,122)
         mat(k,432) = mat(k,432) + rxt(k,456)*y(k,122)
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
         mat(k,384) = -(rxt(k,288)*y(k,230) + rxt(k,290)*y(k,122))
         mat(k,1450) = -rxt(k,288)*y(k,231)
         mat(k,1340) = -rxt(k,290)*y(k,231)
         mat(k,1709) = rxt(k,274)*y(k,230)
         mat(k,1450) = mat(k,1450) + rxt(k,274)*y(k,41)
         mat(k,1152) = -(rxt(k,385)*y(k,224) + rxt(k,386)*y(k,225) + rxt(k,387) &
                      *y(k,230) + rxt(k,388)*y(k,122) + rxt(k,389)*y(k,124))
         mat(k,1275) = -rxt(k,385)*y(k,232)
         mat(k,1857) = -rxt(k,386)*y(k,232)
         mat(k,1502) = -rxt(k,387)*y(k,232)
         mat(k,1386) = -rxt(k,388)*y(k,232)
         mat(k,1806) = -rxt(k,389)*y(k,232)
         mat(k,771) = .600_r8*rxt(k,406)*y(k,241)
         mat(k,1685) = .600_r8*rxt(k,406)*y(k,97)
         mat(k,1230) = -(rxt(k,390)*y(k,224) + rxt(k,391)*y(k,225) + rxt(k,392) &
                      *y(k,230) + rxt(k,394)*y(k,122) + rxt(k,395)*y(k,124))
         mat(k,1278) = -rxt(k,390)*y(k,233)
         mat(k,1860) = -rxt(k,391)*y(k,233)
         mat(k,1505) = -rxt(k,392)*y(k,233)
         mat(k,1389) = -rxt(k,394)*y(k,233)
         mat(k,1809) = -rxt(k,395)*y(k,233)
         mat(k,773) = .400_r8*rxt(k,406)*y(k,241)
         mat(k,1688) = .400_r8*rxt(k,406)*y(k,97)
         mat(k,1197) = -(rxt(k,352)*y(k,224) + rxt(k,353)*y(k,225) + rxt(k,354) &
                      *y(k,230) + rxt(k,355)*y(k,124) + (rxt(k,356) + rxt(k,357) &
                      ) * y(k,122))
         mat(k,1277) = -rxt(k,352)*y(k,234)
         mat(k,1859) = -rxt(k,353)*y(k,234)
         mat(k,1504) = -rxt(k,354)*y(k,234)
         mat(k,1808) = -rxt(k,355)*y(k,234)
         mat(k,1388) = -(rxt(k,356) + rxt(k,357)) * y(k,234)
         mat(k,1106) = .500_r8*rxt(k,359)*y(k,241)
         mat(k,290) = .200_r8*rxt(k,360)*y(k,241)
         mat(k,1175) = rxt(k,373)*y(k,241)
         mat(k,1687) = .500_r8*rxt(k,359)*y(k,105) + .200_r8*rxt(k,360)*y(k,106) &
                      + rxt(k,373)*y(k,111)
         mat(k,627) = -(rxt(k,433)*y(k,230) + rxt(k,434)*y(k,122) + rxt(k,435) &
                      *y(k,123))
         mat(k,1471) = -rxt(k,433)*y(k,235)
         mat(k,1357) = -rxt(k,434)*y(k,235)
         mat(k,1742) = -rxt(k,435)*y(k,235)
         mat(k,1249) = -(rxt(k,361)*y(k,224) + rxt(k,362)*y(k,225) + rxt(k,363) &
                      *y(k,230) + 4._r8*rxt(k,364)*y(k,236) + rxt(k,365)*y(k,122) &
                      + rxt(k,366)*y(k,124) + rxt(k,374)*y(k,123))
         mat(k,1279) = -rxt(k,361)*y(k,236)
         mat(k,1861) = -rxt(k,362)*y(k,236)
         mat(k,1506) = -rxt(k,363)*y(k,236)
         mat(k,1390) = -rxt(k,365)*y(k,236)
         mat(k,1810) = -rxt(k,366)*y(k,236)
         mat(k,1753) = -rxt(k,374)*y(k,236)
         mat(k,1107) = .500_r8*rxt(k,359)*y(k,241)
         mat(k,291) = .500_r8*rxt(k,360)*y(k,241)
         mat(k,1689) = .500_r8*rxt(k,359)*y(k,105) + .500_r8*rxt(k,360)*y(k,106)
         mat(k,782) = -(rxt(k,436)*y(k,230) + rxt(k,437)*y(k,122) + rxt(k,438) &
                      *y(k,123))
         mat(k,1483) = -rxt(k,436)*y(k,237)
         mat(k,1366) = -rxt(k,437)*y(k,237)
         mat(k,1746) = -rxt(k,438)*y(k,237)
         mat(k,588) = -(rxt(k,367)*y(k,230) + rxt(k,368)*y(k,122))
         mat(k,1467) = -rxt(k,367)*y(k,238)
         mat(k,1355) = -rxt(k,368)*y(k,238)
         mat(k,448) = rxt(k,369)*y(k,241)
         mat(k,285) = rxt(k,370)*y(k,241)
         mat(k,1640) = rxt(k,369)*y(k,107) + rxt(k,370)*y(k,108)
         mat(k,1013) = -(rxt(k,465)*y(k,225) + rxt(k,466)*y(k,230) + rxt(k,467) &
                      *y(k,122) + rxt(k,468)*y(k,124))
         mat(k,1850) = -rxt(k,465)*y(k,239)
         mat(k,1494) = -rxt(k,466)*y(k,239)
         mat(k,1379) = -rxt(k,467)*y(k,239)
         mat(k,1798) = -rxt(k,468)*y(k,239)
         mat(k,838) = rxt(k,459)*y(k,124)
         mat(k,812) = rxt(k,462)*y(k,124)
         mat(k,1798) = mat(k,1798) + rxt(k,459)*y(k,5) + rxt(k,462)*y(k,110) &
                      + .500_r8*rxt(k,479)*y(k,206)
         mat(k,330) = rxt(k,469)*y(k,241)
         mat(k,885) = .500_r8*rxt(k,479)*y(k,124)
         mat(k,1676) = rxt(k,469)*y(k,126)
         mat(k,1538) = -(rxt(k,144)*y(k,76) + rxt(k,145)*y(k,251) + rxt(k,148) &
                      *y(k,132) + (rxt(k,226) + rxt(k,227)) * y(k,84) + (rxt(k,249) &
                      + rxt(k,250)) * y(k,80) + rxt(k,255)*y(k,63) + rxt(k,256) &
                      *y(k,64) + rxt(k,294)*y(k,85))
         mat(k,1055) = -rxt(k,144)*y(k,240)
         mat(k,2067) = -rxt(k,145)*y(k,240)
         mat(k,1928) = -rxt(k,148)*y(k,240)
         mat(k,1420) = -(rxt(k,226) + rxt(k,227)) * y(k,240)
         mat(k,723) = -(rxt(k,249) + rxt(k,250)) * y(k,240)
         mat(k,111) = -rxt(k,255)*y(k,240)
         mat(k,142) = -rxt(k,256)*y(k,240)
         mat(k,146) = -rxt(k,294)*y(k,240)
         mat(k,1698) = -(rxt(k,161)*y(k,76) + rxt(k,162)*y(k,78) + rxt(k,163)*y(k,230) &
                      + rxt(k,164)*y(k,131) + rxt(k,165)*y(k,132) + (4._r8*rxt(k,166) &
                      + 4._r8*rxt(k,167)) * y(k,241) + rxt(k,169)*y(k,89) + rxt(k,181) &
                      *y(k,124) + rxt(k,182)*y(k,112) + rxt(k,190)*y(k,123) + rxt(k,191) &
                      *y(k,88) + rxt(k,210)*y(k,59) + (rxt(k,212) + rxt(k,213) &
                      ) * y(k,58) + rxt(k,215)*y(k,84) + rxt(k,218)*y(k,91) + rxt(k,242) &
                      *y(k,18) + rxt(k,244)*y(k,80) + rxt(k,277)*y(k,41) + rxt(k,282) &
                      *y(k,51) + rxt(k,283)*y(k,52) + (rxt(k,285) + rxt(k,295) &
                      ) * y(k,61) + rxt(k,286)*y(k,85) + rxt(k,287)*y(k,86) + rxt(k,297) &
                      *y(k,23) + rxt(k,304)*y(k,25) + rxt(k,305)*y(k,26) + rxt(k,307) &
                      *y(k,27) + rxt(k,309)*y(k,44) + rxt(k,310)*y(k,46) + rxt(k,315) &
                      *y(k,49) + rxt(k,316)*y(k,50) + rxt(k,321)*y(k,73) + rxt(k,322) &
                      *y(k,74) + rxt(k,323)*y(k,136) + rxt(k,324)*y(k,24) + rxt(k,332) &
                      *y(k,29) + rxt(k,333)*y(k,30) + rxt(k,335)*y(k,48) + rxt(k,336) &
                      *y(k,94) + rxt(k,337)*y(k,125) + rxt(k,340)*y(k,145) + rxt(k,344) &
                      *y(k,146) + rxt(k,345)*y(k,28) + rxt(k,346)*y(k,47) + rxt(k,348) &
                      *y(k,15) + rxt(k,351)*y(k,92) + rxt(k,359)*y(k,105) + rxt(k,360) &
                      *y(k,106) + rxt(k,369)*y(k,107) + rxt(k,370)*y(k,108) + rxt(k,371) &
                      *y(k,109) + rxt(k,373)*y(k,111) + rxt(k,376)*y(k,1) + rxt(k,380) &
                      *y(k,2) + rxt(k,381)*y(k,14) + rxt(k,382)*y(k,93) + rxt(k,383) &
                      *y(k,95) + rxt(k,384)*y(k,96) + rxt(k,396)*y(k,98) + rxt(k,397) &
                      *y(k,99) + rxt(k,404)*y(k,101) + rxt(k,406)*y(k,97) + rxt(k,407) &
                      *y(k,102) + rxt(k,408)*y(k,115) + rxt(k,409)*y(k,116) + rxt(k,415) &
                      *y(k,210) + rxt(k,418)*y(k,6) + rxt(k,421)*y(k,7) + rxt(k,422) &
                      *y(k,21) + rxt(k,424)*y(k,22) + rxt(k,428)*y(k,31) + rxt(k,429) &
                      *y(k,65) + rxt(k,441)*y(k,139) + rxt(k,444)*y(k,140) + rxt(k,448) &
                      *y(k,208) + rxt(k,449)*y(k,209) + rxt(k,451)*y(k,211) + rxt(k,454) &
                      *y(k,212) + rxt(k,457)*y(k,213) + rxt(k,458)*y(k,214) + rxt(k,461) &
                      *y(k,5) + rxt(k,464)*y(k,110) + rxt(k,469)*y(k,126) + rxt(k,473) &
                      *y(k,203) + rxt(k,474)*y(k,204) + rxt(k,478)*y(k,205) + rxt(k,480) &
                      *y(k,206) + rxt(k,481)*y(k,207) + rxt(k,483)*y(k,134) + rxt(k,488) &
                      *y(k,147) + rxt(k,493)*y(k,149) + rxt(k,494)*y(k,150) + (rxt(k,497) &
                      + rxt(k,499)) * y(k,66) + rxt(k,498)*y(k,120))
         mat(k,1056) = -rxt(k,161)*y(k,241)
         mat(k,523) = -rxt(k,162)*y(k,241)
         mat(k,1515) = -rxt(k,163)*y(k,241)
         mat(k,2009) = -rxt(k,164)*y(k,241)
         mat(k,1929) = -rxt(k,165)*y(k,241)
         mat(k,360) = -rxt(k,169)*y(k,241)
         mat(k,1819) = -rxt(k,181)*y(k,241)
         mat(k,367) = -rxt(k,182)*y(k,241)
         mat(k,1762) = -rxt(k,190)*y(k,241)
         mat(k,1311) = -rxt(k,191)*y(k,241)
         mat(k,874) = -rxt(k,210)*y(k,241)
         mat(k,1979) = -(rxt(k,212) + rxt(k,213)) * y(k,241)
         mat(k,1421) = -rxt(k,215)*y(k,241)
         mat(k,731) = -rxt(k,218)*y(k,241)
         mat(k,1953) = -rxt(k,242)*y(k,241)
         mat(k,724) = -rxt(k,244)*y(k,241)
         mat(k,1721) = -rxt(k,277)*y(k,241)
         mat(k,703) = -rxt(k,282)*y(k,241)
         mat(k,354) = -rxt(k,283)*y(k,241)
         mat(k,958) = -(rxt(k,285) + rxt(k,295)) * y(k,241)
         mat(k,147) = -rxt(k,286)*y(k,241)
         mat(k,707) = -rxt(k,287)*y(k,241)
         mat(k,242) = -rxt(k,297)*y(k,241)
         mat(k,219) = -rxt(k,304)*y(k,241)
         mat(k,280) = -rxt(k,305)*y(k,241)
         mat(k,246) = -rxt(k,307)*y(k,241)
         mat(k,1045) = -rxt(k,309)*y(k,241)
         mat(k,95) = -rxt(k,310)*y(k,241)
         mat(k,468) = -rxt(k,315)*y(k,241)
         mat(k,437) = -rxt(k,316)*y(k,241)
         mat(k,927) = -rxt(k,321)*y(k,241)
         mat(k,795) = -rxt(k,322)*y(k,241)
         mat(k,406) = -rxt(k,323)*y(k,241)
         mat(k,492) = -rxt(k,324)*y(k,241)
         mat(k,344) = -rxt(k,332)*y(k,241)
         mat(k,101) = -rxt(k,333)*y(k,241)
         mat(k,1121) = -rxt(k,335)*y(k,241)
         mat(k,963) = -rxt(k,336)*y(k,241)
         mat(k,750) = -rxt(k,337)*y(k,241)
         mat(k,464) = -rxt(k,340)*y(k,241)
         mat(k,319) = -rxt(k,344)*y(k,241)
         mat(k,901) = -rxt(k,345)*y(k,241)
         mat(k,866) = -rxt(k,346)*y(k,241)
         mat(k,298) = -rxt(k,348)*y(k,241)
         mat(k,940) = -rxt(k,351)*y(k,241)
         mat(k,1111) = -rxt(k,359)*y(k,241)
         mat(k,293) = -rxt(k,360)*y(k,241)
         mat(k,451) = -rxt(k,369)*y(k,241)
         mat(k,288) = -rxt(k,370)*y(k,241)
         mat(k,483) = -rxt(k,371)*y(k,241)
         mat(k,1181) = -rxt(k,373)*y(k,241)
         mat(k,573) = -rxt(k,376)*y(k,241)
         mat(k,585) = -rxt(k,380)*y(k,241)
         mat(k,196) = -rxt(k,381)*y(k,241)
         mat(k,207) = -rxt(k,382)*y(k,241)
         mat(k,283) = -rxt(k,383)*y(k,241)
         mat(k,105) = -rxt(k,384)*y(k,241)
         mat(k,517) = -rxt(k,396)*y(k,241)
         mat(k,478) = -rxt(k,397)*y(k,241)
         mat(k,326) = -rxt(k,404)*y(k,241)
         mat(k,776) = -rxt(k,406)*y(k,241)
         mat(k,601) = -rxt(k,407)*y(k,241)
         mat(k,315) = -rxt(k,408)*y(k,241)
         mat(k,919) = -rxt(k,409)*y(k,241)
         mat(k,171) = -rxt(k,415)*y(k,241)
         mat(k,128) = -rxt(k,418)*y(k,241)
         mat(k,351) = -rxt(k,421)*y(k,241)
         mat(k,202) = -rxt(k,422)*y(k,241)
         mat(k,275) = -rxt(k,424)*y(k,241)
         mat(k,223) = -rxt(k,428)*y(k,241)
         mat(k,164) = -rxt(k,429)*y(k,241)
         mat(k,137) = -rxt(k,441)*y(k,241)
         mat(k,269) = -rxt(k,444)*y(k,241)
         mat(k,546) = -rxt(k,448)*y(k,241)
         mat(k,159) = -rxt(k,449)*y(k,241)
         mat(k,185) = -rxt(k,451)*y(k,241)
         mat(k,617) = -rxt(k,454)*y(k,241)
         mat(k,190) = -rxt(k,457)*y(k,241)
         mat(k,375) = -rxt(k,458)*y(k,241)
         mat(k,843) = -rxt(k,461)*y(k,241)
         mat(k,817) = -rxt(k,464)*y(k,241)
         mat(k,332) = -rxt(k,469)*y(k,241)
         mat(k,534) = -rxt(k,473)*y(k,241)
         mat(k,555) = -rxt(k,474)*y(k,241)
         mat(k,416) = -rxt(k,478)*y(k,241)
         mat(k,887) = -rxt(k,480)*y(k,241)
         mat(k,951) = -rxt(k,481)*y(k,241)
         mat(k,308) = -rxt(k,483)*y(k,241)
         mat(k,444) = -rxt(k,488)*y(k,241)
         mat(k,1132) = -rxt(k,493)*y(k,241)
         mat(k,757) = -rxt(k,494)*y(k,241)
         mat(k,253) = -(rxt(k,497) + rxt(k,499)) * y(k,241)
         mat(k,91) = -rxt(k,498)*y(k,241)
         mat(k,843) = mat(k,843) + .630_r8*rxt(k,460)*y(k,132)
         mat(k,242) = mat(k,242) + .650_r8*rxt(k,297)*y(k,241)
         mat(k,492) = mat(k,492) + .130_r8*rxt(k,299)*y(k,132)
         mat(k,280) = mat(k,280) + .500_r8*rxt(k,305)*y(k,241)
         mat(k,901) = mat(k,901) + .360_r8*rxt(k,328)*y(k,132)
         mat(k,1721) = mat(k,1721) + rxt(k,276)*y(k,131)
         mat(k,354) = mat(k,354) + .300_r8*rxt(k,283)*y(k,241)
         mat(k,2043) = rxt(k,199)*y(k,230)
         mat(k,697) = rxt(k,253)*y(k,251)
         mat(k,1324) = rxt(k,160)*y(k,132) + 2.000_r8*rxt(k,155)*y(k,230)
         mat(k,1056) = mat(k,1056) + rxt(k,152)*y(k,131) + rxt(k,144)*y(k,240)
         mat(k,523) = mat(k,523) + rxt(k,153)*y(k,131)
         mat(k,724) = mat(k,724) + rxt(k,243)*y(k,131) + rxt(k,249)*y(k,240)
         mat(k,1421) = mat(k,1421) + rxt(k,214)*y(k,131) + rxt(k,226)*y(k,240)
         mat(k,147) = mat(k,147) + rxt(k,294)*y(k,240)
         mat(k,687) = rxt(k,245)*y(k,131)
         mat(k,731) = mat(k,731) + rxt(k,217)*y(k,131)
         mat(k,776) = mat(k,776) + .320_r8*rxt(k,405)*y(k,132)
         mat(k,601) = mat(k,601) + .600_r8*rxt(k,407)*y(k,241)
         mat(k,1111) = mat(k,1111) + .240_r8*rxt(k,358)*y(k,132)
         mat(k,293) = mat(k,293) + .100_r8*rxt(k,360)*y(k,241)
         mat(k,817) = mat(k,817) + .630_r8*rxt(k,463)*y(k,132)
         mat(k,1181) = mat(k,1181) + .360_r8*rxt(k,372)*y(k,132)
         mat(k,1399) = rxt(k,183)*y(k,230)
         mat(k,1819) = mat(k,1819) + rxt(k,178)*y(k,230)
         mat(k,2009) = mat(k,2009) + rxt(k,276)*y(k,41) + rxt(k,152)*y(k,76) &
                      + rxt(k,153)*y(k,78) + rxt(k,243)*y(k,80) + rxt(k,214)*y(k,84) &
                      + rxt(k,245)*y(k,90) + rxt(k,217)*y(k,91) + rxt(k,158)*y(k,230)
         mat(k,1929) = mat(k,1929) + .630_r8*rxt(k,460)*y(k,5) + .130_r8*rxt(k,299) &
                      *y(k,24) + .360_r8*rxt(k,328)*y(k,28) + rxt(k,160)*y(k,75) &
                      + .320_r8*rxt(k,405)*y(k,97) + .240_r8*rxt(k,358)*y(k,105) &
                      + .630_r8*rxt(k,463)*y(k,110) + .360_r8*rxt(k,372)*y(k,111) &
                      + rxt(k,159)*y(k,230)
         mat(k,464) = mat(k,464) + .500_r8*rxt(k,340)*y(k,241)
         mat(k,171) = mat(k,171) + .500_r8*rxt(k,415)*y(k,241)
         mat(k,457) = .400_r8*rxt(k,416)*y(k,230)
         mat(k,1285) = .450_r8*rxt(k,313)*y(k,230)
         mat(k,670) = .400_r8*rxt(k,430)*y(k,230)
         mat(k,1515) = mat(k,1515) + rxt(k,199)*y(k,55) + 2.000_r8*rxt(k,155)*y(k,75) &
                      + rxt(k,183)*y(k,122) + rxt(k,178)*y(k,124) + rxt(k,158) &
                      *y(k,131) + rxt(k,159)*y(k,132) + .400_r8*rxt(k,416)*y(k,217) &
                      + .450_r8*rxt(k,313)*y(k,224) + .400_r8*rxt(k,430)*y(k,226) &
                      + .450_r8*rxt(k,363)*y(k,236) + .400_r8*rxt(k,436)*y(k,237) &
                      + .200_r8*rxt(k,367)*y(k,238) + .150_r8*rxt(k,342)*y(k,244)
         mat(k,1255) = .450_r8*rxt(k,363)*y(k,230)
         mat(k,788) = .400_r8*rxt(k,436)*y(k,230)
         mat(k,593) = .200_r8*rxt(k,367)*y(k,230)
         mat(k,1539) = rxt(k,144)*y(k,76) + rxt(k,249)*y(k,80) + rxt(k,226)*y(k,84) &
                      + rxt(k,294)*y(k,85) + 2.000_r8*rxt(k,145)*y(k,251)
         mat(k,1698) = mat(k,1698) + .650_r8*rxt(k,297)*y(k,23) + .500_r8*rxt(k,305) &
                      *y(k,26) + .300_r8*rxt(k,283)*y(k,52) + .600_r8*rxt(k,407) &
                      *y(k,102) + .100_r8*rxt(k,360)*y(k,106) + .500_r8*rxt(k,340) &
                      *y(k,145) + .500_r8*rxt(k,415)*y(k,210)
         mat(k,1034) = .150_r8*rxt(k,342)*y(k,230)
         mat(k,2068) = rxt(k,253)*y(k,72) + 2.000_r8*rxt(k,145)*y(k,240)
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
         mat(k,391) = -(rxt(k,439)*y(k,230) + rxt(k,440)*y(k,122))
         mat(k,1451) = -rxt(k,439)*y(k,242)
         mat(k,1341) = -rxt(k,440)*y(k,242)
         mat(k,161) = .200_r8*rxt(k,429)*y(k,241)
         mat(k,134) = .140_r8*rxt(k,441)*y(k,241)
         mat(k,266) = rxt(k,444)*y(k,241)
         mat(k,1614) = .200_r8*rxt(k,429)*y(k,65) + .140_r8*rxt(k,441)*y(k,139) &
                      + rxt(k,444)*y(k,140)
         mat(k,674) = -(rxt(k,338)*y(k,230) + rxt(k,339)*y(k,122))
         mat(k,1475) = -rxt(k,338)*y(k,243)
         mat(k,1361) = -rxt(k,339)*y(k,243)
         mat(k,890) = rxt(k,345)*y(k,241)
         mat(k,460) = .500_r8*rxt(k,340)*y(k,241)
         mat(k,1648) = rxt(k,345)*y(k,28) + .500_r8*rxt(k,340)*y(k,145)
         mat(k,1029) = -(rxt(k,341)*y(k,225) + rxt(k,342)*y(k,230) + rxt(k,343) &
                      *y(k,122))
         mat(k,1851) = -rxt(k,341)*y(k,244)
         mat(k,1495) = -rxt(k,342)*y(k,244)
         mat(k,1380) = -rxt(k,343)*y(k,244)
         mat(k,839) = .060_r8*rxt(k,460)*y(k,132)
         mat(k,864) = rxt(k,346)*y(k,241)
         mat(k,813) = .060_r8*rxt(k,463)*y(k,132)
         mat(k,1910) = .060_r8*rxt(k,460)*y(k,5) + .060_r8*rxt(k,463)*y(k,110)
         mat(k,317) = rxt(k,344)*y(k,241)
         mat(k,948) = .150_r8*rxt(k,481)*y(k,241)
         mat(k,1677) = rxt(k,346)*y(k,47) + rxt(k,344)*y(k,146) + .150_r8*rxt(k,481) &
                      *y(k,207)
         mat(k,994) = -(rxt(k,470)*y(k,225) + rxt(k,471)*y(k,230) + rxt(k,472) &
                      *y(k,122))
         mat(k,1849) = -rxt(k,470)*y(k,245)
         mat(k,1493) = -rxt(k,471)*y(k,245)
         mat(k,1378) = -rxt(k,472)*y(k,245)
         mat(k,1797) = .500_r8*rxt(k,479)*y(k,206)
         mat(k,532) = rxt(k,473)*y(k,241)
         mat(k,884) = .500_r8*rxt(k,479)*y(k,124) + rxt(k,480)*y(k,241)
         mat(k,1675) = rxt(k,473)*y(k,203) + rxt(k,480)*y(k,206)
         mat(k,972) = -(rxt(k,475)*y(k,225) + rxt(k,476)*y(k,230) + rxt(k,477) &
                      *y(k,122))
         mat(k,1848) = -rxt(k,475)*y(k,246)
         mat(k,1492) = -rxt(k,476)*y(k,246)
         mat(k,1377) = -rxt(k,477)*y(k,246)
         mat(k,837) = rxt(k,461)*y(k,241)
         mat(k,811) = rxt(k,464)*y(k,241)
         mat(k,414) = rxt(k,478)*y(k,241)
         mat(k,1674) = rxt(k,461)*y(k,5) + rxt(k,464)*y(k,110) + rxt(k,478)*y(k,205)
         mat(k,638) = -(rxt(k,446)*y(k,230) + rxt(k,447)*y(k,122))
         mat(k,1472) = -rxt(k,446)*y(k,247)
         mat(k,1358) = -rxt(k,447)*y(k,247)
         mat(k,541) = rxt(k,448)*y(k,241)
         mat(k,157) = .650_r8*rxt(k,449)*y(k,241)
         mat(k,1645) = rxt(k,448)*y(k,208) + .650_r8*rxt(k,449)*y(k,209)
         mat(k,1068) = -(rxt(k,410)*y(k,224) + rxt(k,411)*y(k,225) + rxt(k,412) &
                      *y(k,230) + rxt(k,413)*y(k,122) + rxt(k,414)*y(k,124))
         mat(k,1271) = -rxt(k,410)*y(k,248)
         mat(k,1853) = -rxt(k,411)*y(k,248)
         mat(k,1498) = -rxt(k,412)*y(k,248)
         mat(k,1382) = -rxt(k,413)*y(k,248)
         mat(k,1801) = -rxt(k,414)*y(k,248)
         mat(k,205) = rxt(k,382)*y(k,241)
         mat(k,282) = rxt(k,383)*y(k,241)
         mat(k,104) = rxt(k,384)*y(k,241)
         mat(k,597) = .400_r8*rxt(k,407)*y(k,241)
         mat(k,170) = .500_r8*rxt(k,415)*y(k,241)
         mat(k,1680) = rxt(k,382)*y(k,93) + rxt(k,383)*y(k,95) + rxt(k,384)*y(k,96) &
                      + .400_r8*rxt(k,407)*y(k,102) + .500_r8*rxt(k,415)*y(k,210)
         mat(k,654) = -(rxt(k,452)*y(k,230) + rxt(k,453)*y(k,122))
         mat(k,1473) = -rxt(k,452)*y(k,249)
         mat(k,1359) = -rxt(k,453)*y(k,249)
         mat(k,181) = .560_r8*rxt(k,451)*y(k,241)
         mat(k,609) = rxt(k,454)*y(k,241)
         mat(k,1646) = .560_r8*rxt(k,451)*y(k,211) + rxt(k,454)*y(k,212)
         mat(k,428) = -(rxt(k,455)*y(k,230) + rxt(k,456)*y(k,122))
         mat(k,1456) = -rxt(k,455)*y(k,250)
         mat(k,1345) = -rxt(k,456)*y(k,250)
         mat(k,188) = .300_r8*rxt(k,457)*y(k,241)
         mat(k,371) = rxt(k,458)*y(k,241)
         mat(k,1620) = .300_r8*rxt(k,457)*y(k,213) + rxt(k,458)*y(k,214)
         mat(k,2078) = -(rxt(k,145)*y(k,240) + rxt(k,253)*y(k,72) + rxt(k,495) &
                      *y(k,151))
         mat(k,1549) = -rxt(k,145)*y(k,251)
         mat(k,700) = -rxt(k,253)*y(k,251)
         mat(k,215) = -rxt(k,495)*y(k,251)
         mat(k,248) = rxt(k,307)*y(k,241)
         mat(k,345) = rxt(k,332)*y(k,241)
         mat(k,102) = rxt(k,333)*y(k,241)
         mat(k,1731) = rxt(k,277)*y(k,241)
         mat(k,1048) = rxt(k,309)*y(k,241)
         mat(k,868) = rxt(k,346)*y(k,241)
         mat(k,1123) = rxt(k,335)*y(k,241)
         mat(k,470) = rxt(k,315)*y(k,241)
         mat(k,440) = rxt(k,316)*y(k,241)
         mat(k,357) = rxt(k,283)*y(k,241)
         mat(k,1331) = rxt(k,156)*y(k,230)
         mat(k,1061) = rxt(k,161)*y(k,241)
         mat(k,526) = rxt(k,162)*y(k,241)
         mat(k,727) = rxt(k,244)*y(k,241)
         mat(k,1431) = (rxt(k,538)+rxt(k,543))*y(k,90) + (rxt(k,531)+rxt(k,537) &
                       +rxt(k,542))*y(k,91) + rxt(k,215)*y(k,241)
         mat(k,708) = rxt(k,287)*y(k,241)
         mat(k,1317) = rxt(k,191)*y(k,241)
         mat(k,363) = rxt(k,169)*y(k,241)
         mat(k,691) = (rxt(k,538)+rxt(k,543))*y(k,84)
         mat(k,735) = (rxt(k,531)+rxt(k,537)+rxt(k,542))*y(k,84) + rxt(k,218)*y(k,241)
         mat(k,1114) = .500_r8*rxt(k,359)*y(k,241)
         mat(k,92) = rxt(k,498)*y(k,241)
         mat(k,466) = rxt(k,340)*y(k,241)
         mat(k,321) = rxt(k,344)*y(k,241)
         mat(k,1525) = rxt(k,156)*y(k,75) + rxt(k,163)*y(k,241)
         mat(k,1708) = rxt(k,307)*y(k,27) + rxt(k,332)*y(k,29) + rxt(k,333)*y(k,30) &
                      + rxt(k,277)*y(k,41) + rxt(k,309)*y(k,44) + rxt(k,346)*y(k,47) &
                      + rxt(k,335)*y(k,48) + rxt(k,315)*y(k,49) + rxt(k,316)*y(k,50) &
                      + rxt(k,283)*y(k,52) + rxt(k,161)*y(k,76) + rxt(k,162)*y(k,78) &
                      + rxt(k,244)*y(k,80) + rxt(k,215)*y(k,84) + rxt(k,287)*y(k,86) &
                      + rxt(k,191)*y(k,88) + rxt(k,169)*y(k,89) + rxt(k,218)*y(k,91) &
                      + .500_r8*rxt(k,359)*y(k,105) + rxt(k,498)*y(k,120) + rxt(k,340) &
                      *y(k,145) + rxt(k,344)*y(k,146) + rxt(k,163)*y(k,230) &
                      + 2.000_r8*rxt(k,166)*y(k,241)
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
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
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
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = mat(k, 99) + lmat(k, 99)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = lmat(k, 107)
         mat(k, 108) = lmat(k, 108)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 112) = lmat(k, 112)
         mat(k, 113) = lmat(k, 113)
         mat(k, 114) = lmat(k, 114)
         mat(k, 115) = lmat(k, 115)
         mat(k, 116) = lmat(k, 116)
         mat(k, 117) = lmat(k, 117)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 138) = lmat(k, 138)
         mat(k, 139) = lmat(k, 139)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 165) = lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = lmat(k, 209)
         mat(k, 210) = lmat(k, 210)
         mat(k, 212) = mat(k, 212) + lmat(k, 212)
         mat(k, 213) = lmat(k, 213)
         mat(k, 214) = lmat(k, 214)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 221) = lmat(k, 221)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 224) = lmat(k, 224)
         mat(k, 225) = lmat(k, 225)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = lmat(k, 227)
         mat(k, 228) = lmat(k, 228)
         mat(k, 229) = lmat(k, 229)
         mat(k, 230) = lmat(k, 230)
         mat(k, 231) = lmat(k, 231)
         mat(k, 232) = lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 255) = lmat(k, 255)
         mat(k, 256) = lmat(k, 256)
         mat(k, 257) = lmat(k, 257)
         mat(k, 258) = lmat(k, 258)
         mat(k, 259) = lmat(k, 259)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 264) = lmat(k, 264)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 267) = lmat(k, 267)
         mat(k, 268) = lmat(k, 268)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 278) = mat(k, 278) + lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 286) = lmat(k, 286)
         mat(k, 287) = lmat(k, 287)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 289) = mat(k, 289) + lmat(k, 289)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 303) = lmat(k, 303)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 313) = lmat(k, 313)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 322) = mat(k, 322) + lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 327) = lmat(k, 327)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 329) = lmat(k, 329)
         mat(k, 331) = lmat(k, 331)
         mat(k, 332) = mat(k, 332) + lmat(k, 332)
         mat(k, 333) = lmat(k, 333)
         mat(k, 334) = lmat(k, 334)
         mat(k, 335) = lmat(k, 335)
         mat(k, 336) = lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = lmat(k, 339)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 347) = lmat(k, 347)
         mat(k, 349) = lmat(k, 349)
         mat(k, 350) = lmat(k, 350)
         mat(k, 351) = mat(k, 351) + lmat(k, 351)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 354) = mat(k, 354) + lmat(k, 354)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 359) = lmat(k, 359)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 362) = lmat(k, 362)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 372) = lmat(k, 372)
         mat(k, 373) = lmat(k, 373)
         mat(k, 374) = lmat(k, 374)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 378) = mat(k, 378) + lmat(k, 378)
         mat(k, 384) = mat(k, 384) + lmat(k, 384)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 388) = lmat(k, 388)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 397) = lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 399) = lmat(k, 399)
         mat(k, 400) = mat(k, 400) + lmat(k, 400)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 405) = lmat(k, 405)
         mat(k, 408) = lmat(k, 408)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 410) = lmat(k, 410)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 412) = lmat(k, 412)
         mat(k, 413) = lmat(k, 413)
         mat(k, 415) = lmat(k, 415)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 437) = mat(k, 437) + lmat(k, 437)
         mat(k, 439) = lmat(k, 439)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 446) = lmat(k, 446)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 449) = lmat(k, 449)
         mat(k, 450) = lmat(k, 450)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = lmat(k, 463)
         mat(k, 464) = mat(k, 464) + lmat(k, 464)
         mat(k, 465) = lmat(k, 465)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 471) = mat(k, 471) + lmat(k, 471)
         mat(k, 476) = lmat(k, 476)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 481) = lmat(k, 481)
         mat(k, 485) = lmat(k, 485)
         mat(k, 487) = mat(k, 487) + lmat(k, 487)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 503) = mat(k, 503) + lmat(k, 503)
         mat(k, 504) = lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 515) = lmat(k, 515)
         mat(k, 520) = mat(k, 520) + lmat(k, 520)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 528) = lmat(k, 528)
         mat(k, 529) = lmat(k, 529)
         mat(k, 530) = lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 533) = lmat(k, 533)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 535) = lmat(k, 535)
         mat(k, 536) = lmat(k, 536)
         mat(k, 537) = lmat(k, 537)
         mat(k, 538) = mat(k, 538) + lmat(k, 538)
         mat(k, 539) = lmat(k, 539)
         mat(k, 543) = lmat(k, 543)
         mat(k, 544) = lmat(k, 544)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 547) = lmat(k, 547)
         mat(k, 548) = lmat(k, 548)
         mat(k, 549) = lmat(k, 549)
         mat(k, 550) = lmat(k, 550)
         mat(k, 551) = mat(k, 551) + lmat(k, 551)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 553) = lmat(k, 553)
         mat(k, 554) = lmat(k, 554)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 559) = mat(k, 559) + lmat(k, 559)
         mat(k, 565) = lmat(k, 565)
         mat(k, 566) = mat(k, 566) + lmat(k, 566)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 570) = mat(k, 570) + lmat(k, 570)
         mat(k, 572) = lmat(k, 572)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 576) = lmat(k, 576)
         mat(k, 577) = mat(k, 577) + lmat(k, 577)
         mat(k, 581) = lmat(k, 581)
         mat(k, 582) = lmat(k, 582)
         mat(k, 584) = lmat(k, 584)
         mat(k, 585) = mat(k, 585) + lmat(k, 585)
         mat(k, 586) = lmat(k, 586)
         mat(k, 588) = mat(k, 588) + lmat(k, 588)
         mat(k, 596) = mat(k, 596) + lmat(k, 596)
         mat(k, 598) = lmat(k, 598)
         mat(k, 599) = lmat(k, 599)
         mat(k, 600) = lmat(k, 600)
         mat(k, 602) = lmat(k, 602)
         mat(k, 603) = lmat(k, 603)
         mat(k, 604) = lmat(k, 604)
         mat(k, 605) = lmat(k, 605)
         mat(k, 606) = lmat(k, 606)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 612) = lmat(k, 612)
         mat(k, 614) = lmat(k, 614)
         mat(k, 616) = lmat(k, 616)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 638) = mat(k, 638) + lmat(k, 638)
         mat(k, 654) = mat(k, 654) + lmat(k, 654)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 674) = mat(k, 674) + lmat(k, 674)
         mat(k, 684) = mat(k, 684) + lmat(k, 684)
         mat(k, 685) = lmat(k, 685)
         mat(k, 687) = mat(k, 687) + lmat(k, 687)
         mat(k, 692) = mat(k, 692) + lmat(k, 692)
         mat(k, 693) = mat(k, 693) + lmat(k, 693)
         mat(k, 699) = lmat(k, 699)
         mat(k, 701) = mat(k, 701) + lmat(k, 701)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 710) = mat(k, 710) + lmat(k, 710)
         mat(k, 720) = mat(k, 720) + lmat(k, 720)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 722) = mat(k, 722) + lmat(k, 722)
         mat(k, 729) = mat(k, 729) + lmat(k, 729)
         mat(k, 731) = mat(k, 731) + lmat(k, 731)
         mat(k, 734) = mat(k, 734) + lmat(k, 734)
         mat(k, 739) = mat(k, 739) + lmat(k, 739)
         mat(k, 747) = mat(k, 747) + lmat(k, 747)
         mat(k, 749) = lmat(k, 749)
         mat(k, 751) = lmat(k, 751)
         mat(k, 752) = mat(k, 752) + lmat(k, 752)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 755) = lmat(k, 755)
         mat(k, 758) = lmat(k, 758)
         mat(k, 766) = mat(k, 766) + lmat(k, 766)
         mat(k, 782) = mat(k, 782) + lmat(k, 782)
         mat(k, 791) = lmat(k, 791)
         mat(k, 792) = mat(k, 792) + lmat(k, 792)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 794) = mat(k, 794) + lmat(k, 794)
         mat(k, 805) = mat(k, 805) + lmat(k, 805)
         mat(k, 831) = mat(k, 831) + lmat(k, 831)
         mat(k, 853) = mat(k, 853) + lmat(k, 853)
         mat(k, 863) = mat(k, 863) + lmat(k, 863)
         mat(k, 865) = lmat(k, 865)
         mat(k, 867) = lmat(k, 867)
         mat(k, 870) = mat(k, 870) + lmat(k, 870)
         mat(k, 871) = mat(k, 871) + lmat(k, 871)
         mat(k, 872) = mat(k, 872) + lmat(k, 872)
         mat(k, 875) = lmat(k, 875)
         mat(k, 876) = mat(k, 876) + lmat(k, 876)
         mat(k, 877) = mat(k, 877) + lmat(k, 877)
         mat(k, 879) = mat(k, 879) + lmat(k, 879)
         mat(k, 881) = mat(k, 881) + lmat(k, 881)
         mat(k, 882) = lmat(k, 882)
         mat(k, 883) = lmat(k, 883)
         mat(k, 886) = lmat(k, 886)
         mat(k, 893) = mat(k, 893) + lmat(k, 893)
         mat(k, 908) = lmat(k, 908)
         mat(k, 912) = mat(k, 912) + lmat(k, 912)
         mat(k, 916) = lmat(k, 916)
         mat(k, 918) = mat(k, 918) + lmat(k, 918)
         mat(k, 921) = lmat(k, 921)
         mat(k, 924) = mat(k, 924) + lmat(k, 924)
         mat(k, 925) = lmat(k, 925)
         mat(k, 926) = mat(k, 926) + lmat(k, 926)
         mat(k, 928) = mat(k, 928) + lmat(k, 928)
         mat(k, 930) = lmat(k, 930)
         mat(k, 931) = lmat(k, 931)
         mat(k, 932) = mat(k, 932) + lmat(k, 932)
         mat(k, 933) = lmat(k, 933)
         mat(k, 934) = lmat(k, 934)
         mat(k, 936) = lmat(k, 936)
         mat(k, 937) = lmat(k, 937)
         mat(k, 938) = lmat(k, 938)
         mat(k, 939) = mat(k, 939) + lmat(k, 939)
         mat(k, 941) = lmat(k, 941)
         mat(k, 942) = lmat(k, 942)
         mat(k, 945) = mat(k, 945) + lmat(k, 945)
         mat(k, 946) = mat(k, 946) + lmat(k, 946)
         mat(k, 947) = mat(k, 947) + lmat(k, 947)
         mat(k, 948) = mat(k, 948) + lmat(k, 948)
         mat(k, 949) = mat(k, 949) + lmat(k, 949)
         mat(k, 950) = mat(k, 950) + lmat(k, 950)
         mat(k, 952) = mat(k, 952) + lmat(k, 952)
         mat(k, 955) = mat(k, 955) + lmat(k, 955)
         mat(k, 959) = mat(k, 959) + lmat(k, 959)
         mat(k, 961) = lmat(k, 961)
         mat(k, 962) = mat(k, 962) + lmat(k, 962)
         mat(k, 964) = lmat(k, 964)
         mat(k, 972) = mat(k, 972) + lmat(k, 972)
         mat(k, 994) = mat(k, 994) + lmat(k, 994)
         mat(k,1013) = mat(k,1013) + lmat(k,1013)
         mat(k,1029) = mat(k,1029) + lmat(k,1029)
         mat(k,1039) = lmat(k,1039)
         mat(k,1040) = mat(k,1040) + lmat(k,1040)
         mat(k,1044) = lmat(k,1044)
         mat(k,1047) = lmat(k,1047)
         mat(k,1051) = mat(k,1051) + lmat(k,1051)
         mat(k,1068) = mat(k,1068) + lmat(k,1068)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1107) = mat(k,1107) + lmat(k,1107)
         mat(k,1108) = mat(k,1108) + lmat(k,1108)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1112) = mat(k,1112) + lmat(k,1112)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1116) = mat(k,1116) + lmat(k,1116)
         mat(k,1117) = mat(k,1117) + lmat(k,1117)
         mat(k,1120) = lmat(k,1120)
         mat(k,1125) = lmat(k,1125)
         mat(k,1126) = mat(k,1126) + lmat(k,1126)
         mat(k,1127) = mat(k,1127) + lmat(k,1127)
         mat(k,1137) = lmat(k,1137)
         mat(k,1152) = mat(k,1152) + lmat(k,1152)
         mat(k,1169) = lmat(k,1169)
         mat(k,1170) = mat(k,1170) + lmat(k,1170)
         mat(k,1174) = mat(k,1174) + lmat(k,1174)
         mat(k,1176) = mat(k,1176) + lmat(k,1176)
         mat(k,1185) = lmat(k,1185)
         mat(k,1197) = mat(k,1197) + lmat(k,1197)
         mat(k,1210) = lmat(k,1210)
         mat(k,1230) = mat(k,1230) + lmat(k,1230)
         mat(k,1236) = mat(k,1236) + lmat(k,1236)
         mat(k,1249) = mat(k,1249) + lmat(k,1249)
         mat(k,1280) = mat(k,1280) + lmat(k,1280)
         mat(k,1294) = mat(k,1294) + lmat(k,1294)
         mat(k,1307) = mat(k,1307) + lmat(k,1307)
         mat(k,1311) = mat(k,1311) + lmat(k,1311)
         mat(k,1312) = lmat(k,1312)
         mat(k,1320) = mat(k,1320) + lmat(k,1320)
         mat(k,1322) = mat(k,1322) + lmat(k,1322)
         mat(k,1338) = mat(k,1338) + lmat(k,1338)
         mat(k,1395) = mat(k,1395) + lmat(k,1395)
         mat(k,1407) = mat(k,1407) + lmat(k,1407)
         mat(k,1417) = mat(k,1417) + lmat(k,1417)
         mat(k,1418) = mat(k,1418) + lmat(k,1418)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1461) = mat(k,1461) + lmat(k,1461)
         mat(k,1513) = mat(k,1513) + lmat(k,1513)
         mat(k,1526) = mat(k,1526) + lmat(k,1526)
         mat(k,1527) = mat(k,1527) + lmat(k,1527)
         mat(k,1529) = mat(k,1529) + lmat(k,1529)
         mat(k,1531) = mat(k,1531) + lmat(k,1531)
         mat(k,1532) = mat(k,1532) + lmat(k,1532)
         mat(k,1534) = mat(k,1534) + lmat(k,1534)
         mat(k,1535) = lmat(k,1535)
         mat(k,1537) = lmat(k,1537)
         mat(k,1538) = mat(k,1538) + lmat(k,1538)
         mat(k,1539) = mat(k,1539) + lmat(k,1539)
         mat(k,1540) = lmat(k,1540)
         mat(k,1543) = lmat(k,1543)
         mat(k,1547) = lmat(k,1547)
         mat(k,1548) = mat(k,1548) + lmat(k,1548)
         mat(k,1573) = lmat(k,1573)
         mat(k,1577) = lmat(k,1577)
         mat(k,1691) = mat(k,1691) + lmat(k,1691)
         mat(k,1696) = mat(k,1696) + lmat(k,1696)
         mat(k,1698) = mat(k,1698) + lmat(k,1698)
         mat(k,1702) = mat(k,1702) + lmat(k,1702)
         mat(k,1707) = mat(k,1707) + lmat(k,1707)
         mat(k,1708) = mat(k,1708) + lmat(k,1708)
         mat(k,1712) = mat(k,1712) + lmat(k,1712)
         mat(k,1713) = lmat(k,1713)
         mat(k,1716) = mat(k,1716) + lmat(k,1716)
         mat(k,1722) = mat(k,1722) + lmat(k,1722)
         mat(k,1756) = mat(k,1756) + lmat(k,1756)
         mat(k,1758) = mat(k,1758) + lmat(k,1758)
         mat(k,1762) = mat(k,1762) + lmat(k,1762)
         mat(k,1764) = mat(k,1764) + lmat(k,1764)
         mat(k,1770) = mat(k,1770) + lmat(k,1770)
         mat(k,1813) = mat(k,1813) + lmat(k,1813)
         mat(k,1815) = mat(k,1815) + lmat(k,1815)
         mat(k,1821) = mat(k,1821) + lmat(k,1821)
         mat(k,1822) = mat(k,1822) + lmat(k,1822)
         mat(k,1827) = mat(k,1827) + lmat(k,1827)
         mat(k,1873) = mat(k,1873) + lmat(k,1873)
         mat(k,1928) = mat(k,1928) + lmat(k,1928)
         mat(k,1934) = mat(k,1934) + lmat(k,1934)
         mat(k,1937) = mat(k,1937) + lmat(k,1937)
         mat(k,1946) = mat(k,1946) + lmat(k,1946)
         mat(k,1959) = mat(k,1959) + lmat(k,1959)
         mat(k,1961) = mat(k,1961) + lmat(k,1961)
         mat(k,1986) = mat(k,1986) + lmat(k,1986)
         mat(k,1987) = mat(k,1987) + lmat(k,1987)
         mat(k,1988) = mat(k,1988) + lmat(k,1988)
         mat(k,2014) = mat(k,2014) + lmat(k,2014)
         mat(k,2017) = mat(k,2017) + lmat(k,2017)
         mat(k,2032) = mat(k,2032) + lmat(k,2032)
         mat(k,2036) = lmat(k,2036)
         mat(k,2040) = mat(k,2040) + lmat(k,2040)
         mat(k,2041) = mat(k,2041) + lmat(k,2041)
         mat(k,2047) = lmat(k,2047)
         mat(k,2052) = mat(k,2052) + lmat(k,2052)
         mat(k,2059) = lmat(k,2059)
         mat(k,2063) = lmat(k,2063)
         mat(k,2067) = mat(k,2067) + lmat(k,2067)
         mat(k,2068) = mat(k,2068) + lmat(k,2068)
         mat(k,2076) = lmat(k,2076)
         mat(k,2078) = mat(k,2078) + lmat(k,2078)
         mat(k, 182) = 0._r8
         mat(k, 183) = 0._r8
         mat(k, 273) = 0._r8
         mat(k, 379) = 0._r8
         mat(k, 382) = 0._r8
         mat(k, 395) = 0._r8
         mat(k, 421) = 0._r8
         mat(k, 425) = 0._r8
         mat(k, 433) = 0._r8
         mat(k, 540) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 564) = 0._r8
         mat(k, 567) = 0._r8
         mat(k, 568) = 0._r8
         mat(k, 571) = 0._r8
         mat(k, 578) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 583) = 0._r8
         mat(k, 608) = 0._r8
         mat(k, 610) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 615) = 0._r8
         mat(k, 621) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 646) = 0._r8
         mat(k, 653) = 0._r8
         mat(k, 655) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 658) = 0._r8
         mat(k, 660) = 0._r8
         mat(k, 663) = 0._r8
         mat(k, 675) = 0._r8
         mat(k, 679) = 0._r8
         mat(k, 682) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 715) = 0._r8
         mat(k, 719) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 759) = 0._r8
         mat(k, 781) = 0._r8
         mat(k, 804) = 0._r8
         mat(k, 806) = 0._r8
         mat(k, 814) = 0._r8
         mat(k, 821) = 0._r8
         mat(k, 830) = 0._r8
         mat(k, 832) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 847) = 0._r8
         mat(k, 851) = 0._r8
         mat(k, 852) = 0._r8
         mat(k, 856) = 0._r8
         mat(k, 859) = 0._r8
         mat(k, 862) = 0._r8
         mat(k, 880) = 0._r8
         mat(k, 895) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 903) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 913) = 0._r8
         mat(k, 914) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 920) = 0._r8
         mat(k, 922) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 943) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 970) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 973) = 0._r8
         mat(k, 974) = 0._r8
         mat(k, 975) = 0._r8
         mat(k, 976) = 0._r8
         mat(k, 977) = 0._r8
         mat(k, 978) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 986) = 0._r8
         mat(k, 995) = 0._r8
         mat(k, 996) = 0._r8
         mat(k, 997) = 0._r8
         mat(k, 998) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1005) = 0._r8
         mat(k,1010) = 0._r8
         mat(k,1011) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1016) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1020) = 0._r8
         mat(k,1025) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1043) = 0._r8
         mat(k,1052) = 0._r8
         mat(k,1057) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1081) = 0._r8
         mat(k,1083) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1087) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1097) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1131) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1171) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1213) = 0._r8
         mat(k,1216) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1220) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1224) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1233) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1248) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1288) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1327) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1394) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1486) = 0._r8
         mat(k,1487) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1503) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1615) = 0._r8
         mat(k,1635) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1649) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1681) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1710) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1720) = 0._r8
         mat(k,1723) = 0._r8
         mat(k,1725) = 0._r8
         mat(k,1726) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1741) = 0._r8
         mat(k,1744) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1749) = 0._r8
         mat(k,1750) = 0._r8
         mat(k,1751) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1759) = 0._r8
         mat(k,1761) = 0._r8
         mat(k,1763) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1772) = 0._r8
         mat(k,1779) = 0._r8
         mat(k,1786) = 0._r8
         mat(k,1792) = 0._r8
         mat(k,1793) = 0._r8
         mat(k,1796) = 0._r8
         mat(k,1799) = 0._r8
         mat(k,1805) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1814) = 0._r8
         mat(k,1816) = 0._r8
         mat(k,1818) = 0._r8
         mat(k,1823) = 0._r8
         mat(k,1824) = 0._r8
         mat(k,1825) = 0._r8
         mat(k,1826) = 0._r8
         mat(k,1828) = 0._r8
         mat(k,1829) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1863) = 0._r8
         mat(k,1864) = 0._r8
         mat(k,1866) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1869) = 0._r8
         mat(k,1872) = 0._r8
         mat(k,1874) = 0._r8
         mat(k,1875) = 0._r8
         mat(k,1877) = 0._r8
         mat(k,1879) = 0._r8
         mat(k,1892) = 0._r8
         mat(k,1896) = 0._r8
         mat(k,1899) = 0._r8
         mat(k,1903) = 0._r8
         mat(k,1906) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1908) = 0._r8
         mat(k,1909) = 0._r8
         mat(k,1912) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1919) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1926) = 0._r8
         mat(k,1939) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1956) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1963) = 0._r8
         mat(k,1973) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1984) = 0._r8
         mat(k,1989) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1997) = 0._r8
         mat(k,2003) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2013) = 0._r8
         mat(k,2019) = 0._r8
         mat(k,2025) = 0._r8
         mat(k,2026) = 0._r8
         mat(k,2029) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2037) = 0._r8
         mat(k,2039) = 0._r8
         mat(k,2042) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2049) = 0._r8
         mat(k,2051) = 0._r8
         mat(k,2053) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2060) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2064) = 0._r8
         mat(k,2065) = 0._r8
         mat(k,2066) = 0._r8
         mat(k,2069) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2073) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2075) = 0._r8
         mat(k,2077) = 0._r8
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
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 89) = mat(k, 89) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 99) = mat(k, 99) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 106) = mat(k, 106) - dti(k)
         mat(k, 109) = mat(k, 109) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
         mat(k, 123) = mat(k, 123) - dti(k)
         mat(k, 129) = mat(k, 129) - dti(k)
         mat(k, 133) = mat(k, 133) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 140) = mat(k, 140) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 153) = mat(k, 153) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 165) = mat(k, 165) - dti(k)
         mat(k, 169) = mat(k, 169) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 186) = mat(k, 186) - dti(k)
         mat(k, 191) = mat(k, 191) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 197) = mat(k, 197) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 203) = mat(k, 203) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 212) = mat(k, 212) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 220) = mat(k, 220) - dti(k)
         mat(k, 224) = mat(k, 224) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 231) = mat(k, 231) - dti(k)
         mat(k, 237) = mat(k, 237) - dti(k)
         mat(k, 243) = mat(k, 243) - dti(k)
         mat(k, 249) = mat(k, 249) - dti(k)
         mat(k, 255) = mat(k, 255) - dti(k)
         mat(k, 260) = mat(k, 260) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 281) = mat(k, 281) - dti(k)
         mat(k, 284) = mat(k, 284) - dti(k)
         mat(k, 289) = mat(k, 289) - dti(k)
         mat(k, 294) = mat(k, 294) - dti(k)
         mat(k, 302) = mat(k, 302) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 322) = mat(k, 322) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 346) = mat(k, 346) - dti(k)
         mat(k, 352) = mat(k, 352) - dti(k)
         mat(k, 358) = mat(k, 358) - dti(k)
         mat(k, 364) = mat(k, 364) - dti(k)
         mat(k, 370) = mat(k, 370) - dti(k)
         mat(k, 378) = mat(k, 378) - dti(k)
         mat(k, 384) = mat(k, 384) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 397) = mat(k, 397) - dti(k)
         mat(k, 400) = mat(k, 400) - dti(k)
         mat(k, 404) = mat(k, 404) - dti(k)
         mat(k, 411) = mat(k, 411) - dti(k)
         mat(k, 420) = mat(k, 420) - dti(k)
         mat(k, 428) = mat(k, 428) - dti(k)
         mat(k, 435) = mat(k, 435) - dti(k)
         mat(k, 441) = mat(k, 441) - dti(k)
         mat(k, 447) = mat(k, 447) - dti(k)
         mat(k, 453) = mat(k, 453) - dti(k)
         mat(k, 459) = mat(k, 459) - dti(k)
         mat(k, 467) = mat(k, 467) - dti(k)
         mat(k, 471) = mat(k, 471) - dti(k)
         mat(k, 479) = mat(k, 479) - dti(k)
         mat(k, 487) = mat(k, 487) - dti(k)
         mat(k, 495) = mat(k, 495) - dti(k)
         mat(k, 503) = mat(k, 503) - dti(k)
         mat(k, 511) = mat(k, 511) - dti(k)
         mat(k, 520) = mat(k, 520) - dti(k)
         mat(k, 527) = mat(k, 527) - dti(k)
         mat(k, 538) = mat(k, 538) - dti(k)
         mat(k, 547) = mat(k, 547) - dti(k)
         mat(k, 551) = mat(k, 551) - dti(k)
         mat(k, 559) = mat(k, 559) - dti(k)
         mat(k, 566) = mat(k, 566) - dti(k)
         mat(k, 577) = mat(k, 577) - dti(k)
         mat(k, 588) = mat(k, 588) - dti(k)
         mat(k, 596) = mat(k, 596) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 627) = mat(k, 627) - dti(k)
         mat(k, 638) = mat(k, 638) - dti(k)
         mat(k, 654) = mat(k, 654) - dti(k)
         mat(k, 665) = mat(k, 665) - dti(k)
         mat(k, 674) = mat(k, 674) - dti(k)
         mat(k, 684) = mat(k, 684) - dti(k)
         mat(k, 693) = mat(k, 693) - dti(k)
         mat(k, 701) = mat(k, 701) - dti(k)
         mat(k, 705) = mat(k, 705) - dti(k)
         mat(k, 710) = mat(k, 710) - dti(k)
         mat(k, 720) = mat(k, 720) - dti(k)
         mat(k, 729) = mat(k, 729) - dti(k)
         mat(k, 739) = mat(k, 739) - dti(k)
         mat(k, 747) = mat(k, 747) - dti(k)
         mat(k, 754) = mat(k, 754) - dti(k)
         mat(k, 766) = mat(k, 766) - dti(k)
         mat(k, 782) = mat(k, 782) - dti(k)
         mat(k, 792) = mat(k, 792) - dti(k)
         mat(k, 805) = mat(k, 805) - dti(k)
         mat(k, 831) = mat(k, 831) - dti(k)
         mat(k, 853) = mat(k, 853) - dti(k)
         mat(k, 863) = mat(k, 863) - dti(k)
         mat(k, 871) = mat(k, 871) - dti(k)
         mat(k, 881) = mat(k, 881) - dti(k)
         mat(k, 893) = mat(k, 893) - dti(k)
         mat(k, 912) = mat(k, 912) - dti(k)
         mat(k, 924) = mat(k, 924) - dti(k)
         mat(k, 932) = mat(k, 932) - dti(k)
         mat(k, 946) = mat(k, 946) - dti(k)
         mat(k, 955) = mat(k, 955) - dti(k)
         mat(k, 959) = mat(k, 959) - dti(k)
         mat(k, 972) = mat(k, 972) - dti(k)
         mat(k, 994) = mat(k, 994) - dti(k)
         mat(k,1013) = mat(k,1013) - dti(k)
         mat(k,1029) = mat(k,1029) - dti(k)
         mat(k,1040) = mat(k,1040) - dti(k)
         mat(k,1051) = mat(k,1051) - dti(k)
         mat(k,1068) = mat(k,1068) - dti(k)
         mat(k,1088) = mat(k,1088) - dti(k)
         mat(k,1104) = mat(k,1104) - dti(k)
         mat(k,1116) = mat(k,1116) - dti(k)
         mat(k,1127) = mat(k,1127) - dti(k)
         mat(k,1152) = mat(k,1152) - dti(k)
         mat(k,1174) = mat(k,1174) - dti(k)
         mat(k,1197) = mat(k,1197) - dti(k)
         mat(k,1230) = mat(k,1230) - dti(k)
         mat(k,1249) = mat(k,1249) - dti(k)
         mat(k,1280) = mat(k,1280) - dti(k)
         mat(k,1294) = mat(k,1294) - dti(k)
         mat(k,1307) = mat(k,1307) - dti(k)
         mat(k,1320) = mat(k,1320) - dti(k)
         mat(k,1395) = mat(k,1395) - dti(k)
         mat(k,1418) = mat(k,1418) - dti(k)
         mat(k,1513) = mat(k,1513) - dti(k)
         mat(k,1538) = mat(k,1538) - dti(k)
         mat(k,1698) = mat(k,1698) - dti(k)
         mat(k,1722) = mat(k,1722) - dti(k)
         mat(k,1764) = mat(k,1764) - dti(k)
         mat(k,1822) = mat(k,1822) - dti(k)
         mat(k,1873) = mat(k,1873) - dti(k)
         mat(k,1934) = mat(k,1934) - dti(k)
         mat(k,1959) = mat(k,1959) - dti(k)
         mat(k,1986) = mat(k,1986) - dti(k)
         mat(k,2017) = mat(k,2017) - dti(k)
         mat(k,2052) = mat(k,2052) - dti(k)
         mat(k,2078) = mat(k,2078) - dti(k)
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
