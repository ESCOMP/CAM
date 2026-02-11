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
         mat(k,661) = -(rxt(k,368)*y(k,221))
         mat(k,1721) = -rxt(k,368)*y(k,1)
         mat(k,2288) = rxt(k,371)*y(k,193)
         mat(k,957) = rxt(k,371)*y(k,123)
         mat(k,627) = -(rxt(k,372)*y(k,221))
         mat(k,1718) = -rxt(k,372)*y(k,2)
         mat(k,956) = rxt(k,369)*y(k,207)
         mat(k,1864) = rxt(k,369)*y(k,193)
         mat(k,936) = -(rxt(k,451)*y(k,125) + rxt(k,452)*y(k,137) + rxt(k,453) &
                      *y(k,221))
         mat(k,2085) = -rxt(k,451)*y(k,5)
         mat(k,1949) = -rxt(k,452)*y(k,5)
         mat(k,1746) = -rxt(k,453)*y(k,5)
         mat(k,168) = -(rxt(k,410)*y(k,221))
         mat(k,1651) = -rxt(k,410)*y(k,6)
         mat(k,412) = -(rxt(k,413)*y(k,221))
         mat(k,1690) = -rxt(k,413)*y(k,7)
         mat(k,490) = rxt(k,411)*y(k,207)
         mat(k,1846) = rxt(k,411)*y(k,195)
         mat(k,169) = .120_r8*rxt(k,410)*y(k,221)
         mat(k,1652) = .120_r8*rxt(k,410)*y(k,6)
         mat(k,934) = .100_r8*rxt(k,452)*y(k,137)
         mat(k,906) = .100_r8*rxt(k,455)*y(k,137)
         mat(k,1938) = .100_r8*rxt(k,452)*y(k,5) + .100_r8*rxt(k,455)*y(k,109)
         mat(k,2276) = .500_r8*rxt(k,412)*y(k,195) + .200_r8*rxt(k,439)*y(k,227) &
                      + .060_r8*rxt(k,445)*y(k,230)
         mat(k,491) = .500_r8*rxt(k,412)*y(k,123)
         mat(k,730) = .200_r8*rxt(k,439)*y(k,123)
         mat(k,746) = .060_r8*rxt(k,445)*y(k,123)
         mat(k,2269) = .200_r8*rxt(k,439)*y(k,227) + .200_r8*rxt(k,445)*y(k,230)
         mat(k,729) = .200_r8*rxt(k,439)*y(k,123)
         mat(k,744) = .200_r8*rxt(k,445)*y(k,123)
         mat(k,2285) = .200_r8*rxt(k,439)*y(k,227) + .150_r8*rxt(k,445)*y(k,230)
         mat(k,731) = .200_r8*rxt(k,439)*y(k,123)
         mat(k,747) = .150_r8*rxt(k,445)*y(k,123)
         mat(k,2270) = .210_r8*rxt(k,445)*y(k,230)
         mat(k,745) = .210_r8*rxt(k,445)*y(k,123)
         mat(k,229) = -(rxt(k,373)*y(k,221))
         mat(k,1662) = -rxt(k,373)*y(k,14)
         mat(k,933) = .050_r8*rxt(k,452)*y(k,137)
         mat(k,905) = .050_r8*rxt(k,455)*y(k,137)
         mat(k,1937) = .050_r8*rxt(k,452)*y(k,5) + .050_r8*rxt(k,455)*y(k,109)
         mat(k,352) = -(rxt(k,339)*y(k,125) + rxt(k,340)*y(k,221))
         mat(k,2078) = -rxt(k,339)*y(k,15)
         mat(k,1682) = -rxt(k,340)*y(k,15)
         mat(k,1446) = -(rxt(k,223)*y(k,41) + rxt(k,224)*y(k,207) + rxt(k,225) &
                      *y(k,137))
         mat(k,1524) = -rxt(k,223)*y(k,16)
         mat(k,1910) = -rxt(k,224)*y(k,16)
         mat(k,1975) = -rxt(k,225)*y(k,16)
         mat(k,1803) = 4.000_r8*rxt(k,226)*y(k,18) + (rxt(k,227)+rxt(k,228))*y(k,58) &
                      + rxt(k,231)*y(k,123) + rxt(k,234)*y(k,133) + rxt(k,481) &
                      *y(k,153) + rxt(k,235)*y(k,221)
         mat(k,144) = rxt(k,213)*y(k,220)
         mat(k,150) = rxt(k,239)*y(k,220)
         mat(k,477) = 2.000_r8*rxt(k,250)*y(k,55) + 2.000_r8*rxt(k,262)*y(k,220) &
                      + 2.000_r8*rxt(k,251)*y(k,221)
         mat(k,590) = rxt(k,252)*y(k,55) + rxt(k,263)*y(k,220) + rxt(k,253)*y(k,221)
         mat(k,383) = 3.000_r8*rxt(k,257)*y(k,55) + 3.000_r8*rxt(k,240)*y(k,220) &
                      + 3.000_r8*rxt(k,258)*y(k,221)
         mat(k,2208) = 2.000_r8*rxt(k,250)*y(k,40) + rxt(k,252)*y(k,42) &
                      + 3.000_r8*rxt(k,257)*y(k,54)
         mat(k,2236) = (rxt(k,227)+rxt(k,228))*y(k,18)
         mat(k,108) = 2.000_r8*rxt(k,241)*y(k,220)
         mat(k,816) = rxt(k,236)*y(k,133) + rxt(k,242)*y(k,220) + rxt(k,237)*y(k,221)
         mat(k,2328) = rxt(k,231)*y(k,18)
         mat(k,2009) = rxt(k,234)*y(k,18) + rxt(k,236)*y(k,80)
         mat(k,1431) = rxt(k,481)*y(k,18)
         mat(k,1612) = rxt(k,213)*y(k,33) + rxt(k,239)*y(k,34) + 2.000_r8*rxt(k,262) &
                      *y(k,40) + rxt(k,263)*y(k,42) + 3.000_r8*rxt(k,240)*y(k,54) &
                      + 2.000_r8*rxt(k,241)*y(k,77) + rxt(k,242)*y(k,80)
         mat(k,1777) = rxt(k,235)*y(k,18) + 2.000_r8*rxt(k,251)*y(k,40) + rxt(k,253) &
                      *y(k,42) + 3.000_r8*rxt(k,258)*y(k,54) + rxt(k,237)*y(k,80)
         mat(k,1797) = rxt(k,229)*y(k,58)
         mat(k,2230) = rxt(k,229)*y(k,18)
         mat(k,2050) = (rxt(k,542)+rxt(k,547))*y(k,90)
         mat(k,769) = (rxt(k,542)+rxt(k,547))*y(k,84)
         mat(k,1810) = -(4._r8*rxt(k,226)*y(k,18) + (rxt(k,227) + rxt(k,228) + rxt(k,229) &
                      ) * y(k,58) + rxt(k,230)*y(k,207) + rxt(k,231)*y(k,123) &
                      + rxt(k,232)*y(k,124) + rxt(k,234)*y(k,133) + rxt(k,235) &
                      *y(k,221) + rxt(k,481)*y(k,153))
         mat(k,2243) = -(rxt(k,227) + rxt(k,228) + rxt(k,229)) * y(k,18)
         mat(k,1918) = -rxt(k,230)*y(k,18)
         mat(k,2336) = -rxt(k,231)*y(k,18)
         mat(k,1577) = -rxt(k,232)*y(k,18)
         mat(k,2017) = -rxt(k,234)*y(k,18)
         mat(k,1785) = -rxt(k,235)*y(k,18)
         mat(k,1435) = -rxt(k,481)*y(k,18)
         mat(k,1451) = rxt(k,225)*y(k,137)
         mat(k,548) = rxt(k,233)*y(k,133)
         mat(k,819) = rxt(k,243)*y(k,220)
         mat(k,773) = rxt(k,238)*y(k,133)
         mat(k,2017) = mat(k,2017) + rxt(k,233)*y(k,19) + rxt(k,238)*y(k,90)
         mat(k,1983) = rxt(k,225)*y(k,16)
         mat(k,1620) = rxt(k,243)*y(k,80)
         mat(k,543) = -(rxt(k,233)*y(k,133))
         mat(k,1998) = -rxt(k,233)*y(k,19)
         mat(k,1799) = rxt(k,232)*y(k,124)
         mat(k,1551) = rxt(k,232)*y(k,18)
         mat(k,241) = -(rxt(k,414)*y(k,221))
         mat(k,1665) = -rxt(k,414)*y(k,21)
         mat(k,2268) = rxt(k,417)*y(k,197)
         mat(k,430) = rxt(k,417)*y(k,123)
         mat(k,339) = -(rxt(k,416)*y(k,221))
         mat(k,1679) = -rxt(k,416)*y(k,22)
         mat(k,431) = rxt(k,415)*y(k,207)
         mat(k,1840) = rxt(k,415)*y(k,197)
         mat(k,294) = -(rxt(k,288)*y(k,55) + rxt(k,289)*y(k,221))
         mat(k,2188) = -rxt(k,288)*y(k,23)
         mat(k,1673) = -rxt(k,289)*y(k,23)
         mat(k,535) = -(rxt(k,290)*y(k,55) + rxt(k,291)*y(k,137) + rxt(k,316)*y(k,221))
         mat(k,2193) = -rxt(k,290)*y(k,24)
         mat(k,1940) = -rxt(k,291)*y(k,24)
         mat(k,1707) = -rxt(k,316)*y(k,24)
         mat(k,268) = -(rxt(k,296)*y(k,221))
         mat(k,1670) = -rxt(k,296)*y(k,25)
         mat(k,856) = .800_r8*rxt(k,292)*y(k,198) + .200_r8*rxt(k,293)*y(k,202)
         mat(k,2134) = .200_r8*rxt(k,293)*y(k,198)
         mat(k,347) = -(rxt(k,297)*y(k,221))
         mat(k,1681) = -rxt(k,297)*y(k,26)
         mat(k,857) = rxt(k,294)*y(k,207)
         mat(k,1841) = rxt(k,294)*y(k,198)
         mat(k,300) = -(rxt(k,298)*y(k,55) + rxt(k,299)*y(k,221))
         mat(k,2189) = -rxt(k,298)*y(k,27)
         mat(k,1674) = -rxt(k,299)*y(k,27)
         mat(k,1023) = -(rxt(k,319)*y(k,125) + rxt(k,320)*y(k,137) + rxt(k,337) &
                      *y(k,221))
         mat(k,2091) = -rxt(k,319)*y(k,28)
         mat(k,1954) = -rxt(k,320)*y(k,28)
         mat(k,1753) = -rxt(k,337)*y(k,28)
         mat(k,836) = .130_r8*rxt(k,397)*y(k,137)
         mat(k,1954) = mat(k,1954) + .130_r8*rxt(k,397)*y(k,97)
         mat(k,406) = -(rxt(k,324)*y(k,221))
         mat(k,1689) = -rxt(k,324)*y(k,29)
         mat(k,803) = rxt(k,322)*y(k,207)
         mat(k,1845) = rxt(k,322)*y(k,199)
         mat(k,110) = -(rxt(k,325)*y(k,221))
         mat(k,1648) = -rxt(k,325)*y(k,30)
         mat(k,272) = -(rxt(k,420)*y(k,221))
         mat(k,1671) = -rxt(k,420)*y(k,31)
         mat(k,618) = rxt(k,418)*y(k,207)
         mat(k,1836) = rxt(k,418)*y(k,200)
         mat(k,100) = -(rxt(k,212)*y(k,220))
         mat(k,1589) = -rxt(k,212)*y(k,32)
         mat(k,142) = -(rxt(k,213)*y(k,220))
         mat(k,1594) = -rxt(k,213)*y(k,33)
         mat(k,147) = -(rxt(k,239)*y(k,220))
         mat(k,1595) = -rxt(k,239)*y(k,34)
         mat(k,114) = -(rxt(k,214)*y(k,220))
         mat(k,1591) = -rxt(k,214)*y(k,35)
         mat(k,152) = -(rxt(k,215)*y(k,220))
         mat(k,1596) = -rxt(k,215)*y(k,36)
         mat(k,118) = -(rxt(k,216)*y(k,220))
         mat(k,1592) = -rxt(k,216)*y(k,37)
         mat(k,157) = -(rxt(k,217)*y(k,220))
         mat(k,1597) = -rxt(k,217)*y(k,38)
         mat(k,122) = -(rxt(k,218)*y(k,220))
         mat(k,1593) = -rxt(k,218)*y(k,39)
         mat(k,476) = -(rxt(k,250)*y(k,55) + rxt(k,251)*y(k,221) + rxt(k,262)*y(k,220))
         mat(k,2192) = -rxt(k,250)*y(k,40)
         mat(k,1699) = -rxt(k,251)*y(k,40)
         mat(k,1607) = -rxt(k,262)*y(k,40)
         mat(k,1528) = -(rxt(k,187)*y(k,55) + rxt(k,223)*y(k,16) + rxt(k,267)*y(k,207) &
                      + rxt(k,268)*y(k,125) + rxt(k,269)*y(k,133) + rxt(k,270) &
                      *y(k,221))
         mat(k,2212) = -rxt(k,187)*y(k,41)
         mat(k,1448) = -rxt(k,223)*y(k,41)
         mat(k,1914) = -rxt(k,267)*y(k,41)
         mat(k,2118) = -rxt(k,268)*y(k,41)
         mat(k,2013) = -rxt(k,269)*y(k,41)
         mat(k,1781) = -rxt(k,270)*y(k,41)
         mat(k,667) = .400_r8*rxt(k,368)*y(k,221)
         mat(k,949) = .340_r8*rxt(k,452)*y(k,137)
         mat(k,356) = .500_r8*rxt(k,339)*y(k,125)
         mat(k,539) = rxt(k,291)*y(k,137)
         mat(k,1032) = .500_r8*rxt(k,320)*y(k,137)
         mat(k,608) = .500_r8*rxt(k,308)*y(k,221)
         mat(k,800) = rxt(k,275)*y(k,221)
         mat(k,452) = .300_r8*rxt(k,276)*y(k,221)
         mat(k,1465) = (rxt(k,284)+rxt(k,285))*y(k,220)
         mat(k,2239) = rxt(k,194)*y(k,202)
         mat(k,1113) = .800_r8*rxt(k,313)*y(k,221)
         mat(k,844) = .910_r8*rxt(k,397)*y(k,137)
         mat(k,581) = .300_r8*rxt(k,388)*y(k,221)
         mat(k,1216) = .800_r8*rxt(k,392)*y(k,202)
         mat(k,1233) = .120_r8*rxt(k,350)*y(k,137)
         mat(k,571) = .500_r8*rxt(k,363)*y(k,221)
         mat(k,921) = .340_r8*rxt(k,455)*y(k,137)
         mat(k,1363) = .600_r8*rxt(k,364)*y(k,137)
         mat(k,2332) = .100_r8*rxt(k,370)*y(k,193) + rxt(k,274)*y(k,202) &
                      + .500_r8*rxt(k,341)*y(k,204) + .500_r8*rxt(k,310)*y(k,206) &
                      + .920_r8*rxt(k,380)*y(k,209) + .250_r8*rxt(k,348)*y(k,213) &
                      + rxt(k,357)*y(k,215) + rxt(k,331)*y(k,223) + rxt(k,335) &
                      *y(k,224) + .340_r8*rxt(k,464)*y(k,225) + .320_r8*rxt(k,469) &
                      *y(k,226) + .250_r8*rxt(k,405)*y(k,229)
         mat(k,2118) = mat(k,2118) + .500_r8*rxt(k,339)*y(k,15) + rxt(k,381)*y(k,209) &
                      + .250_r8*rxt(k,347)*y(k,213) + rxt(k,358)*y(k,215)
         mat(k,1979) = .340_r8*rxt(k,452)*y(k,5) + rxt(k,291)*y(k,24) &
                      + .500_r8*rxt(k,320)*y(k,28) + .910_r8*rxt(k,397)*y(k,97) &
                      + .120_r8*rxt(k,350)*y(k,104) + .340_r8*rxt(k,455)*y(k,109) &
                      + .600_r8*rxt(k,364)*y(k,110)
         mat(k,522) = rxt(k,315)*y(k,221)
         mat(k,1102) = .680_r8*rxt(k,473)*y(k,221)
         mat(k,965) = .100_r8*rxt(k,370)*y(k,123)
         mat(k,862) = .700_r8*rxt(k,293)*y(k,202)
         mat(k,808) = rxt(k,321)*y(k,202)
         mat(k,1417) = rxt(k,304)*y(k,202) + rxt(k,377)*y(k,209) + .250_r8*rxt(k,344) &
                      *y(k,213) + rxt(k,353)*y(k,215) + .250_r8*rxt(k,402)*y(k,229)
         mat(k,2171) = rxt(k,194)*y(k,58) + .800_r8*rxt(k,392)*y(k,100) + rxt(k,274) &
                      *y(k,123) + .700_r8*rxt(k,293)*y(k,198) + rxt(k,321)*y(k,199) &
                      + rxt(k,304)*y(k,201) + (4.000_r8*rxt(k,271)+2.000_r8*rxt(k,272)) &
                      *y(k,202) + 1.500_r8*rxt(k,378)*y(k,209) + .750_r8*rxt(k,383) &
                      *y(k,210) + .880_r8*rxt(k,345)*y(k,213) + 2.000_r8*rxt(k,354) &
                      *y(k,215) + .750_r8*rxt(k,457)*y(k,219) + .800_r8*rxt(k,333) &
                      *y(k,224) + .930_r8*rxt(k,462)*y(k,225) + .950_r8*rxt(k,467) &
                      *y(k,226) + .800_r8*rxt(k,403)*y(k,229)
         mat(k,555) = .500_r8*rxt(k,341)*y(k,123)
         mat(k,783) = .500_r8*rxt(k,310)*y(k,123)
         mat(k,1914) = mat(k,1914) + .450_r8*rxt(k,355)*y(k,215) + .150_r8*rxt(k,334) &
                      *y(k,224)
         mat(k,1282) = .920_r8*rxt(k,380)*y(k,123) + rxt(k,381)*y(k,125) + rxt(k,377) &
                      *y(k,201) + 1.500_r8*rxt(k,378)*y(k,202)
         mat(k,1317) = .750_r8*rxt(k,383)*y(k,202)
         mat(k,1340) = .250_r8*rxt(k,348)*y(k,123) + .250_r8*rxt(k,347)*y(k,125) &
                      + .250_r8*rxt(k,344)*y(k,201) + .880_r8*rxt(k,345)*y(k,202)
         mat(k,1384) = rxt(k,357)*y(k,123) + rxt(k,358)*y(k,125) + rxt(k,353)*y(k,201) &
                      + 2.000_r8*rxt(k,354)*y(k,202) + .450_r8*rxt(k,355)*y(k,207) &
                      + 4.000_r8*rxt(k,356)*y(k,215)
         mat(k,1065) = .750_r8*rxt(k,457)*y(k,202)
         mat(k,1616) = (rxt(k,284)+rxt(k,285))*y(k,53)
         mat(k,1781) = mat(k,1781) + .400_r8*rxt(k,368)*y(k,1) + .500_r8*rxt(k,308) &
                      *y(k,50) + rxt(k,275)*y(k,51) + .300_r8*rxt(k,276)*y(k,52) &
                      + .800_r8*rxt(k,313)*y(k,73) + .300_r8*rxt(k,388)*y(k,98) &
                      + .500_r8*rxt(k,363)*y(k,108) + rxt(k,315)*y(k,142) &
                      + .680_r8*rxt(k,473)*y(k,182)
         mat(k,793) = rxt(k,331)*y(k,123)
         mat(k,1167) = rxt(k,335)*y(k,123) + .800_r8*rxt(k,333)*y(k,202) &
                      + .150_r8*rxt(k,334)*y(k,207)
         mat(k,1153) = .340_r8*rxt(k,464)*y(k,123) + .930_r8*rxt(k,462)*y(k,202)
         mat(k,1131) = .320_r8*rxt(k,469)*y(k,123) + .950_r8*rxt(k,467)*y(k,202)
         mat(k,1191) = .250_r8*rxt(k,405)*y(k,123) + .250_r8*rxt(k,402)*y(k,201) &
                      + .800_r8*rxt(k,403)*y(k,202)
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
         mat(k,589) = -(rxt(k,252)*y(k,55) + rxt(k,253)*y(k,221) + rxt(k,263)*y(k,220))
         mat(k,2194) = -rxt(k,252)*y(k,42)
         mat(k,1713) = -rxt(k,253)*y(k,42)
         mat(k,1608) = -rxt(k,263)*y(k,42)
         mat(k,126) = -(rxt(k,254)*y(k,221))
         mat(k,1649) = -rxt(k,254)*y(k,43)
         mat(k,1047) = -(rxt(k,300)*y(k,125) + rxt(k,301)*y(k,221))
         mat(k,2093) = -rxt(k,300)*y(k,44)
         mat(k,1755) = -rxt(k,301)*y(k,44)
         mat(k,665) = .800_r8*rxt(k,368)*y(k,221)
         mat(k,355) = rxt(k,339)*y(k,125)
         mat(k,269) = rxt(k,296)*y(k,221)
         mat(k,349) = .500_r8*rxt(k,297)*y(k,221)
         mat(k,1025) = .500_r8*rxt(k,320)*y(k,137)
         mat(k,1352) = .100_r8*rxt(k,364)*y(k,137)
         mat(k,2308) = .400_r8*rxt(k,370)*y(k,193) + rxt(k,295)*y(k,198) &
                      + .270_r8*rxt(k,323)*y(k,199) + rxt(k,341)*y(k,204) + rxt(k,360) &
                      *y(k,217) + rxt(k,331)*y(k,223)
         mat(k,2093) = mat(k,2093) + rxt(k,339)*y(k,15)
         mat(k,1956) = .500_r8*rxt(k,320)*y(k,28) + .100_r8*rxt(k,364)*y(k,110)
         mat(k,962) = .400_r8*rxt(k,370)*y(k,123)
         mat(k,860) = rxt(k,295)*y(k,123) + 3.200_r8*rxt(k,292)*y(k,198) &
                      + .800_r8*rxt(k,293)*y(k,202)
         mat(k,806) = .270_r8*rxt(k,323)*y(k,123)
         mat(k,2149) = .800_r8*rxt(k,293)*y(k,198)
         mat(k,553) = rxt(k,341)*y(k,123)
         mat(k,1890) = .200_r8*rxt(k,359)*y(k,217)
         mat(k,673) = rxt(k,360)*y(k,123) + .200_r8*rxt(k,359)*y(k,207)
         mat(k,1755) = mat(k,1755) + .800_r8*rxt(k,368)*y(k,1) + rxt(k,296)*y(k,25) &
                      + .500_r8*rxt(k,297)*y(k,26)
         mat(k,790) = rxt(k,331)*y(k,123)
         mat(k,368) = -(rxt(k,255)*y(k,55) + rxt(k,256)*y(k,221))
         mat(k,2190) = -rxt(k,255)*y(k,45)
         mat(k,1684) = -rxt(k,256)*y(k,45)
         mat(k,103) = -(rxt(k,302)*y(k,221))
         mat(k,1647) = -rxt(k,302)*y(k,46)
         mat(k,971) = -(rxt(k,338)*y(k,221))
         mat(k,1748) = -rxt(k,338)*y(k,47)
         mat(k,664) = .800_r8*rxt(k,368)*y(k,221)
         mat(k,938) = .520_r8*rxt(k,452)*y(k,137)
         mat(k,354) = .500_r8*rxt(k,339)*y(k,125)
         mat(k,910) = .520_r8*rxt(k,455)*y(k,137)
         mat(k,2303) = .250_r8*rxt(k,370)*y(k,193) + .820_r8*rxt(k,323)*y(k,199) &
                      + .500_r8*rxt(k,341)*y(k,204) + .270_r8*rxt(k,464)*y(k,225) &
                      + .040_r8*rxt(k,469)*y(k,226)
         mat(k,2087) = .500_r8*rxt(k,339)*y(k,15)
         mat(k,1951) = .520_r8*rxt(k,452)*y(k,5) + .520_r8*rxt(k,455)*y(k,109)
         mat(k,1095) = .500_r8*rxt(k,473)*y(k,221)
         mat(k,961) = .250_r8*rxt(k,370)*y(k,123)
         mat(k,805) = .820_r8*rxt(k,323)*y(k,123) + .820_r8*rxt(k,321)*y(k,202)
         mat(k,2144) = .820_r8*rxt(k,321)*y(k,199) + .150_r8*rxt(k,462)*y(k,225) &
                      + .025_r8*rxt(k,467)*y(k,226)
         mat(k,552) = .500_r8*rxt(k,341)*y(k,123)
         mat(k,1748) = mat(k,1748) + .800_r8*rxt(k,368)*y(k,1) + .500_r8*rxt(k,473) &
                      *y(k,182)
         mat(k,1143) = .270_r8*rxt(k,464)*y(k,123) + .150_r8*rxt(k,462)*y(k,202)
         mat(k,1120) = .040_r8*rxt(k,469)*y(k,123) + .025_r8*rxt(k,467)*y(k,202)
         mat(k,1240) = -(rxt(k,326)*y(k,125) + rxt(k,327)*y(k,221))
         mat(k,2106) = -rxt(k,326)*y(k,48)
         mat(k,1768) = -rxt(k,327)*y(k,48)
         mat(k,1090) = rxt(k,328)*y(k,221)
         mat(k,1228) = .880_r8*rxt(k,350)*y(k,137)
         mat(k,1355) = .500_r8*rxt(k,364)*y(k,137)
         mat(k,2321) = .170_r8*rxt(k,423)*y(k,203) + .050_r8*rxt(k,386)*y(k,210) &
                      + .250_r8*rxt(k,348)*y(k,213) + .170_r8*rxt(k,429)*y(k,216) &
                      + .400_r8*rxt(k,439)*y(k,227) + .250_r8*rxt(k,405)*y(k,229) &
                      + .540_r8*rxt(k,445)*y(k,230) + .510_r8*rxt(k,448)*y(k,232)
         mat(k,2106) = mat(k,2106) + .050_r8*rxt(k,387)*y(k,210) + .250_r8*rxt(k,347) &
                      *y(k,213) + .250_r8*rxt(k,406)*y(k,229)
         mat(k,851) = rxt(k,329)*y(k,221)
         mat(k,1967) = .880_r8*rxt(k,350)*y(k,104) + .500_r8*rxt(k,364)*y(k,110)
         mat(k,1407) = .250_r8*rxt(k,344)*y(k,213) + .250_r8*rxt(k,402)*y(k,229)
         mat(k,2161) = .240_r8*rxt(k,345)*y(k,213) + .500_r8*rxt(k,333)*y(k,224) &
                      + .100_r8*rxt(k,403)*y(k,229)
         mat(k,763) = .170_r8*rxt(k,423)*y(k,123) + .070_r8*rxt(k,422)*y(k,207)
         mat(k,1902) = .070_r8*rxt(k,422)*y(k,203) + .070_r8*rxt(k,428)*y(k,216)
         mat(k,1309) = .050_r8*rxt(k,386)*y(k,123) + .050_r8*rxt(k,387)*y(k,125)
         mat(k,1334) = .250_r8*rxt(k,348)*y(k,123) + .250_r8*rxt(k,347)*y(k,125) &
                      + .250_r8*rxt(k,344)*y(k,201) + .240_r8*rxt(k,345)*y(k,202)
         mat(k,876) = .170_r8*rxt(k,429)*y(k,123) + .070_r8*rxt(k,428)*y(k,207)
         mat(k,1768) = mat(k,1768) + rxt(k,328)*y(k,94) + rxt(k,329)*y(k,126)
         mat(k,1165) = .500_r8*rxt(k,333)*y(k,202)
         mat(k,739) = .400_r8*rxt(k,439)*y(k,123)
         mat(k,1188) = .250_r8*rxt(k,405)*y(k,123) + .250_r8*rxt(k,406)*y(k,125) &
                      + .250_r8*rxt(k,402)*y(k,201) + .100_r8*rxt(k,403)*y(k,202)
         mat(k,755) = .540_r8*rxt(k,445)*y(k,123)
         mat(k,502) = .510_r8*rxt(k,448)*y(k,123)
         mat(k,688) = -(rxt(k,307)*y(k,221))
         mat(k,1724) = -rxt(k,307)*y(k,49)
         mat(k,1019) = .120_r8*rxt(k,320)*y(k,137)
         mat(k,1942) = .120_r8*rxt(k,320)*y(k,28)
         mat(k,1397) = .100_r8*rxt(k,304)*y(k,202) + .150_r8*rxt(k,305)*y(k,207)
         mat(k,2139) = .100_r8*rxt(k,304)*y(k,201)
         mat(k,1868) = .150_r8*rxt(k,305)*y(k,201) + .150_r8*rxt(k,355)*y(k,215)
         mat(k,1375) = .150_r8*rxt(k,355)*y(k,207)
         mat(k,605) = -(rxt(k,308)*y(k,221))
         mat(k,1715) = -rxt(k,308)*y(k,50)
         mat(k,1396) = .360_r8*rxt(k,305)*y(k,207)
         mat(k,1862) = .360_r8*rxt(k,305)*y(k,201) + .400_r8*rxt(k,355)*y(k,215)
         mat(k,1374) = .400_r8*rxt(k,355)*y(k,207)
         mat(k,799) = -(rxt(k,275)*y(k,221))
         mat(k,1734) = -rxt(k,275)*y(k,51)
         mat(k,1202) = .200_r8*rxt(k,392)*y(k,202)
         mat(k,858) = .300_r8*rxt(k,293)*y(k,202)
         mat(k,2140) = .200_r8*rxt(k,392)*y(k,100) + .300_r8*rxt(k,293)*y(k,198) &
                      + 2.000_r8*rxt(k,272)*y(k,202) + .250_r8*rxt(k,378)*y(k,209) &
                      + .250_r8*rxt(k,383)*y(k,210) + .250_r8*rxt(k,345)*y(k,213) &
                      + .250_r8*rxt(k,457)*y(k,219) + .500_r8*rxt(k,333)*y(k,224) &
                      + .250_r8*rxt(k,462)*y(k,225) + .250_r8*rxt(k,467)*y(k,226) &
                      + .300_r8*rxt(k,403)*y(k,229)
         mat(k,1265) = .250_r8*rxt(k,378)*y(k,202)
         mat(k,1297) = .250_r8*rxt(k,383)*y(k,202)
         mat(k,1328) = .250_r8*rxt(k,345)*y(k,202)
         mat(k,1058) = .250_r8*rxt(k,457)*y(k,202)
         mat(k,1162) = .500_r8*rxt(k,333)*y(k,202)
         mat(k,1142) = .250_r8*rxt(k,462)*y(k,202)
         mat(k,1119) = .250_r8*rxt(k,467)*y(k,202)
         mat(k,1181) = .300_r8*rxt(k,403)*y(k,202)
         mat(k,451) = -(rxt(k,276)*y(k,221))
         mat(k,1695) = -rxt(k,276)*y(k,52)
         mat(k,2137) = rxt(k,273)*y(k,207)
         mat(k,1851) = rxt(k,273)*y(k,202)
         mat(k,1462) = -(rxt(k,188)*y(k,55) + rxt(k,244)*y(k,72) + rxt(k,277)*y(k,221) &
                      + (rxt(k,283) + rxt(k,284) + rxt(k,285)) * y(k,220))
         mat(k,2209) = -rxt(k,188)*y(k,53)
         mat(k,885) = -rxt(k,244)*y(k,53)
         mat(k,1778) = -rxt(k,277)*y(k,53)
         mat(k,1613) = -(rxt(k,283) + rxt(k,284) + rxt(k,285)) * y(k,53)
         mat(k,1030) = .100_r8*rxt(k,320)*y(k,137)
         mat(k,1976) = .100_r8*rxt(k,320)*y(k,28)
         mat(k,382) = -(rxt(k,240)*y(k,220) + rxt(k,257)*y(k,55) + rxt(k,258)*y(k,221))
         mat(k,1606) = -rxt(k,240)*y(k,54)
         mat(k,2191) = -rxt(k,257)*y(k,54)
         mat(k,1685) = -rxt(k,258)*y(k,54)
         mat(k,2224) = -(rxt(k,187)*y(k,41) + rxt(k,188)*y(k,53) + rxt(k,189)*y(k,76) &
                      + rxt(k,190)*y(k,78) + (rxt(k,191) + rxt(k,192)) * y(k,207) &
                      + rxt(k,193)*y(k,137) + rxt(k,200)*y(k,59) + rxt(k,209)*y(k,91) &
                      + rxt(k,250)*y(k,40) + rxt(k,252)*y(k,42) + rxt(k,255)*y(k,45) &
                      + rxt(k,257)*y(k,54) + rxt(k,298)*y(k,27))
         mat(k,1540) = -rxt(k,187)*y(k,55)
         mat(k,1474) = -rxt(k,188)*y(k,55)
         mat(k,1261) = -rxt(k,189)*y(k,55)
         mat(k,603) = -rxt(k,190)*y(k,55)
         mat(k,1926) = -(rxt(k,191) + rxt(k,192)) * y(k,55)
         mat(k,1991) = -rxt(k,193)*y(k,55)
         mat(k,992) = -rxt(k,200)*y(k,55)
         mat(k,828) = -rxt(k,209)*y(k,55)
         mat(k,481) = -rxt(k,250)*y(k,55)
         mat(k,596) = -rxt(k,252)*y(k,55)
         mat(k,374) = -rxt(k,255)*y(k,55)
         mat(k,387) = -rxt(k,257)*y(k,55)
         mat(k,304) = -rxt(k,298)*y(k,55)
         mat(k,1818) = rxt(k,228)*y(k,58)
         mat(k,102) = 4.000_r8*rxt(k,212)*y(k,220)
         mat(k,146) = rxt(k,213)*y(k,220)
         mat(k,117) = 2.000_r8*rxt(k,214)*y(k,220)
         mat(k,156) = 2.000_r8*rxt(k,215)*y(k,220)
         mat(k,121) = 2.000_r8*rxt(k,216)*y(k,220)
         mat(k,161) = rxt(k,217)*y(k,220)
         mat(k,125) = 2.000_r8*rxt(k,218)*y(k,220)
         mat(k,128) = 3.000_r8*rxt(k,254)*y(k,221)
         mat(k,374) = mat(k,374) + rxt(k,256)*y(k,221)
         mat(k,2251) = rxt(k,228)*y(k,18) + (4.000_r8*rxt(k,195)+2.000_r8*rxt(k,197)) &
                      *y(k,58) + rxt(k,199)*y(k,123) + rxt(k,204)*y(k,133) &
                      + rxt(k,482)*y(k,153) + rxt(k,194)*y(k,202) + rxt(k,205) &
                      *y(k,221)
         mat(k,252) = rxt(k,249)*y(k,220)
         mat(k,248) = rxt(k,264)*y(k,220) + rxt(k,259)*y(k,221)
         mat(k,258) = rxt(k,265)*y(k,220) + rxt(k,260)*y(k,221)
         mat(k,311) = rxt(k,266)*y(k,220) + rxt(k,261)*y(k,221)
         mat(k,2069) = rxt(k,207)*y(k,133) + rxt(k,219)*y(k,220) + rxt(k,208)*y(k,221)
         mat(k,2344) = rxt(k,199)*y(k,58)
         mat(k,2025) = rxt(k,204)*y(k,58) + rxt(k,207)*y(k,84)
         mat(k,1440) = rxt(k,482)*y(k,58)
         mat(k,2183) = rxt(k,194)*y(k,58)
         mat(k,1628) = 4.000_r8*rxt(k,212)*y(k,32) + rxt(k,213)*y(k,33) &
                      + 2.000_r8*rxt(k,214)*y(k,35) + 2.000_r8*rxt(k,215)*y(k,36) &
                      + 2.000_r8*rxt(k,216)*y(k,37) + rxt(k,217)*y(k,38) &
                      + 2.000_r8*rxt(k,218)*y(k,39) + rxt(k,249)*y(k,64) + rxt(k,264) &
                      *y(k,81) + rxt(k,265)*y(k,82) + rxt(k,266)*y(k,83) + rxt(k,219) &
                      *y(k,84)
         mat(k,1793) = 3.000_r8*rxt(k,254)*y(k,43) + rxt(k,256)*y(k,45) + rxt(k,205) &
                      *y(k,58) + rxt(k,259)*y(k,81) + rxt(k,260)*y(k,82) + rxt(k,261) &
                      *y(k,83) + rxt(k,208)*y(k,84)
         mat(k,2187) = rxt(k,200)*y(k,59)
         mat(k,2229) = 2.000_r8*rxt(k,196)*y(k,58)
         mat(k,983) = rxt(k,200)*y(k,55) + (rxt(k,540)+rxt(k,545)+rxt(k,550))*y(k,84)
         mat(k,2049) = (rxt(k,540)+rxt(k,545)+rxt(k,550))*y(k,59) + (rxt(k,535) &
                       +rxt(k,541)+rxt(k,546))*y(k,91)
         mat(k,823) = (rxt(k,535)+rxt(k,541)+rxt(k,546))*y(k,84)
         mat(k,2228) = 2.000_r8*rxt(k,221)*y(k,58)
         mat(k,2252) = -(rxt(k,194)*y(k,202) + (4._r8*rxt(k,195) + 4._r8*rxt(k,196) &
                      + 4._r8*rxt(k,197) + 4._r8*rxt(k,221)) * y(k,58) + rxt(k,198) &
                      *y(k,207) + rxt(k,199)*y(k,123) + rxt(k,201)*y(k,124) + rxt(k,204) &
                      *y(k,133) + (rxt(k,205) + rxt(k,206)) * y(k,221) + (rxt(k,227) &
                      + rxt(k,228) + rxt(k,229)) * y(k,18) + rxt(k,482)*y(k,153))
         mat(k,2184) = -rxt(k,194)*y(k,58)
         mat(k,1927) = -rxt(k,198)*y(k,58)
         mat(k,2345) = -rxt(k,199)*y(k,58)
         mat(k,1586) = -rxt(k,201)*y(k,58)
         mat(k,2026) = -rxt(k,204)*y(k,58)
         mat(k,1794) = -(rxt(k,205) + rxt(k,206)) * y(k,58)
         mat(k,1819) = -(rxt(k,227) + rxt(k,228) + rxt(k,229)) * y(k,58)
         mat(k,1441) = -rxt(k,482)*y(k,58)
         mat(k,2225) = rxt(k,209)*y(k,91) + rxt(k,193)*y(k,137) + rxt(k,192)*y(k,207)
         mat(k,993) = rxt(k,202)*y(k,133)
         mat(k,2070) = rxt(k,220)*y(k,220)
         mat(k,829) = rxt(k,209)*y(k,55) + rxt(k,210)*y(k,133) + rxt(k,211)*y(k,221)
         mat(k,2026) = mat(k,2026) + rxt(k,202)*y(k,59) + rxt(k,210)*y(k,91)
         mat(k,1992) = rxt(k,193)*y(k,55)
         mat(k,332) = rxt(k,487)*y(k,153)
         mat(k,1441) = mat(k,1441) + rxt(k,487)*y(k,139)
         mat(k,1927) = mat(k,1927) + rxt(k,192)*y(k,55)
         mat(k,1629) = rxt(k,220)*y(k,84)
         mat(k,1794) = mat(k,1794) + rxt(k,211)*y(k,91)
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
         mat(k,985) = -(rxt(k,200)*y(k,55) + rxt(k,202)*y(k,133) + rxt(k,203)*y(k,221) &
                      + (rxt(k,540) + rxt(k,545) + rxt(k,550)) * y(k,84))
         mat(k,2201) = -rxt(k,200)*y(k,59)
         mat(k,2005) = -rxt(k,202)*y(k,59)
         mat(k,1750) = -rxt(k,203)*y(k,59)
         mat(k,2053) = -(rxt(k,540) + rxt(k,545) + rxt(k,550)) * y(k,59)
         mat(k,2234) = rxt(k,201)*y(k,124)
         mat(k,1561) = rxt(k,201)*y(k,58)
         mat(k,1175) = -(rxt(k,287)*y(k,221))
         mat(k,1764) = -rxt(k,287)*y(k,61)
         mat(k,946) = .230_r8*rxt(k,452)*y(k,137)
         mat(k,1445) = rxt(k,223)*y(k,41)
         mat(k,297) = .350_r8*rxt(k,289)*y(k,221)
         mat(k,538) = .630_r8*rxt(k,291)*y(k,137)
         mat(k,1027) = .560_r8*rxt(k,320)*y(k,137)
         mat(k,1522) = rxt(k,223)*y(k,16) + rxt(k,187)*y(k,55) + rxt(k,268)*y(k,125) &
                      + rxt(k,269)*y(k,133) + rxt(k,270)*y(k,221)
         mat(k,369) = rxt(k,255)*y(k,55)
         mat(k,1239) = rxt(k,326)*y(k,125) + rxt(k,327)*y(k,221)
         mat(k,2205) = rxt(k,187)*y(k,41) + rxt(k,255)*y(k,45)
         mat(k,980) = rxt(k,314)*y(k,221)
         mat(k,837) = .620_r8*rxt(k,397)*y(k,137)
         mat(k,1226) = .650_r8*rxt(k,350)*y(k,137)
         mat(k,918) = .230_r8*rxt(k,455)*y(k,137)
         mat(k,1354) = .560_r8*rxt(k,364)*y(k,137)
         mat(k,2317) = .170_r8*rxt(k,423)*y(k,203) + .220_r8*rxt(k,348)*y(k,213) &
                      + .400_r8*rxt(k,426)*y(k,214) + .350_r8*rxt(k,429)*y(k,216) &
                      + .225_r8*rxt(k,464)*y(k,225) + .250_r8*rxt(k,405)*y(k,229)
         mat(k,2102) = rxt(k,268)*y(k,41) + rxt(k,326)*y(k,48) + .220_r8*rxt(k,347) &
                      *y(k,213) + .500_r8*rxt(k,406)*y(k,229)
         mat(k,2006) = rxt(k,269)*y(k,41) + rxt(k,476)*y(k,140)
         mat(k,1964) = .230_r8*rxt(k,452)*y(k,5) + .630_r8*rxt(k,291)*y(k,24) &
                      + .560_r8*rxt(k,320)*y(k,28) + .620_r8*rxt(k,397)*y(k,97) &
                      + .650_r8*rxt(k,350)*y(k,104) + .230_r8*rxt(k,455)*y(k,109) &
                      + .560_r8*rxt(k,364)*y(k,110)
         mat(k,363) = rxt(k,476)*y(k,133) + rxt(k,477)*y(k,221)
         mat(k,1099) = .700_r8*rxt(k,473)*y(k,221)
         mat(k,1403) = .220_r8*rxt(k,344)*y(k,213) + .250_r8*rxt(k,402)*y(k,229)
         mat(k,2157) = .110_r8*rxt(k,345)*y(k,213) + .125_r8*rxt(k,462)*y(k,225) &
                      + .200_r8*rxt(k,403)*y(k,229)
         mat(k,762) = .170_r8*rxt(k,423)*y(k,123) + .070_r8*rxt(k,422)*y(k,207)
         mat(k,1898) = .070_r8*rxt(k,422)*y(k,203) + .160_r8*rxt(k,425)*y(k,214) &
                      + .140_r8*rxt(k,428)*y(k,216)
         mat(k,1333) = .220_r8*rxt(k,348)*y(k,123) + .220_r8*rxt(k,347)*y(k,125) &
                      + .220_r8*rxt(k,344)*y(k,201) + .110_r8*rxt(k,345)*y(k,202)
         mat(k,725) = .400_r8*rxt(k,426)*y(k,123) + .160_r8*rxt(k,425)*y(k,207)
         mat(k,875) = .350_r8*rxt(k,429)*y(k,123) + .140_r8*rxt(k,428)*y(k,207)
         mat(k,1764) = mat(k,1764) + .350_r8*rxt(k,289)*y(k,23) + rxt(k,270)*y(k,41) &
                      + rxt(k,327)*y(k,48) + rxt(k,314)*y(k,74) + rxt(k,477)*y(k,140) &
                      + .700_r8*rxt(k,473)*y(k,182)
         mat(k,1149) = .225_r8*rxt(k,464)*y(k,123) + .125_r8*rxt(k,462)*y(k,202)
         mat(k,1186) = .250_r8*rxt(k,405)*y(k,123) + .500_r8*rxt(k,406)*y(k,125) &
                      + .250_r8*rxt(k,402)*y(k,201) + .200_r8*rxt(k,403)*y(k,202)
         mat(k,941) = .270_r8*rxt(k,452)*y(k,137)
         mat(k,1024) = .200_r8*rxt(k,320)*y(k,137)
         mat(k,689) = rxt(k,307)*y(k,221)
         mat(k,606) = .500_r8*rxt(k,308)*y(k,221)
         mat(k,1174) = rxt(k,287)*y(k,221)
         mat(k,1109) = .800_r8*rxt(k,313)*y(k,221)
         mat(k,979) = rxt(k,314)*y(k,221)
         mat(k,869) = rxt(k,279)*y(k,221)
         mat(k,568) = .500_r8*rxt(k,363)*y(k,221)
         mat(k,913) = .270_r8*rxt(k,455)*y(k,137)
         mat(k,1351) = .100_r8*rxt(k,364)*y(k,137)
         mat(k,2307) = rxt(k,306)*y(k,201) + .900_r8*rxt(k,464)*y(k,225)
         mat(k,1955) = .270_r8*rxt(k,452)*y(k,5) + .200_r8*rxt(k,320)*y(k,28) &
                      + .270_r8*rxt(k,455)*y(k,109) + .100_r8*rxt(k,364)*y(k,110)
         mat(k,1096) = 1.800_r8*rxt(k,473)*y(k,221)
         mat(k,1400) = rxt(k,306)*y(k,123) + 4.000_r8*rxt(k,303)*y(k,201) &
                      + .900_r8*rxt(k,304)*y(k,202) + .490_r8*rxt(k,305)*y(k,207) &
                      + rxt(k,377)*y(k,209) + 2.000_r8*rxt(k,353)*y(k,215) &
                      + rxt(k,402)*y(k,229)
         mat(k,2148) = .900_r8*rxt(k,304)*y(k,201) + rxt(k,354)*y(k,215) &
                      + .500_r8*rxt(k,462)*y(k,225)
         mat(k,1889) = .490_r8*rxt(k,305)*y(k,201) + .450_r8*rxt(k,355)*y(k,215)
         mat(k,1266) = rxt(k,377)*y(k,201)
         mat(k,1376) = 2.000_r8*rxt(k,353)*y(k,201) + rxt(k,354)*y(k,202) &
                      + .450_r8*rxt(k,355)*y(k,207) + 4.000_r8*rxt(k,356)*y(k,215)
         mat(k,1754) = rxt(k,307)*y(k,49) + .500_r8*rxt(k,308)*y(k,50) + rxt(k,287) &
                      *y(k,61) + .800_r8*rxt(k,313)*y(k,73) + rxt(k,314)*y(k,74) &
                      + rxt(k,279)*y(k,86) + .500_r8*rxt(k,363)*y(k,108) &
                      + 1.800_r8*rxt(k,473)*y(k,182)
         mat(k,1144) = .900_r8*rxt(k,464)*y(k,123) + .500_r8*rxt(k,462)*y(k,202)
         mat(k,1183) = rxt(k,402)*y(k,201)
         mat(k,235) = -(rxt(k,248)*y(k,220))
         mat(k,1600) = -rxt(k,248)*y(k,63)
         mat(k,143) = rxt(k,213)*y(k,220)
         mat(k,148) = rxt(k,239)*y(k,220)
         mat(k,153) = rxt(k,215)*y(k,220)
         mat(k,119) = 2.000_r8*rxt(k,216)*y(k,220)
         mat(k,158) = 2.000_r8*rxt(k,217)*y(k,220)
         mat(k,123) = rxt(k,218)*y(k,220)
         mat(k,107) = 2.000_r8*rxt(k,241)*y(k,220)
         mat(k,253) = rxt(k,265)*y(k,220) + rxt(k,260)*y(k,221)
         mat(k,306) = rxt(k,266)*y(k,220) + rxt(k,261)*y(k,221)
         mat(k,1600) = mat(k,1600) + rxt(k,213)*y(k,33) + rxt(k,239)*y(k,34) &
                      + rxt(k,215)*y(k,36) + 2.000_r8*rxt(k,216)*y(k,37) &
                      + 2.000_r8*rxt(k,217)*y(k,38) + rxt(k,218)*y(k,39) &
                      + 2.000_r8*rxt(k,241)*y(k,77) + rxt(k,265)*y(k,82) + rxt(k,266) &
                      *y(k,83)
         mat(k,1663) = rxt(k,260)*y(k,82) + rxt(k,261)*y(k,83)
         mat(k,249) = -(rxt(k,249)*y(k,220))
         mat(k,1602) = -rxt(k,249)*y(k,64)
         mat(k,115) = rxt(k,214)*y(k,220)
         mat(k,154) = rxt(k,215)*y(k,220)
         mat(k,245) = rxt(k,264)*y(k,220) + rxt(k,259)*y(k,221)
         mat(k,1602) = mat(k,1602) + rxt(k,214)*y(k,35) + rxt(k,215)*y(k,36) &
                      + rxt(k,264)*y(k,81)
         mat(k,1667) = rxt(k,259)*y(k,81)
         mat(k,197) = -(rxt(k,421)*y(k,221))
         mat(k,1656) = -rxt(k,421)*y(k,65)
         mat(k,191) = .180_r8*rxt(k,441)*y(k,221)
         mat(k,1656) = mat(k,1656) + .180_r8*rxt(k,441)*y(k,184)
         mat(k,288) = -(rxt(k,474)*y(k,125) + (rxt(k,475) + rxt(k,489)) * y(k,221))
         mat(k,2076) = -rxt(k,474)*y(k,66)
         mat(k,1672) = -(rxt(k,475) + rxt(k,489)) * y(k,66)
         mat(k,778) = rxt(k,309)*y(k,207)
         mat(k,1834) = rxt(k,309)*y(k,206)
         mat(k,883) = -(rxt(k,244)*y(k,53) + rxt(k,245)*y(k,76) + rxt(k,246)*y(k,233) &
                      + rxt(k,247)*y(k,88))
         mat(k,1458) = -rxt(k,244)*y(k,72)
         mat(k,1250) = -rxt(k,245)*y(k,72)
         mat(k,2351) = -rxt(k,246)*y(k,72)
         mat(k,1477) = -rxt(k,247)*y(k,72)
         mat(k,149) = rxt(k,239)*y(k,220)
         mat(k,159) = rxt(k,217)*y(k,220)
         mat(k,236) = 2.000_r8*rxt(k,248)*y(k,220)
         mat(k,250) = rxt(k,249)*y(k,220)
         mat(k,1610) = rxt(k,239)*y(k,34) + rxt(k,217)*y(k,38) + 2.000_r8*rxt(k,248) &
                      *y(k,63) + rxt(k,249)*y(k,64)
         mat(k,1110) = -(rxt(k,313)*y(k,221))
         mat(k,1760) = -rxt(k,313)*y(k,73)
         mat(k,579) = .700_r8*rxt(k,388)*y(k,221)
         mat(k,563) = .500_r8*rxt(k,389)*y(k,221)
         mat(k,426) = rxt(k,400)*y(k,221)
         mat(k,2313) = .050_r8*rxt(k,386)*y(k,210) + .530_r8*rxt(k,348)*y(k,213) &
                      + .225_r8*rxt(k,464)*y(k,225) + .250_r8*rxt(k,405)*y(k,229)
         mat(k,2098) = .050_r8*rxt(k,387)*y(k,210) + .530_r8*rxt(k,347)*y(k,213) &
                      + .250_r8*rxt(k,406)*y(k,229)
         mat(k,1500) = rxt(k,312)*y(k,205)
         mat(k,1402) = .530_r8*rxt(k,344)*y(k,213) + .250_r8*rxt(k,402)*y(k,229)
         mat(k,2153) = .260_r8*rxt(k,345)*y(k,213) + .125_r8*rxt(k,462)*y(k,225) &
                      + .100_r8*rxt(k,403)*y(k,229)
         mat(k,458) = rxt(k,312)*y(k,134)
         mat(k,1304) = .050_r8*rxt(k,386)*y(k,123) + .050_r8*rxt(k,387)*y(k,125)
         mat(k,1331) = .530_r8*rxt(k,348)*y(k,123) + .530_r8*rxt(k,347)*y(k,125) &
                      + .530_r8*rxt(k,344)*y(k,201) + .260_r8*rxt(k,345)*y(k,202)
         mat(k,1760) = mat(k,1760) + .700_r8*rxt(k,388)*y(k,98) + .500_r8*rxt(k,389) &
                      *y(k,99) + rxt(k,400)*y(k,114)
         mat(k,1146) = .225_r8*rxt(k,464)*y(k,123) + .125_r8*rxt(k,462)*y(k,202)
         mat(k,1185) = .250_r8*rxt(k,405)*y(k,123) + .250_r8*rxt(k,406)*y(k,125) &
                      + .250_r8*rxt(k,402)*y(k,201) + .100_r8*rxt(k,403)*y(k,202)
         mat(k,978) = -(rxt(k,314)*y(k,221))
         mat(k,1749) = -rxt(k,314)*y(k,74)
         mat(k,296) = .650_r8*rxt(k,289)*y(k,221)
         mat(k,1108) = .200_r8*rxt(k,313)*y(k,221)
         mat(k,1006) = rxt(k,401)*y(k,221)
         mat(k,2304) = rxt(k,412)*y(k,195) + .050_r8*rxt(k,386)*y(k,210) &
                      + .400_r8*rxt(k,426)*y(k,214) + .170_r8*rxt(k,429)*y(k,216) &
                      + .700_r8*rxt(k,432)*y(k,222) + .600_r8*rxt(k,439)*y(k,227) &
                      + .250_r8*rxt(k,405)*y(k,229) + .340_r8*rxt(k,445)*y(k,230) &
                      + .170_r8*rxt(k,448)*y(k,232)
         mat(k,2088) = .050_r8*rxt(k,387)*y(k,210) + .250_r8*rxt(k,406)*y(k,229)
         mat(k,494) = rxt(k,412)*y(k,123)
         mat(k,1398) = .250_r8*rxt(k,402)*y(k,229)
         mat(k,2145) = .100_r8*rxt(k,403)*y(k,229)
         mat(k,1887) = .160_r8*rxt(k,425)*y(k,214) + .070_r8*rxt(k,428)*y(k,216)
         mat(k,1300) = .050_r8*rxt(k,386)*y(k,123) + .050_r8*rxt(k,387)*y(k,125)
         mat(k,724) = .400_r8*rxt(k,426)*y(k,123) + .160_r8*rxt(k,425)*y(k,207)
         mat(k,874) = .170_r8*rxt(k,429)*y(k,123) + .070_r8*rxt(k,428)*y(k,207)
         mat(k,1749) = mat(k,1749) + .650_r8*rxt(k,289)*y(k,23) + .200_r8*rxt(k,313) &
                      *y(k,73) + rxt(k,401)*y(k,115)
         mat(k,446) = .700_r8*rxt(k,432)*y(k,123)
         mat(k,737) = .600_r8*rxt(k,439)*y(k,123)
         mat(k,1182) = .250_r8*rxt(k,405)*y(k,123) + .250_r8*rxt(k,406)*y(k,125) &
                      + .250_r8*rxt(k,402)*y(k,201) + .100_r8*rxt(k,403)*y(k,202)
         mat(k,753) = .340_r8*rxt(k,445)*y(k,123)
         mat(k,501) = .170_r8*rxt(k,448)*y(k,123)
         mat(k,2041) = -((rxt(k,147) + rxt(k,148) + rxt(k,149)) * y(k,207) + rxt(k,150) &
                      *y(k,134) + rxt(k,153)*y(k,137))
         mat(k,1922) = -(rxt(k,147) + rxt(k,148) + rxt(k,149)) * y(k,75)
         mat(k,1513) = -rxt(k,150)*y(k,75)
         mat(k,1987) = -rxt(k,153)*y(k,75)
         mat(k,1536) = rxt(k,270)*y(k,221)
         mat(k,1470) = rxt(k,284)*y(k,220)
         mat(k,2220) = rxt(k,189)*y(k,76)
         mat(k,888) = rxt(k,245)*y(k,76)
         mat(k,1257) = rxt(k,189)*y(k,55) + rxt(k,245)*y(k,72) + rxt(k,145)*y(k,133) &
                      + rxt(k,127)*y(k,220) + rxt(k,154)*y(k,221)
         mat(k,821) = rxt(k,243)*y(k,220)
         mat(k,2065) = rxt(k,220)*y(k,220)
         mat(k,686) = rxt(k,175)*y(k,221)
         mat(k,2021) = rxt(k,145)*y(k,76) + rxt(k,157)*y(k,221)
         mat(k,367) = rxt(k,477)*y(k,221)
         mat(k,722) = rxt(k,483)*y(k,221)
         mat(k,1439) = rxt(k,488)*y(k,221)
         mat(k,1624) = rxt(k,284)*y(k,53) + rxt(k,127)*y(k,76) + rxt(k,243)*y(k,80) &
                      + rxt(k,220)*y(k,84)
         mat(k,1789) = rxt(k,270)*y(k,41) + rxt(k,154)*y(k,76) + rxt(k,175)*y(k,111) &
                      + rxt(k,157)*y(k,133) + rxt(k,477)*y(k,140) + rxt(k,483) &
                      *y(k,151) + rxt(k,488)*y(k,153)
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
         mat(k,1251) = -(rxt(k,127)*y(k,220) + rxt(k,145)*y(k,133) + rxt(k,154) &
                      *y(k,221) + rxt(k,189)*y(k,55) + rxt(k,245)*y(k,72))
         mat(k,1611) = -rxt(k,127)*y(k,76)
         mat(k,2007) = -rxt(k,145)*y(k,76)
         mat(k,1769) = -rxt(k,154)*y(k,76)
         mat(k,2206) = -rxt(k,189)*y(k,76)
         mat(k,884) = -rxt(k,245)*y(k,76)
         mat(k,1461) = rxt(k,285)*y(k,220)
         mat(k,2029) = rxt(k,147)*y(k,207)
         mat(k,1903) = rxt(k,147)*y(k,75)
         mat(k,1611) = mat(k,1611) + rxt(k,285)*y(k,53)
         mat(k,106) = -(rxt(k,241)*y(k,220))
         mat(k,1590) = -rxt(k,241)*y(k,77)
         mat(k,598) = -(rxt(k,146)*y(k,133) + rxt(k,155)*y(k,221) + rxt(k,190)*y(k,55))
         mat(k,1999) = -rxt(k,146)*y(k,78)
         mat(k,1714) = -rxt(k,155)*y(k,78)
         mat(k,2195) = -rxt(k,190)*y(k,78)
         mat(k,1861) = 2.000_r8*rxt(k,161)*y(k,207)
         mat(k,1714) = mat(k,1714) + 2.000_r8*rxt(k,160)*y(k,221)
         mat(k,263) = rxt(k,490)*y(k,233)
         mat(k,2348) = rxt(k,490)*y(k,155)
         mat(k,815) = -(rxt(k,236)*y(k,133) + rxt(k,237)*y(k,221) + (rxt(k,242) &
                      + rxt(k,243)) * y(k,220))
         mat(k,2002) = -rxt(k,236)*y(k,80)
         mat(k,1736) = -rxt(k,237)*y(k,80)
         mat(k,1609) = -(rxt(k,242) + rxt(k,243)) * y(k,80)
         mat(k,1444) = rxt(k,223)*y(k,41) + rxt(k,224)*y(k,207)
         mat(k,1519) = rxt(k,223)*y(k,16)
         mat(k,1879) = rxt(k,224)*y(k,16)
         mat(k,244) = -(rxt(k,259)*y(k,221) + rxt(k,264)*y(k,220))
         mat(k,1666) = -rxt(k,259)*y(k,81)
         mat(k,1601) = -rxt(k,264)*y(k,81)
         mat(k,254) = -(rxt(k,260)*y(k,221) + rxt(k,265)*y(k,220))
         mat(k,1668) = -rxt(k,260)*y(k,82)
         mat(k,1603) = -rxt(k,265)*y(k,82)
         mat(k,307) = -(rxt(k,261)*y(k,221) + rxt(k,266)*y(k,220))
         mat(k,1675) = -rxt(k,261)*y(k,83)
         mat(k,1605) = -rxt(k,266)*y(k,83)
         mat(k,2066) = -(rxt(k,207)*y(k,133) + rxt(k,208)*y(k,221) + (rxt(k,219) &
                      + rxt(k,220)) * y(k,220) + (rxt(k,535) + rxt(k,541) + rxt(k,546) &
                      ) * y(k,91) + (rxt(k,540) + rxt(k,545) + rxt(k,550)) * y(k,59) &
                      + (rxt(k,542) + rxt(k,547)) * y(k,90))
         mat(k,2022) = -rxt(k,207)*y(k,84)
         mat(k,1790) = -rxt(k,208)*y(k,84)
         mat(k,1625) = -(rxt(k,219) + rxt(k,220)) * y(k,84)
         mat(k,827) = -(rxt(k,535) + rxt(k,541) + rxt(k,546)) * y(k,84)
         mat(k,990) = -(rxt(k,540) + rxt(k,545) + rxt(k,550)) * y(k,84)
         mat(k,775) = -(rxt(k,542) + rxt(k,547)) * y(k,84)
         mat(k,303) = rxt(k,298)*y(k,55)
         mat(k,480) = rxt(k,250)*y(k,55)
         mat(k,1537) = rxt(k,187)*y(k,55)
         mat(k,594) = rxt(k,252)*y(k,55)
         mat(k,372) = 2.000_r8*rxt(k,255)*y(k,55)
         mat(k,1471) = rxt(k,188)*y(k,55)
         mat(k,386) = rxt(k,257)*y(k,55)
         mat(k,2221) = rxt(k,298)*y(k,27) + rxt(k,250)*y(k,40) + rxt(k,187)*y(k,41) &
                      + rxt(k,252)*y(k,42) + 2.000_r8*rxt(k,255)*y(k,45) + rxt(k,188) &
                      *y(k,53) + rxt(k,257)*y(k,54) + rxt(k,189)*y(k,76) + rxt(k,190) &
                      *y(k,78) + rxt(k,209)*y(k,91) + rxt(k,191)*y(k,207)
         mat(k,2248) = rxt(k,206)*y(k,221)
         mat(k,1258) = rxt(k,189)*y(k,55)
         mat(k,602) = rxt(k,190)*y(k,55)
         mat(k,827) = mat(k,827) + rxt(k,209)*y(k,55)
         mat(k,1923) = rxt(k,191)*y(k,55)
         mat(k,1790) = mat(k,1790) + rxt(k,206)*y(k,58)
         mat(k,185) = -(rxt(k,278)*y(k,221) + rxt(k,286)*y(k,220))
         mat(k,1654) = -rxt(k,278)*y(k,85)
         mat(k,1599) = -rxt(k,286)*y(k,85)
         mat(k,868) = -(rxt(k,279)*y(k,221))
         mat(k,1741) = -rxt(k,279)*y(k,86)
         mat(k,935) = .050_r8*rxt(k,452)*y(k,137)
         mat(k,295) = .350_r8*rxt(k,289)*y(k,221)
         mat(k,537) = .370_r8*rxt(k,291)*y(k,137)
         mat(k,1022) = .120_r8*rxt(k,320)*y(k,137)
         mat(k,835) = .110_r8*rxt(k,397)*y(k,137)
         mat(k,1224) = .330_r8*rxt(k,350)*y(k,137)
         mat(k,907) = .050_r8*rxt(k,455)*y(k,137)
         mat(k,1349) = .120_r8*rxt(k,364)*y(k,137)
         mat(k,2300) = rxt(k,282)*y(k,208)
         mat(k,1946) = .050_r8*rxt(k,452)*y(k,5) + .370_r8*rxt(k,291)*y(k,24) &
                      + .120_r8*rxt(k,320)*y(k,28) + .110_r8*rxt(k,397)*y(k,97) &
                      + .330_r8*rxt(k,350)*y(k,104) + .050_r8*rxt(k,455)*y(k,109) &
                      + .120_r8*rxt(k,364)*y(k,110)
         mat(k,1883) = rxt(k,280)*y(k,208)
         mat(k,439) = rxt(k,282)*y(k,123) + rxt(k,280)*y(k,207)
         mat(k,1741) = mat(k,1741) + .350_r8*rxt(k,289)*y(k,23)
         mat(k,1457) = rxt(k,244)*y(k,72)
         mat(k,882) = rxt(k,244)*y(k,53) + rxt(k,245)*y(k,76) + rxt(k,247)*y(k,88) &
                      + rxt(k,246)*y(k,233)
         mat(k,1249) = rxt(k,245)*y(k,72)
         mat(k,1476) = rxt(k,247)*y(k,72)
         mat(k,2350) = rxt(k,246)*y(k,72)
         mat(k,1480) = -(rxt(k,184)*y(k,221) + rxt(k,247)*y(k,72))
         mat(k,1779) = -rxt(k,184)*y(k,88)
         mat(k,886) = -rxt(k,247)*y(k,88)
         mat(k,1526) = rxt(k,268)*y(k,125)
         mat(k,1050) = rxt(k,300)*y(k,125)
         mat(k,1242) = rxt(k,326)*y(k,125)
         mat(k,986) = (rxt(k,540)+rxt(k,545)+rxt(k,550))*y(k,84)
         mat(k,290) = rxt(k,474)*y(k,125)
         mat(k,2055) = (rxt(k,540)+rxt(k,545)+rxt(k,550))*y(k,59)
         mat(k,1571) = rxt(k,183)*y(k,221)
         mat(k,2116) = rxt(k,268)*y(k,41) + rxt(k,300)*y(k,44) + rxt(k,326)*y(k,48) &
                      + rxt(k,474)*y(k,66)
         mat(k,1779) = mat(k,1779) + rxt(k,183)*y(k,124)
         mat(k,469) = -(rxt(k,162)*y(k,221))
         mat(k,1698) = -rxt(k,162)*y(k,89)
         mat(k,1548) = rxt(k,181)*y(k,207)
         mat(k,1855) = rxt(k,181)*y(k,124)
         mat(k,770) = -(rxt(k,238)*y(k,133) + (rxt(k,542) + rxt(k,547)) * y(k,84))
         mat(k,2001) = -rxt(k,238)*y(k,90)
         mat(k,2051) = -(rxt(k,542) + rxt(k,547)) * y(k,90)
         mat(k,1800) = rxt(k,230)*y(k,207)
         mat(k,1875) = rxt(k,230)*y(k,18)
         mat(k,824) = -(rxt(k,209)*y(k,55) + rxt(k,210)*y(k,133) + rxt(k,211)*y(k,221) &
                      + (rxt(k,535) + rxt(k,541) + rxt(k,546)) * y(k,84))
         mat(k,2197) = -rxt(k,209)*y(k,91)
         mat(k,2003) = -rxt(k,210)*y(k,91)
         mat(k,1737) = -rxt(k,211)*y(k,91)
         mat(k,2052) = -(rxt(k,535) + rxt(k,541) + rxt(k,546)) * y(k,91)
         mat(k,2232) = rxt(k,198)*y(k,207)
         mat(k,984) = rxt(k,203)*y(k,221)
         mat(k,1880) = rxt(k,198)*y(k,58)
         mat(k,1737) = mat(k,1737) + rxt(k,203)*y(k,59)
         mat(k,1075) = -(rxt(k,343)*y(k,221))
         mat(k,1757) = -rxt(k,343)*y(k,92)
         mat(k,577) = .300_r8*rxt(k,388)*y(k,221)
         mat(k,561) = .500_r8*rxt(k,389)*y(k,221)
         mat(k,2310) = rxt(k,342)*y(k,204) + rxt(k,349)*y(k,213)
         mat(k,554) = rxt(k,342)*y(k,123)
         mat(k,1329) = rxt(k,349)*y(k,123)
         mat(k,1757) = mat(k,1757) + .300_r8*rxt(k,388)*y(k,98) + .500_r8*rxt(k,389) &
                      *y(k,99)
         mat(k,224) = -(rxt(k,374)*y(k,221))
         mat(k,1661) = -rxt(k,374)*y(k,93)
         mat(k,1089) = -(rxt(k,328)*y(k,221))
         mat(k,1758) = -rxt(k,328)*y(k,94)
         mat(k,578) = .700_r8*rxt(k,388)*y(k,221)
         mat(k,562) = .500_r8*rxt(k,389)*y(k,221)
         mat(k,569) = .500_r8*rxt(k,363)*y(k,221)
         mat(k,2311) = .050_r8*rxt(k,386)*y(k,210) + .220_r8*rxt(k,348)*y(k,213) &
                      + .250_r8*rxt(k,405)*y(k,229)
         mat(k,2096) = .050_r8*rxt(k,387)*y(k,210) + .220_r8*rxt(k,347)*y(k,213) &
                      + .250_r8*rxt(k,406)*y(k,229)
         mat(k,530) = .500_r8*rxt(k,332)*y(k,221)
         mat(k,1401) = .220_r8*rxt(k,344)*y(k,213) + .250_r8*rxt(k,402)*y(k,229)
         mat(k,2151) = .230_r8*rxt(k,345)*y(k,213) + .200_r8*rxt(k,333)*y(k,224) &
                      + .100_r8*rxt(k,403)*y(k,229)
         mat(k,1303) = .050_r8*rxt(k,386)*y(k,123) + .050_r8*rxt(k,387)*y(k,125)
         mat(k,1330) = .220_r8*rxt(k,348)*y(k,123) + .220_r8*rxt(k,347)*y(k,125) &
                      + .220_r8*rxt(k,344)*y(k,201) + .230_r8*rxt(k,345)*y(k,202)
         mat(k,1758) = mat(k,1758) + .700_r8*rxt(k,388)*y(k,98) + .500_r8*rxt(k,389) &
                      *y(k,99) + .500_r8*rxt(k,363)*y(k,108) + .500_r8*rxt(k,332) &
                      *y(k,149)
         mat(k,1163) = .200_r8*rxt(k,333)*y(k,202)
         mat(k,1184) = .250_r8*rxt(k,405)*y(k,123) + .250_r8*rxt(k,406)*y(k,125) &
                      + .250_r8*rxt(k,402)*y(k,201) + .100_r8*rxt(k,403)*y(k,202)
         mat(k,344) = -(rxt(k,375)*y(k,221))
         mat(k,1680) = -rxt(k,375)*y(k,95)
         mat(k,2271) = .870_r8*rxt(k,386)*y(k,210)
         mat(k,2077) = .950_r8*rxt(k,387)*y(k,210)
         mat(k,1394) = rxt(k,382)*y(k,210)
         mat(k,2135) = .750_r8*rxt(k,383)*y(k,210)
         mat(k,1293) = .870_r8*rxt(k,386)*y(k,123) + .950_r8*rxt(k,387)*y(k,125) &
                      + rxt(k,382)*y(k,201) + .750_r8*rxt(k,383)*y(k,202)
         mat(k,136) = -(rxt(k,376)*y(k,221))
         mat(k,1650) = -rxt(k,376)*y(k,96)
         mat(k,693) = .600_r8*rxt(k,399)*y(k,221)
         mat(k,1650) = mat(k,1650) + .600_r8*rxt(k,399)*y(k,102)
         mat(k,834) = -(rxt(k,390)*y(k,125) + rxt(k,397)*y(k,137) + rxt(k,398) &
                      *y(k,221))
         mat(k,2080) = -rxt(k,390)*y(k,97)
         mat(k,1945) = -rxt(k,397)*y(k,97)
         mat(k,1738) = -rxt(k,398)*y(k,97)
         mat(k,576) = -(rxt(k,388)*y(k,221))
         mat(k,1711) = -rxt(k,388)*y(k,98)
         mat(k,2284) = .080_r8*rxt(k,380)*y(k,209)
         mat(k,1263) = .080_r8*rxt(k,380)*y(k,123)
         mat(k,559) = -(rxt(k,389)*y(k,221))
         mat(k,1709) = -rxt(k,389)*y(k,99)
         mat(k,2283) = .080_r8*rxt(k,386)*y(k,210)
         mat(k,1294) = .080_r8*rxt(k,386)*y(k,123)
         mat(k,1209) = -(rxt(k,391)*y(k,201) + rxt(k,392)*y(k,202) + rxt(k,393) &
                      *y(k,207) + rxt(k,394)*y(k,123) + rxt(k,395)*y(k,125))
         mat(k,1405) = -rxt(k,391)*y(k,100)
         mat(k,2159) = -rxt(k,392)*y(k,100)
         mat(k,1900) = -rxt(k,393)*y(k,100)
         mat(k,2319) = -rxt(k,394)*y(k,100)
         mat(k,2104) = -rxt(k,395)*y(k,100)
         mat(k,838) = rxt(k,390)*y(k,125)
         mat(k,2104) = mat(k,2104) + rxt(k,390)*y(k,97)
         mat(k,400) = -(rxt(k,396)*y(k,221))
         mat(k,1688) = -rxt(k,396)*y(k,101)
         mat(k,1199) = rxt(k,393)*y(k,207)
         mat(k,1844) = rxt(k,393)*y(k,100)
         mat(k,694) = -(rxt(k,399)*y(k,221))
         mat(k,1725) = -rxt(k,399)*y(k,102)
         mat(k,1869) = rxt(k,379)*y(k,209) + rxt(k,384)*y(k,210)
         mat(k,1264) = rxt(k,379)*y(k,207)
         mat(k,1296) = rxt(k,384)*y(k,207)
         mat(k,75) = -(rxt(k,520)*y(k,221))
         mat(k,1642) = -rxt(k,520)*y(k,103)
         mat(k,1227) = -(rxt(k,350)*y(k,137) + rxt(k,351)*y(k,221))
         mat(k,1966) = -rxt(k,350)*y(k,104)
         mat(k,1767) = -rxt(k,351)*y(k,104)
         mat(k,839) = .300_r8*rxt(k,397)*y(k,137)
         mat(k,2320) = .360_r8*rxt(k,380)*y(k,209)
         mat(k,2105) = .400_r8*rxt(k,381)*y(k,209)
         mat(k,1966) = mat(k,1966) + .300_r8*rxt(k,397)*y(k,97)
         mat(k,1406) = .390_r8*rxt(k,377)*y(k,209)
         mat(k,2160) = .310_r8*rxt(k,378)*y(k,209)
         mat(k,1273) = .360_r8*rxt(k,380)*y(k,123) + .400_r8*rxt(k,381)*y(k,125) &
                      + .390_r8*rxt(k,377)*y(k,201) + .310_r8*rxt(k,378)*y(k,202)
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
         mat(k,313) = -(rxt(k,352)*y(k,221))
         mat(k,1676) = -rxt(k,352)*y(k,105)
         mat(k,1837) = rxt(k,346)*y(k,213)
         mat(k,1327) = rxt(k,346)*y(k,207)
         mat(k,507) = -(rxt(k,361)*y(k,221))
         mat(k,1703) = -rxt(k,361)*y(k,106)
         mat(k,2280) = .800_r8*rxt(k,370)*y(k,193)
         mat(k,955) = .800_r8*rxt(k,370)*y(k,123)
         mat(k,318) = -(rxt(k,362)*y(k,221))
         mat(k,1677) = -rxt(k,362)*y(k,107)
         mat(k,1838) = .800_r8*rxt(k,359)*y(k,217)
         mat(k,671) = .800_r8*rxt(k,359)*y(k,207)
         mat(k,567) = -(rxt(k,363)*y(k,221))
         mat(k,1710) = -rxt(k,363)*y(k,108)
         mat(k,1552) = rxt(k,366)*y(k,215)
         mat(k,1373) = rxt(k,366)*y(k,124)
         mat(k,908) = -(rxt(k,454)*y(k,125) + rxt(k,455)*y(k,137) + rxt(k,456) &
                      *y(k,221))
         mat(k,2084) = -rxt(k,454)*y(k,109)
         mat(k,1948) = -rxt(k,455)*y(k,109)
         mat(k,1745) = -rxt(k,456)*y(k,109)
         mat(k,1357) = -(rxt(k,364)*y(k,137) + rxt(k,365)*y(k,221))
         mat(k,1971) = -rxt(k,364)*y(k,110)
         mat(k,1773) = -rxt(k,365)*y(k,110)
         mat(k,842) = .200_r8*rxt(k,397)*y(k,137)
         mat(k,2325) = .560_r8*rxt(k,380)*y(k,209)
         mat(k,2110) = .600_r8*rxt(k,381)*y(k,209)
         mat(k,1971) = mat(k,1971) + .200_r8*rxt(k,397)*y(k,97)
         mat(k,1411) = .610_r8*rxt(k,377)*y(k,209)
         mat(k,2165) = .440_r8*rxt(k,378)*y(k,209)
         mat(k,1277) = .560_r8*rxt(k,380)*y(k,123) + .600_r8*rxt(k,381)*y(k,125) &
                      + .610_r8*rxt(k,377)*y(k,201) + .440_r8*rxt(k,378)*y(k,202)
         mat(k,680) = -(rxt(k,163)*y(k,123) + (rxt(k,164) + rxt(k,165) + rxt(k,166) &
                      ) * y(k,124) + rxt(k,167)*y(k,134) + rxt(k,175)*y(k,221))
         mat(k,2290) = -rxt(k,163)*y(k,111)
         mat(k,1554) = -(rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,111)
         mat(k,1497) = -rxt(k,167)*y(k,111)
         mat(k,1723) = -rxt(k,175)*y(k,111)
         mat(k,259) = -((rxt(k,179) + rxt(k,180)) * y(k,220))
         mat(k,1604) = -(rxt(k,179) + rxt(k,180)) * y(k,112)
         mat(k,679) = rxt(k,164)*y(k,124)
         mat(k,1545) = rxt(k,164)*y(k,111)
         mat(k,1546) = rxt(k,182)*y(k,125)
         mat(k,2075) = rxt(k,182)*y(k,124)
         mat(k,424) = -(rxt(k,400)*y(k,221))
         mat(k,1692) = -rxt(k,400)*y(k,114)
         mat(k,1200) = .200_r8*rxt(k,392)*y(k,202)
         mat(k,2136) = .200_r8*rxt(k,392)*y(k,100)
         mat(k,1007) = -(rxt(k,401)*y(k,221))
         mat(k,1752) = -rxt(k,401)*y(k,115)
         mat(k,1204) = rxt(k,394)*y(k,123) + rxt(k,395)*y(k,125) + rxt(k,391)*y(k,201) &
                      + .800_r8*rxt(k,392)*y(k,202)
         mat(k,2306) = rxt(k,394)*y(k,100)
         mat(k,2090) = rxt(k,395)*y(k,100)
         mat(k,1399) = rxt(k,391)*y(k,100)
         mat(k,2147) = .800_r8*rxt(k,392)*y(k,100)
         mat(k,97) = -(rxt(k,491)*y(k,221))
         mat(k,1646) = -rxt(k,491)*y(k,119)
         mat(k,2346) = -(rxt(k,163)*y(k,111) + rxt(k,172)*y(k,125) + rxt(k,176) &
                      *y(k,207) + rxt(k,177)*y(k,137) + rxt(k,178)*y(k,133) + rxt(k,199) &
                      *y(k,58) + rxt(k,231)*y(k,18) + rxt(k,274)*y(k,202) + rxt(k,282) &
                      *y(k,208) + rxt(k,295)*y(k,198) + rxt(k,306)*y(k,201) + rxt(k,310) &
                      *y(k,206) + rxt(k,323)*y(k,199) + rxt(k,331)*y(k,223) + rxt(k,335) &
                      *y(k,224) + (rxt(k,341) + rxt(k,342)) * y(k,204) + (rxt(k,348) &
                      + rxt(k,349)) * y(k,213) + rxt(k,357)*y(k,215) + rxt(k,360) &
                      *y(k,217) + (rxt(k,370) + rxt(k,371)) * y(k,193) + rxt(k,380) &
                      *y(k,209) + rxt(k,386)*y(k,210) + rxt(k,394)*y(k,100) + rxt(k,405) &
                      *y(k,229) + rxt(k,409)*y(k,192) + rxt(k,412)*y(k,195) + rxt(k,417) &
                      *y(k,197) + rxt(k,419)*y(k,200) + rxt(k,423)*y(k,203) + rxt(k,426) &
                      *y(k,214) + rxt(k,429)*y(k,216) + rxt(k,432)*y(k,222) + rxt(k,439) &
                      *y(k,227) + rxt(k,445)*y(k,230) + rxt(k,448)*y(k,232) + rxt(k,459) &
                      *y(k,219) + rxt(k,464)*y(k,225) + rxt(k,469)*y(k,226))
         mat(k,687) = -rxt(k,163)*y(k,123)
         mat(k,2132) = -rxt(k,172)*y(k,123)
         mat(k,1928) = -rxt(k,176)*y(k,123)
         mat(k,1993) = -rxt(k,177)*y(k,123)
         mat(k,2027) = -rxt(k,178)*y(k,123)
         mat(k,2253) = -rxt(k,199)*y(k,123)
         mat(k,1820) = -rxt(k,231)*y(k,123)
         mat(k,2185) = -rxt(k,274)*y(k,123)
         mat(k,443) = -rxt(k,282)*y(k,123)
         mat(k,867) = -rxt(k,295)*y(k,123)
         mat(k,1425) = -rxt(k,306)*y(k,123)
         mat(k,787) = -rxt(k,310)*y(k,123)
         mat(k,813) = -rxt(k,323)*y(k,123)
         mat(k,797) = -rxt(k,331)*y(k,123)
         mat(k,1172) = -rxt(k,335)*y(k,123)
         mat(k,558) = -(rxt(k,341) + rxt(k,342)) * y(k,123)
         mat(k,1347) = -(rxt(k,348) + rxt(k,349)) * y(k,123)
         mat(k,1392) = -rxt(k,357)*y(k,123)
         mat(k,678) = -rxt(k,360)*y(k,123)
         mat(k,970) = -(rxt(k,370) + rxt(k,371)) * y(k,123)
         mat(k,1290) = -rxt(k,380)*y(k,123)
         mat(k,1325) = -rxt(k,386)*y(k,123)
         mat(k,1223) = -rxt(k,394)*y(k,123)
         mat(k,1198) = -rxt(k,405)*y(k,123)
         mat(k,518) = -rxt(k,409)*y(k,123)
         mat(k,498) = -rxt(k,412)*y(k,123)
         mat(k,437) = -rxt(k,417)*y(k,123)
         mat(k,625) = -rxt(k,419)*y(k,123)
         mat(k,768) = -rxt(k,423)*y(k,123)
         mat(k,728) = -rxt(k,426)*y(k,123)
         mat(k,881) = -rxt(k,429)*y(k,123)
         mat(k,450) = -rxt(k,432)*y(k,123)
         mat(k,743) = -rxt(k,439)*y(k,123)
         mat(k,760) = -rxt(k,445)*y(k,123)
         mat(k,506) = -rxt(k,448)*y(k,123)
         mat(k,1071) = -rxt(k,459)*y(k,123)
         mat(k,1159) = -rxt(k,464)*y(k,123)
         mat(k,1138) = -rxt(k,469)*y(k,123)
         mat(k,687) = mat(k,687) + 2.000_r8*rxt(k,165)*y(k,124) + rxt(k,167)*y(k,134) &
                      + rxt(k,175)*y(k,221)
         mat(k,262) = 2.000_r8*rxt(k,179)*y(k,220)
         mat(k,1587) = 2.000_r8*rxt(k,165)*y(k,111) + rxt(k,168)*y(k,133) + rxt(k,484) &
                      *y(k,153)
         mat(k,2027) = mat(k,2027) + rxt(k,168)*y(k,124)
         mat(k,1516) = rxt(k,167)*y(k,111)
         mat(k,1442) = rxt(k,484)*y(k,124)
         mat(k,1630) = 2.000_r8*rxt(k,179)*y(k,112)
         mat(k,1795) = rxt(k,175)*y(k,111)
         mat(k,1574) = -((rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,111) + (rxt(k,168) &
                      + rxt(k,170)) * y(k,133) + rxt(k,169)*y(k,137) + rxt(k,181) &
                      *y(k,207) + rxt(k,182)*y(k,125) + rxt(k,183)*y(k,221) + rxt(k,201) &
                      *y(k,58) + rxt(k,232)*y(k,18) + rxt(k,317)*y(k,201) + rxt(k,366) &
                      *y(k,215) + rxt(k,424)*y(k,203) + rxt(k,427)*y(k,214) + rxt(k,430) &
                      *y(k,216) + rxt(k,434)*y(k,144) + rxt(k,437)*y(k,192) + rxt(k,484) &
                      *y(k,153))
         mat(k,682) = -(rxt(k,164) + rxt(k,165) + rxt(k,166)) * y(k,124)
         mat(k,2014) = -(rxt(k,168) + rxt(k,170)) * y(k,124)
         mat(k,1980) = -rxt(k,169)*y(k,124)
         mat(k,1915) = -rxt(k,181)*y(k,124)
         mat(k,2119) = -rxt(k,182)*y(k,124)
         mat(k,1782) = -rxt(k,183)*y(k,124)
         mat(k,2240) = -rxt(k,201)*y(k,124)
         mat(k,1807) = -rxt(k,232)*y(k,124)
         mat(k,1418) = -rxt(k,317)*y(k,124)
         mat(k,1385) = -rxt(k,366)*y(k,124)
         mat(k,764) = -rxt(k,424)*y(k,124)
         mat(k,726) = -rxt(k,427)*y(k,124)
         mat(k,877) = -rxt(k,430)*y(k,124)
         mat(k,467) = -rxt(k,434)*y(k,124)
         mat(k,515) = -rxt(k,437)*y(k,124)
         mat(k,1433) = -rxt(k,484)*y(k,124)
         mat(k,668) = rxt(k,368)*y(k,221)
         mat(k,357) = rxt(k,339)*y(k,125)
         mat(k,1807) = mat(k,1807) + rxt(k,231)*y(k,123)
         mat(k,2240) = mat(k,2240) + rxt(k,199)*y(k,123)
         mat(k,471) = rxt(k,162)*y(k,221)
         mat(k,582) = .700_r8*rxt(k,388)*y(k,221)
         mat(k,1217) = rxt(k,394)*y(k,123) + rxt(k,395)*y(k,125)
         mat(k,2333) = rxt(k,231)*y(k,18) + rxt(k,199)*y(k,58) + rxt(k,394)*y(k,100) &
                      + 2.000_r8*rxt(k,172)*y(k,125) + rxt(k,178)*y(k,133) &
                      + rxt(k,177)*y(k,137) + rxt(k,409)*y(k,192) + rxt(k,370) &
                      *y(k,193) + rxt(k,412)*y(k,195) + rxt(k,417)*y(k,197) &
                      + rxt(k,295)*y(k,198) + rxt(k,323)*y(k,199) + rxt(k,419) &
                      *y(k,200) + rxt(k,306)*y(k,201) + rxt(k,274)*y(k,202) &
                      + rxt(k,423)*y(k,203) + rxt(k,341)*y(k,204) + rxt(k,310) &
                      *y(k,206) + rxt(k,176)*y(k,207) + rxt(k,282)*y(k,208) &
                      + .920_r8*rxt(k,380)*y(k,209) + .920_r8*rxt(k,386)*y(k,210) &
                      + rxt(k,348)*y(k,213) + rxt(k,426)*y(k,214) + rxt(k,357) &
                      *y(k,215) + rxt(k,429)*y(k,216) + rxt(k,360)*y(k,217) &
                      + 1.600_r8*rxt(k,459)*y(k,219) + rxt(k,432)*y(k,222) &
                      + rxt(k,331)*y(k,223) + rxt(k,335)*y(k,224) + .900_r8*rxt(k,464) &
                      *y(k,225) + .800_r8*rxt(k,469)*y(k,226) + rxt(k,439)*y(k,227) &
                      + rxt(k,405)*y(k,229) + rxt(k,445)*y(k,230) + rxt(k,448) &
                      *y(k,232)
         mat(k,2119) = mat(k,2119) + rxt(k,339)*y(k,15) + rxt(k,395)*y(k,100) &
                      + 2.000_r8*rxt(k,172)*y(k,123) + rxt(k,173)*y(k,133) &
                      + rxt(k,171)*y(k,207) + rxt(k,381)*y(k,209) + rxt(k,387) &
                      *y(k,210) + rxt(k,347)*y(k,213) + rxt(k,358)*y(k,215) &
                      + 2.000_r8*rxt(k,460)*y(k,219) + rxt(k,174)*y(k,221) &
                      + rxt(k,406)*y(k,229)
         mat(k,854) = rxt(k,329)*y(k,221)
         mat(k,2014) = mat(k,2014) + rxt(k,178)*y(k,123) + rxt(k,173)*y(k,125)
         mat(k,1980) = mat(k,1980) + rxt(k,177)*y(k,123)
         mat(k,615) = rxt(k,466)*y(k,221)
         mat(k,515) = mat(k,515) + rxt(k,409)*y(k,123)
         mat(k,966) = rxt(k,370)*y(k,123)
         mat(k,495) = rxt(k,412)*y(k,123)
         mat(k,434) = rxt(k,417)*y(k,123)
         mat(k,863) = rxt(k,295)*y(k,123)
         mat(k,809) = rxt(k,323)*y(k,123)
         mat(k,621) = rxt(k,419)*y(k,123)
         mat(k,1418) = mat(k,1418) + rxt(k,306)*y(k,123)
         mat(k,2172) = rxt(k,274)*y(k,123) + .500_r8*rxt(k,457)*y(k,219)
         mat(k,764) = mat(k,764) + rxt(k,423)*y(k,123)
         mat(k,556) = rxt(k,341)*y(k,123)
         mat(k,784) = rxt(k,310)*y(k,123)
         mat(k,1915) = mat(k,1915) + rxt(k,176)*y(k,123) + rxt(k,171)*y(k,125)
         mat(k,441) = rxt(k,282)*y(k,123)
         mat(k,1283) = .920_r8*rxt(k,380)*y(k,123) + rxt(k,381)*y(k,125)
         mat(k,1318) = .920_r8*rxt(k,386)*y(k,123) + rxt(k,387)*y(k,125)
         mat(k,1341) = rxt(k,348)*y(k,123) + rxt(k,347)*y(k,125)
         mat(k,726) = mat(k,726) + rxt(k,426)*y(k,123)
         mat(k,1385) = mat(k,1385) + rxt(k,357)*y(k,123) + rxt(k,358)*y(k,125)
         mat(k,877) = mat(k,877) + rxt(k,429)*y(k,123)
         mat(k,675) = rxt(k,360)*y(k,123)
         mat(k,1066) = 1.600_r8*rxt(k,459)*y(k,123) + 2.000_r8*rxt(k,460)*y(k,125) &
                      + .500_r8*rxt(k,457)*y(k,202)
         mat(k,1782) = mat(k,1782) + rxt(k,368)*y(k,1) + rxt(k,162)*y(k,89) &
                      + .700_r8*rxt(k,388)*y(k,98) + rxt(k,174)*y(k,125) + rxt(k,329) &
                      *y(k,126) + rxt(k,466)*y(k,179)
         mat(k,447) = rxt(k,432)*y(k,123)
         mat(k,794) = rxt(k,331)*y(k,123)
         mat(k,1168) = rxt(k,335)*y(k,123)
         mat(k,1154) = .900_r8*rxt(k,464)*y(k,123)
         mat(k,1132) = .800_r8*rxt(k,469)*y(k,123)
         mat(k,740) = rxt(k,439)*y(k,123)
         mat(k,1192) = rxt(k,405)*y(k,123) + rxt(k,406)*y(k,125)
         mat(k,757) = rxt(k,445)*y(k,123)
         mat(k,503) = rxt(k,448)*y(k,123)
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
         mat(k,2128) = -(rxt(k,171)*y(k,207) + rxt(k,172)*y(k,123) + rxt(k,173) &
                      *y(k,133) + rxt(k,174)*y(k,221) + rxt(k,182)*y(k,124) + rxt(k,268) &
                      *y(k,41) + rxt(k,300)*y(k,44) + rxt(k,319)*y(k,28) + rxt(k,326) &
                      *y(k,48) + rxt(k,339)*y(k,15) + rxt(k,347)*y(k,213) + rxt(k,358) &
                      *y(k,215) + rxt(k,381)*y(k,209) + rxt(k,387)*y(k,210) + rxt(k,390) &
                      *y(k,97) + rxt(k,395)*y(k,100) + rxt(k,406)*y(k,229) + rxt(k,451) &
                      *y(k,5) + rxt(k,454)*y(k,109) + rxt(k,460)*y(k,219) + rxt(k,471) &
                      *y(k,181) + rxt(k,474)*y(k,66))
         mat(k,1924) = -rxt(k,171)*y(k,125)
         mat(k,2342) = -rxt(k,172)*y(k,125)
         mat(k,2023) = -rxt(k,173)*y(k,125)
         mat(k,1791) = -rxt(k,174)*y(k,125)
         mat(k,1583) = -rxt(k,182)*y(k,125)
         mat(k,1538) = -rxt(k,268)*y(k,125)
         mat(k,1053) = -rxt(k,300)*y(k,125)
         mat(k,1037) = -rxt(k,319)*y(k,125)
         mat(k,1247) = -rxt(k,326)*y(k,125)
         mat(k,359) = -rxt(k,339)*y(k,125)
         mat(k,1345) = -rxt(k,347)*y(k,125)
         mat(k,1390) = -rxt(k,358)*y(k,125)
         mat(k,1288) = -rxt(k,381)*y(k,125)
         mat(k,1323) = -rxt(k,387)*y(k,125)
         mat(k,848) = -rxt(k,390)*y(k,125)
         mat(k,1221) = -rxt(k,395)*y(k,125)
         mat(k,1196) = -rxt(k,406)*y(k,125)
         mat(k,953) = -rxt(k,451)*y(k,125)
         mat(k,925) = -rxt(k,454)*y(k,125)
         mat(k,1069) = -rxt(k,460)*y(k,125)
         mat(k,1002) = -rxt(k,471)*y(k,125)
         mat(k,293) = -rxt(k,474)*y(k,125)
         mat(k,550) = rxt(k,233)*y(k,133)
         mat(k,2222) = rxt(k,200)*y(k,59)
         mat(k,991) = rxt(k,200)*y(k,55) + rxt(k,202)*y(k,133) + rxt(k,203)*y(k,221)
         mat(k,889) = rxt(k,247)*y(k,88)
         mat(k,1490) = rxt(k,247)*y(k,72) + rxt(k,184)*y(k,221)
         mat(k,575) = .500_r8*rxt(k,363)*y(k,221)
         mat(k,1583) = mat(k,1583) + rxt(k,170)*y(k,133) + rxt(k,169)*y(k,137)
         mat(k,2023) = mat(k,2023) + rxt(k,233)*y(k,19) + rxt(k,202)*y(k,59) &
                      + rxt(k,170)*y(k,124)
         mat(k,1989) = rxt(k,169)*y(k,124)
         mat(k,525) = rxt(k,315)*y(k,221)
         mat(k,1791) = mat(k,1791) + rxt(k,203)*y(k,59) + rxt(k,184)*y(k,88) &
                      + .500_r8*rxt(k,363)*y(k,108) + rxt(k,315)*y(k,142)
         mat(k,850) = -(rxt(k,329)*y(k,221))
         mat(k,1739) = -rxt(k,329)*y(k,126)
         mat(k,1021) = rxt(k,319)*y(k,125)
         mat(k,560) = .500_r8*rxt(k,389)*y(k,221)
         mat(k,402) = rxt(k,396)*y(k,221)
         mat(k,425) = rxt(k,400)*y(k,221)
         mat(k,1004) = rxt(k,401)*y(k,221)
         mat(k,2081) = rxt(k,319)*y(k,28)
         mat(k,1739) = mat(k,1739) + .500_r8*rxt(k,389)*y(k,99) + rxt(k,396)*y(k,101) &
                      + rxt(k,400)*y(k,114) + rxt(k,401)*y(k,115)
         mat(k,388) = -(rxt(k,461)*y(k,221))
         mat(k,1686) = -rxt(k,461)*y(k,127)
         mat(k,1842) = rxt(k,458)*y(k,219)
         mat(k,1056) = rxt(k,458)*y(k,207)
         mat(k,2020) = -(rxt(k,142)*y(k,137) + 4._r8*rxt(k,143)*y(k,133) + rxt(k,144) &
                      *y(k,134) + rxt(k,145)*y(k,76) + rxt(k,146)*y(k,78) + rxt(k,151) &
                      *y(k,207) + rxt(k,157)*y(k,221) + (rxt(k,168) + rxt(k,170) &
                      ) * y(k,124) + rxt(k,173)*y(k,125) + rxt(k,178)*y(k,123) &
                      + rxt(k,202)*y(k,59) + rxt(k,204)*y(k,58) + rxt(k,207)*y(k,84) &
                      + rxt(k,210)*y(k,91) + rxt(k,233)*y(k,19) + rxt(k,234)*y(k,18) &
                      + rxt(k,236)*y(k,80) + rxt(k,238)*y(k,90) + rxt(k,269)*y(k,41) &
                      + rxt(k,476)*y(k,140))
         mat(k,1986) = -rxt(k,142)*y(k,133)
         mat(k,1512) = -rxt(k,144)*y(k,133)
         mat(k,1256) = -rxt(k,145)*y(k,133)
         mat(k,601) = -rxt(k,146)*y(k,133)
         mat(k,1921) = -rxt(k,151)*y(k,133)
         mat(k,1788) = -rxt(k,157)*y(k,133)
         mat(k,1580) = -(rxt(k,168) + rxt(k,170)) * y(k,133)
         mat(k,2125) = -rxt(k,173)*y(k,133)
         mat(k,2339) = -rxt(k,178)*y(k,133)
         mat(k,989) = -rxt(k,202)*y(k,133)
         mat(k,2246) = -rxt(k,204)*y(k,133)
         mat(k,2064) = -rxt(k,207)*y(k,133)
         mat(k,826) = -rxt(k,210)*y(k,133)
         mat(k,549) = -rxt(k,233)*y(k,133)
         mat(k,1813) = -rxt(k,234)*y(k,133)
         mat(k,820) = -rxt(k,236)*y(k,133)
         mat(k,774) = -rxt(k,238)*y(k,133)
         mat(k,1535) = -rxt(k,269)*y(k,133)
         mat(k,366) = -rxt(k,476)*y(k,133)
         mat(k,2040) = rxt(k,149)*y(k,207)
         mat(k,685) = rxt(k,163)*y(k,123) + rxt(k,164)*y(k,124) + rxt(k,167)*y(k,134)
         mat(k,2339) = mat(k,2339) + rxt(k,163)*y(k,111)
         mat(k,1580) = mat(k,1580) + rxt(k,164)*y(k,111)
         mat(k,1512) = mat(k,1512) + rxt(k,167)*y(k,111) + rxt(k,478)*y(k,151) &
                      + rxt(k,485)*y(k,153) + (rxt(k,130)+rxt(k,131))*y(k,220)
         mat(k,1986) = mat(k,1986) + 2.000_r8*rxt(k,133)*y(k,220)
         mat(k,721) = rxt(k,478)*y(k,134)
         mat(k,1438) = rxt(k,485)*y(k,134)
         mat(k,1921) = mat(k,1921) + rxt(k,149)*y(k,75)
         mat(k,1623) = (rxt(k,130)+rxt(k,131))*y(k,134) + 2.000_r8*rxt(k,133)*y(k,137)
         mat(k,1788) = mat(k,1788) + 2.000_r8*rxt(k,159)*y(k,221)
         mat(k,1504) = -(rxt(k,130)*y(k,220) + rxt(k,136)*y(k,135) + rxt(k,144) &
                      *y(k,133) + rxt(k,150)*y(k,75) + rxt(k,167)*y(k,111) + rxt(k,312) &
                      *y(k,205) + rxt(k,478)*y(k,151) + rxt(k,485)*y(k,153))
         mat(k,1615) = -rxt(k,130)*y(k,134)
         mat(k,163) = -rxt(k,136)*y(k,134)
         mat(k,2012) = -rxt(k,144)*y(k,134)
         mat(k,2032) = -rxt(k,150)*y(k,134)
         mat(k,681) = -rxt(k,167)*y(k,134)
         mat(k,459) = -rxt(k,312)*y(k,134)
         mat(k,718) = -rxt(k,478)*y(k,134)
         mat(k,1432) = -rxt(k,485)*y(k,134)
         mat(k,1447) = rxt(k,225)*y(k,137) + rxt(k,224)*y(k,207)
         mat(k,1805) = 2.000_r8*rxt(k,226)*y(k,18) + (rxt(k,228)+rxt(k,229))*y(k,58) &
                      + rxt(k,234)*y(k,133) + rxt(k,230)*y(k,207)
         mat(k,2211) = rxt(k,193)*y(k,137) + rxt(k,191)*y(k,207)
         mat(k,2238) = (rxt(k,228)+rxt(k,229))*y(k,18) + (2.000_r8*rxt(k,195) &
                       +2.000_r8*rxt(k,196))*y(k,58) + rxt(k,204)*y(k,133) &
                      + rxt(k,198)*y(k,207) + rxt(k,206)*y(k,221)
         mat(k,2032) = mat(k,2032) + rxt(k,153)*y(k,137) + rxt(k,147)*y(k,207)
         mat(k,470) = rxt(k,162)*y(k,221)
         mat(k,681) = mat(k,681) + rxt(k,166)*y(k,124)
         mat(k,260) = rxt(k,180)*y(k,220)
         mat(k,2331) = rxt(k,177)*y(k,137)
         mat(k,1572) = rxt(k,166)*y(k,111) + rxt(k,168)*y(k,133) + rxt(k,169)*y(k,137)
         mat(k,2117) = rxt(k,173)*y(k,133) + rxt(k,171)*y(k,207)
         mat(k,2012) = mat(k,2012) + rxt(k,234)*y(k,18) + rxt(k,204)*y(k,58) &
                      + rxt(k,168)*y(k,124) + rxt(k,173)*y(k,125) &
                      + 2.000_r8*rxt(k,143)*y(k,133) + rxt(k,135)*y(k,135) &
                      + 2.000_r8*rxt(k,142)*y(k,137) + rxt(k,151)*y(k,207) &
                      + rxt(k,157)*y(k,221)
         mat(k,1504) = mat(k,1504) + 2.000_r8*rxt(k,136)*y(k,135)
         mat(k,163) = mat(k,163) + rxt(k,135)*y(k,133) + 2.000_r8*rxt(k,136)*y(k,134)
         mat(k,1978) = rxt(k,225)*y(k,16) + rxt(k,193)*y(k,55) + rxt(k,153)*y(k,75) &
                      + rxt(k,177)*y(k,123) + rxt(k,169)*y(k,124) &
                      + 2.000_r8*rxt(k,142)*y(k,133) + rxt(k,480)*y(k,151) &
                      + rxt(k,486)*y(k,153) + 2.000_r8*rxt(k,152)*y(k,207) + ( &
                      + 2.000_r8*rxt(k,132)+rxt(k,133))*y(k,220) + rxt(k,158)*y(k,221)
         mat(k,718) = mat(k,718) + rxt(k,480)*y(k,137)
         mat(k,1432) = mat(k,1432) + rxt(k,486)*y(k,137)
         mat(k,861) = rxt(k,294)*y(k,207)
         mat(k,807) = rxt(k,322)*y(k,207)
         mat(k,2170) = rxt(k,273)*y(k,207)
         mat(k,1913) = rxt(k,224)*y(k,16) + rxt(k,230)*y(k,18) + rxt(k,191)*y(k,55) &
                      + rxt(k,198)*y(k,58) + rxt(k,147)*y(k,75) + rxt(k,171)*y(k,125) &
                      + rxt(k,151)*y(k,133) + 2.000_r8*rxt(k,152)*y(k,137) &
                      + rxt(k,294)*y(k,198) + rxt(k,322)*y(k,199) + rxt(k,273) &
                      *y(k,202) + 2.000_r8*rxt(k,161)*y(k,207) + rxt(k,156)*y(k,221) &
                      + rxt(k,330)*y(k,223)
         mat(k,1615) = mat(k,1615) + rxt(k,180)*y(k,112) + (2.000_r8*rxt(k,132) &
                       +rxt(k,133))*y(k,137)
         mat(k,1780) = rxt(k,206)*y(k,58) + rxt(k,162)*y(k,89) + rxt(k,157)*y(k,133) &
                      + rxt(k,158)*y(k,137) + rxt(k,156)*y(k,207)
         mat(k,792) = rxt(k,330)*y(k,207)
         mat(k,162) = -(rxt(k,135)*y(k,133) + rxt(k,136)*y(k,134))
         mat(k,1995) = -rxt(k,135)*y(k,135)
         mat(k,1494) = -rxt(k,136)*y(k,135)
         mat(k,1041) = rxt(k,137)*y(k,136)
         mat(k,1995) = mat(k,1995) + rxt(k,139)*y(k,136)
         mat(k,1494) = mat(k,1494) + rxt(k,140)*y(k,136)
         mat(k,164) = rxt(k,137)*y(k,62) + rxt(k,139)*y(k,133) + rxt(k,140)*y(k,134) &
                      + rxt(k,141)*y(k,137)
         mat(k,1935) = rxt(k,141)*y(k,136)
         mat(k,165) = -(rxt(k,137)*y(k,62) + rxt(k,139)*y(k,133) + rxt(k,140)*y(k,134) &
                      + rxt(k,141)*y(k,137))
         mat(k,1042) = -rxt(k,137)*y(k,136)
         mat(k,1996) = -rxt(k,139)*y(k,136)
         mat(k,1495) = -rxt(k,140)*y(k,136)
         mat(k,1936) = -rxt(k,141)*y(k,136)
         mat(k,1495) = mat(k,1495) + rxt(k,130)*y(k,220)
         mat(k,1598) = rxt(k,130)*y(k,134)
         mat(k,1985) = -((rxt(k,132) + rxt(k,133)) * y(k,220) + rxt(k,142)*y(k,133) &
                      + rxt(k,152)*y(k,207) + rxt(k,153)*y(k,75) + rxt(k,158)*y(k,221) &
                      + rxt(k,169)*y(k,124) + rxt(k,177)*y(k,123) + rxt(k,193)*y(k,55) &
                      + rxt(k,225)*y(k,16) + rxt(k,291)*y(k,24) + rxt(k,320)*y(k,28) &
                      + rxt(k,350)*y(k,104) + rxt(k,364)*y(k,110) + rxt(k,397)*y(k,97) &
                      + rxt(k,435)*y(k,144) + rxt(k,452)*y(k,5) + rxt(k,455)*y(k,109) &
                      + rxt(k,480)*y(k,151) + rxt(k,486)*y(k,153))
         mat(k,1622) = -(rxt(k,132) + rxt(k,133)) * y(k,137)
         mat(k,2019) = -rxt(k,142)*y(k,137)
         mat(k,1920) = -rxt(k,152)*y(k,137)
         mat(k,2039) = -rxt(k,153)*y(k,137)
         mat(k,1787) = -rxt(k,158)*y(k,137)
         mat(k,1579) = -rxt(k,169)*y(k,137)
         mat(k,2338) = -rxt(k,177)*y(k,137)
         mat(k,2218) = -rxt(k,193)*y(k,137)
         mat(k,1453) = -rxt(k,225)*y(k,137)
         mat(k,542) = -rxt(k,291)*y(k,137)
         mat(k,1036) = -rxt(k,320)*y(k,137)
         mat(k,1236) = -rxt(k,350)*y(k,137)
         mat(k,1367) = -rxt(k,364)*y(k,137)
         mat(k,847) = -rxt(k,397)*y(k,137)
         mat(k,468) = -rxt(k,435)*y(k,137)
         mat(k,952) = -rxt(k,452)*y(k,137)
         mat(k,924) = -rxt(k,455)*y(k,137)
         mat(k,720) = -rxt(k,480)*y(k,137)
         mat(k,1437) = -rxt(k,486)*y(k,137)
         mat(k,2019) = mat(k,2019) + rxt(k,144)*y(k,134)
         mat(k,1511) = rxt(k,144)*y(k,133)
         mat(k,1421) = .150_r8*rxt(k,305)*y(k,207)
         mat(k,1920) = mat(k,1920) + .150_r8*rxt(k,305)*y(k,201) + .150_r8*rxt(k,355) &
                      *y(k,215)
         mat(k,1388) = .150_r8*rxt(k,355)*y(k,207)
         mat(k,328) = -(rxt(k,487)*y(k,153))
         mat(k,1427) = -rxt(k,487)*y(k,139)
         mat(k,1798) = rxt(k,227)*y(k,58)
         mat(k,2231) = rxt(k,227)*y(k,18) + 2.000_r8*rxt(k,197)*y(k,58)
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
         mat(k,360) = -(rxt(k,476)*y(k,133) + rxt(k,477)*y(k,221))
         mat(k,1997) = -rxt(k,476)*y(k,140)
         mat(k,1683) = -rxt(k,477)*y(k,140)
         mat(k,1072) = rxt(k,343)*y(k,221)
         mat(k,2267) = .100_r8*rxt(k,464)*y(k,225)
         mat(k,1664) = rxt(k,343)*y(k,92)
         mat(k,1140) = .100_r8*rxt(k,464)*y(k,123)
         mat(k,519) = -(rxt(k,315)*y(k,221))
         mat(k,1705) = -rxt(k,315)*y(k,142)
         mat(k,1550) = rxt(k,317)*y(k,201)
         mat(k,1395) = rxt(k,317)*y(k,124)
         mat(k,1544) = rxt(k,437)*y(k,192)
         mat(k,512) = rxt(k,437)*y(k,124)
         mat(k,465) = -(rxt(k,434)*y(k,124) + rxt(k,435)*y(k,137))
         mat(k,1547) = -rxt(k,434)*y(k,144)
         mat(k,1939) = -rxt(k,435)*y(k,144)
         mat(k,199) = .070_r8*rxt(k,421)*y(k,221)
         mat(k,2277) = rxt(k,419)*y(k,200)
         mat(k,180) = .060_r8*rxt(k,433)*y(k,221)
         mat(k,220) = .070_r8*rxt(k,449)*y(k,221)
         mat(k,619) = rxt(k,419)*y(k,123)
         mat(k,1697) = .070_r8*rxt(k,421)*y(k,65) + .060_r8*rxt(k,433)*y(k,145) &
                      + .070_r8*rxt(k,449)*y(k,188)
         mat(k,178) = -(rxt(k,433)*y(k,221))
         mat(k,1653) = -rxt(k,433)*y(k,145)
         mat(k,170) = .530_r8*rxt(k,410)*y(k,221)
         mat(k,1653) = mat(k,1653) + .530_r8*rxt(k,410)*y(k,6)
         mat(k,333) = -(rxt(k,436)*y(k,221))
         mat(k,1678) = -rxt(k,436)*y(k,146)
         mat(k,1839) = rxt(k,431)*y(k,222)
         mat(k,444) = rxt(k,431)*y(k,207)
         mat(k,527) = -(rxt(k,332)*y(k,221))
         mat(k,1706) = -rxt(k,332)*y(k,149)
         mat(k,1860) = rxt(k,330)*y(k,223)
         mat(k,788) = rxt(k,330)*y(k,207)
         mat(k,394) = -(rxt(k,336)*y(k,221))
         mat(k,1687) = -rxt(k,336)*y(k,150)
         mat(k,1843) = .850_r8*rxt(k,334)*y(k,224)
         mat(k,1161) = .850_r8*rxt(k,334)*y(k,207)
         mat(k,716) = -(rxt(k,478)*y(k,134) + rxt(k,480)*y(k,137) + rxt(k,483) &
                      *y(k,221))
         mat(k,1498) = -rxt(k,478)*y(k,151)
         mat(k,1943) = -rxt(k,480)*y(k,151)
         mat(k,1727) = -rxt(k,483)*y(k,151)
         mat(k,1430) = -(rxt(k,481)*y(k,18) + rxt(k,482)*y(k,58) + rxt(k,484)*y(k,124) &
                      + rxt(k,485)*y(k,134) + rxt(k,486)*y(k,137) + rxt(k,487) &
                      *y(k,139) + rxt(k,488)*y(k,221))
         mat(k,1802) = -rxt(k,481)*y(k,153)
         mat(k,2235) = -rxt(k,482)*y(k,153)
         mat(k,1568) = -rxt(k,484)*y(k,153)
         mat(k,1502) = -rxt(k,485)*y(k,153)
         mat(k,1974) = -rxt(k,486)*y(k,153)
         mat(k,330) = -rxt(k,487)*y(k,153)
         mat(k,1776) = -rxt(k,488)*y(k,153)
         mat(k,2008) = rxt(k,476)*y(k,140)
         mat(k,1502) = mat(k,1502) + rxt(k,478)*y(k,151)
         mat(k,1974) = mat(k,1974) + rxt(k,480)*y(k,151)
         mat(k,364) = rxt(k,476)*y(k,133)
         mat(k,717) = rxt(k,478)*y(k,134) + rxt(k,480)*y(k,137) + rxt(k,483)*y(k,221)
         mat(k,1776) = mat(k,1776) + rxt(k,483)*y(k,151)
         mat(k,893) = -(rxt(k,479)*y(k,221))
         mat(k,1744) = -rxt(k,479)*y(k,154)
         mat(k,1801) = rxt(k,481)*y(k,153)
         mat(k,2233) = rxt(k,482)*y(k,153)
         mat(k,289) = rxt(k,474)*y(k,125) + (rxt(k,475)+.500_r8*rxt(k,489))*y(k,221)
         mat(k,1559) = rxt(k,484)*y(k,153)
         mat(k,2083) = rxt(k,474)*y(k,66)
         mat(k,1499) = rxt(k,485)*y(k,153)
         mat(k,1947) = rxt(k,486)*y(k,153)
         mat(k,329) = rxt(k,487)*y(k,153)
         mat(k,362) = rxt(k,477)*y(k,221)
         mat(k,1429) = rxt(k,481)*y(k,18) + rxt(k,482)*y(k,58) + rxt(k,484)*y(k,124) &
                      + rxt(k,485)*y(k,134) + rxt(k,486)*y(k,137) + rxt(k,487) &
                      *y(k,139) + rxt(k,488)*y(k,221)
         mat(k,1744) = mat(k,1744) + (rxt(k,475)+.500_r8*rxt(k,489))*y(k,66) &
                      + rxt(k,477)*y(k,140) + rxt(k,488)*y(k,153)
         mat(k,264) = -(rxt(k,490)*y(k,233))
         mat(k,2349) = -rxt(k,490)*y(k,155)
         mat(k,892) = rxt(k,479)*y(k,221)
         mat(k,1669) = rxt(k,479)*y(k,154)
         mat(k,927) = .2202005_r8*rxt(k,508)*y(k,137)
         mat(k,899) = .0508005_r8*rxt(k,524)*y(k,137)
         mat(k,2255) = .1279005_r8*rxt(k,507)*y(k,194) + .0097005_r8*rxt(k,512) &
                      *y(k,196) + .0003005_r8*rxt(k,515)*y(k,211) &
                      + .1056005_r8*rxt(k,519)*y(k,212) + .0245005_r8*rxt(k,523) &
                      *y(k,218) + .0154005_r8*rxt(k,529)*y(k,228) &
                      + .0063005_r8*rxt(k,533)*y(k,231)
         mat(k,1930) = .2202005_r8*rxt(k,508)*y(k,5) + .0508005_r8*rxt(k,524)*y(k,109)
         mat(k,44) = .5931005_r8*rxt(k,526)*y(k,221)
         mat(k,50) = .1279005_r8*rxt(k,507)*y(k,123) + .2202005_r8*rxt(k,506)*y(k,207)
         mat(k,56) = .0097005_r8*rxt(k,512)*y(k,123) + .0023005_r8*rxt(k,511)*y(k,207)
         mat(k,1822) = .2202005_r8*rxt(k,506)*y(k,194) + .0023005_r8*rxt(k,511) &
                      *y(k,196) + .0031005_r8*rxt(k,514)*y(k,211) &
                      + .2381005_r8*rxt(k,518)*y(k,212) + .0508005_r8*rxt(k,522) &
                      *y(k,218) + .1364005_r8*rxt(k,528)*y(k,228) &
                      + .1677005_r8*rxt(k,532)*y(k,231)
         mat(k,62) = .0003005_r8*rxt(k,515)*y(k,123) + .0031005_r8*rxt(k,514)*y(k,207)
         mat(k,68) = .1056005_r8*rxt(k,519)*y(k,123) + .2381005_r8*rxt(k,518)*y(k,207)
         mat(k,76) = .0245005_r8*rxt(k,523)*y(k,123) + .0508005_r8*rxt(k,522)*y(k,207)
         mat(k,1632) = .5931005_r8*rxt(k,526)*y(k,176)
         mat(k,82) = .0154005_r8*rxt(k,529)*y(k,123) + .1364005_r8*rxt(k,528)*y(k,207)
         mat(k,88) = .0063005_r8*rxt(k,533)*y(k,123) + .1677005_r8*rxt(k,532)*y(k,207)
         mat(k,928) = .2067005_r8*rxt(k,508)*y(k,137)
         mat(k,900) = .1149005_r8*rxt(k,524)*y(k,137)
         mat(k,2256) = .1792005_r8*rxt(k,507)*y(k,194) + .0034005_r8*rxt(k,512) &
                      *y(k,196) + .0003005_r8*rxt(k,515)*y(k,211) &
                      + .1026005_r8*rxt(k,519)*y(k,212) + .0082005_r8*rxt(k,523) &
                      *y(k,218) + .0452005_r8*rxt(k,529)*y(k,228) &
                      + .0237005_r8*rxt(k,533)*y(k,231)
         mat(k,1931) = .2067005_r8*rxt(k,508)*y(k,5) + .1149005_r8*rxt(k,524)*y(k,109)
         mat(k,45) = .1534005_r8*rxt(k,526)*y(k,221)
         mat(k,51) = .1792005_r8*rxt(k,507)*y(k,123) + .2067005_r8*rxt(k,506)*y(k,207)
         mat(k,57) = .0034005_r8*rxt(k,512)*y(k,123) + .0008005_r8*rxt(k,511)*y(k,207)
         mat(k,1823) = .2067005_r8*rxt(k,506)*y(k,194) + .0008005_r8*rxt(k,511) &
                      *y(k,196) + .0035005_r8*rxt(k,514)*y(k,211) &
                      + .1308005_r8*rxt(k,518)*y(k,212) + .1149005_r8*rxt(k,522) &
                      *y(k,218) + .0101005_r8*rxt(k,528)*y(k,228) &
                      + .0174005_r8*rxt(k,532)*y(k,231)
         mat(k,63) = .0003005_r8*rxt(k,515)*y(k,123) + .0035005_r8*rxt(k,514)*y(k,207)
         mat(k,69) = .1026005_r8*rxt(k,519)*y(k,123) + .1308005_r8*rxt(k,518)*y(k,207)
         mat(k,77) = .0082005_r8*rxt(k,523)*y(k,123) + .1149005_r8*rxt(k,522)*y(k,207)
         mat(k,1633) = .1534005_r8*rxt(k,526)*y(k,176)
         mat(k,83) = .0452005_r8*rxt(k,529)*y(k,123) + .0101005_r8*rxt(k,528)*y(k,207)
         mat(k,89) = .0237005_r8*rxt(k,533)*y(k,123) + .0174005_r8*rxt(k,532)*y(k,207)
         mat(k,929) = .0653005_r8*rxt(k,508)*y(k,137)
         mat(k,901) = .0348005_r8*rxt(k,524)*y(k,137)
         mat(k,2257) = .0676005_r8*rxt(k,507)*y(k,194) + .1579005_r8*rxt(k,512) &
                      *y(k,196) + .0073005_r8*rxt(k,515)*y(k,211) &
                      + .0521005_r8*rxt(k,519)*y(k,212) + .0772005_r8*rxt(k,523) &
                      *y(k,218) + .0966005_r8*rxt(k,529)*y(k,228) &
                      + .0025005_r8*rxt(k,533)*y(k,231)
         mat(k,1932) = .0653005_r8*rxt(k,508)*y(k,5) + .0348005_r8*rxt(k,524)*y(k,109)
         mat(k,46) = .0459005_r8*rxt(k,526)*y(k,221)
         mat(k,52) = .0676005_r8*rxt(k,507)*y(k,123) + .0653005_r8*rxt(k,506)*y(k,207)
         mat(k,58) = .1579005_r8*rxt(k,512)*y(k,123) + .0843005_r8*rxt(k,511)*y(k,207)
         mat(k,1824) = .0653005_r8*rxt(k,506)*y(k,194) + .0843005_r8*rxt(k,511) &
                      *y(k,196) + .0003005_r8*rxt(k,514)*y(k,211) &
                      + .0348005_r8*rxt(k,518)*y(k,212) + .0348005_r8*rxt(k,522) &
                      *y(k,218) + .0763005_r8*rxt(k,528)*y(k,228) + .086_r8*rxt(k,532) &
                      *y(k,231)
         mat(k,64) = .0073005_r8*rxt(k,515)*y(k,123) + .0003005_r8*rxt(k,514)*y(k,207)
         mat(k,70) = .0521005_r8*rxt(k,519)*y(k,123) + .0348005_r8*rxt(k,518)*y(k,207)
         mat(k,78) = .0772005_r8*rxt(k,523)*y(k,123) + .0348005_r8*rxt(k,522)*y(k,207)
         mat(k,1634) = .0459005_r8*rxt(k,526)*y(k,176)
         mat(k,84) = .0966005_r8*rxt(k,529)*y(k,123) + .0763005_r8*rxt(k,528)*y(k,207)
         mat(k,90) = .0025005_r8*rxt(k,533)*y(k,123) + .086_r8*rxt(k,532)*y(k,207)
         mat(k,930) = .1749305_r8*rxt(k,505)*y(k,125) + .1284005_r8*rxt(k,508) &
                      *y(k,137)
         mat(k,831) = .0590245_r8*rxt(k,513)*y(k,125) + .0033005_r8*rxt(k,516) &
                      *y(k,137)
         mat(k,902) = .1749305_r8*rxt(k,521)*y(k,125) + .0554005_r8*rxt(k,524) &
                      *y(k,137)
         mat(k,2258) = .079_r8*rxt(k,507)*y(k,194) + .0059005_r8*rxt(k,512)*y(k,196) &
                      + .0057005_r8*rxt(k,515)*y(k,211) + .0143005_r8*rxt(k,519) &
                      *y(k,212) + .0332005_r8*rxt(k,523)*y(k,218) &
                      + .0073005_r8*rxt(k,529)*y(k,228) + .011_r8*rxt(k,533)*y(k,231)
         mat(k,2073) = .1749305_r8*rxt(k,505)*y(k,5) + .0590245_r8*rxt(k,513)*y(k,97) &
                      + .1749305_r8*rxt(k,521)*y(k,109)
         mat(k,1933) = .1284005_r8*rxt(k,508)*y(k,5) + .0033005_r8*rxt(k,516)*y(k,97) &
                      + .0554005_r8*rxt(k,524)*y(k,109)
         mat(k,47) = .0085005_r8*rxt(k,526)*y(k,221)
         mat(k,53) = .079_r8*rxt(k,507)*y(k,123) + .1284005_r8*rxt(k,506)*y(k,207)
         mat(k,59) = .0059005_r8*rxt(k,512)*y(k,123) + .0443005_r8*rxt(k,511)*y(k,207)
         mat(k,1825) = .1284005_r8*rxt(k,506)*y(k,194) + .0443005_r8*rxt(k,511) &
                      *y(k,196) + .0271005_r8*rxt(k,514)*y(k,211) &
                      + .0076005_r8*rxt(k,518)*y(k,212) + .0554005_r8*rxt(k,522) &
                      *y(k,218) + .2157005_r8*rxt(k,528)*y(k,228) &
                      + .0512005_r8*rxt(k,532)*y(k,231)
         mat(k,65) = .0057005_r8*rxt(k,515)*y(k,123) + .0271005_r8*rxt(k,514)*y(k,207)
         mat(k,71) = .0143005_r8*rxt(k,519)*y(k,123) + .0076005_r8*rxt(k,518)*y(k,207)
         mat(k,79) = .0332005_r8*rxt(k,523)*y(k,123) + .0554005_r8*rxt(k,522)*y(k,207)
         mat(k,1635) = .0085005_r8*rxt(k,526)*y(k,176)
         mat(k,85) = .0073005_r8*rxt(k,529)*y(k,123) + .2157005_r8*rxt(k,528)*y(k,207)
         mat(k,91) = .011_r8*rxt(k,533)*y(k,123) + .0512005_r8*rxt(k,532)*y(k,207)
         mat(k,931) = .5901905_r8*rxt(k,505)*y(k,125) + .114_r8*rxt(k,508)*y(k,137)
         mat(k,832) = .0250245_r8*rxt(k,513)*y(k,125)
         mat(k,903) = .5901905_r8*rxt(k,521)*y(k,125) + .1278005_r8*rxt(k,524) &
                      *y(k,137)
         mat(k,2259) = .1254005_r8*rxt(k,507)*y(k,194) + .0536005_r8*rxt(k,512) &
                      *y(k,196) + .0623005_r8*rxt(k,515)*y(k,211) &
                      + .0166005_r8*rxt(k,519)*y(k,212) + .130_r8*rxt(k,523)*y(k,218) &
                      + .238_r8*rxt(k,529)*y(k,228) + .1185005_r8*rxt(k,533)*y(k,231)
         mat(k,2074) = .5901905_r8*rxt(k,505)*y(k,5) + .0250245_r8*rxt(k,513)*y(k,97) &
                      + .5901905_r8*rxt(k,521)*y(k,109)
         mat(k,1934) = .114_r8*rxt(k,508)*y(k,5) + .1278005_r8*rxt(k,524)*y(k,109)
         mat(k,48) = .0128005_r8*rxt(k,526)*y(k,221)
         mat(k,54) = .1254005_r8*rxt(k,507)*y(k,123) + .114_r8*rxt(k,506)*y(k,207)
         mat(k,60) = .0536005_r8*rxt(k,512)*y(k,123) + .1621005_r8*rxt(k,511)*y(k,207)
         mat(k,1826) = .114_r8*rxt(k,506)*y(k,194) + .1621005_r8*rxt(k,511)*y(k,196) &
                      + .0474005_r8*rxt(k,514)*y(k,211) + .0113005_r8*rxt(k,518) &
                      *y(k,212) + .1278005_r8*rxt(k,522)*y(k,218) &
                      + .0738005_r8*rxt(k,528)*y(k,228) + .1598005_r8*rxt(k,532) &
                      *y(k,231)
         mat(k,66) = .0623005_r8*rxt(k,515)*y(k,123) + .0474005_r8*rxt(k,514)*y(k,207)
         mat(k,72) = .0166005_r8*rxt(k,519)*y(k,123) + .0113005_r8*rxt(k,518)*y(k,207)
         mat(k,80) = .130_r8*rxt(k,523)*y(k,123) + .1278005_r8*rxt(k,522)*y(k,207)
         mat(k,1636) = .0128005_r8*rxt(k,526)*y(k,176)
         mat(k,86) = .238_r8*rxt(k,529)*y(k,123) + .0738005_r8*rxt(k,528)*y(k,207)
         mat(k,92) = .1185005_r8*rxt(k,533)*y(k,123) + .1598005_r8*rxt(k,532)*y(k,207)
         mat(k,49) = -(rxt(k,526)*y(k,221))
         mat(k,1637) = -rxt(k,526)*y(k,176)
         mat(k,192) = .100_r8*rxt(k,441)*y(k,221)
         mat(k,210) = .230_r8*rxt(k,443)*y(k,221)
         mat(k,1657) = .100_r8*rxt(k,441)*y(k,184) + .230_r8*rxt(k,443)*y(k,186)
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
         mat(k,637) = -(rxt(k,465)*y(k,221))
         mat(k,1719) = -rxt(k,465)*y(k,178)
         mat(k,1865) = rxt(k,463)*y(k,225)
         mat(k,1141) = rxt(k,463)*y(k,207)
         mat(k,612) = -(rxt(k,466)*y(k,221))
         mat(k,1716) = -rxt(k,466)*y(k,179)
         mat(k,2286) = .200_r8*rxt(k,459)*y(k,219) + .200_r8*rxt(k,469)*y(k,226)
         mat(k,2138) = .500_r8*rxt(k,457)*y(k,219)
         mat(k,1057) = .200_r8*rxt(k,459)*y(k,123) + .500_r8*rxt(k,457)*y(k,202)
         mat(k,1118) = .200_r8*rxt(k,469)*y(k,123)
         mat(k,483) = -(rxt(k,470)*y(k,221))
         mat(k,1700) = -rxt(k,470)*y(k,180)
         mat(k,1856) = rxt(k,468)*y(k,226)
         mat(k,1117) = rxt(k,468)*y(k,207)
         mat(k,995) = -(rxt(k,471)*y(k,125) + rxt(k,472)*y(k,221))
         mat(k,2089) = -rxt(k,471)*y(k,181)
         mat(k,1751) = -rxt(k,472)*y(k,181)
         mat(k,940) = .330_r8*rxt(k,452)*y(k,137)
         mat(k,912) = .330_r8*rxt(k,455)*y(k,137)
         mat(k,2305) = .800_r8*rxt(k,459)*y(k,219) + .800_r8*rxt(k,469)*y(k,226)
         mat(k,2089) = mat(k,2089) + rxt(k,460)*y(k,219)
         mat(k,1953) = .330_r8*rxt(k,452)*y(k,5) + .330_r8*rxt(k,455)*y(k,109)
         mat(k,613) = rxt(k,466)*y(k,221)
         mat(k,2146) = .500_r8*rxt(k,457)*y(k,219) + rxt(k,467)*y(k,226)
         mat(k,1059) = .800_r8*rxt(k,459)*y(k,123) + rxt(k,460)*y(k,125) &
                      + .500_r8*rxt(k,457)*y(k,202)
         mat(k,1751) = mat(k,1751) + rxt(k,466)*y(k,179)
         mat(k,1121) = .800_r8*rxt(k,469)*y(k,123) + rxt(k,467)*y(k,202)
         mat(k,1097) = -(rxt(k,473)*y(k,221))
         mat(k,1759) = -rxt(k,473)*y(k,182)
         mat(k,943) = .300_r8*rxt(k,452)*y(k,137)
         mat(k,915) = .300_r8*rxt(k,455)*y(k,137)
         mat(k,2312) = .900_r8*rxt(k,464)*y(k,225)
         mat(k,1959) = .300_r8*rxt(k,452)*y(k,5) + .300_r8*rxt(k,455)*y(k,109)
         mat(k,2152) = rxt(k,462)*y(k,225)
         mat(k,1145) = .900_r8*rxt(k,464)*y(k,123) + rxt(k,462)*y(k,202)
         mat(k,650) = -(rxt(k,440)*y(k,221))
         mat(k,1720) = -rxt(k,440)*y(k,183)
         mat(k,1866) = rxt(k,438)*y(k,227)
         mat(k,732) = rxt(k,438)*y(k,207)
         mat(k,190) = -(rxt(k,441)*y(k,221))
         mat(k,1655) = -rxt(k,441)*y(k,184)
         mat(k,206) = -(rxt(k,407)*y(k,221))
         mat(k,1658) = -rxt(k,407)*y(k,185)
         mat(k,1835) = rxt(k,404)*y(k,229)
         mat(k,1180) = rxt(k,404)*y(k,207)
         mat(k,211) = -(rxt(k,443)*y(k,221))
         mat(k,1659) = -rxt(k,443)*y(k,186)
         mat(k,705) = -(rxt(k,446)*y(k,221))
         mat(k,1726) = -rxt(k,446)*y(k,187)
         mat(k,1870) = rxt(k,444)*y(k,230)
         mat(k,748) = rxt(k,444)*y(k,207)
         mat(k,219) = -(rxt(k,449)*y(k,221))
         mat(k,1660) = -rxt(k,449)*y(k,188)
         mat(k,212) = .150_r8*rxt(k,443)*y(k,221)
         mat(k,1660) = mat(k,1660) + .150_r8*rxt(k,443)*y(k,186)
         mat(k,418) = -(rxt(k,450)*y(k,221))
         mat(k,1691) = -rxt(k,450)*y(k,189)
         mat(k,1847) = rxt(k,447)*y(k,232)
         mat(k,499) = rxt(k,447)*y(k,207)
         mat(k,513) = -(rxt(k,408)*y(k,207) + rxt(k,409)*y(k,123) + rxt(k,437) &
                      *y(k,124))
         mat(k,1859) = -rxt(k,408)*y(k,192)
         mat(k,2281) = -rxt(k,409)*y(k,192)
         mat(k,1549) = -rxt(k,437)*y(k,192)
         mat(k,242) = rxt(k,414)*y(k,221)
         mat(k,1704) = rxt(k,414)*y(k,21)
         mat(k,960) = -(rxt(k,369)*y(k,207) + (rxt(k,370) + rxt(k,371)) * y(k,123))
         mat(k,1885) = -rxt(k,369)*y(k,193)
         mat(k,2302) = -(rxt(k,370) + rxt(k,371)) * y(k,193)
         mat(k,630) = rxt(k,372)*y(k,221)
         mat(k,230) = rxt(k,373)*y(k,221)
         mat(k,1747) = rxt(k,372)*y(k,2) + rxt(k,373)*y(k,14)
         mat(k,55) = -(rxt(k,506)*y(k,207) + rxt(k,507)*y(k,123))
         mat(k,1827) = -rxt(k,506)*y(k,194)
         mat(k,2260) = -rxt(k,507)*y(k,194)
         mat(k,932) = rxt(k,509)*y(k,221)
         mat(k,1638) = rxt(k,509)*y(k,5)
         mat(k,492) = -(rxt(k,411)*y(k,207) + rxt(k,412)*y(k,123))
         mat(k,1857) = -rxt(k,411)*y(k,195)
         mat(k,2278) = -rxt(k,412)*y(k,195)
         mat(k,171) = .350_r8*rxt(k,410)*y(k,221)
         mat(k,414) = rxt(k,413)*y(k,221)
         mat(k,1701) = .350_r8*rxt(k,410)*y(k,6) + rxt(k,413)*y(k,7)
         mat(k,61) = -(rxt(k,511)*y(k,207) + rxt(k,512)*y(k,123))
         mat(k,1828) = -rxt(k,511)*y(k,196)
         mat(k,2261) = -rxt(k,512)*y(k,196)
         mat(k,167) = rxt(k,510)*y(k,221)
         mat(k,1639) = rxt(k,510)*y(k,6)
         mat(k,432) = -(rxt(k,415)*y(k,207) + rxt(k,417)*y(k,123))
         mat(k,1848) = -rxt(k,415)*y(k,197)
         mat(k,2272) = -rxt(k,417)*y(k,197)
         mat(k,340) = rxt(k,416)*y(k,221)
         mat(k,193) = .070_r8*rxt(k,441)*y(k,221)
         mat(k,213) = .060_r8*rxt(k,443)*y(k,221)
         mat(k,1693) = rxt(k,416)*y(k,22) + .070_r8*rxt(k,441)*y(k,184) &
                      + .060_r8*rxt(k,443)*y(k,186)
         mat(k,859) = -(4._r8*rxt(k,292)*y(k,198) + rxt(k,293)*y(k,202) + rxt(k,294) &
                      *y(k,207) + rxt(k,295)*y(k,123))
         mat(k,2143) = -rxt(k,293)*y(k,198)
         mat(k,1882) = -rxt(k,294)*y(k,198)
         mat(k,2299) = -rxt(k,295)*y(k,198)
         mat(k,348) = .500_r8*rxt(k,297)*y(k,221)
         mat(k,301) = rxt(k,298)*y(k,55) + rxt(k,299)*y(k,221)
         mat(k,2198) = rxt(k,298)*y(k,27)
         mat(k,1740) = .500_r8*rxt(k,297)*y(k,26) + rxt(k,299)*y(k,27)
         mat(k,804) = -(rxt(k,321)*y(k,202) + rxt(k,322)*y(k,207) + rxt(k,323) &
                      *y(k,123))
         mat(k,2141) = -rxt(k,321)*y(k,199)
         mat(k,1878) = -rxt(k,322)*y(k,199)
         mat(k,2297) = -rxt(k,323)*y(k,199)
         mat(k,407) = rxt(k,324)*y(k,221)
         mat(k,111) = rxt(k,325)*y(k,221)
         mat(k,1735) = rxt(k,324)*y(k,29) + rxt(k,325)*y(k,30)
         mat(k,620) = -(rxt(k,418)*y(k,207) + rxt(k,419)*y(k,123))
         mat(k,1863) = -rxt(k,418)*y(k,200)
         mat(k,2287) = -rxt(k,419)*y(k,200)
         mat(k,274) = rxt(k,420)*y(k,221)
         mat(k,2287) = mat(k,2287) + rxt(k,409)*y(k,192)
         mat(k,1941) = rxt(k,435)*y(k,144)
         mat(k,466) = rxt(k,435)*y(k,137)
         mat(k,514) = rxt(k,409)*y(k,123) + .400_r8*rxt(k,408)*y(k,207)
         mat(k,1863) = mat(k,1863) + .400_r8*rxt(k,408)*y(k,192)
         mat(k,1717) = rxt(k,420)*y(k,31)
         mat(k,1413) = -(4._r8*rxt(k,303)*y(k,201) + rxt(k,304)*y(k,202) + rxt(k,305) &
                      *y(k,207) + rxt(k,306)*y(k,123) + rxt(k,317)*y(k,124) + rxt(k,344) &
                      *y(k,213) + rxt(k,377)*y(k,209) + rxt(k,382)*y(k,210) + rxt(k,391) &
                      *y(k,100) + rxt(k,402)*y(k,229))
         mat(k,2167) = -rxt(k,304)*y(k,201)
         mat(k,1909) = -rxt(k,305)*y(k,201)
         mat(k,2327) = -rxt(k,306)*y(k,201)
         mat(k,1567) = -rxt(k,317)*y(k,201)
         mat(k,1337) = -rxt(k,344)*y(k,201)
         mat(k,1279) = -rxt(k,377)*y(k,201)
         mat(k,1314) = -rxt(k,382)*y(k,201)
         mat(k,1213) = -rxt(k,391)*y(k,201)
         mat(k,1189) = -rxt(k,402)*y(k,201)
         mat(k,948) = .060_r8*rxt(k,452)*y(k,137)
         mat(k,1049) = rxt(k,300)*y(k,125) + rxt(k,301)*y(k,221)
         mat(k,1241) = rxt(k,326)*y(k,125) + rxt(k,327)*y(k,221)
         mat(k,607) = .500_r8*rxt(k,308)*y(k,221)
         mat(k,843) = .080_r8*rxt(k,397)*y(k,137)
         mat(k,1231) = .100_r8*rxt(k,350)*y(k,137)
         mat(k,920) = .060_r8*rxt(k,455)*y(k,137)
         mat(k,1359) = .280_r8*rxt(k,364)*y(k,137)
         mat(k,2327) = mat(k,2327) + .530_r8*rxt(k,348)*y(k,213) + rxt(k,357)*y(k,215) &
                      + rxt(k,360)*y(k,217) + rxt(k,335)*y(k,224)
         mat(k,2112) = rxt(k,300)*y(k,44) + rxt(k,326)*y(k,48) + .530_r8*rxt(k,347) &
                      *y(k,213) + rxt(k,358)*y(k,215)
         mat(k,1973) = .060_r8*rxt(k,452)*y(k,5) + .080_r8*rxt(k,397)*y(k,97) &
                      + .100_r8*rxt(k,350)*y(k,104) + .060_r8*rxt(k,455)*y(k,109) &
                      + .280_r8*rxt(k,364)*y(k,110)
         mat(k,1100) = .650_r8*rxt(k,473)*y(k,221)
         mat(k,1413) = mat(k,1413) + .530_r8*rxt(k,344)*y(k,213)
         mat(k,2167) = mat(k,2167) + .260_r8*rxt(k,345)*y(k,213) + rxt(k,354)*y(k,215) &
                      + .300_r8*rxt(k,333)*y(k,224)
         mat(k,1909) = mat(k,1909) + .450_r8*rxt(k,355)*y(k,215) + .200_r8*rxt(k,359) &
                      *y(k,217) + .150_r8*rxt(k,334)*y(k,224)
         mat(k,1337) = mat(k,1337) + .530_r8*rxt(k,348)*y(k,123) + .530_r8*rxt(k,347) &
                      *y(k,125) + .530_r8*rxt(k,344)*y(k,201) + .260_r8*rxt(k,345) &
                      *y(k,202)
         mat(k,1381) = rxt(k,357)*y(k,123) + rxt(k,358)*y(k,125) + rxt(k,354)*y(k,202) &
                      + .450_r8*rxt(k,355)*y(k,207) + 4.000_r8*rxt(k,356)*y(k,215)
         mat(k,674) = rxt(k,360)*y(k,123) + .200_r8*rxt(k,359)*y(k,207)
         mat(k,1775) = rxt(k,301)*y(k,44) + rxt(k,327)*y(k,48) + .500_r8*rxt(k,308) &
                      *y(k,50) + .650_r8*rxt(k,473)*y(k,182)
         mat(k,1166) = rxt(k,335)*y(k,123) + .300_r8*rxt(k,333)*y(k,202) &
                      + .150_r8*rxt(k,334)*y(k,207)
         mat(k,2182) = -(rxt(k,194)*y(k,58) + (4._r8*rxt(k,271) + 4._r8*rxt(k,272) &
                      ) * y(k,202) + rxt(k,273)*y(k,207) + rxt(k,274)*y(k,123) &
                      + rxt(k,293)*y(k,198) + rxt(k,304)*y(k,201) + rxt(k,321) &
                      *y(k,199) + rxt(k,333)*y(k,224) + rxt(k,345)*y(k,213) + rxt(k,354) &
                      *y(k,215) + rxt(k,378)*y(k,209) + rxt(k,383)*y(k,210) + rxt(k,392) &
                      *y(k,100) + rxt(k,403)*y(k,229) + rxt(k,457)*y(k,219) + rxt(k,462) &
                      *y(k,225) + rxt(k,467)*y(k,226))
         mat(k,2250) = -rxt(k,194)*y(k,202)
         mat(k,1925) = -rxt(k,273)*y(k,202)
         mat(k,2343) = -rxt(k,274)*y(k,202)
         mat(k,866) = -rxt(k,293)*y(k,202)
         mat(k,1424) = -rxt(k,304)*y(k,202)
         mat(k,812) = -rxt(k,321)*y(k,202)
         mat(k,1171) = -rxt(k,333)*y(k,202)
         mat(k,1346) = -rxt(k,345)*y(k,202)
         mat(k,1391) = -rxt(k,354)*y(k,202)
         mat(k,1289) = -rxt(k,378)*y(k,202)
         mat(k,1324) = -rxt(k,383)*y(k,202)
         mat(k,1222) = -rxt(k,392)*y(k,202)
         mat(k,1197) = -rxt(k,403)*y(k,202)
         mat(k,1070) = -rxt(k,457)*y(k,202)
         mat(k,1158) = -rxt(k,462)*y(k,202)
         mat(k,1137) = -rxt(k,467)*y(k,202)
         mat(k,1038) = .280_r8*rxt(k,320)*y(k,137)
         mat(k,691) = rxt(k,307)*y(k,221)
         mat(k,455) = .700_r8*rxt(k,276)*y(k,221)
         mat(k,1473) = rxt(k,188)*y(k,55) + rxt(k,244)*y(k,72) + rxt(k,283)*y(k,220) &
                      + rxt(k,277)*y(k,221)
         mat(k,2223) = rxt(k,188)*y(k,53)
         mat(k,890) = rxt(k,244)*y(k,53)
         mat(k,849) = .050_r8*rxt(k,397)*y(k,137)
         mat(k,1222) = mat(k,1222) + rxt(k,391)*y(k,201)
         mat(k,2343) = mat(k,2343) + rxt(k,306)*y(k,201) + .830_r8*rxt(k,423)*y(k,203) &
                      + .170_r8*rxt(k,429)*y(k,216)
         mat(k,1990) = .280_r8*rxt(k,320)*y(k,28) + .050_r8*rxt(k,397)*y(k,97)
         mat(k,1424) = mat(k,1424) + rxt(k,391)*y(k,100) + rxt(k,306)*y(k,123) &
                      + 4.000_r8*rxt(k,303)*y(k,201) + .900_r8*rxt(k,304)*y(k,202) &
                      + .490_r8*rxt(k,305)*y(k,207) + rxt(k,377)*y(k,209) + rxt(k,382) &
                      *y(k,210) + rxt(k,344)*y(k,213) + rxt(k,353)*y(k,215) &
                      + rxt(k,402)*y(k,229)
         mat(k,2182) = mat(k,2182) + .900_r8*rxt(k,304)*y(k,201)
         mat(k,767) = .830_r8*rxt(k,423)*y(k,123) + .330_r8*rxt(k,422)*y(k,207)
         mat(k,1925) = mat(k,1925) + .490_r8*rxt(k,305)*y(k,201) + .330_r8*rxt(k,422) &
                      *y(k,203) + .070_r8*rxt(k,428)*y(k,216)
         mat(k,1289) = mat(k,1289) + rxt(k,377)*y(k,201)
         mat(k,1324) = mat(k,1324) + rxt(k,382)*y(k,201)
         mat(k,1346) = mat(k,1346) + rxt(k,344)*y(k,201)
         mat(k,1391) = mat(k,1391) + rxt(k,353)*y(k,201)
         mat(k,880) = .170_r8*rxt(k,429)*y(k,123) + .070_r8*rxt(k,428)*y(k,207)
         mat(k,1627) = rxt(k,283)*y(k,53)
         mat(k,1792) = rxt(k,307)*y(k,49) + .700_r8*rxt(k,276)*y(k,52) + rxt(k,277) &
                      *y(k,53)
         mat(k,1197) = mat(k,1197) + rxt(k,402)*y(k,201)
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
         mat(k,761) = -(rxt(k,422)*y(k,207) + rxt(k,423)*y(k,123) + rxt(k,424) &
                      *y(k,124))
         mat(k,1874) = -rxt(k,422)*y(k,203)
         mat(k,2294) = -rxt(k,423)*y(k,203)
         mat(k,1556) = -rxt(k,424)*y(k,203)
         mat(k,551) = -((rxt(k,341) + rxt(k,342)) * y(k,123))
         mat(k,2282) = -(rxt(k,341) + rxt(k,342)) * y(k,204)
         mat(k,353) = rxt(k,340)*y(k,221)
         mat(k,1708) = rxt(k,340)*y(k,15)
         mat(k,457) = -(rxt(k,312)*y(k,134))
         mat(k,1496) = -rxt(k,312)*y(k,205)
         mat(k,2275) = .750_r8*rxt(k,310)*y(k,206)
         mat(k,779) = .750_r8*rxt(k,310)*y(k,123)
         mat(k,780) = -(rxt(k,309)*y(k,207) + rxt(k,310)*y(k,123))
         mat(k,1876) = -rxt(k,309)*y(k,206)
         mat(k,2295) = -rxt(k,310)*y(k,206)
         mat(k,536) = rxt(k,316)*y(k,221)
         mat(k,1732) = rxt(k,316)*y(k,24)
         mat(k,1919) = -((rxt(k,147) + rxt(k,148) + rxt(k,149)) * y(k,75) + rxt(k,151) &
                      *y(k,133) + rxt(k,152)*y(k,137) + rxt(k,156)*y(k,221) &
                      + 4._r8*rxt(k,161)*y(k,207) + rxt(k,171)*y(k,125) + rxt(k,176) &
                      *y(k,123) + rxt(k,181)*y(k,124) + (rxt(k,191) + rxt(k,192) &
                      ) * y(k,55) + rxt(k,198)*y(k,58) + rxt(k,224)*y(k,16) + rxt(k,230) &
                      *y(k,18) + rxt(k,267)*y(k,41) + rxt(k,273)*y(k,202) + rxt(k,280) &
                      *y(k,208) + rxt(k,294)*y(k,198) + rxt(k,305)*y(k,201) + rxt(k,309) &
                      *y(k,206) + rxt(k,322)*y(k,199) + rxt(k,330)*y(k,223) + rxt(k,334) &
                      *y(k,224) + rxt(k,346)*y(k,213) + rxt(k,355)*y(k,215) + rxt(k,359) &
                      *y(k,217) + rxt(k,369)*y(k,193) + rxt(k,379)*y(k,209) + rxt(k,384) &
                      *y(k,210) + rxt(k,393)*y(k,100) + rxt(k,404)*y(k,229) + rxt(k,408) &
                      *y(k,192) + rxt(k,411)*y(k,195) + rxt(k,415)*y(k,197) + rxt(k,418) &
                      *y(k,200) + rxt(k,422)*y(k,203) + rxt(k,425)*y(k,214) + rxt(k,428) &
                      *y(k,216) + rxt(k,431)*y(k,222) + rxt(k,438)*y(k,227) + rxt(k,444) &
                      *y(k,230) + rxt(k,447)*y(k,232) + rxt(k,458)*y(k,219) + rxt(k,463) &
                      *y(k,225) + rxt(k,468)*y(k,226))
         mat(k,2038) = -(rxt(k,147) + rxt(k,148) + rxt(k,149)) * y(k,207)
         mat(k,2018) = -rxt(k,151)*y(k,207)
         mat(k,1984) = -rxt(k,152)*y(k,207)
         mat(k,1786) = -rxt(k,156)*y(k,207)
         mat(k,2123) = -rxt(k,171)*y(k,207)
         mat(k,2337) = -rxt(k,176)*y(k,207)
         mat(k,1578) = -rxt(k,181)*y(k,207)
         mat(k,2217) = -(rxt(k,191) + rxt(k,192)) * y(k,207)
         mat(k,2244) = -rxt(k,198)*y(k,207)
         mat(k,1452) = -rxt(k,224)*y(k,207)
         mat(k,1811) = -rxt(k,230)*y(k,207)
         mat(k,1533) = -rxt(k,267)*y(k,207)
         mat(k,2176) = -rxt(k,273)*y(k,207)
         mat(k,442) = -rxt(k,280)*y(k,207)
         mat(k,865) = -rxt(k,294)*y(k,207)
         mat(k,1420) = -rxt(k,305)*y(k,207)
         mat(k,786) = -rxt(k,309)*y(k,207)
         mat(k,811) = -rxt(k,322)*y(k,207)
         mat(k,796) = -rxt(k,330)*y(k,207)
         mat(k,1170) = -rxt(k,334)*y(k,207)
         mat(k,1343) = -rxt(k,346)*y(k,207)
         mat(k,1387) = -rxt(k,355)*y(k,207)
         mat(k,677) = -rxt(k,359)*y(k,207)
         mat(k,968) = -rxt(k,369)*y(k,207)
         mat(k,1285) = -rxt(k,379)*y(k,207)
         mat(k,1320) = -rxt(k,384)*y(k,207)
         mat(k,1219) = -rxt(k,393)*y(k,207)
         mat(k,1194) = -rxt(k,404)*y(k,207)
         mat(k,517) = -rxt(k,408)*y(k,207)
         mat(k,497) = -rxt(k,411)*y(k,207)
         mat(k,436) = -rxt(k,415)*y(k,207)
         mat(k,623) = -rxt(k,418)*y(k,207)
         mat(k,766) = -rxt(k,422)*y(k,207)
         mat(k,727) = -rxt(k,425)*y(k,207)
         mat(k,879) = -rxt(k,428)*y(k,207)
         mat(k,449) = -rxt(k,431)*y(k,207)
         mat(k,742) = -rxt(k,438)*y(k,207)
         mat(k,759) = -rxt(k,444)*y(k,207)
         mat(k,505) = -rxt(k,447)*y(k,207)
         mat(k,1068) = -rxt(k,458)*y(k,207)
         mat(k,1156) = -rxt(k,463)*y(k,207)
         mat(k,1134) = -rxt(k,468)*y(k,207)
         mat(k,951) = .570_r8*rxt(k,452)*y(k,137)
         mat(k,173) = .650_r8*rxt(k,410)*y(k,221)
         mat(k,1452) = mat(k,1452) + rxt(k,223)*y(k,41)
         mat(k,1811) = mat(k,1811) + rxt(k,235)*y(k,221)
         mat(k,299) = .350_r8*rxt(k,289)*y(k,221)
         mat(k,541) = .130_r8*rxt(k,291)*y(k,137)
         mat(k,271) = rxt(k,296)*y(k,221)
         mat(k,1035) = .280_r8*rxt(k,320)*y(k,137)
         mat(k,1533) = mat(k,1533) + rxt(k,223)*y(k,16) + rxt(k,187)*y(k,55) &
                      + rxt(k,268)*y(k,125) + rxt(k,269)*y(k,133)
         mat(k,593) = rxt(k,252)*y(k,55) + rxt(k,253)*y(k,221)
         mat(k,371) = rxt(k,255)*y(k,55) + rxt(k,256)*y(k,221)
         mat(k,105) = rxt(k,302)*y(k,221)
         mat(k,802) = rxt(k,275)*y(k,221)
         mat(k,1468) = rxt(k,284)*y(k,220)
         mat(k,2217) = mat(k,2217) + rxt(k,187)*y(k,41) + rxt(k,252)*y(k,42) &
                      + rxt(k,255)*y(k,45) + rxt(k,190)*y(k,78)
         mat(k,2244) = mat(k,2244) + rxt(k,194)*y(k,202) + rxt(k,205)*y(k,221)
         mat(k,1178) = rxt(k,287)*y(k,221)
         mat(k,201) = .730_r8*rxt(k,421)*y(k,221)
         mat(k,292) = .500_r8*rxt(k,489)*y(k,221)
         mat(k,1115) = rxt(k,313)*y(k,221)
         mat(k,982) = rxt(k,314)*y(k,221)
         mat(k,2038) = mat(k,2038) + rxt(k,150)*y(k,134)
         mat(k,600) = rxt(k,190)*y(k,55) + rxt(k,146)*y(k,133) + rxt(k,155)*y(k,221)
         mat(k,188) = rxt(k,278)*y(k,221)
         mat(k,871) = rxt(k,279)*y(k,221)
         mat(k,1085) = rxt(k,343)*y(k,221)
         mat(k,1094) = rxt(k,328)*y(k,221)
         mat(k,846) = .370_r8*rxt(k,397)*y(k,137)
         mat(k,584) = .300_r8*rxt(k,388)*y(k,221)
         mat(k,566) = rxt(k,389)*y(k,221)
         mat(k,1219) = mat(k,1219) + rxt(k,394)*y(k,123) + rxt(k,395)*y(k,125) &
                      + rxt(k,391)*y(k,201) + 1.200_r8*rxt(k,392)*y(k,202)
         mat(k,405) = rxt(k,396)*y(k,221)
         mat(k,1235) = .140_r8*rxt(k,350)*y(k,137)
         mat(k,317) = .200_r8*rxt(k,352)*y(k,221)
         mat(k,574) = .500_r8*rxt(k,363)*y(k,221)
         mat(k,923) = .570_r8*rxt(k,455)*y(k,137)
         mat(k,1366) = .280_r8*rxt(k,364)*y(k,137)
         mat(k,429) = rxt(k,400)*y(k,221)
         mat(k,1016) = rxt(k,401)*y(k,221)
         mat(k,2337) = mat(k,2337) + rxt(k,394)*y(k,100) + rxt(k,370)*y(k,193) &
                      + rxt(k,412)*y(k,195) + rxt(k,417)*y(k,197) + rxt(k,295) &
                      *y(k,198) + rxt(k,323)*y(k,199) + rxt(k,274)*y(k,202) &
                      + .170_r8*rxt(k,423)*y(k,203) + rxt(k,341)*y(k,204) &
                      + .250_r8*rxt(k,310)*y(k,206) + rxt(k,282)*y(k,208) &
                      + .920_r8*rxt(k,380)*y(k,209) + .920_r8*rxt(k,386)*y(k,210) &
                      + .470_r8*rxt(k,348)*y(k,213) + .400_r8*rxt(k,426)*y(k,214) &
                      + .830_r8*rxt(k,429)*y(k,216) + rxt(k,432)*y(k,222) + rxt(k,331) &
                      *y(k,223) + .900_r8*rxt(k,464)*y(k,225) + .800_r8*rxt(k,469) &
                      *y(k,226) + rxt(k,439)*y(k,227) + rxt(k,405)*y(k,229) &
                      + rxt(k,445)*y(k,230) + rxt(k,448)*y(k,232)
         mat(k,2123) = mat(k,2123) + rxt(k,268)*y(k,41) + rxt(k,395)*y(k,100) &
                      + rxt(k,381)*y(k,209) + rxt(k,387)*y(k,210) + .470_r8*rxt(k,347) &
                      *y(k,213) + rxt(k,174)*y(k,221) + rxt(k,406)*y(k,229)
         mat(k,2018) = mat(k,2018) + rxt(k,269)*y(k,41) + rxt(k,146)*y(k,78)
         mat(k,1510) = rxt(k,150)*y(k,75) + rxt(k,312)*y(k,205)
         mat(k,1984) = mat(k,1984) + .570_r8*rxt(k,452)*y(k,5) + .130_r8*rxt(k,291) &
                      *y(k,24) + .280_r8*rxt(k,320)*y(k,28) + .370_r8*rxt(k,397) &
                      *y(k,97) + .140_r8*rxt(k,350)*y(k,104) + .570_r8*rxt(k,455) &
                      *y(k,109) + .280_r8*rxt(k,364)*y(k,110) + rxt(k,158)*y(k,221)
         mat(k,182) = .800_r8*rxt(k,433)*y(k,221)
         mat(k,896) = rxt(k,479)*y(k,221)
         mat(k,1104) = .200_r8*rxt(k,473)*y(k,221)
         mat(k,196) = .280_r8*rxt(k,441)*y(k,221)
         mat(k,218) = .380_r8*rxt(k,443)*y(k,221)
         mat(k,223) = .630_r8*rxt(k,449)*y(k,221)
         mat(k,968) = mat(k,968) + rxt(k,370)*y(k,123)
         mat(k,497) = mat(k,497) + rxt(k,412)*y(k,123)
         mat(k,436) = mat(k,436) + rxt(k,417)*y(k,123)
         mat(k,865) = mat(k,865) + rxt(k,295)*y(k,123) + 2.400_r8*rxt(k,292)*y(k,198) &
                      + rxt(k,293)*y(k,202)
         mat(k,811) = mat(k,811) + rxt(k,323)*y(k,123) + rxt(k,321)*y(k,202)
         mat(k,1420) = mat(k,1420) + rxt(k,391)*y(k,100) + .900_r8*rxt(k,304)*y(k,202) &
                      + rxt(k,377)*y(k,209) + rxt(k,382)*y(k,210) + .470_r8*rxt(k,344) &
                      *y(k,213) + rxt(k,402)*y(k,229)
         mat(k,2176) = mat(k,2176) + rxt(k,194)*y(k,58) + 1.200_r8*rxt(k,392)*y(k,100) &
                      + rxt(k,274)*y(k,123) + rxt(k,293)*y(k,198) + rxt(k,321) &
                      *y(k,199) + .900_r8*rxt(k,304)*y(k,201) + 4.000_r8*rxt(k,271) &
                      *y(k,202) + rxt(k,378)*y(k,209) + rxt(k,383)*y(k,210) &
                      + .730_r8*rxt(k,345)*y(k,213) + rxt(k,354)*y(k,215) &
                      + .500_r8*rxt(k,457)*y(k,219) + .300_r8*rxt(k,333)*y(k,224) &
                      + rxt(k,462)*y(k,225) + rxt(k,467)*y(k,226) + .800_r8*rxt(k,403) &
                      *y(k,229)
         mat(k,766) = mat(k,766) + .170_r8*rxt(k,423)*y(k,123) + .070_r8*rxt(k,422) &
                      *y(k,207)
         mat(k,557) = rxt(k,341)*y(k,123)
         mat(k,461) = rxt(k,312)*y(k,134)
         mat(k,786) = mat(k,786) + .250_r8*rxt(k,310)*y(k,123)
         mat(k,1919) = mat(k,1919) + .070_r8*rxt(k,422)*y(k,203) + .160_r8*rxt(k,425) &
                      *y(k,214) + .330_r8*rxt(k,428)*y(k,216)
         mat(k,442) = mat(k,442) + rxt(k,282)*y(k,123)
         mat(k,1285) = mat(k,1285) + .920_r8*rxt(k,380)*y(k,123) + rxt(k,381)*y(k,125) &
                      + rxt(k,377)*y(k,201) + rxt(k,378)*y(k,202)
         mat(k,1320) = mat(k,1320) + .920_r8*rxt(k,386)*y(k,123) + rxt(k,387)*y(k,125) &
                      + rxt(k,382)*y(k,201) + rxt(k,383)*y(k,202)
         mat(k,1343) = mat(k,1343) + .470_r8*rxt(k,348)*y(k,123) + .470_r8*rxt(k,347) &
                      *y(k,125) + .470_r8*rxt(k,344)*y(k,201) + .730_r8*rxt(k,345) &
                      *y(k,202)
         mat(k,727) = mat(k,727) + .400_r8*rxt(k,426)*y(k,123) + .160_r8*rxt(k,425) &
                      *y(k,207)
         mat(k,1387) = mat(k,1387) + rxt(k,354)*y(k,202)
         mat(k,879) = mat(k,879) + .830_r8*rxt(k,429)*y(k,123) + .330_r8*rxt(k,428) &
                      *y(k,207)
         mat(k,1068) = mat(k,1068) + .500_r8*rxt(k,457)*y(k,202)
         mat(k,1621) = rxt(k,284)*y(k,53)
         mat(k,1786) = mat(k,1786) + .650_r8*rxt(k,410)*y(k,6) + rxt(k,235)*y(k,18) &
                      + .350_r8*rxt(k,289)*y(k,23) + rxt(k,296)*y(k,25) + rxt(k,253) &
                      *y(k,42) + rxt(k,256)*y(k,45) + rxt(k,302)*y(k,46) + rxt(k,275) &
                      *y(k,51) + rxt(k,205)*y(k,58) + rxt(k,287)*y(k,61) &
                      + .730_r8*rxt(k,421)*y(k,65) + .500_r8*rxt(k,489)*y(k,66) &
                      + rxt(k,313)*y(k,73) + rxt(k,314)*y(k,74) + rxt(k,155)*y(k,78) &
                      + rxt(k,278)*y(k,85) + rxt(k,279)*y(k,86) + rxt(k,343)*y(k,92) &
                      + rxt(k,328)*y(k,94) + .300_r8*rxt(k,388)*y(k,98) + rxt(k,389) &
                      *y(k,99) + rxt(k,396)*y(k,101) + .200_r8*rxt(k,352)*y(k,105) &
                      + .500_r8*rxt(k,363)*y(k,108) + rxt(k,400)*y(k,114) + rxt(k,401) &
                      *y(k,115) + rxt(k,174)*y(k,125) + rxt(k,158)*y(k,137) &
                      + .800_r8*rxt(k,433)*y(k,145) + rxt(k,479)*y(k,154) &
                      + .200_r8*rxt(k,473)*y(k,182) + .280_r8*rxt(k,441)*y(k,184) &
                      + .380_r8*rxt(k,443)*y(k,186) + .630_r8*rxt(k,449)*y(k,188)
         mat(k,449) = mat(k,449) + rxt(k,432)*y(k,123)
         mat(k,796) = mat(k,796) + rxt(k,331)*y(k,123)
         mat(k,1170) = mat(k,1170) + .300_r8*rxt(k,333)*y(k,202)
         mat(k,1156) = mat(k,1156) + .900_r8*rxt(k,464)*y(k,123) + rxt(k,462)*y(k,202)
         mat(k,1134) = mat(k,1134) + .800_r8*rxt(k,469)*y(k,123) + rxt(k,467)*y(k,202)
         mat(k,742) = mat(k,742) + rxt(k,439)*y(k,123)
         mat(k,1194) = mat(k,1194) + rxt(k,405)*y(k,123) + rxt(k,406)*y(k,125) &
                      + rxt(k,402)*y(k,201) + .800_r8*rxt(k,403)*y(k,202)
         mat(k,759) = mat(k,759) + rxt(k,445)*y(k,123)
         mat(k,505) = mat(k,505) + rxt(k,448)*y(k,123)
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
         mat(k,438) = -(rxt(k,280)*y(k,207) + rxt(k,282)*y(k,123))
         mat(k,1849) = -rxt(k,280)*y(k,208)
         mat(k,2273) = -rxt(k,282)*y(k,208)
         mat(k,1518) = rxt(k,267)*y(k,207)
         mat(k,1849) = mat(k,1849) + rxt(k,267)*y(k,41)
         mat(k,1275) = -(rxt(k,377)*y(k,201) + rxt(k,378)*y(k,202) + rxt(k,379) &
                      *y(k,207) + rxt(k,380)*y(k,123) + rxt(k,381)*y(k,125))
         mat(k,1408) = -rxt(k,377)*y(k,209)
         mat(k,2162) = -rxt(k,378)*y(k,209)
         mat(k,1904) = -rxt(k,379)*y(k,209)
         mat(k,2322) = -rxt(k,380)*y(k,209)
         mat(k,2107) = -rxt(k,381)*y(k,209)
         mat(k,840) = .600_r8*rxt(k,398)*y(k,221)
         mat(k,1770) = .600_r8*rxt(k,398)*y(k,97)
         mat(k,1310) = -(rxt(k,382)*y(k,201) + rxt(k,383)*y(k,202) + rxt(k,384) &
                      *y(k,207) + rxt(k,386)*y(k,123) + rxt(k,387)*y(k,125))
         mat(k,1409) = -rxt(k,382)*y(k,210)
         mat(k,2163) = -rxt(k,383)*y(k,210)
         mat(k,1905) = -rxt(k,384)*y(k,210)
         mat(k,2323) = -rxt(k,386)*y(k,210)
         mat(k,2108) = -rxt(k,387)*y(k,210)
         mat(k,841) = .400_r8*rxt(k,398)*y(k,221)
         mat(k,1771) = .400_r8*rxt(k,398)*y(k,97)
         mat(k,67) = -(rxt(k,514)*y(k,207) + rxt(k,515)*y(k,123))
         mat(k,1829) = -rxt(k,514)*y(k,211)
         mat(k,2262) = -rxt(k,515)*y(k,211)
         mat(k,833) = rxt(k,517)*y(k,221)
         mat(k,1640) = rxt(k,517)*y(k,97)
         mat(k,73) = -(rxt(k,518)*y(k,207) + rxt(k,519)*y(k,123))
         mat(k,1830) = -rxt(k,518)*y(k,212)
         mat(k,2263) = -rxt(k,519)*y(k,212)
         mat(k,74) = rxt(k,520)*y(k,221)
         mat(k,1641) = rxt(k,520)*y(k,103)
         mat(k,1335) = -(rxt(k,344)*y(k,201) + rxt(k,345)*y(k,202) + rxt(k,346) &
                      *y(k,207) + rxt(k,347)*y(k,125) + (rxt(k,348) + rxt(k,349) &
                      ) * y(k,123))
         mat(k,1410) = -rxt(k,344)*y(k,213)
         mat(k,2164) = -rxt(k,345)*y(k,213)
         mat(k,1906) = -rxt(k,346)*y(k,213)
         mat(k,2109) = -rxt(k,347)*y(k,213)
         mat(k,2324) = -(rxt(k,348) + rxt(k,349)) * y(k,213)
         mat(k,1229) = .500_r8*rxt(k,351)*y(k,221)
         mat(k,314) = .200_r8*rxt(k,352)*y(k,221)
         mat(k,1356) = rxt(k,365)*y(k,221)
         mat(k,1772) = .500_r8*rxt(k,351)*y(k,104) + .200_r8*rxt(k,352)*y(k,105) &
                      + rxt(k,365)*y(k,110)
         mat(k,723) = -(rxt(k,425)*y(k,207) + rxt(k,426)*y(k,123) + rxt(k,427) &
                      *y(k,124))
         mat(k,1871) = -rxt(k,425)*y(k,214)
         mat(k,2291) = -rxt(k,426)*y(k,214)
         mat(k,1555) = -rxt(k,427)*y(k,214)
         mat(k,1380) = -(rxt(k,353)*y(k,201) + rxt(k,354)*y(k,202) + rxt(k,355) &
                      *y(k,207) + 4._r8*rxt(k,356)*y(k,215) + rxt(k,357)*y(k,123) &
                      + rxt(k,358)*y(k,125) + rxt(k,366)*y(k,124))
         mat(k,1412) = -rxt(k,353)*y(k,215)
         mat(k,2166) = -rxt(k,354)*y(k,215)
         mat(k,1908) = -rxt(k,355)*y(k,215)
         mat(k,2326) = -rxt(k,357)*y(k,215)
         mat(k,2111) = -rxt(k,358)*y(k,215)
         mat(k,1566) = -rxt(k,366)*y(k,215)
         mat(k,1230) = .500_r8*rxt(k,351)*y(k,221)
         mat(k,315) = .500_r8*rxt(k,352)*y(k,221)
         mat(k,1774) = .500_r8*rxt(k,351)*y(k,104) + .500_r8*rxt(k,352)*y(k,105)
         mat(k,873) = -(rxt(k,428)*y(k,207) + rxt(k,429)*y(k,123) + rxt(k,430) &
                      *y(k,124))
         mat(k,1884) = -rxt(k,428)*y(k,216)
         mat(k,2301) = -rxt(k,429)*y(k,216)
         mat(k,1558) = -rxt(k,430)*y(k,216)
         mat(k,672) = -(rxt(k,359)*y(k,207) + rxt(k,360)*y(k,123))
         mat(k,1867) = -rxt(k,359)*y(k,217)
         mat(k,2289) = -rxt(k,360)*y(k,217)
         mat(k,508) = rxt(k,361)*y(k,221)
         mat(k,319) = rxt(k,362)*y(k,221)
         mat(k,1722) = rxt(k,361)*y(k,106) + rxt(k,362)*y(k,107)
         mat(k,81) = -(rxt(k,522)*y(k,207) + rxt(k,523)*y(k,123))
         mat(k,1831) = -rxt(k,522)*y(k,218)
         mat(k,2264) = -rxt(k,523)*y(k,218)
         mat(k,904) = rxt(k,525)*y(k,221)
         mat(k,1643) = rxt(k,525)*y(k,109)
         mat(k,1060) = -(rxt(k,457)*y(k,202) + rxt(k,458)*y(k,207) + rxt(k,459) &
                      *y(k,123) + rxt(k,460)*y(k,125))
         mat(k,2150) = -rxt(k,457)*y(k,219)
         mat(k,1891) = -rxt(k,458)*y(k,219)
         mat(k,2309) = -rxt(k,459)*y(k,219)
         mat(k,2094) = -rxt(k,460)*y(k,219)
         mat(k,942) = rxt(k,451)*y(k,125)
         mat(k,914) = rxt(k,454)*y(k,125)
         mat(k,2094) = mat(k,2094) + rxt(k,451)*y(k,5) + rxt(k,454)*y(k,109) &
                      + .500_r8*rxt(k,471)*y(k,181)
         mat(k,390) = rxt(k,461)*y(k,221)
         mat(k,996) = .500_r8*rxt(k,471)*y(k,125)
         mat(k,1756) = rxt(k,461)*y(k,127)
         mat(k,1618) = -(rxt(k,127)*y(k,76) + rxt(k,128)*y(k,233) + (rxt(k,130) &
                      + rxt(k,131)) * y(k,134) + (rxt(k,132) + rxt(k,133)) * y(k,137) &
                      + (rxt(k,179) + rxt(k,180)) * y(k,112) + rxt(k,212)*y(k,32) &
                      + rxt(k,213)*y(k,33) + rxt(k,214)*y(k,35) + rxt(k,215)*y(k,36) &
                      + rxt(k,216)*y(k,37) + rxt(k,217)*y(k,38) + rxt(k,218)*y(k,39) &
                      + (rxt(k,219) + rxt(k,220)) * y(k,84) + rxt(k,239)*y(k,34) &
                      + rxt(k,240)*y(k,54) + rxt(k,241)*y(k,77) + (rxt(k,242) &
                      + rxt(k,243)) * y(k,80) + rxt(k,248)*y(k,63) + rxt(k,249) &
                      *y(k,64) + rxt(k,262)*y(k,40) + rxt(k,263)*y(k,42) + rxt(k,264) &
                      *y(k,81) + rxt(k,265)*y(k,82) + rxt(k,266)*y(k,83) + (rxt(k,283) &
                      + rxt(k,284) + rxt(k,285)) * y(k,53) + rxt(k,286)*y(k,85))
         mat(k,1254) = -rxt(k,127)*y(k,220)
         mat(k,2361) = -rxt(k,128)*y(k,220)
         mat(k,1507) = -(rxt(k,130) + rxt(k,131)) * y(k,220)
         mat(k,1981) = -(rxt(k,132) + rxt(k,133)) * y(k,220)
         mat(k,261) = -(rxt(k,179) + rxt(k,180)) * y(k,220)
         mat(k,101) = -rxt(k,212)*y(k,220)
         mat(k,145) = -rxt(k,213)*y(k,220)
         mat(k,116) = -rxt(k,214)*y(k,220)
         mat(k,155) = -rxt(k,215)*y(k,220)
         mat(k,120) = -rxt(k,216)*y(k,220)
         mat(k,160) = -rxt(k,217)*y(k,220)
         mat(k,124) = -rxt(k,218)*y(k,220)
         mat(k,2059) = -(rxt(k,219) + rxt(k,220)) * y(k,220)
         mat(k,151) = -rxt(k,239)*y(k,220)
         mat(k,384) = -rxt(k,240)*y(k,220)
         mat(k,109) = -rxt(k,241)*y(k,220)
         mat(k,817) = -(rxt(k,242) + rxt(k,243)) * y(k,220)
         mat(k,237) = -rxt(k,248)*y(k,220)
         mat(k,251) = -rxt(k,249)*y(k,220)
         mat(k,478) = -rxt(k,262)*y(k,220)
         mat(k,591) = -rxt(k,263)*y(k,220)
         mat(k,246) = -rxt(k,264)*y(k,220)
         mat(k,256) = -rxt(k,265)*y(k,220)
         mat(k,309) = -rxt(k,266)*y(k,220)
         mat(k,1466) = -(rxt(k,283) + rxt(k,284) + rxt(k,285)) * y(k,220)
         mat(k,186) = -rxt(k,286)*y(k,220)
         mat(k,1784) = -(rxt(k,154)*y(k,76) + rxt(k,155)*y(k,78) + rxt(k,156)*y(k,207) &
                      + rxt(k,157)*y(k,133) + rxt(k,158)*y(k,137) + (4._r8*rxt(k,159) &
                      + 4._r8*rxt(k,160)) * y(k,221) + rxt(k,162)*y(k,89) + rxt(k,174) &
                      *y(k,125) + rxt(k,175)*y(k,111) + rxt(k,183)*y(k,124) + rxt(k,184) &
                      *y(k,88) + rxt(k,203)*y(k,59) + (rxt(k,205) + rxt(k,206) &
                      ) * y(k,58) + rxt(k,208)*y(k,84) + rxt(k,211)*y(k,91) + rxt(k,235) &
                      *y(k,18) + rxt(k,237)*y(k,80) + rxt(k,251)*y(k,40) + rxt(k,253) &
                      *y(k,42) + rxt(k,254)*y(k,43) + rxt(k,256)*y(k,45) + rxt(k,258) &
                      *y(k,54) + rxt(k,259)*y(k,81) + rxt(k,260)*y(k,82) + rxt(k,261) &
                      *y(k,83) + rxt(k,270)*y(k,41) + rxt(k,275)*y(k,51) + rxt(k,276) &
                      *y(k,52) + rxt(k,277)*y(k,53) + rxt(k,278)*y(k,85) + rxt(k,279) &
                      *y(k,86) + rxt(k,287)*y(k,61) + rxt(k,289)*y(k,23) + rxt(k,296) &
                      *y(k,25) + rxt(k,297)*y(k,26) + rxt(k,299)*y(k,27) + rxt(k,301) &
                      *y(k,44) + rxt(k,302)*y(k,46) + rxt(k,307)*y(k,49) + rxt(k,308) &
                      *y(k,50) + rxt(k,313)*y(k,73) + rxt(k,314)*y(k,74) + rxt(k,315) &
                      *y(k,142) + rxt(k,316)*y(k,24) + rxt(k,324)*y(k,29) + rxt(k,325) &
                      *y(k,30) + rxt(k,327)*y(k,48) + rxt(k,328)*y(k,94) + rxt(k,329) &
                      *y(k,126) + rxt(k,332)*y(k,149) + rxt(k,336)*y(k,150) + rxt(k,337) &
                      *y(k,28) + rxt(k,338)*y(k,47) + rxt(k,340)*y(k,15) + rxt(k,343) &
                      *y(k,92) + rxt(k,351)*y(k,104) + rxt(k,352)*y(k,105) + rxt(k,361) &
                      *y(k,106) + rxt(k,362)*y(k,107) + rxt(k,363)*y(k,108) + rxt(k,365) &
                      *y(k,110) + rxt(k,368)*y(k,1) + rxt(k,372)*y(k,2) + rxt(k,373) &
                      *y(k,14) + rxt(k,374)*y(k,93) + rxt(k,375)*y(k,95) + rxt(k,376) &
                      *y(k,96) + rxt(k,388)*y(k,98) + rxt(k,389)*y(k,99) + rxt(k,396) &
                      *y(k,101) + rxt(k,398)*y(k,97) + rxt(k,399)*y(k,102) + rxt(k,400) &
                      *y(k,114) + rxt(k,401)*y(k,115) + rxt(k,407)*y(k,185) + rxt(k,410) &
                      *y(k,6) + rxt(k,413)*y(k,7) + rxt(k,414)*y(k,21) + rxt(k,416) &
                      *y(k,22) + rxt(k,420)*y(k,31) + rxt(k,421)*y(k,65) + rxt(k,433) &
                      *y(k,145) + rxt(k,436)*y(k,146) + rxt(k,440)*y(k,183) + rxt(k,441) &
                      *y(k,184) + rxt(k,443)*y(k,186) + rxt(k,446)*y(k,187) + rxt(k,449) &
                      *y(k,188) + rxt(k,450)*y(k,189) + rxt(k,453)*y(k,5) + rxt(k,456) &
                      *y(k,109) + rxt(k,461)*y(k,127) + rxt(k,465)*y(k,178) + rxt(k,466) &
                      *y(k,179) + rxt(k,470)*y(k,180) + rxt(k,472)*y(k,181) + rxt(k,473) &
                      *y(k,182) + (rxt(k,475) + rxt(k,489)) * y(k,66) + rxt(k,477) &
                      *y(k,140) + rxt(k,479)*y(k,154) + rxt(k,483)*y(k,151) + rxt(k,488) &
                      *y(k,153) + rxt(k,491)*y(k,119))
         mat(k,1255) = -rxt(k,154)*y(k,221)
         mat(k,599) = -rxt(k,155)*y(k,221)
         mat(k,1917) = -rxt(k,156)*y(k,221)
         mat(k,2016) = -rxt(k,157)*y(k,221)
         mat(k,1982) = -rxt(k,158)*y(k,221)
         mat(k,472) = -rxt(k,162)*y(k,221)
         mat(k,2121) = -rxt(k,174)*y(k,221)
         mat(k,684) = -rxt(k,175)*y(k,221)
         mat(k,1576) = -rxt(k,183)*y(k,221)
         mat(k,1485) = -rxt(k,184)*y(k,221)
         mat(k,988) = -rxt(k,203)*y(k,221)
         mat(k,2242) = -(rxt(k,205) + rxt(k,206)) * y(k,221)
         mat(k,2060) = -rxt(k,208)*y(k,221)
         mat(k,825) = -rxt(k,211)*y(k,221)
         mat(k,1809) = -rxt(k,235)*y(k,221)
         mat(k,818) = -rxt(k,237)*y(k,221)
         mat(k,479) = -rxt(k,251)*y(k,221)
         mat(k,592) = -rxt(k,253)*y(k,221)
         mat(k,127) = -rxt(k,254)*y(k,221)
         mat(k,370) = -rxt(k,256)*y(k,221)
         mat(k,385) = -rxt(k,258)*y(k,221)
         mat(k,247) = -rxt(k,259)*y(k,221)
         mat(k,257) = -rxt(k,260)*y(k,221)
         mat(k,310) = -rxt(k,261)*y(k,221)
         mat(k,1531) = -rxt(k,270)*y(k,221)
         mat(k,801) = -rxt(k,275)*y(k,221)
         mat(k,453) = -rxt(k,276)*y(k,221)
         mat(k,1467) = -rxt(k,277)*y(k,221)
         mat(k,187) = -rxt(k,278)*y(k,221)
         mat(k,870) = -rxt(k,279)*y(k,221)
         mat(k,1177) = -rxt(k,287)*y(k,221)
         mat(k,298) = -rxt(k,289)*y(k,221)
         mat(k,270) = -rxt(k,296)*y(k,221)
         mat(k,350) = -rxt(k,297)*y(k,221)
         mat(k,302) = -rxt(k,299)*y(k,221)
         mat(k,1051) = -rxt(k,301)*y(k,221)
         mat(k,104) = -rxt(k,302)*y(k,221)
         mat(k,690) = -rxt(k,307)*y(k,221)
         mat(k,609) = -rxt(k,308)*y(k,221)
         mat(k,1114) = -rxt(k,313)*y(k,221)
         mat(k,981) = -rxt(k,314)*y(k,221)
         mat(k,524) = -rxt(k,315)*y(k,221)
         mat(k,540) = -rxt(k,316)*y(k,221)
         mat(k,409) = -rxt(k,324)*y(k,221)
         mat(k,112) = -rxt(k,325)*y(k,221)
         mat(k,1244) = -rxt(k,327)*y(k,221)
         mat(k,1093) = -rxt(k,328)*y(k,221)
         mat(k,855) = -rxt(k,329)*y(k,221)
         mat(k,532) = -rxt(k,332)*y(k,221)
         mat(k,398) = -rxt(k,336)*y(k,221)
         mat(k,1034) = -rxt(k,337)*y(k,221)
         mat(k,974) = -rxt(k,338)*y(k,221)
         mat(k,358) = -rxt(k,340)*y(k,221)
         mat(k,1084) = -rxt(k,343)*y(k,221)
         mat(k,1234) = -rxt(k,351)*y(k,221)
         mat(k,316) = -rxt(k,352)*y(k,221)
         mat(k,511) = -rxt(k,361)*y(k,221)
         mat(k,322) = -rxt(k,362)*y(k,221)
         mat(k,573) = -rxt(k,363)*y(k,221)
         mat(k,1365) = -rxt(k,365)*y(k,221)
         mat(k,669) = -rxt(k,368)*y(k,221)
         mat(k,635) = -rxt(k,372)*y(k,221)
         mat(k,231) = -rxt(k,373)*y(k,221)
         mat(k,227) = -rxt(k,374)*y(k,221)
         mat(k,346) = -rxt(k,375)*y(k,221)
         mat(k,138) = -rxt(k,376)*y(k,221)
         mat(k,583) = -rxt(k,388)*y(k,221)
         mat(k,565) = -rxt(k,389)*y(k,221)
         mat(k,404) = -rxt(k,396)*y(k,221)
         mat(k,845) = -rxt(k,398)*y(k,221)
         mat(k,699) = -rxt(k,399)*y(k,221)
         mat(k,428) = -rxt(k,400)*y(k,221)
         mat(k,1015) = -rxt(k,401)*y(k,221)
         mat(k,208) = -rxt(k,407)*y(k,221)
         mat(k,172) = -rxt(k,410)*y(k,221)
         mat(k,416) = -rxt(k,413)*y(k,221)
         mat(k,243) = -rxt(k,414)*y(k,221)
         mat(k,342) = -rxt(k,416)*y(k,221)
         mat(k,275) = -rxt(k,420)*y(k,221)
         mat(k,200) = -rxt(k,421)*y(k,221)
         mat(k,181) = -rxt(k,433)*y(k,221)
         mat(k,336) = -rxt(k,436)*y(k,221)
         mat(k,658) = -rxt(k,440)*y(k,221)
         mat(k,195) = -rxt(k,441)*y(k,221)
         mat(k,217) = -rxt(k,443)*y(k,221)
         mat(k,714) = -rxt(k,446)*y(k,221)
         mat(k,222) = -rxt(k,449)*y(k,221)
         mat(k,422) = -rxt(k,450)*y(k,221)
         mat(k,950) = -rxt(k,453)*y(k,221)
         mat(k,922) = -rxt(k,456)*y(k,221)
         mat(k,393) = -rxt(k,461)*y(k,221)
         mat(k,645) = -rxt(k,465)*y(k,221)
         mat(k,616) = -rxt(k,466)*y(k,221)
         mat(k,488) = -rxt(k,470)*y(k,221)
         mat(k,1000) = -rxt(k,472)*y(k,221)
         mat(k,1103) = -rxt(k,473)*y(k,221)
         mat(k,291) = -(rxt(k,475) + rxt(k,489)) * y(k,221)
         mat(k,365) = -rxt(k,477)*y(k,221)
         mat(k,895) = -rxt(k,479)*y(k,221)
         mat(k,719) = -rxt(k,483)*y(k,221)
         mat(k,1434) = -rxt(k,488)*y(k,221)
         mat(k,98) = -rxt(k,491)*y(k,221)
         mat(k,950) = mat(k,950) + .630_r8*rxt(k,452)*y(k,137)
         mat(k,298) = mat(k,298) + .650_r8*rxt(k,289)*y(k,221)
         mat(k,540) = mat(k,540) + .130_r8*rxt(k,291)*y(k,137)
         mat(k,350) = mat(k,350) + .500_r8*rxt(k,297)*y(k,221)
         mat(k,1034) = mat(k,1034) + .360_r8*rxt(k,320)*y(k,137)
         mat(k,1531) = mat(k,1531) + rxt(k,269)*y(k,133)
         mat(k,453) = mat(k,453) + .300_r8*rxt(k,276)*y(k,221)
         mat(k,1467) = mat(k,1467) + rxt(k,283)*y(k,220)
         mat(k,2215) = rxt(k,192)*y(k,207)
         mat(k,887) = rxt(k,246)*y(k,233)
         mat(k,2036) = rxt(k,153)*y(k,137) + 2.000_r8*rxt(k,148)*y(k,207)
         mat(k,1255) = mat(k,1255) + rxt(k,145)*y(k,133) + rxt(k,127)*y(k,220)
         mat(k,599) = mat(k,599) + rxt(k,146)*y(k,133)
         mat(k,818) = mat(k,818) + rxt(k,236)*y(k,133) + rxt(k,242)*y(k,220)
         mat(k,2060) = mat(k,2060) + rxt(k,207)*y(k,133) + rxt(k,219)*y(k,220)
         mat(k,187) = mat(k,187) + rxt(k,286)*y(k,220)
         mat(k,772) = rxt(k,238)*y(k,133)
         mat(k,825) = mat(k,825) + rxt(k,210)*y(k,133)
         mat(k,845) = mat(k,845) + .320_r8*rxt(k,397)*y(k,137)
         mat(k,699) = mat(k,699) + .600_r8*rxt(k,399)*y(k,221)
         mat(k,1234) = mat(k,1234) + .240_r8*rxt(k,350)*y(k,137)
         mat(k,316) = mat(k,316) + .100_r8*rxt(k,352)*y(k,221)
         mat(k,922) = mat(k,922) + .630_r8*rxt(k,455)*y(k,137)
         mat(k,1365) = mat(k,1365) + .360_r8*rxt(k,364)*y(k,137)
         mat(k,2335) = rxt(k,176)*y(k,207)
         mat(k,2121) = mat(k,2121) + rxt(k,171)*y(k,207)
         mat(k,2016) = mat(k,2016) + rxt(k,269)*y(k,41) + rxt(k,145)*y(k,76) &
                      + rxt(k,146)*y(k,78) + rxt(k,236)*y(k,80) + rxt(k,207)*y(k,84) &
                      + rxt(k,238)*y(k,90) + rxt(k,210)*y(k,91) + rxt(k,151)*y(k,207)
         mat(k,1982) = mat(k,1982) + .630_r8*rxt(k,452)*y(k,5) + .130_r8*rxt(k,291) &
                      *y(k,24) + .360_r8*rxt(k,320)*y(k,28) + rxt(k,153)*y(k,75) &
                      + .320_r8*rxt(k,397)*y(k,97) + .240_r8*rxt(k,350)*y(k,104) &
                      + .630_r8*rxt(k,455)*y(k,109) + .360_r8*rxt(k,364)*y(k,110) &
                      + rxt(k,152)*y(k,207)
         mat(k,532) = mat(k,532) + .500_r8*rxt(k,332)*y(k,221)
         mat(k,208) = mat(k,208) + .500_r8*rxt(k,407)*y(k,221)
         mat(k,516) = .400_r8*rxt(k,408)*y(k,207)
         mat(k,1419) = .490_r8*rxt(k,305)*y(k,207)
         mat(k,765) = .400_r8*rxt(k,422)*y(k,207)
         mat(k,1917) = mat(k,1917) + rxt(k,192)*y(k,55) + 2.000_r8*rxt(k,148)*y(k,75) &
                      + rxt(k,176)*y(k,123) + rxt(k,171)*y(k,125) + rxt(k,151) &
                      *y(k,133) + rxt(k,152)*y(k,137) + .400_r8*rxt(k,408)*y(k,192) &
                      + .490_r8*rxt(k,305)*y(k,201) + .400_r8*rxt(k,422)*y(k,203) &
                      + .450_r8*rxt(k,355)*y(k,215) + .400_r8*rxt(k,428)*y(k,216) &
                      + .200_r8*rxt(k,359)*y(k,217) + .150_r8*rxt(k,334)*y(k,224)
         mat(k,1386) = .450_r8*rxt(k,355)*y(k,207)
         mat(k,878) = .400_r8*rxt(k,428)*y(k,207)
         mat(k,676) = .200_r8*rxt(k,359)*y(k,207)
         mat(k,1619) = rxt(k,283)*y(k,53) + rxt(k,127)*y(k,76) + rxt(k,242)*y(k,80) &
                      + rxt(k,219)*y(k,84) + rxt(k,286)*y(k,85) + 2.000_r8*rxt(k,128) &
                      *y(k,233)
         mat(k,1784) = mat(k,1784) + .650_r8*rxt(k,289)*y(k,23) + .500_r8*rxt(k,297) &
                      *y(k,26) + .300_r8*rxt(k,276)*y(k,52) + .600_r8*rxt(k,399) &
                      *y(k,102) + .100_r8*rxt(k,352)*y(k,105) + .500_r8*rxt(k,332) &
                      *y(k,149) + .500_r8*rxt(k,407)*y(k,185)
         mat(k,1169) = .150_r8*rxt(k,334)*y(k,207)
         mat(k,2362) = rxt(k,246)*y(k,72) + 2.000_r8*rxt(k,128)*y(k,220)
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
         mat(k,445) = -(rxt(k,431)*y(k,207) + rxt(k,432)*y(k,123))
         mat(k,1850) = -rxt(k,431)*y(k,222)
         mat(k,2274) = -rxt(k,432)*y(k,222)
         mat(k,198) = .200_r8*rxt(k,421)*y(k,221)
         mat(k,179) = .140_r8*rxt(k,433)*y(k,221)
         mat(k,334) = rxt(k,436)*y(k,221)
         mat(k,1694) = .200_r8*rxt(k,421)*y(k,65) + .140_r8*rxt(k,433)*y(k,145) &
                      + rxt(k,436)*y(k,146)
         mat(k,789) = -(rxt(k,330)*y(k,207) + rxt(k,331)*y(k,123))
         mat(k,1877) = -rxt(k,330)*y(k,223)
         mat(k,2296) = -rxt(k,331)*y(k,223)
         mat(k,1020) = rxt(k,337)*y(k,221)
         mat(k,528) = .500_r8*rxt(k,332)*y(k,221)
         mat(k,1733) = rxt(k,337)*y(k,28) + .500_r8*rxt(k,332)*y(k,149)
         mat(k,1164) = -(rxt(k,333)*y(k,202) + rxt(k,334)*y(k,207) + rxt(k,335) &
                      *y(k,123))
         mat(k,2156) = -rxt(k,333)*y(k,224)
         mat(k,1897) = -rxt(k,334)*y(k,224)
         mat(k,2316) = -rxt(k,335)*y(k,224)
         mat(k,945) = .060_r8*rxt(k,452)*y(k,137)
         mat(k,972) = rxt(k,338)*y(k,221)
         mat(k,917) = .060_r8*rxt(k,455)*y(k,137)
         mat(k,1963) = .060_r8*rxt(k,452)*y(k,5) + .060_r8*rxt(k,455)*y(k,109)
         mat(k,395) = rxt(k,336)*y(k,221)
         mat(k,1098) = .150_r8*rxt(k,473)*y(k,221)
         mat(k,1763) = rxt(k,338)*y(k,47) + rxt(k,336)*y(k,150) + .150_r8*rxt(k,473) &
                      *y(k,182)
         mat(k,1147) = -(rxt(k,462)*y(k,202) + rxt(k,463)*y(k,207) + rxt(k,464) &
                      *y(k,123))
         mat(k,2155) = -rxt(k,462)*y(k,225)
         mat(k,1896) = -rxt(k,463)*y(k,225)
         mat(k,2315) = -rxt(k,464)*y(k,225)
         mat(k,2100) = .500_r8*rxt(k,471)*y(k,181)
         mat(k,642) = rxt(k,465)*y(k,221)
         mat(k,998) = .500_r8*rxt(k,471)*y(k,125) + rxt(k,472)*y(k,221)
         mat(k,1762) = rxt(k,465)*y(k,178) + rxt(k,472)*y(k,181)
         mat(k,1124) = -(rxt(k,467)*y(k,202) + rxt(k,468)*y(k,207) + rxt(k,469) &
                      *y(k,123))
         mat(k,2154) = -rxt(k,467)*y(k,226)
         mat(k,1895) = -rxt(k,468)*y(k,226)
         mat(k,2314) = -rxt(k,469)*y(k,226)
         mat(k,944) = rxt(k,453)*y(k,221)
         mat(k,916) = rxt(k,456)*y(k,221)
         mat(k,486) = rxt(k,470)*y(k,221)
         mat(k,1761) = rxt(k,453)*y(k,5) + rxt(k,456)*y(k,109) + rxt(k,470)*y(k,180)
         mat(k,734) = -(rxt(k,438)*y(k,207) + rxt(k,439)*y(k,123))
         mat(k,1872) = -rxt(k,438)*y(k,227)
         mat(k,2292) = -rxt(k,439)*y(k,227)
         mat(k,652) = rxt(k,440)*y(k,221)
         mat(k,194) = .650_r8*rxt(k,441)*y(k,221)
         mat(k,1729) = rxt(k,440)*y(k,183) + .650_r8*rxt(k,441)*y(k,184)
         mat(k,87) = -(rxt(k,528)*y(k,207) + rxt(k,529)*y(k,123))
         mat(k,1832) = -rxt(k,528)*y(k,228)
         mat(k,2265) = -rxt(k,529)*y(k,228)
         mat(k,189) = rxt(k,527)*y(k,221)
         mat(k,1644) = rxt(k,527)*y(k,184)
         mat(k,1187) = -(rxt(k,402)*y(k,201) + rxt(k,403)*y(k,202) + rxt(k,404) &
                      *y(k,207) + rxt(k,405)*y(k,123) + rxt(k,406)*y(k,125))
         mat(k,1404) = -rxt(k,402)*y(k,229)
         mat(k,2158) = -rxt(k,403)*y(k,229)
         mat(k,1899) = -rxt(k,404)*y(k,229)
         mat(k,2318) = -rxt(k,405)*y(k,229)
         mat(k,2103) = -rxt(k,406)*y(k,229)
         mat(k,226) = rxt(k,374)*y(k,221)
         mat(k,345) = rxt(k,375)*y(k,221)
         mat(k,137) = rxt(k,376)*y(k,221)
         mat(k,695) = .400_r8*rxt(k,399)*y(k,221)
         mat(k,207) = .500_r8*rxt(k,407)*y(k,221)
         mat(k,1765) = rxt(k,374)*y(k,93) + rxt(k,375)*y(k,95) + rxt(k,376)*y(k,96) &
                      + .400_r8*rxt(k,399)*y(k,102) + .500_r8*rxt(k,407)*y(k,185)
         mat(k,750) = -(rxt(k,444)*y(k,207) + rxt(k,445)*y(k,123))
         mat(k,1873) = -rxt(k,444)*y(k,230)
         mat(k,2293) = -rxt(k,445)*y(k,230)
         mat(k,214) = .560_r8*rxt(k,443)*y(k,221)
         mat(k,707) = rxt(k,446)*y(k,221)
         mat(k,1730) = .560_r8*rxt(k,443)*y(k,186) + rxt(k,446)*y(k,187)
         mat(k,93) = -(rxt(k,532)*y(k,207) + rxt(k,533)*y(k,123))
         mat(k,1833) = -rxt(k,532)*y(k,231)
         mat(k,2266) = -rxt(k,533)*y(k,231)
         mat(k,209) = rxt(k,531)*y(k,221)
         mat(k,1645) = rxt(k,531)*y(k,186)
         mat(k,500) = -(rxt(k,447)*y(k,207) + rxt(k,448)*y(k,123))
         mat(k,1858) = -rxt(k,447)*y(k,232)
         mat(k,2279) = -rxt(k,448)*y(k,232)
         mat(k,221) = .300_r8*rxt(k,449)*y(k,221)
         mat(k,419) = rxt(k,450)*y(k,221)
         mat(k,1702) = .300_r8*rxt(k,449)*y(k,188) + rxt(k,450)*y(k,189)
         mat(k,2374) = -(rxt(k,128)*y(k,220) + rxt(k,246)*y(k,72) + rxt(k,490) &
                      *y(k,155))
         mat(k,1631) = -rxt(k,128)*y(k,233)
         mat(k,891) = -rxt(k,246)*y(k,233)
         mat(k,267) = -rxt(k,490)*y(k,233)
         mat(k,305) = rxt(k,299)*y(k,221)
         mat(k,411) = rxt(k,324)*y(k,221)
         mat(k,113) = rxt(k,325)*y(k,221)
         mat(k,482) = rxt(k,251)*y(k,221)
         mat(k,1543) = rxt(k,270)*y(k,221)
         mat(k,597) = rxt(k,253)*y(k,221)
         mat(k,129) = rxt(k,254)*y(k,221)
         mat(k,1055) = rxt(k,301)*y(k,221)
         mat(k,375) = rxt(k,256)*y(k,221)
         mat(k,976) = rxt(k,338)*y(k,221)
         mat(k,1248) = rxt(k,327)*y(k,221)
         mat(k,692) = rxt(k,307)*y(k,221)
         mat(k,611) = rxt(k,308)*y(k,221)
         mat(k,456) = rxt(k,276)*y(k,221)
         mat(k,1475) = rxt(k,277)*y(k,221)
         mat(k,2048) = rxt(k,149)*y(k,207)
         mat(k,1262) = rxt(k,154)*y(k,221)
         mat(k,604) = rxt(k,155)*y(k,221)
         mat(k,822) = rxt(k,237)*y(k,221)
         mat(k,312) = rxt(k,261)*y(k,221)
         mat(k,2072) = (rxt(k,542)+rxt(k,547))*y(k,90) + (rxt(k,535)+rxt(k,541) &
                       +rxt(k,546))*y(k,91) + rxt(k,208)*y(k,221)
         mat(k,872) = rxt(k,279)*y(k,221)
         mat(k,1493) = rxt(k,184)*y(k,221)
         mat(k,475) = rxt(k,162)*y(k,221)
         mat(k,777) = (rxt(k,542)+rxt(k,547))*y(k,84)
         mat(k,830) = (rxt(k,535)+rxt(k,541)+rxt(k,546))*y(k,84) + rxt(k,211)*y(k,221)
         mat(k,1238) = .500_r8*rxt(k,351)*y(k,221)
         mat(k,99) = rxt(k,491)*y(k,221)
         mat(k,534) = rxt(k,332)*y(k,221)
         mat(k,399) = rxt(k,336)*y(k,221)
         mat(k,1929) = rxt(k,149)*y(k,75) + rxt(k,156)*y(k,221)
         mat(k,1796) = rxt(k,299)*y(k,27) + rxt(k,324)*y(k,29) + rxt(k,325)*y(k,30) &
                      + rxt(k,251)*y(k,40) + rxt(k,270)*y(k,41) + rxt(k,253)*y(k,42) &
                      + rxt(k,254)*y(k,43) + rxt(k,301)*y(k,44) + rxt(k,256)*y(k,45) &
                      + rxt(k,338)*y(k,47) + rxt(k,327)*y(k,48) + rxt(k,307)*y(k,49) &
                      + rxt(k,308)*y(k,50) + rxt(k,276)*y(k,52) + rxt(k,277)*y(k,53) &
                      + rxt(k,154)*y(k,76) + rxt(k,155)*y(k,78) + rxt(k,237)*y(k,80) &
                      + rxt(k,261)*y(k,83) + rxt(k,208)*y(k,84) + rxt(k,279)*y(k,86) &
                      + rxt(k,184)*y(k,88) + rxt(k,162)*y(k,89) + rxt(k,211)*y(k,91) &
                      + .500_r8*rxt(k,351)*y(k,104) + rxt(k,491)*y(k,119) + rxt(k,332) &
                      *y(k,149) + rxt(k,336)*y(k,150) + rxt(k,156)*y(k,207) &
                      + 2.000_r8*rxt(k,159)*y(k,221)
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
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 73) = mat(k, 73) + lmat(k, 73)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 106) = mat(k, 106) + lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 139) = lmat(k, 139)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 165) = mat(k, 165) + lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = lmat(k, 225)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = lmat(k, 228)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 232) = lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 238) = lmat(k, 238)
         mat(k, 239) = lmat(k, 239)
         mat(k, 240) = lmat(k, 240)
         mat(k, 241) = mat(k, 241) + lmat(k, 241)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 265) = lmat(k, 265)
         mat(k, 266) = lmat(k, 266)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = lmat(k, 285)
         mat(k, 286) = lmat(k, 286)
         mat(k, 287) = lmat(k, 287)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 306) = mat(k, 306) + lmat(k, 306)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 320) = lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = mat(k, 322) + lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = lmat(k, 325)
         mat(k, 326) = lmat(k, 326)
         mat(k, 327) = lmat(k, 327)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 331) = lmat(k, 331)
         mat(k, 332) = mat(k, 332) + lmat(k, 332)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 335) = lmat(k, 335)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 361) = lmat(k, 361)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 373) = lmat(k, 373)
         mat(k, 374) = mat(k, 374) + lmat(k, 374)
         mat(k, 376) = lmat(k, 376)
         mat(k, 377) = lmat(k, 377)
         mat(k, 378) = lmat(k, 378)
         mat(k, 379) = lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = mat(k, 382) + lmat(k, 382)
         mat(k, 383) = mat(k, 383) + lmat(k, 383)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 389) = lmat(k, 389)
         mat(k, 391) = lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = lmat(k, 397)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 400) = mat(k, 400) + lmat(k, 400)
         mat(k, 401) = lmat(k, 401)
         mat(k, 403) = lmat(k, 403)
         mat(k, 405) = mat(k, 405) + lmat(k, 405)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 408) = lmat(k, 408)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 410) = lmat(k, 410)
         mat(k, 412) = mat(k, 412) + lmat(k, 412)
         mat(k, 413) = lmat(k, 413)
         mat(k, 415) = lmat(k, 415)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 420) = lmat(k, 420)
         mat(k, 421) = lmat(k, 421)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 427) = lmat(k, 427)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 440) = lmat(k, 440)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 452) = mat(k, 452) + lmat(k, 452)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 454) = lmat(k, 454)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 460) = lmat(k, 460)
         mat(k, 461) = mat(k, 461) + lmat(k, 461)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = lmat(k, 463)
         mat(k, 464) = lmat(k, 464)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 471) = mat(k, 471) + lmat(k, 471)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 473) = lmat(k, 473)
         mat(k, 474) = lmat(k, 474)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 483) = mat(k, 483) + lmat(k, 483)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 487) = lmat(k, 487)
         mat(k, 488) = mat(k, 488) + lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 509) = lmat(k, 509)
         mat(k, 510) = lmat(k, 510)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 521) = lmat(k, 521)
         mat(k, 523) = lmat(k, 523)
         mat(k, 525) = mat(k, 525) + lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 529) = lmat(k, 529)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 533) = lmat(k, 533)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 544) = lmat(k, 544)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = lmat(k, 546)
         mat(k, 547) = lmat(k, 547)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 551) = mat(k, 551) + lmat(k, 551)
         mat(k, 559) = mat(k, 559) + lmat(k, 559)
         mat(k, 564) = lmat(k, 564)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 570) = lmat(k, 570)
         mat(k, 572) = lmat(k, 572)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 580) = lmat(k, 580)
         mat(k, 585) = lmat(k, 585)
         mat(k, 586) = lmat(k, 586)
         mat(k, 587) = lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 590) = mat(k, 590) + lmat(k, 590)
         mat(k, 595) = lmat(k, 595)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 599) = mat(k, 599) + lmat(k, 599)
         mat(k, 605) = mat(k, 605) + lmat(k, 605)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 610) = lmat(k, 610)
         mat(k, 612) = mat(k, 612) + lmat(k, 612)
         mat(k, 613) = mat(k, 613) + lmat(k, 613)
         mat(k, 614) = lmat(k, 614)
         mat(k, 615) = mat(k, 615) + lmat(k, 615)
         mat(k, 617) = lmat(k, 617)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 626) = lmat(k, 626)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 631) = lmat(k, 631)
         mat(k, 632) = lmat(k, 632)
         mat(k, 634) = lmat(k, 634)
         mat(k, 635) = mat(k, 635) + lmat(k, 635)
         mat(k, 636) = lmat(k, 636)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 638) = lmat(k, 638)
         mat(k, 639) = lmat(k, 639)
         mat(k, 640) = lmat(k, 640)
         mat(k, 641) = lmat(k, 641)
         mat(k, 643) = lmat(k, 643)
         mat(k, 644) = lmat(k, 644)
         mat(k, 645) = mat(k, 645) + lmat(k, 645)
         mat(k, 646) = lmat(k, 646)
         mat(k, 647) = lmat(k, 647)
         mat(k, 648) = lmat(k, 648)
         mat(k, 649) = lmat(k, 649)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 655) = lmat(k, 655)
         mat(k, 657) = lmat(k, 657)
         mat(k, 658) = mat(k, 658) + lmat(k, 658)
         mat(k, 659) = lmat(k, 659)
         mat(k, 660) = lmat(k, 660)
         mat(k, 661) = mat(k, 661) + lmat(k, 661)
         mat(k, 664) = mat(k, 664) + lmat(k, 664)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 667) = mat(k, 667) + lmat(k, 667)
         mat(k, 668) = mat(k, 668) + lmat(k, 668)
         mat(k, 670) = lmat(k, 670)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 680) = mat(k, 680) + lmat(k, 680)
         mat(k, 688) = mat(k, 688) + lmat(k, 688)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 696) = lmat(k, 696)
         mat(k, 697) = lmat(k, 697)
         mat(k, 698) = lmat(k, 698)
         mat(k, 699) = mat(k, 699) + lmat(k, 699)
         mat(k, 700) = lmat(k, 700)
         mat(k, 701) = lmat(k, 701)
         mat(k, 702) = lmat(k, 702)
         mat(k, 703) = lmat(k, 703)
         mat(k, 704) = lmat(k, 704)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 710) = lmat(k, 710)
         mat(k, 712) = lmat(k, 712)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 715) = lmat(k, 715)
         mat(k, 716) = mat(k, 716) + lmat(k, 716)
         mat(k, 723) = mat(k, 723) + lmat(k, 723)
         mat(k, 734) = mat(k, 734) + lmat(k, 734)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 761) = mat(k, 761) + lmat(k, 761)
         mat(k, 770) = mat(k, 770) + lmat(k, 770)
         mat(k, 771) = lmat(k, 771)
         mat(k, 772) = mat(k, 772) + lmat(k, 772)
         mat(k, 780) = mat(k, 780) + lmat(k, 780)
         mat(k, 789) = mat(k, 789) + lmat(k, 789)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 804) = mat(k, 804) + lmat(k, 804)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 816) = mat(k, 816) + lmat(k, 816)
         mat(k, 821) = mat(k, 821) + lmat(k, 821)
         mat(k, 824) = mat(k, 824) + lmat(k, 824)
         mat(k, 825) = mat(k, 825) + lmat(k, 825)
         mat(k, 828) = mat(k, 828) + lmat(k, 828)
         mat(k, 834) = mat(k, 834) + lmat(k, 834)
         mat(k, 850) = mat(k, 850) + lmat(k, 850)
         mat(k, 852) = lmat(k, 852)
         mat(k, 853) = lmat(k, 853)
         mat(k, 854) = mat(k, 854) + lmat(k, 854)
         mat(k, 859) = mat(k, 859) + lmat(k, 859)
         mat(k, 868) = mat(k, 868) + lmat(k, 868)
         mat(k, 873) = mat(k, 873) + lmat(k, 873)
         mat(k, 883) = mat(k, 883) + lmat(k, 883)
         mat(k, 893) = mat(k, 893) + lmat(k, 893)
         mat(k, 894) = lmat(k, 894)
         mat(k, 897) = lmat(k, 897)
         mat(k, 908) = mat(k, 908) + lmat(k, 908)
         mat(k, 936) = mat(k, 936) + lmat(k, 936)
         mat(k, 960) = mat(k, 960) + lmat(k, 960)
         mat(k, 971) = mat(k, 971) + lmat(k, 971)
         mat(k, 973) = lmat(k, 973)
         mat(k, 975) = lmat(k, 975)
         mat(k, 977) = lmat(k, 977)
         mat(k, 978) = mat(k, 978) + lmat(k, 978)
         mat(k, 980) = mat(k, 980) + lmat(k, 980)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 984) = mat(k, 984) + lmat(k, 984)
         mat(k, 985) = mat(k, 985) + lmat(k, 985)
         mat(k, 986) = mat(k, 986) + lmat(k, 986)
         mat(k, 987) = lmat(k, 987)
         mat(k, 991) = mat(k, 991) + lmat(k, 991)
         mat(k, 992) = mat(k, 992) + lmat(k, 992)
         mat(k, 993) = mat(k, 993) + lmat(k, 993)
         mat(k, 995) = mat(k, 995) + lmat(k, 995)
         mat(k, 997) = lmat(k, 997)
         mat(k, 999) = lmat(k, 999)
         mat(k,1001) = lmat(k,1001)
         mat(k,1003) = lmat(k,1003)
         mat(k,1007) = mat(k,1007) + lmat(k,1007)
         mat(k,1012) = lmat(k,1012)
         mat(k,1014) = lmat(k,1014)
         mat(k,1016) = mat(k,1016) + lmat(k,1016)
         mat(k,1023) = mat(k,1023) + lmat(k,1023)
         mat(k,1043) = lmat(k,1043)
         mat(k,1044) = lmat(k,1044)
         mat(k,1046) = lmat(k,1046)
         mat(k,1047) = mat(k,1047) + lmat(k,1047)
         mat(k,1048) = lmat(k,1048)
         mat(k,1052) = lmat(k,1052)
         mat(k,1054) = lmat(k,1054)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1073) = lmat(k,1073)
         mat(k,1074) = lmat(k,1074)
         mat(k,1075) = mat(k,1075) + lmat(k,1075)
         mat(k,1076) = lmat(k,1076)
         mat(k,1077) = lmat(k,1077)
         mat(k,1079) = lmat(k,1079)
         mat(k,1080) = lmat(k,1080)
         mat(k,1081) = lmat(k,1081)
         mat(k,1082) = lmat(k,1082)
         mat(k,1083) = lmat(k,1083)
         mat(k,1085) = mat(k,1085) + lmat(k,1085)
         mat(k,1089) = mat(k,1089) + lmat(k,1089)
         mat(k,1091) = lmat(k,1091)
         mat(k,1092) = lmat(k,1092)
         mat(k,1094) = mat(k,1094) + lmat(k,1094)
         mat(k,1095) = mat(k,1095) + lmat(k,1095)
         mat(k,1096) = mat(k,1096) + lmat(k,1096)
         mat(k,1097) = mat(k,1097) + lmat(k,1097)
         mat(k,1098) = mat(k,1098) + lmat(k,1098)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1102) = mat(k,1102) + lmat(k,1102)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1111) = lmat(k,1111)
         mat(k,1113) = mat(k,1113) + lmat(k,1113)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1147) = mat(k,1147) + lmat(k,1147)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1175) = mat(k,1175) + lmat(k,1175)
         mat(k,1187) = mat(k,1187) + lmat(k,1187)
         mat(k,1209) = mat(k,1209) + lmat(k,1209)
         mat(k,1226) = mat(k,1226) + lmat(k,1226)
         mat(k,1227) = mat(k,1227) + lmat(k,1227)
         mat(k,1230) = mat(k,1230) + lmat(k,1230)
         mat(k,1231) = mat(k,1231) + lmat(k,1231)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1235) = mat(k,1235) + lmat(k,1235)
         mat(k,1239) = mat(k,1239) + lmat(k,1239)
         mat(k,1240) = mat(k,1240) + lmat(k,1240)
         mat(k,1241) = mat(k,1241) + lmat(k,1241)
         mat(k,1245) = lmat(k,1245)
         mat(k,1251) = mat(k,1251) + lmat(k,1251)
         mat(k,1275) = mat(k,1275) + lmat(k,1275)
         mat(k,1292) = lmat(k,1292)
         mat(k,1310) = mat(k,1310) + lmat(k,1310)
         mat(k,1320) = mat(k,1320) + lmat(k,1320)
         mat(k,1335) = mat(k,1335) + lmat(k,1335)
         mat(k,1350) = lmat(k,1350)
         mat(k,1354) = mat(k,1354) + lmat(k,1354)
         mat(k,1357) = mat(k,1357) + lmat(k,1357)
         mat(k,1359) = mat(k,1359) + lmat(k,1359)
         mat(k,1370) = lmat(k,1370)
         mat(k,1380) = mat(k,1380) + lmat(k,1380)
         mat(k,1413) = mat(k,1413) + lmat(k,1413)
         mat(k,1428) = lmat(k,1428)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1438) = mat(k,1438) + lmat(k,1438)
         mat(k,1446) = mat(k,1446) + lmat(k,1446)
         mat(k,1459) = lmat(k,1459)
         mat(k,1460) = lmat(k,1460)
         mat(k,1461) = mat(k,1461) + lmat(k,1461)
         mat(k,1462) = mat(k,1462) + lmat(k,1462)
         mat(k,1465) = mat(k,1465) + lmat(k,1465)
         mat(k,1467) = mat(k,1467) + lmat(k,1467)
         mat(k,1469) = lmat(k,1469)
         mat(k,1470) = mat(k,1470) + lmat(k,1470)
         mat(k,1473) = mat(k,1473) + lmat(k,1473)
         mat(k,1475) = mat(k,1475) + lmat(k,1475)
         mat(k,1480) = mat(k,1480) + lmat(k,1480)
         mat(k,1483) = lmat(k,1483)
         mat(k,1485) = mat(k,1485) + lmat(k,1485)
         mat(k,1504) = mat(k,1504) + lmat(k,1504)
         mat(k,1507) = mat(k,1507) + lmat(k,1507)
         mat(k,1512) = mat(k,1512) + lmat(k,1512)
         mat(k,1522) = mat(k,1522) + lmat(k,1522)
         mat(k,1523) = lmat(k,1523)
         mat(k,1528) = mat(k,1528) + lmat(k,1528)
         mat(k,1536) = mat(k,1536) + lmat(k,1536)
         mat(k,1571) = mat(k,1571) + lmat(k,1571)
         mat(k,1574) = mat(k,1574) + lmat(k,1574)
         mat(k,1576) = mat(k,1576) + lmat(k,1576)
         mat(k,1580) = mat(k,1580) + lmat(k,1580)
         mat(k,1587) = mat(k,1587) + lmat(k,1587)
         mat(k,1618) = mat(k,1618) + lmat(k,1618)
         mat(k,1623) = mat(k,1623) + lmat(k,1623)
         mat(k,1784) = mat(k,1784) + lmat(k,1784)
         mat(k,1803) = mat(k,1803) + lmat(k,1803)
         mat(k,1810) = mat(k,1810) + lmat(k,1810)
         mat(k,1813) = mat(k,1813) + lmat(k,1813)
         mat(k,1919) = mat(k,1919) + lmat(k,1919)
         mat(k,1929) = mat(k,1929) + lmat(k,1929)
         mat(k,1935) = mat(k,1935) + lmat(k,1935)
         mat(k,1978) = mat(k,1978) + lmat(k,1978)
         mat(k,1981) = mat(k,1981) + lmat(k,1981)
         mat(k,1985) = mat(k,1985) + lmat(k,1985)
         mat(k,1986) = mat(k,1986) + lmat(k,1986)
         mat(k,2020) = mat(k,2020) + lmat(k,2020)
         mat(k,2041) = mat(k,2041) + lmat(k,2041)
         mat(k,2065) = mat(k,2065) + lmat(k,2065)
         mat(k,2066) = mat(k,2066) + lmat(k,2066)
         mat(k,2069) = mat(k,2069) + lmat(k,2069)
         mat(k,2116) = mat(k,2116) + lmat(k,2116)
         mat(k,2117) = mat(k,2117) + lmat(k,2117)
         mat(k,2119) = mat(k,2119) + lmat(k,2119)
         mat(k,2125) = mat(k,2125) + lmat(k,2125)
         mat(k,2128) = mat(k,2128) + lmat(k,2128)
         mat(k,2132) = mat(k,2132) + lmat(k,2132)
         mat(k,2182) = mat(k,2182) + lmat(k,2182)
         mat(k,2224) = mat(k,2224) + lmat(k,2224)
         mat(k,2246) = mat(k,2246) + lmat(k,2246)
         mat(k,2251) = mat(k,2251) + lmat(k,2251)
         mat(k,2252) = mat(k,2252) + lmat(k,2252)
         mat(k,2290) = mat(k,2290) + lmat(k,2290)
         mat(k,2339) = mat(k,2339) + lmat(k,2339)
         mat(k,2346) = mat(k,2346) + lmat(k,2346)
         mat(k,2353) = lmat(k,2353)
         mat(k,2361) = mat(k,2361) + lmat(k,2361)
         mat(k,2362) = mat(k,2362) + lmat(k,2362)
         mat(k,2366) = lmat(k,2366)
         mat(k,2367) = lmat(k,2367)
         mat(k,2374) = mat(k,2374) + lmat(k,2374)
         mat(k, 215) = 0._r8
         mat(k, 216) = 0._r8
         mat(k, 255) = 0._r8
         mat(k, 308) = 0._r8
         mat(k, 341) = 0._r8
         mat(k, 433) = 0._r8
         mat(k, 435) = 0._r8
         mat(k, 448) = 0._r8
         mat(k, 493) = 0._r8
         mat(k, 496) = 0._r8
         mat(k, 504) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 628) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 633) = 0._r8
         mat(k, 651) = 0._r8
         mat(k, 653) = 0._r8
         mat(k, 654) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 662) = 0._r8
         mat(k, 663) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 706) = 0._r8
         mat(k, 708) = 0._r8
         mat(k, 709) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 713) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 735) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 741) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 776) = 0._r8
         mat(k, 781) = 0._r8
         mat(k, 782) = 0._r8
         mat(k, 785) = 0._r8
         mat(k, 791) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 798) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 814) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 909) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 926) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 969) = 0._r8
         mat(k, 994) = 0._r8
         mat(k,1005) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1010) = 0._r8
         mat(k,1011) = 0._r8
         mat(k,1013) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1029) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1033) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1045) = 0._r8
         mat(k,1061) = 0._r8
         mat(k,1062) = 0._r8
         mat(k,1063) = 0._r8
         mat(k,1064) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1087) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1112) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1122) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1129) = 0._r8
         mat(k,1130) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1210) = 0._r8
         mat(k,1211) = 0._r8
         mat(k,1212) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1220) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1267) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1278) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1332) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1348) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1360) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1362) = 0._r8
         mat(k,1364) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1372) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1414) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1416) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1481) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1486) = 0._r8
         mat(k,1487) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1503) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1515) = 0._r8
         mat(k,1517) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1527) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1530) = 0._r8
         mat(k,1532) = 0._r8
         mat(k,1534) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1557) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1562) = 0._r8
         mat(k,1563) = 0._r8
         mat(k,1564) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1569) = 0._r8
         mat(k,1570) = 0._r8
         mat(k,1573) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1581) = 0._r8
         mat(k,1582) = 0._r8
         mat(k,1584) = 0._r8
         mat(k,1585) = 0._r8
         mat(k,1588) = 0._r8
         mat(k,1614) = 0._r8
         mat(k,1617) = 0._r8
         mat(k,1626) = 0._r8
         mat(k,1696) = 0._r8
         mat(k,1712) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1731) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1783) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1808) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1814) = 0._r8
         mat(k,1815) = 0._r8
         mat(k,1816) = 0._r8
         mat(k,1817) = 0._r8
         mat(k,1821) = 0._r8
         mat(k,1852) = 0._r8
         mat(k,1853) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1881) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1888) = 0._r8
         mat(k,1892) = 0._r8
         mat(k,1893) = 0._r8
         mat(k,1894) = 0._r8
         mat(k,1901) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1912) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1944) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1960) = 0._r8
         mat(k,1961) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1965) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1977) = 0._r8
         mat(k,1988) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,2000) = 0._r8
         mat(k,2004) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2011) = 0._r8
         mat(k,2015) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2028) = 0._r8
         mat(k,2030) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2034) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2037) = 0._r8
         mat(k,2042) = 0._r8
         mat(k,2043) = 0._r8
         mat(k,2044) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2054) = 0._r8
         mat(k,2056) = 0._r8
         mat(k,2057) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2067) = 0._r8
         mat(k,2068) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2079) = 0._r8
         mat(k,2082) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2092) = 0._r8
         mat(k,2095) = 0._r8
         mat(k,2097) = 0._r8
         mat(k,2099) = 0._r8
         mat(k,2101) = 0._r8
         mat(k,2113) = 0._r8
         mat(k,2114) = 0._r8
         mat(k,2115) = 0._r8
         mat(k,2120) = 0._r8
         mat(k,2122) = 0._r8
         mat(k,2124) = 0._r8
         mat(k,2126) = 0._r8
         mat(k,2127) = 0._r8
         mat(k,2129) = 0._r8
         mat(k,2130) = 0._r8
         mat(k,2131) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2142) = 0._r8
         mat(k,2168) = 0._r8
         mat(k,2169) = 0._r8
         mat(k,2173) = 0._r8
         mat(k,2174) = 0._r8
         mat(k,2175) = 0._r8
         mat(k,2177) = 0._r8
         mat(k,2178) = 0._r8
         mat(k,2179) = 0._r8
         mat(k,2180) = 0._r8
         mat(k,2181) = 0._r8
         mat(k,2186) = 0._r8
         mat(k,2196) = 0._r8
         mat(k,2199) = 0._r8
         mat(k,2200) = 0._r8
         mat(k,2202) = 0._r8
         mat(k,2203) = 0._r8
         mat(k,2204) = 0._r8
         mat(k,2207) = 0._r8
         mat(k,2210) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2214) = 0._r8
         mat(k,2216) = 0._r8
         mat(k,2219) = 0._r8
         mat(k,2226) = 0._r8
         mat(k,2227) = 0._r8
         mat(k,2237) = 0._r8
         mat(k,2241) = 0._r8
         mat(k,2245) = 0._r8
         mat(k,2247) = 0._r8
         mat(k,2249) = 0._r8
         mat(k,2254) = 0._r8
         mat(k,2298) = 0._r8
         mat(k,2329) = 0._r8
         mat(k,2330) = 0._r8
         mat(k,2334) = 0._r8
         mat(k,2340) = 0._r8
         mat(k,2341) = 0._r8
         mat(k,2347) = 0._r8
         mat(k,2352) = 0._r8
         mat(k,2354) = 0._r8
         mat(k,2355) = 0._r8
         mat(k,2356) = 0._r8
         mat(k,2357) = 0._r8
         mat(k,2358) = 0._r8
         mat(k,2359) = 0._r8
         mat(k,2360) = 0._r8
         mat(k,2363) = 0._r8
         mat(k,2364) = 0._r8
         mat(k,2365) = 0._r8
         mat(k,2368) = 0._r8
         mat(k,2369) = 0._r8
         mat(k,2370) = 0._r8
         mat(k,2371) = 0._r8
         mat(k,2372) = 0._r8
         mat(k,2373) = 0._r8
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
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 81) = mat(k, 81) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 94) = mat(k, 94) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 106) = mat(k, 106) - dti(k)
         mat(k, 110) = mat(k, 110) - dti(k)
         mat(k, 114) = mat(k, 114) - dti(k)
         mat(k, 118) = mat(k, 118) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 126) = mat(k, 126) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 133) = mat(k, 133) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 162) = mat(k, 162) - dti(k)
         mat(k, 165) = mat(k, 165) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 190) = mat(k, 190) - dti(k)
         mat(k, 197) = mat(k, 197) - dti(k)
         mat(k, 202) = mat(k, 202) - dti(k)
         mat(k, 206) = mat(k, 206) - dti(k)
         mat(k, 211) = mat(k, 211) - dti(k)
         mat(k, 219) = mat(k, 219) - dti(k)
         mat(k, 224) = mat(k, 224) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 232) = mat(k, 232) - dti(k)
         mat(k, 235) = mat(k, 235) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 241) = mat(k, 241) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 249) = mat(k, 249) - dti(k)
         mat(k, 254) = mat(k, 254) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 264) = mat(k, 264) - dti(k)
         mat(k, 268) = mat(k, 268) - dti(k)
         mat(k, 272) = mat(k, 272) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 288) = mat(k, 288) - dti(k)
         mat(k, 294) = mat(k, 294) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 307) = mat(k, 307) - dti(k)
         mat(k, 313) = mat(k, 313) - dti(k)
         mat(k, 318) = mat(k, 318) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 344) = mat(k, 344) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 352) = mat(k, 352) - dti(k)
         mat(k, 360) = mat(k, 360) - dti(k)
         mat(k, 368) = mat(k, 368) - dti(k)
         mat(k, 376) = mat(k, 376) - dti(k)
         mat(k, 382) = mat(k, 382) - dti(k)
         mat(k, 388) = mat(k, 388) - dti(k)
         mat(k, 394) = mat(k, 394) - dti(k)
         mat(k, 400) = mat(k, 400) - dti(k)
         mat(k, 406) = mat(k, 406) - dti(k)
         mat(k, 412) = mat(k, 412) - dti(k)
         mat(k, 418) = mat(k, 418) - dti(k)
         mat(k, 424) = mat(k, 424) - dti(k)
         mat(k, 432) = mat(k, 432) - dti(k)
         mat(k, 438) = mat(k, 438) - dti(k)
         mat(k, 445) = mat(k, 445) - dti(k)
         mat(k, 451) = mat(k, 451) - dti(k)
         mat(k, 457) = mat(k, 457) - dti(k)
         mat(k, 462) = mat(k, 462) - dti(k)
         mat(k, 465) = mat(k, 465) - dti(k)
         mat(k, 469) = mat(k, 469) - dti(k)
         mat(k, 476) = mat(k, 476) - dti(k)
         mat(k, 483) = mat(k, 483) - dti(k)
         mat(k, 492) = mat(k, 492) - dti(k)
         mat(k, 500) = mat(k, 500) - dti(k)
         mat(k, 507) = mat(k, 507) - dti(k)
         mat(k, 513) = mat(k, 513) - dti(k)
         mat(k, 519) = mat(k, 519) - dti(k)
         mat(k, 527) = mat(k, 527) - dti(k)
         mat(k, 535) = mat(k, 535) - dti(k)
         mat(k, 543) = mat(k, 543) - dti(k)
         mat(k, 551) = mat(k, 551) - dti(k)
         mat(k, 559) = mat(k, 559) - dti(k)
         mat(k, 567) = mat(k, 567) - dti(k)
         mat(k, 576) = mat(k, 576) - dti(k)
         mat(k, 585) = mat(k, 585) - dti(k)
         mat(k, 589) = mat(k, 589) - dti(k)
         mat(k, 598) = mat(k, 598) - dti(k)
         mat(k, 605) = mat(k, 605) - dti(k)
         mat(k, 612) = mat(k, 612) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 627) = mat(k, 627) - dti(k)
         mat(k, 637) = mat(k, 637) - dti(k)
         mat(k, 650) = mat(k, 650) - dti(k)
         mat(k, 661) = mat(k, 661) - dti(k)
         mat(k, 672) = mat(k, 672) - dti(k)
         mat(k, 680) = mat(k, 680) - dti(k)
         mat(k, 688) = mat(k, 688) - dti(k)
         mat(k, 694) = mat(k, 694) - dti(k)
         mat(k, 705) = mat(k, 705) - dti(k)
         mat(k, 716) = mat(k, 716) - dti(k)
         mat(k, 723) = mat(k, 723) - dti(k)
         mat(k, 734) = mat(k, 734) - dti(k)
         mat(k, 750) = mat(k, 750) - dti(k)
         mat(k, 761) = mat(k, 761) - dti(k)
         mat(k, 770) = mat(k, 770) - dti(k)
         mat(k, 780) = mat(k, 780) - dti(k)
         mat(k, 789) = mat(k, 789) - dti(k)
         mat(k, 799) = mat(k, 799) - dti(k)
         mat(k, 804) = mat(k, 804) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 824) = mat(k, 824) - dti(k)
         mat(k, 834) = mat(k, 834) - dti(k)
         mat(k, 850) = mat(k, 850) - dti(k)
         mat(k, 859) = mat(k, 859) - dti(k)
         mat(k, 868) = mat(k, 868) - dti(k)
         mat(k, 873) = mat(k, 873) - dti(k)
         mat(k, 883) = mat(k, 883) - dti(k)
         mat(k, 893) = mat(k, 893) - dti(k)
         mat(k, 908) = mat(k, 908) - dti(k)
         mat(k, 936) = mat(k, 936) - dti(k)
         mat(k, 960) = mat(k, 960) - dti(k)
         mat(k, 971) = mat(k, 971) - dti(k)
         mat(k, 978) = mat(k, 978) - dti(k)
         mat(k, 985) = mat(k, 985) - dti(k)
         mat(k, 995) = mat(k, 995) - dti(k)
         mat(k,1007) = mat(k,1007) - dti(k)
         mat(k,1023) = mat(k,1023) - dti(k)
         mat(k,1043) = mat(k,1043) - dti(k)
         mat(k,1047) = mat(k,1047) - dti(k)
         mat(k,1060) = mat(k,1060) - dti(k)
         mat(k,1075) = mat(k,1075) - dti(k)
         mat(k,1089) = mat(k,1089) - dti(k)
         mat(k,1097) = mat(k,1097) - dti(k)
         mat(k,1110) = mat(k,1110) - dti(k)
         mat(k,1124) = mat(k,1124) - dti(k)
         mat(k,1147) = mat(k,1147) - dti(k)
         mat(k,1164) = mat(k,1164) - dti(k)
         mat(k,1175) = mat(k,1175) - dti(k)
         mat(k,1187) = mat(k,1187) - dti(k)
         mat(k,1209) = mat(k,1209) - dti(k)
         mat(k,1227) = mat(k,1227) - dti(k)
         mat(k,1240) = mat(k,1240) - dti(k)
         mat(k,1251) = mat(k,1251) - dti(k)
         mat(k,1275) = mat(k,1275) - dti(k)
         mat(k,1310) = mat(k,1310) - dti(k)
         mat(k,1335) = mat(k,1335) - dti(k)
         mat(k,1357) = mat(k,1357) - dti(k)
         mat(k,1380) = mat(k,1380) - dti(k)
         mat(k,1413) = mat(k,1413) - dti(k)
         mat(k,1430) = mat(k,1430) - dti(k)
         mat(k,1446) = mat(k,1446) - dti(k)
         mat(k,1462) = mat(k,1462) - dti(k)
         mat(k,1480) = mat(k,1480) - dti(k)
         mat(k,1504) = mat(k,1504) - dti(k)
         mat(k,1528) = mat(k,1528) - dti(k)
         mat(k,1574) = mat(k,1574) - dti(k)
         mat(k,1618) = mat(k,1618) - dti(k)
         mat(k,1784) = mat(k,1784) - dti(k)
         mat(k,1810) = mat(k,1810) - dti(k)
         mat(k,1919) = mat(k,1919) - dti(k)
         mat(k,1985) = mat(k,1985) - dti(k)
         mat(k,2020) = mat(k,2020) - dti(k)
         mat(k,2041) = mat(k,2041) - dti(k)
         mat(k,2066) = mat(k,2066) - dti(k)
         mat(k,2128) = mat(k,2128) - dti(k)
         mat(k,2182) = mat(k,2182) - dti(k)
         mat(k,2224) = mat(k,2224) - dti(k)
         mat(k,2252) = mat(k,2252) - dti(k)
         mat(k,2346) = mat(k,2346) - dti(k)
         mat(k,2374) = mat(k,2374) - dti(k)
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
