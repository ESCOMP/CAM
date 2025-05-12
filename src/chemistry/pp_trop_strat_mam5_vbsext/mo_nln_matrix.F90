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
         mat(k,688) = -(rxt(k,375)*y(k,253))
         mat(k,1759) = -rxt(k,375)*y(k,1)
         mat(k,1903) = rxt(k,378)*y(k,230)
         mat(k,957) = rxt(k,378)*y(k,129)
         mat(k,677) = -(rxt(k,379)*y(k,253))
         mat(k,1758) = -rxt(k,379)*y(k,2)
         mat(k,956) = rxt(k,376)*y(k,242)
         mat(k,2249) = rxt(k,376)*y(k,230)
         mat(k,1015) = -(rxt(k,458)*y(k,131) + rxt(k,459)*y(k,140) + rxt(k,460) &
                      *y(k,253))
         mat(k,2051) = -rxt(k,458)*y(k,6)
         mat(k,2155) = -rxt(k,459)*y(k,6)
         mat(k,1789) = -rxt(k,460)*y(k,6)
         mat(k,82) = -(rxt(k,513)*y(k,242) + rxt(k,514)*y(k,129))
         mat(k,2207) = -rxt(k,513)*y(k,7)
         mat(k,1869) = -rxt(k,514)*y(k,7)
         mat(k,1007) = rxt(k,516)*y(k,253)
         mat(k,1669) = rxt(k,516)*y(k,6)
         mat(k,207) = -(rxt(k,417)*y(k,253))
         mat(k,1689) = -rxt(k,417)*y(k,8)
         mat(k,113) = -(rxt(k,518)*y(k,242) + rxt(k,519)*y(k,129))
         mat(k,2216) = -rxt(k,518)*y(k,9)
         mat(k,1878) = -rxt(k,519)*y(k,9)
         mat(k,206) = rxt(k,517)*y(k,253)
         mat(k,1679) = rxt(k,517)*y(k,8)
         mat(k,472) = -(rxt(k,420)*y(k,253))
         mat(k,1731) = -rxt(k,420)*y(k,10)
         mat(k,534) = rxt(k,418)*y(k,242)
         mat(k,2234) = rxt(k,418)*y(k,231)
         mat(k,208) = .120_r8*rxt(k,417)*y(k,253)
         mat(k,1690) = .120_r8*rxt(k,417)*y(k,8)
         mat(k,1009) = .100_r8*rxt(k,459)*y(k,140)
         mat(k,934) = .100_r8*rxt(k,462)*y(k,140)
         mat(k,2140) = .100_r8*rxt(k,459)*y(k,6) + .100_r8*rxt(k,462)*y(k,116)
         mat(k,1890) = .500_r8*rxt(k,419)*y(k,231) + .200_r8*rxt(k,446)*y(k,259) &
                      + .060_r8*rxt(k,452)*y(k,261)
         mat(k,535) = .500_r8*rxt(k,419)*y(k,129)
         mat(k,765) = .200_r8*rxt(k,446)*y(k,129)
         mat(k,789) = .060_r8*rxt(k,452)*y(k,129)
         mat(k,1884) = .200_r8*rxt(k,446)*y(k,259) + .200_r8*rxt(k,452)*y(k,261)
         mat(k,764) = .200_r8*rxt(k,446)*y(k,129)
         mat(k,787) = .200_r8*rxt(k,452)*y(k,129)
         mat(k,1900) = .200_r8*rxt(k,446)*y(k,259) + .150_r8*rxt(k,452)*y(k,261)
         mat(k,766) = .200_r8*rxt(k,446)*y(k,129)
         mat(k,790) = .150_r8*rxt(k,452)*y(k,129)
         mat(k,1885) = .210_r8*rxt(k,452)*y(k,261)
         mat(k,788) = .210_r8*rxt(k,452)*y(k,129)
         mat(k,266) = -(rxt(k,380)*y(k,253))
         mat(k,1699) = -rxt(k,380)*y(k,17)
         mat(k,1008) = .050_r8*rxt(k,459)*y(k,140)
         mat(k,933) = .050_r8*rxt(k,462)*y(k,140)
         mat(k,2139) = .050_r8*rxt(k,459)*y(k,6) + .050_r8*rxt(k,462)*y(k,116)
         mat(k,394) = -(rxt(k,346)*y(k,131) + rxt(k,347)*y(k,253))
         mat(k,2041) = -rxt(k,346)*y(k,18)
         mat(k,1720) = -rxt(k,347)*y(k,18)
         mat(k,1462) = -(rxt(k,230)*y(k,44) + rxt(k,231)*y(k,242) + rxt(k,232) &
                      *y(k,140))
         mat(k,1528) = -rxt(k,230)*y(k,19)
         mat(k,2295) = -rxt(k,231)*y(k,19)
         mat(k,2177) = -rxt(k,232)*y(k,19)
         mat(k,1840) = 4.000_r8*rxt(k,233)*y(k,21) + (rxt(k,234)+rxt(k,235))*y(k,61) &
                      + rxt(k,238)*y(k,129) + rxt(k,241)*y(k,139) + rxt(k,488) &
                      *y(k,158) + rxt(k,242)*y(k,253)
         mat(k,182) = rxt(k,220)*y(k,252)
         mat(k,188) = rxt(k,246)*y(k,252)
         mat(k,509) = 2.000_r8*rxt(k,257)*y(k,58) + 2.000_r8*rxt(k,269)*y(k,252) &
                      + 2.000_r8*rxt(k,258)*y(k,253)
         mat(k,627) = rxt(k,259)*y(k,58) + rxt(k,270)*y(k,252) + rxt(k,260)*y(k,253)
         mat(k,461) = 3.000_r8*rxt(k,264)*y(k,58) + 3.000_r8*rxt(k,247)*y(k,252) &
                      + 3.000_r8*rxt(k,265)*y(k,253)
         mat(k,2115) = 2.000_r8*rxt(k,257)*y(k,43) + rxt(k,259)*y(k,45) &
                      + 3.000_r8*rxt(k,264)*y(k,57)
         mat(k,2322) = (rxt(k,234)+rxt(k,235))*y(k,21)
         mat(k,152) = 2.000_r8*rxt(k,248)*y(k,252)
         mat(k,847) = rxt(k,243)*y(k,139) + rxt(k,249)*y(k,252) + rxt(k,244)*y(k,253)
         mat(k,1942) = rxt(k,238)*y(k,21)
         mat(k,1602) = rxt(k,241)*y(k,21) + rxt(k,243)*y(k,83)
         mat(k,1277) = rxt(k,488)*y(k,21)
         mat(k,1643) = rxt(k,220)*y(k,36) + rxt(k,246)*y(k,37) + 2.000_r8*rxt(k,269) &
                      *y(k,43) + rxt(k,270)*y(k,45) + 3.000_r8*rxt(k,247)*y(k,57) &
                      + 2.000_r8*rxt(k,248)*y(k,80) + rxt(k,249)*y(k,83)
         mat(k,1815) = rxt(k,242)*y(k,21) + 2.000_r8*rxt(k,258)*y(k,43) + rxt(k,260) &
                      *y(k,45) + 3.000_r8*rxt(k,265)*y(k,57) + rxt(k,244)*y(k,83)
         mat(k,1834) = rxt(k,236)*y(k,61)
         mat(k,2316) = rxt(k,236)*y(k,21)
         mat(k,2014) = (rxt(k,553)+rxt(k,558))*y(k,93)
         mat(k,822) = (rxt(k,553)+rxt(k,558))*y(k,87)
         mat(k,1848) = -(4._r8*rxt(k,233)*y(k,21) + (rxt(k,234) + rxt(k,235) + rxt(k,236) &
                      ) * y(k,61) + rxt(k,237)*y(k,242) + rxt(k,238)*y(k,129) &
                      + rxt(k,239)*y(k,130) + rxt(k,241)*y(k,139) + rxt(k,242) &
                      *y(k,253) + rxt(k,488)*y(k,158))
         mat(k,2330) = -(rxt(k,234) + rxt(k,235) + rxt(k,236)) * y(k,21)
         mat(k,2304) = -rxt(k,237)*y(k,21)
         mat(k,1951) = -rxt(k,238)*y(k,21)
         mat(k,1580) = -rxt(k,239)*y(k,21)
         mat(k,1611) = -rxt(k,241)*y(k,21)
         mat(k,1824) = -rxt(k,242)*y(k,21)
         mat(k,1282) = -rxt(k,488)*y(k,21)
         mat(k,1468) = rxt(k,232)*y(k,140)
         mat(k,583) = rxt(k,240)*y(k,139)
         mat(k,852) = rxt(k,250)*y(k,252)
         mat(k,827) = rxt(k,245)*y(k,139)
         mat(k,1611) = mat(k,1611) + rxt(k,240)*y(k,22) + rxt(k,245)*y(k,93)
         mat(k,2186) = rxt(k,232)*y(k,19)
         mat(k,1652) = rxt(k,250)*y(k,83)
         mat(k,577) = -(rxt(k,240)*y(k,139))
         mat(k,1592) = -rxt(k,240)*y(k,22)
         mat(k,1836) = rxt(k,239)*y(k,130)
         mat(k,1554) = rxt(k,239)*y(k,21)
         mat(k,275) = -(rxt(k,421)*y(k,253))
         mat(k,1701) = -rxt(k,421)*y(k,24)
         mat(k,1882) = rxt(k,424)*y(k,232)
         mat(k,484) = rxt(k,424)*y(k,129)
         mat(k,371) = -(rxt(k,423)*y(k,253))
         mat(k,1715) = -rxt(k,423)*y(k,25)
         mat(k,485) = rxt(k,422)*y(k,242)
         mat(k,2224) = rxt(k,422)*y(k,232)
         mat(k,321) = -(rxt(k,295)*y(k,58) + rxt(k,296)*y(k,253))
         mat(k,2096) = -rxt(k,295)*y(k,26)
         mat(k,1710) = -rxt(k,296)*y(k,26)
         mat(k,601) = -(rxt(k,297)*y(k,58) + rxt(k,298)*y(k,140) + rxt(k,323)*y(k,253))
         mat(k,2101) = -rxt(k,297)*y(k,27)
         mat(k,2143) = -rxt(k,298)*y(k,27)
         mat(k,1748) = -rxt(k,323)*y(k,27)
         mat(k,306) = -(rxt(k,303)*y(k,253))
         mat(k,1708) = -rxt(k,303)*y(k,28)
         mat(k,862) = .800_r8*rxt(k,299)*y(k,233) + .200_r8*rxt(k,300)*y(k,237)
         mat(k,1961) = .200_r8*rxt(k,300)*y(k,233)
         mat(k,376) = -(rxt(k,304)*y(k,253))
         mat(k,1716) = -rxt(k,304)*y(k,29)
         mat(k,863) = rxt(k,301)*y(k,242)
         mat(k,2225) = rxt(k,301)*y(k,233)
         mat(k,327) = -(rxt(k,305)*y(k,58) + rxt(k,306)*y(k,253))
         mat(k,2097) = -rxt(k,305)*y(k,30)
         mat(k,1711) = -rxt(k,306)*y(k,30)
         mat(k,1048) = -(rxt(k,326)*y(k,131) + rxt(k,327)*y(k,140) + rxt(k,344) &
                      *y(k,253))
         mat(k,2053) = -rxt(k,326)*y(k,31)
         mat(k,2157) = -rxt(k,327)*y(k,31)
         mat(k,1791) = -rxt(k,344)*y(k,31)
         mat(k,881) = .130_r8*rxt(k,404)*y(k,140)
         mat(k,2157) = mat(k,2157) + .130_r8*rxt(k,404)*y(k,100)
         mat(k,436) = -(rxt(k,331)*y(k,253))
         mat(k,1725) = -rxt(k,331)*y(k,32)
         mat(k,835) = rxt(k,329)*y(k,242)
         mat(k,2230) = rxt(k,329)*y(k,234)
         mat(k,154) = -(rxt(k,332)*y(k,253))
         mat(k,1686) = -rxt(k,332)*y(k,33)
         mat(k,310) = -(rxt(k,427)*y(k,253))
         mat(k,1709) = -rxt(k,427)*y(k,34)
         mat(k,668) = rxt(k,425)*y(k,242)
         mat(k,2221) = rxt(k,425)*y(k,235)
         mat(k,141) = -(rxt(k,219)*y(k,252))
         mat(k,1621) = -rxt(k,219)*y(k,35)
         mat(k,180) = -(rxt(k,220)*y(k,252))
         mat(k,1626) = -rxt(k,220)*y(k,36)
         mat(k,185) = -(rxt(k,246)*y(k,252))
         mat(k,1627) = -rxt(k,246)*y(k,37)
         mat(k,158) = -(rxt(k,221)*y(k,252))
         mat(k,1623) = -rxt(k,221)*y(k,38)
         mat(k,190) = -(rxt(k,222)*y(k,252))
         mat(k,1628) = -rxt(k,222)*y(k,39)
         mat(k,162) = -(rxt(k,223)*y(k,252))
         mat(k,1624) = -rxt(k,223)*y(k,40)
         mat(k,195) = -(rxt(k,224)*y(k,252))
         mat(k,1629) = -rxt(k,224)*y(k,41)
         mat(k,166) = -(rxt(k,225)*y(k,252))
         mat(k,1625) = -rxt(k,225)*y(k,42)
         mat(k,508) = -(rxt(k,257)*y(k,58) + rxt(k,258)*y(k,253) + rxt(k,269)*y(k,252))
         mat(k,2100) = -rxt(k,257)*y(k,43)
         mat(k,1736) = -rxt(k,258)*y(k,43)
         mat(k,1638) = -rxt(k,269)*y(k,43)
         mat(k,1532) = -(rxt(k,194)*y(k,58) + rxt(k,230)*y(k,19) + rxt(k,274)*y(k,242) &
                      + rxt(k,275)*y(k,131) + rxt(k,276)*y(k,139) + rxt(k,277) &
                      *y(k,253))
         mat(k,2119) = -rxt(k,194)*y(k,44)
         mat(k,1464) = -rxt(k,230)*y(k,44)
         mat(k,2299) = -rxt(k,274)*y(k,44)
         mat(k,2080) = -rxt(k,275)*y(k,44)
         mat(k,1606) = -rxt(k,276)*y(k,44)
         mat(k,1819) = -rxt(k,277)*y(k,44)
         mat(k,694) = .400_r8*rxt(k,375)*y(k,253)
         mat(k,1025) = .340_r8*rxt(k,459)*y(k,140)
         mat(k,398) = .500_r8*rxt(k,346)*y(k,131)
         mat(k,605) = rxt(k,298)*y(k,140)
         mat(k,1055) = .500_r8*rxt(k,327)*y(k,140)
         mat(k,638) = .500_r8*rxt(k,315)*y(k,253)
         mat(k,832) = rxt(k,282)*y(k,253)
         mat(k,456) = .300_r8*rxt(k,283)*y(k,253)
         mat(k,1480) = (rxt(k,291)+rxt(k,292))*y(k,252)
         mat(k,2325) = rxt(k,201)*y(k,237)
         mat(k,1094) = .800_r8*rxt(k,320)*y(k,253)
         mat(k,889) = .910_r8*rxt(k,404)*y(k,140)
         mat(k,622) = .300_r8*rxt(k,395)*y(k,253)
         mat(k,1243) = .800_r8*rxt(k,399)*y(k,237)
         mat(k,1258) = .120_r8*rxt(k,357)*y(k,140)
         mat(k,657) = .500_r8*rxt(k,370)*y(k,253)
         mat(k,949) = .340_r8*rxt(k,462)*y(k,140)
         mat(k,1384) = .600_r8*rxt(k,371)*y(k,140)
         mat(k,1946) = .100_r8*rxt(k,377)*y(k,230) + rxt(k,281)*y(k,237) &
                      + .500_r8*rxt(k,348)*y(k,239) + .500_r8*rxt(k,317)*y(k,241) &
                      + .920_r8*rxt(k,387)*y(k,244) + .250_r8*rxt(k,355)*y(k,246) &
                      + rxt(k,364)*y(k,248) + rxt(k,338)*y(k,255) + rxt(k,342) &
                      *y(k,256) + .340_r8*rxt(k,471)*y(k,257) + .320_r8*rxt(k,476) &
                      *y(k,258) + .250_r8*rxt(k,412)*y(k,260)
         mat(k,2080) = mat(k,2080) + .500_r8*rxt(k,346)*y(k,18) + rxt(k,388)*y(k,244) &
                      + .250_r8*rxt(k,354)*y(k,246) + rxt(k,365)*y(k,248)
         mat(k,2181) = .340_r8*rxt(k,459)*y(k,6) + rxt(k,298)*y(k,27) &
                      + .500_r8*rxt(k,327)*y(k,31) + .910_r8*rxt(k,404)*y(k,100) &
                      + .120_r8*rxt(k,357)*y(k,111) + .340_r8*rxt(k,462)*y(k,116) &
                      + .600_r8*rxt(k,371)*y(k,118)
         mat(k,572) = rxt(k,322)*y(k,253)
         mat(k,1120) = .680_r8*rxt(k,480)*y(k,253)
         mat(k,964) = .100_r8*rxt(k,377)*y(k,129)
         mat(k,867) = .700_r8*rxt(k,300)*y(k,237)
         mat(k,839) = rxt(k,328)*y(k,237)
         mat(k,1436) = rxt(k,311)*y(k,237) + rxt(k,384)*y(k,244) + .250_r8*rxt(k,351) &
                      *y(k,246) + rxt(k,360)*y(k,248) + .250_r8*rxt(k,409)*y(k,260)
         mat(k,1998) = rxt(k,201)*y(k,61) + .800_r8*rxt(k,399)*y(k,103) + rxt(k,281) &
                      *y(k,129) + .700_r8*rxt(k,300)*y(k,233) + rxt(k,328)*y(k,234) &
                      + rxt(k,311)*y(k,236) + (4.000_r8*rxt(k,278)+2.000_r8*rxt(k,279)) &
                      *y(k,237) + 1.500_r8*rxt(k,385)*y(k,244) + .750_r8*rxt(k,390) &
                      *y(k,245) + .880_r8*rxt(k,352)*y(k,246) + 2.000_r8*rxt(k,361) &
                      *y(k,248) + .750_r8*rxt(k,464)*y(k,251) + .800_r8*rxt(k,340) &
                      *y(k,256) + .930_r8*rxt(k,469)*y(k,257) + .950_r8*rxt(k,474) &
                      *y(k,258) + .800_r8*rxt(k,410)*y(k,260)
         mat(k,613) = .500_r8*rxt(k,348)*y(k,129)
         mat(k,753) = .500_r8*rxt(k,317)*y(k,129)
         mat(k,2299) = mat(k,2299) + .450_r8*rxt(k,362)*y(k,248) + .150_r8*rxt(k,341) &
                      *y(k,256)
         mat(k,1307) = .920_r8*rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131) + rxt(k,384) &
                      *y(k,236) + 1.500_r8*rxt(k,385)*y(k,237)
         mat(k,1340) = .750_r8*rxt(k,390)*y(k,237)
         mat(k,1362) = .250_r8*rxt(k,355)*y(k,129) + .250_r8*rxt(k,354)*y(k,131) &
                      + .250_r8*rxt(k,351)*y(k,236) + .880_r8*rxt(k,352)*y(k,237)
         mat(k,1404) = rxt(k,364)*y(k,129) + rxt(k,365)*y(k,131) + rxt(k,360)*y(k,236) &
                      + 2.000_r8*rxt(k,361)*y(k,237) + .450_r8*rxt(k,362)*y(k,242) &
                      + 4.000_r8*rxt(k,363)*y(k,248)
         mat(k,1107) = .750_r8*rxt(k,464)*y(k,237)
         mat(k,1647) = (rxt(k,291)+rxt(k,292))*y(k,56)
         mat(k,1819) = mat(k,1819) + .400_r8*rxt(k,375)*y(k,1) + .500_r8*rxt(k,315) &
                      *y(k,53) + rxt(k,282)*y(k,54) + .300_r8*rxt(k,283)*y(k,55) &
                      + .800_r8*rxt(k,320)*y(k,76) + .300_r8*rxt(k,395)*y(k,101) &
                      + .500_r8*rxt(k,370)*y(k,115) + rxt(k,322)*y(k,145) &
                      + .680_r8*rxt(k,480)*y(k,217)
         mat(k,816) = rxt(k,338)*y(k,129)
         mat(k,1203) = rxt(k,342)*y(k,129) + .800_r8*rxt(k,340)*y(k,237) &
                      + .150_r8*rxt(k,341)*y(k,242)
         mat(k,1167) = .340_r8*rxt(k,471)*y(k,129) + .930_r8*rxt(k,469)*y(k,237)
         mat(k,1188) = .320_r8*rxt(k,476)*y(k,129) + .950_r8*rxt(k,474)*y(k,237)
         mat(k,1220) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,409)*y(k,236) &
                      + .800_r8*rxt(k,410)*y(k,237)
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
         mat(k,626) = -(rxt(k,259)*y(k,58) + rxt(k,260)*y(k,253) + rxt(k,270)*y(k,252))
         mat(k,2102) = -rxt(k,259)*y(k,45)
         mat(k,1751) = -rxt(k,260)*y(k,45)
         mat(k,1639) = -rxt(k,270)*y(k,45)
         mat(k,170) = -(rxt(k,261)*y(k,253))
         mat(k,1687) = -rxt(k,261)*y(k,46)
         mat(k,1081) = -(rxt(k,307)*y(k,131) + rxt(k,308)*y(k,253))
         mat(k,2055) = -rxt(k,307)*y(k,47)
         mat(k,1793) = -rxt(k,308)*y(k,47)
         mat(k,692) = .800_r8*rxt(k,375)*y(k,253)
         mat(k,397) = rxt(k,346)*y(k,131)
         mat(k,307) = rxt(k,303)*y(k,253)
         mat(k,378) = .500_r8*rxt(k,304)*y(k,253)
         mat(k,1049) = .500_r8*rxt(k,327)*y(k,140)
         mat(k,1374) = .100_r8*rxt(k,371)*y(k,140)
         mat(k,1922) = .400_r8*rxt(k,377)*y(k,230) + rxt(k,302)*y(k,233) &
                      + .270_r8*rxt(k,330)*y(k,234) + rxt(k,348)*y(k,239) + rxt(k,367) &
                      *y(k,250) + rxt(k,338)*y(k,255)
         mat(k,2055) = mat(k,2055) + rxt(k,346)*y(k,18)
         mat(k,2158) = .500_r8*rxt(k,327)*y(k,31) + .100_r8*rxt(k,371)*y(k,118)
         mat(k,962) = .400_r8*rxt(k,377)*y(k,129)
         mat(k,866) = rxt(k,302)*y(k,129) + 3.200_r8*rxt(k,299)*y(k,233) &
                      + .800_r8*rxt(k,300)*y(k,237)
         mat(k,838) = .270_r8*rxt(k,330)*y(k,129)
         mat(k,1976) = .800_r8*rxt(k,300)*y(k,233)
         mat(k,611) = rxt(k,348)*y(k,129)
         mat(k,2275) = .200_r8*rxt(k,366)*y(k,250)
         mat(k,700) = rxt(k,367)*y(k,129) + .200_r8*rxt(k,366)*y(k,242)
         mat(k,1793) = mat(k,1793) + .800_r8*rxt(k,375)*y(k,1) + rxt(k,303)*y(k,28) &
                      + .500_r8*rxt(k,304)*y(k,29)
         mat(k,814) = rxt(k,338)*y(k,129)
         mat(k,402) = -(rxt(k,262)*y(k,58) + rxt(k,263)*y(k,253))
         mat(k,2098) = -rxt(k,262)*y(k,48)
         mat(k,1721) = -rxt(k,263)*y(k,48)
         mat(k,144) = -(rxt(k,309)*y(k,253))
         mat(k,1685) = -rxt(k,309)*y(k,49)
         mat(k,982) = -(rxt(k,345)*y(k,253))
         mat(k,1786) = -rxt(k,345)*y(k,50)
         mat(k,691) = .800_r8*rxt(k,375)*y(k,253)
         mat(k,1012) = .520_r8*rxt(k,459)*y(k,140)
         mat(k,396) = .500_r8*rxt(k,346)*y(k,131)
         mat(k,938) = .520_r8*rxt(k,462)*y(k,140)
         mat(k,1917) = .250_r8*rxt(k,377)*y(k,230) + .820_r8*rxt(k,330)*y(k,234) &
                      + .500_r8*rxt(k,348)*y(k,239) + .270_r8*rxt(k,471)*y(k,257) &
                      + .040_r8*rxt(k,476)*y(k,258)
         mat(k,2048) = .500_r8*rxt(k,346)*y(k,18)
         mat(k,2152) = .520_r8*rxt(k,459)*y(k,6) + .520_r8*rxt(k,462)*y(k,116)
         mat(k,1115) = .500_r8*rxt(k,480)*y(k,253)
         mat(k,961) = .250_r8*rxt(k,377)*y(k,129)
         mat(k,837) = .820_r8*rxt(k,330)*y(k,129) + .820_r8*rxt(k,328)*y(k,237)
         mat(k,1972) = .820_r8*rxt(k,328)*y(k,234) + .150_r8*rxt(k,469)*y(k,257) &
                      + .025_r8*rxt(k,474)*y(k,258)
         mat(k,610) = .500_r8*rxt(k,348)*y(k,129)
         mat(k,1786) = mat(k,1786) + .800_r8*rxt(k,375)*y(k,1) + .500_r8*rxt(k,480) &
                      *y(k,217)
         mat(k,1159) = .270_r8*rxt(k,471)*y(k,129) + .150_r8*rxt(k,469)*y(k,237)
         mat(k,1178) = .040_r8*rxt(k,476)*y(k,129) + .025_r8*rxt(k,474)*y(k,237)
         mat(k,1265) = -(rxt(k,333)*y(k,131) + rxt(k,334)*y(k,253))
         mat(k,2068) = -rxt(k,333)*y(k,51)
         mat(k,1806) = -rxt(k,334)*y(k,51)
         mat(k,1150) = rxt(k,335)*y(k,253)
         mat(k,1254) = .880_r8*rxt(k,357)*y(k,140)
         mat(k,1377) = .500_r8*rxt(k,371)*y(k,140)
         mat(k,1935) = .170_r8*rxt(k,430)*y(k,238) + .050_r8*rxt(k,393)*y(k,245) &
                      + .250_r8*rxt(k,355)*y(k,246) + .170_r8*rxt(k,436)*y(k,249) &
                      + .400_r8*rxt(k,446)*y(k,259) + .250_r8*rxt(k,412)*y(k,260) &
                      + .540_r8*rxt(k,452)*y(k,261) + .510_r8*rxt(k,455)*y(k,262)
         mat(k,2068) = mat(k,2068) + .050_r8*rxt(k,394)*y(k,245) + .250_r8*rxt(k,354) &
                      *y(k,246) + .250_r8*rxt(k,413)*y(k,260)
         mat(k,896) = rxt(k,336)*y(k,253)
         mat(k,2169) = .880_r8*rxt(k,357)*y(k,111) + .500_r8*rxt(k,371)*y(k,118)
         mat(k,1427) = .250_r8*rxt(k,351)*y(k,246) + .250_r8*rxt(k,409)*y(k,260)
         mat(k,1988) = .240_r8*rxt(k,352)*y(k,246) + .500_r8*rxt(k,340)*y(k,256) &
                      + .100_r8*rxt(k,410)*y(k,260)
         mat(k,806) = .170_r8*rxt(k,430)*y(k,129) + .070_r8*rxt(k,429)*y(k,242)
         mat(k,2287) = .070_r8*rxt(k,429)*y(k,238) + .070_r8*rxt(k,435)*y(k,249)
         mat(k,1333) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,1357) = .250_r8*rxt(k,355)*y(k,129) + .250_r8*rxt(k,354)*y(k,131) &
                      + .250_r8*rxt(k,351)*y(k,236) + .240_r8*rxt(k,352)*y(k,237)
         mat(k,921) = .170_r8*rxt(k,436)*y(k,129) + .070_r8*rxt(k,435)*y(k,242)
         mat(k,1806) = mat(k,1806) + rxt(k,335)*y(k,97) + rxt(k,336)*y(k,132)
         mat(k,1201) = .500_r8*rxt(k,340)*y(k,237)
         mat(k,774) = .400_r8*rxt(k,446)*y(k,129)
         mat(k,1218) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,236) + .100_r8*rxt(k,410)*y(k,237)
         mat(k,798) = .540_r8*rxt(k,452)*y(k,129)
         mat(k,546) = .510_r8*rxt(k,455)*y(k,129)
         mat(k,729) = -(rxt(k,314)*y(k,253))
         mat(k,1763) = -rxt(k,314)*y(k,52)
         mat(k,1043) = .120_r8*rxt(k,327)*y(k,140)
         mat(k,2145) = .120_r8*rxt(k,327)*y(k,31)
         mat(k,1417) = .100_r8*rxt(k,311)*y(k,237) + .150_r8*rxt(k,312)*y(k,242)
         mat(k,1966) = .100_r8*rxt(k,311)*y(k,236)
         mat(k,2253) = .150_r8*rxt(k,312)*y(k,236) + .150_r8*rxt(k,362)*y(k,248)
         mat(k,1396) = .150_r8*rxt(k,362)*y(k,242)
         mat(k,635) = -(rxt(k,315)*y(k,253))
         mat(k,1752) = -rxt(k,315)*y(k,53)
         mat(k,1416) = .400_r8*rxt(k,312)*y(k,242)
         mat(k,2246) = .400_r8*rxt(k,312)*y(k,236) + .400_r8*rxt(k,362)*y(k,248)
         mat(k,1394) = .400_r8*rxt(k,362)*y(k,242)
         mat(k,831) = -(rxt(k,282)*y(k,253))
         mat(k,1772) = -rxt(k,282)*y(k,54)
         mat(k,1231) = .200_r8*rxt(k,399)*y(k,237)
         mat(k,864) = .300_r8*rxt(k,300)*y(k,237)
         mat(k,1967) = .200_r8*rxt(k,399)*y(k,103) + .300_r8*rxt(k,300)*y(k,233) &
                      + 2.000_r8*rxt(k,279)*y(k,237) + .250_r8*rxt(k,385)*y(k,244) &
                      + .250_r8*rxt(k,390)*y(k,245) + .250_r8*rxt(k,352)*y(k,246) &
                      + .250_r8*rxt(k,464)*y(k,251) + .500_r8*rxt(k,340)*y(k,256) &
                      + .250_r8*rxt(k,469)*y(k,257) + .250_r8*rxt(k,474)*y(k,258) &
                      + .300_r8*rxt(k,410)*y(k,260)
         mat(k,1291) = .250_r8*rxt(k,385)*y(k,237)
         mat(k,1322) = .250_r8*rxt(k,390)*y(k,237)
         mat(k,1351) = .250_r8*rxt(k,352)*y(k,237)
         mat(k,1100) = .250_r8*rxt(k,464)*y(k,237)
         mat(k,1198) = .500_r8*rxt(k,340)*y(k,237)
         mat(k,1157) = .250_r8*rxt(k,469)*y(k,237)
         mat(k,1177) = .250_r8*rxt(k,474)*y(k,237)
         mat(k,1211) = .300_r8*rxt(k,410)*y(k,237)
         mat(k,454) = -(rxt(k,283)*y(k,253))
         mat(k,1728) = -rxt(k,283)*y(k,55)
         mat(k,1963) = rxt(k,280)*y(k,242)
         mat(k,2233) = rxt(k,280)*y(k,237)
         mat(k,1477) = -(rxt(k,195)*y(k,58) + rxt(k,251)*y(k,75) + rxt(k,284)*y(k,253) &
                      + (rxt(k,290) + rxt(k,291) + rxt(k,292)) * y(k,252))
         mat(k,2116) = -rxt(k,195)*y(k,56)
         mat(k,911) = -rxt(k,251)*y(k,56)
         mat(k,1816) = -rxt(k,284)*y(k,56)
         mat(k,1644) = -(rxt(k,290) + rxt(k,291) + rxt(k,292)) * y(k,56)
         mat(k,1054) = .100_r8*rxt(k,327)*y(k,140)
         mat(k,2178) = .100_r8*rxt(k,327)*y(k,31)
         mat(k,460) = -(rxt(k,247)*y(k,252) + rxt(k,264)*y(k,58) + rxt(k,265)*y(k,253))
         mat(k,1637) = -rxt(k,247)*y(k,57)
         mat(k,2099) = -rxt(k,264)*y(k,57)
         mat(k,1729) = -rxt(k,265)*y(k,57)
         mat(k,2129) = -(rxt(k,194)*y(k,44) + rxt(k,195)*y(k,56) + rxt(k,196)*y(k,79) &
                      + rxt(k,197)*y(k,81) + (rxt(k,198) + rxt(k,199)) * y(k,242) &
                      + rxt(k,200)*y(k,140) + rxt(k,207)*y(k,62) + rxt(k,216)*y(k,94) &
                      + rxt(k,257)*y(k,43) + rxt(k,259)*y(k,45) + rxt(k,262)*y(k,48) &
                      + rxt(k,264)*y(k,57) + rxt(k,305)*y(k,30))
         mat(k,1542) = -rxt(k,194)*y(k,58)
         mat(k,1487) = -rxt(k,195)*y(k,58)
         mat(k,1458) = -rxt(k,196)*y(k,58)
         mat(k,650) = -rxt(k,197)*y(k,58)
         mat(k,2309) = -(rxt(k,198) + rxt(k,199)) * y(k,58)
         mat(k,2191) = -rxt(k,200)*y(k,58)
         mat(k,979) = -rxt(k,207)*y(k,58)
         mat(k,859) = -rxt(k,216)*y(k,58)
         mat(k,513) = -rxt(k,257)*y(k,58)
         mat(k,632) = -rxt(k,259)*y(k,58)
         mat(k,407) = -rxt(k,262)*y(k,58)
         mat(k,465) = -rxt(k,264)*y(k,58)
         mat(k,331) = -rxt(k,305)*y(k,58)
         mat(k,1853) = rxt(k,235)*y(k,61)
         mat(k,143) = 4.000_r8*rxt(k,219)*y(k,252)
         mat(k,184) = rxt(k,220)*y(k,252)
         mat(k,161) = 2.000_r8*rxt(k,221)*y(k,252)
         mat(k,194) = 2.000_r8*rxt(k,222)*y(k,252)
         mat(k,165) = 2.000_r8*rxt(k,223)*y(k,252)
         mat(k,199) = rxt(k,224)*y(k,252)
         mat(k,169) = 2.000_r8*rxt(k,225)*y(k,252)
         mat(k,172) = 3.000_r8*rxt(k,261)*y(k,253)
         mat(k,407) = mat(k,407) + rxt(k,263)*y(k,253)
         mat(k,2335) = rxt(k,235)*y(k,21) + (4.000_r8*rxt(k,202)+2.000_r8*rxt(k,204)) &
                      *y(k,61) + rxt(k,206)*y(k,129) + rxt(k,211)*y(k,139) &
                      + rxt(k,489)*y(k,158) + rxt(k,201)*y(k,237) + rxt(k,212) &
                      *y(k,253)
         mat(k,294) = rxt(k,256)*y(k,252)
         mat(k,290) = rxt(k,271)*y(k,252) + rxt(k,266)*y(k,253)
         mat(k,300) = rxt(k,272)*y(k,252) + rxt(k,267)*y(k,253)
         mat(k,338) = rxt(k,273)*y(k,252) + rxt(k,268)*y(k,253)
         mat(k,2031) = rxt(k,214)*y(k,139) + rxt(k,226)*y(k,252) + rxt(k,215)*y(k,253)
         mat(k,1956) = rxt(k,206)*y(k,61)
         mat(k,1616) = rxt(k,211)*y(k,61) + rxt(k,214)*y(k,87)
         mat(k,1284) = rxt(k,489)*y(k,61)
         mat(k,2008) = rxt(k,201)*y(k,61)
         mat(k,1657) = 4.000_r8*rxt(k,219)*y(k,35) + rxt(k,220)*y(k,36) &
                      + 2.000_r8*rxt(k,221)*y(k,38) + 2.000_r8*rxt(k,222)*y(k,39) &
                      + 2.000_r8*rxt(k,223)*y(k,40) + rxt(k,224)*y(k,41) &
                      + 2.000_r8*rxt(k,225)*y(k,42) + rxt(k,256)*y(k,67) + rxt(k,271) &
                      *y(k,84) + rxt(k,272)*y(k,85) + rxt(k,273)*y(k,86) + rxt(k,226) &
                      *y(k,87)
         mat(k,1829) = 3.000_r8*rxt(k,261)*y(k,46) + rxt(k,263)*y(k,48) + rxt(k,212) &
                      *y(k,61) + rxt(k,266)*y(k,84) + rxt(k,267)*y(k,85) + rxt(k,268) &
                      *y(k,86) + rxt(k,215)*y(k,87)
         mat(k,2095) = rxt(k,207)*y(k,62)
         mat(k,2315) = 2.000_r8*rxt(k,203)*y(k,61)
         mat(k,970) = rxt(k,207)*y(k,58) + (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,87)
         mat(k,2013) = (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,62) + (rxt(k,546) &
                       +rxt(k,552)+rxt(k,557))*y(k,94)
         mat(k,854) = (rxt(k,546)+rxt(k,552)+rxt(k,557))*y(k,87)
         mat(k,2314) = 2.000_r8*rxt(k,228)*y(k,61)
         mat(k,2338) = -(rxt(k,201)*y(k,237) + (4._r8*rxt(k,202) + 4._r8*rxt(k,203) &
                      + 4._r8*rxt(k,204) + 4._r8*rxt(k,228)) * y(k,61) + rxt(k,205) &
                      *y(k,242) + rxt(k,206)*y(k,129) + rxt(k,208)*y(k,130) + rxt(k,211) &
                      *y(k,139) + (rxt(k,212) + rxt(k,213)) * y(k,253) + (rxt(k,234) &
                      + rxt(k,235) + rxt(k,236)) * y(k,21) + rxt(k,489)*y(k,158))
         mat(k,2011) = -rxt(k,201)*y(k,61)
         mat(k,2312) = -rxt(k,205)*y(k,61)
         mat(k,1959) = -rxt(k,206)*y(k,61)
         mat(k,1588) = -rxt(k,208)*y(k,61)
         mat(k,1619) = -rxt(k,211)*y(k,61)
         mat(k,1832) = -(rxt(k,212) + rxt(k,213)) * y(k,61)
         mat(k,1856) = -(rxt(k,234) + rxt(k,235) + rxt(k,236)) * y(k,61)
         mat(k,1287) = -rxt(k,489)*y(k,61)
         mat(k,2132) = rxt(k,216)*y(k,94) + rxt(k,200)*y(k,140) + rxt(k,199)*y(k,242)
         mat(k,980) = rxt(k,209)*y(k,139)
         mat(k,2034) = rxt(k,227)*y(k,252)
         mat(k,860) = rxt(k,216)*y(k,58) + rxt(k,217)*y(k,139) + rxt(k,218)*y(k,253)
         mat(k,1619) = mat(k,1619) + rxt(k,209)*y(k,62) + rxt(k,217)*y(k,94)
         mat(k,2194) = rxt(k,200)*y(k,58)
         mat(k,364) = rxt(k,494)*y(k,158)
         mat(k,1287) = mat(k,1287) + rxt(k,494)*y(k,142)
         mat(k,2312) = mat(k,2312) + rxt(k,199)*y(k,58)
         mat(k,1660) = rxt(k,227)*y(k,87)
         mat(k,1832) = mat(k,1832) + rxt(k,218)*y(k,94)
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
         mat(k,972) = -(rxt(k,207)*y(k,58) + rxt(k,209)*y(k,139) + rxt(k,210)*y(k,253) &
                      + (rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,87))
         mat(k,2107) = -rxt(k,207)*y(k,62)
         mat(k,1598) = -rxt(k,209)*y(k,62)
         mat(k,1785) = -rxt(k,210)*y(k,62)
         mat(k,2017) = -(rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,62)
         mat(k,2320) = rxt(k,208)*y(k,130)
         mat(k,1563) = rxt(k,208)*y(k,61)
         mat(k,1127) = -(rxt(k,294)*y(k,253))
         mat(k,1797) = -rxt(k,294)*y(k,64)
         mat(k,1020) = .230_r8*rxt(k,459)*y(k,140)
         mat(k,1461) = rxt(k,230)*y(k,44)
         mat(k,324) = .350_r8*rxt(k,296)*y(k,253)
         mat(k,604) = .630_r8*rxt(k,298)*y(k,140)
         mat(k,1050) = .560_r8*rxt(k,327)*y(k,140)
         mat(k,1526) = rxt(k,230)*y(k,19) + rxt(k,194)*y(k,58) + rxt(k,275)*y(k,131) &
                      + rxt(k,276)*y(k,139) + rxt(k,277)*y(k,253)
         mat(k,403) = rxt(k,262)*y(k,58)
         mat(k,1264) = rxt(k,333)*y(k,131) + rxt(k,334)*y(k,253)
         mat(k,2112) = rxt(k,194)*y(k,44) + rxt(k,262)*y(k,48)
         mat(k,998) = rxt(k,321)*y(k,253)
         mat(k,882) = .620_r8*rxt(k,404)*y(k,140)
         mat(k,1252) = .650_r8*rxt(k,357)*y(k,140)
         mat(k,944) = .230_r8*rxt(k,462)*y(k,140)
         mat(k,1375) = .560_r8*rxt(k,371)*y(k,140)
         mat(k,1926) = .170_r8*rxt(k,430)*y(k,238) + .220_r8*rxt(k,355)*y(k,246) &
                      + .400_r8*rxt(k,433)*y(k,247) + .350_r8*rxt(k,436)*y(k,249) &
                      + .225_r8*rxt(k,471)*y(k,257) + .250_r8*rxt(k,412)*y(k,260)
         mat(k,2059) = rxt(k,275)*y(k,44) + rxt(k,333)*y(k,51) + .220_r8*rxt(k,354) &
                      *y(k,246) + .500_r8*rxt(k,413)*y(k,260)
         mat(k,1599) = rxt(k,276)*y(k,44) + rxt(k,483)*y(k,143)
         mat(k,2162) = .230_r8*rxt(k,459)*y(k,6) + .630_r8*rxt(k,298)*y(k,27) &
                      + .560_r8*rxt(k,327)*y(k,31) + .620_r8*rxt(k,404)*y(k,100) &
                      + .650_r8*rxt(k,357)*y(k,111) + .230_r8*rxt(k,462)*y(k,116) &
                      + .560_r8*rxt(k,371)*y(k,118)
         mat(k,413) = rxt(k,483)*y(k,139) + rxt(k,484)*y(k,253)
         mat(k,1117) = .700_r8*rxt(k,480)*y(k,253)
         mat(k,1422) = .220_r8*rxt(k,351)*y(k,246) + .250_r8*rxt(k,409)*y(k,260)
         mat(k,1980) = .110_r8*rxt(k,352)*y(k,246) + .125_r8*rxt(k,469)*y(k,257) &
                      + .200_r8*rxt(k,410)*y(k,260)
         mat(k,805) = .170_r8*rxt(k,430)*y(k,129) + .070_r8*rxt(k,429)*y(k,242)
         mat(k,2279) = .070_r8*rxt(k,429)*y(k,238) + .160_r8*rxt(k,432)*y(k,247) &
                      + .140_r8*rxt(k,435)*y(k,249)
         mat(k,1353) = .220_r8*rxt(k,355)*y(k,129) + .220_r8*rxt(k,354)*y(k,131) &
                      + .220_r8*rxt(k,351)*y(k,236) + .110_r8*rxt(k,352)*y(k,237)
         mat(k,760) = .400_r8*rxt(k,433)*y(k,129) + .160_r8*rxt(k,432)*y(k,242)
         mat(k,920) = .350_r8*rxt(k,436)*y(k,129) + .140_r8*rxt(k,435)*y(k,242)
         mat(k,1797) = mat(k,1797) + .350_r8*rxt(k,296)*y(k,26) + rxt(k,277)*y(k,44) &
                      + rxt(k,334)*y(k,51) + rxt(k,321)*y(k,77) + rxt(k,484)*y(k,143) &
                      + .700_r8*rxt(k,480)*y(k,217)
         mat(k,1162) = .225_r8*rxt(k,471)*y(k,129) + .125_r8*rxt(k,469)*y(k,237)
         mat(k,1215) = .250_r8*rxt(k,412)*y(k,129) + .500_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,236) + .200_r8*rxt(k,410)*y(k,237)
         mat(k,1010) = .270_r8*rxt(k,459)*y(k,140)
         mat(k,1045) = .200_r8*rxt(k,327)*y(k,140)
         mat(k,730) = rxt(k,314)*y(k,253)
         mat(k,636) = .500_r8*rxt(k,315)*y(k,253)
         mat(k,1126) = rxt(k,294)*y(k,253)
         mat(k,1090) = .800_r8*rxt(k,320)*y(k,253)
         mat(k,996) = rxt(k,321)*y(k,253)
         mat(k,988) = rxt(k,286)*y(k,253)
         mat(k,654) = .500_r8*rxt(k,370)*y(k,253)
         mat(k,935) = .270_r8*rxt(k,462)*y(k,140)
         mat(k,1371) = .100_r8*rxt(k,371)*y(k,140)
         mat(k,1913) = rxt(k,313)*y(k,236) + .900_r8*rxt(k,471)*y(k,257)
         mat(k,2147) = .270_r8*rxt(k,459)*y(k,6) + .200_r8*rxt(k,327)*y(k,31) &
                      + .270_r8*rxt(k,462)*y(k,116) + .100_r8*rxt(k,371)*y(k,118)
         mat(k,1114) = 1.800_r8*rxt(k,480)*y(k,253)
         mat(k,1418) = rxt(k,313)*y(k,129) + 4.000_r8*rxt(k,310)*y(k,236) &
                      + .900_r8*rxt(k,311)*y(k,237) + rxt(k,384)*y(k,244) &
                      + 2.000_r8*rxt(k,360)*y(k,248) + rxt(k,409)*y(k,260)
         mat(k,1970) = .900_r8*rxt(k,311)*y(k,236) + rxt(k,361)*y(k,248) &
                      + .500_r8*rxt(k,469)*y(k,257)
         mat(k,2267) = .450_r8*rxt(k,362)*y(k,248)
         mat(k,1292) = rxt(k,384)*y(k,236)
         mat(k,1397) = 2.000_r8*rxt(k,360)*y(k,236) + rxt(k,361)*y(k,237) &
                      + .450_r8*rxt(k,362)*y(k,242) + 4.000_r8*rxt(k,363)*y(k,248)
         mat(k,1777) = rxt(k,314)*y(k,52) + .500_r8*rxt(k,315)*y(k,53) + rxt(k,294) &
                      *y(k,64) + .800_r8*rxt(k,320)*y(k,76) + rxt(k,321)*y(k,77) &
                      + rxt(k,286)*y(k,89) + .500_r8*rxt(k,370)*y(k,115) &
                      + 1.800_r8*rxt(k,480)*y(k,217)
         mat(k,1158) = .900_r8*rxt(k,471)*y(k,129) + .500_r8*rxt(k,469)*y(k,237)
         mat(k,1212) = rxt(k,409)*y(k,236)
         mat(k,283) = -(rxt(k,255)*y(k,252))
         mat(k,1632) = -rxt(k,255)*y(k,66)
         mat(k,181) = rxt(k,220)*y(k,252)
         mat(k,186) = rxt(k,246)*y(k,252)
         mat(k,191) = rxt(k,222)*y(k,252)
         mat(k,163) = 2.000_r8*rxt(k,223)*y(k,252)
         mat(k,196) = 2.000_r8*rxt(k,224)*y(k,252)
         mat(k,167) = rxt(k,225)*y(k,252)
         mat(k,151) = 2.000_r8*rxt(k,248)*y(k,252)
         mat(k,295) = rxt(k,272)*y(k,252) + rxt(k,267)*y(k,253)
         mat(k,333) = rxt(k,273)*y(k,252) + rxt(k,268)*y(k,253)
         mat(k,1632) = mat(k,1632) + rxt(k,220)*y(k,36) + rxt(k,246)*y(k,37) &
                      + rxt(k,222)*y(k,39) + 2.000_r8*rxt(k,223)*y(k,40) &
                      + 2.000_r8*rxt(k,224)*y(k,41) + rxt(k,225)*y(k,42) &
                      + 2.000_r8*rxt(k,248)*y(k,80) + rxt(k,272)*y(k,85) + rxt(k,273) &
                      *y(k,86)
         mat(k,1703) = rxt(k,267)*y(k,85) + rxt(k,268)*y(k,86)
         mat(k,291) = -(rxt(k,256)*y(k,252))
         mat(k,1634) = -rxt(k,256)*y(k,67)
         mat(k,159) = rxt(k,221)*y(k,252)
         mat(k,192) = rxt(k,222)*y(k,252)
         mat(k,287) = rxt(k,271)*y(k,252) + rxt(k,266)*y(k,253)
         mat(k,1634) = mat(k,1634) + rxt(k,221)*y(k,38) + rxt(k,222)*y(k,39) &
                      + rxt(k,271)*y(k,84)
         mat(k,1705) = rxt(k,266)*y(k,84)
         mat(k,239) = -(rxt(k,428)*y(k,253))
         mat(k,1694) = -rxt(k,428)*y(k,68)
         mat(k,233) = .180_r8*rxt(k,448)*y(k,253)
         mat(k,1694) = mat(k,1694) + .180_r8*rxt(k,448)*y(k,219)
         mat(k,349) = -(rxt(k,481)*y(k,131) + (rxt(k,482) + rxt(k,496)) * y(k,253))
         mat(k,2039) = -rxt(k,481)*y(k,69)
         mat(k,1713) = -(rxt(k,482) + rxt(k,496)) * y(k,69)
         mat(k,749) = rxt(k,316)*y(k,242)
         mat(k,2219) = rxt(k,316)*y(k,241)
         mat(k,909) = -(rxt(k,251)*y(k,56) + rxt(k,252)*y(k,79) + rxt(k,253)*y(k,263) &
                      + rxt(k,254)*y(k,91))
         mat(k,1474) = -rxt(k,251)*y(k,75)
         mat(k,1447) = -rxt(k,252)*y(k,75)
         mat(k,2344) = -rxt(k,253)*y(k,75)
         mat(k,1491) = -rxt(k,254)*y(k,75)
         mat(k,187) = rxt(k,246)*y(k,252)
         mat(k,197) = rxt(k,224)*y(k,252)
         mat(k,284) = 2.000_r8*rxt(k,255)*y(k,252)
         mat(k,292) = rxt(k,256)*y(k,252)
         mat(k,1641) = rxt(k,246)*y(k,37) + rxt(k,224)*y(k,41) + 2.000_r8*rxt(k,255) &
                      *y(k,66) + rxt(k,256)*y(k,67)
         mat(k,1092) = -(rxt(k,320)*y(k,253))
         mat(k,1794) = -rxt(k,320)*y(k,76)
         mat(k,618) = .700_r8*rxt(k,395)*y(k,253)
         mat(k,595) = .500_r8*rxt(k,396)*y(k,253)
         mat(k,468) = rxt(k,407)*y(k,253)
         mat(k,1923) = .050_r8*rxt(k,393)*y(k,245) + .530_r8*rxt(k,355)*y(k,246) &
                      + .225_r8*rxt(k,471)*y(k,257) + .250_r8*rxt(k,412)*y(k,260)
         mat(k,2056) = .050_r8*rxt(k,394)*y(k,245) + .530_r8*rxt(k,354)*y(k,246) &
                      + .250_r8*rxt(k,413)*y(k,260)
         mat(k,1421) = .530_r8*rxt(k,351)*y(k,246) + .250_r8*rxt(k,409)*y(k,260)
         mat(k,1977) = .260_r8*rxt(k,352)*y(k,246) + .125_r8*rxt(k,469)*y(k,257) &
                      + .100_r8*rxt(k,410)*y(k,260)
         mat(k,1326) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,1352) = .530_r8*rxt(k,355)*y(k,129) + .530_r8*rxt(k,354)*y(k,131) &
                      + .530_r8*rxt(k,351)*y(k,236) + .260_r8*rxt(k,352)*y(k,237)
         mat(k,1794) = mat(k,1794) + .700_r8*rxt(k,395)*y(k,101) + .500_r8*rxt(k,396) &
                      *y(k,102) + rxt(k,407)*y(k,122)
         mat(k,1160) = .225_r8*rxt(k,471)*y(k,129) + .125_r8*rxt(k,469)*y(k,237)
         mat(k,1214) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,236) + .100_r8*rxt(k,410)*y(k,237)
         mat(k,997) = -(rxt(k,321)*y(k,253))
         mat(k,1788) = -rxt(k,321)*y(k,77)
         mat(k,323) = .650_r8*rxt(k,296)*y(k,253)
         mat(k,1091) = .200_r8*rxt(k,320)*y(k,253)
         mat(k,1068) = rxt(k,408)*y(k,253)
         mat(k,1919) = rxt(k,419)*y(k,231) + .050_r8*rxt(k,393)*y(k,245) &
                      + .400_r8*rxt(k,433)*y(k,247) + .170_r8*rxt(k,436)*y(k,249) &
                      + .700_r8*rxt(k,439)*y(k,254) + .600_r8*rxt(k,446)*y(k,259) &
                      + .250_r8*rxt(k,412)*y(k,260) + .340_r8*rxt(k,452)*y(k,261) &
                      + .170_r8*rxt(k,455)*y(k,262)
         mat(k,2050) = .050_r8*rxt(k,394)*y(k,245) + .250_r8*rxt(k,413)*y(k,260)
         mat(k,538) = rxt(k,419)*y(k,129)
         mat(k,1419) = .250_r8*rxt(k,409)*y(k,260)
         mat(k,1973) = .100_r8*rxt(k,410)*y(k,260)
         mat(k,2273) = .160_r8*rxt(k,432)*y(k,247) + .070_r8*rxt(k,435)*y(k,249)
         mat(k,1325) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,759) = .400_r8*rxt(k,433)*y(k,129) + .160_r8*rxt(k,432)*y(k,242)
         mat(k,919) = .170_r8*rxt(k,436)*y(k,129) + .070_r8*rxt(k,435)*y(k,242)
         mat(k,1788) = mat(k,1788) + .650_r8*rxt(k,296)*y(k,26) + .200_r8*rxt(k,320) &
                      *y(k,76) + rxt(k,408)*y(k,123)
         mat(k,500) = .700_r8*rxt(k,439)*y(k,129)
         mat(k,772) = .600_r8*rxt(k,446)*y(k,129)
         mat(k,1213) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,236) + .100_r8*rxt(k,410)*y(k,237)
         mat(k,796) = .340_r8*rxt(k,452)*y(k,129)
         mat(k,545) = .170_r8*rxt(k,455)*y(k,129)
         mat(k,1510) = -((rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,242) + rxt(k,160) &
                      *y(k,140))
         mat(k,2298) = -(rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,78)
         mat(k,2180) = -rxt(k,160)*y(k,78)
         mat(k,1531) = rxt(k,277)*y(k,253)
         mat(k,1479) = rxt(k,291)*y(k,252)
         mat(k,2118) = rxt(k,196)*y(k,79)
         mat(k,913) = rxt(k,252)*y(k,79)
         mat(k,1451) = rxt(k,196)*y(k,58) + rxt(k,252)*y(k,75) + rxt(k,152)*y(k,139) &
                      + rxt(k,144)*y(k,252) + rxt(k,161)*y(k,253)
         mat(k,848) = rxt(k,250)*y(k,252)
         mat(k,2020) = rxt(k,227)*y(k,252)
         mat(k,528) = rxt(k,182)*y(k,253)
         mat(k,1605) = rxt(k,152)*y(k,79) + rxt(k,164)*y(k,253)
         mat(k,415) = rxt(k,484)*y(k,253)
         mat(k,558) = rxt(k,490)*y(k,253)
         mat(k,1278) = rxt(k,495)*y(k,253)
         mat(k,1646) = rxt(k,291)*y(k,56) + rxt(k,144)*y(k,79) + rxt(k,250)*y(k,83) &
                      + rxt(k,227)*y(k,87)
         mat(k,1818) = rxt(k,277)*y(k,44) + rxt(k,161)*y(k,79) + rxt(k,182)*y(k,119) &
                      + rxt(k,164)*y(k,139) + rxt(k,484)*y(k,143) + rxt(k,490) &
                      *y(k,156) + rxt(k,495)*y(k,158)
         mat(k,1448) = -(rxt(k,144)*y(k,252) + rxt(k,152)*y(k,139) + rxt(k,161) &
                      *y(k,253) + rxt(k,196)*y(k,58) + rxt(k,252)*y(k,75))
         mat(k,1642) = -rxt(k,144)*y(k,79)
         mat(k,1601) = -rxt(k,152)*y(k,79)
         mat(k,1814) = -rxt(k,161)*y(k,79)
         mat(k,2114) = -rxt(k,196)*y(k,79)
         mat(k,910) = -rxt(k,252)*y(k,79)
         mat(k,1476) = rxt(k,292)*y(k,252)
         mat(k,1507) = rxt(k,154)*y(k,242)
         mat(k,2294) = rxt(k,154)*y(k,78)
         mat(k,1642) = mat(k,1642) + rxt(k,292)*y(k,56)
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
         mat(k,150) = -(rxt(k,248)*y(k,252))
         mat(k,1622) = -rxt(k,248)*y(k,80)
         mat(k,646) = -(rxt(k,153)*y(k,139) + rxt(k,162)*y(k,253) + rxt(k,197)*y(k,58))
         mat(k,1593) = -rxt(k,153)*y(k,81)
         mat(k,1754) = -rxt(k,162)*y(k,81)
         mat(k,2103) = -rxt(k,197)*y(k,81)
         mat(k,2247) = 2.000_r8*rxt(k,168)*y(k,242)
         mat(k,1754) = mat(k,1754) + 2.000_r8*rxt(k,167)*y(k,253)
         mat(k,301) = rxt(k,497)*y(k,263)
         mat(k,2340) = rxt(k,497)*y(k,160)
         mat(k,846) = -(rxt(k,243)*y(k,139) + rxt(k,244)*y(k,253) + (rxt(k,249) &
                      + rxt(k,250)) * y(k,252))
         mat(k,1595) = -rxt(k,243)*y(k,83)
         mat(k,1774) = -rxt(k,244)*y(k,83)
         mat(k,1640) = -(rxt(k,249) + rxt(k,250)) * y(k,83)
         mat(k,1460) = rxt(k,230)*y(k,44) + rxt(k,231)*y(k,242)
         mat(k,1524) = rxt(k,230)*y(k,19)
         mat(k,2264) = rxt(k,231)*y(k,19)
         mat(k,286) = -(rxt(k,266)*y(k,253) + rxt(k,271)*y(k,252))
         mat(k,1704) = -rxt(k,266)*y(k,84)
         mat(k,1633) = -rxt(k,271)*y(k,84)
         mat(k,296) = -(rxt(k,267)*y(k,253) + rxt(k,272)*y(k,252))
         mat(k,1706) = -rxt(k,267)*y(k,85)
         mat(k,1635) = -rxt(k,272)*y(k,85)
         mat(k,334) = -(rxt(k,268)*y(k,253) + rxt(k,273)*y(k,252))
         mat(k,1712) = -rxt(k,268)*y(k,86)
         mat(k,1636) = -rxt(k,273)*y(k,86)
         mat(k,2029) = -(rxt(k,214)*y(k,139) + rxt(k,215)*y(k,253) + (rxt(k,226) &
                      + rxt(k,227)) * y(k,252) + (rxt(k,546) + rxt(k,552) + rxt(k,557) &
                      ) * y(k,94) + (rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,62) &
                      + (rxt(k,553) + rxt(k,558)) * y(k,93))
         mat(k,1614) = -rxt(k,214)*y(k,87)
         mat(k,1827) = -rxt(k,215)*y(k,87)
         mat(k,1655) = -(rxt(k,226) + rxt(k,227)) * y(k,87)
         mat(k,858) = -(rxt(k,546) + rxt(k,552) + rxt(k,557)) * y(k,87)
         mat(k,977) = -(rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,87)
         mat(k,828) = -(rxt(k,553) + rxt(k,558)) * y(k,87)
         mat(k,330) = rxt(k,305)*y(k,58)
         mat(k,512) = rxt(k,257)*y(k,58)
         mat(k,1540) = rxt(k,194)*y(k,58)
         mat(k,631) = rxt(k,259)*y(k,58)
         mat(k,406) = 2.000_r8*rxt(k,262)*y(k,58)
         mat(k,1485) = rxt(k,195)*y(k,58)
         mat(k,464) = rxt(k,264)*y(k,58)
         mat(k,2127) = rxt(k,305)*y(k,30) + rxt(k,257)*y(k,43) + rxt(k,194)*y(k,44) &
                      + rxt(k,259)*y(k,45) + 2.000_r8*rxt(k,262)*y(k,48) + rxt(k,195) &
                      *y(k,56) + rxt(k,264)*y(k,57) + rxt(k,196)*y(k,79) + rxt(k,197) &
                      *y(k,81) + rxt(k,216)*y(k,94) + rxt(k,198)*y(k,242)
         mat(k,2333) = rxt(k,213)*y(k,253)
         mat(k,1456) = rxt(k,196)*y(k,58)
         mat(k,649) = rxt(k,197)*y(k,58)
         mat(k,858) = mat(k,858) + rxt(k,216)*y(k,58)
         mat(k,2307) = rxt(k,198)*y(k,58)
         mat(k,1827) = mat(k,1827) + rxt(k,213)*y(k,61)
         mat(k,227) = -(rxt(k,285)*y(k,253) + rxt(k,293)*y(k,252))
         mat(k,1692) = -rxt(k,285)*y(k,88)
         mat(k,1631) = -rxt(k,293)*y(k,88)
         mat(k,989) = -(rxt(k,286)*y(k,253))
         mat(k,1787) = -rxt(k,286)*y(k,89)
         mat(k,1013) = .050_r8*rxt(k,459)*y(k,140)
         mat(k,322) = .350_r8*rxt(k,296)*y(k,253)
         mat(k,603) = .370_r8*rxt(k,298)*y(k,140)
         mat(k,1047) = .120_r8*rxt(k,327)*y(k,140)
         mat(k,880) = .110_r8*rxt(k,404)*y(k,140)
         mat(k,1251) = .330_r8*rxt(k,357)*y(k,140)
         mat(k,939) = .050_r8*rxt(k,462)*y(k,140)
         mat(k,1372) = .120_r8*rxt(k,371)*y(k,140)
         mat(k,1918) = rxt(k,289)*y(k,243)
         mat(k,2153) = .050_r8*rxt(k,459)*y(k,6) + .370_r8*rxt(k,298)*y(k,27) &
                      + .120_r8*rxt(k,327)*y(k,31) + .110_r8*rxt(k,404)*y(k,100) &
                      + .330_r8*rxt(k,357)*y(k,111) + .050_r8*rxt(k,462)*y(k,116) &
                      + .120_r8*rxt(k,371)*y(k,118)
         mat(k,2272) = rxt(k,287)*y(k,243)
         mat(k,493) = rxt(k,289)*y(k,129) + rxt(k,287)*y(k,242)
         mat(k,1787) = mat(k,1787) + .350_r8*rxt(k,296)*y(k,26)
         mat(k,1472) = rxt(k,251)*y(k,75)
         mat(k,908) = rxt(k,251)*y(k,56) + rxt(k,252)*y(k,79) + rxt(k,254)*y(k,91) &
                      + rxt(k,253)*y(k,263)
         mat(k,1446) = rxt(k,252)*y(k,75)
         mat(k,1490) = rxt(k,254)*y(k,75)
         mat(k,2342) = rxt(k,253)*y(k,75)
         mat(k,1494) = -(rxt(k,191)*y(k,253) + rxt(k,254)*y(k,75))
         mat(k,1817) = -rxt(k,191)*y(k,91)
         mat(k,912) = -rxt(k,254)*y(k,91)
         mat(k,1530) = rxt(k,275)*y(k,131)
         mat(k,1084) = rxt(k,307)*y(k,131)
         mat(k,1267) = rxt(k,333)*y(k,131)
         mat(k,973) = (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,87)
         mat(k,351) = rxt(k,481)*y(k,131)
         mat(k,2019) = (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,62)
         mat(k,1573) = rxt(k,190)*y(k,253)
         mat(k,2078) = rxt(k,275)*y(k,44) + rxt(k,307)*y(k,47) + rxt(k,333)*y(k,51) &
                      + rxt(k,481)*y(k,69)
         mat(k,1817) = mat(k,1817) + rxt(k,190)*y(k,130)
         mat(k,442) = -(rxt(k,169)*y(k,253))
         mat(k,1726) = -rxt(k,169)*y(k,92)
         mat(k,1549) = rxt(k,188)*y(k,242)
         mat(k,2231) = rxt(k,188)*y(k,130)
         mat(k,823) = -(rxt(k,245)*y(k,139) + (rxt(k,553) + rxt(k,558)) * y(k,87))
         mat(k,1594) = -rxt(k,245)*y(k,93)
         mat(k,2015) = -(rxt(k,553) + rxt(k,558)) * y(k,93)
         mat(k,1837) = rxt(k,237)*y(k,242)
         mat(k,2262) = rxt(k,237)*y(k,21)
         mat(k,855) = -(rxt(k,216)*y(k,58) + rxt(k,217)*y(k,139) + rxt(k,218)*y(k,253) &
                      + (rxt(k,546) + rxt(k,552) + rxt(k,557)) * y(k,87))
         mat(k,2105) = -rxt(k,216)*y(k,94)
         mat(k,1596) = -rxt(k,217)*y(k,94)
         mat(k,1775) = -rxt(k,218)*y(k,94)
         mat(k,2016) = -(rxt(k,546) + rxt(k,552) + rxt(k,557)) * y(k,94)
         mat(k,2318) = rxt(k,205)*y(k,242)
         mat(k,971) = rxt(k,210)*y(k,253)
         mat(k,2265) = rxt(k,205)*y(k,61)
         mat(k,1775) = mat(k,1775) + rxt(k,210)*y(k,62)
         mat(k,1136) = -(rxt(k,350)*y(k,253))
         mat(k,1798) = -rxt(k,350)*y(k,95)
         mat(k,619) = .300_r8*rxt(k,395)*y(k,253)
         mat(k,596) = .500_r8*rxt(k,396)*y(k,253)
         mat(k,1927) = rxt(k,349)*y(k,239) + rxt(k,356)*y(k,246)
         mat(k,612) = rxt(k,349)*y(k,129)
         mat(k,1354) = rxt(k,356)*y(k,129)
         mat(k,1798) = mat(k,1798) + .300_r8*rxt(k,395)*y(k,101) + .500_r8*rxt(k,396) &
                      *y(k,102)
         mat(k,278) = -(rxt(k,381)*y(k,253))
         mat(k,1702) = -rxt(k,381)*y(k,96)
         mat(k,1149) = -(rxt(k,335)*y(k,253))
         mat(k,1799) = -rxt(k,335)*y(k,97)
         mat(k,620) = .700_r8*rxt(k,395)*y(k,253)
         mat(k,597) = .500_r8*rxt(k,396)*y(k,253)
         mat(k,655) = .500_r8*rxt(k,370)*y(k,253)
         mat(k,1928) = .050_r8*rxt(k,393)*y(k,245) + .220_r8*rxt(k,355)*y(k,246) &
                      + .250_r8*rxt(k,412)*y(k,260)
         mat(k,2061) = .050_r8*rxt(k,394)*y(k,245) + .220_r8*rxt(k,354)*y(k,246) &
                      + .250_r8*rxt(k,413)*y(k,260)
         mat(k,588) = .500_r8*rxt(k,339)*y(k,253)
         mat(k,1423) = .220_r8*rxt(k,351)*y(k,246) + .250_r8*rxt(k,409)*y(k,260)
         mat(k,1981) = .230_r8*rxt(k,352)*y(k,246) + .200_r8*rxt(k,340)*y(k,256) &
                      + .100_r8*rxt(k,410)*y(k,260)
         mat(k,1329) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,1355) = .220_r8*rxt(k,355)*y(k,129) + .220_r8*rxt(k,354)*y(k,131) &
                      + .220_r8*rxt(k,351)*y(k,236) + .230_r8*rxt(k,352)*y(k,237)
         mat(k,1799) = mat(k,1799) + .700_r8*rxt(k,395)*y(k,101) + .500_r8*rxt(k,396) &
                      *y(k,102) + .500_r8*rxt(k,370)*y(k,115) + .500_r8*rxt(k,339) &
                      *y(k,154)
         mat(k,1199) = .200_r8*rxt(k,340)*y(k,237)
         mat(k,1216) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,236) + .100_r8*rxt(k,410)*y(k,237)
         mat(k,381) = -(rxt(k,382)*y(k,253))
         mat(k,1717) = -rxt(k,382)*y(k,98)
         mat(k,1886) = .870_r8*rxt(k,393)*y(k,245)
         mat(k,2040) = .950_r8*rxt(k,394)*y(k,245)
         mat(k,1414) = rxt(k,389)*y(k,245)
         mat(k,1962) = .750_r8*rxt(k,390)*y(k,245)
         mat(k,1318) = .870_r8*rxt(k,393)*y(k,129) + .950_r8*rxt(k,394)*y(k,131) &
                      + rxt(k,389)*y(k,236) + .750_r8*rxt(k,390)*y(k,237)
         mat(k,174) = -(rxt(k,383)*y(k,253))
         mat(k,1688) = -rxt(k,383)*y(k,99)
         mat(k,779) = .600_r8*rxt(k,406)*y(k,253)
         mat(k,1688) = mat(k,1688) + .600_r8*rxt(k,406)*y(k,106)
         mat(k,879) = -(rxt(k,397)*y(k,131) + rxt(k,404)*y(k,140) + rxt(k,405) &
                      *y(k,253))
         mat(k,2043) = -rxt(k,397)*y(k,100)
         mat(k,2148) = -rxt(k,404)*y(k,100)
         mat(k,1778) = -rxt(k,405)*y(k,100)
         mat(k,617) = -(rxt(k,395)*y(k,253))
         mat(k,1750) = -rxt(k,395)*y(k,101)
         mat(k,1899) = .080_r8*rxt(k,387)*y(k,244)
         mat(k,1289) = .080_r8*rxt(k,387)*y(k,129)
         mat(k,593) = -(rxt(k,396)*y(k,253))
         mat(k,1747) = -rxt(k,396)*y(k,102)
         mat(k,1897) = .080_r8*rxt(k,393)*y(k,245)
         mat(k,1319) = .080_r8*rxt(k,393)*y(k,129)
         mat(k,1237) = -(rxt(k,398)*y(k,236) + rxt(k,399)*y(k,237) + rxt(k,400) &
                      *y(k,242) + rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131))
         mat(k,1425) = -rxt(k,398)*y(k,103)
         mat(k,1986) = -rxt(k,399)*y(k,103)
         mat(k,2285) = -rxt(k,400)*y(k,103)
         mat(k,1933) = -rxt(k,401)*y(k,103)
         mat(k,2066) = -rxt(k,402)*y(k,103)
         mat(k,883) = rxt(k,397)*y(k,131)
         mat(k,2066) = mat(k,2066) + rxt(k,397)*y(k,100)
         mat(k,424) = -(rxt(k,403)*y(k,253))
         mat(k,1724) = -rxt(k,403)*y(k,104)
         mat(k,1228) = rxt(k,400)*y(k,242)
         mat(k,2229) = rxt(k,400)*y(k,103)
         mat(k,88) = -(rxt(k,521)*y(k,242) + rxt(k,522)*y(k,129))
         mat(k,2208) = -rxt(k,521)*y(k,105)
         mat(k,1870) = -rxt(k,522)*y(k,105)
         mat(k,878) = rxt(k,524)*y(k,253)
         mat(k,1670) = rxt(k,524)*y(k,100)
         mat(k,780) = -(rxt(k,406)*y(k,253))
         mat(k,1768) = -rxt(k,406)*y(k,106)
         mat(k,2258) = rxt(k,386)*y(k,244) + rxt(k,391)*y(k,245)
         mat(k,1290) = rxt(k,386)*y(k,242)
         mat(k,1321) = rxt(k,391)*y(k,242)
         mat(k,71) = -(rxt(k,527)*y(k,253))
         mat(k,1668) = -rxt(k,527)*y(k,107)
         mat(k,69) = -(rxt(k,525)*y(k,242) + rxt(k,526)*y(k,129))
         mat(k,2201) = -rxt(k,525)*y(k,108)
         mat(k,1863) = -rxt(k,526)*y(k,108)
         mat(k,70) = rxt(k,527)*y(k,253)
         mat(k,1667) = rxt(k,527)*y(k,107)
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
         mat(k,107) = -(rxt(k,530)*y(k,253))
         mat(k,1678) = -rxt(k,530)*y(k,109)
         mat(k,105) = -(rxt(k,528)*y(k,242) + rxt(k,529)*y(k,129))
         mat(k,2215) = -rxt(k,528)*y(k,110)
         mat(k,1877) = -rxt(k,529)*y(k,110)
         mat(k,106) = rxt(k,530)*y(k,253)
         mat(k,1677) = rxt(k,530)*y(k,109)
         mat(k,1253) = -(rxt(k,357)*y(k,140) + rxt(k,358)*y(k,253))
         mat(k,2168) = -rxt(k,357)*y(k,111)
         mat(k,1805) = -rxt(k,358)*y(k,111)
         mat(k,884) = .300_r8*rxt(k,404)*y(k,140)
         mat(k,1934) = .360_r8*rxt(k,387)*y(k,244)
         mat(k,2067) = .400_r8*rxt(k,388)*y(k,244)
         mat(k,2168) = mat(k,2168) + .300_r8*rxt(k,404)*y(k,100)
         mat(k,1426) = .390_r8*rxt(k,384)*y(k,244)
         mat(k,1987) = .310_r8*rxt(k,385)*y(k,244)
         mat(k,1299) = .360_r8*rxt(k,387)*y(k,129) + .400_r8*rxt(k,388)*y(k,131) &
                      + .390_r8*rxt(k,384)*y(k,236) + .310_r8*rxt(k,385)*y(k,237)
         mat(k,384) = -(rxt(k,359)*y(k,253))
         mat(k,1718) = -rxt(k,359)*y(k,112)
         mat(k,2226) = rxt(k,353)*y(k,246)
         mat(k,1350) = rxt(k,353)*y(k,242)
         mat(k,551) = -(rxt(k,368)*y(k,253))
         mat(k,1742) = -rxt(k,368)*y(k,113)
         mat(k,1895) = .800_r8*rxt(k,377)*y(k,230)
         mat(k,955) = .800_r8*rxt(k,377)*y(k,129)
         mat(k,389) = -(rxt(k,369)*y(k,253))
         mat(k,1719) = -rxt(k,369)*y(k,114)
         mat(k,2227) = .800_r8*rxt(k,366)*y(k,250)
         mat(k,698) = .800_r8*rxt(k,366)*y(k,242)
         mat(k,653) = -(rxt(k,370)*y(k,253))
         mat(k,1755) = -rxt(k,370)*y(k,115)
         mat(k,1555) = rxt(k,373)*y(k,248)
         mat(k,1395) = rxt(k,373)*y(k,130)
         mat(k,936) = -(rxt(k,461)*y(k,131) + rxt(k,462)*y(k,140) + rxt(k,463) &
                      *y(k,253))
         mat(k,2046) = -rxt(k,461)*y(k,116)
         mat(k,2150) = -rxt(k,462)*y(k,116)
         mat(k,1783) = -rxt(k,463)*y(k,116)
         mat(k,94) = -(rxt(k,532)*y(k,242) + rxt(k,533)*y(k,129))
         mat(k,2209) = -rxt(k,532)*y(k,117)
         mat(k,1871) = -rxt(k,533)*y(k,117)
         mat(k,932) = rxt(k,535)*y(k,253)
         mat(k,1671) = rxt(k,535)*y(k,116)
         mat(k,1379) = -(rxt(k,371)*y(k,140) + rxt(k,372)*y(k,253))
         mat(k,2174) = -rxt(k,371)*y(k,118)
         mat(k,1811) = -rxt(k,372)*y(k,118)
         mat(k,887) = .200_r8*rxt(k,404)*y(k,140)
         mat(k,1939) = .560_r8*rxt(k,387)*y(k,244)
         mat(k,2073) = .600_r8*rxt(k,388)*y(k,244)
         mat(k,2174) = mat(k,2174) + .200_r8*rxt(k,404)*y(k,100)
         mat(k,1431) = .610_r8*rxt(k,384)*y(k,244)
         mat(k,1992) = .440_r8*rxt(k,385)*y(k,244)
         mat(k,1303) = .560_r8*rxt(k,387)*y(k,129) + .600_r8*rxt(k,388)*y(k,131) &
                      + .610_r8*rxt(k,384)*y(k,236) + .440_r8*rxt(k,385)*y(k,237)
         mat(k,527) = -(rxt(k,170)*y(k,129) + (rxt(k,171) + rxt(k,172) + rxt(k,173) &
                      ) * y(k,130) + rxt(k,182)*y(k,253))
         mat(k,1892) = -rxt(k,170)*y(k,119)
         mat(k,1551) = -(rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,119)
         mat(k,1739) = -rxt(k,182)*y(k,119)
         mat(k,224) = -((rxt(k,186) + rxt(k,187)) * y(k,252))
         mat(k,1630) = -(rxt(k,186) + rxt(k,187)) * y(k,120)
         mat(k,526) = rxt(k,171)*y(k,130)
         mat(k,1547) = rxt(k,171)*y(k,119)
         mat(k,1548) = rxt(k,189)*y(k,131)
         mat(k,2038) = rxt(k,189)*y(k,130)
         mat(k,466) = -(rxt(k,407)*y(k,253))
         mat(k,1730) = -rxt(k,407)*y(k,122)
         mat(k,1229) = .200_r8*rxt(k,399)*y(k,237)
         mat(k,1964) = .200_r8*rxt(k,399)*y(k,103)
         mat(k,1069) = -(rxt(k,408)*y(k,253))
         mat(k,1792) = -rxt(k,408)*y(k,123)
         mat(k,1233) = rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131) + rxt(k,398)*y(k,236) &
                      + .800_r8*rxt(k,399)*y(k,237)
         mat(k,1921) = rxt(k,401)*y(k,103)
         mat(k,2054) = rxt(k,402)*y(k,103)
         mat(k,1420) = rxt(k,398)*y(k,103)
         mat(k,1975) = .800_r8*rxt(k,399)*y(k,103)
         mat(k,138) = -(rxt(k,498)*y(k,253))
         mat(k,1684) = -rxt(k,498)*y(k,127)
         mat(k,1952) = -(rxt(k,170)*y(k,119) + rxt(k,179)*y(k,131) + rxt(k,183) &
                      *y(k,242) + rxt(k,184)*y(k,140) + rxt(k,185)*y(k,139) + rxt(k,206) &
                      *y(k,61) + rxt(k,238)*y(k,21) + rxt(k,281)*y(k,237) + rxt(k,289) &
                      *y(k,243) + rxt(k,302)*y(k,233) + rxt(k,313)*y(k,236) + rxt(k,317) &
                      *y(k,241) + rxt(k,330)*y(k,234) + rxt(k,338)*y(k,255) + rxt(k,342) &
                      *y(k,256) + (rxt(k,348) + rxt(k,349)) * y(k,239) + (rxt(k,355) &
                      + rxt(k,356)) * y(k,246) + rxt(k,364)*y(k,248) + rxt(k,367) &
                      *y(k,250) + (rxt(k,377) + rxt(k,378)) * y(k,230) + rxt(k,387) &
                      *y(k,244) + rxt(k,393)*y(k,245) + rxt(k,401)*y(k,103) + rxt(k,412) &
                      *y(k,260) + rxt(k,416)*y(k,229) + rxt(k,419)*y(k,231) + rxt(k,424) &
                      *y(k,232) + rxt(k,426)*y(k,235) + rxt(k,430)*y(k,238) + rxt(k,433) &
                      *y(k,247) + rxt(k,436)*y(k,249) + rxt(k,439)*y(k,254) + rxt(k,446) &
                      *y(k,259) + rxt(k,452)*y(k,261) + rxt(k,455)*y(k,262) + rxt(k,466) &
                      *y(k,251) + rxt(k,471)*y(k,257) + rxt(k,476)*y(k,258))
         mat(k,533) = -rxt(k,170)*y(k,129)
         mat(k,2086) = -rxt(k,179)*y(k,129)
         mat(k,2305) = -rxt(k,183)*y(k,129)
         mat(k,2187) = -rxt(k,184)*y(k,129)
         mat(k,1612) = -rxt(k,185)*y(k,129)
         mat(k,2331) = -rxt(k,206)*y(k,129)
         mat(k,1849) = -rxt(k,238)*y(k,129)
         mat(k,2004) = -rxt(k,281)*y(k,129)
         mat(k,496) = -rxt(k,289)*y(k,129)
         mat(k,870) = -rxt(k,302)*y(k,129)
         mat(k,1440) = -rxt(k,313)*y(k,129)
         mat(k,756) = -rxt(k,317)*y(k,129)
         mat(k,842) = -rxt(k,330)*y(k,129)
         mat(k,819) = -rxt(k,338)*y(k,129)
         mat(k,1206) = -rxt(k,342)*y(k,129)
         mat(k,615) = -(rxt(k,348) + rxt(k,349)) * y(k,129)
         mat(k,1366) = -(rxt(k,355) + rxt(k,356)) * y(k,129)
         mat(k,1408) = -rxt(k,364)*y(k,129)
         mat(k,704) = -rxt(k,367)*y(k,129)
         mat(k,967) = -(rxt(k,377) + rxt(k,378)) * y(k,129)
         mat(k,1311) = -rxt(k,387)*y(k,129)
         mat(k,1344) = -rxt(k,393)*y(k,129)
         mat(k,1247) = -rxt(k,401)*y(k,129)
         mat(k,1224) = -rxt(k,412)*y(k,129)
         mat(k,567) = -rxt(k,416)*y(k,129)
         mat(k,541) = -rxt(k,419)*y(k,129)
         mat(k,490) = -rxt(k,424)*y(k,129)
         mat(k,673) = -rxt(k,426)*y(k,129)
         mat(k,809) = -rxt(k,430)*y(k,129)
         mat(k,762) = -rxt(k,433)*y(k,129)
         mat(k,924) = -rxt(k,436)*y(k,129)
         mat(k,503) = -rxt(k,439)*y(k,129)
         mat(k,777) = -rxt(k,446)*y(k,129)
         mat(k,802) = -rxt(k,452)*y(k,129)
         mat(k,549) = -rxt(k,455)*y(k,129)
         mat(k,1110) = -rxt(k,466)*y(k,129)
         mat(k,1171) = -rxt(k,471)*y(k,129)
         mat(k,1192) = -rxt(k,476)*y(k,129)
         mat(k,533) = mat(k,533) + 2.000_r8*rxt(k,172)*y(k,130) + rxt(k,182)*y(k,253)
         mat(k,226) = 2.000_r8*rxt(k,186)*y(k,252)
         mat(k,1581) = 2.000_r8*rxt(k,172)*y(k,119) + rxt(k,175)*y(k,139) + rxt(k,491) &
                      *y(k,158)
         mat(k,1612) = mat(k,1612) + rxt(k,175)*y(k,130)
         mat(k,1283) = rxt(k,491)*y(k,130)
         mat(k,1653) = 2.000_r8*rxt(k,186)*y(k,120)
         mat(k,1825) = rxt(k,182)*y(k,119)
         mat(k,1576) = -((rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,119) + (rxt(k,175) &
                      + rxt(k,177)) * y(k,139) + rxt(k,176)*y(k,140) + rxt(k,188) &
                      *y(k,242) + rxt(k,189)*y(k,131) + rxt(k,190)*y(k,253) + rxt(k,208) &
                      *y(k,61) + rxt(k,239)*y(k,21) + rxt(k,324)*y(k,236) + rxt(k,373) &
                      *y(k,248) + rxt(k,431)*y(k,238) + rxt(k,434)*y(k,247) + rxt(k,437) &
                      *y(k,249) + rxt(k,441)*y(k,147) + rxt(k,444)*y(k,229) + rxt(k,491) &
                      *y(k,158))
         mat(k,529) = -(rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,130)
         mat(k,1607) = -(rxt(k,175) + rxt(k,177)) * y(k,130)
         mat(k,2182) = -rxt(k,176)*y(k,130)
         mat(k,2300) = -rxt(k,188)*y(k,130)
         mat(k,2081) = -rxt(k,189)*y(k,130)
         mat(k,1820) = -rxt(k,190)*y(k,130)
         mat(k,2326) = -rxt(k,208)*y(k,130)
         mat(k,1844) = -rxt(k,239)*y(k,130)
         mat(k,1437) = -rxt(k,324)*y(k,130)
         mat(k,1405) = -rxt(k,373)*y(k,130)
         mat(k,807) = -rxt(k,431)*y(k,130)
         mat(k,761) = -rxt(k,434)*y(k,130)
         mat(k,922) = -rxt(k,437)*y(k,130)
         mat(k,524) = -rxt(k,441)*y(k,130)
         mat(k,565) = -rxt(k,444)*y(k,130)
         mat(k,1279) = -rxt(k,491)*y(k,130)
         mat(k,695) = rxt(k,375)*y(k,253)
         mat(k,399) = rxt(k,346)*y(k,131)
         mat(k,1844) = mat(k,1844) + rxt(k,238)*y(k,129)
         mat(k,2326) = mat(k,2326) + rxt(k,206)*y(k,129)
         mat(k,443) = rxt(k,169)*y(k,253)
         mat(k,623) = .700_r8*rxt(k,395)*y(k,253)
         mat(k,1244) = rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131)
         mat(k,1947) = rxt(k,238)*y(k,21) + rxt(k,206)*y(k,61) + rxt(k,401)*y(k,103) &
                      + 2.000_r8*rxt(k,179)*y(k,131) + rxt(k,185)*y(k,139) &
                      + rxt(k,184)*y(k,140) + rxt(k,416)*y(k,229) + rxt(k,377) &
                      *y(k,230) + rxt(k,419)*y(k,231) + rxt(k,424)*y(k,232) &
                      + rxt(k,302)*y(k,233) + rxt(k,330)*y(k,234) + rxt(k,426) &
                      *y(k,235) + rxt(k,313)*y(k,236) + rxt(k,281)*y(k,237) &
                      + rxt(k,430)*y(k,238) + rxt(k,348)*y(k,239) + rxt(k,317) &
                      *y(k,241) + rxt(k,183)*y(k,242) + rxt(k,289)*y(k,243) &
                      + .920_r8*rxt(k,387)*y(k,244) + .920_r8*rxt(k,393)*y(k,245) &
                      + rxt(k,355)*y(k,246) + rxt(k,433)*y(k,247) + rxt(k,364) &
                      *y(k,248) + rxt(k,436)*y(k,249) + rxt(k,367)*y(k,250) &
                      + 1.600_r8*rxt(k,466)*y(k,251) + rxt(k,439)*y(k,254) &
                      + rxt(k,338)*y(k,255) + rxt(k,342)*y(k,256) + .900_r8*rxt(k,471) &
                      *y(k,257) + .800_r8*rxt(k,476)*y(k,258) + rxt(k,446)*y(k,259) &
                      + rxt(k,412)*y(k,260) + rxt(k,452)*y(k,261) + rxt(k,455) &
                      *y(k,262)
         mat(k,2081) = mat(k,2081) + rxt(k,346)*y(k,18) + rxt(k,402)*y(k,103) &
                      + 2.000_r8*rxt(k,179)*y(k,129) + rxt(k,180)*y(k,139) &
                      + rxt(k,178)*y(k,242) + rxt(k,388)*y(k,244) + rxt(k,394) &
                      *y(k,245) + rxt(k,354)*y(k,246) + rxt(k,365)*y(k,248) &
                      + 2.000_r8*rxt(k,467)*y(k,251) + rxt(k,181)*y(k,253) &
                      + rxt(k,413)*y(k,260)
         mat(k,899) = rxt(k,336)*y(k,253)
         mat(k,1607) = mat(k,1607) + rxt(k,185)*y(k,129) + rxt(k,180)*y(k,131)
         mat(k,2182) = mat(k,2182) + rxt(k,184)*y(k,129)
         mat(k,665) = rxt(k,473)*y(k,253)
         mat(k,565) = mat(k,565) + rxt(k,416)*y(k,129)
         mat(k,965) = rxt(k,377)*y(k,129)
         mat(k,539) = rxt(k,419)*y(k,129)
         mat(k,488) = rxt(k,424)*y(k,129)
         mat(k,868) = rxt(k,302)*y(k,129)
         mat(k,840) = rxt(k,330)*y(k,129)
         mat(k,671) = rxt(k,426)*y(k,129)
         mat(k,1437) = mat(k,1437) + rxt(k,313)*y(k,129)
         mat(k,1999) = rxt(k,281)*y(k,129) + .500_r8*rxt(k,464)*y(k,251)
         mat(k,807) = mat(k,807) + rxt(k,430)*y(k,129)
         mat(k,614) = rxt(k,348)*y(k,129)
         mat(k,754) = rxt(k,317)*y(k,129)
         mat(k,2300) = mat(k,2300) + rxt(k,183)*y(k,129) + rxt(k,178)*y(k,131)
         mat(k,495) = rxt(k,289)*y(k,129)
         mat(k,1308) = .920_r8*rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131)
         mat(k,1341) = .920_r8*rxt(k,393)*y(k,129) + rxt(k,394)*y(k,131)
         mat(k,1363) = rxt(k,355)*y(k,129) + rxt(k,354)*y(k,131)
         mat(k,761) = mat(k,761) + rxt(k,433)*y(k,129)
         mat(k,1405) = mat(k,1405) + rxt(k,364)*y(k,129) + rxt(k,365)*y(k,131)
         mat(k,922) = mat(k,922) + rxt(k,436)*y(k,129)
         mat(k,702) = rxt(k,367)*y(k,129)
         mat(k,1108) = 1.600_r8*rxt(k,466)*y(k,129) + 2.000_r8*rxt(k,467)*y(k,131) &
                      + .500_r8*rxt(k,464)*y(k,237)
         mat(k,1820) = mat(k,1820) + rxt(k,375)*y(k,1) + rxt(k,169)*y(k,92) &
                      + .700_r8*rxt(k,395)*y(k,101) + rxt(k,181)*y(k,131) + rxt(k,336) &
                      *y(k,132) + rxt(k,473)*y(k,214)
         mat(k,501) = rxt(k,439)*y(k,129)
         mat(k,817) = rxt(k,338)*y(k,129)
         mat(k,1204) = rxt(k,342)*y(k,129)
         mat(k,1168) = .900_r8*rxt(k,471)*y(k,129)
         mat(k,1189) = .800_r8*rxt(k,476)*y(k,129)
         mat(k,775) = rxt(k,446)*y(k,129)
         mat(k,1221) = rxt(k,412)*y(k,129) + rxt(k,413)*y(k,131)
         mat(k,800) = rxt(k,452)*y(k,129)
         mat(k,547) = rxt(k,455)*y(k,129)
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
         mat(k,2089) = -(rxt(k,178)*y(k,242) + rxt(k,179)*y(k,129) + rxt(k,180) &
                      *y(k,139) + rxt(k,181)*y(k,253) + rxt(k,189)*y(k,130) + rxt(k,275) &
                      *y(k,44) + rxt(k,307)*y(k,47) + rxt(k,326)*y(k,31) + rxt(k,333) &
                      *y(k,51) + rxt(k,346)*y(k,18) + rxt(k,354)*y(k,246) + rxt(k,365) &
                      *y(k,248) + rxt(k,388)*y(k,244) + rxt(k,394)*y(k,245) + rxt(k,397) &
                      *y(k,100) + rxt(k,402)*y(k,103) + rxt(k,413)*y(k,260) + rxt(k,458) &
                      *y(k,6) + rxt(k,461)*y(k,116) + rxt(k,467)*y(k,251) + rxt(k,478) &
                      *y(k,216) + rxt(k,481)*y(k,69))
         mat(k,2308) = -rxt(k,178)*y(k,131)
         mat(k,1955) = -rxt(k,179)*y(k,131)
         mat(k,1615) = -rxt(k,180)*y(k,131)
         mat(k,1828) = -rxt(k,181)*y(k,131)
         mat(k,1584) = -rxt(k,189)*y(k,131)
         mat(k,1541) = -rxt(k,275)*y(k,131)
         mat(k,1087) = -rxt(k,307)*y(k,131)
         mat(k,1061) = -rxt(k,326)*y(k,131)
         mat(k,1270) = -rxt(k,333)*y(k,131)
         mat(k,401) = -rxt(k,346)*y(k,131)
         mat(k,1368) = -rxt(k,354)*y(k,131)
         mat(k,1410) = -rxt(k,365)*y(k,131)
         mat(k,1313) = -rxt(k,388)*y(k,131)
         mat(k,1346) = -rxt(k,394)*y(k,131)
         mat(k,892) = -rxt(k,397)*y(k,131)
         mat(k,1249) = -rxt(k,402)*y(k,131)
         mat(k,1226) = -rxt(k,413)*y(k,131)
         mat(k,1031) = -rxt(k,458)*y(k,131)
         mat(k,952) = -rxt(k,461)*y(k,131)
         mat(k,1112) = -rxt(k,467)*y(k,131)
         mat(k,1041) = -rxt(k,478)*y(k,131)
         mat(k,353) = -rxt(k,481)*y(k,131)
         mat(k,584) = rxt(k,240)*y(k,139)
         mat(k,2128) = rxt(k,207)*y(k,62)
         mat(k,978) = rxt(k,207)*y(k,58) + rxt(k,209)*y(k,139) + rxt(k,210)*y(k,253)
         mat(k,916) = rxt(k,254)*y(k,91)
         mat(k,1503) = rxt(k,254)*y(k,75) + rxt(k,191)*y(k,253)
         mat(k,660) = .500_r8*rxt(k,370)*y(k,253)
         mat(k,1584) = mat(k,1584) + rxt(k,177)*y(k,139) + rxt(k,176)*y(k,140)
         mat(k,1615) = mat(k,1615) + rxt(k,240)*y(k,22) + rxt(k,209)*y(k,62) &
                      + rxt(k,177)*y(k,130)
         mat(k,2190) = rxt(k,176)*y(k,130)
         mat(k,576) = rxt(k,322)*y(k,253)
         mat(k,1828) = mat(k,1828) + rxt(k,210)*y(k,62) + rxt(k,191)*y(k,91) &
                      + .500_r8*rxt(k,370)*y(k,115) + rxt(k,322)*y(k,145)
         mat(k,895) = -(rxt(k,336)*y(k,253))
         mat(k,1779) = -rxt(k,336)*y(k,132)
         mat(k,1046) = rxt(k,326)*y(k,131)
         mat(k,594) = .500_r8*rxt(k,396)*y(k,253)
         mat(k,426) = rxt(k,403)*y(k,253)
         mat(k,467) = rxt(k,407)*y(k,253)
         mat(k,1066) = rxt(k,408)*y(k,253)
         mat(k,2044) = rxt(k,326)*y(k,31)
         mat(k,1779) = mat(k,1779) + .500_r8*rxt(k,396)*y(k,102) + rxt(k,403)*y(k,104) &
                      + rxt(k,407)*y(k,122) + rxt(k,408)*y(k,123)
         mat(k,448) = -(rxt(k,468)*y(k,253))
         mat(k,1727) = -rxt(k,468)*y(k,133)
         mat(k,2232) = rxt(k,465)*y(k,251)
         mat(k,1098) = rxt(k,465)*y(k,242)
         mat(k,1608) = -(rxt(k,149)*y(k,140) + 4._r8*rxt(k,150)*y(k,139) + rxt(k,152) &
                      *y(k,79) + rxt(k,153)*y(k,81) + rxt(k,158)*y(k,242) + rxt(k,164) &
                      *y(k,253) + (rxt(k,175) + rxt(k,177)) * y(k,130) + rxt(k,180) &
                      *y(k,131) + rxt(k,185)*y(k,129) + rxt(k,209)*y(k,62) + rxt(k,211) &
                      *y(k,61) + rxt(k,214)*y(k,87) + rxt(k,217)*y(k,94) + rxt(k,240) &
                      *y(k,22) + rxt(k,241)*y(k,21) + rxt(k,243)*y(k,83) + rxt(k,245) &
                      *y(k,93) + rxt(k,276)*y(k,44) + rxt(k,483)*y(k,143))
         mat(k,2183) = -rxt(k,149)*y(k,139)
         mat(k,1452) = -rxt(k,152)*y(k,139)
         mat(k,647) = -rxt(k,153)*y(k,139)
         mat(k,2301) = -rxt(k,158)*y(k,139)
         mat(k,1821) = -rxt(k,164)*y(k,139)
         mat(k,1577) = -(rxt(k,175) + rxt(k,177)) * y(k,139)
         mat(k,2082) = -rxt(k,180)*y(k,139)
         mat(k,1948) = -rxt(k,185)*y(k,139)
         mat(k,975) = -rxt(k,209)*y(k,139)
         mat(k,2327) = -rxt(k,211)*y(k,139)
         mat(k,2023) = -rxt(k,214)*y(k,139)
         mat(k,856) = -rxt(k,217)*y(k,139)
         mat(k,582) = -rxt(k,240)*y(k,139)
         mat(k,1845) = -rxt(k,241)*y(k,139)
         mat(k,849) = -rxt(k,243)*y(k,139)
         mat(k,825) = -rxt(k,245)*y(k,139)
         mat(k,1534) = -rxt(k,276)*y(k,139)
         mat(k,416) = -rxt(k,483)*y(k,139)
         mat(k,1513) = rxt(k,156)*y(k,242)
         mat(k,530) = rxt(k,170)*y(k,129) + rxt(k,171)*y(k,130)
         mat(k,1948) = mat(k,1948) + rxt(k,170)*y(k,119)
         mat(k,1577) = mat(k,1577) + rxt(k,171)*y(k,119)
         mat(k,2301) = mat(k,2301) + rxt(k,156)*y(k,78)
         mat(k,1821) = mat(k,1821) + 2.000_r8*rxt(k,166)*y(k,253)
         mat(k,2192) = -(rxt(k,148)*y(k,252) + rxt(k,149)*y(k,139) + rxt(k,159) &
                      *y(k,242) + rxt(k,160)*y(k,78) + rxt(k,165)*y(k,253) + rxt(k,176) &
                      *y(k,130) + rxt(k,184)*y(k,129) + rxt(k,200)*y(k,58) + rxt(k,232) &
                      *y(k,19) + rxt(k,298)*y(k,27) + rxt(k,327)*y(k,31) + rxt(k,357) &
                      *y(k,111) + rxt(k,371)*y(k,118) + rxt(k,404)*y(k,100) + rxt(k,442) &
                      *y(k,147) + rxt(k,459)*y(k,6) + rxt(k,462)*y(k,116) + rxt(k,487) &
                      *y(k,156) + rxt(k,493)*y(k,158))
         mat(k,1658) = -rxt(k,148)*y(k,140)
         mat(k,1617) = -rxt(k,149)*y(k,140)
         mat(k,2310) = -rxt(k,159)*y(k,140)
         mat(k,1520) = -rxt(k,160)*y(k,140)
         mat(k,1830) = -rxt(k,165)*y(k,140)
         mat(k,1586) = -rxt(k,176)*y(k,140)
         mat(k,1957) = -rxt(k,184)*y(k,140)
         mat(k,2130) = -rxt(k,200)*y(k,140)
         mat(k,1469) = -rxt(k,232)*y(k,140)
         mat(k,607) = -rxt(k,298)*y(k,140)
         mat(k,1062) = -rxt(k,327)*y(k,140)
         mat(k,1261) = -rxt(k,357)*y(k,140)
         mat(k,1391) = -rxt(k,371)*y(k,140)
         mat(k,893) = -rxt(k,404)*y(k,140)
         mat(k,525) = -rxt(k,442)*y(k,140)
         mat(k,1032) = -rxt(k,459)*y(k,140)
         mat(k,953) = -rxt(k,462)*y(k,140)
         mat(k,561) = -rxt(k,487)*y(k,140)
         mat(k,1285) = -rxt(k,493)*y(k,140)
         mat(k,1443) = .150_r8*rxt(k,312)*y(k,242)
         mat(k,2310) = mat(k,2310) + .150_r8*rxt(k,312)*y(k,236) + .150_r8*rxt(k,362) &
                      *y(k,248)
         mat(k,1411) = .150_r8*rxt(k,362)*y(k,242)
         mat(k,360) = -(rxt(k,494)*y(k,158))
         mat(k,1273) = -rxt(k,494)*y(k,142)
         mat(k,1835) = rxt(k,234)*y(k,61)
         mat(k,2317) = rxt(k,234)*y(k,21) + 2.000_r8*rxt(k,204)*y(k,61)
         mat(k,410) = -(rxt(k,483)*y(k,139) + rxt(k,484)*y(k,253))
         mat(k,1590) = -rxt(k,483)*y(k,143)
         mat(k,1722) = -rxt(k,484)*y(k,143)
         mat(k,1131) = rxt(k,350)*y(k,253)
         mat(k,1881) = .100_r8*rxt(k,471)*y(k,257)
         mat(k,1700) = rxt(k,350)*y(k,95)
         mat(k,1155) = .100_r8*rxt(k,471)*y(k,129)
         mat(k,569) = -(rxt(k,322)*y(k,253))
         mat(k,1745) = -rxt(k,322)*y(k,145)
         mat(k,1553) = rxt(k,324)*y(k,236)
         mat(k,1415) = rxt(k,324)*y(k,130)
         mat(k,1546) = rxt(k,444)*y(k,229)
         mat(k,562) = rxt(k,444)*y(k,130)
         mat(k,522) = -(rxt(k,441)*y(k,130) + rxt(k,442)*y(k,140))
         mat(k,1550) = -rxt(k,441)*y(k,147)
         mat(k,2141) = -rxt(k,442)*y(k,147)
         mat(k,241) = .070_r8*rxt(k,428)*y(k,253)
         mat(k,1891) = rxt(k,426)*y(k,235)
         mat(k,219) = .060_r8*rxt(k,440)*y(k,253)
         mat(k,262) = .070_r8*rxt(k,456)*y(k,253)
         mat(k,669) = rxt(k,426)*y(k,129)
         mat(k,1738) = .070_r8*rxt(k,428)*y(k,68) + .060_r8*rxt(k,440)*y(k,148) &
                      + .070_r8*rxt(k,456)*y(k,225)
         mat(k,217) = -(rxt(k,440)*y(k,253))
         mat(k,1691) = -rxt(k,440)*y(k,148)
         mat(k,209) = .530_r8*rxt(k,417)*y(k,253)
         mat(k,1691) = mat(k,1691) + .530_r8*rxt(k,417)*y(k,8)
         mat(k,365) = -(rxt(k,443)*y(k,253))
         mat(k,1714) = -rxt(k,443)*y(k,149)
         mat(k,2223) = rxt(k,438)*y(k,254)
         mat(k,498) = rxt(k,438)*y(k,242)
         mat(k,585) = -(rxt(k,339)*y(k,253))
         mat(k,1746) = -rxt(k,339)*y(k,154)
         mat(k,2245) = rxt(k,337)*y(k,255)
         mat(k,812) = rxt(k,337)*y(k,242)
         mat(k,418) = -(rxt(k,343)*y(k,253))
         mat(k,1723) = -rxt(k,343)*y(k,155)
         mat(k,2228) = .850_r8*rxt(k,341)*y(k,256)
         mat(k,1197) = .850_r8*rxt(k,341)*y(k,242)
         mat(k,556) = -(rxt(k,487)*y(k,140) + rxt(k,490)*y(k,253))
         mat(k,2142) = -rxt(k,487)*y(k,156)
         mat(k,1743) = -rxt(k,490)*y(k,156)
         mat(k,1276) = -(rxt(k,488)*y(k,21) + rxt(k,489)*y(k,61) + rxt(k,491)*y(k,130) &
                      + rxt(k,493)*y(k,140) + rxt(k,494)*y(k,142) + rxt(k,495) &
                      *y(k,253))
         mat(k,1839) = -rxt(k,488)*y(k,158)
         mat(k,2321) = -rxt(k,489)*y(k,158)
         mat(k,1568) = -rxt(k,491)*y(k,158)
         mat(k,2170) = -rxt(k,493)*y(k,158)
         mat(k,362) = -rxt(k,494)*y(k,158)
         mat(k,1807) = -rxt(k,495)*y(k,158)
         mat(k,1600) = rxt(k,483)*y(k,143)
         mat(k,2170) = mat(k,2170) + rxt(k,487)*y(k,156)
         mat(k,414) = rxt(k,483)*y(k,139)
         mat(k,557) = rxt(k,487)*y(k,140) + rxt(k,490)*y(k,253)
         mat(k,1807) = mat(k,1807) + rxt(k,490)*y(k,156)
         mat(k,902) = -(rxt(k,486)*y(k,253))
         mat(k,1780) = -rxt(k,486)*y(k,159)
         mat(k,1838) = rxt(k,488)*y(k,158)
         mat(k,2319) = rxt(k,489)*y(k,158)
         mat(k,350) = rxt(k,481)*y(k,131) + (rxt(k,482)+.500_r8*rxt(k,496))*y(k,253)
         mat(k,1561) = rxt(k,491)*y(k,158)
         mat(k,2045) = rxt(k,481)*y(k,69)
         mat(k,2149) = rxt(k,493)*y(k,158)
         mat(k,361) = rxt(k,494)*y(k,158)
         mat(k,412) = rxt(k,484)*y(k,253)
         mat(k,1275) = rxt(k,488)*y(k,21) + rxt(k,489)*y(k,61) + rxt(k,491)*y(k,130) &
                      + rxt(k,493)*y(k,140) + rxt(k,494)*y(k,142) + rxt(k,495) &
                      *y(k,253)
         mat(k,1780) = mat(k,1780) + (rxt(k,482)+.500_r8*rxt(k,496))*y(k,69) &
                      + rxt(k,484)*y(k,143) + rxt(k,495)*y(k,158)
         mat(k,302) = -(rxt(k,497)*y(k,263))
         mat(k,2341) = -rxt(k,497)*y(k,160)
         mat(k,901) = rxt(k,486)*y(k,253)
         mat(k,1707) = rxt(k,486)*y(k,159)
         mat(k,64) = .1056005_r8*rxt(k,526)*y(k,129) + .2381005_r8*rxt(k,525)*y(k,242)
         mat(k,1858) = .1056005_r8*rxt(k,526)*y(k,108)
         mat(k,114) = .5931005_r8*rxt(k,536)*y(k,253)
         mat(k,2196) = .2381005_r8*rxt(k,525)*y(k,108)
         mat(k,1662) = .5931005_r8*rxt(k,536)*y(k,210)
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
         mat(k,65) = .1026005_r8*rxt(k,526)*y(k,129) + .1308005_r8*rxt(k,525)*y(k,242)
         mat(k,1859) = .1026005_r8*rxt(k,526)*y(k,108)
         mat(k,115) = .1534005_r8*rxt(k,536)*y(k,253)
         mat(k,2197) = .1308005_r8*rxt(k,525)*y(k,108)
         mat(k,1663) = .1534005_r8*rxt(k,536)*y(k,210)
         mat(k,66) = .0521005_r8*rxt(k,526)*y(k,129) + .0348005_r8*rxt(k,525)*y(k,242)
         mat(k,1860) = .0521005_r8*rxt(k,526)*y(k,108)
         mat(k,116) = .0459005_r8*rxt(k,536)*y(k,253)
         mat(k,2198) = .0348005_r8*rxt(k,525)*y(k,108)
         mat(k,1664) = .0459005_r8*rxt(k,536)*y(k,210)
         mat(k,67) = .0143005_r8*rxt(k,526)*y(k,129) + .0076005_r8*rxt(k,525)*y(k,242)
         mat(k,1861) = .0143005_r8*rxt(k,526)*y(k,108)
         mat(k,117) = .0085005_r8*rxt(k,536)*y(k,253)
         mat(k,2199) = .0076005_r8*rxt(k,525)*y(k,108)
         mat(k,1665) = .0085005_r8*rxt(k,536)*y(k,210)
         mat(k,68) = .0166005_r8*rxt(k,526)*y(k,129) + .0113005_r8*rxt(k,525)*y(k,242)
         mat(k,1862) = .0166005_r8*rxt(k,526)*y(k,108)
         mat(k,118) = .0128005_r8*rxt(k,536)*y(k,253)
         mat(k,2200) = .0113005_r8*rxt(k,525)*y(k,108)
         mat(k,1666) = .0128005_r8*rxt(k,536)*y(k,210)
         mat(k,1002) = .2202005_r8*rxt(k,515)*y(k,140)
         mat(k,77) = .1279005_r8*rxt(k,514)*y(k,129) + .2202005_r8*rxt(k,513)*y(k,242)
         mat(k,83) = .0003005_r8*rxt(k,522)*y(k,129) + .0031005_r8*rxt(k,521)*y(k,242)
         mat(k,927) = .0508005_r8*rxt(k,534)*y(k,140)
         mat(k,89) = .0245005_r8*rxt(k,533)*y(k,129) + .0508005_r8*rxt(k,532)*y(k,242)
         mat(k,1864) = .1279005_r8*rxt(k,514)*y(k,7) + .0003005_r8*rxt(k,522)*y(k,105) &
                      + .0245005_r8*rxt(k,533)*y(k,117)
         mat(k,2134) = .2202005_r8*rxt(k,515)*y(k,6) + .0508005_r8*rxt(k,534)*y(k,116)
         mat(k,2202) = .2202005_r8*rxt(k,513)*y(k,7) + .0031005_r8*rxt(k,521)*y(k,105) &
                      + .0508005_r8*rxt(k,532)*y(k,117)
         mat(k,1003) = .2067005_r8*rxt(k,515)*y(k,140)
         mat(k,78) = .1792005_r8*rxt(k,514)*y(k,129) + .2067005_r8*rxt(k,513)*y(k,242)
         mat(k,84) = .0003005_r8*rxt(k,522)*y(k,129) + .0035005_r8*rxt(k,521)*y(k,242)
         mat(k,928) = .1149005_r8*rxt(k,534)*y(k,140)
         mat(k,90) = .0082005_r8*rxt(k,533)*y(k,129) + .1149005_r8*rxt(k,532)*y(k,242)
         mat(k,1865) = .1792005_r8*rxt(k,514)*y(k,7) + .0003005_r8*rxt(k,522)*y(k,105) &
                      + .0082005_r8*rxt(k,533)*y(k,117)
         mat(k,2135) = .2067005_r8*rxt(k,515)*y(k,6) + .1149005_r8*rxt(k,534)*y(k,116)
         mat(k,2203) = .2067005_r8*rxt(k,513)*y(k,7) + .0035005_r8*rxt(k,521)*y(k,105) &
                      + .1149005_r8*rxt(k,532)*y(k,117)
         mat(k,1004) = .0653005_r8*rxt(k,515)*y(k,140)
         mat(k,79) = .0676005_r8*rxt(k,514)*y(k,129) + .0653005_r8*rxt(k,513)*y(k,242)
         mat(k,85) = .0073005_r8*rxt(k,522)*y(k,129) + .0003005_r8*rxt(k,521)*y(k,242)
         mat(k,929) = .0348005_r8*rxt(k,534)*y(k,140)
         mat(k,91) = .0772005_r8*rxt(k,533)*y(k,129) + .0348005_r8*rxt(k,532)*y(k,242)
         mat(k,1866) = .0676005_r8*rxt(k,514)*y(k,7) + .0073005_r8*rxt(k,522)*y(k,105) &
                      + .0772005_r8*rxt(k,533)*y(k,117)
         mat(k,2136) = .0653005_r8*rxt(k,515)*y(k,6) + .0348005_r8*rxt(k,534)*y(k,116)
         mat(k,2204) = .0653005_r8*rxt(k,513)*y(k,7) + .0003005_r8*rxt(k,521)*y(k,105) &
                      + .0348005_r8*rxt(k,532)*y(k,117)
         mat(k,1005) = .1749305_r8*rxt(k,512)*y(k,131) + .1284005_r8*rxt(k,515) &
                      *y(k,140)
         mat(k,80) = .079_r8*rxt(k,514)*y(k,129) + .1284005_r8*rxt(k,513)*y(k,242)
         mat(k,876) = .0590245_r8*rxt(k,520)*y(k,131) + .0033005_r8*rxt(k,523) &
                      *y(k,140)
         mat(k,86) = .0057005_r8*rxt(k,522)*y(k,129) + .0271005_r8*rxt(k,521)*y(k,242)
         mat(k,930) = .1749305_r8*rxt(k,531)*y(k,131) + .0554005_r8*rxt(k,534) &
                      *y(k,140)
         mat(k,92) = .0332005_r8*rxt(k,533)*y(k,129) + .0554005_r8*rxt(k,532)*y(k,242)
         mat(k,1867) = .079_r8*rxt(k,514)*y(k,7) + .0057005_r8*rxt(k,522)*y(k,105) &
                      + .0332005_r8*rxt(k,533)*y(k,117)
         mat(k,2036) = .1749305_r8*rxt(k,512)*y(k,6) + .0590245_r8*rxt(k,520)*y(k,100) &
                      + .1749305_r8*rxt(k,531)*y(k,116)
         mat(k,2137) = .1284005_r8*rxt(k,515)*y(k,6) + .0033005_r8*rxt(k,523)*y(k,100) &
                      + .0554005_r8*rxt(k,534)*y(k,116)
         mat(k,2205) = .1284005_r8*rxt(k,513)*y(k,7) + .0271005_r8*rxt(k,521)*y(k,105) &
                      + .0554005_r8*rxt(k,532)*y(k,117)
         mat(k,1006) = .5901905_r8*rxt(k,512)*y(k,131) + .114_r8*rxt(k,515)*y(k,140)
         mat(k,81) = .1254005_r8*rxt(k,514)*y(k,129) + .114_r8*rxt(k,513)*y(k,242)
         mat(k,877) = .0250245_r8*rxt(k,520)*y(k,131)
         mat(k,87) = .0623005_r8*rxt(k,522)*y(k,129) + .0474005_r8*rxt(k,521)*y(k,242)
         mat(k,931) = .5901905_r8*rxt(k,531)*y(k,131) + .1278005_r8*rxt(k,534) &
                      *y(k,140)
         mat(k,93) = .130_r8*rxt(k,533)*y(k,129) + .1278005_r8*rxt(k,532)*y(k,242)
         mat(k,1868) = .1254005_r8*rxt(k,514)*y(k,7) + .0623005_r8*rxt(k,522)*y(k,105) &
                      + .130_r8*rxt(k,533)*y(k,117)
         mat(k,2037) = .5901905_r8*rxt(k,512)*y(k,6) + .0250245_r8*rxt(k,520)*y(k,100) &
                      + .5901905_r8*rxt(k,531)*y(k,116)
         mat(k,2138) = .114_r8*rxt(k,515)*y(k,6) + .1278005_r8*rxt(k,534)*y(k,116)
         mat(k,2206) = .114_r8*rxt(k,513)*y(k,7) + .0474005_r8*rxt(k,521)*y(k,105) &
                      + .1278005_r8*rxt(k,532)*y(k,117)
         mat(k,108) = .0097005_r8*rxt(k,519)*y(k,129) + .0023005_r8*rxt(k,518) &
                      *y(k,242)
         mat(k,100) = .1056005_r8*rxt(k,529)*y(k,129) + .2381005_r8*rxt(k,528) &
                      *y(k,242)
         mat(k,1872) = .0097005_r8*rxt(k,519)*y(k,9) + .1056005_r8*rxt(k,529)*y(k,110) &
                      + .0154005_r8*rxt(k,540)*y(k,220) + .0063005_r8*rxt(k,544) &
                      *y(k,224)
         mat(k,120) = .5931005_r8*rxt(k,537)*y(k,253)
         mat(k,126) = .0154005_r8*rxt(k,540)*y(k,129) + .1364005_r8*rxt(k,539) &
                      *y(k,242)
         mat(k,132) = .0063005_r8*rxt(k,544)*y(k,129) + .1677005_r8*rxt(k,543) &
                      *y(k,242)
         mat(k,2210) = .0023005_r8*rxt(k,518)*y(k,9) + .2381005_r8*rxt(k,528)*y(k,110) &
                      + .1364005_r8*rxt(k,539)*y(k,220) + .1677005_r8*rxt(k,543) &
                      *y(k,224)
         mat(k,1672) = .5931005_r8*rxt(k,537)*y(k,211)
         mat(k,109) = .0034005_r8*rxt(k,519)*y(k,129) + .0008005_r8*rxt(k,518) &
                      *y(k,242)
         mat(k,101) = .1026005_r8*rxt(k,529)*y(k,129) + .1308005_r8*rxt(k,528) &
                      *y(k,242)
         mat(k,1873) = .0034005_r8*rxt(k,519)*y(k,9) + .1026005_r8*rxt(k,529)*y(k,110) &
                      + .0452005_r8*rxt(k,540)*y(k,220) + .0237005_r8*rxt(k,544) &
                      *y(k,224)
         mat(k,121) = .1534005_r8*rxt(k,537)*y(k,253)
         mat(k,127) = .0452005_r8*rxt(k,540)*y(k,129) + .0101005_r8*rxt(k,539) &
                      *y(k,242)
         mat(k,133) = .0237005_r8*rxt(k,544)*y(k,129) + .0174005_r8*rxt(k,543) &
                      *y(k,242)
         mat(k,2211) = .0008005_r8*rxt(k,518)*y(k,9) + .1308005_r8*rxt(k,528)*y(k,110) &
                      + .0101005_r8*rxt(k,539)*y(k,220) + .0174005_r8*rxt(k,543) &
                      *y(k,224)
         mat(k,1673) = .1534005_r8*rxt(k,537)*y(k,211)
         mat(k,110) = .1579005_r8*rxt(k,519)*y(k,129) + .0843005_r8*rxt(k,518) &
                      *y(k,242)
         mat(k,102) = .0521005_r8*rxt(k,529)*y(k,129) + .0348005_r8*rxt(k,528) &
                      *y(k,242)
         mat(k,1874) = .1579005_r8*rxt(k,519)*y(k,9) + .0521005_r8*rxt(k,529)*y(k,110) &
                      + .0966005_r8*rxt(k,540)*y(k,220) + .0025005_r8*rxt(k,544) &
                      *y(k,224)
         mat(k,122) = .0459005_r8*rxt(k,537)*y(k,253)
         mat(k,128) = .0966005_r8*rxt(k,540)*y(k,129) + .0763005_r8*rxt(k,539) &
                      *y(k,242)
         mat(k,134) = .0025005_r8*rxt(k,544)*y(k,129) + .086_r8*rxt(k,543)*y(k,242)
         mat(k,2212) = .0843005_r8*rxt(k,518)*y(k,9) + .0348005_r8*rxt(k,528)*y(k,110) &
                      + .0763005_r8*rxt(k,539)*y(k,220) + .086_r8*rxt(k,543)*y(k,224)
         mat(k,1674) = .0459005_r8*rxt(k,537)*y(k,211)
         mat(k,111) = .0059005_r8*rxt(k,519)*y(k,129) + .0443005_r8*rxt(k,518) &
                      *y(k,242)
         mat(k,103) = .0143005_r8*rxt(k,529)*y(k,129) + .0076005_r8*rxt(k,528) &
                      *y(k,242)
         mat(k,1875) = .0059005_r8*rxt(k,519)*y(k,9) + .0143005_r8*rxt(k,529)*y(k,110) &
                      + .0073005_r8*rxt(k,540)*y(k,220) + .011_r8*rxt(k,544)*y(k,224)
         mat(k,123) = .0085005_r8*rxt(k,537)*y(k,253)
         mat(k,129) = .0073005_r8*rxt(k,540)*y(k,129) + .2157005_r8*rxt(k,539) &
                      *y(k,242)
         mat(k,135) = .011_r8*rxt(k,544)*y(k,129) + .0512005_r8*rxt(k,543)*y(k,242)
         mat(k,2213) = .0443005_r8*rxt(k,518)*y(k,9) + .0076005_r8*rxt(k,528)*y(k,110) &
                      + .2157005_r8*rxt(k,539)*y(k,220) + .0512005_r8*rxt(k,543) &
                      *y(k,224)
         mat(k,1675) = .0085005_r8*rxt(k,537)*y(k,211)
         mat(k,112) = .0536005_r8*rxt(k,519)*y(k,129) + .1621005_r8*rxt(k,518) &
                      *y(k,242)
         mat(k,104) = .0166005_r8*rxt(k,529)*y(k,129) + .0113005_r8*rxt(k,528) &
                      *y(k,242)
         mat(k,1876) = .0536005_r8*rxt(k,519)*y(k,9) + .0166005_r8*rxt(k,529)*y(k,110) &
                      + .238_r8*rxt(k,540)*y(k,220) + .1185005_r8*rxt(k,544)*y(k,224)
         mat(k,124) = .0128005_r8*rxt(k,537)*y(k,253)
         mat(k,130) = .238_r8*rxt(k,540)*y(k,129) + .0738005_r8*rxt(k,539)*y(k,242)
         mat(k,136) = .1185005_r8*rxt(k,544)*y(k,129) + .1598005_r8*rxt(k,543) &
                      *y(k,242)
         mat(k,2214) = .1621005_r8*rxt(k,518)*y(k,9) + .0113005_r8*rxt(k,528)*y(k,110) &
                      + .0738005_r8*rxt(k,539)*y(k,220) + .1598005_r8*rxt(k,543) &
                      *y(k,224)
         mat(k,1676) = .0128005_r8*rxt(k,537)*y(k,211)
         mat(k,119) = -(rxt(k,536)*y(k,253))
         mat(k,1680) = -rxt(k,536)*y(k,210)
         mat(k,125) = -(rxt(k,537)*y(k,253))
         mat(k,1681) = -rxt(k,537)*y(k,211)
         mat(k,234) = .100_r8*rxt(k,448)*y(k,253)
         mat(k,252) = .230_r8*rxt(k,450)*y(k,253)
         mat(k,1695) = .100_r8*rxt(k,448)*y(k,219) + .230_r8*rxt(k,450)*y(k,222)
         mat(k,706) = -(rxt(k,472)*y(k,253))
         mat(k,1761) = -rxt(k,472)*y(k,213)
         mat(k,2251) = rxt(k,470)*y(k,257)
         mat(k,1156) = rxt(k,470)*y(k,242)
         mat(k,662) = -(rxt(k,473)*y(k,253))
         mat(k,1756) = -rxt(k,473)*y(k,214)
         mat(k,1901) = .200_r8*rxt(k,466)*y(k,251) + .200_r8*rxt(k,476)*y(k,258)
         mat(k,1965) = .500_r8*rxt(k,464)*y(k,251)
         mat(k,1099) = .200_r8*rxt(k,466)*y(k,129) + .500_r8*rxt(k,464)*y(k,237)
         mat(k,1176) = .200_r8*rxt(k,476)*y(k,129)
         mat(k,515) = -(rxt(k,477)*y(k,253))
         mat(k,1737) = -rxt(k,477)*y(k,215)
         mat(k,2240) = rxt(k,475)*y(k,258)
         mat(k,1175) = rxt(k,475)*y(k,242)
         mat(k,1035) = -(rxt(k,478)*y(k,131) + rxt(k,479)*y(k,253))
         mat(k,2052) = -rxt(k,478)*y(k,216)
         mat(k,1790) = -rxt(k,479)*y(k,216)
         mat(k,1016) = .330_r8*rxt(k,459)*y(k,140)
         mat(k,941) = .330_r8*rxt(k,462)*y(k,140)
         mat(k,1920) = .800_r8*rxt(k,466)*y(k,251) + .800_r8*rxt(k,476)*y(k,258)
         mat(k,2052) = mat(k,2052) + rxt(k,467)*y(k,251)
         mat(k,2156) = .330_r8*rxt(k,459)*y(k,6) + .330_r8*rxt(k,462)*y(k,116)
         mat(k,663) = rxt(k,473)*y(k,253)
         mat(k,1974) = .500_r8*rxt(k,464)*y(k,251) + rxt(k,474)*y(k,258)
         mat(k,1101) = .800_r8*rxt(k,466)*y(k,129) + rxt(k,467)*y(k,131) &
                      + .500_r8*rxt(k,464)*y(k,237)
         mat(k,1790) = mat(k,1790) + rxt(k,473)*y(k,214)
         mat(k,1179) = .800_r8*rxt(k,476)*y(k,129) + rxt(k,474)*y(k,237)
         mat(k,1116) = -(rxt(k,480)*y(k,253))
         mat(k,1796) = -rxt(k,480)*y(k,217)
         mat(k,1019) = .300_r8*rxt(k,459)*y(k,140)
         mat(k,943) = .300_r8*rxt(k,462)*y(k,140)
         mat(k,1925) = .900_r8*rxt(k,471)*y(k,257)
         mat(k,2161) = .300_r8*rxt(k,459)*y(k,6) + .300_r8*rxt(k,462)*y(k,116)
         mat(k,1979) = rxt(k,469)*y(k,257)
         mat(k,1161) = .900_r8*rxt(k,471)*y(k,129) + rxt(k,469)*y(k,237)
         mat(k,719) = -(rxt(k,447)*y(k,253))
         mat(k,1762) = -rxt(k,447)*y(k,218)
         mat(k,2252) = rxt(k,445)*y(k,259)
         mat(k,767) = rxt(k,445)*y(k,242)
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
         mat(k,232) = -(rxt(k,448)*y(k,253))
         mat(k,1693) = -rxt(k,448)*y(k,219)
         mat(k,131) = -(rxt(k,539)*y(k,242) + rxt(k,540)*y(k,129))
         mat(k,2217) = -rxt(k,539)*y(k,220)
         mat(k,1879) = -rxt(k,540)*y(k,220)
         mat(k,231) = rxt(k,538)*y(k,253)
         mat(k,1682) = rxt(k,538)*y(k,219)
         mat(k,248) = -(rxt(k,414)*y(k,253))
         mat(k,1696) = -rxt(k,414)*y(k,221)
         mat(k,2220) = rxt(k,411)*y(k,260)
         mat(k,1210) = rxt(k,411)*y(k,242)
         mat(k,253) = -(rxt(k,450)*y(k,253))
         mat(k,1697) = -rxt(k,450)*y(k,222)
         mat(k,738) = -(rxt(k,453)*y(k,253))
         mat(k,1764) = -rxt(k,453)*y(k,223)
         mat(k,2254) = rxt(k,451)*y(k,261)
         mat(k,791) = rxt(k,451)*y(k,242)
         mat(k,137) = -(rxt(k,543)*y(k,242) + rxt(k,544)*y(k,129))
         mat(k,2218) = -rxt(k,543)*y(k,224)
         mat(k,1880) = -rxt(k,544)*y(k,224)
         mat(k,251) = rxt(k,542)*y(k,253)
         mat(k,1683) = rxt(k,542)*y(k,222)
         mat(k,261) = -(rxt(k,456)*y(k,253))
         mat(k,1698) = -rxt(k,456)*y(k,225)
         mat(k,254) = .150_r8*rxt(k,450)*y(k,253)
         mat(k,1698) = mat(k,1698) + .150_r8*rxt(k,450)*y(k,222)
         mat(k,478) = -(rxt(k,457)*y(k,253))
         mat(k,1732) = -rxt(k,457)*y(k,226)
         mat(k,2235) = rxt(k,454)*y(k,262)
         mat(k,543) = rxt(k,454)*y(k,242)
         mat(k,563) = -(rxt(k,415)*y(k,242) + rxt(k,416)*y(k,129) + rxt(k,444) &
                      *y(k,130))
         mat(k,2244) = -rxt(k,415)*y(k,229)
         mat(k,1896) = -rxt(k,416)*y(k,229)
         mat(k,1552) = -rxt(k,444)*y(k,229)
         mat(k,276) = rxt(k,421)*y(k,253)
         mat(k,1744) = rxt(k,421)*y(k,24)
         mat(k,960) = -(rxt(k,376)*y(k,242) + (rxt(k,377) + rxt(k,378)) * y(k,129))
         mat(k,2270) = -rxt(k,376)*y(k,230)
         mat(k,1916) = -(rxt(k,377) + rxt(k,378)) * y(k,230)
         mat(k,680) = rxt(k,379)*y(k,253)
         mat(k,267) = rxt(k,380)*y(k,253)
         mat(k,1784) = rxt(k,379)*y(k,2) + rxt(k,380)*y(k,17)
         mat(k,536) = -(rxt(k,418)*y(k,242) + rxt(k,419)*y(k,129))
         mat(k,2242) = -rxt(k,418)*y(k,231)
         mat(k,1893) = -rxt(k,419)*y(k,231)
         mat(k,210) = .350_r8*rxt(k,417)*y(k,253)
         mat(k,474) = rxt(k,420)*y(k,253)
         mat(k,1740) = .350_r8*rxt(k,417)*y(k,8) + rxt(k,420)*y(k,10)
         mat(k,486) = -(rxt(k,422)*y(k,242) + rxt(k,424)*y(k,129))
         mat(k,2236) = -rxt(k,422)*y(k,232)
         mat(k,1887) = -rxt(k,424)*y(k,232)
         mat(k,372) = rxt(k,423)*y(k,253)
         mat(k,235) = .070_r8*rxt(k,448)*y(k,253)
         mat(k,255) = .060_r8*rxt(k,450)*y(k,253)
         mat(k,1733) = rxt(k,423)*y(k,25) + .070_r8*rxt(k,448)*y(k,219) &
                      + .060_r8*rxt(k,450)*y(k,222)
         mat(k,865) = -(4._r8*rxt(k,299)*y(k,233) + rxt(k,300)*y(k,237) + rxt(k,301) &
                      *y(k,242) + rxt(k,302)*y(k,129))
         mat(k,1969) = -rxt(k,300)*y(k,233)
         mat(k,2266) = -rxt(k,301)*y(k,233)
         mat(k,1912) = -rxt(k,302)*y(k,233)
         mat(k,377) = .500_r8*rxt(k,304)*y(k,253)
         mat(k,328) = rxt(k,305)*y(k,58) + rxt(k,306)*y(k,253)
         mat(k,2106) = rxt(k,305)*y(k,30)
         mat(k,1776) = .500_r8*rxt(k,304)*y(k,29) + rxt(k,306)*y(k,30)
         mat(k,836) = -(rxt(k,328)*y(k,237) + rxt(k,329)*y(k,242) + rxt(k,330) &
                      *y(k,129))
         mat(k,1968) = -rxt(k,328)*y(k,234)
         mat(k,2263) = -rxt(k,329)*y(k,234)
         mat(k,1911) = -rxt(k,330)*y(k,234)
         mat(k,437) = rxt(k,331)*y(k,253)
         mat(k,155) = rxt(k,332)*y(k,253)
         mat(k,1773) = rxt(k,331)*y(k,32) + rxt(k,332)*y(k,33)
         mat(k,670) = -(rxt(k,425)*y(k,242) + rxt(k,426)*y(k,129))
         mat(k,2248) = -rxt(k,425)*y(k,235)
         mat(k,1902) = -rxt(k,426)*y(k,235)
         mat(k,312) = rxt(k,427)*y(k,253)
         mat(k,1902) = mat(k,1902) + rxt(k,416)*y(k,229)
         mat(k,2144) = rxt(k,442)*y(k,147)
         mat(k,523) = rxt(k,442)*y(k,140)
         mat(k,564) = rxt(k,416)*y(k,129) + .400_r8*rxt(k,415)*y(k,242)
         mat(k,2248) = mat(k,2248) + .400_r8*rxt(k,415)*y(k,229)
         mat(k,1757) = rxt(k,427)*y(k,34)
         mat(k,1433) = -(4._r8*rxt(k,310)*y(k,236) + rxt(k,311)*y(k,237) + rxt(k,312) &
                      *y(k,242) + rxt(k,313)*y(k,129) + rxt(k,324)*y(k,130) + rxt(k,351) &
                      *y(k,246) + rxt(k,384)*y(k,244) + rxt(k,389)*y(k,245) + rxt(k,398) &
                      *y(k,103) + rxt(k,409)*y(k,260))
         mat(k,1994) = -rxt(k,311)*y(k,236)
         mat(k,2293) = -rxt(k,312)*y(k,236)
         mat(k,1941) = -rxt(k,313)*y(k,236)
         mat(k,1570) = -rxt(k,324)*y(k,236)
         mat(k,1360) = -rxt(k,351)*y(k,236)
         mat(k,1305) = -rxt(k,384)*y(k,236)
         mat(k,1338) = -rxt(k,389)*y(k,236)
         mat(k,1241) = -rxt(k,398)*y(k,236)
         mat(k,1219) = -rxt(k,409)*y(k,236)
         mat(k,1024) = .060_r8*rxt(k,459)*y(k,140)
         mat(k,1083) = rxt(k,307)*y(k,131) + rxt(k,308)*y(k,253)
         mat(k,1266) = rxt(k,333)*y(k,131) + rxt(k,334)*y(k,253)
         mat(k,637) = .500_r8*rxt(k,315)*y(k,253)
         mat(k,888) = .080_r8*rxt(k,404)*y(k,140)
         mat(k,1257) = .100_r8*rxt(k,357)*y(k,140)
         mat(k,948) = .060_r8*rxt(k,462)*y(k,140)
         mat(k,1381) = .280_r8*rxt(k,371)*y(k,140)
         mat(k,1941) = mat(k,1941) + .530_r8*rxt(k,355)*y(k,246) + rxt(k,364)*y(k,248) &
                      + rxt(k,367)*y(k,250) + rxt(k,342)*y(k,256)
         mat(k,2075) = rxt(k,307)*y(k,47) + rxt(k,333)*y(k,51) + .530_r8*rxt(k,354) &
                      *y(k,246) + rxt(k,365)*y(k,248)
         mat(k,2176) = .060_r8*rxt(k,459)*y(k,6) + .080_r8*rxt(k,404)*y(k,100) &
                      + .100_r8*rxt(k,357)*y(k,111) + .060_r8*rxt(k,462)*y(k,116) &
                      + .280_r8*rxt(k,371)*y(k,118)
         mat(k,1119) = .650_r8*rxt(k,480)*y(k,253)
         mat(k,1433) = mat(k,1433) + .530_r8*rxt(k,351)*y(k,246)
         mat(k,1994) = mat(k,1994) + .260_r8*rxt(k,352)*y(k,246) + rxt(k,361)*y(k,248) &
                      + .300_r8*rxt(k,340)*y(k,256)
         mat(k,2293) = mat(k,2293) + .450_r8*rxt(k,362)*y(k,248) + .200_r8*rxt(k,366) &
                      *y(k,250) + .150_r8*rxt(k,341)*y(k,256)
         mat(k,1360) = mat(k,1360) + .530_r8*rxt(k,355)*y(k,129) + .530_r8*rxt(k,354) &
                      *y(k,131) + .530_r8*rxt(k,351)*y(k,236) + .260_r8*rxt(k,352) &
                      *y(k,237)
         mat(k,1402) = rxt(k,364)*y(k,129) + rxt(k,365)*y(k,131) + rxt(k,361)*y(k,237) &
                      + .450_r8*rxt(k,362)*y(k,242) + 4.000_r8*rxt(k,363)*y(k,248)
         mat(k,701) = rxt(k,367)*y(k,129) + .200_r8*rxt(k,366)*y(k,242)
         mat(k,1813) = rxt(k,308)*y(k,47) + rxt(k,334)*y(k,51) + .500_r8*rxt(k,315) &
                      *y(k,53) + .650_r8*rxt(k,480)*y(k,217)
         mat(k,1202) = rxt(k,342)*y(k,129) + .300_r8*rxt(k,340)*y(k,237) &
                      + .150_r8*rxt(k,341)*y(k,242)
         mat(k,2005) = -(rxt(k,201)*y(k,61) + (4._r8*rxt(k,278) + 4._r8*rxt(k,279) &
                      ) * y(k,237) + rxt(k,280)*y(k,242) + rxt(k,281)*y(k,129) &
                      + rxt(k,300)*y(k,233) + rxt(k,311)*y(k,236) + rxt(k,328) &
                      *y(k,234) + rxt(k,340)*y(k,256) + rxt(k,352)*y(k,246) + rxt(k,361) &
                      *y(k,248) + rxt(k,385)*y(k,244) + rxt(k,390)*y(k,245) + rxt(k,399) &
                      *y(k,103) + rxt(k,410)*y(k,260) + rxt(k,464)*y(k,251) + rxt(k,469) &
                      *y(k,257) + rxt(k,474)*y(k,258))
         mat(k,2332) = -rxt(k,201)*y(k,237)
         mat(k,2306) = -rxt(k,280)*y(k,237)
         mat(k,1953) = -rxt(k,281)*y(k,237)
         mat(k,871) = -rxt(k,300)*y(k,237)
         mat(k,1441) = -rxt(k,311)*y(k,237)
         mat(k,843) = -rxt(k,328)*y(k,237)
         mat(k,1207) = -rxt(k,340)*y(k,237)
         mat(k,1367) = -rxt(k,352)*y(k,237)
         mat(k,1409) = -rxt(k,361)*y(k,237)
         mat(k,1312) = -rxt(k,385)*y(k,237)
         mat(k,1345) = -rxt(k,390)*y(k,237)
         mat(k,1248) = -rxt(k,399)*y(k,237)
         mat(k,1225) = -rxt(k,410)*y(k,237)
         mat(k,1111) = -rxt(k,464)*y(k,237)
         mat(k,1172) = -rxt(k,469)*y(k,237)
         mat(k,1193) = -rxt(k,474)*y(k,237)
         mat(k,1060) = .280_r8*rxt(k,327)*y(k,140)
         mat(k,732) = rxt(k,314)*y(k,253)
         mat(k,458) = .700_r8*rxt(k,283)*y(k,253)
         mat(k,1484) = rxt(k,195)*y(k,58) + rxt(k,251)*y(k,75) + rxt(k,290)*y(k,252) &
                      + rxt(k,284)*y(k,253)
         mat(k,2126) = rxt(k,195)*y(k,56)
         mat(k,915) = rxt(k,251)*y(k,56)
         mat(k,891) = .050_r8*rxt(k,404)*y(k,140)
         mat(k,1248) = mat(k,1248) + rxt(k,398)*y(k,236)
         mat(k,1953) = mat(k,1953) + rxt(k,313)*y(k,236) + .830_r8*rxt(k,430)*y(k,238) &
                      + .170_r8*rxt(k,436)*y(k,249)
         mat(k,2188) = .280_r8*rxt(k,327)*y(k,31) + .050_r8*rxt(k,404)*y(k,100)
         mat(k,1441) = mat(k,1441) + rxt(k,398)*y(k,103) + rxt(k,313)*y(k,129) &
                      + 4.000_r8*rxt(k,310)*y(k,236) + .900_r8*rxt(k,311)*y(k,237) &
                      + .450_r8*rxt(k,312)*y(k,242) + rxt(k,384)*y(k,244) + rxt(k,389) &
                      *y(k,245) + rxt(k,351)*y(k,246) + rxt(k,360)*y(k,248) &
                      + rxt(k,409)*y(k,260)
         mat(k,2005) = mat(k,2005) + .900_r8*rxt(k,311)*y(k,236)
         mat(k,810) = .830_r8*rxt(k,430)*y(k,129) + .330_r8*rxt(k,429)*y(k,242)
         mat(k,2306) = mat(k,2306) + .450_r8*rxt(k,312)*y(k,236) + .330_r8*rxt(k,429) &
                      *y(k,238) + .070_r8*rxt(k,435)*y(k,249)
         mat(k,1312) = mat(k,1312) + rxt(k,384)*y(k,236)
         mat(k,1345) = mat(k,1345) + rxt(k,389)*y(k,236)
         mat(k,1367) = mat(k,1367) + rxt(k,351)*y(k,236)
         mat(k,1409) = mat(k,1409) + rxt(k,360)*y(k,236)
         mat(k,925) = .170_r8*rxt(k,436)*y(k,129) + .070_r8*rxt(k,435)*y(k,242)
         mat(k,1654) = rxt(k,290)*y(k,56)
         mat(k,1826) = rxt(k,314)*y(k,52) + .700_r8*rxt(k,283)*y(k,55) + rxt(k,284) &
                      *y(k,56)
         mat(k,1225) = mat(k,1225) + rxt(k,409)*y(k,236)
         mat(k,804) = -(rxt(k,429)*y(k,242) + rxt(k,430)*y(k,129) + rxt(k,431) &
                      *y(k,130))
         mat(k,2260) = -rxt(k,429)*y(k,238)
         mat(k,1909) = -rxt(k,430)*y(k,238)
         mat(k,1558) = -rxt(k,431)*y(k,238)
         mat(k,609) = -((rxt(k,348) + rxt(k,349)) * y(k,129))
         mat(k,1898) = -(rxt(k,348) + rxt(k,349)) * y(k,239)
         mat(k,395) = rxt(k,347)*y(k,253)
         mat(k,1749) = rxt(k,347)*y(k,18)
         mat(k,1883) = .750_r8*rxt(k,317)*y(k,241)
         mat(k,750) = .750_r8*rxt(k,317)*y(k,129)
         mat(k,751) = -(rxt(k,316)*y(k,242) + rxt(k,317)*y(k,129))
         mat(k,2255) = -rxt(k,316)*y(k,241)
         mat(k,1905) = -rxt(k,317)*y(k,241)
         mat(k,602) = rxt(k,323)*y(k,253)
         mat(k,1765) = rxt(k,323)*y(k,27)
         mat(k,2311) = -((rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,78) + rxt(k,158) &
                      *y(k,139) + rxt(k,159)*y(k,140) + rxt(k,163)*y(k,253) &
                      + 4._r8*rxt(k,168)*y(k,242) + rxt(k,178)*y(k,131) + rxt(k,183) &
                      *y(k,129) + rxt(k,188)*y(k,130) + (rxt(k,198) + rxt(k,199) &
                      ) * y(k,58) + rxt(k,205)*y(k,61) + rxt(k,231)*y(k,19) + rxt(k,237) &
                      *y(k,21) + rxt(k,274)*y(k,44) + rxt(k,280)*y(k,237) + rxt(k,287) &
                      *y(k,243) + rxt(k,301)*y(k,233) + rxt(k,312)*y(k,236) + rxt(k,316) &
                      *y(k,241) + rxt(k,329)*y(k,234) + rxt(k,337)*y(k,255) + rxt(k,341) &
                      *y(k,256) + rxt(k,353)*y(k,246) + rxt(k,362)*y(k,248) + rxt(k,366) &
                      *y(k,250) + rxt(k,376)*y(k,230) + rxt(k,386)*y(k,244) + rxt(k,391) &
                      *y(k,245) + rxt(k,400)*y(k,103) + rxt(k,411)*y(k,260) + rxt(k,415) &
                      *y(k,229) + rxt(k,418)*y(k,231) + rxt(k,422)*y(k,232) + rxt(k,425) &
                      *y(k,235) + rxt(k,429)*y(k,238) + rxt(k,432)*y(k,247) + rxt(k,435) &
                      *y(k,249) + rxt(k,438)*y(k,254) + rxt(k,445)*y(k,259) + rxt(k,451) &
                      *y(k,261) + rxt(k,454)*y(k,262) + rxt(k,465)*y(k,251) + rxt(k,470) &
                      *y(k,257) + rxt(k,475)*y(k,258))
         mat(k,1521) = -(rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,242)
         mat(k,1618) = -rxt(k,158)*y(k,242)
         mat(k,2193) = -rxt(k,159)*y(k,242)
         mat(k,1831) = -rxt(k,163)*y(k,242)
         mat(k,2092) = -rxt(k,178)*y(k,242)
         mat(k,1958) = -rxt(k,183)*y(k,242)
         mat(k,1587) = -rxt(k,188)*y(k,242)
         mat(k,2131) = -(rxt(k,198) + rxt(k,199)) * y(k,242)
         mat(k,2337) = -rxt(k,205)*y(k,242)
         mat(k,1470) = -rxt(k,231)*y(k,242)
         mat(k,1855) = -rxt(k,237)*y(k,242)
         mat(k,1544) = -rxt(k,274)*y(k,242)
         mat(k,2010) = -rxt(k,280)*y(k,242)
         mat(k,497) = -rxt(k,287)*y(k,242)
         mat(k,872) = -rxt(k,301)*y(k,242)
         mat(k,1444) = -rxt(k,312)*y(k,242)
         mat(k,757) = -rxt(k,316)*y(k,242)
         mat(k,844) = -rxt(k,329)*y(k,242)
         mat(k,820) = -rxt(k,337)*y(k,242)
         mat(k,1208) = -rxt(k,341)*y(k,242)
         mat(k,1369) = -rxt(k,353)*y(k,242)
         mat(k,1412) = -rxt(k,362)*y(k,242)
         mat(k,705) = -rxt(k,366)*y(k,242)
         mat(k,969) = -rxt(k,376)*y(k,242)
         mat(k,1315) = -rxt(k,386)*y(k,242)
         mat(k,1348) = -rxt(k,391)*y(k,242)
         mat(k,1250) = -rxt(k,400)*y(k,242)
         mat(k,1227) = -rxt(k,411)*y(k,242)
         mat(k,568) = -rxt(k,415)*y(k,242)
         mat(k,542) = -rxt(k,418)*y(k,242)
         mat(k,491) = -rxt(k,422)*y(k,242)
         mat(k,675) = -rxt(k,425)*y(k,242)
         mat(k,811) = -rxt(k,429)*y(k,242)
         mat(k,763) = -rxt(k,432)*y(k,242)
         mat(k,926) = -rxt(k,435)*y(k,242)
         mat(k,504) = -rxt(k,438)*y(k,242)
         mat(k,778) = -rxt(k,445)*y(k,242)
         mat(k,803) = -rxt(k,451)*y(k,242)
         mat(k,550) = -rxt(k,454)*y(k,242)
         mat(k,1113) = -rxt(k,465)*y(k,242)
         mat(k,1173) = -rxt(k,470)*y(k,242)
         mat(k,1195) = -rxt(k,475)*y(k,242)
         mat(k,1033) = .570_r8*rxt(k,459)*y(k,140)
         mat(k,212) = .650_r8*rxt(k,417)*y(k,253)
         mat(k,1470) = mat(k,1470) + rxt(k,230)*y(k,44)
         mat(k,1855) = mat(k,1855) + rxt(k,242)*y(k,253)
         mat(k,326) = .350_r8*rxt(k,296)*y(k,253)
         mat(k,608) = .130_r8*rxt(k,298)*y(k,140)
         mat(k,309) = rxt(k,303)*y(k,253)
         mat(k,1063) = .280_r8*rxt(k,327)*y(k,140)
         mat(k,1544) = mat(k,1544) + rxt(k,230)*y(k,19) + rxt(k,194)*y(k,58) &
                      + rxt(k,275)*y(k,131) + rxt(k,276)*y(k,139)
         mat(k,633) = rxt(k,259)*y(k,58) + rxt(k,260)*y(k,253)
         mat(k,408) = rxt(k,262)*y(k,58) + rxt(k,263)*y(k,253)
         mat(k,146) = rxt(k,309)*y(k,253)
         mat(k,834) = rxt(k,282)*y(k,253)
         mat(k,1488) = rxt(k,291)*y(k,252)
         mat(k,2131) = mat(k,2131) + rxt(k,194)*y(k,44) + rxt(k,259)*y(k,45) &
                      + rxt(k,262)*y(k,48) + rxt(k,197)*y(k,81)
         mat(k,2337) = mat(k,2337) + rxt(k,201)*y(k,237) + rxt(k,212)*y(k,253)
         mat(k,1130) = rxt(k,294)*y(k,253)
         mat(k,243) = .730_r8*rxt(k,428)*y(k,253)
         mat(k,354) = .500_r8*rxt(k,496)*y(k,253)
         mat(k,1097) = rxt(k,320)*y(k,253)
         mat(k,1001) = rxt(k,321)*y(k,253)
         mat(k,651) = rxt(k,197)*y(k,58) + rxt(k,153)*y(k,139) + rxt(k,162)*y(k,253)
         mat(k,230) = rxt(k,285)*y(k,253)
         mat(k,993) = rxt(k,286)*y(k,253)
         mat(k,1147) = rxt(k,350)*y(k,253)
         mat(k,1154) = rxt(k,335)*y(k,253)
         mat(k,894) = .370_r8*rxt(k,404)*y(k,140)
         mat(k,625) = .300_r8*rxt(k,395)*y(k,253)
         mat(k,600) = rxt(k,396)*y(k,253)
         mat(k,1250) = mat(k,1250) + rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131) &
                      + rxt(k,398)*y(k,236) + 1.200_r8*rxt(k,399)*y(k,237)
         mat(k,429) = rxt(k,403)*y(k,253)
         mat(k,1262) = .140_r8*rxt(k,357)*y(k,140)
         mat(k,388) = .200_r8*rxt(k,359)*y(k,253)
         mat(k,661) = .500_r8*rxt(k,370)*y(k,253)
         mat(k,954) = .570_r8*rxt(k,462)*y(k,140)
         mat(k,1392) = .280_r8*rxt(k,371)*y(k,140)
         mat(k,471) = rxt(k,407)*y(k,253)
         mat(k,1080) = rxt(k,408)*y(k,253)
         mat(k,1958) = mat(k,1958) + rxt(k,401)*y(k,103) + rxt(k,377)*y(k,230) &
                      + rxt(k,419)*y(k,231) + rxt(k,424)*y(k,232) + rxt(k,302) &
                      *y(k,233) + rxt(k,330)*y(k,234) + rxt(k,281)*y(k,237) &
                      + .170_r8*rxt(k,430)*y(k,238) + rxt(k,348)*y(k,239) &
                      + .250_r8*rxt(k,317)*y(k,241) + rxt(k,289)*y(k,243) &
                      + .920_r8*rxt(k,387)*y(k,244) + .920_r8*rxt(k,393)*y(k,245) &
                      + .470_r8*rxt(k,355)*y(k,246) + .400_r8*rxt(k,433)*y(k,247) &
                      + .830_r8*rxt(k,436)*y(k,249) + rxt(k,439)*y(k,254) + rxt(k,338) &
                      *y(k,255) + .900_r8*rxt(k,471)*y(k,257) + .800_r8*rxt(k,476) &
                      *y(k,258) + rxt(k,446)*y(k,259) + rxt(k,412)*y(k,260) &
                      + rxt(k,452)*y(k,261) + rxt(k,455)*y(k,262)
         mat(k,2092) = mat(k,2092) + rxt(k,275)*y(k,44) + rxt(k,402)*y(k,103) &
                      + rxt(k,388)*y(k,244) + rxt(k,394)*y(k,245) + .470_r8*rxt(k,354) &
                      *y(k,246) + rxt(k,181)*y(k,253) + rxt(k,413)*y(k,260)
         mat(k,1618) = mat(k,1618) + rxt(k,276)*y(k,44) + rxt(k,153)*y(k,81)
         mat(k,2193) = mat(k,2193) + .570_r8*rxt(k,459)*y(k,6) + .130_r8*rxt(k,298) &
                      *y(k,27) + .280_r8*rxt(k,327)*y(k,31) + .370_r8*rxt(k,404) &
                      *y(k,100) + .140_r8*rxt(k,357)*y(k,111) + .570_r8*rxt(k,462) &
                      *y(k,116) + .280_r8*rxt(k,371)*y(k,118) + rxt(k,165)*y(k,253)
         mat(k,221) = .800_r8*rxt(k,440)*y(k,253)
         mat(k,906) = rxt(k,486)*y(k,253)
         mat(k,1124) = .200_r8*rxt(k,480)*y(k,253)
         mat(k,238) = .280_r8*rxt(k,448)*y(k,253)
         mat(k,260) = .380_r8*rxt(k,450)*y(k,253)
         mat(k,265) = .630_r8*rxt(k,456)*y(k,253)
         mat(k,969) = mat(k,969) + rxt(k,377)*y(k,129)
         mat(k,542) = mat(k,542) + rxt(k,419)*y(k,129)
         mat(k,491) = mat(k,491) + rxt(k,424)*y(k,129)
         mat(k,872) = mat(k,872) + rxt(k,302)*y(k,129) + 2.400_r8*rxt(k,299)*y(k,233) &
                      + rxt(k,300)*y(k,237)
         mat(k,844) = mat(k,844) + rxt(k,330)*y(k,129) + rxt(k,328)*y(k,237)
         mat(k,1444) = mat(k,1444) + rxt(k,398)*y(k,103) + .900_r8*rxt(k,311)*y(k,237) &
                      + rxt(k,384)*y(k,244) + rxt(k,389)*y(k,245) + .470_r8*rxt(k,351) &
                      *y(k,246) + rxt(k,409)*y(k,260)
         mat(k,2010) = mat(k,2010) + rxt(k,201)*y(k,61) + 1.200_r8*rxt(k,399)*y(k,103) &
                      + rxt(k,281)*y(k,129) + rxt(k,300)*y(k,233) + rxt(k,328) &
                      *y(k,234) + .900_r8*rxt(k,311)*y(k,236) + 4.000_r8*rxt(k,278) &
                      *y(k,237) + rxt(k,385)*y(k,244) + rxt(k,390)*y(k,245) &
                      + .730_r8*rxt(k,352)*y(k,246) + rxt(k,361)*y(k,248) &
                      + .500_r8*rxt(k,464)*y(k,251) + .300_r8*rxt(k,340)*y(k,256) &
                      + rxt(k,469)*y(k,257) + rxt(k,474)*y(k,258) + .800_r8*rxt(k,410) &
                      *y(k,260)
         mat(k,811) = mat(k,811) + .170_r8*rxt(k,430)*y(k,129) + .070_r8*rxt(k,429) &
                      *y(k,242)
         mat(k,616) = rxt(k,348)*y(k,129)
         mat(k,757) = mat(k,757) + .250_r8*rxt(k,317)*y(k,129)
         mat(k,2311) = mat(k,2311) + .070_r8*rxt(k,429)*y(k,238) + .160_r8*rxt(k,432) &
                      *y(k,247) + .330_r8*rxt(k,435)*y(k,249)
         mat(k,497) = mat(k,497) + rxt(k,289)*y(k,129)
         mat(k,1315) = mat(k,1315) + .920_r8*rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131) &
                      + rxt(k,384)*y(k,236) + rxt(k,385)*y(k,237)
         mat(k,1348) = mat(k,1348) + .920_r8*rxt(k,393)*y(k,129) + rxt(k,394)*y(k,131) &
                      + rxt(k,389)*y(k,236) + rxt(k,390)*y(k,237)
         mat(k,1369) = mat(k,1369) + .470_r8*rxt(k,355)*y(k,129) + .470_r8*rxt(k,354) &
                      *y(k,131) + .470_r8*rxt(k,351)*y(k,236) + .730_r8*rxt(k,352) &
                      *y(k,237)
         mat(k,763) = mat(k,763) + .400_r8*rxt(k,433)*y(k,129) + .160_r8*rxt(k,432) &
                      *y(k,242)
         mat(k,1412) = mat(k,1412) + rxt(k,361)*y(k,237)
         mat(k,926) = mat(k,926) + .830_r8*rxt(k,436)*y(k,129) + .330_r8*rxt(k,435) &
                      *y(k,242)
         mat(k,1113) = mat(k,1113) + .500_r8*rxt(k,464)*y(k,237)
         mat(k,1659) = rxt(k,291)*y(k,56)
         mat(k,1831) = mat(k,1831) + .650_r8*rxt(k,417)*y(k,8) + rxt(k,242)*y(k,21) &
                      + .350_r8*rxt(k,296)*y(k,26) + rxt(k,303)*y(k,28) + rxt(k,260) &
                      *y(k,45) + rxt(k,263)*y(k,48) + rxt(k,309)*y(k,49) + rxt(k,282) &
                      *y(k,54) + rxt(k,212)*y(k,61) + rxt(k,294)*y(k,64) &
                      + .730_r8*rxt(k,428)*y(k,68) + .500_r8*rxt(k,496)*y(k,69) &
                      + rxt(k,320)*y(k,76) + rxt(k,321)*y(k,77) + rxt(k,162)*y(k,81) &
                      + rxt(k,285)*y(k,88) + rxt(k,286)*y(k,89) + rxt(k,350)*y(k,95) &
                      + rxt(k,335)*y(k,97) + .300_r8*rxt(k,395)*y(k,101) + rxt(k,396) &
                      *y(k,102) + rxt(k,403)*y(k,104) + .200_r8*rxt(k,359)*y(k,112) &
                      + .500_r8*rxt(k,370)*y(k,115) + rxt(k,407)*y(k,122) + rxt(k,408) &
                      *y(k,123) + rxt(k,181)*y(k,131) + rxt(k,165)*y(k,140) &
                      + .800_r8*rxt(k,440)*y(k,148) + rxt(k,486)*y(k,159) &
                      + .200_r8*rxt(k,480)*y(k,217) + .280_r8*rxt(k,448)*y(k,219) &
                      + .380_r8*rxt(k,450)*y(k,222) + .630_r8*rxt(k,456)*y(k,225)
         mat(k,504) = mat(k,504) + rxt(k,439)*y(k,129)
         mat(k,820) = mat(k,820) + rxt(k,338)*y(k,129)
         mat(k,1208) = mat(k,1208) + .300_r8*rxt(k,340)*y(k,237)
         mat(k,1173) = mat(k,1173) + .900_r8*rxt(k,471)*y(k,129) + rxt(k,469)*y(k,237)
         mat(k,1195) = mat(k,1195) + .800_r8*rxt(k,476)*y(k,129) + rxt(k,474)*y(k,237)
         mat(k,778) = mat(k,778) + rxt(k,446)*y(k,129)
         mat(k,1227) = mat(k,1227) + rxt(k,412)*y(k,129) + rxt(k,413)*y(k,131) &
                      + rxt(k,409)*y(k,236) + .800_r8*rxt(k,410)*y(k,237)
         mat(k,803) = mat(k,803) + rxt(k,452)*y(k,129)
         mat(k,550) = mat(k,550) + rxt(k,455)*y(k,129)
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
         mat(k,492) = -(rxt(k,287)*y(k,242) + rxt(k,289)*y(k,129))
         mat(k,2237) = -rxt(k,287)*y(k,243)
         mat(k,1888) = -rxt(k,289)*y(k,243)
         mat(k,1523) = rxt(k,274)*y(k,242)
         mat(k,2237) = mat(k,2237) + rxt(k,274)*y(k,44)
         mat(k,1301) = -(rxt(k,384)*y(k,236) + rxt(k,385)*y(k,237) + rxt(k,386) &
                      *y(k,242) + rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131))
         mat(k,1428) = -rxt(k,384)*y(k,244)
         mat(k,1989) = -rxt(k,385)*y(k,244)
         mat(k,2288) = -rxt(k,386)*y(k,244)
         mat(k,1936) = -rxt(k,387)*y(k,244)
         mat(k,2070) = -rxt(k,388)*y(k,244)
         mat(k,885) = .600_r8*rxt(k,405)*y(k,253)
         mat(k,1808) = .600_r8*rxt(k,405)*y(k,100)
         mat(k,1334) = -(rxt(k,389)*y(k,236) + rxt(k,390)*y(k,237) + rxt(k,391) &
                      *y(k,242) + rxt(k,393)*y(k,129) + rxt(k,394)*y(k,131))
         mat(k,1429) = -rxt(k,389)*y(k,245)
         mat(k,1990) = -rxt(k,390)*y(k,245)
         mat(k,2289) = -rxt(k,391)*y(k,245)
         mat(k,1937) = -rxt(k,393)*y(k,245)
         mat(k,2071) = -rxt(k,394)*y(k,245)
         mat(k,886) = .400_r8*rxt(k,405)*y(k,253)
         mat(k,1809) = .400_r8*rxt(k,405)*y(k,100)
         mat(k,1358) = -(rxt(k,351)*y(k,236) + rxt(k,352)*y(k,237) + rxt(k,353) &
                      *y(k,242) + rxt(k,354)*y(k,131) + (rxt(k,355) + rxt(k,356) &
                      ) * y(k,129))
         mat(k,1430) = -rxt(k,351)*y(k,246)
         mat(k,1991) = -rxt(k,352)*y(k,246)
         mat(k,2290) = -rxt(k,353)*y(k,246)
         mat(k,2072) = -rxt(k,354)*y(k,246)
         mat(k,1938) = -(rxt(k,355) + rxt(k,356)) * y(k,246)
         mat(k,1255) = .500_r8*rxt(k,358)*y(k,253)
         mat(k,385) = .200_r8*rxt(k,359)*y(k,253)
         mat(k,1378) = rxt(k,372)*y(k,253)
         mat(k,1810) = .500_r8*rxt(k,358)*y(k,111) + .200_r8*rxt(k,359)*y(k,112) &
                      + rxt(k,372)*y(k,118)
         mat(k,758) = -(rxt(k,432)*y(k,242) + rxt(k,433)*y(k,129) + rxt(k,434) &
                      *y(k,130))
         mat(k,2256) = -rxt(k,432)*y(k,247)
         mat(k,1906) = -rxt(k,433)*y(k,247)
         mat(k,1557) = -rxt(k,434)*y(k,247)
         mat(k,1401) = -(rxt(k,360)*y(k,236) + rxt(k,361)*y(k,237) + rxt(k,362) &
                      *y(k,242) + 4._r8*rxt(k,363)*y(k,248) + rxt(k,364)*y(k,129) &
                      + rxt(k,365)*y(k,131) + rxt(k,373)*y(k,130))
         mat(k,1432) = -rxt(k,360)*y(k,248)
         mat(k,1993) = -rxt(k,361)*y(k,248)
         mat(k,2292) = -rxt(k,362)*y(k,248)
         mat(k,1940) = -rxt(k,364)*y(k,248)
         mat(k,2074) = -rxt(k,365)*y(k,248)
         mat(k,1569) = -rxt(k,373)*y(k,248)
         mat(k,1256) = .500_r8*rxt(k,358)*y(k,253)
         mat(k,386) = .500_r8*rxt(k,359)*y(k,253)
         mat(k,1812) = .500_r8*rxt(k,358)*y(k,111) + .500_r8*rxt(k,359)*y(k,112)
         mat(k,918) = -(rxt(k,435)*y(k,242) + rxt(k,436)*y(k,129) + rxt(k,437) &
                      *y(k,130))
         mat(k,2269) = -rxt(k,435)*y(k,249)
         mat(k,1915) = -rxt(k,436)*y(k,249)
         mat(k,1562) = -rxt(k,437)*y(k,249)
         mat(k,699) = -(rxt(k,366)*y(k,242) + rxt(k,367)*y(k,129))
         mat(k,2250) = -rxt(k,366)*y(k,250)
         mat(k,1904) = -rxt(k,367)*y(k,250)
         mat(k,552) = rxt(k,368)*y(k,253)
         mat(k,390) = rxt(k,369)*y(k,253)
         mat(k,1760) = rxt(k,368)*y(k,113) + rxt(k,369)*y(k,114)
         mat(k,1102) = -(rxt(k,464)*y(k,237) + rxt(k,465)*y(k,242) + rxt(k,466) &
                      *y(k,129) + rxt(k,467)*y(k,131))
         mat(k,1978) = -rxt(k,464)*y(k,251)
         mat(k,2277) = -rxt(k,465)*y(k,251)
         mat(k,1924) = -rxt(k,466)*y(k,251)
         mat(k,2057) = -rxt(k,467)*y(k,251)
         mat(k,1018) = rxt(k,458)*y(k,131)
         mat(k,942) = rxt(k,461)*y(k,131)
         mat(k,2057) = mat(k,2057) + rxt(k,458)*y(k,6) + rxt(k,461)*y(k,116) &
                      + .500_r8*rxt(k,478)*y(k,216)
         mat(k,450) = rxt(k,468)*y(k,253)
         mat(k,1036) = .500_r8*rxt(k,478)*y(k,131)
         mat(k,1795) = rxt(k,468)*y(k,133)
         mat(k,1650) = -(rxt(k,144)*y(k,79) + rxt(k,145)*y(k,263) + rxt(k,148) &
                      *y(k,140) + (rxt(k,186) + rxt(k,187)) * y(k,120) + rxt(k,219) &
                      *y(k,35) + rxt(k,220)*y(k,36) + rxt(k,221)*y(k,38) + rxt(k,222) &
                      *y(k,39) + rxt(k,223)*y(k,40) + rxt(k,224)*y(k,41) + rxt(k,225) &
                      *y(k,42) + (rxt(k,226) + rxt(k,227)) * y(k,87) + rxt(k,246) &
                      *y(k,37) + rxt(k,247)*y(k,57) + rxt(k,248)*y(k,80) + (rxt(k,249) &
                      + rxt(k,250)) * y(k,83) + rxt(k,255)*y(k,66) + rxt(k,256) &
                      *y(k,67) + rxt(k,269)*y(k,43) + rxt(k,270)*y(k,45) + rxt(k,271) &
                      *y(k,84) + rxt(k,272)*y(k,85) + rxt(k,273)*y(k,86) + (rxt(k,290) &
                      + rxt(k,291) + rxt(k,292)) * y(k,56) + rxt(k,293)*y(k,88))
         mat(k,1453) = -rxt(k,144)*y(k,252)
         mat(k,2354) = -rxt(k,145)*y(k,252)
         mat(k,2184) = -rxt(k,148)*y(k,252)
         mat(k,225) = -(rxt(k,186) + rxt(k,187)) * y(k,252)
         mat(k,142) = -rxt(k,219)*y(k,252)
         mat(k,183) = -rxt(k,220)*y(k,252)
         mat(k,160) = -rxt(k,221)*y(k,252)
         mat(k,193) = -rxt(k,222)*y(k,252)
         mat(k,164) = -rxt(k,223)*y(k,252)
         mat(k,198) = -rxt(k,224)*y(k,252)
         mat(k,168) = -rxt(k,225)*y(k,252)
         mat(k,2024) = -(rxt(k,226) + rxt(k,227)) * y(k,252)
         mat(k,189) = -rxt(k,246)*y(k,252)
         mat(k,462) = -rxt(k,247)*y(k,252)
         mat(k,153) = -rxt(k,248)*y(k,252)
         mat(k,850) = -(rxt(k,249) + rxt(k,250)) * y(k,252)
         mat(k,285) = -rxt(k,255)*y(k,252)
         mat(k,293) = -rxt(k,256)*y(k,252)
         mat(k,510) = -rxt(k,269)*y(k,252)
         mat(k,628) = -rxt(k,270)*y(k,252)
         mat(k,288) = -rxt(k,271)*y(k,252)
         mat(k,298) = -rxt(k,272)*y(k,252)
         mat(k,336) = -rxt(k,273)*y(k,252)
         mat(k,1482) = -(rxt(k,290) + rxt(k,291) + rxt(k,292)) * y(k,252)
         mat(k,228) = -rxt(k,293)*y(k,252)
         mat(k,1823) = -(rxt(k,161)*y(k,79) + rxt(k,162)*y(k,81) + rxt(k,163)*y(k,242) &
                      + rxt(k,164)*y(k,139) + rxt(k,165)*y(k,140) + (4._r8*rxt(k,166) &
                      + 4._r8*rxt(k,167)) * y(k,253) + rxt(k,169)*y(k,92) + rxt(k,181) &
                      *y(k,131) + rxt(k,182)*y(k,119) + rxt(k,190)*y(k,130) + rxt(k,191) &
                      *y(k,91) + rxt(k,210)*y(k,62) + (rxt(k,212) + rxt(k,213) &
                      ) * y(k,61) + rxt(k,215)*y(k,87) + rxt(k,218)*y(k,94) + rxt(k,242) &
                      *y(k,21) + rxt(k,244)*y(k,83) + rxt(k,258)*y(k,43) + rxt(k,260) &
                      *y(k,45) + rxt(k,261)*y(k,46) + rxt(k,263)*y(k,48) + rxt(k,265) &
                      *y(k,57) + rxt(k,266)*y(k,84) + rxt(k,267)*y(k,85) + rxt(k,268) &
                      *y(k,86) + rxt(k,277)*y(k,44) + rxt(k,282)*y(k,54) + rxt(k,283) &
                      *y(k,55) + rxt(k,284)*y(k,56) + rxt(k,285)*y(k,88) + rxt(k,286) &
                      *y(k,89) + rxt(k,294)*y(k,64) + rxt(k,296)*y(k,26) + rxt(k,303) &
                      *y(k,28) + rxt(k,304)*y(k,29) + rxt(k,306)*y(k,30) + rxt(k,308) &
                      *y(k,47) + rxt(k,309)*y(k,49) + rxt(k,314)*y(k,52) + rxt(k,315) &
                      *y(k,53) + rxt(k,320)*y(k,76) + rxt(k,321)*y(k,77) + rxt(k,322) &
                      *y(k,145) + rxt(k,323)*y(k,27) + rxt(k,331)*y(k,32) + rxt(k,332) &
                      *y(k,33) + rxt(k,334)*y(k,51) + rxt(k,335)*y(k,97) + rxt(k,336) &
                      *y(k,132) + rxt(k,339)*y(k,154) + rxt(k,343)*y(k,155) + rxt(k,344) &
                      *y(k,31) + rxt(k,345)*y(k,50) + rxt(k,347)*y(k,18) + rxt(k,350) &
                      *y(k,95) + rxt(k,358)*y(k,111) + rxt(k,359)*y(k,112) + rxt(k,368) &
                      *y(k,113) + rxt(k,369)*y(k,114) + rxt(k,370)*y(k,115) + rxt(k,372) &
                      *y(k,118) + rxt(k,375)*y(k,1) + rxt(k,379)*y(k,2) + rxt(k,380) &
                      *y(k,17) + rxt(k,381)*y(k,96) + rxt(k,382)*y(k,98) + rxt(k,383) &
                      *y(k,99) + rxt(k,395)*y(k,101) + rxt(k,396)*y(k,102) + rxt(k,403) &
                      *y(k,104) + rxt(k,405)*y(k,100) + rxt(k,406)*y(k,106) + rxt(k,407) &
                      *y(k,122) + rxt(k,408)*y(k,123) + rxt(k,414)*y(k,221) + rxt(k,417) &
                      *y(k,8) + rxt(k,420)*y(k,10) + rxt(k,421)*y(k,24) + rxt(k,423) &
                      *y(k,25) + rxt(k,427)*y(k,34) + rxt(k,428)*y(k,68) + rxt(k,440) &
                      *y(k,148) + rxt(k,443)*y(k,149) + rxt(k,447)*y(k,218) + rxt(k,448) &
                      *y(k,219) + rxt(k,450)*y(k,222) + rxt(k,453)*y(k,223) + rxt(k,456) &
                      *y(k,225) + rxt(k,457)*y(k,226) + rxt(k,460)*y(k,6) + rxt(k,463) &
                      *y(k,116) + rxt(k,468)*y(k,133) + rxt(k,472)*y(k,213) + rxt(k,473) &
                      *y(k,214) + rxt(k,477)*y(k,215) + rxt(k,479)*y(k,216) + rxt(k,480) &
                      *y(k,217) + (rxt(k,482) + rxt(k,496)) * y(k,69) + rxt(k,484) &
                      *y(k,143) + rxt(k,486)*y(k,159) + rxt(k,490)*y(k,156) + rxt(k,495) &
                      *y(k,158) + rxt(k,498)*y(k,127))
         mat(k,1454) = -rxt(k,161)*y(k,253)
         mat(k,648) = -rxt(k,162)*y(k,253)
         mat(k,2303) = -rxt(k,163)*y(k,253)
         mat(k,1610) = -rxt(k,164)*y(k,253)
         mat(k,2185) = -rxt(k,165)*y(k,253)
         mat(k,444) = -rxt(k,169)*y(k,253)
         mat(k,2084) = -rxt(k,181)*y(k,253)
         mat(k,532) = -rxt(k,182)*y(k,253)
         mat(k,1579) = -rxt(k,190)*y(k,253)
         mat(k,1500) = -rxt(k,191)*y(k,253)
         mat(k,976) = -rxt(k,210)*y(k,253)
         mat(k,2329) = -(rxt(k,212) + rxt(k,213)) * y(k,253)
         mat(k,2025) = -rxt(k,215)*y(k,253)
         mat(k,857) = -rxt(k,218)*y(k,253)
         mat(k,1847) = -rxt(k,242)*y(k,253)
         mat(k,851) = -rxt(k,244)*y(k,253)
         mat(k,511) = -rxt(k,258)*y(k,253)
         mat(k,629) = -rxt(k,260)*y(k,253)
         mat(k,171) = -rxt(k,261)*y(k,253)
         mat(k,404) = -rxt(k,263)*y(k,253)
         mat(k,463) = -rxt(k,265)*y(k,253)
         mat(k,289) = -rxt(k,266)*y(k,253)
         mat(k,299) = -rxt(k,267)*y(k,253)
         mat(k,337) = -rxt(k,268)*y(k,253)
         mat(k,1536) = -rxt(k,277)*y(k,253)
         mat(k,833) = -rxt(k,282)*y(k,253)
         mat(k,457) = -rxt(k,283)*y(k,253)
         mat(k,1483) = -rxt(k,284)*y(k,253)
         mat(k,229) = -rxt(k,285)*y(k,253)
         mat(k,992) = -rxt(k,286)*y(k,253)
         mat(k,1129) = -rxt(k,294)*y(k,253)
         mat(k,325) = -rxt(k,296)*y(k,253)
         mat(k,308) = -rxt(k,303)*y(k,253)
         mat(k,379) = -rxt(k,304)*y(k,253)
         mat(k,329) = -rxt(k,306)*y(k,253)
         mat(k,1085) = -rxt(k,308)*y(k,253)
         mat(k,145) = -rxt(k,309)*y(k,253)
         mat(k,731) = -rxt(k,314)*y(k,253)
         mat(k,639) = -rxt(k,315)*y(k,253)
         mat(k,1096) = -rxt(k,320)*y(k,253)
         mat(k,1000) = -rxt(k,321)*y(k,253)
         mat(k,574) = -rxt(k,322)*y(k,253)
         mat(k,606) = -rxt(k,323)*y(k,253)
         mat(k,439) = -rxt(k,331)*y(k,253)
         mat(k,156) = -rxt(k,332)*y(k,253)
         mat(k,1269) = -rxt(k,334)*y(k,253)
         mat(k,1153) = -rxt(k,335)*y(k,253)
         mat(k,900) = -rxt(k,336)*y(k,253)
         mat(k,590) = -rxt(k,339)*y(k,253)
         mat(k,422) = -rxt(k,343)*y(k,253)
         mat(k,1058) = -rxt(k,344)*y(k,253)
         mat(k,985) = -rxt(k,345)*y(k,253)
         mat(k,400) = -rxt(k,347)*y(k,253)
         mat(k,1144) = -rxt(k,350)*y(k,253)
         mat(k,1260) = -rxt(k,358)*y(k,253)
         mat(k,387) = -rxt(k,359)*y(k,253)
         mat(k,555) = -rxt(k,368)*y(k,253)
         mat(k,393) = -rxt(k,369)*y(k,253)
         mat(k,659) = -rxt(k,370)*y(k,253)
         mat(k,1387) = -rxt(k,372)*y(k,253)
         mat(k,696) = -rxt(k,375)*y(k,253)
         mat(k,685) = -rxt(k,379)*y(k,253)
         mat(k,268) = -rxt(k,380)*y(k,253)
         mat(k,281) = -rxt(k,381)*y(k,253)
         mat(k,383) = -rxt(k,382)*y(k,253)
         mat(k,176) = -rxt(k,383)*y(k,253)
         mat(k,624) = -rxt(k,395)*y(k,253)
         mat(k,599) = -rxt(k,396)*y(k,253)
         mat(k,428) = -rxt(k,403)*y(k,253)
         mat(k,890) = -rxt(k,405)*y(k,253)
         mat(k,785) = -rxt(k,406)*y(k,253)
         mat(k,470) = -rxt(k,407)*y(k,253)
         mat(k,1077) = -rxt(k,408)*y(k,253)
         mat(k,250) = -rxt(k,414)*y(k,253)
         mat(k,211) = -rxt(k,417)*y(k,253)
         mat(k,476) = -rxt(k,420)*y(k,253)
         mat(k,277) = -rxt(k,421)*y(k,253)
         mat(k,374) = -rxt(k,423)*y(k,253)
         mat(k,313) = -rxt(k,427)*y(k,253)
         mat(k,242) = -rxt(k,428)*y(k,253)
         mat(k,220) = -rxt(k,440)*y(k,253)
         mat(k,368) = -rxt(k,443)*y(k,253)
         mat(k,727) = -rxt(k,447)*y(k,253)
         mat(k,237) = -rxt(k,448)*y(k,253)
         mat(k,259) = -rxt(k,450)*y(k,253)
         mat(k,747) = -rxt(k,453)*y(k,253)
         mat(k,264) = -rxt(k,456)*y(k,253)
         mat(k,482) = -rxt(k,457)*y(k,253)
         mat(k,1028) = -rxt(k,460)*y(k,253)
         mat(k,951) = -rxt(k,463)*y(k,253)
         mat(k,453) = -rxt(k,468)*y(k,253)
         mat(k,714) = -rxt(k,472)*y(k,253)
         mat(k,666) = -rxt(k,473)*y(k,253)
         mat(k,520) = -rxt(k,477)*y(k,253)
         mat(k,1040) = -rxt(k,479)*y(k,253)
         mat(k,1122) = -rxt(k,480)*y(k,253)
         mat(k,352) = -(rxt(k,482) + rxt(k,496)) * y(k,253)
         mat(k,417) = -rxt(k,484)*y(k,253)
         mat(k,905) = -rxt(k,486)*y(k,253)
         mat(k,560) = -rxt(k,490)*y(k,253)
         mat(k,1281) = -rxt(k,495)*y(k,253)
         mat(k,139) = -rxt(k,498)*y(k,253)
         mat(k,1028) = mat(k,1028) + .630_r8*rxt(k,459)*y(k,140)
         mat(k,325) = mat(k,325) + .650_r8*rxt(k,296)*y(k,253)
         mat(k,606) = mat(k,606) + .130_r8*rxt(k,298)*y(k,140)
         mat(k,379) = mat(k,379) + .500_r8*rxt(k,304)*y(k,253)
         mat(k,1058) = mat(k,1058) + .360_r8*rxt(k,327)*y(k,140)
         mat(k,1536) = mat(k,1536) + rxt(k,276)*y(k,139)
         mat(k,457) = mat(k,457) + .300_r8*rxt(k,283)*y(k,253)
         mat(k,1483) = mat(k,1483) + rxt(k,290)*y(k,252)
         mat(k,2123) = rxt(k,199)*y(k,242)
         mat(k,914) = rxt(k,253)*y(k,263)
         mat(k,1515) = rxt(k,160)*y(k,140) + 2.000_r8*rxt(k,155)*y(k,242)
         mat(k,1454) = mat(k,1454) + rxt(k,152)*y(k,139) + rxt(k,144)*y(k,252)
         mat(k,648) = mat(k,648) + rxt(k,153)*y(k,139)
         mat(k,851) = mat(k,851) + rxt(k,243)*y(k,139) + rxt(k,249)*y(k,252)
         mat(k,2025) = mat(k,2025) + rxt(k,214)*y(k,139) + rxt(k,226)*y(k,252)
         mat(k,229) = mat(k,229) + rxt(k,293)*y(k,252)
         mat(k,826) = rxt(k,245)*y(k,139)
         mat(k,857) = mat(k,857) + rxt(k,217)*y(k,139)
         mat(k,890) = mat(k,890) + .320_r8*rxt(k,404)*y(k,140)
         mat(k,785) = mat(k,785) + .600_r8*rxt(k,406)*y(k,253)
         mat(k,1260) = mat(k,1260) + .240_r8*rxt(k,357)*y(k,140)
         mat(k,387) = mat(k,387) + .100_r8*rxt(k,359)*y(k,253)
         mat(k,951) = mat(k,951) + .630_r8*rxt(k,462)*y(k,140)
         mat(k,1387) = mat(k,1387) + .360_r8*rxt(k,371)*y(k,140)
         mat(k,1950) = rxt(k,183)*y(k,242)
         mat(k,2084) = mat(k,2084) + rxt(k,178)*y(k,242)
         mat(k,1610) = mat(k,1610) + rxt(k,276)*y(k,44) + rxt(k,152)*y(k,79) &
                      + rxt(k,153)*y(k,81) + rxt(k,243)*y(k,83) + rxt(k,214)*y(k,87) &
                      + rxt(k,245)*y(k,93) + rxt(k,217)*y(k,94) + rxt(k,158)*y(k,242)
         mat(k,2185) = mat(k,2185) + .630_r8*rxt(k,459)*y(k,6) + .130_r8*rxt(k,298) &
                      *y(k,27) + .360_r8*rxt(k,327)*y(k,31) + rxt(k,160)*y(k,78) &
                      + .320_r8*rxt(k,404)*y(k,100) + .240_r8*rxt(k,357)*y(k,111) &
                      + .630_r8*rxt(k,462)*y(k,116) + .360_r8*rxt(k,371)*y(k,118) &
                      + rxt(k,159)*y(k,242)
         mat(k,590) = mat(k,590) + .500_r8*rxt(k,339)*y(k,253)
         mat(k,250) = mat(k,250) + .500_r8*rxt(k,414)*y(k,253)
         mat(k,566) = .400_r8*rxt(k,415)*y(k,242)
         mat(k,1439) = .450_r8*rxt(k,312)*y(k,242)
         mat(k,808) = .400_r8*rxt(k,429)*y(k,242)
         mat(k,2303) = mat(k,2303) + rxt(k,199)*y(k,58) + 2.000_r8*rxt(k,155)*y(k,78) &
                      + rxt(k,183)*y(k,129) + rxt(k,178)*y(k,131) + rxt(k,158) &
                      *y(k,139) + rxt(k,159)*y(k,140) + .400_r8*rxt(k,415)*y(k,229) &
                      + .450_r8*rxt(k,312)*y(k,236) + .400_r8*rxt(k,429)*y(k,238) &
                      + .450_r8*rxt(k,362)*y(k,248) + .400_r8*rxt(k,435)*y(k,249) &
                      + .200_r8*rxt(k,366)*y(k,250) + .150_r8*rxt(k,341)*y(k,256)
         mat(k,1407) = .450_r8*rxt(k,362)*y(k,242)
         mat(k,923) = .400_r8*rxt(k,435)*y(k,242)
         mat(k,703) = .200_r8*rxt(k,366)*y(k,242)
         mat(k,1651) = rxt(k,290)*y(k,56) + rxt(k,144)*y(k,79) + rxt(k,249)*y(k,83) &
                      + rxt(k,226)*y(k,87) + rxt(k,293)*y(k,88) + 2.000_r8*rxt(k,145) &
                      *y(k,263)
         mat(k,1823) = mat(k,1823) + .650_r8*rxt(k,296)*y(k,26) + .500_r8*rxt(k,304) &
                      *y(k,29) + .300_r8*rxt(k,283)*y(k,55) + .600_r8*rxt(k,406) &
                      *y(k,106) + .100_r8*rxt(k,359)*y(k,112) + .500_r8*rxt(k,339) &
                      *y(k,154) + .500_r8*rxt(k,414)*y(k,221)
         mat(k,1205) = .150_r8*rxt(k,341)*y(k,242)
         mat(k,2355) = rxt(k,253)*y(k,75) + 2.000_r8*rxt(k,145)*y(k,252)
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
         mat(k,499) = -(rxt(k,438)*y(k,242) + rxt(k,439)*y(k,129))
         mat(k,2238) = -rxt(k,438)*y(k,254)
         mat(k,1889) = -rxt(k,439)*y(k,254)
         mat(k,240) = .200_r8*rxt(k,428)*y(k,253)
         mat(k,218) = .140_r8*rxt(k,440)*y(k,253)
         mat(k,366) = rxt(k,443)*y(k,253)
         mat(k,1734) = .200_r8*rxt(k,428)*y(k,68) + .140_r8*rxt(k,440)*y(k,148) &
                      + rxt(k,443)*y(k,149)
         mat(k,813) = -(rxt(k,337)*y(k,242) + rxt(k,338)*y(k,129))
         mat(k,2261) = -rxt(k,337)*y(k,255)
         mat(k,1910) = -rxt(k,338)*y(k,255)
         mat(k,1044) = rxt(k,344)*y(k,253)
         mat(k,586) = .500_r8*rxt(k,339)*y(k,253)
         mat(k,1771) = rxt(k,344)*y(k,31) + .500_r8*rxt(k,339)*y(k,154)
         mat(k,1200) = -(rxt(k,340)*y(k,237) + rxt(k,341)*y(k,242) + rxt(k,342) &
                      *y(k,129))
         mat(k,1984) = -rxt(k,340)*y(k,256)
         mat(k,2283) = -rxt(k,341)*y(k,256)
         mat(k,1931) = -rxt(k,342)*y(k,256)
         mat(k,1022) = .060_r8*rxt(k,459)*y(k,140)
         mat(k,983) = rxt(k,345)*y(k,253)
         mat(k,946) = .060_r8*rxt(k,462)*y(k,140)
         mat(k,2166) = .060_r8*rxt(k,459)*y(k,6) + .060_r8*rxt(k,462)*y(k,116)
         mat(k,419) = rxt(k,343)*y(k,253)
         mat(k,1118) = .150_r8*rxt(k,480)*y(k,253)
         mat(k,1802) = rxt(k,345)*y(k,50) + rxt(k,343)*y(k,155) + .150_r8*rxt(k,480) &
                      *y(k,217)
         mat(k,1163) = -(rxt(k,469)*y(k,237) + rxt(k,470)*y(k,242) + rxt(k,471) &
                      *y(k,129))
         mat(k,1982) = -rxt(k,469)*y(k,257)
         mat(k,2281) = -rxt(k,470)*y(k,257)
         mat(k,1929) = -rxt(k,471)*y(k,257)
         mat(k,2062) = .500_r8*rxt(k,478)*y(k,216)
         mat(k,712) = rxt(k,472)*y(k,253)
         mat(k,1039) = .500_r8*rxt(k,478)*y(k,131) + rxt(k,479)*y(k,253)
         mat(k,1800) = rxt(k,472)*y(k,213) + rxt(k,479)*y(k,216)
         mat(k,1184) = -(rxt(k,474)*y(k,237) + rxt(k,475)*y(k,242) + rxt(k,476) &
                      *y(k,129))
         mat(k,1983) = -rxt(k,474)*y(k,258)
         mat(k,2282) = -rxt(k,475)*y(k,258)
         mat(k,1930) = -rxt(k,476)*y(k,258)
         mat(k,1021) = rxt(k,460)*y(k,253)
         mat(k,945) = rxt(k,463)*y(k,253)
         mat(k,518) = rxt(k,477)*y(k,253)
         mat(k,1801) = rxt(k,460)*y(k,6) + rxt(k,463)*y(k,116) + rxt(k,477)*y(k,215)
         mat(k,769) = -(rxt(k,445)*y(k,242) + rxt(k,446)*y(k,129))
         mat(k,2257) = -rxt(k,445)*y(k,259)
         mat(k,1907) = -rxt(k,446)*y(k,259)
         mat(k,721) = rxt(k,447)*y(k,253)
         mat(k,236) = .650_r8*rxt(k,448)*y(k,253)
         mat(k,1767) = rxt(k,447)*y(k,218) + .650_r8*rxt(k,448)*y(k,219)
         mat(k,1217) = -(rxt(k,409)*y(k,236) + rxt(k,410)*y(k,237) + rxt(k,411) &
                      *y(k,242) + rxt(k,412)*y(k,129) + rxt(k,413)*y(k,131))
         mat(k,1424) = -rxt(k,409)*y(k,260)
         mat(k,1985) = -rxt(k,410)*y(k,260)
         mat(k,2284) = -rxt(k,411)*y(k,260)
         mat(k,1932) = -rxt(k,412)*y(k,260)
         mat(k,2065) = -rxt(k,413)*y(k,260)
         mat(k,280) = rxt(k,381)*y(k,253)
         mat(k,382) = rxt(k,382)*y(k,253)
         mat(k,175) = rxt(k,383)*y(k,253)
         mat(k,781) = .400_r8*rxt(k,406)*y(k,253)
         mat(k,249) = .500_r8*rxt(k,414)*y(k,253)
         mat(k,1803) = rxt(k,381)*y(k,96) + rxt(k,382)*y(k,98) + rxt(k,383)*y(k,99) &
                      + .400_r8*rxt(k,406)*y(k,106) + .500_r8*rxt(k,414)*y(k,221)
         mat(k,793) = -(rxt(k,451)*y(k,242) + rxt(k,452)*y(k,129))
         mat(k,2259) = -rxt(k,451)*y(k,261)
         mat(k,1908) = -rxt(k,452)*y(k,261)
         mat(k,256) = .560_r8*rxt(k,450)*y(k,253)
         mat(k,740) = rxt(k,453)*y(k,253)
         mat(k,1769) = .560_r8*rxt(k,450)*y(k,222) + rxt(k,453)*y(k,223)
         mat(k,544) = -(rxt(k,454)*y(k,242) + rxt(k,455)*y(k,129))
         mat(k,2243) = -rxt(k,454)*y(k,262)
         mat(k,1894) = -rxt(k,455)*y(k,262)
         mat(k,263) = .300_r8*rxt(k,456)*y(k,253)
         mat(k,479) = rxt(k,457)*y(k,253)
         mat(k,1741) = .300_r8*rxt(k,456)*y(k,225) + rxt(k,457)*y(k,226)
         mat(k,2365) = -(rxt(k,145)*y(k,252) + rxt(k,253)*y(k,75) + rxt(k,497) &
                      *y(k,160))
         mat(k,1661) = -rxt(k,145)*y(k,263)
         mat(k,917) = -rxt(k,253)*y(k,263)
         mat(k,305) = -rxt(k,497)*y(k,263)
         mat(k,332) = rxt(k,306)*y(k,253)
         mat(k,441) = rxt(k,331)*y(k,253)
         mat(k,157) = rxt(k,332)*y(k,253)
         mat(k,514) = rxt(k,258)*y(k,253)
         mat(k,1545) = rxt(k,277)*y(k,253)
         mat(k,634) = rxt(k,260)*y(k,253)
         mat(k,173) = rxt(k,261)*y(k,253)
         mat(k,1089) = rxt(k,308)*y(k,253)
         mat(k,409) = rxt(k,263)*y(k,253)
         mat(k,987) = rxt(k,345)*y(k,253)
         mat(k,1272) = rxt(k,334)*y(k,253)
         mat(k,733) = rxt(k,314)*y(k,253)
         mat(k,641) = rxt(k,315)*y(k,253)
         mat(k,459) = rxt(k,283)*y(k,253)
         mat(k,1489) = rxt(k,284)*y(k,253)
         mat(k,1522) = rxt(k,156)*y(k,242)
         mat(k,1459) = rxt(k,161)*y(k,253)
         mat(k,652) = rxt(k,162)*y(k,253)
         mat(k,853) = rxt(k,244)*y(k,253)
         mat(k,339) = rxt(k,268)*y(k,253)
         mat(k,2035) = (rxt(k,553)+rxt(k,558))*y(k,93) + (rxt(k,546)+rxt(k,552) &
                       +rxt(k,557))*y(k,94) + rxt(k,215)*y(k,253)
         mat(k,994) = rxt(k,286)*y(k,253)
         mat(k,1506) = rxt(k,191)*y(k,253)
         mat(k,447) = rxt(k,169)*y(k,253)
         mat(k,830) = (rxt(k,553)+rxt(k,558))*y(k,87)
         mat(k,861) = (rxt(k,546)+rxt(k,552)+rxt(k,557))*y(k,87) + rxt(k,218)*y(k,253)
         mat(k,1263) = .500_r8*rxt(k,358)*y(k,253)
         mat(k,140) = rxt(k,498)*y(k,253)
         mat(k,592) = rxt(k,339)*y(k,253)
         mat(k,423) = rxt(k,343)*y(k,253)
         mat(k,2313) = rxt(k,156)*y(k,78) + rxt(k,163)*y(k,253)
         mat(k,1833) = rxt(k,306)*y(k,30) + rxt(k,331)*y(k,32) + rxt(k,332)*y(k,33) &
                      + rxt(k,258)*y(k,43) + rxt(k,277)*y(k,44) + rxt(k,260)*y(k,45) &
                      + rxt(k,261)*y(k,46) + rxt(k,308)*y(k,47) + rxt(k,263)*y(k,48) &
                      + rxt(k,345)*y(k,50) + rxt(k,334)*y(k,51) + rxt(k,314)*y(k,52) &
                      + rxt(k,315)*y(k,53) + rxt(k,283)*y(k,55) + rxt(k,284)*y(k,56) &
                      + rxt(k,161)*y(k,79) + rxt(k,162)*y(k,81) + rxt(k,244)*y(k,83) &
                      + rxt(k,268)*y(k,86) + rxt(k,215)*y(k,87) + rxt(k,286)*y(k,89) &
                      + rxt(k,191)*y(k,91) + rxt(k,169)*y(k,92) + rxt(k,218)*y(k,94) &
                      + .500_r8*rxt(k,358)*y(k,111) + rxt(k,498)*y(k,127) + rxt(k,339) &
                      *y(k,154) + rxt(k,343)*y(k,155) + rxt(k,163)*y(k,242) &
                      + 2.000_r8*rxt(k,166)*y(k,253)
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
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 147) = lmat(k, 147)
         mat(k, 148) = lmat(k, 148)
         mat(k, 149) = lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 165) = mat(k, 165) + lmat(k, 165)
         mat(k, 166) = mat(k, 166) + lmat(k, 166)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = lmat(k, 178)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 213) = lmat(k, 213)
         mat(k, 214) = lmat(k, 214)
         mat(k, 215) = lmat(k, 215)
         mat(k, 216) = lmat(k, 216)
         mat(k, 217) = mat(k, 217) + lmat(k, 217)
         mat(k, 222) = lmat(k, 222)
         mat(k, 223) = lmat(k, 223)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 269) = lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 278) = mat(k, 278) + lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 290) = mat(k, 290) + lmat(k, 290)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 295) = mat(k, 295) + lmat(k, 295)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 303) = lmat(k, 303)
         mat(k, 304) = lmat(k, 304)
         mat(k, 306) = mat(k, 306) + lmat(k, 306)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 314) = lmat(k, 314)
         mat(k, 315) = lmat(k, 315)
         mat(k, 316) = lmat(k, 316)
         mat(k, 317) = lmat(k, 317)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 338) = mat(k, 338) + lmat(k, 338)
         mat(k, 340) = lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = lmat(k, 345)
         mat(k, 346) = lmat(k, 346)
         mat(k, 347) = lmat(k, 347)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 355) = lmat(k, 355)
         mat(k, 356) = lmat(k, 356)
         mat(k, 357) = lmat(k, 357)
         mat(k, 358) = lmat(k, 358)
         mat(k, 359) = lmat(k, 359)
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
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 391) = lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 402) = mat(k, 402) + lmat(k, 402)
         mat(k, 405) = lmat(k, 405)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 411) = lmat(k, 411)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 420) = lmat(k, 420)
         mat(k, 421) = lmat(k, 421)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 425) = lmat(k, 425)
         mat(k, 427) = lmat(k, 427)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 430) = lmat(k, 430)
         mat(k, 431) = lmat(k, 431)
         mat(k, 432) = lmat(k, 432)
         mat(k, 433) = lmat(k, 433)
         mat(k, 434) = lmat(k, 434)
         mat(k, 435) = lmat(k, 435)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 438) = lmat(k, 438)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 440) = lmat(k, 440)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 444) = mat(k, 444) + lmat(k, 444)
         mat(k, 445) = lmat(k, 445)
         mat(k, 446) = lmat(k, 446)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 449) = lmat(k, 449)
         mat(k, 451) = lmat(k, 451)
         mat(k, 452) = lmat(k, 452)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 455) = lmat(k, 455)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 461) = mat(k, 461) + lmat(k, 461)
         mat(k, 466) = mat(k, 466) + lmat(k, 466)
         mat(k, 469) = lmat(k, 469)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 473) = lmat(k, 473)
         mat(k, 475) = lmat(k, 475)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 477) = lmat(k, 477)
         mat(k, 478) = mat(k, 478) + lmat(k, 478)
         mat(k, 480) = lmat(k, 480)
         mat(k, 481) = lmat(k, 481)
         mat(k, 482) = mat(k, 482) + lmat(k, 482)
         mat(k, 483) = lmat(k, 483)
         mat(k, 486) = mat(k, 486) + lmat(k, 486)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 494) = lmat(k, 494)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 499) = mat(k, 499) + lmat(k, 499)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 516) = lmat(k, 516)
         mat(k, 517) = lmat(k, 517)
         mat(k, 519) = lmat(k, 519)
         mat(k, 520) = mat(k, 520) + lmat(k, 520)
         mat(k, 521) = lmat(k, 521)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 551) = mat(k, 551) + lmat(k, 551)
         mat(k, 553) = lmat(k, 553)
         mat(k, 554) = lmat(k, 554)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 557) = mat(k, 557) + lmat(k, 557)
         mat(k, 559) = lmat(k, 559)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 570) = lmat(k, 570)
         mat(k, 571) = lmat(k, 571)
         mat(k, 573) = lmat(k, 573)
         mat(k, 575) = lmat(k, 575)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 577) = mat(k, 577) + lmat(k, 577)
         mat(k, 578) = lmat(k, 578)
         mat(k, 579) = lmat(k, 579)
         mat(k, 580) = lmat(k, 580)
         mat(k, 581) = lmat(k, 581)
         mat(k, 583) = mat(k, 583) + lmat(k, 583)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 585) = mat(k, 585) + lmat(k, 585)
         mat(k, 587) = lmat(k, 587)
         mat(k, 589) = lmat(k, 589)
         mat(k, 590) = mat(k, 590) + lmat(k, 590)
         mat(k, 591) = lmat(k, 591)
         mat(k, 593) = mat(k, 593) + lmat(k, 593)
         mat(k, 598) = lmat(k, 598)
         mat(k, 601) = mat(k, 601) + lmat(k, 601)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 621) = lmat(k, 621)
         mat(k, 626) = mat(k, 626) + lmat(k, 626)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 630) = lmat(k, 630)
         mat(k, 635) = mat(k, 635) + lmat(k, 635)
         mat(k, 636) = mat(k, 636) + lmat(k, 636)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 640) = lmat(k, 640)
         mat(k, 642) = lmat(k, 642)
         mat(k, 643) = lmat(k, 643)
         mat(k, 644) = lmat(k, 644)
         mat(k, 645) = lmat(k, 645)
         mat(k, 646) = mat(k, 646) + lmat(k, 646)
         mat(k, 648) = mat(k, 648) + lmat(k, 648)
         mat(k, 653) = mat(k, 653) + lmat(k, 653)
         mat(k, 656) = lmat(k, 656)
         mat(k, 658) = lmat(k, 658)
         mat(k, 662) = mat(k, 662) + lmat(k, 662)
         mat(k, 663) = mat(k, 663) + lmat(k, 663)
         mat(k, 664) = lmat(k, 664)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 667) = lmat(k, 667)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 676) = lmat(k, 676)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 681) = lmat(k, 681)
         mat(k, 682) = lmat(k, 682)
         mat(k, 684) = lmat(k, 684)
         mat(k, 685) = mat(k, 685) + lmat(k, 685)
         mat(k, 686) = lmat(k, 686)
         mat(k, 687) = lmat(k, 687)
         mat(k, 688) = mat(k, 688) + lmat(k, 688)
         mat(k, 691) = mat(k, 691) + lmat(k, 691)
         mat(k, 692) = mat(k, 692) + lmat(k, 692)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 697) = lmat(k, 697)
         mat(k, 699) = mat(k, 699) + lmat(k, 699)
         mat(k, 706) = mat(k, 706) + lmat(k, 706)
         mat(k, 707) = lmat(k, 707)
         mat(k, 708) = lmat(k, 708)
         mat(k, 709) = lmat(k, 709)
         mat(k, 710) = lmat(k, 710)
         mat(k, 711) = lmat(k, 711)
         mat(k, 713) = lmat(k, 713)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 715) = lmat(k, 715)
         mat(k, 716) = lmat(k, 716)
         mat(k, 717) = lmat(k, 717)
         mat(k, 718) = lmat(k, 718)
         mat(k, 719) = mat(k, 719) + lmat(k, 719)
         mat(k, 724) = lmat(k, 724)
         mat(k, 726) = lmat(k, 726)
         mat(k, 727) = mat(k, 727) + lmat(k, 727)
         mat(k, 728) = lmat(k, 728)
         mat(k, 729) = mat(k, 729) + lmat(k, 729)
         mat(k, 734) = lmat(k, 734)
         mat(k, 735) = lmat(k, 735)
         mat(k, 736) = lmat(k, 736)
         mat(k, 737) = lmat(k, 737)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 743) = lmat(k, 743)
         mat(k, 745) = lmat(k, 745)
         mat(k, 747) = mat(k, 747) + lmat(k, 747)
         mat(k, 748) = lmat(k, 748)
         mat(k, 751) = mat(k, 751) + lmat(k, 751)
         mat(k, 758) = mat(k, 758) + lmat(k, 758)
         mat(k, 769) = mat(k, 769) + lmat(k, 769)
         mat(k, 780) = mat(k, 780) + lmat(k, 780)
         mat(k, 782) = lmat(k, 782)
         mat(k, 783) = lmat(k, 783)
         mat(k, 784) = lmat(k, 784)
         mat(k, 785) = mat(k, 785) + lmat(k, 785)
         mat(k, 786) = lmat(k, 786)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 804) = mat(k, 804) + lmat(k, 804)
         mat(k, 813) = mat(k, 813) + lmat(k, 813)
         mat(k, 823) = mat(k, 823) + lmat(k, 823)
         mat(k, 824) = lmat(k, 824)
         mat(k, 826) = mat(k, 826) + lmat(k, 826)
         mat(k, 831) = mat(k, 831) + lmat(k, 831)
         mat(k, 836) = mat(k, 836) + lmat(k, 836)
         mat(k, 846) = mat(k, 846) + lmat(k, 846)
         mat(k, 847) = mat(k, 847) + lmat(k, 847)
         mat(k, 848) = mat(k, 848) + lmat(k, 848)
         mat(k, 855) = mat(k, 855) + lmat(k, 855)
         mat(k, 857) = mat(k, 857) + lmat(k, 857)
         mat(k, 859) = mat(k, 859) + lmat(k, 859)
         mat(k, 865) = mat(k, 865) + lmat(k, 865)
         mat(k, 873) = lmat(k, 873)
         mat(k, 874) = lmat(k, 874)
         mat(k, 875) = lmat(k, 875)
         mat(k, 879) = mat(k, 879) + lmat(k, 879)
         mat(k, 895) = mat(k, 895) + lmat(k, 895)
         mat(k, 897) = lmat(k, 897)
         mat(k, 898) = lmat(k, 898)
         mat(k, 899) = mat(k, 899) + lmat(k, 899)
         mat(k, 902) = mat(k, 902) + lmat(k, 902)
         mat(k, 903) = lmat(k, 903)
         mat(k, 904) = lmat(k, 904)
         mat(k, 909) = mat(k, 909) + lmat(k, 909)
         mat(k, 918) = mat(k, 918) + lmat(k, 918)
         mat(k, 936) = mat(k, 936) + lmat(k, 936)
         mat(k, 960) = mat(k, 960) + lmat(k, 960)
         mat(k, 971) = mat(k, 971) + lmat(k, 971)
         mat(k, 972) = mat(k, 972) + lmat(k, 972)
         mat(k, 973) = mat(k, 973) + lmat(k, 973)
         mat(k, 974) = lmat(k, 974)
         mat(k, 978) = mat(k, 978) + lmat(k, 978)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 980) = mat(k, 980) + lmat(k, 980)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 984) = lmat(k, 984)
         mat(k, 986) = lmat(k, 986)
         mat(k, 989) = mat(k, 989) + lmat(k, 989)
         mat(k, 995) = lmat(k, 995)
         mat(k, 997) = mat(k, 997) + lmat(k, 997)
         mat(k, 998) = mat(k, 998) + lmat(k, 998)
         mat(k,1001) = mat(k,1001) + lmat(k,1001)
         mat(k,1015) = mat(k,1015) + lmat(k,1015)
         mat(k,1035) = mat(k,1035) + lmat(k,1035)
         mat(k,1037) = lmat(k,1037)
         mat(k,1038) = lmat(k,1038)
         mat(k,1042) = lmat(k,1042)
         mat(k,1048) = mat(k,1048) + lmat(k,1048)
         mat(k,1065) = lmat(k,1065)
         mat(k,1069) = mat(k,1069) + lmat(k,1069)
         mat(k,1073) = lmat(k,1073)
         mat(k,1075) = lmat(k,1075)
         mat(k,1080) = mat(k,1080) + lmat(k,1080)
         mat(k,1081) = mat(k,1081) + lmat(k,1081)
         mat(k,1082) = lmat(k,1082)
         mat(k,1086) = lmat(k,1086)
         mat(k,1088) = lmat(k,1088)
         mat(k,1092) = mat(k,1092) + lmat(k,1092)
         mat(k,1093) = lmat(k,1093)
         mat(k,1094) = mat(k,1094) + lmat(k,1094)
         mat(k,1097) = mat(k,1097) + lmat(k,1097)
         mat(k,1102) = mat(k,1102) + lmat(k,1102)
         mat(k,1114) = mat(k,1114) + lmat(k,1114)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1116) = mat(k,1116) + lmat(k,1116)
         mat(k,1117) = mat(k,1117) + lmat(k,1117)
         mat(k,1118) = mat(k,1118) + lmat(k,1118)
         mat(k,1119) = mat(k,1119) + lmat(k,1119)
         mat(k,1120) = mat(k,1120) + lmat(k,1120)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1127) = mat(k,1127) + lmat(k,1127)
         mat(k,1132) = lmat(k,1132)
         mat(k,1133) = lmat(k,1133)
         mat(k,1134) = lmat(k,1134)
         mat(k,1135) = lmat(k,1135)
         mat(k,1136) = mat(k,1136) + lmat(k,1136)
         mat(k,1137) = lmat(k,1137)
         mat(k,1139) = lmat(k,1139)
         mat(k,1140) = lmat(k,1140)
         mat(k,1141) = lmat(k,1141)
         mat(k,1142) = lmat(k,1142)
         mat(k,1147) = mat(k,1147) + lmat(k,1147)
         mat(k,1149) = mat(k,1149) + lmat(k,1149)
         mat(k,1151) = lmat(k,1151)
         mat(k,1152) = lmat(k,1152)
         mat(k,1154) = mat(k,1154) + lmat(k,1154)
         mat(k,1163) = mat(k,1163) + lmat(k,1163)
         mat(k,1184) = mat(k,1184) + lmat(k,1184)
         mat(k,1200) = mat(k,1200) + lmat(k,1200)
         mat(k,1217) = mat(k,1217) + lmat(k,1217)
         mat(k,1237) = mat(k,1237) + lmat(k,1237)
         mat(k,1252) = mat(k,1252) + lmat(k,1252)
         mat(k,1253) = mat(k,1253) + lmat(k,1253)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1257) = mat(k,1257) + lmat(k,1257)
         mat(k,1258) = mat(k,1258) + lmat(k,1258)
         mat(k,1262) = mat(k,1262) + lmat(k,1262)
         mat(k,1264) = mat(k,1264) + lmat(k,1264)
         mat(k,1265) = mat(k,1265) + lmat(k,1265)
         mat(k,1266) = mat(k,1266) + lmat(k,1266)
         mat(k,1271) = lmat(k,1271)
         mat(k,1274) = lmat(k,1274)
         mat(k,1275) = mat(k,1275) + lmat(k,1275)
         mat(k,1276) = mat(k,1276) + lmat(k,1276)
         mat(k,1280) = lmat(k,1280)
         mat(k,1301) = mat(k,1301) + lmat(k,1301)
         mat(k,1317) = lmat(k,1317)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1348) = mat(k,1348) + lmat(k,1348)
         mat(k,1358) = mat(k,1358) + lmat(k,1358)
         mat(k,1373) = lmat(k,1373)
         mat(k,1375) = mat(k,1375) + lmat(k,1375)
         mat(k,1379) = mat(k,1379) + lmat(k,1379)
         mat(k,1381) = mat(k,1381) + lmat(k,1381)
         mat(k,1389) = lmat(k,1389)
         mat(k,1401) = mat(k,1401) + lmat(k,1401)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1448) = mat(k,1448) + lmat(k,1448)
         mat(k,1462) = mat(k,1462) + lmat(k,1462)
         mat(k,1473) = lmat(k,1473)
         mat(k,1475) = lmat(k,1475)
         mat(k,1476) = mat(k,1476) + lmat(k,1476)
         mat(k,1477) = mat(k,1477) + lmat(k,1477)
         mat(k,1479) = mat(k,1479) + lmat(k,1479)
         mat(k,1480) = mat(k,1480) + lmat(k,1480)
         mat(k,1481) = lmat(k,1481)
         mat(k,1483) = mat(k,1483) + lmat(k,1483)
         mat(k,1484) = mat(k,1484) + lmat(k,1484)
         mat(k,1489) = mat(k,1489) + lmat(k,1489)
         mat(k,1494) = mat(k,1494) + lmat(k,1494)
         mat(k,1497) = lmat(k,1497)
         mat(k,1500) = mat(k,1500) + lmat(k,1500)
         mat(k,1510) = mat(k,1510) + lmat(k,1510)
         mat(k,1521) = mat(k,1521) + lmat(k,1521)
         mat(k,1526) = mat(k,1526) + lmat(k,1526)
         mat(k,1527) = lmat(k,1527)
         mat(k,1531) = mat(k,1531) + lmat(k,1531)
         mat(k,1532) = mat(k,1532) + lmat(k,1532)
         mat(k,1573) = mat(k,1573) + lmat(k,1573)
         mat(k,1576) = mat(k,1576) + lmat(k,1576)
         mat(k,1577) = mat(k,1577) + lmat(k,1577)
         mat(k,1579) = mat(k,1579) + lmat(k,1579)
         mat(k,1581) = mat(k,1581) + lmat(k,1581)
         mat(k,1608) = mat(k,1608) + lmat(k,1608)
         mat(k,1617) = mat(k,1617) + lmat(k,1617)
         mat(k,1649) = lmat(k,1649)
         mat(k,1650) = mat(k,1650) + lmat(k,1650)
         mat(k,1823) = mat(k,1823) + lmat(k,1823)
         mat(k,1840) = mat(k,1840) + lmat(k,1840)
         mat(k,1845) = mat(k,1845) + lmat(k,1845)
         mat(k,1848) = mat(k,1848) + lmat(k,1848)
         mat(k,1892) = mat(k,1892) + lmat(k,1892)
         mat(k,1948) = mat(k,1948) + lmat(k,1948)
         mat(k,1952) = mat(k,1952) + lmat(k,1952)
         mat(k,2005) = mat(k,2005) + lmat(k,2005)
         mat(k,2020) = mat(k,2020) + lmat(k,2020)
         mat(k,2029) = mat(k,2029) + lmat(k,2029)
         mat(k,2031) = mat(k,2031) + lmat(k,2031)
         mat(k,2078) = mat(k,2078) + lmat(k,2078)
         mat(k,2081) = mat(k,2081) + lmat(k,2081)
         mat(k,2082) = mat(k,2082) + lmat(k,2082)
         mat(k,2086) = mat(k,2086) + lmat(k,2086)
         mat(k,2089) = mat(k,2089) + lmat(k,2089)
         mat(k,2129) = mat(k,2129) + lmat(k,2129)
         mat(k,2183) = mat(k,2183) + lmat(k,2183)
         mat(k,2184) = mat(k,2184) + lmat(k,2184)
         mat(k,2192) = mat(k,2192) + lmat(k,2192)
         mat(k,2247) = mat(k,2247) + lmat(k,2247)
         mat(k,2311) = mat(k,2311) + lmat(k,2311)
         mat(k,2327) = mat(k,2327) + lmat(k,2327)
         mat(k,2335) = mat(k,2335) + lmat(k,2335)
         mat(k,2338) = mat(k,2338) + lmat(k,2338)
         mat(k,2346) = lmat(k,2346)
         mat(k,2350) = lmat(k,2350)
         mat(k,2353) = lmat(k,2353)
         mat(k,2354) = mat(k,2354) + lmat(k,2354)
         mat(k,2355) = mat(k,2355) + lmat(k,2355)
         mat(k,2365) = mat(k,2365) + lmat(k,2365)
         mat(k, 257) = 0._r8
         mat(k, 258) = 0._r8
         mat(k, 297) = 0._r8
         mat(k, 335) = 0._r8
         mat(k, 373) = 0._r8
         mat(k, 487) = 0._r8
         mat(k, 489) = 0._r8
         mat(k, 502) = 0._r8
         mat(k, 531) = 0._r8
         mat(k, 537) = 0._r8
         mat(k, 540) = 0._r8
         mat(k, 548) = 0._r8
         mat(k, 672) = 0._r8
         mat(k, 674) = 0._r8
         mat(k, 678) = 0._r8
         mat(k, 679) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 693) = 0._r8
         mat(k, 720) = 0._r8
         mat(k, 722) = 0._r8
         mat(k, 723) = 0._r8
         mat(k, 725) = 0._r8
         mat(k, 739) = 0._r8
         mat(k, 741) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 744) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 768) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 773) = 0._r8
         mat(k, 776) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 799) = 0._r8
         mat(k, 801) = 0._r8
         mat(k, 815) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 821) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 841) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 940) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 990) = 0._r8
         mat(k, 991) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1011) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1023) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1029) = 0._r8
         mat(k,1030) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1052) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1056) = 0._r8
         mat(k,1057) = 0._r8
         mat(k,1059) = 0._r8
         mat(k,1064) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1070) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1072) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1079) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1222) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1230) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1240) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1288) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1320) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1327) = 0._r8
         mat(k,1328) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1331) = 0._r8
         mat(k,1332) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1349) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1364) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1370) = 0._r8
         mat(k,1376) = 0._r8
         mat(k,1380) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1390) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1400) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1413) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1442) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1466) = 0._r8
         mat(k,1467) = 0._r8
         mat(k,1471) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1486) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1493) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1502) = 0._r8
         mat(k,1504) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1512) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1516) = 0._r8
         mat(k,1517) = 0._r8
         mat(k,1518) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1556) = 0._r8
         mat(k,1559) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1564) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1566) = 0._r8
         mat(k,1567) = 0._r8
         mat(k,1571) = 0._r8
         mat(k,1572) = 0._r8
         mat(k,1574) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1578) = 0._r8
         mat(k,1582) = 0._r8
         mat(k,1583) = 0._r8
         mat(k,1585) = 0._r8
         mat(k,1589) = 0._r8
         mat(k,1591) = 0._r8
         mat(k,1597) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1609) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1648) = 0._r8
         mat(k,1656) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1770) = 0._r8
         mat(k,1781) = 0._r8
         mat(k,1782) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1841) = 0._r8
         mat(k,1842) = 0._r8
         mat(k,1843) = 0._r8
         mat(k,1846) = 0._r8
         mat(k,1850) = 0._r8
         mat(k,1851) = 0._r8
         mat(k,1852) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1943) = 0._r8
         mat(k,1944) = 0._r8
         mat(k,1945) = 0._r8
         mat(k,1949) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1960) = 0._r8
         mat(k,1971) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1996) = 0._r8
         mat(k,1997) = 0._r8
         mat(k,2000) = 0._r8
         mat(k,2001) = 0._r8
         mat(k,2002) = 0._r8
         mat(k,2003) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2007) = 0._r8
         mat(k,2009) = 0._r8
         mat(k,2012) = 0._r8
         mat(k,2018) = 0._r8
         mat(k,2021) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2026) = 0._r8
         mat(k,2027) = 0._r8
         mat(k,2028) = 0._r8
         mat(k,2030) = 0._r8
         mat(k,2032) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2042) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2049) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2060) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2064) = 0._r8
         mat(k,2069) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2079) = 0._r8
         mat(k,2083) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2088) = 0._r8
         mat(k,2090) = 0._r8
         mat(k,2091) = 0._r8
         mat(k,2093) = 0._r8
         mat(k,2094) = 0._r8
         mat(k,2104) = 0._r8
         mat(k,2108) = 0._r8
         mat(k,2109) = 0._r8
         mat(k,2110) = 0._r8
         mat(k,2111) = 0._r8
         mat(k,2113) = 0._r8
         mat(k,2117) = 0._r8
         mat(k,2120) = 0._r8
         mat(k,2121) = 0._r8
         mat(k,2122) = 0._r8
         mat(k,2124) = 0._r8
         mat(k,2125) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2154) = 0._r8
         mat(k,2159) = 0._r8
         mat(k,2160) = 0._r8
         mat(k,2163) = 0._r8
         mat(k,2164) = 0._r8
         mat(k,2165) = 0._r8
         mat(k,2167) = 0._r8
         mat(k,2171) = 0._r8
         mat(k,2172) = 0._r8
         mat(k,2173) = 0._r8
         mat(k,2175) = 0._r8
         mat(k,2179) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2195) = 0._r8
         mat(k,2222) = 0._r8
         mat(k,2239) = 0._r8
         mat(k,2241) = 0._r8
         mat(k,2268) = 0._r8
         mat(k,2271) = 0._r8
         mat(k,2274) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2280) = 0._r8
         mat(k,2286) = 0._r8
         mat(k,2291) = 0._r8
         mat(k,2296) = 0._r8
         mat(k,2297) = 0._r8
         mat(k,2302) = 0._r8
         mat(k,2323) = 0._r8
         mat(k,2324) = 0._r8
         mat(k,2328) = 0._r8
         mat(k,2334) = 0._r8
         mat(k,2336) = 0._r8
         mat(k,2339) = 0._r8
         mat(k,2343) = 0._r8
         mat(k,2345) = 0._r8
         mat(k,2347) = 0._r8
         mat(k,2348) = 0._r8
         mat(k,2349) = 0._r8
         mat(k,2351) = 0._r8
         mat(k,2352) = 0._r8
         mat(k,2356) = 0._r8
         mat(k,2357) = 0._r8
         mat(k,2358) = 0._r8
         mat(k,2359) = 0._r8
         mat(k,2360) = 0._r8
         mat(k,2361) = 0._r8
         mat(k,2362) = 0._r8
         mat(k,2363) = 0._r8
         mat(k,2364) = 0._r8
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
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 137) = mat(k, 137) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 141) = mat(k, 141) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 150) = mat(k, 150) - dti(k)
         mat(k, 154) = mat(k, 154) - dti(k)
         mat(k, 158) = mat(k, 158) - dti(k)
         mat(k, 162) = mat(k, 162) - dti(k)
         mat(k, 166) = mat(k, 166) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 180) = mat(k, 180) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 190) = mat(k, 190) - dti(k)
         mat(k, 195) = mat(k, 195) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 203) = mat(k, 203) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 217) = mat(k, 217) - dti(k)
         mat(k, 222) = mat(k, 222) - dti(k)
         mat(k, 224) = mat(k, 224) - dti(k)
         mat(k, 227) = mat(k, 227) - dti(k)
         mat(k, 232) = mat(k, 232) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 261) = mat(k, 261) - dti(k)
         mat(k, 266) = mat(k, 266) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 272) = mat(k, 272) - dti(k)
         mat(k, 275) = mat(k, 275) - dti(k)
         mat(k, 278) = mat(k, 278) - dti(k)
         mat(k, 283) = mat(k, 283) - dti(k)
         mat(k, 286) = mat(k, 286) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 296) = mat(k, 296) - dti(k)
         mat(k, 302) = mat(k, 302) - dti(k)
         mat(k, 306) = mat(k, 306) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 314) = mat(k, 314) - dti(k)
         mat(k, 318) = mat(k, 318) - dti(k)
         mat(k, 321) = mat(k, 321) - dti(k)
         mat(k, 327) = mat(k, 327) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 343) = mat(k, 343) - dti(k)
         mat(k, 349) = mat(k, 349) - dti(k)
         mat(k, 355) = mat(k, 355) - dti(k)
         mat(k, 360) = mat(k, 360) - dti(k)
         mat(k, 365) = mat(k, 365) - dti(k)
         mat(k, 371) = mat(k, 371) - dti(k)
         mat(k, 376) = mat(k, 376) - dti(k)
         mat(k, 381) = mat(k, 381) - dti(k)
         mat(k, 384) = mat(k, 384) - dti(k)
         mat(k, 389) = mat(k, 389) - dti(k)
         mat(k, 394) = mat(k, 394) - dti(k)
         mat(k, 402) = mat(k, 402) - dti(k)
         mat(k, 410) = mat(k, 410) - dti(k)
         mat(k, 418) = mat(k, 418) - dti(k)
         mat(k, 424) = mat(k, 424) - dti(k)
         mat(k, 430) = mat(k, 430) - dti(k)
         mat(k, 436) = mat(k, 436) - dti(k)
         mat(k, 442) = mat(k, 442) - dti(k)
         mat(k, 448) = mat(k, 448) - dti(k)
         mat(k, 454) = mat(k, 454) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 466) = mat(k, 466) - dti(k)
         mat(k, 472) = mat(k, 472) - dti(k)
         mat(k, 478) = mat(k, 478) - dti(k)
         mat(k, 486) = mat(k, 486) - dti(k)
         mat(k, 492) = mat(k, 492) - dti(k)
         mat(k, 499) = mat(k, 499) - dti(k)
         mat(k, 505) = mat(k, 505) - dti(k)
         mat(k, 508) = mat(k, 508) - dti(k)
         mat(k, 515) = mat(k, 515) - dti(k)
         mat(k, 522) = mat(k, 522) - dti(k)
         mat(k, 527) = mat(k, 527) - dti(k)
         mat(k, 536) = mat(k, 536) - dti(k)
         mat(k, 544) = mat(k, 544) - dti(k)
         mat(k, 551) = mat(k, 551) - dti(k)
         mat(k, 556) = mat(k, 556) - dti(k)
         mat(k, 563) = mat(k, 563) - dti(k)
         mat(k, 569) = mat(k, 569) - dti(k)
         mat(k, 577) = mat(k, 577) - dti(k)
         mat(k, 585) = mat(k, 585) - dti(k)
         mat(k, 593) = mat(k, 593) - dti(k)
         mat(k, 601) = mat(k, 601) - dti(k)
         mat(k, 609) = mat(k, 609) - dti(k)
         mat(k, 617) = mat(k, 617) - dti(k)
         mat(k, 626) = mat(k, 626) - dti(k)
         mat(k, 635) = mat(k, 635) - dti(k)
         mat(k, 642) = mat(k, 642) - dti(k)
         mat(k, 646) = mat(k, 646) - dti(k)
         mat(k, 653) = mat(k, 653) - dti(k)
         mat(k, 662) = mat(k, 662) - dti(k)
         mat(k, 670) = mat(k, 670) - dti(k)
         mat(k, 677) = mat(k, 677) - dti(k)
         mat(k, 688) = mat(k, 688) - dti(k)
         mat(k, 699) = mat(k, 699) - dti(k)
         mat(k, 706) = mat(k, 706) - dti(k)
         mat(k, 719) = mat(k, 719) - dti(k)
         mat(k, 729) = mat(k, 729) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 751) = mat(k, 751) - dti(k)
         mat(k, 758) = mat(k, 758) - dti(k)
         mat(k, 769) = mat(k, 769) - dti(k)
         mat(k, 780) = mat(k, 780) - dti(k)
         mat(k, 793) = mat(k, 793) - dti(k)
         mat(k, 804) = mat(k, 804) - dti(k)
         mat(k, 813) = mat(k, 813) - dti(k)
         mat(k, 823) = mat(k, 823) - dti(k)
         mat(k, 831) = mat(k, 831) - dti(k)
         mat(k, 836) = mat(k, 836) - dti(k)
         mat(k, 846) = mat(k, 846) - dti(k)
         mat(k, 855) = mat(k, 855) - dti(k)
         mat(k, 865) = mat(k, 865) - dti(k)
         mat(k, 873) = mat(k, 873) - dti(k)
         mat(k, 879) = mat(k, 879) - dti(k)
         mat(k, 895) = mat(k, 895) - dti(k)
         mat(k, 902) = mat(k, 902) - dti(k)
         mat(k, 909) = mat(k, 909) - dti(k)
         mat(k, 918) = mat(k, 918) - dti(k)
         mat(k, 936) = mat(k, 936) - dti(k)
         mat(k, 960) = mat(k, 960) - dti(k)
         mat(k, 972) = mat(k, 972) - dti(k)
         mat(k, 982) = mat(k, 982) - dti(k)
         mat(k, 989) = mat(k, 989) - dti(k)
         mat(k, 997) = mat(k, 997) - dti(k)
         mat(k,1015) = mat(k,1015) - dti(k)
         mat(k,1035) = mat(k,1035) - dti(k)
         mat(k,1048) = mat(k,1048) - dti(k)
         mat(k,1069) = mat(k,1069) - dti(k)
         mat(k,1081) = mat(k,1081) - dti(k)
         mat(k,1092) = mat(k,1092) - dti(k)
         mat(k,1102) = mat(k,1102) - dti(k)
         mat(k,1116) = mat(k,1116) - dti(k)
         mat(k,1127) = mat(k,1127) - dti(k)
         mat(k,1136) = mat(k,1136) - dti(k)
         mat(k,1149) = mat(k,1149) - dti(k)
         mat(k,1163) = mat(k,1163) - dti(k)
         mat(k,1184) = mat(k,1184) - dti(k)
         mat(k,1200) = mat(k,1200) - dti(k)
         mat(k,1217) = mat(k,1217) - dti(k)
         mat(k,1237) = mat(k,1237) - dti(k)
         mat(k,1253) = mat(k,1253) - dti(k)
         mat(k,1265) = mat(k,1265) - dti(k)
         mat(k,1276) = mat(k,1276) - dti(k)
         mat(k,1301) = mat(k,1301) - dti(k)
         mat(k,1334) = mat(k,1334) - dti(k)
         mat(k,1358) = mat(k,1358) - dti(k)
         mat(k,1379) = mat(k,1379) - dti(k)
         mat(k,1401) = mat(k,1401) - dti(k)
         mat(k,1433) = mat(k,1433) - dti(k)
         mat(k,1448) = mat(k,1448) - dti(k)
         mat(k,1462) = mat(k,1462) - dti(k)
         mat(k,1477) = mat(k,1477) - dti(k)
         mat(k,1494) = mat(k,1494) - dti(k)
         mat(k,1510) = mat(k,1510) - dti(k)
         mat(k,1532) = mat(k,1532) - dti(k)
         mat(k,1576) = mat(k,1576) - dti(k)
         mat(k,1608) = mat(k,1608) - dti(k)
         mat(k,1650) = mat(k,1650) - dti(k)
         mat(k,1823) = mat(k,1823) - dti(k)
         mat(k,1848) = mat(k,1848) - dti(k)
         mat(k,1952) = mat(k,1952) - dti(k)
         mat(k,2005) = mat(k,2005) - dti(k)
         mat(k,2029) = mat(k,2029) - dti(k)
         mat(k,2089) = mat(k,2089) - dti(k)
         mat(k,2129) = mat(k,2129) - dti(k)
         mat(k,2192) = mat(k,2192) - dti(k)
         mat(k,2311) = mat(k,2311) - dti(k)
         mat(k,2338) = mat(k,2338) - dti(k)
         mat(k,2365) = mat(k,2365) - dti(k)
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
