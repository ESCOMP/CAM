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
         mat(k,685) = -(rxt(k,375)*y(k,251))
         mat(k,1853) = -rxt(k,375)*y(k,1)
         mat(k,1561) = rxt(k,378)*y(k,228)
         mat(k,935) = rxt(k,378)*y(k,129)
         mat(k,706) = -(rxt(k,379)*y(k,251))
         mat(k,1855) = -rxt(k,379)*y(k,2)
         mat(k,936) = rxt(k,376)*y(k,240)
         mat(k,2168) = rxt(k,376)*y(k,228)
         mat(k,988) = -(rxt(k,458)*y(k,131) + rxt(k,459)*y(k,139) + rxt(k,460) &
                      *y(k,251))
         mat(k,2043) = -rxt(k,458)*y(k,6)
         mat(k,2289) = -rxt(k,459)*y(k,6)
         mat(k,1882) = -rxt(k,460)*y(k,6)
         mat(k,80) = -(rxt(k,513)*y(k,240) + rxt(k,514)*y(k,129))
         mat(k,2123) = -rxt(k,513)*y(k,7)
         mat(k,1527) = -rxt(k,514)*y(k,7)
         mat(k,980) = rxt(k,516)*y(k,251)
         mat(k,1763) = rxt(k,516)*y(k,6)
         mat(k,207) = -(rxt(k,417)*y(k,251))
         mat(k,1783) = -rxt(k,417)*y(k,8)
         mat(k,111) = -(rxt(k,518)*y(k,240) + rxt(k,519)*y(k,129))
         mat(k,2132) = -rxt(k,518)*y(k,9)
         mat(k,1536) = -rxt(k,519)*y(k,9)
         mat(k,206) = rxt(k,517)*y(k,251)
         mat(k,1773) = rxt(k,517)*y(k,8)
         mat(k,421) = -(rxt(k,420)*y(k,251))
         mat(k,1818) = -rxt(k,420)*y(k,10)
         mat(k,531) = rxt(k,418)*y(k,240)
         mat(k,2145) = rxt(k,418)*y(k,229)
         mat(k,208) = .120_r8*rxt(k,417)*y(k,251)
         mat(k,1784) = .120_r8*rxt(k,417)*y(k,8)
         mat(k,982) = .100_r8*rxt(k,459)*y(k,139)
         mat(k,1031) = .100_r8*rxt(k,462)*y(k,139)
         mat(k,2276) = .100_r8*rxt(k,459)*y(k,6) + .100_r8*rxt(k,462)*y(k,116)
         mat(k,1548) = .500_r8*rxt(k,419)*y(k,229) + .200_r8*rxt(k,446)*y(k,257) &
                      + .060_r8*rxt(k,452)*y(k,259)
         mat(k,532) = .500_r8*rxt(k,419)*y(k,129)
         mat(k,768) = .200_r8*rxt(k,446)*y(k,129)
         mat(k,784) = .060_r8*rxt(k,452)*y(k,129)
         mat(k,1542) = .200_r8*rxt(k,446)*y(k,257) + .200_r8*rxt(k,452)*y(k,259)
         mat(k,767) = .200_r8*rxt(k,446)*y(k,129)
         mat(k,782) = .200_r8*rxt(k,452)*y(k,129)
         mat(k,1558) = .200_r8*rxt(k,446)*y(k,257) + .150_r8*rxt(k,452)*y(k,259)
         mat(k,770) = .200_r8*rxt(k,446)*y(k,129)
         mat(k,785) = .150_r8*rxt(k,452)*y(k,129)
         mat(k,1543) = .210_r8*rxt(k,452)*y(k,259)
         mat(k,783) = .210_r8*rxt(k,452)*y(k,129)
         mat(k,273) = -(rxt(k,380)*y(k,251))
         mat(k,1795) = -rxt(k,380)*y(k,17)
         mat(k,981) = .050_r8*rxt(k,459)*y(k,139)
         mat(k,1030) = .050_r8*rxt(k,462)*y(k,139)
         mat(k,2275) = .050_r8*rxt(k,459)*y(k,6) + .050_r8*rxt(k,462)*y(k,116)
         mat(k,399) = -(rxt(k,346)*y(k,131) + rxt(k,347)*y(k,251))
         mat(k,2036) = -rxt(k,346)*y(k,18)
         mat(k,1815) = -rxt(k,347)*y(k,18)
         mat(k,1455) = -(rxt(k,230)*y(k,44) + rxt(k,231)*y(k,240) + rxt(k,232) &
                      *y(k,139))
         mat(k,2012) = -rxt(k,230)*y(k,19)
         mat(k,2212) = -rxt(k,231)*y(k,19)
         mat(k,2312) = -rxt(k,232)*y(k,19)
         mat(k,2094) = 4.000_r8*rxt(k,233)*y(k,21) + (rxt(k,234)+rxt(k,235))*y(k,61) &
                      + rxt(k,238)*y(k,129) + rxt(k,241)*y(k,138) + rxt(k,488) &
                      *y(k,157) + rxt(k,242)*y(k,251)
         mat(k,180) = rxt(k,220)*y(k,250)
         mat(k,186) = rxt(k,246)*y(k,250)
         mat(k,518) = 2.000_r8*rxt(k,257)*y(k,58) + 2.000_r8*rxt(k,269)*y(k,250) &
                      + 2.000_r8*rxt(k,258)*y(k,251)
         mat(k,615) = rxt(k,259)*y(k,58) + rxt(k,270)*y(k,250) + rxt(k,260)*y(k,251)
         mat(k,458) = 3.000_r8*rxt(k,264)*y(k,58) + 3.000_r8*rxt(k,247)*y(k,250) &
                      + 3.000_r8*rxt(k,265)*y(k,251)
         mat(k,2251) = 2.000_r8*rxt(k,257)*y(k,43) + rxt(k,259)*y(k,45) &
                      + 3.000_r8*rxt(k,264)*y(k,57)
         mat(k,1936) = (rxt(k,234)+rxt(k,235))*y(k,21)
         mat(k,170) = 2.000_r8*rxt(k,248)*y(k,250)
         mat(k,850) = rxt(k,243)*y(k,138) + rxt(k,249)*y(k,250) + rxt(k,244)*y(k,251)
         mat(k,1600) = rxt(k,238)*y(k,21)
         mat(k,1988) = rxt(k,241)*y(k,21) + rxt(k,243)*y(k,83)
         mat(k,1273) = rxt(k,488)*y(k,21)
         mat(k,1737) = rxt(k,220)*y(k,36) + rxt(k,246)*y(k,37) + 2.000_r8*rxt(k,269) &
                      *y(k,43) + rxt(k,270)*y(k,45) + 3.000_r8*rxt(k,247)*y(k,57) &
                      + 2.000_r8*rxt(k,248)*y(k,80) + rxt(k,249)*y(k,83)
         mat(k,1909) = rxt(k,242)*y(k,21) + 2.000_r8*rxt(k,258)*y(k,43) + rxt(k,260) &
                      *y(k,45) + 3.000_r8*rxt(k,265)*y(k,57) + rxt(k,244)*y(k,83)
         mat(k,2088) = rxt(k,236)*y(k,61)
         mat(k,1930) = rxt(k,236)*y(k,21)
         mat(k,1955) = (rxt(k,553)+rxt(k,558))*y(k,93)
         mat(k,817) = (rxt(k,553)+rxt(k,558))*y(k,87)
         mat(k,2107) = -(4._r8*rxt(k,233)*y(k,21) + (rxt(k,234) + rxt(k,235) + rxt(k,236) &
                      ) * y(k,61) + rxt(k,237)*y(k,240) + rxt(k,238)*y(k,129) &
                      + rxt(k,239)*y(k,130) + rxt(k,241)*y(k,138) + rxt(k,242) &
                      *y(k,251) + rxt(k,488)*y(k,157))
         mat(k,1949) = -(rxt(k,234) + rxt(k,235) + rxt(k,236)) * y(k,21)
         mat(k,2226) = -rxt(k,237)*y(k,21)
         mat(k,1614) = -rxt(k,238)*y(k,21)
         mat(k,1710) = -rxt(k,239)*y(k,21)
         mat(k,2002) = -rxt(k,241)*y(k,21)
         mat(k,1923) = -rxt(k,242)*y(k,21)
         mat(k,1280) = -rxt(k,488)*y(k,21)
         mat(k,1461) = rxt(k,232)*y(k,139)
         mat(k,581) = rxt(k,240)*y(k,138)
         mat(k,855) = rxt(k,250)*y(k,250)
         mat(k,823) = rxt(k,245)*y(k,138)
         mat(k,2002) = mat(k,2002) + rxt(k,240)*y(k,22) + rxt(k,245)*y(k,93)
         mat(k,2326) = rxt(k,232)*y(k,19)
         mat(k,1751) = rxt(k,250)*y(k,83)
         mat(k,574) = -(rxt(k,240)*y(k,138))
         mat(k,1978) = -rxt(k,240)*y(k,22)
         mat(k,2090) = rxt(k,239)*y(k,130)
         mat(k,1679) = rxt(k,239)*y(k,21)
         mat(k,282) = -(rxt(k,421)*y(k,251))
         mat(k,1797) = -rxt(k,421)*y(k,24)
         mat(k,1540) = rxt(k,424)*y(k,230)
         mat(k,475) = rxt(k,424)*y(k,129)
         mat(k,373) = -(rxt(k,423)*y(k,251))
         mat(k,1810) = -rxt(k,423)*y(k,25)
         mat(k,476) = rxt(k,422)*y(k,240)
         mat(k,2141) = rxt(k,422)*y(k,230)
         mat(k,324) = -(rxt(k,295)*y(k,58) + rxt(k,296)*y(k,251))
         mat(k,2232) = -rxt(k,295)*y(k,26)
         mat(k,1805) = -rxt(k,296)*y(k,26)
         mat(k,598) = -(rxt(k,297)*y(k,58) + rxt(k,298)*y(k,139) + rxt(k,323)*y(k,251))
         mat(k,2237) = -rxt(k,297)*y(k,27)
         mat(k,2279) = -rxt(k,298)*y(k,27)
         mat(k,1842) = -rxt(k,323)*y(k,27)
         mat(k,303) = -(rxt(k,303)*y(k,251))
         mat(k,1802) = -rxt(k,303)*y(k,28)
         mat(k,857) = .800_r8*rxt(k,299)*y(k,231) + .200_r8*rxt(k,300)*y(k,235)
         mat(k,1619) = .200_r8*rxt(k,300)*y(k,231)
         mat(k,378) = -(rxt(k,304)*y(k,251))
         mat(k,1811) = -rxt(k,304)*y(k,29)
         mat(k,858) = rxt(k,301)*y(k,240)
         mat(k,2142) = rxt(k,301)*y(k,231)
         mat(k,330) = -(rxt(k,305)*y(k,58) + rxt(k,306)*y(k,251))
         mat(k,2233) = -rxt(k,305)*y(k,30)
         mat(k,1806) = -rxt(k,306)*y(k,30)
         mat(k,1064) = -(rxt(k,326)*y(k,131) + rxt(k,327)*y(k,139) + rxt(k,344) &
                      *y(k,251))
         mat(k,2047) = -rxt(k,326)*y(k,31)
         mat(k,2293) = -rxt(k,327)*y(k,31)
         mat(k,1886) = -rxt(k,344)*y(k,31)
         mat(k,876) = .130_r8*rxt(k,404)*y(k,139)
         mat(k,2293) = mat(k,2293) + .130_r8*rxt(k,404)*y(k,100)
         mat(k,427) = -(rxt(k,331)*y(k,251))
         mat(k,1819) = -rxt(k,331)*y(k,32)
         mat(k,830) = rxt(k,329)*y(k,240)
         mat(k,2146) = rxt(k,329)*y(k,232)
         mat(k,148) = -(rxt(k,332)*y(k,251))
         mat(k,1780) = -rxt(k,332)*y(k,33)
         mat(k,307) = -(rxt(k,427)*y(k,251))
         mat(k,1803) = -rxt(k,427)*y(k,34)
         mat(k,676) = rxt(k,425)*y(k,240)
         mat(k,2137) = rxt(k,425)*y(k,233)
         mat(k,139) = -(rxt(k,219)*y(k,250))
         mat(k,1715) = -rxt(k,219)*y(k,35)
         mat(k,178) = -(rxt(k,220)*y(k,250))
         mat(k,1720) = -rxt(k,220)*y(k,36)
         mat(k,183) = -(rxt(k,246)*y(k,250))
         mat(k,1721) = -rxt(k,246)*y(k,37)
         mat(k,152) = -(rxt(k,221)*y(k,250))
         mat(k,1716) = -rxt(k,221)*y(k,38)
         mat(k,188) = -(rxt(k,222)*y(k,250))
         mat(k,1722) = -rxt(k,222)*y(k,39)
         mat(k,156) = -(rxt(k,223)*y(k,250))
         mat(k,1717) = -rxt(k,223)*y(k,40)
         mat(k,193) = -(rxt(k,224)*y(k,250))
         mat(k,1723) = -rxt(k,224)*y(k,41)
         mat(k,160) = -(rxt(k,225)*y(k,250))
         mat(k,1718) = -rxt(k,225)*y(k,42)
         mat(k,517) = -(rxt(k,257)*y(k,58) + rxt(k,258)*y(k,251) + rxt(k,269)*y(k,250))
         mat(k,2236) = -rxt(k,257)*y(k,43)
         mat(k,1832) = -rxt(k,258)*y(k,43)
         mat(k,1732) = -rxt(k,269)*y(k,43)
         mat(k,2024) = -(rxt(k,194)*y(k,58) + rxt(k,230)*y(k,19) + rxt(k,274)*y(k,240) &
                      + rxt(k,275)*y(k,131) + rxt(k,276)*y(k,138) + rxt(k,277) &
                      *y(k,251))
         mat(k,2263) = -rxt(k,194)*y(k,44)
         mat(k,1460) = -rxt(k,230)*y(k,44)
         mat(k,2224) = -rxt(k,274)*y(k,44)
         mat(k,2081) = -rxt(k,275)*y(k,44)
         mat(k,2000) = -rxt(k,276)*y(k,44)
         mat(k,1921) = -rxt(k,277)*y(k,44)
         mat(k,693) = .400_r8*rxt(k,375)*y(k,251)
         mat(k,1004) = .340_r8*rxt(k,459)*y(k,139)
         mat(k,405) = .500_r8*rxt(k,346)*y(k,131)
         mat(k,603) = rxt(k,298)*y(k,139)
         mat(k,1076) = .500_r8*rxt(k,327)*y(k,139)
         mat(k,628) = .500_r8*rxt(k,315)*y(k,251)
         mat(k,828) = rxt(k,282)*y(k,251)
         mat(k,455) = .300_r8*rxt(k,283)*y(k,251)
         mat(k,1478) = (rxt(k,291)+rxt(k,292))*y(k,250)
         mat(k,1947) = rxt(k,201)*y(k,235)
         mat(k,1145) = .800_r8*rxt(k,320)*y(k,251)
         mat(k,886) = .910_r8*rxt(k,404)*y(k,139)
         mat(k,664) = .300_r8*rxt(k,395)*y(k,251)
         mat(k,1244) = .800_r8*rxt(k,399)*y(k,235)
         mat(k,1256) = .120_r8*rxt(k,357)*y(k,139)
         mat(k,636) = .500_r8*rxt(k,370)*y(k,251)
         mat(k,1054) = .340_r8*rxt(k,462)*y(k,139)
         mat(k,1382) = .600_r8*rxt(k,371)*y(k,139)
         mat(k,1612) = .100_r8*rxt(k,377)*y(k,228) + rxt(k,281)*y(k,235) &
                      + .500_r8*rxt(k,348)*y(k,237) + .500_r8*rxt(k,317)*y(k,239) &
                      + .920_r8*rxt(k,387)*y(k,242) + .250_r8*rxt(k,355)*y(k,244) &
                      + rxt(k,364)*y(k,246) + rxt(k,338)*y(k,253) + rxt(k,342) &
                      *y(k,254) + .340_r8*rxt(k,471)*y(k,255) + .320_r8*rxt(k,476) &
                      *y(k,256) + .250_r8*rxt(k,412)*y(k,258)
         mat(k,2081) = mat(k,2081) + .500_r8*rxt(k,346)*y(k,18) + rxt(k,388)*y(k,242) &
                      + .250_r8*rxt(k,354)*y(k,244) + rxt(k,365)*y(k,246)
         mat(k,2324) = .340_r8*rxt(k,459)*y(k,6) + rxt(k,298)*y(k,27) &
                      + .500_r8*rxt(k,327)*y(k,31) + .910_r8*rxt(k,404)*y(k,100) &
                      + .120_r8*rxt(k,357)*y(k,111) + .340_r8*rxt(k,462)*y(k,116) &
                      + .600_r8*rxt(k,371)*y(k,118)
         mat(k,572) = rxt(k,322)*y(k,251)
         mat(k,1122) = .680_r8*rxt(k,480)*y(k,251)
         mat(k,947) = .100_r8*rxt(k,377)*y(k,129)
         mat(k,866) = .700_r8*rxt(k,300)*y(k,235)
         mat(k,838) = rxt(k,328)*y(k,235)
         mat(k,1434) = rxt(k,311)*y(k,235) + rxt(k,384)*y(k,242) + .250_r8*rxt(k,351) &
                      *y(k,244) + rxt(k,360)*y(k,246) + .250_r8*rxt(k,409)*y(k,258)
         mat(k,1664) = rxt(k,201)*y(k,61) + .800_r8*rxt(k,399)*y(k,103) + rxt(k,281) &
                      *y(k,129) + .700_r8*rxt(k,300)*y(k,231) + rxt(k,328)*y(k,232) &
                      + rxt(k,311)*y(k,234) + (4.000_r8*rxt(k,278)+2.000_r8*rxt(k,279)) &
                      *y(k,235) + 1.500_r8*rxt(k,385)*y(k,242) + .750_r8*rxt(k,390) &
                      *y(k,243) + .880_r8*rxt(k,352)*y(k,244) + 2.000_r8*rxt(k,361) &
                      *y(k,246) + .750_r8*rxt(k,464)*y(k,249) + .800_r8*rxt(k,340) &
                      *y(k,254) + .930_r8*rxt(k,469)*y(k,255) + .950_r8*rxt(k,474) &
                      *y(k,256) + .800_r8*rxt(k,410)*y(k,258)
         mat(k,612) = .500_r8*rxt(k,348)*y(k,129)
         mat(k,759) = .500_r8*rxt(k,317)*y(k,129)
         mat(k,2224) = mat(k,2224) + .450_r8*rxt(k,362)*y(k,246) + .150_r8*rxt(k,341) &
                      *y(k,254)
         mat(k,1307) = .920_r8*rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131) + rxt(k,384) &
                      *y(k,234) + 1.500_r8*rxt(k,385)*y(k,235)
         mat(k,1339) = .750_r8*rxt(k,390)*y(k,235)
         mat(k,1360) = .250_r8*rxt(k,355)*y(k,129) + .250_r8*rxt(k,354)*y(k,131) &
                      + .250_r8*rxt(k,351)*y(k,234) + .880_r8*rxt(k,352)*y(k,235)
         mat(k,1402) = rxt(k,364)*y(k,129) + rxt(k,365)*y(k,131) + rxt(k,360)*y(k,234) &
                      + 2.000_r8*rxt(k,361)*y(k,235) + .450_r8*rxt(k,362)*y(k,240) &
                      + 4.000_r8*rxt(k,363)*y(k,246)
         mat(k,1110) = .750_r8*rxt(k,464)*y(k,235)
         mat(k,1749) = (rxt(k,291)+rxt(k,292))*y(k,56)
         mat(k,1921) = mat(k,1921) + .400_r8*rxt(k,375)*y(k,1) + .500_r8*rxt(k,315) &
                      *y(k,53) + rxt(k,282)*y(k,54) + .300_r8*rxt(k,283)*y(k,55) &
                      + .800_r8*rxt(k,320)*y(k,76) + .300_r8*rxt(k,395)*y(k,101) &
                      + .500_r8*rxt(k,370)*y(k,115) + rxt(k,322)*y(k,144) &
                      + .680_r8*rxt(k,480)*y(k,215)
         mat(k,814) = rxt(k,338)*y(k,129)
         mat(k,1177) = rxt(k,342)*y(k,129) + .800_r8*rxt(k,340)*y(k,235) &
                      + .150_r8*rxt(k,341)*y(k,240)
         mat(k,1164) = .340_r8*rxt(k,471)*y(k,129) + .930_r8*rxt(k,469)*y(k,235)
         mat(k,960) = .320_r8*rxt(k,476)*y(k,129) + .950_r8*rxt(k,474)*y(k,235)
         mat(k,1221) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,409)*y(k,234) &
                      + .800_r8*rxt(k,410)*y(k,235)
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
         mat(k,614) = -(rxt(k,259)*y(k,58) + rxt(k,260)*y(k,251) + rxt(k,270)*y(k,250))
         mat(k,2238) = -rxt(k,259)*y(k,45)
         mat(k,1844) = -rxt(k,260)*y(k,45)
         mat(k,1733) = -rxt(k,270)*y(k,45)
         mat(k,164) = -(rxt(k,261)*y(k,251))
         mat(k,1781) = -rxt(k,261)*y(k,46)
         mat(k,1125) = -(rxt(k,307)*y(k,131) + rxt(k,308)*y(k,251))
         mat(k,2051) = -rxt(k,307)*y(k,47)
         mat(k,1890) = -rxt(k,308)*y(k,47)
         mat(k,689) = .800_r8*rxt(k,375)*y(k,251)
         mat(k,402) = rxt(k,346)*y(k,131)
         mat(k,304) = rxt(k,303)*y(k,251)
         mat(k,380) = .500_r8*rxt(k,304)*y(k,251)
         mat(k,1065) = .500_r8*rxt(k,327)*y(k,139)
         mat(k,1367) = .100_r8*rxt(k,371)*y(k,139)
         mat(k,1583) = .400_r8*rxt(k,377)*y(k,228) + rxt(k,302)*y(k,231) &
                      + .270_r8*rxt(k,330)*y(k,232) + rxt(k,348)*y(k,237) + rxt(k,367) &
                      *y(k,248) + rxt(k,338)*y(k,253)
         mat(k,2051) = mat(k,2051) + rxt(k,346)*y(k,18)
         mat(k,2296) = .500_r8*rxt(k,327)*y(k,31) + .100_r8*rxt(k,371)*y(k,118)
         mat(k,941) = .400_r8*rxt(k,377)*y(k,129)
         mat(k,861) = rxt(k,302)*y(k,129) + 3.200_r8*rxt(k,299)*y(k,231) &
                      + .800_r8*rxt(k,300)*y(k,235)
         mat(k,833) = .270_r8*rxt(k,330)*y(k,129)
         mat(k,1637) = .800_r8*rxt(k,300)*y(k,231)
         mat(k,608) = rxt(k,348)*y(k,129)
         mat(k,2195) = .200_r8*rxt(k,366)*y(k,248)
         mat(k,718) = rxt(k,367)*y(k,129) + .200_r8*rxt(k,366)*y(k,240)
         mat(k,1890) = mat(k,1890) + .800_r8*rxt(k,375)*y(k,1) + rxt(k,303)*y(k,28) &
                      + .500_r8*rxt(k,304)*y(k,29)
         mat(k,809) = rxt(k,338)*y(k,129)
         mat(k,407) = -(rxt(k,262)*y(k,58) + rxt(k,263)*y(k,251))
         mat(k,2234) = -rxt(k,262)*y(k,48)
         mat(k,1816) = -rxt(k,263)*y(k,48)
         mat(k,142) = -(rxt(k,309)*y(k,251))
         mat(k,1779) = -rxt(k,309)*y(k,49)
         mat(k,962) = -(rxt(k,345)*y(k,251))
         mat(k,1880) = -rxt(k,345)*y(k,50)
         mat(k,688) = .800_r8*rxt(k,375)*y(k,251)
         mat(k,986) = .520_r8*rxt(k,459)*y(k,139)
         mat(k,401) = .500_r8*rxt(k,346)*y(k,131)
         mat(k,1035) = .520_r8*rxt(k,462)*y(k,139)
         mat(k,1576) = .250_r8*rxt(k,377)*y(k,228) + .820_r8*rxt(k,330)*y(k,232) &
                      + .500_r8*rxt(k,348)*y(k,237) + .270_r8*rxt(k,471)*y(k,255) &
                      + .040_r8*rxt(k,476)*y(k,256)
         mat(k,2041) = .500_r8*rxt(k,346)*y(k,18)
         mat(k,2287) = .520_r8*rxt(k,459)*y(k,6) + .520_r8*rxt(k,462)*y(k,116)
         mat(k,1114) = .500_r8*rxt(k,480)*y(k,251)
         mat(k,940) = .250_r8*rxt(k,377)*y(k,129)
         mat(k,832) = .820_r8*rxt(k,330)*y(k,129) + .820_r8*rxt(k,328)*y(k,235)
         mat(k,1631) = .820_r8*rxt(k,328)*y(k,232) + .150_r8*rxt(k,469)*y(k,255) &
                      + .025_r8*rxt(k,474)*y(k,256)
         mat(k,607) = .500_r8*rxt(k,348)*y(k,129)
         mat(k,1880) = mat(k,1880) + .800_r8*rxt(k,375)*y(k,1) + .500_r8*rxt(k,480) &
                      *y(k,215)
         mat(k,1151) = .270_r8*rxt(k,471)*y(k,129) + .150_r8*rxt(k,469)*y(k,235)
         mat(k,953) = .040_r8*rxt(k,476)*y(k,129) + .025_r8*rxt(k,474)*y(k,235)
         mat(k,1261) = -(rxt(k,333)*y(k,131) + rxt(k,334)*y(k,251))
         mat(k,2061) = -rxt(k,333)*y(k,51)
         mat(k,1900) = -rxt(k,334)*y(k,51)
         mat(k,1181) = rxt(k,335)*y(k,251)
         mat(k,1250) = .880_r8*rxt(k,357)*y(k,139)
         mat(k,1370) = .500_r8*rxt(k,371)*y(k,139)
         mat(k,1593) = .170_r8*rxt(k,430)*y(k,236) + .050_r8*rxt(k,393)*y(k,243) &
                      + .250_r8*rxt(k,355)*y(k,244) + .170_r8*rxt(k,436)*y(k,247) &
                      + .400_r8*rxt(k,446)*y(k,257) + .250_r8*rxt(k,412)*y(k,258) &
                      + .540_r8*rxt(k,452)*y(k,259) + .510_r8*rxt(k,455)*y(k,260)
         mat(k,2061) = mat(k,2061) + .050_r8*rxt(k,394)*y(k,243) + .250_r8*rxt(k,354) &
                      *y(k,244) + .250_r8*rxt(k,413)*y(k,258)
         mat(k,891) = rxt(k,336)*y(k,251)
         mat(k,2304) = .880_r8*rxt(k,357)*y(k,111) + .500_r8*rxt(k,371)*y(k,118)
         mat(k,1420) = .250_r8*rxt(k,351)*y(k,244) + .250_r8*rxt(k,409)*y(k,258)
         mat(k,1646) = .240_r8*rxt(k,352)*y(k,244) + .500_r8*rxt(k,340)*y(k,254) &
                      + .100_r8*rxt(k,410)*y(k,258)
         mat(k,801) = .170_r8*rxt(k,430)*y(k,129) + .070_r8*rxt(k,429)*y(k,240)
         mat(k,2204) = .070_r8*rxt(k,429)*y(k,236) + .070_r8*rxt(k,435)*y(k,247)
         mat(k,1327) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,1350) = .250_r8*rxt(k,355)*y(k,129) + .250_r8*rxt(k,354)*y(k,131) &
                      + .250_r8*rxt(k,351)*y(k,234) + .240_r8*rxt(k,352)*y(k,235)
         mat(k,906) = .170_r8*rxt(k,436)*y(k,129) + .070_r8*rxt(k,435)*y(k,240)
         mat(k,1900) = mat(k,1900) + rxt(k,335)*y(k,97) + rxt(k,336)*y(k,132)
         mat(k,1171) = .500_r8*rxt(k,340)*y(k,235)
         mat(k,777) = .400_r8*rxt(k,446)*y(k,129)
         mat(k,1214) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,234) + .100_r8*rxt(k,410)*y(k,235)
         mat(k,793) = .540_r8*rxt(k,452)*y(k,129)
         mat(k,543) = .510_r8*rxt(k,455)*y(k,129)
         mat(k,724) = -(rxt(k,314)*y(k,251))
         mat(k,1857) = -rxt(k,314)*y(k,52)
         mat(k,1059) = .120_r8*rxt(k,327)*y(k,139)
         mat(k,2281) = .120_r8*rxt(k,327)*y(k,31)
         mat(k,1410) = .100_r8*rxt(k,311)*y(k,235) + .150_r8*rxt(k,312)*y(k,240)
         mat(k,1624) = .100_r8*rxt(k,311)*y(k,234)
         mat(k,2170) = .150_r8*rxt(k,312)*y(k,234) + .150_r8*rxt(k,362)*y(k,246)
         mat(k,1389) = .150_r8*rxt(k,362)*y(k,240)
         mat(k,623) = -(rxt(k,315)*y(k,251))
         mat(k,1845) = -rxt(k,315)*y(k,53)
         mat(k,1409) = .400_r8*rxt(k,312)*y(k,240)
         mat(k,2162) = .400_r8*rxt(k,312)*y(k,234) + .400_r8*rxt(k,362)*y(k,246)
         mat(k,1387) = .400_r8*rxt(k,362)*y(k,240)
         mat(k,826) = -(rxt(k,282)*y(k,251))
         mat(k,1866) = -rxt(k,282)*y(k,54)
         mat(k,1227) = .200_r8*rxt(k,399)*y(k,235)
         mat(k,859) = .300_r8*rxt(k,300)*y(k,235)
         mat(k,1625) = .200_r8*rxt(k,399)*y(k,103) + .300_r8*rxt(k,300)*y(k,231) &
                      + 2.000_r8*rxt(k,279)*y(k,235) + .250_r8*rxt(k,385)*y(k,242) &
                      + .250_r8*rxt(k,390)*y(k,243) + .250_r8*rxt(k,352)*y(k,244) &
                      + .250_r8*rxt(k,464)*y(k,249) + .500_r8*rxt(k,340)*y(k,254) &
                      + .250_r8*rxt(k,469)*y(k,255) + .250_r8*rxt(k,474)*y(k,256) &
                      + .300_r8*rxt(k,410)*y(k,258)
         mat(k,1287) = .250_r8*rxt(k,385)*y(k,235)
         mat(k,1317) = .250_r8*rxt(k,390)*y(k,235)
         mat(k,1345) = .250_r8*rxt(k,352)*y(k,235)
         mat(k,1099) = .250_r8*rxt(k,464)*y(k,235)
         mat(k,1168) = .500_r8*rxt(k,340)*y(k,235)
         mat(k,1149) = .250_r8*rxt(k,469)*y(k,235)
         mat(k,951) = .250_r8*rxt(k,474)*y(k,235)
         mat(k,1207) = .300_r8*rxt(k,410)*y(k,235)
         mat(k,451) = -(rxt(k,283)*y(k,251))
         mat(k,1823) = -rxt(k,283)*y(k,55)
         mat(k,1622) = rxt(k,280)*y(k,240)
         mat(k,2149) = rxt(k,280)*y(k,235)
         mat(k,1470) = -(rxt(k,195)*y(k,58) + rxt(k,251)*y(k,75) + rxt(k,284)*y(k,251) &
                      + (rxt(k,290) + rxt(k,291) + rxt(k,292)) * y(k,250))
         mat(k,2252) = -rxt(k,195)*y(k,56)
         mat(k,915) = -rxt(k,251)*y(k,56)
         mat(k,1910) = -rxt(k,284)*y(k,56)
         mat(k,1738) = -(rxt(k,290) + rxt(k,291) + rxt(k,292)) * y(k,56)
         mat(k,1070) = .100_r8*rxt(k,327)*y(k,139)
         mat(k,2313) = .100_r8*rxt(k,327)*y(k,31)
         mat(k,457) = -(rxt(k,247)*y(k,250) + rxt(k,264)*y(k,58) + rxt(k,265)*y(k,251))
         mat(k,1731) = -rxt(k,247)*y(k,57)
         mat(k,2235) = -rxt(k,264)*y(k,57)
         mat(k,1824) = -rxt(k,265)*y(k,57)
         mat(k,2267) = -(rxt(k,194)*y(k,44) + rxt(k,195)*y(k,56) + rxt(k,196)*y(k,79) &
                      + rxt(k,197)*y(k,81) + (rxt(k,198) + rxt(k,199)) * y(k,240) &
                      + rxt(k,200)*y(k,139) + rxt(k,207)*y(k,62) + rxt(k,216)*y(k,94) &
                      + rxt(k,257)*y(k,43) + rxt(k,259)*y(k,45) + rxt(k,262)*y(k,48) &
                      + rxt(k,264)*y(k,57) + rxt(k,305)*y(k,30))
         mat(k,2028) = -rxt(k,194)*y(k,58)
         mat(k,1481) = -rxt(k,195)*y(k,58)
         mat(k,1451) = -rxt(k,196)*y(k,58)
         mat(k,644) = -rxt(k,197)*y(k,58)
         mat(k,2228) = -(rxt(k,198) + rxt(k,199)) * y(k,58)
         mat(k,2328) = -rxt(k,200)*y(k,58)
         mat(k,932) = -rxt(k,207)*y(k,58)
         mat(k,847) = -rxt(k,216)*y(k,58)
         mat(k,522) = -rxt(k,257)*y(k,58)
         mat(k,621) = -rxt(k,259)*y(k,58)
         mat(k,413) = -rxt(k,262)*y(k,58)
         mat(k,462) = -rxt(k,264)*y(k,58)
         mat(k,334) = -rxt(k,305)*y(k,58)
         mat(k,2109) = rxt(k,235)*y(k,61)
         mat(k,141) = 4.000_r8*rxt(k,219)*y(k,250)
         mat(k,182) = rxt(k,220)*y(k,250)
         mat(k,155) = 2.000_r8*rxt(k,221)*y(k,250)
         mat(k,192) = 2.000_r8*rxt(k,222)*y(k,250)
         mat(k,159) = 2.000_r8*rxt(k,223)*y(k,250)
         mat(k,197) = rxt(k,224)*y(k,250)
         mat(k,163) = 2.000_r8*rxt(k,225)*y(k,250)
         mat(k,166) = 3.000_r8*rxt(k,261)*y(k,251)
         mat(k,413) = mat(k,413) + rxt(k,263)*y(k,251)
         mat(k,1951) = rxt(k,235)*y(k,21) + (4.000_r8*rxt(k,202)+2.000_r8*rxt(k,204)) &
                      *y(k,61) + rxt(k,206)*y(k,129) + rxt(k,211)*y(k,138) &
                      + rxt(k,489)*y(k,157) + rxt(k,201)*y(k,235) + rxt(k,212) &
                      *y(k,251)
         mat(k,272) = rxt(k,256)*y(k,250)
         mat(k,268) = rxt(k,271)*y(k,250) + rxt(k,266)*y(k,251)
         mat(k,289) = rxt(k,272)*y(k,250) + rxt(k,267)*y(k,251)
         mat(k,344) = rxt(k,273)*y(k,250) + rxt(k,268)*y(k,251)
         mat(k,1973) = rxt(k,214)*y(k,138) + rxt(k,226)*y(k,250) + rxt(k,215)*y(k,251)
         mat(k,1616) = rxt(k,206)*y(k,61)
         mat(k,2004) = rxt(k,211)*y(k,61) + rxt(k,214)*y(k,87)
         mat(k,1282) = rxt(k,489)*y(k,61)
         mat(k,1668) = rxt(k,201)*y(k,61)
         mat(k,1753) = 4.000_r8*rxt(k,219)*y(k,35) + rxt(k,220)*y(k,36) &
                      + 2.000_r8*rxt(k,221)*y(k,38) + 2.000_r8*rxt(k,222)*y(k,39) &
                      + 2.000_r8*rxt(k,223)*y(k,40) + rxt(k,224)*y(k,41) &
                      + 2.000_r8*rxt(k,225)*y(k,42) + rxt(k,256)*y(k,67) + rxt(k,271) &
                      *y(k,84) + rxt(k,272)*y(k,85) + rxt(k,273)*y(k,86) + rxt(k,226) &
                      *y(k,87)
         mat(k,1925) = 3.000_r8*rxt(k,261)*y(k,46) + rxt(k,263)*y(k,48) + rxt(k,212) &
                      *y(k,61) + rxt(k,266)*y(k,84) + rxt(k,267)*y(k,85) + rxt(k,268) &
                      *y(k,86) + rxt(k,215)*y(k,87)
         mat(k,2231) = rxt(k,207)*y(k,62)
         mat(k,1929) = 2.000_r8*rxt(k,203)*y(k,61)
         mat(k,922) = rxt(k,207)*y(k,58) + (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,87)
         mat(k,1954) = (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,62) + (rxt(k,546) &
                       +rxt(k,552)+rxt(k,557))*y(k,94)
         mat(k,841) = (rxt(k,546)+rxt(k,552)+rxt(k,557))*y(k,87)
         mat(k,1928) = 2.000_r8*rxt(k,228)*y(k,61)
         mat(k,1944) = -(rxt(k,201)*y(k,235) + (4._r8*rxt(k,202) + 4._r8*rxt(k,203) &
                      + 4._r8*rxt(k,204) + 4._r8*rxt(k,228)) * y(k,61) + rxt(k,205) &
                      *y(k,240) + rxt(k,206)*y(k,129) + rxt(k,208)*y(k,130) + rxt(k,211) &
                      *y(k,138) + (rxt(k,212) + rxt(k,213)) * y(k,251) + (rxt(k,234) &
                      + rxt(k,235) + rxt(k,236)) * y(k,21) + rxt(k,489)*y(k,157))
         mat(k,1661) = -rxt(k,201)*y(k,61)
         mat(k,2221) = -rxt(k,205)*y(k,61)
         mat(k,1609) = -rxt(k,206)*y(k,61)
         mat(k,1705) = -rxt(k,208)*y(k,61)
         mat(k,1997) = -rxt(k,211)*y(k,61)
         mat(k,1918) = -(rxt(k,212) + rxt(k,213)) * y(k,61)
         mat(k,2102) = -(rxt(k,234) + rxt(k,235) + rxt(k,236)) * y(k,61)
         mat(k,1278) = -rxt(k,489)*y(k,61)
         mat(k,2260) = rxt(k,216)*y(k,94) + rxt(k,200)*y(k,139) + rxt(k,199)*y(k,240)
         mat(k,928) = rxt(k,209)*y(k,138)
         mat(k,1966) = rxt(k,227)*y(k,250)
         mat(k,844) = rxt(k,216)*y(k,58) + rxt(k,217)*y(k,138) + rxt(k,218)*y(k,251)
         mat(k,1997) = mat(k,1997) + rxt(k,209)*y(k,62) + rxt(k,217)*y(k,94)
         mat(k,2321) = rxt(k,200)*y(k,58)
         mat(k,365) = rxt(k,494)*y(k,157)
         mat(k,1278) = mat(k,1278) + rxt(k,494)*y(k,141)
         mat(k,2221) = mat(k,2221) + rxt(k,199)*y(k,58)
         mat(k,1746) = rxt(k,227)*y(k,87)
         mat(k,1918) = mat(k,1918) + rxt(k,218)*y(k,94)
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
         mat(k,924) = -(rxt(k,207)*y(k,58) + rxt(k,209)*y(k,138) + rxt(k,210)*y(k,251) &
                      + (rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,87))
         mat(k,2243) = -rxt(k,207)*y(k,62)
         mat(k,1984) = -rxt(k,209)*y(k,62)
         mat(k,1877) = -rxt(k,210)*y(k,62)
         mat(k,1958) = -(rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,62)
         mat(k,1934) = rxt(k,208)*y(k,130)
         mat(k,1688) = rxt(k,208)*y(k,61)
         mat(k,1135) = -(rxt(k,294)*y(k,251))
         mat(k,1891) = -rxt(k,294)*y(k,64)
         mat(k,994) = .230_r8*rxt(k,459)*y(k,139)
         mat(k,1454) = rxt(k,230)*y(k,44)
         mat(k,327) = .350_r8*rxt(k,296)*y(k,251)
         mat(k,601) = .630_r8*rxt(k,298)*y(k,139)
         mat(k,1066) = .560_r8*rxt(k,327)*y(k,139)
         mat(k,2010) = rxt(k,230)*y(k,19) + rxt(k,194)*y(k,58) + rxt(k,275)*y(k,131) &
                      + rxt(k,276)*y(k,138) + rxt(k,277)*y(k,251)
         mat(k,408) = rxt(k,262)*y(k,58)
         mat(k,1260) = rxt(k,333)*y(k,131) + rxt(k,334)*y(k,251)
         mat(k,2247) = rxt(k,194)*y(k,44) + rxt(k,262)*y(k,48)
         mat(k,971) = rxt(k,321)*y(k,251)
         mat(k,877) = .620_r8*rxt(k,404)*y(k,139)
         mat(k,1248) = .650_r8*rxt(k,357)*y(k,139)
         mat(k,1043) = .230_r8*rxt(k,462)*y(k,139)
         mat(k,1368) = .560_r8*rxt(k,371)*y(k,139)
         mat(k,1584) = .170_r8*rxt(k,430)*y(k,236) + .220_r8*rxt(k,355)*y(k,244) &
                      + .400_r8*rxt(k,433)*y(k,245) + .350_r8*rxt(k,436)*y(k,247) &
                      + .225_r8*rxt(k,471)*y(k,255) + .250_r8*rxt(k,412)*y(k,258)
         mat(k,2052) = rxt(k,275)*y(k,44) + rxt(k,333)*y(k,51) + .220_r8*rxt(k,354) &
                      *y(k,244) + .500_r8*rxt(k,413)*y(k,258)
         mat(k,1985) = rxt(k,276)*y(k,44) + rxt(k,483)*y(k,142)
         mat(k,2297) = .230_r8*rxt(k,459)*y(k,6) + .630_r8*rxt(k,298)*y(k,27) &
                      + .560_r8*rxt(k,327)*y(k,31) + .620_r8*rxt(k,404)*y(k,100) &
                      + .650_r8*rxt(k,357)*y(k,111) + .230_r8*rxt(k,462)*y(k,116) &
                      + .560_r8*rxt(k,371)*y(k,118)
         mat(k,394) = rxt(k,483)*y(k,138) + rxt(k,484)*y(k,251)
         mat(k,1116) = .700_r8*rxt(k,480)*y(k,251)
         mat(k,1414) = .220_r8*rxt(k,351)*y(k,244) + .250_r8*rxt(k,409)*y(k,258)
         mat(k,1638) = .110_r8*rxt(k,352)*y(k,244) + .125_r8*rxt(k,469)*y(k,255) &
                      + .200_r8*rxt(k,410)*y(k,258)
         mat(k,800) = .170_r8*rxt(k,430)*y(k,129) + .070_r8*rxt(k,429)*y(k,240)
         mat(k,2196) = .070_r8*rxt(k,429)*y(k,236) + .160_r8*rxt(k,432)*y(k,245) &
                      + .140_r8*rxt(k,435)*y(k,247)
         mat(k,1346) = .220_r8*rxt(k,355)*y(k,129) + .220_r8*rxt(k,354)*y(k,131) &
                      + .220_r8*rxt(k,351)*y(k,234) + .110_r8*rxt(k,352)*y(k,235)
         mat(k,763) = .400_r8*rxt(k,433)*y(k,129) + .160_r8*rxt(k,432)*y(k,240)
         mat(k,905) = .350_r8*rxt(k,436)*y(k,129) + .140_r8*rxt(k,435)*y(k,240)
         mat(k,1891) = mat(k,1891) + .350_r8*rxt(k,296)*y(k,26) + rxt(k,277)*y(k,44) &
                      + rxt(k,334)*y(k,51) + rxt(k,321)*y(k,77) + rxt(k,484)*y(k,142) &
                      + .700_r8*rxt(k,480)*y(k,215)
         mat(k,1153) = .225_r8*rxt(k,471)*y(k,129) + .125_r8*rxt(k,469)*y(k,235)
         mat(k,1210) = .250_r8*rxt(k,412)*y(k,129) + .500_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,234) + .200_r8*rxt(k,410)*y(k,235)
         mat(k,983) = .270_r8*rxt(k,459)*y(k,139)
         mat(k,1061) = .200_r8*rxt(k,327)*y(k,139)
         mat(k,725) = rxt(k,314)*y(k,251)
         mat(k,624) = .500_r8*rxt(k,315)*y(k,251)
         mat(k,1134) = rxt(k,294)*y(k,251)
         mat(k,1139) = .800_r8*rxt(k,320)*y(k,251)
         mat(k,969) = rxt(k,321)*y(k,251)
         mat(k,1009) = rxt(k,286)*y(k,251)
         mat(k,631) = .500_r8*rxt(k,370)*y(k,251)
         mat(k,1032) = .270_r8*rxt(k,462)*y(k,139)
         mat(k,1364) = .100_r8*rxt(k,371)*y(k,139)
         mat(k,1571) = rxt(k,313)*y(k,234) + .900_r8*rxt(k,471)*y(k,255)
         mat(k,2283) = .270_r8*rxt(k,459)*y(k,6) + .200_r8*rxt(k,327)*y(k,31) &
                      + .270_r8*rxt(k,462)*y(k,116) + .100_r8*rxt(k,371)*y(k,118)
         mat(k,1113) = 1.800_r8*rxt(k,480)*y(k,251)
         mat(k,1411) = rxt(k,313)*y(k,129) + 4.000_r8*rxt(k,310)*y(k,234) &
                      + .900_r8*rxt(k,311)*y(k,235) + rxt(k,384)*y(k,242) &
                      + 2.000_r8*rxt(k,360)*y(k,246) + rxt(k,409)*y(k,258)
         mat(k,1628) = .900_r8*rxt(k,311)*y(k,234) + rxt(k,361)*y(k,246) &
                      + .500_r8*rxt(k,469)*y(k,255)
         mat(k,2184) = .450_r8*rxt(k,362)*y(k,246)
         mat(k,1288) = rxt(k,384)*y(k,234)
         mat(k,1390) = 2.000_r8*rxt(k,360)*y(k,234) + rxt(k,361)*y(k,235) &
                      + .450_r8*rxt(k,362)*y(k,240) + 4.000_r8*rxt(k,363)*y(k,246)
         mat(k,1871) = rxt(k,314)*y(k,52) + .500_r8*rxt(k,315)*y(k,53) + rxt(k,294) &
                      *y(k,64) + .800_r8*rxt(k,320)*y(k,76) + rxt(k,321)*y(k,77) &
                      + rxt(k,286)*y(k,89) + .500_r8*rxt(k,370)*y(k,115) &
                      + 1.800_r8*rxt(k,480)*y(k,215)
         mat(k,1150) = .900_r8*rxt(k,471)*y(k,129) + .500_r8*rxt(k,469)*y(k,235)
         mat(k,1208) = rxt(k,409)*y(k,234)
         mat(k,290) = -(rxt(k,255)*y(k,250))
         mat(k,1729) = -rxt(k,255)*y(k,66)
         mat(k,179) = rxt(k,220)*y(k,250)
         mat(k,184) = rxt(k,246)*y(k,250)
         mat(k,190) = rxt(k,222)*y(k,250)
         mat(k,157) = 2.000_r8*rxt(k,223)*y(k,250)
         mat(k,194) = 2.000_r8*rxt(k,224)*y(k,250)
         mat(k,161) = rxt(k,225)*y(k,250)
         mat(k,169) = 2.000_r8*rxt(k,248)*y(k,250)
         mat(k,286) = rxt(k,272)*y(k,250) + rxt(k,267)*y(k,251)
         mat(k,339) = rxt(k,273)*y(k,250) + rxt(k,268)*y(k,251)
         mat(k,1729) = mat(k,1729) + rxt(k,220)*y(k,36) + rxt(k,246)*y(k,37) &
                      + rxt(k,222)*y(k,39) + 2.000_r8*rxt(k,223)*y(k,40) &
                      + 2.000_r8*rxt(k,224)*y(k,41) + rxt(k,225)*y(k,42) &
                      + 2.000_r8*rxt(k,248)*y(k,80) + rxt(k,272)*y(k,85) + rxt(k,273) &
                      *y(k,86)
         mat(k,1799) = rxt(k,267)*y(k,85) + rxt(k,268)*y(k,86)
         mat(k,269) = -(rxt(k,256)*y(k,250))
         mat(k,1727) = -rxt(k,256)*y(k,67)
         mat(k,153) = rxt(k,221)*y(k,250)
         mat(k,189) = rxt(k,222)*y(k,250)
         mat(k,265) = rxt(k,271)*y(k,250) + rxt(k,266)*y(k,251)
         mat(k,1727) = mat(k,1727) + rxt(k,221)*y(k,38) + rxt(k,222)*y(k,39) &
                      + rxt(k,271)*y(k,84)
         mat(k,1794) = rxt(k,266)*y(k,84)
         mat(k,237) = -(rxt(k,428)*y(k,251))
         mat(k,1788) = -rxt(k,428)*y(k,68)
         mat(k,231) = .180_r8*rxt(k,448)*y(k,251)
         mat(k,1788) = mat(k,1788) + .180_r8*rxt(k,448)*y(k,217)
         mat(k,318) = -(rxt(k,481)*y(k,131) + (rxt(k,482) + rxt(k,496)) * y(k,251))
         mat(k,2033) = -rxt(k,481)*y(k,69)
         mat(k,1804) = -(rxt(k,482) + rxt(k,496)) * y(k,69)
         mat(k,752) = rxt(k,316)*y(k,240)
         mat(k,2135) = rxt(k,316)*y(k,239)
         mat(k,913) = -(rxt(k,251)*y(k,56) + rxt(k,252)*y(k,79) + rxt(k,253)*y(k,261) &
                      + rxt(k,254)*y(k,91))
         mat(k,1467) = -rxt(k,251)*y(k,75)
         mat(k,1440) = -rxt(k,252)*y(k,75)
         mat(k,2335) = -rxt(k,253)*y(k,75)
         mat(k,1484) = -rxt(k,254)*y(k,75)
         mat(k,185) = rxt(k,246)*y(k,250)
         mat(k,195) = rxt(k,224)*y(k,250)
         mat(k,291) = 2.000_r8*rxt(k,255)*y(k,250)
         mat(k,270) = rxt(k,256)*y(k,250)
         mat(k,1735) = rxt(k,246)*y(k,37) + rxt(k,224)*y(k,41) + 2.000_r8*rxt(k,255) &
                      *y(k,66) + rxt(k,256)*y(k,67)
         mat(k,1142) = -(rxt(k,320)*y(k,251))
         mat(k,1892) = -rxt(k,320)*y(k,76)
         mat(k,658) = .700_r8*rxt(k,395)*y(k,251)
         mat(k,584) = .500_r8*rxt(k,396)*y(k,251)
         mat(k,441) = rxt(k,407)*y(k,251)
         mat(k,1585) = .050_r8*rxt(k,393)*y(k,243) + .530_r8*rxt(k,355)*y(k,244) &
                      + .225_r8*rxt(k,471)*y(k,255) + .250_r8*rxt(k,412)*y(k,258)
         mat(k,2053) = .050_r8*rxt(k,394)*y(k,243) + .530_r8*rxt(k,354)*y(k,244) &
                      + .250_r8*rxt(k,413)*y(k,258)
         mat(k,1415) = .530_r8*rxt(k,351)*y(k,244) + .250_r8*rxt(k,409)*y(k,258)
         mat(k,1639) = .260_r8*rxt(k,352)*y(k,244) + .125_r8*rxt(k,469)*y(k,255) &
                      + .100_r8*rxt(k,410)*y(k,258)
         mat(k,1322) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,1347) = .530_r8*rxt(k,355)*y(k,129) + .530_r8*rxt(k,354)*y(k,131) &
                      + .530_r8*rxt(k,351)*y(k,234) + .260_r8*rxt(k,352)*y(k,235)
         mat(k,1892) = mat(k,1892) + .700_r8*rxt(k,395)*y(k,101) + .500_r8*rxt(k,396) &
                      *y(k,102) + rxt(k,407)*y(k,122)
         mat(k,1154) = .225_r8*rxt(k,471)*y(k,129) + .125_r8*rxt(k,469)*y(k,235)
         mat(k,1211) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,234) + .100_r8*rxt(k,410)*y(k,235)
         mat(k,970) = -(rxt(k,321)*y(k,251))
         mat(k,1881) = -rxt(k,321)*y(k,77)
         mat(k,325) = .650_r8*rxt(k,296)*y(k,251)
         mat(k,1140) = .200_r8*rxt(k,320)*y(k,251)
         mat(k,1084) = rxt(k,408)*y(k,251)
         mat(k,1577) = rxt(k,419)*y(k,229) + .050_r8*rxt(k,393)*y(k,243) &
                      + .400_r8*rxt(k,433)*y(k,245) + .170_r8*rxt(k,436)*y(k,247) &
                      + .700_r8*rxt(k,439)*y(k,252) + .600_r8*rxt(k,446)*y(k,257) &
                      + .250_r8*rxt(k,412)*y(k,258) + .340_r8*rxt(k,452)*y(k,259) &
                      + .170_r8*rxt(k,455)*y(k,260)
         mat(k,2042) = .050_r8*rxt(k,394)*y(k,243) + .250_r8*rxt(k,413)*y(k,258)
         mat(k,535) = rxt(k,419)*y(k,129)
         mat(k,1412) = .250_r8*rxt(k,409)*y(k,258)
         mat(k,1632) = .100_r8*rxt(k,410)*y(k,258)
         mat(k,2190) = .160_r8*rxt(k,432)*y(k,245) + .070_r8*rxt(k,435)*y(k,247)
         mat(k,1320) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,762) = .400_r8*rxt(k,433)*y(k,129) + .160_r8*rxt(k,432)*y(k,240)
         mat(k,904) = .170_r8*rxt(k,436)*y(k,129) + .070_r8*rxt(k,435)*y(k,240)
         mat(k,1881) = mat(k,1881) + .650_r8*rxt(k,296)*y(k,26) + .200_r8*rxt(k,320) &
                      *y(k,76) + rxt(k,408)*y(k,123)
         mat(k,497) = .700_r8*rxt(k,439)*y(k,129)
         mat(k,775) = .600_r8*rxt(k,446)*y(k,129)
         mat(k,1209) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,234) + .100_r8*rxt(k,410)*y(k,235)
         mat(k,791) = .340_r8*rxt(k,452)*y(k,129)
         mat(k,542) = .170_r8*rxt(k,455)*y(k,129)
         mat(k,1503) = -((rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,240) + rxt(k,160) &
                      *y(k,139))
         mat(k,2215) = -(rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,78)
         mat(k,2315) = -rxt(k,160)*y(k,78)
         mat(k,2015) = rxt(k,277)*y(k,251)
         mat(k,1472) = rxt(k,291)*y(k,250)
         mat(k,2254) = rxt(k,196)*y(k,79)
         mat(k,917) = rxt(k,252)*y(k,79)
         mat(k,1444) = rxt(k,196)*y(k,58) + rxt(k,252)*y(k,75) + rxt(k,152)*y(k,138) &
                      + rxt(k,144)*y(k,250) + rxt(k,161)*y(k,251)
         mat(k,851) = rxt(k,250)*y(k,250)
         mat(k,1961) = rxt(k,227)*y(k,250)
         mat(k,507) = rxt(k,182)*y(k,251)
         mat(k,1991) = rxt(k,152)*y(k,79) + rxt(k,164)*y(k,251)
         mat(k,396) = rxt(k,484)*y(k,251)
         mat(k,562) = rxt(k,490)*y(k,251)
         mat(k,1274) = rxt(k,495)*y(k,251)
         mat(k,1740) = rxt(k,291)*y(k,56) + rxt(k,144)*y(k,79) + rxt(k,250)*y(k,83) &
                      + rxt(k,227)*y(k,87)
         mat(k,1912) = rxt(k,277)*y(k,44) + rxt(k,161)*y(k,79) + rxt(k,182)*y(k,119) &
                      + rxt(k,164)*y(k,138) + rxt(k,484)*y(k,142) + rxt(k,490) &
                      *y(k,155) + rxt(k,495)*y(k,157)
         mat(k,1441) = -(rxt(k,144)*y(k,250) + rxt(k,152)*y(k,138) + rxt(k,161) &
                      *y(k,251) + rxt(k,196)*y(k,58) + rxt(k,252)*y(k,75))
         mat(k,1736) = -rxt(k,144)*y(k,79)
         mat(k,1987) = -rxt(k,152)*y(k,79)
         mat(k,1908) = -rxt(k,161)*y(k,79)
         mat(k,2250) = -rxt(k,196)*y(k,79)
         mat(k,914) = -rxt(k,252)*y(k,79)
         mat(k,1469) = rxt(k,292)*y(k,250)
         mat(k,1500) = rxt(k,154)*y(k,240)
         mat(k,2211) = rxt(k,154)*y(k,78)
         mat(k,1736) = mat(k,1736) + rxt(k,292)*y(k,56)
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
         mat(k,168) = -(rxt(k,248)*y(k,250))
         mat(k,1719) = -rxt(k,248)*y(k,80)
         mat(k,639) = -(rxt(k,153)*y(k,138) + rxt(k,162)*y(k,251) + rxt(k,197)*y(k,58))
         mat(k,1979) = -rxt(k,153)*y(k,81)
         mat(k,1847) = -rxt(k,162)*y(k,81)
         mat(k,2239) = -rxt(k,197)*y(k,81)
         mat(k,2163) = 2.000_r8*rxt(k,168)*y(k,240)
         mat(k,1847) = mat(k,1847) + 2.000_r8*rxt(k,167)*y(k,251)
         mat(k,298) = rxt(k,497)*y(k,261)
         mat(k,2331) = rxt(k,497)*y(k,159)
         mat(k,849) = -(rxt(k,243)*y(k,138) + rxt(k,244)*y(k,251) + (rxt(k,249) &
                      + rxt(k,250)) * y(k,250))
         mat(k,1982) = -rxt(k,243)*y(k,83)
         mat(k,1869) = -rxt(k,244)*y(k,83)
         mat(k,1734) = -(rxt(k,249) + rxt(k,250)) * y(k,83)
         mat(k,1453) = rxt(k,230)*y(k,44) + rxt(k,231)*y(k,240)
         mat(k,2008) = rxt(k,230)*y(k,19)
         mat(k,2182) = rxt(k,231)*y(k,19)
         mat(k,264) = -(rxt(k,266)*y(k,251) + rxt(k,271)*y(k,250))
         mat(k,1793) = -rxt(k,266)*y(k,84)
         mat(k,1726) = -rxt(k,271)*y(k,84)
         mat(k,285) = -(rxt(k,267)*y(k,251) + rxt(k,272)*y(k,250))
         mat(k,1798) = -rxt(k,267)*y(k,85)
         mat(k,1728) = -rxt(k,272)*y(k,85)
         mat(k,340) = -(rxt(k,268)*y(k,251) + rxt(k,273)*y(k,250))
         mat(k,1807) = -rxt(k,268)*y(k,86)
         mat(k,1730) = -rxt(k,273)*y(k,86)
         mat(k,1967) = -(rxt(k,214)*y(k,138) + rxt(k,215)*y(k,251) + (rxt(k,226) &
                      + rxt(k,227)) * y(k,250) + (rxt(k,546) + rxt(k,552) + rxt(k,557) &
                      ) * y(k,94) + (rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,62) &
                      + (rxt(k,553) + rxt(k,558)) * y(k,93))
         mat(k,1998) = -rxt(k,214)*y(k,87)
         mat(k,1919) = -rxt(k,215)*y(k,87)
         mat(k,1747) = -(rxt(k,226) + rxt(k,227)) * y(k,87)
         mat(k,845) = -(rxt(k,546) + rxt(k,552) + rxt(k,557)) * y(k,87)
         mat(k,929) = -(rxt(k,551) + rxt(k,556) + rxt(k,561)) * y(k,87)
         mat(k,821) = -(rxt(k,553) + rxt(k,558)) * y(k,87)
         mat(k,333) = rxt(k,305)*y(k,58)
         mat(k,521) = rxt(k,257)*y(k,58)
         mat(k,2022) = rxt(k,194)*y(k,58)
         mat(k,619) = rxt(k,259)*y(k,58)
         mat(k,411) = 2.000_r8*rxt(k,262)*y(k,58)
         mat(k,1476) = rxt(k,195)*y(k,58)
         mat(k,461) = rxt(k,264)*y(k,58)
         mat(k,2261) = rxt(k,305)*y(k,30) + rxt(k,257)*y(k,43) + rxt(k,194)*y(k,44) &
                      + rxt(k,259)*y(k,45) + 2.000_r8*rxt(k,262)*y(k,48) + rxt(k,195) &
                      *y(k,56) + rxt(k,264)*y(k,57) + rxt(k,196)*y(k,79) + rxt(k,197) &
                      *y(k,81) + rxt(k,216)*y(k,94) + rxt(k,198)*y(k,240)
         mat(k,1945) = rxt(k,213)*y(k,251)
         mat(k,1448) = rxt(k,196)*y(k,58)
         mat(k,641) = rxt(k,197)*y(k,58)
         mat(k,845) = mat(k,845) + rxt(k,216)*y(k,58)
         mat(k,2222) = rxt(k,198)*y(k,58)
         mat(k,1919) = mat(k,1919) + rxt(k,213)*y(k,61)
         mat(k,225) = -(rxt(k,285)*y(k,251) + rxt(k,293)*y(k,250))
         mat(k,1786) = -rxt(k,285)*y(k,88)
         mat(k,1725) = -rxt(k,293)*y(k,88)
         mat(k,1010) = -(rxt(k,286)*y(k,251))
         mat(k,1883) = -rxt(k,286)*y(k,89)
         mat(k,989) = .050_r8*rxt(k,459)*y(k,139)
         mat(k,326) = .350_r8*rxt(k,296)*y(k,251)
         mat(k,600) = .370_r8*rxt(k,298)*y(k,139)
         mat(k,1063) = .120_r8*rxt(k,327)*y(k,139)
         mat(k,875) = .110_r8*rxt(k,404)*y(k,139)
         mat(k,1247) = .330_r8*rxt(k,357)*y(k,139)
         mat(k,1037) = .050_r8*rxt(k,462)*y(k,139)
         mat(k,1365) = .120_r8*rxt(k,371)*y(k,139)
         mat(k,1578) = rxt(k,289)*y(k,241)
         mat(k,2290) = .050_r8*rxt(k,459)*y(k,6) + .370_r8*rxt(k,298)*y(k,27) &
                      + .120_r8*rxt(k,327)*y(k,31) + .110_r8*rxt(k,404)*y(k,100) &
                      + .330_r8*rxt(k,357)*y(k,111) + .050_r8*rxt(k,462)*y(k,116) &
                      + .120_r8*rxt(k,371)*y(k,118)
         mat(k,2191) = rxt(k,287)*y(k,241)
         mat(k,490) = rxt(k,289)*y(k,129) + rxt(k,287)*y(k,240)
         mat(k,1883) = mat(k,1883) + .350_r8*rxt(k,296)*y(k,26)
         mat(k,1465) = rxt(k,251)*y(k,75)
         mat(k,912) = rxt(k,251)*y(k,56) + rxt(k,252)*y(k,79) + rxt(k,254)*y(k,91) &
                      + rxt(k,253)*y(k,261)
         mat(k,1439) = rxt(k,252)*y(k,75)
         mat(k,1483) = rxt(k,254)*y(k,75)
         mat(k,2333) = rxt(k,253)*y(k,75)
         mat(k,1487) = -(rxt(k,191)*y(k,251) + rxt(k,254)*y(k,75))
         mat(k,1911) = -rxt(k,191)*y(k,91)
         mat(k,916) = -rxt(k,254)*y(k,91)
         mat(k,2014) = rxt(k,275)*y(k,131)
         mat(k,1128) = rxt(k,307)*y(k,131)
         mat(k,1263) = rxt(k,333)*y(k,131)
         mat(k,925) = (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,87)
         mat(k,320) = rxt(k,481)*y(k,131)
         mat(k,1960) = (rxt(k,551)+rxt(k,556)+rxt(k,561))*y(k,62)
         mat(k,1698) = rxt(k,190)*y(k,251)
         mat(k,2071) = rxt(k,275)*y(k,44) + rxt(k,307)*y(k,47) + rxt(k,333)*y(k,51) &
                      + rxt(k,481)*y(k,69)
         mat(k,1911) = mat(k,1911) + rxt(k,190)*y(k,130)
         mat(k,445) = -(rxt(k,169)*y(k,251))
         mat(k,1822) = -rxt(k,169)*y(k,92)
         mat(k,1674) = rxt(k,188)*y(k,240)
         mat(k,2148) = rxt(k,188)*y(k,130)
         mat(k,818) = -(rxt(k,245)*y(k,138) + (rxt(k,553) + rxt(k,558)) * y(k,87))
         mat(k,1980) = -rxt(k,245)*y(k,93)
         mat(k,1956) = -(rxt(k,553) + rxt(k,558)) * y(k,93)
         mat(k,2091) = rxt(k,237)*y(k,240)
         mat(k,2179) = rxt(k,237)*y(k,21)
         mat(k,842) = -(rxt(k,216)*y(k,58) + rxt(k,217)*y(k,138) + rxt(k,218)*y(k,251) &
                      + (rxt(k,546) + rxt(k,552) + rxt(k,557)) * y(k,87))
         mat(k,2241) = -rxt(k,216)*y(k,94)
         mat(k,1981) = -rxt(k,217)*y(k,94)
         mat(k,1868) = -rxt(k,218)*y(k,94)
         mat(k,1957) = -(rxt(k,546) + rxt(k,552) + rxt(k,557)) * y(k,94)
         mat(k,1932) = rxt(k,205)*y(k,240)
         mat(k,923) = rxt(k,210)*y(k,251)
         mat(k,2181) = rxt(k,205)*y(k,61)
         mat(k,1868) = mat(k,1868) + rxt(k,210)*y(k,62)
         mat(k,1193) = -(rxt(k,350)*y(k,251))
         mat(k,1896) = -rxt(k,350)*y(k,95)
         mat(k,660) = .300_r8*rxt(k,395)*y(k,251)
         mat(k,586) = .500_r8*rxt(k,396)*y(k,251)
         mat(k,1589) = rxt(k,349)*y(k,237) + rxt(k,356)*y(k,244)
         mat(k,609) = rxt(k,349)*y(k,129)
         mat(k,1349) = rxt(k,356)*y(k,129)
         mat(k,1896) = mat(k,1896) + .300_r8*rxt(k,395)*y(k,101) + .500_r8*rxt(k,396) &
                      *y(k,102)
         mat(k,293) = -(rxt(k,381)*y(k,251))
         mat(k,1800) = -rxt(k,381)*y(k,96)
         mat(k,1180) = -(rxt(k,335)*y(k,251))
         mat(k,1895) = -rxt(k,335)*y(k,97)
         mat(k,659) = .700_r8*rxt(k,395)*y(k,251)
         mat(k,585) = .500_r8*rxt(k,396)*y(k,251)
         mat(k,632) = .500_r8*rxt(k,370)*y(k,251)
         mat(k,1588) = .050_r8*rxt(k,393)*y(k,243) + .220_r8*rxt(k,355)*y(k,244) &
                      + .250_r8*rxt(k,412)*y(k,258)
         mat(k,2056) = .050_r8*rxt(k,394)*y(k,243) + .220_r8*rxt(k,354)*y(k,244) &
                      + .250_r8*rxt(k,413)*y(k,258)
         mat(k,593) = .500_r8*rxt(k,339)*y(k,251)
         mat(k,1416) = .220_r8*rxt(k,351)*y(k,244) + .250_r8*rxt(k,409)*y(k,258)
         mat(k,1642) = .230_r8*rxt(k,352)*y(k,244) + .200_r8*rxt(k,340)*y(k,254) &
                      + .100_r8*rxt(k,410)*y(k,258)
         mat(k,1323) = .050_r8*rxt(k,393)*y(k,129) + .050_r8*rxt(k,394)*y(k,131)
         mat(k,1348) = .220_r8*rxt(k,355)*y(k,129) + .220_r8*rxt(k,354)*y(k,131) &
                      + .220_r8*rxt(k,351)*y(k,234) + .230_r8*rxt(k,352)*y(k,235)
         mat(k,1895) = mat(k,1895) + .700_r8*rxt(k,395)*y(k,101) + .500_r8*rxt(k,396) &
                      *y(k,102) + .500_r8*rxt(k,370)*y(k,115) + .500_r8*rxt(k,339) &
                      *y(k,153)
         mat(k,1170) = .200_r8*rxt(k,340)*y(k,235)
         mat(k,1212) = .250_r8*rxt(k,412)*y(k,129) + .250_r8*rxt(k,413)*y(k,131) &
                      + .250_r8*rxt(k,409)*y(k,234) + .100_r8*rxt(k,410)*y(k,235)
         mat(k,383) = -(rxt(k,382)*y(k,251))
         mat(k,1812) = -rxt(k,382)*y(k,98)
         mat(k,1544) = .870_r8*rxt(k,393)*y(k,243)
         mat(k,2035) = .950_r8*rxt(k,394)*y(k,243)
         mat(k,1407) = rxt(k,389)*y(k,243)
         mat(k,1620) = .750_r8*rxt(k,390)*y(k,243)
         mat(k,1313) = .870_r8*rxt(k,393)*y(k,129) + .950_r8*rxt(k,394)*y(k,131) &
                      + rxt(k,389)*y(k,234) + .750_r8*rxt(k,390)*y(k,235)
         mat(k,172) = -(rxt(k,383)*y(k,251))
         mat(k,1782) = -rxt(k,383)*y(k,99)
         mat(k,729) = .600_r8*rxt(k,406)*y(k,251)
         mat(k,1782) = mat(k,1782) + .600_r8*rxt(k,406)*y(k,106)
         mat(k,874) = -(rxt(k,397)*y(k,131) + rxt(k,404)*y(k,139) + rxt(k,405) &
                      *y(k,251))
         mat(k,2038) = -rxt(k,397)*y(k,100)
         mat(k,2284) = -rxt(k,404)*y(k,100)
         mat(k,1872) = -rxt(k,405)*y(k,100)
         mat(k,657) = -(rxt(k,395)*y(k,251))
         mat(k,1849) = -rxt(k,395)*y(k,101)
         mat(k,1557) = .080_r8*rxt(k,387)*y(k,242)
         mat(k,1285) = .080_r8*rxt(k,387)*y(k,129)
         mat(k,582) = -(rxt(k,396)*y(k,251))
         mat(k,1840) = -rxt(k,396)*y(k,102)
         mat(k,1555) = .080_r8*rxt(k,393)*y(k,243)
         mat(k,1314) = .080_r8*rxt(k,393)*y(k,129)
         mat(k,1233) = -(rxt(k,398)*y(k,234) + rxt(k,399)*y(k,235) + rxt(k,400) &
                      *y(k,240) + rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131))
         mat(k,1418) = -rxt(k,398)*y(k,103)
         mat(k,1644) = -rxt(k,399)*y(k,103)
         mat(k,2202) = -rxt(k,400)*y(k,103)
         mat(k,1591) = -rxt(k,401)*y(k,103)
         mat(k,2059) = -rxt(k,402)*y(k,103)
         mat(k,878) = rxt(k,397)*y(k,131)
         mat(k,2059) = mat(k,2059) + rxt(k,397)*y(k,100)
         mat(k,433) = -(rxt(k,403)*y(k,251))
         mat(k,1820) = -rxt(k,403)*y(k,104)
         mat(k,1224) = rxt(k,400)*y(k,240)
         mat(k,2147) = rxt(k,400)*y(k,103)
         mat(k,86) = -(rxt(k,521)*y(k,240) + rxt(k,522)*y(k,129))
         mat(k,2124) = -rxt(k,521)*y(k,105)
         mat(k,1528) = -rxt(k,522)*y(k,105)
         mat(k,873) = rxt(k,524)*y(k,251)
         mat(k,1764) = rxt(k,524)*y(k,100)
         mat(k,730) = -(rxt(k,406)*y(k,251))
         mat(k,1858) = -rxt(k,406)*y(k,106)
         mat(k,2171) = rxt(k,386)*y(k,242) + rxt(k,391)*y(k,243)
         mat(k,1286) = rxt(k,386)*y(k,240)
         mat(k,1316) = rxt(k,391)*y(k,240)
         mat(k,69) = -(rxt(k,527)*y(k,251))
         mat(k,1762) = -rxt(k,527)*y(k,107)
         mat(k,67) = -(rxt(k,525)*y(k,240) + rxt(k,526)*y(k,129))
         mat(k,2117) = -rxt(k,525)*y(k,108)
         mat(k,1521) = -rxt(k,526)*y(k,108)
         mat(k,68) = rxt(k,527)*y(k,251)
         mat(k,1761) = rxt(k,527)*y(k,107)
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
         mat(k,105) = -(rxt(k,530)*y(k,251))
         mat(k,1772) = -rxt(k,530)*y(k,109)
         mat(k,103) = -(rxt(k,528)*y(k,240) + rxt(k,529)*y(k,129))
         mat(k,2131) = -rxt(k,528)*y(k,110)
         mat(k,1535) = -rxt(k,529)*y(k,110)
         mat(k,104) = rxt(k,530)*y(k,251)
         mat(k,1771) = rxt(k,530)*y(k,109)
         mat(k,1249) = -(rxt(k,357)*y(k,139) + rxt(k,358)*y(k,251))
         mat(k,2303) = -rxt(k,357)*y(k,111)
         mat(k,1899) = -rxt(k,358)*y(k,111)
         mat(k,879) = .300_r8*rxt(k,404)*y(k,139)
         mat(k,1592) = .360_r8*rxt(k,387)*y(k,242)
         mat(k,2060) = .400_r8*rxt(k,388)*y(k,242)
         mat(k,2303) = mat(k,2303) + .300_r8*rxt(k,404)*y(k,100)
         mat(k,1419) = .390_r8*rxt(k,384)*y(k,242)
         mat(k,1645) = .310_r8*rxt(k,385)*y(k,242)
         mat(k,1294) = .360_r8*rxt(k,387)*y(k,129) + .400_r8*rxt(k,388)*y(k,131) &
                      + .390_r8*rxt(k,384)*y(k,234) + .310_r8*rxt(k,385)*y(k,235)
         mat(k,386) = -(rxt(k,359)*y(k,251))
         mat(k,1813) = -rxt(k,359)*y(k,112)
         mat(k,2143) = rxt(k,353)*y(k,244)
         mat(k,1344) = rxt(k,353)*y(k,240)
         mat(k,548) = -(rxt(k,368)*y(k,251))
         mat(k,1836) = -rxt(k,368)*y(k,113)
         mat(k,1553) = .800_r8*rxt(k,377)*y(k,228)
         mat(k,934) = .800_r8*rxt(k,377)*y(k,129)
         mat(k,352) = -(rxt(k,369)*y(k,251))
         mat(k,1808) = -rxt(k,369)*y(k,114)
         mat(k,2139) = .800_r8*rxt(k,366)*y(k,248)
         mat(k,716) = .800_r8*rxt(k,366)*y(k,240)
         mat(k,630) = -(rxt(k,370)*y(k,251))
         mat(k,1846) = -rxt(k,370)*y(k,115)
         mat(k,1680) = rxt(k,373)*y(k,246)
         mat(k,1388) = rxt(k,373)*y(k,130)
         mat(k,1039) = -(rxt(k,461)*y(k,131) + rxt(k,462)*y(k,139) + rxt(k,463) &
                      *y(k,251))
         mat(k,2046) = -rxt(k,461)*y(k,116)
         mat(k,2292) = -rxt(k,462)*y(k,116)
         mat(k,1885) = -rxt(k,463)*y(k,116)
         mat(k,92) = -(rxt(k,532)*y(k,240) + rxt(k,533)*y(k,129))
         mat(k,2125) = -rxt(k,532)*y(k,117)
         mat(k,1529) = -rxt(k,533)*y(k,117)
         mat(k,1029) = rxt(k,535)*y(k,251)
         mat(k,1765) = rxt(k,535)*y(k,116)
         mat(k,1372) = -(rxt(k,371)*y(k,139) + rxt(k,372)*y(k,251))
         mat(k,2309) = -rxt(k,371)*y(k,118)
         mat(k,1905) = -rxt(k,372)*y(k,118)
         mat(k,882) = .200_r8*rxt(k,404)*y(k,139)
         mat(k,1597) = .560_r8*rxt(k,387)*y(k,242)
         mat(k,2066) = .600_r8*rxt(k,388)*y(k,242)
         mat(k,2309) = mat(k,2309) + .200_r8*rxt(k,404)*y(k,100)
         mat(k,1424) = .610_r8*rxt(k,384)*y(k,242)
         mat(k,1650) = .440_r8*rxt(k,385)*y(k,242)
         mat(k,1298) = .560_r8*rxt(k,387)*y(k,129) + .600_r8*rxt(k,388)*y(k,131) &
                      + .610_r8*rxt(k,384)*y(k,234) + .440_r8*rxt(k,385)*y(k,235)
         mat(k,506) = -(rxt(k,170)*y(k,129) + (rxt(k,171) + rxt(k,172) + rxt(k,173) &
                      ) * y(k,130) + rxt(k,182)*y(k,251))
         mat(k,1549) = -rxt(k,170)*y(k,119)
         mat(k,1675) = -(rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,119)
         mat(k,1830) = -rxt(k,182)*y(k,119)
         mat(k,222) = -((rxt(k,186) + rxt(k,187)) * y(k,250))
         mat(k,1724) = -(rxt(k,186) + rxt(k,187)) * y(k,120)
         mat(k,505) = rxt(k,171)*y(k,130)
         mat(k,1672) = rxt(k,171)*y(k,119)
         mat(k,1673) = rxt(k,189)*y(k,131)
         mat(k,2034) = rxt(k,189)*y(k,130)
         mat(k,439) = -(rxt(k,407)*y(k,251))
         mat(k,1821) = -rxt(k,407)*y(k,122)
         mat(k,1225) = .200_r8*rxt(k,399)*y(k,235)
         mat(k,1621) = .200_r8*rxt(k,399)*y(k,103)
         mat(k,1085) = -(rxt(k,408)*y(k,251))
         mat(k,1887) = -rxt(k,408)*y(k,123)
         mat(k,1229) = rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131) + rxt(k,398)*y(k,234) &
                      + .800_r8*rxt(k,399)*y(k,235)
         mat(k,1580) = rxt(k,401)*y(k,103)
         mat(k,2048) = rxt(k,402)*y(k,103)
         mat(k,1413) = rxt(k,398)*y(k,103)
         mat(k,1634) = .800_r8*rxt(k,399)*y(k,103)
         mat(k,136) = -(rxt(k,498)*y(k,251))
         mat(k,1778) = -rxt(k,498)*y(k,127)
         mat(k,1604) = -(rxt(k,170)*y(k,119) + rxt(k,179)*y(k,131) + rxt(k,183) &
                      *y(k,240) + rxt(k,184)*y(k,139) + rxt(k,185)*y(k,138) + rxt(k,206) &
                      *y(k,61) + rxt(k,238)*y(k,21) + rxt(k,281)*y(k,235) + rxt(k,289) &
                      *y(k,241) + rxt(k,302)*y(k,231) + rxt(k,313)*y(k,234) + rxt(k,317) &
                      *y(k,239) + rxt(k,330)*y(k,232) + rxt(k,338)*y(k,253) + rxt(k,342) &
                      *y(k,254) + (rxt(k,348) + rxt(k,349)) * y(k,237) + (rxt(k,355) &
                      + rxt(k,356)) * y(k,244) + rxt(k,364)*y(k,246) + rxt(k,367) &
                      *y(k,248) + (rxt(k,377) + rxt(k,378)) * y(k,228) + rxt(k,387) &
                      *y(k,242) + rxt(k,393)*y(k,243) + rxt(k,401)*y(k,103) + rxt(k,412) &
                      *y(k,258) + rxt(k,416)*y(k,227) + rxt(k,419)*y(k,229) + rxt(k,424) &
                      *y(k,230) + rxt(k,426)*y(k,233) + rxt(k,430)*y(k,236) + rxt(k,433) &
                      *y(k,245) + rxt(k,436)*y(k,247) + rxt(k,439)*y(k,252) + rxt(k,446) &
                      *y(k,257) + rxt(k,452)*y(k,259) + rxt(k,455)*y(k,260) + rxt(k,466) &
                      *y(k,249) + rxt(k,471)*y(k,255) + rxt(k,476)*y(k,256))
         mat(k,508) = -rxt(k,170)*y(k,129)
         mat(k,2073) = -rxt(k,179)*y(k,129)
         mat(k,2216) = -rxt(k,183)*y(k,129)
         mat(k,2316) = -rxt(k,184)*y(k,129)
         mat(k,1992) = -rxt(k,185)*y(k,129)
         mat(k,1939) = -rxt(k,206)*y(k,129)
         mat(k,2097) = -rxt(k,238)*y(k,129)
         mat(k,1656) = -rxt(k,281)*y(k,129)
         mat(k,491) = -rxt(k,289)*y(k,129)
         mat(k,862) = -rxt(k,302)*y(k,129)
         mat(k,1429) = -rxt(k,313)*y(k,129)
         mat(k,756) = -rxt(k,317)*y(k,129)
         mat(k,834) = -rxt(k,330)*y(k,129)
         mat(k,811) = -rxt(k,338)*y(k,129)
         mat(k,1173) = -rxt(k,342)*y(k,129)
         mat(k,610) = -(rxt(k,348) + rxt(k,349)) * y(k,129)
         mat(k,1355) = -(rxt(k,355) + rxt(k,356)) * y(k,129)
         mat(k,1397) = -rxt(k,364)*y(k,129)
         mat(k,720) = -rxt(k,367)*y(k,129)
         mat(k,943) = -(rxt(k,377) + rxt(k,378)) * y(k,129)
         mat(k,1302) = -rxt(k,387)*y(k,129)
         mat(k,1334) = -rxt(k,393)*y(k,129)
         mat(k,1239) = -rxt(k,401)*y(k,129)
         mat(k,1216) = -rxt(k,412)*y(k,129)
         mat(k,556) = -rxt(k,416)*y(k,129)
         mat(k,536) = -rxt(k,419)*y(k,129)
         mat(k,479) = -rxt(k,424)*y(k,129)
         mat(k,679) = -rxt(k,426)*y(k,129)
         mat(k,802) = -rxt(k,430)*y(k,129)
         mat(k,764) = -rxt(k,433)*y(k,129)
         mat(k,907) = -rxt(k,436)*y(k,129)
         mat(k,498) = -rxt(k,439)*y(k,129)
         mat(k,778) = -rxt(k,446)*y(k,129)
         mat(k,795) = -rxt(k,452)*y(k,129)
         mat(k,544) = -rxt(k,455)*y(k,129)
         mat(k,1106) = -rxt(k,466)*y(k,129)
         mat(k,1159) = -rxt(k,471)*y(k,129)
         mat(k,956) = -rxt(k,476)*y(k,129)
         mat(k,508) = mat(k,508) + 2.000_r8*rxt(k,172)*y(k,130) + rxt(k,182)*y(k,251)
         mat(k,223) = 2.000_r8*rxt(k,186)*y(k,250)
         mat(k,1700) = 2.000_r8*rxt(k,172)*y(k,119) + rxt(k,175)*y(k,138) + rxt(k,491) &
                      *y(k,157)
         mat(k,1992) = mat(k,1992) + rxt(k,175)*y(k,130)
         mat(k,1275) = rxt(k,491)*y(k,130)
         mat(k,1741) = 2.000_r8*rxt(k,186)*y(k,120)
         mat(k,1913) = rxt(k,182)*y(k,119)
         mat(k,1702) = -((rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,119) + (rxt(k,175) &
                      + rxt(k,177)) * y(k,138) + rxt(k,176)*y(k,139) + rxt(k,188) &
                      *y(k,240) + rxt(k,189)*y(k,131) + rxt(k,190)*y(k,251) + rxt(k,208) &
                      *y(k,61) + rxt(k,239)*y(k,21) + rxt(k,324)*y(k,234) + rxt(k,373) &
                      *y(k,246) + rxt(k,431)*y(k,236) + rxt(k,434)*y(k,245) + rxt(k,437) &
                      *y(k,247) + rxt(k,441)*y(k,146) + rxt(k,444)*y(k,227) + rxt(k,491) &
                      *y(k,157))
         mat(k,509) = -(rxt(k,171) + rxt(k,172) + rxt(k,173)) * y(k,130)
         mat(k,1994) = -(rxt(k,175) + rxt(k,177)) * y(k,130)
         mat(k,2318) = -rxt(k,176)*y(k,130)
         mat(k,2218) = -rxt(k,188)*y(k,130)
         mat(k,2075) = -rxt(k,189)*y(k,130)
         mat(k,1915) = -rxt(k,190)*y(k,130)
         mat(k,1941) = -rxt(k,208)*y(k,130)
         mat(k,2099) = -rxt(k,239)*y(k,130)
         mat(k,1431) = -rxt(k,324)*y(k,130)
         mat(k,1399) = -rxt(k,373)*y(k,130)
         mat(k,804) = -rxt(k,431)*y(k,130)
         mat(k,765) = -rxt(k,434)*y(k,130)
         mat(k,909) = -rxt(k,437)*y(k,130)
         mat(k,515) = -rxt(k,441)*y(k,130)
         mat(k,557) = -rxt(k,444)*y(k,130)
         mat(k,1276) = -rxt(k,491)*y(k,130)
         mat(k,691) = rxt(k,375)*y(k,251)
         mat(k,403) = rxt(k,346)*y(k,131)
         mat(k,2099) = mat(k,2099) + rxt(k,238)*y(k,129)
         mat(k,1941) = mat(k,1941) + rxt(k,206)*y(k,129)
         mat(k,446) = rxt(k,169)*y(k,251)
         mat(k,662) = .700_r8*rxt(k,395)*y(k,251)
         mat(k,1241) = rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131)
         mat(k,1606) = rxt(k,238)*y(k,21) + rxt(k,206)*y(k,61) + rxt(k,401)*y(k,103) &
                      + 2.000_r8*rxt(k,179)*y(k,131) + rxt(k,185)*y(k,138) &
                      + rxt(k,184)*y(k,139) + rxt(k,416)*y(k,227) + rxt(k,377) &
                      *y(k,228) + rxt(k,419)*y(k,229) + rxt(k,424)*y(k,230) &
                      + rxt(k,302)*y(k,231) + rxt(k,330)*y(k,232) + rxt(k,426) &
                      *y(k,233) + rxt(k,313)*y(k,234) + rxt(k,281)*y(k,235) &
                      + rxt(k,430)*y(k,236) + rxt(k,348)*y(k,237) + rxt(k,317) &
                      *y(k,239) + rxt(k,183)*y(k,240) + rxt(k,289)*y(k,241) &
                      + .920_r8*rxt(k,387)*y(k,242) + .920_r8*rxt(k,393)*y(k,243) &
                      + rxt(k,355)*y(k,244) + rxt(k,433)*y(k,245) + rxt(k,364) &
                      *y(k,246) + rxt(k,436)*y(k,247) + rxt(k,367)*y(k,248) &
                      + 1.600_r8*rxt(k,466)*y(k,249) + rxt(k,439)*y(k,252) &
                      + rxt(k,338)*y(k,253) + rxt(k,342)*y(k,254) + .900_r8*rxt(k,471) &
                      *y(k,255) + .800_r8*rxt(k,476)*y(k,256) + rxt(k,446)*y(k,257) &
                      + rxt(k,412)*y(k,258) + rxt(k,452)*y(k,259) + rxt(k,455) &
                      *y(k,260)
         mat(k,2075) = mat(k,2075) + rxt(k,346)*y(k,18) + rxt(k,402)*y(k,103) &
                      + 2.000_r8*rxt(k,179)*y(k,129) + rxt(k,180)*y(k,138) &
                      + rxt(k,178)*y(k,240) + rxt(k,388)*y(k,242) + rxt(k,394) &
                      *y(k,243) + rxt(k,354)*y(k,244) + rxt(k,365)*y(k,246) &
                      + 2.000_r8*rxt(k,467)*y(k,249) + rxt(k,181)*y(k,251) &
                      + rxt(k,413)*y(k,258)
         mat(k,893) = rxt(k,336)*y(k,251)
         mat(k,1994) = mat(k,1994) + rxt(k,185)*y(k,129) + rxt(k,180)*y(k,131)
         mat(k,2318) = mat(k,2318) + rxt(k,184)*y(k,129)
         mat(k,673) = rxt(k,473)*y(k,251)
         mat(k,557) = mat(k,557) + rxt(k,416)*y(k,129)
         mat(k,945) = rxt(k,377)*y(k,129)
         mat(k,537) = rxt(k,419)*y(k,129)
         mat(k,480) = rxt(k,424)*y(k,129)
         mat(k,864) = rxt(k,302)*y(k,129)
         mat(k,836) = rxt(k,330)*y(k,129)
         mat(k,680) = rxt(k,426)*y(k,129)
         mat(k,1431) = mat(k,1431) + rxt(k,313)*y(k,129)
         mat(k,1658) = rxt(k,281)*y(k,129) + .500_r8*rxt(k,464)*y(k,249)
         mat(k,804) = mat(k,804) + rxt(k,430)*y(k,129)
         mat(k,611) = rxt(k,348)*y(k,129)
         mat(k,757) = rxt(k,317)*y(k,129)
         mat(k,2218) = mat(k,2218) + rxt(k,183)*y(k,129) + rxt(k,178)*y(k,131)
         mat(k,492) = rxt(k,289)*y(k,129)
         mat(k,1304) = .920_r8*rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131)
         mat(k,1336) = .920_r8*rxt(k,393)*y(k,129) + rxt(k,394)*y(k,131)
         mat(k,1357) = rxt(k,355)*y(k,129) + rxt(k,354)*y(k,131)
         mat(k,765) = mat(k,765) + rxt(k,433)*y(k,129)
         mat(k,1399) = mat(k,1399) + rxt(k,364)*y(k,129) + rxt(k,365)*y(k,131)
         mat(k,909) = mat(k,909) + rxt(k,436)*y(k,129)
         mat(k,721) = rxt(k,367)*y(k,129)
         mat(k,1108) = 1.600_r8*rxt(k,466)*y(k,129) + 2.000_r8*rxt(k,467)*y(k,131) &
                      + .500_r8*rxt(k,464)*y(k,235)
         mat(k,1915) = mat(k,1915) + rxt(k,375)*y(k,1) + rxt(k,169)*y(k,92) &
                      + .700_r8*rxt(k,395)*y(k,101) + rxt(k,181)*y(k,131) + rxt(k,336) &
                      *y(k,132) + rxt(k,473)*y(k,212)
         mat(k,499) = rxt(k,439)*y(k,129)
         mat(k,812) = rxt(k,338)*y(k,129)
         mat(k,1175) = rxt(k,342)*y(k,129)
         mat(k,1161) = .900_r8*rxt(k,471)*y(k,129)
         mat(k,958) = .800_r8*rxt(k,476)*y(k,129)
         mat(k,779) = rxt(k,446)*y(k,129)
         mat(k,1218) = rxt(k,412)*y(k,129) + rxt(k,413)*y(k,131)
         mat(k,796) = rxt(k,452)*y(k,129)
         mat(k,545) = rxt(k,455)*y(k,129)
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
         mat(k,2082) = -(rxt(k,178)*y(k,240) + rxt(k,179)*y(k,129) + rxt(k,180) &
                      *y(k,138) + rxt(k,181)*y(k,251) + rxt(k,189)*y(k,130) + rxt(k,275) &
                      *y(k,44) + rxt(k,307)*y(k,47) + rxt(k,326)*y(k,31) + rxt(k,333) &
                      *y(k,51) + rxt(k,346)*y(k,18) + rxt(k,354)*y(k,244) + rxt(k,365) &
                      *y(k,246) + rxt(k,388)*y(k,242) + rxt(k,394)*y(k,243) + rxt(k,397) &
                      *y(k,100) + rxt(k,402)*y(k,103) + rxt(k,413)*y(k,258) + rxt(k,458) &
                      *y(k,6) + rxt(k,461)*y(k,116) + rxt(k,467)*y(k,249) + rxt(k,478) &
                      *y(k,214) + rxt(k,481)*y(k,69))
         mat(k,2225) = -rxt(k,178)*y(k,131)
         mat(k,1613) = -rxt(k,179)*y(k,131)
         mat(k,2001) = -rxt(k,180)*y(k,131)
         mat(k,1922) = -rxt(k,181)*y(k,131)
         mat(k,1709) = -rxt(k,189)*y(k,131)
         mat(k,2025) = -rxt(k,275)*y(k,131)
         mat(k,1131) = -rxt(k,307)*y(k,131)
         mat(k,1077) = -rxt(k,326)*y(k,131)
         mat(k,1266) = -rxt(k,333)*y(k,131)
         mat(k,406) = -rxt(k,346)*y(k,131)
         mat(k,1361) = -rxt(k,354)*y(k,131)
         mat(k,1403) = -rxt(k,365)*y(k,131)
         mat(k,1308) = -rxt(k,388)*y(k,131)
         mat(k,1340) = -rxt(k,394)*y(k,131)
         mat(k,887) = -rxt(k,397)*y(k,131)
         mat(k,1245) = -rxt(k,402)*y(k,131)
         mat(k,1222) = -rxt(k,413)*y(k,131)
         mat(k,1005) = -rxt(k,458)*y(k,131)
         mat(k,1055) = -rxt(k,461)*y(k,131)
         mat(k,1111) = -rxt(k,467)*y(k,131)
         mat(k,1022) = -rxt(k,478)*y(k,131)
         mat(k,322) = -rxt(k,481)*y(k,131)
         mat(k,580) = rxt(k,240)*y(k,138)
         mat(k,2264) = rxt(k,207)*y(k,62)
         mat(k,931) = rxt(k,207)*y(k,58) + rxt(k,209)*y(k,138) + rxt(k,210)*y(k,251)
         mat(k,920) = rxt(k,254)*y(k,91)
         mat(k,1496) = rxt(k,254)*y(k,75) + rxt(k,191)*y(k,251)
         mat(k,637) = .500_r8*rxt(k,370)*y(k,251)
         mat(k,1709) = mat(k,1709) + rxt(k,177)*y(k,138) + rxt(k,176)*y(k,139)
         mat(k,2001) = mat(k,2001) + rxt(k,240)*y(k,22) + rxt(k,209)*y(k,62) &
                      + rxt(k,177)*y(k,130)
         mat(k,2325) = rxt(k,176)*y(k,130)
         mat(k,573) = rxt(k,322)*y(k,251)
         mat(k,1922) = mat(k,1922) + rxt(k,210)*y(k,62) + rxt(k,191)*y(k,91) &
                      + .500_r8*rxt(k,370)*y(k,115) + rxt(k,322)*y(k,144)
         mat(k,890) = -(rxt(k,336)*y(k,251))
         mat(k,1873) = -rxt(k,336)*y(k,132)
         mat(k,1062) = rxt(k,326)*y(k,131)
         mat(k,583) = .500_r8*rxt(k,396)*y(k,251)
         mat(k,435) = rxt(k,403)*y(k,251)
         mat(k,440) = rxt(k,407)*y(k,251)
         mat(k,1082) = rxt(k,408)*y(k,251)
         mat(k,2039) = rxt(k,326)*y(k,31)
         mat(k,1873) = mat(k,1873) + .500_r8*rxt(k,396)*y(k,102) + rxt(k,403)*y(k,104) &
                      + rxt(k,407)*y(k,122) + rxt(k,408)*y(k,123)
         mat(k,463) = -(rxt(k,468)*y(k,251))
         mat(k,1825) = -rxt(k,468)*y(k,133)
         mat(k,2150) = rxt(k,465)*y(k,249)
         mat(k,1097) = rxt(k,465)*y(k,240)
         mat(k,1999) = -(rxt(k,149)*y(k,139) + 4._r8*rxt(k,150)*y(k,138) + rxt(k,152) &
                      *y(k,79) + rxt(k,153)*y(k,81) + rxt(k,158)*y(k,240) + rxt(k,164) &
                      *y(k,251) + (rxt(k,175) + rxt(k,177)) * y(k,130) + rxt(k,180) &
                      *y(k,131) + rxt(k,185)*y(k,129) + rxt(k,209)*y(k,62) + rxt(k,211) &
                      *y(k,61) + rxt(k,214)*y(k,87) + rxt(k,217)*y(k,94) + rxt(k,240) &
                      *y(k,22) + rxt(k,241)*y(k,21) + rxt(k,243)*y(k,83) + rxt(k,245) &
                      *y(k,93) + rxt(k,276)*y(k,44) + rxt(k,483)*y(k,142))
         mat(k,2323) = -rxt(k,149)*y(k,138)
         mat(k,1449) = -rxt(k,152)*y(k,138)
         mat(k,642) = -rxt(k,153)*y(k,138)
         mat(k,2223) = -rxt(k,158)*y(k,138)
         mat(k,1920) = -rxt(k,164)*y(k,138)
         mat(k,1707) = -(rxt(k,175) + rxt(k,177)) * y(k,138)
         mat(k,2080) = -rxt(k,180)*y(k,138)
         mat(k,1611) = -rxt(k,185)*y(k,138)
         mat(k,930) = -rxt(k,209)*y(k,138)
         mat(k,1946) = -rxt(k,211)*y(k,138)
         mat(k,1968) = -rxt(k,214)*y(k,138)
         mat(k,846) = -rxt(k,217)*y(k,138)
         mat(k,579) = -rxt(k,240)*y(k,138)
         mat(k,2104) = -rxt(k,241)*y(k,138)
         mat(k,854) = -rxt(k,243)*y(k,138)
         mat(k,822) = -rxt(k,245)*y(k,138)
         mat(k,2023) = -rxt(k,276)*y(k,138)
         mat(k,398) = -rxt(k,483)*y(k,138)
         mat(k,1509) = rxt(k,156)*y(k,240)
         mat(k,512) = rxt(k,170)*y(k,129) + rxt(k,171)*y(k,130)
         mat(k,1611) = mat(k,1611) + rxt(k,170)*y(k,119)
         mat(k,1707) = mat(k,1707) + rxt(k,171)*y(k,119)
         mat(k,2223) = mat(k,2223) + rxt(k,156)*y(k,78)
         mat(k,1920) = mat(k,1920) + 2.000_r8*rxt(k,166)*y(k,251)
         mat(k,2329) = -(rxt(k,148)*y(k,250) + rxt(k,149)*y(k,138) + rxt(k,159) &
                      *y(k,240) + rxt(k,160)*y(k,78) + rxt(k,165)*y(k,251) + rxt(k,176) &
                      *y(k,130) + rxt(k,184)*y(k,129) + rxt(k,200)*y(k,58) + rxt(k,232) &
                      *y(k,19) + rxt(k,298)*y(k,27) + rxt(k,327)*y(k,31) + rxt(k,357) &
                      *y(k,111) + rxt(k,371)*y(k,118) + rxt(k,404)*y(k,100) + rxt(k,442) &
                      *y(k,146) + rxt(k,459)*y(k,6) + rxt(k,462)*y(k,116) + rxt(k,487) &
                      *y(k,155) + rxt(k,493)*y(k,157))
         mat(k,1754) = -rxt(k,148)*y(k,139)
         mat(k,2005) = -rxt(k,149)*y(k,139)
         mat(k,2229) = -rxt(k,159)*y(k,139)
         mat(k,1514) = -rxt(k,160)*y(k,139)
         mat(k,1926) = -rxt(k,165)*y(k,139)
         mat(k,1713) = -rxt(k,176)*y(k,139)
         mat(k,1617) = -rxt(k,184)*y(k,139)
         mat(k,2268) = -rxt(k,200)*y(k,139)
         mat(k,1463) = -rxt(k,232)*y(k,139)
         mat(k,605) = -rxt(k,298)*y(k,139)
         mat(k,1079) = -rxt(k,327)*y(k,139)
         mat(k,1258) = -rxt(k,357)*y(k,139)
         mat(k,1385) = -rxt(k,371)*y(k,139)
         mat(k,889) = -rxt(k,404)*y(k,139)
         mat(k,516) = -rxt(k,442)*y(k,139)
         mat(k,1007) = -rxt(k,459)*y(k,139)
         mat(k,1057) = -rxt(k,462)*y(k,139)
         mat(k,565) = -rxt(k,487)*y(k,139)
         mat(k,1283) = -rxt(k,493)*y(k,139)
         mat(k,1437) = .150_r8*rxt(k,312)*y(k,240)
         mat(k,2229) = mat(k,2229) + .150_r8*rxt(k,312)*y(k,234) + .150_r8*rxt(k,362) &
                      *y(k,246)
         mat(k,1405) = .150_r8*rxt(k,362)*y(k,240)
         mat(k,362) = -(rxt(k,494)*y(k,157))
         mat(k,1269) = -rxt(k,494)*y(k,141)
         mat(k,2089) = rxt(k,234)*y(k,61)
         mat(k,1931) = rxt(k,234)*y(k,21) + 2.000_r8*rxt(k,204)*y(k,61)
         mat(k,391) = -(rxt(k,483)*y(k,138) + rxt(k,484)*y(k,251))
         mat(k,1976) = -rxt(k,483)*y(k,142)
         mat(k,1814) = -rxt(k,484)*y(k,142)
         mat(k,1186) = rxt(k,350)*y(k,251)
         mat(k,1539) = .100_r8*rxt(k,471)*y(k,255)
         mat(k,1796) = rxt(k,350)*y(k,95)
         mat(k,1147) = .100_r8*rxt(k,471)*y(k,129)
         mat(k,566) = -(rxt(k,322)*y(k,251))
         mat(k,1839) = -rxt(k,322)*y(k,144)
         mat(k,1678) = rxt(k,324)*y(k,234)
         mat(k,1408) = rxt(k,324)*y(k,130)
         mat(k,1671) = rxt(k,444)*y(k,227)
         mat(k,553) = rxt(k,444)*y(k,130)
         mat(k,513) = -(rxt(k,441)*y(k,130) + rxt(k,442)*y(k,139))
         mat(k,1676) = -rxt(k,441)*y(k,146)
         mat(k,2277) = -rxt(k,442)*y(k,146)
         mat(k,239) = .070_r8*rxt(k,428)*y(k,251)
         mat(k,1550) = rxt(k,426)*y(k,233)
         mat(k,219) = .060_r8*rxt(k,440)*y(k,251)
         mat(k,260) = .070_r8*rxt(k,456)*y(k,251)
         mat(k,677) = rxt(k,426)*y(k,129)
         mat(k,1831) = .070_r8*rxt(k,428)*y(k,68) + .060_r8*rxt(k,440)*y(k,147) &
                      + .070_r8*rxt(k,456)*y(k,223)
         mat(k,217) = -(rxt(k,440)*y(k,251))
         mat(k,1785) = -rxt(k,440)*y(k,147)
         mat(k,209) = .530_r8*rxt(k,417)*y(k,251)
         mat(k,1785) = mat(k,1785) + .530_r8*rxt(k,417)*y(k,8)
         mat(k,367) = -(rxt(k,443)*y(k,251))
         mat(k,1809) = -rxt(k,443)*y(k,148)
         mat(k,2140) = rxt(k,438)*y(k,252)
         mat(k,495) = rxt(k,438)*y(k,240)
         mat(k,590) = -(rxt(k,339)*y(k,251))
         mat(k,1841) = -rxt(k,339)*y(k,153)
         mat(k,2161) = rxt(k,337)*y(k,253)
         mat(k,807) = rxt(k,337)*y(k,240)
         mat(k,415) = -(rxt(k,343)*y(k,251))
         mat(k,1817) = -rxt(k,343)*y(k,154)
         mat(k,2144) = .850_r8*rxt(k,341)*y(k,254)
         mat(k,1167) = .850_r8*rxt(k,341)*y(k,240)
         mat(k,560) = -(rxt(k,487)*y(k,139) + rxt(k,490)*y(k,251))
         mat(k,2278) = -rxt(k,487)*y(k,155)
         mat(k,1838) = -rxt(k,490)*y(k,155)
         mat(k,1272) = -(rxt(k,488)*y(k,21) + rxt(k,489)*y(k,61) + rxt(k,491)*y(k,130) &
                      + rxt(k,493)*y(k,139) + rxt(k,494)*y(k,141) + rxt(k,495) &
                      *y(k,251))
         mat(k,2093) = -rxt(k,488)*y(k,157)
         mat(k,1935) = -rxt(k,489)*y(k,157)
         mat(k,1693) = -rxt(k,491)*y(k,157)
         mat(k,2305) = -rxt(k,493)*y(k,157)
         mat(k,364) = -rxt(k,494)*y(k,157)
         mat(k,1901) = -rxt(k,495)*y(k,157)
         mat(k,1986) = rxt(k,483)*y(k,142)
         mat(k,2305) = mat(k,2305) + rxt(k,487)*y(k,155)
         mat(k,395) = rxt(k,483)*y(k,138)
         mat(k,561) = rxt(k,487)*y(k,139) + rxt(k,490)*y(k,251)
         mat(k,1901) = mat(k,1901) + rxt(k,490)*y(k,155)
         mat(k,897) = -(rxt(k,486)*y(k,251))
         mat(k,1874) = -rxt(k,486)*y(k,158)
         mat(k,2092) = rxt(k,488)*y(k,157)
         mat(k,1933) = rxt(k,489)*y(k,157)
         mat(k,319) = rxt(k,481)*y(k,131) + (rxt(k,482)+.500_r8*rxt(k,496))*y(k,251)
         mat(k,1686) = rxt(k,491)*y(k,157)
         mat(k,2040) = rxt(k,481)*y(k,69)
         mat(k,2285) = rxt(k,493)*y(k,157)
         mat(k,363) = rxt(k,494)*y(k,157)
         mat(k,393) = rxt(k,484)*y(k,251)
         mat(k,1271) = rxt(k,488)*y(k,21) + rxt(k,489)*y(k,61) + rxt(k,491)*y(k,130) &
                      + rxt(k,493)*y(k,139) + rxt(k,494)*y(k,141) + rxt(k,495) &
                      *y(k,251)
         mat(k,1874) = mat(k,1874) + (rxt(k,482)+.500_r8*rxt(k,496))*y(k,69) &
                      + rxt(k,484)*y(k,142) + rxt(k,495)*y(k,157)
         mat(k,299) = -(rxt(k,497)*y(k,261))
         mat(k,2332) = -rxt(k,497)*y(k,159)
         mat(k,896) = rxt(k,486)*y(k,251)
         mat(k,1801) = rxt(k,486)*y(k,158)
         mat(k,62) = .1056005_r8*rxt(k,526)*y(k,129) + .2381005_r8*rxt(k,525)*y(k,240)
         mat(k,1516) = .1056005_r8*rxt(k,526)*y(k,108)
         mat(k,112) = .5931005_r8*rxt(k,536)*y(k,251)
         mat(k,2112) = .2381005_r8*rxt(k,525)*y(k,108)
         mat(k,1756) = .5931005_r8*rxt(k,536)*y(k,208)
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
         mat(k,63) = .1026005_r8*rxt(k,526)*y(k,129) + .1308005_r8*rxt(k,525)*y(k,240)
         mat(k,1517) = .1026005_r8*rxt(k,526)*y(k,108)
         mat(k,113) = .1534005_r8*rxt(k,536)*y(k,251)
         mat(k,2113) = .1308005_r8*rxt(k,525)*y(k,108)
         mat(k,1757) = .1534005_r8*rxt(k,536)*y(k,208)
         mat(k,64) = .0521005_r8*rxt(k,526)*y(k,129) + .0348005_r8*rxt(k,525)*y(k,240)
         mat(k,1518) = .0521005_r8*rxt(k,526)*y(k,108)
         mat(k,114) = .0459005_r8*rxt(k,536)*y(k,251)
         mat(k,2114) = .0348005_r8*rxt(k,525)*y(k,108)
         mat(k,1758) = .0459005_r8*rxt(k,536)*y(k,208)
         mat(k,65) = .0143005_r8*rxt(k,526)*y(k,129) + .0076005_r8*rxt(k,525)*y(k,240)
         mat(k,1519) = .0143005_r8*rxt(k,526)*y(k,108)
         mat(k,115) = .0085005_r8*rxt(k,536)*y(k,251)
         mat(k,2115) = .0076005_r8*rxt(k,525)*y(k,108)
         mat(k,1759) = .0085005_r8*rxt(k,536)*y(k,208)
         mat(k,66) = .0166005_r8*rxt(k,526)*y(k,129) + .0113005_r8*rxt(k,525)*y(k,240)
         mat(k,1520) = .0166005_r8*rxt(k,526)*y(k,108)
         mat(k,116) = .0128005_r8*rxt(k,536)*y(k,251)
         mat(k,2116) = .0113005_r8*rxt(k,525)*y(k,108)
         mat(k,1760) = .0128005_r8*rxt(k,536)*y(k,208)
         mat(k,975) = .2202005_r8*rxt(k,515)*y(k,139)
         mat(k,75) = .1279005_r8*rxt(k,514)*y(k,129) + .2202005_r8*rxt(k,513)*y(k,240)
         mat(k,81) = .0003005_r8*rxt(k,522)*y(k,129) + .0031005_r8*rxt(k,521)*y(k,240)
         mat(k,1024) = .0508005_r8*rxt(k,534)*y(k,139)
         mat(k,87) = .0245005_r8*rxt(k,533)*y(k,129) + .0508005_r8*rxt(k,532)*y(k,240)
         mat(k,1522) = .1279005_r8*rxt(k,514)*y(k,7) + .0003005_r8*rxt(k,522)*y(k,105) &
                      + .0245005_r8*rxt(k,533)*y(k,117)
         mat(k,2270) = .2202005_r8*rxt(k,515)*y(k,6) + .0508005_r8*rxt(k,534)*y(k,116)
         mat(k,2118) = .2202005_r8*rxt(k,513)*y(k,7) + .0031005_r8*rxt(k,521)*y(k,105) &
                      + .0508005_r8*rxt(k,532)*y(k,117)
         mat(k,976) = .2067005_r8*rxt(k,515)*y(k,139)
         mat(k,76) = .1792005_r8*rxt(k,514)*y(k,129) + .2067005_r8*rxt(k,513)*y(k,240)
         mat(k,82) = .0003005_r8*rxt(k,522)*y(k,129) + .0035005_r8*rxt(k,521)*y(k,240)
         mat(k,1025) = .1149005_r8*rxt(k,534)*y(k,139)
         mat(k,88) = .0082005_r8*rxt(k,533)*y(k,129) + .1149005_r8*rxt(k,532)*y(k,240)
         mat(k,1523) = .1792005_r8*rxt(k,514)*y(k,7) + .0003005_r8*rxt(k,522)*y(k,105) &
                      + .0082005_r8*rxt(k,533)*y(k,117)
         mat(k,2271) = .2067005_r8*rxt(k,515)*y(k,6) + .1149005_r8*rxt(k,534)*y(k,116)
         mat(k,2119) = .2067005_r8*rxt(k,513)*y(k,7) + .0035005_r8*rxt(k,521)*y(k,105) &
                      + .1149005_r8*rxt(k,532)*y(k,117)
         mat(k,977) = .0653005_r8*rxt(k,515)*y(k,139)
         mat(k,77) = .0676005_r8*rxt(k,514)*y(k,129) + .0653005_r8*rxt(k,513)*y(k,240)
         mat(k,83) = .0073005_r8*rxt(k,522)*y(k,129) + .0003005_r8*rxt(k,521)*y(k,240)
         mat(k,1026) = .0348005_r8*rxt(k,534)*y(k,139)
         mat(k,89) = .0772005_r8*rxt(k,533)*y(k,129) + .0348005_r8*rxt(k,532)*y(k,240)
         mat(k,1524) = .0676005_r8*rxt(k,514)*y(k,7) + .0073005_r8*rxt(k,522)*y(k,105) &
                      + .0772005_r8*rxt(k,533)*y(k,117)
         mat(k,2272) = .0653005_r8*rxt(k,515)*y(k,6) + .0348005_r8*rxt(k,534)*y(k,116)
         mat(k,2120) = .0653005_r8*rxt(k,513)*y(k,7) + .0003005_r8*rxt(k,521)*y(k,105) &
                      + .0348005_r8*rxt(k,532)*y(k,117)
         mat(k,978) = .1749305_r8*rxt(k,512)*y(k,131) + .1284005_r8*rxt(k,515) &
                      *y(k,139)
         mat(k,78) = .079_r8*rxt(k,514)*y(k,129) + .1284005_r8*rxt(k,513)*y(k,240)
         mat(k,871) = .0590245_r8*rxt(k,520)*y(k,131) + .0033005_r8*rxt(k,523) &
                      *y(k,139)
         mat(k,84) = .0057005_r8*rxt(k,522)*y(k,129) + .0271005_r8*rxt(k,521)*y(k,240)
         mat(k,1027) = .1749305_r8*rxt(k,531)*y(k,131) + .0554005_r8*rxt(k,534) &
                      *y(k,139)
         mat(k,90) = .0332005_r8*rxt(k,533)*y(k,129) + .0554005_r8*rxt(k,532)*y(k,240)
         mat(k,1525) = .079_r8*rxt(k,514)*y(k,7) + .0057005_r8*rxt(k,522)*y(k,105) &
                      + .0332005_r8*rxt(k,533)*y(k,117)
         mat(k,2031) = .1749305_r8*rxt(k,512)*y(k,6) + .0590245_r8*rxt(k,520)*y(k,100) &
                      + .1749305_r8*rxt(k,531)*y(k,116)
         mat(k,2273) = .1284005_r8*rxt(k,515)*y(k,6) + .0033005_r8*rxt(k,523)*y(k,100) &
                      + .0554005_r8*rxt(k,534)*y(k,116)
         mat(k,2121) = .1284005_r8*rxt(k,513)*y(k,7) + .0271005_r8*rxt(k,521)*y(k,105) &
                      + .0554005_r8*rxt(k,532)*y(k,117)
         mat(k,979) = .5901905_r8*rxt(k,512)*y(k,131) + .114_r8*rxt(k,515)*y(k,139)
         mat(k,79) = .1254005_r8*rxt(k,514)*y(k,129) + .114_r8*rxt(k,513)*y(k,240)
         mat(k,872) = .0250245_r8*rxt(k,520)*y(k,131)
         mat(k,85) = .0623005_r8*rxt(k,522)*y(k,129) + .0474005_r8*rxt(k,521)*y(k,240)
         mat(k,1028) = .5901905_r8*rxt(k,531)*y(k,131) + .1278005_r8*rxt(k,534) &
                      *y(k,139)
         mat(k,91) = .130_r8*rxt(k,533)*y(k,129) + .1278005_r8*rxt(k,532)*y(k,240)
         mat(k,1526) = .1254005_r8*rxt(k,514)*y(k,7) + .0623005_r8*rxt(k,522)*y(k,105) &
                      + .130_r8*rxt(k,533)*y(k,117)
         mat(k,2032) = .5901905_r8*rxt(k,512)*y(k,6) + .0250245_r8*rxt(k,520)*y(k,100) &
                      + .5901905_r8*rxt(k,531)*y(k,116)
         mat(k,2274) = .114_r8*rxt(k,515)*y(k,6) + .1278005_r8*rxt(k,534)*y(k,116)
         mat(k,2122) = .114_r8*rxt(k,513)*y(k,7) + .0474005_r8*rxt(k,521)*y(k,105) &
                      + .1278005_r8*rxt(k,532)*y(k,117)
         mat(k,106) = .0097005_r8*rxt(k,519)*y(k,129) + .0023005_r8*rxt(k,518) &
                      *y(k,240)
         mat(k,98) = .1056005_r8*rxt(k,529)*y(k,129) + .2381005_r8*rxt(k,528)*y(k,240)
         mat(k,1530) = .0097005_r8*rxt(k,519)*y(k,9) + .1056005_r8*rxt(k,529)*y(k,110) &
                      + .0154005_r8*rxt(k,540)*y(k,218) + .0063005_r8*rxt(k,544) &
                      *y(k,222)
         mat(k,118) = .5931005_r8*rxt(k,537)*y(k,251)
         mat(k,124) = .0154005_r8*rxt(k,540)*y(k,129) + .1364005_r8*rxt(k,539) &
                      *y(k,240)
         mat(k,130) = .0063005_r8*rxt(k,544)*y(k,129) + .1677005_r8*rxt(k,543) &
                      *y(k,240)
         mat(k,2126) = .0023005_r8*rxt(k,518)*y(k,9) + .2381005_r8*rxt(k,528)*y(k,110) &
                      + .1364005_r8*rxt(k,539)*y(k,218) + .1677005_r8*rxt(k,543) &
                      *y(k,222)
         mat(k,1766) = .5931005_r8*rxt(k,537)*y(k,209)
         mat(k,107) = .0034005_r8*rxt(k,519)*y(k,129) + .0008005_r8*rxt(k,518) &
                      *y(k,240)
         mat(k,99) = .1026005_r8*rxt(k,529)*y(k,129) + .1308005_r8*rxt(k,528)*y(k,240)
         mat(k,1531) = .0034005_r8*rxt(k,519)*y(k,9) + .1026005_r8*rxt(k,529)*y(k,110) &
                      + .0452005_r8*rxt(k,540)*y(k,218) + .0237005_r8*rxt(k,544) &
                      *y(k,222)
         mat(k,119) = .1534005_r8*rxt(k,537)*y(k,251)
         mat(k,125) = .0452005_r8*rxt(k,540)*y(k,129) + .0101005_r8*rxt(k,539) &
                      *y(k,240)
         mat(k,131) = .0237005_r8*rxt(k,544)*y(k,129) + .0174005_r8*rxt(k,543) &
                      *y(k,240)
         mat(k,2127) = .0008005_r8*rxt(k,518)*y(k,9) + .1308005_r8*rxt(k,528)*y(k,110) &
                      + .0101005_r8*rxt(k,539)*y(k,218) + .0174005_r8*rxt(k,543) &
                      *y(k,222)
         mat(k,1767) = .1534005_r8*rxt(k,537)*y(k,209)
         mat(k,108) = .1579005_r8*rxt(k,519)*y(k,129) + .0843005_r8*rxt(k,518) &
                      *y(k,240)
         mat(k,100) = .0521005_r8*rxt(k,529)*y(k,129) + .0348005_r8*rxt(k,528) &
                      *y(k,240)
         mat(k,1532) = .1579005_r8*rxt(k,519)*y(k,9) + .0521005_r8*rxt(k,529)*y(k,110) &
                      + .0966005_r8*rxt(k,540)*y(k,218) + .0025005_r8*rxt(k,544) &
                      *y(k,222)
         mat(k,120) = .0459005_r8*rxt(k,537)*y(k,251)
         mat(k,126) = .0966005_r8*rxt(k,540)*y(k,129) + .0763005_r8*rxt(k,539) &
                      *y(k,240)
         mat(k,132) = .0025005_r8*rxt(k,544)*y(k,129) + .086_r8*rxt(k,543)*y(k,240)
         mat(k,2128) = .0843005_r8*rxt(k,518)*y(k,9) + .0348005_r8*rxt(k,528)*y(k,110) &
                      + .0763005_r8*rxt(k,539)*y(k,218) + .086_r8*rxt(k,543)*y(k,222)
         mat(k,1768) = .0459005_r8*rxt(k,537)*y(k,209)
         mat(k,109) = .0059005_r8*rxt(k,519)*y(k,129) + .0443005_r8*rxt(k,518) &
                      *y(k,240)
         mat(k,101) = .0143005_r8*rxt(k,529)*y(k,129) + .0076005_r8*rxt(k,528) &
                      *y(k,240)
         mat(k,1533) = .0059005_r8*rxt(k,519)*y(k,9) + .0143005_r8*rxt(k,529)*y(k,110) &
                      + .0073005_r8*rxt(k,540)*y(k,218) + .011_r8*rxt(k,544)*y(k,222)
         mat(k,121) = .0085005_r8*rxt(k,537)*y(k,251)
         mat(k,127) = .0073005_r8*rxt(k,540)*y(k,129) + .2157005_r8*rxt(k,539) &
                      *y(k,240)
         mat(k,133) = .011_r8*rxt(k,544)*y(k,129) + .0512005_r8*rxt(k,543)*y(k,240)
         mat(k,2129) = .0443005_r8*rxt(k,518)*y(k,9) + .0076005_r8*rxt(k,528)*y(k,110) &
                      + .2157005_r8*rxt(k,539)*y(k,218) + .0512005_r8*rxt(k,543) &
                      *y(k,222)
         mat(k,1769) = .0085005_r8*rxt(k,537)*y(k,209)
         mat(k,110) = .0536005_r8*rxt(k,519)*y(k,129) + .1621005_r8*rxt(k,518) &
                      *y(k,240)
         mat(k,102) = .0166005_r8*rxt(k,529)*y(k,129) + .0113005_r8*rxt(k,528) &
                      *y(k,240)
         mat(k,1534) = .0536005_r8*rxt(k,519)*y(k,9) + .0166005_r8*rxt(k,529)*y(k,110) &
                      + .238_r8*rxt(k,540)*y(k,218) + .1185005_r8*rxt(k,544)*y(k,222)
         mat(k,122) = .0128005_r8*rxt(k,537)*y(k,251)
         mat(k,128) = .238_r8*rxt(k,540)*y(k,129) + .0738005_r8*rxt(k,539)*y(k,240)
         mat(k,134) = .1185005_r8*rxt(k,544)*y(k,129) + .1598005_r8*rxt(k,543) &
                      *y(k,240)
         mat(k,2130) = .1621005_r8*rxt(k,518)*y(k,9) + .0113005_r8*rxt(k,528)*y(k,110) &
                      + .0738005_r8*rxt(k,539)*y(k,218) + .1598005_r8*rxt(k,543) &
                      *y(k,222)
         mat(k,1770) = .0128005_r8*rxt(k,537)*y(k,209)
         mat(k,117) = -(rxt(k,536)*y(k,251))
         mat(k,1774) = -rxt(k,536)*y(k,208)
         mat(k,123) = -(rxt(k,537)*y(k,251))
         mat(k,1775) = -rxt(k,537)*y(k,209)
         mat(k,232) = .100_r8*rxt(k,448)*y(k,251)
         mat(k,250) = .230_r8*rxt(k,450)*y(k,251)
         mat(k,1789) = .100_r8*rxt(k,448)*y(k,217) + .230_r8*rxt(k,450)*y(k,220)
         mat(k,695) = -(rxt(k,472)*y(k,251))
         mat(k,1854) = -rxt(k,472)*y(k,211)
         mat(k,2167) = rxt(k,470)*y(k,255)
         mat(k,1148) = rxt(k,470)*y(k,240)
         mat(k,670) = -(rxt(k,473)*y(k,251))
         mat(k,1851) = -rxt(k,473)*y(k,212)
         mat(k,1559) = .200_r8*rxt(k,466)*y(k,249) + .200_r8*rxt(k,476)*y(k,256)
         mat(k,1623) = .500_r8*rxt(k,464)*y(k,249)
         mat(k,1098) = .200_r8*rxt(k,466)*y(k,129) + .500_r8*rxt(k,464)*y(k,235)
         mat(k,950) = .200_r8*rxt(k,476)*y(k,129)
         mat(k,524) = -(rxt(k,477)*y(k,251))
         mat(k,1833) = -rxt(k,477)*y(k,213)
         mat(k,2157) = rxt(k,475)*y(k,256)
         mat(k,949) = rxt(k,475)*y(k,240)
         mat(k,1016) = -(rxt(k,478)*y(k,131) + rxt(k,479)*y(k,251))
         mat(k,2045) = -rxt(k,478)*y(k,214)
         mat(k,1884) = -rxt(k,479)*y(k,214)
         mat(k,990) = .330_r8*rxt(k,459)*y(k,139)
         mat(k,1038) = .330_r8*rxt(k,462)*y(k,139)
         mat(k,1579) = .800_r8*rxt(k,466)*y(k,249) + .800_r8*rxt(k,476)*y(k,256)
         mat(k,2045) = mat(k,2045) + rxt(k,467)*y(k,249)
         mat(k,2291) = .330_r8*rxt(k,459)*y(k,6) + .330_r8*rxt(k,462)*y(k,116)
         mat(k,671) = rxt(k,473)*y(k,251)
         mat(k,1633) = .500_r8*rxt(k,464)*y(k,249) + rxt(k,474)*y(k,256)
         mat(k,1100) = .800_r8*rxt(k,466)*y(k,129) + rxt(k,467)*y(k,131) &
                      + .500_r8*rxt(k,464)*y(k,235)
         mat(k,1884) = mat(k,1884) + rxt(k,473)*y(k,212)
         mat(k,954) = .800_r8*rxt(k,476)*y(k,129) + rxt(k,474)*y(k,235)
         mat(k,1115) = -(rxt(k,480)*y(k,251))
         mat(k,1889) = -rxt(k,480)*y(k,215)
         mat(k,992) = .300_r8*rxt(k,459)*y(k,139)
         mat(k,1041) = .300_r8*rxt(k,462)*y(k,139)
         mat(k,1582) = .900_r8*rxt(k,471)*y(k,255)
         mat(k,2295) = .300_r8*rxt(k,459)*y(k,6) + .300_r8*rxt(k,462)*y(k,116)
         mat(k,1636) = rxt(k,469)*y(k,255)
         mat(k,1152) = .900_r8*rxt(k,471)*y(k,129) + rxt(k,469)*y(k,235)
         mat(k,648) = -(rxt(k,447)*y(k,251))
         mat(k,1848) = -rxt(k,447)*y(k,216)
         mat(k,2164) = rxt(k,445)*y(k,257)
         mat(k,769) = rxt(k,445)*y(k,240)
         mat(k,230) = -(rxt(k,448)*y(k,251))
         mat(k,1787) = -rxt(k,448)*y(k,217)
         mat(k,129) = -(rxt(k,539)*y(k,240) + rxt(k,540)*y(k,129))
         mat(k,2133) = -rxt(k,539)*y(k,218)
         mat(k,1537) = -rxt(k,540)*y(k,218)
         mat(k,229) = rxt(k,538)*y(k,251)
         mat(k,1776) = rxt(k,538)*y(k,217)
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
         mat(k,246) = -(rxt(k,414)*y(k,251))
         mat(k,1790) = -rxt(k,414)*y(k,219)
         mat(k,2136) = rxt(k,411)*y(k,258)
         mat(k,1206) = rxt(k,411)*y(k,240)
         mat(k,251) = -(rxt(k,450)*y(k,251))
         mat(k,1791) = -rxt(k,450)*y(k,220)
         mat(k,741) = -(rxt(k,453)*y(k,251))
         mat(k,1859) = -rxt(k,453)*y(k,221)
         mat(k,2172) = rxt(k,451)*y(k,259)
         mat(k,786) = rxt(k,451)*y(k,240)
         mat(k,135) = -(rxt(k,543)*y(k,240) + rxt(k,544)*y(k,129))
         mat(k,2134) = -rxt(k,543)*y(k,222)
         mat(k,1538) = -rxt(k,544)*y(k,222)
         mat(k,249) = rxt(k,542)*y(k,251)
         mat(k,1777) = rxt(k,542)*y(k,220)
         mat(k,259) = -(rxt(k,456)*y(k,251))
         mat(k,1792) = -rxt(k,456)*y(k,223)
         mat(k,252) = .150_r8*rxt(k,450)*y(k,251)
         mat(k,1792) = mat(k,1792) + .150_r8*rxt(k,450)*y(k,220)
         mat(k,469) = -(rxt(k,457)*y(k,251))
         mat(k,1826) = -rxt(k,457)*y(k,224)
         mat(k,2151) = rxt(k,454)*y(k,260)
         mat(k,540) = rxt(k,454)*y(k,240)
         mat(k,554) = -(rxt(k,415)*y(k,240) + rxt(k,416)*y(k,129) + rxt(k,444) &
                      *y(k,130))
         mat(k,2160) = -rxt(k,415)*y(k,227)
         mat(k,1554) = -rxt(k,416)*y(k,227)
         mat(k,1677) = -rxt(k,444)*y(k,227)
         mat(k,283) = rxt(k,421)*y(k,251)
         mat(k,1837) = rxt(k,421)*y(k,24)
         mat(k,939) = -(rxt(k,376)*y(k,240) + (rxt(k,377) + rxt(k,378)) * y(k,129))
         mat(k,2187) = -rxt(k,376)*y(k,228)
         mat(k,1574) = -(rxt(k,377) + rxt(k,378)) * y(k,228)
         mat(k,709) = rxt(k,379)*y(k,251)
         mat(k,274) = rxt(k,380)*y(k,251)
         mat(k,1878) = rxt(k,379)*y(k,2) + rxt(k,380)*y(k,17)
         mat(k,533) = -(rxt(k,418)*y(k,240) + rxt(k,419)*y(k,129))
         mat(k,2158) = -rxt(k,418)*y(k,229)
         mat(k,1551) = -rxt(k,419)*y(k,229)
         mat(k,210) = .350_r8*rxt(k,417)*y(k,251)
         mat(k,423) = rxt(k,420)*y(k,251)
         mat(k,1834) = .350_r8*rxt(k,417)*y(k,8) + rxt(k,420)*y(k,10)
         mat(k,477) = -(rxt(k,422)*y(k,240) + rxt(k,424)*y(k,129))
         mat(k,2152) = -rxt(k,422)*y(k,230)
         mat(k,1545) = -rxt(k,424)*y(k,230)
         mat(k,374) = rxt(k,423)*y(k,251)
         mat(k,233) = .070_r8*rxt(k,448)*y(k,251)
         mat(k,253) = .060_r8*rxt(k,450)*y(k,251)
         mat(k,1827) = rxt(k,423)*y(k,25) + .070_r8*rxt(k,448)*y(k,217) &
                      + .060_r8*rxt(k,450)*y(k,220)
         mat(k,860) = -(4._r8*rxt(k,299)*y(k,231) + rxt(k,300)*y(k,235) + rxt(k,301) &
                      *y(k,240) + rxt(k,302)*y(k,129))
         mat(k,1627) = -rxt(k,300)*y(k,231)
         mat(k,2183) = -rxt(k,301)*y(k,231)
         mat(k,1570) = -rxt(k,302)*y(k,231)
         mat(k,379) = .500_r8*rxt(k,304)*y(k,251)
         mat(k,331) = rxt(k,305)*y(k,58) + rxt(k,306)*y(k,251)
         mat(k,2242) = rxt(k,305)*y(k,30)
         mat(k,1870) = .500_r8*rxt(k,304)*y(k,29) + rxt(k,306)*y(k,30)
         mat(k,831) = -(rxt(k,328)*y(k,235) + rxt(k,329)*y(k,240) + rxt(k,330) &
                      *y(k,129))
         mat(k,1626) = -rxt(k,328)*y(k,232)
         mat(k,2180) = -rxt(k,329)*y(k,232)
         mat(k,1569) = -rxt(k,330)*y(k,232)
         mat(k,428) = rxt(k,331)*y(k,251)
         mat(k,149) = rxt(k,332)*y(k,251)
         mat(k,1867) = rxt(k,331)*y(k,32) + rxt(k,332)*y(k,33)
         mat(k,678) = -(rxt(k,425)*y(k,240) + rxt(k,426)*y(k,129))
         mat(k,2166) = -rxt(k,425)*y(k,233)
         mat(k,1560) = -rxt(k,426)*y(k,233)
         mat(k,309) = rxt(k,427)*y(k,251)
         mat(k,1560) = mat(k,1560) + rxt(k,416)*y(k,227)
         mat(k,2280) = rxt(k,442)*y(k,146)
         mat(k,514) = rxt(k,442)*y(k,139)
         mat(k,555) = rxt(k,416)*y(k,129) + .400_r8*rxt(k,415)*y(k,240)
         mat(k,2166) = mat(k,2166) + .400_r8*rxt(k,415)*y(k,227)
         mat(k,1852) = rxt(k,427)*y(k,34)
         mat(k,1426) = -(4._r8*rxt(k,310)*y(k,234) + rxt(k,311)*y(k,235) + rxt(k,312) &
                      *y(k,240) + rxt(k,313)*y(k,129) + rxt(k,324)*y(k,130) + rxt(k,351) &
                      *y(k,244) + rxt(k,384)*y(k,242) + rxt(k,389)*y(k,243) + rxt(k,398) &
                      *y(k,103) + rxt(k,409)*y(k,258))
         mat(k,1652) = -rxt(k,311)*y(k,234)
         mat(k,2210) = -rxt(k,312)*y(k,234)
         mat(k,1599) = -rxt(k,313)*y(k,234)
         mat(k,1695) = -rxt(k,324)*y(k,234)
         mat(k,1353) = -rxt(k,351)*y(k,234)
         mat(k,1300) = -rxt(k,384)*y(k,234)
         mat(k,1332) = -rxt(k,389)*y(k,234)
         mat(k,1237) = -rxt(k,398)*y(k,234)
         mat(k,1215) = -rxt(k,409)*y(k,234)
         mat(k,997) = .060_r8*rxt(k,459)*y(k,139)
         mat(k,1127) = rxt(k,307)*y(k,131) + rxt(k,308)*y(k,251)
         mat(k,1262) = rxt(k,333)*y(k,131) + rxt(k,334)*y(k,251)
         mat(k,625) = .500_r8*rxt(k,315)*y(k,251)
         mat(k,883) = .080_r8*rxt(k,404)*y(k,139)
         mat(k,1253) = .100_r8*rxt(k,357)*y(k,139)
         mat(k,1047) = .060_r8*rxt(k,462)*y(k,139)
         mat(k,1374) = .280_r8*rxt(k,371)*y(k,139)
         mat(k,1599) = mat(k,1599) + .530_r8*rxt(k,355)*y(k,244) + rxt(k,364)*y(k,246) &
                      + rxt(k,367)*y(k,248) + rxt(k,342)*y(k,254)
         mat(k,2068) = rxt(k,307)*y(k,47) + rxt(k,333)*y(k,51) + .530_r8*rxt(k,354) &
                      *y(k,244) + rxt(k,365)*y(k,246)
         mat(k,2311) = .060_r8*rxt(k,459)*y(k,6) + .080_r8*rxt(k,404)*y(k,100) &
                      + .100_r8*rxt(k,357)*y(k,111) + .060_r8*rxt(k,462)*y(k,116) &
                      + .280_r8*rxt(k,371)*y(k,118)
         mat(k,1118) = .650_r8*rxt(k,480)*y(k,251)
         mat(k,1426) = mat(k,1426) + .530_r8*rxt(k,351)*y(k,244)
         mat(k,1652) = mat(k,1652) + .260_r8*rxt(k,352)*y(k,244) + rxt(k,361)*y(k,246) &
                      + .300_r8*rxt(k,340)*y(k,254)
         mat(k,2210) = mat(k,2210) + .450_r8*rxt(k,362)*y(k,246) + .200_r8*rxt(k,366) &
                      *y(k,248) + .150_r8*rxt(k,341)*y(k,254)
         mat(k,1353) = mat(k,1353) + .530_r8*rxt(k,355)*y(k,129) + .530_r8*rxt(k,354) &
                      *y(k,131) + .530_r8*rxt(k,351)*y(k,234) + .260_r8*rxt(k,352) &
                      *y(k,235)
         mat(k,1395) = rxt(k,364)*y(k,129) + rxt(k,365)*y(k,131) + rxt(k,361)*y(k,235) &
                      + .450_r8*rxt(k,362)*y(k,240) + 4.000_r8*rxt(k,363)*y(k,246)
         mat(k,719) = rxt(k,367)*y(k,129) + .200_r8*rxt(k,366)*y(k,240)
         mat(k,1907) = rxt(k,308)*y(k,47) + rxt(k,334)*y(k,51) + .500_r8*rxt(k,315) &
                      *y(k,53) + .650_r8*rxt(k,480)*y(k,215)
         mat(k,1172) = rxt(k,342)*y(k,129) + .300_r8*rxt(k,340)*y(k,235) &
                      + .150_r8*rxt(k,341)*y(k,240)
         mat(k,1657) = -(rxt(k,201)*y(k,61) + (4._r8*rxt(k,278) + 4._r8*rxt(k,279) &
                      ) * y(k,235) + rxt(k,280)*y(k,240) + rxt(k,281)*y(k,129) &
                      + rxt(k,300)*y(k,231) + rxt(k,311)*y(k,234) + rxt(k,328) &
                      *y(k,232) + rxt(k,340)*y(k,254) + rxt(k,352)*y(k,244) + rxt(k,361) &
                      *y(k,246) + rxt(k,385)*y(k,242) + rxt(k,390)*y(k,243) + rxt(k,399) &
                      *y(k,103) + rxt(k,410)*y(k,258) + rxt(k,464)*y(k,249) + rxt(k,469) &
                      *y(k,255) + rxt(k,474)*y(k,256))
         mat(k,1940) = -rxt(k,201)*y(k,235)
         mat(k,2217) = -rxt(k,280)*y(k,235)
         mat(k,1605) = -rxt(k,281)*y(k,235)
         mat(k,863) = -rxt(k,300)*y(k,235)
         mat(k,1430) = -rxt(k,311)*y(k,235)
         mat(k,835) = -rxt(k,328)*y(k,235)
         mat(k,1174) = -rxt(k,340)*y(k,235)
         mat(k,1356) = -rxt(k,352)*y(k,235)
         mat(k,1398) = -rxt(k,361)*y(k,235)
         mat(k,1303) = -rxt(k,385)*y(k,235)
         mat(k,1335) = -rxt(k,390)*y(k,235)
         mat(k,1240) = -rxt(k,399)*y(k,235)
         mat(k,1217) = -rxt(k,410)*y(k,235)
         mat(k,1107) = -rxt(k,464)*y(k,235)
         mat(k,1160) = -rxt(k,469)*y(k,235)
         mat(k,957) = -rxt(k,474)*y(k,235)
         mat(k,1072) = .280_r8*rxt(k,327)*y(k,139)
         mat(k,726) = rxt(k,314)*y(k,251)
         mat(k,453) = .700_r8*rxt(k,283)*y(k,251)
         mat(k,1473) = rxt(k,195)*y(k,58) + rxt(k,251)*y(k,75) + rxt(k,290)*y(k,250) &
                      + rxt(k,284)*y(k,251)
         mat(k,2256) = rxt(k,195)*y(k,56)
         mat(k,918) = rxt(k,251)*y(k,56)
         mat(k,884) = .050_r8*rxt(k,404)*y(k,139)
         mat(k,1240) = mat(k,1240) + rxt(k,398)*y(k,234)
         mat(k,1605) = mat(k,1605) + rxt(k,313)*y(k,234) + .830_r8*rxt(k,430)*y(k,236) &
                      + .170_r8*rxt(k,436)*y(k,247)
         mat(k,2317) = .280_r8*rxt(k,327)*y(k,31) + .050_r8*rxt(k,404)*y(k,100)
         mat(k,1430) = mat(k,1430) + rxt(k,398)*y(k,103) + rxt(k,313)*y(k,129) &
                      + 4.000_r8*rxt(k,310)*y(k,234) + .900_r8*rxt(k,311)*y(k,235) &
                      + .450_r8*rxt(k,312)*y(k,240) + rxt(k,384)*y(k,242) + rxt(k,389) &
                      *y(k,243) + rxt(k,351)*y(k,244) + rxt(k,360)*y(k,246) &
                      + rxt(k,409)*y(k,258)
         mat(k,1657) = mat(k,1657) + .900_r8*rxt(k,311)*y(k,234)
         mat(k,803) = .830_r8*rxt(k,430)*y(k,129) + .330_r8*rxt(k,429)*y(k,240)
         mat(k,2217) = mat(k,2217) + .450_r8*rxt(k,312)*y(k,234) + .330_r8*rxt(k,429) &
                      *y(k,236) + .070_r8*rxt(k,435)*y(k,247)
         mat(k,1303) = mat(k,1303) + rxt(k,384)*y(k,234)
         mat(k,1335) = mat(k,1335) + rxt(k,389)*y(k,234)
         mat(k,1356) = mat(k,1356) + rxt(k,351)*y(k,234)
         mat(k,1398) = mat(k,1398) + rxt(k,360)*y(k,234)
         mat(k,908) = .170_r8*rxt(k,436)*y(k,129) + .070_r8*rxt(k,435)*y(k,240)
         mat(k,1742) = rxt(k,290)*y(k,56)
         mat(k,1914) = rxt(k,314)*y(k,52) + .700_r8*rxt(k,283)*y(k,55) + rxt(k,284) &
                      *y(k,56)
         mat(k,1217) = mat(k,1217) + rxt(k,409)*y(k,234)
         mat(k,799) = -(rxt(k,429)*y(k,240) + rxt(k,430)*y(k,129) + rxt(k,431) &
                      *y(k,130))
         mat(k,2177) = -rxt(k,429)*y(k,236)
         mat(k,1567) = -rxt(k,430)*y(k,236)
         mat(k,1683) = -rxt(k,431)*y(k,236)
         mat(k,606) = -((rxt(k,348) + rxt(k,349)) * y(k,129))
         mat(k,1556) = -(rxt(k,348) + rxt(k,349)) * y(k,237)
         mat(k,400) = rxt(k,347)*y(k,251)
         mat(k,1843) = rxt(k,347)*y(k,18)
         mat(k,1541) = .750_r8*rxt(k,317)*y(k,239)
         mat(k,753) = .750_r8*rxt(k,317)*y(k,129)
         mat(k,754) = -(rxt(k,316)*y(k,240) + rxt(k,317)*y(k,129))
         mat(k,2173) = -rxt(k,316)*y(k,239)
         mat(k,1563) = -rxt(k,317)*y(k,239)
         mat(k,599) = rxt(k,323)*y(k,251)
         mat(k,1860) = rxt(k,323)*y(k,27)
         mat(k,2227) = -((rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,78) + rxt(k,158) &
                      *y(k,138) + rxt(k,159)*y(k,139) + rxt(k,163)*y(k,251) &
                      + 4._r8*rxt(k,168)*y(k,240) + rxt(k,178)*y(k,131) + rxt(k,183) &
                      *y(k,129) + rxt(k,188)*y(k,130) + (rxt(k,198) + rxt(k,199) &
                      ) * y(k,58) + rxt(k,205)*y(k,61) + rxt(k,231)*y(k,19) + rxt(k,237) &
                      *y(k,21) + rxt(k,274)*y(k,44) + rxt(k,280)*y(k,235) + rxt(k,287) &
                      *y(k,241) + rxt(k,301)*y(k,231) + rxt(k,312)*y(k,234) + rxt(k,316) &
                      *y(k,239) + rxt(k,329)*y(k,232) + rxt(k,337)*y(k,253) + rxt(k,341) &
                      *y(k,254) + rxt(k,353)*y(k,244) + rxt(k,362)*y(k,246) + rxt(k,366) &
                      *y(k,248) + rxt(k,376)*y(k,228) + rxt(k,386)*y(k,242) + rxt(k,391) &
                      *y(k,243) + rxt(k,400)*y(k,103) + rxt(k,411)*y(k,258) + rxt(k,415) &
                      *y(k,227) + rxt(k,418)*y(k,229) + rxt(k,422)*y(k,230) + rxt(k,425) &
                      *y(k,233) + rxt(k,429)*y(k,236) + rxt(k,432)*y(k,245) + rxt(k,435) &
                      *y(k,247) + rxt(k,438)*y(k,252) + rxt(k,445)*y(k,257) + rxt(k,451) &
                      *y(k,259) + rxt(k,454)*y(k,260) + rxt(k,465)*y(k,249) + rxt(k,470) &
                      *y(k,255) + rxt(k,475)*y(k,256))
         mat(k,1512) = -(rxt(k,154) + rxt(k,155) + rxt(k,156)) * y(k,240)
         mat(k,2003) = -rxt(k,158)*y(k,240)
         mat(k,2327) = -rxt(k,159)*y(k,240)
         mat(k,1924) = -rxt(k,163)*y(k,240)
         mat(k,2084) = -rxt(k,178)*y(k,240)
         mat(k,1615) = -rxt(k,183)*y(k,240)
         mat(k,1711) = -rxt(k,188)*y(k,240)
         mat(k,2266) = -(rxt(k,198) + rxt(k,199)) * y(k,240)
         mat(k,1950) = -rxt(k,205)*y(k,240)
         mat(k,1462) = -rxt(k,231)*y(k,240)
         mat(k,2108) = -rxt(k,237)*y(k,240)
         mat(k,2027) = -rxt(k,274)*y(k,240)
         mat(k,1667) = -rxt(k,280)*y(k,240)
         mat(k,494) = -rxt(k,287)*y(k,240)
         mat(k,867) = -rxt(k,301)*y(k,240)
         mat(k,1436) = -rxt(k,312)*y(k,240)
         mat(k,760) = -rxt(k,316)*y(k,240)
         mat(k,839) = -rxt(k,329)*y(k,240)
         mat(k,815) = -rxt(k,337)*y(k,240)
         mat(k,1178) = -rxt(k,341)*y(k,240)
         mat(k,1362) = -rxt(k,353)*y(k,240)
         mat(k,1404) = -rxt(k,362)*y(k,240)
         mat(k,723) = -rxt(k,366)*y(k,240)
         mat(k,948) = -rxt(k,376)*y(k,240)
         mat(k,1309) = -rxt(k,386)*y(k,240)
         mat(k,1341) = -rxt(k,391)*y(k,240)
         mat(k,1246) = -rxt(k,400)*y(k,240)
         mat(k,1223) = -rxt(k,411)*y(k,240)
         mat(k,559) = -rxt(k,415)*y(k,240)
         mat(k,539) = -rxt(k,418)*y(k,240)
         mat(k,482) = -rxt(k,422)*y(k,240)
         mat(k,682) = -rxt(k,425)*y(k,240)
         mat(k,806) = -rxt(k,429)*y(k,240)
         mat(k,766) = -rxt(k,432)*y(k,240)
         mat(k,911) = -rxt(k,435)*y(k,240)
         mat(k,501) = -rxt(k,438)*y(k,240)
         mat(k,781) = -rxt(k,445)*y(k,240)
         mat(k,798) = -rxt(k,451)*y(k,240)
         mat(k,547) = -rxt(k,454)*y(k,240)
         mat(k,1112) = -rxt(k,465)*y(k,240)
         mat(k,1165) = -rxt(k,470)*y(k,240)
         mat(k,961) = -rxt(k,475)*y(k,240)
         mat(k,1006) = .570_r8*rxt(k,459)*y(k,139)
         mat(k,212) = .650_r8*rxt(k,417)*y(k,251)
         mat(k,1462) = mat(k,1462) + rxt(k,230)*y(k,44)
         mat(k,2108) = mat(k,2108) + rxt(k,242)*y(k,251)
         mat(k,329) = .350_r8*rxt(k,296)*y(k,251)
         mat(k,604) = .130_r8*rxt(k,298)*y(k,139)
         mat(k,306) = rxt(k,303)*y(k,251)
         mat(k,1078) = .280_r8*rxt(k,327)*y(k,139)
         mat(k,2027) = mat(k,2027) + rxt(k,230)*y(k,19) + rxt(k,194)*y(k,58) &
                      + rxt(k,275)*y(k,131) + rxt(k,276)*y(k,138)
         mat(k,620) = rxt(k,259)*y(k,58) + rxt(k,260)*y(k,251)
         mat(k,412) = rxt(k,262)*y(k,58) + rxt(k,263)*y(k,251)
         mat(k,144) = rxt(k,309)*y(k,251)
         mat(k,829) = rxt(k,282)*y(k,251)
         mat(k,1480) = rxt(k,291)*y(k,250)
         mat(k,2266) = mat(k,2266) + rxt(k,194)*y(k,44) + rxt(k,259)*y(k,45) &
                      + rxt(k,262)*y(k,48) + rxt(k,197)*y(k,81)
         mat(k,1950) = mat(k,1950) + rxt(k,201)*y(k,235) + rxt(k,212)*y(k,251)
         mat(k,1138) = rxt(k,294)*y(k,251)
         mat(k,241) = .730_r8*rxt(k,428)*y(k,251)
         mat(k,323) = .500_r8*rxt(k,496)*y(k,251)
         mat(k,1146) = rxt(k,320)*y(k,251)
         mat(k,974) = rxt(k,321)*y(k,251)
         mat(k,643) = rxt(k,197)*y(k,58) + rxt(k,153)*y(k,138) + rxt(k,162)*y(k,251)
         mat(k,228) = rxt(k,285)*y(k,251)
         mat(k,1014) = rxt(k,286)*y(k,251)
         mat(k,1204) = rxt(k,350)*y(k,251)
         mat(k,1185) = rxt(k,335)*y(k,251)
         mat(k,888) = .370_r8*rxt(k,404)*y(k,139)
         mat(k,665) = .300_r8*rxt(k,395)*y(k,251)
         mat(k,589) = rxt(k,396)*y(k,251)
         mat(k,1246) = mat(k,1246) + rxt(k,401)*y(k,129) + rxt(k,402)*y(k,131) &
                      + rxt(k,398)*y(k,234) + 1.200_r8*rxt(k,399)*y(k,235)
         mat(k,438) = rxt(k,403)*y(k,251)
         mat(k,1257) = .140_r8*rxt(k,357)*y(k,139)
         mat(k,390) = .200_r8*rxt(k,359)*y(k,251)
         mat(k,638) = .500_r8*rxt(k,370)*y(k,251)
         mat(k,1056) = .570_r8*rxt(k,462)*y(k,139)
         mat(k,1384) = .280_r8*rxt(k,371)*y(k,139)
         mat(k,444) = rxt(k,407)*y(k,251)
         mat(k,1096) = rxt(k,408)*y(k,251)
         mat(k,1615) = mat(k,1615) + rxt(k,401)*y(k,103) + rxt(k,377)*y(k,228) &
                      + rxt(k,419)*y(k,229) + rxt(k,424)*y(k,230) + rxt(k,302) &
                      *y(k,231) + rxt(k,330)*y(k,232) + rxt(k,281)*y(k,235) &
                      + .170_r8*rxt(k,430)*y(k,236) + rxt(k,348)*y(k,237) &
                      + .250_r8*rxt(k,317)*y(k,239) + rxt(k,289)*y(k,241) &
                      + .920_r8*rxt(k,387)*y(k,242) + .920_r8*rxt(k,393)*y(k,243) &
                      + .470_r8*rxt(k,355)*y(k,244) + .400_r8*rxt(k,433)*y(k,245) &
                      + .830_r8*rxt(k,436)*y(k,247) + rxt(k,439)*y(k,252) + rxt(k,338) &
                      *y(k,253) + .900_r8*rxt(k,471)*y(k,255) + .800_r8*rxt(k,476) &
                      *y(k,256) + rxt(k,446)*y(k,257) + rxt(k,412)*y(k,258) &
                      + rxt(k,452)*y(k,259) + rxt(k,455)*y(k,260)
         mat(k,2084) = mat(k,2084) + rxt(k,275)*y(k,44) + rxt(k,402)*y(k,103) &
                      + rxt(k,388)*y(k,242) + rxt(k,394)*y(k,243) + .470_r8*rxt(k,354) &
                      *y(k,244) + rxt(k,181)*y(k,251) + rxt(k,413)*y(k,258)
         mat(k,2003) = mat(k,2003) + rxt(k,276)*y(k,44) + rxt(k,153)*y(k,81)
         mat(k,2327) = mat(k,2327) + .570_r8*rxt(k,459)*y(k,6) + .130_r8*rxt(k,298) &
                      *y(k,27) + .280_r8*rxt(k,327)*y(k,31) + .370_r8*rxt(k,404) &
                      *y(k,100) + .140_r8*rxt(k,357)*y(k,111) + .570_r8*rxt(k,462) &
                      *y(k,116) + .280_r8*rxt(k,371)*y(k,118) + rxt(k,165)*y(k,251)
         mat(k,221) = .800_r8*rxt(k,440)*y(k,251)
         mat(k,901) = rxt(k,486)*y(k,251)
         mat(k,1123) = .200_r8*rxt(k,480)*y(k,251)
         mat(k,236) = .280_r8*rxt(k,448)*y(k,251)
         mat(k,258) = .380_r8*rxt(k,450)*y(k,251)
         mat(k,263) = .630_r8*rxt(k,456)*y(k,251)
         mat(k,948) = mat(k,948) + rxt(k,377)*y(k,129)
         mat(k,539) = mat(k,539) + rxt(k,419)*y(k,129)
         mat(k,482) = mat(k,482) + rxt(k,424)*y(k,129)
         mat(k,867) = mat(k,867) + rxt(k,302)*y(k,129) + 2.400_r8*rxt(k,299)*y(k,231) &
                      + rxt(k,300)*y(k,235)
         mat(k,839) = mat(k,839) + rxt(k,330)*y(k,129) + rxt(k,328)*y(k,235)
         mat(k,1436) = mat(k,1436) + rxt(k,398)*y(k,103) + .900_r8*rxt(k,311)*y(k,235) &
                      + rxt(k,384)*y(k,242) + rxt(k,389)*y(k,243) + .470_r8*rxt(k,351) &
                      *y(k,244) + rxt(k,409)*y(k,258)
         mat(k,1667) = mat(k,1667) + rxt(k,201)*y(k,61) + 1.200_r8*rxt(k,399)*y(k,103) &
                      + rxt(k,281)*y(k,129) + rxt(k,300)*y(k,231) + rxt(k,328) &
                      *y(k,232) + .900_r8*rxt(k,311)*y(k,234) + 4.000_r8*rxt(k,278) &
                      *y(k,235) + rxt(k,385)*y(k,242) + rxt(k,390)*y(k,243) &
                      + .730_r8*rxt(k,352)*y(k,244) + rxt(k,361)*y(k,246) &
                      + .500_r8*rxt(k,464)*y(k,249) + .300_r8*rxt(k,340)*y(k,254) &
                      + rxt(k,469)*y(k,255) + rxt(k,474)*y(k,256) + .800_r8*rxt(k,410) &
                      *y(k,258)
         mat(k,806) = mat(k,806) + .170_r8*rxt(k,430)*y(k,129) + .070_r8*rxt(k,429) &
                      *y(k,240)
         mat(k,613) = rxt(k,348)*y(k,129)
         mat(k,760) = mat(k,760) + .250_r8*rxt(k,317)*y(k,129)
         mat(k,2227) = mat(k,2227) + .070_r8*rxt(k,429)*y(k,236) + .160_r8*rxt(k,432) &
                      *y(k,245) + .330_r8*rxt(k,435)*y(k,247)
         mat(k,494) = mat(k,494) + rxt(k,289)*y(k,129)
         mat(k,1309) = mat(k,1309) + .920_r8*rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131) &
                      + rxt(k,384)*y(k,234) + rxt(k,385)*y(k,235)
         mat(k,1341) = mat(k,1341) + .920_r8*rxt(k,393)*y(k,129) + rxt(k,394)*y(k,131) &
                      + rxt(k,389)*y(k,234) + rxt(k,390)*y(k,235)
         mat(k,1362) = mat(k,1362) + .470_r8*rxt(k,355)*y(k,129) + .470_r8*rxt(k,354) &
                      *y(k,131) + .470_r8*rxt(k,351)*y(k,234) + .730_r8*rxt(k,352) &
                      *y(k,235)
         mat(k,766) = mat(k,766) + .400_r8*rxt(k,433)*y(k,129) + .160_r8*rxt(k,432) &
                      *y(k,240)
         mat(k,1404) = mat(k,1404) + rxt(k,361)*y(k,235)
         mat(k,911) = mat(k,911) + .830_r8*rxt(k,436)*y(k,129) + .330_r8*rxt(k,435) &
                      *y(k,240)
         mat(k,1112) = mat(k,1112) + .500_r8*rxt(k,464)*y(k,235)
         mat(k,1752) = rxt(k,291)*y(k,56)
         mat(k,1924) = mat(k,1924) + .650_r8*rxt(k,417)*y(k,8) + rxt(k,242)*y(k,21) &
                      + .350_r8*rxt(k,296)*y(k,26) + rxt(k,303)*y(k,28) + rxt(k,260) &
                      *y(k,45) + rxt(k,263)*y(k,48) + rxt(k,309)*y(k,49) + rxt(k,282) &
                      *y(k,54) + rxt(k,212)*y(k,61) + rxt(k,294)*y(k,64) &
                      + .730_r8*rxt(k,428)*y(k,68) + .500_r8*rxt(k,496)*y(k,69) &
                      + rxt(k,320)*y(k,76) + rxt(k,321)*y(k,77) + rxt(k,162)*y(k,81) &
                      + rxt(k,285)*y(k,88) + rxt(k,286)*y(k,89) + rxt(k,350)*y(k,95) &
                      + rxt(k,335)*y(k,97) + .300_r8*rxt(k,395)*y(k,101) + rxt(k,396) &
                      *y(k,102) + rxt(k,403)*y(k,104) + .200_r8*rxt(k,359)*y(k,112) &
                      + .500_r8*rxt(k,370)*y(k,115) + rxt(k,407)*y(k,122) + rxt(k,408) &
                      *y(k,123) + rxt(k,181)*y(k,131) + rxt(k,165)*y(k,139) &
                      + .800_r8*rxt(k,440)*y(k,147) + rxt(k,486)*y(k,158) &
                      + .200_r8*rxt(k,480)*y(k,215) + .280_r8*rxt(k,448)*y(k,217) &
                      + .380_r8*rxt(k,450)*y(k,220) + .630_r8*rxt(k,456)*y(k,223)
         mat(k,501) = mat(k,501) + rxt(k,439)*y(k,129)
         mat(k,815) = mat(k,815) + rxt(k,338)*y(k,129)
         mat(k,1178) = mat(k,1178) + .300_r8*rxt(k,340)*y(k,235)
         mat(k,1165) = mat(k,1165) + .900_r8*rxt(k,471)*y(k,129) + rxt(k,469)*y(k,235)
         mat(k,961) = mat(k,961) + .800_r8*rxt(k,476)*y(k,129) + rxt(k,474)*y(k,235)
         mat(k,781) = mat(k,781) + rxt(k,446)*y(k,129)
         mat(k,1223) = mat(k,1223) + rxt(k,412)*y(k,129) + rxt(k,413)*y(k,131) &
                      + rxt(k,409)*y(k,234) + .800_r8*rxt(k,410)*y(k,235)
         mat(k,798) = mat(k,798) + rxt(k,452)*y(k,129)
         mat(k,547) = mat(k,547) + rxt(k,455)*y(k,129)
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
         mat(k,489) = -(rxt(k,287)*y(k,240) + rxt(k,289)*y(k,129))
         mat(k,2153) = -rxt(k,287)*y(k,241)
         mat(k,1546) = -rxt(k,289)*y(k,241)
         mat(k,2007) = rxt(k,274)*y(k,240)
         mat(k,2153) = mat(k,2153) + rxt(k,274)*y(k,44)
         mat(k,1296) = -(rxt(k,384)*y(k,234) + rxt(k,385)*y(k,235) + rxt(k,386) &
                      *y(k,240) + rxt(k,387)*y(k,129) + rxt(k,388)*y(k,131))
         mat(k,1421) = -rxt(k,384)*y(k,242)
         mat(k,1647) = -rxt(k,385)*y(k,242)
         mat(k,2205) = -rxt(k,386)*y(k,242)
         mat(k,1594) = -rxt(k,387)*y(k,242)
         mat(k,2063) = -rxt(k,388)*y(k,242)
         mat(k,880) = .600_r8*rxt(k,405)*y(k,251)
         mat(k,1902) = .600_r8*rxt(k,405)*y(k,100)
         mat(k,1328) = -(rxt(k,389)*y(k,234) + rxt(k,390)*y(k,235) + rxt(k,391) &
                      *y(k,240) + rxt(k,393)*y(k,129) + rxt(k,394)*y(k,131))
         mat(k,1422) = -rxt(k,389)*y(k,243)
         mat(k,1648) = -rxt(k,390)*y(k,243)
         mat(k,2206) = -rxt(k,391)*y(k,243)
         mat(k,1595) = -rxt(k,393)*y(k,243)
         mat(k,2064) = -rxt(k,394)*y(k,243)
         mat(k,881) = .400_r8*rxt(k,405)*y(k,251)
         mat(k,1903) = .400_r8*rxt(k,405)*y(k,100)
         mat(k,1351) = -(rxt(k,351)*y(k,234) + rxt(k,352)*y(k,235) + rxt(k,353) &
                      *y(k,240) + rxt(k,354)*y(k,131) + (rxt(k,355) + rxt(k,356) &
                      ) * y(k,129))
         mat(k,1423) = -rxt(k,351)*y(k,244)
         mat(k,1649) = -rxt(k,352)*y(k,244)
         mat(k,2207) = -rxt(k,353)*y(k,244)
         mat(k,2065) = -rxt(k,354)*y(k,244)
         mat(k,1596) = -(rxt(k,355) + rxt(k,356)) * y(k,244)
         mat(k,1251) = .500_r8*rxt(k,358)*y(k,251)
         mat(k,387) = .200_r8*rxt(k,359)*y(k,251)
         mat(k,1371) = rxt(k,372)*y(k,251)
         mat(k,1904) = .500_r8*rxt(k,358)*y(k,111) + .200_r8*rxt(k,359)*y(k,112) &
                      + rxt(k,372)*y(k,118)
         mat(k,761) = -(rxt(k,432)*y(k,240) + rxt(k,433)*y(k,129) + rxt(k,434) &
                      *y(k,130))
         mat(k,2174) = -rxt(k,432)*y(k,245)
         mat(k,1564) = -rxt(k,433)*y(k,245)
         mat(k,1682) = -rxt(k,434)*y(k,245)
         mat(k,1394) = -(rxt(k,360)*y(k,234) + rxt(k,361)*y(k,235) + rxt(k,362) &
                      *y(k,240) + 4._r8*rxt(k,363)*y(k,246) + rxt(k,364)*y(k,129) &
                      + rxt(k,365)*y(k,131) + rxt(k,373)*y(k,130))
         mat(k,1425) = -rxt(k,360)*y(k,246)
         mat(k,1651) = -rxt(k,361)*y(k,246)
         mat(k,2209) = -rxt(k,362)*y(k,246)
         mat(k,1598) = -rxt(k,364)*y(k,246)
         mat(k,2067) = -rxt(k,365)*y(k,246)
         mat(k,1694) = -rxt(k,373)*y(k,246)
         mat(k,1252) = .500_r8*rxt(k,358)*y(k,251)
         mat(k,388) = .500_r8*rxt(k,359)*y(k,251)
         mat(k,1906) = .500_r8*rxt(k,358)*y(k,111) + .500_r8*rxt(k,359)*y(k,112)
         mat(k,903) = -(rxt(k,435)*y(k,240) + rxt(k,436)*y(k,129) + rxt(k,437) &
                      *y(k,130))
         mat(k,2186) = -rxt(k,435)*y(k,247)
         mat(k,1573) = -rxt(k,436)*y(k,247)
         mat(k,1687) = -rxt(k,437)*y(k,247)
         mat(k,717) = -(rxt(k,366)*y(k,240) + rxt(k,367)*y(k,129))
         mat(k,2169) = -rxt(k,366)*y(k,248)
         mat(k,1562) = -rxt(k,367)*y(k,248)
         mat(k,549) = rxt(k,368)*y(k,251)
         mat(k,353) = rxt(k,369)*y(k,251)
         mat(k,1856) = rxt(k,368)*y(k,113) + rxt(k,369)*y(k,114)
         mat(k,1101) = -(rxt(k,464)*y(k,235) + rxt(k,465)*y(k,240) + rxt(k,466) &
                      *y(k,129) + rxt(k,467)*y(k,131))
         mat(k,1635) = -rxt(k,464)*y(k,249)
         mat(k,2193) = -rxt(k,465)*y(k,249)
         mat(k,1581) = -rxt(k,466)*y(k,249)
         mat(k,2049) = -rxt(k,467)*y(k,249)
         mat(k,991) = rxt(k,458)*y(k,131)
         mat(k,1040) = rxt(k,461)*y(k,131)
         mat(k,2049) = mat(k,2049) + rxt(k,458)*y(k,6) + rxt(k,461)*y(k,116) &
                      + .500_r8*rxt(k,478)*y(k,214)
         mat(k,465) = rxt(k,468)*y(k,251)
         mat(k,1017) = .500_r8*rxt(k,478)*y(k,131)
         mat(k,1888) = rxt(k,468)*y(k,133)
         mat(k,1744) = -(rxt(k,144)*y(k,79) + rxt(k,145)*y(k,261) + rxt(k,148) &
                      *y(k,139) + (rxt(k,186) + rxt(k,187)) * y(k,120) + rxt(k,219) &
                      *y(k,35) + rxt(k,220)*y(k,36) + rxt(k,221)*y(k,38) + rxt(k,222) &
                      *y(k,39) + rxt(k,223)*y(k,40) + rxt(k,224)*y(k,41) + rxt(k,225) &
                      *y(k,42) + (rxt(k,226) + rxt(k,227)) * y(k,87) + rxt(k,246) &
                      *y(k,37) + rxt(k,247)*y(k,57) + rxt(k,248)*y(k,80) + (rxt(k,249) &
                      + rxt(k,250)) * y(k,83) + rxt(k,255)*y(k,66) + rxt(k,256) &
                      *y(k,67) + rxt(k,269)*y(k,43) + rxt(k,270)*y(k,45) + rxt(k,271) &
                      *y(k,84) + rxt(k,272)*y(k,85) + rxt(k,273)*y(k,86) + (rxt(k,290) &
                      + rxt(k,291) + rxt(k,292)) * y(k,56) + rxt(k,293)*y(k,88))
         mat(k,1446) = -rxt(k,144)*y(k,250)
         mat(k,2345) = -rxt(k,145)*y(k,250)
         mat(k,2319) = -rxt(k,148)*y(k,250)
         mat(k,224) = -(rxt(k,186) + rxt(k,187)) * y(k,250)
         mat(k,140) = -rxt(k,219)*y(k,250)
         mat(k,181) = -rxt(k,220)*y(k,250)
         mat(k,154) = -rxt(k,221)*y(k,250)
         mat(k,191) = -rxt(k,222)*y(k,250)
         mat(k,158) = -rxt(k,223)*y(k,250)
         mat(k,196) = -rxt(k,224)*y(k,250)
         mat(k,162) = -rxt(k,225)*y(k,250)
         mat(k,1964) = -(rxt(k,226) + rxt(k,227)) * y(k,250)
         mat(k,187) = -rxt(k,246)*y(k,250)
         mat(k,459) = -rxt(k,247)*y(k,250)
         mat(k,171) = -rxt(k,248)*y(k,250)
         mat(k,852) = -(rxt(k,249) + rxt(k,250)) * y(k,250)
         mat(k,292) = -rxt(k,255)*y(k,250)
         mat(k,271) = -rxt(k,256)*y(k,250)
         mat(k,519) = -rxt(k,269)*y(k,250)
         mat(k,617) = -rxt(k,270)*y(k,250)
         mat(k,266) = -rxt(k,271)*y(k,250)
         mat(k,287) = -rxt(k,272)*y(k,250)
         mat(k,342) = -rxt(k,273)*y(k,250)
         mat(k,1474) = -(rxt(k,290) + rxt(k,291) + rxt(k,292)) * y(k,250)
         mat(k,226) = -rxt(k,293)*y(k,250)
         mat(k,1917) = -(rxt(k,161)*y(k,79) + rxt(k,162)*y(k,81) + rxt(k,163)*y(k,240) &
                      + rxt(k,164)*y(k,138) + rxt(k,165)*y(k,139) + (4._r8*rxt(k,166) &
                      + 4._r8*rxt(k,167)) * y(k,251) + rxt(k,169)*y(k,92) + rxt(k,181) &
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
                      *y(k,144) + rxt(k,323)*y(k,27) + rxt(k,331)*y(k,32) + rxt(k,332) &
                      *y(k,33) + rxt(k,334)*y(k,51) + rxt(k,335)*y(k,97) + rxt(k,336) &
                      *y(k,132) + rxt(k,339)*y(k,153) + rxt(k,343)*y(k,154) + rxt(k,344) &
                      *y(k,31) + rxt(k,345)*y(k,50) + rxt(k,347)*y(k,18) + rxt(k,350) &
                      *y(k,95) + rxt(k,358)*y(k,111) + rxt(k,359)*y(k,112) + rxt(k,368) &
                      *y(k,113) + rxt(k,369)*y(k,114) + rxt(k,370)*y(k,115) + rxt(k,372) &
                      *y(k,118) + rxt(k,375)*y(k,1) + rxt(k,379)*y(k,2) + rxt(k,380) &
                      *y(k,17) + rxt(k,381)*y(k,96) + rxt(k,382)*y(k,98) + rxt(k,383) &
                      *y(k,99) + rxt(k,395)*y(k,101) + rxt(k,396)*y(k,102) + rxt(k,403) &
                      *y(k,104) + rxt(k,405)*y(k,100) + rxt(k,406)*y(k,106) + rxt(k,407) &
                      *y(k,122) + rxt(k,408)*y(k,123) + rxt(k,414)*y(k,219) + rxt(k,417) &
                      *y(k,8) + rxt(k,420)*y(k,10) + rxt(k,421)*y(k,24) + rxt(k,423) &
                      *y(k,25) + rxt(k,427)*y(k,34) + rxt(k,428)*y(k,68) + rxt(k,440) &
                      *y(k,147) + rxt(k,443)*y(k,148) + rxt(k,447)*y(k,216) + rxt(k,448) &
                      *y(k,217) + rxt(k,450)*y(k,220) + rxt(k,453)*y(k,221) + rxt(k,456) &
                      *y(k,223) + rxt(k,457)*y(k,224) + rxt(k,460)*y(k,6) + rxt(k,463) &
                      *y(k,116) + rxt(k,468)*y(k,133) + rxt(k,472)*y(k,211) + rxt(k,473) &
                      *y(k,212) + rxt(k,477)*y(k,213) + rxt(k,479)*y(k,214) + rxt(k,480) &
                      *y(k,215) + (rxt(k,482) + rxt(k,496)) * y(k,69) + rxt(k,484) &
                      *y(k,142) + rxt(k,486)*y(k,158) + rxt(k,490)*y(k,155) + rxt(k,495) &
                      *y(k,157) + rxt(k,498)*y(k,127))
         mat(k,1447) = -rxt(k,161)*y(k,251)
         mat(k,640) = -rxt(k,162)*y(k,251)
         mat(k,2220) = -rxt(k,163)*y(k,251)
         mat(k,1996) = -rxt(k,164)*y(k,251)
         mat(k,2320) = -rxt(k,165)*y(k,251)
         mat(k,447) = -rxt(k,169)*y(k,251)
         mat(k,2077) = -rxt(k,181)*y(k,251)
         mat(k,511) = -rxt(k,182)*y(k,251)
         mat(k,1704) = -rxt(k,190)*y(k,251)
         mat(k,1492) = -rxt(k,191)*y(k,251)
         mat(k,927) = -rxt(k,210)*y(k,251)
         mat(k,1943) = -(rxt(k,212) + rxt(k,213)) * y(k,251)
         mat(k,1965) = -rxt(k,215)*y(k,251)
         mat(k,843) = -rxt(k,218)*y(k,251)
         mat(k,2101) = -rxt(k,242)*y(k,251)
         mat(k,853) = -rxt(k,244)*y(k,251)
         mat(k,520) = -rxt(k,258)*y(k,251)
         mat(k,618) = -rxt(k,260)*y(k,251)
         mat(k,165) = -rxt(k,261)*y(k,251)
         mat(k,410) = -rxt(k,263)*y(k,251)
         mat(k,460) = -rxt(k,265)*y(k,251)
         mat(k,267) = -rxt(k,266)*y(k,251)
         mat(k,288) = -rxt(k,267)*y(k,251)
         mat(k,343) = -rxt(k,268)*y(k,251)
         mat(k,2020) = -rxt(k,277)*y(k,251)
         mat(k,827) = -rxt(k,282)*y(k,251)
         mat(k,454) = -rxt(k,283)*y(k,251)
         mat(k,1475) = -rxt(k,284)*y(k,251)
         mat(k,227) = -rxt(k,285)*y(k,251)
         mat(k,1012) = -rxt(k,286)*y(k,251)
         mat(k,1136) = -rxt(k,294)*y(k,251)
         mat(k,328) = -rxt(k,296)*y(k,251)
         mat(k,305) = -rxt(k,303)*y(k,251)
         mat(k,381) = -rxt(k,304)*y(k,251)
         mat(k,332) = -rxt(k,306)*y(k,251)
         mat(k,1130) = -rxt(k,308)*y(k,251)
         mat(k,143) = -rxt(k,309)*y(k,251)
         mat(k,727) = -rxt(k,314)*y(k,251)
         mat(k,627) = -rxt(k,315)*y(k,251)
         mat(k,1143) = -rxt(k,320)*y(k,251)
         mat(k,972) = -rxt(k,321)*y(k,251)
         mat(k,571) = -rxt(k,322)*y(k,251)
         mat(k,602) = -rxt(k,323)*y(k,251)
         mat(k,430) = -rxt(k,331)*y(k,251)
         mat(k,150) = -rxt(k,332)*y(k,251)
         mat(k,1264) = -rxt(k,334)*y(k,251)
         mat(k,1183) = -rxt(k,335)*y(k,251)
         mat(k,894) = -rxt(k,336)*y(k,251)
         mat(k,594) = -rxt(k,339)*y(k,251)
         mat(k,418) = -rxt(k,343)*y(k,251)
         mat(k,1074) = -rxt(k,344)*y(k,251)
         mat(k,966) = -rxt(k,345)*y(k,251)
         mat(k,404) = -rxt(k,347)*y(k,251)
         mat(k,1200) = -rxt(k,350)*y(k,251)
         mat(k,1254) = -rxt(k,358)*y(k,251)
         mat(k,389) = -rxt(k,359)*y(k,251)
         mat(k,552) = -rxt(k,368)*y(k,251)
         mat(k,356) = -rxt(k,369)*y(k,251)
         mat(k,635) = -rxt(k,370)*y(k,251)
         mat(k,1380) = -rxt(k,372)*y(k,251)
         mat(k,692) = -rxt(k,375)*y(k,251)
         mat(k,713) = -rxt(k,379)*y(k,251)
         mat(k,275) = -rxt(k,380)*y(k,251)
         mat(k,296) = -rxt(k,381)*y(k,251)
         mat(k,385) = -rxt(k,382)*y(k,251)
         mat(k,174) = -rxt(k,383)*y(k,251)
         mat(k,663) = -rxt(k,395)*y(k,251)
         mat(k,588) = -rxt(k,396)*y(k,251)
         mat(k,437) = -rxt(k,403)*y(k,251)
         mat(k,885) = -rxt(k,405)*y(k,251)
         mat(k,734) = -rxt(k,406)*y(k,251)
         mat(k,443) = -rxt(k,407)*y(k,251)
         mat(k,1093) = -rxt(k,408)*y(k,251)
         mat(k,248) = -rxt(k,414)*y(k,251)
         mat(k,211) = -rxt(k,417)*y(k,251)
         mat(k,425) = -rxt(k,420)*y(k,251)
         mat(k,284) = -rxt(k,421)*y(k,251)
         mat(k,376) = -rxt(k,423)*y(k,251)
         mat(k,310) = -rxt(k,427)*y(k,251)
         mat(k,240) = -rxt(k,428)*y(k,251)
         mat(k,220) = -rxt(k,440)*y(k,251)
         mat(k,370) = -rxt(k,443)*y(k,251)
         mat(k,655) = -rxt(k,447)*y(k,251)
         mat(k,235) = -rxt(k,448)*y(k,251)
         mat(k,257) = -rxt(k,450)*y(k,251)
         mat(k,750) = -rxt(k,453)*y(k,251)
         mat(k,262) = -rxt(k,456)*y(k,251)
         mat(k,473) = -rxt(k,457)*y(k,251)
         mat(k,1002) = -rxt(k,460)*y(k,251)
         mat(k,1052) = -rxt(k,463)*y(k,251)
         mat(k,468) = -rxt(k,468)*y(k,251)
         mat(k,702) = -rxt(k,472)*y(k,251)
         mat(k,674) = -rxt(k,473)*y(k,251)
         mat(k,528) = -rxt(k,477)*y(k,251)
         mat(k,1021) = -rxt(k,479)*y(k,251)
         mat(k,1120) = -rxt(k,480)*y(k,251)
         mat(k,321) = -(rxt(k,482) + rxt(k,496)) * y(k,251)
         mat(k,397) = -rxt(k,484)*y(k,251)
         mat(k,899) = -rxt(k,486)*y(k,251)
         mat(k,563) = -rxt(k,490)*y(k,251)
         mat(k,1277) = -rxt(k,495)*y(k,251)
         mat(k,137) = -rxt(k,498)*y(k,251)
         mat(k,1002) = mat(k,1002) + .630_r8*rxt(k,459)*y(k,139)
         mat(k,328) = mat(k,328) + .650_r8*rxt(k,296)*y(k,251)
         mat(k,602) = mat(k,602) + .130_r8*rxt(k,298)*y(k,139)
         mat(k,381) = mat(k,381) + .500_r8*rxt(k,304)*y(k,251)
         mat(k,1074) = mat(k,1074) + .360_r8*rxt(k,327)*y(k,139)
         mat(k,2020) = mat(k,2020) + rxt(k,276)*y(k,138)
         mat(k,454) = mat(k,454) + .300_r8*rxt(k,283)*y(k,251)
         mat(k,1475) = mat(k,1475) + rxt(k,290)*y(k,250)
         mat(k,2259) = rxt(k,199)*y(k,240)
         mat(k,919) = rxt(k,253)*y(k,261)
         mat(k,1507) = rxt(k,160)*y(k,139) + 2.000_r8*rxt(k,155)*y(k,240)
         mat(k,1447) = mat(k,1447) + rxt(k,152)*y(k,138) + rxt(k,144)*y(k,250)
         mat(k,640) = mat(k,640) + rxt(k,153)*y(k,138)
         mat(k,853) = mat(k,853) + rxt(k,243)*y(k,138) + rxt(k,249)*y(k,250)
         mat(k,1965) = mat(k,1965) + rxt(k,214)*y(k,138) + rxt(k,226)*y(k,250)
         mat(k,227) = mat(k,227) + rxt(k,293)*y(k,250)
         mat(k,820) = rxt(k,245)*y(k,138)
         mat(k,843) = mat(k,843) + rxt(k,217)*y(k,138)
         mat(k,885) = mat(k,885) + .320_r8*rxt(k,404)*y(k,139)
         mat(k,734) = mat(k,734) + .600_r8*rxt(k,406)*y(k,251)
         mat(k,1254) = mat(k,1254) + .240_r8*rxt(k,357)*y(k,139)
         mat(k,389) = mat(k,389) + .100_r8*rxt(k,359)*y(k,251)
         mat(k,1052) = mat(k,1052) + .630_r8*rxt(k,462)*y(k,139)
         mat(k,1380) = mat(k,1380) + .360_r8*rxt(k,371)*y(k,139)
         mat(k,1608) = rxt(k,183)*y(k,240)
         mat(k,2077) = mat(k,2077) + rxt(k,178)*y(k,240)
         mat(k,1996) = mat(k,1996) + rxt(k,276)*y(k,44) + rxt(k,152)*y(k,79) &
                      + rxt(k,153)*y(k,81) + rxt(k,243)*y(k,83) + rxt(k,214)*y(k,87) &
                      + rxt(k,245)*y(k,93) + rxt(k,217)*y(k,94) + rxt(k,158)*y(k,240)
         mat(k,2320) = mat(k,2320) + .630_r8*rxt(k,459)*y(k,6) + .130_r8*rxt(k,298) &
                      *y(k,27) + .360_r8*rxt(k,327)*y(k,31) + rxt(k,160)*y(k,78) &
                      + .320_r8*rxt(k,404)*y(k,100) + .240_r8*rxt(k,357)*y(k,111) &
                      + .630_r8*rxt(k,462)*y(k,116) + .360_r8*rxt(k,371)*y(k,118) &
                      + rxt(k,159)*y(k,240)
         mat(k,594) = mat(k,594) + .500_r8*rxt(k,339)*y(k,251)
         mat(k,248) = mat(k,248) + .500_r8*rxt(k,414)*y(k,251)
         mat(k,558) = .400_r8*rxt(k,415)*y(k,240)
         mat(k,1432) = .450_r8*rxt(k,312)*y(k,240)
         mat(k,805) = .400_r8*rxt(k,429)*y(k,240)
         mat(k,2220) = mat(k,2220) + rxt(k,199)*y(k,58) + 2.000_r8*rxt(k,155)*y(k,78) &
                      + rxt(k,183)*y(k,129) + rxt(k,178)*y(k,131) + rxt(k,158) &
                      *y(k,138) + rxt(k,159)*y(k,139) + .400_r8*rxt(k,415)*y(k,227) &
                      + .450_r8*rxt(k,312)*y(k,234) + .400_r8*rxt(k,429)*y(k,236) &
                      + .450_r8*rxt(k,362)*y(k,246) + .400_r8*rxt(k,435)*y(k,247) &
                      + .200_r8*rxt(k,366)*y(k,248) + .150_r8*rxt(k,341)*y(k,254)
         mat(k,1400) = .450_r8*rxt(k,362)*y(k,240)
         mat(k,910) = .400_r8*rxt(k,435)*y(k,240)
         mat(k,722) = .200_r8*rxt(k,366)*y(k,240)
         mat(k,1745) = rxt(k,290)*y(k,56) + rxt(k,144)*y(k,79) + rxt(k,249)*y(k,83) &
                      + rxt(k,226)*y(k,87) + rxt(k,293)*y(k,88) + 2.000_r8*rxt(k,145) &
                      *y(k,261)
         mat(k,1917) = mat(k,1917) + .650_r8*rxt(k,296)*y(k,26) + .500_r8*rxt(k,304) &
                      *y(k,29) + .300_r8*rxt(k,283)*y(k,55) + .600_r8*rxt(k,406) &
                      *y(k,106) + .100_r8*rxt(k,359)*y(k,112) + .500_r8*rxt(k,339) &
                      *y(k,153) + .500_r8*rxt(k,414)*y(k,219)
         mat(k,1176) = .150_r8*rxt(k,341)*y(k,240)
         mat(k,2346) = rxt(k,253)*y(k,75) + 2.000_r8*rxt(k,145)*y(k,250)
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
         mat(k,496) = -(rxt(k,438)*y(k,240) + rxt(k,439)*y(k,129))
         mat(k,2154) = -rxt(k,438)*y(k,252)
         mat(k,1547) = -rxt(k,439)*y(k,252)
         mat(k,238) = .200_r8*rxt(k,428)*y(k,251)
         mat(k,218) = .140_r8*rxt(k,440)*y(k,251)
         mat(k,368) = rxt(k,443)*y(k,251)
         mat(k,1828) = .200_r8*rxt(k,428)*y(k,68) + .140_r8*rxt(k,440)*y(k,147) &
                      + rxt(k,443)*y(k,148)
         mat(k,808) = -(rxt(k,337)*y(k,240) + rxt(k,338)*y(k,129))
         mat(k,2178) = -rxt(k,337)*y(k,253)
         mat(k,1568) = -rxt(k,338)*y(k,253)
         mat(k,1060) = rxt(k,344)*y(k,251)
         mat(k,591) = .500_r8*rxt(k,339)*y(k,251)
         mat(k,1865) = rxt(k,344)*y(k,31) + .500_r8*rxt(k,339)*y(k,153)
         mat(k,1169) = -(rxt(k,340)*y(k,235) + rxt(k,341)*y(k,240) + rxt(k,342) &
                      *y(k,129))
         mat(k,1641) = -rxt(k,340)*y(k,254)
         mat(k,2199) = -rxt(k,341)*y(k,254)
         mat(k,1587) = -rxt(k,342)*y(k,254)
         mat(k,995) = .060_r8*rxt(k,459)*y(k,139)
         mat(k,963) = rxt(k,345)*y(k,251)
         mat(k,1045) = .060_r8*rxt(k,462)*y(k,139)
         mat(k,2300) = .060_r8*rxt(k,459)*y(k,6) + .060_r8*rxt(k,462)*y(k,116)
         mat(k,416) = rxt(k,343)*y(k,251)
         mat(k,1117) = .150_r8*rxt(k,480)*y(k,251)
         mat(k,1894) = rxt(k,345)*y(k,50) + rxt(k,343)*y(k,154) + .150_r8*rxt(k,480) &
                      *y(k,215)
         mat(k,1155) = -(rxt(k,469)*y(k,235) + rxt(k,470)*y(k,240) + rxt(k,471) &
                      *y(k,129))
         mat(k,1640) = -rxt(k,469)*y(k,255)
         mat(k,2198) = -rxt(k,470)*y(k,255)
         mat(k,1586) = -rxt(k,471)*y(k,255)
         mat(k,2054) = .500_r8*rxt(k,478)*y(k,214)
         mat(k,701) = rxt(k,472)*y(k,251)
         mat(k,1020) = .500_r8*rxt(k,478)*y(k,131) + rxt(k,479)*y(k,251)
         mat(k,1893) = rxt(k,472)*y(k,211) + rxt(k,479)*y(k,214)
         mat(k,952) = -(rxt(k,474)*y(k,235) + rxt(k,475)*y(k,240) + rxt(k,476) &
                      *y(k,129))
         mat(k,1630) = -rxt(k,474)*y(k,256)
         mat(k,2188) = -rxt(k,475)*y(k,256)
         mat(k,1575) = -rxt(k,476)*y(k,256)
         mat(k,985) = rxt(k,460)*y(k,251)
         mat(k,1034) = rxt(k,463)*y(k,251)
         mat(k,525) = rxt(k,477)*y(k,251)
         mat(k,1879) = rxt(k,460)*y(k,6) + rxt(k,463)*y(k,116) + rxt(k,477)*y(k,213)
         mat(k,772) = -(rxt(k,445)*y(k,240) + rxt(k,446)*y(k,129))
         mat(k,2175) = -rxt(k,445)*y(k,257)
         mat(k,1565) = -rxt(k,446)*y(k,257)
         mat(k,651) = rxt(k,447)*y(k,251)
         mat(k,234) = .650_r8*rxt(k,448)*y(k,251)
         mat(k,1862) = rxt(k,447)*y(k,216) + .650_r8*rxt(k,448)*y(k,217)
         mat(k,1213) = -(rxt(k,409)*y(k,234) + rxt(k,410)*y(k,235) + rxt(k,411) &
                      *y(k,240) + rxt(k,412)*y(k,129) + rxt(k,413)*y(k,131))
         mat(k,1417) = -rxt(k,409)*y(k,258)
         mat(k,1643) = -rxt(k,410)*y(k,258)
         mat(k,2201) = -rxt(k,411)*y(k,258)
         mat(k,1590) = -rxt(k,412)*y(k,258)
         mat(k,2058) = -rxt(k,413)*y(k,258)
         mat(k,295) = rxt(k,381)*y(k,251)
         mat(k,384) = rxt(k,382)*y(k,251)
         mat(k,173) = rxt(k,383)*y(k,251)
         mat(k,731) = .400_r8*rxt(k,406)*y(k,251)
         mat(k,247) = .500_r8*rxt(k,414)*y(k,251)
         mat(k,1897) = rxt(k,381)*y(k,96) + rxt(k,382)*y(k,98) + rxt(k,383)*y(k,99) &
                      + .400_r8*rxt(k,406)*y(k,106) + .500_r8*rxt(k,414)*y(k,219)
         mat(k,788) = -(rxt(k,451)*y(k,240) + rxt(k,452)*y(k,129))
         mat(k,2176) = -rxt(k,451)*y(k,259)
         mat(k,1566) = -rxt(k,452)*y(k,259)
         mat(k,254) = .560_r8*rxt(k,450)*y(k,251)
         mat(k,743) = rxt(k,453)*y(k,251)
         mat(k,1863) = .560_r8*rxt(k,450)*y(k,220) + rxt(k,453)*y(k,221)
         mat(k,541) = -(rxt(k,454)*y(k,240) + rxt(k,455)*y(k,129))
         mat(k,2159) = -rxt(k,454)*y(k,260)
         mat(k,1552) = -rxt(k,455)*y(k,260)
         mat(k,261) = .300_r8*rxt(k,456)*y(k,251)
         mat(k,470) = rxt(k,457)*y(k,251)
         mat(k,1835) = .300_r8*rxt(k,456)*y(k,223) + rxt(k,457)*y(k,224)
         mat(k,2356) = -(rxt(k,145)*y(k,250) + rxt(k,253)*y(k,75) + rxt(k,497) &
                      *y(k,159))
         mat(k,1755) = -rxt(k,145)*y(k,261)
         mat(k,921) = -rxt(k,253)*y(k,261)
         mat(k,302) = -rxt(k,497)*y(k,261)
         mat(k,335) = rxt(k,306)*y(k,251)
         mat(k,432) = rxt(k,331)*y(k,251)
         mat(k,151) = rxt(k,332)*y(k,251)
         mat(k,523) = rxt(k,258)*y(k,251)
         mat(k,2030) = rxt(k,277)*y(k,251)
         mat(k,622) = rxt(k,260)*y(k,251)
         mat(k,167) = rxt(k,261)*y(k,251)
         mat(k,1133) = rxt(k,308)*y(k,251)
         mat(k,414) = rxt(k,263)*y(k,251)
         mat(k,967) = rxt(k,345)*y(k,251)
         mat(k,1268) = rxt(k,334)*y(k,251)
         mat(k,728) = rxt(k,314)*y(k,251)
         mat(k,629) = rxt(k,315)*y(k,251)
         mat(k,456) = rxt(k,283)*y(k,251)
         mat(k,1482) = rxt(k,284)*y(k,251)
         mat(k,1515) = rxt(k,156)*y(k,240)
         mat(k,1452) = rxt(k,161)*y(k,251)
         mat(k,645) = rxt(k,162)*y(k,251)
         mat(k,856) = rxt(k,244)*y(k,251)
         mat(k,345) = rxt(k,268)*y(k,251)
         mat(k,1975) = (rxt(k,553)+rxt(k,558))*y(k,93) + (rxt(k,546)+rxt(k,552) &
                       +rxt(k,557))*y(k,94) + rxt(k,215)*y(k,251)
         mat(k,1015) = rxt(k,286)*y(k,251)
         mat(k,1499) = rxt(k,191)*y(k,251)
         mat(k,450) = rxt(k,169)*y(k,251)
         mat(k,825) = (rxt(k,553)+rxt(k,558))*y(k,87)
         mat(k,848) = (rxt(k,546)+rxt(k,552)+rxt(k,557))*y(k,87) + rxt(k,218)*y(k,251)
         mat(k,1259) = .500_r8*rxt(k,358)*y(k,251)
         mat(k,138) = rxt(k,498)*y(k,251)
         mat(k,597) = rxt(k,339)*y(k,251)
         mat(k,420) = rxt(k,343)*y(k,251)
         mat(k,2230) = rxt(k,156)*y(k,78) + rxt(k,163)*y(k,251)
         mat(k,1927) = rxt(k,306)*y(k,30) + rxt(k,331)*y(k,32) + rxt(k,332)*y(k,33) &
                      + rxt(k,258)*y(k,43) + rxt(k,277)*y(k,44) + rxt(k,260)*y(k,45) &
                      + rxt(k,261)*y(k,46) + rxt(k,308)*y(k,47) + rxt(k,263)*y(k,48) &
                      + rxt(k,345)*y(k,50) + rxt(k,334)*y(k,51) + rxt(k,314)*y(k,52) &
                      + rxt(k,315)*y(k,53) + rxt(k,283)*y(k,55) + rxt(k,284)*y(k,56) &
                      + rxt(k,161)*y(k,79) + rxt(k,162)*y(k,81) + rxt(k,244)*y(k,83) &
                      + rxt(k,268)*y(k,86) + rxt(k,215)*y(k,87) + rxt(k,286)*y(k,89) &
                      + rxt(k,191)*y(k,91) + rxt(k,169)*y(k,92) + rxt(k,218)*y(k,94) &
                      + .500_r8*rxt(k,358)*y(k,111) + rxt(k,498)*y(k,127) + rxt(k,339) &
                      *y(k,153) + rxt(k,343)*y(k,154) + rxt(k,163)*y(k,240) &
                      + 2.000_r8*rxt(k,166)*y(k,251)
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
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 70) = lmat(k, 70)
         mat(k, 71) = lmat(k, 71)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 92) = mat(k, 92) + lmat(k, 92)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 129) = mat(k, 129) + lmat(k, 129)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 145) = lmat(k, 145)
         mat(k, 146) = lmat(k, 146)
         mat(k, 147) = lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
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
         mat(k, 166) = mat(k, 166) + lmat(k, 166)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 175) = lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
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
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 230) = mat(k, 230) + lmat(k, 230)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 242) = lmat(k, 242)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 273) = mat(k, 273) + lmat(k, 273)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 289) = mat(k, 289) + lmat(k, 289)
         mat(k, 290) = mat(k, 290) + lmat(k, 290)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 293) = mat(k, 293) + lmat(k, 293)
         mat(k, 294) = lmat(k, 294)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 297) = lmat(k, 297)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 300) = lmat(k, 300)
         mat(k, 301) = lmat(k, 301)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 308) = lmat(k, 308)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 312) = lmat(k, 312)
         mat(k, 313) = lmat(k, 313)
         mat(k, 314) = lmat(k, 314)
         mat(k, 315) = lmat(k, 315)
         mat(k, 316) = lmat(k, 316)
         mat(k, 317) = lmat(k, 317)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 336) = lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 346) = lmat(k, 346)
         mat(k, 347) = lmat(k, 347)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = lmat(k, 349)
         mat(k, 350) = lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = lmat(k, 355)
         mat(k, 356) = mat(k, 356) + lmat(k, 356)
         mat(k, 357) = lmat(k, 357)
         mat(k, 358) = lmat(k, 358)
         mat(k, 359) = lmat(k, 359)
         mat(k, 360) = lmat(k, 360)
         mat(k, 361) = lmat(k, 361)
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
         mat(k, 383) = mat(k, 383) + lmat(k, 383)
         mat(k, 386) = mat(k, 386) + lmat(k, 386)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 409) = lmat(k, 409)
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
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 431) = lmat(k, 431)
         mat(k, 433) = mat(k, 433) + lmat(k, 433)
         mat(k, 434) = lmat(k, 434)
         mat(k, 436) = lmat(k, 436)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 442) = lmat(k, 442)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = lmat(k, 448)
         mat(k, 449) = lmat(k, 449)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 452) = lmat(k, 452)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 458) = mat(k, 458) + lmat(k, 458)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 464) = lmat(k, 464)
         mat(k, 466) = lmat(k, 466)
         mat(k, 467) = lmat(k, 467)
         mat(k, 468) = mat(k, 468) + lmat(k, 468)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 471) = lmat(k, 471)
         mat(k, 472) = lmat(k, 472)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 474) = lmat(k, 474)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 483) = lmat(k, 483)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = lmat(k, 486)
         mat(k, 487) = lmat(k, 487)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = mat(k, 489) + lmat(k, 489)
         mat(k, 493) = lmat(k, 493)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 502) = lmat(k, 502)
         mat(k, 503) = lmat(k, 503)
         mat(k, 504) = lmat(k, 504)
         mat(k, 506) = mat(k, 506) + lmat(k, 506)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 526) = lmat(k, 526)
         mat(k, 527) = lmat(k, 527)
         mat(k, 528) = mat(k, 528) + lmat(k, 528)
         mat(k, 529) = lmat(k, 529)
         mat(k, 530) = lmat(k, 530)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 550) = lmat(k, 550)
         mat(k, 551) = lmat(k, 551)
         mat(k, 554) = mat(k, 554) + lmat(k, 554)
         mat(k, 560) = mat(k, 560) + lmat(k, 560)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 564) = lmat(k, 564)
         mat(k, 566) = mat(k, 566) + lmat(k, 566)
         mat(k, 567) = lmat(k, 567)
         mat(k, 568) = lmat(k, 568)
         mat(k, 569) = lmat(k, 569)
         mat(k, 570) = lmat(k, 570)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 575) = lmat(k, 575)
         mat(k, 576) = lmat(k, 576)
         mat(k, 577) = lmat(k, 577)
         mat(k, 578) = lmat(k, 578)
         mat(k, 580) = mat(k, 580) + lmat(k, 580)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 587) = lmat(k, 587)
         mat(k, 590) = mat(k, 590) + lmat(k, 590)
         mat(k, 592) = lmat(k, 592)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 595) = lmat(k, 595)
         mat(k, 596) = lmat(k, 596)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 614) = mat(k, 614) + lmat(k, 614)
         mat(k, 615) = mat(k, 615) + lmat(k, 615)
         mat(k, 616) = lmat(k, 616)
         mat(k, 623) = mat(k, 623) + lmat(k, 623)
         mat(k, 624) = mat(k, 624) + lmat(k, 624)
         mat(k, 626) = lmat(k, 626)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 630) = mat(k, 630) + lmat(k, 630)
         mat(k, 633) = lmat(k, 633)
         mat(k, 634) = lmat(k, 634)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 646) = lmat(k, 646)
         mat(k, 647) = lmat(k, 647)
         mat(k, 648) = mat(k, 648) + lmat(k, 648)
         mat(k, 649) = lmat(k, 649)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 656) = lmat(k, 656)
         mat(k, 657) = mat(k, 657) + lmat(k, 657)
         mat(k, 661) = lmat(k, 661)
         mat(k, 666) = lmat(k, 666)
         mat(k, 667) = lmat(k, 667)
         mat(k, 668) = lmat(k, 668)
         mat(k, 669) = lmat(k, 669)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 672) = lmat(k, 672)
         mat(k, 673) = mat(k, 673) + lmat(k, 673)
         mat(k, 675) = lmat(k, 675)
         mat(k, 678) = mat(k, 678) + lmat(k, 678)
         mat(k, 684) = lmat(k, 684)
         mat(k, 685) = mat(k, 685) + lmat(k, 685)
         mat(k, 688) = mat(k, 688) + lmat(k, 688)
         mat(k, 689) = mat(k, 689) + lmat(k, 689)
         mat(k, 691) = mat(k, 691) + lmat(k, 691)
         mat(k, 693) = mat(k, 693) + lmat(k, 693)
         mat(k, 694) = lmat(k, 694)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 696) = lmat(k, 696)
         mat(k, 697) = lmat(k, 697)
         mat(k, 698) = lmat(k, 698)
         mat(k, 699) = lmat(k, 699)
         mat(k, 700) = lmat(k, 700)
         mat(k, 702) = mat(k, 702) + lmat(k, 702)
         mat(k, 703) = lmat(k, 703)
         mat(k, 704) = lmat(k, 704)
         mat(k, 705) = lmat(k, 705)
         mat(k, 706) = mat(k, 706) + lmat(k, 706)
         mat(k, 710) = lmat(k, 710)
         mat(k, 711) = lmat(k, 711)
         mat(k, 713) = mat(k, 713) + lmat(k, 713)
         mat(k, 714) = lmat(k, 714)
         mat(k, 715) = lmat(k, 715)
         mat(k, 717) = mat(k, 717) + lmat(k, 717)
         mat(k, 724) = mat(k, 724) + lmat(k, 724)
         mat(k, 730) = mat(k, 730) + lmat(k, 730)
         mat(k, 732) = lmat(k, 732)
         mat(k, 733) = lmat(k, 733)
         mat(k, 734) = mat(k, 734) + lmat(k, 734)
         mat(k, 735) = lmat(k, 735)
         mat(k, 736) = lmat(k, 736)
         mat(k, 737) = lmat(k, 737)
         mat(k, 738) = lmat(k, 738)
         mat(k, 739) = lmat(k, 739)
         mat(k, 740) = lmat(k, 740)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 746) = lmat(k, 746)
         mat(k, 748) = lmat(k, 748)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 751) = lmat(k, 751)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 761) = mat(k, 761) + lmat(k, 761)
         mat(k, 772) = mat(k, 772) + lmat(k, 772)
         mat(k, 788) = mat(k, 788) + lmat(k, 788)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 808) = mat(k, 808) + lmat(k, 808)
         mat(k, 818) = mat(k, 818) + lmat(k, 818)
         mat(k, 819) = lmat(k, 819)
         mat(k, 820) = mat(k, 820) + lmat(k, 820)
         mat(k, 826) = mat(k, 826) + lmat(k, 826)
         mat(k, 831) = mat(k, 831) + lmat(k, 831)
         mat(k, 842) = mat(k, 842) + lmat(k, 842)
         mat(k, 843) = mat(k, 843) + lmat(k, 843)
         mat(k, 847) = mat(k, 847) + lmat(k, 847)
         mat(k, 849) = mat(k, 849) + lmat(k, 849)
         mat(k, 850) = mat(k, 850) + lmat(k, 850)
         mat(k, 851) = mat(k, 851) + lmat(k, 851)
         mat(k, 860) = mat(k, 860) + lmat(k, 860)
         mat(k, 868) = lmat(k, 868)
         mat(k, 869) = lmat(k, 869)
         mat(k, 870) = lmat(k, 870)
         mat(k, 874) = mat(k, 874) + lmat(k, 874)
         mat(k, 890) = mat(k, 890) + lmat(k, 890)
         mat(k, 892) = lmat(k, 892)
         mat(k, 893) = mat(k, 893) + lmat(k, 893)
         mat(k, 895) = lmat(k, 895)
         mat(k, 897) = mat(k, 897) + lmat(k, 897)
         mat(k, 898) = lmat(k, 898)
         mat(k, 900) = lmat(k, 900)
         mat(k, 903) = mat(k, 903) + lmat(k, 903)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 923) = mat(k, 923) + lmat(k, 923)
         mat(k, 924) = mat(k, 924) + lmat(k, 924)
         mat(k, 925) = mat(k, 925) + lmat(k, 925)
         mat(k, 926) = lmat(k, 926)
         mat(k, 928) = mat(k, 928) + lmat(k, 928)
         mat(k, 931) = mat(k, 931) + lmat(k, 931)
         mat(k, 932) = mat(k, 932) + lmat(k, 932)
         mat(k, 939) = mat(k, 939) + lmat(k, 939)
         mat(k, 952) = mat(k, 952) + lmat(k, 952)
         mat(k, 962) = mat(k, 962) + lmat(k, 962)
         mat(k, 964) = lmat(k, 964)
         mat(k, 965) = lmat(k, 965)
         mat(k, 968) = lmat(k, 968)
         mat(k, 970) = mat(k, 970) + lmat(k, 970)
         mat(k, 971) = mat(k, 971) + lmat(k, 971)
         mat(k, 974) = mat(k, 974) + lmat(k, 974)
         mat(k, 988) = mat(k, 988) + lmat(k, 988)
         mat(k,1010) = mat(k,1010) + lmat(k,1010)
         mat(k,1016) = mat(k,1016) + lmat(k,1016)
         mat(k,1018) = lmat(k,1018)
         mat(k,1019) = lmat(k,1019)
         mat(k,1023) = lmat(k,1023)
         mat(k,1039) = mat(k,1039) + lmat(k,1039)
         mat(k,1064) = mat(k,1064) + lmat(k,1064)
         mat(k,1081) = lmat(k,1081)
         mat(k,1085) = mat(k,1085) + lmat(k,1085)
         mat(k,1089) = lmat(k,1089)
         mat(k,1092) = lmat(k,1092)
         mat(k,1096) = mat(k,1096) + lmat(k,1096)
         mat(k,1101) = mat(k,1101) + lmat(k,1101)
         mat(k,1113) = mat(k,1113) + lmat(k,1113)
         mat(k,1114) = mat(k,1114) + lmat(k,1114)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1116) = mat(k,1116) + lmat(k,1116)
         mat(k,1117) = mat(k,1117) + lmat(k,1117)
         mat(k,1118) = mat(k,1118) + lmat(k,1118)
         mat(k,1122) = mat(k,1122) + lmat(k,1122)
         mat(k,1123) = mat(k,1123) + lmat(k,1123)
         mat(k,1125) = mat(k,1125) + lmat(k,1125)
         mat(k,1126) = lmat(k,1126)
         mat(k,1129) = lmat(k,1129)
         mat(k,1132) = lmat(k,1132)
         mat(k,1135) = mat(k,1135) + lmat(k,1135)
         mat(k,1141) = lmat(k,1141)
         mat(k,1142) = mat(k,1142) + lmat(k,1142)
         mat(k,1145) = mat(k,1145) + lmat(k,1145)
         mat(k,1146) = mat(k,1146) + lmat(k,1146)
         mat(k,1155) = mat(k,1155) + lmat(k,1155)
         mat(k,1169) = mat(k,1169) + lmat(k,1169)
         mat(k,1180) = mat(k,1180) + lmat(k,1180)
         mat(k,1182) = lmat(k,1182)
         mat(k,1184) = lmat(k,1184)
         mat(k,1185) = mat(k,1185) + lmat(k,1185)
         mat(k,1187) = lmat(k,1187)
         mat(k,1188) = lmat(k,1188)
         mat(k,1189) = lmat(k,1189)
         mat(k,1190) = lmat(k,1190)
         mat(k,1192) = lmat(k,1192)
         mat(k,1193) = mat(k,1193) + lmat(k,1193)
         mat(k,1195) = lmat(k,1195)
         mat(k,1196) = lmat(k,1196)
         mat(k,1199) = lmat(k,1199)
         mat(k,1202) = lmat(k,1202)
         mat(k,1204) = mat(k,1204) + lmat(k,1204)
         mat(k,1213) = mat(k,1213) + lmat(k,1213)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1248) = mat(k,1248) + lmat(k,1248)
         mat(k,1249) = mat(k,1249) + lmat(k,1249)
         mat(k,1252) = mat(k,1252) + lmat(k,1252)
         mat(k,1253) = mat(k,1253) + lmat(k,1253)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1257) = mat(k,1257) + lmat(k,1257)
         mat(k,1260) = mat(k,1260) + lmat(k,1260)
         mat(k,1261) = mat(k,1261) + lmat(k,1261)
         mat(k,1262) = mat(k,1262) + lmat(k,1262)
         mat(k,1267) = lmat(k,1267)
         mat(k,1270) = lmat(k,1270)
         mat(k,1271) = mat(k,1271) + lmat(k,1271)
         mat(k,1272) = mat(k,1272) + lmat(k,1272)
         mat(k,1279) = lmat(k,1279)
         mat(k,1296) = mat(k,1296) + lmat(k,1296)
         mat(k,1312) = lmat(k,1312)
         mat(k,1328) = mat(k,1328) + lmat(k,1328)
         mat(k,1341) = mat(k,1341) + lmat(k,1341)
         mat(k,1351) = mat(k,1351) + lmat(k,1351)
         mat(k,1366) = lmat(k,1366)
         mat(k,1368) = mat(k,1368) + lmat(k,1368)
         mat(k,1372) = mat(k,1372) + lmat(k,1372)
         mat(k,1374) = mat(k,1374) + lmat(k,1374)
         mat(k,1378) = lmat(k,1378)
         mat(k,1394) = mat(k,1394) + lmat(k,1394)
         mat(k,1426) = mat(k,1426) + lmat(k,1426)
         mat(k,1441) = mat(k,1441) + lmat(k,1441)
         mat(k,1455) = mat(k,1455) + lmat(k,1455)
         mat(k,1466) = lmat(k,1466)
         mat(k,1468) = lmat(k,1468)
         mat(k,1469) = mat(k,1469) + lmat(k,1469)
         mat(k,1470) = mat(k,1470) + lmat(k,1470)
         mat(k,1472) = mat(k,1472) + lmat(k,1472)
         mat(k,1473) = mat(k,1473) + lmat(k,1473)
         mat(k,1475) = mat(k,1475) + lmat(k,1475)
         mat(k,1477) = lmat(k,1477)
         mat(k,1478) = mat(k,1478) + lmat(k,1478)
         mat(k,1482) = mat(k,1482) + lmat(k,1482)
         mat(k,1487) = mat(k,1487) + lmat(k,1487)
         mat(k,1490) = lmat(k,1490)
         mat(k,1492) = mat(k,1492) + lmat(k,1492)
         mat(k,1503) = mat(k,1503) + lmat(k,1503)
         mat(k,1512) = mat(k,1512) + lmat(k,1512)
         mat(k,1549) = mat(k,1549) + lmat(k,1549)
         mat(k,1604) = mat(k,1604) + lmat(k,1604)
         mat(k,1611) = mat(k,1611) + lmat(k,1611)
         mat(k,1657) = mat(k,1657) + lmat(k,1657)
         mat(k,1698) = mat(k,1698) + lmat(k,1698)
         mat(k,1700) = mat(k,1700) + lmat(k,1700)
         mat(k,1702) = mat(k,1702) + lmat(k,1702)
         mat(k,1704) = mat(k,1704) + lmat(k,1704)
         mat(k,1707) = mat(k,1707) + lmat(k,1707)
         mat(k,1744) = mat(k,1744) + lmat(k,1744)
         mat(k,1748) = lmat(k,1748)
         mat(k,1917) = mat(k,1917) + lmat(k,1917)
         mat(k,1944) = mat(k,1944) + lmat(k,1944)
         mat(k,1946) = mat(k,1946) + lmat(k,1946)
         mat(k,1951) = mat(k,1951) + lmat(k,1951)
         mat(k,1961) = mat(k,1961) + lmat(k,1961)
         mat(k,1967) = mat(k,1967) + lmat(k,1967)
         mat(k,1973) = mat(k,1973) + lmat(k,1973)
         mat(k,1999) = mat(k,1999) + lmat(k,1999)
         mat(k,2005) = mat(k,2005) + lmat(k,2005)
         mat(k,2010) = mat(k,2010) + lmat(k,2010)
         mat(k,2011) = lmat(k,2011)
         mat(k,2015) = mat(k,2015) + lmat(k,2015)
         mat(k,2024) = mat(k,2024) + lmat(k,2024)
         mat(k,2071) = mat(k,2071) + lmat(k,2071)
         mat(k,2073) = mat(k,2073) + lmat(k,2073)
         mat(k,2075) = mat(k,2075) + lmat(k,2075)
         mat(k,2080) = mat(k,2080) + lmat(k,2080)
         mat(k,2082) = mat(k,2082) + lmat(k,2082)
         mat(k,2094) = mat(k,2094) + lmat(k,2094)
         mat(k,2104) = mat(k,2104) + lmat(k,2104)
         mat(k,2107) = mat(k,2107) + lmat(k,2107)
         mat(k,2163) = mat(k,2163) + lmat(k,2163)
         mat(k,2227) = mat(k,2227) + lmat(k,2227)
         mat(k,2267) = mat(k,2267) + lmat(k,2267)
         mat(k,2319) = mat(k,2319) + lmat(k,2319)
         mat(k,2323) = mat(k,2323) + lmat(k,2323)
         mat(k,2329) = mat(k,2329) + lmat(k,2329)
         mat(k,2337) = lmat(k,2337)
         mat(k,2341) = lmat(k,2341)
         mat(k,2345) = mat(k,2345) + lmat(k,2345)
         mat(k,2346) = mat(k,2346) + lmat(k,2346)
         mat(k,2349) = lmat(k,2349)
         mat(k,2356) = mat(k,2356) + lmat(k,2356)
         mat(k, 255) = 0._r8
         mat(k, 256) = 0._r8
         mat(k, 341) = 0._r8
         mat(k, 375) = 0._r8
         mat(k, 478) = 0._r8
         mat(k, 481) = 0._r8
         mat(k, 500) = 0._r8
         mat(k, 510) = 0._r8
         mat(k, 534) = 0._r8
         mat(k, 538) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 650) = 0._r8
         mat(k, 652) = 0._r8
         mat(k, 681) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 686) = 0._r8
         mat(k, 687) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 707) = 0._r8
         mat(k, 708) = 0._r8
         mat(k, 712) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 744) = 0._r8
         mat(k, 745) = 0._r8
         mat(k, 747) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 773) = 0._r8
         mat(k, 774) = 0._r8
         mat(k, 776) = 0._r8
         mat(k, 780) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 824) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 955) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 973) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 987) = 0._r8
         mat(k, 993) = 0._r8
         mat(k, 996) = 0._r8
         mat(k, 998) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1000) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1003) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1011) = 0._r8
         mat(k,1013) = 0._r8
         mat(k,1033) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1042) = 0._r8
         mat(k,1044) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1069) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1083) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1087) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1124) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1163) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1219) = 0._r8
         mat(k,1220) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1230) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1318) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1329) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1331) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1363) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1375) = 0._r8
         mat(k,1376) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1442) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1459) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1471) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1485) = 0._r8
         mat(k,1486) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1493) = 0._r8
         mat(k,1494) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1502) = 0._r8
         mat(k,1504) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1510) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1513) = 0._r8
         mat(k,1572) = 0._r8
         mat(k,1601) = 0._r8
         mat(k,1602) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1607) = 0._r8
         mat(k,1610) = 0._r8
         mat(k,1618) = 0._r8
         mat(k,1629) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1654) = 0._r8
         mat(k,1655) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1662) = 0._r8
         mat(k,1663) = 0._r8
         mat(k,1665) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1670) = 0._r8
         mat(k,1681) = 0._r8
         mat(k,1684) = 0._r8
         mat(k,1685) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1690) = 0._r8
         mat(k,1691) = 0._r8
         mat(k,1692) = 0._r8
         mat(k,1696) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1699) = 0._r8
         mat(k,1701) = 0._r8
         mat(k,1703) = 0._r8
         mat(k,1706) = 0._r8
         mat(k,1708) = 0._r8
         mat(k,1712) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1739) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1750) = 0._r8
         mat(k,1829) = 0._r8
         mat(k,1850) = 0._r8
         mat(k,1861) = 0._r8
         mat(k,1864) = 0._r8
         mat(k,1875) = 0._r8
         mat(k,1876) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1937) = 0._r8
         mat(k,1938) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1953) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1963) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1971) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1977) = 0._r8
         mat(k,1983) = 0._r8
         mat(k,1989) = 0._r8
         mat(k,1990) = 0._r8
         mat(k,1993) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2009) = 0._r8
         mat(k,2013) = 0._r8
         mat(k,2016) = 0._r8
         mat(k,2017) = 0._r8
         mat(k,2018) = 0._r8
         mat(k,2019) = 0._r8
         mat(k,2021) = 0._r8
         mat(k,2026) = 0._r8
         mat(k,2029) = 0._r8
         mat(k,2037) = 0._r8
         mat(k,2044) = 0._r8
         mat(k,2050) = 0._r8
         mat(k,2055) = 0._r8
         mat(k,2057) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2069) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2079) = 0._r8
         mat(k,2083) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2095) = 0._r8
         mat(k,2096) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2100) = 0._r8
         mat(k,2103) = 0._r8
         mat(k,2105) = 0._r8
         mat(k,2106) = 0._r8
         mat(k,2110) = 0._r8
         mat(k,2111) = 0._r8
         mat(k,2138) = 0._r8
         mat(k,2155) = 0._r8
         mat(k,2156) = 0._r8
         mat(k,2165) = 0._r8
         mat(k,2185) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2192) = 0._r8
         mat(k,2194) = 0._r8
         mat(k,2197) = 0._r8
         mat(k,2200) = 0._r8
         mat(k,2203) = 0._r8
         mat(k,2208) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2214) = 0._r8
         mat(k,2219) = 0._r8
         mat(k,2240) = 0._r8
         mat(k,2244) = 0._r8
         mat(k,2245) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2248) = 0._r8
         mat(k,2249) = 0._r8
         mat(k,2253) = 0._r8
         mat(k,2255) = 0._r8
         mat(k,2257) = 0._r8
         mat(k,2258) = 0._r8
         mat(k,2262) = 0._r8
         mat(k,2265) = 0._r8
         mat(k,2269) = 0._r8
         mat(k,2282) = 0._r8
         mat(k,2286) = 0._r8
         mat(k,2288) = 0._r8
         mat(k,2294) = 0._r8
         mat(k,2298) = 0._r8
         mat(k,2299) = 0._r8
         mat(k,2301) = 0._r8
         mat(k,2302) = 0._r8
         mat(k,2306) = 0._r8
         mat(k,2307) = 0._r8
         mat(k,2308) = 0._r8
         mat(k,2310) = 0._r8
         mat(k,2314) = 0._r8
         mat(k,2322) = 0._r8
         mat(k,2330) = 0._r8
         mat(k,2334) = 0._r8
         mat(k,2336) = 0._r8
         mat(k,2338) = 0._r8
         mat(k,2339) = 0._r8
         mat(k,2340) = 0._r8
         mat(k,2342) = 0._r8
         mat(k,2343) = 0._r8
         mat(k,2344) = 0._r8
         mat(k,2347) = 0._r8
         mat(k,2348) = 0._r8
         mat(k,2350) = 0._r8
         mat(k,2351) = 0._r8
         mat(k,2352) = 0._r8
         mat(k,2353) = 0._r8
         mat(k,2354) = 0._r8
         mat(k,2355) = 0._r8
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
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 70) = mat(k, 70) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 72) = mat(k, 72) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 80) = mat(k, 80) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 94) = mat(k, 94) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 117) = mat(k, 117) - dti(k)
         mat(k, 123) = mat(k, 123) - dti(k)
         mat(k, 129) = mat(k, 129) - dti(k)
         mat(k, 135) = mat(k, 135) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 145) = mat(k, 145) - dti(k)
         mat(k, 148) = mat(k, 148) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 156) = mat(k, 156) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 172) = mat(k, 172) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 188) = mat(k, 188) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 198) = mat(k, 198) - dti(k)
         mat(k, 201) = mat(k, 201) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 217) = mat(k, 217) - dti(k)
         mat(k, 222) = mat(k, 222) - dti(k)
         mat(k, 225) = mat(k, 225) - dti(k)
         mat(k, 230) = mat(k, 230) - dti(k)
         mat(k, 237) = mat(k, 237) - dti(k)
         mat(k, 242) = mat(k, 242) - dti(k)
         mat(k, 246) = mat(k, 246) - dti(k)
         mat(k, 251) = mat(k, 251) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 264) = mat(k, 264) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 290) = mat(k, 290) - dti(k)
         mat(k, 293) = mat(k, 293) - dti(k)
         mat(k, 299) = mat(k, 299) - dti(k)
         mat(k, 303) = mat(k, 303) - dti(k)
         mat(k, 307) = mat(k, 307) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 315) = mat(k, 315) - dti(k)
         mat(k, 318) = mat(k, 318) - dti(k)
         mat(k, 324) = mat(k, 324) - dti(k)
         mat(k, 330) = mat(k, 330) - dti(k)
         mat(k, 336) = mat(k, 336) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 346) = mat(k, 346) - dti(k)
         mat(k, 352) = mat(k, 352) - dti(k)
         mat(k, 357) = mat(k, 357) - dti(k)
         mat(k, 362) = mat(k, 362) - dti(k)
         mat(k, 367) = mat(k, 367) - dti(k)
         mat(k, 373) = mat(k, 373) - dti(k)
         mat(k, 378) = mat(k, 378) - dti(k)
         mat(k, 383) = mat(k, 383) - dti(k)
         mat(k, 386) = mat(k, 386) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 399) = mat(k, 399) - dti(k)
         mat(k, 407) = mat(k, 407) - dti(k)
         mat(k, 415) = mat(k, 415) - dti(k)
         mat(k, 421) = mat(k, 421) - dti(k)
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
         mat(k, 496) = mat(k, 496) - dti(k)
         mat(k, 502) = mat(k, 502) - dti(k)
         mat(k, 506) = mat(k, 506) - dti(k)
         mat(k, 513) = mat(k, 513) - dti(k)
         mat(k, 517) = mat(k, 517) - dti(k)
         mat(k, 524) = mat(k, 524) - dti(k)
         mat(k, 533) = mat(k, 533) - dti(k)
         mat(k, 541) = mat(k, 541) - dti(k)
         mat(k, 548) = mat(k, 548) - dti(k)
         mat(k, 554) = mat(k, 554) - dti(k)
         mat(k, 560) = mat(k, 560) - dti(k)
         mat(k, 566) = mat(k, 566) - dti(k)
         mat(k, 574) = mat(k, 574) - dti(k)
         mat(k, 582) = mat(k, 582) - dti(k)
         mat(k, 590) = mat(k, 590) - dti(k)
         mat(k, 598) = mat(k, 598) - dti(k)
         mat(k, 606) = mat(k, 606) - dti(k)
         mat(k, 614) = mat(k, 614) - dti(k)
         mat(k, 623) = mat(k, 623) - dti(k)
         mat(k, 630) = mat(k, 630) - dti(k)
         mat(k, 639) = mat(k, 639) - dti(k)
         mat(k, 648) = mat(k, 648) - dti(k)
         mat(k, 657) = mat(k, 657) - dti(k)
         mat(k, 666) = mat(k, 666) - dti(k)
         mat(k, 670) = mat(k, 670) - dti(k)
         mat(k, 678) = mat(k, 678) - dti(k)
         mat(k, 685) = mat(k, 685) - dti(k)
         mat(k, 695) = mat(k, 695) - dti(k)
         mat(k, 706) = mat(k, 706) - dti(k)
         mat(k, 717) = mat(k, 717) - dti(k)
         mat(k, 724) = mat(k, 724) - dti(k)
         mat(k, 730) = mat(k, 730) - dti(k)
         mat(k, 741) = mat(k, 741) - dti(k)
         mat(k, 754) = mat(k, 754) - dti(k)
         mat(k, 761) = mat(k, 761) - dti(k)
         mat(k, 772) = mat(k, 772) - dti(k)
         mat(k, 788) = mat(k, 788) - dti(k)
         mat(k, 799) = mat(k, 799) - dti(k)
         mat(k, 808) = mat(k, 808) - dti(k)
         mat(k, 818) = mat(k, 818) - dti(k)
         mat(k, 826) = mat(k, 826) - dti(k)
         mat(k, 831) = mat(k, 831) - dti(k)
         mat(k, 842) = mat(k, 842) - dti(k)
         mat(k, 849) = mat(k, 849) - dti(k)
         mat(k, 860) = mat(k, 860) - dti(k)
         mat(k, 868) = mat(k, 868) - dti(k)
         mat(k, 874) = mat(k, 874) - dti(k)
         mat(k, 890) = mat(k, 890) - dti(k)
         mat(k, 897) = mat(k, 897) - dti(k)
         mat(k, 903) = mat(k, 903) - dti(k)
         mat(k, 913) = mat(k, 913) - dti(k)
         mat(k, 924) = mat(k, 924) - dti(k)
         mat(k, 939) = mat(k, 939) - dti(k)
         mat(k, 952) = mat(k, 952) - dti(k)
         mat(k, 962) = mat(k, 962) - dti(k)
         mat(k, 970) = mat(k, 970) - dti(k)
         mat(k, 988) = mat(k, 988) - dti(k)
         mat(k,1010) = mat(k,1010) - dti(k)
         mat(k,1016) = mat(k,1016) - dti(k)
         mat(k,1039) = mat(k,1039) - dti(k)
         mat(k,1064) = mat(k,1064) - dti(k)
         mat(k,1085) = mat(k,1085) - dti(k)
         mat(k,1101) = mat(k,1101) - dti(k)
         mat(k,1115) = mat(k,1115) - dti(k)
         mat(k,1125) = mat(k,1125) - dti(k)
         mat(k,1135) = mat(k,1135) - dti(k)
         mat(k,1142) = mat(k,1142) - dti(k)
         mat(k,1155) = mat(k,1155) - dti(k)
         mat(k,1169) = mat(k,1169) - dti(k)
         mat(k,1180) = mat(k,1180) - dti(k)
         mat(k,1193) = mat(k,1193) - dti(k)
         mat(k,1213) = mat(k,1213) - dti(k)
         mat(k,1233) = mat(k,1233) - dti(k)
         mat(k,1249) = mat(k,1249) - dti(k)
         mat(k,1261) = mat(k,1261) - dti(k)
         mat(k,1272) = mat(k,1272) - dti(k)
         mat(k,1296) = mat(k,1296) - dti(k)
         mat(k,1328) = mat(k,1328) - dti(k)
         mat(k,1351) = mat(k,1351) - dti(k)
         mat(k,1372) = mat(k,1372) - dti(k)
         mat(k,1394) = mat(k,1394) - dti(k)
         mat(k,1426) = mat(k,1426) - dti(k)
         mat(k,1441) = mat(k,1441) - dti(k)
         mat(k,1455) = mat(k,1455) - dti(k)
         mat(k,1470) = mat(k,1470) - dti(k)
         mat(k,1487) = mat(k,1487) - dti(k)
         mat(k,1503) = mat(k,1503) - dti(k)
         mat(k,1604) = mat(k,1604) - dti(k)
         mat(k,1657) = mat(k,1657) - dti(k)
         mat(k,1702) = mat(k,1702) - dti(k)
         mat(k,1744) = mat(k,1744) - dti(k)
         mat(k,1917) = mat(k,1917) - dti(k)
         mat(k,1944) = mat(k,1944) - dti(k)
         mat(k,1967) = mat(k,1967) - dti(k)
         mat(k,1999) = mat(k,1999) - dti(k)
         mat(k,2024) = mat(k,2024) - dti(k)
         mat(k,2082) = mat(k,2082) - dti(k)
         mat(k,2107) = mat(k,2107) - dti(k)
         mat(k,2227) = mat(k,2227) - dti(k)
         mat(k,2267) = mat(k,2267) - dti(k)
         mat(k,2329) = mat(k,2329) - dti(k)
         mat(k,2356) = mat(k,2356) - dti(k)
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
