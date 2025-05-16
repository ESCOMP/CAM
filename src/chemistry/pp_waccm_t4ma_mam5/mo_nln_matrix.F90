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
         mat(k,32) = -(rxt(k,341)*y(k,111))
         mat(k,1330) = -rxt(k,341)*y(k,3)
         mat(k,876) = -(rxt(k,215)*y(k,27) + rxt(k,216)*y(k,139) + rxt(k,217)*y(k,108))
         mat(k,1016) = -rxt(k,215)*y(k,4)
         mat(k,1075) = -rxt(k,216)*y(k,4)
         mat(k,1121) = -rxt(k,217)*y(k,4)
         mat(k,930) = 4.000_r8*rxt(k,218)*y(k,6) + (rxt(k,219)+rxt(k,220))*y(k,45) &
                      + rxt(k,223)*y(k,97) + rxt(k,226)*y(k,106) + rxt(k,227)*y(k,111) &
                      + rxt(k,373)*y(k,122)
         mat(k,73) = rxt(k,205)*y(k,144)
         mat(k,79) = rxt(k,231)*y(k,144)
         mat(k,248) = 2.000_r8*rxt(k,242)*y(k,42) + 2.000_r8*rxt(k,243)*y(k,111) &
                      + 2.000_r8*rxt(k,254)*y(k,144)
         mat(k,327) = rxt(k,244)*y(k,42) + rxt(k,245)*y(k,111) + rxt(k,255)*y(k,144)
         mat(k,223) = 3.000_r8*rxt(k,249)*y(k,42) + 3.000_r8*rxt(k,250)*y(k,111) &
                      + 3.000_r8*rxt(k,232)*y(k,144)
         mat(k,1311) = 2.000_r8*rxt(k,242)*y(k,26) + rxt(k,244)*y(k,28) &
                      + 3.000_r8*rxt(k,249)*y(k,41)
         mat(k,956) = (rxt(k,219)+rxt(k,220))*y(k,6)
         mat(k,43) = 2.000_r8*rxt(k,233)*y(k,144)
         mat(k,434) = rxt(k,228)*y(k,106) + rxt(k,229)*y(k,111) + rxt(k,234)*y(k,144)
         mat(k,1493) = rxt(k,223)*y(k,6)
         mat(k,1445) = rxt(k,226)*y(k,6) + rxt(k,228)*y(k,67)
         mat(k,1400) = rxt(k,227)*y(k,6) + 2.000_r8*rxt(k,243)*y(k,26) + rxt(k,245) &
                      *y(k,28) + 3.000_r8*rxt(k,250)*y(k,41) + rxt(k,229)*y(k,67)
         mat(k,842) = rxt(k,373)*y(k,6)
         mat(k,1273) = rxt(k,205)*y(k,19) + rxt(k,231)*y(k,20) + 2.000_r8*rxt(k,254) &
                      *y(k,26) + rxt(k,255)*y(k,28) + 3.000_r8*rxt(k,232)*y(k,41) &
                      + 2.000_r8*rxt(k,233)*y(k,64) + rxt(k,234)*y(k,67)
         mat(k,923) = rxt(k,221)*y(k,45)
         mat(k,949) = rxt(k,221)*y(k,6)
         mat(k,856) = (rxt(k,399)+rxt(k,404))*y(k,75)
         mat(k,385) = (rxt(k,399)+rxt(k,404))*y(k,71)
         mat(k,932) = -(4._r8*rxt(k,218)*y(k,6) + (rxt(k,219) + rxt(k,220) + rxt(k,221) &
                      ) * y(k,45) + rxt(k,222)*y(k,139) + rxt(k,223)*y(k,97) + rxt(k,224) &
                      *y(k,98) + rxt(k,226)*y(k,106) + rxt(k,227)*y(k,111) + rxt(k,373) &
                      *y(k,122))
         mat(k,958) = -(rxt(k,219) + rxt(k,220) + rxt(k,221)) * y(k,6)
         mat(k,1077) = -rxt(k,222)*y(k,6)
         mat(k,1495) = -rxt(k,223)*y(k,6)
         mat(k,994) = -rxt(k,224)*y(k,6)
         mat(k,1447) = -rxt(k,226)*y(k,6)
         mat(k,1402) = -rxt(k,227)*y(k,6)
         mat(k,844) = -rxt(k,373)*y(k,6)
         mat(k,878) = rxt(k,217)*y(k,108)
         mat(k,297) = rxt(k,225)*y(k,106)
         mat(k,435) = rxt(k,235)*y(k,144)
         mat(k,389) = rxt(k,230)*y(k,106)
         mat(k,1447) = mat(k,1447) + rxt(k,225)*y(k,7) + rxt(k,230)*y(k,75)
         mat(k,1123) = rxt(k,217)*y(k,4)
         mat(k,1275) = rxt(k,235)*y(k,67)
         mat(k,294) = -(rxt(k,225)*y(k,106))
         mat(k,1424) = -rxt(k,225)*y(k,7)
         mat(k,925) = rxt(k,224)*y(k,98)
         mat(k,978) = rxt(k,224)*y(k,6)
         mat(k,239) = -(rxt(k,273)*y(k,42) + rxt(k,274)*y(k,108) + rxt(k,298)*y(k,111))
         mat(k,1296) = -rxt(k,273)*y(k,9)
         mat(k,1096) = -rxt(k,274)*y(k,9)
         mat(k,1355) = -rxt(k,298)*y(k,9)
         mat(k,424) = -(4._r8*rxt(k,275)*y(k,10) + rxt(k,276)*y(k,37) + rxt(k,277) &
                      *y(k,139) + rxt(k,278)*y(k,97))
         mat(k,1170) = -rxt(k,276)*y(k,10)
         mat(k,1054) = -rxt(k,277)*y(k,10)
         mat(k,1472) = -rxt(k,278)*y(k,10)
         mat(k,169) = .500_r8*rxt(k,280)*y(k,111)
         mat(k,158) = rxt(k,281)*y(k,42) + rxt(k,282)*y(k,111)
         mat(k,1301) = rxt(k,281)*y(k,13)
         mat(k,1372) = .500_r8*rxt(k,280)*y(k,12) + rxt(k,282)*y(k,13)
         mat(k,126) = -(rxt(k,279)*y(k,111))
         mat(k,1339) = -rxt(k,279)*y(k,11)
         mat(k,421) = .800_r8*rxt(k,275)*y(k,10) + .200_r8*rxt(k,276)*y(k,37)
         mat(k,1163) = .200_r8*rxt(k,276)*y(k,10)
         mat(k,168) = -(rxt(k,280)*y(k,111))
         mat(k,1345) = -rxt(k,280)*y(k,12)
         mat(k,422) = rxt(k,277)*y(k,139)
         mat(k,1038) = rxt(k,277)*y(k,10)
         mat(k,157) = -(rxt(k,281)*y(k,42) + rxt(k,282)*y(k,111))
         mat(k,1293) = -rxt(k,281)*y(k,13)
         mat(k,1343) = -rxt(k,282)*y(k,13)
         mat(k,623) = -(rxt(k,301)*y(k,99) + rxt(k,302)*y(k,108) + rxt(k,319)*y(k,111))
         mat(k,1218) = -rxt(k,301)*y(k,14)
         mat(k,1108) = -rxt(k,302)*y(k,14)
         mat(k,1385) = -rxt(k,319)*y(k,14)
         mat(k,475) = .130_r8*rxt(k,352)*y(k,108)
         mat(k,1108) = mat(k,1108) + .130_r8*rxt(k,352)*y(k,79)
         mat(k,460) = -(rxt(k,303)*y(k,37) + rxt(k,304)*y(k,139) + rxt(k,305)*y(k,97))
         mat(k,1171) = -rxt(k,303)*y(k,15)
         mat(k,1057) = -rxt(k,304)*y(k,15)
         mat(k,1473) = -rxt(k,305)*y(k,15)
         mat(k,33) = 1.670_r8*rxt(k,341)*y(k,111)
         mat(k,205) = rxt(k,306)*y(k,111)
         mat(k,50) = rxt(k,307)*y(k,111)
         mat(k,1375) = 1.670_r8*rxt(k,341)*y(k,3) + rxt(k,306)*y(k,16) + rxt(k,307) &
                      *y(k,17)
         mat(k,203) = -(rxt(k,306)*y(k,111))
         mat(k,1350) = -rxt(k,306)*y(k,16)
         mat(k,458) = rxt(k,304)*y(k,139)
         mat(k,1039) = rxt(k,304)*y(k,15)
         mat(k,49) = -(rxt(k,307)*y(k,111))
         mat(k,1333) = -rxt(k,307)*y(k,17)
         mat(k,38) = -(rxt(k,204)*y(k,144))
         mat(k,1250) = -rxt(k,204)*y(k,18)
         mat(k,71) = -(rxt(k,205)*y(k,144))
         mat(k,1255) = -rxt(k,205)*y(k,19)
         mat(k,76) = -(rxt(k,231)*y(k,144))
         mat(k,1256) = -rxt(k,231)*y(k,20)
         mat(k,53) = -(rxt(k,206)*y(k,144))
         mat(k,1252) = -rxt(k,206)*y(k,21)
         mat(k,81) = -(rxt(k,207)*y(k,144))
         mat(k,1257) = -rxt(k,207)*y(k,22)
         mat(k,57) = -(rxt(k,208)*y(k,144))
         mat(k,1253) = -rxt(k,208)*y(k,23)
         mat(k,86) = -(rxt(k,209)*y(k,144))
         mat(k,1258) = -rxt(k,209)*y(k,24)
         mat(k,61) = -(rxt(k,210)*y(k,144))
         mat(k,1254) = -rxt(k,210)*y(k,25)
         mat(k,246) = -(rxt(k,242)*y(k,42) + rxt(k,243)*y(k,111) + rxt(k,254)*y(k,144))
         mat(k,1297) = -rxt(k,242)*y(k,26)
         mat(k,1356) = -rxt(k,243)*y(k,26)
         mat(k,1267) = -rxt(k,254)*y(k,26)
         mat(k,1021) = -(rxt(k,179)*y(k,42) + rxt(k,215)*y(k,4) + rxt(k,259)*y(k,99) &
                      + rxt(k,260)*y(k,106) + rxt(k,261)*y(k,111))
         mat(k,1316) = -rxt(k,179)*y(k,27)
         mat(k,880) = -rxt(k,215)*y(k,27)
         mat(k,1236) = -rxt(k,259)*y(k,27)
         mat(k,1450) = -rxt(k,260)*y(k,27)
         mat(k,1405) = -rxt(k,261)*y(k,27)
         mat(k,242) = rxt(k,274)*y(k,108)
         mat(k,428) = .700_r8*rxt(k,276)*y(k,37)
         mat(k,631) = .500_r8*rxt(k,302)*y(k,108)
         mat(k,466) = rxt(k,303)*y(k,37)
         mat(k,797) = rxt(k,286)*y(k,37) + .600_r8*rxt(k,347)*y(k,81) &
                      + .250_r8*rxt(k,321)*y(k,84) + rxt(k,330)*y(k,86) &
                      + .250_r8*rxt(k,357)*y(k,134)
         mat(k,321) = .500_r8*rxt(k,290)*y(k,111)
         mat(k,1190) = .700_r8*rxt(k,276)*y(k,10) + rxt(k,303)*y(k,15) + rxt(k,286) &
                      *y(k,32) + (4.000_r8*rxt(k,262)+2.000_r8*rxt(k,263))*y(k,37) &
                      + rxt(k,186)*y(k,45) + 1.200_r8*rxt(k,348)*y(k,81) &
                      + .880_r8*rxt(k,322)*y(k,84) + 2.000_r8*rxt(k,331)*y(k,86) &
                      + rxt(k,265)*y(k,97) + .800_r8*rxt(k,315)*y(k,118) &
                      + .800_r8*rxt(k,358)*y(k,134)
         mat(k,283) = rxt(k,266)*y(k,111)
         mat(k,210) = .300_r8*rxt(k,267)*y(k,111)
         mat(k,1523) = (rxt(k,270)+rxt(k,271))*y(k,144)
         mat(k,961) = rxt(k,186)*y(k,37)
         mat(k,381) = .500_r8*rxt(k,292)*y(k,97)
         mat(k,610) = .800_r8*rxt(k,295)*y(k,111)
         mat(k,482) = .910_r8*rxt(k,352)*y(k,108)
         mat(k,416) = .072_r8*rxt(k,345)*y(k,97) + .072_r8*rxt(k,346)*y(k,99) &
                      + .206_r8*rxt(k,344)*y(k,139)
         mat(k,709) = .600_r8*rxt(k,347)*y(k,32) + 1.200_r8*rxt(k,348)*y(k,37) &
                      + .550_r8*rxt(k,350)*y(k,97) + .600_r8*rxt(k,351)*y(k,99)
         mat(k,660) = .120_r8*rxt(k,327)*y(k,108)
         mat(k,731) = .250_r8*rxt(k,321)*y(k,32) + .880_r8*rxt(k,322)*y(k,37) &
                      + .250_r8*rxt(k,325)*y(k,97) + .250_r8*rxt(k,324)*y(k,99)
         mat(k,769) = rxt(k,330)*y(k,32) + 2.000_r8*rxt(k,331)*y(k,37) &
                      + 4.000_r8*rxt(k,333)*y(k,86) + rxt(k,334)*y(k,97) + rxt(k,335) &
                      *y(k,99) + .450_r8*rxt(k,332)*y(k,139)
         mat(k,307) = .500_r8*rxt(k,336)*y(k,111)
         mat(k,750) = .600_r8*rxt(k,337)*y(k,108)
         mat(k,1498) = rxt(k,265)*y(k,37) + .500_r8*rxt(k,292)*y(k,57) &
                      + .072_r8*rxt(k,345)*y(k,80) + .550_r8*rxt(k,350)*y(k,81) &
                      + .250_r8*rxt(k,325)*y(k,84) + rxt(k,334)*y(k,86) + rxt(k,313) &
                      *y(k,114) + rxt(k,317)*y(k,118) + .250_r8*rxt(k,360)*y(k,134)
         mat(k,1236) = mat(k,1236) + .072_r8*rxt(k,346)*y(k,80) + .600_r8*rxt(k,351) &
                      *y(k,81) + .250_r8*rxt(k,324)*y(k,84) + rxt(k,335)*y(k,86)
         mat(k,1126) = rxt(k,274)*y(k,9) + .500_r8*rxt(k,302)*y(k,14) &
                      + .910_r8*rxt(k,352)*y(k,79) + .120_r8*rxt(k,327)*y(k,83) &
                      + .600_r8*rxt(k,337)*y(k,88)
         mat(k,1405) = mat(k,1405) + .500_r8*rxt(k,290)*y(k,36) + rxt(k,266)*y(k,38) &
                      + .300_r8*rxt(k,267)*y(k,39) + .800_r8*rxt(k,295)*y(k,60) &
                      + .500_r8*rxt(k,336)*y(k,87) + rxt(k,297)*y(k,113)
         mat(k,278) = rxt(k,297)*y(k,111)
         mat(k,406) = rxt(k,313)*y(k,97)
         mat(k,600) = .800_r8*rxt(k,315)*y(k,37) + rxt(k,317)*y(k,97) &
                      + .150_r8*rxt(k,316)*y(k,139)
         mat(k,676) = .250_r8*rxt(k,357)*y(k,32) + .800_r8*rxt(k,358)*y(k,37) &
                      + .250_r8*rxt(k,360)*y(k,97)
         mat(k,1080) = .206_r8*rxt(k,344)*y(k,80) + .450_r8*rxt(k,332)*y(k,86) &
                      + .150_r8*rxt(k,316)*y(k,118)
         mat(k,1278) = (rxt(k,270)+rxt(k,271))*y(k,40)
         mat(k,325) = -(rxt(k,244)*y(k,42) + rxt(k,245)*y(k,111) + rxt(k,255)*y(k,144))
         mat(k,1299) = -rxt(k,244)*y(k,28)
         mat(k,1364) = -rxt(k,245)*y(k,28)
         mat(k,1268) = -rxt(k,255)*y(k,28)
         mat(k,45) = -(rxt(k,246)*y(k,111))
         mat(k,1332) = -rxt(k,246)*y(k,29)
         mat(k,640) = -(rxt(k,283)*y(k,99) + rxt(k,284)*y(k,111))
         mat(k,1219) = -rxt(k,283)*y(k,30)
         mat(k,1386) = -rxt(k,284)*y(k,30)
         mat(k,425) = 3.200_r8*rxt(k,275)*y(k,10) + .800_r8*rxt(k,276)*y(k,37) &
                      + rxt(k,278)*y(k,97)
         mat(k,127) = rxt(k,279)*y(k,111)
         mat(k,170) = .500_r8*rxt(k,280)*y(k,111)
         mat(k,624) = .500_r8*rxt(k,302)*y(k,108)
         mat(k,462) = .270_r8*rxt(k,305)*y(k,97)
         mat(k,1175) = .800_r8*rxt(k,276)*y(k,10)
         mat(k,740) = .100_r8*rxt(k,337)*y(k,108)
         mat(k,1482) = rxt(k,278)*y(k,10) + .270_r8*rxt(k,305)*y(k,15) + rxt(k,313) &
                      *y(k,114)
         mat(k,1109) = .500_r8*rxt(k,302)*y(k,14) + .100_r8*rxt(k,337)*y(k,88)
         mat(k,1386) = mat(k,1386) + rxt(k,279)*y(k,11) + .500_r8*rxt(k,280)*y(k,12)
         mat(k,403) = rxt(k,313)*y(k,97)
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
         mat(k,189) = -(rxt(k,247)*y(k,42) + rxt(k,248)*y(k,111))
         mat(k,1294) = -rxt(k,247)*y(k,31)
         mat(k,1348) = -rxt(k,248)*y(k,31)
         mat(k,793) = -(4._r8*rxt(k,285)*y(k,32) + rxt(k,286)*y(k,37) + rxt(k,287) &
                      *y(k,139) + rxt(k,288)*y(k,97) + rxt(k,299)*y(k,98) + rxt(k,321) &
                      *y(k,84) + rxt(k,347)*y(k,81) + rxt(k,357)*y(k,134))
         mat(k,1184) = -rxt(k,286)*y(k,32)
         mat(k,1071) = -rxt(k,287)*y(k,32)
         mat(k,1491) = -rxt(k,288)*y(k,32)
         mat(k,988) = -rxt(k,299)*y(k,32)
         mat(k,728) = -rxt(k,321)*y(k,32)
         mat(k,706) = -rxt(k,347)*y(k,32)
         mat(k,673) = -rxt(k,357)*y(k,32)
         mat(k,642) = rxt(k,283)*y(k,99) + rxt(k,284)*y(k,111)
         mat(k,793) = mat(k,793) + .530_r8*rxt(k,321)*y(k,84)
         mat(k,684) = rxt(k,308)*y(k,99) + rxt(k,309)*y(k,111)
         mat(k,319) = .500_r8*rxt(k,290)*y(k,111)
         mat(k,1184) = mat(k,1184) + .260_r8*rxt(k,322)*y(k,84) + rxt(k,331)*y(k,86) &
                      + .300_r8*rxt(k,315)*y(k,118)
         mat(k,480) = .080_r8*rxt(k,352)*y(k,108)
         mat(k,658) = .100_r8*rxt(k,327)*y(k,108)
         mat(k,728) = mat(k,728) + .530_r8*rxt(k,321)*y(k,32) + .260_r8*rxt(k,322) &
                      *y(k,37) + .530_r8*rxt(k,325)*y(k,97) + .530_r8*rxt(k,324) &
                      *y(k,99)
         mat(k,766) = rxt(k,331)*y(k,37) + 4.000_r8*rxt(k,333)*y(k,86) + rxt(k,334) &
                      *y(k,97) + rxt(k,335)*y(k,99) + .450_r8*rxt(k,332)*y(k,139)
         mat(k,746) = .280_r8*rxt(k,337)*y(k,108)
         mat(k,1491) = mat(k,1491) + .530_r8*rxt(k,325)*y(k,84) + rxt(k,334)*y(k,86) &
                      + rxt(k,317)*y(k,118)
         mat(k,1228) = rxt(k,283)*y(k,30) + rxt(k,308)*y(k,34) + .530_r8*rxt(k,324) &
                      *y(k,84) + rxt(k,335)*y(k,86)
         mat(k,1118) = .080_r8*rxt(k,352)*y(k,79) + .100_r8*rxt(k,327)*y(k,83) &
                      + .280_r8*rxt(k,337)*y(k,88)
         mat(k,1395) = rxt(k,284)*y(k,30) + rxt(k,309)*y(k,34) + .500_r8*rxt(k,290) &
                      *y(k,36)
         mat(k,598) = .300_r8*rxt(k,315)*y(k,37) + rxt(k,317)*y(k,97) &
                      + .150_r8*rxt(k,316)*y(k,139)
         mat(k,1071) = mat(k,1071) + .450_r8*rxt(k,332)*y(k,86) + .150_r8*rxt(k,316) &
                      *y(k,118)
         mat(k,394) = -(rxt(k,320)*y(k,111))
         mat(k,1370) = -rxt(k,320)*y(k,33)
         mat(k,459) = .820_r8*rxt(k,303)*y(k,37) + .820_r8*rxt(k,305)*y(k,97)
         mat(k,1169) = .820_r8*rxt(k,303)*y(k,15)
         mat(k,1469) = .820_r8*rxt(k,305)*y(k,15)
         mat(k,1370) = mat(k,1370) + .100_r8*rxt(k,365)*y(k,133)
         mat(k,344) = .100_r8*rxt(k,365)*y(k,111)
         mat(k,683) = -(rxt(k,308)*y(k,99) + rxt(k,309)*y(k,111))
         mat(k,1223) = -rxt(k,308)*y(k,34)
         mat(k,1390) = -rxt(k,309)*y(k,34)
         mat(k,788) = .250_r8*rxt(k,321)*y(k,84) + .250_r8*rxt(k,357)*y(k,134)
         mat(k,1179) = .240_r8*rxt(k,322)*y(k,84) + .500_r8*rxt(k,315)*y(k,118) &
                      + .100_r8*rxt(k,358)*y(k,134)
         mat(k,614) = rxt(k,310)*y(k,111)
         mat(k,701) = .020_r8*rxt(k,350)*y(k,97)
         mat(k,655) = .880_r8*rxt(k,327)*y(k,108)
         mat(k,725) = .250_r8*rxt(k,321)*y(k,32) + .240_r8*rxt(k,322)*y(k,37) &
                      + .250_r8*rxt(k,325)*y(k,97) + .250_r8*rxt(k,324)*y(k,99)
         mat(k,742) = .500_r8*rxt(k,337)*y(k,108)
         mat(k,1486) = .020_r8*rxt(k,350)*y(k,81) + .250_r8*rxt(k,325)*y(k,84) &
                      + .250_r8*rxt(k,360)*y(k,134)
         mat(k,1223) = mat(k,1223) + .250_r8*rxt(k,324)*y(k,84) + .250_r8*rxt(k,361) &
                      *y(k,134)
         mat(k,198) = rxt(k,311)*y(k,111)
         mat(k,1113) = .880_r8*rxt(k,327)*y(k,83) + .500_r8*rxt(k,337)*y(k,88)
         mat(k,1390) = mat(k,1390) + rxt(k,310)*y(k,77) + rxt(k,311)*y(k,100)
         mat(k,597) = .500_r8*rxt(k,315)*y(k,37)
         mat(k,672) = .250_r8*rxt(k,357)*y(k,32) + .100_r8*rxt(k,358)*y(k,37) &
                      + .250_r8*rxt(k,360)*y(k,97) + .250_r8*rxt(k,361)*y(k,99)
         mat(k,353) = -(rxt(k,289)*y(k,111))
         mat(k,1366) = -rxt(k,289)*y(k,35)
         mat(k,620) = .120_r8*rxt(k,302)*y(k,108)
         mat(k,781) = .100_r8*rxt(k,286)*y(k,37) + .150_r8*rxt(k,287)*y(k,139)
         mat(k,1167) = .100_r8*rxt(k,286)*y(k,32)
         mat(k,762) = .150_r8*rxt(k,332)*y(k,139)
         mat(k,1098) = .120_r8*rxt(k,302)*y(k,14)
         mat(k,1048) = .150_r8*rxt(k,287)*y(k,32) + .150_r8*rxt(k,332)*y(k,86)
         mat(k,318) = -(rxt(k,290)*y(k,111))
         mat(k,1363) = -rxt(k,290)*y(k,36)
         mat(k,780) = .400_r8*rxt(k,287)*y(k,139)
         mat(k,761) = .400_r8*rxt(k,332)*y(k,139)
         mat(k,1047) = .400_r8*rxt(k,287)*y(k,32) + .400_r8*rxt(k,332)*y(k,86)
         mat(k,1194) = -(rxt(k,186)*y(k,45) + (4._r8*rxt(k,262) + 4._r8*rxt(k,263) &
                      ) * y(k,37) + rxt(k,264)*y(k,139) + rxt(k,265)*y(k,97) + rxt(k,276) &
                      *y(k,10) + rxt(k,286)*y(k,32) + rxt(k,303)*y(k,15) + rxt(k,315) &
                      *y(k,118) + rxt(k,322)*y(k,84) + rxt(k,331)*y(k,86) + rxt(k,348) &
                      *y(k,81) + rxt(k,358)*y(k,134))
         mat(k,965) = -rxt(k,186)*y(k,37)
         mat(k,1084) = -rxt(k,264)*y(k,37)
         mat(k,1502) = -rxt(k,265)*y(k,37)
         mat(k,430) = -rxt(k,276)*y(k,37)
         mat(k,801) = -rxt(k,286)*y(k,37)
         mat(k,468) = -rxt(k,303)*y(k,37)
         mat(k,602) = -rxt(k,315)*y(k,37)
         mat(k,734) = -rxt(k,322)*y(k,37)
         mat(k,773) = -rxt(k,331)*y(k,37)
         mat(k,713) = -rxt(k,348)*y(k,37)
         mat(k,678) = -rxt(k,358)*y(k,37)
         mat(k,634) = .280_r8*rxt(k,302)*y(k,108)
         mat(k,801) = mat(k,801) + 4.000_r8*rxt(k,285)*y(k,32) + .900_r8*rxt(k,286) &
                      *y(k,37) + rxt(k,347)*y(k,81) + rxt(k,321)*y(k,84) + rxt(k,330) &
                      *y(k,86) + rxt(k,288)*y(k,97) + rxt(k,357)*y(k,134) &
                      + .450_r8*rxt(k,287)*y(k,139)
         mat(k,355) = rxt(k,289)*y(k,111)
         mat(k,1194) = mat(k,1194) + .900_r8*rxt(k,286)*y(k,32)
         mat(k,211) = .700_r8*rxt(k,267)*y(k,111)
         mat(k,1527) = rxt(k,180)*y(k,42) + rxt(k,236)*y(k,59) + rxt(k,268)*y(k,111) &
                      + rxt(k,269)*y(k,144)
         mat(k,1320) = rxt(k,180)*y(k,40)
         mat(k,493) = rxt(k,236)*y(k,40)
         mat(k,485) = .050_r8*rxt(k,352)*y(k,108)
         mat(k,713) = mat(k,713) + rxt(k,347)*y(k,32)
         mat(k,734) = mat(k,734) + rxt(k,321)*y(k,32)
         mat(k,773) = mat(k,773) + rxt(k,330)*y(k,32)
         mat(k,1502) = mat(k,1502) + rxt(k,288)*y(k,32)
         mat(k,1130) = .280_r8*rxt(k,302)*y(k,14) + .050_r8*rxt(k,352)*y(k,79)
         mat(k,1409) = rxt(k,289)*y(k,35) + .700_r8*rxt(k,267)*y(k,39) + rxt(k,268) &
                      *y(k,40)
         mat(k,678) = mat(k,678) + rxt(k,357)*y(k,32)
         mat(k,1084) = mat(k,1084) + .450_r8*rxt(k,287)*y(k,32)
         mat(k,1282) = rxt(k,269)*y(k,40)
         mat(k,282) = -(rxt(k,266)*y(k,111))
         mat(k,1359) = -rxt(k,266)*y(k,38)
         mat(k,423) = .300_r8*rxt(k,276)*y(k,37)
         mat(k,1166) = .300_r8*rxt(k,276)*y(k,10) + 2.000_r8*rxt(k,263)*y(k,37) &
                      + .250_r8*rxt(k,348)*y(k,81) + .250_r8*rxt(k,322)*y(k,84) &
                      + .500_r8*rxt(k,315)*y(k,118) + .300_r8*rxt(k,358)*y(k,134)
         mat(k,693) = .250_r8*rxt(k,348)*y(k,37)
         mat(k,719) = .250_r8*rxt(k,322)*y(k,37)
         mat(k,594) = .500_r8*rxt(k,315)*y(k,37)
         mat(k,666) = .300_r8*rxt(k,358)*y(k,37)
         mat(k,209) = -(rxt(k,267)*y(k,111))
         mat(k,1351) = -rxt(k,267)*y(k,39)
         mat(k,1165) = rxt(k,264)*y(k,139)
         mat(k,1040) = rxt(k,264)*y(k,37)
         mat(k,1534) = -(rxt(k,180)*y(k,42) + rxt(k,236)*y(k,59) + rxt(k,268)*y(k,111) &
                      + (rxt(k,269) + rxt(k,270) + rxt(k,271)) * y(k,144))
         mat(k,1327) = -rxt(k,180)*y(k,40)
         mat(k,496) = -rxt(k,236)*y(k,40)
         mat(k,1416) = -rxt(k,268)*y(k,40)
         mat(k,1289) = -(rxt(k,269) + rxt(k,270) + rxt(k,271)) * y(k,40)
         mat(k,638) = .100_r8*rxt(k,302)*y(k,108)
         mat(k,1137) = .100_r8*rxt(k,302)*y(k,14)
         mat(k,221) = -(rxt(k,232)*y(k,144) + rxt(k,249)*y(k,42) + rxt(k,250)*y(k,111))
         mat(k,1266) = -rxt(k,232)*y(k,41)
         mat(k,1295) = -rxt(k,249)*y(k,41)
         mat(k,1353) = -rxt(k,250)*y(k,41)
         mat(k,1323) = -(rxt(k,179)*y(k,27) + rxt(k,180)*y(k,40) + rxt(k,181)*y(k,63) &
                      + rxt(k,182)*y(k,65) + (rxt(k,183) + rxt(k,184)) * y(k,139) &
                      + rxt(k,185)*y(k,108) + rxt(k,192)*y(k,46) + rxt(k,201)*y(k,76) &
                      + rxt(k,242)*y(k,26) + rxt(k,244)*y(k,28) + rxt(k,247)*y(k,31) &
                      + rxt(k,249)*y(k,41) + rxt(k,281)*y(k,13))
         mat(k,1028) = -rxt(k,179)*y(k,42)
         mat(k,1530) = -rxt(k,180)*y(k,42)
         mat(k,832) = -rxt(k,181)*y(k,42)
         mat(k,314) = -rxt(k,182)*y(k,42)
         mat(k,1087) = -(rxt(k,183) + rxt(k,184)) * y(k,42)
         mat(k,1133) = -rxt(k,185)*y(k,42)
         mat(k,579) = -rxt(k,192)*y(k,42)
         mat(k,454) = -rxt(k,201)*y(k,42)
         mat(k,250) = -rxt(k,242)*y(k,42)
         mat(k,331) = -rxt(k,244)*y(k,42)
         mat(k,194) = -rxt(k,247)*y(k,42)
         mat(k,225) = -rxt(k,249)*y(k,42)
         mat(k,160) = -rxt(k,281)*y(k,42)
         mat(k,941) = rxt(k,220)*y(k,45)
         mat(k,40) = 4.000_r8*rxt(k,204)*y(k,144)
         mat(k,75) = rxt(k,205)*y(k,144)
         mat(k,56) = 2.000_r8*rxt(k,206)*y(k,144)
         mat(k,85) = 2.000_r8*rxt(k,207)*y(k,144)
         mat(k,60) = 2.000_r8*rxt(k,208)*y(k,144)
         mat(k,90) = rxt(k,209)*y(k,144)
         mat(k,64) = 2.000_r8*rxt(k,210)*y(k,144)
         mat(k,46) = 3.000_r8*rxt(k,246)*y(k,111)
         mat(k,194) = mat(k,194) + rxt(k,248)*y(k,111)
         mat(k,1197) = rxt(k,186)*y(k,45)
         mat(k,968) = rxt(k,220)*y(k,6) + rxt(k,186)*y(k,37) + (4.000_r8*rxt(k,187) &
                       +2.000_r8*rxt(k,189))*y(k,45) + rxt(k,191)*y(k,97) + rxt(k,196) &
                      *y(k,106) + rxt(k,197)*y(k,111) + rxt(k,374)*y(k,122)
         mat(k,115) = rxt(k,241)*y(k,144)
         mat(k,110) = rxt(k,251)*y(k,111) + rxt(k,256)*y(k,144)
         mat(k,120) = rxt(k,252)*y(k,111) + rxt(k,257)*y(k,144)
         mat(k,145) = rxt(k,253)*y(k,111) + rxt(k,258)*y(k,144)
         mat(k,868) = rxt(k,199)*y(k,106) + rxt(k,200)*y(k,111) + rxt(k,211)*y(k,144)
         mat(k,1505) = rxt(k,191)*y(k,45)
         mat(k,1457) = rxt(k,196)*y(k,45) + rxt(k,199)*y(k,71)
         mat(k,1412) = 3.000_r8*rxt(k,246)*y(k,29) + rxt(k,248)*y(k,31) + rxt(k,197) &
                      *y(k,45) + rxt(k,251)*y(k,68) + rxt(k,252)*y(k,69) + rxt(k,253) &
                      *y(k,70) + rxt(k,200)*y(k,71)
         mat(k,849) = rxt(k,374)*y(k,45)
         mat(k,1285) = 4.000_r8*rxt(k,204)*y(k,18) + rxt(k,205)*y(k,19) &
                      + 2.000_r8*rxt(k,206)*y(k,21) + 2.000_r8*rxt(k,207)*y(k,22) &
                      + 2.000_r8*rxt(k,208)*y(k,23) + rxt(k,209)*y(k,24) &
                      + 2.000_r8*rxt(k,210)*y(k,25) + rxt(k,241)*y(k,51) + rxt(k,256) &
                      *y(k,68) + rxt(k,257)*y(k,69) + rxt(k,258)*y(k,70) + rxt(k,211) &
                      *y(k,71)
         mat(k,1292) = rxt(k,192)*y(k,46)
         mat(k,948) = 2.000_r8*rxt(k,188)*y(k,45)
         mat(k,571) = rxt(k,192)*y(k,42) + (rxt(k,397)+rxt(k,402)+rxt(k,407))*y(k,71)
         mat(k,855) = (rxt(k,397)+rxt(k,402)+rxt(k,407))*y(k,46) + (rxt(k,392) &
                       +rxt(k,398)+rxt(k,403))*y(k,76)
         mat(k,450) = (rxt(k,392)+rxt(k,398)+rxt(k,403))*y(k,71)
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
         mat(k,947) = 2.000_r8*rxt(k,213)*y(k,45)
         mat(k,959) = -(rxt(k,186)*y(k,37) + (4._r8*rxt(k,187) + 4._r8*rxt(k,188) &
                      + 4._r8*rxt(k,189) + 4._r8*rxt(k,213)) * y(k,45) + rxt(k,190) &
                      *y(k,139) + rxt(k,191)*y(k,97) + rxt(k,193)*y(k,98) + rxt(k,196) &
                      *y(k,106) + (rxt(k,197) + rxt(k,198)) * y(k,111) + (rxt(k,219) &
                      + rxt(k,220) + rxt(k,221)) * y(k,6) + rxt(k,374)*y(k,122))
         mat(k,1188) = -rxt(k,186)*y(k,45)
         mat(k,1078) = -rxt(k,190)*y(k,45)
         mat(k,1496) = -rxt(k,191)*y(k,45)
         mat(k,995) = -rxt(k,193)*y(k,45)
         mat(k,1448) = -rxt(k,196)*y(k,45)
         mat(k,1403) = -(rxt(k,197) + rxt(k,198)) * y(k,45)
         mat(k,933) = -(rxt(k,219) + rxt(k,220) + rxt(k,221)) * y(k,45)
         mat(k,845) = -rxt(k,374)*y(k,45)
         mat(k,1314) = rxt(k,201)*y(k,76) + rxt(k,185)*y(k,108) + rxt(k,184)*y(k,139)
         mat(k,575) = rxt(k,194)*y(k,106)
         mat(k,863) = rxt(k,212)*y(k,144)
         mat(k,453) = rxt(k,201)*y(k,42) + rxt(k,202)*y(k,106) + rxt(k,203)*y(k,111)
         mat(k,1448) = mat(k,1448) + rxt(k,194)*y(k,46) + rxt(k,202)*y(k,76)
         mat(k,1124) = rxt(k,185)*y(k,42)
         mat(k,176) = rxt(k,379)*y(k,122)
         mat(k,1403) = mat(k,1403) + rxt(k,203)*y(k,76)
         mat(k,845) = mat(k,845) + rxt(k,379)*y(k,109)
         mat(k,1078) = mat(k,1078) + rxt(k,184)*y(k,42)
         mat(k,1276) = rxt(k,212)*y(k,71)
         mat(k,573) = -(rxt(k,192)*y(k,42) + rxt(k,194)*y(k,106) + rxt(k,195)*y(k,111) &
                      + (rxt(k,397) + rxt(k,402) + rxt(k,407)) * y(k,71))
         mat(k,1303) = -rxt(k,192)*y(k,46)
         mat(k,1439) = -rxt(k,194)*y(k,46)
         mat(k,1380) = -rxt(k,195)*y(k,46)
         mat(k,859) = -(rxt(k,397) + rxt(k,402) + rxt(k,407)) * y(k,46)
         mat(k,953) = rxt(k,193)*y(k,98)
         mat(k,983) = rxt(k,193)*y(k,45)
         mat(k,649) = -(rxt(k,272)*y(k,111))
         mat(k,1387) = -rxt(k,272)*y(k,48)
         mat(k,874) = rxt(k,215)*y(k,27)
         mat(k,241) = .630_r8*rxt(k,274)*y(k,108)
         mat(k,625) = .560_r8*rxt(k,302)*y(k,108)
         mat(k,1012) = rxt(k,215)*y(k,4) + rxt(k,179)*y(k,42) + rxt(k,259)*y(k,99) &
                      + rxt(k,260)*y(k,106) + rxt(k,261)*y(k,111)
         mat(k,190) = rxt(k,247)*y(k,42)
         mat(k,785) = .220_r8*rxt(k,321)*y(k,84) + .250_r8*rxt(k,357)*y(k,134)
         mat(k,682) = rxt(k,308)*y(k,99) + rxt(k,309)*y(k,111)
         mat(k,1176) = .110_r8*rxt(k,322)*y(k,84) + .200_r8*rxt(k,358)*y(k,134)
         mat(k,1306) = rxt(k,179)*y(k,27) + rxt(k,247)*y(k,31)
         mat(k,814) = rxt(k,423)*y(k,148)
         mat(k,366) = rxt(k,296)*y(k,111)
         mat(k,476) = .620_r8*rxt(k,352)*y(k,108)
         mat(k,653) = .650_r8*rxt(k,327)*y(k,108)
         mat(k,723) = .220_r8*rxt(k,321)*y(k,32) + .110_r8*rxt(k,322)*y(k,37) &
                      + .220_r8*rxt(k,325)*y(k,97) + .220_r8*rxt(k,324)*y(k,99)
         mat(k,741) = .560_r8*rxt(k,337)*y(k,108)
         mat(k,1483) = .220_r8*rxt(k,325)*y(k,84) + .250_r8*rxt(k,360)*y(k,134)
         mat(k,1220) = rxt(k,259)*y(k,27) + rxt(k,308)*y(k,34) + .220_r8*rxt(k,324) &
                      *y(k,84) + .500_r8*rxt(k,361)*y(k,134)
         mat(k,1440) = rxt(k,260)*y(k,27) + rxt(k,368)*y(k,110)
         mat(k,1110) = .630_r8*rxt(k,274)*y(k,9) + .560_r8*rxt(k,302)*y(k,14) &
                      + .620_r8*rxt(k,352)*y(k,79) + .650_r8*rxt(k,327)*y(k,83) &
                      + .560_r8*rxt(k,337)*y(k,88)
         mat(k,184) = rxt(k,368)*y(k,106) + rxt(k,369)*y(k,111)
         mat(k,1387) = mat(k,1387) + rxt(k,261)*y(k,27) + rxt(k,309)*y(k,34) &
                      + rxt(k,296)*y(k,61) + rxt(k,369)*y(k,110)
         mat(k,670) = .250_r8*rxt(k,357)*y(k,32) + .200_r8*rxt(k,358)*y(k,37) &
                      + .250_r8*rxt(k,360)*y(k,97) + .500_r8*rxt(k,361)*y(k,99)
         mat(k,506) = rxt(k,423)*y(k,49)
         mat(k,815) = -(rxt(k,423)*y(k,148))
         mat(k,507) = -rxt(k,423)*y(k,49)
         mat(k,628) = .200_r8*rxt(k,302)*y(k,108)
         mat(k,794) = 4.000_r8*rxt(k,285)*y(k,32) + .900_r8*rxt(k,286)*y(k,37) &
                      + 2.000_r8*rxt(k,330)*y(k,86) + rxt(k,288)*y(k,97) + rxt(k,357) &
                      *y(k,134)
         mat(k,354) = rxt(k,289)*y(k,111)
         mat(k,320) = .500_r8*rxt(k,290)*y(k,111)
         mat(k,1185) = .900_r8*rxt(k,286)*y(k,32) + rxt(k,331)*y(k,86)
         mat(k,650) = rxt(k,272)*y(k,111)
         mat(k,609) = .800_r8*rxt(k,295)*y(k,111)
         mat(k,367) = rxt(k,296)*y(k,111)
         mat(k,767) = 2.000_r8*rxt(k,330)*y(k,32) + rxt(k,331)*y(k,37) &
                      + 4.000_r8*rxt(k,333)*y(k,86) + .450_r8*rxt(k,332)*y(k,139)
         mat(k,305) = .500_r8*rxt(k,336)*y(k,111)
         mat(k,747) = .100_r8*rxt(k,337)*y(k,108)
         mat(k,1492) = rxt(k,288)*y(k,32)
         mat(k,1119) = .200_r8*rxt(k,302)*y(k,14) + .100_r8*rxt(k,337)*y(k,88)
         mat(k,1396) = rxt(k,289)*y(k,35) + .500_r8*rxt(k,290)*y(k,36) + rxt(k,272) &
                      *y(k,48) + .800_r8*rxt(k,295)*y(k,60) + rxt(k,296)*y(k,61) &
                      + .500_r8*rxt(k,336)*y(k,87)
         mat(k,674) = rxt(k,357)*y(k,32)
         mat(k,1072) = .450_r8*rxt(k,332)*y(k,86)
         mat(k,101) = -(rxt(k,240)*y(k,144))
         mat(k,1260) = -rxt(k,240)*y(k,50)
         mat(k,72) = rxt(k,205)*y(k,144)
         mat(k,77) = rxt(k,231)*y(k,144)
         mat(k,82) = rxt(k,207)*y(k,144)
         mat(k,58) = 2.000_r8*rxt(k,208)*y(k,144)
         mat(k,87) = 2.000_r8*rxt(k,209)*y(k,144)
         mat(k,62) = rxt(k,210)*y(k,144)
         mat(k,42) = 2.000_r8*rxt(k,233)*y(k,144)
         mat(k,116) = rxt(k,252)*y(k,111) + rxt(k,257)*y(k,144)
         mat(k,141) = rxt(k,253)*y(k,111) + rxt(k,258)*y(k,144)
         mat(k,1335) = rxt(k,252)*y(k,69) + rxt(k,253)*y(k,70)
         mat(k,1260) = mat(k,1260) + rxt(k,205)*y(k,19) + rxt(k,231)*y(k,20) &
                      + rxt(k,207)*y(k,22) + 2.000_r8*rxt(k,208)*y(k,23) &
                      + 2.000_r8*rxt(k,209)*y(k,24) + rxt(k,210)*y(k,25) &
                      + 2.000_r8*rxt(k,233)*y(k,64) + rxt(k,257)*y(k,69) + rxt(k,258) &
                      *y(k,70)
         mat(k,112) = -(rxt(k,241)*y(k,144))
         mat(k,1262) = -rxt(k,241)*y(k,51)
         mat(k,54) = rxt(k,206)*y(k,144)
         mat(k,83) = rxt(k,207)*y(k,144)
         mat(k,108) = rxt(k,251)*y(k,111) + rxt(k,256)*y(k,144)
         mat(k,1337) = rxt(k,251)*y(k,68)
         mat(k,1262) = mat(k,1262) + rxt(k,206)*y(k,21) + rxt(k,207)*y(k,22) &
                      + rxt(k,256)*y(k,68)
         mat(k,135) = -(rxt(k,366)*y(k,99) + (rxt(k,367) + rxt(k,381)) * y(k,111))
         mat(k,1204) = -rxt(k,366)*y(k,52)
         mat(k,1341) = -(rxt(k,367) + rxt(k,381)) * y(k,52)
         mat(k,227) = -(rxt(k,294)*y(k,107))
         mat(k,891) = -rxt(k,294)*y(k,56)
         mat(k,376) = .750_r8*rxt(k,292)*y(k,97)
         mat(k,1465) = .750_r8*rxt(k,292)*y(k,57)
         mat(k,377) = -(rxt(k,291)*y(k,139) + rxt(k,292)*y(k,97))
         mat(k,1049) = -rxt(k,291)*y(k,57)
         mat(k,1468) = -rxt(k,292)*y(k,57)
         mat(k,240) = rxt(k,298)*y(k,111)
         mat(k,1369) = rxt(k,298)*y(k,9)
         mat(k,375) = rxt(k,291)*y(k,139)
         mat(k,1035) = rxt(k,291)*y(k,57)
         mat(k,490) = -(rxt(k,236)*y(k,40) + rxt(k,237)*y(k,63) + rxt(k,238)*y(k,151) &
                      + rxt(k,239)*y(k,73))
         mat(k,1513) = -rxt(k,236)*y(k,59)
         mat(k,825) = -rxt(k,237)*y(k,59)
         mat(k,1561) = -rxt(k,238)*y(k,59)
         mat(k,1141) = -rxt(k,239)*y(k,59)
         mat(k,78) = rxt(k,231)*y(k,144)
         mat(k,88) = rxt(k,209)*y(k,144)
         mat(k,102) = 2.000_r8*rxt(k,240)*y(k,144)
         mat(k,113) = rxt(k,241)*y(k,144)
         mat(k,1270) = rxt(k,231)*y(k,20) + rxt(k,209)*y(k,24) + 2.000_r8*rxt(k,240) &
                      *y(k,50) + rxt(k,241)*y(k,51)
         mat(k,607) = -(rxt(k,295)*y(k,111))
         mat(k,1383) = -rxt(k,295)*y(k,60)
         mat(k,783) = .530_r8*rxt(k,321)*y(k,84) + .250_r8*rxt(k,357)*y(k,134)
         mat(k,1173) = .260_r8*rxt(k,322)*y(k,84) + .100_r8*rxt(k,358)*y(k,134)
         mat(k,228) = rxt(k,294)*y(k,107)
         mat(k,696) = .020_r8*rxt(k,350)*y(k,97)
         mat(k,721) = .530_r8*rxt(k,321)*y(k,32) + .260_r8*rxt(k,322)*y(k,37) &
                      + .530_r8*rxt(k,325)*y(k,97) + .530_r8*rxt(k,324)*y(k,99)
         mat(k,1480) = .020_r8*rxt(k,350)*y(k,81) + .530_r8*rxt(k,325)*y(k,84) &
                      + .250_r8*rxt(k,360)*y(k,134)
         mat(k,1216) = .530_r8*rxt(k,324)*y(k,84) + .250_r8*rxt(k,361)*y(k,134)
         mat(k,904) = rxt(k,294)*y(k,56)
         mat(k,668) = .250_r8*rxt(k,357)*y(k,32) + .100_r8*rxt(k,358)*y(k,37) &
                      + .250_r8*rxt(k,360)*y(k,97) + .250_r8*rxt(k,361)*y(k,99)
         mat(k,365) = -(rxt(k,296)*y(k,111))
         mat(k,1368) = -rxt(k,296)*y(k,61)
         mat(k,782) = .250_r8*rxt(k,357)*y(k,134)
         mat(k,1168) = .100_r8*rxt(k,358)*y(k,134)
         mat(k,606) = .200_r8*rxt(k,295)*y(k,111)
         mat(k,694) = .020_r8*rxt(k,350)*y(k,97)
         mat(k,1466) = .020_r8*rxt(k,350)*y(k,81) + .250_r8*rxt(k,360)*y(k,134)
         mat(k,1209) = .250_r8*rxt(k,361)*y(k,134)
         mat(k,1368) = mat(k,1368) + .200_r8*rxt(k,295)*y(k,60)
         mat(k,667) = .250_r8*rxt(k,357)*y(k,32) + .100_r8*rxt(k,358)*y(k,37) &
                      + .250_r8*rxt(k,360)*y(k,97) + .250_r8*rxt(k,361)*y(k,99)
         mat(k,1556) = -((rxt(k,137) + rxt(k,138) + rxt(k,139)) * y(k,139) + rxt(k,140) &
                      *y(k,107) + rxt(k,143)*y(k,108))
         mat(k,1092) = -(rxt(k,137) + rxt(k,138) + rxt(k,139)) * y(k,62)
         mat(k,921) = -rxt(k,140)*y(k,62)
         mat(k,1138) = -rxt(k,143)*y(k,62)
         mat(k,1033) = rxt(k,261)*y(k,111)
         mat(k,1535) = rxt(k,270)*y(k,144)
         mat(k,1328) = rxt(k,181)*y(k,63)
         mat(k,497) = rxt(k,237)*y(k,63)
         mat(k,836) = rxt(k,181)*y(k,42) + rxt(k,237)*y(k,59) + rxt(k,135)*y(k,106) &
                      + rxt(k,144)*y(k,111) + rxt(k,117)*y(k,144)
         mat(k,439) = rxt(k,235)*y(k,144)
         mat(k,871) = rxt(k,212)*y(k,144)
         mat(k,570) = rxt(k,167)*y(k,111)
         mat(k,1462) = rxt(k,135)*y(k,63) + rxt(k,147)*y(k,111)
         mat(k,188) = rxt(k,369)*y(k,111)
         mat(k,1417) = rxt(k,261)*y(k,27) + rxt(k,144)*y(k,63) + rxt(k,167)*y(k,89) &
                      + rxt(k,147)*y(k,106) + rxt(k,369)*y(k,110) + rxt(k,375) &
                      *y(k,120) + rxt(k,380)*y(k,122)
         mat(k,364) = rxt(k,375)*y(k,111)
         mat(k,853) = rxt(k,380)*y(k,111)
         mat(k,1290) = rxt(k,270)*y(k,40) + rxt(k,117)*y(k,63) + rxt(k,235)*y(k,67) &
                      + rxt(k,212)*y(k,71)
         mat(k,826) = -(rxt(k,117)*y(k,144) + rxt(k,135)*y(k,106) + rxt(k,144) &
                      *y(k,111) + rxt(k,181)*y(k,42) + rxt(k,237)*y(k,59))
         mat(k,1271) = -rxt(k,117)*y(k,63)
         mat(k,1442) = -rxt(k,135)*y(k,63)
         mat(k,1397) = -rxt(k,144)*y(k,63)
         mat(k,1309) = -rxt(k,181)*y(k,63)
         mat(k,491) = -rxt(k,237)*y(k,63)
         mat(k,1516) = rxt(k,271)*y(k,144)
         mat(k,1537) = rxt(k,137)*y(k,139)
         mat(k,1073) = rxt(k,137)*y(k,62)
         mat(k,1271) = mat(k,1271) + rxt(k,271)*y(k,40)
         mat(k,41) = -(rxt(k,233)*y(k,144))
         mat(k,1251) = -rxt(k,233)*y(k,64)
         mat(k,311) = -(rxt(k,136)*y(k,106) + rxt(k,145)*y(k,111) + rxt(k,182)*y(k,42))
         mat(k,1425) = -rxt(k,136)*y(k,65)
         mat(k,1362) = -rxt(k,145)*y(k,65)
         mat(k,1298) = -rxt(k,182)*y(k,65)
         mat(k,1362) = mat(k,1362) + 2.000_r8*rxt(k,150)*y(k,111)
         mat(k,1046) = 2.000_r8*rxt(k,151)*y(k,139)
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
         mat(k,130) = rxt(k,382)*y(k,151)
         mat(k,1558) = rxt(k,382)*y(k,124)
         mat(k,433) = -(rxt(k,228)*y(k,106) + rxt(k,229)*y(k,111) + (rxt(k,234) &
                      + rxt(k,235)) * y(k,144))
         mat(k,1430) = -rxt(k,228)*y(k,67)
         mat(k,1373) = -rxt(k,229)*y(k,67)
         mat(k,1269) = -(rxt(k,234) + rxt(k,235)) * y(k,67)
         mat(k,873) = rxt(k,215)*y(k,27) + rxt(k,216)*y(k,139)
         mat(k,1011) = rxt(k,215)*y(k,4)
         mat(k,1055) = rxt(k,216)*y(k,4)
         mat(k,107) = -(rxt(k,251)*y(k,111) + rxt(k,256)*y(k,144))
         mat(k,1336) = -rxt(k,251)*y(k,68)
         mat(k,1261) = -rxt(k,256)*y(k,68)
         mat(k,117) = -(rxt(k,252)*y(k,111) + rxt(k,257)*y(k,144))
         mat(k,1338) = -rxt(k,252)*y(k,69)
         mat(k,1263) = -rxt(k,257)*y(k,69)
         mat(k,142) = -(rxt(k,253)*y(k,111) + rxt(k,258)*y(k,144))
         mat(k,1342) = -rxt(k,253)*y(k,70)
         mat(k,1265) = -rxt(k,258)*y(k,70)
         mat(k,860) = -(rxt(k,199)*y(k,106) + rxt(k,200)*y(k,111) + (rxt(k,211) &
                      + rxt(k,212)) * y(k,144) + (rxt(k,392) + rxt(k,398) + rxt(k,403) &
                      ) * y(k,76) + (rxt(k,397) + rxt(k,402) + rxt(k,407)) * y(k,46) &
                      + (rxt(k,399) + rxt(k,404)) * y(k,75))
         mat(k,1444) = -rxt(k,199)*y(k,71)
         mat(k,1399) = -rxt(k,200)*y(k,71)
         mat(k,1272) = -(rxt(k,211) + rxt(k,212)) * y(k,71)
         mat(k,452) = -(rxt(k,392) + rxt(k,398) + rxt(k,403)) * y(k,71)
         mat(k,574) = -(rxt(k,397) + rxt(k,402) + rxt(k,407)) * y(k,71)
         mat(k,387) = -(rxt(k,399) + rxt(k,404)) * y(k,71)
         mat(k,159) = rxt(k,281)*y(k,42)
         mat(k,247) = rxt(k,242)*y(k,42)
         mat(k,1015) = rxt(k,179)*y(k,42)
         mat(k,326) = rxt(k,244)*y(k,42)
         mat(k,191) = 2.000_r8*rxt(k,247)*y(k,42)
         mat(k,1517) = rxt(k,180)*y(k,42)
         mat(k,222) = rxt(k,249)*y(k,42)
         mat(k,1310) = rxt(k,281)*y(k,13) + rxt(k,242)*y(k,26) + rxt(k,179)*y(k,27) &
                      + rxt(k,244)*y(k,28) + 2.000_r8*rxt(k,247)*y(k,31) + rxt(k,180) &
                      *y(k,40) + rxt(k,249)*y(k,41) + rxt(k,181)*y(k,63) + rxt(k,182) &
                      *y(k,65) + rxt(k,201)*y(k,76) + rxt(k,183)*y(k,139)
         mat(k,955) = rxt(k,198)*y(k,111)
         mat(k,827) = rxt(k,181)*y(k,42)
         mat(k,312) = rxt(k,182)*y(k,42)
         mat(k,452) = mat(k,452) + rxt(k,201)*y(k,42)
         mat(k,1399) = mat(k,1399) + rxt(k,198)*y(k,45)
         mat(k,1074) = rxt(k,183)*y(k,42)
         mat(k,1512) = rxt(k,236)*y(k,59)
         mat(k,489) = rxt(k,236)*y(k,40) + rxt(k,237)*y(k,63) + rxt(k,239)*y(k,73) &
                      + rxt(k,238)*y(k,151)
         mat(k,824) = rxt(k,237)*y(k,59)
         mat(k,1140) = rxt(k,239)*y(k,59)
         mat(k,1560) = rxt(k,238)*y(k,59)
         mat(k,1152) = -(rxt(k,176)*y(k,111) + rxt(k,239)*y(k,59))
         mat(k,1408) = -rxt(k,176)*y(k,73)
         mat(k,492) = -rxt(k,239)*y(k,73)
         mat(k,1024) = rxt(k,259)*y(k,99)
         mat(k,644) = rxt(k,283)*y(k,99)
         mat(k,687) = rxt(k,308)*y(k,99)
         mat(k,577) = (rxt(k,397)+rxt(k,402)+rxt(k,407))*y(k,71)
         mat(k,138) = rxt(k,366)*y(k,99)
         mat(k,865) = (rxt(k,397)+rxt(k,402)+rxt(k,407))*y(k,46)
         mat(k,1000) = rxt(k,175)*y(k,111)
         mat(k,1239) = rxt(k,259)*y(k,27) + rxt(k,283)*y(k,30) + rxt(k,308)*y(k,34) &
                      + rxt(k,366)*y(k,52)
         mat(k,1408) = mat(k,1408) + rxt(k,175)*y(k,98)
         mat(k,232) = -(rxt(k,152)*y(k,111))
         mat(k,1354) = -rxt(k,152)*y(k,74)
         mat(k,976) = rxt(k,173)*y(k,139)
         mat(k,1043) = rxt(k,173)*y(k,98)
         mat(k,386) = -(rxt(k,230)*y(k,106) + (rxt(k,399) + rxt(k,404)) * y(k,71))
         mat(k,1429) = -rxt(k,230)*y(k,75)
         mat(k,857) = -(rxt(k,399) + rxt(k,404)) * y(k,75)
         mat(k,926) = rxt(k,222)*y(k,139)
         mat(k,1050) = rxt(k,222)*y(k,6)
         mat(k,451) = -(rxt(k,201)*y(k,42) + rxt(k,202)*y(k,106) + rxt(k,203)*y(k,111) &
                      + (rxt(k,392) + rxt(k,398) + rxt(k,403)) * y(k,71))
         mat(k,1302) = -rxt(k,201)*y(k,76)
         mat(k,1432) = -rxt(k,202)*y(k,76)
         mat(k,1374) = -rxt(k,203)*y(k,76)
         mat(k,858) = -(rxt(k,392) + rxt(k,398) + rxt(k,403)) * y(k,76)
         mat(k,951) = rxt(k,190)*y(k,139)
         mat(k,572) = rxt(k,195)*y(k,111)
         mat(k,1374) = mat(k,1374) + rxt(k,195)*y(k,46)
         mat(k,1056) = rxt(k,190)*y(k,45)
         mat(k,613) = -(rxt(k,310)*y(k,111))
         mat(k,1384) = -rxt(k,310)*y(k,77)
         mat(k,784) = .220_r8*rxt(k,321)*y(k,84) + .250_r8*rxt(k,357)*y(k,134)
         mat(k,1174) = .230_r8*rxt(k,322)*y(k,84) + .200_r8*rxt(k,315)*y(k,118) &
                      + .100_r8*rxt(k,358)*y(k,134)
         mat(k,697) = .020_r8*rxt(k,350)*y(k,97)
         mat(k,722) = .220_r8*rxt(k,321)*y(k,32) + .230_r8*rxt(k,322)*y(k,37) &
                      + .220_r8*rxt(k,325)*y(k,97) + .220_r8*rxt(k,324)*y(k,99)
         mat(k,303) = .500_r8*rxt(k,336)*y(k,111)
         mat(k,1481) = .020_r8*rxt(k,350)*y(k,81) + .220_r8*rxt(k,325)*y(k,84) &
                      + .250_r8*rxt(k,360)*y(k,134)
         mat(k,1217) = .220_r8*rxt(k,324)*y(k,84) + .250_r8*rxt(k,361)*y(k,134)
         mat(k,1384) = mat(k,1384) + .500_r8*rxt(k,336)*y(k,87) + .500_r8*rxt(k,314) &
                      *y(k,117)
         mat(k,288) = .500_r8*rxt(k,314)*y(k,111)
         mat(k,596) = .200_r8*rxt(k,315)*y(k,37)
         mat(k,669) = .250_r8*rxt(k,357)*y(k,32) + .100_r8*rxt(k,358)*y(k,37) &
                      + .250_r8*rxt(k,360)*y(k,97) + .250_r8*rxt(k,361)*y(k,99)
         mat(k,178) = -(rxt(k,342)*y(k,111))
         mat(k,1346) = -rxt(k,342)*y(k,78)
         mat(k,778) = .400_r8*rxt(k,347)*y(k,81)
         mat(k,1164) = .300_r8*rxt(k,348)*y(k,81)
         mat(k,691) = .400_r8*rxt(k,347)*y(k,32) + .300_r8*rxt(k,348)*y(k,37) &
                      + .330_r8*rxt(k,350)*y(k,97) + .400_r8*rxt(k,351)*y(k,99)
         mat(k,1464) = .330_r8*rxt(k,350)*y(k,81)
         mat(k,1206) = .400_r8*rxt(k,351)*y(k,81) + rxt(k,355)*y(k,112)
         mat(k,1346) = mat(k,1346) + rxt(k,356)*y(k,112)
         mat(k,583) = rxt(k,355)*y(k,99) + rxt(k,356)*y(k,111)
         mat(k,473) = -(rxt(k,343)*y(k,99) + rxt(k,352)*y(k,108) + rxt(k,353)*y(k,111))
         mat(k,1212) = -rxt(k,343)*y(k,79)
         mat(k,1102) = -rxt(k,352)*y(k,79)
         mat(k,1376) = -rxt(k,353)*y(k,79)
         mat(k,411) = -(rxt(k,344)*y(k,139) + rxt(k,345)*y(k,97) + rxt(k,346)*y(k,99))
         mat(k,1053) = -rxt(k,344)*y(k,80)
         mat(k,1471) = -rxt(k,345)*y(k,80)
         mat(k,1211) = -rxt(k,346)*y(k,80)
         mat(k,472) = rxt(k,343)*y(k,99)
         mat(k,1211) = mat(k,1211) + rxt(k,343)*y(k,79)
         mat(k,702) = -(rxt(k,347)*y(k,32) + rxt(k,348)*y(k,37) + rxt(k,349)*y(k,139) &
                      + rxt(k,350)*y(k,97) + rxt(k,351)*y(k,99))
         mat(k,789) = -rxt(k,347)*y(k,81)
         mat(k,1180) = -rxt(k,348)*y(k,81)
         mat(k,1067) = -rxt(k,349)*y(k,81)
         mat(k,1487) = -rxt(k,350)*y(k,81)
         mat(k,1224) = -rxt(k,351)*y(k,81)
         mat(k,478) = rxt(k,353)*y(k,111)
         mat(k,269) = .200_r8*rxt(k,354)*y(k,111)
         mat(k,1224) = mat(k,1224) + 1.700_r8*rxt(k,363)*y(k,133)
         mat(k,1391) = rxt(k,353)*y(k,79) + .200_r8*rxt(k,354)*y(k,82) &
                      + 1.640_r8*rxt(k,365)*y(k,133)
         mat(k,346) = 1.700_r8*rxt(k,363)*y(k,99) + 1.640_r8*rxt(k,365)*y(k,111)
         mat(k,266) = -(rxt(k,354)*y(k,111))
         mat(k,1357) = -rxt(k,354)*y(k,82)
         mat(k,692) = rxt(k,349)*y(k,139)
         mat(k,1044) = rxt(k,349)*y(k,81)
         mat(k,654) = -(rxt(k,327)*y(k,108) + rxt(k,328)*y(k,111))
         mat(k,1111) = -rxt(k,327)*y(k,83)
         mat(k,1388) = -rxt(k,328)*y(k,83)
         mat(k,786) = .250_r8*rxt(k,347)*y(k,81)
         mat(k,1177) = .190_r8*rxt(k,348)*y(k,81)
         mat(k,477) = .300_r8*rxt(k,352)*y(k,108)
         mat(k,413) = .167_r8*rxt(k,345)*y(k,97) + .167_r8*rxt(k,346)*y(k,99) &
                      + .167_r8*rxt(k,344)*y(k,139)
         mat(k,699) = .250_r8*rxt(k,347)*y(k,32) + .190_r8*rxt(k,348)*y(k,37) &
                      + .230_r8*rxt(k,350)*y(k,97) + .250_r8*rxt(k,351)*y(k,99)
         mat(k,1484) = .167_r8*rxt(k,345)*y(k,80) + .230_r8*rxt(k,350)*y(k,81)
         mat(k,1221) = .167_r8*rxt(k,346)*y(k,80) + .250_r8*rxt(k,351)*y(k,81)
         mat(k,1111) = mat(k,1111) + .300_r8*rxt(k,352)*y(k,79) + 1.122_r8*rxt(k,364) &
                      *y(k,133)
         mat(k,345) = 1.122_r8*rxt(k,364)*y(k,108)
         mat(k,1064) = .167_r8*rxt(k,344)*y(k,80)
         mat(k,726) = -(rxt(k,321)*y(k,32) + rxt(k,322)*y(k,37) + rxt(k,323)*y(k,139) &
                      + rxt(k,324)*y(k,99) + (rxt(k,325) + rxt(k,326)) * y(k,97))
         mat(k,790) = -rxt(k,321)*y(k,84)
         mat(k,1181) = -rxt(k,322)*y(k,84)
         mat(k,1068) = -rxt(k,323)*y(k,84)
         mat(k,1225) = -rxt(k,324)*y(k,84)
         mat(k,1488) = -(rxt(k,325) + rxt(k,326)) * y(k,84)
         mat(k,656) = .500_r8*rxt(k,328)*y(k,111)
         mat(k,164) = .200_r8*rxt(k,329)*y(k,111)
         mat(k,743) = rxt(k,338)*y(k,111)
         mat(k,1392) = .500_r8*rxt(k,328)*y(k,83) + .200_r8*rxt(k,329)*y(k,85) &
                      + rxt(k,338)*y(k,88)
         mat(k,163) = -(rxt(k,329)*y(k,111))
         mat(k,1344) = -rxt(k,329)*y(k,85)
         mat(k,718) = rxt(k,323)*y(k,139)
         mat(k,1037) = rxt(k,323)*y(k,84)
         mat(k,765) = -(rxt(k,330)*y(k,32) + rxt(k,331)*y(k,37) + rxt(k,332)*y(k,139) &
                      + 4._r8*rxt(k,333)*y(k,86) + rxt(k,334)*y(k,97) + rxt(k,335) &
                      *y(k,99) + rxt(k,339)*y(k,98))
         mat(k,792) = -rxt(k,330)*y(k,86)
         mat(k,1183) = -rxt(k,331)*y(k,86)
         mat(k,1070) = -rxt(k,332)*y(k,86)
         mat(k,1490) = -rxt(k,334)*y(k,86)
         mat(k,1227) = -rxt(k,335)*y(k,86)
         mat(k,987) = -rxt(k,339)*y(k,86)
         mat(k,657) = .500_r8*rxt(k,328)*y(k,111)
         mat(k,165) = .500_r8*rxt(k,329)*y(k,111)
         mat(k,1394) = .500_r8*rxt(k,328)*y(k,83) + .500_r8*rxt(k,329)*y(k,85)
         mat(k,302) = -(rxt(k,336)*y(k,111))
         mat(k,1361) = -rxt(k,336)*y(k,87)
         mat(k,760) = rxt(k,339)*y(k,98)
         mat(k,979) = rxt(k,339)*y(k,86)
         mat(k,744) = -(rxt(k,337)*y(k,108) + rxt(k,338)*y(k,111))
         mat(k,1116) = -rxt(k,337)*y(k,88)
         mat(k,1393) = -rxt(k,338)*y(k,88)
         mat(k,791) = .350_r8*rxt(k,347)*y(k,81)
         mat(k,1182) = .260_r8*rxt(k,348)*y(k,81)
         mat(k,479) = .200_r8*rxt(k,352)*y(k,108)
         mat(k,414) = .039_r8*rxt(k,345)*y(k,97) + .039_r8*rxt(k,346)*y(k,99) &
                      + .039_r8*rxt(k,344)*y(k,139)
         mat(k,704) = .350_r8*rxt(k,347)*y(k,32) + .260_r8*rxt(k,348)*y(k,37) &
                      + .320_r8*rxt(k,350)*y(k,97) + .350_r8*rxt(k,351)*y(k,99)
         mat(k,1489) = .039_r8*rxt(k,345)*y(k,80) + .320_r8*rxt(k,350)*y(k,81)
         mat(k,1226) = .039_r8*rxt(k,346)*y(k,80) + .350_r8*rxt(k,351)*y(k,81)
         mat(k,1116) = mat(k,1116) + .200_r8*rxt(k,352)*y(k,79) + .442_r8*rxt(k,364) &
                      *y(k,133)
         mat(k,347) = .442_r8*rxt(k,364)*y(k,108)
         mat(k,1069) = .039_r8*rxt(k,344)*y(k,80)
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
         mat(k,561) = -(rxt(k,155)*y(k,97) + (rxt(k,156) + rxt(k,157) + rxt(k,158) &
                      ) * y(k,98) + rxt(k,159)*y(k,107) + rxt(k,167)*y(k,111) &
                      + rxt(k,420)*y(k,147))
         mat(k,1477) = -rxt(k,155)*y(k,89)
         mat(k,982) = -(rxt(k,156) + rxt(k,157) + rxt(k,158)) * y(k,89)
         mat(k,903) = -rxt(k,159)*y(k,89)
         mat(k,1379) = -rxt(k,167)*y(k,89)
         mat(k,514) = -rxt(k,420)*y(k,89)
         mat(k,1438) = rxt(k,153)*y(k,140) + rxt(k,417)*y(k,143)
         mat(k,903) = mat(k,903) + rxt(k,418)*y(k,143)
         mat(k,528) = 1.100_r8*rxt(k,413)*y(k,141) + .200_r8*rxt(k,411)*y(k,142)
         mat(k,540) = rxt(k,153)*y(k,106)
         mat(k,340) = 1.100_r8*rxt(k,413)*y(k,138)
         mat(k,373) = .200_r8*rxt(k,411)*y(k,138)
         mat(k,447) = rxt(k,417)*y(k,106) + rxt(k,418)*y(k,107)
         mat(k,122) = -((rxt(k,171) + rxt(k,172)) * y(k,144))
         mat(k,1264) = -(rxt(k,171) + rxt(k,172)) * y(k,90)
         mat(k,554) = rxt(k,156)*y(k,98)
         mat(k,974) = rxt(k,156)*y(k,89)
         mat(k,975) = rxt(k,174)*y(k,99)
         mat(k,1205) = rxt(k,174)*y(k,98)
         mat(k,35) = -(rxt(k,383)*y(k,111))
         mat(k,1331) = -rxt(k,383)*y(k,95)
         mat(k,1508) = -(rxt(k,155)*y(k,89) + rxt(k,164)*y(k,99) + rxt(k,168)*y(k,139) &
                      + rxt(k,169)*y(k,108) + rxt(k,170)*y(k,106) + rxt(k,191)*y(k,45) &
                      + rxt(k,223)*y(k,6) + rxt(k,265)*y(k,37) + rxt(k,278)*y(k,10) &
                      + rxt(k,288)*y(k,32) + rxt(k,292)*y(k,57) + rxt(k,305)*y(k,15) &
                      + rxt(k,313)*y(k,114) + rxt(k,317)*y(k,118) + (rxt(k,325) &
                      + rxt(k,326)) * y(k,84) + rxt(k,334)*y(k,86) + rxt(k,345) &
                      *y(k,80) + rxt(k,350)*y(k,81) + rxt(k,360)*y(k,134) + rxt(k,422) &
                      *y(k,147))
         mat(k,569) = -rxt(k,155)*y(k,97)
         mat(k,1246) = -rxt(k,164)*y(k,97)
         mat(k,1090) = -rxt(k,168)*y(k,97)
         mat(k,1136) = -rxt(k,169)*y(k,97)
         mat(k,1460) = -rxt(k,170)*y(k,97)
         mat(k,971) = -rxt(k,191)*y(k,97)
         mat(k,944) = -rxt(k,223)*y(k,97)
         mat(k,1200) = -rxt(k,265)*y(k,97)
         mat(k,432) = -rxt(k,278)*y(k,97)
         mat(k,804) = -rxt(k,288)*y(k,97)
         mat(k,384) = -rxt(k,292)*y(k,97)
         mat(k,470) = -rxt(k,305)*y(k,97)
         mat(k,409) = -rxt(k,313)*y(k,97)
         mat(k,604) = -rxt(k,317)*y(k,97)
         mat(k,737) = -(rxt(k,325) + rxt(k,326)) * y(k,97)
         mat(k,776) = -rxt(k,334)*y(k,97)
         mat(k,420) = -rxt(k,345)*y(k,97)
         mat(k,716) = -rxt(k,350)*y(k,97)
         mat(k,681) = -rxt(k,360)*y(k,97)
         mat(k,518) = -rxt(k,422)*y(k,97)
         mat(k,569) = mat(k,569) + 2.000_r8*rxt(k,157)*y(k,98) + rxt(k,159)*y(k,107) &
                      + rxt(k,167)*y(k,111)
         mat(k,125) = 2.000_r8*rxt(k,171)*y(k,144)
         mat(k,1007) = 2.000_r8*rxt(k,157)*y(k,89) + rxt(k,160)*y(k,106) + rxt(k,376) &
                      *y(k,122)
         mat(k,1460) = mat(k,1460) + rxt(k,160)*y(k,98)
         mat(k,920) = rxt(k,159)*y(k,89) + rxt(k,154)*y(k,140)
         mat(k,1415) = rxt(k,167)*y(k,89)
         mat(k,852) = rxt(k,376)*y(k,98)
         mat(k,546) = rxt(k,154)*y(k,107)
         mat(k,1288) = 2.000_r8*rxt(k,171)*y(k,90)
         mat(k,996) = -((rxt(k,156) + rxt(k,157) + rxt(k,158)) * y(k,89) + (rxt(k,160) &
                      + rxt(k,162)) * y(k,106) + rxt(k,161)*y(k,108) + rxt(k,173) &
                      *y(k,139) + rxt(k,174)*y(k,99) + rxt(k,175)*y(k,111) + rxt(k,193) &
                      *y(k,45) + rxt(k,224)*y(k,6) + rxt(k,299)*y(k,32) + rxt(k,339) &
                      *y(k,86) + rxt(k,376)*y(k,122))
         mat(k,565) = -(rxt(k,156) + rxt(k,157) + rxt(k,158)) * y(k,98)
         mat(k,1449) = -(rxt(k,160) + rxt(k,162)) * y(k,98)
         mat(k,1125) = -rxt(k,161)*y(k,98)
         mat(k,1079) = -rxt(k,173)*y(k,98)
         mat(k,1235) = -rxt(k,174)*y(k,98)
         mat(k,1404) = -rxt(k,175)*y(k,98)
         mat(k,960) = -rxt(k,193)*y(k,98)
         mat(k,934) = -rxt(k,224)*y(k,98)
         mat(k,796) = -rxt(k,299)*y(k,98)
         mat(k,768) = -rxt(k,339)*y(k,98)
         mat(k,846) = -rxt(k,376)*y(k,98)
         mat(k,934) = mat(k,934) + rxt(k,223)*y(k,97)
         mat(k,427) = rxt(k,278)*y(k,97)
         mat(k,465) = rxt(k,305)*y(k,97)
         mat(k,796) = mat(k,796) + rxt(k,288)*y(k,97)
         mat(k,1189) = rxt(k,265)*y(k,97)
         mat(k,960) = mat(k,960) + rxt(k,191)*y(k,97)
         mat(k,380) = rxt(k,292)*y(k,97)
         mat(k,234) = rxt(k,152)*y(k,111)
         mat(k,415) = 1.206_r8*rxt(k,345)*y(k,97) + 1.206_r8*rxt(k,346)*y(k,99) &
                      + .206_r8*rxt(k,344)*y(k,139)
         mat(k,708) = .920_r8*rxt(k,350)*y(k,97) + rxt(k,351)*y(k,99)
         mat(k,730) = rxt(k,325)*y(k,97) + rxt(k,324)*y(k,99)
         mat(k,768) = mat(k,768) + rxt(k,334)*y(k,97) + rxt(k,335)*y(k,99)
         mat(k,1497) = rxt(k,223)*y(k,6) + rxt(k,278)*y(k,10) + rxt(k,305)*y(k,15) &
                      + rxt(k,288)*y(k,32) + rxt(k,265)*y(k,37) + rxt(k,191)*y(k,45) &
                      + rxt(k,292)*y(k,57) + 1.206_r8*rxt(k,345)*y(k,80) &
                      + .920_r8*rxt(k,350)*y(k,81) + rxt(k,325)*y(k,84) + rxt(k,334) &
                      *y(k,86) + 2.000_r8*rxt(k,164)*y(k,99) + rxt(k,170)*y(k,106) &
                      + rxt(k,169)*y(k,108) + rxt(k,313)*y(k,114) + rxt(k,317) &
                      *y(k,118) + rxt(k,360)*y(k,134) + rxt(k,168)*y(k,139)
         mat(k,1235) = mat(k,1235) + 1.206_r8*rxt(k,346)*y(k,80) + rxt(k,351)*y(k,81) &
                      + rxt(k,324)*y(k,84) + rxt(k,335)*y(k,86) + 2.000_r8*rxt(k,164) &
                      *y(k,97) + rxt(k,165)*y(k,106) + rxt(k,166)*y(k,111) &
                      + rxt(k,355)*y(k,112) + rxt(k,363)*y(k,133) + rxt(k,361) &
                      *y(k,134) + rxt(k,163)*y(k,139)
         mat(k,200) = rxt(k,311)*y(k,111)
         mat(k,1449) = mat(k,1449) + rxt(k,170)*y(k,97) + rxt(k,165)*y(k,99)
         mat(k,1125) = mat(k,1125) + rxt(k,169)*y(k,97)
         mat(k,1404) = mat(k,1404) + rxt(k,152)*y(k,74) + rxt(k,166)*y(k,99) &
                      + rxt(k,311)*y(k,100) + .400_r8*rxt(k,356)*y(k,112)
         mat(k,587) = rxt(k,355)*y(k,99) + .400_r8*rxt(k,356)*y(k,111)
         mat(k,405) = rxt(k,313)*y(k,97)
         mat(k,599) = rxt(k,317)*y(k,97)
         mat(k,348) = rxt(k,363)*y(k,99)
         mat(k,675) = rxt(k,360)*y(k,97) + rxt(k,361)*y(k,99)
         mat(k,1079) = mat(k,1079) + .206_r8*rxt(k,344)*y(k,80) + rxt(k,168)*y(k,97) &
                      + rxt(k,163)*y(k,99)
         mat(k,1241) = -(rxt(k,163)*y(k,139) + rxt(k,164)*y(k,97) + rxt(k,165) &
                      *y(k,106) + rxt(k,166)*y(k,111) + rxt(k,174)*y(k,98) + rxt(k,259) &
                      *y(k,27) + rxt(k,283)*y(k,30) + rxt(k,301)*y(k,14) + rxt(k,308) &
                      *y(k,34) + rxt(k,324)*y(k,84) + rxt(k,335)*y(k,86) + rxt(k,343) &
                      *y(k,79) + rxt(k,346)*y(k,80) + rxt(k,351)*y(k,81) + rxt(k,355) &
                      *y(k,112) + rxt(k,361)*y(k,134) + rxt(k,363)*y(k,133) + rxt(k,366) &
                      *y(k,52))
         mat(k,1085) = -rxt(k,163)*y(k,99)
         mat(k,1503) = -rxt(k,164)*y(k,99)
         mat(k,1455) = -rxt(k,165)*y(k,99)
         mat(k,1410) = -rxt(k,166)*y(k,99)
         mat(k,1002) = -rxt(k,174)*y(k,99)
         mat(k,1026) = -rxt(k,259)*y(k,99)
         mat(k,646) = -rxt(k,283)*y(k,99)
         mat(k,635) = -rxt(k,301)*y(k,99)
         mat(k,688) = -rxt(k,308)*y(k,99)
         mat(k,735) = -rxt(k,324)*y(k,99)
         mat(k,774) = -rxt(k,335)*y(k,99)
         mat(k,486) = -rxt(k,343)*y(k,99)
         mat(k,418) = -rxt(k,346)*y(k,99)
         mat(k,714) = -rxt(k,351)*y(k,99)
         mat(k,591) = -rxt(k,355)*y(k,99)
         mat(k,679) = -rxt(k,361)*y(k,99)
         mat(k,350) = -rxt(k,363)*y(k,99)
         mat(k,139) = -rxt(k,366)*y(k,99)
         mat(k,300) = rxt(k,225)*y(k,106)
         mat(k,1321) = rxt(k,192)*y(k,46)
         mat(k,578) = rxt(k,192)*y(k,42) + rxt(k,194)*y(k,106) + rxt(k,195)*y(k,111)
         mat(k,494) = rxt(k,239)*y(k,73)
         mat(k,1154) = rxt(k,239)*y(k,59) + rxt(k,176)*y(k,111)
         mat(k,309) = .500_r8*rxt(k,336)*y(k,111)
         mat(k,1002) = mat(k,1002) + rxt(k,162)*y(k,106) + rxt(k,161)*y(k,108)
         mat(k,1455) = mat(k,1455) + rxt(k,225)*y(k,7) + rxt(k,194)*y(k,46) &
                      + rxt(k,162)*y(k,98)
         mat(k,1131) = rxt(k,161)*y(k,98)
         mat(k,1410) = mat(k,1410) + rxt(k,195)*y(k,46) + rxt(k,176)*y(k,73) &
                      + .500_r8*rxt(k,336)*y(k,87) + rxt(k,297)*y(k,113)
         mat(k,280) = rxt(k,297)*y(k,111)
         mat(k,197) = -(rxt(k,311)*y(k,111))
         mat(k,1349) = -rxt(k,311)*y(k,100)
         mat(k,619) = rxt(k,301)*y(k,99)
         mat(k,1207) = rxt(k,301)*y(k,14)
         mat(k,1459) = -(rxt(k,132)*y(k,108) + 4._r8*rxt(k,133)*y(k,106) + rxt(k,134) &
                      *y(k,107) + rxt(k,135)*y(k,63) + rxt(k,136)*y(k,65) + rxt(k,141) &
                      *y(k,139) + rxt(k,147)*y(k,111) + (rxt(k,160) + rxt(k,162) &
                      ) * y(k,98) + rxt(k,165)*y(k,99) + rxt(k,170)*y(k,97) + rxt(k,194) &
                      *y(k,46) + rxt(k,196)*y(k,45) + rxt(k,199)*y(k,71) + rxt(k,202) &
                      *y(k,76) + rxt(k,225)*y(k,7) + rxt(k,226)*y(k,6) + rxt(k,228) &
                      *y(k,67) + rxt(k,230)*y(k,75) + rxt(k,260)*y(k,27) + rxt(k,368) &
                      *y(k,110) + (rxt(k,415) + rxt(k,416)) * y(k,141) + rxt(k,417) &
                      *y(k,143))
         mat(k,1135) = -rxt(k,132)*y(k,106)
         mat(k,919) = -rxt(k,134)*y(k,106)
         mat(k,834) = -rxt(k,135)*y(k,106)
         mat(k,316) = -rxt(k,136)*y(k,106)
         mat(k,1089) = -rxt(k,141)*y(k,106)
         mat(k,1414) = -rxt(k,147)*y(k,106)
         mat(k,1006) = -(rxt(k,160) + rxt(k,162)) * y(k,106)
         mat(k,1245) = -rxt(k,165)*y(k,106)
         mat(k,1507) = -rxt(k,170)*y(k,106)
         mat(k,581) = -rxt(k,194)*y(k,106)
         mat(k,970) = -rxt(k,196)*y(k,106)
         mat(k,870) = -rxt(k,199)*y(k,106)
         mat(k,456) = -rxt(k,202)*y(k,106)
         mat(k,301) = -rxt(k,225)*y(k,106)
         mat(k,943) = -rxt(k,226)*y(k,106)
         mat(k,438) = -rxt(k,228)*y(k,106)
         mat(k,392) = -rxt(k,230)*y(k,106)
         mat(k,1030) = -rxt(k,260)*y(k,106)
         mat(k,187) = -rxt(k,368)*y(k,106)
         mat(k,342) = -(rxt(k,415) + rxt(k,416)) * y(k,106)
         mat(k,449) = -rxt(k,417)*y(k,106)
         mat(k,1553) = rxt(k,139)*y(k,139)
         mat(k,568) = rxt(k,155)*y(k,97) + rxt(k,156)*y(k,98) + rxt(k,159)*y(k,107) &
                      + rxt(k,420)*y(k,147)
         mat(k,1507) = mat(k,1507) + rxt(k,155)*y(k,89)
         mat(k,1006) = mat(k,1006) + rxt(k,156)*y(k,89)
         mat(k,919) = mat(k,919) + rxt(k,159)*y(k,89) + rxt(k,370)*y(k,120) &
                      + rxt(k,377)*y(k,122) + rxt(k,419)*y(k,143) + (rxt(k,120) &
                       +rxt(k,121))*y(k,144) + rxt(k,426)*y(k,148) + rxt(k,430) &
                      *y(k,149)
         mat(k,1135) = mat(k,1135) + .765_r8*rxt(k,364)*y(k,133) + 2.000_r8*rxt(k,123) &
                      *y(k,144)
         mat(k,1414) = mat(k,1414) + 2.000_r8*rxt(k,149)*y(k,111)
         mat(k,363) = rxt(k,370)*y(k,107)
         mat(k,851) = rxt(k,377)*y(k,107)
         mat(k,352) = .765_r8*rxt(k,364)*y(k,108)
         mat(k,533) = rxt(k,411)*y(k,142) + 1.150_r8*rxt(k,412)*y(k,147)
         mat(k,1089) = mat(k,1089) + rxt(k,139)*y(k,62)
         mat(k,545) = rxt(k,425)*y(k,148)
         mat(k,374) = rxt(k,411)*y(k,138)
         mat(k,449) = mat(k,449) + rxt(k,419)*y(k,107)
         mat(k,1287) = (rxt(k,120)+rxt(k,121))*y(k,107) + 2.000_r8*rxt(k,123)*y(k,108)
         mat(k,517) = rxt(k,420)*y(k,89) + 1.150_r8*rxt(k,412)*y(k,138)
         mat(k,509) = rxt(k,426)*y(k,107) + rxt(k,425)*y(k,140)
         mat(k,265) = rxt(k,430)*y(k,107)
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
         mat(k,909) = -(rxt(k,120)*y(k,144) + rxt(k,126)*y(k,145) + rxt(k,134) &
                      *y(k,106) + rxt(k,140)*y(k,62) + rxt(k,154)*y(k,140) + rxt(k,159) &
                      *y(k,89) + rxt(k,294)*y(k,56) + rxt(k,370)*y(k,120) + rxt(k,377) &
                      *y(k,122) + rxt(k,414)*y(k,141) + (rxt(k,418) + rxt(k,419) &
                      ) * y(k,143) + rxt(k,426)*y(k,148) + rxt(k,430)*y(k,149))
         mat(k,1274) = -rxt(k,120)*y(k,107)
         mat(k,94) = -rxt(k,126)*y(k,107)
         mat(k,1446) = -rxt(k,134)*y(k,107)
         mat(k,1540) = -rxt(k,140)*y(k,107)
         mat(k,543) = -rxt(k,154)*y(k,107)
         mat(k,564) = -rxt(k,159)*y(k,107)
         mat(k,229) = -rxt(k,294)*y(k,107)
         mat(k,360) = -rxt(k,370)*y(k,107)
         mat(k,843) = -rxt(k,377)*y(k,107)
         mat(k,341) = -rxt(k,414)*y(k,107)
         mat(k,448) = -(rxt(k,418) + rxt(k,419)) * y(k,107)
         mat(k,508) = -rxt(k,426)*y(k,107)
         mat(k,264) = -rxt(k,430)*y(k,107)
         mat(k,877) = rxt(k,217)*y(k,108) + rxt(k,216)*y(k,139)
         mat(k,931) = 2.000_r8*rxt(k,218)*y(k,6) + (rxt(k,220)+rxt(k,221))*y(k,45) &
                      + rxt(k,226)*y(k,106) + rxt(k,222)*y(k,139)
         mat(k,426) = rxt(k,277)*y(k,139)
         mat(k,464) = rxt(k,304)*y(k,139)
         mat(k,1186) = rxt(k,264)*y(k,139)
         mat(k,1312) = rxt(k,185)*y(k,108) + rxt(k,183)*y(k,139)
         mat(k,957) = (rxt(k,220)+rxt(k,221))*y(k,6) + (2.000_r8*rxt(k,187) &
                       +2.000_r8*rxt(k,188))*y(k,45) + rxt(k,196)*y(k,106) &
                      + rxt(k,198)*y(k,111) + rxt(k,190)*y(k,139)
         mat(k,1540) = mat(k,1540) + rxt(k,143)*y(k,108) + rxt(k,137)*y(k,139)
         mat(k,233) = rxt(k,152)*y(k,111)
         mat(k,564) = mat(k,564) + rxt(k,158)*y(k,98)
         mat(k,123) = rxt(k,172)*y(k,144)
         mat(k,1494) = rxt(k,169)*y(k,108) + rxt(k,422)*y(k,147)
         mat(k,993) = rxt(k,158)*y(k,89) + rxt(k,160)*y(k,106) + rxt(k,161)*y(k,108)
         mat(k,1232) = rxt(k,165)*y(k,106) + rxt(k,163)*y(k,139)
         mat(k,1446) = mat(k,1446) + rxt(k,226)*y(k,6) + rxt(k,196)*y(k,45) &
                      + rxt(k,160)*y(k,98) + rxt(k,165)*y(k,99) + 2.000_r8*rxt(k,133) &
                      *y(k,106) + 2.000_r8*rxt(k,132)*y(k,108) + rxt(k,147)*y(k,111) &
                      + rxt(k,141)*y(k,139) + rxt(k,125)*y(k,145)
         mat(k,909) = mat(k,909) + 2.000_r8*rxt(k,126)*y(k,145)
         mat(k,1122) = rxt(k,217)*y(k,4) + rxt(k,185)*y(k,42) + rxt(k,143)*y(k,62) &
                      + rxt(k,169)*y(k,97) + rxt(k,161)*y(k,98) + 2.000_r8*rxt(k,132) &
                      *y(k,106) + rxt(k,148)*y(k,111) + rxt(k,372)*y(k,120) &
                      + rxt(k,378)*y(k,122) + 2.000_r8*rxt(k,142)*y(k,139) + ( &
                      + 2.000_r8*rxt(k,122)+rxt(k,123))*y(k,144)
         mat(k,1401) = rxt(k,198)*y(k,45) + rxt(k,152)*y(k,74) + rxt(k,147)*y(k,106) &
                      + rxt(k,148)*y(k,108) + rxt(k,146)*y(k,139)
         mat(k,404) = rxt(k,312)*y(k,139)
         mat(k,360) = mat(k,360) + rxt(k,372)*y(k,108)
         mat(k,843) = mat(k,843) + rxt(k,378)*y(k,108)
         mat(k,1076) = rxt(k,216)*y(k,4) + rxt(k,222)*y(k,6) + rxt(k,277)*y(k,10) &
                      + rxt(k,304)*y(k,15) + rxt(k,264)*y(k,37) + rxt(k,183)*y(k,42) &
                      + rxt(k,190)*y(k,45) + rxt(k,137)*y(k,62) + rxt(k,163)*y(k,99) &
                      + rxt(k,141)*y(k,106) + 2.000_r8*rxt(k,142)*y(k,108) &
                      + rxt(k,146)*y(k,111) + rxt(k,312)*y(k,114) &
                      + 2.000_r8*rxt(k,151)*y(k,139)
         mat(k,1274) = mat(k,1274) + rxt(k,172)*y(k,90) + (2.000_r8*rxt(k,122) &
                       +rxt(k,123))*y(k,108)
         mat(k,94) = mat(k,94) + rxt(k,125)*y(k,106) + 2.000_r8*rxt(k,126)*y(k,107)
         mat(k,515) = rxt(k,422)*y(k,97)
         mat(k,1128) = -((rxt(k,122) + rxt(k,123)) * y(k,144) + rxt(k,132)*y(k,106) &
                      + rxt(k,142)*y(k,139) + rxt(k,143)*y(k,62) + rxt(k,148)*y(k,111) &
                      + rxt(k,161)*y(k,98) + rxt(k,169)*y(k,97) + rxt(k,185)*y(k,42) &
                      + rxt(k,217)*y(k,4) + rxt(k,274)*y(k,9) + rxt(k,302)*y(k,14) &
                      + rxt(k,327)*y(k,83) + rxt(k,337)*y(k,88) + rxt(k,352)*y(k,79) &
                      + rxt(k,364)*y(k,133) + rxt(k,372)*y(k,120) + rxt(k,378) &
                      *y(k,122))
         mat(k,1280) = -(rxt(k,122) + rxt(k,123)) * y(k,108)
         mat(k,1452) = -rxt(k,132)*y(k,108)
         mat(k,1082) = -rxt(k,142)*y(k,108)
         mat(k,1546) = -rxt(k,143)*y(k,108)
         mat(k,1407) = -rxt(k,148)*y(k,108)
         mat(k,999) = -rxt(k,161)*y(k,108)
         mat(k,1500) = -rxt(k,169)*y(k,108)
         mat(k,1318) = -rxt(k,185)*y(k,108)
         mat(k,882) = -rxt(k,217)*y(k,108)
         mat(k,244) = -rxt(k,274)*y(k,108)
         mat(k,633) = -rxt(k,302)*y(k,108)
         mat(k,662) = -rxt(k,327)*y(k,108)
         mat(k,752) = -rxt(k,337)*y(k,108)
         mat(k,484) = -rxt(k,352)*y(k,108)
         mat(k,349) = -rxt(k,364)*y(k,108)
         mat(k,361) = -rxt(k,372)*y(k,108)
         mat(k,848) = -rxt(k,378)*y(k,108)
         mat(k,799) = .150_r8*rxt(k,287)*y(k,139)
         mat(k,771) = .150_r8*rxt(k,332)*y(k,139)
         mat(k,1452) = mat(k,1452) + rxt(k,134)*y(k,107)
         mat(k,915) = rxt(k,134)*y(k,106)
         mat(k,1082) = mat(k,1082) + .150_r8*rxt(k,287)*y(k,32) + .150_r8*rxt(k,332) &
                      *y(k,86)
         mat(k,173) = -(rxt(k,379)*y(k,122))
         mat(k,838) = -rxt(k,379)*y(k,109)
         mat(k,924) = rxt(k,219)*y(k,45)
         mat(k,950) = rxt(k,219)*y(k,6) + 2.000_r8*rxt(k,189)*y(k,45)
         mat(k,181) = -(rxt(k,368)*y(k,106) + rxt(k,369)*y(k,111))
         mat(k,1421) = -rxt(k,368)*y(k,110)
         mat(k,1347) = -rxt(k,369)*y(k,110)
         mat(k,1413) = -(rxt(k,144)*y(k,63) + rxt(k,145)*y(k,65) + rxt(k,146)*y(k,139) &
                      + rxt(k,147)*y(k,106) + rxt(k,148)*y(k,108) + (4._r8*rxt(k,149) &
                      + 4._r8*rxt(k,150)) * y(k,111) + rxt(k,152)*y(k,74) + rxt(k,166) &
                      *y(k,99) + rxt(k,167)*y(k,89) + rxt(k,175)*y(k,98) + rxt(k,176) &
                      *y(k,73) + rxt(k,195)*y(k,46) + (rxt(k,197) + rxt(k,198) &
                      ) * y(k,45) + rxt(k,200)*y(k,71) + rxt(k,203)*y(k,76) + rxt(k,227) &
                      *y(k,6) + rxt(k,229)*y(k,67) + rxt(k,243)*y(k,26) + rxt(k,245) &
                      *y(k,28) + rxt(k,246)*y(k,29) + rxt(k,248)*y(k,31) + rxt(k,250) &
                      *y(k,41) + rxt(k,251)*y(k,68) + rxt(k,252)*y(k,69) + rxt(k,253) &
                      *y(k,70) + rxt(k,261)*y(k,27) + rxt(k,266)*y(k,38) + rxt(k,267) &
                      *y(k,39) + rxt(k,268)*y(k,40) + rxt(k,272)*y(k,48) + rxt(k,279) &
                      *y(k,11) + rxt(k,280)*y(k,12) + rxt(k,282)*y(k,13) + rxt(k,284) &
                      *y(k,30) + rxt(k,289)*y(k,35) + rxt(k,290)*y(k,36) + rxt(k,295) &
                      *y(k,60) + rxt(k,296)*y(k,61) + rxt(k,297)*y(k,113) + rxt(k,298) &
                      *y(k,9) + rxt(k,306)*y(k,16) + rxt(k,307)*y(k,17) + rxt(k,309) &
                      *y(k,34) + rxt(k,310)*y(k,77) + rxt(k,311)*y(k,100) + rxt(k,314) &
                      *y(k,117) + rxt(k,318)*y(k,119) + rxt(k,319)*y(k,14) + rxt(k,320) &
                      *y(k,33) + rxt(k,328)*y(k,83) + rxt(k,329)*y(k,85) + rxt(k,336) &
                      *y(k,87) + rxt(k,338)*y(k,88) + rxt(k,341)*y(k,3) + rxt(k,342) &
                      *y(k,78) + rxt(k,353)*y(k,79) + rxt(k,354)*y(k,82) + rxt(k,356) &
                      *y(k,112) + rxt(k,362)*y(k,135) + rxt(k,365)*y(k,133) + (rxt(k,367) &
                      + rxt(k,381)) * y(k,52) + rxt(k,369)*y(k,110) + rxt(k,371) &
                      *y(k,123) + rxt(k,375)*y(k,120) + rxt(k,380)*y(k,122) + rxt(k,383) &
                      *y(k,95))
         mat(k,833) = -rxt(k,144)*y(k,111)
         mat(k,315) = -rxt(k,145)*y(k,111)
         mat(k,1088) = -rxt(k,146)*y(k,111)
         mat(k,1458) = -rxt(k,147)*y(k,111)
         mat(k,1134) = -rxt(k,148)*y(k,111)
         mat(k,237) = -rxt(k,152)*y(k,111)
         mat(k,1244) = -rxt(k,166)*y(k,111)
         mat(k,567) = -rxt(k,167)*y(k,111)
         mat(k,1005) = -rxt(k,175)*y(k,111)
         mat(k,1157) = -rxt(k,176)*y(k,111)
         mat(k,580) = -rxt(k,195)*y(k,111)
         mat(k,969) = -(rxt(k,197) + rxt(k,198)) * y(k,111)
         mat(k,869) = -rxt(k,200)*y(k,111)
         mat(k,455) = -rxt(k,203)*y(k,111)
         mat(k,942) = -rxt(k,227)*y(k,111)
         mat(k,437) = -rxt(k,229)*y(k,111)
         mat(k,251) = -rxt(k,243)*y(k,111)
         mat(k,332) = -rxt(k,245)*y(k,111)
         mat(k,47) = -rxt(k,246)*y(k,111)
         mat(k,195) = -rxt(k,248)*y(k,111)
         mat(k,226) = -rxt(k,250)*y(k,111)
         mat(k,111) = -rxt(k,251)*y(k,111)
         mat(k,121) = -rxt(k,252)*y(k,111)
         mat(k,146) = -rxt(k,253)*y(k,111)
         mat(k,1029) = -rxt(k,261)*y(k,111)
         mat(k,285) = -rxt(k,266)*y(k,111)
         mat(k,212) = -rxt(k,267)*y(k,111)
         mat(k,1531) = -rxt(k,268)*y(k,111)
         mat(k,652) = -rxt(k,272)*y(k,111)
         mat(k,129) = -rxt(k,279)*y(k,111)
         mat(k,172) = -rxt(k,280)*y(k,111)
         mat(k,161) = -rxt(k,282)*y(k,111)
         mat(k,647) = -rxt(k,284)*y(k,111)
         mat(k,356) = -rxt(k,289)*y(k,111)
         mat(k,323) = -rxt(k,290)*y(k,111)
         mat(k,612) = -rxt(k,295)*y(k,111)
         mat(k,369) = -rxt(k,296)*y(k,111)
         mat(k,281) = -rxt(k,297)*y(k,111)
         mat(k,245) = -rxt(k,298)*y(k,111)
         mat(k,207) = -rxt(k,306)*y(k,111)
         mat(k,51) = -rxt(k,307)*y(k,111)
         mat(k,689) = -rxt(k,309)*y(k,111)
         mat(k,618) = -rxt(k,310)*y(k,111)
         mat(k,202) = -rxt(k,311)*y(k,111)
         mat(k,292) = -rxt(k,314)*y(k,111)
         mat(k,219) = -rxt(k,318)*y(k,111)
         mat(k,636) = -rxt(k,319)*y(k,111)
         mat(k,398) = -rxt(k,320)*y(k,111)
         mat(k,663) = -rxt(k,328)*y(k,111)
         mat(k,167) = -rxt(k,329)*y(k,111)
         mat(k,310) = -rxt(k,336)*y(k,111)
         mat(k,756) = -rxt(k,338)*y(k,111)
         mat(k,34) = -rxt(k,341)*y(k,111)
         mat(k,180) = -rxt(k,342)*y(k,111)
         mat(k,487) = -rxt(k,353)*y(k,111)
         mat(k,273) = -rxt(k,354)*y(k,111)
         mat(k,592) = -rxt(k,356)*y(k,111)
         mat(k,100) = -rxt(k,362)*y(k,111)
         mat(k,351) = -rxt(k,365)*y(k,111)
         mat(k,140) = -(rxt(k,367) + rxt(k,381)) * y(k,111)
         mat(k,186) = -rxt(k,369)*y(k,111)
         mat(k,551) = -rxt(k,371)*y(k,111)
         mat(k,362) = -rxt(k,375)*y(k,111)
         mat(k,850) = -rxt(k,380)*y(k,111)
         mat(k,36) = -rxt(k,383)*y(k,111)
         mat(k,245) = mat(k,245) + .130_r8*rxt(k,274)*y(k,108)
         mat(k,172) = mat(k,172) + .500_r8*rxt(k,280)*y(k,111)
         mat(k,636) = mat(k,636) + .360_r8*rxt(k,302)*y(k,108)
         mat(k,1029) = mat(k,1029) + rxt(k,260)*y(k,106)
         mat(k,803) = .450_r8*rxt(k,287)*y(k,139)
         mat(k,212) = mat(k,212) + .300_r8*rxt(k,267)*y(k,111)
         mat(k,1531) = mat(k,1531) + rxt(k,269)*y(k,144)
         mat(k,1324) = rxt(k,184)*y(k,139)
         mat(k,495) = rxt(k,238)*y(k,151)
         mat(k,1552) = rxt(k,143)*y(k,108) + 2.000_r8*rxt(k,138)*y(k,139)
         mat(k,833) = mat(k,833) + rxt(k,135)*y(k,106) + rxt(k,117)*y(k,144)
         mat(k,315) = mat(k,315) + rxt(k,136)*y(k,106)
         mat(k,437) = mat(k,437) + rxt(k,228)*y(k,106) + rxt(k,234)*y(k,144)
         mat(k,869) = mat(k,869) + rxt(k,199)*y(k,106) + rxt(k,211)*y(k,144)
         mat(k,391) = rxt(k,230)*y(k,106)
         mat(k,455) = mat(k,455) + rxt(k,202)*y(k,106)
         mat(k,487) = mat(k,487) + .320_r8*rxt(k,352)*y(k,108)
         mat(k,419) = .206_r8*rxt(k,344)*y(k,139)
         mat(k,663) = mat(k,663) + .240_r8*rxt(k,327)*y(k,108)
         mat(k,167) = mat(k,167) + .100_r8*rxt(k,329)*y(k,111)
         mat(k,775) = .450_r8*rxt(k,332)*y(k,139)
         mat(k,756) = mat(k,756) + .360_r8*rxt(k,337)*y(k,108)
         mat(k,1506) = rxt(k,168)*y(k,139)
         mat(k,1244) = mat(k,1244) + rxt(k,163)*y(k,139)
         mat(k,1458) = mat(k,1458) + rxt(k,260)*y(k,27) + rxt(k,135)*y(k,63) &
                      + rxt(k,136)*y(k,65) + rxt(k,228)*y(k,67) + rxt(k,199)*y(k,71) &
                      + rxt(k,230)*y(k,75) + rxt(k,202)*y(k,76) + rxt(k,141)*y(k,139)
         mat(k,1134) = mat(k,1134) + .130_r8*rxt(k,274)*y(k,9) + .360_r8*rxt(k,302) &
                      *y(k,14) + rxt(k,143)*y(k,62) + .320_r8*rxt(k,352)*y(k,79) &
                      + .240_r8*rxt(k,327)*y(k,83) + .360_r8*rxt(k,337)*y(k,88) &
                      + 1.156_r8*rxt(k,364)*y(k,133) + rxt(k,142)*y(k,139)
         mat(k,1413) = mat(k,1413) + .500_r8*rxt(k,280)*y(k,12) + .300_r8*rxt(k,267) &
                      *y(k,39) + .100_r8*rxt(k,329)*y(k,85) + .500_r8*rxt(k,314) &
                      *y(k,117) + .500_r8*rxt(k,362)*y(k,135)
         mat(k,292) = mat(k,292) + .500_r8*rxt(k,314)*y(k,111)
         mat(k,603) = .150_r8*rxt(k,316)*y(k,139)
         mat(k,351) = mat(k,351) + 1.156_r8*rxt(k,364)*y(k,108)
         mat(k,100) = mat(k,100) + .500_r8*rxt(k,362)*y(k,111)
         mat(k,1088) = mat(k,1088) + .450_r8*rxt(k,287)*y(k,32) + rxt(k,184)*y(k,42) &
                      + 2.000_r8*rxt(k,138)*y(k,62) + .206_r8*rxt(k,344)*y(k,80) &
                      + .450_r8*rxt(k,332)*y(k,86) + rxt(k,168)*y(k,97) + rxt(k,163) &
                      *y(k,99) + rxt(k,141)*y(k,106) + rxt(k,142)*y(k,108) &
                      + .150_r8*rxt(k,316)*y(k,118)
         mat(k,1286) = rxt(k,269)*y(k,40) + rxt(k,117)*y(k,63) + rxt(k,234)*y(k,67) &
                      + rxt(k,211)*y(k,71) + 2.000_r8*rxt(k,118)*y(k,151)
         mat(k,1579) = rxt(k,238)*y(k,59) + 2.000_r8*rxt(k,118)*y(k,144)
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
         mat(k,584) = -(rxt(k,355)*y(k,99) + rxt(k,356)*y(k,111))
         mat(k,1214) = -rxt(k,355)*y(k,112)
         mat(k,1381) = -rxt(k,356)*y(k,112)
         mat(k,412) = .794_r8*rxt(k,345)*y(k,97) + .794_r8*rxt(k,346)*y(k,99) &
                      + .794_r8*rxt(k,344)*y(k,139)
         mat(k,695) = .080_r8*rxt(k,350)*y(k,97)
         mat(k,720) = .800_r8*rxt(k,326)*y(k,97)
         mat(k,1478) = .794_r8*rxt(k,345)*y(k,80) + .080_r8*rxt(k,350)*y(k,81) &
                      + .800_r8*rxt(k,326)*y(k,84)
         mat(k,1214) = mat(k,1214) + .794_r8*rxt(k,346)*y(k,80)
         mat(k,1058) = .794_r8*rxt(k,344)*y(k,80)
         mat(k,274) = -(rxt(k,297)*y(k,111))
         mat(k,1358) = -rxt(k,297)*y(k,113)
         mat(k,779) = rxt(k,299)*y(k,98)
         mat(k,977) = rxt(k,299)*y(k,32)
         mat(k,401) = -(rxt(k,312)*y(k,139) + rxt(k,313)*y(k,97))
         mat(k,1052) = -rxt(k,312)*y(k,114)
         mat(k,1470) = -rxt(k,313)*y(k,114)
         mat(k,621) = rxt(k,319)*y(k,111)
         mat(k,1371) = rxt(k,319)*y(k,14) + .500_r8*rxt(k,314)*y(k,117)
         mat(k,287) = .500_r8*rxt(k,314)*y(k,111)
         mat(k,286) = -(rxt(k,314)*y(k,111))
         mat(k,1360) = -rxt(k,314)*y(k,117)
         mat(k,400) = rxt(k,312)*y(k,139)
         mat(k,1045) = rxt(k,312)*y(k,114)
         mat(k,595) = -(rxt(k,315)*y(k,37) + rxt(k,316)*y(k,139) + rxt(k,317)*y(k,97))
         mat(k,1172) = -rxt(k,315)*y(k,118)
         mat(k,1059) = -rxt(k,316)*y(k,118)
         mat(k,1479) = -rxt(k,317)*y(k,118)
         mat(k,395) = rxt(k,320)*y(k,111)
         mat(k,1382) = rxt(k,320)*y(k,33) + rxt(k,318)*y(k,119)
         mat(k,216) = rxt(k,318)*y(k,111)
         mat(k,215) = -(rxt(k,318)*y(k,111))
         mat(k,1352) = -rxt(k,318)*y(k,119)
         mat(k,593) = .850_r8*rxt(k,316)*y(k,139)
         mat(k,1041) = .850_r8*rxt(k,316)*y(k,118)
         mat(k,358) = -(rxt(k,370)*y(k,107) + rxt(k,372)*y(k,108) + rxt(k,375) &
                      *y(k,111))
         mat(k,895) = -rxt(k,370)*y(k,120)
         mat(k,1099) = -rxt(k,372)*y(k,120)
         mat(k,1367) = -rxt(k,375)*y(k,120)
         mat(k,841) = -(rxt(k,373)*y(k,6) + rxt(k,374)*y(k,45) + rxt(k,376)*y(k,98) &
                      + rxt(k,377)*y(k,107) + rxt(k,378)*y(k,108) + rxt(k,379) &
                      *y(k,109) + rxt(k,380)*y(k,111))
         mat(k,928) = -rxt(k,373)*y(k,122)
         mat(k,954) = -rxt(k,374)*y(k,122)
         mat(k,990) = -rxt(k,376)*y(k,122)
         mat(k,907) = -rxt(k,377)*y(k,122)
         mat(k,1120) = -rxt(k,378)*y(k,122)
         mat(k,175) = -rxt(k,379)*y(k,122)
         mat(k,1398) = -rxt(k,380)*y(k,122)
         mat(k,1443) = rxt(k,368)*y(k,110)
         mat(k,907) = mat(k,907) + rxt(k,370)*y(k,120)
         mat(k,1120) = mat(k,1120) + rxt(k,372)*y(k,120)
         mat(k,185) = rxt(k,368)*y(k,106)
         mat(k,1398) = mat(k,1398) + rxt(k,375)*y(k,120)
         mat(k,359) = rxt(k,370)*y(k,107) + rxt(k,372)*y(k,108) + rxt(k,375)*y(k,111)
         mat(k,548) = -(rxt(k,371)*y(k,111))
         mat(k,1378) = -rxt(k,371)*y(k,123)
         mat(k,927) = rxt(k,373)*y(k,122)
         mat(k,952) = rxt(k,374)*y(k,122)
         mat(k,136) = rxt(k,366)*y(k,99) + (rxt(k,367)+.500_r8*rxt(k,381))*y(k,111)
         mat(k,981) = rxt(k,376)*y(k,122)
         mat(k,1213) = rxt(k,366)*y(k,52)
         mat(k,902) = rxt(k,377)*y(k,122)
         mat(k,1103) = rxt(k,378)*y(k,122)
         mat(k,174) = rxt(k,379)*y(k,122)
         mat(k,183) = rxt(k,369)*y(k,111)
         mat(k,1378) = mat(k,1378) + (rxt(k,367)+.500_r8*rxt(k,381))*y(k,52) &
                      + rxt(k,369)*y(k,110) + rxt(k,380)*y(k,122)
         mat(k,840) = rxt(k,373)*y(k,6) + rxt(k,374)*y(k,45) + rxt(k,376)*y(k,98) &
                      + rxt(k,377)*y(k,107) + rxt(k,378)*y(k,108) + rxt(k,379) &
                      *y(k,109) + rxt(k,380)*y(k,111)
         mat(k,131) = -(rxt(k,382)*y(k,151))
         mat(k,1559) = -rxt(k,382)*y(k,124)
         mat(k,1340) = rxt(k,371)*y(k,123)
         mat(k,547) = rxt(k,371)*y(k,111)
         mat(k,343) = -(rxt(k,363)*y(k,99) + rxt(k,364)*y(k,108) + rxt(k,365)*y(k,111))
         mat(k,1208) = -rxt(k,363)*y(k,133)
         mat(k,1097) = -rxt(k,364)*y(k,133)
         mat(k,1365) = -rxt(k,365)*y(k,133)
         mat(k,671) = -(rxt(k,357)*y(k,32) + rxt(k,358)*y(k,37) + rxt(k,359)*y(k,139) &
                      + rxt(k,360)*y(k,97) + rxt(k,361)*y(k,99))
         mat(k,787) = -rxt(k,357)*y(k,134)
         mat(k,1178) = -rxt(k,358)*y(k,134)
         mat(k,1065) = -rxt(k,359)*y(k,134)
         mat(k,1485) = -rxt(k,360)*y(k,134)
         mat(k,1222) = -rxt(k,361)*y(k,134)
         mat(k,179) = rxt(k,342)*y(k,111)
         mat(k,268) = .800_r8*rxt(k,354)*y(k,111)
         mat(k,1389) = rxt(k,342)*y(k,78) + .800_r8*rxt(k,354)*y(k,82) &
                      + .500_r8*rxt(k,362)*y(k,135)
         mat(k,99) = .500_r8*rxt(k,362)*y(k,111)
         mat(k,98) = -(rxt(k,362)*y(k,111))
         mat(k,1334) = -rxt(k,362)*y(k,135)
         mat(k,665) = rxt(k,359)*y(k,139)
         mat(k,1036) = rxt(k,359)*y(k,134)
         mat(k,526) = -(rxt(k,411)*y(k,142) + rxt(k,412)*y(k,147) + rxt(k,413) &
                      *y(k,141))
         mat(k,371) = -rxt(k,411)*y(k,138)
         mat(k,512) = -rxt(k,412)*y(k,138)
         mat(k,338) = -rxt(k,413)*y(k,138)
         mat(k,1081) = -((rxt(k,137) + rxt(k,138) + rxt(k,139)) * y(k,62) + rxt(k,141) &
                      *y(k,106) + rxt(k,142)*y(k,108) + rxt(k,146)*y(k,111) &
                      + 4._r8*rxt(k,151)*y(k,139) + rxt(k,163)*y(k,99) + rxt(k,168) &
                      *y(k,97) + rxt(k,173)*y(k,98) + (rxt(k,183) + rxt(k,184) &
                      ) * y(k,42) + rxt(k,190)*y(k,45) + rxt(k,216)*y(k,4) + rxt(k,222) &
                      *y(k,6) + rxt(k,264)*y(k,37) + rxt(k,277)*y(k,10) + rxt(k,287) &
                      *y(k,32) + rxt(k,291)*y(k,57) + rxt(k,304)*y(k,15) + rxt(k,312) &
                      *y(k,114) + rxt(k,316)*y(k,118) + rxt(k,323)*y(k,84) + rxt(k,332) &
                      *y(k,86) + rxt(k,344)*y(k,80) + rxt(k,349)*y(k,81) + rxt(k,359) &
                      *y(k,134))
         mat(k,1545) = -(rxt(k,137) + rxt(k,138) + rxt(k,139)) * y(k,139)
         mat(k,1451) = -rxt(k,141)*y(k,139)
         mat(k,1127) = -rxt(k,142)*y(k,139)
         mat(k,1406) = -rxt(k,146)*y(k,139)
         mat(k,1237) = -rxt(k,163)*y(k,139)
         mat(k,1499) = -rxt(k,168)*y(k,139)
         mat(k,998) = -rxt(k,173)*y(k,139)
         mat(k,1317) = -(rxt(k,183) + rxt(k,184)) * y(k,139)
         mat(k,962) = -rxt(k,190)*y(k,139)
         mat(k,881) = -rxt(k,216)*y(k,139)
         mat(k,936) = -rxt(k,222)*y(k,139)
         mat(k,1191) = -rxt(k,264)*y(k,139)
         mat(k,429) = -rxt(k,277)*y(k,139)
         mat(k,798) = -rxt(k,287)*y(k,139)
         mat(k,382) = -rxt(k,291)*y(k,139)
         mat(k,467) = -rxt(k,304)*y(k,139)
         mat(k,407) = -rxt(k,312)*y(k,139)
         mat(k,601) = -rxt(k,316)*y(k,139)
         mat(k,732) = -rxt(k,323)*y(k,139)
         mat(k,770) = -rxt(k,332)*y(k,139)
         mat(k,417) = -rxt(k,344)*y(k,139)
         mat(k,710) = -rxt(k,349)*y(k,139)
         mat(k,677) = -rxt(k,359)*y(k,139)
         mat(k,881) = mat(k,881) + rxt(k,215)*y(k,27)
         mat(k,936) = mat(k,936) + rxt(k,227)*y(k,111)
         mat(k,243) = .130_r8*rxt(k,274)*y(k,108)
         mat(k,429) = mat(k,429) + 2.400_r8*rxt(k,275)*y(k,10) + rxt(k,276)*y(k,37) &
                      + rxt(k,278)*y(k,97)
         mat(k,128) = rxt(k,279)*y(k,111)
         mat(k,632) = .280_r8*rxt(k,302)*y(k,108)
         mat(k,467) = mat(k,467) + rxt(k,303)*y(k,37) + rxt(k,305)*y(k,97)
         mat(k,1022) = rxt(k,215)*y(k,4) + rxt(k,179)*y(k,42) + rxt(k,259)*y(k,99) &
                      + rxt(k,260)*y(k,106)
         mat(k,328) = rxt(k,244)*y(k,42) + rxt(k,245)*y(k,111)
         mat(k,192) = rxt(k,247)*y(k,42) + rxt(k,248)*y(k,111)
         mat(k,798) = mat(k,798) + .900_r8*rxt(k,286)*y(k,37) + rxt(k,347)*y(k,81) &
                      + .470_r8*rxt(k,321)*y(k,84) + rxt(k,357)*y(k,134)
         mat(k,1191) = mat(k,1191) + rxt(k,276)*y(k,10) + rxt(k,303)*y(k,15) &
                      + .900_r8*rxt(k,286)*y(k,32) + 4.000_r8*rxt(k,262)*y(k,37) &
                      + rxt(k,186)*y(k,45) + rxt(k,348)*y(k,81) + .730_r8*rxt(k,322) &
                      *y(k,84) + rxt(k,331)*y(k,86) + rxt(k,265)*y(k,97) &
                      + .300_r8*rxt(k,315)*y(k,118) + .800_r8*rxt(k,358)*y(k,134)
         mat(k,284) = rxt(k,266)*y(k,111)
         mat(k,1524) = rxt(k,270)*y(k,144)
         mat(k,1317) = mat(k,1317) + rxt(k,179)*y(k,27) + rxt(k,244)*y(k,28) &
                      + rxt(k,247)*y(k,31) + rxt(k,182)*y(k,65)
         mat(k,962) = mat(k,962) + rxt(k,186)*y(k,37) + rxt(k,197)*y(k,111)
         mat(k,651) = rxt(k,272)*y(k,111)
         mat(k,137) = .500_r8*rxt(k,381)*y(k,111)
         mat(k,231) = rxt(k,294)*y(k,107)
         mat(k,382) = mat(k,382) + .250_r8*rxt(k,292)*y(k,97)
         mat(k,611) = rxt(k,295)*y(k,111)
         mat(k,368) = rxt(k,296)*y(k,111)
         mat(k,1545) = mat(k,1545) + rxt(k,140)*y(k,107)
         mat(k,313) = rxt(k,182)*y(k,42) + rxt(k,136)*y(k,106) + rxt(k,145)*y(k,111)
         mat(k,617) = rxt(k,310)*y(k,111)
         mat(k,483) = .370_r8*rxt(k,352)*y(k,108)
         mat(k,417) = mat(k,417) + .794_r8*rxt(k,345)*y(k,97) + .794_r8*rxt(k,346) &
                      *y(k,99)
         mat(k,710) = mat(k,710) + rxt(k,347)*y(k,32) + rxt(k,348)*y(k,37) &
                      + .920_r8*rxt(k,350)*y(k,97) + rxt(k,351)*y(k,99)
         mat(k,661) = .140_r8*rxt(k,327)*y(k,108)
         mat(k,732) = mat(k,732) + .470_r8*rxt(k,321)*y(k,32) + .730_r8*rxt(k,322) &
                      *y(k,37) + .470_r8*rxt(k,325)*y(k,97) + .470_r8*rxt(k,324) &
                      *y(k,99)
         mat(k,166) = .200_r8*rxt(k,329)*y(k,111)
         mat(k,770) = mat(k,770) + rxt(k,331)*y(k,37)
         mat(k,308) = .500_r8*rxt(k,336)*y(k,111)
         mat(k,751) = .280_r8*rxt(k,337)*y(k,108)
         mat(k,1499) = mat(k,1499) + rxt(k,278)*y(k,10) + rxt(k,305)*y(k,15) &
                      + rxt(k,265)*y(k,37) + .250_r8*rxt(k,292)*y(k,57) &
                      + .794_r8*rxt(k,345)*y(k,80) + .920_r8*rxt(k,350)*y(k,81) &
                      + .470_r8*rxt(k,325)*y(k,84) + rxt(k,313)*y(k,114) + rxt(k,360) &
                      *y(k,134)
         mat(k,1237) = mat(k,1237) + rxt(k,259)*y(k,27) + .794_r8*rxt(k,346)*y(k,80) &
                      + rxt(k,351)*y(k,81) + .470_r8*rxt(k,324)*y(k,84) + rxt(k,166) &
                      *y(k,111) + rxt(k,355)*y(k,112) + rxt(k,361)*y(k,134)
         mat(k,1451) = mat(k,1451) + rxt(k,260)*y(k,27) + rxt(k,136)*y(k,65)
         mat(k,914) = rxt(k,294)*y(k,56) + rxt(k,140)*y(k,62)
         mat(k,1127) = mat(k,1127) + .130_r8*rxt(k,274)*y(k,9) + .280_r8*rxt(k,302) &
                      *y(k,14) + .370_r8*rxt(k,352)*y(k,79) + .140_r8*rxt(k,327) &
                      *y(k,83) + .280_r8*rxt(k,337)*y(k,88) + rxt(k,148)*y(k,111)
         mat(k,1406) = mat(k,1406) + rxt(k,227)*y(k,6) + rxt(k,279)*y(k,11) &
                      + rxt(k,245)*y(k,28) + rxt(k,248)*y(k,31) + rxt(k,266)*y(k,38) &
                      + rxt(k,197)*y(k,45) + rxt(k,272)*y(k,48) + .500_r8*rxt(k,381) &
                      *y(k,52) + rxt(k,295)*y(k,60) + rxt(k,296)*y(k,61) + rxt(k,145) &
                      *y(k,65) + rxt(k,310)*y(k,77) + .200_r8*rxt(k,329)*y(k,85) &
                      + .500_r8*rxt(k,336)*y(k,87) + rxt(k,166)*y(k,99) + rxt(k,148) &
                      *y(k,108) + rxt(k,356)*y(k,112) + rxt(k,371)*y(k,123)
         mat(k,589) = rxt(k,355)*y(k,99) + rxt(k,356)*y(k,111)
         mat(k,407) = mat(k,407) + rxt(k,313)*y(k,97)
         mat(k,601) = mat(k,601) + .300_r8*rxt(k,315)*y(k,37)
         mat(k,550) = rxt(k,371)*y(k,111)
         mat(k,677) = mat(k,677) + rxt(k,357)*y(k,32) + .800_r8*rxt(k,358)*y(k,37) &
                      + rxt(k,360)*y(k,97) + rxt(k,361)*y(k,99)
         mat(k,1279) = rxt(k,270)*y(k,40)
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
         mat(k,539) = -(rxt(k,153)*y(k,106) + rxt(k,154)*y(k,107) + rxt(k,425) &
                      *y(k,148))
         mat(k,1436) = -rxt(k,153)*y(k,140)
         mat(k,901) = -rxt(k,154)*y(k,140)
         mat(k,504) = -rxt(k,425)*y(k,140)
         mat(k,1436) = mat(k,1436) + rxt(k,415)*y(k,141)
         mat(k,527) = .900_r8*rxt(k,413)*y(k,141) + .800_r8*rxt(k,411)*y(k,142)
         mat(k,339) = rxt(k,415)*y(k,106) + .900_r8*rxt(k,413)*y(k,138)
         mat(k,372) = .800_r8*rxt(k,411)*y(k,138)
         mat(k,334) = -(rxt(k,413)*y(k,138) + rxt(k,414)*y(k,107) + (rxt(k,415) &
                      + rxt(k,416)) * y(k,106))
         mat(k,521) = -rxt(k,413)*y(k,141)
         mat(k,894) = -rxt(k,414)*y(k,141)
         mat(k,1426) = -(rxt(k,415) + rxt(k,416)) * y(k,141)
         mat(k,370) = -(rxt(k,411)*y(k,138))
         mat(k,522) = -rxt(k,411)*y(k,142)
         mat(k,555) = rxt(k,420)*y(k,147)
         mat(k,1467) = rxt(k,422)*y(k,147)
         mat(k,1428) = rxt(k,415)*y(k,141)
         mat(k,896) = rxt(k,419)*y(k,143)
         mat(k,335) = rxt(k,415)*y(k,106)
         mat(k,441) = rxt(k,419)*y(k,107)
         mat(k,510) = rxt(k,420)*y(k,89) + rxt(k,422)*y(k,97)
         mat(k,442) = -(rxt(k,417)*y(k,106) + (rxt(k,418) + rxt(k,419)) * y(k,107))
         mat(k,1431) = -rxt(k,417)*y(k,143)
         mat(k,897) = -(rxt(k,418) + rxt(k,419)) * y(k,143)
         mat(k,535) = rxt(k,425)*y(k,148)
         mat(k,500) = rxt(k,425)*y(k,140)
         mat(k,1284) = -(rxt(k,117)*y(k,63) + rxt(k,118)*y(k,151) + (rxt(k,120) &
                      + rxt(k,121)) * y(k,107) + (rxt(k,122) + rxt(k,123)) * y(k,108) &
                      + (rxt(k,171) + rxt(k,172)) * y(k,90) + rxt(k,204)*y(k,18) &
                      + rxt(k,205)*y(k,19) + rxt(k,206)*y(k,21) + rxt(k,207)*y(k,22) &
                      + rxt(k,208)*y(k,23) + rxt(k,209)*y(k,24) + rxt(k,210)*y(k,25) &
                      + (rxt(k,211) + rxt(k,212)) * y(k,71) + rxt(k,231)*y(k,20) &
                      + rxt(k,232)*y(k,41) + rxt(k,233)*y(k,64) + (rxt(k,234) &
                      + rxt(k,235)) * y(k,67) + rxt(k,240)*y(k,50) + rxt(k,241) &
                      *y(k,51) + rxt(k,254)*y(k,26) + rxt(k,255)*y(k,28) + rxt(k,256) &
                      *y(k,68) + rxt(k,257)*y(k,69) + rxt(k,258)*y(k,70) + (rxt(k,269) &
                      + rxt(k,270) + rxt(k,271)) * y(k,40))
         mat(k,831) = -rxt(k,117)*y(k,144)
         mat(k,1577) = -rxt(k,118)*y(k,144)
         mat(k,916) = -(rxt(k,120) + rxt(k,121)) * y(k,144)
         mat(k,1132) = -(rxt(k,122) + rxt(k,123)) * y(k,144)
         mat(k,124) = -(rxt(k,171) + rxt(k,172)) * y(k,144)
         mat(k,39) = -rxt(k,204)*y(k,144)
         mat(k,74) = -rxt(k,205)*y(k,144)
         mat(k,55) = -rxt(k,206)*y(k,144)
         mat(k,84) = -rxt(k,207)*y(k,144)
         mat(k,59) = -rxt(k,208)*y(k,144)
         mat(k,89) = -rxt(k,209)*y(k,144)
         mat(k,63) = -rxt(k,210)*y(k,144)
         mat(k,867) = -(rxt(k,211) + rxt(k,212)) * y(k,144)
         mat(k,80) = -rxt(k,231)*y(k,144)
         mat(k,224) = -rxt(k,232)*y(k,144)
         mat(k,44) = -rxt(k,233)*y(k,144)
         mat(k,436) = -(rxt(k,234) + rxt(k,235)) * y(k,144)
         mat(k,103) = -rxt(k,240)*y(k,144)
         mat(k,114) = -rxt(k,241)*y(k,144)
         mat(k,249) = -rxt(k,254)*y(k,144)
         mat(k,330) = -rxt(k,255)*y(k,144)
         mat(k,109) = -rxt(k,256)*y(k,144)
         mat(k,119) = -rxt(k,257)*y(k,144)
         mat(k,144) = -rxt(k,258)*y(k,144)
         mat(k,1529) = -(rxt(k,269) + rxt(k,270) + rxt(k,271)) * y(k,144)
         mat(k,916) = mat(k,916) + rxt(k,154)*y(k,140)
         mat(k,532) = .850_r8*rxt(k,412)*y(k,147)
         mat(k,544) = rxt(k,154)*y(k,107)
         mat(k,516) = .850_r8*rxt(k,412)*y(k,138)
         mat(k,93) = -(rxt(k,125)*y(k,106) + rxt(k,126)*y(k,107))
         mat(k,1419) = -rxt(k,125)*y(k,145)
         mat(k,889) = -rxt(k,126)*y(k,145)
         mat(k,807) = rxt(k,127)*y(k,146)
         mat(k,1419) = mat(k,1419) + rxt(k,129)*y(k,146)
         mat(k,889) = mat(k,889) + rxt(k,130)*y(k,146)
         mat(k,1094) = rxt(k,131)*y(k,146)
         mat(k,95) = rxt(k,127)*y(k,49) + rxt(k,129)*y(k,106) + rxt(k,130)*y(k,107) &
                      + rxt(k,131)*y(k,108)
         mat(k,96) = -(rxt(k,127)*y(k,49) + rxt(k,129)*y(k,106) + rxt(k,130)*y(k,107) &
                      + rxt(k,131)*y(k,108))
         mat(k,808) = -rxt(k,127)*y(k,146)
         mat(k,1420) = -rxt(k,129)*y(k,146)
         mat(k,890) = -rxt(k,130)*y(k,146)
         mat(k,1095) = -rxt(k,131)*y(k,146)
         mat(k,890) = mat(k,890) + rxt(k,120)*y(k,144)
         mat(k,1259) = rxt(k,120)*y(k,107)
         mat(k,511) = -(rxt(k,412)*y(k,138) + rxt(k,420)*y(k,89) + rxt(k,422)*y(k,97))
         mat(k,525) = -rxt(k,412)*y(k,147)
         mat(k,558) = -rxt(k,420)*y(k,147)
         mat(k,1474) = -rxt(k,422)*y(k,147)
         mat(k,810) = rxt(k,423)*y(k,148)
         mat(k,899) = rxt(k,414)*y(k,141) + rxt(k,418)*y(k,143) + rxt(k,426)*y(k,148) &
                      + rxt(k,430)*y(k,149)
         mat(k,337) = rxt(k,414)*y(k,107)
         mat(k,444) = rxt(k,418)*y(k,107)
         mat(k,502) = rxt(k,423)*y(k,49) + rxt(k,426)*y(k,107)
         mat(k,263) = rxt(k,430)*y(k,107)
         mat(k,501) = -(rxt(k,423)*y(k,49) + rxt(k,425)*y(k,140) + rxt(k,426)*y(k,107))
         mat(k,809) = -rxt(k,423)*y(k,148)
         mat(k,536) = -rxt(k,425)*y(k,148)
         mat(k,898) = -rxt(k,426)*y(k,148)
         mat(k,1433) = rxt(k,416)*y(k,141) + rxt(k,417)*y(k,143) + rxt(k,429)*y(k,149) &
                      + rxt(k,435)*y(k,150)
         mat(k,524) = rxt(k,427)*y(k,149) + rxt(k,432)*y(k,150)
         mat(k,336) = rxt(k,416)*y(k,106)
         mat(k,443) = rxt(k,417)*y(k,106)
         mat(k,262) = rxt(k,429)*y(k,106) + rxt(k,427)*y(k,138)
         mat(k,257) = rxt(k,435)*y(k,106) + rxt(k,432)*y(k,138)
         mat(k,260) = -(rxt(k,427)*y(k,138) + rxt(k,429)*y(k,106) + rxt(k,430) &
                      *y(k,107))
         mat(k,520) = -rxt(k,427)*y(k,149)
         mat(k,1423) = -rxt(k,429)*y(k,149)
         mat(k,893) = -rxt(k,430)*y(k,149)
         mat(k,520) = mat(k,520) + rxt(k,431)*y(k,150)
         mat(k,254) = rxt(k,431)*y(k,138)
         mat(k,253) = -((rxt(k,431) + rxt(k,432)) * y(k,138) + rxt(k,435)*y(k,106))
         mat(k,519) = -(rxt(k,431) + rxt(k,432)) * y(k,150)
         mat(k,1422) = -rxt(k,435)*y(k,150)
         mat(k,1584) = -(rxt(k,118)*y(k,144) + rxt(k,238)*y(k,59) + rxt(k,382) &
                      *y(k,124))
         mat(k,1291) = -rxt(k,118)*y(k,151)
         mat(k,498) = -rxt(k,238)*y(k,151)
         mat(k,134) = -rxt(k,382)*y(k,151)
         mat(k,162) = rxt(k,282)*y(k,111)
         mat(k,208) = rxt(k,306)*y(k,111)
         mat(k,52) = rxt(k,307)*y(k,111)
         mat(k,252) = rxt(k,243)*y(k,111)
         mat(k,1034) = rxt(k,261)*y(k,111)
         mat(k,333) = rxt(k,245)*y(k,111)
         mat(k,48) = rxt(k,246)*y(k,111)
         mat(k,648) = rxt(k,284)*y(k,111)
         mat(k,196) = rxt(k,248)*y(k,111)
         mat(k,399) = rxt(k,320)*y(k,111)
         mat(k,690) = rxt(k,309)*y(k,111)
         mat(k,357) = rxt(k,289)*y(k,111)
         mat(k,324) = rxt(k,290)*y(k,111)
         mat(k,214) = rxt(k,267)*y(k,111)
         mat(k,1536) = rxt(k,268)*y(k,111)
         mat(k,1557) = rxt(k,139)*y(k,139)
         mat(k,837) = rxt(k,144)*y(k,111)
         mat(k,317) = rxt(k,145)*y(k,111)
         mat(k,440) = rxt(k,229)*y(k,111)
         mat(k,147) = rxt(k,253)*y(k,111)
         mat(k,872) = (rxt(k,399)+rxt(k,404))*y(k,75) + (rxt(k,392)+rxt(k,398) &
                       +rxt(k,403))*y(k,76) + rxt(k,200)*y(k,111)
         mat(k,1162) = rxt(k,176)*y(k,111)
         mat(k,238) = rxt(k,152)*y(k,111)
         mat(k,393) = (rxt(k,399)+rxt(k,404))*y(k,71)
         mat(k,457) = (rxt(k,392)+rxt(k,398)+rxt(k,403))*y(k,71) + rxt(k,203)*y(k,111)
         mat(k,664) = .500_r8*rxt(k,328)*y(k,111)
         mat(k,37) = rxt(k,383)*y(k,111)
         mat(k,1418) = rxt(k,282)*y(k,13) + rxt(k,306)*y(k,16) + rxt(k,307)*y(k,17) &
                      + rxt(k,243)*y(k,26) + rxt(k,261)*y(k,27) + rxt(k,245)*y(k,28) &
                      + rxt(k,246)*y(k,29) + rxt(k,284)*y(k,30) + rxt(k,248)*y(k,31) &
                      + rxt(k,320)*y(k,33) + rxt(k,309)*y(k,34) + rxt(k,289)*y(k,35) &
                      + rxt(k,290)*y(k,36) + rxt(k,267)*y(k,39) + rxt(k,268)*y(k,40) &
                      + rxt(k,144)*y(k,63) + rxt(k,145)*y(k,65) + rxt(k,229)*y(k,67) &
                      + rxt(k,253)*y(k,70) + rxt(k,200)*y(k,71) + rxt(k,176)*y(k,73) &
                      + rxt(k,152)*y(k,74) + rxt(k,203)*y(k,76) + .500_r8*rxt(k,328) &
                      *y(k,83) + rxt(k,383)*y(k,95) + 2.000_r8*rxt(k,149)*y(k,111) &
                      + rxt(k,314)*y(k,117) + rxt(k,318)*y(k,119) + rxt(k,146) &
                      *y(k,139)
         mat(k,293) = rxt(k,314)*y(k,111)
         mat(k,220) = rxt(k,318)*y(k,111)
         mat(k,1093) = rxt(k,139)*y(k,62) + rxt(k,146)*y(k,111)
      end do
      end subroutine nlnmat08
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
         mat(k, 32) = mat(k, 32) + lmat(k, 32)
         mat(k, 35) = mat(k, 35) + lmat(k, 35)
         mat(k, 38) = mat(k, 38) + lmat(k, 38)
         mat(k, 40) = mat(k, 40) + lmat(k, 40)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 42) = mat(k, 42) + lmat(k, 42)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 69) = lmat(k, 69)
         mat(k, 70) = lmat(k, 70)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 73) = mat(k, 73) + lmat(k, 73)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 82) = mat(k, 82) + lmat(k, 82)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 85) = mat(k, 85) + lmat(k, 85)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 91) = lmat(k, 91)
         mat(k, 92) = lmat(k, 92)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 104) = lmat(k, 104)
         mat(k, 105) = lmat(k, 105)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 148) = lmat(k, 148)
         mat(k, 149) = lmat(k, 149)
         mat(k, 150) = lmat(k, 150)
         mat(k, 151) = lmat(k, 151)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 171) = lmat(k, 171)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 193) = lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 212) = mat(k, 212) + lmat(k, 212)
         mat(k, 213) = lmat(k, 213)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 217) = lmat(k, 217)
         mat(k, 218) = lmat(k, 218)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 230) = lmat(k, 230)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 234) = mat(k, 234) + lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 255) = lmat(k, 255)
         mat(k, 256) = lmat(k, 256)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 258) = lmat(k, 258)
         mat(k, 259) = lmat(k, 259)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 261) = lmat(k, 261)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 267) = lmat(k, 267)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 274) = mat(k, 274) + lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 289) = lmat(k, 289)
         mat(k, 290) = lmat(k, 290)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 295) = lmat(k, 295)
         mat(k, 296) = lmat(k, 296)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 298) = lmat(k, 298)
         mat(k, 299) = lmat(k, 299)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 304) = lmat(k, 304)
         mat(k, 306) = lmat(k, 306)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 322) = lmat(k, 322)
         mat(k, 323) = mat(k, 323) + lmat(k, 323)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 329) = lmat(k, 329)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 343) = mat(k, 343) + lmat(k, 343)
         mat(k, 353) = mat(k, 353) + lmat(k, 353)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 386) = mat(k, 386) + lmat(k, 386)
         mat(k, 388) = lmat(k, 388)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = lmat(k, 397)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 433) = mat(k, 433) + lmat(k, 433)
         mat(k, 434) = mat(k, 434) + lmat(k, 434)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 490) = mat(k, 490) + lmat(k, 490)
         mat(k, 499) = lmat(k, 499)
         mat(k, 501) = mat(k, 501) + lmat(k, 501)
         mat(k, 505) = lmat(k, 505)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 526) = mat(k, 526) + lmat(k, 526)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 549) = lmat(k, 549)
         mat(k, 552) = lmat(k, 552)
         mat(k, 556) = lmat(k, 556)
         mat(k, 559) = lmat(k, 559)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 572) = mat(k, 572) + lmat(k, 572)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 576) = lmat(k, 576)
         mat(k, 577) = mat(k, 577) + lmat(k, 577)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 579) = mat(k, 579) + lmat(k, 579)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 585) = lmat(k, 585)
         mat(k, 587) = mat(k, 587) + lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 590) = lmat(k, 590)
         mat(k, 595) = mat(k, 595) + lmat(k, 595)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 608) = lmat(k, 608)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 611) = mat(k, 611) + lmat(k, 611)
         mat(k, 613) = mat(k, 613) + lmat(k, 613)
         mat(k, 615) = lmat(k, 615)
         mat(k, 616) = lmat(k, 616)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 623) = mat(k, 623) + lmat(k, 623)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 641) = lmat(k, 641)
         mat(k, 643) = lmat(k, 643)
         mat(k, 645) = lmat(k, 645)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 653) = mat(k, 653) + lmat(k, 653)
         mat(k, 654) = mat(k, 654) + lmat(k, 654)
         mat(k, 657) = mat(k, 657) + lmat(k, 657)
         mat(k, 658) = mat(k, 658) + lmat(k, 658)
         mat(k, 660) = mat(k, 660) + lmat(k, 660)
         mat(k, 661) = mat(k, 661) + lmat(k, 661)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 682) = mat(k, 682) + lmat(k, 682)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 684) = mat(k, 684) + lmat(k, 684)
         mat(k, 686) = lmat(k, 686)
         mat(k, 702) = mat(k, 702) + lmat(k, 702)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 739) = lmat(k, 739)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 746) = mat(k, 746) + lmat(k, 746)
         mat(k, 754) = lmat(k, 754)
         mat(k, 765) = mat(k, 765) + lmat(k, 765)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 814) = mat(k, 814) + lmat(k, 814)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 821) = lmat(k, 821)
         mat(k, 826) = mat(k, 826) + lmat(k, 826)
         mat(k, 839) = lmat(k, 839)
         mat(k, 841) = mat(k, 841) + lmat(k, 841)
         mat(k, 851) = mat(k, 851) + lmat(k, 851)
         mat(k, 860) = mat(k, 860) + lmat(k, 860)
         mat(k, 868) = mat(k, 868) + lmat(k, 868)
         mat(k, 871) = mat(k, 871) + lmat(k, 871)
         mat(k, 876) = mat(k, 876) + lmat(k, 876)
         mat(k, 892) = lmat(k, 892)
         mat(k, 893) = mat(k, 893) + lmat(k, 893)
         mat(k, 898) = mat(k, 898) + lmat(k, 898)
         mat(k, 899) = mat(k, 899) + lmat(k, 899)
         mat(k, 900) = lmat(k, 900)
         mat(k, 909) = mat(k, 909) + lmat(k, 909)
         mat(k, 916) = mat(k, 916) + lmat(k, 916)
         mat(k, 919) = mat(k, 919) + lmat(k, 919)
         mat(k, 930) = mat(k, 930) + lmat(k, 930)
         mat(k, 932) = mat(k, 932) + lmat(k, 932)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 959) = mat(k, 959) + lmat(k, 959)
         mat(k, 968) = mat(k, 968) + lmat(k, 968)
         mat(k, 970) = mat(k, 970) + lmat(k, 970)
         mat(k, 996) = mat(k, 996) + lmat(k, 996)
         mat(k,1000) = mat(k,1000) + lmat(k,1000)
         mat(k,1005) = mat(k,1005) + lmat(k,1005)
         mat(k,1006) = mat(k,1006) + lmat(k,1006)
         mat(k,1007) = mat(k,1007) + lmat(k,1007)
         mat(k,1012) = mat(k,1012) + lmat(k,1012)
         mat(k,1014) = lmat(k,1014)
         mat(k,1021) = mat(k,1021) + lmat(k,1021)
         mat(k,1033) = mat(k,1033) + lmat(k,1033)
         mat(k,1081) = mat(k,1081) + lmat(k,1081)
         mat(k,1093) = mat(k,1093) + lmat(k,1093)
         mat(k,1094) = mat(k,1094) + lmat(k,1094)
         mat(k,1122) = mat(k,1122) + lmat(k,1122)
         mat(k,1128) = mat(k,1128) + lmat(k,1128)
         mat(k,1132) = mat(k,1132) + lmat(k,1132)
         mat(k,1135) = mat(k,1135) + lmat(k,1135)
         mat(k,1148) = lmat(k,1148)
         mat(k,1152) = mat(k,1152) + lmat(k,1152)
         mat(k,1157) = mat(k,1157) + lmat(k,1157)
         mat(k,1194) = mat(k,1194) + lmat(k,1194)
         mat(k,1232) = mat(k,1232) + lmat(k,1232)
         mat(k,1235) = mat(k,1235) + lmat(k,1235)
         mat(k,1239) = mat(k,1239) + lmat(k,1239)
         mat(k,1241) = mat(k,1241) + lmat(k,1241)
         mat(k,1245) = mat(k,1245) + lmat(k,1245)
         mat(k,1246) = mat(k,1246) + lmat(k,1246)
         mat(k,1284) = mat(k,1284) + lmat(k,1284)
         mat(k,1287) = mat(k,1287) + lmat(k,1287)
         mat(k,1323) = mat(k,1323) + lmat(k,1323)
         mat(k,1413) = mat(k,1413) + lmat(k,1413)
         mat(k,1422) = mat(k,1422) + lmat(k,1422)
         mat(k,1423) = mat(k,1423) + lmat(k,1423)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1435) = lmat(k,1435)
         mat(k,1459) = mat(k,1459) + lmat(k,1459)
         mat(k,1467) = mat(k,1467) + lmat(k,1467)
         mat(k,1475) = lmat(k,1475)
         mat(k,1477) = mat(k,1477) + lmat(k,1477)
         mat(k,1507) = mat(k,1507) + lmat(k,1507)
         mat(k,1508) = mat(k,1508) + lmat(k,1508)
         mat(k,1514) = lmat(k,1514)
         mat(k,1515) = lmat(k,1515)
         mat(k,1516) = mat(k,1516) + lmat(k,1516)
         mat(k,1523) = mat(k,1523) + lmat(k,1523)
         mat(k,1527) = mat(k,1527) + lmat(k,1527)
         mat(k,1531) = mat(k,1531) + lmat(k,1531)
         mat(k,1532) = lmat(k,1532)
         mat(k,1534) = mat(k,1534) + lmat(k,1534)
         mat(k,1535) = mat(k,1535) + lmat(k,1535)
         mat(k,1536) = mat(k,1536) + lmat(k,1536)
         mat(k,1556) = mat(k,1556) + lmat(k,1556)
         mat(k,1563) = lmat(k,1563)
         mat(k,1577) = mat(k,1577) + lmat(k,1577)
         mat(k,1579) = mat(k,1579) + lmat(k,1579)
         mat(k,1580) = lmat(k,1580)
         mat(k,1583) = lmat(k,1583)
         mat(k,1584) = mat(k,1584) + lmat(k,1584)
         mat(k, 118) = 0._r8
         mat(k, 143) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 379) = 0._r8
         mat(k, 383) = 0._r8
         mat(k, 390) = 0._r8
         mat(k, 402) = 0._r8
         mat(k, 408) = 0._r8
         mat(k, 410) = 0._r8
         mat(k, 431) = 0._r8
         mat(k, 445) = 0._r8
         mat(k, 446) = 0._r8
         mat(k, 461) = 0._r8
         mat(k, 463) = 0._r8
         mat(k, 469) = 0._r8
         mat(k, 471) = 0._r8
         mat(k, 474) = 0._r8
         mat(k, 481) = 0._r8
         mat(k, 488) = 0._r8
         mat(k, 503) = 0._r8
         mat(k, 513) = 0._r8
         mat(k, 523) = 0._r8
         mat(k, 529) = 0._r8
         mat(k, 530) = 0._r8
         mat(k, 531) = 0._r8
         mat(k, 534) = 0._r8
         mat(k, 537) = 0._r8
         mat(k, 538) = 0._r8
         mat(k, 541) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 553) = 0._r8
         mat(k, 557) = 0._r8
         mat(k, 560) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 563) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 582) = 0._r8
         mat(k, 586) = 0._r8
         mat(k, 605) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 627) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 630) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 680) = 0._r8
         mat(k, 685) = 0._r8
         mat(k, 698) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 707) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 712) = 0._r8
         mat(k, 715) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 729) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 745) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 757) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 759) = 0._r8
         mat(k, 763) = 0._r8
         mat(k, 764) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 777) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 805) = 0._r8
         mat(k, 806) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 812) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 817) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 819) = 0._r8
         mat(k, 820) = 0._r8
         mat(k, 822) = 0._r8
         mat(k, 823) = 0._r8
         mat(k, 828) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 830) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 847) = 0._r8
         mat(k, 854) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 862) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 866) = 0._r8
         mat(k, 875) = 0._r8
         mat(k, 879) = 0._r8
         mat(k, 883) = 0._r8
         mat(k, 884) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 886) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 888) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 908) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 912) = 0._r8
         mat(k, 913) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 922) = 0._r8
         mat(k, 929) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 940) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 972) = 0._r8
         mat(k, 973) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 986) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 991) = 0._r8
         mat(k, 992) = 0._r8
         mat(k, 997) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1003) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1010) = 0._r8
         mat(k,1013) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1019) = 0._r8
         mat(k,1020) = 0._r8
         mat(k,1023) = 0._r8
         mat(k,1025) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1032) = 0._r8
         mat(k,1042) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1060) = 0._r8
         mat(k,1061) = 0._r8
         mat(k,1062) = 0._r8
         mat(k,1063) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1083) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1100) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1112) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1129) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1192) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1202) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1210) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1230) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1233) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1240) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1248) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1277) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1329) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1411) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1441) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1461) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1476) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1504) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1510) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1518) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1522) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1526) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1547) = 0._r8
         mat(k,1548) = 0._r8
         mat(k,1549) = 0._r8
         mat(k,1550) = 0._r8
         mat(k,1551) = 0._r8
         mat(k,1554) = 0._r8
         mat(k,1555) = 0._r8
         mat(k,1562) = 0._r8
         mat(k,1564) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1566) = 0._r8
         mat(k,1567) = 0._r8
         mat(k,1568) = 0._r8
         mat(k,1569) = 0._r8
         mat(k,1570) = 0._r8
         mat(k,1571) = 0._r8
         mat(k,1572) = 0._r8
         mat(k,1573) = 0._r8
         mat(k,1574) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1576) = 0._r8
         mat(k,1578) = 0._r8
         mat(k,1581) = 0._r8
         mat(k,1582) = 0._r8
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
         mat(k, 28) = mat(k, 28) - dti(k)
         mat(k, 29) = mat(k, 29) - dti(k)
         mat(k, 32) = mat(k, 32) - dti(k)
         mat(k, 35) = mat(k, 35) - dti(k)
         mat(k, 38) = mat(k, 38) - dti(k)
         mat(k, 41) = mat(k, 41) - dti(k)
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 68) = mat(k, 68) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 81) = mat(k, 81) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 91) = mat(k, 91) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 98) = mat(k, 98) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 104) = mat(k, 104) - dti(k)
         mat(k, 107) = mat(k, 107) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 117) = mat(k, 117) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 126) = mat(k, 126) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 135) = mat(k, 135) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 148) = mat(k, 148) - dti(k)
         mat(k, 151) = mat(k, 151) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 163) = mat(k, 163) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 173) = mat(k, 173) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 181) = mat(k, 181) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 197) = mat(k, 197) - dti(k)
         mat(k, 203) = mat(k, 203) - dti(k)
         mat(k, 209) = mat(k, 209) - dti(k)
         mat(k, 215) = mat(k, 215) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 227) = mat(k, 227) - dti(k)
         mat(k, 232) = mat(k, 232) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 246) = mat(k, 246) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 260) = mat(k, 260) - dti(k)
         mat(k, 266) = mat(k, 266) - dti(k)
         mat(k, 274) = mat(k, 274) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 286) = mat(k, 286) - dti(k)
         mat(k, 294) = mat(k, 294) - dti(k)
         mat(k, 302) = mat(k, 302) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 318) = mat(k, 318) - dti(k)
         mat(k, 325) = mat(k, 325) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 343) = mat(k, 343) - dti(k)
         mat(k, 353) = mat(k, 353) - dti(k)
         mat(k, 358) = mat(k, 358) - dti(k)
         mat(k, 365) = mat(k, 365) - dti(k)
         mat(k, 370) = mat(k, 370) - dti(k)
         mat(k, 377) = mat(k, 377) - dti(k)
         mat(k, 386) = mat(k, 386) - dti(k)
         mat(k, 394) = mat(k, 394) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 411) = mat(k, 411) - dti(k)
         mat(k, 424) = mat(k, 424) - dti(k)
         mat(k, 433) = mat(k, 433) - dti(k)
         mat(k, 442) = mat(k, 442) - dti(k)
         mat(k, 451) = mat(k, 451) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 473) = mat(k, 473) - dti(k)
         mat(k, 490) = mat(k, 490) - dti(k)
         mat(k, 501) = mat(k, 501) - dti(k)
         mat(k, 511) = mat(k, 511) - dti(k)
         mat(k, 526) = mat(k, 526) - dti(k)
         mat(k, 539) = mat(k, 539) - dti(k)
         mat(k, 548) = mat(k, 548) - dti(k)
         mat(k, 561) = mat(k, 561) - dti(k)
         mat(k, 573) = mat(k, 573) - dti(k)
         mat(k, 584) = mat(k, 584) - dti(k)
         mat(k, 595) = mat(k, 595) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 613) = mat(k, 613) - dti(k)
         mat(k, 623) = mat(k, 623) - dti(k)
         mat(k, 640) = mat(k, 640) - dti(k)
         mat(k, 649) = mat(k, 649) - dti(k)
         mat(k, 654) = mat(k, 654) - dti(k)
         mat(k, 671) = mat(k, 671) - dti(k)
         mat(k, 683) = mat(k, 683) - dti(k)
         mat(k, 702) = mat(k, 702) - dti(k)
         mat(k, 726) = mat(k, 726) - dti(k)
         mat(k, 744) = mat(k, 744) - dti(k)
         mat(k, 765) = mat(k, 765) - dti(k)
         mat(k, 793) = mat(k, 793) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 826) = mat(k, 826) - dti(k)
         mat(k, 841) = mat(k, 841) - dti(k)
         mat(k, 860) = mat(k, 860) - dti(k)
         mat(k, 876) = mat(k, 876) - dti(k)
         mat(k, 909) = mat(k, 909) - dti(k)
         mat(k, 932) = mat(k, 932) - dti(k)
         mat(k, 959) = mat(k, 959) - dti(k)
         mat(k, 996) = mat(k, 996) - dti(k)
         mat(k,1021) = mat(k,1021) - dti(k)
         mat(k,1081) = mat(k,1081) - dti(k)
         mat(k,1128) = mat(k,1128) - dti(k)
         mat(k,1152) = mat(k,1152) - dti(k)
         mat(k,1194) = mat(k,1194) - dti(k)
         mat(k,1241) = mat(k,1241) - dti(k)
         mat(k,1284) = mat(k,1284) - dti(k)
         mat(k,1323) = mat(k,1323) - dti(k)
         mat(k,1413) = mat(k,1413) - dti(k)
         mat(k,1459) = mat(k,1459) - dti(k)
         mat(k,1508) = mat(k,1508) - dti(k)
         mat(k,1534) = mat(k,1534) - dti(k)
         mat(k,1556) = mat(k,1556) - dti(k)
         mat(k,1584) = mat(k,1584) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
