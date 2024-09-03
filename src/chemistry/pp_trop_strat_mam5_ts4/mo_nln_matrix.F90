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
         mat(k,39) = -(rxt(k,293)*y(k,137))
         mat(k,1092) = -rxt(k,293)*y(k,3)
         mat(k,800) = -(rxt(k,173)*y(k,25) + rxt(k,174)*y(k,132) + rxt(k,175)*y(k,98))
         mat(k,784) = -rxt(k,173)*y(k,4)
         mat(k,939) = -rxt(k,174)*y(k,4)
         mat(k,1005) = -rxt(k,175)*y(k,4)
         mat(k,961) = 4.000_r8*rxt(k,176)*y(k,6) + (rxt(k,177)+rxt(k,178))*y(k,41) &
                      + rxt(k,181)*y(k,88) + rxt(k,184)*y(k,97) + rxt(k,325)*y(k,110) &
                      + rxt(k,185)*y(k,137)
         mat(k,74) = rxt(k,163)*y(k,136)
         mat(k,52) = rxt(k,189)*y(k,136)
         mat(k,228) = 2.000_r8*rxt(k,194)*y(k,38) + 2.000_r8*rxt(k,206)*y(k,136) &
                      + 2.000_r8*rxt(k,195)*y(k,137)
         mat(k,294) = rxt(k,196)*y(k,38) + rxt(k,207)*y(k,136) + rxt(k,197)*y(k,137)
         mat(k,183) = 3.000_r8*rxt(k,201)*y(k,38) + 3.000_r8*rxt(k,190)*y(k,136) &
                      + 3.000_r8*rxt(k,202)*y(k,137)
         mat(k,854) = 2.000_r8*rxt(k,194)*y(k,24) + rxt(k,196)*y(k,26) &
                      + 3.000_r8*rxt(k,201)*y(k,37)
         mat(k,1232) = (rxt(k,177)+rxt(k,178))*y(k,6)
         mat(k,43) = 2.000_r8*rxt(k,191)*y(k,136)
         mat(k,410) = rxt(k,186)*y(k,97) + rxt(k,192)*y(k,136) + rxt(k,187)*y(k,137)
         mat(k,1274) = rxt(k,181)*y(k,6)
         mat(k,883) = rxt(k,184)*y(k,6) + rxt(k,186)*y(k,59)
         mat(k,604) = rxt(k,325)*y(k,6)
         mat(k,1078) = rxt(k,163)*y(k,17) + rxt(k,189)*y(k,18) + 2.000_r8*rxt(k,206) &
                      *y(k,24) + rxt(k,207)*y(k,26) + 3.000_r8*rxt(k,190)*y(k,37) &
                      + 2.000_r8*rxt(k,191)*y(k,56) + rxt(k,192)*y(k,59)
         mat(k,1163) = rxt(k,185)*y(k,6) + 2.000_r8*rxt(k,195)*y(k,24) + rxt(k,197) &
                      *y(k,26) + 3.000_r8*rxt(k,202)*y(k,37) + rxt(k,187)*y(k,59)
         mat(k,953) = rxt(k,179)*y(k,41)
         mat(k,1223) = rxt(k,179)*y(k,6)
         mat(k,813) = (rxt(k,351)+rxt(k,356))*y(k,67)
         mat(k,351) = (rxt(k,351)+rxt(k,356))*y(k,63)
         mat(k,966) = -(4._r8*rxt(k,176)*y(k,6) + (rxt(k,177) + rxt(k,178) + rxt(k,179) &
                      ) * y(k,41) + rxt(k,180)*y(k,132) + rxt(k,181)*y(k,88) + rxt(k,182) &
                      *y(k,89) + rxt(k,184)*y(k,97) + rxt(k,185)*y(k,137) + rxt(k,325) &
                      *y(k,110))
         mat(k,1237) = -(rxt(k,177) + rxt(k,178) + rxt(k,179)) * y(k,6)
         mat(k,944) = -rxt(k,180)*y(k,6)
         mat(k,1279) = -rxt(k,181)*y(k,6)
         mat(k,1046) = -rxt(k,182)*y(k,6)
         mat(k,888) = -rxt(k,184)*y(k,6)
         mat(k,1168) = -rxt(k,185)*y(k,6)
         mat(k,608) = -rxt(k,325)*y(k,6)
         mat(k,805) = rxt(k,175)*y(k,98)
         mat(k,273) = rxt(k,183)*y(k,97)
         mat(k,412) = rxt(k,193)*y(k,136)
         mat(k,357) = rxt(k,188)*y(k,97)
         mat(k,888) = mat(k,888) + rxt(k,183)*y(k,7) + rxt(k,188)*y(k,67)
         mat(k,1010) = rxt(k,175)*y(k,4)
         mat(k,1083) = rxt(k,193)*y(k,59)
         mat(k,268) = -(rxt(k,183)*y(k,97))
         mat(k,870) = -rxt(k,183)*y(k,7)
         mat(k,955) = rxt(k,182)*y(k,89)
         mat(k,1024) = rxt(k,182)*y(k,6)
         mat(k,220) = -(rxt(k,225)*y(k,38) + rxt(k,226)*y(k,98) + rxt(k,250)*y(k,137))
         mat(k,836) = -rxt(k,225)*y(k,9)
         mat(k,975) = -rxt(k,226)*y(k,9)
         mat(k,1116) = -rxt(k,250)*y(k,9)
         mat(k,111) = -(rxt(k,231)*y(k,137))
         mat(k,1100) = -rxt(k,231)*y(k,10)
         mat(k,366) = .800_r8*rxt(k,227)*y(k,126) + .200_r8*rxt(k,228)*y(k,129)
         mat(k,730) = .200_r8*rxt(k,228)*y(k,126)
         mat(k,147) = -(rxt(k,232)*y(k,137))
         mat(k,1105) = -rxt(k,232)*y(k,11)
         mat(k,367) = rxt(k,229)*y(k,132)
         mat(k,901) = rxt(k,229)*y(k,126)
         mat(k,130) = -(rxt(k,233)*y(k,38) + rxt(k,234)*y(k,137))
         mat(k,833) = -rxt(k,233)*y(k,12)
         mat(k,1102) = -rxt(k,234)*y(k,12)
         mat(k,530) = -(rxt(k,253)*y(k,90) + rxt(k,254)*y(k,98) + rxt(k,271)*y(k,137))
         mat(k,1193) = -rxt(k,253)*y(k,13)
         mat(k,991) = -rxt(k,254)*y(k,13)
         mat(k,1149) = -rxt(k,271)*y(k,13)
         mat(k,432) = .130_r8*rxt(k,304)*y(k,98)
         mat(k,991) = mat(k,991) + .130_r8*rxt(k,304)*y(k,71)
         mat(k,188) = -(rxt(k,258)*y(k,137))
         mat(k,1111) = -rxt(k,258)*y(k,14)
         mat(k,387) = rxt(k,256)*y(k,132)
         mat(k,903) = rxt(k,256)*y(k,127)
         mat(k,69) = -(rxt(k,259)*y(k,137))
         mat(k,1094) = -rxt(k,259)*y(k,15)
         mat(k,48) = -(rxt(k,162)*y(k,136))
         mat(k,1056) = -rxt(k,162)*y(k,16)
         mat(k,73) = -(rxt(k,163)*y(k,136))
         mat(k,1063) = -rxt(k,163)*y(k,17)
         mat(k,51) = -(rxt(k,189)*y(k,136))
         mat(k,1057) = -rxt(k,189)*y(k,18)
         mat(k,54) = -(rxt(k,164)*y(k,136))
         mat(k,1058) = -rxt(k,164)*y(k,19)
         mat(k,57) = -(rxt(k,165)*y(k,136))
         mat(k,1059) = -rxt(k,165)*y(k,20)
         mat(k,60) = -(rxt(k,166)*y(k,136))
         mat(k,1060) = -rxt(k,166)*y(k,21)
         mat(k,63) = -(rxt(k,167)*y(k,136))
         mat(k,1061) = -rxt(k,167)*y(k,22)
         mat(k,66) = -(rxt(k,168)*y(k,136))
         mat(k,1062) = -rxt(k,168)*y(k,23)
         mat(k,227) = -(rxt(k,194)*y(k,38) + rxt(k,195)*y(k,137) + rxt(k,206)*y(k,136))
         mat(k,837) = -rxt(k,194)*y(k,24)
         mat(k,1117) = -rxt(k,195)*y(k,24)
         mat(k,1069) = -rxt(k,206)*y(k,24)
         mat(k,783) = -(rxt(k,137)*y(k,38) + rxt(k,173)*y(k,4) + rxt(k,211)*y(k,90) &
                      + rxt(k,212)*y(k,97) + rxt(k,213)*y(k,137))
         mat(k,853) = -rxt(k,137)*y(k,25)
         mat(k,799) = -rxt(k,173)*y(k,25)
         mat(k,1206) = -rxt(k,211)*y(k,25)
         mat(k,882) = -rxt(k,212)*y(k,25)
         mat(k,1162) = -rxt(k,213)*y(k,25)
         mat(k,223) = rxt(k,226)*y(k,98)
         mat(k,536) = .500_r8*rxt(k,254)*y(k,98)
         mat(k,305) = .500_r8*rxt(k,242)*y(k,137)
         mat(k,249) = rxt(k,218)*y(k,137)
         mat(k,203) = .300_r8*rxt(k,219)*y(k,137)
         mat(k,482) = (rxt(k,222)+rxt(k,223))*y(k,136)
         mat(k,1231) = rxt(k,144)*y(k,129)
         mat(k,461) = .800_r8*rxt(k,247)*y(k,137)
         mat(k,438) = .910_r8*rxt(k,304)*y(k,98)
         mat(k,381) = .072_r8*rxt(k,297)*y(k,88) + .072_r8*rxt(k,298)*y(k,90) &
                      + .206_r8*rxt(k,296)*y(k,132)
         mat(k,565) = .120_r8*rxt(k,279)*y(k,98)
         mat(k,287) = .500_r8*rxt(k,288)*y(k,137)
         mat(k,718) = .600_r8*rxt(k,289)*y(k,98)
         mat(k,1273) = .072_r8*rxt(k,297)*y(k,72) + rxt(k,217)*y(k,129) &
                      + .500_r8*rxt(k,244)*y(k,131) + .550_r8*rxt(k,302)*y(k,133) &
                      + .250_r8*rxt(k,277)*y(k,134) + rxt(k,286)*y(k,135) + rxt(k,265) &
                      *y(k,138) + rxt(k,269)*y(k,139) + .250_r8*rxt(k,312)*y(k,140)
         mat(k,1206) = mat(k,1206) + .072_r8*rxt(k,298)*y(k,72) + .600_r8*rxt(k,303) &
                      *y(k,133) + .250_r8*rxt(k,276)*y(k,134) + rxt(k,287)*y(k,135)
         mat(k,1004) = rxt(k,226)*y(k,9) + .500_r8*rxt(k,254)*y(k,13) &
                      + .910_r8*rxt(k,304)*y(k,71) + .120_r8*rxt(k,279)*y(k,74) &
                      + .600_r8*rxt(k,289)*y(k,77)
         mat(k,256) = rxt(k,249)*y(k,137)
         mat(k,372) = .700_r8*rxt(k,228)*y(k,129)
         mat(k,394) = rxt(k,255)*y(k,129)
         mat(k,698) = rxt(k,238)*y(k,129) + .600_r8*rxt(k,299)*y(k,133) &
                      + .250_r8*rxt(k,273)*y(k,134) + rxt(k,282)*y(k,135) &
                      + .250_r8*rxt(k,309)*y(k,140)
         mat(k,755) = rxt(k,144)*y(k,41) + rxt(k,217)*y(k,88) + .700_r8*rxt(k,228) &
                      *y(k,126) + rxt(k,255)*y(k,127) + rxt(k,238)*y(k,128) + ( &
                      + 4.000_r8*rxt(k,214)+2.000_r8*rxt(k,215))*y(k,129) &
                      + 1.200_r8*rxt(k,300)*y(k,133) + .880_r8*rxt(k,274)*y(k,134) &
                      + 2.000_r8*rxt(k,283)*y(k,135) + .800_r8*rxt(k,267)*y(k,139) &
                      + .800_r8*rxt(k,310)*y(k,140)
         mat(k,330) = .500_r8*rxt(k,244)*y(k,88)
         mat(k,938) = .206_r8*rxt(k,296)*y(k,72) + .450_r8*rxt(k,284)*y(k,135) &
                      + .150_r8*rxt(k,268)*y(k,139)
         mat(k,632) = .550_r8*rxt(k,302)*y(k,88) + .600_r8*rxt(k,303)*y(k,90) &
                      + .600_r8*rxt(k,299)*y(k,128) + 1.200_r8*rxt(k,300)*y(k,129)
         mat(k,653) = .250_r8*rxt(k,277)*y(k,88) + .250_r8*rxt(k,276)*y(k,90) &
                      + .250_r8*rxt(k,273)*y(k,128) + .880_r8*rxt(k,274)*y(k,129)
         mat(k,671) = rxt(k,286)*y(k,88) + rxt(k,287)*y(k,90) + rxt(k,282)*y(k,128) &
                      + 2.000_r8*rxt(k,283)*y(k,129) + .450_r8*rxt(k,284)*y(k,132) &
                      + 4.000_r8*rxt(k,285)*y(k,135)
         mat(k,1077) = (rxt(k,222)+rxt(k,223))*y(k,36)
         mat(k,1162) = mat(k,1162) + .500_r8*rxt(k,242)*y(k,33) + rxt(k,218)*y(k,34) &
                      + .300_r8*rxt(k,219)*y(k,35) + .800_r8*rxt(k,247)*y(k,52) &
                      + .500_r8*rxt(k,288)*y(k,76) + rxt(k,249)*y(k,103)
         mat(k,345) = rxt(k,265)*y(k,88)
         mat(k,497) = rxt(k,269)*y(k,88) + .800_r8*rxt(k,267)*y(k,129) &
                      + .150_r8*rxt(k,268)*y(k,132)
         mat(k,582) = .250_r8*rxt(k,312)*y(k,88) + .250_r8*rxt(k,309)*y(k,128) &
                      + .800_r8*rxt(k,310)*y(k,129)
         mat(k,292) = -(rxt(k,196)*y(k,38) + rxt(k,197)*y(k,137) + rxt(k,207)*y(k,136))
         mat(k,839) = -rxt(k,196)*y(k,26)
         mat(k,1125) = -rxt(k,197)*y(k,26)
         mat(k,1070) = -rxt(k,207)*y(k,26)
         mat(k,77) = -(rxt(k,198)*y(k,137))
         mat(k,1095) = -rxt(k,198)*y(k,27)
         mat(k,550) = -(rxt(k,235)*y(k,90) + rxt(k,236)*y(k,137))
         mat(k,1194) = -rxt(k,235)*y(k,28)
         mat(k,1150) = -rxt(k,236)*y(k,28)
         mat(k,112) = rxt(k,231)*y(k,137)
         mat(k,149) = .500_r8*rxt(k,232)*y(k,137)
         mat(k,531) = .500_r8*rxt(k,254)*y(k,98)
         mat(k,710) = .100_r8*rxt(k,289)*y(k,98)
         mat(k,1262) = rxt(k,230)*y(k,126) + .270_r8*rxt(k,257)*y(k,127) + rxt(k,265) &
                      *y(k,138)
         mat(k,992) = .500_r8*rxt(k,254)*y(k,13) + .100_r8*rxt(k,289)*y(k,77)
         mat(k,370) = rxt(k,230)*y(k,88) + 3.200_r8*rxt(k,227)*y(k,126) &
                      + .800_r8*rxt(k,228)*y(k,129)
         mat(k,391) = .270_r8*rxt(k,257)*y(k,88)
         mat(k,744) = .800_r8*rxt(k,228)*y(k,126)
         mat(k,1150) = mat(k,1150) + rxt(k,231)*y(k,10) + .500_r8*rxt(k,232)*y(k,11)
         mat(k,344) = rxt(k,265)*y(k,88)
         mat(k,168) = -(rxt(k,199)*y(k,38) + rxt(k,200)*y(k,137))
         mat(k,834) = -rxt(k,199)*y(k,29)
         mat(k,1108) = -rxt(k,200)*y(k,29)
         mat(k,360) = -(rxt(k,272)*y(k,137))
         mat(k,1133) = -rxt(k,272)*y(k,30)
         mat(k,1253) = .820_r8*rxt(k,257)*y(k,127)
         mat(k,309) = .100_r8*rxt(k,317)*y(k,137)
         mat(k,388) = .820_r8*rxt(k,257)*y(k,88) + .820_r8*rxt(k,255)*y(k,129)
         mat(k,737) = .820_r8*rxt(k,255)*y(k,127)
         mat(k,1133) = mat(k,1133) + .100_r8*rxt(k,317)*y(k,122)
         mat(k,591) = -(rxt(k,260)*y(k,90) + rxt(k,261)*y(k,137))
         mat(k,1197) = -rxt(k,260)*y(k,31)
         mat(k,1153) = -rxt(k,261)*y(k,31)
         mat(k,504) = rxt(k,262)*y(k,137)
         mat(k,561) = .880_r8*rxt(k,279)*y(k,98)
         mat(k,711) = .500_r8*rxt(k,289)*y(k,98)
         mat(k,1265) = .020_r8*rxt(k,302)*y(k,133) + .250_r8*rxt(k,277)*y(k,134) &
                      + .250_r8*rxt(k,312)*y(k,140)
         mat(k,1197) = mat(k,1197) + .250_r8*rxt(k,276)*y(k,134) + .250_r8*rxt(k,313) &
                      *y(k,140)
         mat(k,195) = rxt(k,263)*y(k,137)
         mat(k,995) = .880_r8*rxt(k,279)*y(k,74) + .500_r8*rxt(k,289)*y(k,77)
         mat(k,691) = .250_r8*rxt(k,273)*y(k,134) + .250_r8*rxt(k,309)*y(k,140)
         mat(k,747) = .240_r8*rxt(k,274)*y(k,134) + .500_r8*rxt(k,267)*y(k,139) &
                      + .100_r8*rxt(k,310)*y(k,140)
         mat(k,625) = .020_r8*rxt(k,302)*y(k,88)
         mat(k,648) = .250_r8*rxt(k,277)*y(k,88) + .250_r8*rxt(k,276)*y(k,90) &
                      + .250_r8*rxt(k,273)*y(k,128) + .240_r8*rxt(k,274)*y(k,129)
         mat(k,1153) = mat(k,1153) + rxt(k,262)*y(k,69) + rxt(k,263)*y(k,91)
         mat(k,494) = .500_r8*rxt(k,267)*y(k,129)
         mat(k,579) = .250_r8*rxt(k,312)*y(k,88) + .250_r8*rxt(k,313)*y(k,90) &
                      + .250_r8*rxt(k,309)*y(k,128) + .100_r8*rxt(k,310)*y(k,129)
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
         mat(k,318) = -(rxt(k,241)*y(k,137))
         mat(k,1128) = -rxt(k,241)*y(k,32)
         mat(k,523) = .120_r8*rxt(k,254)*y(k,98)
         mat(k,978) = .120_r8*rxt(k,254)*y(k,13)
         mat(k,683) = .100_r8*rxt(k,238)*y(k,129) + .150_r8*rxt(k,239)*y(k,132)
         mat(k,734) = .100_r8*rxt(k,238)*y(k,128)
         mat(k,910) = .150_r8*rxt(k,239)*y(k,128) + .150_r8*rxt(k,284)*y(k,135)
         mat(k,663) = .150_r8*rxt(k,284)*y(k,132)
         mat(k,301) = -(rxt(k,242)*y(k,137))
         mat(k,1126) = -rxt(k,242)*y(k,33)
         mat(k,682) = .400_r8*rxt(k,239)*y(k,132)
         mat(k,909) = .400_r8*rxt(k,239)*y(k,128) + .400_r8*rxt(k,284)*y(k,135)
         mat(k,662) = .400_r8*rxt(k,284)*y(k,132)
         mat(k,248) = -(rxt(k,218)*y(k,137))
         mat(k,1120) = -rxt(k,218)*y(k,34)
         mat(k,368) = .300_r8*rxt(k,228)*y(k,129)
         mat(k,733) = .300_r8*rxt(k,228)*y(k,126) + 2.000_r8*rxt(k,215)*y(k,129) &
                      + .250_r8*rxt(k,300)*y(k,133) + .250_r8*rxt(k,274)*y(k,134) &
                      + .500_r8*rxt(k,267)*y(k,139) + .300_r8*rxt(k,310)*y(k,140)
         mat(k,617) = .250_r8*rxt(k,300)*y(k,129)
         mat(k,642) = .250_r8*rxt(k,274)*y(k,129)
         mat(k,491) = .500_r8*rxt(k,267)*y(k,129)
         mat(k,572) = .300_r8*rxt(k,310)*y(k,129)
         mat(k,200) = -(rxt(k,219)*y(k,137))
         mat(k,1113) = -rxt(k,219)*y(k,35)
         mat(k,732) = rxt(k,216)*y(k,132)
         mat(k,904) = rxt(k,216)*y(k,129)
         mat(k,477) = -(rxt(k,138)*y(k,38) + rxt(k,220)*y(k,137) + (rxt(k,221) &
                      + rxt(k,222) + rxt(k,223)) * y(k,136))
         mat(k,845) = -rxt(k,138)*y(k,36)
         mat(k,1144) = -rxt(k,220)*y(k,36)
         mat(k,1072) = -(rxt(k,221) + rxt(k,222) + rxt(k,223)) * y(k,36)
         mat(k,526) = .100_r8*rxt(k,254)*y(k,98)
         mat(k,986) = .100_r8*rxt(k,254)*y(k,13)
         mat(k,182) = -(rxt(k,190)*y(k,136) + rxt(k,201)*y(k,38) + rxt(k,202)*y(k,137))
         mat(k,1068) = -rxt(k,190)*y(k,37)
         mat(k,835) = -rxt(k,201)*y(k,37)
         mat(k,1110) = -rxt(k,202)*y(k,37)
         mat(k,856) = -(rxt(k,137)*y(k,25) + rxt(k,138)*y(k,36) + rxt(k,139)*y(k,55) &
                      + rxt(k,140)*y(k,57) + (rxt(k,141) + rxt(k,142)) * y(k,132) &
                      + rxt(k,143)*y(k,98) + rxt(k,150)*y(k,42) + rxt(k,159)*y(k,68) &
                      + rxt(k,194)*y(k,24) + rxt(k,196)*y(k,26) + rxt(k,199)*y(k,29) &
                      + rxt(k,201)*y(k,37) + rxt(k,233)*y(k,12))
         mat(k,786) = -rxt(k,137)*y(k,38)
         mat(k,484) = -rxt(k,138)*y(k,38)
         mat(k,517) = -rxt(k,139)*y(k,38)
         mat(k,278) = -rxt(k,140)*y(k,38)
         mat(k,941) = -(rxt(k,141) + rxt(k,142)) * y(k,38)
         mat(k,1007) = -rxt(k,143)*y(k,38)
         mat(k,450) = -rxt(k,150)*y(k,38)
         mat(k,403) = -rxt(k,159)*y(k,38)
         mat(k,230) = -rxt(k,194)*y(k,38)
         mat(k,296) = -rxt(k,196)*y(k,38)
         mat(k,172) = -rxt(k,199)*y(k,38)
         mat(k,185) = -rxt(k,201)*y(k,38)
         mat(k,133) = -rxt(k,233)*y(k,38)
         mat(k,963) = rxt(k,178)*y(k,41)
         mat(k,49) = 4.000_r8*rxt(k,162)*y(k,136)
         mat(k,75) = rxt(k,163)*y(k,136)
         mat(k,55) = 3.000_r8*rxt(k,164)*y(k,136)
         mat(k,58) = 3.000_r8*rxt(k,165)*y(k,136)
         mat(k,61) = 2.000_r8*rxt(k,166)*y(k,136)
         mat(k,64) = rxt(k,167)*y(k,136)
         mat(k,67) = 2.000_r8*rxt(k,168)*y(k,136)
         mat(k,78) = 3.000_r8*rxt(k,198)*y(k,137)
         mat(k,172) = mat(k,172) + rxt(k,200)*y(k,137)
         mat(k,1234) = rxt(k,178)*y(k,6) + (4.000_r8*rxt(k,145)+2.000_r8*rxt(k,147)) &
                      *y(k,41) + rxt(k,149)*y(k,88) + rxt(k,154)*y(k,97) + rxt(k,326) &
                      *y(k,110) + rxt(k,144)*y(k,129) + rxt(k,155)*y(k,137)
         mat(k,90) = 2.000_r8*rxt(k,208)*y(k,136) + 2.000_r8*rxt(k,203)*y(k,137)
         mat(k,94) = rxt(k,209)*y(k,136) + rxt(k,204)*y(k,137)
         mat(k,104) = rxt(k,210)*y(k,136) + rxt(k,205)*y(k,137)
         mat(k,821) = rxt(k,157)*y(k,97) + rxt(k,169)*y(k,136) + rxt(k,158)*y(k,137)
         mat(k,1276) = rxt(k,149)*y(k,41)
         mat(k,885) = rxt(k,154)*y(k,41) + rxt(k,157)*y(k,63)
         mat(k,605) = rxt(k,326)*y(k,41)
         mat(k,757) = rxt(k,144)*y(k,41)
         mat(k,1080) = 4.000_r8*rxt(k,162)*y(k,16) + rxt(k,163)*y(k,17) &
                      + 3.000_r8*rxt(k,164)*y(k,19) + 3.000_r8*rxt(k,165)*y(k,20) &
                      + 2.000_r8*rxt(k,166)*y(k,21) + rxt(k,167)*y(k,22) &
                      + 2.000_r8*rxt(k,168)*y(k,23) + 2.000_r8*rxt(k,208)*y(k,60) &
                      + rxt(k,209)*y(k,61) + rxt(k,210)*y(k,62) + rxt(k,169)*y(k,63)
         mat(k,1165) = 3.000_r8*rxt(k,198)*y(k,27) + rxt(k,200)*y(k,29) + rxt(k,155) &
                      *y(k,41) + 2.000_r8*rxt(k,203)*y(k,60) + rxt(k,204)*y(k,61) &
                      + rxt(k,205)*y(k,62) + rxt(k,158)*y(k,63)
         mat(k,832) = rxt(k,150)*y(k,42)
         mat(k,1222) = 2.000_r8*rxt(k,146)*y(k,41)
         mat(k,445) = rxt(k,150)*y(k,38) + (rxt(k,349)+rxt(k,354)+rxt(k,359))*y(k,63)
         mat(k,812) = (rxt(k,349)+rxt(k,354)+rxt(k,359))*y(k,42) + (rxt(k,344) &
                       +rxt(k,350)+rxt(k,355))*y(k,68)
         mat(k,400) = (rxt(k,344)+rxt(k,350)+rxt(k,355))*y(k,63)
         mat(k,1221) = 2.000_r8*rxt(k,171)*y(k,41)
         mat(k,1243) = -(rxt(k,144)*y(k,129) + (4._r8*rxt(k,145) + 4._r8*rxt(k,146) &
                      + 4._r8*rxt(k,147) + 4._r8*rxt(k,171)) * y(k,41) + rxt(k,148) &
                      *y(k,132) + rxt(k,149)*y(k,88) + rxt(k,151)*y(k,89) + rxt(k,154) &
                      *y(k,97) + (rxt(k,155) + rxt(k,156)) * y(k,137) + (rxt(k,177) &
                      + rxt(k,178) + rxt(k,179)) * y(k,6) + rxt(k,326)*y(k,110))
         mat(k,765) = -rxt(k,144)*y(k,41)
         mat(k,950) = -rxt(k,148)*y(k,41)
         mat(k,1285) = -rxt(k,149)*y(k,41)
         mat(k,1052) = -rxt(k,151)*y(k,41)
         mat(k,894) = -rxt(k,154)*y(k,41)
         mat(k,1174) = -(rxt(k,155) + rxt(k,156)) * y(k,41)
         mat(k,972) = -(rxt(k,177) + rxt(k,178) + rxt(k,179)) * y(k,41)
         mat(k,612) = -rxt(k,326)*y(k,41)
         mat(k,865) = rxt(k,159)*y(k,68) + rxt(k,143)*y(k,98) + rxt(k,142)*y(k,132)
         mat(k,455) = rxt(k,152)*y(k,97)
         mat(k,830) = rxt(k,170)*y(k,136)
         mat(k,406) = rxt(k,159)*y(k,38) + rxt(k,160)*y(k,97) + rxt(k,161)*y(k,137)
         mat(k,894) = mat(k,894) + rxt(k,152)*y(k,42) + rxt(k,160)*y(k,68)
         mat(k,1016) = rxt(k,143)*y(k,38)
         mat(k,156) = rxt(k,331)*y(k,110)
         mat(k,612) = mat(k,612) + rxt(k,331)*y(k,100)
         mat(k,950) = mat(k,950) + rxt(k,142)*y(k,38)
         mat(k,1089) = rxt(k,170)*y(k,63)
         mat(k,1174) = mat(k,1174) + rxt(k,161)*y(k,68)
         mat(k,448) = -(rxt(k,150)*y(k,38) + rxt(k,152)*y(k,97) + rxt(k,153)*y(k,137) &
                      + (rxt(k,349) + rxt(k,354) + rxt(k,359)) * y(k,63))
         mat(k,843) = -rxt(k,150)*y(k,42)
         mat(k,877) = -rxt(k,152)*y(k,42)
         mat(k,1141) = -rxt(k,153)*y(k,42)
         mat(k,817) = -(rxt(k,349) + rxt(k,354) + rxt(k,359)) * y(k,42)
         mat(k,1227) = rxt(k,151)*y(k,89)
         mat(k,1030) = rxt(k,151)*y(k,41)
         mat(k,510) = -(rxt(k,224)*y(k,137))
         mat(k,1147) = -rxt(k,224)*y(k,44)
         mat(k,797) = rxt(k,173)*y(k,25)
         mat(k,222) = .630_r8*rxt(k,226)*y(k,98)
         mat(k,528) = .560_r8*rxt(k,254)*y(k,98)
         mat(k,780) = rxt(k,173)*y(k,4) + rxt(k,137)*y(k,38) + rxt(k,211)*y(k,90) &
                      + rxt(k,212)*y(k,97) + rxt(k,213)*y(k,137)
         mat(k,169) = rxt(k,199)*y(k,38)
         mat(k,590) = rxt(k,260)*y(k,90) + rxt(k,261)*y(k,137)
         mat(k,846) = rxt(k,137)*y(k,25) + rxt(k,199)*y(k,29)
         mat(k,337) = rxt(k,248)*y(k,137)
         mat(k,431) = .620_r8*rxt(k,304)*y(k,98)
         mat(k,559) = .650_r8*rxt(k,279)*y(k,98)
         mat(k,708) = .560_r8*rxt(k,289)*y(k,98)
         mat(k,1261) = .220_r8*rxt(k,277)*y(k,134) + .250_r8*rxt(k,312)*y(k,140)
         mat(k,1192) = rxt(k,211)*y(k,25) + rxt(k,260)*y(k,31) + .220_r8*rxt(k,276) &
                      *y(k,134) + .500_r8*rxt(k,313)*y(k,140)
         mat(k,878) = rxt(k,212)*y(k,25) + rxt(k,320)*y(k,101)
         mat(k,989) = .630_r8*rxt(k,226)*y(k,9) + .560_r8*rxt(k,254)*y(k,13) &
                      + .620_r8*rxt(k,304)*y(k,71) + .650_r8*rxt(k,279)*y(k,74) &
                      + .560_r8*rxt(k,289)*y(k,77)
         mat(k,163) = rxt(k,320)*y(k,97) + rxt(k,321)*y(k,137)
         mat(k,688) = .220_r8*rxt(k,273)*y(k,134) + .250_r8*rxt(k,309)*y(k,140)
         mat(k,743) = .110_r8*rxt(k,274)*y(k,134) + .200_r8*rxt(k,310)*y(k,140)
         mat(k,646) = .220_r8*rxt(k,277)*y(k,88) + .220_r8*rxt(k,276)*y(k,90) &
                      + .220_r8*rxt(k,273)*y(k,128) + .110_r8*rxt(k,274)*y(k,129)
         mat(k,1147) = mat(k,1147) + rxt(k,213)*y(k,25) + rxt(k,261)*y(k,31) &
                      + rxt(k,248)*y(k,53) + rxt(k,321)*y(k,101)
         mat(k,577) = .250_r8*rxt(k,312)*y(k,88) + .500_r8*rxt(k,313)*y(k,90) &
                      + .250_r8*rxt(k,309)*y(k,128) + .200_r8*rxt(k,310)*y(k,129)
         mat(k,524) = .200_r8*rxt(k,254)*y(k,98)
         mat(k,319) = rxt(k,241)*y(k,137)
         mat(k,302) = .500_r8*rxt(k,242)*y(k,137)
         mat(k,509) = rxt(k,224)*y(k,137)
         mat(k,457) = .800_r8*rxt(k,247)*y(k,137)
         mat(k,335) = rxt(k,248)*y(k,137)
         mat(k,284) = .500_r8*rxt(k,288)*y(k,137)
         mat(k,707) = .100_r8*rxt(k,289)*y(k,98)
         mat(k,1249) = rxt(k,240)*y(k,128)
         mat(k,979) = .200_r8*rxt(k,254)*y(k,13) + .100_r8*rxt(k,289)*y(k,77)
         mat(k,684) = rxt(k,240)*y(k,88) + 4.000_r8*rxt(k,237)*y(k,128) &
                      + .900_r8*rxt(k,238)*y(k,129) + 2.000_r8*rxt(k,282)*y(k,135) &
                      + rxt(k,309)*y(k,140)
         mat(k,735) = .900_r8*rxt(k,238)*y(k,128) + rxt(k,283)*y(k,135)
         mat(k,911) = .450_r8*rxt(k,284)*y(k,135)
         mat(k,664) = 2.000_r8*rxt(k,282)*y(k,128) + rxt(k,283)*y(k,129) &
                      + .450_r8*rxt(k,284)*y(k,132) + 4.000_r8*rxt(k,285)*y(k,135)
         mat(k,1129) = rxt(k,241)*y(k,32) + .500_r8*rxt(k,242)*y(k,33) + rxt(k,224) &
                      *y(k,44) + .800_r8*rxt(k,247)*y(k,52) + rxt(k,248)*y(k,53) &
                      + .500_r8*rxt(k,288)*y(k,76)
         mat(k,573) = rxt(k,309)*y(k,128)
         mat(k,136) = -(rxt(k,318)*y(k,90) + (rxt(k,319) + rxt(k,333)) * y(k,137))
         mat(k,1178) = -rxt(k,318)*y(k,46)
         mat(k,1103) = -(rxt(k,319) + rxt(k,333)) * y(k,46)
         mat(k,326) = rxt(k,243)*y(k,132)
         mat(k,897) = rxt(k,243)*y(k,131)
         mat(k,459) = -(rxt(k,247)*y(k,137))
         mat(k,1142) = -rxt(k,247)*y(k,52)
         mat(k,1257) = .020_r8*rxt(k,302)*y(k,133) + .530_r8*rxt(k,277)*y(k,134) &
                      + .250_r8*rxt(k,312)*y(k,140)
         mat(k,1188) = .530_r8*rxt(k,276)*y(k,134) + .250_r8*rxt(k,313)*y(k,140)
         mat(k,686) = .530_r8*rxt(k,273)*y(k,134) + .250_r8*rxt(k,309)*y(k,140)
         mat(k,740) = .260_r8*rxt(k,274)*y(k,134) + .100_r8*rxt(k,310)*y(k,140)
         mat(k,619) = .020_r8*rxt(k,302)*y(k,88)
         mat(k,643) = .530_r8*rxt(k,277)*y(k,88) + .530_r8*rxt(k,276)*y(k,90) &
                      + .530_r8*rxt(k,273)*y(k,128) + .260_r8*rxt(k,274)*y(k,129)
         mat(k,575) = .250_r8*rxt(k,312)*y(k,88) + .250_r8*rxt(k,313)*y(k,90) &
                      + .250_r8*rxt(k,309)*y(k,128) + .100_r8*rxt(k,310)*y(k,129)
         mat(k,336) = -(rxt(k,248)*y(k,137))
         mat(k,1131) = -rxt(k,248)*y(k,53)
         mat(k,458) = .200_r8*rxt(k,247)*y(k,137)
         mat(k,1251) = .020_r8*rxt(k,302)*y(k,133) + .250_r8*rxt(k,312)*y(k,140)
         mat(k,1182) = .250_r8*rxt(k,313)*y(k,140)
         mat(k,685) = .250_r8*rxt(k,309)*y(k,140)
         mat(k,736) = .100_r8*rxt(k,310)*y(k,140)
         mat(k,618) = .020_r8*rxt(k,302)*y(k,88)
         mat(k,1131) = mat(k,1131) + .200_r8*rxt(k,247)*y(k,52)
         mat(k,574) = .250_r8*rxt(k,312)*y(k,88) + .250_r8*rxt(k,313)*y(k,90) &
                      + .250_r8*rxt(k,309)*y(k,128) + .100_r8*rxt(k,310)*y(k,129)
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
         mat(k,769) = -((rxt(k,97) + rxt(k,98) + rxt(k,99)) * y(k,132) + rxt(k,103) &
                      *y(k,98))
         mat(k,937) = -(rxt(k,97) + rxt(k,98) + rxt(k,99)) * y(k,54)
         mat(k,1003) = -rxt(k,103)*y(k,54)
         mat(k,782) = rxt(k,213)*y(k,137)
         mat(k,481) = rxt(k,222)*y(k,136)
         mat(k,852) = rxt(k,139)*y(k,55)
         mat(k,515) = rxt(k,139)*y(k,38) + rxt(k,95)*y(k,97) + rxt(k,86)*y(k,136) &
                      + rxt(k,104)*y(k,137)
         mat(k,409) = rxt(k,193)*y(k,136)
         mat(k,818) = rxt(k,170)*y(k,136)
         mat(k,214) = rxt(k,125)*y(k,137)
         mat(k,881) = rxt(k,95)*y(k,55) + rxt(k,107)*y(k,137)
         mat(k,165) = rxt(k,321)*y(k,137)
         mat(k,236) = rxt(k,327)*y(k,137)
         mat(k,603) = rxt(k,332)*y(k,137)
         mat(k,1076) = rxt(k,222)*y(k,36) + rxt(k,86)*y(k,55) + rxt(k,193)*y(k,59) &
                      + rxt(k,170)*y(k,63)
         mat(k,1161) = rxt(k,213)*y(k,25) + rxt(k,104)*y(k,55) + rxt(k,125)*y(k,78) &
                      + rxt(k,107)*y(k,97) + rxt(k,321)*y(k,101) + rxt(k,327)*y(k,108) &
                      + rxt(k,332)*y(k,110)
         mat(k,514) = -(rxt(k,86)*y(k,136) + rxt(k,95)*y(k,97) + rxt(k,104)*y(k,137) &
                      + rxt(k,139)*y(k,38))
         mat(k,1074) = -rxt(k,86)*y(k,55)
         mat(k,879) = -rxt(k,95)*y(k,55)
         mat(k,1148) = -rxt(k,104)*y(k,55)
         mat(k,847) = -rxt(k,139)*y(k,55)
         mat(k,479) = rxt(k,223)*y(k,136)
         mat(k,768) = rxt(k,97)*y(k,132)
         mat(k,926) = rxt(k,97)*y(k,54)
         mat(k,1074) = mat(k,1074) + rxt(k,223)*y(k,36)
         mat(k,42) = -(rxt(k,191)*y(k,136))
         mat(k,1055) = -rxt(k,191)*y(k,56)
         mat(k,276) = -(rxt(k,96)*y(k,97) + rxt(k,105)*y(k,137) + rxt(k,140)*y(k,38))
         mat(k,871) = -rxt(k,96)*y(k,57)
         mat(k,1123) = -rxt(k,105)*y(k,57)
         mat(k,838) = -rxt(k,140)*y(k,57)
         mat(k,908) = 2.000_r8*rxt(k,111)*y(k,132)
         mat(k,1123) = mat(k,1123) + 2.000_r8*rxt(k,110)*y(k,137)
         mat(k,115) = rxt(k,334)*y(k,141)
         mat(k,1288) = rxt(k,334)*y(k,112)
         mat(k,408) = -(rxt(k,186)*y(k,97) + rxt(k,187)*y(k,137) + (rxt(k,192) &
                      + rxt(k,193)) * y(k,136))
         mat(k,874) = -rxt(k,186)*y(k,59)
         mat(k,1137) = -rxt(k,187)*y(k,59)
         mat(k,1071) = -(rxt(k,192) + rxt(k,193)) * y(k,59)
         mat(k,796) = rxt(k,173)*y(k,25) + rxt(k,174)*y(k,132)
         mat(k,778) = rxt(k,173)*y(k,4)
         mat(k,920) = rxt(k,174)*y(k,4)
         mat(k,89) = -(rxt(k,203)*y(k,137) + rxt(k,208)*y(k,136))
         mat(k,1096) = -rxt(k,203)*y(k,60)
         mat(k,1064) = -rxt(k,208)*y(k,60)
         mat(k,93) = -(rxt(k,204)*y(k,137) + rxt(k,209)*y(k,136))
         mat(k,1097) = -rxt(k,204)*y(k,61)
         mat(k,1065) = -rxt(k,209)*y(k,61)
         mat(k,103) = -(rxt(k,205)*y(k,137) + rxt(k,210)*y(k,136))
         mat(k,1099) = -rxt(k,205)*y(k,62)
         mat(k,1067) = -rxt(k,210)*y(k,62)
         mat(k,820) = -(rxt(k,157)*y(k,97) + rxt(k,158)*y(k,137) + (rxt(k,169) &
                      + rxt(k,170)) * y(k,136) + (rxt(k,344) + rxt(k,350) + rxt(k,355) &
                      ) * y(k,68) + (rxt(k,349) + rxt(k,354) + rxt(k,359)) * y(k,42) &
                      + (rxt(k,351) + rxt(k,356)) * y(k,67))
         mat(k,884) = -rxt(k,157)*y(k,63)
         mat(k,1164) = -rxt(k,158)*y(k,63)
         mat(k,1079) = -(rxt(k,169) + rxt(k,170)) * y(k,63)
         mat(k,402) = -(rxt(k,344) + rxt(k,350) + rxt(k,355)) * y(k,63)
         mat(k,449) = -(rxt(k,349) + rxt(k,354) + rxt(k,359)) * y(k,63)
         mat(k,354) = -(rxt(k,351) + rxt(k,356)) * y(k,63)
         mat(k,132) = rxt(k,233)*y(k,38)
         mat(k,229) = rxt(k,194)*y(k,38)
         mat(k,785) = rxt(k,137)*y(k,38)
         mat(k,295) = rxt(k,196)*y(k,38)
         mat(k,171) = 2.000_r8*rxt(k,199)*y(k,38)
         mat(k,483) = rxt(k,138)*y(k,38)
         mat(k,184) = rxt(k,201)*y(k,38)
         mat(k,855) = rxt(k,233)*y(k,12) + rxt(k,194)*y(k,24) + rxt(k,137)*y(k,25) &
                      + rxt(k,196)*y(k,26) + 2.000_r8*rxt(k,199)*y(k,29) + rxt(k,138) &
                      *y(k,36) + rxt(k,201)*y(k,37) + rxt(k,139)*y(k,55) + rxt(k,140) &
                      *y(k,57) + rxt(k,159)*y(k,68) + rxt(k,141)*y(k,132)
         mat(k,1233) = rxt(k,156)*y(k,137)
         mat(k,516) = rxt(k,139)*y(k,38)
         mat(k,277) = rxt(k,140)*y(k,38)
         mat(k,402) = mat(k,402) + rxt(k,159)*y(k,38)
         mat(k,940) = rxt(k,141)*y(k,38)
         mat(k,1164) = mat(k,1164) + rxt(k,156)*y(k,41)
         mat(k,416) = -(rxt(k,134)*y(k,137))
         mat(k,1138) = -rxt(k,134)*y(k,65)
         mat(k,779) = rxt(k,211)*y(k,90)
         mat(k,548) = rxt(k,235)*y(k,90)
         mat(k,589) = rxt(k,260)*y(k,90)
         mat(k,447) = (rxt(k,349)+rxt(k,354)+rxt(k,359))*y(k,63)
         mat(k,137) = rxt(k,318)*y(k,90)
         mat(k,816) = (rxt(k,349)+rxt(k,354)+rxt(k,359))*y(k,42)
         mat(k,1028) = rxt(k,133)*y(k,137)
         mat(k,1185) = rxt(k,211)*y(k,25) + rxt(k,235)*y(k,28) + rxt(k,260)*y(k,31) &
                      + rxt(k,318)*y(k,46)
         mat(k,1138) = mat(k,1138) + rxt(k,133)*y(k,89)
         mat(k,176) = -(rxt(k,112)*y(k,137))
         mat(k,1109) = -rxt(k,112)*y(k,66)
         mat(k,1021) = rxt(k,131)*y(k,132)
         mat(k,902) = rxt(k,131)*y(k,89)
         mat(k,352) = -(rxt(k,188)*y(k,97) + (rxt(k,351) + rxt(k,356)) * y(k,63))
         mat(k,872) = -rxt(k,188)*y(k,67)
         mat(k,814) = -(rxt(k,351) + rxt(k,356)) * y(k,67)
         mat(k,956) = rxt(k,180)*y(k,132)
         mat(k,914) = rxt(k,180)*y(k,6)
         mat(k,401) = -(rxt(k,159)*y(k,38) + rxt(k,160)*y(k,97) + rxt(k,161)*y(k,137) &
                      + (rxt(k,344) + rxt(k,350) + rxt(k,355)) * y(k,63))
         mat(k,842) = -rxt(k,159)*y(k,68)
         mat(k,873) = -rxt(k,160)*y(k,68)
         mat(k,1136) = -rxt(k,161)*y(k,68)
         mat(k,815) = -(rxt(k,344) + rxt(k,350) + rxt(k,355)) * y(k,68)
         mat(k,1225) = rxt(k,148)*y(k,132)
         mat(k,446) = rxt(k,153)*y(k,137)
         mat(k,919) = rxt(k,148)*y(k,41)
         mat(k,1136) = mat(k,1136) + rxt(k,153)*y(k,42)
         mat(k,503) = -(rxt(k,262)*y(k,137))
         mat(k,1146) = -rxt(k,262)*y(k,69)
         mat(k,285) = .500_r8*rxt(k,288)*y(k,137)
         mat(k,1260) = .020_r8*rxt(k,302)*y(k,133) + .220_r8*rxt(k,277)*y(k,134) &
                      + .250_r8*rxt(k,312)*y(k,140)
         mat(k,1191) = .220_r8*rxt(k,276)*y(k,134) + .250_r8*rxt(k,313)*y(k,140)
         mat(k,262) = .500_r8*rxt(k,266)*y(k,137)
         mat(k,687) = .220_r8*rxt(k,273)*y(k,134) + .250_r8*rxt(k,309)*y(k,140)
         mat(k,742) = .230_r8*rxt(k,274)*y(k,134) + .200_r8*rxt(k,267)*y(k,139) &
                      + .100_r8*rxt(k,310)*y(k,140)
         mat(k,621) = .020_r8*rxt(k,302)*y(k,88)
         mat(k,645) = .220_r8*rxt(k,277)*y(k,88) + .220_r8*rxt(k,276)*y(k,90) &
                      + .220_r8*rxt(k,273)*y(k,128) + .230_r8*rxt(k,274)*y(k,129)
         mat(k,1146) = mat(k,1146) + .500_r8*rxt(k,288)*y(k,76) + .500_r8*rxt(k,266) &
                      *y(k,106)
         mat(k,493) = .200_r8*rxt(k,267)*y(k,129)
         mat(k,576) = .250_r8*rxt(k,312)*y(k,88) + .250_r8*rxt(k,313)*y(k,90) &
                      + .250_r8*rxt(k,309)*y(k,128) + .100_r8*rxt(k,310)*y(k,129)
         mat(k,157) = -(rxt(k,294)*y(k,137))
         mat(k,1106) = -rxt(k,294)*y(k,70)
         mat(k,1247) = .330_r8*rxt(k,302)*y(k,133)
         mat(k,1179) = rxt(k,307)*y(k,102) + .400_r8*rxt(k,303)*y(k,133)
         mat(k,465) = rxt(k,307)*y(k,90) + rxt(k,308)*y(k,137)
         mat(k,680) = .400_r8*rxt(k,299)*y(k,133)
         mat(k,731) = .300_r8*rxt(k,300)*y(k,133)
         mat(k,615) = .330_r8*rxt(k,302)*y(k,88) + .400_r8*rxt(k,303)*y(k,90) &
                      + .400_r8*rxt(k,299)*y(k,128) + .300_r8*rxt(k,300)*y(k,129)
         mat(k,1106) = mat(k,1106) + rxt(k,308)*y(k,102)
         mat(k,429) = -(rxt(k,295)*y(k,90) + rxt(k,304)*y(k,98) + rxt(k,305)*y(k,137))
         mat(k,1187) = -rxt(k,295)*y(k,71)
         mat(k,983) = -rxt(k,304)*y(k,71)
         mat(k,1140) = -rxt(k,305)*y(k,71)
         mat(k,377) = -(rxt(k,296)*y(k,132) + rxt(k,297)*y(k,88) + rxt(k,298)*y(k,90))
         mat(k,917) = -rxt(k,296)*y(k,72)
         mat(k,1255) = -rxt(k,297)*y(k,72)
         mat(k,1184) = -rxt(k,298)*y(k,72)
         mat(k,428) = rxt(k,295)*y(k,90)
         mat(k,1184) = mat(k,1184) + rxt(k,295)*y(k,71)
         mat(k,240) = -(rxt(k,306)*y(k,137))
         mat(k,1119) = -rxt(k,306)*y(k,73)
         mat(k,906) = rxt(k,301)*y(k,133)
         mat(k,616) = rxt(k,301)*y(k,132)
         mat(k,560) = -(rxt(k,279)*y(k,98) + rxt(k,280)*y(k,137))
         mat(k,993) = -rxt(k,279)*y(k,74)
         mat(k,1151) = -rxt(k,280)*y(k,74)
         mat(k,433) = .300_r8*rxt(k,304)*y(k,98)
         mat(k,379) = .167_r8*rxt(k,297)*y(k,88) + .167_r8*rxt(k,298)*y(k,90) &
                      + .167_r8*rxt(k,296)*y(k,132)
         mat(k,1263) = .167_r8*rxt(k,297)*y(k,72) + .230_r8*rxt(k,302)*y(k,133)
         mat(k,1195) = .167_r8*rxt(k,298)*y(k,72) + .250_r8*rxt(k,303)*y(k,133)
         mat(k,993) = mat(k,993) + .300_r8*rxt(k,304)*y(k,71) + 1.122_r8*rxt(k,316) &
                      *y(k,122)
         mat(k,310) = 1.122_r8*rxt(k,316)*y(k,98)
         mat(k,689) = .250_r8*rxt(k,299)*y(k,133)
         mat(k,745) = .190_r8*rxt(k,300)*y(k,133)
         mat(k,928) = .167_r8*rxt(k,296)*y(k,72)
         mat(k,623) = .230_r8*rxt(k,302)*y(k,88) + .250_r8*rxt(k,303)*y(k,90) &
                      + .250_r8*rxt(k,299)*y(k,128) + .190_r8*rxt(k,300)*y(k,129)
         mat(k,142) = -(rxt(k,281)*y(k,137))
         mat(k,1104) = -rxt(k,281)*y(k,75)
         mat(k,900) = rxt(k,275)*y(k,134)
         mat(k,641) = rxt(k,275)*y(k,132)
         mat(k,283) = -(rxt(k,288)*y(k,137))
         mat(k,1124) = -rxt(k,288)*y(k,76)
         mat(k,1025) = rxt(k,291)*y(k,135)
         mat(k,661) = rxt(k,291)*y(k,89)
         mat(k,715) = -(rxt(k,289)*y(k,98) + rxt(k,290)*y(k,137))
         mat(k,1001) = -rxt(k,289)*y(k,77)
         mat(k,1159) = -rxt(k,290)*y(k,77)
         mat(k,436) = .200_r8*rxt(k,304)*y(k,98)
         mat(k,380) = .039_r8*rxt(k,297)*y(k,88) + .039_r8*rxt(k,298)*y(k,90) &
                      + .039_r8*rxt(k,296)*y(k,132)
         mat(k,1270) = .039_r8*rxt(k,297)*y(k,72) + .320_r8*rxt(k,302)*y(k,133)
         mat(k,1203) = .039_r8*rxt(k,298)*y(k,72) + .350_r8*rxt(k,303)*y(k,133)
         mat(k,1001) = mat(k,1001) + .200_r8*rxt(k,304)*y(k,71) + .442_r8*rxt(k,316) &
                      *y(k,122)
         mat(k,312) = .442_r8*rxt(k,316)*y(k,98)
         mat(k,696) = .350_r8*rxt(k,299)*y(k,133)
         mat(k,752) = .260_r8*rxt(k,300)*y(k,133)
         mat(k,935) = .039_r8*rxt(k,296)*y(k,72)
         mat(k,630) = .320_r8*rxt(k,302)*y(k,88) + .350_r8*rxt(k,303)*y(k,90) &
                      + .350_r8*rxt(k,299)*y(k,128) + .260_r8*rxt(k,300)*y(k,129)
         mat(k,213) = -(rxt(k,113)*y(k,88) + (rxt(k,114) + rxt(k,115) + rxt(k,116) &
                      ) * y(k,89) + rxt(k,125)*y(k,137))
         mat(k,1248) = -rxt(k,113)*y(k,78)
         mat(k,1022) = -(rxt(k,114) + rxt(k,115) + rxt(k,116)) * y(k,78)
         mat(k,1115) = -rxt(k,125)*y(k,78)
         mat(k,97) = -((rxt(k,129) + rxt(k,130)) * y(k,136))
         mat(k,1066) = -(rxt(k,129) + rxt(k,130)) * y(k,79)
         mat(k,212) = rxt(k,114)*y(k,89)
         mat(k,1019) = rxt(k,114)*y(k,78)
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
         mat(k,1020) = rxt(k,132)*y(k,90)
         mat(k,1177) = rxt(k,132)*y(k,89)
         mat(k,45) = -(rxt(k,335)*y(k,137))
         mat(k,1093) = -rxt(k,335)*y(k,84)
         mat(k,1286) = -(rxt(k,113)*y(k,78) + rxt(k,122)*y(k,90) + rxt(k,126)*y(k,132) &
                      + rxt(k,127)*y(k,98) + rxt(k,128)*y(k,97) + rxt(k,149)*y(k,41) &
                      + rxt(k,181)*y(k,6) + rxt(k,217)*y(k,129) + rxt(k,230)*y(k,126) &
                      + rxt(k,240)*y(k,128) + rxt(k,244)*y(k,131) + rxt(k,257) &
                      *y(k,127) + rxt(k,265)*y(k,138) + rxt(k,269)*y(k,139) + (rxt(k,277) &
                      + rxt(k,278)) * y(k,134) + rxt(k,286)*y(k,135) + rxt(k,297) &
                      *y(k,72) + rxt(k,302)*y(k,133) + rxt(k,312)*y(k,140))
         mat(k,219) = -rxt(k,113)*y(k,88)
         mat(k,1219) = -rxt(k,122)*y(k,88)
         mat(k,951) = -rxt(k,126)*y(k,88)
         mat(k,1017) = -rxt(k,127)*y(k,88)
         mat(k,895) = -rxt(k,128)*y(k,88)
         mat(k,1244) = -rxt(k,149)*y(k,88)
         mat(k,973) = -rxt(k,181)*y(k,88)
         mat(k,766) = -rxt(k,217)*y(k,88)
         mat(k,376) = -rxt(k,230)*y(k,88)
         mat(k,705) = -rxt(k,240)*y(k,88)
         mat(k,334) = -rxt(k,244)*y(k,88)
         mat(k,398) = -rxt(k,257)*y(k,88)
         mat(k,349) = -rxt(k,265)*y(k,88)
         mat(k,501) = -rxt(k,269)*y(k,88)
         mat(k,659) = -(rxt(k,277) + rxt(k,278)) * y(k,88)
         mat(k,678) = -rxt(k,286)*y(k,88)
         mat(k,386) = -rxt(k,297)*y(k,88)
         mat(k,639) = -rxt(k,302)*y(k,88)
         mat(k,588) = -rxt(k,312)*y(k,88)
         mat(k,219) = mat(k,219) + 2.000_r8*rxt(k,115)*y(k,89) + rxt(k,125)*y(k,137)
         mat(k,99) = 2.000_r8*rxt(k,129)*y(k,136)
         mat(k,1053) = 2.000_r8*rxt(k,115)*y(k,78) + rxt(k,118)*y(k,97) + rxt(k,328) &
                      *y(k,110)
         mat(k,895) = mat(k,895) + rxt(k,118)*y(k,89)
         mat(k,613) = rxt(k,328)*y(k,89)
         mat(k,1090) = 2.000_r8*rxt(k,129)*y(k,79)
         mat(k,1175) = rxt(k,125)*y(k,78)
         mat(k,1048) = -((rxt(k,114) + rxt(k,115) + rxt(k,116)) * y(k,78) + (rxt(k,118) &
                      + rxt(k,120)) * y(k,97) + rxt(k,119)*y(k,98) + rxt(k,131) &
                      *y(k,132) + rxt(k,132)*y(k,90) + rxt(k,133)*y(k,137) + rxt(k,151) &
                      *y(k,41) + rxt(k,182)*y(k,6) + rxt(k,251)*y(k,128) + rxt(k,291) &
                      *y(k,135) + rxt(k,328)*y(k,110))
         mat(k,216) = -(rxt(k,114) + rxt(k,115) + rxt(k,116)) * y(k,89)
         mat(k,890) = -(rxt(k,118) + rxt(k,120)) * y(k,89)
         mat(k,1012) = -rxt(k,119)*y(k,89)
         mat(k,946) = -rxt(k,131)*y(k,89)
         mat(k,1214) = -rxt(k,132)*y(k,89)
         mat(k,1170) = -rxt(k,133)*y(k,89)
         mat(k,1239) = -rxt(k,151)*y(k,89)
         mat(k,968) = -rxt(k,182)*y(k,89)
         mat(k,702) = -rxt(k,251)*y(k,89)
         mat(k,675) = -rxt(k,291)*y(k,89)
         mat(k,610) = -rxt(k,328)*y(k,89)
         mat(k,968) = mat(k,968) + rxt(k,181)*y(k,88)
         mat(k,1239) = mat(k,1239) + rxt(k,149)*y(k,88)
         mat(k,178) = rxt(k,112)*y(k,137)
         mat(k,383) = 1.206_r8*rxt(k,297)*y(k,88) + 1.206_r8*rxt(k,298)*y(k,90) &
                      + .206_r8*rxt(k,296)*y(k,132)
         mat(k,1281) = rxt(k,181)*y(k,6) + rxt(k,149)*y(k,41) + 1.206_r8*rxt(k,297) &
                      *y(k,72) + 2.000_r8*rxt(k,122)*y(k,90) + rxt(k,128)*y(k,97) &
                      + rxt(k,127)*y(k,98) + rxt(k,230)*y(k,126) + rxt(k,257)*y(k,127) &
                      + rxt(k,240)*y(k,128) + rxt(k,217)*y(k,129) + rxt(k,244) &
                      *y(k,131) + rxt(k,126)*y(k,132) + .920_r8*rxt(k,302)*y(k,133) &
                      + rxt(k,277)*y(k,134) + rxt(k,286)*y(k,135) + rxt(k,265) &
                      *y(k,138) + rxt(k,269)*y(k,139) + rxt(k,312)*y(k,140)
         mat(k,1214) = mat(k,1214) + 1.206_r8*rxt(k,298)*y(k,72) + 2.000_r8*rxt(k,122) &
                      *y(k,88) + rxt(k,123)*y(k,97) + rxt(k,307)*y(k,102) + rxt(k,315) &
                      *y(k,122) + rxt(k,121)*y(k,132) + rxt(k,303)*y(k,133) &
                      + rxt(k,276)*y(k,134) + rxt(k,287)*y(k,135) + rxt(k,124) &
                      *y(k,137) + rxt(k,313)*y(k,140)
         mat(k,198) = rxt(k,263)*y(k,137)
         mat(k,890) = mat(k,890) + rxt(k,128)*y(k,88) + rxt(k,123)*y(k,90)
         mat(k,1012) = mat(k,1012) + rxt(k,127)*y(k,88)
         mat(k,472) = rxt(k,307)*y(k,90) + .400_r8*rxt(k,308)*y(k,137)
         mat(k,315) = rxt(k,315)*y(k,90)
         mat(k,374) = rxt(k,230)*y(k,88)
         mat(k,396) = rxt(k,257)*y(k,88)
         mat(k,702) = mat(k,702) + rxt(k,240)*y(k,88)
         mat(k,761) = rxt(k,217)*y(k,88)
         mat(k,332) = rxt(k,244)*y(k,88)
         mat(k,946) = mat(k,946) + .206_r8*rxt(k,296)*y(k,72) + rxt(k,126)*y(k,88) &
                      + rxt(k,121)*y(k,90)
         mat(k,636) = .920_r8*rxt(k,302)*y(k,88) + rxt(k,303)*y(k,90)
         mat(k,656) = rxt(k,277)*y(k,88) + rxt(k,276)*y(k,90)
         mat(k,675) = mat(k,675) + rxt(k,286)*y(k,88) + rxt(k,287)*y(k,90)
         mat(k,1170) = mat(k,1170) + rxt(k,112)*y(k,66) + rxt(k,124)*y(k,90) &
                      + rxt(k,263)*y(k,91) + .400_r8*rxt(k,308)*y(k,102)
         mat(k,347) = rxt(k,265)*y(k,88)
         mat(k,499) = rxt(k,269)*y(k,88)
         mat(k,585) = rxt(k,312)*y(k,88) + rxt(k,313)*y(k,90)
         mat(k,1217) = -(rxt(k,121)*y(k,132) + rxt(k,122)*y(k,88) + rxt(k,123)*y(k,97) &
                      + rxt(k,124)*y(k,137) + rxt(k,132)*y(k,89) + rxt(k,211)*y(k,25) &
                      + rxt(k,235)*y(k,28) + rxt(k,253)*y(k,13) + rxt(k,260)*y(k,31) &
                      + rxt(k,276)*y(k,134) + rxt(k,287)*y(k,135) + rxt(k,295)*y(k,71) &
                      + rxt(k,298)*y(k,72) + rxt(k,303)*y(k,133) + rxt(k,307)*y(k,102) &
                      + rxt(k,313)*y(k,140) + rxt(k,315)*y(k,122) + rxt(k,318)*y(k,46))
         mat(k,949) = -rxt(k,121)*y(k,90)
         mat(k,1284) = -rxt(k,122)*y(k,90)
         mat(k,893) = -rxt(k,123)*y(k,90)
         mat(k,1173) = -rxt(k,124)*y(k,90)
         mat(k,1051) = -rxt(k,132)*y(k,90)
         mat(k,794) = -rxt(k,211)*y(k,90)
         mat(k,557) = -rxt(k,235)*y(k,90)
         mat(k,545) = -rxt(k,253)*y(k,90)
         mat(k,597) = -rxt(k,260)*y(k,90)
         mat(k,658) = -rxt(k,276)*y(k,90)
         mat(k,677) = -rxt(k,287)*y(k,90)
         mat(k,443) = -rxt(k,295)*y(k,90)
         mat(k,385) = -rxt(k,298)*y(k,90)
         mat(k,638) = -rxt(k,303)*y(k,90)
         mat(k,474) = -rxt(k,307)*y(k,90)
         mat(k,587) = -rxt(k,313)*y(k,90)
         mat(k,317) = -rxt(k,315)*y(k,90)
         mat(k,141) = -rxt(k,318)*y(k,90)
         mat(k,275) = rxt(k,183)*y(k,97)
         mat(k,864) = rxt(k,150)*y(k,42)
         mat(k,454) = rxt(k,150)*y(k,38) + rxt(k,152)*y(k,97) + rxt(k,153)*y(k,137)
         mat(k,419) = rxt(k,134)*y(k,137)
         mat(k,291) = .500_r8*rxt(k,288)*y(k,137)
         mat(k,1051) = mat(k,1051) + rxt(k,120)*y(k,97) + rxt(k,119)*y(k,98)
         mat(k,893) = mat(k,893) + rxt(k,183)*y(k,7) + rxt(k,152)*y(k,42) + rxt(k,120) &
                      *y(k,89)
         mat(k,1015) = rxt(k,119)*y(k,89)
         mat(k,259) = rxt(k,249)*y(k,137)
         mat(k,1173) = mat(k,1173) + rxt(k,153)*y(k,42) + rxt(k,134)*y(k,65) &
                      + .500_r8*rxt(k,288)*y(k,76) + rxt(k,249)*y(k,103)
         mat(k,194) = -(rxt(k,263)*y(k,137))
         mat(k,1112) = -rxt(k,263)*y(k,91)
         mat(k,522) = rxt(k,253)*y(k,90)
         mat(k,1180) = rxt(k,253)*y(k,13)
         mat(k,886) = -(rxt(k,92)*y(k,98) + 4._r8*rxt(k,93)*y(k,97) + rxt(k,95) &
                      *y(k,55) + rxt(k,96)*y(k,57) + rxt(k,101)*y(k,132) + rxt(k,107) &
                      *y(k,137) + (rxt(k,118) + rxt(k,120)) * y(k,89) + rxt(k,123) &
                      *y(k,90) + rxt(k,128)*y(k,88) + rxt(k,152)*y(k,42) + rxt(k,154) &
                      *y(k,41) + rxt(k,157)*y(k,63) + rxt(k,160)*y(k,68) + rxt(k,183) &
                      *y(k,7) + rxt(k,184)*y(k,6) + rxt(k,186)*y(k,59) + rxt(k,188) &
                      *y(k,67) + rxt(k,212)*y(k,25) + rxt(k,320)*y(k,101))
         mat(k,1008) = -rxt(k,92)*y(k,97)
         mat(k,518) = -rxt(k,95)*y(k,97)
         mat(k,279) = -rxt(k,96)*y(k,97)
         mat(k,942) = -rxt(k,101)*y(k,97)
         mat(k,1166) = -rxt(k,107)*y(k,97)
         mat(k,1044) = -(rxt(k,118) + rxt(k,120)) * y(k,97)
         mat(k,1210) = -rxt(k,123)*y(k,97)
         mat(k,1277) = -rxt(k,128)*y(k,97)
         mat(k,451) = -rxt(k,152)*y(k,97)
         mat(k,1235) = -rxt(k,154)*y(k,97)
         mat(k,822) = -rxt(k,157)*y(k,97)
         mat(k,404) = -rxt(k,160)*y(k,97)
         mat(k,272) = -rxt(k,183)*y(k,97)
         mat(k,964) = -rxt(k,184)*y(k,97)
         mat(k,411) = -rxt(k,186)*y(k,97)
         mat(k,356) = -rxt(k,188)*y(k,97)
         mat(k,787) = -rxt(k,212)*y(k,97)
         mat(k,166) = -rxt(k,320)*y(k,97)
         mat(k,772) = rxt(k,99)*y(k,132)
         mat(k,215) = rxt(k,113)*y(k,88) + rxt(k,114)*y(k,89)
         mat(k,1277) = mat(k,1277) + rxt(k,113)*y(k,78)
         mat(k,1044) = mat(k,1044) + rxt(k,114)*y(k,78)
         mat(k,1008) = mat(k,1008) + .765_r8*rxt(k,316)*y(k,122) + 2.000_r8*rxt(k,91) &
                      *y(k,136)
         mat(k,313) = .765_r8*rxt(k,316)*y(k,98)
         mat(k,942) = mat(k,942) + rxt(k,99)*y(k,54)
         mat(k,1081) = 2.000_r8*rxt(k,91)*y(k,98)
         mat(k,1166) = mat(k,1166) + 2.000_r8*rxt(k,109)*y(k,137)
         mat(k,1011) = -((rxt(k,90) + rxt(k,91)) * y(k,136) + rxt(k,92)*y(k,97) &
                      + rxt(k,102)*y(k,132) + rxt(k,103)*y(k,54) + rxt(k,108)*y(k,137) &
                      + rxt(k,119)*y(k,89) + rxt(k,127)*y(k,88) + rxt(k,143)*y(k,38) &
                      + rxt(k,175)*y(k,4) + rxt(k,226)*y(k,9) + rxt(k,254)*y(k,13) &
                      + rxt(k,279)*y(k,74) + rxt(k,289)*y(k,77) + rxt(k,304)*y(k,71) &
                      + rxt(k,316)*y(k,122) + rxt(k,324)*y(k,108) + rxt(k,330) &
                      *y(k,110))
         mat(k,1084) = -(rxt(k,90) + rxt(k,91)) * y(k,98)
         mat(k,889) = -rxt(k,92)*y(k,98)
         mat(k,945) = -rxt(k,102)*y(k,98)
         mat(k,774) = -rxt(k,103)*y(k,98)
         mat(k,1169) = -rxt(k,108)*y(k,98)
         mat(k,1047) = -rxt(k,119)*y(k,98)
         mat(k,1280) = -rxt(k,127)*y(k,98)
         mat(k,860) = -rxt(k,143)*y(k,98)
         mat(k,806) = -rxt(k,175)*y(k,98)
         mat(k,225) = -rxt(k,226)*y(k,98)
         mat(k,541) = -rxt(k,254)*y(k,98)
         mat(k,568) = -rxt(k,279)*y(k,98)
         mat(k,723) = -rxt(k,289)*y(k,98)
         mat(k,440) = -rxt(k,304)*y(k,98)
         mat(k,314) = -rxt(k,316)*y(k,98)
         mat(k,238) = -rxt(k,324)*y(k,98)
         mat(k,609) = -rxt(k,330)*y(k,98)
         mat(k,701) = .150_r8*rxt(k,239)*y(k,132)
         mat(k,945) = mat(k,945) + .150_r8*rxt(k,239)*y(k,128) + .150_r8*rxt(k,284) &
                      *y(k,135)
         mat(k,674) = .150_r8*rxt(k,284)*y(k,132)
         mat(k,152) = -(rxt(k,331)*y(k,110))
         mat(k,599) = -rxt(k,331)*y(k,100)
         mat(k,954) = rxt(k,177)*y(k,41)
         mat(k,1224) = rxt(k,177)*y(k,6) + 2.000_r8*rxt(k,147)*y(k,41)
         mat(k,160) = -(rxt(k,320)*y(k,97) + rxt(k,321)*y(k,137))
         mat(k,868) = -rxt(k,320)*y(k,101)
         mat(k,1107) = -rxt(k,321)*y(k,101)
         mat(k,467) = -(rxt(k,307)*y(k,90) + rxt(k,308)*y(k,137))
         mat(k,1189) = -rxt(k,307)*y(k,102)
         mat(k,1143) = -rxt(k,308)*y(k,102)
         mat(k,378) = .794_r8*rxt(k,297)*y(k,88) + .794_r8*rxt(k,298)*y(k,90) &
                      + .794_r8*rxt(k,296)*y(k,132)
         mat(k,1258) = .794_r8*rxt(k,297)*y(k,72) + .080_r8*rxt(k,302)*y(k,133) &
                      + .800_r8*rxt(k,278)*y(k,134)
         mat(k,1189) = mat(k,1189) + .794_r8*rxt(k,298)*y(k,72)
         mat(k,922) = .794_r8*rxt(k,296)*y(k,72)
         mat(k,620) = .080_r8*rxt(k,302)*y(k,88)
         mat(k,644) = .800_r8*rxt(k,278)*y(k,88)
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
         mat(k,252) = -(rxt(k,249)*y(k,137))
         mat(k,1121) = -rxt(k,249)*y(k,103)
         mat(k,1023) = rxt(k,251)*y(k,128)
         mat(k,681) = rxt(k,251)*y(k,89)
         mat(k,260) = -(rxt(k,266)*y(k,137))
         mat(k,1122) = -rxt(k,266)*y(k,106)
         mat(k,907) = rxt(k,264)*y(k,138)
         mat(k,341) = rxt(k,264)*y(k,132)
         mat(k,206) = -(rxt(k,270)*y(k,137))
         mat(k,1114) = -rxt(k,270)*y(k,107)
         mat(k,905) = .850_r8*rxt(k,268)*y(k,139)
         mat(k,490) = .850_r8*rxt(k,268)*y(k,132)
         mat(k,234) = -(rxt(k,324)*y(k,98) + rxt(k,327)*y(k,137))
         mat(k,976) = -rxt(k,324)*y(k,108)
         mat(k,1118) = -rxt(k,327)*y(k,108)
         mat(k,602) = -(rxt(k,325)*y(k,6) + rxt(k,326)*y(k,41) + rxt(k,328)*y(k,89) &
                      + rxt(k,330)*y(k,98) + rxt(k,331)*y(k,100) + rxt(k,332)*y(k,137))
         mat(k,959) = -rxt(k,325)*y(k,110)
         mat(k,1228) = -rxt(k,326)*y(k,110)
         mat(k,1034) = -rxt(k,328)*y(k,110)
         mat(k,996) = -rxt(k,330)*y(k,110)
         mat(k,154) = -rxt(k,331)*y(k,110)
         mat(k,1154) = -rxt(k,332)*y(k,110)
         mat(k,880) = rxt(k,320)*y(k,101)
         mat(k,996) = mat(k,996) + rxt(k,324)*y(k,108)
         mat(k,164) = rxt(k,320)*y(k,97)
         mat(k,235) = rxt(k,324)*y(k,98) + rxt(k,327)*y(k,137)
         mat(k,1154) = mat(k,1154) + rxt(k,327)*y(k,108)
         mat(k,422) = -(rxt(k,323)*y(k,137))
         mat(k,1139) = -rxt(k,323)*y(k,111)
         mat(k,958) = rxt(k,325)*y(k,110)
         mat(k,1226) = rxt(k,326)*y(k,110)
         mat(k,138) = rxt(k,318)*y(k,90) + (rxt(k,319)+.500_r8*rxt(k,333))*y(k,137)
         mat(k,1029) = rxt(k,328)*y(k,110)
         mat(k,1186) = rxt(k,318)*y(k,46)
         mat(k,982) = rxt(k,330)*y(k,110)
         mat(k,153) = rxt(k,331)*y(k,110)
         mat(k,162) = rxt(k,321)*y(k,137)
         mat(k,601) = rxt(k,325)*y(k,6) + rxt(k,326)*y(k,41) + rxt(k,328)*y(k,89) &
                      + rxt(k,330)*y(k,98) + rxt(k,331)*y(k,100) + rxt(k,332)*y(k,137)
         mat(k,1139) = mat(k,1139) + (rxt(k,319)+.500_r8*rxt(k,333))*y(k,46) &
                      + rxt(k,321)*y(k,101) + rxt(k,332)*y(k,110)
         mat(k,116) = -(rxt(k,334)*y(k,141))
         mat(k,1289) = -rxt(k,334)*y(k,112)
         mat(k,421) = rxt(k,323)*y(k,137)
         mat(k,1101) = rxt(k,323)*y(k,111)
         mat(k,308) = -(rxt(k,315)*y(k,90) + rxt(k,316)*y(k,98) + rxt(k,317)*y(k,137))
         mat(k,1181) = -rxt(k,315)*y(k,122)
         mat(k,977) = -rxt(k,316)*y(k,122)
         mat(k,1127) = -rxt(k,317)*y(k,122)
         mat(k,100) = -(rxt(k,314)*y(k,137))
         mat(k,1098) = -rxt(k,314)*y(k,123)
         mat(k,898) = rxt(k,311)*y(k,140)
         mat(k,571) = rxt(k,311)*y(k,132)
         mat(k,369) = -(4._r8*rxt(k,227)*y(k,126) + rxt(k,228)*y(k,129) + rxt(k,229) &
                      *y(k,132) + rxt(k,230)*y(k,88))
         mat(k,738) = -rxt(k,228)*y(k,126)
         mat(k,916) = -rxt(k,229)*y(k,126)
         mat(k,1254) = -rxt(k,230)*y(k,126)
         mat(k,148) = .500_r8*rxt(k,232)*y(k,137)
         mat(k,131) = rxt(k,233)*y(k,38) + rxt(k,234)*y(k,137)
         mat(k,841) = rxt(k,233)*y(k,12)
         mat(k,1134) = .500_r8*rxt(k,232)*y(k,11) + rxt(k,234)*y(k,12)
         mat(k,389) = -(rxt(k,255)*y(k,129) + rxt(k,256)*y(k,132) + rxt(k,257)*y(k,88))
         mat(k,739) = -rxt(k,255)*y(k,127)
         mat(k,918) = -rxt(k,256)*y(k,127)
         mat(k,1256) = -rxt(k,257)*y(k,127)
         mat(k,40) = 1.670_r8*rxt(k,293)*y(k,137)
         mat(k,190) = rxt(k,258)*y(k,137)
         mat(k,70) = rxt(k,259)*y(k,137)
         mat(k,1135) = 1.670_r8*rxt(k,293)*y(k,3) + rxt(k,258)*y(k,14) + rxt(k,259) &
                      *y(k,15)
         mat(k,695) = -(4._r8*rxt(k,237)*y(k,128) + rxt(k,238)*y(k,129) + rxt(k,239) &
                      *y(k,132) + rxt(k,240)*y(k,88) + rxt(k,251)*y(k,89) + rxt(k,273) &
                      *y(k,134) + rxt(k,299)*y(k,133) + rxt(k,309)*y(k,140))
         mat(k,751) = -rxt(k,238)*y(k,128)
         mat(k,934) = -rxt(k,239)*y(k,128)
         mat(k,1269) = -rxt(k,240)*y(k,128)
         mat(k,1036) = -rxt(k,251)*y(k,128)
         mat(k,651) = -rxt(k,273)*y(k,128)
         mat(k,629) = -rxt(k,299)*y(k,128)
         mat(k,580) = -rxt(k,309)*y(k,128)
         mat(k,551) = rxt(k,235)*y(k,90) + rxt(k,236)*y(k,137)
         mat(k,592) = rxt(k,260)*y(k,90) + rxt(k,261)*y(k,137)
         mat(k,303) = .500_r8*rxt(k,242)*y(k,137)
         mat(k,435) = .080_r8*rxt(k,304)*y(k,98)
         mat(k,564) = .100_r8*rxt(k,279)*y(k,98)
         mat(k,714) = .280_r8*rxt(k,289)*y(k,98)
         mat(k,1269) = mat(k,1269) + .530_r8*rxt(k,277)*y(k,134) + rxt(k,286)*y(k,135) &
                      + rxt(k,269)*y(k,139)
         mat(k,1202) = rxt(k,235)*y(k,28) + rxt(k,260)*y(k,31) + .530_r8*rxt(k,276) &
                      *y(k,134) + rxt(k,287)*y(k,135)
         mat(k,1000) = .080_r8*rxt(k,304)*y(k,71) + .100_r8*rxt(k,279)*y(k,74) &
                      + .280_r8*rxt(k,289)*y(k,77)
         mat(k,695) = mat(k,695) + .530_r8*rxt(k,273)*y(k,134)
         mat(k,751) = mat(k,751) + .260_r8*rxt(k,274)*y(k,134) + rxt(k,283)*y(k,135) &
                      + .300_r8*rxt(k,267)*y(k,139)
         mat(k,934) = mat(k,934) + .450_r8*rxt(k,284)*y(k,135) + .150_r8*rxt(k,268) &
                      *y(k,139)
         mat(k,651) = mat(k,651) + .530_r8*rxt(k,277)*y(k,88) + .530_r8*rxt(k,276) &
                      *y(k,90) + .530_r8*rxt(k,273)*y(k,128) + .260_r8*rxt(k,274) &
                      *y(k,129)
         mat(k,669) = rxt(k,286)*y(k,88) + rxt(k,287)*y(k,90) + rxt(k,283)*y(k,129) &
                      + .450_r8*rxt(k,284)*y(k,132) + 4.000_r8*rxt(k,285)*y(k,135)
         mat(k,1158) = rxt(k,236)*y(k,28) + rxt(k,261)*y(k,31) + .500_r8*rxt(k,242) &
                      *y(k,33)
         mat(k,495) = rxt(k,269)*y(k,88) + .300_r8*rxt(k,267)*y(k,129) &
                      + .150_r8*rxt(k,268)*y(k,132)
         mat(k,753) = -(rxt(k,144)*y(k,41) + (4._r8*rxt(k,214) + 4._r8*rxt(k,215) &
                      ) * y(k,129) + rxt(k,216)*y(k,132) + rxt(k,217)*y(k,88) &
                      + rxt(k,228)*y(k,126) + rxt(k,238)*y(k,128) + rxt(k,255) &
                      *y(k,127) + rxt(k,267)*y(k,139) + rxt(k,274)*y(k,134) + rxt(k,283) &
                      *y(k,135) + rxt(k,300)*y(k,133) + rxt(k,310)*y(k,140))
         mat(k,1229) = -rxt(k,144)*y(k,129)
         mat(k,936) = -rxt(k,216)*y(k,129)
         mat(k,1271) = -rxt(k,217)*y(k,129)
         mat(k,371) = -rxt(k,228)*y(k,129)
         mat(k,697) = -rxt(k,238)*y(k,129)
         mat(k,393) = -rxt(k,255)*y(k,129)
         mat(k,496) = -rxt(k,267)*y(k,129)
         mat(k,652) = -rxt(k,274)*y(k,129)
         mat(k,670) = -rxt(k,283)*y(k,129)
         mat(k,631) = -rxt(k,300)*y(k,129)
         mat(k,581) = -rxt(k,310)*y(k,129)
         mat(k,534) = .280_r8*rxt(k,254)*y(k,98)
         mat(k,320) = rxt(k,241)*y(k,137)
         mat(k,201) = .700_r8*rxt(k,219)*y(k,137)
         mat(k,480) = rxt(k,138)*y(k,38) + rxt(k,221)*y(k,136) + rxt(k,220)*y(k,137)
         mat(k,851) = rxt(k,138)*y(k,36)
         mat(k,437) = .050_r8*rxt(k,304)*y(k,98)
         mat(k,1271) = mat(k,1271) + rxt(k,240)*y(k,128)
         mat(k,1002) = .280_r8*rxt(k,254)*y(k,13) + .050_r8*rxt(k,304)*y(k,71)
         mat(k,697) = mat(k,697) + rxt(k,240)*y(k,88) + 4.000_r8*rxt(k,237)*y(k,128) &
                      + .900_r8*rxt(k,238)*y(k,129) + .450_r8*rxt(k,239)*y(k,132) &
                      + rxt(k,299)*y(k,133) + rxt(k,273)*y(k,134) + rxt(k,282) &
                      *y(k,135) + rxt(k,309)*y(k,140)
         mat(k,753) = mat(k,753) + .900_r8*rxt(k,238)*y(k,128)
         mat(k,936) = mat(k,936) + .450_r8*rxt(k,239)*y(k,128)
         mat(k,631) = mat(k,631) + rxt(k,299)*y(k,128)
         mat(k,652) = mat(k,652) + rxt(k,273)*y(k,128)
         mat(k,670) = mat(k,670) + rxt(k,282)*y(k,128)
         mat(k,1075) = rxt(k,221)*y(k,36)
         mat(k,1160) = rxt(k,241)*y(k,32) + .700_r8*rxt(k,219)*y(k,35) + rxt(k,220) &
                      *y(k,36)
         mat(k,581) = mat(k,581) + rxt(k,309)*y(k,128)
         mat(k,1246) = .750_r8*rxt(k,244)*y(k,131)
         mat(k,327) = .750_r8*rxt(k,244)*y(k,88)
         mat(k,328) = -(rxt(k,243)*y(k,132) + rxt(k,244)*y(k,88))
         mat(k,912) = -rxt(k,243)*y(k,131)
         mat(k,1250) = -rxt(k,244)*y(k,131)
         mat(k,221) = rxt(k,250)*y(k,137)
         mat(k,1130) = rxt(k,250)*y(k,9)
         mat(k,943) = -((rxt(k,97) + rxt(k,98) + rxt(k,99)) * y(k,54) + rxt(k,101) &
                      *y(k,97) + rxt(k,102)*y(k,98) + rxt(k,106)*y(k,137) &
                      + 4._r8*rxt(k,111)*y(k,132) + rxt(k,121)*y(k,90) + rxt(k,126) &
                      *y(k,88) + rxt(k,131)*y(k,89) + (rxt(k,141) + rxt(k,142) &
                      ) * y(k,38) + rxt(k,148)*y(k,41) + rxt(k,174)*y(k,4) + rxt(k,180) &
                      *y(k,6) + rxt(k,216)*y(k,129) + rxt(k,229)*y(k,126) + rxt(k,239) &
                      *y(k,128) + rxt(k,243)*y(k,131) + rxt(k,256)*y(k,127) + rxt(k,264) &
                      *y(k,138) + rxt(k,268)*y(k,139) + rxt(k,275)*y(k,134) + rxt(k,284) &
                      *y(k,135) + rxt(k,296)*y(k,72) + rxt(k,301)*y(k,133) + rxt(k,311) &
                      *y(k,140))
         mat(k,773) = -(rxt(k,97) + rxt(k,98) + rxt(k,99)) * y(k,132)
         mat(k,887) = -rxt(k,101)*y(k,132)
         mat(k,1009) = -rxt(k,102)*y(k,132)
         mat(k,1167) = -rxt(k,106)*y(k,132)
         mat(k,1211) = -rxt(k,121)*y(k,132)
         mat(k,1278) = -rxt(k,126)*y(k,132)
         mat(k,1045) = -rxt(k,131)*y(k,132)
         mat(k,858) = -(rxt(k,141) + rxt(k,142)) * y(k,132)
         mat(k,1236) = -rxt(k,148)*y(k,132)
         mat(k,804) = -rxt(k,174)*y(k,132)
         mat(k,965) = -rxt(k,180)*y(k,132)
         mat(k,759) = -rxt(k,216)*y(k,132)
         mat(k,373) = -rxt(k,229)*y(k,132)
         mat(k,700) = -rxt(k,239)*y(k,132)
         mat(k,331) = -rxt(k,243)*y(k,132)
         mat(k,395) = -rxt(k,256)*y(k,132)
         mat(k,346) = -rxt(k,264)*y(k,132)
         mat(k,498) = -rxt(k,268)*y(k,132)
         mat(k,655) = -rxt(k,275)*y(k,132)
         mat(k,673) = -rxt(k,284)*y(k,132)
         mat(k,382) = -rxt(k,296)*y(k,132)
         mat(k,634) = -rxt(k,301)*y(k,132)
         mat(k,584) = -rxt(k,311)*y(k,132)
         mat(k,804) = mat(k,804) + rxt(k,173)*y(k,25)
         mat(k,965) = mat(k,965) + rxt(k,185)*y(k,137)
         mat(k,224) = .130_r8*rxt(k,226)*y(k,98)
         mat(k,113) = rxt(k,231)*y(k,137)
         mat(k,540) = .280_r8*rxt(k,254)*y(k,98)
         mat(k,788) = rxt(k,173)*y(k,4) + rxt(k,137)*y(k,38) + rxt(k,211)*y(k,90) &
                      + rxt(k,212)*y(k,97)
         mat(k,297) = rxt(k,196)*y(k,38) + rxt(k,197)*y(k,137)
         mat(k,173) = rxt(k,199)*y(k,38) + rxt(k,200)*y(k,137)
         mat(k,250) = rxt(k,218)*y(k,137)
         mat(k,486) = rxt(k,222)*y(k,136)
         mat(k,858) = mat(k,858) + rxt(k,137)*y(k,25) + rxt(k,196)*y(k,26) &
                      + rxt(k,199)*y(k,29) + rxt(k,140)*y(k,57)
         mat(k,1236) = mat(k,1236) + rxt(k,144)*y(k,129) + rxt(k,155)*y(k,137)
         mat(k,512) = rxt(k,224)*y(k,137)
         mat(k,139) = .500_r8*rxt(k,333)*y(k,137)
         mat(k,463) = rxt(k,247)*y(k,137)
         mat(k,339) = rxt(k,248)*y(k,137)
         mat(k,280) = rxt(k,140)*y(k,38) + rxt(k,96)*y(k,97) + rxt(k,105)*y(k,137)
         mat(k,507) = rxt(k,262)*y(k,137)
         mat(k,439) = .370_r8*rxt(k,304)*y(k,98)
         mat(k,382) = mat(k,382) + .794_r8*rxt(k,297)*y(k,88) + .794_r8*rxt(k,298) &
                      *y(k,90)
         mat(k,567) = .140_r8*rxt(k,279)*y(k,98)
         mat(k,145) = .200_r8*rxt(k,281)*y(k,137)
         mat(k,288) = .500_r8*rxt(k,288)*y(k,137)
         mat(k,722) = .280_r8*rxt(k,289)*y(k,98)
         mat(k,1278) = mat(k,1278) + .794_r8*rxt(k,297)*y(k,72) + rxt(k,230)*y(k,126) &
                      + rxt(k,257)*y(k,127) + rxt(k,217)*y(k,129) + .250_r8*rxt(k,244) &
                      *y(k,131) + .920_r8*rxt(k,302)*y(k,133) + .470_r8*rxt(k,277) &
                      *y(k,134) + rxt(k,265)*y(k,138) + rxt(k,312)*y(k,140)
         mat(k,1211) = mat(k,1211) + rxt(k,211)*y(k,25) + .794_r8*rxt(k,298)*y(k,72) &
                      + rxt(k,307)*y(k,102) + rxt(k,303)*y(k,133) + .470_r8*rxt(k,276) &
                      *y(k,134) + rxt(k,124)*y(k,137) + rxt(k,313)*y(k,140)
         mat(k,887) = mat(k,887) + rxt(k,212)*y(k,25) + rxt(k,96)*y(k,57)
         mat(k,1009) = mat(k,1009) + .130_r8*rxt(k,226)*y(k,9) + .280_r8*rxt(k,254) &
                      *y(k,13) + .370_r8*rxt(k,304)*y(k,71) + .140_r8*rxt(k,279) &
                      *y(k,74) + .280_r8*rxt(k,289)*y(k,77) + rxt(k,108)*y(k,137)
         mat(k,471) = rxt(k,307)*y(k,90) + rxt(k,308)*y(k,137)
         mat(k,425) = rxt(k,323)*y(k,137)
         mat(k,373) = mat(k,373) + rxt(k,230)*y(k,88) + 2.400_r8*rxt(k,227)*y(k,126) &
                      + rxt(k,228)*y(k,129)
         mat(k,395) = mat(k,395) + rxt(k,257)*y(k,88) + rxt(k,255)*y(k,129)
         mat(k,700) = mat(k,700) + .900_r8*rxt(k,238)*y(k,129) + rxt(k,299)*y(k,133) &
                      + .470_r8*rxt(k,273)*y(k,134) + rxt(k,309)*y(k,140)
         mat(k,759) = mat(k,759) + rxt(k,144)*y(k,41) + rxt(k,217)*y(k,88) &
                      + rxt(k,228)*y(k,126) + rxt(k,255)*y(k,127) + .900_r8*rxt(k,238) &
                      *y(k,128) + 4.000_r8*rxt(k,214)*y(k,129) + rxt(k,300)*y(k,133) &
                      + .730_r8*rxt(k,274)*y(k,134) + rxt(k,283)*y(k,135) &
                      + .300_r8*rxt(k,267)*y(k,139) + .800_r8*rxt(k,310)*y(k,140)
         mat(k,331) = mat(k,331) + .250_r8*rxt(k,244)*y(k,88)
         mat(k,634) = mat(k,634) + .920_r8*rxt(k,302)*y(k,88) + rxt(k,303)*y(k,90) &
                      + rxt(k,299)*y(k,128) + rxt(k,300)*y(k,129)
         mat(k,655) = mat(k,655) + .470_r8*rxt(k,277)*y(k,88) + .470_r8*rxt(k,276) &
                      *y(k,90) + .470_r8*rxt(k,273)*y(k,128) + .730_r8*rxt(k,274) &
                      *y(k,129)
         mat(k,673) = mat(k,673) + rxt(k,283)*y(k,129)
         mat(k,1082) = rxt(k,222)*y(k,36)
         mat(k,1167) = mat(k,1167) + rxt(k,185)*y(k,6) + rxt(k,231)*y(k,10) &
                      + rxt(k,197)*y(k,26) + rxt(k,200)*y(k,29) + rxt(k,218)*y(k,34) &
                      + rxt(k,155)*y(k,41) + rxt(k,224)*y(k,44) + .500_r8*rxt(k,333) &
                      *y(k,46) + rxt(k,247)*y(k,52) + rxt(k,248)*y(k,53) + rxt(k,105) &
                      *y(k,57) + rxt(k,262)*y(k,69) + .200_r8*rxt(k,281)*y(k,75) &
                      + .500_r8*rxt(k,288)*y(k,76) + rxt(k,124)*y(k,90) + rxt(k,108) &
                      *y(k,98) + rxt(k,308)*y(k,102) + rxt(k,323)*y(k,111)
         mat(k,346) = mat(k,346) + rxt(k,265)*y(k,88)
         mat(k,498) = mat(k,498) + .300_r8*rxt(k,267)*y(k,129)
         mat(k,584) = mat(k,584) + rxt(k,312)*y(k,88) + rxt(k,313)*y(k,90) &
                      + rxt(k,309)*y(k,128) + .800_r8*rxt(k,310)*y(k,129)
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
         mat(k,626) = -(rxt(k,299)*y(k,128) + rxt(k,300)*y(k,129) + rxt(k,301) &
                      *y(k,132) + rxt(k,302)*y(k,88) + rxt(k,303)*y(k,90))
         mat(k,692) = -rxt(k,299)*y(k,133)
         mat(k,748) = -rxt(k,300)*y(k,133)
         mat(k,931) = -rxt(k,301)*y(k,133)
         mat(k,1266) = -rxt(k,302)*y(k,133)
         mat(k,1199) = -rxt(k,303)*y(k,133)
         mat(k,434) = rxt(k,305)*y(k,137)
         mat(k,243) = .200_r8*rxt(k,306)*y(k,137)
         mat(k,1199) = mat(k,1199) + 1.700_r8*rxt(k,315)*y(k,122)
         mat(k,311) = 1.700_r8*rxt(k,315)*y(k,90) + 1.640_r8*rxt(k,317)*y(k,137)
         mat(k,1155) = rxt(k,305)*y(k,71) + .200_r8*rxt(k,306)*y(k,73) &
                      + 1.640_r8*rxt(k,317)*y(k,122)
         mat(k,649) = -(rxt(k,273)*y(k,128) + rxt(k,274)*y(k,129) + rxt(k,275) &
                      *y(k,132) + rxt(k,276)*y(k,90) + (rxt(k,277) + rxt(k,278) &
                      ) * y(k,88))
         mat(k,693) = -rxt(k,273)*y(k,134)
         mat(k,749) = -rxt(k,274)*y(k,134)
         mat(k,932) = -rxt(k,275)*y(k,134)
         mat(k,1200) = -rxt(k,276)*y(k,134)
         mat(k,1267) = -(rxt(k,277) + rxt(k,278)) * y(k,134)
         mat(k,562) = .500_r8*rxt(k,280)*y(k,137)
         mat(k,143) = .200_r8*rxt(k,281)*y(k,137)
         mat(k,712) = rxt(k,290)*y(k,137)
         mat(k,1156) = .500_r8*rxt(k,280)*y(k,74) + .200_r8*rxt(k,281)*y(k,75) &
                      + rxt(k,290)*y(k,77)
         mat(k,668) = -(rxt(k,282)*y(k,128) + rxt(k,283)*y(k,129) + rxt(k,284) &
                      *y(k,132) + 4._r8*rxt(k,285)*y(k,135) + rxt(k,286)*y(k,88) &
                      + rxt(k,287)*y(k,90) + rxt(k,291)*y(k,89))
         mat(k,694) = -rxt(k,282)*y(k,135)
         mat(k,750) = -rxt(k,283)*y(k,135)
         mat(k,933) = -rxt(k,284)*y(k,135)
         mat(k,1268) = -rxt(k,286)*y(k,135)
         mat(k,1201) = -rxt(k,287)*y(k,135)
         mat(k,1035) = -rxt(k,291)*y(k,135)
         mat(k,563) = .500_r8*rxt(k,280)*y(k,137)
         mat(k,144) = .500_r8*rxt(k,281)*y(k,137)
         mat(k,1157) = .500_r8*rxt(k,280)*y(k,74) + .500_r8*rxt(k,281)*y(k,75)
         mat(k,1086) = -(rxt(k,86)*y(k,55) + rxt(k,87)*y(k,141) + (rxt(k,90) + rxt(k,91) &
                      ) * y(k,98) + (rxt(k,129) + rxt(k,130)) * y(k,79) + rxt(k,162) &
                      *y(k,16) + rxt(k,163)*y(k,17) + rxt(k,164)*y(k,19) + rxt(k,165) &
                      *y(k,20) + rxt(k,166)*y(k,21) + rxt(k,167)*y(k,22) + rxt(k,168) &
                      *y(k,23) + (rxt(k,169) + rxt(k,170)) * y(k,63) + rxt(k,189) &
                      *y(k,18) + rxt(k,190)*y(k,37) + rxt(k,191)*y(k,56) + (rxt(k,192) &
                      + rxt(k,193)) * y(k,59) + rxt(k,206)*y(k,24) + rxt(k,207) &
                      *y(k,26) + rxt(k,208)*y(k,60) + rxt(k,209)*y(k,61) + rxt(k,210) &
                      *y(k,62) + (rxt(k,221) + rxt(k,222) + rxt(k,223)) * y(k,36))
         mat(k,519) = -rxt(k,86)*y(k,136)
         mat(k,1302) = -rxt(k,87)*y(k,136)
         mat(k,1013) = -(rxt(k,90) + rxt(k,91)) * y(k,136)
         mat(k,98) = -(rxt(k,129) + rxt(k,130)) * y(k,136)
         mat(k,50) = -rxt(k,162)*y(k,136)
         mat(k,76) = -rxt(k,163)*y(k,136)
         mat(k,56) = -rxt(k,164)*y(k,136)
         mat(k,59) = -rxt(k,165)*y(k,136)
         mat(k,62) = -rxt(k,166)*y(k,136)
         mat(k,65) = -rxt(k,167)*y(k,136)
         mat(k,68) = -rxt(k,168)*y(k,136)
         mat(k,827) = -(rxt(k,169) + rxt(k,170)) * y(k,136)
         mat(k,53) = -rxt(k,189)*y(k,136)
         mat(k,186) = -rxt(k,190)*y(k,136)
         mat(k,44) = -rxt(k,191)*y(k,136)
         mat(k,413) = -(rxt(k,192) + rxt(k,193)) * y(k,136)
         mat(k,231) = -rxt(k,206)*y(k,136)
         mat(k,298) = -rxt(k,207)*y(k,136)
         mat(k,91) = -rxt(k,208)*y(k,136)
         mat(k,95) = -rxt(k,209)*y(k,136)
         mat(k,105) = -rxt(k,210)*y(k,136)
         mat(k,487) = -(rxt(k,221) + rxt(k,222) + rxt(k,223)) * y(k,136)
         mat(k,1172) = -(rxt(k,104)*y(k,55) + rxt(k,105)*y(k,57) + rxt(k,106)*y(k,132) &
                      + rxt(k,107)*y(k,97) + rxt(k,108)*y(k,98) + (4._r8*rxt(k,109) &
                      + 4._r8*rxt(k,110)) * y(k,137) + rxt(k,112)*y(k,66) + rxt(k,124) &
                      *y(k,90) + rxt(k,125)*y(k,78) + rxt(k,133)*y(k,89) + rxt(k,134) &
                      *y(k,65) + rxt(k,153)*y(k,42) + (rxt(k,155) + rxt(k,156) &
                      ) * y(k,41) + rxt(k,158)*y(k,63) + rxt(k,161)*y(k,68) + rxt(k,185) &
                      *y(k,6) + rxt(k,187)*y(k,59) + rxt(k,195)*y(k,24) + rxt(k,197) &
                      *y(k,26) + rxt(k,198)*y(k,27) + rxt(k,200)*y(k,29) + rxt(k,202) &
                      *y(k,37) + rxt(k,203)*y(k,60) + rxt(k,204)*y(k,61) + rxt(k,205) &
                      *y(k,62) + rxt(k,213)*y(k,25) + rxt(k,218)*y(k,34) + rxt(k,219) &
                      *y(k,35) + rxt(k,220)*y(k,36) + rxt(k,224)*y(k,44) + rxt(k,231) &
                      *y(k,10) + rxt(k,232)*y(k,11) + rxt(k,234)*y(k,12) + rxt(k,236) &
                      *y(k,28) + rxt(k,241)*y(k,32) + rxt(k,242)*y(k,33) + rxt(k,247) &
                      *y(k,52) + rxt(k,248)*y(k,53) + rxt(k,249)*y(k,103) + rxt(k,250) &
                      *y(k,9) + rxt(k,258)*y(k,14) + rxt(k,259)*y(k,15) + rxt(k,261) &
                      *y(k,31) + rxt(k,262)*y(k,69) + rxt(k,263)*y(k,91) + rxt(k,266) &
                      *y(k,106) + rxt(k,270)*y(k,107) + rxt(k,271)*y(k,13) + rxt(k,272) &
                      *y(k,30) + rxt(k,280)*y(k,74) + rxt(k,281)*y(k,75) + rxt(k,288) &
                      *y(k,76) + rxt(k,290)*y(k,77) + rxt(k,293)*y(k,3) + rxt(k,294) &
                      *y(k,70) + rxt(k,305)*y(k,71) + rxt(k,306)*y(k,73) + rxt(k,308) &
                      *y(k,102) + rxt(k,314)*y(k,123) + rxt(k,317)*y(k,122) + (rxt(k,319) &
                      + rxt(k,333)) * y(k,46) + rxt(k,321)*y(k,101) + rxt(k,323) &
                      *y(k,111) + rxt(k,327)*y(k,108) + rxt(k,332)*y(k,110) + rxt(k,335) &
                      *y(k,84))
         mat(k,520) = -rxt(k,104)*y(k,137)
         mat(k,281) = -rxt(k,105)*y(k,137)
         mat(k,948) = -rxt(k,106)*y(k,137)
         mat(k,892) = -rxt(k,107)*y(k,137)
         mat(k,1014) = -rxt(k,108)*y(k,137)
         mat(k,179) = -rxt(k,112)*y(k,137)
         mat(k,1216) = -rxt(k,124)*y(k,137)
         mat(k,218) = -rxt(k,125)*y(k,137)
         mat(k,1050) = -rxt(k,133)*y(k,137)
         mat(k,418) = -rxt(k,134)*y(k,137)
         mat(k,453) = -rxt(k,153)*y(k,137)
         mat(k,1241) = -(rxt(k,155) + rxt(k,156)) * y(k,137)
         mat(k,828) = -rxt(k,158)*y(k,137)
         mat(k,405) = -rxt(k,161)*y(k,137)
         mat(k,970) = -rxt(k,185)*y(k,137)
         mat(k,414) = -rxt(k,187)*y(k,137)
         mat(k,232) = -rxt(k,195)*y(k,137)
         mat(k,299) = -rxt(k,197)*y(k,137)
         mat(k,79) = -rxt(k,198)*y(k,137)
         mat(k,174) = -rxt(k,200)*y(k,137)
         mat(k,187) = -rxt(k,202)*y(k,137)
         mat(k,92) = -rxt(k,203)*y(k,137)
         mat(k,96) = -rxt(k,204)*y(k,137)
         mat(k,106) = -rxt(k,205)*y(k,137)
         mat(k,793) = -rxt(k,213)*y(k,137)
         mat(k,251) = -rxt(k,218)*y(k,137)
         mat(k,204) = -rxt(k,219)*y(k,137)
         mat(k,488) = -rxt(k,220)*y(k,137)
         mat(k,513) = -rxt(k,224)*y(k,137)
         mat(k,114) = -rxt(k,231)*y(k,137)
         mat(k,151) = -rxt(k,232)*y(k,137)
         mat(k,134) = -rxt(k,234)*y(k,137)
         mat(k,556) = -rxt(k,236)*y(k,137)
         mat(k,321) = -rxt(k,241)*y(k,137)
         mat(k,306) = -rxt(k,242)*y(k,137)
         mat(k,464) = -rxt(k,247)*y(k,137)
         mat(k,340) = -rxt(k,248)*y(k,137)
         mat(k,258) = -rxt(k,249)*y(k,137)
         mat(k,226) = -rxt(k,250)*y(k,137)
         mat(k,192) = -rxt(k,258)*y(k,137)
         mat(k,71) = -rxt(k,259)*y(k,137)
         mat(k,596) = -rxt(k,261)*y(k,137)
         mat(k,508) = -rxt(k,262)*y(k,137)
         mat(k,199) = -rxt(k,263)*y(k,137)
         mat(k,266) = -rxt(k,266)*y(k,137)
         mat(k,210) = -rxt(k,270)*y(k,137)
         mat(k,544) = -rxt(k,271)*y(k,137)
         mat(k,364) = -rxt(k,272)*y(k,137)
         mat(k,569) = -rxt(k,280)*y(k,137)
         mat(k,146) = -rxt(k,281)*y(k,137)
         mat(k,290) = -rxt(k,288)*y(k,137)
         mat(k,726) = -rxt(k,290)*y(k,137)
         mat(k,41) = -rxt(k,293)*y(k,137)
         mat(k,159) = -rxt(k,294)*y(k,137)
         mat(k,442) = -rxt(k,305)*y(k,137)
         mat(k,247) = -rxt(k,306)*y(k,137)
         mat(k,473) = -rxt(k,308)*y(k,137)
         mat(k,102) = -rxt(k,314)*y(k,137)
         mat(k,316) = -rxt(k,317)*y(k,137)
         mat(k,140) = -(rxt(k,319) + rxt(k,333)) * y(k,137)
         mat(k,167) = -rxt(k,321)*y(k,137)
         mat(k,426) = -rxt(k,323)*y(k,137)
         mat(k,239) = -rxt(k,327)*y(k,137)
         mat(k,611) = -rxt(k,332)*y(k,137)
         mat(k,46) = -rxt(k,335)*y(k,137)
         mat(k,226) = mat(k,226) + .130_r8*rxt(k,226)*y(k,98)
         mat(k,151) = mat(k,151) + .500_r8*rxt(k,232)*y(k,137)
         mat(k,544) = mat(k,544) + .360_r8*rxt(k,254)*y(k,98)
         mat(k,793) = mat(k,793) + rxt(k,212)*y(k,97)
         mat(k,204) = mat(k,204) + .300_r8*rxt(k,219)*y(k,137)
         mat(k,488) = mat(k,488) + rxt(k,221)*y(k,136)
         mat(k,863) = rxt(k,142)*y(k,132)
         mat(k,776) = rxt(k,103)*y(k,98) + 2.000_r8*rxt(k,98)*y(k,132)
         mat(k,520) = mat(k,520) + rxt(k,95)*y(k,97) + rxt(k,86)*y(k,136)
         mat(k,281) = mat(k,281) + rxt(k,96)*y(k,97)
         mat(k,414) = mat(k,414) + rxt(k,186)*y(k,97) + rxt(k,192)*y(k,136)
         mat(k,828) = mat(k,828) + rxt(k,157)*y(k,97) + rxt(k,169)*y(k,136)
         mat(k,358) = rxt(k,188)*y(k,97)
         mat(k,405) = mat(k,405) + rxt(k,160)*y(k,97)
         mat(k,442) = mat(k,442) + .320_r8*rxt(k,304)*y(k,98)
         mat(k,384) = .206_r8*rxt(k,296)*y(k,132)
         mat(k,569) = mat(k,569) + .240_r8*rxt(k,279)*y(k,98)
         mat(k,146) = mat(k,146) + .100_r8*rxt(k,281)*y(k,137)
         mat(k,726) = mat(k,726) + .360_r8*rxt(k,289)*y(k,98)
         mat(k,1283) = rxt(k,126)*y(k,132)
         mat(k,1216) = mat(k,1216) + rxt(k,121)*y(k,132)
         mat(k,892) = mat(k,892) + rxt(k,212)*y(k,25) + rxt(k,95)*y(k,55) + rxt(k,96) &
                      *y(k,57) + rxt(k,186)*y(k,59) + rxt(k,157)*y(k,63) + rxt(k,188) &
                      *y(k,67) + rxt(k,160)*y(k,68) + rxt(k,101)*y(k,132)
         mat(k,1014) = mat(k,1014) + .130_r8*rxt(k,226)*y(k,9) + .360_r8*rxt(k,254) &
                      *y(k,13) + rxt(k,103)*y(k,54) + .320_r8*rxt(k,304)*y(k,71) &
                      + .240_r8*rxt(k,279)*y(k,74) + .360_r8*rxt(k,289)*y(k,77) &
                      + 1.156_r8*rxt(k,316)*y(k,122) + rxt(k,102)*y(k,132)
         mat(k,266) = mat(k,266) + .500_r8*rxt(k,266)*y(k,137)
         mat(k,316) = mat(k,316) + 1.156_r8*rxt(k,316)*y(k,98)
         mat(k,102) = mat(k,102) + .500_r8*rxt(k,314)*y(k,137)
         mat(k,703) = .450_r8*rxt(k,239)*y(k,132)
         mat(k,948) = mat(k,948) + rxt(k,142)*y(k,38) + 2.000_r8*rxt(k,98)*y(k,54) &
                      + .206_r8*rxt(k,296)*y(k,72) + rxt(k,126)*y(k,88) + rxt(k,121) &
                      *y(k,90) + rxt(k,101)*y(k,97) + rxt(k,102)*y(k,98) &
                      + .450_r8*rxt(k,239)*y(k,128) + .450_r8*rxt(k,284)*y(k,135) &
                      + .150_r8*rxt(k,268)*y(k,139)
         mat(k,676) = .450_r8*rxt(k,284)*y(k,132)
         mat(k,1087) = rxt(k,221)*y(k,36) + rxt(k,86)*y(k,55) + rxt(k,192)*y(k,59) &
                      + rxt(k,169)*y(k,63) + 2.000_r8*rxt(k,87)*y(k,141)
         mat(k,1172) = mat(k,1172) + .500_r8*rxt(k,232)*y(k,11) + .300_r8*rxt(k,219) &
                      *y(k,35) + .100_r8*rxt(k,281)*y(k,75) + .500_r8*rxt(k,266) &
                      *y(k,106) + .500_r8*rxt(k,314)*y(k,123)
         mat(k,500) = .150_r8*rxt(k,268)*y(k,132)
         mat(k,1303) = 2.000_r8*rxt(k,87)*y(k,136)
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
         mat(k,342) = -(rxt(k,264)*y(k,132) + rxt(k,265)*y(k,88))
         mat(k,913) = -rxt(k,264)*y(k,138)
         mat(k,1252) = -rxt(k,265)*y(k,138)
         mat(k,525) = rxt(k,271)*y(k,137)
         mat(k,261) = .500_r8*rxt(k,266)*y(k,137)
         mat(k,1132) = rxt(k,271)*y(k,13) + .500_r8*rxt(k,266)*y(k,106)
         mat(k,492) = -(rxt(k,267)*y(k,129) + rxt(k,268)*y(k,132) + rxt(k,269)*y(k,88))
         mat(k,741) = -rxt(k,267)*y(k,139)
         mat(k,923) = -rxt(k,268)*y(k,139)
         mat(k,1259) = -rxt(k,269)*y(k,139)
         mat(k,361) = rxt(k,272)*y(k,137)
         mat(k,207) = rxt(k,270)*y(k,137)
         mat(k,1145) = rxt(k,272)*y(k,30) + rxt(k,270)*y(k,107)
         mat(k,578) = -(rxt(k,309)*y(k,128) + rxt(k,310)*y(k,129) + rxt(k,311) &
                      *y(k,132) + rxt(k,312)*y(k,88) + rxt(k,313)*y(k,90))
         mat(k,690) = -rxt(k,309)*y(k,140)
         mat(k,746) = -rxt(k,310)*y(k,140)
         mat(k,929) = -rxt(k,311)*y(k,140)
         mat(k,1264) = -rxt(k,312)*y(k,140)
         mat(k,1196) = -rxt(k,313)*y(k,140)
         mat(k,158) = rxt(k,294)*y(k,137)
         mat(k,242) = .800_r8*rxt(k,306)*y(k,137)
         mat(k,101) = .500_r8*rxt(k,314)*y(k,137)
         mat(k,1152) = rxt(k,294)*y(k,70) + .800_r8*rxt(k,306)*y(k,73) &
                      + .500_r8*rxt(k,314)*y(k,123)
         mat(k,1307) = -(rxt(k,87)*y(k,136) + rxt(k,334)*y(k,112))
         mat(k,1091) = -rxt(k,87)*y(k,141)
         mat(k,119) = -rxt(k,334)*y(k,141)
         mat(k,135) = rxt(k,234)*y(k,137)
         mat(k,193) = rxt(k,258)*y(k,137)
         mat(k,72) = rxt(k,259)*y(k,137)
         mat(k,233) = rxt(k,195)*y(k,137)
         mat(k,795) = rxt(k,213)*y(k,137)
         mat(k,300) = rxt(k,197)*y(k,137)
         mat(k,80) = rxt(k,198)*y(k,137)
         mat(k,558) = rxt(k,236)*y(k,137)
         mat(k,175) = rxt(k,200)*y(k,137)
         mat(k,365) = rxt(k,272)*y(k,137)
         mat(k,598) = rxt(k,261)*y(k,137)
         mat(k,322) = rxt(k,241)*y(k,137)
         mat(k,307) = rxt(k,242)*y(k,137)
         mat(k,205) = rxt(k,219)*y(k,137)
         mat(k,489) = rxt(k,220)*y(k,137)
         mat(k,777) = rxt(k,99)*y(k,132)
         mat(k,521) = rxt(k,104)*y(k,137)
         mat(k,282) = rxt(k,105)*y(k,137)
         mat(k,415) = rxt(k,187)*y(k,137)
         mat(k,107) = rxt(k,205)*y(k,137)
         mat(k,831) = (rxt(k,351)+rxt(k,356))*y(k,67) + (rxt(k,344)+rxt(k,350) &
                       +rxt(k,355))*y(k,68) + rxt(k,158)*y(k,137)
         mat(k,420) = rxt(k,134)*y(k,137)
         mat(k,181) = rxt(k,112)*y(k,137)
         mat(k,359) = (rxt(k,351)+rxt(k,356))*y(k,63)
         mat(k,407) = (rxt(k,344)+rxt(k,350)+rxt(k,355))*y(k,63) + rxt(k,161)*y(k,137)
         mat(k,570) = .500_r8*rxt(k,280)*y(k,137)
         mat(k,47) = rxt(k,335)*y(k,137)
         mat(k,267) = rxt(k,266)*y(k,137)
         mat(k,211) = rxt(k,270)*y(k,137)
         mat(k,952) = rxt(k,99)*y(k,54) + rxt(k,106)*y(k,137)
         mat(k,1176) = rxt(k,234)*y(k,12) + rxt(k,258)*y(k,14) + rxt(k,259)*y(k,15) &
                      + rxt(k,195)*y(k,24) + rxt(k,213)*y(k,25) + rxt(k,197)*y(k,26) &
                      + rxt(k,198)*y(k,27) + rxt(k,236)*y(k,28) + rxt(k,200)*y(k,29) &
                      + rxt(k,272)*y(k,30) + rxt(k,261)*y(k,31) + rxt(k,241)*y(k,32) &
                      + rxt(k,242)*y(k,33) + rxt(k,219)*y(k,35) + rxt(k,220)*y(k,36) &
                      + rxt(k,104)*y(k,55) + rxt(k,105)*y(k,57) + rxt(k,187)*y(k,59) &
                      + rxt(k,205)*y(k,62) + rxt(k,158)*y(k,63) + rxt(k,134)*y(k,65) &
                      + rxt(k,112)*y(k,66) + rxt(k,161)*y(k,68) + .500_r8*rxt(k,280) &
                      *y(k,74) + rxt(k,335)*y(k,84) + rxt(k,266)*y(k,106) + rxt(k,270) &
                      *y(k,107) + rxt(k,106)*y(k,132) + 2.000_r8*rxt(k,109)*y(k,137)
      end do
      end subroutine nlnmat07
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
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 42) = mat(k, 42) + lmat(k, 42)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 51) = mat(k, 51) + lmat(k, 51)
         mat(k, 52) = mat(k, 52) + lmat(k, 52)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 66) = mat(k, 66) + lmat(k, 66)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 73) = mat(k, 73) + lmat(k, 73)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 81) = lmat(k, 81)
         mat(k, 82) = lmat(k, 82)
         mat(k, 83) = lmat(k, 83)
         mat(k, 84) = lmat(k, 84)
         mat(k, 85) = lmat(k, 85)
         mat(k, 86) = lmat(k, 86)
         mat(k, 87) = lmat(k, 87)
         mat(k, 88) = lmat(k, 88)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 108) = lmat(k, 108)
         mat(k, 109) = lmat(k, 109)
         mat(k, 110) = lmat(k, 110)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = lmat(k, 117)
         mat(k, 118) = lmat(k, 118)
         mat(k, 120) = lmat(k, 120)
         mat(k, 121) = lmat(k, 121)
         mat(k, 122) = lmat(k, 122)
         mat(k, 123) = lmat(k, 123)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = lmat(k, 128)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 155) = lmat(k, 155)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 170) = lmat(k, 170)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 189) = lmat(k, 189)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 196) = lmat(k, 196)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = mat(k, 198) + lmat(k, 198)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = lmat(k, 209)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 213) = mat(k, 213) + lmat(k, 213)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 234) = mat(k, 234) + lmat(k, 234)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 237) = lmat(k, 237)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 241) = lmat(k, 241)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 253) = lmat(k, 253)
         mat(k, 254) = lmat(k, 254)
         mat(k, 255) = lmat(k, 255)
         mat(k, 257) = lmat(k, 257)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 263) = lmat(k, 263)
         mat(k, 264) = lmat(k, 264)
         mat(k, 265) = lmat(k, 265)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 273) = mat(k, 273) + lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 286) = lmat(k, 286)
         mat(k, 289) = lmat(k, 289)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 293) = lmat(k, 293)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 301) = mat(k, 301) + lmat(k, 301)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 304) = lmat(k, 304)
         mat(k, 306) = mat(k, 306) + lmat(k, 306)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = lmat(k, 325)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 337) = mat(k, 337) + lmat(k, 337)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 362) = lmat(k, 362)
         mat(k, 363) = lmat(k, 363)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 405) = mat(k, 405) + lmat(k, 405)
         mat(k, 408) = mat(k, 408) + lmat(k, 408)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = lmat(k, 424)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 452) = lmat(k, 452)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 460) = lmat(k, 460)
         mat(k, 461) = mat(k, 461) + lmat(k, 461)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 466) = lmat(k, 466)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 468) = lmat(k, 468)
         mat(k, 470) = lmat(k, 470)
         mat(k, 471) = mat(k, 471) + lmat(k, 471)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 476) = lmat(k, 476)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 478) = lmat(k, 478)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 480) = mat(k, 480) + lmat(k, 480)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 482) = mat(k, 482) + lmat(k, 482)
         mat(k, 485) = lmat(k, 485)
         mat(k, 488) = mat(k, 488) + lmat(k, 488)
         mat(k, 489) = mat(k, 489) + lmat(k, 489)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 503) = mat(k, 503) + lmat(k, 503)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 549) = lmat(k, 549)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 552) = lmat(k, 552)
         mat(k, 554) = lmat(k, 554)
         mat(k, 559) = mat(k, 559) + lmat(k, 559)
         mat(k, 560) = mat(k, 560) + lmat(k, 560)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 565) = mat(k, 565) + lmat(k, 565)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 590) = mat(k, 590) + lmat(k, 590)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 592) = mat(k, 592) + lmat(k, 592)
         mat(k, 594) = lmat(k, 594)
         mat(k, 600) = lmat(k, 600)
         mat(k, 601) = mat(k, 601) + lmat(k, 601)
         mat(k, 602) = mat(k, 602) + lmat(k, 602)
         mat(k, 606) = lmat(k, 606)
         mat(k, 626) = mat(k, 626) + lmat(k, 626)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 668) = mat(k, 668) + lmat(k, 668)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 709) = lmat(k, 709)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 716) = lmat(k, 716)
         mat(k, 753) = mat(k, 753) + lmat(k, 753)
         mat(k, 769) = mat(k, 769) + lmat(k, 769)
         mat(k, 773) = mat(k, 773) + lmat(k, 773)
         mat(k, 780) = mat(k, 780) + lmat(k, 780)
         mat(k, 781) = lmat(k, 781)
         mat(k, 782) = mat(k, 782) + lmat(k, 782)
         mat(k, 783) = mat(k, 783) + lmat(k, 783)
         mat(k, 800) = mat(k, 800) + lmat(k, 800)
         mat(k, 818) = mat(k, 818) + lmat(k, 818)
         mat(k, 820) = mat(k, 820) + lmat(k, 820)
         mat(k, 821) = mat(k, 821) + lmat(k, 821)
         mat(k, 856) = mat(k, 856) + lmat(k, 856)
         mat(k, 886) = mat(k, 886) + lmat(k, 886)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 952) = mat(k, 952) + lmat(k, 952)
         mat(k, 961) = mat(k, 961) + lmat(k, 961)
         mat(k, 964) = mat(k, 964) + lmat(k, 964)
         mat(k, 966) = mat(k, 966) + lmat(k, 966)
         mat(k,1008) = mat(k,1008) + lmat(k,1008)
         mat(k,1011) = mat(k,1011) + lmat(k,1011)
         mat(k,1013) = mat(k,1013) + lmat(k,1013)
         mat(k,1028) = mat(k,1028) + lmat(k,1028)
         mat(k,1044) = mat(k,1044) + lmat(k,1044)
         mat(k,1048) = mat(k,1048) + lmat(k,1048)
         mat(k,1050) = mat(k,1050) + lmat(k,1050)
         mat(k,1053) = mat(k,1053) + lmat(k,1053)
         mat(k,1081) = mat(k,1081) + lmat(k,1081)
         mat(k,1086) = mat(k,1086) + lmat(k,1086)
         mat(k,1172) = mat(k,1172) + lmat(k,1172)
         mat(k,1185) = mat(k,1185) + lmat(k,1185)
         mat(k,1210) = mat(k,1210) + lmat(k,1210)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1217) = mat(k,1217) + lmat(k,1217)
         mat(k,1219) = mat(k,1219) + lmat(k,1219)
         mat(k,1234) = mat(k,1234) + lmat(k,1234)
         mat(k,1235) = mat(k,1235) + lmat(k,1235)
         mat(k,1243) = mat(k,1243) + lmat(k,1243)
         mat(k,1248) = mat(k,1248) + lmat(k,1248)
         mat(k,1277) = mat(k,1277) + lmat(k,1277)
         mat(k,1286) = mat(k,1286) + lmat(k,1286)
         mat(k,1291) = lmat(k,1291)
         mat(k,1293) = lmat(k,1293)
         mat(k,1297) = lmat(k,1297)
         mat(k,1302) = mat(k,1302) + lmat(k,1302)
         mat(k,1303) = mat(k,1303) + lmat(k,1303)
         mat(k,1307) = mat(k,1307) + lmat(k,1307)
         mat(k, 217) = 0._r8
         mat(k, 329) = 0._r8
         mat(k, 333) = 0._r8
         mat(k, 338) = 0._r8
         mat(k, 343) = 0._r8
         mat(k, 348) = 0._r8
         mat(k, 350) = 0._r8
         mat(k, 355) = 0._r8
         mat(k, 375) = 0._r8
         mat(k, 390) = 0._r8
         mat(k, 392) = 0._r8
         mat(k, 397) = 0._r8
         mat(k, 399) = 0._r8
         mat(k, 427) = 0._r8
         mat(k, 430) = 0._r8
         mat(k, 441) = 0._r8
         mat(k, 444) = 0._r8
         mat(k, 456) = 0._r8
         mat(k, 462) = 0._r8
         mat(k, 469) = 0._r8
         mat(k, 475) = 0._r8
         mat(k, 502) = 0._r8
         mat(k, 511) = 0._r8
         mat(k, 527) = 0._r8
         mat(k, 529) = 0._r8
         mat(k, 532) = 0._r8
         mat(k, 533) = 0._r8
         mat(k, 535) = 0._r8
         mat(k, 537) = 0._r8
         mat(k, 538) = 0._r8
         mat(k, 539) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 543) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 547) = 0._r8
         mat(k, 553) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 583) = 0._r8
         mat(k, 586) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 595) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 614) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 627) = 0._r8
         mat(k, 628) = 0._r8
         mat(k, 633) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 647) = 0._r8
         mat(k, 650) = 0._r8
         mat(k, 654) = 0._r8
         mat(k, 657) = 0._r8
         mat(k, 660) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 667) = 0._r8
         mat(k, 672) = 0._r8
         mat(k, 679) = 0._r8
         mat(k, 699) = 0._r8
         mat(k, 704) = 0._r8
         mat(k, 706) = 0._r8
         mat(k, 713) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 719) = 0._r8
         mat(k, 720) = 0._r8
         mat(k, 721) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 725) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 728) = 0._r8
         mat(k, 729) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 760) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 763) = 0._r8
         mat(k, 764) = 0._r8
         mat(k, 767) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 775) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 791) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 798) = 0._r8
         mat(k, 801) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 803) = 0._r8
         mat(k, 807) = 0._r8
         mat(k, 808) = 0._r8
         mat(k, 809) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 819) = 0._r8
         mat(k, 823) = 0._r8
         mat(k, 824) = 0._r8
         mat(k, 825) = 0._r8
         mat(k, 826) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 844) = 0._r8
         mat(k, 848) = 0._r8
         mat(k, 849) = 0._r8
         mat(k, 850) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 859) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 862) = 0._r8
         mat(k, 866) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 875) = 0._r8
         mat(k, 876) = 0._r8
         mat(k, 891) = 0._r8
         mat(k, 896) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 921) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 925) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 930) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 960) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 969) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 974) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 987) = 0._r8
         mat(k, 988) = 0._r8
         mat(k, 990) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 997) = 0._r8
         mat(k, 998) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1032) = 0._r8
         mat(k,1033) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1041) = 0._r8
         mat(k,1042) = 0._r8
         mat(k,1043) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1054) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1171) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1212) = 0._r8
         mat(k,1213) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1220) = 0._r8
         mat(k,1230) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1240) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1306) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 2) = mat(k, 2) - dti(k)
         mat(k, 3) = mat(k, 3) - dti(k)
         mat(k, 4) = mat(k, 4) - dti(k)
         mat(k, 5) = mat(k, 5) - dti(k)
         mat(k, 6) = mat(k, 6) - dti(k)
         mat(k, 7) = mat(k, 7) - dti(k)
         mat(k, 8) = mat(k, 8) - dti(k)
         mat(k, 9) = mat(k, 9) - dti(k)
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
         mat(k, 34) = mat(k, 34) - dti(k)
         mat(k, 35) = mat(k, 35) - dti(k)
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 42) = mat(k, 42) - dti(k)
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 48) = mat(k, 48) - dti(k)
         mat(k, 51) = mat(k, 51) - dti(k)
         mat(k, 54) = mat(k, 54) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 60) = mat(k, 60) - dti(k)
         mat(k, 63) = mat(k, 63) - dti(k)
         mat(k, 66) = mat(k, 66) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 81) = mat(k, 81) - dti(k)
         mat(k, 84) = mat(k, 84) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 89) = mat(k, 89) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 120) = mat(k, 120) - dti(k)
         mat(k, 124) = mat(k, 124) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 176) = mat(k, 176) - dti(k)
         mat(k, 182) = mat(k, 182) - dti(k)
         mat(k, 188) = mat(k, 188) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 206) = mat(k, 206) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 220) = mat(k, 220) - dti(k)
         mat(k, 227) = mat(k, 227) - dti(k)
         mat(k, 234) = mat(k, 234) - dti(k)
         mat(k, 240) = mat(k, 240) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 260) = mat(k, 260) - dti(k)
         mat(k, 268) = mat(k, 268) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 283) = mat(k, 283) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 301) = mat(k, 301) - dti(k)
         mat(k, 308) = mat(k, 308) - dti(k)
         mat(k, 318) = mat(k, 318) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 336) = mat(k, 336) - dti(k)
         mat(k, 342) = mat(k, 342) - dti(k)
         mat(k, 352) = mat(k, 352) - dti(k)
         mat(k, 360) = mat(k, 360) - dti(k)
         mat(k, 369) = mat(k, 369) - dti(k)
         mat(k, 377) = mat(k, 377) - dti(k)
         mat(k, 389) = mat(k, 389) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 408) = mat(k, 408) - dti(k)
         mat(k, 416) = mat(k, 416) - dti(k)
         mat(k, 422) = mat(k, 422) - dti(k)
         mat(k, 429) = mat(k, 429) - dti(k)
         mat(k, 448) = mat(k, 448) - dti(k)
         mat(k, 459) = mat(k, 459) - dti(k)
         mat(k, 467) = mat(k, 467) - dti(k)
         mat(k, 477) = mat(k, 477) - dti(k)
         mat(k, 492) = mat(k, 492) - dti(k)
         mat(k, 503) = mat(k, 503) - dti(k)
         mat(k, 510) = mat(k, 510) - dti(k)
         mat(k, 514) = mat(k, 514) - dti(k)
         mat(k, 530) = mat(k, 530) - dti(k)
         mat(k, 550) = mat(k, 550) - dti(k)
         mat(k, 560) = mat(k, 560) - dti(k)
         mat(k, 578) = mat(k, 578) - dti(k)
         mat(k, 591) = mat(k, 591) - dti(k)
         mat(k, 602) = mat(k, 602) - dti(k)
         mat(k, 626) = mat(k, 626) - dti(k)
         mat(k, 649) = mat(k, 649) - dti(k)
         mat(k, 668) = mat(k, 668) - dti(k)
         mat(k, 695) = mat(k, 695) - dti(k)
         mat(k, 715) = mat(k, 715) - dti(k)
         mat(k, 753) = mat(k, 753) - dti(k)
         mat(k, 769) = mat(k, 769) - dti(k)
         mat(k, 783) = mat(k, 783) - dti(k)
         mat(k, 800) = mat(k, 800) - dti(k)
         mat(k, 820) = mat(k, 820) - dti(k)
         mat(k, 856) = mat(k, 856) - dti(k)
         mat(k, 886) = mat(k, 886) - dti(k)
         mat(k, 943) = mat(k, 943) - dti(k)
         mat(k, 966) = mat(k, 966) - dti(k)
         mat(k,1011) = mat(k,1011) - dti(k)
         mat(k,1048) = mat(k,1048) - dti(k)
         mat(k,1086) = mat(k,1086) - dti(k)
         mat(k,1172) = mat(k,1172) - dti(k)
         mat(k,1217) = mat(k,1217) - dti(k)
         mat(k,1243) = mat(k,1243) - dti(k)
         mat(k,1286) = mat(k,1286) - dti(k)
         mat(k,1307) = mat(k,1307) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
