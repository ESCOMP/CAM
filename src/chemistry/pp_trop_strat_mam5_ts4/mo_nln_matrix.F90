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
         mat(k,37) = -(rxt(k,302)*y(k,139))
         mat(k,1192) = -rxt(k,302)*y(k,3)
         mat(k,1099) = -(rxt(k,182)*y(k,25) + rxt(k,183)*y(k,134) + rxt(k,184) &
                      *y(k,102))
         mat(k,815) = -rxt(k,182)*y(k,4)
         mat(k,1036) = -rxt(k,183)*y(k,4)
         mat(k,1082) = -rxt(k,184)*y(k,4)
         mat(k,1325) = 4.000_r8*rxt(k,185)*y(k,6) + (rxt(k,186)+rxt(k,187))*y(k,41) &
                      + rxt(k,190)*y(k,89) + rxt(k,193)*y(k,98) + rxt(k,334)*y(k,114) &
                      + rxt(k,194)*y(k,139)
         mat(k,73) = rxt(k,172)*y(k,138)
         mat(k,50) = rxt(k,198)*y(k,138)
         mat(k,228) = 2.000_r8*rxt(k,203)*y(k,38) + 2.000_r8*rxt(k,215)*y(k,138) &
                      + 2.000_r8*rxt(k,204)*y(k,139)
         mat(k,296) = rxt(k,205)*y(k,38) + rxt(k,216)*y(k,138) + rxt(k,206)*y(k,139)
         mat(k,191) = 3.000_r8*rxt(k,210)*y(k,38) + 3.000_r8*rxt(k,199)*y(k,138) &
                      + 3.000_r8*rxt(k,211)*y(k,139)
         mat(k,979) = 2.000_r8*rxt(k,203)*y(k,24) + rxt(k,205)*y(k,26) &
                      + 3.000_r8*rxt(k,210)*y(k,37)
         mat(k,884) = (rxt(k,186)+rxt(k,187))*y(k,6)
         mat(k,41) = 2.000_r8*rxt(k,200)*y(k,138)
         mat(k,404) = rxt(k,195)*y(k,98) + rxt(k,201)*y(k,138) + rxt(k,196)*y(k,139)
         mat(k,858) = rxt(k,190)*y(k,6)
         mat(k,1302) = rxt(k,193)*y(k,6) + rxt(k,195)*y(k,59)
         mat(k,724) = rxt(k,334)*y(k,6)
         mat(k,1184) = rxt(k,172)*y(k,17) + rxt(k,198)*y(k,18) + 2.000_r8*rxt(k,215) &
                      *y(k,24) + rxt(k,216)*y(k,26) + 3.000_r8*rxt(k,199)*y(k,37) &
                      + 2.000_r8*rxt(k,200)*y(k,56) + rxt(k,201)*y(k,59)
         mat(k,1270) = rxt(k,194)*y(k,6) + 2.000_r8*rxt(k,204)*y(k,24) + rxt(k,206) &
                      *y(k,26) + 3.000_r8*rxt(k,211)*y(k,37) + rxt(k,196)*y(k,59)
         mat(k,1310) = rxt(k,188)*y(k,41)
         mat(k,868) = rxt(k,188)*y(k,6)
         mat(k,930) = (rxt(k,360)+rxt(k,365))*y(k,66)
         mat(k,336) = (rxt(k,360)+rxt(k,365))*y(k,63)
         mat(k,1330) = -(4._r8*rxt(k,185)*y(k,6) + (rxt(k,186) + rxt(k,187) + rxt(k,188) &
                      ) * y(k,41) + rxt(k,189)*y(k,134) + rxt(k,190)*y(k,89) + rxt(k,191) &
                      *y(k,90) + rxt(k,193)*y(k,98) + rxt(k,194)*y(k,139) + rxt(k,334) &
                      *y(k,114))
         mat(k,889) = -(rxt(k,186) + rxt(k,187) + rxt(k,188)) * y(k,6)
         mat(k,1041) = -rxt(k,189)*y(k,6)
         mat(k,863) = -rxt(k,190)*y(k,6)
         mat(k,926) = -rxt(k,191)*y(k,6)
         mat(k,1307) = -rxt(k,193)*y(k,6)
         mat(k,1275) = -rxt(k,194)*y(k,6)
         mat(k,727) = -rxt(k,334)*y(k,6)
         mat(k,1104) = rxt(k,184)*y(k,102)
         mat(k,239) = rxt(k,192)*y(k,98)
         mat(k,408) = rxt(k,202)*y(k,138)
         mat(k,343) = rxt(k,197)*y(k,98)
         mat(k,1307) = mat(k,1307) + rxt(k,192)*y(k,7) + rxt(k,197)*y(k,66)
         mat(k,1087) = rxt(k,184)*y(k,4)
         mat(k,1189) = rxt(k,202)*y(k,59)
         mat(k,232) = -(rxt(k,192)*y(k,98))
         mat(k,1281) = -rxt(k,192)*y(k,7)
         mat(k,1312) = rxt(k,191)*y(k,90)
         mat(k,895) = rxt(k,191)*y(k,6)
         mat(k,211) = -(rxt(k,234)*y(k,38) + rxt(k,235)*y(k,102) + rxt(k,259)*y(k,139))
         mat(k,954) = -rxt(k,234)*y(k,9)
         mat(k,1046) = -rxt(k,235)*y(k,9)
         mat(k,1214) = -rxt(k,259)*y(k,9)
         mat(k,115) = -(rxt(k,240)*y(k,139))
         mat(k,1200) = -rxt(k,240)*y(k,10)
         mat(k,382) = .800_r8*rxt(k,236)*y(k,128) + .200_r8*rxt(k,237)*y(k,131)
         mat(k,754) = .200_r8*rxt(k,237)*y(k,128)
         mat(k,147) = -(rxt(k,241)*y(k,139))
         mat(k,1205) = -rxt(k,241)*y(k,11)
         mat(k,383) = rxt(k,238)*y(k,134)
         mat(k,990) = rxt(k,238)*y(k,128)
         mat(k,130) = -(rxt(k,242)*y(k,38) + rxt(k,243)*y(k,139))
         mat(k,951) = -rxt(k,242)*y(k,12)
         mat(k,1202) = -rxt(k,243)*y(k,12)
         mat(k,541) = -(rxt(k,262)*y(k,91) + rxt(k,263)*y(k,102) + rxt(k,280)*y(k,139))
         mat(k,1123) = -rxt(k,262)*y(k,13)
         mat(k,1060) = -rxt(k,263)*y(k,13)
         mat(k,1248) = -rxt(k,280)*y(k,13)
         mat(k,432) = .130_r8*rxt(k,313)*y(k,102)
         mat(k,1060) = mat(k,1060) + .130_r8*rxt(k,313)*y(k,70)
         mat(k,176) = -(rxt(k,267)*y(k,139))
         mat(k,1209) = -rxt(k,267)*y(k,14)
         mat(k,415) = rxt(k,265)*y(k,134)
         mat(k,991) = rxt(k,265)*y(k,129)
         mat(k,67) = -(rxt(k,268)*y(k,139))
         mat(k,1194) = -rxt(k,268)*y(k,15)
         mat(k,46) = -(rxt(k,171)*y(k,138))
         mat(k,1154) = -rxt(k,171)*y(k,16)
         mat(k,71) = -(rxt(k,172)*y(k,138))
         mat(k,1161) = -rxt(k,172)*y(k,17)
         mat(k,49) = -(rxt(k,198)*y(k,138))
         mat(k,1155) = -rxt(k,198)*y(k,18)
         mat(k,52) = -(rxt(k,173)*y(k,138))
         mat(k,1156) = -rxt(k,173)*y(k,19)
         mat(k,55) = -(rxt(k,174)*y(k,138))
         mat(k,1157) = -rxt(k,174)*y(k,20)
         mat(k,58) = -(rxt(k,175)*y(k,138))
         mat(k,1158) = -rxt(k,175)*y(k,21)
         mat(k,61) = -(rxt(k,176)*y(k,138))
         mat(k,1159) = -rxt(k,176)*y(k,22)
         mat(k,64) = -(rxt(k,177)*y(k,138))
         mat(k,1160) = -rxt(k,177)*y(k,23)
         mat(k,225) = -(rxt(k,203)*y(k,38) + rxt(k,204)*y(k,139) + rxt(k,215)*y(k,138))
         mat(k,955) = -rxt(k,203)*y(k,24)
         mat(k,1216) = -rxt(k,204)*y(k,24)
         mat(k,1168) = -rxt(k,215)*y(k,24)
         mat(k,809) = -(rxt(k,146)*y(k,38) + rxt(k,182)*y(k,4) + rxt(k,220)*y(k,91) &
                      + rxt(k,221)*y(k,98) + rxt(k,222)*y(k,139))
         mat(k,971) = -rxt(k,146)*y(k,25)
         mat(k,1093) = -rxt(k,182)*y(k,25)
         mat(k,1137) = -rxt(k,220)*y(k,25)
         mat(k,1294) = -rxt(k,221)*y(k,25)
         mat(k,1262) = -rxt(k,222)*y(k,25)
         mat(k,214) = rxt(k,235)*y(k,102)
         mat(k,547) = .500_r8*rxt(k,263)*y(k,102)
         mat(k,288) = .500_r8*rxt(k,251)*y(k,139)
         mat(k,241) = rxt(k,227)*y(k,139)
         mat(k,197) = .300_r8*rxt(k,228)*y(k,139)
         mat(k,509) = (rxt(k,231)+rxt(k,232))*y(k,138)
         mat(k,876) = rxt(k,153)*y(k,131)
         mat(k,564) = .800_r8*rxt(k,256)*y(k,139)
         mat(k,439) = .910_r8*rxt(k,313)*y(k,102)
         mat(k,365) = .072_r8*rxt(k,306)*y(k,89) + .072_r8*rxt(k,307)*y(k,91) &
                      + .206_r8*rxt(k,305)*y(k,134)
         mat(k,575) = .120_r8*rxt(k,288)*y(k,102)
         mat(k,279) = .500_r8*rxt(k,297)*y(k,139)
         mat(k,741) = .600_r8*rxt(k,298)*y(k,102)
         mat(k,850) = .072_r8*rxt(k,306)*y(k,71) + rxt(k,226)*y(k,131) &
                      + .500_r8*rxt(k,253)*y(k,133) + .550_r8*rxt(k,311)*y(k,135) &
                      + .250_r8*rxt(k,286)*y(k,136) + rxt(k,295)*y(k,137) + rxt(k,274) &
                      *y(k,140) + rxt(k,278)*y(k,141) + .250_r8*rxt(k,321)*y(k,142)
         mat(k,1137) = mat(k,1137) + .072_r8*rxt(k,307)*y(k,71) + .600_r8*rxt(k,312) &
                      *y(k,135) + .250_r8*rxt(k,285)*y(k,136) + rxt(k,296)*y(k,137)
         mat(k,1074) = rxt(k,235)*y(k,9) + .500_r8*rxt(k,263)*y(k,13) &
                      + .910_r8*rxt(k,313)*y(k,70) + .120_r8*rxt(k,288)*y(k,73) &
                      + .600_r8*rxt(k,298)*y(k,76)
         mat(k,256) = rxt(k,258)*y(k,139)
         mat(k,388) = .700_r8*rxt(k,237)*y(k,131)
         mat(k,422) = rxt(k,264)*y(k,131)
         mat(k,703) = rxt(k,247)*y(k,131) + .600_r8*rxt(k,308)*y(k,135) &
                      + .250_r8*rxt(k,282)*y(k,136) + rxt(k,291)*y(k,137) &
                      + .250_r8*rxt(k,318)*y(k,142)
         mat(k,779) = rxt(k,153)*y(k,41) + rxt(k,226)*y(k,89) + .700_r8*rxt(k,237) &
                      *y(k,128) + rxt(k,264)*y(k,129) + rxt(k,247)*y(k,130) + ( &
                      + 4.000_r8*rxt(k,223)+2.000_r8*rxt(k,224))*y(k,131) &
                      + 1.200_r8*rxt(k,309)*y(k,135) + .880_r8*rxt(k,283)*y(k,136) &
                      + 2.000_r8*rxt(k,292)*y(k,137) + .800_r8*rxt(k,276)*y(k,141) &
                      + .800_r8*rxt(k,319)*y(k,142)
         mat(k,355) = .500_r8*rxt(k,253)*y(k,89)
         mat(k,1028) = .206_r8*rxt(k,305)*y(k,71) + .450_r8*rxt(k,293)*y(k,137) &
                      + .150_r8*rxt(k,277)*y(k,141)
         mat(k,634) = .550_r8*rxt(k,311)*y(k,89) + .600_r8*rxt(k,312)*y(k,91) &
                      + .600_r8*rxt(k,308)*y(k,130) + 1.200_r8*rxt(k,309)*y(k,131)
         mat(k,656) = .250_r8*rxt(k,286)*y(k,89) + .250_r8*rxt(k,285)*y(k,91) &
                      + .250_r8*rxt(k,282)*y(k,130) + .880_r8*rxt(k,283)*y(k,131)
         mat(k,675) = rxt(k,295)*y(k,89) + rxt(k,296)*y(k,91) + rxt(k,291)*y(k,130) &
                      + 2.000_r8*rxt(k,292)*y(k,131) + .450_r8*rxt(k,293)*y(k,134) &
                      + 4.000_r8*rxt(k,294)*y(k,137)
         mat(k,1176) = (rxt(k,231)+rxt(k,232))*y(k,36)
         mat(k,1262) = mat(k,1262) + .500_r8*rxt(k,251)*y(k,33) + rxt(k,227)*y(k,34) &
                      + .300_r8*rxt(k,228)*y(k,35) + .800_r8*rxt(k,256)*y(k,52) &
                      + .500_r8*rxt(k,297)*y(k,75) + rxt(k,258)*y(k,107)
         mat(k,375) = rxt(k,274)*y(k,89)
         mat(k,489) = rxt(k,278)*y(k,89) + .800_r8*rxt(k,276)*y(k,131) &
                      + .150_r8*rxt(k,277)*y(k,134)
         mat(k,597) = .250_r8*rxt(k,321)*y(k,89) + .250_r8*rxt(k,318)*y(k,130) &
                      + .800_r8*rxt(k,319)*y(k,131)
         mat(k,291) = -(rxt(k,205)*y(k,38) + rxt(k,206)*y(k,139) + rxt(k,216)*y(k,138))
         mat(k,957) = -rxt(k,205)*y(k,26)
         mat(k,1224) = -rxt(k,206)*y(k,26)
         mat(k,1169) = -rxt(k,216)*y(k,26)
         mat(k,75) = -(rxt(k,207)*y(k,139))
         mat(k,1195) = -rxt(k,207)*y(k,27)
         mat(k,525) = -(rxt(k,244)*y(k,91) + rxt(k,245)*y(k,139))
         mat(k,1122) = -rxt(k,244)*y(k,28)
         mat(k,1247) = -rxt(k,245)*y(k,28)
         mat(k,116) = rxt(k,240)*y(k,139)
         mat(k,149) = .500_r8*rxt(k,241)*y(k,139)
         mat(k,540) = .500_r8*rxt(k,263)*y(k,102)
         mat(k,731) = .100_r8*rxt(k,298)*y(k,102)
         mat(k,837) = rxt(k,239)*y(k,128) + .270_r8*rxt(k,266)*y(k,129) + rxt(k,274) &
                      *y(k,140)
         mat(k,1059) = .500_r8*rxt(k,263)*y(k,13) + .100_r8*rxt(k,298)*y(k,76)
         mat(k,386) = rxt(k,239)*y(k,89) + 3.200_r8*rxt(k,236)*y(k,128) &
                      + .800_r8*rxt(k,237)*y(k,131)
         mat(k,419) = .270_r8*rxt(k,266)*y(k,89)
         mat(k,766) = .800_r8*rxt(k,237)*y(k,128)
         mat(k,1247) = mat(k,1247) + rxt(k,240)*y(k,10) + .500_r8*rxt(k,241)*y(k,11)
         mat(k,374) = rxt(k,274)*y(k,89)
         mat(k,168) = -(rxt(k,208)*y(k,38) + rxt(k,209)*y(k,139))
         mat(k,952) = -rxt(k,208)*y(k,29)
         mat(k,1208) = -rxt(k,209)*y(k,29)
         mat(k,345) = -(rxt(k,281)*y(k,139))
         mat(k,1230) = -rxt(k,281)*y(k,30)
         mat(k,827) = .820_r8*rxt(k,266)*y(k,129)
         mat(k,301) = .100_r8*rxt(k,326)*y(k,139)
         mat(k,416) = .820_r8*rxt(k,266)*y(k,89) + .820_r8*rxt(k,264)*y(k,131)
         mat(k,760) = .820_r8*rxt(k,264)*y(k,129)
         mat(k,1230) = mat(k,1230) + .100_r8*rxt(k,326)*y(k,126)
         mat(k,607) = -(rxt(k,269)*y(k,91) + rxt(k,270)*y(k,139))
         mat(k,1128) = -rxt(k,269)*y(k,31)
         mat(k,1253) = -rxt(k,270)*y(k,31)
         mat(k,519) = rxt(k,271)*y(k,139)
         mat(k,571) = .880_r8*rxt(k,288)*y(k,102)
         mat(k,734) = .500_r8*rxt(k,298)*y(k,102)
         mat(k,842) = .020_r8*rxt(k,311)*y(k,135) + .250_r8*rxt(k,286)*y(k,136) &
                      + .250_r8*rxt(k,321)*y(k,142)
         mat(k,1128) = mat(k,1128) + .250_r8*rxt(k,285)*y(k,136) + .250_r8*rxt(k,322) &
                      *y(k,142)
         mat(k,183) = rxt(k,272)*y(k,139)
         mat(k,1065) = .880_r8*rxt(k,288)*y(k,73) + .500_r8*rxt(k,298)*y(k,76)
         mat(k,696) = .250_r8*rxt(k,282)*y(k,136) + .250_r8*rxt(k,318)*y(k,142)
         mat(k,771) = .240_r8*rxt(k,283)*y(k,136) + .500_r8*rxt(k,276)*y(k,141) &
                      + .100_r8*rxt(k,319)*y(k,142)
         mat(k,627) = .020_r8*rxt(k,311)*y(k,89)
         mat(k,651) = .250_r8*rxt(k,286)*y(k,89) + .250_r8*rxt(k,285)*y(k,91) &
                      + .250_r8*rxt(k,282)*y(k,130) + .240_r8*rxt(k,283)*y(k,131)
         mat(k,1253) = mat(k,1253) + rxt(k,271)*y(k,68) + rxt(k,272)*y(k,92)
         mat(k,486) = .500_r8*rxt(k,276)*y(k,131)
         mat(k,594) = .250_r8*rxt(k,321)*y(k,89) + .250_r8*rxt(k,322)*y(k,91) &
                      + .250_r8*rxt(k,318)*y(k,130) + .100_r8*rxt(k,319)*y(k,131)
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
         mat(k,319) = -(rxt(k,250)*y(k,139))
         mat(k,1227) = -rxt(k,250)*y(k,32)
         mat(k,535) = .120_r8*rxt(k,263)*y(k,102)
         mat(k,1048) = .120_r8*rxt(k,263)*y(k,13)
         mat(k,688) = .100_r8*rxt(k,247)*y(k,131) + .150_r8*rxt(k,248)*y(k,134)
         mat(k,758) = .100_r8*rxt(k,247)*y(k,130)
         mat(k,1000) = .150_r8*rxt(k,248)*y(k,130) + .150_r8*rxt(k,293)*y(k,137)
         mat(k,667) = .150_r8*rxt(k,293)*y(k,134)
         mat(k,284) = -(rxt(k,251)*y(k,139))
         mat(k,1223) = -rxt(k,251)*y(k,33)
         mat(k,687) = .400_r8*rxt(k,248)*y(k,134)
         mat(k,999) = .400_r8*rxt(k,248)*y(k,130) + .400_r8*rxt(k,293)*y(k,137)
         mat(k,666) = .400_r8*rxt(k,293)*y(k,134)
         mat(k,240) = -(rxt(k,227)*y(k,139))
         mat(k,1217) = -rxt(k,227)*y(k,34)
         mat(k,384) = .300_r8*rxt(k,237)*y(k,131)
         mat(k,757) = .300_r8*rxt(k,237)*y(k,128) + 2.000_r8*rxt(k,224)*y(k,131) &
                      + .250_r8*rxt(k,309)*y(k,135) + .250_r8*rxt(k,283)*y(k,136) &
                      + .500_r8*rxt(k,276)*y(k,141) + .300_r8*rxt(k,319)*y(k,142)
         mat(k,617) = .250_r8*rxt(k,309)*y(k,131)
         mat(k,645) = .250_r8*rxt(k,283)*y(k,131)
         mat(k,483) = .500_r8*rxt(k,276)*y(k,131)
         mat(k,587) = .300_r8*rxt(k,319)*y(k,131)
         mat(k,194) = -(rxt(k,228)*y(k,139))
         mat(k,1212) = -rxt(k,228)*y(k,35)
         mat(k,756) = rxt(k,225)*y(k,134)
         mat(k,992) = rxt(k,225)*y(k,131)
         mat(k,505) = -(rxt(k,147)*y(k,38) + rxt(k,229)*y(k,139) + (rxt(k,230) &
                      + rxt(k,231) + rxt(k,232)) * y(k,138))
         mat(k,963) = -rxt(k,147)*y(k,36)
         mat(k,1245) = -rxt(k,229)*y(k,36)
         mat(k,1172) = -(rxt(k,230) + rxt(k,231) + rxt(k,232)) * y(k,36)
         mat(k,538) = .100_r8*rxt(k,263)*y(k,102)
         mat(k,1057) = .100_r8*rxt(k,263)*y(k,13)
         mat(k,188) = -(rxt(k,199)*y(k,138) + rxt(k,210)*y(k,38) + rxt(k,211)*y(k,139))
         mat(k,1167) = -rxt(k,199)*y(k,37)
         mat(k,953) = -rxt(k,210)*y(k,37)
         mat(k,1211) = -rxt(k,211)*y(k,37)
         mat(k,976) = -(rxt(k,146)*y(k,25) + rxt(k,147)*y(k,36) + rxt(k,148)*y(k,55) &
                      + rxt(k,149)*y(k,57) + (rxt(k,150) + rxt(k,151)) * y(k,134) &
                      + rxt(k,152)*y(k,102) + rxt(k,159)*y(k,42) + rxt(k,168)*y(k,67) &
                      + rxt(k,203)*y(k,24) + rxt(k,205)*y(k,26) + rxt(k,208)*y(k,29) &
                      + rxt(k,210)*y(k,37) + rxt(k,242)*y(k,12))
         mat(k,812) = -rxt(k,146)*y(k,38)
         mat(k,511) = -rxt(k,147)*y(k,38)
         mat(k,498) = -rxt(k,148)*y(k,38)
         mat(k,270) = -rxt(k,149)*y(k,38)
         mat(k,1033) = -(rxt(k,150) + rxt(k,151)) * y(k,38)
         mat(k,1079) = -rxt(k,152)*y(k,38)
         mat(k,453) = -rxt(k,159)*y(k,38)
         mat(k,398) = -rxt(k,168)*y(k,38)
         mat(k,227) = -rxt(k,203)*y(k,38)
         mat(k,294) = -rxt(k,205)*y(k,38)
         mat(k,172) = -rxt(k,208)*y(k,38)
         mat(k,190) = -rxt(k,210)*y(k,38)
         mat(k,133) = -rxt(k,242)*y(k,38)
         mat(k,1322) = rxt(k,187)*y(k,41)
         mat(k,47) = 4.000_r8*rxt(k,171)*y(k,138)
         mat(k,72) = rxt(k,172)*y(k,138)
         mat(k,53) = 3.000_r8*rxt(k,173)*y(k,138)
         mat(k,56) = 3.000_r8*rxt(k,174)*y(k,138)
         mat(k,59) = 2.000_r8*rxt(k,175)*y(k,138)
         mat(k,62) = rxt(k,176)*y(k,138)
         mat(k,65) = 2.000_r8*rxt(k,177)*y(k,138)
         mat(k,76) = 3.000_r8*rxt(k,207)*y(k,139)
         mat(k,172) = mat(k,172) + rxt(k,209)*y(k,139)
         mat(k,881) = rxt(k,187)*y(k,6) + (4.000_r8*rxt(k,154)+2.000_r8*rxt(k,156)) &
                      *y(k,41) + rxt(k,158)*y(k,89) + rxt(k,163)*y(k,98) + rxt(k,335) &
                      *y(k,114) + rxt(k,153)*y(k,131) + rxt(k,164)*y(k,139)
         mat(k,93) = 2.000_r8*rxt(k,217)*y(k,138) + 2.000_r8*rxt(k,212)*y(k,139)
         mat(k,97) = rxt(k,218)*y(k,138) + rxt(k,213)*y(k,139)
         mat(k,104) = rxt(k,219)*y(k,138) + rxt(k,214)*y(k,139)
         mat(k,939) = rxt(k,166)*y(k,98) + rxt(k,178)*y(k,138) + rxt(k,167)*y(k,139)
         mat(k,855) = rxt(k,158)*y(k,41)
         mat(k,1299) = rxt(k,163)*y(k,41) + rxt(k,166)*y(k,63)
         mat(k,721) = rxt(k,335)*y(k,41)
         mat(k,784) = rxt(k,153)*y(k,41)
         mat(k,1181) = 4.000_r8*rxt(k,171)*y(k,16) + rxt(k,172)*y(k,17) &
                      + 3.000_r8*rxt(k,173)*y(k,19) + 3.000_r8*rxt(k,174)*y(k,20) &
                      + 2.000_r8*rxt(k,175)*y(k,21) + rxt(k,176)*y(k,22) &
                      + 2.000_r8*rxt(k,177)*y(k,23) + 2.000_r8*rxt(k,217)*y(k,60) &
                      + rxt(k,218)*y(k,61) + rxt(k,219)*y(k,62) + rxt(k,178)*y(k,63)
         mat(k,1267) = 3.000_r8*rxt(k,207)*y(k,27) + rxt(k,209)*y(k,29) + rxt(k,164) &
                      *y(k,41) + 2.000_r8*rxt(k,212)*y(k,60) + rxt(k,213)*y(k,61) &
                      + rxt(k,214)*y(k,62) + rxt(k,167)*y(k,63)
         mat(k,950) = rxt(k,159)*y(k,42)
         mat(k,867) = 2.000_r8*rxt(k,155)*y(k,41)
         mat(k,446) = rxt(k,159)*y(k,38) + (rxt(k,358)+rxt(k,363)+rxt(k,368))*y(k,63)
         mat(k,929) = (rxt(k,358)+rxt(k,363)+rxt(k,368))*y(k,42) + (rxt(k,353) &
                       +rxt(k,359)+rxt(k,364))*y(k,67)
         mat(k,394) = (rxt(k,353)+rxt(k,359)+rxt(k,364))*y(k,63)
         mat(k,866) = 2.000_r8*rxt(k,180)*y(k,41)
         mat(k,878) = -(rxt(k,153)*y(k,131) + (4._r8*rxt(k,154) + 4._r8*rxt(k,155) &
                      + 4._r8*rxt(k,156) + 4._r8*rxt(k,180)) * y(k,41) + rxt(k,157) &
                      *y(k,134) + rxt(k,158)*y(k,89) + rxt(k,160)*y(k,90) + rxt(k,163) &
                      *y(k,98) + (rxt(k,164) + rxt(k,165)) * y(k,139) + (rxt(k,186) &
                      + rxt(k,187) + rxt(k,188)) * y(k,6) + rxt(k,335)*y(k,114))
         mat(k,781) = -rxt(k,153)*y(k,41)
         mat(k,1030) = -rxt(k,157)*y(k,41)
         mat(k,852) = -rxt(k,158)*y(k,41)
         mat(k,915) = -rxt(k,160)*y(k,41)
         mat(k,1296) = -rxt(k,163)*y(k,41)
         mat(k,1264) = -(rxt(k,164) + rxt(k,165)) * y(k,41)
         mat(k,1319) = -(rxt(k,186) + rxt(k,187) + rxt(k,188)) * y(k,41)
         mat(k,719) = -rxt(k,335)*y(k,41)
         mat(k,973) = rxt(k,168)*y(k,67) + rxt(k,152)*y(k,102) + rxt(k,151)*y(k,134)
         mat(k,450) = rxt(k,161)*y(k,98)
         mat(k,936) = rxt(k,179)*y(k,138)
         mat(k,396) = rxt(k,168)*y(k,38) + rxt(k,169)*y(k,98) + rxt(k,170)*y(k,139)
         mat(k,1296) = mat(k,1296) + rxt(k,161)*y(k,42) + rxt(k,169)*y(k,67)
         mat(k,1076) = rxt(k,152)*y(k,38)
         mat(k,155) = rxt(k,340)*y(k,114)
         mat(k,719) = mat(k,719) + rxt(k,340)*y(k,104)
         mat(k,1030) = mat(k,1030) + rxt(k,151)*y(k,38)
         mat(k,1178) = rxt(k,179)*y(k,63)
         mat(k,1264) = mat(k,1264) + rxt(k,170)*y(k,67)
         mat(k,449) = -(rxt(k,159)*y(k,38) + rxt(k,161)*y(k,98) + rxt(k,162)*y(k,139) &
                      + (rxt(k,358) + rxt(k,363) + rxt(k,368)) * y(k,63))
         mat(k,961) = -rxt(k,159)*y(k,42)
         mat(k,1288) = -rxt(k,161)*y(k,42)
         mat(k,1239) = -rxt(k,162)*y(k,42)
         mat(k,934) = -(rxt(k,358) + rxt(k,363) + rxt(k,368)) * y(k,42)
         mat(k,871) = rxt(k,160)*y(k,90)
         mat(k,901) = rxt(k,160)*y(k,41)
         mat(k,581) = -(rxt(k,233)*y(k,139))
         mat(k,1251) = -rxt(k,233)*y(k,44)
         mat(k,1091) = rxt(k,182)*y(k,25)
         mat(k,213) = .630_r8*rxt(k,235)*y(k,102)
         mat(k,542) = .560_r8*rxt(k,263)*y(k,102)
         mat(k,807) = rxt(k,182)*y(k,4) + rxt(k,146)*y(k,38) + rxt(k,220)*y(k,91) &
                      + rxt(k,221)*y(k,98) + rxt(k,222)*y(k,139)
         mat(k,169) = rxt(k,208)*y(k,38)
         mat(k,606) = rxt(k,269)*y(k,91) + rxt(k,270)*y(k,139)
         mat(k,966) = rxt(k,146)*y(k,25) + rxt(k,208)*y(k,29)
         mat(k,333) = rxt(k,257)*y(k,139)
         mat(k,434) = .620_r8*rxt(k,313)*y(k,102)
         mat(k,570) = .650_r8*rxt(k,288)*y(k,102)
         mat(k,733) = .560_r8*rxt(k,298)*y(k,102)
         mat(k,840) = .220_r8*rxt(k,286)*y(k,136) + .250_r8*rxt(k,321)*y(k,142)
         mat(k,1126) = rxt(k,220)*y(k,25) + rxt(k,269)*y(k,31) + .220_r8*rxt(k,285) &
                      *y(k,136) + .500_r8*rxt(k,322)*y(k,142)
         mat(k,1291) = rxt(k,221)*y(k,25) + rxt(k,329)*y(k,105)
         mat(k,1063) = .630_r8*rxt(k,235)*y(k,9) + .560_r8*rxt(k,263)*y(k,13) &
                      + .620_r8*rxt(k,313)*y(k,70) + .650_r8*rxt(k,288)*y(k,73) &
                      + .560_r8*rxt(k,298)*y(k,76)
         mat(k,163) = rxt(k,329)*y(k,98) + rxt(k,330)*y(k,139)
         mat(k,694) = .220_r8*rxt(k,282)*y(k,136) + .250_r8*rxt(k,318)*y(k,142)
         mat(k,769) = .110_r8*rxt(k,283)*y(k,136) + .200_r8*rxt(k,319)*y(k,142)
         mat(k,649) = .220_r8*rxt(k,286)*y(k,89) + .220_r8*rxt(k,285)*y(k,91) &
                      + .220_r8*rxt(k,282)*y(k,130) + .110_r8*rxt(k,283)*y(k,131)
         mat(k,1251) = mat(k,1251) + rxt(k,222)*y(k,25) + rxt(k,270)*y(k,31) &
                      + rxt(k,257)*y(k,53) + rxt(k,330)*y(k,105)
         mat(k,592) = .250_r8*rxt(k,321)*y(k,89) + .500_r8*rxt(k,322)*y(k,91) &
                      + .250_r8*rxt(k,318)*y(k,130) + .200_r8*rxt(k,319)*y(k,131)
         mat(k,537) = .200_r8*rxt(k,263)*y(k,102)
         mat(k,320) = rxt(k,250)*y(k,139)
         mat(k,285) = .500_r8*rxt(k,251)*y(k,139)
         mat(k,580) = rxt(k,233)*y(k,139)
         mat(k,561) = .800_r8*rxt(k,256)*y(k,139)
         mat(k,332) = rxt(k,257)*y(k,139)
         mat(k,276) = .500_r8*rxt(k,297)*y(k,139)
         mat(k,730) = .100_r8*rxt(k,298)*y(k,102)
         mat(k,833) = rxt(k,249)*y(k,130)
         mat(k,1053) = .200_r8*rxt(k,263)*y(k,13) + .100_r8*rxt(k,298)*y(k,76)
         mat(k,690) = rxt(k,249)*y(k,89) + 4.000_r8*rxt(k,246)*y(k,130) &
                      + .900_r8*rxt(k,247)*y(k,131) + 2.000_r8*rxt(k,291)*y(k,137) &
                      + rxt(k,318)*y(k,142)
         mat(k,763) = .900_r8*rxt(k,247)*y(k,130) + rxt(k,292)*y(k,137)
         mat(k,1010) = .450_r8*rxt(k,293)*y(k,137)
         mat(k,668) = 2.000_r8*rxt(k,291)*y(k,130) + rxt(k,292)*y(k,131) &
                      + .450_r8*rxt(k,293)*y(k,134) + 4.000_r8*rxt(k,294)*y(k,137)
         mat(k,1240) = rxt(k,250)*y(k,32) + .500_r8*rxt(k,251)*y(k,33) + rxt(k,233) &
                      *y(k,44) + .800_r8*rxt(k,256)*y(k,52) + rxt(k,257)*y(k,53) &
                      + .500_r8*rxt(k,297)*y(k,75)
         mat(k,589) = rxt(k,318)*y(k,130)
         mat(k,136) = -(rxt(k,327)*y(k,91) + (rxt(k,328) + rxt(k,342)) * y(k,139))
         mat(k,1108) = -rxt(k,327)*y(k,46)
         mat(k,1203) = -(rxt(k,328) + rxt(k,342)) * y(k,46)
         mat(k,351) = rxt(k,252)*y(k,134)
         mat(k,987) = rxt(k,252)*y(k,133)
         mat(k,562) = -(rxt(k,256)*y(k,139))
         mat(k,1249) = -rxt(k,256)*y(k,52)
         mat(k,838) = .020_r8*rxt(k,311)*y(k,135) + .530_r8*rxt(k,286)*y(k,136) &
                      + .250_r8*rxt(k,321)*y(k,142)
         mat(k,1124) = .530_r8*rxt(k,285)*y(k,136) + .250_r8*rxt(k,322)*y(k,142)
         mat(k,1339) = rxt(k,255)*y(k,132)
         mat(k,692) = .530_r8*rxt(k,282)*y(k,136) + .250_r8*rxt(k,318)*y(k,142)
         mat(k,767) = .260_r8*rxt(k,283)*y(k,136) + .100_r8*rxt(k,319)*y(k,142)
         mat(k,207) = rxt(k,255)*y(k,99)
         mat(k,623) = .020_r8*rxt(k,311)*y(k,89)
         mat(k,648) = .530_r8*rxt(k,286)*y(k,89) + .530_r8*rxt(k,285)*y(k,91) &
                      + .530_r8*rxt(k,282)*y(k,130) + .260_r8*rxt(k,283)*y(k,131)
         mat(k,591) = .250_r8*rxt(k,321)*y(k,89) + .250_r8*rxt(k,322)*y(k,91) &
                      + .250_r8*rxt(k,318)*y(k,130) + .100_r8*rxt(k,319)*y(k,131)
         mat(k,331) = -(rxt(k,257)*y(k,139))
         mat(k,1229) = -rxt(k,257)*y(k,53)
         mat(k,560) = .200_r8*rxt(k,256)*y(k,139)
         mat(k,826) = .020_r8*rxt(k,311)*y(k,135) + .250_r8*rxt(k,321)*y(k,142)
         mat(k,1112) = .250_r8*rxt(k,322)*y(k,142)
         mat(k,689) = .250_r8*rxt(k,318)*y(k,142)
         mat(k,759) = .100_r8*rxt(k,319)*y(k,142)
         mat(k,619) = .020_r8*rxt(k,311)*y(k,89)
         mat(k,1229) = mat(k,1229) + .200_r8*rxt(k,256)*y(k,52)
         mat(k,588) = .250_r8*rxt(k,321)*y(k,89) + .250_r8*rxt(k,322)*y(k,91) &
                      + .250_r8*rxt(k,318)*y(k,130) + .100_r8*rxt(k,319)*y(k,131)
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
         mat(k,794) = -((rxt(k,106) + rxt(k,107) + rxt(k,108)) * y(k,134) + rxt(k,109) &
                      *y(k,99) + rxt(k,112)*y(k,102))
         mat(k,1027) = -(rxt(k,106) + rxt(k,107) + rxt(k,108)) * y(k,54)
         mat(k,1342) = -rxt(k,109)*y(k,54)
         mat(k,1073) = -rxt(k,112)*y(k,54)
         mat(k,808) = rxt(k,222)*y(k,139)
         mat(k,508) = rxt(k,231)*y(k,138)
         mat(k,970) = rxt(k,148)*y(k,55)
         mat(k,496) = rxt(k,148)*y(k,38) + rxt(k,104)*y(k,98) + rxt(k,85)*y(k,138) &
                      + rxt(k,113)*y(k,139)
         mat(k,403) = rxt(k,202)*y(k,138)
         mat(k,935) = rxt(k,179)*y(k,138)
         mat(k,312) = rxt(k,134)*y(k,139)
         mat(k,1293) = rxt(k,104)*y(k,55) + rxt(k,116)*y(k,139)
         mat(k,165) = rxt(k,330)*y(k,139)
         mat(k,326) = rxt(k,336)*y(k,139)
         mat(k,717) = rxt(k,341)*y(k,139)
         mat(k,1175) = rxt(k,231)*y(k,36) + rxt(k,85)*y(k,55) + rxt(k,202)*y(k,59) &
                      + rxt(k,179)*y(k,63)
         mat(k,1261) = rxt(k,222)*y(k,25) + rxt(k,113)*y(k,55) + rxt(k,134)*y(k,77) &
                      + rxt(k,116)*y(k,98) + rxt(k,330)*y(k,105) + rxt(k,336)*y(k,112) &
                      + rxt(k,341)*y(k,114)
         mat(k,495) = -(rxt(k,85)*y(k,138) + rxt(k,104)*y(k,98) + rxt(k,113)*y(k,139) &
                      + rxt(k,148)*y(k,38))
         mat(k,1171) = -rxt(k,85)*y(k,55)
         mat(k,1290) = -rxt(k,104)*y(k,55)
         mat(k,1244) = -rxt(k,113)*y(k,55)
         mat(k,962) = -rxt(k,148)*y(k,55)
         mat(k,504) = rxt(k,232)*y(k,138)
         mat(k,793) = rxt(k,106)*y(k,134)
         mat(k,1013) = rxt(k,106)*y(k,54)
         mat(k,1171) = mat(k,1171) + rxt(k,232)*y(k,36)
         mat(k,40) = -(rxt(k,200)*y(k,138))
         mat(k,1153) = -rxt(k,200)*y(k,56)
         mat(k,268) = -(rxt(k,105)*y(k,98) + rxt(k,114)*y(k,139) + rxt(k,149)*y(k,38))
         mat(k,1282) = -rxt(k,105)*y(k,57)
         mat(k,1221) = -rxt(k,114)*y(k,57)
         mat(k,956) = -rxt(k,149)*y(k,57)
         mat(k,998) = 2.000_r8*rxt(k,120)*y(k,134)
         mat(k,1221) = mat(k,1221) + 2.000_r8*rxt(k,119)*y(k,139)
         mat(k,119) = rxt(k,343)*y(k,143)
         mat(k,1359) = rxt(k,343)*y(k,116)
         mat(k,402) = -(rxt(k,195)*y(k,98) + rxt(k,196)*y(k,139) + (rxt(k,201) &
                      + rxt(k,202)) * y(k,138))
         mat(k,1286) = -rxt(k,195)*y(k,59)
         mat(k,1235) = -rxt(k,196)*y(k,59)
         mat(k,1170) = -(rxt(k,201) + rxt(k,202)) * y(k,59)
         mat(k,1090) = rxt(k,182)*y(k,25) + rxt(k,183)*y(k,134)
         mat(k,804) = rxt(k,182)*y(k,4)
         mat(k,1008) = rxt(k,183)*y(k,4)
         mat(k,92) = -(rxt(k,212)*y(k,139) + rxt(k,217)*y(k,138))
         mat(k,1196) = -rxt(k,212)*y(k,60)
         mat(k,1163) = -rxt(k,217)*y(k,60)
         mat(k,96) = -(rxt(k,213)*y(k,139) + rxt(k,218)*y(k,138))
         mat(k,1197) = -rxt(k,213)*y(k,61)
         mat(k,1164) = -rxt(k,218)*y(k,61)
         mat(k,103) = -(rxt(k,214)*y(k,139) + rxt(k,219)*y(k,138))
         mat(k,1199) = -rxt(k,214)*y(k,62)
         mat(k,1165) = -rxt(k,219)*y(k,62)
         mat(k,938) = -(rxt(k,166)*y(k,98) + rxt(k,167)*y(k,139) + (rxt(k,178) &
                      + rxt(k,179)) * y(k,138) + (rxt(k,353) + rxt(k,359) + rxt(k,364) &
                      ) * y(k,67) + (rxt(k,358) + rxt(k,363) + rxt(k,368)) * y(k,42) &
                      + (rxt(k,360) + rxt(k,365)) * y(k,66))
         mat(k,1298) = -rxt(k,166)*y(k,63)
         mat(k,1266) = -rxt(k,167)*y(k,63)
         mat(k,1180) = -(rxt(k,178) + rxt(k,179)) * y(k,63)
         mat(k,397) = -(rxt(k,353) + rxt(k,359) + rxt(k,364)) * y(k,63)
         mat(k,452) = -(rxt(k,358) + rxt(k,363) + rxt(k,368)) * y(k,63)
         mat(k,338) = -(rxt(k,360) + rxt(k,365)) * y(k,63)
         mat(k,132) = rxt(k,242)*y(k,38)
         mat(k,226) = rxt(k,203)*y(k,38)
         mat(k,811) = rxt(k,146)*y(k,38)
         mat(k,293) = rxt(k,205)*y(k,38)
         mat(k,171) = 2.000_r8*rxt(k,208)*y(k,38)
         mat(k,510) = rxt(k,147)*y(k,38)
         mat(k,189) = rxt(k,210)*y(k,38)
         mat(k,975) = rxt(k,242)*y(k,12) + rxt(k,203)*y(k,24) + rxt(k,146)*y(k,25) &
                      + rxt(k,205)*y(k,26) + 2.000_r8*rxt(k,208)*y(k,29) + rxt(k,147) &
                      *y(k,36) + rxt(k,210)*y(k,37) + rxt(k,148)*y(k,55) + rxt(k,149) &
                      *y(k,57) + rxt(k,168)*y(k,67) + rxt(k,150)*y(k,134)
         mat(k,880) = rxt(k,165)*y(k,139)
         mat(k,497) = rxt(k,148)*y(k,38)
         mat(k,269) = rxt(k,149)*y(k,38)
         mat(k,397) = mat(k,397) + rxt(k,168)*y(k,38)
         mat(k,1032) = rxt(k,150)*y(k,38)
         mat(k,1266) = mat(k,1266) + rxt(k,165)*y(k,41)
         mat(k,410) = -(rxt(k,143)*y(k,139))
         mat(k,1236) = -rxt(k,143)*y(k,64)
         mat(k,805) = rxt(k,220)*y(k,91)
         mat(k,524) = rxt(k,244)*y(k,91)
         mat(k,605) = rxt(k,269)*y(k,91)
         mat(k,448) = (rxt(k,358)+rxt(k,363)+rxt(k,368))*y(k,63)
         mat(k,137) = rxt(k,327)*y(k,91)
         mat(k,933) = (rxt(k,358)+rxt(k,363)+rxt(k,368))*y(k,42)
         mat(k,900) = rxt(k,142)*y(k,139)
         mat(k,1115) = rxt(k,220)*y(k,25) + rxt(k,244)*y(k,28) + rxt(k,269)*y(k,31) &
                      + rxt(k,327)*y(k,46)
         mat(k,1236) = mat(k,1236) + rxt(k,142)*y(k,90)
         mat(k,218) = -(rxt(k,121)*y(k,139))
         mat(k,1215) = -rxt(k,121)*y(k,65)
         mat(k,894) = rxt(k,140)*y(k,134)
         mat(k,995) = rxt(k,140)*y(k,90)
         mat(k,337) = -(rxt(k,197)*y(k,98) + (rxt(k,360) + rxt(k,365)) * y(k,63))
         mat(k,1284) = -rxt(k,197)*y(k,66)
         mat(k,931) = -(rxt(k,360) + rxt(k,365)) * y(k,66)
         mat(k,1313) = rxt(k,189)*y(k,134)
         mat(k,1001) = rxt(k,189)*y(k,6)
         mat(k,395) = -(rxt(k,168)*y(k,38) + rxt(k,169)*y(k,98) + rxt(k,170)*y(k,139) &
                      + (rxt(k,353) + rxt(k,359) + rxt(k,364)) * y(k,63))
         mat(k,960) = -rxt(k,168)*y(k,67)
         mat(k,1285) = -rxt(k,169)*y(k,67)
         mat(k,1234) = -rxt(k,170)*y(k,67)
         mat(k,932) = -(rxt(k,353) + rxt(k,359) + rxt(k,364)) * y(k,67)
         mat(k,870) = rxt(k,157)*y(k,134)
         mat(k,447) = rxt(k,162)*y(k,139)
         mat(k,1007) = rxt(k,157)*y(k,41)
         mat(k,1234) = mat(k,1234) + rxt(k,162)*y(k,42)
         mat(k,518) = -(rxt(k,271)*y(k,139))
         mat(k,1246) = -rxt(k,271)*y(k,68)
         mat(k,277) = .500_r8*rxt(k,297)*y(k,139)
         mat(k,836) = .020_r8*rxt(k,311)*y(k,135) + .220_r8*rxt(k,286)*y(k,136) &
                      + .250_r8*rxt(k,321)*y(k,142)
         mat(k,1121) = .220_r8*rxt(k,285)*y(k,136) + .250_r8*rxt(k,322)*y(k,142)
         mat(k,262) = .500_r8*rxt(k,275)*y(k,139)
         mat(k,691) = .220_r8*rxt(k,282)*y(k,136) + .250_r8*rxt(k,318)*y(k,142)
         mat(k,765) = .230_r8*rxt(k,283)*y(k,136) + .200_r8*rxt(k,276)*y(k,141) &
                      + .100_r8*rxt(k,319)*y(k,142)
         mat(k,622) = .020_r8*rxt(k,311)*y(k,89)
         mat(k,647) = .220_r8*rxt(k,286)*y(k,89) + .220_r8*rxt(k,285)*y(k,91) &
                      + .220_r8*rxt(k,282)*y(k,130) + .230_r8*rxt(k,283)*y(k,131)
         mat(k,1246) = mat(k,1246) + .500_r8*rxt(k,297)*y(k,75) + .500_r8*rxt(k,275) &
                      *y(k,110)
         mat(k,485) = .200_r8*rxt(k,276)*y(k,131)
         mat(k,590) = .250_r8*rxt(k,321)*y(k,89) + .250_r8*rxt(k,322)*y(k,91) &
                      + .250_r8*rxt(k,318)*y(k,130) + .100_r8*rxt(k,319)*y(k,131)
         mat(k,157) = -(rxt(k,303)*y(k,139))
         mat(k,1206) = -rxt(k,303)*y(k,69)
         mat(k,823) = .330_r8*rxt(k,311)*y(k,135)
         mat(k,1109) = rxt(k,316)*y(k,106) + .400_r8*rxt(k,312)*y(k,135)
         mat(k,471) = rxt(k,316)*y(k,91) + rxt(k,317)*y(k,139)
         mat(k,685) = .400_r8*rxt(k,308)*y(k,135)
         mat(k,755) = .300_r8*rxt(k,309)*y(k,135)
         mat(k,616) = .330_r8*rxt(k,311)*y(k,89) + .400_r8*rxt(k,312)*y(k,91) &
                      + .400_r8*rxt(k,308)*y(k,130) + .300_r8*rxt(k,309)*y(k,131)
         mat(k,1206) = mat(k,1206) + rxt(k,317)*y(k,106)
         mat(k,430) = -(rxt(k,304)*y(k,91) + rxt(k,313)*y(k,102) + rxt(k,314)*y(k,139))
         mat(k,1116) = -rxt(k,304)*y(k,70)
         mat(k,1052) = -rxt(k,313)*y(k,70)
         mat(k,1238) = -rxt(k,314)*y(k,70)
         mat(k,361) = -(rxt(k,305)*y(k,134) + rxt(k,306)*y(k,89) + rxt(k,307)*y(k,91))
         mat(k,1004) = -rxt(k,305)*y(k,71)
         mat(k,829) = -rxt(k,306)*y(k,71)
         mat(k,1114) = -rxt(k,307)*y(k,71)
         mat(k,429) = rxt(k,304)*y(k,91)
         mat(k,1114) = mat(k,1114) + rxt(k,304)*y(k,70)
         mat(k,244) = -(rxt(k,315)*y(k,139))
         mat(k,1218) = -rxt(k,315)*y(k,72)
         mat(k,996) = rxt(k,310)*y(k,135)
         mat(k,618) = rxt(k,310)*y(k,134)
         mat(k,569) = -(rxt(k,288)*y(k,102) + rxt(k,289)*y(k,139))
         mat(k,1062) = -rxt(k,288)*y(k,73)
         mat(k,1250) = -rxt(k,289)*y(k,73)
         mat(k,433) = .300_r8*rxt(k,313)*y(k,102)
         mat(k,363) = .167_r8*rxt(k,306)*y(k,89) + .167_r8*rxt(k,307)*y(k,91) &
                      + .167_r8*rxt(k,305)*y(k,134)
         mat(k,839) = .167_r8*rxt(k,306)*y(k,71) + .230_r8*rxt(k,311)*y(k,135)
         mat(k,1125) = .167_r8*rxt(k,307)*y(k,71) + .250_r8*rxt(k,312)*y(k,135)
         mat(k,1062) = mat(k,1062) + .300_r8*rxt(k,313)*y(k,70) + 1.122_r8*rxt(k,325) &
                      *y(k,126)
         mat(k,302) = 1.122_r8*rxt(k,325)*y(k,102)
         mat(k,693) = .250_r8*rxt(k,308)*y(k,135)
         mat(k,768) = .190_r8*rxt(k,309)*y(k,135)
         mat(k,1017) = .167_r8*rxt(k,305)*y(k,71)
         mat(k,624) = .230_r8*rxt(k,311)*y(k,89) + .250_r8*rxt(k,312)*y(k,91) &
                      + .250_r8*rxt(k,308)*y(k,130) + .190_r8*rxt(k,309)*y(k,131)
         mat(k,142) = -(rxt(k,290)*y(k,139))
         mat(k,1204) = -rxt(k,290)*y(k,74)
         mat(k,989) = rxt(k,284)*y(k,136)
         mat(k,644) = rxt(k,284)*y(k,134)
         mat(k,275) = -(rxt(k,297)*y(k,139))
         mat(k,1222) = -rxt(k,297)*y(k,75)
         mat(k,897) = rxt(k,300)*y(k,137)
         mat(k,665) = rxt(k,300)*y(k,90)
         mat(k,738) = -(rxt(k,298)*y(k,102) + rxt(k,299)*y(k,139))
         mat(k,1071) = -rxt(k,298)*y(k,76)
         mat(k,1259) = -rxt(k,299)*y(k,76)
         mat(k,437) = .200_r8*rxt(k,313)*y(k,102)
         mat(k,364) = .039_r8*rxt(k,306)*y(k,89) + .039_r8*rxt(k,307)*y(k,91) &
                      + .039_r8*rxt(k,305)*y(k,134)
         mat(k,847) = .039_r8*rxt(k,306)*y(k,71) + .320_r8*rxt(k,311)*y(k,135)
         mat(k,1134) = .039_r8*rxt(k,307)*y(k,71) + .350_r8*rxt(k,312)*y(k,135)
         mat(k,1071) = mat(k,1071) + .200_r8*rxt(k,313)*y(k,70) + .442_r8*rxt(k,325) &
                      *y(k,126)
         mat(k,304) = .442_r8*rxt(k,325)*y(k,102)
         mat(k,701) = .350_r8*rxt(k,308)*y(k,135)
         mat(k,776) = .260_r8*rxt(k,309)*y(k,135)
         mat(k,1025) = .039_r8*rxt(k,305)*y(k,71)
         mat(k,632) = .320_r8*rxt(k,311)*y(k,89) + .350_r8*rxt(k,312)*y(k,91) &
                      + .350_r8*rxt(k,308)*y(k,130) + .260_r8*rxt(k,309)*y(k,131)
         mat(k,311) = -(rxt(k,122)*y(k,89) + (rxt(k,123) + rxt(k,124) + rxt(k,125) &
                      ) * y(k,90) + rxt(k,126)*y(k,99) + rxt(k,134)*y(k,139))
         mat(k,825) = -rxt(k,122)*y(k,77)
         mat(k,898) = -(rxt(k,123) + rxt(k,124) + rxt(k,125)) * y(k,77)
         mat(k,1336) = -rxt(k,126)*y(k,77)
         mat(k,1226) = -rxt(k,134)*y(k,77)
         mat(k,111) = -((rxt(k,138) + rxt(k,139)) * y(k,138))
         mat(k,1166) = -(rxt(k,138) + rxt(k,139)) * y(k,78)
         mat(k,310) = rxt(k,123)*y(k,90)
         mat(k,892) = rxt(k,123)*y(k,77)
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
         mat(k,893) = rxt(k,141)*y(k,91)
         mat(k,1107) = rxt(k,141)*y(k,90)
         mat(k,43) = -(rxt(k,344)*y(k,139))
         mat(k,1193) = -rxt(k,344)*y(k,84)
         mat(k,851) = -(rxt(k,122)*y(k,77) + rxt(k,131)*y(k,91) + rxt(k,135)*y(k,134) &
                      + rxt(k,136)*y(k,102) + rxt(k,137)*y(k,98) + rxt(k,158)*y(k,41) &
                      + rxt(k,190)*y(k,6) + rxt(k,226)*y(k,131) + rxt(k,239)*y(k,128) &
                      + rxt(k,249)*y(k,130) + rxt(k,253)*y(k,133) + rxt(k,266) &
                      *y(k,129) + rxt(k,274)*y(k,140) + rxt(k,278)*y(k,141) + (rxt(k,286) &
                      + rxt(k,287)) * y(k,136) + rxt(k,295)*y(k,137) + rxt(k,306) &
                      *y(k,71) + rxt(k,311)*y(k,135) + rxt(k,321)*y(k,142))
         mat(k,313) = -rxt(k,122)*y(k,89)
         mat(k,1138) = -rxt(k,131)*y(k,89)
         mat(k,1029) = -rxt(k,135)*y(k,89)
         mat(k,1075) = -rxt(k,136)*y(k,89)
         mat(k,1295) = -rxt(k,137)*y(k,89)
         mat(k,877) = -rxt(k,158)*y(k,89)
         mat(k,1318) = -rxt(k,190)*y(k,89)
         mat(k,780) = -rxt(k,226)*y(k,89)
         mat(k,389) = -rxt(k,239)*y(k,89)
         mat(k,704) = -rxt(k,249)*y(k,89)
         mat(k,356) = -rxt(k,253)*y(k,89)
         mat(k,423) = -rxt(k,266)*y(k,89)
         mat(k,376) = -rxt(k,274)*y(k,89)
         mat(k,490) = -rxt(k,278)*y(k,89)
         mat(k,657) = -(rxt(k,286) + rxt(k,287)) * y(k,89)
         mat(k,676) = -rxt(k,295)*y(k,89)
         mat(k,366) = -rxt(k,306)*y(k,89)
         mat(k,635) = -rxt(k,311)*y(k,89)
         mat(k,598) = -rxt(k,321)*y(k,89)
         mat(k,313) = mat(k,313) + 2.000_r8*rxt(k,124)*y(k,90) + rxt(k,126)*y(k,99) &
                      + rxt(k,134)*y(k,139)
         mat(k,112) = 2.000_r8*rxt(k,138)*y(k,138)
         mat(k,914) = 2.000_r8*rxt(k,124)*y(k,77) + rxt(k,127)*y(k,98) + rxt(k,337) &
                      *y(k,114)
         mat(k,1295) = mat(k,1295) + rxt(k,127)*y(k,90)
         mat(k,1344) = rxt(k,126)*y(k,77)
         mat(k,718) = rxt(k,337)*y(k,90)
         mat(k,1177) = 2.000_r8*rxt(k,138)*y(k,78)
         mat(k,1263) = rxt(k,134)*y(k,77)
         mat(k,916) = -((rxt(k,123) + rxt(k,124) + rxt(k,125)) * y(k,77) + (rxt(k,127) &
                      + rxt(k,129)) * y(k,98) + rxt(k,128)*y(k,102) + rxt(k,140) &
                      *y(k,134) + rxt(k,141)*y(k,91) + rxt(k,142)*y(k,139) + rxt(k,160) &
                      *y(k,41) + rxt(k,191)*y(k,6) + rxt(k,260)*y(k,130) + rxt(k,300) &
                      *y(k,137) + rxt(k,337)*y(k,114))
         mat(k,314) = -(rxt(k,123) + rxt(k,124) + rxt(k,125)) * y(k,90)
         mat(k,1297) = -(rxt(k,127) + rxt(k,129)) * y(k,90)
         mat(k,1077) = -rxt(k,128)*y(k,90)
         mat(k,1031) = -rxt(k,140)*y(k,90)
         mat(k,1140) = -rxt(k,141)*y(k,90)
         mat(k,1265) = -rxt(k,142)*y(k,90)
         mat(k,879) = -rxt(k,160)*y(k,90)
         mat(k,1320) = -rxt(k,191)*y(k,90)
         mat(k,705) = -rxt(k,260)*y(k,90)
         mat(k,677) = -rxt(k,300)*y(k,90)
         mat(k,720) = -rxt(k,337)*y(k,90)
         mat(k,1320) = mat(k,1320) + rxt(k,190)*y(k,89)
         mat(k,879) = mat(k,879) + rxt(k,158)*y(k,89)
         mat(k,219) = rxt(k,121)*y(k,139)
         mat(k,367) = 1.206_r8*rxt(k,306)*y(k,89) + 1.206_r8*rxt(k,307)*y(k,91) &
                      + .206_r8*rxt(k,305)*y(k,134)
         mat(k,853) = rxt(k,190)*y(k,6) + rxt(k,158)*y(k,41) + 1.206_r8*rxt(k,306) &
                      *y(k,71) + 2.000_r8*rxt(k,131)*y(k,91) + rxt(k,137)*y(k,98) &
                      + rxt(k,136)*y(k,102) + rxt(k,239)*y(k,128) + rxt(k,266) &
                      *y(k,129) + rxt(k,249)*y(k,130) + rxt(k,226)*y(k,131) &
                      + rxt(k,253)*y(k,133) + rxt(k,135)*y(k,134) + .920_r8*rxt(k,311) &
                      *y(k,135) + rxt(k,286)*y(k,136) + rxt(k,295)*y(k,137) &
                      + rxt(k,274)*y(k,140) + rxt(k,278)*y(k,141) + rxt(k,321) &
                      *y(k,142)
         mat(k,1140) = mat(k,1140) + 1.206_r8*rxt(k,307)*y(k,71) + 2.000_r8*rxt(k,131) &
                      *y(k,89) + rxt(k,132)*y(k,98) + rxt(k,316)*y(k,106) + rxt(k,324) &
                      *y(k,126) + rxt(k,130)*y(k,134) + rxt(k,312)*y(k,135) &
                      + rxt(k,285)*y(k,136) + rxt(k,296)*y(k,137) + rxt(k,133) &
                      *y(k,139) + rxt(k,322)*y(k,142)
         mat(k,186) = rxt(k,272)*y(k,139)
         mat(k,1297) = mat(k,1297) + rxt(k,137)*y(k,89) + rxt(k,132)*y(k,91)
         mat(k,1077) = mat(k,1077) + rxt(k,136)*y(k,89)
         mat(k,477) = rxt(k,316)*y(k,91) + .400_r8*rxt(k,317)*y(k,139)
         mat(k,305) = rxt(k,324)*y(k,91)
         mat(k,390) = rxt(k,239)*y(k,89)
         mat(k,424) = rxt(k,266)*y(k,89)
         mat(k,705) = mat(k,705) + rxt(k,249)*y(k,89)
         mat(k,782) = rxt(k,226)*y(k,89)
         mat(k,357) = rxt(k,253)*y(k,89)
         mat(k,1031) = mat(k,1031) + .206_r8*rxt(k,305)*y(k,71) + rxt(k,135)*y(k,89) &
                      + rxt(k,130)*y(k,91)
         mat(k,636) = .920_r8*rxt(k,311)*y(k,89) + rxt(k,312)*y(k,91)
         mat(k,658) = rxt(k,286)*y(k,89) + rxt(k,285)*y(k,91)
         mat(k,677) = mat(k,677) + rxt(k,295)*y(k,89) + rxt(k,296)*y(k,91)
         mat(k,1265) = mat(k,1265) + rxt(k,121)*y(k,65) + rxt(k,133)*y(k,91) &
                      + rxt(k,272)*y(k,92) + .400_r8*rxt(k,317)*y(k,106)
         mat(k,377) = rxt(k,274)*y(k,89)
         mat(k,491) = rxt(k,278)*y(k,89)
         mat(k,599) = rxt(k,321)*y(k,89) + rxt(k,322)*y(k,91)
         mat(k,1146) = -(rxt(k,130)*y(k,134) + rxt(k,131)*y(k,89) + rxt(k,132)*y(k,98) &
                      + rxt(k,133)*y(k,139) + rxt(k,141)*y(k,90) + rxt(k,220)*y(k,25) &
                      + rxt(k,244)*y(k,28) + rxt(k,262)*y(k,13) + rxt(k,269)*y(k,31) &
                      + rxt(k,285)*y(k,136) + rxt(k,296)*y(k,137) + rxt(k,304)*y(k,70) &
                      + rxt(k,307)*y(k,71) + rxt(k,312)*y(k,135) + rxt(k,316)*y(k,106) &
                      + rxt(k,322)*y(k,142) + rxt(k,324)*y(k,126) + rxt(k,327)*y(k,46))
         mat(k,1037) = -rxt(k,130)*y(k,91)
         mat(k,859) = -rxt(k,131)*y(k,91)
         mat(k,1303) = -rxt(k,132)*y(k,91)
         mat(k,1271) = -rxt(k,133)*y(k,91)
         mat(k,922) = -rxt(k,141)*y(k,91)
         mat(k,816) = -rxt(k,220)*y(k,91)
         mat(k,531) = -rxt(k,244)*y(k,91)
         mat(k,554) = -rxt(k,262)*y(k,91)
         mat(k,611) = -rxt(k,269)*y(k,91)
         mat(k,660) = -rxt(k,285)*y(k,91)
         mat(k,680) = -rxt(k,296)*y(k,91)
         mat(k,444) = -rxt(k,304)*y(k,91)
         mat(k,369) = -rxt(k,307)*y(k,91)
         mat(k,639) = -rxt(k,312)*y(k,91)
         mat(k,479) = -rxt(k,316)*y(k,91)
         mat(k,601) = -rxt(k,322)*y(k,91)
         mat(k,307) = -rxt(k,324)*y(k,91)
         mat(k,140) = -rxt(k,327)*y(k,91)
         mat(k,237) = rxt(k,192)*y(k,98)
         mat(k,980) = rxt(k,159)*y(k,42)
         mat(k,454) = rxt(k,159)*y(k,38) + rxt(k,161)*y(k,98) + rxt(k,162)*y(k,139)
         mat(k,412) = rxt(k,143)*y(k,139)
         mat(k,282) = .500_r8*rxt(k,297)*y(k,139)
         mat(k,922) = mat(k,922) + rxt(k,129)*y(k,98) + rxt(k,128)*y(k,102)
         mat(k,1303) = mat(k,1303) + rxt(k,192)*y(k,7) + rxt(k,161)*y(k,42) &
                      + rxt(k,129)*y(k,90)
         mat(k,1083) = rxt(k,128)*y(k,90)
         mat(k,258) = rxt(k,258)*y(k,139)
         mat(k,1271) = mat(k,1271) + rxt(k,162)*y(k,42) + rxt(k,143)*y(k,64) &
                      + .500_r8*rxt(k,297)*y(k,75) + rxt(k,258)*y(k,107)
         mat(k,182) = -(rxt(k,272)*y(k,139))
         mat(k,1210) = -rxt(k,272)*y(k,92)
         mat(k,534) = rxt(k,262)*y(k,91)
         mat(k,1110) = rxt(k,262)*y(k,13)
         mat(k,1306) = -(rxt(k,101)*y(k,102) + 4._r8*rxt(k,102)*y(k,98) + rxt(k,103) &
                      *y(k,99) + rxt(k,104)*y(k,55) + rxt(k,105)*y(k,57) + rxt(k,110) &
                      *y(k,134) + rxt(k,116)*y(k,139) + (rxt(k,127) + rxt(k,129) &
                      ) * y(k,90) + rxt(k,132)*y(k,91) + rxt(k,137)*y(k,89) + rxt(k,161) &
                      *y(k,42) + rxt(k,163)*y(k,41) + rxt(k,166)*y(k,63) + rxt(k,169) &
                      *y(k,67) + rxt(k,192)*y(k,7) + rxt(k,193)*y(k,6) + rxt(k,195) &
                      *y(k,59) + rxt(k,197)*y(k,66) + rxt(k,221)*y(k,25) + rxt(k,329) &
                      *y(k,105))
         mat(k,1086) = -rxt(k,101)*y(k,98)
         mat(k,1355) = -rxt(k,103)*y(k,98)
         mat(k,501) = -rxt(k,104)*y(k,98)
         mat(k,273) = -rxt(k,105)*y(k,98)
         mat(k,1040) = -rxt(k,110)*y(k,98)
         mat(k,1274) = -rxt(k,116)*y(k,98)
         mat(k,925) = -(rxt(k,127) + rxt(k,129)) * y(k,98)
         mat(k,1149) = -rxt(k,132)*y(k,98)
         mat(k,862) = -rxt(k,137)*y(k,98)
         mat(k,456) = -rxt(k,161)*y(k,98)
         mat(k,888) = -rxt(k,163)*y(k,98)
         mat(k,946) = -rxt(k,166)*y(k,98)
         mat(k,400) = -rxt(k,169)*y(k,98)
         mat(k,238) = -rxt(k,192)*y(k,98)
         mat(k,1329) = -rxt(k,193)*y(k,98)
         mat(k,407) = -rxt(k,195)*y(k,98)
         mat(k,342) = -rxt(k,197)*y(k,98)
         mat(k,819) = -rxt(k,221)*y(k,98)
         mat(k,167) = -rxt(k,329)*y(k,98)
         mat(k,801) = rxt(k,108)*y(k,134)
         mat(k,317) = rxt(k,122)*y(k,89) + rxt(k,123)*y(k,90) + rxt(k,126)*y(k,99)
         mat(k,862) = mat(k,862) + rxt(k,122)*y(k,77)
         mat(k,925) = mat(k,925) + rxt(k,123)*y(k,77)
         mat(k,1355) = mat(k,1355) + rxt(k,126)*y(k,77) + rxt(k,331)*y(k,112) &
                      + rxt(k,338)*y(k,114) + (rxt(k,88)+rxt(k,89)+rxt(k,90))*y(k,138)
         mat(k,1086) = mat(k,1086) + .765_r8*rxt(k,325)*y(k,126) + 2.000_r8*rxt(k,92) &
                      *y(k,138)
         mat(k,329) = rxt(k,331)*y(k,99)
         mat(k,726) = rxt(k,338)*y(k,99)
         mat(k,309) = .765_r8*rxt(k,325)*y(k,102)
         mat(k,1040) = mat(k,1040) + rxt(k,108)*y(k,54)
         mat(k,1188) = (rxt(k,88)+rxt(k,89)+rxt(k,90))*y(k,99) + 2.000_r8*rxt(k,92) &
                      *y(k,102)
         mat(k,1274) = mat(k,1274) + 2.000_r8*rxt(k,118)*y(k,139)
         mat(k,1357) = -(rxt(k,88)*y(k,138) + rxt(k,95)*y(k,100) + rxt(k,103)*y(k,98) &
                      + rxt(k,109)*y(k,54) + rxt(k,126)*y(k,77) + rxt(k,255)*y(k,132) &
                      + rxt(k,331)*y(k,112) + rxt(k,338)*y(k,114))
         mat(k,1190) = -rxt(k,88)*y(k,99)
         mat(k,88) = -rxt(k,95)*y(k,99)
         mat(k,1308) = -rxt(k,103)*y(k,99)
         mat(k,802) = -rxt(k,109)*y(k,99)
         mat(k,318) = -rxt(k,126)*y(k,99)
         mat(k,210) = -rxt(k,255)*y(k,99)
         mat(k,330) = -rxt(k,331)*y(k,99)
         mat(k,728) = -rxt(k,338)*y(k,99)
         mat(k,1105) = rxt(k,184)*y(k,102) + rxt(k,183)*y(k,134)
         mat(k,1331) = 2.000_r8*rxt(k,185)*y(k,6) + (rxt(k,187)+rxt(k,188))*y(k,41) &
                      + rxt(k,193)*y(k,98) + rxt(k,189)*y(k,134)
         mat(k,985) = rxt(k,152)*y(k,102) + rxt(k,150)*y(k,134)
         mat(k,890) = (rxt(k,187)+rxt(k,188))*y(k,6) + (2.000_r8*rxt(k,154) &
                       +2.000_r8*rxt(k,155))*y(k,41) + rxt(k,163)*y(k,98) + rxt(k,157) &
                      *y(k,134) + rxt(k,165)*y(k,139)
         mat(k,802) = mat(k,802) + rxt(k,112)*y(k,102) + rxt(k,106)*y(k,134)
         mat(k,223) = rxt(k,121)*y(k,139)
         mat(k,318) = mat(k,318) + rxt(k,125)*y(k,90)
         mat(k,114) = rxt(k,139)*y(k,138)
         mat(k,864) = rxt(k,136)*y(k,102)
         mat(k,927) = rxt(k,125)*y(k,77) + rxt(k,127)*y(k,98) + rxt(k,128)*y(k,102)
         mat(k,1151) = rxt(k,132)*y(k,98) + rxt(k,130)*y(k,134)
         mat(k,1308) = mat(k,1308) + rxt(k,193)*y(k,6) + rxt(k,163)*y(k,41) &
                      + rxt(k,127)*y(k,90) + rxt(k,132)*y(k,91) + 2.000_r8*rxt(k,102) &
                      *y(k,98) + rxt(k,94)*y(k,100) + 2.000_r8*rxt(k,101)*y(k,102) &
                      + rxt(k,110)*y(k,134) + rxt(k,116)*y(k,139)
         mat(k,1357) = mat(k,1357) + 2.000_r8*rxt(k,95)*y(k,100)
         mat(k,88) = mat(k,88) + rxt(k,94)*y(k,98) + 2.000_r8*rxt(k,95)*y(k,99)
         mat(k,1088) = rxt(k,184)*y(k,4) + rxt(k,152)*y(k,38) + rxt(k,112)*y(k,54) &
                      + rxt(k,136)*y(k,89) + rxt(k,128)*y(k,90) + 2.000_r8*rxt(k,101) &
                      *y(k,98) + rxt(k,333)*y(k,112) + rxt(k,339)*y(k,114) &
                      + 2.000_r8*rxt(k,111)*y(k,134) + (2.000_r8*rxt(k,91)+rxt(k,92)) &
                      *y(k,138) + rxt(k,117)*y(k,139)
         mat(k,330) = mat(k,330) + rxt(k,333)*y(k,102)
         mat(k,728) = mat(k,728) + rxt(k,339)*y(k,102)
         mat(k,393) = rxt(k,238)*y(k,134)
         mat(k,427) = rxt(k,265)*y(k,134)
         mat(k,791) = rxt(k,225)*y(k,134)
         mat(k,1042) = rxt(k,183)*y(k,4) + rxt(k,189)*y(k,6) + rxt(k,150)*y(k,38) &
                      + rxt(k,157)*y(k,41) + rxt(k,106)*y(k,54) + rxt(k,130)*y(k,91) &
                      + rxt(k,110)*y(k,98) + 2.000_r8*rxt(k,111)*y(k,102) + rxt(k,238) &
                      *y(k,128) + rxt(k,265)*y(k,129) + rxt(k,225)*y(k,131) &
                      + 2.000_r8*rxt(k,120)*y(k,134) + rxt(k,115)*y(k,139) &
                      + rxt(k,273)*y(k,140)
         mat(k,1190) = mat(k,1190) + rxt(k,139)*y(k,78) + (2.000_r8*rxt(k,91) &
                       +rxt(k,92))*y(k,102)
         mat(k,1276) = rxt(k,165)*y(k,41) + rxt(k,121)*y(k,65) + rxt(k,116)*y(k,98) &
                      + rxt(k,117)*y(k,102) + rxt(k,115)*y(k,134)
         mat(k,380) = rxt(k,273)*y(k,134)
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
         mat(k,87) = -(rxt(k,94)*y(k,98) + rxt(k,95)*y(k,99))
         mat(k,1278) = -rxt(k,94)*y(k,100)
         mat(k,1333) = -rxt(k,95)*y(k,100)
         mat(k,458) = rxt(k,96)*y(k,101)
         mat(k,1278) = mat(k,1278) + rxt(k,98)*y(k,101)
         mat(k,1333) = mat(k,1333) + rxt(k,99)*y(k,101)
         mat(k,89) = rxt(k,96)*y(k,45) + rxt(k,98)*y(k,98) + rxt(k,99)*y(k,99) &
                      + rxt(k,100)*y(k,102)
         mat(k,1044) = rxt(k,100)*y(k,101)
         mat(k,90) = -(rxt(k,96)*y(k,45) + rxt(k,98)*y(k,98) + rxt(k,99)*y(k,99) &
                      + rxt(k,100)*y(k,102))
         mat(k,459) = -rxt(k,96)*y(k,101)
         mat(k,1279) = -rxt(k,98)*y(k,101)
         mat(k,1334) = -rxt(k,99)*y(k,101)
         mat(k,1045) = -rxt(k,100)*y(k,101)
         mat(k,1334) = mat(k,1334) + rxt(k,88)*y(k,138)
         mat(k,1162) = rxt(k,88)*y(k,99)
         mat(k,1081) = -((rxt(k,91) + rxt(k,92)) * y(k,138) + rxt(k,101)*y(k,98) &
                      + rxt(k,111)*y(k,134) + rxt(k,112)*y(k,54) + rxt(k,117)*y(k,139) &
                      + rxt(k,128)*y(k,90) + rxt(k,136)*y(k,89) + rxt(k,152)*y(k,38) &
                      + rxt(k,184)*y(k,4) + rxt(k,235)*y(k,9) + rxt(k,263)*y(k,13) &
                      + rxt(k,288)*y(k,73) + rxt(k,298)*y(k,76) + rxt(k,313)*y(k,70) &
                      + rxt(k,325)*y(k,126) + rxt(k,333)*y(k,112) + rxt(k,339) &
                      *y(k,114))
         mat(k,1183) = -(rxt(k,91) + rxt(k,92)) * y(k,102)
         mat(k,1301) = -rxt(k,101)*y(k,102)
         mat(k,1035) = -rxt(k,111)*y(k,102)
         mat(k,798) = -rxt(k,112)*y(k,102)
         mat(k,1269) = -rxt(k,117)*y(k,102)
         mat(k,920) = -rxt(k,128)*y(k,102)
         mat(k,857) = -rxt(k,136)*y(k,102)
         mat(k,978) = -rxt(k,152)*y(k,102)
         mat(k,1098) = -rxt(k,184)*y(k,102)
         mat(k,216) = -rxt(k,235)*y(k,102)
         mat(k,553) = -rxt(k,263)*y(k,102)
         mat(k,577) = -rxt(k,288)*y(k,102)
         mat(k,747) = -rxt(k,298)*y(k,102)
         mat(k,443) = -rxt(k,313)*y(k,102)
         mat(k,306) = -rxt(k,325)*y(k,102)
         mat(k,327) = -rxt(k,333)*y(k,102)
         mat(k,723) = -rxt(k,339)*y(k,102)
         mat(k,1301) = mat(k,1301) + rxt(k,103)*y(k,99)
         mat(k,1350) = rxt(k,103)*y(k,98)
         mat(k,707) = .150_r8*rxt(k,248)*y(k,134)
         mat(k,1035) = mat(k,1035) + .150_r8*rxt(k,248)*y(k,130) + .150_r8*rxt(k,293) &
                      *y(k,137)
         mat(k,679) = .150_r8*rxt(k,293)*y(k,134)
         mat(k,152) = -(rxt(k,340)*y(k,114))
         mat(k,713) = -rxt(k,340)*y(k,104)
         mat(k,1311) = rxt(k,186)*y(k,41)
         mat(k,869) = rxt(k,186)*y(k,6) + 2.000_r8*rxt(k,156)*y(k,41)
         mat(k,160) = -(rxt(k,329)*y(k,98) + rxt(k,330)*y(k,139))
         mat(k,1280) = -rxt(k,329)*y(k,105)
         mat(k,1207) = -rxt(k,330)*y(k,105)
         mat(k,473) = -(rxt(k,316)*y(k,91) + rxt(k,317)*y(k,139))
         mat(k,1119) = -rxt(k,316)*y(k,106)
         mat(k,1242) = -rxt(k,317)*y(k,106)
         mat(k,362) = .794_r8*rxt(k,306)*y(k,89) + .794_r8*rxt(k,307)*y(k,91) &
                      + .794_r8*rxt(k,305)*y(k,134)
         mat(k,834) = .794_r8*rxt(k,306)*y(k,71) + .080_r8*rxt(k,311)*y(k,135) &
                      + .800_r8*rxt(k,287)*y(k,136)
         mat(k,1119) = mat(k,1119) + .794_r8*rxt(k,307)*y(k,71)
         mat(k,1011) = .794_r8*rxt(k,305)*y(k,71)
         mat(k,621) = .080_r8*rxt(k,311)*y(k,89)
         mat(k,646) = .800_r8*rxt(k,287)*y(k,89)
         mat(k,252) = -(rxt(k,258)*y(k,139))
         mat(k,1219) = -rxt(k,258)*y(k,107)
         mat(k,896) = rxt(k,260)*y(k,130)
         mat(k,686) = rxt(k,260)*y(k,90)
         mat(k,260) = -(rxt(k,275)*y(k,139))
         mat(k,1220) = -rxt(k,275)*y(k,110)
         mat(k,997) = rxt(k,273)*y(k,140)
         mat(k,371) = rxt(k,273)*y(k,134)
         mat(k,200) = -(rxt(k,279)*y(k,139))
         mat(k,1213) = -rxt(k,279)*y(k,111)
         mat(k,993) = .850_r8*rxt(k,277)*y(k,141)
         mat(k,482) = .850_r8*rxt(k,277)*y(k,134)
         mat(k,324) = -(rxt(k,331)*y(k,99) + rxt(k,333)*y(k,102) + rxt(k,336)*y(k,139))
         mat(k,1337) = -rxt(k,331)*y(k,112)
         mat(k,1049) = -rxt(k,333)*y(k,112)
         mat(k,1228) = -rxt(k,336)*y(k,112)
         mat(k,716) = -(rxt(k,334)*y(k,6) + rxt(k,335)*y(k,41) + rxt(k,337)*y(k,90) &
                      + rxt(k,338)*y(k,99) + rxt(k,339)*y(k,102) + rxt(k,340)*y(k,104) &
                      + rxt(k,341)*y(k,139))
         mat(k,1316) = -rxt(k,334)*y(k,114)
         mat(k,873) = -rxt(k,335)*y(k,114)
         mat(k,909) = -rxt(k,337)*y(k,114)
         mat(k,1341) = -rxt(k,338)*y(k,114)
         mat(k,1070) = -rxt(k,339)*y(k,114)
         mat(k,154) = -rxt(k,340)*y(k,114)
         mat(k,1258) = -rxt(k,341)*y(k,114)
         mat(k,1292) = rxt(k,329)*y(k,105)
         mat(k,1341) = mat(k,1341) + rxt(k,331)*y(k,112)
         mat(k,1070) = mat(k,1070) + rxt(k,333)*y(k,112)
         mat(k,164) = rxt(k,329)*y(k,98)
         mat(k,325) = rxt(k,331)*y(k,99) + rxt(k,333)*y(k,102) + rxt(k,336)*y(k,139)
         mat(k,1258) = mat(k,1258) + rxt(k,336)*y(k,112)
         mat(k,465) = -(rxt(k,332)*y(k,139))
         mat(k,1241) = -rxt(k,332)*y(k,115)
         mat(k,1315) = rxt(k,334)*y(k,114)
         mat(k,872) = rxt(k,335)*y(k,114)
         mat(k,138) = rxt(k,327)*y(k,91) + (rxt(k,328)+.500_r8*rxt(k,342))*y(k,139)
         mat(k,903) = rxt(k,337)*y(k,114)
         mat(k,1118) = rxt(k,327)*y(k,46)
         mat(k,1338) = rxt(k,338)*y(k,114)
         mat(k,1054) = rxt(k,339)*y(k,114)
         mat(k,153) = rxt(k,340)*y(k,114)
         mat(k,162) = rxt(k,330)*y(k,139)
         mat(k,715) = rxt(k,334)*y(k,6) + rxt(k,335)*y(k,41) + rxt(k,337)*y(k,90) &
                      + rxt(k,338)*y(k,99) + rxt(k,339)*y(k,102) + rxt(k,340)*y(k,104) &
                      + rxt(k,341)*y(k,139)
         mat(k,1241) = mat(k,1241) + (rxt(k,328)+.500_r8*rxt(k,342))*y(k,46) &
                      + rxt(k,330)*y(k,105) + rxt(k,341)*y(k,114)
         mat(k,120) = -(rxt(k,343)*y(k,143))
         mat(k,1360) = -rxt(k,343)*y(k,116)
         mat(k,464) = rxt(k,332)*y(k,139)
         mat(k,1201) = rxt(k,332)*y(k,115)
         mat(k,300) = -(rxt(k,324)*y(k,91) + rxt(k,325)*y(k,102) + rxt(k,326)*y(k,139))
         mat(k,1111) = -rxt(k,324)*y(k,126)
         mat(k,1047) = -rxt(k,325)*y(k,126)
         mat(k,1225) = -rxt(k,326)*y(k,126)
         mat(k,100) = -(rxt(k,323)*y(k,139))
         mat(k,1198) = -rxt(k,323)*y(k,127)
         mat(k,988) = rxt(k,320)*y(k,142)
         mat(k,586) = rxt(k,320)*y(k,134)
         mat(k,385) = -(4._r8*rxt(k,236)*y(k,128) + rxt(k,237)*y(k,131) + rxt(k,238) &
                      *y(k,134) + rxt(k,239)*y(k,89))
         mat(k,761) = -rxt(k,237)*y(k,128)
         mat(k,1006) = -rxt(k,238)*y(k,128)
         mat(k,831) = -rxt(k,239)*y(k,128)
         mat(k,148) = .500_r8*rxt(k,241)*y(k,139)
         mat(k,131) = rxt(k,242)*y(k,38) + rxt(k,243)*y(k,139)
         mat(k,959) = rxt(k,242)*y(k,12)
         mat(k,1233) = .500_r8*rxt(k,241)*y(k,11) + rxt(k,243)*y(k,12)
         mat(k,417) = -(rxt(k,264)*y(k,131) + rxt(k,265)*y(k,134) + rxt(k,266)*y(k,89))
         mat(k,762) = -rxt(k,264)*y(k,129)
         mat(k,1009) = -rxt(k,265)*y(k,129)
         mat(k,832) = -rxt(k,266)*y(k,129)
         mat(k,38) = 1.670_r8*rxt(k,302)*y(k,139)
         mat(k,178) = rxt(k,267)*y(k,139)
         mat(k,68) = rxt(k,268)*y(k,139)
         mat(k,1237) = 1.670_r8*rxt(k,302)*y(k,3) + rxt(k,267)*y(k,14) + rxt(k,268) &
                      *y(k,15)
         mat(k,700) = -(4._r8*rxt(k,246)*y(k,130) + rxt(k,247)*y(k,131) + rxt(k,248) &
                      *y(k,134) + rxt(k,249)*y(k,89) + rxt(k,260)*y(k,90) + rxt(k,282) &
                      *y(k,136) + rxt(k,308)*y(k,135) + rxt(k,318)*y(k,142))
         mat(k,775) = -rxt(k,247)*y(k,130)
         mat(k,1024) = -rxt(k,248)*y(k,130)
         mat(k,846) = -rxt(k,249)*y(k,130)
         mat(k,908) = -rxt(k,260)*y(k,130)
         mat(k,654) = -rxt(k,282)*y(k,130)
         mat(k,631) = -rxt(k,308)*y(k,130)
         mat(k,595) = -rxt(k,318)*y(k,130)
         mat(k,527) = rxt(k,244)*y(k,91) + rxt(k,245)*y(k,139)
         mat(k,608) = rxt(k,269)*y(k,91) + rxt(k,270)*y(k,139)
         mat(k,286) = .500_r8*rxt(k,251)*y(k,139)
         mat(k,436) = .080_r8*rxt(k,313)*y(k,102)
         mat(k,574) = .100_r8*rxt(k,288)*y(k,102)
         mat(k,737) = .280_r8*rxt(k,298)*y(k,102)
         mat(k,846) = mat(k,846) + .530_r8*rxt(k,286)*y(k,136) + rxt(k,295)*y(k,137) &
                      + rxt(k,278)*y(k,141)
         mat(k,1132) = rxt(k,244)*y(k,28) + rxt(k,269)*y(k,31) + .530_r8*rxt(k,285) &
                      *y(k,136) + rxt(k,296)*y(k,137)
         mat(k,1069) = .080_r8*rxt(k,313)*y(k,70) + .100_r8*rxt(k,288)*y(k,73) &
                      + .280_r8*rxt(k,298)*y(k,76)
         mat(k,700) = mat(k,700) + .530_r8*rxt(k,282)*y(k,136)
         mat(k,775) = mat(k,775) + .260_r8*rxt(k,283)*y(k,136) + rxt(k,292)*y(k,137) &
                      + .300_r8*rxt(k,276)*y(k,141)
         mat(k,1024) = mat(k,1024) + .450_r8*rxt(k,293)*y(k,137) + .150_r8*rxt(k,277) &
                      *y(k,141)
         mat(k,654) = mat(k,654) + .530_r8*rxt(k,286)*y(k,89) + .530_r8*rxt(k,285) &
                      *y(k,91) + .530_r8*rxt(k,282)*y(k,130) + .260_r8*rxt(k,283) &
                      *y(k,131)
         mat(k,673) = rxt(k,295)*y(k,89) + rxt(k,296)*y(k,91) + rxt(k,292)*y(k,131) &
                      + .450_r8*rxt(k,293)*y(k,134) + 4.000_r8*rxt(k,294)*y(k,137)
         mat(k,1257) = rxt(k,245)*y(k,28) + rxt(k,270)*y(k,31) + .500_r8*rxt(k,251) &
                      *y(k,33)
         mat(k,487) = rxt(k,278)*y(k,89) + .300_r8*rxt(k,276)*y(k,131) &
                      + .150_r8*rxt(k,277)*y(k,134)
         mat(k,777) = -(rxt(k,153)*y(k,41) + (4._r8*rxt(k,223) + 4._r8*rxt(k,224) &
                      ) * y(k,131) + rxt(k,225)*y(k,134) + rxt(k,226)*y(k,89) &
                      + rxt(k,237)*y(k,128) + rxt(k,247)*y(k,130) + rxt(k,264) &
                      *y(k,129) + rxt(k,276)*y(k,141) + rxt(k,283)*y(k,136) + rxt(k,292) &
                      *y(k,137) + rxt(k,309)*y(k,135) + rxt(k,319)*y(k,142))
         mat(k,874) = -rxt(k,153)*y(k,131)
         mat(k,1026) = -rxt(k,225)*y(k,131)
         mat(k,848) = -rxt(k,226)*y(k,131)
         mat(k,387) = -rxt(k,237)*y(k,131)
         mat(k,702) = -rxt(k,247)*y(k,131)
         mat(k,421) = -rxt(k,264)*y(k,131)
         mat(k,488) = -rxt(k,276)*y(k,131)
         mat(k,655) = -rxt(k,283)*y(k,131)
         mat(k,674) = -rxt(k,292)*y(k,131)
         mat(k,633) = -rxt(k,309)*y(k,131)
         mat(k,596) = -rxt(k,319)*y(k,131)
         mat(k,545) = .280_r8*rxt(k,263)*y(k,102)
         mat(k,321) = rxt(k,250)*y(k,139)
         mat(k,195) = .700_r8*rxt(k,228)*y(k,139)
         mat(k,507) = rxt(k,147)*y(k,38) + rxt(k,230)*y(k,138) + rxt(k,229)*y(k,139)
         mat(k,969) = rxt(k,147)*y(k,36)
         mat(k,438) = .050_r8*rxt(k,313)*y(k,102)
         mat(k,848) = mat(k,848) + rxt(k,249)*y(k,130)
         mat(k,1072) = .280_r8*rxt(k,263)*y(k,13) + .050_r8*rxt(k,313)*y(k,70)
         mat(k,702) = mat(k,702) + rxt(k,249)*y(k,89) + 4.000_r8*rxt(k,246)*y(k,130) &
                      + .900_r8*rxt(k,247)*y(k,131) + .450_r8*rxt(k,248)*y(k,134) &
                      + rxt(k,308)*y(k,135) + rxt(k,282)*y(k,136) + rxt(k,291) &
                      *y(k,137) + rxt(k,318)*y(k,142)
         mat(k,777) = mat(k,777) + .900_r8*rxt(k,247)*y(k,130)
         mat(k,1026) = mat(k,1026) + .450_r8*rxt(k,248)*y(k,130)
         mat(k,633) = mat(k,633) + rxt(k,308)*y(k,130)
         mat(k,655) = mat(k,655) + rxt(k,282)*y(k,130)
         mat(k,674) = mat(k,674) + rxt(k,291)*y(k,130)
         mat(k,1174) = rxt(k,230)*y(k,36)
         mat(k,1260) = rxt(k,250)*y(k,32) + .700_r8*rxt(k,228)*y(k,35) + rxt(k,229) &
                      *y(k,36)
         mat(k,596) = mat(k,596) + rxt(k,318)*y(k,130)
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
         mat(k,206) = -(rxt(k,255)*y(k,99))
         mat(k,1335) = -rxt(k,255)*y(k,132)
         mat(k,824) = .750_r8*rxt(k,253)*y(k,133)
         mat(k,352) = .750_r8*rxt(k,253)*y(k,89)
         mat(k,353) = -(rxt(k,252)*y(k,134) + rxt(k,253)*y(k,89))
         mat(k,1003) = -rxt(k,252)*y(k,133)
         mat(k,828) = -rxt(k,253)*y(k,133)
         mat(k,212) = rxt(k,259)*y(k,139)
         mat(k,1231) = rxt(k,259)*y(k,9)
         mat(k,1034) = -((rxt(k,106) + rxt(k,107) + rxt(k,108)) * y(k,54) + rxt(k,110) &
                      *y(k,98) + rxt(k,111)*y(k,102) + rxt(k,115)*y(k,139) &
                      + 4._r8*rxt(k,120)*y(k,134) + rxt(k,130)*y(k,91) + rxt(k,135) &
                      *y(k,89) + rxt(k,140)*y(k,90) + (rxt(k,150) + rxt(k,151) &
                      ) * y(k,38) + rxt(k,157)*y(k,41) + rxt(k,183)*y(k,4) + rxt(k,189) &
                      *y(k,6) + rxt(k,225)*y(k,131) + rxt(k,238)*y(k,128) + rxt(k,248) &
                      *y(k,130) + rxt(k,252)*y(k,133) + rxt(k,265)*y(k,129) + rxt(k,273) &
                      *y(k,140) + rxt(k,277)*y(k,141) + rxt(k,284)*y(k,136) + rxt(k,293) &
                      *y(k,137) + rxt(k,305)*y(k,71) + rxt(k,310)*y(k,135) + rxt(k,320) &
                      *y(k,142))
         mat(k,797) = -(rxt(k,106) + rxt(k,107) + rxt(k,108)) * y(k,134)
         mat(k,1300) = -rxt(k,110)*y(k,134)
         mat(k,1080) = -rxt(k,111)*y(k,134)
         mat(k,1268) = -rxt(k,115)*y(k,134)
         mat(k,1143) = -rxt(k,130)*y(k,134)
         mat(k,856) = -rxt(k,135)*y(k,134)
         mat(k,919) = -rxt(k,140)*y(k,134)
         mat(k,977) = -(rxt(k,150) + rxt(k,151)) * y(k,134)
         mat(k,882) = -rxt(k,157)*y(k,134)
         mat(k,1097) = -rxt(k,183)*y(k,134)
         mat(k,1323) = -rxt(k,189)*y(k,134)
         mat(k,785) = -rxt(k,225)*y(k,134)
         mat(k,391) = -rxt(k,238)*y(k,134)
         mat(k,706) = -rxt(k,248)*y(k,134)
         mat(k,358) = -rxt(k,252)*y(k,134)
         mat(k,425) = -rxt(k,265)*y(k,134)
         mat(k,378) = -rxt(k,273)*y(k,134)
         mat(k,492) = -rxt(k,277)*y(k,134)
         mat(k,659) = -rxt(k,284)*y(k,134)
         mat(k,678) = -rxt(k,293)*y(k,134)
         mat(k,368) = -rxt(k,305)*y(k,134)
         mat(k,637) = -rxt(k,310)*y(k,134)
         mat(k,600) = -rxt(k,320)*y(k,134)
         mat(k,1097) = mat(k,1097) + rxt(k,182)*y(k,25)
         mat(k,1323) = mat(k,1323) + rxt(k,194)*y(k,139)
         mat(k,215) = .130_r8*rxt(k,235)*y(k,102)
         mat(k,117) = rxt(k,240)*y(k,139)
         mat(k,552) = .280_r8*rxt(k,263)*y(k,102)
         mat(k,813) = rxt(k,182)*y(k,4) + rxt(k,146)*y(k,38) + rxt(k,220)*y(k,91) &
                      + rxt(k,221)*y(k,98)
         mat(k,295) = rxt(k,205)*y(k,38) + rxt(k,206)*y(k,139)
         mat(k,173) = rxt(k,208)*y(k,38) + rxt(k,209)*y(k,139)
         mat(k,242) = rxt(k,227)*y(k,139)
         mat(k,512) = rxt(k,231)*y(k,138)
         mat(k,977) = mat(k,977) + rxt(k,146)*y(k,25) + rxt(k,205)*y(k,26) &
                      + rxt(k,208)*y(k,29) + rxt(k,149)*y(k,57)
         mat(k,882) = mat(k,882) + rxt(k,153)*y(k,131) + rxt(k,164)*y(k,139)
         mat(k,582) = rxt(k,233)*y(k,139)
         mat(k,139) = .500_r8*rxt(k,342)*y(k,139)
         mat(k,565) = rxt(k,256)*y(k,139)
         mat(k,334) = rxt(k,257)*y(k,139)
         mat(k,797) = mat(k,797) + rxt(k,109)*y(k,99)
         mat(k,271) = rxt(k,149)*y(k,38) + rxt(k,105)*y(k,98) + rxt(k,114)*y(k,139)
         mat(k,522) = rxt(k,271)*y(k,139)
         mat(k,442) = .370_r8*rxt(k,313)*y(k,102)
         mat(k,368) = mat(k,368) + .794_r8*rxt(k,306)*y(k,89) + .794_r8*rxt(k,307) &
                      *y(k,91)
         mat(k,576) = .140_r8*rxt(k,288)*y(k,102)
         mat(k,145) = .200_r8*rxt(k,290)*y(k,139)
         mat(k,281) = .500_r8*rxt(k,297)*y(k,139)
         mat(k,746) = .280_r8*rxt(k,298)*y(k,102)
         mat(k,856) = mat(k,856) + .794_r8*rxt(k,306)*y(k,71) + rxt(k,239)*y(k,128) &
                      + rxt(k,266)*y(k,129) + rxt(k,226)*y(k,131) + .250_r8*rxt(k,253) &
                      *y(k,133) + .920_r8*rxt(k,311)*y(k,135) + .470_r8*rxt(k,286) &
                      *y(k,136) + rxt(k,274)*y(k,140) + rxt(k,321)*y(k,142)
         mat(k,1143) = mat(k,1143) + rxt(k,220)*y(k,25) + .794_r8*rxt(k,307)*y(k,71) &
                      + rxt(k,316)*y(k,106) + rxt(k,312)*y(k,135) + .470_r8*rxt(k,285) &
                      *y(k,136) + rxt(k,133)*y(k,139) + rxt(k,322)*y(k,142)
         mat(k,1300) = mat(k,1300) + rxt(k,221)*y(k,25) + rxt(k,105)*y(k,57)
         mat(k,1349) = rxt(k,109)*y(k,54) + rxt(k,255)*y(k,132)
         mat(k,1080) = mat(k,1080) + .130_r8*rxt(k,235)*y(k,9) + .280_r8*rxt(k,263) &
                      *y(k,13) + .370_r8*rxt(k,313)*y(k,70) + .140_r8*rxt(k,288) &
                      *y(k,73) + .280_r8*rxt(k,298)*y(k,76) + rxt(k,117)*y(k,139)
         mat(k,478) = rxt(k,316)*y(k,91) + rxt(k,317)*y(k,139)
         mat(k,467) = rxt(k,332)*y(k,139)
         mat(k,391) = mat(k,391) + rxt(k,239)*y(k,89) + 2.400_r8*rxt(k,236)*y(k,128) &
                      + rxt(k,237)*y(k,131)
         mat(k,425) = mat(k,425) + rxt(k,266)*y(k,89) + rxt(k,264)*y(k,131)
         mat(k,706) = mat(k,706) + .900_r8*rxt(k,247)*y(k,131) + rxt(k,308)*y(k,135) &
                      + .470_r8*rxt(k,282)*y(k,136) + rxt(k,318)*y(k,142)
         mat(k,785) = mat(k,785) + rxt(k,153)*y(k,41) + rxt(k,226)*y(k,89) &
                      + rxt(k,237)*y(k,128) + rxt(k,264)*y(k,129) + .900_r8*rxt(k,247) &
                      *y(k,130) + 4.000_r8*rxt(k,223)*y(k,131) + rxt(k,309)*y(k,135) &
                      + .730_r8*rxt(k,283)*y(k,136) + rxt(k,292)*y(k,137) &
                      + .300_r8*rxt(k,276)*y(k,141) + .800_r8*rxt(k,319)*y(k,142)
         mat(k,209) = rxt(k,255)*y(k,99)
         mat(k,358) = mat(k,358) + .250_r8*rxt(k,253)*y(k,89)
         mat(k,637) = mat(k,637) + .920_r8*rxt(k,311)*y(k,89) + rxt(k,312)*y(k,91) &
                      + rxt(k,308)*y(k,130) + rxt(k,309)*y(k,131)
         mat(k,659) = mat(k,659) + .470_r8*rxt(k,286)*y(k,89) + .470_r8*rxt(k,285) &
                      *y(k,91) + .470_r8*rxt(k,282)*y(k,130) + .730_r8*rxt(k,283) &
                      *y(k,131)
         mat(k,678) = mat(k,678) + rxt(k,292)*y(k,131)
         mat(k,1182) = rxt(k,231)*y(k,36)
         mat(k,1268) = mat(k,1268) + rxt(k,194)*y(k,6) + rxt(k,240)*y(k,10) &
                      + rxt(k,206)*y(k,26) + rxt(k,209)*y(k,29) + rxt(k,227)*y(k,34) &
                      + rxt(k,164)*y(k,41) + rxt(k,233)*y(k,44) + .500_r8*rxt(k,342) &
                      *y(k,46) + rxt(k,256)*y(k,52) + rxt(k,257)*y(k,53) + rxt(k,114) &
                      *y(k,57) + rxt(k,271)*y(k,68) + .200_r8*rxt(k,290)*y(k,74) &
                      + .500_r8*rxt(k,297)*y(k,75) + rxt(k,133)*y(k,91) + rxt(k,117) &
                      *y(k,102) + rxt(k,317)*y(k,106) + rxt(k,332)*y(k,115)
         mat(k,378) = mat(k,378) + rxt(k,274)*y(k,89)
         mat(k,492) = mat(k,492) + .300_r8*rxt(k,276)*y(k,131)
         mat(k,600) = mat(k,600) + rxt(k,321)*y(k,89) + rxt(k,322)*y(k,91) &
                      + rxt(k,318)*y(k,130) + .800_r8*rxt(k,319)*y(k,131)
         mat(k,628) = -(rxt(k,308)*y(k,130) + rxt(k,309)*y(k,131) + rxt(k,310) &
                      *y(k,134) + rxt(k,311)*y(k,89) + rxt(k,312)*y(k,91))
         mat(k,697) = -rxt(k,308)*y(k,135)
         mat(k,772) = -rxt(k,309)*y(k,135)
         mat(k,1021) = -rxt(k,310)*y(k,135)
         mat(k,843) = -rxt(k,311)*y(k,135)
         mat(k,1129) = -rxt(k,312)*y(k,135)
         mat(k,435) = rxt(k,314)*y(k,139)
         mat(k,247) = .200_r8*rxt(k,315)*y(k,139)
         mat(k,1129) = mat(k,1129) + 1.700_r8*rxt(k,324)*y(k,126)
         mat(k,303) = 1.700_r8*rxt(k,324)*y(k,91) + 1.640_r8*rxt(k,326)*y(k,139)
         mat(k,1254) = rxt(k,314)*y(k,70) + .200_r8*rxt(k,315)*y(k,72) &
                      + 1.640_r8*rxt(k,326)*y(k,126)
         mat(k,652) = -(rxt(k,282)*y(k,130) + rxt(k,283)*y(k,131) + rxt(k,284) &
                      *y(k,134) + rxt(k,285)*y(k,91) + (rxt(k,286) + rxt(k,287) &
                      ) * y(k,89))
         mat(k,698) = -rxt(k,282)*y(k,136)
         mat(k,773) = -rxt(k,283)*y(k,136)
         mat(k,1022) = -rxt(k,284)*y(k,136)
         mat(k,1130) = -rxt(k,285)*y(k,136)
         mat(k,844) = -(rxt(k,286) + rxt(k,287)) * y(k,136)
         mat(k,572) = .500_r8*rxt(k,289)*y(k,139)
         mat(k,143) = .200_r8*rxt(k,290)*y(k,139)
         mat(k,735) = rxt(k,299)*y(k,139)
         mat(k,1255) = .500_r8*rxt(k,289)*y(k,73) + .200_r8*rxt(k,290)*y(k,74) &
                      + rxt(k,299)*y(k,76)
         mat(k,672) = -(rxt(k,291)*y(k,130) + rxt(k,292)*y(k,131) + rxt(k,293) &
                      *y(k,134) + 4._r8*rxt(k,294)*y(k,137) + rxt(k,295)*y(k,89) &
                      + rxt(k,296)*y(k,91) + rxt(k,300)*y(k,90))
         mat(k,699) = -rxt(k,291)*y(k,137)
         mat(k,774) = -rxt(k,292)*y(k,137)
         mat(k,1023) = -rxt(k,293)*y(k,137)
         mat(k,845) = -rxt(k,295)*y(k,137)
         mat(k,1131) = -rxt(k,296)*y(k,137)
         mat(k,907) = -rxt(k,300)*y(k,137)
         mat(k,573) = .500_r8*rxt(k,289)*y(k,139)
         mat(k,144) = .500_r8*rxt(k,290)*y(k,139)
         mat(k,1256) = .500_r8*rxt(k,289)*y(k,73) + .500_r8*rxt(k,290)*y(k,74)
         mat(k,1186) = -(rxt(k,85)*y(k,55) + rxt(k,86)*y(k,143) + (rxt(k,88) + rxt(k,89) &
                      + rxt(k,90)) * y(k,99) + (rxt(k,91) + rxt(k,92)) * y(k,102) &
                      + (rxt(k,138) + rxt(k,139)) * y(k,78) + rxt(k,171)*y(k,16) &
                      + rxt(k,172)*y(k,17) + rxt(k,173)*y(k,19) + rxt(k,174)*y(k,20) &
                      + rxt(k,175)*y(k,21) + rxt(k,176)*y(k,22) + rxt(k,177)*y(k,23) &
                      + (rxt(k,178) + rxt(k,179)) * y(k,63) + rxt(k,198)*y(k,18) &
                      + rxt(k,199)*y(k,37) + rxt(k,200)*y(k,56) + (rxt(k,201) &
                      + rxt(k,202)) * y(k,59) + rxt(k,215)*y(k,24) + rxt(k,216) &
                      *y(k,26) + rxt(k,217)*y(k,60) + rxt(k,218)*y(k,61) + rxt(k,219) &
                      *y(k,62) + (rxt(k,230) + rxt(k,231) + rxt(k,232)) * y(k,36))
         mat(k,499) = -rxt(k,85)*y(k,138)
         mat(k,1374) = -rxt(k,86)*y(k,138)
         mat(k,1353) = -(rxt(k,88) + rxt(k,89) + rxt(k,90)) * y(k,138)
         mat(k,1084) = -(rxt(k,91) + rxt(k,92)) * y(k,138)
         mat(k,113) = -(rxt(k,138) + rxt(k,139)) * y(k,138)
         mat(k,48) = -rxt(k,171)*y(k,138)
         mat(k,74) = -rxt(k,172)*y(k,138)
         mat(k,54) = -rxt(k,173)*y(k,138)
         mat(k,57) = -rxt(k,174)*y(k,138)
         mat(k,60) = -rxt(k,175)*y(k,138)
         mat(k,63) = -rxt(k,176)*y(k,138)
         mat(k,66) = -rxt(k,177)*y(k,138)
         mat(k,944) = -(rxt(k,178) + rxt(k,179)) * y(k,138)
         mat(k,51) = -rxt(k,198)*y(k,138)
         mat(k,192) = -rxt(k,199)*y(k,138)
         mat(k,42) = -rxt(k,200)*y(k,138)
         mat(k,405) = -(rxt(k,201) + rxt(k,202)) * y(k,138)
         mat(k,229) = -rxt(k,215)*y(k,138)
         mat(k,297) = -rxt(k,216)*y(k,138)
         mat(k,94) = -rxt(k,217)*y(k,138)
         mat(k,98) = -rxt(k,218)*y(k,138)
         mat(k,105) = -rxt(k,219)*y(k,138)
         mat(k,513) = -(rxt(k,230) + rxt(k,231) + rxt(k,232)) * y(k,138)
         mat(k,1273) = -(rxt(k,113)*y(k,55) + rxt(k,114)*y(k,57) + rxt(k,115)*y(k,134) &
                      + rxt(k,116)*y(k,98) + rxt(k,117)*y(k,102) + (4._r8*rxt(k,118) &
                      + 4._r8*rxt(k,119)) * y(k,139) + rxt(k,121)*y(k,65) + rxt(k,133) &
                      *y(k,91) + rxt(k,134)*y(k,77) + rxt(k,142)*y(k,90) + rxt(k,143) &
                      *y(k,64) + rxt(k,162)*y(k,42) + (rxt(k,164) + rxt(k,165) &
                      ) * y(k,41) + rxt(k,167)*y(k,63) + rxt(k,170)*y(k,67) + rxt(k,194) &
                      *y(k,6) + rxt(k,196)*y(k,59) + rxt(k,204)*y(k,24) + rxt(k,206) &
                      *y(k,26) + rxt(k,207)*y(k,27) + rxt(k,209)*y(k,29) + rxt(k,211) &
                      *y(k,37) + rxt(k,212)*y(k,60) + rxt(k,213)*y(k,61) + rxt(k,214) &
                      *y(k,62) + rxt(k,222)*y(k,25) + rxt(k,227)*y(k,34) + rxt(k,228) &
                      *y(k,35) + rxt(k,229)*y(k,36) + rxt(k,233)*y(k,44) + rxt(k,240) &
                      *y(k,10) + rxt(k,241)*y(k,11) + rxt(k,243)*y(k,12) + rxt(k,245) &
                      *y(k,28) + rxt(k,250)*y(k,32) + rxt(k,251)*y(k,33) + rxt(k,256) &
                      *y(k,52) + rxt(k,257)*y(k,53) + rxt(k,258)*y(k,107) + rxt(k,259) &
                      *y(k,9) + rxt(k,267)*y(k,14) + rxt(k,268)*y(k,15) + rxt(k,270) &
                      *y(k,31) + rxt(k,271)*y(k,68) + rxt(k,272)*y(k,92) + rxt(k,275) &
                      *y(k,110) + rxt(k,279)*y(k,111) + rxt(k,280)*y(k,13) + rxt(k,281) &
                      *y(k,30) + rxt(k,289)*y(k,73) + rxt(k,290)*y(k,74) + rxt(k,297) &
                      *y(k,75) + rxt(k,299)*y(k,76) + rxt(k,302)*y(k,3) + rxt(k,303) &
                      *y(k,69) + rxt(k,314)*y(k,70) + rxt(k,315)*y(k,72) + rxt(k,317) &
                      *y(k,106) + rxt(k,323)*y(k,127) + rxt(k,326)*y(k,126) + (rxt(k,328) &
                      + rxt(k,342)) * y(k,46) + rxt(k,330)*y(k,105) + rxt(k,332) &
                      *y(k,115) + rxt(k,336)*y(k,112) + rxt(k,341)*y(k,114) + rxt(k,344) &
                      *y(k,84))
         mat(k,500) = -rxt(k,113)*y(k,139)
         mat(k,272) = -rxt(k,114)*y(k,139)
         mat(k,1039) = -rxt(k,115)*y(k,139)
         mat(k,1305) = -rxt(k,116)*y(k,139)
         mat(k,1085) = -rxt(k,117)*y(k,139)
         mat(k,222) = -rxt(k,121)*y(k,139)
         mat(k,1148) = -rxt(k,133)*y(k,139)
         mat(k,316) = -rxt(k,134)*y(k,139)
         mat(k,924) = -rxt(k,142)*y(k,139)
         mat(k,413) = -rxt(k,143)*y(k,139)
         mat(k,455) = -rxt(k,162)*y(k,139)
         mat(k,887) = -(rxt(k,164) + rxt(k,165)) * y(k,139)
         mat(k,945) = -rxt(k,167)*y(k,139)
         mat(k,399) = -rxt(k,170)*y(k,139)
         mat(k,1328) = -rxt(k,194)*y(k,139)
         mat(k,406) = -rxt(k,196)*y(k,139)
         mat(k,230) = -rxt(k,204)*y(k,139)
         mat(k,298) = -rxt(k,206)*y(k,139)
         mat(k,77) = -rxt(k,207)*y(k,139)
         mat(k,174) = -rxt(k,209)*y(k,139)
         mat(k,193) = -rxt(k,211)*y(k,139)
         mat(k,95) = -rxt(k,212)*y(k,139)
         mat(k,99) = -rxt(k,213)*y(k,139)
         mat(k,106) = -rxt(k,214)*y(k,139)
         mat(k,818) = -rxt(k,222)*y(k,139)
         mat(k,243) = -rxt(k,227)*y(k,139)
         mat(k,198) = -rxt(k,228)*y(k,139)
         mat(k,514) = -rxt(k,229)*y(k,139)
         mat(k,583) = -rxt(k,233)*y(k,139)
         mat(k,118) = -rxt(k,240)*y(k,139)
         mat(k,151) = -rxt(k,241)*y(k,139)
         mat(k,134) = -rxt(k,243)*y(k,139)
         mat(k,532) = -rxt(k,245)*y(k,139)
         mat(k,322) = -rxt(k,250)*y(k,139)
         mat(k,289) = -rxt(k,251)*y(k,139)
         mat(k,566) = -rxt(k,256)*y(k,139)
         mat(k,335) = -rxt(k,257)*y(k,139)
         mat(k,259) = -rxt(k,258)*y(k,139)
         mat(k,217) = -rxt(k,259)*y(k,139)
         mat(k,180) = -rxt(k,267)*y(k,139)
         mat(k,69) = -rxt(k,268)*y(k,139)
         mat(k,612) = -rxt(k,270)*y(k,139)
         mat(k,523) = -rxt(k,271)*y(k,139)
         mat(k,187) = -rxt(k,272)*y(k,139)
         mat(k,266) = -rxt(k,275)*y(k,139)
         mat(k,204) = -rxt(k,279)*y(k,139)
         mat(k,556) = -rxt(k,280)*y(k,139)
         mat(k,349) = -rxt(k,281)*y(k,139)
         mat(k,578) = -rxt(k,289)*y(k,139)
         mat(k,146) = -rxt(k,290)*y(k,139)
         mat(k,283) = -rxt(k,297)*y(k,139)
         mat(k,750) = -rxt(k,299)*y(k,139)
         mat(k,39) = -rxt(k,302)*y(k,139)
         mat(k,159) = -rxt(k,303)*y(k,139)
         mat(k,445) = -rxt(k,314)*y(k,139)
         mat(k,251) = -rxt(k,315)*y(k,139)
         mat(k,480) = -rxt(k,317)*y(k,139)
         mat(k,102) = -rxt(k,323)*y(k,139)
         mat(k,308) = -rxt(k,326)*y(k,139)
         mat(k,141) = -(rxt(k,328) + rxt(k,342)) * y(k,139)
         mat(k,166) = -rxt(k,330)*y(k,139)
         mat(k,468) = -rxt(k,332)*y(k,139)
         mat(k,328) = -rxt(k,336)*y(k,139)
         mat(k,725) = -rxt(k,341)*y(k,139)
         mat(k,44) = -rxt(k,344)*y(k,139)
         mat(k,217) = mat(k,217) + .130_r8*rxt(k,235)*y(k,102)
         mat(k,151) = mat(k,151) + .500_r8*rxt(k,241)*y(k,139)
         mat(k,556) = mat(k,556) + .360_r8*rxt(k,263)*y(k,102)
         mat(k,818) = mat(k,818) + rxt(k,221)*y(k,98)
         mat(k,198) = mat(k,198) + .300_r8*rxt(k,228)*y(k,139)
         mat(k,514) = mat(k,514) + rxt(k,230)*y(k,138)
         mat(k,982) = rxt(k,151)*y(k,134)
         mat(k,800) = rxt(k,112)*y(k,102) + 2.000_r8*rxt(k,107)*y(k,134)
         mat(k,500) = mat(k,500) + rxt(k,104)*y(k,98) + rxt(k,85)*y(k,138)
         mat(k,272) = mat(k,272) + rxt(k,105)*y(k,98)
         mat(k,406) = mat(k,406) + rxt(k,195)*y(k,98) + rxt(k,201)*y(k,138)
         mat(k,945) = mat(k,945) + rxt(k,166)*y(k,98) + rxt(k,178)*y(k,138)
         mat(k,341) = rxt(k,197)*y(k,98)
         mat(k,399) = mat(k,399) + rxt(k,169)*y(k,98)
         mat(k,445) = mat(k,445) + .320_r8*rxt(k,313)*y(k,102)
         mat(k,370) = .206_r8*rxt(k,305)*y(k,134)
         mat(k,578) = mat(k,578) + .240_r8*rxt(k,288)*y(k,102)
         mat(k,146) = mat(k,146) + .100_r8*rxt(k,290)*y(k,139)
         mat(k,750) = mat(k,750) + .360_r8*rxt(k,298)*y(k,102)
         mat(k,861) = rxt(k,135)*y(k,134)
         mat(k,1148) = mat(k,1148) + rxt(k,130)*y(k,134)
         mat(k,1305) = mat(k,1305) + rxt(k,221)*y(k,25) + rxt(k,104)*y(k,55) &
                      + rxt(k,105)*y(k,57) + rxt(k,195)*y(k,59) + rxt(k,166)*y(k,63) &
                      + rxt(k,197)*y(k,66) + rxt(k,169)*y(k,67) + rxt(k,110)*y(k,134)
         mat(k,1085) = mat(k,1085) + .130_r8*rxt(k,235)*y(k,9) + .360_r8*rxt(k,263) &
                      *y(k,13) + rxt(k,112)*y(k,54) + .320_r8*rxt(k,313)*y(k,70) &
                      + .240_r8*rxt(k,288)*y(k,73) + .360_r8*rxt(k,298)*y(k,76) &
                      + 1.156_r8*rxt(k,325)*y(k,126) + rxt(k,111)*y(k,134)
         mat(k,266) = mat(k,266) + .500_r8*rxt(k,275)*y(k,139)
         mat(k,308) = mat(k,308) + 1.156_r8*rxt(k,325)*y(k,102)
         mat(k,102) = mat(k,102) + .500_r8*rxt(k,323)*y(k,139)
         mat(k,709) = .450_r8*rxt(k,248)*y(k,134)
         mat(k,1039) = mat(k,1039) + rxt(k,151)*y(k,38) + 2.000_r8*rxt(k,107)*y(k,54) &
                      + .206_r8*rxt(k,305)*y(k,71) + rxt(k,135)*y(k,89) + rxt(k,130) &
                      *y(k,91) + rxt(k,110)*y(k,98) + rxt(k,111)*y(k,102) &
                      + .450_r8*rxt(k,248)*y(k,130) + .450_r8*rxt(k,293)*y(k,137) &
                      + .150_r8*rxt(k,277)*y(k,141)
         mat(k,681) = .450_r8*rxt(k,293)*y(k,134)
         mat(k,1187) = rxt(k,230)*y(k,36) + rxt(k,85)*y(k,55) + rxt(k,201)*y(k,59) &
                      + rxt(k,178)*y(k,63) + 2.000_r8*rxt(k,86)*y(k,143)
         mat(k,1273) = mat(k,1273) + .500_r8*rxt(k,241)*y(k,11) + .300_r8*rxt(k,228) &
                      *y(k,35) + .100_r8*rxt(k,290)*y(k,74) + .500_r8*rxt(k,275) &
                      *y(k,110) + .500_r8*rxt(k,323)*y(k,127)
         mat(k,493) = .150_r8*rxt(k,277)*y(k,134)
         mat(k,1375) = 2.000_r8*rxt(k,86)*y(k,138)
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
         mat(k,372) = -(rxt(k,273)*y(k,134) + rxt(k,274)*y(k,89))
         mat(k,1005) = -rxt(k,273)*y(k,140)
         mat(k,830) = -rxt(k,274)*y(k,140)
         mat(k,536) = rxt(k,280)*y(k,139)
         mat(k,261) = .500_r8*rxt(k,275)*y(k,139)
         mat(k,1232) = rxt(k,280)*y(k,13) + .500_r8*rxt(k,275)*y(k,110)
         mat(k,484) = -(rxt(k,276)*y(k,131) + rxt(k,277)*y(k,134) + rxt(k,278)*y(k,89))
         mat(k,764) = -rxt(k,276)*y(k,141)
         mat(k,1012) = -rxt(k,277)*y(k,141)
         mat(k,835) = -rxt(k,278)*y(k,141)
         mat(k,346) = rxt(k,281)*y(k,139)
         mat(k,201) = rxt(k,279)*y(k,139)
         mat(k,1243) = rxt(k,281)*y(k,30) + rxt(k,279)*y(k,111)
         mat(k,593) = -(rxt(k,318)*y(k,130) + rxt(k,319)*y(k,131) + rxt(k,320) &
                      *y(k,134) + rxt(k,321)*y(k,89) + rxt(k,322)*y(k,91))
         mat(k,695) = -rxt(k,318)*y(k,142)
         mat(k,770) = -rxt(k,319)*y(k,142)
         mat(k,1019) = -rxt(k,320)*y(k,142)
         mat(k,841) = -rxt(k,321)*y(k,142)
         mat(k,1127) = -rxt(k,322)*y(k,142)
         mat(k,158) = rxt(k,303)*y(k,139)
         mat(k,246) = .800_r8*rxt(k,315)*y(k,139)
         mat(k,101) = .500_r8*rxt(k,323)*y(k,139)
         mat(k,1252) = rxt(k,303)*y(k,69) + .800_r8*rxt(k,315)*y(k,72) &
                      + .500_r8*rxt(k,323)*y(k,127)
         mat(k,1379) = -(rxt(k,86)*y(k,138) + rxt(k,343)*y(k,116))
         mat(k,1191) = -rxt(k,86)*y(k,143)
         mat(k,123) = -rxt(k,343)*y(k,143)
         mat(k,135) = rxt(k,243)*y(k,139)
         mat(k,181) = rxt(k,267)*y(k,139)
         mat(k,70) = rxt(k,268)*y(k,139)
         mat(k,231) = rxt(k,204)*y(k,139)
         mat(k,822) = rxt(k,222)*y(k,139)
         mat(k,299) = rxt(k,206)*y(k,139)
         mat(k,78) = rxt(k,207)*y(k,139)
         mat(k,533) = rxt(k,245)*y(k,139)
         mat(k,175) = rxt(k,209)*y(k,139)
         mat(k,350) = rxt(k,281)*y(k,139)
         mat(k,615) = rxt(k,270)*y(k,139)
         mat(k,323) = rxt(k,250)*y(k,139)
         mat(k,290) = rxt(k,251)*y(k,139)
         mat(k,199) = rxt(k,228)*y(k,139)
         mat(k,517) = rxt(k,229)*y(k,139)
         mat(k,803) = rxt(k,108)*y(k,134)
         mat(k,502) = rxt(k,113)*y(k,139)
         mat(k,274) = rxt(k,114)*y(k,139)
         mat(k,409) = rxt(k,196)*y(k,139)
         mat(k,107) = rxt(k,214)*y(k,139)
         mat(k,949) = (rxt(k,360)+rxt(k,365))*y(k,66) + (rxt(k,353)+rxt(k,359) &
                       +rxt(k,364))*y(k,67) + rxt(k,167)*y(k,139)
         mat(k,414) = rxt(k,143)*y(k,139)
         mat(k,224) = rxt(k,121)*y(k,139)
         mat(k,344) = (rxt(k,360)+rxt(k,365))*y(k,63)
         mat(k,401) = (rxt(k,353)+rxt(k,359)+rxt(k,364))*y(k,63) + rxt(k,170)*y(k,139)
         mat(k,579) = .500_r8*rxt(k,289)*y(k,139)
         mat(k,45) = rxt(k,344)*y(k,139)
         mat(k,267) = rxt(k,275)*y(k,139)
         mat(k,205) = rxt(k,279)*y(k,139)
         mat(k,1043) = rxt(k,108)*y(k,54) + rxt(k,115)*y(k,139)
         mat(k,1277) = rxt(k,243)*y(k,12) + rxt(k,267)*y(k,14) + rxt(k,268)*y(k,15) &
                      + rxt(k,204)*y(k,24) + rxt(k,222)*y(k,25) + rxt(k,206)*y(k,26) &
                      + rxt(k,207)*y(k,27) + rxt(k,245)*y(k,28) + rxt(k,209)*y(k,29) &
                      + rxt(k,281)*y(k,30) + rxt(k,270)*y(k,31) + rxt(k,250)*y(k,32) &
                      + rxt(k,251)*y(k,33) + rxt(k,228)*y(k,35) + rxt(k,229)*y(k,36) &
                      + rxt(k,113)*y(k,55) + rxt(k,114)*y(k,57) + rxt(k,196)*y(k,59) &
                      + rxt(k,214)*y(k,62) + rxt(k,167)*y(k,63) + rxt(k,143)*y(k,64) &
                      + rxt(k,121)*y(k,65) + rxt(k,170)*y(k,67) + .500_r8*rxt(k,289) &
                      *y(k,73) + rxt(k,344)*y(k,84) + rxt(k,275)*y(k,110) + rxt(k,279) &
                      *y(k,111) + rxt(k,115)*y(k,134) + 2.000_r8*rxt(k,118)*y(k,139)
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
         mat(k, 37) = mat(k, 37) + lmat(k, 37)
         mat(k, 40) = mat(k, 40) + lmat(k, 40)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 47) = mat(k, 47) + lmat(k, 47)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 52) = mat(k, 52) + lmat(k, 52)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 73) = mat(k, 73) + lmat(k, 73)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 79) = lmat(k, 79)
         mat(k, 80) = lmat(k, 80)
         mat(k, 81) = lmat(k, 81)
         mat(k, 82) = lmat(k, 82)
         mat(k, 83) = lmat(k, 83)
         mat(k, 84) = lmat(k, 84)
         mat(k, 85) = lmat(k, 85)
         mat(k, 86) = lmat(k, 86)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 92) = mat(k, 92) + lmat(k, 92)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 108) = lmat(k, 108)
         mat(k, 109) = lmat(k, 109)
         mat(k, 110) = lmat(k, 110)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 121) = lmat(k, 121)
         mat(k, 122) = lmat(k, 122)
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
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 170) = lmat(k, 170)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 196) = lmat(k, 196)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 198) = mat(k, 198) + lmat(k, 198)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 218) = mat(k, 218) + lmat(k, 218)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 220) = lmat(k, 220)
         mat(k, 221) = lmat(k, 221)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 248) = lmat(k, 248)
         mat(k, 249) = lmat(k, 249)
         mat(k, 250) = lmat(k, 250)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 253) = lmat(k, 253)
         mat(k, 254) = lmat(k, 254)
         mat(k, 255) = lmat(k, 255)
         mat(k, 257) = lmat(k, 257)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 263) = lmat(k, 263)
         mat(k, 264) = lmat(k, 264)
         mat(k, 265) = lmat(k, 265)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 278) = lmat(k, 278)
         mat(k, 280) = lmat(k, 280)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 287) = lmat(k, 287)
         mat(k, 289) = mat(k, 289) + lmat(k, 289)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 292) = lmat(k, 292)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 337) = mat(k, 337) + lmat(k, 337)
         mat(k, 340) = lmat(k, 340)
         mat(k, 341) = mat(k, 341) + lmat(k, 341)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 347) = lmat(k, 347)
         mat(k, 348) = lmat(k, 348)
         mat(k, 353) = mat(k, 353) + lmat(k, 353)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 372) = mat(k, 372) + lmat(k, 372)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 402) = mat(k, 402) + lmat(k, 402)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 411) = lmat(k, 411)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 451) = lmat(k, 451)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 460) = lmat(k, 460)
         mat(k, 461) = lmat(k, 461)
         mat(k, 462) = lmat(k, 462)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 466) = lmat(k, 466)
         mat(k, 469) = lmat(k, 469)
         mat(k, 472) = lmat(k, 472)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 474) = lmat(k, 474)
         mat(k, 476) = lmat(k, 476)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 478) = mat(k, 478) + lmat(k, 478)
         mat(k, 484) = mat(k, 484) + lmat(k, 484)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 503) = lmat(k, 503)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 515) = lmat(k, 515)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 520) = lmat(k, 520)
         mat(k, 521) = lmat(k, 521)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 525) = mat(k, 525) + lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 528) = lmat(k, 528)
         mat(k, 530) = lmat(k, 530)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 562) = mat(k, 562) + lmat(k, 562)
         mat(k, 563) = lmat(k, 563)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 565) = mat(k, 565) + lmat(k, 565)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 570) = mat(k, 570) + lmat(k, 570)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 593) = mat(k, 593) + lmat(k, 593)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 608) = mat(k, 608) + lmat(k, 608)
         mat(k, 610) = lmat(k, 610)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 652) = mat(k, 652) + lmat(k, 652)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 700) = mat(k, 700) + lmat(k, 700)
         mat(k, 714) = lmat(k, 714)
         mat(k, 716) = mat(k, 716) + lmat(k, 716)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 732) = lmat(k, 732)
         mat(k, 733) = mat(k, 733) + lmat(k, 733)
         mat(k, 737) = mat(k, 737) + lmat(k, 737)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 739) = lmat(k, 739)
         mat(k, 777) = mat(k, 777) + lmat(k, 777)
         mat(k, 794) = mat(k, 794) + lmat(k, 794)
         mat(k, 806) = lmat(k, 806)
         mat(k, 807) = mat(k, 807) + lmat(k, 807)
         mat(k, 808) = mat(k, 808) + lmat(k, 808)
         mat(k, 809) = mat(k, 809) + lmat(k, 809)
         mat(k, 825) = mat(k, 825) + lmat(k, 825)
         mat(k, 851) = mat(k, 851) + lmat(k, 851)
         mat(k, 862) = mat(k, 862) + lmat(k, 862)
         mat(k, 878) = mat(k, 878) + lmat(k, 878)
         mat(k, 881) = mat(k, 881) + lmat(k, 881)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 900) = mat(k, 900) + lmat(k, 900)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 916) = mat(k, 916) + lmat(k, 916)
         mat(k, 924) = mat(k, 924) + lmat(k, 924)
         mat(k, 925) = mat(k, 925) + lmat(k, 925)
         mat(k, 935) = mat(k, 935) + lmat(k, 935)
         mat(k, 938) = mat(k, 938) + lmat(k, 938)
         mat(k, 939) = mat(k, 939) + lmat(k, 939)
         mat(k, 976) = mat(k, 976) + lmat(k, 976)
         mat(k,1034) = mat(k,1034) + lmat(k,1034)
         mat(k,1043) = mat(k,1043) + lmat(k,1043)
         mat(k,1044) = mat(k,1044) + lmat(k,1044)
         mat(k,1081) = mat(k,1081) + lmat(k,1081)
         mat(k,1084) = mat(k,1084) + lmat(k,1084)
         mat(k,1086) = mat(k,1086) + lmat(k,1086)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1138) = mat(k,1138) + lmat(k,1138)
         mat(k,1140) = mat(k,1140) + lmat(k,1140)
         mat(k,1146) = mat(k,1146) + lmat(k,1146)
         mat(k,1149) = mat(k,1149) + lmat(k,1149)
         mat(k,1151) = mat(k,1151) + lmat(k,1151)
         mat(k,1186) = mat(k,1186) + lmat(k,1186)
         mat(k,1188) = mat(k,1188) + lmat(k,1188)
         mat(k,1273) = mat(k,1273) + lmat(k,1273)
         mat(k,1306) = mat(k,1306) + lmat(k,1306)
         mat(k,1325) = mat(k,1325) + lmat(k,1325)
         mat(k,1329) = mat(k,1329) + lmat(k,1329)
         mat(k,1330) = mat(k,1330) + lmat(k,1330)
         mat(k,1353) = mat(k,1353) + lmat(k,1353)
         mat(k,1355) = mat(k,1355) + lmat(k,1355)
         mat(k,1357) = mat(k,1357) + lmat(k,1357)
         mat(k,1362) = lmat(k,1362)
         mat(k,1364) = lmat(k,1364)
         mat(k,1374) = mat(k,1374) + lmat(k,1374)
         mat(k,1375) = mat(k,1375) + lmat(k,1375)
         mat(k,1376) = lmat(k,1376)
         mat(k,1379) = mat(k,1379) + lmat(k,1379)
         mat(k, 91) = 0._r8
         mat(k, 315) = 0._r8
         mat(k, 339) = 0._r8
         mat(k, 354) = 0._r8
         mat(k, 359) = 0._r8
         mat(k, 360) = 0._r8
         mat(k, 373) = 0._r8
         mat(k, 379) = 0._r8
         mat(k, 381) = 0._r8
         mat(k, 392) = 0._r8
         mat(k, 418) = 0._r8
         mat(k, 420) = 0._r8
         mat(k, 426) = 0._r8
         mat(k, 428) = 0._r8
         mat(k, 431) = 0._r8
         mat(k, 440) = 0._r8
         mat(k, 441) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 463) = 0._r8
         mat(k, 470) = 0._r8
         mat(k, 475) = 0._r8
         mat(k, 481) = 0._r8
         mat(k, 494) = 0._r8
         mat(k, 516) = 0._r8
         mat(k, 529) = 0._r8
         mat(k, 539) = 0._r8
         mat(k, 543) = 0._r8
         mat(k, 544) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 548) = 0._r8
         mat(k, 549) = 0._r8
         mat(k, 550) = 0._r8
         mat(k, 551) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 557) = 0._r8
         mat(k, 558) = 0._r8
         mat(k, 559) = 0._r8
         mat(k, 567) = 0._r8
         mat(k, 568) = 0._r8
         mat(k, 584) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 602) = 0._r8
         mat(k, 603) = 0._r8
         mat(k, 604) = 0._r8
         mat(k, 609) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 614) = 0._r8
         mat(k, 620) = 0._r8
         mat(k, 625) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 630) = 0._r8
         mat(k, 638) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 641) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 643) = 0._r8
         mat(k, 650) = 0._r8
         mat(k, 653) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 662) = 0._r8
         mat(k, 663) = 0._r8
         mat(k, 664) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 670) = 0._r8
         mat(k, 671) = 0._r8
         mat(k, 682) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 684) = 0._r8
         mat(k, 708) = 0._r8
         mat(k, 710) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 712) = 0._r8
         mat(k, 722) = 0._r8
         mat(k, 729) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 744) = 0._r8
         mat(k, 745) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 778) = 0._r8
         mat(k, 783) = 0._r8
         mat(k, 786) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 788) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 799) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 814) = 0._r8
         mat(k, 817) = 0._r8
         mat(k, 820) = 0._r8
         mat(k, 821) = 0._r8
         mat(k, 849) = 0._r8
         mat(k, 854) = 0._r8
         mat(k, 860) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 875) = 0._r8
         mat(k, 883) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 886) = 0._r8
         mat(k, 891) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 912) = 0._r8
         mat(k, 913) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 921) = 0._r8
         mat(k, 923) = 0._r8
         mat(k, 928) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 940) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 943) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 948) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 965) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 972) = 0._r8
         mat(k, 974) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 983) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 986) = 0._r8
         mat(k, 994) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1016) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1020) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1055) = 0._r8
         mat(k,1056) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1061) = 0._r8
         mat(k,1064) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1092) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1096) = 0._r8
         mat(k,1100) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1120) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1141) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1317) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1327) = 0._r8
         mat(k,1332) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1348) = 0._r8
         mat(k,1351) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1363) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1367) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1370) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1372) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1378) = 0._r8
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
         mat(k, 32) = mat(k, 32) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 34) = mat(k, 34) - dti(k)
         mat(k, 37) = mat(k, 37) - dti(k)
         mat(k, 40) = mat(k, 40) - dti(k)
         mat(k, 43) = mat(k, 43) - dti(k)
         mat(k, 46) = mat(k, 46) - dti(k)
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 58) = mat(k, 58) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 79) = mat(k, 79) - dti(k)
         mat(k, 82) = mat(k, 82) - dti(k)
         mat(k, 85) = mat(k, 85) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
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
         mat(k, 211) = mat(k, 211) - dti(k)
         mat(k, 218) = mat(k, 218) - dti(k)
         mat(k, 225) = mat(k, 225) - dti(k)
         mat(k, 232) = mat(k, 232) - dti(k)
         mat(k, 240) = mat(k, 240) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 260) = mat(k, 260) - dti(k)
         mat(k, 268) = mat(k, 268) - dti(k)
         mat(k, 275) = mat(k, 275) - dti(k)
         mat(k, 284) = mat(k, 284) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 319) = mat(k, 319) - dti(k)
         mat(k, 324) = mat(k, 324) - dti(k)
         mat(k, 331) = mat(k, 331) - dti(k)
         mat(k, 337) = mat(k, 337) - dti(k)
         mat(k, 345) = mat(k, 345) - dti(k)
         mat(k, 353) = mat(k, 353) - dti(k)
         mat(k, 361) = mat(k, 361) - dti(k)
         mat(k, 372) = mat(k, 372) - dti(k)
         mat(k, 385) = mat(k, 385) - dti(k)
         mat(k, 395) = mat(k, 395) - dti(k)
         mat(k, 402) = mat(k, 402) - dti(k)
         mat(k, 410) = mat(k, 410) - dti(k)
         mat(k, 417) = mat(k, 417) - dti(k)
         mat(k, 430) = mat(k, 430) - dti(k)
         mat(k, 449) = mat(k, 449) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 465) = mat(k, 465) - dti(k)
         mat(k, 473) = mat(k, 473) - dti(k)
         mat(k, 484) = mat(k, 484) - dti(k)
         mat(k, 495) = mat(k, 495) - dti(k)
         mat(k, 505) = mat(k, 505) - dti(k)
         mat(k, 518) = mat(k, 518) - dti(k)
         mat(k, 525) = mat(k, 525) - dti(k)
         mat(k, 541) = mat(k, 541) - dti(k)
         mat(k, 562) = mat(k, 562) - dti(k)
         mat(k, 569) = mat(k, 569) - dti(k)
         mat(k, 581) = mat(k, 581) - dti(k)
         mat(k, 593) = mat(k, 593) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 628) = mat(k, 628) - dti(k)
         mat(k, 652) = mat(k, 652) - dti(k)
         mat(k, 672) = mat(k, 672) - dti(k)
         mat(k, 700) = mat(k, 700) - dti(k)
         mat(k, 716) = mat(k, 716) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 777) = mat(k, 777) - dti(k)
         mat(k, 794) = mat(k, 794) - dti(k)
         mat(k, 809) = mat(k, 809) - dti(k)
         mat(k, 851) = mat(k, 851) - dti(k)
         mat(k, 878) = mat(k, 878) - dti(k)
         mat(k, 916) = mat(k, 916) - dti(k)
         mat(k, 938) = mat(k, 938) - dti(k)
         mat(k, 976) = mat(k, 976) - dti(k)
         mat(k,1034) = mat(k,1034) - dti(k)
         mat(k,1081) = mat(k,1081) - dti(k)
         mat(k,1099) = mat(k,1099) - dti(k)
         mat(k,1146) = mat(k,1146) - dti(k)
         mat(k,1186) = mat(k,1186) - dti(k)
         mat(k,1273) = mat(k,1273) - dti(k)
         mat(k,1306) = mat(k,1306) - dti(k)
         mat(k,1330) = mat(k,1330) - dti(k)
         mat(k,1357) = mat(k,1357) - dti(k)
         mat(k,1379) = mat(k,1379) - dti(k)
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
