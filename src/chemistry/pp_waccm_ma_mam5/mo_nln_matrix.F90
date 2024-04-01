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
         mat(k,808) = -(rxt(k,191)*y(k,17) + rxt(k,192)*y(k,90) + rxt(k,193)*y(k,71))
         mat(k,486) = -rxt(k,191)*y(k,3)
         mat(k,563) = -rxt(k,192)*y(k,3)
         mat(k,538) = -rxt(k,193)*y(k,3)
         mat(k,927) = 4.000_r8*rxt(k,194)*y(k,5) + (rxt(k,195)+rxt(k,196))*y(k,28) &
                      + rxt(k,199)*y(k,61) + rxt(k,202)*y(k,69) + rxt(k,253)*y(k,79) &
                      + rxt(k,203)*y(k,99)
         mat(k,59) = rxt(k,181)*y(k,95)
         mat(k,64) = rxt(k,207)*y(k,95)
         mat(k,178) = 2.000_r8*rxt(k,218)*y(k,25) + 2.000_r8*rxt(k,230)*y(k,95) &
                      + 2.000_r8*rxt(k,219)*y(k,99)
         mat(k,216) = rxt(k,220)*y(k,25) + rxt(k,231)*y(k,95) + rxt(k,221)*y(k,99)
         mat(k,165) = 3.000_r8*rxt(k,225)*y(k,25) + 3.000_r8*rxt(k,208)*y(k,95) &
                      + 3.000_r8*rxt(k,226)*y(k,99)
         mat(k,786) = 2.000_r8*rxt(k,218)*y(k,16) + rxt(k,220)*y(k,18) &
                      + 3.000_r8*rxt(k,225)*y(k,24)
         mat(k,616) = (rxt(k,195)+rxt(k,196))*y(k,5)
         mat(k,35) = 2.000_r8*rxt(k,209)*y(k,95)
         mat(k,260) = rxt(k,204)*y(k,69) + rxt(k,210)*y(k,95) + rxt(k,205)*y(k,99)
         mat(k,754) = rxt(k,199)*y(k,5)
         mat(k,687) = rxt(k,202)*y(k,5) + rxt(k,204)*y(k,45)
         mat(k,435) = rxt(k,253)*y(k,5)
         mat(k,729) = rxt(k,181)*y(k,9) + rxt(k,207)*y(k,10) + 2.000_r8*rxt(k,230) &
                      *y(k,16) + rxt(k,231)*y(k,18) + 3.000_r8*rxt(k,208)*y(k,24) &
                      + 2.000_r8*rxt(k,209)*y(k,42) + rxt(k,210)*y(k,45)
         mat(k,857) = rxt(k,203)*y(k,5) + 2.000_r8*rxt(k,219)*y(k,16) + rxt(k,221) &
                      *y(k,18) + 3.000_r8*rxt(k,226)*y(k,24) + rxt(k,205)*y(k,45)
         mat(k,909) = rxt(k,197)*y(k,28)
         mat(k,596) = rxt(k,197)*y(k,5)
         mat(k,448) = (rxt(k,275)+rxt(k,280))*y(k,53)
         mat(k,247) = (rxt(k,275)+rxt(k,280))*y(k,49)
         mat(k,931) = -(4._r8*rxt(k,194)*y(k,5) + (rxt(k,195) + rxt(k,196) + rxt(k,197) &
                      ) * y(k,28) + rxt(k,198)*y(k,90) + rxt(k,199)*y(k,61) + rxt(k,200) &
                      *y(k,62) + rxt(k,202)*y(k,69) + rxt(k,203)*y(k,99) + rxt(k,253) &
                      *y(k,79))
         mat(k,620) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,5)
         mat(k,567) = -rxt(k,198)*y(k,5)
         mat(k,758) = -rxt(k,199)*y(k,5)
         mat(k,649) = -rxt(k,200)*y(k,5)
         mat(k,691) = -rxt(k,202)*y(k,5)
         mat(k,861) = -rxt(k,203)*y(k,5)
         mat(k,437) = -rxt(k,253)*y(k,5)
         mat(k,812) = rxt(k,193)*y(k,71)
         mat(k,209) = rxt(k,201)*y(k,69)
         mat(k,262) = rxt(k,211)*y(k,95)
         mat(k,254) = rxt(k,206)*y(k,69)
         mat(k,691) = mat(k,691) + rxt(k,201)*y(k,6) + rxt(k,206)*y(k,53)
         mat(k,541) = rxt(k,193)*y(k,3)
         mat(k,733) = rxt(k,211)*y(k,45)
         mat(k,202) = -(rxt(k,201)*y(k,69))
         mat(k,656) = -rxt(k,201)*y(k,6)
         mat(k,911) = rxt(k,200)*y(k,62)
         mat(k,625) = rxt(k,200)*y(k,5)
         mat(k,29) = -(rxt(k,180)*y(k,95))
         mat(k,693) = -rxt(k,180)*y(k,8)
         mat(k,55) = -(rxt(k,181)*y(k,95))
         mat(k,698) = -rxt(k,181)*y(k,9)
         mat(k,60) = -(rxt(k,207)*y(k,95))
         mat(k,699) = -rxt(k,207)*y(k,10)
         mat(k,36) = -(rxt(k,182)*y(k,95))
         mat(k,695) = -rxt(k,182)*y(k,11)
         mat(k,65) = -(rxt(k,183)*y(k,95))
         mat(k,700) = -rxt(k,183)*y(k,12)
         mat(k,40) = -(rxt(k,184)*y(k,95))
         mat(k,696) = -rxt(k,184)*y(k,13)
         mat(k,70) = -(rxt(k,185)*y(k,95))
         mat(k,701) = -rxt(k,185)*y(k,14)
         mat(k,44) = -(rxt(k,186)*y(k,95))
         mat(k,697) = -rxt(k,186)*y(k,15)
         mat(k,174) = -(rxt(k,218)*y(k,25) + rxt(k,219)*y(k,99) + rxt(k,230)*y(k,95))
         mat(k,763) = -rxt(k,218)*y(k,16)
         mat(k,827) = -rxt(k,219)*y(k,16)
         mat(k,710) = -rxt(k,230)*y(k,16)
         mat(k,475) = -(rxt(k,155)*y(k,25) + rxt(k,191)*y(k,3) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,69) + rxt(k,237)*y(k,99))
         mat(k,775) = -rxt(k,155)*y(k,17)
         mat(k,797) = -rxt(k,191)*y(k,17)
         mat(k,577) = -rxt(k,235)*y(k,17)
         mat(k,676) = -rxt(k,236)*y(k,17)
         mat(k,846) = -rxt(k,237)*y(k,17)
         mat(k,413) = rxt(k,162)*y(k,28) + rxt(k,239)*y(k,61)
         mat(k,158) = .300_r8*rxt(k,240)*y(k,99)
         mat(k,397) = (rxt(k,243)+rxt(k,244))*y(k,95)
         mat(k,605) = rxt(k,162)*y(k,21)
         mat(k,743) = rxt(k,239)*y(k,21)
         mat(k,718) = (rxt(k,243)+rxt(k,244))*y(k,23)
         mat(k,846) = mat(k,846) + .300_r8*rxt(k,240)*y(k,22)
         mat(k,210) = -(rxt(k,220)*y(k,25) + rxt(k,221)*y(k,99) + rxt(k,231)*y(k,95))
         mat(k,764) = -rxt(k,220)*y(k,18)
         mat(k,829) = -rxt(k,221)*y(k,18)
         mat(k,711) = -rxt(k,231)*y(k,18)
         mat(k,48) = -(rxt(k,222)*y(k,99))
         mat(k,814) = -rxt(k,222)*y(k,19)
         mat(k,147) = -(rxt(k,223)*y(k,25) + rxt(k,224)*y(k,99))
         mat(k,761) = -rxt(k,223)*y(k,20)
         mat(k,823) = -rxt(k,224)*y(k,20)
         mat(k,411) = -(rxt(k,162)*y(k,28) + rxt(k,238)*y(k,90) + rxt(k,239)*y(k,61))
         mat(k,601) = -rxt(k,162)*y(k,21)
         mat(k,549) = -rxt(k,238)*y(k,21)
         mat(k,741) = -rxt(k,239)*y(k,21)
         mat(k,156) = .700_r8*rxt(k,240)*y(k,99)
         mat(k,394) = rxt(k,156)*y(k,25) + rxt(k,212)*y(k,39) + rxt(k,242)*y(k,95) &
                      + rxt(k,241)*y(k,99)
         mat(k,772) = rxt(k,156)*y(k,23)
         mat(k,292) = rxt(k,212)*y(k,23)
         mat(k,715) = rxt(k,242)*y(k,23)
         mat(k,842) = .700_r8*rxt(k,240)*y(k,22) + rxt(k,241)*y(k,23)
         mat(k,155) = -(rxt(k,240)*y(k,99))
         mat(k,824) = -rxt(k,240)*y(k,22)
         mat(k,410) = rxt(k,238)*y(k,90)
         mat(k,543) = rxt(k,238)*y(k,21)
         mat(k,393) = -(rxt(k,156)*y(k,25) + rxt(k,212)*y(k,39) + rxt(k,241)*y(k,99) &
                      + (rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,95))
         mat(k,771) = -rxt(k,156)*y(k,23)
         mat(k,291) = -rxt(k,212)*y(k,23)
         mat(k,841) = -rxt(k,241)*y(k,23)
         mat(k,714) = -(rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,23)
         mat(k,161) = -(rxt(k,208)*y(k,95) + rxt(k,225)*y(k,25) + rxt(k,226)*y(k,99))
         mat(k,709) = -rxt(k,208)*y(k,24)
         mat(k,762) = -rxt(k,225)*y(k,24)
         mat(k,825) = -rxt(k,226)*y(k,24)
         mat(k,785) = -(rxt(k,155)*y(k,17) + rxt(k,156)*y(k,23) + rxt(k,157)*y(k,41) &
                      + rxt(k,158)*y(k,43) + (rxt(k,159) + rxt(k,160)) * y(k,90) &
                      + rxt(k,161)*y(k,71) + rxt(k,168)*y(k,29) + rxt(k,177)*y(k,54) &
                      + rxt(k,218)*y(k,16) + rxt(k,220)*y(k,18) + rxt(k,223)*y(k,20) &
                      + rxt(k,225)*y(k,24))
         mat(k,485) = -rxt(k,155)*y(k,25)
         mat(k,405) = -rxt(k,156)*y(k,25)
         mat(k,879) = -rxt(k,157)*y(k,25)
         mat(k,232) = -rxt(k,158)*y(k,25)
         mat(k,562) = -(rxt(k,159) + rxt(k,160)) * y(k,25)
         mat(k,537) = -rxt(k,161)*y(k,25)
         mat(k,314) = -rxt(k,168)*y(k,25)
         mat(k,269) = -rxt(k,177)*y(k,25)
         mat(k,177) = -rxt(k,218)*y(k,25)
         mat(k,215) = -rxt(k,220)*y(k,25)
         mat(k,152) = -rxt(k,223)*y(k,25)
         mat(k,164) = -rxt(k,225)*y(k,25)
         mat(k,926) = rxt(k,196)*y(k,28)
         mat(k,31) = 4.000_r8*rxt(k,180)*y(k,95)
         mat(k,58) = rxt(k,181)*y(k,95)
         mat(k,39) = 2.000_r8*rxt(k,182)*y(k,95)
         mat(k,69) = 2.000_r8*rxt(k,183)*y(k,95)
         mat(k,43) = 2.000_r8*rxt(k,184)*y(k,95)
         mat(k,74) = rxt(k,185)*y(k,95)
         mat(k,47) = 2.000_r8*rxt(k,186)*y(k,95)
         mat(k,49) = 3.000_r8*rxt(k,222)*y(k,99)
         mat(k,152) = mat(k,152) + rxt(k,224)*y(k,99)
         mat(k,419) = rxt(k,162)*y(k,28)
         mat(k,615) = rxt(k,196)*y(k,5) + rxt(k,162)*y(k,21) + (4.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,165))*y(k,28) + rxt(k,167)*y(k,61) + rxt(k,172) &
                      *y(k,69) + rxt(k,254)*y(k,79) + rxt(k,173)*y(k,99)
         mat(k,93) = rxt(k,217)*y(k,95)
         mat(k,88) = rxt(k,232)*y(k,95) + rxt(k,227)*y(k,99)
         mat(k,98) = rxt(k,233)*y(k,95) + rxt(k,228)*y(k,99)
         mat(k,116) = rxt(k,234)*y(k,95) + rxt(k,229)*y(k,99)
         mat(k,462) = rxt(k,175)*y(k,69) + rxt(k,187)*y(k,95) + rxt(k,176)*y(k,99)
         mat(k,753) = rxt(k,167)*y(k,28)
         mat(k,686) = rxt(k,172)*y(k,28) + rxt(k,175)*y(k,49)
         mat(k,434) = rxt(k,254)*y(k,28)
         mat(k,728) = 4.000_r8*rxt(k,180)*y(k,8) + rxt(k,181)*y(k,9) &
                      + 2.000_r8*rxt(k,182)*y(k,11) + 2.000_r8*rxt(k,183)*y(k,12) &
                      + 2.000_r8*rxt(k,184)*y(k,13) + rxt(k,185)*y(k,14) &
                      + 2.000_r8*rxt(k,186)*y(k,15) + rxt(k,217)*y(k,34) + rxt(k,232) &
                      *y(k,46) + rxt(k,233)*y(k,47) + rxt(k,234)*y(k,48) + rxt(k,187) &
                      *y(k,49)
         mat(k,856) = 3.000_r8*rxt(k,222)*y(k,19) + rxt(k,224)*y(k,20) + rxt(k,173) &
                      *y(k,28) + rxt(k,227)*y(k,46) + rxt(k,228)*y(k,47) + rxt(k,229) &
                      *y(k,48) + rxt(k,176)*y(k,49)
         mat(k,760) = rxt(k,168)*y(k,29)
         mat(k,595) = 2.000_r8*rxt(k,164)*y(k,28)
         mat(k,306) = rxt(k,168)*y(k,25) + (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,49)
         mat(k,447) = (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,29) + (rxt(k,268) &
                       +rxt(k,274)+rxt(k,279))*y(k,54)
         mat(k,264) = (rxt(k,268)+rxt(k,274)+rxt(k,279))*y(k,49)
         mat(k,594) = 2.000_r8*rxt(k,189)*y(k,28)
         mat(k,610) = -(rxt(k,162)*y(k,21) + (4._r8*rxt(k,163) + 4._r8*rxt(k,164) &
                      + 4._r8*rxt(k,165) + 4._r8*rxt(k,189)) * y(k,28) + rxt(k,166) &
                      *y(k,90) + rxt(k,167)*y(k,61) + rxt(k,169)*y(k,62) + rxt(k,172) &
                      *y(k,69) + (rxt(k,173) + rxt(k,174)) * y(k,99) + (rxt(k,195) &
                      + rxt(k,196) + rxt(k,197)) * y(k,5) + rxt(k,254)*y(k,79))
         mat(k,416) = -rxt(k,162)*y(k,28)
         mat(k,557) = -rxt(k,166)*y(k,28)
         mat(k,748) = -rxt(k,167)*y(k,28)
         mat(k,639) = -rxt(k,169)*y(k,28)
         mat(k,681) = -rxt(k,172)*y(k,28)
         mat(k,851) = -(rxt(k,173) + rxt(k,174)) * y(k,28)
         mat(k,921) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,28)
         mat(k,430) = -rxt(k,254)*y(k,28)
         mat(k,780) = rxt(k,177)*y(k,54) + rxt(k,161)*y(k,71) + rxt(k,160)*y(k,90)
         mat(k,311) = rxt(k,170)*y(k,69)
         mat(k,458) = rxt(k,188)*y(k,95)
         mat(k,267) = rxt(k,177)*y(k,25) + rxt(k,178)*y(k,69) + rxt(k,179)*y(k,99)
         mat(k,681) = mat(k,681) + rxt(k,170)*y(k,29) + rxt(k,178)*y(k,54)
         mat(k,532) = rxt(k,161)*y(k,25)
         mat(k,137) = rxt(k,259)*y(k,79)
         mat(k,430) = mat(k,430) + rxt(k,259)*y(k,73)
         mat(k,557) = mat(k,557) + rxt(k,160)*y(k,25)
         mat(k,723) = rxt(k,188)*y(k,49)
         mat(k,851) = mat(k,851) + rxt(k,179)*y(k,54)
         mat(k,308) = -(rxt(k,168)*y(k,25) + rxt(k,170)*y(k,69) + rxt(k,171)*y(k,99) &
                      + (rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,49))
         mat(k,768) = -rxt(k,168)*y(k,29)
         mat(k,667) = -rxt(k,170)*y(k,29)
         mat(k,837) = -rxt(k,171)*y(k,29)
         mat(k,451) = -(rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,29)
         mat(k,600) = rxt(k,169)*y(k,62)
         mat(k,628) = rxt(k,169)*y(k,28)
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
         mat(k,282) = -(rxt(k,245)*y(k,99))
         mat(k,834) = -rxt(k,245)*y(k,31)
         mat(k,793) = rxt(k,191)*y(k,17)
         mat(k,470) = rxt(k,191)*y(k,3) + rxt(k,155)*y(k,25) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,69) + rxt(k,237)*y(k,99)
         mat(k,148) = rxt(k,223)*y(k,25)
         mat(k,767) = rxt(k,155)*y(k,17) + rxt(k,223)*y(k,20)
         mat(k,197) = rxt(k,299)*y(k,100)
         mat(k,571) = rxt(k,235)*y(k,17)
         mat(k,665) = rxt(k,236)*y(k,17) + rxt(k,248)*y(k,74)
         mat(k,141) = rxt(k,248)*y(k,69) + rxt(k,249)*y(k,99)
         mat(k,834) = mat(k,834) + rxt(k,237)*y(k,17) + rxt(k,249)*y(k,74)
         mat(k,373) = rxt(k,299)*y(k,32)
         mat(k,196) = -(rxt(k,299)*y(k,100))
         mat(k,370) = -rxt(k,299)*y(k,32)
         mat(k,281) = rxt(k,245)*y(k,99)
         mat(k,828) = rxt(k,245)*y(k,31)
         mat(k,82) = -(rxt(k,216)*y(k,95))
         mat(k,703) = -rxt(k,216)*y(k,33)
         mat(k,56) = rxt(k,181)*y(k,95)
         mat(k,61) = rxt(k,207)*y(k,95)
         mat(k,66) = rxt(k,183)*y(k,95)
         mat(k,41) = 2.000_r8*rxt(k,184)*y(k,95)
         mat(k,71) = 2.000_r8*rxt(k,185)*y(k,95)
         mat(k,45) = rxt(k,186)*y(k,95)
         mat(k,33) = 2.000_r8*rxt(k,209)*y(k,95)
         mat(k,94) = rxt(k,233)*y(k,95) + rxt(k,228)*y(k,99)
         mat(k,112) = rxt(k,234)*y(k,95) + rxt(k,229)*y(k,99)
         mat(k,703) = mat(k,703) + rxt(k,181)*y(k,9) + rxt(k,207)*y(k,10) + rxt(k,183) &
                      *y(k,12) + 2.000_r8*rxt(k,184)*y(k,13) + 2.000_r8*rxt(k,185) &
                      *y(k,14) + rxt(k,186)*y(k,15) + 2.000_r8*rxt(k,209)*y(k,42) &
                      + rxt(k,233)*y(k,47) + rxt(k,234)*y(k,48)
         mat(k,815) = rxt(k,228)*y(k,47) + rxt(k,229)*y(k,48)
         mat(k,90) = -(rxt(k,217)*y(k,95))
         mat(k,705) = -rxt(k,217)*y(k,34)
         mat(k,37) = rxt(k,182)*y(k,95)
         mat(k,67) = rxt(k,183)*y(k,95)
         mat(k,86) = rxt(k,232)*y(k,95) + rxt(k,227)*y(k,99)
         mat(k,705) = mat(k,705) + rxt(k,182)*y(k,11) + rxt(k,183)*y(k,12) &
                      + rxt(k,232)*y(k,46)
         mat(k,817) = rxt(k,227)*y(k,46)
         mat(k,128) = -(rxt(k,246)*y(k,63) + (rxt(k,247) + rxt(k,261)) * y(k,99))
         mat(k,570) = -rxt(k,246)*y(k,35)
         mat(k,821) = -(rxt(k,247) + rxt(k,261)) * y(k,35)
         mat(k,290) = -(rxt(k,212)*y(k,23) + rxt(k,213)*y(k,41) + rxt(k,214)*y(k,103) &
                      + rxt(k,215)*y(k,51))
         mat(k,390) = -rxt(k,212)*y(k,39)
         mat(k,864) = -rxt(k,213)*y(k,39)
         mat(k,936) = -rxt(k,214)*y(k,39)
         mat(k,887) = -rxt(k,215)*y(k,39)
         mat(k,62) = rxt(k,207)*y(k,95)
         mat(k,72) = rxt(k,185)*y(k,95)
         mat(k,83) = 2.000_r8*rxt(k,216)*y(k,95)
         mat(k,91) = rxt(k,217)*y(k,95)
         mat(k,713) = rxt(k,207)*y(k,10) + rxt(k,185)*y(k,14) + 2.000_r8*rxt(k,216) &
                      *y(k,33) + rxt(k,217)*y(k,34)
         mat(k,439) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,90) + rxt(k,116) &
                      *y(k,70) + rxt(k,119)*y(k,71))
         mat(k,550) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,40)
         mat(k,507) = -rxt(k,116)*y(k,40)
         mat(k,527) = -rxt(k,119)*y(k,40)
         mat(k,473) = rxt(k,237)*y(k,99)
         mat(k,395) = rxt(k,243)*y(k,95)
         mat(k,773) = rxt(k,157)*y(k,41)
         mat(k,293) = rxt(k,213)*y(k,41)
         mat(k,867) = rxt(k,157)*y(k,25) + rxt(k,213)*y(k,39) + rxt(k,111)*y(k,69) &
                      + rxt(k,94)*y(k,95) + rxt(k,120)*y(k,99)
         mat(k,257) = rxt(k,211)*y(k,95)
         mat(k,452) = rxt(k,188)*y(k,95)
         mat(k,350) = rxt(k,143)*y(k,99)
         mat(k,674) = rxt(k,111)*y(k,41) + rxt(k,123)*y(k,99)
         mat(k,144) = rxt(k,249)*y(k,99)
         mat(k,237) = rxt(k,255)*y(k,99)
         mat(k,426) = rxt(k,260)*y(k,99)
         mat(k,716) = rxt(k,243)*y(k,23) + rxt(k,94)*y(k,41) + rxt(k,211)*y(k,45) &
                      + rxt(k,188)*y(k,49)
         mat(k,844) = rxt(k,237)*y(k,17) + rxt(k,120)*y(k,41) + rxt(k,143)*y(k,55) &
                      + rxt(k,123)*y(k,69) + rxt(k,249)*y(k,74) + rxt(k,255)*y(k,77) &
                      + rxt(k,260)*y(k,79)
         mat(k,882) = -(rxt(k,94)*y(k,95) + rxt(k,111)*y(k,69) + rxt(k,120)*y(k,99) &
                      + rxt(k,157)*y(k,25) + rxt(k,213)*y(k,39))
         mat(k,731) = -rxt(k,94)*y(k,41)
         mat(k,689) = -rxt(k,111)*y(k,41)
         mat(k,859) = -rxt(k,120)*y(k,41)
         mat(k,788) = -rxt(k,157)*y(k,41)
         mat(k,296) = -rxt(k,213)*y(k,41)
         mat(k,407) = rxt(k,244)*y(k,95)
         mat(k,445) = rxt(k,113)*y(k,90)
         mat(k,565) = rxt(k,113)*y(k,40)
         mat(k,731) = mat(k,731) + rxt(k,244)*y(k,23)
         mat(k,32) = -(rxt(k,209)*y(k,95))
         mat(k,694) = -rxt(k,209)*y(k,42)
         mat(k,228) = -(rxt(k,112)*y(k,69) + rxt(k,121)*y(k,99) + rxt(k,158)*y(k,25))
         mat(k,658) = -rxt(k,112)*y(k,43)
         mat(k,830) = -rxt(k,121)*y(k,43)
         mat(k,765) = -rxt(k,158)*y(k,43)
         mat(k,545) = 2.000_r8*rxt(k,127)*y(k,90)
         mat(k,830) = mat(k,830) + 2.000_r8*rxt(k,126)*y(k,99)
         mat(k,107) = rxt(k,262)*y(k,103)
         mat(k,933) = rxt(k,262)*y(k,81)
         mat(k,256) = -(rxt(k,204)*y(k,69) + rxt(k,205)*y(k,99) + (rxt(k,210) &
                      + rxt(k,211)) * y(k,95))
         mat(k,662) = -rxt(k,204)*y(k,45)
         mat(k,832) = -rxt(k,205)*y(k,45)
         mat(k,712) = -(rxt(k,210) + rxt(k,211)) * y(k,45)
         mat(k,792) = rxt(k,191)*y(k,17) + rxt(k,192)*y(k,90)
         mat(k,469) = rxt(k,191)*y(k,3)
         mat(k,547) = rxt(k,192)*y(k,3)
         mat(k,85) = -(rxt(k,227)*y(k,99) + rxt(k,232)*y(k,95))
         mat(k,816) = -rxt(k,227)*y(k,46)
         mat(k,704) = -rxt(k,232)*y(k,46)
         mat(k,95) = -(rxt(k,228)*y(k,99) + rxt(k,233)*y(k,95))
         mat(k,818) = -rxt(k,228)*y(k,47)
         mat(k,706) = -rxt(k,233)*y(k,47)
         mat(k,113) = -(rxt(k,229)*y(k,99) + rxt(k,234)*y(k,95))
         mat(k,820) = -rxt(k,229)*y(k,48)
         mat(k,708) = -rxt(k,234)*y(k,48)
         mat(k,453) = -(rxt(k,175)*y(k,69) + rxt(k,176)*y(k,99) + (rxt(k,187) &
                      + rxt(k,188)) * y(k,95) + (rxt(k,268) + rxt(k,274) + rxt(k,279) &
                      ) * y(k,54) + (rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,29) &
                      + (rxt(k,275) + rxt(k,280)) * y(k,53))
         mat(k,675) = -rxt(k,175)*y(k,49)
         mat(k,845) = -rxt(k,176)*y(k,49)
         mat(k,717) = -(rxt(k,187) + rxt(k,188)) * y(k,49)
         mat(k,266) = -(rxt(k,268) + rxt(k,274) + rxt(k,279)) * y(k,49)
         mat(k,309) = -(rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,49)
         mat(k,249) = -(rxt(k,275) + rxt(k,280)) * y(k,49)
         mat(k,175) = rxt(k,218)*y(k,25)
         mat(k,474) = rxt(k,155)*y(k,25)
         mat(k,212) = rxt(k,220)*y(k,25)
         mat(k,150) = 2.000_r8*rxt(k,223)*y(k,25)
         mat(k,396) = rxt(k,156)*y(k,25)
         mat(k,162) = rxt(k,225)*y(k,25)
         mat(k,774) = rxt(k,218)*y(k,16) + rxt(k,155)*y(k,17) + rxt(k,220)*y(k,18) &
                      + 2.000_r8*rxt(k,223)*y(k,20) + rxt(k,156)*y(k,23) + rxt(k,225) &
                      *y(k,24) + rxt(k,157)*y(k,41) + rxt(k,158)*y(k,43) + rxt(k,177) &
                      *y(k,54) + rxt(k,159)*y(k,90)
         mat(k,604) = rxt(k,174)*y(k,99)
         mat(k,868) = rxt(k,157)*y(k,25)
         mat(k,229) = rxt(k,158)*y(k,25)
         mat(k,266) = mat(k,266) + rxt(k,177)*y(k,25)
         mat(k,551) = rxt(k,159)*y(k,25)
         mat(k,845) = mat(k,845) + rxt(k,174)*y(k,28)
         mat(k,387) = rxt(k,212)*y(k,39)
         mat(k,289) = rxt(k,212)*y(k,23) + rxt(k,213)*y(k,41) + rxt(k,215)*y(k,51) &
                      + rxt(k,214)*y(k,103)
         mat(k,863) = rxt(k,213)*y(k,39)
         mat(k,886) = rxt(k,215)*y(k,39)
         mat(k,935) = rxt(k,214)*y(k,39)
         mat(k,906) = -(rxt(k,152)*y(k,99) + rxt(k,215)*y(k,39))
         mat(k,860) = -rxt(k,152)*y(k,51)
         mat(k,297) = -rxt(k,215)*y(k,51)
         mat(k,489) = rxt(k,235)*y(k,63)
         mat(k,316) = (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,49)
         mat(k,133) = rxt(k,246)*y(k,63)
         mat(k,466) = (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,29)
         mat(k,648) = rxt(k,151)*y(k,99)
         mat(k,591) = rxt(k,235)*y(k,17) + rxt(k,246)*y(k,35)
         mat(k,860) = mat(k,860) + rxt(k,151)*y(k,62)
         mat(k,167) = -(rxt(k,128)*y(k,99))
         mat(k,826) = -rxt(k,128)*y(k,52)
         mat(k,624) = rxt(k,149)*y(k,90)
         mat(k,544) = rxt(k,149)*y(k,62)
         mat(k,248) = -(rxt(k,206)*y(k,69) + (rxt(k,275) + rxt(k,280)) * y(k,49))
         mat(k,661) = -rxt(k,206)*y(k,53)
         mat(k,449) = -(rxt(k,275) + rxt(k,280)) * y(k,53)
         mat(k,912) = rxt(k,198)*y(k,90)
         mat(k,546) = rxt(k,198)*y(k,5)
         mat(k,265) = -(rxt(k,177)*y(k,25) + rxt(k,178)*y(k,69) + rxt(k,179)*y(k,99) &
                      + (rxt(k,268) + rxt(k,274) + rxt(k,279)) * y(k,49))
         mat(k,766) = -rxt(k,177)*y(k,54)
         mat(k,663) = -rxt(k,178)*y(k,54)
         mat(k,833) = -rxt(k,179)*y(k,54)
         mat(k,450) = -(rxt(k,268) + rxt(k,274) + rxt(k,279)) * y(k,54)
         mat(k,598) = rxt(k,166)*y(k,90)
         mat(k,307) = rxt(k,171)*y(k,99)
         mat(k,548) = rxt(k,166)*y(k,28)
         mat(k,833) = mat(k,833) + rxt(k,171)*y(k,29)
         mat(k,347) = -(rxt(k,131)*y(k,61) + (rxt(k,132) + rxt(k,133) + rxt(k,134) &
                      ) * y(k,62) + rxt(k,135)*y(k,70) + rxt(k,143)*y(k,99) + rxt(k,296) &
                      *y(k,98))
         mat(k,738) = -rxt(k,131)*y(k,55)
         mat(k,629) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,55)
         mat(k,503) = -rxt(k,135)*y(k,55)
         mat(k,838) = -rxt(k,143)*y(k,55)
         mat(k,360) = -rxt(k,296)*y(k,55)
         mat(k,670) = rxt(k,129)*y(k,91) + rxt(k,293)*y(k,94)
         mat(k,503) = mat(k,503) + rxt(k,294)*y(k,94)
         mat(k,325) = 1.100_r8*rxt(k,289)*y(k,92) + .200_r8*rxt(k,287)*y(k,93)
         mat(k,335) = rxt(k,129)*y(k,69)
         mat(k,223) = 1.100_r8*rxt(k,289)*y(k,89)
         mat(k,245) = .200_r8*rxt(k,287)*y(k,89)
         mat(k,276) = rxt(k,293)*y(k,69) + rxt(k,294)*y(k,70)
         mat(k,103) = -((rxt(k,147) + rxt(k,148)) * y(k,95))
         mat(k,707) = -(rxt(k,147) + rxt(k,148)) * y(k,56)
         mat(k,342) = rxt(k,132)*y(k,62)
         mat(k,622) = rxt(k,132)*y(k,55)
         mat(k,623) = rxt(k,150)*y(k,63)
         mat(k,569) = rxt(k,150)*y(k,62)
         mat(k,752) = -(rxt(k,131)*y(k,55) + rxt(k,140)*y(k,63) + rxt(k,144)*y(k,90) &
                      + rxt(k,145)*y(k,71) + rxt(k,146)*y(k,69) + rxt(k,167)*y(k,28) &
                      + rxt(k,199)*y(k,5) + rxt(k,239)*y(k,21) + rxt(k,298)*y(k,98))
         mat(k,355) = -rxt(k,131)*y(k,61)
         mat(k,586) = -rxt(k,140)*y(k,61)
         mat(k,561) = -rxt(k,144)*y(k,61)
         mat(k,536) = -rxt(k,145)*y(k,61)
         mat(k,685) = -rxt(k,146)*y(k,61)
         mat(k,614) = -rxt(k,167)*y(k,61)
         mat(k,925) = -rxt(k,199)*y(k,61)
         mat(k,418) = -rxt(k,239)*y(k,61)
         mat(k,368) = -rxt(k,298)*y(k,61)
         mat(k,355) = mat(k,355) + 2.000_r8*rxt(k,133)*y(k,62) + rxt(k,135)*y(k,70) &
                      + rxt(k,143)*y(k,99)
         mat(k,106) = 2.000_r8*rxt(k,147)*y(k,95)
         mat(k,643) = 2.000_r8*rxt(k,133)*y(k,55) + rxt(k,136)*y(k,69) + rxt(k,256) &
                      *y(k,79)
         mat(k,685) = mat(k,685) + rxt(k,136)*y(k,62)
         mat(k,515) = rxt(k,135)*y(k,55) + rxt(k,130)*y(k,91)
         mat(k,433) = rxt(k,256)*y(k,62)
         mat(k,341) = rxt(k,130)*y(k,70)
         mat(k,727) = 2.000_r8*rxt(k,147)*y(k,56)
         mat(k,855) = rxt(k,143)*y(k,55)
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
         mat(k,640) = -((rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,55) + (rxt(k,136) &
                      + rxt(k,138)) * y(k,69) + rxt(k,137)*y(k,71) + rxt(k,149) &
                      *y(k,90) + rxt(k,150)*y(k,63) + rxt(k,151)*y(k,99) + rxt(k,169) &
                      *y(k,28) + rxt(k,200)*y(k,5) + rxt(k,256)*y(k,79))
         mat(k,352) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,62)
         mat(k,682) = -(rxt(k,136) + rxt(k,138)) * y(k,62)
         mat(k,533) = -rxt(k,137)*y(k,62)
         mat(k,558) = -rxt(k,149)*y(k,62)
         mat(k,583) = -rxt(k,150)*y(k,62)
         mat(k,852) = -rxt(k,151)*y(k,62)
         mat(k,611) = -rxt(k,169)*y(k,62)
         mat(k,922) = -rxt(k,200)*y(k,62)
         mat(k,431) = -rxt(k,256)*y(k,62)
         mat(k,922) = mat(k,922) + rxt(k,199)*y(k,61)
         mat(k,417) = rxt(k,239)*y(k,61)
         mat(k,611) = mat(k,611) + rxt(k,167)*y(k,61)
         mat(k,171) = rxt(k,128)*y(k,99)
         mat(k,749) = rxt(k,199)*y(k,5) + rxt(k,239)*y(k,21) + rxt(k,167)*y(k,28) &
                      + 2.000_r8*rxt(k,140)*y(k,63) + rxt(k,146)*y(k,69) + rxt(k,145) &
                      *y(k,71) + rxt(k,144)*y(k,90)
         mat(k,583) = mat(k,583) + 2.000_r8*rxt(k,140)*y(k,61) + rxt(k,141)*y(k,69) &
                      + rxt(k,139)*y(k,90) + rxt(k,142)*y(k,99)
         mat(k,682) = mat(k,682) + rxt(k,146)*y(k,61) + rxt(k,141)*y(k,63)
         mat(k,533) = mat(k,533) + rxt(k,145)*y(k,61)
         mat(k,558) = mat(k,558) + rxt(k,144)*y(k,61) + rxt(k,139)*y(k,63)
         mat(k,852) = mat(k,852) + rxt(k,128)*y(k,52) + rxt(k,142)*y(k,63)
         mat(k,581) = -(rxt(k,139)*y(k,90) + rxt(k,140)*y(k,61) + rxt(k,141)*y(k,69) &
                      + rxt(k,142)*y(k,99) + rxt(k,150)*y(k,62) + rxt(k,235)*y(k,17) &
                      + rxt(k,246)*y(k,35))
         mat(k,556) = -rxt(k,139)*y(k,63)
         mat(k,747) = -rxt(k,140)*y(k,63)
         mat(k,680) = -rxt(k,141)*y(k,63)
         mat(k,850) = -rxt(k,142)*y(k,63)
         mat(k,638) = -rxt(k,150)*y(k,63)
         mat(k,479) = -rxt(k,235)*y(k,63)
         mat(k,131) = -rxt(k,246)*y(k,63)
         mat(k,204) = rxt(k,201)*y(k,69)
         mat(k,779) = rxt(k,168)*y(k,29)
         mat(k,310) = rxt(k,168)*y(k,25) + rxt(k,170)*y(k,69) + rxt(k,171)*y(k,99)
         mat(k,294) = rxt(k,215)*y(k,51)
         mat(k,896) = rxt(k,215)*y(k,39) + rxt(k,152)*y(k,99)
         mat(k,638) = mat(k,638) + rxt(k,138)*y(k,69) + rxt(k,137)*y(k,71)
         mat(k,680) = mat(k,680) + rxt(k,201)*y(k,6) + rxt(k,170)*y(k,29) + rxt(k,138) &
                      *y(k,62)
         mat(k,531) = rxt(k,137)*y(k,62)
         mat(k,850) = mat(k,850) + rxt(k,171)*y(k,29) + rxt(k,152)*y(k,51)
         mat(k,683) = -(rxt(k,108)*y(k,71) + 4._r8*rxt(k,109)*y(k,69) + rxt(k,110) &
                      *y(k,70) + rxt(k,111)*y(k,41) + rxt(k,112)*y(k,43) + rxt(k,117) &
                      *y(k,90) + rxt(k,123)*y(k,99) + (rxt(k,136) + rxt(k,138) &
                      ) * y(k,62) + rxt(k,141)*y(k,63) + rxt(k,146)*y(k,61) + rxt(k,170) &
                      *y(k,29) + rxt(k,172)*y(k,28) + rxt(k,175)*y(k,49) + rxt(k,178) &
                      *y(k,54) + rxt(k,201)*y(k,6) + rxt(k,202)*y(k,5) + rxt(k,204) &
                      *y(k,45) + rxt(k,206)*y(k,53) + rxt(k,236)*y(k,17) + rxt(k,248) &
                      *y(k,74) + (rxt(k,291) + rxt(k,292)) * y(k,92) + rxt(k,293) &
                      *y(k,94))
         mat(k,534) = -rxt(k,108)*y(k,69)
         mat(k,513) = -rxt(k,110)*y(k,69)
         mat(k,876) = -rxt(k,111)*y(k,69)
         mat(k,231) = -rxt(k,112)*y(k,69)
         mat(k,559) = -rxt(k,117)*y(k,69)
         mat(k,853) = -rxt(k,123)*y(k,69)
         mat(k,641) = -(rxt(k,136) + rxt(k,138)) * y(k,69)
         mat(k,584) = -rxt(k,141)*y(k,69)
         mat(k,750) = -rxt(k,146)*y(k,69)
         mat(k,313) = -rxt(k,170)*y(k,69)
         mat(k,612) = -rxt(k,172)*y(k,69)
         mat(k,460) = -rxt(k,175)*y(k,69)
         mat(k,268) = -rxt(k,178)*y(k,69)
         mat(k,206) = -rxt(k,201)*y(k,69)
         mat(k,923) = -rxt(k,202)*y(k,69)
         mat(k,258) = -rxt(k,204)*y(k,69)
         mat(k,250) = -rxt(k,206)*y(k,69)
         mat(k,482) = -rxt(k,236)*y(k,69)
         mat(k,145) = -rxt(k,248)*y(k,69)
         mat(k,227) = -(rxt(k,291) + rxt(k,292)) * y(k,69)
         mat(k,280) = -rxt(k,293)*y(k,69)
         mat(k,443) = rxt(k,115)*y(k,90)
         mat(k,353) = rxt(k,131)*y(k,61) + rxt(k,132)*y(k,62) + rxt(k,135)*y(k,70) &
                      + rxt(k,296)*y(k,98)
         mat(k,750) = mat(k,750) + rxt(k,131)*y(k,55)
         mat(k,641) = mat(k,641) + rxt(k,132)*y(k,55)
         mat(k,513) = mat(k,513) + rxt(k,135)*y(k,55) + rxt(k,250)*y(k,77) &
                      + rxt(k,257)*y(k,79) + rxt(k,295)*y(k,94) + (rxt(k,97)+rxt(k,98)) &
                      *y(k,95) + rxt(k,302)*y(k,100) + rxt(k,306)*y(k,101)
         mat(k,240) = rxt(k,250)*y(k,70)
         mat(k,432) = rxt(k,257)*y(k,70)
         mat(k,329) = rxt(k,287)*y(k,93) + 1.150_r8*rxt(k,288)*y(k,98)
         mat(k,559) = mat(k,559) + rxt(k,115)*y(k,40)
         mat(k,339) = rxt(k,301)*y(k,100)
         mat(k,246) = rxt(k,287)*y(k,89)
         mat(k,280) = mat(k,280) + rxt(k,295)*y(k,70)
         mat(k,725) = (rxt(k,97)+rxt(k,98))*y(k,70)
         mat(k,366) = rxt(k,296)*y(k,55) + 1.150_r8*rxt(k,288)*y(k,89)
         mat(k,853) = mat(k,853) + 2.000_r8*rxt(k,125)*y(k,99)
         mat(k,383) = rxt(k,302)*y(k,70) + rxt(k,301)*y(k,91)
         mat(k,193) = rxt(k,306)*y(k,70)
         mat(k,508) = -(rxt(k,97)*y(k,95) + rxt(k,102)*y(k,96) + rxt(k,110)*y(k,69) &
                      + rxt(k,116)*y(k,40) + rxt(k,130)*y(k,91) + rxt(k,135)*y(k,55) &
                      + rxt(k,250)*y(k,77) + rxt(k,257)*y(k,79) + rxt(k,290)*y(k,92) &
                      + (rxt(k,294) + rxt(k,295)) * y(k,94) + rxt(k,302)*y(k,100) &
                      + rxt(k,306)*y(k,101))
         mat(k,719) = -rxt(k,97)*y(k,70)
         mat(k,78) = -rxt(k,102)*y(k,70)
         mat(k,677) = -rxt(k,110)*y(k,70)
         mat(k,440) = -rxt(k,116)*y(k,70)
         mat(k,338) = -rxt(k,130)*y(k,70)
         mat(k,351) = -rxt(k,135)*y(k,70)
         mat(k,238) = -rxt(k,250)*y(k,70)
         mat(k,427) = -rxt(k,257)*y(k,70)
         mat(k,226) = -rxt(k,290)*y(k,70)
         mat(k,279) = -(rxt(k,294) + rxt(k,295)) * y(k,70)
         mat(k,380) = -rxt(k,302)*y(k,70)
         mat(k,192) = -rxt(k,306)*y(k,70)
         mat(k,798) = rxt(k,193)*y(k,71) + rxt(k,192)*y(k,90)
         mat(k,917) = 2.000_r8*rxt(k,194)*y(k,5) + (rxt(k,196)+rxt(k,197))*y(k,28) &
                      + rxt(k,202)*y(k,69) + rxt(k,198)*y(k,90)
         mat(k,414) = rxt(k,238)*y(k,90)
         mat(k,776) = rxt(k,161)*y(k,71) + rxt(k,159)*y(k,90)
         mat(k,606) = (rxt(k,196)+rxt(k,197))*y(k,5) + (2.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,164))*y(k,28) + rxt(k,172)*y(k,69) + rxt(k,166) &
                      *y(k,90) + rxt(k,174)*y(k,99)
         mat(k,440) = mat(k,440) + rxt(k,119)*y(k,71) + rxt(k,113)*y(k,90)
         mat(k,168) = rxt(k,128)*y(k,99)
         mat(k,351) = mat(k,351) + rxt(k,134)*y(k,62)
         mat(k,104) = rxt(k,148)*y(k,95)
         mat(k,744) = rxt(k,145)*y(k,71) + rxt(k,298)*y(k,98)
         mat(k,635) = rxt(k,134)*y(k,55) + rxt(k,136)*y(k,69) + rxt(k,137)*y(k,71)
         mat(k,578) = rxt(k,141)*y(k,69) + rxt(k,139)*y(k,90)
         mat(k,677) = mat(k,677) + rxt(k,202)*y(k,5) + rxt(k,172)*y(k,28) + rxt(k,136) &
                      *y(k,62) + rxt(k,141)*y(k,63) + 2.000_r8*rxt(k,109)*y(k,69) &
                      + 2.000_r8*rxt(k,108)*y(k,71) + rxt(k,117)*y(k,90) + rxt(k,101) &
                      *y(k,96) + rxt(k,123)*y(k,99)
         mat(k,508) = mat(k,508) + 2.000_r8*rxt(k,102)*y(k,96)
         mat(k,528) = rxt(k,193)*y(k,3) + rxt(k,161)*y(k,25) + rxt(k,119)*y(k,40) &
                      + rxt(k,145)*y(k,61) + rxt(k,137)*y(k,62) + 2.000_r8*rxt(k,108) &
                      *y(k,69) + rxt(k,252)*y(k,77) + rxt(k,258)*y(k,79) &
                      + 2.000_r8*rxt(k,118)*y(k,90) + 2.000_r8*rxt(k,99)*y(k,95) &
                      + rxt(k,124)*y(k,99)
         mat(k,238) = mat(k,238) + rxt(k,252)*y(k,71)
         mat(k,427) = mat(k,427) + rxt(k,258)*y(k,71)
         mat(k,553) = rxt(k,192)*y(k,3) + rxt(k,198)*y(k,5) + rxt(k,238)*y(k,21) &
                      + rxt(k,159)*y(k,25) + rxt(k,166)*y(k,28) + rxt(k,113)*y(k,40) &
                      + rxt(k,139)*y(k,63) + rxt(k,117)*y(k,69) + 2.000_r8*rxt(k,118) &
                      *y(k,71) + 2.000_r8*rxt(k,127)*y(k,90) + rxt(k,122)*y(k,99)
         mat(k,719) = mat(k,719) + rxt(k,148)*y(k,56) + 2.000_r8*rxt(k,99)*y(k,71)
         mat(k,78) = mat(k,78) + rxt(k,101)*y(k,69) + 2.000_r8*rxt(k,102)*y(k,70)
         mat(k,364) = rxt(k,298)*y(k,61)
         mat(k,847) = rxt(k,174)*y(k,28) + rxt(k,128)*y(k,52) + rxt(k,123)*y(k,69) &
                      + rxt(k,124)*y(k,71) + rxt(k,122)*y(k,90)
         mat(k,529) = -(rxt(k,99)*y(k,95) + rxt(k,108)*y(k,69) + rxt(k,118)*y(k,90) &
                      + rxt(k,119)*y(k,40) + rxt(k,124)*y(k,99) + rxt(k,137)*y(k,62) &
                      + rxt(k,145)*y(k,61) + rxt(k,161)*y(k,25) + rxt(k,193)*y(k,3) &
                      + rxt(k,252)*y(k,77) + rxt(k,258)*y(k,79))
         mat(k,720) = -rxt(k,99)*y(k,71)
         mat(k,678) = -rxt(k,108)*y(k,71)
         mat(k,554) = -rxt(k,118)*y(k,71)
         mat(k,441) = -rxt(k,119)*y(k,71)
         mat(k,848) = -rxt(k,124)*y(k,71)
         mat(k,636) = -rxt(k,137)*y(k,71)
         mat(k,745) = -rxt(k,145)*y(k,71)
         mat(k,777) = -rxt(k,161)*y(k,71)
         mat(k,799) = -rxt(k,193)*y(k,71)
         mat(k,239) = -rxt(k,252)*y(k,71)
         mat(k,428) = -rxt(k,258)*y(k,71)
         mat(k,678) = mat(k,678) + rxt(k,110)*y(k,70)
         mat(k,509) = rxt(k,110)*y(k,69)
         mat(k,134) = -(rxt(k,259)*y(k,79))
         mat(k,422) = -rxt(k,259)*y(k,73)
         mat(k,910) = rxt(k,195)*y(k,28)
         mat(k,597) = rxt(k,195)*y(k,5) + 2.000_r8*rxt(k,165)*y(k,28)
         mat(k,139) = -(rxt(k,248)*y(k,69) + rxt(k,249)*y(k,99))
         mat(k,653) = -rxt(k,248)*y(k,74)
         mat(k,822) = -rxt(k,249)*y(k,74)
         mat(k,235) = -(rxt(k,250)*y(k,70) + rxt(k,252)*y(k,71) + rxt(k,255)*y(k,99))
         mat(k,497) = -rxt(k,250)*y(k,77)
         mat(k,524) = -rxt(k,252)*y(k,77)
         mat(k,831) = -rxt(k,255)*y(k,77)
         mat(k,425) = -(rxt(k,253)*y(k,5) + rxt(k,254)*y(k,28) + rxt(k,256)*y(k,62) &
                      + rxt(k,257)*y(k,70) + rxt(k,258)*y(k,71) + rxt(k,259)*y(k,73) &
                      + rxt(k,260)*y(k,99))
         mat(k,914) = -rxt(k,253)*y(k,79)
         mat(k,602) = -rxt(k,254)*y(k,79)
         mat(k,632) = -rxt(k,256)*y(k,79)
         mat(k,506) = -rxt(k,257)*y(k,79)
         mat(k,526) = -rxt(k,258)*y(k,79)
         mat(k,136) = -rxt(k,259)*y(k,79)
         mat(k,843) = -rxt(k,260)*y(k,79)
         mat(k,673) = rxt(k,248)*y(k,74)
         mat(k,506) = mat(k,506) + rxt(k,250)*y(k,77)
         mat(k,526) = mat(k,526) + rxt(k,252)*y(k,77)
         mat(k,143) = rxt(k,248)*y(k,69)
         mat(k,236) = rxt(k,250)*y(k,70) + rxt(k,252)*y(k,71) + rxt(k,255)*y(k,99)
         mat(k,843) = mat(k,843) + rxt(k,255)*y(k,77)
         mat(k,300) = -(rxt(k,251)*y(k,99))
         mat(k,836) = -rxt(k,251)*y(k,80)
         mat(k,913) = rxt(k,253)*y(k,79)
         mat(k,599) = rxt(k,254)*y(k,79)
         mat(k,129) = rxt(k,246)*y(k,63) + (rxt(k,247)+.500_r8*rxt(k,261))*y(k,99)
         mat(k,627) = rxt(k,256)*y(k,79)
         mat(k,572) = rxt(k,246)*y(k,35)
         mat(k,500) = rxt(k,257)*y(k,79)
         mat(k,525) = rxt(k,258)*y(k,79)
         mat(k,135) = rxt(k,259)*y(k,79)
         mat(k,142) = rxt(k,249)*y(k,99)
         mat(k,424) = rxt(k,253)*y(k,5) + rxt(k,254)*y(k,28) + rxt(k,256)*y(k,62) &
                      + rxt(k,257)*y(k,70) + rxt(k,258)*y(k,71) + rxt(k,259)*y(k,73) &
                      + rxt(k,260)*y(k,99)
         mat(k,836) = mat(k,836) + (rxt(k,247)+.500_r8*rxt(k,261))*y(k,35) &
                      + rxt(k,249)*y(k,74) + rxt(k,260)*y(k,79)
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
         mat(k,108) = -(rxt(k,262)*y(k,103))
         mat(k,934) = -rxt(k,262)*y(k,81)
         mat(k,299) = rxt(k,251)*y(k,99)
         mat(k,819) = rxt(k,251)*y(k,80)
         mat(k,323) = -(rxt(k,287)*y(k,93) + rxt(k,288)*y(k,98) + rxt(k,289)*y(k,92))
         mat(k,243) = -rxt(k,287)*y(k,89)
         mat(k,358) = -rxt(k,288)*y(k,89)
         mat(k,221) = -rxt(k,289)*y(k,89)
         mat(k,555) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,40) + rxt(k,117) &
                      *y(k,69) + rxt(k,118)*y(k,71) + rxt(k,122)*y(k,99) &
                      + 4._r8*rxt(k,127)*y(k,90) + rxt(k,139)*y(k,63) + rxt(k,144) &
                      *y(k,61) + rxt(k,149)*y(k,62) + (rxt(k,159) + rxt(k,160) &
                      ) * y(k,25) + rxt(k,166)*y(k,28) + rxt(k,192)*y(k,3) + rxt(k,198) &
                      *y(k,5) + rxt(k,238)*y(k,21))
         mat(k,442) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,90)
         mat(k,679) = -rxt(k,117)*y(k,90)
         mat(k,530) = -rxt(k,118)*y(k,90)
         mat(k,849) = -rxt(k,122)*y(k,90)
         mat(k,580) = -rxt(k,139)*y(k,90)
         mat(k,746) = -rxt(k,144)*y(k,90)
         mat(k,637) = -rxt(k,149)*y(k,90)
         mat(k,778) = -(rxt(k,159) + rxt(k,160)) * y(k,90)
         mat(k,608) = -rxt(k,166)*y(k,90)
         mat(k,800) = -rxt(k,192)*y(k,90)
         mat(k,919) = -rxt(k,198)*y(k,90)
         mat(k,415) = -rxt(k,238)*y(k,90)
         mat(k,800) = mat(k,800) + rxt(k,191)*y(k,17)
         mat(k,919) = mat(k,919) + rxt(k,203)*y(k,99)
         mat(k,478) = rxt(k,191)*y(k,3) + rxt(k,155)*y(k,25) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,69)
         mat(k,213) = rxt(k,220)*y(k,25) + rxt(k,221)*y(k,99)
         mat(k,151) = rxt(k,223)*y(k,25) + rxt(k,224)*y(k,99)
         mat(k,415) = mat(k,415) + rxt(k,162)*y(k,28) + rxt(k,239)*y(k,61)
         mat(k,399) = rxt(k,243)*y(k,95)
         mat(k,778) = mat(k,778) + rxt(k,155)*y(k,17) + rxt(k,220)*y(k,18) &
                      + rxt(k,223)*y(k,20) + rxt(k,158)*y(k,43)
         mat(k,608) = mat(k,608) + rxt(k,162)*y(k,21) + rxt(k,173)*y(k,99)
         mat(k,286) = rxt(k,245)*y(k,99)
         mat(k,130) = .500_r8*rxt(k,261)*y(k,99)
         mat(k,442) = mat(k,442) + rxt(k,116)*y(k,70)
         mat(k,230) = rxt(k,158)*y(k,25) + rxt(k,112)*y(k,69) + rxt(k,121)*y(k,99)
         mat(k,746) = mat(k,746) + rxt(k,239)*y(k,21)
         mat(k,580) = mat(k,580) + rxt(k,235)*y(k,17) + rxt(k,142)*y(k,99)
         mat(k,679) = mat(k,679) + rxt(k,236)*y(k,17) + rxt(k,112)*y(k,43)
         mat(k,510) = rxt(k,116)*y(k,40)
         mat(k,530) = mat(k,530) + rxt(k,124)*y(k,99)
         mat(k,302) = rxt(k,251)*y(k,99)
         mat(k,721) = rxt(k,243)*y(k,23)
         mat(k,849) = mat(k,849) + rxt(k,203)*y(k,5) + rxt(k,221)*y(k,18) + rxt(k,224) &
                      *y(k,20) + rxt(k,173)*y(k,28) + rxt(k,245)*y(k,31) &
                      + .500_r8*rxt(k,261)*y(k,35) + rxt(k,121)*y(k,43) + rxt(k,142) &
                      *y(k,63) + rxt(k,124)*y(k,71) + rxt(k,251)*y(k,80)
         mat(k,334) = -(rxt(k,129)*y(k,69) + rxt(k,130)*y(k,70) + rxt(k,301)*y(k,100))
         mat(k,669) = -rxt(k,129)*y(k,91)
         mat(k,502) = -rxt(k,130)*y(k,91)
         mat(k,375) = -rxt(k,301)*y(k,91)
         mat(k,669) = mat(k,669) + rxt(k,291)*y(k,92)
         mat(k,324) = .900_r8*rxt(k,289)*y(k,92) + .800_r8*rxt(k,287)*y(k,93)
         mat(k,222) = rxt(k,291)*y(k,69) + .900_r8*rxt(k,289)*y(k,89)
         mat(k,244) = .800_r8*rxt(k,287)*y(k,89)
         mat(k,219) = -(rxt(k,289)*y(k,89) + rxt(k,290)*y(k,70) + (rxt(k,291) &
                      + rxt(k,292)) * y(k,69))
         mat(k,320) = -rxt(k,289)*y(k,92)
         mat(k,496) = -rxt(k,290)*y(k,92)
         mat(k,657) = -(rxt(k,291) + rxt(k,292)) * y(k,92)
         mat(k,242) = -(rxt(k,287)*y(k,89))
         mat(k,321) = -rxt(k,287)*y(k,93)
         mat(k,343) = rxt(k,296)*y(k,98)
         mat(k,735) = rxt(k,298)*y(k,98)
         mat(k,660) = rxt(k,291)*y(k,92)
         mat(k,498) = rxt(k,295)*y(k,94)
         mat(k,220) = rxt(k,291)*y(k,69)
         mat(k,272) = rxt(k,295)*y(k,70)
         mat(k,357) = rxt(k,296)*y(k,55) + rxt(k,298)*y(k,61)
         mat(k,273) = -(rxt(k,293)*y(k,69) + (rxt(k,294) + rxt(k,295)) * y(k,70))
         mat(k,664) = -rxt(k,293)*y(k,94)
         mat(k,499) = -(rxt(k,294) + rxt(k,295)) * y(k,94)
         mat(k,332) = rxt(k,301)*y(k,100)
         mat(k,372) = rxt(k,301)*y(k,91)
         mat(k,726) = -(rxt(k,94)*y(k,41) + rxt(k,95)*y(k,103) + (rxt(k,97) + rxt(k,98) &
                      ) * y(k,70) + rxt(k,99)*y(k,71) + (rxt(k,147) + rxt(k,148) &
                      ) * y(k,56) + rxt(k,180)*y(k,8) + rxt(k,181)*y(k,9) + rxt(k,182) &
                      *y(k,11) + rxt(k,183)*y(k,12) + rxt(k,184)*y(k,13) + rxt(k,185) &
                      *y(k,14) + rxt(k,186)*y(k,15) + (rxt(k,187) + rxt(k,188) &
                      ) * y(k,49) + rxt(k,207)*y(k,10) + rxt(k,208)*y(k,24) + rxt(k,209) &
                      *y(k,42) + (rxt(k,210) + rxt(k,211)) * y(k,45) + rxt(k,216) &
                      *y(k,33) + rxt(k,217)*y(k,34) + rxt(k,230)*y(k,16) + rxt(k,231) &
                      *y(k,18) + rxt(k,232)*y(k,46) + rxt(k,233)*y(k,47) + rxt(k,234) &
                      *y(k,48) + (rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,23))
         mat(k,877) = -rxt(k,94)*y(k,95)
         mat(k,951) = -rxt(k,95)*y(k,95)
         mat(k,514) = -(rxt(k,97) + rxt(k,98)) * y(k,95)
         mat(k,535) = -rxt(k,99)*y(k,95)
         mat(k,105) = -(rxt(k,147) + rxt(k,148)) * y(k,95)
         mat(k,30) = -rxt(k,180)*y(k,95)
         mat(k,57) = -rxt(k,181)*y(k,95)
         mat(k,38) = -rxt(k,182)*y(k,95)
         mat(k,68) = -rxt(k,183)*y(k,95)
         mat(k,42) = -rxt(k,184)*y(k,95)
         mat(k,73) = -rxt(k,185)*y(k,95)
         mat(k,46) = -rxt(k,186)*y(k,95)
         mat(k,461) = -(rxt(k,187) + rxt(k,188)) * y(k,95)
         mat(k,63) = -rxt(k,207)*y(k,95)
         mat(k,163) = -rxt(k,208)*y(k,95)
         mat(k,34) = -rxt(k,209)*y(k,95)
         mat(k,259) = -(rxt(k,210) + rxt(k,211)) * y(k,95)
         mat(k,84) = -rxt(k,216)*y(k,95)
         mat(k,92) = -rxt(k,217)*y(k,95)
         mat(k,176) = -rxt(k,230)*y(k,95)
         mat(k,214) = -rxt(k,231)*y(k,95)
         mat(k,87) = -rxt(k,232)*y(k,95)
         mat(k,97) = -rxt(k,233)*y(k,95)
         mat(k,115) = -rxt(k,234)*y(k,95)
         mat(k,403) = -(rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,95)
         mat(k,514) = mat(k,514) + rxt(k,130)*y(k,91)
         mat(k,330) = .850_r8*rxt(k,288)*y(k,98)
         mat(k,340) = rxt(k,130)*y(k,70)
         mat(k,367) = .850_r8*rxt(k,288)*y(k,89)
         mat(k,77) = -(rxt(k,101)*y(k,69) + rxt(k,102)*y(k,70))
         mat(k,651) = -rxt(k,101)*y(k,96)
         mat(k,492) = -rxt(k,102)*y(k,96)
         mat(k,194) = rxt(k,103)*y(k,97)
         mat(k,651) = mat(k,651) + rxt(k,105)*y(k,97)
         mat(k,492) = mat(k,492) + rxt(k,106)*y(k,97)
         mat(k,522) = rxt(k,107)*y(k,97)
         mat(k,79) = rxt(k,103)*y(k,32) + rxt(k,105)*y(k,69) + rxt(k,106)*y(k,70) &
                      + rxt(k,107)*y(k,71)
         mat(k,80) = -(rxt(k,103)*y(k,32) + rxt(k,105)*y(k,69) + rxt(k,106)*y(k,70) &
                      + rxt(k,107)*y(k,71))
         mat(k,195) = -rxt(k,103)*y(k,97)
         mat(k,652) = -rxt(k,105)*y(k,97)
         mat(k,493) = -rxt(k,106)*y(k,97)
         mat(k,523) = -rxt(k,107)*y(k,97)
         mat(k,493) = mat(k,493) + rxt(k,97)*y(k,95)
         mat(k,702) = rxt(k,97)*y(k,70)
         mat(k,361) = -(rxt(k,288)*y(k,89) + rxt(k,296)*y(k,55) + rxt(k,298)*y(k,61))
         mat(k,326) = -rxt(k,288)*y(k,98)
         mat(k,348) = -rxt(k,296)*y(k,98)
         mat(k,739) = -rxt(k,298)*y(k,98)
         mat(k,198) = rxt(k,299)*y(k,100)
         mat(k,504) = rxt(k,290)*y(k,92) + rxt(k,294)*y(k,94) + rxt(k,302)*y(k,100) &
                      + rxt(k,306)*y(k,101)
         mat(k,224) = rxt(k,290)*y(k,70)
         mat(k,277) = rxt(k,294)*y(k,70)
         mat(k,377) = rxt(k,299)*y(k,32) + rxt(k,302)*y(k,70)
         mat(k,190) = rxt(k,306)*y(k,70)
         mat(k,858) = -(rxt(k,120)*y(k,41) + rxt(k,121)*y(k,43) + rxt(k,122)*y(k,90) &
                      + rxt(k,123)*y(k,69) + rxt(k,124)*y(k,71) + (4._r8*rxt(k,125) &
                      + 4._r8*rxt(k,126)) * y(k,99) + rxt(k,128)*y(k,52) + rxt(k,142) &
                      *y(k,63) + rxt(k,143)*y(k,55) + rxt(k,151)*y(k,62) + rxt(k,152) &
                      *y(k,51) + rxt(k,171)*y(k,29) + (rxt(k,173) + rxt(k,174) &
                      ) * y(k,28) + rxt(k,176)*y(k,49) + rxt(k,179)*y(k,54) + rxt(k,203) &
                      *y(k,5) + rxt(k,205)*y(k,45) + rxt(k,219)*y(k,16) + rxt(k,221) &
                      *y(k,18) + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,20) + rxt(k,226) &
                      *y(k,24) + rxt(k,227)*y(k,46) + rxt(k,228)*y(k,47) + rxt(k,229) &
                      *y(k,48) + rxt(k,237)*y(k,17) + rxt(k,240)*y(k,22) + rxt(k,241) &
                      *y(k,23) + rxt(k,245)*y(k,31) + (rxt(k,247) + rxt(k,261) &
                      ) * y(k,35) + rxt(k,249)*y(k,74) + rxt(k,251)*y(k,80) + rxt(k,255) &
                      *y(k,77) + rxt(k,260)*y(k,79))
         mat(k,881) = -rxt(k,120)*y(k,99)
         mat(k,233) = -rxt(k,121)*y(k,99)
         mat(k,564) = -rxt(k,122)*y(k,99)
         mat(k,688) = -rxt(k,123)*y(k,99)
         mat(k,539) = -rxt(k,124)*y(k,99)
         mat(k,172) = -rxt(k,128)*y(k,99)
         mat(k,589) = -rxt(k,142)*y(k,99)
         mat(k,356) = -rxt(k,143)*y(k,99)
         mat(k,646) = -rxt(k,151)*y(k,99)
         mat(k,904) = -rxt(k,152)*y(k,99)
         mat(k,315) = -rxt(k,171)*y(k,99)
         mat(k,617) = -(rxt(k,173) + rxt(k,174)) * y(k,99)
         mat(k,464) = -rxt(k,176)*y(k,99)
         mat(k,270) = -rxt(k,179)*y(k,99)
         mat(k,928) = -rxt(k,203)*y(k,99)
         mat(k,261) = -rxt(k,205)*y(k,99)
         mat(k,179) = -rxt(k,219)*y(k,99)
         mat(k,217) = -rxt(k,221)*y(k,99)
         mat(k,50) = -rxt(k,222)*y(k,99)
         mat(k,153) = -rxt(k,224)*y(k,99)
         mat(k,166) = -rxt(k,226)*y(k,99)
         mat(k,89) = -rxt(k,227)*y(k,99)
         mat(k,99) = -rxt(k,228)*y(k,99)
         mat(k,117) = -rxt(k,229)*y(k,99)
         mat(k,487) = -rxt(k,237)*y(k,99)
         mat(k,159) = -rxt(k,240)*y(k,99)
         mat(k,406) = -rxt(k,241)*y(k,99)
         mat(k,288) = -rxt(k,245)*y(k,99)
         mat(k,132) = -(rxt(k,247) + rxt(k,261)) * y(k,99)
         mat(k,146) = -rxt(k,249)*y(k,99)
         mat(k,304) = -rxt(k,251)*y(k,99)
         mat(k,241) = -rxt(k,255)*y(k,99)
         mat(k,436) = -rxt(k,260)*y(k,99)
         mat(k,487) = mat(k,487) + rxt(k,236)*y(k,69)
         mat(k,159) = mat(k,159) + .300_r8*rxt(k,240)*y(k,99)
         mat(k,406) = mat(k,406) + rxt(k,242)*y(k,95)
         mat(k,787) = rxt(k,160)*y(k,90)
         mat(k,295) = rxt(k,214)*y(k,103)
         mat(k,444) = rxt(k,119)*y(k,71) + 2.000_r8*rxt(k,114)*y(k,90)
         mat(k,881) = mat(k,881) + rxt(k,111)*y(k,69) + rxt(k,94)*y(k,95)
         mat(k,233) = mat(k,233) + rxt(k,112)*y(k,69)
         mat(k,261) = mat(k,261) + rxt(k,204)*y(k,69) + rxt(k,210)*y(k,95)
         mat(k,464) = mat(k,464) + rxt(k,175)*y(k,69) + rxt(k,187)*y(k,95)
         mat(k,253) = rxt(k,206)*y(k,69)
         mat(k,270) = mat(k,270) + rxt(k,178)*y(k,69)
         mat(k,755) = rxt(k,144)*y(k,90)
         mat(k,589) = mat(k,589) + rxt(k,139)*y(k,90)
         mat(k,688) = mat(k,688) + rxt(k,236)*y(k,17) + rxt(k,111)*y(k,41) &
                      + rxt(k,112)*y(k,43) + rxt(k,204)*y(k,45) + rxt(k,175)*y(k,49) &
                      + rxt(k,206)*y(k,53) + rxt(k,178)*y(k,54) + rxt(k,117)*y(k,90)
         mat(k,539) = mat(k,539) + rxt(k,119)*y(k,40) + rxt(k,118)*y(k,90)
         mat(k,564) = mat(k,564) + rxt(k,160)*y(k,25) + 2.000_r8*rxt(k,114)*y(k,40) &
                      + rxt(k,144)*y(k,61) + rxt(k,139)*y(k,63) + rxt(k,117)*y(k,69) &
                      + rxt(k,118)*y(k,71)
         mat(k,730) = rxt(k,242)*y(k,23) + rxt(k,94)*y(k,41) + rxt(k,210)*y(k,45) &
                      + rxt(k,187)*y(k,49) + 2.000_r8*rxt(k,95)*y(k,103)
         mat(k,858) = mat(k,858) + .300_r8*rxt(k,240)*y(k,22)
         mat(k,955) = rxt(k,214)*y(k,39) + 2.000_r8*rxt(k,95)*y(k,95)
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
         mat(k,378) = -(rxt(k,299)*y(k,32) + rxt(k,301)*y(k,91) + rxt(k,302)*y(k,70))
         mat(k,199) = -rxt(k,299)*y(k,100)
         mat(k,337) = -rxt(k,301)*y(k,100)
         mat(k,505) = -rxt(k,302)*y(k,100)
         mat(k,672) = rxt(k,292)*y(k,92) + rxt(k,293)*y(k,94) + rxt(k,305)*y(k,101) &
                      + rxt(k,311)*y(k,102)
         mat(k,327) = rxt(k,303)*y(k,101) + rxt(k,308)*y(k,102)
         mat(k,225) = rxt(k,292)*y(k,69)
         mat(k,278) = rxt(k,293)*y(k,69)
         mat(k,191) = rxt(k,305)*y(k,69) + rxt(k,303)*y(k,89)
         mat(k,185) = rxt(k,311)*y(k,69) + rxt(k,308)*y(k,89)
         mat(k,188) = -(rxt(k,303)*y(k,89) + rxt(k,305)*y(k,69) + rxt(k,306)*y(k,70))
         mat(k,319) = -rxt(k,303)*y(k,101)
         mat(k,655) = -rxt(k,305)*y(k,101)
         mat(k,495) = -rxt(k,306)*y(k,101)
         mat(k,319) = mat(k,319) + rxt(k,307)*y(k,102)
         mat(k,182) = rxt(k,307)*y(k,89)
         mat(k,181) = -((rxt(k,307) + rxt(k,308)) * y(k,89) + rxt(k,311)*y(k,69))
         mat(k,318) = -(rxt(k,307) + rxt(k,308)) * y(k,102)
         mat(k,654) = -rxt(k,311)*y(k,102)
         mat(k,959) = -(rxt(k,95)*y(k,95) + rxt(k,214)*y(k,39) + rxt(k,262)*y(k,81))
         mat(k,734) = -rxt(k,95)*y(k,103)
         mat(k,298) = -rxt(k,214)*y(k,103)
         mat(k,111) = -rxt(k,262)*y(k,103)
         mat(k,180) = rxt(k,219)*y(k,99)
         mat(k,491) = rxt(k,237)*y(k,99)
         mat(k,218) = rxt(k,221)*y(k,99)
         mat(k,51) = rxt(k,222)*y(k,99)
         mat(k,154) = rxt(k,224)*y(k,99)
         mat(k,160) = rxt(k,240)*y(k,99)
         mat(k,409) = rxt(k,241)*y(k,99)
         mat(k,446) = rxt(k,115)*y(k,90)
         mat(k,885) = rxt(k,120)*y(k,99)
         mat(k,234) = rxt(k,121)*y(k,99)
         mat(k,263) = rxt(k,205)*y(k,99)
         mat(k,118) = rxt(k,229)*y(k,99)
         mat(k,468) = (rxt(k,275)+rxt(k,280))*y(k,53) + (rxt(k,268)+rxt(k,274) &
                       +rxt(k,279))*y(k,54) + rxt(k,176)*y(k,99)
         mat(k,908) = rxt(k,152)*y(k,99)
         mat(k,173) = rxt(k,128)*y(k,99)
         mat(k,255) = (rxt(k,275)+rxt(k,280))*y(k,49)
         mat(k,271) = (rxt(k,268)+rxt(k,274)+rxt(k,279))*y(k,49) + rxt(k,179)*y(k,99)
         mat(k,568) = rxt(k,115)*y(k,40) + rxt(k,122)*y(k,99)
         mat(k,862) = rxt(k,219)*y(k,16) + rxt(k,237)*y(k,17) + rxt(k,221)*y(k,18) &
                      + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,20) + rxt(k,240)*y(k,22) &
                      + rxt(k,241)*y(k,23) + rxt(k,120)*y(k,41) + rxt(k,121)*y(k,43) &
                      + rxt(k,205)*y(k,45) + rxt(k,229)*y(k,48) + rxt(k,176)*y(k,49) &
                      + rxt(k,152)*y(k,51) + rxt(k,128)*y(k,52) + rxt(k,179)*y(k,54) &
                      + rxt(k,122)*y(k,90) + 2.000_r8*rxt(k,125)*y(k,99)
      end do
      end subroutine nlnmat05
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
         mat(k, 29) = mat(k, 29) + lmat(k, 29)
         mat(k, 31) = mat(k, 31) + lmat(k, 31)
         mat(k, 32) = mat(k, 32) + lmat(k, 32)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 35) = mat(k, 35) + lmat(k, 35)
         mat(k, 36) = mat(k, 36) + lmat(k, 36)
         mat(k, 37) = mat(k, 37) + lmat(k, 37)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 40) = mat(k, 40) + lmat(k, 40)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 44) = mat(k, 44) + lmat(k, 44)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 47) = mat(k, 47) + lmat(k, 47)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = lmat(k, 53)
         mat(k, 54) = lmat(k, 54)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 66) = mat(k, 66) + lmat(k, 66)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = lmat(k, 81)
         mat(k, 82) = mat(k, 82) + lmat(k, 82)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 85) = mat(k, 85) + lmat(k, 85)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 100) = lmat(k, 100)
         mat(k, 101) = lmat(k, 101)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 109) = lmat(k, 109)
         mat(k, 110) = lmat(k, 110)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 119) = lmat(k, 119)
         mat(k, 120) = lmat(k, 120)
         mat(k, 121) = lmat(k, 121)
         mat(k, 122) = lmat(k, 122)
         mat(k, 123) = lmat(k, 123)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 138) = lmat(k, 138)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 149) = lmat(k, 149)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 157) = lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 165) = mat(k, 165) + lmat(k, 165)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 169) = lmat(k, 169)
         mat(k, 170) = lmat(k, 170)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = lmat(k, 186)
         mat(k, 187) = lmat(k, 187)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 189) = lmat(k, 189)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = mat(k, 202) + lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 207) = lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 211) = lmat(k, 211)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 252) = lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 273) = mat(k, 273) + lmat(k, 273)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 290) = mat(k, 290) + lmat(k, 290)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 301) = lmat(k, 301)
         mat(k, 303) = lmat(k, 303)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 312) = lmat(k, 312)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 323) = mat(k, 323) + lmat(k, 323)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = lmat(k, 345)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 357) = mat(k, 357) + lmat(k, 357)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 371) = lmat(k, 371)
         mat(k, 376) = lmat(k, 376)
         mat(k, 378) = mat(k, 378) + lmat(k, 378)
         mat(k, 388) = lmat(k, 388)
         mat(k, 389) = lmat(k, 389)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 402) = lmat(k, 402)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 423) = lmat(k, 423)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 452) = mat(k, 452) + lmat(k, 452)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 462) = mat(k, 462) + lmat(k, 462)
         mat(k, 470) = mat(k, 470) + lmat(k, 470)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 475) = mat(k, 475) + lmat(k, 475)
         mat(k, 488) = lmat(k, 488)
         mat(k, 494) = lmat(k, 494)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 501) = lmat(k, 501)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 528) = mat(k, 528) + lmat(k, 528)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 545) = mat(k, 545) + lmat(k, 545)
         mat(k, 555) = mat(k, 555) + lmat(k, 555)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 583) = mat(k, 583) + lmat(k, 583)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 586) = mat(k, 586) + lmat(k, 586)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 612) = mat(k, 612) + lmat(k, 612)
         mat(k, 615) = mat(k, 615) + lmat(k, 615)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 641) = mat(k, 641) + lmat(k, 641)
         mat(k, 643) = mat(k, 643) + lmat(k, 643)
         mat(k, 646) = mat(k, 646) + lmat(k, 646)
         mat(k, 648) = mat(k, 648) + lmat(k, 648)
         mat(k, 654) = mat(k, 654) + lmat(k, 654)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 668) = lmat(k, 668)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 725) = mat(k, 725) + lmat(k, 725)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 735) = mat(k, 735) + lmat(k, 735)
         mat(k, 736) = lmat(k, 736)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 752) = mat(k, 752) + lmat(k, 752)
         mat(k, 785) = mat(k, 785) + lmat(k, 785)
         mat(k, 808) = mat(k, 808) + lmat(k, 808)
         mat(k, 858) = mat(k, 858) + lmat(k, 858)
         mat(k, 882) = mat(k, 882) + lmat(k, 882)
         mat(k, 898) = lmat(k, 898)
         mat(k, 904) = mat(k, 904) + lmat(k, 904)
         mat(k, 906) = mat(k, 906) + lmat(k, 906)
         mat(k, 923) = mat(k, 923) + lmat(k, 923)
         mat(k, 927) = mat(k, 927) + lmat(k, 927)
         mat(k, 931) = mat(k, 931) + lmat(k, 931)
         mat(k, 941) = lmat(k, 941)
         mat(k, 950) = lmat(k, 950)
         mat(k, 951) = mat(k, 951) + lmat(k, 951)
         mat(k, 955) = mat(k, 955) + lmat(k, 955)
         mat(k, 956) = lmat(k, 956)
         mat(k, 959) = mat(k, 959) + lmat(k, 959)
         mat(k, 96) = 0._r8
         mat(k, 114) = 0._r8
         mat(k, 200) = 0._r8
         mat(k, 251) = 0._r8
         mat(k, 274) = 0._r8
         mat(k, 275) = 0._r8
         mat(k, 283) = 0._r8
         mat(k, 284) = 0._r8
         mat(k, 285) = 0._r8
         mat(k, 287) = 0._r8
         mat(k, 305) = 0._r8
         mat(k, 317) = 0._r8
         mat(k, 322) = 0._r8
         mat(k, 328) = 0._r8
         mat(k, 331) = 0._r8
         mat(k, 333) = 0._r8
         mat(k, 336) = 0._r8
         mat(k, 346) = 0._r8
         mat(k, 349) = 0._r8
         mat(k, 354) = 0._r8
         mat(k, 359) = 0._r8
         mat(k, 362) = 0._r8
         mat(k, 363) = 0._r8
         mat(k, 365) = 0._r8
         mat(k, 369) = 0._r8
         mat(k, 374) = 0._r8
         mat(k, 379) = 0._r8
         mat(k, 381) = 0._r8
         mat(k, 382) = 0._r8
         mat(k, 384) = 0._r8
         mat(k, 385) = 0._r8
         mat(k, 386) = 0._r8
         mat(k, 391) = 0._r8
         mat(k, 392) = 0._r8
         mat(k, 398) = 0._r8
         mat(k, 400) = 0._r8
         mat(k, 401) = 0._r8
         mat(k, 404) = 0._r8
         mat(k, 408) = 0._r8
         mat(k, 412) = 0._r8
         mat(k, 420) = 0._r8
         mat(k, 421) = 0._r8
         mat(k, 429) = 0._r8
         mat(k, 438) = 0._r8
         mat(k, 454) = 0._r8
         mat(k, 455) = 0._r8
         mat(k, 456) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 459) = 0._r8
         mat(k, 463) = 0._r8
         mat(k, 465) = 0._r8
         mat(k, 467) = 0._r8
         mat(k, 471) = 0._r8
         mat(k, 472) = 0._r8
         mat(k, 476) = 0._r8
         mat(k, 477) = 0._r8
         mat(k, 480) = 0._r8
         mat(k, 481) = 0._r8
         mat(k, 483) = 0._r8
         mat(k, 484) = 0._r8
         mat(k, 490) = 0._r8
         mat(k, 511) = 0._r8
         mat(k, 512) = 0._r8
         mat(k, 516) = 0._r8
         mat(k, 517) = 0._r8
         mat(k, 518) = 0._r8
         mat(k, 519) = 0._r8
         mat(k, 520) = 0._r8
         mat(k, 521) = 0._r8
         mat(k, 540) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 552) = 0._r8
         mat(k, 560) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 573) = 0._r8
         mat(k, 574) = 0._r8
         mat(k, 575) = 0._r8
         mat(k, 576) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 582) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 587) = 0._r8
         mat(k, 588) = 0._r8
         mat(k, 590) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 603) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 609) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 619) = 0._r8
         mat(k, 621) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 630) = 0._r8
         mat(k, 631) = 0._r8
         mat(k, 633) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 644) = 0._r8
         mat(k, 645) = 0._r8
         mat(k, 647) = 0._r8
         mat(k, 650) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 671) = 0._r8
         mat(k, 684) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 692) = 0._r8
         mat(k, 722) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 732) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 757) = 0._r8
         mat(k, 759) = 0._r8
         mat(k, 769) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 781) = 0._r8
         mat(k, 782) = 0._r8
         mat(k, 783) = 0._r8
         mat(k, 784) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 791) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 801) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 803) = 0._r8
         mat(k, 804) = 0._r8
         mat(k, 805) = 0._r8
         mat(k, 806) = 0._r8
         mat(k, 807) = 0._r8
         mat(k, 809) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 839) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 854) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 866) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 870) = 0._r8
         mat(k, 871) = 0._r8
         mat(k, 872) = 0._r8
         mat(k, 873) = 0._r8
         mat(k, 874) = 0._r8
         mat(k, 875) = 0._r8
         mat(k, 878) = 0._r8
         mat(k, 880) = 0._r8
         mat(k, 883) = 0._r8
         mat(k, 884) = 0._r8
         mat(k, 888) = 0._r8
         mat(k, 889) = 0._r8
         mat(k, 890) = 0._r8
         mat(k, 891) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 895) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 900) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 903) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 916) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 920) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 929) = 0._r8
         mat(k, 930) = 0._r8
         mat(k, 932) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 940) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 943) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 948) = 0._r8
         mat(k, 949) = 0._r8
         mat(k, 952) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 958) = 0._r8
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
         mat(k, 29) = mat(k, 29) - dti(k)
         mat(k, 32) = mat(k, 32) - dti(k)
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 40) = mat(k, 40) - dti(k)
         mat(k, 44) = mat(k, 44) - dti(k)
         mat(k, 48) = mat(k, 48) - dti(k)
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 60) = mat(k, 60) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 70) = mat(k, 70) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 80) = mat(k, 80) - dti(k)
         mat(k, 82) = mat(k, 82) - dti(k)
         mat(k, 85) = mat(k, 85) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 161) = mat(k, 161) - dti(k)
         mat(k, 167) = mat(k, 167) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 181) = mat(k, 181) - dti(k)
         mat(k, 188) = mat(k, 188) - dti(k)
         mat(k, 196) = mat(k, 196) - dti(k)
         mat(k, 202) = mat(k, 202) - dti(k)
         mat(k, 210) = mat(k, 210) - dti(k)
         mat(k, 219) = mat(k, 219) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 235) = mat(k, 235) - dti(k)
         mat(k, 242) = mat(k, 242) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 256) = mat(k, 256) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 290) = mat(k, 290) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 308) = mat(k, 308) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 361) = mat(k, 361) - dti(k)
         mat(k, 378) = mat(k, 378) - dti(k)
         mat(k, 393) = mat(k, 393) - dti(k)
         mat(k, 411) = mat(k, 411) - dti(k)
         mat(k, 425) = mat(k, 425) - dti(k)
         mat(k, 439) = mat(k, 439) - dti(k)
         mat(k, 453) = mat(k, 453) - dti(k)
         mat(k, 475) = mat(k, 475) - dti(k)
         mat(k, 508) = mat(k, 508) - dti(k)
         mat(k, 529) = mat(k, 529) - dti(k)
         mat(k, 555) = mat(k, 555) - dti(k)
         mat(k, 581) = mat(k, 581) - dti(k)
         mat(k, 610) = mat(k, 610) - dti(k)
         mat(k, 640) = mat(k, 640) - dti(k)
         mat(k, 683) = mat(k, 683) - dti(k)
         mat(k, 726) = mat(k, 726) - dti(k)
         mat(k, 752) = mat(k, 752) - dti(k)
         mat(k, 785) = mat(k, 785) - dti(k)
         mat(k, 808) = mat(k, 808) - dti(k)
         mat(k, 858) = mat(k, 858) - dti(k)
         mat(k, 882) = mat(k, 882) - dti(k)
         mat(k, 906) = mat(k, 906) - dti(k)
         mat(k, 931) = mat(k, 931) - dti(k)
         mat(k, 959) = mat(k, 959) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
