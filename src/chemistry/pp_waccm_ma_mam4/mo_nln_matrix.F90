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
         mat(k,767) = -(rxt(k,191)*y(k,17) + rxt(k,192)*y(k,87) + rxt(k,193)*y(k,70))
         mat(k,482) = -rxt(k,191)*y(k,3)
         mat(k,559) = -rxt(k,192)*y(k,3)
         mat(k,534) = -rxt(k,193)*y(k,3)
         mat(k,745) = 4.000_r8*rxt(k,194)*y(k,5) + (rxt(k,195)+rxt(k,196))*y(k,28) &
                      + rxt(k,199)*y(k,61) + rxt(k,202)*y(k,68) + rxt(k,253)*y(k,77) &
                      + rxt(k,203)*y(k,96)
         mat(k,56) = rxt(k,181)*y(k,92)
         mat(k,62) = rxt(k,207)*y(k,92)
         mat(k,167) = 2.000_r8*rxt(k,218)*y(k,25) + 2.000_r8*rxt(k,230)*y(k,92) &
                      + 2.000_r8*rxt(k,219)*y(k,96)
         mat(k,212) = rxt(k,220)*y(k,25) + rxt(k,231)*y(k,92) + rxt(k,221)*y(k,96)
         mat(k,155) = 3.000_r8*rxt(k,225)*y(k,25) + 3.000_r8*rxt(k,208)*y(k,92) &
                      + 3.000_r8*rxt(k,226)*y(k,96)
         mat(k,900) = 2.000_r8*rxt(k,218)*y(k,16) + rxt(k,220)*y(k,18) &
                      + 3.000_r8*rxt(k,225)*y(k,24)
         mat(k,587) = (rxt(k,195)+rxt(k,196))*y(k,5)
         mat(k,33) = 2.000_r8*rxt(k,209)*y(k,92)
         mat(k,258) = rxt(k,204)*y(k,68) + rxt(k,210)*y(k,92) + rxt(k,205)*y(k,96)
         mat(k,679) = rxt(k,199)*y(k,5)
         mat(k,654) = rxt(k,202)*y(k,5) + rxt(k,204)*y(k,45)
         mat(k,431) = rxt(k,253)*y(k,5)
         mat(k,721) = rxt(k,181)*y(k,9) + rxt(k,207)*y(k,10) + 2.000_r8*rxt(k,230) &
                      *y(k,16) + rxt(k,231)*y(k,18) + 3.000_r8*rxt(k,208)*y(k,24) &
                      + 2.000_r8*rxt(k,209)*y(k,42) + rxt(k,210)*y(k,45)
         mat(k,845) = rxt(k,203)*y(k,5) + 2.000_r8*rxt(k,219)*y(k,16) + rxt(k,221) &
                      *y(k,18) + 3.000_r8*rxt(k,226)*y(k,24) + rxt(k,205)*y(k,45)
         mat(k,728) = rxt(k,197)*y(k,28)
         mat(k,568) = rxt(k,197)*y(k,5)
         mat(k,445) = (rxt(k,275)+rxt(k,280))*y(k,53)
         mat(k,244) = (rxt(k,275)+rxt(k,280))*y(k,49)
         mat(k,744) = -(4._r8*rxt(k,194)*y(k,5) + (rxt(k,195) + rxt(k,196) + rxt(k,197) &
                      ) * y(k,28) + rxt(k,198)*y(k,87) + rxt(k,199)*y(k,61) + rxt(k,200) &
                      *y(k,62) + rxt(k,202)*y(k,68) + rxt(k,203)*y(k,96) + rxt(k,253) &
                      *y(k,77))
         mat(k,586) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,5)
         mat(k,558) = -rxt(k,198)*y(k,5)
         mat(k,678) = -rxt(k,199)*y(k,5)
         mat(k,795) = -rxt(k,200)*y(k,5)
         mat(k,653) = -rxt(k,202)*y(k,5)
         mat(k,844) = -rxt(k,203)*y(k,5)
         mat(k,430) = -rxt(k,253)*y(k,5)
         mat(k,766) = rxt(k,193)*y(k,70)
         mat(k,203) = rxt(k,201)*y(k,68)
         mat(k,257) = rxt(k,211)*y(k,92)
         mat(k,248) = rxt(k,206)*y(k,68)
         mat(k,653) = mat(k,653) + rxt(k,201)*y(k,6) + rxt(k,206)*y(k,53)
         mat(k,533) = rxt(k,193)*y(k,3)
         mat(k,720) = rxt(k,211)*y(k,45)
         mat(k,199) = -(rxt(k,201)*y(k,68))
         mat(k,624) = -rxt(k,201)*y(k,6)
         mat(k,730) = rxt(k,200)*y(k,62)
         mat(k,777) = rxt(k,200)*y(k,5)
         mat(k,27) = -(rxt(k,180)*y(k,92))
         mat(k,686) = -rxt(k,180)*y(k,8)
         mat(k,53) = -(rxt(k,181)*y(k,92))
         mat(k,691) = -rxt(k,181)*y(k,9)
         mat(k,58) = -(rxt(k,207)*y(k,92))
         mat(k,692) = -rxt(k,207)*y(k,10)
         mat(k,34) = -(rxt(k,182)*y(k,92))
         mat(k,688) = -rxt(k,182)*y(k,11)
         mat(k,63) = -(rxt(k,183)*y(k,92))
         mat(k,693) = -rxt(k,183)*y(k,12)
         mat(k,38) = -(rxt(k,184)*y(k,92))
         mat(k,689) = -rxt(k,184)*y(k,13)
         mat(k,68) = -(rxt(k,185)*y(k,92))
         mat(k,694) = -rxt(k,185)*y(k,14)
         mat(k,42) = -(rxt(k,186)*y(k,92))
         mat(k,690) = -rxt(k,186)*y(k,15)
         mat(k,164) = -(rxt(k,218)*y(k,25) + rxt(k,219)*y(k,96) + rxt(k,230)*y(k,92))
         mat(k,878) = -rxt(k,218)*y(k,16)
         mat(k,815) = -rxt(k,219)*y(k,16)
         mat(k,703) = -rxt(k,230)*y(k,16)
         mat(k,472) = -(rxt(k,155)*y(k,25) + rxt(k,191)*y(k,3) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,68) + rxt(k,237)*y(k,96))
         mat(k,890) = -rxt(k,155)*y(k,17)
         mat(k,757) = -rxt(k,191)*y(k,17)
         mat(k,602) = -rxt(k,235)*y(k,17)
         mat(k,644) = -rxt(k,236)*y(k,17)
         mat(k,835) = -rxt(k,237)*y(k,17)
         mat(k,410) = rxt(k,162)*y(k,28) + rxt(k,239)*y(k,61)
         mat(k,161) = .300_r8*rxt(k,240)*y(k,96)
         mat(k,394) = (rxt(k,243)+rxt(k,244))*y(k,92)
         mat(k,577) = rxt(k,162)*y(k,21)
         mat(k,669) = rxt(k,239)*y(k,21)
         mat(k,711) = (rxt(k,243)+rxt(k,244))*y(k,23)
         mat(k,835) = mat(k,835) + .300_r8*rxt(k,240)*y(k,22)
         mat(k,207) = -(rxt(k,220)*y(k,25) + rxt(k,221)*y(k,96) + rxt(k,231)*y(k,92))
         mat(k,879) = -rxt(k,220)*y(k,18)
         mat(k,818) = -rxt(k,221)*y(k,18)
         mat(k,704) = -rxt(k,231)*y(k,18)
         mat(k,46) = -(rxt(k,222)*y(k,96))
         mat(k,803) = -rxt(k,222)*y(k,19)
         mat(k,144) = -(rxt(k,223)*y(k,25) + rxt(k,224)*y(k,96))
         mat(k,876) = -rxt(k,223)*y(k,20)
         mat(k,812) = -rxt(k,224)*y(k,20)
         mat(k,408) = -(rxt(k,162)*y(k,28) + rxt(k,238)*y(k,87) + rxt(k,239)*y(k,61))
         mat(k,573) = -rxt(k,162)*y(k,21)
         mat(k,546) = -rxt(k,238)*y(k,21)
         mat(k,667) = -rxt(k,239)*y(k,21)
         mat(k,159) = .700_r8*rxt(k,240)*y(k,96)
         mat(k,391) = rxt(k,156)*y(k,25) + rxt(k,212)*y(k,39) + rxt(k,242)*y(k,92) &
                      + rxt(k,241)*y(k,96)
         mat(k,887) = rxt(k,156)*y(k,23)
         mat(k,289) = rxt(k,212)*y(k,23)
         mat(k,708) = rxt(k,242)*y(k,23)
         mat(k,831) = .700_r8*rxt(k,240)*y(k,22) + rxt(k,241)*y(k,23)
         mat(k,158) = -(rxt(k,240)*y(k,96))
         mat(k,814) = -rxt(k,240)*y(k,22)
         mat(k,407) = rxt(k,238)*y(k,87)
         mat(k,540) = rxt(k,238)*y(k,21)
         mat(k,390) = -(rxt(k,156)*y(k,25) + rxt(k,212)*y(k,39) + rxt(k,241)*y(k,96) &
                      + (rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,92))
         mat(k,886) = -rxt(k,156)*y(k,23)
         mat(k,288) = -rxt(k,212)*y(k,23)
         mat(k,830) = -rxt(k,241)*y(k,23)
         mat(k,707) = -(rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,23)
         mat(k,152) = -(rxt(k,208)*y(k,92) + rxt(k,225)*y(k,25) + rxt(k,226)*y(k,96))
         mat(k,702) = -rxt(k,208)*y(k,24)
         mat(k,877) = -rxt(k,225)*y(k,24)
         mat(k,813) = -rxt(k,226)*y(k,24)
         mat(k,904) = -(rxt(k,155)*y(k,17) + rxt(k,156)*y(k,23) + rxt(k,157)*y(k,41) &
                      + rxt(k,158)*y(k,43) + (rxt(k,159) + rxt(k,160)) * y(k,87) &
                      + rxt(k,161)*y(k,70) + rxt(k,168)*y(k,29) + rxt(k,177)*y(k,54) &
                      + rxt(k,218)*y(k,16) + rxt(k,220)*y(k,18) + rxt(k,223)*y(k,20) &
                      + rxt(k,225)*y(k,24))
         mat(k,486) = -rxt(k,155)*y(k,25)
         mat(k,404) = -rxt(k,156)*y(k,25)
         mat(k,927) = -rxt(k,157)*y(k,25)
         mat(k,221) = -rxt(k,158)*y(k,25)
         mat(k,563) = -(rxt(k,159) + rxt(k,160)) * y(k,25)
         mat(k,537) = -rxt(k,161)*y(k,25)
         mat(k,313) = -rxt(k,168)*y(k,25)
         mat(k,267) = -rxt(k,177)*y(k,25)
         mat(k,169) = -rxt(k,218)*y(k,25)
         mat(k,214) = -rxt(k,220)*y(k,25)
         mat(k,150) = -rxt(k,223)*y(k,25)
         mat(k,157) = -rxt(k,225)*y(k,25)
         mat(k,749) = rxt(k,196)*y(k,28)
         mat(k,29) = 4.000_r8*rxt(k,180)*y(k,92)
         mat(k,57) = rxt(k,181)*y(k,92)
         mat(k,37) = 2.000_r8*rxt(k,182)*y(k,92)
         mat(k,67) = 2.000_r8*rxt(k,183)*y(k,92)
         mat(k,41) = 2.000_r8*rxt(k,184)*y(k,92)
         mat(k,72) = rxt(k,185)*y(k,92)
         mat(k,45) = 2.000_r8*rxt(k,186)*y(k,92)
         mat(k,48) = 3.000_r8*rxt(k,222)*y(k,96)
         mat(k,150) = mat(k,150) + rxt(k,224)*y(k,96)
         mat(k,417) = rxt(k,162)*y(k,28)
         mat(k,591) = rxt(k,196)*y(k,5) + rxt(k,162)*y(k,21) + (4.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,165))*y(k,28) + rxt(k,167)*y(k,61) + rxt(k,172) &
                      *y(k,68) + rxt(k,254)*y(k,77) + rxt(k,173)*y(k,96)
         mat(k,88) = rxt(k,217)*y(k,92)
         mat(k,84) = rxt(k,232)*y(k,92) + rxt(k,227)*y(k,96)
         mat(k,93) = rxt(k,233)*y(k,92) + rxt(k,228)*y(k,96)
         mat(k,114) = rxt(k,234)*y(k,92) + rxt(k,229)*y(k,96)
         mat(k,463) = rxt(k,175)*y(k,68) + rxt(k,187)*y(k,92) + rxt(k,176)*y(k,96)
         mat(k,683) = rxt(k,167)*y(k,28)
         mat(k,658) = rxt(k,172)*y(k,28) + rxt(k,175)*y(k,49)
         mat(k,434) = rxt(k,254)*y(k,28)
         mat(k,725) = 4.000_r8*rxt(k,180)*y(k,8) + rxt(k,181)*y(k,9) &
                      + 2.000_r8*rxt(k,182)*y(k,11) + 2.000_r8*rxt(k,183)*y(k,12) &
                      + 2.000_r8*rxt(k,184)*y(k,13) + rxt(k,185)*y(k,14) &
                      + 2.000_r8*rxt(k,186)*y(k,15) + rxt(k,217)*y(k,34) + rxt(k,232) &
                      *y(k,46) + rxt(k,233)*y(k,47) + rxt(k,234)*y(k,48) + rxt(k,187) &
                      *y(k,49)
         mat(k,849) = 3.000_r8*rxt(k,222)*y(k,19) + rxt(k,224)*y(k,20) + rxt(k,173) &
                      *y(k,28) + rxt(k,227)*y(k,46) + rxt(k,228)*y(k,47) + rxt(k,229) &
                      *y(k,48) + rxt(k,176)*y(k,49)
         mat(k,875) = rxt(k,168)*y(k,29)
         mat(k,567) = 2.000_r8*rxt(k,164)*y(k,28)
         mat(k,303) = rxt(k,168)*y(k,25) + (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,49)
         mat(k,444) = (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,29) + (rxt(k,268) &
                       +rxt(k,274)+rxt(k,279))*y(k,54)
         mat(k,261) = (rxt(k,268)+rxt(k,274)+rxt(k,279))*y(k,49)
         mat(k,566) = 2.000_r8*rxt(k,189)*y(k,28)
         mat(k,581) = -(rxt(k,162)*y(k,21) + (4._r8*rxt(k,163) + 4._r8*rxt(k,164) &
                      + 4._r8*rxt(k,165) + 4._r8*rxt(k,189)) * y(k,28) + rxt(k,166) &
                      *y(k,87) + rxt(k,167)*y(k,61) + rxt(k,169)*y(k,62) + rxt(k,172) &
                      *y(k,68) + (rxt(k,173) + rxt(k,174)) * y(k,96) + (rxt(k,195) &
                      + rxt(k,196) + rxt(k,197)) * y(k,5) + rxt(k,254)*y(k,77))
         mat(k,413) = -rxt(k,162)*y(k,28)
         mat(k,553) = -rxt(k,166)*y(k,28)
         mat(k,673) = -rxt(k,167)*y(k,28)
         mat(k,790) = -rxt(k,169)*y(k,28)
         mat(k,648) = -rxt(k,172)*y(k,28)
         mat(k,839) = -(rxt(k,173) + rxt(k,174)) * y(k,28)
         mat(k,739) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,28)
         mat(k,427) = -rxt(k,254)*y(k,28)
         mat(k,894) = rxt(k,177)*y(k,54) + rxt(k,161)*y(k,70) + rxt(k,160)*y(k,87)
         mat(k,307) = rxt(k,170)*y(k,68)
         mat(k,454) = rxt(k,188)*y(k,92)
         mat(k,264) = rxt(k,177)*y(k,25) + rxt(k,178)*y(k,68) + rxt(k,179)*y(k,96)
         mat(k,648) = mat(k,648) + rxt(k,170)*y(k,29) + rxt(k,178)*y(k,54)
         mat(k,528) = rxt(k,161)*y(k,25)
         mat(k,134) = rxt(k,259)*y(k,77)
         mat(k,427) = mat(k,427) + rxt(k,259)*y(k,71)
         mat(k,553) = mat(k,553) + rxt(k,160)*y(k,25)
         mat(k,715) = rxt(k,188)*y(k,49)
         mat(k,839) = mat(k,839) + rxt(k,179)*y(k,54)
         mat(k,305) = -(rxt(k,168)*y(k,25) + rxt(k,170)*y(k,68) + rxt(k,171)*y(k,96) &
                      + (rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,49))
         mat(k,883) = -rxt(k,168)*y(k,29)
         mat(k,635) = -rxt(k,170)*y(k,29)
         mat(k,826) = -rxt(k,171)*y(k,29)
         mat(k,448) = -(rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,29)
         mat(k,572) = rxt(k,169)*y(k,62)
         mat(k,780) = rxt(k,169)*y(k,28)
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
         mat(k,279) = -(rxt(k,245)*y(k,96))
         mat(k,823) = -rxt(k,245)*y(k,31)
         mat(k,753) = rxt(k,191)*y(k,17)
         mat(k,467) = rxt(k,191)*y(k,3) + rxt(k,155)*y(k,25) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,68) + rxt(k,237)*y(k,96)
         mat(k,145) = rxt(k,223)*y(k,25)
         mat(k,882) = rxt(k,155)*y(k,17) + rxt(k,223)*y(k,20)
         mat(k,194) = rxt(k,299)*y(k,97)
         mat(k,596) = rxt(k,235)*y(k,17)
         mat(k,633) = rxt(k,236)*y(k,17) + rxt(k,248)*y(k,72)
         mat(k,138) = rxt(k,248)*y(k,68) + rxt(k,249)*y(k,96)
         mat(k,823) = mat(k,823) + rxt(k,237)*y(k,17) + rxt(k,249)*y(k,72)
         mat(k,370) = rxt(k,299)*y(k,32)
         mat(k,193) = -(rxt(k,299)*y(k,97))
         mat(k,367) = -rxt(k,299)*y(k,32)
         mat(k,278) = rxt(k,245)*y(k,96)
         mat(k,817) = rxt(k,245)*y(k,31)
         mat(k,94) = -(rxt(k,216)*y(k,92))
         mat(k,699) = -rxt(k,216)*y(k,33)
         mat(k,54) = rxt(k,181)*y(k,92)
         mat(k,59) = rxt(k,207)*y(k,92)
         mat(k,65) = rxt(k,183)*y(k,92)
         mat(k,39) = 2.000_r8*rxt(k,184)*y(k,92)
         mat(k,69) = 2.000_r8*rxt(k,185)*y(k,92)
         mat(k,43) = rxt(k,186)*y(k,92)
         mat(k,31) = 2.000_r8*rxt(k,209)*y(k,92)
         mat(k,90) = rxt(k,233)*y(k,92) + rxt(k,228)*y(k,96)
         mat(k,109) = rxt(k,234)*y(k,92) + rxt(k,229)*y(k,96)
         mat(k,699) = mat(k,699) + rxt(k,181)*y(k,9) + rxt(k,207)*y(k,10) + rxt(k,183) &
                      *y(k,12) + 2.000_r8*rxt(k,184)*y(k,13) + 2.000_r8*rxt(k,185) &
                      *y(k,14) + rxt(k,186)*y(k,15) + 2.000_r8*rxt(k,209)*y(k,42) &
                      + rxt(k,233)*y(k,47) + rxt(k,234)*y(k,48)
         mat(k,807) = rxt(k,228)*y(k,47) + rxt(k,229)*y(k,48)
         mat(k,85) = -(rxt(k,217)*y(k,92))
         mat(k,697) = -rxt(k,217)*y(k,34)
         mat(k,35) = rxt(k,182)*y(k,92)
         mat(k,64) = rxt(k,183)*y(k,92)
         mat(k,81) = rxt(k,232)*y(k,92) + rxt(k,227)*y(k,96)
         mat(k,697) = mat(k,697) + rxt(k,182)*y(k,11) + rxt(k,183)*y(k,12) &
                      + rxt(k,232)*y(k,46)
         mat(k,805) = rxt(k,227)*y(k,46)
         mat(k,125) = -(rxt(k,246)*y(k,63) + (rxt(k,247) + rxt(k,261)) * y(k,96))
         mat(k,595) = -rxt(k,246)*y(k,35)
         mat(k,810) = -(rxt(k,247) + rxt(k,261)) * y(k,35)
         mat(k,287) = -(rxt(k,212)*y(k,23) + rxt(k,213)*y(k,41) + rxt(k,214)*y(k,100) &
                      + rxt(k,215)*y(k,51))
         mat(k,387) = -rxt(k,212)*y(k,39)
         mat(k,908) = -rxt(k,213)*y(k,39)
         mat(k,933) = -rxt(k,214)*y(k,39)
         mat(k,853) = -rxt(k,215)*y(k,39)
         mat(k,60) = rxt(k,207)*y(k,92)
         mat(k,70) = rxt(k,185)*y(k,92)
         mat(k,95) = 2.000_r8*rxt(k,216)*y(k,92)
         mat(k,86) = rxt(k,217)*y(k,92)
         mat(k,706) = rxt(k,207)*y(k,10) + rxt(k,185)*y(k,14) + 2.000_r8*rxt(k,216) &
                      *y(k,33) + rxt(k,217)*y(k,34)
         mat(k,436) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,87) + rxt(k,116) &
                      *y(k,69) + rxt(k,119)*y(k,70))
         mat(k,547) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,40)
         mat(k,504) = -rxt(k,116)*y(k,40)
         mat(k,524) = -rxt(k,119)*y(k,40)
         mat(k,470) = rxt(k,237)*y(k,96)
         mat(k,392) = rxt(k,243)*y(k,92)
         mat(k,888) = rxt(k,157)*y(k,41)
         mat(k,290) = rxt(k,213)*y(k,41)
         mat(k,911) = rxt(k,157)*y(k,25) + rxt(k,213)*y(k,39) + rxt(k,111)*y(k,68) &
                      + rxt(k,94)*y(k,92) + rxt(k,120)*y(k,96)
         mat(k,254) = rxt(k,211)*y(k,92)
         mat(k,449) = rxt(k,188)*y(k,92)
         mat(k,347) = rxt(k,143)*y(k,96)
         mat(k,642) = rxt(k,111)*y(k,41) + rxt(k,123)*y(k,96)
         mat(k,141) = rxt(k,249)*y(k,96)
         mat(k,234) = rxt(k,255)*y(k,96)
         mat(k,423) = rxt(k,260)*y(k,96)
         mat(k,709) = rxt(k,243)*y(k,23) + rxt(k,94)*y(k,41) + rxt(k,211)*y(k,45) &
                      + rxt(k,188)*y(k,49)
         mat(k,833) = rxt(k,237)*y(k,17) + rxt(k,120)*y(k,41) + rxt(k,143)*y(k,55) &
                      + rxt(k,123)*y(k,68) + rxt(k,249)*y(k,72) + rxt(k,255)*y(k,75) &
                      + rxt(k,260)*y(k,77)
         mat(k,928) = -(rxt(k,94)*y(k,92) + rxt(k,111)*y(k,68) + rxt(k,120)*y(k,96) &
                      + rxt(k,157)*y(k,25) + rxt(k,213)*y(k,39))
         mat(k,726) = -rxt(k,94)*y(k,41)
         mat(k,659) = -rxt(k,111)*y(k,41)
         mat(k,850) = -rxt(k,120)*y(k,41)
         mat(k,905) = -rxt(k,157)*y(k,41)
         mat(k,294) = -rxt(k,213)*y(k,41)
         mat(k,405) = rxt(k,244)*y(k,92)
         mat(k,442) = rxt(k,113)*y(k,87)
         mat(k,564) = rxt(k,113)*y(k,40)
         mat(k,726) = mat(k,726) + rxt(k,244)*y(k,23)
         mat(k,30) = -(rxt(k,209)*y(k,92))
         mat(k,687) = -rxt(k,209)*y(k,42)
         mat(k,216) = -(rxt(k,112)*y(k,68) + rxt(k,121)*y(k,96) + rxt(k,158)*y(k,25))
         mat(k,625) = -rxt(k,112)*y(k,43)
         mat(k,819) = -rxt(k,121)*y(k,43)
         mat(k,880) = -rxt(k,158)*y(k,43)
         mat(k,542) = 2.000_r8*rxt(k,127)*y(k,87)
         mat(k,819) = mat(k,819) + 2.000_r8*rxt(k,126)*y(k,96)
         mat(k,104) = rxt(k,262)*y(k,100)
         mat(k,930) = rxt(k,262)*y(k,79)
         mat(k,253) = -(rxt(k,204)*y(k,68) + rxt(k,205)*y(k,96) + (rxt(k,210) &
                      + rxt(k,211)) * y(k,92))
         mat(k,630) = -rxt(k,204)*y(k,45)
         mat(k,821) = -rxt(k,205)*y(k,45)
         mat(k,705) = -(rxt(k,210) + rxt(k,211)) * y(k,45)
         mat(k,752) = rxt(k,191)*y(k,17) + rxt(k,192)*y(k,87)
         mat(k,466) = rxt(k,191)*y(k,3)
         mat(k,544) = rxt(k,192)*y(k,3)
         mat(k,80) = -(rxt(k,227)*y(k,96) + rxt(k,232)*y(k,92))
         mat(k,804) = -rxt(k,227)*y(k,46)
         mat(k,696) = -rxt(k,232)*y(k,46)
         mat(k,89) = -(rxt(k,228)*y(k,96) + rxt(k,233)*y(k,92))
         mat(k,806) = -rxt(k,228)*y(k,47)
         mat(k,698) = -rxt(k,233)*y(k,47)
         mat(k,110) = -(rxt(k,229)*y(k,96) + rxt(k,234)*y(k,92))
         mat(k,809) = -rxt(k,229)*y(k,48)
         mat(k,701) = -rxt(k,234)*y(k,48)
         mat(k,450) = -(rxt(k,175)*y(k,68) + rxt(k,176)*y(k,96) + (rxt(k,187) &
                      + rxt(k,188)) * y(k,92) + (rxt(k,268) + rxt(k,274) + rxt(k,279) &
                      ) * y(k,54) + (rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,29) &
                      + (rxt(k,275) + rxt(k,280)) * y(k,53))
         mat(k,643) = -rxt(k,175)*y(k,49)
         mat(k,834) = -rxt(k,176)*y(k,49)
         mat(k,710) = -(rxt(k,187) + rxt(k,188)) * y(k,49)
         mat(k,263) = -(rxt(k,268) + rxt(k,274) + rxt(k,279)) * y(k,49)
         mat(k,306) = -(rxt(k,273) + rxt(k,278) + rxt(k,283)) * y(k,49)
         mat(k,246) = -(rxt(k,275) + rxt(k,280)) * y(k,49)
         mat(k,165) = rxt(k,218)*y(k,25)
         mat(k,471) = rxt(k,155)*y(k,25)
         mat(k,209) = rxt(k,220)*y(k,25)
         mat(k,147) = 2.000_r8*rxt(k,223)*y(k,25)
         mat(k,393) = rxt(k,156)*y(k,25)
         mat(k,153) = rxt(k,225)*y(k,25)
         mat(k,889) = rxt(k,218)*y(k,16) + rxt(k,155)*y(k,17) + rxt(k,220)*y(k,18) &
                      + 2.000_r8*rxt(k,223)*y(k,20) + rxt(k,156)*y(k,23) + rxt(k,225) &
                      *y(k,24) + rxt(k,157)*y(k,41) + rxt(k,158)*y(k,43) + rxt(k,177) &
                      *y(k,54) + rxt(k,159)*y(k,87)
         mat(k,576) = rxt(k,174)*y(k,96)
         mat(k,912) = rxt(k,157)*y(k,25)
         mat(k,217) = rxt(k,158)*y(k,25)
         mat(k,263) = mat(k,263) + rxt(k,177)*y(k,25)
         mat(k,548) = rxt(k,159)*y(k,25)
         mat(k,834) = mat(k,834) + rxt(k,174)*y(k,28)
         mat(k,384) = rxt(k,212)*y(k,39)
         mat(k,286) = rxt(k,212)*y(k,23) + rxt(k,213)*y(k,41) + rxt(k,215)*y(k,51) &
                      + rxt(k,214)*y(k,100)
         mat(k,907) = rxt(k,213)*y(k,39)
         mat(k,852) = rxt(k,215)*y(k,39)
         mat(k,932) = rxt(k,214)*y(k,39)
         mat(k,871) = -(rxt(k,152)*y(k,96) + rxt(k,215)*y(k,39))
         mat(k,848) = -rxt(k,152)*y(k,51)
         mat(k,293) = -rxt(k,215)*y(k,51)
         mat(k,485) = rxt(k,235)*y(k,63)
         mat(k,312) = (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,49)
         mat(k,130) = rxt(k,246)*y(k,63)
         mat(k,462) = (rxt(k,273)+rxt(k,278)+rxt(k,283))*y(k,29)
         mat(k,799) = rxt(k,151)*y(k,96)
         mat(k,615) = rxt(k,235)*y(k,17) + rxt(k,246)*y(k,35)
         mat(k,848) = mat(k,848) + rxt(k,151)*y(k,62)
         mat(k,171) = -(rxt(k,128)*y(k,96))
         mat(k,816) = -rxt(k,128)*y(k,52)
         mat(k,776) = rxt(k,149)*y(k,87)
         mat(k,541) = rxt(k,149)*y(k,62)
         mat(k,245) = -(rxt(k,206)*y(k,68) + (rxt(k,275) + rxt(k,280)) * y(k,49))
         mat(k,629) = -rxt(k,206)*y(k,53)
         mat(k,446) = -(rxt(k,275) + rxt(k,280)) * y(k,53)
         mat(k,731) = rxt(k,198)*y(k,87)
         mat(k,543) = rxt(k,198)*y(k,5)
         mat(k,262) = -(rxt(k,177)*y(k,25) + rxt(k,178)*y(k,68) + rxt(k,179)*y(k,96) &
                      + (rxt(k,268) + rxt(k,274) + rxt(k,279)) * y(k,49))
         mat(k,881) = -rxt(k,177)*y(k,54)
         mat(k,631) = -rxt(k,178)*y(k,54)
         mat(k,822) = -rxt(k,179)*y(k,54)
         mat(k,447) = -(rxt(k,268) + rxt(k,274) + rxt(k,279)) * y(k,54)
         mat(k,570) = rxt(k,166)*y(k,87)
         mat(k,304) = rxt(k,171)*y(k,96)
         mat(k,545) = rxt(k,166)*y(k,28)
         mat(k,822) = mat(k,822) + rxt(k,171)*y(k,29)
         mat(k,344) = -(rxt(k,131)*y(k,61) + (rxt(k,132) + rxt(k,133) + rxt(k,134) &
                      ) * y(k,62) + rxt(k,135)*y(k,69) + rxt(k,143)*y(k,96) + rxt(k,296) &
                      *y(k,95))
         mat(k,664) = -rxt(k,131)*y(k,55)
         mat(k,781) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,55)
         mat(k,500) = -rxt(k,135)*y(k,55)
         mat(k,827) = -rxt(k,143)*y(k,55)
         mat(k,357) = -rxt(k,296)*y(k,55)
         mat(k,638) = rxt(k,129)*y(k,88) + rxt(k,293)*y(k,91)
         mat(k,500) = mat(k,500) + rxt(k,294)*y(k,91)
         mat(k,322) = 1.100_r8*rxt(k,289)*y(k,89) + .200_r8*rxt(k,287)*y(k,90)
         mat(k,332) = rxt(k,129)*y(k,68)
         mat(k,227) = 1.100_r8*rxt(k,289)*y(k,86)
         mat(k,242) = .200_r8*rxt(k,287)*y(k,86)
         mat(k,273) = rxt(k,293)*y(k,68) + rxt(k,294)*y(k,69)
         mat(k,100) = -((rxt(k,147) + rxt(k,148)) * y(k,92))
         mat(k,700) = -(rxt(k,147) + rxt(k,148)) * y(k,56)
         mat(k,339) = rxt(k,132)*y(k,62)
         mat(k,774) = rxt(k,132)*y(k,55)
         mat(k,775) = rxt(k,150)*y(k,63)
         mat(k,594) = rxt(k,150)*y(k,62)
         mat(k,676) = -(rxt(k,131)*y(k,55) + rxt(k,140)*y(k,63) + rxt(k,144)*y(k,87) &
                      + rxt(k,145)*y(k,70) + rxt(k,146)*y(k,68) + rxt(k,167)*y(k,28) &
                      + rxt(k,199)*y(k,5) + rxt(k,239)*y(k,21) + rxt(k,298)*y(k,95))
         mat(k,350) = -rxt(k,131)*y(k,61)
         mat(k,609) = -rxt(k,140)*y(k,61)
         mat(k,556) = -rxt(k,144)*y(k,61)
         mat(k,531) = -rxt(k,145)*y(k,61)
         mat(k,651) = -rxt(k,146)*y(k,61)
         mat(k,584) = -rxt(k,167)*y(k,61)
         mat(k,742) = -rxt(k,199)*y(k,61)
         mat(k,414) = -rxt(k,239)*y(k,61)
         mat(k,363) = -rxt(k,298)*y(k,61)
         mat(k,350) = mat(k,350) + 2.000_r8*rxt(k,133)*y(k,62) + rxt(k,135)*y(k,69) &
                      + rxt(k,143)*y(k,96)
         mat(k,102) = 2.000_r8*rxt(k,147)*y(k,92)
         mat(k,793) = 2.000_r8*rxt(k,133)*y(k,55) + rxt(k,136)*y(k,68) + rxt(k,256) &
                      *y(k,77)
         mat(k,651) = mat(k,651) + rxt(k,136)*y(k,62)
         mat(k,510) = rxt(k,135)*y(k,55) + rxt(k,130)*y(k,88)
         mat(k,429) = rxt(k,256)*y(k,62)
         mat(k,337) = rxt(k,130)*y(k,69)
         mat(k,718) = 2.000_r8*rxt(k,147)*y(k,56)
         mat(k,842) = rxt(k,143)*y(k,55)
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
         mat(k,797) = -((rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,55) + (rxt(k,136) &
                      + rxt(k,138)) * y(k,68) + rxt(k,137)*y(k,70) + rxt(k,149) &
                      *y(k,87) + rxt(k,150)*y(k,63) + rxt(k,151)*y(k,96) + rxt(k,169) &
                      *y(k,28) + rxt(k,200)*y(k,5) + rxt(k,256)*y(k,77))
         mat(k,352) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,62)
         mat(k,655) = -(rxt(k,136) + rxt(k,138)) * y(k,62)
         mat(k,535) = -rxt(k,137)*y(k,62)
         mat(k,560) = -rxt(k,149)*y(k,62)
         mat(k,613) = -rxt(k,150)*y(k,62)
         mat(k,846) = -rxt(k,151)*y(k,62)
         mat(k,588) = -rxt(k,169)*y(k,62)
         mat(k,746) = -rxt(k,200)*y(k,62)
         mat(k,432) = -rxt(k,256)*y(k,62)
         mat(k,746) = mat(k,746) + rxt(k,199)*y(k,61)
         mat(k,415) = rxt(k,239)*y(k,61)
         mat(k,588) = mat(k,588) + rxt(k,167)*y(k,61)
         mat(k,175) = rxt(k,128)*y(k,96)
         mat(k,680) = rxt(k,199)*y(k,5) + rxt(k,239)*y(k,21) + rxt(k,167)*y(k,28) &
                      + 2.000_r8*rxt(k,140)*y(k,63) + rxt(k,146)*y(k,68) + rxt(k,145) &
                      *y(k,70) + rxt(k,144)*y(k,87)
         mat(k,613) = mat(k,613) + 2.000_r8*rxt(k,140)*y(k,61) + rxt(k,141)*y(k,68) &
                      + rxt(k,139)*y(k,87) + rxt(k,142)*y(k,96)
         mat(k,655) = mat(k,655) + rxt(k,146)*y(k,61) + rxt(k,141)*y(k,63)
         mat(k,535) = mat(k,535) + rxt(k,145)*y(k,61)
         mat(k,560) = mat(k,560) + rxt(k,144)*y(k,61) + rxt(k,139)*y(k,63)
         mat(k,846) = mat(k,846) + rxt(k,128)*y(k,52) + rxt(k,142)*y(k,63)
         mat(k,607) = -(rxt(k,139)*y(k,87) + rxt(k,140)*y(k,61) + rxt(k,141)*y(k,68) &
                      + rxt(k,142)*y(k,96) + rxt(k,150)*y(k,62) + rxt(k,235)*y(k,17) &
                      + rxt(k,246)*y(k,35))
         mat(k,554) = -rxt(k,139)*y(k,63)
         mat(k,674) = -rxt(k,140)*y(k,63)
         mat(k,649) = -rxt(k,141)*y(k,63)
         mat(k,840) = -rxt(k,142)*y(k,63)
         mat(k,791) = -rxt(k,150)*y(k,63)
         mat(k,477) = -rxt(k,235)*y(k,63)
         mat(k,128) = -rxt(k,246)*y(k,63)
         mat(k,201) = rxt(k,201)*y(k,68)
         mat(k,895) = rxt(k,168)*y(k,29)
         mat(k,308) = rxt(k,168)*y(k,25) + rxt(k,170)*y(k,68) + rxt(k,171)*y(k,96)
         mat(k,291) = rxt(k,215)*y(k,51)
         mat(k,863) = rxt(k,215)*y(k,39) + rxt(k,152)*y(k,96)
         mat(k,791) = mat(k,791) + rxt(k,138)*y(k,68) + rxt(k,137)*y(k,70)
         mat(k,649) = mat(k,649) + rxt(k,201)*y(k,6) + rxt(k,170)*y(k,29) + rxt(k,138) &
                      *y(k,62)
         mat(k,529) = rxt(k,137)*y(k,62)
         mat(k,840) = mat(k,840) + rxt(k,171)*y(k,29) + rxt(k,152)*y(k,51)
         mat(k,650) = -(rxt(k,108)*y(k,70) + 4._r8*rxt(k,109)*y(k,68) + rxt(k,110) &
                      *y(k,69) + rxt(k,111)*y(k,41) + rxt(k,112)*y(k,43) + rxt(k,117) &
                      *y(k,87) + rxt(k,123)*y(k,96) + (rxt(k,136) + rxt(k,138) &
                      ) * y(k,62) + rxt(k,141)*y(k,63) + rxt(k,146)*y(k,61) + rxt(k,170) &
                      *y(k,29) + rxt(k,172)*y(k,28) + rxt(k,175)*y(k,49) + rxt(k,178) &
                      *y(k,54) + rxt(k,201)*y(k,6) + rxt(k,202)*y(k,5) + rxt(k,204) &
                      *y(k,45) + rxt(k,206)*y(k,53) + rxt(k,236)*y(k,17) + rxt(k,248) &
                      *y(k,72) + (rxt(k,291) + rxt(k,292)) * y(k,89) + rxt(k,293) &
                      *y(k,91))
         mat(k,530) = -rxt(k,108)*y(k,68)
         mat(k,509) = -rxt(k,110)*y(k,68)
         mat(k,919) = -rxt(k,111)*y(k,68)
         mat(k,219) = -rxt(k,112)*y(k,68)
         mat(k,555) = -rxt(k,117)*y(k,68)
         mat(k,841) = -rxt(k,123)*y(k,68)
         mat(k,792) = -(rxt(k,136) + rxt(k,138)) * y(k,68)
         mat(k,608) = -rxt(k,141)*y(k,68)
         mat(k,675) = -rxt(k,146)*y(k,68)
         mat(k,309) = -rxt(k,170)*y(k,68)
         mat(k,583) = -rxt(k,172)*y(k,68)
         mat(k,456) = -rxt(k,175)*y(k,68)
         mat(k,265) = -rxt(k,178)*y(k,68)
         mat(k,202) = -rxt(k,201)*y(k,68)
         mat(k,741) = -rxt(k,202)*y(k,68)
         mat(k,255) = -rxt(k,204)*y(k,68)
         mat(k,247) = -rxt(k,206)*y(k,68)
         mat(k,478) = -rxt(k,236)*y(k,68)
         mat(k,142) = -rxt(k,248)*y(k,68)
         mat(k,231) = -(rxt(k,291) + rxt(k,292)) * y(k,68)
         mat(k,277) = -rxt(k,293)*y(k,68)
         mat(k,440) = rxt(k,115)*y(k,87)
         mat(k,349) = rxt(k,131)*y(k,61) + rxt(k,132)*y(k,62) + rxt(k,135)*y(k,69) &
                      + rxt(k,296)*y(k,95)
         mat(k,675) = mat(k,675) + rxt(k,131)*y(k,55)
         mat(k,792) = mat(k,792) + rxt(k,132)*y(k,55)
         mat(k,509) = mat(k,509) + rxt(k,135)*y(k,55) + rxt(k,250)*y(k,75) &
                      + rxt(k,257)*y(k,77) + rxt(k,295)*y(k,91) + (rxt(k,97)+rxt(k,98)) &
                      *y(k,92) + rxt(k,302)*y(k,97) + rxt(k,306)*y(k,98)
         mat(k,237) = rxt(k,250)*y(k,69)
         mat(k,428) = rxt(k,257)*y(k,69)
         mat(k,326) = rxt(k,287)*y(k,90) + 1.150_r8*rxt(k,288)*y(k,95)
         mat(k,555) = mat(k,555) + rxt(k,115)*y(k,40)
         mat(k,336) = rxt(k,301)*y(k,97)
         mat(k,243) = rxt(k,287)*y(k,86)
         mat(k,277) = mat(k,277) + rxt(k,295)*y(k,69)
         mat(k,717) = (rxt(k,97)+rxt(k,98))*y(k,69)
         mat(k,362) = rxt(k,296)*y(k,55) + 1.150_r8*rxt(k,288)*y(k,86)
         mat(k,841) = mat(k,841) + 2.000_r8*rxt(k,125)*y(k,96)
         mat(k,379) = rxt(k,302)*y(k,69) + rxt(k,301)*y(k,88)
         mat(k,190) = rxt(k,306)*y(k,69)
         mat(k,505) = -(rxt(k,97)*y(k,92) + rxt(k,102)*y(k,93) + rxt(k,110)*y(k,68) &
                      + rxt(k,116)*y(k,40) + rxt(k,130)*y(k,88) + rxt(k,135)*y(k,55) &
                      + rxt(k,250)*y(k,75) + rxt(k,257)*y(k,77) + rxt(k,290)*y(k,89) &
                      + (rxt(k,294) + rxt(k,295)) * y(k,91) + rxt(k,302)*y(k,97) &
                      + rxt(k,306)*y(k,98))
         mat(k,712) = -rxt(k,97)*y(k,69)
         mat(k,76) = -rxt(k,102)*y(k,69)
         mat(k,645) = -rxt(k,110)*y(k,69)
         mat(k,437) = -rxt(k,116)*y(k,69)
         mat(k,335) = -rxt(k,130)*y(k,69)
         mat(k,348) = -rxt(k,135)*y(k,69)
         mat(k,235) = -rxt(k,250)*y(k,69)
         mat(k,424) = -rxt(k,257)*y(k,69)
         mat(k,230) = -rxt(k,290)*y(k,69)
         mat(k,276) = -(rxt(k,294) + rxt(k,295)) * y(k,69)
         mat(k,377) = -rxt(k,302)*y(k,69)
         mat(k,189) = -rxt(k,306)*y(k,69)
         mat(k,758) = rxt(k,193)*y(k,70) + rxt(k,192)*y(k,87)
         mat(k,736) = 2.000_r8*rxt(k,194)*y(k,5) + (rxt(k,196)+rxt(k,197))*y(k,28) &
                      + rxt(k,202)*y(k,68) + rxt(k,198)*y(k,87)
         mat(k,411) = rxt(k,238)*y(k,87)
         mat(k,891) = rxt(k,161)*y(k,70) + rxt(k,159)*y(k,87)
         mat(k,578) = (rxt(k,196)+rxt(k,197))*y(k,5) + (2.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,164))*y(k,28) + rxt(k,172)*y(k,68) + rxt(k,166) &
                      *y(k,87) + rxt(k,174)*y(k,96)
         mat(k,437) = mat(k,437) + rxt(k,119)*y(k,70) + rxt(k,113)*y(k,87)
         mat(k,172) = rxt(k,128)*y(k,96)
         mat(k,348) = mat(k,348) + rxt(k,134)*y(k,62)
         mat(k,101) = rxt(k,148)*y(k,92)
         mat(k,670) = rxt(k,145)*y(k,70) + rxt(k,298)*y(k,95)
         mat(k,787) = rxt(k,134)*y(k,55) + rxt(k,136)*y(k,68) + rxt(k,137)*y(k,70)
         mat(k,603) = rxt(k,141)*y(k,68) + rxt(k,139)*y(k,87)
         mat(k,645) = mat(k,645) + rxt(k,202)*y(k,5) + rxt(k,172)*y(k,28) + rxt(k,136) &
                      *y(k,62) + rxt(k,141)*y(k,63) + 2.000_r8*rxt(k,109)*y(k,68) &
                      + 2.000_r8*rxt(k,108)*y(k,70) + rxt(k,117)*y(k,87) + rxt(k,101) &
                      *y(k,93) + rxt(k,123)*y(k,96)
         mat(k,505) = mat(k,505) + 2.000_r8*rxt(k,102)*y(k,93)
         mat(k,525) = rxt(k,193)*y(k,3) + rxt(k,161)*y(k,25) + rxt(k,119)*y(k,40) &
                      + rxt(k,145)*y(k,61) + rxt(k,137)*y(k,62) + 2.000_r8*rxt(k,108) &
                      *y(k,68) + rxt(k,252)*y(k,75) + rxt(k,258)*y(k,77) &
                      + 2.000_r8*rxt(k,118)*y(k,87) + 2.000_r8*rxt(k,99)*y(k,92) &
                      + rxt(k,124)*y(k,96)
         mat(k,235) = mat(k,235) + rxt(k,252)*y(k,70)
         mat(k,424) = mat(k,424) + rxt(k,258)*y(k,70)
         mat(k,550) = rxt(k,192)*y(k,3) + rxt(k,198)*y(k,5) + rxt(k,238)*y(k,21) &
                      + rxt(k,159)*y(k,25) + rxt(k,166)*y(k,28) + rxt(k,113)*y(k,40) &
                      + rxt(k,139)*y(k,63) + rxt(k,117)*y(k,68) + 2.000_r8*rxt(k,118) &
                      *y(k,70) + 2.000_r8*rxt(k,127)*y(k,87) + rxt(k,122)*y(k,96)
         mat(k,712) = mat(k,712) + rxt(k,148)*y(k,56) + 2.000_r8*rxt(k,99)*y(k,70)
         mat(k,76) = mat(k,76) + rxt(k,101)*y(k,68) + 2.000_r8*rxt(k,102)*y(k,69)
         mat(k,361) = rxt(k,298)*y(k,61)
         mat(k,836) = rxt(k,174)*y(k,28) + rxt(k,128)*y(k,52) + rxt(k,123)*y(k,68) &
                      + rxt(k,124)*y(k,70) + rxt(k,122)*y(k,87)
         mat(k,526) = -(rxt(k,99)*y(k,92) + rxt(k,108)*y(k,68) + rxt(k,118)*y(k,87) &
                      + rxt(k,119)*y(k,40) + rxt(k,124)*y(k,96) + rxt(k,137)*y(k,62) &
                      + rxt(k,145)*y(k,61) + rxt(k,161)*y(k,25) + rxt(k,193)*y(k,3) &
                      + rxt(k,252)*y(k,75) + rxt(k,258)*y(k,77))
         mat(k,713) = -rxt(k,99)*y(k,70)
         mat(k,646) = -rxt(k,108)*y(k,70)
         mat(k,551) = -rxt(k,118)*y(k,70)
         mat(k,438) = -rxt(k,119)*y(k,70)
         mat(k,837) = -rxt(k,124)*y(k,70)
         mat(k,788) = -rxt(k,137)*y(k,70)
         mat(k,671) = -rxt(k,145)*y(k,70)
         mat(k,892) = -rxt(k,161)*y(k,70)
         mat(k,759) = -rxt(k,193)*y(k,70)
         mat(k,236) = -rxt(k,252)*y(k,70)
         mat(k,425) = -rxt(k,258)*y(k,70)
         mat(k,646) = mat(k,646) + rxt(k,110)*y(k,69)
         mat(k,506) = rxt(k,110)*y(k,68)
         mat(k,131) = -(rxt(k,259)*y(k,77))
         mat(k,419) = -rxt(k,259)*y(k,71)
         mat(k,729) = rxt(k,195)*y(k,28)
         mat(k,569) = rxt(k,195)*y(k,5) + 2.000_r8*rxt(k,165)*y(k,28)
         mat(k,136) = -(rxt(k,248)*y(k,68) + rxt(k,249)*y(k,96))
         mat(k,621) = -rxt(k,248)*y(k,72)
         mat(k,811) = -rxt(k,249)*y(k,72)
         mat(k,232) = -(rxt(k,250)*y(k,69) + rxt(k,252)*y(k,70) + rxt(k,255)*y(k,96))
         mat(k,494) = -rxt(k,250)*y(k,75)
         mat(k,521) = -rxt(k,252)*y(k,75)
         mat(k,820) = -rxt(k,255)*y(k,75)
         mat(k,422) = -(rxt(k,253)*y(k,5) + rxt(k,254)*y(k,28) + rxt(k,256)*y(k,62) &
                      + rxt(k,257)*y(k,69) + rxt(k,258)*y(k,70) + rxt(k,259)*y(k,71) &
                      + rxt(k,260)*y(k,96))
         mat(k,733) = -rxt(k,253)*y(k,77)
         mat(k,574) = -rxt(k,254)*y(k,77)
         mat(k,784) = -rxt(k,256)*y(k,77)
         mat(k,503) = -rxt(k,257)*y(k,77)
         mat(k,523) = -rxt(k,258)*y(k,77)
         mat(k,133) = -rxt(k,259)*y(k,77)
         mat(k,832) = -rxt(k,260)*y(k,77)
         mat(k,641) = rxt(k,248)*y(k,72)
         mat(k,503) = mat(k,503) + rxt(k,250)*y(k,75)
         mat(k,523) = mat(k,523) + rxt(k,252)*y(k,75)
         mat(k,140) = rxt(k,248)*y(k,68)
         mat(k,233) = rxt(k,250)*y(k,69) + rxt(k,252)*y(k,70) + rxt(k,255)*y(k,96)
         mat(k,832) = mat(k,832) + rxt(k,255)*y(k,75)
         mat(k,297) = -(rxt(k,251)*y(k,96))
         mat(k,825) = -rxt(k,251)*y(k,78)
         mat(k,732) = rxt(k,253)*y(k,77)
         mat(k,571) = rxt(k,254)*y(k,77)
         mat(k,126) = rxt(k,246)*y(k,63) + (rxt(k,247)+.500_r8*rxt(k,261))*y(k,96)
         mat(k,779) = rxt(k,256)*y(k,77)
         mat(k,597) = rxt(k,246)*y(k,35)
         mat(k,497) = rxt(k,257)*y(k,77)
         mat(k,522) = rxt(k,258)*y(k,77)
         mat(k,132) = rxt(k,259)*y(k,77)
         mat(k,139) = rxt(k,249)*y(k,96)
         mat(k,421) = rxt(k,253)*y(k,5) + rxt(k,254)*y(k,28) + rxt(k,256)*y(k,62) &
                      + rxt(k,257)*y(k,69) + rxt(k,258)*y(k,70) + rxt(k,259)*y(k,71) &
                      + rxt(k,260)*y(k,96)
         mat(k,825) = mat(k,825) + (rxt(k,247)+.500_r8*rxt(k,261))*y(k,35) &
                      + rxt(k,249)*y(k,72) + rxt(k,260)*y(k,77)
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
         mat(k,105) = -(rxt(k,262)*y(k,100))
         mat(k,931) = -rxt(k,262)*y(k,79)
         mat(k,296) = rxt(k,251)*y(k,96)
         mat(k,808) = rxt(k,251)*y(k,78)
         mat(k,320) = -(rxt(k,287)*y(k,90) + rxt(k,288)*y(k,95) + rxt(k,289)*y(k,89))
         mat(k,240) = -rxt(k,287)*y(k,86)
         mat(k,355) = -rxt(k,288)*y(k,86)
         mat(k,225) = -rxt(k,289)*y(k,86)
         mat(k,552) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,40) + rxt(k,117) &
                      *y(k,68) + rxt(k,118)*y(k,70) + rxt(k,122)*y(k,96) &
                      + 4._r8*rxt(k,127)*y(k,87) + rxt(k,139)*y(k,63) + rxt(k,144) &
                      *y(k,61) + rxt(k,149)*y(k,62) + (rxt(k,159) + rxt(k,160) &
                      ) * y(k,25) + rxt(k,166)*y(k,28) + rxt(k,192)*y(k,3) + rxt(k,198) &
                      *y(k,5) + rxt(k,238)*y(k,21))
         mat(k,439) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,87)
         mat(k,647) = -rxt(k,117)*y(k,87)
         mat(k,527) = -rxt(k,118)*y(k,87)
         mat(k,838) = -rxt(k,122)*y(k,87)
         mat(k,605) = -rxt(k,139)*y(k,87)
         mat(k,672) = -rxt(k,144)*y(k,87)
         mat(k,789) = -rxt(k,149)*y(k,87)
         mat(k,893) = -(rxt(k,159) + rxt(k,160)) * y(k,87)
         mat(k,580) = -rxt(k,166)*y(k,87)
         mat(k,760) = -rxt(k,192)*y(k,87)
         mat(k,738) = -rxt(k,198)*y(k,87)
         mat(k,412) = -rxt(k,238)*y(k,87)
         mat(k,760) = mat(k,760) + rxt(k,191)*y(k,17)
         mat(k,738) = mat(k,738) + rxt(k,203)*y(k,96)
         mat(k,475) = rxt(k,191)*y(k,3) + rxt(k,155)*y(k,25) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,68)
         mat(k,210) = rxt(k,220)*y(k,25) + rxt(k,221)*y(k,96)
         mat(k,148) = rxt(k,223)*y(k,25) + rxt(k,224)*y(k,96)
         mat(k,412) = mat(k,412) + rxt(k,162)*y(k,28) + rxt(k,239)*y(k,61)
         mat(k,396) = rxt(k,243)*y(k,92)
         mat(k,893) = mat(k,893) + rxt(k,155)*y(k,17) + rxt(k,220)*y(k,18) &
                      + rxt(k,223)*y(k,20) + rxt(k,158)*y(k,43)
         mat(k,580) = mat(k,580) + rxt(k,162)*y(k,21) + rxt(k,173)*y(k,96)
         mat(k,283) = rxt(k,245)*y(k,96)
         mat(k,127) = .500_r8*rxt(k,261)*y(k,96)
         mat(k,439) = mat(k,439) + rxt(k,116)*y(k,69)
         mat(k,218) = rxt(k,158)*y(k,25) + rxt(k,112)*y(k,68) + rxt(k,121)*y(k,96)
         mat(k,672) = mat(k,672) + rxt(k,239)*y(k,21)
         mat(k,605) = mat(k,605) + rxt(k,235)*y(k,17) + rxt(k,142)*y(k,96)
         mat(k,647) = mat(k,647) + rxt(k,236)*y(k,17) + rxt(k,112)*y(k,43)
         mat(k,507) = rxt(k,116)*y(k,40)
         mat(k,527) = mat(k,527) + rxt(k,124)*y(k,96)
         mat(k,299) = rxt(k,251)*y(k,96)
         mat(k,714) = rxt(k,243)*y(k,23)
         mat(k,838) = mat(k,838) + rxt(k,203)*y(k,5) + rxt(k,221)*y(k,18) + rxt(k,224) &
                      *y(k,20) + rxt(k,173)*y(k,28) + rxt(k,245)*y(k,31) &
                      + .500_r8*rxt(k,261)*y(k,35) + rxt(k,121)*y(k,43) + rxt(k,142) &
                      *y(k,63) + rxt(k,124)*y(k,70) + rxt(k,251)*y(k,78)
         mat(k,331) = -(rxt(k,129)*y(k,68) + rxt(k,130)*y(k,69) + rxt(k,301)*y(k,97))
         mat(k,637) = -rxt(k,129)*y(k,88)
         mat(k,499) = -rxt(k,130)*y(k,88)
         mat(k,372) = -rxt(k,301)*y(k,88)
         mat(k,637) = mat(k,637) + rxt(k,291)*y(k,89)
         mat(k,321) = .900_r8*rxt(k,289)*y(k,89) + .800_r8*rxt(k,287)*y(k,90)
         mat(k,226) = rxt(k,291)*y(k,68) + .900_r8*rxt(k,289)*y(k,86)
         mat(k,241) = .800_r8*rxt(k,287)*y(k,86)
         mat(k,223) = -(rxt(k,289)*y(k,86) + rxt(k,290)*y(k,69) + (rxt(k,291) &
                      + rxt(k,292)) * y(k,68))
         mat(k,317) = -rxt(k,289)*y(k,89)
         mat(k,493) = -rxt(k,290)*y(k,89)
         mat(k,626) = -(rxt(k,291) + rxt(k,292)) * y(k,89)
         mat(k,239) = -(rxt(k,287)*y(k,86))
         mat(k,318) = -rxt(k,287)*y(k,90)
         mat(k,340) = rxt(k,296)*y(k,95)
         mat(k,661) = rxt(k,298)*y(k,95)
         mat(k,628) = rxt(k,291)*y(k,89)
         mat(k,495) = rxt(k,295)*y(k,91)
         mat(k,224) = rxt(k,291)*y(k,68)
         mat(k,269) = rxt(k,295)*y(k,69)
         mat(k,354) = rxt(k,296)*y(k,55) + rxt(k,298)*y(k,61)
         mat(k,270) = -(rxt(k,293)*y(k,68) + (rxt(k,294) + rxt(k,295)) * y(k,69))
         mat(k,632) = -rxt(k,293)*y(k,91)
         mat(k,496) = -(rxt(k,294) + rxt(k,295)) * y(k,91)
         mat(k,329) = rxt(k,301)*y(k,97)
         mat(k,369) = rxt(k,301)*y(k,88)
         mat(k,719) = -(rxt(k,94)*y(k,41) + rxt(k,95)*y(k,100) + (rxt(k,97) + rxt(k,98) &
                      ) * y(k,69) + rxt(k,99)*y(k,70) + (rxt(k,147) + rxt(k,148) &
                      ) * y(k,56) + rxt(k,180)*y(k,8) + rxt(k,181)*y(k,9) + rxt(k,182) &
                      *y(k,11) + rxt(k,183)*y(k,12) + rxt(k,184)*y(k,13) + rxt(k,185) &
                      *y(k,14) + rxt(k,186)*y(k,15) + (rxt(k,187) + rxt(k,188) &
                      ) * y(k,49) + rxt(k,207)*y(k,10) + rxt(k,208)*y(k,24) + rxt(k,209) &
                      *y(k,42) + (rxt(k,210) + rxt(k,211)) * y(k,45) + rxt(k,216) &
                      *y(k,33) + rxt(k,217)*y(k,34) + rxt(k,230)*y(k,16) + rxt(k,231) &
                      *y(k,18) + rxt(k,232)*y(k,46) + rxt(k,233)*y(k,47) + rxt(k,234) &
                      *y(k,48) + (rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,23))
         mat(k,921) = -rxt(k,94)*y(k,92)
         mat(k,948) = -rxt(k,95)*y(k,92)
         mat(k,511) = -(rxt(k,97) + rxt(k,98)) * y(k,92)
         mat(k,532) = -rxt(k,99)*y(k,92)
         mat(k,103) = -(rxt(k,147) + rxt(k,148)) * y(k,92)
         mat(k,28) = -rxt(k,180)*y(k,92)
         mat(k,55) = -rxt(k,181)*y(k,92)
         mat(k,36) = -rxt(k,182)*y(k,92)
         mat(k,66) = -rxt(k,183)*y(k,92)
         mat(k,40) = -rxt(k,184)*y(k,92)
         mat(k,71) = -rxt(k,185)*y(k,92)
         mat(k,44) = -rxt(k,186)*y(k,92)
         mat(k,457) = -(rxt(k,187) + rxt(k,188)) * y(k,92)
         mat(k,61) = -rxt(k,207)*y(k,92)
         mat(k,154) = -rxt(k,208)*y(k,92)
         mat(k,32) = -rxt(k,209)*y(k,92)
         mat(k,256) = -(rxt(k,210) + rxt(k,211)) * y(k,92)
         mat(k,96) = -rxt(k,216)*y(k,92)
         mat(k,87) = -rxt(k,217)*y(k,92)
         mat(k,166) = -rxt(k,230)*y(k,92)
         mat(k,211) = -rxt(k,231)*y(k,92)
         mat(k,82) = -rxt(k,232)*y(k,92)
         mat(k,91) = -rxt(k,233)*y(k,92)
         mat(k,112) = -rxt(k,234)*y(k,92)
         mat(k,400) = -(rxt(k,242) + rxt(k,243) + rxt(k,244)) * y(k,92)
         mat(k,511) = mat(k,511) + rxt(k,130)*y(k,88)
         mat(k,328) = .850_r8*rxt(k,288)*y(k,95)
         mat(k,338) = rxt(k,130)*y(k,69)
         mat(k,364) = .850_r8*rxt(k,288)*y(k,86)
         mat(k,75) = -(rxt(k,101)*y(k,68) + rxt(k,102)*y(k,69))
         mat(k,619) = -rxt(k,101)*y(k,93)
         mat(k,489) = -rxt(k,102)*y(k,93)
         mat(k,191) = rxt(k,103)*y(k,94)
         mat(k,619) = mat(k,619) + rxt(k,105)*y(k,94)
         mat(k,489) = mat(k,489) + rxt(k,106)*y(k,94)
         mat(k,519) = rxt(k,107)*y(k,94)
         mat(k,77) = rxt(k,103)*y(k,32) + rxt(k,105)*y(k,68) + rxt(k,106)*y(k,69) &
                      + rxt(k,107)*y(k,70)
         mat(k,78) = -(rxt(k,103)*y(k,32) + rxt(k,105)*y(k,68) + rxt(k,106)*y(k,69) &
                      + rxt(k,107)*y(k,70))
         mat(k,192) = -rxt(k,103)*y(k,94)
         mat(k,620) = -rxt(k,105)*y(k,94)
         mat(k,490) = -rxt(k,106)*y(k,94)
         mat(k,520) = -rxt(k,107)*y(k,94)
         mat(k,490) = mat(k,490) + rxt(k,97)*y(k,92)
         mat(k,695) = rxt(k,97)*y(k,69)
         mat(k,358) = -(rxt(k,288)*y(k,86) + rxt(k,296)*y(k,55) + rxt(k,298)*y(k,61))
         mat(k,323) = -rxt(k,288)*y(k,95)
         mat(k,345) = -rxt(k,296)*y(k,95)
         mat(k,665) = -rxt(k,298)*y(k,95)
         mat(k,195) = rxt(k,299)*y(k,97)
         mat(k,501) = rxt(k,290)*y(k,89) + rxt(k,294)*y(k,91) + rxt(k,302)*y(k,97) &
                      + rxt(k,306)*y(k,98)
         mat(k,228) = rxt(k,290)*y(k,69)
         mat(k,274) = rxt(k,294)*y(k,69)
         mat(k,374) = rxt(k,299)*y(k,32) + rxt(k,302)*y(k,69)
         mat(k,187) = rxt(k,306)*y(k,69)
         mat(k,847) = -(rxt(k,120)*y(k,41) + rxt(k,121)*y(k,43) + rxt(k,122)*y(k,87) &
                      + rxt(k,123)*y(k,68) + rxt(k,124)*y(k,70) + (4._r8*rxt(k,125) &
                      + 4._r8*rxt(k,126)) * y(k,96) + rxt(k,128)*y(k,52) + rxt(k,142) &
                      *y(k,63) + rxt(k,143)*y(k,55) + rxt(k,151)*y(k,62) + rxt(k,152) &
                      *y(k,51) + rxt(k,171)*y(k,29) + (rxt(k,173) + rxt(k,174) &
                      ) * y(k,28) + rxt(k,176)*y(k,49) + rxt(k,179)*y(k,54) + rxt(k,203) &
                      *y(k,5) + rxt(k,205)*y(k,45) + rxt(k,219)*y(k,16) + rxt(k,221) &
                      *y(k,18) + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,20) + rxt(k,226) &
                      *y(k,24) + rxt(k,227)*y(k,46) + rxt(k,228)*y(k,47) + rxt(k,229) &
                      *y(k,48) + rxt(k,237)*y(k,17) + rxt(k,240)*y(k,22) + rxt(k,241) &
                      *y(k,23) + rxt(k,245)*y(k,31) + (rxt(k,247) + rxt(k,261) &
                      ) * y(k,35) + rxt(k,249)*y(k,72) + rxt(k,251)*y(k,78) + rxt(k,255) &
                      *y(k,75) + rxt(k,260)*y(k,77))
         mat(k,925) = -rxt(k,120)*y(k,96)
         mat(k,220) = -rxt(k,121)*y(k,96)
         mat(k,561) = -rxt(k,122)*y(k,96)
         mat(k,656) = -rxt(k,123)*y(k,96)
         mat(k,536) = -rxt(k,124)*y(k,96)
         mat(k,176) = -rxt(k,128)*y(k,96)
         mat(k,614) = -rxt(k,142)*y(k,96)
         mat(k,353) = -rxt(k,143)*y(k,96)
         mat(k,798) = -rxt(k,151)*y(k,96)
         mat(k,870) = -rxt(k,152)*y(k,96)
         mat(k,311) = -rxt(k,171)*y(k,96)
         mat(k,589) = -(rxt(k,173) + rxt(k,174)) * y(k,96)
         mat(k,461) = -rxt(k,176)*y(k,96)
         mat(k,266) = -rxt(k,179)*y(k,96)
         mat(k,747) = -rxt(k,203)*y(k,96)
         mat(k,259) = -rxt(k,205)*y(k,96)
         mat(k,168) = -rxt(k,219)*y(k,96)
         mat(k,213) = -rxt(k,221)*y(k,96)
         mat(k,47) = -rxt(k,222)*y(k,96)
         mat(k,149) = -rxt(k,224)*y(k,96)
         mat(k,156) = -rxt(k,226)*y(k,96)
         mat(k,83) = -rxt(k,227)*y(k,96)
         mat(k,92) = -rxt(k,228)*y(k,96)
         mat(k,113) = -rxt(k,229)*y(k,96)
         mat(k,484) = -rxt(k,237)*y(k,96)
         mat(k,162) = -rxt(k,240)*y(k,96)
         mat(k,402) = -rxt(k,241)*y(k,96)
         mat(k,285) = -rxt(k,245)*y(k,96)
         mat(k,129) = -(rxt(k,247) + rxt(k,261)) * y(k,96)
         mat(k,143) = -rxt(k,249)*y(k,96)
         mat(k,301) = -rxt(k,251)*y(k,96)
         mat(k,238) = -rxt(k,255)*y(k,96)
         mat(k,433) = -rxt(k,260)*y(k,96)
         mat(k,484) = mat(k,484) + rxt(k,236)*y(k,68)
         mat(k,162) = mat(k,162) + .300_r8*rxt(k,240)*y(k,96)
         mat(k,402) = mat(k,402) + rxt(k,242)*y(k,92)
         mat(k,902) = rxt(k,160)*y(k,87)
         mat(k,292) = rxt(k,214)*y(k,100)
         mat(k,441) = rxt(k,119)*y(k,70) + 2.000_r8*rxt(k,114)*y(k,87)
         mat(k,925) = mat(k,925) + rxt(k,111)*y(k,68) + rxt(k,94)*y(k,92)
         mat(k,220) = mat(k,220) + rxt(k,112)*y(k,68)
         mat(k,259) = mat(k,259) + rxt(k,204)*y(k,68) + rxt(k,210)*y(k,92)
         mat(k,461) = mat(k,461) + rxt(k,175)*y(k,68) + rxt(k,187)*y(k,92)
         mat(k,250) = rxt(k,206)*y(k,68)
         mat(k,266) = mat(k,266) + rxt(k,178)*y(k,68)
         mat(k,681) = rxt(k,144)*y(k,87)
         mat(k,614) = mat(k,614) + rxt(k,139)*y(k,87)
         mat(k,656) = mat(k,656) + rxt(k,236)*y(k,17) + rxt(k,111)*y(k,41) &
                      + rxt(k,112)*y(k,43) + rxt(k,204)*y(k,45) + rxt(k,175)*y(k,49) &
                      + rxt(k,206)*y(k,53) + rxt(k,178)*y(k,54) + rxt(k,117)*y(k,87)
         mat(k,536) = mat(k,536) + rxt(k,119)*y(k,40) + rxt(k,118)*y(k,87)
         mat(k,561) = mat(k,561) + rxt(k,160)*y(k,25) + 2.000_r8*rxt(k,114)*y(k,40) &
                      + rxt(k,144)*y(k,61) + rxt(k,139)*y(k,63) + rxt(k,117)*y(k,68) &
                      + rxt(k,118)*y(k,70)
         mat(k,723) = rxt(k,242)*y(k,23) + rxt(k,94)*y(k,41) + rxt(k,210)*y(k,45) &
                      + rxt(k,187)*y(k,49) + 2.000_r8*rxt(k,95)*y(k,100)
         mat(k,847) = mat(k,847) + .300_r8*rxt(k,240)*y(k,22)
         mat(k,952) = rxt(k,214)*y(k,39) + 2.000_r8*rxt(k,95)*y(k,92)
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
         mat(k,375) = -(rxt(k,299)*y(k,32) + rxt(k,301)*y(k,88) + rxt(k,302)*y(k,69))
         mat(k,196) = -rxt(k,299)*y(k,97)
         mat(k,334) = -rxt(k,301)*y(k,97)
         mat(k,502) = -rxt(k,302)*y(k,97)
         mat(k,640) = rxt(k,292)*y(k,89) + rxt(k,293)*y(k,91) + rxt(k,305)*y(k,98) &
                      + rxt(k,311)*y(k,99)
         mat(k,324) = rxt(k,303)*y(k,98) + rxt(k,308)*y(k,99)
         mat(k,229) = rxt(k,292)*y(k,68)
         mat(k,275) = rxt(k,293)*y(k,68)
         mat(k,188) = rxt(k,305)*y(k,68) + rxt(k,303)*y(k,86)
         mat(k,182) = rxt(k,311)*y(k,68) + rxt(k,308)*y(k,86)
         mat(k,185) = -(rxt(k,303)*y(k,86) + rxt(k,305)*y(k,68) + rxt(k,306)*y(k,69))
         mat(k,316) = -rxt(k,303)*y(k,98)
         mat(k,623) = -rxt(k,305)*y(k,98)
         mat(k,492) = -rxt(k,306)*y(k,98)
         mat(k,316) = mat(k,316) + rxt(k,307)*y(k,99)
         mat(k,179) = rxt(k,307)*y(k,86)
         mat(k,178) = -((rxt(k,307) + rxt(k,308)) * y(k,86) + rxt(k,311)*y(k,68))
         mat(k,315) = -(rxt(k,307) + rxt(k,308)) * y(k,99)
         mat(k,622) = -rxt(k,311)*y(k,99)
         mat(k,956) = -(rxt(k,95)*y(k,92) + rxt(k,214)*y(k,39) + rxt(k,262)*y(k,79))
         mat(k,727) = -rxt(k,95)*y(k,100)
         mat(k,295) = -rxt(k,214)*y(k,100)
         mat(k,108) = -rxt(k,262)*y(k,100)
         mat(k,170) = rxt(k,219)*y(k,96)
         mat(k,488) = rxt(k,237)*y(k,96)
         mat(k,215) = rxt(k,221)*y(k,96)
         mat(k,49) = rxt(k,222)*y(k,96)
         mat(k,151) = rxt(k,224)*y(k,96)
         mat(k,163) = rxt(k,240)*y(k,96)
         mat(k,406) = rxt(k,241)*y(k,96)
         mat(k,443) = rxt(k,115)*y(k,87)
         mat(k,929) = rxt(k,120)*y(k,96)
         mat(k,222) = rxt(k,121)*y(k,96)
         mat(k,260) = rxt(k,205)*y(k,96)
         mat(k,115) = rxt(k,229)*y(k,96)
         mat(k,465) = (rxt(k,275)+rxt(k,280))*y(k,53) + (rxt(k,268)+rxt(k,274) &
                       +rxt(k,279))*y(k,54) + rxt(k,176)*y(k,96)
         mat(k,874) = rxt(k,152)*y(k,96)
         mat(k,177) = rxt(k,128)*y(k,96)
         mat(k,252) = (rxt(k,275)+rxt(k,280))*y(k,49)
         mat(k,268) = (rxt(k,268)+rxt(k,274)+rxt(k,279))*y(k,49) + rxt(k,179)*y(k,96)
         mat(k,565) = rxt(k,115)*y(k,40) + rxt(k,122)*y(k,96)
         mat(k,851) = rxt(k,219)*y(k,16) + rxt(k,237)*y(k,17) + rxt(k,221)*y(k,18) &
                      + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,20) + rxt(k,240)*y(k,22) &
                      + rxt(k,241)*y(k,23) + rxt(k,120)*y(k,41) + rxt(k,121)*y(k,43) &
                      + rxt(k,205)*y(k,45) + rxt(k,229)*y(k,48) + rxt(k,176)*y(k,49) &
                      + rxt(k,152)*y(k,51) + rxt(k,128)*y(k,52) + rxt(k,179)*y(k,54) &
                      + rxt(k,122)*y(k,87) + 2.000_r8*rxt(k,125)*y(k,96)
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
         mat(k, 27) = mat(k, 27) + lmat(k, 27)
         mat(k, 29) = mat(k, 29) + lmat(k, 29)
         mat(k, 30) = mat(k, 30) + lmat(k, 30)
         mat(k, 31) = mat(k, 31) + lmat(k, 31)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 34) = mat(k, 34) + lmat(k, 34)
         mat(k, 35) = mat(k, 35) + lmat(k, 35)
         mat(k, 37) = mat(k, 37) + lmat(k, 37)
         mat(k, 38) = mat(k, 38) + lmat(k, 38)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 42) = mat(k, 42) + lmat(k, 42)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 50) = lmat(k, 50)
         mat(k, 51) = lmat(k, 51)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 85) = mat(k, 85) + lmat(k, 85)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = lmat(k, 99)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = lmat(k, 107)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 116) = lmat(k, 116)
         mat(k, 117) = lmat(k, 117)
         mat(k, 118) = lmat(k, 118)
         mat(k, 119) = lmat(k, 119)
         mat(k, 120) = lmat(k, 120)
         mat(k, 121) = lmat(k, 121)
         mat(k, 122) = lmat(k, 122)
         mat(k, 123) = lmat(k, 123)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 137) = lmat(k, 137)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 146) = lmat(k, 146)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 160) = lmat(k, 160)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 173) = lmat(k, 173)
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 181) = lmat(k, 181)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = lmat(k, 186)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = mat(k, 201) + lmat(k, 201)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 212) = mat(k, 212) + lmat(k, 212)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 249) = lmat(k, 249)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 298) = lmat(k, 298)
         mat(k, 300) = lmat(k, 300)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 310) = lmat(k, 310)
         mat(k, 312) = mat(k, 312) + lmat(k, 312)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = lmat(k, 342)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 354) = mat(k, 354) + lmat(k, 354)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 368) = lmat(k, 368)
         mat(k, 373) = lmat(k, 373)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 390) = mat(k, 390) + lmat(k, 390)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 392) = mat(k, 392) + lmat(k, 392)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 398) = lmat(k, 398)
         mat(k, 402) = mat(k, 402) + lmat(k, 402)
         mat(k, 405) = mat(k, 405) + lmat(k, 405)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 408) = mat(k, 408) + lmat(k, 408)
         mat(k, 420) = lmat(k, 420)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 470) = mat(k, 470) + lmat(k, 470)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 487) = lmat(k, 487)
         mat(k, 491) = lmat(k, 491)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 498) = lmat(k, 498)
         mat(k, 501) = mat(k, 501) + lmat(k, 501)
         mat(k, 502) = mat(k, 502) + lmat(k, 502)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 525) = mat(k, 525) + lmat(k, 525)
         mat(k, 526) = mat(k, 526) + lmat(k, 526)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 583) = mat(k, 583) + lmat(k, 583)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 603) = mat(k, 603) + lmat(k, 603)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 608) = mat(k, 608) + lmat(k, 608)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 613) = mat(k, 613) + lmat(k, 613)
         mat(k, 615) = mat(k, 615) + lmat(k, 615)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 623) = mat(k, 623) + lmat(k, 623)
         mat(k, 636) = lmat(k, 636)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 661) = mat(k, 661) + lmat(k, 661)
         mat(k, 662) = lmat(k, 662)
         mat(k, 664) = mat(k, 664) + lmat(k, 664)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 676) = mat(k, 676) + lmat(k, 676)
         mat(k, 717) = mat(k, 717) + lmat(k, 717)
         mat(k, 719) = mat(k, 719) + lmat(k, 719)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 745) = mat(k, 745) + lmat(k, 745)
         mat(k, 767) = mat(k, 767) + lmat(k, 767)
         mat(k, 792) = mat(k, 792) + lmat(k, 792)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 797) = mat(k, 797) + lmat(k, 797)
         mat(k, 798) = mat(k, 798) + lmat(k, 798)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 847) = mat(k, 847) + lmat(k, 847)
         mat(k, 869) = lmat(k, 869)
         mat(k, 870) = mat(k, 870) + lmat(k, 870)
         mat(k, 871) = mat(k, 871) + lmat(k, 871)
         mat(k, 904) = mat(k, 904) + lmat(k, 904)
         mat(k, 928) = mat(k, 928) + lmat(k, 928)
         mat(k, 938) = lmat(k, 938)
         mat(k, 946) = lmat(k, 946)
         mat(k, 948) = mat(k, 948) + lmat(k, 948)
         mat(k, 952) = mat(k, 952) + lmat(k, 952)
         mat(k, 955) = lmat(k, 955)
         mat(k, 956) = mat(k, 956) + lmat(k, 956)
         mat(k, 111) = 0._r8
         mat(k, 197) = 0._r8
         mat(k, 251) = 0._r8
         mat(k, 271) = 0._r8
         mat(k, 272) = 0._r8
         mat(k, 280) = 0._r8
         mat(k, 281) = 0._r8
         mat(k, 282) = 0._r8
         mat(k, 284) = 0._r8
         mat(k, 302) = 0._r8
         mat(k, 314) = 0._r8
         mat(k, 319) = 0._r8
         mat(k, 325) = 0._r8
         mat(k, 327) = 0._r8
         mat(k, 330) = 0._r8
         mat(k, 333) = 0._r8
         mat(k, 343) = 0._r8
         mat(k, 346) = 0._r8
         mat(k, 351) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 359) = 0._r8
         mat(k, 360) = 0._r8
         mat(k, 365) = 0._r8
         mat(k, 366) = 0._r8
         mat(k, 371) = 0._r8
         mat(k, 376) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 380) = 0._r8
         mat(k, 381) = 0._r8
         mat(k, 382) = 0._r8
         mat(k, 383) = 0._r8
         mat(k, 388) = 0._r8
         mat(k, 389) = 0._r8
         mat(k, 395) = 0._r8
         mat(k, 397) = 0._r8
         mat(k, 399) = 0._r8
         mat(k, 401) = 0._r8
         mat(k, 403) = 0._r8
         mat(k, 409) = 0._r8
         mat(k, 416) = 0._r8
         mat(k, 418) = 0._r8
         mat(k, 426) = 0._r8
         mat(k, 435) = 0._r8
         mat(k, 451) = 0._r8
         mat(k, 452) = 0._r8
         mat(k, 453) = 0._r8
         mat(k, 455) = 0._r8
         mat(k, 458) = 0._r8
         mat(k, 459) = 0._r8
         mat(k, 460) = 0._r8
         mat(k, 464) = 0._r8
         mat(k, 468) = 0._r8
         mat(k, 469) = 0._r8
         mat(k, 473) = 0._r8
         mat(k, 474) = 0._r8
         mat(k, 476) = 0._r8
         mat(k, 479) = 0._r8
         mat(k, 480) = 0._r8
         mat(k, 481) = 0._r8
         mat(k, 483) = 0._r8
         mat(k, 508) = 0._r8
         mat(k, 512) = 0._r8
         mat(k, 513) = 0._r8
         mat(k, 514) = 0._r8
         mat(k, 515) = 0._r8
         mat(k, 516) = 0._r8
         mat(k, 517) = 0._r8
         mat(k, 518) = 0._r8
         mat(k, 538) = 0._r8
         mat(k, 539) = 0._r8
         mat(k, 549) = 0._r8
         mat(k, 557) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 575) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 582) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 590) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 598) = 0._r8
         mat(k, 599) = 0._r8
         mat(k, 600) = 0._r8
         mat(k, 601) = 0._r8
         mat(k, 604) = 0._r8
         mat(k, 606) = 0._r8
         mat(k, 610) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 612) = 0._r8
         mat(k, 616) = 0._r8
         mat(k, 617) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 627) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 652) = 0._r8
         mat(k, 657) = 0._r8
         mat(k, 660) = 0._r8
         mat(k, 663) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 668) = 0._r8
         mat(k, 677) = 0._r8
         mat(k, 682) = 0._r8
         mat(k, 684) = 0._r8
         mat(k, 685) = 0._r8
         mat(k, 716) = 0._r8
         mat(k, 722) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 735) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 761) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 763) = 0._r8
         mat(k, 764) = 0._r8
         mat(k, 765) = 0._r8
         mat(k, 768) = 0._r8
         mat(k, 769) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 773) = 0._r8
         mat(k, 778) = 0._r8
         mat(k, 782) = 0._r8
         mat(k, 783) = 0._r8
         mat(k, 785) = 0._r8
         mat(k, 786) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 801) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 824) = 0._r8
         mat(k, 828) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 854) = 0._r8
         mat(k, 855) = 0._r8
         mat(k, 856) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 858) = 0._r8
         mat(k, 859) = 0._r8
         mat(k, 860) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 862) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 866) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 868) = 0._r8
         mat(k, 872) = 0._r8
         mat(k, 873) = 0._r8
         mat(k, 884) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 896) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 903) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 909) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 913) = 0._r8
         mat(k, 914) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 916) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 920) = 0._r8
         mat(k, 922) = 0._r8
         mat(k, 923) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 926) = 0._r8
         mat(k, 934) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 936) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 940) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 943) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 949) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 954) = 0._r8
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
         mat(k, 27) = mat(k, 27) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 34) = mat(k, 34) - dti(k)
         mat(k, 38) = mat(k, 38) - dti(k)
         mat(k, 42) = mat(k, 42) - dti(k)
         mat(k, 46) = mat(k, 46) - dti(k)
         mat(k, 50) = mat(k, 50) - dti(k)
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 58) = mat(k, 58) - dti(k)
         mat(k, 63) = mat(k, 63) - dti(k)
         mat(k, 68) = mat(k, 68) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 78) = mat(k, 78) - dti(k)
         mat(k, 80) = mat(k, 80) - dti(k)
         mat(k, 85) = mat(k, 85) - dti(k)
         mat(k, 89) = mat(k, 89) - dti(k)
         mat(k, 94) = mat(k, 94) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 110) = mat(k, 110) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 158) = mat(k, 158) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 171) = mat(k, 171) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 223) = mat(k, 223) - dti(k)
         mat(k, 232) = mat(k, 232) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 245) = mat(k, 245) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 262) = mat(k, 262) - dti(k)
         mat(k, 270) = mat(k, 270) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 297) = mat(k, 297) - dti(k)
         mat(k, 305) = mat(k, 305) - dti(k)
         mat(k, 320) = mat(k, 320) - dti(k)
         mat(k, 331) = mat(k, 331) - dti(k)
         mat(k, 344) = mat(k, 344) - dti(k)
         mat(k, 358) = mat(k, 358) - dti(k)
         mat(k, 375) = mat(k, 375) - dti(k)
         mat(k, 390) = mat(k, 390) - dti(k)
         mat(k, 408) = mat(k, 408) - dti(k)
         mat(k, 422) = mat(k, 422) - dti(k)
         mat(k, 436) = mat(k, 436) - dti(k)
         mat(k, 450) = mat(k, 450) - dti(k)
         mat(k, 472) = mat(k, 472) - dti(k)
         mat(k, 505) = mat(k, 505) - dti(k)
         mat(k, 526) = mat(k, 526) - dti(k)
         mat(k, 552) = mat(k, 552) - dti(k)
         mat(k, 581) = mat(k, 581) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 650) = mat(k, 650) - dti(k)
         mat(k, 676) = mat(k, 676) - dti(k)
         mat(k, 719) = mat(k, 719) - dti(k)
         mat(k, 744) = mat(k, 744) - dti(k)
         mat(k, 767) = mat(k, 767) - dti(k)
         mat(k, 797) = mat(k, 797) - dti(k)
         mat(k, 847) = mat(k, 847) - dti(k)
         mat(k, 871) = mat(k, 871) - dti(k)
         mat(k, 904) = mat(k, 904) - dti(k)
         mat(k, 928) = mat(k, 928) - dti(k)
         mat(k, 956) = mat(k, 956) - dti(k)
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
