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
         mat(k,486) = rxt(k,192)*y(k,24)
         mat(k,348) = rxt(k,192)*y(k,2)
         mat(k,248) = (rxt(k,254)+rxt(k,259))*y(k,45)
         mat(k,111) = (rxt(k,254)+rxt(k,259))*y(k,40)
         mat(k,502) = -(4._r8*rxt(k,189)*y(k,2) + (rxt(k,190) + rxt(k,191) + rxt(k,192) &
                      ) * y(k,24) + rxt(k,193)*y(k,43) + rxt(k,194)*y(k,50) + rxt(k,195) &
                      *y(k,51) + rxt(k,197)*y(k,53) + rxt(k,198)*y(k,69))
         mat(k,365) = -(rxt(k,190) + rxt(k,191) + rxt(k,192)) * y(k,2)
         mat(k,390) = -rxt(k,193)*y(k,2)
         mat(k,314) = -rxt(k,194)*y(k,2)
         mat(k,415) = -rxt(k,195)*y(k,2)
         mat(k,572) = -rxt(k,197)*y(k,2)
         mat(k,533) = -rxt(k,198)*y(k,2)
         mat(k,69) = rxt(k,196)*y(k,53)
         mat(k,189) = rxt(k,206)*y(k,65)
         mat(k,117) = rxt(k,201)*y(k,53)
         mat(k,572) = mat(k,572) + rxt(k,196)*y(k,3) + rxt(k,201)*y(k,45)
         mat(k,631) = rxt(k,188)*y(k,58)
         mat(k,460) = rxt(k,188)*y(k,55)
         mat(k,440) = rxt(k,206)*y(k,36)
         mat(k,64) = -(rxt(k,196)*y(k,53))
         mat(k,544) = -rxt(k,196)*y(k,3)
         mat(k,487) = rxt(k,195)*y(k,51)
         mat(k,398) = rxt(k,195)*y(k,2)
         mat(k,272) = -(rxt(k,150)*y(k,59) + rxt(k,186)*y(k,58) + rxt(k,230)*y(k,52) &
                      + rxt(k,231)*y(k,53) + rxt(k,232)*y(k,69))
         mat(k,329) = -rxt(k,150)*y(k,14)
         mat(k,450) = -rxt(k,186)*y(k,14)
         mat(k,583) = -rxt(k,230)*y(k,14)
         mat(k,562) = -rxt(k,231)*y(k,14)
         mat(k,523) = -rxt(k,232)*y(k,14)
         mat(k,469) = rxt(k,157)*y(k,24) + rxt(k,234)*y(k,50)
         mat(k,35) = .300_r8*rxt(k,235)*y(k,69)
         mat(k,355) = rxt(k,157)*y(k,18)
         mat(k,304) = rxt(k,234)*y(k,18)
         mat(k,523) = mat(k,523) + .300_r8*rxt(k,235)*y(k,19)
         mat(k,478) = -(rxt(k,157)*y(k,24) + rxt(k,233)*y(k,43) + rxt(k,234)*y(k,50))
         mat(k,364) = -rxt(k,157)*y(k,18)
         mat(k,389) = -rxt(k,233)*y(k,18)
         mat(k,313) = -rxt(k,234)*y(k,18)
         mat(k,37) = .700_r8*rxt(k,235)*y(k,69)
         mat(k,532) = .700_r8*rxt(k,235)*y(k,19)
         mat(k,33) = -(rxt(k,235)*y(k,69))
         mat(k,510) = -rxt(k,235)*y(k,19)
         mat(k,466) = rxt(k,233)*y(k,43)
         mat(k,371) = rxt(k,233)*y(k,18)
         mat(k,347) = 2.000_r8*rxt(k,159)*y(k,24)
         mat(k,192) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,40) + rxt(k,163)*y(k,59)
         mat(k,247) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,25) + (rxt(k,247) &
                       +rxt(k,253)+rxt(k,258))*y(k,46)
         mat(k,176) = (rxt(k,247)+rxt(k,253)+rxt(k,258))*y(k,40)
         mat(k,320) = rxt(k,163)*y(k,25)
         mat(k,345) = 2.000_r8*rxt(k,184)*y(k,24)
         mat(k,359) = -(rxt(k,157)*y(k,18) + (4._r8*rxt(k,158) + 4._r8*rxt(k,159) &
                      + 4._r8*rxt(k,160) + 4._r8*rxt(k,184)) * y(k,24) + rxt(k,161) &
                      *y(k,43) + rxt(k,162)*y(k,50) + rxt(k,164)*y(k,51) + rxt(k,167) &
                      *y(k,53) + (rxt(k,168) + rxt(k,169)) * y(k,69) + (rxt(k,190) &
                      + rxt(k,191) + rxt(k,192)) * y(k,2))
         mat(k,473) = -rxt(k,157)*y(k,24)
         mat(k,384) = -rxt(k,161)*y(k,24)
         mat(k,308) = -rxt(k,162)*y(k,24)
         mat(k,409) = -rxt(k,164)*y(k,24)
         mat(k,566) = -rxt(k,167)*y(k,24)
         mat(k,527) = -(rxt(k,168) + rxt(k,169)) * y(k,24)
         mat(k,496) = -(rxt(k,190) + rxt(k,191) + rxt(k,192)) * y(k,24)
         mat(k,199) = rxt(k,165)*y(k,53)
         mat(k,257) = rxt(k,183)*y(k,65)
         mat(k,384) = mat(k,384) + rxt(k,155)*y(k,59)
         mat(k,181) = rxt(k,173)*y(k,53) + rxt(k,172)*y(k,59) + rxt(k,174)*y(k,69)
         mat(k,566) = mat(k,566) + rxt(k,165)*y(k,25) + rxt(k,173)*y(k,46)
         mat(k,625) = rxt(k,156)*y(k,59)
         mat(k,333) = rxt(k,155)*y(k,43) + rxt(k,172)*y(k,46) + rxt(k,156)*y(k,55)
         mat(k,434) = rxt(k,183)*y(k,40)
         mat(k,527) = mat(k,527) + rxt(k,174)*y(k,46)
         mat(k,194) = -(rxt(k,163)*y(k,59) + rxt(k,165)*y(k,53) + rxt(k,166)*y(k,69) &
                      + (rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,40))
         mat(k,324) = -rxt(k,163)*y(k,25)
         mat(k,557) = -rxt(k,165)*y(k,25)
         mat(k,518) = -rxt(k,166)*y(k,25)
         mat(k,251) = -(rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,25)
         mat(k,350) = rxt(k,164)*y(k,51)
         mat(k,401) = rxt(k,164)*y(k,24)
         mat(k,72) = -((rxt(k,237) + rxt(k,241)) * y(k,69))
         mat(k,512) = -(rxt(k,237) + rxt(k,241)) * y(k,27)
         mat(k,267) = rxt(k,230)*y(k,52) + rxt(k,231)*y(k,53) + rxt(k,186)*y(k,58) &
                      + rxt(k,150)*y(k,59) + rxt(k,232)*y(k,69)
         mat(k,579) = rxt(k,230)*y(k,14)
         mat(k,545) = rxt(k,231)*y(k,14)
         mat(k,446) = rxt(k,186)*y(k,14)
         mat(k,321) = rxt(k,150)*y(k,14)
         mat(k,512) = mat(k,512) + rxt(k,232)*y(k,14)
         mat(k,4) = -(rxt(k,211)*y(k,65))
         mat(k,421) = -rxt(k,211)*y(k,29)
         mat(k,17) = -(rxt(k,212)*y(k,65))
         mat(k,423) = -rxt(k,212)*y(k,30)
         mat(k,154) = -(rxt(k,208)*y(k,33) + rxt(k,209)*y(k,73) + rxt(k,210)*y(k,42))
         mat(k,600) = -rxt(k,208)*y(k,31)
         mat(k,236) = -rxt(k,209)*y(k,31)
         mat(k,205) = -rxt(k,210)*y(k,31)
         mat(k,5) = 2.000_r8*rxt(k,211)*y(k,65)
         mat(k,18) = rxt(k,212)*y(k,65)
         mat(k,424) = 2.000_r8*rxt(k,211)*y(k,29) + rxt(k,212)*y(k,30)
         mat(k,288) = -((rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,43) + rxt(k,111) &
                      *y(k,54) + rxt(k,114)*y(k,55))
         mat(k,381) = -(rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,32)
         mat(k,227) = -rxt(k,111)*y(k,32)
         mat(k,622) = -rxt(k,114)*y(k,32)
         mat(k,273) = rxt(k,232)*y(k,69)
         mat(k,73) = rxt(k,241)*y(k,69)
         mat(k,157) = rxt(k,208)*y(k,33)
         mat(k,604) = rxt(k,208)*y(k,31) + rxt(k,106)*y(k,53) + rxt(k,152)*y(k,59) &
                      + rxt(k,89)*y(k,65) + rxt(k,115)*y(k,69)
         mat(k,186) = rxt(k,206)*y(k,65)
         mat(k,255) = rxt(k,183)*y(k,65)
         mat(k,169) = rxt(k,138)*y(k,69)
         mat(k,563) = rxt(k,106)*y(k,33) + rxt(k,118)*y(k,69)
         mat(k,330) = rxt(k,152)*y(k,33)
         mat(k,431) = rxt(k,89)*y(k,33) + rxt(k,206)*y(k,36) + rxt(k,183)*y(k,40)
         mat(k,524) = rxt(k,232)*y(k,14) + rxt(k,241)*y(k,27) + rxt(k,115)*y(k,33) &
                      + rxt(k,138)*y(k,47) + rxt(k,118)*y(k,53)
         mat(k,617) = -(rxt(k,89)*y(k,65) + rxt(k,106)*y(k,53) + rxt(k,115)*y(k,69) &
                      + rxt(k,152)*y(k,59) + rxt(k,208)*y(k,31))
         mat(k,444) = -rxt(k,89)*y(k,33)
         mat(k,576) = -rxt(k,106)*y(k,33)
         mat(k,537) = -rxt(k,115)*y(k,33)
         mat(k,343) = -rxt(k,152)*y(k,33)
         mat(k,161) = -rxt(k,208)*y(k,33)
         mat(k,297) = rxt(k,108)*y(k,43)
         mat(k,394) = rxt(k,108)*y(k,32)
         mat(k,83) = -(rxt(k,107)*y(k,53) + rxt(k,116)*y(k,69) + rxt(k,153)*y(k,59))
         mat(k,547) = -rxt(k,107)*y(k,35)
         mat(k,513) = -rxt(k,116)*y(k,35)
         mat(k,322) = -rxt(k,153)*y(k,35)
         mat(k,373) = 2.000_r8*rxt(k,122)*y(k,43)
         mat(k,513) = mat(k,513) + 2.000_r8*rxt(k,121)*y(k,69)
         mat(k,184) = -(rxt(k,199)*y(k,53) + rxt(k,200)*y(k,69) + (rxt(k,205) &
                      + rxt(k,206)) * y(k,65))
         mat(k,556) = -rxt(k,199)*y(k,36)
         mat(k,517) = -rxt(k,200)*y(k,36)
         mat(k,425) = -(rxt(k,205) + rxt(k,206)) * y(k,36)
         mat(k,268) = rxt(k,186)*y(k,58)
         mat(k,376) = rxt(k,187)*y(k,58)
         mat(k,447) = rxt(k,186)*y(k,14) + rxt(k,187)*y(k,43)
         mat(k,254) = -(rxt(k,170)*y(k,53) + rxt(k,171)*y(k,69) + (rxt(k,182) &
                      + rxt(k,183)) * y(k,65) + (rxt(k,247) + rxt(k,253) + rxt(k,258) &
                      ) * y(k,46) + (rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,25) &
                      + (rxt(k,254) + rxt(k,259)) * y(k,45))
         mat(k,561) = -rxt(k,170)*y(k,40)
         mat(k,522) = -rxt(k,171)*y(k,40)
         mat(k,429) = -(rxt(k,182) + rxt(k,183)) * y(k,40)
         mat(k,179) = -(rxt(k,247) + rxt(k,253) + rxt(k,258)) * y(k,40)
         mat(k,197) = -(rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,40)
         mat(k,114) = -(rxt(k,254) + rxt(k,259)) * y(k,40)
         mat(k,271) = rxt(k,150)*y(k,59)
         mat(k,354) = rxt(k,169)*y(k,69)
         mat(k,603) = rxt(k,152)*y(k,59)
         mat(k,85) = rxt(k,153)*y(k,59)
         mat(k,379) = rxt(k,154)*y(k,59)
         mat(k,179) = mat(k,179) + rxt(k,172)*y(k,59)
         mat(k,328) = rxt(k,150)*y(k,14) + rxt(k,152)*y(k,33) + rxt(k,153)*y(k,35) &
                      + rxt(k,154)*y(k,43) + rxt(k,172)*y(k,46)
         mat(k,522) = mat(k,522) + rxt(k,169)*y(k,24)
         mat(k,153) = rxt(k,208)*y(k,33) + rxt(k,210)*y(k,42) + rxt(k,209)*y(k,73)
         mat(k,599) = rxt(k,208)*y(k,31)
         mat(k,204) = rxt(k,210)*y(k,31)
         mat(k,235) = rxt(k,209)*y(k,31)
         mat(k,206) = -(rxt(k,147)*y(k,69) + rxt(k,210)*y(k,31))
         mat(k,519) = -rxt(k,147)*y(k,42)
         mat(k,155) = -rxt(k,210)*y(k,42)
         mat(k,269) = rxt(k,230)*y(k,52)
         mat(k,195) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,40)
         mat(k,252) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,25)
         mat(k,402) = rxt(k,146)*y(k,69)
         mat(k,580) = rxt(k,230)*y(k,14)
         mat(k,519) = mat(k,519) + rxt(k,146)*y(k,51)
         mat(k,385) = -((rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,32) + rxt(k,112) &
                      *y(k,53) + rxt(k,113)*y(k,55) + rxt(k,117)*y(k,69) &
                      + 4._r8*rxt(k,122)*y(k,43) + rxt(k,134)*y(k,52) + rxt(k,139) &
                      *y(k,50) + rxt(k,144)*y(k,51) + (rxt(k,154) + rxt(k,155) &
                      ) * y(k,59) + rxt(k,161)*y(k,24) + rxt(k,187)*y(k,58) + rxt(k,193) &
                      *y(k,2) + rxt(k,233)*y(k,18))
         mat(k,290) = -(rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,43)
         mat(k,567) = -rxt(k,112)*y(k,43)
         mat(k,626) = -rxt(k,113)*y(k,43)
         mat(k,528) = -rxt(k,117)*y(k,43)
         mat(k,588) = -rxt(k,134)*y(k,43)
         mat(k,309) = -rxt(k,139)*y(k,43)
         mat(k,410) = -rxt(k,144)*y(k,43)
         mat(k,334) = -(rxt(k,154) + rxt(k,155)) * y(k,43)
         mat(k,360) = -rxt(k,161)*y(k,43)
         mat(k,455) = -rxt(k,187)*y(k,43)
         mat(k,497) = -rxt(k,193)*y(k,43)
         mat(k,474) = -rxt(k,233)*y(k,43)
         mat(k,497) = mat(k,497) + rxt(k,198)*y(k,69)
         mat(k,276) = rxt(k,230)*y(k,52) + rxt(k,231)*y(k,53) + rxt(k,186)*y(k,58) &
                      + rxt(k,150)*y(k,59)
         mat(k,474) = mat(k,474) + rxt(k,157)*y(k,24) + rxt(k,234)*y(k,50)
         mat(k,360) = mat(k,360) + rxt(k,157)*y(k,18) + rxt(k,168)*y(k,69)
         mat(k,74) = rxt(k,237)*y(k,69)
         mat(k,290) = mat(k,290) + rxt(k,111)*y(k,54)
         mat(k,87) = rxt(k,107)*y(k,53) + rxt(k,153)*y(k,59) + rxt(k,116)*y(k,69)
         mat(k,309) = mat(k,309) + rxt(k,234)*y(k,18)
         mat(k,588) = mat(k,588) + rxt(k,230)*y(k,14) + rxt(k,137)*y(k,69)
         mat(k,567) = mat(k,567) + rxt(k,231)*y(k,14) + rxt(k,107)*y(k,35)
         mat(k,229) = rxt(k,111)*y(k,32)
         mat(k,626) = mat(k,626) + rxt(k,119)*y(k,69)
         mat(k,455) = mat(k,455) + rxt(k,186)*y(k,14)
         mat(k,334) = mat(k,334) + rxt(k,150)*y(k,14) + rxt(k,153)*y(k,35)
         mat(k,528) = mat(k,528) + rxt(k,198)*y(k,2) + rxt(k,168)*y(k,24) + rxt(k,237) &
                      *y(k,27) + rxt(k,116)*y(k,35) + rxt(k,137)*y(k,52) + rxt(k,119) &
                      *y(k,55)
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
         mat(k,39) = -(rxt(k,123)*y(k,69))
         mat(k,511) = -rxt(k,123)*y(k,44)
         mat(k,372) = rxt(k,144)*y(k,51)
         mat(k,397) = rxt(k,144)*y(k,43)
         mat(k,112) = -(rxt(k,201)*y(k,53) + (rxt(k,254) + rxt(k,259)) * y(k,40))
         mat(k,550) = -rxt(k,201)*y(k,45)
         mat(k,249) = -(rxt(k,254) + rxt(k,259)) * y(k,45)
         mat(k,488) = rxt(k,193)*y(k,43)
         mat(k,374) = rxt(k,193)*y(k,2)
         mat(k,177) = -(rxt(k,172)*y(k,59) + rxt(k,173)*y(k,53) + rxt(k,174)*y(k,69) &
                      + (rxt(k,247) + rxt(k,253) + rxt(k,258)) * y(k,40))
         mat(k,323) = -rxt(k,172)*y(k,46)
         mat(k,555) = -rxt(k,173)*y(k,46)
         mat(k,516) = -rxt(k,174)*y(k,46)
         mat(k,250) = -(rxt(k,247) + rxt(k,253) + rxt(k,258)) * y(k,46)
         mat(k,349) = rxt(k,161)*y(k,43)
         mat(k,193) = rxt(k,166)*y(k,69)
         mat(k,375) = rxt(k,161)*y(k,24)
         mat(k,516) = mat(k,516) + rxt(k,166)*y(k,25)
         mat(k,167) = -(rxt(k,126)*y(k,50) + (rxt(k,127) + rxt(k,128) + rxt(k,129) &
                      ) * y(k,51) + rxt(k,130)*y(k,54) + rxt(k,138)*y(k,69) + rxt(k,275) &
                      *y(k,68))
         mat(k,302) = -rxt(k,126)*y(k,47)
         mat(k,400) = -(rxt(k,127) + rxt(k,128) + rxt(k,129)) * y(k,47)
         mat(k,225) = -rxt(k,130)*y(k,47)
         mat(k,515) = -rxt(k,138)*y(k,47)
         mat(k,123) = -rxt(k,275)*y(k,47)
         mat(k,554) = rxt(k,124)*y(k,61) + rxt(k,272)*y(k,64)
         mat(k,225) = mat(k,225) + rxt(k,273)*y(k,64)
         mat(k,145) = 1.100_r8*rxt(k,268)*y(k,62) + .200_r8*rxt(k,266)*y(k,63)
         mat(k,60) = rxt(k,124)*y(k,53)
         mat(k,96) = 1.100_r8*rxt(k,268)*y(k,60)
         mat(k,131) = .200_r8*rxt(k,266)*y(k,60)
         mat(k,80) = rxt(k,272)*y(k,53) + rxt(k,273)*y(k,54)
         mat(k,396) = rxt(k,145)*y(k,52)
         mat(k,578) = rxt(k,145)*y(k,51)
         mat(k,306) = -(rxt(k,126)*y(k,47) + rxt(k,135)*y(k,52) + rxt(k,139)*y(k,43) &
                      + rxt(k,140)*y(k,55) + rxt(k,141)*y(k,53) + rxt(k,162)*y(k,24) &
                      + rxt(k,194)*y(k,2) + rxt(k,234)*y(k,18) + rxt(k,277)*y(k,68))
         mat(k,170) = -rxt(k,126)*y(k,50)
         mat(k,585) = -rxt(k,135)*y(k,50)
         mat(k,382) = -rxt(k,139)*y(k,50)
         mat(k,623) = -rxt(k,140)*y(k,50)
         mat(k,564) = -rxt(k,141)*y(k,50)
         mat(k,357) = -rxt(k,162)*y(k,50)
         mat(k,494) = -rxt(k,194)*y(k,50)
         mat(k,471) = -rxt(k,234)*y(k,50)
         mat(k,125) = -rxt(k,277)*y(k,50)
         mat(k,170) = mat(k,170) + 2.000_r8*rxt(k,128)*y(k,51) + rxt(k,130)*y(k,54) &
                      + rxt(k,138)*y(k,69)
         mat(k,407) = 2.000_r8*rxt(k,128)*y(k,47) + rxt(k,131)*y(k,53)
         mat(k,564) = mat(k,564) + rxt(k,131)*y(k,51)
         mat(k,228) = rxt(k,130)*y(k,47) + rxt(k,125)*y(k,61)
         mat(k,62) = rxt(k,125)*y(k,54)
         mat(k,525) = rxt(k,138)*y(k,47)
         mat(k,411) = -((rxt(k,127) + rxt(k,128) + rxt(k,129)) * y(k,47) + (rxt(k,131) &
                      + rxt(k,133)) * y(k,53) + rxt(k,132)*y(k,55) + rxt(k,144) &
                      *y(k,43) + rxt(k,145)*y(k,52) + rxt(k,146)*y(k,69) + rxt(k,164) &
                      *y(k,24) + rxt(k,195)*y(k,2))
         mat(k,172) = -(rxt(k,127) + rxt(k,128) + rxt(k,129)) * y(k,51)
         mat(k,568) = -(rxt(k,131) + rxt(k,133)) * y(k,51)
         mat(k,627) = -rxt(k,132)*y(k,51)
         mat(k,386) = -rxt(k,144)*y(k,51)
         mat(k,589) = -rxt(k,145)*y(k,51)
         mat(k,529) = -rxt(k,146)*y(k,51)
         mat(k,361) = -rxt(k,164)*y(k,51)
         mat(k,498) = -rxt(k,195)*y(k,51)
         mat(k,498) = mat(k,498) + rxt(k,194)*y(k,50)
         mat(k,475) = rxt(k,234)*y(k,50)
         mat(k,361) = mat(k,361) + rxt(k,162)*y(k,50)
         mat(k,386) = mat(k,386) + rxt(k,139)*y(k,50) + rxt(k,134)*y(k,52)
         mat(k,43) = rxt(k,123)*y(k,69)
         mat(k,310) = rxt(k,194)*y(k,2) + rxt(k,234)*y(k,18) + rxt(k,162)*y(k,24) &
                      + rxt(k,139)*y(k,43) + 2.000_r8*rxt(k,135)*y(k,52) + rxt(k,141) &
                      *y(k,53) + rxt(k,140)*y(k,55)
         mat(k,589) = mat(k,589) + rxt(k,134)*y(k,43) + 2.000_r8*rxt(k,135)*y(k,50) &
                      + rxt(k,136)*y(k,53) + rxt(k,137)*y(k,69)
         mat(k,568) = mat(k,568) + rxt(k,141)*y(k,50) + rxt(k,136)*y(k,52)
         mat(k,627) = mat(k,627) + rxt(k,140)*y(k,50)
         mat(k,529) = mat(k,529) + rxt(k,123)*y(k,44) + rxt(k,137)*y(k,52)
         mat(k,596) = -(rxt(k,134)*y(k,43) + rxt(k,135)*y(k,50) + rxt(k,136)*y(k,53) &
                      + rxt(k,137)*y(k,69) + rxt(k,145)*y(k,51) + rxt(k,230)*y(k,14))
         mat(k,393) = -rxt(k,134)*y(k,52)
         mat(k,317) = -rxt(k,135)*y(k,52)
         mat(k,575) = -rxt(k,136)*y(k,52)
         mat(k,536) = -rxt(k,137)*y(k,52)
         mat(k,418) = -rxt(k,145)*y(k,52)
         mat(k,284) = -rxt(k,230)*y(k,52)
         mat(k,71) = rxt(k,196)*y(k,53)
         mat(k,203) = rxt(k,165)*y(k,53) + rxt(k,163)*y(k,59) + rxt(k,166)*y(k,69)
         mat(k,160) = rxt(k,210)*y(k,42)
         mat(k,212) = rxt(k,210)*y(k,31) + rxt(k,147)*y(k,69)
         mat(k,418) = mat(k,418) + rxt(k,133)*y(k,53) + rxt(k,132)*y(k,55)
         mat(k,575) = mat(k,575) + rxt(k,196)*y(k,3) + rxt(k,165)*y(k,25) + rxt(k,133) &
                      *y(k,51)
         mat(k,634) = rxt(k,132)*y(k,51)
         mat(k,342) = rxt(k,163)*y(k,25)
         mat(k,536) = mat(k,536) + rxt(k,166)*y(k,25) + rxt(k,147)*y(k,42)
         mat(k,574) = -(rxt(k,103)*y(k,55) + 4._r8*rxt(k,104)*y(k,53) + rxt(k,105) &
                      *y(k,54) + rxt(k,106)*y(k,33) + rxt(k,107)*y(k,35) + rxt(k,112) &
                      *y(k,43) + rxt(k,118)*y(k,69) + (rxt(k,131) + rxt(k,133) &
                      ) * y(k,51) + rxt(k,136)*y(k,52) + rxt(k,141)*y(k,50) + rxt(k,165) &
                      *y(k,25) + rxt(k,167)*y(k,24) + rxt(k,170)*y(k,40) + rxt(k,173) &
                      *y(k,46) + rxt(k,196)*y(k,3) + rxt(k,197)*y(k,2) + rxt(k,199) &
                      *y(k,36) + rxt(k,201)*y(k,45) + rxt(k,231)*y(k,14) + (rxt(k,270) &
                      + rxt(k,271)) * y(k,62) + rxt(k,272)*y(k,64))
         mat(k,633) = -rxt(k,103)*y(k,53)
         mat(k,233) = -rxt(k,105)*y(k,53)
         mat(k,615) = -rxt(k,106)*y(k,53)
         mat(k,89) = -rxt(k,107)*y(k,53)
         mat(k,392) = -rxt(k,112)*y(k,53)
         mat(k,535) = -rxt(k,118)*y(k,53)
         mat(k,417) = -(rxt(k,131) + rxt(k,133)) * y(k,53)
         mat(k,595) = -rxt(k,136)*y(k,53)
         mat(k,316) = -rxt(k,141)*y(k,53)
         mat(k,202) = -rxt(k,165)*y(k,53)
         mat(k,367) = -rxt(k,167)*y(k,53)
         mat(k,264) = -rxt(k,170)*y(k,53)
         mat(k,183) = -rxt(k,173)*y(k,53)
         mat(k,70) = -rxt(k,196)*y(k,53)
         mat(k,504) = -rxt(k,197)*y(k,53)
         mat(k,191) = -rxt(k,199)*y(k,53)
         mat(k,119) = -rxt(k,201)*y(k,53)
         mat(k,283) = -rxt(k,231)*y(k,53)
         mat(k,100) = -(rxt(k,270) + rxt(k,271)) * y(k,53)
         mat(k,82) = -rxt(k,272)*y(k,53)
         mat(k,295) = rxt(k,110)*y(k,43)
         mat(k,392) = mat(k,392) + rxt(k,110)*y(k,32)
         mat(k,175) = rxt(k,126)*y(k,50) + rxt(k,127)*y(k,51) + rxt(k,130)*y(k,54) &
                      + rxt(k,275)*y(k,68)
         mat(k,316) = mat(k,316) + rxt(k,126)*y(k,47)
         mat(k,417) = mat(k,417) + rxt(k,127)*y(k,47)
         mat(k,233) = mat(k,233) + rxt(k,130)*y(k,47) + rxt(k,274)*y(k,64) + ( &
                      + rxt(k,92)+rxt(k,93))*y(k,65) + rxt(k,280)*y(k,70) + rxt(k,284) &
                      *y(k,71)
         mat(k,152) = rxt(k,266)*y(k,63) + 1.150_r8*rxt(k,267)*y(k,68)
         mat(k,135) = rxt(k,266)*y(k,60)
         mat(k,82) = mat(k,82) + rxt(k,274)*y(k,54)
         mat(k,442) = (rxt(k,92)+rxt(k,93))*y(k,54)
         mat(k,127) = rxt(k,275)*y(k,47) + 1.150_r8*rxt(k,267)*y(k,60)
         mat(k,535) = mat(k,535) + 2.000_r8*rxt(k,120)*y(k,69)
         mat(k,110) = rxt(k,280)*y(k,54)
         mat(k,58) = rxt(k,284)*y(k,54)
         mat(k,226) = -(rxt(k,92)*y(k,65) + rxt(k,97)*y(k,66) + rxt(k,105)*y(k,53) &
                      + rxt(k,111)*y(k,32) + rxt(k,125)*y(k,61) + rxt(k,130)*y(k,47) &
                      + rxt(k,269)*y(k,62) + (rxt(k,273) + rxt(k,274)) * y(k,64) &
                      + rxt(k,280)*y(k,70) + rxt(k,284)*y(k,71))
         mat(k,427) = -rxt(k,92)*y(k,54)
         mat(k,11) = -rxt(k,97)*y(k,54)
         mat(k,559) = -rxt(k,105)*y(k,54)
         mat(k,286) = -rxt(k,111)*y(k,54)
         mat(k,61) = -rxt(k,125)*y(k,54)
         mat(k,168) = -rxt(k,130)*y(k,54)
         mat(k,97) = -rxt(k,269)*y(k,54)
         mat(k,81) = -(rxt(k,273) + rxt(k,274)) * y(k,54)
         mat(k,106) = -rxt(k,280)*y(k,54)
         mat(k,57) = -rxt(k,284)*y(k,54)
         mat(k,490) = 2.000_r8*rxt(k,189)*y(k,2) + (rxt(k,191)+rxt(k,192))*y(k,24) &
                      + rxt(k,193)*y(k,43) + rxt(k,197)*y(k,53)
         mat(k,467) = rxt(k,233)*y(k,43)
         mat(k,352) = (rxt(k,191)+rxt(k,192))*y(k,2) + (2.000_r8*rxt(k,158) &
                       +2.000_r8*rxt(k,159))*y(k,24) + rxt(k,161)*y(k,43) + rxt(k,167) &
                      *y(k,53) + rxt(k,169)*y(k,69)
         mat(k,286) = mat(k,286) + rxt(k,108)*y(k,43) + rxt(k,114)*y(k,55)
         mat(k,377) = rxt(k,193)*y(k,2) + rxt(k,233)*y(k,18) + rxt(k,161)*y(k,24) &
                      + rxt(k,108)*y(k,32) + 2.000_r8*rxt(k,122)*y(k,43) + rxt(k,134) &
                      *y(k,52) + rxt(k,112)*y(k,53) + 2.000_r8*rxt(k,113)*y(k,55) &
                      + rxt(k,187)*y(k,58) + rxt(k,154)*y(k,59) + rxt(k,117)*y(k,69)
         mat(k,40) = rxt(k,123)*y(k,69)
         mat(k,168) = mat(k,168) + rxt(k,129)*y(k,51)
         mat(k,303) = rxt(k,140)*y(k,55) + rxt(k,277)*y(k,68)
         mat(k,403) = rxt(k,129)*y(k,47) + rxt(k,131)*y(k,53) + rxt(k,132)*y(k,55)
         mat(k,581) = rxt(k,134)*y(k,43) + rxt(k,136)*y(k,53)
         mat(k,559) = mat(k,559) + rxt(k,197)*y(k,2) + rxt(k,167)*y(k,24) + rxt(k,112) &
                      *y(k,43) + rxt(k,131)*y(k,51) + rxt(k,136)*y(k,52) &
                      + 2.000_r8*rxt(k,104)*y(k,53) + 2.000_r8*rxt(k,103)*y(k,55) &
                      + rxt(k,96)*y(k,66) + rxt(k,118)*y(k,69)
         mat(k,226) = mat(k,226) + 2.000_r8*rxt(k,97)*y(k,66)
         mat(k,621) = rxt(k,114)*y(k,32) + 2.000_r8*rxt(k,113)*y(k,43) + rxt(k,140) &
                      *y(k,50) + rxt(k,132)*y(k,51) + 2.000_r8*rxt(k,103)*y(k,53) &
                      + rxt(k,188)*y(k,58) + rxt(k,156)*y(k,59) + 2.000_r8*rxt(k,94) &
                      *y(k,65) + rxt(k,119)*y(k,69)
         mat(k,448) = rxt(k,187)*y(k,43) + rxt(k,188)*y(k,55)
         mat(k,326) = rxt(k,154)*y(k,43) + rxt(k,156)*y(k,55)
         mat(k,427) = mat(k,427) + 2.000_r8*rxt(k,94)*y(k,55)
         mat(k,11) = mat(k,11) + rxt(k,96)*y(k,53) + 2.000_r8*rxt(k,97)*y(k,54)
         mat(k,124) = rxt(k,277)*y(k,50)
         mat(k,520) = rxt(k,169)*y(k,24) + rxt(k,117)*y(k,43) + rxt(k,123)*y(k,44) &
                      + rxt(k,118)*y(k,53) + rxt(k,119)*y(k,55)
         mat(k,636) = -(rxt(k,94)*y(k,65) + rxt(k,103)*y(k,53) + rxt(k,113)*y(k,43) &
                      + rxt(k,114)*y(k,32) + rxt(k,119)*y(k,69) + rxt(k,132)*y(k,51) &
                      + rxt(k,140)*y(k,50) + rxt(k,156)*y(k,59) + rxt(k,188)*y(k,58))
         mat(k,445) = -rxt(k,94)*y(k,55)
         mat(k,577) = -rxt(k,103)*y(k,55)
         mat(k,395) = -rxt(k,113)*y(k,55)
         mat(k,298) = -rxt(k,114)*y(k,55)
         mat(k,538) = -rxt(k,119)*y(k,55)
         mat(k,420) = -rxt(k,132)*y(k,55)
         mat(k,319) = -rxt(k,140)*y(k,55)
         mat(k,344) = -rxt(k,156)*y(k,55)
         mat(k,465) = -rxt(k,188)*y(k,55)
         mat(k,577) = mat(k,577) + rxt(k,105)*y(k,54)
         mat(k,234) = rxt(k,105)*y(k,53)
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
         mat(k,485) = rxt(k,190)*y(k,24)
         mat(k,346) = rxt(k,190)*y(k,2) + 2.000_r8*rxt(k,160)*y(k,24)
         mat(k,458) = -(rxt(k,186)*y(k,14) + rxt(k,187)*y(k,43) + rxt(k,188)*y(k,55))
         mat(k,279) = -rxt(k,186)*y(k,58)
         mat(k,388) = -rxt(k,187)*y(k,58)
         mat(k,629) = -rxt(k,188)*y(k,58)
         mat(k,500) = 4.000_r8*rxt(k,189)*y(k,2) + (rxt(k,190)+rxt(k,191))*y(k,24) &
                      + rxt(k,194)*y(k,50) + rxt(k,197)*y(k,53) + rxt(k,198)*y(k,69)
         mat(k,363) = (rxt(k,190)+rxt(k,191))*y(k,2)
         mat(k,188) = rxt(k,199)*y(k,53) + rxt(k,205)*y(k,65) + rxt(k,200)*y(k,69)
         mat(k,312) = rxt(k,194)*y(k,2)
         mat(k,570) = rxt(k,197)*y(k,2) + rxt(k,199)*y(k,36)
         mat(k,438) = rxt(k,205)*y(k,36)
         mat(k,531) = rxt(k,198)*y(k,2) + rxt(k,200)*y(k,36)
         mat(k,332) = -(rxt(k,150)*y(k,14) + rxt(k,152)*y(k,33) + rxt(k,153)*y(k,35) &
                      + (rxt(k,154) + rxt(k,155)) * y(k,43) + rxt(k,156)*y(k,55) &
                      + rxt(k,163)*y(k,25) + rxt(k,172)*y(k,46))
         mat(k,274) = -rxt(k,150)*y(k,59)
         mat(k,606) = -rxt(k,152)*y(k,59)
         mat(k,86) = -rxt(k,153)*y(k,59)
         mat(k,383) = -(rxt(k,154) + rxt(k,155)) * y(k,59)
         mat(k,624) = -rxt(k,156)*y(k,59)
         mat(k,198) = -rxt(k,163)*y(k,59)
         mat(k,180) = -rxt(k,172)*y(k,59)
         mat(k,495) = rxt(k,191)*y(k,24)
         mat(k,472) = rxt(k,157)*y(k,24)
         mat(k,358) = rxt(k,191)*y(k,2) + rxt(k,157)*y(k,18) + (4.000_r8*rxt(k,158) &
                       +2.000_r8*rxt(k,160))*y(k,24) + rxt(k,162)*y(k,50) + rxt(k,167) &
                      *y(k,53) + rxt(k,168)*y(k,69)
         mat(k,19) = rxt(k,212)*y(k,65)
         mat(k,256) = rxt(k,170)*y(k,53) + rxt(k,182)*y(k,65) + rxt(k,171)*y(k,69)
         mat(k,307) = rxt(k,162)*y(k,24)
         mat(k,565) = rxt(k,167)*y(k,24) + rxt(k,170)*y(k,40)
         mat(k,433) = rxt(k,212)*y(k,30) + rxt(k,182)*y(k,40)
         mat(k,526) = rxt(k,168)*y(k,24) + rxt(k,171)*y(k,40)
         mat(k,144) = -(rxt(k,266)*y(k,63) + rxt(k,267)*y(k,68) + rxt(k,268)*y(k,62))
         mat(k,130) = -rxt(k,266)*y(k,60)
         mat(k,122) = -rxt(k,267)*y(k,60)
         mat(k,95) = -rxt(k,268)*y(k,60)
         mat(k,59) = -(rxt(k,124)*y(k,53) + rxt(k,125)*y(k,54))
         mat(k,543) = -rxt(k,124)*y(k,61)
         mat(k,218) = -rxt(k,125)*y(k,61)
         mat(k,543) = mat(k,543) + rxt(k,270)*y(k,62)
         mat(k,138) = .900_r8*rxt(k,268)*y(k,62) + .800_r8*rxt(k,266)*y(k,63)
         mat(k,90) = rxt(k,270)*y(k,53) + .900_r8*rxt(k,268)*y(k,60)
         mat(k,128) = .800_r8*rxt(k,266)*y(k,60)
         mat(k,91) = -(rxt(k,268)*y(k,60) + rxt(k,269)*y(k,54) + (rxt(k,270) + rxt(k,271) &
                      ) * y(k,53))
         mat(k,140) = -rxt(k,268)*y(k,62)
         mat(k,220) = -rxt(k,269)*y(k,62)
         mat(k,548) = -(rxt(k,270) + rxt(k,271)) * y(k,62)
         mat(k,129) = -(rxt(k,266)*y(k,60))
         mat(k,143) = -rxt(k,266)*y(k,63)
         mat(k,165) = rxt(k,275)*y(k,68)
         mat(k,300) = rxt(k,277)*y(k,68)
         mat(k,552) = rxt(k,270)*y(k,62)
         mat(k,223) = rxt(k,274)*y(k,64)
         mat(k,94) = rxt(k,270)*y(k,53)
         mat(k,79) = rxt(k,274)*y(k,54)
         mat(k,121) = rxt(k,275)*y(k,47) + rxt(k,277)*y(k,50)
         mat(k,76) = -(rxt(k,272)*y(k,53) + (rxt(k,273) + rxt(k,274)) * y(k,54))
         mat(k,546) = -rxt(k,272)*y(k,64)
         mat(k,219) = -(rxt(k,273) + rxt(k,274)) * y(k,64)
         mat(k,437) = -(rxt(k,89)*y(k,33) + rxt(k,90)*y(k,73) + (rxt(k,92) + rxt(k,93) &
                      ) * y(k,54) + rxt(k,94)*y(k,55) + (rxt(k,182) + rxt(k,183) &
                      ) * y(k,40) + (rxt(k,205) + rxt(k,206)) * y(k,36) + rxt(k,211) &
                      *y(k,29) + rxt(k,212)*y(k,30))
         mat(k,610) = -rxt(k,89)*y(k,65)
         mat(k,241) = -rxt(k,90)*y(k,65)
         mat(k,231) = -(rxt(k,92) + rxt(k,93)) * y(k,65)
         mat(k,628) = -rxt(k,94)*y(k,65)
         mat(k,259) = -(rxt(k,182) + rxt(k,183)) * y(k,65)
         mat(k,187) = -(rxt(k,205) + rxt(k,206)) * y(k,65)
         mat(k,6) = -rxt(k,211)*y(k,65)
         mat(k,20) = -rxt(k,212)*y(k,65)
         mat(k,231) = mat(k,231) + rxt(k,125)*y(k,61)
         mat(k,150) = .850_r8*rxt(k,267)*y(k,68)
         mat(k,63) = rxt(k,125)*y(k,54)
         mat(k,126) = .850_r8*rxt(k,267)*y(k,60)
         mat(k,10) = -(rxt(k,96)*y(k,53) + rxt(k,97)*y(k,54))
         mat(k,539) = -rxt(k,96)*y(k,66)
         mat(k,214) = -rxt(k,97)*y(k,66)
         mat(k,539) = mat(k,539) + rxt(k,100)*y(k,67)
         mat(k,214) = mat(k,214) + rxt(k,101)*y(k,67)
         mat(k,619) = rxt(k,102)*y(k,67)
         mat(k,12) = rxt(k,100)*y(k,53) + rxt(k,101)*y(k,54) + rxt(k,102)*y(k,55)
         mat(k,13) = -(rxt(k,100)*y(k,53) + rxt(k,101)*y(k,54) + rxt(k,102)*y(k,55))
         mat(k,540) = -rxt(k,100)*y(k,67)
         mat(k,215) = -rxt(k,101)*y(k,67)
         mat(k,620) = -rxt(k,102)*y(k,67)
         mat(k,215) = mat(k,215) + rxt(k,92)*y(k,65)
         mat(k,422) = rxt(k,92)*y(k,54)
         mat(k,120) = -(rxt(k,267)*y(k,60) + rxt(k,275)*y(k,47) + rxt(k,277)*y(k,50))
         mat(k,142) = -rxt(k,267)*y(k,68)
         mat(k,164) = -rxt(k,275)*y(k,68)
         mat(k,299) = -rxt(k,277)*y(k,68)
         mat(k,222) = rxt(k,269)*y(k,62) + rxt(k,273)*y(k,64) + rxt(k,280)*y(k,70) &
                      + rxt(k,284)*y(k,71)
         mat(k,93) = rxt(k,269)*y(k,54)
         mat(k,78) = rxt(k,273)*y(k,54)
         mat(k,103) = rxt(k,280)*y(k,54)
         mat(k,56) = rxt(k,284)*y(k,54)
         mat(k,534) = -(rxt(k,115)*y(k,33) + rxt(k,116)*y(k,35) + rxt(k,117)*y(k,43) &
                      + rxt(k,118)*y(k,53) + rxt(k,119)*y(k,55) + (4._r8*rxt(k,120) &
                      + 4._r8*rxt(k,121)) * y(k,69) + rxt(k,123)*y(k,44) + rxt(k,137) &
                      *y(k,52) + rxt(k,138)*y(k,47) + rxt(k,146)*y(k,51) + rxt(k,147) &
                      *y(k,42) + rxt(k,166)*y(k,25) + (rxt(k,168) + rxt(k,169) &
                      ) * y(k,24) + rxt(k,171)*y(k,40) + rxt(k,174)*y(k,46) + rxt(k,198) &
                      *y(k,2) + rxt(k,200)*y(k,36) + rxt(k,232)*y(k,14) + rxt(k,235) &
                      *y(k,19) + (rxt(k,237) + rxt(k,241)) * y(k,27))
         mat(k,614) = -rxt(k,115)*y(k,69)
         mat(k,88) = -rxt(k,116)*y(k,69)
         mat(k,391) = -rxt(k,117)*y(k,69)
         mat(k,573) = -rxt(k,118)*y(k,69)
         mat(k,632) = -rxt(k,119)*y(k,69)
         mat(k,44) = -rxt(k,123)*y(k,69)
         mat(k,594) = -rxt(k,137)*y(k,69)
         mat(k,174) = -rxt(k,138)*y(k,69)
         mat(k,416) = -rxt(k,146)*y(k,69)
         mat(k,211) = -rxt(k,147)*y(k,69)
         mat(k,201) = -rxt(k,166)*y(k,69)
         mat(k,366) = -(rxt(k,168) + rxt(k,169)) * y(k,69)
         mat(k,263) = -rxt(k,171)*y(k,69)
         mat(k,182) = -rxt(k,174)*y(k,69)
         mat(k,503) = -rxt(k,198)*y(k,69)
         mat(k,190) = -rxt(k,200)*y(k,69)
         mat(k,282) = -rxt(k,232)*y(k,69)
         mat(k,38) = -rxt(k,235)*y(k,69)
         mat(k,75) = -(rxt(k,237) + rxt(k,241)) * y(k,69)
         mat(k,282) = mat(k,282) + rxt(k,231)*y(k,53)
         mat(k,38) = mat(k,38) + .300_r8*rxt(k,235)*y(k,69)
         mat(k,159) = rxt(k,209)*y(k,73)
         mat(k,294) = 2.000_r8*rxt(k,109)*y(k,43) + rxt(k,114)*y(k,55)
         mat(k,614) = mat(k,614) + rxt(k,106)*y(k,53) + rxt(k,89)*y(k,65)
         mat(k,88) = mat(k,88) + rxt(k,107)*y(k,53)
         mat(k,190) = mat(k,190) + rxt(k,199)*y(k,53) + rxt(k,205)*y(k,65)
         mat(k,263) = mat(k,263) + rxt(k,170)*y(k,53) + rxt(k,182)*y(k,65)
         mat(k,391) = mat(k,391) + 2.000_r8*rxt(k,109)*y(k,32) + rxt(k,139)*y(k,50) &
                      + rxt(k,134)*y(k,52) + rxt(k,112)*y(k,53) + rxt(k,113)*y(k,55) &
                      + rxt(k,155)*y(k,59)
         mat(k,118) = rxt(k,201)*y(k,53)
         mat(k,182) = mat(k,182) + rxt(k,173)*y(k,53)
         mat(k,315) = rxt(k,139)*y(k,43)
         mat(k,594) = mat(k,594) + rxt(k,134)*y(k,43)
         mat(k,573) = mat(k,573) + rxt(k,231)*y(k,14) + rxt(k,106)*y(k,33) &
                      + rxt(k,107)*y(k,35) + rxt(k,199)*y(k,36) + rxt(k,170)*y(k,40) &
                      + rxt(k,112)*y(k,43) + rxt(k,201)*y(k,45) + rxt(k,173)*y(k,46)
         mat(k,632) = mat(k,632) + rxt(k,114)*y(k,32) + rxt(k,113)*y(k,43)
         mat(k,340) = rxt(k,155)*y(k,43)
         mat(k,441) = rxt(k,89)*y(k,33) + rxt(k,205)*y(k,36) + rxt(k,182)*y(k,40) &
                      + 2.000_r8*rxt(k,90)*y(k,73)
         mat(k,534) = mat(k,534) + .300_r8*rxt(k,235)*y(k,19)
         mat(k,243) = rxt(k,209)*y(k,31) + 2.000_r8*rxt(k,90)*y(k,65)
         mat(k,102) = -(rxt(k,280)*y(k,54))
         mat(k,221) = -rxt(k,280)*y(k,70)
         mat(k,549) = rxt(k,271)*y(k,62) + rxt(k,272)*y(k,64) + rxt(k,283)*y(k,71) &
                      + rxt(k,289)*y(k,72)
         mat(k,141) = rxt(k,281)*y(k,71) + rxt(k,286)*y(k,72)
         mat(k,92) = rxt(k,271)*y(k,53)
         mat(k,77) = rxt(k,272)*y(k,53)
         mat(k,55) = rxt(k,283)*y(k,53) + rxt(k,281)*y(k,60)
         mat(k,50) = rxt(k,289)*y(k,53) + rxt(k,286)*y(k,60)
         mat(k,53) = -(rxt(k,281)*y(k,60) + rxt(k,283)*y(k,53) + rxt(k,284)*y(k,54))
         mat(k,137) = -rxt(k,281)*y(k,71)
         mat(k,542) = -rxt(k,283)*y(k,71)
         mat(k,217) = -rxt(k,284)*y(k,71)
         mat(k,137) = mat(k,137) + rxt(k,285)*y(k,72)
         mat(k,47) = rxt(k,285)*y(k,60)
         mat(k,46) = -((rxt(k,285) + rxt(k,286)) * y(k,60) + rxt(k,289)*y(k,53))
         mat(k,136) = -(rxt(k,285) + rxt(k,286)) * y(k,72)
         mat(k,541) = -rxt(k,289)*y(k,72)
         mat(k,238) = -(rxt(k,90)*y(k,65) + rxt(k,209)*y(k,31))
         mat(k,428) = -rxt(k,90)*y(k,73)
         mat(k,156) = -rxt(k,209)*y(k,73)
         mat(k,270) = rxt(k,232)*y(k,69)
         mat(k,34) = rxt(k,235)*y(k,69)
         mat(k,287) = rxt(k,110)*y(k,43)
         mat(k,602) = rxt(k,115)*y(k,69)
         mat(k,84) = rxt(k,116)*y(k,69)
         mat(k,185) = rxt(k,200)*y(k,69)
         mat(k,253) = (rxt(k,254)+rxt(k,259))*y(k,45) + (rxt(k,247)+rxt(k,253) &
                       +rxt(k,258))*y(k,46) + rxt(k,171)*y(k,69)
         mat(k,207) = rxt(k,147)*y(k,69)
         mat(k,378) = rxt(k,110)*y(k,32) + rxt(k,117)*y(k,69)
         mat(k,41) = rxt(k,123)*y(k,69)
         mat(k,113) = (rxt(k,254)+rxt(k,259))*y(k,40)
         mat(k,178) = (rxt(k,247)+rxt(k,253)+rxt(k,258))*y(k,40) + rxt(k,174)*y(k,69)
         mat(k,521) = rxt(k,232)*y(k,14) + rxt(k,235)*y(k,19) + rxt(k,115)*y(k,33) &
                      + rxt(k,116)*y(k,35) + rxt(k,200)*y(k,36) + rxt(k,171)*y(k,40) &
                      + rxt(k,147)*y(k,42) + rxt(k,117)*y(k,43) + rxt(k,123)*y(k,44) &
                      + rxt(k,174)*y(k,46) + 2.000_r8*rxt(k,120)*y(k,69)
      end do
      end subroutine nlnmat03
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
         mat(k, 4) = mat(k, 4) + lmat(k, 4)
         mat(k, 5) = mat(k, 5) + lmat(k, 5)
         mat(k, 7) = lmat(k, 7)
         mat(k, 8) = lmat(k, 8)
         mat(k, 9) = lmat(k, 9)
         mat(k, 10) = mat(k, 10) + lmat(k, 10)
         mat(k, 11) = mat(k, 11) + lmat(k, 11)
         mat(k, 12) = mat(k, 12) + lmat(k, 12)
         mat(k, 13) = mat(k, 13) + lmat(k, 13)
         mat(k, 14) = lmat(k, 14)
         mat(k, 15) = lmat(k, 15)
         mat(k, 16) = lmat(k, 16)
         mat(k, 17) = mat(k, 17) + lmat(k, 17)
         mat(k, 18) = mat(k, 18) + lmat(k, 18)
         mat(k, 19) = mat(k, 19) + lmat(k, 19)
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
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 35) = mat(k, 35) + lmat(k, 35)
         mat(k, 36) = lmat(k, 36)
         mat(k, 38) = mat(k, 38) + lmat(k, 38)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 42) = lmat(k, 42)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 44) = mat(k, 44) + lmat(k, 44)
         mat(k, 45) = lmat(k, 45)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 47) = mat(k, 47) + lmat(k, 47)
         mat(k, 48) = lmat(k, 48)
         mat(k, 49) = lmat(k, 49)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 51) = lmat(k, 51)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 54) = lmat(k, 54)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 101) = lmat(k, 101)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = lmat(k, 104)
         mat(k, 105) = lmat(k, 105)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 116) = lmat(k, 116)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 129) = mat(k, 129) + lmat(k, 129)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 158) = lmat(k, 158)
         mat(k, 162) = lmat(k, 162)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 198) = mat(k, 198) + lmat(k, 198)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 209) = lmat(k, 209)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 216) = lmat(k, 216)
         mat(k, 217) = mat(k, 217) + lmat(k, 217)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 224) = lmat(k, 224)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 239) = lmat(k, 239)
         mat(k, 241) = mat(k, 241) + lmat(k, 241)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 246) = lmat(k, 246)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 273) = mat(k, 273) + lmat(k, 273)
         mat(k, 285) = lmat(k, 285)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 301) = lmat(k, 301)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 306) = mat(k, 306) + lmat(k, 306)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 332) = mat(k, 332) + lmat(k, 332)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 359) = mat(k, 359) + lmat(k, 359)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 402) = mat(k, 402) + lmat(k, 402)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 427) = mat(k, 427) + lmat(k, 427)
         mat(k, 430) = lmat(k, 430)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 432) = lmat(k, 432)
         mat(k, 433) = mat(k, 433) + lmat(k, 433)
         mat(k, 435) = lmat(k, 435)
         mat(k, 437) = mat(k, 437) + lmat(k, 437)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 439) = lmat(k, 439)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 444) = mat(k, 444) + lmat(k, 444)
         mat(k, 458) = mat(k, 458) + lmat(k, 458)
         mat(k, 478) = mat(k, 478) + lmat(k, 478)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 502) = mat(k, 502) + lmat(k, 502)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 508) = lmat(k, 508)
         mat(k, 509) = lmat(k, 509)
         mat(k, 521) = mat(k, 521) + lmat(k, 521)
         mat(k, 526) = mat(k, 526) + lmat(k, 526)
         mat(k, 528) = mat(k, 528) + lmat(k, 528)
         mat(k, 531) = mat(k, 531) + lmat(k, 531)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 549) = mat(k, 549) + lmat(k, 549)
         mat(k, 553) = lmat(k, 553)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 580) = mat(k, 580) + lmat(k, 580)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 585) = mat(k, 585) + lmat(k, 585)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 595) = mat(k, 595) + lmat(k, 595)
         mat(k, 596) = mat(k, 596) + lmat(k, 596)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 619) = mat(k, 619) + lmat(k, 619)
         mat(k, 621) = mat(k, 621) + lmat(k, 621)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 636) = mat(k, 636) + lmat(k, 636)
         mat(k, 98) = 0._r8
         mat(k, 99) = 0._r8
         mat(k, 107) = 0._r8
         mat(k, 108) = 0._r8
         mat(k, 109) = 0._r8
         mat(k, 115) = 0._r8
         mat(k, 132) = 0._r8
         mat(k, 133) = 0._r8
         mat(k, 134) = 0._r8
         mat(k, 139) = 0._r8
         mat(k, 146) = 0._r8
         mat(k, 147) = 0._r8
         mat(k, 148) = 0._r8
         mat(k, 149) = 0._r8
         mat(k, 151) = 0._r8
         mat(k, 163) = 0._r8
         mat(k, 171) = 0._r8
         mat(k, 173) = 0._r8
         mat(k, 196) = 0._r8
         mat(k, 208) = 0._r8
         mat(k, 210) = 0._r8
         mat(k, 213) = 0._r8
         mat(k, 230) = 0._r8
         mat(k, 232) = 0._r8
         mat(k, 237) = 0._r8
         mat(k, 240) = 0._r8
         mat(k, 242) = 0._r8
         mat(k, 245) = 0._r8
         mat(k, 258) = 0._r8
         mat(k, 260) = 0._r8
         mat(k, 261) = 0._r8
         mat(k, 262) = 0._r8
         mat(k, 265) = 0._r8
         mat(k, 266) = 0._r8
         mat(k, 275) = 0._r8
         mat(k, 277) = 0._r8
         mat(k, 278) = 0._r8
         mat(k, 280) = 0._r8
         mat(k, 281) = 0._r8
         mat(k, 289) = 0._r8
         mat(k, 291) = 0._r8
         mat(k, 292) = 0._r8
         mat(k, 293) = 0._r8
         mat(k, 296) = 0._r8
         mat(k, 305) = 0._r8
         mat(k, 311) = 0._r8
         mat(k, 318) = 0._r8
         mat(k, 325) = 0._r8
         mat(k, 327) = 0._r8
         mat(k, 331) = 0._r8
         mat(k, 335) = 0._r8
         mat(k, 336) = 0._r8
         mat(k, 339) = 0._r8
         mat(k, 341) = 0._r8
         mat(k, 351) = 0._r8
         mat(k, 353) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 362) = 0._r8
         mat(k, 368) = 0._r8
         mat(k, 369) = 0._r8
         mat(k, 370) = 0._r8
         mat(k, 380) = 0._r8
         mat(k, 387) = 0._r8
         mat(k, 399) = 0._r8
         mat(k, 404) = 0._r8
         mat(k, 405) = 0._r8
         mat(k, 406) = 0._r8
         mat(k, 408) = 0._r8
         mat(k, 412) = 0._r8
         mat(k, 413) = 0._r8
         mat(k, 414) = 0._r8
         mat(k, 419) = 0._r8
         mat(k, 426) = 0._r8
         mat(k, 436) = 0._r8
         mat(k, 443) = 0._r8
         mat(k, 449) = 0._r8
         mat(k, 451) = 0._r8
         mat(k, 452) = 0._r8
         mat(k, 453) = 0._r8
         mat(k, 454) = 0._r8
         mat(k, 456) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 459) = 0._r8
         mat(k, 461) = 0._r8
         mat(k, 462) = 0._r8
         mat(k, 463) = 0._r8
         mat(k, 464) = 0._r8
         mat(k, 468) = 0._r8
         mat(k, 470) = 0._r8
         mat(k, 476) = 0._r8
         mat(k, 477) = 0._r8
         mat(k, 479) = 0._r8
         mat(k, 480) = 0._r8
         mat(k, 481) = 0._r8
         mat(k, 482) = 0._r8
         mat(k, 483) = 0._r8
         mat(k, 484) = 0._r8
         mat(k, 489) = 0._r8
         mat(k, 491) = 0._r8
         mat(k, 492) = 0._r8
         mat(k, 493) = 0._r8
         mat(k, 499) = 0._r8
         mat(k, 501) = 0._r8
         mat(k, 505) = 0._r8
         mat(k, 506) = 0._r8
         mat(k, 507) = 0._r8
         mat(k, 514) = 0._r8
         mat(k, 530) = 0._r8
         mat(k, 551) = 0._r8
         mat(k, 558) = 0._r8
         mat(k, 560) = 0._r8
         mat(k, 569) = 0._r8
         mat(k, 571) = 0._r8
         mat(k, 582) = 0._r8
         mat(k, 584) = 0._r8
         mat(k, 586) = 0._r8
         mat(k, 587) = 0._r8
         mat(k, 590) = 0._r8
         mat(k, 591) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 598) = 0._r8
         mat(k, 601) = 0._r8
         mat(k, 605) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 608) = 0._r8
         mat(k, 609) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 612) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 616) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 630) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 4) = mat(k, 4) - dti(k)
         mat(k, 7) = mat(k, 7) - dti(k)
         mat(k, 10) = mat(k, 10) - dti(k)
         mat(k, 13) = mat(k, 13) - dti(k)
         mat(k, 15) = mat(k, 15) - dti(k)
         mat(k, 17) = mat(k, 17) - dti(k)
         mat(k, 21) = mat(k, 21) - dti(k)
         mat(k, 24) = mat(k, 24) - dti(k)
         mat(k, 27) = mat(k, 27) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 46) = mat(k, 46) - dti(k)
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 59) = mat(k, 59) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 72) = mat(k, 72) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 91) = mat(k, 91) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 120) = mat(k, 120) - dti(k)
         mat(k, 129) = mat(k, 129) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 154) = mat(k, 154) - dti(k)
         mat(k, 167) = mat(k, 167) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 184) = mat(k, 184) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 206) = mat(k, 206) - dti(k)
         mat(k, 226) = mat(k, 226) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 254) = mat(k, 254) - dti(k)
         mat(k, 272) = mat(k, 272) - dti(k)
         mat(k, 288) = mat(k, 288) - dti(k)
         mat(k, 306) = mat(k, 306) - dti(k)
         mat(k, 332) = mat(k, 332) - dti(k)
         mat(k, 359) = mat(k, 359) - dti(k)
         mat(k, 385) = mat(k, 385) - dti(k)
         mat(k, 411) = mat(k, 411) - dti(k)
         mat(k, 437) = mat(k, 437) - dti(k)
         mat(k, 458) = mat(k, 458) - dti(k)
         mat(k, 478) = mat(k, 478) - dti(k)
         mat(k, 502) = mat(k, 502) - dti(k)
         mat(k, 534) = mat(k, 534) - dti(k)
         mat(k, 574) = mat(k, 574) - dti(k)
         mat(k, 596) = mat(k, 596) - dti(k)
         mat(k, 617) = mat(k, 617) - dti(k)
         mat(k, 636) = mat(k, 636) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
