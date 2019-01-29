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
         mat(k,353) = rxt(k,192)*y(k,24)
         mat(k,403) = rxt(k,192)*y(k,2)
         mat(k,258) = (rxt(k,254)+rxt(k,259))*y(k,45)
         mat(k,92) = (rxt(k,254)+rxt(k,259))*y(k,40)
         mat(k,363) = -(4._r8*rxt(k,189)*y(k,2) + (rxt(k,190) + rxt(k,191) + rxt(k,192) &
                      ) * y(k,24) + rxt(k,193)*y(k,43) + rxt(k,194)*y(k,50) + rxt(k,195) &
                      *y(k,51) + rxt(k,197)*y(k,53) + rxt(k,198)*y(k,69))
         mat(k,414) = -(rxt(k,190) + rxt(k,191) + rxt(k,192)) * y(k,2)
         mat(k,559) = -rxt(k,193)*y(k,2)
         mat(k,340) = -rxt(k,194)*y(k,2)
         mat(k,503) = -rxt(k,195)*y(k,2)
         mat(k,598) = -rxt(k,197)*y(k,2)
         mat(k,534) = -rxt(k,198)*y(k,2)
         mat(k,63) = rxt(k,196)*y(k,53)
         mat(k,121) = rxt(k,206)*y(k,65)
         mat(k,96) = rxt(k,201)*y(k,53)
         mat(k,598) = mat(k,598) + rxt(k,196)*y(k,3) + rxt(k,201)*y(k,45)
         mat(k,636) = rxt(k,188)*y(k,58)
         mat(k,459) = rxt(k,188)*y(k,55)
         mat(k,439) = rxt(k,206)*y(k,36)
         mat(k,59) = -(rxt(k,196)*y(k,53))
         mat(k,575) = -rxt(k,196)*y(k,3)
         mat(k,354) = rxt(k,195)*y(k,51)
         mat(k,492) = rxt(k,195)*y(k,2)
         mat(k,282) = -(rxt(k,150)*y(k,59) + rxt(k,186)*y(k,58) + rxt(k,230)*y(k,52) &
                      + rxt(k,231)*y(k,53) + rxt(k,232)*y(k,69))
         mat(k,384) = -rxt(k,150)*y(k,14)
         mat(k,455) = -rxt(k,186)*y(k,14)
         mat(k,314) = -rxt(k,230)*y(k,14)
         mat(k,594) = -rxt(k,231)*y(k,14)
         mat(k,530) = -rxt(k,232)*y(k,14)
         mat(k,474) = rxt(k,157)*y(k,24) + rxt(k,234)*y(k,50)
         mat(k,35) = .300_r8*rxt(k,235)*y(k,69)
         mat(k,410) = rxt(k,157)*y(k,18)
         mat(k,336) = rxt(k,234)*y(k,18)
         mat(k,530) = mat(k,530) + .300_r8*rxt(k,235)*y(k,19)
         mat(k,483) = -(rxt(k,157)*y(k,24) + rxt(k,233)*y(k,43) + rxt(k,234)*y(k,50))
         mat(k,419) = -rxt(k,157)*y(k,18)
         mat(k,564) = -rxt(k,233)*y(k,18)
         mat(k,345) = -rxt(k,234)*y(k,18)
         mat(k,37) = .700_r8*rxt(k,235)*y(k,69)
         mat(k,539) = .700_r8*rxt(k,235)*y(k,19)
         mat(k,33) = -(rxt(k,235)*y(k,69))
         mat(k,517) = -rxt(k,235)*y(k,19)
         mat(k,471) = rxt(k,233)*y(k,43)
         mat(k,546) = rxt(k,233)*y(k,18)
         mat(k,402) = 2.000_r8*rxt(k,159)*y(k,24)
         mat(k,202) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,40) + rxt(k,163)*y(k,59)
         mat(k,257) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,25) + (rxt(k,247) &
                       +rxt(k,253)+rxt(k,258))*y(k,46)
         mat(k,110) = (rxt(k,247)+rxt(k,253)+rxt(k,258))*y(k,40)
         mat(k,375) = rxt(k,163)*y(k,25)
         mat(k,400) = 2.000_r8*rxt(k,184)*y(k,24)
         mat(k,416) = -(rxt(k,157)*y(k,18) + (4._r8*rxt(k,158) + 4._r8*rxt(k,159) &
                      + 4._r8*rxt(k,160) + 4._r8*rxt(k,184)) * y(k,24) + rxt(k,161) &
                      *y(k,43) + rxt(k,162)*y(k,50) + rxt(k,164)*y(k,51) + rxt(k,167) &
                      *y(k,53) + (rxt(k,168) + rxt(k,169)) * y(k,69) + (rxt(k,190) &
                      + rxt(k,191) + rxt(k,192)) * y(k,2))
         mat(k,480) = -rxt(k,157)*y(k,24)
         mat(k,561) = -rxt(k,161)*y(k,24)
         mat(k,342) = -rxt(k,162)*y(k,24)
         mat(k,505) = -rxt(k,164)*y(k,24)
         mat(k,600) = -rxt(k,167)*y(k,24)
         mat(k,536) = -(rxt(k,168) + rxt(k,169)) * y(k,24)
         mat(k,365) = -(rxt(k,190) + rxt(k,191) + rxt(k,192)) * y(k,24)
         mat(k,210) = rxt(k,165)*y(k,53)
         mat(k,269) = rxt(k,183)*y(k,65)
         mat(k,561) = mat(k,561) + rxt(k,155)*y(k,59)
         mat(k,115) = rxt(k,173)*y(k,53) + rxt(k,172)*y(k,59) + rxt(k,174)*y(k,69)
         mat(k,600) = mat(k,600) + rxt(k,165)*y(k,25) + rxt(k,173)*y(k,46)
         mat(k,638) = rxt(k,156)*y(k,59)
         mat(k,390) = rxt(k,155)*y(k,43) + rxt(k,172)*y(k,46) + rxt(k,156)*y(k,55)
         mat(k,441) = rxt(k,183)*y(k,40)
         mat(k,536) = mat(k,536) + rxt(k,174)*y(k,46)
         mat(k,204) = -(rxt(k,163)*y(k,59) + rxt(k,165)*y(k,53) + rxt(k,166)*y(k,69) &
                      + (rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,40))
         mat(k,379) = -rxt(k,163)*y(k,25)
         mat(k,589) = -rxt(k,165)*y(k,25)
         mat(k,525) = -rxt(k,166)*y(k,25)
         mat(k,261) = -(rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,25)
         mat(k,405) = rxt(k,164)*y(k,51)
         mat(k,495) = rxt(k,164)*y(k,24)
         mat(k,67) = -((rxt(k,237) + rxt(k,241)) * y(k,69))
         mat(k,519) = -(rxt(k,237) + rxt(k,241)) * y(k,27)
         mat(k,277) = rxt(k,230)*y(k,52) + rxt(k,231)*y(k,53) + rxt(k,186)*y(k,58) &
                      + rxt(k,150)*y(k,59) + rxt(k,232)*y(k,69)
         mat(k,310) = rxt(k,230)*y(k,14)
         mat(k,576) = rxt(k,231)*y(k,14)
         mat(k,451) = rxt(k,186)*y(k,14)
         mat(k,376) = rxt(k,150)*y(k,14)
         mat(k,519) = mat(k,519) + rxt(k,232)*y(k,14)
         mat(k,4) = -(rxt(k,211)*y(k,65))
         mat(k,426) = -rxt(k,211)*y(k,29)
         mat(k,17) = -(rxt(k,212)*y(k,65))
         mat(k,428) = -rxt(k,212)*y(k,30)
         mat(k,102) = -(rxt(k,208)*y(k,33) + rxt(k,209)*y(k,73) + rxt(k,210)*y(k,42))
         mat(k,611) = -rxt(k,208)*y(k,31)
         mat(k,246) = -rxt(k,209)*y(k,31)
         mat(k,215) = -rxt(k,210)*y(k,31)
         mat(k,5) = 2.000_r8*rxt(k,211)*y(k,65)
         mat(k,18) = rxt(k,212)*y(k,65)
         mat(k,429) = 2.000_r8*rxt(k,211)*y(k,29) + rxt(k,212)*y(k,30)
         mat(k,298) = -((rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,43) + rxt(k,111) &
                      *y(k,54) + rxt(k,114)*y(k,55))
         mat(k,556) = -(rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,32)
         mat(k,237) = -rxt(k,111)*y(k,32)
         mat(k,633) = -rxt(k,114)*y(k,32)
         mat(k,283) = rxt(k,232)*y(k,69)
         mat(k,68) = rxt(k,241)*y(k,69)
         mat(k,105) = rxt(k,208)*y(k,33)
         mat(k,615) = rxt(k,208)*y(k,31) + rxt(k,106)*y(k,53) + rxt(k,152)*y(k,59) &
                      + rxt(k,89)*y(k,65) + rxt(k,115)*y(k,69)
         mat(k,120) = rxt(k,206)*y(k,65)
         mat(k,265) = rxt(k,183)*y(k,65)
         mat(k,195) = rxt(k,138)*y(k,69)
         mat(k,595) = rxt(k,106)*y(k,33) + rxt(k,118)*y(k,69)
         mat(k,385) = rxt(k,152)*y(k,33)
         mat(k,436) = rxt(k,89)*y(k,33) + rxt(k,206)*y(k,36) + rxt(k,183)*y(k,40)
         mat(k,531) = rxt(k,232)*y(k,14) + rxt(k,241)*y(k,27) + rxt(k,115)*y(k,33) &
                      + rxt(k,138)*y(k,47) + rxt(k,118)*y(k,53)
         mat(k,628) = -(rxt(k,89)*y(k,65) + rxt(k,106)*y(k,53) + rxt(k,115)*y(k,69) &
                      + rxt(k,152)*y(k,59) + rxt(k,208)*y(k,31))
         mat(k,449) = -rxt(k,89)*y(k,33)
         mat(k,608) = -rxt(k,106)*y(k,33)
         mat(k,544) = -rxt(k,115)*y(k,33)
         mat(k,398) = -rxt(k,152)*y(k,33)
         mat(k,109) = -rxt(k,208)*y(k,33)
         mat(k,307) = rxt(k,108)*y(k,43)
         mat(k,569) = rxt(k,108)*y(k,32)
         mat(k,71) = -(rxt(k,107)*y(k,53) + rxt(k,116)*y(k,69) + rxt(k,153)*y(k,59))
         mat(k,577) = -rxt(k,107)*y(k,35)
         mat(k,520) = -rxt(k,116)*y(k,35)
         mat(k,377) = -rxt(k,153)*y(k,35)
         mat(k,548) = 2.000_r8*rxt(k,122)*y(k,43)
         mat(k,520) = mat(k,520) + 2.000_r8*rxt(k,121)*y(k,69)
         mat(k,118) = -(rxt(k,199)*y(k,53) + rxt(k,200)*y(k,69) + (rxt(k,205) &
                      + rxt(k,206)) * y(k,65))
         mat(k,582) = -rxt(k,199)*y(k,36)
         mat(k,523) = -rxt(k,200)*y(k,36)
         mat(k,430) = -(rxt(k,205) + rxt(k,206)) * y(k,36)
         mat(k,278) = rxt(k,186)*y(k,58)
         mat(k,551) = rxt(k,187)*y(k,58)
         mat(k,452) = rxt(k,186)*y(k,14) + rxt(k,187)*y(k,43)
         mat(k,264) = -(rxt(k,170)*y(k,53) + rxt(k,171)*y(k,69) + (rxt(k,182) &
                      + rxt(k,183)) * y(k,65) + (rxt(k,247) + rxt(k,253) + rxt(k,258) &
                      ) * y(k,46) + (rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,25) &
                      + (rxt(k,254) + rxt(k,259)) * y(k,45))
         mat(k,593) = -rxt(k,170)*y(k,40)
         mat(k,529) = -rxt(k,171)*y(k,40)
         mat(k,434) = -(rxt(k,182) + rxt(k,183)) * y(k,40)
         mat(k,113) = -(rxt(k,247) + rxt(k,253) + rxt(k,258)) * y(k,40)
         mat(k,207) = -(rxt(k,252) + rxt(k,257) + rxt(k,262)) * y(k,40)
         mat(k,95) = -(rxt(k,254) + rxt(k,259)) * y(k,40)
         mat(k,281) = rxt(k,150)*y(k,59)
         mat(k,409) = rxt(k,169)*y(k,69)
         mat(k,614) = rxt(k,152)*y(k,59)
         mat(k,73) = rxt(k,153)*y(k,59)
         mat(k,554) = rxt(k,154)*y(k,59)
         mat(k,113) = mat(k,113) + rxt(k,172)*y(k,59)
         mat(k,383) = rxt(k,150)*y(k,14) + rxt(k,152)*y(k,33) + rxt(k,153)*y(k,35) &
                      + rxt(k,154)*y(k,43) + rxt(k,172)*y(k,46)
         mat(k,529) = mat(k,529) + rxt(k,169)*y(k,24)
         mat(k,101) = rxt(k,208)*y(k,33) + rxt(k,210)*y(k,42) + rxt(k,209)*y(k,73)
         mat(k,610) = rxt(k,208)*y(k,31)
         mat(k,214) = rxt(k,210)*y(k,31)
         mat(k,245) = rxt(k,209)*y(k,31)
         mat(k,216) = -(rxt(k,147)*y(k,69) + rxt(k,210)*y(k,31))
         mat(k,526) = -rxt(k,147)*y(k,42)
         mat(k,103) = -rxt(k,210)*y(k,42)
         mat(k,279) = rxt(k,230)*y(k,52)
         mat(k,205) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,40)
         mat(k,262) = (rxt(k,252)+rxt(k,257)+rxt(k,262))*y(k,25)
         mat(k,496) = rxt(k,146)*y(k,69)
         mat(k,311) = rxt(k,230)*y(k,14)
         mat(k,526) = mat(k,526) + rxt(k,146)*y(k,51)
         mat(k,567) = -((rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,32) + rxt(k,112) &
                      *y(k,53) + rxt(k,113)*y(k,55) + rxt(k,117)*y(k,69) &
                      + 4._r8*rxt(k,122)*y(k,43) + rxt(k,134)*y(k,52) + rxt(k,139) &
                      *y(k,50) + rxt(k,144)*y(k,51) + (rxt(k,154) + rxt(k,155) &
                      ) * y(k,59) + rxt(k,161)*y(k,24) + rxt(k,187)*y(k,58) + rxt(k,193) &
                      *y(k,2) + rxt(k,233)*y(k,18))
         mat(k,305) = -(rxt(k,108) + rxt(k,109) + rxt(k,110)) * y(k,43)
         mat(k,606) = -rxt(k,112)*y(k,43)
         mat(k,644) = -rxt(k,113)*y(k,43)
         mat(k,542) = -rxt(k,117)*y(k,43)
         mat(k,326) = -rxt(k,134)*y(k,43)
         mat(k,348) = -rxt(k,139)*y(k,43)
         mat(k,511) = -rxt(k,144)*y(k,43)
         mat(k,396) = -(rxt(k,154) + rxt(k,155)) * y(k,43)
         mat(k,422) = -rxt(k,161)*y(k,43)
         mat(k,467) = -rxt(k,187)*y(k,43)
         mat(k,371) = -rxt(k,193)*y(k,43)
         mat(k,486) = -rxt(k,233)*y(k,43)
         mat(k,371) = mat(k,371) + rxt(k,198)*y(k,69)
         mat(k,293) = rxt(k,230)*y(k,52) + rxt(k,231)*y(k,53) + rxt(k,186)*y(k,58) &
                      + rxt(k,150)*y(k,59)
         mat(k,486) = mat(k,486) + rxt(k,157)*y(k,24) + rxt(k,234)*y(k,50)
         mat(k,422) = mat(k,422) + rxt(k,157)*y(k,18) + rxt(k,168)*y(k,69)
         mat(k,70) = rxt(k,237)*y(k,69)
         mat(k,305) = mat(k,305) + rxt(k,111)*y(k,54)
         mat(k,76) = rxt(k,107)*y(k,53) + rxt(k,153)*y(k,59) + rxt(k,116)*y(k,69)
         mat(k,348) = mat(k,348) + rxt(k,234)*y(k,18)
         mat(k,326) = mat(k,326) + rxt(k,230)*y(k,14) + rxt(k,137)*y(k,69)
         mat(k,606) = mat(k,606) + rxt(k,231)*y(k,14) + rxt(k,107)*y(k,35)
         mat(k,242) = rxt(k,111)*y(k,32)
         mat(k,644) = mat(k,644) + rxt(k,119)*y(k,69)
         mat(k,467) = mat(k,467) + rxt(k,186)*y(k,14)
         mat(k,396) = mat(k,396) + rxt(k,150)*y(k,14) + rxt(k,153)*y(k,35)
         mat(k,542) = mat(k,542) + rxt(k,198)*y(k,2) + rxt(k,168)*y(k,24) + rxt(k,237) &
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
         mat(k,518) = -rxt(k,123)*y(k,44)
         mat(k,547) = rxt(k,144)*y(k,51)
         mat(k,491) = rxt(k,144)*y(k,43)
         mat(k,93) = -(rxt(k,201)*y(k,53) + (rxt(k,254) + rxt(k,259)) * y(k,40))
         mat(k,580) = -rxt(k,201)*y(k,45)
         mat(k,259) = -(rxt(k,254) + rxt(k,259)) * y(k,45)
         mat(k,355) = rxt(k,193)*y(k,43)
         mat(k,549) = rxt(k,193)*y(k,2)
         mat(k,111) = -(rxt(k,172)*y(k,59) + rxt(k,173)*y(k,53) + rxt(k,174)*y(k,69) &
                      + (rxt(k,247) + rxt(k,253) + rxt(k,258)) * y(k,40))
         mat(k,378) = -rxt(k,172)*y(k,46)
         mat(k,581) = -rxt(k,173)*y(k,46)
         mat(k,522) = -rxt(k,174)*y(k,46)
         mat(k,260) = -(rxt(k,247) + rxt(k,253) + rxt(k,258)) * y(k,46)
         mat(k,404) = rxt(k,161)*y(k,43)
         mat(k,203) = rxt(k,166)*y(k,69)
         mat(k,550) = rxt(k,161)*y(k,24)
         mat(k,522) = mat(k,522) + rxt(k,166)*y(k,25)
         mat(k,193) = -(rxt(k,126)*y(k,50) + (rxt(k,127) + rxt(k,128) + rxt(k,129) &
                      ) * y(k,51) + rxt(k,130)*y(k,54) + rxt(k,138)*y(k,69) + rxt(k,275) &
                      *y(k,68))
         mat(k,334) = -rxt(k,126)*y(k,47)
         mat(k,494) = -(rxt(k,127) + rxt(k,128) + rxt(k,129)) * y(k,47)
         mat(k,235) = -rxt(k,130)*y(k,47)
         mat(k,524) = -rxt(k,138)*y(k,47)
         mat(k,152) = -rxt(k,275)*y(k,47)
         mat(k,588) = rxt(k,124)*y(k,61) + rxt(k,272)*y(k,64)
         mat(k,235) = mat(k,235) + rxt(k,273)*y(k,64)
         mat(k,166) = 1.100_r8*rxt(k,268)*y(k,62) + .200_r8*rxt(k,266)*y(k,63)
         mat(k,179) = rxt(k,124)*y(k,53)
         mat(k,84) = 1.100_r8*rxt(k,268)*y(k,60)
         mat(k,90) = .200_r8*rxt(k,266)*y(k,60)
         mat(k,132) = rxt(k,272)*y(k,53) + rxt(k,273)*y(k,54)
         mat(k,490) = rxt(k,145)*y(k,52)
         mat(k,309) = rxt(k,145)*y(k,51)
         mat(k,339) = -(rxt(k,126)*y(k,47) + rxt(k,135)*y(k,52) + rxt(k,139)*y(k,43) &
                      + rxt(k,140)*y(k,55) + rxt(k,141)*y(k,53) + rxt(k,162)*y(k,24) &
                      + rxt(k,194)*y(k,2) + rxt(k,234)*y(k,18) + rxt(k,277)*y(k,68))
         mat(k,196) = -rxt(k,126)*y(k,50)
         mat(k,317) = -rxt(k,135)*y(k,50)
         mat(k,558) = -rxt(k,139)*y(k,50)
         mat(k,635) = -rxt(k,140)*y(k,50)
         mat(k,597) = -rxt(k,141)*y(k,50)
         mat(k,413) = -rxt(k,162)*y(k,50)
         mat(k,362) = -rxt(k,194)*y(k,50)
         mat(k,477) = -rxt(k,234)*y(k,50)
         mat(k,154) = -rxt(k,277)*y(k,50)
         mat(k,196) = mat(k,196) + 2.000_r8*rxt(k,128)*y(k,51) + rxt(k,130)*y(k,54) &
                      + rxt(k,138)*y(k,69)
         mat(k,502) = 2.000_r8*rxt(k,128)*y(k,47) + rxt(k,131)*y(k,53)
         mat(k,597) = mat(k,597) + rxt(k,131)*y(k,51)
         mat(k,238) = rxt(k,130)*y(k,47) + rxt(k,125)*y(k,61)
         mat(k,182) = rxt(k,125)*y(k,54)
         mat(k,533) = rxt(k,138)*y(k,47)
         mat(k,509) = -((rxt(k,127) + rxt(k,128) + rxt(k,129)) * y(k,47) + (rxt(k,131) &
                      + rxt(k,133)) * y(k,53) + rxt(k,132)*y(k,55) + rxt(k,144) &
                      *y(k,43) + rxt(k,145)*y(k,52) + rxt(k,146)*y(k,69) + rxt(k,164) &
                      *y(k,24) + rxt(k,195)*y(k,2))
         mat(k,198) = -(rxt(k,127) + rxt(k,128) + rxt(k,129)) * y(k,51)
         mat(k,604) = -(rxt(k,131) + rxt(k,133)) * y(k,51)
         mat(k,642) = -rxt(k,132)*y(k,51)
         mat(k,565) = -rxt(k,144)*y(k,51)
         mat(k,324) = -rxt(k,145)*y(k,51)
         mat(k,540) = -rxt(k,146)*y(k,51)
         mat(k,420) = -rxt(k,164)*y(k,51)
         mat(k,369) = -rxt(k,195)*y(k,51)
         mat(k,369) = mat(k,369) + rxt(k,194)*y(k,50)
         mat(k,484) = rxt(k,234)*y(k,50)
         mat(k,420) = mat(k,420) + rxt(k,162)*y(k,50)
         mat(k,565) = mat(k,565) + rxt(k,139)*y(k,50) + rxt(k,134)*y(k,52)
         mat(k,43) = rxt(k,123)*y(k,69)
         mat(k,346) = rxt(k,194)*y(k,2) + rxt(k,234)*y(k,18) + rxt(k,162)*y(k,24) &
                      + rxt(k,139)*y(k,43) + 2.000_r8*rxt(k,135)*y(k,52) + rxt(k,141) &
                      *y(k,53) + rxt(k,140)*y(k,55)
         mat(k,324) = mat(k,324) + rxt(k,134)*y(k,43) + 2.000_r8*rxt(k,135)*y(k,50) &
                      + rxt(k,136)*y(k,53) + rxt(k,137)*y(k,69)
         mat(k,604) = mat(k,604) + rxt(k,141)*y(k,50) + rxt(k,136)*y(k,52)
         mat(k,642) = mat(k,642) + rxt(k,140)*y(k,50)
         mat(k,540) = mat(k,540) + rxt(k,123)*y(k,44) + rxt(k,137)*y(k,52)
         mat(k,316) = -(rxt(k,134)*y(k,43) + rxt(k,135)*y(k,50) + rxt(k,136)*y(k,53) &
                      + rxt(k,137)*y(k,69) + rxt(k,145)*y(k,51) + rxt(k,230)*y(k,14))
         mat(k,557) = -rxt(k,134)*y(k,52)
         mat(k,338) = -rxt(k,135)*y(k,52)
         mat(k,596) = -rxt(k,136)*y(k,52)
         mat(k,532) = -rxt(k,137)*y(k,52)
         mat(k,501) = -rxt(k,145)*y(k,52)
         mat(k,284) = -rxt(k,230)*y(k,52)
         mat(k,62) = rxt(k,196)*y(k,53)
         mat(k,208) = rxt(k,165)*y(k,53) + rxt(k,163)*y(k,59) + rxt(k,166)*y(k,69)
         mat(k,106) = rxt(k,210)*y(k,42)
         mat(k,219) = rxt(k,210)*y(k,31) + rxt(k,147)*y(k,69)
         mat(k,501) = mat(k,501) + rxt(k,133)*y(k,53) + rxt(k,132)*y(k,55)
         mat(k,596) = mat(k,596) + rxt(k,196)*y(k,3) + rxt(k,165)*y(k,25) + rxt(k,133) &
                      *y(k,51)
         mat(k,634) = rxt(k,132)*y(k,51)
         mat(k,386) = rxt(k,163)*y(k,25)
         mat(k,532) = mat(k,532) + rxt(k,166)*y(k,25) + rxt(k,147)*y(k,42)
         mat(k,607) = -(rxt(k,103)*y(k,55) + 4._r8*rxt(k,104)*y(k,53) + rxt(k,105) &
                      *y(k,54) + rxt(k,106)*y(k,33) + rxt(k,107)*y(k,35) + rxt(k,112) &
                      *y(k,43) + rxt(k,118)*y(k,69) + (rxt(k,131) + rxt(k,133) &
                      ) * y(k,51) + rxt(k,136)*y(k,52) + rxt(k,141)*y(k,50) + rxt(k,165) &
                      *y(k,25) + rxt(k,167)*y(k,24) + rxt(k,170)*y(k,40) + rxt(k,173) &
                      *y(k,46) + rxt(k,196)*y(k,3) + rxt(k,197)*y(k,2) + rxt(k,199) &
                      *y(k,36) + rxt(k,201)*y(k,45) + rxt(k,231)*y(k,14) + (rxt(k,270) &
                      + rxt(k,271)) * y(k,62) + rxt(k,272)*y(k,64))
         mat(k,645) = -rxt(k,103)*y(k,53)
         mat(k,243) = -rxt(k,105)*y(k,53)
         mat(k,627) = -rxt(k,106)*y(k,53)
         mat(k,77) = -rxt(k,107)*y(k,53)
         mat(k,568) = -rxt(k,112)*y(k,53)
         mat(k,543) = -rxt(k,118)*y(k,53)
         mat(k,512) = -(rxt(k,131) + rxt(k,133)) * y(k,53)
         mat(k,327) = -rxt(k,136)*y(k,53)
         mat(k,349) = -rxt(k,141)*y(k,53)
         mat(k,213) = -rxt(k,165)*y(k,53)
         mat(k,423) = -rxt(k,167)*y(k,53)
         mat(k,275) = -rxt(k,170)*y(k,53)
         mat(k,117) = -rxt(k,173)*y(k,53)
         mat(k,66) = -rxt(k,196)*y(k,53)
         mat(k,372) = -rxt(k,197)*y(k,53)
         mat(k,125) = -rxt(k,199)*y(k,53)
         mat(k,100) = -rxt(k,201)*y(k,53)
         mat(k,294) = -rxt(k,231)*y(k,53)
         mat(k,86) = -(rxt(k,270) + rxt(k,271)) * y(k,53)
         mat(k,134) = -rxt(k,272)*y(k,53)
         mat(k,306) = rxt(k,110)*y(k,43)
         mat(k,568) = mat(k,568) + rxt(k,110)*y(k,32)
         mat(k,201) = rxt(k,126)*y(k,50) + rxt(k,127)*y(k,51) + rxt(k,130)*y(k,54) &
                      + rxt(k,275)*y(k,68)
         mat(k,349) = mat(k,349) + rxt(k,126)*y(k,47)
         mat(k,512) = mat(k,512) + rxt(k,127)*y(k,47)
         mat(k,243) = mat(k,243) + rxt(k,130)*y(k,47) + rxt(k,274)*y(k,64) + ( &
                      + rxt(k,92)+rxt(k,93))*y(k,65) + rxt(k,281)*y(k,70) + rxt(k,285) &
                      *y(k,71)
         mat(k,173) = rxt(k,266)*y(k,63) + 1.150_r8*rxt(k,267)*y(k,68)
         mat(k,186) = rxt(k,280)*y(k,70)
         mat(k,91) = rxt(k,266)*y(k,60)
         mat(k,134) = mat(k,134) + rxt(k,274)*y(k,54)
         mat(k,448) = (rxt(k,92)+rxt(k,93))*y(k,54)
         mat(k,156) = rxt(k,275)*y(k,47) + 1.150_r8*rxt(k,267)*y(k,60)
         mat(k,543) = mat(k,543) + 2.000_r8*rxt(k,120)*y(k,69)
         mat(k,147) = rxt(k,281)*y(k,54) + rxt(k,280)*y(k,61)
         mat(k,58) = rxt(k,285)*y(k,54)
         mat(k,236) = -(rxt(k,92)*y(k,65) + rxt(k,97)*y(k,66) + rxt(k,105)*y(k,53) &
                      + rxt(k,111)*y(k,32) + rxt(k,125)*y(k,61) + rxt(k,130)*y(k,47) &
                      + rxt(k,269)*y(k,62) + (rxt(k,273) + rxt(k,274)) * y(k,64) &
                      + rxt(k,281)*y(k,70) + rxt(k,285)*y(k,71))
         mat(k,432) = -rxt(k,92)*y(k,54)
         mat(k,11) = -rxt(k,97)*y(k,54)
         mat(k,591) = -rxt(k,105)*y(k,54)
         mat(k,296) = -rxt(k,111)*y(k,54)
         mat(k,180) = -rxt(k,125)*y(k,54)
         mat(k,194) = -rxt(k,130)*y(k,54)
         mat(k,85) = -rxt(k,269)*y(k,54)
         mat(k,133) = -(rxt(k,273) + rxt(k,274)) * y(k,54)
         mat(k,143) = -rxt(k,281)*y(k,54)
         mat(k,57) = -rxt(k,285)*y(k,54)
         mat(k,357) = 2.000_r8*rxt(k,189)*y(k,2) + (rxt(k,191)+rxt(k,192))*y(k,24) &
                      + rxt(k,193)*y(k,43) + rxt(k,197)*y(k,53)
         mat(k,472) = rxt(k,233)*y(k,43)
         mat(k,407) = (rxt(k,191)+rxt(k,192))*y(k,2) + (2.000_r8*rxt(k,158) &
                       +2.000_r8*rxt(k,159))*y(k,24) + rxt(k,161)*y(k,43) + rxt(k,167) &
                      *y(k,53) + rxt(k,169)*y(k,69)
         mat(k,296) = mat(k,296) + rxt(k,108)*y(k,43) + rxt(k,114)*y(k,55)
         mat(k,552) = rxt(k,193)*y(k,2) + rxt(k,233)*y(k,18) + rxt(k,161)*y(k,24) &
                      + rxt(k,108)*y(k,32) + 2.000_r8*rxt(k,122)*y(k,43) + rxt(k,134) &
                      *y(k,52) + rxt(k,112)*y(k,53) + 2.000_r8*rxt(k,113)*y(k,55) &
                      + rxt(k,187)*y(k,58) + rxt(k,154)*y(k,59) + rxt(k,117)*y(k,69)
         mat(k,40) = rxt(k,123)*y(k,69)
         mat(k,194) = mat(k,194) + rxt(k,129)*y(k,51)
         mat(k,335) = rxt(k,140)*y(k,55) + rxt(k,277)*y(k,68)
         mat(k,497) = rxt(k,129)*y(k,47) + rxt(k,131)*y(k,53) + rxt(k,132)*y(k,55)
         mat(k,312) = rxt(k,134)*y(k,43) + rxt(k,136)*y(k,53)
         mat(k,591) = mat(k,591) + rxt(k,197)*y(k,2) + rxt(k,167)*y(k,24) + rxt(k,112) &
                      *y(k,43) + rxt(k,131)*y(k,51) + rxt(k,136)*y(k,52) &
                      + 2.000_r8*rxt(k,104)*y(k,53) + 2.000_r8*rxt(k,103)*y(k,55) &
                      + rxt(k,96)*y(k,66) + rxt(k,118)*y(k,69)
         mat(k,236) = mat(k,236) + 2.000_r8*rxt(k,97)*y(k,66)
         mat(k,632) = rxt(k,114)*y(k,32) + 2.000_r8*rxt(k,113)*y(k,43) + rxt(k,140) &
                      *y(k,50) + rxt(k,132)*y(k,51) + 2.000_r8*rxt(k,103)*y(k,53) &
                      + rxt(k,188)*y(k,58) + rxt(k,156)*y(k,59) + 2.000_r8*rxt(k,94) &
                      *y(k,65) + rxt(k,119)*y(k,69)
         mat(k,453) = rxt(k,187)*y(k,43) + rxt(k,188)*y(k,55)
         mat(k,381) = rxt(k,154)*y(k,43) + rxt(k,156)*y(k,55)
         mat(k,432) = mat(k,432) + 2.000_r8*rxt(k,94)*y(k,55)
         mat(k,11) = mat(k,11) + rxt(k,96)*y(k,53) + 2.000_r8*rxt(k,97)*y(k,54)
         mat(k,153) = rxt(k,277)*y(k,50)
         mat(k,527) = rxt(k,169)*y(k,24) + rxt(k,117)*y(k,43) + rxt(k,123)*y(k,44) &
                      + rxt(k,118)*y(k,53) + rxt(k,119)*y(k,55)
         mat(k,647) = -(rxt(k,94)*y(k,65) + rxt(k,103)*y(k,53) + rxt(k,113)*y(k,43) &
                      + rxt(k,114)*y(k,32) + rxt(k,119)*y(k,69) + rxt(k,132)*y(k,51) &
                      + rxt(k,140)*y(k,50) + rxt(k,156)*y(k,59) + rxt(k,188)*y(k,58))
         mat(k,450) = -rxt(k,94)*y(k,55)
         mat(k,609) = -rxt(k,103)*y(k,55)
         mat(k,570) = -rxt(k,113)*y(k,55)
         mat(k,308) = -rxt(k,114)*y(k,55)
         mat(k,545) = -rxt(k,119)*y(k,55)
         mat(k,514) = -rxt(k,132)*y(k,55)
         mat(k,351) = -rxt(k,140)*y(k,55)
         mat(k,399) = -rxt(k,156)*y(k,55)
         mat(k,470) = -rxt(k,188)*y(k,55)
         mat(k,609) = mat(k,609) + rxt(k,105)*y(k,54)
         mat(k,244) = rxt(k,105)*y(k,53)
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
         mat(k,352) = rxt(k,190)*y(k,24)
         mat(k,401) = rxt(k,190)*y(k,2) + 2.000_r8*rxt(k,160)*y(k,24)
         mat(k,463) = -(rxt(k,186)*y(k,14) + rxt(k,187)*y(k,43) + rxt(k,188)*y(k,55))
         mat(k,289) = -rxt(k,186)*y(k,58)
         mat(k,563) = -rxt(k,187)*y(k,58)
         mat(k,640) = -rxt(k,188)*y(k,58)
         mat(k,367) = 4.000_r8*rxt(k,189)*y(k,2) + (rxt(k,190)+rxt(k,191))*y(k,24) &
                      + rxt(k,194)*y(k,50) + rxt(k,197)*y(k,53) + rxt(k,198)*y(k,69)
         mat(k,418) = (rxt(k,190)+rxt(k,191))*y(k,2)
         mat(k,123) = rxt(k,199)*y(k,53) + rxt(k,205)*y(k,65) + rxt(k,200)*y(k,69)
         mat(k,344) = rxt(k,194)*y(k,2)
         mat(k,602) = rxt(k,197)*y(k,2) + rxt(k,199)*y(k,36)
         mat(k,443) = rxt(k,205)*y(k,36)
         mat(k,538) = rxt(k,198)*y(k,2) + rxt(k,200)*y(k,36)
         mat(k,389) = -(rxt(k,150)*y(k,14) + rxt(k,152)*y(k,33) + rxt(k,153)*y(k,35) &
                      + (rxt(k,154) + rxt(k,155)) * y(k,43) + rxt(k,156)*y(k,55) &
                      + rxt(k,163)*y(k,25) + rxt(k,172)*y(k,46))
         mat(k,286) = -rxt(k,150)*y(k,59)
         mat(k,619) = -rxt(k,152)*y(k,59)
         mat(k,74) = -rxt(k,153)*y(k,59)
         mat(k,560) = -(rxt(k,154) + rxt(k,155)) * y(k,59)
         mat(k,637) = -rxt(k,156)*y(k,59)
         mat(k,209) = -rxt(k,163)*y(k,59)
         mat(k,114) = -rxt(k,172)*y(k,59)
         mat(k,364) = rxt(k,191)*y(k,24)
         mat(k,479) = rxt(k,157)*y(k,24)
         mat(k,415) = rxt(k,191)*y(k,2) + rxt(k,157)*y(k,18) + (4.000_r8*rxt(k,158) &
                       +2.000_r8*rxt(k,160))*y(k,24) + rxt(k,162)*y(k,50) + rxt(k,167) &
                      *y(k,53) + rxt(k,168)*y(k,69)
         mat(k,19) = rxt(k,212)*y(k,65)
         mat(k,268) = rxt(k,170)*y(k,53) + rxt(k,182)*y(k,65) + rxt(k,171)*y(k,69)
         mat(k,341) = rxt(k,162)*y(k,24)
         mat(k,599) = rxt(k,167)*y(k,24) + rxt(k,170)*y(k,40)
         mat(k,440) = rxt(k,212)*y(k,30) + rxt(k,182)*y(k,40)
         mat(k,535) = rxt(k,168)*y(k,24) + rxt(k,171)*y(k,40)
         mat(k,164) = -(rxt(k,266)*y(k,63) + rxt(k,267)*y(k,68) + rxt(k,268)*y(k,62))
         mat(k,88) = -rxt(k,266)*y(k,60)
         mat(k,150) = -rxt(k,267)*y(k,60)
         mat(k,82) = -rxt(k,268)*y(k,60)
         mat(k,178) = -(rxt(k,124)*y(k,53) + rxt(k,125)*y(k,54) + rxt(k,280)*y(k,70))
         mat(k,587) = -rxt(k,124)*y(k,61)
         mat(k,234) = -rxt(k,125)*y(k,61)
         mat(k,141) = -rxt(k,280)*y(k,61)
         mat(k,587) = mat(k,587) + rxt(k,270)*y(k,62)
         mat(k,165) = .900_r8*rxt(k,268)*y(k,62) + .800_r8*rxt(k,266)*y(k,63)
         mat(k,83) = rxt(k,270)*y(k,53) + .900_r8*rxt(k,268)*y(k,60)
         mat(k,89) = .800_r8*rxt(k,266)*y(k,60)
         mat(k,78) = -(rxt(k,268)*y(k,60) + rxt(k,269)*y(k,54) + (rxt(k,270) + rxt(k,271) &
                      ) * y(k,53))
         mat(k,159) = -rxt(k,268)*y(k,62)
         mat(k,228) = -rxt(k,269)*y(k,62)
         mat(k,578) = -(rxt(k,270) + rxt(k,271)) * y(k,62)
         mat(k,87) = -(rxt(k,266)*y(k,60))
         mat(k,160) = -rxt(k,266)*y(k,63)
         mat(k,187) = rxt(k,275)*y(k,68)
         mat(k,330) = rxt(k,277)*y(k,68)
         mat(k,579) = rxt(k,270)*y(k,62)
         mat(k,229) = rxt(k,274)*y(k,64)
         mat(k,79) = rxt(k,270)*y(k,53)
         mat(k,126) = rxt(k,274)*y(k,54)
         mat(k,148) = rxt(k,275)*y(k,47) + rxt(k,277)*y(k,50)
         mat(k,127) = -(rxt(k,272)*y(k,53) + (rxt(k,273) + rxt(k,274)) * y(k,54))
         mat(k,583) = -rxt(k,272)*y(k,64)
         mat(k,230) = -(rxt(k,273) + rxt(k,274)) * y(k,64)
         mat(k,174) = rxt(k,280)*y(k,70)
         mat(k,137) = rxt(k,280)*y(k,61)
         mat(k,442) = -(rxt(k,89)*y(k,33) + rxt(k,90)*y(k,73) + (rxt(k,92) + rxt(k,93) &
                      ) * y(k,54) + rxt(k,94)*y(k,55) + (rxt(k,182) + rxt(k,183) &
                      ) * y(k,40) + (rxt(k,205) + rxt(k,206)) * y(k,36) + rxt(k,211) &
                      *y(k,29) + rxt(k,212)*y(k,30))
         mat(k,621) = -rxt(k,89)*y(k,65)
         mat(k,251) = -rxt(k,90)*y(k,65)
         mat(k,239) = -(rxt(k,92) + rxt(k,93)) * y(k,65)
         mat(k,639) = -rxt(k,94)*y(k,65)
         mat(k,270) = -(rxt(k,182) + rxt(k,183)) * y(k,65)
         mat(k,122) = -(rxt(k,205) + rxt(k,206)) * y(k,65)
         mat(k,6) = -rxt(k,211)*y(k,65)
         mat(k,20) = -rxt(k,212)*y(k,65)
         mat(k,239) = mat(k,239) + rxt(k,125)*y(k,61)
         mat(k,170) = .850_r8*rxt(k,267)*y(k,68)
         mat(k,183) = rxt(k,125)*y(k,54)
         mat(k,155) = .850_r8*rxt(k,267)*y(k,60)
         mat(k,10) = -(rxt(k,96)*y(k,53) + rxt(k,97)*y(k,54))
         mat(k,571) = -rxt(k,96)*y(k,66)
         mat(k,224) = -rxt(k,97)*y(k,66)
         mat(k,571) = mat(k,571) + rxt(k,100)*y(k,67)
         mat(k,224) = mat(k,224) + rxt(k,101)*y(k,67)
         mat(k,630) = rxt(k,102)*y(k,67)
         mat(k,12) = rxt(k,100)*y(k,53) + rxt(k,101)*y(k,54) + rxt(k,102)*y(k,55)
         mat(k,13) = -(rxt(k,100)*y(k,53) + rxt(k,101)*y(k,54) + rxt(k,102)*y(k,55))
         mat(k,572) = -rxt(k,100)*y(k,67)
         mat(k,225) = -rxt(k,101)*y(k,67)
         mat(k,631) = -rxt(k,102)*y(k,67)
         mat(k,225) = mat(k,225) + rxt(k,92)*y(k,65)
         mat(k,427) = rxt(k,92)*y(k,54)
         mat(k,149) = -(rxt(k,267)*y(k,60) + rxt(k,275)*y(k,47) + rxt(k,277)*y(k,50))
         mat(k,163) = -rxt(k,267)*y(k,68)
         mat(k,190) = -rxt(k,275)*y(k,68)
         mat(k,331) = -rxt(k,277)*y(k,68)
         mat(k,232) = rxt(k,269)*y(k,62) + rxt(k,273)*y(k,64) + rxt(k,281)*y(k,70) &
                      + rxt(k,285)*y(k,71)
         mat(k,81) = rxt(k,269)*y(k,54)
         mat(k,129) = rxt(k,273)*y(k,54)
         mat(k,139) = rxt(k,281)*y(k,54)
         mat(k,56) = rxt(k,285)*y(k,54)
         mat(k,541) = -(rxt(k,115)*y(k,33) + rxt(k,116)*y(k,35) + rxt(k,117)*y(k,43) &
                      + rxt(k,118)*y(k,53) + rxt(k,119)*y(k,55) + (4._r8*rxt(k,120) &
                      + 4._r8*rxt(k,121)) * y(k,69) + rxt(k,123)*y(k,44) + rxt(k,137) &
                      *y(k,52) + rxt(k,138)*y(k,47) + rxt(k,146)*y(k,51) + rxt(k,147) &
                      *y(k,42) + rxt(k,166)*y(k,25) + (rxt(k,168) + rxt(k,169) &
                      ) * y(k,24) + rxt(k,171)*y(k,40) + rxt(k,174)*y(k,46) + rxt(k,198) &
                      *y(k,2) + rxt(k,200)*y(k,36) + rxt(k,232)*y(k,14) + rxt(k,235) &
                      *y(k,19) + (rxt(k,237) + rxt(k,241)) * y(k,27))
         mat(k,625) = -rxt(k,115)*y(k,69)
         mat(k,75) = -rxt(k,116)*y(k,69)
         mat(k,566) = -rxt(k,117)*y(k,69)
         mat(k,605) = -rxt(k,118)*y(k,69)
         mat(k,643) = -rxt(k,119)*y(k,69)
         mat(k,44) = -rxt(k,123)*y(k,69)
         mat(k,325) = -rxt(k,137)*y(k,69)
         mat(k,199) = -rxt(k,138)*y(k,69)
         mat(k,510) = -rxt(k,146)*y(k,69)
         mat(k,222) = -rxt(k,147)*y(k,69)
         mat(k,212) = -rxt(k,166)*y(k,69)
         mat(k,421) = -(rxt(k,168) + rxt(k,169)) * y(k,69)
         mat(k,274) = -rxt(k,171)*y(k,69)
         mat(k,116) = -rxt(k,174)*y(k,69)
         mat(k,370) = -rxt(k,198)*y(k,69)
         mat(k,124) = -rxt(k,200)*y(k,69)
         mat(k,292) = -rxt(k,232)*y(k,69)
         mat(k,38) = -rxt(k,235)*y(k,69)
         mat(k,69) = -(rxt(k,237) + rxt(k,241)) * y(k,69)
         mat(k,292) = mat(k,292) + rxt(k,231)*y(k,53)
         mat(k,38) = mat(k,38) + .300_r8*rxt(k,235)*y(k,69)
         mat(k,108) = rxt(k,209)*y(k,73)
         mat(k,304) = 2.000_r8*rxt(k,109)*y(k,43) + rxt(k,114)*y(k,55)
         mat(k,625) = mat(k,625) + rxt(k,106)*y(k,53) + rxt(k,89)*y(k,65)
         mat(k,75) = mat(k,75) + rxt(k,107)*y(k,53)
         mat(k,124) = mat(k,124) + rxt(k,199)*y(k,53) + rxt(k,205)*y(k,65)
         mat(k,274) = mat(k,274) + rxt(k,170)*y(k,53) + rxt(k,182)*y(k,65)
         mat(k,566) = mat(k,566) + 2.000_r8*rxt(k,109)*y(k,32) + rxt(k,139)*y(k,50) &
                      + rxt(k,134)*y(k,52) + rxt(k,112)*y(k,53) + rxt(k,113)*y(k,55) &
                      + rxt(k,155)*y(k,59)
         mat(k,99) = rxt(k,201)*y(k,53)
         mat(k,116) = mat(k,116) + rxt(k,173)*y(k,53)
         mat(k,347) = rxt(k,139)*y(k,43)
         mat(k,325) = mat(k,325) + rxt(k,134)*y(k,43)
         mat(k,605) = mat(k,605) + rxt(k,231)*y(k,14) + rxt(k,106)*y(k,33) &
                      + rxt(k,107)*y(k,35) + rxt(k,199)*y(k,36) + rxt(k,170)*y(k,40) &
                      + rxt(k,112)*y(k,43) + rxt(k,201)*y(k,45) + rxt(k,173)*y(k,46)
         mat(k,643) = mat(k,643) + rxt(k,114)*y(k,32) + rxt(k,113)*y(k,43)
         mat(k,395) = rxt(k,155)*y(k,43)
         mat(k,446) = rxt(k,89)*y(k,33) + rxt(k,205)*y(k,36) + rxt(k,182)*y(k,40) &
                      + 2.000_r8*rxt(k,90)*y(k,73)
         mat(k,541) = mat(k,541) + .300_r8*rxt(k,235)*y(k,19)
         mat(k,254) = rxt(k,209)*y(k,31) + 2.000_r8*rxt(k,90)*y(k,65)
         mat(k,138) = -(rxt(k,280)*y(k,61) + rxt(k,281)*y(k,54))
         mat(k,175) = -rxt(k,280)*y(k,70)
         mat(k,231) = -rxt(k,281)*y(k,70)
         mat(k,584) = rxt(k,271)*y(k,62) + rxt(k,272)*y(k,64) + rxt(k,284)*y(k,71) &
                      + rxt(k,290)*y(k,72)
         mat(k,162) = rxt(k,282)*y(k,71) + rxt(k,287)*y(k,72)
         mat(k,80) = rxt(k,271)*y(k,53)
         mat(k,128) = rxt(k,272)*y(k,53)
         mat(k,55) = rxt(k,284)*y(k,53) + rxt(k,282)*y(k,60)
         mat(k,50) = rxt(k,290)*y(k,53) + rxt(k,287)*y(k,60)
         mat(k,53) = -(rxt(k,282)*y(k,60) + rxt(k,284)*y(k,53) + rxt(k,285)*y(k,54))
         mat(k,158) = -rxt(k,282)*y(k,71)
         mat(k,574) = -rxt(k,284)*y(k,71)
         mat(k,227) = -rxt(k,285)*y(k,71)
         mat(k,158) = mat(k,158) + rxt(k,286)*y(k,72)
         mat(k,47) = rxt(k,286)*y(k,60)
         mat(k,46) = -((rxt(k,286) + rxt(k,287)) * y(k,60) + rxt(k,290)*y(k,53))
         mat(k,157) = -(rxt(k,286) + rxt(k,287)) * y(k,72)
         mat(k,573) = -rxt(k,290)*y(k,72)
         mat(k,248) = -(rxt(k,90)*y(k,65) + rxt(k,209)*y(k,31))
         mat(k,433) = -rxt(k,90)*y(k,73)
         mat(k,104) = -rxt(k,209)*y(k,73)
         mat(k,280) = rxt(k,232)*y(k,69)
         mat(k,34) = rxt(k,235)*y(k,69)
         mat(k,297) = rxt(k,110)*y(k,43)
         mat(k,613) = rxt(k,115)*y(k,69)
         mat(k,72) = rxt(k,116)*y(k,69)
         mat(k,119) = rxt(k,200)*y(k,69)
         mat(k,263) = (rxt(k,254)+rxt(k,259))*y(k,45) + (rxt(k,247)+rxt(k,253) &
                       +rxt(k,258))*y(k,46) + rxt(k,171)*y(k,69)
         mat(k,217) = rxt(k,147)*y(k,69)
         mat(k,553) = rxt(k,110)*y(k,32) + rxt(k,117)*y(k,69)
         mat(k,41) = rxt(k,123)*y(k,69)
         mat(k,94) = (rxt(k,254)+rxt(k,259))*y(k,40)
         mat(k,112) = (rxt(k,247)+rxt(k,253)+rxt(k,258))*y(k,40) + rxt(k,174)*y(k,69)
         mat(k,528) = rxt(k,232)*y(k,14) + rxt(k,235)*y(k,19) + rxt(k,115)*y(k,33) &
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
         mat(k, 60) = lmat(k, 60)
         mat(k, 61) = lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = mat(k, 99) + lmat(k, 99)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 107) = lmat(k, 107)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = lmat(k, 136)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 142) = lmat(k, 142)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 188) = lmat(k, 188)
         mat(k, 191) = lmat(k, 191)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 211) = lmat(k, 211)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 221) = lmat(k, 221)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 249) = lmat(k, 249)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 255) = lmat(k, 255)
         mat(k, 256) = lmat(k, 256)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 295) = lmat(k, 295)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 312) = mat(k, 312) + lmat(k, 312)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 332) = lmat(k, 332)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 372) = mat(k, 372) + lmat(k, 372)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 383) = mat(k, 383) + lmat(k, 383)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 392) = lmat(k, 392)
         mat(k, 393) = lmat(k, 393)
         mat(k, 396) = mat(k, 396) + lmat(k, 396)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 426) = mat(k, 426) + lmat(k, 426)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 435) = lmat(k, 435)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 438) = lmat(k, 438)
         mat(k, 440) = mat(k, 440) + lmat(k, 440)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 444) = lmat(k, 444)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 447) = lmat(k, 447)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 463) = mat(k, 463) + lmat(k, 463)
         mat(k, 483) = mat(k, 483) + lmat(k, 483)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 502) = mat(k, 502) + lmat(k, 502)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 515) = lmat(k, 515)
         mat(k, 516) = lmat(k, 516)
         mat(k, 528) = mat(k, 528) + lmat(k, 528)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 538) = mat(k, 538) + lmat(k, 538)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 586) = lmat(k, 586)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 630) = mat(k, 630) + lmat(k, 630)
         mat(k, 632) = mat(k, 632) + lmat(k, 632)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 645) = mat(k, 645) + lmat(k, 645)
         mat(k, 647) = mat(k, 647) + lmat(k, 647)
         mat(k, 97) = 0._r8
         mat(k, 130) = 0._r8
         mat(k, 131) = 0._r8
         mat(k, 140) = 0._r8
         mat(k, 144) = 0._r8
         mat(k, 145) = 0._r8
         mat(k, 146) = 0._r8
         mat(k, 151) = 0._r8
         mat(k, 161) = 0._r8
         mat(k, 167) = 0._r8
         mat(k, 168) = 0._r8
         mat(k, 169) = 0._r8
         mat(k, 171) = 0._r8
         mat(k, 172) = 0._r8
         mat(k, 176) = 0._r8
         mat(k, 177) = 0._r8
         mat(k, 181) = 0._r8
         mat(k, 184) = 0._r8
         mat(k, 185) = 0._r8
         mat(k, 189) = 0._r8
         mat(k, 192) = 0._r8
         mat(k, 197) = 0._r8
         mat(k, 200) = 0._r8
         mat(k, 206) = 0._r8
         mat(k, 218) = 0._r8
         mat(k, 220) = 0._r8
         mat(k, 223) = 0._r8
         mat(k, 240) = 0._r8
         mat(k, 241) = 0._r8
         mat(k, 247) = 0._r8
         mat(k, 250) = 0._r8
         mat(k, 252) = 0._r8
         mat(k, 253) = 0._r8
         mat(k, 266) = 0._r8
         mat(k, 267) = 0._r8
         mat(k, 271) = 0._r8
         mat(k, 272) = 0._r8
         mat(k, 273) = 0._r8
         mat(k, 276) = 0._r8
         mat(k, 285) = 0._r8
         mat(k, 287) = 0._r8
         mat(k, 288) = 0._r8
         mat(k, 290) = 0._r8
         mat(k, 291) = 0._r8
         mat(k, 299) = 0._r8
         mat(k, 300) = 0._r8
         mat(k, 301) = 0._r8
         mat(k, 302) = 0._r8
         mat(k, 303) = 0._r8
         mat(k, 313) = 0._r8
         mat(k, 315) = 0._r8
         mat(k, 318) = 0._r8
         mat(k, 319) = 0._r8
         mat(k, 320) = 0._r8
         mat(k, 321) = 0._r8
         mat(k, 322) = 0._r8
         mat(k, 323) = 0._r8
         mat(k, 328) = 0._r8
         mat(k, 329) = 0._r8
         mat(k, 333) = 0._r8
         mat(k, 337) = 0._r8
         mat(k, 343) = 0._r8
         mat(k, 350) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 358) = 0._r8
         mat(k, 359) = 0._r8
         mat(k, 360) = 0._r8
         mat(k, 361) = 0._r8
         mat(k, 366) = 0._r8
         mat(k, 368) = 0._r8
         mat(k, 373) = 0._r8
         mat(k, 374) = 0._r8
         mat(k, 380) = 0._r8
         mat(k, 382) = 0._r8
         mat(k, 387) = 0._r8
         mat(k, 388) = 0._r8
         mat(k, 391) = 0._r8
         mat(k, 394) = 0._r8
         mat(k, 397) = 0._r8
         mat(k, 406) = 0._r8
         mat(k, 408) = 0._r8
         mat(k, 411) = 0._r8
         mat(k, 412) = 0._r8
         mat(k, 417) = 0._r8
         mat(k, 424) = 0._r8
         mat(k, 425) = 0._r8
         mat(k, 431) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 445) = 0._r8
         mat(k, 454) = 0._r8
         mat(k, 456) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 458) = 0._r8
         mat(k, 460) = 0._r8
         mat(k, 461) = 0._r8
         mat(k, 462) = 0._r8
         mat(k, 464) = 0._r8
         mat(k, 465) = 0._r8
         mat(k, 466) = 0._r8
         mat(k, 468) = 0._r8
         mat(k, 469) = 0._r8
         mat(k, 473) = 0._r8
         mat(k, 475) = 0._r8
         mat(k, 476) = 0._r8
         mat(k, 478) = 0._r8
         mat(k, 481) = 0._r8
         mat(k, 482) = 0._r8
         mat(k, 485) = 0._r8
         mat(k, 487) = 0._r8
         mat(k, 488) = 0._r8
         mat(k, 489) = 0._r8
         mat(k, 493) = 0._r8
         mat(k, 498) = 0._r8
         mat(k, 499) = 0._r8
         mat(k, 500) = 0._r8
         mat(k, 504) = 0._r8
         mat(k, 506) = 0._r8
         mat(k, 507) = 0._r8
         mat(k, 508) = 0._r8
         mat(k, 513) = 0._r8
         mat(k, 521) = 0._r8
         mat(k, 537) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 590) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 601) = 0._r8
         mat(k, 603) = 0._r8
         mat(k, 612) = 0._r8
         mat(k, 616) = 0._r8
         mat(k, 617) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 620) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 623) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 641) = 0._r8
         mat(k, 646) = 0._r8
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
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 78) = mat(k, 78) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 118) = mat(k, 118) - dti(k)
         mat(k, 127) = mat(k, 127) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 149) = mat(k, 149) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 236) = mat(k, 236) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 264) = mat(k, 264) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 298) = mat(k, 298) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 363) = mat(k, 363) - dti(k)
         mat(k, 389) = mat(k, 389) - dti(k)
         mat(k, 416) = mat(k, 416) - dti(k)
         mat(k, 442) = mat(k, 442) - dti(k)
         mat(k, 463) = mat(k, 463) - dti(k)
         mat(k, 483) = mat(k, 483) - dti(k)
         mat(k, 509) = mat(k, 509) - dti(k)
         mat(k, 541) = mat(k, 541) - dti(k)
         mat(k, 567) = mat(k, 567) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 628) = mat(k, 628) - dti(k)
         mat(k, 647) = mat(k, 647) - dti(k)
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
