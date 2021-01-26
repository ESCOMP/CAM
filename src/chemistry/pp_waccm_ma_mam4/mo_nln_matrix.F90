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
         mat(k,310) = -(rxt(k,191)*y(k,17) + rxt(k,192)*y(k,87) + rxt(k,193)*y(k,70))
         mat(k,722) = -rxt(k,191)*y(k,3)
         mat(k,391) = -rxt(k,192)*y(k,3)
         mat(k,676) = -rxt(k,193)*y(k,3)
         mat(k,700) = 4.000_r8*rxt(k,194)*y(k,5) + (rxt(k,195)+rxt(k,196))*y(k,28) &
                      + rxt(k,199)*y(k,61) + rxt(k,202)*y(k,68) + rxt(k,253)*y(k,77) &
                      + rxt(k,203)*y(k,96)
         mat(k,615) = (rxt(k,195)+rxt(k,196))*y(k,5)
         mat(k,171) = rxt(k,204)*y(k,68) + rxt(k,210)*y(k,92) + rxt(k,205)*y(k,96)
         mat(k,588) = rxt(k,199)*y(k,5)
         mat(k,487) = rxt(k,202)*y(k,5) + rxt(k,204)*y(k,45)
         mat(k,287) = rxt(k,253)*y(k,5)
         mat(k,543) = rxt(k,210)*y(k,45)
         mat(k,652) = rxt(k,203)*y(k,5) + rxt(k,205)*y(k,45)
         mat(k,692) = rxt(k,197)*y(k,28)
         mat(k,607) = rxt(k,197)*y(k,5)
         mat(k,322) = (rxt(k,276)+rxt(k,281))*y(k,53)
         mat(k,143) = (rxt(k,276)+rxt(k,281))*y(k,49)
         mat(k,715) = -(4._r8*rxt(k,194)*y(k,5) + (rxt(k,195) + rxt(k,196) + rxt(k,197) &
                      ) * y(k,28) + rxt(k,198)*y(k,87) + rxt(k,199)*y(k,61) + rxt(k,200) &
                      *y(k,62) + rxt(k,202)*y(k,68) + rxt(k,203)*y(k,96) + rxt(k,253) &
                      *y(k,77))
         mat(k,630) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,5)
         mat(k,406) = -rxt(k,198)*y(k,5)
         mat(k,602) = -rxt(k,199)*y(k,5)
         mat(k,459) = -rxt(k,200)*y(k,5)
         mat(k,502) = -rxt(k,202)*y(k,5)
         mat(k,667) = -rxt(k,203)*y(k,5)
         mat(k,297) = -rxt(k,253)*y(k,5)
         mat(k,318) = rxt(k,193)*y(k,70)
         mat(k,110) = rxt(k,201)*y(k,68)
         mat(k,175) = rxt(k,211)*y(k,92)
         mat(k,150) = rxt(k,206)*y(k,68)
         mat(k,502) = mat(k,502) + rxt(k,201)*y(k,6) + rxt(k,206)*y(k,53)
         mat(k,689) = rxt(k,193)*y(k,3)
         mat(k,558) = rxt(k,211)*y(k,45)
         mat(k,103) = -(rxt(k,201)*y(k,68))
         mat(k,467) = -rxt(k,201)*y(k,6)
         mat(k,694) = rxt(k,200)*y(k,62)
         mat(k,436) = rxt(k,200)*y(k,5)
         mat(k,738) = -(rxt(k,155)*y(k,25) + rxt(k,191)*y(k,3) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,68) + rxt(k,237)*y(k,96))
         mat(k,432) = -rxt(k,155)*y(k,17)
         mat(k,319) = -rxt(k,191)*y(k,17)
         mat(k,382) = -rxt(k,235)*y(k,17)
         mat(k,503) = -rxt(k,236)*y(k,17)
         mat(k,668) = -rxt(k,237)*y(k,17)
         mat(k,358) = rxt(k,162)*y(k,28) + rxt(k,239)*y(k,61)
         mat(k,81) = .300_r8*rxt(k,240)*y(k,96)
         mat(k,631) = rxt(k,162)*y(k,21)
         mat(k,603) = rxt(k,239)*y(k,21)
         mat(k,668) = mat(k,668) + .300_r8*rxt(k,240)*y(k,22)
         mat(k,347) = -(rxt(k,162)*y(k,28) + rxt(k,238)*y(k,87) + rxt(k,239)*y(k,61))
         mat(k,617) = -rxt(k,162)*y(k,21)
         mat(k,393) = -rxt(k,238)*y(k,21)
         mat(k,589) = -rxt(k,239)*y(k,21)
         mat(k,79) = .700_r8*rxt(k,240)*y(k,96)
         mat(k,654) = .700_r8*rxt(k,240)*y(k,22)
         mat(k,77) = -(rxt(k,240)*y(k,96))
         mat(k,638) = -rxt(k,240)*y(k,22)
         mat(k,345) = rxt(k,238)*y(k,87)
         mat(k,384) = rxt(k,238)*y(k,21)
         mat(k,421) = -(rxt(k,155)*y(k,17) + rxt(k,157)*y(k,41) + rxt(k,158)*y(k,43) &
                      + (rxt(k,159) + rxt(k,160)) * y(k,87) + rxt(k,161)*y(k,70) &
                      + rxt(k,168)*y(k,29) + rxt(k,177)*y(k,54))
         mat(k,727) = -rxt(k,155)*y(k,25)
         mat(k,569) = -rxt(k,157)*y(k,25)
         mat(k,127) = -rxt(k,158)*y(k,25)
         mat(k,396) = -(rxt(k,159) + rxt(k,160)) * y(k,25)
         mat(k,679) = -rxt(k,161)*y(k,25)
         mat(k,266) = -rxt(k,168)*y(k,25)
         mat(k,164) = -rxt(k,177)*y(k,25)
         mat(k,705) = rxt(k,196)*y(k,28)
         mat(k,349) = rxt(k,162)*y(k,28)
         mat(k,620) = rxt(k,196)*y(k,5) + rxt(k,162)*y(k,21) + (4.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,165))*y(k,28) + rxt(k,167)*y(k,61) + rxt(k,172) &
                      *y(k,68) + rxt(k,254)*y(k,77) + rxt(k,173)*y(k,96)
         mat(k,39) = rxt(k,217)*y(k,92)
         mat(k,333) = rxt(k,175)*y(k,68) + rxt(k,187)*y(k,92) + rxt(k,176)*y(k,96)
         mat(k,592) = rxt(k,167)*y(k,28)
         mat(k,492) = rxt(k,172)*y(k,28) + rxt(k,175)*y(k,49)
         mat(k,289) = rxt(k,254)*y(k,28)
         mat(k,548) = rxt(k,217)*y(k,34) + rxt(k,187)*y(k,49)
         mat(k,657) = rxt(k,173)*y(k,28) + rxt(k,176)*y(k,49)
         mat(k,409) = rxt(k,168)*y(k,29)
         mat(k,606) = 2.000_r8*rxt(k,164)*y(k,28)
         mat(k,260) = rxt(k,168)*y(k,25) + (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,49)
         mat(k,321) = (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,29) + (rxt(k,269) &
                       +rxt(k,275)+rxt(k,280))*y(k,54)
         mat(k,161) = (rxt(k,269)+rxt(k,275)+rxt(k,280))*y(k,49)
         mat(k,605) = 2.000_r8*rxt(k,189)*y(k,28)
         mat(k,627) = -(rxt(k,162)*y(k,21) + (4._r8*rxt(k,163) + 4._r8*rxt(k,164) &
                      + 4._r8*rxt(k,165) + 4._r8*rxt(k,189)) * y(k,28) + rxt(k,166) &
                      *y(k,87) + rxt(k,167)*y(k,61) + rxt(k,169)*y(k,62) + rxt(k,172) &
                      *y(k,68) + (rxt(k,173) + rxt(k,174)) * y(k,96) + (rxt(k,195) &
                      + rxt(k,196) + rxt(k,197)) * y(k,5) + rxt(k,254)*y(k,77))
         mat(k,355) = -rxt(k,162)*y(k,28)
         mat(k,403) = -rxt(k,166)*y(k,28)
         mat(k,599) = -rxt(k,167)*y(k,28)
         mat(k,456) = -rxt(k,169)*y(k,28)
         mat(k,499) = -rxt(k,172)*y(k,28)
         mat(k,664) = -(rxt(k,173) + rxt(k,174)) * y(k,28)
         mat(k,712) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,28)
         mat(k,294) = -rxt(k,254)*y(k,28)
         mat(k,428) = rxt(k,177)*y(k,54) + rxt(k,161)*y(k,70) + rxt(k,160)*y(k,87)
         mat(k,269) = rxt(k,170)*y(k,68)
         mat(k,339) = rxt(k,188)*y(k,92)
         mat(k,166) = rxt(k,177)*y(k,25) + rxt(k,178)*y(k,68) + rxt(k,179)*y(k,96)
         mat(k,499) = mat(k,499) + rxt(k,170)*y(k,29) + rxt(k,178)*y(k,54)
         mat(k,686) = rxt(k,161)*y(k,25)
         mat(k,68) = rxt(k,259)*y(k,77)
         mat(k,294) = mat(k,294) + rxt(k,259)*y(k,71)
         mat(k,403) = mat(k,403) + rxt(k,160)*y(k,25)
         mat(k,555) = rxt(k,188)*y(k,49)
         mat(k,664) = mat(k,664) + rxt(k,179)*y(k,54)
         mat(k,262) = -(rxt(k,168)*y(k,25) + rxt(k,170)*y(k,68) + rxt(k,171)*y(k,96) &
                      + (rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,49))
         mat(k,413) = -rxt(k,168)*y(k,29)
         mat(k,483) = -rxt(k,170)*y(k,29)
         mat(k,648) = -rxt(k,171)*y(k,29)
         mat(k,325) = -(rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,29)
         mat(k,611) = rxt(k,169)*y(k,62)
         mat(k,440) = rxt(k,169)*y(k,28)
         mat(k,111) = -((rxt(k,242) + rxt(k,246)) * y(k,96))
         mat(k,640) = -(rxt(k,242) + rxt(k,246)) * y(k,31)
         mat(k,307) = rxt(k,191)*y(k,17)
         mat(k,718) = rxt(k,191)*y(k,3) + rxt(k,155)*y(k,25) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,68) + rxt(k,237)*y(k,96)
         mat(k,410) = rxt(k,155)*y(k,17)
         mat(k,362) = rxt(k,235)*y(k,17)
         mat(k,468) = rxt(k,236)*y(k,17) + rxt(k,249)*y(k,72)
         mat(k,70) = rxt(k,249)*y(k,68) + rxt(k,250)*y(k,96)
         mat(k,640) = mat(k,640) + rxt(k,237)*y(k,17) + rxt(k,250)*y(k,72)
         mat(k,27) = -(rxt(k,216)*y(k,92))
         mat(k,536) = -rxt(k,216)*y(k,33)
         mat(k,37) = -(rxt(k,217)*y(k,92))
         mat(k,538) = -rxt(k,217)*y(k,34)
         mat(k,58) = -(rxt(k,247)*y(k,63) + (rxt(k,248) + rxt(k,261)) * y(k,96))
         mat(k,361) = -rxt(k,247)*y(k,35)
         mat(k,636) = -(rxt(k,248) + rxt(k,261)) * y(k,35)
         mat(k,153) = -(rxt(k,213)*y(k,41) + rxt(k,214)*y(k,100) + rxt(k,215)*y(k,51))
         mat(k,562) = -rxt(k,213)*y(k,39)
         mat(k,743) = -rxt(k,214)*y(k,39)
         mat(k,273) = -rxt(k,215)*y(k,39)
         mat(k,28) = 2.000_r8*rxt(k,216)*y(k,92)
         mat(k,38) = rxt(k,217)*y(k,92)
         mat(k,539) = 2.000_r8*rxt(k,216)*y(k,33) + rxt(k,217)*y(k,34)
         mat(k,299) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,87) + rxt(k,116) &
                      *y(k,69) + rxt(k,119)*y(k,70))
         mat(k,390) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,40)
         mat(k,520) = -rxt(k,116)*y(k,40)
         mat(k,675) = -rxt(k,119)*y(k,40)
         mat(k,721) = rxt(k,237)*y(k,96)
         mat(k,415) = rxt(k,157)*y(k,41)
         mat(k,112) = rxt(k,246)*y(k,96)
         mat(k,155) = rxt(k,213)*y(k,41)
         mat(k,564) = rxt(k,157)*y(k,25) + rxt(k,213)*y(k,39) + rxt(k,111)*y(k,68) &
                      + rxt(k,94)*y(k,92) + rxt(k,120)*y(k,96)
         mat(k,170) = rxt(k,211)*y(k,92)
         mat(k,327) = rxt(k,188)*y(k,92)
         mat(k,245) = rxt(k,143)*y(k,96)
         mat(k,486) = rxt(k,111)*y(k,41) + rxt(k,123)*y(k,96)
         mat(k,74) = rxt(k,250)*y(k,96)
         mat(k,133) = rxt(k,255)*y(k,96)
         mat(k,286) = rxt(k,260)*y(k,96)
         mat(k,542) = rxt(k,94)*y(k,41) + rxt(k,211)*y(k,45) + rxt(k,188)*y(k,49)
         mat(k,651) = rxt(k,237)*y(k,17) + rxt(k,246)*y(k,31) + rxt(k,120)*y(k,41) &
                      + rxt(k,143)*y(k,55) + rxt(k,123)*y(k,68) + rxt(k,250)*y(k,72) &
                      + rxt(k,255)*y(k,75) + rxt(k,260)*y(k,77)
         mat(k,574) = -(rxt(k,94)*y(k,92) + rxt(k,111)*y(k,68) + rxt(k,120)*y(k,96) &
                      + rxt(k,157)*y(k,25) + rxt(k,213)*y(k,39))
         mat(k,553) = -rxt(k,94)*y(k,41)
         mat(k,497) = -rxt(k,111)*y(k,41)
         mat(k,662) = -rxt(k,120)*y(k,41)
         mat(k,426) = -rxt(k,157)*y(k,41)
         mat(k,158) = -rxt(k,213)*y(k,41)
         mat(k,303) = rxt(k,113)*y(k,87)
         mat(k,401) = rxt(k,113)*y(k,40)
         mat(k,124) = -(rxt(k,112)*y(k,68) + rxt(k,121)*y(k,96) + rxt(k,158)*y(k,25))
         mat(k,470) = -rxt(k,112)*y(k,43)
         mat(k,641) = -rxt(k,121)*y(k,43)
         mat(k,411) = -rxt(k,158)*y(k,43)
         mat(k,386) = 2.000_r8*rxt(k,127)*y(k,87)
         mat(k,641) = mat(k,641) + 2.000_r8*rxt(k,126)*y(k,96)
         mat(k,47) = rxt(k,263)*y(k,100)
         mat(k,740) = rxt(k,263)*y(k,79)
         mat(k,169) = -(rxt(k,204)*y(k,68) + rxt(k,205)*y(k,96) + (rxt(k,210) &
                      + rxt(k,211)) * y(k,92))
         mat(k,475) = -rxt(k,204)*y(k,45)
         mat(k,645) = -rxt(k,205)*y(k,45)
         mat(k,540) = -(rxt(k,210) + rxt(k,211)) * y(k,45)
         mat(k,308) = rxt(k,191)*y(k,17) + rxt(k,192)*y(k,87)
         mat(k,719) = rxt(k,191)*y(k,3)
         mat(k,389) = rxt(k,192)*y(k,3)
         mat(k,329) = -(rxt(k,175)*y(k,68) + rxt(k,176)*y(k,96) + (rxt(k,187) &
                      + rxt(k,188)) * y(k,92) + (rxt(k,269) + rxt(k,275) + rxt(k,280) &
                      ) * y(k,54) + (rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,29) &
                      + (rxt(k,276) + rxt(k,281)) * y(k,53))
         mat(k,488) = -rxt(k,175)*y(k,49)
         mat(k,653) = -rxt(k,176)*y(k,49)
         mat(k,544) = -(rxt(k,187) + rxt(k,188)) * y(k,49)
         mat(k,163) = -(rxt(k,269) + rxt(k,275) + rxt(k,280)) * y(k,49)
         mat(k,264) = -(rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,49)
         mat(k,146) = -(rxt(k,276) + rxt(k,281)) * y(k,49)
         mat(k,723) = rxt(k,155)*y(k,25)
         mat(k,417) = rxt(k,155)*y(k,17) + rxt(k,157)*y(k,41) + rxt(k,158)*y(k,43) &
                      + rxt(k,177)*y(k,54) + rxt(k,159)*y(k,87)
         mat(k,616) = rxt(k,174)*y(k,96)
         mat(k,565) = rxt(k,157)*y(k,25)
         mat(k,125) = rxt(k,158)*y(k,25)
         mat(k,163) = mat(k,163) + rxt(k,177)*y(k,25)
         mat(k,392) = rxt(k,159)*y(k,25)
         mat(k,653) = mat(k,653) + rxt(k,174)*y(k,28)
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
         mat(k,152) = rxt(k,213)*y(k,41) + rxt(k,215)*y(k,51) + rxt(k,214)*y(k,100)
         mat(k,561) = rxt(k,213)*y(k,39)
         mat(k,272) = rxt(k,215)*y(k,39)
         mat(k,741) = rxt(k,214)*y(k,39)
         mat(k,274) = -(rxt(k,152)*y(k,96) + rxt(k,215)*y(k,39))
         mat(k,649) = -rxt(k,152)*y(k,51)
         mat(k,154) = -rxt(k,215)*y(k,51)
         mat(k,720) = rxt(k,235)*y(k,63)
         mat(k,263) = (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,49)
         mat(k,60) = rxt(k,247)*y(k,63)
         mat(k,326) = (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,29)
         mat(k,441) = rxt(k,151)*y(k,96)
         mat(k,364) = rxt(k,235)*y(k,17) + rxt(k,247)*y(k,35)
         mat(k,649) = mat(k,649) + rxt(k,151)*y(k,62)
         mat(k,83) = -(rxt(k,128)*y(k,96))
         mat(k,639) = -rxt(k,128)*y(k,52)
         mat(k,435) = rxt(k,149)*y(k,87)
         mat(k,385) = rxt(k,149)*y(k,62)
         mat(k,144) = -(rxt(k,206)*y(k,68) + (rxt(k,276) + rxt(k,281)) * y(k,49))
         mat(k,473) = -rxt(k,206)*y(k,53)
         mat(k,323) = -(rxt(k,276) + rxt(k,281)) * y(k,53)
         mat(k,695) = rxt(k,198)*y(k,87)
         mat(k,387) = rxt(k,198)*y(k,5)
         mat(k,162) = -(rxt(k,177)*y(k,25) + rxt(k,178)*y(k,68) + rxt(k,179)*y(k,96) &
                      + (rxt(k,269) + rxt(k,275) + rxt(k,280)) * y(k,49))
         mat(k,412) = -rxt(k,177)*y(k,54)
         mat(k,474) = -rxt(k,178)*y(k,54)
         mat(k,644) = -rxt(k,179)*y(k,54)
         mat(k,324) = -(rxt(k,269) + rxt(k,275) + rxt(k,280)) * y(k,54)
         mat(k,609) = rxt(k,166)*y(k,87)
         mat(k,261) = rxt(k,171)*y(k,96)
         mat(k,388) = rxt(k,166)*y(k,28)
         mat(k,644) = mat(k,644) + rxt(k,171)*y(k,29)
         mat(k,244) = -(rxt(k,131)*y(k,61) + (rxt(k,132) + rxt(k,133) + rxt(k,134) &
                      ) * y(k,62) + rxt(k,135)*y(k,69) + rxt(k,143)*y(k,96) + rxt(k,297) &
                      *y(k,95))
         mat(k,586) = -rxt(k,131)*y(k,55)
         mat(k,438) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,55)
         mat(k,517) = -rxt(k,135)*y(k,55)
         mat(k,646) = -rxt(k,143)*y(k,55)
         mat(k,203) = -rxt(k,297)*y(k,55)
         mat(k,481) = rxt(k,129)*y(k,88) + rxt(k,294)*y(k,91)
         mat(k,517) = mat(k,517) + rxt(k,295)*y(k,91)
         mat(k,217) = 1.100_r8*rxt(k,290)*y(k,89) + .200_r8*rxt(k,288)*y(k,90)
         mat(k,230) = rxt(k,129)*y(k,68)
         mat(k,121) = 1.100_r8*rxt(k,290)*y(k,86)
         mat(k,141) = .200_r8*rxt(k,288)*y(k,86)
         mat(k,183) = rxt(k,294)*y(k,68) + rxt(k,295)*y(k,69)
         mat(k,434) = rxt(k,150)*y(k,63)
         mat(k,360) = rxt(k,150)*y(k,62)
         mat(k,598) = -(rxt(k,131)*y(k,55) + rxt(k,140)*y(k,63) + rxt(k,144)*y(k,87) &
                      + rxt(k,145)*y(k,70) + rxt(k,146)*y(k,68) + rxt(k,167)*y(k,28) &
                      + rxt(k,199)*y(k,5) + rxt(k,239)*y(k,21) + rxt(k,299)*y(k,95))
         mat(k,251) = -rxt(k,131)*y(k,61)
         mat(k,377) = -rxt(k,140)*y(k,61)
         mat(k,402) = -rxt(k,144)*y(k,61)
         mat(k,685) = -rxt(k,145)*y(k,61)
         mat(k,498) = -rxt(k,146)*y(k,61)
         mat(k,626) = -rxt(k,167)*y(k,61)
         mat(k,711) = -rxt(k,199)*y(k,61)
         mat(k,354) = -rxt(k,239)*y(k,61)
         mat(k,207) = -rxt(k,299)*y(k,61)
         mat(k,251) = mat(k,251) + 2.000_r8*rxt(k,133)*y(k,62) + rxt(k,135)*y(k,69) &
                      + rxt(k,143)*y(k,96)
         mat(k,455) = 2.000_r8*rxt(k,133)*y(k,55) + rxt(k,136)*y(k,68) + rxt(k,256) &
                      *y(k,77)
         mat(k,498) = mat(k,498) + rxt(k,136)*y(k,62)
         mat(k,529) = rxt(k,135)*y(k,55) + rxt(k,130)*y(k,88)
         mat(k,293) = rxt(k,256)*y(k,62)
         mat(k,236) = rxt(k,130)*y(k,69)
         mat(k,663) = rxt(k,143)*y(k,55)
         mat(k,450) = -((rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,55) + (rxt(k,136) &
                      + rxt(k,138)) * y(k,68) + rxt(k,137)*y(k,70) + rxt(k,149) &
                      *y(k,87) + rxt(k,150)*y(k,63) + rxt(k,151)*y(k,96) + rxt(k,169) &
                      *y(k,28) + rxt(k,200)*y(k,5) + rxt(k,256)*y(k,77))
         mat(k,247) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,62)
         mat(k,493) = -(rxt(k,136) + rxt(k,138)) * y(k,62)
         mat(k,680) = -rxt(k,137)*y(k,62)
         mat(k,397) = -rxt(k,149)*y(k,62)
         mat(k,372) = -rxt(k,150)*y(k,62)
         mat(k,658) = -rxt(k,151)*y(k,62)
         mat(k,621) = -rxt(k,169)*y(k,62)
         mat(k,706) = -rxt(k,200)*y(k,62)
         mat(k,290) = -rxt(k,256)*y(k,62)
         mat(k,706) = mat(k,706) + rxt(k,199)*y(k,61)
         mat(k,350) = rxt(k,239)*y(k,61)
         mat(k,621) = mat(k,621) + rxt(k,167)*y(k,61)
         mat(k,86) = rxt(k,128)*y(k,96)
         mat(k,593) = rxt(k,199)*y(k,5) + rxt(k,239)*y(k,21) + rxt(k,167)*y(k,28) &
                      + 2.000_r8*rxt(k,140)*y(k,63) + rxt(k,146)*y(k,68) + rxt(k,145) &
                      *y(k,70) + rxt(k,144)*y(k,87)
         mat(k,372) = mat(k,372) + 2.000_r8*rxt(k,140)*y(k,61) + rxt(k,141)*y(k,68) &
                      + rxt(k,139)*y(k,87) + rxt(k,142)*y(k,96)
         mat(k,493) = mat(k,493) + rxt(k,146)*y(k,61) + rxt(k,141)*y(k,63)
         mat(k,680) = mat(k,680) + rxt(k,145)*y(k,61)
         mat(k,397) = mat(k,397) + rxt(k,144)*y(k,61) + rxt(k,139)*y(k,63)
         mat(k,658) = mat(k,658) + rxt(k,128)*y(k,52) + rxt(k,142)*y(k,63)
         mat(k,369) = -(rxt(k,139)*y(k,87) + rxt(k,140)*y(k,61) + rxt(k,141)*y(k,68) &
                      + rxt(k,142)*y(k,96) + rxt(k,150)*y(k,62) + rxt(k,235)*y(k,17) &
                      + rxt(k,247)*y(k,35))
         mat(k,394) = -rxt(k,139)*y(k,63)
         mat(k,590) = -rxt(k,140)*y(k,63)
         mat(k,490) = -rxt(k,141)*y(k,63)
         mat(k,655) = -rxt(k,142)*y(k,63)
         mat(k,447) = -rxt(k,150)*y(k,63)
         mat(k,725) = -rxt(k,235)*y(k,63)
         mat(k,61) = -rxt(k,247)*y(k,63)
         mat(k,107) = rxt(k,201)*y(k,68)
         mat(k,419) = rxt(k,168)*y(k,29)
         mat(k,265) = rxt(k,168)*y(k,25) + rxt(k,170)*y(k,68) + rxt(k,171)*y(k,96)
         mat(k,157) = rxt(k,215)*y(k,51)
         mat(k,277) = rxt(k,215)*y(k,39) + rxt(k,152)*y(k,96)
         mat(k,447) = mat(k,447) + rxt(k,138)*y(k,68) + rxt(k,137)*y(k,70)
         mat(k,490) = mat(k,490) + rxt(k,201)*y(k,6) + rxt(k,170)*y(k,29) + rxt(k,138) &
                      *y(k,62)
         mat(k,677) = rxt(k,137)*y(k,62)
         mat(k,655) = mat(k,655) + rxt(k,171)*y(k,29) + rxt(k,152)*y(k,51)
         mat(k,494) = -(rxt(k,108)*y(k,70) + 4._r8*rxt(k,109)*y(k,68) + rxt(k,110) &
                      *y(k,69) + rxt(k,111)*y(k,41) + rxt(k,112)*y(k,43) + rxt(k,117) &
                      *y(k,87) + rxt(k,123)*y(k,96) + (rxt(k,136) + rxt(k,138) &
                      ) * y(k,62) + rxt(k,141)*y(k,63) + rxt(k,146)*y(k,61) + rxt(k,170) &
                      *y(k,29) + rxt(k,172)*y(k,28) + rxt(k,175)*y(k,49) + rxt(k,178) &
                      *y(k,54) + rxt(k,201)*y(k,6) + rxt(k,202)*y(k,5) + rxt(k,204) &
                      *y(k,45) + rxt(k,206)*y(k,53) + rxt(k,236)*y(k,17) + rxt(k,249) &
                      *y(k,72) + (rxt(k,292) + rxt(k,293)) * y(k,89) + rxt(k,294) &
                      *y(k,91))
         mat(k,681) = -rxt(k,108)*y(k,68)
         mat(k,525) = -rxt(k,110)*y(k,68)
         mat(k,571) = -rxt(k,111)*y(k,68)
         mat(k,128) = -rxt(k,112)*y(k,68)
         mat(k,398) = -rxt(k,117)*y(k,68)
         mat(k,659) = -rxt(k,123)*y(k,68)
         mat(k,451) = -(rxt(k,136) + rxt(k,138)) * y(k,68)
         mat(k,373) = -rxt(k,141)*y(k,68)
         mat(k,594) = -rxt(k,146)*y(k,68)
         mat(k,268) = -rxt(k,170)*y(k,68)
         mat(k,622) = -rxt(k,172)*y(k,68)
         mat(k,335) = -rxt(k,175)*y(k,68)
         mat(k,165) = -rxt(k,178)*y(k,68)
         mat(k,109) = -rxt(k,201)*y(k,68)
         mat(k,707) = -rxt(k,202)*y(k,68)
         mat(k,172) = -rxt(k,204)*y(k,68)
         mat(k,148) = -rxt(k,206)*y(k,68)
         mat(k,729) = -rxt(k,236)*y(k,68)
         mat(k,75) = -rxt(k,249)*y(k,68)
         mat(k,122) = -(rxt(k,292) + rxt(k,293)) * y(k,68)
         mat(k,184) = -rxt(k,294)*y(k,68)
         mat(k,301) = rxt(k,115)*y(k,87)
         mat(k,248) = rxt(k,131)*y(k,61) + rxt(k,132)*y(k,62) + rxt(k,135)*y(k,69) &
                      + rxt(k,297)*y(k,95)
         mat(k,594) = mat(k,594) + rxt(k,131)*y(k,55)
         mat(k,451) = mat(k,451) + rxt(k,132)*y(k,55)
         mat(k,525) = mat(k,525) + rxt(k,135)*y(k,55) + rxt(k,251)*y(k,75) &
                      + rxt(k,257)*y(k,77) + rxt(k,296)*y(k,91) + (rxt(k,97)+rxt(k,98)) &
                      *y(k,92) + rxt(k,303)*y(k,97) + rxt(k,307)*y(k,98)
         mat(k,134) = rxt(k,251)*y(k,69)
         mat(k,291) = rxt(k,257)*y(k,69)
         mat(k,220) = rxt(k,288)*y(k,90) + 1.150_r8*rxt(k,289)*y(k,95)
         mat(k,398) = mat(k,398) + rxt(k,115)*y(k,40)
         mat(k,233) = rxt(k,302)*y(k,97)
         mat(k,142) = rxt(k,288)*y(k,86)
         mat(k,184) = mat(k,184) + rxt(k,296)*y(k,69)
         mat(k,550) = (rxt(k,97)+rxt(k,98))*y(k,69)
         mat(k,204) = rxt(k,297)*y(k,55) + 1.150_r8*rxt(k,289)*y(k,86)
         mat(k,659) = mat(k,659) + 2.000_r8*rxt(k,125)*y(k,96)
         mat(k,196) = rxt(k,303)*y(k,69) + rxt(k,302)*y(k,88)
         mat(k,101) = rxt(k,307)*y(k,69)
         mat(k,526) = -(rxt(k,97)*y(k,92) + rxt(k,102)*y(k,93) + rxt(k,110)*y(k,68) &
                      + rxt(k,116)*y(k,40) + rxt(k,130)*y(k,88) + rxt(k,135)*y(k,55) &
                      + rxt(k,251)*y(k,75) + rxt(k,257)*y(k,77) + rxt(k,291)*y(k,89) &
                      + (rxt(k,295) + rxt(k,296)) * y(k,91) + rxt(k,303)*y(k,97) &
                      + rxt(k,307)*y(k,98))
         mat(k,551) = -rxt(k,97)*y(k,69)
         mat(k,31) = -rxt(k,102)*y(k,69)
         mat(k,495) = -rxt(k,110)*y(k,69)
         mat(k,302) = -rxt(k,116)*y(k,69)
         mat(k,234) = -rxt(k,130)*y(k,69)
         mat(k,249) = -rxt(k,135)*y(k,69)
         mat(k,135) = -rxt(k,251)*y(k,69)
         mat(k,292) = -rxt(k,257)*y(k,69)
         mat(k,123) = -rxt(k,291)*y(k,69)
         mat(k,185) = -(rxt(k,295) + rxt(k,296)) * y(k,69)
         mat(k,197) = -rxt(k,303)*y(k,69)
         mat(k,102) = -rxt(k,307)*y(k,69)
         mat(k,313) = rxt(k,193)*y(k,70) + rxt(k,192)*y(k,87)
         mat(k,708) = 2.000_r8*rxt(k,194)*y(k,5) + (rxt(k,196)+rxt(k,197))*y(k,28) &
                      + rxt(k,202)*y(k,68) + rxt(k,198)*y(k,87)
         mat(k,352) = rxt(k,238)*y(k,87)
         mat(k,424) = rxt(k,161)*y(k,70) + rxt(k,159)*y(k,87)
         mat(k,623) = (rxt(k,196)+rxt(k,197))*y(k,5) + (2.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,164))*y(k,28) + rxt(k,172)*y(k,68) + rxt(k,166) &
                      *y(k,87) + rxt(k,174)*y(k,96)
         mat(k,302) = mat(k,302) + rxt(k,119)*y(k,70) + rxt(k,113)*y(k,87)
         mat(k,87) = rxt(k,128)*y(k,96)
         mat(k,249) = mat(k,249) + rxt(k,134)*y(k,62)
         mat(k,595) = rxt(k,145)*y(k,70) + rxt(k,299)*y(k,95)
         mat(k,452) = rxt(k,134)*y(k,55) + rxt(k,136)*y(k,68) + rxt(k,137)*y(k,70)
         mat(k,374) = rxt(k,141)*y(k,68) + rxt(k,139)*y(k,87)
         mat(k,495) = mat(k,495) + rxt(k,202)*y(k,5) + rxt(k,172)*y(k,28) + rxt(k,136) &
                      *y(k,62) + rxt(k,141)*y(k,63) + 2.000_r8*rxt(k,109)*y(k,68) &
                      + 2.000_r8*rxt(k,108)*y(k,70) + rxt(k,117)*y(k,87) + rxt(k,101) &
                      *y(k,93) + rxt(k,123)*y(k,96)
         mat(k,526) = mat(k,526) + 2.000_r8*rxt(k,102)*y(k,93)
         mat(k,682) = rxt(k,193)*y(k,3) + rxt(k,161)*y(k,25) + rxt(k,119)*y(k,40) &
                      + rxt(k,145)*y(k,61) + rxt(k,137)*y(k,62) + 2.000_r8*rxt(k,108) &
                      *y(k,68) + rxt(k,252)*y(k,75) + rxt(k,258)*y(k,77) &
                      + 2.000_r8*rxt(k,118)*y(k,87) + 2.000_r8*rxt(k,99)*y(k,92) &
                      + rxt(k,124)*y(k,96)
         mat(k,135) = mat(k,135) + rxt(k,252)*y(k,70)
         mat(k,292) = mat(k,292) + rxt(k,258)*y(k,70)
         mat(k,399) = rxt(k,192)*y(k,3) + rxt(k,198)*y(k,5) + rxt(k,238)*y(k,21) &
                      + rxt(k,159)*y(k,25) + rxt(k,166)*y(k,28) + rxt(k,113)*y(k,40) &
                      + rxt(k,139)*y(k,63) + rxt(k,117)*y(k,68) + 2.000_r8*rxt(k,118) &
                      *y(k,70) + 2.000_r8*rxt(k,127)*y(k,87) + rxt(k,122)*y(k,96)
         mat(k,551) = mat(k,551) + 2.000_r8*rxt(k,99)*y(k,70)
         mat(k,31) = mat(k,31) + rxt(k,101)*y(k,68) + 2.000_r8*rxt(k,102)*y(k,69)
         mat(k,205) = rxt(k,299)*y(k,61)
         mat(k,660) = rxt(k,174)*y(k,28) + rxt(k,128)*y(k,52) + rxt(k,123)*y(k,68) &
                      + rxt(k,124)*y(k,70) + rxt(k,122)*y(k,87)
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
         mat(k,688) = -(rxt(k,99)*y(k,92) + rxt(k,108)*y(k,68) + rxt(k,118)*y(k,87) &
                      + rxt(k,119)*y(k,40) + rxt(k,124)*y(k,96) + rxt(k,137)*y(k,62) &
                      + rxt(k,145)*y(k,61) + rxt(k,161)*y(k,25) + rxt(k,193)*y(k,3) &
                      + rxt(k,252)*y(k,75) + rxt(k,258)*y(k,77))
         mat(k,557) = -rxt(k,99)*y(k,70)
         mat(k,501) = -rxt(k,108)*y(k,70)
         mat(k,405) = -rxt(k,118)*y(k,70)
         mat(k,305) = -rxt(k,119)*y(k,70)
         mat(k,666) = -rxt(k,124)*y(k,70)
         mat(k,458) = -rxt(k,137)*y(k,70)
         mat(k,601) = -rxt(k,145)*y(k,70)
         mat(k,430) = -rxt(k,161)*y(k,70)
         mat(k,317) = -rxt(k,193)*y(k,70)
         mat(k,137) = -rxt(k,252)*y(k,70)
         mat(k,296) = -rxt(k,258)*y(k,70)
         mat(k,501) = mat(k,501) + rxt(k,110)*y(k,69)
         mat(k,532) = rxt(k,110)*y(k,68)
         mat(k,64) = -(rxt(k,259)*y(k,77))
         mat(k,282) = -rxt(k,259)*y(k,71)
         mat(k,693) = rxt(k,195)*y(k,28)
         mat(k,608) = rxt(k,195)*y(k,5) + 2.000_r8*rxt(k,165)*y(k,28)
         mat(k,69) = -(rxt(k,249)*y(k,68) + rxt(k,250)*y(k,96))
         mat(k,464) = -rxt(k,249)*y(k,72)
         mat(k,637) = -rxt(k,250)*y(k,72)
         mat(k,131) = -(rxt(k,251)*y(k,69) + rxt(k,252)*y(k,70) + rxt(k,255)*y(k,96))
         mat(k,510) = -rxt(k,251)*y(k,75)
         mat(k,672) = -rxt(k,252)*y(k,75)
         mat(k,642) = -rxt(k,255)*y(k,75)
         mat(k,285) = -(rxt(k,253)*y(k,5) + rxt(k,254)*y(k,28) + rxt(k,256)*y(k,62) &
                      + rxt(k,257)*y(k,69) + rxt(k,258)*y(k,70) + rxt(k,259)*y(k,71) &
                      + rxt(k,260)*y(k,96))
         mat(k,698) = -rxt(k,253)*y(k,77)
         mat(k,613) = -rxt(k,254)*y(k,77)
         mat(k,442) = -rxt(k,256)*y(k,77)
         mat(k,519) = -rxt(k,257)*y(k,77)
         mat(k,674) = -rxt(k,258)*y(k,77)
         mat(k,66) = -rxt(k,259)*y(k,77)
         mat(k,650) = -rxt(k,260)*y(k,77)
         mat(k,485) = rxt(k,249)*y(k,72)
         mat(k,519) = mat(k,519) + rxt(k,251)*y(k,75)
         mat(k,674) = mat(k,674) + rxt(k,252)*y(k,75)
         mat(k,73) = rxt(k,249)*y(k,68)
         mat(k,132) = rxt(k,251)*y(k,69) + rxt(k,252)*y(k,70) + rxt(k,255)*y(k,96)
         mat(k,650) = mat(k,650) + rxt(k,255)*y(k,75)
         mat(k,254) = -(rxt(k,262)*y(k,96))
         mat(k,647) = -rxt(k,262)*y(k,78)
         mat(k,696) = rxt(k,253)*y(k,77)
         mat(k,610) = rxt(k,254)*y(k,77)
         mat(k,59) = rxt(k,247)*y(k,63) + (rxt(k,248)+.500_r8*rxt(k,261))*y(k,96)
         mat(k,439) = rxt(k,256)*y(k,77)
         mat(k,363) = rxt(k,247)*y(k,35)
         mat(k,518) = rxt(k,257)*y(k,77)
         mat(k,673) = rxt(k,258)*y(k,77)
         mat(k,65) = rxt(k,259)*y(k,77)
         mat(k,72) = rxt(k,250)*y(k,96)
         mat(k,284) = rxt(k,253)*y(k,5) + rxt(k,254)*y(k,28) + rxt(k,256)*y(k,62) &
                      + rxt(k,257)*y(k,69) + rxt(k,258)*y(k,70) + rxt(k,259)*y(k,71) &
                      + rxt(k,260)*y(k,96)
         mat(k,647) = mat(k,647) + (rxt(k,248)+.500_r8*rxt(k,261))*y(k,35) &
                      + rxt(k,250)*y(k,72) + rxt(k,260)*y(k,77)
         mat(k,48) = -(rxt(k,263)*y(k,100))
         mat(k,742) = -rxt(k,263)*y(k,79)
         mat(k,253) = rxt(k,262)*y(k,96)
         mat(k,635) = rxt(k,262)*y(k,78)
         mat(k,215) = -(rxt(k,288)*y(k,90) + rxt(k,289)*y(k,95) + rxt(k,290)*y(k,89))
         mat(k,139) = -rxt(k,288)*y(k,86)
         mat(k,201) = -rxt(k,289)*y(k,86)
         mat(k,119) = -rxt(k,290)*y(k,86)
         mat(k,395) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,40) + rxt(k,117) &
                      *y(k,68) + rxt(k,118)*y(k,70) + rxt(k,122)*y(k,96) &
                      + 4._r8*rxt(k,127)*y(k,87) + rxt(k,139)*y(k,63) + rxt(k,144) &
                      *y(k,61) + rxt(k,149)*y(k,62) + (rxt(k,159) + rxt(k,160) &
                      ) * y(k,25) + rxt(k,166)*y(k,28) + rxt(k,192)*y(k,3) + rxt(k,198) &
                      *y(k,5) + rxt(k,238)*y(k,21))
         mat(k,300) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,87)
         mat(k,491) = -rxt(k,117)*y(k,87)
         mat(k,678) = -rxt(k,118)*y(k,87)
         mat(k,656) = -rxt(k,122)*y(k,87)
         mat(k,370) = -rxt(k,139)*y(k,87)
         mat(k,591) = -rxt(k,144)*y(k,87)
         mat(k,448) = -rxt(k,149)*y(k,87)
         mat(k,420) = -(rxt(k,159) + rxt(k,160)) * y(k,87)
         mat(k,619) = -rxt(k,166)*y(k,87)
         mat(k,311) = -rxt(k,192)*y(k,87)
         mat(k,704) = -rxt(k,198)*y(k,87)
         mat(k,348) = -rxt(k,238)*y(k,87)
         mat(k,311) = mat(k,311) + rxt(k,191)*y(k,17)
         mat(k,704) = mat(k,704) + rxt(k,203)*y(k,96)
         mat(k,726) = rxt(k,191)*y(k,3) + rxt(k,155)*y(k,25) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,68)
         mat(k,348) = mat(k,348) + rxt(k,162)*y(k,28) + rxt(k,239)*y(k,61)
         mat(k,420) = mat(k,420) + rxt(k,155)*y(k,17) + rxt(k,158)*y(k,43)
         mat(k,619) = mat(k,619) + rxt(k,162)*y(k,21) + rxt(k,173)*y(k,96)
         mat(k,113) = rxt(k,242)*y(k,96)
         mat(k,62) = .500_r8*rxt(k,261)*y(k,96)
         mat(k,300) = mat(k,300) + rxt(k,116)*y(k,69)
         mat(k,126) = rxt(k,158)*y(k,25) + rxt(k,112)*y(k,68) + rxt(k,121)*y(k,96)
         mat(k,591) = mat(k,591) + rxt(k,239)*y(k,21)
         mat(k,370) = mat(k,370) + rxt(k,235)*y(k,17) + rxt(k,142)*y(k,96)
         mat(k,491) = mat(k,491) + rxt(k,236)*y(k,17) + rxt(k,112)*y(k,43)
         mat(k,522) = rxt(k,116)*y(k,40)
         mat(k,678) = mat(k,678) + rxt(k,124)*y(k,96)
         mat(k,256) = rxt(k,262)*y(k,96)
         mat(k,656) = mat(k,656) + rxt(k,203)*y(k,5) + rxt(k,173)*y(k,28) + rxt(k,242) &
                      *y(k,31) + .500_r8*rxt(k,261)*y(k,35) + rxt(k,121)*y(k,43) &
                      + rxt(k,142)*y(k,63) + rxt(k,124)*y(k,70) + rxt(k,262)*y(k,78)
         mat(k,229) = -(rxt(k,129)*y(k,68) + rxt(k,130)*y(k,69) + rxt(k,302)*y(k,97))
         mat(k,480) = -rxt(k,129)*y(k,88)
         mat(k,516) = -rxt(k,130)*y(k,88)
         mat(k,192) = -rxt(k,302)*y(k,88)
         mat(k,480) = mat(k,480) + rxt(k,292)*y(k,89)
         mat(k,216) = .900_r8*rxt(k,290)*y(k,89) + .800_r8*rxt(k,288)*y(k,90)
         mat(k,120) = rxt(k,292)*y(k,68) + .900_r8*rxt(k,290)*y(k,86)
         mat(k,140) = .800_r8*rxt(k,288)*y(k,86)
         mat(k,115) = -(rxt(k,290)*y(k,86) + rxt(k,291)*y(k,69) + (rxt(k,292) &
                      + rxt(k,293)) * y(k,68))
         mat(k,210) = -rxt(k,290)*y(k,89)
         mat(k,509) = -rxt(k,291)*y(k,89)
         mat(k,469) = -(rxt(k,292) + rxt(k,293)) * y(k,89)
         mat(k,138) = -(rxt(k,288)*y(k,86))
         mat(k,211) = -rxt(k,288)*y(k,90)
         mat(k,238) = rxt(k,297)*y(k,95)
         mat(k,582) = rxt(k,299)*y(k,95)
         mat(k,472) = rxt(k,292)*y(k,89)
         mat(k,511) = rxt(k,296)*y(k,91)
         mat(k,116) = rxt(k,292)*y(k,68)
         mat(k,177) = rxt(k,296)*y(k,69)
         mat(k,199) = rxt(k,297)*y(k,55) + rxt(k,299)*y(k,61)
         mat(k,178) = -(rxt(k,294)*y(k,68) + (rxt(k,295) + rxt(k,296)) * y(k,69))
         mat(k,476) = -rxt(k,294)*y(k,91)
         mat(k,512) = -(rxt(k,295) + rxt(k,296)) * y(k,91)
         mat(k,225) = rxt(k,302)*y(k,97)
         mat(k,188) = rxt(k,302)*y(k,88)
         mat(k,552) = -(rxt(k,94)*y(k,41) + rxt(k,95)*y(k,100) + (rxt(k,97) + rxt(k,98) &
                      ) * y(k,69) + rxt(k,99)*y(k,70) + (rxt(k,187) + rxt(k,188) &
                      ) * y(k,49) + (rxt(k,210) + rxt(k,211)) * y(k,45) + rxt(k,216) &
                      *y(k,33) + rxt(k,217)*y(k,34))
         mat(k,573) = -rxt(k,94)*y(k,92)
         mat(k,756) = -rxt(k,95)*y(k,92)
         mat(k,527) = -(rxt(k,97) + rxt(k,98)) * y(k,92)
         mat(k,683) = -rxt(k,99)*y(k,92)
         mat(k,337) = -(rxt(k,187) + rxt(k,188)) * y(k,92)
         mat(k,173) = -(rxt(k,210) + rxt(k,211)) * y(k,92)
         mat(k,29) = -rxt(k,216)*y(k,92)
         mat(k,40) = -rxt(k,217)*y(k,92)
         mat(k,527) = mat(k,527) + rxt(k,130)*y(k,88)
         mat(k,222) = .850_r8*rxt(k,289)*y(k,95)
         mat(k,235) = rxt(k,130)*y(k,69)
         mat(k,206) = .850_r8*rxt(k,289)*y(k,86)
         mat(k,30) = -(rxt(k,101)*y(k,68) + rxt(k,102)*y(k,69))
         mat(k,462) = -rxt(k,101)*y(k,93)
         mat(k,505) = -rxt(k,102)*y(k,93)
         mat(k,462) = mat(k,462) + rxt(k,105)*y(k,94)
         mat(k,505) = mat(k,505) + rxt(k,106)*y(k,94)
         mat(k,670) = rxt(k,107)*y(k,94)
         mat(k,32) = rxt(k,105)*y(k,68) + rxt(k,106)*y(k,69) + rxt(k,107)*y(k,70)
         mat(k,33) = -(rxt(k,105)*y(k,68) + rxt(k,106)*y(k,69) + rxt(k,107)*y(k,70))
         mat(k,463) = -rxt(k,105)*y(k,94)
         mat(k,506) = -rxt(k,106)*y(k,94)
         mat(k,671) = -rxt(k,107)*y(k,94)
         mat(k,506) = mat(k,506) + rxt(k,97)*y(k,92)
         mat(k,537) = rxt(k,97)*y(k,69)
         mat(k,200) = -(rxt(k,289)*y(k,86) + rxt(k,297)*y(k,55) + rxt(k,299)*y(k,61))
         mat(k,214) = -rxt(k,289)*y(k,95)
         mat(k,241) = -rxt(k,297)*y(k,95)
         mat(k,583) = -rxt(k,299)*y(k,95)
         mat(k,514) = rxt(k,291)*y(k,89) + rxt(k,295)*y(k,91) + rxt(k,303)*y(k,97) &
                      + rxt(k,307)*y(k,98)
         mat(k,118) = rxt(k,291)*y(k,69)
         mat(k,180) = rxt(k,295)*y(k,69)
         mat(k,190) = rxt(k,303)*y(k,69)
         mat(k,100) = rxt(k,307)*y(k,69)
         mat(k,665) = -(rxt(k,120)*y(k,41) + rxt(k,121)*y(k,43) + rxt(k,122)*y(k,87) &
                      + rxt(k,123)*y(k,68) + rxt(k,124)*y(k,70) + (4._r8*rxt(k,125) &
                      + 4._r8*rxt(k,126)) * y(k,96) + rxt(k,128)*y(k,52) + rxt(k,142) &
                      *y(k,63) + rxt(k,143)*y(k,55) + rxt(k,151)*y(k,62) + rxt(k,152) &
                      *y(k,51) + rxt(k,171)*y(k,29) + (rxt(k,173) + rxt(k,174) &
                      ) * y(k,28) + rxt(k,176)*y(k,49) + rxt(k,179)*y(k,54) + rxt(k,203) &
                      *y(k,5) + rxt(k,205)*y(k,45) + rxt(k,237)*y(k,17) + rxt(k,240) &
                      *y(k,22) + (rxt(k,242) + rxt(k,246)) * y(k,31) + (rxt(k,248) &
                      + rxt(k,261)) * y(k,35) + rxt(k,250)*y(k,72) + rxt(k,255) &
                      *y(k,75) + rxt(k,260)*y(k,77) + rxt(k,262)*y(k,78))
         mat(k,577) = -rxt(k,120)*y(k,96)
         mat(k,129) = -rxt(k,121)*y(k,96)
         mat(k,404) = -rxt(k,122)*y(k,96)
         mat(k,500) = -rxt(k,123)*y(k,96)
         mat(k,687) = -rxt(k,124)*y(k,96)
         mat(k,88) = -rxt(k,128)*y(k,96)
         mat(k,379) = -rxt(k,142)*y(k,96)
         mat(k,252) = -rxt(k,143)*y(k,96)
         mat(k,457) = -rxt(k,151)*y(k,96)
         mat(k,280) = -rxt(k,152)*y(k,96)
         mat(k,270) = -rxt(k,171)*y(k,96)
         mat(k,628) = -(rxt(k,173) + rxt(k,174)) * y(k,96)
         mat(k,340) = -rxt(k,176)*y(k,96)
         mat(k,167) = -rxt(k,179)*y(k,96)
         mat(k,713) = -rxt(k,203)*y(k,96)
         mat(k,174) = -rxt(k,205)*y(k,96)
         mat(k,735) = -rxt(k,237)*y(k,96)
         mat(k,80) = -rxt(k,240)*y(k,96)
         mat(k,114) = -(rxt(k,242) + rxt(k,246)) * y(k,96)
         mat(k,63) = -(rxt(k,248) + rxt(k,261)) * y(k,96)
         mat(k,76) = -rxt(k,250)*y(k,96)
         mat(k,136) = -rxt(k,255)*y(k,96)
         mat(k,295) = -rxt(k,260)*y(k,96)
         mat(k,258) = -rxt(k,262)*y(k,96)
         mat(k,735) = mat(k,735) + rxt(k,236)*y(k,68)
         mat(k,80) = mat(k,80) + .300_r8*rxt(k,240)*y(k,96)
         mat(k,429) = rxt(k,160)*y(k,87)
         mat(k,159) = rxt(k,214)*y(k,100)
         mat(k,304) = rxt(k,119)*y(k,70) + 2.000_r8*rxt(k,114)*y(k,87)
         mat(k,577) = mat(k,577) + rxt(k,111)*y(k,68) + rxt(k,94)*y(k,92)
         mat(k,129) = mat(k,129) + rxt(k,112)*y(k,68)
         mat(k,174) = mat(k,174) + rxt(k,204)*y(k,68) + rxt(k,210)*y(k,92)
         mat(k,340) = mat(k,340) + rxt(k,175)*y(k,68) + rxt(k,187)*y(k,92)
         mat(k,149) = rxt(k,206)*y(k,68)
         mat(k,167) = mat(k,167) + rxt(k,178)*y(k,68)
         mat(k,600) = rxt(k,144)*y(k,87)
         mat(k,379) = mat(k,379) + rxt(k,139)*y(k,87)
         mat(k,500) = mat(k,500) + rxt(k,236)*y(k,17) + rxt(k,111)*y(k,41) &
                      + rxt(k,112)*y(k,43) + rxt(k,204)*y(k,45) + rxt(k,175)*y(k,49) &
                      + rxt(k,206)*y(k,53) + rxt(k,178)*y(k,54) + rxt(k,117)*y(k,87)
         mat(k,687) = mat(k,687) + rxt(k,119)*y(k,40) + rxt(k,118)*y(k,87)
         mat(k,404) = mat(k,404) + rxt(k,160)*y(k,25) + 2.000_r8*rxt(k,114)*y(k,40) &
                      + rxt(k,144)*y(k,61) + rxt(k,139)*y(k,63) + rxt(k,117)*y(k,68) &
                      + rxt(k,118)*y(k,70)
         mat(k,556) = rxt(k,94)*y(k,41) + rxt(k,210)*y(k,45) + rxt(k,187)*y(k,49) &
                      + 2.000_r8*rxt(k,95)*y(k,100)
         mat(k,665) = mat(k,665) + .300_r8*rxt(k,240)*y(k,22)
         mat(k,760) = rxt(k,214)*y(k,39) + 2.000_r8*rxt(k,95)*y(k,92)
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
         mat(k,189) = -(rxt(k,302)*y(k,88) + rxt(k,303)*y(k,69))
         mat(k,226) = -rxt(k,302)*y(k,97)
         mat(k,513) = -rxt(k,303)*y(k,97)
         mat(k,477) = rxt(k,293)*y(k,89) + rxt(k,294)*y(k,91) + rxt(k,306)*y(k,98) &
                      + rxt(k,312)*y(k,99)
         mat(k,213) = rxt(k,304)*y(k,98) + rxt(k,309)*y(k,99)
         mat(k,117) = rxt(k,293)*y(k,68)
         mat(k,179) = rxt(k,294)*y(k,68)
         mat(k,99) = rxt(k,306)*y(k,68) + rxt(k,304)*y(k,86)
         mat(k,94) = rxt(k,312)*y(k,68) + rxt(k,309)*y(k,86)
         mat(k,97) = -(rxt(k,304)*y(k,86) + rxt(k,306)*y(k,68) + rxt(k,307)*y(k,69))
         mat(k,209) = -rxt(k,304)*y(k,98)
         mat(k,466) = -rxt(k,306)*y(k,98)
         mat(k,508) = -rxt(k,307)*y(k,98)
         mat(k,209) = mat(k,209) + rxt(k,308)*y(k,99)
         mat(k,91) = rxt(k,308)*y(k,86)
         mat(k,90) = -((rxt(k,308) + rxt(k,309)) * y(k,86) + rxt(k,312)*y(k,68))
         mat(k,208) = -(rxt(k,308) + rxt(k,309)) * y(k,99)
         mat(k,465) = -rxt(k,312)*y(k,99)
         mat(k,764) = -(rxt(k,95)*y(k,92) + rxt(k,214)*y(k,39) + rxt(k,263)*y(k,79))
         mat(k,560) = -rxt(k,95)*y(k,100)
         mat(k,160) = -rxt(k,214)*y(k,100)
         mat(k,51) = -rxt(k,263)*y(k,100)
         mat(k,739) = rxt(k,237)*y(k,96)
         mat(k,82) = rxt(k,240)*y(k,96)
         mat(k,306) = rxt(k,115)*y(k,87)
         mat(k,581) = rxt(k,120)*y(k,96)
         mat(k,130) = rxt(k,121)*y(k,96)
         mat(k,176) = rxt(k,205)*y(k,96)
         mat(k,344) = (rxt(k,276)+rxt(k,281))*y(k,53) + (rxt(k,269)+rxt(k,275) &
                       +rxt(k,280))*y(k,54) + rxt(k,176)*y(k,96)
         mat(k,281) = rxt(k,152)*y(k,96)
         mat(k,89) = rxt(k,128)*y(k,96)
         mat(k,151) = (rxt(k,276)+rxt(k,281))*y(k,49)
         mat(k,168) = (rxt(k,269)+rxt(k,275)+rxt(k,280))*y(k,49) + rxt(k,179)*y(k,96)
         mat(k,408) = rxt(k,115)*y(k,40) + rxt(k,122)*y(k,96)
         mat(k,669) = rxt(k,237)*y(k,17) + rxt(k,240)*y(k,22) + rxt(k,120)*y(k,41) &
                      + rxt(k,121)*y(k,43) + rxt(k,205)*y(k,45) + rxt(k,176)*y(k,49) &
                      + rxt(k,152)*y(k,51) + rxt(k,128)*y(k,52) + rxt(k,179)*y(k,54) &
                      + rxt(k,122)*y(k,87) + 2.000_r8*rxt(k,125)*y(k,96)
      end do
      end subroutine nlnmat04
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
         mat(k, 28) = mat(k, 28) + lmat(k, 28)
         mat(k, 30) = mat(k, 30) + lmat(k, 30)
         mat(k, 31) = mat(k, 31) + lmat(k, 31)
         mat(k, 32) = mat(k, 32) + lmat(k, 32)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 34) = lmat(k, 34)
         mat(k, 35) = lmat(k, 35)
         mat(k, 36) = lmat(k, 36)
         mat(k, 37) = mat(k, 37) + lmat(k, 37)
         mat(k, 38) = mat(k, 38) + lmat(k, 38)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 41) = lmat(k, 41)
         mat(k, 42) = lmat(k, 42)
         mat(k, 43) = lmat(k, 43)
         mat(k, 44) = lmat(k, 44)
         mat(k, 45) = lmat(k, 45)
         mat(k, 46) = lmat(k, 46)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 49) = lmat(k, 49)
         mat(k, 50) = lmat(k, 50)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = lmat(k, 53)
         mat(k, 54) = lmat(k, 54)
         mat(k, 55) = lmat(k, 55)
         mat(k, 56) = lmat(k, 56)
         mat(k, 57) = lmat(k, 57)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 71) = lmat(k, 71)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = lmat(k, 78)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 84) = lmat(k, 84)
         mat(k, 85) = lmat(k, 85)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 92) = lmat(k, 92)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = mat(k, 99) + lmat(k, 99)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = lmat(k, 104)
         mat(k, 105) = lmat(k, 105)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = lmat(k, 108)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 129) = mat(k, 129) + lmat(k, 129)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 145) = lmat(k, 145)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 156) = lmat(k, 156)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 186) = lmat(k, 186)
         mat(k, 187) = lmat(k, 187)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 193) = lmat(k, 193)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 239) = lmat(k, 239)
         mat(k, 242) = lmat(k, 242)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 255) = lmat(k, 255)
         mat(k, 257) = lmat(k, 257)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 267) = lmat(k, 267)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 274) = mat(k, 274) + lmat(k, 274)
         mat(k, 278) = lmat(k, 278)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 283) = lmat(k, 283)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 329) = mat(k, 329) + lmat(k, 329)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 372) = mat(k, 372) + lmat(k, 372)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 374) = mat(k, 374) + lmat(k, 374)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 386) = mat(k, 386) + lmat(k, 386)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 416) = lmat(k, 416)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 418) = lmat(k, 418)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 466) = mat(k, 466) + lmat(k, 466)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 479) = lmat(k, 479)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 515) = lmat(k, 515)
         mat(k, 525) = mat(k, 525) + lmat(k, 525)
         mat(k, 526) = mat(k, 526) + lmat(k, 526)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 538) = mat(k, 538) + lmat(k, 538)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 545) = lmat(k, 545)
         mat(k, 547) = lmat(k, 547)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 551) = mat(k, 551) + lmat(k, 551)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 553) = mat(k, 553) + lmat(k, 553)
         mat(k, 554) = lmat(k, 554)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 559) = lmat(k, 559)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 584) = lmat(k, 584)
         mat(k, 586) = mat(k, 586) + lmat(k, 586)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 633) = lmat(k, 633)
         mat(k, 634) = lmat(k, 634)
         mat(k, 652) = mat(k, 652) + lmat(k, 652)
         mat(k, 654) = mat(k, 654) + lmat(k, 654)
         mat(k, 656) = mat(k, 656) + lmat(k, 656)
         mat(k, 657) = mat(k, 657) + lmat(k, 657)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 681) = mat(k, 681) + lmat(k, 681)
         mat(k, 682) = mat(k, 682) + lmat(k, 682)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 688) = mat(k, 688) + lmat(k, 688)
         mat(k, 700) = mat(k, 700) + lmat(k, 700)
         mat(k, 707) = mat(k, 707) + lmat(k, 707)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 718) = mat(k, 718) + lmat(k, 718)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 732) = lmat(k, 732)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 747) = lmat(k, 747)
         mat(k, 754) = lmat(k, 754)
         mat(k, 756) = mat(k, 756) + lmat(k, 756)
         mat(k, 757) = lmat(k, 757)
         mat(k, 760) = mat(k, 760) + lmat(k, 760)
         mat(k, 764) = mat(k, 764) + lmat(k, 764)
         mat(k, 147) = 0._r8
         mat(k, 181) = 0._r8
         mat(k, 182) = 0._r8
         mat(k, 191) = 0._r8
         mat(k, 194) = 0._r8
         mat(k, 195) = 0._r8
         mat(k, 198) = 0._r8
         mat(k, 202) = 0._r8
         mat(k, 212) = 0._r8
         mat(k, 218) = 0._r8
         mat(k, 219) = 0._r8
         mat(k, 221) = 0._r8
         mat(k, 223) = 0._r8
         mat(k, 224) = 0._r8
         mat(k, 227) = 0._r8
         mat(k, 228) = 0._r8
         mat(k, 231) = 0._r8
         mat(k, 232) = 0._r8
         mat(k, 237) = 0._r8
         mat(k, 240) = 0._r8
         mat(k, 243) = 0._r8
         mat(k, 246) = 0._r8
         mat(k, 250) = 0._r8
         mat(k, 259) = 0._r8
         mat(k, 271) = 0._r8
         mat(k, 275) = 0._r8
         mat(k, 276) = 0._r8
         mat(k, 279) = 0._r8
         mat(k, 288) = 0._r8
         mat(k, 298) = 0._r8
         mat(k, 309) = 0._r8
         mat(k, 312) = 0._r8
         mat(k, 314) = 0._r8
         mat(k, 315) = 0._r8
         mat(k, 316) = 0._r8
         mat(k, 320) = 0._r8
         mat(k, 328) = 0._r8
         mat(k, 330) = 0._r8
         mat(k, 331) = 0._r8
         mat(k, 332) = 0._r8
         mat(k, 334) = 0._r8
         mat(k, 336) = 0._r8
         mat(k, 338) = 0._r8
         mat(k, 341) = 0._r8
         mat(k, 342) = 0._r8
         mat(k, 343) = 0._r8
         mat(k, 346) = 0._r8
         mat(k, 351) = 0._r8
         mat(k, 353) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 357) = 0._r8
         mat(k, 359) = 0._r8
         mat(k, 365) = 0._r8
         mat(k, 366) = 0._r8
         mat(k, 367) = 0._r8
         mat(k, 368) = 0._r8
         mat(k, 371) = 0._r8
         mat(k, 375) = 0._r8
         mat(k, 376) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 380) = 0._r8
         mat(k, 381) = 0._r8
         mat(k, 383) = 0._r8
         mat(k, 400) = 0._r8
         mat(k, 407) = 0._r8
         mat(k, 414) = 0._r8
         mat(k, 422) = 0._r8
         mat(k, 423) = 0._r8
         mat(k, 425) = 0._r8
         mat(k, 427) = 0._r8
         mat(k, 431) = 0._r8
         mat(k, 433) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 443) = 0._r8
         mat(k, 444) = 0._r8
         mat(k, 445) = 0._r8
         mat(k, 446) = 0._r8
         mat(k, 449) = 0._r8
         mat(k, 453) = 0._r8
         mat(k, 454) = 0._r8
         mat(k, 460) = 0._r8
         mat(k, 461) = 0._r8
         mat(k, 471) = 0._r8
         mat(k, 478) = 0._r8
         mat(k, 482) = 0._r8
         mat(k, 484) = 0._r8
         mat(k, 489) = 0._r8
         mat(k, 496) = 0._r8
         mat(k, 504) = 0._r8
         mat(k, 521) = 0._r8
         mat(k, 523) = 0._r8
         mat(k, 524) = 0._r8
         mat(k, 528) = 0._r8
         mat(k, 530) = 0._r8
         mat(k, 531) = 0._r8
         mat(k, 533) = 0._r8
         mat(k, 534) = 0._r8
         mat(k, 535) = 0._r8
         mat(k, 541) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 549) = 0._r8
         mat(k, 563) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 567) = 0._r8
         mat(k, 568) = 0._r8
         mat(k, 570) = 0._r8
         mat(k, 572) = 0._r8
         mat(k, 575) = 0._r8
         mat(k, 576) = 0._r8
         mat(k, 578) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 580) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 587) = 0._r8
         mat(k, 596) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 604) = 0._r8
         mat(k, 612) = 0._r8
         mat(k, 614) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 625) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 643) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 684) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 699) = 0._r8
         mat(k, 701) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 709) = 0._r8
         mat(k, 710) = 0._r8
         mat(k, 714) = 0._r8
         mat(k, 716) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 728) = 0._r8
         mat(k, 730) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 744) = 0._r8
         mat(k, 745) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 759) = 0._r8
         mat(k, 761) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 763) = 0._r8
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
         mat(k, 24) = mat(k, 24) - dti(k)
         mat(k, 27) = mat(k, 27) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 35) = mat(k, 35) - dti(k)
         mat(k, 37) = mat(k, 37) - dti(k)
         mat(k, 41) = mat(k, 41) - dti(k)
         mat(k, 44) = mat(k, 44) - dti(k)
         mat(k, 48) = mat(k, 48) - dti(k)
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 58) = mat(k, 58) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
         mat(k, 124) = mat(k, 124) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 153) = mat(k, 153) - dti(k)
         mat(k, 162) = mat(k, 162) - dti(k)
         mat(k, 169) = mat(k, 169) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 215) = mat(k, 215) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 254) = mat(k, 254) - dti(k)
         mat(k, 262) = mat(k, 262) - dti(k)
         mat(k, 274) = mat(k, 274) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 299) = mat(k, 299) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 329) = mat(k, 329) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 369) = mat(k, 369) - dti(k)
         mat(k, 395) = mat(k, 395) - dti(k)
         mat(k, 421) = mat(k, 421) - dti(k)
         mat(k, 450) = mat(k, 450) - dti(k)
         mat(k, 494) = mat(k, 494) - dti(k)
         mat(k, 526) = mat(k, 526) - dti(k)
         mat(k, 552) = mat(k, 552) - dti(k)
         mat(k, 574) = mat(k, 574) - dti(k)
         mat(k, 598) = mat(k, 598) - dti(k)
         mat(k, 627) = mat(k, 627) - dti(k)
         mat(k, 665) = mat(k, 665) - dti(k)
         mat(k, 688) = mat(k, 688) - dti(k)
         mat(k, 715) = mat(k, 715) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 764) = mat(k, 764) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
