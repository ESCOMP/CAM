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
         mat(k,291) = -(rxt(k,191)*y(k,15) + rxt(k,192)*y(k,68) + rxt(k,193)*y(k,58))
         mat(k,546) = -rxt(k,191)*y(k,1)
         mat(k,372) = -rxt(k,192)*y(k,1)
         mat(k,570) = -rxt(k,193)*y(k,1)
         mat(k,703) = 4.000_r8*rxt(k,194)*y(k,3) + (rxt(k,195)+rxt(k,196))*y(k,26) &
                      + rxt(k,199)*y(k,53) + rxt(k,202)*y(k,56) + rxt(k,253)*y(k,63) &
                      + rxt(k,203)*y(k,77)
         mat(k,443) = (rxt(k,195)+rxt(k,196))*y(k,3)
         mat(k,152) = rxt(k,204)*y(k,56) + rxt(k,210)*y(k,73) + rxt(k,205)*y(k,77)
         mat(k,657) = rxt(k,199)*y(k,3)
         mat(k,415) = rxt(k,202)*y(k,3) + rxt(k,204)*y(k,40)
         mat(k,268) = rxt(k,253)*y(k,3)
         mat(k,524) = rxt(k,210)*y(k,40)
         mat(k,633) = rxt(k,203)*y(k,3) + rxt(k,205)*y(k,40)
         mat(k,695) = rxt(k,197)*y(k,26)
         mat(k,435) = rxt(k,197)*y(k,3)
         mat(k,303) = (rxt(k,276)+rxt(k,281))*y(k,48)
         mat(k,133) = (rxt(k,276)+rxt(k,281))*y(k,44)
         mat(k,719) = -(4._r8*rxt(k,194)*y(k,3) + (rxt(k,195) + rxt(k,196) + rxt(k,197) &
                      ) * y(k,26) + rxt(k,198)*y(k,68) + rxt(k,199)*y(k,53) + rxt(k,200) &
                      *y(k,54) + rxt(k,202)*y(k,56) + rxt(k,203)*y(k,77) + rxt(k,253) &
                      *y(k,63))
         mat(k,459) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,3)
         mat(k,388) = -rxt(k,198)*y(k,3)
         mat(k,672) = -rxt(k,199)*y(k,3)
         mat(k,612) = -rxt(k,200)*y(k,3)
         mat(k,431) = -rxt(k,202)*y(k,3)
         mat(k,649) = -rxt(k,203)*y(k,3)
         mat(k,278) = -rxt(k,253)*y(k,3)
         mat(k,300) = rxt(k,193)*y(k,58)
         mat(k,91) = rxt(k,201)*y(k,56)
         mat(k,156) = rxt(k,211)*y(k,73)
         mat(k,140) = rxt(k,206)*y(k,56)
         mat(k,431) = mat(k,431) + rxt(k,201)*y(k,4) + rxt(k,206)*y(k,48)
         mat(k,584) = rxt(k,193)*y(k,1)
         mat(k,540) = rxt(k,211)*y(k,40)
         mat(k,84) = -(rxt(k,201)*y(k,56))
         mat(k,395) = -rxt(k,201)*y(k,4)
         mat(k,697) = rxt(k,200)*y(k,54)
         mat(k,588) = rxt(k,200)*y(k,3)
         mat(k,556) = -(rxt(k,155)*y(k,23) + rxt(k,191)*y(k,1) + rxt(k,235)*y(k,55) &
                      + rxt(k,236)*y(k,56) + rxt(k,237)*y(k,77))
         mat(k,478) = -rxt(k,155)*y(k,15)
         mat(k,296) = -rxt(k,191)*y(k,15)
         mat(k,357) = -rxt(k,235)*y(k,15)
         mat(k,425) = -rxt(k,236)*y(k,15)
         mat(k,643) = -rxt(k,237)*y(k,15)
         mat(k,334) = rxt(k,162)*y(k,26) + rxt(k,239)*y(k,53)
         mat(k,61) = .300_r8*rxt(k,240)*y(k,77)
         mat(k,453) = rxt(k,162)*y(k,19)
         mat(k,666) = rxt(k,239)*y(k,19)
         mat(k,643) = mat(k,643) + .300_r8*rxt(k,240)*y(k,20)
         mat(k,328) = -(rxt(k,162)*y(k,26) + rxt(k,238)*y(k,68) + rxt(k,239)*y(k,53))
         mat(k,445) = -rxt(k,162)*y(k,19)
         mat(k,374) = -rxt(k,238)*y(k,19)
         mat(k,658) = -rxt(k,239)*y(k,19)
         mat(k,60) = .700_r8*rxt(k,240)*y(k,77)
         mat(k,635) = .700_r8*rxt(k,240)*y(k,20)
         mat(k,58) = -(rxt(k,240)*y(k,77))
         mat(k,619) = -rxt(k,240)*y(k,20)
         mat(k,326) = rxt(k,238)*y(k,68)
         mat(k,365) = rxt(k,238)*y(k,19)
         mat(k,475) = -(rxt(k,155)*y(k,15) + rxt(k,157)*y(k,36) + rxt(k,158)*y(k,38) &
                      + (rxt(k,159) + rxt(k,160)) * y(k,68) + rxt(k,161)*y(k,58) &
                      + rxt(k,168)*y(k,27) + rxt(k,177)*y(k,49))
         mat(k,553) = -rxt(k,155)*y(k,23)
         mat(k,684) = -rxt(k,157)*y(k,23)
         mat(k,100) = -rxt(k,158)*y(k,23)
         mat(k,379) = -(rxt(k,159) + rxt(k,160)) * y(k,23)
         mat(k,575) = -rxt(k,161)*y(k,23)
         mat(k,249) = -rxt(k,168)*y(k,23)
         mat(k,147) = -rxt(k,177)*y(k,23)
         mat(k,710) = rxt(k,196)*y(k,26)
         mat(k,332) = rxt(k,162)*y(k,26)
         mat(k,450) = rxt(k,196)*y(k,3) + rxt(k,162)*y(k,19) + (4.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,165))*y(k,26) + rxt(k,167)*y(k,53) + rxt(k,172) &
                      *y(k,56) + rxt(k,254)*y(k,63) + rxt(k,173)*y(k,77)
         mat(k,20) = rxt(k,217)*y(k,73)
         mat(k,316) = rxt(k,175)*y(k,56) + rxt(k,187)*y(k,73) + rxt(k,176)*y(k,77)
         mat(k,663) = rxt(k,167)*y(k,26)
         mat(k,422) = rxt(k,172)*y(k,26) + rxt(k,175)*y(k,44)
         mat(k,272) = rxt(k,254)*y(k,26)
         mat(k,531) = rxt(k,217)*y(k,32) + rxt(k,187)*y(k,44)
         mat(k,640) = rxt(k,173)*y(k,26) + rxt(k,176)*y(k,44)
         mat(k,461) = rxt(k,168)*y(k,27)
         mat(k,434) = 2.000_r8*rxt(k,164)*y(k,26)
         mat(k,241) = rxt(k,168)*y(k,23) + (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,44)
         mat(k,302) = (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,27) + (rxt(k,269) &
                       +rxt(k,275)+rxt(k,280))*y(k,49)
         mat(k,142) = (rxt(k,269)+rxt(k,275)+rxt(k,280))*y(k,44)
         mat(k,433) = 2.000_r8*rxt(k,189)*y(k,26)
         mat(k,449) = -(rxt(k,162)*y(k,19) + (4._r8*rxt(k,163) + 4._r8*rxt(k,164) &
                      + 4._r8*rxt(k,165) + 4._r8*rxt(k,189)) * y(k,26) + rxt(k,166) &
                      *y(k,68) + rxt(k,167)*y(k,53) + rxt(k,169)*y(k,54) + rxt(k,172) &
                      *y(k,56) + (rxt(k,173) + rxt(k,174)) * y(k,77) + (rxt(k,195) &
                      + rxt(k,196) + rxt(k,197)) * y(k,3) + rxt(k,254)*y(k,63))
         mat(k,331) = -rxt(k,162)*y(k,26)
         mat(k,378) = -rxt(k,166)*y(k,26)
         mat(k,662) = -rxt(k,167)*y(k,26)
         mat(k,602) = -rxt(k,169)*y(k,26)
         mat(k,421) = -rxt(k,172)*y(k,26)
         mat(k,639) = -(rxt(k,173) + rxt(k,174)) * y(k,26)
         mat(k,709) = -(rxt(k,195) + rxt(k,196) + rxt(k,197)) * y(k,26)
         mat(k,271) = -rxt(k,254)*y(k,26)
         mat(k,474) = rxt(k,177)*y(k,49) + rxt(k,161)*y(k,58) + rxt(k,160)*y(k,68)
         mat(k,248) = rxt(k,170)*y(k,56)
         mat(k,315) = rxt(k,188)*y(k,73)
         mat(k,146) = rxt(k,177)*y(k,23) + rxt(k,178)*y(k,56) + rxt(k,179)*y(k,77)
         mat(k,421) = mat(k,421) + rxt(k,170)*y(k,27) + rxt(k,178)*y(k,49)
         mat(k,574) = rxt(k,161)*y(k,23)
         mat(k,49) = rxt(k,259)*y(k,63)
         mat(k,271) = mat(k,271) + rxt(k,259)*y(k,59)
         mat(k,378) = mat(k,378) + rxt(k,160)*y(k,23)
         mat(k,530) = rxt(k,188)*y(k,44)
         mat(k,639) = mat(k,639) + rxt(k,179)*y(k,49)
         mat(k,243) = -(rxt(k,168)*y(k,23) + rxt(k,170)*y(k,56) + rxt(k,171)*y(k,77) &
                      + (rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,44))
         mat(k,465) = -rxt(k,168)*y(k,27)
         mat(k,411) = -rxt(k,170)*y(k,27)
         mat(k,629) = -rxt(k,171)*y(k,27)
         mat(k,306) = -(rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,27)
         mat(k,439) = rxt(k,169)*y(k,54)
         mat(k,592) = rxt(k,169)*y(k,26)
         mat(k,92) = -((rxt(k,242) + rxt(k,246)) * y(k,77))
         mat(k,621) = -(rxt(k,242) + rxt(k,246)) * y(k,29)
         mat(k,288) = rxt(k,191)*y(k,15)
         mat(k,542) = rxt(k,191)*y(k,1) + rxt(k,155)*y(k,23) + rxt(k,235)*y(k,55) &
                      + rxt(k,236)*y(k,56) + rxt(k,237)*y(k,77)
         mat(k,462) = rxt(k,155)*y(k,15)
         mat(k,343) = rxt(k,235)*y(k,15)
         mat(k,396) = rxt(k,236)*y(k,15) + rxt(k,249)*y(k,60)
         mat(k,51) = rxt(k,249)*y(k,56) + rxt(k,250)*y(k,77)
         mat(k,621) = mat(k,621) + rxt(k,237)*y(k,15) + rxt(k,250)*y(k,60)
         mat(k,5) = -(rxt(k,216)*y(k,73))
         mat(k,517) = -rxt(k,216)*y(k,31)
         mat(k,18) = -(rxt(k,217)*y(k,73))
         mat(k,519) = -rxt(k,217)*y(k,32)
         mat(k,33) = -(rxt(k,247)*y(k,55) + (rxt(k,248) + rxt(k,261)) * y(k,77))
         mat(k,341) = -rxt(k,247)*y(k,33)
         mat(k,617) = -(rxt(k,248) + rxt(k,261)) * y(k,33)
         mat(k,125) = -(rxt(k,213)*y(k,36) + rxt(k,214)*y(k,81) + rxt(k,215)*y(k,46))
         mat(k,675) = -rxt(k,213)*y(k,34)
         mat(k,724) = -rxt(k,214)*y(k,34)
         mat(k,254) = -rxt(k,215)*y(k,34)
         mat(k,6) = 2.000_r8*rxt(k,216)*y(k,73)
         mat(k,19) = rxt(k,217)*y(k,73)
         mat(k,520) = 2.000_r8*rxt(k,216)*y(k,31) + rxt(k,217)*y(k,32)
         mat(k,280) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,68) + rxt(k,116) &
                      *y(k,57) + rxt(k,119)*y(k,58))
         mat(k,371) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,35)
         mat(k,501) = -rxt(k,116)*y(k,35)
         mat(k,569) = -rxt(k,119)*y(k,35)
         mat(k,545) = rxt(k,237)*y(k,77)
         mat(k,467) = rxt(k,157)*y(k,36)
         mat(k,93) = rxt(k,246)*y(k,77)
         mat(k,127) = rxt(k,213)*y(k,36)
         mat(k,677) = rxt(k,157)*y(k,23) + rxt(k,213)*y(k,34) + rxt(k,111)*y(k,56) &
                      + rxt(k,94)*y(k,73) + rxt(k,120)*y(k,77)
         mat(k,151) = rxt(k,211)*y(k,73)
         mat(k,308) = rxt(k,188)*y(k,73)
         mat(k,226) = rxt(k,143)*y(k,77)
         mat(k,414) = rxt(k,111)*y(k,36) + rxt(k,123)*y(k,77)
         mat(k,55) = rxt(k,250)*y(k,77)
         mat(k,114) = rxt(k,255)*y(k,77)
         mat(k,267) = rxt(k,260)*y(k,77)
         mat(k,523) = rxt(k,94)*y(k,36) + rxt(k,211)*y(k,40) + rxt(k,188)*y(k,44)
         mat(k,632) = rxt(k,237)*y(k,15) + rxt(k,246)*y(k,29) + rxt(k,120)*y(k,36) &
                      + rxt(k,143)*y(k,50) + rxt(k,123)*y(k,56) + rxt(k,250)*y(k,60) &
                      + rxt(k,255)*y(k,61) + rxt(k,260)*y(k,63)
         mat(k,692) = -(rxt(k,94)*y(k,73) + rxt(k,111)*y(k,56) + rxt(k,120)*y(k,77) &
                      + rxt(k,157)*y(k,23) + rxt(k,213)*y(k,34))
         mat(k,539) = -rxt(k,94)*y(k,36)
         mat(k,430) = -rxt(k,111)*y(k,36)
         mat(k,648) = -rxt(k,120)*y(k,36)
         mat(k,483) = -rxt(k,157)*y(k,36)
         mat(k,131) = -rxt(k,213)*y(k,36)
         mat(k,286) = rxt(k,113)*y(k,68)
         mat(k,387) = rxt(k,113)*y(k,35)
         mat(k,96) = -(rxt(k,112)*y(k,56) + rxt(k,121)*y(k,77) + rxt(k,158)*y(k,23))
         mat(k,397) = -rxt(k,112)*y(k,38)
         mat(k,622) = -rxt(k,121)*y(k,38)
         mat(k,463) = -rxt(k,158)*y(k,38)
         mat(k,367) = 2.000_r8*rxt(k,127)*y(k,68)
         mat(k,622) = mat(k,622) + 2.000_r8*rxt(k,126)*y(k,77)
         mat(k,28) = rxt(k,263)*y(k,81)
         mat(k,721) = rxt(k,263)*y(k,65)
         mat(k,150) = -(rxt(k,204)*y(k,56) + rxt(k,205)*y(k,77) + (rxt(k,210) &
                      + rxt(k,211)) * y(k,73))
         mat(k,403) = -rxt(k,204)*y(k,40)
         mat(k,626) = -rxt(k,205)*y(k,40)
         mat(k,521) = -(rxt(k,210) + rxt(k,211)) * y(k,40)
         mat(k,289) = rxt(k,191)*y(k,15) + rxt(k,192)*y(k,68)
         mat(k,543) = rxt(k,191)*y(k,1)
         mat(k,370) = rxt(k,192)*y(k,1)
         mat(k,310) = -(rxt(k,175)*y(k,56) + rxt(k,176)*y(k,77) + (rxt(k,187) &
                      + rxt(k,188)) * y(k,73) + (rxt(k,269) + rxt(k,275) + rxt(k,280) &
                      ) * y(k,49) + (rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,27) &
                      + (rxt(k,276) + rxt(k,281)) * y(k,48))
         mat(k,416) = -rxt(k,175)*y(k,44)
         mat(k,634) = -rxt(k,176)*y(k,44)
         mat(k,525) = -(rxt(k,187) + rxt(k,188)) * y(k,44)
         mat(k,144) = -(rxt(k,269) + rxt(k,275) + rxt(k,280)) * y(k,44)
         mat(k,245) = -(rxt(k,274) + rxt(k,279) + rxt(k,284)) * y(k,44)
         mat(k,136) = -(rxt(k,276) + rxt(k,281)) * y(k,44)
         mat(k,547) = rxt(k,155)*y(k,23)
         mat(k,469) = rxt(k,155)*y(k,15) + rxt(k,157)*y(k,36) + rxt(k,158)*y(k,38) &
                      + rxt(k,177)*y(k,49) + rxt(k,159)*y(k,68)
         mat(k,444) = rxt(k,174)*y(k,77)
         mat(k,678) = rxt(k,157)*y(k,23)
         mat(k,97) = rxt(k,158)*y(k,23)
         mat(k,144) = mat(k,144) + rxt(k,177)*y(k,23)
         mat(k,373) = rxt(k,159)*y(k,23)
         mat(k,634) = mat(k,634) + rxt(k,174)*y(k,26)
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
         mat(k,124) = rxt(k,213)*y(k,36) + rxt(k,215)*y(k,46) + rxt(k,214)*y(k,81)
         mat(k,674) = rxt(k,213)*y(k,34)
         mat(k,253) = rxt(k,215)*y(k,34)
         mat(k,722) = rxt(k,214)*y(k,34)
         mat(k,255) = -(rxt(k,152)*y(k,77) + rxt(k,215)*y(k,34))
         mat(k,630) = -rxt(k,152)*y(k,46)
         mat(k,126) = -rxt(k,215)*y(k,46)
         mat(k,544) = rxt(k,235)*y(k,55)
         mat(k,244) = (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,44)
         mat(k,35) = rxt(k,247)*y(k,55)
         mat(k,307) = (rxt(k,274)+rxt(k,279)+rxt(k,284))*y(k,27)
         mat(k,593) = rxt(k,151)*y(k,77)
         mat(k,345) = rxt(k,235)*y(k,15) + rxt(k,247)*y(k,33)
         mat(k,630) = mat(k,630) + rxt(k,151)*y(k,54)
         mat(k,64) = -(rxt(k,128)*y(k,77))
         mat(k,620) = -rxt(k,128)*y(k,47)
         mat(k,587) = rxt(k,149)*y(k,68)
         mat(k,366) = rxt(k,149)*y(k,54)
         mat(k,134) = -(rxt(k,206)*y(k,56) + (rxt(k,276) + rxt(k,281)) * y(k,44))
         mat(k,401) = -rxt(k,206)*y(k,48)
         mat(k,304) = -(rxt(k,276) + rxt(k,281)) * y(k,48)
         mat(k,698) = rxt(k,198)*y(k,68)
         mat(k,368) = rxt(k,198)*y(k,3)
         mat(k,143) = -(rxt(k,177)*y(k,23) + rxt(k,178)*y(k,56) + rxt(k,179)*y(k,77) &
                      + (rxt(k,269) + rxt(k,275) + rxt(k,280)) * y(k,44))
         mat(k,464) = -rxt(k,177)*y(k,49)
         mat(k,402) = -rxt(k,178)*y(k,49)
         mat(k,625) = -rxt(k,179)*y(k,49)
         mat(k,305) = -(rxt(k,269) + rxt(k,275) + rxt(k,280)) * y(k,49)
         mat(k,437) = rxt(k,166)*y(k,68)
         mat(k,242) = rxt(k,171)*y(k,77)
         mat(k,369) = rxt(k,166)*y(k,26)
         mat(k,625) = mat(k,625) + rxt(k,171)*y(k,27)
         mat(k,225) = -(rxt(k,131)*y(k,53) + (rxt(k,132) + rxt(k,133) + rxt(k,134) &
                      ) * y(k,54) + rxt(k,135)*y(k,57) + rxt(k,143)*y(k,77) + rxt(k,297) &
                      *y(k,76))
         mat(k,655) = -rxt(k,131)*y(k,50)
         mat(k,590) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,50)
         mat(k,498) = -rxt(k,135)*y(k,50)
         mat(k,627) = -rxt(k,143)*y(k,50)
         mat(k,184) = -rxt(k,297)*y(k,50)
         mat(k,409) = rxt(k,129)*y(k,69) + rxt(k,294)*y(k,72)
         mat(k,498) = mat(k,498) + rxt(k,295)*y(k,72)
         mat(k,198) = 1.100_r8*rxt(k,290)*y(k,70) + .200_r8*rxt(k,288)*y(k,71)
         mat(k,211) = rxt(k,129)*y(k,56)
         mat(k,109) = 1.100_r8*rxt(k,290)*y(k,67)
         mat(k,122) = .200_r8*rxt(k,288)*y(k,67)
         mat(k,164) = rxt(k,294)*y(k,56) + rxt(k,295)*y(k,57)
         mat(k,586) = rxt(k,150)*y(k,55)
         mat(k,342) = rxt(k,150)*y(k,54)
         mat(k,670) = -(rxt(k,131)*y(k,50) + rxt(k,140)*y(k,55) + rxt(k,144)*y(k,68) &
                      + rxt(k,145)*y(k,58) + rxt(k,146)*y(k,56) + rxt(k,167)*y(k,26) &
                      + rxt(k,199)*y(k,3) + rxt(k,239)*y(k,19) + rxt(k,299)*y(k,76))
         mat(k,233) = -rxt(k,131)*y(k,53)
         mat(k,361) = -rxt(k,140)*y(k,53)
         mat(k,386) = -rxt(k,144)*y(k,53)
         mat(k,582) = -rxt(k,145)*y(k,53)
         mat(k,429) = -rxt(k,146)*y(k,53)
         mat(k,457) = -rxt(k,167)*y(k,53)
         mat(k,717) = -rxt(k,199)*y(k,53)
         mat(k,338) = -rxt(k,239)*y(k,53)
         mat(k,188) = -rxt(k,299)*y(k,53)
         mat(k,233) = mat(k,233) + 2.000_r8*rxt(k,133)*y(k,54) + rxt(k,135)*y(k,57) &
                      + rxt(k,143)*y(k,77)
         mat(k,610) = 2.000_r8*rxt(k,133)*y(k,50) + rxt(k,136)*y(k,56) + rxt(k,256) &
                      *y(k,63)
         mat(k,429) = mat(k,429) + rxt(k,136)*y(k,54)
         mat(k,513) = rxt(k,135)*y(k,50) + rxt(k,130)*y(k,69)
         mat(k,277) = rxt(k,256)*y(k,54)
         mat(k,218) = rxt(k,130)*y(k,57)
         mat(k,647) = rxt(k,143)*y(k,50)
         mat(k,608) = -((rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,50) + (rxt(k,136) &
                      + rxt(k,138)) * y(k,56) + rxt(k,137)*y(k,58) + rxt(k,149) &
                      *y(k,68) + rxt(k,150)*y(k,55) + rxt(k,151)*y(k,77) + rxt(k,169) &
                      *y(k,26) + rxt(k,200)*y(k,3) + rxt(k,256)*y(k,63))
         mat(k,231) = -(rxt(k,132) + rxt(k,133) + rxt(k,134)) * y(k,54)
         mat(k,427) = -(rxt(k,136) + rxt(k,138)) * y(k,54)
         mat(k,580) = -rxt(k,137)*y(k,54)
         mat(k,384) = -rxt(k,149)*y(k,54)
         mat(k,359) = -rxt(k,150)*y(k,54)
         mat(k,645) = -rxt(k,151)*y(k,54)
         mat(k,455) = -rxt(k,169)*y(k,54)
         mat(k,715) = -rxt(k,200)*y(k,54)
         mat(k,275) = -rxt(k,256)*y(k,54)
         mat(k,715) = mat(k,715) + rxt(k,199)*y(k,53)
         mat(k,336) = rxt(k,239)*y(k,53)
         mat(k,455) = mat(k,455) + rxt(k,167)*y(k,53)
         mat(k,68) = rxt(k,128)*y(k,77)
         mat(k,668) = rxt(k,199)*y(k,3) + rxt(k,239)*y(k,19) + rxt(k,167)*y(k,26) &
                      + 2.000_r8*rxt(k,140)*y(k,55) + rxt(k,146)*y(k,56) + rxt(k,145) &
                      *y(k,58) + rxt(k,144)*y(k,68)
         mat(k,359) = mat(k,359) + 2.000_r8*rxt(k,140)*y(k,53) + rxt(k,141)*y(k,56) &
                      + rxt(k,139)*y(k,68) + rxt(k,142)*y(k,77)
         mat(k,427) = mat(k,427) + rxt(k,146)*y(k,53) + rxt(k,141)*y(k,55)
         mat(k,580) = mat(k,580) + rxt(k,145)*y(k,53)
         mat(k,384) = mat(k,384) + rxt(k,144)*y(k,53) + rxt(k,139)*y(k,55)
         mat(k,645) = mat(k,645) + rxt(k,128)*y(k,47) + rxt(k,142)*y(k,55)
         mat(k,350) = -(rxt(k,139)*y(k,68) + rxt(k,140)*y(k,53) + rxt(k,141)*y(k,56) &
                      + rxt(k,142)*y(k,77) + rxt(k,150)*y(k,54) + rxt(k,235)*y(k,15) &
                      + rxt(k,247)*y(k,33))
         mat(k,375) = -rxt(k,139)*y(k,55)
         mat(k,659) = -rxt(k,140)*y(k,55)
         mat(k,418) = -rxt(k,141)*y(k,55)
         mat(k,636) = -rxt(k,142)*y(k,55)
         mat(k,599) = -rxt(k,150)*y(k,55)
         mat(k,549) = -rxt(k,235)*y(k,55)
         mat(k,36) = -rxt(k,247)*y(k,55)
         mat(k,88) = rxt(k,201)*y(k,56)
         mat(k,471) = rxt(k,168)*y(k,27)
         mat(k,246) = rxt(k,168)*y(k,23) + rxt(k,170)*y(k,56) + rxt(k,171)*y(k,77)
         mat(k,129) = rxt(k,215)*y(k,46)
         mat(k,258) = rxt(k,215)*y(k,34) + rxt(k,152)*y(k,77)
         mat(k,599) = mat(k,599) + rxt(k,138)*y(k,56) + rxt(k,137)*y(k,58)
         mat(k,418) = mat(k,418) + rxt(k,201)*y(k,4) + rxt(k,170)*y(k,27) + rxt(k,138) &
                      *y(k,54)
         mat(k,571) = rxt(k,137)*y(k,54)
         mat(k,636) = mat(k,636) + rxt(k,171)*y(k,27) + rxt(k,152)*y(k,46)
         mat(k,420) = -(rxt(k,108)*y(k,58) + 4._r8*rxt(k,109)*y(k,56) + rxt(k,110) &
                      *y(k,57) + rxt(k,111)*y(k,36) + rxt(k,112)*y(k,38) + rxt(k,117) &
                      *y(k,68) + rxt(k,123)*y(k,77) + (rxt(k,136) + rxt(k,138) &
                      ) * y(k,54) + rxt(k,141)*y(k,55) + rxt(k,146)*y(k,53) + rxt(k,170) &
                      *y(k,27) + rxt(k,172)*y(k,26) + rxt(k,175)*y(k,44) + rxt(k,178) &
                      *y(k,49) + rxt(k,201)*y(k,4) + rxt(k,202)*y(k,3) + rxt(k,204) &
                      *y(k,40) + rxt(k,206)*y(k,48) + rxt(k,236)*y(k,15) + rxt(k,249) &
                      *y(k,60) + (rxt(k,292) + rxt(k,293)) * y(k,70) + rxt(k,294) &
                      *y(k,72))
         mat(k,573) = -rxt(k,108)*y(k,56)
         mat(k,504) = -rxt(k,110)*y(k,56)
         mat(k,682) = -rxt(k,111)*y(k,56)
         mat(k,99) = -rxt(k,112)*y(k,56)
         mat(k,377) = -rxt(k,117)*y(k,56)
         mat(k,638) = -rxt(k,123)*y(k,56)
         mat(k,601) = -(rxt(k,136) + rxt(k,138)) * y(k,56)
         mat(k,352) = -rxt(k,141)*y(k,56)
         mat(k,661) = -rxt(k,146)*y(k,56)
         mat(k,247) = -rxt(k,170)*y(k,56)
         mat(k,448) = -rxt(k,172)*y(k,56)
         mat(k,314) = -rxt(k,175)*y(k,56)
         mat(k,145) = -rxt(k,178)*y(k,56)
         mat(k,89) = -rxt(k,201)*y(k,56)
         mat(k,708) = -rxt(k,202)*y(k,56)
         mat(k,153) = -rxt(k,204)*y(k,56)
         mat(k,137) = -rxt(k,206)*y(k,56)
         mat(k,551) = -rxt(k,236)*y(k,56)
         mat(k,56) = -rxt(k,249)*y(k,56)
         mat(k,110) = -(rxt(k,292) + rxt(k,293)) * y(k,56)
         mat(k,165) = -rxt(k,294)*y(k,56)
         mat(k,282) = rxt(k,115)*y(k,68)
         mat(k,228) = rxt(k,131)*y(k,53) + rxt(k,132)*y(k,54) + rxt(k,135)*y(k,57) &
                      + rxt(k,297)*y(k,76)
         mat(k,661) = mat(k,661) + rxt(k,131)*y(k,50)
         mat(k,601) = mat(k,601) + rxt(k,132)*y(k,50)
         mat(k,504) = mat(k,504) + rxt(k,135)*y(k,50) + rxt(k,251)*y(k,61) &
                      + rxt(k,257)*y(k,63) + rxt(k,296)*y(k,72) + (rxt(k,97)+rxt(k,98)) &
                      *y(k,73) + rxt(k,303)*y(k,78) + rxt(k,307)*y(k,79)
         mat(k,115) = rxt(k,251)*y(k,57)
         mat(k,270) = rxt(k,257)*y(k,57)
         mat(k,201) = rxt(k,288)*y(k,71) + 1.150_r8*rxt(k,289)*y(k,76)
         mat(k,377) = mat(k,377) + rxt(k,115)*y(k,35)
         mat(k,214) = rxt(k,302)*y(k,78)
         mat(k,123) = rxt(k,288)*y(k,67)
         mat(k,165) = mat(k,165) + rxt(k,296)*y(k,57)
         mat(k,529) = (rxt(k,97)+rxt(k,98))*y(k,57)
         mat(k,185) = rxt(k,297)*y(k,50) + 1.150_r8*rxt(k,289)*y(k,67)
         mat(k,638) = mat(k,638) + 2.000_r8*rxt(k,125)*y(k,77)
         mat(k,177) = rxt(k,303)*y(k,57) + rxt(k,302)*y(k,69)
         mat(k,82) = rxt(k,307)*y(k,57)
         mat(k,507) = -(rxt(k,97)*y(k,73) + rxt(k,102)*y(k,74) + rxt(k,110)*y(k,56) &
                      + rxt(k,116)*y(k,35) + rxt(k,130)*y(k,69) + rxt(k,135)*y(k,50) &
                      + rxt(k,251)*y(k,61) + rxt(k,257)*y(k,63) + rxt(k,291)*y(k,70) &
                      + (rxt(k,295) + rxt(k,296)) * y(k,72) + rxt(k,303)*y(k,78) &
                      + rxt(k,307)*y(k,79))
         mat(k,532) = -rxt(k,97)*y(k,57)
         mat(k,12) = -rxt(k,102)*y(k,57)
         mat(k,423) = -rxt(k,110)*y(k,57)
         mat(k,283) = -rxt(k,116)*y(k,57)
         mat(k,215) = -rxt(k,130)*y(k,57)
         mat(k,229) = -rxt(k,135)*y(k,57)
         mat(k,116) = -rxt(k,251)*y(k,57)
         mat(k,273) = -rxt(k,257)*y(k,57)
         mat(k,111) = -rxt(k,291)*y(k,57)
         mat(k,166) = -(rxt(k,295) + rxt(k,296)) * y(k,57)
         mat(k,178) = -rxt(k,303)*y(k,57)
         mat(k,83) = -rxt(k,307)*y(k,57)
         mat(k,294) = rxt(k,193)*y(k,58) + rxt(k,192)*y(k,68)
         mat(k,711) = 2.000_r8*rxt(k,194)*y(k,3) + (rxt(k,196)+rxt(k,197))*y(k,26) &
                      + rxt(k,202)*y(k,56) + rxt(k,198)*y(k,68)
         mat(k,333) = rxt(k,238)*y(k,68)
         mat(k,476) = rxt(k,161)*y(k,58) + rxt(k,159)*y(k,68)
         mat(k,451) = (rxt(k,196)+rxt(k,197))*y(k,3) + (2.000_r8*rxt(k,163) &
                       +2.000_r8*rxt(k,164))*y(k,26) + rxt(k,172)*y(k,56) + rxt(k,166) &
                      *y(k,68) + rxt(k,174)*y(k,77)
         mat(k,283) = mat(k,283) + rxt(k,119)*y(k,58) + rxt(k,113)*y(k,68)
         mat(k,67) = rxt(k,128)*y(k,77)
         mat(k,229) = mat(k,229) + rxt(k,134)*y(k,54)
         mat(k,664) = rxt(k,145)*y(k,58) + rxt(k,299)*y(k,76)
         mat(k,604) = rxt(k,134)*y(k,50) + rxt(k,136)*y(k,56) + rxt(k,137)*y(k,58)
         mat(k,355) = rxt(k,141)*y(k,56) + rxt(k,139)*y(k,68)
         mat(k,423) = mat(k,423) + rxt(k,202)*y(k,3) + rxt(k,172)*y(k,26) + rxt(k,136) &
                      *y(k,54) + rxt(k,141)*y(k,55) + 2.000_r8*rxt(k,109)*y(k,56) &
                      + 2.000_r8*rxt(k,108)*y(k,58) + rxt(k,117)*y(k,68) + rxt(k,101) &
                      *y(k,74) + rxt(k,123)*y(k,77)
         mat(k,507) = mat(k,507) + 2.000_r8*rxt(k,102)*y(k,74)
         mat(k,576) = rxt(k,193)*y(k,1) + rxt(k,161)*y(k,23) + rxt(k,119)*y(k,35) &
                      + rxt(k,145)*y(k,53) + rxt(k,137)*y(k,54) + 2.000_r8*rxt(k,108) &
                      *y(k,56) + rxt(k,252)*y(k,61) + rxt(k,258)*y(k,63) &
                      + 2.000_r8*rxt(k,118)*y(k,68) + 2.000_r8*rxt(k,99)*y(k,73) &
                      + rxt(k,124)*y(k,77)
         mat(k,116) = mat(k,116) + rxt(k,252)*y(k,58)
         mat(k,273) = mat(k,273) + rxt(k,258)*y(k,58)
         mat(k,380) = rxt(k,192)*y(k,1) + rxt(k,198)*y(k,3) + rxt(k,238)*y(k,19) &
                      + rxt(k,159)*y(k,23) + rxt(k,166)*y(k,26) + rxt(k,113)*y(k,35) &
                      + rxt(k,139)*y(k,55) + rxt(k,117)*y(k,56) + 2.000_r8*rxt(k,118) &
                      *y(k,58) + 2.000_r8*rxt(k,127)*y(k,68) + rxt(k,122)*y(k,77)
         mat(k,532) = mat(k,532) + 2.000_r8*rxt(k,99)*y(k,58)
         mat(k,12) = mat(k,12) + rxt(k,101)*y(k,56) + 2.000_r8*rxt(k,102)*y(k,57)
         mat(k,186) = rxt(k,299)*y(k,53)
         mat(k,641) = rxt(k,174)*y(k,26) + rxt(k,128)*y(k,47) + rxt(k,123)*y(k,56) &
                      + rxt(k,124)*y(k,58) + rxt(k,122)*y(k,68)
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
         mat(k,579) = -(rxt(k,99)*y(k,73) + rxt(k,108)*y(k,56) + rxt(k,118)*y(k,68) &
                      + rxt(k,119)*y(k,35) + rxt(k,124)*y(k,77) + rxt(k,137)*y(k,54) &
                      + rxt(k,145)*y(k,53) + rxt(k,161)*y(k,23) + rxt(k,193)*y(k,1) &
                      + rxt(k,252)*y(k,61) + rxt(k,258)*y(k,63))
         mat(k,535) = -rxt(k,99)*y(k,58)
         mat(k,426) = -rxt(k,108)*y(k,58)
         mat(k,383) = -rxt(k,118)*y(k,58)
         mat(k,284) = -rxt(k,119)*y(k,58)
         mat(k,644) = -rxt(k,124)*y(k,58)
         mat(k,607) = -rxt(k,137)*y(k,58)
         mat(k,667) = -rxt(k,145)*y(k,58)
         mat(k,479) = -rxt(k,161)*y(k,58)
         mat(k,297) = -rxt(k,193)*y(k,58)
         mat(k,117) = -rxt(k,252)*y(k,58)
         mat(k,274) = -rxt(k,258)*y(k,58)
         mat(k,426) = mat(k,426) + rxt(k,110)*y(k,57)
         mat(k,510) = rxt(k,110)*y(k,56)
         mat(k,45) = -(rxt(k,259)*y(k,63))
         mat(k,263) = -rxt(k,259)*y(k,59)
         mat(k,696) = rxt(k,195)*y(k,26)
         mat(k,436) = rxt(k,195)*y(k,3) + 2.000_r8*rxt(k,165)*y(k,26)
         mat(k,50) = -(rxt(k,249)*y(k,56) + rxt(k,250)*y(k,77))
         mat(k,392) = -rxt(k,249)*y(k,60)
         mat(k,618) = -rxt(k,250)*y(k,60)
         mat(k,112) = -(rxt(k,251)*y(k,57) + rxt(k,252)*y(k,58) + rxt(k,255)*y(k,77))
         mat(k,491) = -rxt(k,251)*y(k,61)
         mat(k,566) = -rxt(k,252)*y(k,61)
         mat(k,623) = -rxt(k,255)*y(k,61)
         mat(k,266) = -(rxt(k,253)*y(k,3) + rxt(k,254)*y(k,26) + rxt(k,256)*y(k,54) &
                      + rxt(k,257)*y(k,57) + rxt(k,258)*y(k,58) + rxt(k,259)*y(k,59) &
                      + rxt(k,260)*y(k,77))
         mat(k,701) = -rxt(k,253)*y(k,63)
         mat(k,441) = -rxt(k,254)*y(k,63)
         mat(k,594) = -rxt(k,256)*y(k,63)
         mat(k,500) = -rxt(k,257)*y(k,63)
         mat(k,568) = -rxt(k,258)*y(k,63)
         mat(k,47) = -rxt(k,259)*y(k,63)
         mat(k,631) = -rxt(k,260)*y(k,63)
         mat(k,413) = rxt(k,249)*y(k,60)
         mat(k,500) = mat(k,500) + rxt(k,251)*y(k,61)
         mat(k,568) = mat(k,568) + rxt(k,252)*y(k,61)
         mat(k,54) = rxt(k,249)*y(k,56)
         mat(k,113) = rxt(k,251)*y(k,57) + rxt(k,252)*y(k,58) + rxt(k,255)*y(k,77)
         mat(k,631) = mat(k,631) + rxt(k,255)*y(k,61)
         mat(k,235) = -(rxt(k,262)*y(k,77))
         mat(k,628) = -rxt(k,262)*y(k,64)
         mat(k,699) = rxt(k,253)*y(k,63)
         mat(k,438) = rxt(k,254)*y(k,63)
         mat(k,34) = rxt(k,247)*y(k,55) + (rxt(k,248)+.500_r8*rxt(k,261))*y(k,77)
         mat(k,591) = rxt(k,256)*y(k,63)
         mat(k,344) = rxt(k,247)*y(k,33)
         mat(k,499) = rxt(k,257)*y(k,63)
         mat(k,567) = rxt(k,258)*y(k,63)
         mat(k,46) = rxt(k,259)*y(k,63)
         mat(k,53) = rxt(k,250)*y(k,77)
         mat(k,265) = rxt(k,253)*y(k,3) + rxt(k,254)*y(k,26) + rxt(k,256)*y(k,54) &
                      + rxt(k,257)*y(k,57) + rxt(k,258)*y(k,58) + rxt(k,259)*y(k,59) &
                      + rxt(k,260)*y(k,77)
         mat(k,628) = mat(k,628) + (rxt(k,248)+.500_r8*rxt(k,261))*y(k,33) &
                      + rxt(k,250)*y(k,60) + rxt(k,260)*y(k,63)
         mat(k,29) = -(rxt(k,263)*y(k,81))
         mat(k,723) = -rxt(k,263)*y(k,65)
         mat(k,234) = rxt(k,262)*y(k,77)
         mat(k,616) = rxt(k,262)*y(k,64)
         mat(k,196) = -(rxt(k,288)*y(k,71) + rxt(k,289)*y(k,76) + rxt(k,290)*y(k,70))
         mat(k,120) = -rxt(k,288)*y(k,67)
         mat(k,182) = -rxt(k,289)*y(k,67)
         mat(k,107) = -rxt(k,290)*y(k,67)
         mat(k,376) = -((rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,35) + rxt(k,117) &
                      *y(k,56) + rxt(k,118)*y(k,58) + rxt(k,122)*y(k,77) &
                      + 4._r8*rxt(k,127)*y(k,68) + rxt(k,139)*y(k,55) + rxt(k,144) &
                      *y(k,53) + rxt(k,149)*y(k,54) + (rxt(k,159) + rxt(k,160) &
                      ) * y(k,23) + rxt(k,166)*y(k,26) + rxt(k,192)*y(k,1) + rxt(k,198) &
                      *y(k,3) + rxt(k,238)*y(k,19))
         mat(k,281) = -(rxt(k,113) + rxt(k,114) + rxt(k,115)) * y(k,68)
         mat(k,419) = -rxt(k,117)*y(k,68)
         mat(k,572) = -rxt(k,118)*y(k,68)
         mat(k,637) = -rxt(k,122)*y(k,68)
         mat(k,351) = -rxt(k,139)*y(k,68)
         mat(k,660) = -rxt(k,144)*y(k,68)
         mat(k,600) = -rxt(k,149)*y(k,68)
         mat(k,472) = -(rxt(k,159) + rxt(k,160)) * y(k,68)
         mat(k,447) = -rxt(k,166)*y(k,68)
         mat(k,292) = -rxt(k,192)*y(k,68)
         mat(k,707) = -rxt(k,198)*y(k,68)
         mat(k,329) = -rxt(k,238)*y(k,68)
         mat(k,292) = mat(k,292) + rxt(k,191)*y(k,15)
         mat(k,707) = mat(k,707) + rxt(k,203)*y(k,77)
         mat(k,550) = rxt(k,191)*y(k,1) + rxt(k,155)*y(k,23) + rxt(k,235)*y(k,55) &
                      + rxt(k,236)*y(k,56)
         mat(k,329) = mat(k,329) + rxt(k,162)*y(k,26) + rxt(k,239)*y(k,53)
         mat(k,472) = mat(k,472) + rxt(k,155)*y(k,15) + rxt(k,158)*y(k,38)
         mat(k,447) = mat(k,447) + rxt(k,162)*y(k,19) + rxt(k,173)*y(k,77)
         mat(k,94) = rxt(k,242)*y(k,77)
         mat(k,37) = .500_r8*rxt(k,261)*y(k,77)
         mat(k,281) = mat(k,281) + rxt(k,116)*y(k,57)
         mat(k,98) = rxt(k,158)*y(k,23) + rxt(k,112)*y(k,56) + rxt(k,121)*y(k,77)
         mat(k,660) = mat(k,660) + rxt(k,239)*y(k,19)
         mat(k,351) = mat(k,351) + rxt(k,235)*y(k,15) + rxt(k,142)*y(k,77)
         mat(k,419) = mat(k,419) + rxt(k,236)*y(k,15) + rxt(k,112)*y(k,38)
         mat(k,503) = rxt(k,116)*y(k,35)
         mat(k,572) = mat(k,572) + rxt(k,124)*y(k,77)
         mat(k,237) = rxt(k,262)*y(k,77)
         mat(k,637) = mat(k,637) + rxt(k,203)*y(k,3) + rxt(k,173)*y(k,26) + rxt(k,242) &
                      *y(k,29) + .500_r8*rxt(k,261)*y(k,33) + rxt(k,121)*y(k,38) &
                      + rxt(k,142)*y(k,55) + rxt(k,124)*y(k,58) + rxt(k,262)*y(k,64)
         mat(k,210) = -(rxt(k,129)*y(k,56) + rxt(k,130)*y(k,57) + rxt(k,302)*y(k,78))
         mat(k,408) = -rxt(k,129)*y(k,69)
         mat(k,497) = -rxt(k,130)*y(k,69)
         mat(k,173) = -rxt(k,302)*y(k,69)
         mat(k,408) = mat(k,408) + rxt(k,292)*y(k,70)
         mat(k,197) = .900_r8*rxt(k,290)*y(k,70) + .800_r8*rxt(k,288)*y(k,71)
         mat(k,108) = rxt(k,292)*y(k,56) + .900_r8*rxt(k,290)*y(k,67)
         mat(k,121) = .800_r8*rxt(k,288)*y(k,67)
         mat(k,103) = -(rxt(k,290)*y(k,67) + rxt(k,291)*y(k,57) + (rxt(k,292) &
                      + rxt(k,293)) * y(k,56))
         mat(k,191) = -rxt(k,290)*y(k,70)
         mat(k,490) = -rxt(k,291)*y(k,70)
         mat(k,398) = -(rxt(k,292) + rxt(k,293)) * y(k,70)
         mat(k,119) = -(rxt(k,288)*y(k,67))
         mat(k,192) = -rxt(k,288)*y(k,71)
         mat(k,219) = rxt(k,297)*y(k,76)
         mat(k,651) = rxt(k,299)*y(k,76)
         mat(k,400) = rxt(k,292)*y(k,70)
         mat(k,492) = rxt(k,296)*y(k,72)
         mat(k,104) = rxt(k,292)*y(k,56)
         mat(k,158) = rxt(k,296)*y(k,57)
         mat(k,180) = rxt(k,297)*y(k,50) + rxt(k,299)*y(k,53)
         mat(k,159) = -(rxt(k,294)*y(k,56) + (rxt(k,295) + rxt(k,296)) * y(k,57))
         mat(k,404) = -rxt(k,294)*y(k,72)
         mat(k,493) = -(rxt(k,295) + rxt(k,296)) * y(k,72)
         mat(k,206) = rxt(k,302)*y(k,78)
         mat(k,169) = rxt(k,302)*y(k,69)
         mat(k,533) = -(rxt(k,94)*y(k,36) + rxt(k,95)*y(k,81) + (rxt(k,97) + rxt(k,98) &
                      ) * y(k,57) + rxt(k,99)*y(k,58) + (rxt(k,187) + rxt(k,188) &
                      ) * y(k,44) + (rxt(k,210) + rxt(k,211)) * y(k,40) + rxt(k,216) &
                      *y(k,31) + rxt(k,217)*y(k,32))
         mat(k,686) = -rxt(k,94)*y(k,73)
         mat(k,737) = -rxt(k,95)*y(k,73)
         mat(k,508) = -(rxt(k,97) + rxt(k,98)) * y(k,73)
         mat(k,577) = -rxt(k,99)*y(k,73)
         mat(k,318) = -(rxt(k,187) + rxt(k,188)) * y(k,73)
         mat(k,154) = -(rxt(k,210) + rxt(k,211)) * y(k,73)
         mat(k,7) = -rxt(k,216)*y(k,73)
         mat(k,21) = -rxt(k,217)*y(k,73)
         mat(k,508) = mat(k,508) + rxt(k,130)*y(k,69)
         mat(k,203) = .850_r8*rxt(k,289)*y(k,76)
         mat(k,216) = rxt(k,130)*y(k,57)
         mat(k,187) = .850_r8*rxt(k,289)*y(k,67)
         mat(k,11) = -(rxt(k,101)*y(k,56) + rxt(k,102)*y(k,57))
         mat(k,390) = -rxt(k,101)*y(k,74)
         mat(k,486) = -rxt(k,102)*y(k,74)
         mat(k,390) = mat(k,390) + rxt(k,105)*y(k,75)
         mat(k,486) = mat(k,486) + rxt(k,106)*y(k,75)
         mat(k,564) = rxt(k,107)*y(k,75)
         mat(k,13) = rxt(k,105)*y(k,56) + rxt(k,106)*y(k,57) + rxt(k,107)*y(k,58)
         mat(k,14) = -(rxt(k,105)*y(k,56) + rxt(k,106)*y(k,57) + rxt(k,107)*y(k,58))
         mat(k,391) = -rxt(k,105)*y(k,75)
         mat(k,487) = -rxt(k,106)*y(k,75)
         mat(k,565) = -rxt(k,107)*y(k,75)
         mat(k,487) = mat(k,487) + rxt(k,97)*y(k,73)
         mat(k,518) = rxt(k,97)*y(k,57)
         mat(k,181) = -(rxt(k,289)*y(k,67) + rxt(k,297)*y(k,50) + rxt(k,299)*y(k,53))
         mat(k,195) = -rxt(k,289)*y(k,76)
         mat(k,222) = -rxt(k,297)*y(k,76)
         mat(k,652) = -rxt(k,299)*y(k,76)
         mat(k,495) = rxt(k,291)*y(k,70) + rxt(k,295)*y(k,72) + rxt(k,303)*y(k,78) &
                      + rxt(k,307)*y(k,79)
         mat(k,106) = rxt(k,291)*y(k,57)
         mat(k,161) = rxt(k,295)*y(k,57)
         mat(k,171) = rxt(k,303)*y(k,57)
         mat(k,81) = rxt(k,307)*y(k,57)
         mat(k,646) = -(rxt(k,120)*y(k,36) + rxt(k,121)*y(k,38) + rxt(k,122)*y(k,68) &
                      + rxt(k,123)*y(k,56) + rxt(k,124)*y(k,58) + (4._r8*rxt(k,125) &
                      + 4._r8*rxt(k,126)) * y(k,77) + rxt(k,128)*y(k,47) + rxt(k,142) &
                      *y(k,55) + rxt(k,143)*y(k,50) + rxt(k,151)*y(k,54) + rxt(k,152) &
                      *y(k,46) + rxt(k,171)*y(k,27) + (rxt(k,173) + rxt(k,174) &
                      ) * y(k,26) + rxt(k,176)*y(k,44) + rxt(k,179)*y(k,49) + rxt(k,203) &
                      *y(k,3) + rxt(k,205)*y(k,40) + rxt(k,237)*y(k,15) + rxt(k,240) &
                      *y(k,20) + (rxt(k,242) + rxt(k,246)) * y(k,29) + (rxt(k,248) &
                      + rxt(k,261)) * y(k,33) + rxt(k,250)*y(k,60) + rxt(k,255) &
                      *y(k,61) + rxt(k,260)*y(k,63) + rxt(k,262)*y(k,64))
         mat(k,690) = -rxt(k,120)*y(k,77)
         mat(k,101) = -rxt(k,121)*y(k,77)
         mat(k,385) = -rxt(k,122)*y(k,77)
         mat(k,428) = -rxt(k,123)*y(k,77)
         mat(k,581) = -rxt(k,124)*y(k,77)
         mat(k,69) = -rxt(k,128)*y(k,77)
         mat(k,360) = -rxt(k,142)*y(k,77)
         mat(k,232) = -rxt(k,143)*y(k,77)
         mat(k,609) = -rxt(k,151)*y(k,77)
         mat(k,260) = -rxt(k,152)*y(k,77)
         mat(k,251) = -rxt(k,171)*y(k,77)
         mat(k,456) = -(rxt(k,173) + rxt(k,174)) * y(k,77)
         mat(k,322) = -rxt(k,176)*y(k,77)
         mat(k,148) = -rxt(k,179)*y(k,77)
         mat(k,716) = -rxt(k,203)*y(k,77)
         mat(k,155) = -rxt(k,205)*y(k,77)
         mat(k,559) = -rxt(k,237)*y(k,77)
         mat(k,62) = -rxt(k,240)*y(k,77)
         mat(k,95) = -(rxt(k,242) + rxt(k,246)) * y(k,77)
         mat(k,38) = -(rxt(k,248) + rxt(k,261)) * y(k,77)
         mat(k,57) = -rxt(k,250)*y(k,77)
         mat(k,118) = -rxt(k,255)*y(k,77)
         mat(k,276) = -rxt(k,260)*y(k,77)
         mat(k,239) = -rxt(k,262)*y(k,77)
         mat(k,559) = mat(k,559) + rxt(k,236)*y(k,56)
         mat(k,62) = mat(k,62) + .300_r8*rxt(k,240)*y(k,77)
         mat(k,481) = rxt(k,160)*y(k,68)
         mat(k,130) = rxt(k,214)*y(k,81)
         mat(k,285) = rxt(k,119)*y(k,58) + 2.000_r8*rxt(k,114)*y(k,68)
         mat(k,690) = mat(k,690) + rxt(k,111)*y(k,56) + rxt(k,94)*y(k,73)
         mat(k,101) = mat(k,101) + rxt(k,112)*y(k,56)
         mat(k,155) = mat(k,155) + rxt(k,204)*y(k,56) + rxt(k,210)*y(k,73)
         mat(k,322) = mat(k,322) + rxt(k,175)*y(k,56) + rxt(k,187)*y(k,73)
         mat(k,139) = rxt(k,206)*y(k,56)
         mat(k,148) = mat(k,148) + rxt(k,178)*y(k,56)
         mat(k,669) = rxt(k,144)*y(k,68)
         mat(k,360) = mat(k,360) + rxt(k,139)*y(k,68)
         mat(k,428) = mat(k,428) + rxt(k,236)*y(k,15) + rxt(k,111)*y(k,36) &
                      + rxt(k,112)*y(k,38) + rxt(k,204)*y(k,40) + rxt(k,175)*y(k,44) &
                      + rxt(k,206)*y(k,48) + rxt(k,178)*y(k,49) + rxt(k,117)*y(k,68)
         mat(k,581) = mat(k,581) + rxt(k,119)*y(k,35) + rxt(k,118)*y(k,68)
         mat(k,385) = mat(k,385) + rxt(k,160)*y(k,23) + 2.000_r8*rxt(k,114)*y(k,35) &
                      + rxt(k,144)*y(k,53) + rxt(k,139)*y(k,55) + rxt(k,117)*y(k,56) &
                      + rxt(k,118)*y(k,58)
         mat(k,537) = rxt(k,94)*y(k,36) + rxt(k,210)*y(k,40) + rxt(k,187)*y(k,44) &
                      + 2.000_r8*rxt(k,95)*y(k,81)
         mat(k,646) = mat(k,646) + .300_r8*rxt(k,240)*y(k,20)
         mat(k,741) = rxt(k,214)*y(k,34) + 2.000_r8*rxt(k,95)*y(k,73)
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
         mat(k,170) = -(rxt(k,302)*y(k,69) + rxt(k,303)*y(k,57))
         mat(k,207) = -rxt(k,302)*y(k,78)
         mat(k,494) = -rxt(k,303)*y(k,78)
         mat(k,405) = rxt(k,293)*y(k,70) + rxt(k,294)*y(k,72) + rxt(k,306)*y(k,79) &
                      + rxt(k,312)*y(k,80)
         mat(k,194) = rxt(k,304)*y(k,79) + rxt(k,309)*y(k,80)
         mat(k,105) = rxt(k,293)*y(k,56)
         mat(k,160) = rxt(k,294)*y(k,56)
         mat(k,80) = rxt(k,306)*y(k,56) + rxt(k,304)*y(k,67)
         mat(k,75) = rxt(k,312)*y(k,56) + rxt(k,309)*y(k,67)
         mat(k,78) = -(rxt(k,304)*y(k,67) + rxt(k,306)*y(k,56) + rxt(k,307)*y(k,57))
         mat(k,190) = -rxt(k,304)*y(k,79)
         mat(k,394) = -rxt(k,306)*y(k,79)
         mat(k,489) = -rxt(k,307)*y(k,79)
         mat(k,190) = mat(k,190) + rxt(k,308)*y(k,80)
         mat(k,72) = rxt(k,308)*y(k,67)
         mat(k,71) = -((rxt(k,308) + rxt(k,309)) * y(k,67) + rxt(k,312)*y(k,56))
         mat(k,189) = -(rxt(k,308) + rxt(k,309)) * y(k,80)
         mat(k,393) = -rxt(k,312)*y(k,80)
         mat(k,745) = -(rxt(k,95)*y(k,73) + rxt(k,214)*y(k,34) + rxt(k,263)*y(k,65))
         mat(k,541) = -rxt(k,95)*y(k,81)
         mat(k,132) = -rxt(k,214)*y(k,81)
         mat(k,32) = -rxt(k,263)*y(k,81)
         mat(k,563) = rxt(k,237)*y(k,77)
         mat(k,63) = rxt(k,240)*y(k,77)
         mat(k,287) = rxt(k,115)*y(k,68)
         mat(k,694) = rxt(k,120)*y(k,77)
         mat(k,102) = rxt(k,121)*y(k,77)
         mat(k,157) = rxt(k,205)*y(k,77)
         mat(k,325) = (rxt(k,276)+rxt(k,281))*y(k,48) + (rxt(k,269)+rxt(k,275) &
                       +rxt(k,280))*y(k,49) + rxt(k,176)*y(k,77)
         mat(k,262) = rxt(k,152)*y(k,77)
         mat(k,70) = rxt(k,128)*y(k,77)
         mat(k,141) = (rxt(k,276)+rxt(k,281))*y(k,44)
         mat(k,149) = (rxt(k,269)+rxt(k,275)+rxt(k,280))*y(k,44) + rxt(k,179)*y(k,77)
         mat(k,389) = rxt(k,115)*y(k,35) + rxt(k,122)*y(k,77)
         mat(k,650) = rxt(k,237)*y(k,15) + rxt(k,240)*y(k,20) + rxt(k,120)*y(k,36) &
                      + rxt(k,121)*y(k,38) + rxt(k,205)*y(k,40) + rxt(k,176)*y(k,44) &
                      + rxt(k,152)*y(k,46) + rxt(k,128)*y(k,47) + rxt(k,179)*y(k,49) &
                      + rxt(k,122)*y(k,68) + 2.000_r8*rxt(k,125)*y(k,77)
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
         mat(k, 5) = mat(k, 5) + lmat(k, 5)
         mat(k, 6) = mat(k, 6) + lmat(k, 6)
         mat(k, 8) = lmat(k, 8)
         mat(k, 9) = lmat(k, 9)
         mat(k, 10) = lmat(k, 10)
         mat(k, 11) = mat(k, 11) + lmat(k, 11)
         mat(k, 12) = mat(k, 12) + lmat(k, 12)
         mat(k, 13) = mat(k, 13) + lmat(k, 13)
         mat(k, 14) = mat(k, 14) + lmat(k, 14)
         mat(k, 15) = lmat(k, 15)
         mat(k, 16) = lmat(k, 16)
         mat(k, 17) = lmat(k, 17)
         mat(k, 18) = mat(k, 18) + lmat(k, 18)
         mat(k, 19) = mat(k, 19) + lmat(k, 19)
         mat(k, 20) = mat(k, 20) + lmat(k, 20)
         mat(k, 22) = lmat(k, 22)
         mat(k, 23) = lmat(k, 23)
         mat(k, 24) = lmat(k, 24)
         mat(k, 25) = lmat(k, 25)
         mat(k, 26) = lmat(k, 26)
         mat(k, 27) = lmat(k, 27)
         mat(k, 29) = mat(k, 29) + lmat(k, 29)
         mat(k, 30) = lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 39) = lmat(k, 39)
         mat(k, 40) = lmat(k, 40)
         mat(k, 41) = lmat(k, 41)
         mat(k, 42) = lmat(k, 42)
         mat(k, 43) = lmat(k, 43)
         mat(k, 44) = lmat(k, 44)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 48) = lmat(k, 48)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 51) = mat(k, 51) + lmat(k, 51)
         mat(k, 52) = lmat(k, 52)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = lmat(k, 59)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 82) = mat(k, 82) + lmat(k, 82)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 85) = lmat(k, 85)
         mat(k, 86) = lmat(k, 86)
         mat(k, 87) = lmat(k, 87)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 90) = lmat(k, 90)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 92) = mat(k, 92) + lmat(k, 92)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 128) = lmat(k, 128)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 135) = lmat(k, 135)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 174) = lmat(k, 174)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 220) = lmat(k, 220)
         mat(k, 223) = lmat(k, 223)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 236) = lmat(k, 236)
         mat(k, 238) = lmat(k, 238)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 250) = lmat(k, 250)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 259) = lmat(k, 259)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 264) = lmat(k, 264)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 359) = mat(k, 359) + lmat(k, 359)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 405) = mat(k, 405) + lmat(k, 405)
         mat(k, 407) = lmat(k, 407)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 462) = mat(k, 462) + lmat(k, 462)
         mat(k, 468) = lmat(k, 468)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 470) = lmat(k, 470)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 475) = mat(k, 475) + lmat(k, 475)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = mat(k, 489) + lmat(k, 489)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 496) = lmat(k, 496)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 520) = mat(k, 520) + lmat(k, 520)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 526) = lmat(k, 526)
         mat(k, 528) = lmat(k, 528)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 531) = mat(k, 531) + lmat(k, 531)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 534) = lmat(k, 534)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 538) = lmat(k, 538)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 545) = mat(k, 545) + lmat(k, 545)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 561) = lmat(k, 561)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 577) = mat(k, 577) + lmat(k, 577)
         mat(k, 579) = mat(k, 579) + lmat(k, 579)
         mat(k, 593) = mat(k, 593) + lmat(k, 593)
         mat(k, 601) = mat(k, 601) + lmat(k, 601)
         mat(k, 608) = mat(k, 608) + lmat(k, 608)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 610) = mat(k, 610) + lmat(k, 610)
         mat(k, 614) = lmat(k, 614)
         mat(k, 615) = lmat(k, 615)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 635) = mat(k, 635) + lmat(k, 635)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 646) = mat(k, 646) + lmat(k, 646)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 651) = mat(k, 651) + lmat(k, 651)
         mat(k, 653) = lmat(k, 653)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 661) = mat(k, 661) + lmat(k, 661)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 692) = mat(k, 692) + lmat(k, 692)
         mat(k, 703) = mat(k, 703) + lmat(k, 703)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 719) = mat(k, 719) + lmat(k, 719)
         mat(k, 728) = lmat(k, 728)
         mat(k, 733) = lmat(k, 733)
         mat(k, 737) = mat(k, 737) + lmat(k, 737)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 743) = lmat(k, 743)
         mat(k, 745) = mat(k, 745) + lmat(k, 745)
         mat(k, 138) = 0._r8
         mat(k, 162) = 0._r8
         mat(k, 163) = 0._r8
         mat(k, 172) = 0._r8
         mat(k, 175) = 0._r8
         mat(k, 176) = 0._r8
         mat(k, 179) = 0._r8
         mat(k, 183) = 0._r8
         mat(k, 193) = 0._r8
         mat(k, 199) = 0._r8
         mat(k, 200) = 0._r8
         mat(k, 202) = 0._r8
         mat(k, 204) = 0._r8
         mat(k, 205) = 0._r8
         mat(k, 208) = 0._r8
         mat(k, 209) = 0._r8
         mat(k, 212) = 0._r8
         mat(k, 213) = 0._r8
         mat(k, 217) = 0._r8
         mat(k, 221) = 0._r8
         mat(k, 224) = 0._r8
         mat(k, 227) = 0._r8
         mat(k, 230) = 0._r8
         mat(k, 240) = 0._r8
         mat(k, 252) = 0._r8
         mat(k, 256) = 0._r8
         mat(k, 257) = 0._r8
         mat(k, 261) = 0._r8
         mat(k, 269) = 0._r8
         mat(k, 279) = 0._r8
         mat(k, 290) = 0._r8
         mat(k, 293) = 0._r8
         mat(k, 295) = 0._r8
         mat(k, 298) = 0._r8
         mat(k, 299) = 0._r8
         mat(k, 301) = 0._r8
         mat(k, 309) = 0._r8
         mat(k, 311) = 0._r8
         mat(k, 312) = 0._r8
         mat(k, 313) = 0._r8
         mat(k, 317) = 0._r8
         mat(k, 319) = 0._r8
         mat(k, 320) = 0._r8
         mat(k, 321) = 0._r8
         mat(k, 323) = 0._r8
         mat(k, 324) = 0._r8
         mat(k, 327) = 0._r8
         mat(k, 330) = 0._r8
         mat(k, 335) = 0._r8
         mat(k, 337) = 0._r8
         mat(k, 339) = 0._r8
         mat(k, 340) = 0._r8
         mat(k, 346) = 0._r8
         mat(k, 347) = 0._r8
         mat(k, 348) = 0._r8
         mat(k, 349) = 0._r8
         mat(k, 353) = 0._r8
         mat(k, 354) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 358) = 0._r8
         mat(k, 362) = 0._r8
         mat(k, 363) = 0._r8
         mat(k, 364) = 0._r8
         mat(k, 381) = 0._r8
         mat(k, 382) = 0._r8
         mat(k, 399) = 0._r8
         mat(k, 406) = 0._r8
         mat(k, 410) = 0._r8
         mat(k, 412) = 0._r8
         mat(k, 417) = 0._r8
         mat(k, 424) = 0._r8
         mat(k, 432) = 0._r8
         mat(k, 440) = 0._r8
         mat(k, 442) = 0._r8
         mat(k, 446) = 0._r8
         mat(k, 452) = 0._r8
         mat(k, 454) = 0._r8
         mat(k, 458) = 0._r8
         mat(k, 460) = 0._r8
         mat(k, 466) = 0._r8
         mat(k, 473) = 0._r8
         mat(k, 477) = 0._r8
         mat(k, 480) = 0._r8
         mat(k, 482) = 0._r8
         mat(k, 484) = 0._r8
         mat(k, 485) = 0._r8
         mat(k, 502) = 0._r8
         mat(k, 505) = 0._r8
         mat(k, 506) = 0._r8
         mat(k, 509) = 0._r8
         mat(k, 511) = 0._r8
         mat(k, 512) = 0._r8
         mat(k, 514) = 0._r8
         mat(k, 515) = 0._r8
         mat(k, 516) = 0._r8
         mat(k, 522) = 0._r8
         mat(k, 527) = 0._r8
         mat(k, 536) = 0._r8
         mat(k, 548) = 0._r8
         mat(k, 552) = 0._r8
         mat(k, 554) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 557) = 0._r8
         mat(k, 558) = 0._r8
         mat(k, 560) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 578) = 0._r8
         mat(k, 583) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 589) = 0._r8
         mat(k, 595) = 0._r8
         mat(k, 596) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 598) = 0._r8
         mat(k, 603) = 0._r8
         mat(k, 605) = 0._r8
         mat(k, 606) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 654) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 671) = 0._r8
         mat(k, 673) = 0._r8
         mat(k, 676) = 0._r8
         mat(k, 679) = 0._r8
         mat(k, 680) = 0._r8
         mat(k, 681) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 685) = 0._r8
         mat(k, 687) = 0._r8
         mat(k, 688) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 693) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 704) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 706) = 0._r8
         mat(k, 712) = 0._r8
         mat(k, 713) = 0._r8
         mat(k, 714) = 0._r8
         mat(k, 718) = 0._r8
         mat(k, 720) = 0._r8
         mat(k, 725) = 0._r8
         mat(k, 726) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 729) = 0._r8
         mat(k, 730) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 732) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 735) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 739) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 744) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 2) = mat(k, 2) - dti(k)
         mat(k, 5) = mat(k, 5) - dti(k)
         mat(k, 8) = mat(k, 8) - dti(k)
         mat(k, 11) = mat(k, 11) - dti(k)
         mat(k, 14) = mat(k, 14) - dti(k)
         mat(k, 16) = mat(k, 16) - dti(k)
         mat(k, 18) = mat(k, 18) - dti(k)
         mat(k, 22) = mat(k, 22) - dti(k)
         mat(k, 25) = mat(k, 25) - dti(k)
         mat(k, 29) = mat(k, 29) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 50) = mat(k, 50) - dti(k)
         mat(k, 58) = mat(k, 58) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 78) = mat(k, 78) - dti(k)
         mat(k, 84) = mat(k, 84) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 143) = mat(k, 143) - dti(k)
         mat(k, 150) = mat(k, 150) - dti(k)
         mat(k, 159) = mat(k, 159) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 181) = mat(k, 181) - dti(k)
         mat(k, 196) = mat(k, 196) - dti(k)
         mat(k, 210) = mat(k, 210) - dti(k)
         mat(k, 225) = mat(k, 225) - dti(k)
         mat(k, 235) = mat(k, 235) - dti(k)
         mat(k, 243) = mat(k, 243) - dti(k)
         mat(k, 255) = mat(k, 255) - dti(k)
         mat(k, 266) = mat(k, 266) - dti(k)
         mat(k, 280) = mat(k, 280) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 376) = mat(k, 376) - dti(k)
         mat(k, 420) = mat(k, 420) - dti(k)
         mat(k, 449) = mat(k, 449) - dti(k)
         mat(k, 475) = mat(k, 475) - dti(k)
         mat(k, 507) = mat(k, 507) - dti(k)
         mat(k, 533) = mat(k, 533) - dti(k)
         mat(k, 556) = mat(k, 556) - dti(k)
         mat(k, 579) = mat(k, 579) - dti(k)
         mat(k, 608) = mat(k, 608) - dti(k)
         mat(k, 646) = mat(k, 646) - dti(k)
         mat(k, 670) = mat(k, 670) - dti(k)
         mat(k, 692) = mat(k, 692) - dti(k)
         mat(k, 719) = mat(k, 719) - dti(k)
         mat(k, 745) = mat(k, 745) - dti(k)
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
