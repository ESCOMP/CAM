      module mo_lin_matrix
      use chem_mods, only: veclen
      private
      public :: linmat
      contains
      subroutine linmat01( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
      do k = 1,avec_len
         mat(k,21) = -( rxt(k,26) + het_rates(k,1) )
         mat(k,363) = -( rxt(k,27) + het_rates(k,2) )
         mat(k,63) = rxt(k,28)
         mat(k,59) = -( rxt(k,28) + rxt(k,29) + rxt(k,248) + rxt(k,251) + rxt(k,256) &
                 + het_rates(k,3) )
         mat(k,282) = -( rxt(k,20) + rxt(k,21) + het_rates(k,14) )
         mat(k,35) = rxt(k,22)
         mat(k,435) = rxt(k,239)*y(k,20) + rxt(k,240)*y(k,20)
         mat(k,483) = -( het_rates(k,18) )
         mat(k,393) = rxt(k,151)*y(k,20)
         mat(k,107) = rxt(k,207)*y(k,20)
         mat(k,539) = rxt(k,236)*y(k,20)
         mat(k,444) = rxt(k,238)*y(k,20)
         mat(k,33) = -( rxt(k,22) + het_rates(k,19) )
         mat(k,15) = -( rxt(k,43) + het_rates(k,22) )
         mat(k,1) = -( rxt(k,44) + rxt(k,185) + het_rates(k,23) )
         mat(k,416) = -( rxt(k,45) + het_rates(k,24) )
         mat(k,210) = rxt(k,47)
         mat(k,8) = rxt(k,59)
         mat(k,3) = 2.000_r8*rxt(k,185)
         mat(k,204) = -( rxt(k,46) + rxt(k,47) + rxt(k,250) + rxt(k,255) + rxt(k,261) &
                 + het_rates(k,25) )
         mat(k,67) = -( het_rates(k,27) )
         mat(k,277) = rxt(k,20) + rxt(k,21)
         mat(k,376) = rxt(k,218)*y(k,17)
         mat(k,135) = rxt(k,278)*y(k,28)
         mat(k,4) = -( rxt(k,48) + het_rates(k,29) )
         mat(k,426) = rxt(k,176)*y(k,6) + rxt(k,178)*y(k,9) + 2.000_r8*rxt(k,179)*y(k,10) &
                      + 2.000_r8*rxt(k,180)*y(k,11) + rxt(k,181)*y(k,12) &
                      + rxt(k,202)*y(k,7) + 2.000_r8*rxt(k,204)*y(k,34) &
                      + rxt(k,228)*y(k,38) + rxt(k,229)*y(k,39)
         mat(k,515) = rxt(k,223)*y(k,38) + rxt(k,224)*y(k,39)
         mat(k,17) = -( rxt(k,49) + het_rates(k,30) )
         mat(k,428) = rxt(k,177)*y(k,8) + rxt(k,178)*y(k,9) + rxt(k,227)*y(k,37)
         mat(k,516) = rxt(k,222)*y(k,37)
         mat(k,102) = -( rxt(k,207)*y(k,20) + het_rates(k,31) )
         mat(k,5) = 2.000_r8*rxt(k,48)
         mat(k,18) = rxt(k,49)
         mat(k,25) = rxt(k,56)
         mat(k,429) = rxt(k,180)*y(k,11) + rxt(k,202)*y(k,7)
         mat(k,298) = -( het_rates(k,32) )
         mat(k,249) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,283) = 2.000_r8*rxt(k,20)
         mat(k,36) = rxt(k,22)
         mat(k,120) = rxt(k,51)
         mat(k,265) = rxt(k,55)
         mat(k,26) = rxt(k,56)
         mat(k,436) = rxt(k,239)*y(k,20)
         mat(k,628) = -( het_rates(k,33) )
         mat(k,256) = rxt(k,1)
         mat(k,295) = rxt(k,21)
         mat(k,449) = rxt(k,240)*y(k,20)
         mat(k,71) = -( rxt(k,4) + het_rates(k,35) )
         mat(k,548) = .500_r8*rxt(k,242)
         mat(k,118) = -( rxt(k,51) + het_rates(k,36) )
         mat(k,264) = -( rxt(k,55) + het_rates(k,40) )
         mat(k,383) = rxt(k,151)*y(k,20) + rxt(k,213)*y(k,13) + rxt(k,215)*y(k,15) &
                      + 2.000_r8*rxt(k,218)*y(k,17) + rxt(k,220)*y(k,21)
         mat(k,24) = -( rxt(k,56) + het_rates(k,41) )
         mat(k,101) = rxt(k,207)*y(k,20)
         mat(k,216) = -( rxt(k,9) + het_rates(k,42) )
         mat(k,28) = 2.000_r8*rxt(k,243) + 2.000_r8*rxt(k,246) + 2.000_r8*rxt(k,249) &
                      + 2.000_r8*rxt(k,260)
         mat(k,496) = .500_r8*rxt(k,244)
         mat(k,311) = rxt(k,245)
         mat(k,61) = rxt(k,248) + rxt(k,251) + rxt(k,256)
         mat(k,205) = rxt(k,250) + rxt(k,255) + rxt(k,261)
         mat(k,567) = -( rxt(k,242) + het_rates(k,43) )
         mat(k,45) = rxt(k,11) + rxt(k,148)
         mat(k,396) = rxt(k,215)*y(k,15) + rxt(k,218)*y(k,17)
         mat(k,542) = rxt(k,216)*y(k,15) + rxt(k,219)*y(k,17)
         mat(k,447) = rxt(k,239)*y(k,20)
         mat(k,39) = -( rxt(k,10) + rxt(k,11) + rxt(k,148) + het_rates(k,44) )
         mat(k,93) = -( rxt(k,57) + het_rates(k,45) )
         mat(k,60) = rxt(k,248) + rxt(k,251) + rxt(k,256)
         mat(k,111) = -( rxt(k,58) + het_rates(k,46) )
         mat(k,203) = rxt(k,250) + rxt(k,255) + rxt(k,261)
         mat(k,193) = -( rxt(k,62) + het_rates(k,47) )
         mat(k,334) = rxt(k,15)
         mat(k,142) = rxt(k,279)
         mat(k,27) = -( rxt(k,13) + rxt(k,14) + rxt(k,149) + rxt(k,243) + rxt(k,246) &
                      + rxt(k,249) + rxt(k,260) + het_rates(k,49) )
         mat(k,339) = -( rxt(k,15) + rxt(k,16) + het_rates(k,50) )
         mat(k,30) = rxt(k,14)
         mat(k,502) = rxt(k,17) + .500_r8*rxt(k,244)
         mat(k,317) = rxt(k,19)
         mat(k,154) = rxt(k,276)
         mat(k,51) = rxt(k,289)
         mat(k,438) = 2.000_r8*rxt(k,142)*y(k,48)
         mat(k,509) = -( rxt(k,17) + rxt(k,244) + het_rates(k,51) )
         mat(k,221) = rxt(k,9)
         mat(k,43) = rxt(k,11) + rxt(k,148)
         mat(k,31) = rxt(k,13) + rxt(k,149)
         mat(k,324) = rxt(k,18)
         mat(k,65) = rxt(k,28)
         mat(k,211) = rxt(k,47)
         mat(k,316) = -( rxt(k,18) + rxt(k,19) + rxt(k,245) + het_rates(k,52) )
         mat(k,42) = rxt(k,10)
         mat(k,29) = rxt(k,13) + rxt(k,14) + rxt(k,149)
         mat(k,62) = rxt(k,29)
         mat(k,208) = rxt(k,46)
         mat(k,607) = -( rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76) + het_rates(k,53) )
         mat(k,255) = rxt(k,2)
         mat(k,243) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,79) + rxt(k,81) &
                      + 2.000_r8*rxt(k,82) + 2.000_r8*rxt(k,83) + rxt(k,84) + rxt(k,85) &
                      + rxt(k,86)
         mat(k,645) = rxt(k,8)
         mat(k,32) = rxt(k,14)
         mat(k,349) = rxt(k,15)
         mat(k,512) = rxt(k,17)
         mat(k,327) = rxt(k,18)
         mat(k,372) = rxt(k,27)
         mat(k,423) = rxt(k,45)
         mat(k,9) = rxt(k,59)
         mat(k,448) = rxt(k,91)
         mat(k,58) = rxt(k,283)
         mat(k,52) = rxt(k,288)
         mat(k,236) = -( rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,78) + rxt(k,79) &
                      + rxt(k,80) + rxt(k,81) + rxt(k,82) + rxt(k,83) + rxt(k,84) &
                      + rxt(k,85) + rxt(k,86) + het_rates(k,54) )
         mat(k,632) = rxt(k,8)
         mat(k,312) = rxt(k,19)
         mat(k,11) = rxt(k,87) + rxt(k,95)
         mat(k,14) = rxt(k,88)
         mat(k,432) = rxt(k,143)*y(k,48)
         mat(k,647) = -( rxt(k,7) + rxt(k,8) + het_rates(k,55) )
         mat(k,7) = -( rxt(k,59) + het_rates(k,56) )
         mat(k,463) = -( het_rates(k,58) )
         mat(k,23) = rxt(k,26)
         mat(k,367) = rxt(k,27)
         mat(k,64) = rxt(k,29)
         mat(k,123) = rxt(k,51)
         mat(k,98) = rxt(k,57)
         mat(k,443) = rxt(k,176)*y(k,6) + rxt(k,202)*y(k,7) + 3.000_r8*rxt(k,203)*y(k,21) &
                      + 2.000_r8*rxt(k,204)*y(k,34) + 2.000_r8*rxt(k,225)*y(k,13) &
                      + rxt(k,226)*y(k,15)
         mat(k,392) = 2.000_r8*rxt(k,213)*y(k,13) + rxt(k,215)*y(k,15) &
                      + 3.000_r8*rxt(k,220)*y(k,21)
         mat(k,538) = 2.000_r8*rxt(k,214)*y(k,13) + rxt(k,216)*y(k,15) &
                      + 3.000_r8*rxt(k,221)*y(k,21)
         mat(k,389) = -( rxt(k,151)*y(k,20) + rxt(k,213)*y(k,13) + rxt(k,215)*y(k,15) &
                      + rxt(k,218)*y(k,17) + rxt(k,220)*y(k,21) + het_rates(k,59) )
         mat(k,22) = rxt(k,26)
         mat(k,16) = 2.000_r8*rxt(k,43)
         mat(k,2) = 2.000_r8*rxt(k,44)
         mat(k,415) = rxt(k,45)
         mat(k,209) = rxt(k,46)
         mat(k,19) = rxt(k,49)
         mat(k,268) = rxt(k,55)
         mat(k,114) = rxt(k,58)
         mat(k,440) = 4.000_r8*rxt(k,175)*y(k,5) + rxt(k,176)*y(k,6) &
                      + 2.000_r8*rxt(k,177)*y(k,8) + 2.000_r8*rxt(k,178)*y(k,9) &
                      + 2.000_r8*rxt(k,179)*y(k,10) + rxt(k,180)*y(k,11) &
                      + 2.000_r8*rxt(k,181)*y(k,12) + rxt(k,227)*y(k,37) &
                      + rxt(k,228)*y(k,38) + rxt(k,229)*y(k,39)
         mat(k,535) = 3.000_r8*rxt(k,217)*y(k,16) + rxt(k,219)*y(k,17) &
                      + rxt(k,222)*y(k,37) + rxt(k,223)*y(k,38) + rxt(k,224)*y(k,39)
         mat(k,164) = -( het_rates(k,60) )
         mat(k,332) = rxt(k,16)
         mat(k,191) = rxt(k,62)
         mat(k,586) = rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76)
         mat(k,233) = rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,84) + rxt(k,85) + rxt(k,86)
         mat(k,178) = -( het_rates(k,61) )
         mat(k,78) = -( het_rates(k,62) )
         mat(k,54) = rxt(k,283)
         mat(k,48) = rxt(k,288)
         mat(k,87) = -( het_rates(k,63) )
         mat(k,330) = rxt(k,16)
         mat(k,148) = rxt(k,276)
         mat(k,136) = rxt(k,279)
         mat(k,127) = -( het_rates(k,64) )
         mat(k,188) = rxt(k,62)
         mat(k,49) = rxt(k,289)
         mat(k,442) = -( rxt(k,91) + rxt(k,142)*y(k,48) + rxt(k,143)*y(k,48) &
                      + rxt(k,175)*y(k,5) + rxt(k,176)*y(k,6) + rxt(k,177)*y(k,8) &
                      + rxt(k,178)*y(k,9) + rxt(k,179)*y(k,10) + rxt(k,180)*y(k,11) &
                      + rxt(k,181)*y(k,12) + rxt(k,202)*y(k,7) + rxt(k,203)*y(k,21) &
                      + rxt(k,204)*y(k,34) + rxt(k,225)*y(k,13) + rxt(k,226)*y(k,15) &
                      + rxt(k,227)*y(k,37) + rxt(k,228)*y(k,38) + rxt(k,229)*y(k,39) &
                      + rxt(k,238)*y(k,20) + rxt(k,239)*y(k,20) + rxt(k,240)*y(k,20) &
                 + het_rates(k,65) )
         mat(k,251) = rxt(k,1)
         mat(k,239) = rxt(k,6)
         mat(k,639) = rxt(k,7)
         mat(k,10) = -( rxt(k,87) + rxt(k,95) + het_rates(k,66) )
         mat(k,630) = rxt(k,7)
         mat(k,12) = rxt(k,99) + rxt(k,98)*y(k,28)
         mat(k,13) = -( rxt(k,88) + rxt(k,99) + rxt(k,98)*y(k,28) + het_rates(k,67) )
         mat(k,149) = -( rxt(k,276) + het_rates(k,68) )
         mat(k,232) = rxt(k,78) + rxt(k,80)
         mat(k,139) = rxt(k,278)*y(k,28)
         mat(k,541) = -( rxt(k,214)*y(k,13) + rxt(k,216)*y(k,15) + rxt(k,217)*y(k,16) &
                      + rxt(k,219)*y(k,17) + rxt(k,221)*y(k,21) + rxt(k,222)*y(k,37) &
                      + rxt(k,223)*y(k,38) + rxt(k,224)*y(k,39) + rxt(k,236)*y(k,20) &
                 + het_rates(k,69) )
         mat(k,254) = rxt(k,3)
         mat(k,75) = 2.000_r8*rxt(k,4)
         mat(k,222) = rxt(k,9)
         mat(k,44) = rxt(k,10)
         mat(k,38) = rxt(k,22)
         mat(k,99) = rxt(k,57)
         mat(k,116) = rxt(k,58)
         mat(k,510) = .500_r8*rxt(k,244)
         mat(k,446) = rxt(k,238)*y(k,20)
      end do
      end subroutine linmat01
      subroutine linmat02( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
      do k = 1,avec_len
         mat(k,138) = -( rxt(k,279) + rxt(k,278)*y(k,28) + het_rates(k,70) )
         mat(k,584) = rxt(k,71) + rxt(k,72)
         mat(k,231) = rxt(k,79) + rxt(k,81)
         mat(k,50) = rxt(k,263)
         mat(k,55) = rxt(k,264)
         mat(k,53) = -( rxt(k,264) + rxt(k,283) + het_rates(k,71) )
         mat(k,574) = rxt(k,73) + rxt(k,75)
         mat(k,227) = rxt(k,84) + rxt(k,86)
         mat(k,47) = rxt(k,265)
         mat(k,46) = -( rxt(k,263) + rxt(k,265) + rxt(k,288) + rxt(k,289) &
                 + het_rates(k,72) )
         mat(k,573) = rxt(k,74) + rxt(k,76)
         mat(k,226) = rxt(k,77) + rxt(k,85)
         mat(k,248) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,73) )
         mat(k,528) = rxt(k,214)*y(k,13) + rxt(k,216)*y(k,15) + rxt(k,217)*y(k,16) &
                      + rxt(k,219)*y(k,17) + rxt(k,224)*y(k,39) + rxt(k,236)*y(k,20)
      end do
      end subroutine linmat02
      subroutine linmat( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
      call linmat01( avec_len, mat, y, rxt, het_rates )
      call linmat02( avec_len, mat, y, rxt, het_rates )
      end subroutine linmat
      end module mo_lin_matrix
