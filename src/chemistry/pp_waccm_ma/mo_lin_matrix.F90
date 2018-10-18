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
         mat(k,502) = -( rxt(k,27) + het_rates(k,2) )
         mat(k,69) = rxt(k,28)
         mat(k,64) = -( rxt(k,28) + rxt(k,29) + rxt(k,248) + rxt(k,251) + rxt(k,256) &
                 + het_rates(k,3) )
         mat(k,272) = -( rxt(k,20) + rxt(k,21) + het_rates(k,14) )
         mat(k,35) = rxt(k,22)
         mat(k,430) = rxt(k,239)*y(k,20) + rxt(k,240)*y(k,20)
         mat(k,478) = -( het_rates(k,18) )
         mat(k,338) = rxt(k,151)*y(k,20)
         mat(k,158) = rxt(k,207)*y(k,20)
         mat(k,532) = rxt(k,236)*y(k,20)
         mat(k,439) = rxt(k,238)*y(k,20)
         mat(k,33) = -( rxt(k,22) + het_rates(k,19) )
         mat(k,15) = -( rxt(k,43) + het_rates(k,22) )
         mat(k,1) = -( rxt(k,44) + rxt(k,185) + het_rates(k,23) )
         mat(k,359) = -( rxt(k,45) + het_rates(k,24) )
         mat(k,199) = rxt(k,47)
         mat(k,8) = rxt(k,59)
         mat(k,3) = 2.000_r8*rxt(k,185)
         mat(k,194) = -( rxt(k,46) + rxt(k,47) + rxt(k,250) + rxt(k,255) + rxt(k,261) &
                 + het_rates(k,25) )
         mat(k,72) = -( het_rates(k,27) )
         mat(k,267) = rxt(k,20) + rxt(k,21)
         mat(k,321) = rxt(k,218)*y(k,17)
         mat(k,101) = rxt(k,278)*y(k,28)
         mat(k,4) = -( rxt(k,48) + het_rates(k,29) )
         mat(k,421) = rxt(k,176)*y(k,6) + rxt(k,178)*y(k,9) + 2.000_r8*rxt(k,179)*y(k,10) &
                      + 2.000_r8*rxt(k,180)*y(k,11) + rxt(k,181)*y(k,12) &
                      + rxt(k,202)*y(k,7) + 2.000_r8*rxt(k,204)*y(k,34) &
                      + rxt(k,228)*y(k,38) + rxt(k,229)*y(k,39)
         mat(k,508) = rxt(k,223)*y(k,38) + rxt(k,224)*y(k,39)
         mat(k,17) = -( rxt(k,49) + het_rates(k,30) )
         mat(k,423) = rxt(k,177)*y(k,8) + rxt(k,178)*y(k,9) + rxt(k,227)*y(k,37)
         mat(k,509) = rxt(k,222)*y(k,37)
         mat(k,154) = -( rxt(k,207)*y(k,20) + het_rates(k,31) )
         mat(k,5) = 2.000_r8*rxt(k,48)
         mat(k,18) = rxt(k,49)
         mat(k,25) = rxt(k,56)
         mat(k,424) = rxt(k,180)*y(k,11) + rxt(k,202)*y(k,7)
         mat(k,288) = -( het_rates(k,32) )
         mat(k,239) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,273) = 2.000_r8*rxt(k,20)
         mat(k,36) = rxt(k,22)
         mat(k,186) = rxt(k,51)
         mat(k,255) = rxt(k,55)
         mat(k,26) = rxt(k,56)
         mat(k,431) = rxt(k,239)*y(k,20)
         mat(k,617) = -( het_rates(k,33) )
         mat(k,246) = rxt(k,1)
         mat(k,285) = rxt(k,21)
         mat(k,444) = rxt(k,240)*y(k,20)
         mat(k,83) = -( rxt(k,4) + het_rates(k,35) )
         mat(k,373) = .500_r8*rxt(k,242)
         mat(k,184) = -( rxt(k,51) + het_rates(k,36) )
         mat(k,254) = -( rxt(k,55) + het_rates(k,40) )
         mat(k,328) = rxt(k,151)*y(k,20) + rxt(k,213)*y(k,13) + rxt(k,215)*y(k,15) &
                      + 2.000_r8*rxt(k,218)*y(k,17) + rxt(k,220)*y(k,21)
         mat(k,24) = -( rxt(k,56) + het_rates(k,41) )
         mat(k,153) = rxt(k,207)*y(k,20)
         mat(k,206) = -( rxt(k,9) + het_rates(k,42) )
         mat(k,28) = 2.000_r8*rxt(k,243) + 2.000_r8*rxt(k,246) + 2.000_r8*rxt(k,249) &
                      + 2.000_r8*rxt(k,260)
         mat(k,402) = .500_r8*rxt(k,244)
         mat(k,580) = rxt(k,245)
         mat(k,66) = rxt(k,248) + rxt(k,251) + rxt(k,256)
         mat(k,195) = rxt(k,250) + rxt(k,255) + rxt(k,261)
         mat(k,385) = -( rxt(k,242) + het_rates(k,43) )
         mat(k,42) = rxt(k,11) + rxt(k,148)
         mat(k,334) = rxt(k,215)*y(k,15) + rxt(k,218)*y(k,17)
         mat(k,528) = rxt(k,216)*y(k,15) + rxt(k,219)*y(k,17)
         mat(k,435) = rxt(k,239)*y(k,20)
         mat(k,39) = -( rxt(k,10) + rxt(k,11) + rxt(k,148) + het_rates(k,44) )
         mat(k,112) = -( rxt(k,57) + het_rates(k,45) )
         mat(k,65) = rxt(k,248) + rxt(k,251) + rxt(k,256)
         mat(k,177) = -( rxt(k,58) + het_rates(k,46) )
         mat(k,193) = rxt(k,250) + rxt(k,255) + rxt(k,261)
         mat(k,167) = -( rxt(k,62) + het_rates(k,47) )
         mat(k,302) = rxt(k,15)
         mat(k,105) = rxt(k,279)
         mat(k,27) = -( rxt(k,13) + rxt(k,14) + rxt(k,149) + rxt(k,243) + rxt(k,246) &
                      + rxt(k,249) + rxt(k,260) + het_rates(k,49) )
         mat(k,306) = -( rxt(k,15) + rxt(k,16) + het_rates(k,50) )
         mat(k,29) = rxt(k,14)
         mat(k,407) = rxt(k,17) + .500_r8*rxt(k,244)
         mat(k,585) = rxt(k,19)
         mat(k,125) = rxt(k,276)
         mat(k,51) = rxt(k,288)
         mat(k,432) = 2.000_r8*rxt(k,142)*y(k,48)
         mat(k,411) = -( rxt(k,17) + rxt(k,244) + het_rates(k,51) )
         mat(k,209) = rxt(k,9)
         mat(k,43) = rxt(k,11) + rxt(k,148)
         mat(k,30) = rxt(k,13) + rxt(k,149)
         mat(k,589) = rxt(k,18)
         mat(k,67) = rxt(k,28)
         mat(k,200) = rxt(k,47)
         mat(k,596) = -( rxt(k,18) + rxt(k,19) + rxt(k,245) + het_rates(k,52) )
         mat(k,45) = rxt(k,10)
         mat(k,32) = rxt(k,13) + rxt(k,14) + rxt(k,149)
         mat(k,71) = rxt(k,29)
         mat(k,203) = rxt(k,46)
         mat(k,574) = -( rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76) + het_rates(k,53) )
         mat(k,244) = rxt(k,2)
         mat(k,233) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,79) + rxt(k,81) &
                      + 2.000_r8*rxt(k,82) + 2.000_r8*rxt(k,83) + rxt(k,84) + rxt(k,85) &
                      + rxt(k,86)
         mat(k,633) = rxt(k,8)
         mat(k,31) = rxt(k,14)
         mat(k,316) = rxt(k,15)
         mat(k,417) = rxt(k,17)
         mat(k,595) = rxt(k,18)
         mat(k,504) = rxt(k,27)
         mat(k,367) = rxt(k,45)
         mat(k,9) = rxt(k,59)
         mat(k,442) = rxt(k,91)
         mat(k,58) = rxt(k,282)
         mat(k,52) = rxt(k,287)
         mat(k,226) = -( rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,78) + rxt(k,79) &
                      + rxt(k,80) + rxt(k,81) + rxt(k,82) + rxt(k,83) + rxt(k,84) &
                      + rxt(k,85) + rxt(k,86) + het_rates(k,54) )
         mat(k,621) = rxt(k,8)
         mat(k,581) = rxt(k,19)
         mat(k,11) = rxt(k,87) + rxt(k,95)
         mat(k,14) = rxt(k,88)
         mat(k,427) = rxt(k,143)*y(k,48)
         mat(k,636) = -( rxt(k,7) + rxt(k,8) + het_rates(k,55) )
         mat(k,7) = -( rxt(k,59) + het_rates(k,56) )
         mat(k,458) = -( het_rates(k,58) )
         mat(k,23) = rxt(k,26)
         mat(k,500) = rxt(k,27)
         mat(k,68) = rxt(k,29)
         mat(k,188) = rxt(k,51)
         mat(k,116) = rxt(k,57)
         mat(k,438) = rxt(k,176)*y(k,6) + rxt(k,202)*y(k,7) + 3.000_r8*rxt(k,203)*y(k,21) &
                      + 2.000_r8*rxt(k,204)*y(k,34) + 2.000_r8*rxt(k,225)*y(k,13) &
                      + rxt(k,226)*y(k,15)
         mat(k,337) = 2.000_r8*rxt(k,213)*y(k,13) + rxt(k,215)*y(k,15) &
                      + 3.000_r8*rxt(k,220)*y(k,21)
         mat(k,531) = 2.000_r8*rxt(k,214)*y(k,13) + rxt(k,216)*y(k,15) &
                      + 3.000_r8*rxt(k,221)*y(k,21)
         mat(k,332) = -( rxt(k,151)*y(k,20) + rxt(k,213)*y(k,13) + rxt(k,215)*y(k,15) &
                      + rxt(k,218)*y(k,17) + rxt(k,220)*y(k,21) + het_rates(k,59) )
         mat(k,22) = rxt(k,26)
         mat(k,16) = 2.000_r8*rxt(k,43)
         mat(k,2) = 2.000_r8*rxt(k,44)
         mat(k,358) = rxt(k,45)
         mat(k,198) = rxt(k,46)
         mat(k,19) = rxt(k,49)
         mat(k,256) = rxt(k,55)
         mat(k,180) = rxt(k,58)
         mat(k,433) = 4.000_r8*rxt(k,175)*y(k,5) + rxt(k,176)*y(k,6) &
                      + 2.000_r8*rxt(k,177)*y(k,8) + 2.000_r8*rxt(k,178)*y(k,9) &
                      + 2.000_r8*rxt(k,179)*y(k,10) + rxt(k,180)*y(k,11) &
                      + 2.000_r8*rxt(k,181)*y(k,12) + rxt(k,227)*y(k,37) &
                      + rxt(k,228)*y(k,38) + rxt(k,229)*y(k,39)
         mat(k,526) = 3.000_r8*rxt(k,217)*y(k,16) + rxt(k,219)*y(k,17) &
                      + rxt(k,222)*y(k,37) + rxt(k,223)*y(k,38) + rxt(k,224)*y(k,39)
         mat(k,144) = -( het_rates(k,60) )
         mat(k,301) = rxt(k,16)
         mat(k,166) = rxt(k,62)
         mat(k,553) = rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76)
         mat(k,224) = rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,84) + rxt(k,85) + rxt(k,86)
         mat(k,59) = -( het_rates(k,61) )
         mat(k,91) = -( het_rates(k,62) )
         mat(k,54) = rxt(k,282)
         mat(k,49) = rxt(k,287)
         mat(k,129) = -( het_rates(k,63) )
         mat(k,300) = rxt(k,16)
         mat(k,121) = rxt(k,276)
         mat(k,104) = rxt(k,279)
         mat(k,76) = -( het_rates(k,64) )
         mat(k,162) = rxt(k,62)
         mat(k,48) = rxt(k,288)
         mat(k,437) = -( rxt(k,91) + rxt(k,142)*y(k,48) + rxt(k,143)*y(k,48) &
                      + rxt(k,175)*y(k,5) + rxt(k,176)*y(k,6) + rxt(k,177)*y(k,8) &
                      + rxt(k,178)*y(k,9) + rxt(k,179)*y(k,10) + rxt(k,180)*y(k,11) &
                      + rxt(k,181)*y(k,12) + rxt(k,202)*y(k,7) + rxt(k,203)*y(k,21) &
                      + rxt(k,204)*y(k,34) + rxt(k,225)*y(k,13) + rxt(k,226)*y(k,15) &
                      + rxt(k,227)*y(k,37) + rxt(k,228)*y(k,38) + rxt(k,229)*y(k,39) &
                      + rxt(k,238)*y(k,20) + rxt(k,239)*y(k,20) + rxt(k,240)*y(k,20) &
                 + het_rates(k,65) )
         mat(k,241) = rxt(k,1)
         mat(k,231) = rxt(k,6)
         mat(k,628) = rxt(k,7)
         mat(k,10) = -( rxt(k,87) + rxt(k,95) + het_rates(k,66) )
         mat(k,619) = rxt(k,7)
         mat(k,12) = rxt(k,99) + rxt(k,98)*y(k,28)
         mat(k,13) = -( rxt(k,88) + rxt(k,99) + rxt(k,98)*y(k,28) + het_rates(k,67) )
         mat(k,120) = -( rxt(k,276) + het_rates(k,68) )
         mat(k,222) = rxt(k,78) + rxt(k,80)
         mat(k,103) = rxt(k,278)*y(k,28)
         mat(k,534) = -( rxt(k,214)*y(k,13) + rxt(k,216)*y(k,15) + rxt(k,217)*y(k,16) &
                      + rxt(k,219)*y(k,17) + rxt(k,221)*y(k,21) + rxt(k,222)*y(k,37) &
                      + rxt(k,223)*y(k,38) + rxt(k,224)*y(k,39) + rxt(k,236)*y(k,20) &
                 + het_rates(k,69) )
         mat(k,243) = rxt(k,3)
         mat(k,88) = 2.000_r8*rxt(k,4)
         mat(k,211) = rxt(k,9)
         mat(k,44) = rxt(k,10)
         mat(k,38) = rxt(k,22)
         mat(k,118) = rxt(k,57)
         mat(k,182) = rxt(k,58)
         mat(k,416) = .500_r8*rxt(k,244)
         mat(k,441) = rxt(k,238)*y(k,20)
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
         mat(k,102) = -( rxt(k,279) + rxt(k,278)*y(k,28) + het_rates(k,70) )
         mat(k,549) = rxt(k,71) + rxt(k,72)
         mat(k,221) = rxt(k,79) + rxt(k,81)
         mat(k,50) = rxt(k,263)
         mat(k,55) = rxt(k,264)
         mat(k,53) = -( rxt(k,264) + rxt(k,282) + het_rates(k,71) )
         mat(k,542) = rxt(k,73) + rxt(k,75)
         mat(k,217) = rxt(k,84) + rxt(k,86)
         mat(k,47) = rxt(k,265)
         mat(k,46) = -( rxt(k,263) + rxt(k,265) + rxt(k,287) + rxt(k,288) &
                 + het_rates(k,72) )
         mat(k,541) = rxt(k,74) + rxt(k,76)
         mat(k,216) = rxt(k,77) + rxt(k,85)
         mat(k,238) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,73) )
         mat(k,521) = rxt(k,214)*y(k,13) + rxt(k,216)*y(k,15) + rxt(k,217)*y(k,16) &
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
