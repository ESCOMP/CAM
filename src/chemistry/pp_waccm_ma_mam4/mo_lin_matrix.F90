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
         mat(k,1) = -( het_rates(k,1) )
         mat(k,2) = -( het_rates(k,2) )
         mat(k,310) = -( het_rates(k,3) )
         mat(k,45) = rxt(k,26)
         mat(k,700) = rxt(k,27)
         mat(k,106) = rxt(k,29)
         mat(k,171) = rxt(k,51)
         mat(k,145) = rxt(k,57)
         mat(k,543) = rxt(k,181)*y(k,9) + rxt(k,207)*y(k,10) + 3.000_r8*rxt(k,208)*y(k,24) &
                      + 2.000_r8*rxt(k,209)*y(k,42) + 2.000_r8*rxt(k,230)*y(k,16) &
                      + rxt(k,231)*y(k,18)
         mat(k,416) = 2.000_r8*rxt(k,218)*y(k,16) + rxt(k,220)*y(k,18) &
                      + 3.000_r8*rxt(k,225)*y(k,24)
         mat(k,652) = 2.000_r8*rxt(k,219)*y(k,16) + rxt(k,221)*y(k,18) &
                      + 3.000_r8*rxt(k,226)*y(k,24)
         mat(k,44) = -( rxt(k,26) + het_rates(k,4) )
         mat(k,715) = -( rxt(k,27) + het_rates(k,5) )
         mat(k,110) = rxt(k,28)
         mat(k,103) = -( rxt(k,28) + rxt(k,29) + rxt(k,270) + rxt(k,273) + rxt(k,278) &
                 + het_rates(k,6) )
         mat(k,738) = -( rxt(k,20) + rxt(k,21) + het_rates(k,17) )
         mat(k,81) = rxt(k,22)
         mat(k,559) = rxt(k,244)*y(k,23) + rxt(k,245)*y(k,23)
         mat(k,347) = -( het_rates(k,21) )
         mat(k,418) = rxt(k,156)*y(k,23)
         mat(k,156) = rxt(k,212)*y(k,23)
         mat(k,654) = rxt(k,241)*y(k,23)
         mat(k,545) = rxt(k,243)*y(k,23)
         mat(k,77) = -( rxt(k,22) + het_rates(k,22) )
         mat(k,421) = -( rxt(k,156)*y(k,23) + rxt(k,218)*y(k,16) + rxt(k,220)*y(k,18) &
                      + rxt(k,223)*y(k,20) + rxt(k,225)*y(k,24) + het_rates(k,25) )
         mat(k,46) = rxt(k,26)
         mat(k,36) = 2.000_r8*rxt(k,43)
         mat(k,22) = 2.000_r8*rxt(k,44)
         mat(k,620) = rxt(k,45)
         mat(k,266) = rxt(k,46)
         mat(k,39) = rxt(k,49)
         mat(k,333) = rxt(k,55)
         mat(k,164) = rxt(k,58)
         mat(k,548) = 4.000_r8*rxt(k,180)*y(k,8) + rxt(k,181)*y(k,9) &
                      + 2.000_r8*rxt(k,182)*y(k,11) + 2.000_r8*rxt(k,183)*y(k,12) &
                      + 2.000_r8*rxt(k,184)*y(k,13) + rxt(k,185)*y(k,14) &
                      + 2.000_r8*rxt(k,186)*y(k,15) + rxt(k,232)*y(k,46) &
                      + rxt(k,233)*y(k,47) + rxt(k,234)*y(k,48)
         mat(k,657) = 3.000_r8*rxt(k,222)*y(k,19) + rxt(k,224)*y(k,20) &
                      + rxt(k,227)*y(k,46) + rxt(k,228)*y(k,47) + rxt(k,229)*y(k,48)
         mat(k,35) = -( rxt(k,43) + het_rates(k,26) )
         mat(k,21) = -( rxt(k,44) + rxt(k,190) + het_rates(k,27) )
         mat(k,627) = -( rxt(k,45) + het_rates(k,28) )
         mat(k,269) = rxt(k,47)
         mat(k,68) = rxt(k,59)
         mat(k,23) = 2.000_r8*rxt(k,190)
         mat(k,262) = -( rxt(k,46) + rxt(k,47) + rxt(k,272) + rxt(k,277) + rxt(k,283) &
                 + het_rates(k,29) )
         mat(k,111) = -( het_rates(k,31) )
         mat(k,718) = rxt(k,20) + rxt(k,21)
         mat(k,70) = rxt(k,88)
         mat(k,410) = rxt(k,223)*y(k,20)
         mat(k,186) = rxt(k,300)*y(k,32)
         mat(k,27) = -( rxt(k,48) + het_rates(k,33) )
         mat(k,536) = rxt(k,181)*y(k,9) + rxt(k,183)*y(k,12) + 2.000_r8*rxt(k,184)*y(k,13) &
                      + 2.000_r8*rxt(k,185)*y(k,14) + rxt(k,186)*y(k,15) &
                      + rxt(k,207)*y(k,10) + 2.000_r8*rxt(k,209)*y(k,42) &
                      + rxt(k,233)*y(k,47) + rxt(k,234)*y(k,48)
         mat(k,633) = rxt(k,228)*y(k,47) + rxt(k,229)*y(k,48)
         mat(k,37) = -( rxt(k,49) + het_rates(k,34) )
         mat(k,538) = rxt(k,182)*y(k,11) + rxt(k,183)*y(k,12) + rxt(k,232)*y(k,46)
         mat(k,634) = rxt(k,227)*y(k,46)
         mat(k,58) = -( het_rates(k,35) )
         mat(k,3) = -( het_rates(k,36) )
         mat(k,4) = -( het_rates(k,37) )
         mat(k,5) = -( het_rates(k,38) )
         mat(k,153) = -( rxt(k,212)*y(k,23) + het_rates(k,39) )
         mat(k,28) = 2.000_r8*rxt(k,48)
         mat(k,38) = rxt(k,49)
         mat(k,42) = rxt(k,56)
         mat(k,539) = rxt(k,185)*y(k,14) + rxt(k,207)*y(k,10)
         mat(k,299) = -( het_rates(k,40) )
         mat(k,747) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,721) = 2.000_r8*rxt(k,20)
         mat(k,78) = rxt(k,22)
         mat(k,170) = rxt(k,51)
         mat(k,327) = rxt(k,55)
         mat(k,43) = rxt(k,56)
         mat(k,542) = rxt(k,244)*y(k,23)
         mat(k,574) = -( het_rates(k,41) )
         mat(k,757) = rxt(k,1)
         mat(k,732) = rxt(k,21)
         mat(k,553) = rxt(k,245)*y(k,23)
         mat(k,124) = -( rxt(k,4) + het_rates(k,43) )
         mat(k,386) = .500_r8*rxt(k,264)
         mat(k,24) = -( rxt(k,87) + het_rates(k,44) )
         mat(k,169) = -( rxt(k,51) + het_rates(k,45) )
         mat(k,329) = -( rxt(k,55) + het_rates(k,49) )
         mat(k,417) = rxt(k,156)*y(k,23) + rxt(k,218)*y(k,16) + rxt(k,220)*y(k,18) &
                      + 2.000_r8*rxt(k,223)*y(k,20) + rxt(k,225)*y(k,24)
         mat(k,41) = -( rxt(k,56) + het_rates(k,50) )
         mat(k,152) = rxt(k,212)*y(k,23)
         mat(k,274) = -( rxt(k,9) + het_rates(k,51) )
         mat(k,53) = 2.000_r8*rxt(k,265) + 2.000_r8*rxt(k,268) + 2.000_r8*rxt(k,271) &
                      + 2.000_r8*rxt(k,282)
         mat(k,441) = .500_r8*rxt(k,266)
         mat(k,364) = rxt(k,267)
         mat(k,105) = rxt(k,270) + rxt(k,273) + rxt(k,278)
         mat(k,263) = rxt(k,272) + rxt(k,277) + rxt(k,283)
         mat(k,83) = -( rxt(k,10) + rxt(k,11) + rxt(k,153) + het_rates(k,52) )
         mat(k,144) = -( rxt(k,57) + het_rates(k,53) )
         mat(k,104) = rxt(k,270) + rxt(k,273) + rxt(k,278)
         mat(k,162) = -( rxt(k,58) + het_rates(k,54) )
         mat(k,261) = rxt(k,272) + rxt(k,277) + rxt(k,283)
         mat(k,244) = -( rxt(k,62) + het_rates(k,55) )
         mat(k,586) = rxt(k,15)
         mat(k,193) = rxt(k,301)
         mat(k,52) = -( rxt(k,13) + rxt(k,14) + rxt(k,154) + rxt(k,265) + rxt(k,268) &
                      + rxt(k,271) + rxt(k,282) + het_rates(k,57) )
         mat(k,6) = -( het_rates(k,58) )
         mat(k,7) = -( het_rates(k,59) )
         mat(k,8) = -( het_rates(k,60) )
         mat(k,598) = -( rxt(k,15) + rxt(k,16) + het_rates(k,61) )
         mat(k,57) = rxt(k,14)
         mat(k,455) = rxt(k,17) + .500_r8*rxt(k,266)
         mat(k,377) = rxt(k,19)
         mat(k,207) = rxt(k,298)
         mat(k,96) = rxt(k,311)
         mat(k,554) = 2.000_r8*rxt(k,147)*y(k,56)
         mat(k,450) = -( rxt(k,17) + rxt(k,266) + het_rates(k,62) )
         mat(k,278) = rxt(k,9)
         mat(k,86) = rxt(k,11) + rxt(k,153)
         mat(k,55) = rxt(k,13) + rxt(k,154)
         mat(k,372) = rxt(k,18)
         mat(k,108) = rxt(k,28)
         mat(k,267) = rxt(k,47)
         mat(k,369) = -( rxt(k,18) + rxt(k,19) + rxt(k,267) + het_rates(k,63) )
         mat(k,84) = rxt(k,10)
         mat(k,54) = rxt(k,13) + rxt(k,14) + rxt(k,154)
         mat(k,107) = rxt(k,29)
         mat(k,265) = rxt(k,46)
         mat(k,9) = -( het_rates(k,64) )
         mat(k,10) = -( het_rates(k,65) )
         mat(k,11) = -( het_rates(k,66) )
         mat(k,12) = -( het_rates(k,67) )
         mat(k,494) = -( rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76) + het_rates(k,68) )
         mat(k,754) = rxt(k,2)
         mat(k,525) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,79) + rxt(k,81) &
                      + 2.000_r8*rxt(k,82) + 2.000_r8*rxt(k,83) + rxt(k,84) + rxt(k,85) &
                      + rxt(k,86)
         mat(k,681) = rxt(k,8)
         mat(k,56) = rxt(k,14)
         mat(k,594) = rxt(k,15)
         mat(k,451) = rxt(k,17)
         mat(k,373) = rxt(k,18)
         mat(k,707) = rxt(k,27)
         mat(k,622) = rxt(k,45)
         mat(k,67) = rxt(k,59)
         mat(k,291) = rxt(k,89)
         mat(k,257) = rxt(k,90)
         mat(k,50) = rxt(k,91)
         mat(k,550) = rxt(k,96)
         mat(k,101) = rxt(k,305)
         mat(k,95) = rxt(k,310)
         mat(k,526) = -( rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,78) + rxt(k,79) &
                      + rxt(k,80) + rxt(k,81) + rxt(k,82) + rxt(k,83) + rxt(k,84) &
                      + rxt(k,85) + rxt(k,86) + het_rates(k,69) )
         mat(k,682) = rxt(k,8)
         mat(k,374) = rxt(k,19)
         mat(k,31) = rxt(k,92) + rxt(k,100)
         mat(k,34) = rxt(k,93)
         mat(k,551) = rxt(k,148)*y(k,56)
         mat(k,688) = -( rxt(k,7) + rxt(k,8) + het_rates(k,70) )
         mat(k,64) = -( rxt(k,59) + het_rates(k,71) )
         mat(k,69) = -( rxt(k,88) + het_rates(k,72) )
         mat(k,13) = -( het_rates(k,73) )
         mat(k,14) = -( het_rates(k,74) )
         mat(k,131) = -( het_rates(k,75) )
         mat(k,71) = rxt(k,88)
         mat(k,283) = rxt(k,89)
         mat(k,285) = -( rxt(k,89) + het_rates(k,77) )
         mat(k,255) = rxt(k,90)
         mat(k,254) = -( rxt(k,90) + het_rates(k,78) )
         mat(k,49) = rxt(k,91)
         mat(k,48) = -( rxt(k,91) + het_rates(k,79) )
         mat(k,25) = rxt(k,87)
         mat(k,15) = -( het_rates(k,80) )
         mat(k,16) = -( het_rates(k,81) )
         mat(k,17) = -( het_rates(k,82) )
         mat(k,18) = -( het_rates(k,83) )
         mat(k,19) = -( het_rates(k,84) )
         mat(k,20) = -( het_rates(k,85) )
         mat(k,215) = -( het_rates(k,86) )
         mat(k,584) = rxt(k,16)
         mat(k,242) = rxt(k,62)
         mat(k,479) = rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76)
         mat(k,515) = rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,84) + rxt(k,85) + rxt(k,86)
         mat(k,395) = -( rxt(k,264) + het_rates(k,87) )
         mat(k,85) = rxt(k,11) + rxt(k,153)
         mat(k,420) = rxt(k,220)*y(k,18) + rxt(k,223)*y(k,20)
         mat(k,656) = rxt(k,221)*y(k,18) + rxt(k,224)*y(k,20)
         mat(k,547) = rxt(k,244)*y(k,23)
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
         mat(k,229) = -( het_rates(k,88) )
         mat(k,115) = -( het_rates(k,89) )
         mat(k,98) = rxt(k,305)
         mat(k,92) = rxt(k,310)
         mat(k,138) = -( het_rates(k,90) )
         mat(k,582) = rxt(k,16)
         mat(k,199) = rxt(k,298)
         mat(k,187) = rxt(k,301)
         mat(k,178) = -( het_rates(k,91) )
         mat(k,239) = rxt(k,62)
         mat(k,93) = rxt(k,311)
         mat(k,552) = -( rxt(k,96) + rxt(k,147)*y(k,56) + rxt(k,148)*y(k,56) &
                      + rxt(k,180)*y(k,8) + rxt(k,181)*y(k,9) + rxt(k,182)*y(k,11) &
                      + rxt(k,183)*y(k,12) + rxt(k,184)*y(k,13) + rxt(k,185)*y(k,14) &
                      + rxt(k,186)*y(k,15) + rxt(k,207)*y(k,10) + rxt(k,208)*y(k,24) &
                      + rxt(k,209)*y(k,42) + rxt(k,230)*y(k,16) + rxt(k,231)*y(k,18) &
                      + rxt(k,232)*y(k,46) + rxt(k,233)*y(k,47) + rxt(k,234)*y(k,48) &
                      + rxt(k,243)*y(k,23) + rxt(k,244)*y(k,23) + rxt(k,245)*y(k,23) &
                 + het_rates(k,92) )
         mat(k,756) = rxt(k,1)
         mat(k,527) = rxt(k,6)
         mat(k,683) = rxt(k,7)
         mat(k,30) = -( rxt(k,92) + rxt(k,100) + het_rates(k,93) )
         mat(k,670) = rxt(k,7)
         mat(k,32) = rxt(k,104) + rxt(k,103)*y(k,32)
         mat(k,33) = -( rxt(k,93) + rxt(k,104) + rxt(k,103)*y(k,32) + het_rates(k,94) )
         mat(k,200) = -( rxt(k,298) + het_rates(k,95) )
         mat(k,514) = rxt(k,78) + rxt(k,80)
         mat(k,190) = rxt(k,300)*y(k,32)
         mat(k,665) = -( rxt(k,219)*y(k,16) + rxt(k,221)*y(k,18) + rxt(k,222)*y(k,19) &
                      + rxt(k,224)*y(k,20) + rxt(k,226)*y(k,24) + rxt(k,227)*y(k,46) &
                      + rxt(k,228)*y(k,47) + rxt(k,229)*y(k,48) + rxt(k,241)*y(k,23) &
                 + het_rates(k,96) )
         mat(k,760) = rxt(k,3)
         mat(k,129) = 2.000_r8*rxt(k,4)
         mat(k,280) = rxt(k,9)
         mat(k,88) = rxt(k,10)
         mat(k,80) = rxt(k,22)
         mat(k,149) = rxt(k,57)
         mat(k,167) = rxt(k,58)
         mat(k,457) = .500_r8*rxt(k,266)
         mat(k,556) = rxt(k,243)*y(k,23)
         mat(k,189) = -( rxt(k,301) + rxt(k,300)*y(k,32) + het_rates(k,97) )
         mat(k,477) = rxt(k,73) + rxt(k,74)
         mat(k,513) = rxt(k,79) + rxt(k,81)
         mat(k,94) = rxt(k,285)
         mat(k,99) = rxt(k,286)
         mat(k,97) = -( rxt(k,286) + rxt(k,305) + het_rates(k,98) )
         mat(k,466) = rxt(k,75) + rxt(k,76)
         mat(k,508) = rxt(k,85) + rxt(k,86)
         mat(k,91) = rxt(k,287)
         mat(k,90) = -( rxt(k,285) + rxt(k,287) + rxt(k,310) + rxt(k,311) &
                 + het_rates(k,99) )
         mat(k,465) = rxt(k,71) + rxt(k,72)
         mat(k,507) = rxt(k,77) + rxt(k,84)
         mat(k,764) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,100) )
         mat(k,26) = rxt(k,87)
         mat(k,669) = rxt(k,219)*y(k,16) + rxt(k,221)*y(k,18) + rxt(k,222)*y(k,19) &
                      + rxt(k,224)*y(k,20) + rxt(k,229)*y(k,48) + rxt(k,241)*y(k,23)
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
