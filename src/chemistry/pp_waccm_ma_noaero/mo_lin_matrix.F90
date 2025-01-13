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
         mat(k,291) = -( het_rates(k,1) )
         mat(k,23) = rxt(k,26)
         mat(k,703) = rxt(k,27)
         mat(k,87) = rxt(k,29)
         mat(k,152) = rxt(k,51)
         mat(k,135) = rxt(k,57)
         mat(k,524) = rxt(k,181)*y(k,7) + rxt(k,207)*y(k,8) + 3.000_r8*rxt(k,208)*y(k,22) &
                      + 2.000_r8*rxt(k,209)*y(k,37) + 2.000_r8*rxt(k,230)*y(k,14) &
                      + rxt(k,231)*y(k,16)
         mat(k,468) = 2.000_r8*rxt(k,218)*y(k,14) + rxt(k,220)*y(k,16) &
                      + 3.000_r8*rxt(k,225)*y(k,22)
         mat(k,633) = 2.000_r8*rxt(k,219)*y(k,14) + rxt(k,221)*y(k,16) &
                      + 3.000_r8*rxt(k,226)*y(k,22)
         mat(k,22) = -( rxt(k,26) + het_rates(k,2) )
         mat(k,719) = -( rxt(k,27) + het_rates(k,3) )
         mat(k,91) = rxt(k,28)
         mat(k,84) = -( rxt(k,28) + rxt(k,29) + rxt(k,270) + rxt(k,273) + rxt(k,278) &
                 + het_rates(k,4) )
         mat(k,556) = -( rxt(k,20) + rxt(k,21) + het_rates(k,15) )
         mat(k,61) = rxt(k,22)
         mat(k,534) = rxt(k,244)*y(k,21) + rxt(k,245)*y(k,21)
         mat(k,328) = -( het_rates(k,19) )
         mat(k,470) = rxt(k,156)*y(k,21)
         mat(k,128) = rxt(k,212)*y(k,21)
         mat(k,635) = rxt(k,241)*y(k,21)
         mat(k,526) = rxt(k,243)*y(k,21)
         mat(k,58) = -( rxt(k,22) + het_rates(k,20) )
         mat(k,475) = -( rxt(k,156)*y(k,21) + rxt(k,218)*y(k,14) + rxt(k,220)*y(k,16) &
                      + rxt(k,223)*y(k,18) + rxt(k,225)*y(k,22) + het_rates(k,23) )
         mat(k,24) = rxt(k,26)
         mat(k,17) = 2.000_r8*rxt(k,43)
         mat(k,4) = 2.000_r8*rxt(k,44)
         mat(k,450) = rxt(k,45)
         mat(k,249) = rxt(k,46)
         mat(k,20) = rxt(k,49)
         mat(k,316) = rxt(k,55)
         mat(k,147) = rxt(k,58)
         mat(k,531) = 4.000_r8*rxt(k,180)*y(k,6) + rxt(k,181)*y(k,7) &
                      + 2.000_r8*rxt(k,182)*y(k,9) + 2.000_r8*rxt(k,183)*y(k,10) &
                      + 2.000_r8*rxt(k,184)*y(k,11) + rxt(k,185)*y(k,12) &
                      + 2.000_r8*rxt(k,186)*y(k,13) + rxt(k,232)*y(k,41) &
                      + rxt(k,233)*y(k,42) + rxt(k,234)*y(k,43)
         mat(k,640) = 3.000_r8*rxt(k,222)*y(k,17) + rxt(k,224)*y(k,18) &
                      + rxt(k,227)*y(k,41) + rxt(k,228)*y(k,42) + rxt(k,229)*y(k,43)
         mat(k,16) = -( rxt(k,43) + het_rates(k,24) )
         mat(k,2) = -( rxt(k,44) + rxt(k,190) + het_rates(k,25) )
         mat(k,449) = -( rxt(k,45) + het_rates(k,26) )
         mat(k,248) = rxt(k,47)
         mat(k,49) = rxt(k,59)
         mat(k,3) = 2.000_r8*rxt(k,190)
         mat(k,243) = -( rxt(k,46) + rxt(k,47) + rxt(k,272) + rxt(k,277) + rxt(k,283) &
                 + het_rates(k,27) )
         mat(k,92) = -( het_rates(k,29) )
         mat(k,542) = rxt(k,20) + rxt(k,21)
         mat(k,51) = rxt(k,88)
         mat(k,462) = rxt(k,223)*y(k,18)
         mat(k,167) = rxt(k,300)*y(k,30)
         mat(k,5) = -( rxt(k,48) + het_rates(k,31) )
         mat(k,517) = rxt(k,181)*y(k,7) + rxt(k,183)*y(k,10) + 2.000_r8*rxt(k,184)*y(k,11) &
                      + 2.000_r8*rxt(k,185)*y(k,12) + rxt(k,186)*y(k,13) &
                      + rxt(k,207)*y(k,8) + 2.000_r8*rxt(k,209)*y(k,37) &
                      + rxt(k,233)*y(k,42) + rxt(k,234)*y(k,43)
         mat(k,614) = rxt(k,228)*y(k,42) + rxt(k,229)*y(k,43)
         mat(k,18) = -( rxt(k,49) + het_rates(k,32) )
         mat(k,519) = rxt(k,182)*y(k,9) + rxt(k,183)*y(k,10) + rxt(k,232)*y(k,41)
         mat(k,615) = rxt(k,227)*y(k,41)
         mat(k,33) = -( het_rates(k,33) )
         mat(k,125) = -( rxt(k,212)*y(k,21) + het_rates(k,34) )
         mat(k,6) = 2.000_r8*rxt(k,48)
         mat(k,19) = rxt(k,49)
         mat(k,26) = rxt(k,56)
         mat(k,520) = rxt(k,185)*y(k,12) + rxt(k,207)*y(k,8)
         mat(k,280) = -( het_rates(k,35) )
         mat(k,728) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,545) = 2.000_r8*rxt(k,20)
         mat(k,59) = rxt(k,22)
         mat(k,151) = rxt(k,51)
         mat(k,308) = rxt(k,55)
         mat(k,27) = rxt(k,56)
         mat(k,523) = rxt(k,244)*y(k,21)
         mat(k,692) = -( het_rates(k,36) )
         mat(k,743) = rxt(k,1)
         mat(k,561) = rxt(k,21)
         mat(k,539) = rxt(k,245)*y(k,21)
         mat(k,96) = -( rxt(k,4) + het_rates(k,38) )
         mat(k,8) = -( rxt(k,87) + het_rates(k,39) )
         mat(k,150) = -( rxt(k,51) + het_rates(k,40) )
         mat(k,310) = -( rxt(k,55) + het_rates(k,44) )
         mat(k,469) = rxt(k,156)*y(k,21) + rxt(k,218)*y(k,14) + rxt(k,220)*y(k,16) &
                      + 2.000_r8*rxt(k,223)*y(k,18) + rxt(k,225)*y(k,22)
         mat(k,25) = -( rxt(k,56) + het_rates(k,45) )
         mat(k,124) = rxt(k,212)*y(k,21)
         mat(k,255) = -( rxt(k,9) + het_rates(k,46) )
         mat(k,40) = 2.000_r8*rxt(k,265) + 2.000_r8*rxt(k,268) + 2.000_r8*rxt(k,271) &
                      + 2.000_r8*rxt(k,282)
         mat(k,593) = .500_r8*rxt(k,266)
         mat(k,345) = rxt(k,267)
         mat(k,86) = rxt(k,270) + rxt(k,273) + rxt(k,278)
         mat(k,244) = rxt(k,272) + rxt(k,277) + rxt(k,283)
         mat(k,64) = -( rxt(k,10) + rxt(k,11) + rxt(k,153) + het_rates(k,47) )
         mat(k,134) = -( rxt(k,57) + het_rates(k,48) )
         mat(k,85) = rxt(k,270) + rxt(k,273) + rxt(k,278)
         mat(k,143) = -( rxt(k,58) + het_rates(k,49) )
         mat(k,242) = rxt(k,272) + rxt(k,277) + rxt(k,283)
         mat(k,225) = -( rxt(k,62) + het_rates(k,50) )
         mat(k,655) = rxt(k,15)
         mat(k,174) = rxt(k,301)
         mat(k,39) = -( rxt(k,13) + rxt(k,14) + rxt(k,154) + rxt(k,265) + rxt(k,268) &
                      + rxt(k,271) + rxt(k,282) + het_rates(k,52) )
         mat(k,670) = -( rxt(k,15) + rxt(k,16) + het_rates(k,53) )
         mat(k,44) = rxt(k,14)
         mat(k,610) = rxt(k,17) + .500_r8*rxt(k,266)
         mat(k,361) = rxt(k,19)
         mat(k,188) = rxt(k,298)
         mat(k,77) = rxt(k,311)
         mat(k,538) = 2.000_r8*rxt(k,147)*y(k,51)
         mat(k,608) = -( rxt(k,17) + rxt(k,266) + het_rates(k,54) )
         mat(k,259) = rxt(k,9)
         mat(k,68) = rxt(k,11) + rxt(k,153)
         mat(k,43) = rxt(k,13) + rxt(k,154)
         mat(k,359) = rxt(k,18)
         mat(k,90) = rxt(k,28)
         mat(k,250) = rxt(k,47)
         mat(k,350) = -( rxt(k,18) + rxt(k,19) + rxt(k,267) + het_rates(k,55) )
         mat(k,65) = rxt(k,10)
         mat(k,41) = rxt(k,13) + rxt(k,14) + rxt(k,154)
         mat(k,88) = rxt(k,29)
         mat(k,246) = rxt(k,46)
         mat(k,420) = -( rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76) + het_rates(k,56) )
         mat(k,733) = rxt(k,2)
         mat(k,504) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,79) + rxt(k,81) &
                      + 2.000_r8*rxt(k,82) + 2.000_r8*rxt(k,83) + rxt(k,84) + rxt(k,85) &
                      + rxt(k,86)
         mat(k,573) = rxt(k,8)
         mat(k,42) = rxt(k,14)
         mat(k,661) = rxt(k,15)
         mat(k,601) = rxt(k,17)
         mat(k,352) = rxt(k,18)
         mat(k,708) = rxt(k,27)
         mat(k,448) = rxt(k,45)
         mat(k,48) = rxt(k,59)
         mat(k,270) = rxt(k,89)
         mat(k,238) = rxt(k,90)
         mat(k,31) = rxt(k,91)
         mat(k,529) = rxt(k,96)
         mat(k,82) = rxt(k,305)
         mat(k,76) = rxt(k,310)
         mat(k,507) = -( rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,78) + rxt(k,79) &
                      + rxt(k,80) + rxt(k,81) + rxt(k,82) + rxt(k,83) + rxt(k,84) &
                      + rxt(k,85) + rxt(k,86) + het_rates(k,57) )
         mat(k,576) = rxt(k,8)
         mat(k,355) = rxt(k,19)
         mat(k,12) = rxt(k,92) + rxt(k,100)
         mat(k,15) = rxt(k,93)
         mat(k,532) = rxt(k,148)*y(k,51)
         mat(k,579) = -( rxt(k,7) + rxt(k,8) + het_rates(k,58) )
         mat(k,45) = -( rxt(k,59) + het_rates(k,59) )
         mat(k,50) = -( rxt(k,88) + het_rates(k,60) )
         mat(k,112) = -( het_rates(k,61) )
         mat(k,52) = rxt(k,88)
         mat(k,264) = rxt(k,89)
         mat(k,266) = -( rxt(k,89) + het_rates(k,63) )
         mat(k,236) = rxt(k,90)
         mat(k,235) = -( rxt(k,90) + het_rates(k,64) )
         mat(k,30) = rxt(k,91)
         mat(k,29) = -( rxt(k,91) + het_rates(k,65) )
         mat(k,9) = rxt(k,87)
         mat(k,1) = -( het_rates(k,66) )
         mat(k,196) = -( het_rates(k,67) )
         mat(k,653) = rxt(k,16)
         mat(k,223) = rxt(k,62)
         mat(k,407) = rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + rxt(k,76)
         mat(k,496) = rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,84) + rxt(k,85) + rxt(k,86)
         mat(k,376) = -( rxt(k,264) + het_rates(k,68) )
         mat(k,66) = rxt(k,11) + rxt(k,153)
         mat(k,472) = rxt(k,220)*y(k,16) + rxt(k,223)*y(k,18)
         mat(k,637) = rxt(k,221)*y(k,16) + rxt(k,224)*y(k,18)
         mat(k,528) = rxt(k,244)*y(k,21)
         mat(k,210) = -( het_rates(k,69) )
         mat(k,103) = -( het_rates(k,70) )
         mat(k,79) = rxt(k,305)
         mat(k,73) = rxt(k,310)
         mat(k,119) = -( het_rates(k,71) )
         mat(k,651) = rxt(k,16)
         mat(k,180) = rxt(k,298)
         mat(k,168) = rxt(k,301)
         mat(k,159) = -( het_rates(k,72) )
         mat(k,220) = rxt(k,62)
         mat(k,74) = rxt(k,311)
         mat(k,533) = -( rxt(k,96) + rxt(k,147)*y(k,51) + rxt(k,148)*y(k,51) &
                      + rxt(k,180)*y(k,6) + rxt(k,181)*y(k,7) + rxt(k,182)*y(k,9) &
                      + rxt(k,183)*y(k,10) + rxt(k,184)*y(k,11) + rxt(k,185)*y(k,12) &
                      + rxt(k,186)*y(k,13) + rxt(k,207)*y(k,8) + rxt(k,208)*y(k,22) &
                      + rxt(k,209)*y(k,37) + rxt(k,230)*y(k,14) + rxt(k,231)*y(k,16) &
                      + rxt(k,232)*y(k,41) + rxt(k,233)*y(k,42) + rxt(k,234)*y(k,43) &
                      + rxt(k,243)*y(k,21) + rxt(k,244)*y(k,21) + rxt(k,245)*y(k,21) &
                 + het_rates(k,73) )
         mat(k,737) = rxt(k,1)
         mat(k,508) = rxt(k,6)
         mat(k,577) = rxt(k,7)
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
         mat(k,11) = -( rxt(k,92) + rxt(k,100) + het_rates(k,74) )
         mat(k,564) = rxt(k,7)
         mat(k,13) = rxt(k,104) + rxt(k,103)*y(k,30)
         mat(k,14) = -( rxt(k,93) + rxt(k,104) + rxt(k,103)*y(k,30) + het_rates(k,75) )
         mat(k,181) = -( rxt(k,298) + het_rates(k,76) )
         mat(k,495) = rxt(k,78) + rxt(k,80)
         mat(k,171) = rxt(k,300)*y(k,30)
         mat(k,646) = -( rxt(k,219)*y(k,14) + rxt(k,221)*y(k,16) + rxt(k,222)*y(k,17) &
                      + rxt(k,224)*y(k,18) + rxt(k,226)*y(k,22) + rxt(k,227)*y(k,41) &
                      + rxt(k,228)*y(k,42) + rxt(k,229)*y(k,43) + rxt(k,241)*y(k,21) &
                 + het_rates(k,77) )
         mat(k,741) = rxt(k,3)
         mat(k,101) = 2.000_r8*rxt(k,4)
         mat(k,260) = rxt(k,9)
         mat(k,69) = rxt(k,10)
         mat(k,62) = rxt(k,22)
         mat(k,139) = rxt(k,57)
         mat(k,148) = rxt(k,58)
         mat(k,609) = .500_r8*rxt(k,266)
         mat(k,537) = rxt(k,243)*y(k,21)
         mat(k,170) = -( rxt(k,301) + rxt(k,300)*y(k,30) + het_rates(k,78) )
         mat(k,405) = rxt(k,73) + rxt(k,74)
         mat(k,494) = rxt(k,79) + rxt(k,81)
         mat(k,75) = rxt(k,285)
         mat(k,80) = rxt(k,286)
         mat(k,78) = -( rxt(k,286) + rxt(k,305) + het_rates(k,79) )
         mat(k,394) = rxt(k,75) + rxt(k,76)
         mat(k,489) = rxt(k,85) + rxt(k,86)
         mat(k,72) = rxt(k,287)
         mat(k,71) = -( rxt(k,285) + rxt(k,287) + rxt(k,310) + rxt(k,311) &
                 + het_rates(k,80) )
         mat(k,393) = rxt(k,71) + rxt(k,72)
         mat(k,488) = rxt(k,77) + rxt(k,84)
         mat(k,745) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,81) )
         mat(k,10) = rxt(k,87)
         mat(k,389) = rxt(k,264)
         mat(k,650) = rxt(k,219)*y(k,14) + rxt(k,221)*y(k,16) + rxt(k,222)*y(k,17) &
                      + rxt(k,224)*y(k,18) + rxt(k,229)*y(k,43) + rxt(k,241)*y(k,21)
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
