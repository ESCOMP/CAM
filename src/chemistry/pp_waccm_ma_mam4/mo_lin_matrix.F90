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
         mat(k,44) = -( rxt(k,26) + het_rates(k,3) )
         mat(k,395) = -( rxt(k,27) + het_rates(k,4) )
         mat(k,111) = rxt(k,28)
         mat(k,106) = -( rxt(k,28) + rxt(k,29) + rxt(k,269) + rxt(k,272) + rxt(k,277) &
                 + het_rates(k,5) )
         mat(k,502) = -( rxt(k,20) + rxt(k,21) + het_rates(k,16) )
         mat(k,80) = rxt(k,22)
         mat(k,549) = rxt(k,243)*y(k,22) + rxt(k,244)*y(k,22)
         mat(k,320) = -( het_rates(k,20) )
         mat(k,366) = rxt(k,155)*y(k,22)
         mat(k,162) = rxt(k,211)*y(k,22)
         mat(k,652) = rxt(k,240)*y(k,22)
         mat(k,542) = rxt(k,242)*y(k,22)
         mat(k,77) = -( rxt(k,22) + het_rates(k,21) )
         mat(k,35) = -( rxt(k,43) + het_rates(k,24) )
         mat(k,21) = -( rxt(k,44) + rxt(k,189) + het_rates(k,25) )
         mat(k,479) = -( rxt(k,45) + het_rates(k,26) )
         mat(k,242) = rxt(k,47)
         mat(k,67) = rxt(k,59)
         mat(k,23) = 2.000_r8*rxt(k,189)
         mat(k,236) = -( rxt(k,46) + rxt(k,47) + rxt(k,271) + rxt(k,276) + rxt(k,282) &
                 + het_rates(k,27) )
         mat(k,102) = -( het_rates(k,29) )
         mat(k,489) = rxt(k,20) + rxt(k,21)
         mat(k,70) = rxt(k,87)
         mat(k,358) = rxt(k,222)*y(k,19)
         mat(k,139) = rxt(k,296)*y(k,30)
         mat(k,27) = -( rxt(k,48) + het_rates(k,31) )
         mat(k,533) = rxt(k,180)*y(k,8) + rxt(k,182)*y(k,11) + 2.000_r8*rxt(k,183)*y(k,12) &
                      + 2.000_r8*rxt(k,184)*y(k,13) + rxt(k,185)*y(k,14) &
                      + rxt(k,206)*y(k,9) + 2.000_r8*rxt(k,208)*y(k,40) &
                      + rxt(k,232)*y(k,45) + rxt(k,233)*y(k,46)
         mat(k,631) = rxt(k,227)*y(k,45) + rxt(k,228)*y(k,46)
         mat(k,37) = -( rxt(k,49) + het_rates(k,32) )
         mat(k,535) = rxt(k,181)*y(k,10) + rxt(k,182)*y(k,11) + rxt(k,231)*y(k,44)
         mat(k,632) = rxt(k,226)*y(k,44)
         mat(k,58) = -( het_rates(k,33) )
         mat(k,3) = -( het_rates(k,34) )
         mat(k,4) = -( het_rates(k,35) )
         mat(k,5) = -( het_rates(k,36) )
         mat(k,159) = -( rxt(k,211)*y(k,22) + het_rates(k,37) )
         mat(k,28) = 2.000_r8*rxt(k,48)
         mat(k,38) = rxt(k,49)
         mat(k,42) = rxt(k,56)
         mat(k,536) = rxt(k,184)*y(k,13) + rxt(k,206)*y(k,9)
         mat(k,273) = -( het_rates(k,38) )
         mat(k,716) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,492) = 2.000_r8*rxt(k,20)
         mat(k,78) = rxt(k,22)
         mat(k,212) = rxt(k,51)
         mat(k,287) = rxt(k,55)
         mat(k,43) = rxt(k,56)
         mat(k,539) = rxt(k,243)*y(k,22)
         mat(k,574) = -( het_rates(k,39) )
         mat(k,728) = rxt(k,1)
         mat(k,505) = rxt(k,21)
         mat(k,552) = rxt(k,244)*y(k,22)
         mat(k,114) = -( rxt(k,4) + het_rates(k,41) )
         mat(k,410) = .500_r8*rxt(k,263)
         mat(k,24) = -( rxt(k,86) + het_rates(k,42) )
         mat(k,211) = -( rxt(k,51) + het_rates(k,43) )
         mat(k,288) = -( rxt(k,55) + het_rates(k,47) )
         mat(k,364) = rxt(k,155)*y(k,22) + rxt(k,217)*y(k,15) + rxt(k,219)*y(k,17) &
                      + 2.000_r8*rxt(k,222)*y(k,19) + rxt(k,224)*y(k,23)
         mat(k,41) = -( rxt(k,56) + het_rates(k,48) )
         mat(k,158) = rxt(k,211)*y(k,22)
         mat(k,248) = -( rxt(k,9) + het_rates(k,49) )
         mat(k,53) = 2.000_r8*rxt(k,264) + 2.000_r8*rxt(k,267) + 2.000_r8*rxt(k,270) &
                      + 2.000_r8*rxt(k,281)
         mat(k,440) = .500_r8*rxt(k,265)
         mat(k,337) = rxt(k,266)
         mat(k,108) = rxt(k,269) + rxt(k,272) + rxt(k,277)
         mat(k,237) = rxt(k,271) + rxt(k,276) + rxt(k,282)
         mat(k,83) = -( rxt(k,10) + rxt(k,11) + rxt(k,152) + het_rates(k,50) )
         mat(k,150) = -( rxt(k,57) + het_rates(k,51) )
         mat(k,107) = rxt(k,269) + rxt(k,272) + rxt(k,277)
         mat(k,220) = -( rxt(k,58) + het_rates(k,52) )
         mat(k,235) = rxt(k,271) + rxt(k,276) + rxt(k,282)
         mat(k,202) = -( rxt(k,61) + het_rates(k,53) )
         mat(k,514) = rxt(k,15)
         mat(k,143) = rxt(k,297)
         mat(k,52) = -( rxt(k,13) + rxt(k,14) + rxt(k,153) + rxt(k,264) + rxt(k,267) &
                      + rxt(k,270) + rxt(k,281) + het_rates(k,55) )
         mat(k,6) = -( het_rates(k,56) )
         mat(k,7) = -( het_rates(k,57) )
         mat(k,8) = -( het_rates(k,58) )
         mat(k,525) = -( rxt(k,15) + rxt(k,16) + het_rates(k,59) )
         mat(k,56) = rxt(k,14)
         mat(k,453) = rxt(k,17) + .500_r8*rxt(k,265)
         mat(k,349) = rxt(k,19)
         mat(k,171) = rxt(k,294)
         mat(k,550) = 2.000_r8*rxt(k,146)*y(k,54)
         mat(k,450) = -( rxt(k,17) + rxt(k,265) + het_rates(k,60) )
         mat(k,252) = rxt(k,9)
         mat(k,86) = rxt(k,11) + rxt(k,152)
         mat(k,55) = rxt(k,13) + rxt(k,153)
         mat(k,346) = rxt(k,18)
         mat(k,112) = rxt(k,28)
         mat(k,241) = rxt(k,47)
         mat(k,342) = -( rxt(k,18) + rxt(k,19) + rxt(k,266) + het_rates(k,61) )
         mat(k,84) = rxt(k,10)
         mat(k,54) = rxt(k,13) + rxt(k,14) + rxt(k,153)
         mat(k,110) = rxt(k,29)
         mat(k,239) = rxt(k,46)
         mat(k,9) = -( het_rates(k,62) )
         mat(k,10) = -( het_rates(k,63) )
         mat(k,11) = -( het_rates(k,64) )
         mat(k,12) = -( het_rates(k,65) )
         mat(k,707) = -( rxt(k,70) + rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) &
                      + rxt(k,75) + het_rates(k,66) )
         mat(k,732) = rxt(k,2)
         mat(k,607) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,76) + rxt(k,78) + rxt(k,80) &
                      + 2.000_r8*rxt(k,81) + 2.000_r8*rxt(k,82) + rxt(k,83) + rxt(k,84) &
                      + rxt(k,85)
         mat(k,629) = rxt(k,8)
         mat(k,57) = rxt(k,14)
         mat(k,531) = rxt(k,15)
         mat(k,459) = rxt(k,17)
         mat(k,355) = rxt(k,18)
         mat(k,406) = rxt(k,27)
         mat(k,487) = rxt(k,45)
         mat(k,68) = rxt(k,59)
         mat(k,271) = rxt(k,88)
         mat(k,232) = rxt(k,89)
         mat(k,50) = rxt(k,90)
         mat(k,556) = rxt(k,95)
         mat(k,604) = -( rxt(k,5) + rxt(k,6) + rxt(k,76) + rxt(k,77) + rxt(k,78) &
                      + rxt(k,79) + rxt(k,80) + rxt(k,81) + rxt(k,82) + rxt(k,83) &
                      + rxt(k,84) + rxt(k,85) + het_rates(k,67) )
         mat(k,626) = rxt(k,8)
         mat(k,352) = rxt(k,19)
         mat(k,31) = rxt(k,91) + rxt(k,99)
         mat(k,34) = rxt(k,92)
         mat(k,553) = rxt(k,147)*y(k,54)
         mat(k,627) = -( rxt(k,7) + rxt(k,8) + het_rates(k,68) )
         mat(k,64) = -( rxt(k,59) + het_rates(k,69) )
         mat(k,69) = -( rxt(k,87) + het_rates(k,70) )
         mat(k,13) = -( het_rates(k,71) )
         mat(k,14) = -( het_rates(k,72) )
         mat(k,132) = -( het_rates(k,73) )
         mat(k,71) = rxt(k,87)
         mat(k,257) = rxt(k,88)
         mat(k,259) = -( rxt(k,88) + het_rates(k,74) )
         mat(k,229) = rxt(k,89)
         mat(k,228) = -( rxt(k,89) + het_rates(k,75) )
         mat(k,49) = rxt(k,90)
         mat(k,48) = -( rxt(k,90) + het_rates(k,76) )
         mat(k,25) = rxt(k,86)
         mat(k,15) = -( het_rates(k,77) )
         mat(k,16) = -( het_rates(k,78) )
         mat(k,17) = -( het_rates(k,79) )
         mat(k,18) = -( het_rates(k,80) )
         mat(k,19) = -( het_rates(k,81) )
         mat(k,20) = -( het_rates(k,82) )
         mat(k,307) = -( het_rates(k,83) )
         mat(k,45) = rxt(k,26)
         mat(k,391) = rxt(k,27)
         mat(k,109) = rxt(k,29)
         mat(k,213) = rxt(k,51)
         mat(k,152) = rxt(k,57)
         mat(k,541) = rxt(k,180)*y(k,8) + rxt(k,206)*y(k,9) + 3.000_r8*rxt(k,207)*y(k,23) &
                      + 2.000_r8*rxt(k,208)*y(k,40) + 2.000_r8*rxt(k,229)*y(k,15) &
                      + rxt(k,230)*y(k,17)
         mat(k,365) = 2.000_r8*rxt(k,217)*y(k,15) + rxt(k,219)*y(k,17) &
                      + 3.000_r8*rxt(k,224)*y(k,23)
         mat(k,651) = 2.000_r8*rxt(k,218)*y(k,15) + rxt(k,220)*y(k,17) &
                      + 3.000_r8*rxt(k,225)*y(k,23)
         mat(k,368) = -( rxt(k,155)*y(k,22) + rxt(k,217)*y(k,15) + rxt(k,219)*y(k,17) &
                      + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,23) + het_rates(k,84) )
         mat(k,46) = rxt(k,26)
         mat(k,36) = 2.000_r8*rxt(k,43)
         mat(k,22) = 2.000_r8*rxt(k,44)
         mat(k,475) = rxt(k,45)
         mat(k,240) = rxt(k,46)
         mat(k,39) = rxt(k,49)
         mat(k,292) = rxt(k,55)
         mat(k,222) = rxt(k,58)
         mat(k,544) = 4.000_r8*rxt(k,179)*y(k,7) + rxt(k,180)*y(k,8) &
                      + 2.000_r8*rxt(k,181)*y(k,10) + 2.000_r8*rxt(k,182)*y(k,11) &
                      + 2.000_r8*rxt(k,183)*y(k,12) + rxt(k,184)*y(k,13) &
                      + 2.000_r8*rxt(k,185)*y(k,14) + rxt(k,231)*y(k,44) &
                      + rxt(k,232)*y(k,45) + rxt(k,233)*y(k,46)
         mat(k,654) = 3.000_r8*rxt(k,221)*y(k,18) + rxt(k,223)*y(k,19) &
                      + rxt(k,226)*y(k,44) + rxt(k,227)*y(k,45) + rxt(k,228)*y(k,46)
         mat(k,188) = -( het_rates(k,85) )
         mat(k,513) = rxt(k,16)
         mat(k,201) = rxt(k,61)
         mat(k,682) = rxt(k,70) + rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) &
                      + rxt(k,75)
         mat(k,589) = rxt(k,76) + rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) &
                      + rxt(k,83) + rxt(k,84) + rxt(k,85)
         mat(k,421) = -( rxt(k,263) + het_rates(k,86) )
         mat(k,85) = rxt(k,11) + rxt(k,152)
         mat(k,370) = rxt(k,219)*y(k,17) + rxt(k,222)*y(k,19)
         mat(k,656) = rxt(k,220)*y(k,17) + rxt(k,223)*y(k,19)
         mat(k,546) = rxt(k,243)*y(k,22)
         mat(k,97) = -( het_rates(k,87) )
         mat(k,122) = -( het_rates(k,88) )
         mat(k,176) = -( het_rates(k,89) )
         mat(k,512) = rxt(k,16)
         mat(k,168) = rxt(k,294)
         mat(k,142) = rxt(k,297)
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
         mat(k,90) = -( het_rates(k,90) )
         mat(k,197) = rxt(k,61)
         mat(k,551) = -( rxt(k,95) + rxt(k,146)*y(k,54) + rxt(k,147)*y(k,54) &
                      + rxt(k,179)*y(k,7) + rxt(k,180)*y(k,8) + rxt(k,181)*y(k,10) &
                      + rxt(k,182)*y(k,11) + rxt(k,183)*y(k,12) + rxt(k,184)*y(k,13) &
                      + rxt(k,185)*y(k,14) + rxt(k,206)*y(k,9) + rxt(k,207)*y(k,23) &
                      + rxt(k,208)*y(k,40) + rxt(k,229)*y(k,15) + rxt(k,230)*y(k,17) &
                      + rxt(k,231)*y(k,44) + rxt(k,232)*y(k,45) + rxt(k,233)*y(k,46) &
                      + rxt(k,242)*y(k,22) + rxt(k,243)*y(k,22) + rxt(k,244)*y(k,22) &
                 + het_rates(k,91) )
         mat(k,727) = rxt(k,1)
         mat(k,602) = rxt(k,6)
         mat(k,624) = rxt(k,7)
         mat(k,30) = -( rxt(k,91) + rxt(k,99) + het_rates(k,92) )
         mat(k,609) = rxt(k,7)
         mat(k,32) = rxt(k,103) + rxt(k,102)*y(k,30)
         mat(k,33) = -( rxt(k,92) + rxt(k,103) + rxt(k,102)*y(k,30) + het_rates(k,93) )
         mat(k,167) = -( rxt(k,294) + het_rates(k,94) )
         mat(k,587) = rxt(k,77) + rxt(k,79)
         mat(k,141) = rxt(k,296)*y(k,30)
         mat(k,665) = -( rxt(k,218)*y(k,15) + rxt(k,220)*y(k,17) + rxt(k,221)*y(k,18) &
                      + rxt(k,223)*y(k,19) + rxt(k,225)*y(k,23) + rxt(k,226)*y(k,44) &
                      + rxt(k,227)*y(k,45) + rxt(k,228)*y(k,46) + rxt(k,240)*y(k,22) &
                 + het_rates(k,95) )
         mat(k,731) = rxt(k,3)
         mat(k,118) = 2.000_r8*rxt(k,4)
         mat(k,254) = rxt(k,9)
         mat(k,88) = rxt(k,10)
         mat(k,81) = rxt(k,22)
         mat(k,155) = rxt(k,57)
         mat(k,224) = rxt(k,58)
         mat(k,458) = .500_r8*rxt(k,265)
         mat(k,555) = rxt(k,242)*y(k,22)
         mat(k,140) = -( rxt(k,297) + rxt(k,296)*y(k,30) + het_rates(k,96) )
         mat(k,678) = rxt(k,70) + rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) &
                      + rxt(k,75)
         mat(k,586) = rxt(k,76) + rxt(k,78) + rxt(k,80) + rxt(k,83) + rxt(k,84) &
                      + rxt(k,85)
         mat(k,733) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,97) )
         mat(k,26) = rxt(k,86)
         mat(k,667) = rxt(k,218)*y(k,15) + rxt(k,220)*y(k,17) + rxt(k,221)*y(k,18) &
                      + rxt(k,223)*y(k,19) + rxt(k,228)*y(k,46) + rxt(k,240)*y(k,22)
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
