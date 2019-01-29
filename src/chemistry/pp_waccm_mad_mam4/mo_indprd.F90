      module mo_indprd
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: indprd
      contains
      subroutine indprd( class, prod, nprod, y, extfrc, rxt, chnkpnts )
      use chem_mods, only : gas_pcnst, extcnt, rxntot
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: chnkpnts
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(in) :: extfrc(chnkpnts,extcnt)
      real(r8), intent(inout) :: prod(chnkpnts,nprod)
!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,1) = 0._r8
         prod(:,2) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,12) = 0._r8
         prod(:,13) = 0._r8
         prod(:,14) = 0._r8
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,17) = (rxt(:,63) +rxt(:,114)*y(:,86) +rxt(:,115)*y(:,86) + &
                 rxt(:,116)*y(:,26) +rxt(:,117)*y(:,38) +rxt(:,124)*y(:,49) + &
                 rxt(:,125)*y(:,67) +rxt(:,126)*y(:,68) +rxt(:,168)*y(:,104) + &
                 rxt(:,170)*y(:,102) +rxt(:,186)*y(:,100) +rxt(:,204)*y(:,119) + &
                 rxt(:,221)*y(:,116) +rxt(:,239)*y(:,115) +rxt(:,256)*y(:,126) + &
                 rxt(:,258)*y(:,102) +rxt(:,265)*y(:,104) +rxt(:,280)*y(:,60) + &
                 rxt(:,281)*y(:,61))*y(:,91) + (rxt(:,120)*y(:,61) + &
                 rxt(:,121)*y(:,61) +rxt(:,122)*y(:,60) +rxt(:,123)*y(:,60) + &
                 rxt(:,155)*y(:,126) +rxt(:,158)*y(:,102) +rxt(:,178)*y(:,104) + &
                 rxt(:,196)*y(:,100) +rxt(:,213)*y(:,119) +rxt(:,231)*y(:,116) + &
                 rxt(:,249)*y(:,115) +rxt(:,260)*y(:,102) +rxt(:,261)*y(:,104)) &
                 *y(:,93) + (rxt(:,65) +rxt(:,127)*y(:,86) +rxt(:,128)*y(:,26) + &
                 rxt(:,130)*y(:,47) +rxt(:,132)*y(:,69) +rxt(:,151)*y(:,126) + &
                 rxt(:,174)*y(:,104) +rxt(:,191)*y(:,100) +rxt(:,209)*y(:,119) + &
                 rxt(:,225)*y(:,102) +rxt(:,227)*y(:,116) +rxt(:,244)*y(:,115)) &
                 *y(:,94) + (rxt(:,153)*y(:,126) +rxt(:,176)*y(:,104) + &
                 rxt(:,194)*y(:,100) +rxt(:,211)*y(:,119) +rxt(:,229)*y(:,116) + &
                 rxt(:,246)*y(:,115) +rxt(:,247)*y(:,102) +rxt(:,259)*y(:,104) + &
                 rxt(:,271)*y(:,102))*y(:,92) + (rxt(:,149)*y(:,126) + &
                 rxt(:,172)*y(:,104) +rxt(:,189)*y(:,100) +rxt(:,203)*y(:,102) + &
                 rxt(:,207)*y(:,119) +rxt(:,224)*y(:,116) +rxt(:,242)*y(:,115)) &
                 *y(:,97) + (rxt(:,369) +rxt(:,306)*y(:,95) +rxt(:,307)*y(:,137)) &
                 *y(:,118) + (rxt(:,537)*y(:,131) +rxt(:,541)*y(:,131))*y(:,29)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) =rxt(:,427)*y(:,61)*y(:,54)
         prod(:,23) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,1) = + extfrc(:,5)
         prod(:,2) = + extfrc(:,6)
         prod(:,27) = 0._r8
         prod(:,82) = 0._r8
         prod(:,45) = 0._r8
         prod(:,83) =.180_r8*rxt(:,25)*y(:,22)
         prod(:,66) =rxt(:,40)*y(:,17) +rxt(:,42)*y(:,19) +rxt(:,24)*y(:,22)
         prod(:,35) = 0._r8
         prod(:,24) = 0._r8
         prod(:,21) = 0._r8
         prod(:,112) = 0._r8
         prod(:,62) = 0._r8
         prod(:,44) = (rxt(:,26) +rxt(:,62))*y(:,30) +.380_r8*rxt(:,25)*y(:,22) &
                  + extfrc(:,13)
         prod(:,23) =rxt(:,32)*y(:,8) +rxt(:,33)*y(:,9) +rxt(:,35)*y(:,11) &
                  +2.000_r8*rxt(:,36)*y(:,12) +2.000_r8*rxt(:,37)*y(:,13) +rxt(:,38) &
                 *y(:,14) +2.000_r8*rxt(:,51)*y(:,40) +rxt(:,54)*y(:,45) +rxt(:,55) &
                 *y(:,46)
         prod(:,25) =rxt(:,34)*y(:,10) +rxt(:,35)*y(:,11) +rxt(:,53)*y(:,44)
         prod(:,32) = + extfrc(:,2)
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,53) =rxt(:,33)*y(:,9) +rxt(:,37)*y(:,13)
         prod(:,113) = (rxt(:,24) +.330_r8*rxt(:,25))*y(:,22)
         prod(:,77) =1.440_r8*rxt(:,25)*y(:,22)
         prod(:,46) = 0._r8
         prod(:,22) = 0._r8
         prod(:,58) = 0._r8
         prod(:,98) = 0._r8
         prod(:,26) = 0._r8
         prod(:,101) = 0._r8
         prod(:,39) = 0._r8
         prod(:,54) = 0._r8
         prod(:,57) = 0._r8
         prod(:,50) = 0._r8
         prod(:,65) = (rxt(:,69) +rxt(:,70) +.800_r8*rxt(:,72) +.800_r8*rxt(:,73)) &
                  + extfrc(:,15)
         prod(:,64) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,106) = + extfrc(:,14)
         prod(:,94) = + extfrc(:,3)
         prod(:,97) = 0._r8
         prod(:,9) = + extfrc(:,7)
         prod(:,10) = + extfrc(:,8)
         prod(:,11) = 0._r8
         prod(:,12) = + extfrc(:,9)
         prod(:,104) = (rxt(:,26) +rxt(:,62))*y(:,30) +.180_r8*rxt(:,25)*y(:,22) &
                  + extfrc(:,22)
         prod(:,100) = 0._r8
         prod(:,99) = 0._r8
         prod(:,33) = 0._r8
         prod(:,34) = 0._r8
         prod(:,13) = + extfrc(:,10)
         prod(:,14) = + extfrc(:,11)
         prod(:,49) = 0._r8
         prod(:,70) = 0._r8
         prod(:,60) = + extfrc(:,4)
         prod(:,29) = 0._r8
         prod(:,15) = + extfrc(:,12)
         prod(:,16) = + extfrc(:,1)
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,75) =rxt(:,32)*y(:,8) +rxt(:,33)*y(:,9) +2.000_r8*rxt(:,39)*y(:,15) &
                  +rxt(:,40)*y(:,17) +3.000_r8*rxt(:,43)*y(:,23) +2.000_r8*rxt(:,51) &
                 *y(:,40)
         prod(:,105) =4.000_r8*rxt(:,31)*y(:,7) +rxt(:,32)*y(:,8) +2.000_r8*rxt(:,34) &
                 *y(:,10) +2.000_r8*rxt(:,35)*y(:,11) +2.000_r8*rxt(:,36)*y(:,12) &
                  +rxt(:,37)*y(:,13) +2.000_r8*rxt(:,38)*y(:,14) +3.000_r8*rxt(:,41) &
                 *y(:,18) +rxt(:,42)*y(:,19) +rxt(:,53)*y(:,44) +rxt(:,54)*y(:,45) &
                  +rxt(:,55)*y(:,46)
         prod(:,86) = 0._r8
         prod(:,72) = 0._r8
         prod(:,71) = 0._r8
         prod(:,63) = 0._r8
         prod(:,90) = 0._r8
         prod(:,69) = 0._r8
         prod(:,80) = 0._r8
         prod(:,85) = 0._r8
         prod(:,107) = (rxt(:,67) +rxt(:,68) +rxt(:,69) +rxt(:,70) +rxt(:,71) + &
                 rxt(:,74)) + extfrc(:,20)
         prod(:,36) = 0._r8
         prod(:,67) = 0._r8
         prod(:,88) = 0._r8
         prod(:,48) = 0._r8
         prod(:,95) = 0._r8
         prod(:,30) = 0._r8
         prod(:,111) = 0._r8
         prod(:,31) = 0._r8
         prod(:,108) = 0._r8
         prod(:,51) = 0._r8
         prod(:,55) = (rxt(:,68) +rxt(:,71) +1.200_r8*rxt(:,72) +1.200_r8*rxt(:,73)) &
                  + extfrc(:,16)
         prod(:,47) = (rxt(:,67) +rxt(:,74)) + extfrc(:,17)
         prod(:,91) = 0._r8
         prod(:,73) = 0._r8
         prod(:,87) = 0._r8
         prod(:,79) = 0._r8
         prod(:,81) = 0._r8
         prod(:,76) = 0._r8
         prod(:,78) = 0._r8
         prod(:,92) = 0._r8
         prod(:,93) = 0._r8
         prod(:,37) = 0._r8
         prod(:,41) = 0._r8
         prod(:,96) = 0._r8
         prod(:,40) = 0._r8
         prod(:,56) = (rxt(:,68) +rxt(:,69) +rxt(:,70) +rxt(:,71)) + extfrc(:,21)
         prod(:,84) =rxt(:,13)*y(:,55)
         prod(:,59) = 0._r8
         prod(:,28) = 0._r8
         prod(:,102) = 0._r8
         prod(:,103) = + extfrc(:,23)
         prod(:,52) = 0._r8
         prod(:,74) = 0._r8
         prod(:,38) = 0._r8
         prod(:,68) = 0._r8
         prod(:,89) =.330_r8*rxt(:,25)*y(:,22) + extfrc(:,18)
         prod(:,109) = 0._r8
         prod(:,110) = 0._r8
         prod(:,61) = + extfrc(:,19)
         prod(:,43) = 0._r8
         prod(:,42) = 0._r8
         prod(:,114) =.050_r8*rxt(:,25)*y(:,22)
      end if
      end subroutine indprd
      end module mo_indprd
