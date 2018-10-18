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
         prod(:,17) = (rxt(:,62) +rxt(:,113)*y(:,85) +rxt(:,114)*y(:,85) + &
                 rxt(:,115)*y(:,26) +rxt(:,116)*y(:,38) +rxt(:,123)*y(:,49) + &
                 rxt(:,124)*y(:,67) +rxt(:,125)*y(:,68) +rxt(:,167)*y(:,103) + &
                 rxt(:,169)*y(:,101) +rxt(:,185)*y(:,99) +rxt(:,203)*y(:,118) + &
                 rxt(:,220)*y(:,115) +rxt(:,238)*y(:,114) +rxt(:,255)*y(:,125) + &
                 rxt(:,257)*y(:,101) +rxt(:,264)*y(:,103) +rxt(:,279)*y(:,60) + &
                 rxt(:,280)*y(:,61))*y(:,90) + (rxt(:,119)*y(:,61) + &
                 rxt(:,120)*y(:,61) +rxt(:,121)*y(:,60) +rxt(:,122)*y(:,60) + &
                 rxt(:,154)*y(:,125) +rxt(:,157)*y(:,101) +rxt(:,177)*y(:,103) + &
                 rxt(:,195)*y(:,99) +rxt(:,212)*y(:,118) +rxt(:,230)*y(:,115) + &
                 rxt(:,248)*y(:,114) +rxt(:,259)*y(:,101) +rxt(:,260)*y(:,103)) &
                 *y(:,92) + (rxt(:,64) +rxt(:,126)*y(:,85) +rxt(:,127)*y(:,26) + &
                 rxt(:,129)*y(:,47) +rxt(:,131)*y(:,69) +rxt(:,150)*y(:,125) + &
                 rxt(:,173)*y(:,103) +rxt(:,190)*y(:,99) +rxt(:,208)*y(:,118) + &
                 rxt(:,224)*y(:,101) +rxt(:,226)*y(:,115) +rxt(:,243)*y(:,114)) &
                 *y(:,93) + (rxt(:,152)*y(:,125) +rxt(:,175)*y(:,103) + &
                 rxt(:,193)*y(:,99) +rxt(:,210)*y(:,118) +rxt(:,228)*y(:,115) + &
                 rxt(:,245)*y(:,114) +rxt(:,246)*y(:,101) +rxt(:,258)*y(:,103) + &
                 rxt(:,270)*y(:,101))*y(:,91) + (rxt(:,148)*y(:,125) + &
                 rxt(:,171)*y(:,103) +rxt(:,188)*y(:,99) +rxt(:,202)*y(:,101) + &
                 rxt(:,206)*y(:,118) +rxt(:,223)*y(:,115) +rxt(:,241)*y(:,114)) &
                 *y(:,96) + (rxt(:,368) +rxt(:,305)*y(:,94) +rxt(:,306)*y(:,134)) &
                 *y(:,117) + (rxt(:,536)*y(:,130) +rxt(:,540)*y(:,130))*y(:,29)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) =rxt(:,426)*y(:,61)*y(:,54)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,1) = + extfrc(:,5)
         prod(:,2) = + extfrc(:,6)
         prod(:,27) = 0._r8
         prod(:,80) = 0._r8
         prod(:,45) = 0._r8
         prod(:,81) =.180_r8*rxt(:,25)*y(:,22)
         prod(:,64) =rxt(:,40)*y(:,17) +rxt(:,42)*y(:,19) +rxt(:,24)*y(:,22)
         prod(:,35) = 0._r8
         prod(:,24) = 0._r8
         prod(:,21) = 0._r8
         prod(:,95) = 0._r8
         prod(:,60) = 0._r8
         prod(:,44) = (rxt(:,26) +rxt(:,61))*y(:,30) +.380_r8*rxt(:,25)*y(:,22) &
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
         prod(:,54) =rxt(:,33)*y(:,9) +rxt(:,37)*y(:,13)
         prod(:,98) = (rxt(:,24) +.330_r8*rxt(:,25))*y(:,22)
         prod(:,75) =1.440_r8*rxt(:,25)*y(:,22)
         prod(:,46) = 0._r8
         prod(:,22) = 0._r8
         prod(:,56) = 0._r8
         prod(:,100) = 0._r8
         prod(:,26) = 0._r8
         prod(:,101) = 0._r8
         prod(:,39) = 0._r8
         prod(:,55) = 0._r8
         prod(:,57) = 0._r8
         prod(:,50) = 0._r8
         prod(:,63) = (.800_r8*rxt(:,66) +.800_r8*rxt(:,68) +rxt(:,72) +rxt(:,73)) &
                  + extfrc(:,15)
         prod(:,62) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,106) = + extfrc(:,14)
         prod(:,104) = + extfrc(:,3)
         prod(:,97) = 0._r8
         prod(:,9) = + extfrc(:,7)
         prod(:,10) = + extfrc(:,8)
         prod(:,11) = 0._r8
         prod(:,12) = + extfrc(:,9)
         prod(:,107) = (rxt(:,26) +rxt(:,61))*y(:,30) +.180_r8*rxt(:,25)*y(:,22) &
                  + extfrc(:,22)
         prod(:,99) = 0._r8
         prod(:,91) = 0._r8
         prod(:,33) = 0._r8
         prod(:,34) = 0._r8
         prod(:,13) = + extfrc(:,10)
         prod(:,14) = + extfrc(:,11)
         prod(:,49) = 0._r8
         prod(:,68) = 0._r8
         prod(:,59) = + extfrc(:,4)
         prod(:,29) = 0._r8
         prod(:,15) = + extfrc(:,12)
         prod(:,16) = + extfrc(:,1)
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,74) =rxt(:,32)*y(:,8) +rxt(:,33)*y(:,9) +2.000_r8*rxt(:,39)*y(:,15) &
                  +rxt(:,40)*y(:,17) +3.000_r8*rxt(:,43)*y(:,23) +2.000_r8*rxt(:,51) &
                 *y(:,40)
         prod(:,90) =4.000_r8*rxt(:,31)*y(:,7) +rxt(:,32)*y(:,8) +2.000_r8*rxt(:,34) &
                 *y(:,10) +2.000_r8*rxt(:,35)*y(:,11) +2.000_r8*rxt(:,36)*y(:,12) &
                  +rxt(:,37)*y(:,13) +2.000_r8*rxt(:,38)*y(:,14) +3.000_r8*rxt(:,41) &
                 *y(:,18) +rxt(:,42)*y(:,19) +rxt(:,53)*y(:,44) +rxt(:,54)*y(:,45) &
                  +rxt(:,55)*y(:,46)
         prod(:,84) = 0._r8
         prod(:,70) = 0._r8
         prod(:,69) = 0._r8
         prod(:,61) = 0._r8
         prod(:,88) = 0._r8
         prod(:,67) = 0._r8
         prod(:,78) = 0._r8
         prod(:,83) = 0._r8
         prod(:,89) = (rxt(:,67) +rxt(:,69) +rxt(:,70) +rxt(:,71) +rxt(:,72) + &
                 rxt(:,73)) + extfrc(:,20)
         prod(:,36) = 0._r8
         prod(:,65) = 0._r8
         prod(:,86) = 0._r8
         prod(:,47) = 0._r8
         prod(:,111) = 0._r8
         prod(:,30) = 0._r8
         prod(:,105) = 0._r8
         prod(:,31) = 0._r8
         prod(:,94) = 0._r8
         prod(:,51) = 0._r8
         prod(:,43) = (1.200_r8*rxt(:,66) +1.200_r8*rxt(:,68) +rxt(:,69) +rxt(:,70)) &
                  + extfrc(:,16)
         prod(:,48) = (rxt(:,67) +rxt(:,71)) + extfrc(:,17)
         prod(:,108) = 0._r8
         prod(:,71) = 0._r8
         prod(:,85) = 0._r8
         prod(:,77) = 0._r8
         prod(:,79) = 0._r8
         prod(:,73) = 0._r8
         prod(:,76) = 0._r8
         prod(:,92) = 0._r8
         prod(:,93) = 0._r8
         prod(:,37) = 0._r8
         prod(:,41) = 0._r8
         prod(:,96) = 0._r8
         prod(:,40) = 0._r8
         prod(:,42) = (rxt(:,69) +rxt(:,70) +rxt(:,72) +rxt(:,73)) + extfrc(:,21)
         prod(:,82) =rxt(:,13)*y(:,55)
         prod(:,58) = 0._r8
         prod(:,28) = 0._r8
         prod(:,102) = 0._r8
         prod(:,103) = + extfrc(:,23)
         prod(:,52) = 0._r8
         prod(:,72) = 0._r8
         prod(:,38) = 0._r8
         prod(:,66) = 0._r8
         prod(:,87) =.330_r8*rxt(:,25)*y(:,22) + extfrc(:,18)
         prod(:,109) = 0._r8
         prod(:,110) = 0._r8
         prod(:,53) = + extfrc(:,19)
         prod(:,112) =.050_r8*rxt(:,25)*y(:,22)
      end if
      end subroutine indprd
      end module mo_indprd
