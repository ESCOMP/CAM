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
         prod(:,1) =rxt(:,386)*y(:,111)*y(:,95) +rxt(:,389)*y(:,96)
         prod(:,2) = (.500_r8*rxt(:,339)*y(:,87) +.600_r8*rxt(:,359)*y(:,112)) &
                 *y(:,111) +rxt(:,358)*y(:,112)*y(:,99)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,1) = 0._r8
         prod(:,2) = + extfrc(:,4)
         prod(:,29) = 0._r8
         prod(:,131) = 0._r8
         prod(:,49) = 0._r8
         prod(:,133) = 0._r8
         prod(:,82) = 0._r8
         prod(:,3) = 0._r8
         prod(:,74) = 0._r8
         prod(:,98) = 0._r8
         prod(:,54) = 0._r8
         prod(:,62) = 0._r8
         prod(:,60) = 0._r8
         prod(:,116) = 0._r8
         prod(:,102) = 0._r8
         prod(:,68) = 0._r8
         prod(:,34) = 0._r8
         prod(:,31) = 0._r8
         prod(:,40) = 0._r8
         prod(:,41) = 0._r8
         prod(:,35) = 0._r8
         prod(:,42) = 0._r8
         prod(:,36) = 0._r8
         prod(:,43) = 0._r8
         prod(:,37) = 0._r8
         prod(:,75) = 0._r8
         prod(:,136) = 0._r8
         prod(:,86) = 0._r8
         prod(:,33) = 0._r8
         prod(:,117) = 0._r8
         prod(:,66) = 0._r8
         prod(:,126) = 0._r8
         prod(:,95) = 0._r8
         prod(:,121) = 0._r8
         prod(:,89) = 0._r8
         prod(:,85) = 0._r8
         prod(:,140) = 0._r8
         prod(:,80) = 0._r8
         prod(:,69) = 0._r8
         prod(:,147) = 0._r8
         prod(:,71) = 0._r8
         prod(:,143) = 0._r8
         prod(:,44) = 0._r8
         prod(:,28) = 0._r8
         prod(:,134) = 0._r8
         prod(:,111) = 0._r8
         prod(:,4) = 0._r8
         prod(:,118) = + extfrc(:,1)
         prod(:,127) = 0._r8
         prod(:,48) = 0._r8
         prod(:,51) = 0._r8
         prod(:,56) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,72) = 0._r8
         prod(:,93) = 0._r8
         prod(:,38) = 0._r8
         prod(:,104) = 0._r8
         prod(:,114) = 0._r8
         prod(:,91) = 0._r8
         prod(:,148) = 0._r8
         prod(:,128) = 0._r8
         prod(:,32) = 0._r8
         prod(:,84) = 0._r8
         prod(:,39) = 0._r8
         prod(:,99) = 0._r8
         prod(:,50) = 0._r8
         prod(:,52) = 0._r8
         prod(:,57) = 0._r8
         prod(:,130) = 0._r8
         prod(:,58) = 0._r8
         prod(:,139) = 0._r8
         prod(:,73) = 0._r8
         prod(:,94) = 0._r8
         prod(:,101) = 0._r8
         prod(:,115) = 0._r8
         prod(:,64) = 0._r8
         prod(:,103) = 0._r8
         prod(:,97) = 0._r8
         prod(:,122) = 0._r8
         prod(:,78) = 0._r8
         prod(:,119) = 0._r8
         prod(:,123) = 0._r8
         prod(:,61) = 0._r8
         prod(:,125) = 0._r8
         prod(:,83) = 0._r8
         prod(:,124) = 0._r8
         prod(:,110) = (.800_r8*rxt(:,85) +.800_r8*rxt(:,86) +rxt(:,89) +rxt(:,90)) &
                  + extfrc(:,19)
         prod(:,53) = 0._r8
         prod(:,59) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,30) = 0._r8
         prod(:,11) = 0._r8
         prod(:,146) = + extfrc(:,7)
         prod(:,135) = + extfrc(:,6)
         prod(:,141) = 0._r8
         prod(:,67) = 0._r8
         prod(:,12) = + extfrc(:,8)
         prod(:,13) = + extfrc(:,9)
         prod(:,14) = 0._r8
         prod(:,15) = + extfrc(:,13)
         prod(:,16) = + extfrc(:,12)
         prod(:,145) = 0._r8
         prod(:,132) = 0._r8
         prod(:,138) = 0._r8
         prod(:,63) = 0._r8
         prod(:,65) = 0._r8
         prod(:,144) = + extfrc(:,21)
         prod(:,112) = 0._r8
         prod(:,79) = 0._r8
         prod(:,96) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = + extfrc(:,2)
         prod(:,81) = 0._r8
         prod(:,113) = 0._r8
         prod(:,70) = 0._r8
         prod(:,90) = 0._r8
         prod(:,19) = 0._r8
         prod(:,129) = 0._r8
         prod(:,109) = + extfrc(:,5)
         prod(:,55) = 0._r8
         prod(:,20) = + extfrc(:,10)
         prod(:,21) = + extfrc(:,11)
         prod(:,22) = 0._r8
         prod(:,23) = + extfrc(:,3)
         prod(:,24) = 0._r8
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,27) = 0._r8
         prod(:,88) = 0._r8
         prod(:,120) = 0._r8
         prod(:,47) = 0._r8
         prod(:,107) = (rxt(:,84) +rxt(:,87) +rxt(:,88) +rxt(:,89) +rxt(:,90) + &
                 rxt(:,91)) + extfrc(:,20)
         prod(:,137) = 0._r8
         prod(:,108) = (1.200_r8*rxt(:,85) +1.200_r8*rxt(:,86) +rxt(:,87) +rxt(:,91)) &
                  + extfrc(:,17)
         prod(:,87) = (rxt(:,84) +rxt(:,88)) + extfrc(:,15)
         prod(:,92) = 0._r8
         prod(:,100) = (rxt(:,87) +rxt(:,89) +rxt(:,90) +rxt(:,91)) + extfrc(:,16)
         prod(:,142) = 0._r8
         prod(:,45) = 0._r8
         prod(:,46) = 0._r8
         prod(:,106) = + extfrc(:,14)
         prod(:,105) = + extfrc(:,18)
         prod(:,77) = 0._r8
         prod(:,76) = 0._r8
         prod(:,149) = 0._r8
      end if
      end subroutine indprd
      end module mo_indprd
