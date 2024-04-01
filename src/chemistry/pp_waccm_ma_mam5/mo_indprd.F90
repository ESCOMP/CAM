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
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,1) = + extfrc(:,16)
         prod(:,2) = + extfrc(:,3)
         prod(:,97) = 0._r8
         prod(:,45) = 0._r8
         prod(:,101) = 0._r8
         prod(:,62) = 0._r8
         prod(:,3) = 0._r8
         prod(:,27) = 0._r8
         prod(:,34) = 0._r8
         prod(:,35) = 0._r8
         prod(:,29) = 0._r8
         prod(:,36) = 0._r8
         prod(:,30) = 0._r8
         prod(:,37) = 0._r8
         prod(:,31) = 0._r8
         prod(:,58) = 0._r8
         prod(:,86) = 0._r8
         prod(:,63) = 0._r8
         prod(:,32) = 0._r8
         prod(:,54) = 0._r8
         prod(:,82) = 0._r8
         prod(:,55) = 0._r8
         prod(:,81) = 0._r8
         prod(:,56) = 0._r8
         prod(:,96) = 0._r8
         prod(:,38) = 0._r8
         prod(:,26) = 0._r8
         prod(:,91) = 0._r8
         prod(:,75) = 0._r8
         prod(:,4) = 0._r8
         prod(:,72) = + extfrc(:,12)
         prod(:,61) = 0._r8
         prod(:,41) = 0._r8
         prod(:,43) = 0._r8
         prod(:,51) = + extfrc(:,2)
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,73) = 0._r8
         prod(:,84) = 0._r8
         prod(:,99) = 0._r8
         prod(:,28) = 0._r8
         prod(:,65) = 0._r8
         prod(:,33) = 0._r8
         prod(:,69) = 0._r8
         prod(:,42) = 0._r8
         prod(:,44) = 0._r8
         prod(:,48) = 0._r8
         prod(:,85) = 0._r8
         prod(:,49) = 0._r8
         prod(:,100) = 0._r8
         prod(:,57) = 0._r8
         prod(:,68) = 0._r8
         prod(:,70) = 0._r8
         prod(:,78) = (rxt(:,64) +.800_r8*rxt(:,66) +.800_r8*rxt(:,68) +rxt(:,70)) &
                  + extfrc(:,17)
         prod(:,46) = 0._r8
         prod(:,50) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,95) = + extfrc(:,13)
         prod(:,92) = + extfrc(:,14)
         prod(:,90) = 0._r8
         prod(:,11) = + extfrc(:,4)
         prod(:,12) = + extfrc(:,5)
         prod(:,13) = 0._r8
         prod(:,14) = + extfrc(:,6)
         prod(:,15) = + extfrc(:,7)
         prod(:,93) = 0._r8
         prod(:,87) = 0._r8
         prod(:,88) = 0._r8
         prod(:,52) = 0._r8
         prod(:,53) = 0._r8
         prod(:,16) = + extfrc(:,8)
         prod(:,17) = + extfrc(:,9)
         prod(:,66) = 0._r8
         prod(:,18) = 0._r8
         prod(:,83) = 0._r8
         prod(:,74) = + extfrc(:,15)
         prod(:,47) = 0._r8
         prod(:,19) = + extfrc(:,10)
         prod(:,20) = + extfrc(:,1)
         prod(:,21) = 0._r8
         prod(:,22) = + extfrc(:,11)
         prod(:,23) = 0._r8
         prod(:,24) = 0._r8
         prod(:,25) = 0._r8
         prod(:,76) = (rxt(:,63) +rxt(:,64) +rxt(:,65) +rxt(:,67) +rxt(:,69) + &
                 rxt(:,70)) + extfrc(:,21)
         prod(:,89) = 0._r8
         prod(:,77) = (rxt(:,65) +1.200_r8*rxt(:,66) +1.200_r8*rxt(:,68) +rxt(:,69)) &
                  + extfrc(:,18)
         prod(:,64) = (rxt(:,63) +rxt(:,67)) + extfrc(:,19)
         prod(:,67) = 0._r8
         prod(:,71) = (rxt(:,64) +rxt(:,65) +rxt(:,69) +rxt(:,70)) + extfrc(:,22)
         prod(:,94) = 0._r8
         prod(:,39) = 0._r8
         prod(:,40) = 0._r8
         prod(:,79) = + extfrc(:,23)
         prod(:,98) = + extfrc(:,24)
         prod(:,80) = + extfrc(:,20)
         prod(:,60) = 0._r8
         prod(:,59) = 0._r8
         prod(:,102) = 0._r8
      end if
      end subroutine indprd
      end module mo_indprd
