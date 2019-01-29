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
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) =rxt(:,422)*y(:,52)*y(:,48)
         prod(:,22) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,6) = 0._r8
         prod(:,52) = 0._r8
         prod(:,17) = 0._r8
         prod(:,56) =.180_r8*rxt(:,25)*y(:,20)
         prod(:,39) =rxt(:,40)*y(:,15) +rxt(:,42)*y(:,17) +rxt(:,24)*y(:,20)
         prod(:,12) = 0._r8
         prod(:,4) = 0._r8
         prod(:,1) = 0._r8
         prod(:,81) = 0._r8
         prod(:,34) = 0._r8
         prod(:,25) =.380_r8*rxt(:,25)*y(:,20) + extfrc(:,3)
         prod(:,2) =rxt(:,32)*y(:,6) +rxt(:,33)*y(:,7) +rxt(:,35)*y(:,9) &
                  +2.000_r8*rxt(:,36)*y(:,10) +2.000_r8*rxt(:,37)*y(:,11) +rxt(:,38) &
                 *y(:,12) +2.000_r8*rxt(:,51)*y(:,34) +rxt(:,54)*y(:,38) +rxt(:,55) &
                 *y(:,39)
         prod(:,5) =rxt(:,34)*y(:,8) +rxt(:,35)*y(:,9) +rxt(:,53)*y(:,37)
         prod(:,78) =.440_r8*rxt(:,25)*y(:,20)
         prod(:,30) =rxt(:,33)*y(:,7) +rxt(:,37)*y(:,11)
         prod(:,86) = (rxt(:,24) +.330_r8*rxt(:,25))*y(:,20)
         prod(:,50) =1.440_r8*rxt(:,25)*y(:,20)
         prod(:,18) = 0._r8
         prod(:,33) = 0._r8
         prod(:,80) = 0._r8
         prod(:,7) = 0._r8
         prod(:,71) = 0._r8
         prod(:,61) = 0._r8
         prod(:,14) = 0._r8
         prod(:,28) = 0._r8
         prod(:,32) = 0._r8
         prod(:,24) = 0._r8
         prod(:,40) = (.800_r8*rxt(:,67) +.800_r8*rxt(:,69) +rxt(:,73) +rxt(:,74)) &
                  + extfrc(:,11)
         prod(:,37) = 0._r8
         prod(:,87) = + extfrc(:,2)
         prod(:,72) = + extfrc(:,1)
         prod(:,74) = 0._r8
         prod(:,68) =.180_r8*rxt(:,25)*y(:,20) + extfrc(:,6)
         prod(:,75) = 0._r8
         prod(:,79) = 0._r8
         prod(:,3) = 0._r8
         prod(:,44) =rxt(:,32)*y(:,6) +rxt(:,33)*y(:,7) +2.000_r8*rxt(:,39)*y(:,13) &
                  +rxt(:,40)*y(:,15) +3.000_r8*rxt(:,43)*y(:,21) +2.000_r8*rxt(:,51) &
                 *y(:,34)
         prod(:,85) =4.000_r8*rxt(:,31)*y(:,5) +rxt(:,32)*y(:,6) +2.000_r8*rxt(:,34) &
                 *y(:,8) +2.000_r8*rxt(:,35)*y(:,9) +2.000_r8*rxt(:,36)*y(:,10) &
                  +rxt(:,37)*y(:,11) +2.000_r8*rxt(:,38)*y(:,12) +3.000_r8*rxt(:,41) &
                 *y(:,16) +rxt(:,42)*y(:,17) +rxt(:,53)*y(:,37) +rxt(:,54)*y(:,38) &
                  +rxt(:,55)*y(:,39)
         prod(:,58) = 0._r8
         prod(:,46) = 0._r8
         prod(:,45) = 0._r8
         prod(:,35) = 0._r8
         prod(:,64) = 0._r8
         prod(:,43) = 0._r8
         prod(:,54) = 0._r8
         prod(:,59) = 0._r8
         prod(:,69) = (rxt(:,68) +rxt(:,70) +rxt(:,71) +rxt(:,72) +rxt(:,73) + &
                 rxt(:,74)) + extfrc(:,10)
         prod(:,11) = 0._r8
         prod(:,42) = 0._r8
         prod(:,21) = 0._r8
         prod(:,65) = 0._r8
         prod(:,8) = 0._r8
         prod(:,73) = 0._r8
         prod(:,9) = 0._r8
         prod(:,82) = 0._r8
         prod(:,26) = 0._r8
         prod(:,29) = (1.200_r8*rxt(:,67) +1.200_r8*rxt(:,69) +rxt(:,70) +rxt(:,71)) &
                  + extfrc(:,12)
         prod(:,19) = (rxt(:,68) +rxt(:,72)) + extfrc(:,4)
         prod(:,62) = 0._r8
         prod(:,47) = 0._r8
         prod(:,60) = 0._r8
         prod(:,53) = 0._r8
         prod(:,55) = 0._r8
         prod(:,48) = 0._r8
         prod(:,51) = 0._r8
         prod(:,66) = 0._r8
         prod(:,67) = 0._r8
         prod(:,13) = 0._r8
         prod(:,23) = 0._r8
         prod(:,70) = 0._r8
         prod(:,22) = 0._r8
         prod(:,31) = (rxt(:,70) +rxt(:,71) +rxt(:,73) +rxt(:,74)) + extfrc(:,5)
         prod(:,57) =rxt(:,13)*y(:,49)
         prod(:,36) = 0._r8
         prod(:,10) = 0._r8
         prod(:,76) = 0._r8
         prod(:,77) = + extfrc(:,7)
         prod(:,27) = 0._r8
         prod(:,49) = 0._r8
         prod(:,20) = 0._r8
         prod(:,41) = 0._r8
         prod(:,63) =.330_r8*rxt(:,25)*y(:,20) + extfrc(:,8)
         prod(:,83) = 0._r8
         prod(:,84) = 0._r8
         prod(:,38) = + extfrc(:,9)
         prod(:,16) = 0._r8
         prod(:,15) = 0._r8
         prod(:,88) =.050_r8*rxt(:,25)*y(:,20)
      end if
      end subroutine indprd
      end module mo_indprd
