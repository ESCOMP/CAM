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
         prod(:,17) = (rxt(:,241)*y(:,95) +rxt(:,245)*y(:,95))*y(:,29)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) =rxt(:,131)*y(:,60)*y(:,53)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,1) = + extfrc(:,14)
         prod(:,2) = + extfrc(:,3)
         prod(:,29) = 0._r8
         prod(:,63) = 0._r8
         prod(:,40) = 0._r8
         prod(:,67) =.180_r8*rxt(:,24)*y(:,22)
         prod(:,60) =rxt(:,39)*y(:,17) +rxt(:,41)*y(:,19) +rxt(:,23)*y(:,22)
         prod(:,35) = 0._r8
         prod(:,26) = 0._r8
         prod(:,21) = 0._r8
         prod(:,66) = 0._r8
         prod(:,54) = 0._r8
         prod(:,39) = (rxt(:,25) +rxt(:,60))*y(:,30) +.380_r8*rxt(:,24)*y(:,22) &
                  + extfrc(:,10)
         prod(:,23) =rxt(:,31)*y(:,8) +rxt(:,32)*y(:,9) +rxt(:,34)*y(:,11) &
                  +2.000_r8*rxt(:,35)*y(:,12) +2.000_r8*rxt(:,36)*y(:,13) +rxt(:,37) &
                 *y(:,14) +2.000_r8*rxt(:,50)*y(:,40) +rxt(:,53)*y(:,45) +rxt(:,54) &
                 *y(:,46)
         prod(:,27) =rxt(:,33)*y(:,10) +rxt(:,34)*y(:,11) +rxt(:,52)*y(:,44)
         prod(:,32) = + extfrc(:,2)
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,46) =rxt(:,32)*y(:,9) +rxt(:,36)*y(:,13)
         prod(:,57) = (rxt(:,23) +.330_r8*rxt(:,24))*y(:,22)
         prod(:,70) =1.440_r8*rxt(:,24)*y(:,22)
         prod(:,41) = 0._r8
         prod(:,22) = 0._r8
         prod(:,51) = 0._r8
         prod(:,58) = 0._r8
         prod(:,28) = 0._r8
         prod(:,55) = 0._r8
         prod(:,36) = 0._r8
         prod(:,45) = 0._r8
         prod(:,52) = 0._r8
         prod(:,50) = (.800_r8*rxt(:,63) +rxt(:,65) +.800_r8*rxt(:,67) +rxt(:,69)) &
                  + extfrc(:,15)
         prod(:,31) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,68) = + extfrc(:,11)
         prod(:,65) = + extfrc(:,12)
         prod(:,61) = 0._r8
         prod(:,9) = + extfrc(:,4)
         prod(:,10) = + extfrc(:,5)
         prod(:,11) = 0._r8
         prod(:,12) = + extfrc(:,6)
         prod(:,74) = (rxt(:,25) +rxt(:,60))*y(:,30) +.180_r8*rxt(:,24)*y(:,22)
         prod(:,71) = 0._r8
         prod(:,72) = 0._r8
         prod(:,33) = 0._r8
         prod(:,34) = 0._r8
         prod(:,13) = + extfrc(:,7)
         prod(:,14) = + extfrc(:,8)
         prod(:,43) = 0._r8
         prod(:,56) = 0._r8
         prod(:,53) = + extfrc(:,13)
         prod(:,30) = 0._r8
         prod(:,15) = + extfrc(:,9)
         prod(:,16) = + extfrc(:,1)
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,59) =rxt(:,31)*y(:,8) +rxt(:,32)*y(:,9) +2.000_r8*rxt(:,38)*y(:,15) &
                  +rxt(:,39)*y(:,17) +3.000_r8*rxt(:,42)*y(:,23) +2.000_r8*rxt(:,50) &
                 *y(:,40)
         prod(:,62) =4.000_r8*rxt(:,30)*y(:,7) +rxt(:,31)*y(:,8) +2.000_r8*rxt(:,33) &
                 *y(:,10) +2.000_r8*rxt(:,34)*y(:,11) +2.000_r8*rxt(:,35)*y(:,12) &
                  +rxt(:,36)*y(:,13) +2.000_r8*rxt(:,37)*y(:,14) +3.000_r8*rxt(:,40) &
                 *y(:,18) +rxt(:,41)*y(:,19) +rxt(:,52)*y(:,44) +rxt(:,53)*y(:,45) &
                  +rxt(:,54)*y(:,46)
         prod(:,49) = (rxt(:,62) +rxt(:,64) +rxt(:,65) +rxt(:,66) +rxt(:,68) + &
                 rxt(:,69)) + extfrc(:,19)
         prod(:,64) = 0._r8
         prod(:,38) = (rxt(:,62) +1.200_r8*rxt(:,63) +1.200_r8*rxt(:,67) +rxt(:,68)) &
                  + extfrc(:,16)
         prod(:,42) = (rxt(:,64) +rxt(:,66)) + extfrc(:,17)
         prod(:,48) = 0._r8
         prod(:,37) = (rxt(:,62) +rxt(:,65) +rxt(:,68) +rxt(:,69)) + extfrc(:,20)
         prod(:,69) =rxt(:,12)*y(:,54)
         prod(:,24) = 0._r8
         prod(:,25) = 0._r8
         prod(:,47) = + extfrc(:,21)
         prod(:,73) =.330_r8*rxt(:,24)*y(:,22) + extfrc(:,22)
         prod(:,44) = + extfrc(:,18)
         prod(:,75) =.050_r8*rxt(:,24)*y(:,22)
      end if
      end subroutine indprd
      end module mo_indprd
