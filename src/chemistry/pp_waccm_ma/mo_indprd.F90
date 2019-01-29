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
         prod(:,17) = (rxt(:,237)*y(:,69) +rxt(:,241)*y(:,69))*y(:,27)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) =rxt(:,127)*y(:,51)*y(:,47)
         prod(:,23) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,8) = 0._r8
         prod(:,39) = 0._r8
         prod(:,15) = 0._r8
         prod(:,35) =.180_r8*rxt(:,24)*y(:,20)
         prod(:,44) =rxt(:,39)*y(:,15) +rxt(:,41)*y(:,17) +rxt(:,23)*y(:,20)
         prod(:,11) = 0._r8
         prod(:,6) = 0._r8
         prod(:,1) = 0._r8
         prod(:,41) = 0._r8
         prod(:,30) = 0._r8
         prod(:,16) = (rxt(:,25) +rxt(:,61))*y(:,28) +.380_r8*rxt(:,24)*y(:,20) &
                  + extfrc(:,2)
         prod(:,2) =rxt(:,31)*y(:,6) +rxt(:,32)*y(:,7) +rxt(:,34)*y(:,9) &
                  +2.000_r8*rxt(:,35)*y(:,10) +2.000_r8*rxt(:,36)*y(:,11) +rxt(:,37) &
                 *y(:,12) +2.000_r8*rxt(:,50)*y(:,34) +rxt(:,53)*y(:,38) +rxt(:,54) &
                 *y(:,39)
         prod(:,7) =rxt(:,33)*y(:,8) +rxt(:,34)*y(:,9) +rxt(:,52)*y(:,37)
         prod(:,21) =rxt(:,32)*y(:,7) +rxt(:,36)*y(:,11)
         prod(:,36) = (rxt(:,23) +.330_r8*rxt(:,24))*y(:,20)
         prod(:,49) =1.440_r8*rxt(:,24)*y(:,20)
         prod(:,17) = 0._r8
         prod(:,23) = 0._r8
         prod(:,34) = 0._r8
         prod(:,9) = 0._r8
         prod(:,31) = 0._r8
         prod(:,47) = 0._r8
         prod(:,12) = 0._r8
         prod(:,20) = 0._r8
         prod(:,22) = 0._r8
         prod(:,29) = (.800_r8*rxt(:,64) +rxt(:,66) +.800_r8*rxt(:,68) +rxt(:,70)) &
                  + extfrc(:,7)
         prod(:,10) = 0._r8
         prod(:,38) = + extfrc(:,3)
         prod(:,45) = + extfrc(:,1)
         prod(:,37) = 0._r8
         prod(:,48) = (rxt(:,25) +rxt(:,61))*y(:,28) +.180_r8*rxt(:,24)*y(:,20)
         prod(:,32) = 0._r8
         prod(:,50) = 0._r8
         prod(:,3) = 0._r8
         prod(:,43) =rxt(:,31)*y(:,6) +rxt(:,32)*y(:,7) +2.000_r8*rxt(:,38)*y(:,13) &
                  +rxt(:,39)*y(:,15) +3.000_r8*rxt(:,42)*y(:,21) +2.000_r8*rxt(:,50) &
                 *y(:,34)
         prod(:,40) =4.000_r8*rxt(:,30)*y(:,5) +rxt(:,31)*y(:,6) +2.000_r8*rxt(:,33) &
                 *y(:,8) +2.000_r8*rxt(:,34)*y(:,9) +2.000_r8*rxt(:,35)*y(:,10) &
                  +rxt(:,36)*y(:,11) +2.000_r8*rxt(:,37)*y(:,12) +3.000_r8*rxt(:,40) &
                 *y(:,16) +rxt(:,41)*y(:,17) +rxt(:,52)*y(:,37) +rxt(:,53)*y(:,38) &
                  +rxt(:,54)*y(:,39)
         prod(:,27) = (rxt(:,63) +rxt(:,65) +rxt(:,66) +rxt(:,67) +rxt(:,69) + &
                 rxt(:,70)) + extfrc(:,9)
         prod(:,28) = (rxt(:,63) +1.200_r8*rxt(:,64) +1.200_r8*rxt(:,68) +rxt(:,69)) &
                  + extfrc(:,11)
         prod(:,18) = (rxt(:,65) +rxt(:,67)) + extfrc(:,6)
         prod(:,19) = 0._r8
         prod(:,24) = (rxt(:,63) +rxt(:,66) +rxt(:,69) +rxt(:,70)) + extfrc(:,10)
         prod(:,42) =rxt(:,12)*y(:,48)
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,26) = + extfrc(:,4)
         prod(:,46) =.330_r8*rxt(:,24)*y(:,20) + extfrc(:,8)
         prod(:,25) = + extfrc(:,5)
         prod(:,14) = 0._r8
         prod(:,13) = 0._r8
         prod(:,33) =.050_r8*rxt(:,24)*y(:,20)
      end if
      end subroutine indprd
      end module mo_indprd
