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
         prod(:,17) = (rxt(:,242)*y(:,77) +rxt(:,246)*y(:,77))*y(:,29)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) =rxt(:,132)*y(:,54)*y(:,50)
         prod(:,23) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,41) =rxt(:,31)*y(:,7) +rxt(:,32)*y(:,8) +2.000_r8*rxt(:,38)*y(:,14) &
                  +rxt(:,39)*y(:,16) +3.000_r8*rxt(:,42)*y(:,22) +2.000_r8*rxt(:,50) &
                 *y(:,37)
         prod(:,9) = 0._r8
         prod(:,57) = 0._r8
         prod(:,20) = 0._r8
         prod(:,51) =.180_r8*rxt(:,24)*y(:,21)
         prod(:,43) =rxt(:,39)*y(:,16) +rxt(:,41)*y(:,18) +rxt(:,23)*y(:,21)
         prod(:,16) = 0._r8
         prod(:,48) =4.000_r8*rxt(:,30)*y(:,6) +rxt(:,31)*y(:,7) +2.000_r8*rxt(:,33) &
                 *y(:,9) +2.000_r8*rxt(:,34)*y(:,10) +2.000_r8*rxt(:,35)*y(:,11) &
                  +rxt(:,36)*y(:,12) +2.000_r8*rxt(:,37)*y(:,13) +3.000_r8*rxt(:,40) &
                 *y(:,17) +rxt(:,41)*y(:,18) +rxt(:,52)*y(:,41) +rxt(:,53)*y(:,42) &
                  +rxt(:,54)*y(:,43)
         prod(:,7) = 0._r8
         prod(:,2) = 0._r8
         prod(:,47) = 0._r8
         prod(:,37) = 0._r8
         prod(:,21) = (rxt(:,25) +rxt(:,61))*y(:,30) +.380_r8*rxt(:,24)*y(:,21) &
                  + extfrc(:,2)
         prod(:,3) =rxt(:,31)*y(:,7) +rxt(:,32)*y(:,8) +rxt(:,34)*y(:,10) &
                  +2.000_r8*rxt(:,35)*y(:,11) +2.000_r8*rxt(:,36)*y(:,12) +rxt(:,37) &
                 *y(:,13) +2.000_r8*rxt(:,50)*y(:,37) +rxt(:,53)*y(:,42) +rxt(:,54) &
                 *y(:,43)
         prod(:,8) =rxt(:,33)*y(:,9) +rxt(:,34)*y(:,10) +rxt(:,52)*y(:,41)
         prod(:,12) = + extfrc(:,1)
         prod(:,26) =rxt(:,32)*y(:,8) +rxt(:,36)*y(:,12)
         prod(:,40) = (rxt(:,23) +.330_r8*rxt(:,24))*y(:,21)
         prod(:,56) =1.440_r8*rxt(:,24)*y(:,21)
         prod(:,22) = 0._r8
         prod(:,4) = 0._r8
         prod(:,29) = 0._r8
         prod(:,42) = 0._r8
         prod(:,10) = 0._r8
         prod(:,38) = 0._r8
         prod(:,17) = 0._r8
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,35) = (rxt(:,64) +.800_r8*rxt(:,66) +.800_r8*rxt(:,68) +rxt(:,70)) &
                  + extfrc(:,6)
         prod(:,13) = 0._r8
         prod(:,55) = + extfrc(:,3)
         prod(:,53) = + extfrc(:,4)
         prod(:,44) = 0._r8
         prod(:,46) = (rxt(:,25) +rxt(:,61))*y(:,30) +.180_r8*rxt(:,24)*y(:,21)
         prod(:,49) = 0._r8
         prod(:,52) = 0._r8
         prod(:,14) = 0._r8
         prod(:,15) = 0._r8
         prod(:,24) = 0._r8
         prod(:,39) = 0._r8
         prod(:,36) = + extfrc(:,5)
         prod(:,11) = 0._r8
         prod(:,1) = 0._r8
         prod(:,33) = (rxt(:,63) +rxt(:,64) +rxt(:,65) +rxt(:,67) +rxt(:,69) + &
                 rxt(:,70)) + extfrc(:,10)
         prod(:,45) = 0._r8
         prod(:,34) = (rxt(:,65) +1.200_r8*rxt(:,66) +1.200_r8*rxt(:,68) +rxt(:,69)) &
                  + extfrc(:,7)
         prod(:,23) = (rxt(:,63) +rxt(:,67)) + extfrc(:,8)
         prod(:,25) = 0._r8
         prod(:,30) = (rxt(:,64) +rxt(:,65) +rxt(:,69) +rxt(:,70)) + extfrc(:,11)
         prod(:,50) =rxt(:,12)*y(:,51)
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,32) = + extfrc(:,12)
         prod(:,54) =.330_r8*rxt(:,24)*y(:,21) + extfrc(:,13)
         prod(:,31) = + extfrc(:,9)
         prod(:,19) = 0._r8
         prod(:,18) = 0._r8
         prod(:,58) =.050_r8*rxt(:,24)*y(:,21)
      end if
      end subroutine indprd
      end module mo_indprd
