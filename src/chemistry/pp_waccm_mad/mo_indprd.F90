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
         prod(:,2) =rxt(:,169)*y(:,10)*y(:,8)
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
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,82) = 0._r8
         prod(:,88) =.180_r8*rxt(:,60)*y(:,15) + extfrc(:,12)
         prod(:,57) =rxt(:,5)*y(:,7)
         prod(:,85) = 0._r8
         prod(:,10) = 0._r8
         prod(:,36) = 0._r8
         prod(:,56) =1.440_r8*rxt(:,60)*y(:,15)
         prod(:,27) =.380_r8*rxt(:,60)*y(:,15) + extfrc(:,3)
         prod(:,81) =.440_r8*rxt(:,60)*y(:,15)
         prod(:,40) = (rxt(:,72) +.800_r8*rxt(:,75) +rxt(:,84) +.800_r8*rxt(:,87)) &
                  + extfrc(:,9)
         prod(:,73) = + extfrc(:,1)
         prod(:,87) = + extfrc(:,2)
         prod(:,63) =.330_r8*rxt(:,60)*y(:,15) + extfrc(:,11)
         prod(:,84) = 0._r8
         prod(:,31) = 0._r8
         prod(:,80) = 0._r8
         prod(:,14) = 0._r8
         prod(:,38) = 0._r8
         prod(:,39) =rxt(:,59)*y(:,15) +rxt(:,37)*y(:,43) +rxt(:,48)*y(:,44)
         prod(:,11) = 0._r8
         prod(:,49) =.180_r8*rxt(:,60)*y(:,15)
         prod(:,78) = (rxt(:,59) +.330_r8*rxt(:,60))*y(:,15)
         prod(:,61) = 0._r8
         prod(:,19) = 0._r8
         prod(:,70) =.050_r8*rxt(:,60)*y(:,15)
         prod(:,66) =rxt(:,37)*y(:,43) +2.000_r8*rxt(:,40)*y(:,45) +2.000_r8*rxt(:,41) &
                 *y(:,46) +2.000_r8*rxt(:,42)*y(:,47) +rxt(:,45)*y(:,48) &
                  +4.000_r8*rxt(:,38)*y(:,49) +3.000_r8*rxt(:,39)*y(:,50) +rxt(:,50) &
                 *y(:,52) +rxt(:,46)*y(:,53) +rxt(:,47)*y(:,54) +2.000_r8*rxt(:,43) &
                 *y(:,55) +rxt(:,44)*y(:,56)
         prod(:,4) = 0._r8
         prod(:,86) = 0._r8
         prod(:,2) = 0._r8
         prod(:,1) = 0._r8
         prod(:,76) = 0._r8
         prod(:,33) = 0._r8
         prod(:,35) = 0._r8
         prod(:,6) = 0._r8
         prod(:,55) =rxt(:,48)*y(:,44) +rxt(:,49)*y(:,51) +rxt(:,50)*y(:,52) &
                  +2.000_r8*rxt(:,53)*y(:,57) +2.000_r8*rxt(:,54)*y(:,58) &
                  +3.000_r8*rxt(:,51)*y(:,59) +2.000_r8*rxt(:,52)*y(:,60)
         prod(:,42) = 0._r8
         prod(:,30) = 0._r8
         prod(:,28) = 0._r8
         prod(:,18) = 0._r8
         prod(:,23) = (rxt(:,68) +rxt(:,80)) + extfrc(:,7)
         prod(:,67) = + extfrc(:,5)
         prod(:,20) = (rxt(:,72) +rxt(:,73) +rxt(:,84) +rxt(:,85)) + extfrc(:,6)
         prod(:,34) = + extfrc(:,4)
         prod(:,75) = 0._r8
         prod(:,17) = (rxt(:,73) +1.200_r8*rxt(:,75) +rxt(:,85) +1.200_r8*rxt(:,87)) &
                  + extfrc(:,8)
         prod(:,83) = (rxt(:,68) +rxt(:,72) +rxt(:,73) +rxt(:,80) +rxt(:,84) + &
                 rxt(:,85)) + extfrc(:,10)
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,41) = 0._r8
         prod(:,26) = 0._r8
         prod(:,32) = 0._r8
         prod(:,22) = 0._r8
         prod(:,74) = 0._r8
         prod(:,77) = 0._r8
         prod(:,79) = 0._r8
         prod(:,12) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,69) = 0._r8
         prod(:,71) = 0._r8
         prod(:,13) = 0._r8
         prod(:,25) = 0._r8
         prod(:,24) = 0._r8
         prod(:,68) = 0._r8
         prod(:,65) = 0._r8
         prod(:,50) = 0._r8
         prod(:,21) = 0._r8
         prod(:,64) = 0._r8
         prod(:,59) = 0._r8
         prod(:,62) = 0._r8
         prod(:,60) = 0._r8
         prod(:,72) = 0._r8
         prod(:,43) = 0._r8
         prod(:,53) = 0._r8
         prod(:,44) = 0._r8
         prod(:,47) = 0._r8
         prod(:,54) = 0._r8
         prod(:,52) = 0._r8
         prod(:,51) = 0._r8
         prod(:,48) = 0._r8
         prod(:,58) = 0._r8
         prod(:,37) = 0._r8
         prod(:,46) = 0._r8
         prod(:,45) = 0._r8
         prod(:,3) =rxt(:,41)*y(:,46) +rxt(:,42)*y(:,47) +rxt(:,45)*y(:,48) +rxt(:,49) &
                 *y(:,51) +rxt(:,50)*y(:,52) +rxt(:,47)*y(:,54) +2.000_r8*rxt(:,43) &
                 *y(:,55) +2.000_r8*rxt(:,44)*y(:,56) +rxt(:,53)*y(:,57) &
                  +2.000_r8*rxt(:,54)*y(:,58)
         prod(:,5) =rxt(:,40)*y(:,45) +rxt(:,42)*y(:,47) +rxt(:,46)*y(:,53)
         prod(:,7) = 0._r8
         prod(:,29) =rxt(:,49)*y(:,51) +rxt(:,44)*y(:,56)
      end if
      end subroutine indprd
      end module mo_indprd
