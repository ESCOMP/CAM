      module mo_indprd
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: indprd
      contains
      subroutine indprd( class, prod, nprod, y, extfrc, rxt, ncol )
      use chem_mods, only : gas_pcnst, extcnt, rxntot
      use ppgrid, only : pver
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(ncol,pver,gas_pcnst)
      real(r8), intent(in) :: rxt(ncol,pver,rxntot)
      real(r8), intent(in) :: extfrc(ncol,pver,extcnt)
      real(r8), intent(inout) :: prod(ncol,pver,nprod)
!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,:,1) = 0._r8
         prod(:,:,2) =rxt(:,:,158)*y(:,:,7)*y(:,:,5)
         prod(:,:,3) = 0._r8
         prod(:,:,4) = 0._r8
         prod(:,:,5) = 0._r8
         prod(:,:,6) = 0._r8
         prod(:,:,7) = 0._r8
         prod(:,:,8) = 0._r8
         prod(:,:,9) = 0._r8
         prod(:,:,10) = 0._r8
         prod(:,:,11) = 0._r8
         prod(:,:,12) = 0._r8
         prod(:,:,13) = 0._r8
         prod(:,:,14) = 0._r8
         prod(:,:,15) = 0._r8
         prod(:,:,16) = 0._r8
         prod(:,:,17) = 0._r8
         prod(:,:,18) = 0._r8
         prod(:,:,19) = 0._r8
         prod(:,:,20) = 0._r8
         prod(:,:,21) = (rxt(:,:,238)*y(:,:,86) +rxt(:,:,239)*y(:,:,86))*y(:,:,16)
         prod(:,:,22) = 0._r8
         prod(:,:,23) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,64) = 0._r8
         prod(:,:,75) = (rxt(:,:,58) +rxt(:,:,91))*y(:,:,56) +.180_r8*rxt(:,:,60) &
                 *y(:,:,12)
         prod(:,:,74) =rxt(:,:,5)*y(:,:,4)
         prod(:,:,70) = 0._r8
         prod(:,:,25) = 0._r8
         prod(:,:,24) = 0._r8
         prod(:,:,60) =1.440_r8*rxt(:,:,60)*y(:,:,12)
         prod(:,:,39) = (rxt(:,:,58) +rxt(:,:,91))*y(:,:,56) +.380_r8*rxt(:,:,60) &
                 *y(:,:,12) + extfrc(:,:,3)
         prod(:,:,51) = (rxt(:,:,75) +.800_r8*rxt(:,:,78) +rxt(:,:,87) + &
                 .800_r8*rxt(:,:,90)) + extfrc(:,:,20)
         prod(:,:,65) = + extfrc(:,:,1)
         prod(:,:,66) = + extfrc(:,:,2)
         prod(:,:,67) =.330_r8*rxt(:,:,60)*y(:,:,12) + extfrc(:,:,22)
         prod(:,:,68) = 0._r8
         prod(:,:,55) = 0._r8
         prod(:,:,37) = 0._r8
         prod(:,:,32) = 0._r8
         prod(:,:,72) =rxt(:,:,59)*y(:,:,12) +rxt(:,:,37)*y(:,:,34) +rxt(:,:,48) &
                 *y(:,:,35)
         prod(:,:,35) = 0._r8
         prod(:,:,58) =.180_r8*rxt(:,:,60)*y(:,:,12)
         prod(:,:,57) = (rxt(:,:,59) +.330_r8*rxt(:,:,60))*y(:,:,12)
         prod(:,:,71) = 0._r8
         prod(:,:,41) = 0._r8
         prod(:,:,69) =.050_r8*rxt(:,:,60)*y(:,:,12)
         prod(:,:,61) =rxt(:,:,37)*y(:,:,34) +2.000_r8*rxt(:,:,40)*y(:,:,36) &
                  +2.000_r8*rxt(:,:,41)*y(:,:,37) +2.000_r8*rxt(:,:,42)*y(:,:,38) &
                  +rxt(:,:,45)*y(:,:,39) +4.000_r8*rxt(:,:,38)*y(:,:,40) &
                  +3.000_r8*rxt(:,:,39)*y(:,:,41) +rxt(:,:,50)*y(:,:,43) +rxt(:,:,46) &
                 *y(:,:,44) +rxt(:,:,47)*y(:,:,45) +2.000_r8*rxt(:,:,43)*y(:,:,46) &
                  +rxt(:,:,44)*y(:,:,47)
         prod(:,:,26) = 0._r8
         prod(:,:,62) = 0._r8
         prod(:,:,33) = 0._r8
         prod(:,:,21) = 0._r8
         prod(:,:,73) = 0._r8
         prod(:,:,52) = 0._r8
         prod(:,:,54) = 0._r8
         prod(:,:,28) = 0._r8
         prod(:,:,63) =rxt(:,:,48)*y(:,:,35) +rxt(:,:,49)*y(:,:,42) +rxt(:,:,50) &
                 *y(:,:,43) +2.000_r8*rxt(:,:,53)*y(:,:,48) +2.000_r8*rxt(:,:,54) &
                 *y(:,:,49) +3.000_r8*rxt(:,:,51)*y(:,:,50) +2.000_r8*rxt(:,:,52) &
                 *y(:,:,51)
         prod(:,:,59) = 0._r8
         prod(:,:,50) = 0._r8
         prod(:,:,49) = 0._r8
         prod(:,:,40) = 0._r8
         prod(:,:,42) = (rxt(:,:,71) +rxt(:,:,83)) + extfrc(:,:,18)
         prod(:,:,46) = + extfrc(:,:,16)
         prod(:,:,36) = (rxt(:,:,75) +rxt(:,:,76) +rxt(:,:,87) +rxt(:,:,88)) &
                  + extfrc(:,:,17)
         prod(:,:,44) = + extfrc(:,:,15)
         prod(:,:,47) = 0._r8
         prod(:,:,38) = (rxt(:,:,76) +1.200_r8*rxt(:,:,78) +rxt(:,:,88) + &
                 1.200_r8*rxt(:,:,90)) + extfrc(:,:,19)
         prod(:,:,48) = (rxt(:,:,71) +rxt(:,:,75) +rxt(:,:,76) +rxt(:,:,83) + &
                 rxt(:,:,87) +rxt(:,:,88)) + extfrc(:,:,21)
         prod(:,:,22) =rxt(:,:,41)*y(:,:,37) +rxt(:,:,42)*y(:,:,38) +rxt(:,:,45) &
                 *y(:,:,39) +rxt(:,:,49)*y(:,:,42) +rxt(:,:,50)*y(:,:,43) +rxt(:,:,47) &
                 *y(:,:,45) +2.000_r8*rxt(:,:,43)*y(:,:,46) +2.000_r8*rxt(:,:,44) &
                 *y(:,:,47) +rxt(:,:,53)*y(:,:,48) +2.000_r8*rxt(:,:,54)*y(:,:,49)
         prod(:,:,27) =rxt(:,:,40)*y(:,:,36) +rxt(:,:,42)*y(:,:,38) +rxt(:,:,46) &
                 *y(:,:,44)
         prod(:,:,29) = 0._r8
         prod(:,:,45) =rxt(:,:,49)*y(:,:,42) +rxt(:,:,44)*y(:,:,47)
         prod(:,:,34) = 0._r8
         prod(:,:,43) = 0._r8
         prod(:,:,56) = 0._r8
         prod(:,:,53) = + extfrc(:,:,4)
         prod(:,:,30) = 0._r8
         prod(:,:,23) = 0._r8
         prod(:,:,31) = + extfrc(:,:,5)
         prod(:,:,1) = 0._r8
         prod(:,:,2) = + extfrc(:,:,6)
         prod(:,:,3) = + extfrc(:,:,8)
         prod(:,:,4) = 0._r8
         prod(:,:,5) = + extfrc(:,:,10)
         prod(:,:,6) = 0._r8
         prod(:,:,7) = 0._r8
         prod(:,:,8) = + extfrc(:,:,12)
         prod(:,:,9) = + extfrc(:,:,7)
         prod(:,:,10) = 0._r8
         prod(:,:,11) = 0._r8
         prod(:,:,12) = + extfrc(:,:,13)
         prod(:,:,13) = 0._r8
         prod(:,:,14) = 0._r8
         prod(:,:,15) = 0._r8
         prod(:,:,16) = 0._r8
         prod(:,:,17) = 0._r8
         prod(:,:,18) = + extfrc(:,:,9)
         prod(:,:,19) = + extfrc(:,:,11)
         prod(:,:,20) = + extfrc(:,:,14)
      end if
      end subroutine indprd
      end module mo_indprd
