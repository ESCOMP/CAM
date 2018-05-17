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
         prod(:,:,2) =rxt(:,:,173)*y(:,:,7)*y(:,:,5)
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
         prod(:,:,21) = 0._r8
         prod(:,:,22) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,94) = 0._r8
         prod(:,:,95) =.180_r8*rxt(:,:,60)*y(:,:,12) + extfrc(:,:,20)
         prod(:,:,82) =rxt(:,:,5)*y(:,:,4)
         prod(:,:,97) = 0._r8
         prod(:,:,32) = 0._r8
         prod(:,:,61) = 0._r8
         prod(:,:,80) =1.440_r8*rxt(:,:,60)*y(:,:,12)
         prod(:,:,50) =.380_r8*rxt(:,:,60)*y(:,:,12) + extfrc(:,:,3)
         prod(:,:,102) =.440_r8*rxt(:,:,60)*y(:,:,12)
         prod(:,:,64) = (rxt(:,:,76) +.800_r8*rxt(:,:,79) +rxt(:,:,88) + &
                 .800_r8*rxt(:,:,91)) + extfrc(:,:,21)
         prod(:,:,104) = + extfrc(:,:,1)
         prod(:,:,105) = + extfrc(:,:,2)
         prod(:,:,88) =.330_r8*rxt(:,:,60)*y(:,:,12) + extfrc(:,:,23)
         prod(:,:,107) = 0._r8
         prod(:,:,51) = 0._r8
         prod(:,:,109) = 0._r8
         prod(:,:,39) = 0._r8
         prod(:,:,62) = 0._r8
         prod(:,:,63) =rxt(:,:,59)*y(:,:,12) +rxt(:,:,37)*y(:,:,35) +rxt(:,:,48) &
                 *y(:,:,36)
         prod(:,:,37) = 0._r8
         prod(:,:,73) =.180_r8*rxt(:,:,60)*y(:,:,12)
         prod(:,:,91) = (rxt(:,:,59) +.330_r8*rxt(:,:,60))*y(:,:,12)
         prod(:,:,86) = 0._r8
         prod(:,:,43) = 0._r8
         prod(:,:,100) =.050_r8*rxt(:,:,60)*y(:,:,12)
         prod(:,:,112) =rxt(:,:,37)*y(:,:,35) +2.000_r8*rxt(:,:,40)*y(:,:,37) &
                  +2.000_r8*rxt(:,:,41)*y(:,:,38) +2.000_r8*rxt(:,:,42)*y(:,:,39) &
                  +rxt(:,:,45)*y(:,:,40) +4.000_r8*rxt(:,:,38)*y(:,:,41) &
                  +3.000_r8*rxt(:,:,39)*y(:,:,42) +rxt(:,:,50)*y(:,:,44) +rxt(:,:,46) &
                 *y(:,:,45) +rxt(:,:,47)*y(:,:,46) +2.000_r8*rxt(:,:,43)*y(:,:,47) &
                  +rxt(:,:,44)*y(:,:,48)
         prod(:,:,24) = 0._r8
         prod(:,:,106) = 0._r8
         prod(:,:,33) = 0._r8
         prod(:,:,21) = 0._r8
         prod(:,:,93) = 0._r8
         prod(:,:,56) = 0._r8
         prod(:,:,59) = 0._r8
         prod(:,:,26) = 0._r8
         prod(:,:,81) =rxt(:,:,48)*y(:,:,36) +rxt(:,:,49)*y(:,:,43) +rxt(:,:,50) &
                 *y(:,:,44) +2.000_r8*rxt(:,:,53)*y(:,:,49) +2.000_r8*rxt(:,:,54) &
                 *y(:,:,50) +3.000_r8*rxt(:,:,51)*y(:,:,51) +2.000_r8*rxt(:,:,52) &
                 *y(:,:,52)
         prod(:,:,75) = 0._r8
         prod(:,:,55) = 0._r8
         prod(:,:,52) = 0._r8
         prod(:,:,41) = 0._r8
         prod(:,:,45) = (rxt(:,:,72) +rxt(:,:,84)) + extfrc(:,:,18)
         prod(:,:,110) = + extfrc(:,:,16)
         prod(:,:,38) = (rxt(:,:,76) +rxt(:,:,77) +rxt(:,:,88) +rxt(:,:,89)) &
                  + extfrc(:,:,17)
         prod(:,:,57) = + extfrc(:,:,15)
         prod(:,:,113) = 0._r8
         prod(:,:,40) = (rxt(:,:,77) +1.200_r8*rxt(:,:,79) +rxt(:,:,89) + &
                 1.200_r8*rxt(:,:,91)) + extfrc(:,:,19)
         prod(:,:,65) = 0._r8
         prod(:,:,49) = 0._r8
         prod(:,:,54) = 0._r8
         prod(:,:,44) = 0._r8
         prod(:,:,101) = 0._r8
         prod(:,:,108) = 0._r8
         prod(:,:,103) = 0._r8
         prod(:,:,35) = 0._r8
         prod(:,:,28) = 0._r8
         prod(:,:,29) = 0._r8
         prod(:,:,98) = 0._r8
         prod(:,:,90) = 0._r8
         prod(:,:,36) = 0._r8
         prod(:,:,47) = 0._r8
         prod(:,:,46) = 0._r8
         prod(:,:,99) = 0._r8
         prod(:,:,111) = 0._r8
         prod(:,:,74) = 0._r8
         prod(:,:,42) = 0._r8
         prod(:,:,89) = 0._r8
         prod(:,:,84) = 0._r8
         prod(:,:,87) = 0._r8
         prod(:,:,85) = 0._r8
         prod(:,:,92) = 0._r8
         prod(:,:,66) = 0._r8
         prod(:,:,78) = 0._r8
         prod(:,:,68) = 0._r8
         prod(:,:,71) = 0._r8
         prod(:,:,79) = 0._r8
         prod(:,:,77) = 0._r8
         prod(:,:,76) = 0._r8
         prod(:,:,72) = 0._r8
         prod(:,:,83) = 0._r8
         prod(:,:,60) = 0._r8
         prod(:,:,70) = 0._r8
         prod(:,:,69) = 0._r8
         prod(:,:,96) = (rxt(:,:,72) +rxt(:,:,76) +rxt(:,:,77) +rxt(:,:,84) + &
                 rxt(:,:,88) +rxt(:,:,89)) + extfrc(:,:,22)
         prod(:,:,22) =rxt(:,:,41)*y(:,:,38) +rxt(:,:,42)*y(:,:,39) +rxt(:,:,45) &
                 *y(:,:,40) +rxt(:,:,49)*y(:,:,43) +rxt(:,:,50)*y(:,:,44) +rxt(:,:,47) &
                 *y(:,:,46) +2.000_r8*rxt(:,:,43)*y(:,:,47) +2.000_r8*rxt(:,:,44) &
                 *y(:,:,48) +rxt(:,:,53)*y(:,:,49) +2.000_r8*rxt(:,:,54)*y(:,:,50)
         prod(:,:,25) =rxt(:,:,40)*y(:,:,37) +rxt(:,:,42)*y(:,:,39) +rxt(:,:,46) &
                 *y(:,:,45)
         prod(:,:,27) = 0._r8
         prod(:,:,53) =rxt(:,:,49)*y(:,:,43) +rxt(:,:,44)*y(:,:,48)
         prod(:,:,34) = 0._r8
         prod(:,:,48) = 0._r8
         prod(:,:,67) = 0._r8
         prod(:,:,58) = + extfrc(:,:,4)
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
