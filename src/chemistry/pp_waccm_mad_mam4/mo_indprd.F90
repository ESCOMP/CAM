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
         prod(:,:,2) = 0._r8
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
         prod(:,:,21) =rxt(:,:,427)*y(:,:,63)*y(:,:,56)
         prod(:,:,22) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,1) = + extfrc(:,:,5)
         prod(:,:,2) = + extfrc(:,:,6)
         prod(:,:,75) =rxt(:,:,32)*y(:,:,9) +rxt(:,:,33)*y(:,:,10) &
                  +2.000_r8*rxt(:,:,39)*y(:,:,16) +rxt(:,:,40)*y(:,:,18) &
                  +3.000_r8*rxt(:,:,43)*y(:,:,24) +2.000_r8*rxt(:,:,51)*y(:,:,42)
         prod(:,:,27) = 0._r8
         prod(:,:,82) = 0._r8
         prod(:,:,41) = 0._r8
         prod(:,:,83) =.180_r8*rxt(:,:,25)*y(:,:,23)
         prod(:,:,65) =rxt(:,:,40)*y(:,:,18) +rxt(:,:,42)*y(:,:,20) +rxt(:,:,24) &
                 *y(:,:,23)
         prod(:,:,35) = 0._r8
         prod(:,:,93) =4.000_r8*rxt(:,:,31)*y(:,:,8) +rxt(:,:,32)*y(:,:,9) &
                  +2.000_r8*rxt(:,:,34)*y(:,:,11) +2.000_r8*rxt(:,:,35)*y(:,:,12) &
                  +2.000_r8*rxt(:,:,36)*y(:,:,13) +rxt(:,:,37)*y(:,:,14) &
                  +2.000_r8*rxt(:,:,38)*y(:,:,15) +3.000_r8*rxt(:,:,41)*y(:,:,19) &
                  +rxt(:,:,42)*y(:,:,20) +rxt(:,:,53)*y(:,:,46) +rxt(:,:,54)*y(:,:,47) &
                  +rxt(:,:,55)*y(:,:,48)
         prod(:,:,24) = 0._r8
         prod(:,:,21) = 0._r8
         prod(:,:,113) = 0._r8
         prod(:,:,60) = 0._r8
         prod(:,:,49) =.380_r8*rxt(:,:,25)*y(:,:,23) + extfrc(:,:,13)
         prod(:,:,23) =rxt(:,:,32)*y(:,:,9) +rxt(:,:,33)*y(:,:,10) +rxt(:,:,35) &
                 *y(:,:,12) +2.000_r8*rxt(:,:,36)*y(:,:,13) +2.000_r8*rxt(:,:,37) &
                 *y(:,:,14) +rxt(:,:,38)*y(:,:,15) +2.000_r8*rxt(:,:,51)*y(:,:,42) &
                  +rxt(:,:,54)*y(:,:,47) +rxt(:,:,55)*y(:,:,48)
         prod(:,:,25) =rxt(:,:,34)*y(:,:,11) +rxt(:,:,35)*y(:,:,12) +rxt(:,:,53) &
                 *y(:,:,46)
         prod(:,:,96) =.440_r8*rxt(:,:,25)*y(:,:,23)
         prod(:,:,31) = + extfrc(:,:,2)
         prod(:,:,3) = 0._r8
         prod(:,:,4) = 0._r8
         prod(:,:,5) = 0._r8
         prod(:,:,53) =rxt(:,:,33)*y(:,:,10) +rxt(:,:,37)*y(:,:,14)
         prod(:,:,99) = (rxt(:,:,24) +.330_r8*rxt(:,:,25))*y(:,:,23)
         prod(:,:,77) =1.440_r8*rxt(:,:,25)*y(:,:,23)
         prod(:,:,42) = 0._r8
         prod(:,:,22) = 0._r8
         prod(:,:,57) = 0._r8
         prod(:,:,108) = 0._r8
         prod(:,:,26) = 0._r8
         prod(:,:,112) = 0._r8
         prod(:,:,38) = 0._r8
         prod(:,:,55) = 0._r8
         prod(:,:,58) = 0._r8
         prod(:,:,52) = 0._r8
         prod(:,:,66) = (rxt(:,:,69) +rxt(:,:,70) +.800_r8*rxt(:,:,72) + &
                 .800_r8*rxt(:,:,73)) + extfrc(:,:,15)
         prod(:,:,63) = 0._r8
         prod(:,:,6) = 0._r8
         prod(:,:,7) = 0._r8
         prod(:,:,8) = 0._r8
         prod(:,:,102) = + extfrc(:,:,14)
         prod(:,:,97) = + extfrc(:,:,3)
         prod(:,:,101) = 0._r8
         prod(:,:,9) = + extfrc(:,:,7)
         prod(:,:,10) = + extfrc(:,:,8)
         prod(:,:,11) = 0._r8
         prod(:,:,12) = + extfrc(:,:,9)
         prod(:,:,100) =.180_r8*rxt(:,:,25)*y(:,:,23) + extfrc(:,:,22)
         prod(:,:,84) =rxt(:,:,13)*y(:,:,57)
         prod(:,:,107) = 0._r8
         prod(:,:,105) = 0._r8
         prod(:,:,33) = 0._r8
         prod(:,:,34) = 0._r8
         prod(:,:,13) = + extfrc(:,:,10)
         prod(:,:,14) = + extfrc(:,:,11)
         prod(:,:,48) = 0._r8
         prod(:,:,68) = 0._r8
         prod(:,:,59) = + extfrc(:,:,4)
         prod(:,:,28) = 0._r8
         prod(:,:,15) = + extfrc(:,:,12)
         prod(:,:,16) = + extfrc(:,:,1)
         prod(:,:,17) = 0._r8
         prod(:,:,18) = 0._r8
         prod(:,:,19) = 0._r8
         prod(:,:,20) = 0._r8
         prod(:,:,85) = 0._r8
         prod(:,:,72) = 0._r8
         prod(:,:,71) = 0._r8
         prod(:,:,61) = 0._r8
         prod(:,:,91) = 0._r8
         prod(:,:,70) = 0._r8
         prod(:,:,80) = 0._r8
         prod(:,:,86) = 0._r8
         prod(:,:,92) = (rxt(:,:,67) +rxt(:,:,68) +rxt(:,:,69) +rxt(:,:,70) + &
                 rxt(:,:,71) +rxt(:,:,74)) + extfrc(:,:,20)
         prod(:,:,36) = 0._r8
         prod(:,:,69) = 0._r8
         prod(:,:,88) = 0._r8
         prod(:,:,45) = 0._r8
         prod(:,:,109) = 0._r8
         prod(:,:,29) = 0._r8
         prod(:,:,106) = 0._r8
         prod(:,:,30) = 0._r8
         prod(:,:,114) = 0._r8
         prod(:,:,50) = 0._r8
         prod(:,:,54) = (rxt(:,:,68) +rxt(:,:,71) +1.200_r8*rxt(:,:,72) + &
                 1.200_r8*rxt(:,:,73)) + extfrc(:,:,16)
         prod(:,:,43) = (rxt(:,:,67) +rxt(:,:,74)) + extfrc(:,:,17)
         prod(:,:,89) = 0._r8
         prod(:,:,73) = 0._r8
         prod(:,:,87) = 0._r8
         prod(:,:,79) = 0._r8
         prod(:,:,81) = 0._r8
         prod(:,:,74) = 0._r8
         prod(:,:,78) = 0._r8
         prod(:,:,94) = 0._r8
         prod(:,:,95) = 0._r8
         prod(:,:,37) = 0._r8
         prod(:,:,47) = 0._r8
         prod(:,:,98) = 0._r8
         prod(:,:,46) = 0._r8
         prod(:,:,56) = (rxt(:,:,68) +rxt(:,:,69) +rxt(:,:,70) +rxt(:,:,71)) &
                  + extfrc(:,:,21)
         prod(:,:,62) = 0._r8
         prod(:,:,32) = 0._r8
         prod(:,:,103) = 0._r8
         prod(:,:,104) = + extfrc(:,:,23)
         prod(:,:,51) = 0._r8
         prod(:,:,76) = 0._r8
         prod(:,:,44) = 0._r8
         prod(:,:,67) = 0._r8
         prod(:,:,90) =.330_r8*rxt(:,:,25)*y(:,:,23) + extfrc(:,:,18)
         prod(:,:,110) = 0._r8
         prod(:,:,111) = 0._r8
         prod(:,:,64) = + extfrc(:,:,19)
         prod(:,:,40) = 0._r8
         prod(:,:,39) = 0._r8
         prod(:,:,115) =.050_r8*rxt(:,:,25)*y(:,:,23)
      end if
      end subroutine indprd
      end module mo_indprd
