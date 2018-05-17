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
         prod(:,:,2) =rxt(:,:,154)*y(:,:,10)*y(:,:,8)
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
         prod(:,:,21) = (rxt(:,:,234)*y(:,:,22) +rxt(:,:,235)*y(:,:,22))*y(:,:,19)
         prod(:,:,22) = 0._r8
         prod(:,:,23) = 0._r8
         prod(:,:,24) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,49) = 0._r8
         prod(:,:,44) = (rxt(:,:,58) +rxt(:,:,87))*y(:,:,64) +.180_r8*rxt(:,:,60) &
                 *y(:,:,15)
         prod(:,:,47) =rxt(:,:,5)*y(:,:,7)
         prod(:,:,32) = 0._r8
         prod(:,:,3) = 0._r8
         prod(:,:,2) = 0._r8
         prod(:,:,35) =1.440_r8*rxt(:,:,60)*y(:,:,15)
         prod(:,:,16) = (rxt(:,:,58) +rxt(:,:,87))*y(:,:,64) +.380_r8*rxt(:,:,60) &
                 *y(:,:,15) + extfrc(:,:,3)
         prod(:,:,28) = (rxt(:,:,71) +.800_r8*rxt(:,:,74) +rxt(:,:,83) + &
                 .800_r8*rxt(:,:,86)) + extfrc(:,:,9)
         prod(:,:,48) = + extfrc(:,:,1)
         prod(:,:,40) = + extfrc(:,:,2)
         prod(:,:,43) =.330_r8*rxt(:,:,60)*y(:,:,15) + extfrc(:,:,11)
         prod(:,:,45) = 0._r8
         prod(:,:,31) = 0._r8
         prod(:,:,12) = 0._r8
         prod(:,:,10) = 0._r8
         prod(:,:,50) =rxt(:,:,59)*y(:,:,15) +rxt(:,:,37)*y(:,:,42) +rxt(:,:,48) &
                 *y(:,:,43)
         prod(:,:,11) = 0._r8
         prod(:,:,34) =.180_r8*rxt(:,:,60)*y(:,:,15)
         prod(:,:,38) = (rxt(:,:,59) +.330_r8*rxt(:,:,60))*y(:,:,15)
         prod(:,:,41) = 0._r8
         prod(:,:,18) = 0._r8
         prod(:,:,33) =.050_r8*rxt(:,:,60)*y(:,:,15)
         prod(:,:,42) =rxt(:,:,37)*y(:,:,42) +2.000_r8*rxt(:,:,40)*y(:,:,44) &
                  +2.000_r8*rxt(:,:,41)*y(:,:,45) +2.000_r8*rxt(:,:,42)*y(:,:,46) &
                  +rxt(:,:,45)*y(:,:,47) +4.000_r8*rxt(:,:,38)*y(:,:,48) &
                  +3.000_r8*rxt(:,:,39)*y(:,:,49) +rxt(:,:,50)*y(:,:,51) +rxt(:,:,46) &
                 *y(:,:,52) +rxt(:,:,47)*y(:,:,53) +2.000_r8*rxt(:,:,43)*y(:,:,54) &
                  +rxt(:,:,44)*y(:,:,55)
         prod(:,:,6) = 0._r8
         prod(:,:,37) = 0._r8
         prod(:,:,4) = 0._r8
         prod(:,:,1) = 0._r8
         prod(:,:,46) = 0._r8
         prod(:,:,29) = 0._r8
         prod(:,:,30) = 0._r8
         prod(:,:,8) = 0._r8
         prod(:,:,39) =rxt(:,:,48)*y(:,:,43) +rxt(:,:,49)*y(:,:,50) +rxt(:,:,50) &
                 *y(:,:,51) +2.000_r8*rxt(:,:,53)*y(:,:,56) +2.000_r8*rxt(:,:,54) &
                 *y(:,:,57) +3.000_r8*rxt(:,:,51)*y(:,:,58) +2.000_r8*rxt(:,:,52) &
                 *y(:,:,59)
         prod(:,:,36) = 0._r8
         prod(:,:,27) = 0._r8
         prod(:,:,22) = 0._r8
         prod(:,:,17) = 0._r8
         prod(:,:,20) = (rxt(:,:,67) +rxt(:,:,79)) + extfrc(:,:,7)
         prod(:,:,23) = + extfrc(:,:,5)
         prod(:,:,19) = (rxt(:,:,71) +rxt(:,:,72) +rxt(:,:,83) +rxt(:,:,84)) &
                  + extfrc(:,:,6)
         prod(:,:,21) = + extfrc(:,:,4)
         prod(:,:,24) = 0._r8
         prod(:,:,15) = (rxt(:,:,72) +1.200_r8*rxt(:,:,74) +rxt(:,:,84) + &
                 1.200_r8*rxt(:,:,86)) + extfrc(:,:,8)
         prod(:,:,25) = (rxt(:,:,67) +rxt(:,:,71) +rxt(:,:,72) +rxt(:,:,79) + &
                 rxt(:,:,83) +rxt(:,:,84)) + extfrc(:,:,10)
         prod(:,:,13) = 0._r8
         prod(:,:,14) = 0._r8
         prod(:,:,5) =rxt(:,:,41)*y(:,:,45) +rxt(:,:,42)*y(:,:,46) +rxt(:,:,45) &
                 *y(:,:,47) +rxt(:,:,49)*y(:,:,50) +rxt(:,:,50)*y(:,:,51) +rxt(:,:,47) &
                 *y(:,:,53) +2.000_r8*rxt(:,:,43)*y(:,:,54) +2.000_r8*rxt(:,:,44) &
                 *y(:,:,55) +rxt(:,:,53)*y(:,:,56) +2.000_r8*rxt(:,:,54)*y(:,:,57)
         prod(:,:,7) =rxt(:,:,40)*y(:,:,44) +rxt(:,:,42)*y(:,:,46) +rxt(:,:,46) &
                 *y(:,:,52)
         prod(:,:,9) = 0._r8
         prod(:,:,26) =rxt(:,:,49)*y(:,:,50) +rxt(:,:,44)*y(:,:,55)
      end if
      end subroutine indprd
      end module mo_indprd
