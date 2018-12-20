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
         prod(:,1) = + extfrc(:,17)
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
         prod(:,15) =.100_r8*rxt(:,308)*y(:,134)*y(:,29)
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = (rxt(:,265)*y(:,62) +rxt(:,267)*y(:,87) +rxt(:,275)*y(:,62) + &
                 rxt(:,295)*y(:,50) +.500_r8*rxt(:,296)*y(:,51) + &
                 .800_r8*rxt(:,301)*y(:,74) +rxt(:,302)*y(:,75) + &
                 .500_r8*rxt(:,351)*y(:,109) +1.800_r8*rxt(:,461)*y(:,178))*y(:,212) &
                  + (2.000_r8*rxt(:,291)*y(:,195) +.900_r8*rxt(:,292)*y(:,196) + &
                 rxt(:,294)*y(:,124) +2.000_r8*rxt(:,341)*y(:,207) + &
                 rxt(:,365)*y(:,203) +rxt(:,390)*y(:,219))*y(:,195) &
                  + (.200_r8*rxt(:,308)*y(:,29) +.100_r8*rxt(:,352)*y(:,111) + &
                 .270_r8*rxt(:,440)*y(:,6) +.270_r8*rxt(:,443)*y(:,110))*y(:,134) &
                  + (rxt(:,342)*y(:,196) +.450_r8*rxt(:,343)*y(:,201) + &
                 2.000_r8*rxt(:,344)*y(:,207))*y(:,207) &
                  + (.500_r8*rxt(:,450)*y(:,196) +.900_r8*rxt(:,452)*y(:,124)) &
                 *y(:,216) +rxt(:,37)*y(:,51) +.400_r8*rxt(:,60)*y(:,139) +rxt(:,65) &
                 *y(:,174) +.800_r8*rxt(:,69)*y(:,178)
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) =rxt(:,151)*y(:,125)*y(:,112)
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,29) = 0._r8
         prod(:,30) =rxt(:,478)*y(:,212)*y(:,120) +rxt(:,488)*y(:,121)
         prod(:,31) = (rxt(:,412)*y(:,197) +rxt(:,415)*y(:,206) +rxt(:,418)*y(:,208) + &
                 rxt(:,422)*y(:,141))*y(:,125) +.500_r8*rxt(:,351)*y(:,212)*y(:,109) &
                  +.200_r8*rxt(:,447)*y(:,210)*y(:,124) +.500_r8*rxt(:,459)*y(:,177) &
                 *y(:,126)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,121) = 0._r8
         prod(:,120) = 0._r8
         prod(:,1) = + extfrc(:,10)
         prod(:,2) = + extfrc(:,11)
         prod(:,144) = 0._r8
         prod(:,46) = 0._r8
         prod(:,85) = 0._r8
         prod(:,47) = 0._r8
         prod(:,86) = 0._r8
         prod(:,96) = 0._r8
         prod(:,68) = 0._r8
         prod(:,117) = 0._r8
         prod(:,75) = 0._r8
         prod(:,60) = 0._r8
         prod(:,81) = 0._r8
         prod(:,174) =rxt(:,79)*y(:,34) +rxt(:,80)*y(:,35) +2.000_r8*rxt(:,86)*y(:,41) &
                  +rxt(:,87)*y(:,43) +3.000_r8*rxt(:,90)*y(:,55) +2.000_r8*rxt(:,98) &
                 *y(:,78)
         prod(:,61) = 0._r8
         prod(:,187) = 0._r8
         prod(:,110) = 0._r8
         prod(:,62) = 0._r8
         prod(:,78) = 0._r8
         prod(:,70) = 0._r8
         prod(:,111) = 0._r8
         prod(:,64) = 0._r8
         prod(:,79) = 0._r8
         prod(:,71) = 0._r8
         prod(:,149) = 0._r8
         prod(:,90) = 0._r8
         prod(:,39) = 0._r8
         prod(:,65) = 0._r8
         prod(:,189) =.180_r8*rxt(:,40)*y(:,54)
         prod(:,161) = 0._r8
         prod(:,38) = 0._r8
         prod(:,147) = 0._r8
         prod(:,166) = 0._r8
         prod(:,108) = 0._r8
         prod(:,102) = 0._r8
         prod(:,134) = 0._r8
         prod(:,91) = 0._r8
         prod(:,184) =4.000_r8*rxt(:,78)*y(:,33) +rxt(:,79)*y(:,34) &
                  +2.000_r8*rxt(:,81)*y(:,36) +2.000_r8*rxt(:,82)*y(:,37) &
                  +2.000_r8*rxt(:,83)*y(:,38) +rxt(:,84)*y(:,39) +2.000_r8*rxt(:,85) &
                 *y(:,40) +3.000_r8*rxt(:,88)*y(:,44) +rxt(:,89)*y(:,46) +rxt(:,100) &
                 *y(:,82) +rxt(:,101)*y(:,83) +rxt(:,102)*y(:,84)
         prod(:,45) = 0._r8
         prod(:,36) = 0._r8
         prod(:,183) = 0._r8
         prod(:,148) = 0._r8
         prod(:,155) =.380_r8*rxt(:,40)*y(:,54) +rxt(:,41)*y(:,63) + extfrc(:,9)
         prod(:,40) =rxt(:,79)*y(:,34) +rxt(:,80)*y(:,35) +rxt(:,82)*y(:,37) &
                  +2.000_r8*rxt(:,83)*y(:,38) +2.000_r8*rxt(:,84)*y(:,39) +rxt(:,85) &
                 *y(:,40) +2.000_r8*rxt(:,98)*y(:,78) +rxt(:,101)*y(:,83) +rxt(:,102) &
                 *y(:,84)
         prod(:,50) =rxt(:,81)*y(:,36) +rxt(:,82)*y(:,37) +rxt(:,100)*y(:,82)
         prod(:,52) = 0._r8
         prod(:,69) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,41) = 0._r8
         prod(:,132) =rxt(:,80)*y(:,35) +rxt(:,84)*y(:,39)
         prod(:,151) = 0._r8
         prod(:,142) = 0._r8
         prod(:,176) = (rxt(:,39) +.330_r8*rxt(:,40))*y(:,54)
         prod(:,162) =1.440_r8*rxt(:,40)*y(:,54)
         prod(:,114) = 0._r8
         prod(:,42) = 0._r8
         prod(:,138) = 0._r8
         prod(:,177) = 0._r8
         prod(:,49) = 0._r8
         prod(:,133) = 0._r8
         prod(:,57) = 0._r8
         prod(:,175) = 0._r8
         prod(:,83) = 0._r8
         prod(:,131) = 0._r8
         prod(:,136) = 0._r8
         prod(:,154) = 0._r8
         prod(:,58) = 0._r8
         prod(:,156) = 0._r8
         prod(:,72) = 0._r8
         prod(:,43) = 0._r8
         prod(:,139) = 0._r8
         prod(:,113) = 0._r8
         prod(:,107) = 0._r8
         prod(:,164) = 0._r8
         prod(:,89) = 0._r8
         prod(:,123) = 0._r8
         prod(:,34) = 0._r8
         prod(:,165) = 0._r8
         prod(:,73) = 0._r8
         prod(:,104) = 0._r8
         prod(:,74) = 0._r8
         prod(:,106) = 0._r8
         prod(:,145) = 0._r8
         prod(:,169) = 0._r8
         prod(:,82) = + extfrc(:,16)
         prod(:,67) = 0._r8
         prod(:,84) = 0._r8
         prod(:,152) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,37) = 0._r8
         prod(:,9) = 0._r8
         prod(:,190) = + extfrc(:,2)
         prod(:,185) = + extfrc(:,3)
         prod(:,179) = 0._r8
         prod(:,141) = 0._r8
         prod(:,87) = 0._r8
         prod(:,10) = + extfrc(:,12)
         prod(:,11) = + extfrc(:,13)
         prod(:,12) = 0._r8
         prod(:,13) = + extfrc(:,14)
         prod(:,188) =.180_r8*rxt(:,40)*y(:,54) +rxt(:,41)*y(:,63) + (rxt(:,5) + &
                 2.000_r8*rxt(:,6))
         prod(:,186) = 0._r8
         prod(:,76) = 0._r8
         prod(:,80) = 0._r8
         prod(:,59) = 0._r8
         prod(:,97) = 0._r8
         prod(:,44) = 0._r8
         prod(:,98) = 0._r8
         prod(:,48) = 0._r8
         prod(:,77) = 0._r8
         prod(:,14) = + extfrc(:,6)
         prod(:,15) = + extfrc(:,7)
         prod(:,109) = 0._r8
         prod(:,88) = 0._r8
         prod(:,103) = 0._r8
         prod(:,167) = 0._r8
         prod(:,140) = + extfrc(:,4)
         prod(:,63) = 0._r8
         prod(:,16) = + extfrc(:,8)
         prod(:,17) = + extfrc(:,1)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) = 0._r8
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,29) = 0._r8
         prod(:,30) = 0._r8
         prod(:,31) = 0._r8
         prod(:,32) = 0._r8
         prod(:,33) = 0._r8
         prod(:,35) = + extfrc(:,5)
         prod(:,53) = 0._r8
         prod(:,115) = 0._r8
         prod(:,118) = 0._r8
         prod(:,99) = 0._r8
         prod(:,150) = 0._r8
         prod(:,153) = 0._r8
         prod(:,116) = 0._r8
         prod(:,51) = 0._r8
         prod(:,54) = 0._r8
         prod(:,55) = 0._r8
         prod(:,124) = 0._r8
         prod(:,56) = 0._r8
         prod(:,92) = 0._r8
         prod(:,105) = 0._r8
         prod(:,146) = 0._r8
         prod(:,100) = 0._r8
         prod(:,93) = 0._r8
         prod(:,137) = 0._r8
         prod(:,135) = 0._r8
         prod(:,119) = 0._r8
         prod(:,173) = 0._r8
         prod(:,182) =rxt(:,87)*y(:,43) +rxt(:,89)*y(:,46) +rxt(:,39)*y(:,54)
         prod(:,129) = 0._r8
         prod(:,112) = 0._r8
         prod(:,66) = 0._r8
         prod(:,125) = 0._r8
         prod(:,178) = 0._r8
         prod(:,94) = 0._r8
         prod(:,168) = 0._r8
         prod(:,171) = 0._r8
         prod(:,170) = 0._r8
         prod(:,126) = 0._r8
         prod(:,172) = 0._r8
         prod(:,143) = 0._r8
         prod(:,122) = 0._r8
         prod(:,159) = 0._r8
         prod(:,180) =rxt(:,12)*y(:,113) +rxt(:,5)
         prod(:,181) =.330_r8*rxt(:,40)*y(:,54) + extfrc(:,15)
         prod(:,95) = 0._r8
         prod(:,130) = 0._r8
         prod(:,160) = 0._r8
         prod(:,158) = 0._r8
         prod(:,157) = 0._r8
         prod(:,127) = 0._r8
         prod(:,163) = 0._r8
         prod(:,128) = 0._r8
         prod(:,101) = 0._r8
         prod(:,191) =.050_r8*rxt(:,40)*y(:,54)
      end if
      end subroutine indprd
      end module mo_indprd
