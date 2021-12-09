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
         prod(:,1) = + extfrc(:,16)
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
         prod(:,15) =.100_r8*rxt(:,328)*y(:,133)*y(:,29)
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = (rxt(:,285)*y(:,62) +rxt(:,287)*y(:,87) +rxt(:,295)*y(:,62) + &
                 rxt(:,315)*y(:,50) +.500_r8*rxt(:,316)*y(:,51) + &
                 .800_r8*rxt(:,321)*y(:,74) +rxt(:,322)*y(:,75) + &
                 .500_r8*rxt(:,371)*y(:,110) +1.800_r8*rxt(:,481)*y(:,209))*y(:,243) &
                  + (2.000_r8*rxt(:,311)*y(:,226) +.900_r8*rxt(:,312)*y(:,227) + &
                 rxt(:,314)*y(:,123) +2.000_r8*rxt(:,361)*y(:,238) + &
                 rxt(:,385)*y(:,234) +rxt(:,410)*y(:,250))*y(:,226) &
                  + (.200_r8*rxt(:,328)*y(:,29) +.100_r8*rxt(:,372)*y(:,112) + &
                 .270_r8*rxt(:,460)*y(:,6) +.270_r8*rxt(:,463)*y(:,111))*y(:,133) &
                  + (rxt(:,362)*y(:,227) +.450_r8*rxt(:,363)*y(:,232) + &
                 2.000_r8*rxt(:,364)*y(:,238))*y(:,238) &
                  + (.500_r8*rxt(:,470)*y(:,227) +.900_r8*rxt(:,472)*y(:,123)) &
                 *y(:,247) +rxt(:,37)*y(:,51) +.400_r8*rxt(:,60)*y(:,138) +rxt(:,65) &
                 *y(:,205) +.800_r8*rxt(:,69)*y(:,209)
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) =rxt(:,171)*y(:,124)*y(:,113)
         prod(:,25) = 0._r8
         prod(:,26) =rxt(:,499)*y(:,243)*y(:,121) +rxt(:,508)*y(:,122)
         prod(:,27) = (rxt(:,432)*y(:,228) +rxt(:,435)*y(:,237) +rxt(:,438)*y(:,239) + &
                 rxt(:,442)*y(:,140))*y(:,124) +.500_r8*rxt(:,371)*y(:,243)*y(:,110) &
                  +.200_r8*rxt(:,467)*y(:,241)*y(:,123) +.500_r8*rxt(:,479)*y(:,208) &
                 *y(:,125)
         prod(:,28) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,154) = 0._r8
         prod(:,155) = 0._r8
         prod(:,1) = + extfrc(:,1)
         prod(:,2) = + extfrc(:,2)
         prod(:,179) = 0._r8
         prod(:,79) = 0._r8
         prod(:,122) = 0._r8
         prod(:,80) = 0._r8
         prod(:,120) = 0._r8
         prod(:,130) = 0._r8
         prod(:,101) = 0._r8
         prod(:,151) = 0._r8
         prod(:,106) = 0._r8
         prod(:,92) = 0._r8
         prod(:,114) = 0._r8
         prod(:,208) =rxt(:,79)*y(:,34) +rxt(:,80)*y(:,35) +2.000_r8*rxt(:,86)*y(:,41) &
                  +rxt(:,87)*y(:,43) +3.000_r8*rxt(:,90)*y(:,55) +2.000_r8*rxt(:,98) &
                 *y(:,78)
         prod(:,93) = 0._r8
         prod(:,221) = 0._r8
         prod(:,146) = 0._r8
         prod(:,94) = 0._r8
         prod(:,109) = 0._r8
         prod(:,103) = 0._r8
         prod(:,144) = 0._r8
         prod(:,98) = 0._r8
         prod(:,110) = 0._r8
         prod(:,104) = 0._r8
         prod(:,184) = 0._r8
         prod(:,121) = 0._r8
         prod(:,73) = 0._r8
         prod(:,99) = 0._r8
         prod(:,216) =.180_r8*rxt(:,39)*y(:,54)
         prod(:,195) = 0._r8
         prod(:,71) = 0._r8
         prod(:,181) = 0._r8
         prod(:,200) = 0._r8
         prod(:,141) = 0._r8
         prod(:,136) = 0._r8
         prod(:,167) = 0._r8
         prod(:,123) = 0._r8
         prod(:,224) =4.000_r8*rxt(:,78)*y(:,33) +rxt(:,79)*y(:,34) &
                  +2.000_r8*rxt(:,81)*y(:,36) +2.000_r8*rxt(:,82)*y(:,37) &
                  +2.000_r8*rxt(:,83)*y(:,38) +rxt(:,84)*y(:,39) +2.000_r8*rxt(:,85) &
                 *y(:,40) +3.000_r8*rxt(:,88)*y(:,44) +rxt(:,89)*y(:,46) +rxt(:,100) &
                 *y(:,82) +rxt(:,101)*y(:,83) +rxt(:,102)*y(:,84)
         prod(:,82) = 0._r8
         prod(:,72) = 0._r8
         prod(:,222) = 0._r8
         prod(:,182) = 0._r8
         prod(:,189) =.380_r8*rxt(:,39)*y(:,54) +rxt(:,41)*y(:,63) + extfrc(:,3)
         prod(:,76) =rxt(:,79)*y(:,34) +rxt(:,80)*y(:,35) +rxt(:,82)*y(:,37) &
                  +2.000_r8*rxt(:,83)*y(:,38) +2.000_r8*rxt(:,84)*y(:,39) +rxt(:,85) &
                 *y(:,40) +2.000_r8*rxt(:,98)*y(:,78) +rxt(:,101)*y(:,83) +rxt(:,102) &
                 *y(:,84)
         prod(:,83) =rxt(:,81)*y(:,36) +rxt(:,82)*y(:,37) +rxt(:,100)*y(:,82)
         prod(:,86) = 0._r8
         prod(:,105) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,77) = 0._r8
         prod(:,166) =rxt(:,80)*y(:,35) +rxt(:,84)*y(:,39)
         prod(:,186) = 0._r8
         prod(:,177) = 0._r8
         prod(:,210) = (.330_r8*rxt(:,39) +rxt(:,40))*y(:,54)
         prod(:,196) =1.440_r8*rxt(:,39)*y(:,54)
         prod(:,148) = 0._r8
         prod(:,78) = 0._r8
         prod(:,170) = 0._r8
         prod(:,212) = 0._r8
         prod(:,84) = 0._r8
         prod(:,168) = 0._r8
         prod(:,96) = 0._r8
         prod(:,209) = 0._r8
         prod(:,124) = 0._r8
         prod(:,165) = 0._r8
         prod(:,171) = 0._r8
         prod(:,187) = 0._r8
         prod(:,95) = 0._r8
         prod(:,190) = 0._r8
         prod(:,111) = 0._r8
         prod(:,74) = 0._r8
         prod(:,175) = 0._r8
         prod(:,147) = 0._r8
         prod(:,142) = 0._r8
         prod(:,198) = 0._r8
         prod(:,118) = 0._r8
         prod(:,157) = 0._r8
         prod(:,56) = 0._r8
         prod(:,67) = 0._r8
         prod(:,199) = 0._r8
         prod(:,113) = 0._r8
         prod(:,138) = 0._r8
         prod(:,112) = 0._r8
         prod(:,143) = 0._r8
         prod(:,178) = 0._r8
         prod(:,203) = 0._r8
         prod(:,125) = + extfrc(:,17)
         prod(:,102) = 0._r8
         prod(:,116) = 0._r8
         prod(:,185) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,70) = 0._r8
         prod(:,9) = 0._r8
         prod(:,211) = + extfrc(:,4)
         prod(:,217) = + extfrc(:,5)
         prod(:,218) = 0._r8
         prod(:,173) = 0._r8
         prod(:,119) = 0._r8
         prod(:,10) = + extfrc(:,6)
         prod(:,11) = + extfrc(:,7)
         prod(:,12) = 0._r8
         prod(:,13) = + extfrc(:,8)
         prod(:,223) =.180_r8*rxt(:,39)*y(:,54) +rxt(:,41)*y(:,63) + (rxt(:,5) + &
                 2.000_r8*rxt(:,6))
         prod(:,220) = 0._r8
         prod(:,107) = 0._r8
         prod(:,115) = 0._r8
         prod(:,91) = 0._r8
         prod(:,132) = 0._r8
         prod(:,75) = 0._r8
         prod(:,131) = 0._r8
         prod(:,81) = 0._r8
         prod(:,108) = 0._r8
         prod(:,14) = 0._r8
         prod(:,15) = + extfrc(:,15)
         prod(:,16) = 0._r8
         prod(:,17) = + extfrc(:,14)
         prod(:,140) = 0._r8
         prod(:,117) = 0._r8
         prod(:,137) = 0._r8
         prod(:,201) = 0._r8
         prod(:,174) = + extfrc(:,9)
         prod(:,97) = 0._r8
         prod(:,18) = + extfrc(:,10)
         prod(:,19) = + extfrc(:,11)
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
         prod(:,34) = 0._r8
         prod(:,35) = 0._r8
         prod(:,36) = 0._r8
         prod(:,37) = 0._r8
         prod(:,38) = 0._r8
         prod(:,39) = 0._r8
         prod(:,40) = 0._r8
         prod(:,41) = 0._r8
         prod(:,42) = 0._r8
         prod(:,43) = 0._r8
         prod(:,44) = 0._r8
         prod(:,45) = 0._r8
         prod(:,46) = 0._r8
         prod(:,47) = 0._r8
         prod(:,48) = 0._r8
         prod(:,49) = 0._r8
         prod(:,50) = 0._r8
         prod(:,51) = 0._r8
         prod(:,52) = 0._r8
         prod(:,53) = 0._r8
         prod(:,54) = 0._r8
         prod(:,55) = 0._r8
         prod(:,57) = 0._r8
         prod(:,58) = 0._r8
         prod(:,59) = 0._r8
         prod(:,60) = 0._r8
         prod(:,61) = 0._r8
         prod(:,62) = 0._r8
         prod(:,63) = 0._r8
         prod(:,64) = 0._r8
         prod(:,65) = 0._r8
         prod(:,66) = 0._r8
         prod(:,68) = + extfrc(:,12)
         prod(:,69) = + extfrc(:,13)
         prod(:,87) = 0._r8
         prod(:,149) = 0._r8
         prod(:,152) = 0._r8
         prod(:,133) = 0._r8
         prod(:,183) = 0._r8
         prod(:,188) = 0._r8
         prod(:,150) = 0._r8
         prod(:,85) = 0._r8
         prod(:,88) = 0._r8
         prod(:,89) = 0._r8
         prod(:,158) = 0._r8
         prod(:,90) = 0._r8
         prod(:,126) = 0._r8
         prod(:,139) = 0._r8
         prod(:,180) = 0._r8
         prod(:,134) = 0._r8
         prod(:,127) = 0._r8
         prod(:,172) = 0._r8
         prod(:,169) = 0._r8
         prod(:,153) = 0._r8
         prod(:,207) = 0._r8
         prod(:,219) =rxt(:,87)*y(:,43) +rxt(:,89)*y(:,46) +rxt(:,40)*y(:,54)
         prod(:,163) = 0._r8
         prod(:,145) = 0._r8
         prod(:,100) = 0._r8
         prod(:,159) = 0._r8
         prod(:,213) = 0._r8
         prod(:,128) = 0._r8
         prod(:,202) = 0._r8
         prod(:,205) = 0._r8
         prod(:,204) = 0._r8
         prod(:,160) = 0._r8
         prod(:,206) = 0._r8
         prod(:,176) = 0._r8
         prod(:,156) = 0._r8
         prod(:,193) = 0._r8
         prod(:,214) =rxt(:,12)*y(:,114) +rxt(:,5)
         prod(:,215) =.330_r8*rxt(:,39)*y(:,54) + extfrc(:,18)
         prod(:,129) = 0._r8
         prod(:,164) = 0._r8
         prod(:,194) = 0._r8
         prod(:,192) = 0._r8
         prod(:,191) = 0._r8
         prod(:,161) = 0._r8
         prod(:,197) = 0._r8
         prod(:,162) = 0._r8
         prod(:,135) = 0._r8
         prod(:,225) =.050_r8*rxt(:,39)*y(:,54)
      end if
      end subroutine indprd
      end module mo_indprd
