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
         prod(:,1) = + extfrc(:,15)
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
         prod(:,15) =.100_r8*rxt(:,348)*y(:,135)*y(:,29)
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = (rxt(:,305)*y(:,62) +rxt(:,307)*y(:,87) +rxt(:,315)*y(:,62) + &
                 rxt(:,335)*y(:,50) +.500_r8*rxt(:,336)*y(:,51) + &
                 .800_r8*rxt(:,341)*y(:,74) +rxt(:,342)*y(:,75) + &
                 .500_r8*rxt(:,391)*y(:,109) +1.800_r8*rxt(:,501)*y(:,179))*y(:,221) &
                  + (2.000_r8*rxt(:,331)*y(:,196) +.900_r8*rxt(:,332)*y(:,197) + &
                 rxt(:,334)*y(:,124) +2.000_r8*rxt(:,381)*y(:,209) + &
                 rxt(:,405)*y(:,205) +rxt(:,430)*y(:,229))*y(:,196) &
                  + (.200_r8*rxt(:,348)*y(:,29) +.100_r8*rxt(:,392)*y(:,111) + &
                 .270_r8*rxt(:,480)*y(:,6) +.270_r8*rxt(:,483)*y(:,110))*y(:,135) &
                  + (rxt(:,382)*y(:,197) +.450_r8*rxt(:,383)*y(:,203) + &
                 2.000_r8*rxt(:,384)*y(:,209))*y(:,209) &
                  + (.500_r8*rxt(:,490)*y(:,197) +.900_r8*rxt(:,492)*y(:,124)) &
                 *y(:,226) +rxt(:,38)*y(:,51) +.400_r8*rxt(:,61)*y(:,140) +rxt(:,66) &
                 *y(:,175) +.800_r8*rxt(:,70)*y(:,179)
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) =rxt(:,191)*y(:,125)*y(:,112)
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,29) = 0._r8
         prod(:,30) =rxt(:,518)*y(:,221)*y(:,120) +rxt(:,528)*y(:,121)
         prod(:,31) = (rxt(:,452)*y(:,198) +rxt(:,455)*y(:,208) +rxt(:,458)*y(:,210) + &
                 rxt(:,462)*y(:,142))*y(:,125) +.500_r8*rxt(:,391)*y(:,221)*y(:,109) &
                  +.200_r8*rxt(:,487)*y(:,216)*y(:,124) +.500_r8*rxt(:,499)*y(:,178) &
                 *y(:,126)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,123) = 0._r8
         prod(:,124) = 0._r8
         prod(:,1) = + extfrc(:,13)
         prod(:,2) = + extfrc(:,14)
         prod(:,153) = 0._r8
         prod(:,48) = 0._r8
         prod(:,84) = 0._r8
         prod(:,49) = 0._r8
         prod(:,85) = 0._r8
         prod(:,95) = 0._r8
         prod(:,70) = 0._r8
         prod(:,118) = 0._r8
         prod(:,76) = 0._r8
         prod(:,62) = 0._r8
         prod(:,82) = 0._r8
         prod(:,184) =rxt(:,80)*y(:,34) +rxt(:,81)*y(:,35) +2.000_r8*rxt(:,87)*y(:,41) &
                  +rxt(:,88)*y(:,43) +3.000_r8*rxt(:,91)*y(:,55) +2.000_r8*rxt(:,99) &
                 *y(:,78)
         prod(:,63) = 0._r8
         prod(:,198) = 0._r8
         prod(:,110) = 0._r8
         prod(:,64) = 0._r8
         prod(:,79) = 0._r8
         prod(:,71) = 0._r8
         prod(:,112) = 0._r8
         prod(:,66) = 0._r8
         prod(:,80) = 0._r8
         prod(:,72) = 0._r8
         prod(:,160) = 0._r8
         prod(:,89) = 0._r8
         prod(:,39) = 0._r8
         prod(:,67) = 0._r8
         prod(:,193) =.180_r8*rxt(:,41)*y(:,54)
         prod(:,170) = 0._r8
         prod(:,38) = 0._r8
         prod(:,156) = 0._r8
         prod(:,175) = 0._r8
         prod(:,111) = 0._r8
         prod(:,105) = 0._r8
         prod(:,140) = 0._r8
         prod(:,90) = 0._r8
         prod(:,188) =4.000_r8*rxt(:,79)*y(:,33) +rxt(:,80)*y(:,34) &
                  +2.000_r8*rxt(:,82)*y(:,36) +2.000_r8*rxt(:,83)*y(:,37) &
                  +2.000_r8*rxt(:,84)*y(:,38) +rxt(:,85)*y(:,39) +2.000_r8*rxt(:,86) &
                 *y(:,40) +3.000_r8*rxt(:,89)*y(:,44) +rxt(:,90)*y(:,46) +rxt(:,101) &
                 *y(:,82) +rxt(:,102)*y(:,83) +rxt(:,103)*y(:,84)
         prod(:,47) = 0._r8
         prod(:,36) = 0._r8
         prod(:,200) = 0._r8
         prod(:,157) = 0._r8
         prod(:,165) = (rxt(:,42) +rxt(:,110))*y(:,63) +.380_r8*rxt(:,41)*y(:,54) &
                  + extfrc(:,3)
         prod(:,40) =rxt(:,80)*y(:,34) +rxt(:,81)*y(:,35) +rxt(:,83)*y(:,37) &
                  +2.000_r8*rxt(:,84)*y(:,38) +2.000_r8*rxt(:,85)*y(:,39) +rxt(:,86) &
                 *y(:,40) +2.000_r8*rxt(:,99)*y(:,78) +rxt(:,102)*y(:,83) +rxt(:,103) &
                 *y(:,84)
         prod(:,51) =rxt(:,82)*y(:,36) +rxt(:,83)*y(:,37) +rxt(:,101)*y(:,82)
         prod(:,54) = 0._r8
         prod(:,69) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,41) = 0._r8
         prod(:,136) =rxt(:,81)*y(:,35) +rxt(:,85)*y(:,39)
         prod(:,161) = 0._r8
         prod(:,149) = 0._r8
         prod(:,195) = (rxt(:,40) +.330_r8*rxt(:,41))*y(:,54)
         prod(:,172) =1.440_r8*rxt(:,41)*y(:,54)
         prod(:,115) = 0._r8
         prod(:,42) = 0._r8
         prod(:,145) = 0._r8
         prod(:,183) = 0._r8
         prod(:,52) = 0._r8
         prod(:,141) = 0._r8
         prod(:,59) = 0._r8
         prod(:,196) = 0._r8
         prod(:,99) = 0._r8
         prod(:,134) = 0._r8
         prod(:,146) = 0._r8
         prod(:,162) = 0._r8
         prod(:,60) = 0._r8
         prod(:,164) = 0._r8
         prod(:,73) = 0._r8
         prod(:,43) = 0._r8
         prod(:,148) = 0._r8
         prod(:,119) = 0._r8
         prod(:,108) = 0._r8
         prod(:,173) = 0._r8
         prod(:,88) = 0._r8
         prod(:,127) = 0._r8
         prod(:,34) = 0._r8
         prod(:,174) = 0._r8
         prod(:,74) = 0._r8
         prod(:,107) = 0._r8
         prod(:,75) = 0._r8
         prod(:,114) = 0._r8
         prod(:,151) = 0._r8
         prod(:,179) = 0._r8
         prod(:,144) = (.800_r8*rxt(:,112) +rxt(:,115) +rxt(:,116) + &
                 .800_r8*rxt(:,118)) + extfrc(:,21)
         prod(:,68) = 0._r8
         prod(:,83) = 0._r8
         prod(:,159) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,37) = 0._r8
         prod(:,9) = 0._r8
         prod(:,187) = + extfrc(:,2)
         prod(:,197) = + extfrc(:,1)
         prod(:,194) = 0._r8
         prod(:,147) = 0._r8
         prod(:,86) = 0._r8
         prod(:,10) = + extfrc(:,10)
         prod(:,11) = + extfrc(:,11)
         prod(:,12) = 0._r8
         prod(:,13) = + extfrc(:,12)
         prod(:,192) = (rxt(:,42) +rxt(:,110))*y(:,63) +.180_r8*rxt(:,41)*y(:,54)
         prod(:,186) = 0._r8
         prod(:,191) = 0._r8
         prod(:,77) = 0._r8
         prod(:,81) = 0._r8
         prod(:,61) = 0._r8
         prod(:,97) = 0._r8
         prod(:,44) = 0._r8
         prod(:,98) = 0._r8
         prod(:,50) = 0._r8
         prod(:,78) = 0._r8
         prod(:,14) = + extfrc(:,8)
         prod(:,15) = + extfrc(:,9)
         prod(:,109) = 0._r8
         prod(:,87) = 0._r8
         prod(:,129) = 0._r8
         prod(:,182) = 0._r8
         prod(:,155) = + extfrc(:,4)
         prod(:,65) = 0._r8
         prod(:,16) = + extfrc(:,6)
         prod(:,17) = + extfrc(:,7)
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
         prod(:,55) = 0._r8
         prod(:,116) = 0._r8
         prod(:,121) = 0._r8
         prod(:,100) = 0._r8
         prod(:,158) = 0._r8
         prod(:,163) = 0._r8
         prod(:,117) = 0._r8
         prod(:,53) = 0._r8
         prod(:,56) = 0._r8
         prod(:,57) = 0._r8
         prod(:,128) = 0._r8
         prod(:,58) = 0._r8
         prod(:,91) = 0._r8
         prod(:,104) = 0._r8
         prod(:,154) = 0._r8
         prod(:,101) = 0._r8
         prod(:,92) = 0._r8
         prod(:,152) = 0._r8
         prod(:,143) = 0._r8
         prod(:,122) = 0._r8
         prod(:,181) = 0._r8
         prod(:,185) =rxt(:,88)*y(:,43) +rxt(:,90)*y(:,46) +rxt(:,40)*y(:,54)
         prod(:,133) = 0._r8
         prod(:,139) = (rxt(:,113) +rxt(:,114) +rxt(:,115) +rxt(:,116) +rxt(:,117) + &
                 rxt(:,119)) + extfrc(:,20)
         prod(:,113) = 0._r8
         prod(:,96) = 0._r8
         prod(:,135) = 0._r8
         prod(:,199) = 0._r8
         prod(:,93) = 0._r8
         prod(:,176) = 0._r8
         prod(:,177) = 0._r8
         prod(:,178) = 0._r8
         prod(:,130) = 0._r8
         prod(:,180) = 0._r8
         prod(:,150) = 0._r8
         prod(:,125) = 0._r8
         prod(:,106) = (1.200_r8*rxt(:,112) +rxt(:,113) +rxt(:,117) + &
                 1.200_r8*rxt(:,118)) + extfrc(:,19)
         prod(:,126) = (rxt(:,114) +rxt(:,119)) + extfrc(:,18)
         prod(:,138) = 0._r8
         prod(:,102) = (rxt(:,113) +rxt(:,115) +rxt(:,116) +rxt(:,117)) + extfrc(:,17)
         prod(:,168) = 0._r8
         prod(:,189) =rxt(:,12)*y(:,113)
         prod(:,45) = 0._r8
         prod(:,46) = 0._r8
         prod(:,137) = + extfrc(:,16)
         prod(:,190) =.330_r8*rxt(:,41)*y(:,54) + extfrc(:,22)
         prod(:,120) = + extfrc(:,23)
         prod(:,94) = 0._r8
         prod(:,142) = 0._r8
         prod(:,169) = 0._r8
         prod(:,167) = 0._r8
         prod(:,166) = 0._r8
         prod(:,131) = 0._r8
         prod(:,171) = 0._r8
         prod(:,132) = 0._r8
         prod(:,103) = 0._r8
         prod(:,201) =.050_r8*rxt(:,41)*y(:,54)
      end if
      end subroutine indprd
      end module mo_indprd
