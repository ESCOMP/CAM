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
         prod(:,1) =rxt(:,594)*y(:,250)*y(:,149) +rxt(:,604)*y(:,150)
         prod(:,2) = (rxt(:,506)*y(:,231) +rxt(:,509)*y(:,243) +rxt(:,512)*y(:,245) + &
                 rxt(:,516)*y(:,172))*y(:,154) +.500_r8*rxt(:,445)*y(:,250)*y(:,138) &
                  +.200_r8*rxt(:,541)*y(:,248)*y(:,153) +.500_r8*rxt(:,553)*y(:,209) &
                 *y(:,155)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,176) = 0._r8
         prod(:,178) = 0._r8
         prod(:,1) = + extfrc(:,18)
         prod(:,2) = + extfrc(:,11)
         prod(:,3) = + extfrc(:,12)
         prod(:,200) = 0._r8
         prod(:,81) = 0._r8
         prod(:,135) = 0._r8
         prod(:,82) = 0._r8
         prod(:,138) = 0._r8
         prod(:,147) = 0._r8
         prod(:,110) = 0._r8
         prod(:,172) = 0._r8
         prod(:,118) = 0._r8
         prod(:,100) = 0._r8
         prod(:,127) = 0._r8
         prod(:,257) = 0._r8
         prod(:,148) = 0._r8
         prod(:,159) = 0._r8
         prod(:,119) = 0._r8
         prod(:,249) = 0._r8
         prod(:,194) = 0._r8
         prod(:,4) = 0._r8
         prod(:,101) = 0._r8
         prod(:,120) = 0._r8
         prod(:,87) = 0._r8
         prod(:,111) = 0._r8
         prod(:,165) = 0._r8
         prod(:,103) = 0._r8
         prod(:,105) = 0._r8
         prod(:,122) = 0._r8
         prod(:,113) = 0._r8
         prod(:,211) = 0._r8
         prod(:,137) = 0._r8
         prod(:,114) = 0._r8
         prod(:,106) = 0._r8
         prod(:,61) = 0._r8
         prod(:,77) = 0._r8
         prod(:,78) = 0._r8
         prod(:,67) = 0._r8
         prod(:,79) = 0._r8
         prod(:,68) = 0._r8
         prod(:,80) = 0._r8
         prod(:,69) = 0._r8
         prod(:,151) = 0._r8
         prod(:,70) = 0._r8
         prod(:,115) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,241) = 0._r8
         prod(:,173) = 0._r8
         prod(:,71) = 0._r8
         prod(:,215) = 0._r8
         prod(:,131) = 0._r8
         prod(:,62) = 0._r8
         prod(:,206) = 0._r8
         prod(:,226) = 0._r8
         prod(:,180) = 0._r8
         prod(:,167) = 0._r8
         prod(:,130) = 0._r8
         prod(:,189) = 0._r8
         prod(:,140) = 0._r8
         prod(:,237) = 0._r8
         prod(:,63) = 0._r8
         prod(:,132) = 0._r8
         prod(:,64) = 0._r8
         prod(:,210) = 0._r8
         prod(:,125) = 0._r8
         prod(:,252) = 0._r8
         prod(:,92) = 0._r8
         prod(:,59) = 0._r8
         prod(:,161) = 0._r8
         prod(:,250) = 0._r8
         prod(:,207) = 0._r8
         prod(:,8) = 0._r8
         prod(:,217) = + extfrc(:,4)
         prod(:,190) = 0._r8
         prod(:,124) = 0._r8
         prod(:,98) = 0._r8
         prod(:,96) = 0._r8
         prod(:,89) = 0._r8
         prod(:,208) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,12) = 0._r8
         prod(:,72) = 0._r8
         prod(:,197) = 0._r8
         prod(:,216) = 0._r8
         prod(:,203) = 0._r8
         prod(:,255) = 0._r8
         prod(:,234) = 0._r8
         prod(:,65) = 0._r8
         prod(:,170) = 0._r8
         prod(:,73) = 0._r8
         prod(:,233) = 0._r8
         prod(:,95) = 0._r8
         prod(:,97) = 0._r8
         prod(:,109) = 0._r8
         prod(:,242) = 0._r8
         prod(:,84) = 0._r8
         prod(:,202) = 0._r8
         prod(:,112) = 0._r8
         prod(:,221) = 0._r8
         prod(:,244) = 0._r8
         prod(:,143) = 0._r8
         prod(:,238) = 0._r8
         prod(:,243) = 0._r8
         prod(:,240) = 0._r8
         prod(:,219) = 0._r8
         prod(:,99) = 0._r8
         prod(:,220) = 0._r8
         prod(:,123) = 0._r8
         prod(:,254) = 0._r8
         prod(:,235) = 0._r8
         prod(:,66) = 0._r8
         prod(:,74) = 0._r8
         prod(:,58) = 0._r8
         prod(:,128) = 0._r8
         prod(:,150) = 0._r8
         prod(:,75) = 0._r8
         prod(:,85) = 0._r8
         prod(:,108) = 0._r8
         prod(:,253) = 0._r8
         prod(:,196) = 0._r8
         prod(:,192) = 0._r8
         prod(:,168) = 0._r8
         prod(:,162) = 0._r8
         prod(:,139) = 0._r8
         prod(:,185) = 0._r8
         prod(:,54) = 0._r8
         prod(:,13) = 0._r8
         prod(:,225) = 0._r8
         prod(:,116) = 0._r8
         prod(:,157) = 0._r8
         prod(:,117) = 0._r8
         prod(:,169) = 0._r8
         prod(:,199) = 0._r8
         prod(:,230) = 0._r8
         prod(:,153) = + extfrc(:,17)
         prod(:,86) = 0._r8
         prod(:,126) = 0._r8
         prod(:,133) = 0._r8
         prod(:,214) = 0._r8
         prod(:,14) = 0._r8
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,60) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,246) = + extfrc(:,5)
         prod(:,245) = + extfrc(:,6)
         prod(:,259) = 0._r8
         prod(:,193) = 0._r8
         prod(:,134) = 0._r8
         prod(:,20) = + extfrc(:,13)
         prod(:,21) = + extfrc(:,14)
         prod(:,22) = 0._r8
         prod(:,23) = + extfrc(:,15)
         prod(:,24) = + extfrc(:,16)
         prod(:,256) = (rxt(:,5) +2.000_r8*rxt(:,6))
         prod(:,258) = 0._r8
         prod(:,25) = 0._r8
         prod(:,146) = 0._r8
         prod(:,129) = 0._r8
         prod(:,156) = 0._r8
         prod(:,102) = 0._r8
         prod(:,163) = 0._r8
         prod(:,76) = 0._r8
         prod(:,149) = 0._r8
         prod(:,83) = 0._r8
         prod(:,121) = 0._r8
         prod(:,26) = + extfrc(:,9)
         prod(:,27) = + extfrc(:,10)
         prod(:,164) = 0._r8
         prod(:,136) = 0._r8
         prod(:,158) = 0._r8
         prod(:,28) = 0._r8
         prod(:,236) = 0._r8
         prod(:,209) = + extfrc(:,7)
         prod(:,104) = 0._r8
         prod(:,29) = + extfrc(:,1)
         prod(:,30) = + extfrc(:,2)
         prod(:,31) = 0._r8
         prod(:,32) = + extfrc(:,3)
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
         prod(:,49) = + extfrc(:,8)
         prod(:,90) = 0._r8
         prod(:,177) = 0._r8
         prod(:,174) = 0._r8
         prod(:,152) = 0._r8
         prod(:,205) = 0._r8
         prod(:,213) = 0._r8
         prod(:,171) = 0._r8
         prod(:,88) = 0._r8
         prod(:,91) = 0._r8
         prod(:,93) = 0._r8
         prod(:,181) = 0._r8
         prod(:,94) = 0._r8
         prod(:,141) = 0._r8
         prod(:,160) = 0._r8
         prod(:,201) = 0._r8
         prod(:,50) = 0._r8
         prod(:,154) = 0._r8
         prod(:,51) = 0._r8
         prod(:,142) = 0._r8
         prod(:,191) = 0._r8
         prod(:,195) = 0._r8
         prod(:,175) = 0._r8
         prod(:,232) = 0._r8
         prod(:,239) = 0._r8
         prod(:,187) = 0._r8
         prod(:,166) = 0._r8
         prod(:,107) = 0._r8
         prod(:,182) = 0._r8
         prod(:,251) = 0._r8
         prod(:,144) = 0._r8
         prod(:,227) = 0._r8
         prod(:,228) = 0._r8
         prod(:,224) = 0._r8
         prod(:,52) = 0._r8
         prod(:,53) = 0._r8
         prod(:,229) = 0._r8
         prod(:,183) = 0._r8
         prod(:,231) = 0._r8
         prod(:,198) = 0._r8
         prod(:,179) = 0._r8
         prod(:,55) = 0._r8
         prod(:,212) = 0._r8
         prod(:,247) =rxt(:,5)
         prod(:,248) = 0._r8
         prod(:,145) = 0._r8
         prod(:,188) = 0._r8
         prod(:,222) = 0._r8
         prod(:,218) = 0._r8
         prod(:,204) = 0._r8
         prod(:,184) = 0._r8
         prod(:,56) = 0._r8
         prod(:,223) = 0._r8
         prod(:,186) = 0._r8
         prod(:,57) = 0._r8
         prod(:,155) = 0._r8
         prod(:,260) = 0._r8
      end if
      end subroutine indprd
      end module mo_indprd
