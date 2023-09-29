      module mo_prod_loss
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : veclen
      private
      public :: exp_prod_loss
      public :: imp_prod_loss
      contains
      subroutine exp_prod_loss( ofl, ofu, prod, loss, y, &
                                rxt, het_rates, chnkpnts )
      use chem_mods, only : gas_pcnst,rxntot,clscnt1
      implicit none
!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: ofl, ofu, chnkpnts
      real(r8), dimension(chnkpnts,max(1,clscnt1)), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(in) :: het_rates(chnkpnts,gas_pcnst)
!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      integer :: k
!--------------------------------------------------------------------
! ... loss and production for Explicit method
!--------------------------------------------------------------------
      do k = ofl,ofu
         loss(k,1) = ( + het_rates(k,233))* y(k,233)
         prod(k,1) = 0._r8
         loss(k,2) = ( + het_rates(k,234))* y(k,234)
         prod(k,2) = 0._r8
      end do
      end subroutine exp_prod_loss
      subroutine imp_prod_loss( avec_len, prod, loss, y, &
                                rxt, het_rates )
      use chem_mods, only : gas_pcnst,rxntot,clscnt4
      implicit none
!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), dimension(veclen,clscnt4), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      integer :: k
!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------
      do k = 1,avec_len
         loss(k,168) = (rxt(k,408)* y(k,295) + rxt(k,19) + het_rates(k,1))* y(k,1)
         prod(k,168) =rxt(k,411)*y(k,236)*y(k,147)
         loss(k,165) = (rxt(k,412)* y(k,295) + rxt(k,20) + het_rates(k,2))* y(k,2)
         prod(k,165) =rxt(k,409)*y(k,258)*y(k,236)
         loss(k,1) = ( + het_rates(k,3))* y(k,3)
         prod(k,1) = 0._r8
         loss(k,219) = (rxt(k,584)* y(k,149) +rxt(k,602)* y(k,158) +rxt(k,603) &
                 * y(k,295) + het_rates(k,4))* y(k,4)
         prod(k,219) = 0._r8
         loss(k,2) = ( + het_rates(k,5))* y(k,5)
         prod(k,2) = 0._r8
         loss(k,3) = ( + het_rates(k,6))* y(k,6)
         prod(k,3) = 0._r8
         loss(k,199) = (rxt(k,604)* y(k,149) +rxt(k,622)* y(k,158) +rxt(k,623) &
                 * y(k,295) + het_rates(k,7))* y(k,7)
         prod(k,199) = 0._r8
         loss(k,76) = (rxt(k,543)* y(k,295) + het_rates(k,8))* y(k,8)
         prod(k,76) = 0._r8
         loss(k,124) = (rxt(k,546)* y(k,295) + rxt(k,21) + het_rates(k,9))* y(k,9)
         prod(k,124) =rxt(k,544)*y(k,258)*y(k,243)
         loss(k,77) = ( + rxt(k,22) + het_rates(k,10))* y(k,10)
         prod(k,77) =.120_r8*rxt(k,543)*y(k,295)*y(k,8)
         loss(k,132) = ( + rxt(k,23) + het_rates(k,11))* y(k,11)
         prod(k,132) = (.500_r8*rxt(k,545)*y(k,243) +.200_r8*rxt(k,572)*y(k,314) + &
                 .060_r8*rxt(k,578)*y(k,316))*y(k,147) +.500_r8*rxt(k,21)*y(k,9) &
                  +rxt(k,22)*y(k,10) +.200_r8*rxt(k,115)*y(k,227) +.060_r8*rxt(k,116) &
                 *y(k,230)
         loss(k,103) = ( + rxt(k,24) + het_rates(k,12))* y(k,12)
         prod(k,103) = (.200_r8*rxt(k,572)*y(k,314) +.200_r8*rxt(k,578)*y(k,316)) &
                 *y(k,147) +.200_r8*rxt(k,115)*y(k,227) +.200_r8*rxt(k,116)*y(k,230)
         loss(k,121) = ( + rxt(k,25) + het_rates(k,13))* y(k,13)
         prod(k,121) = (.200_r8*rxt(k,572)*y(k,314) +.150_r8*rxt(k,578)*y(k,316)) &
                 *y(k,147) +.200_r8*rxt(k,115)*y(k,227) +.150_r8*rxt(k,116)*y(k,230)
         loss(k,110) = ( + rxt(k,26) + het_rates(k,14))* y(k,14)
         prod(k,110) =.210_r8*rxt(k,578)*y(k,316)*y(k,147) +.210_r8*rxt(k,116) &
                 *y(k,230)
         loss(k,81) = (rxt(k,413)* y(k,295) + het_rates(k,15))* y(k,15)
         prod(k,81) =.190_r8*rxt(k,642)*y(k,158)*y(k,17)
         loss(k,118) = (rxt(k,374)* y(k,149) +rxt(k,375)* y(k,295) + het_rates(k,16)) &
                 * y(k,16)
         prod(k,118) = 0._r8
         loss(k,202) = (rxt(k,624)* y(k,149) +rxt(k,642)* y(k,158) +rxt(k,643) &
                 * y(k,295) + het_rates(k,17))* y(k,17)
         prod(k,202) = 0._r8
         loss(k,276) = (rxt(k,254)* y(k,43) +rxt(k,256)* y(k,158) +rxt(k,255) &
                 * y(k,258) + het_rates(k,18))* y(k,18)
         prod(k,276) = (rxt(k,119) +2.000_r8*rxt(k,257)*y(k,20) +rxt(k,258)*y(k,60) + &
                 rxt(k,259)*y(k,60) +rxt(k,262)*y(k,147) +rxt(k,265)*y(k,157) + &
                 rxt(k,266)*y(k,295) +rxt(k,800)*y(k,174))*y(k,20) &
                  + (rxt(k,244)*y(k,35) +rxt(k,270)*y(k,36) + &
                 3.000_r8*rxt(k,271)*y(k,56) +2.000_r8*rxt(k,272)*y(k,80) + &
                 rxt(k,273)*y(k,83) +2.000_r8*rxt(k,293)*y(k,42) +rxt(k,294)*y(k,44)) &
                 *y(k,294) + (rxt(k,268)*y(k,83) +2.000_r8*rxt(k,282)*y(k,42) + &
                 rxt(k,284)*y(k,44) +3.000_r8*rxt(k,289)*y(k,56))*y(k,295) &
                  + (2.000_r8*rxt(k,281)*y(k,42) +rxt(k,283)*y(k,44) + &
                 3.000_r8*rxt(k,288)*y(k,56))*y(k,57) + (rxt(k,143) + &
                 rxt(k,267)*y(k,157))*y(k,83) +rxt(k,118)*y(k,19) +rxt(k,121)*y(k,21) &
                  +rxt(k,123)*y(k,35) +rxt(k,124)*y(k,36) +2.000_r8*rxt(k,130)*y(k,42) &
                  +rxt(k,131)*y(k,44) +3.000_r8*rxt(k,134)*y(k,56) &
                  +2.000_r8*rxt(k,142)*y(k,80) +rxt(k,149)*y(k,95)
         loss(k,88) = ( + rxt(k,118) + het_rates(k,19))* y(k,19)
         prod(k,88) = (rxt(k,892)*y(k,95) +rxt(k,897)*y(k,95))*y(k,87) &
                  +rxt(k,260)*y(k,60)*y(k,20)
         loss(k,309) = (2._r8*rxt(k,257)* y(k,20) + (rxt(k,258) +rxt(k,259) + &
                 rxt(k,260))* y(k,60) +rxt(k,262)* y(k,147) +rxt(k,263)* y(k,148) &
                  +rxt(k,265)* y(k,157) +rxt(k,800)* y(k,174) +rxt(k,261)* y(k,258) &
                  +rxt(k,266)* y(k,295) + rxt(k,119) + het_rates(k,20))* y(k,20)
         prod(k,309) = (rxt(k,120) +rxt(k,264)*y(k,157))*y(k,21) +rxt(k,256)*y(k,158) &
                 *y(k,18) +rxt(k,274)*y(k,294)*y(k,83) +rxt(k,269)*y(k,157)*y(k,95)
         loss(k,146) = (rxt(k,264)* y(k,157) + rxt(k,120) + rxt(k,121) + rxt(k,886) &
                  + rxt(k,889) + rxt(k,894) + het_rates(k,21))* y(k,21)
         prod(k,146) =rxt(k,263)*y(k,148)*y(k,20)
         loss(k,4) = ( + het_rates(k,22))* y(k,22)
         prod(k,4) = 0._r8
         loss(k,89) = (rxt(k,547)* y(k,295) + het_rates(k,23))* y(k,23)
         prod(k,89) =rxt(k,27)*y(k,24) +rxt(k,550)*y(k,248)*y(k,147)
         loss(k,113) = (rxt(k,549)* y(k,295) + rxt(k,27) + het_rates(k,24))* y(k,24)
         prod(k,113) =rxt(k,548)*y(k,258)*y(k,248)
         loss(k,105) = (rxt(k,320)* y(k,57) +rxt(k,321)* y(k,295) + het_rates(k,25)) &
                 * y(k,25)
         prod(k,105) = 0._r8
         loss(k,149) = (rxt(k,322)* y(k,57) +rxt(k,323)* y(k,158) +rxt(k,350) &
                 * y(k,295) + het_rates(k,26))* y(k,26)
         prod(k,149) = 0._r8
         loss(k,99) = (rxt(k,328)* y(k,295) + het_rates(k,27))* y(k,27)
         prod(k,99) = (.400_r8*rxt(k,324)*y(k,249) +.200_r8*rxt(k,325)*y(k,253)) &
                 *y(k,249)
         loss(k,114) = (rxt(k,329)* y(k,295) + rxt(k,28) + het_rates(k,28))* y(k,28)
         prod(k,114) =rxt(k,326)*y(k,258)*y(k,249)
         loss(k,106) = (rxt(k,330)* y(k,57) +rxt(k,331)* y(k,295) + het_rates(k,29)) &
                 * y(k,29)
         prod(k,106) = 0._r8
         loss(k,228) = (rxt(k,353)* y(k,149) +rxt(k,354)* y(k,158) +rxt(k,372) &
                 * y(k,295) + het_rates(k,30))* y(k,30)
         prod(k,228) =.700_r8*rxt(k,79)*y(k,132)
         loss(k,123) = (rxt(k,358)* y(k,295) + rxt(k,29) + het_rates(k,31))* y(k,31)
         prod(k,123) =rxt(k,356)*y(k,258)*y(k,250)
         loss(k,63) = (rxt(k,359)* y(k,295) + het_rates(k,32))* y(k,32)
         prod(k,63) = 0._r8
         loss(k,100) = (rxt(k,553)* y(k,295) + rxt(k,30) + het_rates(k,33))* y(k,33)
         prod(k,100) =rxt(k,551)*y(k,258)*y(k,251)
         loss(k,60) = (rxt(k,243)* y(k,294) + rxt(k,122) + het_rates(k,34))* y(k,34)
         prod(k,60) = 0._r8
         loss(k,71) = (rxt(k,244)* y(k,294) + rxt(k,123) + het_rates(k,35))* y(k,35)
         prod(k,71) = 0._r8
         loss(k,72) = (rxt(k,270)* y(k,294) + rxt(k,124) + het_rates(k,36))* y(k,36)
         prod(k,72) = 0._r8
         loss(k,64) = (rxt(k,245)* y(k,294) + rxt(k,125) + het_rates(k,37))* y(k,37)
         prod(k,64) = 0._r8
         loss(k,73) = (rxt(k,246)* y(k,294) + rxt(k,126) + het_rates(k,38))* y(k,38)
         prod(k,73) = 0._r8
         loss(k,65) = (rxt(k,247)* y(k,294) + rxt(k,127) + het_rates(k,39))* y(k,39)
         prod(k,65) = 0._r8
         loss(k,74) = (rxt(k,248)* y(k,294) + rxt(k,128) + het_rates(k,40))* y(k,40)
         prod(k,74) = 0._r8
         loss(k,66) = (rxt(k,249)* y(k,294) + rxt(k,129) + het_rates(k,41))* y(k,41)
         prod(k,66) = 0._r8
         loss(k,137) = (rxt(k,281)* y(k,57) +rxt(k,293)* y(k,294) +rxt(k,282) &
                 * y(k,295) + rxt(k,130) + het_rates(k,42))* y(k,42)
         prod(k,137) = 0._r8
         loss(k,306) = (rxt(k,254)* y(k,18) +rxt(k,218)* y(k,57) +rxt(k,299)* y(k,149) &
                  +rxt(k,300)* y(k,157) +rxt(k,298)* y(k,258) +rxt(k,301)* y(k,295) &
                  + rxt(k,31) + rxt(k,32) + het_rates(k,43))* y(k,43)
         prod(k,306) = (rxt(k,225)*y(k,60) +2.000_r8*rxt(k,302)*y(k,253) + &
                 rxt(k,303)*y(k,253) +rxt(k,305)*y(k,147) + &
                 .700_r8*rxt(k,325)*y(k,249) +rxt(k,336)*y(k,252) + &
                 rxt(k,355)*y(k,250) +.800_r8*rxt(k,368)*y(k,298) + &
                 1.100_r8*rxt(k,382)*y(k,284) +2.000_r8*rxt(k,389)*y(k,286) + &
                 .870_r8*rxt(k,401)*y(k,289) +1.750_r8*rxt(k,425)*y(k,261) + &
                 1.250_r8*rxt(k,431)*y(k,262) +.750_r8*rxt(k,445)*y(k,267) + &
                 .750_r8*rxt(k,449)*y(k,268) +.710_r8*rxt(k,475)*y(k,274) + &
                 .750_r8*rxt(k,492)*y(k,278) +.750_r8*rxt(k,496)*y(k,279) + &
                 .950_r8*rxt(k,587)*y(k,237) +.830_r8*rxt(k,595)*y(k,238) + &
                 .950_r8*rxt(k,607)*y(k,240) +.750_r8*rxt(k,615)*y(k,241) + &
                 .990_r8*rxt(k,627)*y(k,245) +1.400_r8*rxt(k,635)*y(k,246) + &
                 .910_r8*rxt(k,646)*y(k,281) +1.030_r8*rxt(k,655)*y(k,282) + &
                 .980_r8*rxt(k,666)*y(k,290) +.750_r8*rxt(k,675)*y(k,291) + &
                 .750_r8*rxt(k,694)*y(k,301) +rxt(k,702)*y(k,302) + &
                 rxt(k,710)*y(k,303) +rxt(k,720)*y(k,304) +rxt(k,729)*y(k,305) + &
                 3.000_r8*rxt(k,739)*y(k,306) +rxt(k,750)*y(k,307))*y(k,253) &
                  + (.500_r8*rxt(k,342)*y(k,257) +rxt(k,366)*y(k,297) + &
                 rxt(k,370)*y(k,298) +.500_r8*rxt(k,377)*y(k,255) + &
                 rxt(k,392)*y(k,286) +.100_r8*rxt(k,410)*y(k,236) + &
                 rxt(k,505)*y(k,261) +rxt(k,507)*y(k,262) + &
                 .060_r8*rxt(k,513)*y(k,269) +.270_r8*rxt(k,515)*y(k,270) + &
                 rxt(k,517)*y(k,271) +.130_r8*rxt(k,519)*y(k,272) + &
                 .330_r8*rxt(k,521)*y(k,273) +.460_r8*rxt(k,523)*y(k,274) + &
                 .530_r8*rxt(k,525)*y(k,275) +.040_r8*rxt(k,527)*y(k,276) + &
                 .140_r8*rxt(k,535)*y(k,284) +.240_r8*rxt(k,537)*y(k,289) + &
                 .210_r8*rxt(k,597)*y(k,238) +.020_r8*rxt(k,629)*y(k,245) + &
                 .490_r8*rxt(k,637)*y(k,246) +.430_r8*rxt(k,657)*y(k,282) + &
                 .040_r8*rxt(k,669)*y(k,290) +.300_r8*rxt(k,677)*y(k,291) + &
                 .310_r8*rxt(k,688)*y(k,299) +1.820_r8*rxt(k,741)*y(k,306) + &
                 .310_r8*rxt(k,761)*y(k,308))*y(k,147) &
                  + (.150_r8*rxt(k,369)*y(k,298) +.080_r8*rxt(k,383)*y(k,284) + &
                 .490_r8*rxt(k,390)*y(k,286) +.050_r8*rxt(k,402)*y(k,289) + &
                 .060_r8*rxt(k,426)*y(k,261) +.060_r8*rxt(k,432)*y(k,262) + &
                 .030_r8*rxt(k,457)*y(k,269) +.060_r8*rxt(k,461)*y(k,270) + &
                 .600_r8*rxt(k,464)*y(k,271) +.060_r8*rxt(k,467)*y(k,272) + &
                 .100_r8*rxt(k,471)*y(k,273) +.240_r8*rxt(k,476)*y(k,274) + &
                 .170_r8*rxt(k,479)*y(k,275) +.030_r8*rxt(k,482)*y(k,276) + &
                 .080_r8*rxt(k,596)*y(k,238) +.020_r8*rxt(k,628)*y(k,245) + &
                 .030_r8*rxt(k,636)*y(k,246) +.060_r8*rxt(k,656)*y(k,282) + &
                 .020_r8*rxt(k,667)*y(k,290) +.040_r8*rxt(k,676)*y(k,291) + &
                 .080_r8*rxt(k,687)*y(k,299) +1.060_r8*rxt(k,740)*y(k,306) + &
                 .040_r8*rxt(k,760)*y(k,308))*y(k,258) + (rxt(k,306)*y(k,53) + &
                 .300_r8*rxt(k,307)*y(k,54) +.500_r8*rxt(k,311)*y(k,92) + &
                 .500_r8*rxt(k,340)*y(k,52) +.800_r8*rxt(k,345)*y(k,76) + &
                 .110_r8*rxt(k,347)*y(k,89) +rxt(k,348)*y(k,150) + &
                 rxt(k,349)*y(k,163) +.300_r8*rxt(k,363)*y(k,104) + &
                 .400_r8*rxt(k,408)*y(k,1) +.500_r8*rxt(k,419)*y(k,105) + &
                 .400_r8*rxt(k,422)*y(k,107) +.590_r8*rxt(k,423)*y(k,108) + &
                 2.000_r8*rxt(k,718)*y(k,204) +rxt(k,737)*y(k,206))*y(k,295) &
                  + (.140_r8*rxt(k,381)*y(k,284) +rxt(k,388)*y(k,286) + &
                 .250_r8*rxt(k,400)*y(k,289) +rxt(k,424)*y(k,261) + &
                 rxt(k,430)*y(k,262) +.460_r8*rxt(k,474)*y(k,274) + &
                 .270_r8*rxt(k,594)*y(k,238) +.020_r8*rxt(k,626)*y(k,245) + &
                 .650_r8*rxt(k,634)*y(k,246) +.560_r8*rxt(k,654)*y(k,282) + &
                 .040_r8*rxt(k,665)*y(k,290) +.420_r8*rxt(k,674)*y(k,291) + &
                 2.000_r8*rxt(k,738)*y(k,306))*y(k,252) &
                  + (.500_r8*rxt(k,374)*y(k,16) +rxt(k,393)*y(k,286) + &
                 .460_r8*rxt(k,478)*y(k,274) +.270_r8*rxt(k,598)*y(k,238) + &
                 .020_r8*rxt(k,630)*y(k,245) +.650_r8*rxt(k,638)*y(k,246) + &
                 .560_r8*rxt(k,658)*y(k,282) +.040_r8*rxt(k,670)*y(k,290) + &
                 .420_r8*rxt(k,678)*y(k,291) +2.000_r8*rxt(k,742)*y(k,306) + &
                 .440_r8*rxt(k,759)*y(k,212) +.500_r8*rxt(k,764)*y(k,213))*y(k,149) &
                  + (rxt(k,323)*y(k,26) +.500_r8*rxt(k,354)*y(k,30) + &
                 .120_r8*rxt(k,385)*y(k,126) +.600_r8*rxt(k,403)*y(k,132) + &
                 1.010_r8*rxt(k,486)*y(k,109) +.270_r8*rxt(k,602)*y(k,4) + &
                 .080_r8*rxt(k,622)*y(k,7) +.810_r8*rxt(k,642)*y(k,17) + &
                 .330_r8*rxt(k,662)*y(k,125) +.390_r8*rxt(k,682)*y(k,135) + &
                 .620_r8*rxt(k,762)*y(k,212) +.340_r8*rxt(k,767)*y(k,213))*y(k,158) &
                  + (.270_r8*rxt(k,599)*y(k,238) +.020_r8*rxt(k,631)*y(k,245) + &
                 .650_r8*rxt(k,639)*y(k,246) +.560_r8*rxt(k,659)*y(k,282) + &
                 .040_r8*rxt(k,671)*y(k,290) +.420_r8*rxt(k,679)*y(k,291) + &
                 2.000_r8*rxt(k,743)*y(k,306))*y(k,302) &
                  + (.270_r8*rxt(k,600)*y(k,238) +.020_r8*rxt(k,632)*y(k,245) + &
                 .650_r8*rxt(k,640)*y(k,246) +.560_r8*rxt(k,660)*y(k,282) + &
                 .040_r8*rxt(k,672)*y(k,290) +.420_r8*rxt(k,680)*y(k,291) + &
                 2.000_r8*rxt(k,744)*y(k,306))*y(k,304) &
                  + (.270_r8*rxt(k,601)*y(k,238) +.020_r8*rxt(k,633)*y(k,245) + &
                 .650_r8*rxt(k,641)*y(k,246) +.560_r8*rxt(k,661)*y(k,282) + &
                 .040_r8*rxt(k,673)*y(k,290) +.420_r8*rxt(k,681)*y(k,291) + &
                 2.000_r8*rxt(k,745)*y(k,306))*y(k,307) + (.180_r8*rxt(k,39) + &
                 rxt(k,316)*y(k,294) +rxt(k,317)*y(k,294))*y(k,55) + (rxt(k,55) + &
                 rxt(k,56))*y(k,104) +.100_r8*rxt(k,19)*y(k,1) +.100_r8*rxt(k,20) &
                 *y(k,2) +rxt(k,37)*y(k,54) +.500_r8*rxt(k,41)*y(k,68) +rxt(k,43) &
                 *y(k,76) +rxt(k,45)*y(k,89) +rxt(k,46)*y(k,92) +.330_r8*rxt(k,47) &
                 *y(k,97) +rxt(k,52)*y(k,102) +rxt(k,65)*y(k,116) +rxt(k,66)*y(k,117) &
                  +rxt(k,68)*y(k,119) +rxt(k,69)*y(k,120) +rxt(k,71)*y(k,123) &
                  +rxt(k,72)*y(k,126) +.250_r8*rxt(k,74)*y(k,127) +.140_r8*rxt(k,75) &
                 *y(k,128) +.250_r8*rxt(k,80)*y(k,133) +.440_r8*rxt(k,81)*y(k,134) &
                  +rxt(k,83)*y(k,150) +rxt(k,84)*y(k,151) +rxt(k,88)*y(k,170) &
                  +rxt(k,89)*y(k,171) +.040_r8*rxt(k,625)*y(k,245)*y(k,245) &
                  +2.000_r8*rxt(k,343)*y(k,256) +rxt(k,313)*y(k,259) +rxt(k,427) &
                 *y(k,261) +rxt(k,433)*y(k,262) +.160_r8*rxt(k,477)*y(k,274)*y(k,274) &
                  +2.000_r8*rxt(k,391)*y(k,286)*y(k,286) +.060_r8*rxt(k,668)*y(k,290) &
                 *y(k,290)
         loss(k,156) = (rxt(k,283)* y(k,57) +rxt(k,294)* y(k,294) +rxt(k,284) &
                 * y(k,295) + rxt(k,131) + het_rates(k,44))* y(k,44)
         prod(k,156) = 0._r8
         loss(k,67) = (rxt(k,285)* y(k,295) + rxt(k,132) + het_rates(k,45))* y(k,45)
         prod(k,67) = 0._r8
         loss(k,230) = (rxt(k,332)* y(k,149) +rxt(k,333)* y(k,295) + rxt(k,33) &
                  + het_rates(k,46))* y(k,46)
         prod(k,230) = (rxt(k,327)*y(k,249) +.270_r8*rxt(k,357)*y(k,250) + &
                 rxt(k,366)*y(k,297) +rxt(k,377)*y(k,255) +rxt(k,395)*y(k,288) + &
                 .400_r8*rxt(k,410)*y(k,236))*y(k,147) + (rxt(k,328)*y(k,27) + &
                 .500_r8*rxt(k,329)*y(k,28) +.800_r8*rxt(k,408)*y(k,1))*y(k,295) &
                  + (.500_r8*rxt(k,354)*y(k,30) +.100_r8*rxt(k,403)*y(k,132))*y(k,158) &
                  + (1.600_r8*rxt(k,324)*y(k,249) +.800_r8*rxt(k,325)*y(k,253)) &
                 *y(k,249) +.400_r8*rxt(k,19)*y(k,1) +.400_r8*rxt(k,20)*y(k,2) &
                  +rxt(k,374)*y(k,149)*y(k,16) +rxt(k,28)*y(k,28) +.330_r8*rxt(k,47) &
                 *y(k,97) +rxt(k,77)*y(k,130) +rxt(k,88)*y(k,170) &
                  +.200_r8*rxt(k,394)*y(k,288)*y(k,258)
         loss(k,120) = (rxt(k,286)* y(k,57) +rxt(k,287)* y(k,295) + rxt(k,133) &
                  + het_rates(k,47))* y(k,47)
         prod(k,120) = 0._r8
         loss(k,61) = (rxt(k,334)* y(k,295) + het_rates(k,48))* y(k,48)
         prod(k,61) = 0._r8
         loss(k,280) = (rxt(k,373)* y(k,295) + rxt(k,34) + het_rates(k,49))* y(k,49)
         prod(k,280) = (.910_r8*rxt(k,665)*y(k,252) +.740_r8*rxt(k,666)*y(k,253) + &
                 .460_r8*rxt(k,667)*y(k,258) +1.480_r8*rxt(k,668)*y(k,290) + &
                 .850_r8*rxt(k,669)*y(k,147) +.910_r8*rxt(k,670)*y(k,149) + &
                 .910_r8*rxt(k,671)*y(k,302) +.910_r8*rxt(k,672)*y(k,304) + &
                 .910_r8*rxt(k,673)*y(k,307))*y(k,290) &
                  + (.120_r8*rxt(k,594)*y(k,252) +.060_r8*rxt(k,595)*y(k,253) + &
                 .060_r8*rxt(k,596)*y(k,258) +.090_r8*rxt(k,597)*y(k,147) + &
                 .120_r8*rxt(k,598)*y(k,149) +.120_r8*rxt(k,599)*y(k,302) + &
                 .120_r8*rxt(k,600)*y(k,304) +.120_r8*rxt(k,601)*y(k,307))*y(k,238) &
                  + (rxt(k,728)*y(k,252) +rxt(k,729)*y(k,253) + &
                 .150_r8*rxt(k,730)*y(k,258) +.700_r8*rxt(k,731)*y(k,147) + &
                 rxt(k,732)*y(k,149) +rxt(k,733)*y(k,302) +rxt(k,734)*y(k,304) + &
                 rxt(k,735)*y(k,307))*y(k,305) + (.110_r8*rxt(k,634)*y(k,252) + &
                 .080_r8*rxt(k,635)*y(k,253) +.080_r8*rxt(k,637)*y(k,147) + &
                 .110_r8*rxt(k,638)*y(k,149) +.110_r8*rxt(k,639)*y(k,302) + &
                 .110_r8*rxt(k,640)*y(k,304) +.110_r8*rxt(k,641)*y(k,307))*y(k,246) &
                  + (.460_r8*rxt(k,674)*y(k,252) +.050_r8*rxt(k,676)*y(k,258) + &
                 .330_r8*rxt(k,677)*y(k,147) +.460_r8*rxt(k,678)*y(k,149) + &
                 .460_r8*rxt(k,679)*y(k,302) +.460_r8*rxt(k,680)*y(k,304) + &
                 .460_r8*rxt(k,681)*y(k,307))*y(k,291) &
                  + (.820_r8*rxt(k,357)*y(k,250) +.500_r8*rxt(k,377)*y(k,255) + &
                 .250_r8*rxt(k,410)*y(k,236))*y(k,147) + (.250_r8*rxt(k,19) + &
                 .800_r8*rxt(k,408)*y(k,295))*y(k,1) + (.820_r8*rxt(k,355)*y(k,250) + &
                 .100_r8*rxt(k,382)*y(k,284))*y(k,253) +.250_r8*rxt(k,20)*y(k,2) &
                  +.500_r8*rxt(k,374)*y(k,149)*y(k,16) +.820_r8*rxt(k,29)*y(k,31) &
                  +.170_r8*rxt(k,47)*y(k,97) +.250_r8*rxt(k,682)*y(k,158)*y(k,135) &
                  +rxt(k,718)*y(k,295)*y(k,204)
         loss(k,266) = (rxt(k,360)* y(k,149) +rxt(k,361)* y(k,295) + rxt(k,35) &
                  + het_rates(k,50))* y(k,50)
         prod(k,266) = (rxt(k,362)*y(k,102) +.700_r8*rxt(k,363)*y(k,104) + &
                 rxt(k,364)*y(k,151) +.440_r8*rxt(k,405)*y(k,134) + &
                 .380_r8*rxt(k,414)*y(k,98) +.030_r8*rxt(k,415)*y(k,99) + &
                 .460_r8*rxt(k,418)*y(k,103) +.500_r8*rxt(k,419)*y(k,105) + &
                 .400_r8*rxt(k,422)*y(k,107) +.720_r8*rxt(k,456)*y(k,114))*y(k,295) &
                  + (.710_r8*rxt(k,503)*y(k,260) +.140_r8*rxt(k,535)*y(k,284) + &
                 .240_r8*rxt(k,537)*y(k,289) +.120_r8*rxt(k,539)*y(k,293) + &
                 .170_r8*rxt(k,556)*y(k,254) +.170_r8*rxt(k,562)*y(k,287) + &
                 .400_r8*rxt(k,572)*y(k,314) +.540_r8*rxt(k,578)*y(k,316) + &
                 .510_r8*rxt(k,581)*y(k,318))*y(k,147) &
                  + (.880_r8*rxt(k,385)*y(k,126) +.500_r8*rxt(k,403)*y(k,132) + &
                 .170_r8*rxt(k,459)*y(k,115) +.170_r8*rxt(k,469)*y(k,118) + &
                 .170_r8*rxt(k,484)*y(k,121) +.340_r8*rxt(k,501)*y(k,139))*y(k,158) &
                  + (.080_r8*rxt(k,383)*y(k,284) +.050_r8*rxt(k,402)*y(k,289) + &
                 .460_r8*rxt(k,421)*y(k,260) +.100_r8*rxt(k,499)*y(k,293) + &
                 .070_r8*rxt(k,555)*y(k,254) +.070_r8*rxt(k,561)*y(k,287))*y(k,258) &
                  + (.140_r8*rxt(k,381)*y(k,284) +.250_r8*rxt(k,400)*y(k,289)) &
                 *y(k,252) + (.500_r8*rxt(k,368)*y(k,298) + &
                 .120_r8*rxt(k,401)*y(k,289))*y(k,253) +rxt(k,26)*y(k,14) &
                  +.500_r8*rxt(k,41)*y(k,68) +.680_r8*rxt(k,48)*y(k,98) &
                  +.670_r8*rxt(k,49)*y(k,99) +rxt(k,54)*y(k,103) +.500_r8*rxt(k,60) &
                 *y(k,111) +.500_r8*rxt(k,61)*y(k,112) +.720_r8*rxt(k,63)*y(k,114) &
                  +.250_r8*rxt(k,74)*y(k,127) +.140_r8*rxt(k,75)*y(k,128) &
                  +.250_r8*rxt(k,80)*y(k,133) +.440_r8*rxt(k,81)*y(k,134) &
                  +.400_r8*rxt(k,115)*y(k,227) +.540_r8*rxt(k,116)*y(k,230) &
                  +.510_r8*rxt(k,117)*y(k,232)
         loss(k,173) = (rxt(k,339)* y(k,295) + het_rates(k,51))* y(k,51)
         prod(k,173) = (.100_r8*rxt(k,336)*y(k,253) +.150_r8*rxt(k,337)*y(k,258)) &
                 *y(k,252) +.120_r8*rxt(k,354)*y(k,158)*y(k,30) &
                  +.150_r8*rxt(k,390)*y(k,286)*y(k,258)
         loss(k,163) = (rxt(k,340)* y(k,295) + rxt(k,36) + het_rates(k,52))* y(k,52)
         prod(k,163) = (.360_r8*rxt(k,337)*y(k,252) +.360_r8*rxt(k,390)*y(k,286)) &
                 *y(k,258)
         loss(k,238) = (rxt(k,306)* y(k,295) + het_rates(k,53))* y(k,53)
         prod(k,238) = (rxt(k,303)*y(k,253) +.300_r8*rxt(k,325)*y(k,249) + &
                 .500_r8*rxt(k,368)*y(k,298) +.250_r8*rxt(k,401)*y(k,289) + &
                 .250_r8*rxt(k,431)*y(k,262) +.250_r8*rxt(k,445)*y(k,267) + &
                 .250_r8*rxt(k,449)*y(k,268) +.360_r8*rxt(k,475)*y(k,274) + &
                 .250_r8*rxt(k,492)*y(k,278) +.250_r8*rxt(k,496)*y(k,279) + &
                 .050_r8*rxt(k,587)*y(k,237) +.170_r8*rxt(k,595)*y(k,238) + &
                 .050_r8*rxt(k,607)*y(k,240) +.250_r8*rxt(k,615)*y(k,241) + &
                 .030_r8*rxt(k,627)*y(k,245) +.090_r8*rxt(k,646)*y(k,281) + &
                 .250_r8*rxt(k,655)*y(k,282) +.050_r8*rxt(k,666)*y(k,290) + &
                 .250_r8*rxt(k,675)*y(k,291) +.250_r8*rxt(k,694)*y(k,301))*y(k,253)
         loss(k,129) = (rxt(k,307)* y(k,295) + rxt(k,37) + het_rates(k,54))* y(k,54)
         prod(k,129) =rxt(k,304)*y(k,258)*y(k,253)
         loss(k,279) = (rxt(k,219)* y(k,57) +rxt(k,275)* y(k,75) + (rxt(k,315) + &
                 rxt(k,316) +rxt(k,317))* y(k,294) +rxt(k,308)* y(k,295) + rxt(k,38) &
                  + rxt(k,39) + het_rates(k,55))* y(k,55)
         prod(k,279) =.100_r8*rxt(k,354)*y(k,158)*y(k,30)
         loss(k,131) = (rxt(k,288)* y(k,57) +rxt(k,271)* y(k,294) +rxt(k,289) &
                 * y(k,295) + rxt(k,134) + het_rates(k,56))* y(k,56)
         prod(k,131) = 0._r8
         loss(k,315) = (rxt(k,330)* y(k,29) +rxt(k,281)* y(k,42) +rxt(k,218)* y(k,43) &
                  +rxt(k,283)* y(k,44) +rxt(k,286)* y(k,47) +rxt(k,219)* y(k,55) &
                  +rxt(k,288)* y(k,56) +rxt(k,231)* y(k,61) +rxt(k,220)* y(k,79) &
                  +rxt(k,221)* y(k,81) +rxt(k,240)* y(k,96) +rxt(k,224)* y(k,158) &
                  + (rxt(k,222) +rxt(k,223))* y(k,258) + het_rates(k,57))* y(k,57)
         prod(k,315) = (4.000_r8*rxt(k,243)*y(k,34) +rxt(k,244)*y(k,35) + &
                 2.000_r8*rxt(k,245)*y(k,37) +2.000_r8*rxt(k,246)*y(k,38) + &
                 2.000_r8*rxt(k,247)*y(k,39) +rxt(k,248)*y(k,40) + &
                 2.000_r8*rxt(k,249)*y(k,41) +rxt(k,250)*y(k,87) +rxt(k,280)*y(k,66) + &
                 rxt(k,295)*y(k,84) +rxt(k,296)*y(k,85) +rxt(k,297)*y(k,86))*y(k,294) &
                  + (rxt(k,137) +rxt(k,225)*y(k,253) +2.000_r8*rxt(k,226)*y(k,60) + &
                 rxt(k,228)*y(k,60) +rxt(k,230)*y(k,147) +rxt(k,235)*y(k,157) + &
                 rxt(k,236)*y(k,295) +rxt(k,259)*y(k,20) +rxt(k,801)*y(k,174))*y(k,60) &
                  + (rxt(k,239)*y(k,87) +3.000_r8*rxt(k,285)*y(k,45) + &
                 rxt(k,287)*y(k,47) +rxt(k,290)*y(k,84) +rxt(k,291)*y(k,85) + &
                 rxt(k,292)*y(k,86))*y(k,295) + (rxt(k,147) +rxt(k,238)*y(k,157)) &
                 *y(k,87) +rxt(k,118)*y(k,19) +4.000_r8*rxt(k,122)*y(k,34) +rxt(k,123) &
                 *y(k,35) +2.000_r8*rxt(k,125)*y(k,37) +2.000_r8*rxt(k,126)*y(k,38) &
                  +2.000_r8*rxt(k,127)*y(k,39) +rxt(k,128)*y(k,40) &
                  +2.000_r8*rxt(k,129)*y(k,41) +3.000_r8*rxt(k,132)*y(k,45) &
                  +rxt(k,133)*y(k,47) +2.000_r8*rxt(k,135)*y(k,58) &
                  +2.000_r8*rxt(k,136)*y(k,59) +rxt(k,138)*y(k,61) +rxt(k,141)*y(k,66) &
                  +rxt(k,144)*y(k,84) +rxt(k,145)*y(k,85) +rxt(k,146)*y(k,86) &
                  +rxt(k,150)*y(k,96)
         loss(k,75) = ( + rxt(k,135) + het_rates(k,58))* y(k,58)
         prod(k,75) = (rxt(k,885)*y(k,96) +rxt(k,890)*y(k,61) +rxt(k,891)*y(k,96) + &
                 rxt(k,895)*y(k,61) +rxt(k,896)*y(k,96) +rxt(k,900)*y(k,61))*y(k,87) &
                  +rxt(k,231)*y(k,61)*y(k,57) +rxt(k,227)*y(k,60)*y(k,60)
         loss(k,58) = ( + rxt(k,136) + rxt(k,253) + het_rates(k,59))* y(k,59)
         prod(k,58) =rxt(k,252)*y(k,60)*y(k,60)
         loss(k,310) = ((rxt(k,258) +rxt(k,259) +rxt(k,260))* y(k,20) &
                  + 2._r8*(rxt(k,226) +rxt(k,227) +rxt(k,228) +rxt(k,252))* y(k,60) &
                  +rxt(k,230)* y(k,147) +rxt(k,232)* y(k,148) +rxt(k,235)* y(k,157) &
                  +rxt(k,801)* y(k,174) +rxt(k,225)* y(k,253) +rxt(k,229)* y(k,258) &
                  + (rxt(k,236) +rxt(k,237))* y(k,295) + rxt(k,137) + het_rates(k,60)) &
                 * y(k,60)
         prod(k,310) = (rxt(k,223)*y(k,258) +rxt(k,224)*y(k,158) +rxt(k,240)*y(k,96)) &
                 *y(k,57) + (rxt(k,139) +rxt(k,233)*y(k,157))*y(k,61) &
                  + (rxt(k,241)*y(k,157) +rxt(k,242)*y(k,295))*y(k,96) + (rxt(k,151) + &
                 rxt(k,806)*y(k,174))*y(k,160) +2.000_r8*rxt(k,253)*y(k,59) &
                  +rxt(k,251)*y(k,294)*y(k,87)
         loss(k,223) = (rxt(k,231)* y(k,57) + (rxt(k,890) +rxt(k,895) +rxt(k,900)) &
                 * y(k,87) +rxt(k,233)* y(k,157) +rxt(k,234)* y(k,295) + rxt(k,138) &
                  + rxt(k,139) + rxt(k,888) + rxt(k,893) + rxt(k,899) &
                  + het_rates(k,61))* y(k,61)
         prod(k,223) =rxt(k,232)*y(k,148)*y(k,60)
         loss(k,5) = ( + het_rates(k,62))* y(k,62)
         prod(k,5) = 0._r8
         loss(k,271) = (rxt(k,319)* y(k,295) + het_rates(k,63))* y(k,63)
         prod(k,271) = (rxt(k,301)*y(k,43) +.350_r8*rxt(k,321)*y(k,25) + &
                 rxt(k,346)*y(k,77) +.110_r8*rxt(k,347)*y(k,89) +rxt(k,361)*y(k,50) + &
                 rxt(k,376)*y(k,68) +rxt(k,380)*y(k,127) +rxt(k,387)*y(k,128) + &
                 .250_r8*rxt(k,398)*y(k,131) +.500_r8*rxt(k,399)*y(k,133) + &
                 1.560_r8*rxt(k,405)*y(k,134) +1.060_r8*rxt(k,414)*y(k,98) + &
                 .760_r8*rxt(k,415)*y(k,99) +.420_r8*rxt(k,416)*y(k,100) + &
                 .230_r8*rxt(k,417)*y(k,101) +rxt(k,418)*y(k,103) + &
                 1.500_r8*rxt(k,419)*y(k,105) +.350_r8*rxt(k,423)*y(k,108) + &
                 rxt(k,452)*y(k,111) +rxt(k,454)*y(k,112) + &
                 2.000_r8*rxt(k,456)*y(k,114) +.060_r8*rxt(k,460)*y(k,115) + &
                 .040_r8*rxt(k,470)*y(k,118) +.630_r8*rxt(k,502)*y(k,139) + &
                 2.000_r8*rxt(k,718)*y(k,204) +rxt(k,737)*y(k,206) + &
                 rxt(k,757)*y(k,210) +rxt(k,796)*y(k,161))*y(k,295) &
                  + (.650_r8*rxt(k,392)*y(k,286) +.400_r8*rxt(k,503)*y(k,260) + &
                 .550_r8*rxt(k,509)*y(k,267) +.550_r8*rxt(k,511)*y(k,268) + &
                 .550_r8*rxt(k,530)*y(k,278) +.550_r8*rxt(k,533)*y(k,279) + &
                 .860_r8*rxt(k,535)*y(k,284) +.750_r8*rxt(k,539)*y(k,293) + &
                 .170_r8*rxt(k,556)*y(k,254) +.400_r8*rxt(k,559)*y(k,285) + &
                 .350_r8*rxt(k,562)*y(k,287) +.910_r8*rxt(k,741)*y(k,306))*y(k,147) &
                  + (.510_r8*rxt(k,383)*y(k,284) +.320_r8*rxt(k,390)*y(k,286) + &
                 .260_r8*rxt(k,402)*y(k,289) +.260_r8*rxt(k,421)*y(k,260) + &
                 .600_r8*rxt(k,499)*y(k,293) +.070_r8*rxt(k,555)*y(k,254) + &
                 .160_r8*rxt(k,558)*y(k,285) +.140_r8*rxt(k,561)*y(k,287) + &
                 .530_r8*rxt(k,740)*y(k,306))*y(k,258) &
                  + (.900_r8*rxt(k,382)*y(k,284) +.650_r8*rxt(k,389)*y(k,286) + &
                 rxt(k,401)*y(k,289) +.280_r8*rxt(k,445)*y(k,267) + &
                 .280_r8*rxt(k,449)*y(k,268) +.280_r8*rxt(k,492)*y(k,278) + &
                 .280_r8*rxt(k,496)*y(k,279) +rxt(k,739)*y(k,306))*y(k,253) &
                  + (.630_r8*rxt(k,323)*y(k,26) +.560_r8*rxt(k,354)*y(k,30) + &
                 .650_r8*rxt(k,385)*y(k,126) +.560_r8*rxt(k,403)*y(k,132) + &
                 .350_r8*rxt(k,486)*y(k,109) +.300_r8*rxt(k,501)*y(k,139) + &
                 .170_r8*rxt(k,602)*y(k,4))*y(k,158) + (.860_r8*rxt(k,381)*y(k,284) + &
                 .650_r8*rxt(k,388)*y(k,286) +.550_r8*rxt(k,444)*y(k,267) + &
                 .550_r8*rxt(k,448)*y(k,268) +.550_r8*rxt(k,491)*y(k,278) + &
                 .550_r8*rxt(k,495)*y(k,279) +rxt(k,738)*y(k,306))*y(k,252) &
                  + (rxt(k,31) +rxt(k,32) +rxt(k,218)*y(k,57) +rxt(k,254)*y(k,18) + &
                 rxt(k,299)*y(k,149) +rxt(k,300)*y(k,157))*y(k,43) &
                  + (rxt(k,742)*y(k,149) +rxt(k,743)*y(k,302) +rxt(k,744)*y(k,304) + &
                 rxt(k,745)*y(k,307))*y(k,306) + (rxt(k,35) +rxt(k,360)*y(k,149)) &
                 *y(k,50) + (1.500_r8*rxt(k,53) +rxt(k,54))*y(k,103) + (rxt(k,154) + &
                 rxt(k,795)*y(k,157))*y(k,161) + (1.300_r8*rxt(k,391)*y(k,286) + &
                 .650_r8*rxt(k,393)*y(k,149))*y(k,286) +1.500_r8*rxt(k,22)*y(k,10) &
                  +.600_r8*rxt(k,25)*y(k,13) +rxt(k,26)*y(k,14) +rxt(k,33)*y(k,46) &
                  +rxt(k,286)*y(k,57)*y(k,47) +.380_r8*rxt(k,39)*y(k,55) +rxt(k,40) &
                 *y(k,64) +.500_r8*rxt(k,41)*y(k,68) +rxt(k,43)*y(k,76) &
                  +2.000_r8*rxt(k,44)*y(k,77) +rxt(k,45)*y(k,89) +.330_r8*rxt(k,47) &
                 *y(k,97) +1.320_r8*rxt(k,48)*y(k,98) +1.740_r8*rxt(k,49)*y(k,99) &
                  +rxt(k,50)*y(k,100) +rxt(k,51)*y(k,101) +.550_r8*rxt(k,64)*y(k,115) &
                  +.550_r8*rxt(k,67)*y(k,118) +1.650_r8*rxt(k,72)*y(k,126) &
                  +.750_r8*rxt(k,74)*y(k,127) +.860_r8*rxt(k,75)*y(k,128) &
                  +.700_r8*rxt(k,79)*y(k,132) +rxt(k,83)*y(k,150) +1.500_r8*rxt(k,90) &
                 *y(k,199) +rxt(k,93)*y(k,202) +rxt(k,94)*y(k,203) +rxt(k,96)*y(k,205) &
                  +.600_r8*rxt(k,529)*y(k,278) +.600_r8*rxt(k,532)*y(k,279) &
                  +rxt(k,384)*y(k,284) +rxt(k,500)*y(k,293)
         loss(k,247) = ( + rxt(k,40) + het_rates(k,64))* y(k,64)
         prod(k,247) = (2.000_r8*rxt(k,335)*y(k,252) +.900_r8*rxt(k,336)*y(k,253) + &
                 .490_r8*rxt(k,337)*y(k,258) +rxt(k,338)*y(k,147) + &
                 rxt(k,381)*y(k,284) +2.000_r8*rxt(k,388)*y(k,286) + &
                 rxt(k,400)*y(k,289) +rxt(k,424)*y(k,261) +rxt(k,430)*y(k,262) + &
                 rxt(k,444)*y(k,267) +rxt(k,448)*y(k,268) +rxt(k,474)*y(k,274) + &
                 rxt(k,491)*y(k,278) +rxt(k,495)*y(k,279) +rxt(k,586)*y(k,237) + &
                 rxt(k,594)*y(k,238) +rxt(k,606)*y(k,240) +rxt(k,614)*y(k,241) + &
                 rxt(k,626)*y(k,245) +rxt(k,634)*y(k,246) +rxt(k,645)*y(k,281) + &
                 rxt(k,654)*y(k,282) +rxt(k,665)*y(k,290) +rxt(k,674)*y(k,291) + &
                 rxt(k,693)*y(k,301) +2.000_r8*rxt(k,701)*y(k,302) + &
                 rxt(k,709)*y(k,303) +2.000_r8*rxt(k,719)*y(k,304) + &
                 rxt(k,728)*y(k,305) +rxt(k,738)*y(k,306) + &
                 2.000_r8*rxt(k,749)*y(k,307))*y(k,252) + (rxt(k,591)*y(k,237) + &
                 rxt(k,599)*y(k,238) +rxt(k,611)*y(k,240) +rxt(k,619)*y(k,241) + &
                 rxt(k,631)*y(k,245) +rxt(k,639)*y(k,246) +rxt(k,651)*y(k,281) + &
                 rxt(k,659)*y(k,282) +rxt(k,671)*y(k,290) +rxt(k,679)*y(k,291) + &
                 rxt(k,698)*y(k,301) +rxt(k,702)*y(k,253) + &
                 .490_r8*rxt(k,703)*y(k,258) +rxt(k,704)*y(k,147) + &
                 rxt(k,705)*y(k,149) +2.000_r8*rxt(k,706)*y(k,302) + &
                 2.000_r8*rxt(k,707)*y(k,307) +rxt(k,714)*y(k,303) + &
                 2.000_r8*rxt(k,724)*y(k,304) +rxt(k,733)*y(k,305) + &
                 rxt(k,743)*y(k,306))*y(k,302) + (rxt(k,592)*y(k,237) + &
                 rxt(k,600)*y(k,238) +rxt(k,612)*y(k,240) +rxt(k,620)*y(k,241) + &
                 rxt(k,632)*y(k,245) +rxt(k,640)*y(k,246) +rxt(k,652)*y(k,281) + &
                 rxt(k,660)*y(k,282) +rxt(k,672)*y(k,290) +rxt(k,680)*y(k,291) + &
                 rxt(k,699)*y(k,301) +rxt(k,715)*y(k,303) +rxt(k,720)*y(k,253) + &
                 .490_r8*rxt(k,721)*y(k,258) +rxt(k,722)*y(k,147) + &
                 rxt(k,723)*y(k,149) +2.000_r8*rxt(k,725)*y(k,304) + &
                 2.000_r8*rxt(k,726)*y(k,307) +rxt(k,734)*y(k,305) + &
                 rxt(k,744)*y(k,306))*y(k,304) + (rxt(k,593)*y(k,237) + &
                 rxt(k,601)*y(k,238) +rxt(k,613)*y(k,240) +rxt(k,621)*y(k,241) + &
                 rxt(k,633)*y(k,245) +rxt(k,641)*y(k,246) +rxt(k,653)*y(k,281) + &
                 rxt(k,661)*y(k,282) +rxt(k,673)*y(k,290) +rxt(k,681)*y(k,291) + &
                 rxt(k,700)*y(k,301) +rxt(k,716)*y(k,303) +rxt(k,735)*y(k,305) + &
                 rxt(k,745)*y(k,306) +rxt(k,750)*y(k,253) + &
                 .490_r8*rxt(k,751)*y(k,258) +rxt(k,752)*y(k,147) + &
                 rxt(k,753)*y(k,149) +2.000_r8*rxt(k,754)*y(k,307))*y(k,307) &
                  + (rxt(k,310)*y(k,90) +rxt(k,319)*y(k,63) +rxt(k,339)*y(k,51) + &
                 .500_r8*rxt(k,340)*y(k,52) +.800_r8*rxt(k,345)*y(k,76) + &
                 rxt(k,346)*y(k,77) +rxt(k,348)*y(k,150) +.540_r8*rxt(k,414)*y(k,98) + &
                 .540_r8*rxt(k,415)*y(k,99) +.360_r8*rxt(k,418)*y(k,103) + &
                 .190_r8*rxt(k,423)*y(k,108) +.450_r8*rxt(k,502)*y(k,139) + &
                 2.000_r8*rxt(k,718)*y(k,204) +3.000_r8*rxt(k,737)*y(k,206) + &
                 .290_r8*rxt(k,746)*y(k,208) +.290_r8*rxt(k,747)*y(k,209) + &
                 .290_r8*rxt(k,748)*y(k,207))*y(k,295) + (rxt(k,389)*y(k,253) + &
                 .490_r8*rxt(k,390)*y(k,258) +2.000_r8*rxt(k,391)*y(k,286) + &
                 rxt(k,392)*y(k,147) +rxt(k,393)*y(k,149))*y(k,286) &
                  + (.200_r8*rxt(k,354)*y(k,30) +.100_r8*rxt(k,403)*y(k,132) + &
                 .420_r8*rxt(k,486)*y(k,109) +.190_r8*rxt(k,642)*y(k,17))*y(k,158) &
                  +rxt(k,36)*y(k,52) +.440_r8*rxt(k,39)*y(k,55) +.170_r8*rxt(k,48) &
                 *y(k,98) +.280_r8*rxt(k,49)*y(k,99) +rxt(k,54)*y(k,103) &
                  +.400_r8*rxt(k,86)*y(k,163) +rxt(k,98)*y(k,207) +rxt(k,99)*y(k,208) &
                  +rxt(k,100)*y(k,209)
         loss(k,90) = (rxt(k,279)* y(k,294) + rxt(k,140) + het_rates(k,65))* y(k,65)
         prod(k,90) = (rxt(k,244)*y(k,35) +rxt(k,246)*y(k,38) + &
                 2.000_r8*rxt(k,247)*y(k,39) +2.000_r8*rxt(k,248)*y(k,40) + &
                 rxt(k,249)*y(k,41) +rxt(k,270)*y(k,36) +2.000_r8*rxt(k,272)*y(k,80) + &
                 rxt(k,296)*y(k,85) +rxt(k,297)*y(k,86))*y(k,294) + (rxt(k,145) + &
                 rxt(k,291)*y(k,295))*y(k,85) + (rxt(k,146) +rxt(k,292)*y(k,295)) &
                 *y(k,86) +rxt(k,123)*y(k,35) +rxt(k,124)*y(k,36) +rxt(k,126)*y(k,38) &
                  +2.000_r8*rxt(k,127)*y(k,39) +2.000_r8*rxt(k,128)*y(k,40) &
                  +rxt(k,129)*y(k,41) +2.000_r8*rxt(k,142)*y(k,80)
         loss(k,92) = (rxt(k,280)* y(k,294) + rxt(k,141) + het_rates(k,66))* y(k,66)
         prod(k,92) = (rxt(k,144) +rxt(k,290)*y(k,295) +rxt(k,295)*y(k,294))*y(k,84) &
                  + (rxt(k,125) +rxt(k,245)*y(k,294))*y(k,37) + (rxt(k,126) + &
                 rxt(k,246)*y(k,294))*y(k,38)
         loss(k,84) = (rxt(k,554)* y(k,295) + het_rates(k,67))* y(k,67)
         prod(k,84) =.180_r8*rxt(k,574)*y(k,295)*y(k,228)
         loss(k,155) = (rxt(k,376)* y(k,295) + rxt(k,41) + het_rates(k,68))* y(k,68)
         prod(k,155) = (.070_r8*rxt(k,414)*y(k,98) +.170_r8*rxt(k,415)*y(k,99)) &
                 *y(k,295) +.600_r8*rxt(k,529)*y(k,278) +.600_r8*rxt(k,532)*y(k,279)
         loss(k,104) = (rxt(k,793)* y(k,149) + (rxt(k,794) +rxt(k,808))* y(k,295) &
                  + het_rates(k,69))* y(k,69)
         prod(k,104) = 0._r8
         loss(k,6) = ( + het_rates(k,70))* y(k,70)
         prod(k,6) = 0._r8
         loss(k,7) = ( + het_rates(k,71))* y(k,71)
         prod(k,7) = 0._r8
         loss(k,8) = ( + het_rates(k,72))* y(k,72)
         prod(k,8) = 0._r8
         loss(k,9) = ( + rxt(k,901) + het_rates(k,73))* y(k,73)
         prod(k,9) = 0._r8
         loss(k,68) = ( + rxt(k,42) + het_rates(k,74))* y(k,74)
         prod(k,68) =rxt(k,341)*y(k,258)*y(k,257)
         loss(k,217) = (rxt(k,275)* y(k,55) +rxt(k,276)* y(k,79) +rxt(k,278)* y(k,93) &
                  +rxt(k,277)* y(k,319) + het_rates(k,75))* y(k,75)
         prod(k,217) = (rxt(k,248)*y(k,40) +rxt(k,270)*y(k,36) + &
                 2.000_r8*rxt(k,279)*y(k,65) +rxt(k,280)*y(k,66))*y(k,294) +rxt(k,124) &
                 *y(k,36) +rxt(k,128)*y(k,40) +2.000_r8*rxt(k,140)*y(k,65) +rxt(k,141) &
                 *y(k,66) +rxt(k,148)*y(k,91)
         loss(k,254) = (rxt(k,345)* y(k,295) + rxt(k,43) + het_rates(k,76))* y(k,76)
         prod(k,254) = (.570_r8*rxt(k,503)*y(k,260) +.940_r8*rxt(k,513)*y(k,269) + &
                 .730_r8*rxt(k,515)*y(k,270) +.340_r8*rxt(k,521)*y(k,273) + &
                 .400_r8*rxt(k,525)*y(k,275) +.760_r8*rxt(k,537)*y(k,289))*y(k,147) &
                  + (.360_r8*rxt(k,402)*y(k,289) +.370_r8*rxt(k,421)*y(k,260) + &
                 .550_r8*rxt(k,457)*y(k,269) +.460_r8*rxt(k,461)*y(k,270) + &
                 .150_r8*rxt(k,471)*y(k,273) +.280_r8*rxt(k,479)*y(k,275))*y(k,258) &
                  + (.750_r8*rxt(k,400)*y(k,252) +.380_r8*rxt(k,401)*y(k,253)) &
                 *y(k,289) + (rxt(k,488)*y(k,122) +.070_r8*rxt(k,490)*y(k,123)) &
                 *y(k,295) +.330_r8*rxt(k,47)*y(k,97) +.500_r8*rxt(k,53)*y(k,103) &
                  +rxt(k,59)*y(k,110) +.500_r8*rxt(k,60)*y(k,111) +.500_r8*rxt(k,61) &
                 *y(k,112) +rxt(k,62)*y(k,113) +.720_r8*rxt(k,63)*y(k,114) &
                  +.830_r8*rxt(k,459)*y(k,158)*y(k,115) +.500_r8*rxt(k,80)*y(k,133) &
                  +.560_r8*rxt(k,81)*y(k,134) +rxt(k,344)*y(k,256)
         loss(k,235) = (rxt(k,346)* y(k,295) + rxt(k,44) + rxt(k,811) &
                  + het_rates(k,77))* y(k,77)
         prod(k,235) = (.230_r8*rxt(k,503)*y(k,260) +.130_r8*rxt(k,539)*y(k,293) + &
                 rxt(k,545)*y(k,243) +.400_r8*rxt(k,559)*y(k,285) + &
                 .170_r8*rxt(k,562)*y(k,287) +.700_r8*rxt(k,565)*y(k,296) + &
                 .600_r8*rxt(k,572)*y(k,314) +.340_r8*rxt(k,578)*y(k,316) + &
                 .170_r8*rxt(k,581)*y(k,318))*y(k,147) &
                  + (.170_r8*rxt(k,459)*y(k,115) +.170_r8*rxt(k,469)*y(k,118) + &
                 .170_r8*rxt(k,484)*y(k,121) +.660_r8*rxt(k,501)*y(k,139))*y(k,158) &
                  + (.150_r8*rxt(k,421)*y(k,260) +.100_r8*rxt(k,499)*y(k,293) + &
                 .160_r8*rxt(k,558)*y(k,285) +.070_r8*rxt(k,561)*y(k,287))*y(k,258) &
                  + (.650_r8*rxt(k,321)*y(k,25) +.200_r8*rxt(k,345)*y(k,76) + &
                 .890_r8*rxt(k,347)*y(k,89))*y(k,295) +rxt(k,21)*y(k,9) &
                  +.500_r8*rxt(k,60)*y(k,111) +.500_r8*rxt(k,61)*y(k,112) &
                  +.280_r8*rxt(k,63)*y(k,114) +.700_r8*rxt(k,87)*y(k,167) &
                  +.600_r8*rxt(k,115)*y(k,227) +.340_r8*rxt(k,116)*y(k,230) &
                  +.170_r8*rxt(k,117)*y(k,232)
         loss(k,301) = (rxt(k,184)* y(k,158) + (rxt(k,178) +rxt(k,179) +rxt(k,180)) &
                 * y(k,258) + rxt(k,181) + het_rates(k,78))* y(k,78)
         prod(k,301) = (rxt(k,185)*y(k,79) +rxt(k,188)*y(k,157) +rxt(k,206)*y(k,136) + &
                 rxt(k,301)*y(k,43) +rxt(k,796)*y(k,161) +rxt(k,802)*y(k,172) + &
                 rxt(k,807)*y(k,174))*y(k,295) + (rxt(k,168)*y(k,294) + &
                 rxt(k,176)*y(k,157) +rxt(k,220)*y(k,57) +rxt(k,276)*y(k,75))*y(k,79) &
                  + (rxt(k,38) +.330_r8*rxt(k,39) +rxt(k,316)*y(k,294))*y(k,55) &
                  + (rxt(k,143) +rxt(k,274)*y(k,294))*y(k,83) + (rxt(k,147) + &
                 rxt(k,251)*y(k,294))*y(k,87) + (rxt(k,2) +2.000_r8*rxt(k,3))*y(k,319) &
                  +2.000_r8*rxt(k,31)*y(k,43) +rxt(k,37)*y(k,54) +rxt(k,148)*y(k,91)
         loss(k,251) = (rxt(k,220)* y(k,57) +rxt(k,276)* y(k,75) +rxt(k,176)* y(k,157) &
                  +rxt(k,168)* y(k,294) +rxt(k,185)* y(k,295) + het_rates(k,79)) &
                 * y(k,79)
         prod(k,251) = (1.440_r8*rxt(k,39) +rxt(k,317)*y(k,294))*y(k,55) +rxt(k,32) &
                 *y(k,43) +rxt(k,178)*y(k,258)*y(k,78) +rxt(k,1)*y(k,319)
         loss(k,62) = (rxt(k,272)* y(k,294) + rxt(k,142) + het_rates(k,80))* y(k,80)
         prod(k,62) = 0._r8
         loss(k,237) = (rxt(k,221)* y(k,57) +rxt(k,177)* y(k,157) +rxt(k,186) &
                 * y(k,295) + rxt(k,4) + het_rates(k,81))* y(k,81)
         prod(k,237) = (.660_r8*rxt(k,459)*y(k,115) +.660_r8*rxt(k,469)*y(k,118) + &
                 .660_r8*rxt(k,484)*y(k,121) +.030_r8*rxt(k,486)*y(k,109) + &
                 .660_r8*rxt(k,501)*y(k,139) +.220_r8*rxt(k,602)*y(k,4) + &
                 .170_r8*rxt(k,622)*y(k,7) +.320_r8*rxt(k,642)*y(k,17) + &
                 .330_r8*rxt(k,662)*y(k,125) +.020_r8*rxt(k,762)*y(k,212) + &
                 .040_r8*rxt(k,767)*y(k,213))*y(k,158) +rxt(k,192)*y(k,258)*y(k,258) &
                  +rxt(k,191)*y(k,295)*y(k,295)
         loss(k,69) = ( + rxt(k,153) + het_rates(k,82))* y(k,82)
         prod(k,69) =rxt(k,809)*y(k,319)*y(k,176)
         loss(k,208) = (rxt(k,267)* y(k,157) + (rxt(k,273) +rxt(k,274))* y(k,294) &
                  +rxt(k,268)* y(k,295) + rxt(k,143) + het_rates(k,83))* y(k,83)
         prod(k,208) = (rxt(k,254)*y(k,43) +rxt(k,255)*y(k,258))*y(k,18)
         loss(k,91) = (rxt(k,295)* y(k,294) +rxt(k,290)* y(k,295) + rxt(k,144) &
                  + het_rates(k,84))* y(k,84)
         prod(k,91) = 0._r8
         loss(k,93) = (rxt(k,296)* y(k,294) +rxt(k,291)* y(k,295) + rxt(k,145) &
                  + het_rates(k,85))* y(k,85)
         prod(k,93) = 0._r8
         loss(k,108) = (rxt(k,297)* y(k,294) +rxt(k,292)* y(k,295) + rxt(k,146) &
                  + het_rates(k,86))* y(k,86)
         prod(k,108) = 0._r8
         loss(k,304) = ((rxt(k,890) +rxt(k,895) +rxt(k,900))* y(k,61) + (rxt(k,892) + &
                 rxt(k,897))* y(k,95) + (rxt(k,885) +rxt(k,891) +rxt(k,896))* y(k,96) &
                  +rxt(k,238)* y(k,157) + (rxt(k,250) +rxt(k,251))* y(k,294) &
                  +rxt(k,239)* y(k,295) + rxt(k,147) + het_rates(k,87))* y(k,87)
         prod(k,304) = (rxt(k,218)*y(k,43) +rxt(k,219)*y(k,55) +rxt(k,220)*y(k,79) + &
                 rxt(k,221)*y(k,81) +rxt(k,222)*y(k,258) +rxt(k,240)*y(k,96) + &
                 rxt(k,281)*y(k,42) +rxt(k,283)*y(k,44) +2.000_r8*rxt(k,286)*y(k,47) + &
                 rxt(k,288)*y(k,56) +rxt(k,330)*y(k,29))*y(k,57) +rxt(k,237)*y(k,295) &
                 *y(k,60)
         loss(k,79) = (rxt(k,318)* y(k,294) +rxt(k,309)* y(k,295) + het_rates(k,88)) &
                 * y(k,88)
         prod(k,79) = 0._r8
         loss(k,182) = (rxt(k,347)* y(k,295) + rxt(k,45) + het_rates(k,89))* y(k,89)
         prod(k,182) = (.680_r8*rxt(k,482)*y(k,258) +.810_r8*rxt(k,527)*y(k,147)) &
                 *y(k,276) +.700_r8*rxt(k,484)*y(k,158)*y(k,121)
         loss(k,227) = (rxt(k,310)* y(k,295) + het_rates(k,90))* y(k,90)
         prod(k,227) = (.370_r8*rxt(k,323)*y(k,26) +.120_r8*rxt(k,354)*y(k,30) + &
                 .330_r8*rxt(k,385)*y(k,126) +.120_r8*rxt(k,403)*y(k,132) + &
                 .220_r8*rxt(k,486)*y(k,109) +.080_r8*rxt(k,642)*y(k,17) + &
                 .150_r8*rxt(k,762)*y(k,212) +.260_r8*rxt(k,767)*y(k,213))*y(k,158) &
                  + (.500_r8*rxt(k,311)*y(k,92) +.350_r8*rxt(k,321)*y(k,25) + &
                 .400_r8*rxt(k,422)*y(k,107))*y(k,295) &
                  + (.500_r8*rxt(k,312)*y(k,258) +rxt(k,314)*y(k,147))*y(k,259) &
                  +.410_r8*rxt(k,48)*y(k,98)
         loss(k,107) = ( + rxt(k,148) + het_rates(k,91))* y(k,91)
         prod(k,107) = (rxt(k,275)*y(k,55) +rxt(k,276)*y(k,79) +rxt(k,277)*y(k,319) + &
                 rxt(k,278)*y(k,93))*y(k,75)
         loss(k,213) = (rxt(k,311)* y(k,295) + rxt(k,46) + het_rates(k,92))* y(k,92)
         prod(k,213) = (.330_r8*rxt(k,486)*y(k,109) +.110_r8*rxt(k,642)*y(k,17) + &
                 .230_r8*rxt(k,762)*y(k,212) +.400_r8*rxt(k,767)*y(k,213))*y(k,158) &
                  +.500_r8*rxt(k,312)*y(k,259)*y(k,258)
         loss(k,302) = (rxt(k,278)* y(k,75) +rxt(k,215)* y(k,295) + rxt(k,9) &
                  + het_rates(k,93))* y(k,93)
         prod(k,302) = (rxt(k,831) +rxt(k,299)*y(k,43) +rxt(k,332)*y(k,46) + &
                 rxt(k,360)*y(k,50) +rxt(k,708)*y(k,203) +rxt(k,727)*y(k,205) + &
                 rxt(k,755)*y(k,202) +rxt(k,793)*y(k,69))*y(k,149) + (rxt(k,888) + &
                 rxt(k,893) +rxt(k,899) +rxt(k,890)*y(k,87) +rxt(k,895)*y(k,87) + &
                 rxt(k,900)*y(k,87))*y(k,61) + (2.000_r8*rxt(k,827) + &
                 2.000_r8*rxt(k,884) +2.000_r8*rxt(k,887) +2.000_r8*rxt(k,898)) &
                 *y(k,138) + (rxt(k,886) +rxt(k,889) +rxt(k,894))*y(k,21) &
                  + (.500_r8*rxt(k,830) +rxt(k,214)*y(k,295))*y(k,148) +rxt(k,813) &
                 *y(k,97) +rxt(k,816)*y(k,107) +rxt(k,817)*y(k,108) +rxt(k,819) &
                 *y(k,110) +rxt(k,820)*y(k,111) +rxt(k,824)*y(k,115) +rxt(k,825) &
                 *y(k,116) +rxt(k,826)*y(k,118) +rxt(k,818)*y(k,121) +rxt(k,828) &
                 *y(k,139) +rxt(k,832)*y(k,162) +rxt(k,835)*y(k,214) +rxt(k,838) &
                 *y(k,219) +rxt(k,837)*y(k,220) +rxt(k,840)*y(k,223) +rxt(k,839) &
                 *y(k,224)
         loss(k,128) = (rxt(k,193)* y(k,295) + rxt(k,10) + rxt(k,11) + rxt(k,216) &
                  + het_rates(k,94))* y(k,94)
         prod(k,128) =rxt(k,212)*y(k,258)*y(k,148)
         loss(k,195) = ((rxt(k,892) +rxt(k,897))* y(k,87) +rxt(k,269)* y(k,157) &
                  + rxt(k,149) + het_rates(k,95))* y(k,95)
         prod(k,195) = (rxt(k,886) +rxt(k,889) +rxt(k,894))*y(k,21) &
                  +rxt(k,261)*y(k,258)*y(k,20)
         loss(k,209) = (rxt(k,240)* y(k,57) + (rxt(k,885) +rxt(k,891) +rxt(k,896)) &
                 * y(k,87) +rxt(k,241)* y(k,157) +rxt(k,242)* y(k,295) + rxt(k,150) &
                  + het_rates(k,96))* y(k,96)
         prod(k,209) = (rxt(k,888) +rxt(k,893) +rxt(k,899) +rxt(k,234)*y(k,295)) &
                 *y(k,61) +rxt(k,229)*y(k,258)*y(k,60)
         loss(k,193) = (rxt(k,379)* y(k,295) + rxt(k,47) + rxt(k,813) &
                  + het_rates(k,97))* y(k,97)
         prod(k,193) =rxt(k,378)*y(k,255)*y(k,147)
         loss(k,153) = (rxt(k,414)* y(k,295) + rxt(k,48) + het_rates(k,98))* y(k,98)
         prod(k,153) =.250_r8*rxt(k,529)*y(k,278)
         loss(k,154) = (rxt(k,415)* y(k,295) + rxt(k,49) + het_rates(k,99))* y(k,99)
         prod(k,154) =.250_r8*rxt(k,532)*y(k,279)
         loss(k,136) = (rxt(k,416)* y(k,295) + rxt(k,50) + het_rates(k,100))* y(k,100)
         prod(k,136) =.090_r8*rxt(k,489)*y(k,295)*y(k,123) +.150_r8*rxt(k,529) &
                 *y(k,278)
         loss(k,140) = (rxt(k,417)* y(k,295) + rxt(k,51) + het_rates(k,101))* y(k,101)
         prod(k,140) =.090_r8*rxt(k,489)*y(k,295)*y(k,123) +.150_r8*rxt(k,532) &
                 *y(k,279)
         loss(k,258) = (rxt(k,362)* y(k,295) + rxt(k,52) + het_rates(k,102))* y(k,102)
         prod(k,258) = (.500_r8*rxt(k,367)*y(k,170) +.500_r8*rxt(k,380)*y(k,127) + &
                 rxt(k,387)*y(k,128) +.250_r8*rxt(k,398)*y(k,131) + &
                 .220_r8*rxt(k,418)*y(k,103) +.500_r8*rxt(k,419)*y(k,105) + &
                 .190_r8*rxt(k,423)*y(k,108) +.280_r8*rxt(k,456)*y(k,114) + &
                 rxt(k,488)*y(k,122) +.070_r8*rxt(k,490)*y(k,123))*y(k,295) &
                  + (.290_r8*rxt(k,503)*y(k,260) +.730_r8*rxt(k,515)*y(k,270) + &
                 .870_r8*rxt(k,519)*y(k,272) +.330_r8*rxt(k,521)*y(k,273) + &
                 .070_r8*rxt(k,525)*y(k,275) +.860_r8*rxt(k,535)*y(k,284))*y(k,147) &
                  + (.510_r8*rxt(k,383)*y(k,284) +.190_r8*rxt(k,421)*y(k,260) + &
                 .460_r8*rxt(k,461)*y(k,270) +.440_r8*rxt(k,467)*y(k,272) + &
                 .150_r8*rxt(k,471)*y(k,273) +.060_r8*rxt(k,479)*y(k,275))*y(k,258) &
                  + (rxt(k,384) +.860_r8*rxt(k,381)*y(k,252) + &
                 .900_r8*rxt(k,382)*y(k,253))*y(k,284) &
                  + (.830_r8*rxt(k,469)*y(k,118) +.180_r8*rxt(k,682)*y(k,135)) &
                 *y(k,158) +.170_r8*rxt(k,47)*y(k,97) +.500_r8*rxt(k,53)*y(k,103) &
                  +rxt(k,59)*y(k,110) +.500_r8*rxt(k,60)*y(k,111) +.500_r8*rxt(k,61) &
                 *y(k,112) +rxt(k,62)*y(k,113) +.280_r8*rxt(k,63)*y(k,114) &
                  +.500_r8*rxt(k,74)*y(k,127) +.860_r8*rxt(k,75)*y(k,128) &
                  +.200_r8*rxt(k,368)*y(k,298)*y(k,253)
         loss(k,263) = (rxt(k,418)* y(k,295) + rxt(k,53) + rxt(k,54) &
                  + het_rates(k,103))* y(k,103)
         prod(k,263) = (.250_r8*rxt(k,431)*y(k,262) +.470_r8*rxt(k,445)*y(k,267) + &
                 .470_r8*rxt(k,449)*y(k,268) +.470_r8*rxt(k,492)*y(k,278) + &
                 .470_r8*rxt(k,496)*y(k,279))*y(k,253) &
                  + (.450_r8*rxt(k,509)*y(k,267) +.450_r8*rxt(k,511)*y(k,268) + &
                 .450_r8*rxt(k,530)*y(k,278) +.450_r8*rxt(k,533)*y(k,279))*y(k,147) &
                  + (.450_r8*rxt(k,444)*y(k,267) +.450_r8*rxt(k,448)*y(k,268) + &
                 .450_r8*rxt(k,491)*y(k,278) +.450_r8*rxt(k,495)*y(k,279))*y(k,252) &
                  +.450_r8*rxt(k,64)*y(k,115) +.450_r8*rxt(k,67)*y(k,118) &
                  +.130_r8*rxt(k,489)*y(k,295)*y(k,123) +rxt(k,82)*y(k,139)
         loss(k,187) = (rxt(k,363)* y(k,295) + rxt(k,55) + rxt(k,56) &
                  + het_rates(k,104))* y(k,104)
         prod(k,187) = (.500_r8*rxt(k,41) +rxt(k,376)*y(k,295))*y(k,68) &
                  + (.120_r8*rxt(k,482)*y(k,258) +.150_r8*rxt(k,527)*y(k,147)) &
                 *y(k,276) +.150_r8*rxt(k,415)*y(k,295)*y(k,99) &
                  +.130_r8*rxt(k,484)*y(k,158)*y(k,121)
         loss(k,185) = (rxt(k,419)* y(k,295) + rxt(k,814) + het_rates(k,105)) &
                 * y(k,105)
         prod(k,185) = (.080_r8*rxt(k,414)*y(k,98) +.180_r8*rxt(k,415)*y(k,99) + &
                 .580_r8*rxt(k,416)*y(k,100) +.770_r8*rxt(k,417)*y(k,101) + &
                 .190_r8*rxt(k,420)*y(k,106) +.040_r8*rxt(k,502)*y(k,139))*y(k,295) &
                  +rxt(k,57)*y(k,107) +rxt(k,58)*y(k,108)
         loss(k,241) = (rxt(k,420)* y(k,295) + rxt(k,815) + het_rates(k,106)) &
                 * y(k,106)
         prod(k,241) = (.080_r8*rxt(k,460)*y(k,115) +.150_r8*rxt(k,463)*y(k,116) + &
                 .130_r8*rxt(k,466)*y(k,117) +.040_r8*rxt(k,470)*y(k,118) + &
                 .070_r8*rxt(k,485)*y(k,121) +.850_r8*rxt(k,490)*y(k,123))*y(k,295)
         loss(k,138) = (rxt(k,422)* y(k,295) + rxt(k,57) + rxt(k,816) &
                  + het_rates(k,107))* y(k,107)
         prod(k,138) = (.200_r8*rxt(k,422)*y(k,107) +.400_r8*rxt(k,481)*y(k,120)) &
                 *y(k,295)
         loss(k,218) = (rxt(k,423)* y(k,295) + rxt(k,58) + rxt(k,817) &
                  + het_rates(k,108))* y(k,108)
         prod(k,218) = (.060_r8*rxt(k,423)*y(k,108) +.030_r8*rxt(k,472)*y(k,119) + &
                 .200_r8*rxt(k,485)*y(k,121))*y(k,295)
         loss(k,229) = (rxt(k,473)* y(k,149) +rxt(k,486)* y(k,158) +rxt(k,487) &
                 * y(k,295) + het_rates(k,109))* y(k,109)
         prod(k,229) = 0._r8
         loss(k,250) = (rxt(k,453)* y(k,295) + rxt(k,59) + rxt(k,819) &
                  + het_rates(k,110))* y(k,110)
         prod(k,250) = (rxt(k,514)*y(k,269) +rxt(k,516)*y(k,270) + &
                 rxt(k,518)*y(k,271) +rxt(k,520)*y(k,272) +rxt(k,522)*y(k,273) + &
                 rxt(k,524)*y(k,274) +rxt(k,526)*y(k,275) +rxt(k,528)*y(k,276)) &
                 *y(k,147)
         loss(k,204) = (rxt(k,452)* y(k,295) + rxt(k,60) + rxt(k,820) &
                  + het_rates(k,111))* y(k,111)
         prod(k,204) =rxt(k,453)*y(k,295)*y(k,110) +rxt(k,540)*y(k,293)*y(k,147)
         loss(k,260) = (rxt(k,454)* y(k,295) + rxt(k,61) + rxt(k,821) &
                  + het_rates(k,112))* y(k,112)
         prod(k,260) =rxt(k,455)*y(k,295)*y(k,113) +rxt(k,504)*y(k,260)*y(k,147) &
                  +rxt(k,462)*y(k,270) +rxt(k,465)*y(k,271)
         loss(k,234) = (rxt(k,455)* y(k,295) + rxt(k,62) + rxt(k,822) &
                  + het_rates(k,113))* y(k,113)
         prod(k,234) = (.420_r8*rxt(k,457)*y(k,269) +.480_r8*rxt(k,461)*y(k,270) + &
                 .400_r8*rxt(k,464)*y(k,271) +.500_r8*rxt(k,467)*y(k,272) + &
                 .600_r8*rxt(k,471)*y(k,273) +.490_r8*rxt(k,479)*y(k,275) + &
                 .170_r8*rxt(k,482)*y(k,276) +.200_r8*rxt(k,499)*y(k,293))*y(k,258) &
                  +rxt(k,458)*y(k,269) +rxt(k,468)*y(k,272) +rxt(k,480)*y(k,275) &
                  +rxt(k,483)*y(k,276)
         loss(k,170) = (rxt(k,456)* y(k,295) + rxt(k,63) + rxt(k,823) &
                  + het_rates(k,114))* y(k,114)
         prod(k,170) =.080_r8*rxt(k,490)*y(k,295)*y(k,123) &
                  +.350_r8*rxt(k,421)*y(k,260)*y(k,258)
         loss(k,272) = (rxt(k,459)* y(k,158) +rxt(k,460)* y(k,295) + rxt(k,64) &
                  + rxt(k,824) + het_rates(k,115))* y(k,115)
         prod(k,272) = (rxt(k,512)*y(k,268) +rxt(k,534)*y(k,279))*y(k,147) &
                  + (.280_r8*rxt(k,475)*y(k,253) +.530_r8*rxt(k,477)*y(k,274)) &
                 *y(k,274)
         loss(k,159) = (rxt(k,463)* y(k,295) + rxt(k,65) + rxt(k,825) &
                  + het_rates(k,116))* y(k,116)
         prod(k,159) =rxt(k,506)*y(k,261)*y(k,147)
         loss(k,150) = (rxt(k,466)* y(k,295) + rxt(k,66) + het_rates(k,117))* y(k,117)
         prod(k,150) =rxt(k,508)*y(k,262)*y(k,147)
         loss(k,273) = (rxt(k,469)* y(k,158) +rxt(k,470)* y(k,295) + rxt(k,67) &
                  + rxt(k,826) + het_rates(k,118))* y(k,118)
         prod(k,273) = (rxt(k,510)*y(k,267) +rxt(k,531)*y(k,278))*y(k,147) &
                  + (.050_r8*rxt(k,475)*y(k,253) +.090_r8*rxt(k,477)*y(k,274)) &
                 *y(k,274)
         loss(k,166) = (rxt(k,472)* y(k,295) + rxt(k,68) + het_rates(k,119))* y(k,119)
         prod(k,166) = (.070_r8*rxt(k,475)*y(k,253) +.150_r8*rxt(k,477)*y(k,274)) &
                 *y(k,274)
         loss(k,212) = (rxt(k,481)* y(k,295) + rxt(k,69) + het_rates(k,120))* y(k,120)
         prod(k,212) =.230_r8*rxt(k,476)*y(k,274)*y(k,258)
         loss(k,245) = (rxt(k,484)* y(k,158) +rxt(k,485)* y(k,295) + rxt(k,70) &
                  + rxt(k,818) + het_rates(k,121))* y(k,121)
         prod(k,245) =.530_r8*rxt(k,476)*y(k,274)*y(k,258)
         loss(k,186) = (rxt(k,488)* y(k,295) + het_rates(k,122))* y(k,122)
         prod(k,186) = (.250_r8*rxt(k,425)*y(k,261) +.250_r8*rxt(k,431)*y(k,262) + &
                 .250_r8*rxt(k,445)*y(k,267) +.250_r8*rxt(k,449)*y(k,268) + &
                 .250_r8*rxt(k,492)*y(k,278) +.250_r8*rxt(k,496)*y(k,279))*y(k,253)
         loss(k,264) = ((rxt(k,489) +rxt(k,490))* y(k,295) + rxt(k,71) &
                  + het_rates(k,123))* y(k,123)
         prod(k,264) = (.940_r8*rxt(k,426)*y(k,261) +.940_r8*rxt(k,432)*y(k,262) + &
                 rxt(k,446)*y(k,267) +rxt(k,450)*y(k,268) +rxt(k,493)*y(k,278) + &
                 rxt(k,497)*y(k,279))*y(k,258)
         loss(k,53) = (rxt(k,866)* y(k,295) + het_rates(k,124))* y(k,124)
         prod(k,53) = 0._r8
         loss(k,201) = (rxt(k,644)* y(k,149) +rxt(k,662)* y(k,158) +rxt(k,663) &
                 * y(k,295) + het_rates(k,125))* y(k,125)
         prod(k,201) = 0._r8
         loss(k,269) = (rxt(k,385)* y(k,158) +rxt(k,386)* y(k,295) + rxt(k,72) &
                  + rxt(k,73) + het_rates(k,126))* y(k,126)
         prod(k,269) = (.040_r8*rxt(k,474)*y(k,252) +.020_r8*rxt(k,475)*y(k,253) + &
                 .020_r8*rxt(k,476)*y(k,258) +.160_r8*rxt(k,477)*y(k,274) + &
                 .040_r8*rxt(k,478)*y(k,149) +.040_r8*rxt(k,523)*y(k,147))*y(k,274) &
                  + (rxt(k,433) +rxt(k,430)*y(k,252) +.500_r8*rxt(k,431)*y(k,253) + &
                 .060_r8*rxt(k,432)*y(k,258) +rxt(k,507)*y(k,147))*y(k,262) &
                  + (rxt(k,51) +.140_r8*rxt(k,417)*y(k,295))*y(k,101) &
                  +.350_r8*rxt(k,415)*y(k,295)*y(k,99) +.410_r8*rxt(k,486)*y(k,158) &
                 *y(k,109) +rxt(k,66)*y(k,117) +.500_r8*rxt(k,68)*y(k,119) &
                  +.120_r8*rxt(k,69)*y(k,120) +.300_r8*rxt(k,71)*y(k,123)
         loss(k,259) = (rxt(k,380)* y(k,295) + rxt(k,74) + het_rates(k,127))* y(k,127)
         prod(k,259) = (.060_r8*rxt(k,513)*y(k,269) +.270_r8*rxt(k,515)*y(k,270) + &
                 .210_r8*rxt(k,521)*y(k,273) +.490_r8*rxt(k,525)*y(k,275) + &
                 .020_r8*rxt(k,527)*y(k,276) +rxt(k,536)*y(k,284) + &
                 .390_r8*rxt(k,539)*y(k,293))*y(k,147) &
                  + (.030_r8*rxt(k,457)*y(k,269) +.060_r8*rxt(k,461)*y(k,270) + &
                 .060_r8*rxt(k,471)*y(k,273) +.150_r8*rxt(k,479)*y(k,275) + &
                 .020_r8*rxt(k,482)*y(k,276) +.290_r8*rxt(k,499)*y(k,293))*y(k,258) &
                  + (.500_r8*rxt(k,452)*y(k,111) +.250_r8*rxt(k,454)*y(k,112) + &
                 .060_r8*rxt(k,460)*y(k,115) +.240_r8*rxt(k,502)*y(k,139))*y(k,295) &
                  +.510_r8*rxt(k,500)*y(k,293)
         loss(k,233) = (rxt(k,387)* y(k,295) + rxt(k,75) + het_rates(k,128))* y(k,128)
         prod(k,233) = (.550_r8*rxt(k,448)*y(k,252) +.280_r8*rxt(k,449)*y(k,253) + &
                 .550_r8*rxt(k,511)*y(k,147))*y(k,268) &
                  + (.550_r8*rxt(k,495)*y(k,252) +.280_r8*rxt(k,496)*y(k,253) + &
                 .550_r8*rxt(k,533)*y(k,147))*y(k,279) &
                  + (.090_r8*rxt(k,417)*y(k,101) +.250_r8*rxt(k,454)*y(k,112)) &
                 *y(k,295) +.550_r8*rxt(k,64)*y(k,115) +.410_r8*rxt(k,383)*y(k,284) &
                 *y(k,258)
         loss(k,145) = (rxt(k,396)* y(k,295) + rxt(k,76) + het_rates(k,129))* y(k,129)
         prod(k,145) =.800_r8*rxt(k,19)*y(k,1) +.800_r8*rxt(k,20)*y(k,2) &
                  +.800_r8*rxt(k,410)*y(k,236)*y(k,147)
         loss(k,109) = (rxt(k,397)* y(k,295) + rxt(k,77) + het_rates(k,130))* y(k,130)
         prod(k,109) =.800_r8*rxt(k,394)*y(k,288)*y(k,258)
         loss(k,141) = (rxt(k,398)* y(k,295) + rxt(k,78) + rxt(k,407) &
                  + het_rates(k,131))* y(k,131)
         prod(k,141) =rxt(k,406)*y(k,286)*y(k,148)
         loss(k,270) = (rxt(k,403)* y(k,158) +rxt(k,404)* y(k,295) + rxt(k,79) &
                  + het_rates(k,132))* y(k,132)
         prod(k,270) = (rxt(k,427) +rxt(k,424)*y(k,252) +.750_r8*rxt(k,425)*y(k,253) + &
                 .060_r8*rxt(k,426)*y(k,258) +rxt(k,505)*y(k,147))*y(k,261) &
                  + (.420_r8*rxt(k,474)*y(k,252) +.050_r8*rxt(k,475)*y(k,253) + &
                 .220_r8*rxt(k,476)*y(k,258) +.420_r8*rxt(k,478)*y(k,149) + &
                 .420_r8*rxt(k,523)*y(k,147))*y(k,274) + (rxt(k,50) + &
                 .230_r8*rxt(k,416)*y(k,295))*y(k,100) +.350_r8*rxt(k,414)*y(k,295) &
                 *y(k,98) +.170_r8*rxt(k,486)*y(k,158)*y(k,109) +rxt(k,65)*y(k,116) &
                  +.500_r8*rxt(k,68)*y(k,119) +.880_r8*rxt(k,69)*y(k,120) &
                  +.700_r8*rxt(k,71)*y(k,123)
         loss(k,265) = (rxt(k,399)* y(k,295) + rxt(k,80) + het_rates(k,133))* y(k,133)
         prod(k,265) = (rxt(k,517)*y(k,271) +.130_r8*rxt(k,519)*y(k,272) + &
                 .120_r8*rxt(k,521)*y(k,273) +.040_r8*rxt(k,525)*y(k,275) + &
                 .020_r8*rxt(k,527)*y(k,276) +rxt(k,538)*y(k,289) + &
                 .360_r8*rxt(k,539)*y(k,293))*y(k,147) &
                  + (.600_r8*rxt(k,464)*y(k,271) +.060_r8*rxt(k,467)*y(k,272) + &
                 .040_r8*rxt(k,471)*y(k,273) +.020_r8*rxt(k,479)*y(k,275) + &
                 .010_r8*rxt(k,482)*y(k,276) +.310_r8*rxt(k,499)*y(k,293))*y(k,258) &
                  + (.050_r8*rxt(k,423)*y(k,108) +.500_r8*rxt(k,452)*y(k,111) + &
                 .250_r8*rxt(k,454)*y(k,112) +.040_r8*rxt(k,470)*y(k,118) + &
                 .040_r8*rxt(k,502)*y(k,139))*y(k,295) +.490_r8*rxt(k,500)*y(k,293)
         loss(k,239) = (rxt(k,405)* y(k,295) + rxt(k,81) + het_rates(k,134))* y(k,134)
         prod(k,239) = (.550_r8*rxt(k,444)*y(k,252) +.280_r8*rxt(k,445)*y(k,253) + &
                 .550_r8*rxt(k,509)*y(k,147))*y(k,267) &
                  + (.550_r8*rxt(k,491)*y(k,252) +.280_r8*rxt(k,492)*y(k,253) + &
                 .550_r8*rxt(k,530)*y(k,147))*y(k,278) &
                  + (.190_r8*rxt(k,416)*y(k,100) +.250_r8*rxt(k,454)*y(k,112)) &
                 *y(k,295) +.550_r8*rxt(k,67)*y(k,118) +.460_r8*rxt(k,402)*y(k,289) &
                 *y(k,258)
         loss(k,177) = (rxt(k,664)* y(k,149) +rxt(k,682)* y(k,158) +rxt(k,683) &
                 * y(k,295) + het_rates(k,135))* y(k,135)
         prod(k,177) = 0._r8
         loss(k,134) = (rxt(k,194)* y(k,147) + (rxt(k,195) +rxt(k,196) +rxt(k,197)) &
                 * y(k,148) +rxt(k,206)* y(k,295) + rxt(k,198) + het_rates(k,136)) &
                 * y(k,136)
         prod(k,134) =rxt(k,15)*y(k,147)
         loss(k,80) = ((rxt(k,210) +rxt(k,211))* y(k,294) + rxt(k,12) &
                  + het_rates(k,137))* y(k,137)
         prod(k,80) =rxt(k,195)*y(k,148)*y(k,136)
         loss(k,102) = ( + rxt(k,13) + rxt(k,14) + rxt(k,217) + rxt(k,827) &
                  + rxt(k,884) + rxt(k,887) + rxt(k,898) + het_rates(k,138))* y(k,138)
         prod(k,102) =rxt(k,213)*y(k,149)*y(k,148)
         loss(k,274) = (rxt(k,501)* y(k,158) +rxt(k,502)* y(k,295) + rxt(k,82) &
                  + rxt(k,828) + het_rates(k,139))* y(k,139)
         prod(k,274) = (.540_r8*rxt(k,474)*y(k,252) +.530_r8*rxt(k,475)*y(k,253) + &
                 1.070_r8*rxt(k,477)*y(k,274) +.540_r8*rxt(k,478)*y(k,149) + &
                 .540_r8*rxt(k,523)*y(k,147))*y(k,274) &
                  + (.040_r8*rxt(k,460)*y(k,115) +.030_r8*rxt(k,470)*y(k,118) + &
                 .050_r8*rxt(k,472)*y(k,119) +.020_r8*rxt(k,481)*y(k,120) + &
                 .090_r8*rxt(k,485)*y(k,121))*y(k,295) +rxt(k,70)*y(k,121)
         loss(k,10) = ( + het_rates(k,140))* y(k,140)
         prod(k,10) = 0._r8
         loss(k,11) = ( + het_rates(k,141))* y(k,141)
         prod(k,11) = 0._r8
         loss(k,12) = ( + het_rates(k,142))* y(k,142)
         prod(k,12) = 0._r8
         loss(k,59) = (rxt(k,810)* y(k,295) + het_rates(k,143))* y(k,143)
         prod(k,59) = 0._r8
         loss(k,13) = ( + rxt(k,829) + het_rates(k,144))* y(k,144)
         prod(k,13) = 0._r8
         loss(k,14) = ( + rxt(k,903) + het_rates(k,145))* y(k,145)
         prod(k,14) = 0._r8
         loss(k,15) = ( + rxt(k,902) + het_rates(k,146))* y(k,146)
         prod(k,15) = 0._r8
         loss(k,303) = (rxt(k,262)* y(k,20) +rxt(k,230)* y(k,60) +rxt(k,194)* y(k,136) &
                  +rxt(k,203)* y(k,149) +rxt(k,209)* y(k,157) +rxt(k,208)* y(k,158) &
                  +rxt(k,542)* y(k,235) + (rxt(k,410) +rxt(k,411))* y(k,236) &
                  +rxt(k,589)* y(k,237) +rxt(k,597)* y(k,238) +rxt(k,609)* y(k,240) &
                  +rxt(k,617)* y(k,241) +rxt(k,545)* y(k,243) +rxt(k,629)* y(k,245) &
                  +rxt(k,637)* y(k,246) +rxt(k,550)* y(k,248) +rxt(k,327)* y(k,249) &
                  +rxt(k,357)* y(k,250) +rxt(k,552)* y(k,251) +rxt(k,338)* y(k,252) &
                  +rxt(k,305)* y(k,253) +rxt(k,556)* y(k,254) + (rxt(k,377) + &
                 rxt(k,378))* y(k,255) +rxt(k,342)* y(k,257) +rxt(k,207)* y(k,258) &
                  +rxt(k,314)* y(k,259) + (rxt(k,503) +rxt(k,504))* y(k,260) &
                  + (rxt(k,505) +rxt(k,506))* y(k,261) + (rxt(k,507) +rxt(k,508)) &
                 * y(k,262) + (rxt(k,509) +rxt(k,510))* y(k,267) + (rxt(k,511) + &
                 rxt(k,512))* y(k,268) + (rxt(k,513) +rxt(k,514))* y(k,269) &
                  + (rxt(k,515) +rxt(k,516))* y(k,270) + (rxt(k,517) +rxt(k,518)) &
                 * y(k,271) + (rxt(k,519) +rxt(k,520))* y(k,272) + (rxt(k,521) + &
                 rxt(k,522))* y(k,273) + (rxt(k,523) +rxt(k,524))* y(k,274) &
                  + (rxt(k,525) +rxt(k,526))* y(k,275) + (rxt(k,527) +rxt(k,528)) &
                 * y(k,276) + (rxt(k,530) +rxt(k,531))* y(k,278) + (rxt(k,533) + &
                 rxt(k,534))* y(k,279) +rxt(k,649)* y(k,281) +rxt(k,657)* y(k,282) &
                  + (rxt(k,535) +rxt(k,536))* y(k,284) +rxt(k,559)* y(k,285) &
                  +rxt(k,392)* y(k,286) +rxt(k,562)* y(k,287) +rxt(k,395)* y(k,288) &
                  + (rxt(k,537) +rxt(k,538))* y(k,289) +rxt(k,669)* y(k,290) &
                  +rxt(k,677)* y(k,291) + (rxt(k,539) +rxt(k,540))* y(k,293) &
                  +rxt(k,565)* y(k,296) +rxt(k,366)* y(k,297) +rxt(k,370)* y(k,298) &
                  +rxt(k,688)* y(k,299) +rxt(k,692)* y(k,300) +rxt(k,696)* y(k,301) &
                  +rxt(k,704)* y(k,302) +rxt(k,712)* y(k,303) +rxt(k,722)* y(k,304) &
                  +rxt(k,731)* y(k,305) +rxt(k,741)* y(k,306) +rxt(k,752)* y(k,307) &
                  +rxt(k,761)* y(k,308) +rxt(k,766)* y(k,309) +rxt(k,773)* y(k,310) &
                  +rxt(k,777)* y(k,311) +rxt(k,781)* y(k,312) +rxt(k,785)* y(k,313) &
                  +rxt(k,572)* y(k,314) +rxt(k,578)* y(k,316) +rxt(k,581)* y(k,318) &
                  + rxt(k,15) + het_rates(k,147))* y(k,147)
         prod(k,303) = (rxt(k,16) +.500_r8*rxt(k,830) +2.000_r8*rxt(k,196)*y(k,136) + &
                 rxt(k,199)*y(k,157) +rxt(k,803)*y(k,174))*y(k,148) + (rxt(k,198) + &
                 rxt(k,206)*y(k,295))*y(k,136) +2.000_r8*rxt(k,210)*y(k,294)*y(k,137) &
                  +rxt(k,14)*y(k,138) +rxt(k,17)*y(k,149)
         loss(k,312) = (rxt(k,263)* y(k,20) +rxt(k,232)* y(k,60) + (rxt(k,195) + &
                 rxt(k,196) +rxt(k,197))* y(k,136) +rxt(k,213)* y(k,149) &
                  + (rxt(k,199) +rxt(k,201))* y(k,157) +rxt(k,200)* y(k,158) &
                  +rxt(k,567)* y(k,165) +rxt(k,803)* y(k,174) +rxt(k,570)* y(k,235) &
                  +rxt(k,351)* y(k,252) +rxt(k,557)* y(k,254) +rxt(k,212)* y(k,258) &
                  +rxt(k,560)* y(k,285) +rxt(k,406)* y(k,286) +rxt(k,563)* y(k,287) &
                  +rxt(k,214)* y(k,295) +rxt(k,684)* y(k,302) +rxt(k,685)* y(k,304) &
                  +rxt(k,686)* y(k,307) + rxt(k,16) + rxt(k,830) + het_rates(k,148)) &
                 * y(k,148)
         prod(k,312) = (2.000_r8*rxt(k,203)*y(k,149) +rxt(k,207)*y(k,258) + &
                 rxt(k,208)*y(k,158) +rxt(k,209)*y(k,157) +rxt(k,230)*y(k,60) + &
                 rxt(k,262)*y(k,20) +rxt(k,305)*y(k,253) +rxt(k,314)*y(k,259) + &
                 rxt(k,327)*y(k,249) +rxt(k,338)*y(k,252) +rxt(k,342)*y(k,257) + &
                 rxt(k,357)*y(k,250) +rxt(k,366)*y(k,297) +rxt(k,370)*y(k,298) + &
                 rxt(k,377)*y(k,255) +rxt(k,392)*y(k,286) +rxt(k,395)*y(k,288) + &
                 rxt(k,410)*y(k,236) +rxt(k,503)*y(k,260) +rxt(k,505)*y(k,261) + &
                 rxt(k,507)*y(k,262) +rxt(k,509)*y(k,267) +rxt(k,511)*y(k,268) + &
                 rxt(k,513)*y(k,269) +1.730_r8*rxt(k,515)*y(k,270) + &
                 rxt(k,517)*y(k,271) +rxt(k,519)*y(k,272) +rxt(k,521)*y(k,273) + &
                 1.460_r8*rxt(k,523)*y(k,274) +rxt(k,525)*y(k,275) + &
                 rxt(k,527)*y(k,276) +rxt(k,530)*y(k,278) +rxt(k,533)*y(k,279) + &
                 rxt(k,535)*y(k,284) +rxt(k,537)*y(k,289) +rxt(k,539)*y(k,293) + &
                 rxt(k,542)*y(k,235) +rxt(k,545)*y(k,243) +rxt(k,550)*y(k,248) + &
                 rxt(k,552)*y(k,251) +rxt(k,556)*y(k,254) +rxt(k,559)*y(k,285) + &
                 rxt(k,562)*y(k,287) +rxt(k,565)*y(k,296) +rxt(k,572)*y(k,314) + &
                 rxt(k,578)*y(k,316) +rxt(k,581)*y(k,318) + &
                 1.860_r8*rxt(k,589)*y(k,237) +.770_r8*rxt(k,597)*y(k,238) + &
                 1.860_r8*rxt(k,609)*y(k,240) +.700_r8*rxt(k,617)*y(k,241) + &
                 1.390_r8*rxt(k,629)*y(k,245) +.750_r8*rxt(k,637)*y(k,246) + &
                 1.360_r8*rxt(k,649)*y(k,281) +.770_r8*rxt(k,657)*y(k,282) + &
                 1.820_r8*rxt(k,669)*y(k,290) +.710_r8*rxt(k,677)*y(k,291) + &
                 .700_r8*rxt(k,688)*y(k,299) +.700_r8*rxt(k,692)*y(k,300) + &
                 .700_r8*rxt(k,696)*y(k,301) +rxt(k,704)*y(k,302) + &
                 .830_r8*rxt(k,712)*y(k,303) +rxt(k,722)*y(k,304) + &
                 .700_r8*rxt(k,731)*y(k,305) +.910_r8*rxt(k,741)*y(k,306) + &
                 rxt(k,752)*y(k,307) +.700_r8*rxt(k,761)*y(k,308) + &
                 .700_r8*rxt(k,766)*y(k,309) +.700_r8*rxt(k,773)*y(k,310) + &
                 .700_r8*rxt(k,777)*y(k,311) +.700_r8*rxt(k,781)*y(k,312) + &
                 .700_r8*rxt(k,785)*y(k,313))*y(k,147) + (rxt(k,18) + &
                 rxt(k,202)*y(k,258) +rxt(k,204)*y(k,157) +rxt(k,205)*y(k,295) + &
                 rxt(k,374)*y(k,16) +rxt(k,393)*y(k,286) + &
                 1.460_r8*rxt(k,478)*y(k,274) +2.000_r8*rxt(k,590)*y(k,237) + &
                 rxt(k,598)*y(k,238) +2.000_r8*rxt(k,610)*y(k,240) + &
                 rxt(k,618)*y(k,241) +1.500_r8*rxt(k,630)*y(k,245) + &
                 rxt(k,638)*y(k,246) +1.460_r8*rxt(k,650)*y(k,281) + &
                 rxt(k,658)*y(k,282) +1.950_r8*rxt(k,670)*y(k,290) + &
                 rxt(k,678)*y(k,291) +rxt(k,697)*y(k,301) +rxt(k,705)*y(k,302) + &
                 rxt(k,713)*y(k,303) +rxt(k,723)*y(k,304) +rxt(k,732)*y(k,305) + &
                 rxt(k,742)*y(k,306) +rxt(k,753)*y(k,307) +rxt(k,759)*y(k,212) + &
                 .500_r8*rxt(k,764)*y(k,213))*y(k,149) + (rxt(k,193)*y(k,94) + &
                 rxt(k,348)*y(k,150) +rxt(k,364)*y(k,151) + &
                 .500_r8*rxt(k,380)*y(k,127) +rxt(k,408)*y(k,1) + &
                 .400_r8*rxt(k,422)*y(k,107) +.190_r8*rxt(k,423)*y(k,108) + &
                 rxt(k,452)*y(k,111) +.500_r8*rxt(k,454)*y(k,112) + &
                 .080_r8*rxt(k,460)*y(k,115) +.150_r8*rxt(k,463)*y(k,116) + &
                 .130_r8*rxt(k,466)*y(k,117) +.040_r8*rxt(k,470)*y(k,118) + &
                 .070_r8*rxt(k,485)*y(k,121) +.040_r8*rxt(k,502)*y(k,139) + &
                 rxt(k,718)*y(k,204) +rxt(k,737)*y(k,206) +rxt(k,757)*y(k,210) + &
                 rxt(k,769)*y(k,214) +rxt(k,783)*y(k,221) +rxt(k,787)*y(k,223)) &
                 *y(k,295) + (1.640_r8*rxt(k,585)*y(k,237) +rxt(k,586)*y(k,252) + &
                 .820_r8*rxt(k,587)*y(k,253) +.700_r8*rxt(k,588)*y(k,258) + &
                 rxt(k,591)*y(k,302) +rxt(k,592)*y(k,304) +rxt(k,593)*y(k,307)) &
                 *y(k,237) + (1.640_r8*rxt(k,605)*y(k,240) +rxt(k,606)*y(k,252) + &
                 .820_r8*rxt(k,607)*y(k,253) +.500_r8*rxt(k,608)*y(k,258) + &
                 rxt(k,611)*y(k,302) +rxt(k,612)*y(k,304) +rxt(k,613)*y(k,307)) &
                 *y(k,240) + (.940_r8*rxt(k,625)*y(k,245) + &
                 .500_r8*rxt(k,626)*y(k,252) +.360_r8*rxt(k,627)*y(k,253) + &
                 .240_r8*rxt(k,628)*y(k,258) +.500_r8*rxt(k,631)*y(k,302) + &
                 .500_r8*rxt(k,632)*y(k,304) +.500_r8*rxt(k,633)*y(k,307))*y(k,245) &
                  + (.460_r8*rxt(k,645)*y(k,252) +.310_r8*rxt(k,646)*y(k,253) + &
                 .230_r8*rxt(k,647)*y(k,258) +.860_r8*rxt(k,648)*y(k,281) + &
                 .460_r8*rxt(k,651)*y(k,302) +.460_r8*rxt(k,652)*y(k,304) + &
                 .460_r8*rxt(k,653)*y(k,307))*y(k,281) &
                  + (.950_r8*rxt(k,665)*y(k,252) +.770_r8*rxt(k,666)*y(k,253) + &
                 .480_r8*rxt(k,667)*y(k,258) +1.540_r8*rxt(k,668)*y(k,290) + &
                 .950_r8*rxt(k,671)*y(k,302) +.950_r8*rxt(k,672)*y(k,304) + &
                 .950_r8*rxt(k,673)*y(k,307))*y(k,290) &
                  + (.170_r8*rxt(k,459)*y(k,115) +.170_r8*rxt(k,469)*y(k,118) + &
                 .170_r8*rxt(k,484)*y(k,121) +.170_r8*rxt(k,501)*y(k,139))*y(k,158) &
                  + (.460_r8*rxt(k,474)*y(k,252) +.070_r8*rxt(k,475)*y(k,253) + &
                 .240_r8*rxt(k,476)*y(k,258) +.160_r8*rxt(k,477)*y(k,274))*y(k,274) &
                  + (rxt(k,11) +rxt(k,216))*y(k,94) + (rxt(k,78) +rxt(k,407))*y(k,131) &
                  + (rxt(k,13) +rxt(k,217))*y(k,138) + (.600_r8*rxt(k,86) +rxt(k,352)) &
                 *y(k,163) + (rxt(k,95) +rxt(k,790))*y(k,204) + (rxt(k,97) + &
                 rxt(k,791))*y(k,206) + (rxt(k,101) +rxt(k,792))*y(k,210) +rxt(k,19) &
                 *y(k,1) +rxt(k,120)*y(k,21) +rxt(k,139)*y(k,61) +rxt(k,9)*y(k,93) &
                  +rxt(k,47)*y(k,97) +rxt(k,57)*y(k,107) +rxt(k,58)*y(k,108) &
                  +2.000_r8*rxt(k,59)*y(k,110) +2.000_r8*rxt(k,60)*y(k,111) +rxt(k,61) &
                 *y(k,112) +rxt(k,62)*y(k,113) +rxt(k,64)*y(k,115) +rxt(k,65)*y(k,116) &
                  +rxt(k,66)*y(k,117) +rxt(k,67)*y(k,118) +rxt(k,68)*y(k,119) &
                  +rxt(k,69)*y(k,120) +.750_r8*rxt(k,74)*y(k,127) +.750_r8*rxt(k,80) &
                 *y(k,133) +rxt(k,82)*y(k,139) +rxt(k,83)*y(k,150) +rxt(k,84)*y(k,151) &
                  +rxt(k,85)*y(k,162) +rxt(k,575)*y(k,164) +rxt(k,103)*y(k,214) &
                  +.500_r8*rxt(k,105)*y(k,217) +.460_r8*rxt(k,106)*y(k,218) &
                  +rxt(k,107)*y(k,219) +.460_r8*rxt(k,108)*y(k,220) +rxt(k,109) &
                 *y(k,221) +rxt(k,110)*y(k,222) +rxt(k,111)*y(k,223) +rxt(k,112) &
                 *y(k,224) +.460_r8*rxt(k,461)*y(k,270)*y(k,258)
         loss(k,311) = (rxt(k,584)* y(k,4) +rxt(k,604)* y(k,7) +rxt(k,374)* y(k,16) &
                  +rxt(k,624)* y(k,17) +rxt(k,353)* y(k,30) +rxt(k,299)* y(k,43) &
                  +rxt(k,332)* y(k,46) +rxt(k,360)* y(k,50) +rxt(k,793)* y(k,69) &
                  +rxt(k,473)* y(k,109) +rxt(k,644)* y(k,125) +rxt(k,664)* y(k,135) &
                  +rxt(k,203)* y(k,147) +rxt(k,213)* y(k,148) +rxt(k,204)* y(k,157) &
                  +rxt(k,755)* y(k,202) +rxt(k,708)* y(k,203) +rxt(k,727)* y(k,205) &
                  +rxt(k,759)* y(k,212) +rxt(k,764)* y(k,213) +rxt(k,590)* y(k,237) &
                  +rxt(k,598)* y(k,238) +rxt(k,610)* y(k,240) +rxt(k,618)* y(k,241) &
                  +rxt(k,630)* y(k,245) +rxt(k,638)* y(k,246) +rxt(k,202)* y(k,258) &
                  +rxt(k,478)* y(k,274) +rxt(k,650)* y(k,281) +rxt(k,658)* y(k,282) &
                  +rxt(k,393)* y(k,286) +rxt(k,670)* y(k,290) +rxt(k,678)* y(k,291) &
                  +rxt(k,205)* y(k,295) +rxt(k,697)* y(k,301) +rxt(k,705)* y(k,302) &
                  +rxt(k,713)* y(k,303) +rxt(k,723)* y(k,304) +rxt(k,732)* y(k,305) &
                  +rxt(k,742)* y(k,306) +rxt(k,753)* y(k,307) + rxt(k,17) + rxt(k,18) &
                  + rxt(k,831) + het_rates(k,149))* y(k,149)
         prod(k,311) = (rxt(k,138) +rxt(k,231)*y(k,57) +rxt(k,233)*y(k,157) + &
                 rxt(k,234)*y(k,295))*y(k,61) + (rxt(k,13) +rxt(k,14) +rxt(k,217)) &
                 *y(k,138) + (rxt(k,215)*y(k,93) +rxt(k,349)*y(k,163) + &
                 rxt(k,398)*y(k,131))*y(k,295) + (rxt(k,121) +rxt(k,264)*y(k,157)) &
                 *y(k,21) + (rxt(k,200)*y(k,158) +rxt(k,201)*y(k,157))*y(k,148) &
                  +rxt(k,278)*y(k,93)*y(k,75) +rxt(k,10)*y(k,94) +.400_r8*rxt(k,86) &
                 *y(k,163)
         loss(k,232) = (rxt(k,348)* y(k,295) + rxt(k,83) + het_rates(k,150))* y(k,150)
         prod(k,232) = (.870_r8*rxt(k,519)*y(k,272) +.330_r8*rxt(k,521)*y(k,273) + &
                 .070_r8*rxt(k,525)*y(k,275) +.150_r8*rxt(k,527)*y(k,276) + &
                 .120_r8*rxt(k,539)*y(k,293))*y(k,147) &
                  + (.440_r8*rxt(k,467)*y(k,272) +.150_r8*rxt(k,471)*y(k,273) + &
                 .060_r8*rxt(k,479)*y(k,275) +.120_r8*rxt(k,482)*y(k,276) + &
                 .100_r8*rxt(k,499)*y(k,293))*y(k,258) &
                  + (.830_r8*rxt(k,469)*y(k,118) +.130_r8*rxt(k,484)*y(k,121) + &
                 .220_r8*rxt(k,501)*y(k,139))*y(k,158) +.250_r8*rxt(k,80)*y(k,133) &
                  +.100_r8*rxt(k,502)*y(k,295)*y(k,139)
         loss(k,236) = (rxt(k,364)* y(k,295) + rxt(k,84) + het_rates(k,151))* y(k,151)
         prod(k,236) = (.940_r8*rxt(k,513)*y(k,269) +.340_r8*rxt(k,521)*y(k,273) + &
                 .400_r8*rxt(k,525)*y(k,275) +.810_r8*rxt(k,527)*y(k,276) + &
                 .130_r8*rxt(k,539)*y(k,293))*y(k,147) &
                  + (.550_r8*rxt(k,457)*y(k,269) +.150_r8*rxt(k,471)*y(k,273) + &
                 .280_r8*rxt(k,479)*y(k,275) +.680_r8*rxt(k,482)*y(k,276) + &
                 .100_r8*rxt(k,499)*y(k,293))*y(k,258) &
                  + (.500_r8*rxt(k,380)*y(k,127) +.500_r8*rxt(k,399)*y(k,133) + &
                 .350_r8*rxt(k,423)*y(k,108) +.350_r8*rxt(k,502)*y(k,139))*y(k,295) &
                  + (.830_r8*rxt(k,459)*y(k,115) +.700_r8*rxt(k,484)*y(k,121) + &
                 .610_r8*rxt(k,501)*y(k,139))*y(k,158) +rxt(k,353)*y(k,149)*y(k,30) &
                  +.250_r8*rxt(k,74)*y(k,127)
         loss(k,16) = ( + het_rates(k,152))* y(k,152)
         prod(k,16) = 0._r8
         loss(k,17) = ( + het_rates(k,153))* y(k,153)
         prod(k,17) = 0._r8
         loss(k,18) = ( + het_rates(k,154))* y(k,154)
         prod(k,18) = 0._r8
         loss(k,19) = ( + het_rates(k,155))* y(k,155)
         prod(k,19) = 0._r8
         loss(k,20) = ( + het_rates(k,156))* y(k,156)
         prod(k,20) = 0._r8
         loss(k,305) = (rxt(k,265)* y(k,20) +rxt(k,264)* y(k,21) +rxt(k,300)* y(k,43) &
                  +rxt(k,235)* y(k,60) +rxt(k,233)* y(k,61) +rxt(k,176)* y(k,79) &
                  +rxt(k,177)* y(k,81) +rxt(k,267)* y(k,83) +rxt(k,238)* y(k,87) &
                  +rxt(k,269)* y(k,95) +rxt(k,241)* y(k,96) +rxt(k,209)* y(k,147) &
                  + (rxt(k,199) +rxt(k,201))* y(k,148) +rxt(k,204)* y(k,149) &
                  + 2._r8*rxt(k,174)* y(k,157) +rxt(k,173)* y(k,158) +rxt(k,795) &
                 * y(k,161) +rxt(k,182)* y(k,258) +rxt(k,188)* y(k,295) + rxt(k,175) &
                  + het_rates(k,157))* y(k,157)
         prod(k,305) = (rxt(k,198) +rxt(k,194)*y(k,147) +rxt(k,195)*y(k,148))*y(k,136) &
                  + (rxt(k,155) +rxt(k,804))*y(k,174) + (rxt(k,170) +rxt(k,171)) &
                 *y(k,294) +rxt(k,119)*y(k,20) +.180_r8*rxt(k,39)*y(k,55) +rxt(k,137) &
                 *y(k,60) +rxt(k,40)*y(k,64) +rxt(k,180)*y(k,258)*y(k,78) +rxt(k,14) &
                 *y(k,138) +rxt(k,15)*y(k,147) +rxt(k,16)*y(k,148) +rxt(k,18)*y(k,149) &
                  +rxt(k,8)*y(k,158) +rxt(k,151)*y(k,160) +rxt(k,797)*y(k,172) &
                  +rxt(k,156)*y(k,175) +rxt(k,157)*y(k,176) +rxt(k,190)*y(k,295) &
                 *y(k,295) +rxt(k,3)*y(k,319)
         loss(k,313) = (rxt(k,602)* y(k,4) +rxt(k,622)* y(k,7) +rxt(k,642)* y(k,17) &
                  +rxt(k,256)* y(k,18) +rxt(k,323)* y(k,26) +rxt(k,354)* y(k,30) &
                  +rxt(k,224)* y(k,57) +rxt(k,184)* y(k,78) +rxt(k,486)* y(k,109) &
                  +rxt(k,459)* y(k,115) +rxt(k,469)* y(k,118) +rxt(k,484)* y(k,121) &
                  +rxt(k,662)* y(k,125) +rxt(k,385)* y(k,126) +rxt(k,403)* y(k,132) &
                  +rxt(k,682)* y(k,135) +rxt(k,501)* y(k,139) +rxt(k,208)* y(k,147) &
                  +rxt(k,200)* y(k,148) +rxt(k,173)* y(k,157) +rxt(k,568)* y(k,165) &
                  +rxt(k,799)* y(k,172) +rxt(k,805)* y(k,174) +rxt(k,762)* y(k,212) &
                  +rxt(k,767)* y(k,213) +rxt(k,183)* y(k,258) +rxt(k,172)* y(k,294) &
                  +rxt(k,189)* y(k,295) + rxt(k,7) + rxt(k,8) + het_rates(k,158)) &
                 * y(k,158)
         prod(k,313) = (.150_r8*rxt(k,337)*y(k,252) +.150_r8*rxt(k,390)*y(k,286) + &
                 .150_r8*rxt(k,703)*y(k,302) +.150_r8*rxt(k,721)*y(k,304) + &
                 .150_r8*rxt(k,751)*y(k,307))*y(k,258) +rxt(k,175)*y(k,157)
         loss(k,21) = ( + het_rates(k,159))* y(k,159)
         prod(k,21) = 0._r8
         loss(k,111) = (rxt(k,806)* y(k,174) + rxt(k,151) + het_rates(k,160)) &
                 * y(k,160)
         prod(k,111) = (rxt(k,228)*y(k,60) +rxt(k,258)*y(k,20))*y(k,60)
         loss(k,119) = (rxt(k,795)* y(k,157) +rxt(k,796)* y(k,295) + rxt(k,154) &
                  + het_rates(k,161))* y(k,161)
         prod(k,119) = 0._r8
         loss(k,82) = ( + rxt(k,85) + rxt(k,832) + het_rates(k,162))* y(k,162)
         prod(k,82) = (rxt(k,379)*y(k,97) +.500_r8*rxt(k,399)*y(k,133))*y(k,295)
         loss(k,147) = (rxt(k,349)* y(k,295) + rxt(k,86) + rxt(k,352) &
                  + het_rates(k,163))* y(k,163)
         prod(k,147) =rxt(k,351)*y(k,252)*y(k,148)
         loss(k,70) = ( + rxt(k,575) + het_rates(k,164))* y(k,164)
         prod(k,70) =rxt(k,570)*y(k,235)*y(k,148)
         loss(k,135) = (rxt(k,567)* y(k,148) +rxt(k,568)* y(k,158) + het_rates(k,165)) &
                 * y(k,165)
         prod(k,135) = (.070_r8*rxt(k,554)*y(k,67) +.060_r8*rxt(k,566)*y(k,166) + &
                 .070_r8*rxt(k,582)*y(k,231))*y(k,295) +rxt(k,30)*y(k,33) &
                  +rxt(k,552)*y(k,251)*y(k,147)
         loss(k,78) = (rxt(k,566)* y(k,295) + het_rates(k,166))* y(k,166)
         prod(k,78) =.530_r8*rxt(k,543)*y(k,295)*y(k,8)
         loss(k,112) = (rxt(k,569)* y(k,295) + rxt(k,87) + het_rates(k,167))* y(k,167)
         prod(k,112) =rxt(k,564)*y(k,296)*y(k,258)
         loss(k,22) = ( + het_rates(k,168))* y(k,168)
         prod(k,22) = 0._r8
         loss(k,23) = ( + het_rates(k,169))* y(k,169)
         prod(k,23) = 0._r8
         loss(k,148) = (rxt(k,367)* y(k,295) + rxt(k,88) + het_rates(k,170))* y(k,170)
         prod(k,148) =rxt(k,365)*y(k,297)*y(k,258)
         loss(k,122) = (rxt(k,371)* y(k,295) + rxt(k,89) + het_rates(k,171))* y(k,171)
         prod(k,122) =.850_r8*rxt(k,369)*y(k,298)*y(k,258)
         loss(k,143) = (rxt(k,799)* y(k,158) +rxt(k,802)* y(k,295) + rxt(k,797) &
                  + het_rates(k,172))* y(k,172)
         prod(k,143) =rxt(k,154)*y(k,161) +rxt(k,155)*y(k,174)
         loss(k,24) = ( + rxt(k,152) + het_rates(k,173))* y(k,173)
         prod(k,24) = 0._r8
         loss(k,246) = (rxt(k,800)* y(k,20) +rxt(k,801)* y(k,60) +rxt(k,803)* y(k,148) &
                  +rxt(k,805)* y(k,158) +rxt(k,806)* y(k,160) +rxt(k,807)* y(k,295) &
                  + rxt(k,155) + rxt(k,804) + het_rates(k,174))* y(k,174)
         prod(k,246) = (rxt(k,797) +rxt(k,799)*y(k,158) +rxt(k,802)*y(k,295))*y(k,172) &
                  +rxt(k,795)*y(k,161)*y(k,157) +rxt(k,156)*y(k,175)
         loss(k,216) = (rxt(k,798)* y(k,295) + rxt(k,156) + het_rates(k,175)) &
                 * y(k,175)
         prod(k,216) = (rxt(k,804) +rxt(k,800)*y(k,20) +rxt(k,801)*y(k,60) + &
                 rxt(k,803)*y(k,148) +rxt(k,805)*y(k,158) +rxt(k,806)*y(k,160) + &
                 rxt(k,807)*y(k,295))*y(k,174) + (rxt(k,793)*y(k,149) + &
                 rxt(k,794)*y(k,295) +.500_r8*rxt(k,808)*y(k,295))*y(k,69) &
                  +rxt(k,796)*y(k,295)*y(k,161) +rxt(k,157)*y(k,176)
         loss(k,98) = (rxt(k,809)* y(k,319) + rxt(k,157) + het_rates(k,176))* y(k,176)
         prod(k,98) =rxt(k,153)*y(k,82) +rxt(k,798)*y(k,295)*y(k,175)
         loss(k,25) = ( + het_rates(k,177))* y(k,177)
         prod(k,25) = 0._r8
         loss(k,26) = ( + het_rates(k,178))* y(k,178)
         prod(k,26) = 0._r8
         loss(k,27) = ( + het_rates(k,179))* y(k,179)
         prod(k,27) = 0._r8
         loss(k,28) = ( + het_rates(k,180))* y(k,180)
         prod(k,28) = 0._r8
         loss(k,29) = ( + rxt(k,158) + het_rates(k,181))* y(k,181)
         prod(k,29) = 0._r8
         loss(k,30) = ( + rxt(k,159) + het_rates(k,182))* y(k,182)
         prod(k,30) = 0._r8
         loss(k,31) = ( + rxt(k,160) + het_rates(k,183))* y(k,183)
         prod(k,31) = 0._r8
         loss(k,32) = ( + rxt(k,161) + het_rates(k,184))* y(k,184)
         prod(k,32) = 0._r8
         loss(k,33) = ( + rxt(k,162) + het_rates(k,185))* y(k,185)
         prod(k,33) = 0._r8
         loss(k,34) = ( + rxt(k,163) + het_rates(k,186))* y(k,186)
         prod(k,34) = 0._r8
         loss(k,35) = ( + rxt(k,164) + het_rates(k,187))* y(k,187)
         prod(k,35) = 0._r8
         loss(k,36) = ( + rxt(k,165) + het_rates(k,188))* y(k,188)
         prod(k,36) = 0._r8
         loss(k,37) = ( + rxt(k,166) + het_rates(k,189))* y(k,189)
         prod(k,37) = 0._r8
         loss(k,38) = ( + rxt(k,167) + het_rates(k,190))* y(k,190)
         prod(k,38) = 0._r8
         loss(k,39) = ( + het_rates(k,191))* y(k,191)
         prod(k,39) = (.0245005_r8*rxt(k,843)*y(k,239) + &
                 .1279005_r8*rxt(k,848)*y(k,242) +.0097005_r8*rxt(k,853)*y(k,244) + &
                 .0245005_r8*rxt(k,856)*y(k,247) +.0003005_r8*rxt(k,861)*y(k,277) + &
                 .1056005_r8*rxt(k,865)*y(k,280) +.0245005_r8*rxt(k,869)*y(k,283) + &
                 .0245005_r8*rxt(k,874)*y(k,292) +.0154005_r8*rxt(k,880)*y(k,315) + &
                 .0063005_r8*rxt(k,883)*y(k,317))*y(k,147) &
                  + (.0508005_r8*rxt(k,842)*y(k,239) + &
                 .2202005_r8*rxt(k,847)*y(k,242) +.0023005_r8*rxt(k,852)*y(k,244) + &
                 .0508005_r8*rxt(k,855)*y(k,247) +.0031005_r8*rxt(k,860)*y(k,277) + &
                 .2381005_r8*rxt(k,864)*y(k,280) +.0508005_r8*rxt(k,868)*y(k,283) + &
                 .0508005_r8*rxt(k,873)*y(k,292) +.1364005_r8*rxt(k,879)*y(k,315) + &
                 .1677005_r8*rxt(k,882)*y(k,317))*y(k,258) &
                  + (.0508005_r8*rxt(k,844)*y(k,4) +.2202005_r8*rxt(k,849)*y(k,7) + &
                 .0508005_r8*rxt(k,857)*y(k,17) +.0508005_r8*rxt(k,870)*y(k,125) + &
                 .0508005_r8*rxt(k,875)*y(k,135))*y(k,158) +rxt(k,811)*y(k,77) &
                  +.5931005_r8*rxt(k,877)*y(k,295)*y(k,198)
         loss(k,40) = ( + het_rates(k,192))* y(k,192)
         prod(k,40) = (.0082005_r8*rxt(k,843)*y(k,239) + &
                 .1792005_r8*rxt(k,848)*y(k,242) +.0034005_r8*rxt(k,853)*y(k,244) + &
                 .0082005_r8*rxt(k,856)*y(k,247) +.0003005_r8*rxt(k,861)*y(k,277) + &
                 .1026005_r8*rxt(k,865)*y(k,280) +.0082005_r8*rxt(k,869)*y(k,283) + &
                 .0082005_r8*rxt(k,874)*y(k,292) +.0452005_r8*rxt(k,880)*y(k,315) + &
                 .0237005_r8*rxt(k,883)*y(k,317))*y(k,147) &
                  + (.1149005_r8*rxt(k,842)*y(k,239) + &
                 .2067005_r8*rxt(k,847)*y(k,242) +.0008005_r8*rxt(k,852)*y(k,244) + &
                 .1149005_r8*rxt(k,855)*y(k,247) +.0035005_r8*rxt(k,860)*y(k,277) + &
                 .1308005_r8*rxt(k,864)*y(k,280) +.1149005_r8*rxt(k,868)*y(k,283) + &
                 .1149005_r8*rxt(k,873)*y(k,292) +.0101005_r8*rxt(k,879)*y(k,315) + &
                 .0174005_r8*rxt(k,882)*y(k,317))*y(k,258) &
                  + (.1149005_r8*rxt(k,844)*y(k,4) +.2067005_r8*rxt(k,849)*y(k,7) + &
                 .1149005_r8*rxt(k,857)*y(k,17) +.1149005_r8*rxt(k,870)*y(k,125) + &
                 .1149005_r8*rxt(k,875)*y(k,135))*y(k,158) &
                  +.1534005_r8*rxt(k,877)*y(k,295)*y(k,198)
         loss(k,41) = ( + het_rates(k,193))* y(k,193)
         prod(k,41) = (.0772005_r8*rxt(k,843)*y(k,239) + &
                 .0676005_r8*rxt(k,848)*y(k,242) +.1579005_r8*rxt(k,853)*y(k,244) + &
                 .0772005_r8*rxt(k,856)*y(k,247) +.0073005_r8*rxt(k,861)*y(k,277) + &
                 .0521005_r8*rxt(k,865)*y(k,280) +.0772005_r8*rxt(k,869)*y(k,283) + &
                 .0772005_r8*rxt(k,874)*y(k,292) +.0966005_r8*rxt(k,880)*y(k,315) + &
                 .0025005_r8*rxt(k,883)*y(k,317))*y(k,147) &
                  + (.0348005_r8*rxt(k,842)*y(k,239) + &
                 .0653005_r8*rxt(k,847)*y(k,242) +.0843005_r8*rxt(k,852)*y(k,244) + &
                 .0348005_r8*rxt(k,855)*y(k,247) +.0003005_r8*rxt(k,860)*y(k,277) + &
                 .0348005_r8*rxt(k,864)*y(k,280) +.0348005_r8*rxt(k,868)*y(k,283) + &
                 .0348005_r8*rxt(k,873)*y(k,292) +.0763005_r8*rxt(k,879)*y(k,315) + &
                 .086_r8*rxt(k,882)*y(k,317))*y(k,258) &
                  + (.0348005_r8*rxt(k,844)*y(k,4) +.0653005_r8*rxt(k,849)*y(k,7) + &
                 .0348005_r8*rxt(k,857)*y(k,17) +.0348005_r8*rxt(k,870)*y(k,125) + &
                 .0348005_r8*rxt(k,875)*y(k,135))*y(k,158) &
                  +.0459005_r8*rxt(k,877)*y(k,295)*y(k,198)
         loss(k,42) = ( + het_rates(k,194))* y(k,194)
         prod(k,42) = (.0332005_r8*rxt(k,843)*y(k,239) +.079_r8*rxt(k,848)*y(k,242) + &
                 .0059005_r8*rxt(k,853)*y(k,244) +.0332005_r8*rxt(k,856)*y(k,247) + &
                 .0057005_r8*rxt(k,861)*y(k,277) +.0143005_r8*rxt(k,865)*y(k,280) + &
                 .0332005_r8*rxt(k,869)*y(k,283) +.0332005_r8*rxt(k,874)*y(k,292) + &
                 .0073005_r8*rxt(k,880)*y(k,315) +.011_r8*rxt(k,883)*y(k,317)) &
                 *y(k,147) + (.0554005_r8*rxt(k,842)*y(k,239) + &
                 .1284005_r8*rxt(k,847)*y(k,242) +.0443005_r8*rxt(k,852)*y(k,244) + &
                 .0554005_r8*rxt(k,855)*y(k,247) +.0271005_r8*rxt(k,860)*y(k,277) + &
                 .0076005_r8*rxt(k,864)*y(k,280) +.0554005_r8*rxt(k,868)*y(k,283) + &
                 .0554005_r8*rxt(k,873)*y(k,292) +.2157005_r8*rxt(k,879)*y(k,315) + &
                 .0512005_r8*rxt(k,882)*y(k,317))*y(k,258) &
                  + (.1749305_r8*rxt(k,841)*y(k,4) +.1749305_r8*rxt(k,846)*y(k,7) + &
                 .1749305_r8*rxt(k,854)*y(k,17) +.0590245_r8*rxt(k,859)*y(k,109) + &
                 .1749305_r8*rxt(k,867)*y(k,125) +.1749305_r8*rxt(k,872)*y(k,135)) &
                 *y(k,149) + (.0554005_r8*rxt(k,844)*y(k,4) + &
                 .1284005_r8*rxt(k,849)*y(k,7) +.0554005_r8*rxt(k,857)*y(k,17) + &
                 .0033005_r8*rxt(k,862)*y(k,109) +.0554005_r8*rxt(k,870)*y(k,125) + &
                 .0554005_r8*rxt(k,875)*y(k,135))*y(k,158) &
                  +.0085005_r8*rxt(k,877)*y(k,295)*y(k,198)
         loss(k,43) = ( + het_rates(k,195))* y(k,195)
         prod(k,43) = (.130_r8*rxt(k,843)*y(k,239) +.1254005_r8*rxt(k,848)*y(k,242) + &
                 .0536005_r8*rxt(k,853)*y(k,244) +.130_r8*rxt(k,856)*y(k,247) + &
                 .0623005_r8*rxt(k,861)*y(k,277) +.0166005_r8*rxt(k,865)*y(k,280) + &
                 .130_r8*rxt(k,869)*y(k,283) +.130_r8*rxt(k,874)*y(k,292) + &
                 .238_r8*rxt(k,880)*y(k,315) +.1185005_r8*rxt(k,883)*y(k,317)) &
                 *y(k,147) + (.1278005_r8*rxt(k,842)*y(k,239) + &
                 .114_r8*rxt(k,847)*y(k,242) +.1621005_r8*rxt(k,852)*y(k,244) + &
                 .1278005_r8*rxt(k,855)*y(k,247) +.0474005_r8*rxt(k,860)*y(k,277) + &
                 .0113005_r8*rxt(k,864)*y(k,280) +.1278005_r8*rxt(k,868)*y(k,283) + &
                 .1278005_r8*rxt(k,873)*y(k,292) +.0738005_r8*rxt(k,879)*y(k,315) + &
                 .1598005_r8*rxt(k,882)*y(k,317))*y(k,258) &
                  + (.5901905_r8*rxt(k,841)*y(k,4) +.5901905_r8*rxt(k,846)*y(k,7) + &
                 .5901905_r8*rxt(k,854)*y(k,17) +.0250245_r8*rxt(k,859)*y(k,109) + &
                 .5901905_r8*rxt(k,867)*y(k,125) +.5901905_r8*rxt(k,872)*y(k,135)) &
                 *y(k,149) + (.1278005_r8*rxt(k,844)*y(k,4) + &
                 .114_r8*rxt(k,849)*y(k,7) +.1278005_r8*rxt(k,857)*y(k,17) + &
                 .1278005_r8*rxt(k,870)*y(k,125) +.1278005_r8*rxt(k,875)*y(k,135)) &
                 *y(k,158) +.0128005_r8*rxt(k,877)*y(k,295)*y(k,198)
         loss(k,44) = ( + rxt(k,833) + het_rates(k,196))* y(k,196)
         prod(k,44) = (.360_r8*rxt(k,605)*y(k,240) +.180_r8*rxt(k,607)*y(k,253) + &
                 .500_r8*rxt(k,608)*y(k,258) +.070_r8*rxt(k,609)*y(k,147))*y(k,240) &
                  +.300_r8*rxt(k,617)*y(k,241)*y(k,147)
         loss(k,45) = ( + rxt(k,904) + het_rates(k,197))* y(k,197)
         prod(k,45) = 0._r8
         loss(k,46) = (rxt(k,877)* y(k,295) + het_rates(k,198))* y(k,198)
         prod(k,46) = 0._r8
         loss(k,85) = ( + rxt(k,90) + het_rates(k,199))* y(k,199)
         prod(k,85) = (.100_r8*rxt(k,574)*y(k,228) +.230_r8*rxt(k,576)*y(k,229)) &
                 *y(k,295)
         loss(k,282) = (rxt(k,689)* y(k,295) + rxt(k,91) + het_rates(k,200))* y(k,200)
         prod(k,282) = (.140_r8*rxt(k,594)*y(k,252) +.130_r8*rxt(k,595)*y(k,253) + &
                 .250_r8*rxt(k,596)*y(k,258) +.110_r8*rxt(k,597)*y(k,147) + &
                 .140_r8*rxt(k,598)*y(k,149) +.140_r8*rxt(k,599)*y(k,302) + &
                 .140_r8*rxt(k,600)*y(k,304) +.140_r8*rxt(k,601)*y(k,307))*y(k,238) &
                  + (.680_r8*rxt(k,636)*y(k,246) +.900_r8*rxt(k,656)*y(k,282) + &
                 .180_r8*rxt(k,691)*y(k,300) +.900_r8*rxt(k,765)*y(k,309))*y(k,258) &
                  +.700_r8*rxt(k,692)*y(k,300)*y(k,147)
         loss(k,133) = (rxt(k,690)* y(k,295) + rxt(k,92) + het_rates(k,201))* y(k,201)
         prod(k,133) = (.900_r8*rxt(k,616)*y(k,241) +.900_r8*rxt(k,676)*y(k,291)) &
                 *y(k,258)
         loss(k,284) = (rxt(k,755)* y(k,149) +rxt(k,756)* y(k,295) + rxt(k,93) &
                  + het_rates(k,202))* y(k,202)
         prod(k,284) = (1.640_r8*rxt(k,585)*y(k,237) +rxt(k,586)*y(k,252) + &
                 .820_r8*rxt(k,587)*y(k,253) +.700_r8*rxt(k,588)*y(k,258) + &
                 .930_r8*rxt(k,589)*y(k,147) +rxt(k,590)*y(k,149) + &
                 rxt(k,591)*y(k,302) +rxt(k,592)*y(k,304) +rxt(k,593)*y(k,307)) &
                 *y(k,237) + (.390_r8*rxt(k,594)*y(k,252) + &
                 .420_r8*rxt(k,595)*y(k,253) +.290_r8*rxt(k,596)*y(k,258) + &
                 .300_r8*rxt(k,597)*y(k,147) +.390_r8*rxt(k,598)*y(k,149) + &
                 .390_r8*rxt(k,599)*y(k,302) +.390_r8*rxt(k,600)*y(k,304) + &
                 .390_r8*rxt(k,601)*y(k,307))*y(k,238) + (rxt(k,783)*y(k,221) + &
                 rxt(k,787)*y(k,223) +rxt(k,789)*y(k,225))*y(k,295) &
                  +.220_r8*rxt(k,602)*y(k,158)*y(k,4) +.500_r8*rxt(k,105)*y(k,217) &
                  +rxt(k,107)*y(k,219) +rxt(k,109)*y(k,221) +rxt(k,111)*y(k,223) &
                  +rxt(k,113)*y(k,225)
         loss(k,220) = (rxt(k,708)* y(k,149) +rxt(k,717)* y(k,295) + rxt(k,94) &
                  + het_rates(k,203))* y(k,203)
         prod(k,220) =.170_r8*rxt(k,602)*y(k,158)*y(k,4) +rxt(k,757)*y(k,295)*y(k,210) &
                  +.500_r8*rxt(k,694)*y(k,301)*y(k,253)
         loss(k,157) = (rxt(k,718)* y(k,295) + rxt(k,95) + rxt(k,790) &
                  + het_rates(k,204))* y(k,204)
         prod(k,157) =rxt(k,684)*y(k,302)*y(k,148)
         loss(k,255) = (rxt(k,727)* y(k,149) +rxt(k,736)* y(k,295) + rxt(k,96) &
                  + het_rates(k,205))* y(k,205)
         prod(k,255) = (.900_r8*rxt(k,625)*y(k,245) +.480_r8*rxt(k,626)*y(k,252) + &
                 .340_r8*rxt(k,627)*y(k,253) +.220_r8*rxt(k,628)*y(k,258) + &
                 .440_r8*rxt(k,629)*y(k,147) +.480_r8*rxt(k,630)*y(k,149) + &
                 .480_r8*rxt(k,631)*y(k,302) +.480_r8*rxt(k,632)*y(k,304) + &
                 .480_r8*rxt(k,633)*y(k,307))*y(k,245) &
                  + (.350_r8*rxt(k,594)*y(k,252) +.200_r8*rxt(k,595)*y(k,253) + &
                 .270_r8*rxt(k,597)*y(k,147) +.350_r8*rxt(k,598)*y(k,149) + &
                 .350_r8*rxt(k,599)*y(k,302) +.350_r8*rxt(k,600)*y(k,304) + &
                 .350_r8*rxt(k,601)*y(k,307))*y(k,238) &
                  + (.410_r8*rxt(k,634)*y(k,252) +.310_r8*rxt(k,635)*y(k,253) + &
                 .310_r8*rxt(k,637)*y(k,147) +.410_r8*rxt(k,638)*y(k,149) + &
                 .410_r8*rxt(k,639)*y(k,302) +.410_r8*rxt(k,640)*y(k,304) + &
                 .410_r8*rxt(k,641)*y(k,307))*y(k,246) + (rxt(k,759)*y(k,149) + &
                 rxt(k,762)*y(k,158))*y(k,212) + (rxt(k,114) +rxt(k,788)*y(k,295)) &
                 *y(k,226) + (.100_r8*rxt(k,760)*y(k,258) + &
                 .700_r8*rxt(k,761)*y(k,147))*y(k,308)
         loss(k,158) = (rxt(k,737)* y(k,295) + rxt(k,97) + rxt(k,791) &
                  + het_rates(k,206))* y(k,206)
         prod(k,158) =rxt(k,685)*y(k,304)*y(k,148)
         loss(k,171) = (rxt(k,748)* y(k,295) + rxt(k,98) + het_rates(k,207))* y(k,207)
         prod(k,171) = (.010_r8*rxt(k,602)*y(k,4) +.130_r8*rxt(k,622)*y(k,7) + &
                 .010_r8*rxt(k,662)*y(k,125))*y(k,158) +.510_r8*rxt(k,751)*y(k,307) &
                 *y(k,258)
         loss(k,115) = (rxt(k,746)* y(k,295) + rxt(k,99) + het_rates(k,208))* y(k,208)
         prod(k,115) =.510_r8*rxt(k,703)*y(k,302)*y(k,258)
         loss(k,116) = (rxt(k,747)* y(k,295) + rxt(k,100) + het_rates(k,209)) &
                 * y(k,209)
         prod(k,116) =.510_r8*rxt(k,721)*y(k,304)*y(k,258)
         loss(k,125) = (rxt(k,757)* y(k,295) + rxt(k,101) + rxt(k,792) &
                  + het_rates(k,210))* y(k,210)
         prod(k,125) =rxt(k,686)*y(k,307)*y(k,148)
         loss(k,117) = (rxt(k,758)* y(k,295) + rxt(k,102) + rxt(k,834) &
                  + het_rates(k,211))* y(k,211)
         prod(k,117) = (.820_r8*rxt(k,687)*y(k,299) +.820_r8*rxt(k,691)*y(k,300)) &
                 *y(k,258)
         loss(k,293) = (rxt(k,759)* y(k,149) +rxt(k,762)* y(k,158) +rxt(k,763) &
                 * y(k,295) + het_rates(k,212))* y(k,212)
         prod(k,293) = (.460_r8*rxt(k,645)*y(k,252) +.310_r8*rxt(k,646)*y(k,253) + &
                 .230_r8*rxt(k,647)*y(k,258) +.860_r8*rxt(k,648)*y(k,281) + &
                 .430_r8*rxt(k,649)*y(k,147) +.460_r8*rxt(k,650)*y(k,149) + &
                 .460_r8*rxt(k,651)*y(k,302) +.460_r8*rxt(k,652)*y(k,304) + &
                 .460_r8*rxt(k,653)*y(k,307))*y(k,281) &
                  + (.120_r8*rxt(k,594)*y(k,252) +.140_r8*rxt(k,595)*y(k,253) + &
                 .060_r8*rxt(k,596)*y(k,258) +.090_r8*rxt(k,597)*y(k,147) + &
                 .120_r8*rxt(k,598)*y(k,149) +.120_r8*rxt(k,599)*y(k,302) + &
                 .120_r8*rxt(k,600)*y(k,304) +.120_r8*rxt(k,601)*y(k,307))*y(k,238) &
                  + (rxt(k,654)*y(k,252) +rxt(k,655)*y(k,253) + &
                 .100_r8*rxt(k,656)*y(k,258) +.770_r8*rxt(k,657)*y(k,147) + &
                 rxt(k,658)*y(k,149) +rxt(k,659)*y(k,302) +rxt(k,660)*y(k,304) + &
                 rxt(k,661)*y(k,307))*y(k,282) + (.270_r8*rxt(k,634)*y(k,252) + &
                 .370_r8*rxt(k,635)*y(k,253) +.200_r8*rxt(k,637)*y(k,147) + &
                 .270_r8*rxt(k,638)*y(k,149) +.270_r8*rxt(k,639)*y(k,302) + &
                 .270_r8*rxt(k,640)*y(k,304) +.270_r8*rxt(k,641)*y(k,307))*y(k,246) &
                  + (.660_r8*rxt(k,662)*y(k,125) +rxt(k,767)*y(k,213))*y(k,158) &
                  + (.100_r8*rxt(k,765)*y(k,258) +.700_r8*rxt(k,766)*y(k,147)) &
                 *y(k,309) +.500_r8*rxt(k,764)*y(k,213)*y(k,149) +rxt(k,91)*y(k,200) &
                  +.460_r8*rxt(k,106)*y(k,218) +.460_r8*rxt(k,108)*y(k,220) &
                  +rxt(k,110)*y(k,222) +rxt(k,112)*y(k,224)
         loss(k,292) = (rxt(k,764)* y(k,149) +rxt(k,767)* y(k,158) +rxt(k,768) &
                 * y(k,295) + het_rates(k,213))* y(k,213)
         prod(k,292) = (1.640_r8*rxt(k,605)*y(k,240) +rxt(k,606)*y(k,252) + &
                 .820_r8*rxt(k,607)*y(k,253) +.500_r8*rxt(k,608)*y(k,258) + &
                 .930_r8*rxt(k,609)*y(k,147) +rxt(k,610)*y(k,149) + &
                 rxt(k,611)*y(k,302) +rxt(k,612)*y(k,304) +rxt(k,613)*y(k,307)) &
                 *y(k,240) + (.950_r8*rxt(k,665)*y(k,252) + &
                 .770_r8*rxt(k,666)*y(k,253) +.480_r8*rxt(k,667)*y(k,258) + &
                 1.540_r8*rxt(k,668)*y(k,290) +.890_r8*rxt(k,669)*y(k,147) + &
                 .950_r8*rxt(k,670)*y(k,149) +.950_r8*rxt(k,671)*y(k,302) + &
                 .950_r8*rxt(k,672)*y(k,304) +.950_r8*rxt(k,673)*y(k,307))*y(k,290) &
                  + (rxt(k,614)*y(k,252) +rxt(k,615)*y(k,253) + &
                 .100_r8*rxt(k,616)*y(k,258) +.700_r8*rxt(k,617)*y(k,147) + &
                 rxt(k,618)*y(k,149) +rxt(k,619)*y(k,302) +rxt(k,620)*y(k,304) + &
                 rxt(k,621)*y(k,307))*y(k,241) + (rxt(k,674)*y(k,252) + &
                 rxt(k,675)*y(k,253) +.100_r8*rxt(k,676)*y(k,258) + &
                 .710_r8*rxt(k,677)*y(k,147) +rxt(k,678)*y(k,149) + &
                 rxt(k,679)*y(k,302) +rxt(k,680)*y(k,304) +rxt(k,681)*y(k,307)) &
                 *y(k,291) + (.870_r8*rxt(k,622)*y(k,7) +rxt(k,682)*y(k,135))*y(k,158) &
                  +rxt(k,92)*y(k,201)
         loss(k,215) = (rxt(k,769)* y(k,295) + rxt(k,103) + rxt(k,835) &
                  + het_rates(k,214))* y(k,214)
         prod(k,215) = (.070_r8*rxt(k,589)*y(k,237) +.070_r8*rxt(k,629)*y(k,245) + &
                 .070_r8*rxt(k,649)*y(k,281) +.070_r8*rxt(k,669)*y(k,290) + &
                 .300_r8*rxt(k,773)*y(k,310) +.300_r8*rxt(k,777)*y(k,311) + &
                 .300_r8*rxt(k,781)*y(k,312) +.300_r8*rxt(k,785)*y(k,313))*y(k,147)
         loss(k,192) = (rxt(k,770)* y(k,295) + rxt(k,104) + rxt(k,836) &
                  + het_rates(k,215))* y(k,215)
         prod(k,192) = (.010_r8*rxt(k,597)*y(k,238) +.300_r8*rxt(k,688)*y(k,299) + &
                 .300_r8*rxt(k,692)*y(k,300) +.300_r8*rxt(k,761)*y(k,308))*y(k,147) &
                  + (.900_r8*rxt(k,772)*y(k,310) +.900_r8*rxt(k,776)*y(k,311) + &
                 .900_r8*rxt(k,780)*y(k,312) +.900_r8*rxt(k,784)*y(k,313))*y(k,258)
         loss(k,203) = (rxt(k,771)* y(k,295) + het_rates(k,216))* y(k,216)
         prod(k,203) = (.040_r8*rxt(k,625)*y(k,245) +.020_r8*rxt(k,626)*y(k,252) + &
                 .020_r8*rxt(k,627)*y(k,253) +.020_r8*rxt(k,628)*y(k,258) + &
                 .020_r8*rxt(k,629)*y(k,147) +.020_r8*rxt(k,630)*y(k,149) + &
                 .020_r8*rxt(k,631)*y(k,302) +.020_r8*rxt(k,632)*y(k,304) + &
                 .020_r8*rxt(k,633)*y(k,307))*y(k,245) &
                  + (.320_r8*rxt(k,634)*y(k,252) +.320_r8*rxt(k,635)*y(k,253) + &
                 .030_r8*rxt(k,636)*y(k,258) +.240_r8*rxt(k,637)*y(k,147) + &
                 .320_r8*rxt(k,638)*y(k,149) +.320_r8*rxt(k,639)*y(k,302) + &
                 .320_r8*rxt(k,640)*y(k,304) +.320_r8*rxt(k,641)*y(k,307))*y(k,246) &
                  +.510_r8*rxt(k,642)*y(k,158)*y(k,17) +.110_r8*rxt(k,595)*y(k,253) &
                 *y(k,238)
         loss(k,194) = (rxt(k,775)* y(k,295) + rxt(k,105) + het_rates(k,217)) &
                 * y(k,217)
         prod(k,194) = (.450_r8*rxt(k,628)*y(k,245) +.100_r8*rxt(k,772)*y(k,310)) &
                 *y(k,258) +.700_r8*rxt(k,773)*y(k,310)*y(k,147)
         loss(k,160) = (rxt(k,774)* y(k,295) + rxt(k,106) + het_rates(k,218)) &
                 * y(k,218)
         prod(k,160) = (.320_r8*rxt(k,647)*y(k,281) +.360_r8*rxt(k,667)*y(k,290)) &
                 *y(k,258)
         loss(k,207) = (rxt(k,779)* y(k,295) + rxt(k,107) + rxt(k,838) &
                  + het_rates(k,219))* y(k,219)
         prod(k,207) = (.300_r8*rxt(k,588)*y(k,237) +.080_r8*rxt(k,628)*y(k,245) + &
                 .100_r8*rxt(k,776)*y(k,311))*y(k,258) +.700_r8*rxt(k,777)*y(k,311) &
                 *y(k,147)
         loss(k,172) = (rxt(k,778)* y(k,295) + rxt(k,108) + rxt(k,837) &
                  + het_rates(k,220))* y(k,220)
         prod(k,172) = (.180_r8*rxt(k,647)*y(k,281) +.160_r8*rxt(k,667)*y(k,290)) &
                 *y(k,258)
         loss(k,244) = (rxt(k,783)* y(k,295) + rxt(k,109) + het_rates(k,221)) &
                 * y(k,221)
         prod(k,244) = (.920_r8*rxt(k,625)*y(k,245) +.450_r8*rxt(k,626)*y(k,252) + &
                 .560_r8*rxt(k,627)*y(k,253) +.230_r8*rxt(k,628)*y(k,258) + &
                 .420_r8*rxt(k,629)*y(k,147) +.450_r8*rxt(k,630)*y(k,149) + &
                 .450_r8*rxt(k,631)*y(k,302) +.450_r8*rxt(k,632)*y(k,304) + &
                 .450_r8*rxt(k,633)*y(k,307))*y(k,245) &
                  + (.100_r8*rxt(k,597)*y(k,238) +.020_r8*rxt(k,637)*y(k,246) + &
                 .300_r8*rxt(k,696)*y(k,301) +.090_r8*rxt(k,741)*y(k,306) + &
                 .700_r8*rxt(k,781)*y(k,312))*y(k,147) + (rxt(k,103) + &
                 rxt(k,769)*y(k,295))*y(k,214) + (rxt(k,104) +rxt(k,770)*y(k,295)) &
                 *y(k,215) + (.090_r8*rxt(k,585)*y(k,237) + &
                 .090_r8*rxt(k,587)*y(k,253))*y(k,237) +.500_r8*rxt(k,105)*y(k,217) &
                  +.100_r8*rxt(k,780)*y(k,312)*y(k,258)
         loss(k,252) = (rxt(k,782)* y(k,295) + rxt(k,110) + het_rates(k,222)) &
                 * y(k,222)
         prod(k,252) = (.350_r8*rxt(k,645)*y(k,252) +.420_r8*rxt(k,646)*y(k,253) + &
                 .180_r8*rxt(k,647)*y(k,258) +.720_r8*rxt(k,648)*y(k,281) + &
                 .330_r8*rxt(k,649)*y(k,147) +.350_r8*rxt(k,650)*y(k,149) + &
                 .350_r8*rxt(k,651)*y(k,302) +.350_r8*rxt(k,652)*y(k,304) + &
                 .350_r8*rxt(k,653)*y(k,307))*y(k,281) &
                  + (.050_r8*rxt(k,665)*y(k,252) +.140_r8*rxt(k,666)*y(k,253) + &
                 .190_r8*rxt(k,668)*y(k,290) +.040_r8*rxt(k,669)*y(k,147) + &
                 .050_r8*rxt(k,670)*y(k,149) +.050_r8*rxt(k,671)*y(k,302) + &
                 .050_r8*rxt(k,672)*y(k,304) +.050_r8*rxt(k,673)*y(k,307))*y(k,290) &
                  + (.020_r8*rxt(k,597)*y(k,238) +.040_r8*rxt(k,637)*y(k,246) + &
                 .060_r8*rxt(k,657)*y(k,282) +.100_r8*rxt(k,677)*y(k,291) + &
                 .120_r8*rxt(k,766)*y(k,309))*y(k,147) +.500_r8*rxt(k,764)*y(k,213) &
                 *y(k,149) +.540_r8*rxt(k,106)*y(k,218)
         loss(k,242) = (rxt(k,787)* y(k,295) + rxt(k,111) + rxt(k,840) &
                  + het_rates(k,223))* y(k,223)
         prod(k,242) = (.140_r8*rxt(k,625)*y(k,245) +.050_r8*rxt(k,626)*y(k,252) + &
                 .080_r8*rxt(k,627)*y(k,253) +.050_r8*rxt(k,629)*y(k,147) + &
                 .050_r8*rxt(k,630)*y(k,149) +.050_r8*rxt(k,631)*y(k,302) + &
                 .050_r8*rxt(k,632)*y(k,304) +.050_r8*rxt(k,633)*y(k,307))*y(k,245) &
                  + (.050_r8*rxt(k,597)*y(k,238) +.060_r8*rxt(k,637)*y(k,246) + &
                 .170_r8*rxt(k,712)*y(k,303) +.300_r8*rxt(k,731)*y(k,305) + &
                 .700_r8*rxt(k,785)*y(k,313))*y(k,147) &
                  + (.270_r8*rxt(k,585)*y(k,237) +.090_r8*rxt(k,587)*y(k,253)) &
                 *y(k,237) +rxt(k,779)*y(k,295)*y(k,219) +.100_r8*rxt(k,784)*y(k,313) &
                 *y(k,258)
         loss(k,253) = (rxt(k,786)* y(k,295) + rxt(k,112) + rxt(k,839) &
                  + het_rates(k,224))* y(k,224)
         prod(k,253) = (.190_r8*rxt(k,645)*y(k,252) +.270_r8*rxt(k,646)*y(k,253) + &
                 .090_r8*rxt(k,647)*y(k,258) +.420_r8*rxt(k,648)*y(k,281) + &
                 .170_r8*rxt(k,649)*y(k,147) +.190_r8*rxt(k,650)*y(k,149) + &
                 .190_r8*rxt(k,651)*y(k,302) +.190_r8*rxt(k,652)*y(k,304) + &
                 .190_r8*rxt(k,653)*y(k,307))*y(k,281) &
                  + (.050_r8*rxt(k,597)*y(k,238) +.130_r8*rxt(k,637)*y(k,246) + &
                 .170_r8*rxt(k,657)*y(k,282) +.190_r8*rxt(k,677)*y(k,291) + &
                 .180_r8*rxt(k,766)*y(k,309))*y(k,147) &
                  + (.090_r8*rxt(k,666)*y(k,253) +.270_r8*rxt(k,668)*y(k,290)) &
                 *y(k,290) +.540_r8*rxt(k,108)*y(k,220)
         loss(k,161) = (rxt(k,789)* y(k,295) + rxt(k,113) + het_rates(k,225)) &
                 * y(k,225)
         prod(k,161) = (.400_r8*rxt(k,596)*y(k,238) +.290_r8*rxt(k,636)*y(k,246) + &
                 rxt(k,695)*y(k,301) +.620_r8*rxt(k,711)*y(k,303))*y(k,258) &
                  + (rxt(k,102) +rxt(k,758)*y(k,295))*y(k,211)
         loss(k,151) = (rxt(k,788)* y(k,295) + rxt(k,114) + het_rates(k,226)) &
                 * y(k,226)
         prod(k,151) = (.180_r8*rxt(k,687)*y(k,299) +.850_r8*rxt(k,730)*y(k,305) + &
                 .470_r8*rxt(k,740)*y(k,306) +.900_r8*rxt(k,760)*y(k,308))*y(k,258) &
                  +.700_r8*rxt(k,688)*y(k,299)*y(k,147)
         loss(k,167) = (rxt(k,573)* y(k,295) + rxt(k,115) + het_rates(k,227)) &
                 * y(k,227)
         prod(k,167) =rxt(k,571)*y(k,314)*y(k,258)
         loss(k,83) = (rxt(k,574)* y(k,295) + het_rates(k,228))* y(k,228)
         prod(k,83) = 0._r8
         loss(k,86) = (rxt(k,576)* y(k,295) + het_rates(k,229))* y(k,229)
         prod(k,86) = 0._r8
         loss(k,178) = (rxt(k,579)* y(k,295) + rxt(k,116) + het_rates(k,230)) &
                 * y(k,230)
         prod(k,178) =rxt(k,577)*y(k,316)*y(k,258)
         loss(k,87) = (rxt(k,582)* y(k,295) + het_rates(k,231))* y(k,231)
         prod(k,87) =.150_r8*rxt(k,576)*y(k,295)*y(k,229)
         loss(k,126) = (rxt(k,583)* y(k,295) + rxt(k,117) + het_rates(k,232)) &
                 * y(k,232)
         prod(k,126) =rxt(k,580)*y(k,318)*y(k,258)
         loss(k,144) = (rxt(k,542)* y(k,147) +rxt(k,570)* y(k,148) +rxt(k,541) &
                 * y(k,258) + het_rates(k,235))* y(k,235)
         prod(k,144) =rxt(k,547)*y(k,295)*y(k,23) +rxt(k,575)*y(k,164)
         loss(k,211) = ((rxt(k,410) +rxt(k,411))* y(k,147) +rxt(k,409)* y(k,258) &
                  + het_rates(k,236))* y(k,236)
         prod(k,211) = (rxt(k,412)*y(k,2) +rxt(k,413)*y(k,15))*y(k,295)
         loss(k,281) = (rxt(k,589)* y(k,147) +rxt(k,590)* y(k,149) + 2._r8*rxt(k,585) &
                 * y(k,237) +rxt(k,586)* y(k,252) +rxt(k,587)* y(k,253) +rxt(k,588) &
                 * y(k,258) +rxt(k,591)* y(k,302) +rxt(k,592)* y(k,304) +rxt(k,593) &
                 * y(k,307) + het_rates(k,237))* y(k,237)
         prod(k,281) =rxt(k,584)*y(k,149)*y(k,4)
         loss(k,287) = (rxt(k,597)* y(k,147) +rxt(k,598)* y(k,149) +rxt(k,594) &
                 * y(k,252) +rxt(k,595)* y(k,253) +rxt(k,596)* y(k,258) +rxt(k,599) &
                 * y(k,302) +rxt(k,600)* y(k,304) +rxt(k,601)* y(k,307) &
                  + het_rates(k,238))* y(k,238)
         prod(k,287) =rxt(k,603)*y(k,295)*y(k,4)
         loss(k,47) = (rxt(k,843)* y(k,147) +rxt(k,842)* y(k,258) + het_rates(k,239)) &
                 * y(k,239)
         prod(k,47) =rxt(k,845)*y(k,295)*y(k,4)
         loss(k,277) = (rxt(k,609)* y(k,147) +rxt(k,610)* y(k,149) + 2._r8*rxt(k,605) &
                 * y(k,240) +rxt(k,606)* y(k,252) +rxt(k,607)* y(k,253) +rxt(k,608) &
                 * y(k,258) +rxt(k,611)* y(k,302) +rxt(k,612)* y(k,304) +rxt(k,613) &
                 * y(k,307) + het_rates(k,240))* y(k,240)
         prod(k,277) =rxt(k,604)*y(k,149)*y(k,7)
         loss(k,286) = (rxt(k,617)* y(k,147) +rxt(k,618)* y(k,149) +rxt(k,614) &
                 * y(k,252) +rxt(k,615)* y(k,253) +rxt(k,616)* y(k,258) +rxt(k,619) &
                 * y(k,302) +rxt(k,620)* y(k,304) +rxt(k,621)* y(k,307) &
                  + het_rates(k,241))* y(k,241)
         prod(k,286) =rxt(k,623)*y(k,295)*y(k,7)
         loss(k,48) = (rxt(k,848)* y(k,147) +rxt(k,847)* y(k,258) + het_rates(k,242)) &
                 * y(k,242)
         prod(k,48) =rxt(k,850)*y(k,295)*y(k,7)
         loss(k,139) = (rxt(k,545)* y(k,147) +rxt(k,544)* y(k,258) + het_rates(k,243)) &
                 * y(k,243)
         prod(k,139) = (.350_r8*rxt(k,543)*y(k,8) +rxt(k,546)*y(k,9))*y(k,295)
         loss(k,49) = (rxt(k,853)* y(k,147) +rxt(k,852)* y(k,258) + het_rates(k,244)) &
                 * y(k,244)
         prod(k,49) =rxt(k,851)*y(k,295)*y(k,8)
         loss(k,290) = (rxt(k,629)* y(k,147) +rxt(k,630)* y(k,149) + 2._r8*rxt(k,625) &
                 * y(k,245) +rxt(k,626)* y(k,252) +rxt(k,627)* y(k,253) +rxt(k,628) &
                 * y(k,258) +rxt(k,631)* y(k,302) +rxt(k,632)* y(k,304) +rxt(k,633) &
                 * y(k,307) + het_rates(k,245))* y(k,245)
         prod(k,290) =rxt(k,624)*y(k,149)*y(k,17) +rxt(k,775)*y(k,295)*y(k,217)
         loss(k,285) = (rxt(k,637)* y(k,147) +rxt(k,638)* y(k,149) +rxt(k,634) &
                 * y(k,252) +rxt(k,635)* y(k,253) +rxt(k,636)* y(k,258) +rxt(k,639) &
                 * y(k,302) +rxt(k,640)* y(k,304) +rxt(k,641)* y(k,307) &
                  + het_rates(k,246))* y(k,246)
         prod(k,285) =rxt(k,643)*y(k,295)*y(k,17)
         loss(k,50) = (rxt(k,856)* y(k,147) +rxt(k,855)* y(k,258) + het_rates(k,247)) &
                 * y(k,247)
         prod(k,50) =rxt(k,858)*y(k,295)*y(k,17)
         loss(k,127) = (rxt(k,550)* y(k,147) +rxt(k,548)* y(k,258) + het_rates(k,248)) &
                 * y(k,248)
         prod(k,127) = (rxt(k,549)*y(k,24) +.070_r8*rxt(k,574)*y(k,228) + &
                 .060_r8*rxt(k,576)*y(k,229))*y(k,295)
         loss(k,225) = (rxt(k,327)* y(k,147) + 2._r8*rxt(k,324)* y(k,249) +rxt(k,325) &
                 * y(k,253) +rxt(k,326)* y(k,258) + het_rates(k,249))* y(k,249)
         prod(k,225) = (rxt(k,330)*y(k,57) +rxt(k,331)*y(k,295))*y(k,29) &
                  +.500_r8*rxt(k,329)*y(k,295)*y(k,28) +rxt(k,76)*y(k,129)
         loss(k,198) = (rxt(k,357)* y(k,147) +rxt(k,355)* y(k,253) +rxt(k,356) &
                 * y(k,258) + het_rates(k,250))* y(k,250)
         prod(k,198) = (rxt(k,358)*y(k,31) +rxt(k,359)*y(k,32))*y(k,295)
         loss(k,164) = (rxt(k,552)* y(k,147) +rxt(k,551)* y(k,258) + het_rates(k,251)) &
                 * y(k,251)
         prod(k,164) = (.400_r8*rxt(k,541)*y(k,258) +rxt(k,542)*y(k,147))*y(k,235) &
                  +rxt(k,553)*y(k,295)*y(k,33) +rxt(k,568)*y(k,165)*y(k,158)
         loss(k,300) = (rxt(k,338)* y(k,147) +rxt(k,351)* y(k,148) +rxt(k,586) &
                 * y(k,237) +rxt(k,594)* y(k,238) +rxt(k,606)* y(k,240) +rxt(k,614) &
                 * y(k,241) +rxt(k,626)* y(k,245) +rxt(k,634)* y(k,246) &
                  + 2._r8*rxt(k,335)* y(k,252) +rxt(k,336)* y(k,253) +rxt(k,337) &
                 * y(k,258) +rxt(k,424)* y(k,261) +rxt(k,430)* y(k,262) +rxt(k,444) &
                 * y(k,267) +rxt(k,448)* y(k,268) +rxt(k,474)* y(k,274) +rxt(k,491) &
                 * y(k,278) +rxt(k,495)* y(k,279) +rxt(k,645)* y(k,281) +rxt(k,654) &
                 * y(k,282) +rxt(k,381)* y(k,284) +rxt(k,388)* y(k,286) +rxt(k,400) &
                 * y(k,289) +rxt(k,665)* y(k,290) +rxt(k,674)* y(k,291) +rxt(k,693) &
                 * y(k,301) +rxt(k,701)* y(k,302) +rxt(k,709)* y(k,303) +rxt(k,719) &
                 * y(k,304) +rxt(k,728)* y(k,305) +rxt(k,749)* y(k,307) &
                  + het_rates(k,252))* y(k,252)
         prod(k,300) = (rxt(k,333)*y(k,46) +.500_r8*rxt(k,340)*y(k,52) + &
                 rxt(k,361)*y(k,50) +.300_r8*rxt(k,363)*y(k,104) + &
                 .560_r8*rxt(k,405)*y(k,134) +.060_r8*rxt(k,414)*y(k,98) + &
                 .060_r8*rxt(k,415)*y(k,99) +.100_r8*rxt(k,502)*y(k,139) + &
                 2.000_r8*rxt(k,737)*y(k,206))*y(k,295) + (rxt(k,739)*y(k,253) + &
                 .530_r8*rxt(k,740)*y(k,258) +.910_r8*rxt(k,741)*y(k,147) + &
                 rxt(k,742)*y(k,149) +rxt(k,743)*y(k,302) +rxt(k,744)*y(k,304) + &
                 rxt(k,745)*y(k,307))*y(k,306) + (.350_r8*rxt(k,388)*y(k,252) + &
                 .350_r8*rxt(k,389)*y(k,253) +.170_r8*rxt(k,390)*y(k,258) + &
                 .700_r8*rxt(k,391)*y(k,286) +.350_r8*rxt(k,392)*y(k,147) + &
                 .350_r8*rxt(k,393)*y(k,149))*y(k,286) &
                  + (.100_r8*rxt(k,385)*y(k,126) +.280_r8*rxt(k,403)*y(k,132) + &
                 .070_r8*rxt(k,486)*y(k,109) +.040_r8*rxt(k,501)*y(k,139) + &
                 .330_r8*rxt(k,662)*y(k,125))*y(k,158) &
                  + (.750_r8*rxt(k,400)*y(k,252) +.880_r8*rxt(k,401)*y(k,253) + &
                 .490_r8*rxt(k,402)*y(k,258) +.760_r8*rxt(k,537)*y(k,147))*y(k,289) &
                  + (.300_r8*rxt(k,368)*y(k,253) +.150_r8*rxt(k,369)*y(k,258) + &
                 rxt(k,370)*y(k,147))*y(k,298) + (rxt(k,35) +rxt(k,360)*y(k,149)) &
                 *y(k,50) + (rxt(k,55) +rxt(k,56))*y(k,104) + (.600_r8*rxt(k,86) + &
                 rxt(k,352))*y(k,163) + (.200_r8*rxt(k,394)*y(k,258) + &
                 rxt(k,395)*y(k,147))*y(k,288) +rxt(k,26)*y(k,14) +rxt(k,332)*y(k,149) &
                 *y(k,46) +rxt(k,34)*y(k,49) +.330_r8*rxt(k,47)*y(k,97) &
                  +.050_r8*rxt(k,48)*y(k,98) +.070_r8*rxt(k,49)*y(k,99) +rxt(k,52) &
                 *y(k,102) +.500_r8*rxt(k,53)*y(k,103) +.350_r8*rxt(k,72)*y(k,126) &
                  +rxt(k,76)*y(k,129) +rxt(k,77)*y(k,130) +.300_r8*rxt(k,79)*y(k,132) &
                  +.750_r8*rxt(k,80)*y(k,133) +.560_r8*rxt(k,81)*y(k,134) +rxt(k,84) &
                 *y(k,151) +rxt(k,89)*y(k,171) +.500_r8*rxt(k,90)*y(k,199)
         loss(k,308) = (rxt(k,225)* y(k,60) +rxt(k,305)* y(k,147) +rxt(k,587) &
                 * y(k,237) +rxt(k,595)* y(k,238) +rxt(k,607)* y(k,240) +rxt(k,615) &
                 * y(k,241) +rxt(k,627)* y(k,245) +rxt(k,635)* y(k,246) +rxt(k,325) &
                 * y(k,249) +rxt(k,355)* y(k,250) +rxt(k,336)* y(k,252) &
                  + 2._r8*(rxt(k,302) +rxt(k,303))* y(k,253) +rxt(k,304)* y(k,258) &
                  +rxt(k,425)* y(k,261) +rxt(k,431)* y(k,262) +rxt(k,445)* y(k,267) &
                  +rxt(k,449)* y(k,268) +rxt(k,475)* y(k,274) +rxt(k,492)* y(k,278) &
                  +rxt(k,496)* y(k,279) +rxt(k,646)* y(k,281) +rxt(k,655)* y(k,282) &
                  +rxt(k,382)* y(k,284) +rxt(k,389)* y(k,286) +rxt(k,401)* y(k,289) &
                  +rxt(k,666)* y(k,290) +rxt(k,675)* y(k,291) +rxt(k,368)* y(k,298) &
                  +rxt(k,694)* y(k,301) +rxt(k,702)* y(k,302) +rxt(k,710)* y(k,303) &
                  +rxt(k,720)* y(k,304) +rxt(k,729)* y(k,305) +rxt(k,739)* y(k,306) &
                  +rxt(k,750)* y(k,307) + het_rates(k,253))* y(k,253)
         prod(k,308) = (2.000_r8*rxt(k,335)*y(k,252) +.900_r8*rxt(k,336)*y(k,253) + &
                 .490_r8*rxt(k,337)*y(k,258) +rxt(k,338)*y(k,147) + &
                 rxt(k,381)*y(k,284) +1.650_r8*rxt(k,388)*y(k,286) + &
                 rxt(k,400)*y(k,289) +rxt(k,424)*y(k,261) +rxt(k,430)*y(k,262) + &
                 rxt(k,444)*y(k,267) +rxt(k,448)*y(k,268) +rxt(k,474)*y(k,274) + &
                 rxt(k,491)*y(k,278) +rxt(k,495)*y(k,279) +rxt(k,586)*y(k,237) + &
                 rxt(k,594)*y(k,238) +rxt(k,606)*y(k,240) +rxt(k,614)*y(k,241) + &
                 rxt(k,626)*y(k,245) +rxt(k,634)*y(k,246) +rxt(k,645)*y(k,281) + &
                 rxt(k,654)*y(k,282) +rxt(k,665)*y(k,290) +rxt(k,674)*y(k,291) + &
                 rxt(k,693)*y(k,301) +rxt(k,701)*y(k,302) +rxt(k,709)*y(k,303) + &
                 rxt(k,719)*y(k,304) +rxt(k,728)*y(k,305) +rxt(k,738)*y(k,306) + &
                 rxt(k,749)*y(k,307))*y(k,252) + (rxt(k,38) +rxt(k,219)*y(k,57) + &
                 rxt(k,275)*y(k,75) +rxt(k,308)*y(k,295) +rxt(k,315)*y(k,294))*y(k,55) &
                  + (.650_r8*rxt(k,389)*y(k,253) +.320_r8*rxt(k,390)*y(k,258) + &
                 1.300_r8*rxt(k,391)*y(k,286) +.650_r8*rxt(k,392)*y(k,147) + &
                 .650_r8*rxt(k,393)*y(k,149))*y(k,286) + (.700_r8*rxt(k,307)*y(k,54) + &
                 rxt(k,339)*y(k,51) +.060_r8*rxt(k,414)*y(k,98) + &
                 .060_r8*rxt(k,415)*y(k,99))*y(k,295) + (.830_r8*rxt(k,556)*y(k,254) + &
                 .170_r8*rxt(k,562)*y(k,287))*y(k,147) + (.280_r8*rxt(k,354)*y(k,30) + &
                 .210_r8*rxt(k,486)*y(k,109))*y(k,158) &
                  + (.330_r8*rxt(k,555)*y(k,254) +.070_r8*rxt(k,561)*y(k,287)) &
                 *y(k,258) +rxt(k,131)*y(k,44) +rxt(k,33)*y(k,46) +rxt(k,133)*y(k,47) &
                  +rxt(k,34)*y(k,49) +rxt(k,36)*y(k,52) +.040_r8*rxt(k,48)*y(k,98) &
                  +.070_r8*rxt(k,49)*y(k,99) +.650_r8*rxt(k,72)*y(k,126) &
                  +.300_r8*rxt(k,79)*y(k,132) +.400_r8*rxt(k,86)*y(k,163)
         loss(k,184) = (rxt(k,556)* y(k,147) +rxt(k,557)* y(k,148) +rxt(k,555) &
                 * y(k,258) + het_rates(k,254))* y(k,254)
         prod(k,184) =.600_r8*rxt(k,24)*y(k,12)
         loss(k,152) = ((rxt(k,377) +rxt(k,378))* y(k,147) + het_rates(k,255)) &
                 * y(k,255)
         prod(k,152) =rxt(k,375)*y(k,295)*y(k,16)
         loss(k,101) = ( + rxt(k,343) + rxt(k,344) + het_rates(k,256))* y(k,256)
         prod(k,101) =rxt(k,42)*y(k,74) +.750_r8*rxt(k,342)*y(k,257)*y(k,147)
         loss(k,179) = (rxt(k,342)* y(k,147) +rxt(k,341)* y(k,258) + het_rates(k,257)) &
                 * y(k,257)
         prod(k,179) =rxt(k,350)*y(k,295)*y(k,26)
         loss(k,307) = (rxt(k,255)* y(k,18) +rxt(k,261)* y(k,20) +rxt(k,298)* y(k,43) &
                  + (rxt(k,222) +rxt(k,223))* y(k,57) +rxt(k,229)* y(k,60) &
                  + (rxt(k,178) +rxt(k,179) +rxt(k,180))* y(k,78) +rxt(k,207) &
                 * y(k,147) +rxt(k,212)* y(k,148) +rxt(k,202)* y(k,149) +rxt(k,182) &
                 * y(k,157) +rxt(k,183)* y(k,158) +rxt(k,541)* y(k,235) +rxt(k,409) &
                 * y(k,236) +rxt(k,588)* y(k,237) +rxt(k,596)* y(k,238) +rxt(k,608) &
                 * y(k,240) +rxt(k,616)* y(k,241) +rxt(k,544)* y(k,243) +rxt(k,628) &
                 * y(k,245) +rxt(k,636)* y(k,246) +rxt(k,548)* y(k,248) +rxt(k,326) &
                 * y(k,249) +rxt(k,356)* y(k,250) +rxt(k,551)* y(k,251) +rxt(k,337) &
                 * y(k,252) +rxt(k,304)* y(k,253) +rxt(k,555)* y(k,254) +rxt(k,341) &
                 * y(k,257) + 2._r8*rxt(k,192)* y(k,258) +rxt(k,312)* y(k,259) &
                  +rxt(k,421)* y(k,260) +rxt(k,426)* y(k,261) +rxt(k,432)* y(k,262) &
                  +rxt(k,446)* y(k,267) +rxt(k,450)* y(k,268) +rxt(k,457)* y(k,269) &
                  +rxt(k,461)* y(k,270) +rxt(k,464)* y(k,271) +rxt(k,467)* y(k,272) &
                  +rxt(k,471)* y(k,273) +rxt(k,476)* y(k,274) +rxt(k,479)* y(k,275) &
                  +rxt(k,482)* y(k,276) +rxt(k,493)* y(k,278) +rxt(k,497)* y(k,279) &
                  +rxt(k,647)* y(k,281) +rxt(k,656)* y(k,282) +rxt(k,383)* y(k,284) &
                  +rxt(k,558)* y(k,285) +rxt(k,390)* y(k,286) +rxt(k,561)* y(k,287) &
                  +rxt(k,394)* y(k,288) +rxt(k,402)* y(k,289) +rxt(k,667)* y(k,290) &
                  +rxt(k,676)* y(k,291) +rxt(k,499)* y(k,293) +rxt(k,187)* y(k,295) &
                  +rxt(k,564)* y(k,296) +rxt(k,365)* y(k,297) +rxt(k,369)* y(k,298) &
                  +rxt(k,687)* y(k,299) +rxt(k,691)* y(k,300) +rxt(k,695)* y(k,301) &
                  +rxt(k,703)* y(k,302) +rxt(k,711)* y(k,303) +rxt(k,721)* y(k,304) &
                  +rxt(k,730)* y(k,305) +rxt(k,740)* y(k,306) +rxt(k,751)* y(k,307) &
                  +rxt(k,760)* y(k,308) +rxt(k,765)* y(k,309) +rxt(k,772)* y(k,310) &
                  +rxt(k,776)* y(k,311) +rxt(k,780)* y(k,312) +rxt(k,784)* y(k,313) &
                  +rxt(k,571)* y(k,314) +rxt(k,577)* y(k,316) +rxt(k,580)* y(k,318) &
                  + rxt(k,812) + het_rates(k,258))* y(k,258)
         prod(k,307) = (rxt(k,305)*y(k,253) +rxt(k,314)*y(k,259) + &
                 rxt(k,327)*y(k,249) +.250_r8*rxt(k,342)*y(k,257) + &
                 rxt(k,357)*y(k,250) +rxt(k,366)*y(k,297) +rxt(k,377)*y(k,255) + &
                 rxt(k,410)*y(k,236) +rxt(k,503)*y(k,260) +rxt(k,505)*y(k,261) + &
                 rxt(k,507)*y(k,262) +.450_r8*rxt(k,509)*y(k,267) + &
                 .450_r8*rxt(k,511)*y(k,268) +rxt(k,513)*y(k,269) + &
                 .270_r8*rxt(k,515)*y(k,270) +rxt(k,517)*y(k,271) + &
                 rxt(k,519)*y(k,272) +rxt(k,521)*y(k,273) + &
                 .540_r8*rxt(k,523)*y(k,274) +.530_r8*rxt(k,525)*y(k,275) + &
                 .960_r8*rxt(k,527)*y(k,276) +.450_r8*rxt(k,530)*y(k,278) + &
                 .450_r8*rxt(k,533)*y(k,279) +rxt(k,535)*y(k,284) + &
                 .240_r8*rxt(k,537)*y(k,289) +rxt(k,539)*y(k,293) + &
                 rxt(k,545)*y(k,243) +rxt(k,550)*y(k,248) + &
                 .170_r8*rxt(k,556)*y(k,254) +.400_r8*rxt(k,559)*y(k,285) + &
                 .830_r8*rxt(k,562)*y(k,287) +rxt(k,565)*y(k,296) + &
                 rxt(k,572)*y(k,314) +rxt(k,578)*y(k,316) +rxt(k,581)*y(k,318) + &
                 .770_r8*rxt(k,597)*y(k,238) +.700_r8*rxt(k,617)*y(k,241) + &
                 .470_r8*rxt(k,629)*y(k,245) +.750_r8*rxt(k,637)*y(k,246) + &
                 .500_r8*rxt(k,649)*y(k,281) +.770_r8*rxt(k,657)*y(k,282) + &
                 .040_r8*rxt(k,669)*y(k,290) +.710_r8*rxt(k,677)*y(k,291) + &
                 .700_r8*rxt(k,688)*y(k,299) +.700_r8*rxt(k,692)*y(k,300) + &
                 .910_r8*rxt(k,741)*y(k,306) +.700_r8*rxt(k,761)*y(k,308) + &
                 .700_r8*rxt(k,766)*y(k,309) +.700_r8*rxt(k,773)*y(k,310) + &
                 .700_r8*rxt(k,777)*y(k,311) +.700_r8*rxt(k,781)*y(k,312) + &
                 .700_r8*rxt(k,785)*y(k,313))*y(k,147) + (rxt(k,186)*y(k,81) + &
                 rxt(k,189)*y(k,158) +rxt(k,205)*y(k,149) +rxt(k,236)*y(k,60) + &
                 rxt(k,266)*y(k,20) +rxt(k,284)*y(k,44) +rxt(k,287)*y(k,47) + &
                 rxt(k,306)*y(k,53) +rxt(k,309)*y(k,88) +rxt(k,310)*y(k,90) + &
                 .500_r8*rxt(k,311)*y(k,92) +rxt(k,319)*y(k,63) + &
                 .350_r8*rxt(k,321)*y(k,25) +rxt(k,328)*y(k,27) +rxt(k,334)*y(k,48) + &
                 rxt(k,345)*y(k,76) +rxt(k,346)*y(k,77) +.110_r8*rxt(k,347)*y(k,89) + &
                 rxt(k,362)*y(k,102) +rxt(k,379)*y(k,97) + &
                 .500_r8*rxt(k,380)*y(k,127) +rxt(k,399)*y(k,133) + &
                 .440_r8*rxt(k,405)*y(k,134) +.510_r8*rxt(k,414)*y(k,98) + &
                 .410_r8*rxt(k,415)*y(k,99) +.320_r8*rxt(k,418)*y(k,103) + &
                 .190_r8*rxt(k,420)*y(k,106) +.400_r8*rxt(k,423)*y(k,108) + &
                 rxt(k,453)*y(k,110) +rxt(k,455)*y(k,113) + &
                 .040_r8*rxt(k,460)*y(k,115) +.030_r8*rxt(k,470)*y(k,118) + &
                 .050_r8*rxt(k,472)*y(k,119) +rxt(k,488)*y(k,122) + &
                 .180_r8*rxt(k,489)*y(k,123) +.630_r8*rxt(k,502)*y(k,139) + &
                 .650_r8*rxt(k,543)*y(k,8) +.730_r8*rxt(k,554)*y(k,67) + &
                 .800_r8*rxt(k,566)*y(k,166) +.280_r8*rxt(k,574)*y(k,228) + &
                 .380_r8*rxt(k,576)*y(k,229) +.630_r8*rxt(k,582)*y(k,231) + &
                 rxt(k,718)*y(k,204) +rxt(k,737)*y(k,206) +rxt(k,798)*y(k,175) + &
                 .500_r8*rxt(k,808)*y(k,69))*y(k,295) + (rxt(k,225)*y(k,60) + &
                 2.000_r8*rxt(k,302)*y(k,253) +rxt(k,325)*y(k,249) + &
                 .900_r8*rxt(k,336)*y(k,252) +rxt(k,355)*y(k,250) + &
                 .300_r8*rxt(k,368)*y(k,298) +1.500_r8*rxt(k,382)*y(k,284) + &
                 rxt(k,389)*y(k,286) +.620_r8*rxt(k,401)*y(k,289) + &
                 1.500_r8*rxt(k,425)*y(k,261) +rxt(k,431)*y(k,262) + &
                 .720_r8*rxt(k,445)*y(k,267) +.720_r8*rxt(k,449)*y(k,268) + &
                 .400_r8*rxt(k,475)*y(k,274) +.720_r8*rxt(k,492)*y(k,278) + &
                 .720_r8*rxt(k,496)*y(k,279) +.820_r8*rxt(k,587)*y(k,237) + &
                 1.160_r8*rxt(k,595)*y(k,238) +.820_r8*rxt(k,607)*y(k,240) + &
                 rxt(k,615)*y(k,241) +1.100_r8*rxt(k,627)*y(k,245) + &
                 1.500_r8*rxt(k,635)*y(k,246) +1.010_r8*rxt(k,646)*y(k,281) + &
                 rxt(k,655)*y(k,282) +.870_r8*rxt(k,666)*y(k,290) + &
                 rxt(k,675)*y(k,291) +.500_r8*rxt(k,694)*y(k,301) + &
                 rxt(k,702)*y(k,302) +rxt(k,710)*y(k,303) +rxt(k,720)*y(k,304) + &
                 rxt(k,729)*y(k,305) +2.000_r8*rxt(k,739)*y(k,306) + &
                 rxt(k,750)*y(k,307))*y(k,253) + (.200_r8*rxt(k,312)*y(k,259) + &
                 .590_r8*rxt(k,383)*y(k,284) +.180_r8*rxt(k,402)*y(k,289) + &
                 .650_r8*rxt(k,421)*y(k,260) +.060_r8*rxt(k,426)*y(k,261) + &
                 .060_r8*rxt(k,432)*y(k,262) +.580_r8*rxt(k,457)*y(k,269) + &
                 .060_r8*rxt(k,461)*y(k,270) +.600_r8*rxt(k,464)*y(k,271) + &
                 .500_r8*rxt(k,467)*y(k,272) +.400_r8*rxt(k,471)*y(k,273) + &
                 .170_r8*rxt(k,479)*y(k,275) +.800_r8*rxt(k,482)*y(k,276) + &
                 .800_r8*rxt(k,499)*y(k,293) +.070_r8*rxt(k,555)*y(k,254) + &
                 .160_r8*rxt(k,558)*y(k,285) +.330_r8*rxt(k,561)*y(k,287) + &
                 .480_r8*rxt(k,596)*y(k,238) +.100_r8*rxt(k,616)*y(k,241) + &
                 .030_r8*rxt(k,636)*y(k,246) +.270_r8*rxt(k,647)*y(k,281) + &
                 .100_r8*rxt(k,656)*y(k,282) +.100_r8*rxt(k,676)*y(k,291) + &
                 .180_r8*rxt(k,687)*y(k,299) +.180_r8*rxt(k,691)*y(k,300) + &
                 .530_r8*rxt(k,740)*y(k,306) +.100_r8*rxt(k,760)*y(k,308) + &
                 .100_r8*rxt(k,765)*y(k,309) +.100_r8*rxt(k,772)*y(k,310) + &
                 .100_r8*rxt(k,776)*y(k,311) +.100_r8*rxt(k,780)*y(k,312) + &
                 .100_r8*rxt(k,784)*y(k,313))*y(k,258) + (rxt(k,381)*y(k,284) + &
                 .250_r8*rxt(k,400)*y(k,289) +rxt(k,424)*y(k,261) + &
                 rxt(k,430)*y(k,262) +.450_r8*rxt(k,444)*y(k,267) + &
                 .450_r8*rxt(k,448)*y(k,268) +.540_r8*rxt(k,474)*y(k,274) + &
                 .450_r8*rxt(k,491)*y(k,278) +.450_r8*rxt(k,495)*y(k,279) + &
                 rxt(k,594)*y(k,238) +rxt(k,614)*y(k,241) + &
                 .500_r8*rxt(k,626)*y(k,245) +rxt(k,634)*y(k,246) + &
                 .540_r8*rxt(k,645)*y(k,281) +rxt(k,654)*y(k,282) + &
                 .050_r8*rxt(k,665)*y(k,290) +rxt(k,674)*y(k,291) + &
                 rxt(k,738)*y(k,306))*y(k,252) + (rxt(k,299)*y(k,43) + &
                 .540_r8*rxt(k,478)*y(k,274) +rxt(k,598)*y(k,238) + &
                 rxt(k,618)*y(k,241) +.500_r8*rxt(k,630)*y(k,245) + &
                 rxt(k,638)*y(k,246) +.540_r8*rxt(k,650)*y(k,281) + &
                 rxt(k,658)*y(k,282) +.050_r8*rxt(k,670)*y(k,290) + &
                 rxt(k,678)*y(k,291) +rxt(k,742)*y(k,306) + &
                 .500_r8*rxt(k,764)*y(k,213))*y(k,149) + (.130_r8*rxt(k,323)*y(k,26) + &
                 .280_r8*rxt(k,354)*y(k,30) +.140_r8*rxt(k,385)*y(k,126) + &
                 .280_r8*rxt(k,403)*y(k,132) +.170_r8*rxt(k,459)*y(k,115) + &
                 .170_r8*rxt(k,469)*y(k,118) +.420_r8*rxt(k,486)*y(k,109) + &
                 .130_r8*rxt(k,501)*y(k,139) +.170_r8*rxt(k,602)*y(k,4) + &
                 .080_r8*rxt(k,622)*y(k,7) +.630_r8*rxt(k,682)*y(k,135))*y(k,158) &
                  + (rxt(k,599)*y(k,238) +rxt(k,619)*y(k,241) + &
                 .500_r8*rxt(k,631)*y(k,245) +rxt(k,639)*y(k,246) + &
                 .540_r8*rxt(k,651)*y(k,281) +rxt(k,659)*y(k,282) + &
                 .050_r8*rxt(k,671)*y(k,290) +rxt(k,679)*y(k,291) + &
                 rxt(k,743)*y(k,306))*y(k,302) + (rxt(k,600)*y(k,238) + &
                 rxt(k,620)*y(k,241) +.500_r8*rxt(k,632)*y(k,245) + &
                 rxt(k,640)*y(k,246) +.540_r8*rxt(k,652)*y(k,281) + &
                 rxt(k,660)*y(k,282) +.050_r8*rxt(k,672)*y(k,290) + &
                 rxt(k,680)*y(k,291) +rxt(k,744)*y(k,306))*y(k,304) &
                  + (rxt(k,601)*y(k,238) +rxt(k,621)*y(k,241) + &
                 .500_r8*rxt(k,633)*y(k,245) +rxt(k,641)*y(k,246) + &
                 .540_r8*rxt(k,653)*y(k,281) +rxt(k,661)*y(k,282) + &
                 .050_r8*rxt(k,673)*y(k,290) +rxt(k,681)*y(k,291) + &
                 rxt(k,745)*y(k,306))*y(k,307) + (rxt(k,218)*y(k,43) + &
                 rxt(k,221)*y(k,81) +rxt(k,283)*y(k,44) +rxt(k,286)*y(k,47))*y(k,57) &
                  + (rxt(k,254)*y(k,18) +rxt(k,300)*y(k,157))*y(k,43) + (rxt(k,11) + &
                 rxt(k,216))*y(k,94) + (1.500_r8*rxt(k,53) +rxt(k,54))*y(k,103) &
                  + (rxt(k,72) +rxt(k,73))*y(k,126) + (rxt(k,343) +rxt(k,344)) &
                 *y(k,256) +rxt(k,19)*y(k,1) +.900_r8*rxt(k,20)*y(k,2) +rxt(k,21) &
                 *y(k,9) +1.500_r8*rxt(k,22)*y(k,10) +rxt(k,23)*y(k,11) &
                  +.600_r8*rxt(k,24)*y(k,12) +.600_r8*rxt(k,25)*y(k,13) +rxt(k,26) &
                 *y(k,14) +rxt(k,27)*y(k,24) +rxt(k,28)*y(k,28) +rxt(k,29)*y(k,31) &
                  +rxt(k,33)*y(k,46) +rxt(k,35)*y(k,50) +rxt(k,316)*y(k,294)*y(k,55) &
                  +.500_r8*rxt(k,41)*y(k,68) +2.000_r8*rxt(k,43)*y(k,76) &
                  +2.000_r8*rxt(k,44)*y(k,77) +rxt(k,181)*y(k,78) +rxt(k,177)*y(k,157) &
                 *y(k,81) +rxt(k,45)*y(k,89) +.670_r8*rxt(k,47)*y(k,97) &
                  +.620_r8*rxt(k,48)*y(k,98) +.560_r8*rxt(k,49)*y(k,99) +rxt(k,50) &
                 *y(k,100) +rxt(k,51)*y(k,101) +rxt(k,52)*y(k,102) +rxt(k,57)*y(k,107) &
                  +rxt(k,58)*y(k,108) +rxt(k,63)*y(k,114) +.450_r8*rxt(k,64)*y(k,115) &
                  +rxt(k,65)*y(k,116) +rxt(k,66)*y(k,117) +.450_r8*rxt(k,67)*y(k,118) &
                  +rxt(k,68)*y(k,119) +rxt(k,70)*y(k,121) +rxt(k,71)*y(k,123) &
                  +1.250_r8*rxt(k,74)*y(k,127) +rxt(k,75)*y(k,128) +.500_r8*rxt(k,80) &
                 *y(k,133) +.440_r8*rxt(k,81)*y(k,134) +rxt(k,82)*y(k,139) +rxt(k,83) &
                 *y(k,150) +rxt(k,87)*y(k,167) +rxt(k,88)*y(k,170) +rxt(k,90)*y(k,199) &
                  +rxt(k,91)*y(k,200) +rxt(k,92)*y(k,201) +rxt(k,93)*y(k,202) &
                  +rxt(k,94)*y(k,203) +rxt(k,96)*y(k,205) +rxt(k,102)*y(k,211) &
                  +rxt(k,103)*y(k,214) +rxt(k,104)*y(k,215) +.500_r8*rxt(k,105) &
                 *y(k,217) +.540_r8*rxt(k,106)*y(k,218) +.540_r8*rxt(k,108)*y(k,220) &
                  +rxt(k,109)*y(k,221) +rxt(k,110)*y(k,222) +rxt(k,111)*y(k,223) &
                  +rxt(k,112)*y(k,224) +rxt(k,113)*y(k,225) +rxt(k,114)*y(k,226) &
                  +rxt(k,115)*y(k,227) +rxt(k,116)*y(k,230) +rxt(k,117)*y(k,232) &
                  +.940_r8*rxt(k,625)*y(k,245)*y(k,245) +1.200_r8*rxt(k,324)*y(k,249) &
                 *y(k,249) +rxt(k,313)*y(k,259) +rxt(k,458)*y(k,269) +rxt(k,462) &
                 *y(k,270) +rxt(k,465)*y(k,271) +rxt(k,468)*y(k,272) &
                  +.400_r8*rxt(k,477)*y(k,274)*y(k,274) +.400_r8*rxt(k,529)*y(k,278) &
                  +.400_r8*rxt(k,532)*y(k,279) +.990_r8*rxt(k,648)*y(k,281)*y(k,281)
         loss(k,162) = (rxt(k,314)* y(k,147) +rxt(k,312)* y(k,258) + rxt(k,313) &
                  + het_rates(k,259))* y(k,259)
         prod(k,162) =rxt(k,298)*y(k,258)*y(k,43)
         loss(k,221) = ((rxt(k,503) +rxt(k,504))* y(k,147) +rxt(k,421)* y(k,258) &
                  + het_rates(k,260))* y(k,260)
         prod(k,221) = (.320_r8*rxt(k,418)*y(k,103) +.810_r8*rxt(k,420)*y(k,106)) &
                 *y(k,295)
         loss(k,267) = ((rxt(k,505) +rxt(k,506))* y(k,147) +rxt(k,424)* y(k,252) &
                  +rxt(k,425)* y(k,253) +rxt(k,426)* y(k,258) + rxt(k,427) &
                  + rxt(k,428) + rxt(k,429) + het_rates(k,261))* y(k,261)
         prod(k,267) =.530_r8*rxt(k,489)*y(k,295)*y(k,123) +rxt(k,436)*y(k,263) &
                  +rxt(k,438)*y(k,264)
         loss(k,268) = ((rxt(k,507) +rxt(k,508))* y(k,147) +rxt(k,430)* y(k,252) &
                  +rxt(k,431)* y(k,253) +rxt(k,432)* y(k,258) + rxt(k,433) &
                  + rxt(k,434) + rxt(k,435) + het_rates(k,262))* y(k,262)
         prod(k,268) =.160_r8*rxt(k,489)*y(k,295)*y(k,123) +rxt(k,440)*y(k,265) &
                  +rxt(k,442)*y(k,266)
         loss(k,94) = ( + rxt(k,436) + rxt(k,437) + het_rates(k,263))* y(k,263)
         prod(k,94) =.315_r8*rxt(k,487)*y(k,295)*y(k,109) +rxt(k,428)*y(k,261) &
                  +rxt(k,494)*y(k,278)
         loss(k,95) = ( + rxt(k,438) + rxt(k,439) + het_rates(k,264))* y(k,264)
         prod(k,95) =.315_r8*rxt(k,487)*y(k,295)*y(k,109) +rxt(k,429)*y(k,261) &
                  +rxt(k,447)*y(k,267)
         loss(k,96) = ( + rxt(k,440) + rxt(k,441) + het_rates(k,265))* y(k,265)
         prod(k,96) =.259_r8*rxt(k,487)*y(k,295)*y(k,109) +rxt(k,434)*y(k,262) &
                  +rxt(k,498)*y(k,279)
         loss(k,97) = ( + rxt(k,442) + rxt(k,443) + het_rates(k,266))* y(k,266)
         prod(k,97) =.111_r8*rxt(k,487)*y(k,295)*y(k,109) +rxt(k,435)*y(k,262) &
                  +rxt(k,451)*y(k,268)
         loss(k,256) = ((rxt(k,509) +rxt(k,510))* y(k,147) +rxt(k,444)* y(k,252) &
                  +rxt(k,445)* y(k,253) +rxt(k,446)* y(k,258) + rxt(k,447) &
                  + het_rates(k,267))* y(k,267)
         prod(k,256) =rxt(k,439)*y(k,264)
         loss(k,257) = ((rxt(k,511) +rxt(k,512))* y(k,147) +rxt(k,448)* y(k,252) &
                  +rxt(k,449)* y(k,253) +rxt(k,450)* y(k,258) + rxt(k,451) &
                  + het_rates(k,268))* y(k,268)
         prod(k,257) =rxt(k,443)*y(k,266)
         loss(k,205) = ((rxt(k,513) +rxt(k,514))* y(k,147) +rxt(k,457)* y(k,258) &
                  + rxt(k,458) + het_rates(k,269))* y(k,269)
         prod(k,205) =.820_r8*rxt(k,460)*y(k,295)*y(k,115)
         loss(k,210) = ((rxt(k,515) +rxt(k,516))* y(k,147) +rxt(k,461)* y(k,258) &
                  + rxt(k,462) + het_rates(k,270))* y(k,270)
         prod(k,210) =.850_r8*rxt(k,463)*y(k,295)*y(k,116)
         loss(k,200) = ((rxt(k,517) +rxt(k,518))* y(k,147) +rxt(k,464)* y(k,258) &
                  + rxt(k,465) + het_rates(k,271))* y(k,271)
         prod(k,200) =.870_r8*rxt(k,466)*y(k,295)*y(k,117)
         loss(k,206) = ((rxt(k,519) +rxt(k,520))* y(k,147) +rxt(k,467)* y(k,258) &
                  + rxt(k,468) + het_rates(k,272))* y(k,272)
         prod(k,206) =.890_r8*rxt(k,470)*y(k,295)*y(k,118)
         loss(k,231) = ((rxt(k,521) +rxt(k,522))* y(k,147) +rxt(k,471)* y(k,258) &
                  + het_rates(k,273))* y(k,273)
         prod(k,231) =.920_r8*rxt(k,472)*y(k,295)*y(k,119)
         loss(k,275) = ((rxt(k,523) +rxt(k,524))* y(k,147) +rxt(k,478)* y(k,149) &
                  +rxt(k,474)* y(k,252) +rxt(k,475)* y(k,253) +rxt(k,476)* y(k,258) &
                  + 2._r8*rxt(k,477)* y(k,274) + het_rates(k,274))* y(k,274)
         prod(k,275) = (.170_r8*rxt(k,481)*y(k,120) +.070_r8*rxt(k,485)*y(k,121)) &
                 *y(k,295) +rxt(k,473)*y(k,149)*y(k,109)
         loss(k,222) = ((rxt(k,525) +rxt(k,526))* y(k,147) +rxt(k,479)* y(k,258) &
                  + rxt(k,480) + het_rates(k,275))* y(k,275)
         prod(k,222) =.410_r8*rxt(k,481)*y(k,295)*y(k,120)
         loss(k,226) = ((rxt(k,527) +rxt(k,528))* y(k,147) +rxt(k,482)* y(k,258) &
                  + rxt(k,483) + het_rates(k,276))* y(k,276)
         prod(k,226) =.570_r8*rxt(k,485)*y(k,295)*y(k,121)
         loss(k,51) = (rxt(k,861)* y(k,147) +rxt(k,860)* y(k,258) + het_rates(k,277)) &
                 * y(k,277)
         prod(k,51) =rxt(k,863)*y(k,295)*y(k,109)
         loss(k,262) = ((rxt(k,530) +rxt(k,531))* y(k,147) +rxt(k,491)* y(k,252) &
                  +rxt(k,492)* y(k,253) +rxt(k,493)* y(k,258) + rxt(k,494) &
                  + rxt(k,529) + het_rates(k,278))* y(k,278)
         prod(k,262) =rxt(k,437)*y(k,263)
         loss(k,261) = ((rxt(k,533) +rxt(k,534))* y(k,147) +rxt(k,495)* y(k,252) &
                  +rxt(k,496)* y(k,253) +rxt(k,497)* y(k,258) + rxt(k,498) &
                  + rxt(k,532) + het_rates(k,279))* y(k,279)
         prod(k,261) =rxt(k,441)*y(k,265)
         loss(k,52) = (rxt(k,865)* y(k,147) +rxt(k,864)* y(k,258) + het_rates(k,280)) &
                 * y(k,280)
         prod(k,52) =rxt(k,866)*y(k,295)*y(k,124)
         loss(k,289) = (rxt(k,649)* y(k,147) +rxt(k,650)* y(k,149) +rxt(k,645) &
                 * y(k,252) +rxt(k,646)* y(k,253) +rxt(k,647)* y(k,258) &
                  + 2._r8*rxt(k,648)* y(k,281) +rxt(k,651)* y(k,302) +rxt(k,652) &
                 * y(k,304) +rxt(k,653)* y(k,307) + het_rates(k,281))* y(k,281)
         prod(k,289) =rxt(k,644)*y(k,149)*y(k,125)
         loss(k,283) = (rxt(k,657)* y(k,147) +rxt(k,658)* y(k,149) +rxt(k,654) &
                 * y(k,252) +rxt(k,655)* y(k,253) +rxt(k,656)* y(k,258) +rxt(k,659) &
                 * y(k,302) +rxt(k,660)* y(k,304) +rxt(k,661)* y(k,307) &
                  + het_rates(k,282))* y(k,282)
         prod(k,283) =rxt(k,663)*y(k,295)*y(k,125)
         loss(k,54) = (rxt(k,869)* y(k,147) +rxt(k,868)* y(k,258) + het_rates(k,283)) &
                 * y(k,283)
         prod(k,54) =rxt(k,871)*y(k,295)*y(k,125)
         loss(k,243) = ((rxt(k,535) +rxt(k,536))* y(k,147) +rxt(k,381)* y(k,252) &
                  +rxt(k,382)* y(k,253) +rxt(k,383)* y(k,258) + rxt(k,384) &
                  + het_rates(k,284))* y(k,284)
         prod(k,243) =.190_r8*rxt(k,49)*y(k,99) +.550_r8*rxt(k,386)*y(k,295)*y(k,126)
         loss(k,180) = (rxt(k,559)* y(k,147) +rxt(k,560)* y(k,148) +rxt(k,558) &
                 * y(k,258) + het_rates(k,285))* y(k,285)
         prod(k,180) =.600_r8*rxt(k,23)*y(k,11)
         loss(k,248) = (rxt(k,392)* y(k,147) +rxt(k,406)* y(k,148) +rxt(k,393) &
                 * y(k,149) +rxt(k,388)* y(k,252) +rxt(k,389)* y(k,253) +rxt(k,390) &
                 * y(k,258) + 2._r8*rxt(k,391)* y(k,286) + het_rates(k,286))* y(k,286)
         prod(k,248) = (rxt(k,73) +.450_r8*rxt(k,386)*y(k,295))*y(k,126) &
                  + (rxt(k,78) +rxt(k,407))*y(k,131)
         loss(k,188) = (rxt(k,562)* y(k,147) +rxt(k,563)* y(k,148) +rxt(k,561) &
                 * y(k,258) + het_rates(k,287))* y(k,287)
         prod(k,188) =.600_r8*rxt(k,25)*y(k,13)
         loss(k,169) = (rxt(k,395)* y(k,147) +rxt(k,394)* y(k,258) + het_rates(k,288)) &
                 * y(k,288)
         prod(k,169) = (rxt(k,396)*y(k,129) +rxt(k,397)*y(k,130))*y(k,295)
         loss(k,240) = ((rxt(k,537) +rxt(k,538))* y(k,147) +rxt(k,400)* y(k,252) &
                  +rxt(k,401)* y(k,253) +rxt(k,402)* y(k,258) + het_rates(k,289)) &
                 * y(k,289)
         prod(k,240) =.230_r8*rxt(k,48)*y(k,98) +rxt(k,404)*y(k,295)*y(k,132)
         loss(k,291) = (rxt(k,669)* y(k,147) +rxt(k,670)* y(k,149) +rxt(k,665) &
                 * y(k,252) +rxt(k,666)* y(k,253) +rxt(k,667)* y(k,258) &
                  + 2._r8*rxt(k,668)* y(k,290) +rxt(k,671)* y(k,302) +rxt(k,672) &
                 * y(k,304) +rxt(k,673)* y(k,307) + het_rates(k,290))* y(k,290)
         prod(k,291) =rxt(k,664)*y(k,149)*y(k,135)
         loss(k,288) = (rxt(k,677)* y(k,147) +rxt(k,678)* y(k,149) +rxt(k,674) &
                 * y(k,252) +rxt(k,675)* y(k,253) +rxt(k,676)* y(k,258) +rxt(k,679) &
                 * y(k,302) +rxt(k,680)* y(k,304) +rxt(k,681)* y(k,307) &
                  + het_rates(k,291))* y(k,291)
         prod(k,288) =rxt(k,683)*y(k,295)*y(k,135)
         loss(k,55) = (rxt(k,874)* y(k,147) +rxt(k,873)* y(k,258) + het_rates(k,292)) &
                 * y(k,292)
         prod(k,55) =rxt(k,876)*y(k,295)*y(k,135)
         loss(k,249) = ((rxt(k,539) +rxt(k,540))* y(k,147) +rxt(k,499)* y(k,258) &
                  + rxt(k,500) + het_rates(k,293))* y(k,293)
         prod(k,249) = (.400_r8*rxt(k,422)*y(k,107) +.350_r8*rxt(k,423)*y(k,108) + &
                 .230_r8*rxt(k,502)*y(k,139))*y(k,295)
         loss(k,314) = (rxt(k,243)* y(k,34) +rxt(k,244)* y(k,35) +rxt(k,270)* y(k,36) &
                  +rxt(k,245)* y(k,37) +rxt(k,246)* y(k,38) +rxt(k,247)* y(k,39) &
                  +rxt(k,248)* y(k,40) +rxt(k,249)* y(k,41) +rxt(k,293)* y(k,42) &
                  +rxt(k,294)* y(k,44) + (rxt(k,315) +rxt(k,316) +rxt(k,317))* y(k,55) &
                  +rxt(k,271)* y(k,56) +rxt(k,279)* y(k,65) +rxt(k,280)* y(k,66) &
                  +rxt(k,168)* y(k,79) +rxt(k,272)* y(k,80) + (rxt(k,273) +rxt(k,274)) &
                 * y(k,83) +rxt(k,295)* y(k,84) +rxt(k,296)* y(k,85) +rxt(k,297) &
                 * y(k,86) + (rxt(k,250) +rxt(k,251))* y(k,87) +rxt(k,318)* y(k,88) &
                  + (rxt(k,210) +rxt(k,211))* y(k,137) +rxt(k,172)* y(k,158) &
                  +rxt(k,169)* y(k,319) + rxt(k,170) + rxt(k,171) + het_rates(k,294)) &
                 * y(k,294)
         prod(k,314) =rxt(k,12)*y(k,137) +rxt(k,7)*y(k,158) +rxt(k,1)*y(k,319)
         loss(k,316) = (rxt(k,408)* y(k,1) +rxt(k,412)* y(k,2) +rxt(k,603)* y(k,4) &
                  +rxt(k,623)* y(k,7) +rxt(k,543)* y(k,8) +rxt(k,546)* y(k,9) &
                  +rxt(k,413)* y(k,15) +rxt(k,375)* y(k,16) +rxt(k,643)* y(k,17) &
                  +rxt(k,266)* y(k,20) +rxt(k,547)* y(k,23) +rxt(k,549)* y(k,24) &
                  +rxt(k,321)* y(k,25) +rxt(k,350)* y(k,26) +rxt(k,328)* y(k,27) &
                  +rxt(k,329)* y(k,28) +rxt(k,331)* y(k,29) +rxt(k,372)* y(k,30) &
                  +rxt(k,358)* y(k,31) +rxt(k,359)* y(k,32) +rxt(k,553)* y(k,33) &
                  +rxt(k,282)* y(k,42) +rxt(k,301)* y(k,43) +rxt(k,284)* y(k,44) &
                  +rxt(k,285)* y(k,45) +rxt(k,333)* y(k,46) +rxt(k,287)* y(k,47) &
                  +rxt(k,334)* y(k,48) +rxt(k,373)* y(k,49) +rxt(k,361)* y(k,50) &
                  +rxt(k,339)* y(k,51) +rxt(k,340)* y(k,52) +rxt(k,306)* y(k,53) &
                  +rxt(k,307)* y(k,54) +rxt(k,308)* y(k,55) +rxt(k,289)* y(k,56) &
                  + (rxt(k,236) +rxt(k,237))* y(k,60) +rxt(k,234)* y(k,61) +rxt(k,319) &
                 * y(k,63) +rxt(k,554)* y(k,67) + (rxt(k,794) +rxt(k,808))* y(k,69) &
                  +rxt(k,345)* y(k,76) +rxt(k,346)* y(k,77) +rxt(k,185)* y(k,79) &
                  +rxt(k,186)* y(k,81) +rxt(k,268)* y(k,83) +rxt(k,290)* y(k,84) &
                  +rxt(k,291)* y(k,85) +rxt(k,292)* y(k,86) +rxt(k,239)* y(k,87) &
                  +rxt(k,309)* y(k,88) +rxt(k,347)* y(k,89) +rxt(k,310)* y(k,90) &
                  +rxt(k,311)* y(k,92) +rxt(k,215)* y(k,93) +rxt(k,193)* y(k,94) &
                  +rxt(k,242)* y(k,96) +rxt(k,379)* y(k,97) +rxt(k,414)* y(k,98) &
                  +rxt(k,415)* y(k,99) +rxt(k,362)* y(k,102) +rxt(k,418)* y(k,103) &
                  +rxt(k,363)* y(k,104) +rxt(k,420)* y(k,106) +rxt(k,422)* y(k,107) &
                  +rxt(k,423)* y(k,108) +rxt(k,487)* y(k,109) +rxt(k,453)* y(k,110) &
                  +rxt(k,452)* y(k,111) +rxt(k,454)* y(k,112) +rxt(k,455)* y(k,113) &
                  +rxt(k,460)* y(k,115) +rxt(k,463)* y(k,116) +rxt(k,466)* y(k,117) &
                  +rxt(k,470)* y(k,118) +rxt(k,472)* y(k,119) +rxt(k,481)* y(k,120) &
                  +rxt(k,485)* y(k,121) +rxt(k,488)* y(k,122) + (rxt(k,489) + &
                 rxt(k,490))* y(k,123) +rxt(k,663)* y(k,125) +rxt(k,386)* y(k,126) &
                  +rxt(k,380)* y(k,127) +rxt(k,396)* y(k,129) +rxt(k,397)* y(k,130) &
                  +rxt(k,398)* y(k,131) +rxt(k,404)* y(k,132) +rxt(k,399)* y(k,133) &
                  +rxt(k,405)* y(k,134) +rxt(k,683)* y(k,135) +rxt(k,206)* y(k,136) &
                  +rxt(k,502)* y(k,139) +rxt(k,810)* y(k,143) +rxt(k,214)* y(k,148) &
                  +rxt(k,205)* y(k,149) +rxt(k,348)* y(k,150) +rxt(k,364)* y(k,151) &
                  +rxt(k,188)* y(k,157) +rxt(k,189)* y(k,158) +rxt(k,796)* y(k,161) &
                  +rxt(k,349)* y(k,163) +rxt(k,566)* y(k,166) +rxt(k,569)* y(k,167) &
                  +rxt(k,367)* y(k,170) +rxt(k,371)* y(k,171) +rxt(k,802)* y(k,172) &
                  +rxt(k,807)* y(k,174) +rxt(k,798)* y(k,175) +rxt(k,689)* y(k,200) &
                  +rxt(k,690)* y(k,201) +rxt(k,756)* y(k,202) +rxt(k,717)* y(k,203) &
                  +rxt(k,718)* y(k,204) +rxt(k,736)* y(k,205) +rxt(k,737)* y(k,206) &
                  +rxt(k,748)* y(k,207) +rxt(k,746)* y(k,208) +rxt(k,747)* y(k,209) &
                  +rxt(k,757)* y(k,210) +rxt(k,763)* y(k,212) +rxt(k,768)* y(k,213) &
                  +rxt(k,769)* y(k,214) +rxt(k,771)* y(k,216) +rxt(k,775)* y(k,217) &
                  +rxt(k,774)* y(k,218) +rxt(k,778)* y(k,220) +rxt(k,783)* y(k,221) &
                  +rxt(k,782)* y(k,222) +rxt(k,787)* y(k,223) +rxt(k,786)* y(k,224) &
                  +rxt(k,573)* y(k,227) +rxt(k,574)* y(k,228) +rxt(k,576)* y(k,229) &
                  +rxt(k,579)* y(k,230) +rxt(k,582)* y(k,231) +rxt(k,583)* y(k,232) &
                  +rxt(k,187)* y(k,258) + 2._r8*(rxt(k,190) +rxt(k,191))* y(k,295) &
                  + het_rates(k,295))* y(k,295)
         prod(k,316) = (2.000_r8*rxt(k,179)*y(k,78) +rxt(k,182)*y(k,157) + &
                 rxt(k,183)*y(k,158) +rxt(k,202)*y(k,149) +rxt(k,207)*y(k,147) + &
                 rxt(k,223)*y(k,57) +.200_r8*rxt(k,312)*y(k,259) + &
                 .490_r8*rxt(k,337)*y(k,252) +.150_r8*rxt(k,369)*y(k,298) + &
                 .590_r8*rxt(k,383)*y(k,284) +.490_r8*rxt(k,390)*y(k,286) + &
                 .200_r8*rxt(k,394)*y(k,288) +.540_r8*rxt(k,402)*y(k,289) + &
                 .650_r8*rxt(k,421)*y(k,260) +.060_r8*rxt(k,426)*y(k,261) + &
                 .060_r8*rxt(k,432)*y(k,262) +.580_r8*rxt(k,457)*y(k,269) + &
                 .520_r8*rxt(k,461)*y(k,270) +.600_r8*rxt(k,464)*y(k,271) + &
                 .500_r8*rxt(k,467)*y(k,272) +.400_r8*rxt(k,471)*y(k,273) + &
                 .240_r8*rxt(k,476)*y(k,274) +.850_r8*rxt(k,479)*y(k,275) + &
                 .860_r8*rxt(k,482)*y(k,276) +.800_r8*rxt(k,499)*y(k,293) + &
                 .400_r8*rxt(k,541)*y(k,235) +.400_r8*rxt(k,555)*y(k,254) + &
                 .400_r8*rxt(k,561)*y(k,287) +.700_r8*rxt(k,588)*y(k,237) + &
                 .350_r8*rxt(k,596)*y(k,238) +.500_r8*rxt(k,608)*y(k,240) + &
                 .100_r8*rxt(k,616)*y(k,241) +.470_r8*rxt(k,628)*y(k,245) + &
                 .030_r8*rxt(k,636)*y(k,246) +.500_r8*rxt(k,647)*y(k,281) + &
                 .100_r8*rxt(k,656)*y(k,282) +.480_r8*rxt(k,667)*y(k,290) + &
                 .100_r8*rxt(k,676)*y(k,291) +.180_r8*rxt(k,687)*y(k,299) + &
                 .180_r8*rxt(k,691)*y(k,300) +.490_r8*rxt(k,703)*y(k,302) + &
                 .380_r8*rxt(k,711)*y(k,303) +.490_r8*rxt(k,721)*y(k,304) + &
                 .150_r8*rxt(k,730)*y(k,305) +.530_r8*rxt(k,740)*y(k,306) + &
                 .490_r8*rxt(k,751)*y(k,307) +.100_r8*rxt(k,760)*y(k,308) + &
                 .100_r8*rxt(k,765)*y(k,309) +.100_r8*rxt(k,772)*y(k,310) + &
                 .100_r8*rxt(k,776)*y(k,311) +.100_r8*rxt(k,780)*y(k,312) + &
                 .100_r8*rxt(k,784)*y(k,313))*y(k,258) + (.300_r8*rxt(k,307)*y(k,54) + &
                 .500_r8*rxt(k,311)*y(k,92) +.650_r8*rxt(k,321)*y(k,25) + &
                 .500_r8*rxt(k,329)*y(k,28) +.890_r8*rxt(k,347)*y(k,89) + &
                 .700_r8*rxt(k,363)*y(k,104) +.500_r8*rxt(k,367)*y(k,170) + &
                 .430_r8*rxt(k,414)*y(k,98) +.530_r8*rxt(k,415)*y(k,99) + &
                 1.080_r8*rxt(k,418)*y(k,103) +.500_r8*rxt(k,454)*y(k,112) + &
                 .060_r8*rxt(k,460)*y(k,115) +.040_r8*rxt(k,470)*y(k,118) + &
                 .030_r8*rxt(k,472)*y(k,119) +.420_r8*rxt(k,481)*y(k,120) + &
                 .290_r8*rxt(k,485)*y(k,121) +.130_r8*rxt(k,489)*y(k,123) + &
                 .920_r8*rxt(k,490)*y(k,123))*y(k,295) + (rxt(k,184)*y(k,78) + &
                 .130_r8*rxt(k,323)*y(k,26) +.360_r8*rxt(k,354)*y(k,30) + &
                 .240_r8*rxt(k,385)*y(k,126) +.360_r8*rxt(k,403)*y(k,132) + &
                 .340_r8*rxt(k,459)*y(k,115) +.340_r8*rxt(k,469)*y(k,118) + &
                 .510_r8*rxt(k,484)*y(k,121) +.250_r8*rxt(k,486)*y(k,109) + &
                 .340_r8*rxt(k,501)*y(k,139) +.770_r8*rxt(k,602)*y(k,4) + &
                 .080_r8*rxt(k,622)*y(k,7) +.300_r8*rxt(k,642)*y(k,17) + &
                 .660_r8*rxt(k,662)*y(k,125) +.630_r8*rxt(k,682)*y(k,135) + &
                 .090_r8*rxt(k,762)*y(k,212))*y(k,158) + (rxt(k,176)*y(k,79) + &
                 rxt(k,177)*y(k,81) +rxt(k,238)*y(k,87) +rxt(k,241)*y(k,96) + &
                 rxt(k,267)*y(k,83) +rxt(k,269)*y(k,95) +rxt(k,300)*y(k,43))*y(k,157) &
                  + (.550_r8*rxt(k,509)*y(k,267) +.550_r8*rxt(k,511)*y(k,268) + &
                 .470_r8*rxt(k,525)*y(k,275) +.040_r8*rxt(k,527)*y(k,276) + &
                 .550_r8*rxt(k,530)*y(k,278) +.550_r8*rxt(k,533)*y(k,279))*y(k,147) &
                  + (rxt(k,168)*y(k,79) +2.000_r8*rxt(k,169)*y(k,319) + &
                 rxt(k,250)*y(k,87) +rxt(k,273)*y(k,83) +rxt(k,315)*y(k,55) + &
                 rxt(k,318)*y(k,88))*y(k,294) + (.550_r8*rxt(k,444)*y(k,267) + &
                 .550_r8*rxt(k,448)*y(k,268) +.550_r8*rxt(k,491)*y(k,278) + &
                 .550_r8*rxt(k,495)*y(k,279))*y(k,252) &
                  + (.280_r8*rxt(k,445)*y(k,267) +.280_r8*rxt(k,449)*y(k,268) + &
                 .280_r8*rxt(k,492)*y(k,278) +.280_r8*rxt(k,496)*y(k,279))*y(k,253) &
                  + (rxt(k,55) +rxt(k,56))*y(k,104) + (rxt(k,2) +rxt(k,277)*y(k,75)) &
                 *y(k,319) +rxt(k,20)*y(k,2) +rxt(k,21)*y(k,9) +rxt(k,27)*y(k,24) &
                  +rxt(k,28)*y(k,28) +rxt(k,29)*y(k,31) +rxt(k,30)*y(k,33) +rxt(k,36) &
                 *y(k,52) +rxt(k,37)*y(k,54) +.330_r8*rxt(k,39)*y(k,55) &
                  +1.500_r8*rxt(k,41)*y(k,68) +rxt(k,42)*y(k,74) +2.000_r8*rxt(k,4) &
                 *y(k,81) +rxt(k,45)*y(k,89) +2.000_r8*rxt(k,46)*y(k,92) +rxt(k,9) &
                 *y(k,93) +rxt(k,10)*y(k,94) +rxt(k,149)*y(k,95) +rxt(k,150)*y(k,96) &
                  +1.110_r8*rxt(k,48)*y(k,98) +1.180_r8*rxt(k,49)*y(k,99) +rxt(k,50) &
                 *y(k,100) +rxt(k,51)*y(k,101) +3.000_r8*rxt(k,54)*y(k,103) +rxt(k,61) &
                 *y(k,112) +rxt(k,62)*y(k,113) +rxt(k,63)*y(k,114) +.550_r8*rxt(k,64) &
                 *y(k,115) +.550_r8*rxt(k,67)*y(k,118) +rxt(k,69)*y(k,120) +rxt(k,70) &
                 *y(k,121) +rxt(k,71)*y(k,123) +rxt(k,75)*y(k,128) +rxt(k,77)*y(k,130) &
                  +rxt(k,81)*y(k,134) +.500_r8*rxt(k,830)*y(k,148) +rxt(k,87)*y(k,167) &
                  +rxt(k,88)*y(k,170) +rxt(k,89)*y(k,171) +rxt(k,91)*y(k,200) &
                  +rxt(k,92)*y(k,201) +rxt(k,98)*y(k,207) +rxt(k,99)*y(k,208) &
                  +rxt(k,100)*y(k,209) +rxt(k,102)*y(k,211) +rxt(k,104)*y(k,215) &
                  +rxt(k,105)*y(k,217) +rxt(k,106)*y(k,218) +rxt(k,107)*y(k,219) &
                  +rxt(k,108)*y(k,220) +rxt(k,113)*y(k,225) +rxt(k,114)*y(k,226) &
                  +rxt(k,115)*y(k,227) +rxt(k,116)*y(k,230) +rxt(k,117)*y(k,232) &
                  +rxt(k,427)*y(k,261) +rxt(k,433)*y(k,262) +rxt(k,480)*y(k,275) &
                  +rxt(k,483)*y(k,276) +.600_r8*rxt(k,529)*y(k,278) &
                  +.600_r8*rxt(k,532)*y(k,279) +rxt(k,384)*y(k,284) +rxt(k,500) &
                 *y(k,293)
         loss(k,130) = (rxt(k,565)* y(k,147) +rxt(k,564)* y(k,258) + het_rates(k,296)) &
                 * y(k,296)
         prod(k,130) = (.200_r8*rxt(k,554)*y(k,67) +.140_r8*rxt(k,566)*y(k,166) + &
                 rxt(k,569)*y(k,167))*y(k,295)
         loss(k,189) = (rxt(k,366)* y(k,147) +rxt(k,365)* y(k,258) + het_rates(k,297)) &
                 * y(k,297)
         prod(k,189) = (.500_r8*rxt(k,367)*y(k,170) +rxt(k,372)*y(k,30))*y(k,295)
         loss(k,224) = (rxt(k,370)* y(k,147) +rxt(k,368)* y(k,253) +rxt(k,369) &
                 * y(k,258) + het_rates(k,298))* y(k,298)
         prod(k,224) = (rxt(k,371)*y(k,171) +rxt(k,373)*y(k,49))*y(k,295)
         loss(k,190) = (rxt(k,688)* y(k,147) +rxt(k,687)* y(k,258) + het_rates(k,299)) &
                 * y(k,299)
         prod(k,190) =rxt(k,689)*y(k,295)*y(k,200)
         loss(k,196) = (rxt(k,692)* y(k,147) +rxt(k,691)* y(k,258) + het_rates(k,300)) &
                 * y(k,300)
         prod(k,196) =rxt(k,690)*y(k,295)*y(k,201)
         loss(k,294) = (rxt(k,696)* y(k,147) +rxt(k,697)* y(k,149) +rxt(k,693) &
                 * y(k,252) +rxt(k,694)* y(k,253) +rxt(k,695)* y(k,258) +rxt(k,698) &
                 * y(k,302) +rxt(k,699)* y(k,304) + het_rates(k,301))* y(k,301)
         prod(k,294) = (rxt(k,593)*y(k,237) +rxt(k,601)*y(k,238) + &
                 rxt(k,613)*y(k,240) +rxt(k,621)*y(k,241) +rxt(k,633)*y(k,245) + &
                 rxt(k,641)*y(k,246) +rxt(k,653)*y(k,281) +rxt(k,661)*y(k,282) + &
                 rxt(k,673)*y(k,290) +rxt(k,681)*y(k,291) +rxt(k,707)*y(k,302) + &
                 rxt(k,716)*y(k,303) +rxt(k,726)*y(k,304) +rxt(k,735)*y(k,305) + &
                 rxt(k,745)*y(k,306) +rxt(k,749)*y(k,252) +rxt(k,750)*y(k,253) + &
                 .490_r8*rxt(k,751)*y(k,258) +rxt(k,752)*y(k,147) + &
                 rxt(k,753)*y(k,149) +2.000_r8*rxt(k,754)*y(k,307))*y(k,307) &
                  + (rxt(k,98) +.290_r8*rxt(k,748)*y(k,295))*y(k,207) +rxt(k,93) &
                 *y(k,202) +.860_r8*rxt(k,771)*y(k,295)*y(k,216)
         loss(k,297) = (rxt(k,704)* y(k,147) +rxt(k,684)* y(k,148) +rxt(k,705) &
                 * y(k,149) +rxt(k,591)* y(k,237) +rxt(k,599)* y(k,238) +rxt(k,611) &
                 * y(k,240) +rxt(k,619)* y(k,241) +rxt(k,631)* y(k,245) +rxt(k,639) &
                 * y(k,246) +rxt(k,701)* y(k,252) +rxt(k,702)* y(k,253) +rxt(k,703) &
                 * y(k,258) +rxt(k,651)* y(k,281) +rxt(k,659)* y(k,282) +rxt(k,671) &
                 * y(k,290) +rxt(k,679)* y(k,291) +rxt(k,698)* y(k,301) &
                  + 2._r8*rxt(k,706)* y(k,302) +rxt(k,714)* y(k,303) +rxt(k,724) &
                 * y(k,304) +rxt(k,733)* y(k,305) +rxt(k,743)* y(k,306) +rxt(k,707) &
                 * y(k,307) + het_rates(k,302))* y(k,302)
         prod(k,297) = (rxt(k,717)*y(k,203) +.710_r8*rxt(k,746)*y(k,208) + &
                 .140_r8*rxt(k,771)*y(k,216))*y(k,295) + (.270_r8*rxt(k,602)*y(k,4) + &
                 .300_r8*rxt(k,642)*y(k,17))*y(k,158) + (rxt(k,95) +rxt(k,790)) &
                 *y(k,204) +rxt(k,708)*y(k,203)*y(k,149)
         loss(k,296) = (rxt(k,712)* y(k,147) +rxt(k,713)* y(k,149) +rxt(k,709) &
                 * y(k,252) +rxt(k,710)* y(k,253) +rxt(k,711)* y(k,258) +rxt(k,715) &
                 * y(k,304) +rxt(k,716)* y(k,307) + het_rates(k,303))* y(k,303)
         prod(k,296) = (rxt(k,591)*y(k,237) +rxt(k,599)*y(k,238) + &
                 rxt(k,611)*y(k,240) +rxt(k,619)*y(k,241) +rxt(k,631)*y(k,245) + &
                 rxt(k,639)*y(k,246) +rxt(k,651)*y(k,281) +rxt(k,659)*y(k,282) + &
                 rxt(k,671)*y(k,290) +rxt(k,679)*y(k,291) + &
                 2.000_r8*rxt(k,698)*y(k,301) +rxt(k,701)*y(k,252) + &
                 rxt(k,702)*y(k,253) +.490_r8*rxt(k,703)*y(k,258) + &
                 rxt(k,704)*y(k,147) +rxt(k,705)*y(k,149) + &
                 2.000_r8*rxt(k,706)*y(k,302) +rxt(k,707)*y(k,307) + &
                 rxt(k,724)*y(k,304) +rxt(k,733)*y(k,305) +rxt(k,743)*y(k,306)) &
                 *y(k,302) + (rxt(k,693)*y(k,252) +.500_r8*rxt(k,694)*y(k,253) + &
                 .700_r8*rxt(k,696)*y(k,147) +rxt(k,697)*y(k,149) + &
                 rxt(k,699)*y(k,304) +rxt(k,700)*y(k,307))*y(k,301) + (rxt(k,99) + &
                 .290_r8*rxt(k,746)*y(k,295))*y(k,208) +.330_r8*rxt(k,602)*y(k,158) &
                 *y(k,4) +.230_r8*rxt(k,756)*y(k,295)*y(k,202) +rxt(k,94)*y(k,203)
         loss(k,298) = (rxt(k,722)* y(k,147) +rxt(k,685)* y(k,148) +rxt(k,723) &
                 * y(k,149) +rxt(k,592)* y(k,237) +rxt(k,600)* y(k,238) +rxt(k,612) &
                 * y(k,240) +rxt(k,620)* y(k,241) +rxt(k,632)* y(k,245) +rxt(k,640) &
                 * y(k,246) +rxt(k,719)* y(k,252) +rxt(k,720)* y(k,253) +rxt(k,721) &
                 * y(k,258) +rxt(k,652)* y(k,281) +rxt(k,660)* y(k,282) +rxt(k,672) &
                 * y(k,290) +rxt(k,680)* y(k,291) +rxt(k,699)* y(k,301) +rxt(k,724) &
                 * y(k,302) +rxt(k,715)* y(k,303) + 2._r8*rxt(k,725)* y(k,304) &
                  +rxt(k,734)* y(k,305) +rxt(k,744)* y(k,306) +rxt(k,726)* y(k,307) &
                  + het_rates(k,304))* y(k,304)
         prod(k,298) = (.750_r8*rxt(k,736)*y(k,205) +.710_r8*rxt(k,747)*y(k,209) + &
                 .170_r8*rxt(k,763)*y(k,212))*y(k,295) + (rxt(k,97) +rxt(k,791)) &
                 *y(k,206) +.330_r8*rxt(k,662)*y(k,158)*y(k,125) +rxt(k,727)*y(k,205) &
                 *y(k,149)
         loss(k,278) = (rxt(k,731)* y(k,147) +rxt(k,732)* y(k,149) +rxt(k,728) &
                 * y(k,252) +rxt(k,729)* y(k,253) +rxt(k,730)* y(k,258) +rxt(k,733) &
                 * y(k,302) +rxt(k,734)* y(k,304) +rxt(k,735)* y(k,307) &
                  + het_rates(k,305))* y(k,305)
         prod(k,278) = (rxt(k,709)*y(k,252) +rxt(k,710)*y(k,253) + &
                 .380_r8*rxt(k,711)*y(k,258) +.830_r8*rxt(k,712)*y(k,147) + &
                 rxt(k,713)*y(k,149) +rxt(k,714)*y(k,302) +rxt(k,715)*y(k,304) + &
                 rxt(k,716)*y(k,307))*y(k,303)
         loss(k,295) = (rxt(k,741)* y(k,147) +rxt(k,742)* y(k,149) +rxt(k,738) &
                 * y(k,252) +rxt(k,739)* y(k,253) +rxt(k,740)* y(k,258) +rxt(k,743) &
                 * y(k,302) +rxt(k,745)* y(k,307) + het_rates(k,306))* y(k,306)
         prod(k,295) = (rxt(k,592)*y(k,237) +rxt(k,600)*y(k,238) + &
                 rxt(k,612)*y(k,240) +rxt(k,620)*y(k,241) +rxt(k,632)*y(k,245) + &
                 rxt(k,640)*y(k,246) +rxt(k,652)*y(k,281) +rxt(k,660)*y(k,282) + &
                 rxt(k,672)*y(k,290) +rxt(k,680)*y(k,291) +rxt(k,699)*y(k,301) + &
                 rxt(k,715)*y(k,303) +rxt(k,719)*y(k,252) +rxt(k,720)*y(k,253) + &
                 .490_r8*rxt(k,721)*y(k,258) +rxt(k,722)*y(k,147) + &
                 rxt(k,723)*y(k,149) +rxt(k,724)*y(k,302) + &
                 2.000_r8*rxt(k,725)*y(k,304) +rxt(k,726)*y(k,307) + &
                 2.000_r8*rxt(k,734)*y(k,305))*y(k,304) + (rxt(k,728)*y(k,252) + &
                 rxt(k,729)*y(k,253) +.150_r8*rxt(k,730)*y(k,258) + &
                 .700_r8*rxt(k,731)*y(k,147) +rxt(k,732)*y(k,149) + &
                 rxt(k,733)*y(k,302) +rxt(k,735)*y(k,307))*y(k,305) + (rxt(k,96) + &
                 .250_r8*rxt(k,736)*y(k,295))*y(k,205) + (rxt(k,100) + &
                 .290_r8*rxt(k,747)*y(k,295))*y(k,209)
         loss(k,299) = (rxt(k,752)* y(k,147) +rxt(k,686)* y(k,148) +rxt(k,753) &
                 * y(k,149) +rxt(k,593)* y(k,237) +rxt(k,601)* y(k,238) +rxt(k,613) &
                 * y(k,240) +rxt(k,621)* y(k,241) +rxt(k,633)* y(k,245) +rxt(k,641) &
                 * y(k,246) +rxt(k,749)* y(k,252) +rxt(k,750)* y(k,253) +rxt(k,751) &
                 * y(k,258) +rxt(k,653)* y(k,281) +rxt(k,661)* y(k,282) +rxt(k,673) &
                 * y(k,290) +rxt(k,681)* y(k,291) +rxt(k,700)* y(k,301) +rxt(k,707) &
                 * y(k,302) +rxt(k,716)* y(k,303) +rxt(k,726)* y(k,304) +rxt(k,735) &
                 * y(k,305) +rxt(k,745)* y(k,306) + 2._r8*rxt(k,754)* y(k,307) &
                  + het_rates(k,307))* y(k,307)
         prod(k,299) = (rxt(k,755)*y(k,149) +.770_r8*rxt(k,756)*y(k,295))*y(k,202) &
                  + (rxt(k,101) +rxt(k,792))*y(k,210) +.710_r8*rxt(k,748)*y(k,295) &
                 *y(k,207)
         loss(k,174) = (rxt(k,761)* y(k,147) +rxt(k,760)* y(k,258) + het_rates(k,308)) &
                 * y(k,308)
         prod(k,174) =.830_r8*rxt(k,763)*y(k,295)*y(k,212)
         loss(k,191) = (rxt(k,766)* y(k,147) +rxt(k,765)* y(k,258) + het_rates(k,309)) &
                 * y(k,309)
         prod(k,191) =rxt(k,768)*y(k,295)*y(k,213)
         loss(k,214) = (rxt(k,773)* y(k,147) +rxt(k,772)* y(k,258) + het_rates(k,310)) &
                 * y(k,310)
         prod(k,214) =rxt(k,774)*y(k,295)*y(k,218)
         loss(k,197) = (rxt(k,777)* y(k,147) +rxt(k,776)* y(k,258) + het_rates(k,311)) &
                 * y(k,311)
         prod(k,197) =rxt(k,778)*y(k,295)*y(k,220)
         loss(k,175) = (rxt(k,781)* y(k,147) +rxt(k,780)* y(k,258) + het_rates(k,312)) &
                 * y(k,312)
         prod(k,175) =rxt(k,782)*y(k,295)*y(k,222)
         loss(k,176) = (rxt(k,785)* y(k,147) +rxt(k,784)* y(k,258) + het_rates(k,313)) &
                 * y(k,313)
         prod(k,176) =rxt(k,786)*y(k,295)*y(k,224)
         loss(k,181) = (rxt(k,572)* y(k,147) +rxt(k,571)* y(k,258) + het_rates(k,314)) &
                 * y(k,314)
         prod(k,181) = (rxt(k,573)*y(k,227) +.650_r8*rxt(k,574)*y(k,228))*y(k,295)
         loss(k,56) = (rxt(k,880)* y(k,147) +rxt(k,879)* y(k,258) + het_rates(k,315)) &
                 * y(k,315)
         prod(k,56) =rxt(k,878)*y(k,295)*y(k,228)
         loss(k,183) = (rxt(k,578)* y(k,147) +rxt(k,577)* y(k,258) + het_rates(k,316)) &
                 * y(k,316)
         prod(k,183) = (.560_r8*rxt(k,576)*y(k,229) +rxt(k,579)*y(k,230))*y(k,295)
         loss(k,57) = (rxt(k,883)* y(k,147) +rxt(k,882)* y(k,258) + het_rates(k,317)) &
                 * y(k,317)
         prod(k,57) =rxt(k,881)*y(k,295)*y(k,229)
         loss(k,142) = (rxt(k,581)* y(k,147) +rxt(k,580)* y(k,258) + het_rates(k,318)) &
                 * y(k,318)
         prod(k,142) = (.300_r8*rxt(k,582)*y(k,231) +rxt(k,583)*y(k,232))*y(k,295)
         loss(k,317) = (rxt(k,277)* y(k,75) +rxt(k,809)* y(k,176) +rxt(k,169) &
                 * y(k,294) + rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,319)) &
                 * y(k,319)
         prod(k,317) = (rxt(k,185)*y(k,79) +rxt(k,186)*y(k,81) +rxt(k,187)*y(k,258) + &
                 rxt(k,190)*y(k,295) +rxt(k,193)*y(k,94) +rxt(k,215)*y(k,93) + &
                 rxt(k,239)*y(k,87) +rxt(k,242)*y(k,96) +rxt(k,268)*y(k,83) + &
                 rxt(k,282)*y(k,42) +rxt(k,284)*y(k,44) +rxt(k,285)*y(k,45) + &
                 rxt(k,287)*y(k,47) +rxt(k,292)*y(k,86) +rxt(k,301)*y(k,43) + &
                 rxt(k,307)*y(k,54) +rxt(k,308)*y(k,55) +rxt(k,310)*y(k,90) + &
                 rxt(k,311)*y(k,92) +rxt(k,331)*y(k,29) +rxt(k,333)*y(k,46) + &
                 rxt(k,339)*y(k,51) +rxt(k,340)*y(k,52) +rxt(k,358)*y(k,31) + &
                 rxt(k,359)*y(k,32) +rxt(k,361)*y(k,50) +rxt(k,367)*y(k,170) + &
                 rxt(k,371)*y(k,171) +rxt(k,373)*y(k,49) + &
                 .450_r8*rxt(k,386)*y(k,126) +rxt(k,775)*y(k,217) + &
                 rxt(k,779)*y(k,219) +rxt(k,810)*y(k,143))*y(k,295) &
                  + (rxt(k,885)*y(k,96) +rxt(k,891)*y(k,96) +rxt(k,892)*y(k,95) + &
                 rxt(k,896)*y(k,96) +rxt(k,897)*y(k,95))*y(k,87) + (rxt(k,812) + &
                 rxt(k,180)*y(k,78) +.300_r8*rxt(k,312)*y(k,259))*y(k,258) &
                  +.050_r8*rxt(k,39)*y(k,55) +rxt(k,153)*y(k,82)
      end do
      end subroutine imp_prod_loss
      end module mo_prod_loss
