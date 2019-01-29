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
         loss(k,1) = ( + het_rates(k,6))* y(k,6)
         prod(k,1) = 0._r8
         loss(k,2) = (rxt(k,180)* y(k,92) + rxt(k,30) + het_rates(k,7))* y(k,7)
         prod(k,2) = 0._r8
         loss(k,3) = (rxt(k,181)* y(k,92) + rxt(k,31) + het_rates(k,8))* y(k,8)
         prod(k,3) = 0._r8
         loss(k,4) = (rxt(k,207)* y(k,92) + rxt(k,32) + het_rates(k,9))* y(k,9)
         prod(k,4) = 0._r8
         loss(k,5) = (rxt(k,182)* y(k,92) + rxt(k,33) + het_rates(k,10))* y(k,10)
         prod(k,5) = 0._r8
         loss(k,6) = (rxt(k,183)* y(k,92) + rxt(k,34) + het_rates(k,11))* y(k,11)
         prod(k,6) = 0._r8
         loss(k,7) = (rxt(k,184)* y(k,92) + rxt(k,35) + het_rates(k,12))* y(k,12)
         prod(k,7) = 0._r8
         loss(k,8) = (rxt(k,185)* y(k,92) + rxt(k,36) + het_rates(k,13))* y(k,13)
         prod(k,8) = 0._r8
         loss(k,9) = (rxt(k,186)* y(k,92) + rxt(k,37) + het_rates(k,14))* y(k,14)
         prod(k,9) = 0._r8
         loss(k,10) = (rxt(k,218)* y(k,85) +rxt(k,230)* y(k,92) +rxt(k,219)* y(k,96) &
                  + rxt(k,38) + het_rates(k,15))* y(k,15)
         prod(k,10) = 0._r8
         loss(k,11) = (rxt(k,220)* y(k,85) +rxt(k,231)* y(k,92) +rxt(k,221)* y(k,96) &
                  + rxt(k,39) + het_rates(k,17))* y(k,17)
         prod(k,11) = 0._r8
         loss(k,12) = (rxt(k,222)* y(k,96) + rxt(k,40) + het_rates(k,18))* y(k,18)
         prod(k,12) = 0._r8
         loss(k,13) = (rxt(k,223)* y(k,85) +rxt(k,224)* y(k,96) + rxt(k,41) &
                  + het_rates(k,19))* y(k,19)
         prod(k,13) = 0._r8
         loss(k,14) = (rxt(k,212)* y(k,37) +rxt(k,156)* y(k,85) + (rxt(k,243) + &
                 rxt(k,244) +rxt(k,245))* y(k,92) +rxt(k,241)* y(k,96) + rxt(k,23) &
                  + rxt(k,24) + het_rates(k,22))* y(k,22)
         prod(k,14) = 0._r8
         loss(k,15) = (rxt(k,225)* y(k,85) +rxt(k,208)* y(k,92) +rxt(k,226)* y(k,96) &
                  + rxt(k,42) + het_rates(k,23))* y(k,23)
         prod(k,15) = 0._r8
         loss(k,16) = ( + het_rates(k,28))* y(k,28)
         prod(k,16) = 0._r8
         loss(k,17) = (rxt(k,300)* y(k,97) + rxt(k,25) + rxt(k,61) + het_rates(k,30)) &
                 * y(k,30)
         prod(k,17) =.440_r8*rxt(k,24)*y(k,22)
         loss(k,18) = (rxt(k,209)* y(k,92) + rxt(k,50) + het_rates(k,40))* y(k,40)
         prod(k,18) = 0._r8
         loss(k,19) = (rxt(k,232)* y(k,92) +rxt(k,227)* y(k,96) + rxt(k,52) &
                  + het_rates(k,44))* y(k,44)
         prod(k,19) = 0._r8
         loss(k,20) = (rxt(k,233)* y(k,92) +rxt(k,228)* y(k,96) + rxt(k,53) &
                  + het_rates(k,45))* y(k,45)
         prod(k,20) = 0._r8
         loss(k,21) = (rxt(k,234)* y(k,92) +rxt(k,229)* y(k,96) + rxt(k,54) &
                  + het_rates(k,46))* y(k,46)
         prod(k,21) = 0._r8
         loss(k,22) = ((rxt(k,147) +rxt(k,148))* y(k,92) + rxt(k,12) &
                  + het_rates(k,54))* y(k,54)
         prod(k,22) = 0._r8
         loss(k,23) = ( + rxt(k,60) + het_rates(k,74))* y(k,74)
         prod(k,23) = 0._r8
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
         loss(k,1) = ( + het_rates(k,1))* y(k,1)
         prod(k,1) = 0._r8
         loss(k,2) = ( + het_rates(k,2))* y(k,2)
         prod(k,2) = 0._r8
         loss(k,29) = ( + rxt(k,26) + het_rates(k,3))* y(k,3)
         prod(k,29) = (rxt(k,276)*y(k,51) +rxt(k,281)*y(k,51))*y(k,47) &
                  +rxt(k,197)*y(k,26)*y(k,4)
         loss(k,71) = (2._r8*rxt(k,194)* y(k,4) + (rxt(k,195) +rxt(k,196) +rxt(k,197)) &
                 * y(k,26) +rxt(k,199)* y(k,59) +rxt(k,200)* y(k,60) +rxt(k,202) &
                 * y(k,66) +rxt(k,251)* y(k,75) +rxt(k,198)* y(k,87) +rxt(k,203) &
                 * y(k,96) + rxt(k,27) + het_rates(k,4))* y(k,4)
         prod(k,71) = (rxt(k,28) +rxt(k,201)*y(k,66))*y(k,5) +rxt(k,211)*y(k,92) &
                 *y(k,43) +rxt(k,206)*y(k,66)*y(k,51) +rxt(k,193)*y(k,84)*y(k,68)
         loss(k,40) = (rxt(k,201)* y(k,66) + rxt(k,28) + rxt(k,29) + rxt(k,270) &
                  + rxt(k,273) + rxt(k,278) + het_rates(k,5))* y(k,5)
         prod(k,40) =rxt(k,200)*y(k,60)*y(k,4)
         loss(k,68) = (rxt(k,235)* y(k,61) +rxt(k,236)* y(k,66) +rxt(k,191)* y(k,84) &
                  +rxt(k,155)* y(k,85) +rxt(k,237)* y(k,96) + rxt(k,20) + rxt(k,21) &
                  + het_rates(k,16))* y(k,16)
         prod(k,68) = (rxt(k,162)*y(k,26) +rxt(k,239)*y(k,59))*y(k,20) + (rxt(k,22) + &
                 .300_r8*rxt(k,240)*y(k,96))*y(k,21) + (rxt(k,244)*y(k,92) + &
                 rxt(k,245)*y(k,92))*y(k,22)
         loss(k,62) = (rxt(k,162)* y(k,26) +rxt(k,239)* y(k,59) +rxt(k,238)* y(k,87) &
                  + het_rates(k,20))* y(k,20)
         prod(k,62) = (rxt(k,156)*y(k,85) +rxt(k,212)*y(k,37) +rxt(k,241)*y(k,96) + &
                 rxt(k,243)*y(k,92))*y(k,22) +.700_r8*rxt(k,240)*y(k,96)*y(k,21)
         loss(k,35) = (rxt(k,240)* y(k,96) + rxt(k,22) + het_rates(k,21))* y(k,21)
         prod(k,35) =rxt(k,238)*y(k,87)*y(k,20)
         loss(k,26) = ( + rxt(k,43) + het_rates(k,24))* y(k,24)
         prod(k,26) = (rxt(k,269)*y(k,52) +rxt(k,274)*y(k,27) +rxt(k,275)*y(k,52) + &
                 rxt(k,279)*y(k,27) +rxt(k,280)*y(k,52) +rxt(k,284)*y(k,27))*y(k,47) &
                  +rxt(k,164)*y(k,26)*y(k,26) +rxt(k,168)*y(k,85)*y(k,27)
         loss(k,21) = ( + rxt(k,44) + rxt(k,190) + het_rates(k,25))* y(k,25)
         prod(k,21) =rxt(k,189)*y(k,26)*y(k,26)
         loss(k,66) = ((rxt(k,195) +rxt(k,196) +rxt(k,197))* y(k,4) +rxt(k,162) &
                 * y(k,20) + 2._r8*(rxt(k,163) +rxt(k,164) +rxt(k,165) +rxt(k,189)) &
                 * y(k,26) +rxt(k,167)* y(k,59) +rxt(k,169)* y(k,60) +rxt(k,172) &
                 * y(k,66) +rxt(k,252)* y(k,75) +rxt(k,166)* y(k,87) + (rxt(k,173) + &
                 rxt(k,174))* y(k,96) + rxt(k,45) + het_rates(k,26))* y(k,26)
         prod(k,66) = (rxt(k,177)*y(k,85) +rxt(k,178)*y(k,66) +rxt(k,179)*y(k,96)) &
                 *y(k,52) + (rxt(k,47) +rxt(k,170)*y(k,66))*y(k,27) + (rxt(k,59) + &
                 rxt(k,257)*y(k,75))*y(k,69) + (rxt(k,160)*y(k,87) + &
                 rxt(k,161)*y(k,68))*y(k,85) +2.000_r8*rxt(k,190)*y(k,25) &
                  +rxt(k,188)*y(k,92)*y(k,47)
         loss(k,56) = ((rxt(k,274) +rxt(k,279) +rxt(k,284))* y(k,47) +rxt(k,170) &
                 * y(k,66) +rxt(k,168)* y(k,85) +rxt(k,171)* y(k,96) + rxt(k,46) &
                  + rxt(k,47) + rxt(k,272) + rxt(k,277) + rxt(k,283) &
                  + het_rates(k,27))* y(k,27)
         prod(k,56) =rxt(k,169)*y(k,60)*y(k,26)
         loss(k,39) = ((rxt(k,242) +rxt(k,246))* y(k,96) + het_rates(k,29))* y(k,29)
         prod(k,39) = (rxt(k,20) +rxt(k,21) +rxt(k,155)*y(k,85) +rxt(k,191)*y(k,84) + &
                 rxt(k,235)*y(k,61) +rxt(k,236)*y(k,66) +rxt(k,237)*y(k,96))*y(k,16) &
                  + (rxt(k,88) +rxt(k,247)*y(k,66) +rxt(k,248)*y(k,96))*y(k,70) &
                  +rxt(k,223)*y(k,85)*y(k,19) +rxt(k,300)*y(k,97)*y(k,30)
         loss(k,23) = (rxt(k,216)* y(k,92) + rxt(k,48) + het_rates(k,31))* y(k,31)
         prod(k,23) = (rxt(k,181)*y(k,8) +rxt(k,183)*y(k,11) + &
                 2.000_r8*rxt(k,184)*y(k,12) +2.000_r8*rxt(k,185)*y(k,13) + &
                 rxt(k,186)*y(k,14) +rxt(k,207)*y(k,9) +2.000_r8*rxt(k,209)*y(k,40) + &
                 rxt(k,233)*y(k,45) +rxt(k,234)*y(k,46))*y(k,92) &
                  + (rxt(k,228)*y(k,45) +rxt(k,229)*y(k,46))*y(k,96)
         loss(k,27) = (rxt(k,217)* y(k,92) + rxt(k,49) + het_rates(k,32))* y(k,32)
         prod(k,27) = (rxt(k,182)*y(k,10) +rxt(k,183)*y(k,11) +rxt(k,232)*y(k,44)) &
                 *y(k,92) +rxt(k,227)*y(k,96)*y(k,44)
         loss(k,32) = (rxt(k,261)* y(k,61) + (rxt(k,262) +rxt(k,263))* y(k,96) &
                  + het_rates(k,33))* y(k,33)
         prod(k,32) = 0._r8
         loss(k,3) = ( + het_rates(k,34))* y(k,34)
         prod(k,3) = 0._r8
         loss(k,4) = ( + het_rates(k,35))* y(k,35)
         prod(k,4) = 0._r8
         loss(k,5) = ( + het_rates(k,36))* y(k,36)
         prod(k,5) = 0._r8
         loss(k,46) = (rxt(k,212)* y(k,22) +rxt(k,213)* y(k,39) +rxt(k,215)* y(k,49) &
                  +rxt(k,214)* y(k,100) + het_rates(k,37))* y(k,37)
         prod(k,46) = (rxt(k,185)*y(k,13) +rxt(k,207)*y(k,9) + &
                 2.000_r8*rxt(k,216)*y(k,31) +rxt(k,217)*y(k,32))*y(k,92) &
                  +2.000_r8*rxt(k,48)*y(k,31) +rxt(k,49)*y(k,32) +rxt(k,56)*y(k,48)
         loss(k,59) = (rxt(k,116)* y(k,67) +rxt(k,119)* y(k,68) + (rxt(k,113) + &
                 rxt(k,114) +rxt(k,115))* y(k,87) + het_rates(k,38))* y(k,38)
         prod(k,59) = (rxt(k,120)*y(k,39) +rxt(k,123)*y(k,66) +rxt(k,143)*y(k,53) + &
                 rxt(k,237)*y(k,16) +rxt(k,246)*y(k,29) +rxt(k,248)*y(k,70) + &
                 rxt(k,253)*y(k,73) +rxt(k,258)*y(k,75))*y(k,96) &
                  + (rxt(k,94)*y(k,92) +rxt(k,111)*y(k,66) +rxt(k,157)*y(k,85) + &
                 rxt(k,213)*y(k,37))*y(k,39) + (rxt(k,244)*y(k,22) + &
                 rxt(k,188)*y(k,47) +rxt(k,211)*y(k,43))*y(k,92) &
                  + (2.000_r8*rxt(k,2) +rxt(k,3))*y(k,100) +2.000_r8*rxt(k,20)*y(k,16) &
                  +rxt(k,22)*y(k,21) +rxt(k,51)*y(k,43) +rxt(k,55)*y(k,47) +rxt(k,56) &
                 *y(k,48)
         loss(k,70) = (rxt(k,213)* y(k,37) +rxt(k,111)* y(k,66) +rxt(k,157)* y(k,85) &
                  +rxt(k,94)* y(k,92) +rxt(k,120)* y(k,96) + het_rates(k,39))* y(k,39)
         prod(k,70) =rxt(k,21)*y(k,16) +rxt(k,245)*y(k,92)*y(k,22) +rxt(k,113)*y(k,87) &
                 *y(k,38) +rxt(k,1)*y(k,100)
         loss(k,41) = (rxt(k,112)* y(k,66) +rxt(k,158)* y(k,85) +rxt(k,121)* y(k,96) &
                  + rxt(k,4) + het_rates(k,41))* y(k,41)
         prod(k,41) = (.500_r8*rxt(k,264) +rxt(k,127)*y(k,87))*y(k,87) &
                  +rxt(k,126)*y(k,96)*y(k,96)
         loss(k,22) = ( + rxt(k,87) + het_rates(k,42))* y(k,42)
         prod(k,22) =rxt(k,260)*y(k,100)*y(k,77)
         loss(k,49) = (rxt(k,204)* y(k,66) + (rxt(k,210) +rxt(k,211))* y(k,92) &
                  +rxt(k,205)* y(k,96) + rxt(k,51) + het_rates(k,43))* y(k,43)
         prod(k,49) = (rxt(k,191)*y(k,16) +rxt(k,192)*y(k,87))*y(k,84)
         loss(k,61) = ((rxt(k,274) +rxt(k,279) +rxt(k,284))* y(k,27) + (rxt(k,276) + &
                 rxt(k,281))* y(k,51) + (rxt(k,269) +rxt(k,275) +rxt(k,280))* y(k,52) &
                  +rxt(k,175)* y(k,66) + (rxt(k,187) +rxt(k,188))* y(k,92) +rxt(k,176) &
                 * y(k,96) + rxt(k,55) + het_rates(k,47))* y(k,47)
         prod(k,61) = (rxt(k,156)*y(k,22) +rxt(k,218)*y(k,15) +rxt(k,220)*y(k,17) + &
                 2.000_r8*rxt(k,223)*y(k,19) +rxt(k,225)*y(k,23) +rxt(k,155)*y(k,16) + &
                 rxt(k,157)*y(k,39) +rxt(k,158)*y(k,41) +rxt(k,159)*y(k,87) + &
                 rxt(k,177)*y(k,52))*y(k,85) +rxt(k,174)*y(k,96)*y(k,26)
         loss(k,28) = ( + rxt(k,56) + het_rates(k,48))* y(k,48)
         prod(k,28) = (rxt(k,212)*y(k,22) +rxt(k,213)*y(k,39) +rxt(k,214)*y(k,100) + &
                 rxt(k,215)*y(k,49))*y(k,37)
         loss(k,57) = (rxt(k,215)* y(k,37) +rxt(k,152)* y(k,96) + rxt(k,9) &
                  + het_rates(k,49))* y(k,49)
         prod(k,57) = (rxt(k,272) +rxt(k,277) +rxt(k,283) +rxt(k,274)*y(k,47) + &
                 rxt(k,279)*y(k,47) +rxt(k,284)*y(k,47))*y(k,27) &
                  + (2.000_r8*rxt(k,265) +2.000_r8*rxt(k,268) +2.000_r8*rxt(k,271) + &
                 2.000_r8*rxt(k,282))*y(k,55) + (rxt(k,270) +rxt(k,273) +rxt(k,278)) &
                 *y(k,5) + (rxt(k,267) +rxt(k,235)*y(k,16) +rxt(k,261)*y(k,33)) &
                 *y(k,61) + (.500_r8*rxt(k,266) +rxt(k,151)*y(k,96))*y(k,60)
         loss(k,36) = (rxt(k,128)* y(k,96) + rxt(k,10) + rxt(k,11) + rxt(k,153) &
                  + het_rates(k,50))* y(k,50)
         prod(k,36) =rxt(k,149)*y(k,87)*y(k,60)
         loss(k,45) = ((rxt(k,276) +rxt(k,281))* y(k,47) +rxt(k,206)* y(k,66) &
                  + rxt(k,57) + het_rates(k,51))* y(k,51)
         prod(k,45) = (rxt(k,270) +rxt(k,273) +rxt(k,278))*y(k,5) +rxt(k,198)*y(k,87) &
                 *y(k,4)
         loss(k,47) = ((rxt(k,269) +rxt(k,275) +rxt(k,280))* y(k,47) +rxt(k,178) &
                 * y(k,66) +rxt(k,177)* y(k,85) +rxt(k,179)* y(k,96) + rxt(k,58) &
                  + het_rates(k,52))* y(k,52)
         prod(k,47) = (rxt(k,272) +rxt(k,277) +rxt(k,283) +rxt(k,171)*y(k,96))*y(k,27) &
                  +rxt(k,166)*y(k,87)*y(k,26)
         loss(k,54) = (rxt(k,131)* y(k,59) + (rxt(k,132) +rxt(k,133) +rxt(k,134)) &
                 * y(k,60) +rxt(k,135)* y(k,67) +rxt(k,297)* y(k,95) +rxt(k,143) &
                 * y(k,96) + rxt(k,62) + het_rates(k,53))* y(k,53)
         prod(k,54) = (rxt(k,129)*y(k,88) +rxt(k,294)*y(k,91))*y(k,66) &
                  + (.200_r8*rxt(k,288)*y(k,90) +1.100_r8*rxt(k,290)*y(k,89))*y(k,86) &
                  +rxt(k,15)*y(k,59) +rxt(k,295)*y(k,91)*y(k,67) +rxt(k,301)*y(k,97)
         loss(k,31) = ( + rxt(k,13) + rxt(k,14) + rxt(k,154) + rxt(k,265) + rxt(k,268) &
                  + rxt(k,271) + rxt(k,282) + het_rates(k,55))* y(k,55)
         prod(k,31) =rxt(k,150)*y(k,61)*y(k,60)
         loss(k,6) = ( + het_rates(k,56))* y(k,56)
         prod(k,6) = 0._r8
         loss(k,7) = ( + het_rates(k,57))* y(k,57)
         prod(k,7) = 0._r8
         loss(k,8) = ( + het_rates(k,58))* y(k,58)
         prod(k,8) = 0._r8
         loss(k,76) = (rxt(k,199)* y(k,4) +rxt(k,239)* y(k,20) +rxt(k,167)* y(k,26) &
                  +rxt(k,131)* y(k,53) +rxt(k,140)* y(k,61) +rxt(k,146)* y(k,66) &
                  +rxt(k,145)* y(k,68) +rxt(k,144)* y(k,87) +rxt(k,299)* y(k,95) &
                  + rxt(k,15) + rxt(k,16) + het_rates(k,59))* y(k,59)
         prod(k,76) = (rxt(k,17) +.500_r8*rxt(k,266) +2.000_r8*rxt(k,133)*y(k,53) + &
                 rxt(k,136)*y(k,66) +rxt(k,254)*y(k,75))*y(k,60) &
                  + (rxt(k,135)*y(k,67) +rxt(k,143)*y(k,96))*y(k,53) &
                  +2.000_r8*rxt(k,147)*y(k,92)*y(k,54) +rxt(k,14)*y(k,55) +rxt(k,19) &
                 *y(k,61) +rxt(k,130)*y(k,88)*y(k,67) +rxt(k,298)*y(k,95) +rxt(k,311) &
                 *y(k,99)
         loss(k,75) = (rxt(k,200)* y(k,4) +rxt(k,169)* y(k,26) + (rxt(k,132) + &
                 rxt(k,133) +rxt(k,134))* y(k,53) +rxt(k,150)* y(k,61) + (rxt(k,136) + &
                 rxt(k,138))* y(k,66) +rxt(k,137)* y(k,68) +rxt(k,254)* y(k,75) &
                  +rxt(k,149)* y(k,87) +rxt(k,151)* y(k,96) + rxt(k,17) + rxt(k,266) &
                  + het_rates(k,60))* y(k,60)
         prod(k,75) = (2.000_r8*rxt(k,140)*y(k,61) +rxt(k,144)*y(k,87) + &
                 rxt(k,145)*y(k,68) +rxt(k,146)*y(k,66) +rxt(k,167)*y(k,26) + &
                 rxt(k,199)*y(k,4) +rxt(k,239)*y(k,20))*y(k,59) + (rxt(k,18) + &
                 rxt(k,139)*y(k,87) +rxt(k,141)*y(k,66) +rxt(k,142)*y(k,96))*y(k,61) &
                  + (rxt(k,11) +rxt(k,153) +rxt(k,128)*y(k,96))*y(k,50) + (rxt(k,13) + &
                 rxt(k,154))*y(k,55) +rxt(k,28)*y(k,5) +rxt(k,47)*y(k,27) +rxt(k,9) &
                 *y(k,49)
         loss(k,63) = (rxt(k,235)* y(k,16) +rxt(k,261)* y(k,33) +rxt(k,140)* y(k,59) &
                  +rxt(k,150)* y(k,60) +rxt(k,141)* y(k,66) +rxt(k,139)* y(k,87) &
                  +rxt(k,142)* y(k,96) + rxt(k,18) + rxt(k,19) + rxt(k,267) &
                  + het_rates(k,61))* y(k,61)
         prod(k,63) = (rxt(k,46) +rxt(k,168)*y(k,85) +rxt(k,170)*y(k,66) + &
                 rxt(k,171)*y(k,96))*y(k,27) + (rxt(k,13) +rxt(k,14) +rxt(k,154)) &
                 *y(k,55) + (rxt(k,29) +rxt(k,201)*y(k,66))*y(k,5) &
                  + (rxt(k,152)*y(k,96) +rxt(k,215)*y(k,37))*y(k,49) &
                  + (rxt(k,137)*y(k,68) +rxt(k,138)*y(k,66))*y(k,60) +rxt(k,10) &
                 *y(k,50)
         loss(k,9) = ( + het_rates(k,62))* y(k,62)
         prod(k,9) = 0._r8
         loss(k,10) = ( + het_rates(k,63))* y(k,63)
         prod(k,10) = 0._r8
         loss(k,11) = ( + het_rates(k,64))* y(k,64)
         prod(k,11) = 0._r8
         loss(k,12) = ( + het_rates(k,65))* y(k,65)
         prod(k,12) = 0._r8
         loss(k,67) = (rxt(k,202)* y(k,4) +rxt(k,201)* y(k,5) +rxt(k,236)* y(k,16) &
                  +rxt(k,172)* y(k,26) +rxt(k,170)* y(k,27) +rxt(k,111)* y(k,39) &
                  +rxt(k,112)* y(k,41) +rxt(k,204)* y(k,43) +rxt(k,175)* y(k,47) &
                  +rxt(k,206)* y(k,51) +rxt(k,178)* y(k,52) +rxt(k,146)* y(k,59) &
                  + (rxt(k,136) +rxt(k,138))* y(k,60) +rxt(k,141)* y(k,61) &
                  + 2._r8*rxt(k,109)* y(k,66) +rxt(k,110)* y(k,67) +rxt(k,108) &
                 * y(k,68) +rxt(k,247)* y(k,70) +rxt(k,117)* y(k,87) + (rxt(k,292) + &
                 rxt(k,293))* y(k,89) +rxt(k,294)* y(k,91) +rxt(k,123)* y(k,96) &
                  + rxt(k,71) + rxt(k,72) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                  + rxt(k,76) + het_rates(k,66))* y(k,66)
         prod(k,67) = (2.000_r8*rxt(k,5) +rxt(k,6) +rxt(k,77) +rxt(k,79) +rxt(k,81) + &
                 2.000_r8*rxt(k,82) +2.000_r8*rxt(k,83) +rxt(k,84) +rxt(k,85) + &
                 rxt(k,86) +rxt(k,97)*y(k,92) +rxt(k,98)*y(k,92) +rxt(k,135)*y(k,53) + &
                 rxt(k,249)*y(k,73) +rxt(k,255)*y(k,75) +rxt(k,296)*y(k,91) + &
                 rxt(k,303)*y(k,97) +rxt(k,307)*y(k,98))*y(k,67) &
                  + (rxt(k,131)*y(k,59) +rxt(k,132)*y(k,60) +rxt(k,297)*y(k,95)) &
                 *y(k,53) + (rxt(k,288)*y(k,90) +1.150_r8*rxt(k,289)*y(k,95))*y(k,86) &
                  +rxt(k,27)*y(k,4) +rxt(k,45)*y(k,26) +rxt(k,115)*y(k,87)*y(k,38) &
                  +rxt(k,14)*y(k,55) +rxt(k,15)*y(k,59) +rxt(k,17)*y(k,60) +rxt(k,18) &
                 *y(k,61) +rxt(k,8)*y(k,68) +rxt(k,59)*y(k,69) +rxt(k,89)*y(k,75) &
                  +rxt(k,90)*y(k,76) +rxt(k,91)*y(k,77) +rxt(k,302)*y(k,97)*y(k,88) &
                  +rxt(k,96)*y(k,92) +rxt(k,125)*y(k,96)*y(k,96) +rxt(k,305)*y(k,98) &
                  +rxt(k,310)*y(k,99) +rxt(k,2)*y(k,100)
         loss(k,72) = (rxt(k,116)* y(k,38) +rxt(k,135)* y(k,53) +rxt(k,110)* y(k,66) &
                  +rxt(k,249)* y(k,73) +rxt(k,255)* y(k,75) +rxt(k,130)* y(k,88) &
                  +rxt(k,291)* y(k,89) + (rxt(k,295) +rxt(k,296))* y(k,91) +rxt(k,97) &
                 * y(k,92) +rxt(k,102)* y(k,93) +rxt(k,303)* y(k,97) +rxt(k,307) &
                 * y(k,98) + rxt(k,5) + rxt(k,6) + rxt(k,77) + rxt(k,78) + rxt(k,79) &
                  + rxt(k,80) + rxt(k,81) + rxt(k,82) + rxt(k,83) + rxt(k,84) &
                  + rxt(k,85) + rxt(k,86) + het_rates(k,67))* y(k,67)
         prod(k,72) = (rxt(k,8) +2.000_r8*rxt(k,99)*y(k,92) + &
                 2.000_r8*rxt(k,108)*y(k,66) +2.000_r8*rxt(k,118)*y(k,87) + &
                 rxt(k,119)*y(k,38) +rxt(k,124)*y(k,96) +rxt(k,137)*y(k,60) + &
                 rxt(k,145)*y(k,59) +rxt(k,161)*y(k,85) +rxt(k,193)*y(k,84) + &
                 rxt(k,250)*y(k,73) +rxt(k,256)*y(k,75))*y(k,68) &
                  + (rxt(k,113)*y(k,38) +rxt(k,117)*y(k,66) +rxt(k,122)*y(k,96) + &
                 rxt(k,127)*y(k,87) +rxt(k,139)*y(k,61) +rxt(k,159)*y(k,85) + &
                 rxt(k,166)*y(k,26) +rxt(k,192)*y(k,84) +rxt(k,198)*y(k,4) + &
                 rxt(k,238)*y(k,20))*y(k,87) + (rxt(k,101)*y(k,93) + &
                 rxt(k,109)*y(k,66) +rxt(k,123)*y(k,96) +rxt(k,136)*y(k,60) + &
                 rxt(k,141)*y(k,61) +rxt(k,172)*y(k,26) +rxt(k,202)*y(k,4))*y(k,66) &
                  + (rxt(k,163)*y(k,26) +rxt(k,164)*y(k,26) +rxt(k,174)*y(k,96) + &
                 rxt(k,196)*y(k,4) +rxt(k,197)*y(k,4))*y(k,26) + (rxt(k,92) + &
                 rxt(k,100) +2.000_r8*rxt(k,102)*y(k,67))*y(k,93) +rxt(k,194)*y(k,4) &
                 *y(k,4) +rxt(k,128)*y(k,96)*y(k,50) +rxt(k,134)*y(k,60)*y(k,53) &
                  +rxt(k,148)*y(k,92)*y(k,54) +rxt(k,299)*y(k,95)*y(k,59) +rxt(k,19) &
                 *y(k,61) +rxt(k,93)*y(k,94)
         loss(k,74) = (rxt(k,119)* y(k,38) +rxt(k,145)* y(k,59) +rxt(k,137)* y(k,60) &
                  +rxt(k,108)* y(k,66) +rxt(k,250)* y(k,73) +rxt(k,256)* y(k,75) &
                  +rxt(k,193)* y(k,84) +rxt(k,161)* y(k,85) +rxt(k,118)* y(k,87) &
                  +rxt(k,99)* y(k,92) +rxt(k,124)* y(k,96) + rxt(k,7) + rxt(k,8) &
                  + het_rates(k,68))* y(k,68)
         prod(k,74) =rxt(k,110)*y(k,67)*y(k,66)
         loss(k,33) = (rxt(k,257)* y(k,75) + rxt(k,59) + het_rates(k,69))* y(k,69)
         prod(k,33) = (rxt(k,165)*y(k,26) +rxt(k,195)*y(k,4))*y(k,26)
         loss(k,34) = (rxt(k,247)* y(k,66) +rxt(k,248)* y(k,96) + rxt(k,88) &
                  + het_rates(k,70))* y(k,70)
         prod(k,34) = 0._r8
         loss(k,13) = ( + het_rates(k,71))* y(k,71)
         prod(k,13) = 0._r8
         loss(k,14) = ( + het_rates(k,72))* y(k,72)
         prod(k,14) = 0._r8
         loss(k,43) = (rxt(k,249)* y(k,67) +rxt(k,250)* y(k,68) +rxt(k,253)* y(k,96) &
                  + het_rates(k,73))* y(k,73)
         prod(k,43) =rxt(k,88)*y(k,70) +rxt(k,89)*y(k,75)
         loss(k,58) = (rxt(k,251)* y(k,4) +rxt(k,252)* y(k,26) +rxt(k,254)* y(k,60) &
                  +rxt(k,255)* y(k,67) +rxt(k,256)* y(k,68) +rxt(k,257)* y(k,69) &
                  +rxt(k,258)* y(k,96) + rxt(k,89) + het_rates(k,75))* y(k,75)
         prod(k,58) = (rxt(k,249)*y(k,67) +rxt(k,250)*y(k,68) +rxt(k,253)*y(k,96)) &
                 *y(k,73) +rxt(k,247)*y(k,70)*y(k,66) +rxt(k,90)*y(k,76)
         loss(k,55) = (rxt(k,259)* y(k,96) + rxt(k,90) + het_rates(k,76))* y(k,76)
         prod(k,55) = (rxt(k,251)*y(k,4) +rxt(k,252)*y(k,26) +rxt(k,254)*y(k,60) + &
                 rxt(k,255)*y(k,67) +rxt(k,256)*y(k,68) +rxt(k,257)*y(k,69) + &
                 rxt(k,258)*y(k,96))*y(k,75) + (rxt(k,261)*y(k,61) + &
                 rxt(k,262)*y(k,96) +.500_r8*rxt(k,263)*y(k,96))*y(k,33) &
                  +rxt(k,248)*y(k,96)*y(k,70) +rxt(k,91)*y(k,77)
         loss(k,30) = (rxt(k,260)* y(k,100) + rxt(k,91) + het_rates(k,77))* y(k,77)
         prod(k,30) =rxt(k,87)*y(k,42) +rxt(k,259)*y(k,96)*y(k,76)
         loss(k,15) = ( + het_rates(k,78))* y(k,78)
         prod(k,15) = 0._r8
         loss(k,16) = ( + het_rates(k,79))* y(k,79)
         prod(k,16) = 0._r8
         loss(k,17) = ( + het_rates(k,80))* y(k,80)
         prod(k,17) = 0._r8
         loss(k,18) = ( + het_rates(k,81))* y(k,81)
         prod(k,18) = 0._r8
         loss(k,19) = ( + het_rates(k,82))* y(k,82)
         prod(k,19) = 0._r8
         loss(k,20) = ( + het_rates(k,83))* y(k,83)
         prod(k,20) = 0._r8
         loss(k,60) = (rxt(k,191)* y(k,16) +rxt(k,193)* y(k,68) +rxt(k,192)* y(k,87) &
                  + het_rates(k,84))* y(k,84)
         prod(k,60) = (rxt(k,27) +2.000_r8*rxt(k,194)*y(k,4) +rxt(k,195)*y(k,26) + &
                 rxt(k,196)*y(k,26) +rxt(k,199)*y(k,59) +rxt(k,202)*y(k,66) + &
                 rxt(k,203)*y(k,96) +rxt(k,251)*y(k,75))*y(k,4) + (rxt(k,181)*y(k,8) + &
                 rxt(k,207)*y(k,9) +3.000_r8*rxt(k,208)*y(k,23) + &
                 2.000_r8*rxt(k,209)*y(k,40) +2.000_r8*rxt(k,230)*y(k,15) + &
                 rxt(k,231)*y(k,17) +rxt(k,210)*y(k,43))*y(k,92) &
                  + (2.000_r8*rxt(k,219)*y(k,15) +rxt(k,221)*y(k,17) + &
                 3.000_r8*rxt(k,226)*y(k,23) +rxt(k,205)*y(k,43))*y(k,96) &
                  + (2.000_r8*rxt(k,218)*y(k,15) +rxt(k,220)*y(k,17) + &
                 3.000_r8*rxt(k,225)*y(k,23))*y(k,85) + (rxt(k,51) + &
                 rxt(k,204)*y(k,66))*y(k,43) +rxt(k,26)*y(k,3) +rxt(k,29)*y(k,5) &
                  +rxt(k,57)*y(k,51)
         loss(k,65) = (rxt(k,218)* y(k,15) +rxt(k,155)* y(k,16) +rxt(k,220)* y(k,17) &
                  +rxt(k,223)* y(k,19) +rxt(k,156)* y(k,22) +rxt(k,225)* y(k,23) &
                  +rxt(k,168)* y(k,27) +rxt(k,157)* y(k,39) +rxt(k,158)* y(k,41) &
                  +rxt(k,177)* y(k,52) +rxt(k,161)* y(k,68) + (rxt(k,159) +rxt(k,160)) &
                 * y(k,87) + het_rates(k,85))* y(k,85)
         prod(k,65) = (4.000_r8*rxt(k,180)*y(k,7) +rxt(k,181)*y(k,8) + &
                 2.000_r8*rxt(k,182)*y(k,10) +2.000_r8*rxt(k,183)*y(k,11) + &
                 2.000_r8*rxt(k,184)*y(k,12) +rxt(k,185)*y(k,13) + &
                 2.000_r8*rxt(k,186)*y(k,14) +rxt(k,232)*y(k,44) +rxt(k,233)*y(k,45) + &
                 rxt(k,234)*y(k,46) +rxt(k,187)*y(k,47) +rxt(k,217)*y(k,32))*y(k,92) &
                  + (rxt(k,45) +rxt(k,162)*y(k,20) +2.000_r8*rxt(k,163)*y(k,26) + &
                 rxt(k,165)*y(k,26) +rxt(k,167)*y(k,59) +rxt(k,172)*y(k,66) + &
                 rxt(k,173)*y(k,96) +rxt(k,196)*y(k,4) +rxt(k,252)*y(k,75))*y(k,26) &
                  + (3.000_r8*rxt(k,222)*y(k,18) +rxt(k,224)*y(k,19) + &
                 rxt(k,227)*y(k,44) +rxt(k,228)*y(k,45) +rxt(k,229)*y(k,46) + &
                 rxt(k,176)*y(k,47))*y(k,96) + (rxt(k,55) +rxt(k,175)*y(k,66))*y(k,47) &
                  +rxt(k,26)*y(k,3) +2.000_r8*rxt(k,43)*y(k,24) +2.000_r8*rxt(k,44) &
                 *y(k,25) +rxt(k,46)*y(k,27) +rxt(k,49)*y(k,32) +rxt(k,58)*y(k,52)
         loss(k,52) = (rxt(k,290)* y(k,89) +rxt(k,288)* y(k,90) +rxt(k,289)* y(k,95) &
                  + het_rates(k,86))* y(k,86)
         prod(k,52) = (rxt(k,77) +rxt(k,78) +rxt(k,79) +rxt(k,80) +rxt(k,81) + &
                 rxt(k,84) +rxt(k,85) +rxt(k,86))*y(k,67) + (rxt(k,71) +rxt(k,72) + &
                 rxt(k,73) +rxt(k,74) +rxt(k,75) +rxt(k,76))*y(k,66) +rxt(k,62) &
                 *y(k,53) +rxt(k,16)*y(k,59)
         loss(k,64) = (rxt(k,198)* y(k,4) +rxt(k,238)* y(k,20) +rxt(k,166)* y(k,26) &
                  + (rxt(k,113) +rxt(k,114) +rxt(k,115))* y(k,38) +rxt(k,144)* y(k,59) &
                  +rxt(k,149)* y(k,60) +rxt(k,139)* y(k,61) +rxt(k,117)* y(k,66) &
                  +rxt(k,118)* y(k,68) +rxt(k,192)* y(k,84) + (rxt(k,159) +rxt(k,160)) &
                 * y(k,85) + 2._r8*rxt(k,127)* y(k,87) +rxt(k,122)* y(k,96) &
                  + rxt(k,264) + het_rates(k,87))* y(k,87)
         prod(k,64) = (rxt(k,221)*y(k,17) +rxt(k,224)*y(k,19) +rxt(k,121)*y(k,41) + &
                 rxt(k,124)*y(k,68) +rxt(k,142)*y(k,61) +rxt(k,173)*y(k,26) + &
                 rxt(k,203)*y(k,4) +rxt(k,242)*y(k,29) +rxt(k,259)*y(k,76) + &
                 .500_r8*rxt(k,263)*y(k,33))*y(k,96) + (rxt(k,155)*y(k,85) + &
                 rxt(k,191)*y(k,84) +rxt(k,235)*y(k,61) +rxt(k,236)*y(k,66))*y(k,16) &
                  + (rxt(k,220)*y(k,17) +rxt(k,223)*y(k,19) +rxt(k,158)*y(k,41)) &
                 *y(k,85) + (rxt(k,162)*y(k,26) +rxt(k,239)*y(k,59))*y(k,20) &
                  + (rxt(k,11) +rxt(k,153))*y(k,50) +rxt(k,244)*y(k,92)*y(k,22) &
                  +rxt(k,116)*y(k,67)*y(k,38) +rxt(k,112)*y(k,66)*y(k,41)
         loss(k,53) = (rxt(k,129)* y(k,66) +rxt(k,130)* y(k,67) +rxt(k,302)* y(k,97) &
                  + het_rates(k,88))* y(k,88)
         prod(k,53) = (.800_r8*rxt(k,288)*y(k,90) +.900_r8*rxt(k,290)*y(k,89))*y(k,86) &
                  +rxt(k,292)*y(k,89)*y(k,66)
         loss(k,42) = ((rxt(k,292) +rxt(k,293))* y(k,66) +rxt(k,291)* y(k,67) &
                  +rxt(k,290)* y(k,86) + het_rates(k,89))* y(k,89)
         prod(k,42) =rxt(k,305)*y(k,98) +rxt(k,310)*y(k,99)
         loss(k,44) = (rxt(k,288)* y(k,86) + het_rates(k,90))* y(k,90)
         prod(k,44) = (rxt(k,298) +rxt(k,297)*y(k,53) +rxt(k,299)*y(k,59))*y(k,95) &
                  +rxt(k,16)*y(k,59) +rxt(k,292)*y(k,89)*y(k,66) +rxt(k,296)*y(k,91) &
                 *y(k,67) +rxt(k,301)*y(k,97)
         loss(k,48) = (rxt(k,294)* y(k,66) + (rxt(k,295) +rxt(k,296))* y(k,67) &
                  + het_rates(k,91))* y(k,91)
         prod(k,48) =rxt(k,62)*y(k,53) +rxt(k,302)*y(k,97)*y(k,88) +rxt(k,311)*y(k,99)
         loss(k,69) = (rxt(k,180)* y(k,7) +rxt(k,181)* y(k,8) +rxt(k,207)* y(k,9) &
                  +rxt(k,182)* y(k,10) +rxt(k,183)* y(k,11) +rxt(k,184)* y(k,12) &
                  +rxt(k,185)* y(k,13) +rxt(k,186)* y(k,14) +rxt(k,230)* y(k,15) &
                  +rxt(k,231)* y(k,17) + (rxt(k,243) +rxt(k,244) +rxt(k,245))* y(k,22) &
                  +rxt(k,208)* y(k,23) +rxt(k,216)* y(k,31) +rxt(k,217)* y(k,32) &
                  +rxt(k,94)* y(k,39) +rxt(k,209)* y(k,40) + (rxt(k,210) +rxt(k,211)) &
                 * y(k,43) +rxt(k,232)* y(k,44) +rxt(k,233)* y(k,45) +rxt(k,234) &
                 * y(k,46) + (rxt(k,187) +rxt(k,188))* y(k,47) + (rxt(k,147) + &
                 rxt(k,148))* y(k,54) + (rxt(k,97) +rxt(k,98))* y(k,67) +rxt(k,99) &
                 * y(k,68) +rxt(k,95)* y(k,100) + rxt(k,96) + het_rates(k,92)) &
                 * y(k,92)
         prod(k,69) = (rxt(k,6) +rxt(k,130)*y(k,88))*y(k,67) +rxt(k,7)*y(k,68) &
                  +.850_r8*rxt(k,289)*y(k,95)*y(k,86) +rxt(k,1)*y(k,100)
         loss(k,24) = (rxt(k,101)* y(k,66) +rxt(k,102)* y(k,67) + rxt(k,92) &
                  + rxt(k,100) + het_rates(k,93))* y(k,93)
         prod(k,24) = (rxt(k,104) +rxt(k,103)*y(k,30) +rxt(k,105)*y(k,66) + &
                 rxt(k,106)*y(k,67) +rxt(k,107)*y(k,68))*y(k,94) +rxt(k,7)*y(k,68)
         loss(k,25) = (rxt(k,103)* y(k,30) +rxt(k,105)* y(k,66) +rxt(k,106)* y(k,67) &
                  +rxt(k,107)* y(k,68) + rxt(k,93) + rxt(k,104) + het_rates(k,94)) &
                 * y(k,94)
         prod(k,25) =rxt(k,97)*y(k,92)*y(k,67)
         loss(k,51) = (rxt(k,297)* y(k,53) +rxt(k,299)* y(k,59) +rxt(k,289)* y(k,86) &
                  + rxt(k,298) + het_rates(k,95))* y(k,95)
         prod(k,51) = (rxt(k,78) +rxt(k,80) +rxt(k,291)*y(k,89) +rxt(k,295)*y(k,91) + &
                 rxt(k,303)*y(k,97) +rxt(k,307)*y(k,98))*y(k,67) +rxt(k,300)*y(k,97) &
                 *y(k,30)
         loss(k,73) = (rxt(k,203)* y(k,4) +rxt(k,219)* y(k,15) +rxt(k,237)* y(k,16) &
                  +rxt(k,221)* y(k,17) +rxt(k,222)* y(k,18) +rxt(k,224)* y(k,19) &
                  +rxt(k,240)* y(k,21) +rxt(k,241)* y(k,22) +rxt(k,226)* y(k,23) &
                  + (rxt(k,173) +rxt(k,174))* y(k,26) +rxt(k,171)* y(k,27) &
                  + (rxt(k,242) +rxt(k,246))* y(k,29) + (rxt(k,262) +rxt(k,263)) &
                 * y(k,33) +rxt(k,120)* y(k,39) +rxt(k,121)* y(k,41) +rxt(k,205) &
                 * y(k,43) +rxt(k,227)* y(k,44) +rxt(k,228)* y(k,45) +rxt(k,229) &
                 * y(k,46) +rxt(k,176)* y(k,47) +rxt(k,152)* y(k,49) +rxt(k,128) &
                 * y(k,50) +rxt(k,179)* y(k,52) +rxt(k,143)* y(k,53) +rxt(k,151) &
                 * y(k,60) +rxt(k,142)* y(k,61) +rxt(k,123)* y(k,66) +rxt(k,124) &
                 * y(k,68) +rxt(k,248)* y(k,70) +rxt(k,253)* y(k,73) +rxt(k,258) &
                 * y(k,75) +rxt(k,259)* y(k,76) +rxt(k,122)* y(k,87) &
                  + 2._r8*(rxt(k,125) +rxt(k,126))* y(k,96) + het_rates(k,96)) &
                 * y(k,96)
         prod(k,73) = (rxt(k,111)*y(k,39) +rxt(k,112)*y(k,41) +rxt(k,117)*y(k,87) + &
                 rxt(k,175)*y(k,47) +rxt(k,178)*y(k,52) +rxt(k,204)*y(k,43) + &
                 rxt(k,206)*y(k,51) +rxt(k,236)*y(k,16))*y(k,66) &
                  + (2.000_r8*rxt(k,114)*y(k,38) +rxt(k,118)*y(k,68) + &
                 rxt(k,139)*y(k,61) +rxt(k,144)*y(k,59) +rxt(k,160)*y(k,85))*y(k,87) &
                  + (rxt(k,243)*y(k,22) +rxt(k,94)*y(k,39) + &
                 2.000_r8*rxt(k,95)*y(k,100) +rxt(k,187)*y(k,47) +rxt(k,210)*y(k,43)) &
                 *y(k,92) + (rxt(k,22) +.300_r8*rxt(k,240)*y(k,96))*y(k,21) &
                  + (rxt(k,3) +rxt(k,214)*y(k,37))*y(k,100) +rxt(k,119)*y(k,68) &
                 *y(k,38) +2.000_r8*rxt(k,4)*y(k,41) +rxt(k,9)*y(k,49) +rxt(k,10) &
                 *y(k,50) +rxt(k,57)*y(k,51) +rxt(k,58)*y(k,52) +.500_r8*rxt(k,266) &
                 *y(k,60)
         loss(k,50) = (rxt(k,300)* y(k,30) +rxt(k,303)* y(k,67) +rxt(k,302)* y(k,88) &
                  + rxt(k,301) + het_rates(k,97))* y(k,97)
         prod(k,50) = (rxt(k,73) +rxt(k,74) +rxt(k,293)*y(k,89) +rxt(k,294)*y(k,91) + &
                 rxt(k,306)*y(k,98) +rxt(k,312)*y(k,99))*y(k,66) + (rxt(k,79) + &
                 rxt(k,81))*y(k,67) + (rxt(k,304)*y(k,98) +rxt(k,309)*y(k,99))*y(k,86) &
                  +rxt(k,286)*y(k,98) +rxt(k,285)*y(k,99)
         loss(k,38) = (rxt(k,306)* y(k,66) +rxt(k,307)* y(k,67) +rxt(k,304)* y(k,86) &
                  + rxt(k,286) + rxt(k,305) + het_rates(k,98))* y(k,98)
         prod(k,38) = (rxt(k,75) +rxt(k,76))*y(k,66) + (rxt(k,85) +rxt(k,86))*y(k,67) &
                  + (rxt(k,287) +rxt(k,308)*y(k,86))*y(k,99)
         loss(k,37) = (rxt(k,312)* y(k,66) + (rxt(k,308) +rxt(k,309))* y(k,86) &
                  + rxt(k,285) + rxt(k,287) + rxt(k,310) + rxt(k,311) &
                  + het_rates(k,99))* y(k,99)
         prod(k,37) = (rxt(k,71) +rxt(k,72))*y(k,66) + (rxt(k,77) +rxt(k,84))*y(k,67)
         loss(k,77) = (rxt(k,214)* y(k,37) +rxt(k,260)* y(k,77) +rxt(k,95)* y(k,92) &
                  + rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,100))* y(k,100)
         prod(k,77) = (rxt(k,219)*y(k,15) +rxt(k,221)*y(k,17) +rxt(k,222)*y(k,18) + &
                 rxt(k,224)*y(k,19) +rxt(k,229)*y(k,46) +rxt(k,241)*y(k,22) + &
                 rxt(k,120)*y(k,39) +rxt(k,121)*y(k,41) +rxt(k,122)*y(k,87) + &
                 rxt(k,125)*y(k,96) +rxt(k,128)*y(k,50) +rxt(k,152)*y(k,49) + &
                 rxt(k,176)*y(k,47) +rxt(k,179)*y(k,52) +rxt(k,205)*y(k,43) + &
                 rxt(k,237)*y(k,16) +rxt(k,240)*y(k,21))*y(k,96) &
                  + (rxt(k,269)*y(k,52) +rxt(k,275)*y(k,52) +rxt(k,276)*y(k,51) + &
                 rxt(k,280)*y(k,52) +rxt(k,281)*y(k,51))*y(k,47) +rxt(k,115)*y(k,87) &
                 *y(k,38) +rxt(k,87)*y(k,42)
      end do
      end subroutine imp_prod_loss
      end module mo_prod_loss
