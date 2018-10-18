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
         loss(k,1) = ( + het_rates(k,4))* y(k,4)
         prod(k,1) = 0._r8
         loss(k,2) = (rxt(k,470)* y(k,95) + rxt(k,31) + het_rates(k,5))* y(k,5)
         prod(k,2) = 0._r8
         loss(k,3) = (rxt(k,471)* y(k,95) + rxt(k,32) + het_rates(k,6))* y(k,6)
         prod(k,3) = 0._r8
         loss(k,4) = (rxt(k,497)* y(k,95) + rxt(k,33) + het_rates(k,7))* y(k,7)
         prod(k,4) = 0._r8
         loss(k,5) = (rxt(k,472)* y(k,95) + rxt(k,34) + het_rates(k,8))* y(k,8)
         prod(k,5) = 0._r8
         loss(k,6) = (rxt(k,473)* y(k,95) + rxt(k,35) + het_rates(k,9))* y(k,9)
         prod(k,6) = 0._r8
         loss(k,7) = (rxt(k,474)* y(k,95) + rxt(k,36) + het_rates(k,10))* y(k,10)
         prod(k,7) = 0._r8
         loss(k,8) = (rxt(k,475)* y(k,95) + rxt(k,37) + het_rates(k,11))* y(k,11)
         prod(k,8) = 0._r8
         loss(k,9) = (rxt(k,476)* y(k,95) + rxt(k,38) + het_rates(k,12))* y(k,12)
         prod(k,9) = 0._r8
         loss(k,10) = (rxt(k,508)* y(k,60) +rxt(k,520)* y(k,95) +rxt(k,509)* y(k,104) &
                  + rxt(k,39) + het_rates(k,13))* y(k,13)
         prod(k,10) = 0._r8
         loss(k,11) = (rxt(k,510)* y(k,60) +rxt(k,521)* y(k,95) +rxt(k,511)* y(k,104) &
                  + rxt(k,40) + het_rates(k,15))* y(k,15)
         prod(k,11) = 0._r8
         loss(k,12) = (rxt(k,512)* y(k,104) + rxt(k,41) + het_rates(k,16))* y(k,16)
         prod(k,12) = 0._r8
         loss(k,13) = (rxt(k,513)* y(k,60) +rxt(k,514)* y(k,104) + rxt(k,42) &
                  + het_rates(k,17))* y(k,17)
         prod(k,13) = 0._r8
         loss(k,14) = (rxt(k,502)* y(k,31) +rxt(k,446)* y(k,60) + (rxt(k,533) + &
                 rxt(k,534) +rxt(k,535))* y(k,95) +rxt(k,531)* y(k,104) + rxt(k,24) &
                  + rxt(k,25) + het_rates(k,20))* y(k,20)
         prod(k,14) = 0._r8
         loss(k,15) = (rxt(k,515)* y(k,60) +rxt(k,498)* y(k,95) +rxt(k,516)* y(k,104) &
                  + rxt(k,43) + het_rates(k,21))* y(k,21)
         prod(k,15) = 0._r8
         loss(k,16) = ( + het_rates(k,26))* y(k,26)
         prod(k,16) = 0._r8
         loss(k,17) = (rxt(k,499)* y(k,95) + rxt(k,51) + het_rates(k,34))* y(k,34)
         prod(k,17) = 0._r8
         loss(k,18) = (rxt(k,522)* y(k,95) +rxt(k,517)* y(k,104) + rxt(k,53) &
                  + het_rates(k,37))* y(k,37)
         prod(k,18) = 0._r8
         loss(k,19) = (rxt(k,523)* y(k,95) +rxt(k,518)* y(k,104) + rxt(k,54) &
                  + het_rates(k,38))* y(k,38)
         prod(k,19) = 0._r8
         loss(k,20) = (rxt(k,524)* y(k,95) +rxt(k,519)* y(k,104) + rxt(k,55) &
                  + het_rates(k,39))* y(k,39)
         prod(k,20) = 0._r8
         loss(k,21) = ((rxt(k,437) +rxt(k,438))* y(k,95) + rxt(k,13) &
                  + het_rates(k,49))* y(k,49)
         prod(k,21) = 0._r8
         loss(k,22) = ( + rxt(k,61) + het_rates(k,58))* y(k,58)
         prod(k,22) = 0._r8
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
         loss(k,6) = ( + rxt(k,27) + het_rates(k,1))* y(k,1)
         prod(k,6) = (rxt(k,549)*y(k,45) +rxt(k,554)*y(k,45))*y(k,40) &
                  +rxt(k,487)*y(k,24)*y(k,2)
         loss(k,52) = (2._r8*rxt(k,484)* y(k,2) + (rxt(k,485) +rxt(k,486) +rxt(k,487)) &
                 * y(k,24) +rxt(k,488)* y(k,43) +rxt(k,489)* y(k,51) +rxt(k,490) &
                 * y(k,52) +rxt(k,492)* y(k,54) +rxt(k,493)* y(k,104) + rxt(k,28) &
                  + het_rates(k,2))* y(k,2)
         prod(k,52) = (rxt(k,29) +rxt(k,491)*y(k,54))*y(k,3) +rxt(k,501)*y(k,95) &
                 *y(k,36) +rxt(k,496)*y(k,54)*y(k,45) +rxt(k,483)*y(k,59)*y(k,56)
         loss(k,18) = (rxt(k,491)* y(k,54) + rxt(k,29) + rxt(k,30) + rxt(k,543) &
                  + rxt(k,546) + rxt(k,551) + het_rates(k,3))* y(k,3)
         prod(k,18) =rxt(k,490)*y(k,52)*y(k,2)
         loss(k,56) = (rxt(k,525)* y(k,53) +rxt(k,526)* y(k,54) +rxt(k,481)* y(k,59) &
                  +rxt(k,445)* y(k,60) +rxt(k,527)* y(k,104) + rxt(k,21) + rxt(k,22) &
                  + het_rates(k,14))* y(k,14)
         prod(k,56) = (rxt(k,452)*y(k,24) +rxt(k,529)*y(k,51))*y(k,18) + (rxt(k,23) + &
                 .300_r8*rxt(k,530)*y(k,104))*y(k,19) + (rxt(k,534)*y(k,95) + &
                 rxt(k,535)*y(k,95))*y(k,20)
         loss(k,39) = (rxt(k,452)* y(k,24) +rxt(k,528)* y(k,43) +rxt(k,529)* y(k,51) &
                  + het_rates(k,18))* y(k,18)
         prod(k,39) = (rxt(k,446)*y(k,60) +rxt(k,502)*y(k,31) +rxt(k,531)*y(k,104) + &
                 rxt(k,533)*y(k,95))*y(k,20) +.700_r8*rxt(k,530)*y(k,104)*y(k,19)
         loss(k,12) = (rxt(k,530)* y(k,104) + rxt(k,23) + het_rates(k,19))* y(k,19)
         prod(k,12) =rxt(k,528)*y(k,43)*y(k,18)
         loss(k,4) = ( + rxt(k,44) + het_rates(k,22))* y(k,22)
         prod(k,4) = (rxt(k,542)*y(k,46) +rxt(k,547)*y(k,25) +rxt(k,548)*y(k,46) + &
                 rxt(k,552)*y(k,25) +rxt(k,553)*y(k,46) +rxt(k,557)*y(k,25))*y(k,40) &
                  +rxt(k,454)*y(k,24)*y(k,24) +rxt(k,458)*y(k,60)*y(k,25)
         loss(k,1) = ( + rxt(k,45) + rxt(k,480) + het_rates(k,23))* y(k,23)
         prod(k,1) =rxt(k,479)*y(k,24)*y(k,24)
         loss(k,81) = ((rxt(k,485) +rxt(k,486) +rxt(k,487))* y(k,2) +rxt(k,452) &
                 * y(k,18) + 2._r8*(rxt(k,453) +rxt(k,454) +rxt(k,455) +rxt(k,479)) &
                 * y(k,24) +rxt(k,456)* y(k,43) +rxt(k,457)* y(k,51) +rxt(k,459) &
                 * y(k,52) +rxt(k,462)* y(k,54) +rxt(k,111)* y(k,65) +rxt(k,123) &
                 * y(k,68) +rxt(k,281)* y(k,81) +rxt(k,310)* y(k,98) + (rxt(k,463) + &
                 rxt(k,464))* y(k,104) +rxt(k,337)* y(k,105) +rxt(k,346)* y(k,106) &
                  + rxt(k,46) + het_rates(k,24))* y(k,24)
         prod(k,81) = (rxt(k,467)*y(k,60) +rxt(k,468)*y(k,54) +rxt(k,469)*y(k,104)) &
                 *y(k,46) + (rxt(k,48) +rxt(k,460)*y(k,54))*y(k,25) &
                  + (rxt(k,450)*y(k,43) +rxt(k,451)*y(k,56))*y(k,60) &
                  +2.000_r8*rxt(k,480)*y(k,23) +rxt(k,478)*y(k,95)*y(k,40) +rxt(k,60) &
                 *y(k,57)
         loss(k,35) = ((rxt(k,547) +rxt(k,552) +rxt(k,557))* y(k,40) +rxt(k,460) &
                 * y(k,54) +rxt(k,458)* y(k,60) +rxt(k,461)* y(k,104) + rxt(k,47) &
                  + rxt(k,48) + rxt(k,545) + rxt(k,550) + rxt(k,556) &
                  + het_rates(k,25))* y(k,25)
         prod(k,35) =rxt(k,459)*y(k,52)*y(k,24)
         loss(k,27) = ((rxt(k,532) +rxt(k,536))* y(k,104) + het_rates(k,27))* y(k,27)
         prod(k,27) = (rxt(k,21) +rxt(k,22) +rxt(k,445)*y(k,60) +rxt(k,481)*y(k,59) + &
                 rxt(k,525)*y(k,53) +rxt(k,526)*y(k,54) +rxt(k,527)*y(k,104))*y(k,14) &
                  + (rxt(k,26) +rxt(k,62) +rxt(k,573)*y(k,107))*y(k,28) &
                  +rxt(k,513)*y(k,60)*y(k,17)
         loss(k,2) = (rxt(k,506)* y(k,95) + rxt(k,49) + het_rates(k,29))* y(k,29)
         prod(k,2) = (rxt(k,471)*y(k,6) +rxt(k,473)*y(k,9) + &
                 2.000_r8*rxt(k,474)*y(k,10) +2.000_r8*rxt(k,475)*y(k,11) + &
                 rxt(k,476)*y(k,12) +rxt(k,497)*y(k,7) +2.000_r8*rxt(k,499)*y(k,34) + &
                 rxt(k,523)*y(k,38) +rxt(k,524)*y(k,39))*y(k,95) &
                  + (rxt(k,518)*y(k,38) +rxt(k,519)*y(k,39))*y(k,104)
         loss(k,5) = (rxt(k,507)* y(k,95) + rxt(k,50) + het_rates(k,30))* y(k,30)
         prod(k,5) = (rxt(k,472)*y(k,8) +rxt(k,473)*y(k,9) +rxt(k,522)*y(k,37)) &
                 *y(k,95) +rxt(k,517)*y(k,104)*y(k,37)
         loss(k,71) = (rxt(k,363)* y(k,88) +rxt(k,307)* y(k,93) +rxt(k,311)* y(k,98) &
                  +rxt(k,325)* y(k,101) +rxt(k,330)* y(k,102) +rxt(k,338)* y(k,105) &
                  +rxt(k,347)* y(k,106) +rxt(k,573)* y(k,107) + rxt(k,26) + rxt(k,62) &
                  + het_rates(k,28))* y(k,28)
         prod(k,71) = (rxt(k,63) +rxt(k,109)*y(k,60) +rxt(k,110)*y(k,60) + &
                 rxt(k,111)*y(k,24) +rxt(k,112)*y(k,32) +rxt(k,119)*y(k,42) + &
                 rxt(k,120)*y(k,54) +rxt(k,121)*y(k,55) +rxt(k,163)*y(k,77) + &
                 rxt(k,165)*y(k,75) +rxt(k,181)*y(k,73) +rxt(k,199)*y(k,92) + &
                 rxt(k,216)*y(k,89) +rxt(k,234)*y(k,88) +rxt(k,251)*y(k,99) + &
                 rxt(k,253)*y(k,75) +rxt(k,260)*y(k,77) +rxt(k,275)*y(k,51) + &
                 rxt(k,276)*y(k,52))*y(k,65) + (rxt(k,115)*y(k,52) + &
                 rxt(k,116)*y(k,52) +rxt(k,117)*y(k,51) +rxt(k,118)*y(k,51) + &
                 rxt(k,150)*y(k,99) +rxt(k,153)*y(k,75) +rxt(k,173)*y(k,77) + &
                 rxt(k,191)*y(k,73) +rxt(k,208)*y(k,92) +rxt(k,226)*y(k,89) + &
                 rxt(k,244)*y(k,88) +rxt(k,255)*y(k,75) +rxt(k,256)*y(k,77))*y(k,67) &
                  + (rxt(k,65) +rxt(k,122)*y(k,60) +rxt(k,123)*y(k,24) + &
                 rxt(k,125)*y(k,40) +rxt(k,127)*y(k,56) +rxt(k,146)*y(k,99) + &
                 rxt(k,169)*y(k,77) +rxt(k,186)*y(k,73) +rxt(k,204)*y(k,92) + &
                 rxt(k,220)*y(k,75) +rxt(k,222)*y(k,89) +rxt(k,239)*y(k,88))*y(k,68) &
                  + (rxt(k,148)*y(k,99) +rxt(k,171)*y(k,77) +rxt(k,189)*y(k,73) + &
                 rxt(k,206)*y(k,92) +rxt(k,224)*y(k,89) +rxt(k,241)*y(k,88) + &
                 rxt(k,242)*y(k,75) +rxt(k,254)*y(k,77) +rxt(k,266)*y(k,75))*y(k,66) &
                  + (rxt(k,144)*y(k,99) +rxt(k,167)*y(k,77) +rxt(k,184)*y(k,73) + &
                 rxt(k,198)*y(k,75) +rxt(k,202)*y(k,92) +rxt(k,219)*y(k,89) + &
                 rxt(k,237)*y(k,88))*y(k,71) + (rxt(k,364) +rxt(k,301)*y(k,69) + &
                 rxt(k,302)*y(k,110))*y(k,91) + (rxt(k,532)*y(k,104) + &
                 rxt(k,536)*y(k,104))*y(k,27)
         loss(k,31) = (rxt(k,502)* y(k,20) +rxt(k,503)* y(k,33) +rxt(k,505)* y(k,42) &
                  +rxt(k,504)* y(k,110) + het_rates(k,31))* y(k,31)
         prod(k,31) = (rxt(k,475)*y(k,11) +rxt(k,497)*y(k,7) + &
                 2.000_r8*rxt(k,506)*y(k,29) +rxt(k,507)*y(k,30))*y(k,95) &
                  +2.000_r8*rxt(k,49)*y(k,29) +rxt(k,50)*y(k,30) +rxt(k,57)*y(k,41)
         loss(k,86) = ((rxt(k,403) +rxt(k,404) +rxt(k,405))* y(k,43) +rxt(k,406) &
                 * y(k,55) +rxt(k,409)* y(k,56) +rxt(k,100)* y(k,61) +rxt(k,112) &
                 * y(k,65) +rxt(k,124)* y(k,68) +rxt(k,282)* y(k,81) +rxt(k,304) &
                 * y(k,92) +rxt(k,312)* y(k,98) +rxt(k,326)* y(k,101) +rxt(k,339) &
                 * y(k,105) + het_rates(k,32))* y(k,32)
         prod(k,86) = (rxt(k,136)*y(k,69) +rxt(k,142)*y(k,61) +rxt(k,153)*y(k,67) + &
                 rxt(k,157)*y(k,82) +rxt(k,158)*y(k,86) +rxt(k,159)*y(k,62) + &
                 rxt(k,160)*y(k,84) +rxt(k,161)*y(k,81) +rxt(k,165)*y(k,65) + &
                 rxt(k,176)*y(k,63) +rxt(k,198)*y(k,71) +rxt(k,209)*y(k,98) + &
                 rxt(k,220)*y(k,68) +rxt(k,231)*y(k,85) +rxt(k,242)*y(k,66) + &
                 rxt(k,253)*y(k,65) +rxt(k,255)*y(k,67) +rxt(k,257)*y(k,85) + &
                 rxt(k,266)*y(k,66))*y(k,75) + (rxt(k,139)*y(k,69) + &
                 rxt(k,163)*y(k,65) +rxt(k,164)*y(k,63) +rxt(k,167)*y(k,71) + &
                 rxt(k,168)*y(k,98) +rxt(k,169)*y(k,68) +rxt(k,170)*y(k,85) + &
                 rxt(k,171)*y(k,66) +rxt(k,172)*y(k,61) +rxt(k,173)*y(k,67) + &
                 rxt(k,174)*y(k,82) +rxt(k,175)*y(k,86) +rxt(k,177)*y(k,62) + &
                 rxt(k,178)*y(k,84) +rxt(k,179)*y(k,81) +rxt(k,254)*y(k,66) + &
                 rxt(k,256)*y(k,67) +rxt(k,258)*y(k,85) +rxt(k,260)*y(k,65))*y(k,77) &
                  + (rxt(k,181)*y(k,65) +rxt(k,182)*y(k,63) +rxt(k,184)*y(k,71) + &
                 rxt(k,185)*y(k,98) +rxt(k,186)*y(k,68) +rxt(k,188)*y(k,85) + &
                 rxt(k,189)*y(k,66) +rxt(k,190)*y(k,61) +rxt(k,191)*y(k,67) + &
                 rxt(k,192)*y(k,82) +rxt(k,193)*y(k,86) +rxt(k,194)*y(k,62) + &
                 rxt(k,195)*y(k,84) +rxt(k,196)*y(k,81) +rxt(k,378)*y(k,69))*y(k,73) &
                  + (rxt(k,349)*y(k,106) +rxt(k,384)*y(k,95) +rxt(k,401)*y(k,54) + &
                 rxt(k,410)*y(k,104) +rxt(k,447)*y(k,60) +rxt(k,503)*y(k,31))*y(k,33) &
                  + (rxt(k,413)*y(k,54) +rxt(k,433)*y(k,48) +rxt(k,527)*y(k,14) + &
                 rxt(k,536)*y(k,27))*y(k,104) + (rxt(k,133)*y(k,70) + &
                 rxt(k,376)*y(k,78) +rxt(k,377)*y(k,72))*y(k,69) &
                  + (rxt(k,534)*y(k,20) +rxt(k,478)*y(k,40) +rxt(k,501)*y(k,36)) &
                 *y(k,95) + (2.000_r8*rxt(k,2) +rxt(k,3))*y(k,110) +2.000_r8*rxt(k,21) &
                 *y(k,14) +rxt(k,23)*y(k,19) +rxt(k,52)*y(k,36) +rxt(k,56)*y(k,40) &
                  +rxt(k,57)*y(k,41)
         loss(k,50) = (rxt(k,503)* y(k,31) +rxt(k,401)* y(k,54) +rxt(k,447)* y(k,60) &
                  +rxt(k,384)* y(k,95) +rxt(k,410)* y(k,104) + (rxt(k,348) + &
                 rxt(k,349))* y(k,106) + het_rates(k,33))* y(k,33)
         prod(k,50) =rxt(k,22)*y(k,14) +rxt(k,535)*y(k,95)*y(k,20) +rxt(k,403)*y(k,43) &
                 *y(k,32) +rxt(k,1)*y(k,110)
         loss(k,19) = (rxt(k,402)* y(k,54) +rxt(k,448)* y(k,60) +rxt(k,411)* y(k,104) &
                  + rxt(k,4) + het_rates(k,35))* y(k,35)
         prod(k,19) = (.500_r8*rxt(k,537) +rxt(k,417)*y(k,43))*y(k,43) &
                  +rxt(k,416)*y(k,104)*y(k,104)
         loss(k,32) = (rxt(k,494)* y(k,54) + (rxt(k,500) +rxt(k,501))* y(k,95) &
                  +rxt(k,495)* y(k,104) + rxt(k,52) + het_rates(k,36))* y(k,36)
         prod(k,32) = (rxt(k,481)*y(k,14) +rxt(k,482)*y(k,43))*y(k,59)
         loss(k,72) = ((rxt(k,547) +rxt(k,552) +rxt(k,557))* y(k,25) + (rxt(k,549) + &
                 rxt(k,554))* y(k,45) + (rxt(k,542) +rxt(k,548) +rxt(k,553))* y(k,46) &
                  +rxt(k,465)* y(k,54) +rxt(k,103)* y(k,61) +rxt(k,101)* y(k,62) &
                  +rxt(k,125)* y(k,68) +rxt(k,284)* y(k,81) + (rxt(k,271) +rxt(k,293)) &
                 * y(k,83) + (rxt(k,477) +rxt(k,478))* y(k,95) +rxt(k,313)* y(k,98) &
                  +rxt(k,466)* y(k,104) +rxt(k,340)* y(k,105) +rxt(k,351)* y(k,106) &
                  + rxt(k,56) + het_rates(k,40))* y(k,40)
         prod(k,72) = (rxt(k,446)*y(k,20) +rxt(k,508)*y(k,13) +rxt(k,510)*y(k,15) + &
                 2.000_r8*rxt(k,513)*y(k,17) +rxt(k,515)*y(k,21) +rxt(k,445)*y(k,14) + &
                 rxt(k,447)*y(k,33) +rxt(k,448)*y(k,35) +rxt(k,449)*y(k,43) + &
                 rxt(k,467)*y(k,46))*y(k,60) + (rxt(k,381) +rxt(k,164)*y(k,77) + &
                 rxt(k,176)*y(k,75) +rxt(k,182)*y(k,73) +rxt(k,200)*y(k,92) + &
                 rxt(k,217)*y(k,89) +rxt(k,235)*y(k,88) +rxt(k,252)*y(k,99) + &
                 2.000_r8*rxt(k,262)*y(k,75) +2.000_r8*rxt(k,263)*y(k,77))*y(k,63) &
                  + (rxt(k,152)*y(k,99) +rxt(k,158)*y(k,75) +rxt(k,175)*y(k,77) + &
                 rxt(k,193)*y(k,73) +rxt(k,211)*y(k,92) +rxt(k,228)*y(k,89) + &
                 rxt(k,246)*y(k,88) +rxt(k,294)*y(k,42))*y(k,86) &
                  + (rxt(k,100)*y(k,32) +rxt(k,104)*y(k,42))*y(k,61) &
                  +rxt(k,464)*y(k,104)*y(k,24)
         loss(k,7) = ( + rxt(k,57) + het_rates(k,41))* y(k,41)
         prod(k,7) = (rxt(k,502)*y(k,20) +rxt(k,503)*y(k,33) +rxt(k,504)*y(k,110) + &
                 rxt(k,505)*y(k,42))*y(k,31)
         loss(k,73) = (rxt(k,505)* y(k,31) +rxt(k,104)* y(k,61) +rxt(k,119)* y(k,65) &
                  +rxt(k,285)* y(k,81) +rxt(k,295)* y(k,83) +rxt(k,290)* y(k,85) &
                  +rxt(k,294)* y(k,86) +rxt(k,314)* y(k,98) +rxt(k,442)* y(k,104) &
                  +rxt(k,352)* y(k,106) + rxt(k,9) + het_rates(k,42))* y(k,42)
         prod(k,73) = (rxt(k,270) +2.000_r8*rxt(k,141)*y(k,75) + &
                 2.000_r8*rxt(k,162)*y(k,77) +2.000_r8*rxt(k,180)*y(k,73) + &
                 rxt(k,197)*y(k,92) +rxt(k,215)*y(k,89) +rxt(k,233)*y(k,88) + &
                 rxt(k,250)*y(k,99) +2.000_r8*rxt(k,264)*y(k,75) + &
                 2.000_r8*rxt(k,265)*y(k,77))*y(k,87) + (2.000_r8*rxt(k,538) + &
                 2.000_r8*rxt(k,541) +2.000_r8*rxt(k,544) +2.000_r8*rxt(k,555) + &
                 rxt(k,137)*y(k,75) +rxt(k,140)*y(k,77) +rxt(k,288)*y(k,84) + &
                 rxt(k,292)*y(k,85))*y(k,50) + (rxt(k,545) +rxt(k,550) +rxt(k,556) + &
                 rxt(k,547)*y(k,40) +rxt(k,552)*y(k,40) +rxt(k,557)*y(k,40))*y(k,25) &
                  + (rxt(k,166)*y(k,77) +rxt(k,183)*y(k,73) +rxt(k,187)*y(k,75) + &
                 rxt(k,259)*y(k,75) +rxt(k,261)*y(k,77) +rxt(k,293)*y(k,40))*y(k,83) &
                  + (rxt(k,543) +rxt(k,546) +rxt(k,551))*y(k,3) &
                  + (.500_r8*rxt(k,539) +rxt(k,441)*y(k,104))*y(k,52) + (rxt(k,540) + &
                 rxt(k,525)*y(k,14))*y(k,53) + (rxt(k,135)*y(k,74) + &
                 rxt(k,138)*y(k,76))*y(k,110)
         loss(k,61) = (rxt(k,488)* y(k,2) +rxt(k,528)* y(k,18) +rxt(k,456)* y(k,24) &
                  + (rxt(k,403) +rxt(k,404) +rxt(k,405))* y(k,32) + 2._r8*rxt(k,417) &
                 * y(k,43) +rxt(k,434)* y(k,51) +rxt(k,439)* y(k,52) +rxt(k,429) &
                 * y(k,53) +rxt(k,407)* y(k,54) +rxt(k,408)* y(k,56) +rxt(k,482) &
                 * y(k,59) + (rxt(k,449) +rxt(k,450))* y(k,60) +rxt(k,305)* y(k,92) &
                  +rxt(k,412)* y(k,104) + rxt(k,537) + het_rates(k,43))* y(k,43)
         prod(k,61) = (rxt(k,511)*y(k,15) +rxt(k,514)*y(k,17) +rxt(k,411)*y(k,35) + &
                 rxt(k,414)*y(k,56) +rxt(k,432)*y(k,53) +rxt(k,463)*y(k,24) + &
                 rxt(k,493)*y(k,2) +rxt(k,532)*y(k,27))*y(k,104) &
                  + (rxt(k,445)*y(k,60) +rxt(k,481)*y(k,59) +rxt(k,525)*y(k,53) + &
                 rxt(k,526)*y(k,54))*y(k,14) + (rxt(k,510)*y(k,15) + &
                 rxt(k,513)*y(k,17) +rxt(k,448)*y(k,35))*y(k,60) &
                  + (rxt(k,312)*y(k,32) +rxt(k,313)*y(k,40) +rxt(k,314)*y(k,42)) &
                 *y(k,98) + (rxt(k,452)*y(k,24) +rxt(k,529)*y(k,51))*y(k,18) &
                  + (rxt(k,11) +rxt(k,443))*y(k,44) + (rxt(k,342)*y(k,105) + &
                 rxt(k,402)*y(k,35))*y(k,54) +rxt(k,534)*y(k,95)*y(k,20) &
                  +rxt(k,406)*y(k,55)*y(k,32) +rxt(k,125)*y(k,68)*y(k,40)
         loss(k,14) = (rxt(k,418)* y(k,104) + rxt(k,10) + rxt(k,11) + rxt(k,443) &
                  + het_rates(k,44))* y(k,44)
         prod(k,14) =rxt(k,439)*y(k,52)*y(k,43)
         loss(k,30) = ((rxt(k,549) +rxt(k,554))* y(k,40) +rxt(k,496)* y(k,54) &
                  + rxt(k,58) + het_rates(k,45))* y(k,45)
         prod(k,30) = (rxt(k,543) +rxt(k,546) +rxt(k,551))*y(k,3) +rxt(k,488)*y(k,43) &
                 *y(k,2)
         loss(k,33) = ((rxt(k,542) +rxt(k,548) +rxt(k,553))* y(k,40) +rxt(k,468) &
                 * y(k,54) +rxt(k,467)* y(k,60) +rxt(k,469)* y(k,104) + rxt(k,59) &
                  + het_rates(k,46))* y(k,46)
         prod(k,33) = (rxt(k,545) +rxt(k,550) +rxt(k,556) +rxt(k,461)*y(k,104)) &
                 *y(k,25) +rxt(k,456)*y(k,43)*y(k,24)
         loss(k,26) = (rxt(k,335)* y(k,104) + rxt(k,12) + het_rates(k,47))* y(k,47)
         prod(k,26) = (rxt(k,284)*y(k,40) +rxt(k,285)*y(k,42))*y(k,81) &
                  +rxt(k,344)*y(k,104)*y(k,51) +rxt(k,300)*y(k,110)*y(k,90)
         loss(k,40) = (rxt(k,421)* y(k,51) + (rxt(k,422) +rxt(k,423) +rxt(k,424)) &
                 * y(k,52) +rxt(k,425)* y(k,55) +rxt(k,570)* y(k,99) +rxt(k,433) &
                 * y(k,104) + rxt(k,66) + het_rates(k,48))* y(k,48)
         prod(k,40) = (rxt(k,419)*y(k,79) +rxt(k,567)*y(k,94))*y(k,54) &
                  + (.200_r8*rxt(k,561)*y(k,88) +1.100_r8*rxt(k,563)*y(k,80))*y(k,69) &
                  +rxt(k,17)*y(k,51) +rxt(k,568)*y(k,94)*y(k,55) +rxt(k,574)*y(k,107)
         loss(k,38) = (rxt(k,137)* y(k,75) +rxt(k,140)* y(k,77) +rxt(k,288)* y(k,84) &
                  +rxt(k,292)* y(k,85) + rxt(k,14) + rxt(k,15) + rxt(k,444) &
                  + rxt(k,538) + rxt(k,541) + rxt(k,544) + rxt(k,555) &
                  + het_rates(k,50))* y(k,50)
         prod(k,38) =rxt(k,440)*y(k,53)*y(k,52)
         loss(k,69) = (rxt(k,489)* y(k,2) +rxt(k,529)* y(k,18) +rxt(k,457)* y(k,24) &
                  +rxt(k,434)* y(k,43) +rxt(k,421)* y(k,48) +rxt(k,430)* y(k,53) &
                  +rxt(k,436)* y(k,54) +rxt(k,435)* y(k,56) + (rxt(k,106) +rxt(k,107)) &
                 * y(k,64) +rxt(k,275)* y(k,65) + (rxt(k,117) +rxt(k,118))* y(k,67) &
                  +rxt(k,572)* y(k,99) + (rxt(k,267) +rxt(k,274))* y(k,101) &
                  +rxt(k,344)* y(k,104) +rxt(k,131)* y(k,106) + rxt(k,16) + rxt(k,17) &
                  + het_rates(k,51))* y(k,51)
         prod(k,69) = (rxt(k,197)*y(k,87) +rxt(k,199)*y(k,65) +rxt(k,200)*y(k,63) + &
                 rxt(k,201)*y(k,83) +rxt(k,202)*y(k,71) +rxt(k,203)*y(k,98) + &
                 rxt(k,204)*y(k,68) +rxt(k,205)*y(k,85) +rxt(k,206)*y(k,66) + &
                 rxt(k,207)*y(k,61) +rxt(k,208)*y(k,67) +rxt(k,210)*y(k,82) + &
                 rxt(k,211)*y(k,86) +rxt(k,212)*y(k,62) +rxt(k,213)*y(k,84) + &
                 rxt(k,214)*y(k,81) +rxt(k,303)*y(k,69) +rxt(k,304)*y(k,32))*y(k,92) &
                  + (rxt(k,215)*y(k,87) +rxt(k,216)*y(k,65) +rxt(k,217)*y(k,63) + &
                 rxt(k,218)*y(k,83) +rxt(k,219)*y(k,71) +rxt(k,221)*y(k,98) + &
                 rxt(k,222)*y(k,68) +rxt(k,223)*y(k,85) +rxt(k,224)*y(k,66) + &
                 rxt(k,225)*y(k,61) +rxt(k,226)*y(k,67) +rxt(k,227)*y(k,82) + &
                 rxt(k,228)*y(k,86) +rxt(k,229)*y(k,62) +rxt(k,230)*y(k,84) + &
                 rxt(k,232)*y(k,81) +rxt(k,298)*y(k,69))*y(k,89) &
                  + (rxt(k,233)*y(k,87) +rxt(k,234)*y(k,65) +rxt(k,235)*y(k,63) + &
                 rxt(k,236)*y(k,83) +rxt(k,237)*y(k,71) +rxt(k,238)*y(k,98) + &
                 rxt(k,239)*y(k,68) +rxt(k,240)*y(k,85) +rxt(k,241)*y(k,66) + &
                 rxt(k,243)*y(k,61) +rxt(k,244)*y(k,67) +rxt(k,245)*y(k,82) + &
                 rxt(k,246)*y(k,86) +rxt(k,247)*y(k,62) +rxt(k,248)*y(k,84) + &
                 rxt(k,249)*y(k,81))*y(k,88) + (rxt(k,18) +.500_r8*rxt(k,539) + &
                 rxt(k,286)*y(k,81) +2.000_r8*rxt(k,423)*y(k,48) +rxt(k,426)*y(k,54)) &
                 *y(k,52) + (rxt(k,299)*y(k,90) +rxt(k,301)*y(k,91) + &
                 rxt(k,379)*y(k,93))*y(k,69) + (rxt(k,425)*y(k,55) + &
                 rxt(k,433)*y(k,104))*y(k,48) +rxt(k,282)*y(k,81)*y(k,32) +rxt(k,12) &
                 *y(k,47) +2.000_r8*rxt(k,437)*y(k,95)*y(k,49) +rxt(k,15)*y(k,50) &
                  +rxt(k,20)*y(k,53) +rxt(k,420)*y(k,79)*y(k,55) +rxt(k,571)*y(k,99) &
                  +rxt(k,583)*y(k,109)
         loss(k,87) = (rxt(k,490)* y(k,2) +rxt(k,459)* y(k,24) +rxt(k,439)* y(k,43) &
                  + (rxt(k,422) +rxt(k,423) +rxt(k,424))* y(k,48) +rxt(k,440)* y(k,53) &
                  + (rxt(k,426) +rxt(k,428))* y(k,54) +rxt(k,427)* y(k,56) +rxt(k,105) &
                 * y(k,61) +rxt(k,276)* y(k,65) + (rxt(k,115) +rxt(k,116))* y(k,67) &
                  +rxt(k,286)* y(k,81) +rxt(k,315)* y(k,98) + (rxt(k,272) +rxt(k,273)) &
                 * y(k,101) +rxt(k,441)* y(k,104) +rxt(k,341)* y(k,105) +rxt(k,354) &
                 * y(k,106) + rxt(k,18) + rxt(k,539) + het_rates(k,52))* y(k,52)
         prod(k,87) = (rxt(k,107)*y(k,64) +rxt(k,131)*y(k,106) + &
                 2.000_r8*rxt(k,430)*y(k,53) +rxt(k,434)*y(k,43) +rxt(k,435)*y(k,56) + &
                 rxt(k,436)*y(k,54) +rxt(k,457)*y(k,24) +rxt(k,489)*y(k,2) + &
                 rxt(k,529)*y(k,18))*y(k,51) + (rxt(k,75) +rxt(k,156)*y(k,99) + &
                 rxt(k,161)*y(k,75) +rxt(k,179)*y(k,77) +rxt(k,196)*y(k,73) + &
                 rxt(k,214)*y(k,92) +rxt(k,232)*y(k,89) +rxt(k,249)*y(k,88) + &
                 rxt(k,280)*y(k,60))*y(k,81) + (rxt(k,151)*y(k,99) + &
                 rxt(k,157)*y(k,75) +rxt(k,174)*y(k,77) +rxt(k,192)*y(k,73) + &
                 rxt(k,210)*y(k,92) +rxt(k,227)*y(k,89) +rxt(k,245)*y(k,88))*y(k,82) &
                  + (rxt(k,19) +rxt(k,429)*y(k,43) +rxt(k,431)*y(k,54) + &
                 rxt(k,432)*y(k,104))*y(k,53) + (rxt(k,11) +rxt(k,443) + &
                 rxt(k,418)*y(k,104))*y(k,44) + (rxt(k,14) +rxt(k,444))*y(k,50) &
                  + (rxt(k,306)*y(k,92) +rxt(k,335)*y(k,47))*y(k,104) +rxt(k,29) &
                 *y(k,3) +rxt(k,48)*y(k,25) +rxt(k,9)*y(k,42)
         loss(k,79) = (rxt(k,525)* y(k,14) +rxt(k,429)* y(k,43) +rxt(k,430)* y(k,51) &
                  +rxt(k,440)* y(k,52) +rxt(k,431)* y(k,54) +rxt(k,432)* y(k,104) &
                  + rxt(k,19) + rxt(k,20) + rxt(k,540) + het_rates(k,53))* y(k,53)
         prod(k,79) = (rxt(k,147)*y(k,99) +rxt(k,170)*y(k,77) +rxt(k,188)*y(k,73) + &
                 rxt(k,205)*y(k,92) +rxt(k,223)*y(k,89) +rxt(k,231)*y(k,75) + &
                 rxt(k,240)*y(k,88) +rxt(k,257)*y(k,75) +rxt(k,258)*y(k,77))*y(k,85) &
                  + (rxt(k,155)*y(k,99) +rxt(k,160)*y(k,75) +rxt(k,178)*y(k,77) + &
                 rxt(k,195)*y(k,73) +rxt(k,213)*y(k,92) +rxt(k,230)*y(k,89) + &
                 rxt(k,248)*y(k,88))*y(k,84) + (rxt(k,152)*y(k,99) + &
                 rxt(k,158)*y(k,75) +rxt(k,175)*y(k,77) +rxt(k,193)*y(k,73) + &
                 rxt(k,211)*y(k,92) +rxt(k,228)*y(k,89) +rxt(k,246)*y(k,88))*y(k,86) &
                  + (rxt(k,76) +rxt(k,143)*y(k,99) +rxt(k,201)*y(k,92) + &
                 rxt(k,218)*y(k,89) +rxt(k,236)*y(k,88))*y(k,83) + (rxt(k,47) + &
                 rxt(k,458)*y(k,60) +rxt(k,460)*y(k,54) +rxt(k,461)*y(k,104))*y(k,25) &
                  + (rxt(k,197)*y(k,92) +rxt(k,215)*y(k,89) +rxt(k,233)*y(k,88) + &
                 rxt(k,250)*y(k,99))*y(k,87) + (rxt(k,14) +rxt(k,15) +rxt(k,444)) &
                 *y(k,50) + (rxt(k,30) +rxt(k,491)*y(k,54))*y(k,3) &
                  + (rxt(k,442)*y(k,104) +rxt(k,505)*y(k,31))*y(k,42) &
                  + (rxt(k,427)*y(k,56) +rxt(k,428)*y(k,54))*y(k,52) &
                  +rxt(k,281)*y(k,81)*y(k,24) +rxt(k,305)*y(k,92)*y(k,43) +rxt(k,10) &
                 *y(k,44)
         loss(k,74) = (rxt(k,492)* y(k,2) +rxt(k,491)* y(k,3) +rxt(k,526)* y(k,14) &
                  +rxt(k,462)* y(k,24) +rxt(k,460)* y(k,25) +rxt(k,401)* y(k,33) &
                  +rxt(k,402)* y(k,35) +rxt(k,494)* y(k,36) +rxt(k,465)* y(k,40) &
                  +rxt(k,407)* y(k,43) +rxt(k,496)* y(k,45) +rxt(k,468)* y(k,46) &
                  +rxt(k,436)* y(k,51) + (rxt(k,426) +rxt(k,428))* y(k,52) +rxt(k,431) &
                 * y(k,53) + 2._r8*rxt(k,399)* y(k,54) +rxt(k,400)* y(k,55) &
                  +rxt(k,398)* y(k,56) +rxt(k,108)* y(k,64) +rxt(k,120)* y(k,65) &
                  +rxt(k,126)* y(k,68) + (rxt(k,565) +rxt(k,566))* y(k,80) +rxt(k,296) &
                 * y(k,83) +rxt(k,567)* y(k,94) + (rxt(k,319) +rxt(k,320))* y(k,98) &
                  + (rxt(k,328) +rxt(k,329))* y(k,101) +rxt(k,331)* y(k,102) &
                  +rxt(k,333)* y(k,103) +rxt(k,413)* y(k,104) +rxt(k,342)* y(k,105) &
                  +rxt(k,355)* y(k,106) + rxt(k,77) + rxt(k,78) + rxt(k,79) &
                  + rxt(k,80) + rxt(k,81) + rxt(k,82) + het_rates(k,54))* y(k,54)
         prod(k,74) = (2.000_r8*rxt(k,5) +rxt(k,6) +rxt(k,83) +2.000_r8*rxt(k,84) + &
                 rxt(k,85) +rxt(k,87) +2.000_r8*rxt(k,89) +rxt(k,90) +rxt(k,91) + &
                 rxt(k,92) +rxt(k,387)*y(k,95) +rxt(k,388)*y(k,95) + &
                 rxt(k,425)*y(k,48) +rxt(k,569)*y(k,94) +rxt(k,575)*y(k,107) + &
                 rxt(k,579)*y(k,108))*y(k,55) + (rxt(k,109)*y(k,60) + &
                 rxt(k,163)*y(k,77) +rxt(k,165)*y(k,75) +rxt(k,181)*y(k,73) + &
                 rxt(k,199)*y(k,92) +rxt(k,216)*y(k,89) +rxt(k,234)*y(k,88) + &
                 rxt(k,251)*y(k,99) +rxt(k,253)*y(k,75) +rxt(k,260)*y(k,77))*y(k,65) &
                  + (rxt(k,148)*y(k,99) +rxt(k,171)*y(k,77) +rxt(k,189)*y(k,73) + &
                 rxt(k,206)*y(k,92) +rxt(k,224)*y(k,89) +rxt(k,241)*y(k,88) + &
                 rxt(k,242)*y(k,75) +rxt(k,254)*y(k,77) +rxt(k,266)*y(k,75))*y(k,66) &
                  + (rxt(k,150)*y(k,99) +rxt(k,153)*y(k,75) +rxt(k,173)*y(k,77) + &
                 rxt(k,191)*y(k,73) +rxt(k,208)*y(k,92) +rxt(k,226)*y(k,89) + &
                 rxt(k,244)*y(k,88) +rxt(k,255)*y(k,75) +rxt(k,256)*y(k,77))*y(k,67) &
                  + (rxt(k,99) +rxt(k,353) +rxt(k,345)*y(k,60) +rxt(k,354)*y(k,52) + &
                 rxt(k,358)*y(k,56))*y(k,106) + (rxt(k,421)*y(k,51) + &
                 rxt(k,422)*y(k,52) +rxt(k,570)*y(k,99))*y(k,48) + (rxt(k,26) + &
                 rxt(k,62))*y(k,28) + (rxt(k,17) +rxt(k,267)*y(k,101))*y(k,51) &
                  + (rxt(k,561)*y(k,88) +1.150_r8*rxt(k,562)*y(k,99))*y(k,69) &
                  +rxt(k,28)*y(k,2) +rxt(k,46)*y(k,24) +rxt(k,405)*y(k,43)*y(k,32) &
                  +rxt(k,15)*y(k,50) +rxt(k,18)*y(k,52) +rxt(k,19)*y(k,53) +rxt(k,8) &
                 *y(k,56) +rxt(k,60)*y(k,57) +rxt(k,386)*y(k,95) +rxt(k,415)*y(k,104) &
                 *y(k,104) +rxt(k,577)*y(k,108) +rxt(k,582)*y(k,109) +rxt(k,2) &
                 *y(k,110)
         loss(k,85) = (rxt(k,406)* y(k,32) +rxt(k,425)* y(k,48) +rxt(k,400)* y(k,54) &
                  +rxt(k,121)* y(k,65) + (rxt(k,128) +rxt(k,130))* y(k,69) +rxt(k,420) &
                 * y(k,79) +rxt(k,564)* y(k,80) + (rxt(k,568) +rxt(k,569))* y(k,94) &
                  +rxt(k,387)* y(k,95) +rxt(k,392)* y(k,96) +rxt(k,317)* y(k,98) &
                  +rxt(k,359)* y(k,99) +rxt(k,357)* y(k,106) +rxt(k,575)* y(k,107) &
                  +rxt(k,579)* y(k,108) + rxt(k,5) + rxt(k,6) + rxt(k,83) + rxt(k,84) &
                  + rxt(k,85) + rxt(k,86) + rxt(k,87) + rxt(k,88) + rxt(k,89) &
                  + rxt(k,90) + rxt(k,91) + rxt(k,92) + het_rates(k,55))* y(k,55)
         prod(k,85) = (rxt(k,108)*y(k,64) +rxt(k,126)*y(k,68) +rxt(k,296)*y(k,83) + &
                 rxt(k,320)*y(k,98) +2.000_r8*rxt(k,328)*y(k,101) + &
                 rxt(k,329)*y(k,101) +rxt(k,331)*y(k,102) +rxt(k,355)*y(k,106) + &
                 rxt(k,391)*y(k,96) +2.000_r8*rxt(k,398)*y(k,56) +rxt(k,399)*y(k,54) + &
                 rxt(k,407)*y(k,43) +rxt(k,413)*y(k,104) +rxt(k,426)*y(k,52) + &
                 rxt(k,431)*y(k,53) +rxt(k,462)*y(k,24) +rxt(k,492)*y(k,2))*y(k,54) &
                  + (rxt(k,143)*y(k,83) +rxt(k,144)*y(k,71) + &
                 2.000_r8*rxt(k,146)*y(k,68) +rxt(k,147)*y(k,85) +rxt(k,148)*y(k,66) + &
                 rxt(k,149)*y(k,61) +rxt(k,150)*y(k,67) +rxt(k,151)*y(k,82) + &
                 rxt(k,152)*y(k,86) +rxt(k,154)*y(k,62) +rxt(k,155)*y(k,84) + &
                 rxt(k,156)*y(k,81) +rxt(k,250)*y(k,87) +rxt(k,251)*y(k,65) + &
                 rxt(k,252)*y(k,63) +rxt(k,572)*y(k,51))*y(k,99) + (rxt(k,8) + &
                 rxt(k,127)*y(k,68) +rxt(k,129)*y(k,69) +rxt(k,287)*y(k,81) + &
                 2.000_r8*rxt(k,297)*y(k,83) +rxt(k,318)*y(k,98) + &
                 3.000_r8*rxt(k,327)*y(k,101) +2.000_r8*rxt(k,389)*y(k,95) + &
                 2.000_r8*rxt(k,408)*y(k,43) +rxt(k,409)*y(k,32) + &
                 rxt(k,414)*y(k,104) +rxt(k,427)*y(k,52) +rxt(k,435)*y(k,51) + &
                 rxt(k,451)*y(k,60) +rxt(k,483)*y(k,59))*y(k,56) + (rxt(k,93) + &
                 rxt(k,132) +rxt(k,168)*y(k,77) +rxt(k,185)*y(k,73) + &
                 rxt(k,203)*y(k,92) +rxt(k,209)*y(k,75) +rxt(k,221)*y(k,89) + &
                 rxt(k,238)*y(k,88) +rxt(k,309)*y(k,60) +rxt(k,310)*y(k,24) + &
                 rxt(k,315)*y(k,52) +2.000_r8*rxt(k,316)*y(k,96))*y(k,98) &
                  + (rxt(k,111)*y(k,65) +rxt(k,123)*y(k,68) +rxt(k,346)*y(k,106) + &
                 rxt(k,453)*y(k,24) +rxt(k,454)*y(k,24) +rxt(k,456)*y(k,43) + &
                 rxt(k,464)*y(k,104) +rxt(k,486)*y(k,2) +rxt(k,487)*y(k,2))*y(k,24) &
                  + (rxt(k,403)*y(k,32) +rxt(k,412)*y(k,104) +rxt(k,417)*y(k,43) + &
                 rxt(k,429)*y(k,53) +rxt(k,449)*y(k,60) +rxt(k,482)*y(k,59) + &
                 rxt(k,488)*y(k,2) +rxt(k,528)*y(k,18))*y(k,43) &
                  + (rxt(k,122)*y(k,60) +rxt(k,169)*y(k,77) +rxt(k,186)*y(k,73) + &
                 rxt(k,204)*y(k,92) +rxt(k,220)*y(k,75) +rxt(k,222)*y(k,89) + &
                 rxt(k,239)*y(k,88))*y(k,68) + (rxt(k,95) +rxt(k,272)*y(k,52) + &
                 rxt(k,274)*y(k,51) +rxt(k,325)*y(k,28) +rxt(k,326)*y(k,32))*y(k,101) &
                  + (rxt(k,382) +rxt(k,390) +2.000_r8*rxt(k,334)*y(k,103) + &
                 2.000_r8*rxt(k,392)*y(k,55))*y(k,96) + (rxt(k,321)*y(k,69) + &
                 rxt(k,322)*y(k,110) +rxt(k,323)*y(k,110))*y(k,100) + (rxt(k,97) + &
                 rxt(k,330)*y(k,28))*y(k,102) + (rxt(k,332)*y(k,110) + &
                 2.000_r8*rxt(k,375)*y(k,69))*y(k,103) +rxt(k,484)*y(k,2)*y(k,2) &
                  +rxt(k,418)*y(k,104)*y(k,44) +rxt(k,424)*y(k,52)*y(k,48) &
                  +rxt(k,438)*y(k,95)*y(k,49) +rxt(k,20)*y(k,53) +rxt(k,383)*y(k,97)
         loss(k,75) = (rxt(k,409)* y(k,32) +rxt(k,408)* y(k,43) +rxt(k,435)* y(k,51) &
                  +rxt(k,427)* y(k,52) +rxt(k,398)* y(k,54) +rxt(k,483)* y(k,59) &
                  +rxt(k,451)* y(k,60) +rxt(k,127)* y(k,68) +rxt(k,129)* y(k,69) &
                  +rxt(k,287)* y(k,81) +rxt(k,297)* y(k,83) +rxt(k,389)* y(k,95) &
                  +rxt(k,318)* y(k,98) +rxt(k,327)* y(k,101) +rxt(k,414)* y(k,104) &
                  +rxt(k,343)* y(k,105) +rxt(k,358)* y(k,106) + rxt(k,7) + rxt(k,8) &
                  + het_rates(k,56))* y(k,56)
         prod(k,75) = (rxt(k,319)*y(k,98) +rxt(k,333)*y(k,103) +rxt(k,400)*y(k,55)) &
                 *y(k,54) + (rxt(k,96) +rxt(k,273)*y(k,52))*y(k,101) &
                  +rxt(k,356)*y(k,106)*y(k,96)
         loss(k,3) = ( + rxt(k,60) + het_rates(k,57))* y(k,57)
         prod(k,3) = (rxt(k,455)*y(k,24) +rxt(k,485)*y(k,2))*y(k,24)
         loss(k,44) = (rxt(k,481)* y(k,14) +rxt(k,482)* y(k,43) +rxt(k,483)* y(k,56) &
                  + het_rates(k,59))* y(k,59)
         prod(k,44) = (rxt(k,28) +2.000_r8*rxt(k,484)*y(k,2) +rxt(k,485)*y(k,24) + &
                 rxt(k,486)*y(k,24) +rxt(k,489)*y(k,51) +rxt(k,492)*y(k,54) + &
                 rxt(k,493)*y(k,104))*y(k,2) + (rxt(k,471)*y(k,6) +rxt(k,497)*y(k,7) + &
                 3.000_r8*rxt(k,498)*y(k,21) +2.000_r8*rxt(k,499)*y(k,34) + &
                 2.000_r8*rxt(k,520)*y(k,13) +rxt(k,521)*y(k,15) +rxt(k,500)*y(k,36)) &
                 *y(k,95) + (2.000_r8*rxt(k,509)*y(k,13) +rxt(k,511)*y(k,15) + &
                 3.000_r8*rxt(k,516)*y(k,21) +rxt(k,495)*y(k,36))*y(k,104) &
                  + (2.000_r8*rxt(k,508)*y(k,13) +rxt(k,510)*y(k,15) + &
                 3.000_r8*rxt(k,515)*y(k,21))*y(k,60) + (rxt(k,52) + &
                 rxt(k,494)*y(k,54))*y(k,36) +rxt(k,27)*y(k,1) +rxt(k,30)*y(k,3) &
                  +rxt(k,58)*y(k,45)
         loss(k,78) = (rxt(k,508)* y(k,13) +rxt(k,445)* y(k,14) +rxt(k,510)* y(k,15) &
                  +rxt(k,513)* y(k,17) +rxt(k,446)* y(k,20) +rxt(k,515)* y(k,21) &
                  +rxt(k,458)* y(k,25) +rxt(k,447)* y(k,33) +rxt(k,448)* y(k,35) &
                  + (rxt(k,449) +rxt(k,450))* y(k,43) +rxt(k,467)* y(k,46) +rxt(k,451) &
                 * y(k,56) + (rxt(k,109) +rxt(k,110))* y(k,65) +rxt(k,122)* y(k,68) &
                  +rxt(k,280)* y(k,81) +rxt(k,309)* y(k,98) +rxt(k,336)* y(k,105) &
                  +rxt(k,345)* y(k,106) + het_rates(k,60))* y(k,60)
         prod(k,78) = (4.000_r8*rxt(k,470)*y(k,5) +rxt(k,471)*y(k,6) + &
                 2.000_r8*rxt(k,472)*y(k,8) +2.000_r8*rxt(k,473)*y(k,9) + &
                 2.000_r8*rxt(k,474)*y(k,10) +rxt(k,475)*y(k,11) + &
                 2.000_r8*rxt(k,476)*y(k,12) +rxt(k,522)*y(k,37) +rxt(k,523)*y(k,38) + &
                 rxt(k,524)*y(k,39) +rxt(k,477)*y(k,40) +rxt(k,507)*y(k,30))*y(k,95) &
                  + (rxt(k,46) +rxt(k,452)*y(k,18) +2.000_r8*rxt(k,453)*y(k,24) + &
                 rxt(k,455)*y(k,24) +rxt(k,457)*y(k,51) +rxt(k,462)*y(k,54) + &
                 rxt(k,463)*y(k,104) +rxt(k,486)*y(k,2))*y(k,24) &
                  + (rxt(k,105)*y(k,52) +rxt(k,142)*y(k,75) +rxt(k,149)*y(k,99) + &
                 rxt(k,172)*y(k,77) +rxt(k,190)*y(k,73) +rxt(k,207)*y(k,92) + &
                 rxt(k,225)*y(k,89) +rxt(k,243)*y(k,88))*y(k,61) &
                  + (rxt(k,154)*y(k,99) +rxt(k,159)*y(k,75) +rxt(k,177)*y(k,77) + &
                 rxt(k,194)*y(k,73) +rxt(k,212)*y(k,92) +rxt(k,229)*y(k,89) + &
                 rxt(k,247)*y(k,88))*y(k,62) + (rxt(k,164)*y(k,77) + &
                 rxt(k,176)*y(k,75) +rxt(k,182)*y(k,73) +rxt(k,200)*y(k,92) + &
                 rxt(k,217)*y(k,89) +rxt(k,235)*y(k,88) +rxt(k,252)*y(k,99))*y(k,63) &
                  + (3.000_r8*rxt(k,512)*y(k,16) +rxt(k,514)*y(k,17) + &
                 rxt(k,517)*y(k,37) +rxt(k,518)*y(k,38) +rxt(k,519)*y(k,39) + &
                 rxt(k,466)*y(k,40))*y(k,104) + (rxt(k,56) +rxt(k,465)*y(k,54)) &
                 *y(k,40) +rxt(k,27)*y(k,1) +2.000_r8*rxt(k,44)*y(k,22) &
                  +2.000_r8*rxt(k,45)*y(k,23) +rxt(k,47)*y(k,25) +rxt(k,50)*y(k,30) &
                  +rxt(k,59)*y(k,46) +rxt(k,106)*y(k,64)*y(k,51)
         loss(k,58) = (rxt(k,100)* y(k,32) +rxt(k,103)* y(k,40) +rxt(k,104)* y(k,42) &
                  +rxt(k,105)* y(k,52) +rxt(k,190)* y(k,73) +rxt(k,142)* y(k,75) &
                  +rxt(k,172)* y(k,77) +rxt(k,243)* y(k,88) +rxt(k,225)* y(k,89) &
                  +rxt(k,207)* y(k,92) +rxt(k,149)* y(k,99) +rxt(k,102)* y(k,110) &
                  + het_rates(k,61))* y(k,61)
         prod(k,58) = (rxt(k,125)*y(k,68) +rxt(k,284)*y(k,81) +rxt(k,293)*y(k,83) + &
                 rxt(k,313)*y(k,98) +rxt(k,340)*y(k,105) +rxt(k,351)*y(k,106))*y(k,40) &
                  + (rxt(k,109)*y(k,65) +rxt(k,122)*y(k,68) +rxt(k,280)*y(k,81) + &
                 rxt(k,309)*y(k,98) +rxt(k,336)*y(k,105) +rxt(k,345)*y(k,106))*y(k,60) &
                  + (rxt(k,111)*y(k,65) +rxt(k,281)*y(k,81) +rxt(k,346)*y(k,106)) &
                 *y(k,24) + (rxt(k,107)*y(k,51) +rxt(k,108)*y(k,54))*y(k,64) &
                  +rxt(k,380)*y(k,62) +rxt(k,381)*y(k,63)
         loss(k,46) = (rxt(k,101)* y(k,40) +rxt(k,194)* y(k,73) +rxt(k,159)* y(k,75) &
                  +rxt(k,177)* y(k,77) +rxt(k,247)* y(k,88) +rxt(k,229)* y(k,89) &
                  +rxt(k,212)* y(k,92) +rxt(k,154)* y(k,99) + rxt(k,380) &
                  + het_rates(k,62))* y(k,62)
         prod(k,46) =rxt(k,102)*y(k,110)*y(k,61)
         loss(k,45) = (rxt(k,182)* y(k,73) + (rxt(k,176) +rxt(k,262))* y(k,75) &
                  + (rxt(k,164) +rxt(k,263))* y(k,77) +rxt(k,235)* y(k,88) +rxt(k,217) &
                 * y(k,89) +rxt(k,200)* y(k,92) +rxt(k,252)* y(k,99) + rxt(k,381) &
                  + het_rates(k,63))* y(k,63)
         prod(k,45) = (rxt(k,101)*y(k,62) +rxt(k,103)*y(k,61))*y(k,40)
         loss(k,36) = ((rxt(k,106) +rxt(k,107))* y(k,51) +rxt(k,108)* y(k,54) &
                  + het_rates(k,64))* y(k,64)
         prod(k,36) = (rxt(k,123)*y(k,68) +rxt(k,310)*y(k,98) +rxt(k,337)*y(k,105)) &
                 *y(k,24) +rxt(k,110)*y(k,65)*y(k,60)
         loss(k,64) = (rxt(k,111)* y(k,24) +rxt(k,112)* y(k,32) +rxt(k,119)* y(k,42) &
                  +rxt(k,275)* y(k,51) +rxt(k,276)* y(k,52) +rxt(k,120)* y(k,54) &
                  +rxt(k,121)* y(k,55) + (rxt(k,109) +rxt(k,110))* y(k,60) +rxt(k,181) &
                 * y(k,73) + (rxt(k,165) +rxt(k,253))* y(k,75) + (rxt(k,163) + &
                 rxt(k,260))* y(k,77) +rxt(k,234)* y(k,88) +rxt(k,216)* y(k,89) &
                  +rxt(k,199)* y(k,92) +rxt(k,251)* y(k,99) +rxt(k,114)* y(k,110) &
                  + rxt(k,63) + het_rates(k,65))* y(k,65)
         prod(k,64) = (rxt(k,325)*y(k,101) +rxt(k,347)*y(k,106))*y(k,28) &
                  + (rxt(k,64) +rxt(k,278))*y(k,67) + (rxt(k,124)*y(k,32) + &
                 rxt(k,126)*y(k,54))*y(k,68)
         loss(k,43) = (rxt(k,189)* y(k,73) + (rxt(k,242) +rxt(k,266))* y(k,75) &
                  + (rxt(k,171) +rxt(k,254))* y(k,77) +rxt(k,241)* y(k,88) +rxt(k,224) &
                 * y(k,89) +rxt(k,206)* y(k,92) +rxt(k,148)* y(k,99) + rxt(k,279) &
                  + het_rates(k,66))* y(k,66)
         prod(k,43) =rxt(k,113)*y(k,110)*y(k,67)
         loss(k,54) = ((rxt(k,117) +rxt(k,118))* y(k,51) + (rxt(k,115) +rxt(k,116)) &
                 * y(k,52) +rxt(k,191)* y(k,73) + (rxt(k,153) +rxt(k,255))* y(k,75) &
                  + (rxt(k,173) +rxt(k,256))* y(k,77) +rxt(k,244)* y(k,88) +rxt(k,226) &
                 * y(k,89) +rxt(k,208)* y(k,92) +rxt(k,150)* y(k,99) +rxt(k,113) &
                 * y(k,110) + rxt(k,64) + rxt(k,278) + het_rates(k,67))* y(k,67)
         prod(k,54) =rxt(k,114)*y(k,110)*y(k,65) +rxt(k,279)*y(k,66)
         loss(k,59) = (rxt(k,123)* y(k,24) +rxt(k,124)* y(k,32) +rxt(k,125)* y(k,40) &
                  +rxt(k,126)* y(k,54) +rxt(k,127)* y(k,56) +rxt(k,122)* y(k,60) &
                  +rxt(k,186)* y(k,73) +rxt(k,220)* y(k,75) +rxt(k,169)* y(k,77) &
                  +rxt(k,239)* y(k,88) +rxt(k,222)* y(k,89) +rxt(k,204)* y(k,92) &
                  +rxt(k,146)* y(k,99) + rxt(k,65) + het_rates(k,68))* y(k,68)
         prod(k,59) = (rxt(k,311)*y(k,98) +rxt(k,330)*y(k,102))*y(k,28)
         loss(k,68) = ((rxt(k,128) +rxt(k,130))* y(k,55) +rxt(k,129)* y(k,56) &
                  +rxt(k,133)* y(k,70) +rxt(k,377)* y(k,72) +rxt(k,378)* y(k,73) &
                  +rxt(k,136)* y(k,75) +rxt(k,139)* y(k,77) +rxt(k,376)* y(k,78) &
                  +rxt(k,563)* y(k,80) +rxt(k,561)* y(k,88) +rxt(k,298)* y(k,89) &
                  +rxt(k,299)* y(k,90) +rxt(k,301)* y(k,91) +rxt(k,303)* y(k,92) &
                  +rxt(k,379)* y(k,93) +rxt(k,562)* y(k,99) +rxt(k,321)* y(k,100) &
                  +rxt(k,375)* y(k,103) + het_rates(k,69))* y(k,69)
         prod(k,68) = (rxt(k,77) +rxt(k,78) +rxt(k,79) +rxt(k,80) +rxt(k,81) + &
                 rxt(k,82) +rxt(k,319)*y(k,98) +rxt(k,328)*y(k,101) + &
                 rxt(k,342)*y(k,105) +rxt(k,355)*y(k,106))*y(k,54) + (rxt(k,83) + &
                 rxt(k,85) +rxt(k,86) +rxt(k,87) +rxt(k,88) +rxt(k,90) +rxt(k,91) + &
                 rxt(k,92))*y(k,55) + (rxt(k,99) +rxt(k,353) +rxt(k,131)*y(k,51) + &
                 rxt(k,348)*y(k,33) +rxt(k,356)*y(k,96))*y(k,106) + (rxt(k,93) + &
                 rxt(k,132) +rxt(k,312)*y(k,32) +rxt(k,316)*y(k,96))*y(k,98) &
                  + (rxt(k,100)*y(k,61) +rxt(k,339)*y(k,105))*y(k,32) + (rxt(k,96) + &
                 rxt(k,327)*y(k,56))*y(k,101) +rxt(k,66)*y(k,48) +rxt(k,16)*y(k,51) &
                  +rxt(k,75)*y(k,81) +rxt(k,76)*y(k,83) +rxt(k,98)*y(k,105)
         loss(k,11) = (rxt(k,133)* y(k,69) +rxt(k,134)* y(k,110) + het_rates(k,70)) &
                 * y(k,70)
         prod(k,11) =rxt(k,322)*y(k,110)*y(k,100)
         loss(k,42) = (rxt(k,184)* y(k,73) +rxt(k,198)* y(k,75) +rxt(k,167)* y(k,77) &
                  +rxt(k,237)* y(k,88) +rxt(k,219)* y(k,89) +rxt(k,202)* y(k,92) &
                  +rxt(k,144)* y(k,99) + het_rates(k,71))* y(k,71)
         prod(k,42) =rxt(k,338)*y(k,105)*y(k,28)
         loss(k,22) = (rxt(k,377)* y(k,69) +rxt(k,369)* y(k,110) + rxt(k,368) &
                  + het_rates(k,72))* y(k,72)
         prod(k,22) = (rxt(k,134)*y(k,70) +rxt(k,367)*y(k,78))*y(k,110) +rxt(k,370) &
                 *y(k,73)
         loss(k,65) = (rxt(k,190)* y(k,61) +rxt(k,194)* y(k,62) +rxt(k,182)* y(k,63) &
                  +rxt(k,181)* y(k,65) +rxt(k,189)* y(k,66) +rxt(k,191)* y(k,67) &
                  +rxt(k,186)* y(k,68) +rxt(k,378)* y(k,69) +rxt(k,184)* y(k,71) &
                  +rxt(k,196)* y(k,81) +rxt(k,192)* y(k,82) +rxt(k,183)* y(k,83) &
                  +rxt(k,195)* y(k,84) +rxt(k,188)* y(k,85) +rxt(k,193)* y(k,86) &
                  +rxt(k,180)* y(k,87) +rxt(k,185)* y(k,98) +rxt(k,371)* y(k,110) &
                  + rxt(k,370) + het_rates(k,73))* y(k,73)
         prod(k,65) = (rxt(k,300)*y(k,90) +rxt(k,369)*y(k,72))*y(k,110) +rxt(k,372) &
                 *y(k,75)
         loss(k,8) = (rxt(k,135)* y(k,110) + het_rates(k,74))* y(k,74)
         prod(k,8) =rxt(k,137)*y(k,75)*y(k,50)
         loss(k,80) = (rxt(k,137)* y(k,50) +rxt(k,142)* y(k,61) +rxt(k,159)* y(k,62) &
                  + (rxt(k,176) +rxt(k,262))* y(k,63) + (rxt(k,165) +rxt(k,253)) &
                 * y(k,65) + (rxt(k,242) +rxt(k,266))* y(k,66) + (rxt(k,153) + &
                 rxt(k,255))* y(k,67) +rxt(k,220)* y(k,68) +rxt(k,136)* y(k,69) &
                  +rxt(k,198)* y(k,71) +rxt(k,161)* y(k,81) +rxt(k,157)* y(k,82) &
                  + (rxt(k,187) +rxt(k,259))* y(k,83) +rxt(k,160)* y(k,84) &
                  + (rxt(k,231) +rxt(k,257))* y(k,85) +rxt(k,158)* y(k,86) &
                  + (rxt(k,141) +rxt(k,264))* y(k,87) +rxt(k,209)* y(k,98) +rxt(k,373) &
                 * y(k,110) + rxt(k,372) + het_rates(k,75))* y(k,75)
         prod(k,80) = (rxt(k,135)*y(k,74) +rxt(k,371)*y(k,73))*y(k,110) +rxt(k,374) &
                 *y(k,77)
         loss(k,9) = (rxt(k,138)* y(k,110) + het_rates(k,76))* y(k,76)
         prod(k,9) =rxt(k,140)*y(k,77)*y(k,50)
         loss(k,82) = (rxt(k,140)* y(k,50) +rxt(k,172)* y(k,61) +rxt(k,177)* y(k,62) &
                  + (rxt(k,164) +rxt(k,263))* y(k,63) + (rxt(k,163) +rxt(k,260)) &
                 * y(k,65) + (rxt(k,171) +rxt(k,254))* y(k,66) + (rxt(k,173) + &
                 rxt(k,256))* y(k,67) +rxt(k,169)* y(k,68) +rxt(k,139)* y(k,69) &
                  +rxt(k,167)* y(k,71) +rxt(k,179)* y(k,81) +rxt(k,174)* y(k,82) &
                  + (rxt(k,166) +rxt(k,261))* y(k,83) +rxt(k,178)* y(k,84) &
                  + (rxt(k,170) +rxt(k,258))* y(k,85) +rxt(k,175)* y(k,86) &
                  + (rxt(k,162) +rxt(k,265))* y(k,87) +rxt(k,168)* y(k,98) &
                  + rxt(k,374) + het_rates(k,77))* y(k,77)
         prod(k,82) = (rxt(k,138)*y(k,76) +rxt(k,373)*y(k,75))*y(k,110)
         loss(k,28) = (rxt(k,376)* y(k,69) +rxt(k,367)* y(k,110) + het_rates(k,78)) &
                 * y(k,78)
         prod(k,28) = (rxt(k,304)*y(k,32) +rxt(k,305)*y(k,43) +rxt(k,306)*y(k,104)) &
                 *y(k,92) +rxt(k,368)*y(k,72) +rxt(k,323)*y(k,110)*y(k,100)
         loss(k,17) = (rxt(k,419)* y(k,54) +rxt(k,420)* y(k,55) + het_rates(k,79)) &
                 * y(k,79)
         prod(k,17) = (.800_r8*rxt(k,561)*y(k,88) +.900_r8*rxt(k,563)*y(k,80))*y(k,69) &
                  +rxt(k,565)*y(k,80)*y(k,54)
         loss(k,23) = ((rxt(k,565) +rxt(k,566))* y(k,54) +rxt(k,564)* y(k,55) &
                  +rxt(k,563)* y(k,69) + het_rates(k,80))* y(k,80)
         prod(k,23) =rxt(k,577)*y(k,108) +rxt(k,582)*y(k,109)
         loss(k,62) = (rxt(k,281)* y(k,24) +rxt(k,282)* y(k,32) +rxt(k,284)* y(k,40) &
                  +rxt(k,285)* y(k,42) +rxt(k,286)* y(k,52) +rxt(k,287)* y(k,56) &
                  +rxt(k,280)* y(k,60) +rxt(k,196)* y(k,73) +rxt(k,161)* y(k,75) &
                  +rxt(k,179)* y(k,77) +rxt(k,249)* y(k,88) +rxt(k,232)* y(k,89) &
                  +rxt(k,214)* y(k,92) +rxt(k,156)* y(k,99) +rxt(k,283)* y(k,110) &
                  + rxt(k,75) + het_rates(k,81))* y(k,81)
         prod(k,62) = (rxt(k,105)*y(k,61) +rxt(k,273)*y(k,101) +rxt(k,315)*y(k,98) + &
                 rxt(k,341)*y(k,105) +rxt(k,354)*y(k,106))*y(k,52) &
                  + (rxt(k,106)*y(k,64) +rxt(k,118)*y(k,67) +rxt(k,274)*y(k,101) + &
                 rxt(k,275)*y(k,65))*y(k,51) + (rxt(k,296)*y(k,54) + &
                 rxt(k,297)*y(k,56))*y(k,83) +rxt(k,268)*y(k,82)
         loss(k,47) = (rxt(k,192)* y(k,73) +rxt(k,157)* y(k,75) +rxt(k,174)* y(k,77) &
                  +rxt(k,245)* y(k,88) +rxt(k,227)* y(k,89) +rxt(k,210)* y(k,92) &
                  +rxt(k,151)* y(k,99) + rxt(k,268) + het_rates(k,82))* y(k,82)
         prod(k,47) =rxt(k,117)*y(k,67)*y(k,51) +rxt(k,283)*y(k,110)*y(k,81)
         loss(k,60) = ((rxt(k,271) +rxt(k,293))* y(k,40) +rxt(k,295)* y(k,42) &
                  +rxt(k,296)* y(k,54) +rxt(k,297)* y(k,56) +rxt(k,183)* y(k,73) &
                  + (rxt(k,187) +rxt(k,259))* y(k,75) + (rxt(k,166) +rxt(k,261)) &
                 * y(k,77) +rxt(k,236)* y(k,88) +rxt(k,218)* y(k,89) +rxt(k,201) &
                 * y(k,92) +rxt(k,143)* y(k,99) +rxt(k,291)* y(k,110) + rxt(k,76) &
                  + het_rates(k,83))* y(k,83)
         prod(k,60) = (rxt(k,104)*y(k,61) +rxt(k,119)*y(k,65) +rxt(k,285)*y(k,81) + &
                 rxt(k,314)*y(k,98) +rxt(k,352)*y(k,106))*y(k,42) &
                  + (rxt(k,115)*y(k,67) +rxt(k,272)*y(k,101) +rxt(k,276)*y(k,65) + &
                 rxt(k,286)*y(k,81))*y(k,52) +rxt(k,267)*y(k,101)*y(k,51) &
                  +rxt(k,287)*y(k,81)*y(k,56) +rxt(k,277)*y(k,85) +rxt(k,270)*y(k,87)
         loss(k,53) = (rxt(k,288)* y(k,50) +rxt(k,195)* y(k,73) +rxt(k,160)* y(k,75) &
                  +rxt(k,178)* y(k,77) +rxt(k,248)* y(k,88) +rxt(k,230)* y(k,89) &
                  +rxt(k,213)* y(k,92) +rxt(k,155)* y(k,99) + rxt(k,269) &
                  + het_rates(k,84))* y(k,84)
         prod(k,53) =rxt(k,289)*y(k,110)*y(k,85)
         loss(k,55) = (rxt(k,290)* y(k,42) +rxt(k,292)* y(k,50) +rxt(k,188)* y(k,73) &
                  + (rxt(k,231) +rxt(k,257))* y(k,75) + (rxt(k,170) +rxt(k,258)) &
                 * y(k,77) +rxt(k,240)* y(k,88) +rxt(k,223)* y(k,89) +rxt(k,205) &
                 * y(k,92) +rxt(k,147)* y(k,99) +rxt(k,289)* y(k,110) + rxt(k,277) &
                  + het_rates(k,85))* y(k,85)
         prod(k,55) =rxt(k,116)*y(k,67)*y(k,52) +rxt(k,291)*y(k,110)*y(k,83) &
                  +rxt(k,269)*y(k,84)
         loss(k,48) = (rxt(k,294)* y(k,42) +rxt(k,193)* y(k,73) +rxt(k,158)* y(k,75) &
                  +rxt(k,175)* y(k,77) +rxt(k,246)* y(k,88) +rxt(k,228)* y(k,89) &
                  +rxt(k,211)* y(k,92) +rxt(k,152)* y(k,99) + het_rates(k,86)) &
                 * y(k,86)
         prod(k,48) =rxt(k,271)*y(k,83)*y(k,40)
         loss(k,51) = (rxt(k,180)* y(k,73) + (rxt(k,141) +rxt(k,264))* y(k,75) &
                  + (rxt(k,162) +rxt(k,265))* y(k,77) +rxt(k,233)* y(k,88) +rxt(k,215) &
                 * y(k,89) +rxt(k,197)* y(k,92) +rxt(k,250)* y(k,99) + rxt(k,270) &
                  + het_rates(k,87))* y(k,87)
         prod(k,51) = (rxt(k,290)*y(k,85) +rxt(k,294)*y(k,86) +rxt(k,295)*y(k,83)) &
                 *y(k,42) + (rxt(k,288)*y(k,84) +rxt(k,292)*y(k,85))*y(k,50)
         loss(k,66) = (rxt(k,363)* y(k,28) +rxt(k,243)* y(k,61) +rxt(k,247)* y(k,62) &
                  +rxt(k,235)* y(k,63) +rxt(k,234)* y(k,65) +rxt(k,241)* y(k,66) &
                  +rxt(k,244)* y(k,67) +rxt(k,239)* y(k,68) +rxt(k,561)* y(k,69) &
                  +rxt(k,237)* y(k,71) +rxt(k,249)* y(k,81) +rxt(k,245)* y(k,82) &
                  +rxt(k,236)* y(k,83) +rxt(k,248)* y(k,84) +rxt(k,240)* y(k,85) &
                  +rxt(k,246)* y(k,86) +rxt(k,233)* y(k,87) +rxt(k,238)* y(k,98) &
                  +rxt(k,360)* y(k,110) + rxt(k,365) + het_rates(k,88))* y(k,88)
         prod(k,66) = (rxt(k,571) +rxt(k,570)*y(k,48) +rxt(k,572)*y(k,51))*y(k,99) &
                  +rxt(k,16)*y(k,51) +rxt(k,565)*y(k,80)*y(k,54) +rxt(k,569)*y(k,94) &
                 *y(k,55) +rxt(k,364)*y(k,91) +rxt(k,366)*y(k,93) +rxt(k,574)*y(k,107)
         loss(k,67) = (rxt(k,225)* y(k,61) +rxt(k,229)* y(k,62) +rxt(k,217)* y(k,63) &
                  +rxt(k,216)* y(k,65) +rxt(k,224)* y(k,66) +rxt(k,226)* y(k,67) &
                  +rxt(k,222)* y(k,68) +rxt(k,298)* y(k,69) +rxt(k,219)* y(k,71) &
                  +rxt(k,232)* y(k,81) +rxt(k,227)* y(k,82) +rxt(k,218)* y(k,83) &
                  +rxt(k,230)* y(k,84) +rxt(k,223)* y(k,85) +rxt(k,228)* y(k,86) &
                  +rxt(k,215)* y(k,87) +rxt(k,221)* y(k,98) +rxt(k,362)* y(k,110) &
                  + het_rates(k,89))* y(k,89)
         prod(k,67) =rxt(k,361)*y(k,110)*y(k,92)
         loss(k,13) = (rxt(k,299)* y(k,69) +rxt(k,300)* y(k,110) + het_rates(k,90)) &
                 * y(k,90)
         prod(k,13) =rxt(k,362)*y(k,110)*y(k,89)
         loss(k,25) = (rxt(k,301)* y(k,69) +rxt(k,302)* y(k,110) + rxt(k,364) &
                  + het_rates(k,91))* y(k,91)
         prod(k,25) = (rxt(k,307)*y(k,93) +rxt(k,363)*y(k,88))*y(k,28)
         loss(k,70) = (rxt(k,304)* y(k,32) +rxt(k,305)* y(k,43) +rxt(k,207)* y(k,61) &
                  +rxt(k,212)* y(k,62) +rxt(k,200)* y(k,63) +rxt(k,199)* y(k,65) &
                  +rxt(k,206)* y(k,66) +rxt(k,208)* y(k,67) +rxt(k,204)* y(k,68) &
                  +rxt(k,303)* y(k,69) +rxt(k,202)* y(k,71) +rxt(k,214)* y(k,81) &
                  +rxt(k,210)* y(k,82) +rxt(k,201)* y(k,83) +rxt(k,213)* y(k,84) &
                  +rxt(k,205)* y(k,85) +rxt(k,211)* y(k,86) +rxt(k,197)* y(k,87) &
                  +rxt(k,203)* y(k,98) +rxt(k,306)* y(k,104) +rxt(k,361)* y(k,110) &
                  + het_rates(k,92))* y(k,92)
         prod(k,70) = (rxt(k,302)*y(k,91) +rxt(k,308)*y(k,93) +rxt(k,360)*y(k,88)) &
                 *y(k,110)
         loss(k,24) = (rxt(k,307)* y(k,28) +rxt(k,379)* y(k,69) +rxt(k,308)* y(k,110) &
                  + rxt(k,366) + het_rates(k,93))* y(k,93)
         prod(k,24) =rxt(k,365)*y(k,88)
         loss(k,20) = (rxt(k,567)* y(k,54) + (rxt(k,568) +rxt(k,569))* y(k,55) &
                  + het_rates(k,94))* y(k,94)
         prod(k,20) =rxt(k,66)*y(k,48) +rxt(k,583)*y(k,109)
         loss(k,57) = (rxt(k,470)* y(k,5) +rxt(k,471)* y(k,6) +rxt(k,497)* y(k,7) &
                  +rxt(k,472)* y(k,8) +rxt(k,473)* y(k,9) +rxt(k,474)* y(k,10) &
                  +rxt(k,475)* y(k,11) +rxt(k,476)* y(k,12) +rxt(k,520)* y(k,13) &
                  +rxt(k,521)* y(k,15) + (rxt(k,533) +rxt(k,534) +rxt(k,535))* y(k,20) &
                  +rxt(k,498)* y(k,21) +rxt(k,506)* y(k,29) +rxt(k,507)* y(k,30) &
                  +rxt(k,384)* y(k,33) +rxt(k,499)* y(k,34) + (rxt(k,500) +rxt(k,501)) &
                 * y(k,36) +rxt(k,522)* y(k,37) +rxt(k,523)* y(k,38) +rxt(k,524) &
                 * y(k,39) + (rxt(k,477) +rxt(k,478))* y(k,40) + (rxt(k,437) + &
                 rxt(k,438))* y(k,49) + (rxt(k,387) +rxt(k,388))* y(k,55) +rxt(k,389) &
                 * y(k,56) +rxt(k,385)* y(k,110) + rxt(k,386) + het_rates(k,95)) &
                 * y(k,95)
         prod(k,57) = (rxt(k,6) +rxt(k,420)*y(k,79))*y(k,55) +rxt(k,7)*y(k,56) &
                  +.850_r8*rxt(k,562)*y(k,99)*y(k,69) +rxt(k,1)*y(k,110)
         loss(k,37) = (rxt(k,391)* y(k,54) +rxt(k,392)* y(k,55) +rxt(k,316)* y(k,98) &
                  +rxt(k,334)* y(k,103) +rxt(k,356)* y(k,106) + rxt(k,382) &
                  + rxt(k,390) + het_rates(k,96))* y(k,96)
         prod(k,37) = (rxt(k,394) +rxt(k,393)*y(k,28) +rxt(k,395)*y(k,54) + &
                 rxt(k,396)*y(k,55) +rxt(k,397)*y(k,56))*y(k,97) +rxt(k,7)*y(k,56)
         loss(k,10) = (rxt(k,393)* y(k,28) +rxt(k,395)* y(k,54) +rxt(k,396)* y(k,55) &
                  +rxt(k,397)* y(k,56) + rxt(k,383) + rxt(k,394) + het_rates(k,97)) &
                 * y(k,97)
         prod(k,10) =rxt(k,387)*y(k,95)*y(k,55)
         loss(k,76) = (rxt(k,310)* y(k,24) +rxt(k,311)* y(k,28) +rxt(k,312)* y(k,32) &
                  +rxt(k,313)* y(k,40) +rxt(k,314)* y(k,42) +rxt(k,315)* y(k,52) &
                  + (rxt(k,319) +rxt(k,320))* y(k,54) +rxt(k,317)* y(k,55) +rxt(k,318) &
                 * y(k,56) +rxt(k,309)* y(k,60) +rxt(k,185)* y(k,73) +rxt(k,209) &
                 * y(k,75) +rxt(k,168)* y(k,77) +rxt(k,238)* y(k,88) +rxt(k,221) &
                 * y(k,89) +rxt(k,203)* y(k,92) +rxt(k,316)* y(k,96) +rxt(k,145) &
                 * y(k,99) + rxt(k,93) + rxt(k,132) + het_rates(k,98))* y(k,98)
         prod(k,76) = (rxt(k,120)*y(k,65) +rxt(k,329)*y(k,101))*y(k,54) &
                  + (rxt(k,128)*y(k,69) +rxt(k,130)*y(k,69))*y(k,55) +rxt(k,65) &
                 *y(k,68) +rxt(k,97)*y(k,102)
         loss(k,77) = (rxt(k,570)* y(k,48) +rxt(k,572)* y(k,51) +rxt(k,359)* y(k,55) &
                  +rxt(k,149)* y(k,61) +rxt(k,154)* y(k,62) +rxt(k,252)* y(k,63) &
                  +rxt(k,251)* y(k,65) +rxt(k,148)* y(k,66) +rxt(k,150)* y(k,67) &
                  +rxt(k,146)* y(k,68) +rxt(k,562)* y(k,69) +rxt(k,144)* y(k,71) &
                  +rxt(k,156)* y(k,81) +rxt(k,151)* y(k,82) +rxt(k,143)* y(k,83) &
                  +rxt(k,155)* y(k,84) +rxt(k,147)* y(k,85) +rxt(k,152)* y(k,86) &
                  +rxt(k,250)* y(k,87) +rxt(k,145)* y(k,98) +rxt(k,324)* y(k,110) &
                  + rxt(k,571) + het_rates(k,99))* y(k,99)
         prod(k,77) = (rxt(k,86) +rxt(k,88) +rxt(k,564)*y(k,80) +rxt(k,568)*y(k,94) + &
                 rxt(k,575)*y(k,107) +rxt(k,579)*y(k,108))*y(k,55) &
                  + (rxt(k,333)*y(k,54) +rxt(k,334)*y(k,96))*y(k,103) &
                  +rxt(k,573)*y(k,107)*y(k,28) +2.000_r8*rxt(k,145)*y(k,99)*y(k,98) &
                  +rxt(k,94)*y(k,100)
         loss(k,29) = (rxt(k,321)* y(k,69) + (rxt(k,322) +rxt(k,323))* y(k,110) &
                  + rxt(k,94) + het_rates(k,100))* y(k,100)
         prod(k,29) = (rxt(k,324)*y(k,99) +rxt(k,332)*y(k,103))*y(k,110)
         loss(k,49) = (rxt(k,325)* y(k,28) +rxt(k,326)* y(k,32) + (rxt(k,267) + &
                 rxt(k,274))* y(k,51) + (rxt(k,272) +rxt(k,273))* y(k,52) &
                  + (rxt(k,328) +rxt(k,329))* y(k,54) +rxt(k,327)* y(k,56) + rxt(k,95) &
                  + rxt(k,96) + het_rates(k,101))* y(k,101)
         prod(k,49) = (rxt(k,127)*y(k,68) +rxt(k,318)*y(k,98) +rxt(k,343)*y(k,105) + &
                 rxt(k,358)*y(k,106))*y(k,56) + (rxt(k,121)*y(k,65) + &
                 rxt(k,357)*y(k,106))*y(k,55) +rxt(k,331)*y(k,102)*y(k,54)
         loss(k,21) = (rxt(k,330)* y(k,28) +rxt(k,331)* y(k,54) + rxt(k,97) &
                  + het_rates(k,102))* y(k,102)
         prod(k,21) =rxt(k,317)*y(k,98)*y(k,55)
         loss(k,41) = (rxt(k,333)* y(k,54) +rxt(k,375)* y(k,69) +rxt(k,334)* y(k,96) &
                  +rxt(k,332)* y(k,110) + het_rates(k,103))* y(k,103)
         prod(k,41) =rxt(k,359)*y(k,99)*y(k,55)
         loss(k,63) = (rxt(k,493)* y(k,2) +rxt(k,509)* y(k,13) +rxt(k,527)* y(k,14) &
                  +rxt(k,511)* y(k,15) +rxt(k,512)* y(k,16) +rxt(k,514)* y(k,17) &
                  +rxt(k,530)* y(k,19) +rxt(k,531)* y(k,20) +rxt(k,516)* y(k,21) &
                  + (rxt(k,463) +rxt(k,464))* y(k,24) +rxt(k,461)* y(k,25) &
                  + (rxt(k,532) +rxt(k,536))* y(k,27) +rxt(k,410)* y(k,33) +rxt(k,411) &
                 * y(k,35) +rxt(k,495)* y(k,36) +rxt(k,517)* y(k,37) +rxt(k,518) &
                 * y(k,38) +rxt(k,519)* y(k,39) +rxt(k,466)* y(k,40) +rxt(k,442) &
                 * y(k,42) +rxt(k,412)* y(k,43) +rxt(k,418)* y(k,44) +rxt(k,469) &
                 * y(k,46) +rxt(k,335)* y(k,47) +rxt(k,433)* y(k,48) +rxt(k,344) &
                 * y(k,51) +rxt(k,441)* y(k,52) +rxt(k,432)* y(k,53) +rxt(k,413) &
                 * y(k,54) +rxt(k,414)* y(k,56) +rxt(k,306)* y(k,92) &
                  + 2._r8*(rxt(k,415) +rxt(k,416))* y(k,104) + het_rates(k,104)) &
                 * y(k,104)
         prod(k,63) = (rxt(k,401)*y(k,33) +rxt(k,402)*y(k,35) +rxt(k,407)*y(k,43) + &
                 rxt(k,465)*y(k,40) +rxt(k,468)*y(k,46) +rxt(k,494)*y(k,36) + &
                 rxt(k,496)*y(k,45) +rxt(k,526)*y(k,14))*y(k,54) &
                  + (rxt(k,144)*y(k,99) +rxt(k,167)*y(k,77) +rxt(k,184)*y(k,73) + &
                 rxt(k,198)*y(k,75) +rxt(k,202)*y(k,92) +rxt(k,219)*y(k,89) + &
                 rxt(k,237)*y(k,88))*y(k,71) + (rxt(k,3) +rxt(k,134)*y(k,70) + &
                 rxt(k,323)*y(k,100) +rxt(k,350)*y(k,106) + &
                 2.000_r8*rxt(k,385)*y(k,95) +rxt(k,504)*y(k,31))*y(k,110) &
                  + (2.000_r8*rxt(k,404)*y(k,32) +rxt(k,408)*y(k,56) + &
                 rxt(k,429)*y(k,53) +rxt(k,434)*y(k,51) +rxt(k,450)*y(k,60))*y(k,43) &
                  + (rxt(k,98) +rxt(k,336)*y(k,60) +rxt(k,337)*y(k,24) + &
                 rxt(k,341)*y(k,52) +rxt(k,343)*y(k,56))*y(k,105) &
                  + (rxt(k,533)*y(k,20) +rxt(k,384)*y(k,33) +rxt(k,477)*y(k,40) + &
                 rxt(k,500)*y(k,36))*y(k,95) + (rxt(k,9) +rxt(k,119)*y(k,65) + &
                 rxt(k,352)*y(k,106))*y(k,42) + (rxt(k,23) + &
                 .300_r8*rxt(k,530)*y(k,104))*y(k,19) + (rxt(k,124)*y(k,68) + &
                 rxt(k,409)*y(k,56))*y(k,32) +2.000_r8*rxt(k,4)*y(k,35) &
                  +rxt(k,351)*y(k,106)*y(k,40) +rxt(k,10)*y(k,44) +rxt(k,58)*y(k,45) &
                  +rxt(k,59)*y(k,46) +rxt(k,12)*y(k,47) +.500_r8*rxt(k,539)*y(k,52) &
                  +rxt(k,133)*y(k,70)*y(k,69)
         loss(k,83) = (rxt(k,337)* y(k,24) +rxt(k,338)* y(k,28) +rxt(k,339)* y(k,32) &
                  +rxt(k,340)* y(k,40) +rxt(k,341)* y(k,52) +rxt(k,342)* y(k,54) &
                  +rxt(k,343)* y(k,56) +rxt(k,336)* y(k,60) + rxt(k,98) &
                  + het_rates(k,105))* y(k,105)
         prod(k,83) = (rxt(k,112)*y(k,65) +rxt(k,282)*y(k,81) +rxt(k,326)*y(k,101)) &
                 *y(k,32) + (rxt(k,349)*y(k,33) +rxt(k,350)*y(k,110))*y(k,106)
         loss(k,84) = (rxt(k,346)* y(k,24) +rxt(k,347)* y(k,28) + (rxt(k,348) + &
                 rxt(k,349))* y(k,33) +rxt(k,351)* y(k,40) +rxt(k,352)* y(k,42) &
                  +rxt(k,131)* y(k,51) +rxt(k,354)* y(k,52) +rxt(k,355)* y(k,54) &
                  +rxt(k,357)* y(k,55) +rxt(k,358)* y(k,56) +rxt(k,345)* y(k,60) &
                  +rxt(k,356)* y(k,96) +rxt(k,350)* y(k,110) + rxt(k,99) + rxt(k,353) &
                  + het_rates(k,106))* y(k,106)
         prod(k,84) =rxt(k,320)*y(k,98)*y(k,54) +rxt(k,129)*y(k,69)*y(k,56) +rxt(k,63) &
                 *y(k,65) +rxt(k,95)*y(k,101)
         loss(k,34) = (rxt(k,573)* y(k,28) +rxt(k,575)* y(k,55) + rxt(k,574) &
                  + het_rates(k,107))* y(k,107)
         prod(k,34) = (rxt(k,78) +rxt(k,79) +rxt(k,566)*y(k,80) +rxt(k,567)*y(k,94) + &
                 rxt(k,578)*y(k,108) +rxt(k,584)*y(k,109))*y(k,54) + (rxt(k,85) + &
                 rxt(k,87))*y(k,55) + (rxt(k,576)*y(k,108) +rxt(k,581)*y(k,109)) &
                 *y(k,69) +rxt(k,559)*y(k,108) +rxt(k,558)*y(k,109)
         loss(k,16) = (rxt(k,578)* y(k,54) +rxt(k,579)* y(k,55) +rxt(k,576)* y(k,69) &
                  + rxt(k,559) + rxt(k,577) + het_rates(k,108))* y(k,108)
         prod(k,16) = (rxt(k,80) +rxt(k,81))*y(k,54) + (rxt(k,90) +rxt(k,91))*y(k,55) &
                  + (rxt(k,560) +rxt(k,580)*y(k,69))*y(k,109)
         loss(k,15) = (rxt(k,584)* y(k,54) + (rxt(k,580) +rxt(k,581))* y(k,69) &
                  + rxt(k,558) + rxt(k,560) + rxt(k,582) + rxt(k,583) &
                  + het_rates(k,109))* y(k,109)
         prod(k,15) = (rxt(k,77) +rxt(k,82))*y(k,54) + (rxt(k,83) +rxt(k,92))*y(k,55)
         loss(k,88) = (rxt(k,504)* y(k,31) +rxt(k,102)* y(k,61) +rxt(k,114)* y(k,65) &
                  +rxt(k,113)* y(k,67) +rxt(k,134)* y(k,70) +rxt(k,369)* y(k,72) &
                  +rxt(k,371)* y(k,73) +rxt(k,135)* y(k,74) +rxt(k,373)* y(k,75) &
                  +rxt(k,138)* y(k,76) +rxt(k,367)* y(k,78) +rxt(k,283)* y(k,81) &
                  +rxt(k,291)* y(k,83) +rxt(k,289)* y(k,85) +rxt(k,360)* y(k,88) &
                  +rxt(k,362)* y(k,89) +rxt(k,300)* y(k,90) +rxt(k,302)* y(k,91) &
                  +rxt(k,361)* y(k,92) +rxt(k,308)* y(k,93) +rxt(k,385)* y(k,95) &
                  +rxt(k,324)* y(k,99) + (rxt(k,322) +rxt(k,323))* y(k,100) &
                  +rxt(k,332)* y(k,103) +rxt(k,350)* y(k,106) + rxt(k,1) + rxt(k,2) &
                  + rxt(k,3) + het_rates(k,110))* y(k,110)
         prod(k,88) = (rxt(k,372) +4.000_r8*rxt(k,136)*y(k,69) + &
                 4.000_r8*rxt(k,141)*y(k,87) +4.000_r8*rxt(k,142)*y(k,61) + &
                 5.000_r8*rxt(k,153)*y(k,67) +5.000_r8*rxt(k,157)*y(k,82) + &
                 4.000_r8*rxt(k,158)*y(k,86) +5.000_r8*rxt(k,159)*y(k,62) + &
                 6.000_r8*rxt(k,160)*y(k,84) +4.000_r8*rxt(k,161)*y(k,81) + &
                 4.000_r8*rxt(k,165)*y(k,65) +4.000_r8*rxt(k,176)*y(k,63) + &
                 4.000_r8*rxt(k,187)*y(k,83) +4.000_r8*rxt(k,198)*y(k,71) + &
                 4.000_r8*rxt(k,209)*y(k,98) +4.000_r8*rxt(k,220)*y(k,68) + &
                 5.000_r8*rxt(k,231)*y(k,85) +6.000_r8*rxt(k,242)*y(k,66) + &
                 4.000_r8*rxt(k,253)*y(k,65) +5.000_r8*rxt(k,255)*y(k,67) + &
                 5.000_r8*rxt(k,257)*y(k,85) +4.000_r8*rxt(k,259)*y(k,83) + &
                 4.000_r8*rxt(k,262)*y(k,63) +4.000_r8*rxt(k,264)*y(k,87) + &
                 6.000_r8*rxt(k,266)*y(k,66))*y(k,75) + (rxt(k,374) + &
                 5.000_r8*rxt(k,139)*y(k,69) +5.000_r8*rxt(k,162)*y(k,87) + &
                 5.000_r8*rxt(k,163)*y(k,65) +5.000_r8*rxt(k,164)*y(k,63) + &
                 5.000_r8*rxt(k,166)*y(k,83) +5.000_r8*rxt(k,167)*y(k,71) + &
                 5.000_r8*rxt(k,168)*y(k,98) +5.000_r8*rxt(k,169)*y(k,68) + &
                 6.000_r8*rxt(k,170)*y(k,85) +7.000_r8*rxt(k,171)*y(k,66) + &
                 5.000_r8*rxt(k,172)*y(k,61) +6.000_r8*rxt(k,173)*y(k,67) + &
                 6.000_r8*rxt(k,174)*y(k,82) +5.000_r8*rxt(k,175)*y(k,86) + &
                 6.000_r8*rxt(k,177)*y(k,62) +7.000_r8*rxt(k,178)*y(k,84) + &
                 5.000_r8*rxt(k,179)*y(k,81) +7.000_r8*rxt(k,254)*y(k,66) + &
                 6.000_r8*rxt(k,256)*y(k,67) +6.000_r8*rxt(k,258)*y(k,85) + &
                 5.000_r8*rxt(k,260)*y(k,65) +5.000_r8*rxt(k,261)*y(k,83) + &
                 5.000_r8*rxt(k,263)*y(k,63) +5.000_r8*rxt(k,265)*y(k,87))*y(k,77) &
                  + (rxt(k,370) +3.000_r8*rxt(k,180)*y(k,87) + &
                 3.000_r8*rxt(k,181)*y(k,65) +3.000_r8*rxt(k,182)*y(k,63) + &
                 3.000_r8*rxt(k,183)*y(k,83) +3.000_r8*rxt(k,184)*y(k,71) + &
                 3.000_r8*rxt(k,185)*y(k,98) +3.000_r8*rxt(k,186)*y(k,68) + &
                 4.000_r8*rxt(k,188)*y(k,85) +5.000_r8*rxt(k,189)*y(k,66) + &
                 3.000_r8*rxt(k,190)*y(k,61) +4.000_r8*rxt(k,191)*y(k,67) + &
                 4.000_r8*rxt(k,192)*y(k,82) +3.000_r8*rxt(k,193)*y(k,86) + &
                 4.000_r8*rxt(k,194)*y(k,62) +5.000_r8*rxt(k,195)*y(k,84) + &
                 3.000_r8*rxt(k,196)*y(k,81) +3.000_r8*rxt(k,378)*y(k,69))*y(k,73) &
                  + (rxt(k,509)*y(k,13) +rxt(k,511)*y(k,15) +rxt(k,512)*y(k,16) + &
                 rxt(k,514)*y(k,17) +rxt(k,519)*y(k,39) +rxt(k,531)*y(k,20) + &
                 rxt(k,335)*y(k,47) +rxt(k,410)*y(k,33) +rxt(k,411)*y(k,35) + &
                 rxt(k,412)*y(k,43) +rxt(k,415)*y(k,104) +rxt(k,418)*y(k,44) + &
                 rxt(k,442)*y(k,42) +rxt(k,466)*y(k,40) +rxt(k,469)*y(k,46) + &
                 rxt(k,495)*y(k,36) +rxt(k,527)*y(k,14) +rxt(k,530)*y(k,19))*y(k,104) &
                  + (2.000_r8*rxt(k,215)*y(k,87) +2.000_r8*rxt(k,216)*y(k,65) + &
                 2.000_r8*rxt(k,217)*y(k,63) +2.000_r8*rxt(k,218)*y(k,83) + &
                 2.000_r8*rxt(k,219)*y(k,71) +2.000_r8*rxt(k,221)*y(k,98) + &
                 2.000_r8*rxt(k,222)*y(k,68) +3.000_r8*rxt(k,223)*y(k,85) + &
                 4.000_r8*rxt(k,224)*y(k,66) +2.000_r8*rxt(k,225)*y(k,61) + &
                 3.000_r8*rxt(k,226)*y(k,67) +3.000_r8*rxt(k,227)*y(k,82) + &
                 2.000_r8*rxt(k,228)*y(k,86) +3.000_r8*rxt(k,229)*y(k,62) + &
                 4.000_r8*rxt(k,230)*y(k,84) +2.000_r8*rxt(k,232)*y(k,81) + &
                 2.000_r8*rxt(k,298)*y(k,69))*y(k,89) + (rxt(k,197)*y(k,87) + &
                 rxt(k,199)*y(k,65) +rxt(k,200)*y(k,63) +rxt(k,201)*y(k,83) + &
                 rxt(k,202)*y(k,71) +rxt(k,203)*y(k,98) +rxt(k,204)*y(k,68) + &
                 2.000_r8*rxt(k,205)*y(k,85) +3.000_r8*rxt(k,206)*y(k,66) + &
                 rxt(k,207)*y(k,61) +2.000_r8*rxt(k,208)*y(k,67) + &
                 2.000_r8*rxt(k,210)*y(k,82) +rxt(k,211)*y(k,86) + &
                 2.000_r8*rxt(k,212)*y(k,62) +3.000_r8*rxt(k,213)*y(k,84) + &
                 rxt(k,214)*y(k,81) +rxt(k,303)*y(k,69))*y(k,92) &
                  + (rxt(k,101)*y(k,62) +rxt(k,340)*y(k,105) +rxt(k,542)*y(k,46) + &
                 rxt(k,548)*y(k,46) +rxt(k,549)*y(k,45) +rxt(k,553)*y(k,46) + &
                 rxt(k,554)*y(k,45))*y(k,40) + (rxt(k,64) +rxt(k,278) + &
                 rxt(k,115)*y(k,52) +rxt(k,118)*y(k,51) +rxt(k,150)*y(k,99) + &
                 rxt(k,244)*y(k,88))*y(k,67) + (rxt(k,133)*y(k,70) + &
                 3.000_r8*rxt(k,299)*y(k,90) +rxt(k,321)*y(k,100) + &
                 rxt(k,376)*y(k,78) +2.000_r8*rxt(k,377)*y(k,72))*y(k,69) &
                  + (rxt(k,240)*y(k,85) +2.000_r8*rxt(k,241)*y(k,66) + &
                 rxt(k,245)*y(k,82) +rxt(k,247)*y(k,62) +2.000_r8*rxt(k,248)*y(k,84)) &
                 *y(k,88) + (rxt(k,147)*y(k,85) +2.000_r8*rxt(k,148)*y(k,66) + &
                 rxt(k,151)*y(k,82) +rxt(k,154)*y(k,62) +2.000_r8*rxt(k,155)*y(k,84)) &
                 *y(k,99) + (rxt(k,339)*y(k,105) +rxt(k,405)*y(k,43))*y(k,32) &
                  + (rxt(k,269) +rxt(k,288)*y(k,50))*y(k,84) + (rxt(k,277) + &
                 rxt(k,290)*y(k,42))*y(k,85) +rxt(k,348)*y(k,106)*y(k,33) +rxt(k,380) &
                 *y(k,62) +rxt(k,279)*y(k,66) +rxt(k,368)*y(k,72) +rxt(k,268)*y(k,82) &
                  +rxt(k,94)*y(k,100)
      end do
      end subroutine imp_prod_loss
      end module mo_prod_loss
