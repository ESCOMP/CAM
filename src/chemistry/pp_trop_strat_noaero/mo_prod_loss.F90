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
         loss(k,1) = ( + het_rates(k,3))* y(k,3)
         prod(k,1) = 0._r8
         loss(k,2) = ( + het_rates(k,19))* y(k,19)
         prod(k,2) = 0._r8
         loss(k,3) = (rxt(k,189)* y(k,189) + rxt(k,78) + het_rates(k,31))* y(k,31)
         prod(k,3) = 0._r8
         loss(k,4) = (rxt(k,190)* y(k,189) + rxt(k,79) + het_rates(k,32))* y(k,32)
         prod(k,4) = 0._r8
         loss(k,5) = (rxt(k,216)* y(k,189) + rxt(k,80) + het_rates(k,33))* y(k,33)
         prod(k,5) = 0._r8
         loss(k,6) = (rxt(k,191)* y(k,189) + rxt(k,81) + het_rates(k,34))* y(k,34)
         prod(k,6) = 0._r8
         loss(k,7) = (rxt(k,192)* y(k,189) + rxt(k,82) + het_rates(k,35))* y(k,35)
         prod(k,7) = 0._r8
         loss(k,8) = (rxt(k,193)* y(k,189) + rxt(k,83) + het_rates(k,36))* y(k,36)
         prod(k,8) = 0._r8
         loss(k,9) = (rxt(k,194)* y(k,189) + rxt(k,84) + het_rates(k,37))* y(k,37)
         prod(k,9) = 0._r8
         loss(k,10) = (rxt(k,195)* y(k,189) + rxt(k,85) + het_rates(k,38))* y(k,38)
         prod(k,10) = 0._r8
         loss(k,11) = (rxt(k,227)* y(k,54) +rxt(k,239)* y(k,189) +rxt(k,228)* y(k,190) &
                  + rxt(k,86) + het_rates(k,39))* y(k,39)
         prod(k,11) = 0._r8
         loss(k,12) = (rxt(k,229)* y(k,54) +rxt(k,240)* y(k,189) +rxt(k,230)* y(k,190) &
                  + rxt(k,87) + het_rates(k,41))* y(k,41)
         prod(k,12) = 0._r8
         loss(k,13) = (rxt(k,231)* y(k,190) + rxt(k,88) + het_rates(k,42))* y(k,42)
         prod(k,13) = 0._r8
         loss(k,14) = (rxt(k,232)* y(k,54) +rxt(k,233)* y(k,190) + rxt(k,89) &
                  + het_rates(k,44))* y(k,44)
         prod(k,14) = 0._r8
         loss(k,15) = (rxt(k,165)* y(k,54) +rxt(k,221)* y(k,68) + (rxt(k,261) + &
                 rxt(k,262) +rxt(k,263))* y(k,189) +rxt(k,254)* y(k,190) + rxt(k,39) &
                  + rxt(k,40) + het_rates(k,52))* y(k,52)
         prod(k,15) = 0._r8
         loss(k,16) = (rxt(k,234)* y(k,54) +rxt(k,217)* y(k,189) +rxt(k,235)* y(k,190) &
                  + rxt(k,90) + het_rates(k,53))* y(k,53)
         prod(k,16) = 0._r8
         loss(k,17) = ( + het_rates(k,59))* y(k,59)
         prod(k,17) = 0._r8
         loss(k,18) = ( + rxt(k,41) + het_rates(k,61))* y(k,61)
         prod(k,18) =.440_r8*rxt(k,40)*y(k,52)
         loss(k,19) = ( + rxt(k,529) + het_rates(k,66))* y(k,66)
         prod(k,19) = 0._r8
         loss(k,20) = (rxt(k,218)* y(k,189) + rxt(k,98) + het_rates(k,73))* y(k,73)
         prod(k,20) = 0._r8
         loss(k,21) = (rxt(k,241)* y(k,189) +rxt(k,236)* y(k,190) + rxt(k,100) &
                  + het_rates(k,77))* y(k,77)
         prod(k,21) = 0._r8
         loss(k,22) = (rxt(k,242)* y(k,189) +rxt(k,237)* y(k,190) + rxt(k,101) &
                  + het_rates(k,78))* y(k,78)
         prod(k,22) = 0._r8
         loss(k,23) = (rxt(k,243)* y(k,189) +rxt(k,238)* y(k,190) + rxt(k,102) &
                  + het_rates(k,79))* y(k,79)
         prod(k,23) = 0._r8
         loss(k,24) = ((rxt(k,156) +rxt(k,157))* y(k,189) + rxt(k,12) &
                  + het_rates(k,108))* y(k,108)
         prod(k,24) = 0._r8
         loss(k,25) = ( + rxt(k,531) + het_rates(k,114))* y(k,114)
         prod(k,25) = 0._r8
         loss(k,26) = ( + rxt(k,530) + het_rates(k,115))* y(k,115)
         prod(k,26) = 0._r8
         loss(k,27) = ( + het_rates(k,123))* y(k,123)
         prod(k,27) = 0._r8
         loss(k,28) = ( + rxt(k,108) + het_rates(k,135))* y(k,135)
         prod(k,28) = 0._r8
         loss(k,29) = ( + rxt(k,532) + het_rates(k,144))* y(k,144)
         prod(k,29) = 0._r8
         loss(k,30) = ( + het_rates(k,159))* y(k,159)
         prod(k,30) = 0._r8
         loss(k,31) = ( + het_rates(k,160))* y(k,160)
         prod(k,31) = 0._r8
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
         loss(k,101) = (rxt(k,346)* y(k,190) + rxt(k,19) + het_rates(k,1))* y(k,1)
         prod(k,101) =rxt(k,349)*y(k,162)*y(k,116)
         loss(k,100) = (rxt(k,350)* y(k,190) + rxt(k,20) + het_rates(k,2))* y(k,2)
         prod(k,100) =rxt(k,347)*y(k,176)*y(k,162)
         loss(k,124) = (rxt(k,429)* y(k,118) +rxt(k,430)* y(k,122) +rxt(k,431) &
                 * y(k,190) + het_rates(k,4))* y(k,4)
         prod(k,124) = 0._r8
         loss(k,26) = (rxt(k,388)* y(k,190) + het_rates(k,5))* y(k,5)
         prod(k,26) = 0._r8
         loss(k,70) = (rxt(k,391)* y(k,190) + rxt(k,21) + het_rates(k,6))* y(k,6)
         prod(k,70) =rxt(k,389)*y(k,176)*y(k,164)
         loss(k,27) = ( + rxt(k,22) + het_rates(k,7))* y(k,7)
         prod(k,27) =.120_r8*rxt(k,388)*y(k,190)*y(k,5)
         loss(k,71) = ( + rxt(k,23) + het_rates(k,8))* y(k,8)
         prod(k,71) = (.100_r8*rxt(k,430)*y(k,4) +.100_r8*rxt(k,433)*y(k,105)) &
                 *y(k,122)
         loss(k,76) = ( + rxt(k,24) + het_rates(k,9))* y(k,9)
         prod(k,76) = (.500_r8*rxt(k,390)*y(k,164) +.200_r8*rxt(k,417)*y(k,196) + &
                 .060_r8*rxt(k,423)*y(k,199))*y(k,116) +.500_r8*rxt(k,21)*y(k,6) &
                  +rxt(k,22)*y(k,7) +.200_r8*rxt(k,70)*y(k,152) +.060_r8*rxt(k,72) &
                 *y(k,156)
         loss(k,51) = ( + rxt(k,25) + het_rates(k,10))* y(k,10)
         prod(k,51) = (.200_r8*rxt(k,417)*y(k,196) +.200_r8*rxt(k,423)*y(k,199)) &
                 *y(k,116) +.200_r8*rxt(k,70)*y(k,152) +.200_r8*rxt(k,72)*y(k,156)
         loss(k,97) = ( + rxt(k,26) + het_rates(k,11))* y(k,11)
         prod(k,97) = (.200_r8*rxt(k,417)*y(k,196) +.150_r8*rxt(k,423)*y(k,199)) &
                 *y(k,116) +rxt(k,46)*y(k,89) +rxt(k,56)*y(k,111) +.200_r8*rxt(k,70) &
                 *y(k,152) +.150_r8*rxt(k,72)*y(k,156)
         loss(k,59) = ( + rxt(k,27) + het_rates(k,12))* y(k,12)
         prod(k,59) =.210_r8*rxt(k,423)*y(k,199)*y(k,116) +.210_r8*rxt(k,72)*y(k,156)
         loss(k,39) = (rxt(k,351)* y(k,190) + het_rates(k,13))* y(k,13)
         prod(k,39) = (.050_r8*rxt(k,430)*y(k,4) +.050_r8*rxt(k,433)*y(k,105)) &
                 *y(k,122)
         loss(k,61) = (rxt(k,317)* y(k,118) +rxt(k,318)* y(k,190) + het_rates(k,14)) &
                 * y(k,14)
         prod(k,61) = 0._r8
         loss(k,154) = (rxt(k,200)* y(k,40) +rxt(k,202)* y(k,122) +rxt(k,201) &
                 * y(k,176) + het_rates(k,15))* y(k,15)
         prod(k,154) = (rxt(k,75) +2.000_r8*rxt(k,203)*y(k,17) +rxt(k,204)*y(k,57) + &
                 rxt(k,205)*y(k,57) +rxt(k,208)*y(k,116) +rxt(k,211)*y(k,121) + &
                 rxt(k,212)*y(k,190) +rxt(k,458)*y(k,136))*y(k,17) &
                  + (rxt(k,190)*y(k,32) +rxt(k,216)*y(k,33) + &
                 3.000_r8*rxt(k,217)*y(k,53) +2.000_r8*rxt(k,218)*y(k,73) + &
                 2.000_r8*rxt(k,239)*y(k,39) +rxt(k,240)*y(k,41) +rxt(k,219)*y(k,76)) &
                 *y(k,189) + (2.000_r8*rxt(k,228)*y(k,39) +rxt(k,230)*y(k,41) + &
                 3.000_r8*rxt(k,235)*y(k,53) +rxt(k,214)*y(k,76))*y(k,190) &
                  + (2.000_r8*rxt(k,227)*y(k,39) +rxt(k,229)*y(k,41) + &
                 3.000_r8*rxt(k,234)*y(k,53))*y(k,54) + (rxt(k,99) + &
                 rxt(k,213)*y(k,121))*y(k,76) +rxt(k,74)*y(k,16) +rxt(k,77)*y(k,18) &
                  +rxt(k,105)*y(k,86)
         loss(k,42) = ( + rxt(k,74) + het_rates(k,16))* y(k,16)
         prod(k,42) = (rxt(k,520)*y(k,86) +rxt(k,525)*y(k,86))*y(k,80) &
                  +rxt(k,206)*y(k,57)*y(k,17)
         loss(k,160) = (2._r8*rxt(k,203)* y(k,17) + (rxt(k,204) +rxt(k,205) + &
                 rxt(k,206))* y(k,57) +rxt(k,208)* y(k,116) +rxt(k,209)* y(k,117) &
                  +rxt(k,211)* y(k,121) +rxt(k,458)* y(k,136) +rxt(k,207)* y(k,176) &
                  +rxt(k,212)* y(k,190) + rxt(k,75) + het_rates(k,17))* y(k,17)
         prod(k,160) = (rxt(k,76) +rxt(k,210)*y(k,121))*y(k,18) +rxt(k,202)*y(k,122) &
                 *y(k,15) +rxt(k,220)*y(k,189)*y(k,76) +rxt(k,215)*y(k,121)*y(k,86)
         loss(k,86) = (rxt(k,210)* y(k,121) + rxt(k,76) + rxt(k,77) + rxt(k,514) &
                  + rxt(k,517) + rxt(k,522) + het_rates(k,18))* y(k,18)
         prod(k,86) =rxt(k,209)*y(k,117)*y(k,17)
         loss(k,40) = (rxt(k,392)* y(k,190) + het_rates(k,20))* y(k,20)
         prod(k,40) =rxt(k,28)*y(k,21) +rxt(k,395)*y(k,166)*y(k,116)
         loss(k,58) = (rxt(k,394)* y(k,190) + rxt(k,28) + het_rates(k,21))* y(k,21)
         prod(k,58) =rxt(k,393)*y(k,176)*y(k,166)
         loss(k,50) = (rxt(k,266)* y(k,54) +rxt(k,267)* y(k,190) + het_rates(k,22)) &
                 * y(k,22)
         prod(k,50) = 0._r8
         loss(k,90) = (rxt(k,268)* y(k,54) +rxt(k,269)* y(k,122) +rxt(k,294)* y(k,190) &
                  + het_rates(k,23))* y(k,23)
         prod(k,90) = 0._r8
         loss(k,43) = (rxt(k,274)* y(k,190) + het_rates(k,24))* y(k,24)
         prod(k,43) = (.400_r8*rxt(k,270)*y(k,167) +.200_r8*rxt(k,271)*y(k,171)) &
                 *y(k,167)
         loss(k,53) = (rxt(k,275)* y(k,190) + rxt(k,29) + het_rates(k,25))* y(k,25)
         prod(k,53) =rxt(k,272)*y(k,176)*y(k,167)
         loss(k,49) = (rxt(k,276)* y(k,54) +rxt(k,277)* y(k,190) + het_rates(k,26)) &
                 * y(k,26)
         prod(k,49) = 0._r8
         loss(k,129) = (rxt(k,297)* y(k,118) +rxt(k,298)* y(k,122) +rxt(k,315) &
                 * y(k,190) + het_rates(k,27))* y(k,27)
         prod(k,129) =.130_r8*rxt(k,375)*y(k,122)*y(k,93) +.700_r8*rxt(k,55)*y(k,106)
         loss(k,69) = (rxt(k,302)* y(k,190) + rxt(k,30) + het_rates(k,28))* y(k,28)
         prod(k,69) =rxt(k,300)*y(k,176)*y(k,168)
         loss(k,19) = (rxt(k,303)* y(k,190) + het_rates(k,29))* y(k,29)
         prod(k,19) = 0._r8
         loss(k,44) = (rxt(k,398)* y(k,190) + rxt(k,31) + het_rates(k,30))* y(k,30)
         prod(k,44) =rxt(k,396)*y(k,176)*y(k,169)
         loss(k,161) = (rxt(k,200)* y(k,15) +rxt(k,164)* y(k,54) +rxt(k,245)* y(k,118) &
                  +rxt(k,246)* y(k,121) +rxt(k,244)* y(k,176) +rxt(k,247)* y(k,190) &
                  + rxt(k,32) + rxt(k,33) + het_rates(k,40))* y(k,40)
         prod(k,161) = (rxt(k,171)*y(k,57) +2.000_r8*rxt(k,248)*y(k,171) + &
                 rxt(k,249)*y(k,171) +rxt(k,251)*y(k,116) + &
                 .700_r8*rxt(k,271)*y(k,167) +rxt(k,282)*y(k,170) + &
                 rxt(k,299)*y(k,168) +.800_r8*rxt(k,311)*y(k,193) + &
                 .880_r8*rxt(k,323)*y(k,182) +2.000_r8*rxt(k,332)*y(k,184) + &
                 1.500_r8*rxt(k,356)*y(k,178) +.750_r8*rxt(k,361)*y(k,179) + &
                 .800_r8*rxt(k,370)*y(k,96) +.800_r8*rxt(k,381)*y(k,198) + &
                 .750_r8*rxt(k,435)*y(k,188) +.930_r8*rxt(k,440)*y(k,194) + &
                 .950_r8*rxt(k,445)*y(k,195))*y(k,171) &
                  + (.500_r8*rxt(k,288)*y(k,175) +rxt(k,309)*y(k,192) + &
                 rxt(k,313)*y(k,193) +.500_r8*rxt(k,319)*y(k,173) + &
                 .250_r8*rxt(k,326)*y(k,182) +rxt(k,335)*y(k,184) + &
                 .100_r8*rxt(k,348)*y(k,162) +.920_r8*rxt(k,358)*y(k,178) + &
                 .250_r8*rxt(k,383)*y(k,198) +.340_r8*rxt(k,442)*y(k,194) + &
                 .320_r8*rxt(k,447)*y(k,195))*y(k,116) + (rxt(k,252)*y(k,50) + &
                 .300_r8*rxt(k,253)*y(k,51) +.500_r8*rxt(k,286)*y(k,49) + &
                 .800_r8*rxt(k,291)*y(k,69) +rxt(k,293)*y(k,127) + &
                 .500_r8*rxt(k,341)*y(k,104) +.400_r8*rxt(k,346)*y(k,1) + &
                 .300_r8*rxt(k,366)*y(k,94) +.680_r8*rxt(k,451)*y(k,151))*y(k,190) &
                  + (rxt(k,269)*y(k,23) +.500_r8*rxt(k,298)*y(k,27) + &
                 .120_r8*rxt(k,328)*y(k,100) +.600_r8*rxt(k,342)*y(k,106) + &
                 .910_r8*rxt(k,375)*y(k,93) +.340_r8*rxt(k,430)*y(k,4) + &
                 .340_r8*rxt(k,433)*y(k,105))*y(k,122) + (.500_r8*rxt(k,317)*y(k,14) + &
                 .250_r8*rxt(k,325)*y(k,182) +rxt(k,336)*y(k,184) + &
                 rxt(k,359)*y(k,178))*y(k,118) + (.250_r8*rxt(k,322)*y(k,182) + &
                 rxt(k,331)*y(k,184) +rxt(k,355)*y(k,178) + &
                 .250_r8*rxt(k,380)*y(k,198))*y(k,170) + (rxt(k,262)*y(k,189) + &
                 rxt(k,263)*y(k,189))*y(k,52) + (.150_r8*rxt(k,312)*y(k,193) + &
                 .450_r8*rxt(k,333)*y(k,184))*y(k,176) +.100_r8*rxt(k,19)*y(k,1) &
                  +.100_r8*rxt(k,20)*y(k,2) +rxt(k,38)*y(k,51) +rxt(k,43)*y(k,69) &
                  +.330_r8*rxt(k,45)*y(k,88) +rxt(k,47)*y(k,90) +.690_r8*rxt(k,49) &
                 *y(k,98) +1.340_r8*rxt(k,50)*y(k,100) +rxt(k,57)*y(k,119) +rxt(k,62) &
                 *y(k,132) +rxt(k,63)*y(k,133) +.375_r8*rxt(k,65)*y(k,147) &
                  +.400_r8*rxt(k,67)*y(k,149) +.680_r8*rxt(k,69)*y(k,151) &
                  +2.000_r8*rxt(k,289)*y(k,174) +rxt(k,259)*y(k,177) &
                  +2.000_r8*rxt(k,334)*y(k,184)*y(k,184)
         loss(k,134) = (rxt(k,278)* y(k,118) +rxt(k,279)* y(k,190) + rxt(k,34) &
                  + het_rates(k,43))* y(k,43)
         prod(k,134) = (rxt(k,273)*y(k,167) +.270_r8*rxt(k,301)*y(k,168) + &
                 rxt(k,309)*y(k,192) +rxt(k,319)*y(k,173) +rxt(k,338)*y(k,186) + &
                 .400_r8*rxt(k,348)*y(k,162))*y(k,116) + (rxt(k,274)*y(k,24) + &
                 .500_r8*rxt(k,275)*y(k,25) +.800_r8*rxt(k,346)*y(k,1))*y(k,190) &
                  + (.500_r8*rxt(k,298)*y(k,27) +.100_r8*rxt(k,342)*y(k,106))*y(k,122) &
                  + (1.600_r8*rxt(k,270)*y(k,167) +.800_r8*rxt(k,271)*y(k,171)) &
                 *y(k,167) +.400_r8*rxt(k,19)*y(k,1) +.400_r8*rxt(k,20)*y(k,2) &
                  +rxt(k,317)*y(k,118)*y(k,14) +rxt(k,29)*y(k,25) +.330_r8*rxt(k,45) &
                 *y(k,88) +rxt(k,53)*y(k,103) +rxt(k,62)*y(k,132) &
                  +.200_r8*rxt(k,337)*y(k,186)*y(k,176)
         loss(k,16) = (rxt(k,280)* y(k,190) + het_rates(k,45))* y(k,45)
         prod(k,16) = 0._r8
         loss(k,127) = (rxt(k,316)* y(k,190) + rxt(k,35) + het_rates(k,46))* y(k,46)
         prod(k,127) = (.820_r8*rxt(k,301)*y(k,168) +.500_r8*rxt(k,319)*y(k,173) + &
                 .250_r8*rxt(k,348)*y(k,162) +.270_r8*rxt(k,442)*y(k,194) + &
                 .040_r8*rxt(k,447)*y(k,195))*y(k,116) &
                  + (.820_r8*rxt(k,299)*y(k,168) +.150_r8*rxt(k,440)*y(k,194) + &
                 .025_r8*rxt(k,445)*y(k,195))*y(k,171) + (.250_r8*rxt(k,19) + &
                 .800_r8*rxt(k,346)*y(k,190))*y(k,1) + (.520_r8*rxt(k,430)*y(k,4) + &
                 .520_r8*rxt(k,433)*y(k,105))*y(k,122) + (.500_r8*rxt(k,69) + &
                 .500_r8*rxt(k,451)*y(k,190))*y(k,151) +.250_r8*rxt(k,20)*y(k,2) &
                  +.500_r8*rxt(k,317)*y(k,118)*y(k,14) +.820_r8*rxt(k,30)*y(k,28) &
                  +.170_r8*rxt(k,45)*y(k,88) +.300_r8*rxt(k,65)*y(k,147) &
                  +.050_r8*rxt(k,67)*y(k,149)
         loss(k,146) = (rxt(k,304)* y(k,118) +rxt(k,305)* y(k,190) + rxt(k,36) &
                  + het_rates(k,47))* y(k,47)
         prod(k,146) = (.250_r8*rxt(k,326)*y(k,182) +.050_r8*rxt(k,364)*y(k,179) + &
                 .250_r8*rxt(k,383)*y(k,198) +.170_r8*rxt(k,401)*y(k,172) + &
                 .170_r8*rxt(k,407)*y(k,185) +.400_r8*rxt(k,417)*y(k,196) + &
                 .540_r8*rxt(k,423)*y(k,199) +.510_r8*rxt(k,426)*y(k,201))*y(k,116) &
                  + (.250_r8*rxt(k,325)*y(k,182) +.050_r8*rxt(k,365)*y(k,179) + &
                 .250_r8*rxt(k,384)*y(k,198))*y(k,118) &
                  + (.500_r8*rxt(k,311)*y(k,193) +.240_r8*rxt(k,323)*y(k,182) + &
                 .100_r8*rxt(k,381)*y(k,198))*y(k,171) &
                  + (.880_r8*rxt(k,328)*y(k,100) +.500_r8*rxt(k,342)*y(k,106)) &
                 *y(k,122) + (.250_r8*rxt(k,322)*y(k,182) + &
                 .250_r8*rxt(k,380)*y(k,198))*y(k,170) &
                  + (.070_r8*rxt(k,400)*y(k,172) +.070_r8*rxt(k,406)*y(k,185)) &
                 *y(k,176) + (rxt(k,306)*y(k,90) +rxt(k,307)*y(k,119))*y(k,190) &
                  +.180_r8*rxt(k,23)*y(k,8) +rxt(k,27)*y(k,12) +.400_r8*rxt(k,70) &
                 *y(k,152) +.540_r8*rxt(k,72)*y(k,156) +.510_r8*rxt(k,73)*y(k,158)
         loss(k,88) = (rxt(k,285)* y(k,190) + het_rates(k,48))* y(k,48)
         prod(k,88) = (.100_r8*rxt(k,282)*y(k,171) +.150_r8*rxt(k,283)*y(k,176)) &
                 *y(k,170) +.120_r8*rxt(k,298)*y(k,122)*y(k,27) &
                  +.150_r8*rxt(k,333)*y(k,184)*y(k,176)
         loss(k,84) = (rxt(k,286)* y(k,190) + rxt(k,37) + het_rates(k,49))* y(k,49)
         prod(k,84) = (.400_r8*rxt(k,283)*y(k,170) +.400_r8*rxt(k,333)*y(k,184)) &
                 *y(k,176)
         loss(k,113) = (rxt(k,252)* y(k,190) + het_rates(k,50))* y(k,50)
         prod(k,113) = (rxt(k,249)*y(k,171) +.300_r8*rxt(k,271)*y(k,167) + &
                 .500_r8*rxt(k,311)*y(k,193) +.250_r8*rxt(k,323)*y(k,182) + &
                 .250_r8*rxt(k,356)*y(k,178) +.250_r8*rxt(k,361)*y(k,179) + &
                 .200_r8*rxt(k,370)*y(k,96) +.300_r8*rxt(k,381)*y(k,198) + &
                 .250_r8*rxt(k,435)*y(k,188) +.250_r8*rxt(k,440)*y(k,194) + &
                 .250_r8*rxt(k,445)*y(k,195))*y(k,171)
         loss(k,67) = (rxt(k,253)* y(k,190) + rxt(k,38) + het_rates(k,51))* y(k,51)
         prod(k,67) =rxt(k,250)*y(k,176)*y(k,171)
         loss(k,169) = (rxt(k,276)* y(k,26) +rxt(k,227)* y(k,39) +rxt(k,164)* y(k,40) &
                  +rxt(k,229)* y(k,41) +rxt(k,232)* y(k,44) +rxt(k,165)* y(k,52) &
                  +rxt(k,234)* y(k,53) +rxt(k,177)* y(k,58) +rxt(k,166)* y(k,72) &
                  +rxt(k,167)* y(k,74) +rxt(k,186)* y(k,87) +rxt(k,170)* y(k,122) &
                  + (rxt(k,168) +rxt(k,169))* y(k,176) + het_rates(k,54))* y(k,54)
         prod(k,169) = (4.000_r8*rxt(k,189)*y(k,31) +rxt(k,190)*y(k,32) + &
                 2.000_r8*rxt(k,191)*y(k,34) +2.000_r8*rxt(k,192)*y(k,35) + &
                 2.000_r8*rxt(k,193)*y(k,36) +rxt(k,194)*y(k,37) + &
                 2.000_r8*rxt(k,195)*y(k,38) +rxt(k,241)*y(k,77) +rxt(k,242)*y(k,78) + &
                 rxt(k,243)*y(k,79) +rxt(k,196)*y(k,80) +rxt(k,226)*y(k,63))*y(k,189) &
                  + (rxt(k,93) +rxt(k,171)*y(k,171) +2.000_r8*rxt(k,172)*y(k,57) + &
                 rxt(k,174)*y(k,57) +rxt(k,176)*y(k,116) +rxt(k,181)*y(k,121) + &
                 rxt(k,182)*y(k,190) +rxt(k,205)*y(k,17) +rxt(k,459)*y(k,136))*y(k,57) &
                  + (3.000_r8*rxt(k,231)*y(k,42) +rxt(k,233)*y(k,44) + &
                 rxt(k,236)*y(k,77) +rxt(k,237)*y(k,78) +rxt(k,238)*y(k,79) + &
                 rxt(k,185)*y(k,80))*y(k,190) + (rxt(k,103) +rxt(k,184)*y(k,121)) &
                 *y(k,80) +rxt(k,74)*y(k,16) +2.000_r8*rxt(k,91)*y(k,55) &
                  +2.000_r8*rxt(k,92)*y(k,56) +rxt(k,94)*y(k,58) +rxt(k,97)*y(k,63) &
                  +rxt(k,106)*y(k,87)
         loss(k,25) = ( + rxt(k,91) + het_rates(k,55))* y(k,55)
         prod(k,25) = (rxt(k,513)*y(k,87) +rxt(k,518)*y(k,58) +rxt(k,519)*y(k,87) + &
                 rxt(k,523)*y(k,58) +rxt(k,524)*y(k,87) +rxt(k,528)*y(k,58))*y(k,80) &
                  +rxt(k,177)*y(k,58)*y(k,54) +rxt(k,173)*y(k,57)*y(k,57)
         loss(k,17) = ( + rxt(k,92) + rxt(k,199) + het_rates(k,56))* y(k,56)
         prod(k,17) =rxt(k,198)*y(k,57)*y(k,57)
         loss(k,164) = ((rxt(k,204) +rxt(k,205) +rxt(k,206))* y(k,17) &
                  + 2._r8*(rxt(k,172) +rxt(k,173) +rxt(k,174) +rxt(k,198))* y(k,57) &
                  +rxt(k,176)* y(k,116) +rxt(k,178)* y(k,117) +rxt(k,181)* y(k,121) &
                  +rxt(k,459)* y(k,136) +rxt(k,171)* y(k,171) +rxt(k,175)* y(k,176) &
                  + (rxt(k,182) +rxt(k,183))* y(k,190) + rxt(k,93) + het_rates(k,57)) &
                 * y(k,57)
         prod(k,164) = (rxt(k,169)*y(k,176) +rxt(k,170)*y(k,122) +rxt(k,186)*y(k,87)) &
                 *y(k,54) + (rxt(k,95) +rxt(k,179)*y(k,121))*y(k,58) &
                  + (rxt(k,187)*y(k,121) +rxt(k,188)*y(k,190))*y(k,87) + (rxt(k,107) + &
                 rxt(k,464)*y(k,136))*y(k,124) +2.000_r8*rxt(k,199)*y(k,56) &
                  +rxt(k,197)*y(k,189)*y(k,80)
         loss(k,128) = (rxt(k,177)* y(k,54) + (rxt(k,518) +rxt(k,523) +rxt(k,528)) &
                 * y(k,80) +rxt(k,179)* y(k,121) +rxt(k,180)* y(k,190) + rxt(k,94) &
                  + rxt(k,95) + rxt(k,516) + rxt(k,521) + rxt(k,527) &
                  + het_rates(k,58))* y(k,58)
         prod(k,128) =rxt(k,178)*y(k,117)*y(k,57)
         loss(k,135) = ((rxt(k,255) +rxt(k,265))* y(k,190) + het_rates(k,60))* y(k,60)
         prod(k,135) = (rxt(k,32) +rxt(k,33) +rxt(k,164)*y(k,54) +rxt(k,200)*y(k,15) + &
                 rxt(k,245)*y(k,118) +rxt(k,246)*y(k,121) +rxt(k,247)*y(k,190)) &
                 *y(k,40) + (.630_r8*rxt(k,269)*y(k,23) +.560_r8*rxt(k,298)*y(k,27) + &
                 .650_r8*rxt(k,328)*y(k,100) +.560_r8*rxt(k,342)*y(k,106) + &
                 .620_r8*rxt(k,375)*y(k,93) +.230_r8*rxt(k,430)*y(k,4) + &
                 .230_r8*rxt(k,433)*y(k,105))*y(k,122) &
                  + (.220_r8*rxt(k,326)*y(k,182) +.250_r8*rxt(k,383)*y(k,198) + &
                 .170_r8*rxt(k,401)*y(k,172) +.400_r8*rxt(k,404)*y(k,183) + &
                 .350_r8*rxt(k,407)*y(k,185) +.225_r8*rxt(k,442)*y(k,194))*y(k,116) &
                  + (.350_r8*rxt(k,267)*y(k,22) +rxt(k,292)*y(k,70) + &
                 rxt(k,305)*y(k,47) +.700_r8*rxt(k,451)*y(k,151) +rxt(k,455)*y(k,125)) &
                 *y(k,190) + (rxt(k,304)*y(k,47) +.220_r8*rxt(k,325)*y(k,182) + &
                 .500_r8*rxt(k,384)*y(k,198))*y(k,118) &
                  + (.110_r8*rxt(k,323)*y(k,182) +.200_r8*rxt(k,381)*y(k,198) + &
                 .125_r8*rxt(k,440)*y(k,194))*y(k,171) &
                  + (.070_r8*rxt(k,400)*y(k,172) +.160_r8*rxt(k,403)*y(k,183) + &
                 .140_r8*rxt(k,406)*y(k,185))*y(k,176) + (rxt(k,110) + &
                 rxt(k,454)*y(k,121))*y(k,125) + (.220_r8*rxt(k,322)*y(k,182) + &
                 .250_r8*rxt(k,380)*y(k,198))*y(k,170) +1.500_r8*rxt(k,22)*y(k,7) &
                  +.450_r8*rxt(k,23)*y(k,8) +.600_r8*rxt(k,26)*y(k,11) +rxt(k,27) &
                 *y(k,12) +rxt(k,34)*y(k,43) +rxt(k,232)*y(k,54)*y(k,44) +rxt(k,36) &
                 *y(k,47) +rxt(k,43)*y(k,69) +2.000_r8*rxt(k,44)*y(k,70) &
                  +.330_r8*rxt(k,45)*y(k,88) +1.340_r8*rxt(k,51)*y(k,100) &
                  +.700_r8*rxt(k,55)*y(k,106) +1.500_r8*rxt(k,64)*y(k,146) &
                  +.250_r8*rxt(k,65)*y(k,147) +rxt(k,68)*y(k,150) +1.700_r8*rxt(k,69) &
                 *y(k,151)
         loss(k,20) = (rxt(k,225)* y(k,189) + rxt(k,96) + het_rates(k,62))* y(k,62)
         prod(k,20) = (rxt(k,190)*y(k,32) +rxt(k,192)*y(k,35) + &
                 2.000_r8*rxt(k,193)*y(k,36) +2.000_r8*rxt(k,194)*y(k,37) + &
                 rxt(k,195)*y(k,38) +rxt(k,216)*y(k,33) +2.000_r8*rxt(k,218)*y(k,73) + &
                 rxt(k,242)*y(k,78) +rxt(k,243)*y(k,79))*y(k,189) &
                  + (rxt(k,237)*y(k,78) +rxt(k,238)*y(k,79))*y(k,190)
         loss(k,29) = (rxt(k,226)* y(k,189) + rxt(k,97) + het_rates(k,63))* y(k,63)
         prod(k,29) = (rxt(k,191)*y(k,34) +rxt(k,192)*y(k,35) +rxt(k,241)*y(k,77)) &
                 *y(k,189) +rxt(k,236)*y(k,190)*y(k,77)
         loss(k,32) = (rxt(k,399)* y(k,190) + het_rates(k,64))* y(k,64)
         prod(k,32) =.180_r8*rxt(k,419)*y(k,190)*y(k,153)
         loss(k,47) = (rxt(k,452)* y(k,118) + (rxt(k,453) +rxt(k,466))* y(k,190) &
                  + het_rates(k,65))* y(k,65)
         prod(k,47) = 0._r8
         loss(k,21) = ( + rxt(k,42) + het_rates(k,67))* y(k,67)
         prod(k,21) =rxt(k,287)*y(k,176)*y(k,175)
         loss(k,111) = (rxt(k,221)* y(k,52) +rxt(k,222)* y(k,72) +rxt(k,224)* y(k,84) &
                  +rxt(k,223)* y(k,202) + het_rates(k,68))* y(k,68)
         prod(k,111) = (rxt(k,194)*y(k,37) +rxt(k,216)*y(k,33) + &
                 2.000_r8*rxt(k,225)*y(k,62) +rxt(k,226)*y(k,63))*y(k,189) &
                  +2.000_r8*rxt(k,96)*y(k,62) +rxt(k,97)*y(k,63) +rxt(k,104)*y(k,83)
         loss(k,131) = (rxt(k,291)* y(k,190) + rxt(k,43) + het_rates(k,69))* y(k,69)
         prod(k,131) = (.530_r8*rxt(k,326)*y(k,182) +.050_r8*rxt(k,364)*y(k,179) + &
                 .250_r8*rxt(k,383)*y(k,198) +.225_r8*rxt(k,442)*y(k,194))*y(k,116) &
                  + (.530_r8*rxt(k,325)*y(k,182) +.050_r8*rxt(k,365)*y(k,179) + &
                 .250_r8*rxt(k,384)*y(k,198))*y(k,118) &
                  + (.260_r8*rxt(k,323)*y(k,182) +.100_r8*rxt(k,381)*y(k,198) + &
                 .125_r8*rxt(k,440)*y(k,194))*y(k,171) + (.700_r8*rxt(k,366)*y(k,94) + &
                 .500_r8*rxt(k,367)*y(k,95) +rxt(k,378)*y(k,110))*y(k,190) &
                  + (.530_r8*rxt(k,322)*y(k,182) +.250_r8*rxt(k,380)*y(k,198)) &
                 *y(k,170) +.330_r8*rxt(k,45)*y(k,88) +.250_r8*rxt(k,65)*y(k,147) &
                  +rxt(k,290)*y(k,174)
         loss(k,122) = (rxt(k,292)* y(k,190) + rxt(k,44) + rxt(k,470) &
                  + het_rates(k,70))* y(k,70)
         prod(k,122) = (.050_r8*rxt(k,364)*y(k,179) +.250_r8*rxt(k,383)*y(k,198) + &
                 rxt(k,390)*y(k,164) +.400_r8*rxt(k,404)*y(k,183) + &
                 .170_r8*rxt(k,407)*y(k,185) +.700_r8*rxt(k,410)*y(k,191) + &
                 .600_r8*rxt(k,417)*y(k,196) +.340_r8*rxt(k,423)*y(k,199) + &
                 .170_r8*rxt(k,426)*y(k,201))*y(k,116) + (.650_r8*rxt(k,267)*y(k,22) + &
                 .200_r8*rxt(k,291)*y(k,69) +rxt(k,379)*y(k,111))*y(k,190) &
                  + (.250_r8*rxt(k,380)*y(k,170) +.100_r8*rxt(k,381)*y(k,171) + &
                 .250_r8*rxt(k,384)*y(k,118))*y(k,198) &
                  + (.160_r8*rxt(k,403)*y(k,183) +.070_r8*rxt(k,406)*y(k,185)) &
                 *y(k,176) +rxt(k,21)*y(k,6) +.130_r8*rxt(k,23)*y(k,8) &
                  +.050_r8*rxt(k,365)*y(k,179)*y(k,118) +.700_r8*rxt(k,61)*y(k,131) &
                  +.600_r8*rxt(k,70)*y(k,152) +.340_r8*rxt(k,72)*y(k,156) &
                  +.170_r8*rxt(k,73)*y(k,158)
         loss(k,156) = (rxt(k,130)* y(k,122) + (rxt(k,124) +rxt(k,125) +rxt(k,126)) &
                 * y(k,176) + rxt(k,127) + het_rates(k,71))* y(k,71)
         prod(k,156) = (rxt(k,131)*y(k,72) +rxt(k,134)*y(k,121) +rxt(k,152)*y(k,107) + &
                 rxt(k,247)*y(k,40) +rxt(k,265)*y(k,60) +rxt(k,455)*y(k,125) + &
                 rxt(k,460)*y(k,134) +rxt(k,465)*y(k,136))*y(k,190) &
                  + (rxt(k,114)*y(k,189) +rxt(k,122)*y(k,121) +rxt(k,166)*y(k,54) + &
                 rxt(k,222)*y(k,68))*y(k,72) + (rxt(k,262)*y(k,52) + &
                 rxt(k,197)*y(k,80) +rxt(k,220)*y(k,76))*y(k,189) + (rxt(k,2) + &
                 2.000_r8*rxt(k,3))*y(k,202) +2.000_r8*rxt(k,32)*y(k,40) +rxt(k,38) &
                 *y(k,51) +rxt(k,99)*y(k,76) +rxt(k,103)*y(k,80) +rxt(k,104)*y(k,83)
         loss(k,143) = (rxt(k,166)* y(k,54) +rxt(k,222)* y(k,68) +rxt(k,122)* y(k,121) &
                  +rxt(k,114)* y(k,189) +rxt(k,131)* y(k,190) + het_rates(k,72)) &
                 * y(k,72)
         prod(k,143) =rxt(k,33)*y(k,40) +rxt(k,263)*y(k,189)*y(k,52) &
                  +rxt(k,124)*y(k,176)*y(k,71) +rxt(k,1)*y(k,202)
         loss(k,93) = (rxt(k,167)* y(k,54) +rxt(k,123)* y(k,121) +rxt(k,132)* y(k,190) &
                  + rxt(k,4) + het_rates(k,74))* y(k,74)
         prod(k,93) =rxt(k,138)*y(k,176)*y(k,176) +rxt(k,137)*y(k,190)*y(k,190)
         loss(k,22) = ( + rxt(k,109) + het_rates(k,75))* y(k,75)
         prod(k,22) =rxt(k,468)*y(k,202)*y(k,138)
         loss(k,118) = (rxt(k,213)* y(k,121) + (rxt(k,219) +rxt(k,220))* y(k,189) &
                  +rxt(k,214)* y(k,190) + rxt(k,99) + het_rates(k,76))* y(k,76)
         prod(k,118) = (rxt(k,200)*y(k,40) +rxt(k,201)*y(k,176))*y(k,15)
         loss(k,168) = ((rxt(k,518) +rxt(k,523) +rxt(k,528))* y(k,58) + (rxt(k,520) + &
                 rxt(k,525))* y(k,86) + (rxt(k,513) +rxt(k,519) +rxt(k,524))* y(k,87) &
                  +rxt(k,184)* y(k,121) + (rxt(k,196) +rxt(k,197))* y(k,189) &
                  +rxt(k,185)* y(k,190) + rxt(k,103) + het_rates(k,80))* y(k,80)
         prod(k,168) = (rxt(k,165)*y(k,52) +rxt(k,227)*y(k,39) +rxt(k,229)*y(k,41) + &
                 2.000_r8*rxt(k,232)*y(k,44) +rxt(k,234)*y(k,53) +rxt(k,164)*y(k,40) + &
                 rxt(k,166)*y(k,72) +rxt(k,167)*y(k,74) +rxt(k,168)*y(k,176) + &
                 rxt(k,186)*y(k,87) +rxt(k,276)*y(k,26))*y(k,54) +rxt(k,183)*y(k,190) &
                 *y(k,57)
         loss(k,30) = (rxt(k,264)* y(k,189) +rxt(k,256)* y(k,190) + het_rates(k,81)) &
                 * y(k,81)
         prod(k,30) = 0._r8
         loss(k,115) = (rxt(k,257)* y(k,190) + het_rates(k,82))* y(k,82)
         prod(k,115) = (.370_r8*rxt(k,269)*y(k,23) +.120_r8*rxt(k,298)*y(k,27) + &
                 .330_r8*rxt(k,328)*y(k,100) +.120_r8*rxt(k,342)*y(k,106) + &
                 .110_r8*rxt(k,375)*y(k,93) +.050_r8*rxt(k,430)*y(k,4) + &
                 .050_r8*rxt(k,433)*y(k,105))*y(k,122) + (rxt(k,258)*y(k,176) + &
                 rxt(k,260)*y(k,116))*y(k,177) +.350_r8*rxt(k,267)*y(k,190)*y(k,22)
         loss(k,37) = ( + rxt(k,104) + het_rates(k,83))* y(k,83)
         prod(k,37) = (rxt(k,221)*y(k,52) +rxt(k,222)*y(k,72) +rxt(k,223)*y(k,202) + &
                 rxt(k,224)*y(k,84))*y(k,68)
         loss(k,155) = (rxt(k,224)* y(k,68) +rxt(k,161)* y(k,190) + rxt(k,9) &
                  + het_rates(k,84))* y(k,84)
         prod(k,155) = (rxt(k,516) +rxt(k,521) +rxt(k,527) +rxt(k,518)*y(k,80) + &
                 rxt(k,523)*y(k,80) +rxt(k,528)*y(k,80))*y(k,58) + (rxt(k,480) + &
                 rxt(k,245)*y(k,40) +rxt(k,278)*y(k,43) +rxt(k,304)*y(k,47) + &
                 rxt(k,452)*y(k,65))*y(k,118) + (2.000_r8*rxt(k,475) + &
                 2.000_r8*rxt(k,512) +2.000_r8*rxt(k,515) +2.000_r8*rxt(k,526)) &
                 *y(k,109) + (rxt(k,514) +rxt(k,517) +rxt(k,522))*y(k,18) &
                  + (.500_r8*rxt(k,479) +rxt(k,160)*y(k,190))*y(k,117) +rxt(k,472) &
                 *y(k,88) +rxt(k,473)*y(k,94) +rxt(k,474)*y(k,95) +rxt(k,476)*y(k,110) &
                  +rxt(k,477)*y(k,111) +rxt(k,481)*y(k,120) +rxt(k,482)*y(k,126) &
                  +rxt(k,483)*y(k,148)
         loss(k,62) = (rxt(k,139)* y(k,190) + rxt(k,10) + rxt(k,11) + rxt(k,162) &
                  + het_rates(k,85))* y(k,85)
         prod(k,62) =rxt(k,158)*y(k,176)*y(k,117)
         loss(k,112) = ((rxt(k,520) +rxt(k,525))* y(k,80) +rxt(k,215)* y(k,121) &
                  + rxt(k,105) + het_rates(k,86))* y(k,86)
         prod(k,112) = (rxt(k,514) +rxt(k,517) +rxt(k,522))*y(k,18) &
                  +rxt(k,207)*y(k,176)*y(k,17)
         loss(k,116) = (rxt(k,186)* y(k,54) + (rxt(k,513) +rxt(k,519) +rxt(k,524)) &
                 * y(k,80) +rxt(k,187)* y(k,121) +rxt(k,188)* y(k,190) + rxt(k,106) &
                  + het_rates(k,87))* y(k,87)
         prod(k,116) = (rxt(k,516) +rxt(k,521) +rxt(k,527) +rxt(k,180)*y(k,190)) &
                 *y(k,58) +rxt(k,175)*y(k,176)*y(k,57)
         loss(k,136) = (rxt(k,321)* y(k,190) + rxt(k,45) + rxt(k,472) &
                  + het_rates(k,88))* y(k,88)
         prod(k,136) = (rxt(k,320)*y(k,173) +rxt(k,327)*y(k,182))*y(k,116) &
                  + (.300_r8*rxt(k,366)*y(k,94) +.500_r8*rxt(k,367)*y(k,95))*y(k,190)
         loss(k,38) = (rxt(k,352)* y(k,190) + rxt(k,46) + het_rates(k,89))* y(k,89)
         prod(k,38) =rxt(k,363)*y(k,179)
         loss(k,137) = (rxt(k,306)* y(k,190) + rxt(k,47) + het_rates(k,90))* y(k,90)
         prod(k,137) = (.220_r8*rxt(k,322)*y(k,170) +.230_r8*rxt(k,323)*y(k,171) + &
                 .220_r8*rxt(k,325)*y(k,118) +.220_r8*rxt(k,326)*y(k,116))*y(k,182) &
                  + (.500_r8*rxt(k,310)*y(k,132) +.500_r8*rxt(k,341)*y(k,104) + &
                 .700_r8*rxt(k,366)*y(k,94) +.500_r8*rxt(k,367)*y(k,95))*y(k,190) &
                  + (.250_r8*rxt(k,380)*y(k,170) +.100_r8*rxt(k,381)*y(k,171) + &
                 .250_r8*rxt(k,383)*y(k,116) +.250_r8*rxt(k,384)*y(k,118))*y(k,198) &
                  + (.050_r8*rxt(k,364)*y(k,116) +.050_r8*rxt(k,365)*y(k,118)) &
                 *y(k,179) +.170_r8*rxt(k,45)*y(k,88) +.200_r8*rxt(k,311)*y(k,193) &
                 *y(k,171)
         loss(k,52) = (rxt(k,353)* y(k,190) + het_rates(k,91))* y(k,91)
         prod(k,52) = (rxt(k,360)*y(k,170) +.750_r8*rxt(k,361)*y(k,171) + &
                 .870_r8*rxt(k,364)*y(k,116) +.950_r8*rxt(k,365)*y(k,118))*y(k,179)
         loss(k,23) = (rxt(k,354)* y(k,190) + het_rates(k,92))* y(k,92)
         prod(k,23) =.600_r8*rxt(k,377)*y(k,190)*y(k,98)
         loss(k,119) = (rxt(k,368)* y(k,118) +rxt(k,375)* y(k,122) +rxt(k,376) &
                 * y(k,190) + het_rates(k,93))* y(k,93)
         prod(k,119) = 0._r8
         loss(k,94) = (rxt(k,366)* y(k,190) + rxt(k,473) + het_rates(k,94))* y(k,94)
         prod(k,94) =.080_r8*rxt(k,358)*y(k,178)*y(k,116)
         loss(k,89) = (rxt(k,367)* y(k,190) + rxt(k,474) + het_rates(k,95))* y(k,95)
         prod(k,89) =.080_r8*rxt(k,364)*y(k,179)*y(k,116)
         loss(k,144) = (rxt(k,372)* y(k,116) +rxt(k,373)* y(k,118) +rxt(k,369) &
                 * y(k,170) +rxt(k,370)* y(k,171) +rxt(k,371)* y(k,176) &
                  + het_rates(k,96))* y(k,96)
         prod(k,144) =rxt(k,368)*y(k,118)*y(k,93)
         loss(k,63) = (rxt(k,374)* y(k,190) + rxt(k,48) + het_rates(k,97))* y(k,97)
         prod(k,63) =rxt(k,371)*y(k,176)*y(k,96)
         loss(k,107) = (rxt(k,377)* y(k,190) + rxt(k,49) + het_rates(k,98))* y(k,98)
         prod(k,107) = (rxt(k,357)*y(k,178) +rxt(k,362)*y(k,179))*y(k,176) +rxt(k,48) &
                 *y(k,97)
         loss(k,12) = (rxt(k,499)* y(k,190) + het_rates(k,99))* y(k,99)
         prod(k,12) = 0._r8
         loss(k,145) = (rxt(k,328)* y(k,122) +rxt(k,329)* y(k,190) + rxt(k,50) &
                  + rxt(k,51) + het_rates(k,100))* y(k,100)
         prod(k,145) = (.390_r8*rxt(k,355)*y(k,170) +.310_r8*rxt(k,356)*y(k,171) + &
                 .360_r8*rxt(k,358)*y(k,116) +.400_r8*rxt(k,359)*y(k,118))*y(k,178) &
                  +.300_r8*rxt(k,375)*y(k,122)*y(k,93) +.288_r8*rxt(k,49)*y(k,98)
         loss(k,54) = (rxt(k,330)* y(k,190) + het_rates(k,101))* y(k,101)
         prod(k,54) =rxt(k,324)*y(k,182)*y(k,176)
         loss(k,83) = (rxt(k,339)* y(k,190) + rxt(k,52) + het_rates(k,102))* y(k,102)
         prod(k,83) =.800_r8*rxt(k,19)*y(k,1) +.800_r8*rxt(k,20)*y(k,2) &
                  +.800_r8*rxt(k,348)*y(k,162)*y(k,116)
         loss(k,55) = (rxt(k,340)* y(k,190) + rxt(k,53) + het_rates(k,103))* y(k,103)
         prod(k,55) =.800_r8*rxt(k,337)*y(k,186)*y(k,176)
         loss(k,92) = (rxt(k,341)* y(k,190) + rxt(k,54) + rxt(k,345) &
                  + het_rates(k,104))* y(k,104)
         prod(k,92) =rxt(k,344)*y(k,184)*y(k,117)
         loss(k,125) = (rxt(k,432)* y(k,118) +rxt(k,433)* y(k,122) +rxt(k,434) &
                 * y(k,190) + het_rates(k,105))* y(k,105)
         prod(k,125) = 0._r8
         loss(k,149) = (rxt(k,342)* y(k,122) +rxt(k,343)* y(k,190) + rxt(k,55) &
                  + het_rates(k,106))* y(k,106)
         prod(k,149) = (.610_r8*rxt(k,355)*y(k,170) +.440_r8*rxt(k,356)*y(k,171) + &
                 .560_r8*rxt(k,358)*y(k,116) +.600_r8*rxt(k,359)*y(k,118))*y(k,178) &
                  +.200_r8*rxt(k,375)*y(k,122)*y(k,93) +.402_r8*rxt(k,49)*y(k,98)
         loss(k,64) = (rxt(k,140)* y(k,116) + (rxt(k,141) +rxt(k,142) +rxt(k,143)) &
                 * y(k,117) +rxt(k,152)* y(k,190) + rxt(k,144) + het_rates(k,107)) &
                 * y(k,107)
         prod(k,64) =rxt(k,15)*y(k,116)
         loss(k,48) = ( + rxt(k,13) + rxt(k,14) + rxt(k,163) + rxt(k,475) + rxt(k,512) &
                  + rxt(k,515) + rxt(k,526) + het_rates(k,109))* y(k,109)
         prod(k,48) =rxt(k,159)*y(k,118)*y(k,117)
         loss(k,65) = (rxt(k,378)* y(k,190) + rxt(k,476) + het_rates(k,110))* y(k,110)
         prod(k,65) =.200_r8*rxt(k,370)*y(k,171)*y(k,96)
         loss(k,132) = (rxt(k,379)* y(k,190) + rxt(k,56) + rxt(k,477) &
                  + het_rates(k,111))* y(k,111)
         prod(k,132) = (rxt(k,369)*y(k,170) +.800_r8*rxt(k,370)*y(k,171) + &
                 rxt(k,372)*y(k,116) +rxt(k,373)*y(k,118))*y(k,96)
         loss(k,18) = (rxt(k,469)* y(k,190) + het_rates(k,112))* y(k,112)
         prod(k,18) = 0._r8
         loss(k,1) = ( + rxt(k,478) + het_rates(k,113))* y(k,113)
         prod(k,1) = 0._r8
         loss(k,165) = (rxt(k,208)* y(k,17) +rxt(k,176)* y(k,57) +rxt(k,372)* y(k,96) &
                  +rxt(k,140)* y(k,107) +rxt(k,149)* y(k,118) +rxt(k,155)* y(k,121) &
                  +rxt(k,154)* y(k,122) +rxt(k,387)* y(k,161) + (rxt(k,348) + &
                 rxt(k,349))* y(k,162) +rxt(k,390)* y(k,164) +rxt(k,395)* y(k,166) &
                  +rxt(k,273)* y(k,167) +rxt(k,301)* y(k,168) +rxt(k,397)* y(k,169) &
                  +rxt(k,284)* y(k,170) +rxt(k,251)* y(k,171) +rxt(k,401)* y(k,172) &
                  + (rxt(k,319) +rxt(k,320))* y(k,173) +rxt(k,288)* y(k,175) &
                  +rxt(k,153)* y(k,176) +rxt(k,260)* y(k,177) +rxt(k,358)* y(k,178) &
                  +rxt(k,364)* y(k,179) + (rxt(k,326) +rxt(k,327))* y(k,182) &
                  +rxt(k,404)* y(k,183) +rxt(k,335)* y(k,184) +rxt(k,407)* y(k,185) &
                  +rxt(k,338)* y(k,186) +rxt(k,437)* y(k,188) +rxt(k,410)* y(k,191) &
                  +rxt(k,309)* y(k,192) +rxt(k,313)* y(k,193) +rxt(k,442)* y(k,194) &
                  +rxt(k,447)* y(k,195) +rxt(k,417)* y(k,196) +rxt(k,383)* y(k,198) &
                  +rxt(k,423)* y(k,199) +rxt(k,426)* y(k,201) + rxt(k,15) &
                  + het_rates(k,116))* y(k,116)
         prod(k,165) = (rxt(k,16) +.500_r8*rxt(k,479) +2.000_r8*rxt(k,142)*y(k,107) + &
                 rxt(k,145)*y(k,121) +rxt(k,461)*y(k,136))*y(k,117) + (rxt(k,144) + &
                 rxt(k,152)*y(k,190))*y(k,107) +2.000_r8*rxt(k,156)*y(k,189)*y(k,108) &
                  +rxt(k,14)*y(k,109) +rxt(k,17)*y(k,118)
         loss(k,167) = (rxt(k,209)* y(k,17) +rxt(k,178)* y(k,57) + (rxt(k,141) + &
                 rxt(k,142) +rxt(k,143))* y(k,107) +rxt(k,159)* y(k,118) &
                  + (rxt(k,145) +rxt(k,147))* y(k,121) +rxt(k,146)* y(k,122) &
                  +rxt(k,412)* y(k,129) +rxt(k,461)* y(k,136) +rxt(k,415)* y(k,161) &
                  +rxt(k,295)* y(k,170) +rxt(k,402)* y(k,172) +rxt(k,158)* y(k,176) &
                  +rxt(k,405)* y(k,183) +rxt(k,344)* y(k,184) +rxt(k,408)* y(k,185) &
                  +rxt(k,160)* y(k,190) + rxt(k,16) + rxt(k,479) + het_rates(k,117)) &
                 * y(k,117)
         prod(k,167) = (2.000_r8*rxt(k,149)*y(k,118) +rxt(k,153)*y(k,176) + &
                 rxt(k,154)*y(k,122) +rxt(k,155)*y(k,121) +rxt(k,176)*y(k,57) + &
                 rxt(k,208)*y(k,17) +rxt(k,251)*y(k,171) +rxt(k,260)*y(k,177) + &
                 rxt(k,273)*y(k,167) +rxt(k,284)*y(k,170) +rxt(k,288)*y(k,175) + &
                 rxt(k,301)*y(k,168) +rxt(k,309)*y(k,192) +rxt(k,313)*y(k,193) + &
                 rxt(k,319)*y(k,173) +rxt(k,326)*y(k,182) +rxt(k,335)*y(k,184) + &
                 rxt(k,338)*y(k,186) +rxt(k,348)*y(k,162) + &
                 .920_r8*rxt(k,358)*y(k,178) +.920_r8*rxt(k,364)*y(k,179) + &
                 rxt(k,372)*y(k,96) +rxt(k,383)*y(k,198) +rxt(k,387)*y(k,161) + &
                 rxt(k,390)*y(k,164) +rxt(k,395)*y(k,166) +rxt(k,397)*y(k,169) + &
                 rxt(k,401)*y(k,172) +rxt(k,404)*y(k,183) +rxt(k,407)*y(k,185) + &
                 rxt(k,410)*y(k,191) +rxt(k,417)*y(k,196) +rxt(k,423)*y(k,199) + &
                 rxt(k,426)*y(k,201) +1.600_r8*rxt(k,437)*y(k,188) + &
                 .900_r8*rxt(k,442)*y(k,194) +.800_r8*rxt(k,447)*y(k,195))*y(k,116) &
                  + (rxt(k,18) +rxt(k,148)*y(k,176) +rxt(k,150)*y(k,121) + &
                 rxt(k,151)*y(k,190) +rxt(k,317)*y(k,14) +rxt(k,325)*y(k,182) + &
                 rxt(k,336)*y(k,184) +rxt(k,359)*y(k,178) +rxt(k,365)*y(k,179) + &
                 rxt(k,373)*y(k,96) +rxt(k,384)*y(k,198) + &
                 2.000_r8*rxt(k,438)*y(k,188))*y(k,118) + (rxt(k,139)*y(k,85) + &
                 rxt(k,307)*y(k,119) +rxt(k,346)*y(k,1) +.700_r8*rxt(k,366)*y(k,94) + &
                 rxt(k,444)*y(k,148))*y(k,190) + (rxt(k,11) +rxt(k,162))*y(k,85) &
                  + (rxt(k,54) +rxt(k,345))*y(k,104) + (rxt(k,13) +rxt(k,163)) &
                 *y(k,109) + (.600_r8*rxt(k,60) +rxt(k,296))*y(k,127) +rxt(k,19) &
                 *y(k,1) +rxt(k,76)*y(k,18) +rxt(k,95)*y(k,58) +rxt(k,9)*y(k,84) &
                  +rxt(k,45)*y(k,88) +rxt(k,48)*y(k,97) +rxt(k,56)*y(k,111) +rxt(k,57) &
                 *y(k,119) +rxt(k,58)*y(k,120) +rxt(k,59)*y(k,126) +rxt(k,420) &
                 *y(k,128) +rxt(k,66)*y(k,148) +.500_r8*rxt(k,435)*y(k,188)*y(k,171)
         loss(k,170) = (rxt(k,429)* y(k,4) +rxt(k,317)* y(k,14) +rxt(k,297)* y(k,27) &
                  +rxt(k,245)* y(k,40) +rxt(k,278)* y(k,43) +rxt(k,304)* y(k,47) &
                  +rxt(k,452)* y(k,65) +rxt(k,368)* y(k,93) +rxt(k,373)* y(k,96) &
                  +rxt(k,432)* y(k,105) +rxt(k,149)* y(k,116) +rxt(k,159)* y(k,117) &
                  +rxt(k,150)* y(k,121) +rxt(k,449)* y(k,150) +rxt(k,148)* y(k,176) &
                  +rxt(k,359)* y(k,178) +rxt(k,365)* y(k,179) +rxt(k,325)* y(k,182) &
                  +rxt(k,336)* y(k,184) +rxt(k,438)* y(k,188) +rxt(k,151)* y(k,190) &
                  +rxt(k,384)* y(k,198) + rxt(k,17) + rxt(k,18) + rxt(k,480) &
                  + het_rates(k,118))* y(k,118)
         prod(k,170) = (rxt(k,94) +rxt(k,177)*y(k,54) +rxt(k,179)*y(k,121) + &
                 rxt(k,180)*y(k,190))*y(k,58) + (rxt(k,13) +rxt(k,14) +rxt(k,163)) &
                 *y(k,109) + (rxt(k,161)*y(k,84) +rxt(k,293)*y(k,127) + &
                 .500_r8*rxt(k,341)*y(k,104))*y(k,190) + (rxt(k,77) + &
                 rxt(k,210)*y(k,121))*y(k,18) + (rxt(k,146)*y(k,122) + &
                 rxt(k,147)*y(k,121))*y(k,117) +rxt(k,224)*y(k,84)*y(k,68) +rxt(k,10) &
                 *y(k,85) +.400_r8*rxt(k,60)*y(k,127)
         loss(k,121) = (rxt(k,307)* y(k,190) + rxt(k,57) + het_rates(k,119))* y(k,119)
         prod(k,121) = (.500_r8*rxt(k,367)*y(k,95) +rxt(k,374)*y(k,97) + &
                 rxt(k,378)*y(k,110) +rxt(k,379)*y(k,111))*y(k,190) &
                  +rxt(k,297)*y(k,118)*y(k,27)
         loss(k,66) = (rxt(k,439)* y(k,190) + rxt(k,58) + rxt(k,481) &
                  + het_rates(k,120))* y(k,120)
         prod(k,66) =rxt(k,436)*y(k,188)*y(k,176)
         loss(k,162) = (rxt(k,211)* y(k,17) +rxt(k,210)* y(k,18) +rxt(k,246)* y(k,40) &
                  +rxt(k,181)* y(k,57) +rxt(k,179)* y(k,58) +rxt(k,122)* y(k,72) &
                  +rxt(k,123)* y(k,74) +rxt(k,213)* y(k,76) +rxt(k,184)* y(k,80) &
                  +rxt(k,215)* y(k,86) +rxt(k,187)* y(k,87) +rxt(k,155)* y(k,116) &
                  + (rxt(k,145) +rxt(k,147))* y(k,117) +rxt(k,150)* y(k,118) &
                  + 2._r8*rxt(k,120)* y(k,121) +rxt(k,119)* y(k,122) +rxt(k,454) &
                 * y(k,125) +rxt(k,128)* y(k,176) +rxt(k,134)* y(k,190) + rxt(k,121) &
                  + het_rates(k,121))* y(k,121)
         prod(k,162) = (rxt(k,144) +rxt(k,140)*y(k,116) +rxt(k,141)*y(k,117))*y(k,107) &
                  + (rxt(k,111) +rxt(k,462))*y(k,136) + (rxt(k,116) +rxt(k,117)) &
                 *y(k,189) +rxt(k,75)*y(k,17) +rxt(k,93)*y(k,57) +rxt(k,126)*y(k,176) &
                 *y(k,71) +rxt(k,14)*y(k,109) +rxt(k,15)*y(k,116) +rxt(k,16)*y(k,117) &
                  +rxt(k,18)*y(k,118) +rxt(k,8)*y(k,122) +rxt(k,107)*y(k,124) &
                  +rxt(k,456)*y(k,134) +rxt(k,112)*y(k,137) +rxt(k,113)*y(k,138) &
                  +rxt(k,136)*y(k,190)*y(k,190) +rxt(k,3)*y(k,202)
         loss(k,166) = (rxt(k,430)* y(k,4) +rxt(k,202)* y(k,15) +rxt(k,269)* y(k,23) &
                  +rxt(k,298)* y(k,27) +rxt(k,170)* y(k,54) +rxt(k,130)* y(k,71) &
                  +rxt(k,375)* y(k,93) +rxt(k,328)* y(k,100) +rxt(k,433)* y(k,105) &
                  +rxt(k,342)* y(k,106) +rxt(k,154)* y(k,116) +rxt(k,146)* y(k,117) &
                  +rxt(k,119)* y(k,121) +rxt(k,413)* y(k,129) +rxt(k,457)* y(k,134) &
                  +rxt(k,463)* y(k,136) +rxt(k,129)* y(k,176) +rxt(k,118)* y(k,189) &
                  +rxt(k,135)* y(k,190) + rxt(k,7) + rxt(k,8) + het_rates(k,122)) &
                 * y(k,122)
         prod(k,166) = (.150_r8*rxt(k,283)*y(k,170) +.150_r8*rxt(k,333)*y(k,184)) &
                 *y(k,176) +rxt(k,121)*y(k,121)
         loss(k,56) = (rxt(k,464)* y(k,136) + rxt(k,107) + het_rates(k,124))* y(k,124)
         prod(k,56) = (rxt(k,174)*y(k,57) +rxt(k,204)*y(k,17))*y(k,57)
         loss(k,60) = (rxt(k,454)* y(k,121) +rxt(k,455)* y(k,190) + rxt(k,110) &
                  + het_rates(k,125))* y(k,125)
         prod(k,60) = 0._r8
         loss(k,41) = ( + rxt(k,59) + rxt(k,482) + het_rates(k,126))* y(k,126)
         prod(k,41) =rxt(k,321)*y(k,190)*y(k,88) +.100_r8*rxt(k,442)*y(k,194)*y(k,116)
         loss(k,77) = (rxt(k,293)* y(k,190) + rxt(k,60) + rxt(k,296) &
                  + het_rates(k,127))* y(k,127)
         prod(k,77) =rxt(k,295)*y(k,170)*y(k,117)
         loss(k,24) = ( + rxt(k,420) + het_rates(k,128))* y(k,128)
         prod(k,24) =rxt(k,415)*y(k,161)*y(k,117)
         loss(k,78) = (rxt(k,412)* y(k,117) +rxt(k,413)* y(k,122) + het_rates(k,129)) &
                 * y(k,129)
         prod(k,78) = (.070_r8*rxt(k,399)*y(k,64) +.060_r8*rxt(k,411)*y(k,130) + &
                 .070_r8*rxt(k,427)*y(k,157))*y(k,190) +rxt(k,31)*y(k,30) &
                  +rxt(k,397)*y(k,169)*y(k,116)
         loss(k,28) = (rxt(k,411)* y(k,190) + het_rates(k,130))* y(k,130)
         prod(k,28) =.530_r8*rxt(k,388)*y(k,190)*y(k,5)
         loss(k,57) = (rxt(k,414)* y(k,190) + rxt(k,61) + het_rates(k,131))* y(k,131)
         prod(k,57) =rxt(k,409)*y(k,191)*y(k,176)
         loss(k,87) = (rxt(k,310)* y(k,190) + rxt(k,62) + het_rates(k,132))* y(k,132)
         prod(k,87) =rxt(k,308)*y(k,192)*y(k,176)
         loss(k,68) = (rxt(k,314)* y(k,190) + rxt(k,63) + het_rates(k,133))* y(k,133)
         prod(k,68) =.850_r8*rxt(k,312)*y(k,193)*y(k,176)
         loss(k,82) = (rxt(k,457)* y(k,122) +rxt(k,460)* y(k,190) + rxt(k,456) &
                  + het_rates(k,134))* y(k,134)
         prod(k,82) =rxt(k,110)*y(k,125) +rxt(k,111)*y(k,136)
         loss(k,147) = (rxt(k,458)* y(k,17) +rxt(k,459)* y(k,57) +rxt(k,461)* y(k,117) &
                  +rxt(k,463)* y(k,122) +rxt(k,464)* y(k,124) +rxt(k,465)* y(k,190) &
                  + rxt(k,111) + rxt(k,462) + het_rates(k,136))* y(k,136)
         prod(k,147) = (rxt(k,456) +rxt(k,457)*y(k,122) +rxt(k,460)*y(k,190))*y(k,134) &
                  +rxt(k,454)*y(k,125)*y(k,121) +rxt(k,112)*y(k,137)
         loss(k,120) = (rxt(k,467)* y(k,190) + rxt(k,112) + het_rates(k,137)) &
                 * y(k,137)
         prod(k,120) = (rxt(k,462) +rxt(k,458)*y(k,17) +rxt(k,459)*y(k,57) + &
                 rxt(k,461)*y(k,117) +rxt(k,463)*y(k,122) +rxt(k,464)*y(k,124) + &
                 rxt(k,465)*y(k,190))*y(k,136) + (rxt(k,452)*y(k,118) + &
                 rxt(k,453)*y(k,190) +.500_r8*rxt(k,466)*y(k,190))*y(k,65) &
                  +rxt(k,455)*y(k,190)*y(k,125) +rxt(k,113)*y(k,138)
         loss(k,45) = (rxt(k,468)* y(k,202) + rxt(k,113) + het_rates(k,138))* y(k,138)
         prod(k,45) =rxt(k,109)*y(k,75) +rxt(k,467)*y(k,190)*y(k,137)
         loss(k,2) = ( + het_rates(k,139))* y(k,139)
         prod(k,2) = (.1279005_r8*rxt(k,486)*y(k,163) + &
                 .0097005_r8*rxt(k,491)*y(k,165) +.0003005_r8*rxt(k,494)*y(k,180) + &
                 .1056005_r8*rxt(k,498)*y(k,181) +.0245005_r8*rxt(k,502)*y(k,187) + &
                 .0154005_r8*rxt(k,508)*y(k,197) +.0063005_r8*rxt(k,511)*y(k,200)) &
                 *y(k,116) + (.2202005_r8*rxt(k,485)*y(k,163) + &
                 .0023005_r8*rxt(k,490)*y(k,165) +.0031005_r8*rxt(k,493)*y(k,180) + &
                 .2381005_r8*rxt(k,497)*y(k,181) +.0508005_r8*rxt(k,501)*y(k,187) + &
                 .1364005_r8*rxt(k,507)*y(k,197) +.1677005_r8*rxt(k,510)*y(k,200)) &
                 *y(k,176) + (.2202005_r8*rxt(k,487)*y(k,4) + &
                 .0508005_r8*rxt(k,503)*y(k,105))*y(k,122) +rxt(k,470)*y(k,70) &
                  +.5931005_r8*rxt(k,505)*y(k,190)*y(k,145)
         loss(k,3) = ( + het_rates(k,140))* y(k,140)
         prod(k,3) = (.1792005_r8*rxt(k,486)*y(k,163) + &
                 .0034005_r8*rxt(k,491)*y(k,165) +.0003005_r8*rxt(k,494)*y(k,180) + &
                 .1026005_r8*rxt(k,498)*y(k,181) +.0082005_r8*rxt(k,502)*y(k,187) + &
                 .0452005_r8*rxt(k,508)*y(k,197) +.0237005_r8*rxt(k,511)*y(k,200)) &
                 *y(k,116) + (.2067005_r8*rxt(k,485)*y(k,163) + &
                 .0008005_r8*rxt(k,490)*y(k,165) +.0035005_r8*rxt(k,493)*y(k,180) + &
                 .1308005_r8*rxt(k,497)*y(k,181) +.1149005_r8*rxt(k,501)*y(k,187) + &
                 .0101005_r8*rxt(k,507)*y(k,197) +.0174005_r8*rxt(k,510)*y(k,200)) &
                 *y(k,176) + (.2067005_r8*rxt(k,487)*y(k,4) + &
                 .1149005_r8*rxt(k,503)*y(k,105))*y(k,122) &
                  +.1534005_r8*rxt(k,505)*y(k,190)*y(k,145)
         loss(k,4) = ( + het_rates(k,141))* y(k,141)
         prod(k,4) = (.0676005_r8*rxt(k,486)*y(k,163) + &
                 .1579005_r8*rxt(k,491)*y(k,165) +.0073005_r8*rxt(k,494)*y(k,180) + &
                 .0521005_r8*rxt(k,498)*y(k,181) +.0772005_r8*rxt(k,502)*y(k,187) + &
                 .0966005_r8*rxt(k,508)*y(k,197) +.0025005_r8*rxt(k,511)*y(k,200)) &
                 *y(k,116) + (.0653005_r8*rxt(k,485)*y(k,163) + &
                 .0843005_r8*rxt(k,490)*y(k,165) +.0003005_r8*rxt(k,493)*y(k,180) + &
                 .0348005_r8*rxt(k,497)*y(k,181) +.0348005_r8*rxt(k,501)*y(k,187) + &
                 .0763005_r8*rxt(k,507)*y(k,197) +.086_r8*rxt(k,510)*y(k,200)) &
                 *y(k,176) + (.0653005_r8*rxt(k,487)*y(k,4) + &
                 .0348005_r8*rxt(k,503)*y(k,105))*y(k,122) &
                  +.0459005_r8*rxt(k,505)*y(k,190)*y(k,145)
         loss(k,5) = ( + het_rates(k,142))* y(k,142)
         prod(k,5) = (.079_r8*rxt(k,486)*y(k,163) +.0059005_r8*rxt(k,491)*y(k,165) + &
                 .0057005_r8*rxt(k,494)*y(k,180) +.0143005_r8*rxt(k,498)*y(k,181) + &
                 .0332005_r8*rxt(k,502)*y(k,187) +.0073005_r8*rxt(k,508)*y(k,197) + &
                 .011_r8*rxt(k,511)*y(k,200))*y(k,116) &
                  + (.1284005_r8*rxt(k,485)*y(k,163) + &
                 .0443005_r8*rxt(k,490)*y(k,165) +.0271005_r8*rxt(k,493)*y(k,180) + &
                 .0076005_r8*rxt(k,497)*y(k,181) +.0554005_r8*rxt(k,501)*y(k,187) + &
                 .2157005_r8*rxt(k,507)*y(k,197) +.0512005_r8*rxt(k,510)*y(k,200)) &
                 *y(k,176) + (.1749305_r8*rxt(k,484)*y(k,4) + &
                 .0590245_r8*rxt(k,492)*y(k,93) +.1749305_r8*rxt(k,500)*y(k,105)) &
                 *y(k,118) + (.1284005_r8*rxt(k,487)*y(k,4) + &
                 .0033005_r8*rxt(k,495)*y(k,93) +.0554005_r8*rxt(k,503)*y(k,105)) &
                 *y(k,122) +.0085005_r8*rxt(k,505)*y(k,190)*y(k,145)
         loss(k,6) = ( + het_rates(k,143))* y(k,143)
         prod(k,6) = (.1254005_r8*rxt(k,486)*y(k,163) + &
                 .0536005_r8*rxt(k,491)*y(k,165) +.0623005_r8*rxt(k,494)*y(k,180) + &
                 .0166005_r8*rxt(k,498)*y(k,181) +.130_r8*rxt(k,502)*y(k,187) + &
                 .238_r8*rxt(k,508)*y(k,197) +.1185005_r8*rxt(k,511)*y(k,200)) &
                 *y(k,116) + (.114_r8*rxt(k,485)*y(k,163) + &
                 .1621005_r8*rxt(k,490)*y(k,165) +.0474005_r8*rxt(k,493)*y(k,180) + &
                 .0113005_r8*rxt(k,497)*y(k,181) +.1278005_r8*rxt(k,501)*y(k,187) + &
                 .0738005_r8*rxt(k,507)*y(k,197) +.1598005_r8*rxt(k,510)*y(k,200)) &
                 *y(k,176) + (.5901905_r8*rxt(k,484)*y(k,4) + &
                 .0250245_r8*rxt(k,492)*y(k,93) +.5901905_r8*rxt(k,500)*y(k,105)) &
                 *y(k,118) + (.114_r8*rxt(k,487)*y(k,4) + &
                 .1278005_r8*rxt(k,503)*y(k,105))*y(k,122) &
                  +.0128005_r8*rxt(k,505)*y(k,190)*y(k,145)
         loss(k,7) = (rxt(k,505)* y(k,190) + het_rates(k,145))* y(k,145)
         prod(k,7) = 0._r8
         loss(k,33) = ( + rxt(k,64) + het_rates(k,146))* y(k,146)
         prod(k,33) = (.100_r8*rxt(k,419)*y(k,153) +.230_r8*rxt(k,421)*y(k,155)) &
                 *y(k,190)
         loss(k,95) = (rxt(k,443)* y(k,190) + rxt(k,65) + het_rates(k,147))* y(k,147)
         prod(k,95) =rxt(k,441)*y(k,194)*y(k,176)
         loss(k,98) = (rxt(k,444)* y(k,190) + rxt(k,66) + rxt(k,483) &
                  + het_rates(k,148))* y(k,148)
         prod(k,98) = (.200_r8*rxt(k,437)*y(k,188) +.200_r8*rxt(k,447)*y(k,195)) &
                 *y(k,116) +.500_r8*rxt(k,435)*y(k,188)*y(k,171)
         loss(k,79) = (rxt(k,448)* y(k,190) + rxt(k,67) + het_rates(k,149))* y(k,149)
         prod(k,79) =rxt(k,446)*y(k,195)*y(k,176)
         loss(k,130) = (rxt(k,449)* y(k,118) +rxt(k,450)* y(k,190) + rxt(k,68) &
                  + het_rates(k,150))* y(k,150)
         prod(k,130) = (.500_r8*rxt(k,435)*y(k,171) +.800_r8*rxt(k,437)*y(k,116) + &
                 rxt(k,438)*y(k,118))*y(k,188) + (.330_r8*rxt(k,430)*y(k,4) + &
                 .330_r8*rxt(k,433)*y(k,105))*y(k,122) + (rxt(k,66) + &
                 rxt(k,444)*y(k,190))*y(k,148) + (rxt(k,445)*y(k,171) + &
                 .800_r8*rxt(k,447)*y(k,116))*y(k,195) +rxt(k,58)*y(k,120) +rxt(k,67) &
                 *y(k,149)
         loss(k,133) = (rxt(k,451)* y(k,190) + rxt(k,69) + het_rates(k,151))* y(k,151)
         prod(k,133) = (.300_r8*rxt(k,430)*y(k,4) +.300_r8*rxt(k,433)*y(k,105)) &
                 *y(k,122) + (rxt(k,440)*y(k,171) +.900_r8*rxt(k,442)*y(k,116)) &
                 *y(k,194) +rxt(k,65)*y(k,147) +rxt(k,68)*y(k,150)
         loss(k,96) = (rxt(k,418)* y(k,190) + rxt(k,70) + het_rates(k,152))* y(k,152)
         prod(k,96) =rxt(k,416)*y(k,196)*y(k,176)
         loss(k,31) = (rxt(k,419)* y(k,190) + het_rates(k,153))* y(k,153)
         prod(k,31) = 0._r8
         loss(k,34) = (rxt(k,385)* y(k,190) + rxt(k,71) + het_rates(k,154))* y(k,154)
         prod(k,34) =rxt(k,382)*y(k,198)*y(k,176)
         loss(k,35) = (rxt(k,421)* y(k,190) + het_rates(k,155))* y(k,155)
         prod(k,35) = 0._r8
         loss(k,103) = (rxt(k,424)* y(k,190) + rxt(k,72) + het_rates(k,156))* y(k,156)
         prod(k,103) =rxt(k,422)*y(k,199)*y(k,176)
         loss(k,36) = (rxt(k,427)* y(k,190) + het_rates(k,157))* y(k,157)
         prod(k,36) =.150_r8*rxt(k,421)*y(k,190)*y(k,155)
         loss(k,72) = (rxt(k,428)* y(k,190) + rxt(k,73) + het_rates(k,158))* y(k,158)
         prod(k,72) =rxt(k,425)*y(k,201)*y(k,176)
         loss(k,85) = (rxt(k,387)* y(k,116) +rxt(k,415)* y(k,117) +rxt(k,386) &
                 * y(k,176) + het_rates(k,161))* y(k,161)
         prod(k,85) =rxt(k,392)*y(k,190)*y(k,20) +rxt(k,420)*y(k,128)
         loss(k,126) = ((rxt(k,348) +rxt(k,349))* y(k,116) +rxt(k,347)* y(k,176) &
                  + het_rates(k,162))* y(k,162)
         prod(k,126) = (rxt(k,350)*y(k,2) +rxt(k,351)*y(k,13))*y(k,190)
         loss(k,8) = (rxt(k,486)* y(k,116) +rxt(k,485)* y(k,176) + het_rates(k,163)) &
                 * y(k,163)
         prod(k,8) =rxt(k,488)*y(k,190)*y(k,4)
         loss(k,80) = (rxt(k,390)* y(k,116) +rxt(k,389)* y(k,176) + het_rates(k,164)) &
                 * y(k,164)
         prod(k,80) = (.350_r8*rxt(k,388)*y(k,5) +rxt(k,391)*y(k,6))*y(k,190)
         loss(k,9) = (rxt(k,491)* y(k,116) +rxt(k,490)* y(k,176) + het_rates(k,165)) &
                 * y(k,165)
         prod(k,9) =rxt(k,489)*y(k,190)*y(k,5)
         loss(k,73) = (rxt(k,395)* y(k,116) +rxt(k,393)* y(k,176) + het_rates(k,166)) &
                 * y(k,166)
         prod(k,73) = (rxt(k,394)*y(k,21) +.070_r8*rxt(k,419)*y(k,153) + &
                 .060_r8*rxt(k,421)*y(k,155))*y(k,190)
         loss(k,117) = (rxt(k,273)* y(k,116) + 2._r8*rxt(k,270)* y(k,167) +rxt(k,271) &
                 * y(k,171) +rxt(k,272)* y(k,176) + het_rates(k,167))* y(k,167)
         prod(k,117) = (rxt(k,276)*y(k,54) +rxt(k,277)*y(k,190))*y(k,26) &
                  +.500_r8*rxt(k,275)*y(k,190)*y(k,25) +rxt(k,52)*y(k,102)
         loss(k,114) = (rxt(k,301)* y(k,116) +rxt(k,299)* y(k,171) +rxt(k,300) &
                 * y(k,176) + het_rates(k,168))* y(k,168)
         prod(k,114) = (rxt(k,302)*y(k,28) +rxt(k,303)*y(k,29))*y(k,190)
         loss(k,99) = (rxt(k,397)* y(k,116) +rxt(k,396)* y(k,176) + het_rates(k,169)) &
                 * y(k,169)
         prod(k,99) = (.400_r8*rxt(k,386)*y(k,176) +rxt(k,387)*y(k,116))*y(k,161) &
                  +rxt(k,398)*y(k,190)*y(k,30) +rxt(k,413)*y(k,129)*y(k,122)
         loss(k,153) = (rxt(k,369)* y(k,96) +rxt(k,284)* y(k,116) +rxt(k,295) &
                 * y(k,117) + 2._r8*rxt(k,281)* y(k,170) +rxt(k,282)* y(k,171) &
                  +rxt(k,283)* y(k,176) +rxt(k,355)* y(k,178) +rxt(k,360)* y(k,179) &
                  +rxt(k,322)* y(k,182) +rxt(k,380)* y(k,198) + het_rates(k,170)) &
                 * y(k,170)
         prod(k,153) = (.100_r8*rxt(k,328)*y(k,100) +.280_r8*rxt(k,342)*y(k,106) + &
                 .080_r8*rxt(k,375)*y(k,93) +.060_r8*rxt(k,430)*y(k,4) + &
                 .060_r8*rxt(k,433)*y(k,105))*y(k,122) + (rxt(k,332)*y(k,171) + &
                 .450_r8*rxt(k,333)*y(k,176) +2.000_r8*rxt(k,334)*y(k,184) + &
                 rxt(k,335)*y(k,116) +rxt(k,336)*y(k,118))*y(k,184) &
                  + (.530_r8*rxt(k,322)*y(k,170) +.260_r8*rxt(k,323)*y(k,171) + &
                 .530_r8*rxt(k,325)*y(k,118) +.530_r8*rxt(k,326)*y(k,116))*y(k,182) &
                  + (rxt(k,279)*y(k,43) +.500_r8*rxt(k,286)*y(k,49) + &
                 rxt(k,305)*y(k,47) +.650_r8*rxt(k,451)*y(k,151))*y(k,190) &
                  + (.300_r8*rxt(k,311)*y(k,171) +.150_r8*rxt(k,312)*y(k,176) + &
                 rxt(k,313)*y(k,116))*y(k,193) + (rxt(k,36) +rxt(k,304)*y(k,118)) &
                 *y(k,47) + (.600_r8*rxt(k,60) +rxt(k,296))*y(k,127) &
                  + (.200_r8*rxt(k,337)*y(k,176) +rxt(k,338)*y(k,116))*y(k,186) &
                  +.130_r8*rxt(k,23)*y(k,8) +rxt(k,27)*y(k,12) +rxt(k,278)*y(k,118) &
                 *y(k,43) +rxt(k,35)*y(k,46) +.330_r8*rxt(k,45)*y(k,88) +rxt(k,47) &
                 *y(k,90) +1.340_r8*rxt(k,50)*y(k,100) +rxt(k,52)*y(k,102) +rxt(k,53) &
                 *y(k,103) +.300_r8*rxt(k,55)*y(k,106) +rxt(k,57)*y(k,119) +rxt(k,63) &
                 *y(k,133) +.500_r8*rxt(k,64)*y(k,146) +.650_r8*rxt(k,69)*y(k,151)
         loss(k,157) = (rxt(k,171)* y(k,57) +rxt(k,370)* y(k,96) +rxt(k,251)* y(k,116) &
                  +rxt(k,271)* y(k,167) +rxt(k,299)* y(k,168) +rxt(k,282)* y(k,170) &
                  + 2._r8*(rxt(k,248) +rxt(k,249))* y(k,171) +rxt(k,250)* y(k,176) &
                  +rxt(k,356)* y(k,178) +rxt(k,361)* y(k,179) +rxt(k,323)* y(k,182) &
                  +rxt(k,332)* y(k,184) +rxt(k,435)* y(k,188) +rxt(k,311)* y(k,193) &
                  +rxt(k,440)* y(k,194) +rxt(k,445)* y(k,195) +rxt(k,381)* y(k,198) &
                  + het_rates(k,171))* y(k,171)
         prod(k,157) = (2.000_r8*rxt(k,281)*y(k,170) +.900_r8*rxt(k,282)*y(k,171) + &
                 .450_r8*rxt(k,283)*y(k,176) +rxt(k,284)*y(k,116) + &
                 rxt(k,322)*y(k,182) +rxt(k,331)*y(k,184) +rxt(k,355)*y(k,178) + &
                 rxt(k,360)*y(k,179) +rxt(k,369)*y(k,96) +rxt(k,380)*y(k,198)) &
                 *y(k,170) + (rxt(k,165)*y(k,54) +rxt(k,221)*y(k,68) + &
                 rxt(k,254)*y(k,190) +rxt(k,261)*y(k,189))*y(k,52) &
                  + (.830_r8*rxt(k,401)*y(k,172) +.170_r8*rxt(k,407)*y(k,185)) &
                 *y(k,116) + (.280_r8*rxt(k,298)*y(k,27) +.050_r8*rxt(k,375)*y(k,93)) &
                 *y(k,122) + (.330_r8*rxt(k,400)*y(k,172) + &
                 .070_r8*rxt(k,406)*y(k,185))*y(k,176) + (.700_r8*rxt(k,253)*y(k,51) + &
                 rxt(k,285)*y(k,48))*y(k,190) +rxt(k,34)*y(k,43) +rxt(k,35)*y(k,46) &
                  +rxt(k,37)*y(k,49) +.300_r8*rxt(k,55)*y(k,106) +.400_r8*rxt(k,60) &
                 *y(k,127)
         loss(k,109) = (rxt(k,401)* y(k,116) +rxt(k,402)* y(k,117) +rxt(k,400) &
                 * y(k,176) + het_rates(k,172))* y(k,172)
         prod(k,109) =.600_r8*rxt(k,25)*y(k,10)
         loss(k,91) = ((rxt(k,319) +rxt(k,320))* y(k,116) + het_rates(k,173)) &
                 * y(k,173)
         prod(k,91) =rxt(k,318)*y(k,190)*y(k,14)
         loss(k,46) = ( + rxt(k,289) + rxt(k,290) + het_rates(k,174))* y(k,174)
         prod(k,46) =rxt(k,42)*y(k,67) +.750_r8*rxt(k,288)*y(k,175)*y(k,116)
         loss(k,104) = (rxt(k,288)* y(k,116) +rxt(k,287)* y(k,176) + het_rates(k,175)) &
                 * y(k,175)
         prod(k,104) =rxt(k,294)*y(k,190)*y(k,23)
         loss(k,163) = (rxt(k,201)* y(k,15) +rxt(k,207)* y(k,17) +rxt(k,244)* y(k,40) &
                  + (rxt(k,168) +rxt(k,169))* y(k,54) +rxt(k,175)* y(k,57) &
                  + (rxt(k,124) +rxt(k,125) +rxt(k,126))* y(k,71) +rxt(k,371)* y(k,96) &
                  +rxt(k,153)* y(k,116) +rxt(k,158)* y(k,117) +rxt(k,148)* y(k,118) &
                  +rxt(k,128)* y(k,121) +rxt(k,129)* y(k,122) +rxt(k,386)* y(k,161) &
                  +rxt(k,347)* y(k,162) +rxt(k,389)* y(k,164) +rxt(k,393)* y(k,166) &
                  +rxt(k,272)* y(k,167) +rxt(k,300)* y(k,168) +rxt(k,396)* y(k,169) &
                  +rxt(k,283)* y(k,170) +rxt(k,250)* y(k,171) +rxt(k,400)* y(k,172) &
                  +rxt(k,287)* y(k,175) + 2._r8*rxt(k,138)* y(k,176) +rxt(k,258) &
                 * y(k,177) +rxt(k,357)* y(k,178) +rxt(k,362)* y(k,179) +rxt(k,324) &
                 * y(k,182) +rxt(k,403)* y(k,183) +rxt(k,333)* y(k,184) +rxt(k,406) &
                 * y(k,185) +rxt(k,337)* y(k,186) +rxt(k,436)* y(k,188) +rxt(k,133) &
                 * y(k,190) +rxt(k,409)* y(k,191) +rxt(k,308)* y(k,192) +rxt(k,312) &
                 * y(k,193) +rxt(k,441)* y(k,194) +rxt(k,446)* y(k,195) +rxt(k,416) &
                 * y(k,196) +rxt(k,382)* y(k,198) +rxt(k,422)* y(k,199) +rxt(k,425) &
                 * y(k,201) + rxt(k,471) + het_rates(k,176))* y(k,176)
         prod(k,163) = (rxt(k,230)*y(k,41) +rxt(k,233)*y(k,44) +rxt(k,132)*y(k,74) + &
                 rxt(k,135)*y(k,122) +rxt(k,151)*y(k,118) +rxt(k,182)*y(k,57) + &
                 rxt(k,212)*y(k,17) +rxt(k,252)*y(k,50) +rxt(k,255)*y(k,60) + &
                 rxt(k,256)*y(k,81) +rxt(k,257)*y(k,82) +.350_r8*rxt(k,267)*y(k,22) + &
                 rxt(k,274)*y(k,24) +rxt(k,280)*y(k,45) +rxt(k,291)*y(k,69) + &
                 rxt(k,292)*y(k,70) +rxt(k,306)*y(k,90) +rxt(k,321)*y(k,88) + &
                 .200_r8*rxt(k,330)*y(k,101) +.500_r8*rxt(k,341)*y(k,104) + &
                 .300_r8*rxt(k,366)*y(k,94) +rxt(k,367)*y(k,95) +rxt(k,374)*y(k,97) + &
                 rxt(k,378)*y(k,110) +rxt(k,379)*y(k,111) +.650_r8*rxt(k,388)*y(k,5) + &
                 .730_r8*rxt(k,399)*y(k,64) +.800_r8*rxt(k,411)*y(k,130) + &
                 .280_r8*rxt(k,419)*y(k,153) +.380_r8*rxt(k,421)*y(k,155) + &
                 .630_r8*rxt(k,427)*y(k,157) +.200_r8*rxt(k,451)*y(k,151) + &
                 .500_r8*rxt(k,466)*y(k,65) +rxt(k,467)*y(k,137))*y(k,190) &
                  + (rxt(k,251)*y(k,171) +rxt(k,260)*y(k,177) +rxt(k,273)*y(k,167) + &
                 .250_r8*rxt(k,288)*y(k,175) +rxt(k,301)*y(k,168) + &
                 rxt(k,309)*y(k,192) +rxt(k,319)*y(k,173) + &
                 .470_r8*rxt(k,326)*y(k,182) +rxt(k,348)*y(k,162) + &
                 .920_r8*rxt(k,358)*y(k,178) +.920_r8*rxt(k,364)*y(k,179) + &
                 rxt(k,372)*y(k,96) +rxt(k,383)*y(k,198) +rxt(k,390)*y(k,164) + &
                 rxt(k,395)*y(k,166) +.170_r8*rxt(k,401)*y(k,172) + &
                 .400_r8*rxt(k,404)*y(k,183) +.830_r8*rxt(k,407)*y(k,185) + &
                 rxt(k,410)*y(k,191) +rxt(k,417)*y(k,196) +rxt(k,423)*y(k,199) + &
                 rxt(k,426)*y(k,201) +.900_r8*rxt(k,442)*y(k,194) + &
                 .800_r8*rxt(k,447)*y(k,195))*y(k,116) + (rxt(k,171)*y(k,57) + &
                 2.000_r8*rxt(k,248)*y(k,171) +rxt(k,271)*y(k,167) + &
                 .900_r8*rxt(k,282)*y(k,170) +rxt(k,299)*y(k,168) + &
                 .300_r8*rxt(k,311)*y(k,193) +.730_r8*rxt(k,323)*y(k,182) + &
                 rxt(k,332)*y(k,184) +rxt(k,356)*y(k,178) +rxt(k,361)*y(k,179) + &
                 1.200_r8*rxt(k,370)*y(k,96) +.800_r8*rxt(k,381)*y(k,198) + &
                 .500_r8*rxt(k,435)*y(k,188) +rxt(k,440)*y(k,194) + &
                 rxt(k,445)*y(k,195))*y(k,171) + (.130_r8*rxt(k,269)*y(k,23) + &
                 .280_r8*rxt(k,298)*y(k,27) +.140_r8*rxt(k,328)*y(k,100) + &
                 .280_r8*rxt(k,342)*y(k,106) +.370_r8*rxt(k,375)*y(k,93) + &
                 .570_r8*rxt(k,430)*y(k,4) +.570_r8*rxt(k,433)*y(k,105))*y(k,122) &
                  + (rxt(k,245)*y(k,40) +.470_r8*rxt(k,325)*y(k,182) + &
                 rxt(k,359)*y(k,178) +rxt(k,365)*y(k,179) +rxt(k,373)*y(k,96) + &
                 rxt(k,384)*y(k,198))*y(k,118) + (.470_r8*rxt(k,322)*y(k,182) + &
                 rxt(k,355)*y(k,178) +rxt(k,360)*y(k,179) +rxt(k,369)*y(k,96) + &
                 rxt(k,380)*y(k,198))*y(k,170) + (rxt(k,229)*y(k,41) + &
                 rxt(k,232)*y(k,44) +rxt(k,164)*y(k,40) +rxt(k,167)*y(k,74))*y(k,54) &
                  + (.070_r8*rxt(k,400)*y(k,172) +.160_r8*rxt(k,403)*y(k,183) + &
                 .330_r8*rxt(k,406)*y(k,185))*y(k,176) + (rxt(k,200)*y(k,15) + &
                 rxt(k,246)*y(k,121))*y(k,40) + (rxt(k,11) +rxt(k,162))*y(k,85) &
                  + (1.340_r8*rxt(k,50) +.660_r8*rxt(k,51))*y(k,100) + (rxt(k,289) + &
                 rxt(k,290))*y(k,174) +rxt(k,19)*y(k,1) +.900_r8*rxt(k,20)*y(k,2) &
                  +rxt(k,21)*y(k,6) +1.500_r8*rxt(k,22)*y(k,7) +.560_r8*rxt(k,23) &
                 *y(k,8) +rxt(k,24)*y(k,9) +.600_r8*rxt(k,25)*y(k,10) &
                  +.600_r8*rxt(k,26)*y(k,11) +rxt(k,27)*y(k,12) +rxt(k,28)*y(k,21) &
                  +rxt(k,29)*y(k,25) +rxt(k,30)*y(k,28) +rxt(k,34)*y(k,43) +rxt(k,36) &
                 *y(k,47) +rxt(k,262)*y(k,189)*y(k,52) +2.000_r8*rxt(k,43)*y(k,69) &
                  +2.000_r8*rxt(k,44)*y(k,70) +rxt(k,127)*y(k,71) +rxt(k,123)*y(k,121) &
                 *y(k,74) +.670_r8*rxt(k,45)*y(k,88) +rxt(k,46)*y(k,89) +rxt(k,47) &
                 *y(k,90) +rxt(k,48)*y(k,97) +rxt(k,49)*y(k,98) +rxt(k,56)*y(k,111) &
                  +rxt(k,61)*y(k,131) +rxt(k,62)*y(k,132) +rxt(k,64)*y(k,146) &
                  +rxt(k,65)*y(k,147) +rxt(k,66)*y(k,148) +rxt(k,67)*y(k,149) &
                  +rxt(k,68)*y(k,150) +1.200_r8*rxt(k,69)*y(k,151) +rxt(k,70)*y(k,152) &
                  +rxt(k,72)*y(k,156) +rxt(k,73)*y(k,158) &
                  +1.200_r8*rxt(k,270)*y(k,167)*y(k,167) +rxt(k,259)*y(k,177) &
                  +rxt(k,363)*y(k,179)
         loss(k,74) = (rxt(k,260)* y(k,116) +rxt(k,258)* y(k,176) + rxt(k,259) &
                  + het_rates(k,177))* y(k,177)
         prod(k,74) =rxt(k,244)*y(k,176)*y(k,40)
         loss(k,151) = (rxt(k,358)* y(k,116) +rxt(k,359)* y(k,118) +rxt(k,355) &
                 * y(k,170) +rxt(k,356)* y(k,171) +rxt(k,357)* y(k,176) &
                  + het_rates(k,178))* y(k,178)
         prod(k,151) =.600_r8*rxt(k,376)*y(k,190)*y(k,93)
         loss(k,148) = (rxt(k,364)* y(k,116) +rxt(k,365)* y(k,118) +rxt(k,360) &
                 * y(k,170) +rxt(k,361)* y(k,171) +rxt(k,362)* y(k,176) + rxt(k,363) &
                  + het_rates(k,179))* y(k,179)
         prod(k,148) =.400_r8*rxt(k,376)*y(k,190)*y(k,93)
         loss(k,10) = (rxt(k,494)* y(k,116) +rxt(k,493)* y(k,176) + het_rates(k,180)) &
                 * y(k,180)
         prod(k,10) =rxt(k,496)*y(k,190)*y(k,93)
         loss(k,11) = (rxt(k,498)* y(k,116) +rxt(k,497)* y(k,176) + het_rates(k,181)) &
                 * y(k,181)
         prod(k,11) =rxt(k,499)*y(k,190)*y(k,99)
         loss(k,150) = ((rxt(k,326) +rxt(k,327))* y(k,116) +rxt(k,325)* y(k,118) &
                  +rxt(k,322)* y(k,170) +rxt(k,323)* y(k,171) +rxt(k,324)* y(k,176) &
                  + het_rates(k,182))* y(k,182)
         prod(k,150) = (.500_r8*rxt(k,329)*y(k,100) +.200_r8*rxt(k,330)*y(k,101) + &
                 rxt(k,343)*y(k,106))*y(k,190)
         loss(k,105) = (rxt(k,404)* y(k,116) +rxt(k,405)* y(k,117) +rxt(k,403) &
                 * y(k,176) + het_rates(k,183))* y(k,183)
         prod(k,105) =.600_r8*rxt(k,24)*y(k,9)
         loss(k,152) = (rxt(k,335)* y(k,116) +rxt(k,344)* y(k,117) +rxt(k,336) &
                 * y(k,118) +rxt(k,331)* y(k,170) +rxt(k,332)* y(k,171) +rxt(k,333) &
                 * y(k,176) + 2._r8*rxt(k,334)* y(k,184) + het_rates(k,184))* y(k,184)
         prod(k,152) = (.660_r8*rxt(k,50) +.500_r8*rxt(k,329)*y(k,190))*y(k,100) &
                  + (rxt(k,54) +rxt(k,345))*y(k,104) +.500_r8*rxt(k,330)*y(k,190) &
                 *y(k,101)
         loss(k,123) = (rxt(k,407)* y(k,116) +rxt(k,408)* y(k,117) +rxt(k,406) &
                 * y(k,176) + het_rates(k,185))* y(k,185)
         prod(k,123) =.600_r8*rxt(k,26)*y(k,11)
         loss(k,102) = (rxt(k,338)* y(k,116) +rxt(k,337)* y(k,176) + het_rates(k,186)) &
                 * y(k,186)
         prod(k,102) = (rxt(k,339)*y(k,102) +rxt(k,340)*y(k,103))*y(k,190)
         loss(k,13) = (rxt(k,502)* y(k,116) +rxt(k,501)* y(k,176) + het_rates(k,187)) &
                 * y(k,187)
         prod(k,13) =rxt(k,504)*y(k,190)*y(k,105)
         loss(k,140) = (rxt(k,437)* y(k,116) +rxt(k,438)* y(k,118) +rxt(k,435) &
                 * y(k,171) +rxt(k,436)* y(k,176) + het_rates(k,188))* y(k,188)
         prod(k,140) = (rxt(k,429)*y(k,4) +rxt(k,432)*y(k,105) + &
                 .500_r8*rxt(k,449)*y(k,150))*y(k,118) +rxt(k,439)*y(k,190)*y(k,120)
         loss(k,158) = (rxt(k,189)* y(k,31) +rxt(k,190)* y(k,32) +rxt(k,216)* y(k,33) &
                  +rxt(k,191)* y(k,34) +rxt(k,192)* y(k,35) +rxt(k,193)* y(k,36) &
                  +rxt(k,194)* y(k,37) +rxt(k,195)* y(k,38) +rxt(k,239)* y(k,39) &
                  +rxt(k,240)* y(k,41) + (rxt(k,261) +rxt(k,262) +rxt(k,263))* y(k,52) &
                  +rxt(k,217)* y(k,53) +rxt(k,225)* y(k,62) +rxt(k,226)* y(k,63) &
                  +rxt(k,114)* y(k,72) +rxt(k,218)* y(k,73) + (rxt(k,219) +rxt(k,220)) &
                 * y(k,76) +rxt(k,241)* y(k,77) +rxt(k,242)* y(k,78) +rxt(k,243) &
                 * y(k,79) + (rxt(k,196) +rxt(k,197))* y(k,80) +rxt(k,264)* y(k,81) &
                  + (rxt(k,156) +rxt(k,157))* y(k,108) +rxt(k,118)* y(k,122) &
                  +rxt(k,115)* y(k,202) + rxt(k,116) + rxt(k,117) + het_rates(k,189)) &
                 * y(k,189)
         prod(k,158) =rxt(k,7)*y(k,122) +rxt(k,1)*y(k,202)
         loss(k,159) = (rxt(k,346)* y(k,1) +rxt(k,350)* y(k,2) +rxt(k,431)* y(k,4) &
                  +rxt(k,388)* y(k,5) +rxt(k,391)* y(k,6) +rxt(k,351)* y(k,13) &
                  +rxt(k,318)* y(k,14) +rxt(k,212)* y(k,17) +rxt(k,392)* y(k,20) &
                  +rxt(k,394)* y(k,21) +rxt(k,267)* y(k,22) +rxt(k,294)* y(k,23) &
                  +rxt(k,274)* y(k,24) +rxt(k,275)* y(k,25) +rxt(k,277)* y(k,26) &
                  +rxt(k,315)* y(k,27) +rxt(k,302)* y(k,28) +rxt(k,303)* y(k,29) &
                  +rxt(k,398)* y(k,30) +rxt(k,228)* y(k,39) +rxt(k,247)* y(k,40) &
                  +rxt(k,230)* y(k,41) +rxt(k,231)* y(k,42) +rxt(k,279)* y(k,43) &
                  +rxt(k,233)* y(k,44) +rxt(k,280)* y(k,45) +rxt(k,316)* y(k,46) &
                  +rxt(k,305)* y(k,47) +rxt(k,285)* y(k,48) +rxt(k,286)* y(k,49) &
                  +rxt(k,252)* y(k,50) +rxt(k,253)* y(k,51) +rxt(k,254)* y(k,52) &
                  +rxt(k,235)* y(k,53) + (rxt(k,182) +rxt(k,183))* y(k,57) +rxt(k,180) &
                 * y(k,58) + (rxt(k,255) +rxt(k,265))* y(k,60) +rxt(k,399)* y(k,64) &
                  + (rxt(k,453) +rxt(k,466))* y(k,65) +rxt(k,291)* y(k,69) +rxt(k,292) &
                 * y(k,70) +rxt(k,131)* y(k,72) +rxt(k,132)* y(k,74) +rxt(k,214) &
                 * y(k,76) +rxt(k,236)* y(k,77) +rxt(k,237)* y(k,78) +rxt(k,238) &
                 * y(k,79) +rxt(k,185)* y(k,80) +rxt(k,256)* y(k,81) +rxt(k,257) &
                 * y(k,82) +rxt(k,161)* y(k,84) +rxt(k,139)* y(k,85) +rxt(k,188) &
                 * y(k,87) +rxt(k,321)* y(k,88) +rxt(k,352)* y(k,89) +rxt(k,306) &
                 * y(k,90) +rxt(k,353)* y(k,91) +rxt(k,354)* y(k,92) +rxt(k,376) &
                 * y(k,93) +rxt(k,366)* y(k,94) +rxt(k,367)* y(k,95) +rxt(k,374) &
                 * y(k,97) +rxt(k,377)* y(k,98) +rxt(k,329)* y(k,100) +rxt(k,330) &
                 * y(k,101) +rxt(k,339)* y(k,102) +rxt(k,340)* y(k,103) +rxt(k,341) &
                 * y(k,104) +rxt(k,434)* y(k,105) +rxt(k,343)* y(k,106) +rxt(k,152) &
                 * y(k,107) +rxt(k,378)* y(k,110) +rxt(k,379)* y(k,111) +rxt(k,469) &
                 * y(k,112) +rxt(k,160)* y(k,117) +rxt(k,151)* y(k,118) +rxt(k,307) &
                 * y(k,119) +rxt(k,439)* y(k,120) +rxt(k,134)* y(k,121) +rxt(k,135) &
                 * y(k,122) +rxt(k,455)* y(k,125) +rxt(k,293)* y(k,127) +rxt(k,411) &
                 * y(k,130) +rxt(k,414)* y(k,131) +rxt(k,310)* y(k,132) +rxt(k,314) &
                 * y(k,133) +rxt(k,460)* y(k,134) +rxt(k,465)* y(k,136) +rxt(k,467) &
                 * y(k,137) +rxt(k,443)* y(k,147) +rxt(k,444)* y(k,148) +rxt(k,448) &
                 * y(k,149) +rxt(k,450)* y(k,150) +rxt(k,451)* y(k,151) +rxt(k,418) &
                 * y(k,152) +rxt(k,419)* y(k,153) +rxt(k,385)* y(k,154) +rxt(k,421) &
                 * y(k,155) +rxt(k,424)* y(k,156) +rxt(k,427)* y(k,157) +rxt(k,428) &
                 * y(k,158) +rxt(k,133)* y(k,176) + 2._r8*(rxt(k,136) +rxt(k,137)) &
                 * y(k,190) + het_rates(k,190))* y(k,190)
         prod(k,159) = (2.000_r8*rxt(k,125)*y(k,71) +rxt(k,128)*y(k,121) + &
                 rxt(k,129)*y(k,122) +rxt(k,148)*y(k,118) +rxt(k,153)*y(k,116) + &
                 rxt(k,169)*y(k,54) +.450_r8*rxt(k,283)*y(k,170) + &
                 .150_r8*rxt(k,312)*y(k,193) +.450_r8*rxt(k,333)*y(k,184) + &
                 .200_r8*rxt(k,337)*y(k,186) +.400_r8*rxt(k,386)*y(k,161) + &
                 .400_r8*rxt(k,400)*y(k,172) +.400_r8*rxt(k,406)*y(k,185))*y(k,176) &
                  + (rxt(k,130)*y(k,71) +.130_r8*rxt(k,269)*y(k,23) + &
                 .360_r8*rxt(k,298)*y(k,27) +.240_r8*rxt(k,328)*y(k,100) + &
                 .360_r8*rxt(k,342)*y(k,106) +.320_r8*rxt(k,375)*y(k,93) + &
                 .630_r8*rxt(k,430)*y(k,4) +.630_r8*rxt(k,433)*y(k,105))*y(k,122) &
                  + (rxt(k,122)*y(k,72) +rxt(k,123)*y(k,74) +rxt(k,184)*y(k,80) + &
                 rxt(k,187)*y(k,87) +rxt(k,213)*y(k,76) +rxt(k,215)*y(k,86) + &
                 rxt(k,246)*y(k,40))*y(k,121) + (.300_r8*rxt(k,253)*y(k,51) + &
                 .650_r8*rxt(k,267)*y(k,22) +.500_r8*rxt(k,275)*y(k,25) + &
                 .500_r8*rxt(k,310)*y(k,132) +.100_r8*rxt(k,330)*y(k,101) + &
                 .600_r8*rxt(k,377)*y(k,98) +.500_r8*rxt(k,385)*y(k,154))*y(k,190) &
                  + (rxt(k,261)*y(k,52) +rxt(k,114)*y(k,72) + &
                 2.000_r8*rxt(k,115)*y(k,202) +rxt(k,196)*y(k,80) + &
                 rxt(k,219)*y(k,76) +rxt(k,264)*y(k,81))*y(k,189) + (rxt(k,2) + &
                 rxt(k,223)*y(k,68))*y(k,202) +rxt(k,20)*y(k,2) +rxt(k,21)*y(k,6) &
                  +rxt(k,28)*y(k,21) +rxt(k,29)*y(k,25) +rxt(k,30)*y(k,28) +rxt(k,31) &
                 *y(k,30) +rxt(k,37)*y(k,49) +rxt(k,38)*y(k,51) +rxt(k,42)*y(k,67) &
                  +2.000_r8*rxt(k,4)*y(k,74) +rxt(k,9)*y(k,84) +rxt(k,10)*y(k,85) &
                  +rxt(k,105)*y(k,86) +rxt(k,106)*y(k,87) +rxt(k,46)*y(k,89) &
                  +rxt(k,53)*y(k,103) +.500_r8*rxt(k,479)*y(k,117) +rxt(k,58)*y(k,120) &
                  +rxt(k,61)*y(k,131) +rxt(k,62)*y(k,132) +rxt(k,63)*y(k,133) &
                  +rxt(k,65)*y(k,147) +rxt(k,67)*y(k,149) +rxt(k,70)*y(k,152) &
                  +rxt(k,71)*y(k,154) +rxt(k,72)*y(k,156) +rxt(k,73)*y(k,158)
         loss(k,75) = (rxt(k,410)* y(k,116) +rxt(k,409)* y(k,176) + het_rates(k,191)) &
                 * y(k,191)
         prod(k,75) = (.200_r8*rxt(k,399)*y(k,64) +.140_r8*rxt(k,411)*y(k,130) + &
                 rxt(k,414)*y(k,131))*y(k,190)
         loss(k,110) = (rxt(k,309)* y(k,116) +rxt(k,308)* y(k,176) + het_rates(k,192)) &
                 * y(k,192)
         prod(k,110) = (.500_r8*rxt(k,310)*y(k,132) +rxt(k,315)*y(k,27))*y(k,190)
         loss(k,141) = (rxt(k,313)* y(k,116) +rxt(k,311)* y(k,171) +rxt(k,312) &
                 * y(k,176) + het_rates(k,193))* y(k,193)
         prod(k,141) = (rxt(k,314)*y(k,133) +rxt(k,316)*y(k,46) + &
                 .150_r8*rxt(k,451)*y(k,151))*y(k,190) + (.060_r8*rxt(k,430)*y(k,4) + &
                 .060_r8*rxt(k,433)*y(k,105))*y(k,122) +.150_r8*rxt(k,69)*y(k,151)
         loss(k,139) = (rxt(k,442)* y(k,116) +rxt(k,440)* y(k,171) +rxt(k,441) &
                 * y(k,176) + het_rates(k,194))* y(k,194)
         prod(k,139) = (.500_r8*rxt(k,449)*y(k,118) +rxt(k,450)*y(k,190))*y(k,150) &
                  +rxt(k,443)*y(k,190)*y(k,147)
         loss(k,138) = (rxt(k,447)* y(k,116) +rxt(k,445)* y(k,171) +rxt(k,446) &
                 * y(k,176) + het_rates(k,195))* y(k,195)
         prod(k,138) = (rxt(k,431)*y(k,4) +rxt(k,434)*y(k,105) +rxt(k,448)*y(k,149)) &
                 *y(k,190)
         loss(k,106) = (rxt(k,417)* y(k,116) +rxt(k,416)* y(k,176) + het_rates(k,196)) &
                 * y(k,196)
         prod(k,106) = (rxt(k,418)*y(k,152) +.650_r8*rxt(k,419)*y(k,153))*y(k,190)
         loss(k,14) = (rxt(k,508)* y(k,116) +rxt(k,507)* y(k,176) + het_rates(k,197)) &
                 * y(k,197)
         prod(k,14) =rxt(k,506)*y(k,190)*y(k,153)
         loss(k,142) = (rxt(k,383)* y(k,116) +rxt(k,384)* y(k,118) +rxt(k,380) &
                 * y(k,170) +rxt(k,381)* y(k,171) +rxt(k,382)* y(k,176) &
                  + het_rates(k,198))* y(k,198)
         prod(k,142) = (rxt(k,352)*y(k,89) +rxt(k,353)*y(k,91) +rxt(k,354)*y(k,92) + &
                 .400_r8*rxt(k,377)*y(k,98) +.500_r8*rxt(k,385)*y(k,154))*y(k,190)
         loss(k,108) = (rxt(k,423)* y(k,116) +rxt(k,422)* y(k,176) + het_rates(k,199)) &
                 * y(k,199)
         prod(k,108) = (.560_r8*rxt(k,421)*y(k,155) +rxt(k,424)*y(k,156))*y(k,190)
         loss(k,15) = (rxt(k,511)* y(k,116) +rxt(k,510)* y(k,176) + het_rates(k,200)) &
                 * y(k,200)
         prod(k,15) =rxt(k,509)*y(k,190)*y(k,155)
         loss(k,81) = (rxt(k,426)* y(k,116) +rxt(k,425)* y(k,176) + het_rates(k,201)) &
                 * y(k,201)
         prod(k,81) = (.300_r8*rxt(k,427)*y(k,157) +rxt(k,428)*y(k,158))*y(k,190)
         loss(k,171) = (rxt(k,223)* y(k,68) +rxt(k,468)* y(k,138) +rxt(k,115) &
                 * y(k,189) + rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,202)) &
                 * y(k,202)
         prod(k,171) = (rxt(k,228)*y(k,39) +rxt(k,230)*y(k,41) +rxt(k,231)*y(k,42) + &
                 rxt(k,233)*y(k,44) +rxt(k,238)*y(k,79) +rxt(k,254)*y(k,52) + &
                 rxt(k,131)*y(k,72) +rxt(k,132)*y(k,74) +rxt(k,133)*y(k,176) + &
                 rxt(k,136)*y(k,190) +rxt(k,139)*y(k,85) +rxt(k,161)*y(k,84) + &
                 rxt(k,185)*y(k,80) +rxt(k,188)*y(k,87) +rxt(k,214)*y(k,76) + &
                 rxt(k,247)*y(k,40) +rxt(k,253)*y(k,51) +rxt(k,257)*y(k,82) + &
                 rxt(k,277)*y(k,26) +rxt(k,279)*y(k,43) +rxt(k,285)*y(k,48) + &
                 rxt(k,286)*y(k,49) +rxt(k,302)*y(k,28) +rxt(k,303)*y(k,29) + &
                 rxt(k,305)*y(k,47) +rxt(k,310)*y(k,132) +rxt(k,314)*y(k,133) + &
                 rxt(k,316)*y(k,46) +.500_r8*rxt(k,329)*y(k,100) +rxt(k,469)*y(k,112)) &
                 *y(k,190) + (rxt(k,513)*y(k,87) +rxt(k,519)*y(k,87) + &
                 rxt(k,520)*y(k,86) +rxt(k,524)*y(k,87) +rxt(k,525)*y(k,86))*y(k,80) &
                  + (rxt(k,471) +rxt(k,126)*y(k,71))*y(k,176) +rxt(k,109)*y(k,75)
      end do
      end subroutine imp_prod_loss
      end module mo_prod_loss
