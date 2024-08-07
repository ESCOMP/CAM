
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      real(r8)  :: itemp(ncol*pver)
      real(r8)  :: exp_fac(ncol*pver)
      real(r8)  :: ko(ncol*pver)
      real(r8)  :: kinf(ncol*pver)

      rate(:,86) = 1.2e-10_r8
      rate(:,90) = 1.2e-10_r8
      rate(:,91) = 1.2e-10_r8
      rate(:,97) = 6.9e-12_r8
      rate(:,98) = 7.2e-11_r8
      rate(:,99) = 1.6e-12_r8
      rate(:,105) = 1.8e-12_r8
      rate(:,109) = 1.8e-12_r8
      rate(:,121) = 3.5e-12_r8
      rate(:,123) = 1.3e-11_r8
      rate(:,124) = 2.2e-11_r8
      rate(:,125) = 5e-11_r8
      rate(:,160) = 1.7e-13_r8
      rate(:,162) = 2.607e-10_r8
      rate(:,163) = 9.75e-11_r8
      rate(:,164) = 2.07e-10_r8
      rate(:,165) = 2.088e-10_r8
      rate(:,166) = 1.17e-10_r8
      rate(:,167) = 4.644e-11_r8
      rate(:,168) = 1.204e-10_r8
      rate(:,169) = 9.9e-11_r8
      rate(:,170) = 3.3e-12_r8
      rate(:,189) = 4.5e-11_r8
      rate(:,190) = 4.62e-10_r8
      rate(:,191) = 1.2e-10_r8
      rate(:,192) = 9e-11_r8
      rate(:,193) = 3e-11_r8
      rate(:,206) = 2.57e-10_r8
      rate(:,207) = 1.8e-10_r8
      rate(:,208) = 1.794e-10_r8
      rate(:,209) = 1.3e-10_r8
      rate(:,210) = 7.65e-11_r8
      rate(:,221) = 1.31e-10_r8
      rate(:,222) = 3.5e-11_r8
      rate(:,223) = 9e-12_r8
      rate(:,227) = 6.8e-14_r8
      rate(:,228) = 2e-13_r8
      rate(:,242) = 1e-12_r8
      rate(:,246) = 1e-14_r8
      rate(:,247) = 1e-11_r8
      rate(:,248) = 1.15e-11_r8
      rate(:,249) = 4e-14_r8
      rate(:,262) = 3e-12_r8
      rate(:,263) = 6.7e-13_r8
      rate(:,273) = 1.4e-11_r8
      rate(:,276) = 2.4e-12_r8
      rate(:,287) = 5e-12_r8
      rate(:,293) = 3.5e-12_r8
      rate(:,298) = 2.4e-12_r8
      rate(:,299) = 1.4e-11_r8
      rate(:,303) = 2.4e-12_r8
      rate(:,308) = 4.5e-11_r8
      rate(:,313) = 2.4e-12_r8
      rate(:,322) = 2.3e-12_r8
      rate(:,324) = 1.2e-11_r8
      rate(:,325) = 5.7e-11_r8
      rate(:,326) = 2.8e-11_r8
      rate(:,327) = 6.6e-11_r8
      rate(:,328) = 1.4e-11_r8
      rate(:,331) = 1.9e-12_r8
      rate(:,338) = 6.34e-08_r8
      rate(:,342) = 1.157e-05_r8
      rate(:,360) = 1.29e-07_r8
      rate(:,361) = 2.31e-07_r8
      rate(:,362) = 2.31e-06_r8
      rate(:,363) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,87) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,88) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,89) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,92) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,95) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,96) = 1.4e-12_r8 * exp_fac(:)
      rate(:,304) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,101) = 3e-11_r8 * exp_fac(:)
      rate(:,187) = 5.5e-12_r8 * exp_fac(:)
      rate(:,219) = 3.8e-12_r8 * exp_fac(:)
      rate(:,232) = 3.8e-12_r8 * exp_fac(:)
      rate(:,258) = 3.8e-12_r8 * exp_fac(:)
      rate(:,266) = 3.8e-12_r8 * exp_fac(:)
      rate(:,270) = 3.8e-12_r8 * exp_fac(:)
      rate(:,281) = 2.3e-11_r8 * exp_fac(:)
      rate(:,306) = 1.52e-11_r8 * exp_fac(:)
      rate(:,314) = 1.52e-12_r8 * exp_fac(:)
      rate(:,102) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,103) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,104) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,106) = 4.8e-11_r8 * exp_fac(:)
      rate(:,185) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,107) = 1.8e-11_r8 * exp_fac(:)
      rate(:,244) = 4.2e-12_r8 * exp_fac(:)
      rate(:,257) = 4.2e-12_r8 * exp_fac(:)
      rate(:,265) = 4.2e-12_r8 * exp_fac(:)
      rate(:,302) = 4.4e-12_r8 * exp_fac(:)
      rate(:,108) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,112) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,113) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,114) = 2.9e-12_r8 * exp_fac(:)
      rate(:,115) = 1.45e-12_r8 * exp_fac(:)
      rate(:,116) = 1.45e-12_r8 * exp_fac(:)
      rate(:,117) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,118) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,119) = 1.2e-13_r8 * exp_fac(:)
      rate(:,145) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,122) = 1.7e-11_r8 * exp_fac(:)
      rate(:,213) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,126) = 3.44e-12_r8 * exp_fac(:)
      rate(:,178) = 2.3e-12_r8 * exp_fac(:)
      rate(:,181) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,127) = 3e-12_r8 * exp_fac(:)
      rate(:,186) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,129) = 7.26e-11_r8 * exp_fac(:)
      rate(:,130) = 4.64e-11_r8 * exp_fac(:)
      rate(:,137) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,138) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,139) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,140) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,141) = 1.4e-11_r8 * exp_fac(:)
      rate(:,155) = 7.4e-12_r8 * exp_fac(:)
      rate(:,240) = 8.1e-12_r8 * exp_fac(:)
      rate(:,142) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,143) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,144) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,146) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,147) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,148) = 2.6e-12_r8 * exp_fac(:)
      rate(:,149) = 6.4e-12_r8 * exp_fac(:)
      rate(:,179) = 4.1e-13_r8 * exp_fac(:)
      rate(:,150) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,152) = 3.6e-12_r8 * exp_fac(:)
      rate(:,195) = 2e-12_r8 * exp_fac(:)
      rate(:,153) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,154) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,156) = 6e-13_r8 * exp_fac(:)
      rate(:,176) = 1.5e-12_r8 * exp_fac(:)
      rate(:,184) = 1.9e-11_r8 * exp_fac(:)
      rate(:,157) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,158) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,159) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      rate(:,161) = 3e-12_r8 * exp( -500._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,173) = 1.7e-11_r8 * exp_fac(:)
      rate(:,194) = 6.3e-12_r8 * exp_fac(:)
      rate(:,174) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,175) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,177) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,180) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,183) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,188) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,196) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,197) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,198) = 1.64e-12_r8 * exp_fac(:)
      rate(:,289) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,199) = 2.03e-11_r8 * exp_fac(:)
      rate(:,330) = 3.4e-12_r8 * exp_fac(:)
      rate(:,200) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,201) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,202) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,203) = 1.25e-12_r8 * exp_fac(:)
      rate(:,212) = 3.4e-11_r8 * exp_fac(:)
      rate(:,204) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,205) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,211) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,214) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,215) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,216) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,217) = 2.8e-12_r8 * exp_fac(:)
      rate(:,269) = 2.9e-12_r8 * exp_fac(:)
      rate(:,218) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,220) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,226) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,229) = 7.5e-13_r8 * exp_fac(:)
      rate(:,243) = 7.5e-13_r8 * exp_fac(:)
      rate(:,256) = 7.5e-13_r8 * exp_fac(:)
      rate(:,264) = 7.5e-13_r8 * exp_fac(:)
      rate(:,268) = 8.6e-13_r8 * exp_fac(:)
      rate(:,275) = 8e-13_r8 * exp_fac(:)
      rate(:,296) = 8e-13_r8 * exp_fac(:)
      rate(:,301) = 8e-13_r8 * exp_fac(:)
      rate(:,311) = 8e-13_r8 * exp_fac(:)
      rate(:,230) = 2.6e-12_r8 * exp( 365._r8 * itemp(:) )
      rate(:,231) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,233) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,234) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,235) = 1.4e-12_r8 * exp_fac(:)
      rate(:,254) = 6.5e-15_r8 * exp_fac(:)
      rate(:,236) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,237) = 2.9e-12_r8 * exp_fac(:)
      rate(:,238) = 2e-12_r8 * exp_fac(:)
      rate(:,267) = 7.1e-13_r8 * exp_fac(:)
      rate(:,283) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,239) = 4.3e-13_r8 * exp_fac(:)
      rate(:,284) = 4.3e-13_r8 * exp_fac(:)
      rate(:,241) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,245) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,253) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,255) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,259) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      exp_fac(:) = exp( -1860._r8 * itemp(:) )
      rate(:,260) = 1.4e-12_r8 * exp_fac(:)
      rate(:,307) = 1.4e-12_r8 * exp_fac(:)
      rate(:,261) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,274) = 5e-13_r8 * exp_fac(:)
      rate(:,300) = 5e-13_r8 * exp_fac(:)
      rate(:,310) = 5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,277) = 2.7e-12_r8 * exp_fac(:)
      rate(:,278) = 1.3e-13_r8 * exp_fac(:)
      rate(:,280) = 9.6e-12_r8 * exp_fac(:)
      rate(:,286) = 5.3e-12_r8 * exp_fac(:)
      rate(:,297) = 2.7e-12_r8 * exp_fac(:)
      rate(:,312) = 2.7e-12_r8 * exp_fac(:)
      rate(:,279) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,282) = 4.6e-12_r8 * exp_fac(:)
      rate(:,285) = 2.3e-12_r8 * exp_fac(:)
      rate(:,290) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,294) = 1.86e-11_r8 * exp( 175._r8 * itemp(:) )
      rate(:,295) = 3.03e-12_r8 * exp( -446._r8 * itemp(:) )
      rate(:,305) = 2.54e-11_r8 * exp( 410._r8 * itemp(:) )
      rate(:,309) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      rate(:,315) = 1.2e-12_r8 * exp( 490._r8 * itemp(:) )
      rate(:,316) = 6.3e-16_r8 * exp( -580._r8 * itemp(:) )
      rate(:,317) = 1.2e-11_r8 * exp( 440._r8 * itemp(:) )
      rate(:,318) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,319) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,320) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,321) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,329) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,332) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,335) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,100), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,110), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,120), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,128), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,131), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,132), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,133), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,151), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,171), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,182), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,225), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,250), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,251), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,271), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,288), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,291), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,323), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      integer   ::  k
      real(r8)  :: itemp(ncol*kbot)
      real(r8)  :: exp_fac(ncol*kbot)
      real(r8)  :: ko(ncol*kbot)
      real(r8)  :: kinf(ncol*kbot)
      real(r8)  :: wrk(ncol*kbot)
 
      n = ncol*kbot

      rate(:n,97) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,88) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,92) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,101) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,102) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,103) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,106) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,107) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,108) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,113) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,117) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,118) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,126) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,127) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,100) = wrk(:)

















      end subroutine setrxt_hrates

      end module mo_setrxt
