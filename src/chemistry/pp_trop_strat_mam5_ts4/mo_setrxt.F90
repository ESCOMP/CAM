
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

      rate(:,85) = 1.2e-10_r8
      rate(:,91) = 1.2e-10_r8
      rate(:,92) = 1.2e-10_r8
      rate(:,93) = 1e-20_r8
      rate(:,94) = 1.3e-16_r8
      rate(:,96) = 4.2e-13_r8
      rate(:,98) = 8e-14_r8
      rate(:,99) = 3.9e-17_r8
      rate(:,106) = 6.9e-12_r8
      rate(:,107) = 7.2e-11_r8
      rate(:,108) = 1.6e-12_r8
      rate(:,114) = 1.8e-12_r8
      rate(:,118) = 1.8e-12_r8
      rate(:,130) = 3.5e-12_r8
      rate(:,132) = 1.3e-11_r8
      rate(:,133) = 2.2e-11_r8
      rate(:,134) = 5e-11_r8
      rate(:,169) = 1.7e-13_r8
      rate(:,171) = 2.607e-10_r8
      rate(:,172) = 9.75e-11_r8
      rate(:,173) = 2.07e-10_r8
      rate(:,174) = 2.088e-10_r8
      rate(:,175) = 1.17e-10_r8
      rate(:,176) = 4.644e-11_r8
      rate(:,177) = 1.204e-10_r8
      rate(:,178) = 9.9e-11_r8
      rate(:,179) = 3.3e-12_r8
      rate(:,198) = 4.5e-11_r8
      rate(:,199) = 4.62e-10_r8
      rate(:,200) = 1.2e-10_r8
      rate(:,201) = 9e-11_r8
      rate(:,202) = 3e-11_r8
      rate(:,215) = 2.57e-10_r8
      rate(:,216) = 1.8e-10_r8
      rate(:,217) = 1.794e-10_r8
      rate(:,218) = 1.3e-10_r8
      rate(:,219) = 7.65e-11_r8
      rate(:,230) = 1.31e-10_r8
      rate(:,231) = 3.5e-11_r8
      rate(:,232) = 9e-12_r8
      rate(:,236) = 6.8e-14_r8
      rate(:,237) = 2e-13_r8
      rate(:,251) = 1e-12_r8
      rate(:,255) = 1e-14_r8
      rate(:,256) = 1e-11_r8
      rate(:,257) = 1.15e-11_r8
      rate(:,258) = 4e-14_r8
      rate(:,271) = 3e-12_r8
      rate(:,272) = 6.7e-13_r8
      rate(:,282) = 1.4e-11_r8
      rate(:,285) = 2.4e-12_r8
      rate(:,296) = 5e-12_r8
      rate(:,302) = 3.5e-12_r8
      rate(:,307) = 2.4e-12_r8
      rate(:,308) = 1.4e-11_r8
      rate(:,312) = 2.4e-12_r8
      rate(:,317) = 4.5e-11_r8
      rate(:,322) = 2.4e-12_r8
      rate(:,331) = 2.3e-12_r8
      rate(:,333) = 1.2e-11_r8
      rate(:,334) = 5.7e-11_r8
      rate(:,335) = 2.8e-11_r8
      rate(:,336) = 6.6e-11_r8
      rate(:,337) = 1.4e-11_r8
      rate(:,340) = 1.9e-12_r8
      rate(:,347) = 6.34e-08_r8
      rate(:,351) = 1.157e-05_r8
      rate(:,369) = 1.29e-07_r8
      rate(:,370) = 2.31e-07_r8
      rate(:,371) = 2.31e-06_r8
      rate(:,372) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,86) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,87) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,88) = 2.64e-11_r8 * exp_fac(:)
      rate(:,89) = 3.3e-11_r8 * exp_fac(:)
      rate(:,90) = 6.6e-12_r8 * exp_fac(:)
      rate(:,95) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,97) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,100) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,101) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,104) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,105) = 1.4e-12_r8 * exp_fac(:)
      rate(:,313) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,110) = 3e-11_r8 * exp_fac(:)
      rate(:,196) = 5.5e-12_r8 * exp_fac(:)
      rate(:,228) = 3.8e-12_r8 * exp_fac(:)
      rate(:,241) = 3.8e-12_r8 * exp_fac(:)
      rate(:,267) = 3.8e-12_r8 * exp_fac(:)
      rate(:,275) = 3.8e-12_r8 * exp_fac(:)
      rate(:,279) = 3.8e-12_r8 * exp_fac(:)
      rate(:,290) = 2.3e-11_r8 * exp_fac(:)
      rate(:,315) = 1.52e-11_r8 * exp_fac(:)
      rate(:,323) = 1.52e-12_r8 * exp_fac(:)
      rate(:,111) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,112) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,113) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,115) = 4.8e-11_r8 * exp_fac(:)
      rate(:,194) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,116) = 1.8e-11_r8 * exp_fac(:)
      rate(:,253) = 4.2e-12_r8 * exp_fac(:)
      rate(:,274) = 4.2e-12_r8 * exp_fac(:)
      rate(:,311) = 4.4e-12_r8 * exp_fac(:)
      rate(:,117) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,121) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,122) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,123) = 2.9e-12_r8 * exp_fac(:)
      rate(:,124) = 1.45e-12_r8 * exp_fac(:)
      rate(:,125) = 1.45e-12_r8 * exp_fac(:)
      rate(:,126) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,127) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,128) = 1.2e-13_r8 * exp_fac(:)
      rate(:,154) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,131) = 1.7e-11_r8 * exp_fac(:)
      rate(:,222) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,135) = 3.44e-12_r8 * exp_fac(:)
      rate(:,187) = 2.3e-12_r8 * exp_fac(:)
      rate(:,190) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,136) = 3e-12_r8 * exp_fac(:)
      rate(:,195) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,138) = 7.26e-11_r8 * exp_fac(:)
      rate(:,139) = 4.64e-11_r8 * exp_fac(:)
      rate(:,146) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,147) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,148) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,149) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,150) = 1.4e-11_r8 * exp_fac(:)
      rate(:,164) = 7.4e-12_r8 * exp_fac(:)
      rate(:,249) = 8.1e-12_r8 * exp_fac(:)
      rate(:,151) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,152) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,153) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,155) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,156) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,157) = 2.6e-12_r8 * exp_fac(:)
      rate(:,158) = 6.4e-12_r8 * exp_fac(:)
      rate(:,188) = 4.1e-13_r8 * exp_fac(:)
      rate(:,159) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,161) = 3.6e-12_r8 * exp_fac(:)
      rate(:,204) = 2e-12_r8 * exp_fac(:)
      rate(:,162) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,163) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,165) = 6e-13_r8 * exp_fac(:)
      rate(:,185) = 1.5e-12_r8 * exp_fac(:)
      rate(:,193) = 1.9e-11_r8 * exp_fac(:)
      rate(:,166) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,167) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,168) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      rate(:,170) = 3e-12_r8 * exp( -500._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,182) = 1.7e-11_r8 * exp_fac(:)
      rate(:,203) = 6.3e-12_r8 * exp_fac(:)
      rate(:,183) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,184) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,186) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,189) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,192) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,197) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,205) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,206) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,207) = 1.64e-12_r8 * exp_fac(:)
      rate(:,298) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,208) = 2.03e-11_r8 * exp_fac(:)
      rate(:,339) = 3.4e-12_r8 * exp_fac(:)
      rate(:,209) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,210) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,211) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,212) = 1.25e-12_r8 * exp_fac(:)
      rate(:,221) = 3.4e-11_r8 * exp_fac(:)
      rate(:,213) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,214) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,220) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,223) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,224) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,225) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,226) = 2.8e-12_r8 * exp_fac(:)
      rate(:,278) = 2.9e-12_r8 * exp_fac(:)
      rate(:,227) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,229) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,235) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,238) = 7.5e-13_r8 * exp_fac(:)
      rate(:,252) = 7.5e-13_r8 * exp_fac(:)
      rate(:,273) = 7.5e-13_r8 * exp_fac(:)
      rate(:,277) = 8.6e-13_r8 * exp_fac(:)
      rate(:,284) = 8e-13_r8 * exp_fac(:)
      rate(:,305) = 8e-13_r8 * exp_fac(:)
      rate(:,310) = 8e-13_r8 * exp_fac(:)
      rate(:,320) = 8e-13_r8 * exp_fac(:)
      rate(:,239) = 2.6e-12_r8 * exp( 365._r8 * itemp(:) )
      rate(:,240) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,242) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,243) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,244) = 1.4e-12_r8 * exp_fac(:)
      rate(:,263) = 6.5e-15_r8 * exp_fac(:)
      rate(:,245) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,246) = 2.9e-12_r8 * exp_fac(:)
      rate(:,247) = 2e-12_r8 * exp_fac(:)
      rate(:,276) = 7.1e-13_r8 * exp_fac(:)
      rate(:,292) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,248) = 4.3e-13_r8 * exp_fac(:)
      rate(:,293) = 4.3e-13_r8 * exp_fac(:)
      rate(:,250) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,254) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,262) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,264) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      rate(:,265) = 1.41e-13_r8 * exp( 1300._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,266) = 2.7e-12_r8 * exp_fac(:)
      rate(:,286) = 2.7e-12_r8 * exp_fac(:)
      rate(:,287) = 1.3e-13_r8 * exp_fac(:)
      rate(:,289) = 9.6e-12_r8 * exp_fac(:)
      rate(:,295) = 5.3e-12_r8 * exp_fac(:)
      rate(:,306) = 2.7e-12_r8 * exp_fac(:)
      rate(:,321) = 2.7e-12_r8 * exp_fac(:)
      rate(:,268) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      exp_fac(:) = exp( -1860._r8 * itemp(:) )
      rate(:,269) = 1.4e-12_r8 * exp_fac(:)
      rate(:,316) = 1.4e-12_r8 * exp_fac(:)
      rate(:,270) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,283) = 5e-13_r8 * exp_fac(:)
      rate(:,309) = 5e-13_r8 * exp_fac(:)
      rate(:,319) = 5e-13_r8 * exp_fac(:)
      rate(:,288) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,291) = 4.6e-12_r8 * exp_fac(:)
      rate(:,294) = 2.3e-12_r8 * exp_fac(:)
      rate(:,299) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,303) = 1.86e-11_r8 * exp( 175._r8 * itemp(:) )
      rate(:,304) = 3.03e-12_r8 * exp( -446._r8 * itemp(:) )
      rate(:,314) = 2.54e-11_r8 * exp( 410._r8 * itemp(:) )
      rate(:,318) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      rate(:,324) = 1.2e-12_r8 * exp( 490._r8 * itemp(:) )
      rate(:,325) = 6.3e-16_r8 * exp( -580._r8 * itemp(:) )
      rate(:,326) = 1.2e-11_r8 * exp( 440._r8 * itemp(:) )
      rate(:,327) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,328) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,329) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,330) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,338) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,341) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,344) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,109), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,119), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,129), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,137), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,140), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,141), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,142), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,160), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,180), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,191), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,234), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,259), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,260), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,280), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,297), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,300), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,332), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,93) = 1e-20_r8
      rate(:n,94) = 1.3e-16_r8
      rate(:n,98) = 8e-14_r8
      rate(:n,99) = 3.9e-17_r8
      rate(:n,106) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,87) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,88) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,90) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,95) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,97) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,100) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,101) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,110) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,111) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,112) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,115) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,116) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,117) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,122) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,126) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,127) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,135) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,136) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,109) = wrk(:)

















      end subroutine setrxt_hrates

      end module mo_setrxt
