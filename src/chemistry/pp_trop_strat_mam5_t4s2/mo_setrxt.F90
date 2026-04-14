
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

      rate(:,85) = 0.000258_r8
      rate(:,86) = 0.085_r8
      rate(:,87) = 1.2e-10_r8
      rate(:,92) = 1.2e-10_r8
      rate(:,93) = 1.2e-10_r8
      rate(:,94) = 1e-20_r8
      rate(:,95) = 1.3e-16_r8
      rate(:,97) = 4.2e-13_r8
      rate(:,99) = 8e-14_r8
      rate(:,100) = 3.9e-17_r8
      rate(:,107) = 6.9e-12_r8
      rate(:,108) = 7.2e-11_r8
      rate(:,109) = 1.6e-12_r8
      rate(:,115) = 1.8e-12_r8
      rate(:,119) = 1.8e-12_r8
      rate(:,131) = 3.5e-12_r8
      rate(:,133) = 1.3e-11_r8
      rate(:,134) = 2.2e-11_r8
      rate(:,135) = 5e-11_r8
      rate(:,170) = 1.7e-13_r8
      rate(:,172) = 2.607e-10_r8
      rate(:,173) = 9.75e-11_r8
      rate(:,174) = 2.07e-10_r8
      rate(:,175) = 2.088e-10_r8
      rate(:,176) = 1.17e-10_r8
      rate(:,177) = 4.644e-11_r8
      rate(:,178) = 1.204e-10_r8
      rate(:,179) = 9.9e-11_r8
      rate(:,180) = 3.3e-12_r8
      rate(:,199) = 4.5e-11_r8
      rate(:,200) = 4.62e-10_r8
      rate(:,201) = 1.2e-10_r8
      rate(:,202) = 9e-11_r8
      rate(:,203) = 3e-11_r8
      rate(:,216) = 2.57e-10_r8
      rate(:,217) = 1.8e-10_r8
      rate(:,218) = 1.794e-10_r8
      rate(:,219) = 1.3e-10_r8
      rate(:,220) = 7.65e-11_r8
      rate(:,231) = 1.31e-10_r8
      rate(:,232) = 3.5e-11_r8
      rate(:,233) = 9e-12_r8
      rate(:,237) = 6.8e-14_r8
      rate(:,238) = 2e-13_r8
      rate(:,252) = 1e-12_r8
      rate(:,256) = 1e-14_r8
      rate(:,257) = 1e-11_r8
      rate(:,258) = 1.15e-11_r8
      rate(:,259) = 4e-14_r8
      rate(:,272) = 3e-12_r8
      rate(:,273) = 6.7e-13_r8
      rate(:,283) = 1.4e-11_r8
      rate(:,286) = 2.4e-12_r8
      rate(:,297) = 5e-12_r8
      rate(:,303) = 3.5e-12_r8
      rate(:,308) = 2.4e-12_r8
      rate(:,309) = 1.4e-11_r8
      rate(:,313) = 2.4e-12_r8
      rate(:,318) = 4.5e-11_r8
      rate(:,323) = 2.4e-12_r8
      rate(:,332) = 2.3e-12_r8
      rate(:,334) = 1.2e-11_r8
      rate(:,335) = 5.7e-11_r8
      rate(:,336) = 2.8e-11_r8
      rate(:,337) = 6.6e-11_r8
      rate(:,338) = 1.4e-11_r8
      rate(:,341) = 1.9e-12_r8
      rate(:,348) = 6.34e-08_r8
      rate(:,352) = 1.157e-05_r8
      rate(:,370) = 1.29e-07_r8
      rate(:,371) = 2.31e-07_r8
      rate(:,372) = 2.31e-06_r8
      rate(:,373) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,88) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,89) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,90) = 2.64e-11_r8 * exp_fac(:)
      rate(:,91) = 6.6e-12_r8 * exp_fac(:)
      rate(:,96) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,98) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,101) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,102) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,105) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,106) = 1.4e-12_r8 * exp_fac(:)
      rate(:,314) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,111) = 3e-11_r8 * exp_fac(:)
      rate(:,197) = 5.5e-12_r8 * exp_fac(:)
      rate(:,229) = 3.8e-12_r8 * exp_fac(:)
      rate(:,242) = 3.8e-12_r8 * exp_fac(:)
      rate(:,268) = 3.8e-12_r8 * exp_fac(:)
      rate(:,276) = 3.8e-12_r8 * exp_fac(:)
      rate(:,280) = 3.8e-12_r8 * exp_fac(:)
      rate(:,291) = 2.3e-11_r8 * exp_fac(:)
      rate(:,316) = 1.52e-11_r8 * exp_fac(:)
      rate(:,324) = 1.52e-12_r8 * exp_fac(:)
      rate(:,112) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,113) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,114) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,116) = 4.8e-11_r8 * exp_fac(:)
      rate(:,195) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,117) = 1.8e-11_r8 * exp_fac(:)
      rate(:,254) = 4.2e-12_r8 * exp_fac(:)
      rate(:,275) = 4.2e-12_r8 * exp_fac(:)
      rate(:,312) = 4.4e-12_r8 * exp_fac(:)
      rate(:,118) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,122) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,123) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,124) = 2.9e-12_r8 * exp_fac(:)
      rate(:,125) = 1.45e-12_r8 * exp_fac(:)
      rate(:,126) = 1.45e-12_r8 * exp_fac(:)
      rate(:,127) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,128) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,129) = 1.2e-13_r8 * exp_fac(:)
      rate(:,155) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,132) = 1.7e-11_r8 * exp_fac(:)
      rate(:,223) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,136) = 3.44e-12_r8 * exp_fac(:)
      rate(:,188) = 2.3e-12_r8 * exp_fac(:)
      rate(:,191) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,137) = 3e-12_r8 * exp_fac(:)
      rate(:,196) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,139) = 7.26e-11_r8 * exp_fac(:)
      rate(:,140) = 4.64e-11_r8 * exp_fac(:)
      rate(:,147) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,148) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,149) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,150) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,151) = 1.4e-11_r8 * exp_fac(:)
      rate(:,165) = 7.4e-12_r8 * exp_fac(:)
      rate(:,250) = 8.1e-12_r8 * exp_fac(:)
      rate(:,152) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,153) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,154) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,156) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,157) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,158) = 2.6e-12_r8 * exp_fac(:)
      rate(:,159) = 6.4e-12_r8 * exp_fac(:)
      rate(:,189) = 4.1e-13_r8 * exp_fac(:)
      rate(:,160) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,162) = 3.6e-12_r8 * exp_fac(:)
      rate(:,205) = 2e-12_r8 * exp_fac(:)
      rate(:,163) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,164) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,166) = 6e-13_r8 * exp_fac(:)
      rate(:,186) = 1.5e-12_r8 * exp_fac(:)
      rate(:,194) = 1.9e-11_r8 * exp_fac(:)
      rate(:,167) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,168) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,169) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      rate(:,171) = 3e-12_r8 * exp( -500._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,183) = 1.7e-11_r8 * exp_fac(:)
      rate(:,204) = 6.3e-12_r8 * exp_fac(:)
      rate(:,184) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,185) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,187) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,190) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,193) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,198) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,206) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,207) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,208) = 1.64e-12_r8 * exp_fac(:)
      rate(:,299) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,209) = 2.03e-11_r8 * exp_fac(:)
      rate(:,340) = 3.4e-12_r8 * exp_fac(:)
      rate(:,210) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,211) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,212) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,213) = 1.25e-12_r8 * exp_fac(:)
      rate(:,222) = 3.4e-11_r8 * exp_fac(:)
      rate(:,214) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,215) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,221) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,224) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,225) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,226) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,227) = 2.8e-12_r8 * exp_fac(:)
      rate(:,279) = 2.9e-12_r8 * exp_fac(:)
      rate(:,228) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,230) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,236) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,239) = 7.5e-13_r8 * exp_fac(:)
      rate(:,253) = 7.5e-13_r8 * exp_fac(:)
      rate(:,274) = 7.5e-13_r8 * exp_fac(:)
      rate(:,278) = 8.6e-13_r8 * exp_fac(:)
      rate(:,285) = 8e-13_r8 * exp_fac(:)
      rate(:,306) = 8e-13_r8 * exp_fac(:)
      rate(:,311) = 8e-13_r8 * exp_fac(:)
      rate(:,321) = 8e-13_r8 * exp_fac(:)
      rate(:,240) = 2.6e-12_r8 * exp( 365._r8 * itemp(:) )
      rate(:,241) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,243) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,244) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,245) = 1.4e-12_r8 * exp_fac(:)
      rate(:,264) = 6.5e-15_r8 * exp_fac(:)
      rate(:,246) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,247) = 2.9e-12_r8 * exp_fac(:)
      rate(:,248) = 2e-12_r8 * exp_fac(:)
      rate(:,277) = 7.1e-13_r8 * exp_fac(:)
      rate(:,293) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,249) = 4.3e-13_r8 * exp_fac(:)
      rate(:,294) = 4.3e-13_r8 * exp_fac(:)
      rate(:,251) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,255) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,263) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,265) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      rate(:,266) = 1.41e-13_r8 * exp( 1300._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,267) = 2.7e-12_r8 * exp_fac(:)
      rate(:,287) = 2.7e-12_r8 * exp_fac(:)
      rate(:,288) = 1.3e-13_r8 * exp_fac(:)
      rate(:,290) = 9.6e-12_r8 * exp_fac(:)
      rate(:,296) = 5.3e-12_r8 * exp_fac(:)
      rate(:,307) = 2.7e-12_r8 * exp_fac(:)
      rate(:,322) = 2.7e-12_r8 * exp_fac(:)
      rate(:,269) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      exp_fac(:) = exp( -1860._r8 * itemp(:) )
      rate(:,270) = 1.4e-12_r8 * exp_fac(:)
      rate(:,317) = 1.4e-12_r8 * exp_fac(:)
      rate(:,271) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,284) = 5e-13_r8 * exp_fac(:)
      rate(:,310) = 5e-13_r8 * exp_fac(:)
      rate(:,320) = 5e-13_r8 * exp_fac(:)
      rate(:,289) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,292) = 4.6e-12_r8 * exp_fac(:)
      rate(:,295) = 2.3e-12_r8 * exp_fac(:)
      rate(:,300) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,304) = 1.86e-11_r8 * exp( 175._r8 * itemp(:) )
      rate(:,305) = 3.03e-12_r8 * exp( -446._r8 * itemp(:) )
      rate(:,315) = 2.54e-11_r8 * exp( 410._r8 * itemp(:) )
      rate(:,319) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      rate(:,325) = 1.2e-12_r8 * exp( 490._r8 * itemp(:) )
      rate(:,326) = 6.3e-16_r8 * exp( -580._r8 * itemp(:) )
      rate(:,327) = 1.2e-11_r8 * exp( 440._r8 * itemp(:) )
      rate(:,328) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,329) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,330) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,331) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,339) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,342) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,345) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,110), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,120), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,130), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,138), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,141), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,142), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,143), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,161), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,181), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,192), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,235), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,260), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,261), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,281), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,298), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,301), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,333), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,94) = 1e-20_r8
      rate(:n,95) = 1.3e-16_r8
      rate(:n,99) = 8e-14_r8
      rate(:n,100) = 3.9e-17_r8
      rate(:n,107) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,89) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,90) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,91) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,96) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,98) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,101) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,102) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,111) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,112) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,113) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,116) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,117) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,118) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,123) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,127) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,128) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,136) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,137) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,110) = wrk(:)

















      end subroutine setrxt_hrates

      end module mo_setrxt
