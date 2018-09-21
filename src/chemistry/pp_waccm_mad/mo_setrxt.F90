
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

      rate(:,104) = 8.00e-14_r8
      rate(:,105) = 3.90e-17_r8
      rate(:,108) = 4.20e-13_r8
      rate(:,109) = 8.50e-2_r8
      rate(:,110) = 1.30e-16_r8
      rate(:,112) = 1.00e-20_r8
      rate(:,113) = 2.58e-04_r8
      rate(:,120) = 1.20e-10_r8
      rate(:,121) = 2.02e-10_r8
      rate(:,122) = 1.204e-10_r8
      rate(:,123) = 1.50e-10_r8
      rate(:,124) = 9.75e-11_r8
      rate(:,125) = 1.50e-11_r8
      rate(:,126) = 7.20e-11_r8
      rate(:,127) = 1.794e-10_r8
      rate(:,128) = 1.628e-10_r8
      rate(:,129) = 2.84e-10_r8
      rate(:,130) = 1.674e-10_r8
      rate(:,131) = 9.60e-11_r8
      rate(:,132) = 4.10e-11_r8
      rate(:,133) = 1.012e-10_r8
      rate(:,134) = 1.20e-10_r8
      rate(:,135) = 4.49e-10_r8
      rate(:,136) = 2.57e-10_r8
      rate(:,137) = 2.14e-11_r8
      rate(:,138) = 1.90e-10_r8
      rate(:,139) = 1.31e-10_r8
      rate(:,140) = 3.50e-11_r8
      rate(:,141) = 9.00e-12_r8
      rate(:,142) = 1.20e-10_r8
      rate(:,143) = 1.50e-10_r8
      rate(:,144) = 1.20e-10_r8
      rate(:,147) = 7.20e-11_r8
      rate(:,148) = 6.90e-12_r8
      rate(:,149) = 1.60e-12_r8
      rate(:,153) = 1.80e-12_r8
      rate(:,156) = 1.80e-12_r8
      rate(:,164) = 5.00e-12_r8
      rate(:,165) = 7.00e-13_r8
      rate(:,166) = 5.00e-11_r8
      rate(:,183) = 1.00e-11_r8
      rate(:,184) = 2.20e-11_r8
      rate(:,185) = 3.50e-12_r8
      rate(:,210) = 1.70e-13_r8
      rate(:,282) = 9.0e-10_r8
      rate(:,283) = 1.0e-10_r8
      rate(:,284) = 4.4e-10_r8
      rate(:,285) = 4.0e-10_r8
      rate(:,286) = 2.0e-10_r8
      rate(:,287) = 1.0e-12_r8
      rate(:,288) = 6.0e-11_r8
      rate(:,289) = 5.0e-16_r8
      rate(:,293) = 4.8e-10_r8
      rate(:,294) = 1.0e-10_r8
      rate(:,295) = 4.0e-10_r8
      rate(:,298) = 5.0e-12_r8
      rate(:,299) = 7.0e-10_r8
      rate(:,300) = 8.0e-10_r8
      rate(:,302) = 4.7e-2_r8
      rate(:,303) = 1.71e-1_r8
      rate(:,304) = 7.7e-5_r8
      rate(:,309) = 1.5e-6_r8
      rate(:,310) = 1.5e-6_r8
      rate(:,311) = 2.0e-6_r8
      rate(:,312) = 2.0e-6_r8
      rate(:,313) = 2.0e-6_r8
      rate(:,314) = 1.5e-6_r8
      rate(:,318) = 3.6e-6_r8
      rate(:,319) = 5.0e-6_r8
      rate(:,321) = 1e-9_r8
      rate(:,322) = 1e-9_r8
      rate(:,324) = 2.8e-28_r8
      rate(:,325) = 1.7e-9_r8
      rate(:,326) = 1.5e-10_r8
      rate(:,327) = 3.e-10_r8
      rate(:,328) = 9.0e-10_r8
      rate(:,329) = 2.4e-10_r8
      rate(:,330) = 2.0e-9_r8
      rate(:,339) = 4.0e-12_r8
      rate(:,340) = 7.0e-12_r8
      rate(:,341) = 1.0e-9_r8
      rate(:,342) = 1.0e-9_r8
      rate(:,345) = 0.5e-9_r8
      rate(:,346) = 1.e-10_r8
      rate(:,347) = 7.e-12_r8
      rate(:,349) = 7e-11_r8
      rate(:,350) = 1.e-9_r8
      rate(:,355) = 1.9e-10_r8
      rate(:,357) = 3.e-10_r8
      rate(:,358) = 0.5e-12_r8
      rate(:,359) = 5.8e-10_r8
      rate(:,360) = 1.5e-10_r8
      rate(:,361) = 2.e-10_r8
      rate(:,363) = 1.4e-9_r8
      rate(:,364) = 1.e-10_r8
      rate(:,365) = 1.e-10_r8
      rate(:,366) = 2.e-10_r8
      rate(:,367) = 1.4e-9_r8
      rate(:,368) = 8.0e-10_r8
      rate(:,369) = 2.9e-31_r8
      rate(:,370) = 6.0e-13_r8
      rate(:,371) = 1.0e-9_r8
      rate(:,372) = 2.0e-28_r8
      rate(:,373) = 3.2e-11_r8
      rate(:,374) = 3.6e-9_r8
      rate(:,375) = 2.e-9_r8
      rate(:,376) = 1.e-10_r8
      rate(:,377) = 1.e-10_r8
      rate(:,378) = 1.5e-10_r8
      rate(:,379) = 7.8e-10_r8
      rate(:,380) = 9.9e-30_r8
      rate(:,381) = 7.e-10_r8
      rate(:,382) = 3.4e-31_r8
      rate(:,383) = 2.9e-9_r8
      rate(:,384) = 1.6e-9_r8
      rate(:,385) = 1.e-10_r8
      rate(:,386) = 1.e-10_r8
      rate(:,387) = 2.5e-10_r8
      rate(:,388) = 8.4e-10_r8
      rate(:,389) = 5.5e-10_r8
      rate(:,394) = 4.e-10_r8
      rate(:,395) = 4.3e-10_r8
      rate(:,396) = 9.e-10_r8
      rate(:,397) = 1.1e-9_r8
      rate(:,398) = 7.6e-28_r8
      rate(:,399) = 1.e-9_r8
      rate(:,400) = 1.e-10_r8
      rate(:,401) = 1.e-10_r8
      rate(:,402) = 1.1e-10_r8
      rate(:,403) = 6.0e-15_r8
      rate(:,404) = 1.7e-10_r8
      rate(:,407) = 3.51e-10_r8
      rate(:,408) = 1.e-10_r8
      rate(:,409) = 1.e-10_r8
      rate(:,410) = 1.e-11_r8
      rate(:,411) = 1.3e-10_r8
      rate(:,412) = 2.2e-10_r8
      rate(:,413) = 1.4e-10_r8
      rate(:,414) = 1.2e-9_r8
      rate(:,415) = 1.e-10_r8
      rate(:,416) = 1.0e-10_r8
      rate(:,417) = 3.e-10_r8
      rate(:,418) = 2.e-13_r8
      rate(:,419) = 1.2e-10_r8
      rate(:,420) = 1.6e-9_r8
      rate(:,421) = 1.4e-9_r8
      rate(:,422) = 1.0e-10_r8
      rate(:,423) = 1.0e-10_r8
      rate(:,424) = 0.5e-11_r8
      rate(:,425) = 1.e-13_r8
      rate(:,426) = 1.e-12_r8
      rate(:,427) = 9.6e-10_r8
      rate(:,428) = 6.0e-12_r8
      rate(:,429) = 1.6e-9_r8
      rate(:,430) = 2.e-29_r8
      rate(:,431) = 1.0e-27_r8
      rate(:,432) = 2.9e-12_r8
      rate(:,433) = 2.9e-11_r8
      rate(:,434) = 2.0e-10_r8
      rate(:,435) = 1.3e-9_r8
      rate(:,438) = 1.0e-28_r8
      rate(:,439) = 1.6e-28_r8
      rate(:,441) = 3.5e-12_r8
      rate(:,442) = 4.0e-11_r8
      rate(:,444) = 4.0e-11_r8
      rate(:,445) = 1.0e-28_r8
      rate(:,447) = 3.5e-12_r8
      rate(:,449) = 1.6e-28_r8
      rate(:,450) = 1.6e-28_r8
      rate(:,452) = 7.0e-10_r8
      rate(:,454) = 1.45e-26_r8
      rate(:,455) = 7.6e-10_r8
      rate(:,457) = 7.0e-10_r8
      rate(:,458) = 1.6e-9_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,102) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,106) = 1.80e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,107) = 3.50e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,111) = 3.60e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,114) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,115) = 3.135e-11_r8 * exp_fac(:)
      rate(:,116) = 1.65e-12_r8 * exp_fac(:)
      rate(:,117) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,118) = 7.25e-11_r8 * exp_fac(:)
      rate(:,119) = 4.63e-11_r8 * exp_fac(:)
      rate(:,146) = 1.40e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,150) = 1.80e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:,151) = 1.70e-12_r8 * exp( -940._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,152) = 4.80e-11_r8 * exp_fac(:)
      rate(:,220) = 1.70e-11_r8 * exp_fac(:)
      rate(:,155) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:) )
      rate(:,157) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,158) = 3.00e-11_r8 * exp_fac(:)
      rate(:,228) = 5.50e-12_r8 * exp_fac(:)
      rate(:,256) = 3.80e-12_r8 * exp_fac(:)
      rate(:,159) = 1.00e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,161) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:) )
      rate(:,163) = 1.8e-11_r8 * exp( 390._r8 * itemp(:) )
      rate(:,167) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,168) = 2.10e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,169) = 2.90e-12_r8 * exp_fac(:)
      rate(:,170) = 1.45e-12_r8 * exp_fac(:)
      rate(:,171) = 1.45e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,173) = 3.30e-12_r8 * exp_fac(:)
      rate(:,192) = 1.40e-11_r8 * exp_fac(:)
      rate(:,197) = 7.40e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,174) = 3.00e-12_r8 * exp_fac(:)
      rate(:,229) = 5.80e-12_r8 * exp_fac(:)
      rate(:,175) = 5.10e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,177) = 1.20e-13_r8 * exp_fac(:)
      rate(:,203) = 3.00e-11_r8 * exp_fac(:)
      rate(:,182) = 1.50e-11_r8 * exp( 170._r8 * itemp(:) )
      rate(:,187) = 1.30e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,189) = 2.30e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,190) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,191) = 1.10e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,193) = 3.60e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,194) = 8.10e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,195) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:) )
      rate(:,196) = 2.80e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,198) = 6.00e-13_r8 * exp_fac(:)
      rate(:,219) = 1.90e-11_r8 * exp_fac(:)
      rate(:,227) = 1.50e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,199) = 2.60e-12_r8 * exp_fac(:)
      rate(:,201) = 6.40e-12_r8 * exp_fac(:)
      rate(:,226) = 4.10e-13_r8 * exp_fac(:)
      rate(:,200) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,204) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,205) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:) )
      rate(:,208) = 1.80e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,209) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,211) = 3.40e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,212) = 3.00e-12_r8 * exp_fac(:)
      rate(:,233) = 1.40e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,213) = 3.60e-12_r8 * exp_fac(:)
      rate(:,244) = 2.00e-12_r8 * exp_fac(:)
      rate(:,214) = 1.20e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,215) = 6.50e-12_r8 * exp( 135._r8 * itemp(:) )
      rate(:,216) = 1.60e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,217) = 4.80e-12_r8 * exp( -310._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,218) = 1.70e-11_r8 * exp_fac(:)
      rate(:,246) = 6.30e-12_r8 * exp_fac(:)
      rate(:,221) = 4.50e-12_r8 * exp( 460._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,222) = 8.80e-12_r8 * exp_fac(:)
      rate(:,225) = 2.30e-12_r8 * exp_fac(:)
      rate(:,224) = 9.50e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,230) = 1.20e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,231) = 1.90e-11_r8 * exp( 215._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,232) = 1.40e-11_r8 * exp_fac(:)
      rate(:,282) = 9.0e-10_r8 * exp_fac(:)
      rate(:,283) = 1.0e-10_r8 * exp_fac(:)
      rate(:,284) = 4.4e-10_r8 * exp_fac(:)
      rate(:,285) = 4.0e-10_r8 * exp_fac(:)
      rate(:,286) = 2.0e-10_r8 * exp_fac(:)
      rate(:,287) = 1.0e-12_r8 * exp_fac(:)
      rate(:,288) = 6.0e-11_r8 * exp_fac(:)
      rate(:,289) = 5.0e-16_r8 * exp_fac(:)
      rate(:,293) = 4.8e-10_r8 * exp_fac(:)
      rate(:,294) = 1.0e-10_r8 * exp_fac(:)
      rate(:,295) = 4.0e-10_r8 * exp_fac(:)
      rate(:,298) = 5.0e-12_r8 * exp_fac(:)
      rate(:,299) = 7.0e-10_r8 * exp_fac(:)
      rate(:,300) = 8.0e-10_r8 * exp_fac(:)
      rate(:,302) = 4.7e-2_r8 * exp_fac(:)
      rate(:,303) = 1.71e-1_r8 * exp_fac(:)
      rate(:,304) = 7.7e-5_r8 * exp_fac(:)
      rate(:,309) = 1.5e-6_r8 * exp_fac(:)
      rate(:,310) = 1.5e-6_r8 * exp_fac(:)
      rate(:,311) = 2.0e-6_r8 * exp_fac(:)
      rate(:,312) = 2.0e-6_r8 * exp_fac(:)
      rate(:,313) = 2.0e-6_r8 * exp_fac(:)
      rate(:,314) = 1.5e-6_r8 * exp_fac(:)
      rate(:,318) = 3.6e-6_r8 * exp_fac(:)
      rate(:,319) = 5.0e-6_r8 * exp_fac(:)
      rate(:,321) = 1e-9_r8 * exp_fac(:)
      rate(:,322) = 1e-9_r8 * exp_fac(:)
      rate(:,324) = 2.8e-28_r8 * exp_fac(:)
      rate(:,325) = 1.7e-9_r8 * exp_fac(:)
      rate(:,326) = 1.5e-10_r8 * exp_fac(:)
      rate(:,327) = 3.e-10_r8 * exp_fac(:)
      rate(:,328) = 9.0e-10_r8 * exp_fac(:)
      rate(:,329) = 2.4e-10_r8 * exp_fac(:)
      rate(:,330) = 2.0e-9_r8 * exp_fac(:)
      rate(:,339) = 4.0e-12_r8 * exp_fac(:)
      rate(:,340) = 7.0e-12_r8 * exp_fac(:)
      rate(:,341) = 1.0e-9_r8 * exp_fac(:)
      rate(:,342) = 1.0e-9_r8 * exp_fac(:)
      rate(:,345) = 0.5e-9_r8 * exp_fac(:)
      rate(:,346) = 1.e-10_r8 * exp_fac(:)
      rate(:,347) = 7.e-12_r8 * exp_fac(:)
      rate(:,349) = 7e-11_r8 * exp_fac(:)
      rate(:,350) = 1.e-9_r8 * exp_fac(:)
      rate(:,355) = 1.9e-10_r8 * exp_fac(:)
      rate(:,357) = 3.e-10_r8 * exp_fac(:)
      rate(:,358) = 0.5e-12_r8 * exp_fac(:)
      rate(:,359) = 5.8e-10_r8 * exp_fac(:)
      rate(:,360) = 1.5e-10_r8 * exp_fac(:)
      rate(:,361) = 2.e-10_r8 * exp_fac(:)
      rate(:,363) = 1.4e-9_r8 * exp_fac(:)
      rate(:,364) = 1.e-10_r8 * exp_fac(:)
      rate(:,365) = 1.e-10_r8 * exp_fac(:)
      rate(:,366) = 2.e-10_r8 * exp_fac(:)
      rate(:,367) = 1.4e-9_r8 * exp_fac(:)
      rate(:,368) = 8.0e-10_r8 * exp_fac(:)
      rate(:,369) = 2.9e-31_r8 * exp_fac(:)
      rate(:,370) = 6.0e-13_r8 * exp_fac(:)
      rate(:,371) = 1.0e-9_r8 * exp_fac(:)
      rate(:,372) = 2.0e-28_r8 * exp_fac(:)
      rate(:,373) = 3.2e-11_r8 * exp_fac(:)
      rate(:,374) = 3.6e-9_r8 * exp_fac(:)
      rate(:,375) = 2.e-9_r8 * exp_fac(:)
      rate(:,376) = 1.e-10_r8 * exp_fac(:)
      rate(:,377) = 1.e-10_r8 * exp_fac(:)
      rate(:,378) = 1.5e-10_r8 * exp_fac(:)
      rate(:,379) = 7.8e-10_r8 * exp_fac(:)
      rate(:,380) = 9.9e-30_r8 * exp_fac(:)
      rate(:,381) = 7.e-10_r8 * exp_fac(:)
      rate(:,382) = 3.4e-31_r8 * exp_fac(:)
      rate(:,383) = 2.9e-9_r8 * exp_fac(:)
      rate(:,384) = 1.6e-9_r8 * exp_fac(:)
      rate(:,385) = 1.e-10_r8 * exp_fac(:)
      rate(:,386) = 1.e-10_r8 * exp_fac(:)
      rate(:,387) = 2.5e-10_r8 * exp_fac(:)
      rate(:,388) = 8.4e-10_r8 * exp_fac(:)
      rate(:,389) = 5.5e-10_r8 * exp_fac(:)
      rate(:,394) = 4.e-10_r8 * exp_fac(:)
      rate(:,395) = 4.3e-10_r8 * exp_fac(:)
      rate(:,396) = 9.e-10_r8 * exp_fac(:)
      rate(:,397) = 1.1e-9_r8 * exp_fac(:)
      rate(:,398) = 7.6e-28_r8 * exp_fac(:)
      rate(:,399) = 1.e-9_r8 * exp_fac(:)
      rate(:,400) = 1.e-10_r8 * exp_fac(:)
      rate(:,401) = 1.e-10_r8 * exp_fac(:)
      rate(:,402) = 1.1e-10_r8 * exp_fac(:)
      rate(:,403) = 6.0e-15_r8 * exp_fac(:)
      rate(:,404) = 1.7e-10_r8 * exp_fac(:)
      rate(:,407) = 3.51e-10_r8 * exp_fac(:)
      rate(:,408) = 1.e-10_r8 * exp_fac(:)
      rate(:,409) = 1.e-10_r8 * exp_fac(:)
      rate(:,410) = 1.e-11_r8 * exp_fac(:)
      rate(:,411) = 1.3e-10_r8 * exp_fac(:)
      rate(:,412) = 2.2e-10_r8 * exp_fac(:)
      rate(:,413) = 1.4e-10_r8 * exp_fac(:)
      rate(:,414) = 1.2e-9_r8 * exp_fac(:)
      rate(:,415) = 1.e-10_r8 * exp_fac(:)
      rate(:,416) = 1.0e-10_r8 * exp_fac(:)
      rate(:,417) = 3.e-10_r8 * exp_fac(:)
      rate(:,418) = 2.e-13_r8 * exp_fac(:)
      rate(:,419) = 1.2e-10_r8 * exp_fac(:)
      rate(:,420) = 1.6e-9_r8 * exp_fac(:)
      rate(:,421) = 1.4e-9_r8 * exp_fac(:)
      rate(:,422) = 1.0e-10_r8 * exp_fac(:)
      rate(:,423) = 1.0e-10_r8 * exp_fac(:)
      rate(:,424) = 0.5e-11_r8 * exp_fac(:)
      rate(:,425) = 1.e-13_r8 * exp_fac(:)
      rate(:,426) = 1.e-12_r8 * exp_fac(:)
      rate(:,427) = 9.6e-10_r8 * exp_fac(:)
      rate(:,428) = 6.0e-12_r8 * exp_fac(:)
      rate(:,429) = 1.6e-9_r8 * exp_fac(:)
      rate(:,430) = 2.e-29_r8 * exp_fac(:)
      rate(:,431) = 1.0e-27_r8 * exp_fac(:)
      rate(:,432) = 2.9e-12_r8 * exp_fac(:)
      rate(:,433) = 2.9e-11_r8 * exp_fac(:)
      rate(:,434) = 2.0e-10_r8 * exp_fac(:)
      rate(:,435) = 1.3e-9_r8 * exp_fac(:)
      rate(:,438) = 1.0e-28_r8 * exp_fac(:)
      rate(:,439) = 1.6e-28_r8 * exp_fac(:)
      rate(:,441) = 3.5e-12_r8 * exp_fac(:)
      rate(:,442) = 4.0e-11_r8 * exp_fac(:)
      rate(:,444) = 4.0e-11_r8 * exp_fac(:)
      rate(:,445) = 1.0e-28_r8 * exp_fac(:)
      rate(:,447) = 3.5e-12_r8 * exp_fac(:)
      rate(:,449) = 1.6e-28_r8 * exp_fac(:)
      rate(:,450) = 1.6e-28_r8 * exp_fac(:)
      rate(:,452) = 7.0e-10_r8 * exp_fac(:)
      rate(:,454) = 1.45e-26_r8 * exp_fac(:)
      rate(:,455) = 7.6e-10_r8 * exp_fac(:)
      rate(:,457) = 7.0e-10_r8 * exp_fac(:)
      rate(:,458) = 1.6e-9_r8 * exp_fac(:)
      rate(:,234) = 1.60e-10_r8 * exp( -260._r8 * itemp(:) )
      rate(:,235) = 6.00e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,236) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:) )
      rate(:,237) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:) )
      rate(:,238) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,239) = 1.05e-12_r8 * exp_fac(:)
      rate(:,242) = 1.25e-12_r8 * exp_fac(:)
      rate(:,253) = 3.40e-11_r8 * exp_fac(:)
      rate(:,240) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:) )
      rate(:,241) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:) )
      rate(:,243) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,245) = 1.35e-12_r8 * exp( -600._r8 * itemp(:) )
      rate(:,247) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,248) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,251) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,252) = 5.50e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,254) = 2.80e-12_r8 * exp( 300._r8 * itemp(:) )
      rate(:,255) = 4.10e-13_r8 * exp( 750._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.40e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,145), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.90e-31_r8 * itemp(:)**1.0_r8
      kinf(:) = 2.60e-11_r8
      call jpl( rate(:,154), m, 0.6_r8, ko, kinf, n )

      ko(:) = 7e-31_r8 * itemp(:)**2.6_r8
      kinf(:) = 3.6e-11_r8 * itemp(:)**0.1_r8
      call jpl( rate(:,162), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.00e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3.0e-11_r8
      call jpl( rate(:,172), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.50e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,176), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.00e-30_r8 * itemp(:)**4.4_r8
      kinf(:) = 1.4e-12_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,178), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.80e-30_r8 * itemp(:)**3.0_r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,180), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.00e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 2.9e-12_r8 * itemp(:)**1.1_r8
      call jpl( rate(:,186), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.80e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,202), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.60e-32_r8 * itemp(:)**4.5_r8
      kinf(:) = 3.0e-12_r8 * itemp(:)**2.0_r8
      call jpl( rate(:,206), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.20e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,223), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.90e-33_r8 * itemp(:)**1.4_r8
      kinf(:) = 1.10e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,250), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,104) = 8.00e-14_r8
      rate(:n,105) = 3.90e-17_r8
      rate(:n,110) = 1.30e-16_r8
      rate(:n,112) = 1.00e-20_r8
      rate(:n,148) = 6.90e-12_r8
      rate(:n,164) = 5.00e-12_r8
      rate(:n,165) = 7.00e-13_r8
      rate(:n,283) = 1.0e-10_r8
      rate(:n,284) = 4.4e-10_r8
      rate(:n,285) = 4.0e-10_r8
      rate(:n,286) = 2.0e-10_r8
      rate(:n,287) = 1.0e-12_r8
      rate(:n,288) = 6.0e-11_r8
      rate(:n,293) = 4.8e-10_r8
      rate(:n,294) = 1.0e-10_r8
      rate(:n,295) = 4.0e-10_r8
      rate(:n,298) = 5.0e-12_r8
      rate(:n,299) = 7.0e-10_r8
      rate(:n,300) = 8.0e-10_r8
      rate(:n,302) = 4.7e-2_r8
      rate(:n,303) = 1.71e-1_r8
      rate(:n,304) = 7.7e-5_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,102) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,106) = 1.80e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,107) = 3.50e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,111) = 3.60e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,114) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,115) = 3.135e-11_r8 * exp_fac(:)
      rate(:n,116) = 1.65e-12_r8 * exp_fac(:)
      rate(:n,146) = 1.40e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,150) = 1.80e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,151) = 1.70e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,152) = 4.80e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,158) = 3.00e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,159) = 1.00e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,167) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,168) = 2.10e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,173) = 3.30e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,174) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:) )
      rate(:n,175) = 5.10e-12_r8 * exp( 210._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.40e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,145) = wrk(:)












      end subroutine setrxt_hrates

      end module mo_setrxt
