
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,108) = 8.00e-14_r8
      rate(:,:,109) = 3.90e-17_r8
      rate(:,:,112) = 4.20e-13_r8
      rate(:,:,113) = 8.50e-2_r8
      rate(:,:,114) = 1.30e-16_r8
      rate(:,:,116) = 1.00e-20_r8
      rate(:,:,117) = 2.58e-04_r8
      rate(:,:,124) = 1.20e-10_r8
      rate(:,:,125) = 2.02e-10_r8
      rate(:,:,126) = 1.204e-10_r8
      rate(:,:,127) = 1.50e-10_r8
      rate(:,:,128) = 9.75e-11_r8
      rate(:,:,129) = 1.50e-11_r8
      rate(:,:,130) = 7.20e-11_r8
      rate(:,:,131) = 1.794e-10_r8
      rate(:,:,132) = 1.628e-10_r8
      rate(:,:,133) = 2.84e-10_r8
      rate(:,:,134) = 1.674e-10_r8
      rate(:,:,135) = 9.60e-11_r8
      rate(:,:,136) = 4.10e-11_r8
      rate(:,:,137) = 1.012e-10_r8
      rate(:,:,138) = 1.20e-10_r8
      rate(:,:,139) = 4.49e-10_r8
      rate(:,:,140) = 2.57e-10_r8
      rate(:,:,141) = 2.14e-11_r8
      rate(:,:,142) = 1.90e-10_r8
      rate(:,:,143) = 1.31e-10_r8
      rate(:,:,144) = 3.50e-11_r8
      rate(:,:,145) = 9.00e-12_r8
      rate(:,:,146) = 1.20e-10_r8
      rate(:,:,147) = 1.50e-10_r8
      rate(:,:,148) = 1.20e-10_r8
      rate(:,:,151) = 7.20e-11_r8
      rate(:,:,152) = 6.90e-12_r8
      rate(:,:,153) = 1.60e-12_r8
      rate(:,:,157) = 1.80e-12_r8
      rate(:,:,160) = 1.80e-12_r8
      rate(:,:,168) = 5.00e-12_r8
      rate(:,:,169) = 7.00e-13_r8
      rate(:,:,170) = 5.00e-11_r8
      rate(:,:,187) = 1.00e-11_r8
      rate(:,:,188) = 2.20e-11_r8
      rate(:,:,189) = 3.50e-12_r8
      rate(:,:,214) = 1.70e-13_r8
      rate(:,:,267) = 6.60E-11_r8
      rate(:,:,268) = 2.30E-12_r8
      rate(:,:,269) = 1.20E-11_r8
      rate(:,:,273) = 1.40E-11_r8
      rate(:,:,274) = 2.80E-11_r8
      rate(:,:,275) = 5.70E-11_r8
      rate(:,:,276) = 1.90E-12_r8
      rate(:,:,303) = 9.0e-10_r8
      rate(:,:,304) = 1.0e-10_r8
      rate(:,:,305) = 4.4e-10_r8
      rate(:,:,306) = 4.0e-10_r8
      rate(:,:,307) = 2.0e-10_r8
      rate(:,:,308) = 1.0e-12_r8
      rate(:,:,309) = 6.0e-11_r8
      rate(:,:,310) = 5.0e-16_r8
      rate(:,:,318) = 1.5e-6_r8
      rate(:,:,319) = 1.5e-6_r8
      rate(:,:,320) = 2.0e-6_r8
      rate(:,:,321) = 2.0e-6_r8
      rate(:,:,322) = 2.0e-6_r8
      rate(:,:,323) = 1.5e-6_r8
      rate(:,:,327) = 3.6e-6_r8
      rate(:,:,328) = 5.0e-6_r8
      rate(:,:,330) = 1e-9_r8
      rate(:,:,331) = 1e-9_r8
      rate(:,:,333) = 2.8e-28_r8
      rate(:,:,334) = 1.7e-9_r8
      rate(:,:,335) = 1.5e-10_r8
      rate(:,:,336) = 3.e-10_r8
      rate(:,:,337) = 9.0e-10_r8
      rate(:,:,338) = 2.4e-10_r8
      rate(:,:,339) = 2.0e-9_r8
      rate(:,:,348) = 4.0e-12_r8
      rate(:,:,349) = 7.0e-12_r8
      rate(:,:,350) = 1.0e-9_r8
      rate(:,:,351) = 1.0e-9_r8
      rate(:,:,354) = 0.5e-9_r8
      rate(:,:,355) = 1.e-10_r8
      rate(:,:,356) = 7.e-12_r8
      rate(:,:,358) = 7e-11_r8
      rate(:,:,359) = 1.e-9_r8
      rate(:,:,364) = 1.9e-10_r8
      rate(:,:,366) = 3.e-10_r8
      rate(:,:,367) = 0.5e-12_r8
      rate(:,:,368) = 5.8e-10_r8
      rate(:,:,369) = 1.5e-10_r8
      rate(:,:,370) = 2.e-10_r8
      rate(:,:,372) = 1.4e-9_r8
      rate(:,:,373) = 1.e-10_r8
      rate(:,:,374) = 1.e-10_r8
      rate(:,:,375) = 2.e-10_r8
      rate(:,:,376) = 1.4e-9_r8
      rate(:,:,377) = 8.0e-10_r8
      rate(:,:,378) = 2.9e-31_r8
      rate(:,:,379) = 6.0e-13_r8
      rate(:,:,380) = 1.0e-9_r8
      rate(:,:,381) = 2.0e-28_r8
      rate(:,:,382) = 3.2e-11_r8
      rate(:,:,383) = 3.6e-9_r8
      rate(:,:,384) = 2.e-9_r8
      rate(:,:,385) = 1.e-10_r8
      rate(:,:,386) = 1.e-10_r8
      rate(:,:,387) = 1.5e-10_r8
      rate(:,:,388) = 7.8e-10_r8
      rate(:,:,389) = 9.9e-30_r8
      rate(:,:,390) = 7.e-10_r8
      rate(:,:,391) = 3.4e-31_r8
      rate(:,:,392) = 2.9e-9_r8
      rate(:,:,393) = 1.6e-9_r8
      rate(:,:,394) = 1.e-10_r8
      rate(:,:,395) = 1.e-10_r8
      rate(:,:,396) = 2.5e-10_r8
      rate(:,:,397) = 8.4e-10_r8
      rate(:,:,398) = 5.5e-10_r8
      rate(:,:,403) = 4.e-10_r8
      rate(:,:,404) = 4.3e-10_r8
      rate(:,:,405) = 9.e-10_r8
      rate(:,:,406) = 1.1e-9_r8
      rate(:,:,407) = 7.6e-28_r8
      rate(:,:,408) = 1.e-9_r8
      rate(:,:,409) = 1.e-10_r8
      rate(:,:,410) = 1.e-10_r8
      rate(:,:,411) = 1.1e-10_r8
      rate(:,:,412) = 6.0e-15_r8
      rate(:,:,413) = 1.7e-10_r8
      rate(:,:,416) = 3.51e-10_r8
      rate(:,:,417) = 1.e-10_r8
      rate(:,:,418) = 1.e-10_r8
      rate(:,:,419) = 1.e-11_r8
      rate(:,:,420) = 1.3e-10_r8
      rate(:,:,421) = 2.2e-10_r8
      rate(:,:,422) = 1.4e-10_r8
      rate(:,:,423) = 1.2e-9_r8
      rate(:,:,424) = 1.e-10_r8
      rate(:,:,425) = 1.0e-10_r8
      rate(:,:,426) = 3.e-10_r8
      rate(:,:,427) = 2.e-13_r8
      rate(:,:,428) = 1.2e-10_r8
      rate(:,:,429) = 1.6e-9_r8
      rate(:,:,430) = 1.4e-9_r8
      rate(:,:,431) = 1.0e-10_r8
      rate(:,:,432) = 1.0e-10_r8
      rate(:,:,433) = 0.5e-11_r8
      rate(:,:,434) = 1.e-13_r8
      rate(:,:,435) = 1.e-12_r8
      rate(:,:,436) = 9.6e-10_r8
      rate(:,:,437) = 6.0e-12_r8
      rate(:,:,438) = 1.6e-9_r8
      rate(:,:,439) = 2.e-29_r8
      rate(:,:,440) = 1.0e-27_r8
      rate(:,:,441) = 2.9e-12_r8
      rate(:,:,442) = 2.9e-11_r8
      rate(:,:,443) = 2.0e-10_r8
      rate(:,:,444) = 1.3e-9_r8
      rate(:,:,447) = 1.0e-28_r8
      rate(:,:,448) = 1.6e-28_r8
      rate(:,:,450) = 3.5e-12_r8
      rate(:,:,451) = 4.0e-11_r8
      rate(:,:,453) = 4.0e-11_r8
      rate(:,:,454) = 1.0e-28_r8
      rate(:,:,456) = 3.5e-12_r8
      rate(:,:,458) = 1.6e-28_r8
      rate(:,:,459) = 1.6e-28_r8
      rate(:,:,461) = 7.0e-10_r8
      rate(:,:,463) = 1.45e-26_r8
      rate(:,:,464) = 7.6e-10_r8
      rate(:,:,466) = 7.0e-10_r8
      rate(:,:,467) = 1.6e-9_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,106) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,110) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,111) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,115) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,118) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,119) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,120) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,121) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,122) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,123) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,150) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,154) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:,155) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,156) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,224) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,159) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,161) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,162) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,232) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,260) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,163) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,165) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,167) = 1.8e-11_r8 * exp( 390._r8 * itemp(:,:) )
      rate(:,:,171) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,172) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,173) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,174) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,175) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,177) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,196) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,201) = 7.40e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,178) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,233) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,179) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,181) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,207) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,186) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,191) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,193) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,194) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,195) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,197) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,198) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,199) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,200) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,202) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,223) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,231) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,203) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,205) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,230) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,204) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,208) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,209) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,212) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,213) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,215) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,216) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,237) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,217) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,248) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,218) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,219) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,220) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,221) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,222) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,250) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,225) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,226) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,229) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,228) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,234) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,235) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,236) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,267) = 6.60E-11_r8 * exp_fac(:,:)
      rate(:,:,268) = 2.30E-12_r8 * exp_fac(:,:)
      rate(:,:,269) = 1.20E-11_r8 * exp_fac(:,:)
      rate(:,:,273) = 1.40E-11_r8 * exp_fac(:,:)
      rate(:,:,274) = 2.80E-11_r8 * exp_fac(:,:)
      rate(:,:,275) = 5.70E-11_r8 * exp_fac(:,:)
      rate(:,:,276) = 1.90E-12_r8 * exp_fac(:,:)
      rate(:,:,303) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,304) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,305) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,306) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,307) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,308) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,309) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,310) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,318) = 1.5e-6_r8 * exp_fac(:,:)
      rate(:,:,319) = 1.5e-6_r8 * exp_fac(:,:)
      rate(:,:,320) = 2.0e-6_r8 * exp_fac(:,:)
      rate(:,:,321) = 2.0e-6_r8 * exp_fac(:,:)
      rate(:,:,322) = 2.0e-6_r8 * exp_fac(:,:)
      rate(:,:,323) = 1.5e-6_r8 * exp_fac(:,:)
      rate(:,:,327) = 3.6e-6_r8 * exp_fac(:,:)
      rate(:,:,328) = 5.0e-6_r8 * exp_fac(:,:)
      rate(:,:,330) = 1e-9_r8 * exp_fac(:,:)
      rate(:,:,331) = 1e-9_r8 * exp_fac(:,:)
      rate(:,:,333) = 2.8e-28_r8 * exp_fac(:,:)
      rate(:,:,334) = 1.7e-9_r8 * exp_fac(:,:)
      rate(:,:,335) = 1.5e-10_r8 * exp_fac(:,:)
      rate(:,:,336) = 3.e-10_r8 * exp_fac(:,:)
      rate(:,:,337) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,338) = 2.4e-10_r8 * exp_fac(:,:)
      rate(:,:,339) = 2.0e-9_r8 * exp_fac(:,:)
      rate(:,:,348) = 4.0e-12_r8 * exp_fac(:,:)
      rate(:,:,349) = 7.0e-12_r8 * exp_fac(:,:)
      rate(:,:,350) = 1.0e-9_r8 * exp_fac(:,:)
      rate(:,:,351) = 1.0e-9_r8 * exp_fac(:,:)
      rate(:,:,354) = 0.5e-9_r8 * exp_fac(:,:)
      rate(:,:,355) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,356) = 7.e-12_r8 * exp_fac(:,:)
      rate(:,:,358) = 7e-11_r8 * exp_fac(:,:)
      rate(:,:,359) = 1.e-9_r8 * exp_fac(:,:)
      rate(:,:,364) = 1.9e-10_r8 * exp_fac(:,:)
      rate(:,:,366) = 3.e-10_r8 * exp_fac(:,:)
      rate(:,:,367) = 0.5e-12_r8 * exp_fac(:,:)
      rate(:,:,368) = 5.8e-10_r8 * exp_fac(:,:)
      rate(:,:,369) = 1.5e-10_r8 * exp_fac(:,:)
      rate(:,:,370) = 2.e-10_r8 * exp_fac(:,:)
      rate(:,:,372) = 1.4e-9_r8 * exp_fac(:,:)
      rate(:,:,373) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,374) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,375) = 2.e-10_r8 * exp_fac(:,:)
      rate(:,:,376) = 1.4e-9_r8 * exp_fac(:,:)
      rate(:,:,377) = 8.0e-10_r8 * exp_fac(:,:)
      rate(:,:,378) = 2.9e-31_r8 * exp_fac(:,:)
      rate(:,:,379) = 6.0e-13_r8 * exp_fac(:,:)
      rate(:,:,380) = 1.0e-9_r8 * exp_fac(:,:)
      rate(:,:,381) = 2.0e-28_r8 * exp_fac(:,:)
      rate(:,:,382) = 3.2e-11_r8 * exp_fac(:,:)
      rate(:,:,383) = 3.6e-9_r8 * exp_fac(:,:)
      rate(:,:,384) = 2.e-9_r8 * exp_fac(:,:)
      rate(:,:,385) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,386) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,387) = 1.5e-10_r8 * exp_fac(:,:)
      rate(:,:,388) = 7.8e-10_r8 * exp_fac(:,:)
      rate(:,:,389) = 9.9e-30_r8 * exp_fac(:,:)
      rate(:,:,390) = 7.e-10_r8 * exp_fac(:,:)
      rate(:,:,391) = 3.4e-31_r8 * exp_fac(:,:)
      rate(:,:,392) = 2.9e-9_r8 * exp_fac(:,:)
      rate(:,:,393) = 1.6e-9_r8 * exp_fac(:,:)
      rate(:,:,394) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,395) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,396) = 2.5e-10_r8 * exp_fac(:,:)
      rate(:,:,397) = 8.4e-10_r8 * exp_fac(:,:)
      rate(:,:,398) = 5.5e-10_r8 * exp_fac(:,:)
      rate(:,:,403) = 4.e-10_r8 * exp_fac(:,:)
      rate(:,:,404) = 4.3e-10_r8 * exp_fac(:,:)
      rate(:,:,405) = 9.e-10_r8 * exp_fac(:,:)
      rate(:,:,406) = 1.1e-9_r8 * exp_fac(:,:)
      rate(:,:,407) = 7.6e-28_r8 * exp_fac(:,:)
      rate(:,:,408) = 1.e-9_r8 * exp_fac(:,:)
      rate(:,:,409) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,410) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,411) = 1.1e-10_r8 * exp_fac(:,:)
      rate(:,:,412) = 6.0e-15_r8 * exp_fac(:,:)
      rate(:,:,413) = 1.7e-10_r8 * exp_fac(:,:)
      rate(:,:,416) = 3.51e-10_r8 * exp_fac(:,:)
      rate(:,:,417) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,418) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,419) = 1.e-11_r8 * exp_fac(:,:)
      rate(:,:,420) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,421) = 2.2e-10_r8 * exp_fac(:,:)
      rate(:,:,422) = 1.4e-10_r8 * exp_fac(:,:)
      rate(:,:,423) = 1.2e-9_r8 * exp_fac(:,:)
      rate(:,:,424) = 1.e-10_r8 * exp_fac(:,:)
      rate(:,:,425) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,426) = 3.e-10_r8 * exp_fac(:,:)
      rate(:,:,427) = 2.e-13_r8 * exp_fac(:,:)
      rate(:,:,428) = 1.2e-10_r8 * exp_fac(:,:)
      rate(:,:,429) = 1.6e-9_r8 * exp_fac(:,:)
      rate(:,:,430) = 1.4e-9_r8 * exp_fac(:,:)
      rate(:,:,431) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,432) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,433) = 0.5e-11_r8 * exp_fac(:,:)
      rate(:,:,434) = 1.e-13_r8 * exp_fac(:,:)
      rate(:,:,435) = 1.e-12_r8 * exp_fac(:,:)
      rate(:,:,436) = 9.6e-10_r8 * exp_fac(:,:)
      rate(:,:,437) = 6.0e-12_r8 * exp_fac(:,:)
      rate(:,:,438) = 1.6e-9_r8 * exp_fac(:,:)
      rate(:,:,439) = 2.e-29_r8 * exp_fac(:,:)
      rate(:,:,440) = 1.0e-27_r8 * exp_fac(:,:)
      rate(:,:,441) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,442) = 2.9e-11_r8 * exp_fac(:,:)
      rate(:,:,443) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,444) = 1.3e-9_r8 * exp_fac(:,:)
      rate(:,:,447) = 1.0e-28_r8 * exp_fac(:,:)
      rate(:,:,448) = 1.6e-28_r8 * exp_fac(:,:)
      rate(:,:,450) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,451) = 4.0e-11_r8 * exp_fac(:,:)
      rate(:,:,453) = 4.0e-11_r8 * exp_fac(:,:)
      rate(:,:,454) = 1.0e-28_r8 * exp_fac(:,:)
      rate(:,:,456) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,458) = 1.6e-28_r8 * exp_fac(:,:)
      rate(:,:,459) = 1.6e-28_r8 * exp_fac(:,:)
      rate(:,:,461) = 7.0e-10_r8 * exp_fac(:,:)
      rate(:,:,463) = 1.45e-26_r8 * exp_fac(:,:)
      rate(:,:,464) = 7.6e-10_r8 * exp_fac(:,:)
      rate(:,:,466) = 7.0e-10_r8 * exp_fac(:,:)
      rate(:,:,467) = 1.6e-9_r8 * exp_fac(:,:)
      rate(:,:,238) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      rate(:,:,239) = 6.00e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,240) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,241) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,242) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,243) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,246) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,257) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,244) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,245) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,247) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,249) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,251) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,252) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,255) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,256) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,258) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,259) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,265) = 2.10E-11_r8 * exp( -2200.0_r8 * itemp(:,:) )
      rate(:,:,266) = 1.10E-13_r8 * exp( -1200.0_r8 * itemp(:,:) )
      rate(:,:,270) = 2.70E-11_r8 * exp( 335._r8 * itemp(:,:) )
      rate(:,:,271) = 1.25E-13_r8 * exp( -2190.0_r8 * itemp(:,:) )
      rate(:,:,272) = 3.40E-12_r8 * exp( -1100.0_r8 * itemp(:,:) )
      rate(:,:,280) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,281) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,149), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,158), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 7e-31_r8 * itemp(:,:)**2.6_r8
      kinf(:,:) = 3.6e-11_r8 * itemp(:,:)**0.1_r8
      call jpl( rate(1,1,166), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,176), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,180), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,182), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,184), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,190), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,206), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,210), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,227), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,254), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)

      rate(:,:kbot,108) = 8.00e-14_r8
      rate(:,:kbot,109) = 3.90e-17_r8
      rate(:,:kbot,114) = 1.30e-16_r8
      rate(:,:kbot,116) = 1.00e-20_r8
      rate(:,:kbot,152) = 6.90e-12_r8
      rate(:,:kbot,168) = 5.00e-12_r8
      rate(:,:kbot,169) = 7.00e-13_r8
      rate(:,:kbot,304) = 1.0e-10_r8
      rate(:,:kbot,305) = 4.4e-10_r8
      rate(:,:kbot,306) = 4.0e-10_r8
      rate(:,:kbot,307) = 2.0e-10_r8
      rate(:,:kbot,308) = 1.0e-12_r8
      rate(:,:kbot,309) = 6.0e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,106) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,110) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,111) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,115) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,118) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,119) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,120) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,150) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,154) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,155) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,156) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,162) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,163) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,171) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,172) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,177) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,178) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,179) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,149) = wrk(:,:)












      end subroutine setrxt_hrates

      end module mo_setrxt
