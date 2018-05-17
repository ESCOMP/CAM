
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

      rate(:,:,95) = 8.00e-14_r8
      rate(:,:,96) = 3.90e-17_r8
      rate(:,:,99) = 4.20e-13_r8
      rate(:,:,100) = 8.50e-2_r8
      rate(:,:,101) = 1.30e-16_r8
      rate(:,:,103) = 1.00e-20_r8
      rate(:,:,104) = 2.58e-04_r8
      rate(:,:,111) = 1.20e-10_r8
      rate(:,:,112) = 2.02e-10_r8
      rate(:,:,113) = 1.204e-10_r8
      rate(:,:,114) = 1.50e-10_r8
      rate(:,:,115) = 9.75e-11_r8
      rate(:,:,116) = 1.50e-11_r8
      rate(:,:,117) = 7.20e-11_r8
      rate(:,:,118) = 1.794e-10_r8
      rate(:,:,119) = 1.628e-10_r8
      rate(:,:,120) = 2.84e-10_r8
      rate(:,:,121) = 1.674e-10_r8
      rate(:,:,122) = 9.60e-11_r8
      rate(:,:,123) = 4.10e-11_r8
      rate(:,:,124) = 1.012e-10_r8
      rate(:,:,125) = 1.20e-10_r8
      rate(:,:,126) = 4.49e-10_r8
      rate(:,:,127) = 2.57e-10_r8
      rate(:,:,128) = 2.14e-11_r8
      rate(:,:,129) = 1.90e-10_r8
      rate(:,:,130) = 1.31e-10_r8
      rate(:,:,131) = 3.50e-11_r8
      rate(:,:,132) = 9.00e-12_r8
      rate(:,:,133) = 1.20e-10_r8
      rate(:,:,134) = 1.50e-10_r8
      rate(:,:,135) = 1.20e-10_r8
      rate(:,:,138) = 7.20e-11_r8
      rate(:,:,139) = 6.90e-12_r8
      rate(:,:,140) = 1.60e-12_r8
      rate(:,:,144) = 1.80e-12_r8
      rate(:,:,147) = 1.80e-12_r8
      rate(:,:,153) = 5.00e-12_r8
      rate(:,:,154) = 7.00e-13_r8
      rate(:,:,155) = 5.00e-11_r8
      rate(:,:,172) = 1.00e-11_r8
      rate(:,:,173) = 2.20e-11_r8
      rate(:,:,174) = 3.50e-12_r8
      rate(:,:,199) = 1.70e-13_r8
      rate(:,:,252) = 6.60E-11_r8
      rate(:,:,253) = 2.30E-12_r8
      rate(:,:,254) = 1.20E-11_r8
      rate(:,:,258) = 1.40E-11_r8
      rate(:,:,259) = 2.80E-11_r8
      rate(:,:,260) = 5.70E-11_r8
      rate(:,:,261) = 1.90E-12_r8
      rate(:,:,288) = 9.0e-10_r8
      rate(:,:,289) = 1.0e-10_r8
      rate(:,:,290) = 4.4e-10_r8
      rate(:,:,291) = 4.0e-10_r8
      rate(:,:,292) = 2.0e-10_r8
      rate(:,:,293) = 1.0e-12_r8
      rate(:,:,294) = 6.0e-11_r8
      rate(:,:,295) = 5.0e-16_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,93) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,97) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,98) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,102) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,105) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,106) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,107) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,108) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,109) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,110) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,137) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,141) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:,142) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,143) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,209) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,146) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,148) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,149) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,217) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,245) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,150) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,152) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,156) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,157) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,158) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,159) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,160) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,162) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,181) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,186) = 7.40e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,163) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,218) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,164) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,166) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,192) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,171) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,176) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,178) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,179) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,180) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,182) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,183) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,184) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,185) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,187) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,208) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,216) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,188) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,215) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,189) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,193) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,194) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,197) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,198) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,200) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,201) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,222) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,202) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,233) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,203) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,204) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,205) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,206) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,207) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,235) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,210) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,211) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,214) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,213) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,219) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,220) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,221) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,252) = 6.60E-11_r8 * exp_fac(:,:)
      rate(:,:,253) = 2.30E-12_r8 * exp_fac(:,:)
      rate(:,:,254) = 1.20E-11_r8 * exp_fac(:,:)
      rate(:,:,258) = 1.40E-11_r8 * exp_fac(:,:)
      rate(:,:,259) = 2.80E-11_r8 * exp_fac(:,:)
      rate(:,:,260) = 5.70E-11_r8 * exp_fac(:,:)
      rate(:,:,261) = 1.90E-12_r8 * exp_fac(:,:)
      rate(:,:,288) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,289) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,290) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,291) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,292) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,293) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,294) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,295) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,223) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      rate(:,:,224) = 6.00e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,225) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,226) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,227) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,228) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,231) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,242) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,229) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,230) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,232) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,234) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,236) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,237) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,240) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,241) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,243) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,244) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,250) = 2.10E-11_r8 * exp( -2200.0_r8 * itemp(:,:) )
      rate(:,:,251) = 1.10E-13_r8 * exp( -1200.0_r8 * itemp(:,:) )
      rate(:,:,255) = 2.70E-11_r8 * exp( 335._r8 * itemp(:,:) )
      rate(:,:,256) = 1.25E-13_r8 * exp( -2190.0_r8 * itemp(:,:) )
      rate(:,:,257) = 3.40E-12_r8 * exp( -1100.0_r8 * itemp(:,:) )
      rate(:,:,265) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,266) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,136), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,145), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,161), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,165), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,167), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,169), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,175), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,191), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,195), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,212), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,239), m, 0.6_r8, ko, kinf, n )

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

      rate(:,:kbot,95) = 8.00e-14_r8
      rate(:,:kbot,96) = 3.90e-17_r8
      rate(:,:kbot,101) = 1.30e-16_r8
      rate(:,:kbot,103) = 1.00e-20_r8
      rate(:,:kbot,139) = 6.90e-12_r8
      rate(:,:kbot,153) = 5.00e-12_r8
      rate(:,:kbot,154) = 7.00e-13_r8
      rate(:,:kbot,289) = 1.0e-10_r8
      rate(:,:kbot,290) = 4.4e-10_r8
      rate(:,:kbot,291) = 4.0e-10_r8
      rate(:,:kbot,292) = 2.0e-10_r8
      rate(:,:kbot,293) = 1.0e-12_r8
      rate(:,:kbot,294) = 6.0e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,93) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,97) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,98) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,102) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,105) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,106) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,107) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,137) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,141) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,142) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,143) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,149) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,150) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,156) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,157) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,162) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,163) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,164) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,136) = wrk(:,:)











      end subroutine setrxt_hrates

      end module mo_setrxt
