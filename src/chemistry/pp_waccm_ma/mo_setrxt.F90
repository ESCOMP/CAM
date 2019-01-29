
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

      rate(:,87) = 0.000258_r8
      rate(:,88) = 0.085_r8
      rate(:,89) = 1.2e-10_r8
      rate(:,94) = 1.2e-10_r8
      rate(:,95) = 1e-20_r8
      rate(:,96) = 1.3e-16_r8
      rate(:,98) = 4.2e-13_r8
      rate(:,100) = 8e-14_r8
      rate(:,101) = 3.9e-17_r8
      rate(:,108) = 6.9e-12_r8
      rate(:,109) = 7.2e-11_r8
      rate(:,110) = 1.6e-12_r8
      rate(:,116) = 1.8e-12_r8
      rate(:,120) = 1.8e-12_r8
      rate(:,124) = 7e-13_r8
      rate(:,125) = 5e-12_r8
      rate(:,134) = 3.5e-12_r8
      rate(:,136) = 1e-11_r8
      rate(:,137) = 2.2e-11_r8
      rate(:,138) = 5e-11_r8
      rate(:,173) = 1.7e-13_r8
      rate(:,175) = 2.607e-10_r8
      rate(:,176) = 9.75e-11_r8
      rate(:,177) = 2.07e-10_r8
      rate(:,178) = 2.088e-10_r8
      rate(:,179) = 1.17e-10_r8
      rate(:,180) = 4.644e-11_r8
      rate(:,181) = 1.204e-10_r8
      rate(:,182) = 9.9e-11_r8
      rate(:,183) = 3.3e-12_r8
      rate(:,202) = 4.5e-11_r8
      rate(:,203) = 4.62e-10_r8
      rate(:,204) = 1.2e-10_r8
      rate(:,205) = 9e-11_r8
      rate(:,206) = 3e-11_r8
      rate(:,211) = 2.14e-11_r8
      rate(:,212) = 1.9e-10_r8
      rate(:,225) = 2.57e-10_r8
      rate(:,226) = 1.8e-10_r8
      rate(:,227) = 1.794e-10_r8
      rate(:,228) = 1.3e-10_r8
      rate(:,229) = 7.65e-11_r8
      rate(:,238) = 1.31e-10_r8
      rate(:,239) = 3.5e-11_r8
      rate(:,240) = 9e-12_r8
      rate(:,263) = 0.047_r8
      rate(:,264) = 7.7e-05_r8
      rate(:,265) = 0.171_r8
      rate(:,269) = 6e-11_r8
      rate(:,272) = 1e-12_r8
      rate(:,273) = 4e-10_r8
      rate(:,274) = 2e-10_r8
      rate(:,275) = 1e-10_r8
      rate(:,276) = 5e-16_r8
      rate(:,277) = 4.4e-10_r8
      rate(:,278) = 9e-10_r8
      rate(:,280) = 1.3e-10_r8
      rate(:,283) = 8e-10_r8
      rate(:,284) = 5e-12_r8
      rate(:,285) = 7e-10_r8
      rate(:,288) = 4.8e-10_r8
      rate(:,289) = 1e-10_r8
      rate(:,290) = 4e-10_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,90) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,91) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,92) = 2.64e-11_r8 * exp_fac(:)
      rate(:,93) = 6.6e-12_r8 * exp_fac(:)
      rate(:,97) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,99) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,102) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,103) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,106) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      rate(:,107) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,112) = 3e-11_r8 * exp_fac(:)
      rate(:,200) = 5.5e-12_r8 * exp_fac(:)
      rate(:,235) = 3.8e-12_r8 * exp_fac(:)
      rate(:,113) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,114) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,115) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,117) = 4.8e-11_r8 * exp_fac(:)
      rate(:,198) = 1.7e-11_r8 * exp_fac(:)
      rate(:,118) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:,119) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,123) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,126) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,127) = 2.9e-12_r8 * exp_fac(:)
      rate(:,128) = 1.45e-12_r8 * exp_fac(:)
      rate(:,129) = 1.45e-12_r8 * exp_fac(:)
      rate(:,130) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,131) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,132) = 1.2e-13_r8 * exp_fac(:)
      rate(:,158) = 3e-11_r8 * exp_fac(:)
      rate(:,135) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,139) = 3.3e-12_r8 * exp_fac(:)
      rate(:,154) = 1.4e-11_r8 * exp_fac(:)
      rate(:,168) = 7.4e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,140) = 3e-12_r8 * exp_fac(:)
      rate(:,199) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,142) = 7.26e-11_r8 * exp_fac(:)
      rate(:,143) = 4.64e-11_r8 * exp_fac(:)
      rate(:,150) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,151) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,152) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,153) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,155) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,156) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,157) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,159) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,160) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,161) = 2.6e-12_r8 * exp_fac(:)
      rate(:,162) = 6.4e-12_r8 * exp_fac(:)
      rate(:,192) = 4.1e-13_r8 * exp_fac(:)
      rate(:,163) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,165) = 3.6e-12_r8 * exp_fac(:)
      rate(:,214) = 2e-12_r8 * exp_fac(:)
      rate(:,166) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,167) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,169) = 6e-13_r8 * exp_fac(:)
      rate(:,189) = 1.5e-12_r8 * exp_fac(:)
      rate(:,197) = 1.9e-11_r8 * exp_fac(:)
      rate(:,170) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,171) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,172) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,174) = 3e-12_r8 * exp_fac(:)
      rate(:,208) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,186) = 1.7e-11_r8 * exp_fac(:)
      rate(:,213) = 6.3e-12_r8 * exp_fac(:)
      rate(:,187) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,188) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,190) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,191) = 2.3e-12_r8 * exp_fac(:)
      rate(:,194) = 8.8e-12_r8 * exp_fac(:)
      rate(:,193) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,196) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,201) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,207) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,209) = 1.4e-11_r8 * exp_fac(:)
      rate(:,211) = 2.14e-11_r8 * exp_fac(:)
      rate(:,212) = 1.9e-10_r8 * exp_fac(:)
      rate(:,225) = 2.57e-10_r8 * exp_fac(:)
      rate(:,226) = 1.8e-10_r8 * exp_fac(:)
      rate(:,227) = 1.794e-10_r8 * exp_fac(:)
      rate(:,228) = 1.3e-10_r8 * exp_fac(:)
      rate(:,229) = 7.65e-11_r8 * exp_fac(:)
      rate(:,238) = 1.31e-10_r8 * exp_fac(:)
      rate(:,239) = 3.5e-11_r8 * exp_fac(:)
      rate(:,240) = 9e-12_r8 * exp_fac(:)
      rate(:,263) = 0.047_r8 * exp_fac(:)
      rate(:,264) = 7.7e-05_r8 * exp_fac(:)
      rate(:,265) = 0.171_r8 * exp_fac(:)
      rate(:,269) = 6e-11_r8 * exp_fac(:)
      rate(:,272) = 1e-12_r8 * exp_fac(:)
      rate(:,273) = 4e-10_r8 * exp_fac(:)
      rate(:,274) = 2e-10_r8 * exp_fac(:)
      rate(:,275) = 1e-10_r8 * exp_fac(:)
      rate(:,276) = 5e-16_r8 * exp_fac(:)
      rate(:,277) = 4.4e-10_r8 * exp_fac(:)
      rate(:,278) = 9e-10_r8 * exp_fac(:)
      rate(:,280) = 1.3e-10_r8 * exp_fac(:)
      rate(:,283) = 8e-10_r8 * exp_fac(:)
      rate(:,284) = 5e-12_r8 * exp_fac(:)
      rate(:,285) = 7e-10_r8 * exp_fac(:)
      rate(:,288) = 4.8e-10_r8 * exp_fac(:)
      rate(:,289) = 1e-10_r8 * exp_fac(:)
      rate(:,290) = 4e-10_r8 * exp_fac(:)
      rate(:,210) = 6e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,215) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,216) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      rate(:,217) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:) )
      rate(:,218) = 2.03e-11_r8 * exp( -1100._r8 * itemp(:) )
      rate(:,219) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,220) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,221) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,222) = 1.25e-12_r8 * exp_fac(:)
      rate(:,231) = 3.4e-11_r8 * exp_fac(:)
      rate(:,223) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,224) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,230) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,232) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,233) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,234) = 2.8e-12_r8 * exp( 300._r8 * itemp(:) )
      rate(:,236) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,111), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,121), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,133), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,141), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,144), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,145), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,146), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,164), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,184), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,195), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,237), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,95) = 1e-20_r8
      rate(:n,96) = 1.3e-16_r8
      rate(:n,100) = 8e-14_r8
      rate(:n,101) = 3.9e-17_r8
      rate(:n,108) = 6.9e-12_r8
      rate(:n,124) = 7e-13_r8
      rate(:n,125) = 5e-12_r8
      rate(:n,263) = 0.047_r8
      rate(:n,264) = 7.7e-05_r8
      rate(:n,265) = 0.171_r8
      rate(:n,269) = 6e-11_r8
      rate(:n,272) = 1e-12_r8
      rate(:n,273) = 4e-10_r8
      rate(:n,274) = 2e-10_r8
      rate(:n,275) = 1e-10_r8
      rate(:n,277) = 4.4e-10_r8
      rate(:n,280) = 1.3e-10_r8
      rate(:n,283) = 8e-10_r8
      rate(:n,284) = 5e-12_r8
      rate(:n,285) = 7e-10_r8
      rate(:n,288) = 4.8e-10_r8
      rate(:n,289) = 1e-10_r8
      rate(:n,290) = 4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,91) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,92) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,93) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,97) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,99) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,102) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,103) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,112) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,113) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,114) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,117) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,118) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,119) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,126) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,130) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,131) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,139) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,140) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,111) = wrk(:)











      end subroutine setrxt_hrates

      end module mo_setrxt
