
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

      rate(:,91) = 0.000258_r8
      rate(:,92) = 0.085_r8
      rate(:,93) = 1.2e-10_r8
      rate(:,98) = 1.2e-10_r8
      rate(:,99) = 1e-20_r8
      rate(:,100) = 1.3e-16_r8
      rate(:,102) = 4.2e-13_r8
      rate(:,104) = 8e-14_r8
      rate(:,105) = 3.9e-17_r8
      rate(:,112) = 6.9e-12_r8
      rate(:,113) = 7.2e-11_r8
      rate(:,114) = 1.6e-12_r8
      rate(:,120) = 1.8e-12_r8
      rate(:,124) = 1.8e-12_r8
      rate(:,128) = 7e-13_r8
      rate(:,129) = 5e-12_r8
      rate(:,138) = 3.5e-12_r8
      rate(:,140) = 1e-11_r8
      rate(:,141) = 2.2e-11_r8
      rate(:,142) = 5e-11_r8
      rate(:,177) = 1.7e-13_r8
      rate(:,179) = 2.607e-10_r8
      rate(:,180) = 9.75e-11_r8
      rate(:,181) = 2.07e-10_r8
      rate(:,182) = 2.088e-10_r8
      rate(:,183) = 1.17e-10_r8
      rate(:,184) = 4.644e-11_r8
      rate(:,185) = 1.204e-10_r8
      rate(:,186) = 9.9e-11_r8
      rate(:,187) = 3.3e-12_r8
      rate(:,206) = 4.5e-11_r8
      rate(:,207) = 4.62e-10_r8
      rate(:,208) = 1.2e-10_r8
      rate(:,209) = 9e-11_r8
      rate(:,210) = 3e-11_r8
      rate(:,215) = 2.14e-11_r8
      rate(:,216) = 1.9e-10_r8
      rate(:,229) = 2.57e-10_r8
      rate(:,230) = 1.8e-10_r8
      rate(:,231) = 1.794e-10_r8
      rate(:,232) = 1.3e-10_r8
      rate(:,233) = 7.65e-11_r8
      rate(:,242) = 1.31e-10_r8
      rate(:,243) = 3.5e-11_r8
      rate(:,244) = 9e-12_r8
      rate(:,248) = 2.3e-12_r8
      rate(:,249) = 1.2e-11_r8
      rate(:,250) = 5.7e-11_r8
      rate(:,251) = 2.8e-11_r8
      rate(:,252) = 6.6e-11_r8
      rate(:,253) = 1.4e-11_r8
      rate(:,256) = 1.9e-12_r8
      rate(:,287) = 6e-11_r8
      rate(:,290) = 1e-12_r8
      rate(:,291) = 4e-10_r8
      rate(:,292) = 2e-10_r8
      rate(:,293) = 1e-10_r8
      rate(:,294) = 5e-16_r8
      rate(:,295) = 4.4e-10_r8
      rate(:,296) = 9e-10_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,94) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,95) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,96) = 2.64e-11_r8 * exp_fac(:)
      rate(:,97) = 6.6e-12_r8 * exp_fac(:)
      rate(:,101) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,103) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,106) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,107) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,110) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      rate(:,111) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,116) = 3e-11_r8 * exp_fac(:)
      rate(:,204) = 5.5e-12_r8 * exp_fac(:)
      rate(:,239) = 3.8e-12_r8 * exp_fac(:)
      rate(:,117) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,118) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,119) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,121) = 4.8e-11_r8 * exp_fac(:)
      rate(:,202) = 1.7e-11_r8 * exp_fac(:)
      rate(:,122) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:,123) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,127) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,130) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,131) = 2.9e-12_r8 * exp_fac(:)
      rate(:,132) = 1.45e-12_r8 * exp_fac(:)
      rate(:,133) = 1.45e-12_r8 * exp_fac(:)
      rate(:,134) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,135) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,136) = 1.2e-13_r8 * exp_fac(:)
      rate(:,162) = 3e-11_r8 * exp_fac(:)
      rate(:,139) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,143) = 3.3e-12_r8 * exp_fac(:)
      rate(:,158) = 1.4e-11_r8 * exp_fac(:)
      rate(:,172) = 7.4e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,144) = 3e-12_r8 * exp_fac(:)
      rate(:,203) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,146) = 7.26e-11_r8 * exp_fac(:)
      rate(:,147) = 4.64e-11_r8 * exp_fac(:)
      rate(:,154) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,155) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,156) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,157) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,159) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,160) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,161) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,163) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,164) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,165) = 2.6e-12_r8 * exp_fac(:)
      rate(:,166) = 6.4e-12_r8 * exp_fac(:)
      rate(:,196) = 4.1e-13_r8 * exp_fac(:)
      rate(:,167) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,169) = 3.6e-12_r8 * exp_fac(:)
      rate(:,218) = 2e-12_r8 * exp_fac(:)
      rate(:,170) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,171) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,173) = 6e-13_r8 * exp_fac(:)
      rate(:,193) = 1.5e-12_r8 * exp_fac(:)
      rate(:,201) = 1.9e-11_r8 * exp_fac(:)
      rate(:,174) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,175) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,176) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,178) = 3e-12_r8 * exp_fac(:)
      rate(:,212) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,190) = 1.7e-11_r8 * exp_fac(:)
      rate(:,217) = 6.3e-12_r8 * exp_fac(:)
      rate(:,191) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,192) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,194) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,195) = 2.3e-12_r8 * exp_fac(:)
      rate(:,198) = 8.8e-12_r8 * exp_fac(:)
      rate(:,197) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,200) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,205) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,211) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,213) = 1.4e-11_r8 * exp_fac(:)
      rate(:,215) = 2.14e-11_r8 * exp_fac(:)
      rate(:,216) = 1.9e-10_r8 * exp_fac(:)
      rate(:,229) = 2.57e-10_r8 * exp_fac(:)
      rate(:,230) = 1.8e-10_r8 * exp_fac(:)
      rate(:,231) = 1.794e-10_r8 * exp_fac(:)
      rate(:,232) = 1.3e-10_r8 * exp_fac(:)
      rate(:,233) = 7.65e-11_r8 * exp_fac(:)
      rate(:,242) = 1.31e-10_r8 * exp_fac(:)
      rate(:,243) = 3.5e-11_r8 * exp_fac(:)
      rate(:,244) = 9e-12_r8 * exp_fac(:)
      rate(:,248) = 2.3e-12_r8 * exp_fac(:)
      rate(:,249) = 1.2e-11_r8 * exp_fac(:)
      rate(:,250) = 5.7e-11_r8 * exp_fac(:)
      rate(:,251) = 2.8e-11_r8 * exp_fac(:)
      rate(:,252) = 6.6e-11_r8 * exp_fac(:)
      rate(:,253) = 1.4e-11_r8 * exp_fac(:)
      rate(:,256) = 1.9e-12_r8 * exp_fac(:)
      rate(:,287) = 6e-11_r8 * exp_fac(:)
      rate(:,290) = 1e-12_r8 * exp_fac(:)
      rate(:,291) = 4e-10_r8 * exp_fac(:)
      rate(:,292) = 2e-10_r8 * exp_fac(:)
      rate(:,293) = 1e-10_r8 * exp_fac(:)
      rate(:,294) = 5e-16_r8 * exp_fac(:)
      rate(:,295) = 4.4e-10_r8 * exp_fac(:)
      rate(:,296) = 9e-10_r8 * exp_fac(:)
      rate(:,214) = 6e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,219) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,220) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      rate(:,221) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:) )
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,222) = 2.03e-11_r8 * exp_fac(:)
      rate(:,255) = 3.4e-12_r8 * exp_fac(:)
      rate(:,223) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,224) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,225) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,226) = 1.25e-12_r8 * exp_fac(:)
      rate(:,235) = 3.4e-11_r8 * exp_fac(:)
      rate(:,227) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,228) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,234) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,236) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,237) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,238) = 2.8e-12_r8 * exp( 300._r8 * itemp(:) )
      rate(:,240) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,246) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,247) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,254) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,257) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,260) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,261) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,115), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,125), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,137), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,145), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,148), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,149), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,150), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,168), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,188), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,199), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,241), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,99) = 1e-20_r8
      rate(:n,100) = 1.3e-16_r8
      rate(:n,104) = 8e-14_r8
      rate(:n,105) = 3.9e-17_r8
      rate(:n,112) = 6.9e-12_r8
      rate(:n,128) = 7e-13_r8
      rate(:n,129) = 5e-12_r8
      rate(:n,287) = 6e-11_r8
      rate(:n,290) = 1e-12_r8
      rate(:n,291) = 4e-10_r8
      rate(:n,292) = 2e-10_r8
      rate(:n,293) = 1e-10_r8
      rate(:n,295) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,95) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,96) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,97) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,101) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,103) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,106) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,107) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,116) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,117) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,118) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,121) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,122) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,123) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,130) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,134) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,135) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,143) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,144) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,115) = wrk(:)











      end subroutine setrxt_hrates

      end module mo_setrxt
