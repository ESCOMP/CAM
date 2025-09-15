
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

      rate(:,92) = 0.000258_r8
      rate(:,93) = 0.085_r8
      rate(:,94) = 1.2e-10_r8
      rate(:,99) = 1.2e-10_r8
      rate(:,100) = 1e-20_r8
      rate(:,101) = 1.3e-16_r8
      rate(:,103) = 4.2e-13_r8
      rate(:,105) = 8e-14_r8
      rate(:,106) = 3.9e-17_r8
      rate(:,113) = 6.9e-12_r8
      rate(:,114) = 7.2e-11_r8
      rate(:,115) = 1.6e-12_r8
      rate(:,121) = 1.8e-12_r8
      rate(:,125) = 1.8e-12_r8
      rate(:,128) = 1.06e-05_r8
      rate(:,130) = 7e-11_r8
      rate(:,131) = 7e-13_r8
      rate(:,139) = 3.5e-12_r8
      rate(:,141) = 1.3e-11_r8
      rate(:,142) = 2.2e-11_r8
      rate(:,143) = 5e-11_r8
      rate(:,181) = 1.7e-13_r8
      rate(:,183) = 2.607e-10_r8
      rate(:,184) = 9.75e-11_r8
      rate(:,185) = 2.07e-10_r8
      rate(:,186) = 2.088e-10_r8
      rate(:,187) = 1.17e-10_r8
      rate(:,188) = 4.644e-11_r8
      rate(:,189) = 1.204e-10_r8
      rate(:,190) = 9.9e-11_r8
      rate(:,191) = 3.3e-12_r8
      rate(:,210) = 4.5e-11_r8
      rate(:,211) = 4.62e-10_r8
      rate(:,212) = 1.2e-10_r8
      rate(:,213) = 9e-11_r8
      rate(:,214) = 3e-11_r8
      rate(:,219) = 2.14e-11_r8
      rate(:,220) = 1.9e-10_r8
      rate(:,233) = 2.57e-10_r8
      rate(:,234) = 1.8e-10_r8
      rate(:,235) = 1.794e-10_r8
      rate(:,236) = 1.3e-10_r8
      rate(:,237) = 7.65e-11_r8
      rate(:,245) = 1.31e-10_r8
      rate(:,246) = 3.5e-11_r8
      rate(:,247) = 9e-12_r8
      rate(:,253) = 2.3e-12_r8
      rate(:,255) = 1.2e-11_r8
      rate(:,256) = 5.7e-11_r8
      rate(:,257) = 2.8e-11_r8
      rate(:,258) = 6.6e-11_r8
      rate(:,259) = 1.4e-11_r8
      rate(:,262) = 1.9e-12_r8
      rate(:,287) = 0.047_r8
      rate(:,288) = 7.7e-05_r8
      rate(:,289) = 0.171_r8
      rate(:,293) = 6e-11_r8
      rate(:,296) = 1e-12_r8
      rate(:,297) = 4e-10_r8
      rate(:,298) = 2e-10_r8
      rate(:,299) = 1e-10_r8
      rate(:,300) = 5e-16_r8
      rate(:,301) = 4.4e-10_r8
      rate(:,302) = 9e-10_r8
      rate(:,304) = 1.3e-10_r8
      rate(:,307) = 8e-10_r8
      rate(:,308) = 5e-12_r8
      rate(:,309) = 7e-10_r8
      rate(:,312) = 4.8e-10_r8
      rate(:,313) = 1e-10_r8
      rate(:,314) = 4e-10_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,95) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,96) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,97) = 2.64e-11_r8 * exp_fac(:)
      rate(:,98) = 6.6e-12_r8 * exp_fac(:)
      rate(:,102) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,104) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,107) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,108) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,111) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      rate(:,112) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,117) = 3e-11_r8 * exp_fac(:)
      rate(:,208) = 5.5e-12_r8 * exp_fac(:)
      rate(:,243) = 3.8e-12_r8 * exp_fac(:)
      rate(:,118) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,119) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,120) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,122) = 4.8e-11_r8 * exp_fac(:)
      rate(:,206) = 1.7e-11_r8 * exp_fac(:)
      rate(:,123) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:,124) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,129) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,132) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,133) = 2.9e-12_r8 * exp_fac(:)
      rate(:,134) = 1.45e-12_r8 * exp_fac(:)
      rate(:,135) = 1.45e-12_r8 * exp_fac(:)
      rate(:,136) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,137) = 1.2e-13_r8 * exp_fac(:)
      rate(:,166) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,140) = 1.7e-11_r8 * exp_fac(:)
      rate(:,240) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,144) = 3.44e-12_r8 * exp_fac(:)
      rate(:,199) = 2.3e-12_r8 * exp_fac(:)
      rate(:,202) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,145) = 3e-12_r8 * exp_fac(:)
      rate(:,207) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,147) = 7.26e-11_r8 * exp_fac(:)
      rate(:,148) = 4.64e-11_r8 * exp_fac(:)
      rate(:,158) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,159) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,160) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,161) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,162) = 1.4e-11_r8 * exp_fac(:)
      rate(:,176) = 7.4e-12_r8 * exp_fac(:)
      rate(:,163) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,164) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,165) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,167) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,168) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,169) = 2.6e-12_r8 * exp_fac(:)
      rate(:,170) = 6.4e-12_r8 * exp_fac(:)
      rate(:,200) = 4.1e-13_r8 * exp_fac(:)
      rate(:,171) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,173) = 3.6e-12_r8 * exp_fac(:)
      rate(:,222) = 2e-12_r8 * exp_fac(:)
      rate(:,174) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,175) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,177) = 6e-13_r8 * exp_fac(:)
      rate(:,197) = 1.5e-12_r8 * exp_fac(:)
      rate(:,205) = 1.9e-11_r8 * exp_fac(:)
      rate(:,178) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,179) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,180) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,182) = 3e-12_r8 * exp_fac(:)
      rate(:,216) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,194) = 1.7e-11_r8 * exp_fac(:)
      rate(:,221) = 6.3e-12_r8 * exp_fac(:)
      rate(:,195) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,196) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,198) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,201) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,204) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,209) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,215) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,217) = 1.4e-11_r8 * exp_fac(:)
      rate(:,219) = 2.14e-11_r8 * exp_fac(:)
      rate(:,220) = 1.9e-10_r8 * exp_fac(:)
      rate(:,233) = 2.57e-10_r8 * exp_fac(:)
      rate(:,234) = 1.8e-10_r8 * exp_fac(:)
      rate(:,235) = 1.794e-10_r8 * exp_fac(:)
      rate(:,236) = 1.3e-10_r8 * exp_fac(:)
      rate(:,237) = 7.65e-11_r8 * exp_fac(:)
      rate(:,245) = 1.31e-10_r8 * exp_fac(:)
      rate(:,246) = 3.5e-11_r8 * exp_fac(:)
      rate(:,247) = 9e-12_r8 * exp_fac(:)
      rate(:,253) = 2.3e-12_r8 * exp_fac(:)
      rate(:,255) = 1.2e-11_r8 * exp_fac(:)
      rate(:,256) = 5.7e-11_r8 * exp_fac(:)
      rate(:,257) = 2.8e-11_r8 * exp_fac(:)
      rate(:,258) = 6.6e-11_r8 * exp_fac(:)
      rate(:,259) = 1.4e-11_r8 * exp_fac(:)
      rate(:,262) = 1.9e-12_r8 * exp_fac(:)
      rate(:,287) = 0.047_r8 * exp_fac(:)
      rate(:,288) = 7.7e-05_r8 * exp_fac(:)
      rate(:,289) = 0.171_r8 * exp_fac(:)
      rate(:,293) = 6e-11_r8 * exp_fac(:)
      rate(:,296) = 1e-12_r8 * exp_fac(:)
      rate(:,297) = 4e-10_r8 * exp_fac(:)
      rate(:,298) = 2e-10_r8 * exp_fac(:)
      rate(:,299) = 1e-10_r8 * exp_fac(:)
      rate(:,300) = 5e-16_r8 * exp_fac(:)
      rate(:,301) = 4.4e-10_r8 * exp_fac(:)
      rate(:,302) = 9e-10_r8 * exp_fac(:)
      rate(:,304) = 1.3e-10_r8 * exp_fac(:)
      rate(:,307) = 8e-10_r8 * exp_fac(:)
      rate(:,308) = 5e-12_r8 * exp_fac(:)
      rate(:,309) = 7e-10_r8 * exp_fac(:)
      rate(:,312) = 4.8e-10_r8 * exp_fac(:)
      rate(:,313) = 1e-10_r8 * exp_fac(:)
      rate(:,314) = 4e-10_r8 * exp_fac(:)
      rate(:,218) = 6e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,223) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,224) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      rate(:,225) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:) )
      rate(:,226) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,227) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,228) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,229) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,230) = 1.25e-12_r8 * exp_fac(:)
      rate(:,239) = 3.4e-11_r8 * exp_fac(:)
      rate(:,231) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,232) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,238) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,241) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,242) = 2.8e-12_r8 * exp( 300._r8 * itemp(:) )
      rate(:,244) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,249) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,250) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,251) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,252) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,260) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,261) = 3.4e-12_r8 * exp( -1100._r8 * itemp(:) )
      rate(:,263) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,116), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,126), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,138), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,146), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,149), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,150), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,151), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,172), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,192), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,203), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,254), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,100) = 1e-20_r8
      rate(:n,101) = 1.3e-16_r8
      rate(:n,105) = 8e-14_r8
      rate(:n,106) = 3.9e-17_r8
      rate(:n,113) = 6.9e-12_r8
      rate(:n,130) = 7e-11_r8
      rate(:n,131) = 7e-13_r8
      rate(:n,287) = 0.047_r8
      rate(:n,288) = 7.7e-05_r8
      rate(:n,289) = 0.171_r8
      rate(:n,293) = 6e-11_r8
      rate(:n,296) = 1e-12_r8
      rate(:n,297) = 4e-10_r8
      rate(:n,298) = 2e-10_r8
      rate(:n,299) = 1e-10_r8
      rate(:n,301) = 4.4e-10_r8
      rate(:n,304) = 1.3e-10_r8
      rate(:n,307) = 8e-10_r8
      rate(:n,308) = 5e-12_r8
      rate(:n,309) = 7e-10_r8
      rate(:n,312) = 4.8e-10_r8
      rate(:n,313) = 1e-10_r8
      rate(:n,314) = 4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,96) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,97) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,98) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,102) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,104) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,107) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,108) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,117) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,118) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,119) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,122) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,123) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,124) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,132) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,136) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,144) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,145) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,116) = wrk(:)











      end subroutine setrxt_hrates

      end module mo_setrxt
