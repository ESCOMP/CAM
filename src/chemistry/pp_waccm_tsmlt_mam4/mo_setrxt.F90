
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

      rate(:,151) = 0.000258_r8
      rate(:,152) = 0.085_r8
      rate(:,153) = 1.2e-10_r8
      rate(:,158) = 1.2e-10_r8
      rate(:,159) = 1e-20_r8
      rate(:,160) = 1.3e-16_r8
      rate(:,162) = 4.2e-13_r8
      rate(:,164) = 8e-14_r8
      rate(:,165) = 3.9e-17_r8
      rate(:,172) = 6.9e-12_r8
      rate(:,173) = 7.2e-11_r8
      rate(:,174) = 1.6e-12_r8
      rate(:,180) = 1.8e-12_r8
      rate(:,184) = 1.8e-12_r8
      rate(:,188) = 7e-13_r8
      rate(:,189) = 5e-12_r8
      rate(:,198) = 3.5e-12_r8
      rate(:,200) = 1e-11_r8
      rate(:,201) = 2.2e-11_r8
      rate(:,202) = 5e-11_r8
      rate(:,237) = 1.7e-13_r8
      rate(:,239) = 2.607e-10_r8
      rate(:,240) = 9.75e-11_r8
      rate(:,241) = 2.07e-10_r8
      rate(:,242) = 2.088e-10_r8
      rate(:,243) = 1.17e-10_r8
      rate(:,244) = 4.644e-11_r8
      rate(:,245) = 1.204e-10_r8
      rate(:,246) = 9.9e-11_r8
      rate(:,247) = 3.3e-12_r8
      rate(:,266) = 4.5e-11_r8
      rate(:,267) = 4.62e-10_r8
      rate(:,268) = 1.2e-10_r8
      rate(:,269) = 9e-11_r8
      rate(:,270) = 3e-11_r8
      rate(:,275) = 2.14e-11_r8
      rate(:,276) = 1.9e-10_r8
      rate(:,289) = 2.57e-10_r8
      rate(:,290) = 1.8e-10_r8
      rate(:,291) = 1.794e-10_r8
      rate(:,292) = 1.3e-10_r8
      rate(:,293) = 7.65e-11_r8
      rate(:,307) = 4e-13_r8
      rate(:,311) = 1.31e-10_r8
      rate(:,312) = 3.5e-11_r8
      rate(:,313) = 9e-12_r8
      rate(:,320) = 6.8e-14_r8
      rate(:,321) = 2e-13_r8
      rate(:,335) = 7e-13_r8
      rate(:,336) = 1e-12_r8
      rate(:,340) = 1e-14_r8
      rate(:,341) = 1e-11_r8
      rate(:,342) = 1.15e-11_r8
      rate(:,343) = 4e-14_r8
      rate(:,356) = 3e-12_r8
      rate(:,357) = 6.7e-13_r8
      rate(:,367) = 3.5e-13_r8
      rate(:,368) = 5.4e-11_r8
      rate(:,371) = 2e-12_r8
      rate(:,372) = 1.4e-11_r8
      rate(:,375) = 2.4e-12_r8
      rate(:,386) = 5e-12_r8
      rate(:,396) = 1.6e-12_r8
      rate(:,398) = 6.7e-12_r8
      rate(:,401) = 3.5e-12_r8
      rate(:,404) = 1.3e-11_r8
      rate(:,405) = 1.4e-11_r8
      rate(:,409) = 2.4e-12_r8
      rate(:,410) = 1.4e-11_r8
      rate(:,415) = 2.4e-12_r8
      rate(:,416) = 4e-11_r8
      rate(:,417) = 4e-11_r8
      rate(:,419) = 1.4e-11_r8
      rate(:,423) = 2.4e-12_r8
      rate(:,424) = 4e-11_r8
      rate(:,428) = 7e-11_r8
      rate(:,429) = 1e-10_r8
      rate(:,434) = 2.4e-12_r8
      rate(:,449) = 4.7e-11_r8
      rate(:,462) = 2.1e-12_r8
      rate(:,463) = 2.8e-13_r8
      rate(:,471) = 1.7e-11_r8
      rate(:,477) = 8.4e-11_r8
      rate(:,479) = 1.9e-11_r8
      rate(:,480) = 1.2e-14_r8
      rate(:,481) = 2e-10_r8
      rate(:,488) = 2.4e-12_r8
      rate(:,489) = 2e-11_r8
      rate(:,493) = 2.3e-11_r8
      rate(:,494) = 2e-11_r8
      rate(:,498) = 3.3e-11_r8
      rate(:,499) = 1e-12_r8
      rate(:,500) = 5.7e-11_r8
      rate(:,501) = 3.4e-11_r8
      rate(:,504) = 2.3e-12_r8
      rate(:,505) = 1.2e-11_r8
      rate(:,506) = 5.7e-11_r8
      rate(:,507) = 2.8e-11_r8
      rate(:,508) = 6.6e-11_r8
      rate(:,509) = 1.4e-11_r8
      rate(:,512) = 1.9e-12_r8
      rate(:,528) = 6.34e-08_r8
      rate(:,534) = 1.9e-11_r8
      rate(:,535) = 1.2e-14_r8
      rate(:,536) = 2e-10_r8
      rate(:,541) = 1.34e-11_r8
      rate(:,545) = 1.34e-11_r8
      rate(:,547) = 1.7e-11_r8
      rate(:,568) = 6e-11_r8
      rate(:,571) = 1e-12_r8
      rate(:,572) = 4e-10_r8
      rate(:,573) = 2e-10_r8
      rate(:,574) = 1e-10_r8
      rate(:,575) = 5e-16_r8
      rate(:,576) = 4.4e-10_r8
      rate(:,577) = 9e-10_r8
      rate(:,580) = 1.29e-07_r8
      rate(:,581) = 2.31e-07_r8
      rate(:,582) = 2.31e-06_r8
      rate(:,583) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,154) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,155) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,156) = 2.64e-11_r8 * exp_fac(:)
      rate(:,157) = 6.6e-12_r8 * exp_fac(:)
      rate(:,161) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,163) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,166) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,167) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,170) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,171) = 1.4e-12_r8 * exp_fac(:)
      rate(:,425) = 1.05e-14_r8 * exp_fac(:)
      rate(:,539) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,176) = 3e-11_r8 * exp_fac(:)
      rate(:,264) = 5.5e-12_r8 * exp_fac(:)
      rate(:,303) = 3.8e-12_r8 * exp_fac(:)
      rate(:,325) = 3.8e-12_r8 * exp_fac(:)
      rate(:,352) = 3.8e-12_r8 * exp_fac(:)
      rate(:,360) = 3.8e-12_r8 * exp_fac(:)
      rate(:,364) = 3.8e-12_r8 * exp_fac(:)
      rate(:,380) = 2.3e-11_r8 * exp_fac(:)
      rate(:,390) = 3.8e-12_r8 * exp_fac(:)
      rate(:,400) = 3.8e-12_r8 * exp_fac(:)
      rate(:,427) = 1.52e-11_r8 * exp_fac(:)
      rate(:,435) = 1.52e-12_r8 * exp_fac(:)
      rate(:,441) = 3.8e-12_r8 * exp_fac(:)
      rate(:,444) = 3.8e-12_r8 * exp_fac(:)
      rate(:,448) = 3.8e-12_r8 * exp_fac(:)
      rate(:,464) = 3.8e-12_r8 * exp_fac(:)
      rate(:,468) = 3.8e-12_r8 * exp_fac(:)
      rate(:,474) = 3.8e-12_r8 * exp_fac(:)
      rate(:,478) = 3.8e-12_r8 * exp_fac(:)
      rate(:,177) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,178) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,179) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,181) = 4.8e-11_r8 * exp_fac(:)
      rate(:,262) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,182) = 1.8e-11_r8 * exp_fac(:)
      rate(:,338) = 4.2e-12_r8 * exp_fac(:)
      rate(:,351) = 4.2e-12_r8 * exp_fac(:)
      rate(:,359) = 4.2e-12_r8 * exp_fac(:)
      rate(:,388) = 4.2e-12_r8 * exp_fac(:)
      rate(:,408) = 4.4e-12_r8 * exp_fac(:)
      rate(:,414) = 4.4e-12_r8 * exp_fac(:)
      rate(:,487) = 4.2e-12_r8 * exp_fac(:)
      rate(:,492) = 4.2e-12_r8 * exp_fac(:)
      rate(:,497) = 4.2e-12_r8 * exp_fac(:)
      rate(:,183) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,187) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,190) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,191) = 2.9e-12_r8 * exp_fac(:)
      rate(:,192) = 1.45e-12_r8 * exp_fac(:)
      rate(:,193) = 1.45e-12_r8 * exp_fac(:)
      rate(:,194) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,195) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,196) = 1.2e-13_r8 * exp_fac(:)
      rate(:,222) = 3e-11_r8 * exp_fac(:)
      rate(:,199) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,203) = 3.3e-12_r8 * exp_fac(:)
      rate(:,218) = 1.4e-11_r8 * exp_fac(:)
      rate(:,232) = 7.4e-12_r8 * exp_fac(:)
      rate(:,334) = 8.1e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,204) = 3e-12_r8 * exp_fac(:)
      rate(:,263) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,206) = 7.26e-11_r8 * exp_fac(:)
      rate(:,207) = 4.64e-11_r8 * exp_fac(:)
      rate(:,214) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,215) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,216) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,217) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,219) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,220) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,221) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,223) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,224) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,225) = 2.6e-12_r8 * exp_fac(:)
      rate(:,226) = 6.4e-12_r8 * exp_fac(:)
      rate(:,256) = 4.1e-13_r8 * exp_fac(:)
      rate(:,437) = 7.5e-12_r8 * exp_fac(:)
      rate(:,451) = 7.5e-12_r8 * exp_fac(:)
      rate(:,454) = 7.5e-12_r8 * exp_fac(:)
      rate(:,457) = 7.5e-12_r8 * exp_fac(:)
      rate(:,227) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,229) = 3.6e-12_r8 * exp_fac(:)
      rate(:,278) = 2e-12_r8 * exp_fac(:)
      rate(:,230) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,231) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,233) = 6e-13_r8 * exp_fac(:)
      rate(:,253) = 1.5e-12_r8 * exp_fac(:)
      rate(:,261) = 1.9e-11_r8 * exp_fac(:)
      rate(:,234) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,235) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,236) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,238) = 3e-12_r8 * exp_fac(:)
      rate(:,272) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,250) = 1.7e-11_r8 * exp_fac(:)
      rate(:,277) = 6.3e-12_r8 * exp_fac(:)
      rate(:,251) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,252) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,254) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,255) = 2.3e-12_r8 * exp_fac(:)
      rate(:,258) = 8.8e-12_r8 * exp_fac(:)
      rate(:,257) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,260) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,265) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,271) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,273) = 1.4e-11_r8 * exp_fac(:)
      rate(:,275) = 2.14e-11_r8 * exp_fac(:)
      rate(:,276) = 1.9e-10_r8 * exp_fac(:)
      rate(:,289) = 2.57e-10_r8 * exp_fac(:)
      rate(:,290) = 1.8e-10_r8 * exp_fac(:)
      rate(:,291) = 1.794e-10_r8 * exp_fac(:)
      rate(:,292) = 1.3e-10_r8 * exp_fac(:)
      rate(:,293) = 7.65e-11_r8 * exp_fac(:)
      rate(:,307) = 4e-13_r8 * exp_fac(:)
      rate(:,311) = 1.31e-10_r8 * exp_fac(:)
      rate(:,312) = 3.5e-11_r8 * exp_fac(:)
      rate(:,313) = 9e-12_r8 * exp_fac(:)
      rate(:,320) = 6.8e-14_r8 * exp_fac(:)
      rate(:,321) = 2e-13_r8 * exp_fac(:)
      rate(:,335) = 7e-13_r8 * exp_fac(:)
      rate(:,336) = 1e-12_r8 * exp_fac(:)
      rate(:,340) = 1e-14_r8 * exp_fac(:)
      rate(:,341) = 1e-11_r8 * exp_fac(:)
      rate(:,342) = 1.15e-11_r8 * exp_fac(:)
      rate(:,343) = 4e-14_r8 * exp_fac(:)
      rate(:,356) = 3e-12_r8 * exp_fac(:)
      rate(:,357) = 6.7e-13_r8 * exp_fac(:)
      rate(:,367) = 3.5e-13_r8 * exp_fac(:)
      rate(:,368) = 5.4e-11_r8 * exp_fac(:)
      rate(:,371) = 2e-12_r8 * exp_fac(:)
      rate(:,372) = 1.4e-11_r8 * exp_fac(:)
      rate(:,375) = 2.4e-12_r8 * exp_fac(:)
      rate(:,386) = 5e-12_r8 * exp_fac(:)
      rate(:,396) = 1.6e-12_r8 * exp_fac(:)
      rate(:,398) = 6.7e-12_r8 * exp_fac(:)
      rate(:,401) = 3.5e-12_r8 * exp_fac(:)
      rate(:,404) = 1.3e-11_r8 * exp_fac(:)
      rate(:,405) = 1.4e-11_r8 * exp_fac(:)
      rate(:,409) = 2.4e-12_r8 * exp_fac(:)
      rate(:,410) = 1.4e-11_r8 * exp_fac(:)
      rate(:,415) = 2.4e-12_r8 * exp_fac(:)
      rate(:,416) = 4e-11_r8 * exp_fac(:)
      rate(:,417) = 4e-11_r8 * exp_fac(:)
      rate(:,419) = 1.4e-11_r8 * exp_fac(:)
      rate(:,423) = 2.4e-12_r8 * exp_fac(:)
      rate(:,424) = 4e-11_r8 * exp_fac(:)
      rate(:,428) = 7e-11_r8 * exp_fac(:)
      rate(:,429) = 1e-10_r8 * exp_fac(:)
      rate(:,434) = 2.4e-12_r8 * exp_fac(:)
      rate(:,449) = 4.7e-11_r8 * exp_fac(:)
      rate(:,462) = 2.1e-12_r8 * exp_fac(:)
      rate(:,463) = 2.8e-13_r8 * exp_fac(:)
      rate(:,471) = 1.7e-11_r8 * exp_fac(:)
      rate(:,477) = 8.4e-11_r8 * exp_fac(:)
      rate(:,479) = 1.9e-11_r8 * exp_fac(:)
      rate(:,480) = 1.2e-14_r8 * exp_fac(:)
      rate(:,481) = 2e-10_r8 * exp_fac(:)
      rate(:,488) = 2.4e-12_r8 * exp_fac(:)
      rate(:,489) = 2e-11_r8 * exp_fac(:)
      rate(:,493) = 2.3e-11_r8 * exp_fac(:)
      rate(:,494) = 2e-11_r8 * exp_fac(:)
      rate(:,498) = 3.3e-11_r8 * exp_fac(:)
      rate(:,499) = 1e-12_r8 * exp_fac(:)
      rate(:,500) = 5.7e-11_r8 * exp_fac(:)
      rate(:,501) = 3.4e-11_r8 * exp_fac(:)
      rate(:,504) = 2.3e-12_r8 * exp_fac(:)
      rate(:,505) = 1.2e-11_r8 * exp_fac(:)
      rate(:,506) = 5.7e-11_r8 * exp_fac(:)
      rate(:,507) = 2.8e-11_r8 * exp_fac(:)
      rate(:,508) = 6.6e-11_r8 * exp_fac(:)
      rate(:,509) = 1.4e-11_r8 * exp_fac(:)
      rate(:,512) = 1.9e-12_r8 * exp_fac(:)
      rate(:,528) = 6.34e-08_r8 * exp_fac(:)
      rate(:,534) = 1.9e-11_r8 * exp_fac(:)
      rate(:,535) = 1.2e-14_r8 * exp_fac(:)
      rate(:,536) = 2e-10_r8 * exp_fac(:)
      rate(:,541) = 1.34e-11_r8 * exp_fac(:)
      rate(:,545) = 1.34e-11_r8 * exp_fac(:)
      rate(:,547) = 1.7e-11_r8 * exp_fac(:)
      rate(:,568) = 6e-11_r8 * exp_fac(:)
      rate(:,571) = 1e-12_r8 * exp_fac(:)
      rate(:,572) = 4e-10_r8 * exp_fac(:)
      rate(:,573) = 2e-10_r8 * exp_fac(:)
      rate(:,574) = 1e-10_r8 * exp_fac(:)
      rate(:,575) = 5e-16_r8 * exp_fac(:)
      rate(:,576) = 4.4e-10_r8 * exp_fac(:)
      rate(:,577) = 9e-10_r8 * exp_fac(:)
      rate(:,580) = 1.29e-07_r8 * exp_fac(:)
      rate(:,581) = 2.31e-07_r8 * exp_fac(:)
      rate(:,582) = 2.31e-06_r8 * exp_fac(:)
      rate(:,583) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,274) = 6e-12_r8 * exp_fac(:)
      rate(:,373) = 5e-13_r8 * exp_fac(:)
      rate(:,406) = 5e-13_r8 * exp_fac(:)
      rate(:,411) = 5e-13_r8 * exp_fac(:)
      rate(:,420) = 5e-13_r8 * exp_fac(:)
      rate(:,431) = 5e-13_r8 * exp_fac(:)
      rate(:,279) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,280) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,281) = 1.64e-12_r8 * exp_fac(:)
      rate(:,392) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,282) = 2.03e-11_r8 * exp_fac(:)
      rate(:,511) = 3.4e-12_r8 * exp_fac(:)
      rate(:,283) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,284) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,285) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,286) = 1.25e-12_r8 * exp_fac(:)
      rate(:,296) = 3.4e-11_r8 * exp_fac(:)
      rate(:,287) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,288) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,294) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,295) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,297) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,298) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,299) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,300) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,301) = 2.8e-12_r8 * exp_fac(:)
      rate(:,363) = 2.9e-12_r8 * exp_fac(:)
      rate(:,302) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,304) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,308) = 7.5e-13_r8 * exp_fac(:)
      rate(:,322) = 7.5e-13_r8 * exp_fac(:)
      rate(:,337) = 7.5e-13_r8 * exp_fac(:)
      rate(:,350) = 7.5e-13_r8 * exp_fac(:)
      rate(:,358) = 7.5e-13_r8 * exp_fac(:)
      rate(:,362) = 8.6e-13_r8 * exp_fac(:)
      rate(:,374) = 8e-13_r8 * exp_fac(:)
      rate(:,387) = 7.5e-13_r8 * exp_fac(:)
      rate(:,397) = 7.5e-13_r8 * exp_fac(:)
      rate(:,407) = 8e-13_r8 * exp_fac(:)
      rate(:,412) = 8e-13_r8 * exp_fac(:)
      rate(:,421) = 8e-13_r8 * exp_fac(:)
      rate(:,432) = 8e-13_r8 * exp_fac(:)
      rate(:,439) = 7.5e-13_r8 * exp_fac(:)
      rate(:,443) = 7.5e-13_r8 * exp_fac(:)
      rate(:,446) = 7.5e-13_r8 * exp_fac(:)
      rate(:,459) = 7.5e-13_r8 * exp_fac(:)
      rate(:,466) = 7.5e-13_r8 * exp_fac(:)
      rate(:,472) = 7.5e-13_r8 * exp_fac(:)
      rate(:,475) = 7.5e-13_r8 * exp_fac(:)
      rate(:,486) = 7.5e-13_r8 * exp_fac(:)
      rate(:,491) = 7.5e-13_r8 * exp_fac(:)
      rate(:,496) = 7.5e-13_r8 * exp_fac(:)
      rate(:,309) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,310) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,314) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,319) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,323) = 2.6e-12_r8 * exp_fac(:)
      rate(:,440) = 2.6e-12_r8 * exp_fac(:)
      rate(:,445) = 2.6e-12_r8 * exp_fac(:)
      rate(:,447) = 2.6e-12_r8 * exp_fac(:)
      rate(:,460) = 2.6e-12_r8 * exp_fac(:)
      rate(:,467) = 2.6e-12_r8 * exp_fac(:)
      rate(:,473) = 2.6e-12_r8 * exp_fac(:)
      rate(:,476) = 2.6e-12_r8 * exp_fac(:)
      rate(:,324) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,326) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,327) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,328) = 1.4e-12_r8 * exp_fac(:)
      rate(:,348) = 6.5e-15_r8 * exp_fac(:)
      rate(:,329) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      rate(:,330) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,331) = 2.9e-12_r8 * exp_fac(:)
      rate(:,332) = 2e-12_r8 * exp_fac(:)
      rate(:,361) = 7.1e-13_r8 * exp_fac(:)
      rate(:,382) = 2e-12_r8 * exp_fac(:)
      rate(:,485) = 2e-12_r8 * exp_fac(:)
      rate(:,490) = 2e-12_r8 * exp_fac(:)
      rate(:,495) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,333) = 4.3e-13_r8 * exp_fac(:)
      rate(:,383) = 4.3e-13_r8 * exp_fac(:)
      rate(:,436) = 4.3e-13_r8 * exp_fac(:)
      rate(:,450) = 4.3e-13_r8 * exp_fac(:)
      rate(:,453) = 4.3e-13_r8 * exp_fac(:)
      rate(:,456) = 4.3e-13_r8 * exp_fac(:)
      rate(:,339) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,347) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,349) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,353) = 8.7e-12_r8 * exp( -615._r8 * itemp(:) )
      rate(:,354) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,355) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,369) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,370) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,376) = 2.7e-12_r8 * exp_fac(:)
      rate(:,377) = 1.3e-13_r8 * exp_fac(:)
      rate(:,379) = 9.6e-12_r8 * exp_fac(:)
      rate(:,385) = 5.3e-12_r8 * exp_fac(:)
      rate(:,422) = 2.7e-12_r8 * exp_fac(:)
      rate(:,433) = 2.7e-12_r8 * exp_fac(:)
      rate(:,378) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,381) = 4.6e-12_r8 * exp_fac(:)
      rate(:,384) = 2.3e-12_r8 * exp_fac(:)
      rate(:,389) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,393) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,399) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,402) = 1.86e-11_r8 * exp_fac(:)
      rate(:,403) = 1.86e-11_r8 * exp_fac(:)
      rate(:,413) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,418) = 3.03e-12_r8 * exp_fac(:)
      rate(:,538) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,426) = 2.54e-11_r8 * exp_fac(:)
      rate(:,540) = 2.54e-11_r8 * exp_fac(:)
      rate(:,430) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,438) = 2.3e-12_r8 * exp_fac(:)
      rate(:,537) = 2.3e-12_r8 * exp_fac(:)
      rate(:,442) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,461) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,469) = 1.7e-12_r8 * exp_fac(:)
      rate(:,546) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,482) = 1.2e-12_r8 * exp_fac(:)
      rate(:,542) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,483) = 6.3e-16_r8 * exp_fac(:)
      rate(:,543) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,484) = 1.2e-11_r8 * exp_fac(:)
      rate(:,544) = 1.2e-11_r8 * exp_fac(:)
      rate(:,502) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,503) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,510) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,513) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,516) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,517) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )
      rate(:,518) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,175), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,185), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,197), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,205), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,208), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,209), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,210), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,228), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,248), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,259), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,305), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.28e-33_r8
      kinf(:) = 9.3e-15_r8 * itemp(:)**(-4.42_r8)
      call jpl( rate(:,306), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,316), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,317), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,318), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,344), m, 0.48_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,345), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,365), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,391), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,452), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,455), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,458), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,465), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,159) = 1e-20_r8
      rate(:n,160) = 1.3e-16_r8
      rate(:n,164) = 8e-14_r8
      rate(:n,165) = 3.9e-17_r8
      rate(:n,172) = 6.9e-12_r8
      rate(:n,188) = 7e-13_r8
      rate(:n,189) = 5e-12_r8
      rate(:n,568) = 6e-11_r8
      rate(:n,571) = 1e-12_r8
      rate(:n,572) = 4e-10_r8
      rate(:n,573) = 2e-10_r8
      rate(:n,574) = 1e-10_r8
      rate(:n,576) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,155) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,156) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,157) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,161) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,163) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,166) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,167) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,176) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,177) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,178) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,181) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,182) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,183) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,190) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,194) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,195) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,203) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,204) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,175) = wrk(:)























      end subroutine setrxt_hrates

      end module mo_setrxt
