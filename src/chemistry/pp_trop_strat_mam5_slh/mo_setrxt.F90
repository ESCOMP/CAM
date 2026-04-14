
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

      rate(:,153) = 1.2e-10_r8
      rate(:,157) = 1.2e-10_r8
      rate(:,158) = 1.2e-10_r8
      rate(:,164) = 6.9e-12_r8
      rate(:,165) = 7.2e-11_r8
      rate(:,166) = 1.6e-12_r8
      rate(:,172) = 1.8e-12_r8
      rate(:,176) = 1.8e-12_r8
      rate(:,188) = 3.5e-12_r8
      rate(:,190) = 1.3e-11_r8
      rate(:,191) = 2.2e-11_r8
      rate(:,192) = 5e-11_r8
      rate(:,229) = 1.7e-13_r8
      rate(:,231) = 2.607e-10_r8
      rate(:,232) = 9.75e-11_r8
      rate(:,233) = 2.07e-10_r8
      rate(:,234) = 2.088e-10_r8
      rate(:,235) = 1.17e-10_r8
      rate(:,236) = 4.644e-11_r8
      rate(:,237) = 1.204e-10_r8
      rate(:,238) = 9.9e-11_r8
      rate(:,239) = 3.3e-12_r8
      rate(:,262) = 4.5e-11_r8
      rate(:,263) = 4.62e-10_r8
      rate(:,264) = 1.2e-10_r8
      rate(:,265) = 9e-11_r8
      rate(:,266) = 3e-11_r8
      rate(:,268) = 2e-13_r8
      rate(:,269) = 1.5e-12_r8
      rate(:,270) = 1.25e-10_r8
      rate(:,271) = 1.8e-10_r8
      rate(:,272) = 1.44e-11_r8
      rate(:,278) = 1e-10_r8
      rate(:,282) = 2.49e-11_r8
      rate(:,291) = 9e-12_r8
      rate(:,292) = 1.4e-10_r8
      rate(:,293) = 3.6e-16_r8
      rate(:,294) = 1e-10_r8
      rate(:,310) = 2.14e-11_r8
      rate(:,311) = 1.9e-10_r8
      rate(:,314) = 1.3e-12_r8
      rate(:,332) = 1.2e-12_r8
      rate(:,333) = 8e-13_r8
      rate(:,337) = 2.3e-12_r8
      rate(:,343) = 2.57e-10_r8
      rate(:,344) = 1.8e-10_r8
      rate(:,345) = 1.794e-10_r8
      rate(:,346) = 1.3e-10_r8
      rate(:,347) = 7.65e-11_r8
      rate(:,360) = 4e-13_r8
      rate(:,364) = 1.31e-10_r8
      rate(:,365) = 3.5e-11_r8
      rate(:,366) = 9e-12_r8
      rate(:,373) = 6.8e-14_r8
      rate(:,374) = 2e-13_r8
      rate(:,389) = 1e-12_r8
      rate(:,393) = 1e-14_r8
      rate(:,394) = 1e-11_r8
      rate(:,395) = 1.15e-11_r8
      rate(:,396) = 4e-14_r8
      rate(:,409) = 1.45e-10_r8
      rate(:,410) = 3e-12_r8
      rate(:,411) = 6.7e-13_r8
      rate(:,421) = 3.5e-13_r8
      rate(:,422) = 5.4e-11_r8
      rate(:,425) = 2e-12_r8
      rate(:,426) = 1.4e-11_r8
      rate(:,429) = 2.4e-12_r8
      rate(:,440) = 5e-12_r8
      rate(:,450) = 2.2e-12_r8
      rate(:,452) = 6.7e-12_r8
      rate(:,455) = 3.5e-12_r8
      rate(:,458) = 1.3e-11_r8
      rate(:,459) = 1.4e-11_r8
      rate(:,463) = 2.4e-12_r8
      rate(:,464) = 1.4e-11_r8
      rate(:,469) = 2.4e-12_r8
      rate(:,470) = 4e-11_r8
      rate(:,471) = 4e-11_r8
      rate(:,473) = 1.4e-11_r8
      rate(:,477) = 2.4e-12_r8
      rate(:,478) = 4e-11_r8
      rate(:,482) = 7e-11_r8
      rate(:,483) = 1e-10_r8
      rate(:,488) = 2.4e-12_r8
      rate(:,503) = 4.7e-11_r8
      rate(:,516) = 2.1e-12_r8
      rate(:,517) = 2.8e-13_r8
      rate(:,525) = 1.7e-11_r8
      rate(:,531) = 8.4e-11_r8
      rate(:,533) = 1.9e-11_r8
      rate(:,534) = 1.2e-14_r8
      rate(:,535) = 2e-10_r8
      rate(:,542) = 2.4e-12_r8
      rate(:,543) = 2e-11_r8
      rate(:,547) = 2.3e-11_r8
      rate(:,548) = 2e-11_r8
      rate(:,552) = 3.3e-11_r8
      rate(:,553) = 1e-12_r8
      rate(:,554) = 5.7e-11_r8
      rate(:,555) = 3.4e-11_r8
      rate(:,557) = 3.3e-10_r8
      rate(:,564) = 2.3e-12_r8
      rate(:,566) = 1.2e-11_r8
      rate(:,567) = 5.7e-11_r8
      rate(:,568) = 2.8e-11_r8
      rate(:,569) = 6.6e-11_r8
      rate(:,570) = 1.4e-11_r8
      rate(:,573) = 1.9e-12_r8
      rate(:,604) = 6.34e-08_r8
      rate(:,626) = 1.9e-11_r8
      rate(:,629) = 1.2e-14_r8
      rate(:,630) = 2e-10_r8
      rate(:,641) = 1.34e-11_r8
      rate(:,647) = 1.34e-11_r8
      rate(:,652) = 1.7e-11_r8
      rate(:,697) = 1.29e-07_r8
      rate(:,698) = 2.31e-07_r8
      rate(:,699) = 2.31e-06_r8
      rate(:,700) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,154) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,155) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,156) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,159) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,162) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,163) = 1.4e-12_r8 * exp_fac(:)
      rate(:,479) = 1.05e-14_r8 * exp_fac(:)
      rate(:,637) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,168) = 3e-11_r8 * exp_fac(:)
      rate(:,260) = 5.5e-12_r8 * exp_fac(:)
      rate(:,357) = 3.8e-12_r8 * exp_fac(:)
      rate(:,378) = 3.8e-12_r8 * exp_fac(:)
      rate(:,405) = 3.8e-12_r8 * exp_fac(:)
      rate(:,414) = 3.8e-12_r8 * exp_fac(:)
      rate(:,418) = 3.8e-12_r8 * exp_fac(:)
      rate(:,434) = 2.3e-11_r8 * exp_fac(:)
      rate(:,444) = 3.8e-12_r8 * exp_fac(:)
      rate(:,454) = 3.8e-12_r8 * exp_fac(:)
      rate(:,481) = 1.52e-11_r8 * exp_fac(:)
      rate(:,489) = 1.52e-12_r8 * exp_fac(:)
      rate(:,495) = 3.8e-12_r8 * exp_fac(:)
      rate(:,498) = 3.8e-12_r8 * exp_fac(:)
      rate(:,502) = 3.8e-12_r8 * exp_fac(:)
      rate(:,518) = 3.8e-12_r8 * exp_fac(:)
      rate(:,522) = 3.8e-12_r8 * exp_fac(:)
      rate(:,528) = 3.8e-12_r8 * exp_fac(:)
      rate(:,532) = 3.8e-12_r8 * exp_fac(:)
      rate(:,169) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,170) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,171) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,173) = 4.8e-11_r8 * exp_fac(:)
      rate(:,258) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,174) = 1.8e-11_r8 * exp_fac(:)
      rate(:,391) = 4.2e-12_r8 * exp_fac(:)
      rate(:,413) = 4.2e-12_r8 * exp_fac(:)
      rate(:,442) = 4.2e-12_r8 * exp_fac(:)
      rate(:,462) = 4.4e-12_r8 * exp_fac(:)
      rate(:,468) = 4.4e-12_r8 * exp_fac(:)
      rate(:,541) = 4.2e-12_r8 * exp_fac(:)
      rate(:,546) = 4.2e-12_r8 * exp_fac(:)
      rate(:,551) = 4.2e-12_r8 * exp_fac(:)
      rate(:,175) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,179) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,180) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,181) = 2.9e-12_r8 * exp_fac(:)
      rate(:,182) = 1.45e-12_r8 * exp_fac(:)
      rate(:,183) = 1.45e-12_r8 * exp_fac(:)
      rate(:,184) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,185) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,186) = 1.2e-13_r8 * exp_fac(:)
      rate(:,214) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,189) = 1.7e-11_r8 * exp_fac(:)
      rate(:,351) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,193) = 3.44e-12_r8 * exp_fac(:)
      rate(:,249) = 2.3e-12_r8 * exp_fac(:)
      rate(:,252) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,194) = 3e-12_r8 * exp_fac(:)
      rate(:,259) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,196) = 7.26e-11_r8 * exp_fac(:)
      rate(:,197) = 4.64e-11_r8 * exp_fac(:)
      rate(:,204) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,205) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,206) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,207) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,208) = 1.4e-11_r8 * exp_fac(:)
      rate(:,224) = 7.4e-12_r8 * exp_fac(:)
      rate(:,387) = 8.1e-12_r8 * exp_fac(:)
      rate(:,209) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,211) = 2.4e-12_r8 * exp( -1250._r8 * itemp(:) )
      rate(:,212) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,213) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,215) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,216) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,217) = 2.6e-12_r8 * exp_fac(:)
      rate(:,218) = 6.4e-12_r8 * exp_fac(:)
      rate(:,250) = 4.1e-13_r8 * exp_fac(:)
      rate(:,491) = 7.5e-12_r8 * exp_fac(:)
      rate(:,505) = 7.5e-12_r8 * exp_fac(:)
      rate(:,508) = 7.5e-12_r8 * exp_fac(:)
      rate(:,511) = 7.5e-12_r8 * exp_fac(:)
      rate(:,219) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,221) = 3.6e-12_r8 * exp_fac(:)
      rate(:,317) = 2e-12_r8 * exp_fac(:)
      rate(:,222) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,223) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,225) = 6e-13_r8 * exp_fac(:)
      rate(:,247) = 1.5e-12_r8 * exp_fac(:)
      rate(:,257) = 1.9e-11_r8 * exp_fac(:)
      rate(:,226) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,227) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,228) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,230) = 3e-12_r8 * exp_fac(:)
      rate(:,307) = 1.4e-10_r8 * exp_fac(:)
      rate(:,242) = 2.1e-11_r8 * exp( 240._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,243) = 1.7e-11_r8 * exp_fac(:)
      rate(:,316) = 6.3e-12_r8 * exp_fac(:)
      rate(:,244) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,246) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,248) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,251) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,253) = 1.78e-11_r8 * exp_fac(:)
      rate(:,376) = 2.6e-12_r8 * exp_fac(:)
      rate(:,494) = 2.6e-12_r8 * exp_fac(:)
      rate(:,499) = 2.6e-12_r8 * exp_fac(:)
      rate(:,501) = 2.6e-12_r8 * exp_fac(:)
      rate(:,514) = 2.6e-12_r8 * exp_fac(:)
      rate(:,521) = 2.6e-12_r8 * exp_fac(:)
      rate(:,527) = 2.6e-12_r8 * exp_fac(:)
      rate(:,530) = 2.6e-12_r8 * exp_fac(:)
      rate(:,633) = 2.6e-12_r8 * exp_fac(:)
      rate(:,640) = 2.6e-12_r8 * exp_fac(:)
      rate(:,650) = 2.6e-12_r8 * exp_fac(:)
      rate(:,654) = 2.6e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 215._r8 * itemp(:) )
      rate(:,254) = 6.28e-11_r8 * exp_fac(:)
      rate(:,256) = 1.9e-11_r8 * exp_fac(:)
      rate(:,261) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,267) = 1.3e-12_r8 * exp( -1830._r8 * itemp(:) )
      rate(:,273) = 1.5e-11_r8 * exp( -1090._r8 * itemp(:) )
      rate(:,274) = 9.1e-11_r8 * exp( -146._r8 * itemp(:) )
      rate(:,275) = 4.7e-13_r8 * exp( -1670._r8 * itemp(:) )
      rate(:,277) = 1.008e+18_r8 * exp( -13670._r8 * itemp(:) )
      rate(:,279) = 8.4e-11_r8 * exp( -2620._r8 * itemp(:) )
      rate(:,281) = 2.1e-11_r8 * exp( -830._r8 * itemp(:) )
      exp_fac(:) = exp( 510._r8 * itemp(:) )
      rate(:,283) = 3e-12_r8 * exp_fac(:)
      rate(:,284) = 1.2e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 280._r8 * itemp(:) )
      rate(:,285) = 2.585e-12_r8 * exp_fac(:)
      rate(:,286) = 1.175e-12_r8 * exp_fac(:)
      rate(:,287) = 9.4e-13_r8 * exp_fac(:)
      rate(:,288) = 1.3e-11_r8 * exp( 570._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,289) = 7.15e-12_r8 * exp_fac(:)
      rate(:,355) = 2.8e-12_r8 * exp_fac(:)
      rate(:,417) = 2.9e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,295) = 1.6e-11_r8 * exp_fac(:)
      rate(:,538) = 1.2e-11_r8 * exp_fac(:)
      rate(:,646) = 1.2e-11_r8 * exp_fac(:)
      rate(:,296) = 1.1e-12_r8 * exp( 542._r8 * itemp(:) )
      rate(:,306) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,308) = 1.4e-11_r8 * exp_fac(:)
      rate(:,310) = 2.14e-11_r8 * exp_fac(:)
      rate(:,311) = 1.9e-10_r8 * exp_fac(:)
      rate(:,314) = 1.3e-12_r8 * exp_fac(:)
      rate(:,332) = 1.2e-12_r8 * exp_fac(:)
      rate(:,333) = 8e-13_r8 * exp_fac(:)
      rate(:,337) = 2.3e-12_r8 * exp_fac(:)
      rate(:,343) = 2.57e-10_r8 * exp_fac(:)
      rate(:,344) = 1.8e-10_r8 * exp_fac(:)
      rate(:,345) = 1.794e-10_r8 * exp_fac(:)
      rate(:,346) = 1.3e-10_r8 * exp_fac(:)
      rate(:,347) = 7.65e-11_r8 * exp_fac(:)
      rate(:,360) = 4e-13_r8 * exp_fac(:)
      rate(:,364) = 1.31e-10_r8 * exp_fac(:)
      rate(:,365) = 3.5e-11_r8 * exp_fac(:)
      rate(:,366) = 9e-12_r8 * exp_fac(:)
      rate(:,373) = 6.8e-14_r8 * exp_fac(:)
      rate(:,374) = 2e-13_r8 * exp_fac(:)
      rate(:,389) = 1e-12_r8 * exp_fac(:)
      rate(:,393) = 1e-14_r8 * exp_fac(:)
      rate(:,394) = 1e-11_r8 * exp_fac(:)
      rate(:,395) = 1.15e-11_r8 * exp_fac(:)
      rate(:,396) = 4e-14_r8 * exp_fac(:)
      rate(:,409) = 1.45e-10_r8 * exp_fac(:)
      rate(:,410) = 3e-12_r8 * exp_fac(:)
      rate(:,411) = 6.7e-13_r8 * exp_fac(:)
      rate(:,421) = 3.5e-13_r8 * exp_fac(:)
      rate(:,422) = 5.4e-11_r8 * exp_fac(:)
      rate(:,425) = 2e-12_r8 * exp_fac(:)
      rate(:,426) = 1.4e-11_r8 * exp_fac(:)
      rate(:,429) = 2.4e-12_r8 * exp_fac(:)
      rate(:,440) = 5e-12_r8 * exp_fac(:)
      rate(:,450) = 2.2e-12_r8 * exp_fac(:)
      rate(:,452) = 6.7e-12_r8 * exp_fac(:)
      rate(:,455) = 3.5e-12_r8 * exp_fac(:)
      rate(:,458) = 1.3e-11_r8 * exp_fac(:)
      rate(:,459) = 1.4e-11_r8 * exp_fac(:)
      rate(:,463) = 2.4e-12_r8 * exp_fac(:)
      rate(:,464) = 1.4e-11_r8 * exp_fac(:)
      rate(:,469) = 2.4e-12_r8 * exp_fac(:)
      rate(:,470) = 4e-11_r8 * exp_fac(:)
      rate(:,471) = 4e-11_r8 * exp_fac(:)
      rate(:,473) = 1.4e-11_r8 * exp_fac(:)
      rate(:,477) = 2.4e-12_r8 * exp_fac(:)
      rate(:,478) = 4e-11_r8 * exp_fac(:)
      rate(:,482) = 7e-11_r8 * exp_fac(:)
      rate(:,483) = 1e-10_r8 * exp_fac(:)
      rate(:,488) = 2.4e-12_r8 * exp_fac(:)
      rate(:,503) = 4.7e-11_r8 * exp_fac(:)
      rate(:,516) = 2.1e-12_r8 * exp_fac(:)
      rate(:,517) = 2.8e-13_r8 * exp_fac(:)
      rate(:,525) = 1.7e-11_r8 * exp_fac(:)
      rate(:,531) = 8.4e-11_r8 * exp_fac(:)
      rate(:,533) = 1.9e-11_r8 * exp_fac(:)
      rate(:,534) = 1.2e-14_r8 * exp_fac(:)
      rate(:,535) = 2e-10_r8 * exp_fac(:)
      rate(:,542) = 2.4e-12_r8 * exp_fac(:)
      rate(:,543) = 2e-11_r8 * exp_fac(:)
      rate(:,547) = 2.3e-11_r8 * exp_fac(:)
      rate(:,548) = 2e-11_r8 * exp_fac(:)
      rate(:,552) = 3.3e-11_r8 * exp_fac(:)
      rate(:,553) = 1e-12_r8 * exp_fac(:)
      rate(:,554) = 5.7e-11_r8 * exp_fac(:)
      rate(:,555) = 3.4e-11_r8 * exp_fac(:)
      rate(:,557) = 3.3e-10_r8 * exp_fac(:)
      rate(:,564) = 2.3e-12_r8 * exp_fac(:)
      rate(:,566) = 1.2e-11_r8 * exp_fac(:)
      rate(:,567) = 5.7e-11_r8 * exp_fac(:)
      rate(:,568) = 2.8e-11_r8 * exp_fac(:)
      rate(:,569) = 6.6e-11_r8 * exp_fac(:)
      rate(:,570) = 1.4e-11_r8 * exp_fac(:)
      rate(:,573) = 1.9e-12_r8 * exp_fac(:)
      rate(:,604) = 6.34e-08_r8 * exp_fac(:)
      rate(:,626) = 1.9e-11_r8 * exp_fac(:)
      rate(:,629) = 1.2e-14_r8 * exp_fac(:)
      rate(:,630) = 2e-10_r8 * exp_fac(:)
      rate(:,641) = 1.34e-11_r8 * exp_fac(:)
      rate(:,647) = 1.34e-11_r8 * exp_fac(:)
      rate(:,652) = 1.7e-11_r8 * exp_fac(:)
      rate(:,697) = 1.29e-07_r8 * exp_fac(:)
      rate(:,698) = 2.31e-07_r8 * exp_fac(:)
      rate(:,699) = 2.31e-06_r8 * exp_fac(:)
      rate(:,700) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,309) = 6e-12_r8 * exp_fac(:)
      rate(:,427) = 5e-13_r8 * exp_fac(:)
      rate(:,460) = 5e-13_r8 * exp_fac(:)
      rate(:,465) = 5e-13_r8 * exp_fac(:)
      rate(:,474) = 5e-13_r8 * exp_fac(:)
      rate(:,485) = 5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( -990._r8 * itemp(:) )
      rate(:,313) = 4.7e-12_r8 * exp_fac(:)
      rate(:,338) = 3.3e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1150._r8 * itemp(:) )
      rate(:,315) = 1.14e-11_r8 * exp_fac(:)
      rate(:,322) = 1.42e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -880._r8 * itemp(:) )
      rate(:,318) = 2.1e-12_r8 * exp_fac(:)
      rate(:,320) = 1.92e-12_r8 * exp_fac(:)
      rate(:,319) = 7.4e-12_r8 * exp( -910._r8 * itemp(:) )
      rate(:,321) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,323) = 1.64e-12_r8 * exp_fac(:)
      rate(:,446) = 8.5e-16_r8 * exp_fac(:)
      rate(:,324) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,325) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,326) = 2.9e-11_r8 * exp( -1000._r8 * itemp(:) )
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,327) = 2.9e-12_r8 * exp_fac(:)
      rate(:,572) = 3.4e-12_r8 * exp_fac(:)
      rate(:,328) = 9e-13_r8 * exp( -420._r8 * itemp(:) )
      rate(:,329) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,330) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      rate(:,331) = 9.4e-13_r8 * exp( -510._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,334) = 3.92e-13_r8 * exp_fac(:)
      rate(:,335) = 1.68e-13_r8 * exp_fac(:)
      rate(:,361) = 7.5e-13_r8 * exp_fac(:)
      rate(:,375) = 7.5e-13_r8 * exp_fac(:)
      rate(:,390) = 7.5e-13_r8 * exp_fac(:)
      rate(:,412) = 7.5e-13_r8 * exp_fac(:)
      rate(:,416) = 8.6e-13_r8 * exp_fac(:)
      rate(:,428) = 8e-13_r8 * exp_fac(:)
      rate(:,441) = 7.5e-13_r8 * exp_fac(:)
      rate(:,451) = 7.5e-13_r8 * exp_fac(:)
      rate(:,461) = 8e-13_r8 * exp_fac(:)
      rate(:,466) = 8e-13_r8 * exp_fac(:)
      rate(:,475) = 8e-13_r8 * exp_fac(:)
      rate(:,486) = 8e-13_r8 * exp_fac(:)
      rate(:,493) = 7.5e-13_r8 * exp_fac(:)
      rate(:,497) = 7.5e-13_r8 * exp_fac(:)
      rate(:,500) = 7.5e-13_r8 * exp_fac(:)
      rate(:,513) = 7.5e-13_r8 * exp_fac(:)
      rate(:,520) = 7.5e-13_r8 * exp_fac(:)
      rate(:,526) = 7.5e-13_r8 * exp_fac(:)
      rate(:,529) = 7.5e-13_r8 * exp_fac(:)
      rate(:,540) = 7.5e-13_r8 * exp_fac(:)
      rate(:,545) = 7.5e-13_r8 * exp_fac(:)
      rate(:,550) = 7.5e-13_r8 * exp_fac(:)
      rate(:,632) = 7.5e-13_r8 * exp_fac(:)
      rate(:,639) = 7.5e-13_r8 * exp_fac(:)
      rate(:,649) = 7.5e-13_r8 * exp_fac(:)
      rate(:,653) = 7.5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,336) = 4.05e-12_r8 * exp_fac(:)
      rate(:,404) = 2.7e-12_r8 * exp_fac(:)
      rate(:,430) = 2.7e-12_r8 * exp_fac(:)
      rate(:,431) = 1.3e-13_r8 * exp_fac(:)
      rate(:,433) = 9.6e-12_r8 * exp_fac(:)
      rate(:,439) = 5.3e-12_r8 * exp_fac(:)
      rate(:,476) = 2.7e-12_r8 * exp_fac(:)
      rate(:,487) = 2.7e-12_r8 * exp_fac(:)
      rate(:,628) = 2.7e-12_r8 * exp_fac(:)
      rate(:,644) = 2.7e-12_r8 * exp_fac(:)
      rate(:,339) = 2.2e-12_r8 * exp( -920._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,340) = 1.25e-12_r8 * exp_fac(:)
      rate(:,350) = 3.4e-11_r8 * exp_fac(:)
      rate(:,341) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,342) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,348) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,349) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,352) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,353) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,354) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,356) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,358) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,362) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,363) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,367) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,372) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      rate(:,377) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,379) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,380) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,381) = 1.4e-12_r8 * exp_fac(:)
      rate(:,401) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,382) = 4.63e-12_r8 * exp_fac(:)
      rate(:,636) = 2.7e-12_r8 * exp_fac(:)
      rate(:,383) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,384) = 2.9e-12_r8 * exp_fac(:)
      rate(:,385) = 2e-12_r8 * exp_fac(:)
      rate(:,415) = 7.1e-13_r8 * exp_fac(:)
      rate(:,436) = 2e-12_r8 * exp_fac(:)
      rate(:,539) = 2e-12_r8 * exp_fac(:)
      rate(:,544) = 2e-12_r8 * exp_fac(:)
      rate(:,549) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,386) = 4.3e-13_r8 * exp_fac(:)
      rate(:,437) = 4.3e-13_r8 * exp_fac(:)
      rate(:,490) = 4.3e-13_r8 * exp_fac(:)
      rate(:,504) = 4.3e-13_r8 * exp_fac(:)
      rate(:,507) = 4.3e-13_r8 * exp_fac(:)
      rate(:,510) = 4.3e-13_r8 * exp_fac(:)
      rate(:,388) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,392) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,400) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,402) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,403) = 1.41e-13_r8 * exp_fac(:)
      rate(:,627) = 2.75e-13_r8 * exp_fac(:)
      rate(:,635) = 2.12e-13_r8 * exp_fac(:)
      rate(:,643) = 2.6e-13_r8 * exp_fac(:)
      rate(:,406) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,407) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,408) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,423) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,424) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,432) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,435) = 4.6e-12_r8 * exp_fac(:)
      rate(:,438) = 2.3e-12_r8 * exp_fac(:)
      rate(:,443) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,447) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,453) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,456) = 1.86e-11_r8 * exp_fac(:)
      rate(:,457) = 1.86e-11_r8 * exp_fac(:)
      rate(:,467) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,472) = 3.03e-12_r8 * exp_fac(:)
      rate(:,634) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,480) = 2.54e-11_r8 * exp_fac(:)
      rate(:,638) = 2.54e-11_r8 * exp_fac(:)
      rate(:,484) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,492) = 2.3e-12_r8 * exp_fac(:)
      rate(:,631) = 2.3e-12_r8 * exp_fac(:)
      rate(:,496) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,515) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,523) = 1.7e-12_r8 * exp_fac(:)
      rate(:,648) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,536) = 1.2e-12_r8 * exp_fac(:)
      rate(:,642) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,537) = 6.3e-16_r8 * exp_fac(:)
      rate(:,645) = 6.3e-16_r8 * exp_fac(:)
      rate(:,556) = 1e-14_r8 * exp( 950._r8 * itemp(:) )
      rate(:,558) = 9.4e-11_r8 * exp( 190._r8 * itemp(:) )
      rate(:,559) = 3.2e-13_r8 * exp( -925._r8 * itemp(:) )
      rate(:,560) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,561) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,562) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,563) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,571) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,574) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,594) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,167), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,177), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,187), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,195), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,198), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,199), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,200), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**2._r8
      kinf(:) = 1e-10_r8 * itemp(:)
      call jpl( rate(:,210), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,220), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,240), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.2e-31_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.7e-11_r8
      call jpl( rate(:,245), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,255), m, 0.6_r8, ko, kinf, n )

      ko(:) = 3e-31_r8 * itemp(:)**1._r8
      kinf(:) = 6.6e-11_r8
      call jpl( rate(:,276), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-32_r8 * itemp(:)**1._r8
      kinf(:) = 1.7e-11_r8
      call jpl( rate(:,280), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.5e-31_r8 * itemp(:)**3.5_r8
      kinf(:) = 7.6e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,290), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.4e-28_r8 * itemp(:)**8.5_r8
      kinf(:) = 4e-11_r8 * itemp(:)**1.2_r8
      call jpl( rate(:,312), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,359), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,369), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,370), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,371), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,397), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,398), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,419), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,445), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,448), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,506), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,509), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,512), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,519), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,565), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,164) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,155) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,159) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,168) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,169) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,170) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,173) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,174) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,175) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,180) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,184) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,185) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,193) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,194) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,167) = wrk(:)






























      end subroutine setrxt_hrates

      end module mo_setrxt
