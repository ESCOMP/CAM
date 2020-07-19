
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

      rate(:,168) = 1.2e-10_r8
      rate(:,172) = 1.2e-10_r8
      rate(:,178) = 6.9e-12_r8
      rate(:,179) = 7.2e-11_r8
      rate(:,180) = 1.6e-12_r8
      rate(:,186) = 1.8e-12_r8
      rate(:,190) = 1.8e-12_r8
      rate(:,202) = 3.5e-12_r8
      rate(:,204) = 1e-11_r8
      rate(:,205) = 2.2e-11_r8
      rate(:,206) = 5e-11_r8
      rate(:,241) = 1.7e-13_r8
      rate(:,243) = 2.607e-10_r8
      rate(:,244) = 9.75e-11_r8
      rate(:,245) = 2.07e-10_r8
      rate(:,246) = 2.088e-10_r8
      rate(:,247) = 1.17e-10_r8
      rate(:,248) = 4.644e-11_r8
      rate(:,249) = 1.204e-10_r8
      rate(:,250) = 9.9e-11_r8
      rate(:,251) = 3.3e-12_r8
      rate(:,270) = 4.5e-11_r8
      rate(:,271) = 4.62e-10_r8
      rate(:,272) = 1.2e-10_r8
      rate(:,273) = 9e-11_r8
      rate(:,274) = 3e-11_r8
      rate(:,279) = 2.14e-11_r8
      rate(:,280) = 1.9e-10_r8
      rate(:,293) = 2.57e-10_r8
      rate(:,294) = 1.8e-10_r8
      rate(:,295) = 1.794e-10_r8
      rate(:,296) = 1.3e-10_r8
      rate(:,297) = 7.65e-11_r8
      rate(:,311) = 4e-13_r8
      rate(:,316) = 1.31e-10_r8
      rate(:,317) = 3.5e-11_r8
      rate(:,318) = 9e-12_r8
      rate(:,325) = 6.8e-14_r8
      rate(:,326) = 2e-13_r8
      rate(:,340) = 7e-13_r8
      rate(:,341) = 1e-12_r8
      rate(:,345) = 1e-14_r8
      rate(:,346) = 1e-11_r8
      rate(:,347) = 1.15e-11_r8
      rate(:,348) = 3.3e-11_r8
      rate(:,349) = 3.4e-12_r8
      rate(:,350) = 4e-14_r8
      rate(:,363) = 3e-12_r8
      rate(:,364) = 1.2e-11_r8
      rate(:,365) = 6.7e-13_r8
      rate(:,375) = 3.5e-13_r8
      rate(:,376) = 5.4e-11_r8
      rate(:,377) = 3.77e-11_r8
      rate(:,380) = 2e-12_r8
      rate(:,381) = 1.29e-11_r8
      rate(:,383) = 4.5e-14_r8
      rate(:,388) = 3.77e-11_r8
      rate(:,394) = 4e-12_r8
      rate(:,400) = 1.78e-12_r8
      rate(:,402) = 6.1e-13_r8
      rate(:,406) = 4.8e-11_r8
      rate(:,409) = 1.6e-12_r8
      rate(:,411) = 6.7e-12_r8
      rate(:,414) = 3.5e-12_r8
      rate(:,419) = 6.42e-11_r8
      rate(:,426) = 1.6e-13_r8
      rate(:,432) = 1.4e-12_r8
      rate(:,437) = 7.5e-13_r8
      rate(:,438) = 1.4e-13_r8
      rate(:,439) = 7.5e-13_r8
      rate(:,440) = 3.6e-13_r8
      rate(:,441) = 6.5e-13_r8
      rate(:,442) = 2.1e-13_r8
      rate(:,443) = 6.5e-13_r8
      rate(:,444) = 4.9e-13_r8
      rate(:,446) = 1.2e-12_r8
      rate(:,450) = 9.8e-13_r8
      rate(:,453) = 1.85e-11_r8
      rate(:,454) = 1.63e-12_r8
      rate(:,455) = 2.5e-11_r8
      rate(:,456) = 1.1e-11_r8
      rate(:,457) = 3.3e-11_r8
      rate(:,460) = 2.8e-17_r8
      rate(:,461) = 8e-11_r8
      rate(:,464) = 3e-11_r8
      rate(:,467) = 4.2e-11_r8
      rate(:,470) = 2.8e-17_r8
      rate(:,471) = 1.1e-10_r8
      rate(:,473) = 3.9e-11_r8
      rate(:,476) = 1.3e-12_r8
      rate(:,478) = 5e-12_r8
      rate(:,479) = 2.3e-12_r8
      rate(:,482) = 3.9e-11_r8
      rate(:,485) = 2.8e-17_r8
      rate(:,486) = 9.2e-11_r8
      rate(:,489) = 3.85e-11_r8
      rate(:,493) = 1.2e-12_r8
      rate(:,497) = 9.8e-13_r8
      rate(:,502) = 4.4e-18_r8
      rate(:,503) = 3.6e-11_r8
      rate(:,555) = 4.7e-11_r8
      rate(:,568) = 2.1e-12_r8
      rate(:,569) = 2.8e-13_r8
      rate(:,577) = 1.7e-11_r8
      rate(:,583) = 8.4e-11_r8
      rate(:,586) = 5.3e-13_r8
      rate(:,588) = 2e-12_r8
      rate(:,591) = 2.3e-12_r8
      rate(:,596) = 2e-12_r8
      rate(:,599) = 2.3e-12_r8
      rate(:,605) = 1.9e-11_r8
      rate(:,606) = 5.3e-13_r8
      rate(:,608) = 2e-12_r8
      rate(:,611) = 2.3e-12_r8
      rate(:,616) = 2e-12_r8
      rate(:,619) = 2.3e-12_r8
      rate(:,623) = 1.2e-14_r8
      rate(:,624) = 2e-10_r8
      rate(:,625) = 2.5e-12_r8
      rate(:,626) = 5.3e-13_r8
      rate(:,628) = 2e-12_r8
      rate(:,631) = 2.3e-12_r8
      rate(:,636) = 2e-12_r8
      rate(:,639) = 2.3e-12_r8
      rate(:,645) = 1.2e-11_r8
      rate(:,647) = 2e-12_r8
      rate(:,649) = 5.3e-13_r8
      rate(:,651) = 2.3e-12_r8
      rate(:,656) = 2e-12_r8
      rate(:,659) = 2.3e-12_r8
      rate(:,665) = 1.1e-11_r8
      rate(:,667) = 2e-12_r8
      rate(:,669) = 5.3e-13_r8
      rate(:,671) = 2.3e-12_r8
      rate(:,676) = 2e-12_r8
      rate(:,679) = 2.3e-12_r8
      rate(:,684) = 2.1e-10_r8
      rate(:,690) = 8.9e-11_r8
      rate(:,691) = 8.9e-11_r8
      rate(:,695) = 2e-12_r8
      rate(:,698) = 2.3e-12_r8
      rate(:,706) = 4e-12_r8
      rate(:,709) = 2e-14_r8
      rate(:,711) = 2e-12_r8
      rate(:,714) = 2.3e-12_r8
      rate(:,719) = 2.52e-11_r8
      rate(:,724) = 4e-12_r8
      rate(:,728) = 2e-14_r8
      rate(:,730) = 2e-12_r8
      rate(:,733) = 2.3e-12_r8
      rate(:,738) = 1.92e-11_r8
      rate(:,740) = 2e-12_r8
      rate(:,743) = 2.3e-12_r8
      rate(:,747) = 8.8e-12_r8
      rate(:,748) = 8.8e-12_r8
      rate(:,749) = 8.8e-12_r8
      rate(:,754) = 4e-12_r8
      rate(:,756) = 2e-14_r8
      rate(:,758) = 3.66e-12_r8
      rate(:,759) = 2.8e-11_r8
      rate(:,760) = 2.6e-13_r8
      rate(:,763) = 8.3e-18_r8
      rate(:,764) = 1.1e-10_r8
      rate(:,768) = 1.1e-16_r8
      rate(:,770) = 3.64e-12_r8
      rate(:,771) = 2.8e-11_r8
      rate(:,772) = 1.7e-11_r8
      rate(:,775) = 1.1e-10_r8
      rate(:,776) = 9.58e-12_r8
      rate(:,779) = 1.1e-10_r8
      rate(:,780) = 1.23e-11_r8
      rate(:,783) = 1.1e-10_r8
      rate(:,784) = 3.64e-12_r8
      rate(:,787) = 1.1e-10_r8
      rate(:,788) = 5.5e-12_r8
      rate(:,789) = 4.65e-11_r8
      rate(:,790) = 2.8e-11_r8
      rate(:,798) = 2.3e-12_r8
      rate(:,799) = 1.2e-11_r8
      rate(:,800) = 5.7e-11_r8
      rate(:,801) = 2.8e-11_r8
      rate(:,802) = 6.6e-11_r8
      rate(:,803) = 1.4e-11_r8
      rate(:,806) = 1.9e-12_r8
      rate(:,830) = 6.34e-08_r8
      rate(:,847) = 1.9e-11_r8
      rate(:,850) = 1.2e-14_r8
      rate(:,851) = 2e-10_r8
      rate(:,855) = 2.5e-12_r8
      rate(:,867) = 1.34e-11_r8
      rate(:,868) = 1.2e-11_r8
      rate(:,873) = 1.1e-11_r8
      rate(:,877) = 2.1e-10_r8
      rate(:,878) = 1.34e-11_r8
      rate(:,882) = 1.7e-11_r8
      rate(:,902) = 1.29e-07_r8
      rate(:,903) = 2.31e-07_r8
      rate(:,904) = 2.31e-06_r8
      rate(:,905) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,169) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,170) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,171) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,173) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,176) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      rate(:,177) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,182) = 3e-11_r8 * exp_fac(:)
      rate(:,268) = 5.5e-12_r8 * exp_fac(:)
      rate(:,307) = 3.8e-12_r8 * exp_fac(:)
      rate(:,330) = 3.8e-12_r8 * exp_fac(:)
      rate(:,359) = 3.8e-12_r8 * exp_fac(:)
      rate(:,368) = 3.8e-12_r8 * exp_fac(:)
      rate(:,372) = 3.8e-12_r8 * exp_fac(:)
      rate(:,398) = 3.8e-12_r8 * exp_fac(:)
      rate(:,413) = 3.8e-12_r8 * exp_fac(:)
      rate(:,490) = 5.53e-12_r8 * exp_fac(:)
      rate(:,547) = 3.8e-12_r8 * exp_fac(:)
      rate(:,550) = 3.8e-12_r8 * exp_fac(:)
      rate(:,554) = 3.8e-12_r8 * exp_fac(:)
      rate(:,570) = 3.8e-12_r8 * exp_fac(:)
      rate(:,574) = 3.8e-12_r8 * exp_fac(:)
      rate(:,580) = 3.8e-12_r8 * exp_fac(:)
      rate(:,584) = 3.8e-12_r8 * exp_fac(:)
      rate(:,183) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,184) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,185) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,187) = 4.8e-11_r8 * exp_fac(:)
      rate(:,266) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,188) = 1.8e-11_r8 * exp_fac(:)
      rate(:,343) = 4.2e-12_r8 * exp_fac(:)
      rate(:,358) = 4.2e-12_r8 * exp_fac(:)
      rate(:,367) = 4.2e-12_r8 * exp_fac(:)
      rate(:,396) = 4.2e-12_r8 * exp_fac(:)
      rate(:,189) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,193) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,194) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,195) = 2.9e-12_r8 * exp_fac(:)
      rate(:,196) = 1.45e-12_r8 * exp_fac(:)
      rate(:,197) = 1.45e-12_r8 * exp_fac(:)
      rate(:,198) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,199) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,200) = 1.2e-13_r8 * exp_fac(:)
      rate(:,226) = 3e-11_r8 * exp_fac(:)
      rate(:,203) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,207) = 3.3e-12_r8 * exp_fac(:)
      rate(:,222) = 1.4e-11_r8 * exp_fac(:)
      rate(:,236) = 7.4e-12_r8 * exp_fac(:)
      rate(:,339) = 8.1e-12_r8 * exp_fac(:)
      rate(:,393) = 8.1e-12_r8 * exp_fac(:)
      rate(:,705) = 8.1e-12_r8 * exp_fac(:)
      rate(:,723) = 8.1e-12_r8 * exp_fac(:)
      rate(:,753) = 8.1e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,208) = 3e-12_r8 * exp_fac(:)
      rate(:,267) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,210) = 7.26e-11_r8 * exp_fac(:)
      rate(:,211) = 4.64e-11_r8 * exp_fac(:)
      rate(:,218) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      exp_fac(:) = exp( -1270._r8 * itemp(:) )
      rate(:,219) = 7.1e-12_r8 * exp_fac(:)
      rate(:,643) = 1.35e-15_r8 * exp_fac(:)
      rate(:,858) = 1.35e-15_r8 * exp_fac(:)
      rate(:,220) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,221) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,223) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,224) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,225) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,227) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,228) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,229) = 2.6e-12_r8 * exp_fac(:)
      rate(:,230) = 6.4e-12_r8 * exp_fac(:)
      rate(:,260) = 4.1e-13_r8 * exp_fac(:)
      rate(:,543) = 7.5e-12_r8 * exp_fac(:)
      rate(:,557) = 7.5e-12_r8 * exp_fac(:)
      rate(:,560) = 7.5e-12_r8 * exp_fac(:)
      rate(:,563) = 7.5e-12_r8 * exp_fac(:)
      rate(:,231) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,233) = 3.6e-12_r8 * exp_fac(:)
      rate(:,282) = 2e-12_r8 * exp_fac(:)
      rate(:,234) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,235) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,237) = 6e-13_r8 * exp_fac(:)
      rate(:,257) = 1.5e-12_r8 * exp_fac(:)
      rate(:,265) = 1.9e-11_r8 * exp_fac(:)
      rate(:,238) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,239) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,240) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,242) = 3e-12_r8 * exp_fac(:)
      rate(:,276) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,254) = 1.7e-11_r8 * exp_fac(:)
      rate(:,281) = 6.3e-12_r8 * exp_fac(:)
      rate(:,255) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,256) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,258) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,259) = 2.3e-12_r8 * exp_fac(:)
      rate(:,262) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 460._r8 * itemp(:) )
      rate(:,261) = 4.5e-12_r8 * exp_fac(:)
      rate(:,644) = 1.62e-11_r8 * exp_fac(:)
      rate(:,859) = 1.62e-11_r8 * exp_fac(:)
      rate(:,264) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,269) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,275) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,277) = 1.4e-11_r8 * exp_fac(:)
      rate(:,279) = 2.14e-11_r8 * exp_fac(:)
      rate(:,280) = 1.9e-10_r8 * exp_fac(:)
      rate(:,293) = 2.57e-10_r8 * exp_fac(:)
      rate(:,294) = 1.8e-10_r8 * exp_fac(:)
      rate(:,295) = 1.794e-10_r8 * exp_fac(:)
      rate(:,296) = 1.3e-10_r8 * exp_fac(:)
      rate(:,297) = 7.65e-11_r8 * exp_fac(:)
      rate(:,311) = 4e-13_r8 * exp_fac(:)
      rate(:,316) = 1.31e-10_r8 * exp_fac(:)
      rate(:,317) = 3.5e-11_r8 * exp_fac(:)
      rate(:,318) = 9e-12_r8 * exp_fac(:)
      rate(:,325) = 6.8e-14_r8 * exp_fac(:)
      rate(:,326) = 2e-13_r8 * exp_fac(:)
      rate(:,340) = 7e-13_r8 * exp_fac(:)
      rate(:,341) = 1e-12_r8 * exp_fac(:)
      rate(:,345) = 1e-14_r8 * exp_fac(:)
      rate(:,346) = 1e-11_r8 * exp_fac(:)
      rate(:,347) = 1.15e-11_r8 * exp_fac(:)
      rate(:,348) = 3.3e-11_r8 * exp_fac(:)
      rate(:,349) = 3.4e-12_r8 * exp_fac(:)
      rate(:,350) = 4e-14_r8 * exp_fac(:)
      rate(:,363) = 3e-12_r8 * exp_fac(:)
      rate(:,364) = 1.2e-11_r8 * exp_fac(:)
      rate(:,365) = 6.7e-13_r8 * exp_fac(:)
      rate(:,375) = 3.5e-13_r8 * exp_fac(:)
      rate(:,376) = 5.4e-11_r8 * exp_fac(:)
      rate(:,377) = 3.77e-11_r8 * exp_fac(:)
      rate(:,380) = 2e-12_r8 * exp_fac(:)
      rate(:,381) = 1.29e-11_r8 * exp_fac(:)
      rate(:,383) = 4.5e-14_r8 * exp_fac(:)
      rate(:,388) = 3.77e-11_r8 * exp_fac(:)
      rate(:,394) = 4e-12_r8 * exp_fac(:)
      rate(:,400) = 1.78e-12_r8 * exp_fac(:)
      rate(:,402) = 6.1e-13_r8 * exp_fac(:)
      rate(:,406) = 4.8e-11_r8 * exp_fac(:)
      rate(:,409) = 1.6e-12_r8 * exp_fac(:)
      rate(:,411) = 6.7e-12_r8 * exp_fac(:)
      rate(:,414) = 3.5e-12_r8 * exp_fac(:)
      rate(:,419) = 6.42e-11_r8 * exp_fac(:)
      rate(:,426) = 1.6e-13_r8 * exp_fac(:)
      rate(:,432) = 1.4e-12_r8 * exp_fac(:)
      rate(:,437) = 7.5e-13_r8 * exp_fac(:)
      rate(:,438) = 1.4e-13_r8 * exp_fac(:)
      rate(:,439) = 7.5e-13_r8 * exp_fac(:)
      rate(:,440) = 3.6e-13_r8 * exp_fac(:)
      rate(:,441) = 6.5e-13_r8 * exp_fac(:)
      rate(:,442) = 2.1e-13_r8 * exp_fac(:)
      rate(:,443) = 6.5e-13_r8 * exp_fac(:)
      rate(:,444) = 4.9e-13_r8 * exp_fac(:)
      rate(:,446) = 1.2e-12_r8 * exp_fac(:)
      rate(:,450) = 9.8e-13_r8 * exp_fac(:)
      rate(:,453) = 1.85e-11_r8 * exp_fac(:)
      rate(:,454) = 1.63e-12_r8 * exp_fac(:)
      rate(:,455) = 2.5e-11_r8 * exp_fac(:)
      rate(:,456) = 1.1e-11_r8 * exp_fac(:)
      rate(:,457) = 3.3e-11_r8 * exp_fac(:)
      rate(:,460) = 2.8e-17_r8 * exp_fac(:)
      rate(:,461) = 8e-11_r8 * exp_fac(:)
      rate(:,464) = 3e-11_r8 * exp_fac(:)
      rate(:,467) = 4.2e-11_r8 * exp_fac(:)
      rate(:,470) = 2.8e-17_r8 * exp_fac(:)
      rate(:,471) = 1.1e-10_r8 * exp_fac(:)
      rate(:,473) = 3.9e-11_r8 * exp_fac(:)
      rate(:,476) = 1.3e-12_r8 * exp_fac(:)
      rate(:,478) = 5e-12_r8 * exp_fac(:)
      rate(:,479) = 2.3e-12_r8 * exp_fac(:)
      rate(:,482) = 3.9e-11_r8 * exp_fac(:)
      rate(:,485) = 2.8e-17_r8 * exp_fac(:)
      rate(:,486) = 9.2e-11_r8 * exp_fac(:)
      rate(:,489) = 3.85e-11_r8 * exp_fac(:)
      rate(:,493) = 1.2e-12_r8 * exp_fac(:)
      rate(:,497) = 9.8e-13_r8 * exp_fac(:)
      rate(:,502) = 4.4e-18_r8 * exp_fac(:)
      rate(:,503) = 3.6e-11_r8 * exp_fac(:)
      rate(:,555) = 4.7e-11_r8 * exp_fac(:)
      rate(:,568) = 2.1e-12_r8 * exp_fac(:)
      rate(:,569) = 2.8e-13_r8 * exp_fac(:)
      rate(:,577) = 1.7e-11_r8 * exp_fac(:)
      rate(:,583) = 8.4e-11_r8 * exp_fac(:)
      rate(:,586) = 5.3e-13_r8 * exp_fac(:)
      rate(:,588) = 2e-12_r8 * exp_fac(:)
      rate(:,591) = 2.3e-12_r8 * exp_fac(:)
      rate(:,596) = 2e-12_r8 * exp_fac(:)
      rate(:,599) = 2.3e-12_r8 * exp_fac(:)
      rate(:,605) = 1.9e-11_r8 * exp_fac(:)
      rate(:,606) = 5.3e-13_r8 * exp_fac(:)
      rate(:,608) = 2e-12_r8 * exp_fac(:)
      rate(:,611) = 2.3e-12_r8 * exp_fac(:)
      rate(:,616) = 2e-12_r8 * exp_fac(:)
      rate(:,619) = 2.3e-12_r8 * exp_fac(:)
      rate(:,623) = 1.2e-14_r8 * exp_fac(:)
      rate(:,624) = 2e-10_r8 * exp_fac(:)
      rate(:,625) = 2.5e-12_r8 * exp_fac(:)
      rate(:,626) = 5.3e-13_r8 * exp_fac(:)
      rate(:,628) = 2e-12_r8 * exp_fac(:)
      rate(:,631) = 2.3e-12_r8 * exp_fac(:)
      rate(:,636) = 2e-12_r8 * exp_fac(:)
      rate(:,639) = 2.3e-12_r8 * exp_fac(:)
      rate(:,645) = 1.2e-11_r8 * exp_fac(:)
      rate(:,647) = 2e-12_r8 * exp_fac(:)
      rate(:,649) = 5.3e-13_r8 * exp_fac(:)
      rate(:,651) = 2.3e-12_r8 * exp_fac(:)
      rate(:,656) = 2e-12_r8 * exp_fac(:)
      rate(:,659) = 2.3e-12_r8 * exp_fac(:)
      rate(:,665) = 1.1e-11_r8 * exp_fac(:)
      rate(:,667) = 2e-12_r8 * exp_fac(:)
      rate(:,669) = 5.3e-13_r8 * exp_fac(:)
      rate(:,671) = 2.3e-12_r8 * exp_fac(:)
      rate(:,676) = 2e-12_r8 * exp_fac(:)
      rate(:,679) = 2.3e-12_r8 * exp_fac(:)
      rate(:,684) = 2.1e-10_r8 * exp_fac(:)
      rate(:,690) = 8.9e-11_r8 * exp_fac(:)
      rate(:,691) = 8.9e-11_r8 * exp_fac(:)
      rate(:,695) = 2e-12_r8 * exp_fac(:)
      rate(:,698) = 2.3e-12_r8 * exp_fac(:)
      rate(:,706) = 4e-12_r8 * exp_fac(:)
      rate(:,709) = 2e-14_r8 * exp_fac(:)
      rate(:,711) = 2e-12_r8 * exp_fac(:)
      rate(:,714) = 2.3e-12_r8 * exp_fac(:)
      rate(:,719) = 2.52e-11_r8 * exp_fac(:)
      rate(:,724) = 4e-12_r8 * exp_fac(:)
      rate(:,728) = 2e-14_r8 * exp_fac(:)
      rate(:,730) = 2e-12_r8 * exp_fac(:)
      rate(:,733) = 2.3e-12_r8 * exp_fac(:)
      rate(:,738) = 1.92e-11_r8 * exp_fac(:)
      rate(:,740) = 2e-12_r8 * exp_fac(:)
      rate(:,743) = 2.3e-12_r8 * exp_fac(:)
      rate(:,747) = 8.8e-12_r8 * exp_fac(:)
      rate(:,748) = 8.8e-12_r8 * exp_fac(:)
      rate(:,749) = 8.8e-12_r8 * exp_fac(:)
      rate(:,754) = 4e-12_r8 * exp_fac(:)
      rate(:,756) = 2e-14_r8 * exp_fac(:)
      rate(:,758) = 3.66e-12_r8 * exp_fac(:)
      rate(:,759) = 2.8e-11_r8 * exp_fac(:)
      rate(:,760) = 2.6e-13_r8 * exp_fac(:)
      rate(:,763) = 8.3e-18_r8 * exp_fac(:)
      rate(:,764) = 1.1e-10_r8 * exp_fac(:)
      rate(:,768) = 1.1e-16_r8 * exp_fac(:)
      rate(:,770) = 3.64e-12_r8 * exp_fac(:)
      rate(:,771) = 2.8e-11_r8 * exp_fac(:)
      rate(:,772) = 1.7e-11_r8 * exp_fac(:)
      rate(:,775) = 1.1e-10_r8 * exp_fac(:)
      rate(:,776) = 9.58e-12_r8 * exp_fac(:)
      rate(:,779) = 1.1e-10_r8 * exp_fac(:)
      rate(:,780) = 1.23e-11_r8 * exp_fac(:)
      rate(:,783) = 1.1e-10_r8 * exp_fac(:)
      rate(:,784) = 3.64e-12_r8 * exp_fac(:)
      rate(:,787) = 1.1e-10_r8 * exp_fac(:)
      rate(:,788) = 5.5e-12_r8 * exp_fac(:)
      rate(:,789) = 4.65e-11_r8 * exp_fac(:)
      rate(:,790) = 2.8e-11_r8 * exp_fac(:)
      rate(:,798) = 2.3e-12_r8 * exp_fac(:)
      rate(:,799) = 1.2e-11_r8 * exp_fac(:)
      rate(:,800) = 5.7e-11_r8 * exp_fac(:)
      rate(:,801) = 2.8e-11_r8 * exp_fac(:)
      rate(:,802) = 6.6e-11_r8 * exp_fac(:)
      rate(:,803) = 1.4e-11_r8 * exp_fac(:)
      rate(:,806) = 1.9e-12_r8 * exp_fac(:)
      rate(:,830) = 6.34e-08_r8 * exp_fac(:)
      rate(:,847) = 1.9e-11_r8 * exp_fac(:)
      rate(:,850) = 1.2e-14_r8 * exp_fac(:)
      rate(:,851) = 2e-10_r8 * exp_fac(:)
      rate(:,855) = 2.5e-12_r8 * exp_fac(:)
      rate(:,867) = 1.34e-11_r8 * exp_fac(:)
      rate(:,868) = 1.2e-11_r8 * exp_fac(:)
      rate(:,873) = 1.1e-11_r8 * exp_fac(:)
      rate(:,877) = 2.1e-10_r8 * exp_fac(:)
      rate(:,878) = 1.34e-11_r8 * exp_fac(:)
      rate(:,882) = 1.7e-11_r8 * exp_fac(:)
      rate(:,902) = 1.29e-07_r8 * exp_fac(:)
      rate(:,903) = 2.31e-07_r8 * exp_fac(:)
      rate(:,904) = 2.31e-06_r8 * exp_fac(:)
      rate(:,905) = 4.63e-07_r8 * exp_fac(:)
      rate(:,278) = 6e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,283) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,284) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,285) = 1.64e-12_r8 * exp_fac(:)
      rate(:,404) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,286) = 2.03e-11_r8 * exp_fac(:)
      rate(:,805) = 3.4e-12_r8 * exp_fac(:)
      rate(:,287) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,288) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,289) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,290) = 1.25e-12_r8 * exp_fac(:)
      rate(:,300) = 3.4e-11_r8 * exp_fac(:)
      rate(:,291) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,292) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,298) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,299) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,301) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,302) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,303) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,304) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,305) = 2.8e-12_r8 * exp_fac(:)
      rate(:,371) = 2.9e-12_r8 * exp_fac(:)
      rate(:,306) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,308) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,312) = 1.3e-12_r8 * exp_fac(:)
      rate(:,336) = 2.9e-12_r8 * exp_fac(:)
      rate(:,337) = 2e-12_r8 * exp_fac(:)
      rate(:,369) = 7.1e-13_r8 * exp_fac(:)
      rate(:,382) = 2e-12_r8 * exp_fac(:)
      rate(:,389) = 2.9e-12_r8 * exp_fac(:)
      rate(:,390) = 2e-12_r8 * exp_fac(:)
      rate(:,392) = 2.9e-12_r8 * exp_fac(:)
      rate(:,401) = 2e-12_r8 * exp_fac(:)
      rate(:,425) = 2e-12_r8 * exp_fac(:)
      rate(:,431) = 2e-12_r8 * exp_fac(:)
      rate(:,445) = 2e-12_r8 * exp_fac(:)
      rate(:,449) = 2e-12_r8 * exp_fac(:)
      rate(:,475) = 2e-12_r8 * exp_fac(:)
      rate(:,492) = 2e-12_r8 * exp_fac(:)
      rate(:,496) = 2e-12_r8 * exp_fac(:)
      rate(:,587) = 2e-12_r8 * exp_fac(:)
      rate(:,592) = 2e-12_r8 * exp_fac(:)
      rate(:,593) = 2e-12_r8 * exp_fac(:)
      rate(:,594) = 2e-12_r8 * exp_fac(:)
      rate(:,595) = 2e-12_r8 * exp_fac(:)
      rate(:,600) = 2e-12_r8 * exp_fac(:)
      rate(:,601) = 2e-12_r8 * exp_fac(:)
      rate(:,602) = 2e-12_r8 * exp_fac(:)
      rate(:,607) = 2e-12_r8 * exp_fac(:)
      rate(:,612) = 2e-12_r8 * exp_fac(:)
      rate(:,613) = 2e-12_r8 * exp_fac(:)
      rate(:,614) = 2e-12_r8 * exp_fac(:)
      rate(:,615) = 2e-12_r8 * exp_fac(:)
      rate(:,620) = 2e-12_r8 * exp_fac(:)
      rate(:,621) = 2e-12_r8 * exp_fac(:)
      rate(:,622) = 2e-12_r8 * exp_fac(:)
      rate(:,627) = 2e-12_r8 * exp_fac(:)
      rate(:,632) = 2e-12_r8 * exp_fac(:)
      rate(:,633) = 2e-12_r8 * exp_fac(:)
      rate(:,634) = 2e-12_r8 * exp_fac(:)
      rate(:,635) = 2e-12_r8 * exp_fac(:)
      rate(:,640) = 2e-12_r8 * exp_fac(:)
      rate(:,641) = 2e-12_r8 * exp_fac(:)
      rate(:,642) = 2e-12_r8 * exp_fac(:)
      rate(:,646) = 2e-12_r8 * exp_fac(:)
      rate(:,652) = 2e-12_r8 * exp_fac(:)
      rate(:,653) = 2e-12_r8 * exp_fac(:)
      rate(:,654) = 2e-12_r8 * exp_fac(:)
      rate(:,655) = 2e-12_r8 * exp_fac(:)
      rate(:,660) = 2e-12_r8 * exp_fac(:)
      rate(:,661) = 2e-12_r8 * exp_fac(:)
      rate(:,662) = 2e-12_r8 * exp_fac(:)
      rate(:,666) = 2e-12_r8 * exp_fac(:)
      rate(:,672) = 2e-12_r8 * exp_fac(:)
      rate(:,673) = 2e-12_r8 * exp_fac(:)
      rate(:,674) = 2e-12_r8 * exp_fac(:)
      rate(:,675) = 2e-12_r8 * exp_fac(:)
      rate(:,680) = 2e-12_r8 * exp_fac(:)
      rate(:,681) = 2e-12_r8 * exp_fac(:)
      rate(:,682) = 2e-12_r8 * exp_fac(:)
      rate(:,694) = 2e-12_r8 * exp_fac(:)
      rate(:,699) = 2e-12_r8 * exp_fac(:)
      rate(:,700) = 2e-12_r8 * exp_fac(:)
      rate(:,701) = 2e-12_r8 * exp_fac(:)
      rate(:,702) = 2.9e-12_r8 * exp_fac(:)
      rate(:,703) = 2e-12_r8 * exp_fac(:)
      rate(:,707) = 2.9e-12_r8 * exp_fac(:)
      rate(:,708) = 2.9e-12_r8 * exp_fac(:)
      rate(:,710) = 2e-12_r8 * exp_fac(:)
      rate(:,715) = 2e-12_r8 * exp_fac(:)
      rate(:,716) = 2e-12_r8 * exp_fac(:)
      rate(:,717) = 2e-12_r8 * exp_fac(:)
      rate(:,720) = 2.9e-12_r8 * exp_fac(:)
      rate(:,721) = 2e-12_r8 * exp_fac(:)
      rate(:,725) = 2.9e-12_r8 * exp_fac(:)
      rate(:,726) = 2.9e-12_r8 * exp_fac(:)
      rate(:,727) = 2.9e-12_r8 * exp_fac(:)
      rate(:,729) = 2e-12_r8 * exp_fac(:)
      rate(:,734) = 2e-12_r8 * exp_fac(:)
      rate(:,735) = 2e-12_r8 * exp_fac(:)
      rate(:,736) = 2e-12_r8 * exp_fac(:)
      rate(:,739) = 2e-12_r8 * exp_fac(:)
      rate(:,744) = 2e-12_r8 * exp_fac(:)
      rate(:,745) = 2e-12_r8 * exp_fac(:)
      rate(:,746) = 2e-12_r8 * exp_fac(:)
      rate(:,750) = 2.9e-12_r8 * exp_fac(:)
      rate(:,751) = 2e-12_r8 * exp_fac(:)
      rate(:,755) = 2.9e-12_r8 * exp_fac(:)
      rate(:,313) = 5.6e-15_r8 * exp( 2300._r8 * itemp(:) )
      rate(:,314) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,315) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,319) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,324) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,327) = 7.5e-13_r8 * exp_fac(:)
      rate(:,342) = 7.5e-13_r8 * exp_fac(:)
      rate(:,357) = 7.5e-13_r8 * exp_fac(:)
      rate(:,366) = 7.5e-13_r8 * exp_fac(:)
      rate(:,370) = 8.6e-13_r8 * exp_fac(:)
      rate(:,395) = 7.5e-13_r8 * exp_fac(:)
      rate(:,410) = 7.5e-13_r8 * exp_fac(:)
      rate(:,545) = 7.5e-13_r8 * exp_fac(:)
      rate(:,549) = 7.5e-13_r8 * exp_fac(:)
      rate(:,552) = 7.5e-13_r8 * exp_fac(:)
      rate(:,565) = 7.5e-13_r8 * exp_fac(:)
      rate(:,572) = 7.5e-13_r8 * exp_fac(:)
      rate(:,578) = 7.5e-13_r8 * exp_fac(:)
      rate(:,581) = 7.5e-13_r8 * exp_fac(:)
      rate(:,853) = 7.5e-13_r8 * exp_fac(:)
      rate(:,865) = 7.5e-13_r8 * exp_fac(:)
      rate(:,880) = 7.5e-13_r8 * exp_fac(:)
      rate(:,883) = 7.5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,328) = 2.6e-12_r8 * exp_fac(:)
      rate(:,546) = 2.6e-12_r8 * exp_fac(:)
      rate(:,551) = 2.6e-12_r8 * exp_fac(:)
      rate(:,553) = 2.6e-12_r8 * exp_fac(:)
      rate(:,566) = 2.6e-12_r8 * exp_fac(:)
      rate(:,573) = 2.6e-12_r8 * exp_fac(:)
      rate(:,579) = 2.6e-12_r8 * exp_fac(:)
      rate(:,582) = 2.6e-12_r8 * exp_fac(:)
      rate(:,854) = 2.6e-12_r8 * exp_fac(:)
      rate(:,866) = 2.6e-12_r8 * exp_fac(:)
      rate(:,881) = 2.6e-12_r8 * exp_fac(:)
      rate(:,884) = 2.6e-12_r8 * exp_fac(:)
      rate(:,329) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,331) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,332) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,333) = 1.4e-12_r8 * exp_fac(:)
      rate(:,355) = 6.5e-15_r8 * exp_fac(:)
      rate(:,334) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      rate(:,335) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,338) = 4.3e-13_r8 * exp_fac(:)
      rate(:,391) = 4.3e-13_r8 * exp_fac(:)
      rate(:,542) = 4.3e-13_r8 * exp_fac(:)
      rate(:,556) = 4.3e-13_r8 * exp_fac(:)
      rate(:,559) = 4.3e-13_r8 * exp_fac(:)
      rate(:,562) = 4.3e-13_r8 * exp_fac(:)
      rate(:,704) = 4.3e-13_r8 * exp_fac(:)
      rate(:,722) = 4.3e-13_r8 * exp_fac(:)
      rate(:,752) = 4.3e-13_r8 * exp_fac(:)
      rate(:,344) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,354) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,356) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,360) = 8.7e-12_r8 * exp( -615._r8 * itemp(:) )
      rate(:,361) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,362) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,378) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,379) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,384) = 2.11e-13_r8 * exp_fac(:)
      rate(:,403) = 2.11e-13_r8 * exp_fac(:)
      rate(:,422) = 2.38e-13_r8 * exp_fac(:)
      rate(:,427) = 2.12e-13_r8 * exp_fac(:)
      rate(:,433) = 2.12e-13_r8 * exp_fac(:)
      rate(:,447) = 2.12e-13_r8 * exp_fac(:)
      rate(:,451) = 2.12e-13_r8 * exp_fac(:)
      rate(:,458) = 2.6e-13_r8 * exp_fac(:)
      rate(:,462) = 2.6e-13_r8 * exp_fac(:)
      rate(:,465) = 2.6e-13_r8 * exp_fac(:)
      rate(:,468) = 2.6e-13_r8 * exp_fac(:)
      rate(:,472) = 2.6e-13_r8 * exp_fac(:)
      rate(:,477) = 2.47e-13_r8 * exp_fac(:)
      rate(:,480) = 2.64e-13_r8 * exp_fac(:)
      rate(:,483) = 2.64e-13_r8 * exp_fac(:)
      rate(:,494) = 2.12e-13_r8 * exp_fac(:)
      rate(:,498) = 2.12e-13_r8 * exp_fac(:)
      rate(:,500) = 2.6e-13_r8 * exp_fac(:)
      rate(:,589) = 2.71e-13_r8 * exp_fac(:)
      rate(:,597) = 2.6e-13_r8 * exp_fac(:)
      rate(:,609) = 2.78e-13_r8 * exp_fac(:)
      rate(:,617) = 2.75e-13_r8 * exp_fac(:)
      rate(:,629) = 2.71e-13_r8 * exp_fac(:)
      rate(:,637) = 2.6e-13_r8 * exp_fac(:)
      rate(:,648) = 2.71e-13_r8 * exp_fac(:)
      rate(:,657) = 2.6e-13_r8 * exp_fac(:)
      rate(:,668) = 2.71e-13_r8 * exp_fac(:)
      rate(:,677) = 2.6e-13_r8 * exp_fac(:)
      rate(:,688) = 2.71e-13_r8 * exp_fac(:)
      rate(:,692) = 2.71e-13_r8 * exp_fac(:)
      rate(:,696) = 2.54e-13_r8 * exp_fac(:)
      rate(:,712) = 2.62e-13_r8 * exp_fac(:)
      rate(:,731) = 2.66e-13_r8 * exp_fac(:)
      rate(:,741) = 2.51e-13_r8 * exp_fac(:)
      rate(:,761) = 2.68e-13_r8 * exp_fac(:)
      rate(:,766) = 2.47e-13_r8 * exp_fac(:)
      rate(:,773) = 2.76e-13_r8 * exp_fac(:)
      rate(:,777) = 2.76e-13_r8 * exp_fac(:)
      rate(:,781) = 2.75e-13_r8 * exp_fac(:)
      rate(:,785) = 2.75e-13_r8 * exp_fac(:)
      rate(:,843) = 2.6e-13_r8 * exp_fac(:)
      rate(:,848) = 2.75e-13_r8 * exp_fac(:)
      rate(:,856) = 2.6e-13_r8 * exp_fac(:)
      rate(:,861) = 2.12e-13_r8 * exp_fac(:)
      rate(:,869) = 2.6e-13_r8 * exp_fac(:)
      rate(:,874) = 2.6e-13_r8 * exp_fac(:)
      rate(:,385) = 2.9e+07_r8 * exp( -5297._r8 * itemp(:) )
      rate(:,386) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,387) = 9.6e-12_r8 * exp_fac(:)
      rate(:,590) = 2.7e-12_r8 * exp_fac(:)
      rate(:,598) = 2.7e-12_r8 * exp_fac(:)
      rate(:,610) = 2.7e-12_r8 * exp_fac(:)
      rate(:,618) = 2.7e-12_r8 * exp_fac(:)
      rate(:,630) = 2.7e-12_r8 * exp_fac(:)
      rate(:,638) = 2.7e-12_r8 * exp_fac(:)
      rate(:,650) = 2.7e-12_r8 * exp_fac(:)
      rate(:,658) = 2.7e-12_r8 * exp_fac(:)
      rate(:,670) = 2.7e-12_r8 * exp_fac(:)
      rate(:,678) = 2.7e-12_r8 * exp_fac(:)
      rate(:,689) = 2.7e-12_r8 * exp_fac(:)
      rate(:,693) = 2.7e-12_r8 * exp_fac(:)
      rate(:,697) = 2.7e-12_r8 * exp_fac(:)
      rate(:,713) = 2.7e-12_r8 * exp_fac(:)
      rate(:,732) = 2.7e-12_r8 * exp_fac(:)
      rate(:,742) = 2.7e-12_r8 * exp_fac(:)
      rate(:,762) = 2.7e-12_r8 * exp_fac(:)
      rate(:,767) = 2.7e-12_r8 * exp_fac(:)
      rate(:,774) = 2.7e-12_r8 * exp_fac(:)
      rate(:,778) = 2.7e-12_r8 * exp_fac(:)
      rate(:,782) = 2.7e-12_r8 * exp_fac(:)
      rate(:,786) = 2.7e-12_r8 * exp_fac(:)
      rate(:,844) = 2.7e-12_r8 * exp_fac(:)
      rate(:,849) = 2.7e-12_r8 * exp_fac(:)
      rate(:,857) = 2.7e-12_r8 * exp_fac(:)
      rate(:,862) = 2.7e-12_r8 * exp_fac(:)
      rate(:,870) = 2.7e-12_r8 * exp_fac(:)
      rate(:,875) = 2.7e-12_r8 * exp_fac(:)
      rate(:,397) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,405) = 2.7e-12_r8 * exp( 580._r8 * itemp(:) )
      rate(:,412) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 450._r8 * itemp(:) )
      rate(:,415) = 1.17e-11_r8 * exp_fac(:)
      rate(:,416) = 1.17e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 390._r8 * itemp(:) )
      rate(:,417) = 2.2e-11_r8 * exp_fac(:)
      rate(:,418) = 3.5e-11_r8 * exp_fac(:)
      rate(:,488) = 2.7e-11_r8 * exp_fac(:)
      rate(:,491) = 2.08e-11_r8 * exp_fac(:)
      rate(:,769) = 2.7e-11_r8 * exp_fac(:)
      rate(:,864) = 2.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,420) = 9.85e-12_r8 * exp_fac(:)
      rate(:,604) = 1.34e-11_r8 * exp_fac(:)
      rate(:,846) = 1.34e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -400._r8 * itemp(:) )
      rate(:,421) = 4.43e-11_r8 * exp_fac(:)
      rate(:,423) = 4.43e-11_r8 * exp_fac(:)
      rate(:,424) = 3.22e-11_r8 * exp_fac(:)
      rate(:,428) = 1.04e+11_r8 * exp( -9746._r8 * itemp(:) )
      rate(:,429) = 2.24e+15_r8 * exp( -10865._r8 * itemp(:) )
      rate(:,430) = 2.22e+15_r8 * exp( -10355._r8 * itemp(:) )
      rate(:,434) = 1.88e+11_r8 * exp( -9752._r8 * itemp(:) )
      rate(:,435) = 2.49e+15_r8 * exp( -11112._r8 * itemp(:) )
      rate(:,436) = 2.49e+15_r8 * exp( -10890._r8 * itemp(:) )
      rate(:,448) = 1.83e+14_r8 * exp( -8930._r8 * itemp(:) )
      rate(:,452) = 2.08e+14_r8 * exp( -9400._r8 * itemp(:) )
      exp_fac(:) = exp( -10000._r8 * itemp(:) )
      rate(:,459) = 1.256e+13_r8 * exp_fac(:)
      rate(:,463) = 1.875e+13_r8 * exp_fac(:)
      rate(:,466) = 1.875e+13_r8 * exp_fac(:)
      rate(:,469) = 5.092e+12_r8 * exp_fac(:)
      rate(:,481) = 8.72e+12_r8 * exp_fac(:)
      rate(:,484) = 6.55e+12_r8 * exp_fac(:)
      exp_fac(:) = exp( -450._r8 * itemp(:) )
      rate(:,474) = 2.95e-12_r8 * exp_fac(:)
      rate(:,765) = 2.95e-12_r8 * exp_fac(:)
      rate(:,860) = 2.95e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1995._r8 * itemp(:) )
      rate(:,487) = 1.03e-14_r8 * exp_fac(:)
      rate(:,863) = 1.03e-14_r8 * exp_fac(:)
      rate(:,495) = 1.79e+14_r8 * exp( -8830._r8 * itemp(:) )
      rate(:,499) = 1.75e+14_r8 * exp( -9054._r8 * itemp(:) )
      rate(:,501) = 1e+07_r8 * exp( -5000._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,544) = 2.3e-12_r8 * exp_fac(:)
      rate(:,852) = 2.3e-12_r8 * exp_fac(:)
      rate(:,548) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,567) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,575) = 1.7e-12_r8 * exp_fac(:)
      rate(:,879) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,585) = 1.2e-12_r8 * exp_fac(:)
      rate(:,842) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -640._r8 * itemp(:) )
      rate(:,603) = 8.05e-16_r8 * exp_fac(:)
      rate(:,845) = 8.05e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -770._r8 * itemp(:) )
      rate(:,663) = 2.8e-15_r8 * exp_fac(:)
      rate(:,871) = 2.8e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 470._r8 * itemp(:) )
      rate(:,664) = 3.41e-11_r8 * exp_fac(:)
      rate(:,872) = 3.41e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -520._r8 * itemp(:) )
      rate(:,683) = 2.65e-15_r8 * exp_fac(:)
      rate(:,876) = 2.65e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 600._r8 * itemp(:) )
      rate(:,718) = 5.2e-12_r8 * exp_fac(:)
      rate(:,737) = 5.2e-12_r8 * exp_fac(:)
      rate(:,757) = 5.2e-12_r8 * exp_fac(:)
      rate(:,794) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,795) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )
      rate(:,796) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,797) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,804) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,807) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,811) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,181), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,191), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,201), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,209), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,212), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,213), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,214), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,232), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,252), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,263), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,309), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.28e-33_r8
      kinf(:) = 9.3e-15_r8 * itemp(:)**(-4.42_r8)
      call jpl( rate(:,310), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,321), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,322), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,323), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,351), m, 0.48_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,352), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,373), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,399), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,407), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,558), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,561), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,564), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,571), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,685), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,686), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,687), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,178) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,170) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,173) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,182) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,183) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,184) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,187) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,188) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,189) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,194) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,198) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,199) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,207) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,208) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,181) = wrk(:)



























      end subroutine setrxt_hrates

      end module mo_setrxt
