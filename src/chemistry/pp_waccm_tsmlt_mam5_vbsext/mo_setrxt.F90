
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

      rate(:,171) = 0.000258_r8
      rate(:,172) = 0.085_r8
      rate(:,173) = 1.2e-10_r8
      rate(:,178) = 1.2e-10_r8
      rate(:,179) = 1e-20_r8
      rate(:,180) = 1.3e-16_r8
      rate(:,182) = 4.2e-13_r8
      rate(:,184) = 8e-14_r8
      rate(:,185) = 3.9e-17_r8
      rate(:,192) = 6.9e-12_r8
      rate(:,193) = 7.2e-11_r8
      rate(:,194) = 1.6e-12_r8
      rate(:,200) = 1.8e-12_r8
      rate(:,204) = 1.8e-12_r8
      rate(:,207) = 1.06e-05_r8
      rate(:,209) = 7e-11_r8
      rate(:,210) = 7e-13_r8
      rate(:,218) = 3.5e-12_r8
      rate(:,220) = 1.3e-11_r8
      rate(:,221) = 2.2e-11_r8
      rate(:,222) = 5e-11_r8
      rate(:,260) = 1.7e-13_r8
      rate(:,262) = 2.607e-10_r8
      rate(:,263) = 9.75e-11_r8
      rate(:,264) = 2.07e-10_r8
      rate(:,265) = 2.088e-10_r8
      rate(:,266) = 1.17e-10_r8
      rate(:,267) = 4.644e-11_r8
      rate(:,268) = 1.204e-10_r8
      rate(:,269) = 9.9e-11_r8
      rate(:,270) = 3.3e-12_r8
      rate(:,289) = 4.5e-11_r8
      rate(:,290) = 4.62e-10_r8
      rate(:,291) = 1.2e-10_r8
      rate(:,292) = 9e-11_r8
      rate(:,293) = 3e-11_r8
      rate(:,298) = 2.14e-11_r8
      rate(:,299) = 1.9e-10_r8
      rate(:,312) = 2.57e-10_r8
      rate(:,313) = 1.8e-10_r8
      rate(:,314) = 1.794e-10_r8
      rate(:,315) = 1.3e-10_r8
      rate(:,316) = 7.65e-11_r8
      rate(:,329) = 4e-13_r8
      rate(:,333) = 1.31e-10_r8
      rate(:,334) = 3.5e-11_r8
      rate(:,335) = 9e-12_r8
      rate(:,342) = 6.8e-14_r8
      rate(:,343) = 2e-13_r8
      rate(:,358) = 1e-12_r8
      rate(:,362) = 1e-14_r8
      rate(:,363) = 1e-11_r8
      rate(:,364) = 1.15e-11_r8
      rate(:,365) = 4e-14_r8
      rate(:,378) = 1.45e-10_r8
      rate(:,379) = 3e-12_r8
      rate(:,380) = 6.7e-13_r8
      rate(:,390) = 3.5e-13_r8
      rate(:,391) = 5.4e-11_r8
      rate(:,394) = 2e-12_r8
      rate(:,395) = 1.4e-11_r8
      rate(:,398) = 2.4e-12_r8
      rate(:,409) = 5e-12_r8
      rate(:,419) = 2.2e-12_r8
      rate(:,421) = 6.7e-12_r8
      rate(:,424) = 3.5e-12_r8
      rate(:,427) = 1.3e-11_r8
      rate(:,428) = 1.4e-11_r8
      rate(:,432) = 2.4e-12_r8
      rate(:,433) = 1.4e-11_r8
      rate(:,438) = 2.4e-12_r8
      rate(:,439) = 4e-11_r8
      rate(:,440) = 4e-11_r8
      rate(:,442) = 1.4e-11_r8
      rate(:,446) = 2.4e-12_r8
      rate(:,447) = 4e-11_r8
      rate(:,451) = 7e-11_r8
      rate(:,452) = 1e-10_r8
      rate(:,457) = 2.4e-12_r8
      rate(:,472) = 4.7e-11_r8
      rate(:,485) = 2.1e-12_r8
      rate(:,486) = 2.8e-13_r8
      rate(:,494) = 1.7e-11_r8
      rate(:,500) = 8.4e-11_r8
      rate(:,502) = 1.9e-11_r8
      rate(:,503) = 1.2e-14_r8
      rate(:,504) = 2e-10_r8
      rate(:,511) = 2.4e-12_r8
      rate(:,512) = 2e-11_r8
      rate(:,516) = 2.3e-11_r8
      rate(:,517) = 2e-11_r8
      rate(:,521) = 3.3e-11_r8
      rate(:,522) = 1e-12_r8
      rate(:,523) = 5.7e-11_r8
      rate(:,524) = 3.4e-11_r8
      rate(:,529) = 2.3e-12_r8
      rate(:,531) = 1.2e-11_r8
      rate(:,532) = 5.7e-11_r8
      rate(:,533) = 2.8e-11_r8
      rate(:,534) = 6.6e-11_r8
      rate(:,535) = 1.4e-11_r8
      rate(:,538) = 1.9e-12_r8
      rate(:,550) = 6.34e-08_r8
      rate(:,556) = 1.9e-11_r8
      rate(:,559) = 1.2e-14_r8
      rate(:,560) = 2e-10_r8
      rate(:,571) = 1.34e-11_r8
      rate(:,574) = 1.34e-11_r8
      rate(:,580) = 1.34e-11_r8
      rate(:,581) = 1.34e-11_r8
      rate(:,586) = 1.7e-11_r8
      rate(:,609) = 6e-11_r8
      rate(:,612) = 1e-12_r8
      rate(:,613) = 4e-10_r8
      rate(:,614) = 2e-10_r8
      rate(:,615) = 1e-10_r8
      rate(:,616) = 5e-16_r8
      rate(:,617) = 4.4e-10_r8
      rate(:,618) = 9e-10_r8
      rate(:,621) = 1.29e-07_r8
      rate(:,622) = 2.31e-07_r8
      rate(:,623) = 2.31e-06_r8
      rate(:,624) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,174) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,175) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,176) = 2.64e-11_r8 * exp_fac(:)
      rate(:,177) = 6.6e-12_r8 * exp_fac(:)
      rate(:,181) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,183) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,186) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,187) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,190) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,191) = 1.4e-12_r8 * exp_fac(:)
      rate(:,448) = 1.05e-14_r8 * exp_fac(:)
      rate(:,567) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,196) = 3e-11_r8 * exp_fac(:)
      rate(:,287) = 5.5e-12_r8 * exp_fac(:)
      rate(:,326) = 3.8e-12_r8 * exp_fac(:)
      rate(:,347) = 3.8e-12_r8 * exp_fac(:)
      rate(:,374) = 3.8e-12_r8 * exp_fac(:)
      rate(:,383) = 3.8e-12_r8 * exp_fac(:)
      rate(:,387) = 3.8e-12_r8 * exp_fac(:)
      rate(:,403) = 2.3e-11_r8 * exp_fac(:)
      rate(:,413) = 3.8e-12_r8 * exp_fac(:)
      rate(:,423) = 3.8e-12_r8 * exp_fac(:)
      rate(:,450) = 1.52e-11_r8 * exp_fac(:)
      rate(:,458) = 1.52e-12_r8 * exp_fac(:)
      rate(:,464) = 3.8e-12_r8 * exp_fac(:)
      rate(:,467) = 3.8e-12_r8 * exp_fac(:)
      rate(:,471) = 3.8e-12_r8 * exp_fac(:)
      rate(:,487) = 3.8e-12_r8 * exp_fac(:)
      rate(:,491) = 3.8e-12_r8 * exp_fac(:)
      rate(:,497) = 3.8e-12_r8 * exp_fac(:)
      rate(:,501) = 3.8e-12_r8 * exp_fac(:)
      rate(:,197) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,198) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,199) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,201) = 4.8e-11_r8 * exp_fac(:)
      rate(:,285) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,202) = 1.8e-11_r8 * exp_fac(:)
      rate(:,360) = 4.2e-12_r8 * exp_fac(:)
      rate(:,382) = 4.2e-12_r8 * exp_fac(:)
      rate(:,411) = 4.2e-12_r8 * exp_fac(:)
      rate(:,431) = 4.4e-12_r8 * exp_fac(:)
      rate(:,437) = 4.4e-12_r8 * exp_fac(:)
      rate(:,510) = 4.2e-12_r8 * exp_fac(:)
      rate(:,515) = 4.2e-12_r8 * exp_fac(:)
      rate(:,520) = 4.2e-12_r8 * exp_fac(:)
      rate(:,203) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,208) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,211) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,212) = 2.9e-12_r8 * exp_fac(:)
      rate(:,213) = 1.45e-12_r8 * exp_fac(:)
      rate(:,214) = 1.45e-12_r8 * exp_fac(:)
      rate(:,215) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,216) = 1.2e-13_r8 * exp_fac(:)
      rate(:,245) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,219) = 1.7e-11_r8 * exp_fac(:)
      rate(:,320) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,223) = 3.44e-12_r8 * exp_fac(:)
      rate(:,278) = 2.3e-12_r8 * exp_fac(:)
      rate(:,281) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,224) = 3e-12_r8 * exp_fac(:)
      rate(:,286) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,226) = 7.26e-11_r8 * exp_fac(:)
      rate(:,227) = 4.64e-11_r8 * exp_fac(:)
      rate(:,237) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,238) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,239) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,240) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,241) = 1.4e-11_r8 * exp_fac(:)
      rate(:,255) = 7.4e-12_r8 * exp_fac(:)
      rate(:,356) = 8.1e-12_r8 * exp_fac(:)
      rate(:,242) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,243) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,244) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,246) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,247) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,248) = 2.6e-12_r8 * exp_fac(:)
      rate(:,249) = 6.4e-12_r8 * exp_fac(:)
      rate(:,279) = 4.1e-13_r8 * exp_fac(:)
      rate(:,460) = 7.5e-12_r8 * exp_fac(:)
      rate(:,474) = 7.5e-12_r8 * exp_fac(:)
      rate(:,477) = 7.5e-12_r8 * exp_fac(:)
      rate(:,480) = 7.5e-12_r8 * exp_fac(:)
      rate(:,250) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,252) = 3.6e-12_r8 * exp_fac(:)
      rate(:,301) = 2e-12_r8 * exp_fac(:)
      rate(:,253) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,254) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,256) = 6e-13_r8 * exp_fac(:)
      rate(:,276) = 1.5e-12_r8 * exp_fac(:)
      rate(:,284) = 1.9e-11_r8 * exp_fac(:)
      rate(:,257) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,258) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,259) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,261) = 3e-12_r8 * exp_fac(:)
      rate(:,295) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,273) = 1.7e-11_r8 * exp_fac(:)
      rate(:,300) = 6.3e-12_r8 * exp_fac(:)
      rate(:,274) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,275) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,277) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,280) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,283) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,288) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,294) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,296) = 1.4e-11_r8 * exp_fac(:)
      rate(:,298) = 2.14e-11_r8 * exp_fac(:)
      rate(:,299) = 1.9e-10_r8 * exp_fac(:)
      rate(:,312) = 2.57e-10_r8 * exp_fac(:)
      rate(:,313) = 1.8e-10_r8 * exp_fac(:)
      rate(:,314) = 1.794e-10_r8 * exp_fac(:)
      rate(:,315) = 1.3e-10_r8 * exp_fac(:)
      rate(:,316) = 7.65e-11_r8 * exp_fac(:)
      rate(:,329) = 4e-13_r8 * exp_fac(:)
      rate(:,333) = 1.31e-10_r8 * exp_fac(:)
      rate(:,334) = 3.5e-11_r8 * exp_fac(:)
      rate(:,335) = 9e-12_r8 * exp_fac(:)
      rate(:,342) = 6.8e-14_r8 * exp_fac(:)
      rate(:,343) = 2e-13_r8 * exp_fac(:)
      rate(:,358) = 1e-12_r8 * exp_fac(:)
      rate(:,362) = 1e-14_r8 * exp_fac(:)
      rate(:,363) = 1e-11_r8 * exp_fac(:)
      rate(:,364) = 1.15e-11_r8 * exp_fac(:)
      rate(:,365) = 4e-14_r8 * exp_fac(:)
      rate(:,378) = 1.45e-10_r8 * exp_fac(:)
      rate(:,379) = 3e-12_r8 * exp_fac(:)
      rate(:,380) = 6.7e-13_r8 * exp_fac(:)
      rate(:,390) = 3.5e-13_r8 * exp_fac(:)
      rate(:,391) = 5.4e-11_r8 * exp_fac(:)
      rate(:,394) = 2e-12_r8 * exp_fac(:)
      rate(:,395) = 1.4e-11_r8 * exp_fac(:)
      rate(:,398) = 2.4e-12_r8 * exp_fac(:)
      rate(:,409) = 5e-12_r8 * exp_fac(:)
      rate(:,419) = 2.2e-12_r8 * exp_fac(:)
      rate(:,421) = 6.7e-12_r8 * exp_fac(:)
      rate(:,424) = 3.5e-12_r8 * exp_fac(:)
      rate(:,427) = 1.3e-11_r8 * exp_fac(:)
      rate(:,428) = 1.4e-11_r8 * exp_fac(:)
      rate(:,432) = 2.4e-12_r8 * exp_fac(:)
      rate(:,433) = 1.4e-11_r8 * exp_fac(:)
      rate(:,438) = 2.4e-12_r8 * exp_fac(:)
      rate(:,439) = 4e-11_r8 * exp_fac(:)
      rate(:,440) = 4e-11_r8 * exp_fac(:)
      rate(:,442) = 1.4e-11_r8 * exp_fac(:)
      rate(:,446) = 2.4e-12_r8 * exp_fac(:)
      rate(:,447) = 4e-11_r8 * exp_fac(:)
      rate(:,451) = 7e-11_r8 * exp_fac(:)
      rate(:,452) = 1e-10_r8 * exp_fac(:)
      rate(:,457) = 2.4e-12_r8 * exp_fac(:)
      rate(:,472) = 4.7e-11_r8 * exp_fac(:)
      rate(:,485) = 2.1e-12_r8 * exp_fac(:)
      rate(:,486) = 2.8e-13_r8 * exp_fac(:)
      rate(:,494) = 1.7e-11_r8 * exp_fac(:)
      rate(:,500) = 8.4e-11_r8 * exp_fac(:)
      rate(:,502) = 1.9e-11_r8 * exp_fac(:)
      rate(:,503) = 1.2e-14_r8 * exp_fac(:)
      rate(:,504) = 2e-10_r8 * exp_fac(:)
      rate(:,511) = 2.4e-12_r8 * exp_fac(:)
      rate(:,512) = 2e-11_r8 * exp_fac(:)
      rate(:,516) = 2.3e-11_r8 * exp_fac(:)
      rate(:,517) = 2e-11_r8 * exp_fac(:)
      rate(:,521) = 3.3e-11_r8 * exp_fac(:)
      rate(:,522) = 1e-12_r8 * exp_fac(:)
      rate(:,523) = 5.7e-11_r8 * exp_fac(:)
      rate(:,524) = 3.4e-11_r8 * exp_fac(:)
      rate(:,529) = 2.3e-12_r8 * exp_fac(:)
      rate(:,531) = 1.2e-11_r8 * exp_fac(:)
      rate(:,532) = 5.7e-11_r8 * exp_fac(:)
      rate(:,533) = 2.8e-11_r8 * exp_fac(:)
      rate(:,534) = 6.6e-11_r8 * exp_fac(:)
      rate(:,535) = 1.4e-11_r8 * exp_fac(:)
      rate(:,538) = 1.9e-12_r8 * exp_fac(:)
      rate(:,550) = 6.34e-08_r8 * exp_fac(:)
      rate(:,556) = 1.9e-11_r8 * exp_fac(:)
      rate(:,559) = 1.2e-14_r8 * exp_fac(:)
      rate(:,560) = 2e-10_r8 * exp_fac(:)
      rate(:,571) = 1.34e-11_r8 * exp_fac(:)
      rate(:,574) = 1.34e-11_r8 * exp_fac(:)
      rate(:,580) = 1.34e-11_r8 * exp_fac(:)
      rate(:,581) = 1.34e-11_r8 * exp_fac(:)
      rate(:,586) = 1.7e-11_r8 * exp_fac(:)
      rate(:,609) = 6e-11_r8 * exp_fac(:)
      rate(:,612) = 1e-12_r8 * exp_fac(:)
      rate(:,613) = 4e-10_r8 * exp_fac(:)
      rate(:,614) = 2e-10_r8 * exp_fac(:)
      rate(:,615) = 1e-10_r8 * exp_fac(:)
      rate(:,616) = 5e-16_r8 * exp_fac(:)
      rate(:,617) = 4.4e-10_r8 * exp_fac(:)
      rate(:,618) = 9e-10_r8 * exp_fac(:)
      rate(:,621) = 1.29e-07_r8 * exp_fac(:)
      rate(:,622) = 2.31e-07_r8 * exp_fac(:)
      rate(:,623) = 2.31e-06_r8 * exp_fac(:)
      rate(:,624) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,297) = 6e-12_r8 * exp_fac(:)
      rate(:,396) = 5e-13_r8 * exp_fac(:)
      rate(:,429) = 5e-13_r8 * exp_fac(:)
      rate(:,434) = 5e-13_r8 * exp_fac(:)
      rate(:,443) = 5e-13_r8 * exp_fac(:)
      rate(:,454) = 5e-13_r8 * exp_fac(:)
      rate(:,302) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,303) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,304) = 1.64e-12_r8 * exp_fac(:)
      rate(:,415) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,305) = 2.03e-11_r8 * exp_fac(:)
      rate(:,537) = 3.4e-12_r8 * exp_fac(:)
      rate(:,306) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,307) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,308) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,309) = 1.25e-12_r8 * exp_fac(:)
      rate(:,319) = 3.4e-11_r8 * exp_fac(:)
      rate(:,310) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,311) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,317) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,318) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,321) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,322) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,323) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,324) = 2.8e-12_r8 * exp_fac(:)
      rate(:,386) = 2.9e-12_r8 * exp_fac(:)
      rate(:,325) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,327) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,330) = 7.5e-13_r8 * exp_fac(:)
      rate(:,344) = 7.5e-13_r8 * exp_fac(:)
      rate(:,359) = 7.5e-13_r8 * exp_fac(:)
      rate(:,381) = 7.5e-13_r8 * exp_fac(:)
      rate(:,385) = 8.6e-13_r8 * exp_fac(:)
      rate(:,397) = 8e-13_r8 * exp_fac(:)
      rate(:,410) = 7.5e-13_r8 * exp_fac(:)
      rate(:,420) = 7.5e-13_r8 * exp_fac(:)
      rate(:,430) = 8e-13_r8 * exp_fac(:)
      rate(:,435) = 8e-13_r8 * exp_fac(:)
      rate(:,444) = 8e-13_r8 * exp_fac(:)
      rate(:,455) = 8e-13_r8 * exp_fac(:)
      rate(:,462) = 7.5e-13_r8 * exp_fac(:)
      rate(:,466) = 7.5e-13_r8 * exp_fac(:)
      rate(:,469) = 7.5e-13_r8 * exp_fac(:)
      rate(:,482) = 7.5e-13_r8 * exp_fac(:)
      rate(:,489) = 7.5e-13_r8 * exp_fac(:)
      rate(:,495) = 7.5e-13_r8 * exp_fac(:)
      rate(:,498) = 7.5e-13_r8 * exp_fac(:)
      rate(:,509) = 7.5e-13_r8 * exp_fac(:)
      rate(:,514) = 7.5e-13_r8 * exp_fac(:)
      rate(:,519) = 7.5e-13_r8 * exp_fac(:)
      rate(:,562) = 7.5e-13_r8 * exp_fac(:)
      rate(:,569) = 7.5e-13_r8 * exp_fac(:)
      rate(:,572) = 7.5e-13_r8 * exp_fac(:)
      rate(:,583) = 7.5e-13_r8 * exp_fac(:)
      rate(:,587) = 7.5e-13_r8 * exp_fac(:)
      rate(:,331) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,332) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,336) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,341) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,345) = 2.6e-12_r8 * exp_fac(:)
      rate(:,463) = 2.6e-12_r8 * exp_fac(:)
      rate(:,468) = 2.6e-12_r8 * exp_fac(:)
      rate(:,470) = 2.6e-12_r8 * exp_fac(:)
      rate(:,483) = 2.6e-12_r8 * exp_fac(:)
      rate(:,490) = 2.6e-12_r8 * exp_fac(:)
      rate(:,496) = 2.6e-12_r8 * exp_fac(:)
      rate(:,499) = 2.6e-12_r8 * exp_fac(:)
      rate(:,563) = 2.6e-12_r8 * exp_fac(:)
      rate(:,570) = 2.6e-12_r8 * exp_fac(:)
      rate(:,573) = 2.6e-12_r8 * exp_fac(:)
      rate(:,584) = 2.6e-12_r8 * exp_fac(:)
      rate(:,588) = 2.6e-12_r8 * exp_fac(:)
      rate(:,346) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,348) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,349) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,350) = 1.4e-12_r8 * exp_fac(:)
      rate(:,370) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,351) = 4.63e-12_r8 * exp_fac(:)
      rate(:,566) = 2.7e-12_r8 * exp_fac(:)
      rate(:,352) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,353) = 2.9e-12_r8 * exp_fac(:)
      rate(:,354) = 2e-12_r8 * exp_fac(:)
      rate(:,384) = 7.1e-13_r8 * exp_fac(:)
      rate(:,405) = 2e-12_r8 * exp_fac(:)
      rate(:,508) = 2e-12_r8 * exp_fac(:)
      rate(:,513) = 2e-12_r8 * exp_fac(:)
      rate(:,518) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,355) = 4.3e-13_r8 * exp_fac(:)
      rate(:,406) = 4.3e-13_r8 * exp_fac(:)
      rate(:,459) = 4.3e-13_r8 * exp_fac(:)
      rate(:,473) = 4.3e-13_r8 * exp_fac(:)
      rate(:,476) = 4.3e-13_r8 * exp_fac(:)
      rate(:,479) = 4.3e-13_r8 * exp_fac(:)
      rate(:,357) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,361) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,369) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,371) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,372) = 1.41e-13_r8 * exp_fac(:)
      rate(:,557) = 2.75e-13_r8 * exp_fac(:)
      rate(:,565) = 2.12e-13_r8 * exp_fac(:)
      rate(:,576) = 2.6e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,373) = 2.7e-12_r8 * exp_fac(:)
      rate(:,399) = 2.7e-12_r8 * exp_fac(:)
      rate(:,400) = 1.3e-13_r8 * exp_fac(:)
      rate(:,402) = 9.6e-12_r8 * exp_fac(:)
      rate(:,408) = 5.3e-12_r8 * exp_fac(:)
      rate(:,445) = 2.7e-12_r8 * exp_fac(:)
      rate(:,456) = 2.7e-12_r8 * exp_fac(:)
      rate(:,558) = 2.7e-12_r8 * exp_fac(:)
      rate(:,577) = 2.7e-12_r8 * exp_fac(:)
      rate(:,375) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,376) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,377) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,392) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,393) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,401) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,404) = 4.6e-12_r8 * exp_fac(:)
      rate(:,407) = 2.3e-12_r8 * exp_fac(:)
      rate(:,412) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,416) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,422) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,425) = 1.86e-11_r8 * exp_fac(:)
      rate(:,426) = 1.86e-11_r8 * exp_fac(:)
      rate(:,436) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,441) = 3.03e-12_r8 * exp_fac(:)
      rate(:,564) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,449) = 2.54e-11_r8 * exp_fac(:)
      rate(:,568) = 2.54e-11_r8 * exp_fac(:)
      rate(:,453) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,461) = 2.3e-12_r8 * exp_fac(:)
      rate(:,561) = 2.3e-12_r8 * exp_fac(:)
      rate(:,465) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,484) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,492) = 1.7e-12_r8 * exp_fac(:)
      rate(:,582) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,505) = 1.2e-12_r8 * exp_fac(:)
      rate(:,575) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,506) = 6.3e-16_r8 * exp_fac(:)
      rate(:,578) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,507) = 1.2e-11_r8 * exp_fac(:)
      rate(:,579) = 1.2e-11_r8 * exp_fac(:)
      rate(:,525) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,526) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,527) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,528) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,536) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,539) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,542) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,195), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,205), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,217), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,225), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,228), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,229), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,230), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,251), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,271), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,282), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,328), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,338), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,339), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,340), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,366), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,367), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,388), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,414), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,417), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,475), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,478), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,481), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,488), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,530), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,179) = 1e-20_r8
      rate(:n,180) = 1.3e-16_r8
      rate(:n,184) = 8e-14_r8
      rate(:n,185) = 3.9e-17_r8
      rate(:n,192) = 6.9e-12_r8
      rate(:n,209) = 7e-11_r8
      rate(:n,210) = 7e-13_r8
      rate(:n,609) = 6e-11_r8
      rate(:n,612) = 1e-12_r8
      rate(:n,613) = 4e-10_r8
      rate(:n,614) = 2e-10_r8
      rate(:n,615) = 1e-10_r8
      rate(:n,617) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,175) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,176) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,177) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,181) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,183) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,186) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,187) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,196) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,197) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,198) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,201) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,202) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,203) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,211) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,215) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,223) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,224) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,195) = wrk(:)
























      end subroutine setrxt_hrates

      end module mo_setrxt
