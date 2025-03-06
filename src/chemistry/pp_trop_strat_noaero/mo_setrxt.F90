
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

      rate(:,114) = 1.2e-10_r8
      rate(:,118) = 1.2e-10_r8
      rate(:,124) = 6.9e-12_r8
      rate(:,125) = 7.2e-11_r8
      rate(:,126) = 1.6e-12_r8
      rate(:,132) = 1.8e-12_r8
      rate(:,136) = 1.8e-12_r8
      rate(:,148) = 3.5e-12_r8
      rate(:,150) = 1e-11_r8
      rate(:,151) = 2.2e-11_r8
      rate(:,152) = 5e-11_r8
      rate(:,187) = 1.7e-13_r8
      rate(:,189) = 2.607e-10_r8
      rate(:,190) = 9.75e-11_r8
      rate(:,191) = 2.07e-10_r8
      rate(:,192) = 2.088e-10_r8
      rate(:,193) = 1.17e-10_r8
      rate(:,194) = 4.644e-11_r8
      rate(:,195) = 1.204e-10_r8
      rate(:,196) = 9.9e-11_r8
      rate(:,197) = 3.3e-12_r8
      rate(:,216) = 4.5e-11_r8
      rate(:,217) = 4.62e-10_r8
      rate(:,218) = 1.2e-10_r8
      rate(:,219) = 9e-11_r8
      rate(:,220) = 3e-11_r8
      rate(:,225) = 2.14e-11_r8
      rate(:,226) = 1.9e-10_r8
      rate(:,239) = 2.57e-10_r8
      rate(:,240) = 1.8e-10_r8
      rate(:,241) = 1.794e-10_r8
      rate(:,242) = 1.3e-10_r8
      rate(:,243) = 7.65e-11_r8
      rate(:,257) = 4e-13_r8
      rate(:,261) = 1.31e-10_r8
      rate(:,262) = 3.5e-11_r8
      rate(:,263) = 9e-12_r8
      rate(:,270) = 6.8e-14_r8
      rate(:,271) = 2e-13_r8
      rate(:,285) = 7e-13_r8
      rate(:,286) = 1e-12_r8
      rate(:,290) = 1e-14_r8
      rate(:,291) = 1e-11_r8
      rate(:,292) = 1.15e-11_r8
      rate(:,293) = 4e-14_r8
      rate(:,306) = 3e-12_r8
      rate(:,307) = 6.7e-13_r8
      rate(:,317) = 3.5e-13_r8
      rate(:,318) = 5.4e-11_r8
      rate(:,321) = 2e-12_r8
      rate(:,322) = 1.4e-11_r8
      rate(:,325) = 2.4e-12_r8
      rate(:,336) = 5e-12_r8
      rate(:,346) = 1.6e-12_r8
      rate(:,348) = 6.7e-12_r8
      rate(:,351) = 3.5e-12_r8
      rate(:,354) = 1.3e-11_r8
      rate(:,355) = 1.4e-11_r8
      rate(:,359) = 2.4e-12_r8
      rate(:,360) = 1.4e-11_r8
      rate(:,365) = 2.4e-12_r8
      rate(:,366) = 4e-11_r8
      rate(:,367) = 4e-11_r8
      rate(:,369) = 1.4e-11_r8
      rate(:,373) = 2.4e-12_r8
      rate(:,374) = 4e-11_r8
      rate(:,378) = 7e-11_r8
      rate(:,379) = 1e-10_r8
      rate(:,384) = 2.4e-12_r8
      rate(:,399) = 4.7e-11_r8
      rate(:,412) = 2.1e-12_r8
      rate(:,413) = 2.8e-13_r8
      rate(:,421) = 1.7e-11_r8
      rate(:,427) = 8.4e-11_r8
      rate(:,429) = 1.9e-11_r8
      rate(:,430) = 1.2e-14_r8
      rate(:,431) = 2e-10_r8
      rate(:,438) = 2.4e-12_r8
      rate(:,439) = 2e-11_r8
      rate(:,443) = 2.3e-11_r8
      rate(:,444) = 2e-11_r8
      rate(:,448) = 3.3e-11_r8
      rate(:,449) = 1e-12_r8
      rate(:,450) = 5.7e-11_r8
      rate(:,451) = 3.4e-11_r8
      rate(:,456) = 2.3e-12_r8
      rate(:,457) = 1.2e-11_r8
      rate(:,458) = 5.7e-11_r8
      rate(:,459) = 2.8e-11_r8
      rate(:,460) = 6.6e-11_r8
      rate(:,461) = 1.4e-11_r8
      rate(:,464) = 1.9e-12_r8
      rate(:,478) = 6.34e-08_r8
      rate(:,484) = 1.9e-11_r8
      rate(:,487) = 1.2e-14_r8
      rate(:,488) = 2e-10_r8
      rate(:,499) = 1.34e-11_r8
      rate(:,505) = 1.34e-11_r8
      rate(:,509) = 1.7e-11_r8
      rate(:,529) = 1.29e-07_r8
      rate(:,530) = 2.31e-07_r8
      rate(:,531) = 2.31e-06_r8
      rate(:,532) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,115) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,116) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,117) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,119) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,122) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,123) = 1.4e-12_r8 * exp_fac(:)
      rate(:,375) = 1.05e-14_r8 * exp_fac(:)
      rate(:,495) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,128) = 3e-11_r8 * exp_fac(:)
      rate(:,214) = 5.5e-12_r8 * exp_fac(:)
      rate(:,253) = 3.8e-12_r8 * exp_fac(:)
      rate(:,275) = 3.8e-12_r8 * exp_fac(:)
      rate(:,302) = 3.8e-12_r8 * exp_fac(:)
      rate(:,310) = 3.8e-12_r8 * exp_fac(:)
      rate(:,314) = 3.8e-12_r8 * exp_fac(:)
      rate(:,330) = 2.3e-11_r8 * exp_fac(:)
      rate(:,340) = 3.8e-12_r8 * exp_fac(:)
      rate(:,350) = 3.8e-12_r8 * exp_fac(:)
      rate(:,377) = 1.52e-11_r8 * exp_fac(:)
      rate(:,385) = 1.52e-12_r8 * exp_fac(:)
      rate(:,391) = 3.8e-12_r8 * exp_fac(:)
      rate(:,394) = 3.8e-12_r8 * exp_fac(:)
      rate(:,398) = 3.8e-12_r8 * exp_fac(:)
      rate(:,414) = 3.8e-12_r8 * exp_fac(:)
      rate(:,418) = 3.8e-12_r8 * exp_fac(:)
      rate(:,424) = 3.8e-12_r8 * exp_fac(:)
      rate(:,428) = 3.8e-12_r8 * exp_fac(:)
      rate(:,129) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,130) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,131) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,133) = 4.8e-11_r8 * exp_fac(:)
      rate(:,212) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,134) = 1.8e-11_r8 * exp_fac(:)
      rate(:,288) = 4.2e-12_r8 * exp_fac(:)
      rate(:,301) = 4.2e-12_r8 * exp_fac(:)
      rate(:,309) = 4.2e-12_r8 * exp_fac(:)
      rate(:,338) = 4.2e-12_r8 * exp_fac(:)
      rate(:,358) = 4.4e-12_r8 * exp_fac(:)
      rate(:,364) = 4.4e-12_r8 * exp_fac(:)
      rate(:,437) = 4.2e-12_r8 * exp_fac(:)
      rate(:,442) = 4.2e-12_r8 * exp_fac(:)
      rate(:,447) = 4.2e-12_r8 * exp_fac(:)
      rate(:,135) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,139) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,140) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,141) = 2.9e-12_r8 * exp_fac(:)
      rate(:,142) = 1.45e-12_r8 * exp_fac(:)
      rate(:,143) = 1.45e-12_r8 * exp_fac(:)
      rate(:,144) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,145) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,146) = 1.2e-13_r8 * exp_fac(:)
      rate(:,172) = 3e-11_r8 * exp_fac(:)
      rate(:,149) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,153) = 3.3e-12_r8 * exp_fac(:)
      rate(:,168) = 1.4e-11_r8 * exp_fac(:)
      rate(:,182) = 7.4e-12_r8 * exp_fac(:)
      rate(:,284) = 8.1e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,154) = 3e-12_r8 * exp_fac(:)
      rate(:,213) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,156) = 7.26e-11_r8 * exp_fac(:)
      rate(:,157) = 4.64e-11_r8 * exp_fac(:)
      rate(:,164) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,165) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,166) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,167) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,169) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,170) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,171) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,173) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,174) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,175) = 2.6e-12_r8 * exp_fac(:)
      rate(:,176) = 6.4e-12_r8 * exp_fac(:)
      rate(:,206) = 4.1e-13_r8 * exp_fac(:)
      rate(:,387) = 7.5e-12_r8 * exp_fac(:)
      rate(:,401) = 7.5e-12_r8 * exp_fac(:)
      rate(:,404) = 7.5e-12_r8 * exp_fac(:)
      rate(:,407) = 7.5e-12_r8 * exp_fac(:)
      rate(:,177) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,179) = 3.6e-12_r8 * exp_fac(:)
      rate(:,228) = 2e-12_r8 * exp_fac(:)
      rate(:,180) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,181) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,183) = 6e-13_r8 * exp_fac(:)
      rate(:,203) = 1.5e-12_r8 * exp_fac(:)
      rate(:,211) = 1.9e-11_r8 * exp_fac(:)
      rate(:,184) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,185) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,186) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,188) = 3e-12_r8 * exp_fac(:)
      rate(:,222) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,200) = 1.7e-11_r8 * exp_fac(:)
      rate(:,227) = 6.3e-12_r8 * exp_fac(:)
      rate(:,201) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,202) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,204) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,205) = 2.3e-12_r8 * exp_fac(:)
      rate(:,208) = 8.8e-12_r8 * exp_fac(:)
      rate(:,207) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,210) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,215) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,221) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,223) = 1.4e-11_r8 * exp_fac(:)
      rate(:,225) = 2.14e-11_r8 * exp_fac(:)
      rate(:,226) = 1.9e-10_r8 * exp_fac(:)
      rate(:,239) = 2.57e-10_r8 * exp_fac(:)
      rate(:,240) = 1.8e-10_r8 * exp_fac(:)
      rate(:,241) = 1.794e-10_r8 * exp_fac(:)
      rate(:,242) = 1.3e-10_r8 * exp_fac(:)
      rate(:,243) = 7.65e-11_r8 * exp_fac(:)
      rate(:,257) = 4e-13_r8 * exp_fac(:)
      rate(:,261) = 1.31e-10_r8 * exp_fac(:)
      rate(:,262) = 3.5e-11_r8 * exp_fac(:)
      rate(:,263) = 9e-12_r8 * exp_fac(:)
      rate(:,270) = 6.8e-14_r8 * exp_fac(:)
      rate(:,271) = 2e-13_r8 * exp_fac(:)
      rate(:,285) = 7e-13_r8 * exp_fac(:)
      rate(:,286) = 1e-12_r8 * exp_fac(:)
      rate(:,290) = 1e-14_r8 * exp_fac(:)
      rate(:,291) = 1e-11_r8 * exp_fac(:)
      rate(:,292) = 1.15e-11_r8 * exp_fac(:)
      rate(:,293) = 4e-14_r8 * exp_fac(:)
      rate(:,306) = 3e-12_r8 * exp_fac(:)
      rate(:,307) = 6.7e-13_r8 * exp_fac(:)
      rate(:,317) = 3.5e-13_r8 * exp_fac(:)
      rate(:,318) = 5.4e-11_r8 * exp_fac(:)
      rate(:,321) = 2e-12_r8 * exp_fac(:)
      rate(:,322) = 1.4e-11_r8 * exp_fac(:)
      rate(:,325) = 2.4e-12_r8 * exp_fac(:)
      rate(:,336) = 5e-12_r8 * exp_fac(:)
      rate(:,346) = 1.6e-12_r8 * exp_fac(:)
      rate(:,348) = 6.7e-12_r8 * exp_fac(:)
      rate(:,351) = 3.5e-12_r8 * exp_fac(:)
      rate(:,354) = 1.3e-11_r8 * exp_fac(:)
      rate(:,355) = 1.4e-11_r8 * exp_fac(:)
      rate(:,359) = 2.4e-12_r8 * exp_fac(:)
      rate(:,360) = 1.4e-11_r8 * exp_fac(:)
      rate(:,365) = 2.4e-12_r8 * exp_fac(:)
      rate(:,366) = 4e-11_r8 * exp_fac(:)
      rate(:,367) = 4e-11_r8 * exp_fac(:)
      rate(:,369) = 1.4e-11_r8 * exp_fac(:)
      rate(:,373) = 2.4e-12_r8 * exp_fac(:)
      rate(:,374) = 4e-11_r8 * exp_fac(:)
      rate(:,378) = 7e-11_r8 * exp_fac(:)
      rate(:,379) = 1e-10_r8 * exp_fac(:)
      rate(:,384) = 2.4e-12_r8 * exp_fac(:)
      rate(:,399) = 4.7e-11_r8 * exp_fac(:)
      rate(:,412) = 2.1e-12_r8 * exp_fac(:)
      rate(:,413) = 2.8e-13_r8 * exp_fac(:)
      rate(:,421) = 1.7e-11_r8 * exp_fac(:)
      rate(:,427) = 8.4e-11_r8 * exp_fac(:)
      rate(:,429) = 1.9e-11_r8 * exp_fac(:)
      rate(:,430) = 1.2e-14_r8 * exp_fac(:)
      rate(:,431) = 2e-10_r8 * exp_fac(:)
      rate(:,438) = 2.4e-12_r8 * exp_fac(:)
      rate(:,439) = 2e-11_r8 * exp_fac(:)
      rate(:,443) = 2.3e-11_r8 * exp_fac(:)
      rate(:,444) = 2e-11_r8 * exp_fac(:)
      rate(:,448) = 3.3e-11_r8 * exp_fac(:)
      rate(:,449) = 1e-12_r8 * exp_fac(:)
      rate(:,450) = 5.7e-11_r8 * exp_fac(:)
      rate(:,451) = 3.4e-11_r8 * exp_fac(:)
      rate(:,456) = 2.3e-12_r8 * exp_fac(:)
      rate(:,457) = 1.2e-11_r8 * exp_fac(:)
      rate(:,458) = 5.7e-11_r8 * exp_fac(:)
      rate(:,459) = 2.8e-11_r8 * exp_fac(:)
      rate(:,460) = 6.6e-11_r8 * exp_fac(:)
      rate(:,461) = 1.4e-11_r8 * exp_fac(:)
      rate(:,464) = 1.9e-12_r8 * exp_fac(:)
      rate(:,478) = 6.34e-08_r8 * exp_fac(:)
      rate(:,484) = 1.9e-11_r8 * exp_fac(:)
      rate(:,487) = 1.2e-14_r8 * exp_fac(:)
      rate(:,488) = 2e-10_r8 * exp_fac(:)
      rate(:,499) = 1.34e-11_r8 * exp_fac(:)
      rate(:,505) = 1.34e-11_r8 * exp_fac(:)
      rate(:,509) = 1.7e-11_r8 * exp_fac(:)
      rate(:,529) = 1.29e-07_r8 * exp_fac(:)
      rate(:,530) = 2.31e-07_r8 * exp_fac(:)
      rate(:,531) = 2.31e-06_r8 * exp_fac(:)
      rate(:,532) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,224) = 6e-12_r8 * exp_fac(:)
      rate(:,323) = 5e-13_r8 * exp_fac(:)
      rate(:,356) = 5e-13_r8 * exp_fac(:)
      rate(:,361) = 5e-13_r8 * exp_fac(:)
      rate(:,370) = 5e-13_r8 * exp_fac(:)
      rate(:,381) = 5e-13_r8 * exp_fac(:)
      rate(:,229) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,230) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,231) = 1.64e-12_r8 * exp_fac(:)
      rate(:,342) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,232) = 2.03e-11_r8 * exp_fac(:)
      rate(:,463) = 3.4e-12_r8 * exp_fac(:)
      rate(:,233) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,234) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,235) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,236) = 1.25e-12_r8 * exp_fac(:)
      rate(:,246) = 3.4e-11_r8 * exp_fac(:)
      rate(:,237) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,238) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,244) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,245) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,247) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,248) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,249) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,250) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,251) = 2.8e-12_r8 * exp_fac(:)
      rate(:,313) = 2.9e-12_r8 * exp_fac(:)
      rate(:,252) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,254) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,258) = 7.5e-13_r8 * exp_fac(:)
      rate(:,272) = 7.5e-13_r8 * exp_fac(:)
      rate(:,287) = 7.5e-13_r8 * exp_fac(:)
      rate(:,300) = 7.5e-13_r8 * exp_fac(:)
      rate(:,308) = 7.5e-13_r8 * exp_fac(:)
      rate(:,312) = 8.6e-13_r8 * exp_fac(:)
      rate(:,324) = 8e-13_r8 * exp_fac(:)
      rate(:,337) = 7.5e-13_r8 * exp_fac(:)
      rate(:,347) = 7.5e-13_r8 * exp_fac(:)
      rate(:,357) = 8e-13_r8 * exp_fac(:)
      rate(:,362) = 8e-13_r8 * exp_fac(:)
      rate(:,371) = 8e-13_r8 * exp_fac(:)
      rate(:,382) = 8e-13_r8 * exp_fac(:)
      rate(:,389) = 7.5e-13_r8 * exp_fac(:)
      rate(:,393) = 7.5e-13_r8 * exp_fac(:)
      rate(:,396) = 7.5e-13_r8 * exp_fac(:)
      rate(:,409) = 7.5e-13_r8 * exp_fac(:)
      rate(:,416) = 7.5e-13_r8 * exp_fac(:)
      rate(:,422) = 7.5e-13_r8 * exp_fac(:)
      rate(:,425) = 7.5e-13_r8 * exp_fac(:)
      rate(:,436) = 7.5e-13_r8 * exp_fac(:)
      rate(:,441) = 7.5e-13_r8 * exp_fac(:)
      rate(:,446) = 7.5e-13_r8 * exp_fac(:)
      rate(:,490) = 7.5e-13_r8 * exp_fac(:)
      rate(:,497) = 7.5e-13_r8 * exp_fac(:)
      rate(:,507) = 7.5e-13_r8 * exp_fac(:)
      rate(:,510) = 7.5e-13_r8 * exp_fac(:)
      rate(:,259) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,260) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,264) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,269) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,273) = 2.6e-12_r8 * exp_fac(:)
      rate(:,390) = 2.6e-12_r8 * exp_fac(:)
      rate(:,395) = 2.6e-12_r8 * exp_fac(:)
      rate(:,397) = 2.6e-12_r8 * exp_fac(:)
      rate(:,410) = 2.6e-12_r8 * exp_fac(:)
      rate(:,417) = 2.6e-12_r8 * exp_fac(:)
      rate(:,423) = 2.6e-12_r8 * exp_fac(:)
      rate(:,426) = 2.6e-12_r8 * exp_fac(:)
      rate(:,491) = 2.6e-12_r8 * exp_fac(:)
      rate(:,498) = 2.6e-12_r8 * exp_fac(:)
      rate(:,508) = 2.6e-12_r8 * exp_fac(:)
      rate(:,511) = 2.6e-12_r8 * exp_fac(:)
      rate(:,274) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,276) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,277) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,278) = 1.4e-12_r8 * exp_fac(:)
      rate(:,298) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,279) = 4.63e-12_r8 * exp_fac(:)
      rate(:,494) = 2.7e-12_r8 * exp_fac(:)
      rate(:,280) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,281) = 2.9e-12_r8 * exp_fac(:)
      rate(:,282) = 2e-12_r8 * exp_fac(:)
      rate(:,311) = 7.1e-13_r8 * exp_fac(:)
      rate(:,332) = 2e-12_r8 * exp_fac(:)
      rate(:,435) = 2e-12_r8 * exp_fac(:)
      rate(:,440) = 2e-12_r8 * exp_fac(:)
      rate(:,445) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,283) = 4.3e-13_r8 * exp_fac(:)
      rate(:,333) = 4.3e-13_r8 * exp_fac(:)
      rate(:,386) = 4.3e-13_r8 * exp_fac(:)
      rate(:,400) = 4.3e-13_r8 * exp_fac(:)
      rate(:,403) = 4.3e-13_r8 * exp_fac(:)
      rate(:,406) = 4.3e-13_r8 * exp_fac(:)
      rate(:,289) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,297) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,299) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,303) = 8.7e-12_r8 * exp( -615._r8 * itemp(:) )
      rate(:,304) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,305) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,319) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,320) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,326) = 2.7e-12_r8 * exp_fac(:)
      rate(:,327) = 1.3e-13_r8 * exp_fac(:)
      rate(:,329) = 9.6e-12_r8 * exp_fac(:)
      rate(:,335) = 5.3e-12_r8 * exp_fac(:)
      rate(:,372) = 2.7e-12_r8 * exp_fac(:)
      rate(:,383) = 2.7e-12_r8 * exp_fac(:)
      rate(:,486) = 2.7e-12_r8 * exp_fac(:)
      rate(:,502) = 2.7e-12_r8 * exp_fac(:)
      rate(:,328) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,331) = 4.6e-12_r8 * exp_fac(:)
      rate(:,334) = 2.3e-12_r8 * exp_fac(:)
      rate(:,339) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,343) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,349) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,352) = 1.86e-11_r8 * exp_fac(:)
      rate(:,353) = 1.86e-11_r8 * exp_fac(:)
      rate(:,363) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,368) = 3.03e-12_r8 * exp_fac(:)
      rate(:,492) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,376) = 2.54e-11_r8 * exp_fac(:)
      rate(:,496) = 2.54e-11_r8 * exp_fac(:)
      rate(:,380) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,388) = 2.3e-12_r8 * exp_fac(:)
      rate(:,489) = 2.3e-12_r8 * exp_fac(:)
      rate(:,392) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,411) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,419) = 1.7e-12_r8 * exp_fac(:)
      rate(:,506) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,432) = 1.2e-12_r8 * exp_fac(:)
      rate(:,500) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,433) = 6.3e-16_r8 * exp_fac(:)
      rate(:,503) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,434) = 1.2e-11_r8 * exp_fac(:)
      rate(:,504) = 1.2e-11_r8 * exp_fac(:)
      rate(:,452) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,453) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )
      rate(:,454) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,455) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,462) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,465) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,469) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,485) = 2.75e-13_r8 * exp_fac(:)
      rate(:,493) = 2.12e-13_r8 * exp_fac(:)
      rate(:,501) = 2.6e-13_r8 * exp_fac(:)

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,127), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,137), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,147), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,155), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,158), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,159), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,160), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,178), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,198), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,209), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,255), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.28e-33_r8
      kinf(:) = 9.3e-15_r8 * itemp(:)**(-4.42_r8)
      call jpl( rate(:,256), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,266), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,267), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,268), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,294), m, 0.48_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,295), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,315), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,341), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,402), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,405), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,408), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,415), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,124) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,116) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,119) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,128) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,129) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,130) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,133) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,134) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,135) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,140) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,144) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,145) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,153) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,154) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,127) = wrk(:)























      end subroutine setrxt_hrates

      end module mo_setrxt
