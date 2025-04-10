
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
      rate(:,159) = 1.2e-10_r8
      rate(:,160) = 1e-20_r8
      rate(:,161) = 1.3e-16_r8
      rate(:,163) = 4.2e-13_r8
      rate(:,165) = 8e-14_r8
      rate(:,166) = 3.9e-17_r8
      rate(:,173) = 6.9e-12_r8
      rate(:,174) = 7.2e-11_r8
      rate(:,175) = 1.6e-12_r8
      rate(:,181) = 1.8e-12_r8
      rate(:,185) = 1.8e-12_r8
      rate(:,188) = 1.06e-05_r8
      rate(:,190) = 7e-11_r8
      rate(:,191) = 7e-13_r8
      rate(:,199) = 3.5e-12_r8
      rate(:,201) = 1.3e-11_r8
      rate(:,202) = 2.2e-11_r8
      rate(:,203) = 5e-11_r8
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
      rate(:,310) = 4e-13_r8
      rate(:,314) = 1.31e-10_r8
      rate(:,315) = 3.5e-11_r8
      rate(:,316) = 9e-12_r8
      rate(:,323) = 6.8e-14_r8
      rate(:,324) = 2e-13_r8
      rate(:,339) = 1e-12_r8
      rate(:,343) = 1e-14_r8
      rate(:,344) = 1e-11_r8
      rate(:,345) = 1.15e-11_r8
      rate(:,346) = 4e-14_r8
      rate(:,359) = 1.45e-10_r8
      rate(:,360) = 3e-12_r8
      rate(:,361) = 6.7e-13_r8
      rate(:,371) = 3.5e-13_r8
      rate(:,372) = 5.4e-11_r8
      rate(:,375) = 2e-12_r8
      rate(:,376) = 1.4e-11_r8
      rate(:,379) = 2.4e-12_r8
      rate(:,390) = 5e-12_r8
      rate(:,400) = 2.2e-12_r8
      rate(:,402) = 6.7e-12_r8
      rate(:,405) = 3.5e-12_r8
      rate(:,408) = 1.3e-11_r8
      rate(:,409) = 1.4e-11_r8
      rate(:,413) = 2.4e-12_r8
      rate(:,414) = 1.4e-11_r8
      rate(:,419) = 2.4e-12_r8
      rate(:,420) = 4e-11_r8
      rate(:,421) = 4e-11_r8
      rate(:,423) = 1.4e-11_r8
      rate(:,427) = 2.4e-12_r8
      rate(:,428) = 4e-11_r8
      rate(:,432) = 7e-11_r8
      rate(:,433) = 1e-10_r8
      rate(:,438) = 2.4e-12_r8
      rate(:,453) = 4.7e-11_r8
      rate(:,466) = 2.1e-12_r8
      rate(:,467) = 2.8e-13_r8
      rate(:,475) = 1.7e-11_r8
      rate(:,481) = 8.4e-11_r8
      rate(:,483) = 1.9e-11_r8
      rate(:,484) = 1.2e-14_r8
      rate(:,485) = 2e-10_r8
      rate(:,492) = 2.4e-12_r8
      rate(:,493) = 2e-11_r8
      rate(:,497) = 2.3e-11_r8
      rate(:,498) = 2e-11_r8
      rate(:,502) = 3.3e-11_r8
      rate(:,503) = 1e-12_r8
      rate(:,504) = 5.7e-11_r8
      rate(:,505) = 3.4e-11_r8
      rate(:,510) = 2.3e-12_r8
      rate(:,512) = 1.2e-11_r8
      rate(:,513) = 5.7e-11_r8
      rate(:,514) = 2.8e-11_r8
      rate(:,515) = 6.6e-11_r8
      rate(:,516) = 1.4e-11_r8
      rate(:,519) = 1.9e-12_r8
      rate(:,532) = 6.34e-08_r8
      rate(:,538) = 1.9e-11_r8
      rate(:,541) = 1.2e-14_r8
      rate(:,542) = 2e-10_r8
      rate(:,553) = 1.34e-11_r8
      rate(:,559) = 1.34e-11_r8
      rate(:,563) = 1.7e-11_r8
      rate(:,586) = 6e-11_r8
      rate(:,589) = 1e-12_r8
      rate(:,590) = 4e-10_r8
      rate(:,591) = 2e-10_r8
      rate(:,592) = 1e-10_r8
      rate(:,593) = 5e-16_r8
      rate(:,594) = 4.4e-10_r8
      rate(:,595) = 9e-10_r8
      rate(:,598) = 1.29e-07_r8
      rate(:,599) = 2.31e-07_r8
      rate(:,600) = 2.31e-06_r8
      rate(:,601) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,154) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,155) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,156) = 2.64e-11_r8 * exp_fac(:)
      rate(:,157) = 6.6e-12_r8 * exp_fac(:)
      rate(:,162) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,164) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,167) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,168) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,171) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,172) = 1.4e-12_r8 * exp_fac(:)
      rate(:,429) = 1.05e-14_r8 * exp_fac(:)
      rate(:,549) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,177) = 3e-11_r8 * exp_fac(:)
      rate(:,268) = 5.5e-12_r8 * exp_fac(:)
      rate(:,307) = 3.8e-12_r8 * exp_fac(:)
      rate(:,328) = 3.8e-12_r8 * exp_fac(:)
      rate(:,355) = 3.8e-12_r8 * exp_fac(:)
      rate(:,364) = 3.8e-12_r8 * exp_fac(:)
      rate(:,368) = 3.8e-12_r8 * exp_fac(:)
      rate(:,384) = 2.3e-11_r8 * exp_fac(:)
      rate(:,394) = 3.8e-12_r8 * exp_fac(:)
      rate(:,404) = 3.8e-12_r8 * exp_fac(:)
      rate(:,431) = 1.52e-11_r8 * exp_fac(:)
      rate(:,439) = 1.52e-12_r8 * exp_fac(:)
      rate(:,445) = 3.8e-12_r8 * exp_fac(:)
      rate(:,448) = 3.8e-12_r8 * exp_fac(:)
      rate(:,452) = 3.8e-12_r8 * exp_fac(:)
      rate(:,468) = 3.8e-12_r8 * exp_fac(:)
      rate(:,472) = 3.8e-12_r8 * exp_fac(:)
      rate(:,478) = 3.8e-12_r8 * exp_fac(:)
      rate(:,482) = 3.8e-12_r8 * exp_fac(:)
      rate(:,178) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,179) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,180) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,182) = 4.8e-11_r8 * exp_fac(:)
      rate(:,266) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,183) = 1.8e-11_r8 * exp_fac(:)
      rate(:,341) = 4.2e-12_r8 * exp_fac(:)
      rate(:,363) = 4.2e-12_r8 * exp_fac(:)
      rate(:,392) = 4.2e-12_r8 * exp_fac(:)
      rate(:,412) = 4.4e-12_r8 * exp_fac(:)
      rate(:,418) = 4.4e-12_r8 * exp_fac(:)
      rate(:,491) = 4.2e-12_r8 * exp_fac(:)
      rate(:,496) = 4.2e-12_r8 * exp_fac(:)
      rate(:,501) = 4.2e-12_r8 * exp_fac(:)
      rate(:,184) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,189) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,192) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,193) = 2.9e-12_r8 * exp_fac(:)
      rate(:,194) = 1.45e-12_r8 * exp_fac(:)
      rate(:,195) = 1.45e-12_r8 * exp_fac(:)
      rate(:,196) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,197) = 1.2e-13_r8 * exp_fac(:)
      rate(:,226) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,200) = 1.7e-11_r8 * exp_fac(:)
      rate(:,301) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,204) = 3.44e-12_r8 * exp_fac(:)
      rate(:,259) = 2.3e-12_r8 * exp_fac(:)
      rate(:,262) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,205) = 3e-12_r8 * exp_fac(:)
      rate(:,267) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,207) = 7.26e-11_r8 * exp_fac(:)
      rate(:,208) = 4.64e-11_r8 * exp_fac(:)
      rate(:,218) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,219) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,220) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,221) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,222) = 1.4e-11_r8 * exp_fac(:)
      rate(:,236) = 7.4e-12_r8 * exp_fac(:)
      rate(:,337) = 8.1e-12_r8 * exp_fac(:)
      rate(:,223) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,224) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,225) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,227) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,228) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,229) = 2.6e-12_r8 * exp_fac(:)
      rate(:,230) = 6.4e-12_r8 * exp_fac(:)
      rate(:,260) = 4.1e-13_r8 * exp_fac(:)
      rate(:,441) = 7.5e-12_r8 * exp_fac(:)
      rate(:,455) = 7.5e-12_r8 * exp_fac(:)
      rate(:,458) = 7.5e-12_r8 * exp_fac(:)
      rate(:,461) = 7.5e-12_r8 * exp_fac(:)
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
      rate(:,261) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
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
      rate(:,310) = 4e-13_r8 * exp_fac(:)
      rate(:,314) = 1.31e-10_r8 * exp_fac(:)
      rate(:,315) = 3.5e-11_r8 * exp_fac(:)
      rate(:,316) = 9e-12_r8 * exp_fac(:)
      rate(:,323) = 6.8e-14_r8 * exp_fac(:)
      rate(:,324) = 2e-13_r8 * exp_fac(:)
      rate(:,339) = 1e-12_r8 * exp_fac(:)
      rate(:,343) = 1e-14_r8 * exp_fac(:)
      rate(:,344) = 1e-11_r8 * exp_fac(:)
      rate(:,345) = 1.15e-11_r8 * exp_fac(:)
      rate(:,346) = 4e-14_r8 * exp_fac(:)
      rate(:,359) = 1.45e-10_r8 * exp_fac(:)
      rate(:,360) = 3e-12_r8 * exp_fac(:)
      rate(:,361) = 6.7e-13_r8 * exp_fac(:)
      rate(:,371) = 3.5e-13_r8 * exp_fac(:)
      rate(:,372) = 5.4e-11_r8 * exp_fac(:)
      rate(:,375) = 2e-12_r8 * exp_fac(:)
      rate(:,376) = 1.4e-11_r8 * exp_fac(:)
      rate(:,379) = 2.4e-12_r8 * exp_fac(:)
      rate(:,390) = 5e-12_r8 * exp_fac(:)
      rate(:,400) = 2.2e-12_r8 * exp_fac(:)
      rate(:,402) = 6.7e-12_r8 * exp_fac(:)
      rate(:,405) = 3.5e-12_r8 * exp_fac(:)
      rate(:,408) = 1.3e-11_r8 * exp_fac(:)
      rate(:,409) = 1.4e-11_r8 * exp_fac(:)
      rate(:,413) = 2.4e-12_r8 * exp_fac(:)
      rate(:,414) = 1.4e-11_r8 * exp_fac(:)
      rate(:,419) = 2.4e-12_r8 * exp_fac(:)
      rate(:,420) = 4e-11_r8 * exp_fac(:)
      rate(:,421) = 4e-11_r8 * exp_fac(:)
      rate(:,423) = 1.4e-11_r8 * exp_fac(:)
      rate(:,427) = 2.4e-12_r8 * exp_fac(:)
      rate(:,428) = 4e-11_r8 * exp_fac(:)
      rate(:,432) = 7e-11_r8 * exp_fac(:)
      rate(:,433) = 1e-10_r8 * exp_fac(:)
      rate(:,438) = 2.4e-12_r8 * exp_fac(:)
      rate(:,453) = 4.7e-11_r8 * exp_fac(:)
      rate(:,466) = 2.1e-12_r8 * exp_fac(:)
      rate(:,467) = 2.8e-13_r8 * exp_fac(:)
      rate(:,475) = 1.7e-11_r8 * exp_fac(:)
      rate(:,481) = 8.4e-11_r8 * exp_fac(:)
      rate(:,483) = 1.9e-11_r8 * exp_fac(:)
      rate(:,484) = 1.2e-14_r8 * exp_fac(:)
      rate(:,485) = 2e-10_r8 * exp_fac(:)
      rate(:,492) = 2.4e-12_r8 * exp_fac(:)
      rate(:,493) = 2e-11_r8 * exp_fac(:)
      rate(:,497) = 2.3e-11_r8 * exp_fac(:)
      rate(:,498) = 2e-11_r8 * exp_fac(:)
      rate(:,502) = 3.3e-11_r8 * exp_fac(:)
      rate(:,503) = 1e-12_r8 * exp_fac(:)
      rate(:,504) = 5.7e-11_r8 * exp_fac(:)
      rate(:,505) = 3.4e-11_r8 * exp_fac(:)
      rate(:,510) = 2.3e-12_r8 * exp_fac(:)
      rate(:,512) = 1.2e-11_r8 * exp_fac(:)
      rate(:,513) = 5.7e-11_r8 * exp_fac(:)
      rate(:,514) = 2.8e-11_r8 * exp_fac(:)
      rate(:,515) = 6.6e-11_r8 * exp_fac(:)
      rate(:,516) = 1.4e-11_r8 * exp_fac(:)
      rate(:,519) = 1.9e-12_r8 * exp_fac(:)
      rate(:,532) = 6.34e-08_r8 * exp_fac(:)
      rate(:,538) = 1.9e-11_r8 * exp_fac(:)
      rate(:,541) = 1.2e-14_r8 * exp_fac(:)
      rate(:,542) = 2e-10_r8 * exp_fac(:)
      rate(:,553) = 1.34e-11_r8 * exp_fac(:)
      rate(:,559) = 1.34e-11_r8 * exp_fac(:)
      rate(:,563) = 1.7e-11_r8 * exp_fac(:)
      rate(:,586) = 6e-11_r8 * exp_fac(:)
      rate(:,589) = 1e-12_r8 * exp_fac(:)
      rate(:,590) = 4e-10_r8 * exp_fac(:)
      rate(:,591) = 2e-10_r8 * exp_fac(:)
      rate(:,592) = 1e-10_r8 * exp_fac(:)
      rate(:,593) = 5e-16_r8 * exp_fac(:)
      rate(:,594) = 4.4e-10_r8 * exp_fac(:)
      rate(:,595) = 9e-10_r8 * exp_fac(:)
      rate(:,598) = 1.29e-07_r8 * exp_fac(:)
      rate(:,599) = 2.31e-07_r8 * exp_fac(:)
      rate(:,600) = 2.31e-06_r8 * exp_fac(:)
      rate(:,601) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,278) = 6e-12_r8 * exp_fac(:)
      rate(:,377) = 5e-13_r8 * exp_fac(:)
      rate(:,410) = 5e-13_r8 * exp_fac(:)
      rate(:,415) = 5e-13_r8 * exp_fac(:)
      rate(:,424) = 5e-13_r8 * exp_fac(:)
      rate(:,435) = 5e-13_r8 * exp_fac(:)
      rate(:,283) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,284) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,285) = 1.64e-12_r8 * exp_fac(:)
      rate(:,396) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,286) = 2.03e-11_r8 * exp_fac(:)
      rate(:,518) = 3.4e-12_r8 * exp_fac(:)
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
      rate(:,302) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,303) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,304) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,305) = 2.8e-12_r8 * exp_fac(:)
      rate(:,367) = 2.9e-12_r8 * exp_fac(:)
      rate(:,306) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,308) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,311) = 7.5e-13_r8 * exp_fac(:)
      rate(:,325) = 7.5e-13_r8 * exp_fac(:)
      rate(:,340) = 7.5e-13_r8 * exp_fac(:)
      rate(:,362) = 7.5e-13_r8 * exp_fac(:)
      rate(:,366) = 8.6e-13_r8 * exp_fac(:)
      rate(:,378) = 8e-13_r8 * exp_fac(:)
      rate(:,391) = 7.5e-13_r8 * exp_fac(:)
      rate(:,401) = 7.5e-13_r8 * exp_fac(:)
      rate(:,411) = 8e-13_r8 * exp_fac(:)
      rate(:,416) = 8e-13_r8 * exp_fac(:)
      rate(:,425) = 8e-13_r8 * exp_fac(:)
      rate(:,436) = 8e-13_r8 * exp_fac(:)
      rate(:,443) = 7.5e-13_r8 * exp_fac(:)
      rate(:,447) = 7.5e-13_r8 * exp_fac(:)
      rate(:,450) = 7.5e-13_r8 * exp_fac(:)
      rate(:,463) = 7.5e-13_r8 * exp_fac(:)
      rate(:,470) = 7.5e-13_r8 * exp_fac(:)
      rate(:,476) = 7.5e-13_r8 * exp_fac(:)
      rate(:,479) = 7.5e-13_r8 * exp_fac(:)
      rate(:,490) = 7.5e-13_r8 * exp_fac(:)
      rate(:,495) = 7.5e-13_r8 * exp_fac(:)
      rate(:,500) = 7.5e-13_r8 * exp_fac(:)
      rate(:,544) = 7.5e-13_r8 * exp_fac(:)
      rate(:,551) = 7.5e-13_r8 * exp_fac(:)
      rate(:,561) = 7.5e-13_r8 * exp_fac(:)
      rate(:,564) = 7.5e-13_r8 * exp_fac(:)
      rate(:,312) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,313) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,317) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,322) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,326) = 2.6e-12_r8 * exp_fac(:)
      rate(:,444) = 2.6e-12_r8 * exp_fac(:)
      rate(:,449) = 2.6e-12_r8 * exp_fac(:)
      rate(:,451) = 2.6e-12_r8 * exp_fac(:)
      rate(:,464) = 2.6e-12_r8 * exp_fac(:)
      rate(:,471) = 2.6e-12_r8 * exp_fac(:)
      rate(:,477) = 2.6e-12_r8 * exp_fac(:)
      rate(:,480) = 2.6e-12_r8 * exp_fac(:)
      rate(:,545) = 2.6e-12_r8 * exp_fac(:)
      rate(:,552) = 2.6e-12_r8 * exp_fac(:)
      rate(:,562) = 2.6e-12_r8 * exp_fac(:)
      rate(:,565) = 2.6e-12_r8 * exp_fac(:)
      rate(:,327) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,329) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,330) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,331) = 1.4e-12_r8 * exp_fac(:)
      rate(:,351) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,332) = 4.63e-12_r8 * exp_fac(:)
      rate(:,548) = 2.7e-12_r8 * exp_fac(:)
      rate(:,333) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,334) = 2.9e-12_r8 * exp_fac(:)
      rate(:,335) = 2e-12_r8 * exp_fac(:)
      rate(:,365) = 7.1e-13_r8 * exp_fac(:)
      rate(:,386) = 2e-12_r8 * exp_fac(:)
      rate(:,489) = 2e-12_r8 * exp_fac(:)
      rate(:,494) = 2e-12_r8 * exp_fac(:)
      rate(:,499) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,336) = 4.3e-13_r8 * exp_fac(:)
      rate(:,387) = 4.3e-13_r8 * exp_fac(:)
      rate(:,440) = 4.3e-13_r8 * exp_fac(:)
      rate(:,454) = 4.3e-13_r8 * exp_fac(:)
      rate(:,457) = 4.3e-13_r8 * exp_fac(:)
      rate(:,460) = 4.3e-13_r8 * exp_fac(:)
      rate(:,338) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,342) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,350) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,352) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,353) = 1.41e-13_r8 * exp_fac(:)
      rate(:,539) = 2.75e-13_r8 * exp_fac(:)
      rate(:,547) = 2.12e-13_r8 * exp_fac(:)
      rate(:,555) = 2.6e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,354) = 2.7e-12_r8 * exp_fac(:)
      rate(:,380) = 2.7e-12_r8 * exp_fac(:)
      rate(:,381) = 1.3e-13_r8 * exp_fac(:)
      rate(:,383) = 9.6e-12_r8 * exp_fac(:)
      rate(:,389) = 5.3e-12_r8 * exp_fac(:)
      rate(:,426) = 2.7e-12_r8 * exp_fac(:)
      rate(:,437) = 2.7e-12_r8 * exp_fac(:)
      rate(:,540) = 2.7e-12_r8 * exp_fac(:)
      rate(:,556) = 2.7e-12_r8 * exp_fac(:)
      rate(:,356) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,357) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,358) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,373) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,374) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,382) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,385) = 4.6e-12_r8 * exp_fac(:)
      rate(:,388) = 2.3e-12_r8 * exp_fac(:)
      rate(:,393) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,397) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,403) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,406) = 1.86e-11_r8 * exp_fac(:)
      rate(:,407) = 1.86e-11_r8 * exp_fac(:)
      rate(:,417) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,422) = 3.03e-12_r8 * exp_fac(:)
      rate(:,546) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,430) = 2.54e-11_r8 * exp_fac(:)
      rate(:,550) = 2.54e-11_r8 * exp_fac(:)
      rate(:,434) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,442) = 2.3e-12_r8 * exp_fac(:)
      rate(:,543) = 2.3e-12_r8 * exp_fac(:)
      rate(:,446) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,465) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,473) = 1.7e-12_r8 * exp_fac(:)
      rate(:,560) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,486) = 1.2e-12_r8 * exp_fac(:)
      rate(:,554) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,487) = 6.3e-16_r8 * exp_fac(:)
      rate(:,557) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,488) = 1.2e-11_r8 * exp_fac(:)
      rate(:,558) = 1.2e-11_r8 * exp_fac(:)
      rate(:,506) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,507) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,508) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,509) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,517) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,520) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,523) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,176), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,186), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,198), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,206), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,209), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,210), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,211), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,232), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,252), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,263), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,309), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,319), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,320), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,321), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,347), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,348), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,369), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,395), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,398), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,456), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,459), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,462), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,469), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,511), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,160) = 1e-20_r8
      rate(:n,161) = 1.3e-16_r8
      rate(:n,165) = 8e-14_r8
      rate(:n,166) = 3.9e-17_r8
      rate(:n,173) = 6.9e-12_r8
      rate(:n,190) = 7e-11_r8
      rate(:n,191) = 7e-13_r8
      rate(:n,586) = 6e-11_r8
      rate(:n,589) = 1e-12_r8
      rate(:n,590) = 4e-10_r8
      rate(:n,591) = 2e-10_r8
      rate(:n,592) = 1e-10_r8
      rate(:n,594) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,155) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,156) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,157) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,162) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,164) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,167) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,168) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,177) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,178) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,179) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,182) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,183) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,184) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,192) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,196) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,204) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,205) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,176) = wrk(:)
























      end subroutine setrxt_hrates

      end module mo_setrxt
