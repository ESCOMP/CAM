
subroutine hdinti(rearth, deltat)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Time independent initialization for the horizontal diffusion.
! 
! Method: 
! 
! Author: 
! Original version:  D. Williamson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Reviewed:          B. Boville, April 1996
!
!-----------------------------------------------------------------------

   use shr_kind_mod,   only: r8=>shr_kind_r8
   use cam_abortutils, only: endrun
   use pmgrid
   use pspect
   use eul_control_mod
   use cam_logfile,    only: iulog
   implicit none

!------------------------------Arguments--------------------------------

   real(r8), intent(in) :: rearth               ! radius of the earth
   real(r8), intent(in) :: deltat               ! time step

!---------------------------Local workspace-----------------------------

   integer :: k             ! level index
   integer :: n             ! n-wavenumber index
   integer :: iexpon
   real(r8) :: fn
!
!-----------------------------------------------------------------------
!
! Initialize physical constants for courant number based spect truncation
!
   nmaxhd = ptrk
   cnlim  = 0.999_r8          ! maximum allowable Courant number
   cnfac  = deltat*real(nmaxhd,r8)/rearth
!
! Initialize arrays used for courant number based spectral truncation
!
   do k=1,plev
      nindex(k) = 2*nmaxhd
   end do
!
! Set the Del^2 and Del^N diffusion coefficients for each wavenumber
!
   hdfst2(1) = 0._r8
   hdfsd2(1) = 0._r8
!
   hdfstn(1) = 0._r8
   hdfsdn(1) = 0._r8

   iexpon    = hdif_order/2

   do n=2,pnmax

      hdfst2(n) = dif2 * (n*(n-1)  ) / rearth**2
      hdfsd2(n) = dif2 * (n*(n-1)-2) / rearth**2

      fn        = n*(n-1)
      fn        = fn/rearth**2
      fn        = fn**iexpon
      
      hdfstn(n) = hdif_coef * fn
      fn        = 2._r8/rearth**2 
      hdfsdn(n) = hdfstn(n) - hdif_coef * fn**iexpon

   end do
!
   return
end subroutine hdinti

