!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_cldprmc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2008/01/03 21:35:36 $

module rrtmg_sw_cldprmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8

use parrrsw,      only: ngptsw

implicit none

!=========================================================================================
contains
!=========================================================================================

subroutine cldprmc_sw(ncol,nlayers, inflag, iceflag, liqflag, zcldfmc, &
                            ciwpmc, clwpmc, reicmc, dgesmc, relqmc, &
                            ztaormc, ztaucmc, zssacmc, zasmcmc, fsfcmc)

   ! Purpose: Compute the cloud optical properties for each cloudy layer
   ! and g-point interval for use by the McICA method.  
   ! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=2,3 are available;
   ! (Hu & Stamnes, Key, and Fu) are implemented. 

   ! ------- Input -------

   integer, intent(in) :: nlayers         ! total number of layers
   integer, intent(in) :: ncol            ! total number of layers
   integer, intent(in) :: inflag          ! see definitions
   integer, intent(in) :: iceflag         ! see definitions
   integer, intent(in) :: liqflag         ! see definitions

   real(kind=r8), intent(in) :: zcldfmc(ncol,nlayers,ngptsw) ! cloud fraction [mcica]
                                                     !    Dimensions: (ncol,nlayers,ngptsw)
   real(kind=r8), intent(in) :: ciwpmc(ncol,ngptsw,nlayers) ! cloud ice water path [mcica]
                                                     !    Dimensions: (ncol,ngptsw,nlayers)
   real(kind=r8), intent(in) :: clwpmc(ncol,ngptsw,nlayers) ! cloud liquid water path [mcica]
                                                     !    Dimensions: (ncol,ngptsw,nlayers)
   real(kind=r8), intent(in) :: relqmc(ncol,nlayers)  ! cloud liquid particle effective radius (microns)
                                                     !    Dimensions: (ncol,nlayers)
   real(kind=r8), intent(in) :: reicmc(ncol,nlayers) ! cloud ice particle effective radius (microns)
                                                     !    Dimensions: (ncol,nlayers)
   real(kind=r8), intent(in) :: dgesmc(ncol,nlayers) ! cloud ice particle generalized effective size (microns)
                                                     !    Dimensions: (ncol,nlayers)
   real(kind=r8), intent(in) :: fsfcmc(ncol,ngptsw,nlayers) ! cloud forward scattering fraction 
                                                     !    Dimensions: (ncol,ngptsw,nlayers)

   ! ------- Output -------

   real(kind=r8), intent(inout) :: ztaucmc(ncol,nlayers,ngptsw) ! cloud optical depth (delta scaled)
                                                     !    Dimensions: (ncol,nlayers,ngptsw)
   real(kind=r8), intent(inout) :: zssacmc(ncol,nlayers,ngptsw) ! single scattering albedo (delta scaled)
                                                     !    Dimensions: (ncol,nlayers,ngptsw)
   real(kind=r8), intent(inout) :: zasmcmc(ncol,nlayers,ngptsw) ! asymmetry parameter (delta scaled)
                                                     !    Dimensions: (ncol,nlayers,ngptsw)
   real(kind=r8), intent(out) :: ztaormc(ncol,nlayers,ngptsw)   ! cloud optical depth (non-delta scaled)
                                                     !    Dimensions: (ncol,nlayers,ngptsw)

   ! ------- Local -------

   integer                  :: lay, ig, iplon
   real(kind=r8), parameter :: eps = 1.e-06_r8     ! epsilon
   real(kind=r8), parameter :: cldmin = 1.e-80_r8  ! minimum value for cloud quantities
   real(kind=r8) :: cwp(ncol)                            ! total cloud water path

   real(kind=r8), dimension(ncol) :: taucldorig_a, taucloud_a, ssacloud_a, ffp, ffp1, ffpssa
   !----------------------------------------------------------------------------

   ! Main layer loop
   do lay = 1, nlayers

      ! Main g-point interval loop
      do ig = 1, ngptsw 

         do iplon=1, ncol

            ztaormc(iplon,lay,ig) = ztaucmc(iplon,lay,ig)
            cwp(iplon) = ciwpmc(iplon,ig,lay) + clwpmc(iplon,ig,lay)
            if (zcldfmc(iplon,lay,ig) .ge. cldmin .and. &
               (cwp(iplon) .ge. cldmin .or. ztaucmc(iplon,lay,ig) .ge. cldmin)) then

               ! (inflag=0): Cloud optical properties input directly
               if (inflag .eq. 0) then
                  ! Cloud optical properties already defined in ztaucmc, zssacmc, zasmcmc are unscaled;
                  ! Apply delta-M scaling here (using Henyey-Greenstein approximation)
                  taucldorig_a(iplon) = ztaucmc(iplon,lay,ig)
                  ffp(iplon) = fsfcmc(iplon,ig,lay)
                  ffp1(iplon) = 1.0_r8 - ffp(iplon)
                  ffpssa(iplon) = 1.0_r8 - ffp(iplon) * zssacmc(iplon,lay,ig)
                  ssacloud_a(iplon) = ffp1(iplon) * zssacmc(iplon,lay,ig) / ffpssa(iplon)
                  taucloud_a(iplon) = ffpssa(iplon) * taucldorig_a(iplon)

                  ztaormc(iplon,lay,ig) = taucldorig_a(iplon)
                  zssacmc(iplon,lay,ig) = ssacloud_a(iplon)
                  ztaucmc(iplon,lay,ig) = taucloud_a(iplon)
                  zasmcmc(iplon,lay,ig) = (zasmcmc(iplon,lay,ig) - ffp(iplon)) / (ffp1(iplon))

               end if
            end if

         end do
      end do
   end do

end subroutine cldprmc_sw

end module rrtmg_sw_cldprmc
