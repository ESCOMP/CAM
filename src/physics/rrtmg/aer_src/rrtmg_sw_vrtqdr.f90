!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_vrtqdr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:15 $
!
      module rrtmg_sw_vrtqdr

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only: jpim, jprb
      use parrrsw, only: ngptsw

      implicit none

      contains

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw(ncol, klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input
      integer, intent (in) :: ncol
      integer, intent (in) :: klev                   ! number of model layers
      integer, intent (in) :: kw(ncol)                     ! g-point index

      real(kind=r8), intent(in) :: pref(ncol,klev+1)          ! direct beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: prefd(ncol,klev+1)         ! diffuse beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptra(ncol,klev+1)            ! direct beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptrad(ncol,klev+1)         ! diffuse beam transmissivity
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(in) :: pdbt(ncol,klev+1)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptdbt(ncol,klev+1)
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(inout) :: prdnd(ncol,klev+1)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prup(ncol,klev+1)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prupd(ncol,klev+1)
                                                                 !   Dimensions: (nlayers+1)

! Output
      real(kind=r8), intent(out) :: pfd(ncol,klev+1,ngptsw)   ! downwelling flux (W/m2)
                                                                 !   Dimensions: (nlayers+1,ngptsw)
                                                                 ! unadjusted for earth/sun distance or zenith angle
      real(kind=r8), intent(out) :: pfu(ncol,klev+1,ngptsw)   ! upwelling flux (W/m2)
                                                                 !   Dimensions: (nlayers+1,ngptsw)
                                                                 ! unadjusted for earth/sun distance or zenith angle

! Local

      integer :: ikp, ikx, jk
      integer :: icol
      real(kind=r8) :: zreflect
      real(kind=r8) :: ztdn(ncol,klev+1)  

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
             
      do icol=1,ncol
         zreflect = 1._r8 / (1._r8 - prefd(icol,klev+1) * prefd(icol,klev))
         prup(icol,klev) = pref(icol,klev) + (ptrad(icol,klev) * &
                    ((ptra(icol,klev) - pdbt(icol,klev)) * prefd(icol,klev+1) + &
                      pdbt(icol,klev) * pref(icol,klev+1))) * zreflect
         prupd(icol,klev) = prefd(icol,klev) + ptrad(icol,klev) * ptrad(icol,klev) * &
                       prefd(icol,klev+1) * zreflect

! Pass from bottom to top 
      end do
      do jk = 1,klev-1
         do icol=1,ncol
            ikp = klev+1-jk                       
            ikx = ikp-1
            zreflect = 1._r8 / (1._r8 -prupd(icol,ikp) * prefd(icol,ikx))
            prup(icol,ikx) = pref(icol,ikx) + (ptrad(icol,ikx) * &
                      ((ptra(icol,ikx) - pdbt(icol,ikx)) * prupd(icol,ikp) + &
                        pdbt(icol,ikx) * prup(icol,ikp))) * zreflect
            prupd(icol,ikx) = prefd(icol,ikx) + ptrad(icol,ikx) * ptrad(icol,ikx) * &
                         prupd(icol,ikp) * zreflect
         enddo
      enddo
    
! Upper boundary conditions
      do icol=1,ncol
         ztdn(icol,1) = 1._r8
         prdnd(icol,1) = 0._r8
         ztdn(icol,2) = ptra(icol,1)
         prdnd(icol,2) = prefd(icol,1)

! Pass from top to bottom
      end do
      do jk = 2,klev
         do icol=1,ncol
            ikp = jk+1
            zreflect = 1._r8 / (1._r8 - prefd(icol,jk) * prdnd(icol,jk))
            ztdn(icol,ikp) = ptdbt(icol,jk) * ptra(icol,jk) + &
                       (ptrad(icol,jk) * ((ztdn(icol,jk) - ptdbt(icol,jk)) + &
                        ptdbt(icol,jk) * pref(icol,jk) * prdnd(icol,jk))) * zreflect
            prdnd(icol,ikp) = prefd(icol,jk) + ptrad(icol,jk) * ptrad(icol,jk) * &
                         prdnd(icol,jk) * zreflect
         enddo
      end do
! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         do icol=1,ncol
            zreflect = 1._r8 / (1._r8 - prdnd(icol,jk) * prupd(icol,jk))
            pfu(icol,jk,kw(icol)) = (ptdbt(icol,jk) * prup(icol,jk) + &
                         (ztdn(icol,jk) - ptdbt(icol,jk)) * prupd(icol,jk)) * zreflect
            pfd(icol,jk,kw(icol)) = ptdbt(icol,jk) + (ztdn(icol,jk) - ptdbt(icol,jk)+ &
                         ptdbt(icol,jk) * prup(icol,jk) * prdnd(icol,jk)) * zreflect
         enddo
      end do
      end subroutine vrtqdr_sw

      end module rrtmg_sw_vrtqdr
