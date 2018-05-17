!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_spcvmc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:14 $

      module rrtmg_sw_spcvmc

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

      use shr_kind_mod,    only: r8 => shr_kind_r8

      use parrrsw,         only: nbndsw, ngptsw, mxmol, jpband
      use rrsw_tbl,        only: tblint, bpade, od_lo, exp_tbl
      use rrsw_wvn,        only: ngc, ngs
      use rrtmg_sw_reftra, only: reftra_sw
      use rrtmg_sw_taumol, only: taumol_sw
      use rrtmg_sw_vrtqdr, only: vrtqdr_sw

      implicit none

      contains

! ---------------------------------------------------------------------------
      subroutine spcvmc_sw &
            (lchnk, ncol, nlayers, istart, iend, icpr, idelm, iout, &
             pavel, tavel, pz, tz, tbound, palbd, palbp, &
             pcldfmc, ptaucmc, pasycmc, pomgcmc, ptaormc, &
             ptaua, pasya, pomga, prmu0, coldry, wkl, adjflux, &
             laytrop, layswtch, laylow, jp, jt, jt1, &
             co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
             fac00, fac01, fac10, fac11, &
             selffac, selffrac, indself, forfac, forfrac, indfor, &
             pbbfd, pbbfu, pbbcd, pbbcu, puvfd, puvcd, pnifd, pnicd, pnifu, pnicu, &
             pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir, &
             pbbfsu, pbbfsd)
! ---------------------------------------------------------------------------
!
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker and McICA, the Monte-Carlo
!          Independent Column Approximation, for the representation of 
!          sub-grid cloud variability (i.e. cloud overlap).
!
! Interface:  *spcvmc_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:
!
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Modified to use McICA: MJIacono, AER, Nov 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer, intent(in) :: lchnk

      integer, intent(in) :: nlayers
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      integer, intent(in) :: icpr
      integer, intent(in) :: idelm     ! delta-m scaling flag
                                       ! [0 = direct and diffuse fluxes are unscaled]
                                       ! [1 = direct and diffuse fluxes are scaled]
      integer, intent(in) :: iout
      integer, intent(in) :: ncol 

      integer, intent(in) :: laytrop(ncol)
      integer, intent(in) :: layswtch(ncol)
      integer, intent(in) :: laylow(ncol)

      integer, intent(in) :: indfor(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      integer, intent(in) :: indself(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      integer, intent(in) :: jp(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      integer, intent(in) :: jt(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      integer, intent(in) :: jt1(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)

      real(kind=r8), intent(in) :: pavel(ncol,nlayers)           ! layer pressure (hPa, mb) 
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: tavel(ncol,nlayers)           ! layer temperature (K)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: pz(ncol,0:nlayers)            ! level (interface) pressure (hPa, mb)
                                                                 !   Dimensions: (ncol,0:nlayers)
      real(kind=r8), intent(in) :: tz(ncol,0:nlayers)            ! level temperatures (hPa, mb)
                                                                 !   Dimensions: (ncol,0:nlayers)
      real(kind=r8), intent(in) :: tbound(ncol)                  ! surface temperature (K)
      real(kind=r8), intent(in) :: wkl(ncol,mxmol,nlayers)       ! molecular amounts (mol/cm2) 
                                                                 !   Dimensions: (ncol,mxmol,nlayers)
      real(kind=r8), intent(in) :: coldry(ncol,nlayers)          ! dry air column density (mol/cm2)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: colmol(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: adjflux(ncol,jpband)          ! Earth/Sun distance adjustment
                                                                 !   Dimensions: (ncol,jpband)

      real(kind=r8), intent(in) :: palbd(ncol,nbndsw)            ! surface albedo (diffuse)
                                                                 !   Dimensions: (ncol,nbndsw)
      real(kind=r8), intent(in) :: palbp(ncol,nbndsw)            ! surface albedo (direct)
                                                                 !   Dimensions: (ncol, nbndsw)
      real(kind=r8), intent(in) :: prmu0(ncol)                   ! cosine of solar zenith angle
      real(kind=r8), intent(in) :: pcldfmc(ncol,nlayers,ngptsw)  ! cloud fraction [mcica]
                                                                 !   Dimensions: (ncol,nlayers,ngptsw)
      real(kind=r8), intent(in) :: ptaucmc(ncol,nlayers,ngptsw)  ! cloud optical depth [mcica]
                                                                 !   Dimensions: (ncol,nlayers,ngptsw)
      real(kind=r8), intent(in) :: pasycmc(ncol,nlayers,ngptsw)  ! cloud asymmetry parameter [mcica]
                                                                 !   Dimensions: (ncol,nlayers,ngptsw)
      real(kind=r8), intent(in) :: pomgcmc(ncol,nlayers,ngptsw)  ! cloud single scattering albedo [mcica]
                                                                 !   Dimensions: (ncol,nlayers,ngptsw)
      real(kind=r8), intent(in) :: ptaormc(ncol,nlayers,ngptsw)  ! cloud optical depth, non-delta scaled [mcica]
                                                                 !   Dimensions: (ncol,nlayers,ngptsw)
      real(kind=r8), intent(in) :: ptaua(ncol,nlayers,nbndsw)    ! aerosol optical depth
                                                                 !   Dimensions: (ncol,nlayers,nbndsw)
      real(kind=r8), intent(in) :: pasya(ncol,nlayers,nbndsw)    ! aerosol asymmetry parameter
                                                                 !   Dimensions: (ncol,nlayers,nbndsw)
      real(kind=r8), intent(in) :: pomga(ncol,nlayers,nbndsw)    ! aerosol single scattering albedo
                                                                 !   Dimensions: (ncol,nlayers,nbndsw)

      real(kind=r8), intent(in) :: colh2o(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: colco2(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: colch4(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: co2mult(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: colo3(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: colo2(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: coln2o(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: forfac(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: forfrac(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: selffac(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: selffrac(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: fac00(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: fac01(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: fac10(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)
      real(kind=r8), intent(in) :: fac11(ncol,nlayers)
                                                                 !   Dimensions: (ncol,nlayers)

! ------- Output -------
                                                                 !   All Dimensions: (nlayers+1)
      real(kind=r8), intent(out) :: pbbcd(:,:)
      real(kind=r8), intent(out) :: pbbcu(:,:)
      real(kind=r8), intent(out) :: pbbfd(:,:)
      real(kind=r8), intent(out) :: pbbfu(:,:)
      real(kind=r8), intent(out) :: pbbfddir(ncol,nlayers+2)
      real(kind=r8), intent(out) :: pbbcddir(ncol,nlayers+2)

      real(kind=r8), intent(out) :: puvcd(ncol,nlayers+2)
      real(kind=r8), intent(out) :: puvfd(ncol,nlayers+2)
      real(kind=r8), intent(out) :: puvcddir(ncol,nlayers+2)
      real(kind=r8), intent(out) :: puvfddir(:,:)

      real(kind=r8), intent(out) :: pnicd(ncol,nlayers+2)
      real(kind=r8), intent(out) :: pnifd(ncol,nlayers+2)
      real(kind=r8), intent(out) :: pnicddir(ncol,nlayers+2)
      real(kind=r8), intent(out) :: pnifddir(:,:)

! Added for net near-IR flux diagnostic
      real(kind=r8), intent(out) :: pnicu(ncol,nlayers+2)
      real(kind=r8), intent(out) :: pnifu(ncol,nlayers+2)

! Output - inactive                                              !   All Dimensions: (nlayers+1)
!      real(kind=r8), intent(out) :: puvcu(:)
!      real(kind=r8), intent(out) :: puvfu(:)
!      real(kind=r8), intent(out) :: pvscd(:)
!      real(kind=r8), intent(out) :: pvscu(:)
!      real(kind=r8), intent(out) :: pvsfd(:)
!      real(kind=r8), intent(out) :: pvsfu(:)

      real(kind=r8), intent(out)  :: pbbfsu(:,:,:)                 ! shortwave spectral flux up (nswbands,nlayers+1)
      real(kind=r8), intent(out)  :: pbbfsd(:,:,:)                 ! shortwave spectral flux down (nswbands,nlayers+1)


! ------- Local -------

      logical :: lrtchkclr(ncol,nlayers),lrtchkcld(ncol,nlayers)

      integer  :: klev
      integer :: ib1, ib2, ibm, igt, ikl, ikp, ikx
      integer :: jb, jg, jl, jk
      integer :: iw(ncol), iplon
!      integer, parameter :: nuv = ?? 
!      integer, parameter :: nvs = ?? 
      integer :: itind(ncol)

      real(kind=r8) :: tblind(ncol), ze1(ncol)
      real(kind=r8) :: zclear(ncol), zcloud(ncol)
      real(kind=r8) :: zdbt(ncol,nlayers+1), zdbt_nodel(ncol,nlayers+1)
      real(kind=r8) :: zgcc(ncol,nlayers), zgco(ncol,nlayers)
      real(kind=r8) :: zomcc(ncol,nlayers), zomco(ncol,nlayers)
      real(kind=r8) :: zrdnd(ncol,nlayers+1), zrdndc(ncol,nlayers+1)
      real(kind=r8) :: zref(ncol,nlayers+1), zrefc(ncol,nlayers+1), zrefo(ncol,nlayers+1)
      real(kind=r8) :: zrefd(ncol,nlayers+1), zrefdc(ncol,nlayers+1), zrefdo(ncol,nlayers+1)
      real(kind=r8) :: zrup(ncol,nlayers+1), zrupd(ncol,nlayers+1)
      real(kind=r8) :: zrupc(ncol,nlayers+1), zrupdc(ncol,nlayers+1)
      real(kind=r8) :: ztauc(ncol,nlayers), ztauo(ncol,nlayers)
      real(kind=r8) :: ztdbt(ncol,nlayers+1)
      real(kind=r8) :: ztra(ncol,nlayers+1), ztrac(ncol,nlayers+1), ztrao(ncol,nlayers+1)
      real(kind=r8) :: ztrad(ncol,nlayers+1), ztradc(ncol,nlayers+1), ztrado(ncol,nlayers+1)
      real(kind=r8) :: zdbtc(ncol,nlayers+1), ztdbtc(ncol,nlayers+1)
      real(kind=r8) :: zincflx(ncol,ngptsw), zdbtc_nodel(ncol,nlayers+1)
      real(kind=r8) :: ztdbt_nodel(ncol,nlayers+1), ztdbtc_nodel(ncol,nlayers+1)

      real(kind=r8) :: zdbtmc(ncol), zdbtmo(ncol), zf
      real(kind=r8) :: zwf, tauorig(ncol), repclc
!     real(kind=r8) :: zincflux                                   ! inactive

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=r8) :: ztaug(nlayers,16), ztaur(nlayers,16)
!      real(kind=r8) :: zsflxzen(16)
      real(kind=r8) :: ztaug(ncol,nlayers,ngptsw), ztaur(ncol,nlayers,ngptsw)
      real(kind=r8) :: zsflxzen(ncol,ngptsw)

! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=r8) :: zcd(ncol,nlayers+1,ngptsw), zcu(ncol,nlayers+1,ngptsw)
      real(kind=r8) :: zfd(ncol,nlayers+1,ngptsw), zfu(ncol,nlayers+1,ngptsw)

! Inactive arrays
!     real(kind=r8) :: zbbcd(nlayers+1), zbbcu(nlayers+1)
!     real(kind=r8) :: zbbfd(nlayers+1), zbbfu(nlayers+1)
!     real(kind=r8) :: zbbfddir(nlayers+1), zbbcddir(nlayers+1)
! ------------------------------------------------------------------

! Initializations

      ib1 = istart
      ib2 = iend
      klev = nlayers


!      zincflux = 0.0_r8

      repclc = 1.e-12_r8
      pbbcd(1:ncol,1:klev+1)=0._r8
      pbbcu(1:ncol,1:klev+1)=0._r8
      pbbfd(1:ncol,1:klev+1)=0._r8
      pbbfu(1:ncol,1:klev+1)=0._r8
      pbbcddir(1:ncol,1:klev+1)=0._r8
      pbbfddir(1:ncol,1:klev+1)=0._r8
      puvcd(1:ncol,1:klev+1)=0._r8
      puvfd(1:ncol,1:klev+1)=0._r8
      puvcddir(1:ncol,1:klev+1)=0._r8
      puvfddir(1:ncol,1:klev+1)=0._r8
      pnicd(1:ncol,1:klev+1)=0._r8
      pnifd(1:ncol,1:klev+1)=0._r8
      pnicddir(1:ncol,1:klev+1)=0._r8
      pnifddir(1:ncol,1:klev+1)=0._r8
      pnicu(1:ncol,1:klev+1)=0._r8
      pnifu(1:ncol,1:klev+1)=0._r8
      pbbfsu(:,1:ncol,1:klev+1)= 0._r8
      pbbfsd(:,1:ncol,1:klev+1)=0._r8


! Calculate the optical depths for gaseous absorption and Rayleigh scattering

      do iplon=1,ncol
         call taumol_sw(klev, &
                     colh2o(iplon,:), colco2(iplon,:), colch4(iplon,:),&
                     colo2(iplon,:), colo3(iplon,:), colmol(iplon,:), &
                     laytrop(iplon), jp(iplon,:), jt(iplon,:), jt1(iplon,:), &
                     fac00(iplon,:), fac01(iplon,:), fac10(iplon,:),&
                     fac11(iplon,:), &
                     selffac(iplon,:), selffrac(iplon,:),&
                     indself(iplon,:), forfac(iplon,:), forfrac(iplon,:), indfor(iplon,:), &
                     zsflxzen(iplon,:), ztaug(iplon,:,:),&
                     ztaur(iplon,:,:))
      enddo

! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

      jb = ib1-1                  ! ???
      do iplon=1,ncol
         iw(iplon) =0
      end do
      do iplon=1,ncol
! Clear-sky    
!   TOA direct beam    
!   Surface values
         ztdbtc(iplon,1)=1.0_r8
         ztdbtc_nodel(iplon,1)=1.0_r8
         zdbtc(iplon,klev+1) =0.0_r8
         ztrac(iplon,klev+1) =0.0_r8
         ztradc(iplon,klev+1)=0.0_r8
! Cloudy-sky
!   Surface values
         ztrao(iplon,klev+1) =0.0_r8
         ztrado(iplon,klev+1)=0.0_r8
! Total sky    
!   TOA direct beam    
         ztdbt(iplon,1)=1.0_r8
         ztdbt_nodel(iplon,1)=1.0_r8
!   Surface values
         zdbt(iplon,klev+1) =0.0_r8
         ztra(iplon,klev+1) =0.0_r8
         ztrad(iplon,klev+1)=0.0_r8
      enddo
      do jb = ib1, ib2
         ibm = jb-15
         igt = ngc(ibm)

! Reinitialize g-point counter for each band if output for each band is requested.
         do iplon=1,ncol
            if (iout.gt.0.and.ibm.ge.2) iw(iplon)= ngs(ibm-1)
         enddo


!        do jk=1,klev+1
!           zbbcd(jk)=0.0_r8
!           zbbcu(jk)=0.0_r8
!           zbbfd(jk)=0.0_r8
!           zbbfu(jk)=0.0_r8
!        enddo

! Top of g-point interval loop within each band (iw is cumulative counter) 
         do jg = 1,igt
            do iplon=1,ncol
              iw(iplon) = iw(iplon)+1

! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux
              zincflx(iplon,iw(iplon)) = adjflux(iplon,jb) * zsflxzen(iplon,iw(iplon)) * prmu0(iplon)
!             zincflux = zincflux + adjflux(jb) * zsflxzen(iw) * prmu0           ! inactive
            enddo

! Compute layer reflectances and transmittances for direct and diffuse sources, 
! first clear then cloudy

! zrefc(jk)  direct albedo for clear
! zrefo(jk)  direct albedo for cloud
! zrefdc(jk) diffuse albedo for clear
! zrefdo(jk) diffuse albedo for cloud
! ztrac(jk)  direct transmittance for clear
! ztrao(jk)  direct transmittance for cloudy
! ztradc(jk) diffuse transmittance for clear
! ztrado(jk) diffuse transmittance for cloudy
!  
! zref(jk)   direct reflectance
! zrefd(jk)  diffuse reflectance
! ztra(jk)   direct transmittance
! ztrad(jk)  diffuse transmittance
!
! zdbtc(jk)  clear direct beam transmittance
! zdbto(jk)  cloudy direct beam transmittance
! zdbt(jk)   layer mean direct beam transmittance
! ztdbt(jk)  total direct beam transmittance at levels
         do iplon=1,ncol
            zrefc(iplon,klev+1) =palbp(iplon,ibm)
            zrefdc(iplon,klev+1)=palbd(iplon,ibm)
            zrupc(iplon,klev+1) =palbp(iplon,ibm)
            zrupdc(iplon,klev+1)=palbd(iplon,ibm)
            zrefo(iplon,klev+1) =palbp(iplon,ibm)
            zrefdo(iplon,klev+1)=palbd(iplon,ibm)
            zref(iplon,klev+1) =palbp(iplon,ibm)
            zrefd(iplon,klev+1)=palbd(iplon,ibm)
            zrup(iplon,klev+1) =palbp(iplon,ibm)
            zrupd(iplon,klev+1)=palbd(iplon,ibm)
         enddo
! Top of layer loop
            do jk=1,klev
               ikl=klev+1-jk
               do iplon=1,ncol
! Note: two-stream calculations proceed from top to bottom; 
!   RRTMG_SW quantities are given bottom to top and are reversed here


! Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
                  lrtchkclr(iplon,jk)=.true.

!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
                  lrtchkcld(iplon,jk)=.false.
                  lrtchkcld(iplon,jk)=(pcldfmc(iplon,ikl,iw(iplon)) > repclc)

! Clear-sky optical parameters - this section inactive     
!   Original
!               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw)
!               zomcc(jk) = ztaur(ikl,iw) / ztauc(jk)
!               zgcc(jk) = 0.0001_r8
!   Total sky optical parameters        
!               ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaucmc(ikl,iw)
!               zomco(jk) = ptaucmc(ikl,iw) * pomgcmc(ikl,iw) + ztaur(ikl,iw)
!               zgco (jk) = (ptaucmc(ikl,iw) * pomgcmc(ikl,iw) * pasycmc(ikl,iw) + &
!                           ztaur(ikl,iw) * 0.0001_r8) / zomco(jk)
!               zomco(jk) = zomco(jk) / ztauo(jk)

! Clear-sky optical parameters including aerosols
                  if (ztaug(iplon,ikl,iw(iplon)) .lt. 0.0_r8) ztaug(iplon,ikl,iw(iplon)) = 0.0_r8

                  ztauc(iplon,jk) = ztaur(iplon,ikl,iw(iplon)) + ztaug(iplon,ikl,iw(iplon)) + ptaua(iplon,ikl,ibm)
                  zomcc(iplon,jk) = ztaur(iplon,ikl,iw(iplon)) * 1.0_r8 + ptaua(iplon,ikl,ibm) * pomga(iplon,ikl,ibm)
                  zgcc(iplon,jk) = pasya(iplon,ikl,ibm) * pomga(iplon,ikl,ibm) * ptaua(iplon,ikl,ibm) / zomcc(iplon,jk)
                  zomcc(iplon,jk) = zomcc(iplon,jk) / ztauc(iplon,jk)

! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
!   \/\/\/ This block of code is only needed for unscaled direct beam calculation
               enddo   
               if (idelm .eq. 0) then
                  do iplon=1,ncol
!     
                     zclear(iplon) = 1.0_r8 - pcldfmc(iplon,ikl,iw(iplon))
                     zcloud(iplon) = pcldfmc(iplon,ikl,iw(iplon))

! Clear
!                   zdbtmc = exp(-ztauc(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                     ze1(iplon) = ztauc(iplon,jk) / prmu0(iplon)
                  enddo
                  do iplon=1,ncol
                     if (ze1(iplon) .le. od_lo) then
                        zdbtmc(iplon) = 1._r8 - ze1(iplon) + 0.5_r8 * ze1(iplon) * ze1(iplon)
                     else 
                        tblind(iplon) = ze1(iplon) / (bpade + ze1(iplon))
                        itind(iplon) = tblint * tblind(iplon) + 0.5_r8
                        zdbtmc(iplon) = exp_tbl(itind(iplon))
                     endif
                  enddo
                  do iplon=1,ncol

                     zdbtc_nodel(iplon,jk) = zdbtmc(iplon)
                     ztdbtc_nodel(iplon,jk+1) = zdbtc_nodel(iplon,jk) * ztdbtc_nodel(iplon,jk)

! Clear + Cloud
                     tauorig(iplon) = ztauc(iplon,jk) + ptaormc(iplon,ikl,iw(iplon))
!                   zdbtmo = exp(-tauorig / prmu0)

! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                     ze1(iplon) = tauorig(iplon) / prmu0(iplon)
                  enddo
                  do iplon=1,ncol
                     if (ze1(iplon) .le. od_lo) then
                        zdbtmo(iplon) = 1._r8 - ze1(iplon) + 0.5_r8 * ze1(iplon) * ze1(iplon)
                     else
                        tblind(iplon) = ze1(iplon) / (bpade + ze1(iplon))
                        itind(iplon) = tblint * tblind(iplon) + 0.5_r8
                        zdbtmo(iplon) = exp_tbl(itind(iplon))
                     endif
                  enddo
                  do iplon=1,ncol

                     zdbt_nodel(iplon,jk) = zclear(iplon)*zdbtmc(iplon) + zcloud(iplon)*zdbtmo(iplon)
                     ztdbt_nodel(iplon,jk+1) = zdbt_nodel(iplon,jk) * ztdbt_nodel(iplon,jk)

                  enddo
               endif
               do iplon=1,ncol
!   /\/\/\ Above code only needed for unscaled direct beam calculation


! Delta scaling - clear   
                  zf = zgcc(iplon,jk) * zgcc(iplon,jk)
                  zwf = zomcc(iplon,jk) * zf
                  ztauc(iplon,jk) = (1.0_r8 - zwf) * ztauc(iplon,jk)
                  zomcc(iplon,jk) = (zomcc(iplon,jk) - zwf) / (1.0_r8 - zwf)
                  zgcc (iplon,jk) = (zgcc(iplon,jk) - zf) / (1.0_r8 - zf)
               enddo
! Total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
               if (icpr .ge. 1) then
                  do iplon=1,ncol
                     ztauo(iplon,jk) = ztauc(iplon,jk) + ptaucmc(iplon,ikl,iw(iplon))
                     zomco(iplon,jk) = ztauc(iplon,jk) * zomcc(iplon,jk) + ptaucmc(iplon,ikl,iw(iplon)) * pomgcmc(iplon,ikl,iw(iplon)) 
                     zgco (iplon,jk) = (ptaucmc(iplon,ikl,iw(iplon)) * pomgcmc(iplon,ikl,iw(iplon)) * pasycmc(iplon,ikl,iw(iplon)) + &
                                       ztauc(iplon,jk) * zomcc(iplon,jk) * zgcc(iplon,jk)) / zomco(iplon,jk)
                     zomco(iplon,jk) = zomco(iplon,jk) / ztauo(iplon,jk)
                  enddo
! Total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
               elseif (icpr .eq. 0) then
                  do iplon=1,ncol
                     ztauo(iplon,jk) = ztaur(iplon,ikl,iw(iplon)) + ztaug(iplon,ikl,iw(iplon)) + ptaua(iplon,ikl,ibm) + ptaucmc(iplon,ikl,iw(iplon))
                     zomco(iplon,jk) = ptaua(iplon,ikl,ibm) * pomga(iplon,ikl,ibm) + ptaucmc(iplon,ikl,iw(iplon)) * pomgcmc(iplon,ikl,iw(iplon)) + &
                                       ztaur(iplon,ikl,iw(iplon)) * 1.0_r8
                     zgco (iplon,jk) = (ptaucmc(iplon,ikl,iw(iplon)) * pomgcmc(iplon,ikl,iw(iplon)) * pasycmc(iplon,ikl,iw(iplon)) + &
                                       ptaua(iplon,ikl,ibm)*pomga(iplon,ikl,ibm)*pasya(iplon,ikl,ibm)) / zomco(iplon,jk)
                     zomco(iplon,jk) = zomco(iplon,jk) / ztauo(iplon,jk)

! Delta scaling - clouds 
!   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
                     zf = zgco(iplon,jk) * zgco(iplon,jk)
                     zwf = zomco(iplon,jk) * zf
                     ztauo(iplon,jk) = (1._r8 - zwf) * ztauo(iplon,jk)
                     zomco(iplon,jk) = (zomco(iplon,jk) - zwf) / (1.0_r8 - zwf)
                     zgco (iplon,jk) = (zgco(iplon,jk) - zf) / (1.0_r8 - zf)
                  enddo
               endif 

! End of layer loop
            enddo

! Clear sky reflectivities
            call reftra_sw (klev,ncol, &
                            lrtchkclr, zgcc, prmu0, ztauc, zomcc, &
                            zrefc, zrefdc, ztrac, ztradc)

! Total sky reflectivities      
            call reftra_sw (klev, ncol, &
                            lrtchkcld, zgco, prmu0, ztauo, zomco, &
                            zrefo, zrefdo, ztrao, ztrado)

            do jk=1,klev
               do iplon=1,ncol
! Combine clear and cloudy contributions for total sky
                  ikl = klev+1-jk 
                  zclear(iplon) = 1.0_r8 - pcldfmc(iplon,ikl,iw(iplon))
                  zcloud(iplon) = pcldfmc(iplon,ikl,iw(iplon))

                  zref(iplon,jk) = zclear(iplon)*zrefc(iplon,jk) + zcloud(iplon)*zrefo(iplon,jk)
                  zrefd(iplon,jk)= zclear(iplon)*zrefdc(iplon,jk) + zcloud(iplon)*zrefdo(iplon,jk)
                  ztra(iplon,jk) = zclear(iplon)*ztrac(iplon,jk) + zcloud(iplon)*ztrao(iplon,jk)
                  ztrad(iplon,jk)= zclear(iplon)*ztradc(iplon,jk) + zcloud(iplon)*ztrado(iplon,jk)

! Direct beam transmittance        

! Clear
!                zdbtmc(iplon) = exp(-ztauc(iplon,jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
                  ze1(iplon) = ztauc(iplon,jk) / prmu0(iplon)
               enddo
               do iplon=1,ncol
                  if (ze1(iplon) .le. od_lo) then
                     zdbtmc(iplon) = 1._r8 - ze1(iplon) + 0.5_r8 * ze1(iplon) * ze1(iplon)
                  else
                     tblind(iplon) = ze1(iplon) / (bpade + ze1(iplon))
                     itind(iplon) = tblint * tblind(iplon) + 0.5_r8
                     zdbtmc(iplon) = exp_tbl(itind(iplon))
                  endif

                  zdbtc(iplon,jk) = zdbtmc(iplon)
                  ztdbtc(iplon,jk+1) = zdbtc(iplon,jk)*ztdbtc(iplon,jk)

! Clear + Cloud
!                zdbtmo = exp(-ztauo(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
                  ze1(iplon) = ztauo(iplon,jk) / prmu0(iplon)
                  if (ze1(iplon) .le. od_lo) then
                     zdbtmo(iplon) = 1._r8 - ze1(iplon) + 0.5_r8 * ze1(iplon) * ze1(iplon)
                  else
                     tblind(iplon) = ze1(iplon) / (bpade + ze1(iplon))
                     itind(iplon) = tblint * tblind(iplon) + 0.5_r8
                     zdbtmo(iplon) = exp_tbl(itind(iplon))
                  endif

                  zdbt(iplon,jk) = zclear(iplon)*zdbtmc(iplon) + zcloud(iplon)*zdbtmo(iplon)
                  ztdbt(iplon,jk+1) = zdbt(iplon,jk)*ztdbt(iplon,jk)
        
               enddo           
            enddo     
! Vertical quadrature for clear-sky fluxes

            call vrtqdr_sw(ncol,klev, iw, &
                           zrefc, zrefdc, ztrac, ztradc, &
                           zdbtc, zrdndc, zrupc, zrupdc, ztdbtc, &
                           zcd, zcu)
      
! Vertical quadrature for cloudy fluxes

            call vrtqdr_sw(ncol,klev, iw, &
                           zref, zrefd, ztra, ztrad, &
                           zdbt, zrdnd, zrup, zrupd, ztdbt, &
                           zfd, zfu)

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays

            do jk=1,klev+1
               do iplon=1,ncol
                  ikl=klev+2-jk

! Accumulate spectral fluxes over bands - inactive
!               zbbfu(ikl) = zbbfu(ikl) + zincflx(iw)*zfu(jk,iw)  
!               zbbfd(ikl) = zbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
!               zbbcu(ikl) = zbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
!               zbbcd(ikl) = zbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
!               zbbfddir(ikl) = zbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!               zbbcddir(ikl) = zbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)

                  pbbfsu(ibm,iplon,ikl) = pbbfsu(ibm,iplon,ikl) + zincflx(iplon,iw(iplon))*zfu(iplon,jk,iw(iplon))
                  pbbfsd(ibm,iplon,ikl) = pbbfsd(ibm,iplon,ikl) + zincflx(iplon,iw(iplon))*zfd(iplon,jk,iw(iplon))

! Accumulate spectral fluxes over whole spectrum  
                  pbbfu(iplon,ikl) = pbbfu(iplon,ikl) + zincflx(iplon,iw(iplon))*zfu(iplon,jk,iw(iplon))
                  pbbfd(iplon,ikl) = pbbfd(iplon,ikl) + zincflx(iplon,iw(iplon))*zfd(iplon,jk,iw(iplon))
                  pbbcu(iplon,ikl) = pbbcu(iplon,ikl) + zincflx(iplon,iw(iplon))*zcu(iplon,jk,iw(iplon))
                  pbbcd(iplon,ikl) = pbbcd(iplon,ikl) + zincflx(iplon,iw(iplon))*zcd(iplon,jk,iw(iplon))
               enddo
               if (idelm .eq. 0) then
                  do iplon=1,ncol
                     pbbfddir(iplon,ikl) = pbbfddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbt_nodel(iplon,jk)
                     pbbcddir(iplon,ikl) = pbbcddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbtc_nodel(iplon,jk)
                  enddo
               elseif (idelm .eq. 1) then
                  do iplon=1,ncol
                     pbbfddir(iplon,ikl) = pbbfddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbt(iplon,jk)
                     pbbcddir(iplon,ikl) = pbbcddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbtc(iplon,jk)
                  enddo
               endif

! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  do iplon=1,ncol
                     puvcd(iplon,ikl) = puvcd(iplon,ikl) + zincflx(iplon,iw(iplon))*zcd(iplon,jk,iw(iplon))
                     puvfd(iplon,ikl) = puvfd(iplon,ikl) + zincflx(iplon,iw(iplon))*zfd(iplon,jk,iw(iplon))
                  enddo
                  if (idelm .eq. 0) then
                     do iplon=1,ncol
                        puvfddir(iplon,ikl) = puvfddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbt_nodel(iplon,jk)
                        puvcddir(iplon,ikl) = puvcddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbtc_nodel(iplon,jk)
                     enddo
                  elseif (idelm .eq. 1) then
                     do iplon=1,ncol
                        puvfddir(iplon,ikl) = puvfddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbt(iplon,jk)
                        puvcddir(iplon,ikl) = puvcddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbtc(iplon,jk)
                     enddo
                  endif
! band 9 is half-NearIR and half-Visible
               else if (ibm == 9) then  
                  do iplon=1,ncol
                     puvcd(iplon,ikl) = puvcd(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*zcd(iplon,jk,iw(iplon))
                     puvfd(iplon,ikl) = puvfd(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*zfd(iplon,jk,iw(iplon))
                     pnicd(iplon,ikl) = pnicd(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*zcd(iplon,jk,iw(iplon))
                     pnifd(iplon,ikl) = pnifd(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*zfd(iplon,jk,iw(iplon))
                  enddo
                  if (idelm .eq. 0) then
                     do iplon=1,ncol
                        puvfddir(iplon,ikl) = puvfddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbt_nodel(iplon,jk)
                        puvcddir(iplon,ikl) = puvcddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbtc_nodel(iplon,jk)
                        pnifddir(iplon,ikl) = pnifddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbt_nodel(iplon,jk)
                        pnicddir(iplon,ikl) = pnicddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbtc_nodel(iplon,jk)
                     enddo
                  elseif (idelm .eq. 1) then
                     do iplon=1,ncol
                        puvfddir(iplon,ikl) = puvfddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbt(iplon,jk)
                        puvcddir(iplon,ikl) = puvcddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbtc(iplon,jk)
                        pnifddir(iplon,ikl) = pnifddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbt(iplon,jk)
                        pnicddir(iplon,ikl) = pnicddir(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*ztdbtc(iplon,jk)
                     enddo
                  endif
                  do iplon=1,ncol
                     pnicu(iplon,ikl) = pnicu(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*zcu(iplon,jk,iw(iplon))
                     pnifu(iplon,ikl) = pnifu(iplon,ikl) + 0.5_r8*zincflx(iplon,iw(iplon))*zfu(iplon,jk,iw(iplon))
! Accumulate direct fluxes for near-IR bands
                  enddo
               else if (ibm == 14 .or. ibm <= 8) then  
                  do iplon=1,ncol
                     pnicd(iplon,ikl) = pnicd(iplon,ikl) + zincflx(iplon,iw(iplon))*zcd(iplon,jk,iw(iplon))
                     pnifd(iplon,ikl) = pnifd(iplon,ikl) + zincflx(iplon,iw(iplon))*zfd(iplon,jk,iw(iplon))
                  enddo
                  if (idelm .eq. 0) then
                     do iplon=1,ncol
                        pnifddir(iplon,ikl) = pnifddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbt_nodel(iplon,jk)
                        pnicddir(iplon,ikl) = pnicddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbtc_nodel(iplon,jk)
                     enddo
                  elseif (idelm .eq. 1) then
                     do iplon=1,ncol
                        pnifddir(iplon,ikl) = pnifddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbt(iplon,jk)
                        pnicddir(iplon,ikl) = pnicddir(iplon,ikl) + zincflx(iplon,iw(iplon))*ztdbtc(iplon,jk)
                     enddo
                  endif
! Added for net near-IR flux diagnostic
                  do iplon=1,ncol
                     pnicu(iplon,ikl) = pnicu(iplon,ikl) + zincflx(iplon,iw(iplon))*zcu(iplon,jk,iw(iplon))
                     pnifu(iplon,ikl) = pnifu(iplon,ikl) + zincflx(iplon,iw(iplon))*zfu(iplon,jk,iw(iplon))
                  enddo
               endif


! End loop on jg, g-point interval
            enddo             

! End loop on jb, spectral band
         enddo                    
      end do

      end subroutine spcvmc_sw

      end module rrtmg_sw_spcvmc


