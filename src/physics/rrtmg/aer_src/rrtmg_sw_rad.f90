!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/01/03 21:35:35 $
!

module rrtmg_sw_rad

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! --------- Modules ---------

use shr_kind_mod,        only: r8=>shr_kind_r8

use mcica_subcol_gen_sw, only: mcica_subcol_sw
use rrtmg_sw_cldprmc,    only: cldprmc_sw
use rrtmg_sw_setcoef,    only: setcoef_sw
use rrtmg_sw_spcvmc,     only: spcvmc_sw

implicit none

public :: rrtmg_sw

! CAM supplies shortwave cloud optical properties
integer, parameter :: inflag  = 0 ! flag for cloud parameterization method
integer, parameter :: iceflag = 0 ! flag for ice cloud parameterization method
integer, parameter :: liqflag = 0 ! flag for liquid cloud parameterization method

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 10, input total aerosol optical depth, single scattering albedo 
!            and asymmetry parameter (tauaer, ssaaer, asmaer) directly
integer, parameter :: iaer = 10

! Set idelm to select between delta-M scaled or unscaled output direct and diffuse fluxes
! NOTE: total downward fluxes are always delta scaled
! idelm = 0, output direct and diffuse flux components are not delta scaled
!            (direct flux does not include forward scattering peak)
! idelm = 1, output direct and diffuse flux components are delta scaled (default)
!            (direct flux includes part or most of forward scattering peak)
integer, parameter :: idelm = 1

!=========================================================================================
contains
!=========================================================================================

subroutine rrtmg_sw &
            (lchnk   ,ncol    ,nlay    ,icld    ,          &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,solvar, &
             cldfmcl ,taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl, &
             ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  ,ssaaer  ,asmaer  , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc, &
             swuflxs, swdflxs)


! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPRMC_SW to set cloud optical depth for McICA based
!       on input cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVMC to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer
!       using McICA, the Monte-Carlo Independent Column Approximation,
!       to represent sub-grid scale cloud variability
!    f) passes the calculated fluxes and cooling rates back to GCM
!
! *** This version uses McICA ***
!     Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!     JC, 2003) method is applied to the forward model calculation 
!     This method is valid for clear sky or partial cloud conditions.
!
! This call to RRTMG_SW must be preceeded by a call to the module
!     mcica_subcol_gen_sw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngptsw) dimension.
!
! *** This version only allows input of cloud optical properties ***
!     Input cloud fraction, cloud optical depth, single scattering albedo 
!     and asymmetry parameter directly (inflg = 0)
!
! *** This version only allows input of aerosol optical properties ***
!     Input aerosol optical depth, single scattering albedo and asymmetry
!     parameter directly by layer and spectral band (iaer=10)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs. 
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90 
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse 
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays. 
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to output direct and diffuse fluxes either with or without
!   delta scaling based on setting of idelm flag
!     Dec 2008: M. J. Iacono, AER, Inc.

   use parrrsw,  only: nbndsw, ngptsw, mxmol, jpband, jpb1, jpb2
   use rrsw_con, only: heatfac, oneminus, pi


   ! ----- Input -----
   integer, intent(in) :: lchnk                      ! chunk identifier
   integer, intent(in) :: ncol                       ! Number of horizontal columns     
   integer, intent(in) :: nlay                       ! Number of model layers
   integer, intent(in) :: icld                       ! Cloud overlap method
                                                     !    0: Clear only
                                                     !    1: Random
                                                     !    2: Maximum/random
                                                     !    3: Maximum
   real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                     !    Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                     !    Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                     !    Dimensions: (ncol)
   real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: asdir(:)             ! UV/vis surface albedo direct rad
                                                     !    Dimensions: (ncol)
   real(kind=r8), intent(in) :: aldir(:)             ! Near-IR surface albedo direct rad
                                                     !    Dimensions: (ncol)
   real(kind=r8), intent(in) :: asdif(:)             ! UV/vis surface albedo: diffuse rad
                                                     !    Dimensions: (ncol)
   real(kind=r8), intent(in) :: aldif(:)             ! Near-IR surface albedo: diffuse rad
                                                     !    Dimensions: (ncol)

   integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                     !  distance if adjflx not provided)
   real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
   real(kind=r8), intent(in) :: coszen(:)            ! Cosine of solar zenith angle
                                                     !    Dimensions: (ncol)
   real(kind=r8), intent(in) :: solvar(1:nbndsw)     ! Solar constant (Wm-2) scaling per band

   real(kind=r8), intent(in) :: cldfmcl(:,:,:)       ! Cloud fraction
                                                     !    Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: taucmcl(:,:,:)       ! Cloud optical depth
                                                     !    Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: ssacmcl(:,:,:)       ! Cloud single scattering albedo
                                                     !    Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: asmcmcl(:,:,:)       ! Cloud asymmetry parameter
                                                     !    Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: fsfcmcl(:,:,:)       ! Cloud forward scattering parameter
                                                     !    Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                                     !    Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                                     !    Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth (iaer=10 only)
                                                     !    Dimensions: (ncol,nlay,nbndsw)
                                                     ! (non-delta scaled)      
   real(kind=r8), intent(in) :: ssaaer(:,:,:)        ! Aerosol single scattering albedo (iaer=10 only)
                                                     !    Dimensions: (ncol,nlay,nbndsw)
                                                     ! (non-delta scaled)      
   real(kind=r8), intent(in) :: asmaer(:,:,:)        ! Aerosol asymmetry parameter (iaer=10 only)
                                                     !    Dimensions: (ncol,nlay,nbndsw)
                                                     ! (non-delta scaled)      

   ! ----- Output -----

   real(kind=r8), intent(out) :: swuflx(:,:)         ! Total sky shortwave upward flux (W/m2)
                                                     !    Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(out) :: swdflx(:,:)         ! Total sky shortwave downward flux (W/m2)
                                                     !    Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(out) :: swhr(:,:)           ! Total sky shortwave radiative heating rate (K/d)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: swuflxc(:,:)        ! Clear sky shortwave upward flux (W/m2)
                                                     !    Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(out) :: swdflxc(:,:)        ! Clear sky shortwave downward flux (W/m2)
                                                     !    Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(out) :: swhrc(:,:)          ! Clear sky shortwave radiative heating rate (K/d)
                                                     !    Dimensions: (ncol,nlay)

   real(kind=r8), intent(out) :: dirdnuv(:,:)        ! Direct downward shortwave flux, UV/vis
   real(kind=r8), intent(out) :: difdnuv(:,:)        ! Diffuse downward shortwave flux, UV/vis
   real(kind=r8), intent(out) :: dirdnir(:,:)        ! Direct downward shortwave flux, near-IR
   real(kind=r8), intent(out) :: difdnir(:,:)        ! Diffuse downward shortwave flux, near-IR

   real(kind=r8), intent(out) :: ninflx(:,:)         ! Net shortwave flux, near-IR
   real(kind=r8), intent(out) :: ninflxc(:,:)        ! Net clear sky shortwave flux, near-IR

   real(kind=r8), intent(out)  :: swuflxs(:,:,:)     ! shortwave spectral flux up
   real(kind=r8), intent(out)  :: swdflxs(:,:,:)     ! shortwave spectral flux down

   ! ----- Local -----

   ! Control
   integer :: istart                         ! beginning band of calculation
   integer :: iend                           ! ending band of calculation
   integer :: icpr                           ! cldprop/cldprmc use flag
   integer :: iout = 0                       ! output option flag (inactive)
   integer :: isccos                         ! instrumental cosine response flag (inactive)
   integer :: iplon                          ! column loop index
   integer :: i                              ! layer loop index                       ! jk
   integer :: ib                             ! band loop index                        ! jsw
   integer :: ia, ig                         ! indices
   integer :: k                              ! layer loop index
   integer :: ims                            ! value for changing mcica permute seed

   real(kind=r8) :: zepsec, zepzen           ! epsilon
   real(kind=r8) :: zdpgcp                   ! flux to heating conversion ratio

   ! Atmosphere
   real(kind=r8) :: pavel(ncol,nlay)            ! layer pressures (mb)
   real(kind=r8) :: tavel(ncol,nlay)            ! layer temperatures (K)
   real(kind=r8) :: pz(ncol,0:nlay)             ! level (interface) pressures (hPa, mb)
   real(kind=r8) :: tz(ncol,0:nlay)             ! level (interface) temperatures (K)
   real(kind=r8) :: tbound(ncol)                   ! surface temperature (K)
   real(kind=r8) :: pdp(ncol,nlay)              ! layer pressure thickness (hPa, mb)
   real(kind=r8) :: coldry(ncol,nlay)           ! dry air column amount
   real(kind=r8) :: wkl(ncol,mxmol,nlay)        ! molecular amounts (mol/cm-2)

   real(kind=r8) :: cossza(ncol)             ! Cosine of solar zenith angle
   real(kind=r8) :: adjflux(ncol,jpband)     ! adjustment for current Earth/Sun distance
                                                !  default value of 1368.22 Wm-2 at 1 AU
   real(kind=r8) :: albdir(ncol,nbndsw)      ! surface albedo, direct          ! zalbp
   real(kind=r8) :: albdif(ncol,nbndsw)      ! surface albedo, diffuse         ! zalbd

   ! Atmosphere - setcoef
   integer :: laytrop(ncol)                        ! tropopause layer index
   integer :: layswtch(ncol)                       !
   integer :: laylow(ncol)                         !
   integer :: jp(ncol,nlay)                     !
   integer :: jt(ncol,nlay)                     !
   integer :: jt1(ncol,nlay)                    !

   real(kind=r8) :: colh2o(ncol,nlay)           ! column amount (h2o)
   real(kind=r8) :: colco2(ncol,nlay)           ! column amount (co2)
   real(kind=r8) :: colo3(ncol,nlay)            ! column amount (o3)
   real(kind=r8) :: coln2o(ncol,nlay)           ! column amount (n2o)
   real(kind=r8) :: colch4(ncol,nlay)           ! column amount (ch4)
   real(kind=r8) :: colo2(ncol,nlay)            ! column amount (o2)
   real(kind=r8) :: colmol(ncol,nlay)           ! column amount
   real(kind=r8) :: co2mult(ncol,nlay)          ! column amount

   integer :: indself(ncol,nlay)
   integer :: indfor(ncol,nlay)
   real(kind=r8) :: selffac(ncol,nlay)
   real(kind=r8) :: selffrac(ncol,nlay)
   real(kind=r8) :: forfac(ncol,nlay)
   real(kind=r8) :: forfrac(ncol,nlay)

   real(kind=r8) :: fac00(ncol,nlay)
   real(kind=r8) :: fac01(ncol,nlay)
   real(kind=r8) :: fac11(ncol,nlay)
   real(kind=r8) :: fac10(ncol,nlay)

   ! Atmosphere/clouds - cldprmc [mcica]
   real(kind=r8) :: ciwpmc(ncol,ngptsw,nlay)    ! cloud ice water path [mcica]
   real(kind=r8) :: clwpmc(ncol,ngptsw,nlay)    ! cloud liquid water path [mcica]
   real(kind=r8) :: relqmc(ncol,nlay)           ! liquid particle size (microns)
   real(kind=r8) :: reicmc(ncol,nlay)           ! ice particle effective radius (microns)
   real(kind=r8) :: dgesmc(ncol,nlay)           ! ice particle generalized effective size (microns)
   real(kind=r8) :: fsfcmc(ncol,ngptsw,nlay)    ! cloud forward scattering fraction [mcica]

   ! Atmosphere/clouds/aerosol - spcvrt,spcvmc
   real(kind=r8) :: ztaua(ncol,nlay,nbndsw)     ! total aerosol optical depth
   real(kind=r8) :: zasya(ncol,nlay,nbndsw)     ! total aerosol asymmetry parameter
   real(kind=r8) :: zomga(ncol,nlay,nbndsw)     ! total aerosol single scattering albedo
   real(kind=r8) :: zcldfmc(ncol,nlay,ngptsw)   ! cloud fraction [mcica]
   real(kind=r8) :: ztaucmc(ncol,nlay,ngptsw)   ! cloud optical depth [mcica]
   real(kind=r8) :: ztaormc(ncol,nlay,ngptsw)   ! unscaled cloud optical depth [mcica]
   real(kind=r8) :: zasycmc(ncol,nlay,ngptsw)   ! cloud asymmetry parameter [mcica]
   real(kind=r8) :: zomgcmc(ncol,nlay,ngptsw)   ! cloud single scattering albedo [mcica]

   real(kind=r8) :: zbbfddir(ncol,nlay+2)       ! temporary downward direct shortwave flux (w/m2)
   real(kind=r8) :: zbbcddir(ncol,nlay+2)       ! temporary clear sky downward direct shortwave flux (w/m2)
   real(kind=r8) :: zuvfd(ncol,nlay+2)          ! temporary UV downward shortwave flux (w/m2)
   real(kind=r8) :: zuvcd(ncol,nlay+2)          ! temporary clear sky UV downward shortwave flux (w/m2)
   real(kind=r8) :: zuvcddir(ncol,nlay+2)       ! temporary clear sky UV downward direct shortwave flux (w/m2)
   real(kind=r8) :: znifd(ncol,nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
   real(kind=r8) :: znicd(ncol,nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)
   real(kind=r8) :: znicddir(ncol,nlay+2)       ! temporary clear sky near-IR downward direct shortwave flux (w/m2)

   ! Added for near-IR flux diagnostic
   real(kind=r8) :: znifu(ncol,nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
   real(kind=r8) :: znicu(ncol,nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)

   ! Optional output fields 
   real(kind=r8) :: swnflx(nlay+2)         ! Total sky shortwave net flux (W/m2)
   real(kind=r8) :: swnflxc(nlay+2)        ! Clear sky shortwave net flux (W/m2)
   real(kind=r8) :: dirdflux(nlay+2)       ! Direct downward shortwave surface flux
   real(kind=r8) :: difdflux(nlay+2)       ! Diffuse downward shortwave surface flux
   real(kind=r8) :: uvdflx(nlay+2)         ! Total sky downward shortwave flux, UV/vis   
   real(kind=r8) :: nidflx(nlay+2)         ! Total sky downward shortwave flux, near-IR  

   ! Initializations

   zepsec = 1.e-06_r8
   zepzen = 1.e-10_r8
   oneminus = 1.0_r8 - zepsec
   pi = 2._r8 * asin(1._r8)

   istart = jpb1
   iend = jpb2
   icpr = 0
   ims = 2

   ! Prepare atmosphere profile from GCM for use in RRTMG, and define
   ! other input parameters
   call inatm_sw (ncol,nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, &
              h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, adjes, dyofyr, solvar, &
              cldfmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl, ciwpmcl, clwpmcl, &
              reicmcl, relqmcl, tauaer, ssaaer, asmaer, &
              pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
              adjflux, zcldfmc, ztaucmc, &
              zomgcmc, zasycmc, fsfcmc, ciwpmc, clwpmc, reicmc, dgesmc, relqmc, &
              ztaua, zomga, zasya)

   !  Cloud fraction and cloud
   !  optical properties are transferred to rrtmg_sw arrays in cldprop.

   call cldprmc_sw(ncol,nlay, inflag, iceflag, liqflag, zcldfmc, &
                         ciwpmc, clwpmc, reicmc, dgesmc, relqmc, &
                         ztaormc, ztaucmc, zomgcmc, zasycmc, fsfcmc)
   icpr = 1

   ! This is the main longitude/column loop in RRTMG.
   ! Modify to loop over all columns (nlon) or over daylight columns

   do iplon = 1, ncol

      ! Calculate coefficients for the temperature and pressure dependence of the
      ! molecular absorption coefficients by interpolating data from stored
      ! reference atmospheres.

      call setcoef_sw(nlay, pavel(iplon,:), tavel(iplon,:), pz(iplon,:), &
                      tz(iplon,:), tbound(iplon), coldry(iplon,:), wkl(iplon,:,:), &
                      laytrop(iplon), layswtch(iplon), laylow(iplon), &
                      jp(iplon,:), jt(iplon,:), jt1(iplon,:), &
                      co2mult(iplon,:), colch4(iplon,:), colco2(iplon,:),&
                      colh2o(iplon,:), colmol(iplon,:), coln2o(iplon,:), &
                      colo2(iplon,:), colo3(iplon,:), fac00(iplon,:),&
                      fac01(iplon,:), fac10(iplon,:), fac11(iplon,:), &
                      selffac(iplon,:), selffrac(iplon,:), indself(iplon,:),&
                      forfac(iplon,:), forfrac(iplon,:), indfor(iplon,:))
   end do

   ! Cosine of the solar zenith angle
   ! Prevent using value of zero; ideally, SW model is not called from host model when sun
   ! is below horizon

   do iplon = 1, ncol
      cossza(iplon) = coszen(iplon)

      if (cossza(iplon) .lt. zepzen) cossza(iplon) = zepzen
   end do

   ! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer

   ! Surface albedo
   !  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
   !         do ib=1,9
   do ib=1,8
      do iplon = 1, ncol
         albdir(iplon,ib) = aldir(iplon)
         albdif(iplon,ib) = aldif(iplon)
      enddo
   enddo

   do iplon = 1, ncol
      albdir(iplon,nbndsw) = aldir(iplon)
      albdif(iplon,nbndsw) = aldif(iplon)
      !  Set band 24 (or, band 9 counting from 1) to use linear average of UV/visible
      !  and near-IR values, since this band straddles 0.7 microns:
      albdir(iplon,9) = 0.5*(aldir(iplon) + asdir(iplon))
      albdif(iplon,9) = 0.5*(aldif(iplon) + asdif(iplon))
   enddo

   !  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
   do ib=10,13
      do iplon = 1, ncol
         albdir(iplon,ib) = asdir(iplon)
         albdif(iplon,ib) = asdif(iplon)
      enddo
   enddo

   ! Clouds
   if (icld.eq.0) then
      do iplon = 1, ncol
         zcldfmc(iplon,1:nlay,1:ngptsw) = 0._r8
         ztaucmc(iplon,1:nlay,1:ngptsw) = 0._r8
         ztaormc(iplon,1:nlay,1:ngptsw) = 0._r8
         zasycmc(iplon,1:nlay,1:ngptsw) = 0._r8
         zomgcmc(iplon,1:nlay,1:ngptsw) = 1._r8
      enddo
   endif

   ! Aerosol
   ! IAER = 0: no aerosols
   if (iaer.eq.0) then
      do iplon = 1, ncol
         ztaua(iplon,:,:) = 0._r8
         zasya(iplon,:,:) = 0._r8
         zomga(iplon,:,:) = 1._r8
      enddo
   endif

   ! Call the 2-stream radiation transfer model

   call spcvmc_sw &
             (lchnk, ncol, nlay, istart, iend, icpr, idelm, iout, &
              pavel, tavel, pz, tz, tbound, albdif, albdir, &
              zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc, &
              ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &
              laytrop, layswtch, laylow, jp, jt, jt1, &
              co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
              fac00, fac01, fac10, fac11, &
              selffac, selffrac, indself, forfac, forfrac, indfor, &
              swdflx, swuflx, swdflxc, swuflxc, zuvfd, zuvcd, znifd, znicd, znifu, znicu, &
              zbbfddir, zbbcddir, dirdnuv, zuvcddir, dirdnir, znicddir, swuflxs, swdflxs)

   ! Transfer up and down, clear and total sky fluxes to output arrays.
   ! Vertical indexing goes from bottom to top

   do i = 1, nlay+1
      uvdflx(i) = zuvfd(ncol,i)
      nidflx(i) = znifd(ncol,i)

      do iplon = 1, ncol
         !  Direct/diffuse fluxes
         dirdflux(i) = zbbfddir(iplon,i)
         difdflux(i) = swdflx(iplon,i) - dirdflux(i)
         !  UV/visible direct/diffuse fluxes
         difdnuv(iplon,i) = zuvfd(iplon,i) - dirdnuv(iplon,i)
         !  Near-IR direct/diffuse fluxes
         difdnir(iplon,i) = znifd(iplon,i) - dirdnir(iplon,i)
         !  Added for net near-IR diagnostic
         ninflx(iplon,i) = znifd(iplon,i) - znifu(iplon,i)
         ninflxc(iplon,i) = znicd(iplon,i) - znicu(iplon,i)
      end do
   end do

   do iplon = 1, ncol
      !  Total and clear sky net fluxes
      do i = 1, nlay+1
         swnflxc(i) = swdflxc(iplon,i) - swuflxc(iplon,i)
         swnflx(i) = swdflx(iplon,i) - swuflx(iplon,i)
      end do

      !  Total and clear sky heating rates
      !  Heating units are in K/d. Flux units are in W/m2.
      do i = 1, nlay
         zdpgcp = heatfac / pdp(iplon,i)
         swhrc(iplon,i) = (swnflxc(i+1) - swnflxc(i)) * zdpgcp
         swhr(iplon,i) = (swnflx(i+1) - swnflx(i)) * zdpgcp
      end do
      swhrc(iplon,nlay) = 0._r8
      swhr(iplon,nlay) = 0._r8

   end do

end subroutine rrtmg_sw

!=========================================================================================

real(kind=r8) function earth_sun(idn)

   !  Purpose: Function to calculate the correction factor of Earth's orbit
   !  for current day of the year

   !  idn        : Day of the year
   !  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

   ! ------- Modules -------

   use rrsw_con, only : pi

   integer, intent(in) :: idn

   real(kind=r8) :: gamma

   gamma = 2._r8*pi*(idn-1)/365._r8

   ! Use Iqbal's equation 1.2.1

   earth_sun = 1.000110_r8 + .034221_r8 * cos(gamma) + .001289_r8 * sin(gamma) + &
                   .000719_r8 * cos(2._r8*gamma) + .000077_r8 * sin(2._r8*gamma)

end function earth_sun

!=========================================================================================

subroutine inatm_sw (ncol, nlay, icld, iaer, &
            play, plev, tlay, tlev, tsfc, &
            h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, adjes, dyofyr, solvar, &
            cldfmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl, ciwpmcl, clwpmcl, &
            reicmcl, relqmcl, tauaer, ssaaer, asmaer, &
            pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
            adjflux, zcldfmc, ztaucmc, &
            zssacmc, zasmcmc, fsfcmc, ciwpmc, clwpmc, reicmc, dgesmc, relqmc, &
            taua, ssaa, asma)

   ! Input atmospheric profile from GCM, and prepare it for use in RRTMG_SW.
   ! Set other RRTMG_SW input parameters.  

   use parrrsw,  only: nbndsw, ngptsw, nmol, mxmol, &
                       jpband, jpb1, jpb2
   use rrsw_con, only: grav, avogad

   ! ----- Input -----
   integer, intent(in) :: ncol                       ! column end index
   integer, intent(in) :: nlay                       ! number of model layers
   integer, intent(in) :: icld                       ! clear/cloud and cloud overlap flag
   integer, intent(in) :: iaer                       ! aerosol option flag

   real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                     ! Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                     ! Dimensions: (ncol,nlay+1)
   real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                     ! Dimensions: (ncol)
   real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                     ! Dimensions: (ncol,nlay)

   integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                     !  distance if adjflx not provided)
   real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
   real(kind=r8), intent(in) :: solvar(jpb1:jpb2)    ! Solar constant (Wm-2) scaling per band

   real(kind=r8), intent(in) :: cldfmcl(:,:,:)       ! Cloud fraction
                                                     ! Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: taucmcl(:,:,:)       ! Cloud optical depth (optional)
                                                     ! Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: ssacmcl(:,:,:)       ! Cloud single scattering albedo
                                                     ! Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: asmcmcl(:,:,:)       ! Cloud asymmetry parameter
                                                     ! Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: fsfcmcl(:,:,:)       ! Cloud forward scattering fraction
                                                     ! Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                                     ! Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                                     ! Dimensions: (ngptsw,ncol,nlay)
   real(kind=r8), intent(in) :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                     ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                     ! Dimensions: (ncol,nlay)

   real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth
                                                     ! Dimensions: (ncol,nlay,nbndsw)
   real(kind=r8), intent(in) :: ssaaer(:,:,:)        ! Aerosol single scattering albedo
                                                     ! Dimensions: (ncol,nlay,nbndsw)
   real(kind=r8), intent(in) :: asmaer(:,:,:)        ! Aerosol asymmetry parameter
                                                     ! Dimensions: (ncol,nlay,nbndsw)

   ! Atmosphere

   real(kind=r8), intent(out) :: pavel(ncol,nlay)      ! layer pressures (mb)
                                                       ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: tavel(ncol,nlay)      ! layer temperatures (K)
                                                       ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: pz(ncol,0:nlay)       ! level (interface) pressures (hPa, mb)
                                                       ! Dimensions: (ncol,0:nlay)
   real(kind=r8), intent(out) :: tz(ncol,0:nlay)       ! level (interface) temperatures (K)
                                                       ! Dimensions: (ncol,0:nlay)
   real(kind=r8), intent(out) :: tbound(ncol)          ! surface temperature (K)
                                                       ! Dimensions: (ncol)
   real(kind=r8), intent(out) :: pdp(ncol,nlay)        ! layer pressure thickness (hPa, mb)
                                                       ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: coldry(ncol,nlay)     ! dry air column density (mol/cm2)
                                                       ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: wkl(ncol,mxmol,nlay)  ! molecular amounts (mol/cm-2)
                                                       ! Dimensions: (ncol,mxmol,nlay)

   real(kind=r8), intent(out) :: adjflux(ncol,jpband)  ! adjustment for current Earth/Sun distance
                                                       ! Dimensions: (ncol,jpband)
   real(kind=r8), intent(out) :: taua(ncol,nlay,nbndsw) ! Aerosol optical depth
                                                        ! Dimensions: (ncol,nlay,nbndsw)
   real(kind=r8), intent(out) :: ssaa(ncol,nlay,nbndsw) ! Aerosol single scattering albedo
                                                        ! Dimensions: (ncol,nlay,nbndsw)
   real(kind=r8), intent(out) :: asma(ncol,nlay,nbndsw) ! Aerosol asymmetry parameter
                                                        ! Dimensions: (ncol,nlay,nbndsw)

   ! Atmosphere/clouds - cldprop

   real(kind=r8), intent(out) :: zcldfmc(ncol,nlay,ngptsw) ! layer cloud fraction
                                                           ! Dimensions: (ncol,nlay,ngptsw)
   real(kind=r8), intent(out) :: ztaucmc(ncol,nlay,ngptsw) ! cloud optical depth (non-delta scaled)
                                                           ! Dimensions: (ncol,nlay,ngptsw)
   real(kind=r8), intent(out) :: zssacmc(ncol,nlay,ngptsw) ! cloud single scattering albedo (non-delta-scaled)
                                                           ! Dimensions: (ncol,nlay,ngptsw)
   real(kind=r8), intent(out) :: zasmcmc(ncol,nlay,ngptsw) ! cloud asymmetry parameter (non-delta scaled)
   real(kind=r8), intent(out) :: fsfcmc(ncol,ngptsw,nlay) ! cloud forward scattering fraction (non-delta scaled)
                                                          ! Dimensions: (ncol,ngptsw,nlay)
   real(kind=r8), intent(out) :: ciwpmc(ncol,ngptsw,nlay) ! cloud ice water path
                                                          ! Dimensions: (ncol,ngptsw,nlay)
   real(kind=r8), intent(out) :: clwpmc(ncol,ngptsw,nlay) ! cloud liquid water path
                                                          ! Dimensions: (ncol,ngptsw,nlay)
   real(kind=r8), intent(out) :: reicmc(ncol,nlay)        ! cloud ice particle effective radius
                                                          ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: dgesmc(ncol,nlay)        ! cloud ice particle effective radius
                                                          ! Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: relqmc(ncol,nlay)        ! cloud liquid particle size
                                                          ! Dimensions: (ncol,nlay)

   ! ----- Local -----
   real(kind=r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
   real(kind=r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)

   ! Set molecular weight ratios (for converting mmr to vmr)
   !  e.g. h2ovmr = h2ommr * amdw)
   real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
   real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
   real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
   real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
   real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide

   real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)

   integer :: isp, l, ix, n, imol, ib, ig, iplon   ! Loop indices
   real(kind=r8) :: amm, summol                      ! 
   real(kind=r8) :: adjflx                           ! flux adjustment for Earth/Sun distance
   !-----------------------------------------------------------------------------------------

   ! Set flux adjustment for current Earth/Sun distance (two options).
   ! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
   adjflx = adjes

   ! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
   !    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 
   if (dyofyr .gt. 0) then
      adjflx = earth_sun(dyofyr)
   endif

   ! Set incoming solar flux adjustment to include adjustment for
   ! current Earth/Sun distance (ADJFLX) and scaling of default internal
   ! solar constant (rrsw_scon = 1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
   ! to a single scaling factor as needed, or to a different value in each 
   ! band, which may be necessary for paleoclimate simulations. 

   do iplon = 1 ,ncol
      adjflux(iplon,:) = 0._r8
   end do

   do ib = jpb1,jpb2
      do iplon = 1, ncol
         adjflux(iplon,ib) = adjflx * solvar(ib)
      end do
   end do

   do iplon = 1, ncol
      !  Set surface temperature.
      tbound(iplon) = tsfc(iplon)

      !  Install input GCM arrays into RRTMG_SW arrays for pressure, temperature,
      !  and molecular amounts.  
      !  Pressures are input in mb, or are converted to mb here.
      !  Molecular amounts are input in volume mixing ratio, or are converted from 
      !  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
      !  here. These are then converted to molecular amount (molec/cm2) below.  
      !  The dry air column COLDRY (in molec/cm2) is calculated from the level 
      !  pressures, pz (in mb), based on the hydrostatic equation and includes a 
      !  correction to account for h2o in the layer.  The molecular weight of moist 
      !  air (amm) is calculated for each layer.  
      !  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
      !  assumes GCM input fields are also bottom to top. Input layer indexing
      !  from GCM fields should be reversed here if necessary.
      pz(iplon,0) = plev(iplon,nlay+1)
      tz(iplon,0) = tlev(iplon,nlay+1)
   end do

   do l = 1, nlay
      do iplon = 1, ncol
         pavel(iplon,l) = play(iplon,nlay-l+1)
         tavel(iplon,l) = tlay(iplon,nlay-l+1)
         pz(iplon,l) = plev(iplon,nlay-l+1)
         tz(iplon,l) = tlev(iplon,nlay-l+1)
         pdp(iplon,l) = pz(iplon,l-1) - pz(iplon,l)
      end do
   end do

   do iplon = 1, ncol
      do l = 1, nlay

         ! For h2o input in vmr:
         wkl(iplon,1,l) = h2ovmr(iplon,nlay-l+1)
         wkl(iplon,2,l) = co2vmr(iplon,nlay-l+1)
         wkl(iplon,3,l) = o3vmr(iplon,nlay-l+1)
         wkl(iplon,4,l) = n2ovmr(iplon,nlay-l+1)
         wkl(iplon,5,l) = 0._r8
         wkl(iplon,6,l) = ch4vmr(iplon,nlay-l+1)
         wkl(iplon,7,l) = o2vmr(iplon,nlay-l+1) 
         amm = (1._r8 - wkl(iplon,1,l)) * amd + wkl(iplon,1,l) * amw            
         coldry(iplon,l) = (pz(iplon,l-1)-pz(iplon,l)) * 1.e3_r8 * avogad / &
            (1.e2_r8 * grav * amm * (1._r8 + wkl(iplon,1,l)))
      end do

      coldry(iplon,nlay) = (pz(iplon,nlay-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(iplon,1,nlay-1)))

      ! At this point all molecular amounts in wkl are in volume mixing ratio;
      ! convert to molec/cm2 based on coldry for use in rrtm.

      do l = 1, nlay
         do imol = 1, nmol
            wkl(iplon,imol,l) = coldry(iplon,l) * wkl(iplon,imol,l)
         end do
      end do
   end do

   ! Transfer aerosol optical properties to RRTM variables;
   ! modify to reverse layer indexing here if necessary.

   if (iaer .ge. 1) then 
      do l = 1, nlay-1
         do ib = 1, nbndsw
            do iplon = 1, ncol
               taua(iplon,l,ib) = tauaer(iplon,nlay-l,ib)
               ssaa(iplon,l,ib) = ssaaer(iplon,nlay-l,ib)
               asma(iplon,l,ib) = asmaer(iplon,nlay-l,ib)
            end do
         end do
      end do
   end if

   ! Transfer cloud fraction and cloud optical properties to RRTM variables;
   ! modify to reverse layer indexing here if necessary.

   if (icld .ge. 1) then 
      ! Move incoming GCM cloud arrays to RRTMG cloud arrays.
      ! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
      ! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

      do l = 1, nlay-1

         do ig = 1, ngptsw
            do iplon = 1, ncol
               zcldfmc(iplon,l,ig) = cldfmcl(ig,iplon,nlay-l)
               ztaucmc(iplon,l,ig) = taucmcl(ig,iplon,nlay-l)
               zssacmc(iplon,l,ig) = ssacmcl(ig,iplon,nlay-l)
               zasmcmc(iplon,l,ig) = asmcmcl(ig,iplon,nlay-l)

               fsfcmc(iplon,ig,l) = fsfcmcl(ig,iplon,nlay-l)
               ciwpmc(iplon,ig,l) = ciwpmcl(ig,iplon,nlay-l)
               clwpmc(iplon,ig,l) = clwpmcl(ig,iplon,nlay-l)
            end do
         end do

         do iplon = 1, ncol
            reicmc(iplon,l) = reicmcl(iplon,nlay-l)
            if (iceflag .eq. 3) then
               dgesmc(iplon,l) = 1.5396_r8 * reicmcl(iplon,nlay-l)
            end if
            relqmc(iplon,l) = relqmcl(iplon,nlay-l)
         end do
      end do

      ! If an extra layer is being used in RRTMG, set all cloud properties to zero
      ! in the extra layer.
      do iplon = 1, ncol
         zcldfmc(iplon,nlay,:) = 0.0_r8
         ztaucmc(iplon,nlay,:) = 0.0_r8
         zssacmc(iplon,nlay,:) = 1.0_r8
         zasmcmc(iplon,nlay,:) = 0.0_r8
         fsfcmc(iplon,:,nlay) = 0.0_r8
         ciwpmc(iplon,:,nlay) = 0.0_r8
         clwpmc(iplon,:,nlay) = 0.0_r8
         reicmc(iplon,nlay) = 0.0_r8
         dgesmc(iplon,nlay) = 0.0_r8
         relqmc(iplon,nlay) = 0.0_r8
         taua(iplon,nlay,:) = 0.0_r8
         ssaa(iplon,nlay,:) = 1.0_r8
         asma(iplon,nlay,:) = 0.0_r8
      end do
   end if

end subroutine inatm_sw

end module rrtmg_sw_rad


