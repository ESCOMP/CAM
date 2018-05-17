!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/mcica_subcol_gen_lw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.3 $
!     created:   $Date: 2007/08/28 22:38:11 $
!

module mcica_subcol_gen_lw

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! Purpose: Create McICA stochastic arrays for cloud physical or optical properties.
! Two options are possible:
! 1) Input cloud physical properties: cloud fraction, ice and liquid water
!    paths, ice fraction, and particle sizes.  Output will be stochastic
!    arrays of these variables.  (inflag = 1)
! 2) Input cloud optical properties directly: cloud optical depth, single
!    scattering albedo and asymmetry parameter.  Output will be stochastic
!    arrays of these variables.  (inflag = 0; longwave scattering is not
!    yet available, ssac and asmc are for future expansion)

! --------- Modules ----------

use shr_kind_mod,     only: r8 => shr_kind_r8
use cam_abortutils,   only: endrun

use parrrtm, only : nbndlw, ngptlw
use rrlw_wvn, only: ngb

implicit none
private

public :: mcica_subcol_lw

!=========================================================================================
contains
!=========================================================================================

subroutine mcica_subcol_lw(lchnk, ncol, nlay, icld, permuteseed, play, &
                       cldfrac, ciwp, clwp, rei, rel, tauc, cldfmcl, &
                       ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl)

   ! ----- Input -----
   ! Control
   integer, intent(in) :: lchnk           ! chunk identifier
   integer, intent(in) :: ncol            ! number of columns
   integer, intent(in) :: nlay            ! number of model layers
   integer, intent(in) :: icld            ! clear/cloud, cloud overlap flag
   integer, intent(in) :: permuteseed     ! if the cloud generator is called multiple times, 
                                                        ! permute the seed between each call.
                                                        ! between calls for LW and SW, recommended
                                                        ! permuteseed differes by 'ngpt'

   ! Atmosphere
   real(kind=r8), intent(in) :: play(:,:)          ! layer pressures (mb) 
                                                        !    Dimensions: (ncol,nlay)

   ! Atmosphere/clouds - cldprop
   real(kind=r8), intent(in) :: cldfrac(:,:)       ! layer cloud fraction
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: tauc(:,:,:)        ! cloud optical depth
                                                     !    Dimensions: (nbndlw,ncol,nlay)
   real(kind=r8), intent(in) :: ciwp(:,:)          ! cloud ice water path
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: clwp(:,:)          ! cloud liquid water path
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: rei(:,:)           ! cloud ice particle size
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: rel(:,:)           ! cloud liquid particle size
                                                     !    Dimensions: (ncol,nlay)

   ! ----- Output -----
   ! Atmosphere/clouds - cldprmc [mcica]
   real(kind=r8), intent(out) :: cldfmcl(:,:,:)    ! cloud fraction [mcica]
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   real(kind=r8), intent(out) :: ciwpmcl(:,:,:)    ! cloud ice water path [mcica]
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   real(kind=r8), intent(out) :: clwpmcl(:,:,:)    ! cloud liquid water path [mcica]
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   real(kind=r8), intent(out) :: relqmcl(:,:)      ! liquid particle size (microns)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: reicmcl(:,:)      ! ice partcle size (microns)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(out) :: taucmcl(:,:,:)    ! cloud optical depth [mcica]
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   ! ----- Local -----

   ! Stochastic cloud generator variables [mcica]
   integer, parameter :: nsubclw = ngptlw ! number of sub-columns (g-point intervals)
   integer :: km, im, nm                  ! loop indices

   real(kind=r8) :: pmid(ncol, nlay)               ! layer pressures (Pa) 
   !----------------------------------------------------------------------------

   ! Return if clear sky; or stop if icld out of range
   if (icld.eq.0) return
   if (icld.lt.0.or.icld.gt.3) then 
      call endrun('MCICA_SUBCOL: INVALID ICLD')
   end if

   ! NOTE: For GCM mode, permuteseed must be offset between LW and SW by at
   !       least the number of subcolumns

   ! Pass particle sizes to new arrays, no subcolumns for these properties yet
   ! Convert pressures from mb to Pa

   reicmcl(:ncol,:nlay) = rei(:ncol,:nlay)
   relqmcl(:ncol,:nlay) = rel(:ncol,:nlay)
   pmid(:ncol,:nlay)    = play(:ncol,:nlay)*1.e2_r8

   ! Generate the stochastic subcolumns of cloud optical properties for the longwave;
   call generate_stochastic_clouds( &
      ncol, nlay, nsubclw, icld, pmid, &
      cldfrac, clwp, ciwp, tauc, cldfmcl, &
      clwpmcl, ciwpmcl, taucmcl, permuteseed)

end subroutine mcica_subcol_lw

!=========================================================================================

subroutine generate_stochastic_clouds( &
   ncol, nlay, nsubcol, icld, pmid, &
   cld, clwp, ciwp, tauc, cld_stoch, &
   clwp_stoch, ciwp_stoch, tauc_stoch, changeSeed) 

   !----------------------------------------------------------------------------------------------------------------
   ! ---------------------
   ! Contact: Cecile Hannay (hannay@ucar.edu)
   ! 
   ! Original code: Based on Raisanen et al., QJRMS, 2004.
   ! 
   ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
   !   random number generator, which can be changed to the optional kissvec random number generator
   !   with flag 'irnd' below. Some extra functionality has been commented or removed.  
   !   Michael J. Iacono, AER, Inc., February 2007
   !
   ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
   ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
   ! and uniform cloud liquid and cloud ice concentration.
   ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
   ! and obeys an overlap assumption in the vertical.   
   ! 
   ! Overlap assumption:
   !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
   !  The default option is maximum-random (option 3)
   !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
   !  This is set with the variable "overlap" 
   !mji - Exponential overlap option (overlap=4) has been deactivated in this version
   !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
   ! 
   ! Seed:
   !  If the stochastic cloud generator is called several times during the same timestep, 
   !  one should change the seed between the call to insure that the subcolumns are different.
   !  This is done by changing the argument 'changeSeed'
   !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
   !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
   !
   ! PDF assumption:
   !  We can use arbitrary complicated PDFS. 
   !  In the present version, we produce homogeneuous clouds (the simplest case).  
   !  Future developments include using the PDF scheme of Ben Johnson. 
   !
   ! History file:
   !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
   !  nsubcol = number of subcolumns
   !  overlap = overlap type (1-3)
   !  Zo = length scale 
   !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
   !  CLDLIQ_S = mean of the subcolumn cloud water
   !  CLDICE_S = mean of the subcolumn cloud ice 
   !
   ! Note:
   !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
   !   i.e we only have cloud condensate when the cell is cloudy. 
   !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
   !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
   !   without cloud condensate or the opposite).
   !---------------------------------------------------------------------------------------------------------------

   use shr_RandNum_mod, only: ShrIntrinsicRandGen, ShrKissRandGen, &
                              ShrF95MtRandGen, ShrDsfmtRandGen

   type(ShrDsfmtRandGen) :: dsfmt_gen
   type(ShrKissRandGen) :: kiss_gen

   ! -- Arguments

   integer, intent(in) :: ncol            ! number of columns
   integer, intent(in) :: nlay            ! number of layers
   integer, intent(in) :: icld            ! clear/cloud, cloud overlap flag
   integer, intent(in) :: nsubcol         ! number of sub-columns (g-point intervals)
   integer, optional, intent(in) :: changeSeed     ! allows permuting seed

   real(kind=r8), intent(in) :: pmid(:,:)          ! layer pressure (Pa)
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: cld(:,:)           ! cloud fraction 
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: clwp(:,:)          ! cloud liquid water path
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: ciwp(:,:)          ! cloud ice water path
                                                     !    Dimensions: (ncol,nlay)
   real(kind=r8), intent(in) :: tauc(:,:,:)        ! cloud optical depth
                                                     !    Dimensions: (nbndlw,ncol,nlay)

   real(kind=r8), intent(out) :: cld_stoch(:,:,:)  ! subcolumn cloud fraction 
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   real(kind=r8), intent(out) :: clwp_stoch(:,:,:) ! subcolumn cloud liquid water path
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   real(kind=r8), intent(out) :: ciwp_stoch(:,:,:) ! subcolumn cloud ice water path
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   real(kind=r8), intent(out) :: tauc_stoch(:,:,:) ! subcolumn cloud optical depth
                                                     !    Dimensions: (ngptlw,ncol,nlay)
   ! -- Local variables
   real(kind=r8) :: cldf(ncol,nlay)                ! cloud fraction 
    
   ! Constants (min value for cloud fraction and cloud water and ice)
   real(kind=r8), parameter :: cldmin = 1.0e-80_r8     ! min cloud fraction

   ! Variables related to random number and seed 
   integer :: irnd                        ! flag for random number generator
                                                        !  0 = kissvec
                                                        !  1 = Mersenne Twister

   real(kind=r8), dimension(nsubcol, ncol, nlay) :: CDF, CDF2      ! random numbers
   integer, dimension(ncol,4) :: kiss_seed
   real(kind=r8), dimension(ncol,1) :: rand_num_1d ! random number (kissvec)
   real(kind=r8), dimension(ncol,nlay) :: rand_num ! random number (kissvec)
   integer, dimension(ncol) :: iseed                       ! seed to create random number (Mersenne Teister)

   ! Flag to identify cloud fraction in subcolumns
   logical,  dimension(nsubcol, ncol, nlay) :: iscloudy   ! flag that says whether a gridbox is cloudy

   ! Indices
   integer :: ilev, isubcol, i, n         ! indices
   !----------------------------------------------------------------------------

   ! Set randum number generator to use (0 = kissvec; 1 = mersennetwister)
   irnd = 0

   ! ensure that cloud fractions are in bounds 
   cldf(:,:) = cld(:ncol,:nlay)
   where (cldf(:,:) < cldmin)
      cldf(:,:) = 0._r8
   end where

   ! ----- Create seed  --------
   
   ! Advance randum number generator by changeseed values
   if (irnd == 0) then   

      ! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
      ! Must use pmid from bottom four layers. 
      do i=1,ncol
         if (pmid(i,nlay).lt.pmid(i,nlay-1)) then 
            call endrun('MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.')
         end if
         kiss_seed(i,1) = (pmid(i,nlay)   - int(pmid(i,nlay)))    * 1000000000
         kiss_seed(i,2) = (pmid(i,nlay-1) - int(pmid(i,nlay-1)))  * 1000000000
         kiss_seed(i,3) = (pmid(i,nlay-2) - int(pmid(i,nlay-2)))  * 1000000000
         kiss_seed(i,4) = (pmid(i,nlay-3) - int(pmid(i,nlay-3)))  * 1000000000
      end do

      kiss_gen = ShrKissRandGen(kiss_seed)

      do i = 1, changeSeed
         call kiss_gen%random(rand_num_1d)
      end do
   elseif (irnd.eq.1) then

      do i = 1, ncol
         if (pmid(i,nlay) < pmid(i,nlay-1)) then
            call endrun('MCICA_SUBCOL: MT SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.')
         end if
         kiss_seed(i,1) = (pmid(i,nlay)   - int(pmid(i,nlay)))    * 1000000000
         kiss_seed(i,2) = (pmid(i,nlay-1) - int(pmid(i,nlay-1)))  * 1000000000
         kiss_seed(i,3) = (pmid(i,nlay-2) - int(pmid(i,nlay-2)))  * 1000000000
         kiss_seed(i,4) = (pmid(i,nlay-3) - int(pmid(i,nlay-3)))  * 1000000000
      end do

      iseed = kiss_seed(:,1) + kiss_seed(:,2) + kiss_seed(:,3) + kiss_seed(:,4)
      dsfmt_gen =ShrDsfmtRandGen(iseed,1)

   end if

   ! ------ Apply overlap assumption --------

   ! generate the random numbers  

   select case (icld)

   case(1) 
      ! Random overlap
      ! i) pick a random value at every level
  
      if (irnd == 0) then 
         do isubcol = 1,nsubcol
            call kiss_gen%random(rand_num)
            CDF(isubcol,:,:) = rand_num(:,:)
         end do
      else if (irnd == 1) then
         do isubcol = 1, nsubcol
            call dsfmt_gen%random(rand_num)
            CDF(isubcol,:,:) = rand_num(:,:) 
         end do
      end if

   case(2) 
      ! Maximum-Random overlap
      ! i) pick a random number for top layer.
      ! ii) walk down the column: 
      !    - if the layer above is cloudy, we use the same random number than in the layer above
      !    - if the layer above is clear, we use a new random number 

      if (irnd == 0) then 
         do isubcol = 1, nsubcol
            call kiss_gen%random(rand_num)
            CDF(isubcol,:,:) = rand_num(:,:)
         end do
      elseif (irnd == 1) then
         do isubcol = 1, nsubcol
            call dsfmt_gen%random(rand_num)
            CDF(isubcol,:,:) = rand_num(:,:)
         end do
      end if

      do ilev = 2,nlay
         do i = 1, ncol
            do isubcol = 1, nsubcol
               if (CDF(isubcol, i, ilev-1) > 1._r8 - cldf(i,ilev-1) ) then
                  CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev-1) 
               else
                  CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev) * (1._r8 - cldf(i,ilev-1))
               end if
            end do
         end do
      end do
       
   case(3) 
      ! Maximum overlap
      ! i) pick the same random numebr at every level  

      if (irnd.eq.0) then 
         do isubcol = 1,nsubcol
            call kiss_gen%random(rand_num_1d)
            do ilev = 1,nlay
               CDF(isubcol,:,ilev) = rand_num_1d(:,1)
            enddo
         enddo
      elseif (irnd.eq.1) then
         do isubcol = 1, nsubcol
            call dsfmt_gen%random(rand_num_1d)
            do ilev = 1, nlay
               CDF(isubcol,:,ilev) = rand_num_1d(:,1)
            enddo
         enddo
      endif

   end select
 
   ! -- generate subcolumns for homogeneous clouds -----
   do ilev = 1,nlay
      iscloudy(:,:,ilev) = (CDF(:,:,ilev) >= 1._r8 - spread(cldf(:,ilev), dim=1, nCopies=nsubcol) )
   end do

   ! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
   ! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0

   do ilev = 1,nlay
      do i = 1, ncol
         do isubcol = 1, nsubcol
            if (iscloudy(isubcol,i,ilev) ) then
               cld_stoch(isubcol,i,ilev) = 1._r8
            else
               cld_stoch(isubcol,i,ilev) = 0._r8
            endif
         end do
      end do
   end do

   ! where there is a cloud, set the subcolumn cloud properties;
   ! incoming clwp, ciwp and tauc should be in-cloud quantites and not grid-averaged quantities

   do ilev = 1, nlay
      do i = 1, ncol
         do isubcol = 1, nsubcol
            if ( iscloudy(isubcol,i,ilev) .and. (cldf(i,ilev) > 0._r8) ) then
               clwp_stoch(isubcol,i,ilev) = clwp(i,ilev)
               ciwp_stoch(isubcol,i,ilev) = ciwp(i,ilev)
            else
               clwp_stoch(isubcol,i,ilev) = 0._r8
               ciwp_stoch(isubcol,i,ilev) = 0._r8
            end if
         end do
      end do
   end do

   do ilev = 1, nlay
      do i = 1, ncol
         do isubcol = 1,nsubcol
            if ( iscloudy(isubcol,i,ilev) .and. (cldf(i,ilev) > 0._r8) ) then
               n = ngb(isubcol)
               tauc_stoch(isubcol,i,ilev) = tauc(n,i,ilev)
            else
               tauc_stoch(isubcol,i,ilev) = 0._r8
            end if
         end do
      end do
   end do

   if (irnd == 0) then 
      call kiss_gen%finalize()
   else if (irnd == 1) then
      call dsfmt_gen%finalize()
   end if

end subroutine generate_stochastic_clouds

end module mcica_subcol_gen_lw
