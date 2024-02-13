module mcica_subcol_gen

!----------------------------------------------------------------------------------------
! 
! Purpose: Create McICA stochastic arrays for cloud optical properties.
! Input cloud optical properties directly: cloud optical depth, single
! scattering albedo and asymmetry parameter.  Output will be stochastic
! arrays of these variables.  (longwave scattering is not yet available)
!
! Original code: From RRTMG, with the following copyright notice,
! based on Raisanen et al., QJRMS, 2004:
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
! This code is a refactored version of code originally in the files
! mcica_subcol_gen_lw.F90 and mcica_subcol_gen_sw.F90
! 
! Uses the KISS random number generator.
!
! Overlap assumption: maximum-random.
! 
!----------------------------------------------------------------------------------------

use shr_kind_mod,         only: r8 => shr_kind_r8
use ppgrid,               only: pcols, pver
use shr_RandNum_mod,      only: ShrKissRandGen
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

implicit none
private
save

public :: mcica_subcol_lw, mcica_subcol_sw

!========================================================================================
contains
!========================================================================================

subroutine mcica_subcol_lw( &
   kdist, nbnd, ngpt, ncol, nver,           &
   changeseed, pmid, cldfrac, tauc, taucmcl )

   ! Arrays use CAM vertical index convention: index increases from top to bottom.
   ! This index ordering is assumed in the maximum-random overlap algorithm which starts
   ! at the top of a column and marches down, with each layer depending on the state
   ! of the layer above it.
   !
   ! For GCM mode, changeseed must be offset between LW and SW by at least the
   ! number of subcolumns

   ! arguments
   class(ty_gas_optics_rrtmgp), intent(in) :: kdist  ! spectral information
   integer,  intent(in)  :: nbnd                     ! number of spectral bands
   integer,  intent(in)  :: ngpt                     ! number of subcolumns (g-point intervals)
   integer,  intent(in)  :: ncol                     ! number of columns
   integer,  intent(in)  :: nver                     ! number of layers
   integer,  intent(in)  :: changeseed               ! if the subcolumn generator is called multiple times, 
                                                     ! permute the seed between each call.
   real(r8), intent(in)  :: pmid(pcols,pver)         ! layer pressures (Pa)
   real(r8), intent(in)  :: cldfrac(ncol,nver)       ! layer cloud fraction
   real(r8), intent(in)  :: tauc(nbnd,ncol,nver)     ! cloud optical depth
   real(r8), intent(out) :: taucmcl(ngpt,ncol,nver)  ! subcolumn cloud optical depth [mcica]

   ! Local variables

   integer :: i, isubcol, k, n

   real(r8), parameter :: cldmin = 1.0e-80_r8  ! min cloud fraction
   real(r8) :: cldf(ncol,pver)      ! cloud fraction clipped to cldmin

   type(ShrKissRandGen) :: kiss_gen  ! KISS RNG object
   integer  :: kiss_seed(ncol,4)
   real(r8) :: rand_num_1d(ncol,1)   ! random number (kissvec)
   real(r8) :: rand_num(ncol,nver)   ! random number (kissvec)

   real(r8) :: cdf(ngpt,ncol,nver)   ! random numbers
   logical  :: iscloudy(ngpt,ncol,nver)   ! flag that says whether a gridbox is cloudy
   !------------------------------------------------------------------------------------------ 

   ! clip cloud fraction
   cldf(:,:) = cldfrac(:ncol,:)
   where (cldf(:,:) < cldmin)
      cldf(:,:) = 0._r8
   end where

   ! Create a seed that depends on the state of the columns.
   ! Use pmid from bottom four layers. 
   do i = 1, ncol
      kiss_seed(i,1) = (pmid(i,pver)   - int(pmid(i,pver)))    * 1000000000
      kiss_seed(i,2) = (pmid(i,pver-1) - int(pmid(i,pver-1)))  * 1000000000
      kiss_seed(i,3) = (pmid(i,pver-2) - int(pmid(i,pver-2)))  * 1000000000
      kiss_seed(i,4) = (pmid(i,pver-3) - int(pmid(i,pver-3)))  * 1000000000
   end do

   ! create the RNG object
   kiss_gen = ShrKissRandGen(kiss_seed)

   ! Advance randum number generator by changeseed values
   do i = 1, changeSeed
      call kiss_gen%random(rand_num_1d)
   end do

   ! Generate random numbers in each subcolumn at every level
   do isubcol = 1,ngpt
      call kiss_gen%random(rand_num)
      cdf(isubcol,:,:) = rand_num(:,:)
   enddo

   ! Maximum-Random overlap
   ! i) pick a random number for top layer.
   ! ii) walk down the column: 
   !    - if the layer above is cloudy, use the same random number as in the layer above
   !    - if the layer above is clear, use a new random number 

   do k = 2, nver
      do i = 1, ncol
         do isubcol = 1, ngpt
            if (cdf(isubcol,i,k-1) > 1._r8 - cldf(i,k-1) ) then
               cdf(isubcol,i,k) = cdf(isubcol,i,k-1) 
            else
               cdf(isubcol,i,k) = cdf(isubcol,i,k) * (1._r8 - cldf(i,k-1))
            end if
         end do
      end do
   end do
 
   do k = 1, nver
      iscloudy(:,:,k) = (cdf(:,:,k) >= 1._r8 - spread(cldf(:,k), dim=1, nCopies=ngpt) )
   end do

   ! -- generate subcolumns for homogeneous clouds -----
   ! where there is a cloud, set the subcolumn cloud properties;
   ! incoming tauc should be in-cloud quantites and not grid-averaged quantities
   do k = 1,nver
      do i = 1,ncol
         do isubcol = 1,ngpt
            if (iscloudy(isubcol,i,k) .and. (cldf(i,k) > 0._r8) ) then
               n = kdist%convert_gpt2band(isubcol)
               taucmcl(isubcol,i,k) = tauc(n,i,k)
            else
               taucmcl(isubcol,i,k) = 0._r8
            end if
         end do
      end do
   end do

   call kiss_gen%finalize()

end subroutine mcica_subcol_lw

!========================================================================================

subroutine mcica_subcol_sw( &
   kdist, nbnd, ngpt, ncol, nlay, nver, changeseed, &
   pmid, cldfrac, tauc, ssac, asmc,     &
   taucmcl, ssacmcl, asmcmcl)

   ! Arrays use CAM vertical index convention: index increases from top to bottom.
   ! This index ordering is assumed in the maximum-random overlap algorithm which starts
   ! at the top of a column and marches down, with each layer depending on the state
   ! of the layer above it.
   !
   ! For GCM mode, changeseed must be offset between LW and SW by at least the
   ! number of subcolumns

   ! arguments
   class(ty_gas_optics_rrtmgp), intent(in) :: kdist  ! spectral information
   integer,  intent(in)  :: nbnd                    ! number of spectral bands
   integer,  intent(in)  :: ngpt                    ! number of subcolumns (g-point intervals)
   integer,  intent(in)  :: ncol                    ! number of columns
   integer,  intent(in)  :: nlay                    ! number of vertical layers in radiation calc;
                                                    !   may include an "extra layer"
   integer,  intent(in)  :: nver                    ! number of CAM's vertical layers in rad calc
   integer,  intent(in)  :: changeseed              ! if the subcolumn generator is called multiple times, 
                                                    ! permute the seed between each call.
   real(r8), intent(in)  :: pmid(ncol,nlay)         ! layer midpoint pressures (Pa)
   real(r8), intent(in)  :: cldfrac(ncol,nver)      ! layer cloud fraction
   real(r8), intent(in)  :: tauc(nbnd,ncol,nver)    ! cloud optical depth
   real(r8), intent(in)  :: ssac(nbnd,ncol,nver)    ! cloud single scattering albedo (non-delta scaled)
   real(r8), intent(in)  :: asmc(nbnd,ncol,nver)    ! cloud asymmetry parameter (non-delta scaled)


   real(r8), intent(out) :: taucmcl(ngpt,ncol,nver) ! subcolumn cloud optical depth [mcica]
   real(r8), intent(out) :: ssacmcl(ngpt,ncol,nver) ! subcolumn cloud single scattering albedo [mcica]
   real(r8), intent(out) :: asmcmcl(ngpt,ncol,nver) ! subcolumn cloud asymmetry parameter [mcica]

   ! Local vars

   integer :: i, isubcol, k, n

   real(r8), parameter :: cldmin = 1.0e-80_r8  ! min cloud fraction
   real(r8) :: cldf(ncol,nver)      ! cloud fraction clipped to cldmin

   type(ShrKissRandGen) :: kiss_gen  ! KISS RNG object
   integer  :: kiss_seed(ncol,4)
   real(r8) :: rand_num_1d(ncol,1)   ! random number (kissvec)
   real(r8) :: rand_num(ncol,nver)   ! random number (kissvec)

   real(r8) :: cdf(ngpt,ncol,nver)   ! random numbers
   logical  :: iscloudy(ngpt,ncol,nver)   ! flag that says whether a gridbox is cloudy
   !------------------------------------------------------------------------------------------ 

   ! clip cloud fraction
   cldf(:,:) = cldfrac(:ncol,:)
   where (cldf(:,:) < cldmin)
      cldf(:,:) = 0._r8
   end where

   ! Create a seed that depends on the state of the columns.
   ! Use pmid from bottom four layers. 
   do i = 1, ncol
      kiss_seed(i,1) = (pmid(i,nlay)   - int(pmid(i,nlay)))    * 1000000000
      kiss_seed(i,2) = (pmid(i,nlay-1) - int(pmid(i,nlay-1)))  * 1000000000
      kiss_seed(i,3) = (pmid(i,nlay-2) - int(pmid(i,nlay-2)))  * 1000000000
      kiss_seed(i,4) = (pmid(i,nlay-3) - int(pmid(i,nlay-3)))  * 1000000000
   end do

   ! create the RNG object
   kiss_gen = ShrKissRandGen(kiss_seed)

   ! Advance randum number generator by changeseed values
   do i = 1, changeSeed
      call kiss_gen%random(rand_num_1d)
   end do

   ! Generate random numbers in each subcolumn at every level
   do isubcol = 1,ngpt
      call kiss_gen%random(rand_num)
      cdf(isubcol,:,:) = rand_num(:,:)
   enddo

   ! Maximum-Random overlap
   ! i) pick a random number for top layer.
   ! ii) walk down the column: 
   !    - if the layer above is cloudy, use the same random number as in the layer above
   !    - if the layer above is clear, use a new random number 

   do k = 2, nver
      do i = 1, ncol
         do isubcol = 1, ngpt
            if (cdf(isubcol,i,k-1) > 1._r8 - cldf(i,k-1) ) then
               cdf(isubcol,i,k) = cdf(isubcol,i,k-1) 
            else
               cdf(isubcol,i,k) = cdf(isubcol,i,k) * (1._r8 - cldf(i,k-1))
            end if
         end do
      end do
   end do
 
   do k = 1, nver
      iscloudy(:,:,k) = (cdf(:,:,k) >= 1._r8 - spread(cldf(:,k), dim=1, nCopies=ngpt) )
   end do

   ! -- generate subcolumns for homogeneous clouds -----
   ! where there is a cloud, set the subcolumn cloud properties;
   ! incoming tauc should be in-cloud quantites and not grid-averaged quantities
   do k = 1,nver
      do i = 1,ncol
         do isubcol = 1,ngpt
            if (iscloudy(isubcol,i,k) .and. (cldf(i,k) > 0._r8) ) then
               n = kdist%convert_gpt2band(isubcol)
               taucmcl(isubcol,i,k) = tauc(n,i,k)
               ssacmcl(isubcol,i,k) = ssac(n,i,k)
               asmcmcl(isubcol,i,k) = asmc(n,i,k)
            else
               taucmcl(isubcol,i,k) = 0._r8
               ssacmcl(isubcol,i,k) = 1._r8
               asmcmcl(isubcol,i,k) = 0._r8
            end if
         end do
      end do
   end do

   call kiss_gen%finalize()

end subroutine mcica_subcol_sw


end module mcica_subcol_gen

