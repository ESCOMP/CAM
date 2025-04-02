module rrtmgp_lw_mcica_subcol_gen
! PEVERWHEE - dependencies = shr_RandNum_mod

!----------------------------------------------------------------------------------------
! 
! Purpose: Create McICA stochastic arrays for lw cloud optical properties.
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
! rrtmgp_lw_mcica_subcol_gen.F90 and mcica_subcol_gen_sw.F90
! 
! Uses the KISS random number generator.
!
! Overlap assumption: maximum-random.
! 
!----------------------------------------------------------------------------------------

use machine,                only: kind_phys
use shr_RandNum_mod,        only: ShrKissRandGen
use ccpp_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp_ccpp
use ccpp_optical_props,     only: ty_optical_props_1scl_ccpp

implicit none
private
save

public :: rrtmgp_lw_mcica_subcol_gen_run

!========================================================================================
contains
!========================================================================================

!>
!> \section arg_table_rrtmgp_lw_mcica_subcol_gen_run Argument Table
!! \htmlinclude rrtmgp_lw_mcica_subcol_gen_run.html
subroutine rrtmgp_lw_mcica_subcol_gen_run( &
   dolw, ktoprad, kdist, nbnd, ngpt, ncol, pver, nver, &
   changeseed, pmid, cldfrac, tauc, cloud_lw,     &
   errmsg, errflg )

   ! Arrays use CAM vertical index convention: index increases from top to bottom.
   ! This index ordering is assumed in the maximum-random overlap algorithm which starts
   ! at the top of a column and marches down, with each layer depending on the state
   ! of the layer above it.
   !
   ! For GCM mode, changeseed must be offset between LW and SW by at least the
   ! number of subcolumns

   ! arguments
   class(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist        ! Gas optics object
   logical,                          intent(in) :: dolw         ! Flag for whether to perform longwave calculation
   integer,                          intent(in) :: ktoprad      ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
   integer,                          intent(in) :: nbnd         ! Number of spectral bands
   integer,                          intent(in) :: ngpt         ! Number of subcolumns (g-point intervals)
   integer,                          intent(in) :: ncol         ! Number of columns
   integer,                          intent(in) :: pver         ! Number of model layers
   integer,                          intent(in) :: nver         ! Number of layers in radiation calculation
   integer,                          intent(in) :: changeseed   ! If the subcolumn generator is called multiple times, 
                                                                ! permute the seed between each call.
   real(kind_phys), dimension(:,:),   intent(in)  :: pmid       ! Layer pressures at midpoints (Pa)
   real(kind_phys), dimension(:,:),   intent(in)  :: cldfrac    ! Layer cloud fraction
   real(kind_phys), dimension(:,:,:), intent(in)  :: tauc       ! Cloud optical depth
   type(ty_optical_props_1scl_ccpp),  intent(inout) :: cloud_lw ! Cloud optics object
   character(len=*),                  intent(out)   :: errmsg
   integer,                           intent(out)   :: errflg

   ! Local variables

   integer :: idx, isubcol, kdx, ndx

   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys  ! min cloud fraction
   real(kind_phys) :: cldf(ncol,nver)      ! cloud fraction clipped to cldmin

   type(ShrKissRandGen) :: kiss_gen  ! KISS RNG object
   integer  :: kiss_seed(ncol,4)
   real(kind_phys) :: rand_num_1d(ncol,1)   ! random number (kissvec)
   real(kind_phys) :: rand_num(ncol,nver)   ! random number (kissvec)

   real(kind_phys) :: cdf(ngpt,ncol,nver)   ! random numbers
   logical  :: iscloudy(ngpt,ncol,nver)   ! flag that says whether a gridbox is cloudy
   real(kind_phys) :: taucmcl(ngpt,ncol,nver)
   !------------------------------------------------------------------------------------------ 

   ! Set error variables
   errflg = 0
   errmsg = ''

   ! If we're not doing longwave this timestep, no need to proceed
   if (.not. dolw) then
      return
   end if

   ! clip cloud fraction
   cldf(:,:) = cldfrac(:ncol,:)
   where (cldf(:,:) < cldmin)
      cldf(:,:) = 0._kind_phys
   end where

   ! Create a seed that depends on the state of the columns.
   ! Use pmid from bottom four layers. 
   do idx = 1, ncol
      kiss_seed(idx,1) = (pmid(idx,pver)   - int(pmid(idx,pver)))    * 1000000000
      kiss_seed(idx,2) = (pmid(idx,pver-1) - int(pmid(idx,pver-1)))  * 1000000000
      kiss_seed(idx,3) = (pmid(idx,pver-2) - int(pmid(idx,pver-2)))  * 1000000000
      kiss_seed(idx,4) = (pmid(idx,pver-3) - int(pmid(idx,pver-3)))  * 1000000000
   end do

   ! create the RNG object
   kiss_gen = ShrKissRandGen(kiss_seed)

   ! Advance randum number generator by changeseed values
   do idx = 1, changeSeed
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

   do kdx = 2, nver
      do idx = 1, ncol
         do isubcol = 1, ngpt
            if (cdf(isubcol,idx,kdx-1) > 1._kind_phys - cldf(idx,kdx-1) ) then
               cdf(isubcol,idx,kdx) = cdf(isubcol,idx,kdx-1) 
            else
               cdf(isubcol,idx,kdx) = cdf(isubcol,idx,kdx) * (1._kind_phys - cldf(idx,kdx-1))
            end if
         end do
      end do
   end do
 
   do kdx = 1, nver
      iscloudy(:,:,kdx) = (cdf(:,:,kdx) >= 1._kind_phys - spread(cldf(:,kdx), dim=1, nCopies=ngpt) )
   end do

   ! -- generate subcolumns for homogeneous clouds -----
   ! where there is a cloud, set the subcolumn cloud properties;
   ! incoming tauc should be in-cloud quantites and not grid-averaged quantities
   do kdx = 1,nver
      do idx = 1,ncol
         do isubcol = 1,ngpt
            if (iscloudy(isubcol,idx,kdx) .and. (cldf(idx,kdx) > 0._kind_phys) ) then
               ndx = kdist%gas_props%convert_gpt2band(isubcol)
               taucmcl(isubcol,idx,kdx) = tauc(ndx,idx,kdx)
            else
               taucmcl(isubcol,idx,kdx) = 0._kind_phys
            end if
         end do
      end do
   end do

   call kiss_gen%finalize()

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there
   cloud_lw%optical_props%tau = 0.0_kind_phys

   ! Set the properties on g-points
   do idx = 1, ngpt
      cloud_lw%optical_props%tau(:,ktoprad:,idx) = taucmcl(idx,:,:)
   end do

   ! validate checks that: tau > 0
   errmsg = cloud_lw%optical_props%validate()
   if (len_trim(errmsg) > 0) then
       errflg = 1
       return
   end if 

end subroutine rrtmgp_lw_mcica_subcol_gen_run


end module rrtmgp_lw_mcica_subcol_gen

