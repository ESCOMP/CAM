module oldcloud_optics

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use constituents,     only: cnst_get_ind
use physconst,        only: gravit
use radconstants,     only: nlwbands
use ebert_curry_ice_optics, only: scalefactor

use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   oldcloud_init,            &
   oldcloud_lw,              &
   old_liq_get_rad_props_lw, &
   old_ice_get_rad_props_lw

integer :: nmu, nlambda
real(r8), allocatable :: g_mu(:)           ! mu samples on grid
real(r8), allocatable :: g_lambda(:,:)     ! lambda scale samples on grid
real(r8), allocatable :: ext_sw_liq(:,:,:)
real(r8), allocatable :: ssa_sw_liq(:,:,:)
real(r8), allocatable :: asm_sw_liq(:,:,:)
real(r8), allocatable :: abs_lw_liq(:,:,:)

integer :: n_g_d
real(r8), allocatable :: g_d_eff(:)        ! radiative effective diameter samples on grid
real(r8), allocatable :: ext_sw_ice(:,:)
real(r8), allocatable :: ssa_sw_ice(:,:)
real(r8), allocatable :: asm_sw_ice(:,:)
real(r8), allocatable :: abs_lw_ice(:,:)

! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
real(r8), parameter :: cldmin = 1.0e-80_r8

! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
real(r8), parameter :: cldeps = 0.0_r8

! indexes into pbuf
integer :: iciwp_idx   = 0 
integer :: iclwp_idx   = 0 
integer :: cld_idx     = 0 
integer :: rel_idx     = 0 
integer :: rei_idx     = 0 

! indexes into constituents for old optics
integer :: &
   ixcldice, & ! cloud ice water index
   ixcldliq    ! cloud liquid water index


!==============================================================================
contains
!==============================================================================

subroutine oldcloud_init()


   integer :: err

   iciwp_idx  = pbuf_get_index('ICIWP',errcode=err)
   iclwp_idx  = pbuf_get_index('ICLWP',errcode=err)
   cld_idx    = pbuf_get_index('CLD')
   rel_idx    = pbuf_get_index('REL')
   rei_idx    = pbuf_get_index('REI')

   ! old optics
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

end subroutine oldcloud_init

!==============================================================================

subroutine oldcloud_lw(state,pbuf,cld_abs_od,oldwp)

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc),pointer :: pbuf(:)
   real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
   logical,intent(in)              :: oldwp ! use old definition of waterpath


   real(r8) :: gicewp(pcols,pver)
   real(r8) :: gliqwp(pcols,pver)
   real(r8) :: cicewp(pcols,pver)
   real(r8) :: cliqwp(pcols,pver)
   real(r8) :: ficemr(pcols,pver)
   real(r8) :: cwp(pcols,pver)
   real(r8) :: cldtau(pcols,pver)

   real(r8), pointer, dimension(:,:) :: cldn
   real(r8), pointer, dimension(:,:) :: rei 
   integer :: ncol, itim_old, lwband, i, k, lchnk
   real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth

   real(r8) :: kabs, kabsi
   real(r8), parameter :: kabsl = 0.090361_r8 ! longwave liquid absorption coeff (m**2/g)

 
   ncol = state%ncol
   lchnk = state%lchnk

   itim_old  =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, rei_idx,   rei)
   call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   if (oldwp) then
     do k=1,pver
         do i = 1,ncol
            gicewp(i,k) = state%q(i,k,ixcldice)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box ice water path.
            gliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box liquid water path.
            cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
            cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
            ficemr(i,k) = state%q(i,k,ixcldice) /                 &
                 max(1.e-10_r8,(state%q(i,k,ixcldice)+state%q(i,k,ixcldliq)))
         end do
     end do
     cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
   else
      if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
         call endrun('oldcloud_lw: oldwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
      endif
      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
      do k=1,pver
         do i = 1,ncol
            cwp(i,k) = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
            ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
         end do
      end do
   endif

   do k=1,pver
       do i=1,ncol

          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs = kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          cldtau(i,k) = kabs*cwp(i,k)
       end do
   end do
!
   do lwband = 1,nlwbands
      cld_abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
   enddo

end subroutine oldcloud_lw

!==============================================================================

subroutine old_liq_get_rad_props_lw(state, pbuf, abs_od, oldliqwp)

   type(physics_state), intent(in)  :: state
   type(physics_buffer_desc),pointer  :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)
   logical, intent(in) :: oldliqwp

   real(r8) :: gicewp(pcols,pver)
   real(r8) :: gliqwp(pcols,pver)
   real(r8) :: cicewp(pcols,pver)
   real(r8) :: cliqwp(pcols,pver)
   real(r8) :: ficemr(pcols,pver)
   real(r8) :: cwp(pcols,pver)
   real(r8) :: cldtau(pcols,pver)
   
   real(r8), pointer, dimension(:,:) :: cldn
   real(r8), pointer, dimension(:,:) :: rei
   integer :: ncol, itim_old, lwband, i, k, lchnk 

   real(r8) :: kabs, kabsi
   real(r8), parameter :: kabsl = 0.090361_r8 ! longwave liquid absorption coeff (m**2/g)

   real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth

   ncol=state%ncol
   lchnk = state%lchnk

   itim_old  =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, rei_idx,   rei)
   call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   if (oldliqwp) then
     do k=1,pver
         do i = 1,ncol
            gicewp(i,k) = state%q(i,k,ixcldice)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box ice water path.
            gliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box liquid water path.
            cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
            cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
            ficemr(i,k) = state%q(i,k,ixcldice) /                 &
                 max(1.e-10_r8,(state%q(i,k,ixcldice)+state%q(i,k,ixcldliq)))
         end do
     end do
     cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
   else
      if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
         call endrun('old_liq_get_rad_props_lw: oldliqwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
      endif
      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
      do k=1,pver
         do i = 1,ncol
            cwp(i,k) = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
            ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
         end do
      end do
   endif


   do k=1,pver
       do i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use liquid water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs = kabsl*(1._r8-ficemr(i,k)) ! + kabsi*ficemr(i,k)
          cldtau(i,k) = kabs*cwp(i,k)
       end do
   end do

   do lwband = 1,nlwbands
      abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
   enddo


end subroutine old_liq_get_rad_props_lw

!==============================================================================

subroutine old_ice_get_rad_props_lw(state, pbuf, abs_od, oldicewp)

   type(physics_state), intent(in)  :: state
   type(physics_buffer_desc),pointer  :: pbuf(:)
    real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)
   logical, intent(in) :: oldicewp

   real(r8) :: gicewp(pcols,pver)
   real(r8) :: gliqwp(pcols,pver)
   real(r8) :: cicewp(pcols,pver)
   real(r8) :: cliqwp(pcols,pver)
   real(r8) :: ficemr(pcols,pver)
   real(r8) :: cwp(pcols,pver)
   real(r8) :: cldtau(pcols,pver)

   real(r8), pointer, dimension(:,:) :: cldn
   real(r8), pointer, dimension(:,:) :: rei
   integer :: ncol, itim_old, lwband, i, k, lchnk

   real(r8) :: kabs, kabsi
   real(r8), parameter :: kabsl = 0.090361_r8 ! longwave liquid absorption coeff (m**2/g)

   real(r8), pointer, dimension(:,:) :: iclwpth, iciwpth


   ncol = state%ncol
   lchnk = state%lchnk

   itim_old  =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, rei_idx,   rei)
   call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   if(oldicewp) then
     do k=1,pver
         do i = 1,ncol
            gicewp(i,k) = state%q(i,k,ixcldice)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box ice water path.
            gliqwp(i,k) = state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box liquid water path.
            cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
            cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
            ficemr(i,k) = state%q(i,k,ixcldice) /                 &
                 max(1.e-10_r8,(state%q(i,k,ixcldice)+state%q(i,k,ixcldliq)))
         end do
     end do
     cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
   else
      if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
         call endrun('old_ice_get_rad_props_lw: oldicewp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
      endif
      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
      do k=1,pver
         do i = 1,ncol
            cwp(i,k) = 1000.0_r8 *iciwpth(i,k) + 1000.0_r8 *iclwpth(i,k)
            ficemr(i,k) = 1000.0_r8*iciwpth(i,k)/(max(1.e-18_r8,cwp(i,k)))
         end do
      end do
   endif

   do k=1,pver
       do i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use ice water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs =  kabsi*ficemr(i,k) ! kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          cldtau(i,k) = kabs*cwp(i,k)
       end do
   end do

   do lwband = 1,nlwbands
      abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
   enddo

end subroutine old_ice_get_rad_props_lw

!==============================================================================

end module oldcloud_optics
