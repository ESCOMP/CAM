module rrtmgp_inputs_cam

!--------------------------------------------------------------------------------
! Transform data for inputs from CAM's data structures to those used by
! RRTMGP.  Subset the number of model levels if CAM's top exceeds RRTMGP's
! valid domain.  Add an extra layer if CAM's top is below 1 Pa.
! The vertical indexing increases from top to bottom of atmosphere in both
! CAM and RRTMGP arrays.   
!--------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use ppgrid,           only: pcols, pver, pverp

use physconst,        only: stebol, pi

use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc
use camsrfexch,       only: cam_in_t

use radconstants,     only: nradgas, gaslist, nswbands, nlwbands

use rad_constituents, only: rad_cnst_get_gas

use cloud_rad_props,  only: get_liquid_optics_sw, liquid_cloud_get_rad_props_lw, &
                            get_ice_optics_sw,    ice_cloud_get_rad_props_lw,    &
                            get_snow_optics_sw,   snow_cloud_get_rad_props_lw,   &
                            get_grau_optics_sw,   grau_cloud_get_rad_props_lw
                                 
use mcica_subcol_gen, only: mcica_subcol_sw

use aer_rad_props,    only: aer_rad_props_sw, aer_rad_props_lw

use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
use ccpp_optical_props,      only: ty_optical_props_2str_ccpp, ty_optical_props_1scl_ccpp

use cam_history_support,   only: fillvalue
use cam_abortutils,   only: endrun
use error_messages,   only: alloc_err
use radiation_utils, only: get_sw_spectral_boundaries_ccpp

implicit none
private
save

public :: &
   rrtmgp_inputs_cam_init, &
   rrtmgp_get_gas_mmrs, &
   rrtmgp_set_aer_lw,   &
   rrtmgp_set_aer_sw


! This value is to match the arbitrary small value used in RRTMG to decide
! when a quantity is effectively zero.
real(r8), parameter :: tiny = 1.0e-80_r8
real(r8) :: sw_low_bounds(nswbands)
real(r8) :: sw_high_bounds(nswbands)
integer :: ktopcam
integer :: ktoprad
integer :: idx_sw_diag
integer :: idx_nir_diag
integer :: idx_uv_diag
integer :: idx_sw_cloudsim
integer :: idx_lw_diag
integer :: idx_lw_cloudsim

! Mapping from RRTMG shortwave bands to RRTMGP.  Currently needed to continue using
! the SW optics datasets from RRTMG (even thought there is a slight mismatch in the
! band boundaries of the 2 bands that overlap with the LW bands).
integer, parameter, dimension(14) :: rrtmg_to_rrtmgp_swbands = &
   [ 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]

!==================================================================================================
contains
!==================================================================================================

!==================================================================================================
subroutine rrtmgp_inputs_cam_init(ktcam, ktrad, idx_sw_diag_in, idx_nir_diag_in, idx_uv_diag_in, &
               idx_sw_cloudsim_in, idx_lw_diag_in, idx_lw_cloudsim_in)
      
   ! Note that this routine must be called after the calls to set_wavenumber_bands which set
   ! the sw/lw band boundaries in the radconstants module.

   integer, intent(in) :: ktcam
   integer, intent(in) :: ktrad
   integer, intent(in) :: idx_sw_diag_in
   integer, intent(in) :: idx_nir_diag_in
   integer, intent(in) :: idx_uv_diag_in
   integer, intent(in) :: idx_sw_cloudsim_in
   integer, intent(in) :: idx_lw_diag_in
   integer, intent(in) :: idx_lw_cloudsim_in
   character(len=512) :: errmsg
   integer :: errflg

   ktopcam = ktcam
   ktoprad = ktrad
   idx_sw_diag = idx_sw_diag_in
   idx_nir_diag = idx_nir_diag_in
   idx_uv_diag = idx_uv_diag_in
   idx_sw_cloudsim = idx_sw_cloudsim_in
   idx_lw_diag = idx_lw_diag_in
   idx_lw_cloudsim = idx_lw_cloudsim_in

   ! Initialize the module data containing the SW band boundaries.
   call get_sw_spectral_boundaries_ccpp(sw_low_bounds, sw_high_bounds, 'cm^-1', errmsg, errflg)
   if (errflg /= 0) then
      call endrun('rrtmgp_inputs_cam_init: error during get_sw_spectral_boundaries_ccpp - message: '//errmsg)
   end if

end subroutine rrtmgp_inputs_cam_init

!=========================================================================================

subroutine rrtmgp_get_gas_mmrs(icall, state, pbuf, nlay, gas_mmrs)

   ! Retrieve mass mixing ratios for radiatively active gases from rad_constituents

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   real(r8),                    intent(out)   :: gas_mmrs(:,:,:)

   ! local variables
   integer :: i, ncol
   real(r8), pointer :: gas_mmr(:,:)
   character(len=*), parameter :: sub = 'rrtmgp_get_gas_mmrs'
   !--------------------------------------------------------------------------------

   ncol = state%ncol
   do i = 1, nradgas
      call rad_cnst_get_gas(icall, gaslist(i), state, pbuf, gas_mmr)
      gas_mmrs(:,:,i) = gas_mmr(:ncol,:)
   end do
end subroutine rrtmgp_get_gas_mmrs

!==================================================================================================

subroutine rrtmgp_set_aer_lw(icall, state, pbuf, aer_lw)

   ! Load LW aerosol optical properties into the RRTMGP object.

   ! Arguments
   integer,                     intent(in) :: icall
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   type(ty_optical_props_1scl_ccpp), intent(inout) :: aer_lw

   ! Local variables
   integer :: ncol

   ! Aerosol LW absorption optical depth
   real(r8) :: aer_lw_abs (pcols,pver,nlwbands)

   character(len=*), parameter :: sub = 'rrtmgp_set_aer_lw'
   character(len=128) :: errmsg
   !--------------------------------------------------------------------------------

   ncol = state%ncol

   ! Get aerosol longwave optical properties.
   call aer_rad_props_lw(icall, state, pbuf, aer_lw_abs)

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there.
   aer_lw%optical_props%tau = 0.0_r8

   aer_lw%optical_props%tau(:ncol, ktoprad:, :) = aer_lw_abs(:ncol, ktopcam:, :)

   errmsg = aer_lw%optical_props%validate()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: aer_lw%optical_props%validate: '//trim(errmsg))
   end if
end subroutine rrtmgp_set_aer_lw

!==================================================================================================

subroutine rrtmgp_set_aer_sw( &
   icall, state, pbuf, nday, idxday, nnite, idxnite, aer_sw)

   ! Load SW aerosol optical properties into the RRTMGP object.

   ! Arguments
   integer,                     intent(in) :: icall
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   integer,  intent(in) :: nday
   integer,  intent(in) :: idxday(:)
   integer,  intent(in) :: nnite          ! number of night columns
   integer,  intent(in) :: idxnite(pcols) ! indices of night columns in the chunk

   type(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw

   ! local variables
   integer  :: i

   ! The optical arrays dimensioned in the vertical as 0:pver.
   ! The index 0 is for the extra layer used in the radiation
   ! calculation.  The index ktopcam assumes the CAM vertical indices are
   ! in the range 1:pver, so using this index correctly ignores vertical
   ! index 0.  If an "extra" layer is used in the calculations, it is
   ! provided and set in the RRTMGP aerosol object aer_sw.
   real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! extinction optical depth
   real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! single scattering albedo * tau
   real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! asymmetry parameter * w * tau
   real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! forward scattered fraction * w * tau
                                                  ! aer_tau_w_f is not used by RRTMGP.
   character(len=*), parameter :: sub = 'rrtmgp_set_aer_sw'
   !--------------------------------------------------------------------------------

   ! Get aerosol shortwave optical properties.
   ! Make outfld calls for aerosol optical property diagnostics.
   call aer_rad_props_sw( &
      icall, state, pbuf, nnite, idxnite, &
      aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)

   ! The aer_sw object is only initialized if nday > 0.
   if (nday > 0) then

      ! aerosol optical properties need to be re-ordered from the RRTMG spectral bands
      ! (as assumed in the optics datasets) to the RRTMGP band order.
      aer_tau(:,:,:)     = aer_tau(    :,:,rrtmg_to_rrtmgp_swbands)
      aer_tau_w(:,:,:)   = aer_tau_w(  :,:,rrtmg_to_rrtmgp_swbands)
      aer_tau_w_g(:,:,:) = aer_tau_w_g(:,:,rrtmg_to_rrtmgp_swbands)
                  
      ! If there is an extra layer in the radiation then this initialization
      ! will provide default values.
      aer_sw%optical_props%tau = 0.0_r8
      aer_sw%optical_props%ssa = 1.0_r8
      aer_sw%optical_props%g   = 0.0_r8

      ! CAM fields are products tau, tau*ssa, tau*ssa*asy
      ! Fields expected by RRTMGP are computed by
      ! aer_sw%optical_props%tau = aer_tau
      ! aer_sw%optical_props%ssa = aer_tau_w / aer_tau
      ! aer_sw%optical_props%g   = aer_tau_w_g / aer_taw_w
      ! aer_sw arrays have dimensions of (nday,nlay,nswbands)

      do i = 1, nday
         ! set aerosol optical depth, clip to zero
         aer_sw%optical_props%tau(i,ktoprad:,:) = max(aer_tau(idxday(i),ktopcam:,:), 0._r8)
         ! set value of single scattering albedo
         aer_sw%optical_props%ssa(i,ktoprad:,:) = merge(aer_tau_w(idxday(i),ktopcam:,:)/aer_tau(idxday(i),ktopcam:,:), &
                                          1._r8, aer_tau(idxday(i),ktopcam:,:) > 0._r8)
         ! set value of asymmetry
         aer_sw%optical_props%g(i,ktoprad:,:) = merge(aer_tau_w_g(idxday(i),ktopcam:,:)/aer_tau_w(idxday(i),ktopcam:,:), &
                                        0._r8, aer_tau_w(idxday(i),ktopcam:,:) > tiny)
      end do

      ! impose limits on the components
      aer_sw%optical_props%ssa = min(max(aer_sw%optical_props%ssa, 0._r8), 1._r8)
      aer_sw%optical_props%g   = min(max(aer_sw%optical_props%g, -1._r8), 1._r8)

   end if

end subroutine rrtmgp_set_aer_sw

!==================================================================================================

end module rrtmgp_inputs_cam
