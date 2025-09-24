module aer_rad_props

!------------------------------------------------------------------------------------------------
! Converts aerosol masses to bulk optical properties for sw and lw radiation
! computations.
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver
use physconst,        only: rga
use physics_types,    only: physics_state

use physics_buffer,   only: physics_buffer_desc
use radconstants,     only: nswbands, nlwbands, idx_sw_diag
use phys_prop,        only: nrh, ot_length
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                            rad_cnst_get_aer_props
use wv_saturation,    only: qsat
use aerosol_optics_cam,only: aerosol_optics_cam_init, aerosol_optics_cam_sw, aerosol_optics_cam_lw
use cam_history,      only: fieldname_len, addfld, outfld, add_default, horiz_only
use cam_history_support, only : fillvalue

use cam_abortutils,   only: endrun
use aer_vis_diag_mod, only: aer_vis_diag_init, aer_vis_diag_out

implicit none
private
save

public :: &
   aer_rad_props_init,        &
   aer_rad_props_sw,          & ! return SW optical props of aerosols
   aer_rad_props_lw             ! return LW optical props of aerosols

!==============================================================================
contains
!==============================================================================

subroutine aer_rad_props_init()
   use phys_control, only: phys_getopts

   integer                    :: numaerosols  ! number of aerosols
   logical                    :: history_amwg         ! output the variables used by the AMWG diag package
   logical                    :: history_aero_optics  ! Output aerosol optics diagnostics
   logical                    :: history_dust         ! Output dust diagnostics
   integer                    :: nmodes               ! number of aerosol modes
   integer                    :: nbins                ! number of aerosol bins

   !----------------------------------------------------------------------------

   call phys_getopts( history_aero_optics_out    = history_aero_optics, &
                      history_amwg_out           = history_amwg,    &
                      history_dust_out           = history_dust )

   call addfld ('AEROD_v ', horiz_only, 'A', '1', &
      'Total Aerosol Optical Depth in visible band', flag_xyfill=.true.)

   call addfld ('AODvstrt', horiz_only, 'A', '1', &
      'Stratospheric Aerosol Optical Depth in visible band', flag_xyfill=.true.)

   ! Contributions to AEROD_v from individual aerosols (climate species).

   ! number of bulk aerosols in climate list
   call rad_cnst_get_info(0, naero=numaerosols, nmodes=nmodes, nbins=nbins)

   call aer_vis_diag_init()

   ! Determine default fields
   if (history_amwg .or. history_dust .or. history_aero_optics) then
      call add_default ('AEROD_v', 1, ' ')
   endif

   if (numaerosols>0 .or. nmodes>0 .or. nbins>0) then
      call aerosol_optics_cam_init()
   end if

end subroutine aer_rad_props_init

!==============================================================================

subroutine aer_rad_props_sw(list_idx, state, pbuf,  nnite, idxnite, &
                            tau, tau_w, tau_w_g, tau_w_f)

   ! Return bulk layer tau, omega, g, f for all spectral intervals.

   use physics_buffer, only : physics_buffer_desc
   use tropopause,     only : tropopause_find_cam
   ! Arguments
   integer,             intent(in) :: list_idx      ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite                ! number of night columns
   integer,             intent(in) :: idxnite(:)           ! local column indices of night columns

   real(r8), intent(out) :: tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8), intent(out) :: tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8), intent(out) :: tau_w_g(pcols,0:pver,nswbands) ! aerosol asymmetry parameter * tau * w
   real(r8), intent(out) :: tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * tau * w

   ! Local variables

   integer :: ncol
   integer :: lchnk
   integer :: troplev(pcols)

   integer  :: numaerosols     ! number of bulk aerosols in climate/diagnostic list
   integer  :: nmodes          ! number of aerosol modes in climate/diagnostic list
   integer  :: nbins           ! number of aerosol bins in climate/diagnostic list

   !-----------------------------------------------------------------------------

   ncol  = state%ncol
   lchnk = state%lchnk

   ! top layer (ilev = 0) has no aerosol (ie tau = 0)
   ! also initialize rest of layers to accumulate optical depths
   tau    (1:ncol,:,:) = 0._r8
   tau_w  (1:ncol,:,:) = 0._r8
   tau_w_g(1:ncol,:,:) = 0._r8
   tau_w_f(1:ncol,:,:) = 0._r8

   ! get number of bulk aerosols and number of modes in current list
   call rad_cnst_get_info(list_idx, naero=numaerosols, nmodes=nmodes, nbins=nbins)

   ! Contributions from modal and bin aerosols.
   if (numaerosols>0 .or. nmodes>0 .or. nbins>0) then
      call aerosol_optics_cam_sw(list_idx, state, pbuf, nnite, idxnite, &
                                 tau, tau_w, tau_w_g, tau_w_f)
   end if

   !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   troplev = 0
   !REMOVECAM_END
   call tropopause_find_cam(state, troplev)

   ! diagnostic output of total aerosol optical properties
   ! currently implemented for climate list only
   call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, 0, tau(:,:,idx_sw_diag), list_idx, troplev)

end subroutine aer_rad_props_sw

!==============================================================================

subroutine aer_rad_props_lw(list_idx, state, pbuf, odap_aer)

   ! Purpose: Compute aerosol transmissions needed in absorptivity/
  !    emissivity calculations

   use physics_buffer, only : physics_buffer_desc

   ! Arguments
   integer,             intent(in)  :: list_idx                      ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state

   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8),            intent(out) :: odap_aer(pcols,pver,nlwbands) ! [fraction] absorption optical depth, per layer

   ! Local variables

   integer :: numaerosols ! number of bulk aerosols in climate/diagnostic list
   integer :: nmodes      ! number of aerosol modes in climate/diagnostic list
   integer :: nbins       ! number of aerosol bins in climate/diagnostic list

   !-----------------------------------------------------------------------------

   ! get number of bulk aerosols and number of modes in current list
   call rad_cnst_get_info(list_idx, naero=numaerosols, nmodes=nmodes, nbins=nbins)

   odap_aer = 0._r8

   ! Contributions from modal and sectional aerosols.
   if (numaerosols>0 .or. nmodes>0 .or. nbins>0) then
      call aerosol_optics_cam_lw(list_idx, state, pbuf, odap_aer)
   end if

end subroutine aer_rad_props_lw

!==============================================================================

end module aer_rad_props
