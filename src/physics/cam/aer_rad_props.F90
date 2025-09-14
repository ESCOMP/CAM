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
! Placed here due to PGI bug.
use ref_pres,         only: clim_modal_aero_top_lev

use cam_abortutils,   only: endrun
use bam_optics_diags_mod, only: bam_optics_diags_init, bam_optics_diags_out

implicit none
private
save

integer :: top_lev = 1

public :: &
   aer_rad_props_init,        &
   aer_rad_props_sw,          & ! return SW optical props of aerosols
   aer_rad_props_lw             ! return LW optical props of aerosols

!==============================================================================
contains
!==============================================================================

subroutine aer_rad_props_init()
   use phys_control, only: phys_getopts


   integer                    :: i
   integer                    :: numaerosols  ! number of aerosols
   character(len=64), pointer :: aernames(:)  ! aerosol names
   logical                    :: history_amwg         ! output the variables used by the AMWG diag package
   logical                    :: history_aero_optics  ! Output aerosol optics diagnostics
   logical                    :: history_dust         ! Output dust diagnostics
   logical                    :: prog_modal_aero      ! Prognostic modal aerosols present
   integer                    :: nmodes               ! number of aerosol modes
   integer                    :: nbins                ! number of aerosol bins

   !----------------------------------------------------------------------------

   call phys_getopts( history_aero_optics_out    = history_aero_optics, &
                      history_amwg_out           = history_amwg,    &
                      history_dust_out           = history_dust,    &
                      prog_modal_aero_out        = prog_modal_aero )

   ! Limit modal aerosols with top_lev here.
   if (prog_modal_aero) top_lev = clim_modal_aero_top_lev

   call addfld ('AEROD_v ', horiz_only, 'A', '1', &
      'Total Aerosol Optical Depth in visible band', flag_xyfill=.true.)

   call addfld ('AODvstrt', horiz_only, 'A', '1', &
      'Stratospheric Aerosol Optical Depth in visible band', flag_xyfill=.true.)

   ! Contributions to AEROD_v from individual aerosols (climate species).

   ! number of bulk aerosols in climate list
   call rad_cnst_get_info(0, naero=numaerosols)

   ! get names of bulk aerosols
   allocate(aernames(numaerosols))
   call rad_cnst_get_info(0, aernames=aernames, nmodes=nmodes, nbins=nbins)

   call bam_optics_diags_init()

   ! Determine default fields
   if (history_amwg .or. history_dust ) then
      call add_default ('AEROD_v', 1, ' ')
   endif

   if (numaerosols>0 .or. nmodes>0 .or. nbins>0) then
      call aerosol_optics_cam_init()
   end if

   deallocate(aernames)

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
   integer :: k       ! index
   integer :: troplev(pcols)

   ! optical props for each aerosol
   ! hygroscopic
   real(r8), pointer :: h_ext(:,:)
   real(r8), pointer :: h_ssa(:,:)
   real(r8), pointer :: h_asm(:,:)
   ! non-hygroscopic
   real(r8), pointer :: n_ext(:)
   real(r8), pointer :: n_ssa(:)
   real(r8), pointer :: n_asm(:)
   real(r8), pointer :: n_scat(:)
   real(r8), pointer :: n_ascat(:)
   ! radius-dependent
   real(r8), pointer :: r_ext(:,:)    ! radius-dependent mass-specific extinction
   real(r8), pointer :: r_scat(:,:)
   real(r8), pointer :: r_ascat(:,:)
   real(r8), pointer :: r_mu(:)       ! log(radius) domain variable for r_ext, r_scat, r_ascat

   ! radiative properties for each aerosol
   real(r8) :: ta (pcols,pver,nswbands)
   real(r8) :: tw (pcols,pver,nswbands)
   real(r8) :: twf(pcols,pver,nswbands)
   real(r8) :: twg(pcols,pver,nswbands)

   ! aerosol masses
   real(r8), pointer :: aermmr(:,:)    ! mass mixing ratio of aerosols
   real(r8) :: mmr_to_mass(pcols,pver) ! conversion factor for mmr to mass
   real(r8) :: aermass(pcols,pver)     ! mass of aerosols

   ! for table lookup into rh grid
   real(r8) :: es(pcols,pver)     ! saturation vapor pressure
   real(r8) :: qs(pcols,pver)     ! saturation specific humidity
   real(r8) :: rh(pcols,pver)
   real(r8) :: rhtrunc(pcols,pver)
   real(r8) :: wrh(pcols,pver)
   integer  :: krh(pcols,pver)

   integer  :: numaerosols     ! number of bulk aerosols in climate/diagnostic list
   integer  :: nmodes          ! number of aerosol modes in climate/diagnostic list
   integer  :: nbins           ! number of aerosol bins in climate/diagnostic list
   integer  :: iaerosol        ! index into bulk aerosol list

   character(len=ot_length) :: opticstype       ! hygro or nonhygro
   character(len=16) :: pbuf_fld
   !-----------------------------------------------------------------------------

   ncol  = state%ncol
   lchnk = state%lchnk

   ! compute mixing ratio to mass conversion
   do k = 1, pver
      mmr_to_mass(:ncol,k) = rga * state%pdeldry(:ncol,k)
   enddo

   ! initialize to conditions that would cause failure
   tau     (:,:,:) = -100._r8
   tau_w   (:,:,:) = -100._r8
   tau_w_g (:,:,:) = -100._r8
   tau_w_f (:,:,:) = -100._r8

   ! top layer (ilev = 0) has no aerosol (ie tau = 0)
   ! also initialize rest of layers to accumulate od's
   tau    (1:ncol,:,:) = 0._r8
   tau_w  (1:ncol,:,:) = 0._r8
   tau_w_g(1:ncol,:,:) = 0._r8
   tau_w_f(1:ncol,:,:) = 0._r8

   ! calculate relative humidity for table lookup into rh grid
   do k = 1, pver
      call qsat(state%t(1:ncol,k), state%pmid(1:ncol,k), es(1:ncol,k), qs(1:ncol,k), ncol)
   end do
   rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qs(1:ncol,1:pver)

   rhtrunc(1:ncol,1:pver) = min(rh(1:ncol,1:pver),1._r8)
   krh(1:ncol,1:pver) = min(floor( rhtrunc(1:ncol,1:pver) * nrh ) + 1, nrh - 1) ! index into rh mesh
   wrh(1:ncol,1:pver) = rhtrunc(1:ncol,1:pver) * nrh - krh(1:ncol,1:pver)       ! (-) weighting on left side values

   ! get number of bulk aerosols and number of modes in current list
   call rad_cnst_get_info(list_idx, naero=numaerosols, nmodes=nmodes, nbins=nbins)

   ! Contributions from modal and bin aerosols.
   if (numaerosols>0 .or. nmodes>0 .or. nbins>0) then
      call aerosol_optics_cam_sw(list_idx, state, pbuf, nnite, idxnite, &
                                 tau, tau_w, tau_w_g, tau_w_f)
   else
      tau    (1:ncol,:,:) = 0._r8
      tau_w  (1:ncol,:,:) = 0._r8
      tau_w_g(1:ncol,:,:) = 0._r8
      tau_w_f(1:ncol,:,:) = 0._r8
   end if

   !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   troplev = 0
   !REMOVECAM_END
   call tropopause_find_cam(state, troplev)


   ! diagnostic output of total aerosol optical properties
   ! currently implemented for climate list only
   call bam_optics_diags_out(lchnk, ncol, nnite, idxnite, 0, tau(:,:,idx_sw_diag), list_idx, troplev)

end subroutine aer_rad_props_sw

!==============================================================================

subroutine aer_rad_props_lw(list_idx, state, pbuf,  odap_aer)

   ! Purpose: Compute aerosol transmissions needed in absorptivity/
   !    emissivity calculations

   ! lw extinction is the same representation for all
   ! species.  If this changes, this routine will need to do something
   ! similar to the sw with routines like get_hygro_lw_abs

   use physics_buffer, only : pbuf_get_field, pbuf_get_index, physics_buffer_desc

   ! Arguments
   integer,             intent(in)  :: list_idx                      ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state

   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8),            intent(out) :: odap_aer(pcols,pver,nlwbands) ! [fraction] absorption optical depth, per layer

   ! Local variables

   integer :: bnd_idx     ! LW band index
   integer :: i           ! column index
   integer :: k           ! lev index
   integer :: ncol        ! number of columns
   integer :: numaerosols ! number of bulk aerosols in climate/diagnostic list
   integer :: nmodes      ! number of aerosol modes in climate/diagnostic list
   integer :: nbins       ! number of aerosol bins in climate/diagnostic list
   integer :: iaerosol    ! index into bulk aerosol list
   character(len=ot_length) :: opticstype       ! hygro or nonhygro

   ! optical props for each aerosol
   real(r8), pointer :: lw_abs(:)
   real(r8), pointer :: lw_hygro_abs(:,:)
   real(r8), pointer :: geometric_radius(:,:)

   ! volcanic lookup table
   real(r8), pointer :: r_lw_abs(:,:)  ! radius dependent mass-specific absorption coefficient
   real(r8), pointer :: r_mu(:)        ! log(geometric_mean_radius) domain samples of r_lw_abs(:,:)
   integer  :: idx                     ! index to pbuf for geometric radius
   real(r8) :: mu(pcols,pver)          ! log(geometric_radius)
   real(r8) :: r_mu_min, r_mu_max, wmu, mutrunc
   integer  :: nmu, kmu

   ! for table lookup into rh grid
   real(r8) :: es(pcols,pver)     ! saturation vapor pressure
   real(r8) :: qs(pcols,pver)     ! saturation specific humidity
   real(r8) :: rh(pcols,pver)
   real(r8) :: rhtrunc(pcols,pver)
   real(r8) :: wrh(pcols,pver)
   integer  :: krh(pcols,pver)

   ! aerosol (vertical) mass path and extinction
   ! aerosol masses
   real(r8), pointer :: aermmr(:,:)    ! mass mixing ratio of aerosols
   real(r8) :: mmr_to_mass(pcols,pver) ! conversion factor for mmr to mass
   real(r8) :: aermass(pcols,pver)     ! mass of aerosols

   character(len=16) :: pbuf_fld
   !-----------------------------------------------------------------------------

   ncol = state%ncol

   ! get number of bulk aerosols and number of modes in current list
   call rad_cnst_get_info(list_idx, naero=numaerosols, nmodes=nmodes, nbins=nbins)

   ! Contributions from modal and sectional aerosols.
   if (numaerosols>0 .or. nmodes>0 .or. nbins>0) then
      call aerosol_optics_cam_lw(list_idx, state, pbuf, odap_aer)
   else
      odap_aer = 0._r8
   end if

end subroutine aer_rad_props_lw

!==============================================================================
! Private methods
!==============================================================================

subroutine get_hygro_rad_props(ncol, krh, wrh, mass, ext, ssa, asm, &
                               tau, tau_w, tau_w_g, tau_w_f)

   ! Arguments
   integer,  intent(in) :: ncol
   integer,  intent(in) :: krh(pcols,pver)  ! index for linear interpolation of optics on rh
   real(r8), intent(in) :: wrh(pcols,pver)  ! weight for linear interpolation of optics on rh
   real(r8), intent(in) :: mass(pcols,pver)
   real(r8), intent(in) :: ext(:,:)
   real(r8), intent(in) :: ssa(:,:)
   real(r8), intent(in) :: asm(:,:)

   real(r8), intent(out) :: tau    (pcols,pver,nswbands)
   real(r8), intent(out) :: tau_w  (pcols,pver,nswbands)
   real(r8), intent(out) :: tau_w_g(pcols,pver,nswbands)
   real(r8), intent(out) :: tau_w_f(pcols,pver,nswbands)

   ! Local variables
   real(r8) :: ext1, ssa1, asm1
   integer :: icol, ilev, iswband
   !-----------------------------------------------------------------------------

   do iswband = 1, nswbands
      do icol = 1, ncol
         do ilev = 1, pver
            ext1 = (1 + wrh(icol,ilev)) * ext(krh(icol,ilev)+1,iswband) &
                      - wrh(icol,ilev)  * ext(krh(icol,ilev),  iswband)
            ssa1 = (1 + wrh(icol,ilev)) * ssa(krh(icol,ilev)+1,iswband) &
                      - wrh(icol,ilev)  * ssa(krh(icol,ilev),  iswband)
            asm1 = (1 + wrh(icol,ilev)) * asm(krh(icol,ilev)+1,iswband) &
                      - wrh(icol,ilev)  * asm(krh(icol,ilev),  iswband)

            tau    (icol, ilev, iswband) = mass(icol, ilev) * ext1
            tau_w  (icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1
            tau_w_g(icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1 * asm1
            tau_w_f(icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1 * asm1 * asm1
         enddo
      enddo
   enddo

end subroutine get_hygro_rad_props

!==============================================================================

subroutine get_nonhygro_rad_props(ncol, mass, ext, ssa, asm, &
                                  tau, tau_w, tau_w_g, tau_w_f)

   ! Arguments
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: mass(pcols, pver)
   real(r8), intent(in) :: ext(:)
   real(r8), intent(in) :: ssa(:)
   real(r8), intent(in) :: asm(:)

   real(r8), intent(out) :: tau    (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w  (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_g(pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_f(pcols, pver, nswbands)

   ! Local variables
   integer  :: iswband
   real(r8) :: ext1, ssa1, asm1
   !-----------------------------------------------------------------------------

   do iswband = 1, nswbands
      ext1 = ext(iswband)
      ssa1 = ssa(iswband)
      asm1 = asm(iswband)
      tau    (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1
      tau_w  (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1 * ssa1
      tau_w_g(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1 * ssa1 * asm1
      tau_w_f(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1 * ssa1 * asm1 * asm1
   enddo

end subroutine get_nonhygro_rad_props

!==============================================================================

subroutine get_volcanic_radius_rad_props(ncol, mass, pbuf_radius_name, pbuf,  r_ext, r_scat, r_ascat, r_mu, &
                                  tau, tau_w, tau_w_g, tau_w_f)


   use physics_buffer, only : pbuf_get_field, pbuf_get_index

   ! Arguments
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: mass(pcols, pver)
   character(len=*) :: pbuf_radius_name
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in) :: r_ext(:,:)
   real(r8), intent(in) :: r_scat(:,:)
   real(r8), intent(in) :: r_ascat(:,:)
   real(r8), intent(in) :: r_mu(:) ! log(radius) domain of mass-specific optics

   real(r8), intent(out) :: tau    (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w  (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_g(pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_f(pcols, pver, nswbands)

   ! Local variables
   integer  :: iswband
   real(r8) :: g

   integer  :: idx                             ! index to radius in physics buffer
   real(r8), pointer :: geometric_radius(:,:)  ! geometric mean radius of volcanic aerosol
   real(r8) :: mu(pcols,pver)                  ! log(geometric mean radius of volcanic aerosol)
   integer  :: kmu, nmu
   real(r8) :: wmu, mutrunc, r_mu_max, r_mu_min

   ! interpolated values from table
   real(r8) :: ext(nswbands)
   real(r8) :: scat(nswbands)
   real(r8) :: ascat(nswbands)

   integer :: i, k ! column level iterator
   !-----------------------------------------------------------------------------

   tau    =0._r8
   tau_w  =0._r8
   tau_w_g=0._r8
   tau_w_f=0._r8

   ! get microphysical properties for volcanic aerosols
   idx = pbuf_get_index(pbuf_radius_name)
   call pbuf_get_field(pbuf, idx, geometric_radius )

   ! interpolate in radius
   ! caution: clip the table with no warning when outside bounds
   nmu = size(r_mu)
   r_mu_max = r_mu(nmu)
   r_mu_min = r_mu(1)
   do i = 1, ncol
      do k = 1, pver
         if(geometric_radius(i,k) > 0._r8) then
            mu(i,k) = log(geometric_radius(i,k))
         else
            mu(i,k) = 0._r8
         endif
         mutrunc = max(min(mu(i,k),r_mu_max),r_mu_min)
         kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
         wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)
         do iswband = 1, nswbands
            ext(iswband) =  &
               ((1._r8 - wmu) * r_ext(iswband, kmu  ) + &
               (wmu) * r_ext(iswband, kmu+1))
            scat(iswband) =  &
               ((1._r8 - wmu) * r_scat(iswband, kmu  ) + &
               (wmu) * r_scat(iswband, kmu+1))
            ascat(iswband) =  &
               ((1._r8 - wmu) * r_ascat(iswband, kmu  ) + &
               (wmu) * r_ascat(iswband, kmu+1))
            if (scat(iswband).gt.0._r8) then
               g = ascat(iswband)/scat(iswband)
            else
               g=0._r8
            endif
            tau    (i,k,iswband) = mass(i,k) * ext(iswband)
            tau_w  (i,k,iswband) = mass(i,k) * scat(iswband)
            tau_w_g(i,k,iswband) = mass(i,k) * ascat(iswband)
            tau_w_f(i,k,iswband) = mass(i,k) * g * ascat(iswband)
         end do
      enddo
   enddo

end subroutine get_volcanic_radius_rad_props

!==============================================================================

subroutine get_volcanic_rad_props(ncol, mass, ext, scat, ascat, &
                                  tau, tau_w, tau_w_g, tau_w_f)

   ! Arguments
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: mass(pcols, pver)
   real(r8), intent(in) :: ext(:)
   real(r8), intent(in) :: scat(:)
   real(r8), intent(in) :: ascat(:)

   real(r8), intent(out) :: tau    (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w  (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_g(pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_f(pcols, pver, nswbands)

   ! Local variables
   integer  :: iswband
   real(r8) :: g
   !-----------------------------------------------------------------------------

   do iswband = 1, nswbands
      if (scat(iswband).gt.0._r8) then
         g = ascat(iswband)/scat(iswband)
      else
         g=0._r8
      endif
      tau    (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext(iswband)
      tau_w  (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * scat(iswband)
      tau_w_g(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ascat(iswband)
      tau_w_f(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * g * ascat(iswband)
   enddo

end subroutine get_volcanic_rad_props

!==============================================================================
!==============================================================================

end module aer_rad_props
