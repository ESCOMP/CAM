module rrtmgp_inputs

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

use radconstants,     only: nradgas, gaslist, nswbands, nlwbands, nswgpts, nlwgpts,   &
                            get_sw_spectral_boundaries, idx_sw_diag, idx_sw_cloudsim, &
                            idx_lw_cloudsim

use rad_constituents, only: rad_cnst_get_gas

use cloud_rad_props,  only: get_liquid_optics_sw, liquid_cloud_get_rad_props_lw, &
                            get_ice_optics_sw,    ice_cloud_get_rad_props_lw,    &
                            get_snow_optics_sw,   snow_cloud_get_rad_props_lw,   &
                            get_grau_optics_sw,   grau_cloud_get_rad_props_lw
                                 
use mcica_subcol_gen, only: mcica_subcol_sw, mcica_subcol_lw

use aer_rad_props,    only: aer_rad_props_sw, aer_rad_props_lw

use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp 
use mo_optical_props,      only: ty_optical_props_2str, ty_optical_props_1scl

use cam_history_support,   only: fillvalue
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use error_messages,   only: alloc_err

implicit none
private
save

public :: &
   rrtmgp_inputs_init,  &
   rrtmgp_set_state,    &
   rrtmgp_set_gases_lw, &
   rrtmgp_set_gases_sw, &
   rrtmgp_set_cloud_lw, &
   rrtmgp_set_cloud_sw, &
   rrtmgp_set_aer_lw,   &
   rrtmgp_set_aer_sw


! This value is to match the arbitrary small value used in RRTMG to decide
! when a quantity is effectively zero.
real(r8), parameter :: tiny = 1.0e-80_r8

! Indices for copying data between cam and rrtmgp arrays
integer :: ktopcam ! Index in CAM arrays of top level (layer or interface) at which
                   ! RRTMGP is active.
integer :: ktoprad ! Index in RRTMGP arrays of the layer or interface corresponding
                   ! to CAM's top layer or interface

! wavenumber (cm^-1) boundaries of shortwave bands
real(r8) :: sw_low_bounds(nswbands), sw_high_bounds(nswbands)

! Mapping from RRTMG shortwave bands to RRTMGP.  Currently needed to continue using
! the SW optics datasets from RRTMG (even thought there is a slight mismatch in the
! band boundaries of the 2 bands that overlap with the LW bands).
integer, parameter, dimension(14) :: rrtmg_to_rrtmgp_swbands = &
   [ 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]

!==================================================================================================
contains
!==================================================================================================

subroutine rrtmgp_inputs_init(ktcam, ktrad)
      
   ! Note that this routine must be called after the calls to set_wavenumber_bands which set
   ! the sw/lw band boundaries in the radconstants module.

   integer, intent(in) :: ktcam
   integer, intent(in) :: ktrad

   ktopcam = ktcam
   ktoprad = ktrad

   ! Initialize the module data containing the SW band boundaries.
   call get_sw_spectral_boundaries(sw_low_bounds, sw_high_bounds, 'cm^-1')

end subroutine rrtmgp_inputs_init

!=========================================================================================

subroutine rrtmgp_set_state( &
   state, cam_in, ncol, nlay, nday,           &
   idxday, coszrs, kdist_sw, t_sfc, emis_sfc,  &
   t_rad, pmid_rad, pint_rad, t_day, pmid_day, &
   pint_day, coszrs_day, alb_dir, alb_dif) 

   ! arguments
   type(physics_state),         intent(in) :: state     ! CAM physics state
   type(cam_in_t),              intent(in) :: cam_in    ! CAM import state
   integer,                     intent(in) :: ncol      ! # cols in CAM chunk
   integer,                     intent(in) :: nlay      ! # layers in rrtmgp grid
   integer,                     intent(in) :: nday      ! # daylight columns
   integer,                     intent(in) :: idxday(:) ! chunk indicies of daylight columns
   real(r8),                    intent(in) :: coszrs(:) ! cosine of solar zenith angle
   class(ty_gas_optics_rrtmgp), intent(in) :: kdist_sw  ! spectral information

   real(r8), intent(out) :: t_sfc(ncol)              ! surface temperature [K] 
   real(r8), intent(out) :: emis_sfc(nlwbands,ncol)  ! emissivity at surface []
   real(r8), intent(out) :: t_rad(ncol,nlay)         ! layer midpoint temperatures [K]
   real(r8), intent(out) :: pmid_rad(ncol,nlay)      ! layer midpoint pressures [Pa]
   real(r8), intent(out) :: pint_rad(ncol,nlay+1)    ! layer interface pressures [Pa]
   real(r8), intent(out) :: t_day(nday,nlay)         ! layer midpoint temperatures [K]
   real(r8), intent(out) :: pmid_day(nday,nlay)      ! layer midpoint pressure [Pa]
   real(r8), intent(out) :: pint_day(nday,nlay+1)    ! layer interface pressures [Pa]
   real(r8), intent(out) :: coszrs_day(nday)         ! cosine of solar zenith angle
   real(r8), intent(out) :: alb_dir(nswbands,nday)   ! surface albedo, direct radiation
   real(r8), intent(out) :: alb_dif(nswbands,nday)   ! surface albedo, diffuse radiation

   ! local variables
   integer :: i, k, iband

   real(r8) :: tref_min, tref_max

   character(len=*), parameter :: sub='rrtmgp_set_state'
   character(len=512) :: errmsg
   !--------------------------------------------------------------------------------

   t_sfc = sqrt(sqrt(cam_in%lwup(:ncol)/stebol))  ! Surface temp set based on longwave up flux.

   ! Set surface emissivity to 1.0.
   ! The land model *does* have its own surface emissivity, but is not spectrally resolved.
   ! The LW upward flux is calculated with that land emissivity, and the "radiative temperature"
   ! t_sfc is derived from that flux. We assume, therefore, that the emissivity is unity
   ! to be consistent with t_sfc.
   emis_sfc(:,:) = 1._r8

   ! Level ordering is the same for both CAM and RRTMGP (top to bottom)
   t_rad(:,ktoprad:) = state%t(:ncol,ktopcam:)
   pmid_rad(:,ktoprad:) = state%pmid(:ncol,ktopcam:)
   pint_rad(:,ktoprad:) = state%pint(:ncol,ktopcam:)

   ! Add extra layer values if needed.
   if (nlay == pverp) then
      t_rad(:,1)      = state%t(:ncol,1)
      pmid_rad(:,1)   = 0.5_r8 * state%pint(:ncol,1)
      ! The top reference pressure from the RRTMGP coefficients datasets is 1.005183574463 Pa
      ! Set the top of the extra layer just below that.
      pint_rad(:,1) = 1.01_r8
   else
      ! nlay < pverp, thus the 1 Pa level is within a CAM layer.  Assuming the top interface of
      ! this layer is at a pressure < 1 Pa, we need to adjust the top of this layer so that it
      ! is within the valid pressure range of RRTMGP (otherwise RRTMGP issues an error).  Then
      ! set the midpoint pressure halfway between the interfaces.
      pint_rad(:,1) = 1.01_r8
      pmid_rad(:,1) = 0.5_r8 * (pint_rad(:,1) + pint_rad(:,2))
   end if

   ! Limit temperatures to be within the limits of RRTMGP validity.
   tref_min = kdist_sw%get_temp_min()
   tref_max = kdist_sw%get_temp_max()
   t_rad = merge(t_rad, tref_min, t_rad > tref_min)
   t_rad = merge(t_rad, tref_max, t_rad < tref_max)

   ! Construct arrays containing only daylight columns
   do i = 1, nday
      t_day(i,:)    = t_rad(idxday(i),:)
      pmid_day(i,:) = pmid_rad(idxday(i),:)
      pint_day(i,:) = pint_rad(idxday(i),:)
      coszrs_day(i) = coszrs(idxday(i))
   end do
 
   ! Assign albedos to the daylight columns (from E3SM implementation)
   ! Albedos are imported from the surface models as broadband (visible, and near-IR),
   ! and we need to map these to appropriate narrower bands used in RRTMGP. Bands
   ! are categorized broadly as "visible/UV" or "infrared" based on wavenumber.
   ! Loop over bands, and determine for each band whether it is broadly in the
   ! visible or infrared part of the spectrum based on a dividing line of
   ! 0.7 micron, or 14286 cm^-1
   do iband = 1,nswbands
      if (is_visible(sw_low_bounds(iband)) .and. &
         is_visible(sw_high_bounds(iband))) then

         ! Entire band is in the visible
         do i = 1, nday
            alb_dir(iband,i) = cam_in%asdir(idxday(i))
            alb_dif(iband,i) = cam_in%asdif(idxday(i))
         end do

      else if (.not.is_visible(sw_low_bounds(iband)) .and. &
               .not.is_visible(sw_high_bounds(iband))) then
         ! Entire band is in the longwave (near-infrared)
         do i = 1, nday
            alb_dir(iband,i) = cam_in%aldir(idxday(i))
            alb_dif(iband,i) = cam_in%aldif(idxday(i))
         end do
      else
         ! Band straddles the visible to near-infrared transition, so we take
         ! the albedo to be the average of the visible and near-infrared
         ! broadband albedos
         do i = 1, nday
            alb_dir(iband,i) = 0.5_r8 * (cam_in%aldir(idxday(i)) + cam_in%asdir(idxday(i)))
            alb_dif(iband,i) = 0.5_r8 * (cam_in%aldif(idxday(i)) + cam_in%asdif(idxday(i)))
         end do
      end if
   end do

   ! Strictly enforce albedo bounds
   where (alb_dir < 0)
       alb_dir = 0.0_r8
   end where
   where (alb_dir > 1)
       alb_dir = 1.0_r8
   end where
   where (alb_dif < 0)
       alb_dif = 0.0_r8
   end where
   where (alb_dif > 1)
       alb_dif = 1.0_r8
   end where

end subroutine rrtmgp_set_state

!=========================================================================================

pure logical function is_visible(wavenumber)

   ! Wavenumber is in the visible if it is above the visible threshold
   ! wavenumber, and in the infrared if it is below the threshold
   ! This function doesn't distinquish between visible and UV.

   ! wavenumber in inverse cm (cm^-1)
   real(r8), intent(in) :: wavenumber

   ! Set threshold between visible and infrared to 0.7 micron, or 14286 cm^-1
   real(r8), parameter :: visible_wavenumber_threshold = 14286._r8  ! cm^-1

   if (wavenumber > visible_wavenumber_threshold) then
      is_visible = .true.
   else
      is_visible = .false.
   end if

end function is_visible

!=========================================================================================

function get_molar_mass_ratio(gas_name) result(massratio)

   ! return the molar mass ratio of dry air to gas based on gas_name

   character(len=*),intent(in) :: gas_name
   real(r8)                    :: massratio

   ! local variables
   real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
   real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
   real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
   real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
   real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
   real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
   real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
   real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

   character(len=*), parameter :: sub='get_molar_mass_ratio'
   !----------------------------------------------------------------------------

   select case (trim(gas_name)) 
      case ('H2O') 
         massratio = amdw
      case ('CO2')
         massratio = amdc
      case ('O3')
         massratio = amdo
      case ('CH4')
         massratio = amdm
      case ('N2O')
         massratio = amdn
      case ('O2')
         massratio = amdo2
      case ('CFC11')
         massratio = amdc1
      case ('CFC12')
         massratio = amdc2
      case default
         call endrun(sub//": Invalid gas: "//trim(gas_name))
   end select

end function get_molar_mass_ratio

!=========================================================================================

subroutine rad_gas_get_vmr(icall, gas_name, state, pbuf, nlay, numactivecols, gas_concs, idxday)

   ! Set volume mixing ratio in gas_concs object.
   ! The gas_concs%set_vmr method copies data into internally allocated storage.

   integer,                     intent(in) :: icall      ! index of climate/diagnostic radiation call
   character(len=*),            intent(in) :: gas_name
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   integer,                     intent(in) :: nlay           ! number of layers in radiation calculation
   integer,                     intent(in) :: numactivecols  ! number of columns, ncol for LW, nday for SW

   type(ty_gas_concs),       intent(inout) :: gas_concs  ! the result is VRM inside gas_concs

   integer, optional,          intent(in) :: idxday(:)   ! indices of daylight columns in a chunk

   ! Local variables
   integer :: i, idx(numactivecols)
   integer :: istat
   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)
   real(r8), allocatable :: mmr(:,:)
   real(r8) :: massratio

   ! For ozone profile above model
   real(r8) :: P_top, P_int, P_mid, alpha, beta, a, b, chi_mid, chi_0, chi_eff

   character(len=128)          :: errmsg
   character(len=*), parameter :: sub = 'rad_gas_get_vmr'
   !----------------------------------------------------------------------------

   ! set the column indices; when idxday is provided (e.g. daylit columns) use them, otherwise just count.
   do i = 1, numactivecols
      if (present(idxday)) then
         idx(i) = idxday(i)
      else
         idx(i) = i
      end if
   end do

   ! gas_mmr points to a "chunk" in either the state or pbuf objects.  That storage is
   ! dimensioned (pcols,pver).
   call rad_cnst_get_gas(icall, gas_name, state, pbuf, gas_mmr)

   ! Copy into storage for RRTMGP
   allocate(mmr(numactivecols, nlay), stat=istat)
   call alloc_err(istat, sub, 'mmr', numactivecols*nlay)
   allocate(gas_vmr(numactivecols, nlay), stat=istat)
   call alloc_err(istat, sub, 'gas_vmr', numactivecols*nlay)

   do i = 1, numactivecols
      mmr(i,ktoprad:) = gas_mmr(idx(i),ktopcam:)
   end do

   ! If an extra layer is being used, copy mmr from the top layer of CAM to the extra layer.
   if (nlay == pverp) then
      mmr(:,1) = mmr(:,2)
   end if

   ! special case: H2O is specific humidity, not mixing ratio. Use r = q/(1-q):
   if (gas_name == 'H2O') then 
      mmr = mmr / (1._r8 - mmr)
   end if  

   ! convert MMR to VMR, multipy by ratio of dry air molar mas to gas molar mass.
   massratio = get_molar_mass_ratio(gas_name)
   gas_vmr = mmr * massratio

   ! special case: Setting O3 in the extra layer:
   ! 
   ! For the purpose of attenuating solar fluxes above the CAM model top, we assume that ozone 
   ! mixing decreases linearly in each column from the value in the top layer of CAM to zero at 
   ! the pressure level set by P_top. P_top has been set to 50 Pa (0.5 hPa) based on model tuning 
   ! to produce temperatures at the top of CAM that are most consistent with WACCM at similar pressure levels. 

   if ((gas_name == 'O3') .and. (nlay == pverp)) then
      P_top = 50.0_r8
      do i = 1, numactivecols
            P_int = state%pint(idx(i),1) ! pressure (Pa) at upper interface of CAM
            P_mid = state%pmid(idx(i),1) ! pressure (Pa) at midpoint of top layer of CAM
            alpha = log(P_int/P_top)
            beta =  log(P_mid/P_int)/log(P_mid/P_top)
      
            a =  ( (1._r8 + alpha) * exp(-alpha) - 1._r8 ) / alpha
            b =  1._r8 - exp(-alpha)
   
            if (alpha .gt. 0) then             ! only apply where top level is below 80 km
               chi_mid = gas_vmr(i,1)          ! molar mixing ratio of O3 at midpoint of top layer
               chi_0 = chi_mid /  (1._r8 + beta)
               chi_eff = chi_0 * (a + b)
               gas_vmr(i,1) = chi_eff
            end if
      end do
   end if

   errmsg = gas_concs%set_vmr(gas_name, gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR, gas_concs%set_vmr: '//trim(errmsg))
   end if

   deallocate(gas_vmr)
   deallocate(mmr)

end subroutine rad_gas_get_vmr

!==================================================================================================

subroutine rrtmgp_set_gases_lw(icall, state, pbuf, nlay, gas_concs)

   ! Set gas vmr for the gases in the radconstants module's gaslist.

   ! The memory management for the gas_concs object is internal.  The arrays passed to it
   ! are copied to the internally allocated memory.  Each call to the set_vmr method checks
   ! whether the gas already has memory allocated, and if it does that memory is deallocated
   ! and new memory is allocated.

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   integer :: i, ncol
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_lw'
   !--------------------------------------------------------------------------------

   ncol = state%ncol
   do i = 1, nradgas
      call rad_gas_get_vmr(icall, gaslist(i), state, pbuf, nlay, ncol, gas_concs)
   end do
end subroutine rrtmgp_set_gases_lw

!==================================================================================================

subroutine rrtmgp_set_gases_sw( &
   icall, state, pbuf, nlay, nday, &
   idxday, gas_concs)

   ! Return gas_concs with gas volume mixing ratio on DAYLIT columns.
   ! Set all gases in radconstants gaslist.

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   integer,                     intent(in)    :: nday
   integer,                     intent(in)    :: idxday(:)
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   integer :: i
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_sw'
   !----------------------------------------------------------------------------

   ! use the optional argument idxday to specify which columns are sunlit
    do i = 1,nradgas
      call rad_gas_get_vmr(icall, gaslist(i), state, pbuf, nlay, nday, gas_concs, idxday=idxday)
   end do

end subroutine rrtmgp_set_gases_sw

!==================================================================================================

subroutine rrtmgp_set_cloud_lw( &
   state, pbuf, ncol, nlay, nlaycam, &
   cld, cldfsnow, cldfgrau, cldfprime, graupel_in_rad, &
   kdist_lw, cloud_lw, cld_lw_abs_cloudsim, snow_lw_abs_cloudsim, grau_lw_abs_cloudsim)

   ! Compute combined cloud optical properties.
   ! Create MCICA stochastic arrays for cloud LW optical properties.
   ! Initialize optical properties object (cloud_lw) and load with MCICA columns.

   ! arguments
   type(physics_state),         intent(in)  :: state
   type(physics_buffer_desc),   pointer     :: pbuf(:)
   integer,  intent(in) :: ncol           ! number of columns in CAM chunk
   integer,  intent(in) :: nlay           ! number of layers in radiation calculation (may include "extra layer")
   integer,  intent(in) :: nlaycam        ! number of CAM layers in radiation calculation
   real(r8), pointer    :: cld(:,:)       ! cloud fraction (liq+ice)
   real(r8), pointer    :: cldfsnow(:,:)  ! cloud fraction of just "snow clouds"
   real(r8), pointer    :: cldfgrau(:,:)  ! cloud fraction of just "graupel clouds"
   real(r8), intent(in) :: cldfprime(pcols,pver) ! combined cloud fraction

   logical,                     intent(in)  :: graupel_in_rad ! use graupel in radiation code
   class(ty_gas_optics_rrtmgp), intent(in)  :: kdist_lw
   type(ty_optical_props_1scl), intent(out) :: cloud_lw

   ! Diagnostic outputs
   real(r8), intent(out) :: cld_lw_abs_cloudsim(pcols,pver) ! in-cloud liq+ice optical depth (for COSP)
   real(r8), intent(out) :: snow_lw_abs_cloudsim(pcols,pver)! in-cloud snow optical depth (for COSP)
   real(r8), intent(out) :: grau_lw_abs_cloudsim(pcols,pver)! in-cloud Graupel optical depth (for COSP)

   ! Local variables

   integer :: i, k

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: liq_lw_abs(nlwbands,pcols,pver)   ! liquid absorption optics depth (LW)
   real(r8) :: ice_lw_abs(nlwbands,pcols,pver)   ! ice absorption optics depth (LW)
   real(r8) :: cld_lw_abs(nlwbands,pcols,pver)   ! cloud absorption optics depth (LW)
   real(r8) :: snow_lw_abs(nlwbands,pcols,pver)  ! snow absorption optics depth (LW)
   real(r8) :: grau_lw_abs(nlwbands,pcols,pver)  ! graupel absorption optics depth (LW)
   real(r8) :: c_cld_lw_abs(nlwbands,pcols,pver) ! combined cloud absorption optics depth (LW)

   ! Arrays for converting from CAM chunks to RRTMGP inputs.
   real(r8) :: cldf(ncol,nlaycam)
   real(r8) :: tauc(nlwbands,ncol,nlaycam)
   real(r8) :: taucmcl(nlwgpts,ncol,nlaycam)

   character(len=128) :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_cloud_lw'
   !--------------------------------------------------------------------------------

   ! Combine the cloud optical properties.  These calculations are done on CAM "chunks".

   ! gammadist liquid optics
   call liquid_cloud_get_rad_props_lw(state, pbuf, liq_lw_abs)
   ! Mitchell ice optics
   call ice_cloud_get_rad_props_lw(state, pbuf, ice_lw_abs)

   cld_lw_abs(:,:ncol,:) = liq_lw_abs(:,:ncol,:) + ice_lw_abs(:,:ncol,:)

   if (associated(cldfsnow)) then
      ! add in snow
      call snow_cloud_get_rad_props_lw(state, pbuf, snow_lw_abs)
      do i = 1, ncol
         do k = 1, pver
            if (cldfprime(i,k) > 0._r8) then
               c_cld_lw_abs(:,i,k) = ( cldfsnow(i,k)*snow_lw_abs(:,i,k) &
                                            + cld(i,k)*cld_lw_abs(:,i,k) )/cldfprime(i,k)
            else
               c_cld_lw_abs(:,i,k) = 0._r8
            end if
         end do
      end do
   else
      c_cld_lw_abs(:,:ncol,:) = cld_lw_abs(:,:ncol,:)
   end if

   ! add in graupel
   if (associated(cldfgrau) .and. graupel_in_rad) then
      call grau_cloud_get_rad_props_lw(state, pbuf, grau_lw_abs)
      do i = 1, ncol
         do k = 1, pver
            if (cldfprime(i,k) > 0._r8) then
               c_cld_lw_abs(:,i,k) = ( cldfgrau(i,k)*grau_lw_abs(:,i,k) &
                                            + cld(i,k)*c_cld_lw_abs(:,i,k) )/cldfprime(i,k)
            else
               c_cld_lw_abs(:,i,k) = 0._r8
            end if
         end do
      end do
   end if

   ! Cloud optics for COSP
   cld_lw_abs_cloudsim  = cld_lw_abs(idx_lw_cloudsim,:,:)
   snow_lw_abs_cloudsim = snow_lw_abs(idx_lw_cloudsim,:,:)
   grau_lw_abs_cloudsim = grau_lw_abs(idx_lw_cloudsim,:,:)
   
   ! Extract just the layers of CAM where RRTMGP does calculations.

   ! Subset "chunk" data so just the number of CAM layers in the
   ! radiation calculation are used by MCICA to produce subcolumns.
   cldf = cldfprime(:ncol, ktopcam:)
   tauc = c_cld_lw_abs(:, :ncol, ktopcam:)

   ! Enforce tauc >= 0.
   tauc = merge(tauc, 0.0_r8, tauc > 0.0_r8)
   
   ! MCICA uses spectral data (on bands) to construct subcolumns (one per g-point)
   call mcica_subcol_lw( &
      kdist_lw, nlwbands, nlwgpts, ncol, nlaycam, &
      nlwgpts, state%pmid, cldf, tauc, taucmcl    )

   errmsg =cloud_lw%alloc_1scl(ncol, nlay, kdist_lw)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: cloud_lw%alloc_1scalar: '//trim(errmsg))
   end if

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there.
   cloud_lw%tau = 0.0_r8

   ! Set the properties on g-points.
   do i = 1, nlwgpts
      cloud_lw%tau(:,ktoprad:,i) = taucmcl(i,:,:)
   end do

   ! validate checks that: tau > 0
   errmsg = cloud_lw%validate()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: cloud_lw%validate: '//trim(errmsg))
   end if

end subroutine rrtmgp_set_cloud_lw

!==================================================================================================

subroutine rrtmgp_set_cloud_sw( &
   state, pbuf, nlay, nday, idxday, &
   nnite, idxnite, pmid, cld, cldfsnow, &
   cldfgrau, cldfprime, graupel_in_rad, kdist_sw, cloud_sw, &
   tot_cld_vistau, tot_icld_vistau, liq_icld_vistau, ice_icld_vistau, snow_icld_vistau, &
   grau_icld_vistau, cld_tau_cloudsim, snow_tau_cloudsim, grau_tau_cloudsim)

   ! Compute combined cloud optical properties.
   ! Create MCICA stochastic arrays for cloud SW optical properties.
   ! Initialize optical properties object (cloud_sw) and load with MCICA columns.

   ! arguments
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   integer,  intent(in) :: nlay           ! number of layers in radiation calculation (may include "extra layer")
   integer,  intent(in) :: nday           ! number of daylight columns
   integer,  intent(in) :: idxday(pcols)  ! indices of daylight columns in the chunk
   integer,  intent(in) :: nnite          ! number of night columns
   integer,  intent(in) :: idxnite(pcols) ! indices of night columns in the chunk

   real(r8), intent(in) :: pmid(nday,nlay)! pressure at layer midpoints (Pa) used to seed RNG.

   real(r8), pointer    :: cld(:,:)      ! cloud fraction (liq+ice)
   real(r8), pointer    :: cldfsnow(:,:) ! cloud fraction of just "snow clouds"
   real(r8), pointer    :: cldfgrau(:,:) ! cloud fraction of just "graupel clouds"
   real(r8), intent(in) :: cldfprime(pcols,pver) ! combined cloud fraction

   logical,                     intent(in)  :: graupel_in_rad ! graupel in radiation code
   class(ty_gas_optics_rrtmgp), intent(in)  :: kdist_sw  ! shortwave gas optics object
   type(ty_optical_props_2str), intent(out) :: cloud_sw  ! SW cloud optical properties object

   ! Diagnostic outputs
   real(r8), intent(out) :: tot_cld_vistau(pcols,pver)   ! gbx total cloud optical depth
   real(r8), intent(out) :: tot_icld_vistau(pcols,pver)  ! in-cld total cloud optical depth
   real(r8), intent(out) :: liq_icld_vistau(pcols,pver)  ! in-cld liq cloud optical depth
   real(r8), intent(out) :: ice_icld_vistau(pcols,pver)  ! in-cld ice cloud optical depth
   real(r8), intent(out) :: snow_icld_vistau(pcols,pver) ! snow in-cloud visible sw optical depth
   real(r8), intent(out) :: grau_icld_vistau(pcols,pver) ! Graupel in-cloud visible sw optical depth
   real(r8), intent(out) :: cld_tau_cloudsim(pcols,pver) ! in-cloud liq+ice optical depth (for COSP)
   real(r8), intent(out) :: snow_tau_cloudsim(pcols,pver)! in-cloud snow optical depth (for COSP)
   real(r8), intent(out) :: grau_tau_cloudsim(pcols,pver)! in-cloud Graupel optical depth (for COSP)

   ! Local variables

   integer :: i, k, ncol
   integer :: igpt, nver
   integer :: istat
   integer, parameter :: changeseed = 1

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: ice_tau    (nswbands,pcols,pver) ! ice extinction optical depth
   real(r8) :: ice_tau_w  (nswbands,pcols,pver) ! ice single scattering albedo * tau
   real(r8) :: ice_tau_w_g(nswbands,pcols,pver) ! ice asymmetry parameter * tau * w
   real(r8) :: liq_tau    (nswbands,pcols,pver) ! liquid extinction optical depth
   real(r8) :: liq_tau_w  (nswbands,pcols,pver) ! liquid single scattering albedo * tau
   real(r8) :: liq_tau_w_g(nswbands,pcols,pver) ! liquid asymmetry parameter * tau * w
   real(r8) :: cld_tau    (nswbands,pcols,pver) ! cloud extinction optical depth
   real(r8) :: cld_tau_w  (nswbands,pcols,pver) ! cloud single scattering albedo * tau
   real(r8) :: cld_tau_w_g(nswbands,pcols,pver) ! cloud asymmetry parameter * w * tau
   real(r8) :: snow_tau    (nswbands,pcols,pver) ! snow extinction optical depth
   real(r8) :: snow_tau_w  (nswbands,pcols,pver) ! snow single scattering albedo * tau
   real(r8) :: snow_tau_w_g(nswbands,pcols,pver) ! snow asymmetry parameter * tau * w
   real(r8) :: grau_tau    (nswbands,pcols,pver) ! graupel extinction optical depth
   real(r8) :: grau_tau_w  (nswbands,pcols,pver) ! graupel single scattering albedo * tau
   real(r8) :: grau_tau_w_g(nswbands,pcols,pver) ! graupel asymmetry parameter * tau * w
   real(r8) :: c_cld_tau    (nswbands,pcols,pver) ! combined cloud extinction optical depth
   real(r8) :: c_cld_tau_w  (nswbands,pcols,pver) ! combined cloud single scattering albedo * tau
   real(r8) :: c_cld_tau_w_g(nswbands,pcols,pver) ! combined cloud asymmetry parameter * w * tau

   ! RRTMGP does not use this property in its 2-stream calculations.
   real(r8) :: sw_tau_w_f(nswbands,pcols,pver) ! Forward scattered fraction * tau * w.

   ! Arrays for converting from CAM chunks to RRTMGP inputs.
   real(r8), allocatable :: cldf(:,:)
   real(r8), allocatable :: tauc(:,:,:)
   real(r8), allocatable :: ssac(:,:,:)
   real(r8), allocatable :: asmc(:,:,:)
   real(r8), allocatable :: taucmcl(:,:,:)
   real(r8), allocatable :: ssacmcl(:,:,:)
   real(r8), allocatable :: asmcmcl(:,:,:)
   real(r8), allocatable :: day_cld_tau(:,:,:)
   real(r8), allocatable :: day_cld_tau_w(:,:,:)
   real(r8), allocatable :: day_cld_tau_w_g(:,:,:)

   character(len=128) :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_cloud_sw'
   !--------------------------------------------------------------------------------

   ncol = state%ncol

   ! Combine the cloud optical properties.  These calculations are done on CAM "chunks".

   ! gammadist liquid optics
   call get_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, sw_tau_w_f)
   ! Mitchell ice optics
   call get_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, sw_tau_w_f)

   cld_tau(:,:ncol,:)     =  liq_tau(:,:ncol,:)     + ice_tau(:,:ncol,:)
   cld_tau_w(:,:ncol,:)   =  liq_tau_w(:,:ncol,:)   + ice_tau_w(:,:ncol,:)
   cld_tau_w_g(:,:ncol,:) =  liq_tau_w_g(:,:ncol,:) + ice_tau_w_g(:,:ncol,:)

   ! add in snow
   if (associated(cldfsnow)) then
      call get_snow_optics_sw(state, pbuf, snow_tau, snow_tau_w, snow_tau_w_g, sw_tau_w_f)
      do i = 1, ncol
         do k = 1, pver
            if (cldfprime(i,k) > 0._r8) then
               c_cld_tau(:,i,k)     = ( cldfsnow(i,k)*snow_tau(:,i,k) &
                                      + cld(i,k)*cld_tau(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w(:,i,k)   = ( cldfsnow(i,k)*snow_tau_w(:,i,k)  &
                                      + cld(i,k)*cld_tau_w(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w_g(:,i,k) = ( cldfsnow(i,k)*snow_tau_w_g(:,i,k) &
                                      + cld(i,k)*cld_tau_w_g(:,i,k) )/cldfprime(i,k)
            else
               c_cld_tau(:,i,k)     = 0._r8
               c_cld_tau_w(:,i,k)   = 0._r8
               c_cld_tau_w_g(:,i,k) = 0._r8
            end if
         end do
      end do
   else
      c_cld_tau(:,:ncol,:)     = cld_tau(:,:ncol,:)
      c_cld_tau_w(:,:ncol,:)   = cld_tau_w(:,:ncol,:)
      c_cld_tau_w_g(:,:ncol,:) = cld_tau_w_g(:,:ncol,:)
   end if

   ! add in graupel
   if (associated(cldfgrau) .and. graupel_in_rad) then
      call get_grau_optics_sw(state, pbuf, grau_tau, grau_tau_w, grau_tau_w_g, sw_tau_w_f)
      do i = 1, ncol
         do k = 1, pver
            if (cldfprime(i,k) > 0._r8) then
               c_cld_tau(:,i,k)     = ( cldfgrau(i,k)*grau_tau(:,i,k) &
                                      + cld(i,k)*c_cld_tau(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w(:,i,k)   = ( cldfgrau(i,k)*grau_tau_w(:,i,k)  &
                                      + cld(i,k)*c_cld_tau_w(:,i,k) )/cldfprime(i,k)
               c_cld_tau_w_g(:,i,k) = ( cldfgrau(i,k)*grau_tau_w_g(:,i,k) &
                                      + cld(i,k)*c_cld_tau_w_g(:,i,k) )/cldfprime(i,k)
            else
               c_cld_tau(:,i,k)     = 0._r8
               c_cld_tau_w(:,i,k)   = 0._r8
               c_cld_tau_w_g(:,i,k) = 0._r8
            end if
         end do
      end do
   end if

   ! cloud optical properties need to be re-ordered from the RRTMG spectral bands
   ! (assumed in the optics datasets) to RRTMGP's
   ice_tau(:,:ncol,:)       = ice_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   liq_tau(:,:ncol,:)       = liq_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   c_cld_tau(:,:ncol,:)     = c_cld_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   c_cld_tau_w(:,:ncol,:)   = c_cld_tau_w(rrtmg_to_rrtmgp_swbands,:ncol,:)
   c_cld_tau_w_g(:,:ncol,:) = c_cld_tau_w_g(rrtmg_to_rrtmgp_swbands,:ncol,:)
   if (associated(cldfsnow)) then
      snow_tau(:,:ncol,:)   = snow_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   end if
   if (associated(cldfgrau) .and. graupel_in_rad) then
      grau_tau(:,:ncol,:)   = grau_tau(rrtmg_to_rrtmgp_swbands,:ncol,:)
   end if

   ! Set arrays for diagnostic output.
   ! cloud optical depth fields for the visible band
   tot_icld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)
   liq_icld_vistau(:ncol,:) = liq_tau(idx_sw_diag,:ncol,:)
   ice_icld_vistau(:ncol,:) = ice_tau(idx_sw_diag,:ncol,:)
   if (associated(cldfsnow)) then
      snow_icld_vistau(:ncol,:) = snow_tau(idx_sw_diag,:ncol,:)
   endif
   if (associated(cldfgrau) .and. graupel_in_rad) then
      grau_icld_vistau(:ncol,:) = grau_tau(idx_sw_diag,:ncol,:)
   endif

   ! multiply by total cloud fraction to get gridbox value
   tot_cld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)*cldfprime(:ncol,:)

   ! overwrite night columns with fillvalue
   do i = 1, Nnite
      tot_cld_vistau(IdxNite(i),:)   = fillvalue
      tot_icld_vistau(IdxNite(i),:)  = fillvalue
      liq_icld_vistau(IdxNite(i),:)  = fillvalue
      ice_icld_vistau(IdxNite(i),:)  = fillvalue
      if (associated(cldfsnow)) then
         snow_icld_vistau(IdxNite(i),:) = fillvalue
      end if
      if (associated(cldfgrau) .and. graupel_in_rad) then
         grau_icld_vistau(IdxNite(i),:) = fillvalue
      end if
   end do

   ! Cloud optics for COSP
   cld_tau_cloudsim = cld_tau(idx_sw_cloudsim,:,:)
   snow_tau_cloudsim = snow_tau(idx_sw_cloudsim,:,:)
   grau_tau_cloudsim = grau_tau(idx_sw_cloudsim,:,:)

   ! if no daylight columns the cloud_sw object isn't initialized
   if (nday > 0) then

      ! number of CAM's layers in radiation calculation.  Does not include the "extra layer".
      nver = pver - ktopcam + 1

      allocate( &
         cldf(nday,nver),                      &
         day_cld_tau(nswbands,nday,nver),      &
         day_cld_tau_w(nswbands,nday,nver),    &
         day_cld_tau_w_g(nswbands,nday,nver),  &
         tauc(nswbands,nday,nver), taucmcl(nswgpts,nday,nver), &
         ssac(nswbands,nday,nver), ssacmcl(nswgpts,nday,nver), &
         asmc(nswbands,nday,nver), asmcmcl(nswgpts,nday,nver), stat=istat)
      call alloc_err(istat, sub, 'cldf,..,asmcmcl', 9*nswgpts*nday*nver)

      ! Subset "chunk" data so just the daylight columns, and the number of CAM layers in the
      ! radiation calculation are used by MCICA to produce subcolumns.
      cldf            = cldfprime(       idxday(1:nday), ktopcam:)
      day_cld_tau     = c_cld_tau(    :, idxday(1:nday), ktopcam:)
      day_cld_tau_w   = c_cld_tau_w(  :, idxday(1:nday), ktopcam:)
      day_cld_tau_w_g = c_cld_tau_w_g(:, idxday(1:nday), ktopcam:)

      ! Compute the optical properties needed for the 2-stream calculations.  These calculations
      ! are the same as the RRTMG version.

      ! set cloud optical depth, clip @ zero
      tauc = merge(day_cld_tau, 0.0_r8, day_cld_tau > 0.0_r8)
      ! set value of asymmetry
      asmc = merge(day_cld_tau_w_g / max(day_cld_tau_w, tiny), 0.0_r8, day_cld_tau_w > 0.0_r8)
      ! set value of single scattering albedo
      ssac = merge(max(day_cld_tau_w, tiny) / max(tauc, tiny), 1.0_r8 , tauc > 0.0_r8)
      ! set asymmetry to zero when tauc = 0
      asmc = merge(asmc, 0.0_r8, tauc > 0.0_r8)

      ! MCICA uses spectral data (on bands) to construct subcolumns (one per g-point)
      call mcica_subcol_sw( &
         kdist_sw, nswbands, nswgpts, nday, nlay, &
         nver, changeseed, pmid, cldf, tauc,     &
         ssac, asmc, taucmcl, ssacmcl, asmcmcl)
   
      ! Initialize object for SW cloud optical properties.
      errmsg = cloud_sw%alloc_2str(nday, nlay, kdist_sw)
      if (len_trim(errmsg) > 0) then
         call endrun(trim(sub)//': ERROR: cloud_sw%alloc_2str: '//trim(errmsg))
      end if

      ! If there is an extra layer in the radiation then this initialization
      ! will provide the optical properties there.
      cloud_sw%tau = 0.0_r8
      cloud_sw%ssa = 1.0_r8
      cloud_sw%g   = 0.0_r8

      ! Set the properties on g-points.
      do igpt = 1,nswgpts
         cloud_sw%g  (:,ktoprad:,igpt) = asmcmcl(igpt,:,:)
         cloud_sw%ssa(:,ktoprad:,igpt) = ssacmcl(igpt,:,:)
         cloud_sw%tau(:,ktoprad:,igpt) = taucmcl(igpt,:,:)
      end do

      ! validate checks that: tau > 0, ssa is in range [0,1], and g is in range [-1,1].
      errmsg = cloud_sw%validate()
      if (len_trim(errmsg) > 0) then
         call endrun(sub//': ERROR: cloud_sw%validate: '//trim(errmsg))
      end if

      ! delta scaling adjusts for forward scattering
      errmsg = cloud_sw%delta_scale()
      if (len_trim(errmsg) > 0) then
         call endrun(sub//': ERROR: cloud_sw%delta_scale: '//trim(errmsg))
      end if

      ! All information is in cloud_sw, now deallocate local vars.
      deallocate( &
         cldf, tauc, ssac, asmc, &
         taucmcl, ssacmcl, asmcmcl,&
         day_cld_tau, day_cld_tau_w, day_cld_tau_w_g )

   end if

end subroutine rrtmgp_set_cloud_sw

!==================================================================================================

subroutine rrtmgp_set_aer_lw(icall, state, pbuf, aer_lw)

   ! Load LW aerosol optical properties into the RRTMGP object.

   ! Arguments
   integer,                     intent(in) :: icall
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   type(ty_optical_props_1scl), intent(inout) :: aer_lw

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
   aer_lw%tau = 0.0_r8

   aer_lw%tau(:ncol, ktoprad:, :) = aer_lw_abs(:ncol, ktopcam:, :)

   errmsg = aer_lw%validate()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: aer_lw%validate: '//trim(errmsg))
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

   type(ty_optical_props_2str), intent(inout) :: aer_sw

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
      aer_sw%tau = 0.0_r8
      aer_sw%ssa = 1.0_r8
      aer_sw%g   = 0.0_r8

      ! CAM fields are products tau, tau*ssa, tau*ssa*asy
      ! Fields expected by RRTMGP are computed by
      ! aer_sw%tau = aer_tau
      ! aer_sw%ssa = aer_tau_w / aer_tau
      ! aer_sw%g   = aer_tau_w_g / aer_taw_w
      ! aer_sw arrays have dimensions of (nday,nlay,nswbands)

      do i = 1, nday
         ! set aerosol optical depth, clip to zero
         aer_sw%tau(i,ktoprad:,:) = max(aer_tau(idxday(i),ktopcam:,:), 0._r8)
         ! set value of single scattering albedo
         aer_sw%ssa(i,ktoprad:,:) = merge(aer_tau_w(idxday(i),ktopcam:,:)/aer_tau(idxday(i),ktopcam:,:), &
                                          1._r8, aer_tau(idxday(i),ktopcam:,:) > 0._r8)
         ! set value of asymmetry
         aer_sw%g(i,ktoprad:,:) = merge(aer_tau_w_g(idxday(i),ktopcam:,:)/aer_tau_w(idxday(i),ktopcam:,:), &
                                        0._r8, aer_tau_w(idxday(i),ktopcam:,:) > tiny)
      end do

      ! impose limits on the components
      aer_sw%ssa = min(max(aer_sw%ssa, 0._r8), 1._r8)
      aer_sw%g   = min(max(aer_sw%g, -1._r8), 1._r8)

   end if

end subroutine rrtmgp_set_aer_sw

!==================================================================================================

end module rrtmgp_inputs
