module rrtmgp_inputs

!--------------------------------------------------------------------------------
! Transform data for state inputs from CAM's data structures to those used by
! RRTMGP.  Subset the number of model levels if CAM's top exceeds RRTMGP's
! valid domain.
!
!--------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use ppgrid,           only: pcols, pver, pverp

use physconst,        only: stebol, pi

use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc
use camsrfexch,       only: cam_in_t

use radconstants,     only: nswbands, nlwbands, get_sw_spectral_boundaries
use radconstants,     only: nradgas, gaslist

use rad_constituents, only: rad_cnst_get_gas

use mcica_subcol_gen, only: mcica_subcol_sw, mcica_subcol_lw

use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp 
use mo_optical_props,      only: ty_optical_props_2str, ty_optical_props_1scl

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

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

real(r8), parameter :: cldmin = 1.0e-80_r8   ! min cloud fraction

real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

! Indices for copying data between cam and rrtmgp arrays
integer :: ktopcam ! Index in CAM arrays of top level (layer or interface) at which
                   ! RRTMGP is active.
integer :: ktoprad ! Index in RRTMGP arrays of the layer or interface corresponding
                   ! to CAM's top layer or interface

! wavenumber (cm^-1) boundaries of shortwave bands
real(r8) :: sw_low_bounds(nswbands), sw_high_bounds(nswbands)

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
   integer,                     intent(in) :: ncol      ! # cols in chunk
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
   integer :: k, kk, i, iband

   real(r8) :: tref_min, tref_max, tmin, tmax

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

   ! Assume level ordering is the same for both CAM and RRTMGP (top to bottom)
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
   end if

   ! Check that the temperatures are within the limits of RRTMGP validity.
   tref_min = kdist_sw%get_temp_min()
   tref_max = kdist_sw%get_temp_max()
   if ( any(t_rad < tref_min) .or. any(t_rad > tref_max) ) then
      ! Find out of range value and quit.
      do i = 1, ncol
         do k = 1, nlay
            if ( t_rad(i,k) < tref_min .or. t_rad(i,k) > tref_max ) then
               write(errmsg,*) 'temp outside valid range: ', t_rad(i,k), ': column lat=', &
                  state%lat(i)*180._r8/pi, ': column lon=', state%lon(i)*180._r8/pi, ': level idx=',k
               call endrun(sub//': ERROR, '//errmsg)
            end if
         end do
      end do
   end if

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
            alb_dir(iband,i) = 0.5 * (cam_in%aldir(idxday(i)) + cam_in%asdir(idxday(i)))
            alb_dif(iband,i) = 0.5 * (cam_in%aldif(idxday(i)) + cam_in%asdif(idxday(i)))
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

logical function is_visible(wavenumber)

   ! Wavenumber is in the visible if it is above the visible threshold
   ! wavenumber, and in the infrared if it is below the threshold
   ! This function doesn't distinquish between visible and UV.

   ! wavenumber in inverse cm (cm^-1)
   real(r8), intent(in) :: wavenumber

   ! Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1
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

   select case (trim(gas_name)) 
      case ('H2O') 
         massratio = 1.607793_r8
      case ('CO2')
         massratio = 0.658114_r8
      case ('O3')
         massratio = 0.603428_r8
      case ('CH4')
         massratio = 1.805423_r8
      case ('N2O')
         massratio = 0.658090_r8
      case ('O2')
         massratio = 0.905140_r8
      case ('CFC11')
         massratio = 0.210852_r8
      case ('CFC12')
         massratio = 0.239546_r8
      case default
         call endrun("Invalid gas: "//trim(gas_name))
   end select
end function get_molar_mass_ratio

subroutine rad_gas_get_vmr(icall, gas_name, state, pbuf, nlay, numactivecols, gas_concs, indices)
   ! provides volume mixing ratio into gas_concs data structure
   ! Assumes gas_name will be found with rad_cnst_get_gas(). 
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   character(len=*),            intent(in)    :: gas_name
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay           ! number of layers in radiation calculation
   integer,                     intent(in)    :: numactivecols  ! number of columns, ncol for LW, nday for SW

   type(ty_gas_concs), intent(inout)          :: gas_concs      ! the result is VRM inside gas_concs

   integer, intent(in), OPTIONAL :: indices(:) ! this would be idxday, providing the indices of the active columns

   ! local
   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)
   character(len=128)    :: errmsg
   real(r8), allocatable :: mmr(:,:)
   character(len=*), parameter :: sub = 'rad_gas_get_vmr'
   ! -- for ozone profile above model
   real(r8), allocatable :: P_int(:), P_mid(:), alpha(:), beta(:), a(:), b(:), chi_mid(:), chi_0(:), chi_eff(:)
   real(r8) :: P_top
   integer :: idx(numactivecols)
   integer :: i
   real(r8) :: alpha_value
   real(r8) :: amdo !! alpha_value of ozone


   allocate(mmr(numactivecols, nlay))
   allocate(gas_vmr(numactivecols, nlay))

   call rad_cnst_get_gas(icall, gas_name, state, pbuf, gas_mmr)
   ! copy the gas and actually convert to mmr in case of H2O (specific to mixing ratio)

   mmr = gas_mmr
   ! special case: H2O is specific humidity, not mixing ratio. Use r = q/(1-q):
   if (gas_name == 'H2O') then 
      mmr = mmr / (1._r8 - mmr)
   end if  

   ! convert MMR to VMR, multipy by ratio of dry air molar mas to gas molar mass.
   alpha_value = get_molar_mass_ratio(gas_name)

   ! set the column indices; when indices is provided (e.g. daylit columns) use them, otherwise just count.
   do i = 1,numactivecols
      if (present(indices)) then
         idx(i) = indices(i)
      else
         idx(i) = i
      end if
   end do


   if (nlay == pver) then
      do i = 1,numactivecols
         gas_vmr(i, :pver) = mmr(idx(i),:pver) * alpha_value
      end do
   else if (nlay < pver) then ! radiation calculation doesn't go through atmospheric depth
      do i = 1,numactivecols
         gas_vmr(i,nlay+1-pver:) = mmr(idx(i),:pver) * alpha_value
      end do
   else if (nlay > pver) then ! radiation has more layers than atmosphere --> only one extra layer allowed, so could say gas_vmr(:ncol, 2:) = gas_mmr(:ncol, :pver)*amdc
      do i = 1,numactivecols
         gas_vmr(i,nlay+1-pver:) = mmr(idx(i),:pver) * alpha_value
      end do
      if (nlay == pverp) then
         gas_vmr(:,1) = gas_vmr(:,nlay+1-pver) 
      else
         call endrun(sub//': Radiation can not have more than 1 extra layer.')
      end if
   end if

   ! special case: O3
   ! 
   ! """
   ! For the purpose of attenuating solar fluxes above the CAM model top, we assume that ozone 
   ! mixing decreases linearly in each column from the value in the top layer of CAM to zero at 
   ! the pressure level set by P_top. P_top has been set to 50 Pa (0.5 hPa) based on model tuning 
   ! to produce temperatures at the top of CAM that are most consistent with WACCM at similar pressure levels. 
   ! """
   if ((gas_name == 'O3') .and. (nlay == pverp)) then
      allocate(P_int(numactivecols), P_mid(numactivecols), alpha(numactivecols), beta(numactivecols), a(numactivecols), b(numactivecols), chi_mid(numactivecols), chi_0(numactivecols), chi_eff(numactivecols))
      amdo = get_molar_mass_ratio('O3')
      do i = 1, numactivecols
            P_top = 50.0_r8                     ! pressure (Pa) at which we assume O3 = 0 in linear decay from CAM top
            P_int(i) = state%pint(idx(i),1) ! pressure (Pa) at upper interface of CAM
            P_mid(i) = state%pmid(idx(i),1) ! pressure (Pa) at midpoint of top layer of CAM
            alpha(i) = 0.0_r8
            beta(i) = 0.0_r8
            alpha(i) = log(P_int(i)/P_top)
            beta(i) =  log(P_mid(i)/P_int(i))/log(P_mid(i)/P_top)
      
            a(i) =  ( (1._r8 + alpha(i)) * exp(-alpha(i)) - 1._r8 ) / alpha(i)
            b(i) =  1._r8 - exp(-alpha(i))
   
            if (alpha(i) .gt. 0) then              ! only apply where top level is below 80 km
               chi_mid(i) = mmr(i,1)*amdo          ! molar mixing ratio of O3 at midpoint of top layer
               chi_0(i) = chi_mid(i) /  (1._r8 + beta(i))
               chi_eff(i) = chi_0(i) * (a(i) + b(i))
               gas_vmr(i,1) = chi_eff(i)
               chi_eff(i) = chi_eff(i) * P_int(i) / amdo / 9.8_r8 ! O3 column above in kg m-2
               chi_eff(i) = chi_eff(i) / 2.1415e-5_r8             ! O3 column above in DU
            end if
      end do
      deallocate(P_int, P_mid, alpha, beta, a, b, chi_mid, chi_0, chi_eff)
   end if

   ! other special cases: 
   ! N2 and CO: If these are in the gas list, would set them to constants
   ! as in E3SM. Currently, these will abort run because they are not found by rad_cnst_get_gas.
   ! So while RTE-RRTMGP can cope with them, we do not use them for radiation at this time.

   errmsg = gas_concs%set_vmr(gas_name, gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR, gas_concs%set_vmr: '//trim(errmsg))
   end if

   deallocate(gas_vmr)
   deallocate(mmr)

end subroutine rad_gas_get_vmr

!==================================================================================================

subroutine rrtmgp_set_gases_lw(icall, state, pbuf, nlay, gas_concs)

   ! The gases in the LW coefficients file are:
   ! H2O, CO2, O3, N2O, CO, CH4, O2, N2
   ! But we only use the gases in the radconstants module's gaslist.

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
   integer :: ncol

   integer :: lchnk
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_lw'
   integer :: i
   !--------------------------------------------------------------------------------

   ncol = state%ncol
   lchnk = state%lchnk
   do i = 1,nradgas
      call rad_gas_get_vmr(icall, gaslist(i), state, pbuf, nlay, ncol, gas_concs)
   end do
end subroutine rrtmgp_set_gases_lw

!==================================================================================================

subroutine rrtmgp_set_gases_sw( &
   icall, state, pbuf, nlay, nday, &
   idxday, gas_concs)

   ! Return gas_concs with gas volume mixing ratio on DAYLIT columns.

   ! The gases in the SW coefficients file are:
   ! H2O, CO2, O3, N2O, CO, CH4, O2, N2, CCL4, CFC11, CFC12, CFC22, HFC143a,
   ! HFC125, HFC23, HFC32, HFC134a, CF4, NO2
   ! We only use the gases in radconstants gaslist. 

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   integer,                     intent(in)    :: nday
   integer,                     intent(in)    :: idxday(:)
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_sw'
   integer :: i

   ! use the optional argument indices to specify which columns are sunlit
    do i = 1,nradgas
      call rad_gas_get_vmr(icall, gaslist(i), state, pbuf, nlay, nday, gas_concs, indices=idxday)
   end do

end subroutine rrtmgp_set_gases_sw

!==================================================================================================

subroutine rrtmgp_set_cloud_lw(state, nlwbands, cldfrac, c_cld_lw_abs, lwkDist, cloud_lw)

   ! Create MCICA stochastic arrays for cloud LW optical properties.

   ! arguments
   type(physics_state),         intent(in)    :: state
   integer,                     intent(in)    :: nlwbands
   real(r8),                    intent(in)    :: cldfrac(pcols,pver)               ! combined cloud fraction (snow plus regular)
   real(r8),                    intent(in)    :: c_cld_lw_abs(nlwbands,pcols,pver) ! combined cloud absorption optics depth (LW)
   class(ty_gas_optics_rrtmgp), intent(in)    :: lwkDist
   type(ty_optical_props_1scl), intent(inout) :: cloud_lw
   ! local vars
   integer :: i
   integer :: ncol
   integer :: ngptlw
   real(r8), allocatable :: taucmcl(:,:,:) ! cloud optical depth [mcica]
   character(len=32)  :: sub = 'rrtmgp_set_cloud_lw'
   character(len=128) :: errmsg
   !--------------------------------------------------------------------------------
   ncol   = state%ncol
   ngptlw = lwkDist%get_ngpt()

   allocate(taucmcl(ngptlw,ncol,pver))
   
   !***NB*** this code is currently set up to create the subcols for all model layers
   !         not just the ones where the radiation calc is being done.  Need
   !         to subset cldfrac and c_cld_lw_abs to avoid computing unneeded random numbers.
   
   call mcica_subcol_lw( &
      lwkdist,      & ! spectral information
      nlwbands,     & ! number of spectral bands
      ngptlw,       & ! number of subcolumns (g-point intervals)
      ncol,         & ! number of columns
      ngptlw,       & ! changeseed, should be set to number of subcolumns
      state%pmid,   & ! layer pressures (Pa)
      cldfrac,      & ! layer cloud fraction
      c_cld_lw_abs, & ! cloud optical depth
      taucmcl       & ! OUTPUT: subcolumn cloud optical depth [mcica] (ngpt, ncol, nver)
      )

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there.
   cloud_lw%tau = 0.0_r8
   do i = 1, ngptlw
      cloud_lw%tau(:ncol, ktoprad:, i) = taucmcl(i, :ncol, ktopcam:)
   end do
   errmsg = cloud_lw%validate()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: cloud_lw%validate: '//trim(errmsg))
   end if
   deallocate(taucmcl)
end subroutine rrtmgp_set_cloud_lw

!==================================================================================================

subroutine rrtmgp_set_aer_lw(ncol, nlwbands, aer_lw_abs, aer_lw)

   ! Load aerosol optical properties into the RRTMGP object.

   ! arguments
   integer,                intent(in)    :: ncol
   integer,                intent(in)    :: nlwbands
   real(r8),               intent(in)    :: aer_lw_abs(pcols,pver,nlwbands) ! aerosol absorption optics depth (LW)
   type(ty_optical_props_1scl), intent(inout) :: aer_lw
   character(len=32)  :: sub = 'rrtmgp_set_aer_lw'
   character(len=128) :: errmsg

   !--------------------------------------------------------------------------------
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

subroutine rrtmgp_set_cloud_sw( &
   nswbands, nday, nlay, idxday, pmid, cldfrac, &
   c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, c_cld_tau_w_f, kdist_sw, &
   cloud_sw)

   ! Create MCICA stochastic arrays for cloud SW optical properties.

   ! arguments
   integer,                intent(in) :: nswbands
   integer,                intent(in) :: nday
   integer,                intent(in) :: nlay           ! number of layers in rad calc (may include "extra layer")
   integer,                intent(in) :: idxday(:)

   real(r8),               intent(in) :: pmid(nday,nlay)                    ! pressure at layer midpoints (Pa)
   real(r8),               intent(in) :: cldfrac(pcols,pver)                ! combined cloud fraction (snow plus regular)
   real(r8),               intent(in) :: c_cld_tau    (nswbands,pcols,pver) ! combined cloud extinction optical depth
   real(r8),               intent(in) :: c_cld_tau_w  (nswbands,pcols,pver) ! combined cloud single scattering albedo * tau
   real(r8),               intent(in) :: c_cld_tau_w_g(nswbands,pcols,pver) ! combined cloud assymetry parameter * w * tau
   real(r8),               intent(in) :: c_cld_tau_w_f(nswbands,pcols,pver) ! combined cloud forward scattered fraction * w * tau

   class(ty_gas_optics_rrtmgp), intent(in)    :: kdist_sw  ! shortwave gas optics object
   type(ty_optical_props_2str), intent(inout) :: cloud_sw  ! cloud optical properties object

   ! local vars
   integer, parameter :: changeseed = 1

   integer :: i, k, kk, ns, igpt
   integer :: ngptsw
   integer :: nver       ! nver is the number of cam layers in the SW calc.  It
                         ! does not include the "extra layer".

   real(r8), allocatable :: cldf(:,:)
   real(r8), allocatable :: tauc(:,:,:)
   real(r8), allocatable :: ssac(:,:,:)
   real(r8), allocatable :: asmc(:,:,:)
   real(r8), allocatable :: taucmcl(:,:,:)
   real(r8), allocatable :: ssacmcl(:,:,:)
   real(r8), allocatable :: asmcmcl(:,:,:)

   character(len=32)  :: sub = 'rrtmgp_set_cloud_sw'
   character(len=128) :: errmsg
   real(r8) :: small_val = 1.e-80_r8
   real(r8), allocatable :: day_cld_tau(:,:,:)
   real(r8), allocatable :: day_cld_tau_w(:,:,:)
   real(r8), allocatable :: day_cld_tau_w_g(:,:,:)
   !--------------------------------------------------------------------------------
   ngptsw = kdist_sw%get_ngpt()
   nver   = pver - ktopcam + 1 ! number of CAM's layers in radiation calculation. 

   ! Compute the input quantities needed for the 2-stream optical props
   ! object.  Also subset the vertical levels and the daylight columns
   ! here.  But don't reorder the vertical index because the mcica sub-column
   ! generator assumes the CAM vertical indexing.
   allocate( &
      cldf(nday,nver),           &
      tauc(nswbands,nday,nver),  &
      ssac(nswbands,nday,nver),  &
      asmc(nswbands,nday,nver),  &
      taucmcl(ngptsw,nday,nver), &
      ssacmcl(ngptsw,nday,nver), &
      asmcmcl(ngptsw,nday,nver), &
      day_cld_tau(nswbands,nday,nver),     &
      day_cld_tau_w(nswbands,nday,nver),   &
      day_cld_tau_w_g(nswbands,nday,nver))

   ! get daylit arrays on radiation levels, note: expect idxday to be truncated to size nday
   day_cld_tau     = c_cld_tau(    :, idxday(1:nday), ktopcam:)
   day_cld_tau_w   = c_cld_tau_w(  :, idxday(1:nday), ktopcam:)
   day_cld_tau_w_g = c_cld_tau_w_g(:, idxday(1:nday), ktopcam:)
   cldf = cldfrac(idxday(1:nday), ktopcam:)  ! daylit cloud fraction on radiation levels
   tauc = merge(day_cld_tau, 0.0_r8, day_cld_tau > 0.0_r8)  ! start by setting cloud optical depth, clip @ zero
   asmc = merge(day_cld_tau_w_g / max(day_cld_tau_w, small_val), 0.0_r8, day_cld_tau_w > 0.0_r8)  ! set value of asymmetry
   ssac = merge(max(day_cld_tau_w, small_val) / max(tauc, small_val), 1.0_r8 , tauc > 0.0_r8)
   asmc = merge(asmc, 0.0_r8, tauc > 0.0_r8) ! double-check asymmetry; reset when tauc = 0


   ! mcica_subcol_sw converts to gpts (e.g., 224 pts instead of 14 bands)
   ! inputs (pmid, cldf, tauc, ssac, asmc) and outputs (taucmcl, ssacmcl, asmcmcl)
   ! are on the same nver vertical levels
   ! output is shape (ngpt, ncol, nver)
   call mcica_subcol_sw( &
         kdist_sw, nswbands, ngptsw, nday, nlay, nver, changeseed, &
         pmid, cldf, tauc, ssac, asmc,     &
         taucmcl, ssacmcl, asmcmcl) ! 32
   

   ! If there is an extra layer in the radiation then this initialization
   ! will provide the optical properties there.
   ! These are shape (ncol, nlay, ngpt)
   cloud_sw%tau(:,:,:) = 0.0_r8
   cloud_sw%ssa(:,:,:) = 1.0_r8
   cloud_sw%g(:,:,:)   = 0.0_r8
   do igpt = 1,ngptsw
      cloud_sw%g  (:, ktoprad:, igpt) = asmcmcl(igpt, ktopcam:, :)
      cloud_sw%ssa(:, ktoprad:, igpt) = ssacmcl(igpt, ktopcam:, :)
      cloud_sw%tau(:, ktoprad:, igpt) = taucmcl(igpt, ktopcam:, :)
   end do


   errmsg = cloud_sw%validate()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: cloud_sw%validate: '//trim(errmsg))
   end if

   ! delta scaling adjusts for forward scattering
   ! If delta_scale() is applied, cloud_sw%tau differs from RRTMG implementation going into SW calculation.   
   errmsg = cloud_sw%delta_scale()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: cloud_sw%delta_scale: '//trim(errmsg))
   end if

   ! all information is in cloud_sw, now deallocate
   deallocate( &
      cldf, tauc, ssac, asmc, &
      taucmcl, ssacmcl, asmcmcl,&
      day_cld_tau, day_cld_tau_w, day_cld_tau_w_g )

end subroutine rrtmgp_set_cloud_sw

!==================================================================================================

subroutine rrtmgp_set_aer_sw( &
   nswbands, nday, idxday, aer_tau, aer_tau_w, &
   aer_tau_w_g, aer_tau_w_f, aer_sw)

   ! Load aerosol SW optical properties into the RRTMGP object.
   !
   ! *** N.B. *** The input optical arrays from CAM are dimensioned in the vertical
   !              as 0:pver.  The index 0 is for the extra layer used in the radiation
   !              calculation.  


   ! arguments
   integer,   intent(in) :: nswbands
   integer,   intent(in) :: nday
   integer,   intent(in) :: idxday(:)
   real(r8),  intent(in) :: aer_tau    (pcols,0:pver,nswbands) ! extinction optical depth
   real(r8),  intent(in) :: aer_tau_w  (pcols,0:pver,nswbands) ! single scattering albedo * tau
   real(r8),  intent(in) :: aer_tau_w_g(pcols,0:pver,nswbands) ! asymmetry parameter * w * tau
   real(r8),  intent(in) :: aer_tau_w_f(pcols,0:pver,nswbands) ! forward scattered fraction * w * tau
   type(ty_optical_props_2str), intent(inout) :: aer_sw

   ! local variables
   integer  :: ns
   integer  :: k, kk
   integer  :: i
   integer, dimension(nday) :: day_cols
   character(len=32)  :: sub = 'rrtmgp_set_aer_sw'
   character(len=128) :: errmsg
   !--------------------------------------------------------------------------------
   ! If there is an extra layer in the radiation then this initialization
   ! will provide default values there.
   aer_sw%tau = 0.0_r8
   aer_sw%ssa = 1.0_r8
   aer_sw%g   = 0.0_r8
   day_cols = idxday(1:nday)

   ! aer_sw is on RAD grid, aer_tau* is on CAM grid ... to make sure they align, use ktop*
   ! aer_sw has dimensions of (nday, nlay, nswbands)
   aer_sw%tau(1:nday, ktoprad:, :) = max(aer_tau(day_cols, ktopcam:, :), 0._r8)
   aer_sw%ssa(1:nday, ktoprad:, :) = merge( aer_tau_w(day_cols, ktopcam:,:)/aer_tau(day_cols, ktopcam:, :), &
                                            1._r8, aer_tau(day_cols, ktopcam:, :) > 0._r8)
   aer_sw%g(  1:nday, ktoprad:, :) = merge( aer_tau_w_g(day_cols, ktopcam:, :) / aer_tau_w(day_cols, ktopcam:, :), &
                                            0._r8, aer_tau_w(day_cols, ktopcam:, :) > 1.e-80_r8)

   ! impose limits on the components:
   ! aer_sw%tau = max(aer_sw%tau, 0._r) <-- already imposed 
   aer_sw%ssa = min(max(aer_sw%ssa, 0._r8), 1._r8)
   aer_sw%g = min(max(aer_sw%g, -1._r8), 1._r8)
   ! by clamping the values here, the validate method should be guaranteed to succeed,
   ! but we're also saying that any errors in the method to this point are being swept aside. 
   ! We might want to check for out-of-bounds values and report them in the log file.

   errmsg = aer_sw%validate()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: aer_sw%validate: '//trim(errmsg))
   end if
end subroutine rrtmgp_set_aer_sw

!==================================================================================================

end module rrtmgp_inputs
