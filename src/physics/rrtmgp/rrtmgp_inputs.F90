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
   end if

   ! Check that the temperatures are within the limits of RRTMGP validity.
   tref_min = kdist_sw%get_temp_min()
   tref_max = kdist_sw%get_temp_max()
   if ( any(t_rad < tref_min) .or. any(t_rad > tref_max) ) then
      ! Report out of range value and quit.
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

   ! Set volume mixing ratio in gas_concs data structure.

   integer,                     intent(in) :: icall      ! index of climate/diagnostic radiation call
   character(len=*),            intent(in) :: gas_name
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   integer,                     intent(in) :: nlay           ! number of layers in radiation calculation
   integer,                     intent(in) :: numactivecols  ! number of columns, ncol for LW, nday for SW

   type(ty_gas_concs),       intent(inout) :: gas_concs  ! the result is VRM inside gas_concs

   integer, optional,          intent(in) :: idxday(:)   ! indices of daylight columns in a chunk

   ! local
   integer :: i, idx(numactivecols)
   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)
   real(r8), allocatable :: mmr(:,:)
   real(r8) :: massratio

   ! -- for ozone profile above model
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
   allocate(mmr(numactivecols, nlay))
   allocate(gas_vmr(numactivecols, nlay))

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
      do i = 1, numactivecols
            P_top = 50.0_r8
            P_int = state%pint(idx(i),1) ! pressure (Pa) at upper interface of CAM
            P_mid = state%pmid(idx(i),1) ! pressure (Pa) at midpoint of top layer of CAM
            alpha = 0.0_r8
            beta = 0.0_r8
            alpha = log(P_int/P_top)
            beta =  log(P_mid/P_int)/log(P_mid/P_top)
      
            a =  ( (1._r8 + alpha) * exp(-alpha) - 1._r8 ) / alpha
            b =  1._r8 - exp(-alpha)
   
            if (alpha .gt. 0) then             ! only apply where top level is below 80 km
               chi_mid = gas_vmr(i,1)          ! molar mixing ratio of O3 at midpoint of top layer
               chi_0 = chi_mid /  (1._r8 + beta)
               chi_eff = chi_0 * (a + b)
               gas_vmr(i,1) = chi_eff
               chi_eff = chi_eff * P_int / massratio / 9.8_r8 ! O3 column above in kg m-2
               chi_eff = chi_eff / 2.1415e-5_r8               ! O3 column above in DU
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

subroutine rrtmgp_set_cloud_sw( &
   nday, nlay, idxday, pmid, cldfrac,                       &
   c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, kdist_sw, cloud_sw)

   ! Create MCICA stochastic arrays for cloud SW optical properties.
   ! Initialize optical properties object (cloud_sw) and load with MCICA columns.
   !
   ! The input optical properties are on the CAM grid and are represented as products
   ! of the extinction optical depth (tau), single scattering albedo (w) and assymetry
   ! parameter (g).  This routine subsets the input to just the layers and the
   ! daylight columns used in the radiation calculation.  It also computes the
   ! individual properties of tau, w, and g for input to the MCICA routine.

   ! arguments
   integer,  intent(in) :: nday           ! number of daylight columns
   integer,  intent(in) :: nlay           ! number of layers in radiation calculation (may include "extra layer")
   integer,  intent(in) :: idxday(:)      ! indices of daylight columns in the chunk
   real(r8), intent(in) :: pmid(nday,nlay)! pressure at layer midpoints (Pa) used to seed RNG.

   ! cloud fraction and optics are input on the CAM grid
   real(r8), intent(in) :: cldfrac(pcols,pver)                ! combined cloud fraction
   real(r8), intent(in) :: c_cld_tau    (nswbands,pcols,pver) ! combined cloud extinction optical depth
   real(r8), intent(in) :: c_cld_tau_w  (nswbands,pcols,pver) ! combined cloud single scattering albedo * tau
   real(r8), intent(in) :: c_cld_tau_w_g(nswbands,pcols,pver) ! combined cloud assymetry parameter * w * tau

   class(ty_gas_optics_rrtmgp), intent(in)  :: kdist_sw  ! shortwave gas optics object
   type(ty_optical_props_2str), intent(out) :: cloud_sw  ! cloud optical properties object

   ! local vars
   integer, parameter :: changeseed = 1

   integer :: i, k, kk, ns, igpt
   integer :: ngptsw
   integer :: nver

   real(r8), allocatable :: cldf(:,:)
   real(r8), allocatable :: tauc(:,:,:)
   real(r8), allocatable :: ssac(:,:,:)
   real(r8), allocatable :: asmc(:,:,:)
   real(r8), allocatable :: taucmcl(:,:,:)
   real(r8), allocatable :: ssacmcl(:,:,:)
   real(r8), allocatable :: asmcmcl(:,:,:)

   real(r8) :: small_val = 1.e-80_r8
   real(r8), allocatable :: day_cld_tau(:,:,:)
   real(r8), allocatable :: day_cld_tau_w(:,:,:)
   real(r8), allocatable :: day_cld_tau_w_g(:,:,:)

   character(len=128) :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_cloud_sw'
   !--------------------------------------------------------------------------------

   ! number of g-points.  This is the number of subcolumns constructed by MCICA.
   ngptsw = kdist_sw%get_ngpt()

   ! number of CAM's layers in radiation calculation.  Does not include the "extra layer".
   nver   = pver - ktopcam + 1

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

   ! Subset the input data so just the daylight columns, and the number of CAM layers in the
   ! radiation calculation are used by MCICA to produce subcolumns.
   cldf            = cldfrac(         idxday(1:nday), ktopcam:)
   day_cld_tau     = c_cld_tau(    :, idxday(1:nday), ktopcam:)
   day_cld_tau_w   = c_cld_tau_w(  :, idxday(1:nday), ktopcam:)
   day_cld_tau_w_g = c_cld_tau_w_g(:, idxday(1:nday), ktopcam:)

   ! Compute the optical properties needed for the 2-stream calculations.  These calculations
   ! are the same as the RRTMG version.

   ! set cloud optical depth, clip @ zero
   tauc = merge(day_cld_tau, 0.0_r8, day_cld_tau > 0.0_r8)
   ! set value of asymmetry
   asmc = merge(day_cld_tau_w_g / max(day_cld_tau_w, small_val), 0.0_r8, day_cld_tau_w > 0.0_r8)
   ! set value of single scattering albedo
   ssac = merge(max(day_cld_tau_w, small_val) / max(tauc, small_val), 1.0_r8 , tauc > 0.0_r8)
   ! set asymmetry to zero when tauc = 0
   asmc = merge(asmc, 0.0_r8, tauc > 0.0_r8)

   ! MCICA converts from bands to gpts (e.g., 224 g-points instead of 14 bands)
   call mcica_subcol_sw( &
      kdist_sw, nswbands, ngptsw, nday, nlay, &
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
   do igpt = 1,ngptsw
      cloud_sw%g  (:, ktoprad:, igpt) = asmcmcl(igpt, ktopcam:, :)
      cloud_sw%ssa(:, ktoprad:, igpt) = ssacmcl(igpt, ktopcam:, :)
      cloud_sw%tau(:, ktoprad:, igpt) = taucmcl(igpt, ktopcam:, :)
   end do

   ! validate checks the tau > 0, ssa is in range [0,1], and g is in range [-1,1].
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

end subroutine rrtmgp_set_cloud_sw

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

subroutine rrtmgp_set_aer_sw( &
   nday, idxday, aer_tau, aer_tau_w, &
   aer_tau_w_g, aer_tau_w_f, aer_sw)

   ! Load aerosol SW optical properties into the RRTMGP object.
   !
   ! CAM fields are products tau, tau*ssa, tau*ssa*asy, tau*ssa*asy*fsf
   ! Fields expected by RRTMGP are computed by
   ! aer_sw%tau = aer_tau
   ! aer_sw%ssa = aer_tau_w / aer_tau
   ! aer_sw%g   = aer_tau_w_g / aer_taw_w
   !
   ! The input optical arrays from CAM are dimensioned in the vertical
   ! as 0:pver.  The index 0 is for the extra layer used in the radiation
   ! calculation.  The index ktopcam assumes the CAM vertical indices are
   ! in the range 1:pver, so using this index correctly ignores vertical
   ! index 0.  If an "extra" layer is used in the calculations, it is
   ! provided and set in the RRTMGP aerosol object aer_sw.

   ! Arguments
   integer,   intent(in) :: nday
   integer,   intent(in) :: idxday(:)
   real(r8),  intent(in) :: aer_tau    (pcols,0:pver,nswbands) ! extinction optical depth
   real(r8),  intent(in) :: aer_tau_w  (pcols,0:pver,nswbands) ! single scattering albedo * tau
   real(r8),  intent(in) :: aer_tau_w_g(pcols,0:pver,nswbands) ! asymmetry parameter * w * tau
   real(r8),  intent(in) :: aer_tau_w_f(pcols,0:pver,nswbands) ! forward scattered fraction * w * tau
   type(ty_optical_props_2str), intent(inout) :: aer_sw

   ! local variables
   integer  :: i

   ! minimum value for aer_tau_w is the same as used in RRTMG code.
   real(r8), parameter :: tiny = 1.e-80_r8

   character(len=32)  :: sub = 'rrtmgp_set_aer_sw'
   character(len=128) :: errmsg
   !--------------------------------------------------------------------------------

   ! If there is an extra layer in the radiation then this initialization
   ! will provide default values there.
   aer_sw%tau = 0.0_r8
   aer_sw%ssa = 1.0_r8
   aer_sw%g   = 0.0_r8

   do i = 1, nday
      ! aer_sw arrays have dimensions of (nday,nlay,nswbands)
      aer_sw%tau(i,ktoprad:,:) = max(aer_tau(idxday(i),ktopcam:,:), 0._r8)
      aer_sw%ssa(i,ktoprad:,:) = merge(aer_tau_w(idxday(i),ktopcam:,:)/aer_tau(idxday(i),ktopcam:,:), &
                                       1._r8, aer_tau(idxday(i),ktopcam:,:) > 0._r8)
      aer_sw%g(i,ktoprad:,:) = merge(aer_tau_w_g(idxday(i),ktopcam:,:)/aer_tau_w(idxday(i),ktopcam:,:), &
                                     0._r8, aer_tau_w(idxday(i),ktopcam:,:) > tiny)
   end do

   ! impose limits on the components:
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
