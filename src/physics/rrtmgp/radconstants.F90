module radconstants

! This module contains constants that are specific to the radiative transfer
! code used in the RRTMGP model.

use shr_kind_mod,         only: r8 => shr_kind_r8
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
use cam_abortutils,       only: endrun
use radiation_utils,      only: get_sw_spectral_boundaries_ccpp
use radiation_utils,      only: get_lw_spectral_boundaries_ccpp

implicit none
private
save

! Number of bands in SW and LW.  These values must match data in the RRTMGP coefficients datasets.
! But they are needed to allocate space in the physics buffer and need to be available before the
! RRTMGP datasets are read.  So they are set as parameters here and checked in the 
! set_wavenumber_bands subroutine after the datasets are read.
integer, parameter, public :: nswbands = 14
integer, parameter, public :: nlwbands = 16

! Band limits (set from data in RRTMGP coefficient datasets)
real(r8), target :: wavenumber_low_shortwave(nswbands)
real(r8), target :: wavenumber_high_shortwave(nswbands)
real(r8), target :: wavenumber_low_longwave(nlwbands)
real(r8), target :: wavenumber_high_longwave(nlwbands)

logical :: wavenumber_boundaries_set = .true.

! First and last g-point for each band.
integer, public, protected :: band2gpt_sw(2,nswbands)

! These are indices to specific bands for diagnostic output and COSP input.
integer, public, protected :: idx_sw_diag = -1     ! band contains 500-nm wave
integer, public, protected :: idx_nir_diag = -1    ! band contains 1000-nm wave
integer, public, protected :: idx_uv_diag = -1     ! band contains 400-nm wave
integer, public, protected :: idx_lw_diag = -1     ! band contains 1000 cm-1 wave (H20 window)

! GASES TREATED BY RADIATION (line spectra)
! These names are recognized by RRTMGP.  They are in the coefficients files as
! lower case strings.  These upper case names are used by CAM's namelist and 
! rad_constituents module.
integer, public, parameter :: gasnamelength = 5
integer, public, parameter :: nradgas = 8
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/'H2O  ','O3   ', 'O2   ', 'CO2  ', 'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)

! what is the minimum mass mixing ratio that can be supported by radiation implementation?
real(r8), public, parameter :: minmmr(nradgas) = epsilon(1._r8)

public :: &
   radconstants_init,          &
   get_sw_spectral_boundaries, &
   get_lw_spectral_boundaries, &
   get_band_index_by_value,    &
   rad_gas_index

!=========================================================================================
contains
!=========================================================================================
subroutine radconstants_init(idx_sw_diag_in, idx_nir_diag_in, idx_uv_diag_in, idx_lw_diag_in)
   integer, intent(in) :: idx_sw_diag_in
   integer, intent(in) :: idx_nir_diag_in
   integer, intent(in) :: idx_uv_diag_in
   integer, intent(in) :: idx_lw_diag_in

   idx_sw_diag = idx_sw_diag_in
   idx_nir_diag = idx_nir_diag_in
   idx_uv_diag = idx_uv_diag_in
   idx_lw_diag = idx_lw_diag_in

end subroutine radconstants_init
!=========================================================================================

 subroutine get_sw_spectral_boundaries(low_boundaries, high_boundaries, units)

   ! provide spectral boundaries of each shortwave band

   real(r8), dimension(:), intent(out) :: low_boundaries
   real(r8), dimension(:), intent(out) :: high_boundaries
   character(*), intent(in) :: units ! requested units

   character(len=512) :: errmsg
   integer :: errflg
   !----------------------------------------------------------------------------

   call get_sw_spectral_boundaries_ccpp(low_boundaries, high_boundaries, units, errmsg, errflg)
   if (errflg /= 0) then
      call endrun(errmsg)
   end if

 end subroutine get_sw_spectral_boundaries

!=========================================================================================

subroutine get_lw_spectral_boundaries(low_boundaries, high_boundaries, units)

   ! provide spectral boundaries of each longwave band

   real(r8), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
   character(*), intent(in) :: units ! requested units

   character(len=512) :: errmsg
   integer :: errflg
   !----------------------------------------------------------------------------

   call get_lw_spectral_boundaries_ccpp(low_boundaries, high_boundaries, units, errmsg, errflg)
   if (errflg /= 0) then
      call endrun(errmsg)
   end if

end subroutine get_lw_spectral_boundaries
  
!=========================================================================================

integer function rad_gas_index(gasname)

   ! return the index in the gaslist array of the specified gasname

   character(len=*),intent(in) :: gasname
   integer :: igas

   rad_gas_index = -1
   do igas = 1, nradgas
      if (trim(gaslist(igas)).eq.trim(gasname)) then
         rad_gas_index = igas
         return
      endif
   enddo
   call endrun ("rad_gas_index: can not find gas with name "//gasname)
end function rad_gas_index

!=========================================================================================

function get_band_index_by_value(swlw, targetvalue, units) result(ans)

   ! Find band index for requested wavelength/wavenumber.

   character(len=*), intent(in) :: swlw        ! sw or lw bands
   real(r8),         intent(in) :: targetvalue  
   character(len=*), intent(in) :: units       ! units of targetvalue
   integer :: ans

   ! local
   real(r8), pointer, dimension(:) :: lowboundaries, highboundaries
   real(r8) :: tgt
   integer  :: nbnds, i

   character(len=128) :: errmsg
   character(len=*), parameter :: sub = 'get_band_index_by_value'
   !----------------------------------------------------------------------------

   select case (swlw)
   case ('sw','SW','shortwave')
      nbnds = nswbands
      lowboundaries  => wavenumber_low_shortwave
      highboundaries => wavenumber_high_shortwave
   case ('lw', 'LW', 'longwave')
      nbnds = nlwbands
      lowboundaries  => wavenumber_low_longwave
      highboundaries => wavenumber_high_longwave
   case default
      call endrun('radconstants.F90: get_band_index_by_value: type of bands not recognized: '//swlw)
   end select

   ! band info is in cm^-1 but target value may be other units,
   ! so convert targetvalue to cm^-1
   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      tgt = targetvalue
   case('m','meter','meters')
      tgt = 1.0_r8 / (targetvalue * 1.e2_r8)
   case('nm','nanometer','nanometers')
      tgt = 1.0_r8 / (targetvalue * 1.e-7_r8)
   case('um','micrometer','micrometers','micron','microns')
      tgt = 1.0_r8 / (targetvalue * 1.e-4_r8)
   case('cm','centimeter','centimeters')
      tgt = 1._r8/targetvalue
   case default
      call endrun('radconstants.F90: get_band_index_by_value: units not recognized: '//units)
   end select

   ! now just loop through the array
   ans = 0
   do i = 1,nbnds
      if ((tgt > lowboundaries(i)) .and. (tgt <= highboundaries(i))) then
         ans = i
         exit
      end if
   end do

   if (ans == 0) then
      write(errmsg,'(f10.3,a,a)') targetvalue, ' ', trim(units)
      call endrun(sub//': band not found containing wave: '//trim(errmsg))
   end if
   
end function get_band_index_by_value

!=========================================================================================

end module radconstants
