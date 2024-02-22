module radconstants

! This module contains constants that are specific to the radiative transfer
! code used in the RRTMGP model.

use shr_kind_mod,         only: r8 => shr_kind_r8
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
use cam_abortutils,       only: endrun

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

logical :: wavenumber_boundaries_set = .false.

integer, public, protected :: nswgpts  ! number of SW g-points
integer, public, protected :: nlwgpts  ! number of LW g-points

! These are indices to specific bands for diagnostic output and COSP input.
integer, public, protected :: idx_sw_diag = -1     ! band contains 500-nm wave
integer, public, protected :: idx_nir_diag = -1    ! band contains 1000-nm wave
integer, public, protected :: idx_uv_diag = -1     ! band contains 400-nm wave
integer, public, protected :: idx_lw_diag = -1     ! band contains 1000 cm-1 wave (H20 window)
integer, public, protected :: idx_sw_cloudsim = -1 ! band contains 670-nm wave (for COSP)
integer, public, protected :: idx_lw_cloudsim = -1 ! band contains 10.5 micron wave (for COSP)

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
   set_wavenumber_bands,       &
   get_sw_spectral_boundaries, &
   get_lw_spectral_boundaries, &
   get_band_index_by_value,    &
   rad_gas_index

!=========================================================================================
contains
!=========================================================================================

subroutine set_wavenumber_bands(kdist_sw, kdist_lw)

   ! Set the low and high limits of the wavenumber grid for sw and lw.
   ! Values come from RRTMGP coefficients datasets, and are stored in the
   ! kdist objects.
   !
   ! Set band indices for bands containing specific wavelengths.

   ! Arguments
   type(ty_gas_optics_rrtmgp), intent(in) :: kdist_sw
   type(ty_gas_optics_rrtmgp), intent(in) :: kdist_lw

   ! Local variables
   integer :: istat
   real(r8), allocatable :: values(:,:)

   character(len=128) :: errmsg
   character(len=*), parameter :: sub = 'set_wavenumber_bands'
   !----------------------------------------------------------------------------

   ! Check that number of sw/lw bands in gas optics files matches the parameters.
   if (kdist_sw%get_nband() /= nswbands) then
      write(errmsg,'(a,i4,a,i4)') 'number of sw bands in file, ', kdist_sw%get_nband(), &
         ", doesn't match parameter nswbands= ", nswbands
      call endrun(sub//': ERROR: '//trim(errmsg))
   end if
   if (kdist_lw%get_nband() /= nlwbands) then
      write(errmsg,'(a,i4,a,i4)') 'number of lw bands in file, ', kdist_lw%get_nband(), &
         ", doesn't match parameter nlwbands= ", nlwbands
      call endrun(sub//': ERROR: '//trim(errmsg))
   end if

   nswgpts = kdist_sw%get_ngpt()
   nlwgpts = kdist_lw%get_ngpt()

   ! SW band bounds in cm^-1
   allocate( values(2,nswbands), stat=istat )
   if (istat/=0) then
      call endrun(sub//': ERROR allocating array: values(2,nswbands)')
   end if
   values = kdist_sw%get_band_lims_wavenumber()
   wavenumber_low_shortwave = values(1,:)
   wavenumber_high_shortwave = values(2,:)

   ! Indices into specific bands
   idx_sw_diag     = get_band_index_by_value('sw', 500.0_r8, 'nm')
   idx_nir_diag    = get_band_index_by_value('sw', 1000.0_r8, 'nm')
   idx_uv_diag     = get_band_index_by_value('sw', 400._r8, 'nm')
   idx_sw_cloudsim = get_band_index_by_value('sw', 0.67_r8, 'micron')

   deallocate(values)

   ! LW band bounds in cm^-1
   allocate( values(2,nlwbands), stat=istat )
   if (istat/=0) then
      call endrun(sub//': ERROR allocating array: values(2,nlwbands)')
   end if
   values = kdist_lw%get_band_lims_wavenumber()
   wavenumber_low_longwave = values(1,:)
   wavenumber_high_longwave = values(2,:)

   ! Indices into specific bands
   idx_lw_diag     = get_band_index_by_value('lw', 1000.0_r8, 'cm^-1')
   idx_lw_cloudsim = get_band_index_by_value('lw', 10.5_r8, 'micron')

   wavenumber_boundaries_set = .true.

end subroutine set_wavenumber_bands

!=========================================================================================

subroutine get_sw_spectral_boundaries(low_boundaries, high_boundaries, units)

   ! provide spectral boundaries of each shortwave band

   real(r8),    intent(out) :: low_boundaries(nswbands), high_boundaries(nswbands)
   character(*), intent(in) :: units ! requested units

   character(len=*), parameter :: sub = 'get_sw_spectral_boundaries'
   !----------------------------------------------------------------------------

   if (.not. wavenumber_boundaries_set) then
      call endrun(sub//': ERROR, wavenumber boundaries not set. ')
   end if

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries = wavenumber_low_shortwave
      high_boundaries = wavenumber_high_shortwave
   case('m','meter','meters')
      low_boundaries = 1.e-2_r8/wavenumber_high_shortwave
      high_boundaries = 1.e-2_r8/wavenumber_low_shortwave
   case('nm','nanometer','nanometers')
      low_boundaries = 1.e7_r8/wavenumber_high_shortwave
      high_boundaries = 1.e7_r8/wavenumber_low_shortwave
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries = 1.e4_r8/wavenumber_high_shortwave
      high_boundaries = 1.e4_r8/wavenumber_low_shortwave
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._r8/wavenumber_high_shortwave
      high_boundaries = 1._r8/wavenumber_low_shortwave
   case default
      call endrun(sub//': ERROR, requested spectral units not recognized: '//units)
   end select

end subroutine get_sw_spectral_boundaries

!=========================================================================================

subroutine get_lw_spectral_boundaries(low_boundaries, high_boundaries, units)

   ! provide spectral boundaries of each longwave band

   real(r8), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
   character(*), intent(in) :: units ! requested units

   character(len=*), parameter :: sub = 'get_lw_spectral_boundaries'
   !----------------------------------------------------------------------------

   if (.not. wavenumber_boundaries_set) then
      call endrun(sub//': ERROR, wavenumber boundaries not set. ')
   end if

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries  = wavenumber_low_longwave
      high_boundaries = wavenumber_high_longwave
   case('m','meter','meters')
      low_boundaries  = 1.e-2_r8/wavenumber_high_longwave
      high_boundaries = 1.e-2_r8/wavenumber_low_longwave
   case('nm','nanometer','nanometers')
      low_boundaries  = 1.e7_r8/wavenumber_high_longwave
      high_boundaries = 1.e7_r8/wavenumber_low_longwave
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries  = 1.e4_r8/wavenumber_high_longwave
      high_boundaries = 1.e4_r8/wavenumber_low_longwave
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._r8/wavenumber_high_longwave
      high_boundaries = 1._r8/wavenumber_low_longwave
   case default
      call endrun(sub//': ERROR, requested spectral units not recognized: '//units)
   end select

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
