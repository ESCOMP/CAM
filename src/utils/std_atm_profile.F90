module std_atm_profile

!-------------------------------------------------------------------------------
!
! The barometric formula for U.S. Standard Atmosphere is valid up to 86 km.
! see https://en.wikipedia.org/wiki/Barometric_formula.
!
! N.B.  The extension above 86 km is using data from Hanli.  It is not complete
!       since the hardcoded parameter (c1) needs adjustment above 86 km.
!
!-------------------------------------------------------------------------------

use shr_kind_mod,        only: r8 => shr_kind_r8
use cam_logfile,         only: iulog
use cam_abortutils,      only: endrun

implicit none
private
save

public :: &
   std_atm_pres,   & ! compute pressure given height
   std_atm_height, & ! compute height given pressure
   std_atm_temp      ! compute temperature given height

! Parameters for barometric formula for U.S. Standard Atmosphere.

integer, parameter  :: nreg = 15  ! number of regions

real(r8), parameter :: hb(nreg) = & ! height at bottom of layer (m)
     (/0.0_r8, 1.1e4_r8, 2.0e4_r8, 3.2e4_r8, 4.7e4_r8, 5.1e4_r8, 7.1e4_r8, 8.6e4_r8, &
     9.1e4_r8, 1.1e5_r8, 1.2e5_r8, 1.5e5_r8, 2.0e5_r8, 3.0e5_r8, 7.e5_r8/)

real(r8), parameter :: pb(nreg) = & ! standard pressure (Pa)
     (/101325._r8, 22632.1_r8, 5474.89_r8, 868.02_r8, 110.91_r8, 66.94_r8, 3.96_r8, 3.7e-1_r8,  &
     1.5e-1_r8, 7.1e-3_r8, 2.5e-3_r8, 4.5e-4_r8, 8.47e-5_r8, 8.77e-6_r8, 3.19e-8_r8/)

real(r8), parameter :: tb(nreg) = & ! standard temperature (K)
     (/288.15_r8, 216.65_r8, 216.65_r8, 228.65_r8, 270.65_r8, 270.65_r8, 214.65_r8, 186.87_r8,  &
     186.87_r8, 240._r8, 360._r8, 634.39_r8, 854.56_r8, 976.01_r8, 1.e3_r8/)

real(r8), parameter :: lb(nreg) = & ! temperature lapse rate (K/m)
     (/-0.0065_r8, 0.0_r8, 0.001_r8, 0.0028_r8, 0.0_r8, -0.0028_r8, -0.001852_r8, 0.0_r8,       &
     2.796e-3_r8, 0.012_r8, 9.15e-3_r8, 4.4e-3_r8, 1.21e-3_r8, 6.e-5_r8, 0.0_r8/)

real(r8), parameter :: rg = 8.3144598_r8 ! universal gas constant (J/mol/K)
real(r8), parameter :: g0 = 9.80665_r8   ! gravitational acceleration (m/s^2)
real(r8), parameter :: mw = 0.0289644_r8 ! molar mass of dry air (kg/mol)
real(r8), parameter :: c1 = g0*mw/rg
  
!=========================================================================================
CONTAINS
!=========================================================================================

subroutine std_atm_pres(height, pstd)
    
   ! arguments
   real(r8), intent(in)  :: height(:) ! height above sea level in meters
   real(r8), intent(out) :: pstd(:)   ! std pressure in Pa
    
   integer :: i, ii, k, nlev
   character(len=*), parameter :: routine = 'std_atm_pres'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      if (height(k) < 0.0_r8) then
         ! Extrapolate below mean sea level using troposphere lapse rate.
         ii = 1
      else
         ! find region containing height
         find_region: do i = nreg, 1, -1
            if (height(k) >= hb(i)) then
               ii = i
               exit find_region
            end if
         end do find_region
      end if
      
      if (lb(ii) /= 0._r8) then
         pstd(k) = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(height(k) - hb(ii)) ) )**(c1/lb(ii))
      else
         pstd(k) = pb(ii) * exp( -c1*(height(k) - hb(ii))/tb(ii) )
      end if
      
   end do

end subroutine std_atm_pres

!=========================================================================================

subroutine std_atm_height(pstd, height)
    
   ! arguments
   real(r8), intent(in)   :: pstd(:)   ! std pressure in Pa
   real(r8), intent(out)  :: height(:) ! height above sea level in meters
    
   integer :: i, ii, k, nlev
   logical :: found_region
   character(len=*), parameter :: routine = 'std_atm_height'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      
      if (pstd(k) <= pb(nreg)) then
         ii = nreg
      else if (pstd(k) > pb(1)) then
         ii = 1
      else
         ! find region containing pressure
         find_region: do i = 2, nreg
            if (pstd(k) > pb(i)) then
               ii = i - 1
               exit find_region
            end if
         end do find_region
      end if

      if (lb(ii) /= 0._r8) then
         height(k) = hb(ii) + (tb(ii)/lb(ii)) * ( (pb(ii)/pstd(k))**(lb(ii)/c1) - 1._r8 )
      else
         height(k) = hb(ii) + (tb(ii)/c1)*log(pb(ii)/pstd(k))
      end if
   end do

end subroutine std_atm_height

!=========================================================================================

subroutine std_atm_temp(height, temp)
    
   ! arguments
   real(r8), intent(in)   :: height(:) ! std pressure in Pa
   real(r8), intent(out)  :: temp(:)   ! temperature
    
   ! local vars
   integer :: i, ii, k, nlev
   character(len=*), parameter :: routine = 'std_atm_temp'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      if (height(k) < 0.0_r8) then
         ii = 1
      else
         ! find region containing height
         find_region: do i = nreg, 1, -1
            if (height(k) >= hb(i)) then
               ii = i
               exit find_region
            end if
         end do find_region
      end if

      if (lb(ii) /= 0._r8) then
         temp(k) = tb(ii) + lb(ii)*(height(k) - hb(ii))
      else
         temp(k) = tb(ii)
      end if
      
   end do

end subroutine std_atm_temp

end module std_atm_profile
