module radconstants

! provide stubs to allow building with no radiation scheme active

use shr_kind_mod,   only: r8 => shr_kind_r8
use cam_abortutils, only: endrun

implicit none
private
save

integer, parameter, public :: nswbands = 1
integer, parameter, public :: nlwbands = 1
integer, parameter, public :: idx_sw_diag = 1
integer, parameter, public :: idx_lw_diag = 1
integer, parameter, public :: idx_nir_diag = 1
integer, parameter, public :: idx_uv_diag = 1
integer, parameter, public :: nrh = 1
integer, parameter, public :: ot_length = 32

public :: rad_gas_index

integer, public, parameter :: gasnamelength = 1
integer, public, parameter :: nradgas = 1
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/' '/)

!========================================================================================
contains
!========================================================================================

integer function rad_gas_index(gasname)

   character(len=*),intent(in) :: gasname

   call endrun('rad_gas_index: ERROR: this is a stub')

end function rad_gas_index

end module radconstants
