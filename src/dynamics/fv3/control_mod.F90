! This module contains constants and namelist variables used through out the model
! to avoid circular dependancies please do not 'use' any further modules here.
!
module control_mod
  use shr_kind_mod,     only: r8=>shr_kind_r8

  integer, public, parameter :: MAX_STRING_LEN=240
  integer, public, parameter :: MAX_FILE_LEN=240
  integer              , public :: runtype
  logical, public :: disable_diagnostics = .FALSE.

end module control_mod
