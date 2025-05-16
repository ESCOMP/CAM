!===============================================================================
! Dust for CARMA Aerosol Model
!===============================================================================
module dust_model
  use shr_kind_mod,    only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,      only: masterproc
  use cam_abortutils,  only: endrun

  implicit none
  private

  public :: dust_names
  public :: dust_nbin

  integer, parameter :: dust_nbin = 4

  character(len=6), parameter :: dust_names(dust_nbin) &
       = (/'NULL01', 'NULL02', 'NULL03', 'NULL04'/)

end module dust_model
