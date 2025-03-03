!===============================================================================
! Seasalt for CARMA Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,       only: pcols, pver

  implicit none
  private

  public :: seasalt_names
  public :: seasalt_nbin

  integer, parameter :: seasalt_nbin = 4

  character(len=6), parameter :: seasalt_names(seasalt_nbin) &
       = (/'NULL01', 'NULL02', 'NULL03', 'NULL04'/)

end module seasalt_model
