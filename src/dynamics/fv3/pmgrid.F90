module pmgrid

! Initialize grid point resolution parameters

implicit none

integer, parameter :: plon   = PLON       ! number of longitudes
integer, parameter :: plev   = PLEV       ! number of vertical levels
integer, parameter :: plat   = PLAT       ! number of latitudes

integer, parameter :: plevp  = plev+1     ! plev + 1
integer, parameter :: plnlv = plon*plev   ! Length of multilevel field slice

integer :: spmd_on = 0 ! 1 for Spmd, 0 for non-Spmd

end module pmgrid
