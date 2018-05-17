module pmgrid

! Initialize grid point resolution parameters

implicit none

integer, parameter :: plon   = PLON       ! number of longitudes
integer, parameter :: plev   = PLEV       ! number of vertical levels
integer, parameter :: plat   = PLAT       ! number of latitudes

integer, parameter :: plevp  = plev+1     ! plev + 1
integer, parameter :: plnlv = plon*plev   ! Length of multilevel field slice

integer :: spmd_on = 0 ! 1 for Spmd, 0 for non-Spmd

! Staggered grid parameters

integer, parameter :: splon = plon     ! Number of longitudes on the staggered grid
integer, parameter :: splat = plat     ! Number of latitudes on the staggered grid

! Note: In reality, the staggered latitude array for Lin-Rood dynamics only
! uses PLAT-1 latitudes, the first one being ignored.  So ideally the line
! above should read:
!   parameter (splat = plat-1)
! to define the staggered latitude winds with the correct dimension.
! However, the assumption that the staggered latitude grid has one extra
! latitude (making it the same dimension as the non-staggered grid) is
! pervasive throughout the Lin-Rood dynamical core, necessitating the
! extra latitude.

end module pmgrid
