module zonal_mean

use shr_kind_mod, only: r8 => shr_kind_r8
use dynamics_vars, only: T_FVDYCORE_GRID
use pmgrid, only: plon

implicit none
private
save

public :: zonal_mean_3D

real(r8), parameter :: rplon = 1._r8/plon

! External that does parallel sums reproducibly.
interface
   subroutine par_xsum(grid, a, ltot, sum)
     import
     type (T_FVDYCORE_GRID), intent(in) :: grid
     integer, intent(in) :: ltot
     real (r8), intent(in) :: a(grid%ifirstxy:grid%ilastxy,ltot)
     real (r8) sum(ltot)
   end subroutine par_xsum
end interface

contains

subroutine zonal_mean_3D(grid, nlev, fld_orig, fld_zm)

  ! FV dynamics grid
  type(T_FVDYCORE_GRID), intent(in) :: grid
  ! Number of vertical levels
  integer, intent(in) :: nlev
  ! Original field
  real(r8), intent(in)  :: fld_orig(grid%ifirstxy:grid%ilastxy,nlev,grid%jfirstxy:grid%jlastxy)
  ! Zonal mean field
  real(r8), intent(out) :: fld_zm(nlev,grid%jfirstxy:grid%jlastxy)

  integer :: j

  ! Rename grid bounds for convenience.
  associate(beglon => grid%ifirstxy, &
    endlon => grid%ilastxy, &
    beglat => grid%jfirstxy, &
    endlat => grid%jlastxy)

    do j = beglat, endlat
       call par_xsum( grid, fld_orig(beglon:endlon,:,j), nlev, fld_zm(:,j) )
       fld_zm(:,j) = fld_zm(:,j) * rplon
    end do

  end associate

end subroutine zonal_mean_3D

end module zonal_mean
