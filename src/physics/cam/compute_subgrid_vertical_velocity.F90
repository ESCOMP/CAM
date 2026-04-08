! Compute raw subgrid-scale vertical velocity from TKE or KVH.
module compute_subgrid_vertical_velocity

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: compute_subgrid_vertical_velocity_tke_run
  public :: compute_subgrid_vertical_velocity_kvh_run

contains

!> \section arg_table_compute_subgrid_vertical_velocity_tke_run Argument Table
!! \htmlinclude compute_subgrid_vertical_velocity_tke_run.html
  subroutine compute_subgrid_vertical_velocity_tke_run( &
    ncol, pver, top_lev, &
    tke,                 &
    wsub,                &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)  :: ncol    ! number of columns
    integer,          intent(in)  :: pver    ! number of vertical levels
    integer,          intent(in)  :: top_lev ! top vertical level for cloud physics

    real(kind_phys),  intent(in)  :: tke(:, :)  ! turbulent kinetic energy at interfaces [m2/s2] (ncol, pver+1)

    ! Output arguments
    real(kind_phys),  intent(out) :: wsub(:, :) ! subgrid vertical velocity [m/s] (ncol, pver)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    wsub(:, :top_lev-1) = 0._kind_phys

    do k = top_lev, pver
      do i = 1, ncol
        wsub(i,k) = sqrt(0.5_kind_phys * (tke(i,k) + tke(i,k+1)) * (2._kind_phys / 3._kind_phys))
        wsub(i,k) = min(wsub(i,k), 10._kind_phys)
      end do
    end do

  end subroutine compute_subgrid_vertical_velocity_tke_run

!> \section arg_table_compute_subgrid_vertical_velocity_kvh_run Argument Table
!! \htmlinclude compute_subgrid_vertical_velocity_kvh_run.html
  subroutine compute_subgrid_vertical_velocity_kvh_run( &
    ncol, pver, top_lev, &
    kvh,                 &
    wsub,                &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)  :: ncol    ! number of columns
    integer,          intent(in)  :: pver    ! number of vertical levels
    integer,          intent(in)  :: top_lev ! top vertical level for cloud physics

    real(kind_phys),  intent(in)  :: kvh(:, :)  ! eddy diffusivity for heat at interfaces [m2/s] (ncol, pver+1)

    ! Output arguments
    real(kind_phys),  intent(out) :: wsub(:, :) ! subgrid vertical velocity [m/s] (ncol, pver)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    wsub(:, :top_lev-1) = 0._kind_phys

    ! Get sub-grid vertical velocity from diffusion coefficient.
    ! Following Morrison et al. 2005, JAS.
    ! Assume mixing length of 30 m.
    ! Use maximum sub-grid vertical vel of 10 m/s.
    do k = top_lev, pver
      do i = 1, ncol
        wsub(i,k) = (kvh(i,k) + kvh(i,k+1)) / 2._kind_phys / 30._kind_phys
        wsub(i,k) = min(wsub(i,k), 10._kind_phys)
      end do
    end do

  end subroutine compute_subgrid_vertical_velocity_kvh_run

end module compute_subgrid_vertical_velocity
