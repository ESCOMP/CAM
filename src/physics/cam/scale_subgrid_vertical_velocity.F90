! Apply min/max/scale to subgrid vertical velocity and derive wsubi.
module scale_subgrid_vertical_velocity

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: scale_subgrid_vertical_velocity_run

contains

!> \section arg_table_scale_subgrid_vertical_velocity_run Argument Table
!! \htmlinclude scale_subgrid_vertical_velocity_run.html
  subroutine scale_subgrid_vertical_velocity_run( &
    ncol, pver, top_lev,                          &
    wsub_min, wsubi_min,                          &
    wsub_scale, wsubi_scale,                      &
    use_preexisting_ice,                          &
    wsub, wsubi,                                  &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev             ! top vertical level for cloud physics

    real(kind_phys),  intent(in)    :: wsub_min            ! minimum subgrid vertical velocity for droplet nucleation [m s-1]
    real(kind_phys),  intent(in)    :: wsubi_min           ! minimum subgrid vertical velocity for ice nucleation [m s-1]
    real(kind_phys),  intent(in)    :: wsub_scale          ! scaling factor for subgrid vertical velocity (liquid)
    real(kind_phys),  intent(in)    :: wsubi_scale         ! scaling factor for subgrid vertical velocity (ice)

    logical,          intent(in)    :: use_preexisting_ice ! use preexisting ice? [flag]

    ! Input/Output arguments
    real(kind_phys),  intent(inout) :: wsub(:, :)          ! subgrid vertical velocity [m s-1] (ncol, pver)

    ! Output arguments
    real(kind_phys),  intent(out)   :: wsubi(:, :)         ! subgrid vertical velocity for ice nucleation [m s-1] (ncol, pver)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    ! Set minimum values above top_lev
    wsub (:, :top_lev-1) = wsub_min
    wsubi(:, :top_lev-1) = wsubi_min

    do k = top_lev, pver
      do i = 1, ncol
        ! Derive wsubi from raw wsub BEFORE scaling wsub
        wsubi(i,k) = max(wsubi_min, wsub(i,k)) * wsubi_scale
        if (.not. use_preexisting_ice) then
          wsubi(i,k) = min(wsubi(i,k), 0.2_kind_phys)
        end if

        wsub(i,k) = max(wsub_min, wsub(i,k)) * wsub_scale
      end do
    end do

  end subroutine scale_subgrid_vertical_velocity_run

end module scale_subgrid_vertical_velocity
