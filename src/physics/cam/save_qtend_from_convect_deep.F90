! Save water vapor tendency from deep convection
! for use by gravity wave parameterization
!
! This scheme has to be run after all deep convective schemes,
! before the final tendencies are applied.
! Since ZM applies tendencies in each sub-scheme, this scheme has
! to run repeatedly to archive all tendencies from ZM.
module save_qtend_from_convect_deep
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: save_qtend_from_convect_deep_timestep_init
  public :: save_qtend_from_convect_deep_run

contains

!> \section arg_table_save_qtend_from_convect_deep_timestep_init Argument Table
!! \htmlinclude save_qtend_from_convect_deep_timestep_init.html
  subroutine save_qtend_from_convect_deep_timestep_init(ncol, pver, qtend_dp, errmsg, errflg)
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver

    real(kind_phys),    intent(out) :: qtend_dp(:, :)          ! Water vapor tendency from deep convection [s-1]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    qtend_dp(:,:) = 0._kind_phys

  end subroutine save_qtend_from_convect_deep_timestep_init

!> \section arg_table_save_qtend_from_convect_deep_run Argument Table
!! \htmlinclude save_qtend_from_convect_deep_run.html
  subroutine save_qtend_from_convect_deep_run( &
    ncol, pver, &
    qtend,  &
    qtend_dp, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: qtend(:, :)              ! Water vapor tendency from deep convection [s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: qtend_dp(:, :)         ! Water vapor tendency from deep convection [s-1] stored for MCSP
    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    ! Convert to temperature tendency
    do k = 1, pver
      do i = 1, ncol
        qtend_dp(i, k) = qtend_dp(i, k) + qtend(i, k)
      end do
    end do

  end subroutine save_qtend_from_convect_deep_run

end module save_qtend_from_convect_deep
