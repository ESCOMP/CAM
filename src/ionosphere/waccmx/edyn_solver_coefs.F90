module edyn_solver_coefs

  use shr_kind_mod ,only: r8 => shr_kind_r8 ! 8-byte reals

!
! nc(1:9) are pointers to beginning of coefficient blocks at each of
!   levels of resolution:
! nc(1) = nc0, pointer to coefficients for highest resolution.
! nc(2) = nc1, pointer to coefficients at half the resolution of nc0,
!   and so on for nc(3), nc(4), nc(5), etc.
! nc(9) = ncee, the dimension of the entire cee array, containing
!   coefficients for all levels.
!
  integer :: nc(9)

!
! Coefficients are stored in 1-d array cee(ncee)
! cee transmits descretized dynamo PDE coefficients to the multi-grid
!   mudpack solver. (cee was formerly in ceee.h)
!
  real(r8), target, allocatable :: cee(:)
!
! Unmodified coefficients for using modified mudpack:
  real(r8), allocatable :: cofum(:,:,:)

contains

!-----------------------------------------------------------------------
  subroutine ceee(cee,nx,ny,cf)

!
! Called from mudpack solvers to transfer coefficients.
!
! Args:
    integer,intent(in) :: nx,ny
    real(r8),intent(in) :: cee(nx,ny,*)
    real(r8),intent(out) :: cf(nx,ny,*)
!
! Local:
    integer :: i,j,n

    do n = 1,9
      do j = 1,ny
        do i = 1,nx
          cf(i,j,n) = cee(i,j,n)
        enddo
      enddo
    enddo
  end subroutine ceee

end module edyn_solver_coefs
