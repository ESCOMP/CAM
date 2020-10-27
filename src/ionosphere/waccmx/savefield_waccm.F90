module savefield_waccm
  use shr_kind_mod,only: r8 => shr_kind_r8 ! 8-byte reals
  use cam_history ,only: hist_fld_active,outfld  ! Routine to output fields to history files
!
! Save fields to WACCM output history file.
!
  implicit none
  private
  public :: savefld_waccm

contains
!-----------------------------------------------------------------------
  subroutine savefld_waccm(f,name,nlev,i0,i1,j0,j1)
!
! Save field to WACCM history. 
! Call to addfld must be made in edyn_init for each field to be saved.
! Field names must be in user_nl_cam to be written to the file. 
!
! Args:
    integer,intent(in) :: nlev,i0,i1,j0,j1
    real(r8),dimension(nlev,i0:i1,j0:j1),intent(in) :: f
    character(len=*),intent(in) :: name
!
! Local:
    integer :: i,j,k
    real(r8) :: diag_ik(i0:i1,nlev)
!
    if (.not.hist_fld_active(name)) return

    if (nlev /= 1) then
      do j=j0,j1
        do i=i0,i1
          do k=1,nlev
            diag_ik(i,k) = f(nlev-k+1,i,j)
          enddo
        enddo
        call outfld(name,diag_ik,i1-i0+1,j)
      enddo
    else
      do j=j0,j1
        do i=i0,i1
          diag_ik(i,1) = f(1,i,j)
        enddo
        call outfld(name,diag_ik,i1-i0+1,j)
      enddo
    endif
  end subroutine savefld_waccm

end module savefield_waccm
