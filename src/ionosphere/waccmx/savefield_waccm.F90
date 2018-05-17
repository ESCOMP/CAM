module savefield_waccm
  use shr_kind_mod,only: r8 => shr_kind_r8 ! 8-byte reals
  use cam_history ,only: hist_fld_active,outfld  ! Routine to output fields to history files
  use edyn_mpi    ,only: array_ptr_type
!
! Save fields to WACCM output history file.
!
  implicit none
  save
  private
  public savefld_waccm,savefld_waccm_switch
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
            diag_ik(i,k) = f(k,i,j)
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
!-----------------------------------------------------------------------
  subroutine savefld_waccm_switch(f,name,plev,i0,i1,j0,j1)
!
! Copy input array to a local array, associate a pointer to the local array,
! switch the "model format" of the pointer (shift longitude and invert vertical),
! (TIEGCM to WACCM in this case), and save the local array to WACCM history.
! (Input array is unchanged)
!
    use edyn_mpi     ,only: switch_model_format
!
! Args:
    integer,intent(in)  :: plev,i0,i1,j0,j1
    real(r8),intent(in) :: f(plev,i0:i1,j0:j1)
    character(len=*),intent(in) :: name
!
! Local:
    integer :: i,j
    real(r8),target :: ftmp(plev,i0:i1,j0:j1)
    type(array_ptr_type) :: ptr(1)

    if (.not.hist_fld_active(name)) return

!
! Copy input to local array:
    do j=j0,j1
      do i=i0,i1
        ftmp(:,i,j) = f(:,i,j)
      enddo
    enddo
!
! Associate local pointer (lonshift_blocks expects an array_ptr_type)
    ptr(1)%ptr => ftmp
!
! Switch from TIEGCM format to WACCM format:
!
    call switch_model_format(ptr,1,plev,i0,i1,j0,j1,1)
!
! Return data to local array, and save on WACCM history:
!
    do j=j0,j1
      do i=i0,i1
        ftmp(1:plev,i,j) = ptr(1)%ptr(1:plev,i,j)
      enddo
    enddo

    call savefld_waccm(ftmp(:,i0:i1,j0:j1),trim(name),plev,i0,i1,j0,j1)

  end subroutine savefld_waccm_switch
!-----------------------------------------------------------------------
end module savefield_waccm
