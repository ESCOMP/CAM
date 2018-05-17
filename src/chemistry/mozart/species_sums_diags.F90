!--------------------------------------------------------------------------------
! Species summations for history
!--------------------------------------------------------------------------------
module species_sums_diags

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_kind_mod,     only : CL => SHR_KIND_CL
  use cam_history,      only : addfld
  use cam_history,      only : outfld
  use ppgrid,           only : pver
  use spmd_utils,       only : masterproc
  use cam_abortutils,   only : endrun
  use mo_chem_utls,     only : get_spc_ndx
  use sums_utils,       only : sums_grp_t, parse_sums

  implicit none
  private 
  public :: species_sums_init
  public :: species_sums_output
  public :: species_sums_readnl

  integer :: n_vmr_grps = 0
  type(sums_grp_t), allocatable :: vmr_grps(:)  
  integer :: n_mmr_grps = 0
  type(sums_grp_t), allocatable :: mmr_grps(:)  

  integer, parameter :: maxlines = 200
  character(len=CL), allocatable :: vmr_sums(:)
  character(len=CL), allocatable :: mmr_sums(:)

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine species_sums_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils,      only: mpicom, mpi_character, masterprocid

    ! args 
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    namelist /species_sums_nl/ vmr_sums, mmr_sums

    allocate( vmr_sums( maxlines ) )
    vmr_sums(:) = ' '
    allocate( mmr_sums( maxlines ) )
    mmr_sums(:) = ' '

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'species_sums_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, species_sums_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('species_sums_readnl:: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(vmr_sums,len(vmr_sums(1))*maxlines, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(mmr_sums,len(mmr_sums(1))*maxlines, mpi_character, masterprocid, mpicom, ierr)

  end subroutine species_sums_readnl
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine species_sums_init

    integer :: i

    ! parse the terms of the summations
    call parse_sums(vmr_sums, n_vmr_grps, vmr_grps)
    deallocate( vmr_sums )
    call parse_sums(mmr_sums, n_mmr_grps, mmr_grps)
    deallocate( mmr_sums )

    ! add history fields
    do i = 1, n_vmr_grps
       call addfld( vmr_grps(i)%name, (/ 'lev' /),'A', 'mole/mole','summation of species volume mixing ratios')
    enddo
    do i = 1, n_mmr_grps
       call addfld( mmr_grps(i)%name, (/ 'lev' /),'A', 'kg/kg','summation of species mass mixing ratios')
    enddo

  end subroutine species_sums_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine species_sums_output( vmr, mmr, ncol, lchnk )

    real(r8), intent(in)    :: vmr(:,:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    integer,  intent(in)    :: ncol, lchnk

    integer :: i, j, spc_ndx
    real(r8) :: group_sum(ncol,pver)
    character(len=16) :: spc_name
    
    ! output species groups ( or families )
    do i = 1, n_vmr_grps
       ! look up the corresponding species index ...
       group_sum(:,:) = 0._r8
       do j = 1, vmr_grps(i)%nmembers
          spc_name = vmr_grps(i)%term(j)
          spc_ndx = get_spc_ndx( spc_name )
          if ( spc_ndx < 1 ) then
             call endrun('species_sums_output species name not found : '//trim(spc_name))
          endif
          group_sum(:ncol,:) = group_sum(:ncol,:) + vmr_grps(i)%multipler(j)*vmr(:ncol,:,spc_ndx)
       enddo
       call outfld( vmr_grps(i)%name, group_sum(:ncol,:), ncol, lchnk )       
    end do
    do i = 1, n_mmr_grps
       ! look up the corresponding species index ...
       group_sum(:,:) = 0._r8
       do j = 1, mmr_grps(i)%nmembers
          spc_name = mmr_grps(i)%term(j)
          spc_ndx = get_spc_ndx( spc_name )
          if ( spc_ndx < 1 ) then
             call endrun('species_sums_output species name not found : '//trim(spc_name))
          endif
          group_sum(:ncol,:) = group_sum(:ncol,:) + mmr_grps(i)%multipler(j)*mmr(:ncol,:,spc_ndx)
       enddo
       call outfld( mmr_grps(i)%name, group_sum(:ncol,:), ncol, lchnk )       
    end do

  end subroutine species_sums_output

end module species_sums_diags
