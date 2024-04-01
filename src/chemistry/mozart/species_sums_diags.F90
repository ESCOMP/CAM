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
  use cam_logfile,      only : iulog
  use mo_chem_utls,     only : get_spc_ndx
  use shr_expr_parser_mod , only : shr_exp_parse, shr_exp_item_t, shr_exp_list_destroy

  implicit none
  private
  public :: species_sums_init
  public :: species_sums_output
  public :: species_sums_readnl
  public :: species_sums_final

  type(shr_exp_item_t), pointer :: vmr_grps => null()
  type(shr_exp_item_t), pointer :: mmr_grps => null()

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
    integer :: unitn, ierr, n,i

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

    if (masterproc) then

       write(iulog,*) ' '
       write(iulog,*) 'species_sums_readnl -- vmr_sums :'
       n = count(len_trim(vmr_sums)>0)
       do i=1,n
          write(iulog,*) trim(vmr_sums(i))
       end do

       write(iulog,*) ' '
       write(iulog,*) 'species_sums_readnl -- mmr_sums :'
       n = count(len_trim(mmr_sums)>0)
       do i=1,n
          write(iulog,*) trim(mmr_sums(i))
       end do

    end if


  end subroutine species_sums_readnl
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine species_sums_init

    type(shr_exp_item_t), pointer :: grp => null()

    integer :: j

    ! parse the terms of the summations
    vmr_grps => shr_exp_parse( vmr_sums )
    deallocate( vmr_sums )
    mmr_grps => shr_exp_parse( mmr_sums )
    deallocate( mmr_sums )

    ! add history fields

    if (masterproc) write(iulog,*) 'species_sums_init -- VMR SUMS:'
    grp => vmr_grps
    do while(associated(grp))

       if (masterproc) then
          write(iulog,*) ' grp name : ',trim(grp%name)

          do j = 1, grp%n_terms
             write(iulog,'(f12.4,a,a)') grp%coeffs(j),' * ',trim(grp%vars(j))
          end do
       end if

       call addfld( trim(grp%name), (/ 'lev' /),'A', 'mole/mole','summation of species volume mixing ratios')
       grp => grp%next_item
    enddo

    if (masterproc) write(iulog,*) 'species_sums_init -- MMR SUMS:'
    grp => mmr_grps
    do while(associated(grp))

       if (masterproc) then
          write(iulog,*) ' grp name : ',trim(grp%name)

          do j = 1, grp%n_terms
             write(iulog,'(f12.4,a,a)') grp%coeffs(j),' * ',trim(grp%vars(j))
          end do
       end if

       call addfld( trim(grp%name), (/ 'lev' /),'A', 'kg/kg','summation of species mass mixing ratios')
       grp => grp%next_item
    enddo

  end subroutine species_sums_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine species_sums_output( vmr, mmr, ncol, lchnk )

    real(r8), intent(in)    :: vmr(:,:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    integer,  intent(in)    :: ncol, lchnk

    integer :: j, spc_ndx
    real(r8) :: group_sum(ncol,pver)
    character(len=16) :: spc_name
    type(shr_exp_item_t), pointer :: grp

    ! output species groups ( or families )
    grp => vmr_grps
    do while(associated(grp))
       ! look up the corresponding species index ...
       group_sum(:,:) = 0._r8
       do j = 1, grp%n_terms
          spc_name = grp%vars(j)
          spc_ndx = get_spc_ndx( spc_name )
          if ( spc_ndx < 1 ) then
             call endrun('species_sums_output species name not found : '//trim(spc_name))
          endif
          group_sum(:ncol,:) = group_sum(:ncol,:) + grp%coeffs(j)*vmr(:ncol,:,spc_ndx)
       enddo
       call outfld( trim(grp%name), group_sum(:ncol,:), ncol, lchnk )
       grp => grp%next_item
    end do

    grp => mmr_grps
    do while(associated(grp))
       ! look up the corresponding species index ...
       group_sum(:,:) = 0._r8
       do j = 1, grp%n_terms
          spc_name = grp%vars(j)
          spc_ndx = get_spc_ndx( spc_name )
          if ( spc_ndx < 1 ) then
             call endrun('species_sums_output species name not found : '//trim(spc_name))
          endif
          group_sum(:ncol,:) = group_sum(:ncol,:) + grp%coeffs(j)*mmr(:ncol,:,spc_ndx)
       enddo
       call outfld( trim(grp%name), group_sum(:ncol,:), ncol, lchnk )
       grp => grp%next_item
    end do

  end subroutine species_sums_output

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine species_sums_final

    if (associated(vmr_grps)) call shr_exp_list_destroy(vmr_grps)
    if (associated(mmr_grps)) call shr_exp_list_destroy(mmr_grps)

  end subroutine species_sums_final

end module species_sums_diags
