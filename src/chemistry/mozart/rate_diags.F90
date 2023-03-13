!--------------------------------------------------------------------------------
! Manages writing reaction rates to history
!--------------------------------------------------------------------------------
module rate_diags

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_kind_mod,     only : CL => SHR_KIND_CL
  use cam_history,      only : fieldname_len
  use cam_history,      only : addfld, add_default
  use cam_history,      only : outfld
  use chem_mods,        only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
  use ppgrid,           only : pver
  use spmd_utils,       only : masterproc
  use cam_logfile,      only : iulog
  use cam_abortutils,   only : endrun
  use shr_expr_parser_mod , only : shr_exp_parse, shr_exp_item_t, shr_exp_list_destroy

  implicit none
  private
  public :: rate_diags_init
  public :: rate_diags_calc
  public :: rate_diags_readnl
  public :: rate_diags_o3s_loss
  public :: rate_diags_final

  character(len=fieldname_len) :: rate_names(rxt_tag_cnt)

  type(shr_exp_item_t), pointer :: grps_list => null()

  integer, parameter :: maxlines = 200
  character(len=CL), allocatable :: rxn_rate_sums(:)

  integer :: o3_ndx = -1

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine rate_diags_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils,      only: mpicom, mpi_character, masterprocid

    ! args
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    namelist /rxn_rate_diags_nl/ rxn_rate_sums

    allocate( rxn_rate_sums( maxlines ) )
    rxn_rate_sums(:) = ' '

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'rxn_rate_diags_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, rxn_rate_diags_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('rate_diags_readnl:: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(rxn_rate_sums,len(rxn_rate_sums(1))*maxlines, mpi_character, masterprocid, mpicom, ierr)

  end subroutine rate_diags_readnl
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_init
    use phys_control, only : phys_getopts
    use mo_chem_utls, only : get_spc_ndx

    integer :: i,j, len, pos
    character(len=64) :: name
    logical :: history_scwaccm_forcing
    type(shr_exp_item_t), pointer :: grp

    call phys_getopts( history_scwaccm_forcing_out = history_scwaccm_forcing )

    do i = 1,rxt_tag_cnt
       pos = 0
       pos = index(rxt_tag_lst(i),'tag_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'usr_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'cph_')
       if (pos <= 0) pos = index(rxt_tag_lst(i),'ion_')
       if (pos>0) then
          name = 'r_'//trim(rxt_tag_lst(i)(5:))
       else
          name = 'r_'//trim(rxt_tag_lst(i)(1:))
       endif
       len = min(fieldname_len,len_trim(name))
       rate_names(i) = trim(name(1:len))
       call addfld(rate_names(i), (/ 'lev' /),'A', 'molecules/cm3/sec','reaction rate')
       if (history_scwaccm_forcing .and. rate_names(i) == 'r_O1D_H2O') then
          call add_default( rate_names(i), 1, ' ')
       endif
    enddo

    ! parse the terms of the summations
    grps_list => shr_exp_parse( rxn_rate_sums )
    deallocate( rxn_rate_sums )

    if (masterproc) write(iulog,*) 'rate_diags_init :'

    grp => grps_list
    do while(associated(grp))

       if (masterproc) then
          write(iulog,*) ' grp name : ',trim(grp%name)

          do j = 1, grp%n_terms
             write(iulog,'(f12.4,a,a)') grp%coeffs(j),' * ',trim(grp%vars(j))
          end do
       end if

       call addfld( grp%name, (/ 'lev' /),'A', 'molecules/cm3/sec','reaction rate group')

       grp => grp%next_item

    enddo

    o3_ndx = get_spc_ndx('O3')

  end subroutine rate_diags_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_calc( rxt_rates, vmr, m, ncol, lchnk )

    use mo_rxt_rates_conv, only: set_rates

    real(r8), intent(inout) :: rxt_rates(:,:,:) ! 'molec/cm3/sec'
    real(r8), intent(in)    :: vmr(:,:,:)
    real(r8), intent(in)    :: m(:,:)           ! air density (molecules/cm3)
    integer,  intent(in)    :: ncol, lchnk

    integer :: i, j, ndx
    real(r8) :: group_rate(ncol,pver)
    type(shr_exp_item_t), pointer :: grp

    call set_rates( rxt_rates, vmr, ncol )

    ! output individual tagged rates
    do i = 1, rxt_tag_cnt
       ! convert from vmr/sec to molecules/cm3/sec
       rxt_rates(:ncol,:,rxt_tag_map(i)) = rxt_rates(:ncol,:,rxt_tag_map(i)) * m(:ncol,:)
       call outfld( rate_names(i), rxt_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    ! output rate groups ( or families )

    grp => grps_list
    do while(associated(grp))

       group_rate(:,:) = 0._r8
       do j = 1, grp%n_terms
         ndx = lookup_tag_ndx(grp%vars(j))
         group_rate(:ncol,:) = group_rate(:ncol,:) + grp%coeffs(j)*rxt_rates(:ncol,:,ndx)
       enddo
       call outfld( grp%name, group_rate(:ncol,:), ncol, lchnk )

       grp => grp%next_item

    end do

  end subroutine rate_diags_calc

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  function rate_diags_o3s_loss( rxt_rates, vmr, ncol ) result(o3s_loss)
    use mo_rxt_rates_conv, only: set_rates
    use chem_mods,         only: rxntot

    real(r8), intent(in) :: rxt_rates(:,:,:)
    real(r8), intent(in) :: vmr(:,:,:)
    integer,  intent(in) :: ncol

    real(r8) :: o3s_loss(ncol,pver) ! /sec

    integer :: j, ndx
    real(r8) :: group_rate(ncol,pver)
    real(r8) :: lcl_rxt_rates(ncol,pver,rxntot)
    type(shr_exp_item_t), pointer :: grp

    o3s_loss(:,:) = 0._r8

    if (o3_ndx>0) then
       lcl_rxt_rates(:ncol,:,:) = rxt_rates(:ncol,:,:)
       call set_rates( lcl_rxt_rates, vmr, ncol )

       grp => grps_list
       loop: do while(associated(grp))

          if (trim(grp%name)=='O3S_Loss') then
             group_rate(:,:) = 0._r8
             do j = 1, grp%n_terms
                ndx = lookup_tag_ndx(grp%vars(j))
                group_rate(:ncol,:) = group_rate(:ncol,:) + grp%coeffs(j)*lcl_rxt_rates(:ncol,:,ndx)
             enddo
             o3s_loss(:ncol,:) = group_rate(:ncol,:)/vmr(:ncol,:,o3_ndx)
             exit loop
          endif
          grp => grp%next_item

       end do loop
    endif

  end function rate_diags_o3s_loss

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_final

    if (associated(grps_list)) call shr_exp_list_destroy(grps_list)

  end subroutine rate_diags_final

!-------------------------------------------------------------------
! Private routines :
!-------------------------------------------------------------------
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! finds the index corresponging to a given reacton name
!-------------------------------------------------------------------
  function lookup_tag_ndx( name ) result( ndx )
    character(len=*) :: name
    integer :: ndx

    integer :: i

    ndx = -1

    findloop: do i = 1,rxt_tag_cnt
       if (trim(name) .eq. trim(rate_names(i)(3:))) then
          ndx = i
          return
       endif
    end do findloop

    if (ndx<0) then
       call endrun('rate_diags: not able to find rxn tag name: '//trim(name))
    endif

  end function lookup_tag_ndx

end module rate_diags
