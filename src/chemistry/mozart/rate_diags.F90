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
  use cam_abortutils,   only : endrun
  use sums_utils,       only : sums_grp_t, parse_sums

  implicit none
  private 
  public :: rate_diags_init
  public :: rate_diags_calc
  public :: rate_diags_readnl

  character(len=fieldname_len) :: rate_names(rxt_tag_cnt)

  integer :: ngrps = 0
  type(sums_grp_t), allocatable :: grps(:)  

  integer, parameter :: maxlines = 200
  character(len=CL), allocatable :: rxn_rate_sums(:)

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

    integer :: i, len, pos
    character(len=64) :: name
    logical :: history_scwaccm_forcing

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
    call parse_sums(rxn_rate_sums, ngrps, grps)
    deallocate( rxn_rate_sums )

    do i = 1, ngrps
       call addfld( grps(i)%name, (/ 'lev' /),'A', 'molecules/cm3/sec','reaction rate group')
    enddo

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

    call set_rates( rxt_rates, vmr, ncol )

    ! output individual tagged rates    
    do i = 1, rxt_tag_cnt
       ! convert from vmr/sec to molecules/cm3/sec
       rxt_rates(:ncol,:,rxt_tag_map(i)) = rxt_rates(:ncol,:,rxt_tag_map(i)) * m(:ncol,:)
       call outfld( rate_names(i), rxt_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    ! output rate groups ( or families )
    do i = 1, ngrps
       group_rate(:,:) = 0._r8
       do j = 1, grps(i)%nmembers
         ndx = lookup_tag_ndx(grps(i)%term(j))
         group_rate(:ncol,:) = group_rate(:ncol,:) + grps(i)%multipler(j)*rxt_rates(:ncol,:,ndx)
       enddo 
       call outfld( grps(i)%name, group_rate(:ncol,:), ncol, lchnk )       
    end do

  end subroutine rate_diags_calc

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
