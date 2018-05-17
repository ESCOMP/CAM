module beljaars_drag_cam

use shr_kind_mod, only: r8 => shr_kind_r8
use spmd_utils, only: masterproc
use cam_abortutils, only: endrun
use shr_log_mod, only: errMsg => shr_log_errMsg
use cam_logfile, only: iulog
use ppgrid, only: pcols, pver

implicit none
private

public :: beljaars_drag_readnl
public :: beljaars_drag_register
public :: beljaars_drag_init
public :: beljaars_drag_tend

! Is this module on at all?
logical, public, protected :: do_beljaars = .false.

! Tuning parameters for TMS.
real(r8) :: blj_orocnst
real(r8) :: blj_z0fac

! pbuf field indices
integer :: &
     sgh30_idx = -1, &
     dragblj_idx = -1, &
     taubljx_idx = -1, &
     taubljy_idx = -1

contains

subroutine beljaars_drag_readnl(nlfile)
  use namelist_utils, only: find_group_name
  use units, only: getunit, freeunit
  use spmd_utils, only: masterprocid, mpi_logical, mpi_real8, mpicom

  ! filepath for file containing namelist input
  character(len=*), intent(in) :: nlfile

  ! file unit and error code
  integer :: unitn, ierr

  character(len=*), parameter :: subname = "beljaars_drag_readnl"

  namelist /blj_nl/ do_beljaars

  ierr = 0

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'blj_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, blj_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  call mpi_bcast(do_beljaars,      1, mpi_logical, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")

end subroutine beljaars_drag_readnl

subroutine beljaars_drag_register()
  use physics_buffer, only: pbuf_add_field, dtype_r8

  call pbuf_add_field("dragblj", "physpkg", dtype_r8, (/pcols,pver/), dragblj_idx)
  call pbuf_add_field("taubljx", "physpkg", dtype_r8, (/pcols/), taubljx_idx)
  call pbuf_add_field("taubljy", "physpkg", dtype_r8, (/pcols/), taubljy_idx)

end subroutine beljaars_drag_register

subroutine beljaars_drag_init()

  use cam_history, only: addfld, add_default, horiz_only
  use error_messages, only: handle_errmsg
  use phys_control, only: phys_getopts
  use physconst, only: karman, gravit, rair
  use physics_buffer, only: pbuf_get_index
  use beljaars_drag, only: init_blj

  logical :: history_amwg

  character(len=128) :: errstring

  if (.not. do_beljaars) return

  call phys_getopts(history_amwg_out=history_amwg)

  call init_blj( r8, gravit, rair, errstring )
  call handle_errmsg(errstring, subname="init_blj")

  call addfld('DRAGBLJ', (/ 'lev' /) , 'A', '1/s', 'Drag profile from Beljaars SGO              ')
  call addfld('TAUBLJX', horiz_only, 'A', 'N/m2',  'Zonal      integrated drag from Beljaars SGO')
  call addfld('TAUBLJY', horiz_only, 'A', 'N/m2',  'Meridional integrated drag from Beljaars SGO')
  if (history_amwg) then
     call add_default( 'TAUBLJX ', 1, ' ' )
     call add_default( 'TAUBLJY ', 1, ' ' )
  end if

  if (masterproc) then
     write(iulog,*)'Using Beljaars SGO drag module'
  end if

  sgh30_idx = pbuf_get_index("SGH30")

end subroutine beljaars_drag_init

subroutine beljaars_drag_tend(state, pbuf, cam_in)
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field
  use physics_types, only: physics_state
  use camsrfexch, only: cam_in_t
  use cam_history, only: outfld
  use beljaars_drag, only: compute_blj

  type(physics_state), intent(in) :: state
  type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
  type(cam_in_t), intent(in) :: cam_in

  real(r8), pointer :: sgh30(:)
  real(r8), pointer :: dragblj(:,:)
  real(r8), pointer :: taubljx(:), taubljy(:)

  call pbuf_get_field(pbuf, dragblj_idx, dragblj)
  call pbuf_get_field(pbuf, taubljx_idx, taubljx)
  call pbuf_get_field(pbuf, taubljy_idx, taubljy)

  if (.not. do_beljaars) then
     dragblj = 0._r8
     taubljx = 0._r8
     taubljy = 0._r8
     return
  end if

  call pbuf_get_field(pbuf, sgh30_idx, sgh30)

  call compute_blj( pcols    , pver    , state%ncol , &
       state%u    , state%v  , state%t , state%pmid , & 
       state%pdel , state%zm , sgh30   , dragblj    , & 
       taubljx    , taubljy  , cam_in%landfrac )

  call outfld("TAUBLJX", taubljx, pcols, state%lchnk)
  call outfld("TAUBLJY", taubljy, pcols, state%lchnk)
  call outfld("DRAGBLJ", dragblj, pcols, state%lchnk)

end subroutine beljaars_drag_tend

end module beljaars_drag_cam
