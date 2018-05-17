module trb_mtn_stress_cam

use shr_kind_mod, only: r8 => shr_kind_r8
use spmd_utils, only: masterproc
use cam_abortutils, only: endrun
use shr_log_mod, only: errMsg => shr_log_errMsg
use cam_logfile, only: iulog
use ppgrid, only: pcols, pver

implicit none
private

public :: trb_mtn_stress_readnl
public :: trb_mtn_stress_register
public :: trb_mtn_stress_init
public :: trb_mtn_stress_tend

! Is this module on at all?
logical :: do_tms = .false.

! Tuning parameters for TMS.
real(r8) :: tms_orocnst
real(r8) :: tms_z0fac

! pbuf field indices
integer :: &
     sgh30_idx = -1, &
     ksrftms_idx = -1, &
     tautmsx_idx = -1, &
     tautmsy_idx = -1

contains

subroutine trb_mtn_stress_readnl(nlfile)
  use namelist_utils, only: find_group_name
  use units, only: getunit, freeunit
  use spmd_utils, only: masterprocid, mpi_logical, mpi_real8, mpicom

  ! filepath for file containing namelist input
  character(len=*), intent(in) :: nlfile

  ! file unit and error code
  integer :: unitn, ierr

  character(len=*), parameter :: subname = "trb_mtn_stress_readnl"

  namelist /tms_nl/ do_tms, tms_orocnst, tms_z0fac

  ierr = 0

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'tms_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, tms_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  call mpi_bcast(do_tms,      1, mpi_logical, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(tms_orocnst, 1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(tms_z0fac,   1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")

end subroutine trb_mtn_stress_readnl

subroutine trb_mtn_stress_register()
  use physics_buffer, only: pbuf_add_field, dtype_r8

  call pbuf_add_field("ksrftms", "physpkg", dtype_r8, [pcols], ksrftms_idx)
  call pbuf_add_field("tautmsx", "physpkg", dtype_r8, [pcols], tautmsx_idx)
  call pbuf_add_field("tautmsy", "physpkg", dtype_r8, [pcols], tautmsy_idx)

end subroutine trb_mtn_stress_register

subroutine trb_mtn_stress_init()

  use cam_history, only: addfld, add_default, horiz_only
  use error_messages, only: handle_errmsg
  use phys_control, only: phys_getopts
  use physconst, only: karman, gravit, rair
  use physics_buffer, only: pbuf_get_index
  use trb_mtn_stress, only: init_tms

  logical :: history_amwg

  character(len=128) :: errstring

  if (.not. do_tms) return

  call phys_getopts(history_amwg_out=history_amwg)

  call init_tms( r8, tms_orocnst, tms_z0fac, karman, gravit, rair, errstring)
  call handle_errmsg(errstring, subname="init_tms")

  call addfld('TAUTMSX', horiz_only, 'A', 'N/m2', 'Zonal      turbulent mountain surface stress')
  call addfld('TAUTMSY', horiz_only, 'A', 'N/m2', 'Meridional turbulent mountain surface stress')
  if (history_amwg) then
     call add_default( 'TAUTMSX ', 1, ' ' )
     call add_default( 'TAUTMSY ', 1, ' ' )
  end if

  if (masterproc) then
     write(iulog,*)'Using turbulent mountain stress module'
     write(iulog,*)'  tms_orocnst = ',tms_orocnst
     write(iulog,*)'  tms_z0fac = ',tms_z0fac
  end if

  sgh30_idx = pbuf_get_index("SGH30")

end subroutine trb_mtn_stress_init

subroutine trb_mtn_stress_tend(state, pbuf, cam_in)
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field
  use physics_types, only: physics_state
  use camsrfexch, only: cam_in_t
  use cam_history, only: outfld
  use trb_mtn_stress, only: compute_tms

  type(physics_state), intent(in) :: state
  type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
  type(cam_in_t), intent(in) :: cam_in

  real(r8), pointer :: sgh30(:)
  real(r8), pointer :: ksrftms(:)
  real(r8), pointer :: tautmsx(:), tautmsy(:)

  call pbuf_get_field(pbuf, ksrftms_idx, ksrftms)
  call pbuf_get_field(pbuf, tautmsx_idx, tautmsx)
  call pbuf_get_field(pbuf, tautmsy_idx, tautmsy)

  if (.not. do_tms) then
     ksrftms = 0._r8
     tautmsx = 0._r8
     tautmsy = 0._r8
     return
  end if

  call pbuf_get_field(pbuf, sgh30_idx, sgh30)

  call compute_tms( pcols    , pver    , state%ncol , &
       state%u    , state%v  , state%t , state%pmid , & 
       state%exner, state%zm , sgh30   , ksrftms    , & 
       tautmsx    , tautmsy  , cam_in%landfrac )

  call outfld("TAUTMSX", tautmsx, pcols, state%lchnk)
  call outfld("TAUTMSY", tautmsy, pcols, state%lchnk)

end subroutine trb_mtn_stress_tend

end module trb_mtn_stress_cam
