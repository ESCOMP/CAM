module gw_rdg_cam

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,only: masterproc
use cam_abortutils, only: endrun


implicit none
private
save

! Public interface
public :: gw_rdg_cam_readnl

!==========================================================================
contains
!==========================================================================

!!$  subroutine gw_rdg_cam_init(band)
!!$    use gw_rdg, only :: gw_rdg_init
!!$    type(GWBand), intent(in)                    :: band
!!$    logical, intent(in)                         :: use_gw_rdg_beta, use_gw_rdg_gamma
!!$
!!$    call gw_rdg_init( band, &
!!$       prdg_in, &
!!$       rearth_c, &
!!$       effgw_rdg_beta, &
!!$       effgw_rdg_gamma, &
!!$       use_gw_rdg_beta, &
!!$       use_gw_rdg_gamma, &
!!$       bnd_topo_file, &
!!$       bnd_rdg_file, &
!!$       gw_rdg_do_divstream, gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, &
!!$       gw_rdg_Frx0, gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, &
!!$       gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
!!$       gw_rdg_do_backward_compat, gw_rdg_orohmin, gw_rdg_orovmin, &
!!$       gw_rdg_orostratmin, gw_rdg_orom2min, gw_rdg_do_vdiff, &
!!$       masterproc, iulog, errmsg, errflg)
!!$  end subroutine gw_rdg_cam_init
subroutine gw_rdg_cam_readnl(nlfile)
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, mpi_logical
  use gw_rdg,          only: gw_rdg_init
  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gw_rdg_readnl'





  !----------------------------------------------------------------------




end subroutine gw_rdg_cam_readnl

end module gw_rdg_cam
