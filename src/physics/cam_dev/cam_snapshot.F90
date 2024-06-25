module cam_snapshot
!--------------------------------------------------------
! The purpose of this module is to handle taking the "snapshot" of CAM data.
!
! This module writes out ALL the state, tend and pbuf fields.  It also includes the cam_in and cam_out
! fields which are used within CAM
!--------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use cam_history,    only: addfld, add_default, outfld
use cam_history,    only: cam_history_snapshot_deactivate, cam_history_snapshot_activate
use cam_history_support, only: horiz_only
use cam_abortutils, only: endrun
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_get_field_name
use physics_types,  only: physics_state, physics_tend, physics_ptend
use camsrfexch,     only: cam_out_t, cam_in_t
use ppgrid,         only: pcols, begchunk, endchunk
use constituents,   only: pcnst
use phys_control,   only: phys_getopts
use cam_logfile,    only: iulog
use cam_snapshot_common, only: snapshot_type, cam_snapshot_deactivate, cam_snapshot_all_outfld, cam_snapshot_ptend_outfld
use cam_snapshot_common, only: snapshot_type, cam_state_snapshot_init, cam_cnst_snapshot_init, cam_tend_snapshot_init
use cam_snapshot_common, only: cam_ptend_snapshot_init, cam_in_snapshot_init, cam_out_snapshot_init
use cam_snapshot_common, only: cam_pbuf_snapshot_init, snapshot_addfld

implicit  none

private

public :: cam_snapshot_init
public :: cam_snapshot_all_outfld_tphysbc, cam_snapshot_all_outfld_tphysac

private :: cam_tphysbc_snapshot_init, cam_tphysac_snapshot_init

integer :: ntphysbc_var
integer :: ntphysac_var

integer :: cam_snapshot_before_num, cam_snapshot_after_num

! Note the maximum number of variables for each type
type (snapshot_type)    ::  tphysbc_snapshot(30)
type (snapshot_type)    ::  tphysac_snapshot(30)

contains

subroutine cam_snapshot_init(cam_in_arr, cam_out_arr, pbuf, index)


!--------------------------------------------------------
! This subroutine does the addfld calls for ALL state, tend, ptend, and pbuf fields.  It also includes the cam_in and cam_out
! elements which are used within CAM
!--------------------------------------------------------
   type(cam_in_t),      intent(in)    :: cam_in_arr(begchunk:endchunk)
   type(cam_out_t),     intent(in)    :: cam_out_arr(begchunk:endchunk)
   type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:,:)
   integer,             intent(in)    :: index


   call phys_getopts(cam_snapshot_before_num_out = cam_snapshot_before_num, &
                     cam_snapshot_after_num_out  = cam_snapshot_after_num)


   ! Return if not turned on
   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   call cam_state_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)
   call cam_cnst_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)
   call cam_tend_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)
   call cam_ptend_snapshot_init(cam_snapshot_after_num)
   call cam_in_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num, cam_in_arr(index))
   call cam_out_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num, cam_out_arr(index))
   call cam_pbuf_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num, pbuf(:,index))
   call cam_tphysac_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)
   call cam_tphysbc_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

end subroutine cam_snapshot_init

subroutine cam_snapshot_all_outfld_tphysbc(file_num, state, tend, cam_in, cam_out, pbuf, cmfmc, cmfcme, &
        zdu, rliq, rice, dlf, dlf2, rliq2, net_flx)

use time_manager,   only: is_first_step, is_first_restart_step

!--------------------------------------------------------
! This subroutine does the outfld calls for ALL state, tend and pbuf fields for routines in tphysbc.
! It also includes the cam_in and cam_out elements which are used within CAM as well as variables which
! are local to tphysbc.
!--------------------------------------------------------

   integer,                            intent(in) :: file_num
   type(physics_state), intent(in) :: state
   type(physics_tend),  intent(in) :: tend
   type(cam_in_t),      intent(in) :: cam_in
   type(cam_out_t),     intent(in) :: cam_out
   type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
   real(r8),            intent(in) :: cmfmc(:,:)    ! convective mass flux
   real(r8),            intent(in) :: cmfcme(:,:)   ! cmf condensation - evaporation
   real(r8),            intent(in) :: zdu(:,:)      ! detraining mass flux from deep convection
   real(r8),            intent(in) :: rliq(:)       ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8),            intent(in) :: rice(:)       ! vertical integral of ice not yet in q(ixcldice)
   real(r8),            intent(in) :: dlf(:,:)      ! local copy of DLFZM (copy so need to output)
   real(r8),            intent(in) :: dlf2(:,:)     ! Detraining cld H20 from shallow convections
   real(r8),            intent(in) :: rliq2(:)      ! vertical integral of liquid from shallow scheme
   real(r8),            intent(in) :: net_flx(:)

   integer :: lchnk

   ! Return if the first timestep as not all fields may be filled in and this will cause a core dump
   if (is_first_step().or. is_first_restart_step()) return

   ! Return if not turned on
   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   lchnk = state%lchnk

   call outfld('tphysbc_cmfmc', cmfmc, pcols, lchnk)
   call outfld('tphysbc_cmfcme', cmfcme, pcols, lchnk)
   call outfld('tphysbc_zdu', zdu, pcols, lchnk)
   call outfld('tphysbc_rliq', rliq, pcols, lchnk)
   call outfld('tphysbc_rice', rice, pcols, lchnk)
   call outfld('tphysbc_dlf', dlf, pcols, lchnk)
   call outfld('tphysbc_dlf2', dlf2, pcols, lchnk)
   call outfld('tphysbc_rliq2', rliq2, pcols, lchnk)
   call outfld('tphysbc_net_flx', net_flx, pcols, lchnk)

   call cam_snapshot_all_outfld(file_num, state, tend, cam_in, cam_out, pbuf)

end subroutine cam_snapshot_all_outfld_tphysbc

subroutine cam_snapshot_all_outfld_tphysac(file_num, state, tend, cam_in, cam_out, pbuf, fh2o, surfric, obklen, flx_heat, &
       cmfmc, dlf, det_s, det_ice, net_flx)

use time_manager,   only: is_first_step

!--------------------------------------------------------
! This subroutine does the outfld calls for ALL state, tend and pbuf fields for routines in tphysac.
! It also includes the cam_in and cam_out elements which are used within CAM as well as variables which
! are local to tphysac.
!--------------------------------------------------------

   integer,                            intent(in) :: file_num
   type(physics_state), intent(in) :: state
   type(physics_tend),  intent(in) :: tend
   type(cam_in_t),      intent(in) :: cam_in
   type(cam_out_t),     intent(in) :: cam_out
   type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
   real(r8),            intent(in) :: fh2o(:)            ! h2o flux to balance source from methane chemistry
   real(r8),            intent(in) :: surfric(:)         ! surface friction velocity
   real(r8),            intent(in) :: obklen(:)          ! Obukhov length
   real(r8),            intent(in) :: flx_heat(:)        ! Heat flux for check_energy_chng.
   real(r8),            intent(in) :: cmfmc(:,:)    ! convective mass flux
   real(r8),            intent(in) :: dlf(:,:)      ! local copy of DLFZM (copy so need to output)
   real(r8),            intent(in) :: det_s(:)      ! vertical integral of detrained static energy from ice
   real(r8),            intent(in) :: det_ice(:)    ! vertical integral of detrained ice
   real(r8),            intent(in) :: net_flx(:)

   integer :: lchnk

   ! Return if the first timestep as not all fields may be filled in and this will cause a core dump
   if (is_first_step()) return

   ! Return if not turned on
   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   lchnk = state%lchnk

   call outfld('tphysac_fh2o', fh2o, pcols, lchnk)
   call outfld('tphysac_surfric', surfric, pcols, lchnk)
   call outfld('tphysac_obklen', obklen, pcols, lchnk)
   call outfld('tphysac_flx_heat', flx_heat, pcols, lchnk)
   call outfld('tphysac_cmfmc', cmfmc, pcols, lchnk)
   call outfld('tphysac_dlf', dlf, pcols, lchnk)
   call outfld('tphysac_det_s', det_s, pcols, lchnk)
   call outfld('tphysac_det_ice', det_ice, pcols, lchnk)
   call outfld('tphysac_net_flx', net_flx, pcols, lchnk)

   call cam_snapshot_all_outfld(file_num, state, tend, cam_in, cam_out, pbuf)

end subroutine cam_snapshot_all_outfld_tphysac

subroutine cam_tphysbc_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for the misc tphysbc physics variables that are passed individually
! into physics packages
!--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   ntphysbc_var = 0

   !--------------------------------------------------------
   ! Add the misc tphysbc variables to the output
   ! NOTE - flx_heat is added in tphysac
   !--------------------------------------------------------

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'flx',        'tphysbc_flx_heat',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cmfmc',        'tphysbc_cmfmc',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cmfcme',        'tphysbc_cmfcme',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'zdu',        'tphysbc_zdu',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rliq',        'tphysbc_rliq',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rice',        'tphysbc_rice',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'dlf',        'tphysbc_dlf',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'dlf2',        'tphysbc_dlf2',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rliq2',        'tphysbc_rliq2',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'net_flx',        'tphysbc_net_flx',         'unset',              horiz_only)


end subroutine cam_tphysbc_snapshot_init

subroutine cam_tphysac_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for the misc tphysac physics variables that are passed individually
! into physics packages
!--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   ntphysac_var = 0

   !--------------------------------------------------------
   ! Add the misc tphysac variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'fh2o',        'tphysac_fh2o',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'surfric',        'tphysac_surfric',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'obklen',        'tphysac_obklen',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'flx',        'tphysac_flx_heat',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cmfmc',        'tphysac_cmfmc',         'unset',              'lev')

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'zdu',        'tphysac_zdu',         'unset',              'lev')

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rliq',        'tphysac_rliq',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'dlf',        'tphysac_dlf',         'unset',              'lev')

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'dlf2',        'tphysac_dlf2',         'unset',              'lev')

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rliq2',        'tphysac_rliq2',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'det_s',        'tphysac_det_s',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'det_ice',        'tphysac_det_ice',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'net_flx',        'tphysac_net_flx',         'unset',              horiz_only)


end subroutine cam_tphysac_snapshot_init

end module cam_snapshot
