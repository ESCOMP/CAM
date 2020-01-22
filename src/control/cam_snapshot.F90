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

implicit  none

private

public :: cam_snapshot_init, cam_snapshot_deactivate
public :: cam_snapshot_all_outfld, cam_snapshot_all_outfld_tphysbc, cam_snapshot_all_outfld_tphysac
public :: cam_snapshot_ptend_outfld

! This is the number of pbuf fields in the CAM code that are declared with the fieldname as opposed to being data driven.
integer, parameter :: npbuf_all = 327

type snapshot_type
  character(len=40)  :: ddt_string
  character(len=256) :: standard_name
  character(len=20)  :: dim_name
  character(len=8)   :: units
end type snapshot_type

type snapshot_type_nd
  character(len=40)  :: ddt_string
  character(len=256) :: standard_name
  character(len=20)  :: dim_name(6) ! hardwired 6 potential dimensions in pbuf
  character(len=8)   :: units
end type snapshot_type_nd

type pbuf_info_type
  character(len=40)  :: name
  character(len=256) :: standard_name
  character(len=8)   :: units
  character(len=100) :: dim_string(6) ! hardwired 6 potential dimensions in pbuf
end type pbuf_info_type

integer :: nstate_var
integer :: ncnst_var
integer :: ntend_var
integer :: ncam_in_var
integer :: ncam_out_var
integer :: npbuf_var
integer :: ntphysbc_var
integer :: ntphysac_var

integer :: cam_snapshot_before_num, cam_snapshot_after_num

! Note the maximum number of variables for each type
type (snapshot_type)    ::  state_snapshot(27)
type (snapshot_type)    ::  cnst_snapshot(pcnst)
type (snapshot_type)    ::  tend_snapshot(6)
type (snapshot_type)    ::  cam_in_snapshot(30)
type (snapshot_type)    ::  cam_out_snapshot(30)
type (snapshot_type)    ::  tphysbc_snapshot(30)
type (snapshot_type)    ::  tphysac_snapshot(30)
type (snapshot_type_nd) ::  pbuf_snapshot(250)

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

subroutine cam_snapshot_all_outfld_tphysbc(file_num, state, tend, cam_in, cam_out, pbuf, flx_heat, cmfmc, cmfcme, &
        pflx, zdu, rliq, rice, dlf, dlf2, rliq2, det_s, det_ice, net_flx)

use time_manager,   only: is_first_step

!--------------------------------------------------------
! This subroutine does the outfld calls for ALL state, tend and pbuf fields for routines in tphysbc.
! It also includes the cam_in and cam_out elements which are used within CAM as well as variables which
! are local to tphysac.
!--------------------------------------------------------

   integer,                            intent(in) :: file_num
   type(physics_state), intent(in) :: state
   type(physics_tend),  intent(in) :: tend
   type(cam_in_t),      intent(in) :: cam_in
   type(cam_out_t),     intent(in) :: cam_out
   type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
   real(r8),            intent(in) :: flx_heat(:)   ! Heat flux for check_energy_chng.
   real(r8),            intent(in) :: cmfmc(:,:)    ! convective mass flux
   real(r8),            intent(in) :: cmfcme(:,:)   ! cmf condensation - evaporation
   real(r8),            intent(in) :: pflx(:,:)     ! convective rain flux throughout bottom of level
   real(r8),            intent(in) :: zdu(:,:)      ! detraining mass flux from deep convection
   real(r8),            intent(in) :: rliq(:)       ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8),            intent(in) :: rice(:)       ! vertical integral of ice not yet in q(ixcldice)
   real(r8),            intent(in) :: dlf(:,:)      ! local copy of DLFZM (copy so need to output)
   real(r8),            intent(in) :: dlf2(:,:)     ! Detraining cld H20 from shallow convections
   real(r8),            intent(in) :: rliq2(:)      ! vertical integral of liquid from shallow scheme
   real(r8),            intent(in) :: det_s(:)      ! vertical integral of detrained static energy from ice
   real(r8),            intent(in) :: det_ice(:)    ! vertical integral of detrained ice
   real(r8),            intent(in) :: net_flx(:)

   integer :: lchnk

   ! Return if the first timestep as not all fields may be filled in and this will cause a core dump
   if (is_first_step()) return

   ! Return if not turned on 
   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   lchnk = state%lchnk

   call outfld('flx_heat_snapshot', flx_heat, pcols, lchnk)
   call outfld('cmfmc_snapshot', cmfmc, pcols, lchnk)
   call outfld('cmfcme_snapshot', cmfcme, pcols, lchnk)
   call outfld('pflx_snapshot', pflx, pcols, lchnk)
   call outfld('zdu_snapshot', zdu, pcols, lchnk)
   call outfld('rliq_snapshot', rliq, pcols, lchnk)
   call outfld('rice_snapshot', rice, pcols, lchnk)
   call outfld('dlf_snapshot', dlf, pcols, lchnk)
   call outfld('dlf2_snapshot', dlf2, pcols, lchnk)
   call outfld('rliq2_snapshot', rliq2, pcols, lchnk)
   call outfld('det_s_snapshot', det_s, pcols, lchnk)
   call outfld('det_ice_snapshot', det_ice, pcols, lchnk)
   call outfld('net_flx_snapshot', net_flx, pcols, lchnk)

   call cam_snapshot_all_outfld(file_num, state, tend, cam_in, cam_out, pbuf)

end subroutine cam_snapshot_all_outfld_tphysbc

subroutine cam_snapshot_all_outfld_tphysac(file_num, state, tend, cam_in, cam_out, pbuf, fh2o, surfric, obklen, flx_heat)

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

   integer :: lchnk

   ! Return if the first timestep as not all fields may be filled in and this will cause a core dump
   if (is_first_step()) return

   ! Return if not turned on 
   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   lchnk = state%lchnk

   call outfld('fh2o_snapshot', fh2o, pcols, lchnk)
   call outfld('surfric_snapshot', surfric, pcols, lchnk)
   call outfld('obklen_snapshot', obklen, pcols, lchnk)
   call outfld('flx_heat_snapshot', flx_heat, pcols, lchnk)

   call cam_snapshot_all_outfld(file_num, state, tend, cam_in, cam_out, pbuf)

end subroutine cam_snapshot_all_outfld_tphysac

subroutine cam_snapshot_all_outfld(file_num, state, tend, cam_in, cam_out, pbuf)

use time_manager,   only: is_first_step

!--------------------------------------------------------
! This subroutine does the outfld calls for ALL state, tend and pbuf fields.  It also includes the cam_in and cam_out
! elements which are used within CAM
!--------------------------------------------------------

   integer,                            intent(in) :: file_num
   type(physics_state), intent(in) :: state
   type(physics_tend),  intent(in) :: tend
   type(cam_in_t),      intent(in) :: cam_in
   type(cam_out_t),     intent(in) :: cam_out
   type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)


   integer :: lchnk

   ! Return if the first timestep as not all fields may be filled in and this will cause a core dump
   if (is_first_step()) return

   ! Return if not turned on
   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   lchnk = state%lchnk

   ! Write out all the state fields
   call state_snapshot_all_outfld(lchnk, file_num, state)

   ! Write out all the constituent fields
   call cnst_snapshot_all_outfld(lchnk, file_num, state%q)

   ! Write out all the tendency fields
   call tend_snapshot_all_outfld(lchnk, file_num, tend)

   ! Write out all the cam_in fields
   call cam_in_snapshot_all_outfld(lchnk, file_num, cam_in)

   ! Write out all the cam_out fields
   call cam_out_snapshot_all_outfld(lchnk, file_num, cam_out)

   ! Write out all the pbuf fields
   call cam_pbuf_snapshot_all_outfld(lchnk, file_num, pbuf)

end subroutine cam_snapshot_all_outfld

subroutine cam_snapshot_deactivate()

!--------------------------------------------------------
! This subroutine deactivates the printing of the snapshot before and after files
! Note - this needs to be done as add_default has been called to setup the proper
!        outputting of the requested fields. The outfld calls will only write
!        one file at a time (using the same name in both files), hence the writing
!        needs to be turned off for all fields, and will be turned on individaully
!        when needed.
!--------------------------------------------------------
   integer :: i

   ! Return if not turned on
   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   do i=1,nstate_var
      call cam_history_snapshot_deactivate(state_snapshot(i)%standard_name)
   end do

   do i=1,ncnst_var
      call cam_history_snapshot_deactivate(cnst_snapshot(i)%standard_name)
   end do

   do i=1,ntend_var
      call cam_history_snapshot_deactivate(tend_snapshot(i)%standard_name)
   end do

   do i=1,ncam_in_var
      call cam_history_snapshot_deactivate(cam_in_snapshot(i)%standard_name)
   end do

   do i=1,ncam_out_var
      call cam_history_snapshot_deactivate(cam_out_snapshot(i)%standard_name)
   end do

   do i=1,npbuf_var
      call cam_history_snapshot_deactivate(pbuf_snapshot(i)%standard_name)
   end do

end subroutine cam_snapshot_deactivate


subroutine cam_state_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for state
!--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   nstate_var = 0

   !--------------------------------------------------------
   ! Add the state variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%ps',        'ps_snapshot',         'Pa',              horiz_only)

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%psdry',     'psdry_snapshot',      'Pa',              horiz_only)

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%phis',      'phis_snapshot',       'm2/m2',           horiz_only)

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%t',         't_snapshot',          'K',               'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%u',         'u_snapshot',          'm s-1',           'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%v',         'v_snapshot',          'm s-1',           'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%s',         's_snapshot',          ' ',               'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%omega',     'omega_snapshot',      'Pa s-1',           'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pmid',      'pmid_snapshot',       'Pa',               'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pmiddry',   'pmiddry_snapshot',    'Pa',               'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pdel',      'pdel_snapshot',       'Pa',               'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pdeldry',   'pdeldry_snapshot',    'Pa',               'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%rpdel',     'rpdel_snapshot',      'Pa',               'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%rpdeldry',  'rpdeldry_snapshot',   'Pa',              'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpmid',    'lnpmid_snapshot',     'unset',           'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpmiddry', 'lnpmiddry_snapshot',  'unset',           'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%exner',     'exner_snapshot',      'unset',            'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%zm',        'zm_snapshot',         'm',                'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pint',      'pint_snapshot',       'Pa',               'ilev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pintdry',   'pintdry_snapshot',    'Pa',               'ilev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpint',    'lnpint_snapshot',     'unset',            'ilev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpintdry', 'lnpintdry_snapshot',  'unset',            'ilev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%zi',        'zi_snapshot',         'm',                'ilev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%te_ini',    'te_ini_snapshot',     'unset',            horiz_only)

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%te_cur',    'te_cur_snapshot',     'unset',            horiz_only)

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%tw_ini',    'tw_ini_snapshot',     'unset',            horiz_only)

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%tw_cur',    'tw_cur_snapshot',     'unset',            horiz_only)

end subroutine cam_state_snapshot_init

subroutine cam_cnst_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for state constituent (q) fields
!--------------------------------------------------------

   use constituents, only: cnst_name, cnst_longname

   integer, intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   !--------------------------------------------------------
   ! Add the cnst variables to the output
   !--------------------------------------------------------

   ncnst_var = 0 ! Updated inside snapshot_addfld

   do while (ncnst_var < pcnst)
      call snapshot_addfld(ncnst_var, cnst_snapshot, cam_snapshot_before_num, &
           cam_snapshot_after_num, cnst_name(ncnst_var+1),                    &
           trim(cnst_name(ncnst_var+1))//'_snapshot', 'kg kg-1', 'lev')
   end do

end subroutine cam_cnst_snapshot_init

subroutine cam_tend_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for tend fields.
!--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   ntend_var = 0

   !--------------------------------------------------------
   ! Add the physics_tend variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%dtdt',        'dtdt_snapshot',         'K s-1',    'lev')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%dudt',        'dudt_snapshot',         '',    'lev')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%dvdt',        'dvdt_snapshot',         '',    'lev')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%flx_net',        'flx_net_snapshot',   '',    horiz_only)

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%te_tnd',        'te_tnd_snapshot',   '',    horiz_only)

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%tw_tnd',        'tw_tnd_snapshot',   '',    horiz_only)

end subroutine cam_tend_snapshot_init

subroutine cam_ptend_snapshot_init(cam_snapshot_after_num)
   use constituents, only: cnst_name, cnst_longname

   !--------------------------------------------------------
   ! This subroutine does the addfld calls for ptend fields.
   !--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_after_num

   integer            :: mcnst
   character(len=64)  :: fname
   character(len=128) :: lname
   character(len=32)  :: cam_take_snapshot_before
   character(len=32)  :: cam_take_snapshot_after

   call phys_getopts(cam_take_snapshot_before_out = cam_take_snapshot_before, &
        cam_take_snapshot_after_out = cam_take_snapshot_after)

   if (trim(cam_take_snapshot_before) == trim(cam_take_snapshot_after)) then

      !--------------------------------------------------------
      ! Add the physics_ptend variables to the output
      !--------------------------------------------------------

      call addfld('ptend_s_snapshot', (/ 'lev' /), 'I', 'J kg-1 s-1',         &
           'heating rate snapshot')
      call add_default('ptend_s_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_u_snapshot', (/ 'lev' /), 'I', 'm s-1 s-1',          &
           'momentum tendency snapshot')
      call add_default('ptend_u_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_v_snapshot', (/ 'lev' /), 'I', 'm s-1 s-1',          &
           'momentum tendency snapshot')
      call add_default('ptend_v_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_hflux_srf_snapshot', horiz_only, 'I', 'W m-2',       &
           'net zonal stress at surface snapshot')
      call add_default('ptend_hflux_srf_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_hflux_top_snapshot', horiz_only, 'I', 'W m-2',       &
           'net zonal stress at top of model snapshot')
      call add_default('ptend_hflux_top_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_taux_srf_snapshot', horiz_only, 'I', 'Pa',           &
           'net meridional stress at surface snapshot')
      call add_default('ptend_taux_srf_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_taux_top_snapshot', horiz_only, 'I', 'Pa',           &
           'net zonal stress at top of model snapshot')
      call add_default('ptend_taux_top_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_tauy_srf_snapshot', horiz_only, 'I', 'Pa',           &
           'net meridional stress at surface snapshot')
      call add_default('ptend_tauy_srf_snapshot', cam_snapshot_after_num, ' ')

      call addfld('ptend_tauy_top_snapshot', horiz_only, 'I', 'Pa',           &
           'net meridional stress at top of model snapshot')
      call add_default('ptend_tauy_top_snapshot', cam_snapshot_after_num, ' ')

      do mcnst = 1, pcnst
         fname = 'ptend_'//trim(cnst_name(mcnst))//'_snapshot'
         lname = 'tendency of '//trim(cnst_longname(mcnst))//' snapshot'
         call addfld(trim(fname), (/ 'lev' /), 'I', 'kg kg-1 s-1', trim(lname))
         call add_default(trim(fname), cam_snapshot_after_num, ' ')

         fname = 'cflx_srf_'//trim(cnst_name(mcnst))//'_snapshot'
         lname = 'flux of '//trim(cnst_longname(mcnst))//' at surface snapshot'
         call addfld(trim(fname), horiz_only, 'I', 'kg m-2 s-1', trim(lname))
         call add_default(trim(fname), cam_snapshot_after_num, ' ')

         fname = 'cflx_top_'//trim(cnst_name(mcnst))//'_snapshot'
         lname = 'flux of '//trim(cnst_longname(mcnst))//' at top of model snapshot'
         call addfld(trim(fname), horiz_only, 'I', 'kg m-2 s-1', trim(lname))
         call add_default(trim(fname), cam_snapshot_after_num, ' ')
      end do

   end if

end subroutine cam_ptend_snapshot_init

subroutine cam_in_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num, cam_in)

!--------------------------------------------------------
! This subroutine does the addfld calls for cam_in fields
!--------------------------------------------------------

   type(cam_in_t), intent(in) :: cam_in

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   ncam_in_var = 0

   !--------------------------------------------------------
   ! Add the state variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%landfrac',        'landfrac_snapshot',          'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%ocnfrac',         'ocnfrac_snapshot',           'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%snowhland',       'snowhland_snapshot',         'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%ts',              'ts_snapshot',                'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%sst',             'sst_snapshot',               'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%icefrac',         'icefrac_snapshot',           'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%shf',             'shf_snapshot',               'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%cflx',            'cflx_snapshot',              'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%wsx',             'wsx_snapshot',               'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%wsy',             'wsy_snapshot',               'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%asdif',           'asdif_snapshot',             'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%aldif',           'aldif_snapshot',             'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%lwup',            'lwup_snapshot',              'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%asdir',           'asdir_snapshot',             'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%aldir',           'aldir_snapshot',             'unset',          horiz_only)

    if (associated (cam_in%meganflx)) &
    call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
      'cam_in%meganflx',        'meganflx_snapshot',          'unset',          horiz_only)

    if (associated (cam_in%fireflx)) &
    call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
      'cam_in%fireflx',         'fireflx_snapshot',           'unset',          horiz_only)

    if (associated (cam_in%fireztop)) &
    call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
      'cam_in%fireztop',        'fireztop_snapshot',          'unset',          horiz_only)

    if (associated (cam_in%depvel)) &
    call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
      'cam_in%depvel',          'depvel_snapshot',            'unset',          horiz_only)

   call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_in%lhf',             'lhf_snapshot',               'unset',          horiz_only)

    if (associated (cam_in%fv)) &
    call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
      'cam_in%fv',              'fv_snapshot',                'unset',          horiz_only)

    if (associated (cam_in%ram1)) &
    call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
      'cam_in%ram1',            'ram1_snapshot',              'unset',          horiz_only)

    if (associated (cam_in%dstflx)) &
    call snapshot_addfld( ncam_in_var, cam_in_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
      'cam_in%dstflx',          'dstflx_snapshot',            'unset',          horiz_only)

end subroutine cam_in_snapshot_init

subroutine cam_out_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num, cam_out)

!--------------------------------------------------------
! This subroutine does the addfld calls for cam_out fields
!--------------------------------------------------------

   type(cam_out_t),  intent(in) :: cam_out

   integer,          intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   ncam_out_var = 0

   !--------------------------------------------------------
   ! Add the state variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precc',               'precc_snapshot',           'm s-1',          horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precl',               'precl_snapshot',           'm s-1',          horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precsc',               'precsc_snapshot',         'm s-1',          horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precsl',               'precsl_snapshot',         'm s-1',          horiz_only)

   if (associated(cam_out%nhx_nitrogen_flx)) &
   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%nhx_nitrogen_flx',     'nhx_nitro_flx_snapshot',  'kgN m2-1 sec-1', horiz_only)

   if (associated(cam_out%noy_nitrogen_flx)) &
   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%noy_nitrogen_flx',     'noy_nitro_flx_snapshot',  'kgN m2-1 sec-1', horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%bcphodry',             'bcphodry_snapshot',       'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%bcphidry',             'bcphidry_snapshot',       'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%ocphodry',             'ocphodry_snapshot',       'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%ocphidry',             'ocphidry_snapshot',       'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%bcphiwet',             'bcphiwet_snapshot',       'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%ocphiwet',             'ocphiwet_snapshot',       'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet1',             'dstwet1_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet2',             'dstwet2_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet3',             'dstwet3_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet4',             'dstwet4_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry1',             'dstdry1_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry2',             'dstdry2_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry3',             'dstdry3_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry4',             'dstdry4_snapshot',         'kg m-2 s-1',     horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%sols',                'sols_snapshot',            'W m-2',          horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%soll',                'soll_snapshot',            'W m-2',          horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%solsd',               'solsd_snapshot',           'W m-2',          horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%solld',               'solld_snapshot',           'W m-2',          horiz_only)

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%netsw',               'netsw_snapshot',           'unset',          horiz_only)

end subroutine cam_out_snapshot_init

subroutine cam_pbuf_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num, pbuf)

!--------------------------------------------------------
! This subroutine does the addfld calls for pbuf fields.
!--------------------------------------------------------

   use physics_buffer, only: pbuf_get_dim_strings

   integer,                   intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num
   type(physics_buffer_desc), intent(in) :: pbuf(:)

   integer :: i, j, npbuf
   type(pbuf_info_type) :: pbuf_info(size(pbuf))
   character(len=40) :: const_cname(ncnst_var)
   character(len=40) :: dim_strings(size(pbuf),6) ! Hardwired 6 potential dimensions in pbuf

   npbuf = size(pbuf(:))

   !--------------------------------------------------------
   ! fill the name, standard name and units for pbuf_info
   !--------------------------------------------------------

   call fill_pbuf_info(pbuf_info, pbuf, const_cname)

   !--------------------------------------------------------
   ! Determine the indices for the addfld call based on the dimensions in the pbuf
   !--------------------------------------------------------

   call pbuf_get_dim_strings(pbuf, dim_strings)
   do i=1, npbuf
      ! If the second dimension is empty, then this is a horiz_only field
      if (trim(dim_strings(i,2)) == '') then
         pbuf_info(i)%dim_string(1) = horiz_only
      else
         ! The first dimension is the horizontal dimension and should not be used in the addfld call
         do j=2,6
            pbuf_info(i)%dim_string(j-1) = dim_strings(i,j)
         end do
      end if
   end do

   !--------------------------------------------------------
   ! Now that all of the information for the pbufs is stored, call the addfld
   !--------------------------------------------------------
   npbuf_var = 0 ! Updated inside snapshot_addfld

   do while (npbuf_var < npbuf)
      call snapshot_addfld_nd( npbuf_var, pbuf_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
        pbuf_info(npbuf_var+1)%name,   pbuf_info(npbuf_var+1)%standard_name,   pbuf_info(npbuf_var+1)%units,&
        pbuf_info(npbuf_var+1)%dim_string)
   end do

end subroutine cam_pbuf_snapshot_init

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
     'cmfmc',        'cmfmc_snapshot',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cmfcme',        'cmfcme_snapshot',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'pflx',        'pflx_snapshot',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'zdu',        'zdu_snapshot',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rliq',        'rliq_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rice',        'rice_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'dlf',        'dlf_snapshot',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'dlf2',        'dlf2_snapshot',         'unset',              'lev')

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'rliq2',        'rliq2_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'det_s',        'det_s_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'det_ice',        'det_ice_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysbc_var, tphysbc_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'net_flx',        'net_flx_snapshot',         'unset',              horiz_only)


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
     'fh2o',        'fh2o_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'surfric',        'surfric_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'obklen',        'obklen_snapshot',         'unset',              horiz_only)

   call snapshot_addfld( ntphysac_var, tphysac_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'flx',        'flx_heat_snapshot',         'unset',              horiz_only)

end subroutine cam_tphysac_snapshot_init

subroutine snapshot_addfld_nd(nddt_var, ddt_snapshot, cam_snapshot_before_num, cam_snapshot_after_num,&
    ddt_string, standard_name, units, dimension_string)

   integer,                 intent(inout)  :: nddt_var
   type (snapshot_type_nd), intent(inout) ::  ddt_snapshot(:)


   integer,          intent(in) :: cam_snapshot_before_num
   integer,          intent(in) :: cam_snapshot_after_num
   character(len=*), intent(in) :: ddt_string
   character(len=*), intent(in) :: standard_name
   character(len=*), intent(in) :: units
   character(len=*), intent(in) :: dimension_string(:)

   integer :: ndims

   nddt_var=nddt_var+1

   if (nddt_var > size(ddt_snapshot)) &
      call endrun(' ERROR in snapshot_addfld: ddt_snapshot array not allocated large enough')

   ndims = count(dimension_string /= '')

   if (trim(dimension_string(1)) == horiz_only) then
      call addfld(standard_name, horiz_only, 'I', units, standard_name)
   else
      call addfld(standard_name, dimension_string(1:ndims), 'I', units, standard_name)
   end if
   if (cam_snapshot_before_num > 0) call add_default(standard_name, cam_snapshot_before_num, ' ')
   if (cam_snapshot_after_num > 0)  call add_default(standard_name, cam_snapshot_after_num, ' ')

   ddt_snapshot(nddt_var)%ddt_string    = ddt_string
   ddt_snapshot(nddt_var)%standard_name = standard_name
   ddt_snapshot(nddt_var)%dim_name(:)   = dimension_string(:)
   ddt_snapshot(nddt_var)%units         = units


end subroutine snapshot_addfld_nd

subroutine snapshot_addfld(nddt_var, ddt_snapshot, cam_snapshot_before_num, cam_snapshot_after_num,&
    ddt_string, standard_name, units, dimension_string)

   integer,              intent(inout)  :: nddt_var
   type (snapshot_type), intent(inout) ::  ddt_snapshot(:)


   integer,          intent(in) :: cam_snapshot_before_num
   integer,          intent(in) :: cam_snapshot_after_num
   character(len=*), intent(in) :: ddt_string
   character(len=*), intent(in) :: standard_name
   character(len=*), intent(in) :: units
   character(len=*), intent(in) :: dimension_string


   nddt_var=nddt_var+1

   if (nddt_var > size(ddt_snapshot)) &
      call endrun(' ERROR in snapshot_addfld: ddt_snapshot array not allocated large enough')

   call addfld(standard_name, dimension_string, 'I', units, standard_name)
   if (cam_snapshot_before_num > 0) call add_default(standard_name, cam_snapshot_before_num, ' ')
   if (cam_snapshot_after_num > 0)  call add_default(standard_name, cam_snapshot_after_num, ' ')

   ddt_snapshot(nddt_var)%ddt_string    = ddt_string
   ddt_snapshot(nddt_var)%standard_name = standard_name
   ddt_snapshot(nddt_var)%dim_name      = dimension_string
   ddt_snapshot(nddt_var)%units         = units


end subroutine snapshot_addfld

subroutine state_snapshot_all_outfld(lchnk, file_num, state)

   integer,              intent(in)  :: lchnk
   integer,              intent(in)  :: file_num
   type(physics_state),  intent(in)  :: state

   integer :: i

   do i=1, nstate_var

      ! Turn on the writing for only the requested tape (file_num)
      call cam_history_snapshot_activate(trim(state_snapshot(i)%standard_name), file_num)

      ! Select the state field which is being written
      select case(state_snapshot(i)%ddt_string)

      case ('state%ps')
         call outfld(state_snapshot(i)%standard_name, state%ps, pcols, lchnk)

      case ('state%psdry')
         call outfld(state_snapshot(i)%standard_name, state%psdry, pcols, lchnk)

      case ('state%phis')
         call outfld(state_snapshot(i)%standard_name, state%phis, pcols, lchnk)

      case ('state%t')
         call outfld(state_snapshot(i)%standard_name, state%t, pcols, lchnk)

      case ('state%u')
         call outfld(state_snapshot(i)%standard_name, state%u, pcols, lchnk)

      case ('state%v')
         call outfld(state_snapshot(i)%standard_name, state%v, pcols, lchnk)

      case ('state%s')
         call outfld(state_snapshot(i)%standard_name, state%s, pcols, lchnk)

      case ('state%omega')
         call outfld(state_snapshot(i)%standard_name, state%omega, pcols, lchnk)

      case ('state%pmid')
         call outfld(state_snapshot(i)%standard_name, state%pmid, pcols, lchnk)

      case ('state%pmiddry')
         call outfld(state_snapshot(i)%standard_name, state%pmiddry, pcols, lchnk)

      case ('state%pdel')
         call outfld(state_snapshot(i)%standard_name, state%pdel, pcols, lchnk)

      case ('state%pdeldry')
         call outfld(state_snapshot(i)%standard_name, state%pdeldry, pcols, lchnk)

      case ('state%rpdel')
         call outfld(state_snapshot(i)%standard_name, state%rpdel, pcols, lchnk)

      case ('state%rpdeldry')
         call outfld(state_snapshot(i)%standard_name, state%rpdeldry, pcols, lchnk)

      case ('state%lnpmid')
         call outfld(state_snapshot(i)%standard_name, state%lnpmid, pcols, lchnk)

      case ('state%lnpmiddry')
         call outfld(state_snapshot(i)%standard_name, state%lnpmiddry, pcols, lchnk)

      case ('state%exner')
         call outfld(state_snapshot(i)%standard_name, state%exner, pcols, lchnk)

      case ('state%zm')
         call outfld(state_snapshot(i)%standard_name, state%zm, pcols, lchnk)

      case ('state%pint')
         call outfld(state_snapshot(i)%standard_name, state%pint, pcols, lchnk)

      case ('state%pintdry')
         call outfld(state_snapshot(i)%standard_name, state%pintdry, pcols, lchnk)

      case ('state%lnpint')
         call outfld(state_snapshot(i)%standard_name, state%lnpint, pcols, lchnk)

      case ('state%lnpintdry')
         call outfld(state_snapshot(i)%standard_name, state%lnpintdry, pcols, lchnk)

      case ('state%zi')
         call outfld(state_snapshot(i)%standard_name, state%zi, pcols, lchnk)

      case ('state%te_ini')
         call outfld(state_snapshot(i)%standard_name, state%te_ini, pcols, lchnk)

      case ('state%te_cur')
         call outfld(state_snapshot(i)%standard_name, state%te_cur, pcols, lchnk)

      case ('state%tw_ini')
         call outfld(state_snapshot(i)%standard_name, state%tw_ini, pcols, lchnk)

      case ('state%tw_cur')
         call outfld(state_snapshot(i)%standard_name, state%tw_cur, pcols, lchnk)

      case default
         call endrun('ERROR in state_snapshot_all_outfld: no match found for '//trim(state_snapshot(i)%ddt_string))

      end select

      call cam_history_snapshot_deactivate(trim(state_snapshot(i)%standard_name))

   end do

end subroutine state_snapshot_all_outfld

subroutine cam_snapshot_ptend_outfld(ptend, lchnk)

   use constituents, only: cnst_name, cnst_longname
   !--------------------------------------------------------
   ! This subroutine does the outfld calls for ptend fields.
   !--------------------------------------------------------

   type(physics_ptend), intent(in) :: ptend
   integer,             intent(in) :: lchnk

   integer                         :: mcnst
   character(len=128)              :: fname

   !--------------------------------------------------------
   ! Add the physics_ptend variables to the output
   !--------------------------------------------------------

   if (ptend%ls) then
      call outfld('ptend_s_snapshot', ptend%s, pcols, lchnk)

      call outfld('ptend_hflux_srf_snapshot', ptend%hflux_srf, pcols, lchnk)

      call outfld('ptend_hflux_top_snapshot', ptend%hflux_top, pcols, lchnk)
   end if

   if (ptend%lu) then
      call outfld('ptend_u_snapshot', ptend%u, pcols, lchnk)

      call outfld('ptend_taux_srf_snapshot', ptend%taux_srf, pcols, lchnk)

      call outfld('ptend_taux_top_snapshot', ptend%taux_top, pcols, lchnk)
   end if

   if (ptend%lv) then
      call outfld('ptend_v_snapshot', ptend%v, pcols, lchnk)

      call outfld('ptend_tauy_srf_snapshot', ptend%tauy_srf, pcols, lchnk)

      call outfld('ptend_tauy_top_snapshot', ptend%tauy_top, pcols, lchnk)
   end if

   do mcnst = 1, pcnst
      if (ptend%lq(mcnst)) then
         fname = 'ptend_'//trim(cnst_name(mcnst))//'_snapshot'
         call outfld(trim(fname), ptend%q(:,:,mcnst), pcols, lchnk)

         fname = 'cflx_srf_'//trim(cnst_name(mcnst))//'_snapshot'
         call outfld(trim(fname), ptend%cflx_srf(:,mcnst), pcols, lchnk)

         fname = 'cflx_top_'//trim(cnst_name(mcnst))//'_snapshot'
         call outfld(trim(fname), ptend%cflx_top(:,mcnst), pcols, lchnk)
      end if
   end do


end subroutine cam_snapshot_ptend_outfld

subroutine cnst_snapshot_all_outfld(lchnk, file_num, cnst)

   integer,             intent(in)  :: lchnk
   integer,             intent(in)  :: file_num
   real(r8),            intent(in)  :: cnst(:,:,:)

   integer :: i

   do i=1, ncnst_var

      ! Turn on the writing for only the requested tape (file_num)
      call cam_history_snapshot_activate(trim(cnst_snapshot(i)%standard_name), file_num)
      call outfld(cnst_snapshot(i)%standard_name, cnst(:,:,i), pcols, lchnk)

      ! Now that the field has been written, turn off the writing for field
      call cam_history_snapshot_deactivate(trim(cnst_snapshot(i)%standard_name))

   end do

end subroutine cnst_snapshot_all_outfld

subroutine tend_snapshot_all_outfld(lchnk, file_num, tend)

   integer,             intent(in)  :: lchnk
   integer,             intent(in)  :: file_num
   type(physics_tend),  intent(in)  :: tend

   integer :: i

   do i=1, ntend_var

      ! Turn on the writing for only the requested tape (file_num)
      call cam_history_snapshot_activate(trim(tend_snapshot(i)%standard_name), file_num)

      ! Select the tend field which is being written
      select case(tend_snapshot(i)%ddt_string)

      case ('tend%dtdt')
         call outfld(tend_snapshot(i)%standard_name, tend%dtdt, pcols, lchnk)

      case ('tend%dudt')
         call outfld(tend_snapshot(i)%standard_name, tend%dudt, pcols, lchnk)

      case ('tend%dvdt')
         call outfld(tend_snapshot(i)%standard_name, tend%dvdt, pcols, lchnk)

      case ('tend%flx_net')
         call outfld(tend_snapshot(i)%standard_name, tend%flx_net, pcols, lchnk)

      case ('tend%te_tnd')
         call outfld(tend_snapshot(i)%standard_name, tend%te_tnd, pcols, lchnk)

      case ('tend%tw_tnd')
         call outfld(tend_snapshot(i)%standard_name, tend%tw_tnd, pcols, lchnk)

      case default
         call endrun('ERROR in tend_snapshot_all_outfld: no match found for '//trim(tend_snapshot(i)%ddt_string))

      end select

      call cam_history_snapshot_deactivate(trim(tend_snapshot(i)%standard_name))

   end do

end subroutine tend_snapshot_all_outfld

subroutine cam_in_snapshot_all_outfld(lchnk, file_num, cam_in)

   integer,         intent(in)  :: lchnk
   integer,         intent(in)  :: file_num
   type(cam_in_t),  intent(in)  :: cam_in

   integer :: i

   do i=1, ncam_in_var

      ! Turn on the writing for only the requested tape (file_num)
      call cam_history_snapshot_activate(trim(cam_in_snapshot(i)%standard_name), file_num)

      ! Select the cam_in field which is being written
      select case(cam_in_snapshot(i)%ddt_string)

      case ('cam_in%landfrac')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%landfrac, pcols, lchnk)
      case ('cam_in%ocnfrac')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%ocnfrac, pcols, lchnk)
      case ('cam_in%snowhland')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%snowhland, pcols, lchnk)
      case ('cam_in%ts')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%ts, pcols, lchnk)
      case ('cam_in%sst')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%sst, pcols, lchnk)
      case ('cam_in%icefrac')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%icefrac, pcols, lchnk)
      case ('cam_in%shf')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%shf, pcols, lchnk)
      case ('cam_in%cflx')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%cflx, pcols, lchnk)
      case ('cam_in%wsx')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%wsx, pcols, lchnk)
      case ('cam_in%wsy')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%wsy, pcols, lchnk)
      case ('cam_in%asdif')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%asdif, pcols, lchnk)
      case ('cam_in%aldif')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%aldif, pcols, lchnk)
      case ('cam_in%lwup')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%lwup, pcols, lchnk)
      case ('cam_in%asdir')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%asdir, pcols, lchnk)
      case ('cam_in%aldir')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%aldir, pcols, lchnk)
      case ('cam_in%meganflx')
         if (associated (cam_in%meganflx)) &
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%meganflx, pcols, lchnk)
      case ('cam_in%fireflx')
         if (associated (cam_in%fireflx)) &
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%fireflx, pcols, lchnk)
      case ('cam_in%fireztop')
         if (associated (cam_in%fireztop)) &
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%fireztop, pcols, lchnk)
      case ('cam_in%depvel')
         if (associated (cam_in%depvel)) &
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%depvel, pcols, lchnk)
      case ('cam_in%lhf')
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%lhf, pcols, lchnk)
      case ('cam_in%fv')
         if (associated (cam_in%fv)) &
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%fv, pcols, lchnk)
      case ('cam_in%ram1')
         if (associated (cam_in%ram1)) &
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%ram1, pcols, lchnk)
      case ('cam_in%dstflx')
         if (associated (cam_in%dstflx)) &
         call outfld(cam_in_snapshot(i)%standard_name, cam_in%dstflx, pcols, lchnk)

      case default
         call endrun('ERROR in cam_in_snapshot_all_outfld: no match found for '//trim(cam_in_snapshot(i)%ddt_string))

      end select

      call cam_history_snapshot_deactivate(trim(cam_in_snapshot(i)%standard_name))

   end do

end subroutine cam_in_snapshot_all_outfld

subroutine cam_out_snapshot_all_outfld(lchnk, file_num, cam_out)

   integer,         intent(in)  :: lchnk
   integer,         intent(in)  :: file_num
   type(cam_out_t), intent(in)  :: cam_out

   integer :: i

   do i=1, ncam_out_var

      ! Turn on the writing for only the requested tape (file_num)
      call cam_history_snapshot_activate(trim(cam_out_snapshot(i)%standard_name), file_num)

      ! Select the cam_out field which is being written
      select case(cam_out_snapshot(i)%ddt_string)

      case ('cam_out%precc')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%precc, pcols, lchnk)

      case ('cam_out%precl')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%precl, pcols, lchnk)

      case ('cam_out%precsc')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%precsc, pcols, lchnk)

      case ('cam_out%precsl')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%precsl, pcols, lchnk)

      case ('cam_out%nhx_nitrogen_flx')
         if (associated(cam_out%nhx_nitrogen_flx)) &
            call outfld(cam_out_snapshot(i)%standard_name, cam_out%nhx_nitrogen_flx, pcols, lchnk)

      case ('cam_out%noy_nitrogen_flx')
         if (associated(cam_out%noy_nitrogen_flx)) &
            call outfld(cam_out_snapshot(i)%standard_name, cam_out%noy_nitrogen_flx, pcols, lchnk)

      case ('cam_out%bcphodry')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%bcphodry, pcols, lchnk)

      case ('cam_out%bcphidry')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%bcphidry, pcols, lchnk)

      case ('cam_out%ocphodry')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%ocphodry, pcols, lchnk)

      case ('cam_out%ocphidry')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%ocphidry, pcols, lchnk)

      case ('cam_out%bcphiwet')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%bcphiwet, pcols, lchnk)

      case ('cam_out%ocphiwet')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%ocphiwet, pcols, lchnk)

      case ('cam_out%dstwet1')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstwet1, pcols, lchnk)

      case ('cam_out%dstwet2')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstwet2, pcols, lchnk)

      case ('cam_out%dstwet3')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstwet3, pcols, lchnk)

      case ('cam_out%dstwet4')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstwet4, pcols, lchnk)

      case ('cam_out%dstdry1')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstdry1, pcols, lchnk)

      case ('cam_out%dstdry2')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstdry2, pcols, lchnk)

      case ('cam_out%dstdry3')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstdry3, pcols, lchnk)

      case ('cam_out%dstdry4')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%dstdry4, pcols, lchnk)

      case ('cam_out%sols')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%sols, pcols, lchnk)

      case ('cam_out%soll')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%soll, pcols, lchnk)

      case ('cam_out%solsd')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%solsd, pcols, lchnk)

      case ('cam_out%solld')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%solld, pcols, lchnk)

      case ('cam_out%flwds')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%flwds, pcols, lchnk)

      case ('cam_out%netsw')
         call outfld(cam_out_snapshot(i)%standard_name, cam_out%netsw, pcols, lchnk)

      case default
         call endrun('ERROR in cam_out_snapshot_all_outfld: no match found for '//trim(cam_out_snapshot(i)%ddt_string))

      end select

      call cam_history_snapshot_deactivate(trim(cam_out_snapshot(i)%standard_name))

   end do

end subroutine cam_out_snapshot_all_outfld

subroutine cam_pbuf_snapshot_all_outfld(lchnk, file_num, pbuf)
   use physics_buffer, only: pbuf_is_used

   integer,                   intent(in) :: lchnk
   integer,                            intent(in) :: file_num
   type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)

   integer :: i, pbuf_idx, ndims
   real(r8), pointer, dimension(:,:)           :: tmpptr2d
   real(r8), pointer, dimension(:,:,:)         :: tmpptr3d
   real(r8), pointer, dimension(:,:,:,:)       :: tmpptr4d
   real(r8), pointer, dimension(:,:,:,:,:)     :: tmpptr5d


   do i=1, npbuf_var

      pbuf_idx= pbuf_get_index(pbuf_snapshot(i)%ddt_string)

      if (pbuf_is_used(pbuf(pbuf_idx))) then
         ! Turn on the writing for only the requested tape (file_num)
         call cam_history_snapshot_activate(trim(pbuf_snapshot(i)%standard_name), file_num)

         ! Retrieve the pbuf data (dependent on the number of dimensions)
         ndims = count(pbuf_snapshot(i)%dim_name(:) /= '')

         select case (ndims)  ! Note that dimension 5 and 6 do not work with pbuf_get_field, so these are not used here

         case (1)
            call pbuf_get_field(pbuf, pbuf_idx, tmpptr2d)
            call outfld(pbuf_snapshot(i)%standard_name, tmpptr2d, pcols, lchnk)

         case (2)
            call pbuf_get_field(pbuf, pbuf_idx, tmpptr3d)
            call outfld(pbuf_snapshot(i)%standard_name, tmpptr3d, pcols, lchnk)

         case (3)
            call pbuf_get_field(pbuf, pbuf_idx, tmpptr3d)
            call outfld(pbuf_snapshot(i)%standard_name, tmpptr4d, pcols, lchnk)

         case (4)
            call pbuf_get_field(pbuf, pbuf_idx, tmpptr5d)
            call outfld(pbuf_snapshot(i)%standard_name, tmpptr5d, pcols, lchnk)

         end select

         ! Now that the field has been written, turn off the writing for field
         call cam_history_snapshot_deactivate(trim(pbuf_snapshot(i)%standard_name))


      end if

   end do

end subroutine cam_pbuf_snapshot_all_outfld

subroutine fill_pbuf_info(pbuf_info, pbuf, const_cname)

!---------------------------------------------------
! This subroutine exists to link the pbuf name with units.  It can be expanded to include standard_names
! at a later date if needed.  It is a list of all the pbuf fields that are called within CAM with actual
! names.
!---------------------------------------------------

     type(pbuf_info_type),   intent(inout) :: pbuf_info(:)
     type(physics_buffer_desc), intent(in) :: pbuf(:)
     character(len=*),          intent(in) :: const_cname(:)

     logical, dimension(size(pbuf)) ::  found
     character(len=24), dimension(2,npbuf_all) :: pbuf_all
     character(len=24) :: pbuf_name
     integer  :: i, ipbuf

     found(:) = .false.

     pbuf_all(1:2,1:100) = reshape ( (/  &
          'ACCRE_ENHAN            ','unset                  ',&
          'ACGCME                 ','unset                  ',&
          'ACLDY_CEN              ','unset                  ',&
          'ACNUM                  ','unset                  ',&
          'ACPRECL                ','unset                  ',&
          'AIST                   ','unset                  ',&
          'ALST                   ','unset                  ',&
          'am_evp_st              ','unset                  ',&
          'AMIE_efxg              ','mW/m2                  ',&
          'AMIE_kevg              ','keV                    ',&
          'AST                    ','1                      ',&
          'AurIPRateSum           ','unset                  ',&
          'awk_PBL                ','unset                  ',&
          'bprod                  ','unset                  ',&
          'cam3_bcphi             ','unset                  ',&
          'cam3_bcpho             ','unset                  ',&
          'cam3_dust1             ','unset                  ',&
          'cam3_dust2             ','unset                  ',&
          'cam3_dust3             ','unset                  ',&
          'cam3_dust4             ','unset                  ',&
          'cam3_ocphi             ','unset                  ',&
          'cam3_ocpho             ','unset                  ',&
          'cam3_ssam              ','unset                  ',&
          'cam3_sscm              ','unset                  ',&
          'cam3_sul               ','unset                  ',&
          'CC_ni                  ','unset                  ',&
          'CC_nl                  ','unset                  ',&
          'CC_qi                  ','unset                  ',&
          'CC_ql                  ','unset                  ',&
          'CC_qlst                ','unset                  ',&
          'CC_qv                  ','unset                  ',&
          'CC_T                   ','unset                  ',&
          'CICEWP                 ','unset                  ',&
          'CLDBOT                 ','1                      ',&
          'CLDEMIS                ','unset                  ',&
          'CLDFGRAU               ','1                      ',&
          'CLDFSNOW               ','1                      ',&
          'CLD                    ','unset                  ',&
          'CLDICEINI              ','unset                  ',&
          'CLDLIQINI              ','unset                  ',&
          'CLDO                   ','unset                  ',&
          'CLDTAU                 ','unset                  ',&
          'CLDTOP                 ','1                      ',&
          'CLIQWP                 ','unset                  ',&
          'CLOUD_FRAC             ','unset                  ',&
          'CLUBB_BUFFER           ','unset                  ',&
          'CMELIQ                 ','kg/kg/s                ',&
          'CMFMC_SH               ','unset                  ',&
          'cmfr_det               ','kg/m2/s                ',&
          'CONCLD                 ','fraction               ',&
          'CRM_CLD_RAD            ','unset                  ',&
          'CRM_DGNUMWET           ','unset                  ',&
          'CRM_NC                 ','/kg                    ',&
          'CRM_NC_RAD             ','unset                  ',&
          'CRM_NG                 ','/kg                    ',&
          'CRM_NI                 ','/kg                    ',&
          'CRM_NI_RAD             ','unset                  ',&
          'CRM_NR                 ','/kg                    ',&
          'CRM_NS                 ','/kg                    ',&
          'CRM_NS_RAD             ','unset                  ',&
          'CRM_QAERWAT            ','unset                  ',&
          'CRM_QC                 ','kg/kg                  ',&
          'CRM_QC_RAD             ','unset                  ',&
          'CRM_QG                 ','kg/kg                  ',&
          'CRM_QI                 ','kg/kg                  ',&
          'CRM_QI_RAD             ','unset                  ',&
          'CRM_QN                 ','unset                  ',&
          'CRM_QP                 ','kg/kg                  ',&
          'CRM_QRAD               ','unset                  ',&
          'CRM_QR                 ','kg/kg                  ',&
          'CRM_QS                 ','kg/kg                  ',&
          'CRM_QS_RAD             ','unset                  ',&
          'CRM_QT                 ','unset                  ',&
          'CRM_QV_RAD             ','unset                  ',&
          'CRM_T                  ',' K                     ',&
          'CRM_T_RAD              ','unset                  ',&
          'CRM_U                  ','m/s                    ',&
          'CRM_V                  ','m/s                    ',&
          'CRM_W                  ','m/s                    ',&
          'CT                     ','unset                  ',&
          'cu_cmfr                ','kg/m2/s                ',&
          'cuorg                  ','unset                  ',&
          'cu_qir                 ','kg/kg                  ',&
          'cu_qlr                 ','kg/kg                  ',&
          'cu_qtr                 ','kg/kg                  ',&
          'cushavg                ','m                      ',&
          'cush                   ','m                      ',&
          'cu_thlr                ','K                      ',&
          'cu_trr                 ','unset                  ',&
          'cu_ur                  ','m/s                    ',&
          'cu_vr                  ','m/s                    ',&
          'CV_REFFICE             ','micron                 ',&
          'CV_REFFLIQ             ','micron                 ',&
          'DEGRAU                 ','unset                  ',&
          'DEI                    ','unset                  ',&
          'delta_qt_PBL           ','unset                  ',&
          'delta_thl_PBL          ','unset                  ',&
          'delta_tr_PBL           ','unset                  ',&
          'delta_u_PBL            ','unset                  ',&
          'delta_v_PBL            ','unset                  '/) ,  (/2,100/))

     pbuf_all(1:2,101:200) = reshape ( (/  &
          'DES                    ','unset                  ',&
          'DGNUM                  ','unset                  ',&
          'DGNUMWET               ','unset                  ',&
          'DIFZM                  ','kg/kg/s                ',&
          'DLFZM                  ','kg/kg/s                ',&
          'DNIFZM                 ','1/kg/s                 ',&
          'DNLFZM                 ','1/kg/s                 ',&
          'DP_CLDICE              ','unset                  ',&
          'DP_CLDLIQ              ','unset                  ',&
          'DP_FLXPRC              ','unset                  ',&
          'DP_FLXSNW              ','unset                  ',&
          'DP_FRAC                ','unset                  ',&
          'dragblj                ','1/s                    ',&
          'DRYMASS                ','unset                  ',&
          'DRYRAD                 ','unset                  ',&
          'DRYVOL                 ','unset                  ',&
          'DTCORE                 ','K/s                    ',&
          'evprain_st             ','unset                  ',&
          'evpsnow_st             ','unset                  ',&
          'FICE                   ','fraction               ',&
          'FLNS                   ','W/m2                   ',&
          'FLNT                   ','W/m2                   ',&
          'FRACIS                 ','unset                  ',&
          'FRACSOA                ','unset                  ',&
          'FRACSOG                ','unset                  ',&
          'FRONTGA                ','unset                  ',&
          'FRONTGF                ','K^2/M^2/S              ',&
          'FRZCNT                 ','unset                  ',&
          'FRZDEP                 ','unset                  ',&
          'FRZIMM                 ','unset                  ',&
          'FSDS                   ','W/m2                   ',&
          'FSNS                   ','W/m2                   ',&
          'FSNT                   ','W/m2                   ',&
          'HallConduct            ','unset                  ',&
          'HYGRO                  ','unset                  ',&
          'ICCWAT                 ','unset                  ',&
          'ICGRAUWP               ','unset                  ',&
          'ICIWP                  ','unset                  ',&
          'ICIWPST                ','unset                  ',&
          'ICLWP                  ','unset                  ',&
          'ICLWPST                ','unset                  ',&
          'ICSWP                  ','unset                  ',&
          'ICWMRDP                ','kg/kg                  ',&
          'ICWMRSH                ','kg/kg                  ',&
          'IonRates               ','unset                  ',&
          'ipbl                   ','unset                  ',&
          'ISS_FRAC               ','unset                  ',&
          'kpblh                  ','unset                  ',&
          'ksrftms                ','unset                  ',&
          'kvh                    ','m2/s                   ',&
          'kvm                    ','m2/s                   ',&
          'kvt                    ','m2/s                   ',&
          'LAMBDAC                ','unset                  ',&
          'LANDM                  ','unset                  ',&
          'LCWAT                  ','unset                  ',&
          'LD                     ','unset                  ',&
          'LHFLX                  ','W/m2                   ',&
          'LHFLX_RES              ','unset                  ',&
          'LS_FLXPRC              ','kg/m2/s                ',&
          'LS_FLXSNW              ','kg/m2/s                ',&
          'LS_MRPRC               ','unset                  ',&
          'LS_MRSNW               ','unset                  ',&
          'LS_REFFRAIN            ','micron                 ',&
          'LS_REFFSNOW            ','micron                 ',&
          'LU                     ','unset                  ',&
          'MAMH2SO4EQ             ','unset                  ',&
          'MU                     ','Pa/s                   ',&
          'NAAI_HOM               ','unset                  ',&
          'NAAI                   ','unset                  ',&
          'NACON                  ','unset                  ',&
          'NAER                   ','unset                  ',&
          'NEVAPR_DPCU            ','unset                  ',&
          'NEVAPR                 ','unset                  ',&
          'NEVAPR_SHCU            ','unset                  ',&
          'NIWAT                  ','unset                  ',&
          'NLWAT                  ','unset                  ',&
          'NMXRGN                 ','unset                  ',&
          'NPCCN                  ','unset                  ',&
          'NRAIN                  ','m-3                    ',&
          'NSNOW                  ','m-3                    ',&
          'O3                     ','unset                  ',&
          'pblh                   ','m                      ',&
          'PDF_PARAMS             ','unset                  ',&
          'PDF_PARAMS_ZM          ','unset                  ',&
          'PedConduct             ','unset                  ',&
          'PMXRGN                 ','unset                  ',&
          'PRAIN                  ','unset                  ',&
          'PREC_DP                ','unset                  ',&
          'PREC_PCW               ','m/s                    ',&
          'PREC_SED               ','unset                  ',&
          'PREC_SH                ','unset                  ',&
          'PREC_SH                ','unset                  ',&
          'PREC_STR               ','unset                  ',&
          'PRER_EVAP              ','unset                  ',&
          'PSL                    ','Pa                     ',&
          'QAERWAT                ','unset                  ',&
          'QCWAT                  ','unset                  ',&
          'QFLX                   ','kg/m2/s                ',&
          'QFLX_RES               ','unset                  ',&
          'QINI                   ','unset                  '  /),    (/2,100/))

     pbuf_all(1:2,201:300) = reshape ( (/  &
          'qir_det                ','kg/kg                  ',&
          'QIST                   ','unset                  ',&
          'qlr_det                ','kg/kg                  ',&
          'QLST                   ','unset                  ',&
          'QME                    ','unset                  ',&
          'qpert                  ','kg/kg                  ',&
          'QRAIN                  ','kg/kg                  ',&
          'QRL                    ','K/s                    ',&
          'qrlin                  ','unset                  ',&
          'QRS                    ','K/s                    ',&
          'qrsin                  ','unset                  ',&
          'QSATFAC                ','-                      ',&
          'QSNOW                  ','kg/kg                  ',&
          'QTeAur                 ','unset                  ',&
          'qti_flx                ','unset                  ',&
          'qtl_flx                ','unset                  ',&
          'RAD_CLUBB              ','unset                  ',&
          'RATE1_CW2PR_ST         ','unset                  ',&
          'RCM                    ','unset                  ',&
          'RE_ICE                 ','unset                  ',&
          'REI                    ','micron                 ',&
          'RELHUM                 ','percent                ',&
          'REL                    ','micron                 ',&
          'RELVAR                 ','-                      ',&
          'RNDST                  ','unset                  ',&
          'RPRDDP                 ','unset                  ',&
          'RPRDSH                 ','unset                  ',&
          'RPRDTOT                ','unset                  ',&
          'RTM                    ','unset                  ',&
          'rtp2_mc_zt             ','unset                  ',&
          'RTP2_nadv              ','unset                  ',&
          'rtpthlp_mc_zt          ','unset                  ',&
          'RTPTHLP_nadv           ','unset                  ',&
          'RTPTHVP                ','unset                  ',&
          'SADICE                 ','cm2/cm3                ',&
          'SADSNOW                ','cm2/cm3                ',&
          'SADSULF                ','unset                  ',&
          'SD                     ','unset                  ',&
          'SGH30                  ','unset                  ',&
          'SGH                    ','unset                  ',&
          'SH_CLDICE1             ','unset                  ',&
          'SH_CLDICE              ','unset                  ',&
          'SH_CLDLIQ1             ','unset                  ',&
          'SH_CLDLIQ              ','unset                  ',&
          'SH_E_ED_RATIO          ','unset                  ',&
          'SHFLX                  ','W/m2                   ',&
          'SH_FLXPRC              ','unset                  ',&
          'SHFLX_RES              ','unset                  ',&
          'SH_FLXSNW              ','unset                  ',&
          'SH_FRAC                ','unset                  ',&
          'shfrc                  ','unset                  ',&
          'smaw                   ','unset                  ',&
          'SNOW_DP                ','unset                  ',&
          'SNOW_PCW               ','unset                  ',&
          'SNOW_SED               ','unset                  ',&
          'SNOW_SH                ','unset                  ',&
          'SNOW_STR               ','unset                  ',&
          'SO4DRYVOL              ','unset                  ',&
          'SSLTA                  ','kg/kg                  ',&
          'SSLTC                  ','kg/kg                  ',&
          'SU                     ','unset                  ',&
          "taubljx                ",'N/m2                   ',&
          "taubljy                ",'N/m2                   ',&
          'tauresx                ','unset                  ',&
          'tauresy                ','unset                  ',&
          "tautmsx                ",'N/m2                   ',&
          "tautmsy                ",'N/m2                   ',&
          'TAUX                   ','N/m2                   ',&
          'TAUX_RES               ','unset                  ',&
          'TAUY                   ','N/m2                   ',&
          'TAUY_RES               ','unset                  ',&
          'tcorr                  ','unset                  ',&
          'TCWAT                  ','unset                  ',&
          'TElec                  ','K                      ',&
          'TEOUT                  ','J/m2                   ',&
          'THLM                   ','unset                  ',&
          'thlp2_mc_zt            ','unset                  ',&
          'THLP2_nadv             ','unset                  ',&
          'THLPTHVP               ','unset                  ',&
          'TIon                   ','K                      ',&
          'TK_CRM                 ','unset                  ',&
          'tke                    ','m2/s2                  ',&
          'tkes                   ','m2/s2                  ',&
          'TND_NSNOW              ','unset                  ',&
          'TND_QSNOW              ','unset                  ',&
          'tpert                  ','K                      ',&
          'TREFMNAV               ','K                      ',&
          'TREFMXAV               ','K                      ',&
          'tropp                  ','unset                  ',&
          'TSTCPY_SCOL            ','unset                  ',&
          'TTEND_DP               ','unset                  ',&
          'TTEND_SH               ','unset                  ',&
          'T_TTEND                ','unset                  ',&
          'turbtype               ','unset                  ',&
          "UI                     ",'m/s                    ',&
          'UM                     ','unset                  ',&
          'UP2_nadv               ','unset                  ',&
          'UPWP                   ','m^2/s^2                ',&
          'UZM                    ','M/S                    ',&
          'VI                     ','m/s                    '    /),                  (/2,100/))

     pbuf_all(1:2,301:npbuf_all) = reshape ( (/  &
          'VM                     ','m/s                    ',&
          'VOLC_MMR               ','unset                  ',&
          'VOLC_RAD_GEOM          ','unset                  ',&
          'VP2_nadv               ','unset                  ',&
          'VPWP                   ','m^2/s^2                ',&
          'went                   ','m/s                    ',&
          'WETDENS_AP             ','unset                  ',&
          "WI                     ",'m/s                    ',&
          'WP3_nadv               ','unset                  ',&
          'wprtp_mc_zt            ','unset                  ',&
          'WPRTP_nadv             ','unset                  ',&
          'wpthlp_mc_zt           ','unset                  ',&
          'WPTHLP_nadv            ','unset                  ',&
          'WPTHVP                 ','unset                  ',&
          'WSEDL                  ','unset                  ',&
          'wstarPBL               ','unset                  ',&
          'ZM_DP                  ','unset                  ',&
          'ZM_DSUBCLD             ','unset                  ',&
          'ZM_DU                  ','unset                  ',&
          'ZM_ED                  ','unset                  ',&
          'ZM_EU                  ','unset                  ',&
          'ZM_IDEEP               ','unset                  ',&
          'ZM_JT                  ','unset                  ',&
          'ZM_MAXG                ','unset                  ',&
          'ZM_MD                  ','unset                  ',&
          'ZM_MU                  ','unset                  ',&
          'ZTODT                  ','unset                  '  /),                     (/2,27/))

! Fields which are added with pbuf_add_field calls, but are data driven.  These are not
! included in the above list.  This means that these fields will not have proper units
! set for them
!          'CG' // shortname,        'unset',        &
!          'CI' // shortname,        'unset',        &
!          'CL' // shortname,        'unset',        &
!          ghg_names(i),             'unset',        &
!          mmr_name1,                'unset',        &
!          mmr_name2,                'unset',        &
!          mmr_name3,                'unset',        &
!          mmr_name,                 'unset',        &
!          ozone_name,               'unset',        &
!          pbufname,                 'unset',        &
!          pbufname,                 'unset',        &
!          pbuf_names(i),            'unset',        &
!          rad_name1,                'unset',        &
!          rad_name2,                'unset',        &
!          rad_name3,                'unset',        &
!          rad_name,                 'unset',        &
!          sad_name,                 'cm2/cm3',      &
!          volcaero_name,            'kg/kg',        &
!          volcrad_name,             'm',            &
!          xname_massptrcw(l,        'unset',        &
!          xname_numptrcw,           'unset',        &
!          aero_names(mm)
!          cnst_names(iconst)

   do ipbuf = 1, size(pbuf)
     pbuf_name = pbuf_get_field_name(ipbuf)
     i = 1
     do while ((i <= npbuf_all) .and. .not. found(ipbuf))
        if (trim(pbuf_all(1,i)) == trim(pbuf_name)) then
           pbuf_info(ipbuf)%name          = trim(pbuf_all(1,i))
           pbuf_info(ipbuf)%standard_name = trim(pbuf_all(1,i))//'_snapshot'
           pbuf_info(ipbuf)%units         = trim(pbuf_all(2,i))
           pbuf_info(ipbuf)%dim_string(:) = ' '
           found(ipbuf) = .true.
       end if
        i = i+1
      end do
      if (.not. found(ipbuf)) then

         i = 1
         ! Check if variable is a variation of constituent - then use the same units
         do while ((i <= ncnst_var) .and. .not. found(ipbuf))
            if (trim(const_cname(i)) == trim(pbuf_name)) then
               pbuf_info(ipbuf) = pbuf_info_type(trim(const_cname(i)),trim(const_cname(i))//'_snapshot',&
                                                 trim(cnst_snapshot(i)%units),      ' ')
               found(ipbuf) = .true.
            end if
            i = i+1
         end do
      end if

        ! Found a pbuf that has not been added to this routine
      if (.not. found(ipbuf)) then
         write(iulog,*) 'WARNING - no units information for: '//trim(pbuf_name)

         pbuf_info(ipbuf)%name          = trim(pbuf_name)
         pbuf_info(ipbuf)%standard_name = trim(pbuf_name)//'_snapshot'
         pbuf_info(ipbuf)%units         = 'unset'
         pbuf_info(ipbuf)%dim_string(:) = ' '
         found(ipbuf) = .true.
      end if

   end do

end subroutine fill_pbuf_info


end module cam_snapshot
