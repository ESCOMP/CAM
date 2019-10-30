module cam_snapshot
!--------------------------------------------------------
! The purpose of this module is to handle taking the "snapshot" of CAM data.
!
! This module writes out ALL the state, tend and pbuf fields.  It also includes the cam_in and cam_out
! fields which are used within CAM
!--------------------------------------------------------

use cam_history,    only: addfld, add_default, outfld, cam_history_snapshot_deactivate, cam_history_snapshot_activate
use cam_abortutils, only: endrun
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
use physics_types,  only: physics_state, physics_tend
use camsrfexch,     only: cam_out_t, cam_in_t
use ppgrid,         only: pcols
use phys_control,   only: phys_getopts

implicit  none

private

public :: cam_snapshot_init, cam_snapshot_all_outfld, cam_snapshot_deactivate

type state_snapshot_type

  character(len=40)  :: state_string
  character(len=256) :: standard_name
  character(len=20)  :: dim_name
  character(len=8)   :: units

end type state_snapshot_type

integer :: nstate_var 
integer, parameter :: max_nstate_var = 27

type (state_snapshot_type) ::  state_snapshot(max_nstate_var)

contains

subroutine cam_snapshot_init()

!--------------------------------------------------------
! This subroutine does the addfld calls for ALL state, tend and pbuf fields.  It also includes the cam_in and cam_out
! elements which are used within CAM
!--------------------------------------------------------

   integer :: i
   integer :: cam_snapshot_before_num, cam_snapshot_after_num
    
   nstate_var = 0

   call phys_getopts(cam_snapshot_before_num_out=cam_snapshot_before_num, &
                     cam_snapshot_after_num_out=cam_snapshot_after_num)

   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   !--------------------------------------------------------
   ! Add the state variables to the output
   !--------------------------------------------------------

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%ps',        'ps_snapshot',         'Pa',    'horiz_only')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%psdry',     'psdry_snapshot',      'Pa',    'horiz_only')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%phis',      'phis_snapshot',       'm2/m2', 'horiz_only')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%t',         't_snapshot',          'K',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%u',         'u_snapshot',          'm/s',   'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%v',         'v_snapshot',          'm/s',   'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%s',         's_snapshot',          ' ',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%omega',     'omega_snapshot',      'Pa/s',  'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pmid',      'pmid_snapshot',       'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pmiddry',   'pmiddry_snapshot',    'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pdel',      'pdel_snapshot',       'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pdeldry',   'pdeldry_snapshot',    'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%rpdel',     'rpdel_snapshot',      'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%rpdeldry',  'rpdeldry_snapshot',   'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpmid',    'lnpmid_snapshot',     ' ',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpmiddry', 'lnpmiddry_snapshot',  ' ',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%exner',     'exner_snapshot',      ' ',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%zm',        'zm_snapshot',         'm',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pint',      'pint_snapshot',       'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pintdry',   'pintdry_snapshot',    'Pa',    'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpint',    'lnpint_snapshot',     ' ',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpintdry', 'lnpintdry_snapshot',  ' ',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%zi',        'zi_snapshot',         'm',     'lev')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%te_ini',    'te_ini_snapshot',     ' ',     'horiz_only')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%te_cur',    'te_cur_snapshot',     ' ',     'horiz_only')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%tw_ini',    'tw_ini_snapshot',     ' ',     'horiz_only')

   call snapshot_state_addfld( cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%tw_cur',    'tw_cur_snapshot',     ' ',     'horiz_only')

end subroutine cam_snapshot_init

subroutine snapshot_state_addfld(cam_snapshot_before_num, cam_snapshot_after_num,&
    state_string, standard_name, units, dimension_string)
 
   integer,          intent(in) :: cam_snapshot_before_num
   integer,          intent(in) :: cam_snapshot_after_num
   character(len=*), intent(in) :: state_string
   character(len=*), intent(in) :: standard_name
   character(len=*), intent(in) :: units
   character(len=*), intent(in) :: dimension_string
  

   nstate_var=nstate_var+1

   if (nstate_var > size(state_snapshot)) &
      call endrun(' ERROR in cam_snapshot_addfld: state_snapshot array not allocated large enough')

   call addfld(standard_name, dimension_string, 'I', units, standard_name)
   if (cam_snapshot_before_num > 0) call add_default(standard_name, cam_snapshot_before_num, ' ')
   if (cam_snapshot_after_num > 0)  call add_default(standard_name, cam_snapshot_after_num, ' ')

   state_snapshot(nstate_var)%state_string  = state_string
   state_snapshot(nstate_var)%standard_name = standard_name
   state_snapshot(nstate_var)%dim_name      = dimension_string
   state_snapshot(nstate_var)%units         = units


end subroutine snapshot_state_addfld

subroutine cam_snapshot_all_outfld(file_num, state)


!--------------------------------------------------------
! This subroutine does the outfld calls for ALL state, tend and pbuf fields.  It also includes the cam_in and cam_out
! elements which are used within CAM
!--------------------------------------------------------

   integer,             intent(in)  :: file_num
   type(physics_state), intent(in) :: state

   integer :: i
   integer :: lchnk

   lchnk = state%lchnk

   !--------------------------------------------------------
   ! Write out all the state fields
   !--------------------------------------------------------

   call state_snapshot_all_outfld(lchnk, file_num, state)

end subroutine cam_snapshot_all_outfld


subroutine state_snapshot_all_outfld(lchnk, file_num, state)

   integer,             intent(in)  :: lchnk
   integer,             intent(in)  :: file_num
   type(physics_state), intent(in)  :: state

   integer :: i
   integer :: ff
   logical :: found

   do i=1, nstate_var

      ! Turn on the writing for only the requested tape (file_num)
      call cam_history_snapshot_activate(trim(state_snapshot(i)%standard_name), file_num)

      ! Select the state field which is being written
      select case(state_snapshot(i)%state_string)

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
         call endrun('ERROR in state_snapshot_all_outfld: no match found for '//trim(state_snapshot(i)%state_string))
      
      end select

      call cam_history_snapshot_deactivate(trim(state_snapshot(i)%standard_name))

   end do

end subroutine state_snapshot_all_outfld
    
subroutine cam_snapshot_deactivate()

   integer :: i

   do i=1,nstate_var
      call cam_history_snapshot_deactivate(state_snapshot(i)%standard_name)
   end do

end subroutine cam_snapshot_deactivate

end module cam_snapshot
