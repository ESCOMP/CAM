module cam_snapshot
!--------------------------------------------------------
! The purpose of this module is to handle taking the "snapshot" of CAM data.
!
! This module writes out ALL the state, tend and pbuf fields.  It also includes the cam_in and cam_out
! fields which are used within CAM
!--------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
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

type snapshot_type

  character(len=40)  :: ddt_string
  character(len=256) :: standard_name
  character(len=20)  :: dim_name
  character(len=8)   :: units

end type snapshot_type

integer :: nstate_var 
integer :: ncnst_var 
integer :: ntend_var 
integer :: ncam_out_var

! Note the maximum number of variables for each type
type (snapshot_type) ::  state_snapshot(27)
type (snapshot_type) ::  cnst_snapshot(40)
type (snapshot_type) ::  tend_snapshot(6)
type (snapshot_type) ::  cam_out_snapshot(30)

contains

subroutine cam_snapshot_init()

use cam_history, only : cam_history_snapshot_nhtfrq_set

!--------------------------------------------------------
! This subroutine does the addfld calls for ALL state, tend and pbuf fields.  It also includes the cam_in and cam_out
! elements which are used within CAM
!--------------------------------------------------------

   integer  :: cam_snapshot_before_num, cam_snapshot_after_num, cam_snapshot_nhtfrq

   call phys_getopts(cam_snapshot_before_num_out = cam_snapshot_before_num, &
                     cam_snapshot_after_num_out  = cam_snapshot_after_num, &
                     cam_snapshot_nhtfrq_out     = cam_snapshot_nhtfrq)

   if (cam_snapshot_before_num <= 0 .and. cam_snapshot_after_num <= 0) return ! No snapshot files are being requested

   call cam_history_snapshot_nhtfrq_set (cam_snapshot_before_num, cam_snapshot_after_num, cam_snapshot_nhtfrq)

   call cam_state_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)
   call cam_cnst_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)
   call cam_tend_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)
   call cam_out_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

end subroutine cam_snapshot_init

subroutine cam_snapshot_all_outfld(file_num, state, tend, cam_out)

!--------------------------------------------------------
! This subroutine does the outfld calls for ALL state, tend and pbuf fields.  It also includes the cam_in and cam_out
! elements which are used within CAM
!--------------------------------------------------------

   integer,             intent(in) :: file_num
   type(physics_state), intent(in) :: state
   type(physics_tend),  intent(in) :: tend
   type(cam_out_t),     intent(in) :: cam_out

   integer :: i
   integer :: lchnk

   lchnk = state%lchnk

   ! Write out all the state fields
   call state_snapshot_all_outfld(lchnk, file_num, state)

   ! Write out all the constituent fields
   call cnst_snapshot_all_outfld(lchnk, file_num, state%q)

   ! Write out all the tendency fields
   call tend_snapshot_all_outfld(lchnk, file_num, tend)

   ! Write out all the cam_out fields
   call cam_out_snapshot_all_outfld(lchnk, file_num, cam_out)

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

   do i=1,nstate_var
      call cam_history_snapshot_deactivate(state_snapshot(i)%standard_name)
   end do

end subroutine cam_snapshot_deactivate


subroutine cam_state_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for state, tend and pbuf fields.
!--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num
   integer :: i
    
   nstate_var = 0

   !--------------------------------------------------------
   ! Add the state variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%ps',        'ps_snapshot',         'Pa',    'horiz_only')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%psdry',     'psdry_snapshot',      'Pa',    'horiz_only')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%phis',      'phis_snapshot',       'm2/m2', 'horiz_only')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%t',         't_snapshot',          'K',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%u',         'u_snapshot',          'm s-1',   'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%v',         'v_snapshot',          'm s-1',   'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%s',         's_snapshot',          ' ',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%omega',     'omega_snapshot',      'Pa s-1',  'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pmid',      'pmid_snapshot',       'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pmiddry',   'pmiddry_snapshot',    'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pdel',      'pdel_snapshot',       'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pdeldry',   'pdeldry_snapshot',    'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%rpdel',     'rpdel_snapshot',      'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%rpdeldry',  'rpdeldry_snapshot',   'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpmid',    'lnpmid_snapshot',     ' ',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpmiddry', 'lnpmiddry_snapshot',  ' ',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%exner',     'exner_snapshot',      ' ',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%zm',        'zm_snapshot',         'm',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pint',      'pint_snapshot',       'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%pintdry',   'pintdry_snapshot',    'Pa',    'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpint',    'lnpint_snapshot',     ' ',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%lnpintdry', 'lnpintdry_snapshot',  ' ',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%zi',        'zi_snapshot',         'm',     'lev')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%te_ini',    'te_ini_snapshot',     ' ',     'horiz_only')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%te_cur',    'te_cur_snapshot',     ' ',     'horiz_only')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%tw_ini',    'tw_ini_snapshot',     ' ',     'horiz_only')

   call snapshot_addfld( nstate_var, state_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'state%tw_cur',    'tw_cur_snapshot',     ' ',     'horiz_only')

end subroutine cam_state_snapshot_init

subroutine cam_cnst_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for state, tend and pbuf fields., hence the writing
!        needs to be turned off for all fields, and will be turned on individaully
!        when needed.
!--------------------------------------------------------

   use constituents, only: pcnst, cnst_name

   integer, intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num

   integer :: i
    
   !--------------------------------------------------------
   ! Add the cnst variables to the output
   !--------------------------------------------------------

   ncnst_var = 0 ! Updated inside snapshot_addfld

   do while (ncnst_var < pcnst)
      call snapshot_addfld( ncnst_var, cnst_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
        'cnst_name(ncnst_var+1)',        trim(cnst_name(ncnst_var+1))//'_snapshot',         'kg kg-1',    'lev')
   end do

end subroutine cam_cnst_snapshot_init

subroutine cam_tend_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for state, tend and pbuf fields.
!--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num
   integer :: i
    
   ntend_var = 0

   !--------------------------------------------------------
   ! Add the state variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%dtdt',        'dtdt_snapshot',         'K s-1',    'lev')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%dudt',        'dudt_snapshot',         '',    'lev')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%dvdt',        'dvdt_snapshot',         '',    'lev')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%flx_net',        'flx_net_snapshot',   '',    'horiz_only')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%te_tnd',        'te_tnd_snapshot',   '',    'horiz_only')

   call snapshot_addfld( ntend_var, tend_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'tend%tw_tnd',        'tw_tnd_snapshot',   '',    'horiz_only')

end subroutine cam_tend_snapshot_init

subroutine cam_out_snapshot_init(cam_snapshot_before_num, cam_snapshot_after_num)

!--------------------------------------------------------
! This subroutine does the addfld calls for state, tend and pbuf fields.
!--------------------------------------------------------

   integer,intent(in) :: cam_snapshot_before_num, cam_snapshot_after_num
   integer :: i
    
   ncam_out_var = 0

   !--------------------------------------------------------
   ! Add the state variables to the output
   !--------------------------------------------------------

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precc',        'precc_snapshot',         'm s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precl',        'precl_snapshot',         'm s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precsc',        'precsc_snapshot',         'm s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%precsl',        'precsl_snapshot',         'm s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%nhx_nitrogen_flx',        'nhx_nitro_flx_snapshot',         'kgN m2-1 sec-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%noy_nitrogen_flx',        'noy_nitro_flx_snapshot',         'kgN m2-1 sec-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%bcphodry',        'bcphodry_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%bcphidry',        'bcphidry_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%ocphodry',        'ocphodry_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%ocphidry',        'ocphidry_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%bcphiwet',        'bcphiwet_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%ocphiwet',        'ocphiwet_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet1',        'dstwet1_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet2',        'dstwet2_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet3',        'dstwet3_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstwet4',        'dstwet4_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry1',        'dstdry1_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry2',        'dstdry2_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry3',        'dstdry3_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%dstdry4',        'dstdry4_snapshot',         'kg m-2 s-1',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%sols',        'sols_snapshot',         'W m-2',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%soll',        'soll_snapshot',         'W m-2',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%solsd',        'solsd_snapshot',         'W m-2',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%solld',        'solld_snapshot',         'W m-2',    'horiz_only')

   call snapshot_addfld( ncam_out_var, cam_out_snapshot,  cam_snapshot_before_num, cam_snapshot_after_num, &
     'cam_out%netsw',        'netsw_snapshot',         '     ',    'horiz_only')

end subroutine cam_out_snapshot_init

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

subroutine cnst_snapshot_all_outfld(lchnk, file_num, cnst)

   integer,             intent(in)  :: lchnk
   integer,             intent(in)  :: file_num
   real(r8),            intent(in)  :: cnst(:,:,:)

   integer :: i
   integer :: ff
   logical :: found

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
   integer :: ff
   logical :: found

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

subroutine cam_out_snapshot_all_outfld(lchnk, file_num, cam_out)

   integer,             intent(in)  :: lchnk
   integer,             intent(in)  :: file_num
   type(cam_out_t),     intent(in)  :: cam_out

   integer :: i
   integer :: ff
   logical :: found

   do i=1, ncam_out_var

      ! Turn on the writing for only the requested tape (file_num)
      call cam_history_snapshot_activate(trim(state_snapshot(i)%standard_name), file_num)

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

      call cam_history_snapshot_deactivate(trim(state_snapshot(i)%standard_name))

   end do

end subroutine cam_out_snapshot_all_outfld

end module cam_snapshot
