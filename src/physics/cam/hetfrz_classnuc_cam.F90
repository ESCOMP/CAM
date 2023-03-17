module hetfrz_classnuc_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for hetfrz_classnuc module.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: rair, cpair, rh2o, rhoh2o, mwh2o, tmelt, pi
use constituents,   only: cnst_get_ind
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: use_hetfrz_classnuc
use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, add_default, outfld, fieldname_len
use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: svp_water_vect, svp_ice_vect
use cam_logfile,    only: iulog
use error_messages, only: handle_errmsg, alloc_err
use cam_abortutils, only: endrun
use string_utils,   only: int2str
use hetfrz_classnuc,only: hetfrz_classnuc_init, hetfrz_classnuc_calc

use aerosol_properties_mod, only: aerosol_properties, aero_name_len
use aerosol_state_mod, only: aerosol_state

implicit none
private
save

public :: &
   hetfrz_classnuc_cam_readnl,   &
   hetfrz_classnuc_cam_register, &
   hetfrz_classnuc_cam_init,     &
   hetfrz_classnuc_cam_calc

! Namelist variables
logical :: hist_hetfrz_classnuc = .false.
real(r8) :: hetfrz_bc_scalfac = -huge(1._r8) ! scaling factor for BC
real(r8) :: hetfrz_dust_scalfac = -huge(1._r8) ! scaling factor for dust

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numliq_idx = -1, &
   numice_idx = -1

! pbuf indices for fields provided by heterogeneous freezing
integer :: &
   frzimm_idx = -1, &
   frzcnt_idx = -1, &
   frzdep_idx = -1

! pbuf indices for fields needed by heterogeneous freezing
integer :: &
   ast_idx = -1

type index_t
   integer :: bin_ndx
   integer :: spc_ndx
end type index_t

type(index_t),allocatable :: indices(:)
character(len=16),allocatable :: types(:)
character(len=fieldname_len),allocatable :: tot_dens_hnames(:)
character(len=fieldname_len),allocatable :: cld_dens_hnames(:)
character(len=fieldname_len),allocatable :: amb_dens_hnames(:)
character(len=fieldname_len),allocatable :: coated_dens_hnames(:)
character(len=fieldname_len),allocatable :: uncoated_dens_hnames(:)
character(len=fieldname_len),allocatable :: cldfn_dens_hnames(:)
character(len=fieldname_len),allocatable :: coated_frac_hnames(:)
character(len=fieldname_len),allocatable :: radius_hnames(:)
character(len=fieldname_len),allocatable :: wactfac_hnames(:)

integer :: tot_num_bins = 0

!===============================================================================
contains
!===============================================================================

subroutine hetfrz_classnuc_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical, mpi_real8, mpi_success

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'hetfrz_classnuc_cam_readnl'

  namelist /hetfrz_classnuc_nl/ hist_hetfrz_classnuc, hetfrz_bc_scalfac, hetfrz_dust_scalfac

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'hetfrz_classnuc_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, hetfrz_classnuc_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  ! Broadcast namelist variables
  call mpi_bcast(hist_hetfrz_classnuc, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: hist_hetfrz_classnuc")
  call mpi_bcast(hetfrz_bc_scalfac, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: hetfrz_bc_scalfac")
  call mpi_bcast(hetfrz_dust_scalfac, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: hetfrz_dust_scalfac")

  if (masterproc) then
     write(iulog,*) subname,': hist_hetfrz_classnuc = ',hist_hetfrz_classnuc
     write(iulog,*) subname,': hetfrz_bc_scalfac = ',hetfrz_bc_scalfac
     write(iulog,*) subname,': hetfrz_dust_scalfac = ',hetfrz_dust_scalfac
  end if

end subroutine hetfrz_classnuc_cam_readnl

!================================================================================================

subroutine hetfrz_classnuc_cam_register()

   if (.not. use_hetfrz_classnuc) return

   ! pbuf fields provided by hetfrz_classnuc
   call pbuf_add_field('FRZIMM', 'physpkg', dtype_r8, (/pcols,pver/), frzimm_idx)
   call pbuf_add_field('FRZCNT', 'physpkg', dtype_r8, (/pcols,pver/), frzcnt_idx)
   call pbuf_add_field('FRZDEP', 'physpkg', dtype_r8, (/pcols,pver/), frzdep_idx)

end subroutine hetfrz_classnuc_cam_register

!================================================================================================

subroutine hetfrz_classnuc_cam_init(mincld_in, aero_props)

   real(r8), intent(in) :: mincld_in
   class(aerosol_properties), intent(in) :: aero_props

   ! local variables
   integer :: istat, cnt, ibin, ispc
   character(len=42) :: tmpstr
   character(len=aero_name_len) :: species_type
   character(len=*), parameter :: routine = 'hetfrz_classnuc_cam_init'

   !--------------------------------------------------------------------------------------------

   if (.not. use_hetfrz_classnuc) return

   cnt = 0
   do ibin = 1, aero_props%nbins()
      do ispc = 1, aero_props%nspecies(ibin)
         if (aero_props%hetfrz_species(ibin,ispc)) then
            cnt = cnt+1
         end if
      end do
   end do

   tot_num_bins = cnt

   allocate(indices(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'indices', tot_num_bins)
   allocate(types(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'types', tot_num_bins)

   allocate(tot_dens_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'tot_dens_hnames', tot_num_bins)

   allocate(cld_dens_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'cld_dens_hnames', tot_num_bins)

   allocate(cldfn_dens_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'cldfn_dens_hnames', tot_num_bins)

   allocate(amb_dens_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'amb_dens_hnames', tot_num_bins)

   allocate(coated_dens_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'coated_dens_hnames', tot_num_bins)

   allocate(uncoated_dens_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'uncoated_dens_hnames', tot_num_bins)

   allocate(coated_frac_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'coated_frac_hnames', tot_num_bins)

   allocate(radius_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'radius_hnames', tot_num_bins)

   allocate(wactfac_hnames(tot_num_bins), stat=istat)
   call alloc_err(istat, routine, 'wactfac_hnames', tot_num_bins)

   cnt = 0
   do ibin = 1, aero_props%nbins()

      do ispc = 1, aero_props%nspecies(ibin)
         if (aero_props%hetfrz_species(ibin,ispc)) then
            call aero_props%species_type(ibin, ispc, species_type)
            cnt = cnt+1
            indices(cnt)%bin_ndx = ibin
            indices(cnt)%spc_ndx = ispc
            types(cnt) = trim(species_type)
            tmpstr = trim(species_type)//trim(int2str(ibin))

            cldfn_dens_hnames(cnt) = trim(tmpstr)//'_cld_fn'
            tot_dens_hnames(cnt) = trim(tmpstr)//'_tot_num'
            cld_dens_hnames(cnt) = trim(tmpstr)//'_cld_num'
            amb_dens_hnames(cnt) = trim(tmpstr)//'_amb_num'
            coated_dens_hnames(cnt) = trim(tmpstr)//'_coated'
            uncoated_dens_hnames(cnt) = trim(tmpstr)//'_uncoated'
            coated_frac_hnames(cnt) = trim(tmpstr)//'_coated_frac'
            radius_hnames(cnt) = trim(tmpstr)//'_radius'
            wactfac_hnames(cnt) = trim(tmpstr)//'_wactfac'

            call addfld(tot_dens_hnames(cnt),(/ 'lev' /), 'A', '#/cm3', &
                 'total '//trim(tmpstr)//' number density' )
            call addfld(cld_dens_hnames(cnt),(/ 'lev' /), 'A', '#/cm3', &
                 'cloud borne '//trim(tmpstr)//' number density' )
            call addfld(cldfn_dens_hnames(cnt),(/ 'lev' /), 'A', '#/cm3', &
                 'cloud borne '//trim(tmpstr)//' number density derived from fn' )
            call addfld(amb_dens_hnames(cnt),(/ 'lev' /), 'A', '#/cm3', &
                 'ambient '//trim(tmpstr)//' number density' )
            call addfld(coated_dens_hnames(cnt),(/ 'lev' /), 'A', '#/cm3', &
                 'coated '//trim(tmpstr)//' number density' )
            call addfld(uncoated_dens_hnames(cnt),(/ 'lev' /), 'A', '#/cm3', &
                 'uncoated '//trim(tmpstr)//' number density' )
            call addfld(coated_frac_hnames(cnt),(/ 'lev' /), 'A', '#/cm3', &
                 'coated '//trim(tmpstr)//' fraction' )
            call addfld(radius_hnames(cnt),(/ 'lev' /), 'A', 'm', &
                 'ambient '//trim(tmpstr)//' radius' )
            call addfld(wactfac_hnames(cnt),(/ 'lev' /), 'A', ' ', &
                 trim(tmpstr)//' water activity mass factor' )

         end if
      end do

   end do

   mincld = mincld_in

   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   ! pbuf fields used by hetfrz_classnuc
   ast_idx      = pbuf_get_index('AST')

   call addfld('FRZIMM', (/ 'lev' /), 'A', ' ', 'immersion  freezing')
   call addfld('FRZCNT', (/ 'lev' /), 'A', ' ', 'contact    freezing')
   call addfld('FRZDEP', (/ 'lev' /), 'A', ' ', 'deposition freezing')
   call addfld('FREQIMM', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of immersion  freezing')
   call addfld('FREQCNT', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of contact    freezing')
   call addfld('FREQDEP', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of deposition freezing')
   call addfld('FREQMIX', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of mixed-phase clouds' )

   call addfld('DSTFREZIMM', (/ 'lev' /), 'A', 'm-3s-1', 'dust immersion  freezing rate')
   call addfld('DSTFREZCNT', (/ 'lev' /), 'A', 'm-3s-1', 'dust contact    freezing rate')
   call addfld('DSTFREZDEP', (/ 'lev' /), 'A', 'm-3s-1', 'dust deposition freezing rate')

   call addfld('BCFREZIMM', (/ 'lev' /), 'A', 'm-3s-1', 'bc immersion  freezing rate')
   call addfld('BCFREZCNT', (/ 'lev' /), 'A', 'm-3s-1', 'bc contact    freezing rate')
   call addfld('BCFREZDEP', (/ 'lev' /), 'A', 'm-3s-1', 'bc deposition freezing rate')

   call addfld('NIMIX_IMM', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to het immersion freezing in Mixed Clouds')
   call addfld('NIMIX_CNT', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to het contact freezing in Mixed Clouds')
   call addfld('NIMIX_DEP', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to het deposition freezing in Mixed Clouds')

   call addfld('DSTNIDEP', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to dst dep freezing in Mixed Clouds')
   call addfld('DSTNICNT', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to dst cnt freezing in Mixed Clouds')
   call addfld('DSTNIIMM', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to dst imm freezing in Mixed Clouds')

   call addfld('BCNIDEP', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to bc dep freezing in Mixed Clouds')
   call addfld('BCNICNT', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to bc cnt freezing in Mixed Clouds')
   call addfld('BCNIIMM', (/ 'lev' /), 'A', '#/m3', &
               'Activated Ice Number Concentration due to bc imm freezing in Mixed Clouds')

   call addfld('NUMICE10s', (/ 'lev' /), 'A', '#/m3', &
               'Ice Number Concentration due to het freezing in Mixed Clouds during 10-s period')
   call addfld('NUMIMM10sDST', (/ 'lev' /), 'A', '#/m3', &
               'Ice Number Concentration due to imm freezing by dst in Mixed Clouds during 10-s period')
   call addfld('NUMIMM10sBC', (/ 'lev' /), 'A', '#/m3', &
               'Ice Number Concentration due to imm freezing by bc in Mixed Clouds during 10-s period')

   call addfld('hetfrz_supersatice',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_fn1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_fn2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_fn3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_fn4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_r3lx',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_icnlx',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_radius1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_radius2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_radius3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_radius4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awcam1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awcam2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awcam3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awcam4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awfacm1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awfacm2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awfacm3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
!!$   call addfld('hetfrz_awfacm4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_dstcoat1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_dstcoat2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_dstcoat3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_dstcoat4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totaernum1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totaernum2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totaernum3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totaernum4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_uncoatednum1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_uncoatednum2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_uncoatednum3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_uncoatednum4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totinum1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totinum2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totinum3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totinum4',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totcnum1',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totcnum2',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totcnum3',(/ 'lev' /), 'A', ' ','hetfrz debug field')
   call addfld('hetfrz_totcnum4',(/ 'lev' /), 'A', ' ','hetfrz debug field')


   call add_default('hetfrz_supersatice', 3, ' ')
   call add_default('hetfrz_fn1', 3, ' ')
   call add_default('hetfrz_fn2', 3, ' ')
   call add_default('hetfrz_fn3', 3, ' ')
   call add_default('hetfrz_fn4', 3, ' ')
   call add_default('hetfrz_r3lx', 3, ' ')
   call add_default('hetfrz_icnlx', 3, ' ')
   call add_default('hetfrz_radius1', 3, ' ')
   call add_default('hetfrz_radius2', 3, ' ')
   call add_default('hetfrz_radius3', 3, ' ')
   call add_default('hetfrz_radius4', 3, ' ')
!!$   call add_default('hetfrz_awcam1', 3, ' ')
!!$   call add_default('hetfrz_awcam2', 3, ' ')
!!$   call add_default('hetfrz_awcam3', 3, ' ')
!!$   call add_default('hetfrz_awcam4', 3, ' ')
!!$   call add_default('hetfrz_awfacm1', 3, ' ')
!!$   call add_default('hetfrz_awfacm2', 3, ' ')
!!$   call add_default('hetfrz_awfacm3', 3, ' ')
!!$   call add_default('hetfrz_awfacm4', 3, ' ')
   call add_default('hetfrz_dstcoat1', 3, ' ')
   call add_default('hetfrz_dstcoat2', 3, ' ')
   call add_default('hetfrz_dstcoat3', 3, ' ')
   call add_default('hetfrz_dstcoat4', 3, ' ')
   call add_default('hetfrz_totaernum1', 3, ' ')
   call add_default('hetfrz_totaernum2', 3, ' ')
   call add_default('hetfrz_totaernum3', 3, ' ')
   call add_default('hetfrz_totaernum4', 3, ' ')
   call add_default('hetfrz_uncoatednum1', 3, ' ')
   call add_default('hetfrz_uncoatednum2', 3, ' ')
   call add_default('hetfrz_uncoatednum3', 3, ' ')
   call add_default('hetfrz_uncoatednum4', 3, ' ')
   call add_default('hetfrz_totinum1', 3, ' ')
   call add_default('hetfrz_totinum2', 3, ' ')
   call add_default('hetfrz_totinum3', 3, ' ')
   call add_default('hetfrz_totinum4', 3, ' ')
   call add_default('hetfrz_totcnum1', 3, ' ')
   call add_default('hetfrz_totcnum2', 3, ' ')
   call add_default('hetfrz_totcnum3', 3, ' ')
   call add_default('hetfrz_totcnum4', 3, ' ')

   if (hist_hetfrz_classnuc) then

      call add_default('FREQIMM', 1, ' ')
      call add_default('FREQCNT', 1, ' ')
      call add_default('FREQDEP', 1, ' ')
      call add_default('FREQMIX', 1, ' ')

      call add_default('DSTFREZIMM', 1, ' ')
      call add_default('DSTFREZCNT', 1, ' ')
      call add_default('DSTFREZDEP', 1, ' ')

      call add_default('BCFREZIMM', 1, ' ')
      call add_default('BCFREZCNT', 1, ' ')
      call add_default('BCFREZDEP', 1, ' ')

      call add_default('NIMIX_IMM', 1, ' ')
      call add_default('NIMIX_CNT', 1, ' ')
      call add_default('NIMIX_DEP', 1, ' ')

      call add_default('DSTNIDEP', 1, ' ')
      call add_default('DSTNICNT', 1, ' ')
      call add_default('DSTNIIMM', 1, ' ')

      call add_default('BCNIDEP', 1, ' ')
      call add_default('BCNICNT', 1, ' ')
      call add_default('BCNIIMM', 1, ' ')

      call add_default('NUMICE10s', 1, ' ')
      call add_default('NUMIMM10sDST', 1, ' ')
      call add_default('NUMIMM10sBC', 1, ' ')

   end if

   call hetfrz_classnuc_init(rair, cpair, rh2o, rhoh2o, mwh2o, tmelt, iulog, &
        hetfrz_bc_scalfac, hetfrz_dust_scalfac )

end subroutine hetfrz_classnuc_cam_init

!================================================================================================

subroutine hetfrz_classnuc_cam_calc( aero_props, aero_state, &
   state, deltatin, factnum, pbuf)

   ! arguments
   class(aerosol_properties),optional, intent(in) :: aero_props
   class(aerosol_state),optional, intent(in) :: aero_state
   type(physics_state), target, intent(in)    :: state
   real(r8),                    intent(in)    :: deltatin       ! time step (s)
   real(r8),                    intent(in)    :: factnum(:,:,:) ! activation fraction for aerosol number
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! local workspace

   ! outputs shared with the microphysics via the pbuf
   real(r8), pointer :: frzimm(:,:)
   real(r8), pointer :: frzcnt(:,:)
   real(r8), pointer :: frzdep(:,:)

   integer :: itim_old
   integer :: i, k

   real(r8) :: rho(pcols,pver)          ! air density (kg m-3)

   real(r8), pointer :: ast(:,:)

   real(r8) :: lcldm(pcols,pver)

   real(r8) :: esi(pcols), esl(pcols)
   real(r8) :: con1, r3lx, supersatice

   real(r8) :: qcic
   real(r8) :: ncic

   real(r8) :: frzbcimm(pcols,pver), frzduimm(pcols,pver)
   real(r8) :: frzbccnt(pcols,pver), frzducnt(pcols,pver)
   real(r8) :: frzbcdep(pcols,pver), frzdudep(pcols,pver)

   real(r8) :: freqimm(pcols,pver), freqcnt(pcols,pver), freqdep(pcols,pver), freqmix(pcols,pver)

   real(r8) :: nnuccc_bc(pcols,pver), nnucct_bc(pcols,pver), nnudep_bc(pcols,pver)
   real(r8) :: nnuccc_dst(pcols,pver), nnucct_dst(pcols,pver), nnudep_dst(pcols,pver)
   real(r8) :: niimm_bc(pcols,pver), nicnt_bc(pcols,pver), nidep_bc(pcols,pver)
   real(r8) :: niimm_dst(pcols,pver), nicnt_dst(pcols,pver), nidep_dst(pcols,pver)
   real(r8) :: numice10s(pcols,pver)
   real(r8) :: numice10s_imm_dst(pcols,pver)
   real(r8) :: numice10s_imm_bc(pcols,pver)

   real(r8) :: coated(pcols,pver,tot_num_bins)
   real(r8) :: aer_radius(pcols,pver,tot_num_bins)
   real(r8) :: aer_wactfac(pcols,pver,tot_num_bins)

   real(r8) :: coated_amb_aer_num(pcols,pver,tot_num_bins)
   real(r8) :: uncoated_amb_aer_num(pcols,pver,tot_num_bins)
   real(r8) :: amb_aer_num(pcols,pver,tot_num_bins)
   real(r8) :: cld_aer_num(pcols,pver,tot_num_bins)
   real(r8) :: tot_aer_num(pcols,pver,tot_num_bins)
   real(r8) :: fn_cld_aer_num(pcols,pver)
   real(r8) :: fraction_activated(pcols,pver,tot_num_bins)

   character(128) :: errstring   ! Error status
   !-------------------------------------------------------------------------------

   real(r8) :: supersatice_tmp(pcols,pver)
   real(r8) :: r3lx_tmp(pcols,pver)
   real(r8) :: icnlx_tmp(pcols,pver)

   supersatice_tmp = 0.0_r8
   r3lx_tmp = 0.0_r8
   icnlx_tmp = 0.0_r8

   associate( &
      lchnk => state%lchnk,             &
      ncol  => state%ncol,              &
      t     => state%t,                 &
      qc    => state%q(:pcols,:pver,cldliq_idx), &
      nc    => state%q(:pcols,:pver,numliq_idx), &
      pmid  => state%pmid               )

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   rho(:,:) = 0._r8

   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

   do k = top_lev, pver
      do i = 1, ncol
         lcldm(i,k) = max(ast(i,k), mincld)
      end do
   end do

   do i = 1,tot_num_bins

      call aero_state%get_amb_species_numdens( indices(i)%bin_ndx, ncol, pver, types(i), aero_props, rho, amb_aer_num(:,:,i))
      call aero_state%get_cld_species_numdens( indices(i)%bin_ndx, ncol, pver, types(i), aero_props, rho, cld_aer_num(:,:,i))

      tot_aer_num(:ncol,:,i) = cld_aer_num(:ncol,:,i) + amb_aer_num(:ncol,:,i)

      call outfld(tot_dens_hnames(i), tot_aer_num(:,:,i), pcols, lchnk)
      call outfld(amb_dens_hnames(i), amb_aer_num(:,:,i), pcols, lchnk)
      call outfld(cld_dens_hnames(i), cld_aer_num(:,:,i), pcols, lchnk)

      aer_radius(:ncol,:,i) = aero_state%mass_mean_radius( indices(i)%bin_ndx, indices(i)%spc_ndx,  ncol, pver, aero_props, rho )

      coated(:ncol,:,i) = aero_state%coated_frac( indices(i)%bin_ndx, types(i), ncol, pver, aero_props, aer_radius(:,:,i) )

      call outfld(coated_frac_hnames(i), coated(:,:,i), pcols, lchnk)

      coated_amb_aer_num(:ncol,:,i) = amb_aer_num(:ncol,:,i)*coated(:ncol,:,i)
      uncoated_amb_aer_num(:ncol,:,i) = amb_aer_num(:ncol,:,i)*(1._r8-coated(:ncol,:,i))

      call outfld(coated_dens_hnames(i), coated_amb_aer_num(:,:,i), pcols, lchnk)
      call outfld(uncoated_dens_hnames(i), uncoated_amb_aer_num(:,:,i), pcols, lchnk)
      call outfld(radius_hnames(i), aer_radius(:ncol,:,i), ncol, lchnk)

      call aero_state%watact_mfactor(indices(i)%bin_ndx, types(i), ncol, pver, aero_props, rho, aer_wactfac(:ncol,:,i))
      call outfld(wactfac_hnames(i), aer_wactfac(:,:,i), pcols, lchnk)

      fn_cld_aer_num(:ncol,:) = tot_aer_num(:ncol,:,i)*factnum(:ncol,:,indices(i)%bin_ndx)
      call outfld(cldfn_dens_hnames(i), fn_cld_aer_num, pcols, lchnk)

      fraction_activated(:ncol,:,i) = factnum(:ncol,:,indices(i)%bin_ndx)

   end do

   ! frzimm, frzcnt, frzdep are the outputs of this parameterization used by the microphysics
   call pbuf_get_field(pbuf, frzimm_idx, frzimm)
   call pbuf_get_field(pbuf, frzcnt_idx, frzcnt)
   call pbuf_get_field(pbuf, frzdep_idx, frzdep)

   frzimm(:ncol,:) = 0._r8
   frzcnt(:ncol,:) = 0._r8
   frzdep(:ncol,:) = 0._r8

   frzbcimm(:ncol,:) = 0._r8
   frzduimm(:ncol,:) = 0._r8
   frzbccnt(:ncol,:) = 0._r8
   frzducnt(:ncol,:) = 0._r8
   frzbcdep(:ncol,:) = 0._r8
   frzdudep(:ncol,:) = 0._r8

   freqimm(:ncol,:) = 0._r8
   freqcnt(:ncol,:) = 0._r8
   freqdep(:ncol,:) = 0._r8
   freqmix(:ncol,:) = 0._r8

   numice10s(:ncol,:)         = 0._r8
   numice10s_imm_dst(:ncol,:) = 0._r8
   numice10s_imm_bc(:ncol,:)  = 0._r8

   nnuccc_bc(:,:) = 0._r8
   nnucct_bc(:,:) = 0._r8
   nnudep_bc(:,:) = 0._r8

   nnuccc_dst(:,:) = 0._r8
   nnucct_dst(:,:) = 0._r8
   nnudep_dst(:,:) = 0._r8

   niimm_bc(:,:) = 0._r8
   nicnt_bc(:,:) = 0._r8
   nidep_bc(:,:) = 0._r8

   niimm_dst(:,:) = 0._r8
   nicnt_dst(:,:) = 0._r8
   nidep_dst(:,:) = 0._r8

   do k = top_lev, pver
      call svp_water_vect(t(1:ncol,k), esl(1:ncol), ncol)
      call svp_ice_vect(t(1:ncol,k), esi(1:ncol), ncol)
      do i = 1, ncol

         if (t(i,k) > 235.15_r8 .and. t(i,k) < 269.15_r8) then
            qcic = min(qc(i,k)/lcldm(i,k), 5.e-3_r8)
            ncic = max(nc(i,k)/lcldm(i,k), 0._r8)

            con1 = 1._r8/(1.333_r8*pi)**0.333_r8
            r3lx = con1*(rho(i,k)*qcic/(rhoh2o*max(ncic*rho(i,k), 1.0e6_r8)))**0.333_r8 ! in m
            r3lx = max(4.e-6_r8, r3lx)
            supersatice = esl(i)/esi(i)

            supersatice_tmp(i,k) = supersatice

            r3lx_tmp(i,k) = r3lx
            icnlx_tmp(i,k) = ncic*rho(i,k)*1.0e-6_r8

            call hetfrz_classnuc_calc( tot_num_bins, types, &
               deltatin,  t(i,k),  pmid(i,k),  supersatice,   &
               fraction_activated(i,k,:),  r3lx,  ncic*rho(i,k)*1.0e-6_r8,  frzbcimm(i,k),  frzduimm(i,k),   &
               frzbccnt(i,k),  frzducnt(i,k),  frzbcdep(i,k),  frzdudep(i,k),  aer_radius(i,k,:), &
               aer_wactfac(i,k,:), coated(i,k,:), tot_aer_num(i,k,:),  &
               uncoated_amb_aer_num(i,k,:), amb_aer_num(i,k,:), &
               cld_aer_num(i,k,:), errstring)

            call handle_errmsg(errstring, subname="hetfrz_classnuc_calc")

            frzimm(i,k) = frzbcimm(i,k) + frzduimm(i,k)
            frzcnt(i,k) = frzbccnt(i,k) + frzducnt(i,k)
            frzdep(i,k) = frzbcdep(i,k) + frzdudep(i,k)

            if (frzimm(i,k) > 0._r8) freqimm(i,k) = 1._r8
            if (frzcnt(i,k) > 0._r8) freqcnt(i,k) = 1._r8
            if (frzdep(i,k) > 0._r8) freqdep(i,k) = 1._r8
            if ((frzimm(i,k) + frzcnt(i,k) + frzdep(i,k)) > 0._r8) freqmix(i,k) = 1._r8

         else
            frzimm(i,k) = 0._r8
            frzcnt(i,k) = 0._r8
            frzdep(i,k) = 0._r8
         end if

         nnuccc_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*ast(i,k)
         nnucct_bc(i,k) = frzbccnt(i,k)*1.0e6_r8*ast(i,k)
         nnudep_bc(i,k) = frzbcdep(i,k)*1.0e6_r8*ast(i,k)

         nnuccc_dst(i,k) = frzduimm(i,k)*1.0e6_r8*ast(i,k)
         nnucct_dst(i,k) = frzducnt(i,k)*1.0e6_r8*ast(i,k)
         nnudep_dst(i,k) = frzdudep(i,k)*1.0e6_r8*ast(i,k)

         niimm_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*deltatin
         nicnt_bc(i,k) = frzbccnt(i,k)*1.0e6_r8*deltatin
         nidep_bc(i,k) = frzbcdep(i,k)*1.0e6_r8*deltatin

         niimm_dst(i,k) = frzduimm(i,k)*1.0e6_r8*deltatin
         nicnt_dst(i,k) = frzducnt(i,k)*1.0e6_r8*deltatin
         nidep_dst(i,k) = frzdudep(i,k)*1.0e6_r8*deltatin

         numice10s(i,k) = (frzimm(i,k)+frzcnt(i,k)+frzdep(i,k))*1.0e6_r8*deltatin*(10._r8/deltatin)
         numice10s_imm_dst(i,k) = frzduimm(i,k)*1.0e6_r8*deltatin*(10._r8/deltatin)
         numice10s_imm_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*deltatin*(10._r8/deltatin)
      end do
   end do

   call outfld('hetfrz_supersatice', supersatice_tmp(:ncol,:), ncol, lchnk)
   call outfld('hetfrz_fn1', fraction_activated(:ncol,:,1), ncol, lchnk)
   call outfld('hetfrz_fn2', fraction_activated(:ncol,:,2), ncol, lchnk)
   call outfld('hetfrz_fn3', fraction_activated(:ncol,:,3), ncol, lchnk)
   call outfld('hetfrz_fn4', fraction_activated(:ncol,:,4), ncol, lchnk)
   call outfld('hetfrz_r3lx', r3lx_tmp(:ncol,:), ncol, lchnk)
   call outfld('hetfrz_icnlx', icnlx_tmp(:ncol,:), ncol, lchnk)

   call outfld('hetfrz_radius1', aer_radius(:ncol,:,1), ncol, lchnk)
   call outfld('hetfrz_radius2', aer_radius(:ncol,:,2), ncol, lchnk)
   call outfld('hetfrz_radius3', aer_radius(:ncol,:,3), ncol, lchnk)
   call outfld('hetfrz_radius4', aer_radius(:ncol,:,4), ncol, lchnk)

!!$
!!$   call outfld('hetfrz_awcam1', awcam(:ncol,:,1), ncol, lchnk)
!!$   call outfld('hetfrz_awcam2', awcam(:ncol,:,2), ncol, lchnk)
!!$   call outfld('hetfrz_awcam3', awcam(:ncol,:,3), ncol, lchnk)
!!$   call outfld('hetfrz_awcam4', awcam(:ncol,:,4), ncol, lchnk)
!!$
!!$   call outfld('hetfrz_awfacm1', awfacm(:ncol,:,1), ncol, lchnk)
!!$   call outfld('hetfrz_awfacm2', awfacm(:ncol,:,2), ncol, lchnk)
!!$   call outfld('hetfrz_awfacm3', awfacm(:ncol,:,3), ncol, lchnk)
!!$   call outfld('hetfrz_awfacm4', awfacm(:ncol,:,4), ncol, lchnk)

   call outfld('hetfrz_dstcoat1', coated(:ncol,:,1), ncol, lchnk)
   call outfld('hetfrz_dstcoat2', coated(:ncol,:,2), ncol, lchnk)
   call outfld('hetfrz_dstcoat3', coated(:ncol,:,3), ncol, lchnk)
   call outfld('hetfrz_dstcoat4', coated(:ncol,:,4), ncol, lchnk)

   call outfld('hetfrz_totaernum1', tot_aer_num(:ncol,:,1), ncol, lchnk)
   call outfld('hetfrz_totaernum2', tot_aer_num(:ncol,:,2), ncol, lchnk)
   call outfld('hetfrz_totaernum3', tot_aer_num(:ncol,:,3), ncol, lchnk)
   call outfld('hetfrz_totaernum4', tot_aer_num(:ncol,:,4), ncol, lchnk)

   call outfld('hetfrz_uncoatednum1', uncoated_amb_aer_num(:ncol,:,1), ncol, lchnk)
   call outfld('hetfrz_uncoatednum2', uncoated_amb_aer_num(:ncol,:,2), ncol, lchnk)
   call outfld('hetfrz_uncoatednum3', uncoated_amb_aer_num(:ncol,:,3), ncol, lchnk)
   call outfld('hetfrz_uncoatednum4', uncoated_amb_aer_num(:ncol,:,4), ncol, lchnk)

   call outfld('hetfrz_totinum1', amb_aer_num(:ncol,:,1), ncol, lchnk)
   call outfld('hetfrz_totinum2', amb_aer_num(:ncol,:,2), ncol, lchnk)
   call outfld('hetfrz_totinum3', amb_aer_num(:ncol,:,3), ncol, lchnk)
   call outfld('hetfrz_totinum4', amb_aer_num(:ncol,:,4), ncol, lchnk)

   call outfld('hetfrz_totcnum1', cld_aer_num(:ncol,:,1), ncol, lchnk)
   call outfld('hetfrz_totcnum2', cld_aer_num(:ncol,:,2), ncol, lchnk)
   call outfld('hetfrz_totcnum3', cld_aer_num(:ncol,:,3), ncol, lchnk)
   call outfld('hetfrz_totcnum4', cld_aer_num(:ncol,:,4), ncol, lchnk)


   call outfld('FRZIMM', frzimm, pcols, lchnk)
   call outfld('FRZCNT', frzcnt, pcols, lchnk)
   call outfld('FRZDEP', frzdep, pcols, lchnk)

   call outfld('FREQIMM', freqimm, pcols, lchnk)
   call outfld('FREQCNT', freqcnt, pcols, lchnk)
   call outfld('FREQDEP', freqdep, pcols, lchnk)
   call outfld('FREQMIX', freqmix, pcols, lchnk)

   call outfld('DSTFREZIMM', nnuccc_dst, pcols, lchnk)
   call outfld('DSTFREZCNT', nnucct_dst, pcols, lchnk)
   call outfld('DSTFREZDEP', nnudep_dst, pcols, lchnk)

   call outfld('BCFREZIMM', nnuccc_bc, pcols, lchnk)
   call outfld('BCFREZCNT', nnucct_bc, pcols, lchnk)
   call outfld('BCFREZDEP', nnudep_bc, pcols, lchnk)

   call outfld('NIMIX_IMM', niimm_bc+niimm_dst, pcols, lchnk)
   call outfld('NIMIX_CNT', nicnt_bc+nicnt_dst, pcols, lchnk)
   call outfld('NIMIX_DEP', nidep_bc+nidep_dst, pcols, lchnk)

   call outfld('DSTNICNT', nicnt_dst, pcols, lchnk)
   call outfld('DSTNIDEP', nidep_dst, pcols, lchnk)
   call outfld('DSTNIIMM', niimm_dst, pcols, lchnk)

   call outfld('BCNICNT', nicnt_bc, pcols, lchnk)
   call outfld('BCNIDEP', nidep_bc, pcols, lchnk)
   call outfld('BCNIIMM', niimm_bc, pcols, lchnk)

   call outfld('NUMICE10s', numice10s, pcols, lchnk)
   call outfld('NUMIMM10sDST', numice10s_imm_dst, pcols, lchnk)
   call outfld('NUMIMM10sBC', numice10s_imm_bc, pcols, lchnk)

   end associate

end subroutine hetfrz_classnuc_cam_calc

!====================================================================================================

end module hetfrz_classnuc_cam
