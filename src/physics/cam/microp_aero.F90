module microp_aero

!---------------------------------------------------------------------------------
! Purpose:
!   CAM driver layer for aerosol activation processes.
!
! ***N.B.*** This module is currently hardcoded to recognize only the aerosols/modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!
! Author: Andrew Gettelman
! Based on code from: Hugh Morrison, Xiaohong Liu and Steve Ghan
! May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Andrew Gettelman  (andrew@ucar.edu)
! Modifications: A. Gettelman Nov 2010  - changed to support separation of
!                  microphysics and macrophysics and concentrate aerosol information here
!                B. Eaton, Sep 2014 - Refactored to move CAM interface code into the CAM
!                  interface modules and preserve just the driver layer functionality here.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use ref_pres,         only: top_lev => trop_cloud_top_lev
use physconst,        only: rair, gravit, tmelt, cpair, rh2o, rhoh2o, latvap, &
                            r_universal, mwh2o
use constituents,     only: cnst_get_ind
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init, physics_ptend_sum, &
                            physics_state_copy, physics_update
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,     only: phys_getopts, use_hetfrz_classnuc
use aerosol_instances_mod, only: aerosol_instances_get_num_models, &
                                 aerosol_instances_is_active, &
                                 aerosol_instances_get_props, &
                                 aerosol_instances_create_states, &
                                 aerosol_instances_destroy_states, &
                                 aero_state_entry_t

use nucleate_ice_cam, only: use_preexisting_ice, nucleate_ice_cam_readnl, nucleate_ice_cam_register, &
                            nucleate_ice_cam_init, nucleate_ice_cam_calc

use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_calc, &
                            aername, psat, ccn_name

use hetfrz_classnuc_cam, only: hetfrz_classnuc_cam_readnl, hetfrz_classnuc_cam_register, hetfrz_classnuc_cam_init, &
                               hetfrz_classnuc_cam_calc

use compute_subgrid_vertical_velocity, only: compute_subgrid_vertical_velocity_tke_run, &
                                             compute_subgrid_vertical_velocity_kvh_run, &
                                             compute_subgrid_vertical_velocity_clubb_run
use scale_subgrid_vertical_velocity,   only: scale_subgrid_vertical_velocity_run

use cam_history,      only: addfld, add_default, outfld
use cam_logfile,      only: iulog
use cam_abortutils,       only: endrun

use aerosol_properties_mod, only: aerosol_properties

use aerosol_state_mod, only: aerosol_state

implicit none
private
save

public :: microp_aero_init, microp_aero_run, microp_aero_readnl, microp_aero_register
public :: microp_aero_final
public :: aerosol_state_object
public :: aerosol_properties_object

! Private module data
character(len=16)   :: eddy_scheme
real(r8), parameter :: unset_r8 = huge(1.0_r8)

! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

! Namelist parameters
real(r8) :: npccn_scale   ! scaling for activated number
real(r8) :: wsub_scale    ! scaling for sub-grid vertical velocity (liquid)
real(r8) :: wsubi_scale   ! scaling for sub-grid vertical velocity (ice)
real(r8) :: wsub_min      ! minimum sub-grid vertical velocity (liquid) before scale factor
real(r8) :: wsub_min_asf  ! minimum sub-grid vertical velocity (liquid) after scale factor
real(r8) :: wsubi_min     ! minimum sub-grid vertical velocity (ice)

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter :: mincld = 0.0001_r8

! indices in state%q and pbuf structures
integer :: cldliq_idx = -1
integer :: cldice_idx = -1
integer :: numliq_idx = -1
integer :: numice_idx = -1
integer :: kvh_idx = -1
integer :: tke_idx = -1
integer :: wp2_idx = -1
integer :: ast_idx = -1
integer :: cldo_idx = -1
integer :: dgnumwet_idx = -1

! carma aerosols
logical :: clim_carma_aero

! modal aerosols
logical :: clim_modal_aero

integer :: iaermod_clim_modal_carma

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode
integer :: coarse_so4_idx = -1  ! index of sulfate in coarse mode

integer :: naer_all = 0
integer :: npccn_idx, rndst_idx, nacon_idx

logical  :: separate_dust = .false.

class(aerosol_properties), pointer :: aero_props_obj=>null()

!=========================================================================================
contains
!=========================================================================================

subroutine microp_aero_register
   !-----------------------------------------------------------------------
   !
   ! Purpose:
   ! Register pbuf fields for aerosols needed by microphysics
   !
   ! Author: Cheryl Craig October 2012
   !
   !-----------------------------------------------------------------------
   use ppgrid,         only: pcols
   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)

   call pbuf_add_field('RNDST',      'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
   call pbuf_add_field('NACON',      'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)

   call nucleate_ice_cam_register()
   call hetfrz_classnuc_cam_register()

end subroutine microp_aero_register

!=========================================================================================

subroutine microp_aero_init(phys_state,pbuf2d)

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   ! Initialize constants for aerosols needed by microphysics
   !
   ! Author: Andrew Gettelman May 2010
   !
   !-----------------------------------------------------------------------

   type(physics_state), pointer       :: phys_state(:)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   integer  :: iaer
   integer  :: m, n, nspec
   integer  :: iaermod

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg

   character(len=512) :: errmsg
   integer            :: errflg

   class(aerosol_properties), pointer :: aero_props_bulk => null()

   !-----------------------------------------------------------------------

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out          = eddy_scheme,  &
                     history_amwg_out         = history_amwg )

   ! Access the physical properties of the aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   ! get indices into state and pbuf structures
   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   select case(trim(eddy_scheme))
   case ('diag_TKE')
      tke_idx      = pbuf_get_index('tke')
   case ('CLUBB_SGS')
      wp2_idx = pbuf_get_index('WP2_nadv')
   case default
      kvh_idx      = pbuf_get_index('kvh')
   end select

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   clim_modal_aero = aerosol_instances_is_active('modal')
   clim_carma_aero = aerosol_instances_is_active('carma')
   iaermod_clim_modal_carma = -1

   ast_idx = pbuf_get_index('AST')

   if (clim_modal_aero .or. clim_carma_aero) then
      cldo_idx = pbuf_get_index('CLDO')
      ! Get modal/CARMA properties object from factory (factory owns the object)
      do iaermod = 1, aerosol_instances_get_num_models()
         aero_props_obj => aerosol_instances_get_props(iaermod, 0)
         if (associated(aero_props_obj)) then
            if (aero_props_obj%model_is('modal') .or. aero_props_obj%model_is('CARMA')) then
               ! store idx for providing to dycore via aerosol_state_object...
               iaermod_clim_modal_carma = iaermod
               exit
            end if
         end if
      end do
      call ndrop_init(aero_props_obj)
   end if

   if (clim_modal_aero) then

      dgnumwet_idx = pbuf_get_index('DGNUMWET')

      ! Init indices for specific modes/species

      ! mode index for specified mode types
      do m = 1, aero_props_obj%nbins()
         str32 = aero_props_obj%bin_name(m)
         select case (trim(str32))
         case ('accum')
            mode_accum_idx = m
         case ('aitken')
            mode_aitken_idx = m
         case ('coarse')
            mode_coarse_idx = m
         case ('coarse_dust')
            mode_coarse_dst_idx = m
         case ('coarse_seasalt')
            mode_coarse_slt_idx = m
         end select
      end do

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode
      if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
      if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! species indices for specified types
      ! find indices for the dust and seasalt species in the coarse mode
      nspec = aero_props_obj%nspecies(mode_coarse_dst_idx)
      do n = 1, nspec
         call aero_props_obj%species_type(mode_coarse_dst_idx, n, str32)
         select case (trim(str32))
         case ('dust')
            coarse_dust_idx = n
         end select
      end do
      nspec = aero_props_obj%nspecies(mode_coarse_slt_idx)
      do n = 1, nspec
         call aero_props_obj%species_type(mode_coarse_slt_idx, n, str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do
      if (mode_coarse_idx>0) then
         nspec = aero_props_obj%nspecies(mode_coarse_idx)
         do n = 1, nspec
            call aero_props_obj%species_type(mode_coarse_idx, n, str32)
            select case (trim(str32))
            case ('sulfate')
               coarse_so4_idx = n
            end select
         end do
      endif

      ! Check that required mode specie types were found
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

   else if (.not.clim_carma_aero) then

      ! Find bulk properties object from factory
      aero_props_bulk => null()
      do iaermod = 1, aerosol_instances_get_num_models()
         aero_props_bulk => aerosol_instances_get_props(iaermod, 0)
         if (associated(aero_props_bulk)) then
            if (aero_props_bulk%model_is('BAM')) exit
         end if
         aero_props_bulk => null()
      end do

      call ndrop_bam_init(masterproc, iulog, &
           mwh2o=mwh2o, r_universal=r_universal, tmelt=tmelt, rhoh2o=rhoh2o, &
           naer_all_out=naer_all, errmsg=errmsg, errflg=errflg)
      if (errflg /= 0) then
         call endrun(routine//': ndrop_bam_init: '//trim(errmsg))
      end if

      ! Set module-level props object for BAM (used by nucleate_ice_cam)
      aero_props_obj => aero_props_bulk

   end if

   call addfld('LCLOUD', (/ 'lev' /), 'A', ' ',   'Liquid cloud fraction used in stratus activation', sampled_on_subcycle=.true.)

   call addfld('WSUB',   (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity',            sampled_on_subcycle=.true.)
   call addfld('WSUBI',  (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity for ice',    sampled_on_subcycle=.true.)

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
   end if

   ! BAM-specific history fields (moved from ndrop_bam for CCPP portability)
   if (.not.clim_modal_aero .and. .not.clim_carma_aero) then
      do iaer = 1, naer_all
         call addfld(trim(aername(iaer))//'_m3', (/ 'lev' /), 'A', 'm-3', &
              'aerosol number concentration', sampled_on_subcycle=.true.)
      end do
      call addfld ('CCN1', (/ 'lev' /), 'A', '#/cm3', 'CCN concentration at S=0.02%', sampled_on_subcycle=.true.)
      call addfld ('CCN2', (/ 'lev' /), 'A', '#/cm3', 'CCN concentration at S=0.05%', sampled_on_subcycle=.true.)
      call addfld ('CCN3', (/ 'lev' /), 'A', '#/cm3', 'CCN concentration at S=0.1%',  sampled_on_subcycle=.true.)
      call addfld ('CCN4', (/ 'lev' /), 'A', '#/cm3', 'CCN concentration at S=0.2%',  sampled_on_subcycle=.true.)
      call addfld ('CCN5', (/ 'lev' /), 'A', '#/cm3', 'CCN concentration at S=0.5%',  sampled_on_subcycle=.true.)
      call addfld ('CCN6', (/ 'lev' /), 'A', '#/cm3', 'CCN concentration at S=1.0%',  sampled_on_subcycle=.true.)
      if (history_amwg) then
         call add_default('CCN3', 1, ' ')
      end if
   end if

   if (associated(aero_props_obj)) then
      call nucleate_ice_cam_init(mincld, pbuf2d, aero_props=aero_props_obj)
   else
      ! this path is used for aquaplanet compsets with two-moment microphysics
      ! where nucleatei still runs Meyers nucleation deposition even without aerosol.
      call nucleate_ice_cam_init(mincld, pbuf2d)
   end if
   if (use_hetfrz_classnuc) then
      if (associated(aero_props_obj)) then
         call hetfrz_classnuc_cam_init(mincld, aero_props_obj)
      else
         call endrun(routine//': cannot use hetfrz_classnuc without prognostic aerosols')
      endif
   endif

end subroutine microp_aero_init

!=========================================================================================
! returns a pointer to an aerosol state object for a given chunk index
! compatibility: for use by the dycore
function aerosol_state_object(lchnk) result(obj)
   use aerosol_instances_mod, only: aerosol_instances_get_state

  integer,intent(in) :: lchnk ! local chunk index
  class(aerosol_state), pointer :: obj ! aerosol state object pointer for local chunk

  if (iaermod_clim_modal_carma > 0) then
    obj => aerosol_instances_get_state(iaermod_clim_modal_carma, list_idx=0, lchnk=lchnk)
  else
    obj => null()
  end if

end function aerosol_state_object

!=========================================================================================
! returns a pointer to an aerosol properties object
function aerosol_properties_object() result(obj)

  class(aerosol_properties), pointer :: obj ! aerosol properties object pointer

  obj => aero_props_obj

end function aerosol_properties_object

!=========================================================================================

subroutine microp_aero_final

  integer :: c

  ! aerosol_instances_mod owns the props obj, so just nullify the pointer.
  nullify(aero_props_obj)

end subroutine microp_aero_final

!=========================================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_npccn_scale = unset_r8  ! prescribed aerosol bulk sulfur scale factor
   real(r8) :: microp_aero_wsub_scale = unset_r8   ! subgrid vertical velocity (liquid) scale factor
   real(r8) :: microp_aero_wsubi_scale = unset_r8  ! subgrid vertical velocity (ice) scale factor
   real(r8) :: microp_aero_wsub_min = unset_r8     ! subgrid vertical velocity (liquid) minimum (before scale factor)
   real(r8) :: microp_aero_wsub_min_asf = unset_r8 ! subgrid vertical velocity (liquid) minimum (after scale factor)
   real(r8) :: microp_aero_wsubi_min = unset_r8    ! subgrid vertical velocity (ice) minimum

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_npccn_scale, microp_aero_wsub_min, &
                             microp_aero_wsubi_min, microp_aero_wsub_scale, microp_aero_wsubi_scale, microp_aero_wsub_min_asf
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'microp_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, microp_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(microp_aero_npccn_scale, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: microp_aero_npccn_scale")
   call mpi_bcast(microp_aero_wsub_scale, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: microp_aero_wsub_scale")
   call mpi_bcast(microp_aero_wsubi_scale, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: microp_aero_wsubi_scale")
   call mpi_bcast(microp_aero_wsub_min, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: microp_aero_wsub_min")
   call mpi_bcast(microp_aero_wsub_min_asf, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: microp_aero_wsub_min_asf")
   call mpi_bcast(microp_aero_wsubi_min, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: microp_aero_wsubi_min")

   ! set local variables
   npccn_scale = microp_aero_npccn_scale
   wsub_scale = microp_aero_wsub_scale
   wsubi_scale = microp_aero_wsubi_scale
   wsub_min = microp_aero_wsub_min
   wsub_min_asf = microp_aero_wsub_min_asf
   wsubi_min = microp_aero_wsubi_min

   if(npccn_scale == unset_r8) call endrun(subname//": FATAL: npccn_scale is not set")
   if(wsub_scale == unset_r8) call endrun(subname//": FATAL: wsub_scale is not set")
   if(wsubi_scale == unset_r8) call endrun(subname//": FATAL: wsubi_scale is not set")
   if(wsub_min == unset_r8) call endrun(subname//": FATAL: wsub_min is not set")
   if(wsub_min_asf == unset_r8) call endrun(subname//": FATAL: wsub_min_asf is not set")
   if(wsubi_min == unset_r8) call endrun(subname//": FATAL: wsubi_min is not set")

   call nucleate_ice_cam_readnl(nlfile)
   call hetfrz_classnuc_cam_readnl(nlfile)

end subroutine microp_aero_readnl

!=========================================================================================

subroutine microp_aero_run ( &
   state, ptend_all, deltatin, pbuf)

   ! input arguments
   type(physics_state),         intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend_all
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! local workspace
   ! all units mks unless otherwise stated

   integer :: i, k, m
   integer :: itim_old

   type(physics_state), target :: state1                ! Local copy of state variable
   type(physics_ptend) :: ptend_loc

   real(r8), pointer :: ast(:,:)

   real(r8), pointer :: npccn(:,:)      ! number of CCN (liquid activated)

   real(r8), pointer :: rndst(:,:,:)    ! radius of 4 dust bins for contact freezing
   real(r8), pointer :: nacon(:,:,:)    ! number in 4 dust bins for contact freezing

   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl
   real(r8), pointer :: coarse_so4(:,:)  ! mass m.r. of coarse sulfate

   real(r8), pointer :: kvh(:,:)        ! vertical eddy diff coef (m2 s-1)
   real(r8), pointer :: tke(:,:)        ! TKE from the UW PBL scheme (m2 s-2)
   real(r8), pointer :: wp2(:,:)        ! CLUBB vertical velocity variance

   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cldo(:,:)       ! old cloud fraction

   real(r8), pointer :: dgnumwet(:,:,:) ! aerosol mode diameter

   real(r8) :: rho(pcols,pver)     ! air density (kg m-3)

   real(r8) :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: cldliqf(pcols,pver) ! fractional of total cloud that is liquid
   real(r8) :: qcld                ! total cloud water
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: dmc, ssmc, so4mc    ! variables for modal scheme.

   ! BAM diagnostic output from ndrop_bam_calc
   real(r8), allocatable :: ccn_bam(:,:,:)    ! CCN at 6 supersaturations (#/cm3)
   real(r8), allocatable :: naer2_bam(:,:,:)  ! aerosol number conc for diagnostics

   real(r8) :: wp2_full(pcols,pverp) ! CLUBB wp2 expanded onto the full interface grid (m2/s2)

   real(r8) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8) :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)

   real(r8) :: wght

   integer :: l, lchnk, ncol, astat

   real(r8), allocatable :: factnum(:,:,:) ! activation fraction for aerosol number

   class(aerosol_state), pointer :: aero_state1_obj
   type(aero_state_entry_t), allocatable :: aero_states1(:)
   integer :: nstates1, iaermod
   class(aerosol_properties), pointer :: props_tmp

   character(len=512)   :: errmsg
   integer              :: errflg

   !-------------------------------------------------------------------------------

   nullify(aero_state1_obj)

   call physics_state_copy(state,state1)

   lchnk = state1%lchnk
   ncol  = state1%ncol

   itim_old = pbuf_old_tim_idx()

   call pbuf_get_field(pbuf, npccn_idx, npccn)

   call pbuf_get_field(pbuf, nacon_idx, nacon)
   call pbuf_get_field(pbuf, rndst_idx, rndst)

   call physics_ptend_init(ptend_all, state%psetcols, 'microp_aero')

   !REMOVECAM: when microp_aero is brought into SIMA intermediate state1 updates should split into separate
   ! physics schemes, run tendency updaters, then the aerosol state is updated, so no need for factory pattern.
   ! create aerosol state objects via factory
   call aerosol_instances_create_states(list_idx=0, state=state1, pbuf=pbuf, aero_states=aero_states1, nstates=nstates1)
   !REMOVECAM_END

   ! find the appropriate state object for the active aerosol model
   do iaermod = 1, nstates1
      props_tmp => aerosol_instances_get_props(iaermod, 0)
      if (clim_modal_aero .and. props_tmp%model_is('modal')) then
         aero_state1_obj => aero_states1(iaermod)%obj
         exit
      else if (clim_carma_aero .and. props_tmp%model_is('CARMA')) then
         aero_state1_obj => aero_states1(iaermod)%obj
         exit
      else if (.not.clim_modal_aero .and. .not.clim_carma_aero .and. props_tmp%model_is('BAM')) then
         aero_state1_obj => aero_states1(iaermod)%obj
         exit
      end if
   end do

   if (clim_modal_aero.or.clim_carma_aero) then

      itim_old = pbuf_old_tim_idx()

      call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   end if

   if (clim_modal_aero) then
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet)
   end if

   ! initialize output
   npccn(1:ncol,1:pver)    = 0._r8

   nacon(1:ncol,1:pver,:)  = 0._r8

   ! set default or fixed dust bins for contact freezing
   rndst(1:ncol,1:pver,1) = rn_dst1
   rndst(1:ncol,1:pver,2) = rn_dst2
   rndst(1:ncol,1:pver,3) = rn_dst3
   rndst(1:ncol,1:pver,4) = rn_dst4

   ! initialize time-varying parameters
   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = state1%pmid(i,k)/(rair*state1%t(i,k))
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! More refined computation of sub-grid vertical velocity
   ! Set to be zero at the surface by initialization.
   wsub(:,:) = 0._r8
   select case (trim(eddy_scheme))
   case ('diag_TKE')
      call pbuf_get_field(pbuf, tke_idx, tke)
      call compute_subgrid_vertical_velocity_tke_run( &
           ncol    = ncol,    &
           pver    = pver,    &
           top_lev = top_lev, &
           tke     = tke(:ncol,:pverp), &
           wsub    = wsub(:ncol,:pver), &
           errmsg  = errmsg,  &
           errflg  = errflg)
      if(errflg /= 0) call endrun('compute_subgrid_vertical_velocity_tke_run: ' // errmsg)
   case ('CLUBB_SGS')
      itim_old = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, wp2_idx, wp2)

      ! The WP2_nadv pbuf field is dimensioned on the CLUBB momentum subgrid
      ! (nzm_clubb = pverp + 1 - top_lev), whose index 1 is CAM interface top_lev.
      ! The scheme expects wp2 on the full interface grid, so expand it here.
      wp2_full(:ncol, :top_lev-1)    = 0._r8
      wp2_full(:ncol, top_lev:pverp) = wp2(:ncol, 1:pverp-top_lev+1)

      call compute_subgrid_vertical_velocity_clubb_run( &
           ncol    = ncol,    &
           pver    = pver,    &
           pverp   = pverp,   &
           top_lev = top_lev, &
           wp2     = wp2_full(:ncol,:), &
           wsub    = wsub(:ncol,:pver), &
           errmsg  = errmsg,  &
           errflg  = errflg)
      if(errflg /= 0) call endrun('compute_subgrid_vertical_velocity_clubb_run: ' // errmsg)
   case default
      call pbuf_get_field(pbuf, kvh_idx, kvh)
      call compute_subgrid_vertical_velocity_kvh_run( &
           ncol    = ncol,    &
           pver    = pver,    &
           top_lev = top_lev, &
           kvh     = kvh(:ncol,:pverp), &
           wsub    = wsub(:ncol,:pver), &
           errmsg  = errmsg,  &
           errflg  = errflg)
      if(errflg /= 0) call endrun('compute_subgrid_vertical_velocity_kvh_run: ' // errmsg)
   end select

   ! Apply min/max/scale and derive wsubi (for ice nucleation)
   call scale_subgrid_vertical_velocity_run( &
        ncol                = ncol,                &
        pver                = pver,                &
        top_lev             = top_lev,             &
        wsub_min            = wsub_min,            &
        wsubi_min           = wsubi_min,           &
        wsub_scale          = wsub_scale,          &
        wsubi_scale         = wsubi_scale,         &
        use_preexisting_ice = use_preexisting_ice, &
        wsub                = wsub(:ncol,:pver),   &
        wsubi               = wsubi(:ncol,:pver),  &
        errmsg              = errmsg,              &
        errflg              = errflg)
   if(errflg /= 0) call endrun('scale_subgrid_vertical_velocity_run: ' // errmsg)

   call outfld('WSUB',   wsub, pcols, lchnk)
   call outfld('WSUBI', wsubi, pcols, lchnk)


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

   if (associated(aero_props_obj) .and. associated(aero_state1_obj)) then
      call nucleate_ice_cam_calc(state1, wsubi, pbuf, deltatin, ptend_loc, aero_props_obj, aero_state1_obj)
   else
      call nucleate_ice_cam_calc(state1, wsubi, pbuf, deltatin, ptend_loc)
   end if

   call physics_ptend_sum(ptend_loc, ptend_all, ncol)
   call physics_update(state1, ptend_loc, deltatin)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Droplet Activation
   if (clim_modal_aero .or. clim_carma_aero) then

      ! for modal or carma aerosol

      ! partition cloud fraction into liquid water part
      lcldn = 0._r8
      lcldo = 0._r8
      cldliqf = 0._r8
      do k = top_lev, pver
         do i = 1, ncol
            qcld = state1%q(i,k,cldliq_idx) + state1%q(i,k,cldice_idx)
            if (qcld > qsmall) then
               lcldn(i,k)   = cldn(i,k)*state1%q(i,k,cldliq_idx)/qcld
               lcldo(i,k)   = cldo(i,k)*state1%q(i,k,cldliq_idx)/qcld
               cldliqf(i,k) = state1%q(i,k,cldliq_idx)/qcld
            end if
         end do
      end do

      call outfld('LCLOUD', lcldn, pcols, lchnk)

      allocate(factnum(pcols,pver,aero_props_obj%nbins()),stat=astat)
      if (astat/=0) then
         call endrun('microp_aero_run: not able to allocate factnum')
      endif

      ! If not using preexsiting ice, then only use cloudbourne aerosol for the
      ! liquid clouds. This is the same behavior as CAM5.
      !
      ! ptend_loc is initialized inside dropmixnuc
      if (use_preexisting_ice) then
         call dropmixnuc( aero_props_obj, aero_state1_obj, &
              state1, ptend_loc, deltatin, pbuf, wsub, wsub_min_asf, &
              cldn, cldo, cldliqf, nctend_mixnuc, factnum)
      else
         cldliqf = 1._r8
         call dropmixnuc( aero_props_obj, aero_state1_obj, &
              state1, ptend_loc, deltatin, pbuf, wsub, wsub_min_asf, &
              lcldn, lcldo, cldliqf, nctend_mixnuc, factnum)
      end if

      npccn(:ncol,:) = nctend_mixnuc(:ncol,:)

      ! this scale is ONLY applied to modal/carma
      npccn(:ncol,:) = npccn(:ncol,:) * npccn_scale

   else

      ! ndrop_bam has no tendencies; this is initialized unconditionally in order
      ! for ptend_loc accummulation to happen the same as dropmixnuc above,
      ! including when no aerosols are present at all (aquaplanet compsets).
      call physics_ptend_init(ptend_loc, state1%psetcols, 'none')

      ! for bulk aerosol: activation, contact freezing, CCN diagnostics
      ! do not run for aquaplanet compsets which also gets in this path.
      if (associated(aero_state1_obj)) then
         ! get liquid cloud fraction. scaling is done within the ndrop_bam scheme.
         call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

         allocate(ccn_bam(pcols, pver, psat))
         allocate(naer2_bam(pcols, pver, naer_all))

         ! zero output to pcols
         npccn(:,:) = 0._r8
         nacon(:,:,:) = 0._r8
         ccn_bam(:,:,:) = 0._r8
         naer2_bam(:,:,:) = 0._r8

         ! Run CCPPized subroutine for droplet activation and contact freezing
         call ndrop_bam_calc( &
              aero_state = aero_state1_obj,                  &
              aero_props = aero_props_obj,                   &
              ncol       = ncol,                             &
              pver       = pver,                             &
              top_lev    = top_lev,                          &
              gravit     = gravit,                           &
              rair       = rair,                             &
              tmelt      = tmelt,                            &
              cpair      = cpair,                            &
              rh2o       = rh2o,                             &
              rhoh2o     = rhoh2o,                           &
              latvap     = latvap,                           &
              rho        = rho(:ncol,:),                     &
              tair       = state1%t(:ncol,:),                &
              wsub       = wsub(:ncol,:),                    &
              qcld       = state1%q(:ncol,:,cldliq_idx),     &
              qsmall_in  = qsmall,                           &
              ast        = ast(:ncol,:),                     &
              numliq     = state1%q(:ncol,:,numliq_idx),     &
              deltatin   = deltatin,                         &
              npccn      = npccn(:ncol,:),                   &
              nacon      = nacon(:ncol,:,:),                 &
              ccn        = ccn_bam(:ncol,:,:),               &
              naer2_diag = naer2_bam(:ncol,:,:),             &
              errmsg     = errmsg, &
              errflg     = errflg)
         if(errflg /= 0) then
            call endrun('ndrop_bam_calc: ' // errmsg)
         end if
      end if

   end if

   call physics_ptend_sum(ptend_loc, ptend_all, ncol)
   call physics_update(state1, ptend_loc, deltatin)


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
   ! estimate rndst and nacon for 4 dust bins here to pass to MG microphysics
   !
   ! For bulk aerosol contact freezing is handled inside ndrop_bam_calc
   ! (nacon is output from ndrop_bam_calc)
   ! Below for modal aerosol only:
   if(clim_modal_aero) then
      ! For modal aerosols:
      ! mode number mixing ratios
      call aero_state1_obj%get_ambient_num(mode_coarse_dst_idx, num_coarse)

      ! mode specie mass m.r.
      call aero_state1_obj%get_ambient_mmr(species_ndx=coarse_dust_idx, bin_ndx=mode_coarse_dst_idx, mmr=coarse_dust)
      call aero_state1_obj%get_ambient_mmr(species_ndx=coarse_nacl_idx, bin_ndx=mode_coarse_slt_idx, mmr=coarse_nacl)
      if (mode_coarse_idx>0) then
         call aero_state1_obj%get_ambient_mmr(species_ndx=coarse_so4_idx, bin_ndx=mode_coarse_idx, mmr=coarse_so4)
      endif


      !  use size '3' for dust coarse mode...
      !  scale by dust fraction in coarse mode
      do k = top_lev, pver
         do i = 1, ncol
            if (state1%t(i,k) < 269.15_r8) then
               dmc  = coarse_dust(i,k)
               ssmc = coarse_nacl(i,k)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting
                  wght = 1._r8
               else
                  so4mc = coarse_so4(i,k)
                  ! 3-mode -- needs weighting for dust since dust, seasalt, and sulfate  are combined in the "coarse" mode type
                  wght = dmc/(ssmc + dmc + so4mc)
               endif

               if (dmc > 0.0_r8) then
                  nacon(i,k,3) = wght*num_coarse(i,k)*rho(i,k)
               else
                  nacon(i,k,3) = 0._r8
               end if

               ! also redefine parameters based on size...
               rndst(i,k,3) = 0.5_r8*dgnumwet(i,k,mode_coarse_dst_idx)
               if (rndst(i,k,3) <= 0._r8) then
                  rndst(i,k,3) = rn_dst3
               end if
            end if
         end do
      end do
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! BAM diagnostic output (CCN concentration and aerosol number)
   if ((.not. clim_modal_aero) .and. (.not.clim_carma_aero) .and. &    ! not modal / carma
       allocated(ccn_bam) .and. allocated(naer2_bam)) then             ! and bulk aerosols are active (not aquap)
      do l = 1, psat
         call outfld(ccn_name(l), ccn_bam(1,1,l), pcols, lchnk)
      end do
      do l = 1, naer_all
         call outfld(trim(aername(l))//'_m3', naer2_bam(1,1,l), pcols, lchnk)
      end do
      deallocate(ccn_bam, naer2_bam)
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! heterogeneous freezing
   if (use_hetfrz_classnuc) then

      call hetfrz_classnuc_cam_calc(aero_props_obj, aero_state1_obj, state1, deltatin, factnum, pbuf)

   end if

   if (clim_modal_aero.or.clim_carma_aero) then
      deallocate(factnum)
   end if

   ! destroy all aerosol state objects created for this chunk
   nullify(aero_state1_obj)
   call aerosol_instances_destroy_states(aero_states1)

 end subroutine microp_aero_run

!=========================================================================================

end module microp_aero
