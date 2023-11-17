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
use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
use ref_pres,         only: top_lev => trop_cloud_top_lev
use physconst,        only: rair
use constituents,     only: cnst_get_ind
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init, physics_ptend_sum, &
                            physics_state_copy, physics_update
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                            pbuf_get_chunk
use phys_control,     only: phys_getopts, use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num

use nucleate_ice_cam, only: use_preexisting_ice, nucleate_ice_cam_readnl, nucleate_ice_cam_register, &
                            nucleate_ice_cam_init, nucleate_ice_cam_calc

use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn

use hetfrz_classnuc_cam, only: hetfrz_classnuc_cam_readnl, hetfrz_classnuc_cam_register, hetfrz_classnuc_cam_init, &
                               hetfrz_classnuc_cam_calc

use cam_history,      only: addfld, add_default, outfld
use cam_logfile,      only: iulog
use cam_abortutils,       only: endrun

use aerosol_properties_mod, only: aerosol_properties
use modal_aerosol_properties_mod, only: modal_aerosol_properties

use aerosol_state_mod, only: aerosol_state
use modal_aerosol_state_mod, only: modal_aerosol_state

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
real(r8) :: bulk_scale    ! prescribed aerosol bulk sulfur scale factor
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

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4

! modal aerosols
logical :: clim_modal_aero

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode
integer :: coarse_so4_idx = -1  ! index of sulfate in coarse mode

integer :: npccn_idx, rndst_idx, nacon_idx

logical  :: separate_dust = .false.

type aero_state_t
   class(aerosol_state), pointer :: obj=>null()
end type aero_state_t

class(aerosol_properties), pointer :: aero_props_obj=>null()
type(aero_state_t), pointer :: aero_state(:) => null()

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
   integer  :: iaer, ierr
   integer  :: m, n, nmodes, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer :: c

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
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ast_idx = pbuf_get_index('AST')

   if (clim_modal_aero) then

      cldo_idx = pbuf_get_index('CLDO')
      dgnumwet_idx = pbuf_get_index('DGNUMWET')

      aero_props_obj => modal_aerosol_properties()
      if (.not.associated(aero_props_obj)) then
         call endrun('ma_convproc_init: construction of modal_aerosol_properties object failed')
      end if
      call ndrop_init(aero_props_obj)
      call nucleate_ice_cam_init(mincld, bulk_scale, pbuf2d, aero_props=aero_props_obj)

      allocate(aero_state(begchunk:endchunk))
      do c = begchunk,endchunk
         pbuf => pbuf_get_chunk(pbuf2d, c)
         aero_state(c)%obj => modal_aerosol_state( phys_state(c), pbuf )
         if (.not.associated(aero_state(c)%obj)) then
            call endrun('microp_aero_init: construction of modal_aerosol_state object failed')
         end if
      end do

      ! Init indices for specific modes/species

      ! mode index for specified mode types
      do m = 1, nmodes
         call rad_cnst_get_info(0, m, mode_type=str32)
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
      call rad_cnst_get_info(0, mode_coarse_dst_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            coarse_dust_idx = n
         end select
      end do
      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do
      if (mode_coarse_idx>0) then
         call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
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

   else

      ! Props needed for BAM number concentration calcs.

      call rad_cnst_get_info(0, naero=naer_all)
      allocate( &
         aername(naer_all),        &
         num_to_mass_aer(naer_all) )

      do iaer = 1, naer_all
         call rad_cnst_get_aer_props(0, iaer, &
            aername         = aername(iaer), &
            num_to_mass_aer = num_to_mass_aer(iaer) )

         ! Look for sulfate, dust, and soot in this list (Bulk aerosol only)
         if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
      end do

      call ndrop_bam_init()
      call nucleate_ice_cam_init(mincld, bulk_scale, pbuf2d)

   end if

   call addfld('LCLOUD', (/ 'lev' /), 'A', ' ',   'Liquid cloud fraction used in stratus activation')

   call addfld('WSUB',   (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity'                   )
   call addfld('WSUBI',  (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity for ice'           )

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
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
function aerosol_state_object(lchnk) result(obj)

  integer,intent(in) :: lchnk ! local chunk index
  class(aerosol_state), pointer :: obj ! aerosol state object pointer for local chunk

  obj => aero_state(lchnk)%obj

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

  if (associated(aero_props_obj)) then
     deallocate(aero_props_obj)
  end if
  nullify(aero_props_obj)

  if (associated(aero_state)) then
     do c = begchunk,endchunk
        deallocate(aero_state(c)%obj)
     end do
     deallocate(aero_state)
     nullify(aero_state)
  end if

end subroutine microp_aero_final

!=========================================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_bulk_scale = unset_r8   ! prescribed aerosol bulk sulfur scale factor
   real(r8) :: microp_aero_npccn_scale = unset_r8  ! prescribed aerosol bulk sulfur scale factor
   real(r8) :: microp_aero_wsub_scale = unset_r8   ! subgrid vertical velocity (liquid) scale factor
   real(r8) :: microp_aero_wsubi_scale = unset_r8  ! subgrid vertical velocity (ice) scale factor
   real(r8) :: microp_aero_wsub_min = unset_r8     ! subgrid vertical velocity (liquid) minimum (before scale factor)
   real(r8) :: microp_aero_wsub_min_asf = unset_r8 ! subgrid vertical velocity (liquid) minimum (after scale factor)
   real(r8) :: microp_aero_wsubi_min = unset_r8    ! subgrid vertical velocity (ice) minimum

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_bulk_scale, microp_aero_npccn_scale, microp_aero_wsub_min, &
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
   call mpi_bcast(microp_aero_bulk_scale, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: microp_aero_bulk_scale")
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
   bulk_scale = microp_aero_bulk_scale
   npccn_scale = microp_aero_npccn_scale
   wsub_scale = microp_aero_wsub_scale
   wsubi_scale = microp_aero_wsubi_scale
   wsub_min = microp_aero_wsub_min
   wsub_min_asf = microp_aero_wsub_min_asf
   wsubi_min = microp_aero_wsubi_min

   if(bulk_scale == unset_r8) call endrun(subname//": FATAL: bulk_scale is not set")
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

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio

   real(r8) :: rho(pcols,pver)     ! air density (kg m-3)

   real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction

   real(r8) :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: cldliqf(pcols,pver) ! fractional of total cloud that is liquid
   real(r8) :: qcld                ! total cloud water
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: dum, dum2           ! temporary dummy variable
   real(r8) :: dmc, ssmc, so4mc    ! variables for modal scheme.

   ! bulk aerosol variables
   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8) :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)
   real(r8) :: nucboas

   real(r8) :: wght

   integer :: lchnk, ncol, astat

   real(r8), allocatable :: factnum(:,:,:) ! activation fraction for aerosol number

   class(aerosol_state), pointer :: aero_state1_obj

   !-------------------------------------------------------------------------------

   nullify(aero_state1_obj)

   call physics_state_copy(state,state1)

   lchnk = state1%lchnk
   ncol  = state1%ncol

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, npccn_idx, npccn)

   call pbuf_get_field(pbuf, nacon_idx, nacon)
   call pbuf_get_field(pbuf, rndst_idx, rndst)

   call physics_ptend_init(ptend_all, state%psetcols, 'microp_aero')

   if (clim_modal_aero) then
      ! create an aerosol state object specifically for cam state1
      aero_state1_obj => modal_aerosol_state( state1, pbuf )
      if (.not.associated(aero_state1_obj)) then
         call endrun('microp_aero_run: construction of aero_state1_obj modal_aerosol_state object failed')
      end if

      itim_old = pbuf_old_tim_idx()

      call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
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

   if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', state1, pbuf, num_coarse)

      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state1, pbuf, coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state1, pbuf, coarse_nacl)
      if (mode_coarse_idx>0) then
         call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_so4_idx, 'a', state1, pbuf, coarse_so4)
      endif

   else
      ! init number/mass arrays for bulk aerosols
      allocate( &
         naer2(pcols,pver,naer_all), &
         maerosol(pcols,pver,naer_all))

      do m = 1, naer_all
         call rad_cnst_get_aer_mmr(0, m, state1, pbuf, aer_mmr)
         maerosol(:ncol,:,m) = aer_mmr(:ncol,:)*rho(:ncol,:)

         if (m .eq. idxsul) then
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)*bulk_scale
         else
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)
         end if
      end do
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! More refined computation of sub-grid vertical velocity
   ! Set to be zero at the surface by initialization.

   select case (trim(eddy_scheme))
   case ('diag_TKE')
      call pbuf_get_field(pbuf, tke_idx, tke)
   case ('CLUBB_SGS')
      itim_old = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim_old/),kount=(/pcols,pverp,1/))
      allocate(tke(pcols,pverp))
      tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)

   case default
      call pbuf_get_field(pbuf, kvh_idx, kvh)
   end select

   ! Set minimum values above top_lev.
   wsub(:ncol,:top_lev-1)  = wsub_min
   wsubi(:ncol,:top_lev-1) = wsubi_min

   do k = top_lev, pver
      do i = 1, ncol

         select case (trim(eddy_scheme))
         case ('diag_TKE', 'CLUBB_SGS')
            wsub(i,k) = sqrt(0.5_r8*(tke(i,k) + tke(i,k+1))*(2._r8/3._r8))
            wsub(i,k) = min(wsub(i,k),10._r8)
         case default
            ! get sub-grid vertical velocity from diff coef.
            ! following morrison et al. 2005, JAS
            ! assume mixing length of 30 m
            dum = (kvh(i,k) + kvh(i,k+1))/2._r8/30._r8
            ! use maximum sub-grid vertical vel of 10 m/s
            dum = min(dum, 10._r8)
            ! set wsub to value at current vertical level
            wsub(i,k)  = dum
         end select

         wsubi(i,k) = max(wsubi_min, wsub(i,k)) * wsubi_scale
         if (.not. use_preexisting_ice) then
            wsubi(i,k) = min(wsubi(i,k), 0.2_r8)
         endif

         wsub(i,k)  = max(wsub_min, wsub(i,k)) * wsub_scale

      end do
   end do

   call outfld('WSUB',   wsub, pcols, lchnk)
   call outfld('WSUBI', wsubi, pcols, lchnk)

   if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

   if (associated(aero_props_obj).and.associated(aero_state1_obj)) then
      call nucleate_ice_cam_calc(state1, wsubi, pbuf, deltatin, ptend_loc, aero_props_obj, aero_state1_obj)
   else
      call nucleate_ice_cam_calc(state1, wsubi, pbuf, deltatin, ptend_loc)
   end if

   call physics_ptend_sum(ptend_loc, ptend_all, ncol)
   call physics_update(state1, ptend_loc, deltatin)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! get liquid cloud fraction, check for minimum

   do k = top_lev, pver
      do i = 1, ncol
         lcldm(i,k) = max(ast(i,k), mincld)
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Droplet Activation

   if (clim_modal_aero) then

      ! for modal aerosol

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

      npccn(:ncol,:) = npccn(:ncol,:) * npccn_scale

   else

      ! for bulk aerosol

      ! no tendencies returned from ndrop_bam_run, so just init ptend here
      call physics_ptend_init(ptend_loc, state1%psetcols, 'none')

      do k = top_lev, pver
         do i = 1, ncol

            if (state1%q(i,k,cldliq_idx) >= qsmall) then

               ! get droplet activation rate

               call ndrop_bam_run( &
                  wsub(i,k), state1%t(i,k), rho(i,k), naer2(i,k,:), naer_all, &
                  naer_all, maerosol(i,k,:),  &
                  dum2)
               dum = dum2
            else
               dum = 0._r8
            end if

            npccn(i,k) = (dum*lcldm(i,k) - state1%q(i,k,numliq_idx))/deltatin
         end do
      end do

   end if

   call physics_ptend_sum(ptend_loc, ptend_all, ncol)
   call physics_update(state1, ptend_loc, deltatin)


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
   ! estimate rndst and nanco for 4 dust bins here to pass to MG microphysics

   do k = top_lev, pver
      do i = 1, ncol

         if (state1%t(i,k) < 269.15_r8) then

            if (clim_modal_aero) then

               ! For modal aerosols:
               !  use size '3' for dust coarse mode...
               !  scale by dust fraction in coarse mode

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

               !also redefine parameters based on size...

               rndst(i,k,3) = 0.5_r8*dgnumwet(i,k,mode_coarse_dst_idx)
               if (rndst(i,k,3) <= 0._r8) then
                  rndst(i,k,3) = rn_dst3
               end if

            else

               !For Bulk Aerosols: set equal to aerosol number for dust for bins 2-4 (bin 1=0)

               if (idxdst2 > 0) then
                  nacon(i,k,2) = naer2(i,k,idxdst2)
               end if
               if (idxdst3 > 0) then
                  nacon(i,k,3) = naer2(i,k,idxdst3)
               end if
               if (idxdst4 > 0) then
                  nacon(i,k,4) = naer2(i,k,idxdst4)
               end if
            end if

         end if
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !bulk aerosol ccn concentration (modal does it in ndrop, from dropmixnuc)

   if (.not. clim_modal_aero) then

      ! ccn concentration as diagnostic
      call ndrop_bam_ccn(lchnk, ncol, maerosol, naer2)

      deallocate( &
         naer2,    &
         maerosol)

   end if

   ! heterogeneous freezing
   if (use_hetfrz_classnuc) then

      call hetfrz_classnuc_cam_calc(aero_props_obj, aero_state1_obj, state1, deltatin, factnum, pbuf)

   end if

   if (clim_modal_aero) then
      deallocate(factnum)
   end if

   if (associated(aero_state1_obj)) then
      ! destroy the aerosol state object
      deallocate(aero_state1_obj)
      nullify(aero_state1_obj)
   endif

 end subroutine microp_aero_run

!=========================================================================================

end module microp_aero
