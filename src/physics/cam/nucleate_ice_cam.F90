module nucleate_ice_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for nucleate_ice module.
!
!  B. Eaton - Sept 2014
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: pi, rair, tmelt
use constituents,   only: pcnst, cnst_get_ind
use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer, only: physics_buffer_desc
use phys_control,   only: use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field, &
                          pbuf_set_field
use cam_history,    only: addfld, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: qsat_water

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

use nucleate_ice,   only: nucleati_init, nucleati

use aerosol_properties_mod, only: aerosol_properties
use aerosol_state_mod, only: aerosol_state

use phys_control,   only: cam_physpkg_is

implicit none
private
save

public :: &
   nucleate_ice_cam_readnl,   &
   nucleate_ice_cam_register, &
   nucleate_ice_cam_init,     &
   nucleate_ice_cam_calc

! Namelist variables
logical, public, protected :: use_preexisting_ice = .false.
logical                    :: hist_preexisting_ice = .false.
logical                    :: nucleate_ice_incloud = .false.
logical                    :: nucleate_ice_use_troplev = .false.
real(r8)                   :: nucleate_ice_subgrid = -1._r8
real(r8)                   :: nucleate_ice_subgrid_strat = -1._r8
real(r8)                   :: nucleate_ice_strat = 0.0_r8

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction
real(r8) :: bulk_scale  ! prescribed aerosol bulk sulfur scale factor

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numice_idx = -1

integer :: &
   naai_idx = -1,     &
   naai_hom_idx = -1

integer :: &
   ast_idx   = -1

integer :: &
    qsatfac_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all = -1 ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst1  = -1 ! index in aerosol list for dust1
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4
integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHIL)

! modal aerosols
logical :: clim_modal_aero = .false.
logical :: prog_modal_aero = .false.

logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies

integer, allocatable :: aer_cnst_idx(:,:)

!===============================================================================
contains
!===============================================================================

subroutine nucleate_ice_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, masterprocid, mpi_logical, mpi_real8

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'nucleate_ice_cam_readnl'

  namelist /nucleate_ice_nl/ use_preexisting_ice, hist_preexisting_ice, &
       nucleate_ice_subgrid, nucleate_ice_subgrid_strat, nucleate_ice_strat, &
       nucleate_ice_incloud, nucleate_ice_use_troplev

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'nucleate_ice_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, nucleate_ice_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

  ! Broadcast namelist variables
  call mpi_bcast(use_preexisting_ice,  1, mpi_logical,masterprocid, mpicom, ierr)
  call mpi_bcast(hist_preexisting_ice, 1, mpi_logical,masterprocid, mpicom, ierr)
  call mpi_bcast(nucleate_ice_subgrid, 1, mpi_real8,  masterprocid, mpicom, ierr)
  call mpi_bcast(nucleate_ice_subgrid_strat, 1, mpi_real8,  masterprocid, mpicom, ierr)
  call mpi_bcast(nucleate_ice_strat,   1, mpi_real8,  masterprocid, mpicom, ierr)
  call mpi_bcast(nucleate_ice_incloud, 1, mpi_logical,masterprocid, mpicom, ierr)
  call mpi_bcast(nucleate_ice_use_troplev, 1, mpi_logical,masterprocid, mpicom, ierr)

end subroutine nucleate_ice_cam_readnl

!================================================================================================

subroutine nucleate_ice_cam_register()

   ! global scope for NAAI needed when clubb_do_icesuper=.true.
   call pbuf_add_field('NAAI',     'global', dtype_r8, (/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

end subroutine nucleate_ice_cam_register

!================================================================================================

subroutine nucleate_ice_cam_init(mincld_in, bulk_scale_in, pbuf2d, aero_props)
   use phys_control, only: phys_getopts
   use time_manager, only: is_first_step

   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: bulk_scale_in
   class(aerosol_properties), optional, intent(in) :: aero_props

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   integer :: iaer
   integer :: ierr
   integer :: ispc, ibin
   integer :: idxtmp
   integer :: nmodes

   character(len=*), parameter :: routine = 'nucleate_ice_cam_init'
   logical :: history_cesm_forcing

   character(len=32) :: tmpname

   !--------------------------------------------------------------------------------------------
   call phys_getopts(prog_modal_aero_out = prog_modal_aero, history_cesm_forcing_out = history_cesm_forcing)

   mincld     = mincld_in
   bulk_scale = bulk_scale_in

   lq(:) = .false.

   if (prog_modal_aero.and.use_preexisting_ice) then

      if (.not. present(aero_props)) then
         call endrun(routine//' :  aero_props must be present')
      end if

      ! constituent tendencies are calculated only if use_preexisting_ice is TRUE
      ! set lq for constituent tendencies --

      allocate(aer_cnst_idx(aero_props%nbins(),0:maxval(aero_props%nspecies())), stat=ierr)
      if( ierr /= 0 ) then
         call endrun(routine//': aer_cnst_idx allocation failed')
      end if
      aer_cnst_idx = -1

      do ibin = 1, aero_props%nbins()
         if (aero_props%icenuc_updates_num(ibin)) then

            ! constituents of this bin will need to be updated

            if (aero_props%icenuc_updates_mmr(ibin,0)) then ! species 0 indicates bin MMR
               call aero_props%amb_mmr_name( ibin, 0, tmpname)
            else
               call aero_props%amb_num_name( ibin, tmpname)
            end if

            call cnst_get_ind(tmpname, idxtmp, abort=.false.)
            aer_cnst_idx(ibin,0) = idxtmp
            if (idxtmp>0) then
               lq(idxtmp) = .true.
            end if

            ! iterate over the species within the bin
            do ispc = 1, aero_props%nspecies(ibin)
               if (aero_props%icenuc_updates_mmr(ibin,ispc)) then
                  ! this aerosol constituent will be updated
                  call aero_props%amb_mmr_name( ibin, ispc, tmpname)
                  call cnst_get_ind(tmpname, idxtmp, abort=.false.)
                  aer_cnst_idx(ibin,ispc) = idxtmp
                  if (idxtmp>0) then
                     lq(idxtmp) = .true.
                  end if
               end if
            end do

         end if
      end do

   end if

   ! Initialize naai.
   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, naai_idx, 0.0_r8)
   end if

   if( masterproc ) then
      write(iulog,*) 'nucleate_ice parameters:'
      write(iulog,*) '  mincld                     = ', mincld_in
      write(iulog,*) '  bulk_scale                 = ', bulk_scale_in
      write(iulog,*) '  use_preexisiting_ice       = ', use_preexisting_ice
      write(iulog,*) '  hist_preexisiting_ice      = ', hist_preexisting_ice
      write(iulog,*) '  nucleate_ice_subgrid       = ', nucleate_ice_subgrid
      write(iulog,*) '  nucleate_ice_subgrid_strat = ', nucleate_ice_subgrid_strat
      write(iulog,*) '  nucleate_ice_strat         = ', nucleate_ice_strat
      write(iulog,*) '  nucleate_ice_incloud       = ', nucleate_ice_incloud
      write(iulog,*) '  nucleate_ice_use_troplev   = ', nucleate_ice_use_troplev
   end if

   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMICE', numice_idx)
   qsatfac_idx  = pbuf_get_index('QSATFAC', ierr)

   if (((nucleate_ice_subgrid .eq. -1._r8) .or. (nucleate_ice_subgrid_strat .eq. -1._r8)) .and. (qsatfac_idx .eq. -1)) then
     call endrun(routine//': ERROR qsatfac is required when subgrid = -1 or subgrid_strat = -1')
   end if

   if (cam_physpkg_is("cam_dev")) then
      ! Updates for PUMAS v1.21+
      call addfld('NIHFTEN',  (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to homogenous freezing')
      call addfld('NIDEPTEN', (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to deposition nucleation')
      call addfld('NIIMMTEN', (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to immersion freezing')
      call addfld('NIMEYTEN', (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to meyers deposition')
   else
      call addfld('NIHF',  (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to homogenous freezing')
      call addfld('NIDEP', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to deposition nucleation')
      call addfld('NIIMM', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to immersion freezing')
      call addfld('NIMEY', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to meyers deposition')
   endif

   call addfld('NIREGM',(/ 'lev' /), 'A', 'C', 'Ice Nucleation Temperature Threshold for Regime')
   call addfld('NISUBGRID',(/ 'lev' /), 'A', '', 'Ice Nucleation subgrid saturation factor')
   call addfld('NITROP_PD',(/ 'lev' /), 'A', '', 'Chemical Tropopause probability')
   if ( history_cesm_forcing ) then
      call add_default('NITROP_PD',8,' ')
   endif

   if (use_preexisting_ice) then
      call addfld('fhom',      (/ 'lev' /), 'A','fraction', 'Fraction of cirrus where homogeneous freezing occur'   )
      call addfld ('WICE',     (/ 'lev' /), 'A','m/s','Vertical velocity Reduction caused by preexisting ice'  )
      call addfld ('WEFF',     (/ 'lev' /), 'A','m/s','Effective Vertical velocity for ice nucleation' )

      if (cam_physpkg_is("cam_dev")) then
         ! Updates for PUMAS v1.21+
         call addfld ('INnso4TEN',   (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency so4 (in) to ice_nucleation')
         call addfld ('INnbcTEN',    (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency bc  (in) to ice_nucleation')
         call addfld ('INndustTEN',  (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency dust (in) ice_nucleation')
         call addfld ('INondustTEN',  (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency dust (out) from ice_nucleation')
         call addfld ('INhetTEN',    (/ 'lev' /), 'A','1/m3/s', &
              'Tendency for contribution for in-cloud ice number density increase by het nucleation in ice cloud')
         call addfld ('INhomTEN',    (/ 'lev' /), 'A','1/m3/s', &
              'Tendency for contribution for in-cloud ice number density increase by hom nucleation in ice cloud')
      else
         call addfld ('INnso4',   (/ 'lev' /), 'A','1/m3','Number Concentration so4 (in) to ice_nucleation')
         call addfld ('INnbc',    (/ 'lev' /), 'A','1/m3','Number Concentration bc  (in) to ice_nucleation')
         call addfld ('INndust',  (/ 'lev' /), 'A','1/m3','Number Concentration dust (in) ice_nucleation')
         call addfld ('INondust',  (/ 'lev' /), 'A','1/m3','Number Concentration dust (out) from ice_nucleation')
         call addfld ('INhet',    (/ 'lev' /), 'A','1/m3', &
              'contribution for in-cloud ice number density increase by het nucleation in ice cloud')
         call addfld ('INhom',    (/ 'lev' /), 'A','1/m3', &
              'contribution for in-cloud ice number density increase by hom nucleation in ice cloud')
      endif

      call addfld ('INFrehom', (/ 'lev' /), 'A','frequency','hom IN frequency ice cloud')
      call addfld ('INFreIN',  (/ 'lev' /), 'A','frequency','frequency of ice nucleation occur')

      if (hist_preexisting_ice) then
         call add_default ('WSUBI   ', 1, ' ')  ! addfld/outfld calls are in microp_aero

         call add_default ('fhom    ', 1, ' ')
         call add_default ('WICE    ', 1, ' ')
         call add_default ('WEFF    ', 1, ' ')
         call add_default ('INnso4  ', 1, ' ')
         call add_default ('INnbc   ', 1, ' ')
         call add_default ('INndust ', 1, ' ')
         call add_default ('INhet   ', 1, ' ')
         call add_default ('INhom   ', 1, ' ')
         call add_default ('INFrehom', 1, ' ')
         call add_default ('INFreIN ', 1, ' ')
      end if
   end if

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)

   clim_modal_aero = (nmodes > 0)

   if (.not. clim_modal_aero) then

      ! Props needed for BAM number concentration calcs.

      call rad_cnst_get_info(0, naero=naer_all)
      allocate( &
         aername(naer_all),        &
         num_to_mass_aer(naer_all) )

      do iaer = 1, naer_all
         call rad_cnst_get_aer_props(0, iaer, &
            aername         = aername(iaer), &
            num_to_mass_aer = num_to_mass_aer(iaer))
         ! Look for sulfate, dust, and soot in this list (Bulk aerosol only)
         if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
         if (trim(aername(iaer)) == 'DUST1') idxdst1 = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
         if (trim(aername(iaer)) == 'BCPHI') idxbcphi = iaer
      end do
   end if


   call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc, nucleate_ice_incloud, iulog, pi, &
        mincld)

   ! get indices for fields in the physics buffer
   ast_idx = pbuf_get_index('AST')

end subroutine nucleate_ice_cam_init

!================================================================================================

subroutine nucleate_ice_cam_calc( &
   state, wsubi, pbuf, dtime, ptend, aero_props, aero_state )

   use tropopause,     only: tropopause_findChemTrop

   ! arguments
   type(physics_state), target, intent(in)    :: state
   real(r8),                    intent(in)    :: wsubi(:,:)
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   real(r8),                    intent(in)    :: dtime
   type(physics_ptend),         intent(out)   :: ptend
   class(aerosol_properties),optional, intent(in) :: aero_props
   class(aerosol_state),optional, intent(in) :: aero_state

   ! local workspace

   ! naai and naai_hom are the outputs shared with the microphysics
   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only)

   integer :: lchnk, ncol
   integer :: itim_old
   integer :: i, k, l, m

   character(len=32) :: spectype

   real(r8), pointer :: t(:,:)          ! input temperature (K)
   real(r8), pointer :: qn(:,:)         ! input water vapor mixing ratio (kg/kg)
   real(r8), pointer :: qc(:,:)         ! cloud water mixing ratio (kg/kg)
   real(r8), pointer :: qi(:,:)         ! cloud ice mixing ratio (kg/kg)
   real(r8), pointer :: ni(:,:)         ! cloud ice number conc (1/kg)
   real(r8), pointer :: pmid(:,:)       ! pressure at layer midpoints (pa)

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio

   real(r8), pointer :: ast(:,:)
   real(r8) :: icecldf(pcols,pver)  ! ice cloud fraction
   real(r8), pointer :: qsatfac(:,:)      ! Subgrid cloud water saturation scaling factor.

   real(r8) :: rho(pcols,pver)      ! air density (kg m-3)

   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
   real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
   real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water
   integer  :: troplev(pcols)       ! tropopause level

   real(r8) :: relhum(pcols,pver)  ! relative humidity
   real(r8) :: icldm(pcols,pver)   ! ice cloud fraction

   real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)
   real(r8) :: dso4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8) :: wght
   real(r8) :: oso4_num
   real(r8) :: odst_num
   real(r8) :: osoot_num
   real(r8) :: so4_num_st_cr_tot
   real(r8) :: ramp

   real(r8) :: subgrid(pcols,pver)
   real(r8) :: trop_pd(pcols,pver)

   ! For pre-existing ice
   real(r8) :: fhom(pcols,pver)    ! how much fraction of cloud can reach Shom
   real(r8) :: wice(pcols,pver)    ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
   real(r8) :: weff(pcols,pver)    ! effective Vertical velocity for ice nucleation (m/s); weff=wsubi-wice
   real(r8) :: INnso4(pcols,pver)   ! #/m3, so4 aerosol number used for ice nucleation
   real(r8) :: INnbc(pcols,pver)    ! #/m3, bc aerosol number used for ice nucleation
   real(r8) :: INndust(pcols,pver)  ! #/m3, dust aerosol number used for ice nucleation
   real(r8) :: INondust(pcols,pver)  ! #/m3, dust aerosol number used for ice nucleation
   real(r8) :: INhet(pcols,pver)    ! #/m3, ice number from het freezing
   real(r8) :: INhom(pcols,pver)    ! #/m3, ice number from hom freezing
   real(r8) :: INFrehom(pcols,pver) !  hom freezing occurence frequency.  1 occur, 0 not occur.
   real(r8) :: INFreIN(pcols,pver)  !  ice nucleation occerence frequency.   1 occur, 0 not occur.

   ! history output for ice nucleation
   real(r8) :: nihf(pcols,pver)  !output number conc of ice nuclei due to heterogenous freezing (1/m3)
   real(r8) :: niimm(pcols,pver) !output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
   real(r8) :: nidep(pcols,pver) !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
   real(r8) :: nimey(pcols,pver) !output number conc of ice nuclei due to meyers deposition (1/m3)
   real(r8) :: regm(pcols,pver)  !output temperature thershold for nucleation regime

   real(r8) :: size_wghts(pcols,pver)
   real(r8) :: type_wghts(pcols,pver)
   real(r8), pointer :: num_col(:,:)
   real(r8) :: dust_num_col(pcols,pver)
   real(r8) :: sulf_num_col(pcols,pver)
   real(r8) :: soot_num_col(pcols,pver)
   real(r8) :: sulf_num_tot_col(pcols,pver)

   integer :: idxtmp
   real(r8), pointer :: amb_num(:,:)
   real(r8), pointer :: amb_mmr(:,:)
   real(r8), pointer :: cld_num(:,:)
   real(r8), pointer :: cld_mmr(:,:)

   real(r8) :: delmmr, delmmr_sum
   real(r8) :: delnum, delnum_sum

   real(r8), parameter :: per_cm3 = 1.e-6_r8 ! factor for m-3 to cm-3 conversions

   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   t     => state%t
   qn    => state%q(:,:,1)
   qc    => state%q(:,:,cldliq_idx)
   qi    => state%q(:,:,cldice_idx)
   ni    => state%q(:,:,numice_idx)
   pmid  => state%pmid

   rho(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))

   if (clim_modal_aero) then

      call physics_ptend_init(ptend, state%psetcols, 'nucleatei', lq=lq)

   else
      ! init number/mass arrays for bulk aerosols
      allocate( &
           naer2(pcols,pver,naer_all), &
           maerosol(pcols,pver,naer_all))

      do m = 1, naer_all
         call rad_cnst_get_aer_mmr(0, m, state, pbuf, aer_mmr)
         maerosol(:ncol,:,m) = aer_mmr(:ncol,:)*rho(:ncol,:)

         if (m .eq. idxsul) then
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)*bulk_scale
         else
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)
         end if
      end do

      call physics_ptend_init(ptend, state%psetcols, 'nucleatei')
   end if

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   icecldf(:ncol,:pver) = ast(:ncol,:pver)

   ! naai and naai_hom are the outputs from this parameterization
   call pbuf_get_field(pbuf, naai_idx, naai)
   call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)
   naai(1:ncol,1:pver)     = 0._r8
   naai_hom(1:ncol,1:pver) = 0._r8

   ! Use the same criteria that is used in chemistry and in CLUBB (for cloud fraction)
   ! to determine whether to use tropospheric or stratospheric settings. Include the
   ! tropopause level so that the cold point tropopause will use the stratospheric values.
   call tropopause_findChemTrop(state, troplev)

   if ((nucleate_ice_subgrid .eq. -1._r8) .or. (nucleate_ice_subgrid_strat .eq. -1._r8)) then
      call pbuf_get_field(pbuf, qsatfac_idx, qsatfac)
   end if

   trop_pd(:,:) = 0._r8

   do k = top_lev, pver
      do i = 1, ncol
         trop_pd(i, troplev(i)) = 1._r8

         if (k <= troplev(i)) then
            if (nucleate_ice_subgrid_strat .eq. -1._r8) then
               subgrid(i, k) = 1._r8 / qsatfac(i, k)
            else
               subgrid(i, k) = nucleate_ice_subgrid_strat
            end if
         else
            if (nucleate_ice_subgrid .eq. -1._r8) then
               subgrid(i, k) = 1._r8 / qsatfac(i, k)
            else
               subgrid(i, k) = nucleate_ice_subgrid
            end if
         end if
      end do
   end do


   ! initialize history output fields for ice nucleation
   nihf(1:ncol,1:pver)  = 0._r8
   niimm(1:ncol,1:pver) = 0._r8
   nidep(1:ncol,1:pver) = 0._r8
   nimey(1:ncol,1:pver) = 0._r8

   regm(1:ncol,1:pver) = 0._r8

   if (use_preexisting_ice) then
      fhom(:,:)     = 0.0_r8
      wice(:,:)     = 0.0_r8
      weff(:,:)     = 0.0_r8
      INnso4(:,:)   = 0.0_r8
      INnbc(:,:)    = 0.0_r8
      INndust(:,:)  = 0.0_r8
      INondust(:,:)  = 0.0_r8
      INhet(:,:)    = 0.0_r8
      INhom(:,:)    = 0.0_r8
      INFrehom(:,:) = 0.0_r8
      INFreIN(:,:)  = 0.0_r8
   endif

   do k = top_lev, pver
      ! Get humidity and saturation vapor pressures
      call qsat_water(t(1:ncol,k), pmid(1:ncol,k), es(1:ncol), qs(1:ncol), ncol, gam=gammas(1:ncol))

      do i = 1, ncol

         relhum(i,k) = qn(i,k)/qs(i)

         ! get cloud fraction, check for minimum
         icldm(i,k) = max(icecldf(i,k), mincld)

      end do
   end do

   dust_num_col = 0._r8
   sulf_num_col = 0._r8
   sulf_num_tot_col = 0._r8
   soot_num_col = 0._r8

   if (clim_modal_aero) then

      if (.not.(present(aero_props).and.present(aero_state))) then
         call endrun('nucleate_ice_cam_calc: aero_props and aero_state must be present')
      end if

      ! collect number densities (#/cm^3) for dust, sulfate, and soot
      call aero_state%nuclice_get_numdens( aero_props, use_preexisting_ice, ncol, pver, rho, &
                                           dust_num_col, sulf_num_col, soot_num_col, sulf_num_tot_col )

   else
      ! for bulk model
      dust_num_col(:ncol,:) = naer2(:ncol,:,idxdst1)/25._r8 * per_cm3 & ! #/cm3
                            + naer2(:ncol,:,idxdst2)/25._r8 * per_cm3 &
                            + naer2(:ncol,:,idxdst3)/25._r8 * per_cm3 &
                            + naer2(:ncol,:,idxdst4)/25._r8 * per_cm3
      sulf_num_col(:ncol,:) = naer2(:ncol,:,idxsul)/25._r8 * per_cm3
      soot_num_col(:ncol,:) = naer2(:ncol,:,idxbcphi)/25._r8 * per_cm3
   endif

   kloop: do k = top_lev, pver
      iloop: do i = 1, ncol

         so4_num_st_cr_tot = 0._r8

         freezing: if (t(i,k) < tmelt - 5._r8) then

            ! set aerosol number for so4, soot, and dust with units #/cm^3
            so4_num = sulf_num_col(i,k)
            dst_num = dust_num_col(i,k)
            so4_num_st_cr_tot=sulf_num_tot_col(i,k)

            ! *** Turn off soot nucleation ***
            soot_num = 0.0_r8

            if (cam_physpkg_is("cam_dev")) then

               call nucleati( &
                    wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
                    qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
                    so4_num, dst_num, soot_num, subgrid(i,k),                 &
                    naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
                    wice(i,k), weff(i,k), fhom(i,k), regm(i,k),               &
                    oso4_num, odst_num, osoot_num, &
                    call_frm_zm_in = .false., add_preexisting_ice_in = .false.)

            else

               call nucleati( &
                    wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
                    qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
                    so4_num, dst_num, soot_num, subgrid(i,k),                 &
                    naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
                    wice(i,k), weff(i,k), fhom(i,k), regm(i,k),               &
                    oso4_num, odst_num, osoot_num)

            end if

            ! Move aerosol used for nucleation from interstial to cloudborne,
            ! otherwise the same coarse mode aerosols will be available again
            ! in the next timestep and will supress homogeneous freezing.


            if (prog_modal_aero .and. use_preexisting_ice) then

               ! compute tendencies for transported aerosol constituents
               ! and update not-transported constituents

               do m = 1, aero_props%nbins()

                  if (aero_props%icenuc_updates_num(m)) then

                     ! constituents of this bin will need to be updated

                     call aero_state%get_ambient_num(m, amb_num)
                     call aero_state%get_cldbrne_num(m, cld_num)

                     if (amb_num(i,k)>0._r8) then
                        delmmr_sum = 0._r8
                        delnum_sum = 0._r8

                        ! iterate over the species within the bin
                        do l = 1, aero_props%nspecies(m)
                           if (aero_props%icenuc_updates_mmr(m,l)) then

                              call aero_props%species_type(m, l, spectype)
                              call aero_state%icenuc_size_wght( m, i,k, spectype, use_preexisting_ice, wght)

                              if (wght>0._r8) then

                                 ! this aerosol constituent will be updated

                                 idxtmp = aer_cnst_idx(m,l)

                                 call aero_state%get_ambient_mmr(l,m,amb_mmr)
                                 call aero_state%get_cldbrne_mmr(l,m,cld_mmr)

                                 ! determine change in aerosol mass
                                 delmmr = 0._r8
                                 delnum = 0._r8
                                 if (trim(spectype)=='dust') then
                                    if (dst_num>0._r8) then
                                       delmmr = (odst_num / dst_num) * icldm(i,k) * amb_mmr(i,k) * wght
                                       delnum = (odst_num * icldm(i,k)) /rho(i,k)/per_cm3
                                    endif
                                 elseif (trim(spectype)=='sulfate') then
                                    if (so4_num>0._r8) then
                                       delmmr = (oso4_num / so4_num) * icldm(i,k) * amb_mmr(i,k) * wght
                                       delnum = (oso4_num * icldm(i,k)) /rho(i,k)/per_cm3
                                    endif
                                 endif

                                 if (idxtmp>0) then
                                    ! constituent tendency (for transported species)
                                    ptend%q(i,k,idxtmp) = -delmmr/dtime
                                 else
                                    ! apply change of mass to not-transported species
                                    amb_mmr(i,k) = amb_mmr(i,k) - delmmr
                                 endif
                                 cld_mmr(i,k) = cld_mmr(i,k) + delmmr

                                 delmmr_sum = delmmr_sum + delmmr
                                 delnum_sum = delnum_sum + delnum
                              end if
                           end if
                        end do

                        idxtmp = aer_cnst_idx(m,0)

                        ! update aerosol state bin and tendency for grid box i,k
                        call aero_state%update_bin( m,i,k, delmmr_sum, delnum_sum, idxtmp, dtime, ptend%q )

                     end if

                  end if
               end do

            end if


            ! Liu&Penner does not generate enough nucleation in the polar winter
            ! stratosphere, which affects surface area density, dehydration and
            ! ozone chemistry. Part of this is that there are a larger number of
            ! particles in the accumulation mode than in the Aitken mode. In volcanic
            ! periods, the coarse mode may also be important. As a short
            ! term work around, include the accumulation and coarse mode particles
            ! and assume a larger fraction of the sulfates nucleate in the polar
            ! stratosphere.
            !
            ! Do not include the tropopause level, as stratospheric aerosols
            ! only exist above the tropopause level.
            !
            ! NOTE: This may still not represent the proper particles that
            ! participate in nucleation, because it doesn't include STS and NAT
            ! particles. It may not represent the proper saturation threshold for
            ! nucleation, and wsubi from CLUBB is probably not representative of
            ! wave driven varaibility in the polar stratosphere.
            if (nucleate_ice_use_troplev .and. clim_modal_aero) then
               if ((k < troplev(i)) .and. (nucleate_ice_strat > 0._r8) .and. (oso4_num > 0._r8)) then
                  dso4_num = max(0._r8, (nucleate_ice_strat*so4_num_st_cr_tot - oso4_num) * 1e6_r8 / rho(i,k))
                  naai(i,k) = naai(i,k) + dso4_num
                  nihf(i,k) = nihf(i,k) + dso4_num
               endif
            else
               ! This maintains backwards compatibility with the previous version.
               if (pmid(i,k) <= 12500._r8 .and. pmid(i,k) > 100._r8 .and. abs(state%lat(i)) >= 60._r8 * pi / 180._r8) then
                  ramp = 1._r8 - min(1._r8, max(0._r8, (pmid(i,k) - 10000._r8) / 2500._r8))

                  if (oso4_num > 0._r8) then
                     dso4_num = (max(oso4_num, ramp * nucleate_ice_strat * so4_num) - oso4_num) * 1e6_r8 / rho(i,k)
                     naai(i,k) = naai(i,k) + dso4_num
                     nihf(i,k) = nihf(i,k) + dso4_num
                  end if
               end if
            end if

            if (cam_physpkg_is("cam_dev")) then
               !Updates for pumas v1.21+

               naai_hom(i,k) = nihf(i,k)/dtime
               naai(i,k)= naai(i,k)/dtime

               ! output activated ice (convert from #/kg -> #/m3/s)
               nihf(i,k)     = nihf(i,k) *rho(i,k)/dtime
               niimm(i,k)    = niimm(i,k)*rho(i,k)/dtime
               nidep(i,k)    = nidep(i,k)*rho(i,k)/dtime
               nimey(i,k)    = nimey(i,k)*rho(i,k)/dtime

               if (use_preexisting_ice) then
                  INnso4(i,k) =so4_num*1e6_r8/dtime  ! (convert from #/cm3 -> #/m3/s)
                  INnbc(i,k)  =soot_num*1e6_r8/dtime
                  INndust(i,k)=dst_num*1e6_r8/dtime
                  INondust(i,k)=odst_num*1e6_r8/dtime
                  INFreIN(i,k)=1.0_r8          ! 1,ice nucleation occur
                  INhet(i,k) = (niimm(i,k) + nidep(i,k))   ! #/m3/s, nimey not in cirrus
                  INhom(i,k) = nihf(i,k)                 ! #/m3/s
                  if (INhom(i,k).gt.1e3_r8)   then ! > 1/L
                     INFrehom(i,k)=1.0_r8       ! 1, hom freezing occur
                  endif

                  ! exclude  no ice nucleaton
                  if ((INFrehom(i,k) < 0.5_r8) .and. (INhet(i,k) < 1.0_r8))   then
                     INnso4(i,k) =0.0_r8
                     INnbc(i,k)  =0.0_r8
                     INndust(i,k)=0.0_r8
                     INondust(i,k)=0.0_r8
                     INFreIN(i,k)=0.0_r8
                     INhet(i,k) = 0.0_r8
                     INhom(i,k) = 0.0_r8
                     INFrehom(i,k)=0.0_r8
                     wice(i,k) = 0.0_r8
                     weff(i,k) = 0.0_r8
                     fhom(i,k) = 0.0_r8
                  endif
               endif

            else ! Not cam_dev

               naai_hom(i,k) = nihf(i,k)

               ! output activated ice (convert from #/kg -> #/m3/s)
               nihf(i,k)     = nihf(i,k) *rho(i,k)
               niimm(i,k)    = niimm(i,k)*rho(i,k)
               nidep(i,k)    = nidep(i,k)*rho(i,k)
               nimey(i,k)    = nimey(i,k)*rho(i,k)

               if (use_preexisting_ice) then
                  INnso4(i,k) =so4_num*1e6_r8 ! (convert from #/cm3 -> #/m3/s)
                  INnbc(i,k)  =soot_num*1e6_r8
                  INndust(i,k)=dst_num*1e6_r8
                  INondust(i,k)=odst_num*1e6_r8
                  INFreIN(i,k)=1.0_r8          ! 1,ice nucleation occur
                  INhet(i,k) = (niimm(i,k) + nidep(i,k))   ! #/m3, nimey not in cirrus
                  INhom(i,k) = nihf(i,k)                 ! #/m3
                  if (INhom(i,k).gt.1e3_r8)   then ! > 1/L
                     INFrehom(i,k)=1.0_r8       ! 1, hom freezing occur
                  endif

                  ! exclude  no ice nucleaton
                  if ((INFrehom(i,k) < 0.5_r8) .and. (INhet(i,k) < 1.0_r8))   then
                     INnso4(i,k) =0.0_r8
                     INnbc(i,k)  =0.0_r8
                     INndust(i,k)=0.0_r8
                     INondust(i,k)=0.0_r8
                     INFreIN(i,k)=0.0_r8
                     INhet(i,k) = 0.0_r8
                     INhom(i,k) = 0.0_r8
                     INFrehom(i,k)=0.0_r8
                     wice(i,k) = 0.0_r8
                     weff(i,k) = 0.0_r8
                     fhom(i,k) = 0.0_r8
                  endif
               end if

            end if ! cam_dev
         end if freezing
      end do iloop
   end do kloop

   if (.not. clim_modal_aero) then
      deallocate( &
           naer2, &
           maerosol)
   end if

   if (cam_physpkg_is("cam_dev")) then
      ! Updates for PUMAS v1.21+
      call outfld('NIHFTEN',   nihf, pcols, lchnk)
      call outfld('NIIMMTEN', niimm, pcols, lchnk)
      call outfld('NIDEPTEN', nidep, pcols, lchnk)
      call outfld('NIMEYTEN', nimey, pcols, lchnk)
   else
      call outfld('NIHF',   nihf, pcols, lchnk)
      call outfld('NIIMM', niimm, pcols, lchnk)
      call outfld('NIDEP', nidep, pcols, lchnk)
      call outfld('NIMEY', nimey, pcols, lchnk)
   end if
   call outfld('NIREGM', regm, pcols, lchnk)
   call outfld('NISUBGRID', subgrid, pcols, lchnk)
   call outfld('NITROP_PD', trop_pd, pcols, lchnk)

   if (use_preexisting_ice) then
      call outfld( 'fhom' , fhom, pcols, lchnk)
      call outfld( 'WICE' , wice, pcols, lchnk)
      call outfld( 'WEFF' , weff, pcols, lchnk)
      if (cam_physpkg_is("cam_dev")) then
         ! Updates for PUMAS v1.21+
         call outfld('INnso4TEN',INnso4 , pcols,lchnk)
         call outfld('INnbcTEN',INnbc  , pcols,lchnk)
         call outfld('INndustTEN',INndust, pcols,lchnk)
         call outfld('INondustTEN',INondust, pcols,lchnk)
         call outfld('INhetTEN',INhet  , pcols,lchnk)
         call outfld('INhomTEN',INhom  , pcols,lchnk)
      else
         call outfld('INnso4  ',INnso4 , pcols,lchnk)
         call outfld('INnbc   ',INnbc  , pcols,lchnk)
         call outfld('INndust ',INndust, pcols,lchnk)
         call outfld('INondust ',INondust, pcols,lchnk)
         call outfld('INhet   ',INhet  , pcols,lchnk)
         call outfld('INhom   ',INhom  , pcols,lchnk)
      end if
      call outfld('INFrehom',INFrehom,pcols,lchnk)
      call outfld('INFreIN ',INFreIN, pcols,lchnk)
   end if

end subroutine nucleate_ice_cam_calc

!================================================================================================

end module nucleate_ice_cam
