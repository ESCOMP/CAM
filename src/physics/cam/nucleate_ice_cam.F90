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
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props, rad_cnst_get_mode_num_idx, &
                            rad_cnst_get_mam_mmr_idx

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: qsat_water, svp_water, svp_ice
use shr_spfn_mod,   only: erf => shr_spfn_erf

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

use nucleate_ice,   only: nucleati_init, nucleati


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
   naai_idx,     &
   naai_hom_idx

integer :: &
   ast_idx   = -1, &
   dgnum_idx = -1

integer :: &
    qsatfac_idx
    
! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst1  = -1 ! index in aerosol list for dust1
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4
integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHIL)

! modal aerosols
logical :: clim_modal_aero
logical :: prog_modal_aero

integer :: nmodes = -1
integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode
integer :: coarse_so4_idx = -1   ! index of sulfate in coarse mode

logical  :: separate_dust = .false.
real(r8) :: sigmag_aitken
real(r8) :: sigmag_accum

logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
integer :: cnum_idx, cdst_idx, cso4_idx

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

   call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

end subroutine nucleate_ice_cam_register

!================================================================================================

subroutine nucleate_ice_cam_init(mincld_in, bulk_scale_in)
   use phys_control, only: phys_getopts

   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: bulk_scale_in

   ! local variables
   integer  :: iaer
   integer :: ierr
   integer  :: m, n, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'nucleate_ice_cam_init'
   logical :: history_cesm_forcing
   !--------------------------------------------------------------------------------------------
   call phys_getopts(prog_modal_aero_out = prog_modal_aero, history_cesm_forcing_out = history_cesm_forcing)

   mincld     = mincld_in
   bulk_scale = bulk_scale_in

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
   
   call addfld('NIHF',  (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to homogenous freezing')
   call addfld('NIDEP', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to deposition nucleation')
   call addfld('NIIMM', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to immersion freezing')
   call addfld('NIMEY', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to meyers deposition')

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
      call addfld ('INnso4',   (/ 'lev' /), 'A','1/m3','Number Concentation so4 (in) to ice_nucleation')
      call addfld ('INnbc',    (/ 'lev' /), 'A','1/m3','Number Concentation bc  (in) to ice_nucleation')
      call addfld ('INndust',  (/ 'lev' /), 'A','1/m3','Number Concentation dust (in) ice_nucleation')
      call addfld ('INondust',  (/ 'lev' /), 'A','1/m3','Number Concentation dust (out) from ice_nucleation')
      call addfld ('INhet',    (/ 'lev' /), 'A','1/m3', &
                'contribution for in-cloud ice number density increase by het nucleation in ice cloud')
      call addfld ('INhom',    (/ 'lev' /), 'A','1/m3', &
                'contribution for in-cloud ice number density increase by hom nucleation in ice cloud')
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

   if (clim_modal_aero) then

      dgnum_idx    = pbuf_get_index('DGNUM' )

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
      if (mode_coarse_dst_idx < 0) mode_coarse_dst_idx = mode_coarse_idx
      if (mode_coarse_slt_idx < 0) mode_coarse_slt_idx = mode_coarse_idx

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! species indices for specified types
      ! find indices for the dust, seasalt and sulfate species in the coarse mode
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


      ! get specific mode properties
      call rad_cnst_get_mode_props(0, mode_aitken_idx, sigmag=sigmag_aitken)
      call rad_cnst_get_mode_props(0, mode_accum_idx, sigmag=sigmag_accum)

      if (prog_modal_aero) then
         call rad_cnst_get_mode_num_idx(mode_coarse_dst_idx, cnum_idx)
         call rad_cnst_get_mam_mmr_idx(mode_coarse_dst_idx, coarse_dust_idx, cdst_idx)
         if (mode_coarse_idx>0) then
            call rad_cnst_get_mam_mmr_idx(mode_coarse_idx, coarse_so4_idx, cso4_idx)
         end if
         lq(cnum_idx) = .true.
         lq(cdst_idx) = .true.
      endif

   else

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
         if (trim(aername(iaer)) == 'BCPHIL') idxbcphi = iaer
      end do
   end if


   call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc, nucleate_ice_incloud, iulog, pi, &
        mincld)

   ! get indices for fields in the physics buffer
   ast_idx      = pbuf_get_index('AST')

end subroutine nucleate_ice_cam_init

!================================================================================================

subroutine nucleate_ice_cam_calc( &
   state, wsubi, pbuf, dtime, ptend)

   use tropopause,     only: tropopause_findChemTrop

   ! arguments
   type(physics_state), target, intent(in)    :: state
   real(r8),                    intent(in)    :: wsubi(:,:)
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   real(r8),                    intent(in)    :: dtime
   type(physics_ptend),         intent(out)   :: ptend
 
   ! local workspace

   ! naai and naai_hom are the outputs shared with the microphysics
   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation 
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only)

   integer :: lchnk, ncol
   integer :: itim_old
   integer :: i, k, m

   real(r8), pointer :: t(:,:)          ! input temperature (K)
   real(r8), pointer :: qn(:,:)         ! input water vapor mixing ratio (kg/kg)
   real(r8), pointer :: qc(:,:)         ! cloud water mixing ratio (kg/kg)
   real(r8), pointer :: qi(:,:)         ! cloud ice mixing ratio (kg/kg)
   real(r8), pointer :: ni(:,:)         ! cloud ice number conc (1/kg)
   real(r8), pointer :: pmid(:,:)       ! pressure at layer midpoints (pa)

   real(r8), pointer :: num_accum(:,:)  ! number m.r. of accumulation mode
   real(r8), pointer :: num_aitken(:,:) ! number m.r. of aitken mode
   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl
   real(r8), pointer :: coarse_so4(:,:) ! mass m.r. of coarse sulfate
   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio
   real(r8), pointer :: dgnum(:,:,:)    ! mode dry radius
   real(r8), pointer :: cld_num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: cld_coarse_dust(:,:) ! mass m.r. of coarse dust

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

   real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
   real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)
   real(r8) :: wght
   real(r8) :: dmc
   real(r8) :: ssmc
   real(r8) :: so4mc
   real(r8) :: oso4_num
   real(r8) :: odst_num
   real(r8) :: osoot_num
   real(r8) :: dso4_num
   real(r8) :: so4_num_ac
   real(r8) :: so4_num_cr
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


   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   t     => state%t
   qn    => state%q(:,:,1)
   qc    => state%q(:,:,cldliq_idx)
   qi    => state%q(:,:,cldice_idx)
   ni    => state%q(:,:,numice_idx)
   pmid  => state%pmid

   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

   if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_accum_idx,  'a', state, pbuf, num_accum)
      call rad_cnst_get_mode_num(0, mode_aitken_idx, 'a', state, pbuf, num_aitken)
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', state, pbuf, num_coarse)

      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state, pbuf, coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state, pbuf, coarse_nacl)
      if (mode_coarse_idx>0) then
         call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_so4_idx, 'a', state, pbuf, coarse_so4)
      endif

      ! Get the cloudbourne coarse mode fields, so aerosol used for nucleated
      ! can be moved from interstial to cloudbourne.
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'c', state, pbuf, cld_num_coarse)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'c', state, pbuf, cld_coarse_dust)

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

   if (clim_modal_aero) then
      call pbuf_get_field(pbuf, dgnum_idx, dgnum)
   end if

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
      call qsat_water(t(:ncol,k), pmid(:ncol,k), &
           es(:ncol), qs(:ncol), gam=gammas(:ncol))

      do i = 1, ncol

         relhum(i,k) = qn(i,k)/qs(i)

         ! get cloud fraction, check for minimum
         icldm(i,k) = max(icecldf(i,k), mincld)

      end do
   end do


   do k = top_lev, pver
      do i = 1, ncol

         if (t(i,k) < tmelt - 5._r8) then

            ! compute aerosol number for so4, soot, and dust with units #/cm^3
            so4_num  = 0._r8
            soot_num = 0._r8
            dst1_num = 0._r8
            dst2_num = 0._r8
            dst3_num = 0._r8
            dst4_num = 0._r8
            dst_num  = 0._r8
            so4_num_cr = 0._r8

            if (clim_modal_aero) then
               !For modal aerosols, assume for the upper troposphere:
               ! soot = accumulation mode
               ! sulfate = aiken mode
               ! dust = coarse mode
               ! since modal has internal mixtures.
               soot_num = num_accum(i,k)*rho(i,k)*1.0e-6_r8
               dmc  = coarse_dust(i,k)*rho(i,k)
               ssmc = coarse_nacl(i,k)*rho(i,k)

               if (dmc > 0._r8) then
                  if ( separate_dust ) then
                     ! 7-mode -- has separate dust and seasalt mode types and
                     !           no need for weighting 
                     wght = 1._r8
                  else
                     ! 3-mode -- needs weighting for dust since dust, seasalt,
                     !           and sulfate are combined in the "coarse" mode type
                     so4mc    = coarse_so4(i,k)*rho(i,k)
                     wght = dmc/(ssmc + dmc + so4mc)
                  endif
                  dst_num = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
               else 
                  dst_num = 0.0_r8
               end if

               if ( separate_dust ) then
                  ! 7-mode -- the 7 mode scheme does not support
                  ! stratospheric sulfates, and the sulfates are mixed in
                  ! with the separate soot and dust modes, so just ignore
                  ! for now.
                  so4_num_cr = 0.0_r8
               else
                  ! 3-mode -- needs weighting for dust since dust, seasalt,
                  !           and sulfate are combined in the "coarse" mode
                  !           type
                  so4mc    = coarse_so4(i,k)*rho(i,k)
                  
                  if (so4mc > 0._r8) then
                    wght = so4mc/(ssmc + dmc + so4mc)
                    so4_num_cr = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
                  else
                    so4_num_cr = 0.0_r8
                  end if
               endif

               so4_num = 0.0_r8 
               if (.not. use_preexisting_ice) then
                  if (dgnum(i,k,mode_aitken_idx) > 0._r8) then
                     ! only allow so4 with D>0.1 um in ice nucleation
                     so4_num = so4_num + max(0._r8, num_aitken(i,k)*rho(i,k)*1.0e-6_r8 &
                        * (0.5_r8 - 0.5_r8*erf(log(0.1e-6_r8/dgnum(i,k,mode_aitken_idx))/  &
                        (2._r8**0.5_r8*log(sigmag_aitken)))))
                  end if
               else
                  ! all so4 from aitken
                  so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8
               end if

            else

               if (idxsul > 0) then 
                  so4_num = naer2(i,k,idxsul)/25._r8 *1.0e-6_r8
               end if
               if (idxbcphi > 0) then 
                  soot_num = naer2(i,k,idxbcphi)/25._r8 *1.0e-6_r8
               end if
               if (idxdst1 > 0) then 
                  dst1_num = naer2(i,k,idxdst1)/25._r8 *1.0e-6_r8
               end if
               if (idxdst2 > 0) then 
                  dst2_num = naer2(i,k,idxdst2)/25._r8 *1.0e-6_r8
               end if
               if (idxdst3 > 0) then 
                  dst3_num = naer2(i,k,idxdst3)/25._r8 *1.0e-6_r8
               end if
               if (idxdst4 > 0) then 
                  dst4_num = naer2(i,k,idxdst4)/25._r8 *1.0e-6_r8
               end if
               dst_num = dst1_num + dst2_num + dst3_num + dst4_num

            end if

            ! *** Turn off soot nucleation ***
            soot_num = 0.0_r8

            call nucleati( &
               wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
               qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
               so4_num, dst_num, soot_num, subgrid(i,k),                 &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
               wice(i,k), weff(i,k), fhom(i,k), regm(i,k),               &
               oso4_num, odst_num, osoot_num)

            ! Move aerosol used for nucleation from interstial to cloudborne, 
            ! otherwise the same coarse mode aerosols will be available again
            ! in the next timestep and will supress homogeneous freezing.
            if (prog_modal_aero .and. use_preexisting_ice) then
               if (separate_dust) then
                  call endrun('nucleate_ice_cam: use_preexisting_ice is not supported in separate_dust mode (MAM7)')
               endif
               ptend%q(i,k,cnum_idx) = -(odst_num * icldm(i,k))/rho(i,k)/1e-6_r8/dtime
               cld_num_coarse(i,k)   = cld_num_coarse(i,k) + (odst_num * icldm(i,k))/rho(i,k)/1e-6_r8

               ptend%q(i,k,cdst_idx) = - odst_num / dst_num * icldm(i,k) * coarse_dust(i,k) / dtime
               cld_coarse_dust(i,k) = cld_coarse_dust(i,k) + odst_num / dst_num *icldm(i,k) * coarse_dust(i,k)
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
              if ((k < troplev(i)) .and. (nucleate_ice_strat > 0._r8)) then
                 if (oso4_num > 0._r8) then
                    so4_num_ac = num_accum(i,k)*rho(i,k)*1.0e-6_r8
                    dso4_num = max(0._r8, (nucleate_ice_strat * (so4_num_cr + so4_num_ac)) - oso4_num) * 1e6_r8 / rho(i,k)
                    naai(i,k) = naai(i,k) + dso4_num
                    nihf(i,k) = nihf(i,k) + dso4_num
                 end if
              end if
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

            naai_hom(i,k) = nihf(i,k)

            ! output activated ice (convert from #/kg -> #/m3)
            nihf(i,k)     = nihf(i,k) *rho(i,k)
            niimm(i,k)    = niimm(i,k)*rho(i,k)
            nidep(i,k)    = nidep(i,k)*rho(i,k)
            nimey(i,k)    = nimey(i,k)*rho(i,k)

            if (use_preexisting_ice) then
               INnso4(i,k) =so4_num*1e6_r8  ! (convert from #/cm3 -> #/m3)
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

         end if
      end do
   end do

   if (.not. clim_modal_aero) then

      deallocate( &
         naer2,    &
         maerosol)

   end if

   call outfld('NIHF',   nihf, pcols, lchnk)
   call outfld('NIIMM', niimm, pcols, lchnk)
   call outfld('NIDEP', nidep, pcols, lchnk)
   call outfld('NIMEY', nimey, pcols, lchnk)
   call outfld('NIREGM', regm, pcols, lchnk)
   call outfld('NISUBGRID', subgrid, pcols, lchnk)
   call outfld('NITROP_PD', trop_pd, pcols, lchnk)

   if (use_preexisting_ice) then
      call outfld( 'fhom' , fhom, pcols, lchnk)
      call outfld( 'WICE' , wice, pcols, lchnk)
      call outfld( 'WEFF' , weff, pcols, lchnk)
      call outfld('INnso4  ',INnso4 , pcols,lchnk)
      call outfld('INnbc   ',INnbc  , pcols,lchnk)
      call outfld('INndust ',INndust, pcols,lchnk)
      call outfld('INondust ',INondust, pcols,lchnk)
      call outfld('INhet   ',INhet  , pcols,lchnk)
      call outfld('INhom   ',INhom  , pcols,lchnk)
      call outfld('INFrehom',INFrehom,pcols,lchnk)
      call outfld('INFreIN ',INFreIN, pcols,lchnk)
   end if

end subroutine nucleate_ice_cam_calc

!================================================================================================

end module nucleate_ice_cam
