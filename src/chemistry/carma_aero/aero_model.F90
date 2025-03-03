!===============================================================================
! CAMRA Aerosol Model
!===============================================================================
module aero_model
  use physics_buffer,    only: physics_buffer_desc, pbuf_get_index, pbuf_add_field, dtype_r8
  use shr_kind_mod,      only: r8 => shr_kind_r8
  use constituents,      only: pcnst, cnst_name, cnst_get_ind
  use perf_mod,          only: t_startf, t_stopf
  use ppgrid,            only: pcols, pver, pverp
  use phys_control,      only: phys_getopts, cam_physpkg_is
  use cam_abortutils,    only: endrun
  use cam_logfile,       only: iulog
  use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
  use camsrfexch,        only: cam_in_t, cam_out_t
  use physics_buffer,    only: pbuf_get_field, pbuf_set_field, dtype_r8
  use physconst,         only: gravit, rair, rhoh2o
  use spmd_utils,        only: masterproc
  use cam_history,       only: outfld
  use chem_mods,         only: gas_pcnst, adv_mass
  use mo_tracname,       only: solsym
  use infnan,            only: nan, assignment(=)
  use rad_constituents,  only: rad_cnst_get_info, rad_cnst_get_info_by_bin, &
                               rad_cnst_get_info_by_bin_spec, rad_cnst_get_bin_props_by_idx, &
                               rad_cnst_get_bin_mmr_by_idx
  use mo_setsox,         only: setsox, has_sox
  use carma_aerosol_properties_mod, only: carma_aerosol_properties

  use carma_intr, only: carma_get_group_by_name, carma_get_dry_radius, carma_get_wet_radius, carma_get_bin_rmass
  use carma_intr, only: carma_get_sad

  use aerosol_properties_mod, only: aero_name_len

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep     ! aerosol dry deposition and sediment
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: aero_model_surfarea    ! tropospheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea   ! stub

   ! Misc private data
  character(len=32), allocatable :: fieldname(:)    ! names for interstitial output fields
  character(len=32), allocatable :: fieldname_cw(:)    ! names for cloud_borne output fields

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: rprddp_idx          = 0
  integer :: rprdsh_idx          = 0
  integer :: nevapr_shcu_idx     = 0
  integer :: nevapr_dpcu_idx     = 0
  integer :: nh3_ndx    = 0
  integer :: nh4_ndx    = 0
  integer :: h2so4_ndx  = 0

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12


  ! description of bin aerosols
  integer, public, protected :: nspec_max = 0
  integer, public, protected :: nbins = 0
  integer, public, protected, allocatable :: nspec(:)

  ! local indexing for bins
  integer, allocatable :: bin_idx(:,:) ! table for local indexing of modal aero number and mmr
  integer :: ncnst_tot                  ! total number of mode number conc + mode species
  integer :: ncnst_extd                  ! twiece total number of mode number conc + mode species

  ! Indices for CARMA species in the ptend%q array.  Needed for prognostic aerosol case.
  logical, allocatable :: bin_cnst_lq(:,:)
  integer, allocatable :: bin_cnst_idx(:,:)


  ! ptr2d_t is used to create arrays of pointers to 2D fields
  type ptr2d_t
    real(r8), pointer :: fld(:,:) => null()
  end type ptr2d_t

  logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
                                 ! in the ptend object

  ! Namelist variables
  real(r8)          :: sol_facti_cloud_borne   = 1._r8
  real(r8)          :: sol_factb_interstitial  = 0.1_r8
  real(r8)          :: sol_factic_interstitial = 0.4_r8

  logical :: convproc_do_aer

  type(carma_aerosol_properties), pointer :: aero_props =>null()

contains

  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use aero_wetdep_cam, only: aero_wetdep_readnl
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! Namelist variables
    namelist /aerosol_nl/ sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aerosol_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(sol_facti_cloud_borne, 1,                         mpir8,   0, mpicom)
    call mpibcast(sol_factb_interstitial, 1,                        mpir8,   0, mpicom)
    call mpibcast(sol_factic_interstitial, 1,                       mpir8,   0, mpicom)
#endif

    call aero_wetdep_readnl(nlfile)

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register()

    use carma_flags_mod,   only: carma_model

    integer :: m, l, i
    integer :: nsoa_vbs
    character(len=32) :: num_name
    character(len=32) :: num_name_cw
    character(len=32) :: spec_name_cw

    integer :: idx, ierr

    call rad_cnst_get_info( 0, nbins=nbins)
    allocate( nspec(nbins), stat=ierr )
    if (ierr/=0) call endrun('aero_model_register: allocate error')

    ! add pbuf fields for interstitial (cloud borne) aerosols in CARMA
    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, num_name=num_name, num_name_cw=num_name_cw, nspec=nspec(m))
       call pbuf_add_field(num_name,'global',dtype_r8,(/pcols,pver/), idx)
       call pbuf_add_field(num_name_cw,'global',dtype_r8,(/pcols,pver/), idx)
       do l = 1, nspec(m)
          call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name_cw=spec_name_cw)
          call pbuf_add_field(spec_name_cw,'global',dtype_r8,(/pcols,pver/),idx)
       enddo
    enddo

  ! SOA information
  ! Define number of VBS bins (nsoa) based on number of SOAG chemistry species
    nsoa_vbs = 0
    do i = 1, pcnst
       if (cnst_name(i)(:4) == 'SOAG') then
          nsoa_vbs = nsoa_vbs + 1
       end if
    end do
    if (masterproc) then
       write(iulog,*) 'nsoa_vbs  = ', nsoa_vbs
    endif

   ! Define pbuf field for soa_fraction
    call pbuf_add_field('FRACVBS','global',dtype_r8,(/pcols,pver,nbins,nsoa_vbs/), idx)

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d )

    use mo_chem_utls,    only: get_inv_ndx
    use cam_history,     only: addfld, add_default, horiz_only
    use mo_chem_utls,    only: get_rxt_ndx, get_spc_ndx
    use aero_wetdep_cam, only: aero_wetdep_init
    use mo_setsox,       only: sox_inti
    use carma_aero_gasaerexch, only: carma_aero_gasaerexch_init

    use time_manager,    only: is_first_step
    use constituents,    only: cnst_set_convtran2
    use aero_deposition_cam, only: aero_deposition_cam_init

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    ! local vars
    character(len=*), parameter :: subrname = 'aero_model_init'
    integer :: m, n, ii, mm
    integer :: idxtmp    = -1

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_chemistry, history_cesm_forcing

    integer :: l

    character(len=2)  :: unit_basename  ! Units 'kg' or '1'
    character(len=32) :: num_name
    character(len=32) :: num_name_cw
    character(len=32) :: spec_name_cw

    integer :: idx, ierr
    real(r8) :: nanval

    aero_props => carma_aerosol_properties()
    call aero_deposition_cam_init(aero_props)

    if (is_first_step()) then
       do m = 1, nbins
          call rad_cnst_get_info_by_bin(0, m, num_name=num_name, num_name_cw=num_name_cw)
          idx = pbuf_get_index(num_name)
          call pbuf_set_field(pbuf2d, idx, 0.0_r8)
          idx = pbuf_get_index(num_name_cw)
          call pbuf_set_field(pbuf2d, idx, 0.0_r8)
          do l = 1, nspec(m)
             call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name_cw=spec_name_cw)
             idx = pbuf_get_index(spec_name_cw)
             call pbuf_set_field(pbuf2d, idx, 0.0_r8)
          enddo
       enddo
    endif

    ! define pbuf field for soa_fraction
    if (is_first_step()) then
       nanval = nan
       idx = pbuf_get_index('FRACVBS')
       call pbuf_set_field(pbuf2d, idx, nanval)
    end if

    ! aqueous chem initialization
    call sox_inti()

    h2so4_ndx = get_spc_ndx('H2SO4')
    nh3_ndx = get_spc_ndx('NH3')
    nh4_ndx = get_spc_ndx('NH4')

    fracis_idx      = pbuf_get_index('FRACIS')
    prain_idx       = pbuf_get_index('PRAIN')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    rprdsh_idx      = pbuf_get_index('RPRDSH')
    nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')

    call phys_getopts(history_aerosol_out = history_aerosol, &
                      history_chemistry_out=history_chemistry, &
                      history_cesm_forcing_out=history_cesm_forcing, &
                      convproc_do_aer_out = convproc_do_aer)

    call carma_aero_gasaerexch_init

    nspec_max = maxval(nspec)

    ncnst_tot = nspec(1)
    do m = 2, nbins
      ncnst_tot = ncnst_tot + nspec(m)
    end do
    ncnst_extd = 2*ncnst_tot

    allocate( &
      bin_idx(nbins,nspec_max),      &
      bin_cnst_lq(nbins,nspec_max), &
      bin_cnst_idx(nbins,nspec_max), &
      fieldname_cw(ncnst_tot), &
      fieldname(ncnst_tot), stat=ierr  )
    if (ierr/=0) call endrun(subrname//' : allocate error')

    ii = 0
    do m = 1, nbins
      do l = 1, nspec(m) ! loop through species
         ii = ii + 1
         bin_idx(m,l) = ii

         if (l <= nspec(m) ) then   ! species
            call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=fieldname(ii), spec_name_cw=fieldname_cw(ii))
         else  !number
            call rad_cnst_get_info_by_bin(0, m, num_name=fieldname(ii), num_name_cw=fieldname_cw(ii))
         end if

         call cnst_get_ind(fieldname(ii), idxtmp, abort=.false.)
          if (idxtmp.gt.0) then
             bin_cnst_lq(m,l) = .true.
             bin_cnst_idx(m,l) = idxtmp
             lq(idxtmp) = .true.
             call cnst_set_convtran2(idxtmp, .not.convproc_do_aer)
          else
             bin_cnst_lq(m,l) = .false.
             bin_cnst_idx(m,l) = 0
          end if

         mm = ii

         unit_basename = 'kg'
         if (l == nspec(m) + 2) then   ! number
          unit_basename = ' 1'
         end if


          call addfld( fieldname_cw(mm),                (/ 'lev' /), 'A', unit_basename//'/kg ',   &
               trim(fieldname_cw(mm))//' in cloud water')
          call addfld (trim(fieldname_cw(mm))//'DDF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' dry deposition flux at bottom (grav + turb)')
          call addfld (trim(fieldname_cw(mm))//'TBF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' turbulent dry deposition flux')
          call addfld (trim(fieldname_cw(mm))//'GVF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' gravitational dry deposition flux')

          if ( history_aerosol.or. history_chemistry ) then
             call add_default( fieldname_cw(mm), 1, ' ' )
          endif
          if ( history_aerosol ) then
             call add_default (trim(fieldname_cw(mm))//'GVF', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'TBF', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'DDF', 1, ' ')
          endif
       enddo
    enddo

    do m = 1,gas_pcnst

       unit_basename = 'kg'  ! Units 'kg' or '1'

       call addfld( 'GS_'//trim(solsym(m)),horiz_only, 'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' gas chemistry/wet removal (for gas species)')
       call addfld( 'AQ_'//trim(solsym(m)),horiz_only, 'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' aqueous chemistry (for gas species)')
       if ( history_aerosol ) then
          call add_default( 'AQ_'//trim(solsym(m)), 1, ' ')
       endif

    enddo

    if (has_sox) then
       do n = 1, nbins
          do l = 1, nspec(n)   ! not for total mass or number
             mm = bin_idx(n, l)
             call addfld (&
                  trim(fieldname_cw(mm))//'AQSO4',horiz_only,  'A','kg/m2/s', &
                  trim(fieldname_cw(mm))//' aqueous phase chemistry')
             call addfld (&
                  trim(fieldname_cw(mm))//'AQH2SO4',horiz_only,  'A','kg/m2/s', &
                  trim(fieldname_cw(mm))//' aqueous phase chemistry')
             if ( history_aerosol ) then
                call add_default (trim(fieldname_cw(mm))//'AQSO4', 1, ' ')
                call add_default (trim(fieldname_cw(mm))//'AQH2SO4', 1, ' ')
             endif
          end do
       end do

       call addfld( 'XPH_LWC',    (/ 'lev' /), 'A','kg/kg',   'pH value multiplied by lwc')
       call addfld ('AQSO4_H2O2', horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to H2O2')
       call addfld ('AQSO4_O3',   horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to O3')

       if ( history_aerosol ) then
          call add_default ('XPH_LWC', 1, ' ')
          call add_default ('AQSO4_H2O2', 1, ' ')
          call add_default ('AQSO4_O3', 1, ' ')
       endif
    endif

    call aero_wetdep_init()

  end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    ! args
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt        ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  endsubroutine aero_model_drydep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    use aero_wetdep_cam, only: aero_wetdep_tend

    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    call aero_wetdep_tend(state, dt, dlf, cam_out, ptend, pbuf)

  end subroutine aero_model_wetdep

  !-------------------------------------------------------------------------
  ! provides wet tropospheric aerosol surface area info for sectional aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  state, mmr, radmean, relhum, pmid, temp, strato_sad, sulfate,  m, ltrop, &
                  dlat, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_trop, reff_trop )

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: m(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_trop(:,:)  ! aerosol surface area density (cm2/cm3), zeroed above the tropopause
    real(r8), intent(out)   :: reff_trop(:,:) ! aerosol effective radius (cm), zeroed above the tropopause

    ! local vars
    integer :: beglev(ncol)
    integer :: endlev(ncol)

    beglev(:ncol)=ltrop(:ncol)+1
    endlev(:ncol)=pver
    call surf_area_dens( state, pbuf, ncol, mmr, beglev, endlev, sad_trop, reff_trop, sfc=sfc, dm_aer=dm_aer )

  end subroutine aero_model_surfarea

  !-------------------------------------------------------------------------
  ! provides wet stratospheric aerosol surface area info for sectional aerosols
  ! called from mo_gas_phase_chemdr.F90
  !-------------------------------------------------------------------------
  subroutine aero_model_strat_surfarea( state, ncol, mmr, pmid, temp, ltrop, pbuf, strato_sad, reff_strat )

    use ref_pres, only: clim_modal_aero_top_lev

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:) ! aerosol surface area density (cm2/cm3), zeroed below the tropopause
    real(r8), intent(out)   :: reff_strat(:,:) ! aerosol effective radius (cm), zeroed below the tropopause

    ! local vars
    integer :: beglev(ncol)
    integer :: endlev(ncol)

    beglev(:ncol) = clim_modal_aero_top_lev
    endlev(:ncol) = ltrop(:ncol)

    call surf_area_dens( state, pbuf, ncol, mmr, beglev, endlev, strato_sad, reff_strat )

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( state, loffset, ncol, lchnk, troplev, delt, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    use carma_aero_gasaerexch, only : carma_aero_gasaerexch_sub
    use time_manager,          only : get_nstep
    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    type(physics_state), intent(in)    :: state    ! Physics state variables
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    integer,  intent(in) :: troplev(:)
    real(r8), intent(in) :: delt                   ! time step size (sec)
    real(r8), intent(in) :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in) :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in) :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in) :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: relhum(:,:)            ! relative humidity
    real(r8), intent(in) :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:)
    real(r8), intent(in) :: zm(:,:)
    real(r8), intent(in) :: qh2o(:,:)
    real(r8), intent(in) :: cwat(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), intent(in) :: cldfr(:,:)
    real(r8), intent(in) :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8), intent(in) :: vmr0(:,:,:)       ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)         ! mixing ratios ( vmr )

    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars

    integer :: n, m, mm
    integer :: i,k,l
    integer :: nstep

    type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
    type(ptr2d_t), allocatable :: qqcw(:)

    real(r8) :: del_h2so4_aeruptk(ncol,pver)

    real(r8), pointer :: pblh(:)                    ! pbl height (m)

    real(r8), dimension(ncol) :: wrk
    character(len=32)         :: name
    real(r8) :: dvmrcwdt(ncol,pver,ncnst_tot)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: delta_so4mass(ncol,pver,ncnst_tot)
    real(r8) :: wetr_n(pcols,pver,nbins)       ! wet radius from CARMA for different bin
    real(r8) :: vmrcw(ncol,pver,ncnst_tot)            ! cloud-borne aerosol (vmr)
    real(r8) :: mmrcw(ncol,pver,ncnst_tot)            ! cloud-borne aerosol (mmr)
    real(r8) :: raervmr(ncol,pver,ncnst_tot)            ! cloud-borne aerosol (vmr)

    real(r8) ::  aqso4(ncol,ncnst_tot)               ! aqueous phase chemistry
    real(r8) ::  aqh2so4(ncol,ncnst_tot)             ! aqueous phase chemistry
    real(r8) ::  aqso4_h2o2(ncol)                     ! SO4 aqueous phase chemistry due to H2O2
    real(r8) ::  aqso4_o3(ncol)                       ! SO4 aqueous phase chemistry due to O3
    real(r8) ::  xphlwc(ncol,pver)                    ! pH value multiplied by lwc
    real(r8) ::  nh3_beg(ncol,pver)
    real(r8) ::  mw_carma(ncnst_tot)
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: sulfeq(:,:,:)
    real(r8) :: wetr(pcols,pver)   ! CARMA wet radius in cm
    real(r8) :: wetrho(pcols,pver)   ! CARMA wet dens
    real(r8), allocatable :: rmass(:)     ! CARMA rmass

    real(r8) :: old_total_mass
    real(r8) :: new_total_mass
    real(r8) :: old_total_number

    character(len=32) :: spectype

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr, ierr
    character(len=*), parameter :: subname = 'aero_model_gasaerexch'

!
! ... initialize nh3
!
    if ( nh3_ndx > 0 ) then
      nh3_beg = vmr(1:ncol,:,nh3_ndx)
    end if
!
! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    nstep = get_nstep()

    ! calculate tendency due to gas phase chemistry and processes
    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0.0_r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m)*adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'GS_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

!
! Aerosol processes ...
!
    allocate( &
      rmass(nbins), &
      raer(ncnst_tot), &
      qqcw(ncnst_tot), stat=ierr )
    if (ierr /= 0) call endrun(subname//': allocate error')

    mw_carma(:) = 0.0_r8
    do m = 1, nbins      ! main loop over aerosol bins
       ! dryr is the dry bin radius
       ! wetr is the dry bin radius
       ! Note: taken here from CARMA pbuf field which may be not any more consistent with changed fields after carma was applied
       ! Need to add new code that recalcuates dryr and wetr
       ! get bin info
       call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m), bin_name=bin_name)

       nchr = len_trim(bin_name)-2
       shortname = bin_name(:nchr)

       call carma_get_group_by_name(shortname, igroup, rc)
       if (rc/=0) then
          call endrun(subname//': ERROR in carma_get_group_by_name')
       end if

       read(bin_name(nchr+1:),*) ibin

       call carma_get_wet_radius(state, igroup, ibin, wetr, wetrho, rc) ! m
       if (rc/=0) then
          call endrun(subname//': ERROR in carma_get_wet_radius')
       end if
       wetr(:ncol,:) = wetr(:ncol,:) * 1.e2_r8 ! cm

       call carma_get_bin_rmass(igroup, ibin, rmass(m), rc) ! grams
       if (rc/=0) then
          call endrun(subname//': ERROR in carma_get_bin_rmass')
       end if

       wetr_n(:,:,m) = wetr(:,:)

       ! Init pointers to mode number and specie mass mixing ratios in
       ! intersitial and cloud borne phases.
       do l = 1, nspec(m)
          mm = bin_idx(m, l)
          if (l <= nspec(m)) then
             call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
             call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'a', state, pbuf, raer(mm)%fld)
             call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
             if (trim(spectype) == 'sulfate') then
                mw_carma(mm) = 96._r8
             end if
             if (trim(spectype) == 'black-c') then
                mw_carma(mm) = 12._r8
             end if
             if (trim(spectype) == 'p-organic') then
                mw_carma(mm) = 12._r8
             end if
             if (trim(spectype) == 's-organic') then
                mw_carma(mm) = 250._r8
             end if
             if (trim(spectype) == 'dust') then
                mw_carma(mm) = 12._r8
             end if
             if (trim(spectype) == 'seasalt') then
                mw_carma(mm) = 57._r8
             end if
          end if
          mmrcw(:ncol,:,mm) = qqcw(mm)%fld(:ncol,:)
          vmrcw(:ncol,:,mm) = qqcw(mm)%fld(:ncol,:)
          raervmr(:ncol,:,mm) = raer(mm)%fld(:ncol,:)
       end do
    end do

    ! qqcw2vrm is different from what is done in MAM, here we pass in the fields set by the qqcw and raer pointer
    ! for all the CARMA aerosols, species, mmr, and number, vmrcw (kg/kg) -> vmr
    call mmr2vmr_carma ( lchnk, vmrcw, mbar, mw_carma, ncol, loffset, rmass )

    dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)   ! all adveced species no aerosols
    dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)  ! cloud borne carma aerosol species
    ! aqueous chemistry ...

    if( has_sox ) then
         call setsox( state,  &
              ncol,     &
              lchnk,    &
              loffset,  &
              delt,     &
              pmid,     &
              pdel,     &
              tfld,     &
              mbar,     &
              cwat,     &
              cldfr,    &
              cldnum,   &
              airdens,  &
              invariants, &
              vmrcw,    &
              vmr,      &
              xphlwc,   &
              aqso4,    &
              aqh2so4,  &
              aqso4_h2o2, &
              aqso4_o3  &
              )

          do n = 1, nbins
            do l = 1, nspec(n)   ! not for total mass or number
                mm = bin_idx(n, l)
                call outfld( trim(fieldname_cw(mm))//'AQSO4',   aqso4(:ncol,mm),   ncol, lchnk)
                call outfld( trim(fieldname_cw(mm))//'AQH2SO4', aqh2so4(:ncol,mm), ncol, lchnk)
             end do
          end do

          call outfld( 'AQSO4_H2O2', aqso4_h2o2(:ncol), ncol, lchnk)
          call outfld( 'AQSO4_O3',   aqso4_o3(:ncol),   ncol, lchnk)
          call outfld( 'XPH_LWC',    xphlwc(:ncol,:),   ncol, lchnk )

    endif

!   Tendency due to aqueous chemistry
    dvmrdt = (vmr - dvmrdt) / delt
    dvmrcwdt = (vmrcw - dvmrcwdt) / delt

    do m = 1, gas_pcnst
       wrk(:) = 0.0_r8
       do k = 1,pver
          wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m) * adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
       end do
       name = 'AQ_'//trim(solsym(m))
       call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    if (h2so4_ndx > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,h2so4_ndx)
    else
       del_h2so4_aeruptk(:,:) = 0.0_r8
    endif

    ! need to transform raer to raervmr from CARMA, routine requires vmr, note number wil not be changed here
    call mmr2vmr_carma ( lchnk, raervmr, mbar, mw_carma, ncol, loffset, rmass)

    call carma_aero_gasaerexch_sub( state, &
          pbuf, lchnk,    ncol,     nstep,      &
          loffset,            delt, mbar ,      &
          tfld,     pmid,     pdel,             &
          qh2o,               troplev,          &
          vmr,                raervmr,          &
          wetr_n     )

    ! note vmr2qqcw does not change qqcw pointer (different than in MAM)
    call vmr2mmr_carma ( lchnk, vmrcw, mbar, mw_carma, ncol, loffset, rmass )

    !vmrcw in kg/kg
    ! change pointer value for total mmr and number. In order to do this correctly
    ! only mass has to be added to each bin (not number). This will require redistributing
    ! mass to different bins. Here, we change both mass and number until we have a better
    ! solution.
    delta_so4mass(:,:,:) = 0.0_r8
    do m = 1, nbins
       do l = 1, nspec(m)  ! for sulfate only
          mm = bin_idx(m, l)
         ! sulfate mass that needs to be added to the total mass
          call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
          if (trim(spectype) == 'sulfate') then
              ! only do loop if vmrcw has changed
              do k=1,pver
                 do i=1,ncol
                  if (vmrcw(i,k,mm) .gt. mmrcw(i,k,mm) .and. mmrcw(i,k,mm) /= 0.0_r8)  then
                   delta_so4mass(i,k,mm) = ( vmrcw(i,k,mm) - mmrcw(i,k,mm) )
                  else
                    delta_so4mass(i,k,mm) = 0.0_r8
                  end if
                 end do
              end do
         end if
       end do
    end do

    do m = 1, nbins
       do l = 1, nspec(m) ! for sulfate only
          mm = bin_idx(m, l)
          qqcw(mm)%fld(:ncol,:) = vmrcw(:ncol,:,mm)
          call outfld( trim(fieldname_cw(mm)), qqcw(mm)%fld(:ncol,:), ncol, lchnk)
       end do
    end do


  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

  end subroutine aero_model_emissions


  !===============================================================================
  !===============================================================================
  ! private methods


  !=============================================================================
  !=============================================================================
  subroutine surf_area_dens( state, pbuf, ncol, mmr, beglev, endlev, sad, reff, sfc, dm_aer )
    use mo_constants, only: pi
    use carma_intr,   only: carma_effecitive_radius

    ! dummy args
    type(physics_state),    intent(in) :: state           ! Physics state variables
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: mmr(:,:,:)
    integer,  intent(in)  :: beglev(:)
    integer,  intent(in)  :: endlev(:)
    real(r8), intent(out) :: sad(:,:)    ! bulk surface area density in cm2/cm3 from beglev to endlev, zero elsewhere
    real(r8), intent(out) :: reff(:,:)   ! bulk effective radius in cm from beglev to endlev, zero elsewhere
    real(r8), optional, intent(out) :: sfc(:,:,:) ! surface area density per bin
    real(r8), optional, intent(out) :: dm_aer(:,:,:) ! diameter per bin

    ! local vars
    real(r8) :: reffaer(pcols,pver) ! bulk effective radius in cm

    real(r8) :: sad_bin(pcols,pver,nbins)
    integer  :: icol, ilev, ibin, ispec !!, reff_pbf_ndx
    real(r8) :: chm_mass, tot_mass
    character(len=32) :: spectype
    real(r8) :: wetr(pcols,pver)      ! CARMA bin wet radius in cm
    real(r8) :: wetrho(pcols,pver)    ! CARMA bin wet density
    real(r8) :: sad_carma(pcols,pver) ! CARMA bin wet surface area density in cm2/cm3
    real(r8), pointer :: aer_bin_mmr(:,:)

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, indxbin, rc, nchr

    sad = 0._r8
    reff = 0._r8

    !
    ! Compute surface aero for each bin.
    ! Total over all bins as the surface area for chemical reactions.
    !

    reffaer = carma_effecitive_radius(state)

    sad = 0._r8
    sad_bin = 0._r8
    reff = 0._r8

    do ibin=1,nbins ! loop over aerosol bins
      call rad_cnst_get_info_by_bin(0, ibin, bin_name=bin_name)

      nchr = len_trim(bin_name)-2
      shortname = bin_name(:nchr)

      call carma_get_group_by_name(shortname, igroup, rc)

      read(bin_name(nchr+1:),*) indxbin

      call carma_get_wet_radius(state, igroup, indxbin, wetr, wetrho, rc) ! m
      wetr(:ncol,:) = wetr(:ncol,:) * 1.e2_r8 ! cm
      call carma_get_sad(state, igroup, indxbin, sad_carma, rc)

      if (present(dm_aer)) then
         dm_aer(:ncol,:,ibin) = 2._r8 * wetr(:ncol,:) ! convert wet radius (cm) to wet diameter (cm)
      endif
      sad_bin(:ncol,:,ibin) = sad_carma(:ncol,:) ! cm^2/cm^3
    end do

    do icol = 1,ncol
      do ilev = beglev(icol),endlev(icol)
        do ibin=1,nbins ! loop over aerosol bins
          !
          ! compute a mass weighting of the number
          !
          tot_mass = 0._r8
          chm_mass = 0._r8
          do ispec=1,nspec(ibin)

             call rad_cnst_get_bin_mmr_by_idx(0, ibin, ispec, 'a', state, pbuf, aer_bin_mmr)

             tot_mass = tot_mass + aer_bin_mmr(icol,ilev)

             call rad_cnst_get_bin_props_by_idx(0, ibin, ispec, spectype=spectype)

             if ( trim(spectype) == 'sulfate'   .or. &
                trim(spectype) == 's-organic' .or. &
                trim(spectype) == 'p-organic' .or. &
                trim(spectype) == 'black-c'   .or. &
                trim(spectype) == 'ammonium') then
                chm_mass = chm_mass + aer_bin_mmr(icol,ilev)
             end if

          end do
          if ( tot_mass > 0._r8 ) then
         ! surface area density
            sad_bin(icol,ilev,ibin) = chm_mass / tot_mass * sad_bin(icol,ilev,ibin) ! cm^2/cm^3
          else
            sad_bin(icol,ilev,ibin) = 0._r8
          end if
        end do
        sad(icol,ilev) = sum(sad_bin(icol,ilev,:))
        reff(icol,ilev) = reffaer(icol,ilev)

       end do
    end do

    if (present(sfc)) then
       sfc(:,:,:) = sad_bin(:,:,:)
    endif

  end subroutine surf_area_dens

  !=============================================================================
  subroutine mmr2vmr_carma(lchnk, vmr, mbar, mw_carma, ncol, im, rmass)
    !-----------------------------------------------------------------
    !   ... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !   ... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: rmass(nbins)
    real(r8), intent(in)    :: mw_carma(ncnst_tot)
    real(r8), intent(inout) :: vmr(ncol,pver,ncnst_tot)

    !-----------------------------------------------------------------
    !   ... Local variables
    !-----------------------------------------------------------------
    integer :: k, m, mm, l

    do m = 1, nbins
       do l = 1, nspec(m)   ! for each species, not total mmr or number, information of mw are missing
          mm = bin_idx(m, l)
          do k=1,pver
             vmr(:ncol,k,mm) = mbar(:ncol,k) * vmr(:ncol,k,mm) / mw_carma(mm)
          end do
       end do
    end do

  end subroutine mmr2vmr_carma
  !=============================================================================

  !=============================================================================
  subroutine vmr2mmr_carma ( lchnk, vmr, mbar, mw_carma, ncol, im, rmass )
    !-----------------------------------------------------------------
    !   ... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !   ... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: rmass(nbins)
    real(r8), intent(inout)    :: vmr(ncol,pver,ncnst_tot)
    real(r8), intent(in)    :: mw_carma(ncnst_tot)

    !-----------------------------------------------------------------
    !   ... Local variables
    !-----------------------------------------------------------------
    integer :: k, m, mm, l
    !-----------------------------------------------------------------
    !   ... The non-group species
    !-----------------------------------------------------------------
    do m = 1, nbins
       do l = 1, nspec(m)   ! for each species, not total mmr or number, information of mw are missing
          mm = bin_idx(m, l)
          do k=1,pver
             vmr(:ncol,k,mm) = mw_carma(mm) * vmr(:ncol,k,mm) / mbar(:ncol,k)
          end do
       end do
    end do

  end subroutine vmr2mmr_carma

end module aero_model
