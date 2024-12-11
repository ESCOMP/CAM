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
  use carma_intr, only: carma_get_total_mmr, carma_get_sad

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

  ! number of modes
  integer :: pblh_idx            = 0
  integer :: wetdens_ap_idx      = 0

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: rprddp_idx          = 0
  integer :: rprdsh_idx          = 0
  integer :: nevapr_shcu_idx     = 0
  integer :: nevapr_dpcu_idx     = 0

  integer :: sulfeq_idx = -1

  integer :: nh3_ndx    = 0
  integer :: nh4_ndx    = 0
  integer :: h2so4_ndx  = 0

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
  real(r8) :: dlndg_nimptblgrow
  real(r8),allocatable :: scavimptblnum(:,:)
  real(r8),allocatable :: scavimptblvol(:,:)


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
  real(r8)          :: seasalt_emis_scale

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
    !st character(len=16) :: aer_wetdep_list(pcnst) = ' '
    !st character(len=16) :: aer_drydep_list(pcnst) = ' '

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
    !st call mpibcast(aer_wetdep_list,   len(aer_wetdep_list(1))*pcnst, mpichar, 0, mpicom)
    !st call mpibcast(aer_drydep_list,   len(aer_drydep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(sol_facti_cloud_borne, 1,                         mpir8,   0, mpicom)
    call mpibcast(sol_factb_interstitial, 1,                        mpir8,   0, mpicom)
    call mpibcast(sol_factic_interstitial, 1,                       mpir8,   0, mpicom)
    !st call mpibcast(modal_strat_sulfate,     1,                       mpilog,  0, mpicom)
    !st call mpibcast(seasalt_emis_scale, 1,                            mpir8,   0, mpicom)
    !st call mpibcast(modal_accum_coarse_exch, 1,                       mpilog,  0, mpicom)
#endif

    call aero_wetdep_readnl(nlfile)

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register()

    use carma_flags_mod,   only: carma_model

    integer :: m, l, i
    integer :: nsoa_vbs
    character(len=32) :: spectype
    character(len=32) :: num_name
    character(len=32) :: num_name_cw
    character(len=32) :: spec_name_cw
    character(len=32) :: soag_name
    character(len=32) :: soa_name

    integer :: idx

    call rad_cnst_get_info( 0, nbins=nbins)
    allocate( nspec(nbins) )

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
    !st use modal_aero_data, only: cnst_name_cw
    !st use modal_aero_data, only: modal_aero_data_init
    !st use dust_model,      only: dust_init, dust_names, dust_active, dust_nbin, dust_nnum
    !st use seasalt_model,   only: seasalt_init, seasalt_names, seasalt_active,seasalt_nbin
    !st use drydep_mod,      only: inidrydep
    use aero_wetdep_cam, only: aero_wetdep_init
    use mo_setsox,       only: sox_inti

    !st use modal_aero_calcsize,   only: modal_aero_calcsize_init
    !st use modal_aero_coag,       only: modal_aero_coag_init
    !st use modal_aero_deposition, only: modal_aero_deposition_init
    use carma_aero_gasaerexch, only: carma_aero_gasaerexch_init
    !st use modal_aero_newnuc,     only: modal_aero_newnuc_init
    !st use modal_aero_rename,     only: modal_aero_rename_init

    use time_manager,    only: is_first_step
    use constituents,    only: cnst_set_convtran2
    use aero_deposition_cam, only: aero_deposition_cam_init

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    ! local vars
    character(len=*), parameter :: subrname = 'aero_model_init'
    integer :: m, n, id, ii, mm
    integer :: lptr      = -1
    integer :: idxtmp    = -1
    character(len=20) :: dummy

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_chemistry, history_cesm_forcing, history_dust

    integer :: l
    character(len=6) :: test_name
    character(len=64) :: errmes

    character(len=2)  :: unit_basename  ! Units 'kg' or '1'
    integer :: errcode
    !st character(len=fieldname_len) :: field_name

    character(len=32) :: spectype
    character(len=32) :: num_name
    character(len=32) :: num_name_cw
    character(len=32) :: spec_name_cw

    integer :: idx
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
    !st sulfeq_idx      = pbuf_get_index('MAMH2SO4EQ',errcode)

    !st not sure if this is needed
    call phys_getopts(history_aerosol_out = history_aerosol, &
                      history_chemistry_out=history_chemistry, &
                      history_cesm_forcing_out=history_cesm_forcing, &
                      convproc_do_aer_out = convproc_do_aer)

!!$    call carma_aero_bcscavcoef_init(pbuf2d)

    !st  call modal_aero_rename_init( modal_accum_coarse_exch )
    !   calcsize call must follow rename call
    !st call modal_aero_calcsize_init( pbuf2d )
    call carma_aero_gasaerexch_init
    !   coag call must follow gasaerexch call
    !st call modal_aero_coag_init
    !st call modal_aero_newnuc_init

    ! call modal_aero_deposition_init only if the user has not specified
    ! prescribed aerosol deposition fluxes
    !st if (.not.aerodep_flx_prescribed()) then
    !st   call modal_aero_deposition_init
    !stendif


    !st all CARMA species are deposited, therefore the following is not used
    !st nwetdep = 0
    !st ndrydep = 0

    !st count_species: do m = 1,pcnst
    !st    if ( len_trim(wetdep_list(m)) /= 0 ) then
    !st       nwetdep = nwetdep+1
    !st    endif
    !st    if ( len_trim(drydep_list(m)) /= 0 ) then
    !st       ndrydep = ndrydep+1
    !st    endif
    !st enddo count_species

    ! add plus one to include number, total mmr and nspec
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
      fieldname(ncnst_tot) )

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

    !st real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:), wetdens(:,:,:)
    real(r8), pointer :: pblh(:)                    ! pbl height (m)

    real(r8), dimension(ncol) :: wrk
    character(len=32)         :: name
    real(r8) :: dvmrcwdt(ncol,pver,ncnst_tot)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: delta_so4mass(ncol,pver,ncnst_tot)
    real(r8) :: wetr_n(pcols,pver,nbins)       ! wet radius from CARMA for different bin
    !st real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)
    !st vmrcw is going only through CARMA aerosols (ncnst_tot)
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
    logical :: is_spcam_m2005

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr
    character(len=*), parameter :: subname = 'aero_model_gasaerexch'

!
! ... initialize nh3
!
    if ( nh3_ndx > 0 ) then
      nh3_beg = vmr(1:ncol,:,nh3_ndx)
    end if
!
    is_spcam_m2005   = cam_physpkg_is('spcam_m2005')

    !st call pbuf_get_field(pbuf, dgnum_idx,      dgnum)
    !st call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet )
    !st call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens )
    !st call pbuf_get_field(pbuf, pblh_idx,       pblh)

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
      rmass(nbins),                &
      raer(ncnst_tot),                &
      qqcw(ncnst_tot)                 )

    mw_carma(:) = 0.0_r8
    do m = 1, nbins      ! main loop over aerosol bins
       !st can we move this part to init???
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
          !write(iulog,*) 'qqcw(mm)%fld) ', qqcw(mm)%fld(:ncol,:)
       end do
    end do
    !write(iulog,*) 'vmrcw(:,:,1) start', maxval(vmrcw(:ncol,:,1) )

    !write(iulog,*) 'mm start vmrcw, raervmr'
    ! qqcw2vrm is different from what is done in MAM, here we pass in the fields set by the qqcw and raer pointer
    ! for all the CARMA aerosols, species, mmr, and number, vmrcw (kg/kg) -> vmr
    call mmr2vmr_carma ( lchnk, vmrcw, mbar, mw_carma, ncol, loffset, rmass )
    !write(iulog,*) 'vmrcw(:,:,1) mmr', maxval(vmrcw(:,:,1))

    if (.not. is_spcam_m2005) then  ! regular CAM
       dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)   ! all adveced species no aerosols
       dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)  ! cloud borne carma aerosol species
    ! aqueous chemistry ...
    ! write(iulog,*) 'start has_sox'

    if( has_sox ) then
         call setsox( state,  &
              pbuf,     &
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

      !write(iulog,*) 'done with has_sox'
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

    else if (is_spcam_m2005) then  ! SPCAM ECPP
! when ECPP is used, aqueous chemistry is done in ECPP,
! and not updated here.
! Minghuai Wang, 2010-02 (Minghuai.Wang@pnl.gov)

      dvmrdt = 0.0_r8
      dvmrcwdt = 0.0_r8
    endif

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    if (h2so4_ndx > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,h2so4_ndx)
    else
       del_h2so4_aeruptk(:,:) = 0.0_r8
    endif


    !call t_startf('modal_gas-aer_exchng')

    !if ( sulfeq_idx>0 ) then
    !   call pbuf_get_field( pbuf, sulfeq_idx, sulfeq )
    !else
    !   nullify( sulfeq )
    !endif
    !write(iulog,*) 'start carma_aero_gasaerexch_sub'
    ! need to transform raer to raervmr from CARMA, routine requires vmr, note number wil not be changed here
    call mmr2vmr_carma ( lchnk, raervmr, mbar, mw_carma, ncol, loffset, rmass)
    !write(iulog,*) 'mm start raervmr done'

    call carma_aero_gasaerexch_sub( state, &
          pbuf, lchnk,    ncol,     nstep,      &
          loffset,            delt, mbar ,      &
          tfld,     pmid,     pdel,             &
          qh2o,               troplev,          &
          vmr,                raervmr,          &
          wetr_n     )

    !if (h2so4_ndx > 0) then
    !   del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,h2so4_ndx) - del_h2so4_aeruptk(1:ncol,:)
    !endif

    !call t_stopf('modal_gas-aer_exchng')


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

    ! Is the loop here needed?
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

    real(r8), pointer, dimension(:,:) :: cmass,tmass ! carma element chemical and total mass
    real(r8) :: sad_bin(pcols,pver,nbins)
    integer  :: err, icol, ilev, ibin, ispec !!, reff_pbf_ndx
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

!!$  !===============================================================================
!!$  !===============================================================================
!!$  subroutine carma_aero_bcscavcoef_init ( pbuf2d )
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    ! Purpose:
!!$    ! Computes lookup table for aerosol impaction/interception scavenging rates
!!$    !
!!$    ! Authors: R. Easter
!!$    ! Simone Tilmes Nov 2021
!!$    ! added modifications for bin model, assuming sigma = 1.
!!$    !
!!$    !-----------------------------------------------------------------------
!!$
!!$    use shr_kind_mod,    only: r8 => shr_kind_r8
!!$    use cam_abortutils,  only: endrun
!!$    use mo_constants, only:  pi
!!$    use ppgrid,          only: begchunk
!!$
!!$    implicit none
!!$
!!$    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
!!$
!!$    !   local variables
!!$    integer nnfit_maxd
!!$    parameter (nnfit_maxd=27)
!!$
!!$    integer m, i, l, jgrow, jdens, jpress, jtemp, nnfit
!!$    integer lunerr
!!$
!!$    character(len=32) :: bin_name
!!$    character(len=32) :: spectype
!!$
!!$    real(r8) dg0, dg0_cgs, press, dg0_base, &
!!$         rhodryaero, rhowetaero, rhowetaero_cgs, rmserr, &
!!$         scavratenum, scavratevol, sigmag,                &
!!$         temp, wetdiaratio, wetvolratio
!!$    real(r8) :: specdens
!!$    real(r8) aafitnum(1), xxfitnum(1,nnfit_maxd), yyfitnum(nnfit_maxd)
!!$    real(r8) aafitvol(1), xxfitvol(1,nnfit_maxd), yyfitvol(nnfit_maxd)
!!$
!!$
!!$    allocate(scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, nbins))
!!$    allocate(scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, nbins))
!!$
!!$    lunerr = iulog
!!$    dlndg_nimptblgrow = log( 1.25_r8 )
!!$
!!$    ! bin model: main loop over aerosol bins
!!$
!!$    modeloop: do m = 1, nbins
!!$       !write(*,*) 'mloop start ',m
!!$       ! r(m) is the dry bin radius
!!$       ! taken here from CARMA pbuf field
!!$       ! get bin info
!!$       call rad_cnst_get_info_by_bin(0, m, bin_name=bin_name)
!!$
!!$       !   for setting up the lookup table, use the dry density of the first
!!$       ! get specdens from sulfate (check)
!!$       do l = 1, nspec(m)
!!$          call aero_props%species_type(m,l, spectype)
!!$          if (trim(spectype) == 'sulfate') then
!!$             call aero_props%get(m,l,density=rhodryaero)
!!$          end if
!!$       end do
!!$
!!$       dg0_base = 2._r8 * aero_props%scav_radius(m)
!!$
!!$           !sigmag = sigmag_amode(mode)
!!$           !dg0_base = dcen_sect(m,n)*exp( -1.5*((log(sigmag))**2) )
!!$           ! for bin approach sigma assumed to be 1., dg0_base equal dry radius
!!$       sigmag = 1._r8
!!$
!!$
!!$       !st rhodryaero = specdens_amode(1,mode)
!!$
!!$       growloop: do jgrow = nimptblgrow_mind, nimptblgrow_maxd
!!$
!!$          wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
!!$          !dg0 = dgnum_amode(mode)*wetdiaratio
!!$          dg0 = dg0_base*wetdiaratio
!!$          !st write(*,*) 'm,l,dg0 ',m,l,dg0
!!$
!!$          wetvolratio = exp( jgrow*dlndg_nimptblgrow*3._r8 )
!!$          rhowetaero = 1.0_r8 + (rhodryaero-1.0_r8)/wetvolratio
!!$          rhowetaero = min( rhowetaero, rhodryaero )
!!$
!!$          !
!!$          !   compute impaction scavenging rates at 1 temp-press pair and save
!!$          !
!!$          nnfit = 0
!!$
!!$          temp = 273.16_r8
!!$          press = 0.75e6_r8   ! dynes/cm2
!!$          rhowetaero = rhodryaero
!!$
!!$          ! CARMA dry radius is in  cm
!!$          !dg0_cgs = dg0*1.0e2_r8   ! m to cm
!!$          dg0_cgs = dg0   ! CARMA  radius / diameter is already in cm
!!$
!!$          rhowetaero_cgs = rhowetaero*1.0e-3_r8   ! kg/m3 to g/cm3
!!$
!!$
!!$          call calc_1_impact_rate( &
!!$               dg0_cgs, sigmag, rhowetaero_cgs, temp, press, &
!!$               scavratenum, scavratevol, lunerr )
!!$
!!$
!!$          nnfit = nnfit + 1
!!$          if (nnfit > nnfit_maxd) then
!!$             write(lunerr,9110)
!!$             call endrun()
!!$          end if
!!$9110      format( '*** subr. carma_aero_bcscavcoef_init -- nnfit too big' )
!!$
!!$          xxfitnum(1,nnfit) = 1._r8
!!$          yyfitnum(nnfit) = log( scavratenum )
!!$
!!$          xxfitvol(1,nnfit) = 1._r8
!!$          yyfitvol(nnfit) = log( scavratevol )
!!$
!!$          !
!!$          ! skip mlinfit stuff because scav table no longer has dependencies on
!!$          !    air temp, air press, and particle wet density
!!$          ! just load the log( scavrate--- ) values
!!$          !
!!$          !!
!!$          !!   do linear regression
!!$          !!    log(scavrate) = a1 + a2*log(wetdens)
!!$          !!
!!$          !     call mlinft( xxfitnum, yyfitnum, aafitnum, nnfit, 1, 1, rmserr )
!!$          !     call mlinft( xxfitvol, yyfitvol, aafitvol, nnfit, 1, 1, rmserr )
!!$          !
!!$          !     scavimptblnum(jgrow,mode) = aafitnum(1)
!!$          !     scavimptblvol(jgrow,mode) = aafitvol(1)
!!$
!!$         !depends on both bins and different species
!!$          scavimptblnum(jgrow,m) = yyfitnum(1)
!!$          scavimptblvol(jgrow,m) = yyfitvol(1)
!!$
!!$       enddo growloop
!!$    enddo modeloop
!!$
!!$    return
!!$  end subroutine carma_aero_bcscavcoef_init
!!$
!!$  !===============================================================================
!!$  !===============================================================================
!!$
!!$
!!$  !===============================================================================
!!$  subroutine carma_aero_bcscavcoef_get( m, ncol, isprx, wetr, dryr, scavcoefnum, scavcoefvol, pbuf )
!!$    !  need to go through both bins and species
!!$    ! need dry radius and wet radius
!!$
!!$    !-----------------------------------------------------------------------
!!$
!!$    use mo_constants, only:  pi
!!$
!!$    implicit none
!!$
!!$    integer,intent(in) :: m, ncol
!!$    logical,intent(in):: isprx(pcols,pver)
!!$    ! wet radius per bin dgn_awet -> wetr
!!$    real(r8), intent(in) :: dryr(pcols,pver)
!!$    real(r8), intent(in) :: wetr(pcols,pver)
!!$    real(r8), intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)
!!$    type(physics_buffer_desc), pointer :: pbuf(:)
!!$
!!$    integer i, k, jgrow, l
!!$    real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum, dg0_base, specdens, rhodryaero
!!$
!!$    character(len=32) :: spectype
!!$    character(len=aero_name_len) :: bin_name, shortname
!!$    integer :: igroup, ibin, rc, nchr
!!$
!!$    real(r8), allocatable :: rmass(:)     ! CARMA rmass
!!$    character(len=*), parameter :: subname = 'carma_aero_bcscavcoef_get'
!!$
!!$    allocate ( rmass(nbins) )
!!$    ! bin model: main loop over aerosol bins
!!$
!!$    ! get bin info
!!$     call rad_cnst_get_info_by_bin(0, m, bin_name=bin_name)
!!$
!!$     nchr = len_trim(bin_name)-2
!!$     shortname = bin_name(:nchr)
!!$
!!$     call carma_get_group_by_name(shortname, igroup, rc)
!!$
!!$     read(bin_name(nchr+1:),*) ibin
!!$
!!$     call carma_get_bin_rmass(igroup, ibin, rmass(m), rc)
!!$     if (rc/=0) then
!!$        call endrun(subname//': ERROR in carma_get_bin_rmass')
!!$     end if
!!$
!!$    ! get rmass and specdens for sulfate
!!$    do l = 1, nspec(m)
!!$        call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype, density_aer=specdens)
!!$
!!$       !   chemical component of the aerosol type (which currently will be so4)
!!$       !   For  CARMA, rmass per bin stays the same, while dry radius varies when the particle density varies
!!$       !   rmass = 4/3 * Pi * density * dry radius
!!$       !   We assume a fixed specie density
!!$        if (trim(spectype) == 'sulfate') then
!!$           rhodryaero = specdens
!!$        end if
!!$    end do
!!$    dg0_base = 2._r8 * (0.75_r8*rmass(m) / pi  / (1.0e-3_r8*rhodryaero)) **(0.33_r8)    ! specdens kg/m3 to g/cm3, convert from radiust to diameter
!!$    !rg0_base = (0.75_r8*rmass(m) / pi  / (1.0e-3_r8*specdens)) **(0.33_r8)    ! specdens kg/m3 to g/cm3
!!$
!!$    do k = 1, pver
!!$       do i = 1, ncol
!!$
!!$          ! do only if no precip
!!$          if ( isprx(i,k) .and. dryr(i,k).gt.0._r8) then
!!$             !
!!$             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)
!!$
!!$             ! dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)
!!$             ! dgnum_amode(m) is the rg0_base radius.
!!$
!!$             dumdgratio = wetr(i,k)/dg0_base
!!$
!!$              if ((dumdgratio >= 0.99_r8) .and. (dumdgratio <= 1.01_r8)) then
!!$                 scavimpvol = scavimptblvol(0,m)
!!$                 scavimpnum = scavimptblnum(0,m)
!!$              else
!!$                xgrow = log( dumdgratio ) / dlndg_nimptblgrow
!!$                jgrow = int( xgrow )
!!$                if (xgrow < 0._r8) jgrow = jgrow - 1
!!$                if (jgrow < nimptblgrow_mind) then
!!$                   jgrow = nimptblgrow_mind
!!$                   xgrow = jgrow
!!$                else
!!$                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
!!$                end if
!!$
!!$                dumfhi = xgrow - jgrow
!!$                dumflo = 1._r8 - dumfhi
!!$
!!$                scavimpvol = dumflo*scavimptblvol(jgrow,m) + &
!!$                     dumfhi*scavimptblvol(jgrow+1,m)
!!$                scavimpnum = dumflo*scavimptblnum(jgrow,m) + &
!!$                     dumfhi*scavimptblnum(jgrow+1,m)
!!$
!!$             end if
!!$
!!$             ! impaction scavenging removal amount for volume
!!$             scavcoefvol(i,k) = exp( scavimpvol )
!!$             ! impaction scavenging removal amount to number
!!$             scavcoefnum(i,k) = exp( scavimpnum )
!!$
!!$             ! scavcoef = impaction scav rate (1/h) for precip = 1 mm/h
!!$             ! scavcoef = impaction scav rate (1/s) for precip = pfx_inrain
!!$             ! (scavcoef/3600) = impaction scav rate (1/s) for precip = 1 mm/h
!!$             ! (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
!!$             ! impactrate = (scavcoef/3600) * (pfx_inrain*3600)
!!$          else
!!$             scavcoefvol(i,k) = 0._r8
!!$             scavcoefnum(i,k) = 0._r8
!!$          end if
!!$
!!$       end do
!!$    end do
!!$
!!$    return
!!$  end subroutine carma_aero_bcscavcoef_get

  !===============================================================================
  subroutine calc_1_impact_rate(             &
       dg0, sigmag, rhoaero, temp, press, &
       scavratenum, scavratevol, lunerr )
    !
    !   routine computes a single impaction scavenging rate
    !    for precipitation rate of 1 mm/h
    !
    !   dg0 = geometric mean diameter of aerosol number size distrib. (for CARMA it is the dry radius) (cm)
    !   sigmag = geometric standard deviation of size distrib.
    !   rhoaero = density of aerosol particles (g/cm^3)
    !   temp = temperature (K)
    !   press = pressure (dyne/cm^2)
    !   scavratenum = number scavenging rate (1/h)
    !   scavratevol = volume or mass scavenging rate (1/h)
    !   lunerr = logical unit for error message
    !
    use shr_kind_mod, only: r8 => shr_kind_r8
    use mo_constants, only: boltz_cgs, pi, rhowater => rhoh2o_cgs, &
                           gravity => gravity_cgs, rgas => rgas_cgs

   implicit none

   !   subr. parameters
   integer lunerr
   real(r8) dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol

   !   local variables
   integer nrainsvmax
   parameter (nrainsvmax=50)
   real(r8) rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
        vfallrainsv(nrainsvmax)

   integer naerosvmax
   parameter (naerosvmax=51)
   real(r8) aaerosv(naerosvmax), &
        ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

   integer i, ja, jr, na, nr
   real(r8) a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
   real(r8) anumsum, avolsum, cair, chi
   real(r8) d, dr, dum, dumfuchs, dx
   real(r8) ebrown, eimpact, eintercept, etotal, freepath
   real(r8) precip, precipmmhr, precipsum
   real(r8) r, rainsweepout, reynolds, rhi, rhoair, rlo, rnumsum
   real(r8) scavsumnum, scavsumnumbb
   real(r8) scavsumvol, scavsumvolbb
   real(r8) schmidt, sqrtreynolds, sstar, stokes, sx
   real(r8) taurelax, vfall, vfallstp
   real(r8) x, xg0, xg3, xhi, xlo, xmuwaterair


   rlo = .005_r8
   rhi = .250_r8
   dr = 0.005_r8
   nr = 1 + nint( (rhi-rlo)/dr )
   if (nr > nrainsvmax) then
      write(lunerr,9110)
      call endrun()
   end if

9110 format( '*** subr. calc_1_impact_rate -- nr > nrainsvmax' )

   precipmmhr = 1.0_r8
   precip = precipmmhr/36000._r8

! if dg0 the diameter, than ag0 equals the radius
   ag0 = dg0/2._r8
  if (sigmag.ne.1._r8) then
   sx = log( sigmag )
   xg0 = log( ag0 )
   xg3 = xg0 + 3._r8*sx*sx

   xlo = xg3 - 4._r8*sx
   xhi = xg3 + 4._r8*sx
   dx = 0.2_r8*sx

   dx = max( 0.2_r8*sx, 0.01_r8 )
   xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
   xhi = xg3 + max( 4._r8*sx, 2._r8*dx )

   na = 1 + nint( (xhi-xlo)/dx )
   if (na > naerosvmax) then
      write(lunerr,9120)
      call endrun()
   end if
  else
   na = 1
   a = ag0
  end if

9120 format( '*** subr. calc_1_impact_rate -- na > naerosvmax' )

   !   air molar density
   cair = press/(rgas*temp)
   !   air mass density
   rhoair = 28.966_r8*cair
   !   molecular freepath
   freepath = 2.8052e-10_r8/cair
   !   air dynamic viscosity
   airdynvisc = 1.8325e-4_r8 * (416.16_r8/(temp+120._r8)) *    &
        ((temp/296.16_r8)**1.5_r8)
   !   air kinemaic viscosity
   airkinvisc = airdynvisc/rhoair
   !   ratio of water viscosity to air viscosity (from Slinn)
   xmuwaterair = 60.0_r8

   !
   !   compute rain drop number concentrations
   !    rrainsv = raindrop radius (cm)
   !    xnumrainsv = raindrop number concentration (#/cm^3)
   !            (number in the bin, not number density)
   !    vfallrainsv = fall velocity (cm/s)
   !
   precipsum = 0._r8
   do i = 1, nr
      r = rlo + (i-1)*dr
      rrainsv(i) = r
      xnumrainsv(i) = exp( -r/2.7e-2_r8 )

      d = 2._r8*r
      if (d <= 0.007_r8) then
         vfallstp = 2.88e5_r8 * d**2._r8
      else if (d <= 0.025_r8) then
         vfallstp = 2.8008e4_r8 * d**1.528_r8
      else if (d <= 0.1_r8) then
         vfallstp = 4104.9_r8 * d**1.008_r8
      else if (d <= 0.25_r8) then
         vfallstp = 1812.1_r8 * d**0.638_r8
      else
         vfallstp = 1069.8_r8 * d**0.235_r8
      end if

      vfall = vfallstp * sqrt(1.204e-3_r8/rhoair)
      vfallrainsv(i) = vfall
      precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
   end do
   precipsum = precipsum*pi*1.333333_r8

   rnumsum = 0._r8
   do i = 1, nr
      xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
      rnumsum = rnumsum + xnumrainsv(i)
   end do

   !
   !   compute aerosol concentrations
   !    aaerosv = particle radius (cm)
   !    fnumaerosv = fraction of total number in the bin (--)
   !    fvolaerosv = fraction of total volume in the bin (--)
   !


   anumsum = 0._r8
   avolsum = 0._r8
   ynumaerosv(:) = 1._r8
   yvolaerosv(:) = 1._r8
   aaerosv(:) = a
  if (na.ne.1) then
   do i = 1, na
      x = xlo + (i-1)*dx
      a = exp( x )
      aaerosv(i) = a
      dum = (x - xg0)/sx
      ynumaerosv(i) = exp( -0.5_r8*dum*dum )
      yvolaerosv(i) = ynumaerosv(i)*1.3333_r8*pi*a*a*a
      anumsum = anumsum + ynumaerosv(i)
      avolsum = avolsum + yvolaerosv(i)
   end do

   do i = 1, na
      ynumaerosv(i) = ynumaerosv(i)/anumsum
      yvolaerosv(i) = yvolaerosv(i)/avolsum
   end do
  end if


   !
   !   compute scavenging
   !
   scavsumnum = 0._r8
   scavsumvol = 0._r8
   !
   !   outer loop for rain drop radius
   !
   jr_loop: do jr = 1, nr

      r = rrainsv(jr)
      vfall = vfallrainsv(jr)

      reynolds = r * vfall / airkinvisc
      sqrtreynolds = sqrt( reynolds )

      !
      !   inner loop for aerosol particle radius
      !
      scavsumnumbb = 0._r8
      scavsumvolbb = 0._r8

      ja_loop: do ja = 1, na

         a = aaerosv(ja)

         chi = a/r

         dum = freepath/a
         dumfuchs = 1._r8 + 1.246_r8*dum + 0.42_r8*dum*exp(-0.87_r8/dum)
         taurelax = 2._r8*rhoaero*a*a*dumfuchs/(9._r8*rhoair*airkinvisc)


         aeromass = 4._r8*pi*a*a*a*rhoaero/3._r8
         aerodiffus = boltz_cgs*temp*taurelax/aeromass

         schmidt = airkinvisc/aerodiffus
         stokes = vfall*taurelax/r

         ebrown = 4._r8*(1._r8 + 0.4_r8*sqrtreynolds*(schmidt**0.3333333_r8)) /  &
              (reynolds*schmidt)

         dum = (1._r8 + 2._r8*xmuwaterair*chi) /         &
              (1._r8 + xmuwaterair/sqrtreynolds)
         eintercept = 4._r8*chi*(chi + dum)

         dum = log( 1._r8 + reynolds )
         sstar = (1.2_r8 + dum/12._r8) / (1._r8 + dum)
         eimpact = 0._r8
         if (stokes > sstar) then
            dum = stokes - sstar
            eimpact = (dum/(dum+0.6666667_r8)) ** 1.5_r8
         end if

         etotal = ebrown + eintercept + eimpact
         etotal = min( etotal, 1.0_r8 )

         rainsweepout = xnumrainsv(jr)*4._r8*pi*r*r*vfall

         scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
         scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

      enddo ja_loop

      scavsumnum = scavsumnum + scavsumnumbb
      scavsumvol = scavsumvol + scavsumvolbb

   enddo jr_loop

   scavratenum = scavsumnum*3600._r8
   scavratevol = scavsumvol*3600._r8

   return
 end subroutine calc_1_impact_rate

  !=============================================================================
  subroutine mmr2vmr_carma(lchnk, vmr, mbar, mw_carma, ncol, im, rmass)
    !-----------------------------------------------------------------
    !   ... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    !st use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !   ... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: rmass(nbins)
    real(r8), intent(in)    :: mw_carma(ncnst_tot)
    real(r8), intent(inout) :: vmr(ncol,pver,ncnst_tot)
    real(r8)                :: vmr_total(ncol,pver)

    !-----------------------------------------------------------------
    !   ... Local variables
    !-----------------------------------------------------------------
    integer :: k, m, mm, l

    do m = 1, nbins
       vmr_total(:ncol,:) = 0._r8
       do l = 1, nspec(m)   ! for each species, not total mmr or number, information of mw are missing
             mm = bin_idx(m, l)
             do k=1,pver
                vmr(:ncol,k,mm) = mbar(:ncol,k) * vmr(:ncol,k,mm) / mw_carma(mm)
             end do
             vmr_total(:ncol,:) = vmr_total(:ncol,:) +  vmr(:ncol,:,mm)
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
    real(r8)                :: vmr_total(ncol,pver)

    !-----------------------------------------------------------------
    !   ... Local variables
    !-----------------------------------------------------------------
    integer :: k, m, mm, l
    !-----------------------------------------------------------------
    !   ... The non-group species
    !-----------------------------------------------------------------
    do m = 1, nbins
       vmr_total(:ncol,:) = 0._r8
       do l = 1, nspec(m)   ! for each species, not total mmr or number, information of mw are missing
             mm = bin_idx(m, l)
             do k=1,pver
                vmr(:ncol,k,mm) = mw_carma(mm) * vmr(:ncol,k,mm) / mbar(:ncol,k)
             end do
             vmr_total(:ncol,:) = vmr_total(:ncol,:) +  vmr(:ncol,:,mm)
       end do
    end do

   end subroutine vmr2mmr_carma

end module aero_model
