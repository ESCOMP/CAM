!===============================================================================
! Modal Aerosol Model
!===============================================================================
module aero_model
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,         only: pcols, pver, pverp
  use phys_control,   only: phys_getopts, cam_physpkg_is
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use perf_mod,       only: t_startf, t_stopf
  use camsrfexch,     only: cam_in_t, cam_out_t
  use aerodep_flx,    only: aerodep_flx_prescribed
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
  use physics_buffer, only: pbuf_get_field, pbuf_get_index, pbuf_set_field
  use physconst,      only: gravit, rair, rhoh2o
  use spmd_utils,     only: masterproc
  use infnan,         only: nan, assignment(=)

  use cam_history,    only: outfld, fieldname_len
  use chem_mods,      only: gas_pcnst, adv_mass
  use mo_tracname,    only: solsym

  use modal_aero_data,only: cnst_name_cw, lptr_so4_cw_amode
  use modal_aero_data,only: ntot_amode, modename_amode, nspec_max
  use modal_aero_data,only: modal_strat_sulfate

  use ref_pres,       only: top_lev => clim_modal_aero_top_lev

  use mo_setsox_cam,          only: setsox, has_sox
  use aerosol_properties_mod, only: aerosol_properties
  use aerosol_state_mod, only: aerosol_state
  use aerosol_instances_mod, only: aerosol_instances_get_props, &
       aerosol_instances_get_state, aerosol_instances_get_num_models

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep     ! aerosol dry deposition and sediment
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: aero_model_surfarea  ! tropopspheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea ! stratospheric aerosol wet surface area for chemistry
  public :: calc_1_impact_rate
  public :: nimptblgrow_mind, nimptblgrow_maxd

  ! Accessor functions
  public ::  get_scavimptblvol, get_scavimptblnum, get_dlndg_nimptblgrow

  ! Misc private data

  ! number of modes
  integer :: nmodes
  integer :: pblh_idx            = 0
  integer :: dgnum_idx           = 0
  integer :: dgnumwet_idx        = 0
  integer :: rate1_cw2pr_st_idx  = 0

  integer :: wetdens_ap_idx      = 0
  integer :: qaerwat_idx         = 0

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: rprddp_idx          = 0
  integer :: rprdsh_idx          = 0
  integer :: nevapr_shcu_idx     = 0
  integer :: nevapr_dpcu_idx     = 0

  integer :: sulfeq_idx = -1

  integer :: nh3_ndx    = 0
  integer :: nh4_ndx    = 0

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
  real(r8) :: dlndg_nimptblgrow
  real(r8),allocatable :: scavimptblnum(:,:)
  real(r8),allocatable :: scavimptblvol(:,:)

  ! for surf_area_dens
  integer,allocatable :: num_idx(:)

  integer :: ndx_h2so4
  character(len=fieldname_len), allocatable :: dgnum_name(:), dgnumwet_name(:)

  ! Namelist variables
  character(len=16) :: drydep_list(pcnst) = ' '
  real(r8)          :: seasalt_emis_scale

  integer, parameter :: max_sad_spec = 16
  character(len=32) :: sad_chem_spec_types(max_sad_spec) = ' '
  character(len=32) :: sad_seasalt_spec_types(max_sad_spec) = ' '
  character(len=32) :: sad_strat_spec_types(max_sad_spec) = ' '

  ! sfc/dm_aer slots mo_usrrxt must reserve beyond the aerosol bins; all modal
  ! surfaces come from the aerosol representation, so no extra slots are needed
  integer, parameter, public :: n_supplemental_sad = 0

  ! Mode types excluded from SAD: primary_carbon should not contribute
  integer, parameter :: num_sad_exclude_modes = 1
  character(len=32), parameter :: sad_exclude_mode_types(num_sad_exclude_modes) = (/ &
     'primary_carbon                  ' /)

  integer :: ndrydep = 0
  integer,allocatable :: drydep_indices(:)
  logical :: drydep_lq(pcnst)

  logical :: modal_accum_coarse_exch = .false.

  class(aerosol_properties), pointer :: aero_props=>null()
  integer :: iaermod_ = -1

  integer :: n_coarse_dust=-1 ! dmleung added n_coarse_dust to determine the index for the
                              ! coarse dust mode for different MAM versions. 29 Oct 2025
  integer :: ncnst_tot = -1
  integer :: chem_map_ndx(gas_pcnst) = -1

contains

  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use aero_wetdep_cam, only: aero_wetdep_readnl
    use dust_model,      only: dust_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! Namelist variables
    character(len=16) :: aer_drydep_list(pcnst) = ' '

    namelist /aerosol_nl/ aer_drydep_list, modal_strat_sulfate, modal_accum_coarse_exch, seasalt_emis_scale, &
       sad_chem_spec_types, sad_seasalt_spec_types, sad_strat_spec_types

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
    call mpibcast(aer_drydep_list,   len(aer_drydep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(modal_strat_sulfate,     1,                       mpilog,  0, mpicom)
    call mpibcast(seasalt_emis_scale, 1,                            mpir8,   0, mpicom)
    call mpibcast(modal_accum_coarse_exch, 1,                       mpilog,  0, mpicom)
    call mpibcast(sad_chem_spec_types,    len(sad_chem_spec_types(1))*max_sad_spec,    mpichar, 0, mpicom)
    call mpibcast(sad_seasalt_spec_types, len(sad_seasalt_spec_types(1))*max_sad_spec, mpichar, 0, mpicom)
    call mpibcast(sad_strat_spec_types,  len(sad_strat_spec_types(1))*max_sad_spec,   mpichar, 0, mpicom)
#endif

    drydep_list = aer_drydep_list

    call aero_wetdep_readnl(nlfile)
    call dust_readnl(nlfile)

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register()
    use modal_aero_data, only: modal_aero_data_reg

    call modal_aero_data_reg()

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d )

    use mo_chem_utls,    only: get_inv_ndx
    use cam_history,     only: addfld, add_default, horiz_only
    use mo_chem_utls,    only: get_spc_ndx
    use modal_aero_data, only: cnst_name_cw
    use modal_aero_data, only: modal_aero_data_init
    use radiative_aerosol,only: rad_aer_get_info
    use dust_model,      only: dust_init, dust_names, dust_active, dust_nbin, dust_nnum
    use seasalt_model,   only: seasalt_init, seasalt_names, seasalt_active,seasalt_nbin
    use aero_wetdep_cam, only: aero_wetdep_init
    use mo_setsox_cam,   only: sox_inti

    use modal_aero_calcsize_cam, only: modal_aero_calcsize_init
    use modal_aero_coag_cam,   only: modal_aero_coag_cam_init
    use aero_deposition_cam, only: aero_deposition_cam_init
    use modal_aero_gasaerexch_cam, only: modal_aero_gasaerexch_cam_init
    use modal_aero_newnuc_cam, only: modal_aero_newnuc_cam_init
    use modal_aero_rename_cam, only: modal_aero_rename_cam_init
    use aerosol_spec_utils,    only: spec_type_in_list

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    character(len=*), parameter :: subrname = 'aero_model_init'
    integer :: m, n, id
    character(len=20) :: dummy

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_chemistry, history_cesm_forcing, history_dust

    integer :: l
    character(len=6) :: test_name
    character(len=64) :: errmes

    character(len=2)  :: unit_basename  ! Units 'kg' or '1'
    integer :: errcode
    character(len=fieldname_len) :: field_name

    character(len=32) :: spec_name
    character(len=32) :: spec_type
    character(len=32) :: mode_type
    integer :: nspec
    integer :: mm
    character(len=32) :: name_a, name_c

    do iaermod_ = 1, aerosol_instances_get_num_models()
       aero_props => aerosol_instances_get_props(iaermod_, 0)
       if (aero_props%model_is('MAM')) exit
    end do
    ncnst_tot = aero_props%ncnst_tot()

    ! aqueous chem initialization
    call sox_inti(aero_props)

    do m = 1,aero_props%nbins()
       do l = 0,aero_props%nspecies(m)
          mm = aero_props%indexer(m,l)
          if (l==0) then
             call aero_props%num_names(m, name_a, name_c)
          else
             call aero_props%mmr_names(m,l, name_a, name_c)
          end if
          chem_map_ndx(mm) = get_spc_ndx( name_a )
       end do
    end do

    dgnum_idx       = pbuf_get_index('DGNUM')
    dgnumwet_idx    = pbuf_get_index('DGNUMWET')
    fracis_idx      = pbuf_get_index('FRACIS')
    prain_idx       = pbuf_get_index('PRAIN')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    rprdsh_idx      = pbuf_get_index('RPRDSH')
    nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
    sulfeq_idx      = pbuf_get_index('MAMH2SO4EQ',errcode)

    call phys_getopts(history_aerosol_out = history_aerosol, &
                      history_chemistry_out=history_chemistry, &
                      history_cesm_forcing_out=history_cesm_forcing, &
                      history_dust_out=history_dust)

    call rad_aer_get_info(0, nmodes=nmodes)

    call modal_aero_data_init(pbuf2d)
    call modal_aero_bcscavcoef_init()

    call modal_aero_rename_cam_init( modal_accum_coarse_exch )
    !   calcsize call must follow rename call
    call modal_aero_calcsize_init( pbuf2d )
    call modal_aero_gasaerexch_cam_init()
    !   coag call must follow gasaerexch call
    call modal_aero_coag_cam_init
    call modal_aero_newnuc_cam_init

    ! call aero_deposition_cam_init only if the user has not specified
    ! prescribed aerosol deposition fluxes
    if (.not.aerodep_flx_prescribed()) then
       call aero_deposition_cam_init(aero_props)
    endif

    call dust_init()
    call seasalt_init(seasalt_emis_scale)

    ndrydep = 0

    count_species: do m = 1,pcnst
       if ( len_trim(drydep_list(m)) /= 0 ) then
          ndrydep = ndrydep+1
       endif
    enddo count_species

    if (ndrydep>0) &
         allocate(drydep_indices(ndrydep))

    do m = 1,ndrydep
       call cnst_get_ind ( drydep_list(m), id, abort=.false. )
       if (id>0) then
          drydep_indices(m) = id
       else
          call endrun(subrname//': invalid drydep species: '//trim(drydep_list(m)) )
       endif

       if (masterproc) then
          write(iulog,*) subrname//': '//drydep_list(m)//' will have drydep applied'
       endif
    enddo

    if (ndrydep>0) then

       dummy = 'RAM1'
       call addfld (dummy,horiz_only, 'A','frac','RAM1')
       if ( history_aerosol ) then
          call add_default (dummy, 1, ' ')
       endif
       dummy = 'airFV'
       call addfld (dummy,horiz_only, 'A','frac','FV')
       if ( history_aerosol ) then
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (dust_active) then
       ! emissions diagnostics ....

       do m = 1, dust_nbin+dust_nnum
          dummy = trim(dust_names(m)) // 'SF'
          call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(dust_names(m))//' dust surface emission')
          if (history_aerosol.or.history_chemistry) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

       dummy = 'DSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol .or. history_dust) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'LND_MBL'
       call addfld (dummy,horiz_only, 'A','frac','Soil erodibility factor')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (seasalt_active) then

       dummy = 'SSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       do m = 1, seasalt_nbin
          dummy = trim(seasalt_names(m)) // 'SF'
          call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(seasalt_names(m))//' seasalt surface emission')
          if (history_aerosol.or.history_chemistry) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

    endif


    ! set flags for drydep tendencies
    drydep_lq(:) = .false.
    do m=1,ndrydep
       id = drydep_indices(m)
       drydep_lq(id) =  .true.
    enddo

    wetdens_ap_idx = pbuf_get_index('WETDENS_AP')
    qaerwat_idx    = pbuf_get_index('QAERWAT')
    pblh_idx       = pbuf_get_index('pblh')

    rate1_cw2pr_st_idx  = pbuf_get_index('RATE1_CW2PR_ST')
    call pbuf_set_field(pbuf2d, rate1_cw2pr_st_idx, 0.0_r8)

    do m = 1,ndrydep

       ! units
       if (drydep_list(m)(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'
       endif

       call addfld (trim(drydep_list(m))//'DDF', horiz_only,  'A',unit_basename//'/m2/s ', &
            trim(drydep_list(m))//' dry deposition flux at bottom (grav + turb)')
       call addfld (trim(drydep_list(m))//'TBF', horiz_only,  'A',unit_basename//'/m2/s',  &
            trim(drydep_list(m))//' turbulent dry deposition flux')
       call addfld (trim(drydep_list(m))//'GVF', horiz_only,  'A',unit_basename//'/m2/s ', &
            trim(drydep_list(m))//' gravitational dry deposition flux')
       call addfld (trim(drydep_list(m))//'DTQ', (/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(drydep_list(m))//' dry deposition')
       call addfld (trim(drydep_list(m))//'DDV', (/ 'lev' /), 'A','m/s',                   &
            trim(drydep_list(m))//' deposition velocity')

       if ( history_aerosol.or.history_chemistry ) then
          call add_default (trim(drydep_list(m))//'DDF', 1, ' ')
       endif
       if ( history_aerosol ) then
          call add_default (trim(drydep_list(m))//'TBF', 1, ' ')
          call add_default (trim(drydep_list(m))//'GVF', 1, ' ')
       endif

    enddo

    do m = 1,gas_pcnst

       if  ( solsym(m)(1:3) == 'num') then
          unit_basename = ' 1'  ! Units 'kg' or '1'
       else
          unit_basename = 'kg'  ! Units 'kg' or '1'
       end if

       call addfld( 'GS_'//trim(solsym(m)),horiz_only, 'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' gas chemistry/wet removal (for gas species)')
       call addfld( 'AQ_'//trim(solsym(m)),horiz_only, 'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' aqueous chemistry (for gas species)')
       if ( history_aerosol ) then
          call add_default( 'AQ_'//trim(solsym(m)), 1, ' ')
       endif

    enddo
    do n = 1,pcnst
       if( .not. (cnst_name_cw(n) == ' ') ) then

          if (cnst_name_cw(n)(1:3) == 'num') then
             unit_basename = ' 1'
          else
             unit_basename = 'kg'
          endif

          call addfld( cnst_name_cw(n),                (/ 'lev' /), 'A', unit_basename//'/kg ',   &
               trim(cnst_name_cw(n))//' in cloud water')
          call addfld (trim(cnst_name_cw(n))//'DDF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' dry deposition flux at bottom (grav + turb)')
          call addfld (trim(cnst_name_cw(n))//'TBF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' turbulent dry deposition flux')
          call addfld (trim(cnst_name_cw(n))//'GVF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(cnst_name_cw(n))//' gravitational dry deposition flux')

          if ( history_aerosol.or. history_chemistry ) then
             call add_default( cnst_name_cw(n), 1, ' ' )
          endif
          if ( history_aerosol ) then
             call add_default (trim(cnst_name_cw(n))//'GVF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'TBF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'DDF', 1, ' ')
          endif
       endif
    enddo

    allocate(dgnum_name(ntot_amode), dgnumwet_name(ntot_amode))
    do n=1,ntot_amode
       dgnum_name(n) = ' '
       dgnumwet_name(n) = ' '
       write(dgnum_name(n),fmt='(a,i1)') 'dgnum',n
       write(dgnumwet_name(n),fmt='(a,i1)') 'dgnumwet',n
       call addfld( dgnum_name(n),    (/ 'lev' /), 'I', 'm', 'Aerosol mode dry diameter' )
       call addfld( dgnumwet_name(n), (/ 'lev' /), 'I', 'm', 'Aerosol mode wet diameter' )
       if ( history_aerosol ) then
          call add_default( dgnum_name(n), 1, ' ' )
          call add_default( dgnumwet_name(n), 1, ' ' )
       endif
       if ( history_cesm_forcing ) then
          call add_default( dgnumwet_name(n), 8, ' ' )
       endif

       if (modal_strat_sulfate) then
          field_name = ' '
          write(field_name,fmt='(a,i1)') 'wtpct_a',n
          call addfld( field_name, (/ 'lev' /), 'I', '%', 'Aerosol mode weight percent H2SO4' )
          if ( history_aerosol ) then
             call add_default (field_name, 0, 'I')
          endif

          field_name = ' '
          write(field_name,fmt='(a,i1)') 'sulfeq_a',n
          call addfld( field_name, (/ 'lev' /), 'I', 'kg/kg', 'H2SO4 equilibrium mixing ratio' )
          if ( history_aerosol ) then
             call add_default (field_name, 0, 'I')
          endif

          field_name = ' '
          write(field_name,fmt='(a,i1)') 'sulden_a',n
          call addfld( field_name, (/ 'lev' /), 'I', 'g/cm3', 'Sulfate aerosol particle mass density' )
          if ( history_aerosol ) then
             call add_default (field_name, 0, 'I')
          endif

       end if
    end do

    ndx_h2so4 = get_spc_ndx('H2SO4')
    nh3_ndx = get_spc_ndx('NH3')
    nh4_ndx = get_spc_ndx('NH4')

    allocate(num_idx(ntot_amode))
    num_idx = -1

    ! for aero_model_surfarea called from mo_usrrxt
    do l=1,ntot_amode
       test_name = ' '
       write(test_name,fmt='(a5,i1)') 'num_a',l
       num_idx(l) = get_spc_ndx( trim(test_name) )
       if (num_idx(l) < 0) then
          write(errmes,fmt='(a,i1)') 'usrrxt_inti: cannot find MAM num_idx ',l
          write(iulog,*) errmes
          call endrun(errmes)
       endif
    end do

    ! determine coarse dust mode number
    do n = 1,nmodes
       call rad_aer_get_info(0, n, mode_type=mode_type, nspec=nspec)
       if (mode_type=='coarse' .or. mode_type=='coarse_dust') then
          n_coarse_dust = n
       end if
    enddo

    if (masterproc) then
       write(iulog,*) 'SAD chemistry spec_types:'
       do l = 1, max_sad_spec
          if (len_trim(sad_chem_spec_types(l)) > 0) then
             write(iulog,*) '  ', trim(sad_chem_spec_types(l))
          end if
       end do
       write(iulog,*) 'SAD seasalt spec_types:'
       do l = 1, max_sad_spec
          if (len_trim(sad_seasalt_spec_types(l)) > 0) then
             write(iulog,*) '  ', trim(sad_seasalt_spec_types(l))
          end if
       end do
       write(iulog,*) 'SAD stratospheric spec_types:'
       do l = 1, max_sad_spec
          if (len_trim(sad_strat_spec_types(l)) > 0) then
             write(iulog,*) '  ', trim(sad_strat_spec_types(l))
          end if
       end do
    end if

    if (has_sox) then
       do m = 1, ntot_amode

          l = lptr_so4_cw_amode(m)
          if (l > 0) then
             call addfld (&
                  trim(cnst_name_cw(l))//'AQSO4',horiz_only,  'A','kg/m2/s', &
                  trim(cnst_name_cw(l))//' aqueous phase chemistry')
             call addfld (&
                  trim(cnst_name_cw(l))//'AQH2SO4',horiz_only,  'A','kg/m2/s', &
                  trim(cnst_name_cw(l))//' aqueous phase chemistry')
             if ( history_aerosol ) then
                call add_default (trim(cnst_name_cw(l))//'AQSO4', 1, ' ')
                call add_default (trim(cnst_name_cw(l))//'AQH2SO4', 1, ' ')
             endif
          end if

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

    use dust_sediment_mod, only: dust_sediment_tend
    use aero_drydep_core,  only: modal_aero_depvel_part, calcram
    use mo_drydep,         only: n_land_type, fraction_landuse
    use physconst,         only: pi, boltz
    use modal_aero_data,   only: qqcw_get_field
    use modal_aero_data,   only: cnst_name_cw
    use modal_aero_data,   only: alnsg_amode
    use modal_aero_data,   only: sigmag_amode
    use modal_aero_data,   only: nspec_amode
    use modal_aero_data,   only: numptr_amode
    use modal_aero_data,   only: numptrcw_amode
    use modal_aero_data,   only: lmassptr_amode
    use modal_aero_data,   only: lmassptrcw_amode
    use aero_deposition_cam,only: aero_deposition_cam_setdry

  ! args
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt             ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  ! local vars
    real(r8), pointer :: landfrac(:) ! land fraction
    real(r8), pointer :: icefrac(:)  ! ice fraction
    real(r8), pointer :: ocnfrac(:)  ! ocean fraction
    real(r8), pointer :: fvin(:)     !
    real(r8), pointer :: ram1in(:)   ! for dry dep velocities from land model for progseasalts

    real(r8) :: fv(pcols)            ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)          ! for dry dep velocities, from land modified over ocean & ice

    integer :: lchnk                   ! chunk identifier
    integer :: ncol                    ! number of atmospheric columns
    integer :: jvlc                    ! index for last dimension of vlc_xxx arrays
    integer :: lphase                  ! index for interstitial / cloudborne aerosol
    integer :: lspec                   ! index for aerosol number / chem-mass / water-mass
    integer :: m                       ! aerosol mode index
    integer :: mm                      ! tracer index
    integer :: i

    real(r8) :: rho(pcols,pver)      ! air density in kg/m3
    real(r8) :: sflx(pcols)          ! deposition flux
    real(r8) :: dep_trb(pcols)       !kg/m2/s
    real(r8) :: dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8) :: pvmzaer(pcols,pverp) ! sedimentation velocity in Pa
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species

    real(r8) :: rad_drop(pcols,pver)
    real(r8) :: dens_drop(pcols,pver)
    real(r8) :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)

    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity
    real(r8) :: aerdepdryis(pcols,pcnst)  ! aerosol dry deposition (interstitial)
    real(r8) :: aerdepdrycw(pcols,pcnst)  ! aerosol dry deposition (cloud water)
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: dgncur_awet(:,:,:)
    real(r8), pointer :: wetdens(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)

    logical :: aspherical

    character(len=512) :: errmsg_local
    integer :: errflg_local

    landfrac => cam_in%landfrac(:)
    icefrac  => cam_in%icefrac(:)
    ocnfrac  => cam_in%ocnfrac(:)
    fvin     => cam_in%fv(:)
    ram1in   => cam_in%ram1(:)

    lchnk = state%lchnk
    ncol  = state%ncol

    ! calc ram and fv over ocean and sea ice ...
    call calcram( ncol,landfrac,icefrac,ocnfrac,obklen,&
                  ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
                  state%pdel(:,pver),fvin,fv,rair,gravit)

    call outfld( 'airFV', fv(:), pcols, lchnk )
    call outfld( 'RAM1', ram1(:), pcols, lchnk )

    ! note that tendencies are not only in sfc layer (because of sedimentation)
    ! and that ptend is updated within each subroutine for different species

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep', lq=drydep_lq)

    call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )

    rho(:ncol,:)=  state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

!
! calc settling/deposition velocities for cloud droplets (and cloud-borne aerosols)
!
! *** mean drop radius should eventually be computed from ndrop and qcldwtr
    rad_drop(:,:) = 5.0e-6_r8
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_r8
    jvlc = 3    ! dmleung: jvlc = 3, moment = 0 => dry dep velocity for number of cloud-borne aerosols
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 0, &
                     pver, top_lev, n_land_type, fraction_landuse(:,:,lchnk), &
                     pi, boltz, gravit, rair)
    jvlc = 4    ! jvlc = 4, moment = 3 => dry dep velocity for vol/mass of cloud-borne aerosols
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 3, &
                     pver, top_lev, n_land_type, fraction_landuse(:,:,lchnk), &
                     pi, boltz, gravit, rair)

    do m = 1, ntot_amode   ! main loop over aerosol modes

       do lphase = 1, 2   ! loop over interstitial / cloud-borne forms

          if (lphase == 1) then   ! interstial aerosol - calc settling/dep velocities of mode

! rad_aer = volume mean wet radius (m)
! dgncur_awet = geometric mean wet diameter for number distribution (m)
             rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,m)   &
                                 *exp(1.5_r8*(alnsg_amode(m)**2))
! dens_aer(1:ncol,:) = wet density (kg/m3)
             dens_aer(1:ncol,:) = wetdens(1:ncol,:,m)
             sg_aer(1:ncol,:) = sigmag_amode(m)

             ! dmleung 20 Oct 2025 ++
             ! dmleung: adding asphericity effect on slowing down gravitational settling velocity
             ! for internally mixed coarse-mode aerosols (Yue Huang et al., 2020)
             ! Huang et al. (2020) showed that aspherical dust has reduced gravitational settling by 15-20 %.
             ! Since (1) MAM modes are internally mixed, and (2) sea spray aerosols are also aspherical,
             ! for now dmleung applies asphericity correction to grav. set. velocity for the whole coarse mode.

             aspherical = (m == n_coarse_dust)

             jvlc = 1   ! dmleung: jvlc = 1, moment = 0 => dry dep velocity for number of interstitial aerosols
             call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  &
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 0, &
                        pver, top_lev, n_land_type, fraction_landuse(:,:,lchnk), &
                        pi, boltz, gravit, rair, aspherical=aspherical)
             jvlc = 2   ! jvlc = 2, moment = 3 => dry dep velocity for vol/mass of interstitial aerosols
             call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  &
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, &
                        pver, top_lev, n_land_type, fraction_landuse(:,:,lchnk), &
                        pi, boltz, gravit, rair, aspherical=aspherical)

          end if

          do lspec = 0, nspec_amode(m)+1   ! loop over number + constituents + water

             if (lspec == 0) then   ! number
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jvlc = 1
                else
                   mm = numptrcw_amode(m)
                   jvlc = 3
                endif
             else if (lspec <= nspec_amode(m)) then   ! non-water mass
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jvlc = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jvlc = 4
                endif
             else   ! water mass
!   bypass dry deposition of aerosol water
                cycle
                if (lphase == 1) then
                   mm = 0
!                  mm = lwaterptr_amode(m)
                   jvlc = 2
                else
                   mm = 0
                   jvlc = 4
                endif
             endif


          if (mm <= 0) cycle

!         if (lphase == 1) then
          if ((lphase == 1) .and. (lspec <= nspec_amode(m))) then
             ptend%lq(mm) = .TRUE.

             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

             call outfld( trim(cnst_name(mm))//'DDV', pvmzaer(:,2:pverp), pcols, lchnk )

             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pdel, &
                     state%q(:,:,mm),  pvmzaer,  ptend%q(:,:,mm), sflx, &
                     pver,             gravit,   errmsg_local,    errflg_local )
                if (errflg_local /= 0) call endrun('aero_model_drydep: '//trim(errmsg_local))

             ! apportion dry deposition into turb and gravitational settling for tapes
             dep_trb = 0._r8
             dep_grv = 0._r8
             do i=1,ncol
                if (vlc_dry(i,pver,jvlc) /= 0._r8) then
                   dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                   dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                end if
             enddo

             call outfld( trim(cnst_name(mm))//'DDF', sflx, pcols, lchnk)
             call outfld( trim(cnst_name(mm))//'TBF', dep_trb, pcols, lchnk )
             call outfld( trim(cnst_name(mm))//'GVF', dep_grv, pcols, lchnk )
             call outfld( trim(cnst_name(mm))//'DTQ', ptend%q(:,:,mm), pcols, lchnk)
             aerdepdryis(:ncol,mm) = sflx(:ncol)

          else if ((lphase == 1) .and. (lspec == nspec_amode(m)+1)) then  ! aerosol water
             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pdel, &
                     qaerwat(:,:,mm),  pvmzaer,  dqdt_tmp(:,:), sflx, &
                     pver,             gravit,   errmsg_local,  errflg_local )
                if (errflg_local /= 0) call endrun('aero_model_drydep: '//trim(errmsg_local))

             ! apportion dry deposition into turb and gravitational settling for tapes
             dep_trb = 0._r8
             dep_grv = 0._r8
             do i=1,ncol
                if (vlc_dry(i,pver,jvlc) /= 0._r8) then
                   dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                   dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                end if
             enddo

             qaerwat(1:ncol,:,mm) = qaerwat(1:ncol,:,mm) + dqdt_tmp(1:ncol,:) * dt

          else  ! lphase == 2
             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)
             fldcw => qqcw_get_field(pbuf, mm,lchnk)

             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pdel, &
                     fldcw(:,:),  pvmzaer,  dqdt_tmp(:,:), sflx, &
                     pver,        gravit,   errmsg_local,  errflg_local )
                if (errflg_local /= 0) call endrun('aero_model_drydep: '//trim(errmsg_local))

             ! apportion dry deposition into turb and gravitational settling for tapes
             dep_trb = 0._r8
             dep_grv = 0._r8
             do i=1,ncol
                if (vlc_dry(i,pver,jvlc) /= 0._r8) then
                   dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                   dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                end if
             enddo

             fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

             call outfld( trim(cnst_name_cw(mm))//'DDF', sflx, pcols, lchnk)
             call outfld( trim(cnst_name_cw(mm))//'TBF', dep_trb, pcols, lchnk )
             call outfld( trim(cnst_name_cw(mm))//'GVF', dep_grv, pcols, lchnk )
             aerdepdrycw(:ncol,mm) = sflx(:ncol)

          endif

          enddo   ! lspec = 0, nspec_amode(m)+1
       enddo   ! lphase = 1, 2
    enddo   ! m = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not.aerodep_flx_prescribed()) then
       call aero_deposition_cam_setdry(aerdepdryis, aerdepdrycw, cam_out)
    endif

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
  ! provides wet tropospheric aerosol surface area info for modal aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  state, relhum, pmid, temp, ltrop, &
                  sfc, dm_aer, sad_trop, reff_trop, sad_ssa )

    use mo_constants, only : pi

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:)
    real(r8), intent(in)    :: relhum(:,:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_trop(:,:)
    real(r8), intent(out)   :: reff_trop(:,:)
    real(r8), intent(out)   :: sad_ssa(:,:)

    ! local vars
    integer :: beglev(pcols)
    integer :: endlev(pcols)
    integer :: lchnk, ncol
    real(r8) :: reff_ssa(pcols,pver)

    class(aerosol_state), pointer :: aero_state

    sad_ssa = -huge(1._r8)

    lchnk = state%lchnk
    ncol = state%ncol

    beglev(:ncol)=ltrop(:ncol)+1
    endlev(:ncol)=pver

    aero_state => aerosol_instances_get_state(iaermod_, 0, lchnk)

    if (len_trim(sad_seasalt_spec_types(1))>0) then
       call aero_state%surf_area_dens(aero_props, sad_seasalt_spec_types, ncol, pver, beglev, endlev, &
            relhum, pmid, temp, pi, sad_ssa, reff_ssa )
    end if

    call aero_state%surf_area_dens(aero_props, sad_chem_spec_types, ncol, pver, beglev, endlev, &
         relhum, pmid, temp, pi, sad_trop, reff_trop, sfc, dm_aer )

  end subroutine aero_model_surfarea

  !-------------------------------------------------------------------------
  ! provides WET stratospheric aerosol surface area info for modal aerosols
  ! if modal_strat_sulfate = TRUE -- called from mo_gas_phase_chemdr
  !-------------------------------------------------------------------------
  subroutine aero_model_strat_surfarea( state, pmid, temp, ltrop, strato_sad, reff_strat )

    use mo_constants, only : pi

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    real(r8), intent(out)   :: strato_sad(:,:)
    real(r8), intent(out)   :: reff_strat(:,:)

    ! local vars
    integer :: i,k, lchnk, ncol

    real(r8) :: relhum(state%ncol,pver)

    class(aerosol_state), pointer :: aero_state
    integer :: beglev(pcols)
    integer :: endlev(pcols)

    reff_strat = 0._r8
    strato_sad = 0._r8

    if (.not.modal_strat_sulfate) return

    relhum = huge(1._r8)

    lchnk = state%lchnk
    ncol = state%ncol
    beglev(:ncol) = top_lev
    endlev(:ncol) = ltrop(:ncol)

    aero_state => aerosol_instances_get_state(iaermod_, 0, lchnk)
    call aero_state%surf_area_dens(aero_props, sad_strat_spec_types, ncol, pver, beglev, endlev, &
         relhum, pmid, temp, pi, strato_sad, reff_strat)

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( state, loffset, ncol, lchnk, troplev, delt, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    use time_manager,          only : get_nstep
    use modal_aero_coag,       only : modal_aero_coag_run
    use modal_aero_gasaerexch, only : modal_aero_gasaerexch_run, modefrm_pcage
    use modal_aero_rename,     only : modal_aero_rename_run
    use modal_aero_rename_cam, only : npair_renamexf, modefrm_renamexf, modetoo_renamexf, &
                                      nspecfrm_renamexf, lspecfrma_renamexf, lspecfrmc_renamexf, &
                                      lspectooa_renamexf, lspectooc_renamexf, &
                                      igrow_shrink_renamexf, ixferable_all_renamexf, &
                                      ixferable_a_renamexf, ixferable_c_renamexf, strat_only_renamexf
    use modal_aero_newnuc,     only : modal_aero_newnuc_run
    use modal_aero_data,       only : cnst_name_cw, qqcw_get_field, &
                                      nsoa, lptr2_soa_a_amode, &
                                      nspec_amode, &
                                      alnsg_amode, voltonumblo_amode, voltonumbhi_amode, &
                                      dgnum_amode, specmw_amode, specdens_amode, &
                                      lmassptr_amode, lmassptrcw_amode, numptr_amode, &
                                      numptrcw_amode, modeptr_accum, modeptr_coarse, &
                                      modeptr_stracoar
    use mo_constants,          only : pi
    use mo_chem_utls,          only : get_spc_ndx
    use constituents,          only : pcnst, cnst_name
    use physconst,             only : mwdry

    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    type(physics_state),target,intent(in)  :: state ! Physics state variables
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    integer,  intent(in) :: troplev(pcols)
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

    integer :: n, m
    integer :: i,k,l
    integer :: nstep

    real(r8) :: del_h2so4_aeruptk(ncol,pver)

    real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:), wetdens(:,:,:)
    real(r8), pointer :: pblh(:)                    ! pbl height (m)

    real(r8), dimension(ncol) :: wrk
    character(len=32)         :: name
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)

    real(r8) ::  aqso4(ncol,ntot_amode)               ! aqueous phase chemistry
    real(r8) ::  aqh2so4(ncol,ntot_amode)             ! aqueous phase chemistry
    real(r8) ::  aqso4_h2o2(ncol)                     ! SO4 aqueous phase chemistry due to H2O2
    real(r8) ::  aqso4_o3(ncol)                       ! SO4 aqueous phase chemistry due to O3
    real(r8) ::  xphlwc(ncol,pver)                    ! pH value multiplied by lwc
    real(r8) ::  nh3_beg(ncol,pver)
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: sulfeq(:,:,:)

    real(r8) :: qqcw(ncol,pver,ncnst_tot)

    integer :: mm
    character(len=32) :: specname
    class(aerosol_state), pointer :: aero_state

    ! Local arrays for refactored gasaerexch call
    real(r8) :: dqdt_gaex(ncol,pver,gas_pcnst)
    real(r8) :: dqdt_gaex_conden(ncol,pver,gas_pcnst)  ! conden-only snapshot (pre-rename) for diagnostics
    real(r8) :: dqdt_rnpos_unused(ncol,pver,gas_pcnst) ! required rename output, unused by CAM
    logical  :: dotend_gaex(gas_pcnst)
    real(r8) :: dqqcwdt_gaex(ncol,pver,gas_pcnst)
    logical  :: dotendrn(gas_pcnst), dotendqqcwrn(gas_pcnst)
    logical  :: is_dorename_atik, dorename_atik(ncol,pver)
    integer, parameter :: jsrflx_gaexch = 1
    integer, parameter :: jsrflx_rename = 2
    integer, parameter :: nsrflx = 2
    real(r8) :: qsrflx(pcols,gas_pcnst,nsrflx)
    real(r8) :: qqcwsrflx(pcols,gas_pcnst,nsrflx)
    real(r8) :: qsrflx_gaexch_out(ncol,gas_pcnst)     ! column-integrated gaexch source/sink from the scheme

    ! Local arrays for refactored newnuc call
    real(r8) :: dqdt_nnuc(ncol,pver,gas_pcnst)
    logical  :: dotend_nnuc(gas_pcnst)
    real(r8) :: qsrflx_nnuc(pcols,gas_pcnst,1)        ! column-integrated nucleation source/sink

    ! Local arrays for refactored coag call
    ! dqdt_coag is diagnostic-only, vmr is updated in-place
    real(r8) :: dqdt_coag(ncol,pver,gas_pcnst)
    logical  :: dotend_coag(gas_pcnst)
    real(r8) :: qsrflx_coag(pcols)                    ! column-integrated coagulation source/sink
    character(len=fieldname_len+3) :: fieldname
    integer  :: jac, jsrf, jsoa, lb
    logical  :: use_sulfeq

    character(len=512) :: errmsg_local
    integer :: errflg_local
    ! Zero-initialized dummy array for intent(in) placeholders (e.g. sulfeq
    ! when use_sulfeq=.false.)
    real(r8) :: dummy_3d(ncol,pver,ntot_amode)

    ! SOA condensation/evaporation diagnostics
    real(r8) :: qconff(pcols,pver),qevapff(pcols,pver)
    real(r8) :: qconbb(pcols,pver),qevapbb(pcols,pver)
    real(r8) :: qconbg(pcols,pver),qevapbg(pcols,pver)
    real(r8) :: qcon(pcols,pver),qevap(pcols,pver)
    real(r8) :: dqdt_soa_val
    integer  :: l_soa

    aero_state => aerosol_instances_get_state(iaermod_, 0, lchnk)

    dummy_3d(:,:,:) = 0.0_r8
!
! ... initialize nh3
!
    if ( nh3_ndx > 0 ) then
      nh3_beg = vmr(1:ncol,:,nh3_ndx)
    end if
!

    call pbuf_get_field(pbuf, dgnum_idx,      dgnum)
    call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens )
    call pbuf_get_field(pbuf, pblh_idx,       pblh)

    do n=1,ntot_amode
       call outfld(dgnum_name(n), dgnum(1:ncol,1:pver,n), ncol, lchnk )
       call outfld(dgnumwet_name(n), dgnumwet(1:ncol,1:pver,n), ncol, lchnk )
    end do

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    nstep = get_nstep()

    ! calculate tendency due to gas phase chemistry and processes
    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0._r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m)*adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'GS_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

!
! Aerosol processes ...
!
    call qqcw2vmr( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)
    dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)

    ! aqueous chemistry ...

    if( has_sox ) then

       ! Temperary code to map cloud-borne aerosol VMRs to aerosol only array (qqcw)
       ! needed for setsox interface.  When refactoring aero_model_gasaerexch
       ! with modal_aero_gasaerexch_sub, this mapping should go away.
       do m = 1,aero_props%nbins()
          do l = 0,aero_props%nspecies(m)
             mm = aero_props%indexer(m,l)
             qqcw(:,:,mm) = vmrcw(:,:,chem_map_ndx(mm))
          end do
       end do

       call setsox( aero_state, state, &
              pbuf,     &
              ncol,     &
              delt,     &
              pmid,     &
              pdel,     &
              tfld,     &
              mbar,     &
              cwat,     &
              cldfr,    &
              cldnum,   &
              invariants, &
              qqcw,     &
              vmr,      &
              xphlwc,   &
              aqso4,    &
              aqh2so4,  &
              aqso4_h2o2, &
              aqso4_o3  &
              )

       ! Map back to all-species chemistry VMR array
       do m = 1,aero_props%nbins()
          do l = 0,aero_props%nspecies(m)
             mm = aero_props%indexer(m,l)
             vmrcw(:,:,chem_map_ndx(mm)) = qqcw(:,:,mm)
          end do
       end do

       do n = 1, ntot_amode
          l = lptr_so4_cw_amode(n)
          if (l > 0) then
             call outfld( trim(cnst_name_cw(l))//'AQSO4',   aqso4(:ncol,n),   ncol, lchnk)
             call outfld( trim(cnst_name_cw(l))//'AQH2SO4', aqh2so4(:ncol,n), ncol, lchnk)
          end if
       end do

       call outfld( 'AQSO4_H2O2', aqso4_h2o2(:ncol), ncol, lchnk)
       call outfld( 'AQSO4_O3',   aqso4_o3(:ncol),   ncol, lchnk)
       call outfld( 'XPH_LWC',    xphlwc(:ncol,:),   ncol, lchnk )

    endif

    ! Tendency due to aqueous chemistry
    dvmrdt = (vmr - dvmrdt) / delt
    dvmrcwdt = (vmrcw - dvmrcwdt) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0._r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m) * adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'AQ_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    call t_startf('modal_gas-aer_exchng')

    if ( sulfeq_idx>0 ) then
       call pbuf_get_field( pbuf, sulfeq_idx, sulfeq )
       use_sulfeq = .true.
    else
       nullify( sulfeq )
       use_sulfeq = .false.
    endif

    ! Call portable gasaerexch_run to get tendencies
    dqdt_gaex(:,:,:) = 0.0_r8
    dotend_gaex(:) = .false.
    if (use_sulfeq) then
       call modal_aero_gasaerexch_run(                       &
            ncol      = ncol,                                &
            pver      = pver,                                &
            deltat    = delt,                                &
            top_lev   = top_lev,                             &
            loffset   = loffset,                             &
            t         = tfld(:ncol,:),                       &
            pmid      = pmid(:ncol,:),                       &
            pdel      = pdel(:ncol,:),                       &
            gravit    = gravit,                              &
            troplev   = troplev(:ncol),                      &
            dgncur_a  = dgnum(:ncol,:,:),                    &
            dgncur_awet = dgnumwet(:ncol,:,:),               &
            use_sulfeq = .true.,                             &
            sulfeq    = sulfeq(:ncol,:,:),                   &
            num_q     = gas_pcnst,                           &
            q         = vmr(:ncol,:,:),                      &
            dqdt      = dqdt_gaex,                           &
            dotend    = dotend_gaex,                         &
            qsrflx_gaexch = qsrflx_gaexch_out,               &
            errmsg    = errmsg_local,                        &
            errflg    = errflg_local)
    else
       call modal_aero_gasaerexch_run(                       &
            ncol      = ncol,                                &
            pver      = pver,                                &
            deltat    = delt,                                &
            top_lev   = top_lev,                             &
            loffset   = loffset,                             &
            t         = tfld(:ncol,:),                       &
            pmid      = pmid(:ncol,:),                       &
            pdel      = pdel(:ncol,:),                       &
            gravit    = gravit,                              &
            troplev   = troplev(:ncol),                      &
            dgncur_a  = dgnum(:ncol,:,:),                    &
            dgncur_awet = dgnumwet(:ncol,:,:),               &
            use_sulfeq = .false.,                            &
            sulfeq    = dummy_3d,                            &
            num_q     = gas_pcnst,                           &
            q         = vmr(:ncol,:,:),                      &
            dqdt      = dqdt_gaex,                           &
            dotend    = dotend_gaex,                         &
            qsrflx_gaexch = qsrflx_gaexch_out,               &
            errmsg    = errmsg_local,                        &
            errflg    = errflg_local)
    end if

    if (errflg_local /= 0) then
       call endrun('aero_model_gasaerexch: ' // trim(errmsg_local))
    end if

    ! Save conden-only tendencies before modal_aero_rename_run
    ! adds mode-transfer tendencies into dqdt_gaex
    ! since they are used by _sfgaex1 and SOA cond/evap diagnostics
    dqdt_gaex_conden(:,:,:) = dqdt_gaex(:,:,:)

    if (ndx_h2so4 > 0) then
       ! Save h2so4 vmr before applying tendencies.
       ! del_h2so4_aeruptk is computed by (vmr_after - vmr_before) after
       ! the apply loop for bfb.
       !
       ! A clearer formulation is
       ! del_h2so4_aeruptk(1:ncol,:) = dqdt_gaex(1:ncol,:,ndx_h2so4) * delt
       ! but is not bfb and the difference propagates down to newnuc.
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4)
    else
       del_h2so4_aeruptk(:,:) = 0.0_r8
    end if

    ! Call portable mode merging (renaming):
    dqqcwdt_gaex(:,:,:) = 0.0_r8
    dotendrn(:) = .false.
    dotendqqcwrn(:) = .false.
    dorename_atik(1:ncol,:) = .true.
    is_dorename_atik = .true.

    ! zero to pcols:
    qsrflx(:,:,:)    = 0.0_r8
    qqcwsrflx(:,:,:) = 0.0_r8
    call modal_aero_rename_run( &
       ncol                    = ncol,                          &
       loffset                 = loffset,                       &
       deltat                  = delt,                          &
       pdel                    = pdel(:ncol,:),                 &
       troplev                 = troplev(:ncol),                &
       dotendrn                = dotendrn,                      &
       q                       = vmr(:ncol,:,:),                &
       dqdt                    = dqdt_gaex(:ncol,:,:),          &
       dqdt_other              = dvmrdt(:ncol,:,:),             &
       dotendqqcwrn            = dotendqqcwrn,                  &
       qqcw                    = vmrcw(:ncol,:,:),              &
       dqqcwdt                 = dqqcwdt_gaex(:ncol,:,:),       &
       dqqcwdt_other           = dvmrcwdt(:ncol,:,:),           &
       is_dorename_atik        = is_dorename_atik,              &
       dorename_atik           = dorename_atik(:ncol,:),        &
       jsrflx_rename           = jsrflx_rename,                 &
       qsrflx                  = qsrflx(:ncol,:,:),             &
       qqcwsrflx               = qqcwsrflx(:ncol,:,:),          &
       dqdt_rnpos              = dqdt_rnpos_unused,             &
       ntot_amode              = ntot_amode,                    &
       npair_renamexf          = npair_renamexf,                &
       modefrm_renamexf        = modefrm_renamexf,              &
       modetoo_renamexf        = modetoo_renamexf,              &
       nspecfrm_renamexf       = nspecfrm_renamexf,             &
       lspecfrma_renamexf      = lspecfrma_renamexf,            &
       lspecfrmc_renamexf      = lspecfrmc_renamexf,            &
       lspectooa_renamexf      = lspectooa_renamexf,            &
       lspectooc_renamexf      = lspectooc_renamexf,            &
       alnsg_amode             = alnsg_amode,                   &
       voltonumblo_amode       = voltonumblo_amode,             &
       voltonumbhi_amode       = voltonumbhi_amode,             &
       dgnum_amode             = dgnum_amode,                   &
       nspec_amode             = nspec_amode,                   &
       specmw_amode            = specmw_amode,                  &
       specdens_amode          = specdens_amode,                &
       lmassptr_amode          = lmassptr_amode,                &
       lmassptrcw_amode        = lmassptrcw_amode,              &
       numptr_amode            = numptr_amode,                  &
       numptrcw_amode          = numptrcw_amode,                &
       pi                      = pi,                            &
       modeptr_accum           = modeptr_accum,                 &
       modeptr_coarse          = modeptr_coarse,                &
       modeptr_stracoar        = modeptr_stracoar,              &
       igrow_shrink_renamexf   = igrow_shrink_renamexf,         &
       ixferable_all_renamexf  = ixferable_all_renamexf,        &
       ixferable_a_renamexf    = ixferable_a_renamexf,          &
       ixferable_c_renamexf    = ixferable_c_renamexf,          &
       strat_only_renamexf     = strat_only_renamexf,           &
       modal_accum_coarse_exch = modal_accum_coarse_exch,       &
       pver                    = pver,                          &
       gravit                  = gravit,                        &
       errmsg                  = errmsg_local,                  &
       errflg                  = errflg_local)

    if (errflg_local /= 0) then
       call endrun('aero_model_gasaerexch (rename): ' // trim(errmsg_local))
    end if

    ! Apply tendencies to vmr and vmrcw
    do l = 1, gas_pcnst
       if ( dotend_gaex(l) .or. dotendrn(l) ) then
          do k = top_lev, pver
             do i = 1, ncol
                vmr(i,k,l) = vmr(i,k,l) + dqdt_gaex(i,k,l)*delt
             end do
          end do
       end if
       if ( dotendqqcwrn(l) ) then
          do k = top_lev, pver
             do i = 1, ncol
                vmrcw(i,k,l) = vmrcw(i,k,l) + dqqcwdt_gaex(i,k,l)*delt
             end do
          end do
       end if
    end do

    ! Get del_h2so4_aeruptk = vmr_after - vmr_before:
    if (ndx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - del_h2so4_aeruptk(1:ncol,:)
    end if

    ! Diagnostics: column tendencies for gas-aerosol exchange and renaming.
    ! The gaexch column source/sink (qsrflx, jsrflx_gaexch)
    ! is accumulated inside the scheme:
    ! the per-mode and primary-carbon-aging contributions
    ! are sum'd term by term for bfb.
    qsrflx(:ncol,:,jsrflx_gaexch) = qsrflx_gaexch_out(:ncol,:)

    ! Output history fields
    do l = 1, gas_pcnst
       lb = l + loffset
       do jsrf = 1, 2
          do jac = 1, 2
             if (jac == 1) then
               if (jsrf == jsrflx_gaexch) then
                  if ( .not. dotend_gaex(l) ) cycle
                  fieldname = trim(cnst_name(lb)) // '_sfgaex1'
               else if (jsrf == jsrflx_rename) then
                  if ( .not. dotendrn(l) ) cycle
                  fieldname = trim(cnst_name(lb)) // '_sfgaex2'
               else
                  cycle
               end if
               do i = 1, ncol
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do
               call outfld( fieldname, qsrflx(:,l,jsrf), pcols, lchnk )
             else
               if (jsrf == jsrflx_gaexch) then
                  cycle
               else if (jsrf == jsrflx_rename) then
                  if ( .not. dotendqqcwrn(l) ) cycle
                  fieldname = trim(cnst_name_cw(lb)) // '_sfgaex2'
               else
                  cycle
               end if
               do i = 1, ncol
                  qqcwsrflx(i,l,jsrf) = qqcwsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do
               call outfld( fieldname, qqcwsrflx(:,l,jsrf), pcols, lchnk )
             end if
          end do ! jac = ...
       end do ! jsrf = ...
    end do ! l = ...

    ! SOA condensation/evaporation diagnostics
    ! Reconstruct from the pre-rename conden tendencies (dqdt_gaex_conden).
    ! NOTE: for the accumulation mode this is not exactly b4b with the original,
    ! which used the per-mode conden tendency dqdt_soa(n,jsoa);
    ! the species-indexed tendency here also absorbs primary-carbon-aged SOA
    qconff(:,:) = 0.0_r8
    qevapff(:,:) = 0.0_r8
    qconbb(:,:) = 0.0_r8
    qevapbb(:,:) = 0.0_r8
    qconbg(:,:) = 0.0_r8
    qevapbg(:,:) = 0.0_r8
    qcon(:,:) = 0.0_r8
    qevap(:,:) = 0.0_r8

    do n = 1, ntot_amode
       do jsoa = 1, nsoa
          l_soa = lptr2_soa_a_amode(n,jsoa) - loffset
          if ((l_soa <= 0) .or. (l_soa > gas_pcnst)) cycle
          ! Skip pcage from-mode: only accumulated for ido_soaa==1
          if (modefrm_pcage > 0 .and. n == modefrm_pcage) cycle
          do k = top_lev, pver
             do i = 1, ncol
                dqdt_soa_val = dqdt_gaex_conden(i,k,l_soa)
                if (nsoa.eq.15) then !check for current SOA package
                   if(jsoa.ge.1.and.jsoa.le.5) then ! Fossil SOA species
                      if (dqdt_soa_val.ge.0.0_r8) then
                         qconff(i,k)=qconff(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      elseif(dqdt_soa_val.lt.0.0_r8) then
                         qevapff(i,k)=qevapff(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      endif

                   elseif(jsoa.ge.6.and.jsoa.le.10) then ! Biomass SOA species
                      if (dqdt_soa_val.ge.0.0_r8) then
                         qconbb(i,k)=qconbb(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      elseif(dqdt_soa_val.lt.0.0_r8) then
                         qevapbb(i,k)=qevapbb(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      endif

                   elseif(jsoa.ge.11.and.jsoa.le.15) then ! Biogenic SOA species
                      if (dqdt_soa_val.ge.0.0_r8) then
                         qconbg(i,k)=qconbg(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      elseif(dqdt_soa_val.lt.0.0_r8) then
                         qevapbg(i,k)=qevapbg(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      endif

                   endif ! jsoa
                endif !nsoa
                if (nsoa.eq.5) then !check for current SOA package
                      if (dqdt_soa_val.ge.0.0_r8) then
                         qcon(i,k)=qcon(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      elseif(dqdt_soa_val.lt.0.0_r8) then
                         qevap(i,k)=qevap(i,k)+dqdt_soa_val*(adv_mass(l_soa)/mwdry)
                      endif
                endif !nsoa
             end do ! i
          end do ! k
       end do ! jsoa
    end do ! n

    if (nsoa.eq.5) then
       call outfld(trim('qcon_gaex'), qcon(:,:), pcols, lchnk )
       call outfld(trim('qevap_gaex'), qevap(:,:), pcols, lchnk )
    endif
    if (nsoa.eq.15) then
       call outfld(trim('qconff_gaex'), qconff(:,:), pcols, lchnk )
       call outfld(trim('qevapff_gaex'), qevapff(:,:), pcols, lchnk )
       call outfld(trim('qconbb_gaex'), qconbb(:,:), pcols, lchnk )
       call outfld(trim('qevapbb_gaex'), qevapbb(:,:), pcols, lchnk )
       call outfld(trim('qconbg_gaex'), qconbg(:,:), pcols, lchnk )
       call outfld(trim('qevapbg_gaex'), qevapbg(:,:), pcols, lchnk )
    endif

    call t_stopf('modal_gas-aer_exchng')

    call t_startf('modal_nucl')

    ! do aerosol nucleation (new particle formation)
    qsrflx_nnuc(:,:,:) = 0.0_r8
    call modal_aero_newnuc_run(                          &
         ncol      = ncol,                               &
         pver      = pver,                               &
         top_lev   = top_lev,                            &
         num_q     = gas_pcnst,                          &
         loffset   = loffset,                            &
         deltat    = delt,                               &
         t         = tfld(:ncol,:),                      &
         pmid      = pmid(:ncol,:),                      &
         pdel      = pdel(:ncol,:),                      &
         zm        = zm(:ncol,:),                        &
         pblh      = pblh(:ncol),                        &
         qv        = qh2o(:ncol,:),                      &
         cld       = cldfr(:ncol,:),                     &
         q         = vmr(:ncol,:,:),                     &
         gravit    = gravit,                             &
         del_h2so4_gasprod = del_h2so4_gasprod(:ncol,:), &
         del_h2so4_aeruptk = del_h2so4_aeruptk(:ncol,:), &
         dqdt      = dqdt_nnuc,                          &
         dotend    = dotend_nnuc,                        &
         qsrflx    = qsrflx_nnuc(:ncol,:,:),             &
         errmsg    = errmsg_local,                       &
         errflg    = errflg_local )

    if (errflg_local /= 0) then
       call endrun('aero_model_gasaerexch (newnuc): ' // trim(errmsg_local))
    end if

    ! Apply nucleation tendencies to vmr (was applied in place by the scheme)
    do l = 1, gas_pcnst
       if ( dotend_nnuc(l) ) then
          do k = top_lev, pver
             do i = 1, ncol
                vmr(i,k,l) = vmr(i,k,l) + dqdt_nnuc(i,k,l)*delt
             end do
          end do
       end if
    end do

    ! do history file column-tendency fields
    do l = 1, gas_pcnst
       if ( .not. dotend_nnuc(l) ) cycle
       lb = l + loffset
       do i = 1, ncol
          qsrflx_nnuc(i,l,1) = qsrflx_nnuc(i,l,1)*(adv_mass(l)/mwdry)
       end do
       fieldname = trim(cnst_name(lb)) // '_sfnnuc1'
       call outfld( fieldname, qsrflx_nnuc(:,l,1), pcols, lchnk )
    end do ! l = ...

    call t_stopf('modal_nucl')

    call t_startf('modal_coag')

    ! do aerosol coagulation
    ! vmr is updated in place by the scheme.
    ! dqdt_coag is returned for the history diagnostics only
    ! (dqdt*delt is not bit-identical to the stored change so it cannot be applied directly)
    call modal_aero_coag_run(                            &
         ncol      = ncol,                               &
         pver      = pver,                               &
         top_lev   = top_lev,                            &
         loffset   = loffset,                            &
         nstep     = nstep,                              &
         deltat_main = delt,                             &
         t         = tfld(:ncol,:),                      &
         pmid      = pmid(:ncol,:),                      &
         q         = vmr(:ncol,:,:),                     &
         dgncur_a  = dgnum(:ncol,:,:),                   &
         dgncur_awet = dgnumwet(:ncol,:,:),              &
         wetdens_a = wetdens(:ncol,:,:),                 &
         dqdt      = dqdt_coag,                          &
         dotend    = dotend_coag,                        &
         errmsg    = errmsg_local,                       &
         errflg    = errflg_local )

    if (errflg_local /= 0) then
       call endrun('aero_model_gasaerexch (coag): ' // trim(errmsg_local))
    end if

    ! do history file column-tendency fields
    do l = 1, gas_pcnst
       if ( .not. dotend_coag(l) ) cycle
       lb = l + loffset

       qsrflx_coag(:) = 0.0_r8
       do k = top_lev, pver
       do i = 1, ncol
          qsrflx_coag(i) = qsrflx_coag(i) + dqdt_coag(i,k,l)*pdel(i,k)
       end do
       end do
       qsrflx_coag(:) = qsrflx_coag(:)*(adv_mass(l)/(gravit*mwdry))
       fieldname = trim(cnst_name(lb)) // '_sfcoag1'
       call outfld( fieldname, qsrflx_coag, pcols, lchnk )
    end do ! l = ...

    call t_stopf('modal_coag')

    call vmr2qqcw( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    ! diagnostics for cloud-borne aerosols...
    do n = 1,pcnst
       fldcw => qqcw_get_field(pbuf,n,lchnk,errorhandle=.true.)
       if(associated(fldcw)) then
          call outfld( cnst_name_cw(n), fldcw(:,:), pcols, lchnk )
       endif
    end do
!
! ... put missing NH3 into NH4
!
    if ( nh3_ndx > 0 .and. nh4_ndx > 0 ) then
      vmr(1:ncol,:,nh4_ndx) = vmr(1:ncol,:,nh4_ndx) + (nh3_beg-vmr(1:ncol,:,nh3_ndx))
      vmr(1:ncol,:,nh4_ndx) = max(0._r8,vmr(1:ncol,:,nh4_ndx))
    end if

  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
    use seasalt_model, only: seasalt_emis, seasalt_names, seasalt_indices, seasalt_active,seasalt_nbin
    use dust_model,    only: dust_emis, dust_names, dust_indices, dust_active,dust_nbin, dust_nnum
    use physics_types, only: physics_state

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    ! local vars

    integer :: lchnk, ncol
    integer :: m, mm
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: sflx(pcols)   ! accumulate over all bins for output

    lchnk = state%lchnk
    ncol = state%ncol

    if (dust_active) then

       call dust_emis( ncol, lchnk, cam_in%dstflx, cam_in%cflx, soil_erod_tmp )

       ! some dust emis diagnostics ...
       sflx(:)=0._r8
       do m=1,dust_nbin+dust_nnum
          mm = dust_indices(m)
          if (m<=dust_nbin) sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,mm),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif

    if (seasalt_active) then
       sflx(:)=0._r8

       call seasalt_emis( state%u(:ncol,pver), state%v(:ncol,pver), state%zm(:ncol,pver), &
                          cam_in%sst, cam_in%ocnfrac, ncol, cam_in%cflx )

       do m=1,seasalt_nbin
          mm = seasalt_indices(m)
          sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(seasalt_names(m))//'SF',cam_in%cflx(:,mm),pcols,lchnk)
       enddo
       call outfld('SSTSFMBL',sflx(:),pcols,lchnk)
    endif

  end subroutine aero_model_emissions

  !===============================================================================
  ! private methods

  !===============================================================================
  !===============================================================================
  subroutine modal_aero_bcscavcoef_init
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Computes lookup table for aerosol impaction/interception scavenging rates
    !
    ! Authors: R. Easter
    !
    !-----------------------------------------------------------------------

    use shr_kind_mod,    only: r8 => shr_kind_r8
    use modal_aero_data
    use cam_abortutils,  only: endrun

    implicit none


    !   local variables
    integer nnfit_maxd
    parameter (nnfit_maxd=27)

    integer i, jgrow, jdens, jpress, jtemp, mode, nnfit
    integer lunerr

    real(r8) dg0, dg0_cgs, press, &
         rhodryaero, rhowetaero, rhowetaero_cgs, rmserr, &
         scavratenum, scavratevol, sigmag,                &
         temp, wetdiaratio, wetvolratio
    real(r8) aafitnum(1), xxfitnum(1,nnfit_maxd), yyfitnum(nnfit_maxd)
    real(r8) aafitvol(1), xxfitvol(1,nnfit_maxd), yyfitvol(nnfit_maxd)

    allocate(scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode))
    allocate(scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode))

    lunerr = 6
    dlndg_nimptblgrow = log( 1.25_r8 )

    modeloop: do mode = 1, ntot_amode

       sigmag = sigmag_amode(mode)

       rhodryaero = specdens_amode(1,mode)

       growloop: do jgrow = nimptblgrow_mind, nimptblgrow_maxd

          wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
          dg0 = dgnum_amode(mode)*wetdiaratio

          wetvolratio = exp( jgrow*dlndg_nimptblgrow*3._r8 )
          rhowetaero = 1.0_r8 + (rhodryaero-1.0_r8)/wetvolratio
          rhowetaero = min( rhowetaero, rhodryaero )

          !
          !   compute impaction scavenging rates at 1 temp-press pair and save
          !
          nnfit = 0

          temp = 273.16_r8
          press = 0.75e6_r8   ! dynes/cm2
          rhowetaero = rhodryaero

          dg0_cgs = dg0*1.0e2_r8   ! m to cm
          rhowetaero_cgs = rhowetaero*1.0e-3_r8   ! kg/m3 to g/cm3
          call calc_1_impact_rate( &
               dg0_cgs, sigmag, rhowetaero_cgs, temp, press, &
               scavratenum, scavratevol, lunerr )

          nnfit = nnfit + 1
          if (nnfit > nnfit_maxd) then
             write(lunerr,9110)
             call endrun()
          end if
9110      format( '*** subr. modal_aero_bcscavcoef_init -- nnfit too big' )

          xxfitnum(1,nnfit) = 1._r8
          yyfitnum(nnfit) = log( scavratenum )

          xxfitvol(1,nnfit) = 1._r8
          yyfitvol(nnfit) = log( scavratevol )

5900      continue

          !
          ! skip mlinfit stuff because scav table no longer has dependencies on
          !    air temp, air press, and particle wet density
          ! just load the log( scavrate--- ) values
          !
          !!
          !!   do linear regression
          !!	log(scavrate) = a1 + a2*log(wetdens)
          !!
          !	call mlinft( xxfitnum, yyfitnum, aafitnum, nnfit, 1, 1, rmserr )
          !	call mlinft( xxfitvol, yyfitvol, aafitvol, nnfit, 1, 1, rmserr )
          !
          !	scavimptblnum(jgrow,mode) = aafitnum(1)
          !	scavimptblvol(jgrow,mode) = aafitvol(1)

          scavimptblnum(jgrow,mode) = yyfitnum(1)
          scavimptblvol(jgrow,mode) = yyfitvol(1)

       enddo growloop
    enddo modeloop
    return
  end subroutine modal_aero_bcscavcoef_init


  !===============================================================================
  subroutine modal_aero_bcscavcoef_get( m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol )

    use modal_aero_data
    !-----------------------------------------------------------------------
    implicit none

    integer,intent(in) :: m, ncol
    logical,intent(in):: isprx(pcols,pver)
    real(r8), intent(in) :: dgn_awet(pcols,pver,ntot_amode)
    real(r8), intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)

    integer i, k, jgrow
    real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum


    do k = 1, pver
       do i = 1, ncol

          ! do only if no precip
          if ( isprx(i,k) ) then
             !
             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)

             dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)

             if ((dumdgratio >= 0.99_r8) .and. (dumdgratio <= 1.01_r8)) then
                scavimpvol = scavimptblvol(0,m)
                scavimpnum = scavimptblnum(0,m)
             else
                xgrow = log( dumdgratio ) / dlndg_nimptblgrow
                jgrow = int( xgrow )
                if (xgrow < 0._r8) jgrow = jgrow - 1
                if (jgrow < nimptblgrow_mind) then
                   jgrow = nimptblgrow_mind
                   xgrow = jgrow
                else
                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
                end if

                dumfhi = xgrow - jgrow
                dumflo = 1._r8 - dumfhi

                scavimpvol = dumflo*scavimptblvol(jgrow,m) + &
                     dumfhi*scavimptblvol(jgrow+1,m)
                scavimpnum = dumflo*scavimptblnum(jgrow,m) + &
                     dumfhi*scavimptblnum(jgrow+1,m)

             end if

             ! impaction scavenging removal amount for volume
             scavcoefvol(i,k) = exp( scavimpvol )
             ! impaction scavenging removal amount to number
             scavcoefnum(i,k) = exp( scavimpnum )

             ! scavcoef = impaction scav rate (1/h) for precip = 1 mm/h
             ! scavcoef = impaction scav rate (1/s) for precip = pfx_inrain
             ! (scavcoef/3600) = impaction scav rate (1/s) for precip = 1 mm/h
             ! (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
             ! impactrate = (scavcoef/3600) * (pfx_inrain*3600)
          else
             scavcoefvol(i,k) = 0._r8
             scavcoefnum(i,k) = 0._r8
          end if

       end do
    end do

    return
  end subroutine modal_aero_bcscavcoef_get

  !===============================================================================
	subroutine calc_1_impact_rate(             &
     		dg0, sigmag, rhoaero, temp, press, &
     		scavratenum, scavratevol, lunerr )
   !
   !   routine computes a single impaction scavenging rate
   !	for precipitation rate of 1 mm/h
   !
   !   dg0 = geometric mean diameter of aerosol number size distrib. (cm)
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

   ag0 = dg0/2._r8
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
   !	rrainsv = raindrop radius (cm)
   !	xnumrainsv = raindrop number concentration (#/cm^3)
   !		(number in the bin, not number density)
   !	vfallrainsv = fall velocity (cm/s)
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
   !	aaerosv = particle radius (cm)
   !	fnumaerosv = fraction of total number in the bin (--)
   !	fvolaerosv = fraction of total volume in the bin (--)
   !
   anumsum = 0._r8
   avolsum = 0._r8
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
  !=============================================================================
  subroutine qqcw2vmr(lchnk, vmr, mbar, ncol, im, pbuf)
    use modal_aero_data, only : qqcw_get_field
    use physics_buffer, only : physics_buffer_desc
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)

    do m=1,gas_pcnst
       if( adv_mass(m) /= 0._r8 ) then
          fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
          if(associated(fldcw)) then
             do k=1,pver
                vmr(:ncol,k,m) = mbar(:ncol,k) * fldcw(:ncol,k) / adv_mass(m)
             end do
          else
             vmr(:,:,m) = 0.0_r8
          end if
       end if
    end do
  end subroutine qqcw2vmr


  !=============================================================================
  !=============================================================================
  subroutine vmr2qqcw( lchnk, vmr, mbar, ncol, im, pbuf )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    use m_spc_id
    use chem_mods,       only : adv_mass, gas_pcnst
    use modal_aero_data, only : qqcw_get_field
    use physics_buffer,  only : physics_buffer_desc

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)
    !-----------------------------------------------------------------
    !	... The non-group species
    !-----------------------------------------------------------------
    do m = 1,gas_pcnst
       fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
       if( adv_mass(m) /= 0._r8 .and. associated(fldcw)) then
          do k = 1,pver
             fldcw(:ncol,k) = adv_mass(m) * vmr(:ncol,k,m) / mbar(:ncol,k)
          end do
       end if
    end do

  end subroutine vmr2qqcw

  function get_dlndg_nimptblgrow() result (dlndg_nimptblgrow_ret)
    real(r8) ::  dlndg_nimptblgrow_ret
    dlndg_nimptblgrow_ret =  dlndg_nimptblgrow
  end function get_dlndg_nimptblgrow

  function get_scavimptblvol() result (scavimptblvol_ret)
    real(r8) :: scavimptblvol_ret(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)
    scavimptblvol_ret = scavimptblvol
  end function get_scavimptblvol

  function get_scavimptblnum() result (scavimptblnum_ret)
    real(r8) :: scavimptblnum_ret(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)
    scavimptblnum_ret = scavimptblnum
  end function get_scavimptblnum

end module aero_model
