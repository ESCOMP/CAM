module chemistry

!---------------------------------------------------------------------------------
! "Interactive" gas phase module
!---------------------------------------------------------------------------------

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use ppgrid,           only : pcols, pver, begchunk, endchunk
  use physconst,        only : gravit
  use constituents,     only : pcnst, cnst_add, cnst_name, cnst_fixed_ubc
  use chem_mods,        only : gas_pcnst
  use cam_history,      only : fieldname_len
  use physics_types,    only : physics_state, physics_ptend, physics_ptend_init
  use spmd_utils,       only : masterproc
  use cam_logfile,      only : iulog
  use mo_gas_phase_chemdr, only : map2chm
  use shr_megan_mod,    only : shr_megan_mechcomps, shr_megan_mechcomps_n 
  use srf_field_check,  only : active_Fall_flxvoc
  use tracer_data,      only : MAXTRCRS
  use gcr_ionization,   only : gcr_ionization_readnl, gcr_ionization_init, gcr_ionization_adv
  use epp_ionization,   only : epp_ionization_readnl, epp_ionization_adv
  use mo_apex,          only : mo_apex_readnl
  use ref_pres,         only : do_molec_diff, ptop_ref
  use phys_control,     only : waccmx_is   ! WACCM-X switch query function

  implicit none
  private
  save

!---------------------------------------------------------------------------------
! Public interfaces
!---------------------------------------------------------------------------------
  public :: chem_is                        ! identify which chemistry is being used
  public :: chem_register                  ! register consituents
  public :: chem_readnl                    ! read chem namelist 
  public :: chem_is_active                 ! returns true
  public :: chem_implements_cnst           ! returns true if consituent is implemented by this package
  public :: chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public :: chem_init                      ! initialize (history) variables
  public :: chem_timestep_init             ! per timestep initializations
  public :: chem_timestep_tend             ! interface to tendency computation
  public :: chem_final
  public :: chem_write_restart
  public :: chem_read_restart
  public :: chem_init_restart
  public :: chem_emissions

  integer, public :: imozart = -1       ! index of 1st constituent
  
  ! Namelist variables
  
  ! control
  
  integer :: chem_freq = 1 ! time steps

  ! ghg

  character(len=shr_kind_cl) :: bndtvg = ' ' ! pathname for greenhouse gas loss rate
  character(len=shr_kind_cl) :: h2orates = ' ' ! pathname for greenhouse gas (lyman-alpha H2O loss)

  ! lightning

  real(r8)           :: lght_no_prd_factor = 1._r8

  ! photolysis

  logical            :: xactive_prates = .false.
  character(len=shr_kind_cl) :: rsf_file = 'rsf_file'
  character(len=shr_kind_cl) :: exo_coldens_file = ''
  character(len=shr_kind_cl) :: tuv_xsect_file = 'tuv_xsect_file'
  character(len=shr_kind_cl) :: o2_xsect_file = 'o2_xsect_file'
  character(len=shr_kind_cl) :: xs_coef_file = 'xs_coef_file'
  character(len=shr_kind_cl) :: xs_short_file = 'xs_short_file'
  character(len=shr_kind_cl) :: xs_long_file = 'xs_long_file'
  character(len=shr_kind_cl) :: electron_file = 'electron_file'
  character(len=shr_kind_cl) :: euvac_file = 'NONE'

  ! solar / geomag data

  character(len=shr_kind_cl) :: photon_file = 'photon_file'

  ! dry dep
  
  character(len=shr_kind_cl) :: depvel_file = 'depvel_file'
  character(len=shr_kind_cl) :: depvel_lnd_file = 'depvel_lnd_file'
  character(len=shr_kind_cl) :: clim_soilw_file = 'clim_soilw_file'
  character(len=shr_kind_cl) :: season_wes_file = 'season_wes_file'

  ! emis

  character(len=shr_kind_cl) :: airpl_emis_file = '' ! airplane emissions
  character(len=shr_kind_cl) :: srf_emis_specifier(pcnst) = ''
  character(len=shr_kind_cl) :: ext_frc_specifier(pcnst) = ''

  character(len=24)  :: srf_emis_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer            :: srf_emis_cycle_yr  = 0
  integer            :: srf_emis_fixed_ymd = 0
  integer            :: srf_emis_fixed_tod = 0

  character(len=24)  :: ext_frc_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer            :: ext_frc_cycle_yr  = 0
  integer            :: ext_frc_fixed_ymd = 0
  integer            :: ext_frc_fixed_tod = 0

  ! fixed stratosphere
  
  character(len=shr_kind_cl) :: fstrat_file = 'fstrat_file'
  character(len=16)  :: fstrat_list(pcnst)  = ''

! for linoz
  character(len=shr_kind_cl) :: chlorine_loading_file = ''
  character(len=8)   :: chlorine_loading_type = 'SERIAL' ! "FIXED" or "SERIAL"
  integer            :: chlorine_loading_fixed_ymd = 0         ! YYYYMMDD for "FIXED" type
  integer            :: chlorine_loading_fixed_tod = 0         ! seconds of day for "FIXED" type

!---------------------------------------------------------------------------------
! dummy values for specific heats at constant pressure
!---------------------------------------------------------------------------------
  real(r8), parameter   :: cptmp = 666._r8

  character(len=fieldname_len) :: srcnam(gas_pcnst) ! names of source/sink tendencies

  integer :: ixcldliq, ixcldice                     ! indicies of liquid and ice cloud water
  integer :: ndx_cld
  integer :: ndx_cmfdqr
  integer :: ndx_nevapr
  integer :: ndx_prain
  integer :: ndx_cldtop
  integer :: h2o_ndx
  integer :: ixndrop             ! cloud droplet number index
  integer :: ndx_pblh
  integer :: ndx_fsds

  logical :: ghg_chem = .false.      ! .true. => use ghg chem package
  logical :: chem_step = .true.
  logical :: is_active = .false.

  character(len=32) :: chem_name = 'NONE'
  logical :: chem_rad_passive = .false.
  
  ! for MEGAN emissions
  integer, allocatable :: megan_indices_map(:) 
  real(r8),allocatable :: megan_wght_factors(:)

  logical :: chem_use_chemtrop = .false.

!================================================================================================
contains
!================================================================================================

logical function chem_is (name)
   use phys_control,     only : cam_chempkg_is

   character(len=*), intent(in) :: name
   chem_is = cam_chempkg_is(name)

end function chem_is

!================================================================================================

  subroutine chem_register
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents and physics buffer fields
! 
!-----------------------------------------------------------------------

    use mo_sim_dat,          only : set_sim_dat
    use chem_mods,           only : gas_pcnst, adv_mass
    use mo_tracname,         only : solsym
    use mo_chem_utls,        only : get_spc_ndx
    use short_lived_species, only : slvd_index, short_lived_map=>map, register_short_lived_species
    use cfc11star,           only : register_cfc11star
    use mo_photo,            only : photo_register
    use mo_aurora,           only : aurora_register
    use aero_model,          only : aero_model_register

    implicit none

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer  :: m, n                                ! tracer index
    real(r8) :: qmin                                ! min value
    logical  :: ic_from_cam2                        ! wrk variable for initial cond input
    logical  :: has_fixed_ubc                       ! wrk variable for upper bndy cond
    logical  :: has_fixed_ubflx                     ! wrk variable for upper bndy flux
    integer  :: ch4_ndx, n2o_ndx, o3_ndx 
    integer  :: cfc11_ndx, cfc12_ndx, o2_1s_ndx, o2_1d_ndx, o2_ndx
    integer  :: n_ndx, no_ndx, h_ndx, h2_ndx, o_ndx, e_ndx, np_ndx
    integer  :: op_ndx, o1d_ndx, n2d_ndx, nop_ndx, n2p_ndx, o2p_ndx
    integer  :: hf_ndx, f_ndx

    character(len=128) :: lng_name                  ! variable long name
    logical :: cam_outfld
    character(len=128) :: mixtype
    character(len=128) :: molectype
    integer :: islvd

!-----------------------------------------------------------------------
! Set the simulation chemistry variables
!-----------------------------------------------------------------------
    call set_sim_dat

    o3_ndx    = get_spc_ndx('O3')
    ch4_ndx   = get_spc_ndx('CH4')
    n2o_ndx   = get_spc_ndx('N2O')

    cfc11_ndx = get_spc_ndx('CFC11')
    cfc12_ndx = get_spc_ndx('CFC12')
    o2_1s_ndx = get_spc_ndx('O2_1S')
    o2_1d_ndx = get_spc_ndx('O2_1D')
    o2_ndx    = get_spc_ndx('O2')
    n_ndx     = get_spc_ndx('N')
    no_ndx    = get_spc_ndx('NO')
    h_ndx     = get_spc_ndx('H')
    h2_ndx    = get_spc_ndx('H2')
    o_ndx     = get_spc_ndx('O')
    e_ndx     = get_spc_ndx('e')
    np_ndx    = get_spc_ndx('Np')
    op_ndx    = get_spc_ndx('Op')
    o1d_ndx   = get_spc_ndx('O1D')
    n2d_ndx   = get_spc_ndx('N2D')
    n2p_ndx   = get_spc_ndx('N2p')
    nop_ndx   = get_spc_ndx('NOp')
    h2o_ndx   = get_spc_ndx('H2O')
    o2p_ndx   = get_spc_ndx('O2p')

    f_ndx     = get_spc_ndx('F')
    hf_ndx    = get_spc_ndx('HF')


    !-----------------------------------------------------------------------
    ! Set names of diffused variable tendencies and declare them as history variables
    !-----------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    ! For WACCM-X, change variable has_fixed_ubc from .true. to .false. which is a flag
    ! used later to check for a fixed upper boundary condition for species. 
    !----------------------------------------------------------------------------------
     do m = 1,gas_pcnst
     ! setting of these variables is for registration of transported species
       ic_from_cam2  = .true.
       has_fixed_ubc = .false.
       has_fixed_ubflx = .false.
       lng_name      = trim( solsym(m) )
       molectype = 'minor'

       qmin = 1.e-36_r8
       
       if ( lng_name(1:5) .eq. 'num_a' ) then ! aerosol number density
          qmin = 1.e-5_r8
       else if ( m == o3_ndx ) then
          qmin = 1.e-12_r8
       else if ( m == ch4_ndx ) then
          qmin = 1.e-12_r8
          if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
            has_fixed_ubc = .false.   ! diffusive equilibrium at UB
          else
            has_fixed_ubc = .true.
          endif
       else if ( m == n2o_ndx ) then
          qmin = 1.e-15_r8
       else if( m == cfc11_ndx .or. m == cfc12_ndx ) then
          qmin = 1.e-20_r8
       else if( m == o2_1s_ndx .or. m == o2_1d_ndx ) then
          ic_from_cam2 = .false.
          if( m == o2_1d_ndx ) then
             lng_name = 'O2(1-delta)'
          else
             lng_name = 'O2(1-sigma)'
          end if
       else if ( m==o2_ndx .or. m==n_ndx .or. m==no_ndx .or. m==h_ndx .or. m==h2_ndx .or. m==o_ndx .or. m==hf_ndx &
               .or. m==f_ndx ) then
         if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
           has_fixed_ubc = .false.   ! diffusive equilibrium at UB
           if ( m == h_ndx ) has_fixed_ubflx = .true. ! fixed flux value for H at UB
           if ( m == o2_ndx .or. m == o_ndx ) molectype = 'major'
         else
           has_fixed_ubc = .true.
         endif
       else if( m == e_ndx ) then
          lng_name = 'electron concentration'
       else if( m == np_ndx ) then
          lng_name = 'N+'
       else if( m == op_ndx ) then
          lng_name = 'O+'
       else if( m == o1d_ndx ) then
          lng_name = 'O(1D)'
       else if( m == n2d_ndx ) then
          lng_name = 'N(2D)'
       else if( m == o2p_ndx ) then
          lng_name = 'O2+'
       else if( m == n2p_ndx ) then
          lng_name = 'N2+'
       else if( m == nop_ndx ) then
          lng_name = 'NO+'
       else if( m == h2o_ndx ) then
          map2chm(1) = m
          cycle
       endif

       cam_outfld=.false.
       is_active = .true.
       mixtype = 'dry'

       islvd = slvd_index(solsym(m))

       if ( islvd > 0 ) then
          short_lived_map(islvd) = m
       else
          call cnst_add( solsym(m), adv_mass(m), cptmp, qmin, n, readiv=ic_from_cam2, cam_outfld=cam_outfld, &
                         mixtype=mixtype, molectype=molectype, fixed_ubc=has_fixed_ubc, fixed_ubflx=has_fixed_ubflx, &
                         longname=trim(lng_name) )

          if( imozart == -1 ) then
             imozart = n
          end if
          map2chm(n) = m
       endif

    end do
    
    call register_short_lived_species()
    call register_cfc11star()

    if ( waccmx_is('ionosphere') ) then 
       call photo_register()
       call aurora_register()
    endif
    
    ! add fields to pbuf needed by aerosol models
    call aero_model_register()

  end subroutine chem_register

!================================================================================================
  
  subroutine chem_readnl(nlfile)

    ! Read chem namelist group.

    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    use linoz_data,       only: linoz_data_defaultopts,  linoz_data_setopts
    use tracer_cnst,      only: tracer_cnst_defaultopts, tracer_cnst_setopts
    use tracer_srcs,      only: tracer_srcs_defaultopts, tracer_srcs_setopts
    use aero_model,       only: aero_model_readnl
    use dust_model,       only: dust_readnl
    use gas_wetdep_opts,  only: gas_wetdep_readnl
    use upper_bc,         only: ubc_defaultopts, ubc_setopts
    use mo_drydep,        only: drydep_srf_file
    use noy_ubc,          only: noy_ubc_readnl
    use mo_sulf,          only: sulf_readnl
    use species_sums_diags,only: species_sums_readnl
    use ocean_emis,       only: ocean_emis_readnl

    ! args

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    ! linoz data
    character(len=shr_kind_cl) :: linoz_data_file               ! prescribed data file
    character(len=shr_kind_cl) :: linoz_data_filelist           ! list of prescribed data files (series of files)
    character(len=shr_kind_cl) :: linoz_data_path               ! absolute path of prescribed data files 
    character(len=24)  :: linoz_data_type               ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)
    logical            :: linoz_data_rmfile             ! remove data file from local disk (default .false.)
    integer            :: linoz_data_cycle_yr
    integer            :: linoz_data_fixed_ymd
    integer            :: linoz_data_fixed_tod

    ! trop_mozart prescribed constituent concentratons
    character(len=shr_kind_cl) :: tracer_cnst_file              ! prescribed data file
    character(len=shr_kind_cl) :: tracer_cnst_filelist          ! list of prescribed data files (series of files)
    character(len=shr_kind_cl) :: tracer_cnst_datapath          ! absolute path of prescribed data files 
    character(len=24)  :: tracer_cnst_type              ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)
    character(len=shr_kind_cl) :: tracer_cnst_specifier(MAXTRCRS) ! string array where each 
    logical            :: tracer_cnst_rmfile            ! remove data file from local disk (default .false.)
    integer            :: tracer_cnst_cycle_yr
    integer            :: tracer_cnst_fixed_ymd
    integer            :: tracer_cnst_fixed_tod

    ! trop_mozart prescribed constituent sourrces/sinks
    character(len=shr_kind_cl) :: tracer_srcs_file              ! prescribed data file
    character(len=shr_kind_cl) :: tracer_srcs_filelist          ! list of prescribed data files (series of files)
    character(len=shr_kind_cl) :: tracer_srcs_datapath          ! absolute path of prescribed data files 
    character(len=24)  :: tracer_srcs_type              ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)
    character(len=shr_kind_cl) :: tracer_srcs_specifier(MAXTRCRS) ! string array where each 
    logical            :: tracer_srcs_rmfile            ! remove data file from local disk (default .false.)
    integer            :: tracer_srcs_cycle_yr
    integer            :: tracer_srcs_fixed_ymd
    integer            :: tracer_srcs_fixed_tod

    ! Upper boundary conditions
    character(len=shr_kind_cl) :: tgcm_ubc_file
    integer            :: tgcm_ubc_cycle_yr
    integer            :: tgcm_ubc_fixed_ymd
    integer            :: tgcm_ubc_fixed_tod
    character(len=32)  :: tgcm_ubc_data_type
    character(len=shr_kind_cl) :: snoe_ubc_file
    ! Upper boundary conditions
    real(r8)           :: t_pert_ubc   ! temperature perturbation at ubc
    real(r8)           :: no_xfac_ubc  ! no multiplicative factor at ubc

    namelist /chem_inparm/ chem_freq, airpl_emis_file, &
         euvac_file, photon_file, electron_file, &
         depvel_file, xs_coef_file, xs_short_file, &
         exo_coldens_file, tuv_xsect_file, o2_xsect_file, &
         xs_long_file, rsf_file, &
         lght_no_prd_factor, xactive_prates, &
         depvel_lnd_file, clim_soilw_file, season_wes_file, drydep_srf_file, &
         srf_emis_type, srf_emis_cycle_yr, srf_emis_fixed_ymd, srf_emis_fixed_tod, srf_emis_specifier,  &
         fstrat_file, fstrat_list, &
         ext_frc_specifier, ext_frc_type, ext_frc_cycle_yr, ext_frc_fixed_ymd, ext_frc_fixed_tod

    namelist /chem_inparm/ chem_rad_passive

    ! ghg chem

    namelist /chem_inparm/ bndtvg, h2orates, ghg_chem

    ! linoz inputs

    namelist /chem_inparm/ &
         linoz_data_file, linoz_data_filelist, linoz_data_path, &
         linoz_data_type, &
         linoz_data_rmfile, linoz_data_cycle_yr, linoz_data_fixed_ymd, linoz_data_fixed_tod
    namelist /chem_inparm/ &
         chlorine_loading_file, chlorine_loading_type, chlorine_loading_fixed_ymd, chlorine_loading_fixed_tod

    ! prescribed chem tracers

    namelist /chem_inparm/ &
         tracer_cnst_file, tracer_cnst_filelist, tracer_cnst_datapath, &
         tracer_cnst_type, tracer_cnst_specifier, &
         tracer_srcs_file, tracer_srcs_filelist, tracer_srcs_datapath, &
         tracer_srcs_type, tracer_srcs_specifier, &
         tracer_cnst_rmfile, tracer_cnst_cycle_yr, tracer_cnst_fixed_ymd, tracer_cnst_fixed_tod, &
         tracer_srcs_rmfile, tracer_srcs_cycle_yr, tracer_srcs_fixed_ymd, tracer_srcs_fixed_tod 
    
    ! upper boundary conditions
    namelist /chem_inparm/ tgcm_ubc_file, tgcm_ubc_data_type, tgcm_ubc_cycle_yr, tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod, &
                           snoe_ubc_file, t_pert_ubc, no_xfac_ubc

    ! tropopause level control
    namelist /chem_inparm/ chem_use_chemtrop

    ! get the default settings

    call linoz_data_defaultopts( &
         linoz_data_file_out      = linoz_data_file,      &
         linoz_data_filelist_out  = linoz_data_filelist,  &
         linoz_data_path_out      = linoz_data_path,      &
         linoz_data_type_out      = linoz_data_type,      &
         linoz_data_rmfile_out    = linoz_data_rmfile,    &
         linoz_data_cycle_yr_out  = linoz_data_cycle_yr,  &
         linoz_data_fixed_ymd_out = linoz_data_fixed_ymd, &
         linoz_data_fixed_tod_out = linoz_data_fixed_tod  ) 
    call tracer_cnst_defaultopts( &
         tracer_cnst_file_out      = tracer_cnst_file,      &
         tracer_cnst_filelist_out  = tracer_cnst_filelist,  &
         tracer_cnst_datapath_out  = tracer_cnst_datapath,  &
         tracer_cnst_type_out      = tracer_cnst_type,      &
         tracer_cnst_specifier_out = tracer_cnst_specifier, &
         tracer_cnst_rmfile_out    = tracer_cnst_rmfile,    &
         tracer_cnst_cycle_yr_out  = tracer_cnst_cycle_yr,  &
         tracer_cnst_fixed_ymd_out = tracer_cnst_fixed_ymd, &
         tracer_cnst_fixed_tod_out = tracer_cnst_fixed_tod  ) 
    call tracer_srcs_defaultopts( &
         tracer_srcs_file_out      = tracer_srcs_file,      &
         tracer_srcs_filelist_out  = tracer_srcs_filelist,  &
         tracer_srcs_datapath_out  = tracer_srcs_datapath,  &
         tracer_srcs_type_out      = tracer_srcs_type,      &
         tracer_srcs_specifier_out = tracer_srcs_specifier, &
         tracer_srcs_rmfile_out    = tracer_srcs_rmfile,    &
         tracer_srcs_cycle_yr_out  = tracer_srcs_cycle_yr,  &
         tracer_srcs_fixed_ymd_out = tracer_srcs_fixed_ymd, &
         tracer_srcs_fixed_tod_out = tracer_srcs_fixed_tod  )

    ! Upper boundary conditions
    call ubc_defaultopts( &
         snoe_ubc_file_out =snoe_ubc_file, &
         t_pert_ubc_out    =t_pert_ubc, &
         no_xfac_ubc_out   =no_xfac_ubc, &
         tgcm_ubc_file_out      = tgcm_ubc_file, &
         tgcm_ubc_data_type_out = tgcm_ubc_data_type, &
         tgcm_ubc_cycle_yr_out  = tgcm_ubc_cycle_yr, &
         tgcm_ubc_fixed_ymd_out = tgcm_ubc_fixed_ymd, &
         tgcm_ubc_fixed_tod_out = tgcm_ubc_fixed_tod )

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'chem_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, chem_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('chem_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables

    ! control

    call mpibcast (chem_freq,         1,                               mpiint,  0, mpicom)

    call mpibcast (chem_rad_passive,  1,                               mpilog,  0, mpicom)

    ! ghg

    call mpibcast (ghg_chem,          1,                               mpilog,  0, mpicom)
    call mpibcast (bndtvg,            len(bndtvg),                     mpichar, 0, mpicom)
    call mpibcast (h2orates,          len(h2orates),                   mpichar, 0, mpicom)

    ! lightning

    call mpibcast (lght_no_prd_factor,1,                               mpir8,   0, mpicom)

    ! photolysis

    call mpibcast (rsf_file,          len(rsf_file),                   mpichar, 0, mpicom)
    call mpibcast (exo_coldens_file,  len(exo_coldens_file),           mpichar, 0, mpicom)
    call mpibcast (tuv_xsect_file,    len(tuv_xsect_file),             mpichar, 0, mpicom)
    call mpibcast (o2_xsect_file,     len(o2_xsect_file),              mpichar, 0, mpicom)
    call mpibcast (xs_coef_file,      len(xs_coef_file),               mpichar, 0, mpicom)
    call mpibcast (xs_short_file,     len(xs_short_file),              mpichar, 0, mpicom)
    call mpibcast (xs_long_file,      len(xs_long_file),               mpichar, 0, mpicom)
    call mpibcast (xactive_prates,    1,                               mpilog,  0, mpicom)
    call mpibcast (electron_file,     len(electron_file),              mpichar, 0, mpicom)
    call mpibcast (euvac_file,        len(euvac_file),                 mpichar, 0, mpicom)

    ! solar / geomag data

    call mpibcast (photon_file,       len(photon_file),                mpichar, 0, mpicom)

    ! dry dep

    call mpibcast (depvel_lnd_file,   len(depvel_lnd_file),            mpichar, 0, mpicom)
    call mpibcast (depvel_file,       len(depvel_file),                mpichar, 0, mpicom)
    call mpibcast (clim_soilw_file,   len(clim_soilw_file),            mpichar, 0, mpicom)
    call mpibcast (season_wes_file,   len(season_wes_file),            mpichar, 0, mpicom)
    call mpibcast (drydep_srf_file,   len(drydep_srf_file),            mpichar, 0, mpicom)

    ! emis

    call mpibcast (airpl_emis_file,   len(airpl_emis_file),            mpichar, 0, mpicom)
    call mpibcast (srf_emis_specifier,len(srf_emis_specifier(1))*pcnst,mpichar, 0, mpicom)
    call mpibcast (srf_emis_type,     len(srf_emis_type),              mpichar, 0, mpicom)
    call mpibcast (srf_emis_cycle_yr, 1,                               mpiint,  0, mpicom)
    call mpibcast (srf_emis_fixed_ymd,1,                               mpiint,  0, mpicom)
    call mpibcast (srf_emis_fixed_tod,1,                               mpiint,  0, mpicom)
    call mpibcast (ext_frc_specifier, len(ext_frc_specifier(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast (ext_frc_type,      len(ext_frc_type),               mpichar, 0, mpicom)
    call mpibcast (ext_frc_cycle_yr,  1,                               mpiint,  0, mpicom)
    call mpibcast (ext_frc_fixed_ymd, 1,                               mpiint,  0, mpicom)
    call mpibcast (ext_frc_fixed_tod, 1,                               mpiint,  0, mpicom)


    ! fixed stratosphere

    call mpibcast (fstrat_file,       len(fstrat_file),                mpichar, 0, mpicom)
    call mpibcast (fstrat_list,       len(fstrat_list(1))*pcnst,       mpichar, 0, mpicom)

    ! upper boundary
    call mpibcast (tgcm_ubc_file,      len(tgcm_ubc_file),     mpichar, 0, mpicom)
    call mpibcast (tgcm_ubc_data_type, len(tgcm_ubc_data_type),mpichar, 0, mpicom)
    call mpibcast (tgcm_ubc_cycle_yr,  1,                      mpiint,  0, mpicom)
    call mpibcast (tgcm_ubc_fixed_ymd, 1,                      mpiint,  0, mpicom)
    call mpibcast (tgcm_ubc_fixed_tod, 1,                      mpiint,  0, mpicom)

    call mpibcast (snoe_ubc_file, len(snoe_ubc_file), mpichar, 0, mpicom)
    call mpibcast (t_pert_ubc,    1,                  mpir8,   0, mpicom)
    call mpibcast (no_xfac_ubc,   1,                  mpir8,   0, mpicom)

    ! linoz data

    call mpibcast (linoz_data_file,      len(linoz_data_file),                  mpichar, 0, mpicom)
    call mpibcast (linoz_data_filelist,  len(linoz_data_filelist),              mpichar, 0, mpicom)
    call mpibcast (linoz_data_path,      len(linoz_data_path),              mpichar, 0, mpicom)
    call mpibcast (linoz_data_type,      len(linoz_data_type),                  mpichar, 0, mpicom)
    call mpibcast (linoz_data_rmfile,    1,                                     mpilog, 0,  mpicom)
    call mpibcast (linoz_data_cycle_yr,       1,                                     mpiint, 0,  mpicom)
    call mpibcast (linoz_data_fixed_ymd,       1,                                     mpiint, 0,  mpicom)
    call mpibcast (linoz_data_fixed_tod,       1,                                     mpiint, 0,  mpicom)

    call mpibcast (chlorine_loading_file,len(chlorine_loading_file), mpichar, 0, mpicom)
    call mpibcast (chlorine_loading_type,len(chlorine_loading_type), mpichar, 0, mpicom)
    call mpibcast (chlorine_loading_fixed_ymd, 1,                    mpiint,  0, mpicom)
    call mpibcast (chlorine_loading_fixed_tod, 1,                    mpiint,  0, mpicom)

    ! prescribed chemical tracers

    call mpibcast (tracer_cnst_specifier, len(tracer_cnst_specifier(1))*MAXTRCRS, mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_file,      len(tracer_cnst_file),                  mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_filelist,  len(tracer_cnst_filelist),              mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_datapath,  len(tracer_cnst_datapath),              mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_type,      len(tracer_cnst_type),                  mpichar, 0, mpicom)
    call mpibcast (tracer_cnst_rmfile,    1,                                      mpilog, 0,  mpicom)
    call mpibcast (tracer_cnst_cycle_yr,  1,                                      mpiint, 0,  mpicom)
    call mpibcast (tracer_cnst_fixed_ymd, 1,                                      mpiint, 0,  mpicom)
    call mpibcast (tracer_cnst_fixed_tod, 1,                                      mpiint, 0,  mpicom)

    call mpibcast (tracer_srcs_specifier, len(tracer_srcs_specifier(1))*MAXTRCRS, mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_file,      len(tracer_srcs_file),                  mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_filelist,  len(tracer_srcs_filelist),              mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_datapath,  len(tracer_srcs_datapath),              mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_type,      len(tracer_srcs_type),                  mpichar, 0, mpicom)
    call mpibcast (tracer_srcs_rmfile,    1,                                      mpilog,  0, mpicom)
    call mpibcast (tracer_srcs_cycle_yr,  1,                                      mpiint,  0, mpicom)
    call mpibcast (tracer_srcs_fixed_ymd, 1,                                      mpiint,  0, mpicom)
    call mpibcast (tracer_srcs_fixed_tod, 1,                                      mpiint,  0, mpicom)

    call mpibcast (chem_use_chemtrop,1,                                    mpilog,  0, mpicom)

#endif

    ! set the options

   call linoz_data_setopts( &
        linoz_data_file_in      = linoz_data_file,      &
        linoz_data_filelist_in  = linoz_data_filelist,  &
        linoz_data_path_in      = linoz_data_path,      &
        linoz_data_type_in      = linoz_data_type,      &
        linoz_data_rmfile_in    = linoz_data_rmfile,    &
        linoz_data_cycle_yr_in  = linoz_data_cycle_yr,  &
        linoz_data_fixed_ymd_in = linoz_data_fixed_ymd, &
        linoz_data_fixed_tod_in = linoz_data_fixed_tod )
   call tracer_cnst_setopts( &
        tracer_cnst_file_in      = tracer_cnst_file,      &
        tracer_cnst_filelist_in  = tracer_cnst_filelist,  &
        tracer_cnst_datapath_in  = tracer_cnst_datapath,  &
        tracer_cnst_type_in      = tracer_cnst_type,      &
        tracer_cnst_specifier_in = tracer_cnst_specifier, &
        tracer_cnst_rmfile_in    = tracer_cnst_rmfile,    &
        tracer_cnst_cycle_yr_in  = tracer_cnst_cycle_yr,  &
        tracer_cnst_fixed_ymd_in = tracer_cnst_fixed_ymd, &
        tracer_cnst_fixed_tod_in = tracer_cnst_fixed_tod )
   call tracer_srcs_setopts( &
        tracer_srcs_file_in      = tracer_srcs_file,      &
        tracer_srcs_filelist_in  = tracer_srcs_filelist,  &
        tracer_srcs_datapath_in  = tracer_srcs_datapath,  &
        tracer_srcs_type_in      = tracer_srcs_type,      &
        tracer_srcs_specifier_in = tracer_srcs_specifier, &
        tracer_srcs_rmfile_in    = tracer_srcs_rmfile,    &
        tracer_srcs_cycle_yr_in  = tracer_srcs_cycle_yr,  &
        tracer_srcs_fixed_ymd_in = tracer_srcs_fixed_ymd, &
        tracer_srcs_fixed_tod_in = tracer_srcs_fixed_tod )

   ! Upper boundary conditions
   call ubc_setopts( &
        snoe_ubc_file_in =snoe_ubc_file, &
        t_pert_ubc_in    =t_pert_ubc, &
        no_xfac_ubc_in   =no_xfac_ubc, &
        tgcm_ubc_file_in =tgcm_ubc_file, &
        tgcm_ubc_data_type_in = tgcm_ubc_data_type, &
        tgcm_ubc_cycle_yr_in = tgcm_ubc_cycle_yr, &
        tgcm_ubc_fixed_ymd_in = tgcm_ubc_fixed_ymd, &
        tgcm_ubc_fixed_tod_in = tgcm_ubc_fixed_tod )

   call aero_model_readnl(nlfile)
   call dust_readnl(nlfile)     
!
   call gas_wetdep_readnl(nlfile)
   call gcr_ionization_readnl(nlfile)
   call epp_ionization_readnl(nlfile)
   call mo_apex_readnl(nlfile)
   call noy_ubc_readnl(nlfile)
   call sulf_readnl(nlfile)
   call species_sums_readnl(nlfile)
   call ocean_emis_readnl(nlfile)

 end subroutine chem_readnl

!================================================================================================

function chem_is_active()
!----------------------------------------------------------------------- 
! Purpose: return true if this package is active
!-----------------------------------------------------------------------
   logical :: chem_is_active
!-----------------------------------------------------------------------
   chem_is_active = is_active
end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
    use chem_mods,       only : gas_pcnst, inv_lst, nfs
    use mo_tracname,     only : solsym

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value
!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer :: m
    
    chem_implements_cnst = .false.
    do m = 1,gas_pcnst
       if( trim(name) /= 'H2O' ) then
          if( trim(name) == solsym(m) ) then
             chem_implements_cnst = .true.
             exit
          end if
       end if
    end do
    do m = 1,nfs
       if( trim(name) /= 'H2O' ) then
          if( trim(name) == inv_lst(m) ) then
             chem_implements_cnst = .true.
             exit
          end if
       endif
    enddo

  end function chem_implements_cnst

  subroutine chem_init(phys_state, pbuf2d)

!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterized greenhouse gas chemistry
!          (declare history variables)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------
    use physics_buffer,      only : physics_buffer_desc, pbuf_get_index
    
    use constituents,        only : cnst_get_ind
    use cam_history,         only : addfld, add_default, horiz_only, fieldname_len
    use chem_mods,           only : gas_pcnst
    use mo_chemini,          only : chemini
    use mo_ghg_chem,         only : ghg_chem_init
    use mo_tracname,         only : solsym
    use llnl_O1D_to_2OH_adj, only : O1D_to_2OH_adj_init
    use lin_strat_chem,      only : lin_strat_chem_inti
    use chlorine_loading_data, only : chlorine_loading_init
    use cfc11star,             only : init_cfc11star
    use phys_control,          only : phys_getopts
    use chem_mods,             only : adv_mass
    use infnan,                only : nan, assignment(=)
    use mo_chem_utls,          only : get_spc_ndx
    use cam_abortutils,        only : endrun
    use aero_model,            only : aero_model_init
    use mo_setsox,             only : sox_inti
    use constituents,          only : sflxnam
    use noy_ubc,             only : noy_ubc_init
    use fire_emissions,      only : fire_emissions_init
    use short_lived_species, only : short_lived_species_initic
    use ocean_emis,          only : ocean_emis_init
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)

    
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer :: m                                ! tracer indicies
    character(len=fieldname_len) :: spc_name
    integer :: n, ii
    logical :: history_aerosol
    logical :: history_chemistry
    logical :: history_cesm_forcing

    character(len=2)  :: unit_basename  ! Units 'kg' or '1' 
    logical :: history_budget                 ! output tendencies and state variables for CAM
                                              ! temperature, water vapor, cloud ice and cloud
                                              ! liquid budgets.
    integer :: history_budget_histfile_num    ! output history file number for budget fields

    call phys_getopts( cam_chempkg_out=chem_name, &
                       history_aerosol_out=history_aerosol , &
                       history_chemistry_out=history_chemistry , &
                       history_budget_out = history_budget , &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       history_cesm_forcing_out = history_cesm_forcing )

    ! aqueous chem initialization
    call sox_inti()

    ! Initialize aerosols
    call aero_model_init( pbuf2d )

!-----------------------------------------------------------------------
! Get liq and ice cloud water indicies
!-----------------------------------------------------------------------
    call cnst_get_ind( 'CLDLIQ', ixcldliq )
    call cnst_get_ind( 'CLDICE', ixcldice )
    call cnst_get_ind( 'NUMLIQ', ixndrop, abort=.false.  )

!-----------------------------------------------------------------------
! get pbuf indicies
!-----------------------------------------------------------------------
    ndx_cld    = pbuf_get_index('CLD')
    ndx_cmfdqr = pbuf_get_index('RPRDTOT')
    ndx_nevapr = pbuf_get_index('NEVAPR')
    ndx_prain  = pbuf_get_index('PRAIN')
    ndx_cldtop = pbuf_get_index('CLDTOP')
    ndx_pblh   = pbuf_get_index('pblh')
    ndx_fsds   = pbuf_get_index('FSDS')

    call addfld( 'HEIGHT',     (/ 'ilev' /),'A','m',       'geopotential height above surface at interfaces (m)' )
    call addfld( 'CT_H2O_GHG', (/ 'lev' /), 'A','kg/kg/s', 'ghg-chem h2o source/sink' )

!-----------------------------------------------------------------------
! Set names of chemistry variable tendencies and declare them as history variables
!-----------------------------------------------------------------------
    do m = 1,gas_pcnst
       spc_name = solsym(m)
       srcnam(m) = 'CT_' // spc_name ! chem tendancy (source/sink)

       call addfld( srcnam(m), (/ 'lev' /), 'A', 'kg/kg/s', trim(spc_name)//' source/sink' )
       call cnst_get_ind(solsym(m), n, abort=.false. ) 
       if ( n > 0 ) then

          if (sflxnam(n)(3:5) == 'num') then  ! name is in the form of "SF****"
             unit_basename = ' 1'
          else
             unit_basename = 'kg'  
          endif

          call addfld (sflxnam(n),horiz_only,    'A',  unit_basename//'/m2/s',trim(solsym(m))//' surface flux')
          if ( history_aerosol .or. history_chemistry ) then 
             call add_default( sflxnam(n), 1, ' ' )
          endif

          if ( history_cesm_forcing ) then
             if ( spc_name == 'NO' .or. spc_name == 'NH3' ) then
                call add_default( sflxnam(n), 1, ' ' )
             endif
          endif
  
       endif
    end do

    ! Add chemical tendency of water vapor to water budget output
    if ( history_budget ) then 
      call add_default ('CT_H2O'  , history_budget_histfile_num, ' ')
    endif

    !-----------------------------------------------------------------------
    ! BAB: 2004-09-01 kludge to define a fixed ubc for water vapor
    !      required because water vapor is not declared by chemistry but
    !      has a fixed ubc only if WACCM chemistry is running.
    !-----------------------------------------------------------------------
    ! this is moved out of chem_register because we need to know where (what pressure) 
    ! the upper boundary is to determine if this is a high top configuration -- after
    ! initialization of ref_pres ...
    if ( 1.e-2_r8 >= ptop_ref .and. ptop_ref > 1.e-5_r8 ) then ! around waccm top, below top of waccmx
       cnst_fixed_ubc(1) = .true.
    else if ( 1.e1_r8 > ptop_ref .and. ptop_ref > 1.e-2_r8 ) then ! well above top of cam and below top of waccm
       call endrun('chem_init: do not know how to set water vapor upper boundary when model top is near mesopause')
    endif

    if ( masterproc ) write(iulog,*) 'chem_init: addfld done'

!-----------------------------------------------------------------------
! Initialize chemistry modules
!-----------------------------------------------------------------------
    call chemini &
       ( euvac_file &
       , photon_file &
       , electron_file &
       , airpl_emis_file &
       , depvel_file &
       , depvel_lnd_file &
       , clim_soilw_file &
       , season_wes_file &
       , xs_coef_file &
       , xs_short_file &
       , xs_long_file &
       , rsf_file &
       , fstrat_file &
       , fstrat_list &
       , srf_emis_specifier &
       , srf_emis_type &
       , srf_emis_cycle_yr &
       , srf_emis_fixed_ymd &
       , srf_emis_fixed_tod &
       , ext_frc_specifier &
       , ext_frc_type &
       , ext_frc_cycle_yr &
       , ext_frc_fixed_ymd &
       , ext_frc_fixed_tod &
       , xactive_prates &
       , exo_coldens_file &
       , tuv_xsect_file &
       , o2_xsect_file &
       , lght_no_prd_factor &
       , pbuf2d &
       )

     if ( ghg_chem ) then
        call ghg_chem_init(phys_state, bndtvg, h2orates)
     endif
     
     call O1D_to_2OH_adj_init()

     call lin_strat_chem_inti(phys_state)
     call chlorine_loading_init( chlorine_loading_file, &
                                 type = chlorine_loading_type, &
                                 ymd = chlorine_loading_fixed_ymd, &
                                 tod = chlorine_loading_fixed_tod )

     call init_cfc11star(pbuf2d)
     
     ! MEGAN emissions initialize
     if (shr_megan_mechcomps_n>0) then

        allocate( megan_indices_map(shr_megan_mechcomps_n) )
        allocate( megan_wght_factors(shr_megan_mechcomps_n) )
        megan_wght_factors(:) = nan

        do n=1,shr_megan_mechcomps_n
           call cnst_get_ind (shr_megan_mechcomps(n)%name,  megan_indices_map(n), abort=.false.)
           ii = get_spc_ndx(shr_megan_mechcomps(n)%name)
           if (ii>0) then
              megan_wght_factors(n) = adv_mass(ii)*1.e-3_r8 ! kg/moles (to convert moles/m2/sec to kg/m2/sec)
           else
              call endrun( 'gas_phase_chemdr_inti: MEGAN compound not in chemistry mechanism : '&
                           //trim(shr_megan_mechcomps(n)%name))
           endif

           ! MEGAN  history fields
           call addfld( 'MEG_'//trim(shr_megan_mechcomps(n)%name),horiz_only,'A','kg/m2/sec',&
                trim(shr_megan_mechcomps(n)%name)//' MEGAN emissions flux')
           if (history_chemistry) then
              call add_default('MEG_'//trim(shr_megan_mechcomps(n)%name), 1, ' ')
           endif

        enddo
     endif
     
     call noy_ubc_init()

     ! Galatic Cosmic Rays ...
     call gcr_ionization_init()

     ! Fire emissions ...
     call fire_emissions_init()

     call short_lived_species_initic()

     call ocean_emis_init()
     
  end subroutine chem_init

!================================================================================
!================================================================================
  subroutine chem_emissions( state, cam_in )
    use aero_model,       only: aero_model_emissions
    use camsrfexch,       only: cam_in_t     
    use constituents,     only: sflxnam
    use cam_history,      only: outfld
    use mo_srf_emissions, only: set_srf_emissions
    use fire_emissions,   only: fire_emissions_srf
    use ocean_emis,       only: ocean_emis_getflux

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    ! local vars

    integer :: lchnk, ncol
    integer :: i, m,n 

    real(r8) :: sflx(pcols,gas_pcnst)
    real(r8) :: megflx(pcols)

    lchnk = state%lchnk
    ncol = state%ncol
    
    ! initialize chemistry constituent surface fluxes to zero
    do m = 2,pcnst
       n = map2chm(m)
       if (n>0) cam_in%cflx(:,m) = 0._r8 
    enddo

    ! aerosol emissions ...
    call aero_model_emissions( state, cam_in )

   ! MEGAN emissions ...
 
    if ( active_Fall_flxvoc .and. shr_megan_mechcomps_n>0 ) then

       ! set MEGAN fluxes 
       do n = 1,shr_megan_mechcomps_n
          do i =1,ncol
             megflx(i) = -cam_in%meganflx(i,n) * megan_wght_factors(n)
             cam_in%cflx(i,megan_indices_map(n)) = cam_in%cflx(i,megan_indices_map(n)) + megflx(i)
          enddo

          ! output MEGAN emis fluxes to history
          call outfld('MEG_'//trim(shr_megan_mechcomps(n)%name), megflx(:ncol), ncol, lchnk)
       enddo

    endif

   ! prescribed emissions from file ...

    !-----------------------------------------------------------------------      
    !        ... Set surface emissions
    !-----------------------------------------------------------------------      
    call set_srf_emissions( lchnk, ncol, sflx(:,:) )

    do m = 1,pcnst
       n = map2chm(m)
       if ( n /= h2o_ndx .and. n > 0 ) then
          cam_in%cflx(:ncol,m) = cam_in%cflx(:ncol,m) + sflx(:ncol,n)
          call outfld( sflxnam(m), cam_in%cflx(:ncol,m), ncol,lchnk )
       endif
    enddo

    ! fire surface emissions if not elevated forcing
    call fire_emissions_srf( lchnk, ncol, cam_in%fireflx, cam_in%cflx )

    ! air-sea exchange of trace gases
    call ocean_emis_getflux(lchnk, ncol, state, cam_in%u10, cam_in%sst, cam_in%ocnfrac, cam_in%icefrac, cam_in%cflx)

  end subroutine chem_emissions

!================================================================================

  subroutine chem_init_cnst( name, latvals, lonvals, mask, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Specify initial mass mixing ratios
! 
!-----------------------------------------------------------------------

    use chem_mods, only : inv_lst

    use physconst,     only : mwdry, mwch4, mwn2o, mwf11, mwf12
    use chem_surfvals, only : chem_surfvals_get

    implicit none

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    character(len=*), intent(in)  :: name       !  constituent name
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    
    real(r8) :: rmwn2o != mwn2o/mwdry ! ratio of mol weight n2o   to dry air
    real(r8) :: rmwch4 != mwch4/mwdry ! ratio of mol weight ch4   to dry air
    real(r8) :: rmwf11 != mwf11/mwdry ! ratio of mol weight cfc11 to dry air
    real(r8) :: rmwf12 != mwf12/mwdry ! ratio of mol weight cfc12 to dry air
    integer  :: ilev, nlev

!-----------------------------------------------------------------------
! initialize local variables
!-----------------------------------------------------------------------

    rmwn2o = mwn2o/mwdry 
    rmwch4 = mwch4/mwdry 
    rmwf11 = mwf11/mwdry 
    rmwf12 = mwf12/mwdry 

!-----------------------------------------------------------------------
! Get initial mixing ratios
!-----------------------------------------------------------------------
    nlev = size(q, 2)
    if ( any( inv_lst .eq. name ) ) then
       do ilev = 1, nlev
          where(mask)
             q(:,ilev) = 0.0_r8
          end where
       end do
    else
       do ilev = 1, nlev
          where(mask)
             q(:,ilev) = 1.e-38_r8
          end where
       end do
    endif

    if ( ghg_chem ) then
       do ilev = 1, nlev
          select case (name)
          case ('N2O')
             where(mask)
                q(:,ilev) = rmwn2o * chem_surfvals_get('N2OVMR')
             end where
          case ('CH4')
             where(mask)
                q(:,ilev) = rmwch4 * chem_surfvals_get('CH4VMR')
             end where
          case ('CFC11')
             where(mask)
                q(:,ilev) = rmwf11 * chem_surfvals_get('F11VMR')
             end where
          case ('CFC12')
             where(mask)
                q(:,ilev) = rmwf12 * chem_surfvals_get('F12VMR')
             end where
          end select
       end do
    end if

  end subroutine chem_init_cnst

  subroutine chem_timestep_init(phys_state,pbuf2d)

    use time_manager,      only : get_nstep
    use time_manager,      only : get_curr_calday
    use mo_srf_emissions,  only : set_srf_emissions_time
    use mo_sulf,           only : set_sulf_time
    use mo_extfrc,         only : extfrc_timestep_init
    use mo_flbc,           only : flbc_chk
    use tracer_cnst,       only : tracer_cnst_adv
    use tracer_srcs,       only : tracer_srcs_adv
    use mo_ghg_chem,       only : ghg_chem_timestep_init

    use mo_aurora,         only : aurora_timestep_init
    use mo_photo,          only : photo_timestep_init
    use linoz_data,        only : linoz_data_adv
    use chlorine_loading_data, only : chlorine_loading_advance
    use noy_ubc,           only : noy_ubc_advance

    use cfc11star,         only : update_cfc11star
    use physics_buffer,    only : physics_buffer_desc
    use ocean_emis,        only : ocean_emis_advance

    implicit none

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    real(r8) :: calday
    integer :: nstep

    nstep = get_nstep()
    chem_step = mod( nstep, chem_freq ) == 0

    if ( .not. chem_step ) return

    !-----------------------------------------------------------------------
    ! get current calendar day of year
    !-----------------------------------------------------------------------
    calday = get_curr_calday( )

    !-----------------------------------------------------------------------
    ! Set emissions timing factors
    !-----------------------------------------------------------------------
    call set_srf_emissions_time( pbuf2d, phys_state )

    !-----------------------------------------------------------------------
    ! Set external forcings timing factors
    !-----------------------------------------------------------------------
    call extfrc_timestep_init( pbuf2d, phys_state )

    !-----------------------------------------------------------------------
    ! Set sulf timing factors
    !-----------------------------------------------------------------------
    call set_sulf_time( pbuf2d, phys_state  )

    !-----------------------------------------------------------------------
    ! Set fixed lower boundary timing factors
    !-----------------------------------------------------------------------
    call flbc_chk

    !-----------------------------------------------------------------------
    ! NOy upper boundary conditions for low top model
    !-----------------------------------------------------------------------
    call noy_ubc_advance(pbuf2d, phys_state)

    !-----------------------------------------------------------------------
    ! Set fixed offline tracers
    !-----------------------------------------------------------------------
    call tracer_cnst_adv(pbuf2d, phys_state)

    !-----------------------------------------------------------------------
    ! Set fixed offline tracer sources
    !-----------------------------------------------------------------------
    call tracer_srcs_adv(pbuf2d, phys_state)

    !-----------------------------------------------------------------------
    ! Advance the linoz data
    !-----------------------------------------------------------------------
    call linoz_data_adv(pbuf2d, phys_state)
    call chlorine_loading_advance()

    if ( ghg_chem ) then
       call ghg_chem_timestep_init(phys_state)
    endif
    
    !-----------------------------------------------------------------------
    ! Set up aurora
    !-----------------------------------------------------------------------
    call aurora_timestep_init

    !-----------------------------------------------------------------------------
    !   ... setup the time interpolation for mo_photo
    !-----------------------------------------------------------------------------
    call photo_timestep_init( calday )

    call update_cfc11star( pbuf2d, phys_state )
    
    ! Galatic Cosmic Rays ...
    call gcr_ionization_adv( pbuf2d, phys_state )
    call epp_ionization_adv()

    call ocean_emis_advance( pbuf2d, phys_state )

  end subroutine chem_timestep_init

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf,  fh2o)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to parameterized greenhouse gas chemisty (source/sink).
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------

    use physics_buffer,      only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use cam_history,         only : outfld
    use time_manager,        only : get_curr_calday
    use mo_gas_phase_chemdr, only : gas_phase_chemdr
    use camsrfexch,          only : cam_in_t, cam_out_t     
    use perf_mod,            only : t_startf, t_stopf
    use tropopause,          only : tropopause_findChemTrop, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE, tropopause_find
    use mo_drydep,           only : drydep_update
    use mo_neu_wetdep,       only : neu_wetdep_tend
    use aerodep_flx,         only : aerodep_flx_prescribed
    use short_lived_species, only : short_lived_species_writeic
    
    implicit none

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    real(r8),            intent(in)    :: dt              ! time step
    type(physics_state), intent(in)    :: state           ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend           ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(inout) :: cam_out
    real(r8),            intent(out)   :: fh2o(pcols)     ! h2o flux to balance source from chemistry
    

    type(physics_buffer_desc), pointer :: pbuf(:)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer  :: i, k, m, n                         ! indicies
    integer  :: lchnk                              ! chunk identifier
    integer  :: ncol                               ! number of atmospheric columns
    real(r8) :: calday                             ! current calendar day of year
    real(r8) :: cldw(pcols,pver)                   ! cloud water (kg/kg)
    real(r8) :: chem_dt              ! time step
    real(r8) :: drydepflx(pcols,pcnst)             ! dry deposition fluxes (kg/m2/s)
    real(r8) :: wetdepflx(pcols,pcnst)             ! wet deposition fluxes (kg/m2/s)
    integer  :: tropLev(pcols), tropLevChem(pcols)
    real(r8) :: ncldwtr(pcols,pver)                ! droplet number concentration (#/kg)
    real(r8), pointer :: fsds(:)     ! longwave down at sfc
    real(r8), pointer :: pblh(:)
    real(r8), pointer :: prain(:,:)
    real(r8), pointer :: cldfr(:,:)
    real(r8), pointer :: cmfdqr(:,:)
    real(r8), pointer :: nevapr(:,:)
    real(r8), pointer :: cldtop(:)
    real(r8) :: nhx_nitrogen_flx(pcols)
    real(r8) :: noy_nitrogen_flx(pcols)

    integer :: tim_ndx

    logical :: lq(pcnst)

    if ( .not. chem_step ) return

    chem_dt = chem_freq*dt

    lchnk = state%lchnk
    ncol  = state%ncol

    call short_lived_species_writeic( lchnk, pbuf )

    lq(:) = .false.
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          lq(n) = .true.
       end if
    end do
    if ( ghg_chem ) lq(1) = .true.

    call physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)
    
    call drydep_update( state, cam_in )

!-----------------------------------------------------------------------
! get current calendar day of year
!-----------------------------------------------------------------------
    calday = get_curr_calday()

!-----------------------------------------------------------------------
! get tropopause level
!-----------------------------------------------------------------------
    if (.not.chem_use_chemtrop) then
       call tropopause_find(state,tropLev)
       tropLevChem=tropLev
    else
       call tropopause_find(state,tropLev)
       call tropopause_findChemTrop(state, tropLevChem)
    endif

    tim_ndx = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ndx_fsds,       fsds)
    call pbuf_get_field(pbuf, ndx_pblh,       pblh)
    call pbuf_get_field(pbuf, ndx_prain,      prain,  start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cld,        cldfr,  start=(/1,1,tim_ndx/), kount=(/ncol,pver,1/) )
    call pbuf_get_field(pbuf, ndx_cmfdqr,     cmfdqr, start=(/1,1/),         kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_nevapr,     nevapr, start=(/1,1/),         kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldtop,     cldtop )

!-----------------------------------------------------------------------
! call Neu wet dep scheme
!-----------------------------------------------------------------------
    call neu_wetdep_tend(lchnk,ncol,state%q,state%pmid,state%pdel,state%zi,state%t,dt, &
         prain, nevapr, cldfr, cmfdqr, ptend%q, wetdepflx)

!-----------------------------------------------------------------------
! compute tendencies and surface fluxes
!-----------------------------------------------------------------------
    call t_startf( 'chemdr' )
    do k = 1,pver
       cldw(:ncol,k) = state%q(:ncol,k,ixcldliq) + state%q(:ncol,k,ixcldice)
       if (ixndrop>0) &
            ncldwtr(:ncol,k) = state%q(:ncol,k,ixndrop)
    end do

    call gas_phase_chemdr(lchnk, ncol, imozart, state%q, &
                          state%phis, state%zm, state%zi, calday, &
                          state%t, state%pmid, state%pdel, state%pint, &
                          cldw, tropLev, tropLevChem, ncldwtr, state%u, state%v, &
                          chem_dt, state%ps, xactive_prates, &
                          fsds, cam_in%ts, cam_in%asdir, cam_in%ocnfrac, cam_in%icefrac, &
                          cam_out%precc, cam_out%precl, cam_in%snowhland, ghg_chem, state%latmapback, &
                          drydepflx, wetdepflx, cam_in%cflx, cam_in%fireflx, cam_in%fireztop, &
                          nhx_nitrogen_flx, noy_nitrogen_flx, ptend%q, pbuf )
    if (associated(cam_out%nhx_nitrogen_flx)) then
       cam_out%nhx_nitrogen_flx(:ncol) = nhx_nitrogen_flx(:ncol)
    endif
    if (associated(cam_out%noy_nitrogen_flx)) then
       cam_out%noy_nitrogen_flx(:ncol) = noy_nitrogen_flx(:ncol)
    endif

    call t_stopf( 'chemdr' )

!-----------------------------------------------------------------------
! set flags for tracer tendencies (water and gas phase constituents)
! record tendencies on history files
!-----------------------------------------------------------------------
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          call outfld( srcnam(m), ptend%q(:,:,n), pcols, lchnk )
       end if

       ! if the user has specified prescribed aerosol dep fluxes then 
       ! do not set cam_out dep fluxes according to the prognostic aerosols
       if (.not.aerodep_flx_prescribed()) then
          ! set deposition fluxes in the export state
          select case (trim(cnst_name(n)))
          case('CB1')
             do i = 1, ncol
                cam_out%bcphodry(i) = max(drydepflx(i,n), 0._r8)
             end do
          case('CB2')
             do i = 1, ncol
                cam_out%bcphidry(i) = max(drydepflx(i,n), 0._r8)
             end do
          case('OC1')
             do i = 1, ncol
                cam_out%ocphodry(i) = max(drydepflx(i,n), 0._r8)
             end do
          case('OC2')
             do i = 1, ncol
                cam_out%ocphidry(i) = max(drydepflx(i,n), 0._r8)
             end do
          end select
       endif
    end do
    if ( ghg_chem ) then
       ptend%lq(1) = .true.
       call outfld( 'CT_H2O_GHG', ptend%q(:,:,1), pcols, lchnk )
    endif

    call outfld( 'HEIGHT', state%zi(:ncol,:),  ncol, lchnk )

!-----------------------------------------------------------------------
!  turn off water vapor tendency if radiatively passive
!-----------------------------------------------------------------------
    if (chem_rad_passive) then
       ptend%lq(1) = .false.
       ptend%q(:ncol,:,1) = 0._r8
    endif

!-----------------------------------------------------------------------
! Compute water vapor flux required to make conservation check
!-----------------------------------------------------------------------
    fh2o(:ncol) = 0._r8
    do k = 1,pver
       fh2o(:ncol) = fh2o(:ncol) + ptend%q(:ncol,k,1)*state%pdel(:ncol,k)/gravit
    end do
    
  end subroutine chem_timestep_tend

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine chem_final
  end subroutine chem_final

!-------------------------------------------------------------------
!-------------------------------------------------------------------

  subroutine chem_init_restart( File )
    use pio, only : file_desc_t
    use tracer_cnst,      only: init_tracer_cnst_restart
    use tracer_srcs,      only: init_tracer_srcs_restart
    use linoz_data, only : init_linoz_data_restart
    implicit none
    type(file_desc_t),intent(inout) :: File     ! pio File pointer

    !
    ! data for offline tracers
    !
    call init_tracer_cnst_restart(File)
    call init_tracer_srcs_restart(File)
    call init_linoz_data_restart(File)
  end subroutine chem_init_restart
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine chem_write_restart( File )
    use tracer_cnst, only: write_tracer_cnst_restart
    use tracer_srcs, only: write_tracer_srcs_restart
    use linoz_data,  only: write_linoz_data_restart
    use pio, only : file_desc_t
    implicit none
    type(file_desc_t) :: File

    !
    ! data for offline tracers
    !
    call write_tracer_cnst_restart(File)
    call write_tracer_srcs_restart(File)
    call write_linoz_data_restart(File)
  end subroutine chem_write_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine chem_read_restart( File )
    use tracer_cnst, only: read_tracer_cnst_restart
    use tracer_srcs, only: read_tracer_srcs_restart
    use linoz_data,  only: read_linoz_data_restart

    use pio, only : file_desc_t
    implicit none
    type(file_desc_t) :: File

    !
    ! data for offline tracers
    !
    call read_tracer_cnst_restart(File)
    call read_tracer_srcs_restart(File)
    call read_linoz_data_restart(File)
  end subroutine chem_read_restart

end module chemistry
