
module mo_chemini

  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmd_utils,   only : masterproc
  use cam_logfile,  only : iulog

  implicit none

  private
  public :: chemini

contains

  subroutine chemini &
       ( euvac_file &
       , photon_file &
       , electron_file &
       , airpl_emis_file &
       , depvel_lnd_file &
       , xs_coef_file &
       , xs_short_file &
       , xs_long_file &
       , photo_max_zen &
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
       , exo_coldens_file &
       , use_hemco &
       , pbuf2d &
       )

    !-----------------------------------------------------------------------
    ! 	... Chemistry module intialization
    !-----------------------------------------------------------------------

    use mo_airplane,       only : airpl_src
    use mo_srf_emissions,  only : srf_emissions_inti
    use mo_sulf,           only : sulf_inti
    use mo_photo,          only : photo_inti
    use mo_tuvx,           only : tuvx_init, tuvx_active
    use mo_drydep,         only : drydep_inti
    use mo_imp_sol,        only : imp_slv_inti
    use mo_exp_sol,        only : exp_sol_inti
    use spmd_utils,        only : iam
    use mo_fstrat,         only : fstrat_inti
    use mo_sethet,         only : sethet_inti
    use mo_usrrxt,         only : usrrxt_inti
    use hco_cc_emissions,  only : hco_extfrc_inti
    use mo_extfrc,         only : extfrc_inti
    use mo_setext,         only : setext_inti
    use mo_setinv,         only : setinv_inti
    use mo_gas_phase_chemdr,only: gas_phase_chemdr_inti

    use tracer_cnst,       only : tracer_cnst_init
    use tracer_srcs,       only : tracer_srcs_init
    use mo_airglow,        only : init_airglow
    use mo_mean_mass,      only : init_mean_mass
    use mo_mass_xforms,    only : init_mass_xforms
    use mo_strato_rates,   only : init_strato_rates
    use mo_cph,            only : init_cph
    use mo_sad,            only : sad_inti
    use euvac,             only : euvac_init
    use mo_heatnirco2,     only : heatnirco2_init
    use mo_waccm_hrates,   only : init_hrates
    use mo_aurora,         only : aurora_inti
    use clybry_fam,        only : clybry_fam_init
    use mo_neu_wetdep,     only : neu_wetdep_init
    use physics_buffer,    only : physics_buffer_desc
    use cam_abortutils,    only : endrun

    character(len=*), intent(in) :: euvac_file
    character(len=*), intent(in) :: photon_file
    character(len=*), intent(in) :: electron_file

    character(len=*), intent(in) :: airpl_emis_file
    character(len=*), intent(in) :: depvel_lnd_file
    character(len=*), intent(in) :: xs_coef_file
    character(len=*), intent(in) :: xs_short_file
    character(len=*), intent(in) :: xs_long_file
    real(r8),         intent(in) :: photo_max_zen
    character(len=*), intent(in) :: rsf_file
    character(len=*), intent(in) :: fstrat_file
    character(len=*), intent(in) :: fstrat_list(:)
    character(len=*), dimension(:), intent(in) :: srf_emis_specifier
    character(len=*), dimension(:), intent(in) :: ext_frc_specifier
    character(len=*), intent(in) :: exo_coldens_file
    character(len=*), intent(in) :: ext_frc_type
    integer,          intent(in) :: ext_frc_cycle_yr
    integer,          intent(in) :: ext_frc_fixed_ymd
    integer,          intent(in) :: ext_frc_fixed_tod
    character(len=*), intent(in) :: srf_emis_type
    integer,          intent(in) :: srf_emis_cycle_yr
    integer,          intent(in) :: srf_emis_fixed_ymd
    integer,          intent(in) :: srf_emis_fixed_tod
    logical,          intent(in) :: use_hemco

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !	... initialize the implicit solver
    !-----------------------------------------------------------------------
    call imp_slv_inti()
    call exp_sol_inti()

    call gas_phase_chemdr_inti()

    call init_mean_mass
    call init_mass_xforms

    call setinv_inti()
    call sethet_inti()
    call usrrxt_inti()
    call init_hrates
    call init_airglow

    call init_strato_rates
    call init_cph

    !-----------------------------------------------------------------------
    ! 	... initialize tracer modules
    !-----------------------------------------------------------------------
    call tracer_cnst_init()
    call tracer_srcs_init()

    !-----------------------------------------------------------------------
    ! 	... read time-independent airplane emissions
    !-----------------------------------------------------------------------
    call airpl_src(airpl_emis_file)
    if (masterproc) write(iulog,*) 'chemini: after airpl_src on node ',iam

    !-----------------------------------------------------------------------
    ! 	... read time-dependent surface flux dataset
    !-----------------------------------------------------------------------
    call srf_emissions_inti ( srf_emis_specifier, srf_emis_type, srf_emis_cycle_yr, srf_emis_fixed_ymd, srf_emis_fixed_tod)

    if (masterproc) write(iulog,*) 'chemini: after srf_emissions_inti on node ',iam

    !-----------------------------------------------------------------------
    ! 	... initialize external forcings module
    !-----------------------------------------------------------------------
    call setext_inti()

    if ( use_hemco ) then
        ! Initialize HEMCO version of extfrc_inti
        call hco_extfrc_inti()
        if (masterproc) write(iulog,*) 'chemini: after hco_extfrc_inti on node ',iam
    else
        call extfrc_inti(ext_frc_specifier, ext_frc_type, ext_frc_cycle_yr, ext_frc_fixed_ymd, ext_frc_fixed_tod)
        if (masterproc) write(iulog,*) 'chemini: after extfrc_inti on node ',iam
    endif

    call sulf_inti()
    if (masterproc) write(iulog,*) 'chemini: after sulf_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... initialize the sad module
    !-----------------------------------------------------------------------
    call sad_inti(pbuf2d)
    if (masterproc) write(iulog,*) 'chemini: after sad_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... initialize the dry deposition module
    !-----------------------------------------------------------------------
    call drydep_inti(depvel_lnd_file)

    if (masterproc) write(iulog,*) 'chemini: after drydep_inti on node ',iam

    !-----------------------------------------------------------------------
    !	... Initialize the upper boundary module
    !-----------------------------------------------------------------------
    call fstrat_inti( fstrat_file, fstrat_list )
    if (masterproc) write(iulog,*) 'chemini: after fstrat_inti on node ',iam

    !-----------------------------------------------------------------------
    ! 	... initialize the co2 nir heating module
    !-----------------------------------------------------------------------
    call heatnirco2_init

    !-----------------------------------------------------------------------
    ! 	... initialize photorate module
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! 	... initialize the euvac etf module
    !-----------------------------------------------------------------------
    call euvac_init (euvac_file)

    call photo_inti( xs_coef_file, xs_short_file, xs_long_file, rsf_file, &
           photon_file, electron_file, exo_coldens_file, photo_max_zen )
    if (masterproc) write(iulog,*) 'chemini: after photo_inti on node ',iam

    !-----------------------------------------------------------------------
    ! 	... initialize the TUV-x photolysis rate constant calculator
    !-----------------------------------------------------------------------
    if( tuvx_active ) then
      call tuvx_init( photon_file, electron_file, photo_max_zen, pbuf2d )
      if (masterproc) write(iulog,*) 'chemini: after tuvx_init on node ',iam
    end if

    !-----------------------------------------------------------------------
    !	... initialize ion production
    !-----------------------------------------------------------------------
    call aurora_inti(pbuf2d)
    if (masterproc) write(iulog,*) 'chemini: after aurora_inti'

    call neu_wetdep_init()
    if (masterproc) write(iulog,*) 'chemini: after wetdep_init'

    call clybry_fam_init()

    if (masterproc) write(iulog,*) 'chemini: finished on node ',iam

  end subroutine chemini

end module mo_chemini
