!> \file rrtmgp_lw_gas_optics_data.F90
!!

!> This module contains an init routine to initialize the gas optics object
!>  with data read in from file on the host side
module rrtmgp_lw_gas_optics_data
  use machine,                 only: kind_phys
  use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
  use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
  use radiation_tools,         only: check_error_msg

  implicit none


contains
!> \section arg_table_rrtmgp_lw_gas_optics_data_init Argument Table
!! \htmlinclude rrtmgp_lw_gas_optics_data_init.html
!!
  subroutine rrtmgp_lw_gas_optics_data_init(kdist, available_gases, gas_names,    &
                  key_species, band2gpt, band_lims_wavenum, press_ref, press_ref_trop, &
                  temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower,     &
                  kminor_upper, gas_minor, identifier_minor, minor_gases_lower,        &
                  minor_gases_upper, minor_limits_gpt_lower, minor_limits_gpt_upper,   &
                  minor_scales_with_density_lower, minor_scales_with_density_upper,    &
                  scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower,     &
                  scale_by_complement_upper, kminor_start_lower, kminor_start_upper,   &
                  totplnk, planck_frac, rayl_lower, rayl_upper, optimal_angle_fit,     &
                  errmsg, errflg)

    ! Inputs
    class(ty_gas_concs_ccpp),            intent(in) :: available_gases                    ! Gas concentrations object
    character(len=*),      dimension(:), intent(in) :: gas_names                          ! Names of absorbing gases
    character(len=*),      dimension(:), intent(in) :: gas_minor                          ! Name of absorbing minor gas
    character(len=*),      dimension(:), intent(in) :: identifier_minor                   ! Unique string identifying minor gas
    character(len=*),      dimension(:), intent(in) :: minor_gases_lower                  ! Names of minor absorbing gases in lower atmosphere
    character(len=*),      dimension(:), intent(in) :: minor_gases_upper                  ! Names of minor absorbing gases in upper atmosphere
    character(len=*),      dimension(:), intent(in) :: scaling_gas_lower                  ! Absorption also depends on the concentration of this gas in the lower atmosphere
    character(len=*),      dimension(:), intent(in) :: scaling_gas_upper                  ! Absorption also depends on the concentration of this gas in the upper atmosphere
    integer,           dimension(:,:,:), intent(in) :: key_species                        ! Key species pair for each band
    integer,             dimension(:,:), intent(in) :: band2gpt                           ! Array for converting shortwave band limits to g-points
    integer,             dimension(:,:), intent(in) :: minor_limits_gpt_lower             ! Beginning and ending gpoint for each minor interval in lower atmosphere
    integer,             dimension(:,:), intent(in) :: minor_limits_gpt_upper             ! Beginning and ending gpoint for each minor interval in upper atmosphere
    integer,               dimension(:), intent(in) :: kminor_start_lower                 ! Starting index in the [1,nContributors] vector for a contributor given by "minor_gases_lower"
    integer,               dimension(:), intent(in) :: kminor_start_upper                 ! Starting index in the [1,nContributors] vector for a contributor given by "minor_gases_upper"
    logical,               dimension(:), intent(in) :: minor_scales_with_density_lower    ! Density scaling is applied to minor absorption coefficients in the lower atmosphere
    logical,               dimension(:), intent(in) :: scale_by_complement_lower          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T) in the lower atmosphere
    logical,               dimension(:), intent(in) :: minor_scales_with_density_upper    ! Density scaling is applied to minor absorption coefficients in the upper atmosphere
    logical,               dimension(:), intent(in) :: scale_by_complement_upper          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T) in the upper atmosphere
    real(kind_phys), dimension(:,:,:,:), intent(in) :: kmajor                             ! Stored absorption coefficients due to major absorbing gases
    real(kind_phys), dimension(:,:,:,:), intent(in) :: planck_frac                        ! Fraction of band-integrated Planck energy associated with each g-point
    real(kind_phys),   dimension(:,:,:), intent(in) :: kminor_lower                       ! Transformed from [nTemp x nEta x nGpt x nAbsorbers] array to [nTemp x nEta x nContributors] array
    real(kind_phys),   dimension(:,:,:), intent(in) :: kminor_upper                       ! Transformed from [nTemp x nEta x nGpt x nAbsorbers] array to [nTemp x nEta x nContributors] array
    real(kind_phys),   dimension(:,:,:), intent(in) :: vmr_ref                            ! Volume mixing ratios for reference atmosphere
    real(kind_phys),     dimension(:,:), intent(in) :: band_lims_wavenum                  ! Beginning and ending wavenumber for each band [cm-1]
    real(kind_phys),     dimension(:,:), intent(in) :: totplnk                            ! Integrated Planck function by band
    real(kind_phys),     dimension(:,:), intent(in) :: optimal_angle_fit                  ! Coefficients for linear fit used in longwave optimal angle RT calculation
    real(kind_phys),       dimension(:), intent(in) :: press_ref                          ! Pressures for reference atmosphere [Pa]
    real(kind_phys),       dimension(:), intent(in) :: temp_ref                           ! Temperatures for reference atmosphere [K]
    real(kind_phys),                     intent(in) :: press_ref_trop                     ! Reference pressure separating the lower and upper atmosphere [Pa]
    real(kind_phys),                     intent(in) :: temp_ref_p                         ! Standard spectroscopic reference pressure [Pa]
    real(kind_phys),                     intent(in) :: temp_ref_t                         ! Standard spectroscopic reference temperature [K]
    real(kind_phys), dimension(:,:,:), allocatable, intent(in) :: rayl_lower              ! Stored coefficients due to rayleigh scattering contribution in lower part of atmosphere
    real(kind_phys), dimension(:,:,:), allocatable, intent(in) :: rayl_upper              ! Stored coefficients due to rayleigh scattering contribution in upper part of atmosphere
 
    ! Outputs
    class(ty_gas_optics_rrtmgp_ccpp),  intent(inout) :: kdist                             ! RRTMGP gas optics object
    character(len=*),                    intent(out) :: errmsg                            ! CCPP error message
    integer,                             intent(out) :: errflg                            ! CCPP error code

    ! Initialize error variables
    errmsg = ''
    errflg = 0

    ! Initialize the gas optics object with data.
    errmsg = kdist%gas_props%load( &
         available_gases%gas_concs, gas_names, key_species,    &
         band2gpt, band_lims_wavenum,                          &
         press_ref, press_ref_trop, temp_ref,                  &
         temp_ref_p, temp_ref_t, vmr_ref,                      &
         kmajor, kminor_lower, kminor_upper,                   &
         gas_minor, identifier_minor,                          &
         minor_gases_lower, minor_gases_upper,                 &
         minor_limits_gpt_lower, minor_limits_gpt_upper,       &
         minor_scales_with_density_lower,                      &
         minor_scales_with_density_upper,                      &
         scaling_gas_lower, scaling_gas_upper,                 &
         scale_by_complement_lower, scale_by_complement_upper, &
         kminor_start_lower, kminor_start_upper,               &
         totplnk, planck_frac, rayl_lower, rayl_upper,         &
         optimal_angle_fit)

    if (len_trim(errmsg) > 0) then
       errflg = 1
    end if
    call check_error_msg('rrtmgp_lw_gas_optics_init_load', errmsg)

  end subroutine rrtmgp_lw_gas_optics_data_init

end module rrtmgp_lw_gas_optics_data
