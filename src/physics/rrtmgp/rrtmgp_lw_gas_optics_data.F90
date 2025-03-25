!> \file rrtmgp_lw_gas_optics_data.F90
!!

!> This module contains an init routine to initialize the gas optics object
!>  with data read in from file on the host side
module rrtmgp_lw_gas_optics_data
  use machine,               only: kind_phys
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs  
!  use radiation_tools,       only: check_error_msg

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
    class(ty_gas_concs),                 intent(in) :: available_gases
    character(len=*),      dimension(:), intent(in) :: gas_names
    character(len=*),      dimension(:), intent(in) :: gas_minor
    character(len=*),      dimension(:), intent(in) :: identifier_minor
    character(len=*),      dimension(:), intent(in) :: minor_gases_lower
    character(len=*),      dimension(:), intent(in) :: minor_gases_upper
    character(len=*),      dimension(:), intent(in) :: scaling_gas_lower
    character(len=*),      dimension(:), intent(in) :: scaling_gas_upper
    integer,           dimension(:,:,:), intent(in) :: key_species
    integer,             dimension(:,:), intent(in) :: band2gpt
    integer,             dimension(:,:), intent(in) :: minor_limits_gpt_lower
    integer,             dimension(:,:), intent(in) :: minor_limits_gpt_upper
    integer,               dimension(:), intent(in) :: kminor_start_lower
    integer,               dimension(:), intent(in) :: kminor_start_upper
    logical,               dimension(:), intent(in) :: minor_scales_with_density_lower
    logical,               dimension(:), intent(in) :: scale_by_complement_lower
    logical,               dimension(:), intent(in) :: minor_scales_with_density_upper
    logical,               dimension(:), intent(in) :: scale_by_complement_upper
    real(kind_phys), dimension(:,:,:,:), intent(in) :: kmajor
    real(kind_phys), dimension(:,:,:,:), intent(in) :: planck_frac
    real(kind_phys),   dimension(:,:,:), intent(in) :: kminor_lower
    real(kind_phys),   dimension(:,:,:), intent(in) :: kminor_upper
    real(kind_phys),   dimension(:,:,:), intent(in) :: vmr_ref
    real(kind_phys),   dimension(:,:,:), allocatable, intent(in) :: rayl_lower
    real(kind_phys),   dimension(:,:,:), allocatable, intent(in) :: rayl_upper
    real(kind_phys),     dimension(:,:), intent(in) :: band_lims_wavenum
    real(kind_phys),     dimension(:,:), intent(in) :: totplnk
    real(kind_phys),     dimension(:,:), intent(in) :: optimal_angle_fit
    real(kind_phys),       dimension(:), intent(in) :: press_ref
    real(kind_phys),       dimension(:), intent(in) :: temp_ref
    real(kind_phys),                     intent(in) :: press_ref_trop
    real(kind_phys),                     intent(in) :: temp_ref_p
    real(kind_phys),                     intent(in) :: temp_ref_t
 
    ! Outputs
    class(ty_gas_optics_rrtmgp), intent(inout) :: kdist  !< RRTMGP gas optics object
    character(len=*), intent(out) :: errmsg              !< CCPP error message
    integer,          intent(out) :: errflg              !< CCPP error code

    ! Initialize error variables
    errmsg = ''
    errflg = 0

    ! Initialize the gas optics object with data.
    errmsg = kdist%load( &
         available_gases, gas_names, key_species,              &
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

  end subroutine rrtmgp_lw_gas_optics_data_init

end module rrtmgp_lw_gas_optics_data
