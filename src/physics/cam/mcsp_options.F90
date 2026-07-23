! Reads MCSP namelist options without run phase
module mcsp_options

 use ccpp_kinds, only: kind_phys
 implicit none
 private


 public :: mcsp_options_init

contains
!> \section arg_table_mcsp_options_init Argument Table
!! \htmlinclude mcsp_conv_options_init.html
 subroutine mcsp_options_init(masterproc, iulog, MCSP_heat_coeff, MCSP_moisture_coeff, MCSP_uwind_coeff, MCSP_vwind_coeff, &
                 MCSP_storm_speed_pref, MCSP_conv_depth_min, mcsp_shear_min)
       

  logical, intent(in) :: masterproc
  integer, intent(in) :: iulog

  real(kind_phys), intent(in) :: MCSP_heat_coeff
  real(kind_phys), intent(in) :: MCSP_moisture_coeff
  real(kind_phys), intent(in) :: MCSP_uwind_coeff
  real(kind_phys), intent(in) :: MCSP_vwind_coeff
  real(kind_phys), intent(in) :: MCSP_storm_speed_pref
  real(kind_phys), intent(in) :: MCSP_conv_depth_min
  real(kind_phys), intent(in) :: mcsp_shear_min

  if ( masterproc ) then
    write(iulog,*) 'tuning parameters mcsp_options_init: MCSP_heat_coeff', MCSP_heat_coeff
    write(iulog,*) 'tuning parameters mcsp_options_init: MCSP_moisture_coeff', MCSP_moisture_coeff
    write(iulog,*) 'tuning parameters mcsp_options_init: MCSP_uwind_coeff', MCSP_uwind_coeff
    write(iulog,*) 'tuning parameters mcsp_options_init: MCSP_vwind_coeff', MCSP_vwind_coeff
    write(iulog,*) 'tuning parameters mcsp_options_init: MCSP_storm_speed_pref', MCSP_storm_speed_pref
    write(iulog,*) 'tuning parameters mcsp_options_init: MCSP_conv_depth_min', MCSP_conv_depth_min
    write(iulog,*) 'tuning parameters mcsp_options_init: mcsp_shear_min', mcsp_shear_min
  end if


 end subroutine mcsp_options_init

end module mcsp_options
