module cam_thermo_formula

   implicit none
   private
   save

   ! energy_formula options for use by CCPPized check_energy
   integer, public, parameter :: ENERGY_FORMULA_DYCORE_FV   = 0  ! vc_moist_pressure
   integer, public, parameter :: ENERGY_FORMULA_DYCORE_SE   = 1  ! vc_dry_pressure
   integer, public, parameter :: ENERGY_FORMULA_DYCORE_MPAS = 2  ! vc_height

   !REMOVECAM: in CAM, energy_formula_physics and energy_formula_dycore still uses vc_physics
   ! and vc_dycore in dyn_tests_utils. The values are the same.
end module cam_thermo_formula
