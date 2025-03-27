module ccpp_gas_concentrations
   ! CCPP wrapper for ty_gas_concs DDT from RRTMGP
   use mo_gas_concentrations, only: ty_gas_concs

   !> \section arg_table_ty_gas_concs_ccpp Argument Table
   !! \htmlinclude ty_gas_concs_ccpp.html
   type, public :: ty_gas_concs_ccpp
      type(ty_gas_concs) :: gas_concs
   end type

end module ccpp_gas_concentrations
