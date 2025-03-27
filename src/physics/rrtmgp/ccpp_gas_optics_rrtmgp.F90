module ccpp_gas_optics_rrtmgp
   ! CCPP wrapper for ty_gas_optics_rrtmgp DDT from RRTMGP
   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

   !> \section arg_table_ty_gas_optics_rrtmgp_ccpp Argument Table
   !! \htmlinclude ty_gas_optics_rrtmgp_ccpp.html
   type, public, extends(ty_gas_optics_rrtmgp) :: ty_gas_optics_rrtmgp_ccpp
   end type

end module ccpp_gas_optics_rrtmgp
