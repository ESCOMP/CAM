module ccpp_source_functions
   ! CCPP wrapper for ty_source_func_lw DDT from RRTMGP
   use mo_source_functions, only: ty_source_func_lw

   !> \section arg_table_ty_source_func_lw_ccpp Argument Table
   !! \htmlinclude ty_source_func_lw_ccpp.html
   type, public :: ty_source_func_lw_ccpp
      type(ty_source_func_lw) :: sources
   end type

end module ccpp_source_functions
