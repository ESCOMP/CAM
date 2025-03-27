module ccpp_optical_props
   ! CCPP wrapper for ty_optical_props_* DDTs from RRTMGP
   use mo_optical_props, only: ty_optical_props_1scl
   use mo_optical_props, only: ty_optical_props_2str
   use mo_optical_props, only: ty_optical_props_arry

   !> \section arg_table_ty_optical_props_1scl_ccpp Argument Table
   !! \htmlinclude ty_optical_props_1scl_ccpp.html
   type, public, extends(ty_optical_props_1scl) :: ty_optical_props_1scl_ccpp
   end type

   !> \section arg_table_ty_optical_props_2str_ccpp Argument Table
   !! \htmlinclude ty_optical_props_2str_ccpp.html
   type, public, extends(ty_optical_props_2str) :: ty_optical_props_2str_ccpp
   end type

   !> \section arg_table_ty_optical_props_arry_ccpp Argument Table
   !! \htmlinclude ty_optical_props_arry_ccpp.html
   type, public, abstract, extends(ty_optical_props_arry) :: ty_optical_props_arry_ccpp
   end type

end module ccpp_optical_props
