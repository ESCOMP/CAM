module ccpp_fluxes
   ! CCPP wrapper for ty_fluxes DDT from RRTMGP
   use mo_fluxes, only: ty_fluxes
   use mo_fluxes, only: ty_fluxes_broadband

   !> \section arg_table_ty_fluxes_broadband_ccpp Argument Table
   !! \htmlinclude ty_fluxes_broadband_ccpp.html
   type, public :: ty_fluxes_broadband_ccpp
      type(ty_fluxes_broadband) :: fluxes
   end type

end module ccpp_fluxes
