module ccpp_fluxes_byband
   ! CCPP wrapper for ty_fluxes_byband DDT from RRTMGP
   use mo_fluxes_byband, only: ty_fluxes_byband
   use ccpp_fluxes,      only: ty_fluxes_broadband_ccpp

   !> \section arg_table_ty_fluxes_byband_ccpp Argument Table
   !! \htmlinclude ty_fluxes_byband_ccpp.html
   type, public :: ty_fluxes_byband_ccpp
      type(ty_fluxes_byband) :: fluxes
   end type

end module ccpp_fluxes_byband
