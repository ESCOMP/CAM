!> \file rrtmgp_lw_gas_optics.F90
!!

!> This module contains a run routine to compute gas optics during the radiation subcycle
module rrtmgp_lw_gas_optics
  use machine,               only: kind_phys
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_1scl
  use mo_source_functions,   only: ty_source_func_lw
  use radiation_tools,       only: check_error_msg

  implicit none

  public :: rrtmgp_lw_gas_optics_run
contains

!> \section arg_table_rrtmgp_lw_gas_optics_run Argument Table
!! \htmlinclude rrtmgp_lw_gas_optics_run.html
!!
  subroutine rrtmgp_lw_gas_optics_run(dolw, iter_num, ncol, rrtmgp_phys_blksz, p_lay, p_lev, t_lay, tsfg, &
             gas_concs, lw_optical_props_clrsky, sources, t_lev, include_interface_temp, lw_gas_props,    &
             errmsg, errflg)
   ! Inputs
   logical,                         intent(in) :: dolw
   logical,                         intent(in) :: include_interface_temp
   integer,                         intent(in) :: iter_num
   integer,                         intent(in) :: ncol
   integer,                         intent(in) :: rrtmgp_phys_blksz
   real(kind_phys), dimension(:,:), intent(in) :: p_lay
   real(kind_phys), dimension(:,:), intent(in) :: p_lev
   real(kind_phys), dimension(:,:), intent(in) :: t_lay
   real(kind_phys), dimension(:),   intent(in) :: tsfg
   real(kind_phys), dimension(:,:), intent(in) :: t_lev
   type(ty_gas_concs),           intent(in)    :: gas_concs                !< RRTMGP gas concentrations object

   ! Outputs
   !type(ty_gas_concs),           intent(in)    :: gas_concs                !< RRTMGP gas concentrations object
   type(ty_optical_props_1scl),  intent(inout) :: lw_optical_props_clrsky  !< Clearsky optical properties
   type(ty_source_func_lw),      intent(inout) :: sources
   character(len=*),             intent(out)   :: errmsg
   integer,                      intent(out)   :: errflg
   type(ty_gas_optics_rrtmgp),   intent(inout) :: lw_gas_props             !< RRTMGP gas optics object

   ! Local variables
   integer :: iCol, iCol2

   ! Set error variables
   errmsg = ''
   errflg = 0

   if (.not. dolw) then
      return
   end if

   iCol = ((iter_num - 1) * rrtmgp_phys_blksz) + 1
   iCol2= min(iCol + rrtmgp_phys_blksz - 1, ncol)
   if (include_interface_temp) then
      call check_error_msg('rrtmgp_lw_main_gas_optics',lw_gas_props%gas_optics(&
            p_lay(iCol:iCol2,:),              & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(iCol:iCol2,:),              & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(iCol:iCol2,:),              & ! IN  - Temperature @ layer-centers (K)
            tsfg(iCol:iCol2),                 & ! IN  - Skin-temperature (K)
            gas_concs,                        & ! IN  - RRTMGP DDT: trace gas volume mixing-ratios
            lw_optical_props_clrsky,          & ! OUT - RRTMGP DDT: longwave optical properties
            sources,                          & ! OUT - RRTMGP DDT: source functions
            tlev=t_lev(iCol:iCol2,:)))          ! IN  - Temperature @ layer-interfaces (K) (optional)
   else
      call check_error_msg('rrtmgp_lw_main_gas_optics',lw_gas_props%gas_optics(&
            p_lay(iCol:iCol2,:),              & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(iCol:iCol2,:),              & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(iCol:iCol2,:),              & ! IN  - Temperature @ layer-centers (K)
            tsfg(iCol:iCol2),                 & ! IN  - Skin-temperature (K)
            gas_concs,                        & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            lw_optical_props_clrsky,          & ! OUT - RRTMGP DDT: longwave optical properties
            sources))                           ! OUT - RRTMGP DDT: source functions
   end if

  end subroutine rrtmgp_lw_gas_optics_run

end module rrtmgp_lw_gas_optics
