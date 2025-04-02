!> \file rrtmgp_lw_gas_optics.F90
!!

!> This module contains a run routine to compute gas optics during the radiation subcycle
module rrtmgp_lw_gas_optics
  use machine,                 only: kind_phys
  use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
  use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
  use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp
  use ccpp_source_functions,   only: ty_source_func_lw_ccpp
  use radiation_tools,         only: check_error_msg

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
   logical,                           intent(in) :: dolw                        !< Flag for whether to perform longwave calculation
   logical,                           intent(in) :: include_interface_temp      !< Flag for whether to include interface temperature in calculation
   integer,                           intent(in) :: iter_num                    !< Subcycle iteration number
   integer,                           intent(in) :: ncol                        !< Total number of columns
   integer,                           intent(in) :: rrtmgp_phys_blksz           !< Number of horizontal points to process at once
   real(kind_phys), dimension(:,:),   intent(in) :: p_lay                       !< Air pressure at midpoints [Pa]
   real(kind_phys), dimension(:,:),   intent(in) :: p_lev                       !< Air pressure at interfaces [Pa]
   real(kind_phys), dimension(:,:),   intent(in) :: t_lay                       !< Air temperature at midpoints [K]
   real(kind_phys), dimension(:),     intent(in) :: tsfg                        !< Surface skin temperature [K]
   real(kind_phys), dimension(:,:),   intent(in) :: t_lev                       !< Air temperature at interfaces [K]
   type(ty_gas_concs_ccpp),           intent(in) :: gas_concs                   !< RRTMGP gas concentrations object

   ! Outputs
   type(ty_optical_props_1scl_ccpp),  intent(inout) :: lw_optical_props_clrsky  !< Clearsky optical properties
   type(ty_gas_optics_rrtmgp_ccpp),   intent(inout) :: lw_gas_props             !< RRTMGP gas optics object
   type(ty_source_func_lw_ccpp),      intent(inout) :: sources                  !< Longwave sources object
   character(len=*),                  intent(out)   :: errmsg
   integer,                           intent(out)   :: errflg

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
      errmsg = lw_gas_props%gas_props%gas_optics(&
            p_lay(iCol:iCol2,:),                   & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(iCol:iCol2,:),                   & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(iCol:iCol2,:),                   & ! IN  - Temperature @ layer-centers (K)
            tsfg(iCol:iCol2),                      & ! IN  - Skin-temperature (K)
            gas_concs%gas_concs,                   & ! IN  - RRTMGP DDT: trace gas volume mixing-ratios
            lw_optical_props_clrsky%optical_props, & ! OUT - RRTMGP DDT: longwave optical properties
            sources%sources,                       & ! OUT - RRTMGP DDT: source functions
            tlev=t_lev(iCol:iCol2,:))                ! IN  - Temperature @ layer-interfaces (K) (optional)
      call check_error_msg('rrtmgp_lw_main_gas_optics', errmsg)
      if (len_trim(errmsg) /= 0) then
         errflg = 1
      end if
   else
      errmsg = lw_gas_props%gas_props%gas_optics(&
            p_lay(iCol:iCol2,:),                   & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(iCol:iCol2,:),                   & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(iCol:iCol2,:),                   & ! IN  - Temperature @ layer-centers (K)
            tsfg(iCol:iCol2),                      & ! IN  - Skin-temperature (K)
            gas_concs%gas_concs,                   & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            lw_optical_props_clrsky%optical_props, & ! OUT - RRTMGP DDT: longwave optical properties
            sources%sources)                         ! OUT - RRTMGP DDT: source functions
      call check_error_msg('rrtmgp_lw_main_gas_optics', errmsg)
      if (len_trim(errmsg) /= 0) then
         errflg = 1
      end if
   end if

  end subroutine rrtmgp_lw_gas_optics_run

end module rrtmgp_lw_gas_optics
