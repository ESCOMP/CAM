!> \file rrtmgp_lw_main.F90
!! This file contains the core longwave RRTMGP radiation calcuation

!> This module contains the call to the RRTMGP-LW radiation routine
module rrtmgp_lw_main
  use machine,                  only: kind_phys
  use mo_rte_lw,                only: rte_lw
  use ccpp_fluxes_byband,       only: ty_fluxes_byband_ccpp
  use ccpp_optical_props,       only: ty_optical_props_1scl_ccpp
  use ccpp_fluxes_byband,       only: ty_fluxes_byband_ccpp
  use ccpp_fluxes,              only: ty_fluxes_broadband_ccpp
  use ccpp_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp_ccpp
  use ccpp_source_functions,    only: ty_source_func_lw_ccpp
  use radiation_tools,          only: check_error_msg
  implicit none

  public rrtmgp_lw_main_run
contains

!> \section arg_table_rrtmgp_lw_main_run Argument Table
!! \htmlinclude rrtmgp_lw_main_run.html
!!
   subroutine rrtmgp_lw_main_run(doLWrad, doLWclrsky, doGP_lwscat, use_LW_jacobian, use_LW_optimal_angles,   &
                                 nGauss_angles,  nCol, iter_num, rrtmgp_phys_blksz, lw_optical_props_clrsky, &
                                 lw_optical_props_clouds, top_at_1, sources, sfc_emiss_byband, lw_gas_props, &
                                 aerlw, fluxlwUP_jac, lw_Ds, flux_clrsky, flux_allsky, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: doLWrad               !< Flag to perform longwave calculation
    logical, intent(in) :: doLWclrsky            !< Flag to compute clear-sky fluxes
    logical, intent(in) :: doGP_lwscat           !< Flag to include scattering in clouds
    logical, intent(in) :: use_LW_jacobian       !< Flag to compute Jacobian
    logical, intent(in) :: use_LW_optimal_angles !< Flag to compute and use optimal angles
    logical, intent(in) :: top_at_1              !< Flag for vertical ordering convention

    integer, intent(in) :: nGauss_angles         !< Number of gaussian quadrature angles used
    integer, intent(in) :: nCol                  !< Number of horizontal points
    integer, intent(in) :: iter_num              !< Radiation subcycle iteration number
    integer, intent(in) :: rrtmgp_phys_blksz     !< Number of horizontal points to process at once

    real(kind_phys), dimension(:,:),   intent(in) :: sfc_emiss_byband           !< Surface emissivity by band
    class(ty_source_func_lw_ccpp),     intent(in) :: sources                    !< Longwave sources object

    ! Outputs
    real(kind_phys), dimension(:,:),   intent(inout) :: fluxlwUP_jac            !< Surface temperature flux Jacobian [W m-2 K-1]
    class(ty_fluxes_byband_ccpp),      intent(inout) :: flux_allsky             !< All-sky flux [W m-2]
    class(ty_fluxes_broadband_ccpp),   intent(inout) :: flux_clrsky             !< Clear-sky flux [W m-2]
    class(ty_optical_props_1scl_ccpp), intent(inout) :: aerlw                   !< Aerosol optical properties object
    class(ty_optical_props_1scl_ccpp), intent(inout) :: lw_optical_props_clrsky !< Clear-sky optical properties object
    class(ty_optical_props_1scl_ccpp), intent(inout) :: lw_optical_props_clouds !< Cloud optical properties object

    class(ty_gas_optics_rrtmgp_ccpp),  intent(inout) :: lw_gas_props            !< Gas optical properties object

    real(kind_phys), dimension(:,:),   intent(out) :: lw_Ds                     !< 1/cos of transport angle per column, g-point
    character(len=*), intent(out) :: errmsg                                     !< CCPP error message
    integer,          intent(out) :: errflg                                     !< CCPP error flag

    ! Local variables
    integer :: iCol, iCol2 

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    iCol = ((iter_num - 1) * rrtmgp_phys_blksz) + 1
    iCol2 = min(iCol + rrtmgp_phys_blksz - 1, nCol)

    ! ###################################################################################
    !
    ! Compute clear-sky fluxes (gaseous+aerosol) (optional)
    !
    ! ###################################################################################
    ! Increment
    errmsg = aerlw%optical_props%increment(lw_optical_props_clrsky%optical_props)
    call check_error_msg('rrtmgp_lw_main_increment_aerosol_to_clrsky', errmsg)
    if (len_trim(errmsg) =/ 0) then
        errflg = 1
        return
    end if

    ! Call RTE solver
    if (doLWclrsky) then
       if (nGauss_angles .gt. 1) then
          errmsg = rte_lw(           &
               lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
               top_at_1,                              & ! IN  - vertical ordering flag
               sources%sources,                       & ! IN  - source function
               sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
               flux_clrsky%fluxes,                    & ! OUT - Fluxes
               n_gauss_angles = nGauss_angles)          ! IN  - Number of angles in Gaussian quadrature
          call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky', errmsg)
          if (len_trim(errmsg) =/ 0) then
             errflg = 1
             return
          end if
       else
          if (use_lw_optimal_angles) then
             errmsg = lw_gas_props%gas_props%compute_optimal_angles(lw_optical_props_clrsky%optical_props,lw_Ds)
             call check_error_msg('rrtmgp_lw_main_opt_angle', errmsg)
             if (len_trim(errmsg) /= 0) then
                errflg = 1
                return
             end if
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_clrsky%fluxes,                    & ! OUT - Fluxes
                  lw_Ds = lw_Ds)
             call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
                return
             end if
          else
            errmsg = rte_lw(           &
                 lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                 top_at_1,                              & ! IN  - vertical ordering flag
                 sources%sources,                       & ! IN  - source function
                 sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                 flux_clrsky%fluxes)                      ! OUT - Fluxes
             call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
                return
             end if
          end if
       endif
    end if

    ! ###################################################################################
    !
    ! All-sky fluxes (clear-sky + clouds + precipitation)
    ! *Note* CCPP does not allow for polymorphic types, they are ambiguous to the CCPP
    ! framework. rte-rrtmgp uses polymorphic types extensively, for example, querying the
    ! type to determine physics configuration/pathway/etc...
    !
    ! The logic in the code below is to satisfy the polymorphishm in the rte-rrtmgp code.
    ! The rte-rrtmgp "increment" procedures are utilized to provide the correct type to the
    ! rte solver (rte_lw). Rte_lw quieries the type determine if scattering is to be 
    ! included in the calculation. The increment procedures are called so that the correct
    ! optical properties are inherited.
    ! 
    ! ###################################################################################

    ! Include LW cloud-scattering?
    if (doGP_lwscat) then 
       ! Increment
       errmsg = lw_optical_props_clrsky%optical_props%increment(lw_optical_props_clouds%optical_props)
       call check_error_msg('rrtmgp_lw_main_increment_clrsky_to_clouds', errmsg)
       if (len_trim(errmsg) =/ 0) then
           errflg = 1
           return
       end if
          
       if (use_LW_jacobian) then
          if (nGauss_angles .gt. 1) then
             ! Compute LW Jacobians; use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles,        & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
             end if
          else
             ! Compute LW Jacobians; don't use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
             end if
          end if
       else
          if (nGauss_angles .gt. 1) then
             ! Compute LW Jacobians; use Gaussian angles
             ! Don't compute LW Jacobians; use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles)          ! IN  - Number of angles in Gaussian quadrature
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
             end if
          else
             ! Don't compute LW Jacobians; don't use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes)                      ! OUT - Fluxes
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
             end if
          end if
       end if
    ! No scattering in LW clouds.   
    else
       ! Increment
       errmsg = lw_optical_props_clouds%optical_props%increment(lw_optical_props_clrsky%optical_props)
       call check_error_msg('rrtmgp_lw_main_increment_clouds_to_clrsky', errmsg)
       if (len_trim(errmsg) =/ 0) then
           errflg = 1
           return
       end if
          
       if (use_LW_jacobian) then
          if (nGauss_angles .gt. 1) then
             ! Compute LW Jacobians; use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles,        & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
             end if
          else
             ! Compute LW Jacobians; don't use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
             end if
          end if
       else
          if (nGauss_angles .gt. 1) then
             ! Don't compute LW Jacobians; use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles)          ! IN  - Number of angles in Gaussian quadrature
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) =/ 0) then
                errflg = 1
             end if
          else
             ! Don't compute LW Jacobians; don't use Gaussian angles
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  top_at_1,                              & ! IN  - vertical ordering flag
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes)                      ! OUT - Fluxes
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky', errmsg)
             if (len_trim(errmsg) /= 0) then
                errflg = 1
             end if
          end if
       end if
     end if

  end subroutine rrtmgp_lw_main_run
end module rrtmgp_lw_main
