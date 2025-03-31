!> \file rrtmgp_lw_main.F90
!! This file contains the longwave RRTMGP radiation scheme.

!> This module contains the call to the RRTMGP-LW radiation scheme
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
    logical, intent(in) :: doLWrad               ! Flag to perform longwave calculation
    logical, intent(in) :: doLWclrsky            ! Flag to compute clear-sky fluxes
    logical, intent(in) :: doGP_lwscat           ! Flag to include scattering in clouds
    logical, intent(in) :: use_LW_jacobian       ! Flag to compute Jacobian
    logical, intent(in) :: use_LW_optimal_angles ! Flag to compute and use optimal angles
    logical, intent(in) :: top_at_1              ! Flag for vertical ordering convention

    integer, intent(in) :: nGauss_angles         ! Number of gaussian quadrature angles used
    integer, intent(in) :: nCol                  ! Number of horizontal points
    integer, intent(in) :: iter_num              ! RRTMGP iteration number
    integer, intent(in) :: rrtmgp_phys_blksz     ! Number of horizontal points to process at once

    real(kind_phys), dimension(:,:), intent(out) :: lw_Ds
    real(kind_phys), dimension(:,:), intent(in) :: sfc_emiss_byband

    class(ty_source_func_lw_ccpp), intent(in) :: sources

    ! Outputs
    real(kind_phys), dimension(:,:),   intent(inout) :: fluxlwUP_jac
    class(ty_fluxes_byband_ccpp),             intent(inout) :: flux_allsky
    class(ty_fluxes_broadband_ccpp),             intent(inout) :: flux_clrsky
    class(ty_optical_props_1scl_ccpp), intent(inout) :: aerlw
    class(ty_optical_props_1scl_ccpp), intent(inout) :: lw_optical_props_clrsky
    class(ty_optical_props_1scl_ccpp), intent(inout) :: lw_optical_props_clouds

    class(ty_gas_optics_rrtmgp_ccpp),   intent(inout) :: lw_gas_props

    character(len=*), intent(out) :: errmsg     ! CCPP error message
    integer,          intent(out) :: errflg     ! CCPP error flag

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
    call check_error_msg('rrtmgp_lw_main_increment_aerosol_to_clrsky',&
         aerlw%optical_props%increment(lw_optical_props_clrsky%optical_props))

    ! Call RTE solver
    if (doLWclrsky) then
       if (nGauss_angles .gt. 1) then
          call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky',rte_lw(           &
               lw_optical_props_clrsky%optical_props,         & ! IN  - optical-properties
               top_at_1,                        & ! IN  - veritcal ordering flag
               sources%sources,                 & ! IN  - source function
               sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
               flux_clrsky%fluxes,              & ! OUT - Fluxes
               n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature
       else
          if (use_lw_optimal_angles) then
             call check_error_msg('rrtmgp_lw_main_opt_angle',&
                  lw_gas_props%gas_props%compute_optimal_angles(lw_optical_props_clrsky%optical_props,lw_Ds))
             call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky',rte_lw(           &
                  lw_optical_props_clrsky%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_clrsky%fluxes,              & ! OUT - Fluxes
                  lw_Ds = lw_Ds))
          else
            call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky',rte_lw(           &
                 lw_optical_props_clrsky%optical_props,         & ! IN  - optical-properties
                 top_at_1,                        & ! IN  - veritcal ordering flag
                 sources%sources,                 & ! IN  - source function
                 sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                 flux_clrsky%fluxes))               ! OUT - Fluxes
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
       call check_error_msg('rrtmgp_lw_main_increment_clrsky_to_clouds',&
            lw_optical_props_clrsky%optical_props%increment(lw_optical_props_clouds%optical_props))
          
       if (use_LW_jacobian) then
          if (nGauss_angles .gt. 1) then
             ! Compute LW Jacobians; use Gaussian angles
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky',rte_lw(           &
                  lw_optical_props_clouds%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,              & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles,  & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac))    ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
          else
             ! Compute LW Jacobians; don't use Gaussian angles
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky',rte_lw(           &
                  lw_optical_props_clouds%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,              & ! OUT - Fluxes
                  flux_up_Jac    = fluxlwUP_jac))    ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
          end if
       else
          if (nGauss_angles .gt. 1) then
             ! Compute LW Jacobians; use Gaussian angles
             ! Don't compute LW Jacobians; use Gaussian angles
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky',rte_lw(           &
                  lw_optical_props_clouds%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,              & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature
          else
             ! Don't compute LW Jacobians; don't use Gaussian angles
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky',rte_lw(           &
                  lw_optical_props_clouds%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes))               ! OUT - Fluxes
          end if
       end if
    ! No scattering in LW clouds.   
    else
       ! Increment
       call check_error_msg('rrtmgp_lw_main_increment_clouds_to_clrsky', &
            lw_optical_props_clouds%optical_props%increment(lw_optical_props_clrsky%optical_props))
          
       if (use_LW_jacobian) then
          if (nGauss_angles .gt. 1) then
             ! Compute LW Jacobians; use Gaussian angles
             call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
                  lw_optical_props_clrsky%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,              & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles,  & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac))    ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
          else
             ! Compute LW Jacobians; don't use Gaussian angles
             call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
                  lw_optical_props_clrsky%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,              & ! OUT - Fluxes
                  flux_up_Jac    = fluxlwUP_jac))    ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
          end if
       else
          if (nGauss_angles .gt. 1) then
             ! Don't compute LW Jacobians; use Gaussian angles
             call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
                  lw_optical_props_clrsky%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,              & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature
          else
             ! Don't compute LW Jacobians; don't use Gaussian angles
             call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
                  lw_optical_props_clrsky%optical_props,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources%sources,                 & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes))               ! OUT - Fluxes
          end if
       end if
     end if

  end subroutine rrtmgp_lw_main_run
end module rrtmgp_lw_main
