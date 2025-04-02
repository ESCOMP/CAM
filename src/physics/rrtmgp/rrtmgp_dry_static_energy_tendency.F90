module rrtmgp_dry_static_energy_tendency
!-----------------------------------------------------------------------
!
! Purpose:  Provide an interface to convert shortwave and longwave
!           radiative heating terms into net heating.
!
!           This module provides a hook to allow incorporating additional
!           radiative terms (eUV heating and nonLTE longwave cooling).
!
! Original version: B.A. Boville
!-----------------------------------------------------------------------

use ccpp_kinds,    only: kind_phys

implicit none
private
save

! Public interfaces
public :: rrtmgp_dry_static_energy_tendency_run

!===============================================================================
contains
!===============================================================================

!> \section arg_table_rrtmgp_dry_static_energy_tendency_run Argument Table
!! \htmlinclude rrtmgp_dry_static_energy_tendency_run.html
!!
subroutine rrtmgp_dry_static_energy_tendency_run(ncol, pdel, calc_sw_heat, calc_lw_heat, &
                qrs, qrl, errmsg, errflg)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------

   ! Arguments
   integer,                         intent(in)    :: ncol          ! Number of columns
   real(kind_phys), dimension(:,:), intent(in)    :: pdel          ! Layer thickness
   logical,                         intent(in)    :: calc_sw_heat  ! Flag to calculate net shortwave heating
   logical,                         intent(in)    :: calc_lw_heat  ! Flag to calculate net longwave heating
   real(kind_phys), dimension(:,:), intent(inout) :: qrs           ! shortwave heating rate (J kg-1 s-1)
   real(kind_phys), dimension(:,:), intent(inout) :: qrl           ! longwave heating rate  (J kg-1 s-1)
   character(len=*),                intent(out)   :: errmsg
   integer,                         intent(out)   :: errflg


   !-----------------------------------------------------------------------
   ! Set error variables
   errmsg = ''
   errflg = 0

   if (calc_sw_heat) then
      qrs(:ncol,:) = qrs(:ncol,:) / pdel(:ncol,:)
   end if

   if (calc_lw_heat) then
      qrl(:ncol,:) = qrl(:ncol,:) / pdel(:ncol,:)
   end if

end subroutine rrtmgp_dry_static_energy_tendency_run

!================================================================================================
end module rrtmgp_dry_static_energy_tendency
