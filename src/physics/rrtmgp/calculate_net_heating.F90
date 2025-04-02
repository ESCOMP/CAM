module calculate_net_heating
! PEVERWHEE - this should go in schemes/rrtmgp/utils
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
public :: calculate_net_heating_run

!===============================================================================
contains
!===============================================================================

!> \section arg_table_calculate_net_heating_run Argument Table
!! \htmlinclude calculate_net_heating_run.html
!!
subroutine calculate_net_heating_run(ncol, rad_heat, qrl, qrs, fsns, fsnt, flns, flnt, &
                is_offline_dyn, net_flx, errmsg, errflg)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------

    ! Arguments
   integer,                    intent(in)  :: ncol           ! horizontal dimension
   real(kind_phys),            intent(in)  :: qrl(:,:)       ! longwave heating [J kg-1 s-1]
   real(kind_phys),            intent(in)  :: qrs(:,:)       ! shortwave heating [J kg-1 s-1]
   real(kind_phys),            intent(in)  :: fsns(:)        ! Surface solar absorbed flux [W m-2]
   real(kind_phys),            intent(in)  :: fsnt(:)        ! Net column abs solar flux at model top [W m-2]
   real(kind_phys),            intent(in)  :: flns(:)        ! Srf longwave cooling (up-down) flux [W m-2]
   real(kind_phys),            intent(in)  :: flnt(:)        ! Net outgoing lw flux at model top [W m-2]
   logical,                    intent(in)  :: is_offline_dyn ! is offline dycore
   real(kind_phys),            intent(out) :: rad_heat(:,:)  ! radiative heating [J kg-1 s-1]
   real(kind_phys),            intent(out) :: net_flx(:)     ! net boundary flux [W m-2]
   character(len=*),           intent(out) :: errmsg
   integer,                    intent(out) :: errflg


   ! Local variables
   integer :: idx
   !-----------------------------------------------------------------------
   ! Set error variables
   errmsg = ''
   errflg = 0
   if (.not. is_offline_dyn) then
      rad_heat(:ncol,:) = (qrs(:ncol,:) + qrl(:ncol,:))
   end if

   do idx = 1, ncol
      net_flx(idx) = fsnt(idx) - fsns(idx) - flnt(idx) + flns(idx)
   end do

end subroutine calculate_net_heating_run

!================================================================================================
end module calculate_net_heating
