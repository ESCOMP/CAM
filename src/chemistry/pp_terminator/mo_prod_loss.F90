      module mo_prod_loss
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: exp_prod_loss
      public :: imp_prod_loss
      contains
      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )
      use ppgrid, only : pver
      implicit none
!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:,:,:)
      real(r8), intent(in) :: rxt(:,:,:)
      real(r8), intent(in) :: het_rates(:,:,:)
      end subroutine exp_prod_loss
      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )
      use ppgrid, only : pver
      implicit none
!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:)
      real(r8), intent(in) :: rxt(:)
      real(r8), intent(in) :: het_rates(:)
!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------
         loss(3) = (2._r8*rxt(2)* y(1) + het_rates(1))* y(1)
         prod(3) =2.000_r8*rxt(1)*y(2)
         loss(2) = ( + rxt(1) + het_rates(2))* y(2)
         prod(2) =rxt(2)*y(1)*y(1)
         loss(1) = ( + het_rates(3))* y(3)
         prod(1) = 0._r8
      end subroutine imp_prod_loss
      end module mo_prod_loss
