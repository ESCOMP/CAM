      module mo_adjrxt
      private
      public :: adjrxt
      contains
      subroutine adjrxt( rate, inv, m, ncol, nlev )
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : nfs, rxntot
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in) :: inv(ncol,nlev,nfs)
      real(r8), intent(in) :: m(ncol,nlev)
      real(r8), intent(inout) :: rate(ncol,nlev,rxntot)
      rate(:,:, 2) = rate(:,:, 2) * m(:,:)
      end subroutine adjrxt
      end module mo_adjrxt
