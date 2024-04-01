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
!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      real(r8) :: im(ncol,nlev)
      im(:,:) = 1._r8 / m(:,:)
      rate(:,:, 5) = rate(:,:, 5) * inv(:,:, 5)
      rate(:,:, 11) = rate(:,:, 11) * inv(:,:, 6)
      rate(:,:, 12) = rate(:,:, 12) * inv(:,:, 5)
      rate(:,:, 14) = rate(:,:, 14) * inv(:,:, 5)
      rate(:,:, 6) = rate(:,:, 6) * inv(:,:, 4) * inv(:,:, 4) * im(:,:)
      rate(:,:, 13) = rate(:,:, 13) * inv(:,:, 5) * inv(:,:, 1)
      end subroutine adjrxt
      end module mo_adjrxt
