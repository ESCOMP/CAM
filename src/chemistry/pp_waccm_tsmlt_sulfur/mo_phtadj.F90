      module mo_phtadj
      private
      public :: phtadj
      contains
      subroutine phtadj( p_rate, inv, m, ncol, nlev )
      use chem_mods, only : nfs, phtcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in) :: inv(ncol,nlev,max(1,nfs))
      real(r8), intent(in) :: m(ncol,nlev)
      real(r8), intent(inout) :: p_rate(ncol,nlev,max(1,phtcnt))
!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      integer :: k
      real(r8) :: im(ncol,nlev)
      do k = 1,nlev
         im(:ncol,k) = 1._r8 / m(:ncol,k)
         p_rate(:,k, 97) = p_rate(:,k, 97) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,101) = p_rate(:,k,101) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,102) = p_rate(:,k,102) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,104) = p_rate(:,k,104) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,109) = p_rate(:,k,109) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,113) = p_rate(:,k,113) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,114) = p_rate(:,k,114) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,116) = p_rate(:,k,116) * inv(:,k, 2) * im(:,k)
      end do
      end subroutine phtadj
      end module mo_phtadj
