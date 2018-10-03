module trsolv_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  
contains
!-----------------------------------------------------------------------
  subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1)
!
! Tri-diagonal solver.
!   a(k,i)*x(k-1,i) + b(k,i)*x(k,i) + c(k,i)*x(k+1,i) = f(k,i)
!
    implicit none
!
! Args:
    integer,intent(in) :: lev0,lev1,k1,k2,lon0,lon1
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in) :: &
      a, & ! input coefficients
      b, & ! input coefficients
      c, & ! input coefficients
      f    ! input RHS
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: &
      x  ! output
!
! Local:
    integer :: k,kk,i
    real(r8),dimension(lev0:lev1,lon0:lon1) :: w1,w2,w3  ! work arrays

!
! Lower boundary (W(K1)=B(K1):
    do i=lon0,lon1
      w1(lev0,i) = b(lev0,i) 
    enddo
!
! Set up work arrays:
    do i=lon0,lon1
      do k=k1+1,k2
!
! W(KF+K-1)=C(K-1)/W(K-1):
        w2(k-1,i) = c(k-1,i) / w1(k-1,i)
!
! W(K)=A(K)*W(KF+K-1)
        w1(k,i) = a(k,i) * w2(k-1,i)
!
! W(K)=B(K)-W(K)
        w1(k,i) = b(k,i) - w1(k,i)
      enddo ! k=k1+1,k2
    enddo ! i=lon0,lon1
!
! Lower boundary (W(2*KF+K1)=F(K1)/W(K1)):
    do i=lon0,lon1
      w3(k1,i) = f(k1,i) / w1(k1,i)
    enddo
!
    do i=lon0,lon1
      do k=k1+1,k2
!
! W(2*KF+K)=A(K)*W(2*KF+K-1)
        w3(k,i) = a(k,i) * w3(k-1,i)
!
! W(2*KF+K)=F(K)-W(2*KF+K)
        w3(k,i) = f(k,i) - w3(k,i)         
!
! W(2*KF+K)=W(2*KF+K)/W(K)
        w3(k,i) = w3(k,i) / w1(k,i)
      enddo ! k=k1+1,k2
    enddo ! i=lon0,lon1
!
! Upper boundary (X(K2)=W(2*KF+K2)):
    do i=lon0,lon1
      x(k2,i) = w3(k2,i)       
    enddo
!

! Back substitution:
    do i=lon0,lon1
      do kk=k1+1,k2          
        k = k1+k2-kk ! k2-1,k1,-1
!
! X(K)=W(KF+K)*X(K+1)
        x(k,i) = w2(k,i) * x(k+1,i)
!
! X(K)=W(2*KF+K)-X(K)
        x(k,i) = w3(k,i) - x(k,i)
      enddo ! k=k1+1,k2
    enddo
  end subroutine trsolv
!-----------------------------------------------------------------------
end module trsolv_mod
