
!**************************************************************************************
!
! fv_mapz contains vertical remapping algorithms that come from the FV3 dycore.
! They have been minimally modified for use in CAM.
!
! The following license statement is from the original code.
!
!**************************************************************************************
  
!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it 
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be 
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.  
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module fv_mapz
  
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use cam_abortutils,         only: endrun
  
  implicit none
  
  public map_scalar, map1_ppm, mapn_tracer
  
  real(kind=r8), parameter::  r3 = 1._r8/3._r8, r23 = 2._r8/3._r8, r12 = 1._r8/12._r8
contains
  
  subroutine map_scalar( km,   pe1,    q1,   qs,           &
       kn,   pe2,    q2,   i1, i2,       &
       j,  ibeg, iend, jbeg, jend, iv,  kord, q_min)
    ! iv=1
    integer, intent(in) :: i1                !< Starting longitude
    integer, intent(in) :: i2                !< Finishing longitude
    integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == temp 2 == remap temp with cs scheme
    integer, intent(in) :: kord              !< Method order
    integer, intent(in) :: j                 !< Current latitude
    integer, intent(in) :: ibeg, iend, jbeg, jend
    integer, intent(in) :: km                !< Original vertical dimension
    integer, intent(in) :: kn                !< Target vertical dimension
    real(kind=r8), intent(in) ::   qs(i1:i2)       !< bottom BC
    real(kind=r8), intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
    real(kind=r8), intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
    real(kind=r8), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
    ! INPUT/OUTPUT PARAMETERS:
    real(kind=r8), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output
    real(kind=r8), intent(in):: q_min
    
    ! DESCRIPTION:
    ! IV = 0: constituents
    ! pe1: pressure at layer edges (from model top to bottom surface)
    !      in the original vertical coordinate
    ! pe2: pressure at layer edges (from model top to bottom surface)
    !      in the new vertical coordinate
    ! LOCAL VARIABLES:
    real(kind=r8)    dp1(i1:i2,km)
    real(kind=r8)   q4(4,i1:i2,km)
    real(kind=r8)    pl, pr, qsum, dp, esl
    integer i, k, l, m, k0
    
    do k=1,km
      do i=i1,i2
        dp1(i,k) = pe1(i,k+1) - pe1(i,k)
        q4(1,i,k) = q1(i,j,k)
      enddo
    enddo
    
    ! Compute vertical subgrid distribution
    if ( kord >7 ) then
      call scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
    else
      call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
    endif
    
    do i=i1,i2
      k0 = 1
      do 555 k=1,kn
        do l=k0,km
          ! locate the top edge: pe2(i,k)
          if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
            pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
            if( pe2(i,k+1) <= pe1(i,l+1) ) then
              ! entire new grid is within the original grid
              pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
              q2(i,j,k) = q4(2,i,l) + 0.5_r8*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                   *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              goto 555
            else
              ! Fractional area...
              qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5_r8*(q4(4,i,l)+   &
                   q4(3,i,l)-q4(2,i,l))*(1._r8+pl)-q4(4,i,l)*           &
                   (r3*(1._r8+pl*(1._r8+pl))))
              do m=l+1,km
                ! locate the bottom edge: pe2(i,k+1)
                if( pe2(i,k+1) > pe1(i,m+1) ) then
                  ! Whole layer
                  qsum = qsum + dp1(i,m)*q4(1,i,m)
                else
                  dp = pe2(i,k+1)-pe1(i,m)
                  esl = dp / dp1(i,m)
                  qsum = qsum + dp*(q4(2,i,m)+0.5_r8*esl*               &
                       (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1._r8-r23*esl)))
                  k0 = m
                  goto 123
                endif
              enddo
              goto 123
            endif
          endif
        enddo
123     q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555     continue
      enddo
  end subroutine map_scalar
    
  
  subroutine mapn_tracer(nq, km, pe1, pe2, q1, dp2, kord, j,     &
       i1, i2, isd, ied, jsd, jed, q_min, fill)
    ! INPUT PARAMETERS:
    integer, intent(in):: km                !< vertical dimension
    integer, intent(in):: j, nq, i1, i2
    integer, intent(in):: isd, ied, jsd, jed
    integer, intent(in):: kord(nq)
    real(kind=r8), intent(in)::  pe1(i1:i2,km+1)     !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
    real(kind=r8), intent(in)::  pe2(i1:i2,km+1)     !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
    real(kind=r8), intent(in)::  dp2(i1:i2,km)
    real(kind=r8), intent(in)::  q_min
    logical, intent(in):: fill
    real(kind=r8), intent(inout):: q1(isd:ied,jsd:jed,km,nq) ! Field input
    ! LOCAL VARIABLES:
    real(kind=r8):: q4(4,i1:i2,km,nq)
    real(kind=r8):: q2(i1:i2,km,nq) !< Field output
    real(kind=r8):: qsum(nq)
    real(kind=r8):: dp1(i1:i2,km)
    real(kind=r8):: qs(i1:i2)
    real(kind=r8):: pl, pr, dp, esl, fac1, fac2
    integer:: i, k, l, m, k0, iq
    
    do k=1,km
      do i=i1,i2
        dp1(i,k) = pe1(i,k+1) - pe1(i,k)
      enddo
    enddo
    
    do iq=1,nq
      do k=1,km
        do i=i1,i2
          q4(1,i,k,iq) = q1(i,j,k,iq)
        enddo
      enddo
      call scalar_profile( qs, q4(1,i1,1,iq), dp1, km, i1, i2, 0, kord(iq), q_min )
    enddo
    ! Mapping
    do 1000 i=i1,i2
      k0 = 1
      do 555 k=1,km
        do 100 l=k0,km
          ! locate the top edge: pe2(i,k)
          if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
            pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
            if(pe2(i,k+1) <= pe1(i,l+1)) then
              ! entire new grid is within the original grid
              pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
              fac1 = pr + pl
              fac2 = r3*(pr*fac1 + pl*pl) 
              fac1 = 0.5_r8*fac1
              do iq=1,nq
                q2(i,k,iq) = q4(2,i,l,iq) + (q4(4,i,l,iq)+q4(3,i,l,iq)-q4(2,i,l,iq))*fac1  &
                     -  q4(4,i,l,iq)*fac2
              enddo
              k0 = l
              goto 555
            else
              ! Fractional area...
              dp = pe1(i,l+1) - pe2(i,k)
              fac1 = 1.0_r8 + pl
              fac2 = r3*(1.0_r8+pl*fac1)
              fac1 = 0.5_r8*fac1
              do iq=1,nq
                qsum(iq) = dp*(q4(2,i,l,iq) + (q4(4,i,l,iq)+   &
                     q4(3,i,l,iq) - q4(2,i,l,iq))*fac1 - q4(4,i,l,iq)*fac2)
              enddo
              do m=l+1,km
                ! locate the bottom edge: pe2(i,k+1)
                if(pe2(i,k+1) > pe1(i,m+1) ) then
                  ! Whole layer..
                  do iq=1,nq
                    qsum(iq) = qsum(iq) + dp1(i,m)*q4(1,i,m,iq)
                  enddo
                else
                  dp = pe2(i,k+1)-pe1(i,m)
                  esl = dp / dp1(i,m)
                  fac1 = 0.5_r8*esl
                  fac2 = 1.0_r8-r23*esl
                  do iq=1,nq
                    qsum(iq) = qsum(iq) + dp*( q4(2,i,m,iq) + fac1*(         &
                         q4(3,i,m,iq)-q4(2,i,m,iq)+q4(4,i,m,iq)*fac2 ) )
                  enddo
                  k0 = m
                  goto 123
                endif
              enddo
              goto 123
            endif
          endif
100     continue
123   continue
      do iq=1,nq
        q2(i,k,iq) = qsum(iq) / dp2(i,k)
      enddo
555   continue
1000  continue
          
      if (fill) call fillz(i2-i1+1, km, nq, q2, dp2)
      
      do iq=1,nq
        !    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
        do k=1,km
          do i=i1,i2
            q1(i,j,k,iq) = q2(i,k,iq)
          enddo
        enddo
      enddo

    end subroutine mapn_tracer


    subroutine map1_ppm( km,   pe1,    q1,   qs,           &
         kn,   pe2,    q2,   i1, i2,       &
         j,    ibeg, iend, jbeg, jend, iv,  kord)
      integer, intent(in) :: i1                !< Starting longitude
      integer, intent(in) :: i2                !< Finishing longitude
      integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == ??? 2 == remap temp with cs scheme
      integer, intent(in) :: kord              !< Method order
      integer, intent(in) :: j                 !< Current latitude
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: km                !< Original vertical dimension
      integer, intent(in) :: kn                !< Target vertical dimension
      real(kind=r8), intent(in) ::   qs(i1:i2)       !< bottom BC
      real(kind=r8), intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
      real(kind=r8), intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
      real(kind=r8), intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
      ! INPUT/OUTPUT PARAMETERS:
      real(kind=r8), intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output
      
      ! DESCRIPTION:
      ! IV = 0: constituents
      ! pe1: pressure at layer edges (from model top to bottom surface)
      !      in the original vertical coordinate
      ! pe2: pressure at layer edges (from model top to bottom surface)
      !      in the new vertical coordinate
      
      ! LOCAL VARIABLES:
      real(kind=r8)    dp1(i1:i2,km)
      real(kind=r8)   q4(4,i1:i2,km)
      real(kind=r8)    pl, pr, qsum, dp, esl
      integer i, k, l, m, k0
      
      do k=1,km
        do i=i1,i2
          dp1(i,k) = pe1(i,k+1) - pe1(i,k)
          q4(1,i,k) = q1(i,j,k)
        enddo
      enddo
      
      ! Compute vertical subgrid distribution
      if ( kord >7 ) then
        call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
      else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
      endif
      
      do i=i1,i2
        k0 = 1
        do 555 k=1,kn
          do l=k0,km
            ! locate the top edge: pe2(i,k)
            if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
              pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
              if( pe2(i,k+1) <= pe1(i,l+1) ) then
                ! entire new grid is within the original grid
                pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
                q2(i,j,k) = q4(2,i,l) + 0.5_r8*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                     *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
                k0 = l
                goto 555
              else
                ! Fractional area...
                qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5_r8*(q4(4,i,l)+   &
                     q4(3,i,l)-q4(2,i,l))*(1.0_r8+pl)-q4(4,i,l)*           &
                     (r3*(1.0_r8+pl*(1.0_r8+pl))))
                do m=l+1,km
                  ! locate the bottom edge: pe2(i,k+1)
                  if( pe2(i,k+1) > pe1(i,m+1) ) then
                    ! Whole layer
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                  else
                    dp = pe2(i,k+1)-pe1(i,m)
                    esl = dp / dp1(i,m)
                    qsum = qsum + dp*(q4(2,i,m)+0.5_r8*esl*               &
                         (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.0_r8-r23*esl)))
                    k0 = m
                    goto 123
                  endif
                enddo
                goto 123
              endif
            endif
          enddo
123       q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555       continue
        enddo

    end subroutine map1_ppm

    subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)
      
      ! INPUT PARAMETERS:
      integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others iv = 2: temp (if remap_t) and w (iv=-2)
      integer, intent(in):: i1      !< Starting longitude
      integer, intent(in):: i2      !< Finishing longitude
      integer, intent(in):: km      !< Vertical dimension
      integer, intent(in):: kord    !< Order (or more accurately method no.):
      ! 
      real(kind=r8) , intent(in):: delp(i1:i2,km)     !< Layer pressure thickness
      
      ! !INPUT/OUTPUT PARAMETERS:
      real(kind=r8) , intent(inout):: a4(4,i1:i2,km)  !< Interpolated values
      
      ! DESCRIPTION:
      !
      !   Perform the piecewise parabolic reconstruction
      ! 
      ! !REVISION HISTORY: 
      ! S.-J. Lin   revised at GFDL 2007
      !-----------------------------------------------------------------------
      ! local arrays:
      real(kind=r8)    dc(i1:i2,km)
      real(kind=r8)    h2(i1:i2,km)
      real(kind=r8)  delq(i1:i2,km)
      real(kind=r8)   df2(i1:i2,km)
      real(kind=r8)    d4(i1:i2,km)
      
      ! local scalars:
      integer i, k, km1, lmt, it
      real(kind=r8)  fac
      real(kind=r8)  a1, a2, c1, c2, c3, d1, d2
      real(kind=r8)  qm, dq, lac, qmp, pmp
      
      km1 = km - 1
      it = i2 - i1 + 1
      
      do k=2,km
        do i=i1,i2
          delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
          d4(i,k  ) = delp(i,k-1) + delp(i,k)
        enddo
      enddo
      
      do k=2,km1
        do i=i1,i2
          c1  = (delp(i,k-1)+0.5_r8*delp(i,k))/d4(i,k+1)
          c2  = (delp(i,k+1)+0.5_r8*delp(i,k))/d4(i,k)
          df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
               (d4(i,k)+delp(i,k+1))
          dc(i,k) = sign( min(abs(df2(i,k)),              &
               max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
               a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
        enddo
      enddo
      
      !-----------------------------------------------------------
      ! 4th order interpolation of the provisional cell edge value
      !-----------------------------------------------------------
      
      do k=3,km1
        do i=i1,i2
          c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
          a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
          a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
          a4(2,i,k) = a4(1,i,k-1) + c1 + 2.0_r8/(d4(i,k-1)+d4(i,k+1)) *    &
               ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
               delp(i,k-1)*a1*dc(i,k  ) )
        enddo
      enddo
      
      !     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)
      
      ! Area preserving cubic with 2nd deriv. = 0 at the boundaries
      ! Top
      do i=i1,i2
        d1 = delp(i,1)
        d2 = delp(i,2)
        qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
        dq = 2.0_r8*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
        c1 = 4.0_r8*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.0_r8*d2*d2+d1*(d2+3.0_r8*d1)) )
        c3 = dq - 0.5_r8*c1*(d2*(5.0_r8*d1+d2)-3.0_r8*d1*d1)
        a4(2,i,2) = qm - 0.25_r8*c1*d1*d2*(d2+3.0_r8*d1)
        ! Top edge:
        !-------------------------------------------------------
        a4(2,i,1) = d1*(2.0_r8*c1*d1**2-c3) + a4(2,i,2)
        !-------------------------------------------------------
        !        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
        !-------------------------------------------------------
        ! No over- and undershoot condition
        a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
        a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
        dc(i,1) =  0.5_r8*(a4(2,i,2) - a4(1,i,1))
      enddo
      
      ! Enforce monotonicity  within the top layer
      
      if( iv==0 ) then
        do i=i1,i2
          a4(2,i,1) = max(0.0_r8, a4(2,i,1))
          a4(2,i,2) = max(0.0_r8, a4(2,i,2))
        enddo
      elseif( iv==-1 ) then
        do i=i1,i2
          if ( a4(2,i,1)*a4(1,i,1) <= 0.0_r8 ) a4(2,i,1) = 0.0_r8
        enddo
      elseif( abs(iv)==2 ) then
        do i=i1,i2
          a4(2,i,1) = a4(1,i,1)
          a4(3,i,1) = a4(1,i,1)
        enddo
      endif
      
      ! Bottom
      ! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
        d1 = delp(i,km)
        d2 = delp(i,km1)
        qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
        dq = 2.0_r8*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
        c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.0_r8*d2*d2+d1*(d2+3.0_r8*d1)))
        c3 = dq - 2.0_r8*c1*(d2*(5._r8*d1+d2)-3.0_r8*d1*d1)
        a4(2,i,km) = qm - c1*d1*d2*(d2+3.0_r8*d1)
        ! Bottom edge:
        !-----------------------------------------------------
        a4(3,i,km) = d1*(8.0_r8*c1*d1**2-c3) + a4(2,i,km)
        !        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
        !-----------------------------------------------------
        !        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
        ! No over- and under-shoot condition
        a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
        a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
        dc(i,km) = 0.5_r8*(a4(1,i,km) - a4(2,i,km))
      enddo
      
      
      ! Enforce constraint on the "slope" at the surface
      
#ifdef BOT_MONO
      do i=i1,i2
        a4(4,i,km) = 0
        if( a4(3,i,km) * a4(1,i,km) <= 0.0_r8 ) a4(3,i,km) = 0.0_r8
        d1 = a4(1,i,km) - a4(2,i,km)
        d2 = a4(3,i,km) - a4(1,i,km)
        if ( d1*d2 < 0.0_r8 ) then
          a4(2,i,km) = a4(1,i,km)
          a4(3,i,km) = a4(1,i,km)
        else
          dq = sign(min(abs(d1),abs(d2),0.5_r8*abs(delq(i,km-1))), d1)
          a4(2,i,km) = a4(1,i,km) - dq
          a4(3,i,km) = a4(1,i,km) + dq
        endif
      enddo
#else
      if( iv==0 ) then
        do i=i1,i2
          a4(2,i,km) = max(0.0_r8,a4(2,i,km))
          a4(3,i,km) = max(0.0_r8,a4(3,i,km))
        enddo
      elseif( iv<0 ) then
        do i=i1,i2
          if( a4(1,i,km)*a4(3,i,km) <= 0.0_r8 )  a4(3,i,km) = 0.0_r8
        enddo
      endif
#endif
      
      do k=1,km1
        do i=i1,i2
          a4(3,i,k) = a4(2,i,k+1)
        enddo
      enddo
      
      !-----------------------------------------------------------
      ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
      !-----------------------------------------------------------
      ! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
        do i=i1,i2
          a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
        call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
      
      if(kord >= 7) then
        !-----------------------
        ! Huynh's 2nd constraint
        !-----------------------
        do k=2,km1
          do i=i1,i2
            ! Method#1
            !           h2(i,k) = delq(i,k) - delq(i,k-1)
            ! Method#2 - better
            h2(i,k) = 2.0_r8*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                 / ( delp(i,k)+0.5_r8*(delp(i,k-1)+delp(i,k+1)) )        &
                 * delp(i,k)**2 
            ! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
          enddo
        enddo
        
        fac = 1.5_r8           ! original quasi-monotone
        
        do k=3,km-2
          do i=i1,i2
            ! Right edges
            !        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
            !        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
            !
            pmp   = 2.0_r8*dc(i,k)
            qmp   = a4(1,i,k) + pmp
            lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
            a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                 max(a4(1,i,k), qmp, lac) )
            ! Left  edges
            !        qmp   = a4(1,i,k) - 2.0*delq(i,k)
            !        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
            !
            qmp   = a4(1,i,k) - pmp
            lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
            a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                 max(a4(1,i,k), qmp, lac))
            !-------------
            ! Recompute A6
            !-------------
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
          ! Additional constraint to ensure positivity when kord=7
          if (iv == 0 .and. kord >= 6 )                      &
               call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
        enddo
        
      else
        
        lmt = kord - 3
        lmt = max(0, lmt)
        if (iv == 0) lmt = min(2, lmt)
        
        do k=3,km-2
          if( kord /= 4) then
            do i=i1,i2
              a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
            enddo
          endif
          if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
        enddo
      endif
      
      do k=km1,km
        do i=i1,i2
          a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
        call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

    end subroutine ppm_profile

    subroutine ppm_limiters(dm, a4, itot, lmt)
      
      ! INPUT PARAMETERS:
      real(kind=r8) , intent(in):: dm(*)     !< Linear slope
      integer, intent(in) :: itot      !< Total Longitudes
      integer, intent(in) :: lmt       !< 0: Standard PPM constraint 1: Improved full monotonicity constraint
      !< (Lin) 2: Positive definite constraint 
      !< 3: do nothing (return immediately)
      ! INPUT/OUTPUT PARAMETERS:
      real(kind=r8) , intent(inout) :: a4(4,*)   !< PPM array AA <-- a4(1,i) AL <-- a4(2,i) AR <-- a4(3,i) A6 <-- a4(4,i)
      ! LOCAL VARIABLES:
      real(kind=r8)  qmp
      real(kind=r8)  da1, da2, a6da
      real(kind=r8)  fmin
      integer i
      
      ! Developer: S.-J. Lin
      
      if ( lmt == 3 ) return
      
      if(lmt == 0) then
        ! Standard PPM constraint
        do i=1,itot
          if(dm(i) == 0.0_r8) then
            a4(2,i) = a4(1,i)
            a4(3,i) = a4(1,i)
            a4(4,i) = 0.0_r8
          else
            da1  = a4(3,i) - a4(2,i)
            da2  = da1**2
            a6da = a4(4,i)*da1
            if(a6da < -da2) then
              a4(4,i) = 3.0_r8*(a4(2,i)-a4(1,i))
              a4(3,i) = a4(2,i) - a4(4,i)
            elseif(a6da > da2) then
              a4(4,i) = 3.0_r8*(a4(3,i)-a4(1,i))
              a4(2,i) = a4(3,i) - a4(4,i)
            endif
          endif
        enddo
        
      elseif (lmt == 1) then
        
        ! Improved full monotonicity constraint (Lin 2004)
        ! Note: no need to provide first guess of A6 <-- a4(4,i)
        do i=1, itot
          qmp = 2.0_r8*dm(i)
          a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
          a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
          a4(4,i) = 3.0_r8*( 2.0_r8*a4(1,i) - (a4(2,i)+a4(3,i)) )
        enddo
        
      elseif (lmt == 2) then
        
        ! Positive definite constraint
        do i=1,itot
          if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
            fmin = a4(1,i)+0.25_r8*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
            if( fmin < 0.0_r8 ) then
              if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
                a4(3,i) = a4(1,i)
                a4(2,i) = a4(1,i)
                a4(4,i) = 0.0_r8
              elseif(a4(3,i) > a4(2,i)) then
                a4(4,i) = 3.0_r8*(a4(2,i)-a4(1,i))
                a4(3,i) = a4(2,i) - a4(4,i)
              else
                a4(4,i) = 3.0_r8*(a4(3,i)-a4(1,i))
                a4(2,i) = a4(3,i) - a4(4,i)
              endif
            endif
          endif
        enddo
        
      endif
      
    end subroutine ppm_limiters
 
 
    subroutine scalar_profile(qs, a4, delp, km, i1, i2, iv, kord, qmin)
      ! Optimized vertical profile reconstruction:
      ! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
      integer, intent(in):: i1, i2
      integer, intent(in):: km      !< vertical dimension
      integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others
      integer, intent(in):: kord
      real(kind=r8), intent(in)   ::   qs(i1:i2)
      real(kind=r8), intent(in)   :: delp(i1:i2,km)     !< Layer pressure thickness
      real(kind=r8), intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
      real(kind=r8), intent(in):: qmin
      !-----------------------------------------------------------------------
      logical, dimension(i1:i2,km):: extm, ext5, ext6
      real(kind=r8)  gam(i1:i2,km)
      real(kind=r8)    q(i1:i2,km+1)
      real(kind=r8)   d4(i1:i2)
      real(kind=r8)   bet, a_bot, grat 
      real(kind=r8)   pmp_1, lac_1, pmp_2, lac_2, x0, x1
      integer i, k, im
      
      if ( iv .eq. -2 ) then
        do i=i1,i2
          gam(i,2) = 0.5_r8
          q(i,1) = 1.5_r8*a4(1,i,1)
        enddo
        do k=2,km-1
          do i=i1, i2
            grat = delp(i,k-1) / delp(i,k)
            bet =  2.0_r8 + grat + grat - gam(i,k)
            q(i,k) = (3.0_r8*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
          enddo
        enddo
        do i=i1,i2
          grat = delp(i,km-1) / delp(i,km) 
          q(i,km) = (3.0_r8*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
               (2.0_r8 + grat + grat - gam(i,km))
          q(i,km+1) = qs(i)
        enddo
        do k=km-1,1,-1
          do i=i1,i2
            q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
          enddo
        enddo
      else
        do i=i1,i2
          grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5_r8)
          q(i,1) = ( (grat+grat)*(grat+1.0_r8)*a4(1,i,1) + a4(1,i,2) ) / bet
          gam(i,1) = ( 1.0_r8 + grat*(grat+1.5_r8) ) / bet
        enddo
        
        do k=2,km
          do i=i1,i2
            d4(i) = delp(i,k-1) / delp(i,k)
            bet =  2.0_r8 + d4(i) + d4(i) - gam(i,k-1)
            q(i,k) = ( 3.0_r8*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
            gam(i,k) = d4(i) / bet
          enddo
        enddo
        
        do i=i1,i2
          a_bot = 1.0_r8 + d4(i)*(d4(i)+1.5_r8)
          q(i,km+1) = (2.0_r8*d4(i)*(d4(i)+1.0_r8)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5_r8) - a_bot*gam(i,km) )
        enddo
        
        do k=km,1,-1
          do i=i1,i2
            q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
          enddo
        enddo
      endif
      
      !----- Perfectly linear scheme --------------------------------
      if ( abs(kord) > 16 ) then
        do k=1,km
          do i=i1,i2
            a4(2,i,k) = q(i,k  )
            a4(3,i,k) = q(i,k+1)
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        enddo
        return
      endif
      !----- Perfectly linear scheme --------------------------------
      !------------------
      ! Apply constraints
      !------------------
      im = i2 - i1 + 1
      
      ! Apply *large-scale* constraints 
      do i=i1,i2
        q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
        q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
      enddo
      
      do k=2,km
        do i=i1,i2
          gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
        enddo
      enddo
      
      ! Interior:
      do k=3,km-1
        do i=i1,i2
          if ( gam(i,k-1)*gam(i,k+1)>0.0_r8 ) then
            ! Apply large-scale constraint to ALL fields if not local max/min
            q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
            q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
          else
            if ( gam(i,k-1) > 0.0_r8 ) then
              ! There exists a local max
              q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
            else
              ! There exists a local min
              q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
              if ( iv==0 ) q(i,k) = max(0.0_r8, q(i,k))
            endif
          endif
        enddo
      enddo
      
      ! Bottom:
      do i=i1,i2
        q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
        q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
      enddo
      
      do k=1,km
        do i=i1,i2
          a4(2,i,k) = q(i,k  )
          a4(3,i,k) = q(i,k+1)
        enddo
      enddo
      
      do k=1,km
        if ( k==1 .or. k==km ) then
          do i=i1,i2
            extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.0_r8
          enddo
        else
          do i=i1,i2
            extm(i,k) = gam(i,k)*gam(i,k+1) < 0.0_r8
          enddo
        endif
        if ( abs(kord) > 9 ) then
          do i=i1,i2
            x0 = 2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
            x1 = abs(a4(2,i,k)-a4(3,i,k))
            a4(4,i,k) = 3.0_r8*x0
            ext5(i,k) = abs(x0) > x1
            ext6(i,k) = abs(a4(4,i,k)) > x1
          enddo
        endif
      enddo
      
      !---------------------------
      ! Apply subgrid constraints:
      !---------------------------
      ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
      ! Top 2 and bottom 2 layers always use monotonic mapping
      
      if ( iv==0 ) then
        do i=i1,i2
          a4(2,i,1) = max(0.0_r8, a4(2,i,1))
        enddo
      elseif ( iv==-1 ) then 
        do i=i1,i2
          if ( a4(2,i,1)*a4(1,i,1) <= 0.0_r8 ) a4(2,i,1) = 0.0_r8
        enddo
      elseif ( iv==2 ) then
        do i=i1,i2
          a4(2,i,1) = a4(1,i,1)
          a4(3,i,1) = a4(1,i,1)
          a4(4,i,1) = 0.0_r8
        enddo
      endif
      
      if ( iv/=2 ) then
        do i=i1,i2
          a4(4,i,1) = 3.0_r8*(2.0_r8*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
        enddo
        call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
      endif
      
      ! k=2
      do i=i1,i2
        a4(4,i,2) = 3.0_r8*(2.0_r8*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
      enddo
      call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)
      
      !-------------------------------------
      ! Huynh's 2nd constraint for interior:
      !-------------------------------------
      do k=3,km-2
        if ( abs(kord)<9 ) then
          do i=i1,i2
            ! Left  edges
            pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
            lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
            a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                 max(a4(1,i,k), pmp_1, lac_1) )
            ! Right edges
            pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
            lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
            a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                 max(a4(1,i,k), pmp_2, lac_2) )
            
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
          
        elseif ( abs(kord)==9 ) then
          do i=i1,i2
            if ( extm(i,k) .and. extm(i,k-1) ) then
              ! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else if ( extm(i,k) .and. extm(i,k+1) ) then
              ! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else if ( extm(i,k) .and. a4(1,i,k)<qmin ) then
              ! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else
              a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              ! Check within the smooth region if subgrid profile is non-monotonic
              if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
                a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              endif
            endif
          enddo
        elseif ( abs(kord)==10 ) then
          do i=i1,i2
            if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                a4(2,i,k) = a4(1,i,k)
                a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
              endif
            elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
              endif
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        elseif ( abs(kord)==12 ) then
          do i=i1,i2
            if( extm(i,k) ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else        ! not a local extremum
              a4(4,i,k) = 6.0_r8*a4(1,i,k) - 3.0_r8*(a4(2,i,k)+a4(3,i,k))
              ! Check within the smooth region if subgrid profile is non-monotonic
              if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
                a4(4,i,k) = 6.0_r8*a4(1,i,k) - 3.0_r8*(a4(2,i,k)+a4(3,i,k))
              endif
            endif
          enddo
        elseif ( abs(kord)==13 ) then
          do i=i1,i2
            if( ext6(i,k) ) then
              if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
                ! grid-scale 2-delta-z wave detected
                a4(2,i,k) = a4(1,i,k)
                a4(3,i,k) = a4(1,i,k)
              endif
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        elseif ( abs(kord)==14 ) then
          
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
          
        elseif ( abs(kord)==15 ) then   ! Revised abs(kord)=9 scheme
          do i=i1,i2
            if ( ext5(i,k) .and. ext5(i,k-1) ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
            else if ( ext5(i,k) .and. ext5(i,k+1) ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
            else if ( ext5(i,k) .and. a4(1,i,k)<qmin ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
            elseif( ext6(i,k) ) then
              pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
              lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                   max(a4(1,i,k), pmp_1, lac_1) )
              pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
              lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                   max(a4(1,i,k), pmp_2, lac_2) )
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        elseif ( abs(kord)==16 ) then
          do i=i1,i2
            if( ext5(i,k) ) then
              if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                a4(2,i,k) = a4(1,i,k)
                a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                ! Left  edges
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                     max(a4(1,i,k), pmp_1, lac_1) )
                ! Right edges
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                     max(a4(1,i,k), pmp_2, lac_2) )
              endif
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        else      ! kord = 11, 13
          do i=i1,i2
            if ( ext5(i,k) .and. (ext5(i,k-1).or.ext5(i,k+1).or.a4(1,i,k)<qmin) ) then
              ! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else
              a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
            endif
          enddo
        endif
        
        ! Additional constraint to ensure positivity
        if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)
        
      enddo      ! k-loop
      
      !----------------------------------
      ! Bottom layer subgrid constraints:
      !----------------------------------
      if ( iv==0 ) then
        do i=i1,i2
          a4(3,i,km) = max(0.0_r8, a4(3,i,km))
        enddo
      elseif ( iv .eq. -1 ) then 
        do i=i1,i2
          if ( a4(3,i,km)*a4(1,i,km) <= 0.0_r8 )  a4(3,i,km) = 0.0_r8
        enddo
      endif
      
      do k=km-1,km
        do i=i1,i2
          a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
        if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
        if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
      enddo

    end subroutine scalar_profile

 
    subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
      ! Optimized vertical profile reconstruction:
      ! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
      integer, intent(in):: i1, i2
      integer, intent(in):: km      !< vertical dimension
      integer, intent(in):: iv      !< iv =-1: winds
      !< iv = 0: positive definite scalars
      !< iv = 1: others
      integer, intent(in):: kord
      real(kind=r8), intent(in)   ::   qs(i1:i2)
      real(kind=r8), intent(in)   :: delp(i1:i2,km)     !< layer pressure thickness
      real(kind=r8), intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
      !-----------------------------------------------------------------------
      logical, dimension(i1:i2,km):: extm, ext5, ext6
      real(kind=r8)  gam(i1:i2,km)
      real(kind=r8)    q(i1:i2,km+1)
      real(kind=r8)   d4(i1:i2)
      real(kind=r8)   bet, a_bot, grat 
      real(kind=r8)   pmp_1, lac_1, pmp_2, lac_2, x0, x1
      integer i, k, im
      
      if ( iv .eq. -2 ) then
        do i=i1,i2
          gam(i,2) = 0.5_r8
          q(i,1) = 1.5_r8*a4(1,i,1)
        enddo
        do k=2,km-1
          do i=i1, i2
            grat = delp(i,k-1) / delp(i,k)
            bet =  2.0_r8 + grat + grat - gam(i,k)
            q(i,k) = (3.0_r8*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
          enddo
        enddo
        do i=i1,i2
          grat = delp(i,km-1) / delp(i,km) 
          q(i,km) = (3.0_r8*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
               (2.0_r8 + grat + grat - gam(i,km))
          q(i,km+1) = qs(i)
        enddo
        do k=km-1,1,-1
          do i=i1,i2
            q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
          enddo
        enddo
      else
        do i=i1,i2
          grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5_r8)
          q(i,1) = ( (grat+grat)*(grat+1.0_r8)*a4(1,i,1) + a4(1,i,2) ) / bet
          gam(i,1) = ( 1.0_r8 + grat*(grat+1.5_r8) ) / bet
        enddo
        
        do k=2,km
          do i=i1,i2
            d4(i) = delp(i,k-1) / delp(i,k)
            bet =  2.0_r8 + d4(i) + d4(i) - gam(i,k-1)
            q(i,k) = ( 3.0_r8*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
            gam(i,k) = d4(i) / bet
          enddo
        enddo
        
        do i=i1,i2
          a_bot = 1.0_r8 + d4(i)*(d4(i)+1.5_r8)
          q(i,km+1) = (2.0_r8*d4(i)*(d4(i)+1.0_r8)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5_r8) - a_bot*gam(i,km) )
        enddo
        
        do k=km,1,-1
          do i=i1,i2
            q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
          enddo
        enddo
      endif
      !----- Perfectly linear scheme --------------------------------
      if ( abs(kord) > 16 ) then
        do k=1,km
          do i=i1,i2
            a4(2,i,k) = q(i,k  )
            a4(3,i,k) = q(i,k+1)
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        enddo
        return
      endif
      !----- Perfectly linear scheme --------------------------------
      
      !------------------
      ! Apply constraints
      !------------------
      im = i2 - i1 + 1
      
      ! Apply *large-scale* constraints 
      do i=i1,i2
        q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
        q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
      enddo
      
      do k=2,km
        do i=i1,i2
          gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
        enddo
      enddo
      
      ! Interior:
      do k=3,km-1
        do i=i1,i2
          if ( gam(i,k-1)*gam(i,k+1)>0.0_r8 ) then
            ! Apply large-scale constraint to ALL fields if not local max/min
            q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
            q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
          else
            if ( gam(i,k-1) > 0.0_r8 ) then
              ! There exists a local max
              q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
            else
              ! There exists a local min
              q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
              if ( iv==0 ) q(i,k) = max(0.0_r8, q(i,k))
            endif
          endif
        enddo
      enddo
      
      ! Bottom:
      do i=i1,i2
        q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
        q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
      enddo
      
      do k=1,km
        do i=i1,i2
          a4(2,i,k) = q(i,k  )
          a4(3,i,k) = q(i,k+1)
        enddo
      enddo
      
      do k=1,km
        if ( k==1 .or. k==km ) then
          do i=i1,i2
            extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.0_r8
          enddo
        else
          do i=i1,i2
            extm(i,k) = gam(i,k)*gam(i,k+1) < 0.0_r8
          enddo
        endif
        if ( abs(kord) > 9 ) then
          do i=i1,i2
            x0 = 2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
            x1 = abs(a4(2,i,k)-a4(3,i,k))
            a4(4,i,k) = 3.0_r8*x0
            ext5(i,k) = abs(x0) > x1
            ext6(i,k) = abs(a4(4,i,k)) > x1
          enddo
        endif
      enddo
      
      !---------------------------
      ! Apply subgrid constraints:
      !---------------------------
      ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
      ! Top 2 and bottom 2 layers always use monotonic mapping
      
      if ( iv==0 ) then
        do i=i1,i2
          a4(2,i,1) = max(0.0_r8, a4(2,i,1))
        enddo
      elseif ( iv==-1 ) then 
        do i=i1,i2
          if ( a4(2,i,1)*a4(1,i,1) <= 0.0_r8 ) a4(2,i,1) = 0.0_r8
        enddo
      elseif ( iv==2 ) then
        do i=i1,i2
          a4(2,i,1) = a4(1,i,1)
          a4(3,i,1) = a4(1,i,1)
          a4(4,i,1) = 0.0_r8
        enddo
      endif
      
      if ( iv/=2 ) then
        do i=i1,i2
          a4(4,i,1) = 3.0_r8*(2.0_r8*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
        enddo
        call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
      endif
      
      do i=i1,i2
        a4(4,i,2) = 3.0_r8*(2.0_r8*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
      enddo
      call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)
      
      !-------------------------------------
      ! Huynh's 2nd constraint for interior:
      !-------------------------------------
      do k=3,km-2
        if ( abs(kord)<9 ) then
          do i=i1,i2
            ! Left  edges
            pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
            lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
            a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                 max(a4(1,i,k), pmp_1, lac_1) )
            ! Right edges
            pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
            lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
            a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                 max(a4(1,i,k), pmp_2, lac_2) )
            
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
          
        elseif ( abs(kord)==9 ) then
          do i=i1,i2
            if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
              ! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
              ! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else
              a4(4,i,k) = 6.0_r8*a4(1,i,k) - 3.0_r8*(a4(2,i,k)+a4(3,i,k))
              ! Check within the smooth region if subgrid profile is non-monotonic
              if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
                a4(4,i,k) = 6.0_r8*a4(1,i,k) - 3.0_r8*(a4(2,i,k)+a4(3,i,k))
              endif
            endif
          enddo
        elseif ( abs(kord)==10 ) then
          do i=i1,i2
            if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                a4(2,i,k) = a4(1,i,k)
                a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
              endif
            elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
              endif
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        elseif ( abs(kord)==12 ) then
          do i=i1,i2
            if( extm(i,k) ) then
              ! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else        ! not a local extremum
              a4(4,i,k) = 6.0_r8*a4(1,i,k) - 3.0_r8*(a4(2,i,k)+a4(3,i,k))
              ! Check within the smooth region if subgrid profile is non-monotonic
              if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                     max(a4(1,i,k), pmp_1, lac_1) )
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                     max(a4(1,i,k), pmp_2, lac_2) )
                a4(4,i,k) = 6.0_r8*a4(1,i,k) - 3.0_r8*(a4(2,i,k)+a4(3,i,k))
              endif
            endif
          enddo
        elseif ( abs(kord)==13 ) then
          do i=i1,i2
            if( ext6(i,k) ) then
              if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
                ! grid-scale 2-delta-z wave detected
                a4(2,i,k) = a4(1,i,k)
                a4(3,i,k) = a4(1,i,k)
              endif
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        elseif ( abs(kord)==14 ) then
          
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
          
        elseif ( abs(kord)==15 ) then   ! revised kord=9 scehem
          do i=i1,i2
            if ( ext5(i,k) ) then  ! c90_mp122
              if ( ext5(i,k-1) .or. ext5(i,k+1) ) then  ! c90_mp122
                ! grid-scale 2-delta-z wave detected
                a4(2,i,k) = a4(1,i,k)
                a4(3,i,k) = a4(1,i,k)
              endif
            elseif( ext6(i,k) ) then
              ! Check within the smooth region if subgrid profile is non-monotonic
              pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
              lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                   max(a4(1,i,k), pmp_1, lac_1) )
              pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
              lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                   max(a4(1,i,k), pmp_2, lac_2) )
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        elseif ( abs(kord)==16 ) then
          do i=i1,i2
            if( ext5(i,k) ) then
              if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                a4(2,i,k) = a4(1,i,k)
                a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                ! Left  edges
                pmp_1 = a4(1,i,k) - 2.0_r8*gam(i,k+1)
                lac_1 = pmp_1 + 1.5_r8*gam(i,k+2)
                a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                     max(a4(1,i,k), pmp_1, lac_1) )
                ! Right edges
                pmp_2 = a4(1,i,k) + 2.0_r8*gam(i,k)
                lac_2 = pmp_2 - 1.5_r8*gam(i,k-1)
                a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                     max(a4(1,i,k), pmp_2, lac_2) )
              endif
            endif
          enddo
          do i=i1,i2
            a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          enddo
        else      ! kord = 11
          do i=i1,i2
            if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
              ! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.0_r8
            else
              a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
            endif
          enddo
        endif
        
        ! Additional constraint to ensure positivity
        if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)
        
      enddo      ! k-loop
      
      !----------------------------------
      ! Bottom layer subgrid constraints:
      !----------------------------------
      if ( iv==0 ) then
        do i=i1,i2
          a4(3,i,km) = max(0.0_r8, a4(3,i,km))
        enddo
      elseif ( iv .eq. -1 ) then 
        do i=i1,i2
          if ( a4(3,i,km)*a4(1,i,km) <= 0.0_r8 )  a4(3,i,km) = 0.0_r8
        enddo
      endif
      
      do k=km-1,km
        do i=i1,i2
          a4(4,i,k) = 3.0_r8*(2.0_r8*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
        if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
        if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
      enddo
      
    end subroutine cs_profile
    
    subroutine cs_limiters(im, extm, a4, iv)
      integer, intent(in) :: im
      integer, intent(in) :: iv
      logical, intent(in) :: extm(im)
      real(kind=r8) , intent(inout) :: a4(4,im)   !< PPM array
      ! LOCAL VARIABLES:
      real(kind=r8)  da1, da2, a6da
      integer i
      
      if ( iv==0 ) then
        ! Positive definite constraint
        do i=1,im
          if( a4(1,i)<=0.0_r8) then
            a4(2,i) = a4(1,i)
            a4(3,i) = a4(1,i)
            a4(4,i) = 0.0_r8
          else
            if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
              if( (a4(1,i)+0.25_r8*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0.0_r8 ) then
                ! local minimum is negative
                if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                  a4(3,i) = a4(1,i)
                  a4(2,i) = a4(1,i)
                  a4(4,i) = 0.0_r8
                elseif( a4(3,i) > a4(2,i) ) then
                  a4(4,i) = 3.0_r8*(a4(2,i)-a4(1,i))
                  a4(3,i) = a4(2,i) - a4(4,i)
                else
                  a4(4,i) = 3.0_r8*(a4(3,i)-a4(1,i))
                  a4(2,i) = a4(3,i) - a4(4,i)
                endif
              endif
            endif
          endif
        enddo
      elseif ( iv==1 ) then
        do i=1,im
          if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0.0_r8 ) then
            a4(2,i) = a4(1,i)
            a4(3,i) = a4(1,i)
            a4(4,i) = 0.0_r8
          else
            da1  = a4(3,i) - a4(2,i)
            da2  = da1**2
            a6da = a4(4,i)*da1
            if(a6da < -da2) then
              a4(4,i) = 3.0_r8*(a4(2,i)-a4(1,i))
              a4(3,i) = a4(2,i) - a4(4,i)
            elseif(a6da > da2) then
              a4(4,i) = 3.0_r8*(a4(3,i)-a4(1,i))
              a4(2,i) = a4(3,i) - a4(4,i)
            endif
          endif
        enddo
      else
        ! Standard PPM constraint
        do i=1,im
          if( extm(i) ) then
            a4(2,i) = a4(1,i)
            a4(3,i) = a4(1,i)
            a4(4,i) = 0.0_r8
          else
            da1  = a4(3,i) - a4(2,i)
            da2  = da1**2
            a6da = a4(4,i)*da1
            if(a6da < -da2) then
              a4(4,i) = 3.0_r8*(a4(2,i)-a4(1,i))
              a4(3,i) = a4(2,i) - a4(4,i)
            elseif(a6da > da2) then
              a4(4,i) = 3.0_r8*(a4(3,i)-a4(1,i))
              a4(2,i) = a4(3,i) - a4(4,i)
            endif
          endif
        enddo
      endif
    end subroutine cs_limiters
    
    
    subroutine fillz(im, km, nq, q, dp)
      integer,  intent(in):: im                !< No. of longitudes
      integer,  intent(in):: km                !< No. of levels
      integer,  intent(in):: nq                !< Total number of tracers
      real(kind=r8) , intent(in)::  dp(im,km)           !< pressure thickness
      real(kind=r8) , intent(inout) :: q(im,km,nq)      !< tracer mixing ratio
      ! LOCAL VARIABLES:
      logical:: zfix(im)
      real(kind=r8) ::  dm(km)
      integer i, k, ic, k1
      real(kind=r8)  qup, qly, dup, dq, sum0, sum1, fac
      
      do ic=1,nq
#ifdef DEV_GFS_PHYS
        ! Bottom up:
        do k=km,2,-1
          k1 = k-1
          do i=1,im
            if( q(i,k,ic) < 0.0_r8 ) then
              q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
              q(i,k ,ic) = 0.0_r8
            endif
          enddo
        enddo
        ! Top down:
        do k=1,km-1
          k1 = k+1
          do i=1,im
            if( q(i,k,ic) < 0.0_r8 ) then
              q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
              q(i,k ,ic) = 0.0_r8
            endif
          enddo
        enddo
#else
        ! Top layer
        do i=1,im
          if( q(i,1,ic) < 0.0_r8 ) then
            q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
            q(i,1,ic) = 0.0_r8
          endif
        enddo
        
        ! Interior
        zfix(:) = .false.
        do k=2,km-1
          do i=1,im
            if( q(i,k,ic) < 0.0_r8 ) then
              zfix(i) = .true.
              if ( q(i,k-1,ic) > 0.0_r8 ) then
                ! Borrow from above
                dq = min ( q(i,k-1,ic)*dp(i,k-1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k-1,ic) = q(i,k-1,ic) - dq/dp(i,k-1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
              endif
              if ( q(i,k,ic)<0.0_r8 .and. q(i,k+1,ic)>0.0_r8 ) then
                ! Borrow from below:
                dq = min ( q(i,k+1,ic)*dp(i,k+1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k+1,ic) = q(i,k+1,ic) - dq/dp(i,k+1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
              endif
            endif
          enddo
        enddo
        
        ! Bottom layer
        k = km
        do i=1,im
          if( q(i,k,ic)<0.0_r8 .and. q(i,k-1,ic)>0.0_r8) then
            zfix(i) = .true.
            ! Borrow from above
            qup =  q(i,k-1,ic)*dp(i,k-1)
            qly = -q(i,k  ,ic)*dp(i,k  )
            dup =  min(qly, qup)
            q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
            q(i,k,  ic) = q(i,k,  ic) + dup/dp(i,k  )
          endif
        enddo
        
        ! Perform final check and non-local fix if needed
        do i=1,im
          if ( zfix(i) ) then
            sum0 = 0.0_r8
            do k=2,km
              dm(k) = q(i,k,ic)*dp(i,k)
              sum0 = sum0 + dm(k)
            enddo
            
            if ( sum0 > 0.0_r8 ) then
              sum1 = 0.0_r8
              do k=2,km
                sum1 = sum1 + max(0.0_r8, dm(k))
              enddo
              fac = sum0 / sum1
              do k=2,km
                q(i,k,ic) = max(0.0_r8, fac*dm(k)/dp(i,k))
              enddo
            endif
            
          endif
        enddo
#endif
        
      enddo
    end subroutine fillz
  end module fv_mapz
