!=======================================================================================================
! extracted from TIEGCM model
! see section 6.3 in https://www.hao.ucar.edu/modeling/tgcm/doc/description/model_description.pdf
!=======================================================================================================
module tei_mod
  !
  ! Calculate electron and ion temperatures.
  !
  use shr_kind_mod,  only : r8 => shr_kind_r8
  use shr_const_mod, only : pi => shr_const_pi           ! Boltzmann constant and pi
  use mo_constants,  only : gask => rgas_cgs
  use mo_constants,  only : boltz=>boltz_cgs
  use physconst,     only : avogad ! molecules/kmole

  implicit none
  private
  public :: settei

  ! g = 8.7 m/s at 400 km altitude
  real(r8), parameter :: grav   = 870._r8      ! (cm/s^2)
  real(r8), parameter :: dipmin = 0.24_r8      ! minimum mag dip angle (2.5 deg horizontal res)
  real(r8), parameter :: rtd = 180._r8/pi
  real(r8), parameter :: evergs = 1.602e-12_r8 ! 1 eV = 1.602e-12 ergs

  real(r8), parameter :: rmassinv_o2 = 1._r8/32._r8
  real(r8), parameter :: rmassinv_o1 = 1._r8/16._r8
  real(r8), parameter :: rmassinv_n2 = 1._r8/28._r8
  real(r8), parameter :: avo = avogad*1.e-3_r8 ! molecules/mole
contains

  ! 
  subroutine settei(tn,o2,o1,n2,ne,te,ti,op,o2p,nop,barm,qji_ti, &
       f107, chi, qtot, qteaur, rlatm, dipmag, pmid, pint, lev0, lev1,ncol, &
       te_out, ti_out, qtotal )

    use trsolv_mod, only : trsolv
    !
    ! Args:
    integer,intent(in) :: lev0,lev1,ncol
    real(r8),dimension(lev0:lev1,ncol),intent(in) :: &
         tn,    &   ! neutral temperature (deg K) 
         o2,    &   ! molecular oxygen (mmr) 
         o1,    &   ! atomic oxygen (mmr)
         n2,    &   ! molecular nitrogen (mmr)
         ne,    &   ! electron density (cm3)
         te,    &   ! electron temperature (from previous time step) (K)
         ti,    &   ! ion temperature (from previous time step) (K)
         op,    &   ! O+ number dens (/cm^3) 
         o2p,   &   ! O2+ number dens (/cm^3) 
         nop,   &   ! NO+ number dens (/cm^3) 
         barm,  &   ! mean molecular weight (g/mole)
         qji_ti     ! joule heating from qjoule_ti (used ui,vi)  ergs/s/g

    real(r8), intent(in) :: f107
    real(r8), intent(in) :: chi(ncol)              ! solar zenith angle (radians)
    real(r8), intent(in) :: qtot(lev0:lev1,ncol)   ! total ionization rate  (s-1 cm-3)
    real(r8), intent(in) :: qteaur(ncol)           ! (units ?? )
    real(r8), intent(in) :: rlatm(ncol)            ! geo mag latitude (radians)
    real(r8), intent(in) :: dipmag(ncol)           ! geo mag dip angle (radians)
    real(r8), intent(in) :: pmid(lev0:lev1,ncol)   ! mid-level press (dyne/cm^2)
    real(r8), intent(in) :: pint(lev0:lev1+1,ncol) ! interface press (dyne/cm^2)

    !
    ! Output args:
    real(r8),dimension(lev0:lev1,ncol),intent(out) :: &
         te_out, & ! output electron temperature (deg K) 
         ti_out    ! output ion temperature (deg K)
    real(r8), intent(out) :: qtotal(lev0:lev1,ncol)! (ergs/sec/gm)
    
    !
    ! Local:
    integer :: i,k,n,ier
    integer :: nk,nkm1
    real(r8),dimension(lev0:lev1,ncol) :: &
         xnmbar,    & ! p0*e(-z)*barm/kT  (minpoints or interfaces)
         p_coef,    & ! coefficient for trisolv     (s1)
         q_coef,    & ! coefficient for trisolv     (s2)
         r_coef,    & ! coefficient for trisolv     (s3)
         rhs          ! right-hand-side for trisolv (s4)
    real(r8),dimension(lev0:lev1) :: &
         te_int,    & ! electron temperature (interfaces)
         tn_int,    & ! neutral temperature (interfaces)
         o2n,       & ! O2 number density (interfaces)
         o1n,       & ! O1 number density (interfaces)
         n2n,       & ! N2 number density (interfaces)
         root_te,   & ! sqrt(te)
         root_tn,   & ! sqrt(tn)
         root_ne,   & ! sqrt(ne)
         tek0,      & ! ke/te**2.5 (s15)
         h_mid,h_int,&
         barm_int,&
         fki,       & ! work array
         qe,        & ! source term                 (s10)
         q_eni,     & ! heating from electron/neutral and electron/ion collisions
         coll_en2v, & ! electron/N2vib collision    (s9)
                                !
                                ! Cooling rates (heat loss):
         loss_en2v, & ! electron/N2vib loss term    (s10)
         loss_en2,  & ! electron/N2 loss
         loss_eo2,  & ! electron/O2 loss
         loss_eo1d, & ! electron/O(1d) loss
         loss_eo1,  & ! electron/O loss
         loss_xen,  & ! L0*(E,N) (s8)
         loss_en      ! electrons/neutrals loss     (s11)
    real(r8),dimension(lev0:lev1,ncol) :: &
         loss_ei,   & ! electron/ion loss           (s10)
         loss_in      ! ion/neutral loss            (s9)
    real(r8),parameter :: &
         fpolar = -3.0e+9_r8, & ! polar te flux
         del    = 1.e-6_r8  , & ! some small value
         root2 = sqrt(2._r8), &
   !
   ! Correction factors for neutral heating due to L(E,O1D)
         alam = 0.0069_r8   , &
         ad   = 0.0091_r8   , &
         sd   = 2.3e-11_r8
    real(r8) ::  &
         f107te   ! solar flux
    !
    ! a,fed,fen,fe,sindipmag have a z dimension only for diagnostic plotting:
    real(r8) :: &
         a,fed,fen,& ! day/night
         fe,       & ! heat flux at upper boundary
         sindipmag   ! sin(dipmag)
    !
    real(r8) :: dellnp(lev0:lev1) ! delta of log interface pressures

    qtotal(:,:) = 0._r8

    col_loop: do i=1,ncol
       do k=lev0,lev1
          dellnp(k) = abs(log(pint(k,i)) - log(pint(k+1,i)))
       end do
       !
       f107te = f107
       if (f107te > 235._r8) f107te = 235._r8
       nk = lev1-lev0+1
       nkm1 = nk-1
       !
       if (abs(rlatm(i))-pi/4.5_r8 >= 0._r8) then
          a = 1._r8
       else
          a = .5_r8*(1._r8+sin(pi*(abs(rlatm(i))-pi/9._r8)/(pi/4.5_r8)))
       endif
       !
       ! Increased heat flux for TE fom protonosphere.
       fed = ( -5.0e+7_r8*f107te*a-4.0e+7_r8*f107te)*1.2_r8
       fen = fed/2._r8
       fed = fed+qteaur(i)     ! t4
       fen = fen+qteaur(i)     ! t5
       if (chi(i)-.5_r8*pi >= 0._r8) then  ! chi==t2
          fe = fen                ! t1
       else
          fe = fed
       endif
       if ((chi(i)*rtd-80._r8)*(chi(i)*rtd-100._r8)>=0._r8) then
          fe = fe*evergs
       else
          fe = (.5_r8*(fed+fen)+.5_r8*(fed-fen)* &
               cos(pi*(chi(i)*rtd-80._r8)/20._r8))*evergs
       endif
       !
       ! Add fpolar if magnetic latitude >= 60 degrees:
       if (abs(rlatm(i)-pi/3._r8)>=0._r8) then
          fe = fe+fpolar*evergs
       end if

       !
       ! te,o2,o,n2,tn at interfaces: 
       do k=lev0+1,lev1-1
          te_int(k) = .5_r8*(te(k,i)+te(k-1,i))
          o2n(k)    = .5_r8*(o2(k,i)+o2(k-1,i))
          o1n(k)    = .5_r8*(o1(k,i)+o1(k-1,i))
          n2n(k)    = .5_r8*(n2(k,i)+n2(k-1,i))
          tn_int(k) = .5_r8*(tn(k,i)+tn(k-1,i))
          barm_int(k) = .5_r8*(barm(k,i)+barm(k-1,i))
       enddo ! k=lev0+1,lev1-2
       !
       ! Bottom:
       te_int(lev0) = 1.5_r8*te(lev0,i)-.5_r8*te(lev0+1,i)
       o2n(lev0)    = 1.5_r8*o2(lev0,i)-.5_r8*o2(lev0+1,i)
       o1n(lev0)    = 1.5_r8*o1(lev0,i)-.5_r8*o1(lev0+1,i)
       n2n(lev0)    = 1.5_r8*n2(lev0,i)-.5_r8*n2(lev0+1,i)
       tn_int(lev0) = 1.5_r8*tn(lev0,i)-.5_r8*tn(lev0+1,i)
       barm_int(lev0) = 1.5_r8*barm(lev0,i)-.5_r8*barm(lev0+1,i)
       !
       ! Top:
       te_int(lev1) = 1.5_r8*te(lev1-1,i)-.5_r8*te(lev1-2,i)
       o2n(lev1)    = 1.5_r8*o2(lev1-1,i)-.5_r8*o2(lev1-2,i)
       o1n(lev1)    = 1.5_r8*o1(lev1-1,i)-.5_r8*o1(lev1-2,i)
       n2n(lev1)    = 1.5_r8*n2(lev1-1,i)-.5_r8*n2(lev1-2,i)
       tn_int(lev1) = 1.5_r8*tn(lev1-1,i)-.5_r8*tn(lev1-2,i)
       barm_int(lev1) = 1.5_r8*barm(lev1-1,i)-.5_r8*barm(lev1-2,i)
       !
       ! N2:
       do k=lev0,lev1
          if (n2n(k) < 0._r8) n2n(k) = 0._r8
       enddo ! k=lev0,lev1

       !
       ! Convert o2,o,n2 to number density (interfaces):

       do k=lev0,lev1
          ! xnmbar at interfaces:
          xnmbar(k,i) = pint(k,i)*barm_int(k)/(boltz*tn_int(k)) ! s8
          o2n(k) = xnmbar(k,i)*o2n(k)*rmassinv_o2       ! s13
          o1n(k) = xnmbar(k,i)*o1n(k)*rmassinv_o1       ! s12
          n2n(k) = xnmbar(k,i)*n2n(k)*rmassinv_n2       ! s11
          root_te(k) = sqrt(te_int(k))
          !
          tek0(k) = 7.5e5_r8/ &
               (1._r8+3.22e4_r8*te_int(k)**2/ &
               ne(k,i)*(root_te(k)* &
               (2.82e-17_r8 - 3.41e-21_r8   * te_int (k))*n2n(k)+ &
               (2.20e-16_r8 + 7.92e-18_r8   * root_te(k))*o2n(k)+ &
               1.10e-16_r8 * (1._r8+5.7e-4_r8 * te_int (k))*o1n(k)))*evergs

       enddo ! k=lev0,lev1

       do k=lev0,lev1-1
          h_mid(k) = gask*tn(k,i)/(barm(k,i)*grav) ! s7
       enddo ! k=lev0,lev1-1
       do k=lev0,lev1
          h_int(k) = gask*tn_int(k)/(barm_int(k)*grav)              ! s6
       enddo ! k=lev0,lev1

       if (abs(dipmag(i)) >= dipmin) then
          sindipmag = (sin(dipmag(i)))**2 ! t2,s2
       else
          sindipmag = (sin(dipmin))**2
       endif
       if (sindipmag < .10_r8) sindipmag = .10_r8
       !
       ! Start coefficients and rhs for trsolv:
       do k=lev0,lev1-1
          p_coef(k,i) = 2._r8/7._r8*sindipmag/(h_mid(k)*dellnp(k)**2) ! s1
          r_coef(k,i) = p_coef(k,i)*tek0(k+1)/h_int(k+1)  ! s3
          p_coef(k,i) = p_coef(k,i)*tek0(k  )/h_int(k  )  ! s1
          q_coef(k,i) = -(p_coef(k,i)+r_coef(k,i))            ! s2
          rhs(k,i) = 0._r8                                       ! s4
       enddo ! k=lev0,lev1-1
       !
       ! Bottom boundary:
       q_coef(lev0,i) = q_coef(lev0,i)-p_coef(lev0,i)
       rhs(lev0,i) = rhs(lev0,i)-2._r8*p_coef(lev0,i)*tn_int(lev0)**3.5_r8
       p_coef(lev0,i) = 0._r8
       !
       ! Upper boundary:
       q_coef(lev1-1,i) = q_coef(lev1-1,i)+r_coef(lev1-1,i)
       rhs(lev1-1,i) = rhs(lev1-1,i)+r_coef(lev1-1,i)*dellnp(lev1-1)*3.5_r8* &
            h_int(lev1)*fe/tek0(lev1)
       r_coef(lev1-1,i) = 0._r8

       !
       !
       ! Set Ne (midpoints "(K+1/2)"):
       !
       do k=lev0,lev1-1
          root_ne(k) = ne(k,i)*ne(k+1,i)
          if (root_ne(k) < 1.e4_r8) root_ne(k) = 1.e4_r8
          root_ne(k) = sqrt(root_ne(k))
       enddo ! k=lev0,lev1-1

       !
       ! Set up o2,o,n2 number densities at midpoints:
       !

       do k=lev0,lev1-1
          xnmbar(k,i) = pmid(k,i)*barm(k,i)/(boltz*tn(k,i))
          o2n(k) = xnmbar(k,i)*o2(k,i)*rmassinv_o2  ! s14
          o1n(k) = xnmbar(k,i)*o1(k,i)*rmassinv_o1  ! s13
          n2n(k) = n2(k,i)
          if (n2n(k) < 0._r8) n2n(k) = 0._r8
          n2n(k) = xnmbar(k,i)*n2n(k)*rmassinv_n2 ! s12
          !
          ! Calculate source term qe (s10)
          ! Comment from earlier version (maybe the *1.0 below was once *2.0):
          !   "Correction facor of 2 increase in TE heating rate"
          !
          qe(k) = log(root_ne(k)/(o2n(k)+n2n(k)+0.1_r8*o1n(k)))
          qe(k) = exp(-((((0.001996_r8*qe(k)+0.08034_r8)*qe(k)+1.166_r8)* &
               qe(k)+6.941_r8)*qe(k)+12.75_r8))*1.0_r8
          !
          ! Subtract qe from right-hand-side:
          rhs(k,i) = rhs(k,i)-qe(k)*qtot(k,i)*evergs
          root_te(k) = sqrt(te(k,i))
          !
          ! Electron/N2 collision A(E,N2,VIB) (s9):
          !
          if (te(k,i) >= 1000._r8) then
             coll_en2v(k) = 2.e-7_r8*exp(-4605.2_r8/te(k,i))
          else
             coll_en2v(k) = 5.71e-8_r8*exp(-3352.6_r8/te(k,i))
          endif
          if (te(k,i) > 2000._r8) then
             coll_en2v(k) = 2.53e-6_r8*root_te(k)*exp(-17620._r8/te(k,i))
          end if
          !
          ! Loss due to electron/n2 collision L0(E,N2,VIB)/(NE*N(N2)) (s10)
          !
          loss_en2v(k) = 3200._r8*(1._r8/te(k,i)-1._r8/tn(k,i))
          loss_en2v(k) = sign(abs(loss_en2v(k))+del,loss_en2v(k)) ! avoid divide by zero in the next line when te = tn
                                                                  ! by adding a small value to loss_en2v
                                                                  ! this assumes te >= tn 
          loss_en2v(k) = -3200._r8/(te(k,i)*tn(k,i))*(1._r8-exp(loss_en2v(k)))/loss_en2v(k)
          loss_en2v(k) = 1.3e-4_r8*loss_en2v(k)*coll_en2v(k)
       enddo ! k=lev0,lev1-1

       !
       ! Calculate and sum cooling rates (heat loss) due to interactions between
       ! electrons/neutrals, electrons/ions, ions/neutrals
       !
       do k=lev0,lev1-1
          !
          ! Electron/N2 loss rate:
          ! loss_en2 = (L0(E,N2)+L0(E,N2,ROT)+L0(E,N2,VIB))/NE (s11)
          !
          loss_en2(k) = n2n(k)*(1.77E-19_r8*(1._r8-1.21E-4_r8*te(k,i))* te(k,i) + 2.9e-14_r8/root_te(k) + loss_en2v(k))
          !
          ! Start total of electron/neutral loss rate (s11):
          !
          loss_en(k) = loss_en2(k)
          !
          ! Electron/O2 loss rates: (L0(E,O2)+L0(E,O2,ROT)+L0(E,O2,VIB)/NE
          !
          loss_eo2(k) = o2n(k)*(1.21e-18_r8*(1._r8+3.6e-2_r8*root_te(k))* &
               root_te(k)+6.9e-14_r8/root_te(k)+3.125e-21_r8*te(k,i)**2)
          loss_en(k) = loss_en(k)+loss_eo2(k)
          !
          ! Electron/O(1d) loss rates: L0(E,O,1D)/(NE*N(O))
          !
          loss_eo1d(k) = 22713._r8*(1._r8/te(k,i)-1._r8/tn(k,i))
          loss_eo1d(k) = sign(abs(loss_eo1d(k))+del,loss_eo1d(k))
          loss_eo1d(k) = 22713._r8/(te(k,i)*tn(k,i))*(1._r8-exp(loss_eo1d(k)))/loss_eo1d(k)
          !
          ! loss_eo1d function often fails here with bad argument to exp()
          ! due to high te and/or high loss_eo1d from above.
          !         loss_eo1d(k) = 1.57e-12*exp((2.4e4+0.3*(te(k)-1500.)-
          !                        1.947e-5*(te(k)-1500.)*(te(k)-4000.))*(te(k)-3000.)/
          !                       (3000.*te(k)))*loss_eo1d(k)
          loss_eo1d(k) = 0._r8
          !
          ! Electron/O1 loss rates: (L0(E,O)+L0(E,O,F))/NE
          !
          loss_eo1(k) = o1n(k)*(7.9e-19_r8*(1._r8+5.7e-4_r8*te(k,i))* &
               root_te(k)+3.4e-12_r8*(1._r8-7.e-5_r8*te(k,i))/tn(k,i)* &
               (150._r8/te(k,i)+0.4_r8))

          loss_en(k) = loss_en(k)+loss_eo1(k)
          !
          ! loss_xen = L0*(E,N) (s8)
          !
          loss_xen(k) = (loss_en(k)+o1n(k)*(1._r8-alam/(ad+sd* &
               n2n(k)))*loss_eo1d(k))*root_ne(k)*evergs
          !
          ! Complete total electron/neutral loss rate L0(E,N) (s11):
          !
          loss_en(k) = (loss_en(k)+o1n(k)*loss_eo1d(k))* &
               root_ne(k)*evergs
          !
          ! Calculate L0(E) = L(E)/(TE-TI), where L(E) is loss due to
          ! interactions between electrons and ions.
          !
          loss_ei(k,i) = 3.2e-8_r8*root_ne(k)/(root_te(k)*te(k,i))* &
               15._r8*(op(k,i)+0.5_r8*o2p(k,i)+0.53_r8*nop(k,i))*evergs

          root_tn(k) = sqrt(tn(k,i))
          ! 
          ! loss_in = ion/neutral cooling = L0(I,N) =L(I,N)/(TI-TN)
          !
          loss_in(k,i) = ((6.6e-14_r8*n2n(k)+5.8e-14_r8*o2n(k)+0.21e-14_r8* &
               o1n(k)*root2*root_tn(k))*op(k,i)+(5.45e-14_r8*o2n(k)+ &
               5.9e-14_r8*n2n(k)+4.5e-14_r8*o1n(k))*nop(k,i)+(5.8e-14_r8* &
               n2n(k)+4.4e-14_r8*o1n(k)+0.14e-14_r8*o2n(k)*root_tn(k))* &
               o2p(k,i))*evergs
          !
          ! Complete tridiagonal matrix coefficients and rhs:
          !
          ! q_coef = q_coef-(L0(E,N)+L0(E,I))/TE**2.5 = Q
          !
          q_coef(k,i) = q_coef(k,i)-(loss_en(k)+loss_ei(k,i))/te(k,i)**2.5_r8
          !          
          ! rhs = rhs-L0(E,N)*TN-L0(E)*TI
          !
          rhs(k,i) = rhs(k,i)-loss_en(k)*tn(k,i)-loss_ei(k,i)*ti(k,i)

       enddo ! k=lev0,lev1-1

       ! Calculate heating due to electron/neutral and electron/ion collisions
       ! (ergs/sec/gm):
       !
       do k=lev0,lev1-1
          if (te(k,i)-ti(k,i) >= 0._r8) then
             q_eni(k)=loss_ei(k,i)*(te(k,i)-ti(k,i))
          else
             q_eni(k) = 0._r8
          endif
          q_eni(k) = (loss_xen(k)*(te(k,i)-tn(k,i))+q_eni(k)) * avo/xnmbar(k,i)
       enddo ! k=lev0,lev1-1

       !
       ! Add collisional heating to Q for use in thermodynamic equation.
       do k=lev0,lev1-2
          qtotal(k+1,i) = qtotal(k+1,i)+.5_r8*(q_eni(k)+q_eni(k+1))
       enddo ! k=lev0,lev1-2
       !
       ! Upper and lower boundaries:
       qtotal(lev0,i) = qtotal(lev0,i)+1.5_r8*q_eni(lev0)-0.5_r8*q_eni(lev0+1)
       qtotal(lev1,i) = qtotal(lev1,i)+1.5_r8*q_eni(lev1-1)-0.5_r8*q_eni(lev1-2)
    end do col_loop

    !
    ! Solve tridiagonal system:
    !
    call trsolv(p_coef,q_coef,r_coef,rhs,te_out,  lev0,lev1, lev0,lev1-1, 1,ncol )

    col_loop2: do i=1,ncol
       !
       ! Te = Te**(2./7.):
       do k=lev0,lev1-1
          te_out(k,i) = te_out(k,i)**(2._r8/7._r8)
       enddo
       !
       ! 10/21/03 btf: make this check after te*(2/7), rather than before.
       !
       ! Te must be >= Tn:
       do k=lev0,lev1-1
          if (te_out(k,i) < tn(k,i)) te_out(k,i) = tn(k,i)
       enddo
       !
       ! 1/9/08 btf: put spval in top level of te:
       !    te_out(lev1) = spval
       te_out(lev1,i) = te_out(lev1-1,i)
       !
       ! Set ion temperature output. Use joule heating qji_ti from sub 
       ! qjoule_ti (see qjoule.F). lev1 not calculated.
       !
       do k=lev0,lev1-1
          ti_out(k,i) = (qji_ti(k,i)*(xnmbar(k,i)/avo)+ &
               loss_ei(k,i)*te_out(k,i)+loss_in(k,i)*tn(k,i))/&
               (loss_ei(k,i)+loss_in(k,i))
          !
          ! ti must be at least as large as tn:
          if (ti_out(k,i) < tn(k,i)) ti_out(k,i) = tn(k,i)
       enddo
       !
       ! 1/9/08 btf: put spval in top level of ti:
       !    ti_out(lev1) = spval
       ti_out(lev1,i) = ti_out(lev1-1,i)
    end do col_loop2

  end subroutine settei

end module tei_mod

