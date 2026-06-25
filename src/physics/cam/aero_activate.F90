module aero_activate

! Portable (CCPP-ready) Abdul-Razzak & Ghan aerosol activation kernel
! (activate_aerosol), extracted from the CAM ndrop module.  The derived
! constants aten and sqrt(pi) are computed once by aero_activate_init from host
! physical constants; activate_aerosol receives the remaining host physical
! constants as arguments and returns errmsg/errflg instead of aborting.  The
! polymorphic aerosol_properties abstraction is deliberately host-portable.
! CAM interface / callers: ndrop, aero_convproc.

use shr_kind_mod,     only: r8 => shr_kind_r8
use wv_saturation,    only: qsat
use shr_spfn_mod,     only: erf => shr_spfn_erf

use aerosol_properties_mod, only: aerosol_properties

implicit none
private

public :: aero_activate_init, activate_aerosol

! mathematical constants
real(r8), parameter :: zero     = 0._r8
real(r8), parameter :: third    = 1._r8/3._r8
real(r8), parameter :: twothird = 2._r8*third
real(r8), parameter :: sixth    = 1._r8/6._r8
real(r8), parameter :: sq2      = sqrt(2._r8)
real(r8), parameter :: tmelt    = 273._r8

! derived constants (set by aero_activate_init)
real(r8) :: aten
real(r8) :: sqpi

!===============================================================================
contains
!===============================================================================

subroutine aero_activate_init(mwh2o, r_universal, rhoh2o, pi)

   ! Compute the derived activation constants from host physical constants.
   ! surften (surface tension of water) is a fixed property of the activation
   ! parameterization, so it is set here rather than threaded from the host.

   real(r8), intent(in) :: mwh2o        ! molecular weight of water (kg/kmol)
   real(r8), intent(in) :: r_universal  ! universal gas constant (J/K/kmol)
   real(r8), intent(in) :: rhoh2o       ! density of liquid water (kg/m3)
   real(r8), intent(in) :: pi           ! pi

   real(r8) :: surften                  ! surface tension of water w/respect to air (N/m)

   surften = 0.076_r8
   aten = 2._r8*mwh2o*surften/(r_universal*tmelt*rhoh2o)
   sqpi = sqrt(pi)

end subroutine aero_activate_init

!===============================================================================

subroutine activate_aerosol(wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,  &
   na, nbins, volume, hygro, aero_props, &
   fn, fm, fluxn, fluxm, flux_fullact, &
   pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, errmsg, errflg, &
   smax_prescribed, in_cloud_in, smax_f)

   !      calculates number, surface, and mass fraction of aerosols activated as CCN
   !      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
   !      assumes an internal mixture within each of up to nbin multiple aerosol bins
   !      a gaussiam spectrum of updrafts can be treated.

   !      mks units

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   !      input

   real(r8), intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: sigw          ! subgrid standard deviation of vertical vel (m/s)
   real(r8), intent(in) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
   real(r8), intent(in) :: wminf         ! minimum updraft velocity for integration (m/s)
   real(r8), intent(in) :: wmaxf         ! maximum updraft velocity for integration (m/s)
   real(r8), intent(in) :: tair          ! air temperature (K)
   real(r8), intent(in) :: rhoair        ! air density (kg/m3)
   real(r8), intent(in) :: na(:)      ! aerosol number concentration (/m3)
   integer,  intent(in) :: nbins      ! number of aerosol bins
   real(r8), intent(in) :: volume(:)  ! aerosol volume concentration (m3/m3)
   real(r8), intent(in) :: hygro(:)   ! hygroscopicity of aerosol mode

   class(aerosol_properties), intent(in) :: aero_props

   !      output

   real(r8), intent(out) :: fn(:)      ! number fraction of aerosols activated
   real(r8), intent(out) :: fm(:)      ! mass fraction of aerosols activated
   real(r8), intent(out) :: fluxn(:)   ! flux of activated aerosol number fraction into cloud (cm/s)
   real(r8), intent(out) :: fluxm(:)   ! flux of activated aerosol mass fraction into cloud (cm/s)
   real(r8), intent(out) :: flux_fullact   ! flux of activated aerosol fraction assuming 100% activation (cm/s)
   !    rce-comment
   !    used for consistency check -- this should match (ekd(k)*zs(k))
   !    also, fluxm/flux_fullact gives fraction of aerosol mass flux
   !       that is activated

   !      host physical constants (from physconst) + error handling
   real(r8), intent(in) :: pi           ! pi
   real(r8), intent(in) :: rhoh2o       ! density of liquid water (kg/m3)
   real(r8), intent(in) :: rh2o         ! water vapor gas constant (J/K/kg)
   real(r8), intent(in) :: gravit       ! gravitational acceleration (m/s2)
   real(r8), intent(in) :: latvap       ! latent heat of vaporization (J/kg)
   real(r8), intent(in) :: cpair        ! specific heat of dry air (J/K/kg)
   real(r8), intent(in) :: rair         ! dry air gas constant (J/K/kg)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   !      optional
   real(r8), optional, intent(in) :: smax_prescribed  ! prescribed max. supersaturation for secondary activation
   logical,  optional, intent(in) :: in_cloud_in      ! switch to modify calculations when above cloud base
   real(r8), optional, intent(in) :: smax_f           ! droplet and rain size distr factor in the smax calculation
                                                      ! used when in_cloud=.true.

   !      local

   integer, parameter:: nx=200
   real(r8) integ,integf
   real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure (Pa)
   real(r8) pres ! pressure (Pa)
   real(r8) diff0,conduct0
   real(r8) es ! saturation vapor pressure
   real(r8) qs ! water vapor saturation mixing ratio
   real(r8) dqsdt ! change in qs with temperature
   real(r8) g ! thermodynamic function (m2/s)
   real(r8) zeta(nbins), eta(nbins)
   real(r8) alpha
   real(r8) gamma
   real(r8) beta
   real(r8) sqrtg
   real(r8) :: amcube(nbins) ! cube of dry bin radius (m)
   real(r8) smc(nbins) ! critical supersaturation for number bin radius
   real(r8) sumflx_fullact
   real(r8) sumflxn(nbins)
   real(r8) sumflxm(nbins)
   real(r8) sumfn(nbins)
   real(r8) sumfm(nbins)
   real(r8) fnold(nbins)   ! number fraction activated
   real(r8) fmold(nbins)   ! mass fraction activated
   real(r8) wold,gold
   real(r8) wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb
   real(r8) dfmin,dfmax,fnew,fold,fnmin,fnbar,fmbar
   real(r8) alw,sqrtalw
   real(r8) smax
   real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
   real(r8) etafactor1,etafactor2(nbins),etafactor2max
   real(r8) grow
   character(len=*), parameter :: subname='activate_aerosol'

   logical :: in_cloud
   integer m,n
   !      numerical integration parameters
   real(r8), parameter :: eps=0.3_r8,fmax=0.99_r8,sds=3._r8

   real(r8), parameter :: namin=1.e6_r8   ! minimum aerosol number concentration (/m3)

   errmsg = ''
   errflg = 0

   if (present(in_cloud_in)) then
      if (.not. present(smax_f)) then
         errmsg = subname//' error: smax_f must be supplied when in_cloud is used'
         errflg = 1
         return
      end if
      in_cloud = in_cloud_in
   else
      in_cloud = .false.
   end if

   fn(:)=0._r8
   fm(:)=0._r8
   fluxn(:)=0._r8
   fluxm(:)=0._r8
   flux_fullact=0._r8

   if(nbins.eq.1.and.na(1).lt.1.e-20_r8)return

   if(sigw.le.1.e-5_r8.and.wbar.le.0._r8)return

   if ( present( smax_prescribed ) ) then
      if (smax_prescribed <= 0.0_r8) return
   end if

   pres=rair*rhoair*tair
   diff0=0.211e-4_r8*(p0/pres)*(tair/tmelt)**1.94_r8
   conduct0=(5.69_r8+0.017_r8*(tair-tmelt))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg
   call qsat(tair, pres, es, qs)
   dqsdt=latvap/(rh2o*tair*tair)*qs
   alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1._r8/(rair*tair))
   gamma=(1.0_r8+latvap/cpair*dqsdt)/(rhoair*qs)
   etafactor2max=1.e10_r8/(alpha*wmaxf)**1.5_r8 ! this should make eta big if na is very small.

   grow  = 1._r8/(rhoh2o/(diff0*rhoair*qs)  &
           + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._r8))
   sqrtg = sqrt(grow)
   beta  = 2._r8*pi*rhoh2o*grow*gamma

   do m=1,nbins

      if(volume(m).gt.1.e-39_r8.and.na(m).gt.1.e-39_r8)then
         ! number mode radius (m)
         amcube(m)=aero_props%amcube(m, volume(m),na(m))
         ! growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
         ! should depend on mean radius of mode to account for gas kinetic effects
         ! see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
         ! for approriate size to use for effective diffusivity.
         etafactor2(m)=1._r8/(na(m)*beta*sqrtg)
         if(hygro(m).gt.1.e-10_r8)then
            smc(m)=2._r8*aten*sqrt(aten/(27._r8*hygro(m)*amcube(m))) ! only if variable size dist
         else
            smc(m)=100._r8
         endif
      else
         smc(m)=1._r8
         etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
      endif

   enddo

   if(sigw.gt.1.e-5_r8)then ! spectrum of updrafts

      wmax=min(wmaxf,wbar+sds*sigw)
      wmin=max(wminf,-wdiab)
      wmin=max(wmin,wbar-sds*sigw)
      w=wmin
      dwmax=eps*sigw
      dw=dwmax
      dfmax=0.2_r8
      dfmin=0.1_r8
      if (wmax <= w) return
      do m=1,nbins
         sumflxn(m)=0._r8
         sumfn(m)=0._r8
         fnold(m)=0._r8
         sumflxm(m)=0._r8
         sumfm(m)=0._r8
         fmold(m)=0._r8
      enddo
      sumflx_fullact=0._r8

      fold=0._r8
      wold=0._r8
      gold=0._r8

      dwmin = min( dwmax, 0.01_r8 )
      do n = 1, nx

100      wnuc=w+wdiab
         !           write(iulog,*)'wnuc=',wnuc
         alw=alpha*wnuc
         sqrtalw=sqrt(alw)
         etafactor1=alw*sqrtalw

         do m=1,nbins
            eta(m)=etafactor1*etafactor2(m)
            zeta(m)=twothird*sqrtalw*aten/sqrtg
         enddo

         if ( present( smax_prescribed ) ) then
            smax = smax_prescribed
         else
            smax = aero_props%maxsat(zeta,eta,smc)
         endif

         call aero_props%actfracs( nbins, smc(nbins), smax, fnew, fm(nbins) )

         dwnew = dw
         if(fnew-fold.gt.dfmax.and.n.gt.1)then
            !              reduce updraft increment for greater accuracy in integration
            if (dw .gt. 1.01_r8*dwmin) then
               dw=0.7_r8*dw
               dw=max(dw,dwmin)
               w=wold+dw
               go to 100
            else
               dwnew = dwmin
            endif
         endif

         if(fnew-fold.lt.dfmin)then
            !              increase updraft increment to accelerate integration
            dwnew=min(1.5_r8*dw,dwmax)
         endif
         fold=fnew

         z=(w-wbar)/(sigw*sq2)
         g=exp(-z*z)
         fnmin=1._r8

         do m=1,nbins
            !              modal
            call aero_props%actfracs( m, smc(m), smax, fn(m), fm(m) )
            fnmin=min(fn(m),fnmin)
            !               integration is second order accurate
            !               assumes linear variation of f*g with w
            fnbar=(fn(m)*g+fnold(m)*gold)
            fmbar=(fm(m)*g+fmold(m)*gold)
            wb=(w+wold)
            if(w.gt.0._r8)then
               sumflxn(m)=sumflxn(m)+sixth*(wb*fnbar           &
                  +(fn(m)*g*w+fnold(m)*gold*wold))*dw
               sumflxm(m)=sumflxm(m)+sixth*(wb*fmbar           &
                  +(fm(m)*g*w+fmold(m)*gold*wold))*dw
            endif
            sumfn(m)=sumfn(m)+0.5_r8*fnbar*dw
            fnold(m)=fn(m)
            sumfm(m)=sumfm(m)+0.5_r8*fmbar*dw
            fmold(m)=fm(m)
         enddo
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact &
            + sixth*(wb*(g+gold) + (g*w+gold*wold))*dw
         gold=g
         wold=w
         dw=dwnew
         if (n > 1 .and. (w > wmax .or. fnmin > fmax)) exit
         w=w+dw
         if (n == nx) then
            errmsg = subname//' -- do loop is too short in activate'
            errflg = 1
            return
         end if

      enddo

      if(w.lt.wmaxf)then

         !            contribution from all updrafts stronger than wmax
         !            assuming constant f (close to fmax)
         wnuc=w+wdiab

         z1=(w-wbar)/(sigw*sq2)
         z2=(wmaxf-wbar)/(sigw*sq2)
         g=exp(-z1*z1)
         integ=sigw*0.5_r8*sq2*sqpi*(erf(z2)-erf(z1))
         !            consider only upward flow into cloud base when estimating flux
         wf1=max(w,zero)
         zf1=(wf1-wbar)/(sigw*sq2)
         gf1=exp(-zf1*zf1)
         wf2=max(wmaxf,zero)
         zf2=(wf2-wbar)/(sigw*sq2)
         gf2=exp(-zf2*zf2)
         gf=(gf1-gf2)
         integf=wbar*sigw*0.5_r8*sq2*sqpi*(erf(zf2)-erf(zf1))+sigw*sigw*gf

         do m=1,nbins
            sumflxn(m)=sumflxn(m)+integf*fn(m)
            sumfn(m)=sumfn(m)+fn(m)*integ
            sumflxm(m)=sumflxm(m)+integf*fm(m)
            sumfm(m)=sumfm(m)+fm(m)*integ
         enddo
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact + integf
         !            sumg=sumg+integ
      endif


      do m=1,nbins
         fn(m)=sumfn(m)/(sq2*sqpi*sigw)
         !            fn(m)=sumfn(m)/(sumg)
         if(fn(m).gt.1.01_r8)then
            errmsg = 'activate -- fn > 1'
            errflg = 1
            return
         endif
         fluxn(m)=sumflxn(m)/(sq2*sqpi*sigw)
         fm(m)=sumfm(m)/(sq2*sqpi*sigw)
         !            fm(m)=sumfm(m)/(sumg)
         fluxm(m)=sumflxm(m)/(sq2*sqpi*sigw)
      enddo
      !        same form as fluxm
      flux_fullact = sumflx_fullact/(sq2*sqpi*sigw)

   else

      !        single updraft
      wnuc=wbar+wdiab

      if(wnuc.gt.0._r8)then

         w=wbar

         if(in_cloud) then

            if (smax_f > 0._r8) then
               smax = alpha*w/(2.0_r8*pi*rhoh2o*grow*gamma*smax_f)
            else
               smax = 1.e-20_r8
            end if

         else ! at cloud base
            alw        = alpha*wnuc
            sqrtalw    = sqrt(alw)
            etafactor1 = alw*sqrtalw

            do m = 1, nbins
               eta(m)  = etafactor1*etafactor2(m)
               zeta(m) = twothird*sqrtalw*aten/sqrtg
            end do
            if ( present(smax_prescribed) ) then
               smax = smax_prescribed
            else
               smax = aero_props%maxsat(zeta,eta,smc)
            end if
         end if

         do m=1,nbins

            call aero_props%actfracs( m, smc(m), smax, fn(m), fm(m) )

            if(wbar.gt.0._r8)then
               fluxn(m)=fn(m)*w
               fluxm(m)=fm(m)*w
            endif
         enddo
         flux_fullact = w
      endif

   endif

end subroutine activate_aerosol

end module aero_activate
