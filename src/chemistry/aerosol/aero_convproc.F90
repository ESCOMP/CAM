module aero_convproc
!---------------------------------------------------------------------------------
! Purpose:
!
! Portable aerosol/trace-gas convective cloud processing scheme
!
! currently these routines assume stratiform and convective clouds only interact
! through the detrainment of convective cloudborne material into stratiform clouds
!
! thus the stratiform-cloudborne aerosols (in the qqcw array) are not processed
! by the convective up/downdrafts, but are affected by the detrainment
!
! Author: R. C. Easter
!
!---------------------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8

use aerosol_properties_mod, only: aerosol_properties

implicit none
private

public :: aero_convproc_run

logical, parameter, public :: use_cwaer_for_activate_maxsat = .false.

!  method1_activate_nlayers = number of layers (including cloud base) where activation is applied
integer, parameter, public :: method1_activate_nlayers = 2
!  method2_activate_smaxmax = the uniform or peak supersat value (as 0-1 fraction = percent*0.01)
real(r8), parameter, public :: method2_activate_smaxmax = 0.003_r8

!  method_reduce_actfrac = 1 -- multiply activation fractions by factor_reduce_actfrac
!                               (this works ok with convproc_method_activate = 1 but not for ... = 2)
!                        = 2 -- do 2 iterations to get an overall reduction by factor_reduce_actfrac
!                               (this works ok with convproc_method_activate = 1 or 2)
!                        = other -- do nothing involving reduce_actfrac
integer, parameter, public :: method_reduce_actfrac = 0
real(r8), parameter, public :: factor_reduce_actfrac = 0.5_r8

!  convproc_method_activate - 1=apply abdulrazzak-ghan to entrained aerosols for lowest nlayers
!                             2=do secondary activation with prescribed supersat
integer, parameter, public :: convproc_method_activate = 2

contains

subroutine aero_convproc_run( aero_props, convtype, lchnk, dt,  &
                     t,          pmid,       q,  du,     eu,         &
                     ed,         dp,         dpdry,      jt,         &
                     mx,         ideep,      il1g,       il2g,       &
                     cldfrac,    icwmr,      rprd,       evapc,      &
                     fracice,    dqdt,       nsrflx,     qsrflx,     &
                     xx_mfup_max, xx_wcldbase, xx_kcldbase,          &
                     dcondt_resusp3d, conu2,  dcondt2,              &
                     ncol,       pver,       ncnstaer,   nbins,      &
                     pi,         rhoh2o,     rh2o,       gravit,     &
                     latvap,     cpair,      rair,                   &
                     convproc_do_evaprain_atonce,                    &
                     convproc_pom_spechygro,                         &
                     errmsg,     errflg )

!-----------------------------------------------------------------------
!
! Purpose:
! Convective transport of trace species.
! The trace species need not be conservative, and source/sink terms for
!    activation, resuspension, aqueous chemistry and gas uptake, and
!    wet removal are all applied.
! Currently this works with the ZM deep convection, but we should be able
!    to adapt it for both Hack and McCaa shallow convection
!
! Compare to subr convproc which does conservative trace species.
!
! Method:
! Computes tracer mixing ratios in updraft and downdraft "cells" in a
! Lagrangian manner, with source/sinks applied in the updraft other.
! Then computes grid-cell-mean tendencies by considering
!    updraft and downdraft fluxes across layer boundaries
!    environment subsidence/lifting fluxes across layer boundaries
!    sources and sinks in the updraft
!    resuspension of activated species in the grid-cell as a whole
!
! Note1:  A better estimate or calculation of either the updraft velocity
!         or fractional area is needed.
! Note2:  If updraft area is a small fraction of over cloud area,
!         then aqueous chemistry is underestimated.  These are both
!         research areas.
!
! Authors: O. Seland and R. Easter, based on convtran by P. Rasch
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Input arguments
!
   class(aerosol_properties), intent(in) :: aero_props

   integer,  intent(in) :: ncol             ! number of atmospheric columns
   integer,  intent(in) :: pver             ! number of vertical layers
   integer,  intent(in) :: ncnstaer         ! number of aerosol constituents (extended, = 2*ncnst)
   integer,  intent(in) :: nbins            ! number of aerosol bins/modes
   real(r8), intent(in) :: pi               ! ratio of circle circumference to diameter
   real(r8), intent(in) :: rhoh2o           ! density of liquid water (STP) (kg/m3)
   real(r8), intent(in) :: rh2o             ! gas constant for water vapor (J/K/kg)
   real(r8), intent(in) :: gravit           ! gravitational acceleration (m/s2)
   real(r8), intent(in) :: latvap           ! latent heat of vaporization (J/kg)
   real(r8), intent(in) :: cpair            ! specific heat of dry air (J/K/kg)
   real(r8), intent(in) :: rair             ! dry air gas constant (J/K/kg)
   logical,  intent(in) :: convproc_do_evaprain_atonce ! resuspend only when rain fully evaporates in a layer
   real(r8), intent(in) :: convproc_pom_spechygro      ! prescribed p-organic hygroscopicity (<0 = use default)

   character(len=*), intent(in) :: convtype  ! identifies the type of
                                             ! convection ("deep", "shcu")
   integer,  intent(in) :: lchnk             ! chunk identifier
   real(r8), intent(in) :: dt                ! Model timestep
   real(r8), intent(in) :: t(:,:)            ! Temperature
   real(r8), intent(in) :: pmid(:,:)         ! Pressure at model levels
   real(r8), intent(in) :: q(:,:,:)          ! Tracer array including moisture

   real(r8), intent(in) :: du(:,:)           ! Mass detrain rate from updraft
   real(r8), intent(in) :: eu(:,:)           ! Mass entrain rate into updraft
   real(r8), intent(in) :: ed(:,:)           ! Mass entrain rate into downdraft
! *** note1 - mu, md, eu, ed, du, dp, dpdry are GATHERED ARRAYS ***
! *** note2 - mu and md units are (mb/s), which is used in the zm_conv code
!           - eventually these should be changed to (kg/m2/s)
! *** note3 - eu, ed, du are "d(massflux)/dp" (with dp units = mb), and are all >= 0

   real(r8), intent(in) :: dp(:,:)           ! Delta pressure between interfaces (mb)
   real(r8), intent(in) :: dpdry(:,:)        ! Delta dry-pressure (mb)
   integer,  intent(in) :: jt(:)             ! Index of cloud top for each column
   integer,  intent(in) :: mx(:)             ! Index of cloud bottom for each column
   integer,  intent(in) :: ideep(:)          ! Gathering array indices
   integer,  intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer,  intent(in) :: il2g              ! Gathered max lon indices over which to operate
! *** note4 -- for il1g <= i <= il2g,  icol = ideep(i) is the "normal" chunk column index

   real(r8), intent(in) :: cldfrac(:,:)      ! Convective cloud fractional area
   real(r8), intent(in) :: icwmr(:,:)        ! Convective cloud water from zhang
   real(r8), intent(in) :: rprd(:,:)         ! Convective precipitation formation rate
   real(r8), intent(in) :: evapc(:,:)        ! Convective precipitation evaporation rate
   real(r8), intent(in) :: fracice(:,:)      ! Ice fraction of cloud droplets

   real(r8), intent(out):: dqdt(:,:,:)       ! Tracer tendency array
   integer,  intent(in) :: nsrflx            ! last dimension of qsrflx
   real(r8), intent(out):: qsrflx(:,:,:)
                              ! process-specific column tracer tendencies
                              ! (1=activation,  2=resuspension, 3=aqueous rxn,
                              !  4=wet removal, 5=renaming)
   real(r8), intent(out) :: xx_mfup_max(:)
   real(r8), intent(out) :: xx_wcldbase(:)
   real(r8), intent(out) :: xx_kcldbase(:)
   real(r8), intent(inout) :: dcondt_resusp3d(:,:,:)
   real(r8), intent(out) :: conu2(:,:,:,:)   ! updraft interface TMR diagnostic (WETC/CONU history)
   real(r8), intent(out) :: dcondt2(:,:,:,:) ! wet-deposition TMR tendency diagnostic (WETC/CONU history)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!--------------------------Local Variables------------------------------

! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = pcnst*2

   integer :: i, icol         ! Work index
   integer :: iconvtype       ! 1=deep, 2=uw shallow
   integer :: iflux_method    ! 1=as in convtran (deep), 2=simpler
   integer :: ipass_calc_updraft
   integer :: jtsub           ! Work index
   integer :: k               ! Work index
   integer :: kactcnt         ! Counter for no. of levels having activation
   integer :: kactcntb        ! Counter for activation diagnostic output
   integer :: kactfirst       ! Lowest layer with activation (= cloudbase)
   integer :: kbot            ! Cloud-flux bottom layer for current i (=mx(i))
   integer :: kbot_prevap     ! Lowest layer for doing resuspension from evaporating precip
   integer :: ktop            ! Cloud-flux top    layer for current i (=jt(i))
                              ! Layers between kbot,ktop have mass fluxes
                              !    but not all have cloud water, because the
                              !    updraft starts below the cloud base
   integer :: km1, km1x       ! Work index
   integer :: kp1, kp1x       ! Work index
   integer :: l, mm           ! Work index
   integer :: m, n            ! Work index
   integer :: nerr            ! number of errors for entire run
   integer :: nerrmax         ! maximum number of errors to report
   integer :: npass_calc_updraft
   integer :: ntsub           !

   logical  do_act_this_lev             ! flag for doing activation at current level

   real(r8) aqfrac(2,ncnstaer)       ! aqueous fraction of constituent in updraft
   real(r8) cldfrac_i(pver)          ! cldfrac at current i (with adjustments)

   real(r8) chat(2,ncnstaer,pver+1)   ! mix ratio in env at interfaces
   real(r8) cond(2,ncnstaer,pver+1)   ! mix ratio in downdraft at interfaces
   real(r8) const(2,ncnstaer,pver)   ! gathered tracer array
   real(r8) conu(2,ncnstaer,pver+1)   ! mix ratio in updraft at interfaces

   real(r8) dcondt(2,ncnstaer,pver)  ! grid-average TMR tendency for current column
   real(r8) dcondt_prevap(2,ncnstaer,pver) ! portion of dcondt from precip evaporation
   real(r8) dcondt_resusp(2,ncnstaer,pver) ! portion of dcondt from resuspension

   real(r8) dcondt_wetdep(2,ncnstaer,pver) ! portion of dcondt from wet deposition
   real(r8) dconudt_activa(2,ncnstaer,pver+1) ! d(conu)/dt by activation
   real(r8) dconudt_aqchem(2,ncnstaer,pver+1) ! d(conu)/dt by aqueous chem
   real(r8) dconudt_wetdep(2,ncnstaer,pver+1) ! d(conu)/dt by wet removal

   real(r8) maxflux(2,ncnstaer)      ! maximum (over layers) of fluxin and fluxout
   real(r8) maxflux2(2,ncnstaer)     ! ditto but computed using method-2 fluxes
   real(r8) maxprevap(2,ncnstaer)    ! maximum (over layers) of dcondt_prevap*dp
   real(r8) maxresusp(2,ncnstaer)    ! maximum (over layers) of dcondt_resusp*dp
   real(r8) maxsrce(2,ncnstaer)      ! maximum (over layers) of netsrce

   real(r8) sumflux(2,ncnstaer)      ! sum (over layers) of netflux
   real(r8) sumflux2(2,ncnstaer)     ! ditto but computed using method-2 fluxes
   real(r8) sumsrce(2,ncnstaer)      ! sum (over layers) of dp*netsrce
   real(r8) sumchng(2,ncnstaer)      ! sum (over layers) of dp*dcondt
   real(r8) sumchng3(2,ncnstaer)     ! ditto but after call to resusp_conv
   real(r8) sumprevap(2,ncnstaer)    ! sum (over layers) of dp*dcondt_prevap
   real(r8) sumwetdep(2,ncnstaer)    ! sum (over layers) of dp*dconudt_wetdep

   real(r8) cabv                 ! mix ratio of constituent above
   real(r8) cbel                 ! mix ratio of constituent below
   real(r8) cdifr                ! normalized diff between cabv and cbel
   real(r8) cdt(pver)            ! (in-updraft first order wet removal rate) * dt
   real(r8) clw_cut              ! threshold clw value for doing updraft
                                 ! transformation and removal
   real(r8) courantmax           ! maximum courant no.
   real(r8) dddp(pver)           ! dd(i,k)*dp(i,k) at current i
   real(r8) dp_i(pver)           ! dp(i,k) at current i
   real(r8) dt_u(pver)           ! lagrangian transport time in the updraft
   real(r8) dudp(pver)           ! du(i,k)*dp(i,k) at current i
   real(r8) dqdt_i(pver,ncnstaer)   ! dqdt(i,k,m) at current i
   real(r8) dtsub                ! dt/ntsub
   real(r8) dz                   ! working layer thickness (m)
   real(r8) eddp(pver)           ! ed(i,k)*dp(i,k) at current i
   real(r8) eudp(pver)           ! eu(i,k)*dp(i,k) at current i
   real(r8) expcdtm1             ! a work variable
   real(r8) fa_u(pver)           ! fractional area of in the updraft
   real(r8) fa_u_dp              ! current fa_u(k)*dp_i(k)
   real(r8) f_ent                ! fraction of the "before-detrainment" updraft
                                 ! massflux at k/k-1 interface resulting from
                                 ! entrainment of level k air
   real(r8) fluxin               ! a work variable
   real(r8) fluxout              ! a work variable
   real(r8) maxc                 ! a work variable
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) minc                 ! a work variable
   real(r8) md_m_eddp            ! a work variable
   real(r8) md_i(pver+1)         ! md(i,k) at current i (note pverp dimension)
   real(r8) md_x(pver+1)         ! md(i,k) at current i (note pverp dimension)
   real(r8) mu_i(pver+1)         ! mu(i,k) at current i (note pverp dimension)
   real(r8) mu_x(pver+1)         ! mu(i,k) at current i (note pverp dimension)
   ! md_i, md_x, mu_i, mu_x are all "dry" mass fluxes
   ! the mu_x/md_x are initially calculated from the incoming mu/md by applying dp/dpdry
   ! the mu_i/md_i are next calculated by applying the mbsth threshold
   real(r8) mu_p_eudp(pver)      ! = mu_i(kp1) + eudp(k)
   real(r8) netflux              ! a work variable
   real(r8) netsrce              ! a work variable
   real(r8) q_i(pver,ncnstaer)      ! q(i,k,m) at current i
   real(r8) qsrflx_i(ncnstaer,nsrflx) ! qsrflx(i,m,n) at current i
   real(r8) rhoair_i(pver)       ! air density at current i
   real(r8) small                ! a small number
   real(r8) tmpa                 ! work variables
   real(r8) tmpf                 ! work variables
   real(r8) xinv_ntsub           ! 1.0/ntsub
   real(r8) wup(pver)            ! working updraft velocity (m/s)
   real(r8) hund_ovr_g           ! = 100.0_r8/gravit
!  used with zm_conv mass fluxes and delta-p
!     for mu = [mbar/s],   mu*hund_ovr_g = [kg/m2/s]
!     for dp = [mbar] and q = [kg/kg],   q*dp*hund_ovr_g = [kg/m2]

   !Fractional area of ensemble mean updrafts in ZM scheme set to 0.01
   !Chosen to reproduce vertical velocities in GATEIII GIGALES (Khairoutdinov etal 2009, JAMES)
   real(r8), parameter :: zm_areafrac = 0.01_r8

!-----------------------------------------------------------------------
!
   errmsg = ''
   errflg = 0

   hund_ovr_g = 100.0_r8/gravit

   iconvtype = -1
   iflux_method = -1

   if (convtype == 'deep') then
      iconvtype = 1
      iflux_method = 1
   else if (convtype == 'uwsh') then
      iconvtype = 2
      iflux_method = 2
   else
      errmsg = '*** aero_convproc_run -- convtype is not |deep| or |uwsh|'
      errflg = 1
      return
   end if

   nerr = 0
   nerrmax = 99

   dcondt_resusp3d(:,:,:) = 0._r8

   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

   qsrflx(:,:,:) = 0.0_r8
   dqdt(:,:,:) = 0.0_r8
   xx_mfup_max(:) = 0.0_r8
   xx_wcldbase(:) = 0.0_r8
   xx_kcldbase(:) = 0.0_r8

   wup(:) = 0.0_r8

   dcondt2 = 0.0_r8
   conu2 = 0.0_r8
   aqfrac = 0.0_r8

! inititialize aqfrac to 1.0 for activated aerosol species, 0.0 otherwise
   do m = 1, aero_props%nbins()
      do l = 0, aero_props%nmasses(m)
         mm = aero_props%indexer(m,l)
         aqfrac(2,mm) = 1.0_r8
      enddo
   enddo

! Loop ever each column that has convection
! *** i is index to gathered arrays; ideep(i) is index to "normal" chunk arrays
i_loop_main_aa: &
   do i = il1g, il2g
   icol = ideep(i)


   if ( (jt(i) <= 0) .and. (mx(i) <= 0) .and. (iconvtype /= 1) ) then
! shallow conv case with jt,mx <= 0, which means there is no shallow conv
! in this column -- skip this column
      cycle i_loop_main_aa

   else if ( (jt(i) < 1) .or. (mx(i) > pver) .or. (jt(i) > mx(i)) ) then
! invalid cloudtop and cloudbase indices -- skip this column
      write(*,9010) 'illegal jt, mx', convtype, lchnk, icol, i,    &
                                      jt(i), mx(i)
9010  format( '*** aero_convproc_run error -- ', a, 5x, 'convtype = ', a /   &
              '*** lchnk, icol, il, jt, mx = ', 5(1x,i10) )
      cycle i_loop_main_aa

   else if (jt(i) == mx(i)) then
! cloudtop = cloudbase (1 layer cloud) -- skip this column
      write(*,9010) 'jt == mx', convtype, lchnk, icol, i, jt(i), mx(i)
      cycle i_loop_main_aa

   end if


!
! cloudtop and cloudbase indices are valid so proceed with calculations
!

! Load dp_i and cldfrac_i, and calc rhoair_i
      do k = 1, pver
         dp_i(k) = dpdry(i,k)
         cldfrac_i(k) = cldfrac(icol,k)
         rhoair_i(k) = pmid(icol,k)/(rair*t(icol,k))
      end do

! Calc dry mass fluxes
!    This is approximate because the updraft air is has different temp and qv than
!    the grid mean, but the whole convective parameterization is highly approximate
      mu_x(:) = 0.0_r8
      md_x(:) = 0.0_r8
! (eu-du) = d(mu)/dp -- integrate upwards, multiplying by dpdry
      do k = pver, 1, -1
         mu_x(k) = mu_x(k+1) + (eu(i,k)-du(i,k))*dp_i(k)
         xx_mfup_max(icol) = max( xx_mfup_max(icol), mu_x(k) )
      end do
! (ed) = d(md)/dp -- integrate downwards, multiplying by dpdry
      do k = 2, pver
         md_x(k) = md_x(k-1) - ed(i,k-1)*dp_i(k-1)
      end do

! Load mass fluxes over cloud layers
! (Note - use of arrays dimensioned k=1,pver+1 simplifies later coding)
! Zero out values below threshold
! Zero out values at "top of cloudtop", "base of cloudbase"
      ktop = jt(i)
      kbot = mx(i)
! usually the updraft ( & downdraft) start ( & end ) at kbot=pver, but sometimes kbot < pver
! transport, activation, resuspension, and wet removal only occur between kbot >= k >= ktop
! resuspension from evaporating precip can occur at k > kbot when kbot < pver
      kbot_prevap = pver
      mu_i(:) = 0.0_r8
      md_i(:) = 0.0_r8
      do k = ktop+1, kbot
         mu_i(k) = mu_x(k)
         if (mu_i(k) <= mbsth) mu_i(k) = 0.0_r8
         md_i(k) = md_x(k)
         if (md_i(k) >= -mbsth) md_i(k) = 0.0_r8
      end do
      mu_i(ktop) = 0.0_r8
      md_i(ktop) = 0.0_r8
      mu_i(kbot+1) = 0.0_r8
      md_i(kbot+1) = 0.0_r8

!  Compute updraft and downdraft "entrainment*dp" from eu and ed
!  Compute "detrainment*dp" from mass conservation
      eudp(:) = 0.0_r8
      dudp(:) = 0.0_r8
      eddp(:) = 0.0_r8
      dddp(:) = 0.0_r8
      courantmax = 0.0_r8
      do k = ktop, kbot
         if ((mu_i(k) > 0) .or. (mu_i(k+1) > 0)) then
            if (du(i,k) <= 0.0_r8) then
               eudp(k) = mu_i(k) - mu_i(k+1)
            else
               eudp(k) = max( eu(i,k)*dp_i(k), 0.0_r8 )
               dudp(k) = (mu_i(k+1) + eudp(k)) - mu_i(k)
               if (dudp(k) < 1.0e-12_r8*eudp(k)) then
                  eudp(k) = mu_i(k) - mu_i(k+1)
                  dudp(k) = 0.0_r8
               end if
            end if
         end if
         if ((md_i(k) < 0) .or. (md_i(k+1) < 0)) then
            eddp(k) = max( ed(i,k)*dp_i(k), 0.0_r8 )
            dddp(k) = (md_i(k+1) + eddp(k)) - md_i(k)
            if (dddp(k) < 1.0e-12_r8*eddp(k)) then
               eddp(k) = md_i(k) - md_i(k+1)
               dddp(k) = 0.0_r8
            end if
         end if
         courantmax = max( courantmax, ( mu_i(k+1)+eudp(k)-md_i(k)+eddp(k) )*dt/dp_i(k) )
      end do ! k

! number of time substeps needed to maintain "courant number" <= 1
      ntsub = 1
      if (courantmax > (1.0_r8 + 1.0e-6_r8)) then
         ntsub = 1 + int( courantmax )
      end if
      xinv_ntsub = 1.0_r8/ntsub
      dtsub = dt*xinv_ntsub
      courantmax = courantmax*xinv_ntsub

!  load tracer mixing ratio array, which will be updated at the end of each jtsub interation
      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            q_i(1:pver,mm) = q(icol,1:pver,mm)
            conu2(icol,1:pver,1,mm) = q(icol,1:pver,mm)
         end do
      end do

!
!   when method_reduce_actfrac = 2, need to do the updraft calc twice
!   (1st to get non-adjusted activation amount, 2nd to apply reduction factor)
      npass_calc_updraft = 1
      if ( (method_reduce_actfrac == 2)      .and. &
           (factor_reduce_actfrac >= 0.0_r8) .and. &
           (factor_reduce_actfrac <= 1.0_r8) ) npass_calc_updraft = 2


jtsub_loop_main_aa: &
      do jtsub = 1, ntsub


ipass_calc_updraft_loop: &
      do ipass_calc_updraft = 1, npass_calc_updraft

      qsrflx_i(:,:) = 0.0_r8
      dqdt_i(:,:) = 0.0_r8

      const = 0.0_r8 ! zero cloud-phase species
      chat = 0.0_r8 ! zero cloud-phase species
      conu = 0.0_r8
      cond = 0.0_r8

      dcondt = 0.0_r8
      dcondt_resusp = 0.0_r8
      dcondt_wetdep = 0.0_r8
      dcondt_prevap = 0.0_r8
      dconudt_aqchem = 0.0_r8
      dconudt_wetdep = 0.0_r8

! only initialize the activation tendency on ipass=1
      if (ipass_calc_updraft == 1) dconudt_activa = 0.0_r8

      ! initialize mixing ratio arrays (chat, const, conu, cond)
      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)

            const(1,mm,:) = q_i(:,mm)

            ! From now on work only with gathered data
            ! Interpolate environment tracer values to interfaces
            do k = 1,pver
               km1 = max(1,k-1)
               minc = min(const(1,mm,km1),const(1,mm,k))
               maxc = max(const(1,mm,km1),const(1,mm,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(1,mm,k)-const(1,mm,km1))/max(maxc,small)
               endif

               ! If the two layers differ significantly use a geometric averaging procedure
               ! But only do that for deep convection.  For shallow, use the simple
               ! averaging which is used in subr cmfmca
               if (iconvtype /= 1) then
                  chat(1,mm,k) = 0.5_r8* (const(1,mm,k)+const(1,mm,km1))
               else if (cdifr > 1.E-6_r8) then
                  cabv = max(const(1,mm,km1),maxc*1.e-12_r8)
                  cbel = max(const(1,mm,k),maxc*1.e-12_r8)
                  chat(1,mm,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel
               else             ! Small diff, so just arithmetic mean
                  chat(1,mm,k) = 0.5_r8* (const(1,mm,k)+const(1,mm,km1))
               end if

               ! Set provisional up and down draft values, and tendencies
               conu(1,mm,k) = chat(1,mm,k)
               cond(1,mm,k) = chat(1,mm,k)
            end do ! k

            ! Values at surface inferface == values in lowest layer
            chat(1,mm,pver+1) = const(1,mm,pver)
            conu(1,mm,pver+1) = const(1,mm,pver)
            cond(1,mm,pver+1) = const(1,mm,pver)
         end do ! l
      end do ! m



! Compute updraft mixing ratios from cloudbase to cloudtop
! No special treatment is needed at k=pver because arrays
!    are dimensioned 1:pver+1
! A time-split approach is used.  First, entrainment is applied to produce
!    an initial conu(m,k) from conu(m,k+1).  Next, chemistry/physics are
!    applied to the initial conu(m,k) to produce a final conu(m,k).
!    Detrainment from the updraft uses this final conu(m,k).
! Note that different time-split approaches would give somewhat different
!    results
      kactcnt = 0 ; kactcntb = 0 ; kactfirst = 1
k_loop_main_bb: &
      do k = kbot, ktop, -1
         kp1 = k+1

! cldfrac = conv cloud fractional area.  This could represent anvil cirrus area,
!    and may not useful for aqueous chem and wet removal calculations
         cldfrac_i(k) = max( cldfrac_i(k), 0.005_r8 )
! mu_p_eudp(k) = updraft massflux at k, without detrainment between kp1,k
         mu_p_eudp(k) = mu_i(kp1) + eudp(k)

         fa_u(k) = 0.0_r8 !BSINGH(10/15/2014): Initialized so that it has a value if the following "if" check yeilds .false.
         if (mu_p_eudp(k) > mbsth) then
! if (mu_p_eudp(k) <= mbsth) the updraft mass flux is negligible at base and top
!    of current layer,
! so current layer is a "gap" between two unconnected updrafts,
! so essentially skip all the updraft calculations for this layer

! First apply changes from entrainment
            f_ent = eudp(k)/mu_p_eudp(k)
            f_ent = max( 0.0_r8, min( 1.0_r8, f_ent ) )
            tmpa = 1.0_r8 - f_ent
            do n = 1,2 ! phase
               do m = 1, aero_props%nbins()
                  do l = 0, aero_props%nmasses(m)
                     mm = aero_props%indexer(m,l)
                     conu(n,mm,k) = tmpa*conu(n,mm,kp1) + f_ent*const(n,mm,k)
                  end do
               end do
            end do

! estimate updraft velocity (wup)
            if (iconvtype /= 1) then
! shallow - wup = (mup in kg/m2/s) / [rhoair * (updraft area)]
               wup(k) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                      / (rhoair_i(k) * (cldfrac_i(k)*0.5_r8))
            else
! deep - as in shallow, but assumed constant updraft_area with height zm_areafrac
               wup(k) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                      / (rhoair_i(k) * zm_areafrac)
            end if

! compute lagrangian transport time (dt_u) and updraft fractional area (fa_u)
! *** these must obey    dt_u(k)*mu_p_eudp(k) = dp_i(k)*fa_u(k)
            dz = dp_i(k)*hund_ovr_g/rhoair_i(k)
            dt_u(k) = dz/wup(k)
            dt_u(k) = min( dt_u(k), dt )
            fa_u(k) = dt_u(k)*(mu_p_eudp(k)/dp_i(k))


! Now apply transformation and removal changes
!    Skip levels where icwmr(icol,k) <= clw_cut (= 1.0e-6) to eliminate
!    occasional very small icwmr values from the ZM module
            clw_cut = 1.0e-6_r8


            if (convproc_method_activate <= 1) then
! aerosol activation - method 1
!    skip levels that are completely glaciated (fracice(icol,k) == 1.0)
!    when kactcnt=1 (first/lowest layer with cloud water) apply
!       activatation to the entire updraft
!    when kactcnt>1 apply activatation to the amount entrained at this level
               if ((icwmr(icol,k) > clw_cut) .and. (fracice(icol,k) < 1.0_r8)) then
                  kactcnt = kactcnt + 1

                  if ((kactcnt == 1) .or. (f_ent > 0.0_r8)) then
                     kactcntb = kactcntb + 1
                  end if

                  if (kactcnt == 1) then
                     ! diagnostic fields
                     ! xx_wcldbase = w at first cloudy layer, estimated from mu and cldfrac
                     xx_wcldbase(icol) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                         / (rhoair_i(k) * (cldfrac_i(k)*0.5_r8))
                     xx_kcldbase(icol) = k

                     kactfirst = k
                     tmpa = 1.0_r8
                     call activate_convproc( aero_props, &
                        conu(:,:,k), dconudt_activa(:,:,k), conu(:,:,k),  &
                        tmpa,       dt_u(k),            wup(k),           &
                        t(icol,k),  rhoair_i(k), ipass_calc_updraft,     &
                        ncnstaer,   nbins,                               &
                        pi, rhoh2o, rh2o, gravit, latvap, cpair, rair,   &
                        errmsg,     errflg )
                     if (errflg /= 0) return
                  else if (f_ent > 0.0_r8) then
                     ! current layer is above cloud base (=first layer with activation)
                     !    only allow activation at k = kactfirst thru kactfirst-(method1_activate_nlayers-1)
                     if (k >= kactfirst-(method1_activate_nlayers-1)) then
                        call activate_convproc( aero_props, &
                           conu(:,:,k),  dconudt_activa(:,:,k), const(:,:,k), &
                           f_ent,      dt_u(k),             wup(k),           &
                           t(icol,k),  rhoair_i(k), ipass_calc_updraft,      &
                           ncnstaer,   nbins,                                &
                           pi, rhoh2o, rh2o, gravit, latvap, cpair, rair,    &
                           errmsg,     errflg )
                        if (errflg /= 0) return
                     end if
                  end if
! the following was for cam2 shallow convection (hack),
! but is not appropriate for cam5 (uwshcu)
!                 else if ((kactcnt > 0) .and. (iconvtype /= 1)) then
! !    for shallow conv, when you move from activation occuring to
! !       not occuring, reset kactcnt=0, because the hack scheme can
! !       produce multiple "1.5 layer clouds" separated by clear air
!                    kactcnt = 0
!                 end if
               end if ! ((icwmr(icol,k) > clw_cut) .and. (fracice(icol,k) < 1.0)) then

            else ! (convproc_method_activate >= 2)
! aerosol activation - method 2
!    skip levels that are completely glaciated (fracice(icol,k) == 1.0)
!    when kactcnt=1 (first/lowest layer with cloud water)
!       apply "primary" activatation to the entire updraft
!    when kactcnt>1
!       apply secondary activatation to the entire updraft
!       do this for all levels above cloud base (even if completely glaciated)
!          (this is something for sensitivity testing)
               do_act_this_lev = .false.
               if (kactcnt <= 0) then
                  if (icwmr(icol,k) > clw_cut) then
                     do_act_this_lev = .true.
                     kactcnt = 1
                     kactfirst = k
                     ! diagnostic fields
                     ! xx_wcldbase = w at first cloudy layer, estimated from mu and cldfrac
                     xx_wcldbase(icol) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                         / (rhoair_i(k) * (cldfrac_i(k)*0.5_r8))
                     xx_kcldbase(icol) = k
                  end if
               else
!                 if ((icwmr(icol,k) > clw_cut) .and. (fracice(icol,k) < 1.0)) then
                     do_act_this_lev = .true.
                     kactcnt = kactcnt + 1
!                 end if
               end if

               if ( do_act_this_lev ) then
                  kactcntb = kactcntb + 1

                  call activate_convproc_method2( aero_props, &
                     conu(:,:,k),  dconudt_activa(:,:,k),     &
                     f_ent,      dt_u(k),   wup(k),           &
                     t(icol,k),  rhoair_i(k), k,              &
                     kactfirst,  ipass_calc_updraft,          &
                     ncnstaer,   nbins,      convproc_pom_spechygro, &
                     pi, rhoh2o, rh2o, gravit, latvap, cpair, rair,  &
                     errmsg,     errflg )
                  if (errflg /= 0) return

               end if
               conu2(icol,k,:,:) = conu(:,:,k)

            end if ! (convproc_method_activate <= 1)

! aqueous chemistry
!    do glaciated levels as aqchem_conv will eventually do acid vapor uptake
!    to ice, and aqchem_conv module checks fracice before doing liquid wtr stuff
!            if (icwmr(icol,k) > clw_cut) then
!              call aqchem_conv( conu(1,k), dconudt_aqchem(1,k), aqfrac,  &
!                 t(icol,k), fracice(icol,k), icwmr(icol,k), rhoair_i(k), &
!                 lh2o2(icol,k), lo3(icol,k), dt_u(k)                     )
!            end if

! wet removal
!
! mirage2
!    rprd               = precip formation as a grid-cell average (kgW/kgA/s)
!    icwmr              = cloud water MR within updraft area (kgW/kgA)
!    fupdr              = updraft fractional area (--)
!    A = rprd/fupdr     = precip formation rate within updraft area (kgW/kgA/s)
!    B = A/icwmr = rprd/(icwmr*fupdr)
!                       = first-order removal rate (1/s)
!    C = dp/(mup/fupdr) = updraft air residence time in the layer (s)
!
!    fraction removed = (1.0 - exp(-cdt)) where
!                 cdt = B*C = (dp/mup)*rprd/icwmr
!
!    Note1:  fupdr cancels out in cdt, so need not be specified
!    Note2:  dp & mup units need only be consistent (e.g., mb & mb/s)
!    Note3:  for shallow conv, cdt = 1-beta (beta defined in Hack scheme)
!    Note4:  the "dp" in C above and code below should be the moist dp
!
! cam5
!    clw_preloss = cloud water MR before loss to precip
!                = icwmr + dt*(rprd/fupdr)
!    B = A/clw_preloss  = (rprd/fupdr)/(icwmr + dt*rprd/fupdr)
!                       = rprd/(fupdr*icwmr + dt*rprd)
!                       = first-order removal rate (1/s)
!
!    fraction removed = (1.0 - exp(-cdt)) where
!                 cdt = B*C = (fupdr*dp/mup)*[rprd/(fupdr*icwmr + dt*rprd)]
!
!    Note1:  *** cdt is now sensitive to fupdr, which we do not really know,
!                and is not the same as the convective cloud fraction
!    Note2:  dt is appropriate in the above cdt expression, not dtsub
!
!    Apply wet removal at levels where
!       icwmr(icol,k) > clw_cut  AND  rprd(icol,k) > 0.0
!    as wet removal occurs in both liquid and ice clouds
!
            cdt(k) = 0.0_r8
            if ((icwmr(icol,k) > clw_cut) .and. (rprd(icol,k) > 0.0_r8)) then
!              if (iconvtype == 1) then
                  tmpf = 0.5_r8*cldfrac_i(k)
                  cdt(k) = (tmpf*dp(i,k)/mu_p_eudp(k)) * rprd(icol,k) / &
                        (tmpf*icwmr(icol,k) + dt*rprd(icol,k))
!              else if (k < pver) then
!                 if (eudp(k+1) > 0) cdt(k) =   &
!                       rprd(icol,k)*dp(i,k)/(icwmr(icol,k)*eudp(k+1))
!              end if
            end if
            if (cdt(k) > 0.0_r8) then
               expcdtm1 = exp(-cdt(k)) - 1.0_r8

               do m = 1, aero_props%nbins()
                  do l = 0, aero_props%nmasses(m)
                     mm = aero_props%indexer(m,l)
                     do n = 1,2
                        dconudt_wetdep(n,mm,k) = conu(n,mm,k)*aqfrac(n,mm)*expcdtm1
                        conu(n,mm,k) = conu(n,mm,k) + dconudt_wetdep(n,mm,k)
                        dconudt_wetdep(n,mm,k) = dconudt_wetdep(n,mm,k) / dt_u(k)
                        conu2(icol,k,n,mm) = conu(n,mm,k)
                     enddo
                  enddo
               enddo

            end if

         end if    ! "(mu_p_eudp(k) > mbsth)"
      end do k_loop_main_bb ! "k = kbot, ktop, -1"

! when doing updraft calcs twice, only need to go this far on the first pass
      if ( (ipass_calc_updraft == 1) .and. &
           (npass_calc_updraft == 2) ) cycle ipass_calc_updraft_loop


! Compute downdraft mixing ratios from cloudtop to cloudbase
! No special treatment is needed at k=2
! No transformation or removal is applied in the downdraft
      do k = ktop, kbot
         kp1 = k + 1
! md_m_eddp = downdraft massflux at kp1, without detrainment between k,kp1
         md_m_eddp = md_i(k) - eddp(k)
         if (md_m_eddp < -mbsth) then

            do m = 1, aero_props%nbins()
               do l = 0, aero_props%nmasses(m)
                  mm = aero_props%indexer(m,l)
                  do n = 1,2
                     cond(n,mm,kp1) = ( md_i(k)*cond(n,mm,k) &
                                      - eddp(k)*const(n,mm,k) ) / md_m_eddp
                  end do
               end do
            end do
         end if
      end do ! k


! Now computes fluxes and tendencies
! NOTE:  The approach used in convtran applies to inert tracers and
!        must be modified to include source and sink terms
      sumflux = 0.0_r8
      sumflux2 = 0.0_r8
      sumsrce = 0.0_r8
      sumchng = 0.0_r8
      sumchng3 = 0.0_r8
      sumwetdep = 0.0_r8
      sumprevap = 0.0_r8

      maxflux = 0.0_r8
      maxflux2 = 0.0_r8
      maxresusp = 0.0_r8
      maxsrce = 0.0_r8
      maxprevap = 0.0_r8

k_loop_main_cc: &
      do k = ktop, kbot
         kp1 = k+1
         km1 = k-1
         kp1x = min( kp1, pver )
         km1x = max( km1, 1 )
         fa_u_dp = fa_u(k)*dp_i(k)

         do m = 1, aero_props%nbins()
            do l = 0, aero_props%nmasses(m)
               mm = aero_props%indexer(m,l)
               do n = 1,2

                  ! First compute fluxes using environment subsidence/lifting and
                  ! entrainment/detrainment into up/downdrafts,
                  ! to provide an additional mass balance check
                  ! (this could be deleted after the code is well tested)
                  fluxin  = mu_i(k)*min(chat(n,mm,k),const(n,mm,km1x))       &
                          - md_i(kp1)*min(chat(n,mm,kp1),const(n,mm,kp1x))   &
                          + dudp(k)*conu(n,mm,k) + dddp(k)*cond(n,mm,kp1)
                  fluxout = mu_i(kp1)*min(chat(n,mm,kp1),const(n,mm,k))      &
                          - md_i(k)*min(chat(n,mm,k),const(n,mm,k))          &
                          + (eudp(k) + eddp(k))*const(n,mm,k)

                  netflux = fluxin - fluxout

                  sumflux2(n,mm) = sumflux2(n,mm) + netflux
                  maxflux2(n,mm) = max( maxflux2(n,mm), abs(fluxin), abs(fluxout) )

                  ! Now compute fluxes as in convtran, and also source/sink terms
                  ! (version 3 limit fluxes outside convection to mass in appropriate layer
                  ! (these limiters are probably only safe for positive definite quantitities
                  ! (it assumes that mu and md already satify a courant number limit of 1)
                  if (iflux_method /= 2) then
                     fluxin  =     mu_i(kp1)*conu(n,mm,kp1)                     &
                                 + mu_i(k  )*min(chat(n,mm,k  ),const(n,mm,km1x))  &
                               - ( md_i(k  )*cond(n,mm,k)                       &
                                 + md_i(kp1)*min(chat(n,mm,kp1),const(n,mm,kp1x)) )
                     fluxout =     mu_i(k  )*conu(n,mm,k)                       &
                                 + mu_i(kp1)*min(chat(n,mm,kp1),const(n,mm,k   ))  &
                               - ( md_i(kp1)*cond(n,mm,kp1)                     &
                                 + md_i(k  )*min(chat(n,mm,k  ),const(n,mm,k   )) )
                  else
                     fluxin  =     mu_i(kp1)*conu(n,mm,kp1)                     &
                               - ( md_i(k  )*cond(n,mm,k) )
                     fluxout =     mu_i(k  )*conu(n,mm,k)                       &
                               - ( md_i(kp1)*cond(n,mm,kp1) )

                     ! new method -- simple upstream method for the env subsidence
                     ! tmpa = net env mass flux (positive up) at top of layer k
                     tmpa = -( mu_i(k  ) + md_i(k  ) )
                     if (tmpa <= 0.0_r8) then
                        fluxin  = fluxin  - tmpa*const(n,mm,km1x)
                     else
                        fluxout = fluxout + tmpa*const(n,mm,k   )
                     end if
                     ! tmpa = net env mass flux (positive up) at base of layer k
                     tmpa = -( mu_i(kp1) + md_i(kp1) )
                     if (tmpa >= 0.0_r8) then
                        fluxin  = fluxin  + tmpa*const(n,mm,kp1x)
                     else
                        fluxout = fluxout - tmpa*const(n,mm,k   )
                     end if
                  end if

                  netflux = fluxin - fluxout
                  netsrce = fa_u_dp*(dconudt_aqchem(n,mm,k) + &
                       dconudt_activa(n,mm,k) + dconudt_wetdep(n,mm,k))
                  dcondt(n,mm,k) = (netflux+netsrce)/dp_i(k)

                  dcondt_wetdep(n,mm,k) = fa_u_dp*dconudt_wetdep(n,mm,k)/dp_i(k)
                  sumwetdep(n,mm) = sumwetdep(n,mm) + fa_u_dp*dconudt_wetdep(n,mm,k)

                  dcondt2(icol,k,n,mm) = dcondt(n,mm,k)

               end do
            end do

         end do
      end do k_loop_main_cc ! "k = ktop, kbot"

! calculate effects of precipitation evaporation
      call precpevap_convproc( aero_props, dcondt, dcondt_wetdep,  dcondt_prevap,   &
                                  rprd,   evapc,          dp_i,            &
                                  icol,   ktop,   pver,   ncnstaer,        &
                                  convproc_do_evaprain_atonce )

! make adjustments to dcondt for activated & unactivated aerosol species
!    pairs to account any (or total) resuspension of convective-cloudborne aerosol
      call resuspend_convproc( aero_props, dcondt, dcondt_resusp, ktop, kbot_prevap, &
                                  pver,   ncnstaer,   convproc_do_evaprain_atonce )

      ! Do resuspension of aerosols from rain only when the rain has
      ! totally evaporated.
      if (convproc_do_evaprain_atonce) then

         do m = 1, aero_props%nbins()
            do l = 0, aero_props%nmasses(m)
               mm = aero_props%indexer(m,l)
               dcondt_resusp3d(mm,icol,:) = dcondt_resusp(2,mm,:)
            end do
         end do

         dcondt_resusp(2,:,:) = 0._r8
      end if

! calculate new column-tendency variables
      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            do n = 1,2
               do k = ktop, kbot_prevap
                  sumprevap(n,mm) = sumprevap(n,mm) + dcondt_prevap(n,mm,k)*dp_i(k)
               end do
            end do
         end do
      end do

!
! note again the aero_convproc_run does not apply convective cloud processing
!    to the stratiform-cloudborne aerosol
! within this routine, cloudborne aerosols are convective-cloudborne
!
! before tendencies (dcondt, which is loaded into dqdt) are returned,
!    the convective-cloudborne aerosol tendencies must be combined
!    with the interstitial tendencies
! resuspend_convproc has already done this for the dcondt
!
! the individual process column tendencies (sumwetdep, sumprevap, ...)
!    are just diagnostic fields that can be written to history
! tendencies for interstitial and convective-cloudborne aerosol could
!    both be passed back and output, if desired
! currently, however, the interstitial and convective-cloudborne tendencies
!    are combined (in the next code block) before being passed back (in qsrflx)
!

      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            sumwetdep(1,mm) = sumwetdep(1,mm) + sumwetdep(2,mm)
            sumprevap(1,mm) = sumprevap(1,mm) + sumprevap(2,mm)
         enddo
      enddo

!
! scatter overall tendency back to full array
!
      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            do k = ktop, kbot_prevap
               dqdt_i(k,mm) = dcondt(1,mm,k)
               dqdt(icol,k,mm) = dqdt(icol,k,mm) + dqdt_i(k,mm)*xinv_ntsub
            end do

         end do
      end do ! m

! scatter column burden tendencies for various processes to qsrflx
      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            qsrflx_i(mm,4) = sumwetdep(1,mm)*hund_ovr_g
            qsrflx_i(mm,5) = sumprevap(1,mm)*hund_ovr_g
            qsrflx(icol,mm,1:5) = qsrflx(icol,mm,1:5) + qsrflx_i(mm,1:5)*xinv_ntsub
         end do
      end do

      if (jtsub < ntsub) then
         ! update the q_i for the next interation of the jtsub loop
         do m = 1, aero_props%nbins()
            do l = 0, aero_props%nmasses(m)
               mm = aero_props%indexer(m,l)
               do k = ktop, kbot_prevap
                  q_i(k,mm) = max( (q_i(k,mm) + dqdt_i(k,mm)*dtsub), 0.0_r8 )
               end do
            end do
         end do
      end if

      end do ipass_calc_updraft_loop

      end do jtsub_loop_main_aa  ! of the main "do jtsub = 1, ntsub" loop


   end do i_loop_main_aa  ! of the main "do i = il1g, il2g" loop

! conu2/dcondt2 are returned as out-args; the WETC/CONU history diagnostics
!    are written by the CAM host layer (aero_convproc_cam)

end subroutine aero_convproc_run

!=========================================================================================
   subroutine precpevap_convproc(  aero_props,      &
              dcondt,  dcondt_wetdep, dcondt_prevap,           &
              rprd,    evapc,         dp_i,                    &
              icol,    ktop,          pver,      ncnstaer,     &
              convproc_do_evaprain_atonce )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate resuspension of wet-removed aerosol species resulting
!    from precip evaporation
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! arguments
! (note:  TMR = tracer mixing ratio)

   class(aerosol_properties), intent(in) :: aero_props
   integer,  intent(in)    :: pver              ! number of vertical layers
   integer,  intent(in)    :: ncnstaer          ! number of aerosol constituents (extended)
   logical,  intent(in)    :: convproc_do_evaprain_atonce
   real(r8), intent(inout) :: dcondt(2,ncnstaer,pver)
                              ! overall TMR tendency from convection
   real(r8), intent(in)    :: dcondt_wetdep(2,ncnstaer,pver)
                              ! portion of TMR tendency due to wet removal
   real(r8), intent(inout) :: dcondt_prevap(2,ncnstaer,pver)
                              ! portion of TMR tendency due to precip evaporation
                              ! (actually, due to the adjustments made here)
                              ! (on entry, this is 0.0)

   real(r8), intent(in)    :: rprd(:,:)  ! conv precip production  rate (gathered)
   real(r8), intent(in)    :: evapc(:,:)  ! conv precip evaporation rate (gathered)
   real(r8), intent(in)    :: dp_i(pver) ! pressure thickness of level (in mb)

   integer,  intent(in)    :: icol  ! normal (ungathered) i index for current column
   integer,  intent(in)    :: ktop  ! index of top cloud level for current column

!-----------------------------------------------------------------------
! local variables
   integer  :: k, l, m, mm, n
   real(r8) :: del_pr_flux_prod      ! change to precip flux from production  [(kg/kg/s)*mb]
   real(r8) :: del_pr_flux_evap      ! change to precip flux from evaporation [(kg/kg/s)*mb]
   real(r8) :: del_wd_flux_evap      ! change to wet deposition flux from evaporation [(kg/kg/s)*mb]
   real(r8) :: fdel_pr_flux_evap     ! fractional change to precip flux from evaporation
   real(r8) :: pr_flux               ! precip flux at base of current layer [(kg/kg/s)*mb]
   real(r8) :: pr_flux_old
   real(r8) :: tmpdp                 ! delta-pressure (mb)
   real(r8) :: wd_flux(2,ncnstaer)   ! tracer wet deposition flux at base of current layer [(kg/kg/s)*mb]
!-----------------------------------------------------------------------

   pr_flux = 0.0_r8
   wd_flux = 0.0_r8

   do k = ktop, pver
      tmpdp = dp_i(k)

      pr_flux_old = pr_flux
      del_pr_flux_prod = tmpdp*max(0.0_r8, rprd(icol,k))
      pr_flux = pr_flux_old + del_pr_flux_prod

      del_pr_flux_evap = min( pr_flux, tmpdp*max(0.0_r8, evapc(icol,k)) )

      ! Do resuspension of aerosols from rain only when the rain has
      ! totally evaporated in one layer.
      if (convproc_do_evaprain_atonce .and. &
          (del_pr_flux_evap.ne.pr_flux)) del_pr_flux_evap = 0._r8

      fdel_pr_flux_evap = del_pr_flux_evap / max(pr_flux, 1.0e-35_r8)

      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            do n = 1,2

               ! use -dcondt_wetdep(m,k) as it is negative (or zero)
               wd_flux(n,mm) = wd_flux(n,mm) + tmpdp*max(0.0_r8, -dcondt_wetdep(n,mm,k))
               del_wd_flux_evap = wd_flux(n,mm)*fdel_pr_flux_evap

               dcondt_prevap(n,mm,k) = del_wd_flux_evap/tmpdp

            end do
         end do
      end do

      ! resuspension --> create larger aerosols
      if (convproc_do_evaprain_atonce) then
         call aero_props%resuspension_resize( dcondt_prevap(1,:,k) )
      endif

      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)
            do n = 1,2
               dcondt(n,mm,k) = dcondt(n,mm,k) + dcondt_prevap(n,mm,k)
            end do
         end do
      end do

      pr_flux = max( 0.0_r8, pr_flux-del_pr_flux_evap )

   end do ! k

   end subroutine precpevap_convproc

!=========================================================================================
   subroutine activate_convproc( aero_props,    &
              conu,       dconudt,   conent,    &
              f_ent,      dt_u,      wup,       &
              tair,       rhoair, ipass_calc_updraft, &
              ncnstaer,   nbins,                 &
              pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
              errmsg,     errflg )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate activation of aerosol species in convective updraft
! for a single column and level
!
! Method:
! conu(l)    = Updraft TMR (tracer mixing ratio) at k/k-1 interface
! conent(l)  = TMR of air that is entrained into the updraft from level k
! f_ent      = Fraction of the "before-detrainment" updraft massflux at
!              k/k-1 interface" resulting from entrainment of level k air
!              (where k is the current level in subr aero_convproc_run)
!
! On entry to this routine, the conu(l) represents the updraft TMR
! after entrainment, but before chemistry/physics and detrainment,
! and is equal to
!    conu(l) = f_ent*conent(l) + (1.0-f_ent)*conu_below(l)
! where
!    conu_below(l) = updraft TMR at the k+1/k interface, and
!    f_ent   = (eudp/mu_p_eudp) is the fraction of the updraft massflux
!              from level k entrainment
!
! This routine applies aerosol activation to the entrained tracer,
! then adjusts the conu so that on exit,
!   conu(la) = conu_incoming(la) - f_ent*conent(la)*f_act(la)
!   conu(lc) = conu_incoming(lc) + f_ent*conent(la)*f_act(la)
! where
!   la, lc   = indices for an unactivated/activated aerosol component pair
!   f_act    = fraction of conent(la) that is activated.  The f_act are
!              calculated with the Razzak-Ghan activation parameterization.
!              The f_act differ for each mode, and for number/surface/mass.
!
! Note:  At the lowest layer with cloud water, subr convproc calls this
! routine with conent==conu and f_ent==1.0, with the result that
! activation is applied to the entire updraft tracer flux
!
! *** The updraft velocity used for activation calculations is rather
!     uncertain and needs more work.  However, an updraft of 1-3 m/s
!     will activate essentially all of accumulation and coarse mode particles.
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

   use aero_activate, only: activate_aerosol

!-----------------------------------------------------------------------
! arguments  (note:  TMR = tracer mixing ratio)

   class(aerosol_properties), intent(in) :: aero_props

   integer,  intent(in)    :: ncnstaer ! number of aerosol constituents (extended)
   integer,  intent(in)    :: nbins    ! number of aerosol bins/modes
   real(r8), intent(in)    :: pi       ! ratio of circle circumference to diameter
   real(r8), intent(in)    :: rhoh2o   ! density of liquid water (STP) (kg/m3)
   real(r8), intent(in)    :: rh2o     ! gas constant for water vapor (J/K/kg)
   real(r8), intent(in)    :: gravit   ! gravitational acceleration (m/s2)
   real(r8), intent(in)    :: latvap   ! latent heat of vaporization (J/kg)
   real(r8), intent(in)    :: cpair    ! specific heat of dry air (J/K/kg)
   real(r8), intent(in)    :: rair     ! dry air gas constant (J/K/kg)

   ! conu = tracer mixing ratios in updraft at top of this (current) level
   !        The conu are changed by activation
   real(r8), intent(inout) :: conu(2,ncnstaer)
   ! conent = TMRs in the entrained air at this level
   real(r8), intent(in)    :: conent(2,ncnstaer)
   real(r8), intent(inout) :: dconudt(2,ncnstaer) ! TMR tendencies due to activation

   real(r8), intent(in)    :: f_ent  ! fraction of updraft massflux that was
                                     ! entrained across this layer == eudp/mu_p_eudp
   real(r8), intent(in)    :: dt_u   ! lagrangian transport time (s) in the
                                     ! updraft at current level
   real(r8), intent(in)    :: wup    ! mean updraft vertical velocity (m/s)
                                     ! at current level updraft

   real(r8), intent(in)    :: tair   ! Temperature in Kelvin
   real(r8), intent(in)    :: rhoair ! air density (kg/m3)

   integer,  intent(in)    :: ipass_calc_updraft

!-----------------------------------------------------------------------
! local variables
   integer  :: l, m, mm

   real(r8) :: delact           ! working variable
   real(r8) :: dt_u_inv         ! 1.0/dt_u
   real(r8) :: fluxm(nbins)     ! to understand this, see subr activate_aerosol
   real(r8) :: fluxn(nbins)     ! to understand this, see subr activate_aerosol
   real(r8) :: flux_fullact     ! to understand this, see subr activate_aerosol
   real(r8) :: fm(nbins)        ! mass fraction of aerosols activated
   real(r8) :: fn(nbins)        ! number fraction of aerosols activated
   real(r8) :: hygro(nbins)     ! current hygroscopicity for int+act
   real(r8) :: naerosol(nbins)  ! interstitial+activated number conc (#/m3)
   real(r8) :: sigw             ! standard deviation of updraft velocity (cm/s)
   real(r8) :: tmp_fact         ! working variable
   real(r8) :: vaerosol(nbins)  ! int+act volume (m3/m3)
   real(r8) :: wbar             ! mean updraft velocity (cm/s)
   real(r8) :: wdiab            ! diabatic vertical velocity (cm/s)
   real(r8) :: wminf, wmaxf     ! limits for integration over updraft spectrum (cm/s)

   real(r8) :: spec_hygro
   real(r8) :: spec_dens
   character(len=32) :: spec_type

   real(r8) :: tmpa, tmpb, tmpc ! working variable
   real(r8) :: naerosol_a(1,1)    ! number conc (1/m3)
   real(r8) :: vaerosol_a(1,1)    ! volume conc (m3/m3)

   character(len=*), intent(out) :: errmsg   ! error message from activate_aerosol
   integer,          intent(out) :: errflg   ! error flag from activate_aerosol

!-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

! when ipass_calc_updraft == 2, apply the activation tendencies
!    from pass 1, but multiplied by factor_reduce_actfrac
! (can only have ipass_calc_updraft == 2 when method_reduce_actfrac = 2)
   if (ipass_calc_updraft == 2) then

      dt_u_inv = 1.0_r8/dt_u

      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)

            delact = dconudt(2,mm)*dt_u * factor_reduce_actfrac
            delact = min( delact, conu(1,mm) )
            delact = max( delact, 0.0_r8 )
            conu(1,mm) = conu(1,mm) - delact
            conu(2,mm) = conu(2,mm) + delact
            dconudt(1,mm) = -delact*dt_u_inv
            dconudt(2,mm) =  delact*dt_u_inv

         end do
      end do

      return

   end if ! (ipass_calc_updraft == 2)

! check f_ent > 0
   if (f_ent <= 0.0_r8) return

   hygro = 0.0_r8
   vaerosol = 0.0_r8
   naerosol = 0.0_r8

   do m = 1, nbins
! compute a (or a+cw) volume and hygroscopicity
      tmpa = 0.0_r8
      tmpb = 0.0_r8
      do l = 1, aero_props%nmasses(m)

         mm = aero_props%indexer(m,l)

         call aero_props%get(m, l, spectype=spec_type, density=spec_dens, hygro=spec_hygro)

         tmpc = max( conent(1,mm), 0.0_r8 )
         if ( use_cwaer_for_activate_maxsat ) &
         tmpc = tmpc + max( conent(2,mm), 0.0_r8 )
         tmpc = tmpc / spec_dens
         tmpa = tmpa + tmpc
         tmpb = tmpb + tmpc * spec_hygro
      end do
      vaerosol(m) = tmpa * rhoair
      if (tmpa < 1.0e-35_r8) then
         hygro(m) = 0.2_r8
      else
         hygro(m) = tmpb/tmpa
      end if

! load a (or a+cw) number and bound it
      tmpa = max( conent(1,mm), 0.0_r8 )
      if ( use_cwaer_for_activate_maxsat ) &
      tmpa = tmpa + max( conent(2,mm), 0.0_r8 )
      naerosol(m) = tmpa * rhoair

      naerosol_a(1,1) = naerosol(m)
      vaerosol_a(1,1) = vaerosol(m)

      call aero_props%apply_number_limits( naerosol_a, vaerosol_a, 1, 1, m )

      naerosol(m) = naerosol_a(1,1)
   end do

! call Razzak-Ghan activation routine with single updraft
   wbar = max( wup, 0.5_r8 )  ! force wbar >= 0.5 m/s for now
   sigw = 0.0_r8
   wdiab = 0.0_r8
   wminf = wbar
   wmaxf = wbar

   call activate_aerosol(                                                    &
         wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,                    &
         naerosol, nbins, vaerosol, hygro, aero_props,            &
         fn, fm, fluxn, fluxm, flux_fullact,                              &
         pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, errmsg, errflg    )
   if (errflg /= 0) return

! apply the activation fractions to the updraft aerosol mixing ratios
   dt_u_inv = 1.0_r8/dt_u

   do m = 1, aero_props%nbins()
      do l = 0, aero_props%nmasses(m)
         mm = aero_props%indexer(m,l)

         if ( (method_reduce_actfrac == 1)      .and. &
              (factor_reduce_actfrac >= 0.0_r8) .and. &
              (factor_reduce_actfrac <  1.0_r8) )     &
              tmp_fact = tmp_fact * factor_reduce_actfrac

         delact = min( conent(1,mm)*tmp_fact*f_ent, conu(1,mm) )
         delact = max( delact, 0.0_r8 )
         conu(1,mm) = conu(1,mm) - delact
         conu(2,mm) = conu(2,mm) + delact
         dconudt(1,mm) = -delact*dt_u_inv
         dconudt(2,mm) =  delact*dt_u_inv
      end do
   end do

   end subroutine activate_convproc

!=========================================================================================
   subroutine activate_convproc_method2( aero_props, &
              conu,       dconudt,             &
              f_ent,      dt_u,     wup,       &
              tair,       rhoair,   k,         &
              kactfirst,  ipass_calc_updraft,  &
              ncnstaer,   nbins,    convproc_pom_spechygro, &
              pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
              errmsg,     errflg )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate activation of aerosol species in convective updraft
! for a single column and level
!
! Method:
! conu(l)    = Updraft TMR (tracer mixing ratio) at k/k-1 interface
! f_ent      = Fraction of the "before-detrainment" updraft massflux at
!              k/k-1 interface" resulting from entrainment of level k air
!              (where k is the current level in subr aero_convproc_run)
!
! On entry to this routine, the conu(l) represents the updraft TMR
! after entrainment, but before chemistry/physics and detrainment.
!
! This routine applies aerosol activation to the conu tracer mixing ratios,
! then adjusts the conu so that on exit,
!   conu(la) = conu_incoming(la) - conu(la)*f_act(la)
!   conu(lc) = conu_incoming(lc) + conu(la)*f_act(la)
! where
!   la, lc   = indices for an unactivated/activated aerosol component pair
!   f_act    = fraction of conu(la) that is activated.  The f_act are
!              calculated with the Razzak-Ghan activation parameterization.
!              The f_act differ for each mode, and for number/surface/mass.
!
! At cloud base (k==kactfirst), primary activation is done using the
! "standard" code in subr activate do diagnose maximum supersaturation.
! Above cloud base, secondary activation is done using a
! prescribed supersaturation.
!
! *** The updraft velocity used for activation calculations is rather
!     uncertain and needs more work.  However, an updraft of 1-3 m/s
!     will activate essentially all of accumulation and coarse mode particles.
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

   use aero_activate, only: activate_aerosol

!-----------------------------------------------------------------------
! arguments  (note:  TMR = tracer mixing ratio)

   class(aerosol_properties), intent(in) :: aero_props

   integer,  intent(in)    :: ncnstaer ! number of aerosol constituents (extended)
   integer,  intent(in)    :: nbins    ! number of aerosol bins/modes
   real(r8), intent(in)    :: convproc_pom_spechygro ! prescribed p-organic hygroscopicity (<0 = use default)
   real(r8), intent(in)    :: pi       ! ratio of circle circumference to diameter
   real(r8), intent(in)    :: rhoh2o   ! density of liquid water (STP) (kg/m3)
   real(r8), intent(in)    :: rh2o     ! gas constant for water vapor (J/K/kg)
   real(r8), intent(in)    :: gravit   ! gravitational acceleration (m/s2)
   real(r8), intent(in)    :: latvap   ! latent heat of vaporization (J/kg)
   real(r8), intent(in)    :: cpair    ! specific heat of dry air (J/K/kg)
   real(r8), intent(in)    :: rair     ! dry air gas constant (J/K/kg)

   ! conu = tracer mixing ratios in updraft at top of this (current) level
   !        The conu are changed by activation
   real(r8), intent(inout) :: conu(2,ncnstaer)
   real(r8), intent(inout) :: dconudt(2,ncnstaer) ! TMR tendencies due to activation

   real(r8), intent(in)    :: f_ent  ! fraction of updraft massflux that was
                                     ! entrained across this layer == eudp/mu_p_eudp
   real(r8), intent(in)    :: dt_u   ! lagrangian transport time (s) in the
                                     ! updraft at current level
   real(r8), intent(in)    :: wup    ! mean updraft vertical velocity (m/s)
                                     ! at current level updraft

   real(r8), intent(in)    :: tair   ! Temperature in Kelvin
   real(r8), intent(in)    :: rhoair ! air density (kg/m3)
                                     ! used as in-cloud wet removal rate
   integer,  intent(in)    :: k      ! level index
   integer,  intent(in)    :: kactfirst ! k at cloud base
   integer,  intent(in)    :: ipass_calc_updraft

!-----------------------------------------------------------------------
! local variables
   integer  :: l, m, mm

   real(r8) :: delact           ! working variable
   real(r8) :: dt_u_inv         ! 1.0/dt_u
   real(r8) :: fluxm(nbins)     ! to understand this, see subr activate_aerosol
   real(r8) :: fluxn(nbins)     ! to understand this, see subr activate_aerosol
   real(r8) :: flux_fullact     ! to understand this, see subr activate_aerosol
   real(r8) :: fm(nbins)        ! mass fraction of aerosols activated
   real(r8) :: fn(nbins)        ! number fraction of aerosols activated
   real(r8) :: hygro(nbins)     ! current hygroscopicity for int+act
   real(r8) :: naerosol(nbins)  ! interstitial+activated number conc (#/m3)
   real(r8) :: sigw             ! standard deviation of updraft velocity (cm/s)
   real(r8) :: smax_prescribed  ! prescribed supersaturation for secondary activation (0-1 fraction)
   real(r8) :: tmp_fact         ! working variable
   real(r8) :: vaerosol(nbins)  ! int+act volume (m3/m3)
   real(r8) :: wbar             ! mean updraft velocity (cm/s)
   real(r8) :: wdiab            ! diabatic vertical velocity (cm/s)
   real(r8) :: wminf, wmaxf     ! limits for integration over updraft spectrum (cm/s)

   real(r8) :: spec_hygro
   real(r8) :: spec_dens
   character(len=32) :: spec_type

   real(r8) :: tmpa, tmpb, tmpc ! working variable
   real(r8) :: naerosol_a(1,1)    ! number conc (1/m3)
   real(r8) :: vaerosol_a(1,1)    ! volume conc (m3/m3)

   character(len=*), intent(out) :: errmsg   ! error message from activate_aerosol
   integer,          intent(out) :: errflg   ! error flag from activate_aerosol

!-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

! when ipass_calc_updraft == 2, apply the activation tendencies
!    from pass 1, but multiplied by factor_reduce_actfrac
! (can only have ipass_calc_updraft == 2 when method_reduce_actfrac = 2)

   if (ipass_calc_updraft == 2) then

      dt_u_inv = 1.0_r8/dt_u

      do m = 1, aero_props%nbins()
         do l = 0, aero_props%nmasses(m)
            mm = aero_props%indexer(m,l)

            delact = dconudt(2,mm)*dt_u * factor_reduce_actfrac
            delact = min( delact, conu(1,mm) )
            delact = max( delact, 0.0_r8 )
            conu(1,mm) = conu(1,mm) - delact
            conu(2,mm) = conu(2,mm) + delact
            dconudt(1,mm) = -delact*dt_u_inv
            dconudt(2,mm) =  delact*dt_u_inv
         end do
      end do   ! "n = 1, ntot_amode"
      return

   end if ! (ipass_calc_updraft == 2)

! check f_ent > 0
   if (f_ent <= 0.0_r8) return

   hygro = 0.0_r8
   vaerosol = 0.0_r8
   naerosol = 0.0_r8

   do m = 1, nbins
! compute a (or a+cw) volume and hygroscopicity
      tmpa = 0.0_r8
      tmpb = 0.0_r8
      do l = 1, aero_props%nspecies(m)

         mm = aero_props%indexer(m,l)

         call aero_props%get(m, l, spectype=spec_type, density=spec_dens, hygro=spec_hygro)

         tmpc = max( conu(1,mm), 0.0_r8 )
         if ( use_cwaer_for_activate_maxsat ) &
         tmpc = tmpc + max( conu(2,mm), 0.0_r8 )
         tmpc = tmpc / spec_dens
         tmpa = tmpa + tmpc

         ! Change the hygroscopicity of POM based on the discussion with Prof.
         ! Xiaohong Liu. Some observational studies found that the primary organic
         ! material from biomass burning emission shows very high hygroscopicity.
         ! Also, found that BC mass will be overestimated if all the aerosols in
         ! the primary mode are free to be removed. Therefore, set the hygroscopicity
         ! of POM here as 0.2 to enhance the wet scavenge of primary BC and POM.

         if (spec_type=='p-organic' .and. convproc_pom_spechygro>0._r8) then
            tmpb = tmpb + tmpc * convproc_pom_spechygro
         else
            tmpb = tmpb + tmpc * spec_hygro
         end if
      end do
      vaerosol(m) = tmpa * rhoair
      if (tmpa < 1.0e-35_r8) then
         hygro(m) = 0.2_r8
      else
         hygro(m) = tmpb/tmpa
      end if

      mm = aero_props%indexer(m,0)

! load a (or a+cw) number and bound it
      tmpa = max( conu(1,mm), 0.0_r8 )
      if ( use_cwaer_for_activate_maxsat ) &
      tmpa = tmpa + max( conu(2,mm), 0.0_r8 )
      naerosol(m) = tmpa * rhoair

      naerosol_a(1,1) = naerosol(m)
      vaerosol_a(1,1) = vaerosol(m)

      call aero_props%apply_number_limits( naerosol_a, vaerosol_a, 1, 1, m )

      naerosol(m) = naerosol_a(1,1)

   end do

! call Razzak-Ghan activation routine with single updraft
   wbar = max( wup, 0.5_r8 )  ! force wbar >= 0.5 m/s for now
   sigw = 0.0_r8
   wdiab = 0.0_r8
   wminf = wbar
   wmaxf = wbar

   if (k == kactfirst) then

      call activate_aerosol(                                                 &
         wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,                    &
         naerosol, nbins, vaerosol, hygro, aero_props,            &
         fn, fm, fluxn, fluxm, flux_fullact,                              &
         pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, errmsg, errflg    )
      if (errflg /= 0) return

   else
! above cloud base - do secondary activation with prescribed supersat
! that is constant with height
      smax_prescribed = method2_activate_smaxmax
      call activate_aerosol(                                                 &
         wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,                    &
         naerosol, nbins, vaerosol, hygro, aero_props,            &
         fn, fm, fluxn, fluxm, flux_fullact,                              &
         pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, errmsg, errflg,   &
         smax_prescribed                                                  )
      if (errflg /= 0) return
   end if

! apply the activation fractions to the updraft aerosol mixing ratios
   dt_u_inv = 1.0_r8/dt_u

   do m = 1, aero_props%nbins()
      do l = 0, aero_props%nmasses(m)
         mm = aero_props%indexer(m,l)
         if (l==0) then
            tmp_fact = fn(m)
         else
            tmp_fact = fm(m)
         end if

         if ( (method_reduce_actfrac == 1)      .and. &
              (factor_reduce_actfrac >= 0.0_r8) .and. &
              (factor_reduce_actfrac <  1.0_r8) )     &
              tmp_fact = tmp_fact * factor_reduce_actfrac

         delact = min( conu(1,mm)*tmp_fact, conu(1,mm) )
         delact = max( delact, 0.0_r8 )
         conu(1,mm) = conu(1,mm) - delact
         conu(2,mm) = conu(2,mm) + delact
         dconudt(1,mm) = -delact*dt_u_inv
         dconudt(2,mm) =  delact*dt_u_inv
      end do
   end do

   end subroutine activate_convproc_method2

!=========================================================================================
   subroutine resuspend_convproc( aero_props, &
              dcondt,  dcondt_resusp, ktop,  kbot_prevap, &
              pver,    ncnstaer,      convproc_do_evaprain_atonce )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate resuspension of activated aerosol species resulting from both
!    detrainment from updraft and downdraft into environment
!    subsidence and lifting of environment, which may move air from
!       levels with large-scale cloud to levels with no large-scale cloud
!
! Method:
! Three possible approaches were considered:
!
! 1. Ad-hoc #1 approach.  At each level, adjust dcondt for the activated
!    and unactivated portions of a particular aerosol species so that the
!    ratio of dcondt (activated/unactivate) is equal to the ratio of the
!    mixing ratios before convection.
!    THIS WAS IMPLEMENTED IN MIRAGE2
!
! 2. Ad-hoc #2 approach.  At each level, adjust dcondt for the activated
!    and unactivated portions of a particular aerosol species so that the
!    change to the activated portion is minimized (zero if possible).  The
!    would minimize effects of convection on the large-scale cloud.
!    THIS IS CURRENTLY IMPLEMENTED IN CAM5 where we assume that convective
!    clouds have no impact on the stratiform-cloudborne aerosol
!
! 3. Mechanistic approach that treats the details of interactions between
!    the large-scale and convective clouds.  (Something for the future.)
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! arguments
! (note:  TMR = tracer mixing ratio)

   class(aerosol_properties), intent(in) :: aero_props
   integer,  intent(in)    :: pver              ! number of vertical layers
   integer,  intent(in)    :: ncnstaer          ! number of aerosol constituents (extended)
   logical,  intent(in)    :: convproc_do_evaprain_atonce
   real(r8), intent(inout) :: dcondt(2,ncnstaer,pver)
                              ! overall TMR tendency from convection
   real(r8), intent(inout) :: dcondt_resusp(2,ncnstaer,pver)
                              ! portion of TMR tendency due to resuspension
                              ! (actually, due to the adjustments made here)
   integer,  intent(in)    :: ktop, kbot_prevap ! indices of top and bottom cloud levels

!-----------------------------------------------------------------------
! local variables
   integer  :: k, l, m, mm
   real(r8) :: qdota, qdotc, qdotac  ! working variables (MR tendencies)
   !-----------------------------------------------------------------------

   ! apply adjustments to dcondt for pairs of unactivated and
   ! activated aerosol species
   do m = 1, aero_props%nbins()
      do l = 0, aero_props%nmasses(m)
         mm = aero_props%indexer(m,l)

         do k = ktop, kbot_prevap
            if (convproc_do_evaprain_atonce) then
               dcondt_resusp(1,mm,k) = dcondt(1,mm,k)
               dcondt_resusp(2,mm,k) = dcondt(2,mm,k)
            else
               qdota = dcondt(1,mm,k)
               qdotc = dcondt(2,mm,k)
               qdotac = qdota + qdotc

               dcondt(1,mm,k) = qdotac
               dcondt(2,mm,k) = 0.0_r8

               dcondt_resusp(1,mm,k) = (dcondt(1,mm,k) - qdota)
               dcondt_resusp(2,mm,k) = (dcondt(2,mm,k) - qdotc)
            end if
         end do

      end do
   end do

   end subroutine resuspend_convproc

!=========================================================================================

end module aero_convproc
