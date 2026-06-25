module ndrop

!---------------------------------------------------------------------------------
! Purpose:
!   Droplet activation and vertical mixing by modal (or bin) aerosols
!   (dropmixnuc).  Portable science routines split from the CAM interface,
!   which now lives in microp_aero.F90.  Aerosol access is polymorphic through
!   aerosol_properties/aerosol_state; host physical constants are provided by
!   ndrop_init; array sizing is by runtime ncol/pver.
!
! ***N.B.*** This module is currently hardcoded to recognize only the modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8, shr_kind_cs
use shr_spfn_mod,     only: erf => shr_spfn_erf

use aerosol_properties_mod, only: aerosol_properties
use aerosol_state_mod, only: aerosol_state, ptr2d_t

use aero_activate,    only: aero_activate_init, activate_aerosol

implicit none
private
save

public ndrop_init, dropmixnuc
public psat   ! needed by the CAM interface to size the ccn diagnostic

! mathematical constants
real(r8), parameter :: zero     = 0._r8
real(r8), parameter :: third    = 1._r8/3._r8
real(r8), parameter :: twothird = 2._r8*third
real(r8), parameter :: sq2      = sqrt(2._r8)
real(r8), parameter :: surften  = 0.076_r8

! CCN diagnostic fields
integer,  parameter :: psat=6    ! number of supersaturations to calc ccn concentration
real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
                       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)

! host physical constants (set by ndrop_init)
real(r8) :: pi           ! pi
real(r8) :: rhoh2o       ! density of liquid water (kg/m3)
real(r8) :: mwh2o        ! molecular weight of water (kg/kmol)
real(r8) :: r_universal  ! universal gas constant (J/K/kmol)
real(r8) :: rh2o         ! water vapor gas constant (J/K/kg)
real(r8) :: gravit       ! gravitational acceleration (m/s2)
real(r8) :: latvap       ! latent heat of vaporization (J/kg)
real(r8) :: cpair        ! specific heat of dry air (J/K/kg)
real(r8) :: rair         ! dry air gas constant (J/K/kg)

real(r8) :: sq2pi        ! sqrt(2*pi), derived from host pi by ndrop_init

integer :: nbin  ! number of bins

!===============================================================================
contains
!===============================================================================

subroutine ndrop_init(aero_props, pi_in, rhoh2o_in, mwh2o_in, r_universal_in, &
                      rh2o_in, gravit_in, latvap_in, cpair_in, rair_in)

   class(aerosol_properties), intent(in) :: aero_props

   ! host physical constants
   real(r8), intent(in) :: pi_in           ! pi
   real(r8), intent(in) :: rhoh2o_in       ! density of liquid water (kg/m3)
   real(r8), intent(in) :: mwh2o_in        ! molecular weight of water (kg/kmol)
   real(r8), intent(in) :: r_universal_in  ! universal gas constant (J/K/kmol)
   real(r8), intent(in) :: rh2o_in         ! water vapor gas constant (J/K/kg)
   real(r8), intent(in) :: gravit_in       ! gravitational acceleration (m/s2)
   real(r8), intent(in) :: latvap_in       ! latent heat of vaporization (J/kg)
   real(r8), intent(in) :: cpair_in        ! specific heat of dry air (J/K/kg)
   real(r8), intent(in) :: rair_in         ! dry air gas constant (J/K/kg)

   !-------------------------------------------------------------------------------

   pi          = pi_in
   rhoh2o      = rhoh2o_in
   mwh2o       = mwh2o_in
   r_universal = r_universal_in
   rh2o        = rh2o_in
   gravit      = gravit_in
   latvap      = latvap_in
   cpair       = cpair_in
   rair        = rair_in

   sq2pi = sqrt(2._r8*pi)

   call aero_activate_init(mwh2o, r_universal, rhoh2o, pi)

   nbin = aero_props%nbins()

end subroutine ndrop_init

!===============================================================================

subroutine dropmixnuc( aero_props, aero_state, &
   ncol, pver, top_lev, dtmicro, &
   temp, pmid, pint, pdel, rpdel, zm, kvh, ncldwtr, &
   wsub, wmixmin, cldn, cldo, cldliqf, &
   dotend, raertend_out, tendnd, factnum, &
   wtke, nsource, ndropmix, ndropcol, &
   ccn, coltend, coltend_cw, &
   errmsg, errflg)

   ! vertical diffusion and nucleation of cloud droplets
   ! assume cloud presence controlled by cloud fraction
   ! doesn't distinguish between warm, cold clouds

   ! arguments
   class(aerosol_properties), intent(in) :: aero_props
   class(aerosol_state), intent(in) :: aero_state

   integer,  intent(in)    :: ncol        ! number of columns
   integer,  intent(in)    :: pver        ! number of vertical layers
   integer,  intent(in)    :: top_lev     ! top level for cloud physics
   real(r8), intent(in)    :: dtmicro     ! time step for microphysics (s)
   real(r8), intent(in)    :: temp(:,:)   ! temperature (K)
   real(r8), intent(in)    :: pmid(:,:)   ! mid-level pressure (Pa)
   real(r8), intent(in)    :: pint(:,:)   ! pressure at layer interfaces (Pa)
   real(r8), intent(in)    :: pdel(:,:)   ! pressure thickess of layer (Pa)
   real(r8), intent(in)    :: rpdel(:,:)  ! inverse of pressure thickess of layer (/Pa)
   real(r8), intent(in)    :: zm(:,:)     ! geopotential height of level (m)
   real(r8), intent(in)    :: kvh(:,:)    ! vertical diffusivity (m2/s), interfaces
   real(r8), intent(in)    :: ncldwtr(:,:)! droplet number concentration (#/kg)
   real(r8), intent(in)    :: wsub(:,:)   ! subgrid vertical velocity
   real(r8), intent(in)    :: wmixmin     ! minimum turbulence vertical velocity (m/s)
   real(r8), intent(in)    :: cldn(:,:)   ! cloud fraction
   real(r8), intent(in)    :: cldo(:,:)   ! cloud fraction on previous time step
   real(r8), intent(in)    :: cldliqf(:,:)! liquid cloud fraction (liquid / (liquid + ice))
   logical,  intent(in)    :: dotend(:)   ! (nele_tot) true for aerosol elements resolving to
                                          ! advected constituents: tendency returned in raertend_out.
                                          ! false elements are updated in place through the
                                          ! aero_state interstitial pointers.

   ! output arguments
   real(r8), intent(out) :: raertend_out(:,:,:) ! (ncol,pver,nele_tot) tendency of interstitial aerosol
                                                ! mass, number mixing ratios, only where dotend is true
   real(r8), intent(out) :: tendnd(:,:)   ! change in droplet number concentration (#/kg/s)
   real(r8), intent(out) :: factnum(:,:,:)     ! activation fraction for aerosol number
   real(r8), intent(out) :: wtke(:,:)     ! turbulent vertical velocity at base of layer k (m/s)
   real(r8), intent(out) :: nsource(:,:)  ! droplet number source (#/kg/s)
   real(r8), intent(out) :: ndropmix(:,:) ! droplet number mixing (#/kg/s)
   real(r8), intent(out) :: ndropcol(:)   ! column droplet number (#/m2)
   real(r8), intent(out) :: ccn(:,:,:)    ! (ncol,pver,psat) number conc of aerosols activated at supersat (#/cm3)
   real(r8), intent(out) :: coltend(:,:)  ! (ncol,nele_tot) column tendency for diagnostic output
   real(r8), intent(out) :: coltend_cw(:,:) ! (ncol,nele_tot) column tendency
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg
   !--------------------Local storage-------------------------------------

   integer  :: nele_tot            ! total number of aerosol elements

   type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
   type(ptr2d_t), allocatable :: qqcw(:)
   real(r8) :: raertend(pver)  ! tendency of aerosol mass, number mixing ratios
   real(r8) :: qqcwtend(pver)  ! tendency of cloudborne aerosol mass, number mixing ratios

   real(r8), parameter :: zkmin = 0.01_r8, zkmax = 100._r8
   integer  :: i, k, l, m, mm, n
   integer  :: km1, kp1
   integer  :: nnew, nsav, ntemp
   integer  :: nsubmix, nsubmix_bnd
   integer, save :: count_submix(100)
   integer  :: phase ! phase of aerosol

   real(r8) :: arg
   real(r8) :: dtinv

   real(r8) :: dtmin, tinv, dtt
   real(r8) :: lcldn(ncol,pver)
   real(r8) :: lcldo(ncol,pver)

   real(r8) :: zs(pver) ! inverse of distance between levels (m)
   real(r8) :: qcld(pver) ! cloud droplet number mixing ratio (#/kg)
   real(r8) :: qncld(pver)     ! droplet number nucleated on cloud boundaries
   real(r8) :: srcn(pver)       ! droplet source rate (/s)
   real(r8) :: cs(ncol,pver)      ! air density (kg/m3)
   real(r8) :: csbot(pver)       ! air density at bottom (interface) of layer (kg/m3)
   real(r8) :: csbot_cscen(pver) ! csbot(i)/cs(i,k)
   real(r8) :: dz(ncol,pver)      ! geometric thickness of layers (m)

   real(r8) :: wtke_cen(ncol,pver) ! turbulent vertical velocity at center of layer k (m/s)
   real(r8) :: wbar, wmix, wmin, wmax

   real(r8) :: zn(pver)   ! g/pdel (m2/g) for layer
   real(r8) :: flxconv    ! convergence of flux into lowest layer

   real(r8) :: wdiab           ! diabatic vertical velocity
   real(r8) :: ekd(pver)       ! diffusivity for droplets (m2/s)
   real(r8) :: ekk(0:pver)     ! density*diffusivity for droplets (kg/m3 m2/s)
   real(r8) :: ekkp(pver)      ! zn*zs*density*diffusivity
   real(r8) :: ekkm(pver)      ! zn*zs*density*diffusivity

   real(r8) :: dum, dumc
   real(r8) :: tmpa
   real(r8) :: dact
   real(r8) :: fluxntot         ! (#/cm2/s)
   real(r8) :: dtmix
   real(r8) :: alogarg
   real(r8) :: overlapp(pver), overlapm(pver) ! cloud overlap

   real(r8) :: cldo_tmp, cldn_tmp
   real(r8) :: tau_cld_regenerate
   real(r8) :: taumix_internal_pver_inv ! 1/(internal mixing time scale for k=pver) (1/s)


   real(r8), allocatable :: nact(:,:)  ! fractional aero. number  activation rate (/s)
   real(r8), allocatable :: mact(:,:)  ! fractional aero. mass    activation rate (/s)

   real(r8), allocatable :: raercol(:,:,:)    ! single column of aerosol mass, number mixing ratios
   real(r8), allocatable :: raercol_cw(:,:,:) ! same as raercol but for cloud-borne phase


   real(r8) :: na(ncol,pver,nbin), va(ncol,pver,nbin), hy(ncol,pver,nbin)
   real(r8), allocatable :: naermod(:)  ! (1/m3)
   real(r8), allocatable :: hygro(:)    ! hygroscopicity of aerosol mode
   real(r8), allocatable :: vaerosol(:) ! interstit+activated aerosol volume conc (cm3/cm3)

   real(r8) :: source(pver)

   real(r8), allocatable :: fn(:)              ! activation fraction for aerosol number
   real(r8), allocatable :: fm(:)              ! activation fraction for aerosol mass

   real(r8), allocatable :: fluxn(:)           ! number  activation fraction flux (cm/s)
   real(r8), allocatable :: fluxm(:)           ! mass    activation fraction flux (cm/s)
   real(r8)              :: flux_fullact(pver) ! 100%    activation fraction flux (cm/s)
   !     note:  activation fraction fluxes are defined as
   !     fluxn = [flux of activated aero. number into cloud (#/cm2/s)]
   !           / [aero. number conc. in updraft, just below cloudbase (#/cm3)]


   integer :: errnum
   character(len=shr_kind_cs) :: errstr
   !-------------------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   nele_tot  = aero_props%ncnst_tot()

   ! Create the liquid weighted cloud fractions that were passsed in
   ! before. This doesn't seem like the best variable, since the cloud could
   ! have liquid condensate, but the part of it that is changing could be the
   ! ice portion; however, this is what was done before.
   lcldo(:ncol,:)  = cldo(:ncol,:)  * cldliqf(:ncol,:)
   lcldn(:ncol,:) = cldn(:ncol,:) * cldliqf(:ncol,:)


   arg = 1.0_r8
   if (abs(0.8427_r8 - erf(arg))/0.8427_r8 > 0.001_r8) then
      write(errmsg,*) 'dropmixnuc: Error function error, erf(1.0) = ',ERF(arg)
      errflg = 1
      return
   endif
   arg = 0.0_r8
   if (erf(arg) /= 0.0_r8) then
      write(errmsg,*) 'dropmixnuc: Error function error, erf(0.0) = ',erf(arg)
      errflg = 1
      return
   endif

   dtinv = 1._r8/dtmicro

   allocate( &
      nact(pver,nbin),             &
      mact(pver,nbin),             &
      raer(nele_tot),              &
      qqcw(nele_tot),              &
      raercol(pver,nele_tot,2),    &
      raercol_cw(pver,nele_tot,2), &
      naermod(nbin),               &
      hygro(nbin),                 &
      vaerosol(nbin),              &
      fn(nbin),                    &
      fm(nbin),                    &
      fluxn(nbin),                 &
      fluxm(nbin)               )

   ! Init pointers to mode number and specie mass mixing ratios in
   ! intersitial and cloud borne phases.
   call aero_state%get_states( aero_props, raer, qqcw )

   factnum = 0._r8
   wtke = 0._r8
   nsource = 0._r8
   ndropmix = 0._r8
   ndropcol = 0._r8
   tendnd = 0._r8
   raertend_out = 0._r8

   ! air density (kg/m3)
   cs(:ncol,:)  = pmid(:ncol,:)/(rair*temp(:ncol,:))

   phase = 1 ! interstitial
   do m = 1, nbin
      call aero_state%loadaer( aero_props, &
           ncol, pver, &
           m, cs, phase, na(:,:,m), va(:,:,m), &
           hy(:,:,m), errnum, errstr)
      if (errnum/=0) then
         errmsg = 'dropmixnuc : '//trim(errstr)
         errflg = 1
         return
      end if
   end do

   ! overall_main_i_loop
   do i = 1, ncol

      do k = top_lev, pver-1
         zs(k) = 1._r8/(zm(i,k) - zm(i,k+1))
      end do
      zs(pver) = zs(pver-1)

      ! load number nucleated into qcld on cloud boundaries

      do k = top_lev, pver

         qcld(k)  = ncldwtr(i,k)
         qncld(k) = 0._r8
         srcn(k)  = 0._r8
         dz(i,k)  = 1._r8/(cs(i,k)*gravit*rpdel(i,k)) ! layer thickness in m

         do m = 1, nbin
            nact(k,m) = 0._r8
            mact(k,m) = 0._r8
         end do

         zn(k) = gravit*rpdel(i,k)

         if (k < pver) then
            ekd(k)   = kvh(i,k+1)
            ekd(k)   = max(ekd(k), zkmin)
            ekd(k)   = min(ekd(k), zkmax)
            csbot(k) = 2.0_r8*pint(i,k+1)/(rair*(temp(i,k) + temp(i,k+1)))
            csbot_cscen(k) = csbot(k)/cs(i,k)
         else
            ekd(k)   = 0._r8
            csbot(k) = cs(i,k)
            csbot_cscen(k) = 1.0_r8
         end if

         ! rce-comment - define wtke at layer centers for new-cloud activation
         !    and at layer boundaries for old-cloud activation
         wtke_cen(i,k) = wsub(i,k)
         wtke(i,k)     = wsub(i,k)
         wtke_cen(i,k) = max(wtke_cen(i,k), wmixmin)
         wtke(i,k)     = max(wtke(i,k), wmixmin)

         nsource(i,k) = 0._r8

      end do

      nsav = 1
      nnew = 2

      do mm = 1,nele_tot
         raercol_cw(:,mm,nsav) = 0.0_r8
         raercol(:,mm,nsav)    = 0.0_r8
         raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(i,top_lev:pver)
         raercol(top_lev:pver,mm,nsav)    = raer(mm)%fld(i,top_lev:pver)
      end do

      ! droplet nucleation/aerosol activation

      ! tau_cld_regenerate = time scale for regeneration of cloudy air
      !    by (horizontal) exchange with clear air
      tau_cld_regenerate = 3600.0_r8 * 3.0_r8

      ! k-loop for growing/shrinking cloud calcs .............................
      ! grow_shrink_main_k_loop: &
      do k = top_lev, pver

         ! This code was designed for liquid clouds, but the cloudbourne
         ! aerosol can be either from liquid or ice clouds. For the ice clouds,
         ! we do not do regeneration, but as cloud fraction decreases the
         ! aerosols should be returned interstitial. The lack of a liquid cloud
         ! should not mean that all of the aerosol is realease. Therefor a
         ! section has been added for shrinking ice clouds and checks were added
         ! to protect ice cloudbourne aerosols from being released when no
         ! liquid cloud is present.

         ! shrinking ice cloud ......................................................
         cldo_tmp = cldo(i,k) * (1._r8 - cldliqf(i,k))
         cldn_tmp = cldn(i,k) * (1._r8 - cldliqf(i,k))

         if (cldn_tmp < cldo_tmp) then

            ! convert activated aerosol to interstitial in decaying cloud

            dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * (1._r8 - cldliqf(i,k))
            do mm = 1,nele_tot
               dact   = raercol_cw(k,mm,nsav)*dumc
               raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact   ! cloud-borne aerosol
               raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
            end do

         end if

         ! shrinking liquid cloud ......................................................
         !    treat the reduction of cloud fraction from when cldn(i,k) < cldo(i,k)
         !    and also dissipate the portion of the cloud that will be regenerated
         cldo_tmp = lcldo(i,k)
         cldn_tmp = lcldn(i,k) * exp( -dtmicro/tau_cld_regenerate )
         !    alternate formulation
         !    cldn_tmp = cldn(i,k) * max( 0.0_r8, (1.0_r8-dtmicro/tau_cld_regenerate) )

         ! fraction is also provided.
         if (cldn_tmp < cldo_tmp) then
            !  droplet loss in decaying cloud
            !++ sungsup
            nsource(i,k) = nsource(i,k) + qcld(k)*(cldn_tmp - cldo_tmp)/cldo_tmp*cldliqf(i,k)*dtinv
            qcld(k)      = qcld(k)*(1._r8 + (cldn_tmp - cldo_tmp)/cldo_tmp)
            !-- sungsup

            ! convert activated aerosol to interstitial in decaying cloud

            dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * cldliqf(i,k)
            do mm = 1,nele_tot
               dact   = raercol_cw(k,mm,nsav)*dumc
               raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact   ! cloud-borne aerosol
               raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
            end do

         end if

         ! growing liquid cloud ......................................................
         !    treat the increase of cloud fraction from when cldn(i,k) > cldo(i,k)
         !    and also regenerate part of the cloud
         cldo_tmp = cldn_tmp
         cldn_tmp = lcldn(i,k)

         if (cldn_tmp-cldo_tmp > 0.01_r8) then

            ! rce-comment - use wtke at layer centers for new-cloud activation
            wbar  = wtke_cen(i,k)
            wmix  = 0._r8
            wmin  = 0._r8
            wmax  = 10._r8
            wdiab = 0._r8

            ! load aerosol properties, assuming external mixtures

            do m = 1, nbin
               naermod(m)  = na(i,k,m)
               vaerosol(m) = va(i,k,m)
               hygro(m)    = hy(i,k,m)
            end do

            call activate_aerosol( &
               wbar, wmix, wdiab, wmin, wmax,                       &
               temp(i,k), cs(i,k), naermod, nbin,             &
               vaerosol, hygro, aero_props, fn, fm, fluxn,     &
               fluxm,flux_fullact(k),                          &
               pi, rhoh2o, rh2o, gravit, latvap, cpair, rair,  &
               errmsg, errflg)
            if (errflg /= 0) return

            factnum(i,k,:) = fn

            dumc = (cldn_tmp - cldo_tmp)

            do m = 1, nbin
               mm = aero_props%indexer(m,0)
               dact   = dumc*fn(m)*raer(mm)%fld(i,k) ! interstitial only
               qcld(k) = qcld(k) + dact
               nsource(i,k) = nsource(i,k) + dact*dtinv
               raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
               raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               dum = dumc*fm(m)
               do l = 1,aero_props%nmasses(m)
                  mm = aero_props%indexer(m,l)
                  dact    = dum*raer(mm)%fld(i,k) ! interstitial only
                  raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               enddo
            enddo

         endif

      enddo  ! grow_shrink_main_k_loop
      ! end of k-loop for growing/shrinking cloud calcs ......................

      ! ......................................................................
      ! start of k-loop for calc of old cloud activation tendencies ..........
      !
      ! rce-comment
      !    changed this part of code to use current cloud fraction (cldn) exclusively
      !    consider case of cldo(:)=0, cldn(k)=1, cldn(k+1)=0
      !    previous code (which used cldo below here) would have no cloud-base activation
      !       into layer k.  however, activated particles in k mix out to k+1,
      !       so they are incorrectly depleted with no replacement

      ! old_cloud_main_k_loop
      do k = top_lev, pver
         kp1 = min0(k+1, pver)
         taumix_internal_pver_inv = 0.0_r8

         if (lcldn(i,k) > 0.01_r8) then

            wdiab = 0._r8
            wmix  = 0._r8                       ! single updraft
            wbar  = wtke(i,k)                   ! single updraft
            if (k == pver) wbar = wtke_cen(i,k) ! single updraft
            wmax  = 10._r8
            wmin  = 0._r8

            if (lcldn(i,k) - lcldn(i,kp1) > 0.01_r8 .or. k == pver) then

               ! cloud base

               ! ekd(k) = wtke(i,k)*dz(i,k)/sq2pi
               ! rce-comments
               !   first, should probably have 1/zs(k) here rather than dz(i,k) because
               !      the turbulent flux is proportional to ekd(k)*zs(k),
               !      while the dz(i,k) is used to get flux divergences
               !      and mixing ratio tendency/change
               !   second and more importantly, using a single updraft velocity here
               !      means having monodisperse turbulent updraft and downdrafts.
               !      The sq2pi factor assumes a normal draft spectrum.
               !      The fluxn/fluxm from activate must be consistent with the
               !      fluxes calculated in explmix.
               ekd(k) = wbar/zs(k)

               alogarg = max(1.e-20_r8, 1._r8/lcldn(i,k) - 1._r8)
               wmin    = wbar + wmix*0.25_r8*sq2pi*log(alogarg)

               do m = 1, nbin
                  ! rce-comment - use kp1 here as old-cloud activation involves
                  !   aerosol from layer below
                  naermod(m)  = na(i,kp1,m)
                  vaerosol(m) = va(i,kp1,m)
                  hygro(m)    = hy(i,kp1,m)
               end do

               call activate_aerosol( &
                  wbar, wmix, wdiab, wmin, wmax,                       &
                  temp(i,k), cs(i,k), naermod, nbin,             &
                  vaerosol, hygro, aero_props, fn, fm, fluxn,     &
                  fluxm, flux_fullact(k),                        &
                  pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
                  errmsg, errflg)
               if (errflg /= 0) return

               factnum(i,k,:) = fn

               if (k < pver) then
                  dumc = lcldn(i,k) - lcldn(i,kp1)
               else
                  dumc = lcldn(i,k)
               endif

               fluxntot = 0

               ! rce-comment 1
               !    flux of activated mass into layer k (in kg/m2/s)
               !       = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(k)
               !    source of activated mass (in kg/kg/s) = flux divergence
               !       = actmassflux/(cs(i,k)*dz(i,k))
               !    so need factor of csbot_cscen = csbot(k)/cs(i,k)
               !                   dum=1./(dz(i,k))
               dum=csbot_cscen(k)/(dz(i,k))

               ! rce-comment 2
               !    code for k=pver was changed to use the following conceptual model
               !    in k=pver, there can be no cloud-base activation unless one considers
               !       a scenario such as the layer being partially cloudy,
               !       with clear air at bottom and cloudy air at top
               !    assume this scenario, and that the clear/cloudy portions mix with
               !       a timescale taumix_internal = dz(i,pver)/wtke_cen(i,pver)
               !    in the absence of other sources/sinks, qact (the activated particle
               !       mixratio) attains a steady state value given by
               !          qact_ss = fcloud*fact*qtot
               !       where fcloud is cloud fraction, fact is activation fraction,
               !       qtot=qact+qint, qint is interstitial particle mixratio
               !    the activation rate (from mixing within the layer) can now be
               !       written as
               !          d(qact)/dt = (qact_ss - qact)/taumix_internal
               !                     = qtot*(fcloud*fact*wtke/dz) - qact*(wtke/dz)
               !    note that (fcloud*fact*wtke/dz) is equal to the nact/mact
               !    also, d(qact)/dt can be negative.  in the code below
               !       it is forced to be >= 0
               !
               ! steve --
               !    you will likely want to change this.  i did not really understand
               !       what was previously being done in k=pver
               !    in the cam3_5_3 code, wtke(i,pver) appears to be equal to the
               !       droplet deposition velocity which is quite small
               !    in the cam3_5_37 version, wtke is done differently and is much
               !       larger in k=pver, so the activation is stronger there
               !
               if (k == pver) then
                  taumix_internal_pver_inv = flux_fullact(k)/dz(i,k)
               end if

               do m = 1, nbin
                  mm = aero_props%indexer(m,0)
                  fluxn(m) = fluxn(m)*dumc
                  fluxm(m) = fluxm(m)*dumc
                  nact(k,m) = nact(k,m) + fluxn(m)*dum
                  mact(k,m) = mact(k,m) + fluxm(m)*dum
                  if (k < pver) then
                     ! note that kp1 is used here
                     fluxntot = fluxntot &
                        + fluxn(m)*raercol(kp1,mm,nsav)*cs(i,k)
                  else
                     tmpa = raercol(kp1,mm,nsav)*fluxn(m) &
                          + raercol_cw(kp1,mm,nsav)*(fluxn(m) &
                          - taumix_internal_pver_inv*dz(i,k))
                     fluxntot = fluxntot + max(0.0_r8, tmpa)*cs(i,k)
                  end if
               end do
               srcn(k)      = srcn(k) + fluxntot/(cs(i,k)*dz(i,k))
               nsource(i,k) = nsource(i,k) + fluxntot/(cs(i,k)*dz(i,k))

            endif  ! (cldn(i,k) - cldn(i,kp1) > 0.01 .or. k == pver)

         else

            ! no liquid cloud
            nsource(i,k) = nsource(i,k) - qcld(k)*dtinv
            qcld(k)      = 0

            if (cldn(i,k) < 0.01_r8) then
               ! no ice cloud either

               ! convert activated aerosol to interstitial in decaying cloud

               do mm = 1,nele_tot
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) + raercol_cw(k,mm,nsav)  ! cloud-borne aerosol
                  raercol_cw(k,mm,nsav) = 0._r8
               end do

            end if
         end if

      end do  ! old_cloud_main_k_loop

      ! switch nsav, nnew so that nnew is the updated aerosol
      ntemp = nsav
      nsav  = nnew
      nnew  = ntemp

      ! load new droplets in layers above, below clouds

      dtmin     = dtmicro
      ekk(top_lev-1)    = 0.0_r8
      ekk(pver) = 0.0_r8
      do k = top_lev, pver-1
         ! rce-comment -- ekd(k) is eddy-diffusivity at k/k+1 interface
         !   want ekk(k) = ekd(k) * (density at k/k+1 interface)
         !   so use pint(i,k+1) as pint is 1:pverp
         !           ekk(k)=ekd(k)*2.*pint(i,k)/(rair*(temp(i,k)+temp(i,k+1)))
         !           ekk(k)=ekd(k)*2.*pint(i,k+1)/(rair*(temp(i,k)+temp(i,k+1)))
         ekk(k) = ekd(k)*csbot(k)
      end do

      do k = top_lev, pver
         km1     = max0(k-1, top_lev)
         ekkp(k) = zn(k)*ekk(k)*zs(k)
         ekkm(k) = zn(k)*ekk(k-1)*zs(km1)
         tinv    = ekkp(k) + ekkm(k)

         ! rce-comment -- tinv is the sum of all first-order-loss-rates
         !    for the layer.  for most layers, the activation loss rate
         !    (for interstitial particles) is accounted for by the loss by
         !    turb-transfer to the layer above.
         !    k=pver is special, and the loss rate for activation within
         !    the layer must be added to tinv.  if not, the time step
         !    can be too big, and explmix can produce negative values.
         !    the negative values are reset to zero, resulting in an
         !    artificial source.
         if (k == pver) tinv = tinv + taumix_internal_pver_inv

         if (tinv .gt. 1.e-6_r8) then
            dtt   = 1._r8/tinv
            dtmin = min(dtmin, dtt)
         end if
      end do

      dtmix   = 0.9_r8*dtmin
      nsubmix = int(dtmicro/dtmix) + 1
      if (nsubmix > 100) then
         nsubmix_bnd = 100
      else
         nsubmix_bnd = nsubmix
      end if
      count_submix(nsubmix_bnd) = count_submix(nsubmix_bnd) + 1
      dtmix = dtmicro/nsubmix

      do k = top_lev, pver
         kp1 = min(k+1, pver)
         km1 = max(k-1, top_lev)
         ! maximum overlap assumption
         if (cldn(i,kp1) > 1.e-10_r8) then
            overlapp(k) = min(cldn(i,k)/cldn(i,kp1), 1._r8)
         else
            overlapp(k) = 1._r8
         end if
         if (cldn(i,km1) > 1.e-10_r8) then
            overlapm(k) = min(cldn(i,k)/cldn(i,km1), 1._r8)
         else
            overlapm(k) = 1._r8
         end if
      end do


      ! rce-comment
      !    the activation source(k) = mact(k,m)*raercol(kp1,lmass)
      !       should not exceed the rate of transfer of unactivated particles
      !       from kp1 to k which = ekkp(k)*raercol(kp1,lmass)
      !    however it might if things are not "just right" in subr activate
      !    the following is a safety measure to avoid negatives in explmix
      do k = top_lev, pver-1
         do m = 1, nbin
            nact(k,m) = min( nact(k,m), ekkp(k) )
            mact(k,m) = min( mact(k,m), ekkp(k) )
         end do
      end do


      ! old_cloud_nsubmix_loop
      do n = 1, nsubmix
         qncld(:) = qcld(:)
         ! switch nsav, nnew so that nsav is the updated aerosol
         ntemp   = nsav
         nsav    = nnew
         nnew    = ntemp
         srcn(:) = 0.0_r8

         do m = 1, nbin
            mm = aero_props%indexer(m,0)

            ! update droplet source
            ! rce-comment- activation source in layer k involves particles from k+1
            !	       srcn(:)=srcn(:)+nact(:,m)*(raercol(:,mm,nsav))
            srcn(top_lev:pver-1) = srcn(top_lev:pver-1) + nact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))

            ! rce-comment- new formulation for k=pver
            !              srcn(  pver  )=srcn(  pver  )+nact(  pver  ,m)*(raercol(  pver,mm,nsav))
            tmpa = raercol(pver,mm,nsav)*nact(pver,m) &
                 + raercol_cw(pver,mm,nsav)*(nact(pver,m) - taumix_internal_pver_inv)
            srcn(pver) = srcn(pver) + max(0.0_r8,tmpa)
         end do
         call explmix(  &
            qcld, srcn, ekkp, ekkm, overlapp,  &
            overlapm, qncld, zero, zero, pver, &
            top_lev, dtmix, .false.)

         ! rce-comment
         !    the interstitial particle mixratio is different in clear/cloudy portions
         !    of a layer, and generally higher in the clear portion.  (we have/had
         !    a method for diagnosing the the clear/cloudy mixratios.)  the activation
         !    source terms involve clear air (from below) moving into cloudy air (above).
         !    in theory, the clear-portion mixratio should be used when calculating
         !    source terms
         do m = 1, nbin
            mm = aero_props%indexer(m,0)
            ! rce-comment -   activation source in layer k involves particles from k+1
            !	              source(:)= nact(:,m)*(raercol(:,mm,nsav))
            source(top_lev:pver-1) = nact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))
            ! rce-comment - new formulation for k=pver
            !               source(  pver  )= nact(  pver,  m)*(raercol(  pver,mm,nsav))
            tmpa = raercol(pver,mm,nsav)*nact(pver,m) &
                 + raercol_cw(pver,mm,nsav)*(nact(pver,m) - taumix_internal_pver_inv)
            source(pver) = max(0.0_r8, tmpa)
            flxconv = 0._r8

            call explmix( &
                 raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                 overlapm, raercol_cw(:,mm,nsav), zero, zero, pver,   &
                 top_lev, dtmix, .false.)

            call explmix( &
                 raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                 overlapm, raercol(:,mm,nsav), zero, flxconv, pver, &
                 top_lev, dtmix, .true., raercol_cw(:,mm,nsav))

            do l = 1,aero_props%nmasses(m)
               mm = aero_props%indexer(m,l)
               ! rce-comment -   activation source in layer k involves particles from k+1
               !	          source(:)= mact(:,m)*(raercol(:,mm,nsav))
               source(top_lev:pver-1) = mact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))
               ! rce-comment- new formulation for k=pver
               !                 source(  pver  )= mact(  pver  ,m)*(raercol(  pver,mm,nsav))
               tmpa = raercol(pver,mm,nsav)*mact(pver,m) &
                    + raercol_cw(pver,mm,nsav)*(mact(pver,m) - taumix_internal_pver_inv)
               source(pver) = max(0.0_r8, tmpa)
               flxconv = 0._r8

               call explmix( &
                    raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                    overlapm, raercol_cw(:,mm,nsav), zero, zero, pver,   &
                    top_lev, dtmix, .false.)

               call explmix( &
                    raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                    overlapm, raercol(:,mm,nsav), zero, flxconv, pver, &
                    top_lev, dtmix, .true., raercol_cw(:,mm,nsav))

            end do
         end do

      end do ! old_cloud_nsubmix_loop

      ! evaporate particles again if no cloud (either ice or liquid)

      do k = top_lev, pver
         if (cldn(i,k) == 0._r8) then
            ! no ice or liquid cloud
            qcld(k)=0._r8

            ! convert activated aerosol to interstitial in decaying cloud
            do mm = 1,nele_tot
               raercol(k,mm,nnew)    = raercol(k,mm,nnew) + raercol_cw(k,mm,nnew)
               raercol_cw(k,mm,nnew) = 0._r8
            end do

         end if
      end do

      ! droplet number

      ndropcol(i) = 0._r8
      do k = top_lev, pver
         ndropmix(i,k) = (qcld(k) - ncldwtr(i,k))*dtinv - nsource(i,k)
         tendnd(i,k)   = (max(qcld(k), 1.e-6_r8) - ncldwtr(i,k))*dtinv
         ndropcol(i)   = ndropcol(i) + ncldwtr(i,k)*pdel(i,k)
      end do
      ndropcol(i) = ndropcol(i)/gravit

      raertend = 0._r8
      qqcwtend = 0._r8

      do m = 1, nbin
         do l = 0, aero_props%nmasses(m)

            mm = aero_props%indexer(m,l)

            raertend(top_lev:pver) = (raercol(top_lev:pver,mm,nnew) - raer(mm)%fld(i,top_lev:pver))*dtinv
            qqcwtend(top_lev:pver) = (raercol_cw(top_lev:pver,mm,nnew) - qqcw(mm)%fld(i,top_lev:pver))*dtinv

            coltend(i,mm)    = sum( pdel(i,:)*raertend )/gravit
            coltend_cw(i,mm) = sum( pdel(i,:)*qqcwtend )/gravit

            ! check for advected aerosol constituents
            if (dotend(mm)) then ! advected aerosol parts
               raertend_out(i,:,mm) = 0.0_r8
               raertend_out(i,top_lev:pver,mm) = raertend(top_lev:pver)      ! set tendencies for interstitial aerosol
            else
               raer(mm)%fld(i,:) = 0.0_r8
               raer(mm)%fld(i,top_lev:pver)  = raercol(top_lev:pver,mm,nnew) ! update non-advected interstitial aerosol (pbuf)
            end if

            qqcw(mm)%fld(i,:) = 0.0_r8
            qqcw(mm)%fld(i,top_lev:pver) = raercol_cw(top_lev:pver,mm,nnew)  ! update cloud-borne aerosol

         end do
      end do

   end do  ! overall_main_i_loop
   ! end of main loop over i/longitude ....................................

   call ccncalc(aero_state, aero_props, ncol, pver, top_lev, temp, cs, ccn, errmsg, errflg)
   if (errflg /= 0) return

   deallocate( &
      nact,       &
      mact,       &
      raer,       &
      qqcw,       &
      raercol,    &
      raercol_cw, &
      naermod,    &
      hygro,      &
      vaerosol,   &
      fn,         &
      fm,         &
      fluxn,      &
      fluxm       )

end subroutine dropmixnuc

!===============================================================================

subroutine explmix( q, src, ekkp, ekkm, overlapp, overlapm, &
   qold, surfrate, flxconv, pver, top_lev, dt, is_unact, qactold )

   !  explicit integration of droplet/aerosol mixing
   !     with source due to activation/nucleation


   integer, intent(in) :: pver ! number of levels
   integer, intent(in) :: top_lev ! top level for cloud physics
   real(r8), intent(out) :: q(pver) ! mixing ratio to be updated
   real(r8), intent(in) :: qold(pver) ! mixing ratio from previous time step
   real(r8), intent(in) :: src(pver) ! source due to activation/nucleation (/s)
   real(r8), intent(in) :: ekkp(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
   ! below layer k  (k,k+1 interface)
   real(r8), intent(in) :: ekkm(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
   ! above layer k  (k,k+1 interface)
   real(r8), intent(in) :: overlapp(pver) ! cloud overlap below
   real(r8), intent(in) :: overlapm(pver) ! cloud overlap above
   real(r8), intent(in) :: surfrate ! surface exchange rate (/s)
   real(r8), intent(in) :: flxconv ! convergence of flux from surface
   real(r8), intent(in) :: dt ! time step (s)
   logical, intent(in) :: is_unact ! true if this is an unactivated species
   real(r8), intent(in),optional :: qactold(pver)
   ! mixing ratio of ACTIVATED species from previous step
   ! *** this should only be present
   !     if the current species is unactivated number/sfc/mass

   integer k,kp1,km1

   if ( is_unact ) then
      !     the qactold*(1-overlap) terms are resuspension of activated material
      do k=top_lev,pver
         kp1=min(k+1,pver)
         km1=max(k-1,top_lev)
         q(k) = qold(k) + dt*( - src(k) + ekkp(k)*(qold(kp1) - qold(k) +       &
            qactold(kp1)*(1.0_r8-overlapp(k)))               &
            + ekkm(k)*(qold(km1) - qold(k) +     &
            qactold(km1)*(1.0_r8-overlapm(k))) )
         !        force to non-negative
         !        if(q(k)<-1.e-30)then
         !           write(iulog,*)'q=',q(k),' in explmix'
         q(k)=max(q(k),0._r8)
         !        endif
      end do

      !     diffusion loss at base of lowest layer
      q(pver)=q(pver)-surfrate*qold(pver)*dt+flxconv*dt
      !        force to non-negative
      !        if(q(pver)<-1.e-30)then
      !           write(iulog,*)'q=',q(pver),' in explmix'
      q(pver)=max(q(pver),0._r8)
      !        endif
   else
      do k=top_lev,pver
         kp1=min(k+1,pver)
         km1=max(k-1,top_lev)
         q(k) = qold(k) + dt*(src(k) + ekkp(k)*(overlapp(k)*qold(kp1)-qold(k)) +      &
            ekkm(k)*(overlapm(k)*qold(km1)-qold(k)) )
         !        force to non-negative
         !        if(q(k)<-1.e-30)then
         !           write(iulog,*)'q=',q(k),' in explmix'
         q(k)=max(q(k),0._r8)
         !        endif
      end do
      !     diffusion loss at base of lowest layer
      q(pver)=q(pver)-surfrate*qold(pver)*dt+flxconv*dt
      !        force to non-negative
      !        if(q(pver)<-1.e-30)then
      !           write(iulog,*)'q=',q(pver),' in explmix'
      q(pver)=max(q(pver),0._r8)

   end if

end subroutine explmix

!===============================================================================

subroutine ccncalc(aero_state, aero_props, ncol, pver, top_lev, tair, cs, ccn, errmsg, errflg)

   ! calculates number concentration of aerosols activated as CCN at
   ! supersaturation supersat.
   ! assumes an internal mixture of a multiple externally-mixed aerosol modes
   ! cgs units

   ! Ghan et al., Atmos. Res., 1993, 198-221.

   ! arguments
   class(aerosol_state), intent(in) :: aero_state
   class(aerosol_properties), intent(in) :: aero_props

   integer,  intent(in)  :: ncol        ! number of columns
   integer,  intent(in)  :: pver        ! number of vertical layers
   integer,  intent(in)  :: top_lev     ! top level for cloud physics
   real(r8), intent(in)  :: tair(:,:)   ! air temperature (K)
   real(r8), intent(in)  :: cs(:,:)     ! air density (kg/m3)
   real(r8), intent(out) :: ccn(:,:,:)  ! (ncol,pver,psat) number conc of aerosols activated at supersat (#/m3)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! local

   real(r8) naerosol(ncol,pver,nbin) ! interstit+activated aerosol number conc (/m3)
   real(r8) vaerosol(ncol,pver,nbin) ! interstit+activated aerosol volume conc (m3/m3)

   real(r8) amcube(ncol)
   real(r8), allocatable :: argfactor(:)
   real(r8) surften_coef
   real(r8) a(ncol) ! surface tension parameter
   real(r8) hygro(ncol,pver,nbin)  ! aerosol hygroscopicity
   real(r8) sm(ncol)  ! critical supersaturation at mode radius
   real(r8) arg(ncol)
   integer l,m,i,k, astat
   real(r8) smcoef(ncol)
   integer phase ! phase of aerosol

   integer :: errnum
   character(len=shr_kind_cs) :: errstr

   !     mathematical constants
   real(r8), parameter :: super(psat) = supersat(:psat)*0.01_r8
   real(r8), parameter :: smcoefcoef  = 2._r8/sqrt(27._r8)

   !-------------------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   allocate( argfactor(nbin), stat=astat )
   if (astat/=0) then
      errmsg = 'ndrop::ccncalc : not able to allocate argfactor'
      errflg = 1
      return
   end if

   surften_coef=2._r8*mwh2o*surften/(r_universal*rhoh2o)

   do m=1,nbin
      argfactor(m)=twothird/(sq2*aero_props%alogsig(m))
   end do

   phase=3 ! interstitial+cloudborne

   do m = 1, nbin
      call aero_state%loadaer( aero_props, &
           ncol, pver, &
           m, cs, phase, naerosol(:,:,m), vaerosol(:,:,m), &
           hygro(:,:,m), errnum, errstr)
      if (errnum/=0) then
         errmsg = 'ccncalc : '//trim(errstr)
         errflg = 1
         return
      end if
   end do

   ccn = 0._r8
   do k=top_lev,pver

      do i=1,ncol
         a(i)=surften_coef/tair(i,k)
         smcoef(i)=smcoefcoef*a(i)*sqrt(a(i))
      end do

      do m=1,nbin

         where(naerosol(:ncol,k,m)>1.e-3_r8 .and. hygro(:ncol,k,m)>1.e-10_r8)
            amcube(:ncol)=aero_props%amcube(m, vaerosol(:ncol,k,m), naerosol(:ncol,k,m) )
            sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol,k,m)*amcube(:ncol)) ! critical supersaturation
         elsewhere
            sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
         endwhere
         do l=1,psat
            do i=1,ncol
               arg(i)=argfactor(m)*log(sm(i)/super(l))
               ccn(i,k,l)=ccn(i,k,l)+naerosol(i,k,m)*0.5_r8*(1._r8-erf(arg(i)))
            enddo
         enddo
      enddo
   enddo
   ccn(:ncol,:,:)=ccn(:ncol,:,:)*1.e-6_r8 ! convert from #/m3 to #/cm3

   deallocate( argfactor )

end subroutine ccncalc

!===============================================================================
end module ndrop
