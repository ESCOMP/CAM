module wetdep_cam

!-----------------------------------------------------------------------
!
! CAM host interface for the portable wetdep module
! Calls portable wetdepa.
! Retains gas-phase Henry's law scavenging (wetdepg).
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, pver
   use physconst,    only: gravit, rair, tmelt
   use phys_control, only: cam_physpkg_is

   use wetdep,       only: clddiag

   implicit none
   private

   public :: wetdepg     ! scavenging of gas phase constituents by henry's law
   public :: wetdep_inputs_t
   public :: wetdep_init
   public :: wetdep_inputs_set

   real(r8), parameter :: cmftau = 3600._r8
   real(r8), parameter :: rhoh2o = 1000._r8            ! density of water
   real(r8), parameter :: molwta = 28.97_r8            ! molecular weight dry air gm/mole
   real(r8), parameter :: omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero

   type wetdep_inputs_t
      real(r8), pointer :: cldt(:,:) => null()  ! cloud fraction
      real(r8), pointer :: prain(:,:) => null()
      real(r8), pointer :: bergso(:,:) => null()
      real(r8), pointer :: evapr(:,:) => null()
      real(r8) :: cldcu(pcols,pver)     ! convective cloud fraction, currently empty
      real(r8) :: evapc(pcols,pver)     ! Evaporation rate of convective precipitation
      real(r8) :: cmfdqr(pcols,pver)    ! convective production of rain
      real(r8) :: conicw(pcols,pver)    ! convective in-cloud water
      real(r8) :: totcond(pcols, pver)  ! total condensate
      real(r8) :: cldv(pcols,pver)      ! cloudy volume undergoing wet chem and scavenging
      real(r8) :: cldvcu(pcols,pver)    ! Convective precipitation area at the top interface of current layer
      real(r8) :: cldvst(pcols,pver)    ! Stratiform precipitation area at the top interface of current layer
   end type wetdep_inputs_t

   integer :: cld_idx             = 0
   integer :: prain_idx           = 0
   integer :: bergso_idx          = 0
   integer :: nevapr_idx          = 0

   integer :: icwmrdp_idx         = 0
   integer :: icwmrsh_idx         = 0
   integer :: rprddp_idx          = 0
   integer :: rprdsh_idx          = 0
   integer :: sh_frac_idx         = 0
   integer :: dp_frac_idx         = 0
   integer :: nevapr_shcu_idx     = 0
   integer :: nevapr_dpcu_idx     = 0
   integer :: ixcldice, ixcldliq

contains

   subroutine wetdep_init()
      use physics_buffer, only: pbuf_get_index
      use constituents,   only: cnst_get_ind

      integer :: ierr

      cld_idx             = pbuf_get_index('CLD')
      prain_idx           = pbuf_get_index('PRAIN')
      bergso_idx          = pbuf_get_index('BERGSO', errcode=ierr )
      nevapr_idx          = pbuf_get_index('NEVAPR')

      icwmrdp_idx         = pbuf_get_index('ICWMRDP')
      rprddp_idx          = pbuf_get_index('RPRDDP')
      icwmrsh_idx         = pbuf_get_index('ICWMRSH')
      rprdsh_idx          = pbuf_get_index('RPRDSH')
      sh_frac_idx         = pbuf_get_index('SH_FRAC' )
      dp_frac_idx         = pbuf_get_index('DP_FRAC')
      nevapr_shcu_idx     = pbuf_get_index('NEVAPR_SHCU')
      nevapr_dpcu_idx     = pbuf_get_index('NEVAPR_DPCU')

      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('CLDLIQ', ixcldliq)

   end subroutine wetdep_init

!==============================================================================
! gathers up the inputs needed for the wetdepa routines
!==============================================================================
   subroutine wetdep_inputs_set( state, pbuf, inputs )
      use physics_types,  only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

      ! args

      type(physics_state),  intent(in )  :: state           !! physics state
      type(physics_buffer_desc), pointer :: pbuf(:)         !! physics buffer
      type(wetdep_inputs_t), intent(out) :: inputs          !! collection of wetdepa inputs

      ! local vars

      real(r8), pointer :: icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
      real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
      real(r8), pointer :: icwmrsh(:,:)    ! in cloud water mixing ratio, deep convection
      real(r8), pointer :: rprdsh(:,:)     ! rain production, deep convection
      real(r8), pointer :: sh_frac(:,:)    ! Shallow convective cloud fraction
      real(r8), pointer :: dp_frac(:,:)    ! Deep convective cloud fraction
      real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
      real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.

      real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume
      real(r8) :: cldst(pcols,pver)        ! Stratiform cloud fraction

      integer :: itim, ncol

      ncol = state%ncol
      itim = pbuf_old_tim_idx()

      call pbuf_get_field(pbuf, cld_idx,         inputs%cldt, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, prain_idx,       inputs%prain   )
      call pbuf_get_field(pbuf, nevapr_idx,      inputs%evapr   )
      call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp )
      call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh )
      call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
      call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
      call pbuf_get_field(pbuf, sh_frac_idx,     sh_frac )
      call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac )
      call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
      call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )

      if (bergso_idx>0) then
         call pbuf_get_field(pbuf, bergso_idx, inputs%bergso )
      else
         if (.not. associated(inputs%bergso)) then
            allocate(inputs%bergso(pcols,pver))
            inputs%bergso(:,:) = 0.0_r8
         endif
      endif

      inputs%cldcu(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)
      cldst(:ncol,:)          = inputs%cldt(:ncol,:) - inputs%cldcu(:ncol,:)       ! Stratiform cloud fraction
      inputs%evapc(:ncol,:)  = evapcsh(:ncol,:) + evapcdp(:ncol,:)
      inputs%cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)

      ! sum deep and shallow convection contributions
      if (cam_physpkg_is('cam5') .or. cam_physpkg_is('cam6')) then
         ! Dec.29.2009. Sungsu
         inputs%conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                                  max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
      else
         inputs%conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
      end if

      inputs%totcond(:ncol,:) = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)

      call clddiag( state%t,     state%pmid,   state%pdel,   inputs%cmfdqr, inputs%evapc, &
                   inputs%cldt,  inputs%cldcu,       cldst,  inputs%evapr, &
                   inputs%prain, inputs%cldv, inputs%cldvcu, inputs%cldvst,       rainmr, &
                   state%ncol, pver, gravit, tmelt, rair )

   end subroutine wetdep_inputs_set

!==============================================================================

! wetdepg is currently being used for both CAM4 and CAM5 by making use of the
! cam_physpkg_is method.

   subroutine wetdepg( t, p, q, pdel, &
                      cldt, cldc, cmfdqr, evapc, precs, evaps, &
                      rain, cwat, tracer, deltat, molwt, &
                      solconst, scavt, iscavt, cldv, icwmr1, &
                      icwmr2, fracis, ncol )

      !-----------------------------------------------------------------------
      ! Purpose:
      ! scavenging of gas phase constituents by henry's law
      !
      ! Author: P. Rasch
      !-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         t(pcols,pver),        &! temperature
         p(pcols,pver),        &! pressure
         q(pcols,pver),        &! moisture
         pdel(pcols,pver),     &! pressure thikness
         cldt(pcols,pver),     &! total cloud fraction
         cldc(pcols,pver),     &! convective cloud fraction
         cmfdqr(pcols,pver),   &! rate of production of convective precip
         rain (pcols,pver),    &! total rainwater mixing ratio
         cwat(pcols,pver),     &! cloud water amount
         precs(pcols,pver),    &! rate of production of stratiform precip
         evaps(pcols,pver),    &! rate of evaporation of precip
         ! Sungsu
         evapc(pcols,pver),    &! Rate of evaporation of convective precipitation
         ! Sungsu
         cldv(pcols,pver),     &! estimate of local volume occupied by clouds
         icwmr1 (pcols,pver),  &! in cloud water mixing ration for zhang scheme
         icwmr2 (pcols,pver),  &! in cloud water mixing ration for hack  scheme
         deltat,               &! time step
         tracer(pcols,pver),   &! trace species
         molwt                  ! molecular weights

      integer, intent(in) :: ncol

      real(r8) &
         solconst(pcols,pver)   ! Henry's law coefficient

      real(r8), intent(out) ::&
         scavt(pcols,pver),    &! scavenging tend
         iscavt(pcols,pver),   &! incloud scavenging tends
         fracis(pcols, pver)    ! fraction of constituent that is insoluble

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cwatc                ! local convective total water amount
      real(r8) cwats                ! local stratiform total water amount
      real(r8) cwatl                ! local cloud liq water amount
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) cwatpl               ! local water amount falling from above precip (liq)
      real(r8) cwatt                ! local sum of strat + conv total water amount
      real(r8) cwatti               ! cwatt/cldv = cloudy grid volume mixing ratio
      real(r8) fracev               ! fraction of precip from above that is evaporating
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
      ! in terms of mixing ratios
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/gravit)
      real(r8) precab(pcols)        ! precip from above (work array)
      real(r8) precbl               ! precip work variable
      real(r8) precxx               ! precip work variable
      real(r8) precxx2               !
      real(r8) precic               ! precip work variable
      real(r8) rat                  ! ratio of amount available to amount removed
      real(r8) scavab(pcols)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(pcols)       ! scavenged tracer flux from above (work array)

      real(r8) scavmax              ! an estimate of the max tracer avail for removal
      real(r8) scavbl               ! flux removed at bottom of layer
      real(r8) fins                 ! in cloud fraction removed by strat rain
      real(r8) finc                 ! in cloud fraction removed by conv rain
      real(r8) rate                 ! max removal rate estimate
      real(r8) scavlimt             ! limiting value 1
      real(r8) scavt1               ! limiting value 2
      real(r8) scavin               ! scavenging by incloud processes
      real(r8) scavbc               ! scavenging by below cloud processes
      real(r8) tc
      real(r8) weight               ! ice fraction
      real(r8) wtpl                 ! work variable
      real(r8) cldmabs(pcols)       ! maximum cloud at or above this level
      real(r8) cldmabc(pcols)       ! maximum cloud at or above this level
      !-----------------------------------------------------------

      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! zero accumulators
      do i = 1,pcols
         precab(i) = 1.e-36_r8
         scavab(i) = 0._r8
         cldmabs(i) = 0._r8
      end do

      do k = 1,pver
         do i = 1,ncol

            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

            cldmabs(i) = max(cldmabs(i),cldt(i,k))

            ! partitioning coefs for gas and aqueous phase
            !              take as a cloud water amount, the sum of the stratiform amount
            !              plus the convective rain water amount

            ! convective amnt is just the local precip rate from the hack scheme
            !              since there is no storage of water, this ignores that falling from above
            !            cwatc = cmfdqr(i,k)*deltat/adjfac
            !++mcb -- test cwatc
            cwatc = (icwmr1(i,k) + icwmr2(i,k)) * (1._r8-weight)
            !--mcb

            ! strat cloud water amount and also ignore the part falling from above
            cwats = cwat(i,k)

            ! cloud water as liq
            !++mcb -- add cwatc later (in cwatti)
            !            cwatl = (1.-weight)*(cwatc+cwats)
            cwatl = (1._r8-weight)*cwats
            ! cloud water as ice
            !*not used        cwati = weight*(cwatc+cwats)

            ! total suspended condensate as liquid
            cwatt = cwatl + rain(i,k)

            ! incloud version
            !++mcb -- add cwatc here
            cwatti = cwatt/max(cldv(i,k), 0.00001_r8) + cwatc

            ! partitioning terms
            patm = p(i,k)/1.013e5_r8 ! pressure in atmospheres
            hconst = molwta*patm*solconst(i,k)*cwatti/rhoh2o
            aqfrac = hconst/(1._r8+hconst)
            gafrac = 1/(1._r8+hconst)
            fracis(i,k) = gafrac

            ! partial pressure of the tracer in the gridbox in atmospheres
            part = patm*gafrac*tracer(i,k)*molwta/molwt

            ! use henrys law to give moles tracer /liter of water
            ! in this volume
            ! then convert to kg tracer /liter of water (kg tracer / kg water)
            mplb = solconst(i,k)*part*molwt/1000._r8

            pdog = pdel(i,k)/gravit

            ! this part of precip will be carried downward but at a new molarity of mpl
            precic = pdog*(precs(i,k) + cmfdqr(i,k))

            ! we cant take out more than entered, plus that available in the cloud
            !                  scavmax = scavab(i)+tracer(i,k)*cldt(i,k)/deltat*pdog
            scavmax = scavab(i)+tracer(i,k)*cldv(i,k)/deltat*pdog

            ! flux of tracer by incloud processes
            scavin = precic*(1._r8-weight)*mplb

            ! fraction of precip which entered above that leaves below
            if (cam_physpkg_is('cam5') .or. cam_physpkg_is('cam6')) then
               ! Sungsu added evaporation of convective precipitation below.
               precxx = precab(i)-pdog*(evaps(i,k)+evapc(i,k))
            else
               precxx = precab(i)-pdog*evaps(i,k)
            end if
            precxx = max (precxx,0.0_r8)

            ! flux of tracer by below cloud processes
            !++mcb -- removed wtpl because it is now not assigned and previously
            !          when it was assigned it was unnecessary:  if(tc.gt.0)wtpl=1
            if (tc>0.0_r8) then
               !               scavbc = precxx*wtpl*mplb ! if liquid
               scavbc = precxx*mplb ! if liquid
            else
               precxx2=max(precxx,1.e-36_r8)
               scavbc = scavab(i)*precxx2/(precab(i)) ! if ice
            endif

            scavbl = min(scavbc + scavin, scavmax)

            ! first guess assuming that henries law works
            scavt1 = (scavab(i)-scavbl)/pdog*omsm

            ! pjr this should not be required, but we put it in to make sure we cant remove too much
            ! remember, scavt1 is generally negative (indicating removal)
            scavt1 = max(scavt1,-tracer(i,k)*cldv(i,k)/deltat)

            !++mcb -- remove this limitation for gas species
            !c use the dana and hales or balkanski limit on scavenging
            !c            rate = precab(i)*0.1
            !            rate = (precic + precxx)*0.1
            !            scavlimt = -tracer(i,k)*cldv(i,k)
            !     $           *rate/(1.+rate*deltat)

            !            scavt(i,k) = max(scavt1, scavlimt)

            ! instead just set scavt to scavt1
            scavt(i,k) = scavt1
            !--mcb

            ! now update the amount leaving the layer
            scavbl = scavab(i) - scavt(i,k)*pdog

            ! in cloud amount is that formed locally over the total flux out bottom
            fins = scavin/(scavin + scavbc + 1.e-36_r8)
            iscavt(i,k) = scavt(i,k)*fins

            scavab(i) = scavbl
            precab(i) = max(precxx + precic,1.e-36_r8)

         end do
      end do

   end subroutine wetdepg

end module wetdep_cam
