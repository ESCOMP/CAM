module  zm_microphysics

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for cumulus microphysics
! 
! Author: Xialiang Song and Guang Jun Zhang, June 2010  
!---------------------------------------------------------------------------------

use shr_kind_mod,      only: r8=>shr_kind_r8
use spmd_utils,        only: masterproc    
use ppgrid,            only: pcols, pver, pverp
use physconst,         only: gravit, rair, tmelt, cpair, rh2o, r_universal, mwh2o, rhoh2o
use physconst,         only: latvap, latice
!use activate_drop_mam, only: actdrop_mam_calc
use ndrop,             only: activate_modal
use ndrop_bam,         only: ndrop_bam_run
use nucleate_ice,      only: nucleati
use shr_spfn_mod,     only: erf => shr_spfn_erf 
use shr_spfn_mod,     only: gamma => shr_spfn_gamma
use wv_saturation,  only: svp_water, svp_ice
use cam_logfile,       only: iulog
use cam_abortutils,        only: endrun
use micro_mg_utils, only:ice_autoconversion, snow_self_aggregation, accrete_cloud_water_snow, &
                         secondary_ice_production, accrete_rain_snow, heterogeneous_rain_freezing, &
                         accrete_cloud_water_rain, self_collection_rain, accrete_cloud_ice_snow

implicit none
private
save

public :: &
   zm_mphyi, &
   zm_mphy,  &
   zm_conv_t,&
   zm_aero_t

! Private module data

! constants remaped
real(r8) :: g      ! gravity
real(r8) :: mw     ! molecular weight of water
real(r8) :: r      ! Dry air Gas constant
real(r8) :: rv     ! water vapor gas contstant
real(r8) :: rr     ! universal gas constant
real(r8) :: cpp    ! specific heat of dry air
real(r8) :: rhow   ! density of liquid water
real(r8) :: xlf    ! latent heat of freezing

!from 'microconstants'
real(r8) :: rhosn  ! bulk density snow
real(r8) :: rhoi   ! bulk density ice

real(r8) :: ac,bc,as,bs,ai,bi,ar,br  !fall speed parameters 
real(r8) :: ci,di    !ice mass-diameter relation parameters
real(r8) :: cs,ds    !snow mass-diameter relation parameters
real(r8) :: cr,dr    !drop mass-diameter relation parameters
real(r8) :: Eii      !collection efficiency aggregation of ice
real(r8) :: Ecc      !collection efficiency
real(r8) :: Ecr      !collection efficiency cloud droplets/rain
real(r8) :: DCS      !autoconversion size threshold
real(r8) :: bimm,aimm !immersion freezing
real(r8) :: rhosu     !typical 850mn air density
real(r8) :: mi0       ! new crystal mass
real(r8) :: rin       ! radius of contact nuclei
real(r8) :: pi       ! pi

! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8


type, public :: ptr2d
   real(r8), pointer :: val(:,:)
end type ptr2d

! Aerosols
type :: zm_aero_t

   ! Aerosol treatment
   character(len=5) :: scheme  ! either 'bulk' or 'modal'

   ! Bulk aerosols
   integer :: nbulk    =  0 ! number of bulk aerosols affecting climate
   integer :: idxsul   = -1 ! index in aerosol list for sulfate
   integer :: idxdst1  = -1 ! index in aerosol list for dust1
   integer :: idxdst2  = -1 ! index in aerosol list for dust2
   integer :: idxdst3  = -1 ! index in aerosol list for dust3
   integer :: idxdst4  = -1 ! index in aerosol list for dust4
   integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHI)

   real(r8),    allocatable :: num_to_mass_aer(:)  ! conversion of mmr to number conc for bulk aerosols
   type(ptr2d), allocatable :: mmr_bulk(:)         ! array of pointers to bulk aerosol mmr
   real(r8),    allocatable :: mmrg_bulk(:,:,:)    ! gathered bulk aerosol mmr

   ! Modal aerosols
   integer                  :: nmodes = 0      ! number of modes
   integer,     allocatable :: nspec(:)        ! number of species in each mode
   type(ptr2d), allocatable :: num_a(:)        ! number mixing ratio of modes (interstitial phase)
   type(ptr2d), allocatable :: mmr_a(:,:)      ! species mmr in each mode (interstitial phase)
   real(r8),    allocatable :: numg_a(:,:,:)   ! gathered number mixing ratio of modes (interstitial phase)
   real(r8),    allocatable :: mmrg_a(:,:,:,:) ! gathered species mmr in each mode (interstitial phase)
   real(r8),    allocatable :: voltonumblo(:)  ! volume to number conversion (lower bound) for each mode
   real(r8),    allocatable :: voltonumbhi(:)  ! volume to number conversion (upper bound) for each mode
   real(r8),    allocatable :: specdens(:,:)   ! density of modal species
   real(r8),    allocatable :: spechygro(:,:)  ! hygroscopicity of modal species

   integer :: mode_accum_idx  = -1  ! index of accumulation mode
   integer :: mode_aitken_idx = -1  ! index of aitken mode
   integer :: mode_coarse_idx = -1  ! index of coarse mode
   integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
   integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

   type(ptr2d), allocatable :: dgnum(:)        ! mode dry radius
   real(r8),    allocatable :: dgnumg(:,:,:)   ! gathered mode dry radius

   real(r8) :: sigmag_aitken

end type zm_aero_t

type :: zm_conv_t

   real(r8),    allocatable :: qi(:,:)             ! wg grid slice of cloud ice.
   real(r8),    allocatable :: qliq(:,:)           ! convective cloud liquid water.
   real(r8),    allocatable :: qice(:,:)           ! convective cloud ice.
   real(r8),    allocatable :: wu(:,:)             ! vertical velocity
   real(r8),    allocatable :: sprd(:,:)           ! rate of production of snow at that layer
   real(r8),    allocatable :: qrain(:,:)          ! convective rain water.
   real(r8),    allocatable :: qsnow(:,:)          ! convective snow.
   real(r8),    allocatable :: qnl(:,:)            ! convective cloud liquid water num concen.
   real(r8),    allocatable :: qni(:,:)            ! convective cloud ice num concen.
   real(r8),    allocatable :: qnr(:,:)            ! convective rain water num concen.
   real(r8),    allocatable :: qns(:,:)            ! convective snow num concen.
   real(r8),    allocatable :: frz(:,:)            ! heating rate due to freezing
   real(r8),    allocatable :: autolm(:,:)    !mass tendency due to autoconversion of droplets to rain
   real(r8),    allocatable :: accrlm(:,:)    !mass tendency due to accretion of droplets by rain
   real(r8),    allocatable :: bergnm(:,:)    !mass tendency due to Bergeron process
   real(r8),    allocatable :: fhtimm(:,:)    !mass tendency due to immersion freezing
   real(r8),    allocatable :: fhtctm(:,:)    !mass tendency due to contact freezing
   real(r8),    allocatable :: fhmlm (:,:)    !mass tendency due to homogeneous freezing
   real(r8),    allocatable :: hmpim (:,:)    !mass tendency due to HM process
   real(r8),    allocatable :: accslm(:,:)    !mass tendency due to accretion of droplets by snow
   real(r8),    allocatable :: dlfm  (:,:)    !mass tendency due to detrainment of droplet
   real(r8),    allocatable :: autoln(:,:)    !num tendency due to autoconversion of droplets to rain
   real(r8),    allocatable :: accrln(:,:)    !num tendency due to accretion of droplets by rain
   real(r8),    allocatable :: bergnn(:,:)    !num tendency due to Bergeron process
   real(r8),    allocatable :: fhtimn(:,:)    !num tendency due to immersion freezing
   real(r8),    allocatable :: fhtctn(:,:)    !num tendency due to contact freezing
   real(r8),    allocatable :: fhmln (:,:)    !num tendency due to homogeneous freezing
   real(r8),    allocatable :: accsln(:,:)    !num tendency due to accretion of droplets by snow
   real(r8),    allocatable :: activn(:,:)    !num tendency due to droplets activation
   real(r8),    allocatable :: dlfn  (:,:)    !num tendency due to detrainment of droplet
   real(r8),    allocatable :: autoim(:,:)    !mass tendency due to autoconversion of cloud ice to snow
   real(r8),    allocatable :: accsim(:,:)    !mass tendency due to accretion of cloud ice by snow
   real(r8),    allocatable :: difm  (:,:)    !mass tendency due to detrainment of cloud ice
   real(r8),    allocatable :: nuclin(:,:)    !num tendency due to ice nucleation
   real(r8),    allocatable :: autoin(:,:)    !num tendency due to autoconversion of cloud ice to snow
   real(r8),    allocatable :: accsin(:,:)    !num tendency due to accretion of cloud ice by snow
   real(r8),    allocatable :: hmpin (:,:)    !num tendency due to HM process
   real(r8),    allocatable :: difn  (:,:)    !num tendency due to detrainment of cloud ice
   real(r8),    allocatable :: cmel  (:,:)    !mass tendency due to condensation
   real(r8),    allocatable :: cmei  (:,:)    !mass tendency due to deposition
   real(r8),    allocatable :: trspcm(:,:)    !LWC tendency due to convective transport
   real(r8),    allocatable :: trspcn(:,:)    !droplet num tendency due to convective transport
   real(r8),    allocatable :: trspim(:,:)    !IWC tendency due to convective transport
   real(r8),    allocatable :: trspin(:,:)    !ice crystal num tendency due to convective transport
   real(r8),    allocatable :: dcape(:)       ! CAPE change due to freezing heating
   real(r8),    allocatable :: lambdadpcu(:,:)! slope of cloud liquid size distr
   real(r8),    allocatable :: mudpcu(:,:)    ! width parameter of droplet size distr
   real(r8),    allocatable :: di(:,:)
   real(r8),    allocatable :: dnl(:,:)
   real(r8),    allocatable :: dni(:,:)
   real(r8),    allocatable :: qide(:,:)     ! cloud ice mixing ratio for detrainment (kg/kg)
   real(r8),    allocatable :: qncde(:,:)    ! cloud water number concentration for detrainment (1/kg)
   real(r8),    allocatable :: qnide(:,:)    ! cloud ice number concentration for detrainment (1/kg)


end type zm_conv_t

real(r8), parameter :: dcon  = 25.e-6_r8
real(r8), parameter :: mucon = 5.3_r8
real(r8), parameter :: lambdadpcu = (mucon + 1._r8)/dcon

!===============================================================================
contains
!===============================================================================

subroutine zm_mphyi

!----------------------------------------------------------------------- 
! 
! Purpose:
! initialize constants for the cumulus microphysics
! called from zm_conv_init() in zm_conv_intr.F90
!
! Author: Xialiang Song, June 2010
! 
!-----------------------------------------------------------------------

!NOTE:
! latent heats should probably be fixed with temperature 
! for energy conservation with the rest of the model
! (this looks like a +/- 3 or 4% effect, but will mess up energy balance)

   xlf = latice          ! latent heat freezing

! from microconstants

! parameters below from Reisner et al. (1998)
! density parameters (kg/m3)

      rhosn = 100._r8    ! bulk density snow
      rhoi = 500._r8     ! bulk density ice
      rhow = 1000._r8    ! bulk density liquid

! fall speed parameters, V = aD^b
! V is in m/s

! droplets
	ac = 3.e7_r8
	bc = 2._r8

! snow
	as = 11.72_r8
	bs = 0.41_r8

! cloud ice
	ai = 700._r8
	bi = 1._r8

! rain
	ar = 841.99667_r8
	br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d

        pi= 3.14159265358979323846_r8     

! cloud ice mass-diameter relationship

	ci = rhoi*pi/6._r8
	di = 3._r8

! snow mass-diameter relationship

	cs = rhosn*pi/6._r8
	ds = 3._r8

! drop mass-diameter relationship

	cr = rhow*pi/6._r8
	dr = 3._r8

! collection efficiency, aggregation of cloud ice and snow

	Eii = 0.1_r8

! collection efficiency, accretion of cloud water by rain

	Ecr = 1.0_r8

! autoconversion size threshold for cloud ice to snow (m)

        Dcs = 150.e-6_r8
! immersion freezing parameters, bigg 1953

	bimm = 100._r8
	aimm = 0.66_r8

! typical air density at 850 mb

	rhosu = 85000._r8/(rair * tmelt)

! mass of new crystal due to aerosol freezing and growth (kg)

	mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)

! radius of contact nuclei aerosol (m)

        rin = 0.1e-6_r8

end subroutine zm_mphyi

!===============================================================================

subroutine zm_mphy(su,    qu,   mu,   du,   eu,    cmel,  cmei,  zf,   pm,   te,   qe,        &
                   eps0,  jb,   jt,   jlcl, msg,   il2g,  grav,  cp,   rd,   aero, gamhat,    &
                   qc,    qi,   nc,   ni,   qcde,  qide,  ncde,  nide, rprd, sprd, frz,       &
                   wu,    qr,   qni,  nr,   ns,    autolm, accrlm, bergnm, fhtimm, fhtctm,    &
                   fhmlm,  hmpim,  accslm, dlfm,   autoln, accrln, bergnn, fhtimn, fhtctn,    &
                   fhmln,  accsln, activn, dlfn,   autoim, accsim, difm,   nuclin, autoin,    &
                   accsin, hmpin,  difn,   trspcm, trspcn, trspim, trspin, lamc,   pgam  )
   

! Purpose:
! microphysic parameterization for Zhang-McFarlane convection scheme
! called from cldprp() in zm_conv.F90
!
! Author: Xialiang Song, June 2010

  use time_manager,    only: get_step_size

! variable declarations

  implicit none

! input variables
  real(r8), intent(in) :: su(pcols,pver)        ! normalized dry stat energy of updraft
  real(r8), intent(in) :: qu(pcols,pver)        ! spec hum of updraft
  real(r8), intent(in) :: mu(pcols,pver)        ! updraft mass flux
  real(r8), intent(in) :: du(pcols,pver)        ! detrainement rate of updraft
  real(r8), intent(in) :: eu(pcols,pver)        ! entrainment rate of updraft
  real(r8), intent(in) :: cmel(pcols,pver)      ! condensation rate of updraft
  real(r8), intent(in) :: cmei(pcols,pver)      ! condensation rate of updraft
  real(r8), intent(in) :: zf(pcols,pverp)       ! height of interfaces
  real(r8), intent(in) :: pm(pcols,pver)        ! pressure of env
  real(r8), intent(in) :: te(pcols,pver)        ! temp of env
  real(r8), intent(in) :: qe(pcols,pver)        ! spec. humidity of env
  real(r8), intent(in) :: eps0(pcols)
  real(r8), intent(in) :: gamhat(pcols,pver)    ! gamma=L/cp(dq*/dT) at interface

  integer, intent(in) :: jb(pcols)              ! updraft base level
  integer, intent(in) :: jt(pcols)              ! updraft plume top
  integer, intent(in) :: jlcl(pcols)            ! updraft lifting cond level
  integer, intent(in) :: msg                    ! missing moisture vals
  integer, intent(in) :: il2g                   ! number of columns in gathered arrays

  type(zm_aero_t), intent(in) :: aero           ! aerosol object

  real(r8) grav                                 ! gravity
  real(r8) cp                                   ! heat capacity of dry air
  real(r8) rd                                   ! gas constant for dry air

! output variables
  real(r8), intent(out) :: qc(pcols,pver)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(out) :: qi(pcols,pver)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(out) :: nc(pcols,pver)       ! cloud water number conc (1/kg)
  real(r8), intent(out) :: ni(pcols,pver)       ! cloud ice number conc (1/kg)
  real(r8), intent(out) :: qcde(pcols,pver)     ! cloud water mixing ratio for detrainment(kg/kg)
  real(r8), intent(out) :: qide(pcols,pver)     ! cloud ice mixing ratio for detrainment (kg/kg)
  real(r8), intent(out) :: ncde(pcols,pver)     ! cloud water number conc for detrainment (1/kg)
  real(r8), intent(out) :: nide(pcols,pver)     ! cloud ice number conc for detrainment (1/kg)
  real(r8), intent(out) :: wu(pcols,pver)
  real(r8), intent(out) :: qni(pcols,pver)      ! snow mixing ratio
  real(r8), intent(out) :: qr(pcols,pver)       ! rain mixing ratio
  real(r8), intent(out) :: ns(pcols,pver)       ! snow number conc
  real(r8), intent(out) :: nr(pcols,pver)       ! rain number conc
  real(r8), intent(out) :: rprd(pcols,pver)     ! rate of production of precip at that layer
  real(r8), intent(out) :: sprd(pcols,pver)     ! rate of production of snow at that layer
  real(r8), intent(out) :: frz(pcols,pver)      ! rate of freezing 


  real(r8), intent(inout) :: lamc(pcols,pver)     ! slope of cloud liquid size distr
  real(r8), intent(inout) :: pgam(pcols,pver)     ! spectral width parameter of droplet size distr

! tendency for output
  real(r8),intent(out) :: autolm(pcols,pver)    !mass tendency due to autoconversion of droplets to rain
  real(r8),intent(out) :: accrlm(pcols,pver)    !mass tendency due to accretion of droplets by rain
  real(r8),intent(out) :: bergnm(pcols,pver)    !mass tendency due to Bergeron process
  real(r8),intent(out) :: fhtimm(pcols,pver)    !mass tendency due to immersion freezing
  real(r8),intent(out) :: fhtctm(pcols,pver)    !mass tendency due to contact freezing
  real(r8),intent(out) :: fhmlm (pcols,pver)    !mass tendency due to homogeneous freezing
  real(r8),intent(out) :: hmpim (pcols,pver)    !mass tendency due to HM process
  real(r8),intent(out) :: accslm(pcols,pver)    !mass tendency due to accretion of droplets by snow
  real(r8),intent(out) :: dlfm  (pcols,pver)    !mass tendency due to detrainment of droplet
  real(r8),intent(out) :: trspcm(pcols,pver)    !mass tendency of droplets due to convective transport

  real(r8),intent(out) :: autoln(pcols,pver)    !num tendency due to autoconversion of droplets to rain
  real(r8),intent(out) :: accrln(pcols,pver)    !num tendency due to accretion of droplets by rain
  real(r8),intent(out) :: bergnn(pcols,pver)    !num tendency due to Bergeron process
  real(r8),intent(out) :: fhtimn(pcols,pver)    !num tendency due to immersion freezing
  real(r8),intent(out) :: fhtctn(pcols,pver)    !num tendency due to contact freezing
  real(r8),intent(out) :: fhmln (pcols,pver)    !num tendency due to homogeneous freezing
  real(r8),intent(out) :: accsln(pcols,pver)    !num tendency due to accretion of droplets by snow
  real(r8),intent(out) :: activn(pcols,pver)    !num tendency due to droplets activation
  real(r8),intent(out) :: dlfn  (pcols,pver)    !num tendency due to detrainment of droplet
  real(r8),intent(out) :: trspcn(pcols,pver)    !num tendency of droplets due to convective transport

  real(r8),intent(out) :: autoim(pcols,pver)    !mass tendency due to autoconversion of cloud ice to snow
  real(r8),intent(out) :: accsim(pcols,pver)    !mass tendency due to accretion of cloud ice by snow
  real(r8),intent(out) :: difm  (pcols,pver)    !mass tendency due to detrainment of cloud ice
  real(r8),intent(out) :: trspim(pcols,pver)    !mass tendency of ice crystal due to convective transport

  real(r8),intent(out) :: nuclin(pcols,pver)    !num tendency due to ice nucleation
  real(r8),intent(out) :: autoin(pcols,pver)    !num tendency due to autoconversion of cloud ice to snow
  real(r8),intent(out) :: accsin(pcols,pver)    !num tendency due to accretion of cloud ice by snow
  real(r8),intent(out) :: hmpin (pcols,pver)    !num tendency due to HM process
  real(r8),intent(out) :: difn  (pcols,pver)    !num tendency due to detrainment of cloud ice
  real(r8),intent(out) :: trspin(pcols,pver)    !num tendency of ice crystal due to convective transport

!................................................................................
! local workspace
! all units mks unless otherwise stated
  real(r8) :: deltat                ! time step (s)
  real(r8) :: omsm                  ! number near unity for round-off issues
  real(r8) :: dum                   ! temporary dummy variable
  real(r8) :: dum1                  ! temporary dummy variable 
  real(r8) :: dum2                  ! temporary dummy variable

  real(r8) :: q(pcols,pver)         ! water vapor mixing ratio (kg/kg)
  real(r8) :: t(pcols,pver)         ! temperature (K)
  real(r8) :: rho(pcols,pver)       ! air density (kg m-3)
  real(r8) :: dz(pcols,pver)        ! height difference across model vertical level

  real(r8) :: qcic(pcols,pver)      ! in-cloud cloud liquid mixing ratio
  real(r8) :: qiic(pcols,pver)      ! in-cloud cloud ice mixing ratio
  real(r8) :: qniic(pcols,pver)     ! in-precip snow mixing ratio
  real(r8) :: qric(pcols,pver)      ! in-precip rain mixing ratio
  real(r8) :: ncic(pcols,pver)      ! in-cloud droplet number conc
  real(r8) :: niic(pcols,pver)      ! in-cloud cloud ice number conc
  real(r8) :: nsic(pcols,pver)      ! in-precip snow number conc
  real(r8) :: nric(pcols,pver)      ! in-precip rain number conc

  real(r8) :: lami(pver)            ! slope of cloud ice size distr
  real(r8) :: n0i(pver)             ! intercept of cloud ice size distr
  real(r8) :: n0c(pver)             ! intercept of cloud liquid size distr
  real(r8) :: lams(pver)            ! slope of snow size distr
  real(r8) :: n0s(pver)             ! intercept of snow size distr
  real(r8) :: lamr(pver)            ! slope of rain size distr
  real(r8) :: n0r(pver)             ! intercept of rain size distr
  real(r8) :: cdist1(pver)          ! size distr parameter to calculate droplet freezing
  real(r8) :: lammax                ! maximum allowed slope of size distr
  real(r8) :: lammin                ! minimum allowed slope of size distr

  real(r8) :: mnuccc(pver)          ! mixing ratio tendency due to freezing of cloud water
  real(r8) :: nnuccc(pver)          ! number conc tendency due to freezing of cloud water
  real(r8) :: mnucct(pver)          ! mixing ratio tendency due to contact freezing of cloud water
  real(r8) :: nnucct(pver)          ! number conc tendency due to contact freezing of cloud water
  real(r8) :: msacwi(pver)          ! mixing ratio tendency due to HM ice multiplication
  real(r8) :: nsacwi(pver)          ! number conc tendency due to HM ice multiplication
  real(r8) :: prf(pver)             ! mixing ratio tendency due to fallout of rain
  real(r8) :: psf(pver)             ! mixing ratio tendency due to fallout of snow
  real(r8) :: pnrf(pver)            ! number conc tendency due to fallout of rain
  real(r8) :: pnsf(pver)            ! number conc tendency due to fallout of snow
  real(r8) :: prc(pver)             ! mixing ratio tendency due to autoconversion of cloud droplets
  real(r8) :: nprc(pver)            ! number conc tendency due to autoconversion of cloud droplets
  real(r8) :: nprc1(pver)           ! qr tendency due to autoconversion of cloud droplets
  real(r8) :: nsagg(pver)           ! ns tendency due to self-aggregation of snow
  real(r8) :: dc0                   ! mean size droplet size distr
  real(r8) :: ds0                   ! mean size snow size distr (area weighted)
  real(r8) :: eci                   ! collection efficiency for riming of snow by droplets
  real(r8) :: dv(pcols,pver)        ! diffusivity of water vapor in air
  real(r8) :: mua(pcols,pver)       ! viscocity of air
  real(r8) :: psacws(pver)          ! mixing rat tendency due to collection of droplets by snow
  real(r8) :: npsacws(pver)         ! number conc tendency due to collection of droplets by snow
  real(r8) :: pracs(pver)           ! mixing rat tendency due to collection of rain	by snow
  real(r8) :: npracs(pver)          ! number conc tendency due to collection of rain by snow
  real(r8) :: mnuccr(pver)          ! mixing rat tendency due to freezing of rain
  real(r8) :: nnuccr(pver)          ! number conc tendency due to freezing of rain
  real(r8) :: pra(pver)             ! mixing rat tendnency due to accretion of droplets by rain
  real(r8) :: npra(pver)            ! nc tendnency due to accretion of droplets by rain
  real(r8) :: nragg(pver)           ! nr tendency due to self-collection of rain
  real(r8) :: prci(pver)            ! mixing rat tendency due to autoconversion of cloud ice to snow
  real(r8) :: nprci(pver)           ! number conc tendency due to autoconversion of cloud ice to snow
  real(r8) :: prai(pver)            ! mixing rat tendency due to accretion of cloud ice by snow
  real(r8) :: nprai(pver)           ! number conc tendency due to accretion of cloud ice by snow
  real(r8) :: prb(pver)             ! rain mixing rat tendency due to Bergeron process
  real(r8) :: nprb(pver)            ! number conc tendency due to Bergeron process
  real(r8) :: fhmrm (pcols,pver)    !mass tendency due to homogeneous freezing of rain

! fall speed
  real(r8) :: arn(pcols,pver)       ! air density corrected rain fallspeed parameter
  real(r8) :: asn(pcols,pver)       ! air density corrected snow fallspeed parameter
  real(r8) :: acn(pcols,pver)       ! air density corrected cloud droplet fallspeed parameter
  real(r8) :: ain(pcols,pver)       ! air density corrected cloud ice fallspeed parameter
  real(r8) :: uns(pver)             ! number-weighted snow fallspeed
  real(r8) :: ums(pver)             ! mass-weighted snow fallspeed
  real(r8) :: unr(pver)             ! number-weighted rain fallspeed
  real(r8) :: umr(pver)             ! mass-weighted rain fallspeed

! conservation check
  real(r8) :: qce                   ! dummy qc for conservation check
  real(r8) :: qie                   ! dummy qi for conservation check
  real(r8) :: nce                   ! dummy nc for conservation check
  real(r8) :: nie                   ! dummy ni for conservation check
  real(r8) :: qre                   ! dummy qr for conservation check
  real(r8) :: nre                   ! dummy nr for conservation check
  real(r8) :: qnie                  ! dummy qni for conservation check
  real(r8) :: nse                   ! dummy ns for conservation check      
  real(r8) :: ratio                 ! parameter for conservation check

! sum of source/sink terms for cloud hydrometeor
  real(r8) :: qctend(pcols,pver)    ! microphysical tendency qc (1/s)
  real(r8) :: qitend(pcols,pver)    ! microphysical tendency qi (1/s)
  real(r8) :: nctend(pcols,pver)    ! microphysical tendency nc (1/(kg*s))
  real(r8) :: nitend(pcols,pver)    ! microphysical tendency ni (1/(kg*s))
  real(r8) :: qnitend(pcols,pver)   ! snow mixing ratio source/sink term
  real(r8) :: nstend(pcols,pver)    ! snow number concentration source/sink term
  real(r8) :: qrtend(pcols,pver)    ! rain mixing ratio source/sink term
  real(r8) :: nrtend(pcols,pver)    ! rain number concentration source/sink term

! terms for Bergeron process
  real(r8) :: bergtsf               !bergeron timescale to remove all liquid
  real(r8) :: plevap                ! cloud liquid water evaporation rate

! variables for droplet activation by modal aerosols
  real(r8) :: wmix, wmin, wmax, wdiab
  real(r8) :: vol, nlsrc
  real(r8), allocatable :: vaerosol(:), hygro(:), naermod(:)
  real(r8), allocatable :: fn(:)      ! number fraction of aerosols activated
  real(r8), allocatable :: fm(:)      ! mass fraction of aerosols activated
  real(r8), allocatable :: fluxn(:)   ! flux of activated aerosol number fraction into cloud (cm/s)
  real(r8), allocatable :: fluxm(:)   ! flux of activated aerosol mass fraction into cloud (cm/s)
  real(r8) :: flux_fullact            ! flux of activated aerosol fraction assuming 100% activation (cm/s)
  real(r8) :: dmc
  real(r8) :: ssmc
  real(r8) :: dgnum_aitken

! bulk aerosol variables
  real(r8), allocatable :: naer2(:,:,:)    ! new aerosol number concentration (/m3)
  real(r8), allocatable :: naer2h(:,:,:)   ! new aerosol number concentration (/m3) 
  real(r8), allocatable :: maerosol(:)     ! aerosol mass conc (kg/m3)
  real(r8) :: so4_num
  real(r8) :: soot_num
  real(r8) :: dst1_num
  real(r8) :: dst2_num
  real(r8) :: dst3_num
  real(r8) :: dst4_num
  real(r8) :: dst_num

! droplet activation
  logical  :: in_cloud              ! true when above cloud base layer (k > jb)
  real(r8) :: smax_f                ! droplet and rain size distr factor used in the
                                    ! in-cloud smax calculation
  real(r8) :: dum2l(pcols,pver)     ! number conc of CCN (1/kg)
  real(r8) :: npccn(pver)           ! droplet activation rate
  real(r8) :: ncmax
  real(r8) :: mtimec                ! factor to account for droplet activation timescale

! ice nucleation
  real(r8) :: dum2i(pcols,pver)     ! number conc of ice nuclei available (1/kg)
  real(r8) :: qs(pcols,pver)        ! liquid-ice weighted sat mixing rat (kg/kg)
  real(r8) :: es(pcols,pver)        ! sat vapor press (pa) over water
  real(r8) :: relhum(pcols,pver)    ! relative humidity
  real(r8) :: esi(pcols,pver)       ! sat vapor press (pa) over ice
  real(r8) :: nnuccd(pver)          ! ice nucleation rate from deposition/cond.-freezing
  real(r8) :: mnuccd(pver)          ! mass tendency from ice nucleation
  real(r8) :: mtime                 ! factor to account for ice nucleation timescale

! output for ice nucleation
  real(r8) :: nimey(pcols,pver)     !number conc of ice nuclei due to meyers deposition (1/m3)
  real(r8) :: nihf(pcols,pver)      !number conc of ice nuclei due to heterogenous freezing (1/m3)
  real(r8) :: nidep(pcols,pver)     !number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
  real(r8) :: niimm(pcols,pver)     !number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)

  real(r8) :: wpice, weff, fhom      ! unused dummies  

! loop array variables
  integer i,k, n, l
  integer ii,kk, m

! loop variables for iteration solution
  integer iter,it,ltrue(pcols)

! used in contact freezing via dust particles
  real(r8)  tcnt, viscosity, mfp
  real(r8)  slip1, slip2, slip3, slip4
  real(r8)  dfaer1, dfaer2, dfaer3, dfaer4
  real(r8)  nacon1,nacon2,nacon3,nacon4

! used in immersion freezing via soot
  real(r8) ttend(pver)
  real(r8) naimm
  real(r8) :: ntaer(pcols,pver)
  real(r8) :: ntaerh(pcols,pver)

! used in homogeneous freezing
  real(r8) :: fholm (pcols,pver)    !mass tendency due to homogeneous freezing
  real(r8) :: fholn (pcols,pver)    !number conc tendency due to homogeneous freezing

! used in secondary ice production
  real(r8) ni_secp

! used in vertical velocity calculation
  real(r8) th(pcols,pver)
  real(r8) qh(pcols,pver)
  real(r8) zkine(pcols,pver) 
  real(r8) zbuo(pcols,pver)    
  real(r8) zfacbuo, cwdrag, cwifrac, retv,  zbuoc
  real(r8) zbc, zbe,  zdkbuo, zdken
  real(r8) arcf(pcols,pver)
  real(r8) p(pcols,pver)
  real(r8) ph(pcols,pver)

! used in vertical integreation
  logical qcimp(pver)             ! true to solve qc with implicit formula
  logical ncimp(pver)             ! true to solve nc with implicit formula
  logical qiimp(pver)             ! true to solve qi with implicit formula
  logical niimp(pver)             ! true to solve ni with implicit formula

! tendency due to adjustment
  real(r8) :: ncadj(pcols,pver)     !droplet num tendency due to adjustment
  real(r8) :: niadj(pcols,pver)     !ice crystal num tendency due to adjustment
  real(r8) :: ncorg, niorg, total

  real(r8) :: rhoh(pcols,pver)    ! air density (kg m-3) at interface 
  real(r8) :: rhom(pcols,pver)    ! air density (kg m-3) at mid-level
  real(r8) :: tu(pcols,pver)      ! temperature in updraft (K)

  integer  kqi(pcols),kqc(pcols)
  logical  lcbase(pcols), libase(pcols)

  real(r8) :: nai_bcphi, nai_dst1, nai_dst2, nai_dst3, nai_dst4

  real(r8) flxrm, mvtrm, flxrn, mvtrn, flxsm, mvtsm, flxsn, mvtsn
  integer  nlr, nls 

  real(r8)  rmean, beta6, beta66, r6, r6c
  real(r8)  temp1, temp2, temp3, temp4   ! variable to store output which is not required by this routine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (aero%scheme == 'modal') then

     allocate(vaerosol(aero%nmodes), hygro(aero%nmodes), naermod(aero%nmodes), &
        fn(aero%nmodes), fm(aero%nmodes), fluxn(aero%nmodes), fluxm(aero%nmodes))

  else if (aero%scheme == 'bulk') then

     allocate( &
        naer2(pcols,pver,aero%nbulk),  &
        naer2h(pcols,pver,aero%nbulk), &
        maerosol(aero%nbulk))

  end if
	
  deltat= get_step_size()      !for FV dynamical core
        
  ! parameters for scheme
  omsm=0.99999_r8
  zfacbuo = 0.5_r8/(1._r8+0.5_r8)
  cwdrag  = 1.875_r8*0.506_r8
  cwifrac = 0.5_r8
  retv    = 0.608_r8
  bergtsf = 1800._r8
 
  ! initialize multi-level fields
  do i=1,il2g
     do k=1,pver
        q(i,k) = qu(i,k)         
        tu(i,k)= su(i,k) - grav/cp*zf(i,k)
        t(i,k) = su(i,k) - grav/cp*zf(i,k)
        p(i,k) = 100._r8*pm(i,k)
        wu(i,k)   = 0._r8
        zkine(i,k)= 0._r8
        arcf(i,k) = 0._r8
        zbuo(i,k) = 0._r8
        nc(i,k)   = 0._r8
        ni(i,k)   = 0._r8
        qc(i,k)   = 0._r8
        qi(i,k)   = 0._r8
        ncde(i,k) = 0._r8
        nide(i,k) = 0._r8
        qcde(i,k) = 0._r8
        qide(i,k) = 0._r8
        rprd(i,k) = 0._r8
        sprd(i,k) = 0._r8
        frz(i,k)  = 0._r8
        qcic(i,k) = 0._r8
        qiic(i,k) = 0._r8
        ncic(i,k) = 0._r8
        niic(i,k) = 0._r8 
        qr(i,k)   = 0._r8
        qni(i,k)  = 0._r8
        nr(i,k)   = 0._r8
        ns(i,k)   = 0._r8
        qric(i,k)  = 0._r8
        qniic(i,k) = 0._r8
        nric(i,k)  = 0._r8
        nsic(i,k)  = 0._r8
        nimey(i,k) = 0._r8
        nihf(i,k)  = 0._r8
        nidep(i,k) = 0._r8
        niimm(i,k) = 0._r8  
        fhmrm(i,k) = 0._r8

        autolm(i,k) = 0._r8
        accrlm(i,k) = 0._r8
        bergnm(i,k) = 0._r8
        fhtimm(i,k) = 0._r8
        fhtctm(i,k) = 0._r8
        fhmlm (i,k) = 0._r8
        fholm (i,k) = 0._r8
        hmpim (i,k) = 0._r8
        accslm(i,k) = 0._r8
        dlfm  (i,k) = 0._r8

        autoln(i,k) = 0._r8
        accrln(i,k) = 0._r8
        bergnn(i,k) = 0._r8
        fhtimn(i,k) = 0._r8
        fhtctn(i,k) = 0._r8
        fhmln (i,k) = 0._r8
        fholn (i,k) = 0._r8
        accsln(i,k) = 0._r8
        activn(i,k) = 0._r8
        dlfn  (i,k) = 0._r8

        autoim(i,k) = 0._r8
        accsim(i,k) = 0._r8
        difm  (i,k) = 0._r8

        nuclin(i,k) = 0._r8
        autoin(i,k) = 0._r8
        accsin(i,k) = 0._r8
        hmpin (i,k) = 0._r8
        difn  (i,k) = 0._r8

        trspcm(i,k) = 0._r8
        trspcn(i,k) = 0._r8
        trspim(i,k) = 0._r8
        trspin(i,k) = 0._r8

        ncadj (i,k) = 0._r8
        niadj (i,k) = 0._r8
     end do
  end do

  ! initialize time-varying parameters
  do k=1,pver
     do i=1,il2g
        if (k .eq.1) then
           rhoh(i,k) = p(i,k)/(t(i,k)*rd)
           rhom(i,k) = p(i,k)/(t(i,k)*rd)
           th (i,k) = te(i,k)
           qh (i,k) = qe(i,k)
           dz (i,k)  = zf(i,k) - zf(i,k+1)
           ph(i,k)   = p(i,k)
        else 
           rhoh(i,k) = 0.5_r8*(p(i,k)+p(i,k-1))/(t(i,k)*rd)
           if (k .eq. pver) then
              rhom(i,k) = p(i,k)/(rd*t(i,k))   
           else
              rhom(i,k) = 2.0_r8*p(i,k)/(rd*(t(i,k)+t(i,k+1)))
           end if
           th (i,k) = 0.5_r8*(te(i,k)+te(i,k-1))
           qh (i,k) = 0.5_r8*(qe(i,k)+qe(i,k-1))
           dz(i,k)  = zf(i,k-1) - zf(i,k)
           ph(i,k)  = 0.5_r8*(p(i,k) + p(i,k-1))
        end if
        dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/ph(i,k)
        mua(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/ &
           (t(i,k)+120._r8)

        rho(i,k) = rhoh(i,k)    

        ! air density adjustment for fallspeed parameters
        ! add air density correction factor to the power of 
        ! 0.54 following Heymsfield and Bansemer 2006

        arn(i,k)=ar*(rhosu/rho(i,k))**0.54_r8
        asn(i,k)=as*(rhosu/rho(i,k))**0.54_r8
        acn(i,k)=ac*(rhosu/rho(i,k))**0.54_r8
        ain(i,k)=ai*(rhosu/rho(i,k))**0.54_r8

     end do
  end do

  if (aero%scheme == 'modal') then

     wmix  = 0._r8
     wmin  = 0._r8
     wmax  = 10._r8
     wdiab = 0._r8

     do k=1,pver
        do i=1,il2g
           dum2l(i,k)=0._r8
           dum2i(i,k)=0._r8
           ntaer(i,k) = 0.0_r8
           ntaerh(i,k) = 0.0_r8
           do m = 1, aero%nmodes
              ntaer(i,k) = ntaer(i,k) + aero%numg_a(i,k,m)*rhom(i,k)
           enddo
        end do
     end do

  else if (aero%scheme == 'bulk') then

     ! initialize aerosol number
     do k=1,pver
        do i=1,il2g
           naer2(i,k,:)=0._r8
           naer2h(i,k,:)=0._r8      
           dum2l(i,k)=0._r8
           dum2i(i,k)=0._r8
        end do
     end do

     do k=1,pver
        do i=1,il2g
           ntaer(i,k) = 0.0_r8
           ntaerh(i,k) = 0.0_r8
           do m = 1, aero%nbulk
              maerosol(m) = aero%mmrg_bulk(i,k,m)*rhom(i,k)

              ! set number nucleated for sulfate based on Lohmann et al. 2000 (JGR) Eq.2
              !    Na=340.*(massSO4)^0.58  where Na=cm-3 and massSO4=ug/m3
              ! convert units to Na [m-3] and SO4 [kgm-3]
              !    Na(m-3)= 1.e6 cm3 m-3 Na(cm-3)=340. *(massSO4[kg/m3]*1.e9ug/kg)^0.58
              !    or Na(m-3)= 1.e6* 340.*(1.e9ug/kg)^0.58 * (massSO4[kg/m3])^0.58

              if (m .eq. aero%idxsul) then
                 naer2(i,k,m)= 5.64259e13_r8 * maerosol(m)**0.58_r8
              else
                 naer2(i,k,m)=maerosol(m)*aero%num_to_mass_aer(m)
              end if
              ntaer(i,k) = ntaer(i,k) + naer2(i,k,m) 
           end do
        end do
     end do

  end if

  do i=1,il2g
     ltrue(i)=0
     do k=1,pver
        if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or.cmel(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
     end do
  end do

  ! skip microphysical calculations if no cloud water
  do i=1,il2g
     if (ltrue(i).eq.0) then
        do k=1,pver
           qctend(i,k)=0._r8
           qitend(i,k)=0._r8
           qnitend(i,k)=0._r8
           qrtend(i,k)=0._r8
           nctend(i,k)=0._r8
           nitend(i,k)=0._r8
           nrtend(i,k)=0._r8
           nstend(i,k)=0._r8
           qniic(i,k)=0._r8
           qric(i,k)=0._r8
           nsic(i,k)=0._r8
           nric(i,k)=0._r8
           qni(i,k)=0._r8
           qr(i,k)=0._r8
           ns(i,k)=0._r8
           nr(i,k)=0._r8
           qc(i,k) = 0._r8
           qi(i,k) = 0._r8
           nc(i,k) = 0._r8
           ni(i,k) = 0._r8
           qcde(i,k) = 0._r8
           qide(i,k) = 0._r8
           ncde(i,k) = 0._r8
           nide(i,k) = 0._r8
           rprd(i,k) = 0._r8
           sprd(i,k) = 0._r8  
           frz(i,k)  = 0._r8
        end do
        goto 300
     end if

     kqc(i) = 1
     kqi(i) = 1
     lcbase(i) = .true.
     libase(i) = .true. 

     ! assign number of steps for iteration
     ! use 2 steps following Song and Zhang, 2011, J. Clim.
     iter = 2

     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !  iteration
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     do it=1,iter

        ! initialize sub-step microphysical tendencies
        do k=1,pver
           qctend(i,k)=0._r8
           qitend(i,k)=0._r8
           qnitend(i,k)=0._r8
           qrtend(i,k)=0._r8
           nctend(i,k)=0._r8
           nitend(i,k)=0._r8
           nrtend(i,k)=0._r8
           nstend(i,k)=0._r8
           rprd(i,k) = 0._r8
           sprd(i,k) = 0._r8
           frz(i,k)  = 0._r8
           qniic(i,k)=0._r8
           qric(i,k)=0._r8
           nsic(i,k)=0._r8
           nric(i,k)=0._r8
           qiic(i,k)=0._r8
           qcic(i,k)=0._r8
           niic(i,k)=0._r8
           ncic(i,k)=0._r8           
           qcimp(k) = .false.
           ncimp(k) = .false.
           qiimp(k) = .false.
           niimp(k) = .false.
           dum2l(i,k)=0._r8
           dum2i(i,k)=0._r8
           autolm(i,k) = 0._r8
           accrlm(i,k) = 0._r8
           bergnm(i,k) = 0._r8
           fhtimm(i,k) = 0._r8
           fhtctm(i,k) = 0._r8
           fhmlm (i,k) = 0._r8
           fholm (i,k) = 0._r8
           hmpim (i,k) = 0._r8
           accslm(i,k) = 0._r8
           dlfm  (i,k) = 0._r8

           autoln(i,k) = 0._r8
           accrln(i,k) = 0._r8
           bergnn(i,k) = 0._r8
           fhtimn(i,k) = 0._r8
           fhtctn(i,k) = 0._r8
           fhmln (i,k) = 0._r8
           fholn (i,k) = 0._r8
           accsln(i,k) = 0._r8
           activn(i,k) = 0._r8
           dlfn  (i,k) = 0._r8
           ncadj (i,k) = 0._r8

           autoim(i,k) = 0._r8
           accsim(i,k) = 0._r8
           difm  (i,k) = 0._r8

           nuclin(i,k) = 0._r8
           autoin(i,k) = 0._r8
           accsin(i,k) = 0._r8
           hmpin (i,k) = 0._r8
           difn  (i,k) = 0._r8
           niadj (i,k) = 0._r8

           trspcm(i,k) = 0._r8
           trspcn(i,k) = 0._r8
           trspim(i,k) = 0._r8
           trspin(i,k) = 0._r8

           fhmrm (i,k) = 0._r8
        end do

        do k = pver,msg+2,-1
  
           if (k > jt(i) .and. k <= jb(i) .and. eps0(i) > 0._r8      &
              .and.mu(i,k).gt.0._r8 .and. mu(i,k-1).gt.0._r8) then

              ! initialize precip fallspeeds to zero
              if (it.eq.1) then
                ums(k)=0._r8 
                uns(k)=0._r8 
                umr(k)=0._r8 
                unr(k)=0._r8
                prf(k)=0._r8
                pnrf(k)=0._r8
                psf(k) =0._r8
                pnsf(k) = 0._r8
              end if
              ttend(k)=0._r8
              nnuccd(k)=0._r8
              npccn(k)=0._r8

              !************************************************************************************
              ! obtain values of cloud water/ice mixing ratios and number concentrations in updraft
              ! for microphysical process calculations
              ! units are kg/kg for mixing ratio, 1/kg for number conc
              !************************************************************************************


              if (it.eq.1) then
                 qcic(i,k) = qc(i,k)
                 qiic(i,k) = qi(i,k)
                 ncic(i,k) = nc(i,k)
                 niic(i,k) = ni(i,k)
                 qniic(i,k)= qni(i,k)
                 qric(i,k) = qr(i,k)
                 nsic(i,k) = ns(i,k)
                 nric(i,k) = nr(i,k)
              else 
                 if (k.le.kqc(i)) then   
                    qcic(i,k) = qc(i,k)
                    ncic(i,k) = nc(i,k)

                    ! consider rain falling from above
                    flxrm = 0._r8
                    mvtrm = 0._r8
                    flxrn = 0._r8
                    mvtrn = 0._r8
                    nlr = 0
                    
                    do kk= k,jt(i)+3,-1
                       if (qr(i,kk-1) .gt. 0._r8) then
                           nlr = nlr + 1
                           flxrm = flxrm + umr(kk-1)*qr(i,kk-1)*arcf(i,kk-1)
                           flxrn = flxrn + unr(kk-1)*nr(i,kk-1)*arcf(i,kk-1)
                           mvtrm = mvtrm + umr(kk-1)*arcf(i,kk-1)
                           mvtrn = mvtrn + unr(kk-1)*arcf(i,kk-1)
                       end if
                    end do
                    if (mvtrm.gt.0) then
                       qric(i,k) = (qr(i,k)*mu(i,k)+flxrm)/(mu(i,k)+mvtrm)
                    else
                       qric(i,k) = qr(i,k)
                    end if
                    if (mvtrn.gt.0) then
                       nric(i,k) = (nr(i,k)*mu(i,k)+flxrn)/(mu(i,k)+mvtrn)
                    else
                       nric(i,k) = nr(i,k)
                    end if

                 end if
                 if (k.eq.kqc(i)) then
                     qcic(i,k) = qc(i,k-1)
                     ncic(i,k) = nc(i,k-1)
                 end if
                 if(k.le.kqi(i)) then
                   qiic(i,k) = qi(i,k)
                   niic(i,k) = ni(i,k)
! consider snow falling from above
                   flxsm = 0._r8
                   mvtsm = 0._r8
                   flxsn = 0._r8
                   mvtsn = 0._r8
                   nls = 0
  
                   do kk= k,jt(i)+3,-1
                     if (qni(i,kk-1) .gt. 0._r8) then
                       nls = nls + 1
                       flxsm = flxsm + ums(kk-1)*qni(i,kk-1)*arcf(i,kk-1)
                       mvtsm = mvtsm + ums(kk-1)*arcf(i,kk-1)
                       flxsn = flxsn + uns(kk-1)*ns(i,kk-1)*arcf(i,kk-1)
                       mvtsn = mvtsn + uns(kk-1)*arcf(i,kk-1)
                     end if
                   end do

                   if (mvtsm.gt.0) then
                      qniic(i,k) = (qni(i,k)*mu(i,k)+flxsm)/(mu(i,k)+mvtsm)
                   else
                      qniic(i,k) = qni(i,k)
                   end if
                   if (mvtsn.gt.0) then
                      nsic(i,k) = (ns(i,k)*mu(i,k)+flxsn)/(mu(i,k)+mvtsn)
                   else
                      nsic(i,k) = ns(i,k)
                   end if
                 end if
                 if(k.eq.kqi(i)) then
                    qiic(i,k) = qi(i,k-1)
                    niic(i,k) = ni(i,k-1)
                 end if
              end if

              !**********************************************************************
              ! boundary condition for cloud liquid water and cloud ice
              !***********************************************************************

              ! boundary condition for provisional cloud water
              if (cmel(i,k-1).gt.qsmall .and. lcbase(i) .and. it.eq.1 ) then
                 kqc(i) = k
                 lcbase(i) = .false.
                 qcic(i,k) = dz(i,k)*cmel(i,k-1)/(mu(i,k-1)+dz(i,k)*du(i,k-1))
                 ncic(i,k) = qcic(i,k)/(4._r8/3._r8*pi*10.e-6_r8**3*rhow)
              end if

              ! boundary condition for provisional cloud ice
              if (qiic(i,k).gt.qsmall .and. libase(i) .and. it.eq.1 ) then
                 kqi(i) = k
                 libase(i) = .false.
              else if ( cmei(i,k-1).gt.qsmall .and.   &
                 cmei(i,k).lt.qsmall .and. k.le.jb(i) .and. libase(i) .and. it.eq.1 ) then
                 kqi(i)=k
                 libase(i) = .false.
                 qiic(i,k) = dz(i,k)*cmei(i,k-1)/(mu(i,k-1)+dz(i,k)*du(i,k-1))
                 niic(i,k) = qiic(i,k)/(4._r8/3._r8*pi*25.e-6_r8**3*rhoi)               
              end if

              !***************************************************************************
              ! get size distribution parameters based on in-cloud cloud water/ice 
              ! these calculations also ensure consistency between number and mixing ratio
              !***************************************************************************
              ! cloud ice
              if (qiic(i,k).ge.qsmall) then

                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 niic(i,k)=min(niic(i,k),qiic(i,k)*1.e20_r8)
                 lami(k) = (gamma(1._r8+di)*ci* &
                    niic(i,k)/qiic(i,k))**(1._r8/di)
                 n0i(k) = niic(i,k)*lami(k)

                 ! check for slope
                 lammax = 1._r8/10.e-6_r8
                 lammin = 1._r8/(2._r8*dcs)

                 ! adjust vars
                 if (lami(k).lt.lammin) then
                    lami(k) = lammin
                    n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*gamma(1._r8+di))
                    niic(i,k) = n0i(k)/lami(k)
                 else if (lami(k).gt.lammax) then
                    lami(k) = lammax
                    n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*gamma(1._r8+di))
                    niic(i,k) = n0i(k)/lami(k)
                 end if
              else
                 lami(k) = 0._r8
                 n0i(k) = 0._r8
              end if

              ! cloud water
              if (qcic(i,k).ge.qsmall) then

                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 ncic(i,k)=min(ncic(i,k),qcic(i,k)*1.e20_r8)

                 ! get pgam from fit to observations of martin et al. 1994

                 pgam(i,k)=0.0005714_r8*(ncic(i,k)/1.e6_r8/rho(i,k))+0.2714_r8
                 pgam(i,k)=1._r8/(pgam(i,k)**2)-1._r8
                 pgam(i,k)=max(pgam(i,k),2._r8)
                 pgam(i,k)=min(pgam(i,k),15._r8)

                 ! calculate lamc
                 lamc(i,k) = (pi/6._r8*rhow*ncic(i,k)*gamma(pgam(i,k)+4._r8)/ &
                    (qcic(i,k)*gamma(pgam(i,k)+1._r8)))**(1._r8/3._r8)

                 ! lammin, 50 micron diameter max mean size
                 lammin = (pgam(i,k)+1._r8)/40.e-6_r8
                 lammax = (pgam(i,k)+1._r8)/1.e-6_r8 

                 if (lamc(i,k).lt.lammin) then
                    lamc(i,k) = lammin
                    ncic(i,k) = 6._r8*lamc(i,k)**3*qcic(i,k)* &
                       gamma(pgam(i,k)+1._r8)/ &
                       (pi*rhow*gamma(pgam(i,k)+4._r8))
                 else if (lamc(i,k).gt.lammax) then
                    lamc(i,k) = lammax
                    ncic(i,k) = 6._r8*lamc(i,k)**3*qcic(i,k)* &
                       gamma(pgam(i,k)+1._r8)/ &
                       (pi*rhow*gamma(pgam(i,k)+4._r8))
                 end if

                 ! parameter to calculate droplet freezing

                 cdist1(k) = ncic(i,k)/gamma(pgam(i,k)+1._r8) 
              else
                 lamc(i,k) = 0._r8
                 cdist1(k) = 0._r8
              end if

              ! boundary condition for cloud liquid water
              if ( kqc(i) .eq. k  ) then
                 qc(i,k) =  0._r8
                 nc(i,k) = 0._r8
              end if

              ! boundary condition for cloud ice
              if (kqi(i).eq.k  ) then
                 qi(i,k) = 0._r8
                 ni(i,k) = 0._r8
              end if
              
              !**************************************************************************
              ! begin micropysical process calculations 
              !**************************************************************************

              !.................................................................
              ! autoconversion of cloud liquid water to rain
              ! formula from Khrouditnov and Kogan (2000)
              ! minimum qc of 1 x 10^-8 prevents floating point error

              if (qcic(i,k).ge.1.e-8_r8) then

                 ! nprc is increase in rain number conc due to autoconversion
                 ! nprc1 is decrease in cloud droplet conc due to autoconversion
                 ! Khrouditnov and Kogan (2000) 
!                 prc(k) = 1350._r8*qcic(i,k)**2.47_r8*    &
!                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)

                ! Liu and Daum(2004)(modified), Wood(2005)
                rmean = 1.e6_r8*((qcic(i,k)/ncic(i,k))/(4._r8/3._r8*pi*rhow))**(1._r8/3._r8)

                if (rmean .ge. 15._r8) then 
 
                  beta6  = (1._r8+3._r8/rmean)**(1._r8/3._r8)
                  beta66 = (1._r8+3._r8/rmean)**2._r8
                  r6     = beta6*rmean
                  r6c    = 7.5_r8/(r6**0.5_r8*(qcic(i,k)*rho(i,k))**(1._r8/6._r8))
                  prc(k) = 1.3e9_r8*beta66*(qcic(i,k)*rho(i,k))**3._r8/            &
                           (ncic(i,k)*rho(i,k))*max(0._r8,r6-r6c)/rho(i,k)

                  nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
                  nprc(k) = nprc1(k)*0.5_r8
                else
                  prc(k)=0._r8
                  nprc(k)=0._r8
                  nprc1(k)=0._r8
                end if
              else
                 prc(k)=0._r8
                 nprc(k)=0._r8
                 nprc1(k)=0._r8
              end if
 
              ! provisional rain mixing ratio and number concentration (qric and nric)
              ! at boundary are estimated via autoconversion

              if (k.eq.kqc(i) .and. it.eq.1) then
                 qric(i,k) = prc(k)*dz(i,k)/0.55_r8
                 nric(i,k) = nprc(k)*dz(i,k)/0.55_r8 
                 qr(i,k) = 0.0_r8
                 nr(i,k) = 0.0_r8
              end if

              !.......................................................................
              ! Autoconversion of cloud ice to snow
              ! similar to Ferrier (1994)

              call ice_autoconversion(t(i,k), qiic(i,k), lami(k), n0i(k), dcs, prci(k), nprci(k), 1)  

              ! provisional snow mixing ratio and number concentration (qniic and nsic) 
              ! at boundary are estimated via autoconversion

              if (k.eq.kqi(i) .and. it.eq.1) then
                 qniic(i,k)= prci(k)*dz(i,k)*0.25_r8
                 nsic(i,k)= nprci(k)*dz(i,k)*0.25_r8
                 qni(i,k)= 0.0_r8
                 ns(i,k)= 0.0_r8
              end if

              ! if precip mix ratio is zero so should number concentration
              if (qniic(i,k).lt.qsmall) then
                 qniic(i,k)=0._r8
                 nsic(i,k)=0._r8
              end if
              if (qric(i,k).lt.qsmall) then
                 qric(i,k)=0._r8
                 nric(i,k)=0._r8
              end if

              ! make sure number concentration is a positive number to avoid 
              ! taking root of negative later
              nric(i,k)=max(nric(i,k),0._r8)
              nsic(i,k)=max(nsic(i,k),0._r8)

              !**********************************************************************
              ! get size distribution parameters for precip
              !**********************************************************************
              ! rain

              if (qric(i,k).ge.qsmall) then
                 lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
                 n0r(k) = nric(i,k)*lamr(k)

                 ! check for slope
                 lammax = 1._r8/150.e-6_r8
                 lammin = 1._r8/3000.e-6_r8

                 ! adjust vars
                 if (lamr(k).lt.lammin) then
                    lamr(k) = lammin
                    n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                    nric(i,k) = n0r(k)/lamr(k)
                 else if (lamr(k).gt.lammax) then
                    lamr(k) = lammax
                    n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                    nric(i,k) = n0r(k)/lamr(k)
                 end if

                 ! provisional rain number and mass weighted mean fallspeed (m/s)
                 ! Eq.18 of Morrison and Gettelman, 2008, J. Climate
                 unr(k) = min(arn(i,k)*gamma(1._r8+br)/lamr(k)**br,10._r8)
                 umr(k) = min(arn(i,k)*gamma(4._r8+br)/(6._r8*lamr(k)**br),10._r8)
              else
                 lamr(k) = 0._r8
                 n0r(k) = 0._r8
                 umr(k) = 0._r8
                 unr(k) = 0._r8
              end if

              !......................................................................
              ! snow
              if (qniic(i,k).ge.qsmall) then
                 lams(k) = (gamma(1._r8+ds)*cs*nsic(i,k)/ &
                    qniic(i,k))**(1._r8/ds)
                 n0s(k) = nsic(i,k)*lams(k)

                 ! check for slope
                 lammax = 1._r8/10.e-6_r8
                 lammin = 1._r8/2000.e-6_r8

                 ! adjust vars
                 if (lams(k).lt.lammin) then
                    lams(k) = lammin
                    n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*gamma(1._r8+ds))
                    nsic(i,k) = n0s(k)/lams(k)
                 else if (lams(k).gt.lammax) then
                    lams(k) = lammax
                    n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*gamma(1._r8+ds))
                    nsic(i,k) = n0s(k)/lams(k)
                 end if

                 ! provisional snow number and mass weighted mean fallspeed (m/s)
                 ums(k) = min(asn(i,k)*gamma(4._r8+bs)/(6._r8*lams(k)**bs),3.6_r8)
                 uns(k) = min(asn(i,k)*gamma(1._r8+bs)/lams(k)**bs,3.6_r8) 
              else
                 lams(k) = 0._r8
                 n0s(k) = 0._r8
                 ums(k) = 0._r8
                 uns(k) = 0._r8
              end if

              !.......................................................................
              ! snow self-aggregation from passarelli, 1978, used by Reisner(1998,Eq.A.35)
              ! this is hard-wired for bs = 0.4 for now
              ! ignore self-collection of cloud ice

              call snow_self_aggregation(t(i,k), rho(i,k), asn(i,k), rhosn, qniic(i,k), nsic(i,k), nsagg(k), 1)

              !.......................................................................
              ! accretion of cloud droplets onto snow/graupel
              ! here use continuous collection equation with
              ! simple gravitational collection kernel
              ! ignore collisions between droplets/cloud ice

              ! ignore collision of snow with droplets above freezing

              call accrete_cloud_water_snow(t(i,k), rho(i,k), asn(i,k), uns(k), mua(i,k),  &
                     qcic(i,k), ncic(i,k), qniic(i,k), pgam(i,k), lamc(i,k), lams(k), n0s(k), &
                     psacws(k), npsacws(k), 1) 

              ! secondary ice production due to accretion of droplets by snow 
              ! (Hallet-Mossop process) (from Cotton et al., 1986)

              call secondary_ice_production(t(i,k), psacws(k), msacwi(k), nsacwi(k), 1)

              !.......................................................................
              ! accretion of rain water by snow
              ! formula from ikawa and saito, 1991, used by reisner et al., 1998

              call accrete_rain_snow(t(i,k), rho(i,k), umr(k), ums(k), unr(k), uns(k), qric(i,k),  &
                        qniic(i,k), lamr(k), n0r(k), lams(k), n0s(k), pracs(k), npracs(k), 1 )

              !.......................................................................
              ! heterogeneous freezing of rain drops
              ! follows from Bigg (1953)

              call heterogeneous_rain_freezing(t(i,k), qric(i,k), nric(i,k), lamr(k), mnuccr(k), nnuccr(k), 1)

              !.......................................................................
              ! accretion of cloud liquid water by rain
              ! formula from Khrouditnov and Kogan (2000)
              ! gravitational collection kernel, droplet fall speed neglected

              call accrete_cloud_water_rain(.true., qric(i,k), qcic(i,k), ncic(i,k), [1._r8], [0._r8], pra(k), npra(k), 1) 

              !.......................................................................
              ! Self-collection of rain drops
              ! from Beheng(1994)

              call self_collection_rain(rho(i,k), qric(i,k), nric(i,k), nragg(k), 1)

              !.......................................................................
              ! Accretion of cloud ice by snow
              ! For this calculation, it is assumed that the Vs >> Vi
              ! and Ds >> Di for continuous collection

              call accrete_cloud_ice_snow(t(i,k), rho(i,k), asn(i,k), qiic(i,k), niic(i,k),    &
                         qniic(i,k), lams(k), n0s(k), prai(k), nprai(k), 1)

              !.......................................................................
              ! fallout term
              prf(k)  = -umr(k)*qric(i,k)/dz(i,k)
              pnrf(k) = -unr(k)*nric(i,k)/dz(i,k)
              psf(k)  = -ums(k)*qniic(i,k)/dz(i,k)
              pnsf(k) = -uns(k)*nsic(i,k)/dz(i,k)

              !........................................................................
              ! calculate vertical velocity in cumulus updraft
              
              if (k.eq.jb(i)) then
                 zkine(i,jb(i)) = 0.5_r8
                 wu   (i,jb(i)) = 1._r8
                 zbuo (i,jb(i)) = (tu(i,jb(i))*(1._r8+retv*qu(i,jb(i)))-    &
                    th(i,jb(i))*(1._r8+retv*qh(i,jb(i))))/   &
                    (th(i,jb(i))*(1._r8+retv*qh(i,jb(i))))
              else
                 if (.true.) then
                    ! ECMWF formula
                    zbc = tu(i,k)*(1._r8+retv*qu(i,k)-qr(i,k)-qni(i,k)-qi(i,k)-qc(i,k))
                    zbe = th(i,k)*(1._r8+retv*qh(i,k))
                    zbuo(i,k) = (zbc-zbe)/zbe
                    zbuoc= (zbuo(i,k)+zbuo(i,k+1))*0.5_r8
                    zdkbuo = dz(i,k+1)*grav*zfacbuo*zbuoc
                    zdken = min(.99_r8,(1._r8+cwdrag)*max(du(i,k),eu(i,k))*dz(i,k+1)/      &
                       max(1.e-10_r8,mu(i,k+1)))
                    zkine(i,k) = (zkine(i,k+1)*(1._r8-zdken)+zdkbuo)/      &
                       (1._r8+zdken)
                 else
                    ! Gregory formula
                    zbc = tu(i,k)*(1._r8+retv*qu(i,k))
                    zbe = th(i,k)*(1._r8+retv*qh(i,k))
                    zbuo(i,k) = (zbc-zbe)/zbe-qr(i,k)-qni(i,k)-qi(i,k)-qc(i,k)
                    zbuoc= (zbuo(i,k)+zbuo(i,k+1))*0.5_r8
                    zdkbuo = dz(i,k+1)*grav*zbuoc*(1.0_r8-0.25_r8)/6._r8
                    zdken = du(i,k)*dz(i,k+1)/max(1.e-10_r8,mu(i,k+1))
                    zkine(i,k) = (zkine(i,k+1)*(1._r8-zdken)+zdkbuo)/      &
                       (1._r8+zdken)
                 end if
                 wu(i,k) = min(15._r8,sqrt(2._r8*max(0.1_r8,zkine(i,k) )))
              end if

              arcf(i,k)= mu(i,k)/wu(i,k)

              !............................................................................
              ! droplet activation
              ! calculate potential for droplet activation if cloud water is present
              ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98

              if (aero%scheme == 'bulk') then
                 naer2h(i,k,:) = 0.5_r8*(naer2(i,k,:) + naer2(i,k-1,:))
              end if

              ntaerh(i,k)   = 0.5_r8*(ntaer(i,k) + ntaer(i,k-1))

              if (qcic(i,k).ge.qsmall ) then

                 if (aero%scheme == 'modal') then

                    nlsrc = 0._r8

                    do m = 1, aero%nmodes
                       vaerosol(m) = 0._r8
                       hygro(m)    = 0._r8
                       do l = 1, aero%nspec(m)
                          vol   = max(0.5_r8*(aero%mmrg_a(i,k,l,m)+aero%mmrg_a(i,k-1,l,m)) , 0._r8)/aero%specdens(l,m)
                          vaerosol(m) = vaerosol(m) + vol
                          hygro(m)    = hygro(m) + vol*aero%spechygro(l,m)
                       end do
                       if (vaerosol(m) > 1.0e-30_r8) then
                          hygro(m)    = hygro(m)/(vaerosol(m))
                          vaerosol(m) = vaerosol(m)*rho(i,k)
                       else
                          hygro(m)    = 0.0_r8
                          vaerosol(m) = 0.0_r8
                       endif
                       naermod(m) = 0.5_r8*(aero%numg_a(i,k,m)+aero%numg_a(i,k-1,m))*rho(i,k)
                       naermod(m) = max(naermod(m), vaerosol(m)*aero%voltonumbhi(m))
                       naermod(m) = min(naermod(m), vaerosol(m)*aero%voltonumblo(m))
                    end do

                    in_cloud = (k < jb(i))
                    smax_f = 0.0_r8
                    if (in_cloud) then
                       if ( qcic(i,k).ge.qsmall )           &
                          smax_f = ncic(i,k)/lamc(i,k) * gamma(2.0_r8 + pgam(i,k))/gamma(1.0_r8 + pgam(i,k))
                       if ( qric(i,k).ge.qsmall)  smax_f = smax_f + nric(i,k)/lamr(k)

                    end if

                    call activate_modal(  &
                       wu(i,k), wmix, wdiab, wmin, wmax,                 &
                       t(i,k), rho(i,k), naermod, aero%nmodes, vaerosol, &
                       hygro, fn, fm,                  &
                       fluxn, fluxm, flux_fullact, in_cloud_in=in_cloud, smax_f=smax_f)

                    do m = 1, aero%nmodes
                       nlsrc = nlsrc + fn(m)*naermod(m)  !  number nucleated
                    end do

                    if (nlsrc .ne. nlsrc) then
                      write(iulog,*) "nlsrc=",nlsrc,"wu(i,k)=",wu(i,k)
                      write(iulog,*) "fn(m)=",fn,"naermod(m)=",naermod,"aero%specdens(l,m)=",aero%specdens
                      write(iulog,*) "vaerosol(m)=",vaerosol,"aero%voltonumbhi(m)=",aero%voltonumbhi
                      write(iulog,*) "aero%voltonumblo(m)=",aero%voltonumblo,"k=",k,"i=",i
                      write(iulog,*) "aero%numg_a(i,k,m)=",aero%numg_a(i,k,:),"rho(i,k)=",rho(i,k)
                      write(iulog,*) "aero%mmrg_a(i,k,l,m)=",aero%mmrg_a(i,k,:,:)                                    
                    end if

                    dum2l(i,k) = nlsrc

                 else if (aero%scheme == 'bulk') then

                    call ndrop_bam_run( &
                       wu(i,k), t(i,k), rho(i,k), naer2h(i,k,:), aero%nbulk, &
                       aero%nbulk, maerosol, dum2)

                    dum2l(i,k) = dum2

                 end if

              else
                 dum2l(i,k) = 0._r8
              end if

              ! get droplet activation rate
              if (qcic(i,k).ge.qsmall .and. t(i,k).gt.238.15_r8 .and. k.gt.jt(i)+2  ) then

                 ! assume aerosols already activated are equal number of existing droplets for simplicity
                 if (k.eq.kqc(i))  then
                    npccn(k) = dum2l(i,k)/deltat
                 else  
                    npccn(k) = (dum2l(i,k)-ncic(i,k))/deltat
                 end if

                 ! make sure number activated > 0
                 npccn(k) = max(0._r8,npccn(k))
                 ncmax = dum2l(i,k)
              else
                 npccn(k)=0._r8
                 ncmax = 0._r8
              end if

              !..............................................................................
              !ice nucleation
              es(i,k)  = svp_water(t(i,k))     ! over water in mixed clouds
              esi(i,k) = svp_ice(t(i,k))     ! over ice
              qs(i,k) = 0.622_r8*es(i,k)/(ph(i,k) - (1.0_r8-0.622_r8)*es(i,k))
              qs(i,k) = min(1.0_r8,qs(i,k))
              if (qs(i,k) < 0.0_r8)  qs(i,k) = 1.0_r8

              relhum(i,k)= 1.0_r8

              if (t(i,k).lt.tmelt ) then

                 ! compute aerosol number for so4, soot, and dust with units #/cm^3
                 so4_num  = 0._r8
                 soot_num = 0._r8
                 dst1_num = 0._r8
                 dst2_num = 0._r8
                 dst3_num = 0._r8
                 dst4_num = 0._r8

                 if (aero%scheme == 'modal') then          

                    !For modal aerosols, assume for the upper troposphere:
                    ! soot = accumulation mode
                    ! sulfate = aiken mode
                    ! dust = coarse mode
                    ! since modal has internal mixtures.
                    soot_num = 0.5_r8*(aero%numg_a(i,k-1,aero%mode_accum_idx)             &
                                      +aero%numg_a(i,k,aero%mode_accum_idx))*rho(i,k)*1.0e-6_r8
                    dmc  = 0.5_r8*(aero%mmrg_a(i,k-1,aero%coarse_dust_idx,aero%mode_coarse_idx)      &
                                  +aero%mmrg_a(i,k,aero%coarse_dust_idx,aero%mode_coarse_idx))
                    ssmc = 0.5_r8*(aero%mmrg_a(i,k-1,aero%coarse_nacl_idx,aero%mode_coarse_idx)      &
                                  +aero%mmrg_a(i,k,aero%coarse_nacl_idx,aero%mode_coarse_idx))  
                    if (dmc > 0._r8) then
                        dst_num = dmc/(ssmc + dmc) *(aero%numg_a(i,k-1,aero%mode_coarse_idx)         & 
                                  + aero%numg_a(i,k,aero%mode_coarse_idx))*0.5_r8*rho(i,k)*1.0e-6_r8
                    else 
                       dst_num = 0.0_r8
                    end if
                    dgnum_aitken = 0.5_r8*(aero%dgnumg(i,k,aero%mode_aitken_idx)+   &
                                           aero%dgnumg(i,k-1,aero%mode_aitken_idx))     
                    if (dgnum_aitken > 0._r8) then
                       ! only allow so4 with D>0.1 um in ice nucleation
                       so4_num  = 0.5_r8*(aero%numg_a(i,k-1,aero%mode_aitken_idx)+          &
                          aero%numg_a(i,k,aero%mode_aitken_idx))*rho(i,k)*1.0e-6_r8         &
                          * (0.5_r8 - 0.5_r8*erf(log(0.1e-6_r8/dgnum_aitken)/  &
                          (2._r8**0.5_r8*log(aero%sigmag_aitken))))
                    else 
                       so4_num = 0.0_r8 
                    end if
                    so4_num = max(0.0_r8, so4_num)

                 else if (aero%scheme == 'bulk') then          

                    if (aero%idxsul > 0) then 
                       so4_num = naer2h(i,k,aero%idxsul)/25._r8 *1.0e-6_r8
                    end if
                    if (aero%idxbcphi > 0) then 
                       soot_num = naer2h(i,k,aero%idxbcphi)/25._r8 *1.0e-6_r8
                    end if
                    if (aero%idxdst1 > 0) then 
                       dst1_num = naer2h(i,k,aero%idxdst1)/25._r8 *1.0e-6_r8
                    end if
                    if (aero%idxdst2 > 0) then 
                       dst2_num = naer2h(i,k,aero%idxdst2)/25._r8 *1.0e-6_r8
                    end if
                    if (aero%idxdst3 > 0) then 
                       dst3_num = naer2h(i,k,aero%idxdst3)/25._r8 *1.0e-6_r8
                    end if
                    if (aero%idxdst4 > 0) then 
                       dst4_num = naer2h(i,k,aero%idxdst4)/25._r8 *1.0e-6_r8
                    end if
                    dst_num = dst1_num + dst2_num + dst3_num + dst4_num

                 end if

                 ! *** Turn off soot nucleation ***
                 soot_num = 0.0_r8

                 ! Liu et al.,J. climate, 2007
                 if ( wu(i,k) .lt. 4.0_r8) then
                    call nucleati( &
                       wu(i,k), t(i,k), ph(i,k), relhum(i,k), 1.0_r8, qcic(i,k), &
                       1.0e-20_r8, 0.0_r8, rho(i,k), so4_num, dst_num, soot_num, 1.0_r8, &
                       dum2i(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k),   &
                       wpice, weff, fhom, temp1, temp2, temp3, temp4, .true.   )
                 end if   
                 nihf(i,k)=nihf(i,k)*rho(i,k)           !  convert from #/kg -> #/m3)
                 niimm(i,k)=niimm(i,k)*rho(i,k)
                 nidep(i,k)=nidep(i,k)*rho(i,k)
                 nimey(i,k)=nimey(i,k)*rho(i,k)

                 if (.false.) then
                    ! cooper curve (factor of 1000 is to convert from L-1 to m-3)
                    !dum2i(i,k)=0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8

                    ! put limit on number of nucleated crystals, set to number at T=-30 C
                    ! cooper (limit to value at -35 C)
                    !dum2i(i,k)=min(dum2i(i,k),208.9e3_r8)/rho(i,k) ! convert from m-3 to kg-1
                 end if

              else
                 dum2i(i,k)=0._r8
              end if

              ! ice nucleation if activated nuclei exist at t<0C 

              if (dum2i(i,k).gt.0._r8.and.t(i,k).lt.tmelt.and. &
                 relhum(i,k)*es(i,k)/esi(i,k).gt. 1.05_r8  .and. k.gt.jt(i)+1) then

                 if (k.eq.kqi(i)) then
                    nnuccd(k)=dum2i(i,k)/deltat
                 else
                    nnuccd(k)=(dum2i(i,k)-niic(i,k))/deltat
                 end if
                 nnuccd(k)=max(nnuccd(k),0._r8)

                 !Calc mass of new particles using new crystal mass...
                 !also this will be multiplied by mtime as nnuccd is...

                 mnuccd(k) = nnuccd(k) * mi0
              else
                 nnuccd(k)=0._r8
                 mnuccd(k) = 0._r8
              end if

              !................................................................................
              ! Bergeron process
              ! If 0C< T <-40C and both ice and liquid exist
         
              if (t(i,k).le.273.15_r8 .and. t(i,k).gt.233.15_r8 .and.  &
                 qiic(i,k).gt.0.5e-6_r8 .and. qcic(i,k).gt. qsmall)  then
                 plevap = qcic(i,k)/bergtsf
                 prb(k) = max(0._r8,plevap) 
                 nprb(k) = prb(k)/(qcic(i,k)/ncic(i,k))
              else
                 prb(k)=0._r8
                 nprb(k)=0._r8
              end if

              !................................................................................
              ! heterogeneous freezing of cloud water (-5C < T < -35C)

              if (qcic(i,k).ge.qsmall .and.ncic(i,k).gt.0._r8 .and. ntaerh(i,k).gt.0._r8 .and.  &
                 t(i,k).le.268.15_r8 .and. t(i,k).gt.238.15_r8 ) then

                 if (aero%scheme == 'bulk')  then
                    ! immersion freezing (Diehl and Wurzler, 2004)
                    ttend(k) = -grav*wu(i,k)/cp/(1.0_r8+gamhat(i,k))

                    nai_bcphi = 0.0_r8
                    nai_dst1  = 0.0_r8
                    nai_dst2  = 0.0_r8
                    nai_dst3  = 0.0_r8
                    nai_dst4  = 0.0_r8

                    if (aero%idxbcphi > 0) nai_bcphi = naer2h(i,k,aero%idxbcphi)
                    if (aero%idxdst1  > 0) nai_dst1  = naer2h(i,k,aero%idxdst1)
                    if (aero%idxdst2  > 0) nai_dst2  = naer2h(i,k,aero%idxdst2)
                    if (aero%idxdst3  > 0) nai_dst3  = naer2h(i,k,aero%idxdst3)
                    if (aero%idxdst4  > 0) nai_dst4  = naer2h(i,k,aero%idxdst4)

                    naimm = (0.00291_r8*nai_bcphi + 32.3_r8*(nai_dst1 + nai_dst2 + &
                       nai_dst3 + nai_dst4))/ntaerh(i,k)             !m-3
                    if (ttend(k) .lt. 0._r8) then
                       nnuccc(k) = -naimm*exp(273.15_r8-t(i,k))*ttend(k)*qcic(i,k)/rhow   ! kg-1s-1                        
                       mnuccc(k) = nnuccc(k)*qcic(i,k)/ncic(i,k)
                    end if
                 else   
                    if (.false.) then
                       ! immersion freezing (Diehl and Wurzler, 2004)
                       ttend(k) = -grav*wu(i,k)/cp/(1.0_r8+gamhat(i,k))
                       naimm = (0.00291_r8*soot_num + 32.3_r8*dst_num )*1.0e6_r8/ntaerh(i,k)             !m-3
                       if (ttend(k) .lt. 0._r8) then
                          nnuccc(k) = -naimm*exp(273.15_r8-t(i,k))*ttend(k)*qcic(i,k)/rhow   ! kg-1s-1
                          mnuccc(k) = nnuccc(k)*qcic(i,k)/ncic(i,k)
                       end if
                    else
                       ! immersion freezing (Bigg, 1953)
                       mnuccc(k) = pi*pi/36._r8*rhow* &
                           cdist1(k)*gamma(7._r8+pgam(i,k))* &
                           bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
                           lamc(i,k)**3/lamc(i,k)**3

                       nnuccc(k) = pi/6._r8*cdist1(k)*gamma(pgam(i,k)+4._r8) &
                           *bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(i,k)**3
                    end if
                 end if

                 ! contact freezing (Young, 1974) with hooks into simulated dust

                 tcnt=(270.16_r8-t(i,k))**1.3_r8
                 viscosity=1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)          
                 mfp=2.0_r8*viscosity/(ph(i,k)  &                  ! Mean free path (m)
                    *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k))))

                 slip1=1.0_r8+(mfp/rn_dst1)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst1/mfp))))! Slip correction factor
                 slip2=1.0_r8+(mfp/rn_dst2)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst2/mfp))))
                 slip3=1.0_r8+(mfp/rn_dst3)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst3/mfp))))
                 slip4=1.0_r8+(mfp/rn_dst4)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst4/mfp))))
                 
                 dfaer1=1.381e-23_r8*t(i,k)*slip1/(6._r8*pi*viscosity*rn_dst1)  ! aerosol diffusivity (m2/s)
                 dfaer2=1.381e-23_r8*t(i,k)*slip2/(6._r8*pi*viscosity*rn_dst2)
                 dfaer3=1.381e-23_r8*t(i,k)*slip3/(6._r8*pi*viscosity*rn_dst3)
                 dfaer4=1.381e-23_r8*t(i,k)*slip4/(6._r8*pi*viscosity*rn_dst4)

                 nacon1=0.0_r8
                 nacon2=0.0_r8
                 nacon3=0.0_r8
                 nacon4=0.0_r8

                 if (aero%scheme == 'modal') then

                    ! For modal aerosols:
                    !  use size '3' for dust coarse mode...
                    !  scale by dust fraction in coarse mode

                    dmc  = 0.5_r8*(aero%mmrg_a(i,k,aero%coarse_dust_idx,aero%mode_coarse_idx)     &
                                  +aero%mmrg_a(i,k-1,aero%coarse_dust_idx,aero%mode_coarse_idx))
                    ssmc = 0.5_r8*(aero%mmrg_a(i,k,aero%coarse_nacl_idx,aero%mode_coarse_idx)     &
                                  +aero%mmrg_a(i,k-1,aero%coarse_nacl_idx,aero%mode_coarse_idx)) 
                    if (dmc > 0.0_r8) then
                        nacon3 = dmc/(ssmc + dmc) * (aero%numg_a(i,k,aero%mode_coarse_idx)     &
                                 + aero%numg_a(i,k-1,aero%mode_coarse_idx))*0.5_r8*rho(i,k)
                    end if

                 else if (aero%scheme == 'bulk') then

                    if (aero%idxdst1.gt.0) then
                       nacon1=naer2h(i,k,aero%idxdst1)*tcnt *0.0_r8
                    endif
                    if (aero%idxdst2.gt.0) then
                       nacon2=naer2h(i,k,aero%idxdst2)*tcnt ! 1/m3
                    endif
                    if (aero%idxdst3.gt.0) then
                       nacon3=naer2h(i,k,aero%idxdst3)*tcnt
                    endif
                    if (aero%idxdst4.gt.0) then
                       nacon4=naer2h(i,k,aero%idxdst4)*tcnt
                    endif
                 end if

                 mnucct(k) = (dfaer1*nacon1+dfaer2*nacon2+dfaer3*nacon3+dfaer4*nacon4)*pi*pi/3._r8*rhow* &
                    cdist1(k)*gamma(pgam(i,k)+5._r8)/lamc(i,k)**4

                 nnucct(k) = (dfaer1*nacon1+dfaer2*nacon2+dfaer3*nacon3+dfaer4*nacon4)*2._r8*pi*  &
                    cdist1(k)*gamma(pgam(i,k)+2._r8)/lamc(i,k)

                 !              if (nnuccc(k).gt.nnuccd(k)) then
                 !                 dum=nnuccd(k)/nnuccc(k)
                 ! scale mixing ratio of droplet freezing with limit
                 !                 mnuccc(k)=mnuccc(k)*dum
                 !                 nnuccc(k)=nnuccd(k)
                 !              end if
 
              else
                 mnuccc(k) = 0._r8
                 nnuccc(k) = 0._r8
                 mnucct(k) = 0._r8
                 nnucct(k) = 0._r8
              end if

              ! freeze cloud liquid water homogeneously at -40 C
              if (t(i,k) < 233.15_r8 .and. qc(i,k) > 0._r8) then

                 ! make sure freezing rain doesn't increase temperature above
                 ! threshold
                 dum = xlf/cp*qc(i,k)
                 if (t(i,k)+dum.gt.233.15_r8) then
                    dum = -(t(i,k)-233.15_r8)*cp/xlf
                    dum = dum/qc(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if
                 fholm(i,k) = mu(i,k)*dum*qc(i,k)
                 fholn(i,k) = mu(i,k)*dum*nc(i,k)
               end if


              !****************************************************************************************
              ! conservation to ensure no negative values of cloud water/precipitation
              ! in case microphysical process rates are large
              ! note: for check on conservation, processes are multiplied by omsm
              ! to prevent problems due to round off error

              ! since activation/nucleation processes are fast, need to take into account
              ! factor mtime = mixing timescale in cloud / model time step
              ! for now mixing timescale is assumed to be 15 min
              !*****************************************************************************************

              mtime=deltat/900._r8
              mtimec=deltat/900._r8

              ! conservation of qc
              ! ice mass production from ice nucleation(deposition/cond.-freezing), mnuccd, 
              ! is considered as a part of cmei.
                            
              qce = mu(i,k)*qc(i,k)-fholm(i,k) +dz(i,k)*cmel(i,k-1)
              dum = arcf(i,k)*(pra(k)+prc(k)+prb(k)+mnuccc(k)+mnucct(k)+msacwi(k)+   &
                               psacws(k))*dz(i,k)
              if( qce.lt.0._r8)  then
                 qcimp(k) = .true.              
                 prc(k) = 0._r8
                 pra(k) = 0._r8
                 prb(k) = 0._r8
                 mnuccc(k) = 0._r8
                 mnucct(k) = 0._r8
                 msacwi(k) = 0._r8
                 psacws(k) = 0._r8
              else  if (dum.gt.qce) then
                 ratio = qce/dum*omsm
                 prc(k) = prc(k)*ratio
                 pra(k) = pra(k)*ratio
                 prb(k) = prb(k)*ratio
                 mnuccc(k) = mnuccc(k)*ratio
                 mnucct(k) = mnucct(k)*ratio  
                 msacwi(k) = msacwi(k)*ratio  
                 psacws(k) = psacws(k)*ratio 
              end if

              ! conservation of nc
              nce = mu(i,k)*nc(i,k)-fholn(i,k) + (arcf(i,k)*npccn(k)*mtimec)*dz(i,k)
              dum = arcf(i,k)*dz(i,k)*(nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
                 npsacws(k)+ nprb(k) )              
              if (nce.lt.0._r8) then
                 ncimp(k) = .true.
                 nprc1(k) = 0._r8
                 npra(k) = 0._r8
                 nnuccc(k) = 0._r8
                 nnucct(k) = 0._r8
                 npsacws(k) = 0._r8
                 nprb(k) = 0._r8                                  
              else if (dum.gt.nce) then
                 ratio = nce/dum*omsm 
                 nprc1(k) = nprc1(k)*ratio
                 npra(k) = npra(k)*ratio
                 nnuccc(k) = nnuccc(k)*ratio
                 nnucct(k) = nnucct(k)*ratio  
                 npsacws(k) = npsacws(k)*ratio
                 nprb(k) = nprb(k)*ratio                           
              end if
 
              ! conservation of qi
              qie = mu(i,k)*qi(i,k)+fholm(i,k) +dz(i,k)*(cmei(i,k-1) +  &
                    ( mnuccc(k)+mnucct(k)+msacwi(k)+prb(k))*arcf(i,k) )
              dum = arcf(i,k)*(prci(k)+ prai(k))*dz(i,k) 
              if (qie.lt.0._r8) then
                 qiimp(k) = .true.
                 prci(k) = 0._r8
                 prai(k) = 0._r8
              else if (dum.gt.qie) then
                 ratio = qie/dum*omsm
                 prci(k) = prci(k)*ratio
                 prai(k) = prai(k)*ratio
              end if

              ! conservation of ni
              nie = mu(i,k)*ni(i,k)+fholn(i,k) +dz(i,k)*(nnuccd(k)*mtime*arcf(i,k)   &
                    +(nnuccc(k)+ nnucct(k))*arcf(i,k) )
              dum = arcf(i,k)*dz(i,k)*(-nsacwi(k)+nprci(k)+ nprai(k))
              if( nie.lt.0._r8) then
                 niimp(k) = .true.
                 nsacwi(k)= 0._r8
                 nprci(k) = 0._r8
                 nprai(k) = 0._r8
              else  if (dum.gt.nie) then
                 ratio = nie/dum*omsm
                 nsacwi(k)= nsacwi(k)*ratio
                 nprci(k) = nprci(k)*ratio
                 nprai(k) = nprai(k)*ratio
              end if

              ! conservation of qr

              qre = mu(i,k)*qr(i,k)+dz(i,k)*(pra(k)+prc(k))*arcf(i,k)
              dum = arcf(i,k)*dz(i,k)*(pracs(k)+ mnuccr(k)-prf(k))
              if (qre.lt.0._r8) then
                 prf(k) = 0._r8
                 pracs(k) = 0._r8
                 mnuccr(k) = 0._r8
              else if (dum.gt.qre) then
                 ratio = qre/dum*omsm
                 prf(k) = prf(k)*ratio
                 pracs(k) = pracs(k)*ratio
                 mnuccr(k) = mnuccr(k)*ratio
              end if

              ! conservation of nr
              nre = mu(i,k)*nr(i,k) + nprc(k)*arcf(i,k)*dz(i,k)
              dum = arcf(i,k)*dz(i,k)*(npracs(k)+nnuccr(k) &
                 -nragg(k)-pnrf(k))
              if(nre.lt.0._r8) then
                 npracs(k)= 0._r8
                 nnuccr(k)= 0._r8
                 nragg(k) = 0._r8
                 pnrf(k) = 0._r8
              else if (dum.gt.nre) then
                 ratio = nre/dum*omsm
                 npracs(k)= npracs(k)*ratio
                 nnuccr(k)= nnuccr(k)*ratio
                 nragg(k) = nragg(k)*ratio
                 pnrf(k) = pnrf(k)*ratio
              end if

              ! conservation of qni

              qnie = mu(i,k)*qni(i,k)+dz(i,k)*( (prai(k)+psacws(k)+prci(k)+     &
                 pracs(k)+mnuccr(k))*arcf(i,k) )
              dum = arcf(i,k)*dz(i,k)*(-psf(k))

              if(qnie.lt.0._r8) then
                 psf(k) = 0._r8
              else if (dum.gt.qnie) then
                 ratio = qnie/dum*omsm
                 psf(k) = psf(k)*ratio
              end if

              ! conservation of ns
              nse = mu(i,k)*ns(i,k)+dz(i,k)*(nprci(k)+nnuccr(k))*arcf(i,k)
              dum = arcf(i,k)*dz(i,k)*(-nsagg(k)-pnsf(k))
              if (nse.lt.0._r8) then
                 nsagg(k) = 0._r8
                 pnsf(k) = 0._r8
              else if (dum.gt.nse) then
                 ratio = nse/dum*omsm
                 nsagg(k) = nsagg(k)*ratio
                 pnsf(k) = pnsf(k)*ratio
              end if

              !*****************************************************************************
              ! get tendencies due to microphysical conversion processes
              !*****************************************************************************

              if (k.le.kqc(i))   then
                 qctend(i,k) = (-pra(k)-prc(k)-prb(k)-mnuccc(k)-mnucct(k)-msacwi(k)- &   
                                psacws(k))

                 qitend(i,k) = (prb(k)+mnuccc(k)+mnucct(k)+msacwi(k)-prci(k)- prai(k))

                 qrtend(i,k) = (pra(k)+prc(k))+(-pracs(k)- mnuccr(k))

                 qnitend(i,k) = (prai(k)+psacws(k)+prci(k))+(pracs(k)+mnuccr(k))

                 ! multiply activation/nucleation by mtime to account for fast timescale

                 nctend(i,k) = npccn(k)*mtimec+(-nnuccc(k)-nnucct(k)-npsacws(k) &    
                               -npra(k)-nprc1(k)-nprb(k))                           

                 nitend(i,k) = nnuccd(k)*mtime+(nnuccc(k)+ nnucct(k)+nsacwi(k)-nprci(k)- &
                               nprai(k))

                 nstend(i,k) = nsagg(k)+nnuccr(k) + nprci(k)

                 nrtend(i,k) = nprc(k)+(-npracs(k)-nnuccr(k) +nragg(k))

                 ! for output
                 ! cloud liquid water-------------

                 autolm(i,k-1) = -prc(k)*arcf(i,k)
                 accrlm(i,k-1) = -pra(k)*arcf(i,k)
                 bergnm(i,k-1) = -prb(k)*arcf(i,k)
                 fhtimm(i,k-1) = -mnuccc(k)*arcf(i,k)
                 fhtctm(i,k-1) = -mnucct(k)*arcf(i,k)
                 hmpim (i,k-1) = -msacwi(k)*arcf(i,k)
                 accslm(i,k-1) = -psacws(k)*arcf(i,k)
                 fhmlm(i,k-1)  = -fholm(i,k)/dz(i,k)

                 autoln(i,k-1) = -nprc1(k)*arcf(i,k)
                 accrln(i,k-1) = -npra(k)*arcf(i,k)
                 bergnn(i,k-1) = -nprb(k)*arcf(i,k)
                 fhtimn(i,k-1) = -nnuccc(k)*arcf(i,k)
                 fhtctn(i,k-1) = -nnucct(k)*arcf(i,k)
                 accsln(i,k-1) = -npsacws(k)*arcf(i,k)
                 activn(i,k-1) = npccn(k)*mtimec*arcf(i,k)
                 fhmln(i,k-1)  = -fholn(i,k)/dz(i,k)
                 
                 !cloud ice------------------------

                 autoim(i,k-1) = -prci(k)*arcf(i,k)
                 accsim(i,k-1) = -prai(k)*arcf(i,k)

                 nuclin(i,k-1) = nnuccd(k)*mtime*arcf(i,k)
                 autoin(i,k-1) = -nprci(k)*arcf(i,k)
                 accsin(i,k-1) = -nprai(k)*arcf(i,k)
                 hmpin (i,k-1) = nsacwi(k)*arcf(i,k)

              else
                 qctend(i,k) = 0._r8
                 qitend(i,k) = 0._r8
                 qrtend(i,k) = 0._r8
                 qnitend(i,k) = 0._r8
                 nctend(i,k) = 0._r8
                 nitend(i,k) = 0._r8
                 nstend(i,k) = 0._r8
                 nrtend(i,k) = 0._r8
              end if

              !********************************************************************************
              !  vertical integration
              !********************************************************************************
              ! snow
              if ( k.le.kqi(i) ) then
                 qni(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                    (mu(i,k)*qni(i,k)+dz(i,k)*(qnitend(i,k)+psf(k))*arcf(i,k) )
                 
                 ns(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                    (mu(i,k)*ns(i,k)+dz(i,k)*(nstend(i,k)+pnsf(k))*arcf(i,k) )

              else
                 qni(i,k-1)=0._r8
                 ns(i,k-1)=0._r8
              end if

              if (qni(i,k-1).le.0._r8) then
                 qni(i,k-1)=0._r8
                 ns(i,k-1)=0._r8
              end if

              ! rain
              if (k.le.kqc(i) ) then
                 qr(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                    (mu(i,k)*qr(i,k)+dz(i,k)*(qrtend(i,k)+prf(k))*arcf(i,k) )

                 nr(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                    (mu(i,k)*nr(i,k)+dz(i,k)*(nrtend(i,k)+pnrf(k))*arcf(i,k) )

              else
                 qr(i,k-1)=0._r8
                 nr(i,k-1)=0._r8
              end if

              if( qr(i,k-1) .le. 0._r8) then
                 qr(i,k-1)=0._r8
                 nr(i,k-1)=0._r8
              end if

              ! freeze rain homogeneously at -40 C

              if (t(i,k-1) < 233.15_r8 .and. qr(i,k-1) > 0._r8) then

                 ! make sure freezing rain doesn't increase temperature above threshold
                 dum = xlf/cp*qr(i,k-1) 
                 if (t(i,k-1)+dum.gt.233.15_r8) then
                    dum = -(t(i,k-1)-233.15_r8)*cp/xlf
                    dum = dum/qr(i,k-1)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if
                 qni(i,k-1)=qni(i,k-1)+dum*qr(i,k-1)
                 ns(i,k-1)=ns(i,k-1)+dum*nr(i,k-1)
                 qr(i,k-1)=(1._r8-dum)*qr(i,k-1)
                 nr(i,k-1)=(1._r8-dum)*nr(i,k-1)
                 fhmrm(i,k-1) = -mu(i,k-1)*dum*qr(i,k-1)/dz(i,k)
              end if


              ! cloud water
              if ( k.le.kqc(i) ) then
                  qc(i,k-1) = (mu(i,k)*qc(i,k)-fholm(i,k)+dz(i,k)*qctend(i,k)*arcf(i,k)      &
                               +dz(i,k)*cmel(i,k-1) )/(mu(i,k-1)+dz(i,k)*du(i,k-1))

                  qcde(i,k) = qc(i,k-1)

                  nc(i,k-1) = (mu(i,k)*nc(i,k) -fholn(i,k) +dz(i,k)*nctend(i,k)*arcf(i,k) )   & 
                              /(mu(i,k-1)+dz(i,k)*du(i,k-1)) 

                  ncde(i,k) = nc(i,k-1)
              else
                 qc(i,k-1)=0._r8
                 nc(i,k-1)=0._r8
              end if

              if (qc(i,k-1).lt.0._r8) write(iulog,*) "negative qc(i,k-1)=",qc(i,k-1)
              dlfm(i,k-1) = -du(i,k-1)*qcde(i,k)
              dlfn(i,k-1) = -du(i,k-1)*ncde(i,k)

              if (qc(i,k-1).le. 0._r8) then
                 qc(i,k-1)=0._r8
                 nc(i,k-1)=0._r8
              end if

              if (nc(i,k-1).lt. 0._r8) then
                 write(iulog,*) "nc(i,k-1)=",nc(i,k-1),"k-1=",k-1,"arcf(i,k)=",arcf(i,k)
                 write(iulog,*) "mu(i,k-1)=",mu(i,k-1),"mu(i,k)=",mu(i,k),"nc(i,k)=",ni(i,k)
                 write(iulog,*) "dz(i,k)=",dz(i,k),"du(i,k-1)=",du(i,k-1),"nctend(i,k)=",nctend(i,k)
                 write(iulog,*) "eu(i,k-1)=",eu(i,k-1)
              end if              
  
              ! cloud ice
              if( k.le.kqi(i)) then  
                  qi(i,k-1) = (mu(i,k)*qi(i,k)+fholm(i,k) +dz(i,k)*qitend(i,k)*arcf(i,k)      &
                               +dz(i,k)*cmei(i,k-1) )/(mu(i,k-1)+dz(i,k)*du(i,k-1))

                  qide(i,k) = qi(i,k-1)
                
                  ni(i,k-1) = (mu(i,k)*ni(i,k)+fholn(i,k)+dz(i,k)*nitend(i,k)*arcf(i,k) )   &
                              /(mu(i,k-1)+dz(i,k)*du(i,k-1))

                  nide(i,k) = ni(i,k-1)
              else
                 qi(i,k-1)=0._r8
                 ni(i,k-1)=0._r8
              end if

              if (qi(i,k-1).lt.0._r8) write(iulog,*) "negative qi(i,k-1)=",qi(i,k-1)
              difm(i,k-1) = -du(i,k-1)*qide(i,k)
              difn(i,k-1) = -du(i,k-1)*nide(i,k)         

              if (qi(i,k-1).le. 0._r8) then
                 qi(i,k-1)=0._r8
                 ni(i,k-1)=0._r8
              end if


              if (ni(i,k-1).lt. 0._r8) then  
                 write(iulog,*) "ni(i,k-1)=",ni(i,k-1),"k-1=",k-1,"arcf(i,k)=",arcf(i,k)
                 write(iulog,*) "mu(i,k-1)=",mu(i,k-1),"mu(i,k)=",mu(i,k),"ni(i,k)=",ni(i,k)
                 write(iulog,*) "dz(i,k)=",dz(i,k),"du(i,k-1)=",du(i,k-1),"nitend(i,k)=",nitend(i,k)
                 write(iulog,*) "eu(i,k-1)=",eu(i,k-1)
              end if            


              frz(i,k-1) = cmei(i,k-1) + arcf(i,k)*(prb(k)+mnuccc(k)+mnucct(k)+msacwi(k)+   &
                     pracs(k)+mnuccr(k)+psacws(k) )-fhmlm(i,k-1)-fhmrm(i,k-1)


              !******************************************************************************
              ! get size distribution parameters based on in-cloud cloud water/ice
              ! these calculations also ensure consistency between number and mixing ratio

              ! following equation(2,3,4) of Morrison and Gettelman, 2008, J. Climate.
              ! Gamma(n)= (n-1)! 
              ! lamc <-> lambda for cloud liquid water
              ! pgam <-> meu    for cloud liquid water
              ! meu=0 for ice,rain and snow         
              !*******************************************************************************

              ! cloud ice
              niorg = ni(i,k-1)
              if (qi(i,k-1).ge.qsmall) then

                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 ni(i,k-1)=min(ni(i,k-1),qi(i,k-1)*1.e20_r8)
                 ! ni should be non-negative
                 ! ni(i,k-1) = max(ni(i,k-1), 0._r8)
                 if (ni(i,k-1).lt. 0._r8) write(iulog,*) "ni(i,k-1)=",ni(i,k-1)

                 lami(k-1) = (gamma(1._r8+di)*ci* &
                    ni(i,k-1)/qi(i,k-1))**(1._r8/di)
                 n0i(k-1) = ni(i,k-1)*lami(k-1)

                 ! check for slope
                 lammax = 1._r8/10.e-6_r8
                 lammin = 1._r8/(2._r8*dcs)

                 ! adjust vars
                 if (lami(k-1).lt.lammin) then
                    lami(k-1) = lammin
                    n0i(k-1) = lami(k-1)**(di+1._r8)*qi(i,k-1)/(ci*gamma(1._r8+di))
                    ni(i,k-1) = n0i(k-1)/lami(k-1)
                 else if (lami(k-1).gt.lammax) then
                    lami(k-1) = lammax
                    n0i(k-1) = lami(k-1)**(di+1._r8)*qi(i,k-1)/(ci*gamma(1._r8+di))
                    ni(i,k-1) = n0i(k-1)/lami(k-1)
                 end if
              else
                 lami(k-1) = 0._r8
                 n0i(k-1) = 0._r8
              end if

              nide(i,k)   = ni(i,k-1)
              difn(i,k-1) = -du(i,k-1)*nide(i,k)

              niadj(i,k-1)= (ni(i,k-1)- niorg)*mu(i,k-1)/dz(i,k)

              if (niadj(i,k-1) .lt. 0._r8) then
                 total = nuclin(i,k-1)-fhtimn(i,k-1)-fhtctn(i,k-1)-fhmln(i,k-1)+ hmpin (i,k-1)
                 if (total .ne. 0._r8) then
                    nuclin(i,k-1) = nuclin(i,k-1) + nuclin(i,k-1)*niadj(i,k-1)/total
                    fhtimn(i,k-1) = fhtimn(i,k-1) + fhtimn(i,k-1)*niadj(i,k-1)/total
                    fhtctn(i,k-1) = fhtctn(i,k-1) + fhtctn(i,k-1)*niadj(i,k-1)/total 
                    fhmln (i,k-1) = fhmln (i,k-1) + fhmln (i,k-1)*niadj(i,k-1)/total
                    hmpin (i,k-1) = hmpin (i,k-1) + hmpin (i,k-1)*niadj(i,k-1)/total
                 else
                    total = 5._r8
                    nuclin(i,k-1) = nuclin(i,k-1) + niadj(i,k-1)/total
                    fhtimn(i,k-1) = fhtimn(i,k-1) + niadj(i,k-1)/total
                    fhtctn(i,k-1) = fhtctn(i,k-1) + niadj(i,k-1)/total
                    fhmln (i,k-1) = fhmln (i,k-1) + niadj(i,k-1)/total
                    hmpin (i,k-1) = hmpin (i,k-1) + niadj(i,k-1)/total
                 end if
              else if (niadj(i,k-1) .gt. 0._r8)   then
                 total = autoin(i,k-1)+accsin(i,k-1)
                 if (total .ne. 0._r8) then
                    autoin(i,k-1) = autoin(i,k-1) + autoin(i,k-1)*niadj(i,k-1)/total
                    accsin(i,k-1) = accsin(i,k-1) + accsin(i,k-1)*niadj(i,k-1)/total 
                 else
                    total = 2._r8
                    autoin(i,k-1) = autoin(i,k-1) + niadj(i,k-1)/total
                    accsin(i,k-1) = accsin(i,k-1) + niadj(i,k-1)/total  
                 end if  
              end if

              !................................................................................
              !cloud water
              ncorg = nc(i,k-1)
              if (qc(i,k-1).ge.qsmall) then

                 ! add upper limit to in-cloud number concentration to prevent numerical error
                 nc(i,k-1)=min(nc(i,k-1),qc(i,k-1)*1.e20_r8)
                 ! and make sure it's non-negative
                 ! nc(i,k-1) = max(nc(i,k-1), 0._r8)
                 if (nc(i,k-1).lt. 0._r8) write(iulog,*) "nc(i,k-1)=",nc(i,k-1) 

                 ! get pgam from fit to observations of martin et al. 1994

                 pgam(i,k-1)=0.0005714_r8*(nc(i,k-1)/1.e6_r8/rho(i,k-1))+0.2714_r8
                 pgam(i,k-1)=1._r8/(pgam(i,k-1)**2)-1._r8
                 pgam(i,k-1)=max(pgam(i,k-1),2._r8)
                 pgam(i,k-1)=min(pgam(i,k-1),15._r8)
                 ! calculate lamc

                 lamc(i,k-1) = (pi/6._r8*rhow*nc(i,k-1)*gamma(pgam(i,k-1)+4._r8)/ &
                    (qc(i,k-1)*gamma(pgam(i,k-1)+1._r8)))**(1._r8/3._r8)

                 ! lammin, 50 micron diameter max mean size
                 lammin = (pgam(i,k-1)+1._r8)/40.e-6_r8
                 lammax = (pgam(i,k-1)+1._r8)/1.e-6_r8

                 if (lamc(i,k-1).lt.lammin) then
                    lamc(i,k-1) = lammin
                    nc(i,k-1) = 6._r8*lamc(i,k-1)**3*qc(i,k-1)* &
                       gamma(pgam(i,k-1)+1._r8)/ &
                       (pi*rhow*gamma(pgam(i,k-1)+4._r8))
                 else if (lamc(i,k-1).gt.lammax) then
                    lamc(i,k-1) = lammax
                    nc(i,k-1) = 6._r8*lamc(i,k-1)**3*qc(i,k-1)* &
                       gamma(pgam(i,k-1)+1._r8)/ &
                       (pi*rhow*gamma(pgam(i,k-1)+4._r8))
                 end if

                 ! parameter to calculate droplet freezing

                 cdist1(k-1) = nc(i,k-1)/gamma(pgam(i,k-1)+1._r8)
              else
                 lamc(i,k-1) = 0._r8
                 cdist1(k-1) = 0._r8
              end if

              ncde(i,k)   = nc(i,k-1)
              dlfn(i,k-1) = -du(i,k-1)*ncde(i,k)
              
              ncadj(i,k-1) = (nc(i,k-1)- ncorg)*mu(i,k-1)/dz(i,k)
              if (ncadj(i,k-1) .lt. 0._r8) then 
                 activn(i,k-1) = activn(i,k-1) + ncadj(i,k-1)
              else if (ncadj(i,k-1) .gt. 0._r8) then
                total = autoln(i,k-1)+accrln(i,k-1)+bergnn(i,k-1)+accsln(i,k-1)
                if (total .ne. 0._r8) then
                    autoln(i,k-1) = autoln(i,k-1) + autoln(i,k-1)*ncadj(i,k-1)/total
                    accrln(i,k-1) = accrln(i,k-1) + accrln(i,k-1)*ncadj(i,k-1)/total
                    bergnn(i,k-1) = bergnn(i,k-1) + bergnn(i,k-1)*ncadj(i,k-1)/total
                    accsln(i,k-1) = accsln(i,k-1) + accsln(i,k-1)*ncadj(i,k-1)/total
                else
                    total = 4._r8
                    autoln(i,k-1) = autoln(i,k-1) + ncadj(i,k-1)/total
                    accrln(i,k-1) = accrln(i,k-1) + ncadj(i,k-1)/total
                    bergnn(i,k-1) = bergnn(i,k-1) + ncadj(i,k-1)/total
                    accsln(i,k-1) = accsln(i,k-1) + ncadj(i,k-1)/total
                end if
              end if

              trspcm(i,k-1) = (mu(i,k)*qc(i,k) - mu(i,k-1)*qc(i,k-1))/dz(i,k)
              trspcn(i,k-1) = (mu(i,k)*nc(i,k) - mu(i,k-1)*nc(i,k-1))/dz(i,k)
              trspim(i,k-1) = (mu(i,k)*qi(i,k) - mu(i,k-1)*qi(i,k-1))/dz(i,k)
              trspin(i,k-1) = (mu(i,k)*ni(i,k) - mu(i,k-1)*ni(i,k-1))/dz(i,k)

              if (k-1 .eq. jt(i)+1)  then
                 trspcm(i,k-2) =  mu(i,k-1)*qc(i,k-1)/dz(i,k-1)
                 trspcn(i,k-2) =  mu(i,k-1)*nc(i,k-1)/dz(i,k-1)
                 trspim(i,k-2) =  mu(i,k-1)*qi(i,k-1)/dz(i,k-1)
                 trspin(i,k-2) =  mu(i,k-1)*ni(i,k-1)/dz(i,k-1)
                 qcde(i,k-1)   = qc(i,k-1)
                 ncde(i,k-1)   = nc(i,k-1)
                 qide(i,k-1)   = qi(i,k-1)
                 nide(i,k-1)   = ni(i,k-1)
                 dlfm  (i,k-2) = -du(i,k-2)*qcde(i,k-1)
                 dlfn  (i,k-2) = -du(i,k-2)*ncde(i,k-1)
                 difm  (i,k-2) = -du(i,k-2)*qide(i,k-1)
                 difn  (i,k-2) = -du(i,k-2)*nide(i,k-1)
              end if


              !.......................................................................
              ! get size distribution parameters for precip
              !......................................................................
              ! rain
              if (qr(i,k-1).ge.qsmall) then

                 lamr(k-1) = (pi*rhow*nr(i,k-1)/qr(i,k-1))**(1._r8/3._r8)
                 n0r(k-1) = nr(i,k-1)*lamr(k-1)

                 ! check for slope
                 lammax = 1._r8/150.e-6_r8
                 lammin = 1._r8/3000.e-6_r8
                 ! adjust vars
                 if (lamr(k-1).lt.lammin) then
                    lamr(k-1) = lammin
                    n0r(k-1) = lamr(k-1)**4*qr(i,k-1)/(pi*rhow)
                    nr(i,k-1) = n0r(k-1)/lamr(k-1)
                 else if (lamr(k-1).gt.lammax) then
                    lamr(k-1) = lammax
                    n0r(k-1) = lamr(k-1)**4*qr(i,k-1)/(pi*rhow)
                    nr(i,k-1) = n0r(k-1)/lamr(k-1)
                 end if

                 unr(k-1) = min(arn(i,k-1)*gamma(1._r8+br)/lamr(k-1)**br,10._r8)
                 umr(k-1) = min(arn(i,k-1)*gamma(4._r8+br)/(6._r8*lamr(k-1)**br),10._r8)                   
              else
                 lamr(k-1) = 0._r8
                 n0r(k-1) = 0._r8
                 umr(k-1) = 0._r8
                 unr(k-1) = 0._r8
              end if

              !......................................................................
              ! snow
              if (qni(i,k-1).ge.qsmall) then
                 lams(k-1) = (gamma(1._r8+ds)*cs*ns(i,k-1)/ &
                    qni(i,k-1))**(1._r8/ds)
                 n0s(k-1) = ns(i,k-1)*lams(k-1)

                 ! check for slope
                 lammax = 1._r8/10.e-6_r8
                 lammin = 1._r8/2000.e-6_r8

                 ! adjust vars
                 if (lams(k-1).lt.lammin) then
                    lams(k-1) = lammin
                    n0s(k-1) = lams(k-1)**(ds+1._r8)*qni(i,k-1)/(cs*gamma(1._r8+ds))
                    ns(i,k-1) = n0s(k-1)/lams(k-1)
                 else if (lams(k-1).gt.lammax) then
                    lams(k-1) = lammax
                    n0s(k-1) = lams(k-1)**(ds+1._r8)*qni(i,k-1)/(cs*gamma(1._r8+ds))
                    ns(i,k-1) = n0s(k-1)/lams(k-1)
                 end if
                 ums(k-1) = min(asn(i,k-1)*gamma(4._r8+bs)/(6._r8*lams(k-1)**bs),3.6_r8)
                 uns(k-1) = min(asn(i,k-1)*gamma(1._r8+bs)/lams(k-1)**bs,3.6_r8)
              else
                 lams(k-1) = 0._r8
                 n0s(k-1) = 0._r8
                 ums(k-1) = 0._r8
                 uns(k-1) = 0._r8
              end if

              rprd(i,k-1)= (qnitend(i,k) + qrtend(i,k))*arcf(i,k)
              sprd(i,k-1)=  qnitend(i,k) *arcf(i,k) -fhmrm(i,k-1)       

           end if  ! k<jlcl

           ! if rain/snow mix ratio is zero so should number concentration

           if (qni(i,k-1).lt.qsmall) then
              qni(i,k-1)=0._r8
              ns(i,k-1)=0._r8
           end if

           if (qr(i,k-1).lt.qsmall) then
              qr(i,k-1)=0._r8
              nr(i,k-1)=0._r8
           end if

           if (qi(i,k-1).lt.qsmall) then
              qi(i,k-1)=0._r8
              ni(i,k-1)=0._r8
           end if

           if (qc(i,k-1).lt.qsmall) then
              qc(i,k-1)=0._r8
              nc(i,k-1)=0._r8
           end if

           ! make sure number concentration is a positive number to avoid 
           ! taking root of negative

           nr(i,k-1)=max(nr(i,k-1),0._r8)
           ns(i,k-1)=max(ns(i,k-1),0._r8)
           ni(i,k-1)=max(ni(i,k-1),0._r8)
           nc(i,k-1)=max(nc(i,k-1),0._r8)

           !!......................................................................
        end do ! k loop

     end do ! it loop, iteration

300  continue  ! continue if no cloud water

  end do ! i loop

  !........................................................................

  if (aero%scheme == 'modal') then
     deallocate( &
        vaerosol, &
        hygro,    &
        naermod,  &
        fn,       &
        fm,       &
        fluxn,    &
        fluxm     )
  else if (aero%scheme == 'bulk') then
     deallocate( &
        naer2,    &
        naer2h,   &
        maerosol  )
  end if

end subroutine zm_mphy

!##############################################################################

end module zm_microphysics
