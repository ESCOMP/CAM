module crmx_vars

use crmx_grid
#ifdef CRM 
#ifdef MODAL_AERO 
use modal_aero_data,   only: ntot_amode
#endif
#endif

implicit none
!--------------------------------------------------------------------
! prognostic variables:

real u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
real v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
real w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
real t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! liquid/ice water static energy 

!--------------------------------------------------------------------
! diagnostic variables:

real p      (0:nx, (1-YES3D):ny, nzm)     ! perturbation pressure (from Poison eq)
real tabs   (nx, ny, nzm)                 ! temperature
real qv      (nx, ny, nzm)                ! water vapor
real qcl     (nx, ny, nzm)                ! liquid water  (condensate)
real qpl     (nx, ny, nzm)                ! liquid water  (precipitation)
real qci     (nx, ny, nzm)                ! ice water  (condensate)
real qpi     (nx, ny, nzm)                ! ice water  (precipitation)

real tke2(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
real tk2     (0:nxp1, (1-YES3D):nyp1, nzm) ! SGS eddyviscosity
        
!--------------------------------------------------------------------
! time-tendencies for prognostic variables

real dudt   (nxp1, ny, nzm, 3)
real dvdt   (nx, nyp1, nzm, 3)
real dwdt   (nx, ny, nz,  3)

!----------------------------------------------------------------
! Temporary storage array:

	real misc(nx, ny, nz)
!------------------------------------------------------------------
! fluxes at the top and bottom of the domain:

real fluxbu (nx, ny), fluxbv (nx, ny), fluxbt (nx, ny)
real fluxbq (nx, ny), fluxtu (nx, ny), fluxtv (nx, ny)
real fluxtt (nx, ny), fluxtq (nx, ny), fzero  (nx, ny)
real precsfc(nx,ny) ! surface precip. rate
real precssfc(nx,ny) ! surface ice precip. rate
                
!-----------------------------------------------------------------
! profiles 

real   t0(nzm), q0(nzm), qv0(nzm), tabs0(nzm), tl0(nzm), &
       tv0(nzm), u0(nzm), v0(nzm), &
       tg0(nzm), qg0(nzm), ug0(nzm), vg0(nzm), p0(nzm), &
       tke0(nzm), t01(nzm), q01(nzm), qp0(nzm), qn0(nzm)
!----------------------------------------------------------------
! "observed" (read from snd file) surface characteristics 

real  sstobs, lhobs, shobs
!----------------------------------------------------------------
!  Domain top stuff:

real   gamt0    ! gradient of t() at the top,K/m
real   gamq0    ! gradient of q() at the top,g/g/m

!-----------------------------------------------------------------
! reference vertical profiles:
 
real   prespot(nzm)  ! (1000./pres)**R/cp
real   rho(nzm)	  ! air density at pressure levels,kg/m3 
real   rhow(nz)   ! air density at vertical velocity levels,kg/m3
real   bet(nzm)	  ! = ggr/tv0
real   gamaz(nzm) ! ggr/cp*z
real   wsub(nz)   ! Large-scale subsidence velocity,m/s
real   qtend(nzm) ! Large-scale tendency for total water
real   ttend(nzm) ! Large-scale tendency for temp.
real   utend(nzm) ! Large-scale tendency for u
real   vtend(nzm) ! Large-scale tendency for v


!---------------------------------------------------------------------
! Large-scale and surface forcing:

integer nlsf	! number of large-scale forcing profiles
integer nrfc	! number of radiative forcing profiles
integer nsfc	! number of surface forcing profiles
integer nsnd	! number of observed soundings
integer nzlsf	! number of large-scale forcing profiles
integer nzrfc	! number of radiative forcing profiles
integer nzsnd	! number of observed soundings

real, allocatable :: dqls(:,:) ! Large-scale tendency for total water
real, allocatable :: dtls(:,:) ! Large-scale tendency for temp.
real, allocatable :: ugls(:,:) ! Large-scale wind in X-direction
real, allocatable :: vgls(:,:) ! Large-scale wind in Y-direction
real, allocatable :: wgls(:,:) ! Large-scale subsidence velocity,m/s
real, allocatable :: pres0ls(:)! Surface pressure, mb
real, allocatable :: zls(:,:)  ! Height
real, allocatable :: pls(:,:)  ! Pressure
real, allocatable :: dayls(:)  ! Large-scale forcing arrays time (days) 
real, allocatable :: dtrfc(:,:)! Radiative tendency for pot. temp.
real, allocatable :: dayrfc(:) ! Radiative forcing arrays time (days) 
real, allocatable :: prfc(:,:) ! Pressure/Height
real, allocatable :: sstsfc(:) ! SSTs
real, allocatable :: shsfc(:)   ! Sensible heat flux,W/m2
real, allocatable :: lhsfc(:)  ! Latent heat flux,W/m2
real, allocatable :: tausfc(:) ! Surface drag,m2/s2
real, allocatable :: daysfc(:) ! Surface forcing arrays time (days) 
real, allocatable :: usnd(:,:) ! Observed zonal wind
real, allocatable :: vsnd(:,:) ! Observed meriod wind
real, allocatable :: tsnd(:,:) ! Observed Abs. temperature
real, allocatable :: qsnd(:,:) ! Observed Moisture
real, allocatable :: zsnd(:,:) ! Height
real, allocatable :: psnd(:,:) ! Pressure
real, allocatable :: daysnd(:) ! number of sounding samples
 
!---------------------------------------------------------------------
!  Horizontally varying stuff (as a function of xy)
!
real sstxy(0:nx,(1-YES3D):ny)	!  surface temperature xy-distribution
real fcory(0:ny)      !  Coriolis parameter xy-distribution
real fcorzy(ny)      !  z-Coriolis parameter xy-distribution
real latitude(nx,ny)	     ! latitude (degrees)
real longitude(nx,ny)	     ! longitude(degrees)
real prec_xy(nx,ny) ! mean precip. rate for outout
real shf_xy(nx,ny) ! mean precip. rate for outout
real lhf_xy(nx,ny) ! mean precip. rate for outout
real lwns_xy(nx,ny) ! mean net lw at SFC
real swns_xy(nx,ny) ! mean net sw at SFC
real lwnsc_xy(nx,ny) ! clear-sky mean net lw at SFC
real swnsc_xy(nx,ny) ! clear-sky mean net sw at SFC
real lwnt_xy(nx,ny) ! mean net lw at TOA
real swnt_xy(nx,ny) ! mean net sw at TOA
real lwntc_xy(nx,ny) ! clear-sky mean net lw at TOA
real swntc_xy(nx,ny) ! clear-sky mean net sw at TOA
real solin_xy(nx,ny) ! solar TOA insolation
real pw_xy(nx,ny)   ! precipitable water
real cw_xy(nx,ny)   ! cloud water path
real iw_xy(nx,ny)   ! ice water path
real cld_xy(nx,ny)   ! cloud frequency
real u200_xy(nx,ny) ! u-wind at 200 mb
real usfc_xy(nx,ny) ! u-wind at at the surface
real v200_xy(nx,ny) ! v-wind at 200 mb
real vsfc_xy(nx,ny) ! v-wind at the surface
real w500_xy(nx,ny) ! w at 500 mb
real qocean_xy(nx,ny) ! ocean cooling in W/m2

!----------------------------------------------------------------------
!	Vertical profiles of quantities sampled for statitistics purposes:

real &
    twle(nz), twsb(nz), precflux(nz), &
    uwle(nz), uwsb(nz), vwle(nz), vwsb(nz), &
    radlwup(nz), radlwdn(nz), radswup(nz), radswdn(nz), &
    radqrlw(nz), radqrsw(nz), w_max, u_max, s_acld, s_acldcold, s_ar, s_arthr, s_sst, &
    s_acldl, s_acldm, s_acldh,  ncmn, nrmn, z_inv, z_cb, z_ct, z_cbmn, z_ctmn, &
    z2_inv, z2_cb, z2_ct, cwpmean, cwp2, precmean, prec2, precmax, nrainy, ncloudy, &
    s_acldisccp, s_acldlisccp, s_acldmisccp, s_acldhisccp, s_ptopisccp, &
    s_acldmodis, s_acldlmodis, s_acldmmodis, s_acldhmodis, s_ptopmodis, &
    s_acldmisr, s_ztopmisr, s_relmodis, s_reimodis, s_lwpmodis, s_iwpmodis, &
    s_tbisccp, s_tbclrisccp, s_acldliqmodis, s_acldicemodis, &
    s_cldtauisccp,s_cldtaumodis,s_cldtaulmodis,s_cldtauimodis,s_cldalbisccp, &
    s_flns,s_flnt,s_flntoa,s_flnsc,s_flntoac,s_flds,s_fsns, &
    s_fsnt,s_fsntoa,s_fsnsc,s_fsntoac,s_fsds,s_solin, & 
    tkeleadv(nz), tkelepress(nz), tkelediss(nz), tkelediff(nz),tkelebuoy(nz), &
    t2leadv(nz),t2legrad(nz),t2lediff(nz),t2leprec(nz),t2lediss(nz), &
    q2leadv(nz),q2legrad(nz),q2lediff(nz),q2leprec(nz),q2lediss(nz), &
    twleadv(nz),twlediff(nz),twlepres(nz),twlebuoy(nz),twleprec(nz), &
    qwleadv(nz),qwlediff(nz),qwlepres(nz),qwlebuoy(nz),qwleprec(nz), &
    momleadv(nz,3),momlepress(nz,3),momlebuoy(nz,3), &
    momlediff(nz,3),tadv(nz),tdiff(nz),tlat(nz), tlatqi(nz),qifall(nz),qpfall(nz)
real tdiff_xy(nz), tdiff_z(nz), ttest0(nzm), ttest1(nz), ttest2(nz, 10)  !+++mhwang test


! register functions:


real, external :: esatw_crm,esati_crm,dtesatw_crm,dtesati_crm
real, external :: qsatw_crm,qsati_crm,dtqsatw_crm,dtqsati_crm
integer, external :: lenstr, bytes_in_rec

! energy conservation diagnostics:
 
  real(kind=selected_real_kind(12)) total_water_before, total_water_after
  real(kind=selected_real_kind(12)) total_water_evap, total_water_prec, total_water_ls
!#ifdef CLUBB_CRM
  real(kind=selected_real_kind(12)) total_water_clubb
  real(kind=selected_real_kind(12)) total_energy_before, total_energy_after
  real(kind=selected_real_kind(12)) total_energy_evap, total_energy_prec, total_energy_ls
  real(kind=selected_real_kind(12)) total_energy_clubb, total_energy_rad
!#endif
  real(kind=selected_real_kind(12))  qtotmicro(5)  ! total water for water conservation test in microphysics +++mhwang

!===========================================================================
! UW ADDITIONS

! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
integer :: ncondavg, icondavg_cld, icondavg_cor, icondavg_cordn, &
     icondavg_satdn, icondavg_satup, icondavg_env
real, allocatable :: condavg_factor(:,:) ! replaces cloud_factor, core_factor
real, allocatable :: condavg_mask(:,:,:,:) ! indicator array for various conditional averages
character(LEN=8), allocatable :: condavgname(:) ! array of short names
character(LEN=25), allocatable :: condavglongname(:) ! array of long names

real   qlsvadv(nzm) ! Large-scale vertical advection tendency for total water
real   tlsvadv(nzm) ! Large-scale vertical advection tendency for temperature
real   ulsvadv(nzm) ! Large-scale vertical advection tendency for zonal velocity
real   vlsvadv(nzm) ! Large-scale vertical advection tendency for meridional velocity

real   qnudge(nzm) ! Nudging of horiz.-averaged total water profile
real   tnudge(nzm) ! Nudging of horiz.-averaged temperature profile
real   unudge(nzm) ! Nudging of horiz.-averaged zonal velocity
real   vnudge(nzm) ! Nudging of horiz.-averaged meridional velocity

real   qstor(nzm) ! Storage of horiz.-averaged total water profile
real   tstor(nzm) ! Storage of horiz.-averaged temperature profile
real   ustor(nzm) ! Storage of horiz.-averaged zonal velocity
real   vstor(nzm) ! Storage of horiz.-averaged meridional velocity
real   qtostor(nzm) ! Storage of horiz.-averaged total water profile (vapor + liquid)

real   utendcor(nzm) ! coriolis acceleration of zonal velocity
real   vtendcor(nzm) ! coriolis acceleration of meridional velocity

real   CF3D(1:nx, 1:ny, 1:nzm)  ! Cloud fraction 
                                ! =1.0 when there is no fractional cloudiness scheme
                                ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

! 850 mbar horizontal winds
real u850_xy(nx,ny) ! zonal velocity at 850 mb
real v850_xy(nx,ny) ! meridional velocity at 850 mb

! Surface pressure
real psfc_xy(nx,ny) ! pressure (in millibar) at lowest grid point

! Saturated water vapor path, useful for computing column relative humidity
real swvp_xy(nx,ny)  ! saturated water vapor path (wrt water)

! Cloud and echo top heights, and cloud top temperature (instantaneous)
real cloudtopheight(nx,ny), echotopheight(nx,ny), cloudtoptemp(nx,ny)

! END UW ADDITIONS
!===========================================================================
! Initial bubble parameters. Activated when perturb_type = 2
  real bubble_x0  
  real bubble_y0 
  real bubble_z0 
  real bubble_radius_hor 
  real bubble_radius_ver 
  real bubble_dtemp 
  real bubble_dq 
  real, allocatable ::  naer(:,:)     ! Aerosol number concentration [/m3]
  real, allocatable ::  vaer(:,:)     ! aerosol volume concentration [m3/m3]
  real, allocatable ::  hgaer(:,:)    ! hygroscopicity of aerosol mode

end module crmx_vars
