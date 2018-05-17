module crmx_params

use crmx_grid, only: nzm
#ifdef CLUBB_CRM
! Use the CLUBB values for these constants for consistency
use crmx_constants_clubb, only: Cp_clubb => Cp, grav_clubb => grav, Lv_clubb => Lv, Lf_clubb => Lf, &
  Ls_clubb => Ls, Rv_clubb => Rv, Rd_clubb => Rd, pi_clubb => pi
#else

#ifdef CRM
use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
                           shr_const_latice, shr_const_latsub, shr_const_rgas, &
                           shr_const_mwwv, shr_const_stebol, shr_const_tkfrz, &
                           shr_const_mwdair, shr_const_g, shr_const_karman, &
                           shr_const_rhofw
#endif /*CRM*/

#endif

implicit none

!   Constants:

#ifdef CLUBB_CRM
! Define Cp, ggr, etc. in module constants_clubb
real, parameter :: cp =  Cp_clubb 
real, parameter :: ggr =  grav_clubb 
real, parameter :: lcond =  Lv_clubb 
real, parameter :: lfus =  Lf_clubb 
real, parameter :: lsub =  Ls_clubb 
real, parameter :: rv =  Rv_clubb 
real, parameter :: rgas=  Rd_clubb 
#else
#ifndef CRM
real, parameter :: cp = 1004.             ! Specific heat of air, J/kg/K
real, parameter :: ggr = 9.81             ! Gravity acceleration, m/s2
real, parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
real, parameter :: lfus = 0.3336e+06      ! Latent heat of fusion, J/kg
real, parameter :: lsub = 2.8440e+06      ! Latent heat of sublimation, J/kg
real, parameter :: rv = 461.              ! Gas constant for water vapor, J/kg/K
real, parameter :: rgas = 287.            ! Gas constant for dry air, J/kg/K
#else
real, parameter :: cp = shr_const_cpdair
real, parameter :: ggr = shr_const_g
real, parameter :: lcond = shr_const_latvap 
real, parameter :: lfus = shr_const_latice 
real, parameter :: lsub = lcond + lfus
real, parameter :: rv = shr_const_rgas/shr_const_mwwv 
real, parameter :: rgas = shr_const_rdair
#endif
#endif
real, parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
real, parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
real, parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

real, parameter :: fac_cond = lcond/cp 
real, parameter :: fac_fus = lfus/cp
real, parameter :: fac_sub = lsub/cp

#ifdef CLUBB_CRM
real, parameter ::  pi =  pi_clubb 
#else
real, parameter ::  pi = 3.141592653589793
#endif

!
! internally set parameters:

real   epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
logical:: dosubsidence = .false.
real fcorz      ! Vertical Coriolis parameter
real coszrs

!----------------------------------------------
! Parameters set by PARAMETERS namelist:
! Initialized to default values.
!----------------------------------------------

real:: ug = 0.        ! Velocity of the Domain's drift in x direction
real:: vg	= 0.        ! Velocity of the Domain's drift in y direction
real:: fcor = -999.   ! Coriolis parameter	
real:: longitude0 = 0.    ! latitude of the domain's center 
real:: latitude0  = 0.    ! longitude of the domain's center 
real:: nxco2 = 1         ! factor to modify co2 concentration
logical:: doradlat = .false.
logical:: doradlon = .false.

real(kind=selected_real_kind(12)):: tabs_s =0.	! surface temperature,K
real:: delta_sst = 0.   ! amplitude of sin-pattern of sst about tabs_s (ocean_type=1)
real:: depth_slab_ocean = 2. ! thickness of the slab-ocean (m)
real:: Szero = 0.  ! mean ocean transport (W/m2)
real:: deltaS = 0. ! amplitude of linear variation of ocean transport (W/m2)
real:: timesimpleocean = 0. ! time to start simple ocean

real::   fluxt0 =0.  ! surface sensible flux, Km/s
real::   fluxq0 =0.  ! surface latent flux, m/s
real::   tau0   =0.  ! surface stress, m2/s2
real::   z0     =0.035	! roughness length
real::   soil_wetness =1.! wetness coeff for soil (from 0 to 1.)
integer:: ocean_type =0 ! type of SST forcing
logical:: cem =.false.    ! flag for Cloud Ensemble Model
logical:: les =.false.    ! flag for Large-Eddy Simulation
logical:: ocean =.false.  ! flag indicating that surface is water
logical:: land =.false.   ! flag indicating that surface is land
logical:: sfc_flx_fxd =.false. ! surface sensible flux is fixed
logical:: sfc_tau_fxd =.false.! surface drag is fixed

real:: timelargescale =0. ! time to start large-scale forcing

! nudging boundaries (between z1 and z2, where z2 > z1): 
real:: nudging_uv_z1 =-1., nudging_uv_z2 = 1000000.
real:: nudging_t_z1 =-1., nudging_t_z2 = 1000000.
real:: nudging_q_z1 =-1., nudging_q_z2 = 1000000.
real:: tauls = 99999999.    ! nudging-to-large-scaler-profile time-scale
real:: tautqls = 99999999.! nudging-to-large-scaler-profile time-scale for scalars

logical:: dodamping = .false.
logical:: doupperbound = .false. 
logical:: docloud = .false. 
logical:: doclubb = .false. ! Enabled the CLUBB parameterization (interactively)
logical:: doclubb_sfc_fluxes = .false. ! Apply the surface fluxes within the CLUBB code rather than SAM
logical:: doclubbnoninter = .false. ! Enable the CLUBB parameterization (non-interactively)
logical:: docam_sfc_fluxes = .false.   ! Apply the surface fluxes within CAM
logical:: doprecip = .false.
logical:: dolongwave = .false. 
logical:: doshortwave = .false.
logical:: dosgs = .false.
logical:: docoriolis = .false. 
logical:: docoriolisz = .false. 
logical:: dofplane = .true.
logical:: dosurface = .false. 
logical:: dolargescale = .false. 
logical:: doradforcing = .false.
logical:: dosfcforcing = .false. 
logical:: doradsimple = .false. 
logical:: donudging_uv = .false. 
logical:: donudging_tq = .false.
logical:: donudging_t = .false. 
logical:: donudging_q = .false.
logical:: doensemble = .false. 
logical:: dowallx = .false. 
logical:: dowally = .false. 
logical:: docolumn = .false. 
logical:: docup = .false.
logical:: doperpetual = .false. 
logical:: doseasons = .false. 
logical:: doradhomo = .false. 
logical:: dosfchomo = .false.
logical:: dossthomo = .false. 
logical:: dodynamicocean = .false. 
logical:: dosolarconstant = .false.
logical:: dotracers = .false. 
logical:: dosmoke = .false. 
logical:: notracegases = .false.

! Specify solar constant and zenith angle for perpetual insolation.
! Based onn Tompkins and Graig (1998)
! Note that if doperpetual=.true. and dosolarconstant=.false.
! the insolation will be set to the daily-averaged value on day0.
real:: solar_constant = 685. ! solar constant (in W/m2)
real:: zenith_angle = 51.7   ! zenith angle (in degrees)

integer:: nensemble =0   ! the number of subensemble set of perturbations
integer:: perturb_type  = 0 ! type of initial noise in setperturb()
integer:: nclubb = 1 ! SAM timesteps per CLUBB timestep
! Initial bubble parameters. Activated when perturb_type = 2
  real:: bubble_x0 = 0.
  real:: bubble_y0 = 0.
  real:: bubble_z0 = 0.
  real:: bubble_radius_hor = 0.
  real:: bubble_radius_ver = 0.
  real:: bubble_dtemp = 0.
  real:: bubble_dq = 0.

real uhl        ! current large-scale velocity in x near sfc
real vhl        ! current large-scale velocity in y near sfc
real ::  taux0 = 0.   ! surface stress in x, m2/s2
real ::  tauy0 = 0.   ! surface stress in y, m2/s2

end module crmx_params
