!===============================================================================
! Aerosol dry deposition
!   Portable science routines split from modal_aero/aero_model.F90 and
!   aer_drydep_mod.F90: surface deposition velocities of particles
!   (Zhang et al. 2001) and the aerodynamic resistance / friction velocity
!   patch over ocean and sea ice. Host constants and the landuse fractions
!   are passed as arguments; array sizing is by ncol/pver runtime arguments.
!===============================================================================
module aero_drydep_core

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: modal_aero_depvel_part
  public :: calcram

contains

  !=============================================================================
  !=============================================================================
  subroutine modal_aero_depvel_part( ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv,  &
                                     radius_part, density_part, sig_part, moment, &
                                     pver, top_lev, n_land_type, fraction_landuse, &
                                     pi, boltz, gravit, rair, aspherical )   ! dmleung added aspherical flag 20 Oct 2025

!    calculates surface deposition velocity of particles
!    L. Zhang, S. Gong, J. Padro, and L. Barrie
!    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
!    Atmospheric Environment, 35, 549-560, 2001.
!
!    Authors: X. Liu

    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8), intent(in) :: t(:,:)       !atm temperature (K)
    real(r8), intent(in) :: pmid(:,:)    !atm pressure (Pa)
    real(r8), intent(in) :: fv(:)           !friction velocity (m/s)
    real(r8), intent(in) :: ram1(:)         !aerodynamical resistance (s/m)
    real(r8), intent(in) :: radius_part(:,:)    ! mean (volume/number) particle radius (m)
    real(r8), intent(in) :: density_part(:,:)   ! density of particle material (kg/m3)
    real(r8), intent(in) :: sig_part(:,:)       ! geometric standard deviation of particles
    integer,  intent(in) :: moment ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: pver                ! number of vertical levels
    integer,  intent(in) :: top_lev             ! top level for modal aerosols
    integer,  intent(in) :: n_land_type         ! number of land use types
    real(r8), intent(in) :: fraction_landuse(:,:)  ! land use fractions (ncol, n_land_type)
    real(r8), intent(in) :: pi                  ! host model constants
    real(r8), intent(in) :: boltz               ! Boltzmann constant (J/K)
    real(r8), intent(in) :: gravit              ! gravitational acceleration (m/s2)
    real(r8), intent(in) :: rair                ! gas constant for dry air (J/K/kg)

    real(r8), intent(out) :: vlc_trb(:)       !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(:,:)       !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(:,:)       !dry deposn velocity (m/s)
    logical,  intent(in), OPTIONAL :: aspherical   ! dmleung: asphericity is strong for coarse-mode interstitial
    ! aerosols only, mostly dust and seasalt. For coarse mode aerosols, asphericity reduces coarse-mode gravitational
    ! settling velocity by 20 % following Fig. 4 of Yue Huang et al. (2020).
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k,ix                !indices
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8) :: vsc_dyn_atm(ncol,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(ncol,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(ncol,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: slp_crc(ncol,pver) ![frc] Slip correction factor
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: rss_lmn       ![s m-1] Quasi-laminar layer resistance
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface
    real(r8) :: radius_moment(ncol,pver) ! median radius (m) for moment
    real(r8) :: lnsig         ! ln(sig_part)
    real(r8) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
                              ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    integer  :: lt
    real(r8) :: lnd_frc
    real(r8) :: wrk1, wrk2, wrk3

    ! constants

     real(r8), parameter :: asphericaldust_drydep = 0.8_r8 ! dmleung added 20 Oct 2025: aspherical dust reduces
     ! gravitational settling velocity by 15-20 %. Yue Huang et al. (2020)
     ! Climate Models and Remote Sensing Retrievals Neglect Substantial Desert Dust Asphericity

    real(r8) gamma(11)      ! exponent of schmidt number
!   data gamma/0.54d+00,  0.56d+00,  0.57d+00,  0.54d+00,  0.54d+00, &
!              0.56d+00,  0.54d+00,  0.54d+00,  0.54d+00,  0.56d+00, &
!              0.50d+00/
    data gamma/0.56e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.56e+00_r8,  0.56e+00_r8, &
               0.56e+00_r8,  0.50e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.54e+00_r8, &
               0.54e+00_r8/
    save gamma

    real(r8) alpha(11)      ! parameter for impaction
!   data alpha/50.00d+00,  0.95d+00,  0.80d+00,  1.20d+00,  1.30d+00, &
!               0.80d+00, 50.00d+00, 50.00d+00,  2.00d+00,  1.50d+00, &
!             100.00d+00/
    data alpha/1.50e+00_r8,   1.20e+00_r8,  1.20e+00_r8,  0.80e+00_r8,  1.00e+00_r8, &
               0.80e+00_r8, 100.00e+00_r8, 50.00e+00_r8,  2.00e+00_r8,  1.20e+00_r8, &
              50.00e+00_r8/
    save alpha

    real(r8) radius_collector(11) ! radius (m) of surface collectors
!   data radius_collector/-1.00d+00,  5.10d-03,  3.50d-03,  3.20d-03, 10.00d-03, &
!                          5.00d-03, -1.00d+00, -1.00d+00, 10.00d-03, 10.00d-03, &
!                         -1.00d+00/
    data radius_collector/10.00e-03_r8,  3.50e-03_r8,  3.50e-03_r8,  5.10e-03_r8,  2.00e-03_r8, &
                           5.00e-03_r8, -1.00e+00_r8, -1.00e+00_r8, 10.00e-03_r8,  3.50e-03_r8, &
                          -1.00e+00_r8/
    save radius_collector

    integer            :: iwet(11) ! flag for wet surface = 1, otherwise = -1
!   data iwet/1,   -1,   -1,   -1,   -1,  &
!            -1,   -1,   -1,    1,   -1,  &
!             1/
    data iwet/-1,  -1,   -1,   -1,   -1,  &
              -1,   1,   -1,    1,   -1,  &
              -1/
    save iwet


    vlc_trb = 0._r8
    vlc_grv = 0._r8
    vlc_dry = 0._r8

    !------------------------------------------------------------------------
    do k=top_lev,pver ! radius_part is not defined above top_lev
       do i=1,ncol

          lnsig = log(sig_part(i,k))
! use a maximum radius of 50 microns when calculating deposition velocity
          radius_moment(i,k) = min(50.0e-6_r8,radius_part(i,k))*   &
                          exp((float(moment)-1.5_r8)*lnsig*lnsig)
          dispersion = exp(2._r8*lnsig*lnsig)

          rho=pmid(i,k)/rair/t(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air

          slp_crc(i,k) = 1.0_r8 + mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*radius_moment(i,k)/(mfp_atm(i,k)))) / &
                  radius_moment(i,k)   ![frc] Slip correction factor SeP97 p. 464
          vlc_grv(i,k) = (4.0_r8/18.0_r8) * radius_moment(i,k)*radius_moment(i,k)*density_part(i,k)* &
                  gravit*slp_crc(i,k) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(i,k) = vlc_grv(i,k) * dispersion

          ! dmleung edited 20 Oct 2025 based on Longlei Li's edits ++
          ! asphericity reduces gravitational settling velocity of coarse-mode aerosols by 20 %.
          ! scale flag is only true for coarse mode (m == n_coarse_dust).
          if (present(aspherical)) then
             if(aspherical) then
                vlc_grv(i,k) = vlc_grv(i,k) * asphericaldust_drydep
             end if
          end if
          ! dmleung --

          vlc_dry(i,k)=vlc_grv(i,k)
       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do i=1,ncol
       dff_aer = boltz * t(i,k) * slp_crc(i,k) / &    ![m2 s-1]
                 (6.0_r8*pi*vsc_dyn_atm(i,k)*radius_moment(i,k)) !SeP97 p.474
       shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972

       wrk2 = 0._r8
       wrk3 = 0._r8
       do lt = 1,n_land_type
          lnd_frc = fraction_landuse(i,lt)
          if ( lnd_frc /= 0._r8 ) then
             brownian = shm_nbr**(-gamma(lt))
             if (radius_collector(lt) > 0.0_r8) then
!       vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) / (gravit*radius_collector(lt))
                interception = 2.0_r8*(radius_moment(i,k)/radius_collector(lt))**2.0_r8
             else
!       non-vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))  ![frc] SeP97 p.965
                interception = 0.0_r8
             endif
             impaction = (stk_nbr/(alpha(lt)+stk_nbr))**2.0_r8

             if (iwet(lt) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = exp(-sqrt(stk_nbr))
                if (stickfrac < 1.0e-10_r8) stickfrac = 1.0e-10_r8
             endif
             rss_lmn = 1.0_r8 / (3.0_r8 * fv(i) * stickfrac * (brownian+interception+impaction))
             rss_trb = ram1(i) + rss_lmn + ram1(i)*rss_lmn*vlc_grv(i,k)

             wrk1 = 1.0_r8 / rss_trb
             wrk2 = wrk2 + lnd_frc*( wrk1 )
             wrk3 = wrk3 + lnd_frc*( wrk1 + vlc_grv(i,k) )
          endif
       enddo  ! n_land_type
       vlc_trb(i) = wrk2
       vlc_dry(i,k) = wrk3
    enddo !ncol

    return
  end subroutine modal_aero_depvel_part

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Calcram
!
! !INTERFACE:
!

      subroutine  calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
           ustar,ram1in,ram1,t,pmid,&
           pdel,fvin,fv,rair,gravit)
        !
        ! !DESCRIPTION:
        !
        ! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
        ! from Seinfeld and Pandis, p.963.
        !
        ! Author: Natalie Mahowald
        !
        implicit none
        integer, intent(in) :: ncol
        real(r8),intent(in) :: ram1in(:)         !aerodynamical resistance (s/m)
        real(r8),intent(in) :: fvin(:)                 ! sfc frc vel from land
        real(r8),intent(out) :: ram1(:)         !aerodynamical resistance (s/m)
        real(r8),intent(out) :: fv(:)                 ! sfc frc vel from land
        real(r8), intent(in) :: obklen(:)                 ! obklen
        real(r8), intent(in) :: ustar(:)                  ! sfc fric vel
        real(r8), intent(in) :: landfrac(:)               ! land fraction
        real(r8), intent(in) :: icefrac(:)                ! ice fraction
        real(r8), intent(in) :: ocnfrac(:)                ! ocean fraction
        real(r8), intent(in) :: t(:)       !atm temperature (K)
        real(r8), intent(in) :: pmid(:)    !atm pressure (Pa)
        real(r8), intent(in) :: pdel(:)    !atm pressure (Pa)
        real(r8), intent(in) :: rair       ! gas constant for dry air (J/K/kg)
        real(r8), intent(in) :: gravit     ! gravitational acceleration (m/s2)
        real(r8), parameter :: zzocen = 0.0001_r8   ! Ocean aerodynamic roughness length
        real(r8), parameter :: zzsice = 0.0400_r8   ! Sea ice aerodynamic roughness length
        real(r8), parameter :: xkar   = 0.4_r8      ! Von Karman constant

        ! local variables
        real(r8) :: z,psi,psi0,nu,nu0,temp,ram
        integer :: i
        !    write(iulog,*) rair,zzsice,zzocen,gravit,xkar


        do i=1,ncol
           z=pdel(i)*rair*t(i)/pmid(i)/gravit/2.0_r8   !use half the layer height like Ganzefeld and Lelieveld, 1995
           if(obklen(i).eq.0) then
              psi=0._r8
              psi0=0._r8
           else
              psi=min(max(z/obklen(i),-1.0_r8),1.0_r8)
              psi0=min(max(zzocen/obklen(i),-1.0_r8),1.0_r8)
           endif
           temp=z/zzocen
           if(icefrac(i) > 0.5_r8) then 
              if(obklen(i).gt.0) then 
                 psi0=min(max(zzsice/obklen(i),-1.0_r8),1.0_r8)
              else
                 psi0=0.0_r8
              endif
              temp=z/zzsice
	   endif
           if(psi> 0._r8) then
              ram=1/xkar/ustar(i)*(log(temp)+4.7_r8*(psi-psi0))
           else
              nu=(1.00_r8-15.000_r8*psi)**(.25_r8)
              nu0=(1.000_r8-15.000_r8*psi0)**(.25_r8)
              if(ustar(i).ne.0._r8) then
                 ram=1/xkar/ustar(i)*(log(temp) &
                      +log(((nu0**2+1.00_r8)*(nu0+1.0_r8)**2)/((nu**2+1.0_r8)*(nu+1.00_r8)**2)) &
                      +2.0_r8*(atan(nu)-atan(nu0)))
              else
	         ram=0._r8
              endif
           endif
           if(landfrac(i) < 0.000000001_r8) then
              fv(i)=ustar(i)
              ram1(i)=ram
           else
              fv(i)=fvin(i)
              ram1(i)=ram1in(i)
           endif
           !          write(iulog,*) i,pdel(i),t(i),pmid(i),gravit,obklen(i),psi,psi0,icefrac(i),nu,nu0,ram,ustar(i),&
           !             log(((nu0**2+1.00)*(nu0+1.0)**2)/((nu**2+1.0)*(nu+1.00)**2)),2.0*(atan(nu)-atan(nu0))

        enddo

        ! fvitt -- fv == 0 causes a floating point exception in 
        ! dry dep of sea salts and dust
        where ( fv(:ncol) == 0._r8 ) 
           fv(:ncol) = 1.e-12_r8
        endwhere

        return
      end subroutine calcram

end module aero_drydep_core
