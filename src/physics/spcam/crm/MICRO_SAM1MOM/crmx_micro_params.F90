module crmx_micro_params

use crmx_grid, only: nzm

implicit none

!  Microphysics stuff:

! Densities of hydrometeors

real, parameter :: rhor = 1000. ! Density of water, kg/m3
real, parameter :: rhos = 100.  ! Density of snow, kg/m3
real, parameter :: rhog = 400.  ! Density of graupel, kg/m3
!real, parameter :: rhog = 917.  ! hail - Lin 1983    

! Temperatures limits for various hydrometeors

real, parameter :: tbgmin = 253.16    ! Minimum temperature for cloud water., K
real, parameter :: tbgmax = 273.16    ! Maximum temperature for cloud ice, K
real, parameter :: tprmin = 268.16    ! Minimum temperature for rain, K
real, parameter :: tprmax = 283.16    ! Maximum temperature for snow+graupel, K
real, parameter :: tgrmin = 223.16    ! Minimum temperature for snow, K
real, parameter :: tgrmax = 283.16    ! Maximum temperature for graupel, K

! Terminal velocity coefficients

real, parameter :: a_rain = 842. ! Coeff.for rain term vel 
real, parameter :: b_rain = 0.8  ! Fall speed exponent for rain
real, parameter :: a_snow = 4.84 ! Coeff.for snow term vel
real, parameter :: b_snow = 0.25 ! Fall speed exponent for snow
!real, parameter :: a_grau = 40.7! Krueger (1994) ! Coef. for graupel term vel
real, parameter :: a_grau = 94.5 ! Lin (1983) (rhog=400)
!real, parameter :: a_grau = 127.94! Lin (1983) (rhog=917)
real, parameter :: b_grau = 0.5  ! Fall speed exponent for graupel

! Autoconversion
#ifdef CLUBB_CRM   /*microphysical tuning for CLUBB*/
real, parameter :: qcw0 = 0.6e-3      ! Threshold for water autoconversion, g/g  
real, parameter :: qci0 = 1.e-4     ! Threshold for ice autoconversion, g/g
real, parameter :: alphaelq = 10.e-3  ! autoconversion of cloud water rate coef
real, parameter :: betaelq = 6.0e-3   ! autoconversion of cloud ice rate coef
#else 
real, parameter :: qcw0 = 1.e-3      ! Threshold for water autoconversion, g/g  
real, parameter :: qci0 = 1.e-4     ! Threshold for ice autoconversion, g/g
real, parameter :: alphaelq = 1.e-3  ! autoconversion of cloud water rate coef
real, parameter :: betaelq = 1.e-3   ! autoconversion of cloud ice rate coef
#endif /*CLUBB_CRM*/

! Accretion

real, parameter :: erccoef = 1.0   ! Rain/Cloud water collection efficiency
real, parameter :: esccoef = 1.0   ! Snow/Cloud water collection efficiency
real, parameter :: esicoef = 0.1   ! Snow/cloud ice collection efficiency
real, parameter :: egccoef = 1.0   ! Graupel/Cloud water collection efficiency
real, parameter :: egicoef = 0.1   ! Graupel/Cloud ice collection efficiency

! Interseption parameters for exponential size spectra

real, parameter :: nzeror = 8.e6   ! Intercept coeff. for rain  
real, parameter :: nzeros = 3.e6   ! Intersept coeff. for snow
real, parameter :: nzerog = 4.e6   ! Intersept coeff. for graupel
!real, parameter :: nzerog = 4.e4   ! hail - Lin 1993 

real, parameter :: qp_threshold = 1.e-8 ! minimal rain/snow water content


! Misc. microphysics variables

real*4 gam3       ! Gamma function of 3
real*4 gams1      ! Gamma function of (3 + b_snow)
real*4 gams2      ! Gamma function of (5 + b_snow)/2
real*4 gams3      ! Gamma function of (4 + b_snow)
real*4 gamg1      ! Gamma function of (3 + b_grau)
real*4 gamg2      ! Gamma function of (5 + b_grau)/2
real*4 gamg3      ! Gamma function of (4 + b_grau)
real*4 gamr1      ! Gamma function of (3 + b_rain)
real*4 gamr2      ! Gamma function of (5 + b_rain)/2
real*4 gamr3      ! Gamma function of (4 + b_rain)
      
real accrsc(nzm),accrsi(nzm),accrrc(nzm),coefice(nzm)
real accrgc(nzm),accrgi(nzm)
real evaps1(nzm),evaps2(nzm),evapr1(nzm),evapr2(nzm)
real evapg1(nzm),evapg2(nzm)
            
real a_bg, a_pr, a_gr 


end module crmx_micro_params
