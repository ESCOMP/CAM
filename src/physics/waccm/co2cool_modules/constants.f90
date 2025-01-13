!=======================================================================!
! CONSTANTS ---  module with parameters                                 !
!-----------------------------------------------------------------------!
! STATUS: bf 23-Feb-2024                      CREATED: bf 23-Feb-2024   !
!-----------------------------------------------------------------------!
! COPYRIGHT: (C) 2023-2024 Instituto de Astrofisica de Andalucia (CSIC) !
!-----------------------------------------------------------------------!
! LICENSE: GNU Lesser General Public License (LGPL) version 2.1         !
!          (see file LICENSE_LGPLv2.1)                                  !
!-----------------------------------------------------------------------!
! CONTAINS:                                                             !
! public:                                                               !
!                                                                       !
!-----------------------------------------------------------------------!
! MAINTENANCE HISTORY:                                                  !
! MLP: 14-March-2024: Added some notes                                  !
!=======================================================================!

module constants
  use precision

! constants
!----------
real(dp), parameter :: kboltz  = 1.38064853e-23_dp ! Boltzmann constant [J K-1]
real(dp), parameter :: hplanck = 6.62606876e-19_dp ! Planck constant [J * s]
real(dp), parameter :: clight  = 2.99792458e+08_dp ! vacuum speed of light [m * s-1]
real(dp), parameter :: Nav     = 6.02214076e+23_dp ! Avogadro constant [mol-1]
real(db),parameter  :: rgasuni = 8.31441_dp        ! molar gas constant R [J/(mol*K)]
real(dp), parameter :: E_fun   = 667.3799_dp       ! Energy of CO2 fundamental (01101-00001) band
real(dp), parameter :: cp_0    = 1.005e7_dp        ! specific enthalpy dry air - erg g-1 K-1
real(dp), parameter :: sperday = 86400.0_dp        ! seconds per day
real(dp), parameter :: mo2 = 32.0_db, mo = 16.0_dp, mn2 = 28.0_dp !molecular masses

! derived constants
!------------------
real(dp), parameter :: kb      = kboltz*1e4_dp
real(dp), parameter :: kbc     = kboltz/(hplanck*clight*1e-13_dp) !  0.695035
real(dp), parameter :: Edkbc   = E_fun/kbc
real(dp), parameter :: Nave1   = 1e1_dp * Nav
real(dp), parameter :: rgasue7 = rgasuni * 1e7_db

!computational parameters
!------------------------
real(dp), parameter :: numfac    = 2.55520997e11_dp
real(dp), parameter :: lmbfac    = 1.5988_dp
real(dp), parameter :: onedthree = one/three
real(db), parameter :: p0        = 1e3_dp

! Collisional rate constants in the form of a*sqrt(T)+b*exp(-g*T**(-1./3))
!-------------------------------------------------------------------------
! Nominal values: Funke et al., GRANADA model, JQSRT, 2012
real(dp), parameter :: a_zo  = 3.5e-13_dp, b_zo  = 2.32e-09_dp, g_zo  = 76.75_dp 
real(dp), parameter :: a_zn2 = 7.0e-17_dp, b_zn2 = 6.70e-10_dp, g_zn2 = 83.80_dp 
real(dp), parameter :: a_zo2 = 7.0e-17_dp, b_zo2 = 1.00e-09_dp, g_zo2 = 83.80_dp

! dimensions
!-----------
integer, parameter :: n_co2prof  = 8  ! number of co2 vmr profiles
integer, parameter :: i_co2ref   = 3  ! reference co2 vmr level for overhead column calculation
integer, parameter :: n_alpha    = 15 ! = n_lev_ru - n_lev_rf + 1
integer, parameter :: n_Lesc     = 61 ! number of points of Lesc(uco2)
integer, parameter :: n_meanco2  = 40 ! upper level for average co2
integer, parameter :: n_lev_rf   = 52 ! lower level recurrence formula (has to be smaller than n_lev_cm)
integer, parameter :: n_lev_cm   = 55 ! upper level curtis matrix
integer, parameter :: n_lev_ru   = 66 ! upper level region 2 of recurrence formula
integer, parameter :: n_lev_cs   = 81 ! lower level cooling to space
integer, parameter :: n_lev_mx   = 83 ! max level

end module constants
