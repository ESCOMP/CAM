!WRF:MODEL_LAYER:CHEMICS
!
    MODULE module_data_radm2

      use shr_kind_mod,    only: r8 => shr_kind_r8

      IMPLICIT NONE
!     REAL(r8), PARAMETER ::  epsilc = 1.E-16_r8
      REAL(r8), PARAMETER ::  epsilc = 1.E-12_r8

!--- for radm solver
! .. Parameters ..
      INTEGER, PARAMETER :: ldiag = 18, lpred = 39, lss = 2, &
        lump = 4, naqre = 70,  nreacj = 21, nreack = 140, &
        ntroe = 7, numchem_radm = 41
      INTEGER, PARAMETER :: lspec = lpred + lss
      INTEGER, DIMENSION(1:NTROE) :: itroe = (/11, 22, 10, 15, 21, 24, 28/)
!
!
!
      INTEGER, PARAMETER :: lso2=1
      INTEGER, PARAMETER :: lsulf=2
      INTEGER, PARAMETER :: lno2=3
      INTEGER, PARAMETER :: lno=4
      INTEGER, PARAMETER :: lo3=5
      INTEGER, PARAMETER :: lhno3=6
      INTEGER, PARAMETER :: lh2o2=7
      INTEGER, PARAMETER :: lald=8
      INTEGER, PARAMETER :: lhcho=9
      INTEGER, PARAMETER :: lop1=10
      INTEGER, PARAMETER :: lop2=11
      INTEGER, PARAMETER :: lpaa=12
      INTEGER, PARAMETER :: lora1=13
      
      INTEGER, PARAMETER :: lora2=14
      INTEGER, PARAMETER :: lnh3=15
      INTEGER, PARAMETER :: ln2o5=16
      INTEGER, PARAMETER :: lno3=17
      INTEGER, PARAMETER :: lpan=18
      INTEGER, PARAMETER :: lhc3=19
      INTEGER, PARAMETER :: lhc5=20
      INTEGER, PARAMETER :: lhc8=21
      
      INTEGER, PARAMETER :: leth=22
      INTEGER, PARAMETER :: lco=23
      INTEGER, PARAMETER :: lol2=24
      INTEGER, PARAMETER :: lolt=25
      INTEGER, PARAMETER :: loli=26
      INTEGER, PARAMETER :: ltol=27
      INTEGER, PARAMETER :: lxyl=28
      INTEGER, PARAMETER :: laco3=29
      
      INTEGER, PARAMETER :: ltpan=30
      INTEGER, PARAMETER :: lhono=31
      INTEGER, PARAMETER :: lhno4=32
      INTEGER, PARAMETER :: lket=33
      INTEGER, PARAMETER :: lgly=34
      INTEGER, PARAMETER :: lmgly=35
      INTEGER, PARAMETER :: ldcb=36
      INTEGER, PARAMETER :: lonit=37
      
      INTEGER, PARAMETER :: lcsl=38
      INTEGER, PARAMETER :: liso=39
      INTEGER, PARAMETER :: lho=40
      INTEGER, PARAMETER :: lho2=41
! parameters for timestep, integration
      INTEGER, DIMENSION(1:lpred) :: intgrt = (/1, 1, 1, 0, 1, &
                                 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, &
                                 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, &
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                                   1, 1, 1, 1 /)
!     INTEGER, DIMENSION(1:lspec) ::  qdtc  = (/0, 0, 1, 0, 1, &
!                                0, 1, 0, 1, 0, 0, 0, 1, 0, 0, &
!                                1, 1, 1, 0, 0, 0, 0, 0, 0, 0, &
!                                0, 0, 0, 1, 1, 1, 1, 0, 0, 0, &
!                                            0, 0, 0, 0, 0, 0 /)
      INTEGER, DIMENSION(1:lspec) ::  qdtc  = (/1, 1, 1, 0, 1, &
                                 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, &
                                 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, &
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                             1, 1, 1, 1, 0, 0 /)
! max, min values,
      INTEGER :: itrdu
!
      REAL(r8), DIMENSION(1:lspec) :: cmin =(/(1.E-16_r8,itrdu=1,lspec)/)
!
      REAL(r8), DIMENSION(1:lspec) :: cmax=(/1._r8,  1._r8,  1._r8,   1._r8, .2_r8, &
                3._r8, .05_r8, .01_r8, .01_r8, .01_r8, .05_r8, .01_r8, .05_r8,  .05_r8,.05_r8, &
                1._r8,  1._r8,  1._r8,  1._r8,  1._r8,  1._r8,  1._r8,  1._r8,   1._r8, 1._r8, &
                1._r8,  1._r8,  1._r8,  1._r8,  1._r8,  1._r8,  1._r8,  1._r8,.0001_r8, .1_r8, &
                                    1._r8, .001_r8, .01_r8, .01_r8, .01_r8, .01_r8/) 

!
!
!
      INTEGER, PARAMETER :: lo3p=1
      INTEGER, PARAMETER :: lo1d=2
      INTEGER, PARAMETER :: ltco3=3
      INTEGER, PARAMETER :: lhc3p=4
      INTEGER, PARAMETER :: lhc5p=5
      INTEGER, PARAMETER :: lhc8p=6
       
      INTEGER, PARAMETER :: lol2p=7
      INTEGER, PARAMETER :: loltp=8
      INTEGER, PARAMETER :: lolip=9
      INTEGER, PARAMETER :: ltolp=10
      INTEGER, PARAMETER :: lxylp=11
      INTEGER, PARAMETER :: lethp=12
      INTEGER, PARAMETER :: lketp=13
      INTEGER, PARAMETER :: loln=14
       
      INTEGER, PARAMETER :: lxo2=15
      INTEGER, PARAMETER :: lxno2=16
      INTEGER, PARAMETER :: lxho=17
      INTEGER, PARAMETER :: lmo2=18
!
!
      INTEGER, PARAMETER ::  lnox=1
      INTEGER, PARAMETER ::  lhox=2
      INTEGER, PARAMETER ::  lpao3=3
      INTEGER, PARAMETER ::  ln2n3=4
! ..
      REAL(r8), PARAMETER ::  ch4=1.7_r8
      REAL(r8), PARAMETER ::  co2=350._r8
      REAL(r8), PARAMETER ::  n2=7.81E5_r8
      REAL(r8), PARAMETER ::  o2=2.09E5_r8
      REAL(r8), PARAMETER ::  pi=3.141592654_r8

! ..
      REAL(r8) :: afac(2), &
        bfac(2),  const(3), eor(nreack), &
        thafac(nreack), &
        xk0300(ntroe), &
        xkf300(ntroe), xmtroe(ntroe), xntroe(ntroe)

! ..
! .. Data Statements ..
      DATA thafac/0.00_r8, 6.50E-12_r8, 1.80E-11_r8, 3.20E-11_r8, 2.20E-10_r8, 2.00E-12_r8, &
        1.60E-12_r8, 1.10E-14_r8, 3.70E-12_r8, 4*0.00_r8, 3.30E-12_r8, 0.00_r8, 3.30E-19_r8, &
        1.40E-13_r8, 1.70E-11_r8, 2.50E-14_r8, 2.50E-12_r8, 2*0.00_r8, 2.00E-21_r8, 2*0.00_r8, &
        1.30E-12_r8, 4.60E-11_r8, 2*0.00_r8, 6.95E-18_r8, 1.37E-17_r8, 1.59E-11_r8, 1.73E-11_r8, &
        3.64E-11_r8, 2.15E-12_r8, 5.32E-12_r8, 1.07E-11_r8, 2.10E-12_r8, 1.89E-11_r8, 4.00E-11_r8, &
        9.00E-12_r8, 6.87E-12_r8, 1.20E-11_r8, 1.15E-11_r8, 1.70E-11_r8, 2.80E-11_r8, 1.00E-11_r8, &
        1.00E-11_r8, 1.00E-11_r8, 6.85E-18_r8, 1.55E-11_r8, 2.55E-11_r8, 2.80E-12_r8, 1.95E+16_r8, &
        4.70E-12_r8, 1.95E+16_r8, 4.20E-12_r8, 4.20E-12_r8, 0.00_r8, 4.20E-12_r8, 0.00_r8, &
        4.20E-12_r8, 0.00_r8, 10*4.20E-12_r8, 6.00E-13_r8, 1.40E-12_r8, 6.00E-13_r8, 1.40E-12_r8, &
        1.40E-12_r8, 2.20E-11_r8, 2.00E-12_r8, 1.00E-11_r8, 3.23E-11_r8, 5.81E-13_r8, 1.20E-14_r8, &
        1.32E-14_r8, 7.29E-15_r8, 1.23E-14_r8, 14*7.70E-14_r8, 1.90E-13_r8, 1.40E-13_r8, &
        4.20E-14_r8, 3.40E-14_r8, 2.90E-14_r8, 1.40E-13_r8, 1.40E-13_r8, 1.70E-14_r8, 1.70E-14_r8, &
        9.60E-13_r8, 1.70E-14_r8, 1.70E-14_r8, 9.60E-13_r8, 3.40E-13_r8, 1.00E-13_r8, 8.40E-14_r8, &
        7.20E-14_r8, 3.40E-13_r8, 3.40E-13_r8, 4.20E-14_r8, 4.20E-14_r8, 1.19E-12_r8, 4.20E-14_r8, &
        4.20E-14_r8, 1.19E-12_r8, 7.70E-14_r8, 1.70E-14_r8, 4.20E-14_r8, 3.60E-16_r8, 4.20E-12_r8, &
        4.20E-12_r8, 7.70E-14_r8, 1.70E-14_r8, 4.20E-14_r8, 3.60E-16_r8, 0.00_r8, 1.70E-14_r8, &
        4.20E-14_r8, 3.60E-16_r8/
! ..
! constants for RADM2 rate coefficients
      DATA eor/0._r8, -120._r8, -110._r8, -70._r8, 0._r8, 1400._r8, 940._r8, 500._r8, -240._r8, 0._r8, 0._r8, &
        0._r8, 0._r8, 200._r8, 0._r8, -530._r8, 2500._r8, -150._r8, 1230._r8, 0._r8, 0._r8, 0._r8, 0._r8, 0._r8, 0._r8, &
        -380._r8, -230._r8, 0._r8, 0._r8, 1280._r8, 444._r8, 540._r8, 380._r8, 380._r8, -411._r8, -504._r8, &
        -549._r8, -322._r8, -116._r8, 0._r8, 0._r8, -256._r8, 745._r8, 0._r8, 0._r8, 0._r8, 0._r8, 0._r8, 0._r8, &
        444._r8, 540._r8, -409._r8, -181._r8, 13543._r8, 0._r8, 13543._r8, -180._r8, -180._r8, 0._r8, -180._r8, &
        0._r8, -180._r8, 0._r8, -180._r8, -180._r8, -180._r8, -180._r8, -180._r8, -180._r8, -180._r8, -180._r8, &
        -180._r8, -180._r8, 2058._r8, 1900._r8, 2058._r8, 1900._r8, 1900._r8, 0._r8, 2923._r8, 1895._r8, &
        975._r8, 0._r8, 2633._r8, 2105._r8, 1136._r8, 2013._r8, -1300._r8, -1300._r8, -1300._r8, -1300._r8, &
        -1300._r8, -1300._r8, -1300._r8, -1300._r8, -1300._r8, -1300._r8, -1300._r8, -1300._r8, &
        -1300._r8, -1300._r8, 25* -220._r8, -1300._r8, -220._r8, -220._r8, -220._r8, -180._r8, -180._r8, &
        -1300._r8, -220._r8, -220._r8, 0._r8, 0._r8, -220._r8, -220._r8, -220._r8/

      DATA xk0300/1.8E-31_r8, 2.2E-30_r8, 1.8E-31_r8, 7.E-31_r8, 2.2E-30_r8, 2.6E-30_r8, 3.E-31_r8/
      DATA xntroe/3.2_r8, 4.3_r8, 3.2_r8, 2.6_r8, 4.3_r8, 3.2_r8, 3.3_r8/
      DATA xkf300/4.7E-12_r8, 1.5E-12_r8, 4.7E-12_r8, 1.5E-11_r8, 1.5E-12_r8, 2.4E-11_r8, &
        1.5E-12_r8/
      DATA xmtroe/1.4_r8, 0.5_r8, 1.4_r8, 2*.5_r8, 1.3_r8, 0._r8/
      DATA afac/2.1E-27_r8, 1.1E-27_r8/
      DATA bfac/10900._r8, 11200._r8/
      DATA const/7.34E21_r8, 4.4E17_r8, 3.23E33_r8/

    END MODULE module_data_radm2
