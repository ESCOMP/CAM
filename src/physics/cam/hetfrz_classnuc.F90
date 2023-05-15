module hetfrz_classnuc

!-----------------------------------------------------------------------
!
! Purpose: Calculate heterogeneous freezing rates from classical nucleation theory
!
! Public interfaces:
!
! hetfrz_classnuc_init
! hetfrz_classnuc_calc
!
! Author:
!   Corinna Hoose, UiO, May 2009
!   Yong Wang and Xiaohong Liu, UWyo, 12/2012,
!   implement in CAM5 and constrain uncertain parameters using natural dust and
!   BC(soot) datasets.
!   Yong Wang and Xiaohong Liu, UWyo, 05/2013, implement the PDF-contact angle approach:
!   Y. Wang et al., Atmos. Chem. Phys., 2014. https://doi.org/10.5194/acp-14-10411-2014
!   Jack Chen, NCAR, 09/2015, modify calculation of dust activation fraction.
!
!-----------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use wv_saturation, only: svp_water
use shr_spfn_mod,  only: erf => shr_spfn_erf

use physconst,     only: pi, planck, boltz, mwso4, amu, pstd

implicit none
private

public :: hetfrz_classnuc_init, hetfrz_classnuc_calc

real(r8) :: rair
real(r8) :: cpair
real(r8) :: rh2o
real(r8) :: rhoh2o
real(r8) :: mwh2o
real(r8) :: tmelt

!*****************************************************************************
!                PDF theta model
!*****************************************************************************
! some variables for PDF theta model
! immersion freezing
!
! With the original value of pdf_n_theta set to 101 the dust activation
! fraction between -15 and 0 C could be overestimated.  This problem was
! eliminated by increasing pdf_n_theta to 301.  To reduce the expense of
! computing the dust activation fraction the integral is only evaluated
! where dim_theta is non-zero.  This was determined to be between
! dim_theta index values of 53 through 113.  These loop bounds are
! hardcoded in the variables i1 and i2.

logical, parameter :: pdf_imm_in = .true.
integer, parameter :: pdf_n_theta = 301
integer, parameter :: i1 = 53
integer, parameter :: i2 = 113

real(r8) :: dim_theta(pdf_n_theta) = -huge(1._r8)
real(r8) :: pdf_imm_theta(pdf_n_theta) = 0.0_r8
real(r8) :: pdf_d_theta = -huge(1._r8)
real(r8) :: dim_f_imm(pdf_n_theta) = 0.0_r8

integer  :: iulog

real(r8), parameter :: n1 = 1.e19_r8           ! number of water molecules in contact with unit area of substrate [m-2]
real(r8), parameter :: rhplanck = 1._r8/planck
real(r8), parameter :: nus = 1.e13_r8          ! frequ. of vibration [s-1] higher freq. (as in P&K, consistent with Anupam's data)
real(r8), parameter :: rhwincloud = 0.98_r8    ! 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)

logical, parameter :: tot_in = .false.

real(r8) :: bc_limfac = -huge(1._r8)   ! soot ice nucleating fraction
real(r8) :: dust_limfac = -huge(1._r8) ! dust ice nucleating fraction

!===================================================================================================
contains
!===================================================================================================

subroutine hetfrz_classnuc_init( &
   rair_in, cpair_in, rh2o_in, rhoh2o_in, mwh2o_in, &
   tmelt_in, iulog_in, bc_limfac_in, dust_limfac_in)

   real(r8), intent(in) :: rair_in
   real(r8), intent(in) :: cpair_in
   real(r8), intent(in) :: rh2o_in
   real(r8), intent(in) :: rhoh2o_in
   real(r8), intent(in) :: mwh2o_in
   real(r8), intent(in) :: tmelt_in
   integer,  intent(in) :: iulog_in
   real(r8), intent(in) :: bc_limfac_in
   real(r8), intent(in) :: dust_limfac_in

   rair   = rair_in
   cpair  = cpair_in
   rh2o   = rh2o_in
   rhoh2o = rhoh2o_in
   mwh2o  = mwh2o_in
   tmelt  = tmelt_in
   iulog  = iulog_in

   bc_limfac = bc_limfac_in
   dust_limfac = dust_limfac_in

   ! Initialize all the PDF theta variables:
   if (pdf_imm_in) then
      call hetfrz_classnuc_init_pdftheta()
   end if

end subroutine hetfrz_classnuc_init

!===================================================================================================

subroutine hetfrz_classnuc_calc(ntypes, types,&
   deltat, T, p, supersatice,                 &
   fn,                                        &
   r3lx, icnlx,                               &
   frzbcimm, frzduimm,                        &
   frzbccnt, frzducnt,                        &
   frzbcdep, frzdudep,                        &
   hetraer, wact_factor, dstcoat,             &
   total_aer_num, uncoated_aer_num,           &
   total_interstitial_aer_num, total_cloudborne_aer_num, errstring)

   integer, intent(in) :: ntypes
   character(len=*), intent(in) :: types(ntypes)
   real(r8), intent(in) :: deltat            ! timestep [s]
   real(r8), intent(in) :: T                 ! temperature [K]
   real(r8), intent(in) :: p                 ! pressure [Pa]
   real(r8), intent(in) :: supersatice       ! supersaturation ratio wrt ice at 100%rh over water [ ]
   real(r8), intent(in) :: r3lx              ! volume mean drop radius [m]
   real(r8), intent(in) :: icnlx             ! in-cloud droplet concentration [cm-3]

   real(r8), intent(in) :: fn(ntypes)               ! fraction activated [ ] for cloud borne aerosol number
                                                    ! index values are 1:bc, 2:dust_a1, 3:dust_a3
   real(r8), intent(in) :: hetraer(ntypes)          ! bc and dust mass mean radius [m]
   real(r8), intent(in) :: wact_factor(ntypes)      ! water activity factor -- density*(1.-(OC+BC)/(OC+BC+SO4)) [mug m-3]
   real(r8), intent(in) :: dstcoat(ntypes)          ! coated fraction
   real(r8), intent(in) :: total_aer_num(ntypes)    ! total bc and dust number concentration(interstitial+cloudborne) [#/cm^3]
   real(r8), intent(in) :: uncoated_aer_num(ntypes) ! uncoated bc and dust number concentration(interstitial)
   real(r8), intent(in) :: total_interstitial_aer_num(ntypes) ! total bc and dust concentration(interstitial)
   real(r8), intent(in) :: total_cloudborne_aer_num(ntypes)   ! total bc and dust concentration(cloudborne)

   real(r8), target, intent(out) :: frzbcimm           ! het. frz by BC immersion nucleation [cm-3 s-1]
   real(r8), target, intent(out) :: frzduimm           ! het. frz by dust immersion nucleation [cm-3 s-1]
   real(r8), target, intent(out) :: frzbccnt           ! het. frz by BC contact nucleation [cm-3 s-1]
   real(r8), target, intent(out) :: frzducnt           ! het. frz by dust contact nucleation [cm-3 s-1]
   real(r8), target, intent(out) :: frzbcdep           ! het. frz by BC deposition nucleation [cm-3 s-1]
   real(r8), target, intent(out) :: frzdudep           ! het. frz by dust deposition nucleation [cm-3 s-1]

   character(len=*), intent(out) :: errstring

   ! local variables

   real(r8) :: tc
   real(r8) :: vwice
   real(r8) :: rhoice
   real(r8) :: sigma_iw                        ! [J/m2]
   real(r8) :: sigma_iv                        ! [J/m2]
   real(r8) :: eswtr                           ! [Pa]

   real(r8) :: rgimm    ! critical germ size
   real(r8) :: rgdep
   real(r8) :: dg0dep   ! homogeneous energy of germ formation
   real(r8) :: dg0cnt
   real(r8) :: Adep     ! prefactors
   real(r8) :: Acnt

   !********************************************************
   ! Hoose et al., 2010 fitting parameters
   !********************************************************
   !freezing parameters for immersion freezing
   !real(r8),parameter :: theta_imm_bc = 40.17         ! contact angle [deg], converted to rad later
   !real(r8),parameter :: dga_imm_bc = 14.4E-20        ! activation energy [J]
   !real(r8),parameter :: theta_imm_dust = 30.98       ! contact angle [deg], converted to rad later
   !real(r8),parameter :: dga_imm_dust = 15.7E-20      ! activation energy [J]
   !freezing parameters for deposition nucleation
   !real(r8),parameter :: theta_dep_dust = 12.7       ! contact angle [deg], converted to rad later !Zimmermann et al (2008), illite
   !real(r8),parameter :: dga_dep_dust = -6.21E-21    ! activation energy [J]
   !real(r8),parameter :: theta_dep_bc = 28.          ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
   !real(r8),parameter :: dga_dep_bc = -2.E-19        ! activation energy [J]
   !********************************************************
   ! Wang et al., 2014 fitting parameters
   !********************************************************
   ! freezing parameters for immersion freezing
   real(r8),parameter :: theta_imm_bc = 48.0_r8            ! contact angle [deg], converted to rad later !DeMott et al (1990)
   real(r8),parameter :: dga_imm_bc = 14.15E-20_r8         ! activation energy [J]
   real(r8),parameter :: theta_imm_dust = 46.0_r8          ! contact angle [deg], converted to rad later !DeMott et al (2011) SD
   real(r8),parameter :: dga_imm_dust = 14.75E-20_r8       ! activation energy [J]
   ! freezing parameters for deposition nucleation
   real(r8),parameter :: theta_dep_dust = 20.0_r8          ! contact angle [deg], converted to rad later !Koehler et al (2010) SD
   real(r8),parameter :: dga_dep_dust = -8.1E-21_r8        ! activation energy [J]
   real(r8),parameter :: theta_dep_bc = 28._r8             ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
   real(r8),parameter :: dga_dep_bc = -2.E-19_r8           ! activation energy [J]

   ! form factor
   ! only consider flat surfaces due to uncertainty of curved surfaces
   real(r8),parameter :: m_depcnt_bc = COS(theta_dep_bc*pi/180._r8)
   real(r8),parameter :: f_depcnt_bc = (2+m_depcnt_bc)*(1-m_depcnt_bc)**2/4._r8

   real(r8),parameter :: m_depcnt_dst = COS(theta_dep_dust*pi/180._r8)
   real(r8),parameter :: f_depcnt_dust = (2+m_depcnt_dst)*(1-m_depcnt_dst)**2/4._r8

   real(r8),parameter :: m_imm_bc = COS(theta_imm_bc*pi/180._r8)
   real(r8),parameter :: f_imm_bc = (2+m_imm_bc)*(1-m_imm_bc)**2/4._r8

   real(r8),parameter :: m_imm_dust = COS(theta_imm_dust*pi/180._r8)
   real(r8),parameter :: f_imm_dust = (2+m_imm_dust)*(1-m_imm_dust)**2/4._r8

   real(r8) :: f_dep, f_cnt, f_imm
   real(r8) :: dga_dep, dga_imm
   real(r8) :: limfac
   real(r8) :: frzimm, frzcnt, frzdep
   real(r8), pointer :: frzimm_ptr, frzcnt_ptr, frzdep_ptr

   logical :: pdf_imm

   integer :: ispc

   real(r8) :: ktherm(ntypes), kcoll(ntypes)

   real(r8), parameter :: Ktherm_bc = 4.2_r8   ! black carbon thermal conductivity [J/(m s K)]
   real(r8), parameter :: Ktherm_dst = 0.72_r8 ! clay thermal conductivity [J/(m s K)]

   !------------------------------------------------------------------------------------------------

   errstring = ' '

   nullify(frzimm_ptr)
   nullify(frzcnt_ptr)
   nullify(frzdep_ptr)

   frzbcimm = 0._r8
   frzbccnt= 0._r8
   frzbcdep = 0._r8
   frzduimm = 0._r8
   frzducnt= 0._r8
   frzdudep = 0._r8

   ! get saturation vapor pressure
   eswtr = svp_water(t)  ! 0 for liquid

   tc = T - tmelt
   rhoice = 916.7_r8-0.175_r8*tc-5.e-4_r8*tc**2
   vwice = mwh2o*amu/rhoice
   sigma_iw = (28.5_r8+0.25_r8*tc)*1E-3_r8
   sigma_iv = (76.1_r8-0.155_r8*tc + 28.5_r8+0.25_r8*tc)*1E-3_r8

   ! critical germ size
   rgimm = 2*vwice*sigma_iw/(boltz*T*LOG(supersatice))

   ! critical germ size
   ! assume 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)
   rgdep=2*vwice*sigma_iv/(boltz*T*LOG(rhwincloud*supersatice))

   ! homogeneous energy of germ formation
   dg0dep = 4*pi/3._r8*sigma_iv*rgdep**2

   ! prefactor
   ! attention: division of small numbers
   Adep = (rhwincloud*eswtr)**2*(vwice/(mwh2o*amu))/(boltz*T*nus)*SQRT(sigma_iv/(boltz*T))

   ! homogeneous energy of germ formation
   dg0cnt = 4*pi/3._r8*sigma_iv*rgimm**2

   ! prefactor
   ! attention: division of small numbers
   Acnt = rhwincloud*eswtr*4*pi/(nus*SQRT(2*pi*mwh2o*amu*boltz*T))

   do ispc = 1, ntypes

      select case (trim(types(ispc)))
      case ('black-c')
         ktherm(ispc) = ktherm_bc
      case ('dust')
         ktherm(ispc) = ktherm_dst
      case default
         errstring = 'hetfrz_classnuc_calc ERROR: unrecognized aerosol type: '//trim(types(ispc))
         return
      end select
   end do

   call collkernel(T, p, eswtr, rhwincloud, r3lx, hetraer, Ktherm, Kcoll)

   do ispc = 1, ntypes

      select case (trim(types(ispc)))
      case ('black-c')
         f_dep = f_depcnt_bc
         f_cnt = f_depcnt_bc
         f_imm = f_imm_bc
         dga_dep = dga_dep_bc
         dga_imm = dga_imm_bc
         pdf_imm = .false.
         limfac = bc_limfac
         frzimm_ptr => frzbcimm
         frzcnt_ptr => frzbccnt
         frzdep_ptr => frzbcdep
      case ('dust')
         f_dep = f_depcnt_dust
         f_cnt = f_depcnt_dust
         f_imm = f_imm_dust
         dga_dep = dga_dep_dust
         dga_imm = dga_imm_dust
         pdf_imm = .true.
         limfac = dust_limfac
         frzimm_ptr => frzduimm
         frzcnt_ptr => frzducnt
         frzdep_ptr => frzdudep
      case default
         errstring = 'hetfrz_classnuc_calc ERROR: unrecognized aerosol type: '//trim(types(ispc))
         return
      end select

      call hetfrz_classnuc_calc_rates( f_dep, f_cnt, f_imm, dga_dep, dga_imm, pdf_imm, limfac, &
           kcoll(ispc), hetraer(ispc), icnlx, r3lx, T, supersatice, sigma_iw, &
           rgimm, rgdep, dg0dep, Adep, dg0cnt, Acnt, vwice, deltat, &
           fn(ispc), wact_factor(ispc), dstcoat(ispc), &
           total_aer_num(ispc), total_interstitial_aer_num(ispc), total_cloudborne_aer_num(ispc), uncoated_aer_num(ispc), &
           frzimm, frzcnt, frzdep, errstring )

      ! accumulate dust and bc frz rates
      frzimm_ptr = frzimm_ptr + frzimm
      frzcnt_ptr = frzcnt_ptr + frzcnt
      frzdep_ptr = frzdep_ptr + frzdep

   end do

 end subroutine  hetfrz_classnuc_calc

 subroutine  hetfrz_classnuc_calc_rates( f_dep, f_cnt, f_imm, dga_dep, dga_imm, pdf_imm, limfac, &
      kcoll, mradius, icnlx, r3lx, T, supersatice, sigma_iw, &
      rgimm, rgdep, dg0dep, Adep, dg0cnt, Acnt, vwice, deltat, &
      fn, wact_factor, dstcoat, &
      total_aer_num, total_interstitial_aer_num, total_cloudborne_aer_num, uncoated_aer_num, &
      frzimm, frzcnt, frzdep, errstring )

   ! input
   real(r8), intent(in) :: f_dep                      ! deposition form factor
   real(r8), intent(in) :: f_cnt                      ! contact form factor
   real(r8), intent(in) :: f_imm                      ! immersion form factor
   real(r8), intent(in) :: dga_dep                    ! deposition activation energy [J]
   real(r8), intent(in) :: dga_imm                    ! immersion activation energy [J]
   logical,  intent(in) :: pdf_imm                    ! PDF theta model switch (TRUE for dust)
   real(r8), intent(in) :: limfac                     ! Limit to 1% of available potential IN (for BC), no limit for dust
   real(r8), intent(in) :: kcoll                      ! collision kernel [cm3 s-1]
   real(r8), intent(in) :: mradius                    ! mass mean radius [m]
   real(r8), intent(in) :: icnlx                      ! in-cloud droplet concentration [cm-3]
   real(r8), intent(in) :: r3lx                       ! volume mean drop radius [m]
   real(r8), intent(in) :: T                          ! temperature [K]
   real(r8), intent(in) :: supersatice                ! supersaturation ratio wrt ice at 100%rh over water [ ]
   real(r8), intent(in) :: sigma_iw                   ! [J/m2]
   real(r8), intent(in) :: rgimm                      ! critical germ size
   real(r8), intent(in) :: rgdep                      ! critical germ size
   real(r8), intent(in) :: dg0dep                     ! homogeneous energy of germ formation
   real(r8), intent(in) :: Adep                       ! deposition nucleation prefactor
   real(r8), intent(in) :: dg0cnt                     ! homogeneous energy of germ formation
   real(r8), intent(in) :: Acnt                       ! contact nucleation prefactor

   real(r8), intent(in) :: vwice
   real(r8), intent(in) :: deltat                     ! timestep [s]
   real(r8), intent(in) :: fn                         ! fraction activated [ ] for cloud borne aerosol number
   real(r8), intent(in) :: wact_factor                ! water activity factor -- density*(1.-(OC+BC)/(OC+BC+SO4)) [mug m-3]
   real(r8), intent(in) :: dstcoat                    ! coated fraction
   real(r8), intent(in) :: total_aer_num              ! total bc and dust number concentration(interstitial+cloudborne) [#/cm^3]
   real(r8), intent(in) :: total_interstitial_aer_num ! total bc and dust concentration(interstitial)
   real(r8), intent(in) :: total_cloudborne_aer_num   ! total bc and dust concentration(cloudborne)
   real(r8), intent(in) :: uncoated_aer_num           ! uncoated bc and dust number concentration(interstitial)
   ! output
   real(r8), intent(out) :: frzimm ! het. frz by immersion nucleation [cm-3 s-1]
   real(r8), intent(out) :: frzcnt ! het. frz by contact nucleation [cm-3 s-1]
   real(r8), intent(out) :: frzdep ! het. frz by deposition nucleation [cm-3 s-1]

   character(len=*), intent(out) :: errstring

   ! local vars
   real(r8) :: aw                          ! water activity [ ]
   real(r8) :: molal                       ! molality [moles/kg]

   real(r8) :: Aimm
   real(r8) :: Jdep
   real(r8) :: Jimm
   real(r8) :: Jcnt
   real(r8) :: dg0imm
   real(r8) :: rgimm_aer
   real(r8) :: sum_imm
   real(r8) :: dim_Jimm(pdf_n_theta)

   logical :: do_frz
   logical :: do_imm

   integer :: i

   !*****************************************************************************
   !                take water activity into account
   !*****************************************************************************
   !   solute effect
   aw = 1._r8
   molal = 0._r8

   ! The heterogeneous ice freezing temperatures of all IN generally decrease with
   ! increasing total solute mole fraction. Therefore, the large solution concentration
   ! will cause the freezing point depression and the ice freezing temperatures of all
   ! IN will get close to the homogeneous ice freezing temperatures. Since we take into
   ! account water activity for three heterogeneous freezing modes(immersion, deposition,
   ! and contact), we utilize interstitial aerosols(not cloudborne aerosols) to calculate
   ! water activity.
   ! If the index of IN is 0, it means three freezing modes of this aerosol are depressed.

   !calculate molality
   if ( total_interstitial_aer_num > 0._r8 ) then
      molal = (1.e-6_r8*wact_factor/(mwso4*total_interstitial_aer_num*1.e6_r8))/ &
           (4*pi/3*rhoh2o*(MAX(r3lx,4.e-6_r8))**3)
      aw = 1._r8/(1._r8+2.9244948e-2_r8*molal+2.3141243e-3_r8*molal**2+7.8184854e-7_r8*molal**3)
   end if

   !*****************************************************************************
   !                immersion freezing begin
   !*****************************************************************************

   frzimm = 0._r8
   frzcnt = 0._r8
   frzdep = 0._r8

   ! take solute effect into account
   rgimm_aer = rgimm

   ! if aw*Si<=1, the freezing point depression is strong enough to prevent freezing

   do_frz = aw*supersatice > 1._r8
   if (do_frz) then
      rgimm_aer = 2*vwice*sigma_iw/(boltz*T*LOG(aw*supersatice))
   else
      return
   endif

   do_imm = T <= 263.15_r8 ! temperature threshold for immersion freezing (-10 C)

   if (do_imm) then
      ! homogeneous energy of germ formation
      dg0imm = 4*pi/3._r8*sigma_iw*rgimm_aer**2

      ! prefactor
      Aimm = n1*((vwice*rhplanck)/(rgimm_aer**3)*SQRT(3._r8/pi*boltz*T*dg0imm))

      ! nucleation rate per particle

      if (pdf_imm) then
         dim_Jimm(:) = 0._r8
         do i = i1,i2
            ! 1/sqrt(f)
            dim_Jimm(i) = Aimm*mradius**2/SQRT(dim_f_imm(i))*EXP((-dga_imm-dim_f_imm(i)*dg0imm)/(boltz*T))
            dim_Jimm(i) = max(dim_Jimm(i), 0._r8)
         end do

         sum_imm  = 0._r8
         do i = i1,i2-1
            sum_imm = sum_imm + 0.5_r8*((pdf_imm_theta(i  )*exp(-dim_Jimm(i  )*deltat)+ &
                 pdf_imm_theta(i+1)*exp(-dim_Jimm(i+1)*deltat)))*pdf_d_theta
         end do
         if (sum_imm > 0.99_r8) then
            sum_imm = 1.0_r8
         end if
      else
         Jimm = Aimm*mradius**2/SQRT(f_imm)*EXP(( -dga_imm - f_imm*dg0imm )/(boltz*T))
         sum_imm = exp(-Jimm*deltat)
      end if
   end if

   !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   Deposition nucleation
   !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! nucleation rate per particle
   if (rgdep > 0) then
      Jdep = Adep*mradius**2/SQRT(f_dep)*EXP((-dga_dep-f_dep*dg0dep)/(boltz*T))
   else
      Jdep = 0._r8
   end if

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! contact nucleation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! nucleation rate per particle
   Jcnt = Acnt*mradius**2*EXP((-dga_dep-f_cnt*dg0cnt)/(boltz*T))*Kcoll*icnlx

   ! Limit to 1% of available potential IN (for BC), no limit for dust
   if (tot_in) then
      if (do_imm) then
         frzimm = MIN(limfac*fn*total_aer_num/deltat, fn*total_aer_num/deltat*(1._r8-sum_imm))
      end if
      frzdep = MIN(limfac*(1._r8-fn)*(1._r8-dstcoat)*total_aer_num/deltat, &
                   (1._r8-fn)*(1._r8-dstcoat)*total_aer_num/deltat*(1._r8-exp(-Jdep*deltat)))
      frzcnt = MIN(limfac*(1._r8-fn)*(1._r8-dstcoat)*total_aer_num/deltat, &
                   (1._r8-fn)*(1._r8-dstcoat)*total_aer_num/deltat*(1._r8-exp(-Jcnt*deltat)))
   else
      if (do_imm) then
         frzimm = MIN(limfac*total_cloudborne_aer_num /deltat, total_cloudborne_aer_num/deltat*(1._r8-sum_imm))
      end if
      frzdep = MIN(limfac*uncoated_aer_num/deltat, uncoated_aer_num/deltat*(1._r8-exp(-Jdep*deltat)))
      frzcnt = MIN(limfac*uncoated_aer_num/deltat, uncoated_aer_num/deltat*(1._r8-exp(-Jcnt*deltat)))
   end if

   if (frzcnt <= -1._r8) then
      write(iulog,*) 'hetfrz_classnuc_calc: frzcnt, Jcnt, Kcoll: ', frzcnt, Jcnt, Kcoll
      errstring = 'ERROR in hetfrz_classnuc_calc::frzducnt'
      return
   end if

end subroutine  hetfrz_classnuc_calc_rates

!===================================================================================================

!-----------------------------------------------------------------------
!
! Purpose: calculate collision kernels as a function of environmental parameters and aerosol/droplet sizes
!
! Author: Corinna Hoose, UiO, October 2009
!
! Modifications: Yong Wang and Xiaohong Liu, UWyo, 12/2012
!
! "Seinfeld & Pandis" referenced in several places in this routine is:
! Atmospheric Chemistry and Physics: From Air Pollution to Climate Change, 3rd Edition
! John H. Seinfeld, Spyros N. Pandis  ISBN: 978-1-118-94740-1
!-----------------------------------------------------------------------

subroutine collkernel( temp, pres, eswtr, rhwincloud, r3lx,  rad, Ktherm, Kcoll )

   real(r8), intent(in) :: temp       ! temperature [K]
   real(r8), intent(in) :: pres       ! pressure [Pa]
   real(r8), intent(in) :: eswtr      ! saturation vapor pressure of water [Pa]
   real(r8), intent(in) :: r3lx       ! volume mean drop radius [m]
   real(r8), intent(in) :: rhwincloud ! in-cloud relative humidity over water [ ]
   real(r8), intent(in) :: rad(:)     ! aerosol radius [m]
   real(r8), intent(in) :: Ktherm(:)  ! thermal conductivity of aerosol [J/(m s K)]
   real(r8), intent(out) :: Kcoll(:)  ! collision kernel [cm3 s-1]

   ! local variables
   real(r8) :: a, b, c, a_f, b_f, c_f, f
   real(r8) :: tc          ! temperature [deg C]
   real(r8) :: rho_air     ! air density [kg m-3]
   real(r8) :: viscos_air  ! dynamic viscosity of air [kg m-1 s-1]
   real(r8) :: Ktherm_air  ! thermal conductivity of air [J/(m s K)]
   real(r8) :: lambda      ! mean free path [m]
   real(r8) :: Kn          ! Knudsen number [ ]
   real(r8) :: Re          ! Reynolds number [ ]
   real(r8) :: Pr          ! Prandtl number [ ]
   real(r8) :: Sc          ! Schmidt number [ ]
   real(r8) :: vterm       ! terminal velocity [m s-1]
   real(r8) :: Dvap        ! water vapor diffusivity [m2 s-1]
   real(r8) :: Daer        ! aerosol diffusivity [m2 s-1]
   real(r8) :: latvap      ! latent heat of vaporization [J kg-1]
   real(r8) :: G           ! thermodynamic function in Cotton et al. [kg m-1 s-1]
   real(r8) :: f_t         ! factor by Waldmann & Schmidt [ ]
   real(r8) :: Q_heat      ! heat flux [J m-2 s-1]
   real(r8) :: Tdiff_cotton ! temperature difference between droplet and environment [K]
   real(r8) :: K_brownian,K_thermo_cotton,K_diffusio_cotton   ! collision kernels [m3 s-1]

   integer :: ntot, idx

   !------------------------------------------------------------------------------------------------

   ntot = size(ktherm)

   Kcoll(:) = 0._r8

   tc = temp - tmelt

   ! air viscosity for tc<0, from depvel_part.F90
   viscos_air = (1.718_r8+0.0049_r8*tc-1.2e-5_r8*tc*tc)*1.e-5_r8
   ! air density
   rho_air = pres/(rair*temp)
   ! mean free path: Seinfeld & Pandis 8.6 (Book: ISBN: 978-1-118-94740-1)
   lambda = 2*viscos_air/(pres*SQRT(8/(pi*rair*temp)))
   ! latent heat of vaporization, varies with T
   latvap = 1000*(-0.0000614342_r8*tc**3 + 0.00158927_r8*tc**2 - 2.36418_r8*tc + 2500.79_r8)
   ! droplet terminal velocity after Chen & Liu, QJRMS 2004 (https://doi-org.cuucar.idm.oclc.org/10.1256/qj.03.41)
   a = 8.8462e2_r8
   b = 9.7593e7_r8
   c = -3.4249e-11_r8
   a_f = 3.1250e-1_r8
   b_f = 1.0552e-3_r8
   c_f = -2.4023_r8
   f = EXP(EXP(a_f + b_f*(LOG(r3lx))**3 + c_f*rho_air**1.5_r8))
   vterm = (a+ (b + c*r3lx)*r3lx)*r3lx*f

   ! Reynolds number
   Re = 2*vterm*r3lx*rho_air/viscos_air
   ! thermal conductivity of air: Seinfeld & Pandis eq. 15.75 (Book: ISBN: 978-1-118-94740-1)
   Ktherm_air = 1.e-3_r8*(4.39_r8+0.071_r8*temp)  !J/(m s K)
   ! Prandtl number
   Pr = viscos_air*cpair/Ktherm_air
   ! water vapor diffusivity: Pruppacher & Klett 13-3 (https://link.springer.com/book/10.1007/978-0-306-48100-0)
   Dvap = 0.211e-4_r8*(temp/tmelt)*(pstd/pres)
   ! G-factor = rhoh2o*Xi in Rogers & Yau, p. 104
   G = rhoh2o/((latvap/(rh2o*temp) - 1)*latvap*rhoh2o/(Ktherm_air*temp) &
       + rhoh2o*rh2o*temp/(Dvap*eswtr))

   do idx = 1,ntot
      if (rad(idx)>0._r8) then
         ! Knudsen number (Seinfeld & Pandis 8.1) (Book: ISBN: 978-1-118-94740-1)
         Kn = lambda/rad(idx)
         ! aerosol diffusivity
         Daer = boltz*temp*(1 + Kn)/(6*pi*rad(idx)*viscos_air)

         ! Schmidt number
         Sc = viscos_air/(Daer*rho_air)

         ! Young (1974) first equ. on page 771 (doi: 10.1175/1520-0469(1974)031<0768:TROCNI>2.0.CO;2)
         K_brownian = 4*pi*r3lx*Daer*(1 + 0.3_r8*Re**0.5_r8*Sc**0.33_r8)

         ! thermal conductivities from Seinfeld & Pandis, Table 8.6 (Book: ISBN: 978-1-118-94740-1)
         ! form factor
         f_t = 0.4_r8*(1._r8 + 1.45_r8*Kn + 0.4_r8*Kn*EXP(-1._r8/Kn))      &
              *(Ktherm_air + 2.5_r8*Kn*Ktherm(idx))                      &
              /((1._r8 + 3._r8*Kn)*(2._r8*Ktherm_air + 5._r8*Kn*Ktherm(idx)+Ktherm(idx)))

         ! calculate T-Tc as in Cotton et al.
         Tdiff_cotton = -G*(rhwincloud - 1._r8)*latvap/Ktherm_air
         Q_heat = Ktherm_air/r3lx*(1._r8 + 0.3_r8*Re**0.5_r8*Pr**0.33_r8)*Tdiff_cotton
         K_thermo_cotton = 4._r8*pi*r3lx*r3lx*f_t*Q_heat/pres
         K_diffusio_cotton = -(1._r8/f_t)*(rh2o*temp/latvap)*K_thermo_cotton
         Kcoll(idx) = 1.e6_r8*(K_brownian + K_thermo_cotton + K_diffusio_cotton)  ! convert m3/s -> cm3/s

         ! set K to 0 if negative
         if (Kcoll(idx) < 0._r8) Kcoll(idx) = 0._r8
      else
         Kcoll(idx) = 0._r8
      endif
   end do

end subroutine collkernel

!===================================================================================================

subroutine hetfrz_classnuc_init_pdftheta()

   ! Local variables:
   real(r8) :: theta_min, theta_max
   real(r8) :: x1_imm, x2_imm
   real(r8) :: norm_theta_imm
   real(r8) :: imm_dust_mean_theta
   real(r8) :: imm_dust_var_theta
   integer  :: i
   real(r8) :: m
   real(r8) :: temp
   !----------------------------------------------------------------------------

   theta_min           = pi/180._r8
   theta_max           = 179._r8/180._r8*pi
   imm_dust_mean_theta = 46.0_r8/180.0_r8*pi
   imm_dust_var_theta  = 0.01_r8

   pdf_d_theta = (179._r8-1._r8)/180._r8*pi/(pdf_n_theta-1)

   x1_imm = (LOG(theta_min) - LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
   x2_imm = (LOG(theta_max) - LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
   norm_theta_imm = (ERF(x2_imm) - ERF(x1_imm))*0.5_r8
   dim_theta      = 0.0_r8
   pdf_imm_theta  = 0.0_r8
   do i = i1, i2
      dim_theta(i)     = 1._r8/180._r8*pi + (i-1)*pdf_d_theta
      pdf_imm_theta(i) = exp(-((LOG(dim_theta(i)) - LOG(imm_dust_mean_theta))**2._r8) / &
                              (2._r8*imm_dust_var_theta**2._r8) ) /                     &
                         (dim_theta(i)*imm_dust_var_theta*SQRT(2*pi))/norm_theta_imm
   end do

   do i = i1, i2
     m = cos(dim_theta(i))
     temp = (2+m)*(1-m)**2/4._r8
     dim_f_imm(i) = temp
   end do

end subroutine hetfrz_classnuc_init_pdftheta

!===================================================================================================

end module hetfrz_classnuc
