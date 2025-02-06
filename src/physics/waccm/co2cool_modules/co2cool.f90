!======================================================================!
! MODULE CO2COOL ---                                                   !
!======================================================================!
! This parameterization calculates the CO2 15 micron cooling under     !
! non-LTE.                                                             !
! It is an improved and extended parameterization of the CO2 15 μm     !
! cooling in the middle/upper atmosphere originally published by       !
! Fomichev et al. (1998). The major improvement is its extension to    !
! cope with CO2 vmrs from ∼0.5 to over 10 times the CO2 pre-industrial !
! value of 284 ppmv (i.e., 150 ppmv to 3000 ppmv). Furthermore, it     !
! incorporates a more contemporary CO2 line list and the collisional   !
! rates that affect the CO2 15 non-LTE cooling rates.                  !
! It has been specifically designed for being used in GCM climate      !
! models and it is very fast.                                          !
! Authors:                                                             !
! Manuel López-Puertas, Federico Fabiano, Victor Fomichev, Bernd Funke,!
! and Daniel R. Marsh                                                  !
! Please reference the parameterization (to be updated) by:            !
! López-Puertas, M., Fabiano, F., Fomichev, V., Funke, B., and Marsh,  !
! D. R.: An improved and extended parameterization of the CO2 15 µm    !
! cooling in the middle/upper atmosphere, EGUsphere [preprint],        !
! https://doi.org/10.5194/egusphere-2023-2424, 2023.                   !
!----------------------------------------------------------------------!
! STATUS: BF 20-Feb-2024                      CREATED: BF 20-Feb-2024  !
!----------------------------------------------------------------------!
! COPYRIGHT: (C) 2024-2024 Instituto de Astrofisica de Andalucia (IAA) !
!----------------------------------------------------------------------!
! LICENSE: GNU Lesser General Public License (LGPL) version 2.1        !
!          (see file LICENSE_LGPLv2.1)                                 !
!----------------------------------------------------------------------!
! CONTAINS:                                                            !
! public:                                                              !
!             NEW_PARAM                                                !
!----------------------------------------------------------------------!
! MAINTENANCE HISTORY:                                                 !
! bf  24-02-24:   adapted from python code                             !
! Notes added by MLP on 17 march 2024                                  !
!======================================================================!

module CO2COOL
  use precision

  implicit NONE
  private
  public :: CO2_NLTE_COOL

CONTAINS

!======================================================================!
!  CO2_NLTE_COOL --- subroutine for external call of parameterization  !
!----------------------------------------------------------------------!
! STATUS: bf 20-Feb-2024                      CREATED: bf 20-Feb-2024  !
!----------------------------------------------------------------------!
! INPUT:                                                               !
!       temp   --- vector input temperatures                           !
!       co2vmr --- vector input co2 vmr                                !
!       ovmr   --- vector input o vmr                                  !
!       o2vmr  --- vector input o2 vmr                                 !
!       n2vmr  --- vector input n2 vmr                                 !
!       surf_temp  --- real surface temperature                        !
!       lev0       --- integer maximum pressure level of input grid    !
!                      to be calculated                                !
!----------------------------------------------------------------------!
! OUTPUT:                                                              !
!       hr    --- computed heating rate on input grid (K day-1)        !
!----------------------------------------------------------------------!
! NOTE: This routine computes the heating rates. It set up the internal!
! reference grid in x, determines the ranges of calculation, calculated!
! the cooling in its internal grid and interpolatess them back to the  !
! pressure grid povided by the user.                                   !
!
!----------------------------------------------------------------------!
! MAINTENANCE HISTORY:                                                 !
! bf  24-02-24:   adapted from python code                             !
! MLP 17-03-24:   Added some notes                                     !
!======================================================================!

subroutine CO2_NLTE_COOL (temp, pres, co2vmr, ovmr, o2vmr, n2vmr, lev0, &
                      surf_temp, hr)

use constants
use varsub, only : error

   implicit NONE
   real(dp), dimension(:), intent(IN)    :: temp, pres, co2vmr, ovmr, &
                                            o2vmr, n2vmr
   real(dp), dimension(:), intent(INOUT) :: hr
   integer,               intent(IN)     :: lev0
   real(dp),              intent(IN)     :: surf_temp

   real(dp), allocatable :: xatm(:), hrout(:), xref(:)
   real(dp) :: c_int
   integer :: n_inlev, i, max_lev, min_lev, grd_lev, istat, ilev
   character(*), parameter :: routine = 'NEW_PARAM'

   n_inlev = size(pres)
   if (size(temp) /= n_inlev) &
      call error(routine, 'inconsistent temperature grid!')
   if (size(co2vmr) /= n_inlev) &
      call error(routine, 'inconsistent CO2 vmr grid!')
   if (size(ovmr) /= n_inlev) &
      call error(routine, 'inconsistent O vmr grid!')
   if (size(o2vmr) /= n_inlev) &
      call error(routine, 'inconsistent O2 vmr grid!')
   if (size(n2vmr) /= n_inlev) &
      call error(routine, 'inconsistent N2 vmr grid!')
   if (size(hr) /= n_inlev) &
      call error(routine, 'inconsistent hr grid!')
   if (lev0 > n_inlev) &
      call error(routine, 'lev0 out of range')

! Convert pressure to log scale
!------------------------------
   allocate(xatm(n_inlev),stat=istat)
   if (istat/=0) call error(routine ,'Allocation of xatm failed')
   xatm = log(p0/pres)

! Set reference grid
!-------------------
   max_lev = ceiling((maxval(xatm)-0.125_dp) / 0.25_dp) + 1
   allocate(xref(max_lev),stat=istat)
   if (istat/=0) call error(routine ,'Allocation of xref failed')
   xref    = (/ (i*0.25_db, i=1,max_lev) /) - 0.125_dp

! Determine ranges
!-----------------
   if (xatm(1) > xatm(n_inlev)) then
      call HUNT (xref,max_lev,xatm(n_inlev),grd_lev)
      if (grd_lev == 0) grd_lev=1
      call HUNT (xref,max_lev,xatm(lev0),min_lev)
      if (min_lev == 0) min_lev=1
   else
      call HUNT (xref,max_lev,xatm(1),grd_lev)
      if (grd_lev == 0) grd_lev=1
      call HUNT (xref,max_lev,xatm(lev0),min_lev)
      if (min_lev == 0) min_lev=1
   endif

! Calculate cooling
!------------------
   allocate(hrout(max_lev),stat=istat)
   if (istat/=0) call error(routine ,'Allocation of hrout failed')

   hrout = CALC_COOL(xatm, temp, co2vmr, ovmr, o2vmr, n2vmr, surf_temp, &
           xref, n_inlev, max_lev, min_lev, grd_lev)

! Interpolate back to input grid at pressures lower than press(lev0)
!-------------------------------------------------------------------
   do ilev = 1,n_inlev
      if (pres(ilev) <= pres(lev0)) then
         call HUNT (xref,max_lev,xatm(ilev),i)
         if (i == 0) i = 1
         c_int = (xatm(ilev) - xref(i)) / (xref(i+1)-xref(i))
         hr(ilev) = hrout(i) + (hrout(i+1) - hrout(i)) * c_int
      endif
   end do

   deallocate(xatm, hrout, xref, stat=istat)
   if (istat/=0) call error(routine ,'Deallocation of xatm, hrout, xref failed')

end subroutine CO2_NLTE_COOL

!======================================================================!
! CALC_COOL --- function                                               !
! This function performs the calculation of the heating rate           !
! Interpolates the user input reference atmosphere into its internal   !
! x- grid and perform all calculations in the LTE and different non-LTE!
! regions.                                                             !
!                                                                      !
!----------------------------------------------------------------------!
! STATUS: bf 20-Feb-2024                      CREATED: bf 20-Feb-2024  !
!----------------------------------------------------------------------!
! INPUT:                                                               !
!       xatm   --- vector input log pressure heights                   !
!       temp   --- vector input temperatures                           !
!       co2vmr --- vector input co2 vmr                                !
!       ovmr   --- vector input o vmr                                  !
!       o2vmr  --- vector input o2 vmr                                 !
!       n2vmr  --- vector input n2 vmr                                 !
!       surf_temp  --- surface temperature                             !
!       n_inlev    --- integer input grid dimension                    !
!       n_lev      --- integer internal max level to be computed       !
!       min_lev    --- integer internal min level to be computed       !
!       grd_lev    --- integer internal surface level                  !
!----------------------------------------------------------------------!
! RESULT:                                                              !
!       hr    --- computed heating rate (K day-1) on internal grid     !
!----------------------------------------------------------------------!
! NOTE:                                                                !
!----------------------------------------------------------------------!
! MAINTENANCE HISTORY:                                                 !
! bf  24-02-24:   adapted from python code                             !
! MLP 17-March-2024:  added some notes                                 !
! MLP 18-March-2024:  Tabs removed                                     !
!======================================================================!

function CALC_COOL (xatm, temp, co2vmr, ovmr, o2vmr, n2vmr, surf_temp, &
                    xref, n_inlev,n_lev, min_lev, grd_lev) result (hr)

use constants
use coedat

   implicit NONE
   real(dp), dimension(n_inlev), intent(IN) :: temp, co2vmr, ovmr, &
                                               o2vmr, n2vmr, xatm
   real(db),dimension(n_lev), intent(IN)   :: xref
   real(dp), intent(IN)    :: surf_temp
   integer,  intent(IN)    :: n_inlev, n_lev, min_lev, grd_lev

   real(db),dimension(n_alpha)           :: alphaf
   real(db),dimension(n_Lesc)            :: Lescf
   real(db),dimension(n_lev_cm)          :: asurff, bsurff
   real(db),dimension(n_lev_cm,n_lev_cm) :: acoeff, bcoeff
   real(db),dimension(n_lev) :: temp_x, pres_x, co2vmr_x, ovmr_x, o2vmr_x, &
                                n2vmr_x, phi_fun, hr, co2col, n_dens, &
                                molmass, lamb, cp, fac, dj, hb, Lesc_x, &
                                alpha_ok, eps_gn
   real(db),dimension(n_co2prof) :: co2p
   real(dp) :: co2m, c_int, phi_fun_g, kost, h, t13, tsq, zo, &
               zo2, zn2, eps125, Djj, Djjm1, Fj, Fjm1, Phi_165
   integer  :: ilev, ico2, i, mxc_lev, mic_lev, n_lev_rui, n_lev_csi


! Interpolate atmospheric input onto reference grid
!--------------------------------------------------
   do ilev = 1, n_lev
      call HUNT (xatm,n_inlev,xref(ilev),i)
      if (i == 0) i = 1
      if (i == n_inlev) i = n_inlev-1
      c_int = (xref(ilev) - xatm(i)) / (xatm(i+1)-xatm(i))
      temp_x(ilev)   = temp(i) + (temp(i+1) - temp(i)) * c_int
      co2vmr_x(ilev) = co2vmr(i) + (co2vmr(i+1) - co2vmr(i)) * c_int
      if (ilev >=n_lev_rf-1) then
         o2vmr_x(ilev)  = o2vmr(i) + (o2vmr(i+1) - o2vmr(i)) * c_int
         ovmr_x(ilev)   = ovmr(i) * (ovmr(i+1) / ovmr(i))**c_int
         n2vmr_x(ilev)  = n2vmr(i) + (n2vmr(i+1) - n2vmr(i)) * c_int
         pres_x(ilev)   = p0 * exp(-xref(ilev))
      endif
   end do

! Set phi functions
!------------------
   phi_fun   = exp(-Edkbc / temp_x)
   phi_fun_g = exp(-Edkbc / surf_temp)

! Interpolate CM coefficients onto actual CO2
!--------------------------------------------
   mic_lev = min(min_lev, n_lev_rf-1)
   mxc_lev = min(n_lev, n_lev_cm)

   do ilev = mic_lev, mxc_lev
      co2p = co2profs(:,ilev)
      call HUNT (co2p,n_co2prof,co2vmr_x(ilev),ico2)
      if (ico2 == 0) ico2 = 1
      if (ico2 == n_co2prof) ico2 = n_co2prof-1
      c_int = (co2vmr_x(ilev) - co2p(ico2)) / (co2p(ico2+1)-co2p(ico2))
      asurff(ilev)   = asurf(ico2,ilev) * (asurf(ico2+1,ilev)/ &
                       asurf(ico2,ilev))**c_int
      bsurff(ilev)   = bsurf(ico2,ilev) * (bsurf(ico2+1,ilev)/ &
                       bsurf(ico2,ilev))**c_int
      acoeff(:,ilev) = acoef(ico2,:,ilev) * (acoef(ico2+1,:,ilev)/ &
                       acoef(ico2,:,ilev))**c_int
      bcoeff(:,ilev) = bcoef(ico2,:,ilev) * (bcoef(ico2+1,:,ilev)/ &
                       bcoef(ico2,:,ilev))**c_int
   end do

! Calculate the cooling in the LTE region with Curtis matrix
!----------------------------------------
   do ilev = mic_lev, mxc_lev
      hr(ilev) = sum((acoeff(grd_lev:mxc_lev,ilev) + &
                 bcoeff(grd_lev:mxc_lev,ilev) * &
                 phi_fun(ilev)) * phi_fun(grd_lev:mxc_lev)) ! only consider contributions from >= grd_lev
      hr(ilev) = hr(ilev) + (asurff(ilev) + &
                 bsurff(ilev) * phi_fun(ilev)) * phi_fun_g
   end do

   if (n_lev < n_lev_rf) then  ! no non-LTE region touched
      return
   endif

! Adjust level boundaries
!------------------------
   n_lev_rui = min(n_lev_ru-1, n_lev)+1
   n_lev_csi = min(n_lev_cs, n_lev)

! Interpolate Lesc (the escape function) onto actual CO2
!-----------------------------------------
   co2m = sum(co2vmr_x(1:n_meanco2)) / (n_meanco2*one)
   call hunt (co2mean, n_co2prof, co2m, ico2)
   if (ico2 == 0) ico2 = 1
   if (ico2 == n_co2prof) ico2 = n_co2prof-1
   c_int = (co2m - co2mean(ico2)) / (co2mean(ico2+1)-co2mean(ico2))
   Lescf  = Lesc(ico2,:) + (Lesc(ico2+1,:) - Lesc(ico2,:)) * c_int

! Interpolate the alpha coefficient (correction factor) onto actual CO2
!-----------------------------------------
   do ilev = n_lev_rf,n_lev_rui
      co2p = co2profs(:,ilev)
      call HUNT (co2p, n_co2prof, co2vmr_x(ilev), ico2)
      if (ico2 == 0) ico2 = 1
      if (ico2 == n_co2prof) ico2 = n_co2prof-1
      c_int = (co2vmr_x(ilev) - co2p(ico2)) / (co2p(ico2+1)-co2p(ico2))
      alphaf(ilev-n_lev_rf+1)  = alpha(ico2,ilev-n_lev_rf+1) + &
                                 (alpha(ico2+1,ilev-n_lev_rf+1) - &
                                  alpha(ico2,ilev-n_lev_rf+1)) * c_int
   end do

! Calculate density and molecular mass, required to express the heating in K/day
!------------------------------
   do ilev = n_lev_rf-1,n_lev
      n_dens(ilev)  = pres_x(ilev) / (kb * temp_x(ilev))
      molmass(ilev) = (n2vmr_x(ilev)*mn2 + o2vmr_x(ilev)*mo2 + ovmr_x(ilev)*mo)/ &
                      (n2vmr_x(ilev) + o2vmr_x(ilev) + ovmr_x(ilev))
   end do

! Calculate CO2 column
!---------------------
   co2col(n_lev_csi) = co2col_ovh(n_lev_csi) * &
                       co2vmr_x(n_lev_csi) / co2profs(i_co2ref,n_lev_csi)
   hb(n_lev_csi) = co2vmr_x(n_lev_csi)/molmass(n_lev_csi) ! need hb at lev n_lev_csi
   do ilev=n_lev_csi-1,n_lev_rf-1,-1
      kost    = -Nave1 / g_grav(ilev)
      hb(ilev) = co2vmr_x(ilev)/molmass(ilev)
      co2col(ilev) = 0.5_db*(hb(ilev+1) + hb(ilev)) * (pres_x(ilev+1) - &
                    pres_x(ilev)) * kost + co2col(ilev+1)
   end do

! Interpolate Lescape in log-log
!-------------------------------
   do ilev = n_lev_rf-1, n_lev_csi
      h = log(co2col(ilev))
      call HUNT (uco2, n_Lesc, h, i)
      if (i == 0) i = 1
      if (i == n_Lesc) i = n_Lesc-1
      c_int = (h - uco2(i)) / (uco2(i+1)-uco2(i))
      Lesc_x(ilev) = min(max(Lescf(i) * (Lescf(i+1) / Lescf(i))**c_int, zero),one)
   end do

! Calculate the collisional non-LTE rates at atmos. temperatures
!----------------------------
   do ilev = n_lev_rf-1,n_lev
      t13  = temp_x(ilev)**(-onedthree)
      tsq  = sqrt(temp_x(ilev))
      zo   = a_zo  * tsq + b_zo  * exp(-g_zo *t13)
      zo2  = a_zo2 * tsq + b_zo2 * exp(-g_zo2*t13)
      zn2  = a_zn2 * tsq + b_zn2 * exp(-g_zn2*t13)
      lamb(ilev) = lmbfac / (lmbfac + n_dens(ilev) * (n2vmr_x(ilev)*zn2 + &
                   o2vmr_x(ilev)*zo2 + ovmr_x(ilev)*zo))
      cp(ilev) = rgasue7 / molmass(ilev) * &
                 (3.5_dp*(one - ovmr_x(ilev)) + 2.5_dp*ovmr_x(ilev))
      fac(ilev) = (numfac * co2vmr_x(ilev) * (one-lamb(ilev))) / molmass(ilev)
   end do

! Calculation using the recurrence formula (Eq. 8)
!-------------------
   alpha_ok = one
   eps_gn   = zero

   alpha_ok(n_lev_rf-1:n_lev_ru-1) = alphaf
   dj(n_lev_rf-1:n_lev) = Lesc_x(n_lev_rf-1:n_lev)*alpha_ok(n_lev_rf-1:n_lev)
   eps125 = hr(n_lev_rf-1) * cp(n_lev_rf-1) / sperday
   eps_gn(n_lev_rf-1) = eps125 / fac(n_lev_rf-1)

! Eq. 10
!-------------------------
   do ilev = n_lev_rf, n_lev_csi
      Djj   = 0.25_dp * (dj(ilev-1) + three * dj(ilev))
      Djjm1 = 0.25_dp * (dj(ilev) + three * dj(ilev-1))
      Fj    = one - lamb(ilev) * (one - Djj)
      Fjm1  = one - lamb(ilev-1) * (one - Djjm1)
      eps_gn(ilev) = (Fjm1 * eps_gn(ilev-1) + Djjm1 * &
                  phi_fun(ilev-1) - Djj * phi_fun(ilev)) / Fj
   end do

! Cooling to space
!-----------------
   if (n_lev > n_lev_cs) then
        Phi_165 = eps_gn(n_lev_cs) + phi_fun(n_lev_cs)
        eps_gn(n_lev_cs:n_lev) = Phi_165 - phi_fun(n_lev_cs:n_lev)
   endif

! Eq. 7
!-------------------------
   hr(n_lev_rf:n_lev) = fac(n_lev_rf:n_lev) * eps_gn(n_lev_rf:n_lev)

! Convert back the heating to K/day
!----------------------
   hr(n_lev_rf:n_lev) = hr(n_lev_rf:n_lev) * sperday / cp(n_lev_rf:n_lev)

end function CALC_COOL

!======================================================================!
! HUNT ---  find interval in xx where x is located in                  !
!----------------------------------------------------------------------!
! STATUS: bf 24-Feb-2024                     CREATED: gra 06-Sep-1999  !
!----------------------------------------------------------------------!
! INPUT:                                                               !
!       xx  --- vector, MONOTONIC real(dp) values                      !
!       n   --- number of values in xx, > 1                            !
!       x   --- real(dp) value to find interval for                    !
!       jlo --- last found interval index (see notes)                  !
!----------------------------------------------------------------------!
! OUTPUT:                                                              !
!       jlo --- lower index of interval (see notes)                    !
!----------------------------------------------------------------------!
! NOTES: This routine has the same result as locate,                   !
!        but does an additional hunting phase for                      !
!        searches for consecutive values. jlo should                   !
!        be set to a initial guess for the interval                    !
!        index, to 0 or n if hunting is not needed.                    !
!        See locate for further details.                               !
!----------------------------------------------------------------------!
! MAINTENANCE HISTORY:                                                 !
! bf  24-02-24: bi-directional hunting included                        !
! gra 06-09-99: adapted from Numerical Recipes                         !
!======================================================================!
subroutine HUNT (xx,n,x,jlo)
   implicit NONE

   real(dp), dimension(:), intent   (IN) :: xx
   integer               , intent   (IN) :: n
   real(dp)              , intent   (IN) :: x
   integer               , intent(INOUT) :: jlo

   integer :: inc,jhi,jm
   logical :: ascnd

!----------------!
!     hunt phase !
!----------------!

   ascnd = xx(n) >= xx(1)

   if (jlo <= 0 .or. jlo >= n) then
      jlo = 0
      jhi = n+1
   else
!
! --> hunt up
!
      inc = 1
      if (x > xx(jlo) .eqv. ascnd) then
         do
            jhi = jlo+inc
            if (jhi > n) then
               jhi = n+1
               exit
            else
               if (x <= xx(jhi) .eqv. ascnd) exit
               jlo = jhi
               inc = inc+inc
            end if
         end do
      else
!
! --> hunt down
!
         jhi = jlo
         do
            jlo = jhi-inc
            if (jlo < 1) then
               jlo = 0
               exit
            else
               if (x > xx(jlo) .eqv. ascnd) exit
               jhi = jlo
               inc = inc+inc
            end if
         end do
      end if
   end if
!---------------------!
!     bisection phase !
!---------------------!
   do while (jhi-jlo > 1)
      jm = (jlo+jhi)/2
      if (x > xx(jm) .eqv. ascnd) then
         jlo = jm
      else
         jhi = jm
      endif
   enddo

   if (x == xx(1)) jlo = 1

end subroutine HUNT

end module co2cool
