
      module mag_parms

      use shr_kind_mod,   only : r8 => shr_kind_r8
      use solar_parms_data, only : wkp=>solar_parms_kp, wf107=>solar_parms_f107
      use cam_abortutils, only : endrun

      implicit none

      private
      public :: get_mag_parms

      contains

      subroutine get_mag_parms( by, bz, hpower, ctpoten )
!---------------------------------------------------------------
!       ... retrieve magnetic field parmaters
!---------------------------------------------------------------

      implicit none

!---------------------------------------------------------------
!       ... dummy arguments
!---------------------------------------------------------------
      real(r8), optional, intent(out) :: by
      real(r8), optional, intent(out) :: bz
      real(r8), optional, intent(out) :: hpower
      real(r8), optional, intent(out) :: ctpoten

      if( present( by ) ) then
         by  =  0._r8
      end if
      if( present( bz ) ) then
         bz = .433726_r8 - wkp*(.0849999_r8*wkp + .0810363_r8) &
              + wf107*(.00793738_r8 - .00219316_r8*wkp)
      end if
! modified by LQIAN, 2008
! for wkp<=7: formula given by Zhang Yongliang based on TIMED/GUVI
! for wkp>7: power is 153.13 when wkp=7 from Zhang's formula, 
! assume power is 300.(based on NOAA satellites) when wkp=9 
! do linear interporation in between 
      if( present( hpower ) ) then
         if (wkp <=7._r8) hpower = 16.82_r8*exp(0.32_r8*wkp)-4.86_r8
         if (wkp > 7._r8) hpower = 153.13_r8+(wkp-7._r8)/ &
                               (9._r8-7._r8)*(300._r8-153.13_r8)

      end if
!
! modified by LQIAN, 2008
! formula given by Wenbin based on data fitting 
!
! 9/18/15 btf: If Weimer model was used for high-latitude potential,
!   then use ctpoten that was returned by weimer (see sub weimer05 in
!   wei05sc.F90, called by dpie_coupling.F90)
!
      if( present( ctpoten ) ) then
         ctpoten = 15._r8+15._r8*wkp + 0.8_r8*wkp**2
      end if

      end subroutine get_mag_parms

      end module mag_parms
