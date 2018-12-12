! NCLFORTSTART
!-----------------------------------------------------------------------
!     subroutine read_namelist(Nudge_Model,Nudge_Path,                 &
!                        Nudge_File_Template,Nudge_Times_Per_Day,      &
!                        Model_Times_Per_Day,                          &
!                        Nudge_Ucoef ,Nudge_Uprof,                     &
!                        Nudge_Vcoef ,Nudge_Vprof,                     &
!                        Nudge_Qcoef ,Nudge_Qprof,                     &
!                        Nudge_Tcoef ,Nudge_Tprof,                     &
!                        Nudge_PScoef,Nudge_PSprof,                    &
!                        Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day, &
!                        Nudge_End_Year,Nudge_End_Month,Nudge_End_Day, &
!                        Nudge_Hwin_Invert,                            &
!                        Nudge_Vwin_Invert,                            &
!                        Nudge_Hwin_lat0,Nudge_Hwin_lon0,              &
!                        Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,      &
!                        Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,      &
!                        Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,          &
!                        Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta )
      subroutine read_namelist(Nudge_Hwin_Invert,                      &
                               Nudge_Vwin_Invert,                      &
                               Nudge_Hwin_lat0,Nudge_Hwin_lon0,        &
                               Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,&
                               Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,&
                               Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,    &
                               Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta )
      !---------------------------------------------------------------
      ! Read Processing Namelist
      !---------------------------------------------------------------
      implicit none

      !     INPUTS
      !---------------------------------------------------------------

      !     OUTPUT
      !---------------------------------------------------------------
      logical::         Nudge_Model       
      character*80      Nudge_Path
      character*80      Nudge_File_Template
      integer           Nudge_Force_Opt
      integer           Nudge_TimeScale_Opt
      integer           Nudge_Times_Per_Day
      integer           Model_Times_Per_Day
      real              Nudge_Ucoef,Nudge_Vcoef
      integer           Nudge_Uprof,Nudge_Vprof
      real              Nudge_Qcoef,Nudge_Tcoef
      integer           Nudge_Qprof,Nudge_Tprof
      real              Nudge_PScoef
      integer           Nudge_PSprof
      integer           Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day
      integer           Nudge_End_Year,Nudge_End_Month,Nudge_End_Day
      logical           Nudge_Hwin_Invert
      real              Nudge_Hwin_lat0
      real              Nudge_Hwin_latWidth
      real              Nudge_Hwin_latDelta
      real              Nudge_Hwin_lon0
      real              Nudge_Hwin_lonWidth
      real              Nudge_Hwin_lonDelta
      logical           Nudge_Vwin_Invert
      real              Nudge_Vwin_Hindex
      real              Nudge_Vwin_Hdelta
      real              Nudge_Vwin_Lindex
      real              Nudge_Vwin_Ldelta

! NCLEND
!-----------------------------------------------------------------------
      namelist /nudging_nl/ Nudge_Model,Nudge_Path,                       &
                            Nudge_File_Template,Nudge_Times_Per_Day,      &
                            Model_Times_Per_Day,                          &
                            Nudge_Force_Opt,                              &
                            Nudge_TimeScale_Opt,                          &
                            Nudge_Ucoef ,Nudge_Uprof,                     &
                            Nudge_Vcoef ,Nudge_Vprof,                     &
                            Nudge_Qcoef ,Nudge_Qprof,                     &
                            Nudge_Tcoef ,Nudge_Tprof,                     &
                            Nudge_PScoef,Nudge_PSprof,                    &
                            Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day, &
                            Nudge_End_Year,Nudge_End_Month,Nudge_End_Day, &
                            Nudge_Hwin_Invert,                            &
                            Nudge_Vwin_Invert,                            &
                            Nudge_Hwin_lat0,Nudge_Hwin_lon0,              &
                            Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,      &
                            Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,      &
                            Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,          &
                            Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta

!     namelist /nudging_nl/ Nudge_Hwin_Invert,                        &
!                           Nudge_Vwin_Invert,                        &
!                           Nudge_Hwin_lat0,Nudge_Hwin_lon0,          &
!                           Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,  &
!                           Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,  &
!                           Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,      &
!                           Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta

      ! open namelist file and read in values
      !-------------------------------------------
      open(unit=31,file='./user_nl_cam',status='old')
      read(31,nudging_nl)
      close(unit=31)

      ! End Routine
      !-------------
      return
      end
