module get_ana_dynfrc_4scam

  use spmd_utils,       only: masterproc
  use cam_logfile,      only: iulog
  use shr_kind_mod,     only: r8 => shr_kind_r8, i8 => shr_kind_i8, &
                          cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use shr_const_mod,    only:  rearth => shr_const_rearth , &  !  =6.37122e6_R8 meters
                               pi => shr_const_pi         , &
                               OOmega => shr_const_omega   , &
                               rdair => shr_const_rdair   , &
                               cpair => shr_const_cpdair

  use scamMod,          only:  use_scm_ana_frc, &
                               scm_ana_frc_path, &
                               scm_ana_frc_file_template, &
                               scm_ana_x_plevels, &
                               scm_ana_direct_omega, &
                               scm_ana_t_react, &
                               scm_ana_q_react, &
                               scm_ana_u_react, &
                               scm_ana_v_react, &
                               scm_ana_upwind, &
                               scm_ana_direct_ttend



  !    shr_const_mod is in ${CESMROOT}/cime/src/share/util/

  implicit none
  private
  save

  public get_ana_dynfrc_fv
!
! Private module data
!

    real(r8) , save , allocatable :: T_1(:,:,:) , U_1(:,:,:), V_1(:,:,:), Q_1(:,:,:),PS_1(:,:),PHIS_1(:,:)
    real(r8) , save , allocatable :: T_2(:,:,:) , U_2(:,:,:), V_2(:,:,:), Q_2(:,:,:),PS_2(:,:),PHIS_2(:,:)
    real(r8) , save , allocatable :: UTCORE_1(:,:,:) , UTCORE_2(:,:,:)
    real(r8) , save , allocatable :: VTCORE_1(:,:,:) , VTCORE_2(:,:,:)
    real(r8) , save , allocatable :: TTCORE_1(:,:,:) , TTCORE_2(:,:,:)
    real(r8) , save , allocatable :: OGCORE_1(:,:,:) , OGCORE_2(:,:,:)
    real(r8) , save , allocatable :: lat_ana(:),lon_ana(:),lev_ana(:)
    integer  , save               :: nlev_ana, nlon_ana, nlat_ana

    real(r8) , save , allocatable :: To_1(:,:,:) , Uo_1(:,:,:), Vo_1(:,:,:), Qo_1(:,:,:),PSo_1(:,:),PHISo_1(:,:)
    real(r8) , save , allocatable :: To_2(:,:,:) , Uo_2(:,:,:), Vo_2(:,:,:), Qo_2(:,:,:),PSo_2(:,:),PHISo_2(:,:)
    real(r8) , save , allocatable :: UTCOREo_1(:,:,:) , UTCOREo_2(:,:,:), UTCOREo_X(:,:,:)
    real(r8) , save , allocatable :: VTCOREo_1(:,:,:) , VTCOREo_2(:,:,:), VTCOREo_X(:,:,:)
    real(r8) , save , allocatable :: TTCOREo_1(:,:,:) , TTCOREo_2(:,:,:), TTCOREo_X(:,:,:)
    real(r8) , save , allocatable :: OGCOREo_1(:,:,:) , OGCOREo_2(:,:,:), OGCOREo_X(:,:,:)



    real(r8) , save , allocatable :: ETAD_X(:,:,:) , OMG_X(:,:,:)
    real(r8) , save , allocatable :: ZETA_X(:)
    real(r8) , save , allocatable :: KEh_X(:,:,:)
    real(r8) , save , allocatable :: Tv_X(:,:,:)

    real(r8) , save , allocatable :: pke_X(:,:,:),pko_X(:,:,:),phik_X(:,:,:),Thv_X(:,:,:)
    real(r8) , save , allocatable :: ple_X(:,:,:) , plo_X(:,:,:), phi_X(:,:,:)

    real(r8) , save , allocatable :: To_X(:,:,:) , Uo_X(:,:,:), Vo_X(:,:,:), Qo_X(:,:,:),PSo_X(:,:),PHISo_X(:,:)


!=======================================================================
contains
!=======================================================================

subroutine get_ana_dynfrc_fv ( scmlon, scmlat ,  &
                               omega_ana, etad_ana, zeta_ana, &
                               t_ana , tv_ana ,  &
                               q_ana ,           &
                               u_ana ,           &
                               v_ana ,           &
                               ps_ana ,          &
                               uten_hadv_ana ,   &
                               vten_hadv_ana ,   &
                               uten_pfrc_ana ,   &
                               vten_pfrc_ana ,   &
                               uten_vort_ana ,   &
                               vten_vort_ana ,   &
                               qten_hadv_ana ,   &
                               tten_hadv_ana ,   &
                               uten_vadv_ana ,   &
                               vten_vadv_ana ,   &
                               tten_vadv_ana ,   &
                               qten_vadv_ana ,   &
                               tten_comp_ana ,   &
                               uten_keg_ana  ,   &
                               uten_phig_ana ,   &
                               uten_prg_ana  ,   &
                               uten_dycore_ana , &
                               vten_dycore_ana , &
                               tten_dycore_ana , &
                               omega_dycore_ana , &
                               omega_recalc_ana , &
                               u_scm, v_scm, t_scm, q_scm,  &
                               u_ana_diag, v_ana_diag, t_ana_diag, q_ana_diag   )
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! US and VS are input (D-grid velocities)
!--------------------------------------------
!  ub(i,j,L)= 0.5*(us(i-1,j,L) + us(i,j,L))
!  vb(i,j,L)= 0.5*(vs(i,j,L)   + vs(i,j+1,L))
!
!  uc(i,j,L)= 0.5*(ub(i,j,L)   + ub(i,j-1,L))
!  vc(i,j,L)= 0.5*(vb(i,j-1,L) + vb(i+1,j-1,L))
!---------------------------------------------
!   Grid arrangement in FV latlon h,i-files
!---------------------------------------------
!    J=NY
!                   ...
!
!     ub,vb(I,J)           us(I,J),vc(I,J+1)
!
!
!     vs(I,J),uc(I,J)      ua,va,T,p(I,J)       vs(I+1,J),uc(I+1,J)
!
!
!                          vc(I,J)
!
!
!                          ua,va,T,p(I,J-1)
!
!                    ...
!    J=1 ...
!----------------------------------------------

     use pmgrid,         only : plev, plat, plevp, plon
     use hycoef,         only: hyai, hybi, ps0, hyam, hybm
     use filenames,      only: interpret_filename_spec
     use time_manager,   only: timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size,is_first_step
     use netcdf
     use cam_abortutils, only: endrun
     use ref_pres,       only: pref_mid  ! In Pascal

    real(r8), intent(in)   :: scmlon, scmlat
    real(r8), intent(out)  :: omega_ana( plev )
    real(r8), intent(out)  :: etad_ana(plev)
    real(r8), intent(out)  :: t_ana(plev) , tv_ana(plev)
    real(r8), intent(out)  :: zeta_ana(plev)
    real(r8), intent(out)  :: u_ana(plev)
    real(r8), intent(out)  :: v_ana(plev)
    real(r8), intent(out)  :: q_ana(plev)
    real(r8), intent(out)  :: ps_ana
    real(r8), intent(out)  :: uten_hadv_ana( plev )
    real(r8), intent(out)  :: vten_hadv_ana( plev )
    real(r8), intent(out)  :: uten_pfrc_ana( plev )
    real(r8), intent(out)  :: vten_pfrc_ana( plev )
    real(r8), intent(out)  :: qten_hadv_ana( plev )
    real(r8), intent(out)  :: tten_hadv_ana( plev )
    real(r8), intent(out)  :: qten_vadv_ana( plev )
    real(r8), intent(out)  :: tten_vadv_ana( plev )
    real(r8), intent(out)  :: uten_vadv_ana( plev )
    real(r8), intent(out)  :: vten_vadv_ana( plev )

    real(r8), intent(out)  :: tten_comp_ana( plev )

    real(r8), intent(out)  :: uten_keg_ana( plev )
    real(r8), intent(out)  :: uten_prg_ana( plev )
    real(r8), intent(out)  :: uten_phig_ana( plev )
    real(r8), intent(out)  :: uten_vort_ana( plev )
    real(r8), intent(out)  :: vten_vort_ana( plev )
    real(r8), intent(out)  :: uten_dycore_ana( plev )
    real(r8), intent(out)  :: vten_dycore_ana( plev )
    real(r8), intent(out)  :: tten_dycore_ana( plev )
    real(r8), intent(out)  :: omega_recalc_ana( plev )
    real(r8), intent(out)  :: omega_dycore_ana( plev )

    real(r8), intent(in) :: u_scm(plev)
    real(r8), intent(in) :: v_scm(plev)
    real(r8), intent(in) :: t_scm(plev)
    real(r8), intent(in) :: q_scm(plev)

    real(r8), intent(out) :: u_ana_diag(plev)
    real(r8), intent(out) :: v_ana_diag(plev)
    real(r8), intent(out) :: t_ana_diag(plev)
    real(r8), intent(out) :: q_ana_diag(plev)

    integer, save :: iax, jax
    integer, save :: Read_year2, Read_month2, Read_day2, Read_sec2, Read_YMD2
    integer, save :: nlev_alc, nlon_alc, nlat_alc

    !!logical , parameter   :: l_vectinv = .FALSE.
    !!real(r8) :: tv_ana(plev)
    real(r8) :: rho_ana( plev ), plo_ana(plev)



    real(r8) :: scmlonx

    real(r8) :: ana_wgt1 , ana_wgt2 , dx0, dy, darea

    integer :: nx, ny,i,j,k,L,LM, iav(1),jav(1),iac,jac

    real(r8) , allocatable :: rlats(:),rlons(:)
    real(r8) :: zeta(plev),absvo(plev)
    ! Horz. gradient profiles (1=X, 2=Y)
    real(r8) :: kehg_ana(plev,2),kehg_X(plev,2)
    real(r8) :: phig_ana(plev,2),phig_X(plev,2)
    real(r8) :: plog_ana(plev,2),plog_X(plev,2)
    real(r8) :: teg_ana(plev,2), teg_X(plev,2)
    real(r8) :: qg_ana(plev,2),  qg_X(plev,2)
    real(r8) :: ug_ana(plev,2),  ug_X(plev,2)
    real(r8) :: vg_ana(plev,2),  vg_X(plev,2)
    real(r8) :: lin_pfc_ana(plev,2) , lin_pfc_X(plev,2)

    real(r8) :: omega_ana_x(plev)
    real(r8) :: alpha_react(plev)

    real(r8) :: lat_alc(3) , lon_alc(3)
    real(r8) :: aalc(3,3,plev)


    character(len=CL):: Ana_File_Template,Ana_file1,Ana_file2,Ana_Path


    integer :: dyn_year,dyn_month,dyn_day,dyn_sec,year,month,day,sec
    integer :: dyn_step,ymd1,ymd2,curr_sec,next_sec,curr_year,curr_month,curr_day,curr_ymd

    integer :: analysis_step
    integer :: ana_year1, ana_month1, ana_day1, ana_sec1
    integer :: ana_year2, ana_month2, ana_day2, ana_sec2

    logical :: l_Read_next_Ana, Alarm_Read_ana, Alarm_Bump_ana, initialize

    character(len=19) :: subname='get_ana_dynfrc_fv: '

    write(iulog,*) subname//" version 07 of get_ana_dynfrc_4scam ... "


    Alarm_Read_Ana = .FALSE.
    Alarm_Bump_Ana = .FALSE.

    if ( scmlon < 0 ) then
       scmlonx = scmlon + 360._r8
    else
       scmlonx = scmlon
    end if

    ! Default to 6 hour steps between ana
    analysis_step = 6 * 3600


    Ana_path = trim(scm_ana_frc_path)
    Ana_File_Template = trim(Ana_path)//trim(scm_ana_frc_file_template)


    call get_curr_date(Year,Month,Day,Sec)

    curr_ymd = (Year*10000) + (Month*100) + Day
    curr_sec = Sec

    ana_sec1   = ( Sec / analysis_step ) * analysis_step
    ana_day1   = Day
    ana_month1 = Month
    ana_year1  = Year

     YMD1=(Ana_Year1*10000) + (Ana_Month1*100) + Ana_Day1


     call timemgr_time_inc(YMD1,Ana_Sec1,              &
                           YMD2,Ana_Sec2,Analysis_Step,0,0)

     Ana_Year2   = YMD2 / 10000
     Ana_Month2  = (YMD2 - Ana_Year2*10000)/100
     Ana_Day2    = YMD2 - Ana_Year2*10000 - Ana_Month2*100

     Ana_File1   = interpret_filename_spec(Ana_File_Template      , &
                                       yr_spec=Ana_Year1 , &
                                      mon_spec=Ana_Month1, &
                                      day_spec=Ana_Day1  , &
                                      sec_spec=Ana_Sec1    )

     Ana_File2   = interpret_filename_spec(Ana_File_Template      , &
                                       yr_spec=Ana_Year2 , &
                                      mon_spec=Ana_Month2, &
                                      day_spec=Ana_Day2  , &
                                      sec_spec=Ana_Sec2    )


       l_Read_next_Ana = .FALSE.
       ! On first time step, read in 2 analysis files
       if (is_first_step().and.masterproc) then
              write(iulog,*) subname//" It's now (First time step):" , curr_YMD, curr_sec
              write(iulog,*) "Read Initial ana files "
              write(iulog,*) Ana_file1
              write(iulog,*) Ana_file2
              Alarm_Read_Ana = .TRUE.
              Alarm_Bump_Ana = .FALSE.
       else
       ! On subsequent steps test to see if "Curr" date is later or same as "Read".
       ! If it is, then l_read_next_ana=.TRUE.
             call timemgr_time_ge(Read_ymd2, Read_Sec2, curr_YMD, curr_Sec, l_Read_next_ana )
       endif

       if (l_Read_next_Ana) then
              Alarm_Read_Ana = .TRUE.
              Alarm_Bump_Ana = .TRUE.
       endif

       ! Aloocate space for analysis fields.
       ! Read in both Initial Analysis files. Nothing to bump yet
       if (  (Alarm_Read_Ana ) .AND. .NOT.(Alarm_Bump_Ana) ) then
           initialize=.TRUE.
           call read_netcdf_ana_fv_ini ( Ana_File1, nlon_ana, nlat_ana, nlev_ana ,iax, jax   )

              if ( plev /= nlev_ana) then
                 call endrun ("SCAM plev NE nlev_ana")
              end if

              ! Full global fields
              allocate( lat_ana(nlat_ana) , lon_ana(nlon_ana), lev_ana(nlev_ana) )
              allocate( U_1(nlon_ana, nlat_ana, nlev_ana), V_1(nlon_ana, nlat_ana, nlev_ana), T_1(nlon_ana, nlat_ana, nlev_ana), &
                        Q_1(nlon_ana, nlat_ana, nlev_ana), PS_1 (nlon_ana, nlat_ana ), PHIS_1 (nlon_ana, nlat_ana ) )
              allocate( U_2(nlon_ana, nlat_ana, nlev_ana), V_2(nlon_ana, nlat_ana, nlev_ana), T_2(nlon_ana, nlat_ana, nlev_ana), &
                        Q_2(nlon_ana, nlat_ana, nlev_ana), PS_2 (nlon_ana, nlat_ana ), PHIS_2 (nlon_ana, nlat_ana ) )

              allocate( UTCORE_1(nlon_ana, nlat_ana, nlev_ana), UTCORE_2(nlon_ana, nlat_ana, nlev_ana) )
              allocate( VTCORE_1(nlon_ana, nlat_ana, nlev_ana), VTCORE_2(nlon_ana, nlat_ana, nlev_ana) )
              allocate( TTCORE_1(nlon_ana, nlat_ana, nlev_ana), TTCORE_2(nlon_ana, nlat_ana, nlev_ana) )
              allocate( OGCORE_1(nlon_ana, nlat_ana, nlev_ana), OGCORE_2(nlon_ana, nlat_ana, nlev_ana) )

              ! SCM "patches"
              nlon_alc=3
              nlat_alc=3
              nlev_alc=nlev_ana



              ! Patches of full global fields
              allocate( Uo_1(nlon_alc, nlat_alc, nlev_alc), Vo_1(nlon_alc, nlat_alc, nlev_alc), To_1(nlon_alc, nlat_alc, nlev_alc), &
                        Qo_1(nlon_alc, nlat_alc, nlev_alc), PSo_1 (nlon_alc, nlat_alc ), PHISo_1 (nlon_alc, nlat_alc ) )
              allocate( Uo_2(nlon_alc, nlat_alc, nlev_alc), Vo_2(nlon_alc, nlat_alc, nlev_alc), To_2(nlon_alc, nlat_alc, nlev_alc), &
                        Qo_2(nlon_alc, nlat_alc, nlev_alc), PSo_2 (nlon_alc, nlat_alc ), PHISo_2 (nlon_alc, nlat_alc ) )

              allocate( UTCOREo_1(nlon_alc, nlat_alc, nlev_alc), UTCOREo_2(nlon_alc, nlat_alc, nlev_alc), UTCOREo_X(nlon_alc, nlat_alc, nlev_alc) )
              allocate( VTCOREo_1(nlon_alc, nlat_alc, nlev_alc), VTCOREo_2(nlon_alc, nlat_alc, nlev_alc), VTCOREo_X(nlon_alc, nlat_alc, nlev_alc) )
              allocate( TTCOREo_1(nlon_alc, nlat_alc, nlev_alc), TTCOREo_2(nlon_alc, nlat_alc, nlev_alc), TTCOREo_X(nlon_alc, nlat_alc, nlev_alc) )
              allocate( OGCOREo_1(nlon_alc, nlat_alc, nlev_alc), OGCOREo_2(nlon_alc, nlat_alc, nlev_alc), OGCOREo_X(nlon_alc, nlat_alc, nlev_alc) )

              allocate( Uo_X(nlon_alc, nlat_alc, nlev_alc), Vo_X(nlon_alc, nlat_alc, nlev_alc), To_X(nlon_alc, nlat_alc, nlev_alc), &
                        Qo_X(nlon_alc, nlat_alc, nlev_alc), PSo_X (nlon_alc, nlat_alc ), PHISo_X (nlon_alc, nlat_alc ) )
              allocate( ETAD_X(nlon_alc,nlat_alc,nlev_alc) )
              allocate( OMG_X(nlon_alc,nlat_alc,nlev_alc)  )
              allocate( ple_X(nlon_alc, nlat_alc, nlev_alc+1), plo_X(nlon_alc, nlat_alc, nlev_alc), phi_X(nlon_alc, nlat_alc, nlev_alc+1) )
              allocate( pke_X(nlon_alc, nlat_alc, nlev_alc+1), pko_X(nlon_alc, nlat_alc, nlev_alc), phik_X(nlon_alc, nlat_alc, nlev_alc+1) )
              allocate( THv_X(nlon_alc, nlat_alc, nlev_alc ) )
              allocate( zeta_X(nlev_alc) )
              allocate( KEh_X(nlon_alc, nlat_alc, nlev_alc ) )
              allocate( Tv_X(nlon_alc, nlat_alc, nlev_alc ) )

           call read_netcdf_ana_fv ( Ana_File1, nlon_ana, nlat_ana, nlev_ana, &
                                   U_1, V_1, &
                                   T_1, Q_1, PS_1, PHIS_1, &
                                   lon_ana, lat_ana, lev_ana &
                               ,   utcore_1, vtcore_1, ttcore_1, ogcore_1 &
                                                )

           call read_netcdf_ana_fv ( Ana_File2, nlon_ana, nlat_ana, nlev_ana, &
                                   U_2, V_2, &
                                   T_2, Q_2, PS_2, PHIS_2, &
                                   lon_ana, lat_ana, lev_ana &
                                ,  utcore_2, vtcore_2, ttcore_2, ogcore_2 &
                                                )

           ! Make patches
           Uo_1      = U_1(iax-1:iax+1,jax-1:jax+1,:)
           Vo_1      = V_1(iax-1:iax+1,jax-1:jax+1,:)
           To_1      = T_1(iax-1:iax+1,jax-1:jax+1,:)
           Qo_1      = Q_1(iax-1:iax+1,jax-1:jax+1,:)
           PSo_1     = PS_1(iax-1:iax+1,jax-1:jax+1 )
           PHISo_1   = PHIS_1(iax-1:iax+1,jax-1:jax+1 )
           UTCOREo_1 = UTCORE_1(iax-1:iax+1,jax-1:jax+1,:)
           VTCOREo_1 = VTCORE_1(iax-1:iax+1,jax-1:jax+1,:)
           TTCOREo_1 = TTCORE_1(iax-1:iax+1,jax-1:jax+1,:)
           OGCOREo_1 = OGCORE_1(iax-1:iax+1,jax-1:jax+1,:)

           Uo_2      = U_2(iax-1:iax+1,jax-1:jax+1,:)
           Vo_2      = V_2(iax-1:iax+1,jax-1:jax+1,:)
           To_2      = T_2(iax-1:iax+1,jax-1:jax+1,:)
           Qo_2      = Q_2(iax-1:iax+1,jax-1:jax+1,:)
           PSo_2     = PS_2(iax-1:iax+1,jax-1:jax+1 )
           PHISo_2   = PHIS_2(iax-1:iax+1,jax-1:jax+1 )
           UTCOREo_2 = UTCORE_2(iax-1:iax+1,jax-1:jax+1,:)
           VTCOREo_2 = VTCORE_2(iax-1:iax+1,jax-1:jax+1,:)
           TTCOREo_2 = TTCORE_2(iax-1:iax+1,jax-1:jax+1,:)
           OGCOREo_2 = OGCORE_2(iax-1:iax+1,jax-1:jax+1,:)


           ! Mark Ana date as read
           Read_year2    = Ana_year2
           Read_month2   = Ana_month2
           Read_day2     = Ana_day2
           Read_sec2     = Ana_sec2
           Read_YMD2     =(Ana_Year2*10000) + (Ana_Month2*100) + Ana_Day2

       end if

       ! Bump second analysis to first postion, and read in next analysis
       if (  (Alarm_Read_Ana ) .AND. (Alarm_Bump_Ana) ) then

           Uo_1      = Uo_2
           Vo_1      = Vo_2
           To_1      = To_2
           Qo_1      = Qo_2
           PSo_1     = PSo_2
           PHISo_1   = PHISo_2
           UTCOREo_1 = UTCOREo_2
           VTCOREo_1 = VTCOREo_2
           TTCOREo_1 = TTCOREo_2

           call read_netcdf_ana_fv ( Ana_File2, nlon_ana, nlat_ana, nlev_ana, &
                                   U_2, V_2, &
                                   T_2, Q_2, PS_2, PHIS_2, &
                                   lon_ana, lat_ana, lev_ana &
                                ,  utcore_2, vtcore_2, ttcore_2, ogcore_2 &
                                   )


           ! Make patches
           Uo_2      = U_2(iax-1:iax+1,jax-1:jax+1,:)
           Vo_2      = V_2(iax-1:iax+1,jax-1:jax+1,:)
           To_2      = T_2(iax-1:iax+1,jax-1:jax+1,:)
           Qo_2      = Q_2(iax-1:iax+1,jax-1:jax+1,:)
           PSo_2     = PS_2(iax-1:iax+1,jax-1:jax+1 )
           PHISo_2   = PHIS_2(iax-1:iax+1,jax-1:jax+1 )
           UTCOREo_2 = UTCORE_2(iax-1:iax+1,jax-1:jax+1,:)
           VTCOREo_2 = VTCORE_2(iax-1:iax+1,jax-1:jax+1,:)
           TTCOREo_2 = TTCORE_2(iax-1:iax+1,jax-1:jax+1,:)
           OGCOREo_2 = OGCORE_2(iax-1:iax+1,jax-1:jax+1,:)


           ! Mark Ana date as read
           Read_year2    = Ana_year2
           Read_month2   = Ana_month2
           Read_day2     = Ana_day2
           Read_sec2     = Ana_sec2
           Read_YMD2=(Ana_Year2*10000) + (Ana_Month2*100) + Ana_Day2
       end if

       Alarm_Read_Ana = .FALSE.
       Alarm_Bump_Ana = .FALSE.




#if 0
     call dynfrc_timewgts(      &
                         (/ Ana_Year1, Ana_Month1, Ana_day1, Ana_sec1 /) , &
                         (/ Ana_Year2, Ana_Month2, Ana_day2, Ana_sec2 /) , &
                         ana_wgt1 , ana_wgt2   )
#else
     ana_wgt1 = 0._r8  ! 0=all weight on t+1
     ana_wgt2 = 1._r8 - ana_wgt1
#endif
     if (masterproc) write(iulog,*) subname//" Ana forcing time wgts ",ana_wgt1,ana_wgt2

          iac=2
          jac=2



          Uo_X      = ana_wgt1 * Uo_1      +   ana_wgt2 * Uo_2
          Vo_X      = ana_wgt1 * Vo_1      +   ana_wgt2 * Vo_2
          To_X      = ana_wgt1 * To_1      +   ana_wgt2 * To_2
          Qo_X      = ana_wgt1 * Qo_1      +   ana_wgt2 * Qo_2
          PSo_X     = ana_wgt1 * PSo_1     +   ana_wgt2 * PSo_2
          PHISo_X   = ana_wgt1 * PHISo_1   +   ana_wgt2 * PHISo_2

          UTCOREo_X = ana_wgt1 * UTCOREo_1 +   ana_wgt2 * UTCOREo_2
          VTCOREo_X = ana_wgt1 * VTCOREo_1 +   ana_wgt2 * VTCOREo_2
          TTCOREo_X = ana_wgt1 * TTCOREo_1 +   ana_wgt2 * TTCOREo_2
          OGCOREo_X = ana_wgt1 * OGCOREo_1 +   ana_wgt2 * OGCOREo_2

          lon_alc = lon_ana(iax-1:iax+1)
          lat_alc = lat_ana(jax-1:jax+1)

          if(masterproc) write(iulog,*) subname//"         SCM lon lat: ",scmlonx,scmlat
          if(masterproc) write(iulog,*) subname//" Closest Ana lon lat: ",lon_ana( iax ) , lat_ana( jax )


          ! Save off analysis fields for diagnostics and
          ! other purposes
          T_ana_diag(:) = To_X( iac, jac, :)
          Q_ana_diag(:) = Qo_X( iac, jac, :)
          U_ana_diag(:) = Uo_X( iac, jac, :)
          V_ana_diag(:) = Vo_X( iac, jac, :)

          !================================================
          ! Patch in SCM profiles here if wanted.
          ! This acts as "dynamical nudging", since
          ! horizontal advective tendencies will become
          ! stronger if SCM state drifts away from re-ana.
          ! Note, this will only be effective w/ upwind
          ! scheme, since 2nd order cntrd skips over central
          ! point in stencil.
          !----
          ! For stability it turns out may be good to scale
          ! with pressure so that high-velocity strato winds
          ! don't lead to CFL violations. So, as a bad, dirty,
          ! dirty short term solution, weight "reaction" by
          ! pref_mid.  Clearly, better soln would be to
          ! sub-step this part of the dynamics as is done
          ! for the other "dycores".
          !=================================================
          ! Calculate "reaction coefficient"
          !---------------------------------
          alpha_react(:)=1.0_r8 !1._r8

          ! Adjust central profiles in stencils
          !------------------------------------
          if (scm_ana_t_react) then
              To_X( iac, jac, :) = alpha_react(:) * T_scm(:)   &
                                + ( 1._r8-alpha_react(:) ) * To_X( iac, jac, :)
              if(masterproc) write(iulog,*) subname//" REACTING to SCM T-state ..... "
          else
              if(masterproc) write(iulog,*) subname//" No reaction to SCM T-state ..... "
          endif
          if (scm_ana_q_react) then
              Qo_X( iac, jac, :) = alpha_react(:) * Q_scm(:) &
                                + ( 1._r8-alpha_react(:) ) * Qo_X( iac, jac, :)
              if(masterproc) write(iulog,*) subname//" REACTING to SCM Q-state ..... "
          else
              if(masterproc) write(iulog,*) subname//" No reaction to SCM Q-state ..... "
          endif
          if (scm_ana_u_react) then
              Uo_X( iac, jac, :) = alpha_react(:) * U_scm(:) &
                                + ( 1._r8-alpha_react(:) ) * Uo_X( iac, jac, :)
              if(masterproc) write(iulog,*) subname//" REACTING to SCM U-state ..... "
          else
              if(masterproc) write(iulog,*) subname//" No reaction to SCM U-state ..... "
          endif
          if (scm_ana_v_react) then
              Vo_X( iac, jac, :) = alpha_react(:) * V_scm(:) &
                                + ( 1._r8-alpha_react(:) ) * Vo_X( iac, jac, :)
              if(masterproc) write(iulog,*) subname//" REACTING to SCM V-state ..... "
          else
              if(masterproc) write(iulog,*) subname//" No reaction to SCM V-state ..... "
          endif



          !=========================================

           call virtual_t(  nlon_alc,nlat_alc,nlev_alc, &
                             To_X , Qo_X , Tv_X )

           call makepr_fv(  nlon_alc,nlat_alc,nlev_alc, &
                            tv_X , pso_X , phiso_X , &
                            plo_X, ple_X, phi_X )
           call etadot_fv ( nlon_alc , nlat_alc , nlev_alc , lon_alc , lat_alc , &
                            uo_X ,  &
                            vo_X , &
                            plo_X, ple_X , etad_X , omg_X )
           call zeta_fv( nlon_alc,nlat_alc,nlev_alc, &
                       lon_alc ,lat_alc , &
                       uo_X , vo_X , zeta_X )

           call makepk_fv( nlon_alc,nlat_alc,nlev_alc, &
                         To_X , Qo_X ,  &
                           pso_X , phiso_X , &
                            pko_X, pke_X, phik_X, thv_X )

           KEh_X  = 0.5 * ( Uo_X**2 + Vo_X**2 )


          if (scm_ana_x_plevels) then
             call patch_eta_x_plv (  nlon_alc , nlat_alc ,  nlev_alc, &
                                  iac, jac, uo_X , plo_X )
             call patch_eta_x_plv (  nlon_alc , nlat_alc ,  nlev_alc, &
                                  iac, jac, vo_X , plo_X )
             call patch_eta_x_plv (  nlon_alc , nlat_alc ,  nlev_alc, &
                                  iac, jac, to_X , plo_X )
             call patch_eta_x_plv (  nlon_alc , nlat_alc ,  nlev_alc, &
                                  iac, jac, qo_X , plo_X )
             call patch_eta_x_plv (  nlon_alc , nlat_alc ,  nlev_alc, &
                                  iac, jac, tv_X , plo_X )
             !Retain p-frc calculation on eta???
             !call patch_eta_x_plv (  nlon_alc , nlat_alc ,  nlev_alc+1, &
             !                     iac, jac, phi_X , ple_X )
             if(masterproc) write(iulog,*) subname//" calcs on PRESSURE levels "
          else
             if(masterproc) write(iulog,*) subname//" calcs on ETA levels "
          end if


          zeta_ana         = zeta_X
          omega_recalc_ana = omg_X( iac,jac,:)
          etad_ana         = etad_X( iac,jac,:)
          plo_ana   = plo_X( iac,jac,:)
          t_ana     = To_X( iac,jac,:)
          tv_ana    = Tv_X( iac,jac,:)
          q_ana     = Qo_X( iac,jac,:)
          ps_ana    = PSo_X( iac,jac )

          u_ana     = Uo_X( iac,jac,:)
          v_ana     = Vo_X( iac,jac,:)

          rho_ana   = plo_ana / ( Rdair * tv_ana )

          uten_dycore_ana     = UTCOREo_X( iac,jac,:)
          vten_dycore_ana     = VTCOREo_X( iac,jac,:)
          tten_dycore_ana     = TTCOREo_X( iac,jac,:)
          omega_dycore_ana    = OGCOREo_X( iac,jac,:)


          ! Horz. gradient calcs

            kehg_X = grad_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, KEh_X )

                 ! T_x, T_y  should be straight T (not virtual)
            teg_X  =  grad_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, Tv_X ) !test 05-31-21

            qg_X  =  grad_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, Qo_X )

            ug_X  =  grad_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, Uo_X )

            vg_X  =  grad_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, Vo_X )

            aalc   =  0.5*( PHI_X( :, :, 2:nlev_alc+1) +  PHI_X(: , : ,1:nlev_alc) )
            phig_X = grad_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, aalc )

            plog_X = grad_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, plo_X(:,:,1:nlev_alc) )


            lin_pfc_X = lin_pfc_fv( nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, ple_X, phi_X )

            kehg_ana  = kehg_X
            plog_ana  = plog_X
            phig_ana  = phig_X
            teg_ana   = teg_X
            qg_ana    = qg_X
            ug_ana    = ug_X
            vg_ana    = vg_X
            lin_pfc_ana  = lin_pfc_X

            !put together pieces for u*grad(u) form of U and V adv tendencies

            if ( scm_ana_upwind .OR. scm_ana_u_react ) then
              uten_hadv_ana =  upwind_hadv(nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, u_ana, v_ana, Uo_X )
            else
              uten_hadv_ana  =  -u_ana * ug_ana(:,1) - v_ana * ug_ana(:,2)
            end if
            if ( scm_ana_upwind .OR. scm_ana_v_react ) then
              vten_hadv_ana =  upwind_hadv(nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, u_ana, v_ana, Vo_X )
            else
              vten_hadv_ana  =  -u_ana * vg_ana(:,1) - v_ana * vg_ana(:,2)
            end if

            ! Coriolis terms
            !======================================
            absvo     = 2._r8 * OOmega * sin( lat_ana(jax) * PI/180._r8 )
            !Allow Coriolis to react to SCM winds
            uten_vort_ana  =   absvo * v_ana
            vten_vort_ana  =  -absvo * u_ana
            ! Force Coriolis to ALWAYS be calc w/ analysis winds
            !!uten_vort_ana  =   absvo * v_ana_diag
            !!vten_vort_ana  =  -absvo * u_ana_diag
            !  -----  Diags for VI form (0-out)
            uten_keg_ana  =    0._r8 ! fill with 0

          if (.FALSE.) then  ! No horz. p-gradient in p-coords
               uten_pfrc_ana  = - phig_ana(:,1)
               vten_pfrc_ana  = - phig_ana(:,2)
          else
              !put together pieces for Pressure and Phi gradient tencency terms
              uten_pfrc_ana  = -(1._r8/rho_ana) * plog_ana(:,1) - phig_ana(:,1)
              vten_pfrc_ana  = -(1._r8/rho_ana) * plog_ana(:,2) - phig_ana(:,2)
          end if


          if ( scm_ana_upwind .OR. scm_ana_t_react ) then
             tten_hadv_ana =  upwind_hadv(nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, u_ana, v_ana, Tv_X )
          else
             tten_hadv_ana  = -u_ana * teg_ana(:,1) - v_ana * teg_ana(:,2) ! should be straight T (not virtual)
          end if
          if ( scm_ana_upwind .OR. scm_ana_q_react ) then
             qten_hadv_ana =  upwind_hadv(nlon_alc,nlat_alc,nlev_alc,iac,jac,lon_alc,lat_alc, u_ana, v_ana, Qo_X )
          else
             qten_hadv_ana  = -u_ana * qg_ana(:,1) - v_ana * qg_ana(:,2)
          end if

         if (.not.(scm_ana_direct_omega)) then
            omega_ana = omega_recalc_ana ! use reconstructed omega
            if(masterproc) write(iulog,*) subname//" Omega recalc from ana U,V etc."
          else
            omega_ana = omega_dycore_ana ! use direct omega from dycore/ana
            if(masterproc) write(iulog,*) subname//" Omega direct from ana"
          end if


          if (.not.(scm_ana_x_plevels)) then
          !Tendencies due to vertical advection (etadot * D_eta ... )
            uten_vadv_ana = vadv_fv( nlev_alc, etad_ana, u_ana )
            vten_vadv_ana = vadv_fv( nlev_alc, etad_ana, v_ana )
            tten_vadv_ana = vadv_fv( nlev_alc, etad_ana, tv_ana ) ! should be straight T (not virtual)
            qten_vadv_ana = vadv_fv( nlev_alc, etad_ana, q_ana )
          else
          !Tendencies due to vertical advection (Omega * D_p ... )
            uten_vadv_ana = vadv_fv_press( nlev_alc, omega_ana, plo_ana, u_ana )
            vten_vadv_ana = vadv_fv_press( nlev_alc, omega_ana, plo_ana, v_ana )
            tten_vadv_ana = vadv_fv_press( nlev_alc, omega_ana, plo_ana, t_ana ) ! should be straight T (not virtual)
            qten_vadv_ana = vadv_fv_press( nlev_alc, omega_ana, plo_ana, q_ana )
          end if

            tten_comp_ana = (1./cpair)*( omega_ana / rho_ana )

          !DIags for pressure/geop grad forces
            uten_phig_ana =  - phig_ana(:,1)
            uten_prg_ana  =  - (1._r8/rho_ana) * plog_ana(:,1)

    end subroutine get_ana_dynfrc_fv

!-----------------------------------------------------
!    Stuff ... useful ojala
!-----------------------------------------------------
    !-------------------------
       function vadv_fv( nlev, etad, aa ) result( tend )
        use hycoef,           only: hyai, hybi, ps0, hyam, hybm
            integer,  intent(in)  :: nlev
            real(r8), intent(in)  :: etad(nlev) , aa(nlev)
            real(r8)              :: tend(nlev)
            real(r8)              :: eta(nlev)
            integer               :: L

            eta = hybm+hyam

            do L=2,nlev-1
               tend(L) = etad(L)* ( aa(L+1) - aa(L-1) ) / ( eta(L+1) - eta(L-1) )
            end do
            L=1
               tend(L) = etad(L)* ( aa(L+1) - aa(L) ) / ( eta(L+1) - eta(L) )
            L=nlev
               tend(L) = etad(L)* ( aa(L) - aa(L-1) ) / ( eta(L) - eta(L-1) )

            tend = -1._r8*tend ! for RHS consistency

         end function vadv_fv
!---------------------------
    !-------------------------
       function vadv_fv_press( nlev, omega, plo, aa ) result( tend )
            integer,  intent(in)  :: nlev
            real(r8), intent(in)  :: omega(nlev) , aa(nlev),plo(nlev)
            real(r8)              :: tend(nlev)
            integer               :: L

            do L=2,nlev-1
               tend(L) = omega(L)* ( aa(L+1) - aa(L-1) ) / ( plo(L+1) - plo(L-1) )
            end do
            L=1
               tend(L) = omega(L)* ( aa(L+1) - aa(L) ) / ( plo(L+1) - plo(L) )
            L=nlev
               tend(L) = omega(L)* ( aa(L) - aa(L-1) ) / ( plo(L) - plo(L-1) )

            tend = -1._r8*tend ! for RHS consistency

         end function vadv_fv_press
!---------------------------
      function lin_pfc_fv( nlon,nlat,nlev,iax,jax,lons,lats, pre, phi ) result( pfc )
         !use shr_kind_mod,  only:  r8 => shr_kind_r8
         !use shr_const_mod, only:  rearth => shr_const_rearth , &  !  =6.37122e6_R8 meters
         !                          pi => shr_const_pi         , &
         !                          omega => shr_const_omega

            integer,  intent(in)  :: nlon,nlat,nlev,iax,jax
            real(r8), intent(in)  :: pre(nlon,nlat,nlev+1),phi(nlon,nlat,nlev+1)
            real(r8), intent(in)  :: lats(nlat),lons(nlon)
            real(r8)              :: pfc(nlev,2)
            real(r8)              :: pfxW(nlev) , pfxE(nlev)
            real(r8)              :: pfyS(nlev) , pfyN(nlev)
            real(r8)              :: rlats(nlat),rlons(nlon),dx,dy,ds
            real(r8)              :: pr1,pr2,pr3,pr4, ph1,ph2,ph3,ph4
            integer               :: L , igg

         ! Begin
            rlons(:) = lons(:) * PI/180._r8
            rlats(:) = lats(:) * PI/180._r8

            dx=( rlons(2)-rlons(1) ) * Rearth
            dy=( rlats(2)-rlats(1) ) * Rearth

            ds  = MAX( dx*cos(rlats(jax)) , .1_r8 )
            igg = iax
            do L=1,nlev
               pr1 = pre(igg-1,jax,L+1)
               pr2 = pre(igg  ,jax,L+1)
               pr3 = pre(igg  ,jax,L  )
               pr4 = pre(igg-1,jax,L  )
               ph1 = phi(igg-1,jax,L+1)
               ph2 = phi(igg  ,jax,L+1)
               ph3 = phi(igg  ,jax,L  )
               ph4 = phi(igg-1,jax,L  )
               pfxW(L) = ( (pr2-pr4)*(ph1-ph3) + (pr1-pr3)*(ph4-ph2) ) /( ds * ( (pr2-pr4) + (pr1-pr3) ) )
            end do
            igg = iax +1
            do L=1,nlev
               pr1 = pre(igg-1,jax,L+1)
               pr2 = pre(igg  ,jax,L+1)
               pr3 = pre(igg  ,jax,L  )
               pr4 = pre(igg-1,jax,L  )
               ph1 = phi(igg-1,jax,L+1)
               ph2 = phi(igg  ,jax,L+1)
               ph3 = phi(igg  ,jax,L  )
               ph4 = phi(igg-1,jax,L  )
               pfxE(L) = ( (pr2-pr4)*(ph1-ph3) + (pr1-pr3)*(ph4-ph2) ) /( ds * ( (pr2-pr4) + (pr1-pr3) ) )
            end do
            ds  = dy
            igg = jax
            do L=1,nlev
               pr1 = pre(iax,igg-1,L+1)
               pr2 = pre(iax,igg  ,L+1)
               pr3 = pre(iax,igg  ,L  )
               pr4 = pre(iax,igg-1,L  )
               ph1 = phi(iax,igg-1,L+1)
               ph2 = phi(iax,igg  ,L+1)
               ph3 = phi(iax,igg  ,L  )
               ph4 = phi(iax,igg-1,L  )
               pfyS(L) = ( (pr2-pr4)*(ph1-ph3) + (pr1-pr3)*(ph4-ph2) ) /( ds * ( (pr2-pr4) + (pr1-pr3) ) )
            end do
            igg = jax +1
            do L=1,nlev
               pr1 = pre(iax,igg-1,L+1)
               pr2 = pre(iax,igg  ,L+1)
               pr3 = pre(iax,igg  ,L  )
               pr4 = pre(iax,igg-1,L  )
               ph1 = phi(iax,igg-1,L+1)
               ph2 = phi(iax,igg  ,L+1)
               ph3 = phi(iax,igg  ,L  )
               ph4 = phi(iax,igg-1,L  )
               pfyN(L) = ( (pr2-pr4)*(ph1-ph3) + (pr1-pr3)*(ph4-ph2) ) /( ds * ( (pr2-pr4) + (pr1-pr3) ) )
            end do


            do L=1,nlev
              pfc(L,1)  = 0.5_r8*( pfxW(L) + pfxE(L) )
              pfc(L,2)  = 0.5_r8*( pfyS(L) + pfyN(L) )
            end do



         end function lin_pfc_fv
   !-------------------------
       function grad_fv( nlon,nlat,nlev,iax,jax,lons,lats, aa ) result( ga )
         !use shr_kind_mod,  only:  r8 => shr_kind_r8
         !use shr_const_mod, only:  rearth => shr_const_rearth , &  !  =6.37122e6_R8 meters
         !                          pi => shr_const_pi         , &
         !                          omega => shr_const_omega

            integer,  intent(in)  :: nlon,nlat,nlev,iax,jax
            real(r8), intent(in)  :: aa(nlon,nlat,nlev)
            real(r8), intent(in)  :: lats(nlat),lons(nlon)
            real(r8)              :: ga(nlev,2)
            real(r8)              :: rlats(nlat),rlons(nlon),dx,dy
            integer               :: L

         ! Begin
            rlons(:) = lons(:) * PI/180._r8
            rlats(:) = lats(:) * PI/180._r8

            dx=( rlons(2)-rlons(1) ) * Rearth
            dy=( rlats(2)-rlats(1) ) * Rearth

            do L=1,nlev
              ga(L,1)  = (aa(iax+1,jax,L) - aa(iax-1,jax,L))/( 2._r8*dx*cos(rlats(jax)) + 0.1_r8 )
              ga(L,2)  = (aa(iax,jax+1,L) - aa(iax,jax-1,L))/( 2._r8*dy )
            end do



         end function grad_fv
   !-------------------------
       function upwind_hadv( nlon,nlat,nlev,iax,jax,lons,lats,u,v, aa ) result( hadv_tend )
         !use shr_kind_mod,  only:  r8 => shr_kind_r8
         !use shr_const_mod, only:  rearth => shr_const_rearth , &  !  =6.37122e6_R8 meters
         !                          pi => shr_const_pi         , &
         !                          omega => shr_const_omega

            integer,  intent(in)  :: nlon,nlat,nlev,iax,jax
            real(r8), intent(in)  :: aa(nlon,nlat,nlev)
            real(r8), intent(in)  :: lats(nlat),lons(nlon),u(nlev),v(nlev)
            real(r8)              :: hadv_tend(nlev)
            real(r8)              :: rlats(nlat),rlons(nlon),dx,dy,xten(nlev),yten(nlev)
            integer               :: L

         ! Begin
            rlons(:) = lons(:) * PI/180._r8
            rlats(:) = lats(:) * PI/180._r8

            dx=( rlons(2)-rlons(1) ) * Rearth
            dy=( rlats(2)-rlats(1) ) * Rearth

            do L=1,nlev
              if ( u(L) >= 0._r8 ) then
                 xten(L)  = u(L) * ( aa(iax,jax,L) - aa(iax-1,jax,L))/( dx*cos(rlats(jax)) + 0.1_r8 )
              else
                 xten(L)  = u(L) * ( aa(iax+1,jax,L) - aa(iax,jax,L))/( dx*cos(rlats(jax)) + 0.1_r8 )
              end if
            end do
            do L=1,nlev
              if ( v(L) >= 0._r8 ) then
                 yten(L)  = v(L) * ( aa(iax,jax,L) - aa(iax,jax-1,L))/( dy )
              else
                 yten(L)  = v(L) * ( aa(iax,jax+1,L) - aa(iax,jax,L))/( dy )
              end if
            end do

            hadv_tend(:) = -1._r8 * ( xten(:) + yten(:) )


         end function upwind_hadv
!=========================================
      subroutine makepk_fv( nlon,nlat,nlev, t, q, ps, phis, pko, pke, phi, th )
           use hycoef,           only: hyai, hybi, ps0, hyam, hybm
           integer,  intent(in)  :: nlon,nlat,nlev
            real(r8), intent(in)  :: t(nlon,nlat,nlev),q(nlon,nlat,nlev),ps(nlon,nlat),phis(nlon,nlat)
            real(r8), intent(out) :: pko(nlon,nlat,nlev),th(nlon,nlat,nlev),pke(nlon,nlat,nlev+1), phi(nlon,nlat,nlev+1)
            real(r8) :: ple(nlon,nlat,nlev+1),plo(nlon,nlat,nlev),rv(nlon,nlat,nlev)
            real(r8) :: kappa, p00
            integer :: L

           do L=1,nlev+1
              ple(:,:,L) = hyai(L)*ps0 + hybi(L)*ps(:,:)
           end do
           do L=1,nlev
              plo(:,:,L) = hyam(L)*ps0 + hybm(L)*ps(:,:)
           end do

           kappa=rdair/cpair

           pko = plo**kappa
           pke = ple**kappa

           p00 =    100000._r8
           th  = ( ( p00 / plo)**kappa ) * t

           rv = 1._r8/(1._r8 - q) - 1._r8
           th = th*(1._r8  + 0.61_r8 * rv )

           phi(:,:,nlev+1) = phis(:,:)
           do L=nlev,1,-1
             phi(:,:,L) = phi(:,:,L+1) - ( CpAir * Th(:,:,L) ) * (  pke(:,:,L) - pke(:,:,L+1) ) / (p00**kappa )
          end do


       end subroutine makepk_fv

!=============================================================================
      subroutine makepr_fv( nlon,nlat,nlev, t, ps, phis, plo, ple, phi )
           use hycoef,           only: hyai, hybi, ps0, hyam, hybm
           use shr_const_mod,    only: rdair => shr_const_rdair
           integer,  intent(in)  :: nlon,nlat,nlev
            real(r8), intent(in)  :: t(nlon,nlat,nlev),ps(nlon,nlat),phis(nlon,nlat)
            real(r8), intent(out) :: plo(nlon,nlat,nlev), ple(nlon,nlat,nlev+1), phi(nlon,nlat,nlev+1)
            real(r8) :: lnple(nlon,nlat,nlev+1)
            integer :: L

           do L=1,nlev+1
              ple(:,:,L) = hyai(L)*ps0 + hybi(L)*ps(:,:)
           end do
           do L=1,nlev
              plo(:,:,L) = hyam(L)*ps0 + hybm(L)*ps(:,:)
           end do

           lnple  = log( ple )
           phi(:,:,nlev+1) = phis(:,:)
           do L=nlev,1,-1
             phi(:,:,L) = phi(:,:,L+1) - (RdAir * T(:,:,L) ) * (  lnple(:,:,L) - lnple(:,:,L+1) )
             !phi(:,:,L) = phi(:,:,L+1) - (RdAir * T(:,:,L) / plo(:,:,L) ) * (  ple(:,:,L) - ple(:,:,L+1) )
          end do

      end subroutine makepr_fv

!=============================================================================
      subroutine virtual_t( nlon,nlat,nlev, t, q, tv )
           use hycoef,           only: hyai, hybi, ps0, hyam, hybm
           use shr_const_mod,    only: rdair => shr_const_rdair
           integer,  intent(in)  :: nlon,nlat,nlev
            real(r8), intent(in)  :: t(nlon,nlat,nlev),q(nlon,nlat,nlev)
            real(r8), intent(out) :: tv(nlon,nlat,nlev)
            real(r8)              :: rv(nlon,nlat,nlev)
            integer :: L


            rv = 1._r8/(1._r8 - q) - 1._r8
            tv = t*(1._r8  + 0.61_r8 * rv )


      end subroutine virtual_t

    !-------------------------
       subroutine zeta_fv( nlon,nlat,nlev,lons,lats, u,v, zeta )
         !use shr_kind_mod,  only:  r8 => shr_kind_r8
         !use shr_const_mod, only:  rearth => shr_const_rearth , &  !  =6.37122e6_R8 meters
         !                          pi => shr_const_pi         , &
         !                          omega => shr_const_omega

            integer,  intent(in)  :: nlon,nlat,nlev
            real(r8), intent(in)  :: u(nlon,nlat,nlev),v(nlon,nlat,nlev)
            real(r8), intent(out) :: zeta(nlev)
            real(r8), intent(in) :: lats(nlat),lons(nlon)
            real(r8) :: rlats(nlat),rlons(nlon)
            real(r8) :: dy,dx0,dx,darea,voo,voo2

            integer  :: iap,jap,iam,jam,i,j,L,iax,jax

            iax=2
            jax=2

            rlons(:) = lons(:) * PI/180._r8
            rlats(:) = lats(:) * PI/180._r8

            dx0 = rearth* ( rlons(2)-rlons(1) )
            dy  = rearth* ( rlats(2)-rlats(1) )

            darea = dy*dx0*cos( rlats(jax) )


             do L =1,nlev
               zeta(L) =                                                &
                           ( V(iax+1,jax, L)   - V(iax-1,jax,L) ) / ( 2._r8*dx0*cos( rlats(jax) ) )  &
                         - ( U(iax,jax+1, L)   - U(iax,jax-1,L) ) / ( 2._r8*dy )
            end do


  end subroutine zeta_fv
!================================================================
  subroutine etadot_fv ( nlon, nlat, nlev, lons, lats, u, v, plo, ple, etadot , omega )
  use shr_kind_mod,     only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use shr_const_mod,    only: rearth => shr_const_rearth , &  !  =6.37122e6_R8 meters
                               pi => shr_const_pi
  use hycoef,           only: hyai, hybi, ps0, hyam, hybm

  integer,  intent(in)  :: nlon,nlat,nlev
  real(r8), intent(in)  :: lons(nlon),lats(nlat)
  real(r8), intent(in)  :: u(nlon,nlat,nlev)  , v(nlon,nlat,nlev) , plo( nlon,nlat,nlev) , ple( nlon,nlat,nlev+1)
  real(r8), intent(out) :: etadot( nlon,nlat,nlev) ,omega(nlon,nlat,nlev)
  !real(r8), intent(in) :: uc(:,:,:)  , vc(:,:,:) , ple(:,:,:)

  ! Local variables
  real(r8),allocatable :: div(:,:,:)
  real(r8),allocatable :: mass(:,:,:), fuc(:,:,:),fvc(:,:,:)
  real(r8) :: rlats(nlat), rlons(nlon), rcos1, eta(nlev+1) , dx,dy! radians
  real(r8), allocatable :: etadot_t1(:,:), etadot_t2(:,:,:)
  integer :: i,j,L,im1,jm1,ip1,jp1
  real    :: uc_ijL , vc_ijL

   allocate (  div(nlon,nlat,nlev) )
   allocate (  mass(nlon,nlat,nlev), fuc(nlon,nlat,nlev),fvc(nlon,nlat+1,nlev) )
   allocate (  etadot_t1(nlon,nlat), etadot_t2(nlon,nlat,nlev) )

   div  = 0._r8
   fuc  = 0._r8
   fvc  = 0._r8
   mass = 0._r8
   etadot    = 0._r8
   etadot_t1 = 0._r8
   etadot_t2 = 0._r8

   rlons(:) = lons(:) * PI/180._r8
   rlats(:) = lats(:) * PI/180._r8

   do L=1,nlev+1
      eta(L) = hyai(L) + hybi(L)  ! 1._r8*L/(nlev+1)
   end do
   do L=1,nlev
      mass(:,:,L) = ( ple(:,:,L+1)-ple(:,:,L) )/( eta(L+1)-eta(L) )
   end do

   ! calculate mass fluxes at gridbox edges, using upwind algorithm
   do L=1,nlev
      do j=1,nlat
         do i=2,nlon
            im1=i-1
            !if ( i == 1) im1=nlon
            uc_ijL = 0.5_r8*( u(im1,j,L) + u(i,j,L) )
            if ( uc_ijL < 0._r8  ) fuc(i,j,L)= uc_ijL * mass(i,j,L)
            if ( uc_ijL >= 0._r8 ) fuc(i,j,L)= uc_ijL * mass(im1,j,L)
         end do
      end do
   end do
             ! Note: cos(lat) term incorporated into fluxes
   do L=1,nlev
      do j=2,nlat
         do i=1,nlon
            jm1=j-1
            vc_ijL = 0.5_r8 * ( v(i,jm1,L)+v(i,j,L) )
            if ( vc_ijL  < 0._r8  ) fvc(i,j,L)= vc_ijL * mass(i,j,L)    *cos( rlats(j) )
            if ( vc_ijL >= 0._r8  ) fvc(i,j,L)= vc_ijL * mass(i,jm1,L)  *cos( rlats(jm1) )
         end do
      end do
   end do


  ! now calculate HORZ divergence of (FUC,FVC). Note coslat term already
  ! incorporated in FVC.
  do L=1,nlev
      do j=1,nlat-1
         do i=1,nlon-1
            ip1=i+1
            jp1=j+1
            rcos1 = 1._r8 /( Rearth*cos( rlats(j) ) )
            div(i,j,L) = rcos1 * ( FUC(ip1,j,L)-FUC(i,j,L) ) / (rlons(ip1)-rlons(i) )  &
                      + rcos1 * ( FVC(i,jp1,L)-FVC(i,j,L) ) / (rlats(jp1)-rlats(j) )
         end do
      end do
   end do


   etadot_t1(:,:)=0._r8
   etadot_t2(:,:,:)=0._r8
   do L=1,nlev
      etadot_t1(:,:) = etadot_t1(:,:) + div(:,:,L)*(eta(L+1)-eta(L))
   end do
   do L=2,nlev
      etadot_t2(:,:,L) = etadot_t2(:,:,L-1) + div(:,:,L)*(eta(L+1)-eta(L))
   end do
   do L=1,nlev
      etadot(:,:,L) = ( hybm(L)*etadot_t1(:,:)  - etadot_t2(:,:,L) )  / mass(:,:,L)
   end do

  dx=( rlons(2)-rlons(1) ) * Rearth
  dy=( rlats(2)-rlats(1) ) * Rearth
  omega = 0._r8

#if 1
  do L=1,nlev
    do j=2,nlat-1
      do i=2,nlon-1
      omega(i,j,L)  = u(i,j,L) * (plo(i+1,j,L)-plo(i-1,j,L))/( 2._r8*dx*cos(rlats(j)) + 0.1_r8 ) &
                    + v(i,j,L) * (plo(i,j+1,L)-plo(i,j-1,L))/( 2._r8*dy ) &
                      - etadot_t2(i,j,L)
     end do
   end do
  end do
#else
  do L=1,nlev
    do j=2,nlat-1
      do i=2,nlon-1
      omega(i,j,L)  = etadot(i,j,L)*mass(i,j,L)
     end do
   end do
  end do
#endif


end subroutine etadot_fv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Reading netcdf files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!================================================================
  subroutine read_netcdf_ana_fv_ini( anal_file , nlon, nlat, nlev,lonidx,latidx )
   !
   ! READ_NETCDF_ANAL_INI:
   !                 Open the given analyses data file. Query dimesnisons.
   !                 Close.
   !===============================================================
   use cam_abortutils,   only : endrun
   use netcdf
   use shr_kind_mod,     only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use hycoef,           only: hyai, hybi, ps0, hyam, hybm
   use shr_const_mod,    only: rdair => shr_const_rdair
   use scammod,          only: scmlon,scmlat
   use shr_scam_mod,     only: shr_scam_getCloseLatLon  ! Standardized system subroutines

   !-------------
   character(len=*),intent(in):: anal_file

   integer, intent(out) ::  nlon,nlat,nlev,latidx,lonidx

   ! Local values
   !-------------
   integer :: ncid,varid,istat
   integer :: ilat,ilon,ilev
   integer :: i,j,L

   real(r8) :: closelon,closelat

   logical :: l_have_us , l_have_vs

   character(len=24) :: subname='read_netcdf_ana_fv_ini: '

   l_have_us = .FALSE.
   l_have_vs = .FALSE.

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlev)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     call shr_scam_getCloseLatLon(ncid ,scmlat,scmlon,closelat,closelon,latidx,lonidx)

     ! Close the analyses file and exit
     !-----------------------
     istat=nf90_close(ncid)
     if(istat /= NF90_NOERR) then
       write(iulog,*) subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

    endif ! (masterproc) then


 end subroutine read_netcdf_ana_fv_ini

!================================================================
  subroutine read_netcdf_ana_fv( anal_file , nlon, nlat, nlev, &
                                 u, v, &
                                 t, q, ps, phis, &
                                 lons, lats, levs &
                               , utcore, vtcore, ttcore, ogcore &
                                    )
   !
   ! READ_NETCDF_ANAL :
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values as well as Lons, Lats.
   !===============================================================
   use cam_abortutils,   only : endrun
   use netcdf
   use shr_kind_mod,     only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use hycoef,           only: hyai, hybi, ps0, hyam, hybm
   use shr_const_mod,    only: rdair => shr_const_rdair
   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   integer,  intent(in   ) ::  nlon,nlat,nlev
   real(r8), intent(out) ::  U(nlon,nlat,nlev), V(nlon,nlat,nlev)
   real(r8), intent(out) ::  T(nlon,nlat,nlev), Q(nlon,nlat,nlev)
   real(r8), intent(out) ::  PS(nlon,nlat), PHIS(nlon,nlat)
   !real(r8), intent(out) ::  PHI(nlon,nlat,nlev+1),PLE(nlon,nlat,nlev+1),PLO(nlon,nlat,nlev)
   real(r8), intent(out) ::  Lats(nlat),Lons(nlon),Levs(nlev)

   real(r8), intent(out) ::  UTCORE(nlon,nlat,nlev), VTCORE(nlon,nlat,nlev), TTCORE(nlon,nlat,nlev)
   real(r8), intent(out) ::  OGCORE(nlon,nlat,nlev)

   ! Local values
   !-------------
   integer :: ncid,varid,istat
   integer :: ilat,ilon,ilev
   integer :: i,j,L

   logical :: l_have_us , l_have_vs

   character(len=20) :: subname='read_netcdf_ana_fv: '

   l_have_us = .FALSE.
   l_have_vs = .FALSE.

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) subname//'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
    end if
   end if



   if(masterproc) then

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lons)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lats)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lev',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Levs)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif


     ! Read in, transpose lat/lev indices,
     ! and scatter data arrays
     !----------------------------------
     !  First block reads U
     !----------------------------------
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid, U )
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid, V )
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid, T )
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid, Q )
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,PS )
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif

    istat=nf90_inq_varid(ncid,'PHIS',varid)
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PHIS )
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif

    istat=nf90_inq_varid(ncid,'UTEND_CORE',varid)
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//"No UTEND_CORE on file: "
      write(iulog,*) trim(anal_file)
      utcore(:,:,:)=-9999._r8
    else
      istat=nf90_get_var(ncid,varid,utcore )
        if(istat /= NF90_NOERR) then
          write(iulog,*)  subname//nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
       end if
    end if

    istat=nf90_inq_varid(ncid,'VTEND_CORE',varid)
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//"No VTEND_CORE on file: "
      write(iulog,*) trim(anal_file)
      vtcore(:,:,:)=-9999._r8
    else
      istat=nf90_get_var(ncid,varid,vtcore )
        if(istat /= NF90_NOERR) then
          write(iulog,*)  subname//nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
        end if
    end if

    istat=nf90_inq_varid(ncid,'DTCORE',varid)
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//"No TTEND_CORE on file: "
      write(iulog,*) trim(anal_file)
      ttcore(:,:,:)=-9999._r8
    else
      istat=nf90_get_var(ncid,varid,ttcore )
        if(istat /= NF90_NOERR) then
          write(iulog,*)  subname//nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
        end if
    end if

    istat=nf90_inq_varid(ncid,'OMEGA',varid)
    if(istat /= NF90_NOERR) then
      write(iulog,*)  subname//"No OMEGA (core) on file: "
      write(iulog,*) trim(anal_file)
      ogcore(:,:,:)=-9999._r8
    else
      istat=nf90_get_var(ncid,varid,ogcore )
        if(istat /= NF90_NOERR) then
          write(iulog,*)  subname//nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
        end if
    end if


     ! Close the analysis file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat /= NF90_NOERR) then
       write(iulog,*)  subname//nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

   end if ! (masterproc) then
   !------------

   return
  end subroutine read_netcdf_ana_fv
!================================================================
!================================================================
  subroutine dynfrc_timewgts (  &
                                ana_prev_date, ana_next_date, &
                                wgt1 , wgt2 )


  use shr_kind_mod,     only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use ESMF
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size

  integer,  intent(in)  :: ana_prev_date(4), ana_next_date(4)
  real(r8) , intent(out) :: wgt1,wgt2

   type(ESMF_Time)         :: Date1,Date2,Date0
   type(ESMF_TimeInterval) :: DateDiff2,DateDiff0,DateDiff, AnaDiff
   integer                 :: DeltaT0, DeltaT2 , YMD, Year,Month,Day,Sec, Ana_interval, rc

   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

       call ESMF_TimeSet(Date0,YY=Ana_prev_date(1), MM=Ana_prev_date(2) , &
                               DD= Ana_prev_date(3)  , S= Ana_prev_date(4)  )
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)

       call ESMF_TimeSet(Date2,YY=Ana_next_date(1), MM=Ana_next_date(2) , &
                               DD= Ana_next_date(3)  , S= Ana_next_date(4)  )
       AnaDiff =Date2-Date0
       call ESMF_TimeIntervalGet(AnaDiff,S=Ana_interval ,rc=rc)

       DateDiff2 =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff2,S=DeltaT2,rc=rc)
       DateDiff0 =Date1-Date0
       call ESMF_TimeIntervalGet(DateDiff0,S=DeltaT0,rc=rc)

            wgt1 = 1._r8 - ( 1._r8 * DeltaT0 ) / Ana_interval
            wgt2 = 1._r8 - ( 1._r8 * DeltaT2 ) / Ana_interval

end subroutine dynfrc_timewgts


!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine patch_eta_x_plv (  nx , ny, nL,ix, jx, aa, plo )
   integer, intent(in)     :: nx,ny,nl,ix,jx
   real(r8), intent(in)    :: plo(nx,ny,nL)
   real(r8), intent(inout) :: aa(nx,ny,nL)

   real(r8) :: plx(nL),plq(nL),aax(nL),aaq(nL),aat(nx,ny,nL)
   real(r8) :: dp,dpk,dpk1,wtk,wtk1
   integer  :: i,j,L,k


   plx(:) = plo(ix,jx,:) ! target pressures

   do j=1,ny
   do i=1,nx
      plq(:) =  plo(i,j,:)
      aaq(:) =  aa(i,j,:)
      do L=1,nl
      do k=2,nl
         if ( ( plx(L) <= plq(k) ).AND.(plx(L) > plq(k-1) ) ) then
            dp  = plq(k)-plq(k-1)
            dpk1 = plx(L)-plq(k-1)
            dpk  = plq(k)-plx(L)
            wtk1 = 1._r8 - dpk1 / dp
            wtk  = 1._r8 - dpk  / dp
            aax(L) = wtk * aaq(k) + wtk1 * aaq(k-1)
         end if
      end do
      if (  plx(L) <= plq(1)  ) aax(L)=aaq(1)
      if (  plx(L) >  plq(NL) ) aax(L)=aaq(NL)
      end do

      aat(i,j,:)=aax(:)
   end do
   end do

   aa=aat

  end subroutine patch_eta_x_plv


end module get_ana_dynfrc_4scam
