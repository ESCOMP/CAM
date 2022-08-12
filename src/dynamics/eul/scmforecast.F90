#define SCAMNUDGERUN
module scmforecast
   ! --------------------------------------------------------------------------- ! 
   !                                                                             !
   ! Compute Time-Marched 'T, u, v, q' for SCAM by summing the 'physics',        !
   ! 'horizontal advection', and 'vertical advection' tendencies.                ! 
   ! This module is used only for SCAM.                                          ! 
   !                                                                             ! 
   ! --------------------------------------------------------------------------- ! 
  use spmd_utils,       only: masterproc
  use cam_logfile,      only: iulog
  use cam_control_mod,  only: adiabatic
!++jtb
#ifdef SCAMNUDGERUN
  use get_ana_dynfrc_4scam,    only: get_ana_dynfrc_fv
#endif
!--jtb
  implicit none
  private
  save

  public forecast
!
! Private module data
!

!======================================================================= 
contains
!======================================================================= 


   subroutine forecast( lat       , nlon      , ztodt     ,           & 
                        psm1      , psm2      , ps        ,           &
                        u3        , u3m1      , u3m2      ,           &
                        v3        , v3m1      , v3m2      ,           &
                        t3        , t3m1      , t3m2      ,           &
                        q3        , q3m1      , q3m2      ,           & 
                        tten_phys , uten_phys , vten_phys ,           &
                        qminus    , qfcst     )

   ! --------------------------------------------------------------------------- ! 
   !                                                                             !
   ! Compute Time-Marched 'T, u, v, q' for SCAM by summing the 'physics',        !
   ! 'horizontal advection', and 'vertical advection' tendencies.                ! 
   ! This module is used only for SCAM.                                          ! 
   !                                                                             ! 
   ! Author : Sungsu Park. 2010. Sep.                                            !
   !                                                                             !
   ! --------------------------------------------------------------------------- !

     use shr_kind_mod,    only : r8 => shr_kind_r8, i8 => shr_kind_i8
     use pmgrid,          only : plev, plat, plevp, plon
     use cam_history,     only : outfld
     use constituents,    only : pcnst, cnst_get_ind, cnst_name
     use physconst,       only : rair, cpair, gravit, rga
     use scammod,         only : divq,divq3d,divt,divu,divt3d,divu3d,have_divv, &
                                 divv,divv3d,have_aldif,have_aldir,have_asdif,have_asdir, &
                                 have_cld,have_cldice,have_cldliq,have_clwp,have_divq,have_divq3d, &
                                 have_divt,have_divt3d,have_divu,have_divu3d,have_divv3d,have_numice, &
                                 have_numliq,have_omega,have_phis,have_prec,have_ps,have_ptend, &
                                 have_q,have_q1,have_q2,have_t,have_u,have_v, &
                                 have_vertdivq,have_vertdivt,have_vertdivu,have_vertdivv,qdiff,qobs, &
                                 scm_relax_bot_p,scm_relax_linear,scm_relax_tau_bot_sec, &
                                 scm_relax_tau_sec,scm_relax_tau_top_sec,scm_relax_top_p, &
                                 scm_relaxation,scm_use_obs_qv,scm_use_obs_t,scm_use_obs_uv,scm_zadv_q,scm_zadv_t, &
                                 scm_zadv_uv,tdiff,tobs,uobs,use_3dfrc,use_camiop,vertdivq, &
                                 vertdivt,vertdivu,vertdivv,vobs,wfld,qinitobs,scm_relax_fincl, &
!++jtb
                                 scmlon,scmlat, & 
                                 scm_ana_direct_ttend, &
                                 scm_use_ana_iop
!--jtb
     use time_manager,     only : get_curr_calday, get_nstep, get_step_size, is_first_step
     use cam_abortutils,   only : endrun
     use string_utils,     only: to_upper
!++jtb
     use shr_const_mod,    only:  rearth => shr_const_rearth , &  !  =6.37122e6_R8 meters
                               pi => shr_const_pi         , & 
                               OOmega => shr_const_omega  
!--jtb

     implicit none

   ! ---------------------- !
   ! Parameters             !
   ! ---------------------- !

    character(len=*), parameter ::  subname = "forecast"
    real(r8),parameter :: hugebad=9.99e12_r8

   ! --------------------------------------------------- !
   ! x = t, u, v, q                                      !
   ! x3m1 : state variable used for computing 'forcing'  !
   ! x3m2 : initial state variable before time-marching  !
   ! x3   : final   state variable after  time-marching  !
   ! --------------------------------------------------- !

     integer,  intent(in)    :: lat                     
     integer,  intent(in)    :: nlon                     
     real(r8), intent(in)    :: ztodt                   ! Twice time step unless nstep = 0 [ s ]

     real(r8), intent(inout) :: ps(plon)                ! Surface pressure [ Pa ]
     real(r8), intent(inout) :: psm1(plon)              ! Surface pressure [ Pa ]
     real(r8), intent(inout) :: psm2(plon)              ! Surface pressure [ Pa ]

     real(r8), intent(in)    :: t3m1(plev)              ! Temperature [ K ]
     real(r8), intent(inout)    :: t3m2(plev)              ! Temperature [ K ]
     real(r8), intent(in)    :: u3m1(plev)              ! Zonal wind [ m/s ]
     real(r8), intent(inout)    :: u3m2(plev)              ! Zonal wind [ m/s ]
     real(r8), intent(in)    :: v3m1(plev)              ! Meridional wind [ m/s ]
     real(r8), intent(inout)    :: v3m2(plev)              ! Meridional wind [ m/s ]
     real(r8), intent(inout) :: q3m1(plev,pcnst)        ! Tracers [ kg/kg, #/kg ]
     real(r8), intent(inout) :: q3m2(plev,pcnst)        ! Tracers [ kg/kg, #/kg ]

     real(r8), intent(inout) :: tten_phys(plev)         ! Tendency of T by the 'physics' [ K/s ]
     real(r8), intent(inout) :: uten_phys(plev)         ! Tendency of u by the sum of 'physics + geostrophic forcing' [ m/s/s ]
     real(r8), intent(inout) :: vten_phys(plev)         ! Tendency of v by the sum of 'physics + geostrophic forcing' [ m/s/s ]
     real(r8)                   qten_phys(plev,pcnst)   ! Tendency of q by the 'physics' [ #/kg/s, kg/kg/s ]
     real(r8), intent(in)    :: qminus(plon,plev,pcnst) ! ( qminus - q3m2 ) / ztodt = Tendency of tracers by the 'physics' [ #/kg/s, kg/kg/s ]  

     real(r8), intent(out)   :: t3(plev)                ! Temperature [ K ]
     real(r8), intent(out)   :: u3(plev)                ! Zonal wind [ m/s ]
     real(r8), intent(out)   :: v3(plev)                ! Meridional wind [ m/s ]
     real(r8), intent(inout) :: q3(plev,pcnst)          ! Tracers [ #/kg, kg/kg ]
     real(r8), intent(inout) :: qfcst(plon,plev,pcnst)  ! ( Input qfcst - q3m2 ) / ztodt = Tendency of q by the sum of 'physics' + 'SLT vertical advection' [ #/kg/s, kg/kg/s ]


   ! --------------- !
   ! Local Variables !
   ! --------------- !

     integer  dummy
     integer  dummy_dyndecomp
     integer  i, k, m              
     integer  ixcldliq, ixcldice, ixnumliq, ixnumice
     real(r8) weight, fac
     real(r8) pmidm1(plev)       
     real(r8) pintm1(plevp)      
     real(r8) pdelm1(plev)       
     real(r8) wfldint(plevp)     
     real(r8) pdelb(plon,plev)   
     real(r8) tfcst(plev)             ! (       tfcst - t3m2 ) / ztodt = Tendency of T by the sum of 'physics' + 'SLT/EUL/XXX vertical advection' [ K/s ]
     real(r8) ufcst(plev)             ! (       ufcst - u3m2 ) / ztodt = Tendency of u by the sum of 'physics' + 'SLT/EUL/XXX vertical advection' [ m/s/s ] 
     real(r8) vfcst(plev)             ! (       vfcst - u3m2 ) / ztodt = Tendency of v by the sum of 'physics' + 'SLT/EUL/XXX vertical advection' [ m/s/s ]
     logical scm_fincl_empty
   ! ----------------------------------------------- !
   ! Centered Eulerian vertical advective tendencies !
   ! ----------------------------------------------- !

     real(r8) tten_zadv_EULc(plev)              ! Vertical advective forcing of t [ K/s ]
     real(r8) uten_zadv_EULc(plev)              ! Vertical advective forcing of u [ m/s/s ] 
     real(r8) vten_zadv_EULc(plev)              ! Vertical advective forcing of v [ m/s/s ] 
     real(r8) qten_zadv_EULc(plev,pcnst)        ! Vertical advective forcing of tracers [ #/kg/s, kg/kg/s ]

   ! --------------------------------- !
   ! SLT vertical advective tendencies !
   ! --------------------------------- !
     real(r8) qten_zadv_SLT(plev,pcnst)         ! Vertical advective forcing of tracers [ #/kg/s, kg/kg/s ]

   ! ---------------------------- !
   ! Eulerian compression heating !
   ! ---------------------------- !

     real(r8) tten_comp_EUL(plev)               ! Compression heating by vertical advection [ K/s ]    
 
   ! ----------------------------------- !
   ! Final vertical advective tendencies !
   ! ----------------------------------- !  

     real(r8) tten_zadv(plev)                   ! Vertical advective forcing of t [ K/s ]
     real(r8) uten_zadv(plev)                   ! Vertical advective forcing of u [ m/s/s ] 
     real(r8) vten_zadv(plev)                   ! Vertical advective forcing of v [ m/s/s ] 
     real(r8) qten_zadv(plev,pcnst)             ! Vertical advective forcing of tracers [ #/kg/s, kg/kg/s ]


   ! --------------------------- !
   ! For 'scm_relaxation' switch !
   ! --------------------------- !

     real(r8) rtau(plev)
     real(r8) relax_T(plev)
     real(r8) relax_u(plev)
     real(r8) relax_v(plev)
     real(r8) relax_q(plev,pcnst)
  ! +++BPM: allow linear relaxation profile
     real(r8) rslope ! [optional] slope for linear relaxation profile
     real(r8) rycept ! [optional] y-intercept for linear relaxtion profile


!++jtb
! ------------------------------------ !
! Quantities derived from Analyses     !
! ------------------------------------ !
!======================================! 
     real(r8) dynfrcp(plev)                     ! Scaling factor for ana-derived tends
     logical  l_vectinv
     real(r8) omega_ana(plev)                   ! Vertical pressure velocity [ Pa/s ]
     real(r8) etad_ana(plev)                    ! "Eta dot" velocity [ Pa/s ]
     real(r8) T_ana(plev),  Q_ana(plev)  , Tv_ana(plev)       ! 
     real(r8) u_ana(plev),  v_ana(plev)         ! 
     real(r8) zeta_ana(plev)                    ! 
     real(r8) ps_ana 
     real(r8) T_ana_diag(plev),  Q_ana_diag(plev)        ! 
     real(r8) u_ana_diag(plev),  v_ana_diag(plev)         ! 
     ! ----------------------------------- !
     ! vertical advective tendencies !
     ! ----------------------------------- !  
     real(r8) tten_vadv_ana(plev)                   ! Vertical advective forcing of t [ K/s ]
     real(r8) uten_vadv_ana(plev)                   ! Vertical advective forcing of u [ m/s/s ] 
     real(r8) vten_vadv_ana(plev)                   ! Vertical advective forcing of v [ m/s/s ] 
     real(r8) qten_vadv_ana(plev)                   ! Vertical advective forcing of tracers [ #/kg/s, kg/kg/s ]
     ! ------------------------------------- !
     ! Horizontal advective/other tendencies !
     ! ------------------------------------- !  
     real(r8) uten_hadv_ana(plev)                    ! of u [ m/s/s ] 
     real(r8) vten_hadv_ana(plev)                    ! of v [ m/s/s ] 
     real(r8) uten_pfrc_ana(plev)                    ! of u [ m/s/s ] 
     real(r8) vten_pfrc_ana(plev)                    ! of v [ m/s/s ] 
     real(r8) uten_vort_ana(plev)                    ! of u [ m/s/s ] 
     real(r8) vten_vort_ana(plev)                    ! of v [ m/s/s ] 
     real(r8) tten_hadv_ana(plev)                   ! of t [ K/s ]
     real(r8) qten_hadv_ana(plev)             ! of tracers [ #/kg/s, kg/kg/s ]

     !---------------------------------!
     ! Adiabatic compression tendency  !
     !---------------------------------!
     real(r8) tten_comp_ana(plev)                   ! of t [ K/s ]


     real(r8) uten_keg_ana(plev)                     ! of u [ m/s/s ] 
     real(r8) uten_prg_ana(plev)                     ! of u [ m/s/s ] 
     real(r8) uten_phig_ana(plev)                    ! of u [ m/s/s ] 
     ! ------------------------------------------ !
     ! Direct dycore or ana tendencies or quants  !
     ! Not recalculated.                          !
     ! (not usually available,                    !
     !  set=-9999 if missing )                    !
     ! ------------------------------------------ !  
     real(r8) tten_dycore_ana(plev)                   ! Total direct Ana forcing of t [ K/s ]
     real(r8) vten_dycore_ana(plev)                   ! Total direct Ana forcing of v [ m/s/s ]
     real(r8) uten_dycore_ana(plev)                   ! Total direct Ana forcing of u [ m/s/s ] 
     real(r8) omega_dycore_ana(plev)                  ! Omega direct from Ana/dycore (not recalc) [ Pa/s ] 
     ! ----------------------------------- !
     ! total recalc. "dycore" tendencies !
     ! ----------------------------------- !  
     real(r8) omega_recalc_ana(plev)                  ! Omega from Ana/dycore (recalculated) [ Pa/s ] 
     real(r8) tten_totdyn_ana(plev)                   ! Total Ana forcing of t [ K/s ]
     real(r8) uten_totdyn_ana(plev)                   ! Total Ana forcing of u [ m/s/s ] 
     real(r8) vten_totdyn_ana(plev)                   ! Total Ana forcing of v [ m/s/s ] 
     real(r8) qten_totdyn_ana(plev)                   ! Total Ana forcing of tracers [ #/kg/s, kg/kg/s ]
     real(r8) fcoriol,uten_coriol(plev),vten_coriol(plev)      
     real(r8) ufcstm2(plev),vfcstm2(plev)
     real(r8) ufcor_0(plev),vfcor_0(plev)
     real(r8) uten_totdyn_anax(plev)                   ! Total Ana forcing of u [ m/s/s ] 
     real(r8) vten_totdyn_anax(plev)                   ! Total Ana forcing of v [ m/s/s ] 
     real(r8) tfw0, tfw1, tfw2, tftotw,ztodtn,AA
     integer  nsubdyn,nt,nstep_curr

!+++ARH
     !logical use_ana_iop
!---ARH
     logical l_use_reconst_ttend  ! use reconstructed T-tendency based on analysis
     logical l_use_direct_ttend   ! use T-tendency direct from dycore


      l_use_reconst_ttend = .NOT.( scm_ana_direct_ttend )
      l_use_direct_ttend  = .NOT.( l_use_reconst_ttend )

!--jtb

!+++ BPM check what we have:
     if (masterproc .and. is_first_step()) write(iulog,*) 'SCAM FORECAST REPORT: ' ,  &
          'have_divq     ',  have_divq      ,  &
          'have_divt     ',  have_divt      ,  &
          'have_divq3d   ',  have_divq3d    ,  &
          'have_vertdivt ',  have_vertdivt  ,  &
          'have_vertdivu ',  have_vertdivu  ,  &
          'have_vertdivv ',  have_vertdivv  ,  &
          'have_vertdivq ',  have_vertdivq  ,  &
          'have_divt3d   ',  have_divt3d    ,  &
          'have_divu3d   ',  have_divu3d    ,  &
          'have_divv3d   ',  have_divv3d    ,  &
          'have_divu     ',  have_divu      ,  &
          'have_divv     ',  have_divv      ,  &
          'have_omega    ',  have_omega     ,  &
          'have_phis     ',  have_phis      ,  &
          'have_ptend    ',  have_ptend     ,  &
          'have_ps       ',  have_ps        ,  &
          'have_q        ',  have_q         ,  &
          'have_q1       ',  have_q1        ,  &
          'have_q2       ',  have_q2        ,  &
          'have_prec     ',  have_prec      ,  &
          'have_t        ',  have_t         ,  &
          'have_u        ',  have_u         ,  &
          'have_v        ',  have_v         ,  &
          'have_cld      ',  have_cld       ,  &
          'have_cldliq   ',  have_cldliq    ,  &
          'have_cldice   ',  have_cldice    ,  &
          'have_numliq   ',  have_numliq    ,  &
          'have_numice   ',  have_numice    ,  &
          'have_clwp     ',  have_clwp      ,  &
          'have_aldir    ',  have_aldir     ,  &
          'have_aldif    ',  have_aldif     ,  &
          'have_asdir    ',  have_asdir     ,  &
          'have_asdif    ',  have_asdif     ,  &
          'use_camiop    ',  use_camiop     ,  &
          'use_obs_uv    ',  scm_use_obs_uv     ,  &
          'use_obs_qv    ',  scm_use_obs_qv     ,  &
          'use_obs_T     ',  scm_use_obs_T      ,  &
          'relaxation    ',  scm_relaxation     ,  &
          'use_3dfrc     ',  use_3dfrc
     
     !---BPM


   ! ---------------------------- !
   !                              ! 
   ! Main Computation Begins Here !
   !                              !
   ! ---------------------------- !

     dummy = 2
     dummy_dyndecomp = 1


   ! ------------------------------------------------------------ !
   ! Calculate midpoint pressure levels                           !
   ! ------------------------------------------------------------ !
     call plevs0( nlon, plon, plev, psm1, pintm1, pmidm1, pdelm1 )

     call cnst_get_ind( 'CLDLIQ', ixcldliq, abort=.false. )
     call cnst_get_ind( 'CLDICE', ixcldice, abort=.false. )
     call cnst_get_ind( 'NUMLIQ', ixnumliq, abort=.false. )
     call cnst_get_ind( 'NUMICE', ixnumice, abort=.false. )

   ! ------------------------------------------------------------ !
   ! Extract physical tendencies of tracers q.                    !
   ! Note 'tten_phys, uten_phys, vten_phys' are already input.    !
   ! ------------------------------------------------------------ !

     qten_phys(:plev,:pcnst)     = ( qminus(1,:plev,:pcnst) - q3m2(:plev,:pcnst) ) / ztodt  

   ! ----------------------------------------------------- !
   ! Extract SLT-transported vertical advective tendencies !
   ! TODO : Add in SLT transport of t u v as well            !
   ! ----------------------------------------------------- !

     qten_zadv_SLT(:plev,:pcnst) = ( qfcst(1,:plev,:pcnst) - qminus(1,:plev,:pcnst) ) / ztodt  

   ! ------------------------------------------------------- !
   ! use_camiop = .true.  : Use  CAM-generated 3D   IOP file ! 
   !            = .false. : Use User-generated SCAM IOP file ! 
   ! ------------------------------------------------------- !  

#ifdef SCAMNUDGERUN
    !!! use_ana_iop needs to get into namelist!!  !!!!
!+++ARH
      !use_ana_iop=.TRUE.
      !!use_ana_iop=.FALSE.
!---ARH
      l_vectinv =.FALSE.

!+++ARH
      !if (use_ana_iop) then
      if (scm_use_ana_iop) then
!---ARH
      call get_ana_dynfrc_fv ( scmlon, scmlat ,  &
                               omega_ana, etad_ana, zeta_ana, &
                               t_ana ,  tv_ana , &
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
                               u3m2, v3m2, t3m2, q3m2(:,1), &
                               u_ana_diag, v_ana_diag, t_ana_diag, q_ana_diag   )
          else
                          ! set these to a "bad" value
                               omega_ana = HugeBad 
                               etad_ana = HugeBad 
                               zeta_ana = HugeBad 
                               t_ana = HugeBad  
                               tv_ana = HugeBad 
                               q_ana = HugeBad           
                               u_ana = HugeBad           
                               v_ana = HugeBad           
                               t_ana_diag = HugeBad  
                               q_ana_diag = HugeBad           
                               u_ana_diag = HugeBad           
                               v_ana_diag = HugeBad           
                               ps_ana = HugeBad           
                               uten_hadv_ana = HugeBad   
                               vten_hadv_ana = HugeBad   
                               uten_pfrc_ana = HugeBad   
                               vten_pfrc_ana = HugeBad   
                               uten_vort_ana = HugeBad   
                               vten_vort_ana = HugeBad   
                               qten_hadv_ana = HugeBad   
                               tten_hadv_ana = HugeBad   
                               uten_vadv_ana = HugeBad   
                               vten_vadv_ana = HugeBad   
                               tten_vadv_ana = HugeBad   
                               qten_vadv_ana = HugeBad   
                               tten_comp_ana = HugeBad   
                               uten_keg_ana   = HugeBad
                               uten_phig_ana = HugeBad    
                               uten_prg_ana = HugeBad
                               uten_dycore_ana = HugeBad 
                               vten_dycore_ana = HugeBad    
                               tten_dycore_ana = HugeBad    
                               omega_dycore_ana = HugeBad    
                               omega_recalc_ana = HugeBad    
          endif
 
   ! -------------------------------------------------------------- !
   ! Re-Calculate midpoint pressure levels if PS_ANA is reasonable  !
   ! -------------------------------------------------------------- !
     if (ps_ana < 500000._r8 ) then
        psm1=ps_ana
        call plevs0( nlon, plon, plev, psm1, pintm1, pmidm1, pdelm1 )
     end if
          if(l_vectinv) then
          uten_totdyn_ana  = uten_hadv_ana + uten_pfrc_ana + uten_vadv_ana
          vten_totdyn_ana  = vten_hadv_ana + vten_pfrc_ana + vten_vadv_ana
          uten_totdyn_anax = uten_hadv_ana + uten_pfrc_ana + uten_vadv_ana
          vten_totdyn_anax = vten_hadv_ana + vten_pfrc_ana + vten_vadv_ana
          else
          uten_totdyn_ana  = uten_hadv_ana + uten_vort_ana + uten_pfrc_ana + uten_vadv_ana
          vten_totdyn_ana  = vten_hadv_ana + vten_vort_ana + vten_pfrc_ana + vten_vadv_ana
          uten_totdyn_anax = uten_hadv_ana + uten_vort_ana + uten_pfrc_ana + uten_vadv_ana
          vten_totdyn_anax = vten_hadv_ana + vten_vort_ana + vten_pfrc_ana + vten_vadv_ana
          endif

          tten_totdyn_ana = tten_hadv_ana + tten_vadv_ana + tten_comp_ana
          qten_totdyn_ana = qten_hadv_ana + qten_vadv_ana
#else
!+++ARH
      !use_ana_iop=.FALSE.
!---ARH
#endif

!++jtb
   ! Need 3rd option 'use_ana_iop'
   !    - suboption: use {u,v,t,q}ten_vadv_ana OR recalculate with etad_ana
   !    - what about other species in q?
   !    - we might want to calculate fu,fv using evolving (local) u's and v's
   !      to allow geostrophic adjustment.
!--jtb

if( use_camiop ) then
       do k = 1, plev
          tfcst(k) = t3m2(k) + ztodt * tten_phys(k) + ztodt * divt3d(k)
          ufcst(k) = u3m2(k) + ztodt * uten_phys(k) + ztodt * divu3d(k)
          vfcst(k) = v3m2(k) + ztodt * vten_phys(k) + ztodt * divv3d(k)
          do m = 1, pcnst
           ! Below two lines are identical but in order to reproduce the bit-by-bit results 
           ! of CAM-3D simulation, I simply rewrite the 'original' into the 'expanded' one.
           ! Below is the 'original' one.
           ! qfcst(1,k,m) = q3m2(k,m) + ztodt * ( qten_phys(k,m) + divq3d(k,m) )
           ! Below is the 'expanded' one.
             qfcst(1,k,m) = qminus(1,k,m) + ztodt * divq3d(k,m)
          enddo
       enddo

else  ! when use_camiop =.FALSE.
!+++ARH
     !if( .NOT.(use_ana_iop) ) then
     if( .NOT.(scm_use_ana_iop) ) then
!---ARH
     ! ---------------------------------------------------------------------------- !
     ! Compute 'omega'( wfldint ) at the interface from the value at the mid-point. ! 
     ! SCAM-IOP file must provide omega at the mid-point not at the interface.      !
     ! ---------------------------------------------------------------------------- !
 
       wfldint(1) = 0._r8
       do k = 2, plev
          weight = ( pintm1(k) - pmidm1(k-1) ) / ( pmidm1(k) - pmidm1(k-1) )
          wfldint(k) = ( 1._r8 - weight ) * wfld(k-1) + weight * wfld(k)
       enddo
       wfldint(plevp) = 0._r8
   
     ! ------------------------------------------------------------ ! 
     ! Compute Eulerian compression heating due to vertical motion. !
     ! ------------------------------------------------------------ !

       do k = 1, plev
          tten_comp_EUL(k) = wfld(k) * t3m1(k) * rair / ( cpair * pmidm1(k) )
       enddo

     ! ---------------------------------------------------------------------------- !
     ! Compute Centered Eulerian vertical advective tendencies for all 't, u, v, q' ! 
     ! ---------------------------------------------------------------------------- ! 

       do k = 2, plev - 1
          fac = 1._r8 / ( 2.0_r8 * pdelm1(k) )
          tten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( t3m1(k+1) - t3m1(k) ) + wfldint(k) * ( t3m1(k) - t3m1(k-1) ) )
          vten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( v3m1(k+1) - v3m1(k) ) + wfldint(k) * ( v3m1(k) - v3m1(k-1) ) ) 
          uten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( u3m1(k+1) - u3m1(k) ) + wfldint(k) * ( u3m1(k) - u3m1(k-1) ) )
          do m = 1, pcnst
             qten_zadv_EULc(k,m) = -fac * ( wfldint(k+1) * ( q3m1(k+1,m) - q3m1(k,m) ) + wfldint(k) * ( q3m1(k,m) - q3m1(k-1,m) ) )
          end do
       end do

       k = 1
       fac = 1._r8 / ( 2.0_r8 * pdelm1(k) )
       tten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( t3m1(k+1) - t3m1(k) ) )
       vten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( v3m1(k+1) - v3m1(k) ) )
       uten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( u3m1(k+1) - u3m1(k) ) )
       do m = 1, pcnst
          qten_zadv_EULc(k,m) = -fac * ( wfldint(k+1) * ( q3m1(k+1,m) - q3m1(k,m) ) )
       end do

       k = plev
       fac = 1._r8 / ( 2.0_r8 * pdelm1(k) )
       tten_zadv_EULc(k) = -fac * ( wfldint(k) * ( t3m1(k) - t3m1(k-1) ) )
       vten_zadv_EULc(k) = -fac * ( wfldint(k) * ( v3m1(k) - v3m1(k-1) ) )
       uten_zadv_EULc(k) = -fac * ( wfldint(k) * ( u3m1(k) - u3m1(k-1) ) )
       do m = 1, pcnst
          qten_zadv_EULc(k,m) = -fac * ( wfldint(k) * ( q3m1(k,m) - q3m1(k-1,m) ) )
       end do

     ! ------------------------------------- !
     ! Manupulate individual forcings before ! 
     ! computing the final forecasted state  !
     ! ------------------------------------- !

     ! -------------------------------------------------------------- !
     ! Select the type of vertical advection :  EULc,IOP,OFF supported!
     ! -------------------------------------------------------------- !

     select case (scm_zadv_T)
     case ('iop')
           if (have_vertdivt) then
              tten_zadv(:plev) = vertdivt(:plev)
           else
              call endrun( subname//':: user set scm_zadv_tten to iop but vertdivt not on file')
           end if
     case ('eulc')
           tten_zadv(:) = tten_zadv_EULc(:) + tten_comp_EUL(:)
     case ('off')
           tten_zadv(:) = 0._r8
     end select

     select case (scm_zadv_uv)
     case ('iop')
           if (have_vertdivu .and. have_vertdivv) then
              uten_zadv(:) = vertdivu(:)
              vten_zadv(:) = vertdivv(:)
           else
              call endrun( subname//':: user set scm_zadv_uv to iop but vertdivu/v not on file')
           end if
     case ('eulc')
           uten_zadv(:) = uten_zadv_EULc(:)
           vten_zadv(:) = vten_zadv_EULc(:)
     case ('off')
           uten_zadv(:) = 0._r8
           vten_zadv(:) = 0._r8
     end select

     select case (scm_zadv_q)
     case ('iop')
           if (have_vertdivq) then
              qten_zadv(:plev,:pcnst) = vertdivq(:plev,:pcnst)
           else
              call endrun( subname//':: user set scm_zadv_qten to iop but vertdivq not on file')
           end if
     case ('eulc')
           qten_zadv(:plev,:pcnst) = qten_zadv_EULc(:plev,:pcnst)
     case ('slt')
        qten_zadv = qten_zadv_SLT
     case ('off')
        qten_zadv = 0._r8
     end select

     ! -------------------------------------------------------------- !
     ! Check horizontal advection u,v,t,q                             !
     ! -------------------------------------------------------------- !
     if (.not. have_divu) divu=0._r8 
     if (.not. have_divv) divv=0._r8 
     if (.not. have_divt) divt=0._r8 
     if (.not. have_divq) divq=0._r8 

     ! ----------------------------------- !
     !                                     ! 
     ! Compute the final forecasted states !
     !                                     !
     ! ----------------------------------- ! 
     ! make sure we have everything        !
     ! ----------------------------------- ! 

       if( .not. scm_use_obs_uv .and. .not. have_divu .and. .not. have_divv ) then 
              call endrun( subname//':: divu and divv not on the iop Unable to forecast Wind Set &
                                     scm_use_obs_uv=true to use observed u and v')
       end if
       if( .not. scm_use_obs_T .and. .not. have_divt) then
              call endrun( subname//':: divt not on the dataset. Unable to forecast Temperature. Stopping')
       end if
       if( .not. scm_use_obs_qv .and. .not. have_divq) then
              call endrun( subname//':: divq not on the dataset. Unable to forecast Humidity. Stopping')
       end if

           
       nstep_curr  = get_nstep()
  
          do k = 1, plev
             tfcst(k) = t3m2(k) + ztodt * ( tten_phys(k) + divt(k) + tten_zadv(k) )
             ufcst(k) = u3m2(k) + ztodt * ( uten_phys(k) + divu(k) + uten_zadv(k) )
             vfcst(k) = v3m2(k) + ztodt * ( vten_phys(k) + divv(k) + vten_zadv(k) )
             do m = 1, pcnst
                qfcst(1,k,m) = q3m2(k,m) + ztodt * ( qten_phys(k,m) + divq(k,m) + qten_zadv(k,m) ) 
             enddo
          enddo

       else
       !-------------------------------------
       ! This is the use_ana_iop=.TRUE. block
       !-------------------------------------

           nstep_curr  = get_nstep()
  
            if (is_first_step()) then
              u3m2 = u_ana
              v3m2 = v_ana
              t3m2 = t_ana
              q3m2(:,1) = q_ana
              psm2 = ps_ana
           endif
              

           ! -----------------------------------------------------
           !  Applied tendencies are in two 
           !  categories: 1) physics (includes nudging);  
           !  and 2) dynamics.   Dynamics tendencies are
           !  grouped and then scaled by dynfrcp. This is
           !  to allow removal of unreliable analysis driven 
           !  dynamics tendencies above some pressure, 
           !  typically <~ 10Pa.
           !------------------------------------------------------
           dynfrcp(:) = 1._r8
           where( pmidm1 < 10._r8)  ! changed from 10. Test.
              dynfrcp = 0._r8
           end where
           !------------------------------------------------------
           fcoriol     =  2._r8 * OOmega * sin( scmlat * PI/180._r8 )
           uten_coriol =  fcoriol * v3m2  
           vten_coriol = -fcoriol * u3m2 
           nsubdyn = 1
           vfcst = v3m2
           ufcst = u3m2
           ztodtn = ztodt/nsubdyn
           do nt= 1, nsubdyn
           do k = 1, plev
              ufcst(k) = ufcst(k) + ztodtn * ( uten_phys(k)                         & 
                                            + dynfrcp(k) *                          &
                                            ( uten_hadv_ana(k) + uten_vadv_ana(k)   & 
                                            + uten_vort_ana(k)                        & 
                                        !!  + fcoriol * vfcstm2(k)                    & 
                                            + uten_pfrc_ana(k) )  )
              vfcst(k) = vfcst(k) + ztodtn * ( vten_phys(k)                          & 
                                            + dynfrcp(k) *                          &
                                            ( vten_hadv_ana(k) + vten_vadv_ana(k)   & 
                                            + vten_vort_ana(k)                        & 
                                        !!  - fcoriol * ufcstm2(k)            & 
                                            + vten_pfrc_ana(k) )  )
           end do
           ufcstm2 = ufcst
           vfcstm2 = vfcst
           end do


           ufcor_0 = ufcst
           vfcor_0 = vfcst

#if 0
           !  Implicit formulation of Coriolis terms
           nsubdyn = 1
           ztodtn = ztodt/nsubdyn
           AA = 1._r8/(1._r8 + (ztodtn*fcoriol)**2 )
           do nt= 1, nsubdyn
           do k = 1, plev
              ufcst(k) = dynfrcp(k) * AA * ( ufcstm2(k) + ztodtn*fcoriol*vfcstm2(k) ) &
                            + (1._r8 - dynfrcp(k) )*ufcst(k)
              vfcst(k) = dynfrcp(k) * AA * ( vfcstm2(k) - ztodtn*fcoriol*ufcstm2(k) ) &
                            + (1._r8 - dynfrcp(k) )*vfcst(k)
           end do
           ufcstm2 = ufcst
           vfcstm2 = vfcst
           end do

           uten_vort_ana   = (ufcst - ufcor_0 )/ztodt
           vten_vort_ana   = (vfcst - vfcor_0 )/ztodt
#endif 

           uten_totdyn_ana = uten_hadv_ana + uten_vort_ana + uten_pfrc_ana + uten_vadv_ana
           vten_totdyn_ana = vten_hadv_ana + vten_vort_ana + vten_pfrc_ana + vten_vadv_ana

#if 1
           !----------------------------
           ! Calculate "usual" T-tendencies from complete IOP-file anyway
           !----------------------------
              ! ---------------------------------------------------------------------------- !
              ! Compute 'omega'( wfldint ) at the interface from the value at the mid-point. ! 
              ! SCAM-IOP file must provide omega at the mid-point not at the interface.      !
              ! ---------------------------------------------------------------------------- !
              wfldint(1) = 0._r8
              do k = 2, plev
                 weight = ( pintm1(k) - pmidm1(k-1) ) / ( pmidm1(k) - pmidm1(k-1) )
                 wfldint(k) = ( 1._r8 - weight ) * wfld(k-1) + weight * wfld(k)
              enddo
              wfldint(plevp) = 0._r8
              ! ------------------------------------------------------------ ! 
              ! Compute Eulerian compression heating due to vertical motion. !
              ! ------------------------------------------------------------ !
              do k = 1, plev
                 tten_comp_EUL(k) = wfld(k) * t3m1(k) * rair / ( cpair * pmidm1(k) )
              enddo
              ! ---------------------------------------------------------------------------- !
              ! Compute Centered Eulerian vertical advective tendencies for all 't, u, v, q' ! 
              ! ---------------------------------------------------------------------------- ! 
              do k = 2, plev - 1
                fac = 1._r8 / ( 2.0_r8 * pdelm1(k) )
                tten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( t3m1(k+1) - t3m1(k) ) + wfldint(k) * ( t3m1(k) - t3m1(k-1) ) )
              end do
              k = 1
              fac = 1._r8 / ( 2.0_r8 * pdelm1(k) )
              tten_zadv_EULc(k) = -fac * ( wfldint(k+1) * ( t3m1(k+1) - t3m1(k) ) )
              k = plev
              fac = 1._r8 / ( 2.0_r8 * pdelm1(k) )
              tten_zadv_EULc(k) = -fac * ( wfldint(k) * ( t3m1(k) - t3m1(k-1) ) )
              !----------------------------------------
              ! Replace ERA-derived T-tendencies with 
              ! IOP-file derived T-tendencies
              !----------------------------------------
              !!tten_vadv_ana(:) = tten_zadv_EULc(:)
              !!tten_comp_ana(:) = tten_comp_EUL(:)
              !!tten_hadv_ana(:) = divt(:)
              !-------------------
              ! For output
              !--------------------
              tten_zadv(:)     =  tten_zadv_EULc(:)
           !----------------------------
           ! End of Calculate "usual" T-tendencies from complete IOP-file anyway
           !----------------------------
#endif



           if (l_use_reconst_ttend) then
           do k=1,plev 
              tfcst(k) = t3m2(k) + ztodt * ( tten_phys(k)                          &
                                          + dynfrcp(k) *                          &
                                          ( tten_hadv_ana(k) + tten_vadv_ana(k)   & 
                                          + tten_comp_ana(k)  ) )
           end do
           end if

           if (l_use_direct_ttend) then
           do k=1,plev 
              tfcst(k) = t3m2(k) + ztodt * ( tten_phys(k)                          &
                                          + dynfrcp(k) *                          &
                                          ( tten_dycore_ana(k) )    )
           end do
           end if

           do k=1,plev 
              do m = 1, 1
                 qfcst(1,k,m) = q3m2(k,m) + ztodt * ( qten_phys(k,m)                       & 
                                            + dynfrcp(k) *                          &
                                            ( qten_hadv_ana(k) + qten_vadv_ana(k) ) ) 
              enddo
            enddo

            ps = ps_ana
               
            write(*,*) " Nstep "  ,nstep_curr
            if (mod( nstep_curr,10)==0) then
               !ufcst = 0.5*(ufcst+u3m1)
               !vfcst = 0.5*(vfcst+v3m1)
            endif

            ! Zero-out NON ana_iop diagnostics
            !    ???? 

      end if  ! END use_ana_iop IF block
     
     ! This code is executed regardless of use_ana_iop value
     ! ------------------ !
     ! Diagnostic Outputs !
     ! ------------------ !
       call outfld( 'TOBS' ,       tobs,                  plon, dummy_dyndecomp )
       call outfld( 'UOBS' ,       uobs,                  plon, dummy_dyndecomp )
       call outfld( 'VOBS' ,       vobs,                  plon, dummy_dyndecomp )
       call outfld( 'TTEN_XYADV' , divt,                  plon, dummy_dyndecomp )
       call outfld( 'UTEN_XYADV' , divu,                  plon, dummy_dyndecomp )
       call outfld( 'VTEN_XYADV' , divv,                  plon, dummy_dyndecomp )
       call outfld( 'QVTEN_XYADV', divq(:,1),             plon, dummy_dyndecomp )
       if (.not.adiabatic) then
          call outfld( 'QLTEN_XYADV', divq(:,ixcldliq),      plon, dummy_dyndecomp )
          call outfld( 'QITEN_XYADV', divq(:,ixcldice),      plon, dummy_dyndecomp )
          call outfld( 'NLTEN_XYADV', divq(:,ixnumliq),      plon, dummy_dyndecomp )
          call outfld( 'NITEN_XYADV', divq(:,ixnumice),      plon, dummy_dyndecomp )
          call outfld( 'QLTEN_ZADV' , qten_zadv(:,ixcldliq), plon, dummy_dyndecomp )
          call outfld( 'QITEN_ZADV' , qten_zadv(:,ixcldice), plon, dummy_dyndecomp )
          call outfld( 'NLTEN_ZADV' , qten_zadv(:,ixnumliq), plon, dummy_dyndecomp )
          call outfld( 'NITEN_ZADV' , qten_zadv(:,ixnumice), plon, dummy_dyndecomp )
          call outfld( 'QLTEN_PHYS' , qten_phys(:,ixcldliq), plon, dummy )
          call outfld( 'QITEN_PHYS' , qten_phys(:,ixcldice), plon, dummy )
          call outfld( 'NLTEN_PHYS' , qten_phys(:,ixnumliq), plon, dummy )
          call outfld( 'NITEN_PHYS' , qten_phys(:,ixnumice), plon, dummy )
       end if
       call outfld( 'TTEN_ZADV'  , tten_zadv,             plon, dummy_dyndecomp )
       call outfld( 'UTEN_ZADV'  , uten_zadv,             plon, dummy_dyndecomp )
       call outfld( 'VTEN_ZADV'  , vten_zadv,             plon, dummy_dyndecomp )
       call outfld( 'QVTEN_ZADV' , qten_zadv(:,1),        plon, dummy_dyndecomp )
       !call outfld( 'TTEN_ZADV'  , vertdivt,             plon, dummy_dyndecomp )
       !call outfld( 'QVTEN_ZADV' , vertdivq(:,1),        plon, dummy_dyndecomp )

       call outfld( 'TTEN_COMP_IOP', tten_comp_eul,      plon, dummy_dyndecomp )

       call outfld( 'TTEN_PHYS'  , tten_phys,             plon, dummy_dyndecomp )
       call outfld( 'UTEN_PHYS'  , uten_phys,             plon, dummy_dyndecomp )
       call outfld( 'VTEN_PHYS'  , vten_phys,             plon, dummy_dyndecomp )
       call outfld( 'QVTEN_PHYS' , qten_phys(:,1),        plon, dummy_dyndecomp )

    end if  ! END of use_camiop IF BLOCK

!!!!#if 0
!+++ARH
  !if( .NOT.(use_ana_iop) ) then
  if( .NOT.(scm_use_ana_iop) ) then
!---ARH
  ! ---------------------------------------------------------------- !
  ! Used the SCAM-IOP-specified state instead of forecasted state    !
  ! at each time step if specified by the switch.                    !
  ! If SCAM-IOP has 't,u,v,q' profile at a single initial time step. !
  ! ---------------------------------------------------------------- ! 
  if( scm_use_obs_T .and. have_t ) then 
     do k = 1, plev
        tfcst(k) = tobs(k)
     enddo
  endif
  
  if( scm_use_obs_uv .and. have_u .and. have_v  ) then 
     do k = 1, plev
        ufcst(k) = uobs(k)
        vfcst(k) = vobs(k)
     enddo
  endif
  
  if( scm_use_obs_qv .and. have_q ) then 
     do k = 1, plev
        qfcst(1,k,1) = qobs(k)
     enddo
  endif
  
  ! ------------------------------------------------------------------- !
  ! Relaxation to the observed or specified state                       !
  ! We should specify relaxation time scale ( rtau ) and                !
  ! target-relaxation state ( in the current case, either 'obs' or 0 )  !
  ! ------------------------------------------------------------------- !
  
  relax_T(:)             = 0._r8
  relax_u(:)             = 0._r8
  relax_v(:)             = 0._r8
  relax_q(:plev,:pcnst)  = 0._r8
  ! +++BPM: allow linear relaxation profile
  ! scm_relaxation is a logical from scamMod
  ! scm_relax_tau_top_sec and scm_relax_tau_bot_sec are the relaxation times at top and bottom of layer
  ! also defined in scamMod
  if ( scm_relaxation.and.scm_relax_linear ) then
     rslope = (scm_relax_top_p - scm_relax_bot_p)/(scm_relax_tau_top_sec - scm_relax_tau_bot_sec)
     rycept = scm_relax_tau_top_sec - (rslope*scm_relax_top_p)
  endif

  !    prepare scm_relax_fincl for comparison in scmforecast.F90
  scm_fincl_empty=.true.
  do i=1,pcnst
     if (len_trim(scm_relax_fincl(i)) > 0) then
        scm_fincl_empty=.false.
        scm_relax_fincl(i)=trim(to_upper(scm_relax_fincl(i)))
     end if
  end do

  do k = 1, plev
     if( scm_relaxation ) then
        if ( pmidm1(k).le.scm_relax_bot_p.and.pmidm1(k).ge.scm_relax_top_p ) then ! inside layer
           if (scm_relax_linear) then
              rtau(k) = rslope*pmidm1(k) + rycept ! linear regime
           else
              rtau(k)         = max( ztodt, scm_relax_tau_sec ) ! constant for whole layer / no relax outside
           endif
        else if  (scm_relax_linear .and. pmidm1(k).le.scm_relax_top_p ) then ! not linear => do nothing / linear => use upper value
           rtau(k) = scm_relax_tau_top_sec ! above layer keep rtau equal to the top
        endif
        ! +BPM: this can't be the best way...
        ! I put this in because if rtau doesn't get set above, then I don't want to do any relaxation in that layer.
        ! maybe the logic of this whole loop needs to be re-thinked. 
        if (rtau(k).ne.0) then
           relax_T(k)      = -  ( tfcst(k)     - tobs(k) )    / rtau(k)
           relax_u(k)      = -  ( ufcst(k)     - uobs(k) )    / rtau(k)
           relax_v(k)      = -  ( vfcst(k)     - vobs(k) )    / rtau(k)         
           relax_q(k,1)    = -  ( qfcst(1,k,1) - qobs(k) )    / rtau(k)
           do m = 2, pcnst
              relax_q(k,m) = -  ( qfcst(1,k,m) - qinitobs(k,m)   )    / rtau(k)
           enddo
           if (scm_fincl_empty .or. ANY(scm_relax_fincl(:).eq.'T')) &
                tfcst(k)        =      tfcst(k)     + relax_T(k)   * ztodt
           if (scm_fincl_empty .or.ANY(scm_relax_fincl(:).eq.'U')) &
                ufcst(k)        =      ufcst(k)     + relax_u(k)   * ztodt
           if (scm_fincl_empty .or. ANY(scm_relax_fincl(:).eq.'V')) &
                vfcst(k)        =      vfcst(k)     + relax_v(k)   * ztodt
           do m = 1, pcnst
              if (scm_fincl_empty .or. ANY(scm_relax_fincl(:) .eq. trim(to_upper(cnst_name(m)))) ) then
                 qfcst(1,k,m) =      qfcst(1,k,m) + relax_q(k,m) * ztodt
              end if
           enddo
        end if
     endif
  enddo
  call outfld( 'TRELAX'   , relax_T           , plon, dummy )
  call outfld( 'QRELAX'   , relax_q(1:plev,1) , plon, dummy )
  call outfld( 'TAURELAX' , rtau              , plon, dummy )
!!!#endif 
  end if ! END of 2nd use_ana_iop BLOCK (exec for use_ana_iop=.F.) 

  ! --------------------------------------------------------- !
  ! Assign the final forecasted state to the output variables !
  ! --------------------------------------------------------- !
  
  t3(1:plev)         = tfcst(1:plev)
  u3(1:plev)         = ufcst(1:plev)
  v3(1:plev)         = vfcst(1:plev)

!+++ARH  
  !if (use_ana_iop) then 
  if (scm_use_ana_iop) then
!---ARH
     q3(1:plev,1:1) = qfcst(1,1:plev,1:1)
  else
     q3(1:plev,1:pcnst) = qfcst(1,1:plev,1:pcnst)
  endif

  tdiff(1:plev)  =   t3(1:plev)   - tobs(1:plev)
  qdiff(1:plev)  =   q3(1:plev,1) - qobs(1:plev)
 

  call outfld( 'QDIFF'  , qdiff,             plon, dummy_dyndecomp )
  call outfld( 'TDIFF'  , tdiff,             plon, dummy_dyndecomp )

#ifdef SCAMNUDGERUN
  call outfld( 'OMEGA_IOP' , wfld,           plon, dummy_dyndecomp )
  call outfld( 'OMEGA_ANA' , omega_ana,      plon, dummy_dyndecomp )
  call outfld( 'ETAD_ANA'  ,  etad_ana,      plon, dummy_dyndecomp )
  call outfld( 'ZETA_ANA'  ,  zeta_ana,      plon, dummy_dyndecomp )
  call outfld( 'T_ANA'     ,  T_ana_diag,    plon, dummy_dyndecomp )
  call outfld( 'Q_ANA'     ,  Q_ana_diag,    plon, dummy_dyndecomp )
  call outfld( 'TV_ANA'    ,  Tv_ana,        plon, dummy_dyndecomp )
  call outfld( 'U_ANA'     ,  U_ana_diag,    plon, dummy_dyndecomp )
  call outfld( 'V_ANA'     ,  V_ana_diag,    plon, dummy_dyndecomp )

  call outfld( 'UTEN_CORIOL' ,  uten_coriol,  plon, dummy_dyndecomp )
  call outfld( 'VTEN_CORIOL' ,  vten_coriol,  plon, dummy_dyndecomp )

  call outfld( 'UTEN_TOTDYN_ANA' ,  uten_totdyn_ana,  plon, dummy_dyndecomp )
  call outfld( 'VTEN_TOTDYN_ANA' ,  vten_totdyn_ana,  plon, dummy_dyndecomp )
  call outfld( 'TTEN_TOTDYN_ANA' ,  tten_totdyn_ana,  plon, dummy_dyndecomp )
  call outfld( 'QTEN_TOTDYN_ANA' ,  qten_totdyn_ana,  plon, dummy_dyndecomp )

  call outfld( 'UTEN_TOTDYN_ANAR' ,  uten_totdyn_anax,  plon, dummy_dyndecomp )
  call outfld( 'VTEN_TOTDYN_ANAR' ,  vten_totdyn_anax,  plon, dummy_dyndecomp )

  call outfld( 'UTEN_DYCORE_ANA' ,  uten_dycore_ana,  plon, dummy_dyndecomp )
  call outfld( 'VTEN_DYCORE_ANA' ,  vten_dycore_ana,  plon, dummy_dyndecomp )
  call outfld( 'TTEN_DYCORE_ANA' ,  tten_dycore_ana,  plon, dummy_dyndecomp )
  call outfld( 'OMEGA_DYCORE_ANA',  omega_dycore_ana, plon, dummy_dyndecomp )
  call outfld( 'OMEGA_RECALC_ANA',  omega_recalc_ana, plon, dummy_dyndecomp )

  call outfld( 'UTEN_HADV_ANA' ,  uten_hadv_ana,  plon, dummy_dyndecomp )
  call outfld( 'UTEN_VADV_ANA' ,  uten_vadv_ana,  plon, dummy_dyndecomp )
  call outfld( 'UTEN_VORT_ANA' ,  uten_vort_ana,  plon, dummy_dyndecomp )
  call outfld( 'UTEN_KEG_ANA' ,   uten_keg_ana,  plon, dummy_dyndecomp )
  call outfld( 'UTEN_PFRC_ANA' ,  uten_pfrc_ana,  plon, dummy_dyndecomp )
  call outfld( 'UTEN_PRG_ANA'  ,  uten_prg_ana,   plon, dummy_dyndecomp )
  call outfld( 'UTEN_PHIG_ANA' ,  uten_phig_ana,  plon, dummy_dyndecomp )

  call outfld( 'VTEN_HADV_ANA' ,  vten_hadv_ana,  plon, dummy_dyndecomp )
  call outfld( 'VTEN_VADV_ANA' ,  vten_vadv_ana,  plon, dummy_dyndecomp )
  call outfld( 'VTEN_VORT_ANA' ,  vten_vort_ana,  plon, dummy_dyndecomp )
  call outfld( 'VTEN_PFRC_ANA' ,  vten_pfrc_ana,  plon, dummy_dyndecomp )

  call outfld( 'TTEN_HADV_ANA' ,  tten_hadv_ana,  plon, dummy_dyndecomp )
  call outfld( 'TTEN_VADV_ANA' ,  tten_vadv_ana,  plon, dummy_dyndecomp )
  call outfld( 'TTEN_COMP_ANA' ,  tten_comp_ana,  plon, dummy_dyndecomp )

  call outfld( 'QTEN_HADV_ANA' ,  qten_hadv_ana,  plon, dummy_dyndecomp )
  call outfld( 'QTEN_VADV_ANA' ,  qten_vadv_ana,  plon, dummy_dyndecomp )

  if (have_u) call outfld( 'U_IOP' ,  uobs,  plon, dummy_dyndecomp )
  if (have_u) call outfld( 'V_IOP' ,  vobs,  plon, dummy_dyndecomp )

#endif
  
   return

   end subroutine forecast


end module scmforecast
