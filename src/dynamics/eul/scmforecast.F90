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
                                 vertdivt,vertdivu,vertdivv,vobs,wfld,qinitobs,scm_relax_fincl
     use time_manager,     only : get_curr_calday, get_nstep, get_step_size, is_first_step
     use cam_abortutils,   only : endrun
     use string_utils,     only: to_upper

     implicit none

   ! ---------------------- !
   ! Parameters             !
   ! ---------------------- !

    character(len=*), parameter ::  subname = "forecast"

   ! --------------------------------------------------- !
   ! x = t, u, v, q                                      !
   ! x3m1 : state variable used for computing 'forcing'  !
   ! x3m2 : initial state variable before time-marching  !
   ! x3   : final   state variable after  time-marching  !
   ! --------------------------------------------------- !

     integer,  intent(in)    :: lat                     
     integer,  intent(in)    :: nlon                     
     real(r8), intent(in)    :: ztodt                   ! Twice time step unless nstep = 0 [ s ]

     real(r8), intent(in)    :: ps(plon)                ! Surface pressure [ Pa ]
     real(r8), intent(in)    :: psm1(plon)              ! Surface pressure [ Pa ]
     real(r8), intent(in)    :: psm2(plon)              ! Surface pressure [ Pa ]

     real(r8), intent(in)    :: t3m1(plev)              ! Temperature [ K ]
     real(r8), intent(in)    :: t3m2(plev)              ! Temperature [ K ]
     real(r8), intent(in)    :: u3m1(plev)              ! Zonal wind [ m/s ]
     real(r8), intent(in)    :: u3m2(plev)              ! Zonal wind [ m/s ]
     real(r8), intent(in)    :: v3m1(plev)              ! Meridional wind [ m/s ]
     real(r8), intent(in)    :: v3m2(plev)              ! Meridional wind [ m/s ]
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

   else

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

       do k = 1, plev
          tfcst(k) = t3m2(k) + ztodt * ( tten_phys(k) + divt(k) + tten_zadv(k) )
          ufcst(k) = u3m2(k) + ztodt * ( uten_phys(k) + divu(k) + uten_zadv(k) )
          vfcst(k) = v3m2(k) + ztodt * ( vten_phys(k) + divv(k) + vten_zadv(k) )
          do m = 1, pcnst
             qfcst(1,k,m) = q3m2(k,m) + ztodt * ( qten_phys(k,m) + divq(k,m) + qten_zadv(k,m) ) 
          enddo
       enddo

     ! ------------------ !
     ! Diagnostic Outputs !
     ! ------------------ !

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
       call outfld( 'TTEN_ZADV'  , vertdivt,             plon, dummy_dyndecomp )
       call outfld( 'QVTEN_ZADV' , vertdivq(:,1),        plon, dummy_dyndecomp )

       call outfld( 'TTEN_PHYS'  , tten_phys,             plon, dummy )
       call outfld( 'UTEN_PHYS'  , uten_phys,             plon, dummy )
       call outfld( 'VTEN_PHYS'  , vten_phys,             plon, dummy )
       call outfld( 'QVTEN_PHYS' , qten_phys(:,1),        plon, dummy )

  endif

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
  
  if( scm_use_obs_uv .and. have_u .and. have_v ) then 
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
  
  ! --------------------------------------------------------- !
  ! Assign the final forecasted state to the output variables !
  ! --------------------------------------------------------- !
  
  t3(1:plev)         = tfcst(1:plev)
  u3(1:plev)         = ufcst(1:plev)
  v3(1:plev)         = vfcst(1:plev)
  q3(1:plev,1:pcnst) = qfcst(1,1:plev,1:pcnst)
  
  tdiff(1:plev)  =   t3(1:plev)   - tobs(1:plev)
  qdiff(1:plev)  =   q3(1:plev,1) - qobs(1:plev)

  call outfld( 'QDIFF'  , qdiff,             plon, dummy_dyndecomp )
  call outfld( 'TDIFF'  , tdiff,             plon, dummy_dyndecomp )
  
   return

   end subroutine forecast
 end module scmforecast
