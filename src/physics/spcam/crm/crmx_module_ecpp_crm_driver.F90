module  crmx_module_ecpp_crm_driver
#ifdef ECPP
!------------------------------------------------------------------------
! F90 module to prepare CRM output for ECPP module in the MMF model. 
!
! This code was written originally by William Gustafson, and is adopted into 
! the MMF model by Minghuai Wang (minghuai.wang@pnl.gov), November, 2009.
!
! Assumptiont built into this code:
!
! Open issues:
!  - The mask for determining a "moving" or limited spatial average
!    is not implemented.
!  - The dependencies in Makefile don't work. If a compile fails,
!    try "make clean; make" instead to clear out the module files.
!  - For uv_in/out, a simple time average is being done and one can
!    argue that it should be a weighted average since the number of in
!    and out points changes with each time step. The affect is probably
!    small for short time averages though.
!  - When calculating the standard deviation of vertical velocity,
!    each cell is treated equally and the std. dev. is over the 3 dims
!    below the cloud tops. We may want to consider weighting each cell
!    by either its volume or mass.
!  - To get cloud values at vertical cell interface, a simple average
!    is being done when an interpolation should technically be done.
!    This only affects quiescent cloudy/clear categories.
!  - Ditto for getting the density at the vertical cell interface (rho8w).
! 
! Differences between the methodology here and in Ferret:
!  - When calculating wup_bar and wdown_bar, points with w==0 are ignored
!    here and were included in wup in Ferret.
!  - Clear fluxes are no longer chopped off at the cloud top.
!  - When calculating the std. dev. in and below the cloud, the level
!    just above the cloud top is now included so we include w out the
!    cloud top.
!  - When determining "cloudyother" in Ferret the cloud above the
!    interface was used. Now, the average of the cloud above and below
!    is used.
!
! William.Gustafson@pnl.gov; 20-Jul-2006
!  v2.0 - Added two-level time averaging, one for the stats and a longer
!         period for output.
!  v2.1 - 25-Jul-2006; Fixed sign bug with uv_in/out.
!
!  v3.0 - aug-sep-2006 - many changes by r.easter and s.ghan
!         major change is option for multiple up and downdraft classes
!
!  v3.1 - 02-nov-2006 r.easter - replaced uv_in/outsum with u_in/outsum
!                                & v_in/outsum
!
!  v4.0 - 25-Jan-2007, wig; 
!         - Added areaavgtype switch to output final areas either as
!           instantaneous, averaged over the last ntavg1 period of each
!           ntavg2 avg, or as averaged over ntavg2.
!         - Output areas as average over ntavg2 and also just at end
!           of it.
!         - Added entrainment averages to output (do not divide by dz).
!
!  postproc_wrfout_bb.f90 from postproc_wrfout.f90 - 15-nov-2007, rce;
!	- do multiple processings
!
!  v5.0 - Nov-2008, wig
!         - Major rewrite to include combinations of cloud, precipitation,
!           and transport classes
!         - Output format changes to multi-dimensional variables based
!           on the classes instead of outputting each class separately
!
! 14-Apr-2009, wig: Fixed bug with mode_updnthresh at model top for
!           bad calculation of w thresholds.
!
! 16-Apr-2009, wig: Added qcloud weighting to qlsink averages
!
!----------------------------------------------------------------------------------------
  use crmx_ecppvars
  use crmx_ecppvars,  only: QUI, UP1, DN1, NCLASS_TR, NCLASS_CL, CLR, CLD, NCLASS_PR, PRN, PRY
  use cam_abortutils, only: endrun

  public ecpp_crm_stat
  public ecpp_crm_init
  public ecpp_crm_cleanup

  integer, public :: ntavg1_ss, ntavg2_ss

  private 
  save

  integer ::  nxstag, nystag, nzstag
  integer :: itavg1, itavg2,  &
      ntavg1, ntavg2

  integer :: mode_updnthresh
  integer :: areaavgtype
                         ! Methodology to compute final area averages:
                         ! 0 = area categories based on instantaneous
                         !     values at last time step of ntavg2
                         ! 1 = area cat. based on last ntavg1 avgeraging
                         !     period of each ntavg2 period
                         ! 2 = area cat. based on average of full ntavg2
                         !     period
  integer :: plumetype
                         ! 1 = single plume
                         ! 2 = two plumes, core and weak
                         ! 3 = multi-plume, number based on setting of
                         !     allcomb
  logical :: allcomb
                         ! true if updrafts and downdrafts have all
                         ! combinations of bases and tops.
  real :: cloudthresh, &
          prcpthresh, &
          downthresh, downthresh2, &
          upthresh, upthresh2

  real :: cloudthresh_trans,   &   !  the threshold total cloud water for updraft or downdraft
          precthresh_trans         !  the threshold total rain, snow and graupel for clear, updraft or downdraft 

  integer, dimension(:),allocatable :: &
          updraftbase, updrafttop, dndrafttop, dndraftbase
  integer :: nupdraft, ndndraft
  integer :: ndraft_max, nupdraft_max, ndndraft_max

contains

!========================================================================================
subroutine ecpp_crm_init()

  use crmx_grid, only: nx, ny, nzm, dt
  use crmx_module_ecpp_stats, only: zero_out_sums1, zero_out_sums2
  use module_ecpp_ppdriver2, only: nupdraft_in, ndndraft_in, ncls_ecpp_in
  implicit none

  integer :: kbase, ktop
  integer :: m
  integer :: nup, ndn
  character(len=100) :: msg

  nxstag = nx+1
  nystag = ny+1
  nzstag = nzm+1

! ntavg1_ss and ntavg1_ss are defined in crm.F90 in the MMF model.
! ntavg1_ss = dt_gl   ! GCM time step
! ntavg1_ss = number of seconds to average between computing categories.
!  ntavg2_ss = dt_gl   ! GCM time step
! ntavg2_ss = number of seconds to average between outputs.
!    This must be a multiple of ntavgt1_ss.

  mode_updnthresh = 16
!     1 = method originally implemented by Bill G
!         wup_thresh   =  wup_stddev*abs(upthresh)
!         wdown_thresh = -wdown_stddev*abs(downthresh)
!     2 = similar to 1, but include the mean wup and wdown
!         wup_thresh   = wup_bar   + wup_stddev*abs(upthresh)
!         wdown_thresh = wdown_bar - wdown_stddev*abs(downthresh)
!     3 = user specifies an absolute threshold
!         wup_thresh   =  abs(upthresh)
!         wdown_thresh = -abs(downthresh)
!     4 = similar to 1, but do
!         wup_thresh   =  wup_rms*abs(upthresh)
!         wdown_thresh = -wdown_rms*abs(downthresh)
!
!     5     = see description  in module_ecpp_stats.f90
!     6,  7 = see descriptions in module_ecpp_stats.f90
!     8,  9 = see descriptions in module_ecpp_stats.f90
!    10, 11 = see descriptions in module_ecpp_stats.f90
!    12, 13 = see descriptions in module_ecpp_stats.f90

  upthresh    = 1.    !Multiples of std. dev. to classify as updraft
  downthresh  = 1.    !Multiples of std. dev. to classify as downdraft
  upthresh2   = 0.5   ! ...ditto, except for weaker 2nd draft type when plumetype=2
  downthresh2 = 0.5

#ifdef CLUBB_CRM
  cloudthresh = 2e-7  !Cloud mixing ratio beyond which cell is "cloudy(liquid)" (kg/kg)
                      ! As now fractional cloudiness is used for classifying cloudy vs. clear, 
                      ! reduce it from 1.0e-6 to 2.0e-7
#else 
  cloudthresh = 1e-6  !Cloud mixing ratio beyond which cell is "cloudy(liquid)" (kg/kg)
#endif

  prcpthresh  = 1e-6  !Preciptation rate (precr) beyond which cell is raining (kg/m2/s)
                      ! this is used to classify precipitating vs. nonprecipitating class for wet scavenging. 

!+++mhwang
! high thresholds are used to classify transport classes (following Xu et al., 2002, Q.J.R.M.S. 
!
  cloudthresh_trans = 1e-5  !Cloud mixing ratio beyond which cell is "cloudy" to classify transport classes (kg/kg)   +++mhwang
                            ! the maxium of cloudthres_trans and 0.01*qvs is used to classify transport class  
  precthresh_trans  = 1e-4  !Preciptation mixing ratio beyond which cell is raining to classify transport classes (kg/kg)  !+++mwhang
!---mhwang

  areaavgtype= 1    !final area avg over 0=instantaneous, 1=ntavg1, 2=ntavg2
  plumetype  = 1    !1 for single plume, 2 for core and weak plumes, 3 for multiple plumes
  allcomb = .false. !true for all combinations of plume bases and tops, false for 1 plume per base

!----------------------------------------------------------------------------------
! Sanity check...
!----------------------------------------------------------------------------------

  if(plumetype>3)then
     msg = 'ecpp_crm, plumetype must be <=3'
     call endrun(trim(msg)) 
  endif

  if(plumetype<3 .and. allcomb)then
     msg='ecpp_crm, allcomb=true requires plumetype=3'
     call endrun(trim(msg))
  endif

  if(areaavgtype>2)then
     msg='ecpp_crm, areaavgtype must be <=2'
     call endrun(trim(msg))
  endif

  if ((mode_updnthresh < 1) .or. (mode_updnthresh > 17)) then
     msg='ecpp_crm, error - must have   1 <= mode_updnthresh <= 17'
     call endrun(trim(msg))
  endif

  if( abs(upthresh2) > 0.90*abs(upthresh) ) then
     msg='ecpp_crm, error - upthresh2 must be < 0.90*upthresh'
     call endrun(trim(msg))
  end if

  if( abs(downthresh2) > 0.90*abs(downthresh) ) then
     msg='ecpp_crm, error - downthresh2 must be < 0.90*downthresh'
     call endrun(trim(msg))
  end if

! determine number of updrafts and downdrafts
!
! updraft kbase & ktop definition:
!    ww(i,j,k ) >  wup_thresh for k=kbase+1:ktop
!               <= wup_thresh at  k=kbase and k=ktop+1
! they identify the "T-points" which enclose the updraft "W-points"
!    and are affected by the subgrid transport of this updraft
!
! downdraft kbase & ktop definition:
!    ww(i,j,k ) <  wdown_thresh for k=kbase+1:ktop
!               >= wdown_thresh at  k=kbase and k=ktop+1
! they identify the "T-points" which enclose the downdraft "W-points"
!    and are affected by the subgrid transport of this downdraft
!
! for both updrafts and downdrafts,
!    1 <= kbase < ktop < nzstag

  nupdraft = 0
  ndndraft = 0
  nupdraft_max = 0
  ndndraft_max = 0

  select case (plumetype)
  case (1) !single plume
     nupdraft = 1
     ndndraft = 1
  case (2) !core and weak plumes
     nupdraft = 2
     ndndraft = 2
  case (3)
     do kbase=1,nzm-1
        if(allcomb)then   ! all possible tops
           nupdraft=nupdraft+nzm-kbase
        else              ! one top per base
           nupdraft=nupdraft+1
        endif
     enddo
     do ktop=nzm,2,-1
        if(allcomb)then   ! all possible bases
           ndndraft=ndndraft+ktop-1
        else              ! one base per top
           ndndraft=ndndraft+1
        endif
     enddo
  end select

  nupdraft_max = max( nupdraft_max, nupdraft )
  ndndraft_max = max( ndndraft_max, ndndraft )

  DN1 = nupdraft + 2 !Setup index of first downdraft class
  NCLASS_TR = nupdraft + ndndraft + 1

  ndraft_max = 1 + nupdraft_max + ndndraft_max

  if(NCLASS_TR.ne.ncls_ecpp_in) then
     call endrun('NCLASS_TR should be equal to ncls_ecpp_in')
  end if 
  if((nupdraft.ne.nupdraft_in) .or. (ndndraft.ne.ndndraft_in)) then
     call endrun('nupdraft or ndndraft is not set correctly')
  end if

  allocate (updraftbase(nupdraft_max), &
            updrafttop( nupdraft_max) )
  allocate (dndraftbase(ndndraft_max), &
            dndrafttop( ndndraft_max) )

  select case (plumetype)
  case (1) !single plume
     updraftbase(1)=1
     updrafttop( 1)=nzm
     dndrafttop( 1)=nzm
     dndraftbase(1)=1
  case (2)
     updraftbase(1:2)=1
     updrafttop( 1:2)=nzm
     dndrafttop( 1:2)=nzm
     dndraftbase(1:2)=1
  case (3)
     m=0
     do kbase=1,nzm-1
        if(allcomb)then     ! loop over all possible tops.
           do ktop=kbase+1,nzm
              m=m+1
              updraftbase(m)=kbase
              updrafttop( m)=ktop
           enddo
        else                ! only one top per base
           m=m+1
           updraftbase(m)=kbase
           updrafttop( m)=nzm
        endif
     enddo

     m=0
     do ktop=nzm,2,-1
        if(allcomb)then    ! loop over all possible bases.
           do kbase=ktop-1,1,-1
              m=m+1
              dndrafttop( m)=ktop
              dndraftbase(m)=kbase
           enddo
        else               ! only one base per top
           m=m+1
           dndrafttop( m)=ktop
           dndraftbase(m)=1
        endif
     enddo
  end select

!---------------------------------------------------------------------------
! Allocate arrays
!---------------------------------------------------------------------------
  allocate( qlsink(nx,ny,nzm), precr(nx,ny,nzm), precsolid(nx,ny,nzm), rh(nx, ny, nzm), qvs(nx, ny, nzm))

  allocate( qlsink_bf(nx, ny, nzm), prain(nx, ny, nzm), qcloud_bf(nx, ny, nzm))

  allocate( qcloudsum1(nx,ny,nzm), qcloud_bfsum1(nx,ny,nzm), qrainsum1(nx,ny,nzm), &
            qicesum1(nx,ny,nzm), qsnowsum1(nx,ny,nzm), qgraupsum1(nx,ny,nzm), &
            qlsinksum1(nx,ny,nzm), precrsum1(nx,ny,nzm), &
            precsolidsum1(nx,ny,nzm), precallsum1(nx,ny,nzm), &
            altsum1(nx,ny,nzm), rhsum1(nx,ny,nzm), cf3dsum1(nx,ny,nzm), &
            wwsum1(nx,ny,nzstag), wwsqsum1(nx,ny,nzstag), &
            tkesgssum1(nx, ny, nzm), qlsink_bfsum1(nx, ny, nzm), prainsum1(nx, ny, nzm), qvssum1(nx, ny, nzm) )

  allocate(           &
            xkhvsum(nzm) )

  allocate( wwqui_cen_sum(nzm), wwqui_bnd_sum(nzm+1),  &
            wwqui_cloudy_cen_sum(nzm), wwqui_cloudy_bnd_sum(nzm+1))

  allocate( wup_thresh(nzm+1), wdown_thresh(nzm+1))

  allocate( area_bnd_final( nzstag,NCLASS_CL,ndraft_max,NCLASS_PR), &
            area_bnd_sum(   nzstag,NCLASS_CL,ndraft_max,NCLASS_PR), &
            area_cen_final( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            area_cen_sum(   nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            mass_bnd_final( nzstag,NCLASS_CL,ndraft_max,NCLASS_PR), &
            mass_bnd_sum(   nzstag,NCLASS_CL,ndraft_max,NCLASS_PR), &
            mass_cen_final( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            mass_cen_sum(   nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            ent_bnd_sum(    nzstag,NCLASS_CL,ndraft_max,NCLASS_PR), &
            rh_cen_sum(     nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qcloud_cen_sum( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qcloud_bf_cen_sum( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qrain_cen_sum(  nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qice_cen_sum(   nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qsnow_cen_sum(  nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qgraup_cen_sum( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qlsink_cen_sum( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            precr_cen_sum(  nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            precsolid_cen_sum(nzm  ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            precall_cen_sum(nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qlsink_bf_cen_sum( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            qlsink_avg_cen_sum( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR), &
            prain_cen_sum( nzm    ,NCLASS_CL,ndraft_max,NCLASS_PR)  )

! Initialize the running sums.
  call zero_out_sums1( qcloudsum1(:,:,:), qcloud_bfsum1(:,:,:), qrainsum1(:,:,:), &
                       qicesum1(:,:,:), qsnowsum1(:,:,:), qgraupsum1(:,:,:), &
                       qlsinksum1(:,:,:), precrsum1(:,:,:), &
                       precsolidsum1(:,:,:), precallsum1(:,:,:), &
                       altsum1(:,:,:), rhsum1(:,:,:), cf3dsum1(:,:,:), &
                       wwsum1(:,:,:), wwsqsum1(:,:,:), tkesgssum1(:,:,:), &
                       qlsink_bfsum1(:,:,:), prainsum1(:,:,:), qvssum1(:,:,:) )
  ndn = ndndraft ; nup = nupdraft
  call zero_out_sums2( &
       xkhvsum(:), &
       wwqui_cen_sum(:), wwqui_bnd_sum(:), wwqui_cloudy_cen_sum(:), wwqui_cloudy_bnd_sum(:),  &
       area_bnd_final(:,:,1:1+nup+ndn,:), area_bnd_sum(:,:,1:1+nup+ndn,:), &
       area_cen_final(:,:,1:1+nup+ndn,:), area_cen_sum(:,:,1:1+nup+ndn,:), &
       mass_bnd_final(:,:,1:1+nup+ndn,:), mass_bnd_sum(:,:,1:1+nup+ndn,:), &
       mass_cen_final(:,:,1:1+nup+ndn,:), mass_cen_sum(:,:,1:1+nup+ndn,:), &
       ent_bnd_sum(:,:,1:1+nup+ndn,:), &
       rh_cen_sum(:,:,1:1+nup+ndn,:), &
       qcloud_cen_sum(:,:,1:1+nup+ndn,:), qcloud_bf_cen_sum(:,:,1:1+nup+ndn,:), qrain_cen_sum(:,:,1:1+nup+ndn,:), &
       qice_cen_sum(:,:,1:1+nup+ndn,:),  qsnow_cen_sum(:,:,1:1+nup+ndn,:), &
       qgraup_cen_sum(:,:,1:1+nup+ndn,:), &
       qlsink_cen_sum(:,:,1:1+nup+ndn,:), precr_cen_sum(:,:,1:1+nup+ndn,:), &
       precsolid_cen_sum(:,:,1:1+nup+ndn,:), precall_cen_sum(:,:,1:1+nup+ndn,:), &
       qlsink_bf_cen_sum(:,:,1:1+nup+ndn,:), qlsink_avg_cen_sum(:,:,1:1+nup+ndn,:), & 
       prain_cen_sum(:,:,1:1+nup+ndn,:)  )

   wup_thresh(:) = 0.0
   wdown_thresh(:) = 0.0

   ntavg1 = ntavg1_ss  / dt
   ntavg2 = ntavg2_ss  / dt
   itavg1 = 0
   itavg2 = 0

end subroutine ecpp_crm_init
!---------------------------------------------------------------------------------------

!=======================================================================================
subroutine ecpp_crm_cleanup ()

! deallocate variables
  deallocate (updraftbase, &
            updrafttop )
  deallocate (dndraftbase, &
            dndrafttop )
 
  deallocate( qlsink, precr, precsolid, rh, qvs)

  deallocate( qlsink_bf, prain, qcloud_bf)

  deallocate( qcloudsum1, qcloud_bfsum1, qrainsum1, &
            qicesum1, qsnowsum1, qgraupsum1, &
            qlsinksum1, precrsum1, &
            precsolidsum1, precallsum1, &
            altsum1, rhsum1, cf3dsum1, &
            wwsum1, wwsqsum1, tkesgssum1, &
            qlsink_bfsum1, prainsum1, qvssum1 )

  deallocate(           &
            xkhvsum, wup_thresh, wdown_thresh )

  deallocate(wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum)

  deallocate( area_bnd_final, &
            area_bnd_sum, &
            area_cen_final, &
            area_cen_sum, &
            mass_bnd_final, &
            mass_bnd_sum, &
            mass_cen_final, &
            mass_cen_sum, &
            ent_bnd_sum, &
            rh_cen_sum, &
            qcloud_cen_sum, &
            qcloud_bf_cen_sum, &
            qrain_cen_sum, &
            qice_cen_sum, &
            qsnow_cen_sum, &
            qgraup_cen_sum, &
            qlsink_cen_sum, &
            precr_cen_sum, &
            precsolid_cen_sum, &
            precall_cen_sum, &
            qlsink_bf_cen_sum, &
            qlsink_avg_cen_sum, &
            prain_cen_sum  )

end subroutine ecpp_crm_cleanup
!---------------------------------------------------------------------------------------

!========================================================================================
subroutine ecpp_crm_stat()

  use crmx_module_ecpp_stats
  use module_data_ecpp1, only: afrac_cut 
  use crmx_grid,  only: nx, ny, nzm, pres
  use crmx_vars,  only: w, tabs, p, CF3D
  use crmx_sgs, only: tke, tk
  use crmx_microphysics, only: micro_field, iqv, iqci, iqr, iqs, iqg, cloudliq
  use crmx_module_mp_GRAUPEL, only: POLYSVP
#ifdef CLUBB_CRM
  use crmx_clubbvars, only: wp2
  use crmx_sgs, only: tk_clubb
#endif 
  implicit none

  integer :: i, ierr, i_tidx, j, &
      ncnt1, ncnt2

  integer :: nup, ndn
  integer :: kbase, ktop, m
  integer :: ii, jj, kk
  integer :: icl, icls, ipr

  real,dimension(nx, ny, nzm) :: &
       qcloud, qrain, qice, qsnow, qgraup, &
       precall, alt, xkhv
  real, dimension(nx, ny, nzstag) :: ww, wwsq

  real :: EVS

!------------------------------------------------------------------------
! Main code section...
!------------------------------------------------------------------------

     ndn = ndndraft ; nup = nupdraft

     itavg1 = itavg1 + 1
     itavg2 = itavg2 + 1
     ndn = ndndraft ; nup = nupdraft

! Get values from SAM cloud fields
     qcloud(1:nx,1:ny,1:nzm) = cloudliq(1:nx,1:ny,1:nzm)
     qrain(1:nx,1:ny,1:nzm)  = micro_field(1:nx,1:ny,1:nzm,iqr)
     qice(1:nx,1:ny,1:nzm)   = micro_field(1:nx,1:ny,1:nzm,iqci)
     qsnow(1:nx,1:ny,1:nzm)  = micro_field(1:nx,1:ny,1:nzm,iqs) 
     qgraup(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqg)

     precall(:,:,:)= precr(:,:,:) + precsolid(:,:,:)
    
     do ii=1, nx
      do jj=1, ny
       do kk=1, nzm
         EVS = POLYSVP(tabs(ii,jj,kk),0)   ! saturation water vapor pressure (PA)
         qvs(ii,jj,kk) = .622*EVS/(pres(kk)*100.-EVS)  ! pres(kk) with unit of hPa
!         rh(ii,jj,kk) = micro_field(ii,jj,kk,iqv)/QVS ! unit 0-1
!         rh(ii,jj,kk) = min(1.0, rh(ii,jj,kk))    ! RH is diagnosed in microphysics
         alt(ii,jj,kk) =  287.*tabs(ii,jj,kk)/(100.*pres(kk))

       end do
      end do
     end do
     
     ww(:,:,:)     = w(1:nx,1:ny,1:nzstag)
#ifdef CLUBB_CRM
     wwsq(:,:,:)  = sqrt(wp2(1:nx, 1:ny, 1:nzstag))
#else
     wwsq(:,:,:)   = 0.  ! subgrid vertical velocity is not used in the current version of ECPP. 
#endif

#ifdef CLUBB_CRM
     xkhv(:,:,:)   = tk_clubb(1:nx,1:ny,1:nzm)  ! eddy viscosity m2/s
#else
     xkhv(:,:,:)   = tk(1:nx,1:ny,1:nzm)  ! eddy viscosity m2/s
#endif

!+++mhwangtest
!     do ii=1, nx
!       do jj=1, ny
!        do kk=1, nzm
!          if(prain(ii,jj,kk).gt.1.0e-15) then
!            if(qcloud_bf(ii,jj,kk)*qlsink_bf(ii,jj,kk)/prain(ii,jj,kk)  .lt. 0.90) then
!             write(0, *) 'qcloud_bf*qlsink_bf/prain, qlsink_bf, qlsink, qlcoud_bf, qcloud, prain', qcloud_bf(ii,jj,kk)*qlsink_bf(ii,jj,kk)/prain(ii,jj,kk),  &
!                 qlsink_bf(ii, jj, kk) * 86400, qlsink(ii, jj, kk)*86400, qcloud_bf(ii, jj, kk), qcloud(ii, jj, kk), prain(ii, jj, kk) 
!            end if
!          end if
!        end do
!       end do
!     end do
!---mhwangest     


! Increment the 3-D running sums for averaging period 1.
     call rsums1( qcloud,    qcloudsum1(:,:,:),    &
                  qcloud_bf,    qcloud_bfsum1(:,:,:),    &
                  qrain,     qrainsum1(:,:,:),     &
                  qice,      qicesum1(:,:,:),      &
                  qsnow,     qsnowsum1(:,:,:),     &
                  qgraup,    qgraupsum1(:,:,:),    &
                  qlsink,    qlsinksum1(:,:,:),    &
                  precr,     precrsum1(:,:,:),     &
                  precsolid, precsolidsum1(:,:,:), &
                  precall,   precallsum1(:,:,:),   &
                  alt,       altsum1(:,:,:),       &
                  rh,        rhsum1(:,:,:),        &
                  CF3D,      cf3dsum1(:,:,:),       & 
                  ww,        wwsum1(:,:,:),        &
                  wwsq,      wwsqsum1(:,:,:),      &
                  tke(1:nx,1:ny,1:nzm),       tkesgssum1(:,:,:),    & 
                  qlsink_bf, qlsink_bfsum1(:,:,:), &
                  prain,     prainsum1(:,:,:),     &
                  qvs,       qvssum1(:,:,:)   )

! Increment the running sums for the level two variables that are not
! already incremented. Consolidate from 3-D to 1-D columns.
     call rsums2( &
          nx, ny, nzm, &
          xkhv,      xkhvsum(:) )

! Check if we have reached the end of the level 1 time averaging period.
     if( mod(itavg1,ntavg1) == 0 ) then

! Turn the running sums into averages.
        if( itavg1 /= 0 ) then
           ncnt1 = ntavg1
        else
           ncnt1 = 1
        end if
        call rsums1ToAvg( ncnt1, qcloudsum1(:,:,:), qcloud_bfsum1(:,:,:), qrainsum1(:,:,:), &
                                 qicesum1(:,:,:), qsnowsum1(:,:,:), &
                                 qgraupsum1(:,:,:), &
                                 qlsinksum1(:,:,:), precrsum1(:,:,:), &
                                 precsolidsum1(:,:,:), precallsum1(:,:,:), &
                                 altsum1(:,:,:), rhsum1(:,:,:), cf3dsum1(:,:,:), &
                                 wwsum1(:,:,:), wwsqsum1(:,:,:), &
                                 tkesgssum1(:,:,:), qlsink_bfsum1(:,:,:), &
                                 prainsum1(:,:,:), qvssum1(:,:,:)  )

! Determine draft categories and get running sums of them.
        call categorization_stats( .true., &
             nx, ny, nzm, nupdraft, ndndraft, ndraft_max, &
             mode_updnthresh, upthresh, downthresh, &
             upthresh2, downthresh2, cloudthresh, prcpthresh, &
             cloudthresh_trans, precthresh_trans,  &
             qvssum1(:,:,:),          & 
             plumetype, allcomb, &
             updraftbase(1:nupdraft), updrafttop(1:nupdraft), &
             dndraftbase(1:ndndraft), dndrafttop(1:ndndraft), &
             qcloudsum1(:,:,:), qcloud_bfsum1(:,:,:), qrainsum1(:,:,:), &
             qicesum1(:,:,:), qsnowsum1(:,:,:), qgraupsum1(:,:,:), &
             qlsinksum1(:,:,:), precrsum1(:,:,:), &
             precsolidsum1(:,:,:), precallsum1(:,:,:), &
             altsum1(:,:,:), rhsum1(:,:,:), cf3dsum1(:,:,:), &
             wwsum1(:,:,:), wwsqsum1(:,:,:), tkesgssum1(:,:,:),  &
             qlsink_bfsum1(:,:,:), prainsum1(:,:,:), &
             area_bnd_final(:,:,1:1+ndn+nup,:), area_cen_final(:,:,1:1+ndn+nup,:), &
             area_bnd_sum(:,:,1:1+ndn+nup,:), area_cen_sum(:,:,1:1+ndn+nup,:), &
             ent_bnd_sum(:,:,1:1+ndn+nup,:), mass_bnd_sum(:,:,1:1+ndn+nup,:), &
             rh_cen_sum(:,:,1:1+ndn+nup,:), &
             qcloud_cen_sum(:,:,1:1+ndn+nup,:), qcloud_bf_cen_sum(:,:,1:1+ndn+nup,:), qrain_cen_sum(:,:,1:1+ndn+nup,:), &
             qice_cen_sum(:,:,1:1+ndn+nup,:), qsnow_cen_sum(:,:,1:1+ndn+nup,:), &
             qgraup_cen_sum(:,:,1:1+ndn+nup,:), &
             qlsink_cen_sum(:,:,1:1+ndn+nup,:), precr_cen_sum(:,:,1:1+ndn+nup,:), &
             precsolid_cen_sum(:,:,1:1+nup+ndn,:), precall_cen_sum(:,:,1:1+nup+ndn,:),  &
             qlsink_bf_cen_sum(:,:,1:1+nup+ndn,:), prain_cen_sum(:,:,1:1+nup+ndn,:),  &
             wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum, &  
             wup_thresh, wdown_thresh )

! If we want final area categories based on the last avg1 period in each
! avg2 then we need to zero out the running sum just created for the areas
! if it is not the last block of time in ntavg2
        if( areaavgtype==1 .and. .not. mod(itavg2,ntavg2)==0 ) then
           call zero_out_areas( &
                area_bnd_final(:,:,1:1+ndn+nup,:), &
                area_cen_final(:,:,1:1+ndn+nup,:) )
        end if

! Done with time level one averages so zero them out for next period.
        call zero_out_sums1( qcloudsum1(:,:,:), qcloud_bfsum1(:,:,:), qrainsum1(:,:,:), &
             qicesum1(:,:,:), qsnowsum1(:,:,:), qgraupsum1(:,:,:), &
             qlsinksum1(:,:,:), precrsum1(:,:,:), &
             precsolidsum1(:,:,:), precallsum1(:,:,:), &
             altsum1(:,:,:), rhsum1(:,:,:), cf3dsum1(:,:,:), &
             wwsum1(:,:,:), wwsqsum1(:,:,:), tkesgssum1(:,:,:), &
             qlsink_bfsum1(:,:,:), prainsum1(:,:,:), qvssum1(:,:,:) )

     end if !End of time level one averaging period

! Check if we have reached the end of a level 2 averaging period.
     if( mod(itavg2,ntavg2) == 0 ) then

! Turn the running sums into averages. ncnt1 in this case is the number
! of calls to categorization_stats during the level 2 averaging period,
! which increment the bnd/cen arrays.
        if( itavg2 /= 0 ) then
           ncnt1 = ntavg2_ss/ntavg1_ss
           ncnt2 = ntavg2
        else
           ncnt1 = 1
           ncnt2 = 1
        end if

        call rsums2ToAvg( areaavgtype, nx, ny, ncnt1, ncnt2, &
             xkhvsum(:), &
             wwqui_cen_sum(:), wwqui_bnd_sum(:), wwqui_cloudy_cen_sum(:), wwqui_cloudy_bnd_sum(:),  &
             area_bnd_final(:,:,1:1+ndn+nup,:), area_bnd_sum(:,:,1:1+ndn+nup,:), &
             area_cen_final(:,:,1:1+ndn+nup,:), area_cen_sum(:,:,1:1+ndn+nup,:), &
             mass_bnd_final(:,:,1:1+ndn+nup,:), mass_bnd_sum(:,:,1:1+ndn+nup,:), &
             mass_cen_final(:,:,1:1+ndn+nup,:), mass_cen_sum(:,:,1:1+ndn+nup,:), &
             ent_bnd_sum(:,:,1:1+ndn+nup,:), &
             rh_cen_sum(:,:,1:1+ndn+nup,:), &
             qcloud_cen_sum(:,:,1:1+ndn+nup,:), qcloud_bf_cen_sum(:,:,1:1+ndn+nup,:), qrain_cen_sum(:,:,1:1+ndn+nup,:), &
             qice_cen_sum(:,:,1:1+ndn+nup,:), qsnow_cen_sum(:,:,1:1+ndn+nup,:), &
             qgraup_cen_sum(:,:,1:1+ndn+nup,:), &
             qlsink_cen_sum(:,:,1:1+ndn+nup,:), precr_cen_sum(:,:,1:1+ndn+nup,:), &
             precsolid_cen_sum(:,:,1:1+ndn+nup,:), precall_cen_sum(:,:,1:1+ndn+nup,:), &
             qlsink_bf_cen_sum(:,:,1:1+ndn+nup,:), prain_cen_sum(:,:,1:1+ndn+nup,:) )

! get in-cloud value for rh, qcloud, qrain, qice, qsnow, qgraup, 
! percr, precsolid, and precall. (qlsink is already in-cloud values)
        do kk=1, nzm
         do icl=1, NCLASS_CL
           do icls=1, ncls_ecpp_in
             do ipr=1, NCLASS_PR
               if(area_cen_sum(kk, icl, icls, ipr).gt.afrac_cut) then
                 rh_cen_sum(kk,icl,icls,ipr) = rh_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 qcloud_cen_sum(kk,icl,icls,ipr) = qcloud_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 qcloud_bf_cen_sum(kk,icl,icls,ipr) = qcloud_bf_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 qrain_cen_sum(kk,icl,icls,ipr) = qrain_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 qice_cen_sum(kk,icl,icls,ipr) = qice_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 qsnow_cen_sum(kk,icl,icls,ipr) = qsnow_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 qgraup_cen_sum(kk,icl,icls,ipr) = qgraup_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 precr_cen_sum(kk,icl,icls,ipr) = precr_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 precsolid_cen_sum(kk,icl,icls,ipr) = precsolid_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 precall_cen_sum(kk,icl,icls,ipr) = precall_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 prain_cen_sum(kk,icl,icls,ipr) = prain_cen_sum(kk,icl,icls,ipr)/area_cen_sum(kk,icl,icls,ipr)
                 if(qcloud_bf_cen_sum(kk,icl,icls,ipr).gt.1.0e-10) then
                   qlsink_avg_cen_sum(kk,icl,icls,ipr) = min(1.0/ntavg2_ss, &
                                                         prain_cen_sum(kk,icl,icls,ipr)/qcloud_bf_cen_sum(kk,icl,icls,ipr))
                 else
                   qlsink_avg_cen_sum(kk,icl,icls,ipr) = 0.0
                 end if
                 qlsink_bf_cen_sum(kk,icl,icls,ipr) = min(1.0/ntavg2_ss, qlsink_bf_cen_sum(kk,icl,icls,ipr))
                 qlsink_cen_sum(kk,icl,icls,ipr) = min(1.0/ntavg2_ss, qlsink_cen_sum(kk,icl,icls,ipr))
               else
                 rh_cen_sum(kk,icl,icls,ipr) = 0.0
                 qcloud_cen_sum(kk,icl,icls,ipr) = 0.0
                 qcloud_bf_cen_sum(kk,icl,icls,ipr) = 0.0
                 qrain_cen_sum(kk,icl,icls,ipr) = 0.0
                 qice_cen_sum(kk,icl,icls,ipr) = 0.0
                 qsnow_cen_sum(kk,icl,icls,ipr) = 0.0
                 qgraup_cen_sum(kk,icl,icls,ipr) = 0.0
                 precr_cen_sum(kk,icl,icls,ipr) = 0.0
                 precsolid_cen_sum(kk,icl,icls,ipr) = 0.0
                 precall_cen_sum(kk,icl,icls,ipr) = 0.0
                 qlsink_bf_cen_sum(kk,icl,icls,ipr) = 0.0
                 prain_cen_sum(kk,icl,icls,ipr) = 0.0
                 qlsink_avg_cen_sum(kk,icl,icls,ipr) = 0.0
                 qlsink_bf_cen_sum(kk,icl,icls,ipr) = 0.0
                 qlsink_cen_sum(kk,icl,icls,ipr) = 0.0
               end if
             end do
           end do
         end do
!
! calculate vertical velocity variance for quiescent class
          if(sum(area_cen_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR)).gt.afrac_cut) then  
             wwqui_cen_sum(kk) = wwqui_cen_sum(kk) / sum(area_cen_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR))
          else
             wwqui_cen_sum(kk) = 0.0
          end if
          if(sum(area_cen_sum(kk,CLD, QUI, 1:NCLASS_PR)).gt.afrac_cut) then
             wwqui_cloudy_cen_sum(kk) = wwqui_cloudy_cen_sum(kk) / sum(area_cen_sum(kk, CLD, QUI, 1:NCLASS_PR))
          else
             wwqui_cloudy_cen_sum(kk) = 0.0
          end if
     
        end do  ! kk
!
! calcualte vertical velocity variance for quiescent calss at lay boundary
        do kk=1, nzm+1
          if(sum(area_bnd_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR)).gt.afrac_cut) then
             wwqui_bnd_sum(kk) = wwqui_bnd_sum(kk) / sum(area_bnd_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR))
          else
             wwqui_bnd_sum(kk) = 0.0
          end if
          if(sum(area_bnd_sum(kk,CLD, QUI, 1:NCLASS_PR)).gt.afrac_cut) then
             wwqui_cloudy_bnd_sum(kk) = wwqui_cloudy_bnd_sum(kk) / sum(area_bnd_sum(kk, CLD, QUI, 1:NCLASS_PR))
          else
             wwqui_cloudy_bnd_sum(kk) = 0.0
          end if
        end do

     end if !End of level two time averaging period

end subroutine ecpp_crm_stat

#endif /*ECPP*/
end module  crmx_module_ecpp_crm_driver
