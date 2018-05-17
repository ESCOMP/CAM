module crmx_crm_module
!---------------------------------------------------------------
!  Super-parameterization's main driver 
!  Marat Khairoutdinov, 2001-2009
!---------------------------------------------------------------

use crmx_setparm_mod, only : setparm

contains

subroutine crm        (lchnk, icol, &
                       tl, ql, qccl, qiil, ul, vl, &
                       ps, pmid, pdel, phis, &
                       zmid, zint, dt_gl, plev, &
                       qltend, qcltend, qiltend, sltend, &
                       u_crm, v_crm, w_crm, t_crm, micro_fields_crm, &
                       qrad_crm, &
                       qc_crm, qi_crm, qpc_crm, qpi_crm, prec_crm, &
                       t_rad, qv_rad, qc_rad, qi_rad, cld_rad, cld3d_crm, &
#ifdef m2005
                       nc_rad, ni_rad, qs_rad, ns_rad, wvar_crm,  &
! hm 7/26/11 new output
                       aut_crm, acc_crm, evpc_crm, evpr_crm, mlt_crm, &
                       sub_crm, dep_crm, con_crm, &
! hm 8/31/11 new output for gcm-grid and time-step avg process rates
                       aut_crm_a, acc_crm_a, evpc_crm_a, evpr_crm_a, mlt_crm_a, &
                       sub_crm_a, dep_crm_a, con_crm_a, &
#endif
                       precc, precl, precsc, precsl, &
                       cltot, clhgh, clmed, cllow, cld, cldtop, &
                       gicewp, gliqwp, &
                       mc, mcup, mcdn, mcuup, mcudn, &
                       crm_qc, crm_qi, crm_qs, crm_qg, crm_qr, &
#ifdef m2005
                       crm_nc, crm_ni, crm_ns, crm_ng, crm_nr, &
#ifdef MODAL_AERO
                       naermod, vaerosol, hygro,     &
#endif 
#endif
#ifdef SPCAM_CLUBB_SGS
                       clubb_buffer,                 &
                       crm_cld,                      &
                       clubb_tk, clubb_tkh,          &
                       relvar, accre_enhan, qclvar,  &
#endif
                       crm_tk, crm_tkh,              &
                       mu_crm, md_crm, du_crm, eu_crm, ed_crm, jt_crm, mx_crm,    &
#ifdef ECPP
                       abnd, abnd_tf, massflxbnd, acen, acen_tf,           &
                       rhcen, qcloudcen, qicecen, qlsinkcen, precrcen, precsolidcen,  & 
                       qlsink_bfcen, qlsink_avgcen, praincen,     &
                       wupthresh_bnd, wdownthresh_bnd,   &
                       wwqui_cen, wwqui_bnd, wwqui_cloudy_cen, wwqui_cloudy_bnd,   &
#endif
                       tkez, tkesgsz, tkz, flux_u, flux_v, flux_qt, fluxsgs_qt,flux_qp, &
                       pflx, qt_ls, qt_trans, qp_trans, qp_fall, &
                       qp_evp, qp_src, t_ls, prectend, precstend, &
                       ocnfrac, wndls, tau00, bflxls, &
                       fluxu00, fluxv00, fluxt00, fluxq00,    &
                       taux_crm, tauy_crm, z0m, timing_factor, qtot)   

!            dolong, doshort, nrad0, &
!            latitude00, longitude00, day00, pres00, tabs_s0, case0, &
!            radlwup0, radlwdn0, radswup0, radswdn0, radqrlw0, radqrsw0, &
!            lwnsxy,swnsxy,lwntxy,swntxy,solinxy,lwnscxy,swnscxy,lwntcxy,swntcxy,lwdsxy,swdsxy)


!---------------------------------------------------------------

        use shr_kind_mod, only: r8 => shr_kind_r8
#ifdef SPCAM_CLUBB_SGS
        use crmdims, only: nclubbvars
#endif
        use phys_grid, only: get_rlon_p, get_rlat_p, get_gcol_all_p
        use ppgrid, only: pcols
        use crmx_vars
        use crmx_params
        use crmx_microphysics
        use crmx_sgs
        use crmx_crmtracers
#ifdef MODAL_AERO
        use modal_aero_data,   only: ntot_amode
#endif
#ifdef SPCAM_CLUBB_SGS
        use crmx_clubb_sgs, only: advance_clubb_sgs, clubb_sgs_setup, clubb_sgs_cleanup, &
	       apply_clubb_sgs_tndcy, apply_clubb_sgs_tndcy_scalars, apply_clubb_sgs_tndcy_mom,   & ! Subroutines
	       t2thetal                 ! Functions 
        use crmx_clubb_sgs, only: total_energy
	    use crmx_clubbvars, only: edsclr_dim, sclr_dim, rho_ds_zt, rho_ds_zm, &
	        rtm_spurious_source, thlm_spurious_source
        use crmx_clubb_precision, only: time_precision
        use crmx_clubbvars,  only:  up2, vp2, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp, &
                             upwp, vpwp, cloud_frac, t_tndcy, qc_tndcy, qv_tndcy, u_tndcy, v_tndcy, lrestart_clubb  
        use crmx_clubbvars, only: rho_ds_zt, rho_ds_zm, thv_ds_zt, thv_ds_zm, &
         invrs_rho_ds_zt, invrs_rho_ds_zm
	use crmx_clubbvars, only: tracer_tndcy, sclrp2, sclrprtp, sclrpthlp, wpsclrp
	use crmx_fill_holes, only: vertical_integral ! Function
	use crmx_numerical_check, only: calculate_spurious_source
	use crmx_grid_class, only: gr   ! Variable
       use crmx_clubb_precision, only: core_rknd ! Constants
       use crmx_clubbvars, only: relvarg, accre_enhang, qclvarg
#endif /*CLUBB_SGS*/
#ifdef ECPP
        use crmx_ecppvars, only: qlsink, precr, precsolid, &
                    area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
                    mass_bnd_final, mass_bnd_sum, rh_cen_sum, qcloud_cen_sum, qice_cen_sum, &
                    qlsink_cen_sum, precr_cen_sum, precsolid_cen_sum, xkhvsum, wup_thresh, wdown_thresh, &
                    wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum, &
                    qlsink_bf_cen_sum, qlsink_avg_cen_sum, prain_cen_sum, qlsink_bf, prain
        use crmx_module_ecpp_crm_driver, only: ecpp_crm_stat, ecpp_crm_init, ecpp_crm_cleanup, ntavg1_ss, ntavg2_ss
        use crmx_ecppvars,  only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif /*ECPP*/

        use cam_abortutils,  only: endrun
        use time_manager,    only: get_nstep

        implicit none

!        integer, parameter :: r8 = 8

!  Input:

         integer, intent(in) :: lchnk    ! chunk identifier
         integer, intent(in) :: icol     ! column identifier
         integer, intent(in) :: plev     ! number of levels
         real(r8), intent(in) :: ps ! Global grid surface pressure (Pa)
         real(r8), intent(in) :: pmid(plev) ! Global grid pressure (Pa)
         real(r8), intent(in) :: pdel(plev) ! Layer's pressure thickness (Pa)
         real(r8), intent(in) :: phis ! Global grid surface geopotential (m2/s2)
         real(r8), intent(in) :: zmid(plev) ! Global grid height (m)
         real(r8), intent(in) :: zint(plev+1)! Global grid interface height (m)
         real(r8), intent(in) :: qrad_crm(crm_nx, crm_ny, crm_nz) ! CRM rad. heating
         real(r8), intent(in) :: dt_gl ! global model's time step
         real(r8), intent(in) :: ocnfrac ! area fraction of the ocean
         real(r8), intent(in) :: tau00  ! large-scale surface stress (N/m2)
         real(r8), intent(in) :: wndls  ! large-scale surface wind (m/s)
         real(r8), intent(in) :: bflxls  ! large-scale surface buoyancy flux (K m/s)
         real(r8), intent(in) :: fluxu00  ! surface momenent fluxes [N/m2]
         real(r8), intent(in) :: fluxv00  ! surface momenent fluxes [N/m2]
         real(r8), intent(in) :: fluxt00  ! surface sensible heat fluxes [K Kg/ (m2 s)]
         real(r8), intent(in) :: fluxq00  ! surface latent heat fluxes [ kg/(m2 s)]
!         logical, intent(in)  :: doshort ! compute shortwave radiation
!         logical, intent(in)  :: dolong ! compute longwave radiation
!         real(r8), intent(in) :: day00 ! initial day
!         real(r8), intent(in) :: latitude00
!         real(r8), intent(in) :: longitude00
!         real(r8), intent(in) :: pres00
!         real(r8), intent(in) :: tabs_s0
!         integer , intent(in) :: nrad0
!         character *40 case0  ! 8-symbol id-string to identify a case-name


! tl, ql, qccl, qiil, ul, vl are not updated in this subroutine, and set to intent(in), but
! not intent(inout). +++mhwang 
         real(r8), intent(in) :: tl(plev) ! Global grid temperature (K)
         real(r8), intent(in) :: ql(plev) ! Global grid water vapor (g/g)
         real(r8), intent(in) :: qccl(plev)! Global grid cloud liquid water (g/g)
         real(r8), intent(in) :: qiil(plev)! Global grid cloud ice (g/g)
         real(r8), intent(in) :: ul(plev) ! Global grid u (m/s)
         real(r8), intent(in) :: vl(plev) ! Global grid v (m/s)

!  Input/Output:
#ifdef SPCAM_CLUBB_SGS
         real(r8), intent(inout), target :: clubb_buffer(crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
         real(r8), intent(inout)  :: crm_cld(crm_nx, crm_ny, crm_nz+1)
         real(r8), intent(inout)  :: clubb_tk(crm_nx, crm_ny, crm_nz)
         real(r8), intent(inout)  :: clubb_tkh(crm_nx, crm_ny, crm_nz)
         real(r8), intent(inout)  :: relvar(crm_nx, crm_ny, crm_nz) 
         real(r8), intent(inout)  :: accre_enhan(crm_nx, crm_ny, crm_nz)
         real(r8), intent(inout)  :: qclvar(crm_nx, crm_ny, crm_nz)
#endif
         real(r8), intent(inout)  :: crm_tk(crm_nx, crm_ny, crm_nz)
         real(r8), intent(inout)  :: crm_tkh(crm_nx, crm_ny, crm_nz)

         real(r8), intent(inout) :: cltot ! shaded cloud fraction
         real(r8), intent(inout) :: clhgh ! shaded cloud fraction
         real(r8), intent(inout) :: clmed ! shaded cloud fraction
         real(r8), intent(inout) :: cllow ! shaded cloud fraction

         
!  Output
         
         real(r8), intent(inout) :: sltend(plev) ! tendency of static energy
!         real(r8), intent(inout) :: u_crm  (:,:,:) ! CRM v-wind component
!         real(r8), intent(inout) :: v_crm  (:,:,:) ! CRM v-wind component
!         real(r8), intent(inout) :: w_crm  (:,:,:) ! CRM w-wind component
!         real(r8), intent(inout) :: t_crm  (:,:,:) ! CRM temperuture
         real(r8), intent(inout) :: u_crm  (crm_nx,crm_ny,crm_nz) ! CRM v-wind component
         real(r8), intent(inout) :: v_crm  (crm_nx,crm_ny,crm_nz) ! CRM v-wind component
         real(r8), intent(inout) :: w_crm  (crm_nx,crm_ny,crm_nz) ! CRM w-wind component
         real(r8), intent(inout) :: t_crm  (crm_nx,crm_ny,crm_nz) ! CRM temperuture
!         real(r8), intent(inout) :: micro_fields_crm  (:,:,:,:) ! CRM total water
         real(r8), intent(inout) :: micro_fields_crm  (crm_nx,crm_ny,crm_nz,nmicro_fields+1) ! CRM total water
         real(r8), intent(inout) :: qltend(plev) ! tendency of water vapor
         real(r8), intent(inout) :: qcltend(plev)! tendency of cloud liquid water
         real(r8), intent(inout) :: qiltend(plev)! tendency of cloud ice
         real(r8), intent(inout) :: t_rad (crm_nx, crm_ny, crm_nz) ! rad temperuture
         real(r8), intent(inout) :: qv_rad(crm_nx, crm_ny, crm_nz) ! rad vapor
         real(r8), intent(inout) :: qc_rad(crm_nx, crm_ny, crm_nz) ! rad cloud water
         real(r8), intent(inout) :: qi_rad(crm_nx, crm_ny, crm_nz) ! rad cloud ice
         real(r8), intent(inout) :: cld_rad(crm_nx, crm_ny, crm_nz) ! rad cloud fraction 
         real(r8), intent(inout) :: cld3d_crm(crm_nx, crm_ny, crm_nz) ! instant 3D cloud fraction
#ifdef m2005
         real(r8), intent(inout) :: nc_rad(crm_nx, crm_ny, crm_nz) ! rad cloud droplet number (#/kg) 
         real(r8), intent(inout) :: ni_rad(crm_nx, crm_ny, crm_nz) ! rad cloud ice crystal number (#/kg)
         real(r8), intent(inout) :: qs_rad(crm_nx, crm_ny, crm_nz) ! rad cloud snow (kg/kg)
         real(r8), intent(inout) :: ns_rad(crm_nx, crm_ny, crm_nz) ! rad cloud snow crystal number (#/kg)
         real(r8), intent(inout) :: wvar_crm(crm_nx, crm_ny, crm_nz) ! vertical velocity variance (m/s)
! hm 7/26/11 new output
         real(r8), intent(inout) :: aut_crm(crm_nx, crm_ny, crm_nz) ! cloud water autoconversion (1/s)
         real(r8), intent(inout) :: acc_crm(crm_nx, crm_ny, crm_nz) ! cloud water accretion (1/s)
         real(r8), intent(inout) :: evpc_crm(crm_nx, crm_ny, crm_nz) ! cloud water evaporation (1/s)
         real(r8), intent(inout) :: evpr_crm(crm_nx, crm_ny, crm_nz) ! rain evaporation (1/s)
         real(r8), intent(inout) :: mlt_crm(crm_nx, crm_ny, crm_nz) ! ice, snow, graupel melting (1/s)
         real(r8), intent(inout) :: sub_crm(crm_nx, crm_ny, crm_nz) ! ice, snow, graupel sublimation (1/s)
         real(r8), intent(inout) :: dep_crm(crm_nx, crm_ny, crm_nz) ! ice, snow, graupel deposition (1/s)
         real(r8), intent(inout) :: con_crm(crm_nx, crm_ny, crm_nz) ! cloud water condensation(1/s)
! hm 8/31/11 new output, gcm-grid and time step-avg
         real(r8), intent(inout) :: aut_crm_a(plev) ! cloud water autoconversion (1/s)
         real(r8), intent(inout) :: acc_crm_a(plev) ! cloud water accretion (1/s)
         real(r8), intent(inout) :: evpc_crm_a(plev) ! cloud water evaporation (1/s)
         real(r8), intent(inout) :: evpr_crm_a(plev) ! rain evaporation (1/s)
         real(r8), intent(inout) :: mlt_crm_a(plev) ! ice, snow, graupel melting (1/s)
         real(r8), intent(inout) :: sub_crm_a(plev) ! ice, snow, graupel sublimation (1/s)
         real(r8), intent(inout) :: dep_crm_a(plev) ! ice, snow, graupel deposition (1/s)
         real(r8), intent(inout) :: con_crm_a(plev) ! cloud water condensation(1/s)
#endif
         real(r8), intent(inout) :: precc ! convective precip rate (m/s)
         real(r8), intent(inout) :: precl ! stratiform precip rate (m/s)
         real(r8), intent(inout) :: cld(plev)  ! cloud fraction
         real(r8), intent(inout) :: cldtop(plev)  ! cloud top pdf
         real(r8), intent(inout) :: gicewp(plev)  ! ice water path
         real(r8), intent(inout) :: gliqwp(plev)  ! ice water path
         real(r8), intent(inout) :: mc(plev)   ! cloud mass flux
         real(r8), intent(inout) :: mcup(plev) ! updraft cloud mass flux
         real(r8), intent(inout) :: mcdn(plev) ! downdraft cloud mass flux
         real(r8), intent(inout) :: mcuup(plev) ! unsat updraft cloud mass flux
         real(r8), intent(inout) :: mcudn(plev) ! unsat downdraft cloud mass flux
         real(r8), intent(inout) :: crm_qc(plev)  ! mean cloud water
         real(r8), intent(inout) :: crm_qi(plev)  ! mean cloud ice
         real(r8), intent(inout) :: crm_qs(plev)  ! mean snow
         real(r8), intent(inout) :: crm_qg(plev)  ! mean graupel
         real(r8), intent(inout) :: crm_qr(plev)  ! mean rain
#ifdef m2005
         real(r8), intent(inout) :: crm_nc(plev)  ! mean cloud water  (#/kg)
         real(r8), intent(inout) :: crm_ni(plev)  ! mean cloud ice    (#/kg)
         real(r8), intent(inout) :: crm_ns(plev)  ! mean snow         (#/kg)
         real(r8), intent(inout) :: crm_ng(plev)  ! mean graupel      (#/kg)
         real(r8), intent(inout) :: crm_nr(plev)  ! mean rain         (#/kg)
#ifdef MODAL_AERO
         real(r8), intent(in)  :: naermod(plev, ntot_amode)     ! Aerosol number concentration [/m3]
         real(r8), intent(in)  :: vaerosol(plev, ntot_amode)    ! aerosol volume concentration [m3/m3]
         real(r8), intent(in)  :: hygro(plev, ntot_amode)       ! hygroscopicity of aerosol mode 
#endif 
#endif
         real(r8), intent(inout) :: mu_crm (plev)             ! mass flux up
         real(r8), intent(inout) :: md_crm (plev)             ! mass flux down
         real(r8), intent(inout) :: du_crm (plev)             ! mass detrainment from updraft
         real(r8), intent(inout) :: eu_crm (plev)             ! mass entrainment from updraft
         real(r8), intent(inout) :: ed_crm (plev)             ! mass detrainment from downdraft
         real(r8)              :: dd_crm (plev)             ! mass entraiment from downdraft
         real(r8), intent(inout) :: jt_crm                    ! index of cloud (convection) top 
         real(r8), intent(inout) :: mx_crm                    ! index of cloud (convection) bottom
         real(r8)              :: mui_crm (plev+1)             ! mass flux up at the interface
         real(r8)              :: mdi_crm (plev+1)             ! mass flux down at the interface

         real(r8), intent(inout) :: flux_qt(plev) ! nonprecipitating water flux           [kg/m2/s]
         real(r8), intent(inout) :: fluxsgs_qt(plev) ! sgs nonprecipitating water flux    [kg/m2/s]
         real(r8), intent(inout) :: tkez(plev) ! tke profile               [kg/m/s2]
         real(r8), intent(inout) :: tkesgsz(plev) ! sgs tke profile        [kg/m/s2]
         real(r8), intent(inout) :: tkz(plev)  ! tk profile                [m2/s]
         real(r8), intent(inout) :: flux_u(plev) ! x-momentum flux          [m2/s2]
         real(r8), intent(inout) :: flux_v(plev) ! y-momentum flux          [m2/s2]
         real(r8), intent(inout) :: flux_qp(plev) ! precipitating water flux [kg/m2/s or mm/s]
         real(r8), intent(inout) :: pflx(plev)    ! precipitation flux      [m/s]
         real(r8), intent(inout) :: qt_ls(plev) ! tendency of nonprec water due to large-scale  [kg/kg/s]
         real(r8), intent(inout) :: qt_trans(plev)! tendency of nonprec water due to transport  [kg/kg/s]
         real(r8), intent(inout) :: qp_trans(plev) ! tendency of prec water due to transport [kg/kg/s]
         real(r8), intent(inout) :: qp_fall(plev) ! tendency of prec water due to fall-out   [kg/kg/s]
         real(r8), intent(inout) :: qp_src(plev) ! tendency of prec water due to conversion  [kg/kg/s]
         real(r8), intent(inout) :: qp_evp(plev) ! tendency of prec water due to evp         [kg/kg/s]
         real(r8), intent(inout) :: t_ls(plev) ! tendency of lwse  due to large-scale        [kg/kg/s] ???
         real(r8), intent(inout) :: prectend ! column integrated tendency in precipitating water+ice (kg/m2/s)
         real(r8), intent(inout) :: precstend ! column integrated tendency in precipitating ice (kg/m2/s)
         real(r8), intent(inout) :: precsc ! convective snow rate (m/s)
         real(r8), intent(inout) :: precsl ! stratiform snow rate (m/s)
         real(r8), intent(inout):: taux_crm  ! zonal CRM surface stress perturbation (N/m2)
         real(r8), intent(inout):: tauy_crm  ! merid CRM surface stress perturbation (N/m2)
         real(r8), intent(inout):: z0m ! surface stress (N/m2)
         real(r8), intent(inout):: timing_factor ! crm cpu efficiency
         real(r8), intent(inout) :: qc_crm (crm_nx, crm_ny, crm_nz)! CRM cloud water
         real(r8), intent(inout) :: qi_crm (crm_nx, crm_ny, crm_nz)! CRM cloud ice
         real(r8), intent(inout) :: qpc_crm(crm_nx, crm_ny, crm_nz)! CRM precip water
         real(r8), intent(inout) :: qpi_crm(crm_nx, crm_ny, crm_nz)! CRM precip ice
         real(r8), intent(inout) :: prec_crm(crm_nx, crm_ny)! CRM precipiation rate
#ifdef ECPP
! at layer center
         real(r8), intent(inout) :: acen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! cloud fraction for each sub-sub class for full time period
         real(r8), intent(inout) :: acen_tf(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! cloud fraction for end-portion of time period
         real(r8), intent(inout) :: rhcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! relative humidity (0-1)
         real(r8), intent(inout) :: qcloudcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water (kg/kg)
         real(r8), intent(inout) :: qicecen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! cloud ice (kg/kg)
         real(r8), intent(inout) :: qlsinkcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (/s??)
         real(r8), intent(inout) :: precrcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! liquid (rain) precipitation rate (kg/m2/s)
         real(r8), intent(inout) :: precsolidcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! solid (rain) precipitation rate (kg/m2/s)
         real(r8), intent(inout) :: qlsink_bfcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated 
                                                                                      ! cloud water before precipitatinog (/s)
         real(r8), intent(inout) :: qlsink_avgcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated 
                                                                                      ! from praincen and qlcoudcen averaged over 
                                                                                      ! ntavg1_ss time step (/s??)
         real(r8), intent(inout) :: praincen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (kg/kg/s)
         real(r8), intent(inout) :: wwqui_cen(plev)                                ! vertical velocity variance in quiescent class (m2/s2)
         real(r8), intent(inout) :: wwqui_cloudy_cen(plev)                         ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
! at layer boundary
         real(r8), intent(inout) :: abnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! cloud fraction for each sub-sub class for full time period
         real(r8), intent(inout) :: abnd_tf(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! cloud fraction for end-portion of time period
         real(r8), intent(inout) :: massflxbnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
         real(r8), intent(inout) :: wupthresh_bnd(plev+1)             ! vertical velocity threshold for updraft (m/s)
         real(r8), intent(inout) :: wdownthresh_bnd(plev+1)           ! vertical velocity threshold for downdraft (m/s)
         real(r8), intent(inout) :: wwqui_bnd(plev+1)                                ! vertical velocity variance in quiescent class (m2/s2)
         real(r8), intent(inout) :: wwqui_cloudy_bnd(plev+1)                         ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
#endif

!  Local space:
        real dummy(nz), t00(nz)
        real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
        real tln(plev), qln(plev), qccln(plev), qiiln(plev), uln(plev), vln(plev)
        real cwp(nx,ny), cwph(nx,ny), cwpm(nx,ny), cwpl(nx,ny)
        real(r8) factor_xy, idt_gl
        real tmp1, tmp2
        real u2z,v2z,w2z
        integer i,j,k,l,ptop,nn,icyc, nstatsteps
        integer kx
        real(r8), parameter :: umax = 0.5*crm_dx/crm_dt ! maxumum ampitude of the l.s. wind
        real(r8), parameter :: wmin = 2.   ! minimum up/downdraft velocity for stat
        real, parameter :: cwp_threshold = 0.001 ! threshold for cloud condensate for shaded fraction calculation
        logical flag_top(nx,ny)
        real ustar, bflx, wnd, z0_est, qsat, omg
        real colprec,colprecs
        real(r8) zs ! surface elevation
        integer igstep    ! GCM time steps
        integer iseed   ! seed for random perturbation
        integer gcolindex(pcols)  ! array of global latitude indices

#ifdef SPCAM_CLUBB_SGS
!Array indicies for spurious RTM check

real(kind=core_rknd) :: &
  rtm_integral_before(nx,ny), rtm_integral_after(nx,ny), rtm_flux_top, rtm_flux_sfc
real(kind=core_rknd) :: &
  thlm_integral_before(nx,ny), thlm_integral_after(nx,ny), thlm_before(nzm), thlm_after(nzm), &
  thlm_flux_top, thlm_flux_sfc

real(kind=core_rknd), dimension(nzm) :: &
  rtm_column ! Total water (vapor + liquid)     [kg/kg]
#endif

        real  cltemp(nx,ny), cmtemp(nx,ny), chtemp(nx, ny), cttemp(nx, ny)

        real(r8), intent(inout) :: qtot(20)
        real ntotal_step

!-----------------------------------------------

        dostatis = .false.    ! no statistics are collected. 
        idt_gl = 1._r8/dt_gl
        ptop = plev-nzm+1
        factor_xy = 1._r8/dble(nx*ny)
        dummy = 0.
        t_rad = 0.
        qv_rad = 0.
        qc_rad = 0.
        qi_rad = 0.
        cld_rad = 0.
#ifdef m2005
        nc_rad = 0.0
        ni_rad = 0.0
        qs_rad = 0.0
        ns_rad = 0.0
#endif
        zs=phis/ggr
        bflx = bflxls
        wnd = wndls

!-----------------------------------------
        igstep = get_nstep() 

#ifdef SPCAM_CLUBB_SGS
        if(igstep == 1) then
          lrestart_clubb = .false.
        else
         lrestart_clubb = .true.
        endif
#endif

        call task_init ()

        call setparm()

!        doshortwave = doshort
!        dolongwave = dolong
!        day0 = day00-dt_gl/86400.
!        latitude = latitude00
!        longitude = longitude00
!        pres0 = pres00
!        tabs_s = tabs_s0
!        case = case0

        latitude0 = get_rlat_p(lchnk, icol)*57.296_r8
        longitude0 = get_rlon_p(lchnk, icol)*57.296_r8
!        pi = acos(-1.)
        if(fcor.eq.-999.) fcor= 4*pi/86400.*sin(latitude0*pi/180.)
        fcorz = sqrt(4.*(2*pi/(3600.*24.))**2-fcor**2)
        fcory(:) = fcor
        fcorzy(:) = fcorz
        do j=1,ny
          do i=1,nx
            latitude(i,j) = latitude0
            longitude(i,j) = longitude0
          end do
        end do

        if(ocnfrac.gt.0.5) then
           OCEAN = .true.
        else
           LAND = .true.
        end if

!        create CRM vertical grid and initialize some vertical reference arrays:
!
        do k = 1, nzm

           z(k) = zmid(plev-k+1) - zint(plev+1)
           zi(k) = zint(plev-k+2)- zint(plev+1)
           pres(k) = pmid(plev-k+1)/100.
           prespot(k)=(1000./pres(k))**(rgas/cp)
           bet(k) = ggr/tl(plev-k+1)
           gamaz(k)=ggr/cp*z(k)

        end do ! k
!        zi(nz) =  zint(plev-nz+2)
        zi(nz) = zint(plev-nz+2)-zint(plev+1) !+++mhwang, 2012-02-04

        dz = 0.5*(z(1)+z(2))
        do k=2,nzm
           adzw(k) = (z(k)-z(k-1))/dz
        end do
        adzw(1) = 1.
        adzw(nz) = adzw(nzm)
!        adz(1) = 1.
!        do k=2,nzm-1
!          adz(k) = 0.5*(z(k+1)-z(k-1))/dz
!        end do
!        adz(nzm) = adzw(nzm)
!+++mhwang fix the adz bug. (adz needs to be consistent with zi)
!2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
        do k=1, nzm
           adz(k)=(zi(k+1)-zi(k))/dz
        end do
        
        do k = 1,nzm
          rho(k) = pdel(plev-k+1)/ggr/(adz(k)*dz)
        end do
        do k=2,nzm
!          rhow(k) = 0.5*(rho(k)+rho(k-1))
!+++mhwang fix the rhow bug (rhow needes to be consistent with pmid)
!2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
           rhow(k) = (pmid(plev-k+2)-pmid(plev-k+1))/ggr/(adzw(k)*dz)
        end do
        rhow(1) = 2*rhow(2) - rhow(3)
#ifdef SPCAM_CLUBB_SGS /* Fix extropolation for 30 point grid */
        if (  2*rhow(nzm) - rhow(nzm-1) > 0. ) then
           rhow(nz)= 2*rhow(nzm) - rhow(nzm-1)
        else
           rhow(nz)= sqrt( rhow(nzm) )
        endif
#else
        rhow(nz)= 2*rhow(nzm) - rhow(nzm-1)
#endif /*CLUBB_SGS*/
        colprec=0
        colprecs=0

!  
!  Initialize:
!       
        

! limit the velocity at the very first step:
        
        if(u_crm(1,1,1).eq.u_crm(2,1,1).and.u_crm(3,1,2).eq.u_crm(4,1,2)) then
         do k=1,nzm
          do j=1,ny
           do i=1,nx
             u_crm(i,j,k) = min( umax, max(-umax,u_crm(i,j,k)) )
             v_crm(i,j,k) = min( umax, max(-umax,v_crm(i,j,k)) )*YES3D
           end do
          end do
         end do
        
        end if

        u(1:nx,1:ny,1:nzm) = u_crm(1:nx,1:ny,1:nzm)
        v(1:nx,1:ny,1:nzm) = v_crm(1:nx,1:ny,1:nzm)*YES3D
        w(1:nx,1:ny,1:nzm) = w_crm(1:nx,1:ny,1:nzm)
        tabs(1:nx,1:ny,1:nzm) = t_crm(1:nx,1:ny,1:nzm)
        micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields) = micro_fields_crm(1:nx,1:ny,1:nzm,1:nmicro_fields)
#ifdef sam1mom
        qn(1:nx,1:ny,1:nzm) =  micro_fields_crm(1:nx,1:ny,1:nzm,3)
#endif

#ifdef m2005
        cloudliq(1:nx,1:ny,1:nzm) = micro_fields_crm(1:nx,1:ny,1:nzm,11)
#endif

#ifdef m2005
        do k=1, nzm
#ifdef MODAL_AERO
! set aerosol data
         l=plev-k+1
         naer(k, 1:ntot_amode) = naermod(l, 1:ntot_amode)
         vaer(k, 1:ntot_amode) = vaerosol(l, 1:ntot_amode)
         hgaer(k, 1:ntot_amode) = hygro(l, 1:ntot_amode)
#endif
         do j=1, ny
          do i=1, nx
!            if(micro_field(i,j,k,iqcl).gt.0) then
            if(cloudliq(i,j,k).gt.0) then
              if(dopredictNc) then 
               if( micro_field(i,j,k,incl).eq.0) micro_field(i,j,k,incl) = 1.0e6*Nc0/rho(k)
              endif
            end if
          enddo
         enddo
        enddo
#endif 

        w(:,:,nz)=0.
        wsub (:) = 0.      !used in clubb, +++mhwang
        dudt(:,:,:,1:3) = 0.
        dvdt(:,:,:,1:3) = 0.
        dwdt(1:nx,1:ny,1:nz,1:3) = 0.
        tke(1:nx,1:ny,1:nzm) = 0.
        tk(1:nx,1:ny,1:nzm) = 0.
        tkh(1:nx,1:ny,1:nzm) = 0.
        p(1:nx,1:ny,1:nzm) = 0.

        CF3D(1:nx,1:ny,1:nzm) = 1.

        call micro_init

! initialize sgs fields
        call sgs_init
        
        do k=1,nzm
          
          u0(k)=0.
          v0(k)=0.
          t0(k)=0.
          t00(k)=0.
          tabs0(k)=0.
          q0(k)=0.
          qv0(k)=0.
!+++mhwang these are not initialized ??
          qn0(k) = 0.0
          qp0(k) = 0.0
          tke0(k) = 0.0
!---mhwang
          do j=1,ny
           do i=1,nx
            
            t(i,j,k) = tabs(i,j,k)+gamaz(k) &
                        -fac_cond*qcl(i,j,k)-fac_sub*qci(i,j,k) &
                        -fac_cond*qpl(i,j,k)-fac_sub*qpi(i,j,k)

            colprec=colprec+(qpl(i,j,k)+qpi(i,j,k))*pdel(plev-k+1)
            colprecs=colprecs+qpi(i,j,k)*pdel(plev-k+1)
            u0(k)=u0(k)+u(i,j,k)
            v0(k)=v0(k)+v(i,j,k)
            t0(k)=t0(k)+t(i,j,k)
            t00(k)=t00(k)+t(i,j,k)+fac_cond*qpl(i,j,k)+fac_sub*qpi(i,j,k)
            tabs0(k)=tabs0(k)+tabs(i,j,k)
            q0(k)=q0(k)+qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)
            qv0(k) = qv0(k) + qv(i,j,k)
            qn0(k) = qn0(k) + qcl(i,j,k) + qci(i,j,k)
            qp0(k) = qp0(k) + qpl(i,j,k) + qpi(i,j,k)
            tke0(k)=tke0(k)+tke(i,j,k)

           end do
          end do

          u0(k) = u0(k) * factor_xy
          v0(k) = v0(k) * factor_xy
          t0(k) = t0(k) * factor_xy
          t00(k) = t00(k) * factor_xy
          tabs0(k) = tabs0(k) * factor_xy
          q0(k) = q0(k) * factor_xy
          qv0(k) = qv0(k) * factor_xy
          qn0(k) = qn0(k) * factor_xy
          qp0(k) = qp0(k) * factor_xy
          tke0(k) = tke0(k) * factor_xy

#ifdef SPCAM_CLUBB_SGS
 ! Update thetav for CLUBB.  This is needed when we have a higher model top 
 ! than is in the sounding, because we subsequently use tv0 to initialize 
 ! thv_ds_zt/zm, which appear in CLUBB's anelastic buoyancy terms. 
 ! -dschanen UWM 11 Feb 2010
          tv0(k) = tabs0(k)*prespot(k)*(1.+epsv*q0(k))
#endif

          l = plev-k+1
          uln(l) = min( umax, max(-umax,ul(l)) )
          vln(l) = min( umax, max(-umax,vl(l)) )*YES3D
          ttend(k) = (tl(l)+gamaz(k)- &
               fac_cond*(qccl(l)+qiil(l))-fac_fus*qiil(l)-t00(k))*idt_gl
          qtend(k) = (ql(l)+qccl(l)+qiil(l)-q0(k))*idt_gl
          utend(k) = (uln(l)-u0(k))*idt_gl
          vtend(k) = (vln(l)-v0(k))*idt_gl
          ug0(k) = uln(l)
          vg0(k) = vln(l)
          tg0(k) = tl(l)+gamaz(k)-fac_cond*qccl(l)-fac_sub*qiil(l)
          qg0(k) = ql(l)+qccl(l)+qiil(l)

        end do ! k

        uhl = u0(1)
        vhl = v0(1)

! estimate roughness length assuming logarithmic profile of velocity near the surface:

        ustar = sqrt(tau00/rho(1))
        z0 = z0_est(z(1),bflx,wnd,ustar)
        z0 = max(0.00001,min(1.,z0))

        timing_factor = 0.

        prectend=colprec
        precstend=colprecs

#ifdef SPCAM_CLUBB_SGS
        if(doclubb) then
          fluxbu(:, :) = fluxu00/rhow(1)
          fluxbv(:, :) = fluxv00/rhow(1)
          fluxbt(:, :) = fluxt00/rhow(1)
          fluxbq(:, :) = fluxq00/rhow(1)
        else
          fluxbu(:, :) = 0.
          fluxbv(:, :) = 0.
          fluxbt(:, :) = 0.
          fluxbq(:, :) = 0.
        end if
#else 
        fluxbu=0.
        fluxbv=0.
        fluxbt=0.
        fluxbq=0.
#endif /*CLUBB_SGS*/

        fluxtu=0.
        fluxtv=0.
        fluxtt=0.
        fluxtq=0.
        fzero =0.
        precsfc=0.
        precssfc=0.

!---------------------------------------------------
        cld = 0.
        cldtop = 0.
        gicewp=0
        gliqwp=0
        mc = 0.
        mcup = 0.
        mcdn = 0.
        mcuup = 0.
        mcudn = 0.
        crm_qc = 0.
        crm_qi = 0.
        crm_qs = 0.
        crm_qg = 0.
        crm_qr = 0.
#ifdef m2005
        crm_nc = 0.
        crm_ni = 0.
        crm_ns = 0.
        crm_ng = 0.
        crm_nr = 0.
! hm 8/31/11 add new variables
        aut_crm_a = 0.
        acc_crm_a = 0.
        evpc_crm_a = 0.
        evpr_crm_a = 0.
        mlt_crm_a = 0.
        sub_crm_a = 0.
        dep_crm_a = 0.
        con_crm_a = 0.

! hm 8/31/11 add new output
! these are increments added to calculate gcm-grid and time-step avg
! note - these values are also averaged over the icycle loop following
! the approach for precsfc
        aut1a = 0.
        acc1a = 0.
        evpc1a = 0.
        evpr1a = 0.
        mlt1a = 0.
        sub1a = 0.
        dep1a = 0.
        con1a = 0.

#endif 

        mu_crm = 0.
        md_crm = 0.
        eu_crm = 0.
        du_crm = 0.
        ed_crm = 0.
        dd_crm = 0.
        jt_crm = 0.
        mx_crm = 0.

        mui_crm = 0.
        mdi_crm = 0.

        flux_qt = 0.
        flux_u = 0.
        flux_v = 0.
        fluxsgs_qt = 0.
        tkez = 0.
        tkesgsz = 0.
        tkz = 0.
        flux_qp = 0.
        pflx = 0.
        qt_trans = 0.
        qp_trans = 0.
        qp_fall = 0.
        qp_evp = 0.
        qp_src = 0.
        qt_ls = 0.
        t_ls = 0.

        uwle = 0.
        uwsb = 0.
        vwle = 0.
        vwsb = 0.
        qpsrc = 0.
        qpevp = 0.
        qpfall = 0.
        precflux = 0.

       prec_xy = 0.0
       total_water_evap = 0.0
       total_water_prec = 0.0
       tlat = 0.0
       pw_xy = 0.0; cw_xy=0.0; iw_xy = 0.0
       usfc_xy = 0.0; vsfc_xy =0.0; u200_xy =0.0; v200_xy = 0.0; w500_xy = 0.0
       swvp_xy = 0.0; psfc_xy = 0.0; u850_xy = 0.0; v850_xy = 0.0

!--------------------------------------------------
#ifdef sam1mom
     if(doprecip) call precip_init()
#endif

        call get_gcol_all_p(lchnk, pcols, gcolindex)
        iseed = gcolindex(icol)
        if(u(1,1,1).eq.u(2,1,1).and.u(3,1,2).eq.u(4,1,2)) &
                    call setperturb(iseed)

#ifndef SPCAM_CLUBB_SGS
!--------------------------
! do a CLUBB sanity check
        if ( doclubb .or. doclubbnoninter ) then
           write(0,*) "Cannot call CLUBB if -DCLUBB is not in FFLAGS"
           call endrun('crm main')
        end if
#endif /*CLUBB_SGS*/
#ifdef SPCAM_CLUBB_SGS
!------------------------------------------------------------------
! Do initialization for UWM CLUBB
!------------------------------------------------------------------
        up2(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 1)
        vp2(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 2)
        wprtp(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 3)
        wpthlp(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 4)
        wp2(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 5)
        wp3(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 6)
        rtp2(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 7)
        thlp2(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 8)
        rtpthlp(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 9)
        upwp(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 10)
        vpwp(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 11)
        cloud_frac(1:nx, 1:ny, 1:nz) = clubb_buffer(1:nx, 1:ny, 1:nz, 12)
        t_tndcy(1:nx, 1:ny, 1:nzm) = clubb_buffer(1:nx, 1:ny, 1:nzm, 13)
        qc_tndcy(1:nx, 1:ny, 1:nzm) = clubb_buffer(1:nx, 1:ny, 1:nzm, 14)
        qv_tndcy(1:nx, 1:ny, 1:nzm) = clubb_buffer(1:nx, 1:ny, 1:nzm, 15)
        u_tndcy(1:nx, 1:ny, 1:nzm) = clubb_buffer(1:nx, 1:ny, 1:nzm, 16)
        v_tndcy(1:nx, 1:ny, 1:nzm) = clubb_buffer(1:nx, 1:ny, 1:nzm, 17)

!
! since no tracer is carried in the current version of MMF, these 
! tracer-related restart varialbes are set to zero. +++mhwang, 2011-08
        tracer_tndcy = 0.0
        sclrp2 = 0.0
        sclrprtp = 0.0
        sclrpthlp = 0.0
        wpsclrp =0.0

        if((doclubb.and.docloud).or.(.not.doclubb .and. .not.docloud)) then
          write(0, *) 'doclubb and docloud can not both be true or be false'
          call endrun('crm_clubb2') 
        end if
        if((doclubb_sfc_fluxes.and.docam_sfc_fluxes)) then
          write(0, *) 'doclubb_sfc_fluxes and dosam_sfc_fluxes can not both be true'
          call endrun('crm_clubb_fluxes')
        end if

        if ( doclubb .or. doclubbnoninter ) then
	  call clubb_sgs_setup( real( dt*real( nclubb ), kind=time_precision), &
	    latitude, longitude, z, rho, zi, rhow, tv0, tke )
        end if
#endif  /*CLUBB_SGS*/

#ifdef ECPP
!        ntavg1_ss = dt_gl/3   ! one third of GCM time step, 10 minutes
        ntavg1_ss = min(600._r8, dt_gl)       ! 10 minutes  or the GCM timestep, whichever smaller
              ! ntavg1_ss = number of seconds to average between computing categories.
        ntavg2_ss = dt_gl   ! GCM time step
              ! ntavg2_ss = number of seconds to average between outputs.
              !    This must be a multiple of ntavgt1_ss.
!
!  ecpp_crm_init has to be called after ntavg1_ss and ntavg2_ss are set for
!  their values are used in ecpp_crm_init. 
        call ecpp_crm_init()

        qlsink = 0.0
        qlsink_bf = 0.0
        prain = 0.0
        precr = 0.0
        precsolid = 0.0
#endif /*ECPP*/

!+++mhwangtest
! test water conservtion problem
        ntotal_step = 0.0
        qtot(:) = 0.0
        qtotmicro(:) = 0.0
        do k=1, nzm
         l=plev-k+1
         do j=1, ny
          do i=1, nx
#ifdef m2005
            qtot(1) = qtot(1)+((micro_field(i,j,k,iqr)+micro_field(i,j,k,iqs)+micro_field(i,j,k,iqg)) * pdel(l)/ggr)/(nx*ny) 
#endif
#ifdef sam1mom
            qtot(1) = qtot(1)+(qpl(i,j,k)+qpi(i,j,k)) * pdel(l)/ggr/(nx*ny)
#endif
          enddo
         enddo
         qtot(1) = qtot(1) + (ql(l)+qccl(l)+qiil(l)) * pdel(l)/ggr
        enddo
!---mhwangtest

        nstop = dt_gl/dt
        dt = dt_gl/nstop
        nsave3D = nint(60/dt)
!       if(nint(nsave3D*dt).ne.60)then
!          print *,'CRM: time step=',dt,' is not divisible by 60 seconds'
!          print *,'this is needed for output every 60 seconds'
!          stop
!       endif
        nstep = 0
        nprint = 1
        ncycle = 0
!        nrad = nstop/nrad0
        day=day0

!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------

do while(nstep.lt.nstop) 
        
  nstep = nstep + 1
  time = time + dt
  day = day0 + time/86400.
  timing_factor = timing_factor+1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()

  do icyc=1,ncycle

     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

!-----------------------------------------------------------
!       Buoyancy term:
	     
     call buoyancy()

!+++mhwangtest
! test water conservtion problem
        ntotal_step = ntotal_step + 1.
!---mhwangtest 

!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()

     do k=1,nzm
      do j=1,ny
        do i=1,nx
          t(i,j,k) = t(i,j,k) + qrad_crm(i,j,k)*dtn
        end do
      end do
     end do

!----------------------------------------------------------
!   	suppress turbulence near the upper boundary (spange):

     if(dodamping) call damping()

!---------------------------------------------------------
!   Ice fall-out

#ifdef SPCAM_CLUBB_SGS
      if ( docloud .or. doclubb ) then
        call ice_fall()
      end if
#else
      if(docloud) then
          call ice_fall()
      end if
#endif  /*CLUBB_SGS*/ 

!----------------------------------------------------------
!     Update scalar boundaries after large-scale processes:

     call boundaries(3)

!---------------------------------------------------------
!     Update boundaries for velocities:

      call boundaries(0)

!-----------------------------------------------
!     surface fluxes:

     if(dosurface) call crmsurface(bflx)

!-----------------------------------------------------------
!  SGS physics:

     if (dosgs) call sgs_proc()
        
#ifdef CLUBB_CRM_OLD   
!----------------------------------------------------------
!     Do a timestep with CLUBB if enabled:
!     -dschanen UWM 16 May 2008

      if ( doclubb .or. doclubbnoninter ) then
        ! In case of ice fall, we recompute qci here for the 
        ! single-moment scheme.  Also, subsidence, diffusion and advection have
        ! been applied to micro_field but not qv/qcl so they must be updated.
        call micro_update()
      end if ! doclubb .or. doclubbnoninter

      if ( doclubb ) then
        ! Calculate the vertical integrals for RTM and THLM so we can later
        ! calculate whether CLUBB is a spurious source or sink of either.
        ! - nielsenb UWM 4 Jun 2010
        do i = 1,nx
          do j = 1,ny
            rtm_column = qv(i,j,1:nzm) + qcl(i,j,1:nzm)
            rtm_integral_before(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), & 
                                         rtm_column, gr%invrs_dzt(2:nz) )

            thlm_before = t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                                 qcl(i,j,1:nzm), qpl(i,j,1:nzm), &
                                 qci(i,j,1:nzm), qpi(i,j,1:nzm), &
                                 prespot(1:nzm) )

            thlm_integral_before(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                                          thlm_before(1:nzm), gr%invrs_dzt(2:nz) )
          end do
        end do
        ! End vertical integral

      end if ! doclubb

      if ( doclubb .or. doclubbnoninter ) then

        ! We call CLUBB here because adjustments to the wind
        ! must occur prior to adams() -dschanen 26 Aug 2008
        ! Here we call clubb only if nstep divides the current timestep,
        ! or we're on the very first timestep
        if ( nstep == 1 .or. mod( nstep, nclubb ) == 0 ) then

          call advance_clubb_sgs &
               ( real( dtn*real( nclubb ), kind=time_precision), & ! in
                 real( 0., kind=time_precision ),         & ! in
                 real( time, kind=time_precision ),       & ! in
                 rho, rhow, wsub, u, v, w, qpl, qci, qpi, & ! in
                 t, qv, qcl ) ! in
        end if ! nstep == 1 .or. mod( nstep, nclubb) == 0

      end if ! doclubb .or. doclubbnoninter

#endif  /*CLUBB_CRM_OLD*/
!----------------------------------------------------------
!     Fill boundaries for SGS diagnostic fields:

     call boundaries(4)
!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!----------------------------------------------------------
!	SGS effects on momentum:

     if(dosgs) call sgs_mom()
#ifdef CLUBB_CRM_OLD
     if ( doclubb ) then
!          call apply_clubb_sgs_tndcy_mom &
!               ( dudt, dvdt ) ! in/out
     endif
#endif /*CLUBB_CRM_OLD*/

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
!  Note that at the end of the call, the velocities are in nondimensional form.
	 
     call adams()

!----------------------------------------------------------
!     Update boundaries for all prognostic scalar fields for advection:

     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :

     call advect_all_scalars()

!-----------------------------------------------------------
!    Convert velocity back from nondimensional form:

      call uvw()

!----------------------------------------------------------
!     Update boundaries for scalars to prepare for SGS effects:

     call boundaries(3)

!---------------------------------------------------------
!      SGS effects on scalars :

     if (dosgs) call sgs_scalars()

#ifdef CLUBB_CRM_OLD
      ! Re-compute q/qv/qcl based on values computed in CLUBB
     if ( doclubb ) then

      ! Recalculate q, qv, qcl based on new micro_fields (updated by horizontal
      ! diffusion)
       call micro_update()

      ! Then Re-compute q/qv/qcl based on values computed in CLUBB
       call apply_clubb_sgs_tndcy_scalars &
            ( real( dtn, kind=time_precision), & ! in
              t, qv, qcl) ! in/out

       call micro_adjust( qv, qcl ) ! in

       ! Calculate the vertical integrals for RTM and THLM again so
       ! calculate whether CLUBB is a spurious source or sink of either.
       ! - nielsenb UWM 4 Jun 2010
       do i = 1,nx
         do j = 1,ny
           rtm_flux_top = rho_ds_zm(nz) * wprtp(i,j,nz)
           rtm_flux_sfc = rho_ds_zm(1) * fluxbq(i,j)
           rtm_column = qv(i,j,1:nzm) + qcl(i,j,1:nzm)
           rtm_integral_after(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), & 
                                         rtm_column, gr%invrs_dzt(2:nz) )
                                          
           rtm_spurious_source(i,j) = calculate_spurious_source( rtm_integral_after(i,j), &
                                                      rtm_integral_before(i,j), &
                                                      rtm_flux_top, rtm_flux_sfc, &
                                                      0.0_core_rknd, real( dtn, kind=core_rknd) )

           thlm_flux_top = rho_ds_zm(nz) * wpthlp(i,j,nz)
           thlm_flux_sfc = rho_ds_zm(1) * fluxbt(i,j)

           thlm_after = t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                                  qcl(i,j,1:nzm), qpl(i,j,1:nzm), &
                                  qci(i,j,1:nzm), qpi(i,j,1:nzm), &
                                  prespot(1:nzm) )

           thlm_integral_after(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                                      thlm_after(1:nzm), gr%invrs_dzt(2:nz))
                                         
           thlm_spurious_source(i,j) = calculate_spurious_source( thlm_integral_after(i,j), &
                                                          thlm_integral_before(i,j), &
                                                          thlm_flux_top, thlm_flux_sfc, &
                                                          0.0_core_rknd, real( dtn, kind=core_rknd ))
         end do
       end do
          ! End spurious source calculation

     end if! doclubb
#endif /*CLUBB_CRM_OLD*/

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:
#ifdef SPCAM_CLUBB_SGS
      if(docloud.or.dosmoke.or.doclubb) call micro_proc()
#else
      if(docloud.or.dosmoke) call micro_proc()
#endif /*CLUBB_SGS*/

!-----------------------------------------------------------
!    Compute diagnostics fields:

      call diagnose()

!----------------------------------------------------------
! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

   end do ! icycle	
          
!----------------------------------------------------------
!----------------------------------------------------------
#ifdef ECPP
! Here ecpp_crm_stat is called every CRM time step (dt), not every subcycle time step (dtn). 
! This is what the original MMF model did (t_rad, qv_rad, ...). Do we want to call ecpp_crm_stat
! every subcycle time step??? +++mhwang
       call ecpp_crm_stat()   
#endif /*ECPP*/

        cwp = 0.
        cwph = 0.
        cwpm = 0.
        cwpl = 0.

        flag_top(:,:) = .true.

        cltemp = 0.0; cmtemp = 0.0
        chtemp = 0.0; cttemp = 0.0

        do k=1,nzm
         l = plev-k+1
         do j=1,ny
          do i=1,nx

! hm modify 9/7/11 for end of timestep, GCM-grid scale hydrometeor output
! instead of time-step-averaged
! I also modified this for all q and N variables as well as for sam1mom
! for consistency
!hm           crm_qc(l) = crm_qc(l) + qcl(i,j,k)
!hm           crm_qi(l) = crm_qi(l) + qci(i,j,k)
!hm           crm_qr(l) = crm_qr(l) + qpl(i,j,k)
!hm#ifdef sam1mom
!hm           omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
!hm           crm_qg(l) = crm_qg(l) + qpi(i,j,k)*omg
!hm           crm_qs(l) = crm_qs(l) + qpi(i,j,k)*(1.-omg)
!hm#else
!           crm_qg(l) = crm_qg(l) + qpi(i,j,k)
!           crm_qs(l) = crm_qs(l) + 0.     ! temporerary solution
!hm           crm_qg(l) = crm_qg(l) + micro_field(i,j,k,iqg)
!hm           crm_qs(l) = crm_qs(l) + micro_field(i,j,k,iqs)   

!hm           crm_nc(l) = crm_nc(l) + micro_field(i,j,k,incl)
!hm           crm_ni(l) = crm_ni(l) + micro_field(i,j,k,inci)
!hm           crm_nr(l) = crm_nr(l) + micro_field(i,j,k,inr)
!hm           crm_ng(l) = crm_ng(l) + micro_field(i,j,k,ing)
!hm           crm_ns(l) = crm_ns(l) + micro_field(i,j,k,ins)

!hm#endif

           tmp1 = rho(nz-k)*adz(nz-k)*dz*(qcl(i,j,nz-k)+qci(i,j,nz-k))
           cwp(i,j) = cwp(i,j)+tmp1
           cttemp(i,j) = max(CF3D(i,j,nz-k), cttemp(i,j))
           if(cwp(i,j).gt.cwp_threshold.and.flag_top(i,j)) then
               cldtop(k) = cldtop(k) + 1
               flag_top(i,j) = .false.
           end if
           if(pres(nz-k).ge.700.) then
               cwpl(i,j) = cwpl(i,j)+tmp1
               cltemp(i,j) = max(CF3D(i,j,nz-k), cltemp(i,j))
           else if(pres(nz-k).lt.400.) then
               cwph(i,j) = cwph(i,j)+tmp1
               chtemp(i,j) = max(CF3D(i,j,nz-k), chtemp(i,j))
           else
               cwpm(i,j) = cwpm(i,j)+tmp1
               cmtemp(i,j) = max(CF3D(i,j,nz-k), cmtemp(i,j))
           end if

      !     qsat = qsatw_crm(tabs(i,j,k),pres(k))
      !     if(qcl(i,j,k)+qci(i,j,k).gt.min(1.e-5,0.01*qsat)) then
           tmp1 = rho(k)*adz(k)*dz
           if(tmp1*(qcl(i,j,k)+qci(i,j,k)).gt.cwp_threshold) then
                cld(l) = cld(l) + CF3D(i,j,k)
                if(w(i,j,k+1)+w(i,j,k).gt.2*wmin) then
                  mcup(l) = mcup(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * CF3D(i,j,k)
                  mcuup(l) = mcuup(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * (1.0 - CF3D(i,j,k))
                end if
                if(w(i,j,k+1)+w(i,j,k).lt.-2*wmin) then
                  mcdn(l) = mcdn(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * CF3D(i,j,k)
                  mcudn(l) = mcudn(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * (1. - CF3D(i,j,k))
                end if
           else 
                if(w(i,j,k+1)+w(i,j,k).gt.2*wmin) then
                  mcuup(l) = mcuup(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                end if
                if(w(i,j,k+1)+w(i,j,k).lt.-2*wmin) then
                  mcudn(l) = mcudn(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                end if
           end if
           
           t_rad (i,j,k) = t_rad (i,j,k)+tabs(i,j,k)
           qv_rad(i,j,k) = qv_rad(i,j,k)+max(0.,qv(i,j,k))
           qc_rad(i,j,k) = qc_rad(i,j,k)+qcl(i,j,k)
           qi_rad(i,j,k) = qi_rad(i,j,k)+qci(i,j,k)
           cld_rad(i,j,k) = cld_rad(i,j,k) +  CF3D(i,j,k)
#ifdef m2005
           nc_rad(i,j,k) = nc_rad(i,j,k)+micro_field(i,j,k,incl)
           ni_rad(i,j,k) = ni_rad(i,j,k)+micro_field(i,j,k,inci)
           qs_rad(i,j,k) = qs_rad(i,j,k)+micro_field(i,j,k,iqs)
           ns_rad(i,j,k) = ns_rad(i,j,k)+micro_field(i,j,k,ins)
#endif 
           gliqwp(l)=gliqwp(l)+qcl(i,j,k)
           gicewp(l)=gicewp(l)+qci(i,j,k)
          
          end do
         end do
        end do

! Diagnose mass fluxes to drive CAM's convective transport of tracers. 
! definition of mass fluxes is taken from Xu et al., 2002, QJRMS. 
        do k=1, nzm+1
         l=plev+1-k+1
         do j=1, ny
          do i=1, nx
           if(w(i,j,k).gt.0.) then
             kx=max(1, k-1)
             qsat = qsatw_crm(tabs(i,j,kx),pres(kx))
             if(qcl(i,j,kx)+qci(i,j,kx).gt.min(1.e-5,0.01*qsat)) then
              mui_crm(l) = mui_crm(l)+rhow(k)*w(i,j,k)
             end if
           else if (w(i,j,k).lt.0.) then 
             kx=min(k+1, nzm)
             qsat = qsatw_crm(tabs(i,j,kx),pres(kx))
             if(qcl(i,j,kx)+qci(i,j,kx).gt.min(1.e-5,0.01*qsat)) then
               mdi_crm(l) = mdi_crm(l)+rhow(k)*w(i,j,k)
             else if(qpl(i,j,kx)+qpi(i,j,kx).gt.1.0e-4) then 
               mdi_crm(l) = mdi_crm(l)+rhow(k)*w(i,j,k) 
             end if
           end if
          end do
         end do
        end do

!        do k=1,nzm
!         radlwup0(k)=radlwup0(k)+radlwup(k)
!         radlwdn0(k)=radlwdn0(k)+radlwdn(k)
!         radqrlw0(k)=radqrlw0(k)+radqrlw(k)
!         radswup0(k)=radswup0(k)+radswup(k)
!         radswdn0(k)=radswdn0(k)+radswdn(k)
!         radqrsw0(k)=radqrsw0(k)+radqrsw(k)
!        end do
        
        do j=1,ny
         do i=1,nx
!           if(cwp(i,j).gt.cwp_threshold) cltot = cltot + 1.
!           if(cwph(i,j).gt.cwp_threshold) clhgh = clhgh + 1.
!           if(cwpm(i,j).gt.cwp_threshold) clmed = clmed + 1.
!           if(cwpl(i,j).gt.cwp_threshold) cllow = cllow + 1.
!  use maxmimum cloud overlap to calcluate cltot, clhgh, 
!  cldmed, and cldlow   +++ mhwang
           if(cwp(i,j).gt.cwp_threshold) cltot = cltot + cttemp(i,j) 
           if(cwph(i,j).gt.cwp_threshold) clhgh = clhgh + chtemp(i,j) 
           if(cwpm(i,j).gt.cwp_threshold) clmed = clmed + cmtemp(i,j) 
           if(cwpl(i,j).gt.cwp_threshold) cllow = cllow + cltemp(i,j) 
         end do
        end do

!        call stepout()
!----------------------------------------------------------
        end do ! main loop
!----------------------------------------------------------

        tmp1 = 1._r8/ dble(nstop)
        t_rad = t_rad * tmp1
        qv_rad = qv_rad * tmp1
        qc_rad = qc_rad * tmp1
        qi_rad = qi_rad * tmp1
        cld_rad = cld_rad * tmp1
#ifdef m2005
        nc_rad = nc_rad * tmp1
        ni_rad = ni_rad * tmp1
        qs_rad = qs_rad * tmp1
        ns_rad = ns_rad * tmp1
#endif

! no CRM tendencies above its top
        
        tln(1:ptop-1) = tl(1:ptop-1)
        qln(1:ptop-1) = ql(1:ptop-1)
        qccln(1:ptop-1)= qccl(1:ptop-1)
        qiiln(1:ptop-1)= qiil(1:ptop-1)
        uln(1:ptop-1) = ul(1:ptop-1)
        vln(1:ptop-1) = vl(1:ptop-1)

!  Compute tendencies due to CRM:
        
        tln(ptop:plev) = 0.
        qln(ptop:plev) = 0.
        qccln(ptop:plev)= 0.
        qiiln(ptop:plev)= 0.
        uln(ptop:plev) = 0.
        vln(ptop:plev) = 0.
        
        colprec=0
        colprecs=0
        do k = 1,nzm
         l = plev-k+1
         do i=1,nx
          do j=1,ny
             colprec=colprec+(qpl(i,j,k)+qpi(i,j,k))*pdel(plev-k+1)
             colprecs=colprecs+qpi(i,j,k)*pdel(plev-k+1)
             tln(l) = tln(l)+tabs(i,j,k)
             qln(l) = qln(l)+qv(i,j,k)
             qccln(l)= qccln(l)+qcl(i,j,k)
             qiiln(l)= qiiln(l)+qci(i,j,k)
             uln(l) = uln(l)+u(i,j,k)
             vln(l) = vln(l)+v(i,j,k)
          end do ! k
         end do
        end do ! i


        tln(ptop:plev) = tln(ptop:plev) * factor_xy
        qln(ptop:plev) = qln(ptop:plev) * factor_xy
        qccln(ptop:plev) = qccln(ptop:plev) * factor_xy
        qiiln(ptop:plev) = qiiln(ptop:plev) * factor_xy
        uln(ptop:plev) = uln(ptop:plev) * factor_xy
        vln(ptop:plev) = vln(ptop:plev) * factor_xy

        sltend = cp * (tln - tl) * idt_gl
        qltend = (qln - ql) * idt_gl
        qcltend = (qccln - qccl) * idt_gl
        qiltend = (qiiln - qiil) * idt_gl
        prectend=(colprec-prectend)/ggr*factor_xy * idt_gl
        precstend=(colprecs-precstend)/ggr*factor_xy * idt_gl

! don't use CRM tendencies from two crm top levels
        sltend(ptop:ptop+1) = 0.
        qltend(ptop:ptop+1) = 0.
        qcltend(ptop:ptop+1) = 0.
        qiltend(ptop:ptop+1) = 0.
!-------------------------------------------------------------
! 
! Save the last step to the permanent core:
        
        u_crm  (1:nx,1:ny,1:nzm) = u   (1:nx,1:ny,1:nzm)
        v_crm  (1:nx,1:ny,1:nzm) = v   (1:nx,1:ny,1:nzm)
        w_crm  (1:nx,1:ny,1:nzm) = w   (1:nx,1:ny,1:nzm)
        t_crm  (1:nx,1:ny,1:nzm) = tabs(1:nx,1:ny,1:nzm)
        micro_fields_crm(1:nx,1:ny,1:nzm,1:nmicro_fields) = micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields)
#ifdef sam1mom
        micro_fields_crm(1:nx,1:ny,1:nzm,3) = qn(1:nx,1:ny,1:nzm)
#endif
#ifdef m2005
        micro_fields_crm(1:nx,1:ny,1:nzm,11) = cloudliq(1:nx,1:ny,1:nzm)
#endif
       crm_tk(1:nx,1:ny,1:nzm) = tk(1:nx, 1:ny, 1:nzm)
       crm_tkh(1:nx,1:ny,1:nzm) = tkh(1:nx, 1:ny, 1:nzm)
       cld3d_crm(1:nx, 1:ny, 1:nzm) = CF3D(1:nx, 1:ny, 1:nzm)
#ifdef SPCAM_CLUBB_SGS
       clubb_buffer(1:nx, 1:ny, 1:nz, 1) = up2(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 2) = vp2(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 3) = wprtp(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 4) = wpthlp(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 5) = wp2(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 6) = wp3(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 7) = rtp2(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 8) = thlp2(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 9) = rtpthlp(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 10) = upwp(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 11) = vpwp(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nz, 12) = cloud_frac(1:nx, 1:ny, 1:nz)
       clubb_buffer(1:nx, 1:ny, 1:nzm, 13) = t_tndcy(1:nx, 1:ny, 1:nzm)
       clubb_buffer(1:nx, 1:ny, 1:nzm, 14) = qc_tndcy(1:nx, 1:ny, 1:nzm)
       clubb_buffer(1:nx, 1:ny, 1:nzm, 15) = qv_tndcy(1:nx, 1:ny, 1:nzm)
       clubb_buffer(1:nx, 1:ny, 1:nzm, 16) = u_tndcy(1:nx, 1:ny, 1:nzm)
       clubb_buffer(1:nx, 1:ny, 1:nzm, 17) = v_tndcy(1:nx, 1:ny, 1:nzm)

       crm_cld(1:nx, 1:ny, 1:nz) = cloud_frac(1:nx, 1:ny, 1:nz)
       clubb_tk(1:nx,1:ny,1:nzm) = tk_clubb(1:nx, 1:ny, 1:nzm)
       clubb_tkh(1:nx,1:ny,1:nzm) = tkh_clubb(1:nx, 1:ny, 1:nzm)
       relvar(1:nx, 1:ny, 1:nzm) = relvarg(1:nx, 1:ny, 1:nzm)
       accre_enhan(1:nx, 1:ny, 1:nzm) = accre_enhang(1:nx, 1:ny, 1:nzm) 
       qclvar(1:nx, 1:ny, 1:nzm) = qclvarg(1:nx, 1:ny, 1:nzm)
#endif

        do k=1,nzm
         do j=1,ny
          do i=1,nx
            qc_crm(i,j,k) = qcl(i,j,k)
            qi_crm(i,j,k) = qci(i,j,k)
            qpc_crm(i,j,k) = qpl(i,j,k)
            qpi_crm(i,j,k) = qpi(i,j,k)
#ifdef m2005
            wvar_crm(i,j,k) = wvar(i,j,k)
! hm 7/26/11, new output
            aut_crm(i,j,k) = aut1(i,j,k)
            acc_crm(i,j,k) = acc1(i,j,k)
            evpc_crm(i,j,k) = evpc1(i,j,k)
            evpr_crm(i,j,k) = evpr1(i,j,k)
            mlt_crm(i,j,k) = mlt1(i,j,k)
            sub_crm(i,j,k) = sub1(i,j,k)
            dep_crm(i,j,k) = dep1(i,j,k)
            con_crm(i,j,k) = con1(i,j,k)
#endif
          end do
         end do
        end do
        z0m = z0 
        taux_crm = taux0 / dble(nstop)
        tauy_crm = tauy0 / dble(nstop)

!---------------------------------------------------------------
!
!  Diagnostics:

! hm add 9/7/11, change from GCM-time step avg to end-of-timestep

        do k=1,nzm
         l = plev-k+1
         do j=1,ny
          do i=1,nx

           crm_qc(l) = crm_qc(l) + qcl(i,j,k)
           crm_qi(l) = crm_qi(l) + qci(i,j,k)
           crm_qr(l) = crm_qr(l) + qpl(i,j,k)
#ifdef sam1mom
           omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
           crm_qg(l) = crm_qg(l) + qpi(i,j,k)*omg
           crm_qs(l) = crm_qs(l) + qpi(i,j,k)*(1.-omg)
#else
!           crm_qg(l) = crm_qg(l) + qpi(i,j,k)
!           crm_qs(l) = crm_qs(l) + 0.     ! temporerary solution
           crm_qg(l) = crm_qg(l) + micro_field(i,j,k,iqg)
           crm_qs(l) = crm_qs(l) + micro_field(i,j,k,iqs)

           crm_nc(l) = crm_nc(l) + micro_field(i,j,k,incl)
           crm_ni(l) = crm_ni(l) + micro_field(i,j,k,inci)
           crm_nr(l) = crm_nr(l) + micro_field(i,j,k,inr)
           crm_ng(l) = crm_ng(l) + micro_field(i,j,k,ing)
           crm_ns(l) = crm_ns(l) + micro_field(i,j,k,ins)
#endif

	  end do
	 end do
	end do

        cld = min(1._r8,cld/float(nstop)*factor_xy)
        cldtop = min(1._r8,cldtop/float(nstop)*factor_xy)
        gicewp(:)=gicewp*pdel(:)*1000./ggr/float(nstop)*factor_xy
        gliqwp(:)=gliqwp*pdel(:)*1000./ggr/float(nstop)*factor_xy
        mcup = mcup / float(nstop) * factor_xy
        mcdn = mcdn / float(nstop) * factor_xy
        mcuup = mcuup / float(nstop) * factor_xy
        mcudn = mcudn / float(nstop) * factor_xy
        mc = mcup + mcdn + mcuup + mcudn
! hm 9/7/11 modify for end-of-timestep instead of timestep-avg output
!hm        crm_qc = crm_qc / float(nstop) * factor_xy
!hm        crm_qi = crm_qi / float(nstop) * factor_xy
!hm        crm_qs = crm_qs / float(nstop) * factor_xy
!hm        crm_qg = crm_qg / float(nstop) * factor_xy
!hm        crm_qr = crm_qr / float(nstop) * factor_xy
!hm#ifdef m2005
!hm        crm_nc = crm_nc / float(nstop) * factor_xy
!hm        crm_ni = crm_ni / float(nstop) * factor_xy
!hm        crm_ns = crm_ns / float(nstop) * factor_xy
!hm        crm_ng = crm_ng / float(nstop) * factor_xy
!hm        crm_nr = crm_nr / float(nstop) * factor_xy

        crm_qc = crm_qc * factor_xy
        crm_qi = crm_qi * factor_xy
        crm_qs = crm_qs * factor_xy
        crm_qg = crm_qg * factor_xy
        crm_qr = crm_qr * factor_xy
#ifdef m2005
        crm_nc = crm_nc * factor_xy
        crm_ni = crm_ni * factor_xy
        crm_ns = crm_ns * factor_xy
        crm_ng = crm_ng * factor_xy
        crm_nr = crm_nr * factor_xy


! hm 8/31/11 new output, gcm-grid- and time-step avg
! add loop over i,j do get horizontal avg, and flip vertical array
        do k=1,nzm
         l = plev-k+1
         do j=1,ny
          do i=1,nx
           aut_crm_a(l) = aut_crm_a(l) + aut1a(i,j,k)
           acc_crm_a(l) = acc_crm_a(l) + acc1a(i,j,k)
           evpc_crm_a(l) = evpc_crm_a(l) + evpc1a(i,j,k)
           evpr_crm_a(l) = evpr_crm_a(l) + evpr1a(i,j,k)
           mlt_crm_a(l) = mlt_crm_a(l) + mlt1a(i,j,k)
           sub_crm_a(l) = sub_crm_a(l) + sub1a(i,j,k)
           dep_crm_a(l) = dep_crm_a(l) + dep1a(i,j,k)
           con_crm_a(l) = con_crm_a(l) + con1a(i,j,k)
	  end do
	 end do
	end do

! note, rates are divded by dt to get mean rate over step
        aut_crm_a = aut_crm_a / dble(nstop) * factor_xy / dt
        acc_crm_a = acc_crm_a / dble(nstop) * factor_xy / dt
        evpc_crm_a = evpc_crm_a / dble(nstop) * factor_xy / dt
        evpr_crm_a = evpr_crm_a / dble(nstop) * factor_xy / dt
        mlt_crm_a = mlt_crm_a / dble(nstop) * factor_xy / dt
        sub_crm_a = sub_crm_a / dble(nstop) * factor_xy / dt
        dep_crm_a = dep_crm_a / dble(nstop) * factor_xy / dt
        con_crm_a = con_crm_a / dble(nstop) * factor_xy / dt

#endif
        precc = 0.
        precl = 0.
        precsc = 0.
        precsl = 0.
        do j=1,ny
         do i=1,nx
#ifdef sam1mom
          precsfc(i,j) = precsfc(i,j)*dz/dt/dble(nstop)
          precssfc(i,j) = precssfc(i,j)*dz/dt/dble(nstop)
#endif
#ifdef m2005 
! precsfc and precssfc from the subroutine of micro_proc in M2005 have a unit mm/s/dz
!          precsfc(i,j) = precsfc(i,j)*dz/dble(nstop)     !mm/s/dz --> mm/s
!          precssfc(i,j) = precssfc(i,j)*dz/dble(nstop)   !mm/s/dz --> mm/s
! precsfc and precssfc from the subroutine of micro_proc in M2005 have a unit mm/dz
          precsfc(i,j) = precsfc(i,j)*dz/dt/dble(nstop)     !mm/s/dz --> mm/s
          precssfc(i,j) = precssfc(i,j)*dz/dt/dble(nstop)   !mm/s/dz --> mm/s

#endif
          if(precsfc(i,j).gt.10./86400.) then
             precc = precc + precsfc(i,j)
             precsc = precsc + precssfc(i,j)
          else
             precl = precl + precsfc(i,j)
             precsl = precsl + precssfc(i,j)
          end if
         end do
        end do
        prec_crm = precsfc/1000.           !mm/s --> m/s
        precc = precc*factor_xy/1000.     
        precl = precl*factor_xy/1000.
        precsc = precsc*factor_xy/1000.
        precsl = precsl*factor_xy/1000.

!+++mhwangtest
! test water conservtion problem
      do k=1, nzm
         l=plev-k+1
         do j=1, ny
          do i=1, nx
#ifdef m2005
            qtot(9) = qtot(9)+((micro_field(i,j,k,iqr)+micro_field(i,j,k,iqs)+micro_field(i,j,k,iqg)) * pdel(l)/ggr)/(nx*ny)
            qtot(9) = qtot(9)+((micro_field(i,j,k,iqv)+micro_field(i,j,k,iqci)) * pdel(l)/ggr)/(nx*ny)
#endif
#ifdef sam1mom
            qtot(9) = qtot(9)+((micro_field(i,j,k,1)+micro_field(i,j,k,2)) * pdel(l)/ggr)/(nx*ny)
#endif
          enddo
         enddo
        enddo
        qtot(9) = qtot(9) + (precc+precl)*1000 * dt_gl

        if(abs(qtot(9)-qtot(1))/qtot(1).gt.1.0e-6) then
!           write(0, *) 'in crm water middle       ', igstep, lchnk, icol, qtot(2:8)/ntotal_step, (qtot(5)-qtot(4)) * ntotal_step/qtot(4),  &
!                                                     (qtot(6)+(precc+precl)*1000 * dt_gl-qtot(5))*ntotal_step/qtot(5)
!           write(0, *) 'in crm water middle2       ', igstep, lchnk, icol, qtot(2:8)/ntotal_step, (qtot(8)-qtot(7)) * ntotal_step/qtot(7) 
!           write(0, *) 'total water (liquid+vapor)', qtot(16:19)/nstop, (qtot(17)-qtot(16)) * ntotal_step/qtot(16), &
!                                                     (qtot(18)-qtot(19)) * ntotal_step/qtot(19),
!           call endrun('water conservation in crm.F90')
        end if
!---mhwangtest
        
        cltot = cltot *factor_xy/nstop
        clhgh = clhgh *factor_xy/nstop
        clmed = clmed *factor_xy/nstop
        cllow = cllow *factor_xy/nstop

        jt_crm = plev * 1.0
        mx_crm = 1.0
        do k=1, plev 
         mu_crm(k)=0.5*(mui_crm(k)+mui_crm(k+1))
         md_crm(k)=0.5*(mdi_crm(k)+mdi_crm(k+1))
         mu_crm(k)=mu_crm(k)*ggr/100.          !kg/m2/s --> mb/s
         md_crm(k)=md_crm(k)*ggr/100.          !kg/m2/s --> mb/s
         eu_crm(k) = 0.
         if(mui_crm(k)-mui_crm(k+1).gt.0) then
           eu_crm(k)=(mui_crm(k)-mui_crm(k+1))*ggr/pdel(k)    !/s
         else
           du_crm(k)=-1.0*(mui_crm(k)-mui_crm(k+1))*ggr/pdel(k)   !/s
         end if
         if(mdi_crm(k+1)-mdi_crm(k).lt.0) then
           ed_crm(k)=(mdi_crm(k)-mdi_crm(k+1))*ggr/pdel(k) ! /s
         else
           dd_crm(k)=-1.*(mdi_crm(k)-mdi_crm(k+1))*ggr/pdel(k)   !/s 
         end if
         if(abs(mu_crm(k)).gt.1.0e-15.or.abs(md_crm(k)).gt.1.0e-15) then
           jt_crm = min(k*1.0_r8, jt_crm)
           mx_crm = max(k*1.0_r8, mx_crm)
         end if
        end do
        
!-------------------------------------------------------------
!       Fluxes and other stat:
!-------------------------------------------------------------
        do k=1,nzm
          u2z = 0.
          v2z = 0.
          w2z = 0.
          do j=1,ny
           do i=1,nx
             u2z = u2z+(u(i,j,k)-u0(k))**2
             v2z = v2z+(v(i,j,k)-v0(k))**2
             w2z = w2z+0.5*(w(i,j,k+1)**2+w(i,j,k)**2)
           end do
          end do

!+++mhwang
! mkwsb, mkle, mkadv, mkdiff (also flux_u, flux_v) seem not calculted correclty in the spcam3.5 codes. 
! Only values at the last time step are calculated, but is averaged over the entire GCM 
! time step. 
!---mhwang

          tmp1 = dz/rhow(k)
          tmp2 = tmp1/dtn                        ! dtn is calculated inside of the icyc loop. 
                                                 ! It seems wrong to use it here ???? +++mhwang
          mkwsb(k,:) = mkwsb(k,:) * tmp1*rhow(k) * factor_xy/nstop     !kg/m3/s --> kg/m2/s
          mkwle(k,:) = mkwle(k,:) * tmp2*rhow(k) * factor_xy/nstop     !kg/m3   --> kg/m2/s  
          mkadv(k,:) = mkadv(k,:) * factor_xy*idt_gl     ! kg/kg  --> kg/kg/s
          mkdiff(k,:) = mkdiff(k,:) * factor_xy*idt_gl   ! kg/kg  --> kg/kg/s

! qpsrc, qpevp, qpfall in M2005 are calculated in micro_flux. 
          qpsrc(k) = qpsrc(k) * factor_xy*idt_gl
          qpevp(k) = qpevp(k) * factor_xy*idt_gl
          qpfall(k) = qpfall(k) * factor_xy*idt_gl   ! kg/kg in M2005 ---> kg/kg/s 
          precflux(k) = precflux(k) * factor_xy*dz/dt/nstop  !kg/m2/dz in M2005 -->kg/m2/s or mm/s (idt_gl=1/dt/nstop)

          l = plev-k+1
          flux_u(l) = (uwle(k) + uwsb(k))*tmp1*factor_xy/nstop
          flux_v(l) = (vwle(k) + vwsb(k))*tmp1*factor_xy/nstop
#ifdef sam1mom
          flux_qt(l) = mkwle(k,1) + mkwsb(k,1)
          fluxsgs_qt(l) =  mkwsb(k,1)
          flux_qp(l) = mkwle(k,2) + mkwsb(k,2)
          qt_trans(l) = mkadv(k,1) + mkdiff(k,1)
          qp_trans(l) = mkadv(k,2) + mkdiff(k,2)
#endif
#ifdef m2005
          flux_qt(l) = mkwle(k,1) + mkwsb(k,1) +  &
                   mkwle(k,iqci) + mkwsb(k,iqci)
          fluxsgs_qt(l) =  mkwsb(k,1) + mkwsb(k,iqci)
          flux_qp(l) = mkwle(k,iqr) + mkwsb(k,iqr) +  &
                   mkwle(k,iqs) + mkwsb(k,iqs) + mkwle(k,iqg) + mkwsb(k,iqg)
          qt_trans(l) = mkadv(k,1) + mkadv(k,iqci) + &
                   mkdiff(k,1) + mkdiff(k,iqci) 
          qp_trans(l) = mkadv(k,iqr) + mkadv(k,iqs) + mkadv(k,iqg) + &
                   mkdiff(k,iqr) + mkdiff(k,iqs) + mkdiff(k,iqg) 
#endif
          tkesgsz(l)= rho(k)*sum(tke(1:nx,1:ny,k))*factor_xy
          tkez(l)= rho(k)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + tkesgsz(l)
          tkz(l) = sum(tk(1:nx, 1:ny, k)) * factor_xy
          pflx(l) = precflux(k)/1000.       !mm/s  -->m/s

          qp_fall(l) = qpfall(k)
          qp_evp(l) = qpevp(k)
          qp_src(l) = qpsrc(k)

          qt_ls(l) = qtend(k)
          t_ls(l) = ttend(k)
        end do

#ifdef ECPP
        abnd=0.0
        abnd_tf=0.0
        massflxbnd=0.0
        acen=0.0
        acen_tf=0.0
        rhcen=0.0
        qcloudcen=0.0
        qicecen=0.0
        qlsinkcen=0.0
        precrcen=0.0
        precsolidcen=0.0
        wupthresh_bnd = 0.0
        wdownthresh_bnd = 0.0
        wwqui_cen = 0.0
        wwqui_bnd = 0.0
        wwqui_cloudy_cen = 0.0
        wwqui_cloudy_bnd = 0.0
        qlsink_bfcen = 0.0
        qlsink_avgcen = 0.0
        praincen = 0.0
! default is clear, non-precipitating, and quiescent class
        abnd(:,1,1,1)=1.0 
        abnd_tf(:,1,1,1)=1.0
        acen(:,1,1,1)=1.0
        acen_tf(:,1,1,1)=1.0

        do k=1, nzm
          l=plev-k+1
          acen(l,:,:,:)=area_cen_sum(k,:,1:ncls_ecpp_in,:)
          acen_tf(l,:,:,:)=area_cen_final(k,:,1:ncls_ecpp_in,:)
          rhcen(l,:,:,:)=rh_cen_sum(k,:,1:ncls_ecpp_in,:)
          qcloudcen(l,:,:,:)=qcloud_cen_sum(k,:,1:ncls_ecpp_in,:)
          qicecen(l,:,:,:)=qice_cen_sum(k,:,1:ncls_ecpp_in,:)
          qlsinkcen(l,:,:,:)=qlsink_cen_sum(k,:,1:ncls_ecpp_in,:)
          precrcen(l,:,:,:)=precr_cen_sum(k,:,1:ncls_ecpp_in,:)
          precsolidcen(l,:,:,:)=precsolid_cen_sum(k,:,1:ncls_ecpp_in,:)
          wwqui_cen(l) = wwqui_cen_sum(k)
          wwqui_cloudy_cen(l) = wwqui_cloudy_cen_sum(k)
          qlsink_bfcen(l,:,:,:)=qlsink_bf_cen_sum(k,:,1:ncls_ecpp_in,:)
          qlsink_avgcen(l,:,:,:)=qlsink_avg_cen_sum(k,:,1:ncls_ecpp_in,:)
          praincen(l,:,:,:)=prain_cen_sum(k,:,1:ncls_ecpp_in,:)
        end do
        do k=1, nzm+1
          l=plev+1-k+1
          abnd(l,:,:,:)=area_bnd_sum(k,:,1:ncls_ecpp_in,:)
          abnd_tf(l,:,:,:)=area_bnd_final(k,:,1:ncls_ecpp_in,:)
          massflxbnd(l,:,:,:)=mass_bnd_sum(k,:,1:ncls_ecpp_in,:)
          wupthresh_bnd(l)=wup_thresh(k)
          wdownthresh_bnd(l)=wdown_thresh(k)
          wwqui_bnd(l) = wwqui_bnd_sum(k)
          wwqui_cloudy_bnd(l) = wwqui_cloudy_bnd_sum(k)
        end do
#endif /*ECPP*/
        
        timing_factor = timing_factor / nstop

#ifdef SPCAM_CLUBB_SGS
! Deallocate CLUBB variables, etc.
! -UWM
        if ( doclubb .or. doclubbnoninter ) call clubb_sgs_cleanup( )
#endif
#ifdef ECPP
! Deallocate ECPP variables
       call ecpp_crm_cleanup ()
#endif /*ECPP*/
        
end subroutine crm
end module crmx_crm_module
