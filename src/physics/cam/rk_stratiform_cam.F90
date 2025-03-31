module rk_stratiform_cam

!-------------------------------------------------------------------------------------------------------
!
! Provides the CAM interface to the Rasch and Kristjansson (RK)
! prognostic cloud microphysics, and the cam4 macrophysics.
!
!-------------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use ppgrid,         only: pcols, pver, pverp
use physconst,      only: gravit, latvap, latice
use phys_control,   only: phys_getopts
use constituents,   only: pcnst
use spmd_utils,     only: masterproc

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use perf_mod

implicit none
private
save

public :: rk_stratiform_cam_register, rk_stratiform_cam_init_cnst, rk_stratiform_cam_implements_cnst
public :: rk_stratiform_cam_init
public :: rk_stratiform_cam_tend
public :: rk_stratiform_cam_readnl

! Physics buffer indices
integer  ::  landm_idx          = 0

integer  ::  qcwat_idx          = 0
integer  ::  lcwat_idx          = 0
integer  ::  tcwat_idx          = 0

integer  ::  cld_idx            = 0
integer  ::  ast_idx            = 0
integer  ::  concld_idx         = 0
integer  ::  fice_idx           = 0

integer  ::  qme_idx            = 0
integer  ::  prain_idx          = 0
integer  ::  nevapr_idx         = 0

integer  ::  wsedl_idx          = 0

integer  ::  rei_idx            = 0
integer  ::  rel_idx            = 0

integer  ::  shfrc_idx          = 0
integer  ::  cmfmc_sh_idx       = 0

integer  ::  prec_str_idx       = 0
integer  ::  snow_str_idx       = 0
integer  ::  prec_sed_idx       = 0
integer  ::  snow_sed_idx       = 0
integer  ::  prec_pcw_idx       = 0
integer  ::  snow_pcw_idx       = 0

integer  ::  ls_flxprc_idx      = 0
integer  ::  ls_flxsnw_idx      = 0

! Physics buffer indices for convective_cloud_cover
integer :: sh_frac_idx   = 0
integer :: dp_frac_idx   = 0

integer, parameter :: ncnst = 2                    ! Number of constituents
character(len=8), dimension(ncnst), parameter &    ! Constituent names
                   :: cnst_names = (/'CLDLIQ', 'CLDICE'/)
logical            :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc

logical            :: do_cnst = .false. ! True when this module has registered constituents.

integer :: &
   ixcldliq,     &! cloud liquid amount index
   ixcldice       ! cloud ice amount index

real(r8), parameter :: unset_r8 = huge(1.0_r8)
logical :: do_psrhmin

!===============================================================================
contains
!===============================================================================
  subroutine rk_stratiform_cam_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit

   use prognostic_cloud_water, only: prognostic_cloud_water_init
   use physconst,       only: tmelt, rhodair, pi

   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'rk_stratiform_readnl'

   character(len=512) :: errmsg
   integer            :: errflg

   ! Namelist variables
   real(r8) :: rk_strat_icritw  = unset_r8    !   icritw  = threshold for autoconversion of warm ice
   real(r8) :: rk_strat_icritc  = unset_r8    !   icritc  = threshold for autoconversion of cold ice
   real(r8) :: rk_strat_conke   = unset_r8    !   conke   = tunable constant for evaporation of precip
   real(r8) :: rk_strat_r3lcrit = unset_r8    !   r3lcrit = critical radius where liq conversion begins
   real(r8) :: rk_strat_polstrat_rhmin = unset_r8 ! condensation threadhold in polar stratosphere

   namelist /rk_stratiform_nl/ rk_strat_icritw, rk_strat_icritc, rk_strat_conke, rk_strat_r3lcrit
   namelist /rk_stratiform_nl/ rk_strat_polstrat_rhmin

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rk_stratiform_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rk_stratiform_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(rk_strat_icritw,  1, mpir8, 0, mpicom)
   call mpibcast(rk_strat_icritc,  1, mpir8, 0, mpicom)
   call mpibcast(rk_strat_conke,   1, mpir8, 0, mpicom)
   call mpibcast(rk_strat_r3lcrit, 1, mpir8, 0, mpicom)
   call mpibcast(rk_strat_polstrat_rhmin, 1, mpir8, 0, mpicom)
#endif

   do_psrhmin = rk_strat_polstrat_rhmin .ne. unset_r8

   call prognostic_cloud_water_init(&
      amIRoot = masterproc, &
      iulog   = iulog, &
      tmelt   = tmelt, &
      rhodair = rhodair, &
      pi      = pi, &
      icritc_in = rk_strat_icritc, &
      icritw_in = rk_strat_icritw, &
      conke_in  = rk_strat_conke,  &
      r3lcrit_in = rk_strat_r3lcrit, &
      do_psrhmin_in = do_psrhmin, &
      psrhmin_in = rk_strat_polstrat_rhmin, &
      errmsg = errmsg, &
      errflg = errflg)

   if(errflg /= 0) then
      call endrun(subname // errmsg)
   endif

   ! cloud_particle_sedimentation is initialized in pkg_cld_sediment,
   ! which has its own namelist read.

end subroutine rk_stratiform_cam_readnl

subroutine rk_stratiform_cam_register

   !---------------------------------------------------------------------- !
   !                                                                       !
   ! Register the constituents (cloud liquid and cloud ice) and the fields !
   ! in the physics buffer.                                                !
   !                                                                       !
   !---------------------------------------------------------------------- !

   use constituents, only: cnst_add, pcnst
   use physconst,    only: mwh2o, cpair

   use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls

   !-----------------------------------------------------------------------

   ! Take note of the fact that we are registering constituents.
   do_cnst = .true.

   ! Register cloud water and save indices.
   call cnst_add(cnst_names(1), mwh2o, cpair, 0._r8, ixcldliq, &
      longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
   call cnst_add(cnst_names(2), mwh2o, cpair, 0._r8, ixcldice, &
      longname='Grid box averaged cloud ice amount', is_convtran1=.true.)

   call pbuf_add_field('QCWAT',  'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), qcwat_idx)
   call pbuf_add_field('LCWAT',  'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), lcwat_idx)
   call pbuf_add_field('TCWAT',  'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), tcwat_idx)

   call pbuf_add_field('CLD',    'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), cld_idx) ! cloud_area_fraction
   call pbuf_add_field('AST',    'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), ast_idx) ! stratiform_cloud_area_fraction
   call pbuf_add_field('CONCLD', 'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), concld_idx) ! convective_cloud_area_fraction

   call pbuf_add_field('FICE',   'physpkg', dtype_r8, (/pcols,pver/), fice_idx) ! mass_fraction_of_ice_content_within_stratiform_cloud

   call pbuf_add_field('QME',       'physpkg', dtype_r8, (/pcols,pver/), qme_idx) ! net_condensation_rate_due_to_microphysics
   call pbuf_add_field('PRAIN',     'physpkg', dtype_r8, (/pcols,pver/), prain_idx) ! precipitation_production_due_to_microphysics
   call pbuf_add_field('NEVAPR',    'physpkg', dtype_r8, (/pcols,pver/), nevapr_idx) ! precipitation_evaporation_due_to_microphysics

   call pbuf_add_field('WSEDL',     'physpkg', dtype_r8, (/pcols,pver/), wsedl_idx) ! sedimentation velocity of liquid stratus cloud droplet [m s-1]

   call pbuf_add_field('REI',       'physpkg', dtype_r8, (/pcols,pver/), rei_idx) ! effective_radius_of_stratiform_cloud_ice_particle
   call pbuf_add_field('REL',       'physpkg', dtype_r8, (/pcols,pver/), rel_idx) ! effective_radius_of_stratiform_cloud_liquid_water_particle

   call pbuf_add_field('LS_FLXPRC', 'physpkg', dtype_r8, (/pcols,pverp/), ls_flxprc_idx) ! effective_radius_of_stratiform_cloud_liquid_water_particle
   call pbuf_add_field('LS_FLXSNW', 'physpkg', dtype_r8, (/pcols,pverp/), ls_flxsnw_idx) ! stratiform_snow_flux_at_interface

end subroutine rk_stratiform_cam_register

!===============================================================================

function rk_stratiform_cam_implements_cnst(name)

  !----------------------------------------------------------------------------- !
  !                                                                              !
  ! Return true if specified constituent is implemented by this package          !
  !                                                                              !
  !----------------------------------------------------------------------------- !

   character(len=*), intent(in) :: name      ! constituent name
   logical :: rk_stratiform_cam_implements_cnst     ! return value

   !-----------------------------------------------------------------------

   rk_stratiform_cam_implements_cnst = (do_cnst .and. any(name == cnst_names))

end function rk_stratiform_cam_implements_cnst

!===============================================================================

subroutine rk_stratiform_cam_init_cnst(name, latvals, lonvals, mask, q)

   !----------------------------------------------------------------------- !
   !                                                                        !
   ! Initialize the cloud water mixing ratios (liquid and ice), if they are !
   ! not read from the initial file                                         !
   !                                                                        !
   !----------------------------------------------------------------------- !

   character(len=*), intent(in)  :: name       ! constituent name
   real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
   real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
   logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev
   !-----------------------------------------------------------------------
   integer :: k

   if (any(name == cnst_names)) then
     do k = 1, size(q, 2)
       where(mask)
         q(:, k) = 0.0_r8
       end where
     end do
   end if

end subroutine rk_stratiform_cam_init_cnst

!===============================================================================

subroutine rk_stratiform_cam_init()

   !-------------------------------------------- !
   !                                             !
   ! Initialize the cloud water parameterization !
   !                                             !
   !-------------------------------------------- !

   use physics_buffer,  only: physics_buffer_desc, pbuf_get_index
   use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, sflxnam, apcnst, bpcnst
   use cam_history,     only: addfld, add_default, horiz_only
   use convect_shallow, only: convect_shallow_use_shfrc
   use phys_control,    only: cam_physpkg_is
   use physconst,       only: tmelt, rhodair, rh2o

   integer :: m, mm
   logical :: history_amwg         ! output the variables used by the AMWG diag package
   logical :: history_aerosol      ! Output the MAM aerosol tendencies
   logical :: history_budget       ! Output tendencies and state variables for CAM4
                                   ! temperature, water vapor, cloud ice and cloud
                                   ! liquid budgets.
   integer :: history_budget_histfile_num ! output history file number for budget fields
   !-----------------------------------------------------------------------

   call phys_getopts( history_aerosol_out        = history_aerosol      , &
                      history_amwg_out   = history_amwg                 , &
                      history_budget_out         = history_budget       , &
                      history_budget_histfile_num_out = history_budget_histfile_num)

   landm_idx = pbuf_get_index("LANDM")

   ! Find out whether shfrc from convect_shallow will be used in cldfrc
   if( convect_shallow_use_shfrc() ) then
      use_shfrc = .true.
      shfrc_idx = pbuf_get_index('shfrc')
   else
      use_shfrc = .false.
   endif

   ! Register history variables

   do m = 1, ncnst
      call cnst_get_ind( cnst_names(m), mm )
      call addfld( cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg',   cnst_longname(mm)                    )
      call addfld( sflxnam  (mm), horiz_only,  'A', 'kg/m2/s', trim(cnst_name(mm))//' surface flux' )
      if (history_amwg) then
         call add_default( cnst_name(mm), 1, ' ' )
         call add_default( sflxnam  (mm), 1, ' ' )
      endif
   enddo

   call addfld (apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
   call addfld (apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
   call addfld (bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
   call addfld (bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )

   if( history_budget) then
      call add_default (cnst_name(ixcldliq), history_budget_histfile_num, ' ')
      call add_default (cnst_name(ixcldice), history_budget_histfile_num, ' ')
      call add_default (apcnst   (ixcldliq), history_budget_histfile_num, ' ')
      call add_default (apcnst   (ixcldice), history_budget_histfile_num, ' ')
      call add_default (bpcnst   (ixcldliq), history_budget_histfile_num, ' ')
      call add_default (bpcnst   (ixcldice), history_budget_histfile_num, ' ')
   end if

   call addfld ('FWAUT',      (/ 'lev' /),   'A', 'fraction', 'Relative importance of liquid autoconversion'            )
   call addfld ('FSAUT',      (/ 'lev' /),   'A', 'fraction', 'Relative importance of ice autoconversion'               )
   call addfld ('FRACW',      (/ 'lev' /),   'A', 'fraction', 'Relative importance of rain accreting liquid'            )
   call addfld ('FSACW',      (/ 'lev' /),   'A', 'fraction', 'Relative importance of snow accreting liquid'            )
   call addfld ('FSACI',      (/ 'lev' /),   'A', 'fraction', 'Relative importance of snow accreting ice'               )
   call addfld ('CME',        (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of cond-evap within the cloud'                      )
   call addfld ('CMEICE',     (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of cond-evap of ice within the cloud'               )
   call addfld ('CMELIQ',     (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of cond-evap of liq within the cloud'               )
   call addfld ('ICE2PR',     (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of conversion of ice to precip'                     )
   call addfld ('LIQ2PR',     (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of conversion of liq to precip'                     )
   call addfld ('ZMDLF',      (/ 'lev' /),   'A', 'kg/kg/s' , 'Detrained liquid water from ZM convection'               )
   call addfld ('SHDLF',      (/ 'lev' /),   'A', 'kg/kg/s' , 'Detrained liquid water from shallow convection'          )

   call addfld ('PRODPREC',   (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of conversion of condensate to precip'              )
   call addfld ('EVAPPREC',   (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of evaporation of falling precip'                   )
   call addfld ('EVAPSNOW',   (/ 'lev' /),   'A', 'kg/kg/s' , 'Rate of evaporation of falling snow'                     )
   call addfld ('HPROGCLD',   (/ 'lev' /),   'A', 'W/kg'    , 'Heating from prognostic clouds'                          )
   call addfld ('HCME',       (/ 'lev' /),   'A', 'W/kg'    , 'Heating from cond-evap within the cloud'                 )
   call addfld ('HEVAP',      (/ 'lev' /),   'A', 'W/kg'    , 'Heating from evaporation of falling precip'              )
   call addfld ('HFREEZ',     (/ 'lev' /),   'A', 'W/kg'    , 'Heating rate due to freezing of precip'                  )
   call addfld ('HMELT',      (/ 'lev' /),   'A', 'W/kg'    , 'Heating from snow melt'                                  )
   call addfld ('HREPART',    (/ 'lev' /),   'A', 'W/kg'    , 'Heating from cloud ice/liquid repartitioning'            )
   call addfld ('REPARTICE',  (/ 'lev' /),   'A', 'kg/kg/s' , 'Cloud ice tendency from cloud ice/liquid repartitioning' )
   call addfld ('REPARTLIQ',  (/ 'lev' /),   'A', 'kg/kg/s' , 'Cloud liq tendency from cloud ice/liquid repartitioning' )
   call addfld ('FICE',       (/ 'lev' /),   'A', 'fraction', 'Fractional ice content within cloud'                     )
   call addfld ('ICWMR',      (/ 'lev' /),   'A', 'kg/kg'   , 'Prognostic in-cloud water mixing ratio'                  )
   call addfld ('ICIMR',      (/ 'lev' /),   'A', 'kg/kg'   , 'Prognostic in-cloud ice mixing ratio'                    )
   call addfld ('PCSNOW',     horiz_only   , 'A', 'm/s'     , 'Snow fall from prognostic clouds'                        )

   call addfld ('DQSED',      (/ 'lev' /),   'A', 'kg/kg/s' , 'Water vapor tendency from cloud sedimentation'           )
   call addfld ('DLSED',      (/ 'lev' /),   'A', 'kg/kg/s' , 'Cloud liquid tendency from sedimentation'                )
   call addfld ('DISED',      (/ 'lev' /),   'A', 'kg/kg/s' , 'Cloud ice tendency from sedimentation'                   )
   call addfld ('HSED',       (/ 'lev' /),   'A', 'W/kg'    , 'Heating from cloud sediment evaporation'                 )
   call addfld ('SNOWSED',    horiz_only,    'A', 'm/s'     , 'Snow from cloud ice sedimentation'                       )
   call addfld ('RAINSED',    horiz_only,    'A', 'm/s'     , 'Rain from cloud liquid sedimentation'                    )
   call addfld ('PRECSED',    horiz_only,    'A', 'm/s'     , 'Precipitation from cloud sedimentation'                  )

   call addfld ('CLDST',      (/ 'lev' /),   'A', 'fraction', 'Stratus cloud fraction'                                  )
   call addfld ('CONCLD',     (/ 'lev' /),   'A', 'fraction', 'Convective cloud cover'                                  )

   call addfld ('AST',        (/ 'lev' /),   'A', 'fraction', 'Stratus cloud fraction'                                  )
   call addfld ('LIQCLDF',    (/ 'lev' /),   'A', 'fraction', 'Stratus Liquid cloud fraction'                           )
   call addfld ('ICECLDF',    (/ 'lev' /),   'A', 'fraction', 'Stratus ICE cloud fraction'                              )
   call addfld ('IWC',        (/ 'lev' /),   'A', 'kg/m3'   , 'Grid box average ice water content'                      )
   call addfld ('LWC',        (/ 'lev' /),   'A', 'kg/m3'   , 'Grid box average liquid water content'                   )
   call addfld ('REL',        (/ 'lev' /),   'A', 'micron'  , 'effective liquid drop radius'                            )
   call addfld ('REI',        (/ 'lev' /),   'A', 'micron'  , 'effective ice particle radius'                           )

   if ( history_budget ) then

      call add_default ('EVAPSNOW ', history_budget_histfile_num, ' ')
      call add_default ('EVAPPREC ', history_budget_histfile_num, ' ')
      call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')

      if( cam_physpkg_is('cam4') ) then

         call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')
         call add_default ('CME      ', history_budget_histfile_num, ' ')
         call add_default ('DQSED    ', history_budget_histfile_num, ' ')
         call add_default ('DISED    ', history_budget_histfile_num, ' ')
         call add_default ('DLSED    ', history_budget_histfile_num, ' ')
         call add_default ('HSED     ', history_budget_histfile_num, ' ')
         call add_default ('CMEICE   ', history_budget_histfile_num, ' ')
         call add_default ('LIQ2PR   ', history_budget_histfile_num, ' ')
         call add_default ('ICE2PR   ', history_budget_histfile_num, ' ')
         call add_default ('HCME     ', history_budget_histfile_num, ' ')
         call add_default ('HEVAP    ', history_budget_histfile_num, ' ')
         call add_default ('HFREEZ   ', history_budget_histfile_num, ' ')
         call add_default ('HMELT    ', history_budget_histfile_num, ' ')
         call add_default ('HREPART  ', history_budget_histfile_num, ' ')
         call add_default ('HPROGCLD ', history_budget_histfile_num, ' ')
         call add_default ('REPARTLIQ', history_budget_histfile_num, ' ')
         call add_default ('REPARTICE', history_budget_histfile_num, ' ')

      end if

   end if

   if (history_amwg) then
      call add_default ('ICWMR', 1, ' ')
      call add_default ('ICIMR', 1, ' ')
      call add_default ('CONCLD  ', 1, ' ')
      call add_default ('FICE    ', 1, ' ')
   endif

   ! History Variables for COSP/CFMIP
   call addfld ('LS_FLXPRC', (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface rain+snow flux')
   call addfld ('LS_FLXSNW', (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface snow flux')
   call addfld ('PRACWO',    (/ 'lev' /),  'A', '1/s', 'Accretion of cloud water by rain')
   call addfld ('PSACWO',    (/ 'lev' /),  'A', '1/s', 'Accretion of cloud water by snow')
   call addfld ('PSACIO',    (/ 'lev' /),  'A', '1/s', 'Accretion of cloud ice by snow')

   call addfld ('CLDLIQSTR', (/ 'lev' /),  'A', 'kg/kg', 'Stratiform CLDLIQ')
   call addfld ('CLDICESTR', (/ 'lev' /),  'A', 'kg/kg', 'Stratiform CLDICE')
   call addfld ('CLDLIQCON', (/ 'lev' /),  'A', 'kg/kg', 'Convective CLDLIQ')
   call addfld ('CLDICECON', (/ 'lev' /),  'A', 'kg/kg', 'Convective CLDICE')

   cmfmc_sh_idx = pbuf_get_index('CMFMC_SH') ! atmosphere_convective_mass_flux_due_to_shallow_convection
   prec_str_idx = pbuf_get_index('PREC_STR') ! stratiform_rain_and_snow_surface_flux_due_to_sedimentation
   snow_str_idx = pbuf_get_index('SNOW_STR') ! lwe_snow_and_cloud_ice_precipitation_rate_at_surface_due_to_microphysics
   prec_pcw_idx = pbuf_get_index('PREC_PCW') ! lwe_stratiform_precipitation_rate_at_surface
   snow_pcw_idx = pbuf_get_index('SNOW_PCW') ! lwe_snow_precipitation_rate_at_surface_due_to_microphysics
   prec_sed_idx = pbuf_get_index('PREC_SED') ! stratiform_cloud_water_surface_flux_due_to_sedimentation
   snow_sed_idx = pbuf_get_index('SNOW_SED') ! lwe_cloud_ice_sedimentation_rate_at_surface_due_to_microphysics

   sh_frac_idx  = pbuf_get_index('SH_FRAC')
   dp_frac_idx  = pbuf_get_index('DP_FRAC')

end subroutine rk_stratiform_cam_init

!===============================================================================

subroutine rk_stratiform_cam_tend( &
   state, ptend_all, pbuf, dtime, icefrac, &
   landfrac, ocnfrac, snowh, dlf,   &
   dlf2, rliq, cmfmc, ts,                  &
   sst, zdu)

   !-------------------------------------------------------- !
   !                                                         !
   ! Interface to sedimentation, detrain, cloud fraction and !
   !        cloud macro - microphysics subroutines           !
   !                                                         !
   !-------------------------------------------------------- !

   use cloud_fraction_fice,  only: cloud_fraction_fice_run
   use physics_types,    only: physics_state, physics_ptend
   use physics_types,    only: physics_ptend_init, physics_update
   use physics_types,    only: physics_ptend_sum,  physics_state_copy
   use physics_types,    only: physics_state_dealloc
   use cam_history,      only: outfld
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

   ! ccppized version of the RK scheme
   use rk_stratiform,    only: rk_stratiform_sedimentation_run
   use cloud_particle_sedimentation, only: cloud_particle_sedimentation_run
   use rk_stratiform,    only: rk_stratiform_detrain_convective_condensate_run
   use compute_cloud_fraction, only: compute_cloud_fraction_run
   use rk_stratiform,    only: rk_stratiform_cloud_fraction_perturbation_run
   use rk_stratiform,    only: rk_stratiform_condensate_repartioning_run
   use rk_stratiform,    only: rk_stratiform_external_forcings_run
   use prognostic_cloud_water, only: prognostic_cloud_water_run
   use rk_stratiform,    only: rk_stratiform_prognostic_cloud_water_tendencies_run
   use convective_cloud_cover, only: convective_cloud_cover_run
   use rk_stratiform,    only: rk_stratiform_save_qtlcwat_run
   use rk_stratiform,    only: rk_stratiform_cloud_optical_properties_run

   use prognostic_cloud_water, only: icritc !REMOVECAM no longer need to be public after CAM is retired
   use physconst,        only: gravit, rhoh2o, epsilo, latvap, latice, cpair, tmelt, cappa
   use physconst,        only: rair, rh2o, pi
   use ref_pres,         only: trop_cloud_top_lev

   use cloud_optical_properties,    only: cldefr
   use phys_control,     only: cam_physpkg_is
   use tropopause,       only: tropopause_find_cam
   use phys_grid,        only: get_lat_all_p
   use constituents,     only: cnst_get_ind

   ! Arguments
   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend_all   ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep
   real(r8), intent(in)  :: icefrac (pcols)          ! Sea ice fraction (fraction)
   real(r8), intent(in)  :: landfrac(pcols)          ! Land fraction (fraction)
   real(r8), intent(in)  :: ocnfrac (pcols)          ! Ocean fraction (fraction)
   real(r8), intent(in)  :: snowh(pcols)             ! Snow depth over land, water equivalent (m)

   real(r8), intent(in)  :: dlf(pcols,pver)          ! Detrained water from convection schemes
   real(r8), intent(in)  :: dlf2(pcols,pver)         ! Detrained water from shallow convection scheme
   real(r8), intent(in)  :: rliq(pcols)              ! Vertical integral of liquid not yet in q(ixcldliq)
   real(r8), intent(in)  :: cmfmc(pcols,pverp)       ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]

   real(r8), intent(in)  :: ts(pcols)                ! Surface temperature
   real(r8), intent(in)  :: sst(pcols)               ! Sea surface temperature
   real(r8), intent(in)  :: zdu(pcols,pver)          ! Detrainment rate from deep convection

  ! Local variables

   type(physics_state)   :: state1                   ! Local copy of the state variable
   type(physics_ptend)   :: ptend_loc                ! Package tendencies

   integer :: i, k, m
   integer :: lchnk                                  ! Chunk identifier
   integer :: ncol                                   ! Number of atmospheric columns
   integer :: itim_old

   real(r8), parameter :: rdair = 287.15_r8

   ! Physics buffer fields
   real(r8), pointer :: landm(:)             ! Land fraction ramped over water

   real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ]
   real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
   real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
   real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation
   real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
   real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]
   real(r8), pointer, dimension(:,:) :: deepcu      ! deep convection cloud fraction
   real(r8), pointer, dimension(:,:) :: shallowcu   ! shallow convection cloud fraction

   real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
   real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
   real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
   real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
   real(r8), pointer, dimension(:,:) :: fice         ! Cloud ice/water partitioning ratio.
   real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
   real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction
   real(r8), pointer, dimension(:,:) :: qme          ! rate of cond-evap of condensate (1/s)
   real(r8), pointer, dimension(:,:) :: prain        ! Total precipitation (rain + snow)
   real(r8), pointer, dimension(:,:) :: nevapr       ! Evaporation of total precipitation (rain + snow)
   real(r8), pointer, dimension(:,:) :: rel          ! Liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:) :: rei          ! Ice effective drop size (microns)
   real(r8), pointer, dimension(:,:) :: wsedl        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
   real(r8), pointer, dimension(:,:) :: shfrc        ! Cloud fraction from shallow convection scheme
   real(r8), pointer, dimension(:,:) :: cmfmc_sh     ! Shallow convective mass flux (pcols,pverp) [ kg/s/m^2 ]

   real(r8), target :: shfrc_local(pcols,pver)

   ! physics buffer fields for COSP simulator (RK only)
   real(r8), pointer, dimension(:,:) :: rkflxprc     ! RK grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: rkflxsnw     ! RK grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)

   ! Local variables for stratiform_sediment
   real(r8) :: rain(pcols)                           ! Surface flux of cloud liquid
   real(r8) :: pvliq(pcols,pverp)                    ! Vertical velocity of cloud liquid drops (Pa/s)
   real(r8) :: pvice(pcols,pverp)                    ! Vertical velocity of cloud ice particles (Pa/s)

   ! Local variables for cldfrc

   real(r8) :: cldst(pcols,pver)                       ! Stratus cloud fraction
   real(r8) :: rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
   real(r8) :: relhum(pcols,pver)                      ! RH, output to determine drh/da
   real(r8) :: rhu00(pcols,pver)
   real(r8) :: rhdfda(pcols,pver)
   real(r8) :: icecldf(pcols,pver)                     ! Ice cloud fraction
   real(r8) :: liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)

   ! Local variables for microphysics

   real(r8) :: qtend(pcols,pver)                       ! Moisture tendencies
   real(r8) :: ttend(pcols,pver)                       ! Temperature tendencies
   real(r8) :: ltend(pcols,pver)                       ! Cloud liquid water tendencies
   real(r8) :: evapheat(pcols,pver)                    ! Heating rate due to evaporation of precip
   real(r8) :: evapsnow(pcols,pver)                    ! Local evaporation of snow
   real(r8) :: prfzheat(pcols,pver)                    ! Heating rate due to freezing of precip (W/kg)
   real(r8) :: meltheat(pcols,pver)                    ! Heating rate due to phase change of precip
   real(r8) :: cmeheat (pcols,pver)                    ! Heating rate due to phase change of precip
   real(r8) :: prodsnow(pcols,pver)                    ! Local production of snow
   real(r8) :: totcw(pcols,pver)                       ! Total cloud water mixing ratio
   real(r8) :: fsnow(pcols,pver)                       ! Fractional snow production
   real(r8) :: repartht(pcols,pver)                    ! Heating rate due to phase repartition of input precip
   real(r8) :: icimr(pcols,pver)                       ! In cloud ice mixing ratio
   real(r8) :: icwmr(pcols,pver)                       ! In cloud water mixing ratio
   real(r8) :: fwaut(pcols,pver)
   real(r8) :: fsaut(pcols,pver)
   real(r8) :: fracw(pcols,pver)
   real(r8) :: fsacw(pcols,pver)
   real(r8) :: fsaci(pcols,pver)
   real(r8) :: cmeice(pcols,pver)                      ! Rate of cond-evap of ice within the cloud
   real(r8) :: cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud
   real(r8) :: ice2pr(pcols,pver)                      ! Rate of conversion of ice to precip
   real(r8) :: liq2pr(pcols,pver)                      ! Rate of conversion of liquid to precip
   real(r8) :: liq2snow(pcols,pver)                    ! Rate of conversion of liquid to snow

   ! Local variables for CFMIP calculations
   real(r8) :: mr_lsliq(pcols,pver)                     ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
   real(r8) :: mr_lsice(pcols,pver)                     ! mixing_ratio_large_scale_cloud_ice (kg/kg)
   real(r8) :: mr_ccliq(pcols,pver)                     ! mixing_ratio_convective_cloud_liquid (kg/kg)
   real(r8) :: mr_ccice(pcols,pver)                     ! mixing_ratio_convective_cloud_ice (kg/kg)

   real(r8) :: pracwo(pcols,pver)                       ! RK accretion of cloud water by rain (1/s)
   real(r8) :: psacwo(pcols,pver)                       ! RK accretion of cloud water by snow (1/s)
   real(r8) :: psacio(pcols,pver)                       ! RK accretion of cloud ice by snow (1/s)

   real(r8) :: iwc(pcols,pver)                         ! Grid box average ice water content
   real(r8) :: lwc(pcols,pver)                         ! Grid box average liquid water content

   logical  :: lq(pcnst)
   integer  :: troplev(pcols)
   real(r8) :: dlat(pcols)

   integer  :: top_lev
   integer  :: ixq

   character(len=512)   :: errmsg
   integer              :: errflg


   ! ======================================================================

   lchnk = state%lchnk
   ncol  = state%ncol

   call cnst_get_ind('Q', ixq) ! water vapor index in const array.

   call physics_state_copy(state,state1)             ! Copy state to local state1.

   ! Associate pointers with physics buffer fields

   call pbuf_get_field(pbuf, landm_idx,   landm) ! smoothed_land_area_fraction

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, qcwat_idx,   qcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, tcwat_idx,   tcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, lcwat_idx,   lcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, fice_idx,    fice)

   call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)

   call pbuf_get_field(pbuf, prec_str_idx, prec_str)
   call pbuf_get_field(pbuf, snow_str_idx, snow_str)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

   call pbuf_get_field(pbuf, qme_idx,    qme )
   call pbuf_get_field(pbuf, prain_idx,  prain)
   call pbuf_get_field(pbuf, nevapr_idx, nevapr)

   call pbuf_get_field(pbuf, rel_idx,    rel)
   call pbuf_get_field(pbuf, rei_idx,    rei)

   call pbuf_get_field(pbuf, wsedl_idx,  wsedl)

   ! check that qcwat and tcwat were initialized; if not then do it now.
   if (qcwat(1,1) == huge(1._r8)) then
      qcwat(:ncol,:) = state%q(:ncol,:,1)
   end if
   if (tcwat(1,1) == huge(1._r8)) then
      tcwat(:ncol,:) = state%t(:ncol,:)
   end if

   if ( do_psrhmin ) then
      !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
      troplev(:) = 0
      !REMOVECAM_END
      call tropopause_find_cam(state, troplev)
   endif

   ! dlat required for CCPPized scheme
   call get_lat_all_p(lchnk,ncol,dlat)

   ! ------------- !
   ! Sedimentation !
   ! ------------- !

   ! Allow the cloud liquid drops and ice particles to sediment.
   ! This is done before adding convectively detrained cloud water,
   ! because the phase of the detrained water is unknown.

   call t_startf('stratiform_sediment')

   lq(:)        = .FALSE.
   lq(ixq)        = .TRUE.
   lq(ixcldice) = .TRUE.
   lq(ixcldliq) = .TRUE.
   call physics_ptend_init(ptend_loc, state%psetcols, 'pcwsediment', ls=.true., lq=lq)! Initialize local ptend type

!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   pvliq(:,:) = 0._r8
   pvice(:,:) = 0._r8
   rain(:) = 0._r8
   snow_sed(:) = 0._r8
!REMOVECAM_END
   ! Call CCPP-ized subroutine for cld_sediment_vel and cld_sediment_tend combined
   call cloud_particle_sedimentation_run( &
     ncol                = ncol,                           &
     pver                = pver,                           &
     pverp               = pverp,                          &
     dtime               = dtime,                          &
     tmelt               = tmelt,                          &
     gravit              = gravit,                         &
     latvap              = latvap,                         &
     latice              = latice,                         &
     rair                = rair,                           &
     rhoh2o              = rhoh2o,                         &
     icritc              = icritc,                         &
     pint                = state1%pint(:ncol,:),           &
     pmid                = state1%pmid(:ncol,:),           &
     pdel                = state1%pdel(:ncol,:),           &
     t                   = state1%t(:ncol,:),              &
     cloud               = cld(:ncol,:),                   &
     icefrac             = icefrac(:ncol),                 &
     landfrac            = landfrac(:ncol),                &
     ocnfrac             = ocnfrac(:ncol),                 &
     cldliq              = state1%q(:ncol,:,ixcldliq),     &
     cldice              = state1%q(:ncol,:,ixcldice),     &
     snowh               = snowh(:ncol),                    &
     landm               = landm(:ncol),                    &
     ! Outputs
     pvliq               = pvliq(:ncol,:),                 &
     pvice               = pvice(:ncol,:),                 &
     liqtend             = ptend_loc%q(:ncol,:,ixcldliq),  &
     icetend             = ptend_loc%q(:ncol,:,ixcldice),  &
     wvtend              = ptend_loc%q(:ncol,:,ixq),       &
     htend               = ptend_loc%s(:ncol,:),           &
     sfliq               = rain(:ncol),                    & ! mass units (not precip units)
     sfice               = snow_sed(:ncol),                & ! precip units (as stored in pbuf)
     errmsg              = errmsg,                         &
     errflg              = errflg)

   if(errflg /= 0) then
      call endrun('rk_stratiform_tend:' // errmsg)
   endif

   wsedl(:ncol,:pver) = pvliq(:ncol,:pver)/gravit/(state1%pmid(:ncol,:pver)/(rdair*state1%t(:ncol,:pver)))

!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   prec_sed(:) = 0._r8
   prec_str(:) = 0._r8
   snow_str(:) = 0._r8
!REMOVECAM_END
   ! Call the CCPPized scheme for accumulation of total stratiform fluxes [m s-1]
   call rk_stratiform_sedimentation_run( &
     ncol                = ncol,                           &
     sfliq               = rain(:ncol),                    & ! in from sedimentation
     snow_sed            = snow_sed(:ncol),                & ! in from sedimentation
     prec_sed            = prec_sed(:ncol),                & ! out -- converted to precip units
     prec_str            = prec_str(:ncol),                & ! out -- initial assignment for accumulation
     snow_str            = snow_str(:ncol),                & ! out -- initial assignment for accumulation
     errmsg              = errmsg,                         &
     errflg              = errflg)

   if(errflg /= 0) then
      call endrun('rk_stratiform_tend:' // errmsg)
   endif

   ! Record history variables
   call outfld( 'DQSED'   ,ptend_loc%q(:,:,ixq)     , pcols,lchnk )
   call outfld( 'DISED'   ,ptend_loc%q(:,:,ixcldice), pcols,lchnk )
   call outfld( 'DLSED'   ,ptend_loc%q(:,:,ixcldliq), pcols,lchnk )
   call outfld( 'HSED'    ,ptend_loc%s              , pcols,lchnk )
   call outfld( 'PRECSED' ,prec_sed                 , pcols,lchnk )
   call outfld( 'SNOWSED' ,snow_sed                 , pcols,lchnk )
   call outfld( 'RAINSED' ,rain/1000._r8            , pcols,lchnk ) ! convert from kg m-2 s-1 to m s-1 (precip units) for output

   ! Add tendency from this process to tend from other processes here
   call physics_ptend_init(ptend_all, state%psetcols, 'stratiform')
   call physics_ptend_sum( ptend_loc, ptend_all, ncol )

   ! Update physics state type state1 with ptend_loc
   call physics_update( state1, ptend_loc, dtime )

   call t_stopf('stratiform_sediment')

   ! ----------------------------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment or stratiform cloud !
   ! ----------------------------------------------------------------------------- !

   ! Put all of the detraining cloud water from convection into the large scale cloud.
   ! It all goes in liquid for the moment.
   ! Strictly speaking, this approach is detraining all the cconvective water into
   ! the environment, not the large-scale cloud.

   lq(:)        = .FALSE.
   lq(ixcldliq) = .TRUE.
   call physics_ptend_init( ptend_loc, state1%psetcols, 'pcwdetrain', lq=lq)

   call rk_stratiform_detrain_convective_condensate_run( &
      ncol        = ncol,                          &
      dlf         = dlf(:ncol,:),                  &
      rliq        = rliq(:ncol),                   &
      prec_str    = prec_str(:ncol),               &
      tend_cldliq = ptend_loc%q(:ncol,:,ixcldliq), &
      errmsg      = errmsg, &
      errflg      = errflg)

   if(errflg /= 0) then
      call endrun('rk_stratiform_tend:' // errmsg)
   endif

   call outfld( 'ZMDLF', dlf, pcols, lchnk )
   call outfld( 'SHDLF', dlf2, pcols, lchnk )

   ! Add hie detrainment tendency to tend from the other prior processes

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_update( state1, ptend_loc, dtime )

   ! -------------------------------------- !
   ! Computation of Various Cloud Fractions !
   ! -------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Treatment of cloud fraction in CAM4 and CAM5 differs                          !
   ! (1) CAM4                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( empirical fcn of mass flux )        !
   !     . Stratus AMT = max( RH stratus AMT, Stability Stratus AMT )              !
   !     . Cumulus and Stratus are 'minimally' overlapped without hierarchy.       !
   !     . Cumulus LWC,IWC is assumed to be the same as Stratus LWC,IWC            !
   ! (2) CAM5                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) !
   !     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     !
   !     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  !
   !     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              !
   ! ----------------------------------------------------------------------------- !

   if( use_shfrc ) then
       call pbuf_get_field(pbuf, shfrc_idx, shfrc )
   else
       shfrc=>shfrc_local
       shfrc(:,:) = 0._r8
   endif

   call pbuf_get_field(pbuf, sh_frac_idx, shallowcu )
   call pbuf_get_field(pbuf, dp_frac_idx, deepcu )

   !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   shallowcu(:,:) = 0._r8
   deepcu(:,:) = 0._r8
   concld(:,:) = 0._r8
   !REMOVECAM_END

   ! compute convective cloud fraction using CCPP-ized subroutine
   call convective_cloud_cover_run( &
      ncol        = ncol, &
      pver        = pver, &
      top_lev_cloudphys = 1, & ! CAM4 macrophysics.
      use_shfrc   = use_shfrc, &
      shfrc       = shfrc(:ncol,:), &
      cmfmc_total = cmfmc(:ncol,:), &
      cmfmc_sh    = cmfmc_sh(:ncol,:), &
      shallowcu   = shallowcu(:ncol,:), &
      deepcu      = deepcu(:ncol,:), &
      concld      = concld(:ncol,:), &
      errmsg      = errmsg, &
      errflg      = errflg)

   ! write out convective cloud fraction diagnostic.
   call outfld( 'SH_CLD  ', shallowcu   , pcols, lchnk )
   call outfld( 'DP_CLD  ', deepcu      , pcols, lchnk )
   call outfld( 'CONCLD  ', concld      , pcols, lchnk )

   ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
   ! will be computed using this updated 'concld' in the stratiform macrophysics
   ! scheme (mmacro_pcond) later below.

   call t_startf("cldfrc")

   !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   cld(:,:) = 0._r8
   rhcloud(:,:) = 0._r8
   cldst(:,:) = 0._r8
   rhu00(:,:) = 0._r8
   icecldf(:,:) = 0._r8
   liqcldf(:,:) = 0._r8
   relhum(:,:) = 0._r8
   !REMOVECAM_END

   ! call CCPPized cloud fraction scheme
   call compute_cloud_fraction_run( &
      ncol = ncol, &
      pver = pver, &
      cappa = cappa, &
      gravit = gravit, &
      rair = rair, &
      tmelt = tmelt, &
      top_lev_cloudphys = 1, & ! CAM4 macrophysics.
      pmid = state1%pmid(:ncol,:), &
      ps = state1%pint(:,pverp), &
      temp = state1%t(:ncol,:), &
      sst = sst(:ncol), &
      q = state1%q(:ncol,:,ixq), &
      cldice = state1%q(:ncol,:,ixcldice), &
      phis = state1%phis(:ncol), &
      shallowcu = shallowcu(:ncol,:), &
      deepcu = deepcu(:ncol,:), &
      concld = concld(:ncol,:), &
      landfrac = landfrac(:ncol), &
      ocnfrac = ocnfrac(:ncol), &
      snowh = snowh(:ncol), &
      rhpert_flag = .false., & ! below output:
      cloud = cld(:ncol, :), &
      rhcloud = rhcloud(:ncol, :), &
      cldst = cldst(:ncol,:), &
      rhu00 = rhu00(:ncol,:), &
      icecldf = icecldf(:ncol,:), &
      liqcldf = liqcldf(:ncol,:), &
      relhum = relhum(:ncol,:), &
      errmsg = errmsg, &
      errflg = errflg)

   ! The AST history output is from cld from this call of compute_cloud_fraction,
   ! and not in the subsequent calls. (hplin, 3/6/25)
   call outfld( 'AST',      cld,    pcols, lchnk )

   ! Re-calculate cloud with perturbed rh add call cldfrc to estimate rhdfda.
   !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   rhdfda(:,:) = 0._r8
   !REMOVECAM_END

   ! Re-calculate cloud with perturbed rh add call cldfrc to estimate rhdfda.
   call rk_stratiform_cloud_fraction_perturbation_run( &
      ncol = ncol, &
      pver = pver, &
      cappa = cappa, &
      gravit = gravit, &
      rair = rair, &
      tmelt = tmelt, &
      top_lev_cloudphys = 1, & ! CAM4 macrophysics.
      pmid = state1%pmid(:ncol,:), &
      ps = state1%pint(:,pverp), &
      temp = state1%t(:ncol,:), &
      sst = sst(:ncol), &
      q_wv = state1%q(:ncol,:,ixq), &
      cldice = state1%q(:ncol,:,ixcldice), &
      phis = state1%phis(:ncol), &
      shallowcu = shallowcu(:ncol,:), &
      deepcu = deepcu(:ncol,:), &
      concld = concld(:ncol,:), &
      landfrac = landfrac(:ncol), &
      ocnfrac = ocnfrac(:ncol), &
      snowh = snowh(:ncol), &
      cloud = cld(:ncol,:), & ! unperturbed input
      relhum = relhum(:ncol,:), & ! unperturbed input
      rhu00 = rhu00(:ncol,:), & ! unperturbed input - output
      rhdfda = rhdfda(:ncol,:), & ! output for prognostic_cloud_water
      errmsg = errmsg, &
      errflg = errflg)

   call t_stopf("cldfrc")

   ! ---------------------------------------------- !
   ! Stratiform Cloud Macrophysics and Microphysics !
   ! ---------------------------------------------- !

   call t_startf('stratiform_microphys')

   !-------------------------------------------------------
   ! Non-micro and non-macrophysical external advective forcings to compute net condensation rate.
   ! Note that advective forcing of condensate is aggregated into liquid phase.
   !-------------------------------------------------------

!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   qtend(:,:) = 0._r8
   ttend(:,:) = 0._r8
   ltend(:,:) = 0._r8
!REMOVECAM_END

   call rk_stratiform_external_forcings_run( &
      ncol = ncol, &
      pver = pver, &
      dtime = dtime, &
      t = state1%t(:ncol,:), &
      q_wv = state1%q(:ncol,:pver,ixq), &
      cldice = state1%q(:ncol,:pver,ixcldice), & ! using before-repartitioning values for bit-for-bit.
      cldliq = state1%q(:ncol,:pver,ixcldliq), & ! using before-repartitioning values for bit-for-bit.
      qcwat = qcwat(:ncol,:), &
      tcwat = tcwat(:ncol,:), &
      lcwat = lcwat(:ncol,:), &
      qtend = qtend(:ncol,:), &
      ttend = ttend(:ncol,:), &
      ltend = ltend(:ncol,:), &
      errmsg = errmsg, &
      errflg = errflg)

   !-------------------------------------------------------
   ! Define fractional amount of stratus condensate and precipitation in ice phase.
   ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ).
   ! The ramp within convective cloud may be different
   !
   ! This ice fraction is used both in the repartitioning, and also in prognostic_cloud_water.
   !-------------------------------------------------------

!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   fice(:,:) = 0._r8
   fsnow(:,:) = 0._r8
!REMOVECAM_END
   top_lev = 1
   call cloud_fraction_fice_run(ncol, state1%t(:ncol,:), tmelt, top_lev, pver, fice(:ncol,:), fsnow(:ncol,:), errmsg, errflg)

   !-------------------------------------------------------
   ! Compute Stratiform Macro-Microphysical Tendencies
   !-------------------------------------------------------

   ! Add rain and snow fluxes as output variables from pcond, and into physics buffer
   call pbuf_get_field(pbuf, ls_flxprc_idx, rkflxprc)
   call pbuf_get_field(pbuf, ls_flxsnw_idx, rkflxsnw)

   call t_startf('pcond')

!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   qme(:,:) = 0._r8
   prain(:,:) = 0._r8
   prodsnow(:,:) = 0._r8
   nevapr(:,:) = 0._r8
   evapsnow(:,:) = 0._r8
   evapheat(:,:) = 0._r8
   prfzheat(:,:) = 0._r8
   meltheat(:,:) = 0._r8
   prec_pcw(:) = 0._r8
   snow_pcw(:) = 0._r8
   ice2pr(:,:) = 0._r8
   liq2pr(:,:) = 0._r8
   liq2snow(:,:) = 0._r8
   rkflxprc(:,:) = 0._r8
   rkflxsnw(:,:) = 0._r8
   pracwo(:,:) = 0._r8
   psacwo(:,:) = 0._r8
   psacio(:,:) = 0._r8
!REMOVECAM_END
   ! call CCPP-ized subroutine
   call prognostic_cloud_water_run( &
      ncol            = ncol,                  &
      pver            = pver,                  &
      top_lev         = trop_cloud_top_lev,    &
      deltat          = dtime,                 &
      iulog           = iulog,                 &
      pi              = pi,                    &
      gravit          = gravit,                &
      rh2o            = rh2o,                  &
      epsilo          = epsilo,                &
      latvap          = latvap,                &
      latice          = latice,                &
      cpair           = cpair,                 &
      dlat            = dlat(:ncol),           &
      pmid            = state1%pmid(:ncol,:),  &
      pdel            = state1%pdel(:ncol,:),  &
      zi              = state1%zi(:ncol,:),    &
      troplev         = troplev(:ncol),        &
      ttend           = ttend(:ncol,:),        &
      tn              = state1%t(:ncol,:),     &
      qtend           = qtend(:ncol,:),        &
      qn              = state1%q(:ncol,:,ixq), &
      ltend           = ltend(:ncol,:),        &
      cldice          = state1%q(:ncol,:pver,ixcldice), & ! for bit-for-bit, use values pre-repartitioning.
      cldliq          = state1%q(:ncol,:pver,ixcldliq), & ! for bit-for-bit, use values pre-repartitioning.
      !cwat            = totcw(:ncol,:),        &
      omega           = state1%omega(:ncol,:), &
      cldn            = cld(:ncol,:),          &
      fice            = fice(:ncol,:),         &
      fsnow           = fsnow(:ncol,:),        &
      rhdfda          = rhdfda(:ncol,:),       &
      rhu00           = rhu00(:ncol,:),        &
      landm           = landm(:ncol),          &
      seaicef         = icefrac(:ncol),        &
      snowh           = snowh(:ncol),          &
      ! below output
      qme             = qme(:ncol,:),          &
      prodprec        = prain(:ncol,:),        &
      prodsnow        = prodsnow(:ncol,:),     & ! ignored in rk.
      evapprec        = nevapr(:ncol,:),       &
      evapsnow        = evapsnow(:ncol,:),     &
      evapheat        = evapheat(:ncol,:),     &
      prfzheat        = prfzheat(:ncol,:),     &
      meltheat        = meltheat(:ncol,:),     &
      precip          = prec_pcw(:ncol),       &
      snowab          = snow_pcw(:ncol),       &
      ice2pr          = ice2pr(:ncol,:),       &
      liq2pr          = liq2pr(:ncol,:),       &
      liq2snow        = liq2snow(:ncol,:),     &
      lsflxprc        = rkflxprc(:ncol,:),     &
      lsflxsnw        = rkflxsnw(:ncol,:),     &
      pracwo          = pracwo(:ncol,:),       &
      psacwo          = psacwo(:ncol,:),       &
      psacio          = psacio(:ncol,:),       &
      fwaut           = fwaut(:ncol,:),        &
      fsaut           = fsaut(:ncol,:),        &
      fracw           = fracw(:ncol,:),        &
      fsacw           = fsacw(:ncol,:),        &
      fsaci           = fsaci(:ncol,:),        &
      errmsg          = errmsg,                &
      errflg          = errflg)
   if(errflg /= 0) then
      call endrun('rk_stratiform_tend:' // errmsg)
   endif

   call t_stopf('pcond')

   !-------------------------------------------------------
   ! Perform repartitioning of stratiform condensate.
   ! Corresponding heating tendency will be added later.
   !-------------------------------------------------------

   lq(:)        = .FALSE.
   lq(ixcldice) = .true.
   lq(ixcldliq) = .true.
   call physics_ptend_init( ptend_loc, state1%psetcols, 'cldwat-repartition', lq=lq )
!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   repartht(:,:) = 0._r8
!REMOVECAM_END

   call rk_stratiform_condensate_repartioning_run( &
      ncol = ncol, &
      pver = pver, &
      dtime = dtime, &
      latice = latice, &
      cldice = state1%q(:ncol,:,ixcldice), &
      cldliq = state1%q(:ncol,:,ixcldliq), &
      fice = fice(:ncol,:), & ! below output
      tend_cldice = ptend_loc%q(:ncol,:,ixcldice), &
      tend_cldliq = ptend_loc%q(:ncol,:,ixcldliq), &
      repartht = repartht(:ncol,:), &
      errmsg = errmsg, &
      errflg = errflg)

   call outfld( 'REPARTICE', ptend_loc%q(:,:,ixcldice), pcols, lchnk )
   call outfld( 'REPARTLIQ', ptend_loc%q(:,:,ixcldliq), pcols, lchnk )

   ! note: to restore old way of calculating repartht
   ! because apparently it causes errors greater than roundoff when accumulated
   ! in certain compsets of the CAM regression tests.
   ! in any case, I do not believe this is a good representation. It sets repartht
   ! to an intermediate value with the wrong units...
   ! repartht(:ncol,:pver)  = state1%q(:ncol,:pver,ixcldice)
   ! still bit differences:
   ! repartht(:ncol,:pver)  = (latice/dtime) * state1%q(:ncol,:pver,ixcldice)
   ! // note

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_update( state1, ptend_loc, dtime )

   ! note: second part of restoring old way of calculating repartht
   ! ...then does the difference here to determine repartition heating
   ! repartht(:ncol,:pver) = (latice/dtime) * ( state1%q(:ncol,:pver,ixcldice) - repartht(:ncol,:pver) )
   ! still bit differences if we switch to maintain proper units:
   ! repartht(:ncol,:pver) = (latice/dtime) * state1%q(:ncol,:pver,ixcldice) - repartht(:ncol,:pver)


   !-------------------------------------------------------
   ! APPLY Stratiform Macro-Microphysical Tendencies
   !-------------------------------------------------------

   lq(:)        = .FALSE.
   lq(ixq)        = .true.
   lq(ixcldice) = .true.
   lq(ixcldliq) = .true.
   call physics_ptend_init( ptend_loc, state1%psetcols, 'cldwat', ls=.true., lq=lq)

!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   cmeheat(:,:) = 0._r8
   cmeice(:,:) = 0._r8
   cmeliq(:,:) = 0._r8
!REMOVECAM_END
   ! call CCPP-ized subroutine to compute tendencies from prognostic_cloud_water
   ! (to heating rate, q, cldice, cldliq)
   call rk_stratiform_prognostic_cloud_water_tendencies_run( &
      ncol            = ncol,                  &
      pver            = pver,                  &
      dtime           = dtime,                 &
      latvap          = latvap,                &
      latice          = latice,                &
      qme             = qme(:ncol,:),          &
      fice            = fice(:ncol,:),         &
      evapheat        = evapheat(:ncol,:),     &
      prfzheat        = prfzheat(:ncol,:),     &
      meltheat        = meltheat(:ncol,:),     &
      repartht        = repartht(:ncol,:),     &
      evapprec        = nevapr(:ncol,:),       &
      ice2pr          = ice2pr(:ncol,:),       &
      liq2pr          = liq2pr(:ncol,:),       &
      cmeheat         = cmeheat(:ncol,:),      &
      cmeice          = cmeice(:ncol,:),       &
      cmeliq          = cmeliq(:ncol,:),       &
      prec_pcw        = prec_pcw(:ncol),       &
      snow_pcw        = snow_pcw(:ncol),       &
      prec_str        = prec_str(:ncol),       &
      snow_str        = snow_str(:ncol),       &
      tend_s          = ptend_loc%s(:ncol,:),  &
      tend_q          = ptend_loc%q(:ncol,:,ixq), &
      tend_cldice     = ptend_loc%q(:ncol,:,ixcldice), &
      tend_cldliq     = ptend_loc%q(:ncol,:,ixcldliq), &
      errmsg          = errmsg,                &
      errflg          = errflg)

   ! Record history variables
   call outfld( 'CME'     , qme,         pcols, lchnk )
   call outfld( 'PRODPREC' , prain,       pcols, lchnk )
   call outfld( 'EVAPPREC' , nevapr,      pcols, lchnk )
   call outfld( 'EVAPSNOW' , evapsnow,    pcols, lchnk )

   call outfld( 'FWAUT'   , fwaut,       pcols, lchnk )
   call outfld( 'FSAUT'   , fsaut,       pcols, lchnk )
   call outfld( 'FRACW'   , fracw,       pcols, lchnk )
   call outfld( 'FSACW'   , fsacw,       pcols, lchnk )
   call outfld( 'FSACI'   , fsaci,       pcols, lchnk )

   call outfld( 'PCSNOW'  , snow_pcw,    pcols, lchnk )
   call outfld( 'FICE'    , fice,        pcols, lchnk )
   call outfld( 'CMEICE'  , cmeice,      pcols, lchnk )
   call outfld( 'CMELIQ'  , cmeliq,      pcols, lchnk )
   call outfld( 'ICE2PR'  , ice2pr,      pcols, lchnk )
   call outfld( 'LIQ2PR'  , liq2pr,      pcols, lchnk )
   call outfld( 'HPROGCLD', ptend_loc%s, pcols, lchnk )
   call outfld( 'HEVAP   ', evapheat,    pcols, lchnk )
   call outfld( 'HMELT'   , meltheat,    pcols, lchnk )
   call outfld( 'HCME'    , cmeheat ,    pcols, lchnk )
   call outfld( 'HFREEZ'  , prfzheat,    pcols, lchnk )
   call outfld( 'HREPART' , repartht,    pcols, lchnk )
   call outfld('LS_FLXPRC', rkflxprc,    pcols, lchnk )
   call outfld('LS_FLXSNW', rkflxsnw,    pcols, lchnk )
   call outfld('PRACWO'   , pracwo,      pcols, lchnk )
   call outfld('PSACWO'   , psacwo,      pcols, lchnk )
   call outfld('PSACIO'   , psacio,      pcols, lchnk )

   ! initialize local variables for CFMIP diagnostics
   mr_ccliq(1:ncol,1:pver) = 0._r8
   mr_ccice(1:ncol,1:pver) = 0._r8
   mr_lsliq(1:ncol,1:pver) = 0._r8
   mr_lsice(1:ncol,1:pver) = 0._r8

   do k=1,pver
      do i=1,ncol
         if (cld(i,k) .gt. 0._r8) then
            mr_ccliq(i,k) = (state%q(i,k,ixcldliq)/cld(i,k))*concld(i,k)
            mr_ccice(i,k) = (state%q(i,k,ixcldice)/cld(i,k))*concld(i,k)
            mr_lsliq(i,k) = (state%q(i,k,ixcldliq)/cld(i,k))*(cld(i,k)-concld(i,k))
            mr_lsice(i,k) = (state%q(i,k,ixcldice)/cld(i,k))*(cld(i,k)-concld(i,k))
         else
            mr_ccliq(i,k) = 0._r8
            mr_ccice(i,k) = 0._r8
            mr_lsliq(i,k) = 0._r8
            mr_lsice(i,k) = 0._r8
         end if
      end do
   end do

   call outfld( 'CLDLIQSTR  ', mr_lsliq,    pcols, lchnk )
   call outfld( 'CLDICESTR  ', mr_lsice,    pcols, lchnk )
   call outfld( 'CLDLIQCON  ', mr_ccliq,    pcols, lchnk )
   call outfld( 'CLDICECON  ', mr_ccice,    pcols, lchnk )

   ! ------------------------------- !
   ! Update microphysical tendencies !
   ! ------------------------------- !

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_update( state1, ptend_loc, dtime )

   call t_startf("cldfrc")

   !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   cld(:,:) = 0._r8
   rhcloud(:,:) = 0._r8
   cldst(:,:) = 0._r8
   rhu00(:,:) = 0._r8
   icecldf(:,:) = 0._r8
   liqcldf(:,:) = 0._r8
   relhum(:,:) = 0._r8
   !REMOVECAM_END

   ! call CCPPized cloud fraction scheme
   call compute_cloud_fraction_run( &
      ncol = ncol, &
      pver = pver, &
      cappa = cappa, &
      gravit = gravit, &
      rair = rair, &
      tmelt = tmelt, &
      top_lev_cloudphys = 1, & ! CAM4 macrophysics.
      pmid = state1%pmid(:ncol,:), &
      ps = state1%pint(:,pverp), &
      temp = state1%t(:ncol,:), &
      sst = sst(:ncol), &
      q = state1%q(:ncol,:,ixq), &
      cldice = state1%q(:ncol,:,ixcldice), &
      phis = state1%phis(:ncol), &
      shallowcu = shallowcu(:ncol,:), &
      deepcu = deepcu(:ncol,:), &
      concld = concld(:ncol,:), &
      landfrac = landfrac(:ncol), &
      ocnfrac = ocnfrac(:ncol), &
      snowh = snowh(:ncol), &
      rhpert_flag = .false., & ! below output:
      cloud = cld(:ncol, :), &
      rhcloud = rhcloud(:ncol, :), &
      cldst = cldst(:ncol,:), &
      rhu00 = rhu00(:ncol,:), &
      icecldf = icecldf(:ncol,:), &
      liqcldf = liqcldf(:ncol,:), &
      relhum = relhum(:ncol,:), &
      errmsg = errmsg, &
      errflg = errflg)
   call t_stopf("cldfrc")

   call outfld( 'CLDST   ', cldst,  pcols, lchnk )

   do k = 1, pver
      do i = 1, ncol
         iwc(i,k)   = state1%q(i,k,ixcldice)*state1%pmid(i,k)/(rdair*state1%t(i,k))
         lwc(i,k)   = state1%q(i,k,ixcldliq)*state1%pmid(i,k)/(rdair*state1%t(i,k))
         icimr(i,k) = state1%q(i,k,ixcldice) / max(0.01_r8,rhcloud(i,k))
         icwmr(i,k) = state1%q(i,k,ixcldliq) / max(0.01_r8,rhcloud(i,k))
      end do
   end do

   call outfld( 'IWC'      , iwc,         pcols, lchnk )
   call outfld( 'LWC'      , lwc,         pcols, lchnk )
   call outfld( 'ICIMR'    , icimr,       pcols, lchnk )
   call outfld( 'ICWMR'    , icwmr,       pcols, lchnk )

   call t_stopf('stratiform_microphys')

   ! Save variables for use in the macrophysics at the next time step
   call rk_stratiform_save_qtlcwat_run( &
      ncol = ncol, &
      pver = pver, &
      t = state1%t(:ncol,:), &
      q_wv = state1%q(:ncol,:,ixq), &
      cldice = state1%q(:ncol,:,ixcldice), &
      cldliq = state1%q(:ncol,:,ixcldliq), &
      qcwat = qcwat(:ncol,:), &
      tcwat = tcwat(:ncol,:), &
      lcwat = lcwat(:ncol,:), &
      errmsg = errmsg, &
      errflg = errflg)

   ! Cloud water and ice particle sizes, saved in physics buffer for radiation
   !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   rel(:,:) = 0._r8
   rei(:,:) = 0._r8
   !REMOVECAM_END

   call rk_stratiform_cloud_optical_properties_run( &
      ncol = ncol, &
      pver = pver, &
      tmelt = tmelt, &
      landfrac = landfrac(:ncol), &
      icefrac = icefrac(:ncol), &
      snowh = snowh(:ncol), &
      landm = landm(:ncol), &
      t = state1%t(:ncol,:), &
      ps = state1%ps(:ncol), &
      pmid = state1%pmid(:ncol,:), & ! below output:
      rel = rel(:ncol,:), &
      rei = rei(:ncol,:), &
      errmsg = errmsg, &
      errflg = errflg)

   call outfld('REL', rel, pcols, lchnk)
   call outfld('REI', rei, pcols, lchnk)

   call physics_state_dealloc(state1)

end subroutine rk_stratiform_cam_tend

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_1(state1,ptend,i,k, &
        dtime,qme,fice,snow_pcw,prec_pcw, &
        prain,nevapr,prodsnow, evapsnow, &
        ice2pr,liq2pr,liq2snow)

     use physics_types, only: physics_state, physics_ptend
     use physconst,     only: tmelt

     implicit none

     integer, intent(in) :: i,k
     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     type(physics_ptend), intent(in) :: ptend  ! local copy of the ptend variable
     real(r8), intent(in)  :: dtime                ! timestep
     real(r8), intent(in) :: qme(pcols,pver)          ! local condensation - evaporation of cloud water

     real(r8), intent(in) :: prain(pcols,pver)          ! local production of precipitation
     real(r8), intent(in) :: nevapr(pcols,pver)          ! local evaporation of precipitation
     real(r8), intent(in) :: prodsnow(pcols,pver)          ! local production of snow
     real(r8), intent(in) :: evapsnow(pcols,pver)          ! local evaporation of snow
     real(r8), intent(in) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
     real(r8), intent(in) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
     real(r8), intent(in) :: liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
     real(r8), intent(in) :: fice    (pcols,pver)          ! Fractional ice content within cloud
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: prec_pcw(pcols)

     real(r8) hs1, qv1, ql1, qi1, qs1, qr1, fice2, pr1, w1, w2, w3, fliq, res
     real(r8) w4, wl, wv, wi, wlf, wvf, wif, qif, qlf, qvf

     pr1 = 0
     hs1 = 0
     qv1 = 0
     ql1 = 0
     qi1 = 0
     qs1 = 0
     qr1 = 0
     w1 = 0
     wl = 0
     wv = 0
     wi = 0
     wlf = 0
     wvf = 0
     wif = 0


     write(iulog,*)
     write(iulog,*) ' input state, t, q, l, i ', k, state1%t(i,k), state1%q(i,k,1), state1%q(i,k,ixcldliq),  state1%q(i,k,ixcldice)
     write(iulog,*) ' rain, snow, total from components before accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip before accumulation                      ', k, pr1

     wv = wv + state1%q(i,k,1       )*state1%pdel(i,k)/gravit
     wl = wl + state1%q(i,k,ixcldliq)*state1%pdel(i,k)/gravit
     wi = wi + state1%q(i,k,ixcldice)*state1%pdel(i,k)/gravit

     qvf = state1%q(i,k,1) + ptend%q(i,k,1)*dtime
     qlf = state1%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)*dtime
     qif = state1%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dtime

     if (qvf.lt.0._r8) then
        write(iulog,*) ' qvf is negative *******', qvf
     endif
     if (qlf.lt.0._r8) then
        write(iulog,*) ' qlf is negative *******', qlf
     endif
     if (qif.lt.0._r8) then
        write(iulog,*) ' qif is negative *******', qif
     endif
     write(iulog,*) ' qvf, qlf, qif ', qvf, qlf, qif

     wvf = wvf + qvf*state1%pdel(i,k)/gravit
     wlf = wlf + qlf*state1%pdel(i,k)/gravit
     wif = wif + qif*state1%pdel(i,k)/gravit

     hs1 = hs1 + ptend%s(i,k)*state1%pdel(i,k)/gravit
     pr1 = pr1 + state1%pdel(i,k)/gravit*(prain(i,k)-nevapr(i,k))
     qv1 = qv1 - (qme(i,k)-nevapr(i,k))*state1%pdel(i,k)/gravit    ! vdot
     w1  = w1  + (qme(i,k)-prain(i,k))*state1%pdel(i,k)/gravit    ! cdot
     qi1 = qi1 + ((qme(i,k))*fice(i,k)        -ice2pr(i,k) )*state1%pdel(i,k)/gravit   ! idot
     ql1 = ql1 + ((qme(i,k))*(1._r8-fice(i,k))-liq2pr(i,k) )*state1%pdel(i,k)/gravit   ! ldot

     qr1 = qr1 &
          + ( liq2pr(i,k)-liq2snow(i,k)   &     ! production of rain
          -(nevapr(i,k)-evapsnow(i,k)) &     ! rain evaporation
          )*state1%pdel(i,k)/gravit
     qs1 = qs1 &
          + ( ice2pr(i,k) + liq2snow(i,k) &     ! production of snow.Note last term has phase change
          -evapsnow(i,k)               &     ! snow evaporation
          )*state1%pdel(i,k)/gravit

     if (state1%t(i,k).gt.tmelt) then
        qr1 = qr1 + qs1
        qs1 = 0._r8
     endif
     write(iulog,*) ' rain, snow, total after accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip after accumulation      ', k, pr1
     write(iulog,*)
     write(iulog,*) ' layer prain, nevapr, pdel ', prain(i,k), nevapr(i,k), state1%pdel(i,k)
     write(iulog,*) ' layer prodsnow, ice2pr+liq2snow ', prodsnow(i,k), ice2pr(i,k)+liq2snow(i,k)
     write(iulog,*) ' layer prain-prodsnow, liq2pr-liq2snow ', prain(i,k)-prodsnow(i,k), liq2pr(i,k)-liq2snow(i,k)
     write(iulog,*) ' layer evapsnow, evaprain ', k, evapsnow(i,k), nevapr(i,k)-evapsnow(i,k)
     write(iulog,*) ' layer ice2pr, liq2pr, liq2snow ', ice2pr(i,k), liq2pr(i,k), liq2snow(i,k)
     write(iulog,*) ' layer ice2pr+liq2pr, prain ', ice2pr(i,k)+liq2pr(i,k), prain(i,k)
     write(iulog,*)
     write(iulog,*) ' qv1 vapor removed from col after accum  (vdot)   ', k, qv1
     write(iulog,*) ' - (precip produced - vapor removed) after accum  ', k, -pr1-qv1
     write(iulog,*) ' condensate produce after accum                   ', k, w1
     write(iulog,*) ' liq+ice tends accum                              ', k, ql1+qi1
     write(iulog,*) ' change in total water after accum                ', k, qv1+ql1+qi1
     write(iulog,*) ' imbalance in colum after accum                   ', k, qs1+qr1+qv1+ql1+qi1
     write(iulog,*) ' fice at this lev ', fice(i,k)
     write(iulog,*)

     res = abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1),abs(ql1),abs(qi1),abs(qs1),abs(qr1),1.e-36_r8))
     write(iulog,*) ' relative residual in column method 1             ', k, res

     write(iulog,*) ' relative residual in column method 2             ',&
	 k, abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1+ql1+qi1),1.e-36_r8))
     !            if (abs((qs1+qr1+qv1+ql1+qi1)/(qs1+qr1+1.e-36)).gt.1.e-14) then
     if (res.gt.1.e-14_r8) then
        call endrun ('STRATIFORM_TEND')
     endif

     !             w3  = qme(i,k) * (latvap + latice*fice(i,k)) &
     !               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k)

     res = qs1+qr1-pr1
     w4 = max(abs(qs1),abs(qr1),abs(pr1))
     if (w4.gt.0._r8)  then
        if (res/w4.gt.1.e-14_r8) then
           write(iulog,*) ' imbalance in precips calculated two ways '
           write(iulog,*) ' res/w4, pr1, qr1, qs1, qr1+qs1 ', &
                res/w4, pr1, qr1, qs1, qr1+qs1
           !                   call endrun()
        endif
     endif
     if (k.eq.pver) then
        write(iulog,*) ' pcond returned precip, rain and snow rates ', prec_pcw(i), prec_pcw(i)-snow_pcw(i), snow_pcw(i)
        write(iulog,*) ' I calculate ', pr1, qr1, qs1
        !               call endrun
        write(iulog,*) ' byrons water check ', wv+wl+wi-pr1*dtime, wvf+wlf+wif
     endif
     write(iulog,*)


   end subroutine debug_microphys_1

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_2(state1,&
        snow_pcw,fsaut,fsacw ,fsaci, meltheat)

     use ppgrid,        only: pver
     use physconst,     only: tmelt
     use physics_types, only: physics_state

     implicit none

     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: fsaut(pcols,pver)
     real(r8), intent(in) :: fsacw(pcols,pver)
     real(r8), intent(in) :: fsaci(pcols,pver)
     real(r8), intent(in) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip


     integer  i,ncol,lchnk


     ncol = state1%ncol
     lchnk = state1%lchnk

     do i = 1,ncol
        if (snow_pcw(i) .gt. 0.01_r8/8.64e4_r8  .and.  state1%t(i,pver) .gt. tmelt) then
           write(iulog,*) ' stratiform: snow, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif

        if (snow_pcw(i)*8.64e4_r8 .lt. -1.e-5_r8) then
           write(iulog,*) ' neg snow ', snow_pcw(i)*8.64e4_r8
           write(iulog,*) ' stratiform: snow_pcw, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
     end do

   end subroutine debug_microphys_2

  end module rk_stratiform_cam
