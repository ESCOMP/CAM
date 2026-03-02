! CAM interface for UW shallow convection.
! The actual computation has been converted to CCPP-compliant physics for CAM-SIMA
! and is in src/atmos_phys.
  module uwshcu

  use cam_logfile,    only: iulog
  use ppgrid,         only: pcols, pver, pverp
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc

  implicit none
  private
  save

  public &
     uwshcu_readnl,      &
     init_uwshcu,        &
     uwshcu_cam
  
  integer , parameter :: r8 = selected_real_kind(12)    !  8 byte real

  ! Tuning parameters set via namelist
  real(r8) :: rpen          !  For penetrative entrainment efficiency

!===============================================================================
contains
!===============================================================================
  
subroutine uwshcu_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'uwshcu_readnl'

   ! Namelist variables
   real(r8), parameter :: unset_r8 = huge(1.0_r8)
   real(r8) :: uwshcu_rpen =  unset_r8    !  For penetrative entrainment efficiency

   namelist /uwshcu_nl/ uwshcu_rpen
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'uwshcu_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, uwshcu_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(uwshcu_rpen,            1, mpir8,  0, mpicom)
#endif
   
   rpen=uwshcu_rpen
  

end subroutine uwshcu_readnl

!===============================================================================

  subroutine init_uwshcu( kind, xlv_in, cp_in, xlf_in, zvir_in, r_in, g_in, ep2_in )

    !------------------------------------------------------------- ! 
    ! Purpose:                                                     !
    ! Initialize key constants for the shallow convection package. !
    !------------------------------------------------------------- !

    use cam_history,   only: addfld, horiz_only

    use uw_convect_shallow, only: uw_convect_shallow_init

    integer , intent(in) :: kind       !  kind of reals being passed in
    real(r8), intent(in) :: xlv_in     !  Latent heat of vaporization
    real(r8), intent(in) :: xlf_in     !  Latent heat of fusion
    real(r8), intent(in) :: cp_in      !  Specific heat of dry air
    real(r8), intent(in) :: zvir_in    !  rh2o/rair - 1
    real(r8), intent(in) :: r_in       !  Gas constant for dry air
    real(r8), intent(in) :: g_in       !  Gravitational constant
    real(r8), intent(in) :: ep2_in     !  mol wgt water vapor / mol wgt dry air 

    character(len=*), parameter :: subname = 'init_uwshcu'

    character(len=512)   :: errmsg
    integer              :: errflg

    ! ------------------------- !
    ! Internal Output Variables !
    ! ------------------------- !

    call addfld( 'qtflx_Cu'       , (/ 'ilev' /),  'A', 'kg/m2/s' , 'Convective qt flux'         )
    call addfld( 'slflx_Cu'       , (/ 'ilev' /),  'A', 'J/m2/s'  , 'Convective sl flux'         )
    call addfld( 'uflx_Cu'        , (/ 'ilev' /),  'A', 'kg/m/s2' , 'Convective  u flux'         )
    call addfld( 'vflx_Cu'        , (/ 'ilev' /),  'A', 'kg/m/s2' , 'Convective  v flux'         )

    call addfld( 'qtten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'qt tendency by convection'  )
    call addfld( 'slten_Cu'       , (/ 'lev' /),   'A', 'J/kg/s'  , 'sl tendency by convection'  )
    call addfld( 'uten_Cu'        , (/ 'lev' /),   'A', 'm/s2'    , ' u tendency by convection'  )
    call addfld( 'vten_Cu'        , (/ 'lev' /),   'A', 'm/s2'    , ' v tendency by convection'  )
    call addfld( 'qvten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'qv tendency by convection'  )
    call addfld( 'qlten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'ql tendency by convection'  )
    call addfld( 'qiten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'qi tendency by convection'  )

    call addfld( 'cbmf_Cu'        , horiz_only,    'A', 'kg/m2/s' , 'Cumulus base mass flux'                               )
    call addfld( 'ufrcinvbase_Cu' , horiz_only,    'A', 'fraction', 'Cumulus fraction at PBL top'                          )
    call addfld( 'ufrclcl_Cu'     , horiz_only,    'A', 'fraction', 'Cumulus fraction at LCL'                              )
    call addfld( 'winvbase_Cu'    , horiz_only,    'A', 'm/s'     , 'Cumulus vertical velocity at PBL top'                 )
    call addfld( 'wlcl_Cu'        , horiz_only,    'A', 'm/s'     , 'Cumulus vertical velocity at LCL'                     )
    call addfld( 'plcl_Cu'        , horiz_only,    'A', 'Pa'      , 'LCL of source air'                                    )
    call addfld( 'pinv_Cu'        , horiz_only,    'A', 'Pa'      , 'PBL top pressure'                                     )
    call addfld( 'plfc_Cu'        , horiz_only,    'A', 'Pa'      , 'LFC of source air'                                    )
    call addfld( 'pbup_Cu'        , horiz_only,    'A', 'Pa'      , 'Highest interface level of positive cumulus buoyancy' )
    call addfld( 'ppen_Cu'        , horiz_only,    'A', 'Pa'      , 'Highest level where cumulus w is 0'                   )
    call addfld( 'qtsrc_Cu'       , horiz_only,    'A', 'kg/kg'   , 'Cumulus source air qt'                                )
    call addfld( 'thlsrc_Cu'      , horiz_only,    'A', 'K'       , 'Cumulus source air thl'                               )
    call addfld( 'thvlsrc_Cu'     , horiz_only,    'A', 'K'       , 'Cumulus source air thvl'                              )
    call addfld( 'emfkbup_Cu'     , horiz_only,    'A', 'kg/m2/s' ,  'Penetrative mass flux at kbup'                       )
    call addfld( 'cin_Cu'         , horiz_only,    'A', 'J/kg'    , 'CIN upto LFC'                                         )
    call addfld( 'cinlcl_Cu'      , horiz_only,    'A', 'J/kg'    , 'CIN upto LCL'                                         )
    call addfld( 'cbmflimit_Cu'   , horiz_only,    'A', 'kg/m2/s' , 'cbmflimiter'                                          )
    call addfld( 'tkeavg_Cu'      , horiz_only,    'A', 'm2/s2'   , 'Average tke within PBL for convection scheme'         )
    call addfld( 'zinv_Cu'        , horiz_only,    'A', 'm'       , 'PBL top height'                                       )
    call addfld( 'rcwp_Cu'        , horiz_only,    'A', 'kg/m2'   , 'Cumulus LWP+IWP'                                      )
    call addfld( 'rlwp_Cu'        , horiz_only,    'A', 'kg/m2'   , 'Cumulus LWP'                                          )
    call addfld( 'riwp_Cu'        , horiz_only,    'A', 'kg/m2'   , 'Cumulus IWP'                                          )
    call addfld( 'tophgt_Cu'      , horiz_only,    'A', 'm'       , 'Cumulus top height'                                   )

    call addfld( 'wu_Cu'          , (/ 'ilev' /),  'A', 'm/s'     , 'Convective updraft vertical velocity'             )
    call addfld( 'ufrc_Cu'        , (/ 'ilev' /),  'A', 'fraction', 'Convective updraft fractional area'               )
    call addfld( 'qtu_Cu'         , (/ 'ilev' /),  'A', 'kg/kg'   , 'Cumulus updraft qt'                               )
    call addfld( 'thlu_Cu'        , (/ 'ilev' /),  'A', 'K'       , 'Cumulus updraft thl'                              )
    call addfld( 'thvu_Cu'        , (/ 'ilev' /),  'A', 'K'       , 'Cumulus updraft thv'                              )
    call addfld( 'uu_Cu'          , (/ 'ilev' /),  'A', 'm/s'     , 'Cumulus updraft uwnd'                             )
    call addfld( 'vu_Cu'          , (/ 'ilev' /),  'A', 'm/s'     , 'Cumulus updraft vwnd'                             )
    call addfld( 'qtu_emf_Cu'     , (/ 'ilev' /),  'A', 'kg/kg'   , 'qt of penatratively entrained air'                )
    call addfld( 'thlu_emf_Cu'    , (/ 'ilev' /),  'A', 'K'       , 'thl of penatratively entrained air'               )
    call addfld( 'uu_emf_Cu'      , (/ 'ilev' /),  'A', 'm/s'     , 'uwnd of penatratively entrained air'              )
    call addfld( 'vu_emf_Cu'      , (/ 'ilev' /),  'A', 'm/s'     , 'vwnd of penatratively entrained air'              )
    call addfld( 'umf_Cu'         , (/ 'ilev' /),  'A', 'kg/m2/s' , 'Cumulus updraft mass flux'                        )
    call addfld( 'uemf_Cu'        , (/ 'ilev' /),  'A', 'kg/m2/s' , 'Cumulus net ( updraft + entrainment ) mass flux'  )
    call addfld( 'qcu_Cu'         , (/ 'lev' /),   'A', 'kg/kg'   , 'Cumulus updraft LWC+IWC'                          )
    call addfld( 'qlu_Cu'         , (/ 'lev' /),   'A', 'kg/kg'   , 'Cumulus updraft LWC'                              )
    call addfld( 'qiu_Cu'         , (/ 'lev' /),   'A', 'kg/kg'   , 'Cumulus updraft IWC'                              )
    call addfld( 'cufrc_Cu'       , (/ 'lev' /),   'A', 'fraction', 'Cumulus cloud fraction'                           )
    call addfld( 'fer_Cu'         , (/ 'lev' /),   'A', '1/m'     , 'Cumulus lateral fractional entrainment rate'      )
    call addfld( 'fdr_Cu'         , (/ 'lev' /),   'A', '1/m'     , 'Cumulus lateral fractional detrainment Rate'      )

    call addfld( 'dwten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'Expellsion rate of cumulus cloud water to env.'   )
    call addfld( 'diten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'Expellsion rate of cumulus ice water to env.'     )
    call addfld( 'qrten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'Production rate of rain by cumulus'               )
    call addfld( 'qsten_Cu'       , (/ 'lev' /),   'A', 'kg/kg/s' , 'Production rate of snow by cumulus'               )
    call addfld( 'flxrain_Cu'     , (/ 'ilev' /),  'A', 'kg/m2/s' , 'Rain flux induced by Cumulus'                     )
    call addfld( 'flxsnow_Cu'     , (/ 'ilev' /),  'A', 'kg/m2/s' , 'Snow flux induced by Cumulus'                     )
    call addfld( 'ntraprd_Cu'     , (/ 'lev' /),   'A', 'kg/kg/s' , 'Net production rate of rain by Cumulus'           )
    call addfld( 'ntsnprd_Cu'     , (/ 'lev' /),   'A', 'kg/kg/s' , 'Net production rate of snow by Cumulus'           )

    call addfld( 'excessu_Cu'     , (/ 'lev' /),   'A', 'no'      , 'Updraft saturation excess'                        )
    call addfld( 'excess0_Cu'     , (/ 'lev' /),   'A', 'no'      , 'Environmental saturation excess'                  )
    call addfld( 'xc_Cu'          , (/ 'lev' /),   'A', 'no'      , 'Critical ncoling ratio'                            )
    call addfld( 'aquad_Cu'       , (/ 'lev' /),   'A', 'no'      , 'aquad'                                            )
    call addfld( 'bquad_Cu'       , (/ 'lev' /),   'A', 'no'      , 'bquad'                                            )
    call addfld( 'cquad_Cu'       , (/ 'lev' /),   'A', 'no'      , 'cquad'                                            )
    call addfld( 'bogbot_Cu'      , (/ 'lev' /),   'A', 'no'      , 'Cloud buoyancy at the bottom interface'           )
    call addfld( 'bogtop_Cu'      , (/ 'lev' /),   'A', 'no'      , 'Cloud buoyancy at the top interface'              )

    call addfld('exit_UWCu_Cu'    , horiz_only,    'A', 'no' , 'exit_UWCu'     )
    call addfld('exit_conden_Cu'  , horiz_only,    'A', 'no' , 'exit_conden'   )
    call addfld('exit_klclpver_Cu' , horiz_only,    'A', 'no' , 'exit_klclpver'  )
    call addfld('exit_klfcpver_Cu' , horiz_only,    'A', 'no' , 'exit_klfcpver'  )
    call addfld('exit_ufrc_Cu'    , horiz_only,    'A', 'no' , 'exit_ufrc'     )
    call addfld('exit_wtw_Cu'     , horiz_only,    'A', 'no' , 'exit_wtw'      )
    call addfld('exit_drycore_Cu' , horiz_only,    'A', 'no' , 'exit_drycore'  )
    call addfld('exit_wu_Cu'      , horiz_only,    'A', 'no' , 'exit_wu'       )
    call addfld('exit_cufilter_Cu', horiz_only,    'A', 'no' , 'exit_cufilter' )
    call addfld('exit_kinv1_Cu'   , horiz_only,    'A', 'no' , 'exit_kinv1'    )
    call addfld('exit_rei_Cu'     , horiz_only,    'A', 'no' , 'exit_rei'      )

    call addfld('limit_shcu_Cu'   , horiz_only,    'A', 'no' , 'limit_shcu'    )
    call addfld('limit_negcon_Cu' , horiz_only,    'A', 'no' , 'limit_negcon'  )
    call addfld('limit_ufrc_Cu'   , horiz_only,    'A', 'no' , 'limit_ufrc'    )
    call addfld('limit_ppen_Cu'   , horiz_only,    'A', 'no' , 'limit_ppen'    )
    call addfld('limit_emf_Cu'    , horiz_only,    'A', 'no' , 'limit_emf'     )
    call addfld('limit_cinlcl_Cu' , horiz_only,    'A', 'no' , 'limit_cinlcl'  )
    call addfld('limit_cin_Cu'    , horiz_only,    'A', 'no' , 'limit_cin'     )
    call addfld('limit_cbmf_Cu'   , horiz_only,    'A', 'no' , 'limit_cbmf'    )
    call addfld('limit_rei_Cu'    , horiz_only,    'A', 'no' , 'limit_rei'     )
    call addfld('ind_delcin_Cu'   , horiz_only,    'A', 'no' , 'ind_delcin'    )

    if( kind .ne. r8 ) then
        write(iulog,*) subname//': ERROR -- real KIND does not match internal specification.'
        call endrun(subname//': ERROR -- real KIND does not match internal specification.')
    endif

    ! call the underlying CCPPized subroutine
    call uw_convect_shallow_init(masterproc, iulog, rpen, xlv_in, cp_in, xlf_in, zvir_in, r_in, g_in, ep2_in, errmsg, errflg)

    if(errflg /= 0) then
      call endrun(subname//': '//errmsg)
    end if

  end subroutine init_uwshcu

  subroutine uwshcu_cam(pcols, ncol      , pver      , ncnst     , dt       ,  &
                        ps0_inv  , zs0_inv    , p0_inv        , z0_inv    , dp0_inv  ,  &
                        u0_inv   , v0_inv     , qv0_inv       , ql0_inv   , qi0_inv  ,  &
                        t0_inv   , s0_inv     , tr0_inv       ,                         &
                        tke_inv  , cldfrct_inv, concldfrct_inv, pblh      , cush     ,  &
                        umf_inv  , slflx_inv  , qtflx_inv     ,                         &
                        flxprc1_inv, flxsnow1_inv,                 &
                        sten_inv , uten_inv   , vten_inv      , trten_inv ,             &
                        qrten_inv, qsten_inv  , precip        , snow      , evapc_inv,  &
                        cufrc_inv, qcu_inv    , qlu_inv       , qiu_inv   ,             &
                        cbmf     , qc_inv     , rliq          ,                         &
                        cnt_inv  , cnb_inv    , lchnk         , dpdry0_inv,             &
                        sh_e_ed_ratio)
    use cam_history,     only : outfld
    use ccpp_constituent_prop_mod, only: ccpp_const_props

    use uw_convect_shallow, only: uw_convect_shallow_run

    integer , intent(in)    :: pcols
    integer , intent(in)    :: lchnk
    integer , intent(in)    :: ncol
    integer , intent(in)    :: pver
    integer , intent(in)    :: ncnst
    real(r8), intent(in)    :: dt                       !  Time step : 2*delta_t [ s ]
    real(r8), intent(in)    :: ps0_inv(pcols,pver+1)       !  Environmental pressure at the interfaces [ Pa ]
    real(r8), intent(in)    :: zs0_inv(pcols,pver+1)       !  Environmental height at the interfaces   [ m ]
    real(r8), intent(in)    :: p0_inv(pcols,pver)          !  Environmental pressure at the layer mid-point [ Pa ]
    real(r8), intent(in)    :: z0_inv(pcols,pver)          !  Environmental height at the layer mid-point [ m ]
    real(r8), intent(in)    :: dp0_inv(pcols,pver)         !  Environmental layer pressure thickness [ Pa ] > 0.
    real(r8), intent(in)    :: dpdry0_inv(pcols,pver)      !  Environmental dry layer pressure thickness [ Pa ]
    real(r8), intent(in)    :: u0_inv(pcols,pver)          !  Environmental zonal wind [ m/s ]
    real(r8), intent(in)    :: v0_inv(pcols,pver)          !  Environmental meridional wind [ m/s ]
    real(r8), intent(in)    :: qv0_inv(pcols,pver)         !  Environmental water vapor specific humidity [ kg/kg ]
    real(r8), intent(in)    :: ql0_inv(pcols,pver)         !  Environmental liquid water specific humidity [ kg/kg ]
    real(r8), intent(in)    :: qi0_inv(pcols,pver)         !  Environmental ice specific humidity [ kg/kg ]
    real(r8), intent(in)    :: t0_inv(pcols,pver)          !  Environmental temperature [ K ]
    real(r8), intent(in)    :: s0_inv(pcols,pver)          !  Environmental dry static energy [ J/kg ]
    real(r8), intent(in)    :: tr0_inv(pcols,pver,ncnst)   !  Environmental tracers [ #, kg/kg ]
    real(r8), intent(in)    :: tke_inv(pcols,pver+1)       !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    real(r8), intent(in)    :: cldfrct_inv(pcols,pver)     !  Total cloud fraction at the previous time step [ fraction ]
    real(r8), intent(in)    :: concldfrct_inv(pcols,pver)  !  Total convective ( shallow + deep ) cloud fraction
                                                        !  at the previous time step [ fraction ]
    real(r8), intent(in)    :: pblh(pcols)                !  Height of PBL [ m ]
    real(r8), intent(inout) :: cush(pcols)                !  Convective scale height [ m ]
    real(r8), intent(out)   :: umf_inv(pcols,pver+1)       !  Updraft mass flux at the interfaces [ kg/m2/s ]
    real(r8), intent(out)   :: sten_inv(pcols,pver)        !  Tendency of dry static energy [ J/kg/s ]
    real(r8), intent(out)   :: uten_inv(pcols,pver)        !  Tendency of zonal wind [ m/s2 ]
    real(r8), intent(out)   :: vten_inv(pcols,pver)        !  Tendency of meridional wind [ m/s2 ]
    real(r8), intent(out)   :: trten_inv(pcols,pver,ncnst) !  Tendency of tracers [ #/s, kg/kg/s ]
    real(r8), intent(out)   :: qrten_inv(pcols,pver)       !  Tendency of rain water specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: qsten_inv(pcols,pver)       !  Tendency of snow specific humidity [ kg/kg/s ]
    real(r8), intent(out)   :: precip(pcols)              !  Precipitation ( rain + snow ) flux at the surface [ m/s ]
    real(r8), intent(out)   :: snow(pcols)                !  Snow flux at the surface [ m/s ]
    real(r8), intent(out)   :: evapc_inv(pcols,pver)       !  Evaporation of precipitation [ kg/kg/s ]
    real(r8), intent(out)   :: rliq(pcols)                !  Vertical integral of tendency of detrained cloud condensate qc [ m/s ]
    real(r8), intent(out)   :: slflx_inv(pcols,pver+1)     !  Updraft liquid static energy flux [ J/kg * kg/m2/s ]
    real(r8), intent(out)   :: qtflx_inv(pcols,pver+1)     !  Updraft total water flux [ kg/kg * kg/m2/s ]
    real(r8), intent(out)   :: flxprc1_inv(pcols,pver+1)   ! uw grid-box mean rain+snow flux (kg m^-2 s^-1)
                                                        ! for physics buffer calls in convect_shallow.F90
    real(r8), intent(out)   :: flxsnow1_inv(pcols,pver+1)  ! uw grid-box mean snow flux (kg m^-2 s^-1)
                                                        ! for physics buffer calls in convect_shallow.F90

    real(r8), intent(out)   :: cufrc_inv(pcols,pver)       !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    real(r8), intent(out)   :: qcu_inv(pcols,pver)         !  Liquid+ice specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qlu_inv(pcols,pver)         !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qiu_inv(pcols,pver)         !  Ice specific humidity within cumulus updraft [ kg/kg ]
    real(r8), intent(out)   :: qc_inv(pcols,pver)          !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
    real(r8), intent(out)   :: cbmf(pcols)                !  Cumulus base mass flux [ kg/m2/s ]
    real(r8), intent(out)   :: cnt_inv(pcols)             !  Cumulus top  interface index, cnt = kpen [ no ]
    real(r8), intent(out)   :: cnb_inv(pcols)             !  Cumulus base interface index, cnb = krel - 1 [ no ]

    real(r8), intent(out)   :: sh_e_ed_ratio(pcols,pver)   !  shallow conv [ent/(ent+det)] ratio

    character(len=512)   :: errmsg
    integer              :: errflg

    ! zero pcols before subsetting
    umf_inv(:, :)       = 0.0_r8
    slflx_inv(:, :)     = 0.0_r8
    qtflx_inv(:, :)     = 0.0_r8
    flxprc1_inv(:, :)   = 0.0_r8
    flxsnow1_inv(:, :)  = 0.0_r8
    sten_inv(:, :)      = 0.0_r8
    uten_inv(:, :)      = 0.0_r8
    vten_inv(:, :)      = 0.0_r8
    trten_inv(:, :, :)  = 0.0_r8
    qrten_inv(:, :)     = 0.0_r8
    qsten_inv(:, :)     = 0.0_r8
    precip(:)           = 0.0_r8
    snow(:)             = 0.0_r8
    evapc_inv(:, :)     = 0.0_r8
    cufrc_inv(:, :)     = 0.0_r8
    qcu_inv(:, :)       = 0.0_r8
    qlu_inv(:, :)       = 0.0_r8
    qiu_inv(:, :)       = 0.0_r8
    qc_inv(:, :)        = 0.0_r8
    cbmf(:)             = 0.0_r8
    rliq(:)             = 0.0_r8
    cnt_inv(:)          = 0.0_r8
    cnb_inv(:)          = 0.0_r8
    sh_e_ed_ratio(:, :) = 0.0_r8

    ! call the underlying CCPPized subroutine (dechunkized)
    call uw_convect_shallow_run( &
      ncol          = ncol,                             &
      pver          = pver,                             &
      ncnst         = ncnst,                            &
      dt            = dt,                               &
      const_props   = ccpp_const_props,                 &
      pint          = ps0_inv(:ncol, :pver+1),        & ! interfaces (ncol, pver+1)
      zi            = zs0_inv(:ncol, :pver+1),        & ! interfaces (ncol, pver+1)
      pmid          = p0_inv(:ncol, :pver),            & ! midpoints  (ncol, pver)
      zm            = z0_inv(:ncol, :pver),            & ! midpoints  (ncol, pver)
      pdel          = dp0_inv(:ncol, :pver),           & ! midpoints  (ncol, pver)
      pdeldry       = dpdry0_inv(:ncol, :pver),        & ! midpoints  (ncol, pver)
      u             = u0_inv(:ncol, :pver),            & ! midpoints  (ncol, pver)
      v             = v0_inv(:ncol, :pver),            & ! midpoints  (ncol, pver)
      qv0           = qv0_inv(:ncol, :pver),           & ! midpoints  (ncol, pver)
      ql0           = ql0_inv(:ncol, :pver),           & ! midpoints  (ncol, pver)
      qi0           = qi0_inv(:ncol, :pver),           & ! midpoints  (ncol, pver)
      t             = t0_inv(:ncol, :pver),            & ! midpoints  (ncol, pver)
      s             = s0_inv(:ncol, :pver),            & ! midpoints  (ncol, pver)
      tr0           = tr0_inv(:ncol, :pver, 1:ncnst),  & ! midpoints  (ncol, pver, ncnst)
      tke           = tke_inv(:ncol, :pver+1),         & ! interfaces (ncol, pver+1)
      cld           = cldfrct_inv(:ncol, :pver),       & ! midpoints  (ncol, pver)
      concld        = concldfrct_inv(:ncol, :pver),    & ! midpoints  (ncol, pver)
      pblh          = pblh(:ncol),                      &
      cush          = cush(:ncol),                      &
      cmfmc_sh      = umf_inv(:ncol, :pver+1),        & ! interfaces (ncol, pver+1)
      slflx         = slflx_inv(:ncol, :pver+1),      & ! interfaces (ncol, pver+1)
      qtflx         = qtflx_inv(:ncol, :pver+1),      & ! interfaces (ncol, pver+1)
      flxprc_sh     = flxprc1_inv(:ncol, :pver+1),    & ! interfaces (ncol, pver+1)
      flxsnw_sh     = flxsnow1_inv(:ncol, :pver+1),   & ! interfaces (ncol, pver+1)
      sten          = sten_inv(:ncol, :pver),          & ! midpoints  (ncol, pver)
      uten          = uten_inv(:ncol, :pver),          & ! midpoints  (ncol, pver)
      vten          = vten_inv(:ncol, :pver),          & ! midpoints  (ncol, pver)
      trten         = trten_inv(:ncol, :pver, 1:ncnst),& ! midpoints  (ncol, pver, ncnst)
      qrten         = qrten_inv(:ncol, :pver),         & ! midpoints  (ncol, pver)
      qsten         = qsten_inv(:ncol, :pver),         & ! midpoints  (ncol, pver)
      precip_sh     = precip(:ncol),                    &
      snow_sh       = snow(:ncol),                      &
      evapc_sh      = evapc_inv(:ncol, :pver),        & ! midpoints  (ncol, pver)
      shfrc         = cufrc_inv(:ncol, :pver),         & ! midpoints  (ncol, pver)
      qcu           = qcu_inv(:ncol, :pver),           & ! midpoints  (ncol, pver)
      qlu           = qlu_inv(:ncol, :pver),           & ! midpoints  (ncol, pver)
      qiu           = qiu_inv(:ncol, :pver),           & ! midpoints  (ncol, pver)
      qc            = qc_inv(:ncol, :pver),            & ! midpoints  (ncol, pver)
      cbmf          = cbmf(:ncol),                      &
      rliq          = rliq(:ncol),                      &
      cnt           = cnt_inv(:ncol),                   &
      cnb           = cnb_inv(:ncol),                   &
      sh_e_ed_ratio = sh_e_ed_ratio(:ncol, :pver),    & ! midpoints  (ncol, pver)
      errmsg        = errmsg,                            &
      errflg        = errflg)

     ! ---------------------------------------- !
     ! Writing main diagnostic output variables !
     ! ---------------------------------------- !

     ! call outfld( 'qtflx_Cu'        , qtflx_out(:,pver:0:-1),    ncol,    lchnk )
     ! call outfld( 'slflx_Cu'        , slflx_out(:,pver:0:-1),    ncol,    lchnk )
     ! call outfld( 'uflx_Cu'         , uflx_out,                 ncol,    lchnk )
     ! call outfld( 'vflx_Cu'         , vflx_out,                 ncol,    lchnk )

     ! call outfld( 'qtten_Cu'        , qtten_out,                ncol,    lchnk )
     ! call outfld( 'slten_Cu'        , slten_out,                ncol,    lchnk )
     ! call outfld( 'uten_Cu'         , uten_out(:,pver:1:-1),     ncol,    lchnk )
     ! call outfld( 'vten_Cu'         , vten_out(:,pver:1:-1),     ncol,    lchnk )
     ! call outfld( 'qvten_Cu'        , qvten_out(:,pver:1:-1),    ncol,    lchnk )
     ! call outfld( 'qlten_Cu'        , qlten_out(:,pver:1:-1),    ncol,    lchnk )
     ! call outfld( 'qiten_Cu'        , qiten_out(:,pver:1:-1),    ncol,    lchnk )

     ! call outfld( 'cbmf_Cu'         , cbmf_out,                 ncol,    lchnk )
     ! call outfld( 'ufrcinvbase_Cu'  , ufrcinvbase_out,          ncol,    lchnk )
     ! call outfld( 'ufrclcl_Cu'      , ufrclcl_out,              ncol,    lchnk )
     ! call outfld( 'winvbase_Cu'     , winvbase_out,             ncol,    lchnk )
     ! call outfld( 'wlcl_Cu'         , wlcl_out,                 ncol,    lchnk )
     ! call outfld( 'plcl_Cu'         , plcl_out,                 ncol,    lchnk )
     ! call outfld( 'pinv_Cu'         , pinv_out,                 ncol,    lchnk )
     ! call outfld( 'plfc_Cu'         , plfc_out,                 ncol,    lchnk )
     ! call outfld( 'pbup_Cu'         , pbup_out,                 ncol,    lchnk )
     ! call outfld( 'ppen_Cu'         , ppen_out,                 ncol,    lchnk )
     ! call outfld( 'qtsrc_Cu'        , qtsrc_out,                ncol,    lchnk )
     ! call outfld( 'thlsrc_Cu'       , thlsrc_out,               ncol,    lchnk )
     ! call outfld( 'thvlsrc_Cu'      , thvlsrc_out,              ncol,    lchnk )
     ! call outfld( 'emfkbup_Cu'      , emfkbup_out,              ncol,    lchnk )
     ! call outfld( 'cin_Cu'          , cinh_out,                 ncol,    lchnk )
     ! call outfld( 'cinlcl_Cu'       , cinlclh_out,              ncol,    lchnk )
     ! call outfld( 'cbmflimit_Cu'    , cbmflimit_out,            ncol,    lchnk )
     ! call outfld( 'tkeavg_Cu'       , tkeavg_out,               ncol,    lchnk )
     ! call outfld( 'zinv_Cu'         , zinv_out,                 ncol,    lchnk )
     ! call outfld( 'rcwp_Cu'         , rcwp_out,                 ncol,    lchnk )
     ! call outfld( 'rlwp_Cu'         , rlwp_out,                 ncol,    lchnk )
     ! call outfld( 'riwp_Cu'         , riwp_out,                 ncol,    lchnk )
     ! call outfld( 'tophgt_Cu'       , cush_inout,               ncol,    lchnk )

     ! call outfld( 'wu_Cu'           , wu_out,                   ncol,    lchnk )
     ! call outfld( 'ufrc_Cu'         , ufrc_out,                 ncol,    lchnk )
     ! call outfld( 'qtu_Cu'          , qtu_out,                  ncol,    lchnk )
     ! call outfld( 'thlu_Cu'         , thlu_out,                 ncol,    lchnk )
     ! call outfld( 'thvu_Cu'         , thvu_out,                 ncol,    lchnk )
     ! call outfld( 'uu_Cu'           , uu_out,                   ncol,    lchnk )
     ! call outfld( 'vu_Cu'           , vu_out,                   ncol,    lchnk )
     ! call outfld( 'qtu_emf_Cu'      , qtu_emf_out,              ncol,    lchnk )
     ! call outfld( 'thlu_emf_Cu'     , thlu_emf_out,             ncol,    lchnk )
     ! call outfld( 'uu_emf_Cu'       , uu_emf_out,               ncol,    lchnk )
     ! call outfld( 'vu_emf_Cu'       , vu_emf_out,               ncol,    lchnk )
     ! call outfld( 'umf_Cu'          , umf_out(:,pver:0:-1),      ncol,    lchnk )
     ! call outfld( 'uemf_Cu'         , uemf_out,                 ncol,    lchnk )
     ! call outfld( 'qcu_Cu'          , qcu_out(:,pver:1:-1),      ncol,    lchnk )
     ! call outfld( 'qlu_Cu'          , qlu_out(:,pver:1:-1),      ncol,    lchnk )
     ! call outfld( 'qiu_Cu'          , qiu_out(:,pver:1:-1),      ncol,    lchnk )
     ! call outfld( 'cufrc_Cu'        , cufrc_out(:,pver:1:-1),    ncol,    lchnk )
     ! call outfld( 'fer_Cu'          , fer_out,                  ncol,    lchnk )
     ! call outfld( 'fdr_Cu'          , fdr_out,                  ncol,    lchnk )

     ! call outfld( 'dwten_Cu'        , dwten_out,                ncol,    lchnk )
     ! call outfld( 'diten_Cu'        , diten_out,                ncol,    lchnk )
     ! call outfld( 'qrten_Cu'        , qrten_out(:,pver:1:-1),    ncol,    lchnk )
     ! call outfld( 'qsten_Cu'        , qsten_out(:,pver:1:-1),    ncol,    lchnk )
     ! call outfld( 'flxrain_Cu'      , flxrain_out,              ncol,    lchnk )
     ! call outfld( 'flxsnow_Cu'      , flxsnow_out,              ncol,    lchnk )
     ! call outfld( 'ntraprd_Cu'      , ntraprd_out,              ncol,    lchnk )
     ! call outfld( 'ntsnprd_Cu'      , ntsnprd_out,              ncol,    lchnk )

     ! call outfld( 'excessu_Cu'      , excessu_arr_out,          ncol,    lchnk )
     ! call outfld( 'excess0_Cu'      , excess0_arr_out,          ncol,    lchnk )
     ! call outfld( 'xc_Cu'           , xc_arr_out,               ncol,    lchnk )
     ! call outfld( 'aquad_Cu'        , aquad_arr_out,            ncol,    lchnk )
     ! call outfld( 'bquad_Cu'        , bquad_arr_out,            ncol,    lchnk )
     ! call outfld( 'cquad_Cu'        , cquad_arr_out,            ncol,    lchnk )
     ! call outfld( 'bogbot_Cu'       , bogbot_arr_out,           ncol,    lchnk )
     ! call outfld( 'bogtop_Cu'       , bogtop_arr_out,           ncol,    lchnk )

     ! call outfld( 'exit_UWCu_Cu'    , exit_UWCu,                ncol,    lchnk )
     ! call outfld( 'exit_conden_Cu'  , exit_conden,              ncol,    lchnk )
     ! call outfld( 'exit_klclpver_Cu' , exit_klclpver,             ncol,    lchnk )
     ! call outfld( 'exit_klfcpver_Cu' , exit_klfcpver,             ncol,    lchnk )
     ! call outfld( 'exit_ufrc_Cu'    , exit_ufrc,                ncol,    lchnk )
     ! call outfld( 'exit_wtw_Cu'     , exit_wtw,                 ncol,    lchnk )
     ! call outfld( 'exit_drycore_Cu' , exit_drycore,             ncol,    lchnk )
     ! call outfld( 'exit_wu_Cu'      , exit_wu,                  ncol,    lchnk )
     ! call outfld( 'exit_cufilter_Cu', exit_cufilter,            ncol,    lchnk )
     ! call outfld( 'exit_kinv1_Cu'   , exit_kinv1,               ncol,    lchnk )
     ! call outfld( 'exit_rei_Cu'     , exit_rei,                 ncol,    lchnk )

     ! call outfld( 'limit_shcu_Cu'   , limit_shcu,               ncol,    lchnk )
     ! call outfld( 'limit_negcon_Cu' , limit_negcon,             ncol,    lchnk )
     ! call outfld( 'limit_ufrc_Cu'   , limit_ufrc,               ncol,    lchnk )
     ! call outfld( 'limit_ppen_Cu'   , limit_ppen,               ncol,    lchnk )
     ! call outfld( 'limit_emf_Cu'    , limit_emf,                ncol,    lchnk )
     ! call outfld( 'limit_cinlcl_Cu' , limit_cinlcl,             ncol,    lchnk )
     ! call outfld( 'limit_cin_Cu'    , limit_cin,                ncol,    lchnk )
     ! call outfld( 'limit_cbmf_Cu'   , limit_cbmf,               ncol,    lchnk )
     ! call outfld( 'limit_rei_Cu'    , limit_rei,                ncol,    lchnk )
     ! call outfld( 'ind_delcin_Cu'   , ind_delcin,               ncol,    lchnk )
  end subroutine uwshcu_cam

  end module uwshcu

