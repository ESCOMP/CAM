! PEVERWHEE - dependencies = interpolate_data
!> \file rrtmgp_lw_cloud_optics.F90
!!

!> This module contains two routines: The first initializes data and functions
!! needed to compute the longwave cloud radiative properteis in RRTMGP. The second routine
!! is a ccpp scheme within the "radiation loop", where the shortwave optical properties
!! (optical-depth, single-scattering albedo, asymmetry parameter) are computed for ALL
!! cloud types visible to RRTMGP.
module rrtmgp_lw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
!  use mo_cloud_optics,          only: ty_cloud_optics
!  use rrtmgp_lw_gas_optics,     only: lw_gas_props
  use interpolate_data,         only: interp_type, lininterp_init, &
                                      lininterp, extrap_method_bndry, &
                                      lininterp_finish
  use radiation_utils,          only: get_mu_lambda_weights_ccpp
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp 
  use mo_optical_props,         only: ty_optical_props_1scl

  implicit none
  public :: rrtmgp_lw_cloud_optics_run

  real(kind_phys), allocatable :: abs_lw_liq(:,:,:)
  real(kind_phys), allocatable :: abs_lw_ice(:,:)
  real(kind_phys), allocatable :: g_mu(:)
  real(kind_phys), allocatable :: g_d_eff(:)
  real(kind_phys), allocatable :: g_lambda(:,:)
  real(kind_phys) :: tiny
  integer :: nmu
  integer :: nlambda
  integer :: n_g_d


contains

  ! ######################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_init()
  ! ######################################################################################
!> \section arg_table_rrtmgp_lw_cloud_optics_init Argument Table
!! \htmlinclude rrtmgp_lw_cloud_optics_init.html
!!
  subroutine rrtmgp_lw_cloud_optics_init(nmu_in, nlambda_in, n_g_d_in, &
                  abs_lw_liq_in, abs_lw_ice_in, nlwbands, g_mu_in, g_lambda_in,  &
                  g_d_eff_in, tiny_in, errmsg, errflg)
    ! Inputs
    integer,          intent(in)  :: nmu_in
    integer,          intent(in)  :: nlambda_in
    integer,          intent(in)  :: n_g_d_in
    integer,          intent(in)  :: nlwbands
    real(kind_phys), dimension(:,:,:), intent(in) :: abs_lw_liq_in
    real(kind_phys), dimension(:,:),   intent(in) :: abs_lw_ice_in
    real(kind_phys), dimension(:,:),   intent(in) :: g_lambda_in
    real(kind_phys), dimension(:),     intent(in) :: g_mu_in
    real(kind_phys), dimension(:),     intent(in) :: g_d_eff_in
    real(kind_phys),                   intent(in) :: tiny_in

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    character(len=256) :: alloc_errmsg
    character(len=*), parameter :: sub = 'rrtmgp_lw_cloud_optics_init'

    ! Set error variables
    errmsg = ''
    errflg = 0

    ! Set module-level variables
    nmu = nmu_in
    nlambda = nlambda_in
    n_g_d = n_g_d_in
    tiny = tiny_in
    ! Allocate module-level-variables
    allocate(abs_lw_liq(nmu,nlambda,nlwbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating abs_lw_liq, message: ', alloc_errmsg
       return
    end if
    allocate(abs_lw_ice(n_g_d,nlwbands), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating abs_lw_ice, message: ', alloc_errmsg
       return
    end if
    allocate(g_mu(nmu), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating g_mu, message: ', alloc_errmsg
       return
    end if
    allocate(g_lambda(nmu,nlambda), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating g_lambda, message: ', alloc_errmsg
       return
    end if
    allocate(g_d_eff(n_g_d), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating g_d_eff, message: ', alloc_errmsg
       return
    end if

    abs_lw_liq = abs_lw_liq_in
    abs_lw_ice = abs_lw_ice_in
    g_mu       = g_mu_in
    g_lambda   = g_lambda_in
    g_d_eff    = g_d_eff_in

  end subroutine rrtmgp_lw_cloud_optics_init

  ! ######################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_run()
  ! ######################################################################################
!> \section arg_table_rrtmgp_lw_cloud_optics_run Argument Table
!! \htmlinclude rrtmgp_lw_cloud_optics_run.html
!!
  subroutine rrtmgp_lw_cloud_optics_run(ncol, nlay, nlaycam, cld, cldfsnow, cldfgrau,       &
             cldfprime, graupel_in_rad, kdist_lw, cloud_lw, lamc, pgam, iclwpth, iciwpth,  &
             dei, icswpth, des, icgrauwpth, degrau, nlwbands, do_snow, &
             do_graupel, cld_lw_abs_cloudsim, snow_lw_abs_cloudsim, pver, ktopcam,               &
             grau_lw_abs_cloudsim, idx_lw_cloudsim, tauc, cldf, errmsg, errflg)
    ! Compute combined cloud optical properties
    ! Create MCICA stochastic arrays for cloud LW optical properties
    ! Initialize optical properties object (cloud_lw) and load with MCICA columns

    ! Inputs
    integer,                           intent(in) :: ncol
    integer,                           intent(in) :: nlay
    integer,                           intent(in) :: nlaycam
    integer,                           intent(in) :: nlwbands
    integer,                           intent(in) :: pver
    integer,                           intent(in) :: ktopcam
    integer,                           intent(in) :: idx_lw_cloudsim
    real(kind_phys), dimension(:,:),   intent(in) :: cld
    real(kind_phys), dimension(:,:),   intent(in) :: cldfsnow
    real(kind_phys), dimension(:,:),   intent(in) :: cldfgrau
    real(kind_phys), dimension(:,:),   intent(in) :: cldfprime
    real(kind_phys), dimension(:,:),   intent(in) :: lamc
    real(kind_phys), dimension(:,:),   intent(in) :: pgam
    real(kind_phys), dimension(:,:),   intent(in) :: iclwpth
    real(kind_phys), dimension(:,:),   intent(in) :: iciwpth
    real(kind_phys), dimension(:,:),   intent(in) :: icswpth
    real(kind_phys), dimension(:,:),   intent(in) :: icgrauwpth
    real(kind_phys), dimension(:,:),   intent(in) :: dei
    real(kind_phys), dimension(:,:),   intent(in) :: des
    real(kind_phys), dimension(:,:),   intent(in) :: degrau
    logical,                           intent(in) :: graupel_in_rad
    logical,                           intent(in) :: do_snow
    logical,                           intent(in) :: do_graupel
    class(ty_gas_optics_rrtmgp),       intent(in) :: kdist_lw

    ! Outputs
    type(ty_optical_props_1scl),       intent(out) :: cloud_lw
    real(kind_phys), dimension(:,:),   intent(out) :: cld_lw_abs_cloudsim
    real(kind_phys), dimension(:,:),   intent(out) :: snow_lw_abs_cloudsim
    real(kind_phys), dimension(:,:),   intent(out) :: grau_lw_abs_cloudsim
    real(kind_phys), dimension(:,:),   intent(out) :: cldf
    real(kind_phys), dimension(:,:,:), intent(out) :: tauc
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    ! Local variables
    integer :: idx, kdx

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(kind_phys) :: liq_lw_abs(nlwbands, ncol, pver)   ! liquid absorption optics depth (LW)
    real(kind_phys) :: ice_lw_abs(nlwbands, ncol, pver)   ! ice absorption optics depth (LW)
    real(kind_phys) :: cld_lw_abs(nlwbands, ncol, pver)   ! cloud absorption optics depth (LW)
    real(kind_phys) :: snow_lw_abs(nlwbands, ncol, pver)  ! snow absorption optics depth (LW)
    real(kind_phys) :: grau_lw_abs(nlwbands, ncol, pver)  ! graupel absorption optics depth (LW)
    real(kind_phys) :: c_cld_lw_abs(nlwbands, ncol, pver) ! combined cloud absorption optics depth (LW)

    character(len=*), parameter :: sub = 'rrtmgp_set_cloud_lw'
    !--------------------------------------------------------------------------------

    ! Combine the cloud optical properties.

    ! gammadist liquid optics
    call liquid_cloud_get_rad_props_lw(ncol, pver, nmu, nlambda, nlwbands, lamc, pgam, g_mu, g_lambda, iclwpth, &
            abs_lw_liq, liq_lw_abs, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    ! Mitchell ice optics
    call ice_cloud_get_rad_props_lw(ncol, pver, nlwbands, iciwpth, dei, n_g_d, g_d_eff, abs_lw_ice, ice_lw_abs, &
            errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    cld_lw_abs(:,:,:) = liq_lw_abs(:,:,:) + ice_lw_abs(:,:,:)

    if (do_snow) then
       ! add in snow
       call snow_cloud_get_rad_props_lw(ncol, pver, nlwbands, icswpth, des, n_g_d, g_d_eff, abs_lw_ice, &
               snow_lw_abs, errmsg, errflg)
       if (errflg /= 0) then
          return
       end if
       do idx = 1, ncol
          do kdx = 1, pver
             if (cldfprime(idx,kdx) > 0._kind_phys) then
                c_cld_lw_abs(:,idx,kdx) = ( cldfsnow(idx,kdx)*snow_lw_abs(:,idx,kdx) &
                                                 + cld(idx,kdx)*cld_lw_abs(:,idx,kdx) )/cldfprime(idx,kdx)
             else
                c_cld_lw_abs(:,idx,kdx) = 0._kind_phys
             end if
          end do
       end do
    else
       c_cld_lw_abs(:,:,:) = cld_lw_abs(:,:,:)
    end if

    ! add in graupel
    if (do_graupel .and. graupel_in_rad) then
       call grau_cloud_get_rad_props_lw(ncol, pver, nlwbands, icgrauwpth, degrau, n_g_d, g_d_eff, abs_lw_ice, &
               grau_lw_abs, errmsg, errflg)
       if (errflg /= 0) then
          return
       end if
       do idx = 1, ncol
          do kdx = 1, pver
             if (cldfprime(idx,kdx) > 0._kind_phys) then
                c_cld_lw_abs(:,idx,kdx) = ( cldfgrau(idx,kdx)*grau_lw_abs(:,idx,kdx) &
                                                 + cld(idx,kdx)*c_cld_lw_abs(:,idx,kdx) )/cldfprime(idx,kdx)
             else
                c_cld_lw_abs(:,idx,kdx) = 0._kind_phys
             end if
          end do
       end do
    end if

    ! Cloud optics for COSP
    cld_lw_abs_cloudsim = cld_lw_abs(idx_lw_cloudsim,:,:)
    snow_lw_abs_cloudsim = snow_lw_abs(idx_lw_cloudsim,:,:)
    grau_lw_abs_cloudsim = grau_lw_abs(idx_lw_cloudsim,:,:)

    ! Extract just the layers of CAM where RRTMGP does calculations

    ! Subset "chunk" data so just the number of CAM layers in the
    ! radiation calculation are used by MCICA to produce subcolumns
    cldf = cldfprime(:, ktopcam:)
    tauc = c_cld_lw_abs(:, :, ktopcam:)

    ! Enforce tauc >= 0.
    tauc = merge(tauc, 0.0_kind_phys, tauc > 0.0_kind_phys)

    errmsg =cloud_lw%alloc_1scl(ncol, nlay, kdist_lw)
    if (len_trim(errmsg) > 0) then
       errflg = 1
       return
    end if

  end subroutine rrtmgp_lw_cloud_optics_run

!==============================================================================

  subroutine liquid_cloud_get_rad_props_lw(ncol, pver, nmu, nlambda, nlwbands, lamc, pgam, &
                  g_mu, g_lambda, iclwpth, abs_lw_liq, abs_od, errmsg, errflg)
    ! Inputs
    integer,                           intent(in) :: ncol
    integer,                           intent(in) :: pver
    integer,                           intent(in) :: nmu
    integer,                           intent(in) :: nlambda
    integer,                           intent(in) :: nlwbands
    real(kind_phys), dimension(:,:),   intent(in) :: lamc
    real(kind_phys), dimension(:,:),   intent(in) :: pgam
    real(kind_phys), dimension(:,:,:), intent(in) :: abs_lw_liq
    real(kind_phys), dimension(:),     intent(in) :: g_mu
    real(kind_phys), dimension(:,:),   intent(in) :: g_lambda
    real(kind_phys), dimension(:,:),   intent(in) :: iclwpth
    ! Outputs
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer lwband, idx, kdx

    ! Set error variables
    errflg = 0
    errmsg = ''

    abs_od = 0._kind_phys

    do kdx = 1,pver
       do idx = 1,ncol
          if(lamc(idx,kdx) > 0._kind_phys) then ! This seems to be the clue for no cloud from microphysics formulation
             call gam_liquid_lw(nlwbands, nmu, nlambda, iclwpth(idx,kdx), lamc(idx,kdx), pgam(idx,kdx), abs_lw_liq, &
                     g_mu, g_lambda, abs_od(1:nlwbands,idx,kdx), errmsg, errflg)
          else
             abs_od(1:nlwbands,idx,kdx) = 0._kind_phys
          endif
       enddo
    enddo

  end subroutine liquid_cloud_get_rad_props_lw

!==============================================================================

  subroutine gam_liquid_lw(nlwbands, nmu, nlambda, clwptn, lamc, pgam, abs_lw_liq, g_mu, g_lambda, abs_od, errmsg, errflg)
    ! Inputs
    integer,         intent(in) :: nlwbands
    integer,         intent(in) :: nmu
    integer,         intent(in) :: nlambda
    real(kind_phys), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
    real(kind_phys), intent(in) :: lamc   ! prognosed value of lambda for cloud
    real(kind_phys), intent(in) :: pgam   ! prognosed value of mu for cloud
    real(kind_phys), dimension(:,:,:), intent(in) :: abs_lw_liq
    real(kind_phys), dimension(:),     intent(in) :: g_mu
    real(kind_phys), dimension(:,:)  , intent(in) :: g_lambda
    ! Outputs
    real(kind_phys), dimension(:), intent(out) :: abs_od
    integer,                       intent(out) :: errflg
    character(len=*),              intent(out) :: errmsg
    
    integer :: lwband ! sw band index

    type(interp_type) :: mu_wgts
    type(interp_type) :: lambda_wgts

    if (clwptn < tiny) then
      abs_od = 0._kind_phys
      return
    endif

    call get_mu_lambda_weights_ccpp(nmu, nlambda, g_mu, g_lambda, lamc, pgam, mu_wgts, lambda_wgts, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    do lwband = 1, nlwbands
       call lininterp(abs_lw_liq(:,:,lwband), nmu, nlambda, &
            abs_od(lwband:lwband), 1, mu_wgts, lambda_wgts)
    enddo

    abs_od = clwptn * abs_od

    call lininterp_finish(mu_wgts)
    call lininterp_finish(lambda_wgts)

  end subroutine gam_liquid_lw

!==============================================================================

 subroutine ice_cloud_get_rad_props_lw(ncol, pver, nlwbands, iciwpth, dei, &
                 n_g_d, g_d_eff, abs_lw_ice, abs_od, errmsg, errflg)
    integer,                           intent(in)  :: ncol
    integer,                           intent(in)  :: pver
    integer,                           intent(in)  :: n_g_d
    integer,                           intent(in)  :: nlwbands
    real(kind_phys), dimension(:),     intent(in)  :: g_d_eff
    real(kind_phys), dimension(:,:),   intent(in)  :: iciwpth
    real(kind_phys), dimension(:,:),   intent(in)  :: dei
    real(kind_phys), dimension(:,:),   intent(in)  :: abs_lw_ice
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    ! Set error variables
    errflg = 0
    errmsg = ''

    call interpolate_ice_optics_lw(ncol, pver, nlwbands, iciwpth, dei, &
            n_g_d, g_d_eff, abs_lw_ice, abs_od, errmsg, errflg)

 end subroutine ice_cloud_get_rad_props_lw

!==============================================================================

 subroutine snow_cloud_get_rad_props_lw(ncol, pver, nlwbands, icswpth, des,  &
                 n_g_d, g_d_eff, abs_lw_ice, abs_od, errmsg, errflg)
    integer,                           intent(in)  :: ncol
    integer,                           intent(in)  :: pver
    integer,                           intent(in)  :: n_g_d
    integer,                           intent(in)  :: nlwbands
    real(kind_phys), dimension(:),     intent(in)  :: g_d_eff
    real(kind_phys), dimension(:,:),   intent(in)  :: icswpth
    real(kind_phys), dimension(:,:),   intent(in)  :: des
    real(kind_phys), dimension(:,:),   intent(in)  :: abs_lw_ice
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    errflg = 0
    errmsg = ''

    call interpolate_ice_optics_lw(ncol, pver, nlwbands, icswpth, des, &
            n_g_d, g_d_eff, abs_lw_ice, abs_od, errmsg, errflg)

 end subroutine snow_cloud_get_rad_props_lw

!==============================================================================

 subroutine grau_cloud_get_rad_props_lw(ncol, pver, nlwbands, icgrauwpth, degrau, &
                 n_g_d, g_d_eff, abs_lw_ice, abs_od, errmsg, errflg)
    integer,                           intent(in)  :: ncol
    integer,                           intent(in)  :: pver
    integer,                           intent(in)  :: n_g_d
    integer,                           intent(in)  :: nlwbands
    real(kind_phys), dimension(:),     intent(in)  :: g_d_eff
    real(kind_phys), dimension(:,:),   intent(in)  :: icgrauwpth
    real(kind_phys), dimension(:,:),   intent(in)  :: degrau
    real(kind_phys), dimension(:,:),   intent(in)  :: abs_lw_ice
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    ! This does the same thing as ice_cloud_get_rad_props_lw, except with a
    ! different water path and effective diameter.
    call interpolate_ice_optics_lw(ncol, pver, nlwbands, icgrauwpth, degrau, n_g_d, &
            g_d_eff, abs_lw_ice, abs_od, errmsg, errflg)
 
 end subroutine grau_cloud_get_rad_props_lw

!==============================================================================

  subroutine interpolate_ice_optics_lw(ncol, pver, nlwbands, iciwpth, dei, &
                  n_g_d, g_d_eff, abs_lw_ice, abs_od, errmsg, errflg)

    integer,           intent(in)                  :: ncol
    integer,           intent(in)                  :: n_g_d
    integer,           intent(in)                  :: pver
    integer,           intent(in)                  :: nlwbands
    real(kind_phys), dimension(:),     intent(in)  :: g_d_eff
    real(kind_phys), dimension(:,:),   intent(in)  :: iciwpth
    real(kind_phys), dimension(:,:),   intent(in)  :: dei
    real(kind_phys), dimension(:,:),   intent(in)  :: abs_lw_ice
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    type(interp_type) :: dei_wgts

    integer :: idx, kdx, lwband
    real(kind_phys) :: absor(nlwbands)

    ! Set error variables
    errflg = 0
    errmsg = ''

    do kdx = 1,pver
       do idx = 1,ncol
          ! if ice water path is too small, OD := 0
          if( iciwpth(idx,kdx) < tiny .or. dei(idx,kdx) == 0._kind_phys) then
             abs_od (:,idx,kdx) = 0._kind_phys
          else
             ! for each cell interpolate to find weights in g_d_eff grid.
             call lininterp_init(g_d_eff, n_g_d, dei(idx:idx,kdx), 1, &
                  extrap_method_bndry, dei_wgts)
             ! interpolate into grid and extract radiative properties
             do lwband = 1, nlwbands
                call lininterp(abs_lw_ice(:,lwband), n_g_d, &
                     absor(lwband:lwband), 1, dei_wgts)
             enddo
             abs_od(:,idx,kdx) = iciwpth(idx,kdx) * absor
             where(abs_od(:,idx,kdx) > 50.0_kind_phys) abs_od(:,idx,kdx) = 50.0_kind_phys
             call lininterp_finish(dei_wgts)
          endif
       enddo
    enddo

  end subroutine interpolate_ice_optics_lw

!==============================================================================

end module rrtmgp_lw_cloud_optics
