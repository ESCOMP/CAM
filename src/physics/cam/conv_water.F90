  module conv_water

   ! --------------------------------------------------------------------- ! 
   ! Purpose:                                                              !
   ! Computes grid-box average liquid (and ice) from stratus and cumulus   !
   ! These values used by both the radiation and the COSP diagnostics.     !
   !                                                                       ! 
   ! Method:                                                               !
   ! Extract information about deep+shallow liquid and cloud fraction from !
   ! the physics buffer.                                                   !
   !                                                                       !
   ! Author: Rich Neale, August 2006                                       !
   !         October 2006: Allow averaging of liquid to give a linear      !
   !                       average in emissivity.                          !
   !         Andrew Gettelman October 2010  Separate module                !
   !---------------------------------------------------------------------- !

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use spmd_utils,     only: masterproc
  use ppgrid,         only: pcols, pver, pverp
  use physconst,      only: gravit
  use cam_abortutils, only: endrun

  use perf_mod
  use cam_logfile,    only: iulog

  implicit none
  private
  save

  public :: &
     conv_water_readnl,   &
     conv_water_register, &
     conv_water_init,     &
     conv_water_4rad,     &
     conv_water_in_rad

! pbuf indices

  integer :: icwmrsh_idx, icwmrdp_idx, sh_frac_idx, dp_frac_idx, &
             ast_idx, rei_idx

  integer :: ixcldice, ixcldliq
  integer :: gb_totcldliqmr_idx, gb_totcldicemr_idx

! Namelist
integer, parameter :: unset_int = huge(1)

integer  :: conv_water_in_rad = unset_int  ! 0==> No; 1==> Yes-Arithmetic average;
                                           ! 2==> Yes-Average in emissivity.
integer  :: conv_water_mode
real(r8) :: frac_limit 

!=============================================================================================
contains
!=============================================================================================

subroutine conv_water_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'conv_water_readnl'

   real(r8) :: conv_water_frac_limit = 0._r8

   namelist /conv_water_nl/ conv_water_in_rad, conv_water_frac_limit
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'conv_water_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, conv_water_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(conv_water_in_rad,     1, mpiint, 0, mpicom)
   call mpibcast(conv_water_frac_limit, 1, mpir8,  0, mpicom)
#endif

   conv_water_mode = conv_water_in_rad
   frac_limit      = conv_water_frac_limit

end subroutine conv_water_readnl

!=============================================================================================

  subroutine conv_water_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the fields in the physics buffer.                            !
  !                                                                       !
  !---------------------------------------------------------------------- !

    use constituents, only: cnst_add, pcnst
    use physconst,    only: mwdry, cpair
    
    use physics_buffer, only : pbuf_add_field, dtype_r8

  !-----------------------------------------------------------------------

    ! grid box total cloud liquid water mixing ratio (kg/kg)
    call pbuf_add_field('GB_TOTCLDLIQMR', 'physpkg', dtype_r8, (/pcols,pver/), gb_totcldliqmr_idx)  
    ! grid box total cloud ice water mixing ratio (kg/kg)
    call pbuf_add_field('GB_TOTCLDICEMR', 'physpkg', dtype_r8, (/pcols,pver/), gb_totcldicemr_idx)  

  end subroutine conv_water_register


  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine conv_water_init()
   ! --------------------------------------------------------------------- ! 
   ! Purpose:                                                              !
   !   Initializes the pbuf indices required by conv_water
   ! --------------------------------------------------------------------- ! 

   
   use physics_buffer, only : pbuf_get_index
   use cam_history,    only : addfld

   use constituents,  only: cnst_get_ind

   implicit none

   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)
 
   icwmrsh_idx  = pbuf_get_index('ICWMRSH')
   icwmrdp_idx  = pbuf_get_index('ICWMRDP')
   sh_frac_idx  = pbuf_get_index('SH_FRAC')
   dp_frac_idx  = pbuf_get_index('DP_FRAC')
   ast_idx      = pbuf_get_index('AST')
   rei_idx      = pbuf_get_index('REI')

   ! Convective cloud water variables.
   call addfld ('ICIMRCU',  (/ 'lev' /), 'A', 'kg/kg', 'Convection in-cloud ice mixing ratio '   )
   call addfld ('ICLMRCU',  (/ 'lev' /), 'A', 'kg/kg', 'Convection in-cloud liquid mixing ratio ')
   call addfld ('ICIMRTOT', (/ 'lev' /), 'A', 'kg/kg', 'Total in-cloud ice mixing ratio '        )
   call addfld ('ICLMRTOT', (/ 'lev' /), 'A', 'kg/kg', 'Total in-cloud liquid mixing ratio '     )

   call addfld ('GCLMRDP',  (/ 'lev' /), 'A', 'kg/kg', 'Grid-mean deep convective LWC'           )
   call addfld ('GCIMRDP',  (/ 'lev' /), 'A', 'kg/kg', 'Grid-mean deep convective IWC'           )
   call addfld ('GCLMRSH',  (/ 'lev' /), 'A', 'kg/kg', 'Grid-mean shallow convective LWC'        )
   call addfld ('GCIMRSH',  (/ 'lev' /), 'A', 'kg/kg', 'Grid-mean shallow convective IWC'        )
   call addfld ('FRESH',    (/ 'lev' /), 'A', '1', 'Fractional occurrence of shallow cumulus with condensate')
   call addfld ('FREDP',    (/ 'lev' /), 'A', '1', 'Fractional occurrence of deep cumulus with condensate')
   call addfld ('FRECU',    (/ 'lev' /), 'A', '1', 'Fractional occurrence of cumulus with condensate')
   call addfld ('FRETOT',   (/ 'lev' /), 'A', '1', 'Fractional occurrence of cloud with condensate')

   end subroutine conv_water_init

   subroutine conv_water_4rad(state, pbuf)

   ! --------------------------------------------------------------------- !
   ! Purpose:                                                              !
   ! Computes grid-box average liquid (and ice) from stratus and cumulus   !
   ! Just for the purposes of radiation.                                   !
   !                                                                       !
   ! Method:                                                               !
   ! Extract information about deep+shallow liquid and cloud fraction from !
   ! the physics buffer; the science lives in the portable core            !
   ! (convective_cloud_water.F90).                                         !
   !                                                                       !
   ! Author: Rich Neale, August 2006                                       !
   !         October 2006: Allow averaging of liquid to give a linear      !
   !                       average in emissivity.                          !
   !                                                                       !
   !---------------------------------------------------------------------- !


   use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

   use physics_types,   only: physics_state
   use cam_history,     only: outfld
   use phys_control,    only: phys_getopts
   use convective_cloud_water, only: convective_cloud_water_run

   implicit none

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !


   type(physics_state), target, intent(in) :: state        ! state variables
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   ! --------------- !
   ! Local Workspace !
   ! --------------- !

   ! Physics buffer fields
   real(r8), pointer, dimension(:,:) ::  ast      ! Physical liquid+ice stratus cloud fraction
   real(r8), pointer, dimension(:,:) ::  sh_frac  ! Shallow convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  dp_frac  ! Deep convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  rei      ! Ice effective drop size (microns)

   real(r8), pointer, dimension(:,:) ::  dp_icwmr ! Deep conv. cloud water
   real(r8), pointer, dimension(:,:) ::  sh_icwmr ! Shallow conv. cloud water

   real(r8), pointer, dimension(:,:) :: totg_ice  ! Grid box total cloud ice mixing ratio
   real(r8), pointer, dimension(:,:) :: totg_liq  ! Grid box total cloud liquid mixing ratio

   real(r8) :: conv_ice(pcols,pver)               ! Convective contributions to IC cloud ice
   real(r8) :: conv_liq(pcols,pver)               ! Convective contributions to IC cloud liquid
   real(r8) :: tot_ice(pcols,pver)                ! Total IC ice
   real(r8) :: tot_liq(pcols,pver)                ! Total IC liquid

   integer  :: itim_old                           ! Time index buff stuff.

   real(r8) :: totg_ice_sh(pcols,pver)   ! Grid-mean IWP from shallow convective cloud
   real(r8) :: totg_liq_sh(pcols,pver)   ! Grid-mean LWP from shallow convective cloud
   real(r8) :: totg_ice_dp(pcols,pver)   ! Grid-mean IWP from deep convective cloud
   real(r8) :: totg_liq_dp(pcols,pver)   ! Grid-mean LWP from deep convective cloud
   real(r8) :: fresh(pcols,pver)         ! Fractional occurrence of shallow cumulus
   real(r8) :: fredp(pcols,pver)         ! Fractional occurrence of deep cumulus
   real(r8) :: frecu(pcols,pver)         ! Fractional occurrence of cumulus
   real(r8) :: fretot(pcols,pver)        ! Fractional occurrence of cloud

   character(len=512) :: errmsg
   integer :: errflg

   integer :: lchnk
   integer :: ncol

   character(len=16) :: microp_scheme

   ncol  = state%ncol
   lchnk = state%lchnk

 ! Get microphysics option
   call phys_getopts( microp_scheme_out = microp_scheme )

 ! Get convective in-cloud water.

   call pbuf_get_field(pbuf, icwmrsh_idx, sh_icwmr )
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr )

 ! Get convective in-cloud fraction

   call pbuf_get_field(pbuf, sh_frac_idx,  sh_frac )
   call pbuf_get_field(pbuf, dp_frac_idx,  dp_frac )
   call pbuf_get_field(pbuf, rei_idx,      rei )

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,  ast,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   ! Fields computed below and stored in pbuf.
   call pbuf_get_field(pbuf, gb_totcldicemr_idx, totg_ice)
   call pbuf_get_field(pbuf, gb_totcldliqmr_idx, totg_liq)

   ! Zero the diagnostic arrays so columns beyond ncol hold no garbage.
   conv_ice(:,:)    = 0._r8
   conv_liq(:,:)    = 0._r8
   tot_ice(:,:)     = 0._r8
   tot_liq(:,:)     = 0._r8
   totg_ice_sh(:,:) = 0._r8
   totg_liq_sh(:,:) = 0._r8
   totg_ice_dp(:,:) = 0._r8
   totg_liq_dp(:,:) = 0._r8
   fresh(:,:)  = 0._r8
   fredp(:,:)  = 0._r8
   frecu(:,:)  = 0._r8
   fretot(:,:) = 0._r8

   call convective_cloud_water_run(                       &
      ncol              = ncol,                           &
      pver              = pver,                           &
      conv_water_in_rad = conv_water_mode,                &
      frac_limit        = frac_limit,                     &
      one_mom_clouds    = (microp_scheme == 'RK'),        &
      gravit            = gravit,                         &
      pdel              = state%pdel(:ncol,:),            &
      ls_liq            = state%q(:ncol,:,ixcldliq),      &
      ls_ice            = state%q(:ncol,:,ixcldice),      &
      sh_icwmr          = sh_icwmr(:ncol,:),              &
      dp_icwmr          = dp_icwmr(:ncol,:),              &
      sh_frac           = sh_frac(:ncol,:),               &
      dp_frac           = dp_frac(:ncol,:),               &
      ast               = ast(:ncol,:),                   &
      rei               = rei(:ncol,:),                   &
      totg_liq          = totg_liq(:ncol,:),              &
      totg_ice          = totg_ice(:ncol,:),              &
      conv_liq          = conv_liq(:ncol,:),              &
      conv_ice          = conv_ice(:ncol,:),              &
      tot_liq           = tot_liq(:ncol,:),               &
      tot_ice           = tot_ice(:ncol,:),               &
      totg_liq_sh       = totg_liq_sh(:ncol,:),           &
      totg_liq_dp       = totg_liq_dp(:ncol,:),           &
      totg_ice_sh       = totg_ice_sh(:ncol,:),           &
      totg_ice_dp       = totg_ice_dp(:ncol,:),           &
      fresh             = fresh(:ncol,:),                 &
      fredp             = fredp(:ncol,:),                 &
      frecu             = frecu(:ncol,:),                 &
      fretot            = fretot(:ncol,:),                &
      errmsg            = errmsg,                         &
      errflg            = errflg)
   if (errflg /= 0) then
      call endrun('conv_water_4rad: '//trim(errmsg))
   end if

  ! Output convective IC WMRs

   call outfld( 'ICLMRCU ', conv_liq  , pcols, lchnk )
   call outfld( 'ICIMRCU ', conv_ice  , pcols, lchnk )
   call outfld( 'ICLMRTOT', tot_liq   , pcols, lchnk )
   call outfld( 'ICIMRTOT', tot_ice   , pcols, lchnk )

   call outfld('GCLMRDP', totg_liq_dp, pcols, lchnk)
   call outfld('GCIMRDP', totg_ice_dp, pcols, lchnk)
   call outfld('GCLMRSH', totg_liq_sh, pcols, lchnk)
   call outfld('GCIMRSH', totg_ice_sh, pcols, lchnk)
   call outfld('FRESH',   fresh,       pcols, lchnk)
   call outfld('FREDP',   fredp,       pcols, lchnk)
   call outfld('FRECU',   frecu,       pcols, lchnk)
   call outfld('FRETOT',  fretot,      pcols, lchnk)

  end subroutine conv_water_4rad

end module conv_water
