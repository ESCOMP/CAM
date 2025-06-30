module cloud_rad_props

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
use constituents,     only: cnst_get_ind
use radconstants,     only: nswbands, nlwbands, idx_sw_diag
use rad_constituents, only: iceopticsfile, liqopticsfile
use oldcloud_optics,  only: oldcloud_init, oldcloud_lw, &
                            old_liq_get_rad_props_lw, old_ice_get_rad_props_lw
                            
use slingo_liq_optics,      only: slingo_rad_props_init
use ebert_curry_ice_optics, only: ec_rad_props_init, scalefactor

use interpolate_data, only: interp_type, lininterp_init, lininterp, &
                            extrap_method_bndry, lininterp_finish

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun


implicit none
private
save

public :: &
   cloud_rad_props_init,          &
   cloud_rad_props_get_lw,        & ! return LW optical props for old cloud optics
   get_ice_optics_sw,             & ! return Mitchell SW ice radiative properties
   ice_cloud_get_rad_props_lw,    & ! return Mitchell LW ice radiative properties
   get_liquid_optics_sw,          & ! return Conley SW radiative properties
   liquid_cloud_get_rad_props_lw, & ! return Conley LW radiative properties
   get_snow_optics_sw,            &
   snow_cloud_get_rad_props_lw,   &
   get_grau_optics_sw,            &
   get_mu_lambda_weights,         &
   grau_cloud_get_rad_props_lw


integer :: nmu, nlambda
real(r8), allocatable :: g_mu(:)           ! mu samples on grid
real(r8), allocatable :: g_lambda(:,:)     ! lambda scale samples on grid
real(r8), allocatable :: ext_sw_liq(:,:,:)
real(r8), allocatable :: ssa_sw_liq(:,:,:)
real(r8), allocatable :: asm_sw_liq(:,:,:)
real(r8), allocatable :: abs_lw_liq(:,:,:)

integer :: n_g_d
real(r8), allocatable :: g_d_eff(:)        ! radiative effective diameter samples on grid
real(r8), allocatable :: ext_sw_ice(:,:)
real(r8), allocatable :: ssa_sw_ice(:,:)
real(r8), allocatable :: asm_sw_ice(:,:)
real(r8), allocatable :: abs_lw_ice(:,:)

! indexes into pbuf for optical parameters of MG clouds
integer :: i_dei=0
integer :: i_mu=0
integer :: i_lambda=0
integer :: i_iciwp=0
integer :: i_iclwp=0
integer :: i_des=0
integer :: i_icswp=0
integer :: i_degrau=0
integer :: i_icgrauwp=0

! indexes into constituents for old optics
integer :: &
   ixcldice,           & ! cloud ice water index
   ixcldliq              ! cloud liquid water index

real(r8), parameter :: tiny = 1.e-80_r8

!==============================================================================
contains
!==============================================================================

subroutine cloud_rad_props_init(abs_lw_liq_out, abs_lw_ice_out, ext_sw_liq_out,   &
                  ssa_sw_liq_out, asm_sw_liq_out, ext_sw_ice_out, ssa_sw_ice_out, &
                  asm_sw_ice_out, g_mu_out, g_lambda_out, g_d_eff_out, tiny_out)
   use netcdf
   use spmd_utils,                only: masterproc
   use ioFileMod,                 only: getfil
   use error_messages,            only: handle_ncerr
   use cam_abortutils,            only: handle_allocate_error
   use rrtmgp_cloud_optics_setup, only: rrtmgp_cloud_optics_setup_init
#if ( defined SPMD )
   use mpishorthand
#endif
   real(r8), allocatable, intent(out) :: abs_lw_liq_out(:,:,:)
   real(r8), allocatable, intent(out) :: ext_sw_liq_out(:,:,:)
   real(r8), allocatable, intent(out) :: asm_sw_liq_out(:,:,:)
   real(r8), allocatable, intent(out) :: ssa_sw_liq_out(:,:,:)
   real(r8), allocatable, intent(out) :: abs_lw_ice_out(:,:)
   real(r8), allocatable, intent(out) :: ext_sw_ice_out(:,:)
   real(r8), allocatable, intent(out) :: asm_sw_ice_out(:,:)
   real(r8), allocatable, intent(out) :: ssa_sw_ice_out(:,:)
   real(r8), allocatable, intent(out) :: g_mu_out(:)
   real(r8), allocatable, intent(out) :: g_lambda_out(:,:)
   real(r8), allocatable, intent(out) :: g_d_eff_out(:)
   real(r8),              intent(out) :: tiny_out

   character(len=256) :: liquidfile
   character(len=256) :: icefile
   character(len=256) :: locfn

   integer :: ncid, dimid, f_nlwbands, f_nswbands, ierr
   integer :: vdimids(NF90_MAX_VAR_DIMS), ndims, templen
   ! liquid clouds
   integer :: mudimid, lambdadimid
   integer :: mu_id, lambda_id, ext_sw_liq_id, ssa_sw_liq_id, asm_sw_liq_id, abs_lw_liq_id

   ! ice clouds
   integer :: d_dimid ! diameters
   integer :: d_id, ext_sw_ice_id, ssa_sw_ice_id, asm_sw_ice_id, abs_lw_ice_id

   integer :: err
   character(len=*), parameter :: sub = 'cloud_rad_props_init'
   character(len=512) :: errmsg

   liquidfile = liqopticsfile
   icefile = iceopticsfile

   call slingo_rad_props_init
   call ec_rad_props_init
   call oldcloud_init

   i_dei    = pbuf_get_index('DEI',errcode=err)
   i_mu     = pbuf_get_index('MU',errcode=err)
   i_lambda = pbuf_get_index('LAMBDAC',errcode=err)
   i_iciwp  = pbuf_get_index('ICIWP',errcode=err)
   i_iclwp  = pbuf_get_index('ICLWP',errcode=err)
   i_des    = pbuf_get_index('DES',errcode=err)
   i_icswp  = pbuf_get_index('ICSWP',errcode=err)
   i_icgrauwp  = pbuf_get_index('ICGRAUWP',errcode=err) ! Available when using MG3
   i_degrau    = pbuf_get_index('DEGRAU',errcode=err)   ! Available when using MG3

   ! old optics
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   call rrtmgp_cloud_optics_setup_init(liqopticsfile, iceopticsfile, abs_lw_liq,       &
                  abs_lw_ice, ext_sw_liq, ext_sw_ice, ssa_sw_liq, ssa_sw_ice, asm_sw_liq, &
                  asm_sw_ice, g_lambda, g_mu, g_d_eff, errmsg, err)
   if (err /= 0) then
      call endrun(sub//': '//errmsg)
   end if

   ! Set output variables
   nmu = size(g_mu)
   nlambda = size(g_lambda, 2)
   n_g_d = size(g_d_eff)
   tiny_out = tiny
   allocate(abs_lw_liq_out(nmu,nlambda,nlwbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'abs_lw_liq_out')
   abs_lw_liq_out = abs_lw_liq
   allocate(ext_sw_liq_out(nmu,nlambda,nswbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'ext_sw_liq_out')
   ext_sw_liq_out = ext_sw_liq
   allocate(ssa_sw_liq_out(nmu,nlambda,nswbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'ssa_sw_liq_out')
   ssa_sw_liq_out = ssa_sw_liq
   allocate(asm_sw_liq_out(nmu,nlambda,nswbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'asm_sw_liq_out')
   asm_sw_liq_out = asm_sw_liq
   allocate(abs_lw_ice_out(n_g_d,nlwbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'abs_lw_ice_out')
   abs_lw_ice_out = abs_lw_ice
   allocate(ext_sw_ice_out(n_g_d,nswbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'ext_sw_ice_out')
   ext_sw_ice_out = ext_sw_ice
   allocate(ssa_sw_ice_out(n_g_d,nswbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'ssa_sw_ice_out')
   ssa_sw_ice_out = ssa_sw_ice
   allocate(asm_sw_ice_out(n_g_d,nswbands), stat=ierr)
   call handle_allocate_error(ierr, sub, 'asm_sw_ice_out')
   asm_sw_ice_out = asm_sw_ice
   allocate(g_mu_out(nmu), stat=ierr)
   call handle_allocate_error(ierr, sub, 'g_mu_out')
   g_mu_out = g_mu
   allocate(g_lambda_out(nmu,nlambda), stat=ierr)
   call handle_allocate_error(ierr, sub, 'g_lambda_out')
   g_lambda_out = g_lambda
   allocate(g_d_eff_out(n_g_d), stat=ierr)
   call handle_allocate_error(ierr, sub, 'g_d_eff_out')
   g_d_eff_out = g_d_eff
   return

end subroutine cloud_rad_props_init

!==============================================================================

subroutine cloud_rad_props_get_lw(state, pbuf, cld_abs_od, oldliq, oldice, oldcloud)

   ! Purpose: Compute cloud longwave absorption optical depth

   ! Arguments
   type(physics_state), intent(in)  :: state
   type(physics_buffer_desc),pointer:: pbuf(:)
   real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
   logical, optional,   intent(in)  :: oldliq  ! use old liquid optics
   logical, optional,   intent(in)  :: oldice  ! use old ice optics
   logical, optional,   intent(in)  :: oldcloud  ! use old optics for both (b4b)

   ! Local variables
   integer :: ncol        ! number of columns

   ! rad properties for liquid clouds
   real(r8) :: liq_tau_abs_od(nlwbands,pcols,pver) ! liquid cloud absorption optical depth

   ! rad properties for ice clouds
   real(r8) :: ice_tau_abs_od(nlwbands,pcols,pver) ! ice cloud absorption optical depth
   !-----------------------------------------------------------------------------

   ncol = state%ncol

   cld_abs_od = 0._r8

   if(present(oldcloud))then
      if(oldcloud) then
         call oldcloud_lw(state,pbuf,cld_abs_od,oldwp=.false.)
         return
      endif
   endif

   if(present(oldliq))then
      if(oldliq) then
         call old_liq_get_rad_props_lw(state, pbuf, liq_tau_abs_od, oldliqwp=.false.)
      else
         call liquid_cloud_get_rad_props_lw(state, pbuf, liq_tau_abs_od)
      endif
   else
      call liquid_cloud_get_rad_props_lw(state, pbuf, liq_tau_abs_od)
   endif

   if(present(oldice))then
      if(oldice) then
         call old_ice_get_rad_props_lw(state, pbuf, ice_tau_abs_od, oldicewp=.false.)
      else
         call ice_cloud_get_rad_props_lw(state, pbuf, ice_tau_abs_od)
      endif
   else
      call ice_cloud_get_rad_props_lw(state, pbuf, ice_tau_abs_od)
   endif

   cld_abs_od(:,1:ncol,:) = liq_tau_abs_od(:,1:ncol,:) + ice_tau_abs_od(:,1:ncol,:)

end subroutine cloud_rad_props_get_lw

!==============================================================================

subroutine get_ice_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! asymmetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer :: iciwpth(:,:), dei(:,:)

   ! Get relevant pbuf fields, and interpolate optical properties from
   ! the lookup tables.
   call pbuf_get_field(pbuf, i_iciwp, iciwpth)
   call pbuf_get_field(pbuf, i_dei,   dei)

   call interpolate_ice_optics_sw(state%ncol, iciwpth, dei, tau, tau_w, &
        tau_w_g, tau_w_f)

end subroutine get_ice_optics_sw

!==============================================================================

subroutine get_snow_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! asymmetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer :: icswpth(:,:), des(:,:)

   ! This does the same thing as get_ice_optics_sw, except with a different
   ! water path and effective diameter.
   call pbuf_get_field(pbuf, i_icswp, icswpth)
   call pbuf_get_field(pbuf, i_des,   des)

   call interpolate_ice_optics_sw(state%ncol, icswpth, des, tau, tau_w, &
        tau_w_g, tau_w_f)

end subroutine get_snow_optics_sw

!==============================================================================

subroutine get_grau_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! asymmetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer :: icgrauwpth(:,:), degrau(:,:)

   integer :: i,k
   character(len=*), parameter :: sub = 'get_grau_optics_sw'

   ! This does the same thing as get_ice_optics_sw, except with a different
   ! water path and effective diameter.
   if((i_icgrauwp > 0) .and. (i_degrau > 0)) then

      call pbuf_get_field(pbuf, i_icgrauwp, icgrauwpth)
      call pbuf_get_field(pbuf, i_degrau,   degrau)

      call interpolate_ice_optics_sw(state%ncol, icgrauwpth, degrau, tau, tau_w, &
           tau_w_g, tau_w_f)
      do i = 1, pcols
         do k = 1, pver
            if (tau(idx_sw_diag,i,k).gt.100._r8) then
               write(iulog,*) 'WARNING: SW Graupel Tau > 100  (i,k,icgrauwpth,degrau,tau):'
               write(iulog,*) i,k,icgrauwpth(i,k), degrau(i,k), tau(idx_sw_diag,i,k)
            end if
         enddo
      enddo

   else
      call endrun(sub//': ERROR: Get_grau_optics_sw called when graupel properties not supported')
   end if

end subroutine get_grau_optics_sw

!==============================================================================

subroutine get_liquid_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! asymmetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer, dimension(:,:) :: lamc, pgam, iclwpth
   real(r8), dimension(pcols,pver) :: kext
   integer i,k,swband, ncol

   ncol = state%ncol


   call pbuf_get_field(pbuf, i_lambda,  lamc)
   call pbuf_get_field(pbuf, i_mu,      pgam)
   call pbuf_get_field(pbuf, i_iclwp,   iclwpth)

   do k = 1,pver
      do i = 1,ncol
         if(lamc(i,k) > 0._r8) then ! This seems to be clue from microphysics of no cloud
            call gam_liquid_sw(iclwpth(i,k), lamc(i,k), pgam(i,k), &
                tau(1:nswbands,i,k), tau_w(1:nswbands,i,k), tau_w_g(1:nswbands,i,k), tau_w_f(1:nswbands,i,k))
         else
            tau(1:nswbands,i,k) = 0._r8
            tau_w(1:nswbands,i,k) = 0._r8
            tau_w_g(1:nswbands,i,k) = 0._r8
            tau_w_f(1:nswbands,i,k) = 0._r8
         endif
      enddo
   enddo

end subroutine get_liquid_optics_sw

!==============================================================================

subroutine liquid_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc),pointer  :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   integer :: ncol
   real(r8), pointer, dimension(:,:) :: lamc, pgam, iclwpth

   integer lwband, i, k

   abs_od = 0._r8

   ncol = state%ncol

   call pbuf_get_field(pbuf, i_lambda,  lamc)
   call pbuf_get_field(pbuf, i_mu,      pgam)
   call pbuf_get_field(pbuf, i_iclwp,   iclwpth)

   do k = 1,pver
      do i = 1,ncol
         if(lamc(i,k) > 0._r8) then ! This seems to be the clue for no cloud from microphysics formulation
            call gam_liquid_lw(iclwpth(i,k), lamc(i,k), pgam(i,k), abs_od(1:nlwbands,i,k))
         else
            abs_od(1:nlwbands,i,k) = 0._r8
         endif
      enddo
   enddo

end subroutine liquid_cloud_get_rad_props_lw
!==============================================================================

subroutine snow_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   real(r8), pointer :: icswpth(:,:), des(:,:)

   ! This does the same thing as ice_cloud_get_rad_props_lw, except with a
   ! different water path and effective diameter.
   call pbuf_get_field(pbuf, i_icswp, icswpth)
   call pbuf_get_field(pbuf, i_des,   des)

   call interpolate_ice_optics_lw(state%ncol,icswpth, des, abs_od)

end subroutine snow_cloud_get_rad_props_lw


!==============================================================================

subroutine grau_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   real(r8), pointer :: icgrauwpth(:,:), degrau(:,:)
   character(len=*), parameter :: sub = 'grau_cloud_get_rad_props_lw'

   ! This does the same thing as ice_cloud_get_rad_props_lw, except with a
   ! different water path and effective diameter.
   if((i_icgrauwp > 0) .and. (i_degrau > 0)) then
      call pbuf_get_field(pbuf, i_icgrauwp, icgrauwpth)
      call pbuf_get_field(pbuf, i_degrau,   degrau)

      call interpolate_ice_optics_lw(state%ncol,icgrauwpth, degrau, abs_od)
   else
      call endrun(sub//': ERROR: Grau_cloud_get_rad_props_lw called when graupel &
           &properties not supported')
   end if

end subroutine grau_cloud_get_rad_props_lw

!==============================================================================

subroutine ice_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)     :: state
   type(physics_buffer_desc), pointer  :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   real(r8), pointer :: iciwpth(:,:), dei(:,:)

   ! Get relevant pbuf fields, and interpolate optical properties from
   ! the lookup tables.
   call pbuf_get_field(pbuf, i_iciwp, iciwpth)
   call pbuf_get_field(pbuf, i_dei,   dei)

   call interpolate_ice_optics_lw(state%ncol,iciwpth, dei, abs_od)

end subroutine ice_cloud_get_rad_props_lw

!==============================================================================
! Private methods
!==============================================================================

subroutine interpolate_ice_optics_sw(ncol, iciwpth, dei, tau, tau_w, &
     tau_w_g, tau_w_f)

  integer, intent(in) :: ncol
  real(r8), intent(in) :: iciwpth(pcols,pver)
  real(r8), intent(in) :: dei(pcols,pver)

  real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
  real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
  real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! asymmetry parameter * tau * w
  real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

  type(interp_type) :: dei_wgts

  integer :: i, k, swband
  real(r8) :: ext(nswbands), ssa(nswbands), asm(nswbands)

  do k = 1,pver
     do i = 1,ncol
        if( iciwpth(i,k) < tiny .or. dei(i,k) == 0._r8) then
           ! if ice water path is too small, OD := 0
           tau    (:,i,k) = 0._r8
           tau_w  (:,i,k) = 0._r8
           tau_w_g(:,i,k) = 0._r8
           tau_w_f(:,i,k) = 0._r8
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call lininterp_init(g_d_eff, n_g_d, dei(i:i,k), 1, &
                extrap_method_bndry, dei_wgts)
           ! interpolate into grid and extract radiative properties
           do swband = 1, nswbands
              call lininterp(ext_sw_ice(:,swband), n_g_d, &
                   ext(swband:swband), 1, dei_wgts)
              call lininterp(ssa_sw_ice(:,swband), n_g_d, &
                   ssa(swband:swband), 1, dei_wgts)
              call lininterp(asm_sw_ice(:,swband), n_g_d, &
                   asm(swband:swband), 1, dei_wgts)
           end do
           tau    (:,i,k) = iciwpth(i,k) * ext
           tau_w  (:,i,k) = tau(:,i,k) * ssa
           tau_w_g(:,i,k) = tau_w(:,i,k) * asm
           tau_w_f(:,i,k) = tau_w_g(:,i,k) * asm
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine interpolate_ice_optics_sw

!==============================================================================

subroutine interpolate_ice_optics_lw(ncol, iciwpth, dei, abs_od)

  integer, intent(in) :: ncol
  real(r8), intent(in) :: iciwpth(pcols,pver)
  real(r8), intent(in) :: dei(pcols,pver)

  real(r8),intent(out) :: abs_od(nlwbands,pcols,pver)

  type(interp_type) :: dei_wgts

  integer :: i, k, lwband
  real(r8) :: absor(nlwbands)

  do k = 1,pver
     do i = 1,ncol
        ! if ice water path is too small, OD := 0
        if( iciwpth(i,k) < tiny .or. dei(i,k) == 0._r8) then
           abs_od (:,i,k) = 0._r8
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call lininterp_init(g_d_eff, n_g_d, dei(i:i,k), 1, &
                extrap_method_bndry, dei_wgts)
           ! interpolate into grid and extract radiative properties
           do lwband = 1, nlwbands
              call lininterp(abs_lw_ice(:,lwband), n_g_d, &
                   absor(lwband:lwband), 1, dei_wgts)
           enddo
           abs_od(:,i,k) = iciwpth(i,k) * absor
           where(abs_od(:,i,k) > 50.0_r8) abs_od(:,i,k) = 50.0_r8
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine interpolate_ice_optics_lw

!==============================================================================

subroutine gam_liquid_lw(clwptn, lamc, pgam, abs_od)
  real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  real(r8), intent(out) :: abs_od(1:nlwbands)

  integer :: lwband ! sw band index

  type(interp_type) :: mu_wgts
  type(interp_type) :: lambda_wgts

  if (clwptn < tiny) then
    abs_od = 0._r8
    return
  endif

  call get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts)

  do lwband = 1, nlwbands
     call lininterp(abs_lw_liq(:,:,lwband), nmu, nlambda, &
          abs_od(lwband:lwband), 1, mu_wgts, lambda_wgts)
  enddo

  abs_od = clwptn * abs_od

  call lininterp_finish(mu_wgts)
  call lininterp_finish(lambda_wgts)

end subroutine gam_liquid_lw

!==============================================================================

subroutine gam_liquid_sw(clwptn, lamc, pgam, tau, tau_w, tau_w_g, tau_w_f)
  real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  real(r8), intent(out) :: tau(1:nswbands), tau_w(1:nswbands), tau_w_f(1:nswbands), tau_w_g(1:nswbands)

  integer :: swband ! sw band index

  real(r8) :: ext(nswbands), ssa(nswbands), asm(nswbands)

  type(interp_type) :: mu_wgts
  type(interp_type) :: lambda_wgts

  if (clwptn < tiny) then
    tau = 0._r8
    tau_w = 0._r8
    tau_w_g = 0._r8
    tau_w_f = 0._r8
    return
  endif

  call get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts)

  do swband = 1, nswbands
     call lininterp(ext_sw_liq(:,:,swband), nmu, nlambda, &
          ext(swband:swband), 1, mu_wgts, lambda_wgts)
     call lininterp(ssa_sw_liq(:,:,swband), nmu, nlambda, &
          ssa(swband:swband), 1, mu_wgts, lambda_wgts)
     call lininterp(asm_sw_liq(:,:,swband), nmu, nlambda, &
          asm(swband:swband), 1, mu_wgts, lambda_wgts)
  enddo

  ! compute radiative properties
  tau = clwptn * ext
  tau_w = tau * ssa
  tau_w_g = tau_w * asm
  tau_w_f = tau_w_g * asm

  call lininterp_finish(mu_wgts)
  call lininterp_finish(lambda_wgts)

end subroutine gam_liquid_sw

!==============================================================================

subroutine get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts)
  use radiation_utils, only: get_mu_lambda_weights_ccpp
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  ! Output interpolation weights. Caller is responsible for freeing these.
  type(interp_type), intent(out) :: mu_wgts
  type(interp_type), intent(out) :: lambda_wgts

  character(len=512) :: errmsg
  integer            :: errflg

  call get_mu_lambda_weights_ccpp(nmu, nlambda, g_mu, g_lambda, lamc, pgam, mu_wgts, &
          lambda_wgts, errmsg, errflg)
  if (errflg /= 0) then
     call endrun('get_mu_lambda_weights: ERROR message: '//errmsg)
  end if

end subroutine get_mu_lambda_weights

!==============================================================================

end module cloud_rad_props
