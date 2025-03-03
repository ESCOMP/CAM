! CAM interface for Hack shallow moist convection
module hk_conv
   use shr_kind_mod,     only: r8 => shr_kind_r8
   use cam_logfile,      only: iulog
   use spmd_utils,       only: masterproc
   use cam_abortutils,   only: endrun

   implicit none
   private
   save

   public :: hkconv_readnl
   public :: mfinti
   public :: cmfmca_cam

   !
   ! Private data used for Hack shallow convection
   !
   real(r8), parameter :: unset_r8 = huge(1.0_r8)

   ! Namelist variables
   real(r8) :: hkconv_c0 = unset_r8
   real(r8) :: hkconv_cmftau = unset_r8

   real(r8) :: cmftau      ! characteristic adjustment time scale set from namelist input hkconv_cmftau
   real(r8) :: c0          ! rain water autoconversion coefficient set from namelist input hkconv_c0

contains
   subroutine hkconv_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'hkconv_readnl'

      namelist /hkconv_nl/ hkconv_cmftau, hkconv_c0
      !-----------------------------------------------------------------------------

      if (masterproc) then
         open( newunit=unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'hkconv_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, hkconv_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)

         ! set local variables
         cmftau = hkconv_cmftau
         c0     = hkconv_c0

      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(cmftau,            1, mpir8,  0, mpicom)
      call mpibcast(c0,                1, mpir8,  0, mpicom)
#endif

   end subroutine hkconv_readnl

   subroutine mfinti (rair    ,cpair   ,gravit  ,latvap  ,rhowtr,pref_edge )
      use spmd_utils, only: masterproc
      use hack_convect_shallow,     only: hack_convect_shallow_init
      use ppgrid,     only: pver, pcols

      real(r8), intent(in) :: rair              ! gas constant for dry air
      real(r8), intent(in) :: cpair             ! specific heat of dry air
      real(r8), intent(in) :: gravit            ! acceleration due to gravity
      real(r8), intent(in) :: latvap            ! latent heat of vaporization
      real(r8), intent(in) :: rhowtr            ! density of liquid water (STP)
      real(r8), intent(in)  :: pref_edge(:) ! reference pressures at interface [Pa]

      character(len=512)   :: errmsg
      integer              :: errflg

      ! dummy parameters output from CCPPized scheme.
      ! the dimensions do not matter as the whole array is assigned zero
      logical              :: dummy_use_shfrc
      real(r8)             :: dummy_shfrc_out(pcols, pver)
      integer              :: dummy_top_lev_out

      ! Initialize free parameters for moist convective mass flux procedure
      ! cmftau - characteristic adjustment time scale
      ! c0     - rain water autoconversion coeff (1/m)
      !
      ! call CCPPized subroutine
      call hack_convect_shallow_init( &
         pver = pver, &
         amIRoot = masterproc, &
         iulog = iulog, &
         cmftau_in = cmftau, &
         c0_in = c0, &
         rair = rair, &
         cpair = cpair, &
         gravit = gravit, &
         latvap = latvap, &
         rhoh2o_in = rhowtr, &
         pref_edge = pref_edge, &
         use_shfrc = dummy_use_shfrc, &
         shfrc = dummy_shfrc_out, &
         top_lev = dummy_top_lev_out, &
         errmsg = errmsg, &
         errflg = errflg &
      )
   end subroutine mfinti

   subroutine cmfmca_cam(lchnk   ,ncol    , &
                  nstep   ,ztodt     ,pmid    ,pdel    , &
                  rpdel   ,zm      ,tpert   ,qpert   ,phis    , &
                  pblh    ,t       ,q       ,cmfdt   ,dq      , &
                  cmfmc   ,cmfdqr  ,cmfsl   ,cmflq   ,precc   , &
                  qc      ,cnt     ,cnb     ,icwmr   ,rliq    , &
                  pmiddry ,pdeldry ,rpdeldry)
      use ppgrid,                    only: pcols, pver, pverp
      use constituents,              only: pcnst
      use ccpp_constituent_prop_mod, only: ccpp_const_props

      use hack_convect_shallow, only: hack_convect_shallow_run
      !
      ! Input arguments
      !
      integer, intent(in) :: lchnk                ! chunk identifier
      integer, intent(in) :: ncol                 ! number of atmospheric columns
      integer, intent(in) :: nstep                ! current time step index

      real(r8), intent(in) :: ztodt               ! physics timestep (seconds)
      real(r8), intent(in) :: pmid(pcols,pver)    ! pressure
      real(r8), intent(in) :: pdel(pcols,pver)    ! delta-p
      real(r8), intent(in) :: pmiddry(pcols,pver)    ! dry air pressure
      real(r8), intent(in) :: pdeldry(pcols,pver)    ! dry air delta-p
      real(r8), intent(in) :: rpdel(pcols,pver)   ! 1./pdel
      real(r8), intent(in) :: rpdeldry(pcols,pver)   ! 1./pdeldry
      real(r8), intent(in) :: zm(pcols,pver)      ! height abv sfc at midpoints
      real(r8), intent(in) :: tpert(pcols)        ! PBL perturbation theta
      real(r8), intent(in) :: qpert(pcols)        ! PBL perturbation specific humidity
      real(r8), intent(in) :: phis(pcols)         ! surface geopotential
      real(r8), intent(in) :: pblh(pcols)         ! PBL height (provided by PBL routine)
      real(r8), intent(in) :: t(pcols,pver)       ! temperature (t bar)
      real(r8), intent(in) :: q(pcols,pver,pcnst) ! specific humidity (sh bar)
      !
      ! Output arguments
      !
      real(r8), intent(out) :: cmfdt(pcols,pver)   ! dt/dt due to moist convection
      real(r8), intent(out) :: cmfmc(pcols,pverp)  ! moist convection cloud mass flux
      real(r8), intent(out) :: cmfdqr(pcols,pver)  ! dq/dt due to convective rainout
      real(r8), intent(out) :: cmfsl(pcols,pver )  ! convective lw static energy flux
      real(r8), intent(out) :: cmflq(pcols,pver )  ! convective total water flux
      real(r8), intent(out) :: precc(pcols)        ! convective precipitation rate
      ! JJH mod to explicitly export cloud water
      real(r8), intent(out) :: qc(pcols,pver)      ! dq/dt due to export of cloud water
      real(r8), intent(out) :: cnt(pcols)          ! top level of convective activity
      real(r8), intent(out) :: cnb(pcols)          ! bottom level of convective activity
      real(r8), intent(out) :: dq(pcols,pver,pcnst) ! constituent tendencies
      real(r8), intent(out) :: icwmr(pcols,pver)
      real(r8), intent(out) :: rliq(pcols)

      ! local variables
      character(len=512) :: errmsg
      integer :: errflg
      character(len=64) :: dummy_scheme_name  ! dummy scheme name for CCPP-ized scheme
      real(r8) :: flx_cnd(pcols) ! dummy flx_cnd for CCPP-ized scheme (actual check_energy flux computed in convect_shallow)

      ! integer output arguments for CCPP-ized scheme to be later converted to real
      integer :: cnt_out(ncol)
      integer :: cnb_out(ncol)

      errmsg = ''
      errflg = 0

      ! zero output arguments to pcols
      cmfdt(:,:) = 0.0_r8
      cmfmc(:,:) = 0.0_r8
      cmfdqr(:,:) = 0.0_r8
      cmfsl(:,:) = 0.0_r8
      cmflq(:,:) = 0.0_r8
      precc(:) = 0.0_r8
      qc(:,:) = 0.0_r8
      cnt(:) = 0.0_r8
      cnb(:) = 0.0_r8
      dq(:,:,:) = 0.0_r8
      icwmr(:,:) = 0.0_r8
      rliq(:) = 0.0_r8

      ! Call the CCPPized subroutine with subsetting to ncol
      call hack_convect_shallow_run( &
         ncol = ncol, &
         pver = pver, &
         pcnst = pcnst, &
         iulog = iulog, &
         const_props = ccpp_const_props, &
         ztodt = ztodt, &
         pmid = pmid(:ncol,:), &
         pmiddry = pmiddry(:ncol,:), &
         pdel = pdel(:ncol,:), &
         pdeldry = pdeldry(:ncol,:), &
         rpdel = rpdel(:ncol,:), &
         rpdeldry = rpdeldry(:ncol,:), &
         zm = zm(:ncol,:), &
         ! tpert is always zero
         qpert_in = qpert(:ncol), &
         phis = phis(:ncol), &
         pblh = pblh(:ncol), &
         t = t(:ncol,:), &
         q = q(:ncol,:,:), &
         dq = dq(:ncol,:,:), &
         qc_sh = qc(:ncol,:), &
         cmfdt = cmfdt(:ncol,:), &
         cmfmc_sh = cmfmc(:ncol,:), &
         cmfdqr = cmfdqr(:ncol,:), &
         cmfsl = cmfsl(:ncol,:), &
         cmflq = cmflq(:ncol,:), &
         precc = precc(:ncol), &
         cnt_sh = cnt_out(:ncol), &
         cnb_sh = cnb_out(:ncol), &
         icwmr = icwmr(:ncol,:), &
         rliq_sh = rliq(:ncol), &
         scheme_name = dummy_scheme_name, &
         flx_cnd = flx_cnd(:ncol), &
         errmsg = errmsg, &
         errflg = errflg &
      )

      ! convert back to real for diagnostics
      cnt(:ncol) = real(cnt_out(:ncol), r8)
      cnb(:ncol) = real(cnb_out(:ncol), r8)
   end subroutine cmfmca_cam

end module hk_conv
