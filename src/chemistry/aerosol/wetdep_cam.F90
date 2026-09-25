module wetdep_cam

!-----------------------------------------------------------------------
!
! CAM host interface for the portable wetdep module
! Calls portable wetdepa.
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, pver
   use physconst,    only: gravit, rair, tmelt
   use phys_control, only: cam_physpkg_is

   use wetdep,       only: clddiag

   implicit none
   private

   public :: wetdep_inputs_t
   public :: wetdep_init
   public :: wetdep_inputs_set

   type wetdep_inputs_t
      real(r8), pointer :: cldt(:,:) => null()  ! cloud fraction
      real(r8), pointer :: prain(:,:) => null()
      real(r8), pointer :: bergso(:,:) => null()
      real(r8), pointer :: evapr(:,:) => null()
      real(r8) :: cldcu(pcols,pver)     ! convective cloud fraction, currently empty
      real(r8) :: evapc(pcols,pver)     ! Evaporation rate of convective precipitation
      real(r8) :: cmfdqr(pcols,pver)    ! convective production of rain
      real(r8) :: conicw(pcols,pver)    ! convective in-cloud water
      real(r8) :: totcond(pcols, pver)  ! total condensate
      real(r8) :: cldv(pcols,pver)      ! cloudy volume undergoing wet chem and scavenging
      real(r8) :: cldvcu(pcols,pver)    ! Convective precipitation area at the top interface of current layer
      real(r8) :: cldvst(pcols,pver)    ! Stratiform precipitation area at the top interface of current layer
   end type wetdep_inputs_t

   integer :: cld_idx             = 0
   integer :: prain_idx           = 0
   integer :: bergso_idx          = 0
   integer :: nevapr_idx          = 0

   integer :: icwmrdp_idx         = 0
   integer :: icwmrsh_idx         = 0
   integer :: rprddp_idx          = 0
   integer :: rprdsh_idx          = 0
   integer :: sh_frac_idx         = 0
   integer :: dp_frac_idx         = 0
   integer :: nevapr_shcu_idx     = 0
   integer :: nevapr_dpcu_idx     = 0
   integer :: ixcldice, ixcldliq

contains

   subroutine wetdep_init()
      use physics_buffer, only: pbuf_get_index
      use constituents,   only: cnst_get_ind

      integer :: ierr

      cld_idx             = pbuf_get_index('CLD')
      prain_idx           = pbuf_get_index('PRAIN')
      bergso_idx          = pbuf_get_index('BERGSO', errcode=ierr )
      nevapr_idx          = pbuf_get_index('NEVAPR')

      icwmrdp_idx         = pbuf_get_index('ICWMRDP')
      rprddp_idx          = pbuf_get_index('RPRDDP')
      icwmrsh_idx         = pbuf_get_index('ICWMRSH')
      rprdsh_idx          = pbuf_get_index('RPRDSH')
      sh_frac_idx         = pbuf_get_index('SH_FRAC' )
      dp_frac_idx         = pbuf_get_index('DP_FRAC')
      nevapr_shcu_idx     = pbuf_get_index('NEVAPR_SHCU')
      nevapr_dpcu_idx     = pbuf_get_index('NEVAPR_DPCU')

      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('CLDLIQ', ixcldliq)

   end subroutine wetdep_init

!==============================================================================
! gathers up the inputs needed for the wetdepa routines
!==============================================================================
   subroutine wetdep_inputs_set( state, pbuf, inputs )
      use physics_types,  only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

      ! args

      type(physics_state),  intent(in )  :: state           !! physics state
      type(physics_buffer_desc), pointer :: pbuf(:)         !! physics buffer
      type(wetdep_inputs_t), intent(out) :: inputs          !! collection of wetdepa inputs

      ! local vars

      real(r8), pointer :: icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
      real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
      real(r8), pointer :: icwmrsh(:,:)    ! in cloud water mixing ratio, deep convection
      real(r8), pointer :: rprdsh(:,:)     ! rain production, deep convection
      real(r8), pointer :: sh_frac(:,:)    ! Shallow convective cloud fraction
      real(r8), pointer :: dp_frac(:,:)    ! Deep convective cloud fraction
      real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
      real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.

      real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume
      real(r8) :: cldst(pcols,pver)        ! Stratiform cloud fraction

      integer :: itim, ncol

      ncol = state%ncol
      itim = pbuf_old_tim_idx()

      call pbuf_get_field(pbuf, cld_idx,         inputs%cldt, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, prain_idx,       inputs%prain   )
      call pbuf_get_field(pbuf, nevapr_idx,      inputs%evapr   )
      call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp )
      call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh )
      call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
      call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
      call pbuf_get_field(pbuf, sh_frac_idx,     sh_frac )
      call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac )
      call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
      call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )

      if (bergso_idx>0) then
         call pbuf_get_field(pbuf, bergso_idx, inputs%bergso )
      else
         if (.not. associated(inputs%bergso)) then
            allocate(inputs%bergso(pcols,pver))
            inputs%bergso(:,:) = 0.0_r8
         endif
      endif

      inputs%cldcu(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)
      cldst(:ncol,:)          = inputs%cldt(:ncol,:) - inputs%cldcu(:ncol,:)       ! Stratiform cloud fraction
      inputs%evapc(:ncol,:)  = evapcsh(:ncol,:) + evapcdp(:ncol,:)
      inputs%cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)

      ! sum deep and shallow convection contributions
      if (cam_physpkg_is('cam5') .or. cam_physpkg_is('cam6')) then
         ! Dec.29.2009. Sungsu
         inputs%conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                                  max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
      else
         inputs%conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
      end if

      inputs%totcond(:ncol,:) = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)

      call clddiag( state%t,     state%pmid,   state%pdel,   inputs%cmfdqr, inputs%evapc, &
                   inputs%cldt,  inputs%cldcu,       cldst,  inputs%evapr, &
                   inputs%prain, inputs%cldv, inputs%cldvcu, inputs%cldvst,       rainmr, &
                   state%ncol, pver, gravit, tmelt, rair )

   end subroutine wetdep_inputs_set

end module wetdep_cam
