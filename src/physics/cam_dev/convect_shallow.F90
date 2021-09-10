   module convect_shallow

   !----------------------------------------------- !
   ! Purpose:                                       !
   !                                                !
   ! CAM interface to the shallow convection scheme !
   !                                                !
   ! Author: D.B. Coleman                           !
   !         Sungsu Park. Jan. 2010.                !
   ! STUB version: A. Herrington. Sep. 2021 !
   !----------------------------------------------- !

   use shr_kind_mod,      only : r8=>shr_kind_r8
   use ppgrid,            only : pver, pcols, pverp
   use cam_history,       only : outfld, addfld, horiz_only
!+++ARH
!   use cam_logfile,       only : iulog
!---ARH
   use phys_control,      only : phys_getopts
   implicit none
   private                 
   save

   public :: &
             convect_shallow_register,       & ! Register fields in physics buffer
             convect_shallow_init,           & ! Initialize shallow module
             convect_shallow_tend,           & ! Return tendencies
             convect_shallow_use_shfrc         ! 

   character(len=16) :: shallow_scheme      ! Default set in phys_control.F90, use namelist to change

   ! Physics buffer indices 
   integer    ::     icwmrsh_idx    = 0  
   integer    ::      rprdsh_idx    = 0 
   integer    ::     rprdtot_idx    = 0 
   integer    ::      cldtop_idx    = 0 
   integer    ::      cldbot_idx    = 0 
   integer    :: nevapr_shcu_idx    = 0
   integer    ::      rprddp_idx    = 0
   integer    ::     prec_sh_idx    = 0
   integer    ::     snow_sh_idx    = 0
   integer    ::    cmfmc_sh_idx    = 0

   contains

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine convect_shallow_register

  !-------------------------------------------------- !
  ! Purpose : Register fields with the physics buffer !
  !-------------------------------------------------- !
!+++ARH
!  use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls
!---ARH
  use physics_buffer, only : pbuf_add_field, dtype_r8

  call phys_getopts( shallow_scheme_out = shallow_scheme )

  call pbuf_add_field('ICWMRSH',    'physpkg' ,dtype_r8,(/pcols,pver/),       icwmrsh_idx )
  call pbuf_add_field('RPRDSH',     'physpkg' ,dtype_r8,(/pcols,pver/),       rprdsh_idx )
  call pbuf_add_field('RPRDTOT',    'physpkg' ,dtype_r8,(/pcols,pver/),       rprdtot_idx )
  call pbuf_add_field('CLDTOP',     'physpkg' ,dtype_r8,(/pcols,1/),          cldtop_idx )
  call pbuf_add_field('CLDBOT',     'physpkg' ,dtype_r8,(/pcols,1/),          cldbot_idx )
  call pbuf_add_field('NEVAPR_SHCU','physpkg' ,dtype_r8,(/pcols,pver/),       nevapr_shcu_idx )
  call pbuf_add_field('PREC_SH',    'physpkg' ,dtype_r8,(/pcols/),            prec_sh_idx )
  call pbuf_add_field('SNOW_SH',    'physpkg' ,dtype_r8,(/pcols/),            snow_sh_idx )
  ! Updraft mass flux by shallow convection [ kg/s/m2 ]
  call pbuf_add_field('CMFMC_SH',   'physpkg' ,dtype_r8,(/pcols,pverp/),      cmfmc_sh_idx )

  end subroutine convect_shallow_register

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !


  subroutine convect_shallow_init(pref_edge, pbuf2d)

  !------------------------------------------------------------------------------- !
  ! Purpose : Declare output fields, and initialize variables needed by convection !
  !------------------------------------------------------------------------------- !
!+++ARH
!  use cam_history,       only : addfld, add_default
!---ARH
  use cam_history,       only : addfld
  use pmgrid,            only : plev, plevp
!+++ARH
!  use spmd_utils,        only : masterproc
!  use cam_abortutils,    only : endrun
!---ARH
  use phys_control,      only : cam_physpkg_is
  
  use physics_buffer,    only : pbuf_get_index, physics_buffer_desc, pbuf_set_field

  real(r8),                  intent(in) :: pref_edge(plevp)  ! Reference pressures at interfaces
  type(physics_buffer_desc), pointer    :: pbuf2d(:,:)

  call addfld( 'CMFMC',      (/ 'ilev' /), 'A', 'kg/m2/s',  'Moist convection (deep+shallow) mass flux'                 )
  call addfld( 'CLDTOP',     horiz_only,   'I', '1',        'Vertical index of cloud top'                               )
  call addfld( 'CLDBOT',     horiz_only,   'I', '1',        'Vertical index of cloud base'                              )
  call addfld( 'PCLDTOP',    horiz_only,   'A', '1',        'Pressure of cloud top'                                     )
  call addfld( 'PCLDBOT',    horiz_only,   'A', '1',        'Pressure of cloud base'                                    )

  rprddp_idx   = pbuf_get_index('RPRDDP')

  end subroutine convect_shallow_init

!==================================================================================================
  function convect_shallow_use_shfrc()
  !-------------------------------------------------------------- !
  ! Return true if cloud fraction should use shallow convection   !
  !          calculated convective clouds.                        !
  !-------------------------------------------------------------- !
     implicit none
     logical :: convect_shallow_use_shfrc     ! Return value

     if (shallow_scheme .eq. 'UW' .or. shallow_scheme .eq. 'UNICON') then
          convect_shallow_use_shfrc = .true.
     else
          convect_shallow_use_shfrc = .false.
     endif

     return

  end function convect_shallow_use_shfrc
  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine convect_shallow_tend( ztodt  , cmfmc   , &
                                   qc     , qc2     , rliq     , rliq2    , & 
                                   state  , ptend_all, pbuf, cam_in)
!+++ARH
!   use physics_buffer,  only : physics_buffer_desc, pbuf_get_field, pbuf_set_field, pbuf_old_tim_idx
!---ARH
   use physics_buffer,  only : physics_buffer_desc, pbuf_get_field, pbuf_set_field
   use cam_history,     only : outfld
   use physics_types,   only : physics_state, physics_ptend
   use camsrfexch,      only : cam_in_t
   implicit none

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state), intent(in)    :: state                           ! Physics state variables
   real(r8),            intent(in)    :: ztodt                           ! 2 delta-t  [ s ]

   type(physics_ptend), intent(out)   :: ptend_all                       ! Indivdual parameterization tendencies
   real(r8),            intent(out)   :: rliq2(pcols)                    ! Vertically-integrated reserved cloud condensate [ m/s ]
   real(r8),            intent(out)   :: qc2(pcols,pver)                 ! Same as qc but only from shallow convection scheme

   

   real(r8),            intent(inout) :: cmfmc(pcols,pverp)    ! Moist deep + shallow convection cloud mass flux [ kg/s/m2 ]
   real(r8),            intent(inout) :: qc(pcols,pver)        ! dq/dt due to export of cloud water into environment by shallow
                                                               ! and deep convection [ kg/kg/s ]
   real(r8),            intent(inout) :: rliq(pcols)           ! Vertical integral of qc [ m/s ]

   type(cam_in_t),      intent(in)    :: cam_in

!+++ARH
!   integer  :: i, k, m
!---ARH
   integer  :: i
   integer  :: lchnk                                           ! Chunk identifier
   integer  :: ncol                                            ! Number of atmospheric columns

   real(r8),  pointer   :: precc(:)                            ! Shallow convective precipitation (rain+snow) rate at surface [ m/s ]
   real(r8),  pointer   :: snow(:)                             ! Shallow convective snow rate at surface [ m/s ]

   real(r8) :: cnt2(pcols)                                     ! Top level of shallow convective activity
   real(r8) :: cnb2(pcols)                                     ! Bottom level of convective activity
   real(r8) :: pcnt(pcols)                                     ! Top    pressure level of shallow + deep convective activity
   real(r8) :: pcnb(pcols)                                     ! Bottom pressure level of shallow + deep convective activity
!+++ARH
!   integer itim_old, ifld
!---ARH
   real(r8), pointer, dimension(:,:) :: icwmr                  ! In cloud water + ice mixing ratio
   real(r8), pointer, dimension(:,:) :: rprddp                 ! dq/dt due to deep convective rainout
   real(r8), pointer, dimension(:,:) :: rprdsh                 ! dq/dt due to deep and shallow convective rainout
   real(r8), pointer, dimension(:,:) :: evapcsh                ! Evaporation of shallow convective precipitation >= 0.
   real(r8), pointer, dimension(:)   :: cnt
   real(r8), pointer, dimension(:)   :: cnb
   real(r8), pointer, dimension(:,:) :: cmfmc2                 ! (pcols,pverp) Updraft mass flux by shallow convection [ kg/s/m2 ]

   ! ----------------------- !
   ! Main Computation Begins ! 
   ! ----------------------- !

   lchnk = state%lchnk
   ncol  = state%ncol

   call pbuf_get_field(pbuf, icwmrsh_idx,     icwmr)

   call pbuf_get_field(pbuf, rprddp_idx,      rprddp )

   call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh )

   call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh  )

   call pbuf_get_field(pbuf, cldtop_idx,      cnt )

   call pbuf_get_field(pbuf, cldbot_idx,      cnb )

   call pbuf_get_field(pbuf, prec_sh_idx,   precc )

   call pbuf_get_field(pbuf, snow_sh_idx,    snow )

   call pbuf_get_field(pbuf, cmfmc_sh_idx,  cmfmc2)

   ! Initialization

   cmfmc2      = 0._r8
   rprdsh      = 0._r8
   precc       = 0._r8
   icwmr       = 0._r8
   rliq2       = 0._r8
   qc2         = 0._r8
   cnt2        = pver
   cnb2        = 1._r8
   evapcsh     = 0._r8
   snow        = 0._r8

   ! ------------------------------------------------------------------------------ !
   ! Merge shallow convection output with prior results from deep convection scheme !
   ! ------------------------------------------------------------------------------ !

   ! ----------------------------------------------------------------------- !
   ! Combine cumulus updraft mass flux : 'cmfmc2'(shallow) + 'cmfmc'(deep)   !
   ! ----------------------------------------------------------------------- !

   cmfmc(:ncol,:) = cmfmc(:ncol,:) + cmfmc2(:ncol,:)

   ! -------------------------------------------------------------- !
   ! 'cnt2' & 'cnb2' are from shallow, 'cnt' & 'cnb' are from deep  !
   ! 'cnt2' & 'cnb2' are the interface indices of cloud top & base: ! 
   !        cnt2 = float(kpen)                                      !
   !        cnb2 = float(krel - 1)                                  !
   ! Note that indices decreases with height.                       !
   ! -------------------------------------------------------------- !

   do i = 1, ncol
      if( cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
      if( cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
      pcnt(i) = state%pmid(i,int(cnt(i)))
      pcnb(i) = state%pmid(i,int(cnb(i)))     
   end do
   
   ! ----------------------------------------------- !
   ! This quantity was previously known as CMFDQR.   !
   ! Now CMFDQR is the shallow rain production only. !
   ! ----------------------------------------------- !

   call pbuf_set_field(pbuf, rprdtot_idx, rprdsh(:ncol,:pver) + rprddp(:ncol,:pver), start=(/1,1/), kount=(/ncol,pver/))
 
   ! ----------------------------------------------------------------------- ! 
   ! Add shallow reserved cloud condensate to deep reserved cloud condensate !
   !     qc [ kg/kg/s] , rliq [ m/s ]                                        !
   ! ----------------------------------------------------------------------- !

   qc(:ncol,:pver) = qc(:ncol,:pver) + qc2(:ncol,:pver)
   rliq(:ncol)     = rliq(:ncol) + rliq2(:ncol)    

   ! ---------------------------------------------------------------------------- !
   ! Output new partition of cloud condensate variables, as well as precipitation !
   ! ---------------------------------------------------------------------------- ! 

   call outfld( 'CMFMC'  , cmfmc                     , pcols   , lchnk )
   call outfld( 'CLDTOP' , cnt                       , pcols   , lchnk )
   call outfld( 'CLDBOT' , cnb                       , pcols   , lchnk )
   call outfld( 'PCLDTOP', pcnt                      , pcols   , lchnk )
   call outfld( 'PCLDBOT', pcnb                      , pcols   , lchnk )  

  end subroutine convect_shallow_tend

  end module convect_shallow
