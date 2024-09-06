module iop
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: iop
!
! !DESCRIPTION:
! iop specific routines
!
! !USES:
!
  use cam_abortutils,   only: endrun
  use constituents,     only: pcnst
  use eul_control_mod,  only: eul_nsplit
  use pmgrid,           only: beglat,endlat,plon,plev
  use shr_kind_mod,     only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none


  private

  real(r8), allocatable,target :: dqfx3sav(:,:,:,:)
  real(r8), allocatable,target :: t2sav(:,:,:)
  real(r8), allocatable,target :: fusav(:,:,:)
  real(r8), allocatable,target :: fvsav(:,:,:)
  real(r8), allocatable,target :: divq3dsav(:,:,:,:)
  real(r8), allocatable,target :: divt3dsav(:,:,:)
  real(r8), allocatable,target :: divu3dsav(:,:,:)
  real(r8), allocatable,target :: divv3dsav(:,:,:)
  real(r8), allocatable,target :: betasav(:)

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_iop_fields
  public :: iop_update_prognostics
! !PUBLIC DATA:
  public betasav, &
         dqfx3sav, divq3dsav, divt3dsav,divu3dsav,divv3dsav,t2sav,fusav,fvsav

!
! !REVISION HISTORY:
! Created by John Truesdale
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!-----------------------------------------------------------------------

contains
   subroutine init_iop_fields()
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
   implicit none
   character(len=*), parameter ::  sub = "init_iop_fields"
!-----------------------------------------------------------------------
   if (eul_nsplit>1) then
      call endrun('iop module cannot be used with eul_nsplit>1')
   endif

   if(.not.allocated(betasav)) then
      allocate (betasav(beglat:endlat))
      betasav(:)=0._r8
   endif

   if(.not.allocated(dqfx3sav)) then
      allocate (dqfx3sav(plon,plev,pcnst,beglat:endlat))
      dqfx3sav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divq3dsav)) then
      allocate (divq3dsav(plon,plev,pcnst,beglat:endlat))
      divq3dsav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divt3dsav)) then
      allocate (divt3dsav(plon,plev,beglat:endlat))
      divt3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(divu3dsav)) then
      allocate (divu3dsav(plon,plev,beglat:endlat))
      divu3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(divv3dsav)) then
      allocate (divv3dsav(plon,plev,beglat:endlat))
      divv3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(t2sav)) then
      allocate (t2sav(plon,plev,beglat:endlat))  ! temp tendency
      t2sav(:,:,:)=0._r8
   endif
   if(.not.allocated(fusav)) then
      allocate (fusav(plon,plev,beglat:endlat))  ! U wind tendency
      fusav(:,:,:)=0._r8
   endif
   if(.not.allocated(fvsav)) then
      allocate (fvsav(plon,plev,beglat:endlat))  ! v wind tendency
      fvsav(:,:,:)=0._r8
   endif
  end subroutine init_iop_fields

  subroutine iop_update_prognostics(timelevel,ps,t3,u3,v3,q3)
!------------------------------------------------------------------------------
! Copy IOP forcing fields into prognostics which for Eulerian is just PS
!------------------------------------------------------------------------------
   use scamMod,             only: tobs,uobs,vobs,qobs,psobs
   implicit none

   !-----------------------------------------------------------------------

   integer,  intent(in)   :: timelevel
   real(r8), optional, intent(inout) :: q3(:,:,:,:,:)
   real(r8), optional, intent(inout) :: u3(:,:,:,:)
   real(r8), optional, intent(inout) :: v3(:,:,:,:)
   real(r8), optional, intent(inout) :: t3(:,:,:,:)
   real(r8), optional, intent(inout) :: ps(:,:,:)

!---------------------------Local workspace-----------------------------
   integer                        :: ioptop
   character(len=*), parameter    :: sub = "iop_update_prognostics"
!-----------------------------------------------------------------------
   ! set prognostics from iop
   ! Find level where tobs is no longer zero
   ioptop = minloc(tobs(:), 1, BACK=.true.)+1
   if (present(ps)) ps(1,1,timelevel)           = psobs
   if (present(t3)) t3(1,ioptop:,1,timelevel)   = tobs(ioptop:)
   if (present(q3)) q3(1,ioptop:,1,1,timelevel) = qobs(ioptop:)
   if (present(u3)) u3(1,ioptop:,1,timelevel)   = uobs(ioptop:)
   if (present(v3)) v3(1,ioptop:,1,timelevel)   = vobs(ioptop:)

 end subroutine iop_update_prognostics

end module iop
