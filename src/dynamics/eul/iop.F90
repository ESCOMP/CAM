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

  integer :: closelatidx,closelonidx,latid,lonid,levid,timeid

  real(r8):: closelat,closelon

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_iop_fields
!  public :: scam_use_iop_srf
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


end module iop

