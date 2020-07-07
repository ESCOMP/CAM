module dimensions_mod
  use shr_kind_mod, only: r8=>shr_kind_r8

  implicit none
  private


  !These are convenience variables for local use only, and are set to values in Atm%
  integer, public :: npx, npy, ntiles

  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1

  !
  ! The variables below hold indices of water vapor and condensate loading tracers as well as
  ! associated heat capacities (initialized in dyn_init):
  !
  !   qsize_condensate_loading_idx     = FV3 index of water tracers included in condensate loading according to FV3 dynamics
  !   qsize_condensate_loading_idx_gll = CAM index of water tracers included in condensate loading terms given FV3 index
  !
  integer,            allocatable, public :: qsize_tracer_idx_cam2dyn(:)
  character(len=16),  allocatable, public :: cnst_name_ffsl(:)     ! constituent names for FV3 tracers
  character(len=128), allocatable, public :: cnst_longname_ffsl(:) ! long name of FV3 tracers
  !
  !moist cp in energy conversion term
  !
  ! .false.: force dycore to use cpd (cp dry) instead of moist cp
  ! .true. : use moist cp in dycore
  !
  logical           , public :: fv3_lcp_moist   = .false. 
  logical           , public :: fv3_lcv_moist   = .false. 
  logical           , public :: fv3_scale_ttend = .false. 

end module dimensions_mod

