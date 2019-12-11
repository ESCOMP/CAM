module dimensions_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use constituents, only: qsize_d=>pcnst ! _EXTERNAL

  implicit none
  private


  !These are convenience variables for local use only, and are set to values in Atm%
  integer, public :: npx, npy, npz, ncnst, pnats, dnats
  integer, public :: nq                       ! transported tracers

  integer,               public :: qsize_condensate_loading = 1 !how many water variables to include in full density
  !
  ! The variables below hold indices of water vapor and condensate loading tracers as well as
  ! associated heat capacities (initialized in dyn_init):
  !
  !   qsize_condensate_loading_idx     = FV3 index of water tracers included in condensate loading according to FV3 dynamics
  !   qsize_condensate_loading_idx_gll = CAM index of water tracers included in condensate loading terms given FV3 index
  !
  integer,            allocatable, public :: qsize_tracer_idx_cam2dyn(:)
  integer,            allocatable, public :: qsize_condensate_loading_idx(:)    
  integer,            allocatable, public :: qsize_condensate_loading_idx_gll(:)
  real(r8),           allocatable, public :: qsize_condensate_loading_cp(:)
  real(r8),           allocatable, public :: qsize_condensate_loading_cv(:)
  character(len=16),  allocatable, public :: cnst_name_ffsl(:)     ! constituent names for FV3 tracers
  character(len=128), allocatable, public :: cnst_longname_ffsl(:) ! long name of FV3 tracers
  !
  !moist cp in energy conversion term
  !
  ! .false.: force dycore to use cpd (cp dry) instead of moist cp
  ! .true. : use moist cp in dycore
  !
  logical           , public :: fv3_lcp_moist = .false. 
  logical           , public :: fv3_lcv_moist = .false. 

  logical           , public :: fv3_scale_ttend = .false. 
  
   integer         :: qsize = 0 !qsize is set in dyn_comp

  ! fvm dimensions:
  logical, public :: lprint!for debugging


  integer, parameter :: plon  = 1      ! number of longitudes
  integer, parameter :: plev  = PLEV   ! number of vertical levels
  integer, parameter :: plat  = 1      ! number of latitudes
  
  integer, parameter :: plevp = plev + 1  ! plev + 1
  
  logical :: dyndecomp_set = .false.      ! flag indicates dynamics grid has been set
  
  integer :: beglat     ! beg. index for latitudes owned by a given proc
  integer :: endlat     ! end. index for latitudes owned by a given proc
  integer :: beglon     ! beg. index for longitudes owned by a given proc
  integer :: endlon     ! end. index for longitudes owned by a given proc
  integer :: numlats    ! number of latitudes owned by a given proc
  integer :: numlons    ! number of longitudes owned by a given proc
  
#if ( ! defined SPMD )
  parameter (beglat   = 1)
  parameter (endlat   = plat)
  parameter (beglon   = 1)
  parameter (endlon   = plon)
  parameter (numlats  = plat)
  parameter (numlons  = plon)
#endif
  
  public :: qsize

end module dimensions_mod

