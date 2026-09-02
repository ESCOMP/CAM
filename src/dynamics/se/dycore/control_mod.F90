! This module contains constants and namelist variables used through out the model
! to avoid circular dependancies please do not 'use' any further modules here.
!
module control_mod
  use shr_kind_mod,     only: r8=>shr_kind_r8

  integer, public, parameter :: MAX_STRING_LEN=240
  integer, public, parameter :: MAX_FILE_LEN=240
!  character(len=MAX_STRING_LEN)    , public :: integration    ! time integration (only one currently supported is "explicit")

!shallow water advection tests:
!kmass points to a level with density.  other levels contain test tracers

  integer, public  :: tstep_type= 0                           ! 0 = leapfrog
                                                              ! 1 = RK (foward-in-time)
  integer, public  :: ftype = 2                               ! Forcing Type
  integer, public  :: ftype_conserve = 1  !conserve momentum (dp*u)
  integer, public  :: dribble_in_rsplit_loop = 0
  integer, public  :: statediag_numtrac = 3

  integer, public :: qsplit = 1           ! ratio of dynamics tsteps to tracer tsteps
  integer, public :: rsplit =-1           ! for vertically lagrangian dynamics, apply remap
                                          ! every rsplit tracer timesteps

  logical, public :: refined_mesh

  integer, public :: cubed_sphere_map = -1  ! -1 = chosen at run time
                                            !  0 = equi-angle Gnomonic (default)
                                            !  1 = equi-spaced Gnomonic (not yet coded)
                                            !  2 = element-local projection  (for var-res)
                                            !  3 = parametric (not yet coded)

  integer              , public :: limiter_option = 0

  integer              , public :: partmethod     ! partition methods
  character(len=MAX_STRING_LEN)    , public :: topology       ! options: "cube" is supported
  integer              , public :: statefreq      ! output frequency of synopsis of system state (steps)
  integer              , public :: runtype
  integer              , public :: numnodes
  integer              , public :: multilevel

  character(len=MAX_STRING_LEN)    , public :: columnpackage

  real (kind=r8), public :: tol            ! solver tolerance (convergence criteria)

  integer              , public :: fine_ne = -1              ! set for refined exodus meshes (variable viscosity)
  real (kind=r8), public :: max_hypervis_courant = 1d99 ! upper bound for Courant number
                                                               ! (only used for variable viscosity, recommend 1.9 in namelist)
  real (kind=r8), public :: nu      = 7.0D5           ! viscosity (momentum equ)
  real (kind=r8), public :: nu_div  = -1              ! viscsoity (momentum equ, div component)
  real (kind=r8), public :: nu_t    = -1              ! default = nu   T equ. viscosity
  real (kind=r8), public :: nu_q    = -1              ! default = nu   tracer viscosity
  real (kind=r8), public :: nu_p    = 0.0D5           ! default = 0    ps equ. viscosity
  real (kind=r8), public :: nu_top  = 0.0D5           ! top-of-the-model viscosity

  !
  ! Del4 sponge layer diffusion
  !
  ! Divergence damping hyperviscosity coefficient nu_div [m^4/s] for u,v is increased to
  ! nu_div*sponge_del4_nu_div_fac following a hyperbolic tangent function
  ! centered around pressure at vertical index sponge_del4_lev
  !
  ! Similar for sponge_del4_nu_fac
  !
  real(r8), public :: sponge_del4_nu_fac
  real(r8), public :: sponge_del4_nu_div_fac
  integer , public :: sponge_del4_lev


  integer, public :: hypervis_subcycle=1    ! number of subcycles for hyper viscsosity timestep
  integer, public :: hypervis_subcycle_sponge=1    ! number of subcycles for hyper viscsosity timestep in sponge
  integer, public :: hypervis_subcycle_q=1  ! number of subcycles for hyper viscsosity timestep on TRACERS
  logical, public :: cslam_q_filter = .false. ! mass-conservative del4 damping of water vapor on the CSLAM grid
  real(r8), public :: cslam_q_filter_nu_fac = 0.5_r8 ! factor on nu_p for the CSLAM-grid del4
  logical, public :: del4_cslam_qgll = .true. ! apply del4 hyperviscosity to water vapor on GLL after cslam2gll
  logical, public :: gll_advect_q = .false.   ! GLL double advection of thermodynamic tracers under CSLAM
                                              ! (cslam2gll overwrite not used)
  integer, public :: cslam_q_filter_nsub=1  ! number of subcycles for the CSLAM-grid del4 Q filter
                                            ! (apply_cslam_q_filter_del4); auto-set in print_cfl from the
                                            ! 2D del4 von Neumann stability bound
  integer, public :: hypervis_subcycle_cslam_q=2 ! number of subcycles for the GLL-side del4 on water vapor
                                            ! after cslam2gll (hypervis_Qdp); auto-set in print_cfl when
                                            ! the CSLAM Q filter is active, else the operational value 2

  real (kind=r8), public :: hypervis_power=0     ! if not 0, use variable hyperviscosity based on element area
  real (kind=r8), public :: hypervis_scaling=0      ! use tensor hyperviscosity

!
!three types of hyper viscosity are supported right now:
! (1) const hv:    nu * del^2 del^2
! (2) scalar hv:   nu(lat,lon) * del^2 del^2
! (3) tensor hv,   nu * ( \div * tensor * \grad ) * del^2
!
! (1) default:  hypervis_power=0, hypervis_scaling=0
! (2) Original version for var-res grids. (M. Levy)
!            scalar coefficient within each element
!            hypervisc_scaling=0
!            set hypervis_power>0 and set fine_ne, max_hypervis_courant
! (3) tensor HV var-res grids
!            tensor within each element:
!            set hypervis_scaling > 0 (typical values would be 3.2 or 4.0)
!            hypervis_power=0
!            (\div * tensor * \grad) operator uses cartesian laplace
!

  real (kind=r8), public :: initial_global_ave_dry_ps = 0._r8 ! scale dry surface pressure to initial_global_ave_dry_ps

  integer, public, parameter :: west  = 1
  integer, public, parameter :: east  = 2
  integer, public, parameter :: south = 3
  integer, public, parameter :: north = 4

  integer, public, parameter :: swest = 5
  integer, public, parameter :: seast = 6
  integer, public, parameter :: nwest = 7
  integer, public, parameter :: neast = 8

  !
  ! molecular diffusion
  !
  real(r8), public :: molecular_diff = -1.0_r8

  integer, public  :: vert_remap_uvTq_alg, vert_remap_tracer_alg

  integer, public :: pgf_formulation = -1 !PGF formulation - see prim_advance_mod.F90

  real(r8), public :: min_temperature = 0._r8

end module control_mod
