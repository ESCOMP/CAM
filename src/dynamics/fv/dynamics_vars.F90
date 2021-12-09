module dynamics_vars

!-----------------------------------------------------------------------
! CAM fvcore internal variables
!
! !REVISION HISTORY:
!   01.06.06   Sawyer     Consolidated from various code snippets
!   03.06.25   Sawyer     Cleaned up, used ParPatternCopy (Create)
!   03.08.05   Sawyer     Removed rayf_init and hswf_init, related vars
!   03.10.22   Sawyer     pmgrid removed (now spmd_dyn)
!   03.11.18   Sawyer     Removed set_eta (ak, bk, now read from restart)
!   03.12.04   Sawyer     Moved T_FVDYCORE_GRID here (removed some vars)
!   04.08.25   Sawyer     Removed all module data members, now GRID only
!   04.10.06   Sawyer     Added spmd_dyn vars here; ESMF transpose vars
!   05.04.12   Sawyer     Added support for r4/r8 tracers
!   05.05.24   Sawyer     CAM/GEOS5 merge (removed GEOS_mod dependencies)
!   05.06.10   Sawyer     Scaled down version for CAM (no ESMF)
!   05.11.10   Sawyer     Removed dyn_interface (now in dyn_comp)
!   06.03.01   Sawyer     Removed m_ttrans, q_to_qxy, qxy_to_q, etc.
!   06.05.09   Sawyer     Added CONSV to dyn_state (conserve energy)
!   06.08.27   Sawyer     Removed unused ESMF code for RouteHandle
!-----------------------------------------------------------------------

use shr_kind_mod,       only: r8=>shr_kind_r8
use pmgrid,             only: plon, plat, plev

use decompmodule,       only: decomptype
use ghostmodule,        only: ghosttype

use cam_logfile,        only: iulog
use cam_abortutils,     only: endrun

#if defined(SPMD)
use parutilitiesmodule, only: parpatterntype, REAL4, INT4
#endif

implicit none
private
save

public :: &
   t_fvdycore_vars,      &
   t_fvdycore_grid,      &
   t_fvdycore_constants, &
   t_fvdycore_state,     &
   grid_vars_init,       &
   dynamics_clean

#ifdef SPMD
public :: spmd_vars_init
#endif

! T_FVDYCORE_VARS contains the prognostic variables for FVdycore

type T_FVDYCORE_VARS
   real(r8), dimension(:,:,:  ), pointer     :: U      ! U winds (D-grid)
   real(r8), dimension(:,:,:  ), pointer     :: V      ! V winds (D-grid)
   real(r8), dimension(:,:,:  ), pointer     :: PT     ! scaled virtual pot. temp.
   real(r8), dimension(:,:,:  ), pointer     :: PE     ! Pressure at layer edges
   real(r8), dimension(:,:,:  ), pointer     :: PKZ    ! P^kappa mean
   real(r8), dimension(:,:,:,:), pointer     :: tracer ! Tracers
end type T_FVDYCORE_VARS

! T_FVDYCORE_GRID contains information about the horizontal and vertical
! discretization and decompositions.
 
type T_FVDYCORE_GRID

   ! PILGRIM communication information

   integer :: twod_decomp = 0  ! 1 for multi-2D decompositions, 0 otherwise
              ! To assure that the latitudinal decomposition operates
              ! as efficiently as before, a separate parameter "twod_decomp" has
              ! been defined; a value of 1 refers to the multi-2D decomposition with
              ! transposes; a value of 0 means that the decomposition is effectively
              ! one-dimensional, thereby enabling the transpose logic to be skipped;
              ! there is an option to force computation of transposes even for case
              ! where decomposition is effectively 1-D.

   integer :: npes_xy= 1    ! number of PEs for XY decomposition
   integer :: npes_yz= 1    ! number of PEs for YZ decomposition
   integer :: myid_y = 0    ! subdomain index (0-based) in latitude (y)
   integer :: myid_z = 0    ! subdomain index (0 based) in level (z)
   integer :: npr_y  = 1    ! number of subdomains in y
   integer :: npr_z  = 1    ! number of subdomains in z

   integer :: myidxy_x = 0  ! subdomain index (0-based) in longitude (x) (second. decomp.)
   integer :: myidxy_y = 0  ! subdomain index (0 based) in latitude (y) (second. decomp.)
   integer :: nprxy_x = 1   ! number of subdomains in x (second. decomp.)
   integer :: nprxy_y = 1   ! number of subdomains in y (second. decomp.)
   integer :: iam = 0       ! 

   integer :: mod_method = 0  ! 1 for mpi derived types with transposes, 0 for contiguous buffers
   integer :: mod_geopk = 0   ! 1 for mpi derived types with transposes, 0 for contiguous buffers
   integer :: mod_gatscat = 0 ! 1 for mpi derived types with transposes, 0 for contiguous buffers

   type(decomptype) :: strip2d, strip2dx, strip3dxyz, strip3dxzy,          &
                       strip3dxyzp, strip3zaty, strip3dxzyp,               &
                       strip3yatz, strip3yatzp, strip3zatypt,              &
                       strip3kxyz, strip3kxzy, strip3kxyzp, strip3kxzyp,   &
                       strip3dyz, checker3kxy

   integer :: commdyn           ! communicator for all dynamics
   integer :: commxy            ! communicator for XY decomposition
   integer :: commyz            ! communicator for YZ decomposition
   integer :: commnyz           ! communicator for multiple YZ decomposition

   integer :: comm_y            ! communicator in latitude
   integer :: comm_z            ! communicator in vertical
   integer :: commxy_x          ! communicator in longitude (xy second. decomp.)
   integer :: commxy_y          ! communicator in latitude (xy second. decomp.)
   logical :: geopkdist         ! use distributed method for geopotential calculation 
                                !  with 2D decomp.
   logical :: geopk16byte       ! use Z-parallel distributed method for geopotential 
                                !  calculation with 2D decomp.; otherwise use Z-serial
                                !  pipeline algorithm when using distributed algoritm
   integer :: geopkblocks       ! number of stages to use in Z-serial pipeline 
                                !  (non-transpose) geopotential algorithm
   integer :: modc_dynrun(4)    ! 1: mod_comm irregular underlying communication method for dyn_run/misc
                                ! 2: mod_comm irregular communication handshaking for dyn_run/misc
                                ! 3: mod_comm irregular communication send protocol for dyn_run/misc
                                ! 4: mod_comm irregular communication nonblocking request throttle for dyn_run/misc
   integer :: modc_cdcore(4)    ! 1: mod_comm irregular underlying communication method for cd_core/geopk
                                ! 2: mod_comm irregular communication handshaking for cd_core/geopk
                                ! 3: geopk_d and mod_comm irregular communication send protocol for cd_core/geopk
                                ! 4: mod_comm irregular communication nonblocking request throttle for cd_core/geopk
   integer :: modc_gather(4)    ! 1: mod_comm irregular underlying communication method for gather
                                ! 2: mod_comm irregular communication handshaking for gather
                                ! 3: mod_comm irregular communication send protocol for gather
                                ! 4: mod_comm irregular communication nonblocking request throttle for gather
   integer :: modc_scatter(4)   ! 1: mod_comm irregular underlying communication method for scatter
                                ! 2: mod_comm irregular communication handshaking for scatter
                                ! 3: mod_comm irregular communication send protocol for scatter
                                ! 4: mod_comm irregular communication nonblocking request throttle for scatter 
   integer :: modc_tracer(4)    ! 1: mod_comm irregular underlying communication method for multiple tracers
                                ! 2: mod_comm irregular communication handshaking for multiple tracers
                                ! 3: mod_comm irregular communication send protocol for multiple tracers
                                ! 4: mod_comm irregular communication nonblocking request throttle for multiple tracers 
   integer :: modc_onetwo       ! one or two simultaneous mod_comm irregular communications (excl. tracers)
   integer :: modc_tracers      ! max number of tracers for simultaneous mod_comm irregular communications

#if defined(SPMD)
   type (ghosttype)        :: ghostu_yz, ghostv_yz, ghostpt_yz,                   &
                              ghostpe_yz, ghostpkc_yz
   type (parpatterntype)   :: u_to_uxy, uxy_to_u, v_to_vxy, vxy_to_v,             &
                              ikj_yz_to_xy, ikj_xy_to_yz,                         &
                              ijk_yz_to_xy, ijk_xy_to_yz,                         &
                              pe_to_pexy, pexy_to_pe,                             &
                              pt_to_ptxy, ptxy_to_pt, pkxy_to_pkc,                &
                              r4_xy_to_yz, r4_yz_to_xy, q3_to_qxy3, qxy3_to_q3,   &
                              xy2d_to_yz2d, yz2d_to_xy2d, scatter_3d, gather_3d,  &
                              g_2dxy_r8, g_2dxy_r4, g_2dxy_i4,                    &
                              s_2dxy_r8, s_2dxy_r4, s_2dxy_i4,                    &
                              g_3dxyz_r8, g_3dxyz_r4, g_3dxyzp_r8, g_3dxyzp_r4,   &
                              s_3dxyz_r8, s_3dxyz_r4, s_3dxyzp_r8, s_3dxyzp_r4
#endif

   ! END PILGRIM communication information


   integer  :: JFIRST = 1       ! Start latitude (exclusive)
   integer  :: JLAST  = plat    ! End latitude (exclusive)
 
   integer  :: ng_c = 0         ! Ccore ghosting
   integer  :: ng_d = 0         ! Dcore ghosting
   integer  :: ng_s = 0         ! Staggered grid ghosting for
                                ! certain arrays, max(ng_c+1,ng_d)
   ! For 2D decomposition

   integer  :: IFIRSTXY = 1     ! Start longitude (exclusive)
   integer  :: ILASTXY  = plon  ! End longitude (exclusive)
   integer  :: JFIRSTXY = 1     ! Start latitude (exclusive)
   integer  :: JLASTXY  = plat  ! End latitude (exclusive)

   integer  :: IM               ! Full longitude dim
   integer  :: JM               ! Full latitude dim (including poles)

   real(r8) :: DL
   real(r8) :: DP
   real(r8) :: ACAP
   real(r8) :: RCAP

   real(r8), dimension(:), pointer :: COSP             ! Cosine of lat angle -- volume mean
   real(r8), dimension(:), pointer :: SINP             ! Sine of lat angle -- volume mean
   real(r8), dimension(:), pointer :: COSE             ! Cosine at finite volume edge
   real(r8), dimension(:), pointer :: SINE             ! Sine at finite volume edge
   real(r8), dimension(:), pointer :: ACOSP            ! Reciprocal of cosine of lat angle

   real(r8), dimension(:), pointer :: ACOSU            ! Reciprocal of cosine of lat angle (staggered)

   real(r8), dimension(:), pointer :: COSLON           ! Cosine of longitudes - volume center
   real(r8), dimension(:), pointer :: SINLON           ! Sine of longitudes - volume center
   real(r8), dimension(:), pointer :: COSL5            ! Cosine of longitudes - volume center
   real(r8), dimension(:), pointer :: SINL5            ! Sine of longitudes - volume center
 
   ! Variables which are used repeatedly in CD_CORE

   integer ::       js2g0
   integer ::       jn2g0
   integer ::       jn1g1

   real(r8), pointer :: trigs(:)
   real(r8), pointer :: fc(:), f0(:)
   real(r8), pointer :: dc(:,:), de(:,:), sc(:), se(:)
   real(r8), pointer :: cdx(:,:), cdy(:,:)
   real(r8), pointer :: cdx4(:,:), cdy4(:,:)                           ! for div4 damping
   real(r8), pointer :: cdxde(:,:), cdxdp(:,:),cdyde(:,:),cdydp(:,:)   ! for del2 damping
   real(r8), pointer :: cdxdiv(:,:), cdydiv(:,:), cdtau4(:,:)          ! for del2 damping

   real(r8), pointer :: dcdiv4(:,:), dediv4(:,:), scdiv4(:), sediv4(:) ! for div4 damping

   real(r8), pointer :: dtdx(:), dtdxe(:), txe5(:), dtxe5(:)
   real(r8), pointer :: dyce(:),   dx(:) ,  rdx(:),    cy(:)
   real(r8), pointer :: dtdx2(:), dtdx4(:),  dxdt(:), dxe(:)
   real(r8), pointer :: cye(:),    dycp(:),  rdxe(:)

   real(r8) :: rdy, dtdy, dydt, dtdy5, tdy5
   real(r8) :: dt0 = 0

   integer  :: ifax(13)

   real(r8) ::  zt_c
   real(r8) ::  zt_d

   ! This part refers to the vertical grid

   integer                         :: KM              ! Numer of levels
   integer                         :: KMAX            ! KM+1 (?)

   ! For 2D decomposition

   integer                         :: KFIRST = 1      ! Start level (exclusive)
   integer                         :: KLAST  = plev   ! End level (exclusive)
   integer                         :: KLASTP = plev+1 ! klast+1, except km+1 when klastp=km+1

   integer                         :: KORD            ! monotonicity order for mapping (te_map)
   integer                         :: KS              ! Number of true pressure levels (out of KM+1)
   real(r8)                        :: PTOP            ! pressure at top (ak(1))
   real(r8)                        :: PINT            ! initial pressure (ak(km+1))
   real(r8), dimension(:), pointer :: AK              ! Sigma mapping
   real(r8), dimension(:), pointer :: BK              ! Sigma mapping

   ! Tracers

   integer                         :: NQ              ! Number of advected tracers
   integer                         :: NTOTQ           ! Total number of tracers (NQ <= NC)

   integer  :: ct_overlap    ! nonzero for overlap of cd_core and trac2d, 0 otherwise
   integer  :: trac_decomp   ! size of tracer domain decomposition for trac2d

   ! Extra subdomain bounds for cd_core/trac2d overlap and trac2d decomposition
   ! Relevant for secondary yz decomposition only; refers back to primary yz decomposition
   integer                         :: JFIRSTCT          ! jfirst
   integer                         :: JLASTCT           ! jlast
   integer                         :: KFIRSTCT          ! kfirst
   integer                         :: KLASTCT           ! klast

   ! Bounds for tracer decomposition
   integer, dimension(:), pointer  :: ktloa             ! lower tracer index (global map)
   integer, dimension(:), pointer  :: kthia             ! upper tracer index (global map)
   integer                         :: ktlo              ! lower tracer index (local)
   integer                         :: kthi              ! upper tracer index (local)

   logical :: high_alt  ! high-altitude physics parameters switch

end type T_FVDYCORE_GRID

! Constants used by fvcore
type T_FVDYCORE_CONSTANTS
   real(r8)                             :: pi
   real(r8)                             :: omega    ! angular velocity of earth's rotation  
   real(r8)                             :: cp       ! heat capacity of air at constant pressure
   real(r8)                             :: ae       ! radius of the earth (m)
   real(r8)                             :: rair     ! Gas constant of the air
   real(r8)                             :: cappa    ! Cappa?
   real(r8)                             :: zvir     ! RWV/RAIR-1
end type T_FVDYCORE_CONSTANTS

type t_fvdycore_state
   type (t_fvdycore_vars)      :: vars
   type (t_fvdycore_grid )     :: grid
   type (t_fvdycore_constants) :: constants
   real(r8) :: DT            ! Large time step
   real(r8) :: CHECK_DT      ! Time step to check maxmin
   integer  :: ICD, JCD      ! Algorithm orders (C Grid)
   integer  :: IORD, JORD    ! Algorithm orders (D Grid)
   integer  :: KORD          ! Vertical order
   integer  :: TE_METHOD     ! method for total energy mapping (te_map)
   logical  :: CONSV         ! dycore conserves tot. en.
   integer  :: NSPLIT
   integer  :: NSPLTRAC
   integer  :: NSPLTVRM
   integer  :: FILTCW        ! filter c-grid winds if positive
   integer  :: fft_flt       ! 0 => FFT/algebraic filter; 1 => FFT filter
   integer  :: div24del2flag ! 2 for 2nd order div damping, 4 for 4th order div damping,
                             ! 42 for 4th order div damping plus 2nd order velocity damping
   real(r8) :: del2coef      ! strength of 2nd order velocity damping
   logical  :: high_order_top! use normal 4-order PPM calculation near the model top
   logical  :: am_geom_crrct ! apply correction for angular momentum (AM) conservation in geometry
   logical  :: am_correction ! apply correction for angular momentum (AM) conservation in SW eqns
   logical  :: am_fixer      ! apply global fixer to conserve AM
   logical  :: am_fix_lbl    ! apply global AM fixer level by level
   logical  :: am_diag       ! turns on an AM diagnostic calculations
end type t_fvdycore_state

!========================================================================================
contains
!========================================================================================

#if defined(SPMD)

subroutine spmd_vars_init(imxy, jmxy, jmyz, kmyz, grid)

   ! Initialize SPMD related variables.
   ! !REVISION HISTORY: 
   !   02.11.08    Sawyer       Creation
   !   03.05.07    Sawyer       Use ParPatternCopy for q_to_qxy, etc.
   !   03.07.23    Sawyer       Removed dependency on constituents module
   !   03.09.10    Sawyer       Reactivated u_to_uxy, etc, redefined pe2pexy
   !   03.11.19    Sawyer       Merged in CAM code with mod_method
   !   04.08.25    Sawyer       Added GRID as argument
   
   use decompmodule, only: decompcreate, decompfree
   use ghostmodule, only : ghostcreate, ghostfree
   use parutilitiesmodule, only : gid, parpatterncreate, parsplit
   use mpishorthand, only: mpiint

   ! Arguments

   integer, dimension(:), intent(in) :: imxy
   integer, dimension(:), intent(in) :: jmxy
   integer, dimension(:), intent(in) :: jmyz
   integer, dimension(:), intent(in) :: kmyz

   type( t_fvdycore_grid ), intent(inout) :: grid

   ! local variables:

   type(decomptype) :: global2d, local2d

   integer :: im, jm, km         !  Global dims
   integer :: nq

   integer :: nprxy_x    ! XY decomp - Nr in X
   integer :: nprxy_y    ! XY decomp - Nr in Y
   integer :: npryz_y    ! YZ decomp - Nr in Y
   integer :: npryz_z    ! YZ decomp - Nr in Z
   integer :: npes_xy    ! XY decomp - Total nr
   integer :: npes_yz    ! YZ decomp - Total nr

   integer :: commxy     ! Communicator for XY decomp
   integer :: commyz     ! Communicator for YZ decomp
   integer :: commnyz    ! Communicator for multiple YZ decomp

   integer :: jfirstxy, jlastxy
   integer :: jfirst, jlast
   integer :: kfirst, klast

   integer :: ng_s, ng_d   !  Ghost widths

   integer :: rank_y, rank_z, rankxy_x, rankxy_y  ! Currently not used
   integer :: size_y, size_z, sizexy_x, sizexy_y  ! Currently not used

   integer :: xdist(1), ydistk(1), zdist1(1), zdistxy(1) ! non-distributed dims
   integer, allocatable :: xdist_global(:), ydist_global(:) 
   integer, allocatable :: zdist(:) ! number of levels per subdomain

   integer :: ig1, ig2, jg1, jg2
   integer :: jg1d, jg2d, jg1s, jg2s
   integer :: kg1, kg2, kg2p

   integer :: ct_overlap
   integer :: trac_decomp
   integer :: ktmod, ml
   integer :: myidmod
   integer :: ictstuff(4)
   integer :: kquot, krem, krun, kt, mlt
   !---------------------------------------------------------------------------

   im = grid%im
   jm = grid%jm
   km = grid%km
   nq = grid%nq

   nprxy_x = grid%nprxy_x
   nprxy_y = grid%nprxy_y
   npryz_y = grid%npr_y
   npryz_z = grid%npr_z
   npes_xy = grid%npes_xy
   npes_yz = grid%npes_yz

   commxy  = grid%commxy
   commyz  = grid%commyz
   commnyz = grid%commnyz

   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy

   jfirst = grid%jfirst
   jlast  = grid%jlast
   kfirst = grid%kfirst
   klast  = grid%klast

   ng_s   = grid%ng_s
   ng_d   = grid%ng_d

   ! Split communicators
   call parsplit(commyz, grid%myid_z, gid, grid%comm_y, rank_y, size_y)
   call parsplit(commyz, grid%myid_y, gid, grid%comm_z, rank_z, size_z)
   call parsplit(commxy, grid%myidxy_y, gid, grid%commxy_x, rankxy_x, sizexy_x)
   call parsplit(commxy, grid%myidxy_x, gid, grid%commxy_y, rankxy_y, sizexy_y)


   ! create decompositions for CAM data structures

   allocate(xdist_global(nprxy_x))
   allocate(ydist_global(nprxy_y))
   allocate(zdist (npryz_z))
   xdist(1) = im

   ! Create PILGRIM decompositions (see decompmodule)

   if (gid < npes_xy) then
      xdist_global = 0
      ydist_global = 0
      xdist_global(1) = im
      ydist_global(1) = jm
      call decompcreate(nprxy_x, nprxy_y, xdist_global,            &
                        ydist_global, global2d)
      call decompcreate(nprxy_x, nprxy_y, imxy, jmxy, local2d )

      ! Decompositions needed on xy decomposition for parpatterncreate

      call decompcreate( 1, npryz_y, xdist, jmyz, grid%strip2d )
      call decompcreate( 1, npryz_y, npryz_z, xdist,                &
                         jmyz, kmyz, grid%strip3dxyz )
      call decompcreate( "xzy", 1, npryz_z, grid%npr_y, xdist,      &
                           kmyz, jmyz, grid%strip3dxzy )

      ! For y communication within z subdomain (klast version)
      ! Use myidmod to have valid index for inactive processes
      !  for smaller yz decomposition
      myidmod = mod(grid%myid_z, grid%npr_z)  ! = myid_z for active yz process
      zdist1(1) = kmyz(myidmod+1)
      call decompcreate( 1, npryz_y, 1, xdist, jmyz, zdist1,        &
                         grid%strip3yatz )

      ! For z communication within y subdomain

      ydistk(1) = jmyz(grid%myid_y+1)
      call decompcreate( 1, 1, npryz_z, xdist, ydistk, kmyz,        &
                         grid%strip3zaty )

      ! Arrays dimensioned plev+1

      zdist(:) = kmyz(:)
      zdist(npryz_z) = kmyz(npryz_z) + 1
      call decompcreate( 1, npryz_y, npryz_z, xdist, jmyz, zdist,&
                         grid%strip3dxyzp )
      call decompcreate( "xzy", 1, npryz_z, npryz_y,                &
                         xdist, zdist, jmyz, grid%strip3dxzyp )

      ! Arrays dimensioned plev+1, within y subdomain

      ydistk(1) = jmyz(grid%myid_y+1)
      call decompcreate( "xzy", 1, npryz_z, 1, xdist, zdist, ydistk,   &
                         grid%strip3zatypt )

      ! For y communication within z subdomain (klast+1 version)
      ! Use myidmod to have valid index for inactive processes
      !  for smaller yz decomposition
      myidmod = mod(grid%myid_z, grid%npr_z)  ! = myid_z for active yz process
      zdist1(1) = kmyz(myidmod+1)
      call decompcreate( 1, npryz_y, 1, xdist, jmyz, zdist1,        &
                         grid%strip3yatzp )

      ! For the 2D XY-YZ data transfer, we need a short 3D array
      zdist(:) = 1    ! One copy on each z PE set
      call decompcreate( 1, npryz_y, npryz_z,                       &
                         xdist, jmyz, zdist, grid%strip3dyz )
   end if

   ! Secondary xy decomposition

   if (grid%twod_decomp == 1) then
      if (gid < npes_xy) then
         zdistxy(1) = npryz_z     ! All npr_z copies on 1 PE
         call decompcreate( nprxy_x, nprxy_y, 1,                     &
                            imxy, jmxy, zdistxy, grid%checker3kxy )
         zdistxy(1) = km
         call decompcreate( nprxy_x, nprxy_y, 1,                     &
                            imxy, jmxy, zdistxy, grid%strip3kxyz )
         call decompcreate( "xzy", nprxy_x, 1, nprxy_y,              &
                            imxy, zdistxy, jmxy, grid%strip3kxzy )

         zdistxy(1) = zdistxy(1) + 1
         call decompcreate( nprxy_x, nprxy_y, 1,                     &
                            imxy, jmxy, zdistxy, grid%strip3kxyzp )
         call decompcreate( "xzy", nprxy_x, 1, nprxy_y,              &
                            imxy, zdistxy, jmxy, grid%strip3kxzyp )
         zdistxy(1) = jlastxy - jfirstxy + 1
         call decompcreate( nprxy_x, 1, imxy, zdistxy, grid%strip2dx )
      end if
   end if

   deallocate(zdist)
   deallocate(ydist_global)
   deallocate(xdist_global)

   if ( grid%twod_decomp == 1 ) then

      ! Initialize ghost regions

      ! Set limits for ghostcreate
      ig1 = 1
      ig2 = im
      jg1 = jfirst
      jg2 = jlast
      jg1d = jfirst-ng_d
      jg1s = jfirst-ng_s
      jg2d = jlast+ng_d
      jg2s = jlast+ng_s
      kg1 = kfirst
      kg2 = klast
      kg2p = klast+1

      ! Call ghostcreate with null ranges for non-yz processes
      if (gid >= npes_yz) then
         ig1 = im/2
         ig2 = ig1 - 1
         jg1 = (jfirst+jlast)/2
         jg2 = jg1 - 1
         jg1d = jg1
         jg1s = jg1
         jg2d = jg2
         jg2s = jg2
         kg1 = (kfirst+klast)/2
         kg2 = kg1 - 1
         kg2p = kg2
      end if

      ! Ghosted decompositions needed on xy decomposition for parpatterncreate
      if (gid < npes_xy) then
         call ghostcreate( grid%strip3dxyz, gid, im, ig1, ig2, .true., &
                           jm, jg1d, jg2s, .false., &
                           km, kg1, kg2, .false., grid%ghostu_yz )
         call ghostcreate( grid%strip3dxyz, gid, im, ig1, ig2, .true., &
                           jm, jg1s, jg2d, .false., &
                           km, kg1, kg2, .false., grid%ghostv_yz )
         call ghostcreate( grid%strip3dxyz, gid, im, ig1, ig2, .true., &
                           jm, jg1d, jg2d, .false., &
                           km, kg1, kg2, .false., grid%ghostpt_yz )
         call ghostcreate( grid%strip3dxzyp, gid, im, ig1, ig2, .true., &
                           km+1, kg1, kg2p, .false., &
                           jm, jg1, jg2, .false., grid%ghostpe_yz)
         call ghostcreate( grid%strip3dxyzp, gid, im, ig1, ig2, .true., &
                           jm, jg1, jg2, .false.,       &
                           km+1, kg1, kg2p, .false., grid%ghostpkc_yz)
      end if

      ! Initialize transposes

      if (gid < npes_xy) then
         call parpatterncreate(commxy, grid%ghostu_yz, grid%strip3kxyz, &
                               grid%u_to_uxy, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3kxyz,grid%ghostu_yz, &
                               grid%uxy_to_u, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%ghostv_yz, grid%strip3kxyz, &
                               grid%v_to_vxy, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3kxyz, grid%ghostv_yz, &
                               grid%vxy_to_v, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3dxyz, grid%strip3kxyz,&
                               grid%ijk_yz_to_xy, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3kxyz, grid%strip3dxyz,&
                               grid%ijk_xy_to_yz, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3dxzy, grid%strip3kxzy,&
                               grid%ikj_yz_to_xy, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3kxzy, grid%strip3dxzy,&
                               grid%ikj_xy_to_yz, mod_method=grid%mod_method)

         ! Note PE <-> PEXY has been redefined for PEXY ijk, but PE ikj

         call parpatterncreate(commxy, grid%ghostpe_yz, grid%strip3kxzyp, &
                               grid%pe_to_pexy, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3kxzyp, grid%ghostpe_yz, &
                               grid%pexy_to_pe, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%ghostpt_yz, grid%strip3kxyz,  &
                               grid%pt_to_ptxy, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3kxyz, grid%ghostpt_yz,  &
                               grid%ptxy_to_pt, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3dxyz, grid%strip3kxyz,  &
                               grid%r4_yz_to_xy, mod_method=grid%mod_method,  &
                               T = REAL4 )
         call parpatterncreate(commxy, grid%strip3kxyz, grid%strip3dxyz,  &
                               grid%r4_xy_to_yz, mod_method=grid%mod_method,  &
                               T = REAL4 )
         call parpatterncreate(commxy, grid%strip3kxyzp, grid%ghostpkc_yz, &
                               grid%pkxy_to_pkc, mod_method=grid%mod_method)

         ! These are for 'transposing' 2D arrays from XY YZ
         call parpatterncreate(commxy, grid%checker3kxy, grid%strip3dyz, &
                               grid%xy2d_to_yz2d, mod_method=grid%mod_method)
         call parpatterncreate(commxy, grid%strip3dyz, grid%checker3kxy, &
                               grid%yz2d_to_xy2d, mod_method=grid%mod_method)
      end if

      ! Free unneeded decompositions

      call decompfree(grid%strip3dxzyp)
      call decompfree(grid%strip3dyz)
      call decompfree(grid%strip3yatz)
      call decompfree(grid%strip3yatzp)
      call decompfree(grid%strip3zaty)
      call decompfree(grid%strip3zatypt)
      call decompfree(grid%strip3kxyz)
      call decompfree(grid%strip3kxzy)
      call decompfree(grid%strip3kxyzp)
      call decompfree(grid%strip3kxzyp)
      call decompfree(grid%checker3kxy)

      call ghostfree(grid%ghostu_yz)
      call ghostfree(grid%ghostv_yz)
      call ghostfree(grid%ghostpt_yz)
      call ghostfree(grid%ghostpe_yz)
      call ghostfree(grid%ghostpkc_yz)

   end if

   ! Define scatter and gather patterns for 2D and 3D unghosted arrays

   if (gid < npes_xy) then
      call parpatterncreate( commxy, global2d, local2d, grid%s_2dxy_r8, &
                             mod_method=grid%mod_gatscat )
      call parpatterncreate( commxy, local2d, global2d, grid%g_2dxy_r8,  &
                             mod_method=grid%mod_gatscat )

      call parpatterncreate( commxy, global2d, local2d, grid%s_2dxy_r4, &
                             mod_method=grid%mod_gatscat, T = REAL4 )
      call parpatterncreate( commxy, local2d, global2d, grid%g_2dxy_r4,  &
                             mod_method=grid%mod_gatscat, T = REAL4 )

      call parpatterncreate( commxy, global2d, local2d, grid%s_2dxy_i4, &
                             mod_method=grid%mod_gatscat, T = INT4 )
      call parpatterncreate( commxy, local2d, global2d, grid%g_2dxy_i4,  &
                             mod_method=grid%mod_gatscat, T = INT4 )

      ! 3D XYZ patterns, will replace XZY patterns eventually

      call parpatterncreate( commxy, grid%s_2dxy_r8, grid%s_3dxyz_r8, km )
      call parpatterncreate( commxy, grid%g_2dxy_r8, grid%g_3dxyz_r8, km )
      call parpatterncreate( commxy, grid%s_2dxy_r8, grid%s_3dxyzp_r8, km+1 )
      call parpatterncreate( commxy, grid%g_2dxy_r8, grid%g_3dxyzp_r8, km+1 )

      call parpatterncreate( commxy, grid%s_2dxy_r4, grid%s_3dxyz_r4, km )
      call parpatterncreate( commxy, grid%g_2dxy_r4, grid%g_3dxyz_r4, km )
      call parpatterncreate( commxy, grid%s_2dxy_r4, grid%s_3dxyzp_r4, km+1 )
      call parpatterncreate( commxy, grid%g_2dxy_r4, grid%g_3dxyzp_r4, km+1 )

      call decompfree( global2d )
      call decompfree( local2d )
   end if

   ! Secondary subdomain limits for cd_core/trac2d overlap and trac2d decomposition

   ct_overlap    = grid%ct_overlap
   trac_decomp   = grid%trac_decomp
   grid%jfirstct = grid%jfirst
   grid%jlastct  = grid%jlast
   grid%kfirstct = grid%kfirst
   grid%klastct  = grid%klast
   if (ct_overlap > 0) then
      mlt = 2
   elseif (trac_decomp .gt. 1) then
      mlt = trac_decomp
   else
      mlt = 1
   end if

   if (mlt > 1) then
      if (gid < npes_yz) then
         ictstuff(1) = grid%jfirstct
         ictstuff(2) = grid%jlastct
         ictstuff(3) = grid%kfirstct
         ictstuff(4) = grid%klastct
         do ml = 2, mlt
            call mpisend(ictstuff, 4, mpiint, gid+(ml-1)*npes_yz, gid+(ml-1)*npes_yz, commnyz)
         enddo
      elseif (gid < mlt*npes_yz) then
         ktmod = gid/npes_yz
         call mpirecv(ictstuff, 4, mpiint, gid-ktmod*npes_yz, gid, commnyz)
         grid%jfirstct = ictstuff(1)
         grid%jlastct  = ictstuff(2)
         grid%kfirstct = ictstuff(3)
         grid%klastct  = ictstuff(4)
      end if
   end if

   if (trac_decomp .gt. 1) then
      kquot = nq / trac_decomp
      krem = nq - kquot * trac_decomp
      krun = 0
      do kt = 1, krem
         grid%ktloa(kt) = krun + 1
         krun = krun + kquot + 1
         grid%kthia(kt) = krun
      enddo
      do kt = krem+1, trac_decomp
         grid%ktloa(kt) = krun + 1
         krun = krun + kquot
         grid%kthia(kt) = krun
      enddo
      ktmod = gid/npes_yz + 1
      ktmod = min(ktmod, trac_decomp)
      grid%ktlo = grid%ktloa(ktmod)
      grid%kthi = grid%kthia(ktmod)
   endif

end subroutine spmd_vars_init

#endif

!========================================================================================

subroutine grid_vars_init(pi, ae, om, dt, fft_flt, &
                          am_geom_crrct, grid)

   ! Initialize FV specific GRID vars
   ! 
   ! !REVISION HISTORY: 
   !   00.01.10    Grant        Creation using code from SJ Lin
   !   01.06.06    Sawyer       Modified for dynamics_vars
   !   04.08.25    Sawyer       Now updates GRID
   !   05.06.30    Sawyer       Added initializations from cd_core
   !   06.09.15    Sawyer       PI now passed as argument

   use pft_module, only: pftinit, pft_cf

   ! Arguments
   real(r8), intent(in) :: pi
   real(r8), intent(in) :: ae     ! radius of the earth (m)
   real(r8), intent(in) :: om     ! angular velocity of earth's rotation 
   real(r8), intent(in) :: dt
   integer,  intent(in) :: fft_flt
   logical,  intent(in) :: am_geom_crrct

   type( T_FVDYCORE_GRID ), intent(inout) :: grid

   ! local variables:
   integer :: im
   integer :: jm

   real(r8) :: ph5      ! This is to ensure 64-bit for any choice of r8

   integer  :: i, j, imh
   real(r8) :: zam5, zamda
   integer  :: js2g0, jn2g0, jn1g1, js2gc, jn1gc
   integer  :: js2gs, jn2gd, jn1gs

   real(r8), pointer :: cosp(:), sinp(:), cose(:), sine(:), acosp(:), acosu(:)
   real(r8), pointer :: coslon(:), sinlon(:), cosl5(:), sinl5(:)

   real(r8) :: rat, ycrit
   real(r8) :: dt5

   character(len=*), parameter :: sub='grid_vars_init'
   !---------------------------------------------------------------------------

   im = grid%im
   jm = grid%jm

   grid%dl = (pi+pi)/im
   grid%dp = pi/(jm-1)

   allocate(grid%cosp(jm))
   allocate(grid%sinp(jm))
   allocate(grid%cose(jm))
   allocate(grid%sine(jm))
   allocate(grid%acosp(jm))
   allocate(grid%acosu(jm))

   allocate(grid%coslon(im))
   allocate(grid%sinlon(im))
   allocate(grid%cosl5(im))
   allocate(grid%sinl5(im))

   cosp => grid%cosp
   sinp => grid%sinp
   cose => grid%cose
   sine => grid%sine
   acosp  => grid%acosp
   acosu  => grid%acosu

   coslon => grid%coslon
   sinlon => grid%sinlon
   cosl5  => grid%cosl5
   sinl5  => grid%sinl5

   ! philosophy below: edge values = local true values; centred values = area averages

   do j = 2, jm
      ph5  = -0.5_r8*pi + ((j-1)-0.5_r8)*(pi/(jm-1))
      sine(j) = sin(ph5)
   end do

   ! cos(theta) at cell center distretized as
   !
   ! cos(theta) = d(sin(theta))/d(theta)

   cosp( 1) =  0._r8
   cosp(jm) =  0._r8
   do j = 2, jm-1
      cosp(j) = (sine(j+1)-sine(j)) / grid%dp
   end do

   ! Define cosine at edges..
   
   if (am_geom_crrct) then 
      do j = 2, jm
         ph5     = -0.5_r8*pi + ((j-1)-0.5_r8)*(pi/(jm-1._r8))
         cose(j) = cos(ph5)
      end do
   else
      do j = 2, jm
         cose(j) = 0.5_r8 * (cosp(j-1) + cosp(j))  ! dsine/dpe between j+1 and j-1
      end do
   end if
   cose(1) = cose(2)

   do j = 2, jm-1
      acosu(j) = 2._r8 / (cose(j) + cose(j+1))
   end do

   sinp( 1) = -1._r8
   sinp(jm) =  1._r8
   if (am_geom_crrct) then 
      do j = 2, jm-1
         sinp(j) =        (cose(j) - cose(j+1))/grid%dp  ! sqrt(cosp^2+sinp^2)=1
      end do
   else
      do j = 2, jm-1
         sinp(j) = 0.5_r8 * (sine(j) + sine(j+1))        ! 2*sinp*cosp*dp=d(cose^2)
      end do
   end if

   ! Pole cap area and inverse
   grid%acap = im*(1._r8+sine(2)) / grid%dp
   grid%rcap = 1._r8 / grid%acap
 
   imh = im/2
   if (im /= 2*imh) then
      write(iulog,*) sub//': ERROR: im must be an even integer'
      call endrun(sub//': ERROR: im must be an even integer')
   end if
 
   ! Define logitude at the center of the volume
   ! i=1, Zamda = -pi
 
   do i = 1, imh
      zam5          = ((i-1)-0.5_r8) * grid%dl
      cosl5(i)      =  cos(zam5)
      cosl5(i+imh)  = -cosl5(i)
      sinl5(i)      =  sin(zam5)
      sinl5(i+imh)  = -sinl5(i)
      zamda         = (i-1)*grid%dl
      coslon(i)     =  cos(zamda)
      coslon(i+imh) = -coslon(i)
      sinlon(i)     =  sin(zamda)
      sinlon(i+imh) = -sinlon(i)
   end do

   acosp( 1) = grid%rcap * im
   acosp(jm) = grid%rcap * im
   do j = 2, jm-1
      acosp(j) = 1._r8 / cosp(j)
   enddo

   ! cd_core initializations

   allocate(grid%dtdx(jm))
   allocate(grid%dtdx2(jm))
   allocate(grid%dtdx4(jm))
   allocate(grid%dtdxe(jm))
   allocate(grid%dxdt(jm))
   allocate(grid%dxe(jm))
   allocate(grid%cye(jm))
   allocate(grid%dycp(jm))
   allocate(grid%rdxe(jm))
   allocate(grid%txe5(jm))
   allocate(grid%dtxe5(jm))
   allocate(grid%dyce(jm))
   allocate(grid%dx(jm))
   allocate(grid%rdx(jm))
   allocate(grid%cy(jm))

   js2g0  = max(2,grid%jfirst)
   jn2g0  = min(jm-1,grid%jlast)
   jn1g1  = min(jm,grid%jlast+1)
   js2gc  = max(2,grid%jfirst-grid%ng_c) ! NG lats on S (starting at 2)
   jn1gc  = min(jm,grid%jlast+grid%ng_c) ! ng_c lats on N (ending at jm)

   grid%js2g0  = js2g0
   grid%jn2g0  = jn2g0
   grid%jn1g1  = jn1g1

   js2gs = max(2,grid%jfirst-grid%ng_s)
   jn2gd = min(jm-1,grid%jlast+grid%ng_d)
   jn1gs = min(jm,grid%jlast+grid%ng_s)

   allocate(grid%sc(js2g0:jn2g0))
   allocate(grid%se(js2g0:jn1g1))
   allocate(grid%dc(im,js2g0:jn2g0))
   allocate(grid%de(im,js2g0:jn1g1))

   allocate(grid%scdiv4(js2gs:jn2gd))   !for filtering of u and v in div4 damping 
   allocate(grid%sediv4(js2gs:jn1gs))   !for filtering of u and v in div4 damping 
   allocate(grid%dcdiv4(im,js2gs:jn2gd))!for filtering of u and v in div4 damping 
   allocate(grid%dediv4(im,js2gs:jn1gs))!for filtering of u and v in div4 damping 

   call pftinit(im, fft_flt)

   ! Determine ycrit such that effective DX >= DY
   rat = real(im,r8)/real(2*(jm-1),r8)
   ycrit = acos( min(0.81_r8, rat) ) * (180._r8/pi)

   call pft_cf(im, jm, js2g0, jn2g0, jn1g1, &
               grid%sc, grid%se, grid%dc, grid%de,  &
               grid%cosp, grid%cose, ycrit)

   !for filtering of u and v in div4 damping 
   !(needs larger halo than cam3.5 code)
   call pft_cf(im, jm, js2gs, jn2gd, jn1gs,                             & 
               grid%scdiv4, grid%sediv4, grid%dcdiv4, grid%dediv4,      & 
               grid%cosp, grid%cose, ycrit)                               

   allocate( grid%cdx   (js2g0:jn1g1,grid%kfirst:grid%klast) )
   allocate( grid%cdy   (js2g0:jn1g1,grid%kfirst:grid%klast) )

   allocate( grid%cdx4  (js2g0:jn1g1,grid%kfirst:grid%klast) )!for div4 damping
   allocate( grid%cdy4  (js2g0:jn1g1,grid%kfirst:grid%klast) )!for div4 damping

   allocate( grid%cdxde (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping
   allocate( grid%cdxdp (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping
   allocate( grid%cdyde (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping
   allocate( grid%cdydp (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping

   allocate( grid%cdxdiv(jm,grid%kfirst:grid%klast) )!for div4 damping
   allocate( grid%cdydiv(jm,grid%kfirst:grid%klast) )!for div4 damping
   allocate( grid%cdtau4(js2g0:jn1g1,grid%kfirst:grid%klast) )!for div4 damping

   allocate( grid%f0(grid%jfirst-grid%ng_s-1:grid%jlast+grid%ng_d) )
   allocate( grid%fc(js2gc:jn1gc) )

   do j = max(1,grid%jfirst-grid%ng_s-1), min(jm,grid%jlast+grid%ng_d)
      grid%f0(j) = (om+om)*grid%sinp(j)
   end do

   ! Compute coriolis parameter at cell corners.

   if (am_geom_crrct) then
      do j = js2gc, jn1gc
         grid%fc(j) = (om+om)*grid%sine(j)
      end do
   else
      do j = js2gc, jn1gc                    ! Not the issue with ng_c = ng_d
         grid%fc(j) = 0.5_r8*(grid%f0(j) + grid%f0(j-1))
      end do
   end if

   grid%dt0 = 0._r8
   dt5      = 0.5_r8*dt

   grid%rdy   = 1._r8/(ae*grid%dp)
   grid%dtdy  = dt *grid%rdy
   grid%dtdy5 = dt5*grid%rdy
   grid%dydt  = (ae*grid%dp) / dt
   grid%tdy5  = 0.5_r8/grid%dtdy

end subroutine grid_vars_init

!========================================================================================

subroutine dynamics_clean(grid)

   ! Arguments
   type(t_fvdycore_grid), intent(inout)  :: grid

   ! Temporary data structures

    if(associated(GRID%SINLON          )) deallocate(GRID%SINLON)
    if(associated(GRID%COSLON          )) deallocate(GRID%COSLON)
    if(associated(GRID%SINL5           )) deallocate(GRID%SINL5)
    if(associated(GRID%COSL5           )) deallocate(GRID%COSL5)

    if(associated(GRID%ACOSP           )) deallocate(GRID%ACOSP)
    if(associated(GRID%ACOSU           )) deallocate(GRID%ACOSU)
    if(associated(GRID%SINP            )) deallocate(GRID%SINP)
    if(associated(GRID%COSP            )) deallocate(GRID%COSP)
    if(associated(GRID%SINE            )) deallocate(GRID%SINE)
    if(associated(GRID%COSE            )) deallocate(GRID%COSE)
    if(associated(GRID%AK              )) deallocate(GRID%AK)
    if(associated(GRID%BK              )) deallocate(GRID%BK)

    ! cd_core variables

    if(associated( grid%dtdx  )) deallocate(grid%dtdx)
    if(associated( grid%dtdx2 )) deallocate(grid%dtdx2)
    if(associated( grid%dtdx4 )) deallocate(grid%dtdx4)
    if(associated( grid%dtdxe )) deallocate(grid%dtdxe)
    if(associated( grid%dxdt  )) deallocate(grid%dxdt)
    if(associated( grid%dxe   )) deallocate(grid%dxe)
    if(associated( grid%cye   )) deallocate(grid%cye)
    if(associated( grid%dycp  )) deallocate(grid%dycp)
    if(associated( grid%rdxe  )) deallocate(grid%rdxe)
    if(associated( grid%txe5  )) deallocate(grid%txe5)
    if(associated( grid%dtxe5 )) deallocate(grid%dtxe5)
    if(associated( grid%dyce  )) deallocate(grid%dyce)
    if(associated( grid%dx    )) deallocate(grid%dx)
    if(associated( grid%rdx   )) deallocate(grid%rdx)
    if(associated( grid%cy    )) deallocate(grid%cy)

    if(associated( grid%sc    )) deallocate(grid%sc)
    if(associated( grid%se    )) deallocate(grid%se)
    if(associated( grid%dc    )) deallocate(grid%dc)
    if(associated( grid%de    )) deallocate(grid%de)

    if(associated( grid%cdx   )) deallocate(grid%cdx)
    if(associated( grid%cdy   )) deallocate(grid%cdy)
    if(associated( grid%cdx4  )) deallocate(grid%cdx4)  
    if(associated( grid%cdy4  )) deallocate(grid%cdy4)  
    if(associated( grid%cdxde )) deallocate(grid%cdxde) 
    if(associated( grid%cdxdp )) deallocate(grid%cdxdp) 
    if(associated( grid%cdydp )) deallocate(grid%cdydp) 
    if(associated( grid%cdyde )) deallocate(grid%cdyde) 
    if(associated( grid%cdxdiv)) deallocate(grid%cdxdiv)
    if(associated( grid%cdydiv)) deallocate(grid%cdydiv)
    if(associated( grid%cdtau4)) deallocate(grid%cdtau4)

    if(associated( grid%scdiv4)) deallocate(grid%scdiv4)
    if(associated( grid%sediv4)) deallocate(grid%sediv4)
    if(associated( grid%dcdiv4)) deallocate(grid%dcdiv4)
    if(associated( grid%dediv4)) deallocate(grid%dediv4)

    if(associated( grid%f0    )) deallocate(grid%f0)
    if(associated( grid%fc    )) deallocate(grid%fc)

#if defined(SPMD)
   call spmd_vars_clean(grid)
#endif

end subroutine dynamics_clean

!========================================================================================

#if defined(SPMD)
subroutine spmd_vars_clean(grid)

   use parutilitiesmodule, only : parpatternfree

   ! args
   type (T_FVDYCORE_GRID), intent(inout) :: grid
   !-----------------------------------------------------------------------

   if ( grid%twod_decomp == 1 ) then

      ! Clean transposes

      call parpatternfree(grid%commxy, grid%u_to_uxy)
      call parpatternfree(grid%commxy, grid%uxy_to_u)
      call parpatternfree(grid%commxy, grid%v_to_vxy)
      call parpatternfree(grid%commxy, grid%vxy_to_v)
      call parpatternfree(grid%commxy, grid%ijk_yz_to_xy)
      call parpatternfree(grid%commxy, grid%ijk_xy_to_yz)
      call parpatternfree(grid%commxy, grid%ikj_xy_to_yz)
      call parpatternfree(grid%commxy, grid%ikj_yz_to_xy)
      call parpatternfree(grid%commxy, grid%pe_to_pexy)
      call parpatternfree(grid%commxy, grid%pexy_to_pe)
      call parpatternfree(grid%commxy, grid%pt_to_ptxy)
      call parpatternfree(grid%commxy, grid%ptxy_to_pt)
      call parpatternfree(grid%commxy, grid%r4_xy_to_yz)
      call parpatternfree(grid%commxy, grid%r4_yz_to_xy)
      call parpatternfree(grid%commxy, grid%pkxy_to_pkc)
      call parpatternfree(grid%commxy, grid%xy2d_to_yz2d)
      call parpatternfree(grid%commxy, grid%yz2d_to_xy2d)
   endif

end subroutine spmd_vars_clean
#endif

!========================================================================================

end module dynamics_vars

