module dyn_grid

    ! Purpose: Definition of dynamics computational grid.
    !
    ! Method: Variables are private; interface routines used to extract
    !         information for use in user code. Global column index range
    !         defined using full (unreduced) grid.
    !
    ! Entry points:
    !      get_block_bounds_d       get first and last indices in global
    !                               block ordering
    !      get_block_gcol_d         get column indices for given block
    !      get_block_gcol_cnt_d     get number of columns in given block
    !      get_block_lvl_cnt_d      get number of vertical levels in column
    !      get_block_levels_d       get vertical levels in column
    !      get_gcol_block_d         get global block indices and local columns
    !                               index for given global column index
    !      get_gcol_block_cnt_d     get number of blocks containing data
    !                               from a given global column index
    !      get_block_owner_d        get process "owning" given block
    !      get_horiz_grid_d         get horizontal grid coordinates
    !      get_horiz_grid_dim_d     get horizontal dimensions of dynamics grid
    !      dyn_grid_get_elem_coords get coordinates of a specified block element
    !                               of the dynamics grid
    !      dyn_grid_get_colndx      get block/column and MPI process indices
    !                               corresponding to a specified global column index
    !
    ! Author: Jim Edwards and Patrick Worley

    use physconst,        only: rearth
    use cam_abortutils,   only: endrun
    use cam_grid_support, only: iMap
    use cam_logfile,      only: iulog
    use pio,              only: file_desc_t
    use shr_kind_mod,     only: r8 => shr_kind_r8
    use spmd_utils,       only: masterprocid,mpicom
!    use fv_grid_tools_mod,only: gindex,bindex,lindex
    use fv_arrays_mod,    only: fv_atmos_type
    use mpp_domains_mod,    only: domain2d,mpp_domains_set_stack_size
    use time_manager_mod,   only: time_type, get_time
    use fms_mod,            only: open_namelist_file, file_exist, check_nml_error,  &
                                 error_mesg, fms_init, fms_end, &
                                 write_version_number, uppercase, FATAL, NOTE
    use fv_control_mod,     only: ngrids,fv_init

    use dimensions_mod,     only: npx, npy, npz, ncnst, pnats, dnats,nq
    use mpp_mod,            only: mpp_pe, mpp_root_pe
    use fv_mp_mod,          only: mp_bcst
    use constituents,       only: pcnst
    implicit none
    private
    save

    ! The FV3 dynamics grids
    real(r8), parameter :: pi = 3.1415926535897931_r8
    real(r8), parameter :: rad2deg = 180._r8/pi
    integer, parameter, public :: dyn_decomp   = 101
    integer, parameter, public :: dyn_decomp_ew = 102
    integer, parameter, public :: dyn_decomp_ns = 103
    integer, parameter, public :: dyn_decomp_z = 104
    integer, parameter, public :: dyn_decomp_ew_rst = 105
    integer, parameter, public :: dyn_decomp_ns_rst = 106
    integer, public ::  ntiles   = -999
    integer, public ::  nest_pes = 0
    integer, public ::  p_split  = 1

    integer, parameter, public :: ptimelevels = 2  ! number of time levels in the dycore

    !These are convenience variables for local use only, and are set to values in Atm%
    real(r8)    :: dt_atmos
    real(r8)    :: zvir
    integer :: sec, seconds, days
    logical :: cold_start = .false.       ! read in initial condition
    
    integer, dimension(:), allocatable :: id_tracerdt_dyn
    integer :: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel  !condensate species
    
    integer :: mytile = 1
    integer, allocatable :: pelist(:)
    logical, allocatable, public :: grids_on_this_pe(:)
    type(fv_atmos_type), allocatable, target :: Atm(:)
    
    integer :: id_udt_dyn, id_vdt_dyn
    
    real(r8), parameter:: w0_big = 60.  ! to prevent negative w-tracer diffusion
    real(r8), allocatable :: block_extents_g(:,:)
!-----------------------------------------------------------------------
! zlj, 2014.10.13
! Calculate Global Index

integer, allocatable, target, dimension(:,:) :: mygindex,mygindex_ew,mygindex_ns
integer, allocatable, target, dimension(:,:) :: mygindexdups,mygindexdups_ew,mygindexdups_ns
integer, allocatable, target, dimension(:,:) :: mygindex_tiles,mygindex_tiles_ew,mygindex_tiles_ns
integer, allocatable, target, dimension(:,:) :: mylindex,mylindex_ew,mylindex_ns
integer                                      :: mybindex,mybindex_ew,mybindex_ns
real(r8), allocatable, target, dimension(:,:,:) ::   locidx_g
real(r8), allocatable, target, dimension(:,:,:) ::   blkidx_g
real(r8), allocatable, target, dimension(:,:,:) ::   gindex_g
real(r8), allocatable, target, dimension(:,:) ::   gblidx_g

integer, public :: uniqpts_glob = 0     ! number of dynamics columns
integer, public :: uniqpts_glob_ew = 0     ! number of dynamics columns for Dgrid ew
integer, public :: uniqpts_glob_ns = 0     ! number of dynamics columns for Dgrid ew
integer, public :: uniqpts_loc = 0     ! number of dynamics columns
integer, public :: uniqpts_loc_ew = 0     ! number of dynamics columns for Dgrid ew
integer, public :: uniqpts_loc_ns = 0     ! number of dynamics columns for Dgrid ew

!real(r8), allocatable, dimension(:,:,:) :: grid_ew, grid_ns
real(r8), pointer, dimension(:,:,:) :: grid_ew, grid_ns

! Calculate Global Index
! zlj, 2014.10.13

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! zlj, 2014.10.13
! Calculate Global Index

public :: mygindex,mygindex_ew,mygindex_ns
public :: mygindexdups, mygindexdups_ew,mygindexdups_ns
public :: mygindex_tiles, mygindex_tiles_ew,mygindex_tiles_ns
public :: mylindex,mylindex_ew,mylindex_ns
public :: mybindex,mybindex_ew,mybindex_ns
public :: grid_ew,grid_ns
! Calculate Global Index
! zlj, 2014.10.13
!-----------------------------------------------------------------------

    public :: get_block_bounds_d, get_block_gcol_d, get_block_gcol_cnt_d, &
              get_block_lvl_cnt_d, get_block_levels_d, get_block_owner_d, &
              get_gcol_block_d, get_gcol_block_cnt_d, &
              get_horiz_grid_d, get_horiz_grid_dim_d
    public :: get_dyn_grid_parm, get_block_ldof_d
    public :: get_dyn_grid_parm_real1d, get_dyn_grid_parm_real2d
    public :: dyn_grid_get_elem_coords
    public :: dyn_grid_get_colndx
    public :: define_cam_grids
    public :: physgrid_copy_attributes_d
    public :: dyn_grid_init
    public Atm, mytile

    real(r8), public, pointer :: w(:) => null()        ! weights
    real(r8), pointer :: clat(:) => null()     ! model latitudes (radians)
    real(r8), pointer :: clon(:,:) => null()   ! model longitudes (radians)
    real(r8), pointer :: latdeg(:) => null()   ! model latitudes (degrees)
    real(r8), pointer :: londeg(:,:) => null() ! model longitudes (degrees)

!=======================================================================
contains
!=======================================================================

subroutine dyn_grid_init()

   ! Initialize FV grid, decomposition

   use hycoef,             only: hycoef_init, hyai, hybi, hypi, hypm, nprlev
   use ref_pres,           only: ref_pres_init
   use pio,                only: file_desc_t, pio_seterrorhandling, pio_bcast_error, &
                                 pio_internal_error, pio_noerr, pio_inq_dimid,       &
                                 pio_inq_dimlen
   use cam_initfiles,      only: initial_file_get_id
   use mpp_mod,            only: mpp_init
   use namelist_utils,     only: find_group_name
   use units,              only: getunit, freeunit
   use block_control_mod,  only: block_control_type, define_blocks_packed
   use constants_mod,      only: rdgas, rvgas
   use constants_mod,      only: constants_init
   use fms_mod,            only: close_file
   use fms_io_mod,         only: set_domain, nullify_domain
   use fv_mp_mod,          only: switch_current_Atm,mp_gather, mp_bcst
   use fv_restart_mod,     only: fv_restart
   use memutils_mod,       only: print_memuse_stats
   use mpp_mod,            only: mpp_init, mpp_npes, mpp_get_current_pelist,mpp_gather
   use pmgrid,             only: plev
   use sat_vapor_pres_mod, only: sat_vapor_pres_init
   use time_manager,       only: get_curr_date,get_step_size
   use tracer_manager_mod, only: get_tracer_index
   use xgrid_mod,          only: grid_box_type
   use fv_eta_mod,          only: set_eta

   ! Local variables

   integer :: i, k, lat

   real(r8) :: dt
   real(r8) :: dp
   real(r8) :: sum

   character(len=*), parameter :: sub='dyn_grid_init'

   integer :: m
   integer :: unit, io
   character(len = 256) :: err_msg
   character(len = 9) :: month
   character(len = 17) :: calendar,fv_calendar = 'NOLEAP           '
   
   type(file_desc_t),      pointer :: fh_ini
   real(r8) :: dt_atmos_real = 0.
   real(r8),allocatable :: rtmp(:)
   real(r8) :: block_extents(5)
   integer, allocatable :: be_size(:)
! ----- namelist -----
   integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
   character(len=17) :: fvcalendar = '                 '
   integer :: months=0, days=0, hours=0, minutes=0, seconds=0
   integer :: dt_atmos = 0
   integer :: dt_ocean = 0
   integer :: restart_days = 0
   integer :: restart_secs = 0
   integer :: atmos_nthreads = 1
   logical :: memuse_verbose = .false.
   logical :: use_hyper_thread = .false.
   logical :: debug_affinity = .false.
   integer :: ncores_per_node = 0
   integer :: restart_interval(6) = 0
   integer :: memuse_interval = 0

   integer :: unitn               ! File unit number
   integer :: ierr                ! Error code 
   integer :: n                   ! index
   integer :: nlat,nlon,mlat,mlon
namelist /main_nml/ current_date, fv_calendar, &
     months, days, hours, minutes, seconds,  &
     dt_atmos,atmos_nthreads, memuse_verbose, & 
     use_hyper_thread, ncores_per_node, debug_affinity, &
     restart_secs, restart_days, restart_interval, &
     memuse_interval

!-----------------------------------------------------------------------
integer :: ncolid,ncollen
integer :: blocksize    = -1
logical :: chksum_debug = .false.
logical :: dycore_only  = .false.
logical :: debug        = .false.
logical :: sync         = .false.
logical, save :: block_message = .true.

namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
type (block_control_type), target   :: Atm_block
integer :: is,isc,ie,iec,js,jsc,je,jec,npes,tsize,ssize

   !-----------------------------------------------------------------------

   ! Get file handle for initial file and first consistency check
   fh_ini => initial_file_get_id()

   call pio_seterrorhandling(fh_ini, pio_bcast_error)
   ierr = pio_inq_dimid(fh_ini, 'ncol', ncolid)
   call pio_seterrorhandling(fh_ini, pio_internal_error)

   if (ierr /= pio_noerr) then
      call endrun(sub//': ERROR: initial dataset not on unstructured grid')
   else      
      ierr = pio_inq_dimlen(fh_ini, ncolid, ncollen)
   end if

   !-----------------------------------------------------------------------
   !  from couple_main initialize atm structure - initializes fv3 grid
   !-----------------------------------------------------------------------

   !!********From Coupler_main and coupler_init.F90
   call fms_init(mpicom)
   call mpp_init()
   
   call fms_init
   call constants_init
   call sat_vapor_pres_init

   ! =======================
   ! Read namelist variables
   ! =======================
   
      unitn=getunit()
      open( unitn, file='input.nml', status='old' )
      
      ! Look for fv3_nl group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'main_nml', status=ierr)
      if ( ierr == 0 ) then
         read (unitn,main_nml,iostat=ierr)
         if (ierr /= 0) then
            call endrun( 'dyn_grid_init:: namelist fv3_nl read returns an'// &
                 ' end of file or end of record condition' )
         end if
      end if
      close( unitn )
      call freeunit( unitn )
      
      ! Tests for valid input, and set some defaults
      
      if (dt_atmos .eq. 0) then
         call error_mesg('subroutine dyn_grid_init', 'dt_atmos has not been specified', FATAL)
      end if
 
!-----------------------------------------------------------------------
! initialize atmospheric model -----

   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)

   zvir = rvgas/rdgas - 1.

!---- compute physics/atmos time step in seconds ----

   dt_atmos_real = dt_atmos

!----- initialize FV dynamical core -----

   call fv_init( Atm, dt_atmos_real, grids_on_this_pe, p_split)  ! allocates Atm components

   do n=1,ngrids
      if (grids_on_this_pe(n)) mytile = n
   enddo

!----- write version and namelist to log file -----
   call write_version_number ( version, tagname )

   n = mytile
   call switch_current_Atm(Atm(n)) 

!---------- inline atmosphere_init -------
   IF ( file_exist('input.nml')) THEN
#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
#else
      unitn = open_namelist_file ( )
      ierr=1
      do while (ierr /= 0)
         read  (unitn, nml=atmos_model_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'atmos_model_nml')
      enddo
 10     call close_file (unitn)
#endif
   endif

!!  set up dimensions_mod convenience variables. 

   is = Atm(mytile)%bd%is
   ie = Atm(mytile)%bd%ie
   js = Atm(mytile)%bd%js
   je = Atm(mytile)%bd%je
   isc = Atm(mytile)%bd%isc
   iec = Atm(mytile)%bd%iec
   jsc = Atm(mytile)%bd%jsc
   jec = Atm(mytile)%bd%jec
   ncnst  = Atm(mytile)%flagstruct%ncnst
   npx   = Atm(mytile)%flagstruct%npx
   npy   = Atm(mytile)%flagstruct%npy
   npz    = Atm(mytile)%flagstruct%npz
   nq     = Atm(mytile)%flagstruct%ncnst - Atm(mytile)%flagstruct%pnats
   pnats  = Atm(mytile)%flagstruct%pnats
   dnats = Atm(mytile)%flagstruct%dnats ! Number of non-advected consituents (as seen by dynamics)

!-----------------------------------------------------------------------
!--- get block extents for each task/pe
!-----------------------------------------------------------------------
   npes=mpp_npes()
   allocate(block_extents_g(5,npes))
   allocate(rtmp(5*npes))
   allocate(be_size(npes))
   ssize=5
   be_size(:)=ssize
   block_extents(1)=is;block_extents(2)=ie;block_extents(3)=js;block_extents(4)=je;block_extents(5)=Atm(mytile)%tile
   call mpp_gather(block_extents,ssize,rtmp,be_size)
   call mp_bcst(rtmp,5*npes)
   block_extents_g=reshape(rtmp,(/5,npes/))

   deallocate(rtmp)
   deallocate(be_size)
   call define_blocks_packed ('atmos_model', Atm_block,  isc, iec, &
                              jsc, jec, npz, &
                              blocksize, block_message)

 ! Initialize hybrid coordinate arrays
 call hycoef_init(fh_ini)

 ! Initialize reference pressures
 call ref_pres_init(hypi, hypm, nprlev)

   ! Hybrid coordinate info for FV grid object
   Atm(mytile)%ks = plev
   do k = 1, plev+1
      Atm(mytile)%ak(k) = hyai(k) * 1.e5_r8
      Atm(mytile)%bk(k) = hybi(k)
      if ( Atm(mytile)%bk(k) == 0._r8)  Atm(mytile)%ks = k-1
   end do
   Atm(mytile)%ptop = Atm(mytile)%ak(1)

   ! Define the CAM grids
   call define_cam_grids(Atm)

!  This is adjusted via fms_nml variable domain_stack_size (for cam namelist use fv3_domain_stack_size)
!   call mpp_domains_set_stack_size(3000000)

end subroutine dyn_grid_init

!=======================================================================

subroutine get_block_bounds_d(block_first, block_last)

  use spmd_utils,    only : npes

    implicit none

    integer, intent(out) :: block_first  ! first (global) index used for blocks
    integer, intent(out) :: block_last   ! last (global) index used for blocks

    block_first = 1
    block_last = npes

end subroutine get_block_bounds_d

!=======================================================================

subroutine get_block_gcol_d(blockid, size, cdex)

    implicit none

    integer, intent(in) :: blockid      ! global block id
    integer, intent(in) :: size         ! array size
    integer, intent(out):: cdex(size)   ! global column indices

    integer :: i, j, n,is,ie,js,je,tile

    is=block_extents_g(1,blockid)
    ie=block_extents_g(2,blockid)
    js=block_extents_g(3,blockid)
    je=block_extents_g(4,blockid)
    tile=block_extents_g(5,blockid)

    if (size .ne. (ie - is + 1) * (je - js + 1)) then
        call endrun ('get_block_gcol_d: block sizes are not consistent.')
    end if

    do j = js, je
        do i = is, ie
            n=locidx_g(i,j,tile)
            cdex(n) = gindex_g(i,j,tile)
        end do
    end do

end subroutine get_block_gcol_d

!=======================================================================

integer function get_block_gcol_cnt_d(blockid)

    implicit none

    integer, intent(in) :: blockid

    get_block_gcol_cnt_d=count(blkidx_g.eq.blockid)

end function get_block_gcol_cnt_d

!=======================================================================

integer function get_block_lvl_cnt_d(blockid, bcid)

    use pmgrid, only: plevp

    implicit none

    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid     ! column index within block

    get_block_lvl_cnt_d = plevp

end function get_block_lvl_cnt_d

!=======================================================================

subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

    use pmgrid, only: plev

    implicit none

    integer, intent(in)  :: blockid        ! global block id
    integer, intent(in)  :: bcid           ! column index within block
    integer, intent(in)  :: lvlsiz         ! dimension of levels array
    integer, intent(out) :: levels(lvlsiz) ! levels indices for block

    integer :: k

    if (lvlsiz < plev + 1) then
        write(iulog,*)'GET_BLOCK_LEVELS_D: levels array not large enough (',lvlsiz,' < ',plev + 1,')'
        call endrun
    else
        do k = 0, plev
            levels(k+1) = k
        enddo
        do k = plev + 2, lvlsiz
            levels(k) = -1
        enddo
    end if

end subroutine get_block_levels_d

!=======================================================================

integer function get_block_owner_d(blockid)

    implicit none

    integer, intent(in) :: blockid  ! global block id

    get_block_owner_d = blockid - 1

end function get_block_owner_d

!=======================================================================

subroutine get_gcol_block_d(gcol, cnt, blockid, bcid, localblockid)
 
    use dimensions_mod,     only: npx, npy
    implicit none
    integer, intent(in)  :: gcol         ! global column index
    integer, intent(in)  :: cnt          ! size of blockid and bcid arrays
    integer, intent(out) :: blockid(cnt) ! block index
    integer, intent(out) :: bcid(cnt)    ! column index within block
    integer, intent(out), optional :: localblockid(cnt)

    integer :: ind(2),tot
    integer :: tmpblkid(npx-1,npy-1,6),tmplcid(npx-1,npy-1,6)

    if (cnt .ne. 1) then
       call endrun ('get_gcol_block_d: cnt is not equal to 1:.')
    end if
    tot=(npx-1)*(npy-1)*6
    if (gcol.lt.1.or.gcol.gt.tot) then
       call endrun ('get_gcol_block_d: global column number is out of bounds')
    else
       blockid(1) = gblidx_g(gcol,1)
       bcid(1) = gblidx_g(gcol,2)
    end if

    if (present(localblockid)) then
        localblockid(cnt) = 1
    end if

end subroutine get_gcol_block_d

!=======================================================================

integer function get_gcol_block_cnt_d(gcol)

    implicit none

    integer, intent(in) :: gcol     ! global column index

    get_gcol_block_cnt_d = 1

end function get_gcol_block_cnt_d

!=======================================================================

subroutine get_horiz_grid_d(nxy, clat_d_out, clon_d_out, area_d_out, wght_d_out, lat_d_out, lon_d_out)

    use fv_mp_mod,         only: mp_gather, mp_bcst
    implicit none

    integer, intent(in)   :: nxy                       ! array sizes
    real(r8), intent(out), optional :: clat_d_out(nxy) ! column latitudes
    real(r8), intent(out), optional :: clon_d_out(nxy) ! column longitudes
    real(r8), intent(out), optional :: area_d_out(nxy) ! column surface area
    real(r8), intent(out), optional :: wght_d_out(nxy) ! column integration
    real(r8), intent(out), optional :: lat_d_out(nxy)  ! column degree latitudes
    real(r8), intent(out), optional :: lon_d_out(nxy)  ! column degree longitudes

    integer :: i, j, k

    real(r8) :: clon_d_mpi(nxy)
    real(r8) :: clat_d_mpi(nxy)
    real(r8) :: area_d_mpi(nxy)
    real(r8) :: wght_d_mpi(nxy)

    real(r8), allocatable :: clon_d(:,:,:)
    real(r8), allocatable :: clat_d(:,:,:)
    real(r8), allocatable :: area_d(:,:,:)
    real(r8),  allocatable :: indx_d(:,:,:)

    real(r8), pointer, dimension(:,:,:) :: agrid
    real(r8), pointer, dimension(:,:)   :: area
    integer,  pointer                   :: ntiles_g,tile
    integer                             :: is,ie,js,je
    integer                             :: masterproc,gid
!
    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    masterproc = mpp_root_pe()
    gid = mpp_pe()

    area  => Atm(mytile)%gridstruct%area_64
    agrid => Atm(mytile)%gridstruct%agrid_64
    ntiles_g => Atm(mytile)%gridstruct%ntiles_g
    tile=>  Atm(mytile)%tile

    clon_d_mpi = -999.0_r8
    clat_d_mpi = -999.0_r8
    area_d_mpi = -999.0_r8
    wght_d_mpi = -999.0_r8

    allocate(clon_d(npx-1, npy-1, ntiles_g))
    allocate(clat_d(npx-1, npy-1, ntiles_g))
    allocate(area_d(npx-1, npy-1, ntiles_g))
    allocate(indx_d(npx-1, npy-1, ntiles_g))

    do j = js, je
        do i = is, ie
            clon_d(i,j,tile) = agrid (i,j,1)
            clat_d(i,j,tile) = agrid (i,j,2)
            area_d(i,j,tile) = area  (i,j) / (rearth * rearth)
            indx_d(i,j,tile) = mygindex  (i,j)
        end do
    end do

    call mp_gather(indx_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
    if (gid == masterproc) then
        if (maxval(indx_d) .gt. nxy .or. minval(indx_d) .lt. 1) then
            write(iulog,*)'GET_HORIZ_GRID_D: indx_d is out of bound. nxy = ',nxy
            call endrun
        end if
    end if

    if (present(clat_d_out)) then
        call mp_gather(clat_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (gid == masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                        clat_d_mpi(indx_d(i,j,k)) = clat_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clat_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clat_d_mpi .eq. -999.0_r8)) then
            clat_d_out = clat_d_mpi
        else
            call endrun('clat_d_mpi is not full filled')
        end if
    end if

    if (present(clon_d_out)) then
        call mp_gather(clon_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (gid == masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                        clon_d_mpi(indx_d(i,j,k)) = clon_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clon_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clon_d_mpi .eq. -999.0_r8)) then
            clon_d_out = clon_d_mpi
        else
            call endrun('clon_d_mpi is not full filled')
        end if
    end if

    if (present(lat_d_out)) then
        call mp_gather(clat_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (gid == masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                        clat_d_mpi(indx_d(i,j,k)) = clat_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clat_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clat_d_mpi .eq. -999.0_r8)) then
            lat_d_out = clat_d_mpi * rad2deg
        else
            call endrun('clat_d_mpi is not full filled')
        end if
    end if

    if (present(lon_d_out)) then
        call mp_gather(clon_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (gid == masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                        clon_d_mpi(indx_d(i,j,k)) = clon_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clon_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clon_d_mpi .eq. -999.0_r8)) then
            lon_d_out = clon_d_mpi * rad2deg
        else
            call endrun('clon_d_mpi is not full filled')
        end if
    end if

    if (present(area_d_out)) then
        call mp_gather(area_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (gid == masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                        area_d_mpi(indx_d(i,j,k)) = area_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(area_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(area_d_mpi .eq. -999.0_r8)) then
            area_d_out = area_d_mpi
        else
            call endrun('area_d_mpi is not full filled')
        end if
    end if

    if (present(wght_d_out)) then
        call mp_gather(area_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (gid == masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                        wght_d_mpi(indx_d(i,j,k)) = area_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(wght_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(wght_d_mpi .eq. -999.0_r8)) then
            wght_d_out = wght_d_mpi
        else
            call endrun('wght_d_mpi is not full filled')
        end if
    end if

    if (present(clat_d_out) .and. present(clon_d_out)) then
        if (.not. associated(clat)) then
            allocate(clat(nxy), clon(nxy,1))
        end if
        clat(:  ) = clat_d_out
        clon(:,1) = clon_d_out
    end if

    deallocate(clon_d)
    deallocate(clat_d)
    deallocate(area_d)
    deallocate(indx_d)

end subroutine get_horiz_grid_d

!=======================================================================

subroutine get_horiz_grid_dim_d(hdim1_d, hdim2_d)
    use spmd_utils,    only : npes
    implicit none

    integer, intent(out) :: hdim1_d                 ! first horizontal dimension
    integer, intent(out), optional :: hdim2_d       ! second horizontal dimension
    integer is,ie,js,je

    hdim1_d = count(gindex_g.gt.0)
    if (present(hdim2_d)) hdim2_d = 1

end subroutine get_horiz_grid_dim_d

!=======================================================================


subroutine define_cam_grids(Atm)
  
  use cam_grid_support,  only: horiz_coord_t, horiz_coord_create, cam_grid_get_local_size
  use cam_grid_support,  only: cam_grid_register, cam_grid_attribute_register
  use fv_grid_utils_mod, only: mid_pt_sphere
  use mpp_mod,           only: mpp_pe, mpp_npes
  use fv_mp_mod,         only: mp_gather, mp_bcst, mp_barrier
  use spmd_utils,        only : npes
  use physconst,         only: rearth
  implicit none
  
  type(fv_atmos_type), target, intent(in) :: Atm(:)


  type(horiz_coord_t), pointer :: lat_coord
  type(horiz_coord_t), pointer :: lon_coord
  
  integer(iMap), pointer :: grid_map(:,:),grid_map_dups(:,:)

  integer :: n, i, j, mapind,is,ie,js,je,isd,ied,jsd,jed,ierr,tile,nregions
  real(r8), pointer, dimension(:,:,:) :: agrid
  real(r8), pointer, dimension(:,:,:) :: grid
  real(r8), pointer, dimension(:,:)   :: area
  real(r8), pointer :: area_ffsl(:)                   !fv3 cell centered grid area in sq radians
  real(r8), pointer :: pelon_deg(:)
  real(r8), pointer :: pelat_deg(:)
  real(r8), pointer :: pelon_deg_ew(:)
  logical, pointer  :: mask_ew(:)
  logical, pointer  :: mask_ns(:)
  real(r8), pointer :: pelat_deg_ew(:)
  real(r8), pointer :: pelon_deg_ns(:)
  real(r8), pointer :: pelat_deg_ns(:)
  real(r8), pointer :: coord_map_ns(:)
  real(r8), pointer :: coord_map_ew(:)
  real(r8)               :: pt(2),lonrad,latrad
  integer(iMap), pointer :: pemap(:)
  integer(iMap), pointer :: pemap_ew(:)
  integer(iMap), pointer :: pemap_ns(:)
  integer(iMap), pointer :: pemap_dups_ew(:)
  integer(iMap), pointer :: pemap_dups_ns(:)
  character(12)          :: filename
  character(2)          :: cgid
  integer               :: masterproc,gid
  integer :: ncols_glob = 0     ! number of dynamics columns
  integer :: ncols_glob_ew = 0     ! number of dynamics columns for Dgrid ew
  integer :: ncols_glob_ns = 0     ! number of dynamics columns for Dgrid ns
  integer :: ncols_loc = 0     ! number of dynamics columns
  integer :: ncols_loc_ew = 0     ! number of dynamics columns for Dgrid ew
  integer :: ncols_loc_ns = 0     ! number of dynamics columns for Dgrid ew
  integer :: glob_pts_rst         ! number of dynamics columns for Dgrid restarts (includes dup points.)
  integer :: nx,ny

  area  => Atm(mytile)%gridstruct%area_64
  agrid => Atm(mytile)%gridstruct%agrid_64
  grid => Atm(mytile)%gridstruct%grid_64
  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  isd = Atm(mytile)%bd%isd
  ied = Atm(mytile)%bd%ied
  jsd = Atm(mytile)%bd%jsd
  jed = Atm(mytile)%bd%jed
  tile = Atm(mytile)%tile
  nregions = Atm(mytile)%gridstruct%ntiles_g
  ncols_glob= (npx-1)*(npy-1)*nregions
  ncols_glob_ew= (npx)*(npy-1)*nregions
  ncols_glob_ns= (npx-1)*(npy)*nregions

  allocate(area_ffsl((ie-is+1)*(je-js+1)))
  allocate(grid_ew(isd:ied+1,jsd:jed,2))
  allocate(grid_ns(isd:ied,jsd:jed+1,2))
  allocate(pelon_deg((ie-is+1)*(je-js+1)))
  allocate(mask_ns((ie-is+1)*(je-js+2)))
  allocate(pelon_deg_ns((ie-is+1)*(je-js+2)))
  allocate(mask_ew((ie-is+2)*(je-js+1)))
  allocate(pelon_deg_ew((ie-is+2)*(je-js+1)))
  allocate(pelat_deg((ie-is+1)*(je-js+1)))
  allocate(pelat_deg_ew((ie-is+2)*(je-js+1)))
  allocate(coord_map_ew((ie-is+2)*(je-js+1)))
  allocate(pelat_deg_ns((ie-is+1)*(je-js+2)))
  allocate(coord_map_ns((ie-is+1)*(je-js+2)))
  allocate(pemap((ie-is+1)*(je-js+1)))
  allocate(pemap_ew((ie-is+2)*(je-js+1)))
  allocate(pemap_ns((ie-is+1)*(je-js+2)))
  allocate(pemap_dups_ew((ie-is+2)*(je-js+1)))
  allocate(pemap_dups_ns((ie-is+1)*(je-js+2)))


  do j=jsd,jed
     do i=isd,ied+1
        call mid_pt_sphere(grid(i,  j,1:2), grid(i,  j+1,1:2), grid_ew(i,j,:))
     end do
  end do
  
  do j=jsd,jed+1
     do i=isd,ied
        call mid_pt_sphere(grid(i,j  ,1:2), grid(i+1,j  ,1:2), grid_ns(i,j,:))
     end do
  end do

  allocate(locidx_g(npx-1,npy-1,nregions))
  allocate(blkidx_g(npx-1,npy-1,nregions))
  allocate(gindex_g(npx-1,npy-1,nregions))
  allocate(gblidx_g( (npx-1)*(npy-1)*nregions,2))
  allocate(mygindex(is:ie,js:je))
  allocate(mygindex_ew(is:ie+1,js:je))
  allocate(mygindex_ns(is:ie,js:je+1))

  allocate(mygindexdups(is:ie,js:je))
  allocate(mygindexdups_ew(is:ie+1,js:je))
  allocate(mygindexdups_ns(is:ie,js:je+1))

  allocate(mygindex_tiles(is:ie,js:je))
  allocate(mygindex_tiles_ew(is:ie+1,js:je))
  allocate(mygindex_tiles_ns(is:ie,js:je+1))

  allocate(mylindex(is:ie,js:je))
  allocate(mylindex_ew(is:ie+1,js:je))
  allocate(mylindex_ns(is:ie,js:je+1))

  masterproc = mpp_root_pe()
  gid = mpp_pe()

  call calc_global_indexjt(is,ie  ,js,je   ,npx-1 ,npy-1 ,tile, nregions, mygindex   ,mylindex   ,mybindex   ,mygindexdups, mygindex_tiles, 'gridew.txt' ,uniqpts_loc    ,uniqpts_glob)
  call calc_global_indexjt(is,ie+1,js,je   ,npx   ,npy-1 ,tile, nregions, mygindex_ew,mylindex_ew,mybindex_ew,mygindexdups_ew, mygindex_tiles_ew, 'gridew.txt' ,uniqpts_loc_ew ,uniqpts_glob_ew)
  call calc_global_indexjt(is,ie  ,js,je+1 ,npx-1 ,npy   ,tile, nregions, mygindex_ns,mylindex_ns,mybindex_ns,mygindexdups_ns, mygindex_tiles_ns, 'gridns.txt' ,uniqpts_loc_ns ,uniqpts_glob_ns)


  blkidx_g(is:ie,js:je,tile)=mybindex
  locidx_g(is:ie,js:je,tile)=mylindex(is:ie,js:je)
  gindex_g(is:ie,js:je,tile)=mygindex(is:ie,js:je)

  nx=npx-1
  ny=npy-1
  call mp_gather(locidx_g, is, ie, js, je, nx ,ny, nregions)
  call mp_gather(blkidx_g, is, ie, js, je, nx ,ny, nregions)
  call mp_gather(gindex_g, is, ie, js, je, nx, ny, nregions)
  call mp_bcst(locidx_g, nx, ny, nregions)
  call mp_bcst(blkidx_g, nx, ny, nregions)
  call mp_bcst(gindex_g, nx, ny, nregions)

  do n = 1, nregions
     do j=1,nx
        do i=1,ny
           gblidx_g(gindex_g(i,j,n),1)=blkidx_g(i,j,n)
           gblidx_g(gindex_g(i,j,n),2)=locidx_g(i,j,n)
        end do
     end do
  end do

  !-----------------------
  ! Create FFSL grid object
  !-----------------------
  
  ! Calculate the mapping between FFSL points and file order (tile1 thru tile6)
  do j = js, je
     do i = is, ie
        n = mylindex(i,j)
        pelon_deg(n) = agrid(i,j,1) * rad2deg
        pelat_deg(n) = agrid(i,j,2) * rad2deg
        area_ffsl(n) = area(i,j)/(rearth*rearth)
        pemap(n)     = mygindex(i,j)
     end do
  end do

  pelon_deg_ew = 0
  pelat_deg_ew = 0
  do j = js, je
     do i = is, ie+1
        n=mylindex_ew(i,j)
        lonrad=grid_ew(i,j,1)
        latrad=grid_ew(i,j,2)
        pelon_deg_ew(n) = lonrad * rad2deg
        pelat_deg_ew(n) = latrad * rad2deg
        pemap_ew(n)     = mygindex_ew(i,j)
        pemap_dups_ew(n) = mygindexdups_ew(i,j)
     end do
  end do
  mask_ew=pemap_ew.ne.0
  pelon_deg_ns = 0
  pelat_deg_ns = 0

  do j = js, je+1
     do i = is, ie
        n=mylindex_ns(i,j)
        lonrad=grid_ns(i,j,1)
        latrad=grid_ns(i,j,2)
        pelon_deg_ns(n) = lonrad * rad2deg
        pelat_deg_ns(n) = latrad * rad2deg
        pemap_ns(n)     = mygindex_ns(i,j)
        pemap_dups_ns(n) = mygindexdups_ns(i,j)
     end do
  end do
  mask_ns=pemap_ns.ne.0
    
  allocate(grid_map(3, (ie-is+1)*(je-js+1)))
  grid_map = 0
  do j = js, je
     do i = is, ie
        mapind=mylindex(i,j)
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap(mapind)
     end do
  end do
  
  lat_coord => horiz_coord_create('lat_d', 'ncol_d', uniqpts_glob, 'latitude',      &
       'degrees_north', 1, size(pelat_deg), pelat_deg, map=pemap)
  lon_coord => horiz_coord_create('lon_d', 'ncol_d', uniqpts_glob, 'longitude',     &
       'degrees_east', 1, size(pelon_deg), pelon_deg, map=pemap)

  call cam_grid_register('FFSL', dyn_decomp, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL', 'cell', '', 1)
  call cam_grid_attribute_register('FFSL', 'area_d', 'FFSL grid areas', &
         'ncol_d', area_ffsl, map=pemap)         
  nullify(grid_map)
  nullify(lat_coord)
  nullify(lon_coord)
  nullify(area_ffsl)

  lat_coord => horiz_coord_create('lat_d_ew', 'ncol_d_ew', uniqpts_glob_ew, 'latitude',      &
       'degrees_north', 1, size(pelat_deg_ew), pelat_deg_ew, map=pemap_ew)
  lon_coord => horiz_coord_create('lon_d_ew', 'ncol_d_ew',  uniqpts_glob_ew, 'longitude',     &
       'degrees_east', 1, size(pelon_deg_ew), pelon_deg_ew, map=pemap_ew)
  
  allocate(grid_map(3, (ie-is+2)*(je-js+1)))
  grid_map = 0
  do j = js, je
     do i = is, ie+1
        mapind = mylindex_ew(i,j)
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap_ew(mapind)
     end do
  end do
  
  call cam_grid_register('FFSL_EW', dyn_decomp_ew, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL_EW', 'cell', '', 1)
  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  nullify(lat_coord)         ! Belongs to grid
  nullify(lon_coord)         ! Belongs to grid
  
! create ew grid of all pts for restarts (number of tile points for ew ns grids is npx-1*(npy)*6 or npx*(npy-1)*6, both are equivalent since npx=npy

  lat_coord => horiz_coord_create('lat_d_ew_rst', 'ncol_d_ew_rst', uniqpts_glob_ew, 'latitude',      &
       'degrees_north', 1, size(pelat_deg_ew), pelat_deg_ew, map=pemap_dups_ew)
  lon_coord => horiz_coord_create('lon_d_ew_rst', 'ncol_d_ew_rst', uniqpts_glob_ew, 'longitude',     &
       'degrees_east', 1, size(pelon_deg_ew), pelon_deg_ew, map=pemap_dups_ew)
  
  allocate(grid_map(3, (ie-is+2)*(je-js+1)))
  grid_map = 0
  do j = js, je
     do i = is, ie+1
        mapind = mylindex_ew(i,j)
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap_dups_ew(mapind)
     end do
  end do
  
  call cam_grid_register('FFSL_EW_RST', dyn_decomp_ew_rst, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL_EW_RST', 'cell', '', 1)
  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  nullify(lat_coord)         ! Belongs to grid
  nullify(lon_coord)         ! Belongs to grid
  

  lat_coord => horiz_coord_create('lat_d_ns', 'ncol_d_ns',  uniqpts_glob_ns, 'latitude',      &
       'degrees_north', 1, size(pelat_deg_ns), pelat_deg_ns, map=pemap_ns)
  lon_coord => horiz_coord_create('lon_d_ns', 'ncol_d_ns',  uniqpts_glob_ns, 'longitude',     &
       'degrees_east', 1, size(pelon_deg_ns), pelon_deg_ns, map=pemap_ns)
  
  allocate(grid_map(3, (ie-is+1)*(je-js+2)))
  grid_map = 0
  mapind = 1
  do j = js, je+1
     do i = is, ie
        mapind = mylindex_ns(i,j)
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap_ns(mapind)
     end do
  end do
  
  call cam_grid_register('FFSL_NS', dyn_decomp_ns, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL_NS', 'cell', '', 1)

  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  nullify(lat_coord)         ! Belongs to grid
  nullify(lon_coord)         ! Belongs to grid

! Create nw non uniq grid (includes duplicate points) for restarts

  lat_coord => horiz_coord_create('lat_d_ns_rst', 'ncol_d_ns_rst',  uniqpts_glob_ns, 'latitude',      &
       'degrees_north', 1, size(pelat_deg_ns), pelat_deg_ns, map=pemap_dups_ns)
  lon_coord => horiz_coord_create('lon_d_ns_rst', 'ncol_d_ns_rst',  uniqpts_glob_ns, 'longitude',     &
       'degrees_east', 1, size(pelon_deg_ns), pelon_deg_ns, map=pemap_dups_ns)

  allocate(grid_map(3, (ie-is+1)*(je-js+2)))
  grid_map = 0
  mapind = 1
  do j = js, je+1
     do i = is, ie
        mapind = mylindex_ns(i,j)
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap_dups_ns(mapind)
     end do
  end do
  
  call cam_grid_register('FFSL_NS_RST', dyn_decomp_ns_rst, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL_NS_RST', 'cell', '', 1)

  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  nullify(lat_coord)         ! Belongs to grid
  nullify(lon_coord)         ! Belongs to grid


  deallocate(pelon_deg)
  deallocate(pelat_deg)
  deallocate(pelon_deg_ns)
  deallocate(pelat_deg_ns)
  deallocate(pelon_deg_ew)
  deallocate(pelat_deg_ew)
  deallocate(pemap)
  deallocate(pemap_ew)
  deallocate(pemap_ns)
  deallocate(pemap_dups_ew)
  deallocate(pemap_dups_ns)
  
end subroutine define_cam_grids

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
  use cam_grid_support, only: max_hcoordname_len

  ! Dummy arguments
  character(len=max_hcoordname_len),          intent(out) :: gridname
  character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)

  gridname = 'FFSL'
  allocate(grid_attribute_names(1))
  ! For standard CAM-FV3, we need to copy the area attribute.
  ! For physgrid, the physics grid will create area
  grid_attribute_names(1) = 'cell'

end subroutine physgrid_copy_attributes_d

!=======================================================================

integer function get_dyn_grid_parm(name) result(ival)

    use pmgrid,  only: plon, plev, plat, plevp

    implicit none

    character(len=*), intent(in) :: name
    integer is,ie,js,je

    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je
    
    if(name.eq.'ne') then
        ival = -1
    else if (name.eq.'np') then
        ival = -1
    else if (name .eq. 'plat') then
        ival = plat
    else if (name .eq. 'plon') then
        ival = (je-js+1)*(ie-is+1)
    else if (name .eq. 'plev') then
        ival = plev
    else if (name .eq. 'plevp') then
        ival = plevp
    else if (name .eq. 'beglat') then
        ival = 0
    else if (name .eq. 'endlat') then
        ival = 0
    else if (name .eq. 'beglon') then
        ival = 0
    else if (name .eq. 'endlon') then
        ival = 0
    else if (name .eq. 'beglonxy') then
        ival = 1
    else if (name .eq. 'endlonxy') then
        ival = (ie - is + 1) * (je - js + 1)  ! number of column per block
    else if (name .eq. 'beglatxy') then
        ival = 1
    else if (name .eq. 'endlatxy') then
        ival = 1                              ! number of block per MPI task
    else if (name .eq. 'splon') then
        ival = -1
    else
        call endrun('get_dyn_grid_parm: undefined name: '//adjustl(trim(name)))
    end if

end function get_dyn_grid_parm

!=======================================================================

subroutine get_block_ldof_d(nlev, ldof)

    use pio, only: pio_offset

    implicit none

    integer, intent(in) :: nlev
    integer(kind=pio_offset), pointer :: ldof(:)

    call endrun('get_block_ldof_d: currently not avaliable.')

end subroutine get_block_ldof_d

!=======================================================================

function get_dyn_grid_parm_real2d(name) result(rval)

    implicit none

    character(len=*), intent(in) :: name
    real(r8), pointer :: rval(:,:)

    if (name .eq. 'clon') then
       rval => clon
    else if (name .eq. 'londeg') then
       rval => londeg
    else
       nullify(rval)
    end if

end function get_dyn_grid_parm_real2d

!=======================================================================

function get_dyn_grid_parm_real1d(name) result(rval)

    implicit none

    character(len=*), intent(in) :: name
    real(r8), pointer :: rval(:)

    if (name .eq. 'clat') then
        rval => clat
    else if (name .eq. 'latdeg') then
        rval => latdeg
    else if (name .eq. 'w') then
        rval => w
    else
        nullify(rval)
    end if

end function get_dyn_grid_parm_real1d

!#######################################################################
subroutine dyn_grid_get_colndx( igcol, ncols, owners, indx, jndx)
  use spmd_utils,       only: iam
   ! For each global column index return the owning task.  If the column is owned
   ! by this task, then also return the MPI process indicies for that column
   !
   ! NOTE: this routine needs to be updated for the physgrid

  integer, intent(in)  :: ncols
  integer, intent(in)  :: igcol(ncols)
  integer, intent(out) :: owners(ncols)
  integer, intent(out) :: indx(ncols)
  integer, intent(out) :: jndx(ncols)

  integer  :: i,is,ie,js,je
  integer  :: blockid(1), bcid(1), lclblockid(1), ind(2)

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je

  do i = 1,ncols

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))

     if ( iam==owners(i) ) then
        if (minval(abs(bcid(1)-mylindex)) .eq. 0) then
           ind = minloc(abs(bcid(1)-mylindex))
           indx(i) = is+ind(1)-1
           jndx(i) = js+ind(2)-1
        end if
     else
        indx(i) = -1
        jndx(i) = -1
     endif

  end do

end subroutine dyn_grid_get_colndx

!=======================================================================

subroutine dyn_grid_get_elem_coords(ie, rlon, rlat, cdex)

    implicit none

    integer, intent(in) :: ie ! block element index
    real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the element
    real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the element
    integer, optional, intent(out) :: cdex(:) ! global column index

    call endrun('dyn_grid_get_elem_coords: currently not avaliable.')

end subroutine dyn_grid_get_elem_coords

!=======================================================================
subroutine calc_global_index(ilocs,iloce,jlocs,jloce,agrid,nxtile,nytile,gindex,lindex,bindex,gindexdups,name)

  use fv_mp_mod,          only: mp_gather, mp_bcst
  
  implicit none
  character(10), intent(in)                                            :: name
  integer, intent(in)                                                 :: ilocs,iloce,jlocs,jloce
  real(r8), pointer, dimension(:,:,:), intent(in)                     :: agrid
  integer, intent(in)                                                 :: nxtile,nytile
  integer, dimension(ilocs:iloce,jlocs:jloce), intent(inout)          :: gindex,lindex,gindexdups
  integer, intent(inout)                                              :: bindex
 
  integer :: s, t, l, m, alen, tmp_i,n,mm
  integer :: nregions,i,j
  integer, pointer   ::   tile
  real(r8) :: tmp_r,lonrad,latrad
  integer gid,masterproc
  
  integer, allocatable :: alatlon_1d_idx(:)
  integer :: numuniq
  real(r8), allocatable :: alatlon   (:,:,:),alatlonnew   (:,:,:)
  real(r8), allocatable :: alatlon_id(:,:,:),alatlon_id_dups(:,:,:)
  logical, allocatable :: alatlon_1d_uniqmask(:),alatlon_1d_sorted_uniqmask(:)
  real(r8), allocatable :: alatlon_1d(:),alatlon_1d_sorted(:)
  real(r8), allocatable :: alatlon_1d_sorted_uniq(:)
  
! Calculate Global Index
! zlj, 2014.10.13
!-----------------------------------------------------------------------
  
  nregions = atm(mytile)%gridstruct%ntiles_g

  masterproc = mpp_root_pe()
  gid = mpp_pe()

  !-----------------------------------------------------------------------
  ! zlj, 2014.10.13
  ! Calculate Global Index

  gindex(:,:) = 0
  gindexdups(:,:) = 0
  lindex(:,:) = 0

  allocate(alatlon   ( nxtile,  nytile, nregions))
  allocate(alatlonnew   ( nxtile,  nytile, nregions))
  allocate(alatlon_1d_idx((nxtile)*(nytile)*nregions))
  allocate(alatlon_1d((nxtile)*(nytile)*nregions))
  allocate(alatlon_1d_sorted((nxtile)*(nytile)*nregions))
  allocate(alatlon_id( nxtile,  nytile, nregions))
  allocate(alatlon_id_dups( nxtile,  nytile, nregions))
  allocate(alatlon_1d_sorted_uniqmask((nxtile)*(nytile)*nregions))
  allocate(alatlon_1d_uniqmask((nxtile)*(nytile)*nregions))

  alatlon = 0._r8
  alatlonnew = 0._r8
  alatlon_1d_idx = 0._r8
  alatlon_1d = 0._r8
  alatlon_1d_sorted = 0._r8
  alatlon_id = 0._r8
  alatlon_id_dups = 0._r8
  alatlon_1d_sorted_uniqmask=.true.
  alatlon_1d_uniqmask=.true.
! there are 2 cases where multiplying lat*lon will not identify a uniq point
! 1) when one of lat/lon is 0
! 2) when 2 points have lat/lon values which are swapped ie (pt1 (lat=2,lon=10), pt2(lat=10,lon=2)) 
  do j = jlocs, jloce
     do i = ilocs, iloce
        !take care of case 1 and 2 - make identical lats and lons uniq by adding eps to one and subtracting eps to the other
        !There will be no identical 0's left after adding/subtracting eps.
        lonrad = agrid(i,j,1)+1.e-14_r8
        latrad = agrid(i,j,2)+1.e-13_r8
        alatlon(i,j,Atm(mytile)%tile) = lonrad * latrad
!        lonrad = agrid(i,j,1)
!        latrad = agrid(i,j,2)
!        if (agrid(i,j,1).eq.0._r8) lonrad=1.e-14_r8
!        if (agrid(i,j,2).eq.0._r8) latrad=1.e-14_r8
!        alatlon(i,j,Atm(mytile)%tile) = lonrad*latrad
        alatlonnew(i,j,Atm(mytile)%tile) = agrid(i,j,1)*agrid(i,j,2)
     end do
  end do

  call mp_gather(alatlon, ilocs, iloce, jlocs, jloce, nxtile, nytile, nregions)
  call mp_gather(alatlonnew, ilocs, iloce, jlocs, jloce, nxtile, nytile, nregions)

  if (gid == masterproc) then
     alen = 0
     do n = 1, nregions
        do j = 1, nytile
           do i = 1, nxtile
              alen             = alen + 1
              alatlon_1d(alen) = alatlon(i,j,n)
              alatlon_1d_idx(alen) = alen
           end do
        end do
     end do

     alatlon_1d_sorted=alatlon_1d
     do s = 1, alen - 1
        do t = s + 1, alen
           if (alatlon_1d_sorted(s) .ge. alatlon_1d_sorted(t)) then
              tmp_r         = alatlon_1d_sorted(s)
              alatlon_1d_sorted(s) = alatlon_1d_sorted(t)
              alatlon_1d_sorted(t) = tmp_r
              tmp_i         = alatlon_1d_idx(s)
              alatlon_1d_idx(s) = alatlon_1d_idx(t)
              alatlon_1d_idx(t) = tmp_i
           end if
        end do
     end do

     do s = 1, alen - 1
        if (alatlon_1d_sorted(s) .eq. alatlon_1d_sorted(s+1)) then 
           alatlon_1d_sorted_uniqmask(s)=.false.
           alatlon_1d_uniqmask(alatlon_1d_idx(s))=.false.
        end if
     end do

     numuniq=count(alatlon_1d_sorted_uniqmask)
     allocate(alatlon_1d_sorted_uniq(numuniq))
     alatlon_1d_sorted_uniq=-999.9
     alatlon_1d_sorted_uniq=pack(alatlon_1d_sorted,alatlon_1d_sorted_uniqmask)
     
     m = 0
     do n = 1, nregions
        do j = 1, nytile
           do i = 1, nxtile
              m = m + 1
              if (alatlon_1d_uniqmask(m)) then
                 do l = 1, numuniq
                    if (alatlon(i,j,n).eq.alatlon_1d_sorted_uniq(l)) then
                       alatlon_id(i,j,n) = l
                       exit
                    end if
                 end do
              else
                 alatlon_id(i,j,n) = 0
              end if
              do l = 1, numuniq
                 if (alatlon(i,j,n).eq.alatlon_1d_sorted_uniq(l)) then
                    alatlon_id_dups(i,j,n) = l
                    exit
                 end if
              end do

           end do
        end do
     end do

     deallocate(alatlon_1d_sorted_uniq)

  end if

  call mp_bcst(alatlon_id, nxtile, nytile, nregions)
  call mp_bcst(alatlon_id_dups, nxtile, nytile, nregions)

  n = 1
  do j = jlocs, jloce
     do i = ilocs, iloce
        gindex(i,j) = nint(alatlon_id(i,j,Atm(mytile)%tile))
        lindex(i,j) = n
        gindexdups(i,j) = nint(alatlon_id_dups(i,j,Atm(mytile)%tile))
        n=n+1
     end do
  end do

!
! Have taken care of intertile duplicates for ew and ns grids. Need to 
! account for duplicates of differing tasks within a tile. ie and je 
! for every task can be 0'd out since they are repeated versions of 
! is,js on a neighboring task.  ie and je already accounted for on 
! tile boundaries after executing above code, don't want to touch those. 
!
! For ew grids 0 out ie task column (east most column of task) if not already 
! eastern tile boundary
  if (nxtile.gt.nytile .and. iloce .ne. nytile) then
     do j = jlocs, jloce
        gindex(iloce,j) = 0
     end do
  end if

! For ns grids 0 out je row if not already on north tile edge
  if (nytile.gt.nxtile .and. jloce .ne. nytile) then
     do i = ilocs, iloce
        gindex(i,jloce) = 0
     end do
  end if

  deallocate(alatlon)
  deallocate(alatlonnew)
  deallocate(alatlon_1d_idx)
  deallocate(alatlon_1d)
  deallocate(alatlon_1d_sorted)
  deallocate(alatlon_id)
  deallocate(alatlon_id_dups)
  deallocate(alatlon_1d_sorted_uniqmask)
  deallocate(alatlon_1d_uniqmask)
  bindex = mpp_pe() + 1

  ! Calculate Global Index
  ! zlj, 2014.10.13
  !-----------------------------------------------------------------------
end subroutine calc_global_index

subroutine calc_global_indexjt(ilocs,iloce,jlocs,jloce,nxtile,nytile,tile,nregions,locgindex,lindex,bindex,locgindex_justdups,locgindex_tile, name,uniq_pts_loc,uniq_pts_glob)

  use fv_mp_mod,          only: mp_gather, mp_bcst
  
  implicit none
  character(10), intent(in)                                            :: name
  integer, intent(in)                                                 :: ilocs,iloce,jlocs,jloce
  integer, intent(in)                                                 :: tile,nregions
  integer, intent(in)                                                 :: nxtile,nytile
  integer, intent(inout)                                              :: uniq_pts_loc,uniq_pts_glob
  integer, dimension(ilocs:iloce,jlocs:jloce), intent(inout)          :: locgindex,lindex,locgindex_tile,locgindex_justdups
  integer, intent(inout)                                              :: bindex
 
  integer :: i,j,k,n,alen,tile_alen
  integer gid,masterproc
  real(r8), allocatable :: gindex_zeros4dups(:,:,:)
  real(r8), allocatable :: gindex_tiles(:,:,:)
  integer, dimension(ilocs:iloce,jlocs:jloce)                         :: locgindex_wdups

! Calculate Global Index
!-----------------------------------------------------------------------
  
  masterproc = mpp_root_pe()
  gid = mpp_pe()+1

  !-----------------------------------------------------------------------
  ! zlj, 2014.10.13
  ! Calculate Global Index

  locgindex(:,:) = 0
  locgindex_tile(:,:) = 0
  locgindex_wdups(:,:) = 0
  lindex(:,:) = 0          
  allocate(gindex_zeros4dups    ( nxtile,  nytile, nregions))
  allocate(gindex_tiles    ( nxtile,  nytile, nregions))
  gindex_zeros4dups=0
  gindex_tiles=0

  bindex = mpp_pe() + 1
  
  !                             5 6
  !Tiles are laid out:        3 4
  !                         1 2
  ! 
  ! pattern repeats so tile 1 also sits on top of 6

  ! build the indicies.  All grid types will use these counts for the
  ! local version of gindex. Each grid type will just zero out the appropropriate
  ! rows or columns depending on the grid type tile number.


  if (nytile.gt.nxtile) then ! create  locgindex for NS grid

     ! for ns grids 
     ! if tile 2,4,6 then don't count north boundary points
     ! lindex ordering by rows

     gindex_zeros4dups=0
     alen=1
     tile_alen=1
     do n = 1, nregions
        do j = 1, nytile
           do i = 1, nxtile
              gindex_tiles(i,j,n)=tile_alen
              tile_alen=tile_alen+1
              if (n.eq.2.or.n.eq.4.or.n.eq.6) then
                 if (j.ne.nytile) then
                    gindex_zeros4dups(i,j,n)=alen
                    alen=alen+1
                 else
                    gindex_zeros4dups(i,j,n)=0
                 end if
              else
                 gindex_zeros4dups(i,j,n)=alen
                 alen=alen+1
              end if
           end do
        end do
     end do

     locgindex(ilocs:iloce,jlocs:jloce)=gindex_zeros4dups(ilocs:iloce,jlocs:jloce,tile)

     ! appropriate tiles boundarys (nytile) already 0'd from gindex_zeros4dups need to 
     ! zero inner tile je boundarys (These are also repeated points between tasks in ns direction))
     
     if (jloce.ne.nytile) then
        locgindex(ilocs:iloce,jloce)=0
     end if

     ! output local and global uniq points
     uniq_pts_glob=count(gindex_zeros4dups.ne.0)
     uniq_pts_loc=count(locgindex(ilocs:iloce,jlocs:jloce).ne.0)

     ! locgindex_wdups eq gindex_zeros4dups with all northern boundaries 0's filled from adjacent row on next tile
     locgindex_wdups(ilocs:iloce,jlocs:jloce)=gindex_zeros4dups(ilocs:iloce,jlocs:jloce,tile)
     if (tile.eq.2.or.tile.eq.4) then !need to pull row 1 from next tile.
        if (jloce.eq.nytile) locgindex_wdups(ilocs:iloce,jloce)=gindex_zeros4dups(ilocs:iloce,1,tile+1)
     end if
     if (tile.eq.6) then ! need to pull 1st row of tile 1
        if (jloce.eq.nytile) locgindex_wdups(ilocs:iloce,jloce)=gindex_zeros4dups(ilocs:iloce,1,1)
     endif

     locgindex_justdups=0
     where (locgindex(ilocs:iloce,jlocs:jloce).eq.0)
        locgindex_justdups=locgindex_wdups
     end where

     !lindex ordered by row
     alen=1
     do j = jlocs, jloce
        do i = ilocs, iloce
           lindex(i,j)=alen
           locgindex_tile(i,j)=gindex_tiles(i,j,tile)
           alen=alen+1
        end do
     end do

  end if

  if (nxtile.gt.nytile) then ! create locgindex for EW grid

     ! for ew grids 
     ! if tile 1,3,5 then don't count East boundary points
     ! lindex ordering by row

     gindex_zeros4dups=0
     alen=1
     tile_alen=1
     do n = 1, nregions
        do j = 1, nytile
           do i = 1, nxtile
              gindex_tiles(i,j,n)=tile_alen
              tile_alen=tile_alen+1
              if (n.eq.1.or.n.eq.3.or.n.eq.5) then
                 if (i.ne.nxtile) then
                    gindex_zeros4dups(i,j,n)=alen
                    alen=alen+1
                 else
                    gindex_zeros4dups(i,j,n)=0
                 end if
              else
                 gindex_zeros4dups(i,j,n)=alen
                 alen=alen+1
              end if
           end do
        end do
     end do

     locgindex(ilocs:iloce,jlocs:jloce)=gindex_zeros4dups(ilocs:iloce,jlocs:jloce,tile)

     ! appropriate tiles boundarys (nxtile) already 0'd from gindex_zeros4dups need to 
     ! zero inner tile ie boundarys (These are also repeated points between tasks in ew direction)
     if (iloce.ne.nxtile) then
        locgindex(iloce,jlocs:jloce)=0
     end if

     ! output local and global uniq points
     uniq_pts_glob=count(gindex_zeros4dups.ne.0)
     uniq_pts_loc=count(locgindex(ilocs:iloce,jlocs:jloce).ne.0)

     ! locgindex_wdups eq gindex_zeros4dups with eastern boundaries of tiles 1,3,5 (zeros) 
     ! filled from adjacent column on next tile
     locgindex_wdups(ilocs:iloce,jlocs:jloce)=gindex_zeros4dups(ilocs:iloce,jlocs:jloce,tile)
     if (tile.eq.1.or.tile.eq.3.or.tile.eq.5) then !need to pull row 1 from next tile.
        if (iloce.eq.nxtile) locgindex_wdups(iloce,jlocs:jloce)=gindex_zeros4dups(1,jlocs:jloce,tile+1)
     end if
     
     locgindex_justdups=0
     where (locgindex(ilocs:iloce,jlocs:jloce).eq.0)
        locgindex_justdups=locgindex_wdups
     end where

     alen=1
     do j = jlocs, jloce
        do i = ilocs, iloce
           lindex(i,j)=alen
           locgindex_tile(i,j)=gindex_tiles(i,j,tile)
           alen=alen+1
        end do
     end do

  end if
  
  if (nytile.eq.nxtile) then ! create  locgindex for A

     ! for a-grids just sequential ordering
     ! lindex ordering by rows
     gindex_zeros4dups=0
     alen=1
     do n = 1, nregions
        do j = 1, nytile
           do i = 1, nxtile
              gindex_tiles(i,j,n)=alen
              gindex_zeros4dups(i,j,n)=alen
              alen=alen+1
           end do
        end do
     end do

     locgindex_wdups(ilocs:iloce,jlocs:jloce)=gindex_zeros4dups(ilocs:iloce,jlocs:jloce,tile)
     locgindex(ilocs:iloce,jlocs:jloce)=gindex_zeros4dups(ilocs:iloce,jlocs:jloce,tile)
     locgindex_justdups=0
     where (locgindex(ilocs:iloce,jlocs:jloce).eq.0)
        locgindex_justdups=locgindex_wdups
     end where

     ! output local and global uniq points
     uniq_pts_glob=count(gindex_zeros4dups.ne.0)
     uniq_pts_loc=count(locgindex(ilocs:iloce,jlocs:jloce).ne.0)
     !lindex ordered by row
     alen=1
     do j = jlocs, jloce
        do i = ilocs, iloce
           lindex(i,j)=alen
           locgindex_tile(i,j)=gindex_tiles(i,j,tile)
           alen=alen+1
        end do
     end do

  end if

  deallocate(gindex_zeros4dups)
  deallocate(gindex_tiles)

  !-----------------------------------------------------------------------
end subroutine calc_global_indexjt

end module dyn_grid
