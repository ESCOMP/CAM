module spmd_dyn

!-----------------------------------------------------------------------
! Subroutines to initialize SPMD implementation of FV
!
! !REVISION HISTORY:
!   00.09.30  Sawyer             Alterations for LR SPMD mode
!   01.05.09  Mirin              2-D yz decomposition
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.12.20  Sawyer             Changed index order of Q3 decomposition
!   03.05.07  Sawyer             Removed unneeded decompositions
!   06.03.01  Sawyer             Removed tracertrans-related variables
!-----------------------------------------------------------------------

use spmd_utils,         only: iam, masterproc, npes, mpicom

use pmgrid,             only: plat, plon, plev, spmd_on
use constituents,       only: pcnst

use dynamics_vars,      only: t_fvdycore_grid
use dyn_internal_state, only: get_dyn_state_grid

use cam_abortutils,     only: endrun
use cam_logfile,        only: iulog

implicit none
private
save

public :: &
   spmd_readnl,       &
   spmdinit_dyn,      &
   compute_gsfactors, &
   spmdbuf

public :: &
   local_dp_map, &
   block_buf_nrecs, &
   chunk_buf_nrecs, &
   proc,            &
   lonrangexy,      &
   latrangexy

! local_dp_map, block_buf_nrecs, chunk_buf_nrecs belong somewhere else.  They are just
! stored here without being set or used
logical :: local_dp_map=.false.    ! flag indicates that mapping between dynamics 
                                   !  and physics decompositions does not require 
                                   !  interprocess communication
integer :: block_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                   !  in dynamics decomposition (including level 0)
integer :: chunk_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                   !  in physics decomposition (including level 0)

! used by dyn_grid::get_block_owner_d
integer :: proc(plat)              ! processor id associated with a given lat.

! used by dyn_grid:: get_gcol_block_d, get_block_gcol_cnt_d, get_block_gcol_d
integer, allocatable :: lonrangexy(:,:)   ! global xy-longitude subdomain index
integer, allocatable :: latrangexy(:,:)   ! global xy-latitude subdomain index

integer :: force_2d = 0                 !option to force transpose computation for 1D decomp.
integer :: geopkblocks = 1              !number of stages to use in Z-serial non-transpose
                                               ! geopotential method (routine geopk_d)
                                               ! with 2D decomp.
logical :: geopkdist  = .false.         !use a distributed method for geopotential calculation 
                                               ! with 2D decomp.
logical :: geopk16byte   = .false.      !use Z-parallel distributed method for geopotential 
                                               ! calculation with 2D decomp.; otherwise, use Z-serial 
                                               ! pipeline algorithm
integer :: geopktrans = 0               
integer :: npr_yz(4)                    !yz and xy decompositions
integer :: modcomm_transpose = 0        !mod_comm transpose method
                                               !   0 for temporary contiguous buffers
                                               !   1 for mpi derived types
integer :: modcomm_geopk = 0            !mod_comm geopk method
                                               !   0 for temporary contiguous buffers
                                               !   1 for mpi derived types
integer :: modcomm_gatscat = 0          !mod_comm gather/scatter method
                                               !   0 for temporary contiguous buffers
                                               !   1 for mpi derived types
integer :: modc_sw_dynrun = 0           !mod_comm irregular underlying communication method for dyn_run/misc
                                               !  0 for original mp_sendirr/mp_recvirr
                                               !  1 for mp_swapirr and point-to-point communications
                                               !  2 for mp_swapirr and all-to-all communications
logical :: modc_hs_dynrun = .true.      !mod_comm irreg comm handshaking for dyn_run/misc
logical :: modc_send_dynrun = .true.    ! true for mod_comm irregular communication blocking send for
                                               ! dyn_run/misc, false for nonblocking send
integer :: modc_mxreq_dynrun = -1       !maximum number of nonblocking communication requests to allow
                                               ! when using mp_swapirr and point-to-point communications for
                                               ! dyn_run/misc
                                               ! < 0 implies no limits
integer :: modc_sw_cdcore = 0           !mod_comm irregular underlying communication method for cd_core/geopk
                                               !  0 for original mp_sendirr/mp_recvirr
                                               !  1 for mp_swapirr and point-to-point communications
                                               !  2 for mp_swapirr and all-to-all communications
logical :: modc_hs_cdcore = .true.      ! true for mod_comm irregular communication handshaking for cd_core/geopk
logical :: modc_send_cdcore  = .true.   ! true for geopk_d or mod_comm irregular communication blocking send for
                                               !  cd_core/geopk, false for nonblocking send
integer :: modc_mxreq_cdcore = -1       ! maximum number of nonblocking communication requests to allow
                                               !  when using mp_swapirr and point-to-point communications for
                                               !  cd_core/geopk
                                               !  < 0 implies no limits
integer :: modc_sw_gather = 1           ! mod_comm irregular underlying communication method for gather
                                               !  0 for original mp_sendirr/mp_recvirr
                                               !  1 for mp_swapirr and point-to-point communications
                                               !  2 for mp_swapirr and all-to-all communications
logical :: modc_hs_gather = .true.      ! true for mod_comm irregular communication handshaking for gather
logical :: modc_send_gather = .true.    ! true for mod_comm irregular communication blocking send for
                                               !  gather, false for nonblocking send
integer :: modc_mxreq_gather = 64       ! maximum number of nonblocking communication requests to allow
                                               !  when using mp_swapirr and point-to-point communications for
                                               !  gather
                                               !  < 0 implies no limits
integer :: modc_sw_scatter = 0          ! mod_comm irregular underlying communication method for scatter
                                               !  0 for original mp_sendirr/mp_recvirr
                                               !  1 for mp_swapirr and point-to-point communications
                                               !  2 for mp_swapirr and all-to-all communications
logical :: modc_hs_scatter = .false.    ! true for mod_comm irregular communication handshaking for scatter
logical :: modc_send_scatter = .true.   ! true for mod_comm irregular communication blocking send for
                                               !  scatter, false for nonblocking send
integer :: modc_mxreq_scatter = -1      ! maximum number of nonblocking communication requests to allow
                                               !  when using mp_swapirr and point-to-point communications for
                                               !  scatter
                                               !  < 0 implies no limits
integer :: modc_sw_tracer = 0           ! mod_comm irregular underlying communication method for multiple tracers
                                               !  0 for original mp_sendirr/mp_recvirr
                                               !  1 for mp_swapirr and point-to-point communications
                                               !  2 for mp_swapirr and all-to-all communications
logical :: modc_hs_tracer = .true.      ! true for mod_comm irregular communication handshaking for multiple tracers
logical :: modc_send_tracer = .true.    ! true for mod_comm irregular communication blocking send for
                                               !  multiple tracers, false for nonblocking send
integer :: modc_mxreq_tracer = -1       ! maximum number of nonblocking communication requests to allow
                                               !  when using mp_swapirr and point-to-point communications for
                                               !  multiple tracers
                                               !  < 0 implies no limits
integer :: modc_onetwo = 2              !one or two simultaneous mod_comm irregular communications 
                                               !   (excl. tracers)
integer :: modc_tracers = 3             ! max number of tracers for simultaneous mod_comm irregular communications 
                                               !  0 for original mp_sendirr/mp_recvirr communications
                                               !  positive for special tracer routines

integer :: fv_ct_overlap  = 0           ! nonzero for overlap of cd_core and trac2d, 0 otherwise
integer :: fv_trac_decomp = 1           ! size of tracer domain decomposition for trac2d

! these vars used in beglat and endlat determination, and in compute_gsfactors
integer, allocatable :: nlat_p(:)  ! number of latitudes per subdomain in YZ decomp
integer, allocatable :: cut(:,:)   ! latitude partition for MPI tasks in YZ decomp

integer :: npr_y
integer :: npr_z
integer :: nprxy_x
integer :: nprxy_y
integer :: npes_yz
integer :: npes_xy
integer :: myid_y
integer :: myid_z
integer :: myidxy_x
integer :: myidxy_y

integer :: numlats

!========================================================================================
contains
!========================================================================================

subroutine spmd_readnl(nlfilename)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mstrid=>masterprocid, mpi_integer, mpi_logical,&
                              mpi_success
   ! args
   character(len=*), intent(in) :: nlfilename

   ! Local variables
   integer :: ierr           ! error code
   integer :: unitn          ! namelist unit number

   namelist /spmd_fv_inparm/ npr_yz,             &
        geopktrans, geopkblocks,                 &
        force_2d, modcomm_transpose,             &
        modcomm_geopk, modcomm_gatscat,          &
        modc_sw_dynrun, modc_hs_dynrun,          &
        modc_send_dynrun, modc_mxreq_dynrun,     &
        modc_sw_cdcore, modc_hs_cdcore,          &
        modc_send_cdcore, modc_mxreq_cdcore,     &
        modc_sw_gather, modc_hs_gather,          &
        modc_send_gather, modc_mxreq_gather,     &
        modc_sw_scatter, modc_hs_scatter,        &
        modc_send_scatter, modc_mxreq_scatter,   &
        modc_sw_tracer, modc_hs_tracer,          &
        modc_send_tracer, modc_mxreq_tracer,     &
        modc_onetwo, modc_tracers,               &
        fv_ct_overlap, fv_trac_decomp
   
   character(len=*), parameter ::  sub = "spmd_readnl"

   type(t_fvdycore_grid), pointer :: grid

   integer :: color, ierror, ntemp
   integer :: twod_decomp
   integer :: mpicom_yz         ! communicator for yz decomposition
   integer :: mpicom_nyz        ! communicator for multiple yz decomposition
   integer :: mpicom_xy         ! communicator for xy decomposition
   !----------------------------------------------------------------------

   ! Default 1D domain decomposition
   npr_yz(1) = npes
   npr_yz(2) = 1
   npr_yz(3) = 1
   npr_yz(4) = npes

   if (masterproc) then
      write(iulog,*) sub//': Read in spmd_fv_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      ! Look for spmd_fv_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'spmd_fv_inparm', status=ierr)
      if (ierr == 0) then  ! found spmd_fv_inparm
         read(unitn, spmd_fv_inparm, iostat=ierr)  ! read the spmd_fv_inparm namelist group
         if (ierr /= 0) then
            call endrun(sub//': ERROR reading namelist spmd_fv_inparm')
         end if
      end if
      close( unitn )
      call freeunit( unitn )
   endif

   call mpi_bcast(npr_yz, 4, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: npr_yz")

   call mpi_bcast(geopktrans, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: geopktrans")

   call mpi_bcast(geopkblocks, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: geopkblocks")

   call mpi_bcast(force_2d, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: force_2d")

   call mpi_bcast(modcomm_transpose, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modcomm_transpose")

   call mpi_bcast(modcomm_geopk, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modcomm_geopk")

   call mpi_bcast(modcomm_gatscat, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modcomm_gatscat")

   call mpi_bcast(modc_sw_dynrun, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_sw_dynrun")

   call mpi_bcast(modc_hs_dynrun, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_hs_dynrun")

   call mpi_bcast(modc_send_dynrun, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_send_dynrun")

   call mpi_bcast(modc_mxreq_dynrun, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_mxreq_dynrun")

   call mpi_bcast(modc_sw_cdcore, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_sw_cdcore")

   call mpi_bcast(modc_hs_cdcore, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_hs_cdcore")

   call mpi_bcast(modc_send_cdcore, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_send_cdcore")

   call mpi_bcast(modc_mxreq_cdcore, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_mxreq_cdcore")

   call mpi_bcast(modc_sw_gather, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_sw_gather")

   call mpi_bcast(modc_hs_gather, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_hs_gather")

   call mpi_bcast(modc_send_gather, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_send_gather")

   call mpi_bcast(modc_mxreq_gather, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_mxreq_gather")

   call mpi_bcast(modc_sw_scatter, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_sw_scatter")

   call mpi_bcast(modc_hs_scatter, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_hs_scatter")

   call mpi_bcast(modc_send_scatter, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_send_scatter")

   call mpi_bcast(modc_mxreq_scatter, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_mxreq_scatter")

   call mpi_bcast(modc_sw_tracer, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_sw_tracer")

   call mpi_bcast(modc_hs_tracer, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_hs_tracer")

   call mpi_bcast(modc_send_tracer, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_send_tracer")

   call mpi_bcast(modc_mxreq_tracer, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_mxreq_tracer")

   call mpi_bcast(modc_onetwo, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_onetwo")

   call mpi_bcast(modc_tracers, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: modc_tracers")

   call mpi_bcast(fv_ct_overlap, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_ct_overlap")

   call mpi_bcast(fv_trac_decomp, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_trac_decomp")

   ! Put namelist input into the grid object
   grid => get_dyn_state_grid()

   npr_y   = npr_yz(1)
   npr_z   = npr_yz(2)
   nprxy_x = npr_yz(3)
   nprxy_y = npr_yz(4)
   npes_yz = npr_y*npr_z
   npes_xy = nprxy_x*nprxy_y
   if (npes_yz < 1) then
      call endrun(sub//': ERROR: yz domain decomposition must have at least 1 subdomain')
   endif
   if (npes_yz > npes) then
      call endrun(sub//': ERROR: incorrect yz domain decomposition')
   endif
   if (npes_xy > npes) then
      call endrun(sub//': ERROR: incorrect xy domain decomposition')
   endif

   grid%npr_y    = npr_y
   grid%npr_z    = npr_z
   grid%nprxy_x  = nprxy_x
   grid%nprxy_y  = nprxy_y
   grid%npes_xy  = npes_xy
   grid%npes_yz  = npes_yz

   grid%ct_overlap  = fv_ct_overlap
   grid%trac_decomp = fv_trac_decomp

   if (fv_ct_overlap .ne. 0 .and. npes .lt. 2*npes_yz) then
      call endrun(sub//': ERROR: Not enough processes to overlap cd_core and trac2d')
   end if

   if (fv_trac_decomp .le. 0) then
      call endrun(sub//': ERROR: fv_trac_decomp improperly initialized')
   end if

   if (npes .lt. fv_trac_decomp*npes_yz) then
      call endrun(sub//': ERROR: Not enough processes to decompose tracers')
   endif

   if (fv_ct_overlap .gt. 0 .and. fv_trac_decomp .gt. 1) then
      call endrun(sub//': ERROR: Cannot simultaneously overlap cd_core/trac2d and decompose tracers')
   endif

   ! Tracer decomposition limits
   allocate(grid%ktloa(fv_trac_decomp), grid%kthia(fv_trac_decomp))
   grid%ktloa(:) = 1
   grid%kthia(:) = pcnst
   grid%ktlo = 1
   grid%kthi = pcnst

   grid%commdyn  = mpicom
   grid%iam      = iam

#ifdef SPMD
   myid_z   = iam/npr_y
   myid_y   = iam - myid_z*npr_y
   color = iam/npes_yz
   call mpi_comm_split(mpicom, color, iam, mpicom_yz, ierror)
   if (ierror /= mpi_success) then
      write(iulog,*) sub//': ERROR: mpi_comm_split_yz failed with IER=', ierror
      call endrun(sub//': ERROR: mpi_comm_split_yz failed')
   end if
   call mpi_comm_size(mpicom_yz, ntemp, ierror)
   if (iam .lt. npes_yz .and. ntemp .ne. npes_yz) then
      write(iulog,*) sub//': ERROR: mpicom_yz has incorrect size of ', ntemp
      call endrun(sub//': ERROR: mpicom_yz has incorrect size')
   end if

   if (fv_ct_overlap .gt. 0 .or. fv_trac_decomp .gt. 1) then
      ! These are mutually exclusive options
      if ((fv_ct_overlap .gt. 0 .and. iam .lt. 2*npes_yz) .or.         &
           (fv_trac_decomp .gt. 1 .and. iam .lt. fv_trac_decomp*npes_yz)) then
         color = 1
      else
         color = 0
      endif
      call mpi_comm_split(mpicom, color, iam, mpicom_nyz, ierror)
      if (ierror /= mpi_success) then
         write (iulog,*) sub//': ERROR: mpi_comm_split_nyz failed with IER=', ierror
         call endrun(sub//': ERROR: mpi_comm_split_nyz failed')
      endif
   else
      mpicom_nyz = mpicom_yz
   endif

   myidxy_y = iam/nprxy_x
   myidxy_x = iam - myidxy_y*nprxy_x
   color = iam/npes_xy
   call mpi_comm_split(mpicom, color, iam, mpicom_xy, ierror)
   if (ierror /= mpi_success) then
      write(iulog,*) sub//': ERROR: mpi_comm_split_xy failed with IER=', ierror
      call endrun(sub//': ERROR: mpi_comm_split_xy failed')
   endif
   call mpi_comm_size(mpicom_xy, ntemp, ierror)
   if (iam .lt. npes_xy .and. ntemp .ne. npes_xy) then
      write(iulog,*) sub//': ERROR: mpicom_xy has incorrect size of ', ntemp
      call endrun(sub//': ERROR: mpicom_xy has incorrect size')
   endif

   grid%myid_z = myid_z
   grid%myid_y = myid_y
   grid%myidxy_y = myidxy_y
   grid%myidxy_x = myidxy_x

   grid%commyz   = mpicom_yz
   grid%commnyz  = mpicom_nyz
   grid%commxy   = mpicom_xy
#endif

   twod_decomp = 0
   if (npr_z > 1 .or. nprxy_x > 1 .or. force_2d .eq. 1) then
      twod_decomp = 1
   endif
   grid%twod_decomp = twod_decomp

   if (geopktrans .ne. 0) geopkdist   = .true.
   if (geopktrans .eq. 1) geopk16byte = .true.
#ifdef NO_CRAY_POINTERS
   if (geopk16byte) then
      call endrun (sub//': ERROR: cannot use geopk16 unless compiler supports cray pointers')
   end if
#endif

   geopkblocks = max(1,geopkblocks)

   grid%geopkdist   = geopkdist
   grid%geopk16byte = geopk16byte
   grid%geopkblocks = geopkblocks
   grid%mod_method  = modcomm_transpose
   grid%mod_geopk   = modcomm_geopk
   grid%mod_gatscat = modcomm_gatscat

   if (modc_sw_dynrun > 0 .and. modcomm_transpose > 0) then
      call endrun(sub//': ERROR: modc_sw_dynrun and modcomm_transpose are inconsistent')
   endif

   grid%modc_dynrun(1)   = modc_sw_dynrun
   grid%modc_dynrun(2:3) = 0
   if (modc_hs_dynrun)   grid%modc_dynrun(2) = 1
   if (modc_send_dynrun) grid%modc_dynrun(3) = 1
   grid%modc_dynrun(4)   = modc_mxreq_dynrun

   if (modc_sw_cdcore > 0 .and. (modcomm_transpose > 0  .or. &
                                (modcomm_geopk > 0 .and. geopk16byte))) then
      call endrun(sub//': ERROR: modc_sw_cdcore and modcomm_transpose are inconsistent')
   endif

   grid%modc_cdcore(1)   = modc_sw_cdcore
   grid%modc_cdcore(2:3) = 0
   if (modc_hs_cdcore)   grid%modc_cdcore(2) = 1
   if (modc_send_cdcore) grid%modc_cdcore(3) = 1
   grid%modc_cdcore(4)   = modc_mxreq_cdcore

   if (modc_sw_gather > 0 .and. modcomm_gatscat > 0) then
      call endrun(sub//': ERROR: modc_sw_gather and modcomm_gatscat are inconsistent')
   endif

   grid%modc_gather(1)   = modc_sw_gather
   grid%modc_gather(2:3) = 0
   if (modc_hs_gather)   grid%modc_gather(2) = 1
   if (modc_send_gather) grid%modc_gather(3) = 1
   grid%modc_gather(4)   = modc_mxreq_gather

   if (modc_sw_scatter > 0 .and. modcomm_gatscat > 0) then
      call endrun(sub//': ERROR: modc_sw_scatter and modcomm_gatscat are inconsistent')
   endif

   grid%modc_scatter(1)   = modc_sw_scatter
   grid%modc_scatter(2:3) = 0
   if (modc_hs_scatter)   grid%modc_scatter(2) = 1
   if (modc_send_scatter) grid%modc_scatter(3) = 1
   grid%modc_scatter(4)    = modc_mxreq_scatter

   if (modc_sw_tracer > 0 .and. modcomm_transpose > 0) then
      call endrun(sub//': ERROR: modc_sw_tracer and modcomm_transpose are inconsistent')
   endif

   grid%modc_tracer(1)   = modc_sw_tracer
   grid%modc_tracer(2:3) = 0
   if (modc_hs_tracer)   grid%modc_tracer(2) = 1
   if (modc_send_tracer) grid%modc_tracer(3) = 1
   grid%modc_tracer(4)   = modc_mxreq_tracer

   if (modc_tracers < 0) then
      call endrun(sub//': ERROR: inadmissable value of modc_tracers')
   endif

   grid%modc_onetwo       = modc_onetwo
   grid%modc_tracers      = modc_tracers

   if (masterproc) then
      write(iulog,*)'FV grid decomposition:'
      write(iulog,*)'  npr_y = ', npr_y, '  npr_z = ', npr_z
      write(iulog,*)'  nprxy_x = ', nprxy_x, '  nprxy_y = ', nprxy_y
      write(iulog,*)'  npes = ', npes, '  npes_yz= ', npes_yz, '  npes_xy = ', npes_xy

      if (npes > 1) then

         if (fv_ct_overlap .ne. 0) then
            write(iulog,*)'  Overlapping tracer and dynamics subcycles'
         endif

         write(iulog,*)'  Decomposing tracers into ', fv_trac_decomp, ' groups'

         if (twod_decomp == 0) then
            write(iulog,*)'  decomposition is effectively 1D - skipping transposes'
         else
            write(iulog,*)'  using multi-2d decomposition methodology'
         end if

         write(iulog,*)'  non-transpose geopk communication method = ', geopkdist
         write(iulog,*)'  Z-parallel non-transpose geopk communication method = ', geopk16byte

         if (geopkdist .and. (.not. geopk16byte)) then
            write(iulog,*)'  number of stages in Z-serial non-transpose geopk method = ', geopkblocks
         endif

         write(iulog,*)'  modcomm transpose method = ', modcomm_transpose
         write(iulog,*)'  modcomm geopk method = ', modcomm_geopk
         write(iulog,*)'  modcomm gatscat method = ', modcomm_gatscat

         write(iulog,*)'  modc_sw_dynrun = ', modc_sw_dynrun
         write(iulog,*)'  modc_hs_dynrun = ', modc_hs_dynrun
         write(iulog,*)'  modc_send_dynrun = ', modc_send_dynrun
         write(iulog,*)'  modc_mxreq_dynrun = ', modc_mxreq_dynrun

         write(iulog,*)'  modc_sw_cdcore = ', modc_sw_cdcore
         write(iulog,*)'  modc_hs_cdcore = ', modc_hs_cdcore
         write(iulog,*)'  modc_send_cdcore = ', modc_send_cdcore
         write(iulog,*)'  modc_mxreq_cdcore = ', modc_mxreq_cdcore

         write(iulog,*)'  modc_sw_gather = ', modc_sw_gather
         write(iulog,*)'  modc_hs_gather = ', modc_hs_gather
         write(iulog,*)'  modc_send_gather = ', modc_send_gather
         write(iulog,*)'  modc_mxreq_gather = ', modc_mxreq_gather

         write(iulog,*)'  modc_sw_scatter = ', modc_sw_scatter
         write(iulog,*)'  modc_hs_scatter = ', modc_hs_scatter
         write(iulog,*)'  modc_send_scatter = ', modc_send_scatter
         write(iulog,*)'  modc_mxreq_scatter = ', modc_mxreq_scatter

         write(iulog,*)'  modc_sw_tracer = ', modc_sw_tracer
         write(iulog,*)'  modc_hs_tracer = ', modc_hs_tracer
         write(iulog,*)'  modc_send_tracer = ', modc_send_tracer
         write(iulog,*)'  modc_mxreq_tracer = ', modc_mxreq_tracer

         write(iulog,*)'  modc_onetwo = ', modc_onetwo
         write(iulog,*)'  modc_tracers = ', modc_tracers

      end if

   end if

end subroutine spmd_readnl

!========================================================================================

subroutine spmdinit_dyn(jord, grid)

   !-----------------------------------------------------------------------
   ! Initialize grid decomposition.
   !
   ! !REVISION HISTORY:
   !   00.09.30  Sawyer             Added LR-specific initialization
   !   01.06.27  Mirin              Secondary 2-D xy decomposition
   !   01.10.16  Sawyer             Added Y at each Z decompositions
   !   03.07.22  Sawyer             Removed decomps used by highp2
   !-----------------------------------------------------------------------

#ifdef SPMD
   use parutilitiesmodule, only: parinit, gid, parcollective, maxop
   use dynamics_vars,      only: spmd_vars_init
#endif

   ! arguments
   integer,            intent(in) :: jord
   type(t_fvdycore_grid), pointer :: grid

   ! local variables:

   integer, parameter :: numbnd = 0        ! no.of latitudes passed N and S of forecast lat

   integer :: beglat
   integer :: endlat
   integer :: beglev
   integer :: endlev

   integer :: mod_maxirr

   integer :: procid    ! processor id
   integer :: procids   ! processor id
   integer :: procidn   ! processor id

   integer :: j, k
   integer :: lat
   integer :: vert
   integer :: lonn
   integer :: workleft  ! amount of work still to be parcelled out

   integer :: isum      ! running total of work parcelled out
   integer :: smostlat  ! southern-most latitude index
   integer :: nmostlat  ! northern-most latitude index

   integer, allocatable :: ydist(:) ! number of lats per subdomain in YZ decomp
   integer, allocatable :: zdist(:) ! number of levels per subdomain in YZ decomp

   integer :: beglonxy, endlonxy
   integer :: beglatxy, endlatxy
   integer, allocatable :: xdistxy(:) ! number of xy-longs per subdomain
   integer, allocatable :: ydistxy(:) ! number of xy-lats per subdomain

   integer, allocatable :: tmp(:)
   integer, allocatable :: jmyz(:), kmyz(:) ! used for nonblocking receive
   integer, allocatable :: imxy(:), jmxy(:) ! used for nonblocking receive

   character(len=*), parameter ::  sub = "spmdinit_dyn"
   !---------------------------------------------------------------------------

   ! Default YZ decomposition is 1D
   beglev = 1
   endlev = plev

   mod_maxirr = max(modc_onetwo, modc_tracers)

#ifdef SPMD

   spmd_on = 1

   ! Initialize PILGRIM library
   call parinit(comm=mpicom, &
                npryzxy = (/ npr_y, npr_z, nprxy_x, nprxy_y /), &
                mod_method  = modcomm_transpose, &
                mod_geopk   = modcomm_geopk,     &
                mod_maxirr  = mod_maxirr,    &
                mod_gatscat = modcomm_gatscat )

   ! Compute Y partition for YZ decomposition

   allocate(ydist  (npr_y))
   allocate(nlat_p (0:npes-1))
   allocate(cut    (2,0:npes-1))

   ydist(:) = 0
   nlat_p(:) = 0
   cut(1,:) = -1
   cut(2,:) = -2

   lat = plat / npr_y
   workleft = plat - lat * npr_y
   if (lat < 3) then
      call endrun(sub//': ERROR: less than 3 latitudes per subdomain')
   end if

   ! Be careful:  ydist is 1-based.  CAMs arrays, e.g., cut,  are 0-based

   do procid = 1, npr_y
      ydist(procid) = lat
   end do

   if ( workleft /= 0 ) then
      procids = (npr_y+1) / 2
      procidn = procids + 1
      do while ( workleft /= 0 )
         if ( procids == 1 ) procids = npr_y
         ydist(procids) = ydist(procids) + 1
         workleft = workleft - 1  
         if ( workleft /= 0 ) then
            ydist(procidn) = ydist(procidn) + 1
            workleft = workleft - 1
         end if
         procidn = procidn + 1
         procids = procids - 1
      end do
   end if

   ! Safety checks:
   if (sum(ydist) /= plat) then
      write(iulog,*) sub//': ERROR: sum(ydist)=', sum(ydist),' not equal to plat'
      call endrun(sub//': ERROR: sum(ydist) not equal to plat')
   end if
   if (workleft/=0) then
      write(iulog,*) sub//': ERROR: workleft for ydist not zero. workleft=', workleft
      call endrun(sub//': ERROR: workleft for ydist not zero')
   end if

   ! Set the CAM data structures
   lat  = 0
   do procid = 0, npr_y-1
      cut(1,procid) = lat + 1
      lat = lat + ydist(procid+1)
      cut(2,procid) = lat
      nlat_p(procid) = ydist(procid+1)

      if (masterproc) then
         write(iulog,*) 'nlat_p(',procid,') = ', nlat_p(procid)
      end if

      if (myid_y == procid) then
         beglat  = cut(1,myid_y)
         endlat  = cut(2,myid_y)
         numlats = ydist(procid+1)
      end if
   end do

   do k = 1, npr_z-1
      do j = 0, npr_y-1
         procid = j + k*npr_y
         cut(1,procid) = cut(1,j)
         cut(2,procid) = cut(2,j)
         nlat_p(procid) = nlat_p(j)
      end do
   end do

   proc(:) = 0
   do procid = 0, npr_y*npr_z-1

      ! Determine which processor is responsible for the defined latitudes
      do lat = cut(1,procid), cut(2,procid)
         proc(lat) = procid
      end do
   end do

   ! Compute Z partition for YZ decomposition

   allocate(zdist((npes-1)/npr_y+1))

   zdist(:) = 0

   vert = plev / npr_z
   workleft = plev - vert * npr_z
   if (vert < 1) then
      call endrun(sub//': ERROR: less than 1 vertical levels per subdomain')
   end if

   do procid = 1, npr_z
      zdist(procid) = vert
   end do

   if (workleft /= 0) then
      procids = (npr_z+1) / 2
      procidn = procids + 1
      do while ( workleft /= 0 )
         if (procids == 1) procids = npr_z
         zdist(procids) = zdist(procids) + 1
         workleft = workleft - 1  
         if (workleft /= 0) then
            zdist(procidn) = zdist(procidn) + 1
            workleft = workleft - 1
         end if
         procidn = procidn + 1
         procids = procids - 1
      end do
   end if

   ! Safety checks:
   if (sum(zdist) /= plev) then
      write(iulog,*) sub//': ERROR: sum(zdist)=', sum(zdist),' not equal to plev'
      call endrun(sub//': ERROR: sum(zdist) not equal to plev')
   endif
   if (workleft/=0) then
      write(iulog,*) sub//': ERROR: workleft for zdist not zero. workleft=', workleft
      call endrun(sub//': ERROR: workleft for zdist not zero')
   end if

   ! Compute local limits
   call locallimits(myid_z, zdist, beglev, endlev)

   ! Auxiliary processes only
   if (iam >= npes_yz) then
      beglat   = 1
      endlat   = 0
      numlats  = 0
      beglev   = 1
      endlev   = 0
   end if

   grid%jfirst   = beglat
   grid%jlast    = endlat
   grid%kfirst   = beglev
   grid%klast    = endlev
   if (endlev == plev) then
      grid%klastp = plev+1
   else
      grid%klastp = endlev
   end if

   ! Compute X partition for XY decomposition

   allocate(xdistxy(nprxy_x))
   xdistxy(:) = 0

   lonn = plon / nprxy_x
   workleft = plon - lonn * nprxy_x
   if (lonn < 3) then
      call endrun(sub//': ERROR: less than 3 longitudes per XY subdomain')
   end if

   do procid = 1, nprxy_x
      xdistxy(procid) = lonn
   enddo

   if (workleft /= 0) then
      procids = (nprxy_x+1) / 2
      procidn = procids + 1
      do while (workleft /= 0)
         if (procids == 1) procids = nprxy_x
         xdistxy(procids) = xdistxy(procids) + 1
         workleft = workleft - 1  
         if (workleft /= 0) then
            xdistxy(procidn) = xdistxy(procidn) + 1
            workleft = workleft - 1
         end if
         procidn = procidn + 1
         procids = procids - 1
      end do
   end if

   ! Safety check:
   if ( sum(xdistxy) /= plon ) then
      write(iulog,*) sub//': ERROR: sum(xdistxy)=', sum(xdistxy),' not equal to plon'
      call endrun(sub//': ERROR: sum(xdistxy) not equal to plon ')
   end if
   if (workleft/=0) then
      write(iulog,*) sub//': ERROR: workleft for xdistxy not zero.  workleft=',workleft
      call endrun(sub//': ERROR: workleft for xdistxy not zero')
   end if

   ! Compute local limits
   call locallimits(myidxy_x, xdistxy, beglonxy, endlonxy)

   ! Compute global array for use in dyn_grid module
   allocate (lonrangexy(2,nprxy_x))
   lonrangexy(1,1) = 1
   lonrangexy(2,1) = xdistxy(1)
   do procid = 2, nprxy_x
      lonrangexy(1,procid) = lonrangexy(2,procid-1) + 1
      lonrangexy(2,procid) = lonrangexy(1,procid) + xdistxy(procid) - 1
   end do

   ! Compute Y partition for XY decomposition

   allocate(ydistxy((npes-1)/nprxy_x+1))
   ydistxy(:) = 0

   lat = plat / nprxy_y
   workleft = plat - lat * nprxy_y
   if (lat < 3) then
      call endrun(sub//': ERROR: less than 3 latitudes per XY subdomain')
   end if

   do procid = 1, nprxy_y
      ydistxy(procid) = lat
   end do

   if (workleft /= 0) then
      procids = (nprxy_y+1) / 2
      procidn = procids + 1
      do while (workleft /= 0)
         if (procids == 1) procids = nprxy_y
         ydistxy(procids) = ydistxy(procids) + 1
         workleft = workleft - 1  
         if (workleft /= 0) then
            ydistxy(procidn) = ydistxy(procidn) + 1
            workleft = workleft - 1
         end if
         procidn = procidn + 1
         procids = procids - 1
      end do
   end if

   ! Safety check:
   if (sum(ydistxy) /= plat) then
      write(iulog,*) sub//': ERROR: sum(ydistxy)=', sum(ydistxy),' not equal to plat'
      call endrun(sub//': ERROR: sum(ydistxy) not equal to plat')
   end if
   if (workleft/=0) then
      write(iulog,*) sub//': ERROR: workleft for ydistxy not zero.  workleft=',workleft
      call endrun(sub//': ERROR: workleft for ydistxy not zero')
   end if

   ! Compute local limits
   call locallimits(myidxy_y, ydistxy, beglatxy, endlatxy)

   ! Auxiliary processes only
   if (iam >= npes_xy) then
      beglonxy = 1
      endlonxy = 0
      beglatxy = 1
      endlatxy = 0
   end if

   grid%ifirstxy = beglonxy
   grid%ilastxy  = endlonxy
   grid%jfirstxy = beglatxy
   grid%jlastxy  = endlatxy

   ! Compute global array for use in dyn_grid module
   allocate (latrangexy(2,nprxy_y))
   latrangexy(1,1) = 1
   latrangexy(2,1) = ydistxy(1)
   do procid = 2, nprxy_y
      latrangexy(1,procid) = latrangexy(2,procid-1) + 1
      latrangexy(2,procid) = latrangexy(1,procid) + ydistxy(procid) - 1
   end do

   deallocate(ydist)
   deallocate(zdist)

   deallocate(xdistxy)
   deallocate(ydistxy)

   ! Calculate the ghost region sizes for the SPMD version (tricky stuff)
   grid%ng_c = 2                    ! Avoid the case where ng_c = 1
   grid%ng_d = min( abs(jord), 3)   ! SJL: number of max ghost latitudes
   grid%ng_d = max( grid%ng_d, 2)
   grid%ng_s = max( grid%ng_c+1, grid%ng_d )

   ! Define imxy, jmxy, jmyz, kmyz from beglonxy, endlonxy, etc.
   allocate(tmp(npes), imxy(nprxy_x), jmxy(nprxy_y), jmyz(npr_y), kmyz(npr_z))

   tmp = 0
   tmp(gid+1) = endlonxy - beglonxy + 1
   call parcollective( mpicom, maxop, npes, tmp )
   imxy(1:nprxy_x) = tmp(1:nprxy_x)

   tmp = 0
   tmp(gid+1) = endlatxy - beglatxy + 1
   call parcollective( mpicom, maxop, npes, tmp )
   do k = 1, nprxy_y
      jmxy(k) = tmp((k-1)*nprxy_x+1)
   end do

   tmp = 0
   tmp(gid+1)   = grid%jlast - grid%jfirst + 1
   call parcollective( mpicom, maxop, npes, tmp )
   jmyz(1:grid%npr_y) = tmp(1:grid%npr_y)

   tmp = 0
   tmp(gid+1)   = grid%klast - grid%kfirst + 1
   call parcollective( mpicom, maxop, npes, tmp )
   do k = 1, grid%npr_z
      kmyz(k) = tmp((k-1)*grid%npr_y+1)
   end do

   ! set up the PILGRIM communications
   call spmd_vars_init(imxy, jmxy, jmyz, kmyz, grid)

   deallocate(tmp, imxy, jmxy, jmyz, kmyz)

#endif
end subroutine spmdinit_dyn

!========================================================================================

subroutine spmdbuf

end subroutine spmdbuf

!========================================================================================

subroutine compute_gsfactors(numperlat, numtot, numperproc, displs)

   ! Compute arguments for gatherv, scatterv

   ! arguments
   integer, intent(in)  :: numperlat            ! number of elements per latitude
   integer, intent(out) :: numtot               ! total number of elements (to send or recv)
   integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
   integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements

#ifdef SPMD
   ! Local variables
   integer :: p
   !---------------------------------------------------------------------------
   
   numtot = numperlat*numlats
   
   do p = 0, npes-1
      numperproc(p) = numperlat*nlat_p(p)
   end do
     
   displs(:) = 0
   do p = 1, npr_y-1
      displs(p) = displs(p-1) + numperproc(p-1)
   end do

   if (npr_z > 1) then
      do p = 1, npr_z-1
         displs(p*npr_y:(p+1)*npr_y-1) = displs(0:npr_y-1)
      enddo
   endif
#endif

end subroutine compute_gsfactors

!========================================================================================

subroutine locallimits(myidxy, distxy, begdimxy, enddimxy)

   integer, intent(in)  :: myidxy
   integer, intent(in)  :: distxy(:)
   integer, intent(out) :: begdimxy
   integer, intent(out) :: enddimxy
      
   integer :: procid
      
   begdimxy = 1
   enddimxy = distxy(1)
      
   do procid = 1, myidxy
      begdimxy = enddimxy + 1
      enddimxy = begdimxy + distxy(procid+1) - 1
   end do
end subroutine locallimits

!========================================================================================

end module spmd_dyn

