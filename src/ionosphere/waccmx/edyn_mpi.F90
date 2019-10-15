module edyn_mpi
  use shr_kind_mod  ,only: r8 => shr_kind_r8
  use cam_logfile   ,only: iulog
  use cam_abortutils,only: endrun

  use edyn_geogrid  ,only: nlon,nlat
  use edyn_maggrid  ,only: nmlonp1,nmlat,nmlath,nmlev  ! note nmlev is not a parameter
  use spmd_utils    ,only: masterproc
  use mpi	    ,only: mpi_comm_size, mpi_comm_rank, MPI_PROC_NULL, mpi_comm_split, &
			   MPI_INTEGER, MPI_STATUS_SIZE, mpi_wait, &
			   MPI_REAL8, MPI_SUCCESS, MPI_SUM, &
			   MPI_Comm_rank

  implicit none
  private
 
  public :: array_ptr_type,switch_model_format,mp_geo_halos,mp_pole_halos,mlon0,mlon1,omlon1, &
            mlat0,mlat1,mlev0,mlev1,mytid,lon0,lon1,lat0,lat1,lev0,lev1,mp_mag_halos,mp_scatter_phim, &
            mp_mageq,mp_mageq_jpm1,mp_magpole_2d,mp_mag_foldhem,mp_mag_periodic_f2d,mp_gather_edyn, &
            mp_mageq_jpm3,mp_mag_jslot,mp_magpoles,ixfind,mp_magpole_3d,ntask,ntaski,ntaskj,tasks, &
            nmagtaski,nmagtaskj,setpoles, mp_gatherlons_f3d, mytidi, mp_scatterlons_f3d, mp_exchange_tasks, &
            mp_distribute_mag, mp_distribute_geo, mp_init



!
! Number of MPI tasks and current task id (geo or mag):
!
  integer :: &
    ntask,   & ! number of mpi tasks
    mytid      ! my task id
!
! Geographic subdomains for current task:
!
  integer ::     &
    ntaski,      & ! number of tasks in lon dimension
    ntaskj,      & ! number of tasks in lat dimension
    mytidi,      & ! i coord for current task in task table
    mytidj,      & ! j coord for current task in task table
    lat0,lat1,   & ! first and last lats for each task
    lon0,lon1,   & ! first and last lons for each task
    lev0,lev1,   & ! first and last levs for each task (not distributed)
    mxlon,mxlat    ! max number of subdomain lon,lat points among all tasks
!
! Magnetic subdomains for current task:
!
  integer ::     &
    nmagtaski,   & ! number of tasks in mag lon dimension
    nmagtaskj,   & ! number of tasks in mag lat dimension
    magtidi,     & ! i coord for current task in task table
    magtidj,     & ! j coord for current task in task table
    mlat0,mlat1, & ! first and last mag lats for each task
    mlon0,mlon1, & ! first and last mag lons for each task
    omlon1,      & ! last mag lons for each task to remove periodic point from outputs
    mlev0,mlev1, & ! first and last mag levs (not distributed)
    mxmaglon,    & ! max number of mag subdomain lon points among all tasks
    mxmaglat       ! max number of mag subdomain lat points among all tasks

  integer,allocatable,save :: &
    itask_table_geo(:,:),     & ! 2d table of tasks on geographic grid (i,j)
    itask_table_mag(:,:)        ! 2d table of tasks on mag grid (i,j)

  integer :: cols_comm  ! communicators for each task column
  integer :: rows_comm  ! communicators for each task row
!
! Task type: subdomain information for all tasks, known by all tasks:
!
  type task
    integer :: mytid          ! task id
!
! Geographic subdomains in task structure:
    integer :: mytidi         ! task coord in longitude dimension of task table
    integer :: mytidj         ! task coord in latitude  dimension of task table
    integer :: nlats          ! number of latitudes  calculated by this task
    integer :: nlons          ! number of longitudes calculated by this task
    integer :: lat0,lat1      ! first and last latitude  indices
    integer :: lon0,lon1      ! first and last longitude indices
!
! Magnetic subdomains in task structure:
    integer :: magtidi        ! task coord in mag longitude dimension of task table
    integer :: magtidj        ! task coord in mag latitude  dimension of task table
    integer :: nmaglats       ! number of mag latitudes  calculated by this task
    integer :: nmaglons       ! number of mag longitudes calculated by this task
    integer :: mlat0,mlat1    ! first and last latitude  indices
    integer :: mlon0,mlon1    ! first and last longitude indices
  end type task
!
! type(task) :: tasks(ntask) will be made available to all tasks
! (so each task has information about all tasks)
!
  type(task),allocatable,save :: tasks(:)
!
! Conjugate points in mag subdomains, for mp_mag_foldhem
!
  integer,allocatable,dimension(:),save ::   & ! (ntask)
    nsend_south,       & ! number of south lats to send to north (each task)
    nrecv_north          ! number of north lats to send to south (each task)
  integer,allocatable,dimension(:,:),save :: & ! (mxlats,ntask)
    send_south_coords, & ! south j lats to send to north
    recv_north_coords    ! north j lats to recv from south

  type array_ptr_type
    real(r8),pointer :: ptr(:,:,:) ! (k,i,j)
  end type array_ptr_type

  integer, protected :: mpi_comm_edyn = -9999

  logical, parameter :: debug = .false.

  contains
!-----------------------------------------------------------------------
  subroutine mp_init( mpi_comm )
!
! Initialize MPI, and allocate task table.
!
    integer, intent(in) :: mpi_comm

    integer :: ier

    mpi_comm_edyn = mpi_comm

    call mpi_comm_size(mpi_comm_edyn,ntask,ier)
    call mpi_comm_rank(mpi_comm_edyn,mytid,ier)
!
! Allocate array of task structures:
!
    allocate(tasks(0:ntask-1),stat=ier)
    if (ier /= 0) then
      write(iulog,"('>>> mp_init: error allocating tasks(',i3,')')") ntask
      call endrun('edyn_mpi mp_init')
    endif
  end subroutine mp_init
!-----------------------------------------------------------------------
  subroutine mp_distribute_geo(lonndx0,lonndx1,latndx0,latndx1,levndx0,levndx1, ntaski_in,ntaskj_in)
!
! Args:
    integer, intent(in) :: lonndx0,lonndx1,latndx0,latndx1,levndx0,levndx1, ntaski_in,ntaskj_in
!
! Local:
    integer :: i,j,n,irank,ier,tidrow,nj,ni
!
! Define all task structures with current task values
! (redundant for alltoall):
! Use WACCM subdomains:
!
    lon0 = lonndx0 ; lon1 = lonndx1
    lat0 = latndx0 ; lat1 = latndx1
    lev0 = levndx0 ; lev1 = levndx1

    ntaski = ntaski_in
    ntaskj = ntaskj_in
!
! Allocate and set 2d table of tasks:
!
    allocate(itask_table_geo(-1:ntaski,-1:ntaskj),stat=ier)
    if (ier /= 0) then
      write(iulog,"('>>> Error allocating itable: ntaski,j=',2i4)") ntaski,ntaskj
      call endrun('itask_table_geo')
    endif
    itask_table_geo(:,:) = MPI_PROC_NULL

    irank = 0
    do j =  0,ntaskj-1
      do i = 0,ntaski-1
        itask_table_geo(i,j) = irank
        if (mytid == irank) then
          mytidi = i
          mytidj = j
        endif
        irank = irank+1
      enddo
!
! Tasks are periodic in longitude:
! (this is not done in tiegcm, but here sub mp_geo_halos depends on it)
!
      itask_table_geo(-1,j) = itask_table_geo(ntaski-1,j)
      itask_table_geo(ntaski,j) = itask_table_geo(0,j)

    enddo ! j=0,ntaskj-1

    if (debug.and.masterproc) then
      write(iulog,"('mp_distribute_geo: mytid=',i4,' ntaski,j=',2i4,' mytidi,j=',2i4,&
                    ' lon0,1=',2i4,' lat0,1=',2i4,' lev0,1=',2i4)") &
                  mytid,ntaski,ntaskj,mytidi,mytidj,lon0,lon1,lat0,lat1,lev0,lev1
!
! Print table to stdout, including -1,ntaski:
!
      write(iulog,"(/,'ntask=',i3,' ntaski=',i2,' ntaskj=',i2,' Geo Task Table:')") &
        ntask,ntaski,ntaskj
      do j=-1,ntaskj
        write(iulog,"('j=',i3,' itask_table_geo(:,j)=',100i3)") j,itask_table_geo(:,j)
      enddo
    endif
!
! Calculate start and end indices in lon,lat dimensions for each task:
! For WACCM: do not call distribute_1d - lon0,1, lat0,1 are set from
!   waccm grid above.
!
!     call distribute_1d(1,nlon,ntaski,mytidi,lon0,lon1)
!     call distribute_1d(1,nlat,ntaskj,mytidj,lat0,lat1)

      nj = lat1-lat0+1 ! number of latitudes for this task
      ni = lon1-lon0+1 ! number of longitudes for this task
!
! Report my stats to stdout:
!     write(iulog,"(/,'mytid=',i3,' mytidi,j=',2i3,' lat0,1=',2i3,' (',i2,') lon0,1=',2i3,' (',i2,') ncells=',i4)") &
!       mytid,mytidi,mytidj,lat0,lat1,nj,lon0,lon1,ni
!
! Define all task structures with current task values
! (redundant for alltoall):
!
    do n=0,ntask-1
      tasks(n)%mytid  = mytid
      tasks(n)%mytidi = mytidi
      tasks(n)%mytidj = mytidj
      tasks(n)%nlats  = nj
      tasks(n)%nlons  = ni
      tasks(n)%lat0   = lat0
      tasks(n)%lat1   = lat1
      tasks(n)%lon0   = lon0
      tasks(n)%lon1   = lon1
    enddo
!
! All tasks must have at least 4 longitudes:
!
    do n=0,ntask-1

    if (debug.and.masterproc) then
      write(iulog,"('mp_distribute_geo: n=',i3,' tasks(n)%nlons=',i3,' tasks(n)%nlats=',i3)") &
        n,tasks(n)%nlons,tasks(n)%nlats
    endif

      if (tasks(n)%nlons < 4) then
        write(iulog,"('>>> mp_distribute_geo: each task must carry at least 4 longitudes. task=',i4,' nlons=',i4)") &
          n,tasks(n)%nlons
        call endrun('edyn_mpi: nlons per task')
      endif
    enddo
!
! Create sub-communicators for each task row (used by mp_geopole_3d):
!
!   call mpi_comm_split(mpi_comm_edyn,mod(mytid,ntaskj),mytid,rows_comm,ier)
!   call MPI_Comm_rank(rows_comm,tidrow,ier)

    call mpi_comm_split(mpi_comm_edyn,mytidj,mytid,rows_comm,ier)
    call MPI_Comm_rank(rows_comm,tidrow,ier)

    if (debug.and.masterproc) then
      write(iulog,"('mp_distribute_geo: ntaskj=',i3,' tidrow=',i3)") &
        ntaskj,tidrow
    endif

  end subroutine mp_distribute_geo
!-----------------------------------------------------------------------
  subroutine mp_distribute_mag
!
! Local:
    integer :: i,j,n,irank,ier,tidcol,nj,ni,ncells
!
! Number of tasks in mag lon,lat same as geo grid:
! Also true for WACCM processor distribution.
!
    nmagtaski = ntaski
    nmagtaskj = ntaskj
!
! Vertical dimension is not distributed:
    mlev0 = 1
    mlev1 = nmlev
!
! Allocate and set 2d table of tasks:
    allocate(itask_table_mag(-1:nmagtaski,-1:nmagtaskj),stat=ier)
    if (ier /= 0) then
      write(iulog,"('>>> Error allocating itable: nmagtaski,j=',2i3)") &
        nmagtaski,nmagtaskj
      call endrun('itask_table_mag')
    endif
    itask_table_mag(:,:) = MPI_PROC_NULL
    irank = 0
    do j =  0,nmagtaskj-1
      do i = 0,nmagtaski-1
        itask_table_mag(i,j) = irank
        if (mytid == irank) then
          magtidi = i
          magtidj = j
        endif
        irank = irank+1
      enddo
!
! Tasks are periodic in longitude:
!
      itask_table_mag(-1,j) = itask_table_mag(nmagtaski-1,j)
      itask_table_mag(nmagtaski,j) = itask_table_mag(0,j)
    enddo

    if (debug.and.masterproc) then
!
! Print table to stdout:
      write(iulog,"(/,'ntask=',i3,' nmagtaski=',i2,' nmagtaskj=',i2,' Mag Task Table:')") &
        ntask,nmagtaski,nmagtaskj
      do j=-1,nmagtaskj
        write(iulog,"('j=',i3,' itask_table_mag(:,j)=',100i3)") j,itask_table_mag(:,j)
      enddo
    endif
!
! Calculate start and end indices in mag lon,lat dimensions for each task:
!
    call distribute_1d(1,nmlonp1,nmagtaski,magtidi,mlon0,mlon1)
    call distribute_1d(1,nmlat  ,nmagtaskj,magtidj,mlat0,mlat1)

    omlon1=mlon1
    if (omlon1 == nmlonp1) omlon1=omlon1-1

    nj = mlat1-mlat0+1 ! number of mag latitudes for this task
    ni = mlon1-mlon0+1 ! number of mag longitudes for this task
    ncells = nj*ni     ! total number of grid cells for this task

    if (debug.and.masterproc) then
!
! Report my stats to stdout:
      write(iulog,"(/,'mytid=',i3,' magtidi,j=',2i3,' mlat0,1=',2i3,' (',i2,') mlon0,1=',2i3,' (',i2,') ncells=',i4)") &
        mytid,magtidi,magtidj,mlat0,mlat1,nj,mlon0,mlon1,ni,ncells
    endif
!
! Define all task structures with current task values
! (redundant for alltoall):
!
    do n=0,ntask-1
      tasks(n)%magtidi = magtidi
      tasks(n)%magtidj = magtidj
      tasks(n)%nmaglats  = nj
      tasks(n)%nmaglons  = ni
      tasks(n)%mlat0   = mlat0
      tasks(n)%mlat1   = mlat1
      tasks(n)%mlon0   = mlon0
      tasks(n)%mlon1   = mlon1
    enddo
!
! All tasks must have at least 4 longitudes:
    do n=0,ntask-1
      if (tasks(n)%nmaglons < 4) then
        write(iulog,"('>>> mp_distribute_mag: each task must carry at least 4 longitudes. task=',i4,' nmaglons=',i4)") &
          n,tasks(n)%nmaglons
        call endrun('edyn_mpi: nmaglons per task')
      endif
    enddo
!
! Create subgroup communicators for each task column:
! These communicators will be used by sub mp_mag_jslot (mpi.F).
!
    call mpi_comm_split(mpi_comm_edyn,mod(mytid,nmagtaski),mytid,cols_comm,ier)
    call MPI_Comm_rank(cols_comm,tidcol,ier)

    if (debug.and.masterproc) then
      write(iulog,"('mp_distribute_mag: nmagtaski=',i3,' mod(mytid,nmagtaski)=',i3,' tidcol=',i3)") &
        nmagtaski,mod(mytid,nmagtaski),tidcol
    endif
   
  end subroutine mp_distribute_mag
!-----------------------------------------------------------------------
  subroutine distribute_1d(n1,n2,nprocs,myrank,istart,iend)
!
! Distribute work across a 1d vector(n1->n2) to nprocs.
! Return start and end indices for proc myrank.
!
! Args:
    integer,intent(in)  :: n1,n2,nprocs,myrank
    integer,intent(out) :: istart,iend
!
! Local:
    integer :: lenproc,iremain,n
!
    n = n2-n1+1
    lenproc = n/nprocs
    iremain = mod(n,nprocs)
    istart = n1 + myrank*lenproc + min(myrank,iremain)
    iend = istart+lenproc-1
    if (iremain > myrank) iend = iend+1
  end subroutine distribute_1d
!-----------------------------------------------------------------------
  subroutine mp_exchange_tasks(iprint)
!
! Args:
    integer,intent(in) :: iprint
!
! Local:
! itasks_send(len_task_type,ntask) will be used to send tasks(:) info
!   to all tasks (directly passing mpi derived data types is reportedly
!   not stable, or not available until MPI 2.x).
!
    integer :: n,ier
    integer,parameter :: len_task_type = 17 ! see type task above
    integer,allocatable,save :: &
      itasks_send(:,:), & ! send buffer
      itasks_recv(:,:)    ! send buffer
!
! Pack tasks(mytid) into itasks_send:
    allocate(itasks_send(len_task_type,0:ntask-1),stat=ier)
    if (ier /= 0) then
      write(iulog,"(i4,i4)") '>>> Error allocating itasks_send: len_task_type=',&
        len_task_type,' ntask=',ntask
    endif
    allocate(itasks_recv(len_task_type,0:ntask-1),stat=ier)
    if (ier /= 0) then
      write(iulog,"(i4,i4)") '>>> Error allocating itasks_recv: len_task_type=',&
        len_task_type,' ntask=',ntask
    endif
    do n=0,ntask-1
      itasks_send(1,n) = tasks(mytid)%mytid

      itasks_send(2,n) = tasks(mytid)%mytidi
      itasks_send(3,n) = tasks(mytid)%mytidj
      itasks_send(4,n) = tasks(mytid)%nlats
      itasks_send(5,n) = tasks(mytid)%nlons
      itasks_send(6,n) = tasks(mytid)%lat0
      itasks_send(7,n) = tasks(mytid)%lat1
      itasks_send(8,n) = tasks(mytid)%lon0
      itasks_send(9,n) = tasks(mytid)%lon1

      itasks_send(10,n) = tasks(mytid)%magtidi
      itasks_send(11,n) = tasks(mytid)%magtidj
      itasks_send(12,n) = tasks(mytid)%nmaglats
      itasks_send(13,n) = tasks(mytid)%nmaglons
      itasks_send(14,n) = tasks(mytid)%mlat0
      itasks_send(15,n) = tasks(mytid)%mlat1
      itasks_send(16,n) = tasks(mytid)%mlon0
      itasks_send(17,n) = tasks(mytid)%mlon1
    enddo
!
! Send itasks_send and receive itasks_recv:
    call mpi_alltoall(itasks_send,len_task_type,MPI_INTEGER,&
                      itasks_recv,len_task_type,MPI_INTEGER,&
                      mpi_comm_edyn,ier)
    if (ier /= 0) &
      call handle_mpi_err(ier,'edyn_mpi: mpi_alltoall to send/recv itasks')
!
! Unpack itasks_recv into tasks(n)
!
    do n=0,ntask-1
      tasks(n)%mytid  = itasks_recv(1,n)

      tasks(n)%mytidi = itasks_recv(2,n)
      tasks(n)%mytidj = itasks_recv(3,n)
      tasks(n)%nlats  = itasks_recv(4,n)
      tasks(n)%nlons  = itasks_recv(5,n)
      tasks(n)%lat0   = itasks_recv(6,n)
      tasks(n)%lat1   = itasks_recv(7,n)
      tasks(n)%lon0   = itasks_recv(8,n)
      tasks(n)%lon1   = itasks_recv(9,n)

      tasks(n)%magtidi   = itasks_recv(10,n)
      tasks(n)%magtidj   = itasks_recv(11,n)
      tasks(n)%nmaglats  = itasks_recv(12,n)
      tasks(n)%nmaglons  = itasks_recv(13,n)
      tasks(n)%mlat0   = itasks_recv(14,n)
      tasks(n)%mlat1   = itasks_recv(15,n)
      tasks(n)%mlon0   = itasks_recv(16,n)
      tasks(n)%mlon1   = itasks_recv(17,n)
!
! Report to stdout:
!
      if (n==mytid.and.iprint > 0) then
        write(iulog,"(/,'Task ',i3,':')") n
        write(iulog,"(/,'Subdomain on geographic grid:')")
        write(iulog,"('tasks(',i3,')%mytid =',i3)") n,tasks(n)%mytid
        write(iulog,"('tasks(',i3,')%mytidi=',i3)") n,tasks(n)%mytidi
        write(iulog,"('tasks(',i3,')%mytidj=',i3)") n,tasks(n)%mytidj
        write(iulog,"('tasks(',i3,')%nlats =',i3)") n,tasks(n)%nlats
        write(iulog,"('tasks(',i3,')%nlons =',i3)") n,tasks(n)%nlons
        write(iulog,"('tasks(',i3,')%lat0  =',i3)") n,tasks(n)%lat0
        write(iulog,"('tasks(',i3,')%lat1  =',i3)") n,tasks(n)%lat1
        write(iulog,"('tasks(',i3,')%lon0  =',i3)") n,tasks(n)%lon0
        write(iulog,"('tasks(',i3,')%lon1  =',i3)") n,tasks(n)%lon1
        write(iulog,"('Number of geo subdomain grid points = ',i6)") &
          tasks(n)%nlons * tasks(n)%nlats
        write(iulog,"(/,'Subdomain on geomagnetic grid:')")
        write(iulog,"('tasks(',i3,')%magtidi=',i3)") n,tasks(n)%magtidi
        write(iulog,"('tasks(',i3,')%magtidj=',i3)") n,tasks(n)%magtidj
        write(iulog,"('tasks(',i3,')%nmaglats =',i3)") n,tasks(n)%nmaglats
        write(iulog,"('tasks(',i3,')%nmaglons =',i3)") n,tasks(n)%nmaglons
        write(iulog,"('tasks(',i3,')%mlat0  =',i3)") n,tasks(n)%mlat0
        write(iulog,"('tasks(',i3,')%mlat1  =',i3)") n,tasks(n)%mlat1
        write(iulog,"('tasks(',i3,')%mlon0  =',i3)") n,tasks(n)%mlon0
        write(iulog,"('tasks(',i3,')%mlon1  =',i3)") n,tasks(n)%mlon1
        write(iulog,"('Number of mag subdomain grid points = ',i6)") &
          tasks(n)%nmaglons * tasks(n)%nmaglats
      endif
    enddo
!
! Release locally allocated space:
    deallocate(itasks_send)
    deallocate(itasks_recv)
!
! mxlon,mxlat are maximum number of lons,lats owned by all tasks:
    mxlon = -9999
    do n=0,ntask-1
      if (tasks(n)%nlons > mxlon) mxlon = tasks(n)%nlons
    enddo
    mxlat = -9999
    do n=0,ntask-1
      if (tasks(n)%nlats > mxlat) mxlat = tasks(n)%nlats
    enddo
!
! mxmaglon,mxmaglat are maximum number of mag lons,lats owned by all tasks:
    mxmaglon = -9999
    do n=0,ntask-1
      if (tasks(n)%nmaglons > mxmaglon) mxmaglon = tasks(n)%nmaglons
    enddo
    mxmaglat = -9999
    do n=0,ntask-1
      if (tasks(n)%nmaglats > mxmaglat) mxmaglat = tasks(n)%nmaglats
    enddo
!
! Find conjugate points for folding hemispheres:
    call conjugate_points

  end subroutine mp_exchange_tasks
!-----------------------------------------------------------------------
  subroutine mp_mageq(fin,fout,nf,mlon0,mlon1,mlat0,mlat1,nmlev)
!
! Each task needs values of conductivities and adotv1,2 fields at the
! at the mag equator for its longitude subdomain (and all levels), for
! the fieldline integrations.
!
! On input, fin is ped_mag, hal_mag, adotv1_mag, adotv2_mag
!   on (i,j,k) magnetic subdomain.
! On output, fout(mlon0:mlon1,nmlev,nf) is ped_meq, hal_meq, adotv1_meq,
!   adotv2_meq at mag equator at longitude subdomain and all levels.
!
! Args:
    integer :: mlon0,mlon1,mlat0,mlat1,nmlev,nf
    real(r8),intent(in)  :: fin (mlon0:mlon1,mlat0:mlat1,nmlev,nf)
    real(r8),intent(out) :: fout(mlon0:mlon1,nmlev,nf)
!
! Local:
    real(r8) :: & ! mpi buffers
      sndbuf(mxmaglon,nmlev,nf),             & ! mxmaglon,nmlev,nf
      rcvbuf(mxmaglon,nmlev,nf)                ! mxmaglon,nmlev,nf
    integer :: i,j,n,itask,ier,len,jlateq,ireqsend,ireqrecv
    integer :: irstat(MPI_STATUS_SIZE)      ! mpi receive status
    logical :: have_eq

    sndbuf = 0._r8
    rcvbuf = 0._r8
    len = mxmaglon*nmlev*nf
!
! If mag equator is in current subdomain, load it into sndbuf
! and send to other tasks in my task column (mytidi)
!
    jlateq = (nmlat+1)/2   ! lat index of mag equator (49)
    have_eq = .false.
    do j=mlat0,mlat1
      if (j == jlateq) then ! load send buffer w/ data at equator
        have_eq = .true.
        do i=mlon0,mlon1
          sndbuf(i-mlon0+1,:,:) = fin(i,j,:,:)
        enddo
!
! Send mag equator data to other tasks in my task column (mytidi):
        do itask=0,ntask-1
          if (itask /= mytid.and.tasks(itask)%mytidi==mytidi) then
            call mpi_isend(sndbuf,len,MPI_REAL8,itask,1, &
              mpi_comm_edyn,ireqsend,ier)
            if (ier /= 0) call handle_mpi_err(ier,'mp_mageq isend')
            call mpi_wait(ireqsend,irstat,ier)
          endif ! another task in mytidi
        enddo ! itask=0,ntask-1
      endif ! j==jlateq
    enddo ! j=mlat0,mlat1
!
! Receive by other tasks in the sending task's column:
    fout = 0._r8
    if (.not.have_eq) then ! find task to receive from
      do itask=0,ntask-1
        do j=tasks(itask)%mlat0,tasks(itask)%mlat1
          if (j == jlateq.and.tasks(itask)%mytidi==mytidi) then
            call mpi_irecv(rcvbuf,len,MPI_REAL8,itask,1, &
              mpi_comm_edyn,ireqrecv,ier)
            if (ier /= 0) call handle_mpi_err(ier,'mp_mageq irecv')
            call mpi_wait(ireqrecv,irstat,ier)
            do n=1,nf
              do i=mlon0,mlon1
                fout(i,:,n) = rcvbuf(i-mlon0+1,:,n)
              enddo
            enddo
          endif ! itask has mag eq and is in my column (sending task)
        enddo ! scan itask latitudes
      enddo ! task table search
!
! If I am the sending task, set fout to equator values of input array:
    else
      do n=1,nf
        do i=mlon0,mlon1
          fout(i,:,n) = fin(i,jlateq,:,n)
        enddo
      enddo
    endif ! I am receiving or sending task
  end subroutine mp_mageq
!-----------------------------------------------------------------------
  subroutine mp_mageq_jpm1(f,mlon0,mlon1,mlat0,mlat1,nmlonp1,feq_jpm1,nf)
!
! All tasks need data at mag latitudes equator-1, equator+1 at global
!   longitudes.
! On input: f is 6 fields on mag subdomains: zigm11,zigm22,zigmc,zigm2,rim1,rim2
! On output: feq_jpm1(nmlonp1,2,nf)
!
! Args:
    integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nmlonp1,nf
    real(r8),intent(in) :: f(mlon0:mlon1,mlat0:mlat1,nf)
    real(r8),intent(out) :: feq_jpm1(nmlonp1,2,nf) ! eq-1,eq+1
!
! Local:
    integer :: j,ier,len,jlateq
    real(r8) :: sndbuf(nmlonp1,2,nf)

    sndbuf = 0._r8
    feq_jpm1 = 0._r8
    len = nmlonp1*2*nf
!
! Load send buffer w/ eq +/- 1 for current subdomain
! (redundant to all tasks for alltoall)
!
    jlateq = (nmlat+1)/2
    do j=mlat0,mlat1
       if (j == jlateq+1) then     ! equator+1
          sndbuf(mlon0:mlon1,1,:) = f(mlon0:mlon1,j,:)
       elseif (j == jlateq-1) then ! equator-1
          sndbuf(mlon0:mlon1,2,:) = f(mlon0:mlon1,j,:)
       endif ! j==jlateq
    enddo ! j=mlat0,mlat1
!
! Do the exchange:
!
    call mpi_allreduce( sndbuf(:,:,1:nf), feq_jpm1(:,:,1:nf), len, MPI_REAL8, MPI_SUM, mpi_comm_edyn, ier )
    if ( ier .ne. MPI_SUCCESS ) call handle_mpi_err(ier,'mp_mageq_jpm1 call mpi_allreduce')

!
! Periodic point:
    feq_jpm1(nmlonp1,:,:) = feq_jpm1(1,:,:)

  end subroutine mp_mageq_jpm1
!-----------------------------------------------------------------------
  subroutine mp_mageq_jpm3(f,mlon0,mlon1,mlat0,mlat1,nmlonp1,feq_jpm3,nf)
!
! All tasks need global longitudes at mag latitudes equator,
!   and equator +/- 1,2,3
! On input: f is nf fields on mag subdomains
! On output: feq_jpm3(nmlonp1,-3:3,nf) has global lons at eq, eq +/- 1,2,3
!   2nd dimension of feq_jpm3 (and send/recv buffers) is as follows:
!   +3: eq+3
!   +2: eq+2
!   +1: eq+1
!    0: eq
!   -1: eq-1
!   -2: eq-2
!   -3: eq-3
!
! Args:
    integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nmlonp1,nf
    real(r8),intent(in) :: f(mlon0:mlon1,mlat0:mlat1,nf)
    real(r8),intent(out) :: feq_jpm3(nmlonp1,-3:3,nf)
!
! Local:
    integer :: j,ier,len,jlateq
    integer,parameter :: mxnf=6

    real(r8) :: sndbuf(nmlonp1,-3:3,mxnf)

    if (nf > mxnf) then
      write(iulog,"('>>> mp_mageq_jpm3: nf=',i4,' but cannot be called with greater than mxnf=',i4)") &
        nf,mxnf
      call endrun('mp_mageq_jpm3')
    endif
 
    sndbuf = 0._r8
    feq_jpm3 = 0._r8
    len = nmlonp1*7*nf
!
! Load send buffer w/ eq +/- 3 for current subdomain
!
    jlateq = (nmlat+1)/2
    do j=mlat0,mlat1
       if (j == jlateq-3) then       ! equator-3
          sndbuf(mlon0:mlon1,-3,1:nf) = f(mlon0:mlon1,j,:)
       elseif (j == jlateq-2) then   ! equator-2
          sndbuf(mlon0:mlon1,-2,1:nf) = f(mlon0:mlon1,j,:)
       elseif (j == jlateq-1) then   ! equator-1
          sndbuf(mlon0:mlon1,-1,1:nf) = f(mlon0:mlon1,j,:)
       elseif (j == jlateq) then     ! equator
          sndbuf(mlon0:mlon1,0,1:nf) = f(mlon0:mlon1,j,:)
       elseif (j == jlateq+1) then   ! equator+1
          sndbuf(mlon0:mlon1,1,1:nf) = f(mlon0:mlon1,j,:)
       elseif (j == jlateq+2) then   ! equator+2
          sndbuf(mlon0:mlon1,2,1:nf) = f(mlon0:mlon1,j,:)
       elseif (j == jlateq+3) then   ! equator+3
          sndbuf(mlon0:mlon1,3,1:nf) = f(mlon0:mlon1,j,:)
       endif ! j==jlateq
    enddo ! j=mlat0,mlat1
!
! Do the exchange:
!
    call mpi_allreduce( sndbuf(:,:,1:nf), feq_jpm3(:,:,1:nf), len, MPI_REAL8, MPI_SUM, mpi_comm_edyn, ier )
    if ( ier .ne. MPI_SUCCESS ) call handle_mpi_err(ier,'mp_mageq_jpm3 call mpi_allreduce')

!
! Periodic point:
    feq_jpm3(nmlonp1,:,:) = feq_jpm3(1,:,:)

  end subroutine mp_mageq_jpm3
!-----------------------------------------------------------------------
  subroutine mp_magpole_2d(f,ilon0,ilon1,ilat0,ilat1, &
    nglblon,jspole,jnpole,fpole_jpm2,nf)
!
! Return fpole_jpm2(nglblon,1->4,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2)
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2)
! This can be called with different number of fields nf, but cannot
! be called w/ > mxnf fields. 
!
! Args:
    integer,intent(in) :: ilon0,ilon1,ilat0,ilat1,nglblon,jspole,jnpole,nf
    real(r8),intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nf)
    real(r8),intent(out) :: fpole_jpm2(nglblon,4,nf)
!
! Local:
    integer :: j,ier,len
    integer,parameter :: mxnf=6
    real(r8) :: sndbuf(nglblon,4,mxnf)

    if (nf > mxnf) then
      write(iulog,"('>>> mp_magpole_2d: nf=',i4,' but cannot be called with greater than mxnf=',i4)") &
        nf,mxnf
      call endrun('mp_magpole_2d')
    endif

    sndbuf = 0._r8
    fpole_jpm2 = 0._r8
    len = nglblon*4*nf
!
! Load send buffer with values at poles +/- 2 for current subdomain
!
    do j=ilat0,ilat1
       if (j==jspole+1) then             ! south pole +1
          sndbuf(ilon0:ilon1,1,1:nf) = f(ilon0:ilon1,j,:)
       elseif (j==jspole+2) then         ! south pole +2
          sndbuf(ilon0:ilon1,2,1:nf) = f(ilon0:ilon1,j,:)
       elseif (j==jnpole-1) then         ! north pole -1
          sndbuf(ilon0:ilon1,3,1:nf) = f(ilon0:ilon1,j,:)
       elseif (j==jnpole-2) then         ! north pole -2
          sndbuf(ilon0:ilon1,4,1:nf) = f(ilon0:ilon1,j,:)
       endif
    enddo

!
! Do the exchange:
!
    call mpi_allreduce( sndbuf(:,:,1:nf), fpole_jpm2(:,:,1:nf), len, MPI_REAL8, MPI_SUM, mpi_comm_edyn, ier )
    if ( ier .ne. MPI_SUCCESS ) call handle_mpi_err(ier,'mp_magpole_2d call mpi_allreduce')

  end subroutine mp_magpole_2d
!-----------------------------------------------------------------------
  subroutine mp_magpole_3d(f,ilon0,ilon1,ilat0,ilat1,nlev, nglblon,jspole,jnpole,fpole_jpm2,nf)
!
! Return fpole_jpm2(nglblon,1->4,nlev,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2)
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2)
! This can be called with different number of fields nf, but cannot
! be called w/ > mxnf fields. 
!
! Args:
  integer,intent(in) :: ilon0,ilon1,ilat0,ilat1,nglblon,&
    jspole,jnpole,nf,nlev
  real(r8),intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nlev,nf)
  real(r8),intent(out) :: fpole_jpm2(nglblon,4,nlev,nf)
!
! Local:
  integer :: j,k,ier,len
  integer,parameter :: mxnf=6
  real(r8) :: sndbuf(nglblon,4,nlev,mxnf)

  if (nf > mxnf) then
    write(iulog,"('>>> mp_magpole_3d: nf=',i4,' but cannot be called with greater than mxnf=',i4)") &
      nf,mxnf
    call endrun('mp_magpole_3d')
  endif

  sndbuf = 0._r8
  fpole_jpm2 = 0._r8
  len = nglblon*4*nlev*nf
!
! Load send buffer with values at poles +/- 2 for current subdomain
!
  do j=ilat0,ilat1
     do k=1,nlev
        if (j==jspole+1) then             ! south pole +1
           sndbuf(ilon0:ilon1,1,k,1:nf) = f(ilon0:ilon1,j,k,:)
        elseif (j==jspole+2) then         ! south pole +2
           sndbuf(ilon0:ilon1,2,k,1:nf) = f(ilon0:ilon1,j,k,:)
        elseif (j==jnpole-1) then         ! north pole -1
           sndbuf(ilon0:ilon1,3,k,1:nf) = f(ilon0:ilon1,j,k,:)
        elseif (j==jnpole-2) then         ! north pole -2
           sndbuf(ilon0:ilon1,4,k,1:nf) = f(ilon0:ilon1,j,k,:)
        endif
     enddo
  enddo

!
! Do the exchange:
!
    call mpi_allreduce( sndbuf(:,:,:,1:nf), fpole_jpm2(:,:,:,1:nf), len, MPI_REAL8, MPI_SUM, mpi_comm_edyn, ier )
    if ( ier .ne. MPI_SUCCESS ) call handle_mpi_err(ier,'mp_magpole_3d call mpi_allreduce')

  end subroutine mp_magpole_3d
!-----------------------------------------------------------------------
  subroutine mp_magpoles(f,ilon0,ilon1,ilat0,ilat1,nglblon, jspole,jnpole,fpoles,nf)
!
! Similiar to mp_magpole_2d, but returns global longitudes for
! j==1 and j==nmlat (not for poles +/- 2)
! Return fpoles(nglblon,2,nf) as:
!   1: j = jspole (spole)
!   2: j = jnpole (npole)
! This can be called with different number of fields nf, but cannot
! be called w/ > mxnf fields.
!
! Args:
    integer,intent(in) :: ilon0,ilon1,ilat0,ilat1,nglblon, jspole,jnpole,nf
    real(r8),intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nf)
    real(r8),intent(out) :: fpoles(nglblon,2,nf)
!
! Local:
    integer :: j,ier,len
    real(r8) :: sndbuf(nglblon,2,nf)

    sndbuf = 0._r8
    fpoles = 0._r8
    len = nglblon*2*nf
!
! Load send buffer with values at poles +/- 2 for current subdomain
!
    do j=ilat0,ilat1
       if (j==jspole) then             ! south pole
          sndbuf(ilon0:ilon1,1,1:nf) = f(ilon0:ilon1,j,:)
       elseif (j==jnpole) then         ! npole pole
          sndbuf(ilon0:ilon1,2,1:nf) = f(ilon0:ilon1,j,:)
       endif
    enddo

!
! Do the exchange:
!
    call mpi_allreduce( sndbuf(:,:,1:nf), fpoles(:,:,1:nf), len, MPI_REAL8, MPI_SUM, mpi_comm_edyn, ier )
    if ( ier .ne. MPI_SUCCESS ) call handle_mpi_err(ier,'mp_magpoles call mpi_allreduce')

  end subroutine mp_magpoles
!-----------------------------------------------------------------------
  integer function getpe(ix,jx)
    integer,intent(in) :: ix,jx
    integer :: it

    getpe = -1
    do it=0,ntask-1
      if ((tasks(it)%lon0 <= ix .and. tasks(it)%lon1 >= ix).and.&
          (tasks(it)%lat0 <= jx .and. tasks(it)%lat1 >= jx)) then
        getpe = it
        exit
      endif
    enddo
    if (getpe < 0) then
      write(iulog,"('getpe: pe with ix=',i4,' not found.')") ix
      call endrun('getpe')
    endif
  end function getpe
!-----------------------------------------------------------------------
  subroutine mp_pole_halos(f,lev0,lev1,lon0,lon1,lat0,lat1,nf,polesign)
!
! Set latitude halo points over the poles.
!
! Args:
    integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,nf
    real(r8),intent(in) :: polesign(nf)
    type(array_ptr_type) :: f(nf) ! (plev,i0-2:i1+2,j0-2:j1+2)
!
! Local:
    integer :: if,i,j,k,ihalo,it,i0,i1,j0,j1,itask

!   real(r8) :: fglblon(lev0:lev1,nlon,lat0-2:lat1+2,nf)
    type(array_ptr_type) :: pglblon(nf) ! (lev0:lev1,nlon,lat0-2:lat1+2)

    if (mytidj /= 0 .and. mytidj /= ntaskj-1) return

!   fglblon = 0._r8 ! init
!
! Allocate local fields with global longitudes:
    do if=1,nf
      allocate(pglblon(if)%ptr(lev0:lev1,nlon,lat0-2:lat1+2))
    enddo
!
! Define my subdomain in local fglblon, which has global lon dimension:
!
    do if=1,nf
      do j=lat0-2,lat1+2
        do i=lon0,lon1
          pglblon(if)%ptr(lev0:lev1,i,j) = f(if)%ptr(lev0:lev1,i,j)
        enddo
      enddo
    enddo
!
! Gather longitude data to westernmost processors (far north and south):
!
    call mp_gatherlons_f3d(pglblon,lev0,lev1,lon0,lon1,lat0-2,lat1+2,nf)
!
! Loop over tasks in my latitude row (far north or far south),
! including myself, and set halo points over the poles.
!
    if (mytidi==0) then
      do it=0,ntaski-1
        itask = tasks(itask_table_geo(it,mytidj))%mytid
        i0 = tasks(itask)%lon0
        i1 = tasks(itask)%lon1
        j0 = tasks(itask)%lat0
        j1 = tasks(itask)%lat1
        do if=1,nf
          if (j0==1) then ! south
            do i=i0,i1
              ihalo = 1+mod(i-1+nlon/2,nlon)
              pglblon(if)%ptr(lev0:lev1,i,j0-2) = pglblon(if)%ptr(lev0:lev1,ihalo,j0+2) ! get lat -1 from lat 3
              pglblon(if)%ptr(lev0:lev1,i,j0-1) = pglblon(if)%ptr(lev0:lev1,ihalo,j0+1) ! get lat 0 from lat 2
            enddo
          else            ! north
            do i=i0,i1
              ihalo = 1+mod(i-1+nlon/2,nlon)
              pglblon(if)%ptr(lev0:lev1,i,j1+1) = pglblon(if)%ptr(lev0:lev1,ihalo,j1-1) ! get lat plat+1 from plat-1
              pglblon(if)%ptr(lev0:lev1,i,j1+2) = pglblon(if)%ptr(lev0:lev1,ihalo,j1-2) ! get lat plat+2 from plat-2
            enddo
          endif
        enddo ! if=1,nf
      enddo ! it=0,ntaski-1
    endif ! mytidi==0
!
! Scatter data back out to processors in my latitude row:
!
    call mp_scatterlons_f3d(pglblon,lev0,lev1,lon0,lon1,lat0-2,lat1+2,nf)
!
! Finally, define halo points in data arrays from local global lon array,
! changing sign if necessary (winds):
!
    if (lat0==1) then ! south
      do if=1,nf
        do j=lat0-2,lat0-1
          do k=lev0,lev1
            f(if)%ptr(k,lon0:lon1,j) = pglblon(if)%ptr(k,lon0:lon1,j)*polesign(if)
          enddo
        enddo
      enddo
    else              ! north
      do if=1,nf
        do j=lat1+1,lat1+2
          do k=lev0,lev1
            f(if)%ptr(k,lon0:lon1,j) = pglblon(if)%ptr(k,lon0:lon1,j)*polesign(if)
          enddo
        enddo
      enddo
    endif

    do if=1,nf
      deallocate(pglblon(if)%ptr)
    enddo
  end subroutine mp_pole_halos
!-----------------------------------------------------------------------
  subroutine conjugate_points
  use edyn_maggrid,only: gmlat
!
! Local:
    integer :: ier,j,js,jn,itask,jj
!
! nsend_south(ntask): number of lats in south to send north
! nrecv_north(ntask): number of lats in north to recv from south
!
    allocate(nsend_south(0:ntask-1),stat=ier)
    allocate(nrecv_north(0:ntask-1),stat=ier)
!
! send_south_coords: south j lats to send north
! recv_north_coords: north j lats to recv from south
!
    allocate(send_south_coords(mxmaglat,0:ntask-1),stat=ier)
    allocate(recv_north_coords(mxmaglat,0:ntask-1),stat=ier)

    nsend_south(:) = 0
    nrecv_north(:) = 0
    send_south_coords(:,:) = 0
    recv_north_coords(:,:) = 0

    magloop: do j=mlat0,mlat1
!
! In north hem: find tasks w/ conjugate points in south to recv:
! (nmlath is in params module)
      if (gmlat(j) > 0._r8) then  ! in north hem of current task
        js = nmlath-(j-nmlath) ! j index to south conjugate point (should be -j)
        do itask=0,ntask-1
          do jj = tasks(itask)%mlat0,tasks(itask)%mlat1
!
! Receive these north coords from the south:
            if (jj==js.and.mlon0==tasks(itask)%mlon0.and. &
                           mlon1==tasks(itask)%mlon1) then
              nrecv_north(itask) = nrecv_north(itask)+1
              recv_north_coords(nrecv_north(itask),itask) = j
            endif
          enddo ! jj of remote task
        enddo ! itask=0,ntask-1
        if (all(nrecv_north==0)) &
          write(iulog,"(2a,i4,a,f8.2)") '>>> WARNING: could not find north conjugate',&
            ' points corresponding to south latitude js=',js,' gmlat(js)=',gmlat(js)
!
! In south hem: find tasks w/ conjugate points in north to send:
      elseif (gmlat(j) < 0._r8.and.j /= nmlath) then ! in south hem
        jn = nmlath+(nmlath-j) ! j index of north conjugate point
        do itask=0,ntask-1
          do jj = tasks(itask)%mlat0,tasks(itask)%mlat1
            if (jj==jn.and.mlon0==tasks(itask)%mlon0.and. &
                           mlon1==tasks(itask)%mlon1) then
              nsend_south(itask) = nsend_south(itask)+1
! Send these south coords to the north:
              send_south_coords(nsend_south(itask),itask) = j
            endif
          enddo ! jj of remote task
        enddo ! itask=0,ntask-1
        if (all(nsend_south==0)) &
          write(iulog,"(2a,i4,a,f8.2)") '>>> WARNING: could not find south conjugate',&
            ' points corresponding to north latitude jn=',jn,' gmlat(jn)=',gmlat(jn)
      endif ! in north or south hem
    enddo magloop ! j=mlat0,mlat1
  end subroutine conjugate_points
!-----------------------------------------------------------------------
  subroutine mp_mag_foldhem(f,mlon0,mlon1,mlat0,mlat1,nf)
!
! For each point in northern hemisphere (if any) of the current task
!   subdomain, receive data from conjugate point in the south (from the
!   south task that owns it), and sum it to the north point data.
!   Do this for nf fields. Conjugate point indices to send/recv to/from
!   each task were determined by sub conjugate_points (this module).
! nsend_south,       ! number of south lats to send to north (each task)
! nrecv_north        ! number of north lats to send to south (each task)
!
! This routine is called from edynamo at every timestep.
! Sub conjugate_points is called once per run, from mp_distribute.
!
! Args:
    integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nf
    real(r8),intent(inout) :: f(mlon0:mlon1,mlat0:mlat1,nf)
!
! Local:
    integer :: j,n,len,itask,ifld,ier,nmlons
    real(r8) :: sndbuf(mxmaglon,mxmaglat,nf,0:ntask-1)
    real(r8) :: rcvbuf(mxmaglon,mxmaglat,nf,0:ntask-1)
    integer :: jsend(0:ntask-1),jrecv(0:ntask-1)
    integer :: irstat(MPI_STATUS_SIZE)      ! mpi receive status

!
    sndbuf = 0._r8 ; rcvbuf = 0._r8
    jsend  = 0     ; jrecv  = 0
    len = mxmaglon*mxmaglat*nf
    nmlons = mlon1-mlon0+1
!
! Send south data to north itask:
! (To avoid deadlock, do not send if north task is also myself. This will
!  happen when there is an odd number of tasks in the latitude dimension,
!  e.g., ntask == 12, 30, etc)
!
    do itask=0,ntask-1

! Attempt to fetch from allocatable variable NSEND_SOUTH when it is not allocated

      if (nsend_south(itask) > 0 .and. itask /= mytid) then
        do ifld = 1,nf
          do n=1,nsend_south(itask)
            sndbuf(1:nmlons,n,ifld,itask) = &
              f(:,send_south_coords(n,itask),ifld)
          enddo
        enddo ! ifld=1,nf
        call mpi_isend(sndbuf(1,1,1,itask),len,MPI_REAL8, &
          itask,1,mpi_comm_edyn,jsend(itask),ier)
        call mpi_wait(jsend(itask),irstat,ier)
      endif ! nsend_south(itask) > 0
    enddo ! itask=0,ntask-1
!
! Receive north data from south itask and add to north,
! i.e., north = north+south. (do not receive if south task is
! also myself, but do add south data to my north points, see below)
!
    do itask=0,ntask-1
      if (nrecv_north(itask) > 0 .and. itask /= mytid) then
        call mpi_irecv(rcvbuf(1,1,1,itask),len,MPI_REAL8, &
          itask,1,mpi_comm_edyn,jrecv(itask),ier)
        call mpi_wait(jrecv(itask),irstat,ier)
        do ifld=1,nf
          do n=1,nrecv_north(itask)
!
! Receive lats in reverse order:
            f(mlon0:mlon1, &
              recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) = &
            f(mlon0:mlon1, &
              recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) + &
            rcvbuf(1:nmlons,n,ifld,itask)
          enddo ! n=1,nrecv_north(itask)
        enddo ! ifld=1,nf
!
! If I am send *and* receive task, simply add my south data to my north points:
      elseif (nrecv_north(itask) > 0 .and. itask == mytid) then
        do ifld=1,nf
          do n=1,nrecv_north(itask)
            f(mlon0:mlon1, &
              recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) = &
            f(mlon0:mlon1, &
              recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) + &
            f(mlon0:mlon1,send_south_coords(n,itask),ifld)
          enddo ! n=1,nrecv_north(itask)
        enddo ! ifld=1,nf
      endif ! nrecv_north(itask) > 0
    enddo ! itask=0,ntask-1
!
! Mag equator is also "folded", but not included in conjugate points,
! so double it here:
    do j=mlat0,mlat1
      if (j==nmlath) then
        do ifld=1,nf
          f(:,j,ifld) = f(:,j,ifld)+f(:,j,ifld)
        enddo
      endif
    enddo

  end subroutine mp_mag_foldhem
!-----------------------------------------------------------------------
  subroutine mp_mag_periodic_f2d(f,mlon0,mlon1,mlat0,mlat1,nf)
!
! Args:
    integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nf
    real(r8),intent(inout) :: f(mlon0:mlon1,mlat0:mlat1,nf)
!
! Local:
    integer :: j,ier,idest,isrc,len,ireqsend,ireqrecv,msgtag
    real(r8) :: sndbuf(mxmaglat,nf),rcvbuf(mxmaglat,nf)
    integer  :: irstat(MPI_STATUS_SIZE)      ! mpi receive status

    if (ntaski>1) then
       len = mxmaglat*nf
       !
       ! I am a western-most task. Send lon 1 to eastern-most tasks:
       if (mytidi==0) then
          idest = itask_table_mag(ntaski-1,mytidj)
          do j=mlat0,mlat1
             sndbuf(j-mlat0+1,:) = f(1,j,:)
          enddo
          msgtag = mytid
          call mpi_isend(sndbuf,len,MPI_REAL8,idest,msgtag,mpi_comm_edyn, ireqsend,ier)
          if (ier /= 0) call handle_mpi_err(ier,'mp_mag_periodic_f2d send to idest')
          call mpi_wait(ireqsend,irstat,ier)
          if (ier /= 0) call handle_mpi_err(ier,'mp_mag_periodic_f2d wait for send')
          !
          ! I am eastern-most task. Receive lon 1 from western-most tasks,
          ! and assign to nmlonp1:
       elseif (mytidi==ntaski-1) then
          isrc = itask_table_mag(0,mytidj)
          msgtag = isrc
          call mpi_irecv(rcvbuf,len,MPI_REAL8,isrc,msgtag,mpi_comm_edyn, ireqrecv,ier)
          if (ier /= 0) call handle_mpi_err(ier,'mp_mag_periodic_f2d recv from isrc')
          call mpi_wait(ireqrecv,irstat,ier)
          if (ier /= 0) call handle_mpi_err(ier,'mp_mag_periodic_f2d wait for recv')

          do j=mlat0,mlat1
             f(nmlonp1,j,:) = rcvbuf(j-mlat0+1,:)
          enddo
       endif ! mytidi == 0 or ntaski-1
    else
       do j=mlat0,mlat1
          f(nmlonp1,j,:) = f(1,j,:)
       enddo
    endif

  end subroutine mp_mag_periodic_f2d
!-----------------------------------------------------------------------
  subroutine mp_mag_halos(fmsub,mlon0,mlon1,mlat0,mlat1,nf)
!
! Exchange halo/ghost points between magnetic grid subdomains for nf fields.
! Only a single halo point is required in both lon and lat dimensions.
! Note that all tasks in any row of the task matrix have the same
!   mlat0,mlat1, and that all tasks in any column of the task matrix
!   have the same mlon0,mlon1.
! Longitude halos are done first, exchanging mlat0:mlat1, then latitude
!   halos are done, exchanging mlon0-1:mlon1+1 (i.e., including the
!   longitude halos that were defined first).
!
! Args:
    integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nf
    real(r8),intent(inout) :: fmsub(mlon0-1:mlon1+1,mlat0-1:mlat1+1,nf)
!
! Local:
    integer :: ifld,west,east,north,south,len,isend0,isend1, &
      irecv0,irecv1,ier,nmlats,istat(MPI_STATUS_SIZE,4),ireq(4),nmlons
    real(r8),dimension(mlat1-mlat0+1,nf)::sndlon0,sndlon1,rcvlon0,rcvlon1
    real(r8),dimension((mlon1+1)-(mlon0-1)+1,nf) :: &
      sndlat0,sndlat1,rcvlat0,rcvlat1

!
! Init send/recv buffers for lon halos:
    sndlon0 = 0._r8 ; rcvlon0 = 0._r8
    sndlon1 = 0._r8 ; rcvlon1 = 0._r8
!
! Identify east and west neightbors:
    west  = itask_table_mag(mytidi-1,mytidj)
    east  = itask_table_mag(mytidi+1,mytidj)
!
! Exchange mlat0:mlat1 (lat halos are not yet defined):
    nmlats = mlat1-mlat0+1
    len = nmlats*nf
!
! Send mlon0 to the west neighbor, and mlon1 to the east.
! However, tasks are periodic in longitude (see itask_table_mag),
!   and far west tasks send mlon0+1, and far east tasks send mlon1-1
!
    do ifld=1,nf
! Far west tasks send mlon0+1 to far east (periodic) tasks:
      if (mytidi==0) then
        sndlon0(:,ifld) = fmsub(mlon0+1,mlat0:mlat1,ifld)
! Interior tasks send mlon0 to west neighbor:
      else
        sndlon0(:,ifld) = fmsub(mlon0,mlat0:mlat1,ifld)
      endif

! Far east tasks send mlon1-1 to far west (periodic) tasks:
      if (mytidi==nmagtaski-1) then
        sndlon1(:,ifld) = fmsub(mlon1-1,mlat0:mlat1,ifld)
! Interior tasks send mlon1 to east neighbor:
      else
        sndlon1(:,ifld) = fmsub(mlon1,mlat0:mlat1,ifld)
      endif
    enddo ! ifld=1,nf
!
! Send mlon0 to the west:
    call mpi_isend(sndlon0,len,MPI_REAL8,west,1,mpi_comm_edyn,isend0,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos send mlon0 to west')
!
! Send mlon1 to the east:
    call mpi_isend(sndlon1,len,MPI_REAL8,east,1,mpi_comm_edyn,isend1,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos send mlon1 to east')
!
! Recv mlon0-1 from west:
    call mpi_irecv(rcvlon0,len,MPI_REAL8,west,1,mpi_comm_edyn,irecv0,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos recv mlon0 from west')
!
! Recv mlon1+1 from east:
    call mpi_irecv(rcvlon1,len,MPI_REAL8,east,1,mpi_comm_edyn,irecv1,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos recv mlon1 from east')
!
! Wait for completions:
    ireq = (/isend0,isend1,irecv0,irecv1/)
    istat = 0
    call mpi_waitall(4,ireq,istat,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos waitall for lons')
!
! Copy mlon0-1 from rcvlon0, and mlon1+1 from rcvlon1:
    do ifld=1,nf
      fmsub(mlon0-1,mlat0:mlat1,ifld) = rcvlon0(:,ifld)
      fmsub(mlon1+1,mlat0:mlat1,ifld) = rcvlon1(:,ifld)
!
! Fix special case of 2 tasks in longitude dimension:
      if (east == west) then
        fmsub(mlon0-1,mlat0:mlat1,ifld) = rcvlon1(:,ifld)
        fmsub(mlon1+1,mlat0:mlat1,ifld) = rcvlon0(:,ifld)
      endif
    enddo ! ifld=1,nf
!
! Now exchange latitudes:
    sndlat0 = 0._r8 ; rcvlat0 = 0._r8
    sndlat1 = 0._r8 ; rcvlat1 = 0._r8

    south = itask_table_mag(mytidi,mytidj-1)  ! neighbor to south
    north = itask_table_mag(mytidi,mytidj+1)  ! neighbor to north
!
! Include halo longitudes that were defined by the exchanges above:
    nmlons = (mlon1+1)-(mlon0-1)+1
    len = nmlons*nf
!
! Send mlat0 to south neighbor, and mlat1 to north:
    do ifld=1,nf
      sndlat0(:,ifld) = fmsub(:,mlat0,ifld)
      sndlat1(:,ifld) = fmsub(:,mlat1,ifld)
    enddo
!
! Send mlat0 to south:
    call mpi_isend(sndlat0,len,MPI_REAL8,south,1,mpi_comm_edyn,isend0,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos send mlat0 to south')
!
! Send mlat1 to north:
    call mpi_isend(sndlat1,len,MPI_REAL8,north,1,mpi_comm_edyn,isend1,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos send mlat1 to north')
!
! Recv mlat0-1 from south:
    call mpi_irecv(rcvlat0,len,MPI_REAL8,south,1,mpi_comm_edyn,irecv0,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos recv mlat0-1 from south')
!
! Recv mlat1+1 from north:
    call mpi_irecv(rcvlat1,len,MPI_REAL8,north,1,mpi_comm_edyn,irecv1,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos recv mlat1+1 from north')
!
! Wait for completions:
    ireq = (/isend0,isend1,irecv0,irecv1/)
    istat = 0
    call mpi_waitall(4,ireq,istat,ier)
    if (ier /= 0) call handle_mpi_err(ier,'mp_mag_halos waitall for lats')
!
! Copy mlat0-1 from rcvlat0, and mlat1+1 from rcvlat1:
    do ifld=1,nf
      fmsub(:,mlat0-1,ifld) = rcvlat0(:,ifld)
      fmsub(:,mlat1+1,ifld) = rcvlat1(:,ifld)
    enddo ! ifld=1,nf

  end subroutine mp_mag_halos
!-----------------------------------------------------------------------
  subroutine mp_geo_halos(fmsub,lev0,lev1,lon0,lon1,lat0,lat1,nf)
!
! Exchange halo/ghost points between geographic grid subdomains for nf fields.
! Two halo points are set in both lon and lat dimensions.
! Longitude halos are done first, then latitude halos are done, including
!   longitude halos that were defined first).
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,nf
      type(array_ptr_type) :: fmsub(nf) ! (lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)

!
! Local:
      integer :: k,i,ifld,west,east,north,south,len,isend0,isend1, &
        irecv0,irecv1,ier,nlats,istat(MPI_STATUS_SIZE,4),ireq(4),nlons
      real(r8),dimension(lev0:lev1,2,lat1-lat0+1,nf) :: &
        sndlon0,sndlon1,rcvlon0,rcvlon1
      real(r8),dimension(lev0:lev1,2,(lon1+2)-(lon0-2)+1,nf) :: &
        sndlat0,sndlat1,rcvlat0,rcvlat1

!     if (mpi_timing) starttime = mpi_wtime()
!
! Init send/recv buffers for lon halos:
      sndlon0 = 0._r8 ; rcvlon0 = 0._r8
      sndlon1 = 0._r8 ; rcvlon1 = 0._r8
!
! Identify east and west neighbors:
      west  = itask_table_geo(mytidi-1,mytidj)
      east  = itask_table_geo(mytidi+1,mytidj)
!
! Exchange lat0:lat1 (lat halos are not yet defined):
      nlats = lat1-lat0+1
      len = (lev1-lev0+1)*2*nlats*nf
!
! Send lon0:lon0+1 to the west neighbor, and lon1-1:lon1 to the east.
!
      do ifld=1,nf
        do i=1,2
          do k=lev0,lev1
            sndlon0(k,i,:,ifld) = fmsub(ifld)%ptr(k,lon0+i-1,lat0:lat1) ! lon0, lon0+1
            sndlon1(k,i,:,ifld) = fmsub(ifld)%ptr(k,lon1+i-2,lat0:lat1) ! lon1-1, lon1
          enddo
        enddo
      enddo ! ifld=1,nf
!
! Send lon0:lon0+1 to the west:
      call mpi_isend(sndlon0,len,MPI_REAL8,west,1,mpi_comm_edyn,isend0,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
        'mp_geo_halos send lon0:lon0+1 to west')
!
! Send lon1-1:lon1 to the east:
      call mpi_isend(sndlon1,len,MPI_REAL8,east,1,mpi_comm_edyn,isend1,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
        'mp_geo_halos send lon1-1:lon1 to east')
!
! Recv lon0-2:lon0-1 from west:
      call mpi_irecv(rcvlon0,len,MPI_REAL8,west,1,mpi_comm_edyn,irecv0,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
        'mp_geo_halos recv lon0-2:lon0-1 from west')
!
! Recv lon1+1:lon1+2 from east:
      call mpi_irecv(rcvlon1,len,MPI_REAL8,east,1,mpi_comm_edyn,irecv1,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
          'mp_geo_halos recv lon1+1:lon1+2 from east')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
        'mp_geo_halos waitall for lons')
!
! Copy lon0-2:lon0-1 from rcvlon0, and lon1+1:lon1+2 from rcvlon1:
      do ifld=1,nf
        if (east /= west) then
          do i=1,2
            do k=lev0,lev1
              fmsub(ifld)%ptr(k,lon0-3+i,lat0:lat1) = rcvlon0(k,i,:,ifld) ! lon0-2, lon0-1 
              fmsub(ifld)%ptr(k,lon1+i  ,lat0:lat1) = rcvlon1(k,i,:,ifld) ! lon1+1, lon1+2 
            enddo
          enddo ! i=1,2
!
! Fix special case of 2 tasks in longitude dimension:
        else ! east==west
          do i=1,2
            do k=lev0,lev1
              fmsub(ifld)%ptr(k,lon0-3+i,lat0:lat1) = rcvlon1(k,i,:,ifld) ! lon0-2, lon0-1 
              fmsub(ifld)%ptr(k,lon1+i  ,lat0:lat1) = rcvlon0(k,i,:,ifld) ! lon1+1, lon1+2 
            enddo
          enddo
        endif ! east==west
      enddo ! ifld=1,nf
!
! Now exchange latitudes:
      sndlat0 = 0._r8 ; rcvlat0 = 0._r8
      sndlat1 = 0._r8 ; rcvlat1 = 0._r8

      south = itask_table_geo(mytidi,mytidj-1)  ! neighbor to south
      north = itask_table_geo(mytidi,mytidj+1)  ! neighbor to north
!
! Include halo longitudes that were defined by the exchanges above:
      nlons = (lon1+2)-(lon0-2)+1 
      len = (lev1-lev0+1)*2*nlons*nf
!
! Send lat0:lat0+1 to south neighbor, and lat1-1:lat1 to north:
      do ifld=1,nf
        do k=lev0,lev1
          sndlat0(k,1,:,ifld) = fmsub(ifld)%ptr(k,:,lat0  ) ! send lat0   to south
          sndlat0(k,2,:,ifld) = fmsub(ifld)%ptr(k,:,lat0+1) ! send lat0+1 to south

          sndlat1(k,1,:,ifld) = fmsub(ifld)%ptr(k,:,lat1  ) ! send lat1   to north 
          sndlat1(k,2,:,ifld) = fmsub(ifld)%ptr(k,:,lat1-1) ! send lat1-1 to north 
        enddo
      enddo
!
! Send lat0:lat0+1 to south (matching recv is lat1+1:lat1+2):
      call mpi_isend(sndlat0,len,MPI_REAL8,south,100,mpi_comm_edyn,isend0,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
          'mp_geo_halos send lat0:lat0+1 to south')
!
! Send lat1-1:lat1 to north (matching recv is lat0-2:lat0-1):
      call mpi_isend(sndlat1,len,MPI_REAL8,north,101,mpi_comm_edyn,isend1,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
          'mp_geo_halos send lat1-1:lat1 to north')
!
! Recv lat0-2:lat0-1 from south:
      call mpi_irecv(rcvlat0,len,MPI_REAL8,south,101,mpi_comm_edyn,irecv0,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
          'mp_geo_halos recv lat0-2:lat0-1 from south')
!
! Recv lat1+1:lat1+2 from north:
      call mpi_irecv(rcvlat1,len,MPI_REAL8,north,100,mpi_comm_edyn,irecv1,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
          'mp_geo_halos recv lat1+1:lat1+2 from north')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) call handle_mpi_err(ier, &
        'mp_geo_halos waitall for lats')
!
! Copy lat0-2:lat0-1 from rcvlat0, and lat1+1:lat1+2 from rcvlat1:
      do ifld=1,nf
        do k=lev0,lev1
          fmsub(ifld)%ptr(k,:,lat0-1) = rcvlat0(k,1,:,ifld) ! recv lat0-1 from south
          fmsub(ifld)%ptr(k,:,lat0-2) = rcvlat0(k,2,:,ifld) ! recv lat0-2 from south

          fmsub(ifld)%ptr(k,:,lat1+1) = rcvlat1(k,1,:,ifld) ! recv lat1+1 from north
          fmsub(ifld)%ptr(k,:,lat1+2) = rcvlat1(k,2,:,ifld) ! recv lat1+2 from north
        enddo
!
! Fix special case of 2 tasks in latitude dimension:
! Not sure if this will happen in WACCM:
!
        if (north == south) then
          call endrun('mp_geo_halos: north==south')
        endif
      enddo ! ifld=1,nf

      end subroutine mp_geo_halos
!-----------------------------------------------------------------------
  subroutine mp_gather_edyn(fmsub,mlon0,mlon1,mlat0,mlat1,fmglb,nmlonp1,nmlat,nf)
!
! Gather fields on mag subdomains to root task, so root task can
! complete non-parallel portion of dynamo (starting after rhspde)
!
! Args:
    integer,intent(in)   :: mlon0,mlon1,mlat0,mlat1,nmlonp1,nmlat,nf
    real(r8),intent(in)  :: fmsub(mlon0:mlon1,mlat0:mlat1,nf)
    real(r8),intent(out) :: fmglb(nmlonp1,nmlat,nf)
!
! Local:
    integer :: len,i,j,ifld,ier
    real(r8),dimension(nmlonp1,nmlat,nf) :: sndbuf

    sndbuf = 0._r8
    fmglb = 0._r8

    len = nmlonp1*nmlat*nf
!
! Load send buffer with my subdomain:
    do ifld=1,nf
       do j=mlat0,mlat1
          do i=mlon0, mlon1
             sndbuf(i,j,ifld) = fmsub(i,j,ifld)
          enddo
       enddo
    enddo

!
! Gather to root by using scalable reduce method:

    call mpi_reduce(sndbuf, fmglb, len, MPI_REAL8, MPI_SUM, 0, mpi_comm_edyn, ier )
    if (ier /= 0) call handle_mpi_err(ier,'mp_gather_edyn: mpi_gather to root')

  end subroutine mp_gather_edyn
!-----------------------------------------------------------------------
  subroutine mp_scatter_phim(phim_glb,phim)
    real(r8),intent(in)  :: phim_glb(nmlonp1,nmlat)
    real(r8),intent(out) :: phim(mlon0:mlon1,mlat0:mlat1)
!
! Local:
    integer :: ier,len,i,j

!   if (mpi_timing) starttime = mpi_wtime()
!
! Broadcast global phim (from pdynamo phim(nmlonp1,nmlat)):
    len = nmlat*nmlonp1
    call mpi_bcast(phim_glb,len,MPI_REAL8,0,mpi_comm_edyn,ier)
    if (ier /= 0) &
      call handle_mpi_err(ier,'mp_scatter_phim: bcast global phim')
!
! Define subdomains:
    do j=mlat0,mlat1
      do i=mlon0,mlon1
        phim(i,j) = phim_glb(i,j)
      enddo
    enddo

  end subroutine mp_scatter_phim
!-----------------------------------------------------------------------
  subroutine mp_mag_jslot(fin,mlon00,mlon11,mlat00,mlat11, &
    fout,jneed,mxneed,nf)
!
! Current task needs to receive (from other tasks) field f at (non-zero)
! latitude indices in jneed, at all longitudes in the current subdomain.
! Note subdomains include halo points mlon0-1 and mlat1+1. Data in f also
! includes halo points (will need the lat data at halo-longitudes)
!
! Args:
    integer,intent(in) :: mlon00,mlon11,mlat00,mlat11 ! subdomains w/ halos
    integer,intent(in) :: nf ! number of fields
    integer,intent(in) ::  mxneed ! max number of needed lats (nmlat+2)
    integer,intent(in) :: jneed(mxneed) ! j-indices of needed lats (where /= -1)
    real(r8),intent(in) :: fin(mlon00:mlon11,mlat00:mlat11,nf) ! data at current subdomain
    real(r8),intent(out) :: fout(mlon00:mlon11,mxneed,nf) ! returned data at needed lats
    !
    ! Local:
    integer,parameter :: sndbuf_cntr_max = 20 ! Maximum number of ibsend from one mpi task
    integer :: ier,njneed,i,j,n,nj,idest, &
      icount,len,nlons,isrc,msgid,ifld,sndbuf_cntr
    integer :: tij ! rank in cols_comm (0 to nmagtaskj-1)
    integer :: jhave(mxneed),njhave,wid
    integer :: peersneed(mxneed,0:nmagtaskj-1)
    integer :: jneedall (mxneed,0:nmagtaskj-1)
    real(r8) :: sndbuf(mxmaglon+2,mxneed,nf,sndbuf_cntr_max)
    real(r8) :: rcvbuf(mxmaglon+2,mxneed,nf)
    real(r8) :: buffer((mxmaglon+2)*mxneed*nf*sndbuf_cntr_max)
    integer :: irstat(MPI_STATUS_SIZE)          ! mpi receive status
    integer :: isstat(MPI_STATUS_SIZE,sndbuf_cntr_max) !mpi_ibsend wait status
    integer :: ibsend_requests(sndbuf_cntr_max) !array of ibsend requests

    sndbuf = 0._r8
    rcvbuf = 0._r8
    njneed = 0
    ibsend_requests = 0
    sndbuf_cntr = 0
    do j=1,mxneed
      if (jneed(j) /= -1) njneed=njneed+1
    enddo
    if (any(jneed(1:njneed)==-1)) call endrun('mp_mag_jslot jneed')
    !
    call MPI_Comm_rank(cols_comm,tij,ier)
    call MPI_buffer_attach(buffer,(mxmaglon+2)*mxneed*nf*sndbuf_cntr_max,ier)
    if (ier /= 0) &
      call handle_mpi_err(ier,'mp_mag_jslot call mpi_buffer_attach')

    !
    ! Send needed lat indices to all tasks in my column:
    ! (redundant for alltoall)
    do n=0,nmagtaskj-1
      jneedall(:,n) = jneed(:)
    enddo

    call mpi_alltoall(jneedall,mxneed,MPI_INTEGER, &
      peersneed,mxneed,MPI_INTEGER,cols_comm,ier)
    if (ier /= 0) &
      call handle_mpi_err(ier,'mp_mag_jslot call mpi_alltoall')
    !
    ! Check if I have any needed lats, and who to send to:
    do n=0,nmagtaskj-1
      if (n==tij) cycle
      njhave = 0
      do j=1,mxneed
        if (peersneed(j,n) >= mlat00.and.peersneed(j,n) <= mlat11)then
          njhave = njhave+1
          jhave(njhave) = peersneed(j,n)
          idest = n
          wid = itask_table_geo(mytidi,idest)
        endif
      enddo
      if (njhave > 0) then

        sndbuf_cntr = sndbuf_cntr + 1
        if (sndbuf_cntr > sndbuf_cntr_max) call endrun('sndbuf_cntr exceeded sndbuf_cntr_max')

        !
        ! Load send buffer:
        nlons = mlon11-mlon00+1
        do ifld=1,nf
          do j=1,njhave
            do i=mlon00,mlon11
              sndbuf(i-mlon00+1,j,ifld,sndbuf_cntr) = fin(i,jhave(j),ifld)
            enddo
          enddo
        enddo
        len = nlons*njhave*nf
        msgid = mytid+wid*10000
        call mpi_ibsend(sndbuf(1:nlons,1:njhave,:,sndbuf_cntr),len,MPI_REAL8, &
          idest,msgid,cols_comm,ibsend_requests(sndbuf_cntr),ier)
        if (ier /= 0) &
          call handle_mpi_err(ier,'mp_mag_jslot call mpi_ibsend')
      endif
    enddo ! n=0,nmagtaskj-1

    call MPI_waitall(sndbuf_cntr,ibsend_requests,isstat,ier)
    if (ier /= 0) &
      call handle_mpi_err(ier,'mp_mag_jslot call mpi_waitall')
    call MPI_buffer_detach(buffer,(mxmaglon+2)*mxneed*nf*sndbuf_cntr_max,ier)
    if (ier /= 0) &
      call handle_mpi_err(ier,'mp_mag_jslot call mpi_buffer_detach')

    !
    ! Determine which tasks to receive which lats from. Task to
    ! receive from must be in same task column magtidi as I am.
    if (njneed > 0) then
      njhave = 0
      jhave(:) = -1
      do n=0,ntask-1
        njhave = 0
        do j=1,njneed
          if (jneed(j) >= tasks(n)%mlat0-1 .and. &
              jneed(j) <= tasks(n)%mlat1+1) then
            njhave = njhave+1
            jhave(njhave) = jneed(j)
          endif
        enddo
        if (njhave > 0 .and. tasks(n)%magtidi==magtidi) then
          isrc = tasks(n)%magtidj ! task id in cols_comm to recv from
          nlons = mlon11-mlon00+1
          len = nlons*njhave*nf
          msgid = mytid*10000+n
          rcvbuf = 0._r8
          call mpi_recv(rcvbuf(1:nlons,1:njhave,:),len,MPI_REAL8, &
            isrc,msgid,cols_comm,irstat,ier)
          if (ier /= 0) &
            call handle_mpi_err(ier,'mp_mag_jslot call mpi_recv')
          !
          ! Get data from receive buffer:
          ! real,intent(out) :: fout(mlon00:mlon11,mxneed) ! returned data at needed lats
          do ifld=1,nf
            do j=1,njhave
              nj = ixfind(jneed,mxneed,jhave(j),icount)
              if (nj==0) call endrun('jhave(j) not in jneed')
              do i=mlon00,mlon11
                fout(i,nj,ifld) = rcvbuf(i-mlon00+1,j,ifld)
              enddo
            enddo ! j=1,njhave
          enddo ! ifld=1,nf
        endif ! jhave > 0
      enddo ! n=0,ntask-1
    endif ! njneed > 0

  end subroutine mp_mag_jslot
!-----------------------------------------------------------------------
  subroutine mp_gatherlons_f3d(f,k0,k1,i0,i1,j0,j1,nflds)
!
! Gather longitude data in a row of tasks to leftmost task in the row.
! On entry f(k0:k1,i0:i1,j0:j1,nflds) is defined for current task.
! On exit f(k0:k1,nlonp4,j0:j1,nflds) is defined for task with mytidi==0.
!

!
! Args:
!
    integer,intent(in) :: k0,k1,i0,i1,j0,j1,nflds
!   real(r8),intent(inout) :: f(k0:k1,nlon,j0:j1,nflds)
    type(array_ptr_type) :: f(nflds) ! f(n)%ptr(k0:k1,nlon,j0:j1)
!
! Local:
!
    integer :: irstat(MPI_STATUS_SIZE)      ! mpi receive status
    integer :: j,n,nlons,nlonrecv,nlevs,len,idest,isrc,ier, &
                                        isend,irecv,itask,lonrecv0,lonrecv1,mtag
    real(r8) :: &
      sndbuf(k0:k1,mxlon,mxlat+4,nflds), & ! send buffer
      rcvbuf(k0:k1,mxlon,mxlat+4,nflds)    ! recv buffer
!
! Exec:
!
    nlons = i1-i0+1
    nlevs = k1-k0+1

    sndbuf = 0._r8 
    rcvbuf = 0._r8
    len = nlevs*mxlon*(mxlat+4)*nflds ! +4 is for when this is called from mp_pole_halos
!
! If mytidi==0, receive from other tasks in my row (mytidi>0,mytidj):
    if (mytidi == 0) then
      do itask=1,ntaski-1
        isrc = itask_table_geo(itask,mytidj)
        mtag = isrc+mytid
        call mpi_irecv(rcvbuf,len,MPI_REAL8,isrc,mtag,mpi_comm_edyn,irecv,ier)
        if (ier /= 0) &
          call handle_mpi_err(ier,'mp_gatherlons_f3d recv fm isrc')
        call mpi_wait(irecv,irstat,ier)
        if (ier /= 0) &
           call handle_mpi_err(ier,'mp_gatherlons_f3d wait for recv0')
!
! Copy data from receive buffer:
        lonrecv0 = tasks(isrc)%lon0
        lonrecv1 = tasks(isrc)%lon1
        nlonrecv = lonrecv1-lonrecv0+1
        do n=1,nflds
          do j=j0,j1
            f(n)%ptr(k0:k1,lonrecv0:lonrecv1,j) = rcvbuf(k0:k1,1:nlonrecv,j-j0+1,n)
          enddo ! j=j0,j1
        enddo ! n=1,nflds
      enddo ! itask=1,ntaski-1
!
! If mytidi > 0, load send buffer, and send to task (0,mytidj):
    else ! mytidi /= 0
      idest = itask_table_geo(0,mytidj)
      do n=1,nflds
        do j=j0,j1
          sndbuf(:,1:nlons,j-j0+1,n) = f(n)%ptr(k0:k1,i0:i1,j)
        enddo ! j=j0,j1
      enddo ! n=1,nflds
      mtag = idest+mytid
      call mpi_isend(sndbuf,len,MPI_REAL8,idest,mtag,mpi_comm_edyn,isend,ier)
      if (ier /= 0) &
        call handle_mpi_err(ier,'mp_gatherlons_f3d send0 to idest')
      call mpi_wait(isend,irstat,ier)
      if (ier /= 0) &
        call handle_mpi_err(ier,'mp_gatherlons_f3d wait for send0')
    endif ! mytidi==0
  end subroutine mp_gatherlons_f3d
!-----------------------------------------------------------------------
  subroutine mp_scatterlons_f3d(f,k0,k1,i0,i1,j0,j1,nflds)
!
! Redistribute longitudes from left most task in j-row to other tasks
!   in the row.
! On input, f(:,nlonp4,j0:j1,nflds) is defined for tasks with mytidi==0.
! On output, f(:,i0:i1,j0:j1,nflds) is defined for all tasks.
!
! Args:
!
    integer,intent(in) :: k0,k1,i0,i1,j0,j1,nflds
    type(array_ptr_type) :: f(nflds) ! f(n)%ptr(k0:k1,nlon,j0:j1)
!
! Local:
!
    integer :: irstat(MPI_STATUS_SIZE)      ! mpi receive status
    integer :: j,n,nlevs,nlons,nlonsend,len,idest,isrc,ier, &
               isend,irecv,itask,lonsend0,lonsend1,mtag
    real(r8) :: &
      sndbuf(k0:k1,mxlon,mxlat+4,nflds), & ! send buffer
      rcvbuf(k0:k1,mxlon,mxlat+4,nflds)    ! recv buffer
!
! Exec:
!
    nlons = i1-i0+1
    nlevs = k1-k0+1

    sndbuf = 0._r8 ; rcvbuf = 0._r8
    len = nlevs*mxlon*(mxlat+4)*nflds ! +4 is for when this is called from mp_pole_halos
!
! If mytidi==0, send to other tasks in my row (mytidi>0,mytidj):
    if (mytidi == 0) then
      do itask=1,ntaski-1
        idest = itask_table_geo(itask,mytidj)
        lonsend0 = tasks(idest)%lon0
        lonsend1 = tasks(idest)%lon1
        nlonsend = lonsend1-lonsend0+1
        mtag = idest+mytid
        do n=1,nflds
          do j=j0,j1
            sndbuf(:,1:nlonsend,j-j0+1,n) = f(n)%ptr(:,lonsend0:lonsend1,j) 
          enddo ! j=j0,j1
        enddo ! n=1,nflds
        mtag = idest+mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,mtag,mpi_comm_edyn,isend,ier)
        if (ier /= 0) call handle_mpi_err(ier,'mp_scatterlons_f3d send to idest')
        call mpi_wait(isend,irstat,ier)
        if (ier /= 0) call handle_mpi_err(ier,'mp_scatterlons_f3d wait for send')
      enddo ! itask=1,ntaski-1
!
! If mytidi > 0, receive from task (0,mytidj):
    else
      isrc = itask_table_geo(0,mytidj)
      mtag = isrc+mytid
      call mpi_irecv(rcvbuf,len,MPI_REAL8,isrc,mtag,mpi_comm_edyn,irecv,ier)
      if (ier /= 0) &
         call handle_mpi_err(ier,'mp_scatterlons_f3d recv fm isrc')
      call mpi_wait(irecv,irstat,ier)
      if (ier /= 0) & 
         call handle_mpi_err(ier,'mp_scatterlons_f3d wait for recv')
      do n=1,nflds
        do j=j0,j1 
          f(n)%ptr(:,i0:i1,j) = rcvbuf(:,1:nlons,j-j0+1,n)
        enddo ! j=j0,j1
      enddo ! n=1,nflds
    endif
  end subroutine mp_scatterlons_f3d
!-----------------------------------------------------------------------
  subroutine handle_mpi_err(ierrcode,string)
!
! Args:
    integer,intent(in) :: ierrcode
    character(len=*) :: string
!
! Local:
    character(len=80) :: errstring
    integer :: len_errstring, ierr
! 
    call mpi_error_string(ierrcode,errstring,len_errstring, ierr)
    write(iulog,"(/,'>>> mpi error: ',a)") trim(string)
    write(iulog,"('  ierrcode=',i3,': ',a)") trim(errstring)
  end subroutine handle_mpi_err
!-----------------------------------------------------------------------
  integer function ixfind(iarray,idim,itarget,icount)
!
! Search iarray(idim) for itarget, returning first index in iarray
! where iarray(idim)==target. Also return number of elements of
! iarray that == itarget in icount.
!
! Args:
    integer,intent(in) :: idim,itarget
    integer,intent(in) :: iarray(idim)
    integer,intent(out) :: icount
!
! Local:
    integer :: i
!
    ixfind = 0
    icount = 0
    if (.not.any(iarray==itarget)) return
    icount = count(iarray==itarget)
    do i=1,idim
      if (iarray(i)==itarget) then
        ixfind = i
        exit
      endif
    enddo
  end function ixfind

!-----------------------------------------------------------------------
  subroutine setpoles(f,k0,k1,i0,i1,j0,j1)
!
! Args:
    integer,intent(in) :: k0,k1,i0,i1,j0,j1
    real(r8),intent(inout) :: f(k0:k1,i0:i1,j0:j1)
!
! Local:
    integer :: i,j,k,lon0,lon1,it,itask
    type(array_ptr_type) :: ptr(1)
    real(r8) :: fave(k0:k1)
    real(r8) :: rnlon

    if (j0 /= 1 .and. j1 /= nlat) return ! subdomain does not include poles

    rnlon = dble(nlon)
    allocate(ptr(1)%ptr(k0:k1,nlon,j0:j1))
!
! Define subdomains in global longitude dimension of ptmp:
!
    do j=j0,j1
      do i=i0,i1
        ptr(1)%ptr(k0:k1,i,j) = f(k0:k1,i,j)
      enddo
    enddo
!
! Get values for global longitudes at the latitude below each pole,
! average them at each level, and assign the average redundantly
! to all lons at each pole.
!
    call mp_gatherlons_f3d(ptr,k0,k1,i0,i1,j0,j1,1)
!
    if (mytidi==0) then ! only westernmost tasks have global longitudes

      if (j0 == 1) then ! subdomain includes south pole
        fave(:) = 0._r8
!
! Find average of all lons at each level, at first lat equatorward of south pole.
!
        do k=k0,k1
          do i=1,nlon
            fave(k) = fave(k)+ptr(1)%ptr(k,i,j0+1) 
          enddo
          fave(k) = fave(k) / rnlon
        enddo  
        if (debug.and.masterproc) write(iulog,"('setpoles: spole ave(k0:k1)=',/,(8es12.4))") fave
!
! Define south pole in ptmp on subdomains for each tasks in my latitude row 
! (I am SW corner task):
!
        do it=0,ntaski-1
          itask = tasks(itask_table_geo(it,mytidj))%mytid
          lon0 = tasks(itask)%lon0
          lon1 = tasks(itask)%lon1
          do k=k0,k1
            ptr(1)%ptr(k,lon0:lon1,j0) = fave(k) ! all lons get the average
          enddo
        enddo
      endif ! south pole

      if (j1 == nlat) then ! subdomain includes north pole
        fave(:) = 0._r8
!
! Find average of all lons at each level, at first lat equatorward of north pole.
!
        do k=k0,k1
          do i=1,nlon
            fave(k) = fave(k)+ptr(1)%ptr(k,i,j1-1)
          enddo
          fave(k) = fave(k) / rnlon
        enddo  
        if (debug.and.masterproc) write(iulog,"('setpoles: npole fave(k0:k1)=',/,(8es12.4))") fave
!
! Define north pole in ptmp on subdomains for each tasks in my latitude row 
! (I am NW corner task):
!
        do it=0,ntaski-1
          itask = tasks(itask_table_geo(it,mytidj))%mytid
          lon0 = tasks(itask)%lon0
          lon1 = tasks(itask)%lon1
          do k=k0,k1
            ptr(1)%ptr(k,lon0:lon1,j1) = fave(k)
          enddo
        enddo
      endif ! north pole
    endif ! mytidj==0
!
! Scatter to tasks in my latitude row:
!
    call mp_scatterlons_f3d(ptr,k0,k1,i0,i1,j0,j1,1)
!
! Define poles on current subdomain inout arg array:
!
    if (j0==1) then
      do i=i0,i1
        do k=k0,k1
          f(k,i,j0) = ptr(1)%ptr(k,i,j0)
        enddo
      enddo
    endif
    if (j1==nlat) then
      do i=i0,i1
        do k=k0,k1
          f(k,i,j1) = ptr(1)%ptr(k,i,j1)
        enddo
      enddo
    endif
    deallocate(ptr(1)%ptr)
  end subroutine setpoles
!-----------------------------------------------------------------------
  subroutine lonshift_blocks(f,k0,k1,i0,i1,j0,j1,nfields)
!
! On input, field(s) f are in subdomains 
! On output, field(s) f subdomain longitudes are shifted by 180 degrees
!   (either 0->360 to -180->+180, or the reverse)
!
    use edyn_geogrid ,only: nlon
!
! Args:
    integer :: k0,k1,i0,i1,j0,j1,nfields
    type(array_ptr_type) :: f(nfields) ! f(n)%ptr(k0:k1,i0:i1,j0:j1)
!
! Local variables
!
    integer :: i,j,k,ifield
    integer :: midpoint            ! middle point of longitude dimension
    real(r8) :: flons(nlon)        ! fields at global longitudes
    type(array_ptr_type) :: pglblon(nfields) ! pglblon(n)%ptr(k0:k1,nlon,j0:j1)
!
! Shift longitude grid from 0 to 360 to -180 to 180 for edynamo
!  Check for compatible geographic longitude dimension and quit if not compatible
!
    if (nlon /= 288 .and. nlon /= 144 .and. nlon /= 80 .and. nlon /= 72 .and. nlon /= 24) then
      write(iulog,"('ERROR lonshift_blocks: incompatible nlon = ',i5,' i0,i1=',2i4)") nlon,i0,i1
      call endrun
    end if
!
! Load subdomains into local global longitude pointer:
    do ifield=1,nfields
      allocate(pglblon(ifield)%ptr(k0:k1,nlon,j0:j1)) 
      do j=j0,j1
        do i=i0,i1
          pglblon(ifield)%ptr(k0:k1,i,j) = f(ifield)%ptr(k0:k1,i,j)
        enddo
      enddo
    enddo

    call mp_gatherlons_f3d(pglblon,k0,k1,i0,i1,j0,j1,nfields)
!
! Only leftmost tasks (mytidi=0) at each latitude does the longitude shift for that latitude
!
    if (mytidi==0) then
      do j=j0,j1
        midpoint = nlon/2
        do ifield = 1,nfields
          do k = k0,k1
            flons(:) = pglblon(ifield)%ptr(k,1:nlon,j)
            flons = cshift(flons,midpoint)
            pglblon(ifield)%ptr(k,1:nlon,j) = flons(:)
          enddo ! k0,k1
        enddo ! nfields
      enddo ! j=j0,j1
    endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
    call mp_scatterlons_f3d(pglblon,k0,k1,i0,i1,j0,j1,nfields)
!
! Update fields argument: 
    do ifield=1,nfields
      do j=j0,j1
        do i=i0,i1
          f(ifield)%ptr(k0:k1,i,j) = pglblon(ifield)%ptr(k0:k1,i,j)
        enddo
      enddo
    enddo

    do ifield=1,nfields
      deallocate(pglblon(ifield)%ptr)
    enddo
  end subroutine lonshift_blocks
!-----------------------------------------------------------------------
  subroutine switch_model_format(fptr,k0,k1,i0,i1,j0,j1,nfields)
!
! fptr is array of pointer structures to nfields fields. Convert these 
! fields in "model format", i.e., phase shift longitude data by 180 degrees,
! and invert the vertical dimension. This may be converting from WACCM to 
! TIEGCM, or the reverse. It is up to the calling routine to keep track of 
! which model format the data is being converted from/to.
! (This routine does not do unit conversion on the fields)
!
! Args:
    integer,intent(in)  :: k0,k1,i0,i1,j0,j1,nfields
!
! Pointer structures to each field:
    type(array_ptr_type) :: fptr(nfields) ! (fptr(n)%ptr(k0:k1,i0:i1,j0:j1))
!
! Local:
    integer :: ifield
!
! Phase shift longitudes by 180 degrees:
!
    call lonshift_blocks(fptr,k0,k1,i0,i1,j0,j1,nfields)
!
! Invert vertical dimension:
!
    do ifield=1,nfields
      fptr(ifield)%ptr(k0:k1,i0:i1,j0:j1) = fptr(ifield)%ptr(k1:k0:-1,i0:i1,j0:j1)
    enddo
  end subroutine switch_model_format
!-----------------------------------------------------------------------
end module edyn_mpi
