!-------------------------------------------------------------------------------
! Encapsulates an ESMF regular longitude / latitude grid
!-------------------------------------------------------------------------------
module esmf_lonlat_grid_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils, only: masterproc, mpicom
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun

  use ESMF, only: ESMF_Grid, ESMF_GridCreate1PeriDim, ESMF_GridAddCoord
  use ESMF, only: ESMF_GridGetCoord, ESMF_GridDestroy
  use ESMF, only: ESMF_KIND_R8, ESMF_INDEX_GLOBAL, ESMF_STAGGERLOC_CENTER
  use esmf_check_error_mod, only: check_esmf_error

  implicit none

  public

  type(ESMF_Grid), protected :: lonlat_grid

  integer, protected :: nlon = 0
  integer, protected :: nlat = 0

  integer, protected :: lon_beg = -1
  integer, protected :: lon_end = -1
  integer, protected :: lat_beg = -1
  integer, protected :: lat_end = -1

  real(r8), allocatable, protected :: glats(:)
  real(r8), allocatable, protected :: glons(:)

  integer, protected :: zonal_comm ! zonal direction MPI communicator

contains

  subroutine esmf_lonlat_grid_init(nlats_in)
    use phys_grid, only: get_grid_dims
    use mpi, only: mpi_comm_size, mpi_comm_rank, MPI_PROC_NULL, MPI_INTEGER

    integer, intent(in) :: nlats_in

    real(r8) :: delx, dely

    integer :: npes, ierr, mytid, irank, mytidi, mytidj
    integer :: i,j, n
    integer :: ntasks_lon, ntasks_lat
    integer :: lons_per_task, lons_overflow, lats_per_task, lats_overflow
    integer :: task_cnt
    integer :: mynlats, mynlons

    integer, allocatable :: mytidi_send(:)
    integer, allocatable :: mytidj_send(:)
    integer, allocatable :: mytidi_recv(:)
    integer, allocatable :: mytidj_recv(:)

    integer, allocatable :: nlons_send(:)
    integer, allocatable :: nlats_send(:)
    integer, allocatable :: nlons_recv(:)
    integer, allocatable :: nlats_recv(:)

    integer, allocatable :: nlons_task(:)
    integer, allocatable :: nlats_task(:)

    integer, parameter :: minlats_per_pe = 2
    integer, parameter :: minlons_per_pe = 2

    integer, allocatable :: petmap(:,:,:)
    integer :: petcnt

    integer                       :: lbnd_lat, ubnd_lat, lbnd_lon, ubnd_lon
    integer                       :: lbnd(1), ubnd(1)
    real(ESMF_KIND_R8), pointer   :: coordX(:), coordY(:)

    character(len=*), parameter :: subname  = 'esmf_lonlat_grid_init: '

    ! create reg lon lat grid

    nlat = nlats_in
    dely = 180._r8/nlat

    nlon = 2*nlat
    delx = 360._r8/nlon

    allocate(glons(nlon))
    allocate(glats(nlat))

    glons(1) = 0._r8
    glats(1) = -90._r8 + 0.5_r8 * dely

    do i = 2,nlon
       glons(i) = glons(i-1) + delx
    end do
    do i = 2,nlat
       glats(i) = glats(i-1) + dely
    end do

    ! decompose the grid across mpi tasks ...

    call mpi_comm_size(mpicom, npes, ierr)
    call mpi_comm_rank(mpicom, mytid, ierr)

    decomp_loop: do ntasks_lon = 1,nlon
       ntasks_lat = npes/ntasks_lon
       if ( (minlats_per_pe*ntasks_lat<nlat) .and. (ntasks_lat*ntasks_lon==npes) ) then
          exit decomp_loop
       endif
    end do decomp_loop

    if (masterproc) then
       write(iulog,'(a,3i6)') subname//' npes, nlon, nlat : ',npes, nlon, nlat
       write(iulog,'(a,2i6)') subname//' ntasks_lon,ntasks_lat : ',ntasks_lon,ntasks_lat
    endif

    if (ntasks_lat*ntasks_lon/=npes) then
       call endrun(subname//'ntasks_lat*ntasks_lon/=npes')
    endif

    ! dermine the starting and ending coordinates
    lons_per_task = nlon / ntasks_lon
    lons_overflow = MOD(nlon, ntasks_lon)
    lats_per_task = nlat / ntasks_lat
    lats_overflow = MOD(nlat, ntasks_lat)
    lon_beg = 1
    lon_end = 0
    lat_beg = 1
    lat_end = 0
    task_cnt= 0
    if (mytid<npes) then
       jloop: do j = 0,ntasks_lat-1
          lat_beg = lat_end + 1
          lat_end = lat_beg + lats_per_task - 1
          if (j<lats_overflow) then
             lat_end = lat_end + 1
          end if
          lon_end = 0
          do i = 0,ntasks_lon-1
             lon_beg = lon_end + 1
             lon_end = lon_beg + lons_per_task - 1
             if (i<lons_overflow) then
                lon_end = lon_end + 1
             end if
             task_cnt = task_cnt+1
             if (task_cnt>mytid) exit jloop
          end do
       enddo jloop
    endif

    mynlats = lat_end-lat_beg+1
    mynlons = lon_end-lon_beg+1

    if (mynlats<minlats_per_pe) then
       call endrun(subname//'mynlats < minlats_per_pe')
    end if
    if (mynlons<minlons_per_pe) then
       call endrun(subname//'mynlons < minlons_per_pe')
    end if

    irank = 0
    mytidi = -1
    mytidj = -1
    do j =  0, ntasks_lat-1
       do i = 0, ntasks_lon-1
          if (mytid == irank) then
             mytidi = i
             mytidj = j
          end if
          irank = irank+1
       end do

    end do ! j=0,ntaskj-1

    call mpi_comm_split(mpicom,mytidj,mytid,zonal_comm,ierr)

    allocate(mytidi_send(npes))
    allocate(mytidj_send(npes))
    allocate(mytidi_recv(npes))
    allocate(mytidj_recv(npes))

    mytidi_send = mytidi
    mytidj_send = mytidj
    call mpi_alltoall(mytidi_send, 1, MPI_INTEGER, mytidi_recv, 1, MPI_INTEGER, mpicom, ierr)
    call mpi_alltoall(mytidj_send, 1, MPI_INTEGER, mytidj_recv, 1, MPI_INTEGER, mpicom, ierr)

    allocate(nlons_send(npes))
    allocate(nlats_send(npes))
    allocate(nlons_recv(npes))
    allocate(nlats_recv(npes))

    nlons_send(:) = mynlons
    nlats_send(:) = mynlats

    call mpi_alltoall(nlons_send, 1, MPI_INTEGER, nlons_recv, 1, MPI_INTEGER, mpicom, ierr)
    call mpi_alltoall(nlats_send, 1, MPI_INTEGER, nlats_recv, 1, MPI_INTEGER, mpicom, ierr)

    deallocate(nlons_send)
    deallocate(nlats_send)

    allocate(nlons_task(ntasks_lon))
    allocate(nlats_task(ntasks_lat))

    do i = 1, ntasks_lon
       loop1: do n = 1, npes
          if (mytidi_recv(n) == i-1) then
             nlons_task(i) = nlons_recv(n)
             exit loop1
          end if
       end do loop1
    end do

    do j = 1, ntasks_lat
       loop2: do n = 1, npes
          if (mytidj_recv(n) == j-1) then
             nlats_task(j) = nlats_recv(n)
             exit loop2
          end if
       end do loop2
    end do

    deallocate(mytidi_send)
    deallocate(mytidj_send)
    deallocate(mytidi_recv)
    deallocate(mytidj_recv)

    deallocate(nlons_recv)
    deallocate(nlats_recv)


    ! set up 2D ESMF lon lat grid

    allocate(petmap(ntasks_lon,ntasks_lat,1))

    petcnt = 0
    do j = 1,ntasks_lat
       do i = 1,ntasks_lon
          petmap(i,j,1) = petcnt
          petcnt = petcnt+1
       end do
    end do


    ! Create 2d lon/lat grid
    lonlat_grid = ESMF_GridCreate1PeriDim(                       &
         countsPerDEDim1=nlons_task, coordDep1=(/1/),         &
         countsPerDEDim2=nlats_task, coordDep2=(/2/), petmap=petmap, &
         indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/), rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_GridCreate1PeriDim ERROR')


    ! Set coordinates:

    call ESMF_GridAddCoord(lonlat_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_GridAddCoord ERROR')

    if (mytid<npes) then
       call ESMF_GridGetCoord(lonlat_grid, coordDim=1, &
            computationalLBound=lbnd, computationalUBound=ubnd,  &
            farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
       call check_esmf_error(ierr, subname//'ESMF_GridGetCoord for longitude coords ERROR')

       lbnd_lon = lbnd(1)
       ubnd_lon = ubnd(1)
       do i = lbnd_lon, ubnd_lon
          coordX(i) = glons(i)
       end do

       call ESMF_GridGetCoord(lonlat_grid, coordDim=2, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
       call check_esmf_error(ierr, subname//'ESMF_GridGetCoord for latitude coords ERROR')

       lbnd_lat = lbnd(1)
       ubnd_lat = ubnd(1)
       do i = lbnd_lat, ubnd_lat
          coordY(i) = glats(i)
       end do
    end if

    deallocate(nlons_task)
    deallocate(nlats_task)

  end subroutine esmf_lonlat_grid_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_lonlat_grid_destroy()

    integer :: rc
    character(len=*), parameter :: subname = 'esmf_lonlat_grid_destroy: '

    call ESMF_GridDestroy(lonlat_grid, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_GridDestroy lonlat_grid')

    deallocate(glats)
    deallocate(glons)

  end subroutine esmf_lonlat_grid_destroy


end module esmf_lonlat_grid_mod
