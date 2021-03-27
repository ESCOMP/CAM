module edyn_geogrid
!
! Global geographic grid.
! See sub set_geogrid (edyn_init.F90)
!
  use shr_kind_mod, only: r8 => shr_kind_r8 ! 8-byte reals
  use cam_logfile,  only: iulog
  use cam_abortutils, only: endrun

  implicit none
  private
  save

  integer, public, protected :: & ! dimensions
    nlat,    & ! number of latitudes
    nlon,    & ! number of longitudes
    nlev,    & ! number of midpoint levels
    nilev,   & ! number of interface levels
    npes       ! number of PEs in geogrid

  real(r8), public, protected, allocatable , dimension(:) :: & ! coordinate vars
    glat,    & ! latitude coordinates (degrees)
    glon,    & ! longitude coordinates (degrees)
    ylatg,   & ! latitudes (radians)
    ylong,   & ! longitudes (radians)
    zlev,    & ! midpoint vertical coordinates
    zilev      ! interface vertical coordinates

  real(r8), public, allocatable, protected :: cs(:)   ! cos(phi) (0:nlat+1)

  integer, public, protected :: & ! model independent (set by sub get_geogrid)
    nlonp1,  & ! nlon+1
    nlonp2,  & ! nlon+2
    nlatp1     ! nlat+1

  real(r8), public, protected :: dphi
  real(r8), public, protected :: dlamda
!
! Using p0 in microbars, as in TIEGCM.
  real(r8), parameter, public :: p0 = 5.0e-4_r8  ! standard pressure (microbars)

  integer, public, protected :: & ! model dependent (set by subs read_tgcm, read_waccm)
    jspole,  & ! latitude index to geographic south pole
    jnpole     ! latitude index to geographic north pole

  ! set_geogrid sets up a distributed finite-volume lat/lon grid
  public :: set_geogrid

  logical :: debug = .false. ! set true for prints to stdout at each call

contains

   !-----------------------------------------------------------------------
   subroutine set_geogrid(nlon_g, nlat_g, nlev_in, npes_in, iam, pres_mid_in, pres_edge_in, min_lat_pe_in)
      use shr_const_mod,  only: pi => shr_const_pi
      use edyn_params,    only: kbotdyn, pbotdyn
      use edyn_mpi,       only: mp_distribute_geo
      use spmd_utils,     only: masterproc
      use edyn_maggrid,   only: nmlat

      ! Dummy Args
      integer,            intent(in) :: nlon_g        ! Global num longitudes
      integer,            intent(in) :: nlat_g        ! Global num latitudes
      integer,            intent(in) :: nlev_in       ! Num levels
      integer,            intent(in) :: npes_in
      integer,            intent(in) :: iam
      real(r8),           intent(in) :: pres_mid_in(:)
      real(r8),           intent(in) :: pres_edge_in(:)
      integer,  optional, intent(in) :: min_lat_pe_in ! Min # lats / PE
      !
      ! Local:
      integer                        :: latind, lonind, js, k
      integer                        :: lon_beg, lon_end, lat_beg, lat_end
      integer                        :: lons_per_task, lats_per_task
      integer                        :: lons_overflow, lats_overflow
      integer                        :: ntasks_lat, ntasks_lon
      integer                        :: task_cnt, i,j
      integer                        :: minlats_per_pe
      integer                        :: ierr
      real(r8)                       :: phi
      real(r8)                       :: delta ! Coordinate spacing
      real(r8), parameter            :: eps = 1.e-6_r8

      real(r8)                       :: pmid(nlev_in)

      nlon = nlon_g
      nlat = nlat_g
      nlev = nlev_in
      npes = npes_in

      nilev = nlev+1

      nlonp1 = nlon + 1
      nlonp2 = nlon + 2
      nlatp1 = nlat + 1

      jspole = 1
      jnpole = nlat

      if (present(min_lat_pe_in)) then
         minlats_per_pe = min_lat_pe_in
      else
         minlats_per_pe = 2
      end if

      dphi   = pi / real(nlat,r8)
      dlamda = 2._r8*pi / real(nlon,r8)

      !
      ! Allocate coordinate variables:
      !
      allocate(glon(nlon))
      allocate(glat(nlat))
      !
      ! Create a finite-volume coordinate grid (in degrees)
      !
      delta = 360.0_r8 / real(nlon, r8)
      do lonind = 1, nlon
         glon(lonind) = -180.0_r8 + ((lonind - 1) * delta)
      end do
      delta = 180.0_r8 / real((nlat - 1), r8)
      ! Set the poles exactly (they might be checked later)
      glat(1) = -90.0_r8
      glat(nlat) = 90.0_r8
      do latind = 2, nlat - 1
         glat(latind) = -90.0_r8 + ((latind - 1) * delta)
      end do

      if (masterproc.and.debug) then
        write(iulog,*) 'set_geogrid glon : ',glon(:)
        write(iulog,*) 'set_geogrid glat : ',glat(:)
      end if

      allocate(zlev(nlev))
      allocate(zilev(nilev))
      !
      ! Hybrid-sigma levels from ref_pres module:
      !
      zlev(:nlev)  = pres_mid_in(:)  ! midpoints vertical coord (top down)
      zilev(:nilev) = pres_edge_in(:nilev)  ! interfaces vertical coord

      ! do bottom up search for kbotdyn
      pmid(:nlev) = zlev(nlev:1:-1)
      kloop: do k = 1, nlev
         if ( pmid(k) <= pbotdyn) then
            kbotdyn = k
            exit kloop
         end if
      end do kloop
      if ( kbotdyn < 1 ) then
         call endrun('set_geogrid: kbotdyn is not set')
      endif
      if (debug) then
         write(iulog,"('set_geogrid: kbotdyn=',i4,' pmid(kbotdyn)=',es12.4)") kbotdyn,pmid(kbotdyn)
      endif

      !
      ! Setup a decomposition for the geogrid
      !
      ! First, try using a 1-D latitude decomposition

      do ntasks_lon = 1,nlon_g
         ntasks_lat = npes/ntasks_lon
         if ( minlats_per_pe*ntasks_lat<nmlat .and. ntasks_lat*ntasks_lon==npes ) then
            exit
         endif
      end do
      if (masterproc) then
         write(iulog,'(a,3i6)') 'set_geogrid: npes,nlon_g,nlat_g: ',npes,nlon_g,nlat_g
         write(iulog,'(a,2i6)') 'set_geogrid: ntasks_lon,ntasks_lat (oplus and edyn grids): ',ntasks_lon,ntasks_lat
      endif

      if (ntasks_lat*ntasks_lon/=npes) then
         call endrun('set_geogrid: ntasks_lat*ntasks_lon/=npes')
      endif

      ! Now, figure the starting and ending coordinates
      lons_per_task = nlon / ntasks_lon
      lons_overflow = MOD(nlon, ntasks_lon)
      lats_per_task = nlat / ntasks_lat
      lats_overflow = MOD(nlat, ntasks_lat)
      lon_beg = 1
      lon_end = 0
      lat_beg = 1
      lat_end = 0
      task_cnt= 0
      if (iam<npes) then
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
               if (task_cnt>iam) exit jloop
            end do
         enddo jloop
      endif

      call mp_distribute_geo(lon_beg, lon_end, lat_beg, lat_end, 1, nlev, ntasks_lon, ntasks_lat)

      !
      ! Set horizontal geographic grid in radians (for apex code):
      !
      allocate(ylatg(nlat))   ! waccm grid includes poles
      allocate(ylong(nlonp1)) ! single periodic point
      ylatg(1)    = -pi/2._r8+eps ! south pole
      ylatg(nlat) =  pi/2._r8-eps ! north pole
      do latind = 2, nlat-1
         ylatg(latind) = -0.5_r8*(pi-dphi)+real(latind-1,r8)*dphi
      end do
      do lonind = 1, nlonp1
         ylong(lonind) = -pi+real(lonind-1,r8)*dlamda
      end do
      !
      ! Calculate cosine of latitude
      !
      allocate(cs(0:nlat+1))
      js = -(nlat/2)
      do latind = 1, nlat
         phi = (latind + js - .5_r8) * dphi
         cs(latind) = cos(phi)
      end do
      cs(0) = -cs(1)
      cs(nlat+1) = -cs(nlat)

   end subroutine set_geogrid

end module edyn_geogrid
