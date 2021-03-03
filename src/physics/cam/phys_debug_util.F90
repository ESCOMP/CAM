module phys_debug_util

!----------------------------------------------------------------------------------------

! Module to facilitate debugging of physics parameterizations.
!
! The user requests a location for debugging in lat/lon coordinates
! (degrees).  The initialization routine does a global search to find the
! column in the physics grid closest to the requested location.  The local
! indices of that column in the physics decomposition are stored as module
! data.  The user code then passes the local chunk index of the chunked
! data into the subroutine that will write diagnostic information for the
! column.  The function phys_debug_col returns the local column index if
! the column of interest is contained in the chunk, and zero otherwise.
! Printing is done only if a column index >0 is returned.
!
! Phil Rasch, B. Eaton, Feb 2008
!----------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc, iam, mpicom, npes
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

implicit none
private
save

real(r8), parameter :: uninitr8 = huge(1._r8)

! Public methods
public phys_debug_readnl  ! read namelist input
public phys_debug_init    ! initialize the method to a chunk and column
public phys_debug_col     ! return local column index in debug chunk

! Namelist variables
real(r8) :: phys_debug_lat = uninitr8 ! latitude of requested debug column location in degrees
real(r8) :: phys_debug_lon = uninitr8 ! longitude of requested debug column location in degrees


integer :: debchunk = -999            ! local index of the chuck we will debug
integer :: debcol   = -999            ! the column within the chunk we will debug

!================================================================================
contains
!================================================================================

subroutine phys_debug_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'phys_debug_readnl'

   namelist /phys_debug_nl/ phys_debug_lat, phys_debug_lon
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'phys_debug_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, phys_debug_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
      ! Check inputs
      if (phys_debug_lat /= uninitr8) then
         if (abs(phys_debug_lat) > 90.0_r8) then
            write(iulog, *) subname, ': phys_debug_lat out of range [-90., 90.]'
            call endrun(subname//': phys_debug_lat out of range [-90., 90.]')
         end if
      else
         write(iulog, *) subname, ': phys_debug_lat = ', phys_debug_lat
      end if
      if (phys_debug_lon /= uninitr8) then
         if ((phys_debug_lon < 0.0_r8) .or. (phys_debug_lon > 360.0_r8)) then
            write(iulog, *) subname, ': phys_debug_lon out of range [0., 360.]'
            call endrun(subname//': phys_debug_lon out of range [0., 360.]')
         end if
      else
         write(iulog, *) subname, ': phys_debug_lon = ', phys_debug_lon
      end if
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(phys_debug_lat, 1, mpir8, 0, mpicom)
   call mpibcast(phys_debug_lon, 1, mpir8, 0, mpicom)
#endif

end subroutine phys_debug_readnl

!==============================================================================

subroutine phys_debug_init()
   use mpi,       only: mpi_real8, mpi_integer, mpi_min, mpi_max
   use physconst, only: pi
   use ppgrid,    only: begchunk, endchunk
   use phys_grid, only: get_ncols_p, get_rlat_p, get_rlon_p

   integer             :: owner, lchunk, icol, ncol
   integer             :: lchunk_min, icol_min, minlondist
   real(r8)            :: deblat, deblon
   real(r8)            :: latmin, lonmin
   real(r8)            :: lat, lon, dist, temp1, temp2
   real(r8)            :: mindist
   real(r8), parameter :: maxangle = pi / 4.0_r8
   real(r8), parameter :: rad2deg = 180.0_r8 / pi
   real(r8), parameter :: deg2rad = pi / 180.0_r8
   real(r8), parameter :: maxtol = 0.99999_r8 ! max cos value
   real(r8), parameter :: maxlat = pi * maxtol / 2.0_r8
   !---------------------------------------------------------------------------

   ! If no debug column specified then do nothing
   if ((phys_debug_lat == uninitr8) .or. (phys_debug_lon == uninitr8)) then
      return
   end if

   ! User has specified a column location for debugging.  Find the closest
   ! column in the physics grid.
   mindist = 2.0_r8 * pi
   deblat = pi
   deblon = 3.0_r8 * pi
   latmin = phys_debug_lat * deg2rad
   lonmin = phys_debug_lon * deg2rad
   lchunk_min = -1
   icol_min = -1
   do lchunk = begchunk, endchunk
      ncol = get_ncols_p(lchunk)
      do icol = 1, ncol
         lat = get_rlat_p(lchunk, icol)
         lon = get_rlon_p(lchunk, icol)
         if ( (abs(lat - latmin) <= maxangle) .and.                           &
              (abs(lon - lonmin) <= maxangle)) then
            ! maxangle could be pi but why waste all those trig functions?
            if ((lat == latmin) .and. (lon == lonmin)) then
               dist = 0.0_r8
            else
               temp1 = (sin(latmin) * sin(lat)) +                             &
                    (cos(latmin) * cos(lat) * cos(lon - lonmin))
               if (temp1 > maxtol) then
                  ! Use haversine formula
                  temp1 = sin((latmin - lat) / 2.0_r8)
                  temp2 = sin((lonmin - lon) / 2.0_r8)
                  dist = (temp1 * temp1) +                                    &
                       (cos(latmin)* cos(lat) * temp2 * temp2)
                  dist = 2.0_r8 * asin(sqrt(dist))
               else
                  dist = acos(temp1)
               end if
            end if
            if ( (dist < mindist) .or.                                        &
                 ((dist == mindist) .and.                                     &
                 (abs(lon - lonmin) < abs(deblon - lonmin)))) then
               lchunk_min = lchunk
               icol_min = icol
               mindist = dist
               deblon = lon
               deblat = lat
               if (dist == 0.0_r8) then
                  exit
               end if
            end if
         end if
      end do
   end do
   ! We need to find the minimum mindist and use only that value
   dist = mindist
   call MPI_allreduce(dist, mindist, 1, mpi_real8, mpi_min, mpicom, icol)
   ! Special case for pole points
   if (deblat > pi / 2.0_r8) then
      temp1 = 0.0_r8
   else
      temp1 = abs(deblat)
   end if
   call MPI_allreduce(temp1, lat, 1, mpi_real8, mpi_max, mpicom, icol)
   if ((abs(latmin) > maxlat) .or. (lat > maxlat)) then
      if (dist == mindist) then
         ! Only distance winners can compete
         lon = abs(deblon - lonmin)
      else
         lon = 3.0_r8 * pi
      end if
      call MPI_allreduce(lon, minlondist, 1, mpi_real8, mpi_min, mpicom, icol)
      ! Kill the losers
      if (lon /= minlondist) then
         dist = dist + 1.0_r8
      end if
   end if
   ! Now, only task(s) which have real minimum distance should be owner
   if (dist == mindist) then
      lchunk = iam
   else
      lchunk = npes + 2
   end if
   call MPI_allreduce(lchunk, owner, 1, mpi_integer, mpi_min, mpicom, icol)
   ! If the column is owned by this process then save its local indices
   if (iam == owner) then
      debchunk         = lchunk_min
      debcol           = icol_min
      deblat           = get_rlat_p(lchunk_min, icol_min) * rad2deg
      deblon           = get_rlon_p(lchunk_min, icol_min) * rad2deg
      write(iulog,*) 'phys_debug_init: debugging column at lat=', deblat, '  lon=', deblon
   end if

end subroutine phys_debug_init

!================================================================================

integer function phys_debug_col(chunk)

   integer,  intent(in) :: chunk
   !-----------------------------------------------------------------------------

   if (chunk == debchunk) then
      phys_debug_col = debcol
   else
      phys_debug_col = 0
   endif

end function phys_debug_col

!================================================================================

end module phys_debug_util
