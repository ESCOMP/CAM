module edyn_geogrid
!
! Global geographic grid. 
! See sub set_geogrid (edyn_init.F90)
!
  use shr_kind_mod,   only : r8 => shr_kind_r8            ! 8-byte reals
  implicit none
  save

  integer :: & ! dimensions
    nlat,    & ! number of latitudes
    nlon,    & ! number of longitudes
    nlev,    & ! number of midpoint levels
    nilev,   & ! number of interface latitudes
    ntime      ! number of times on input file

  real(r8),allocatable,dimension(:) :: & ! coordinate vars
    glat,    & ! latitude coordinates (degrees)
    glon,    & ! longitude coordinates (degrees)
    ylatg,   & ! latitudes (radians)
    ylong,   & ! longitudes (radians)
    zlev,    & ! midpoint vertical coordinates
    zilev,   & ! interface vertical coordinates
    time       ! times (histories) on input file

  real(r8),allocatable,dimension(:) :: & 
    cs,      & ! cos(phi) (0:nlat+1)
    zp,      & ! log pressure (as in tiegcm lev(nlev))
    expz       ! exp(-zp)

  integer :: & ! model independent (set by sub get_geogrid)
    nlonp1,  & ! nlon+1
    nlonp2,  & ! nlon+2
    nlatp1     ! nlat+1

  real(r8) :: dlatg,dlong
  real(r8) :: dphi,dlamda
!
! Using p0 in microbars, as in TIEGCM.
  real(r8),parameter :: p0 = 5.0e-4_r8  ! standard pressure (microbars)

  integer :: & ! model dependent (set by subs read_tgcm, read_waccm)
    jspole,  & ! latitude index to geographic south pole
    jnpole     ! latitude index to geographic north pole
!
! lev_sequence is a string indicating ordering of the vertical 
! coordinates lev and ilev, and of the field arrays along the 
! vertical dimension. lev_sequence can have 1 of 2 values:
!
!  'bottom2top' means lev(1) is the bottom boundary, lev(nlev) is the top boundary 
!  'top2bottom' means lev(1) is the top boundary, lev(nlev) is the bottom boundary 
!
! For example, TIMEGCM history files are bottom2top, whereas 
! WACCM files are top2bottom. The edynamo code assumes bottom2top, 
! so WACCM input fields are reversed to be bottom2top for the edynamo 
! calculations, then reversed back to the native WACCM sequence 
! (top2bottom) before writing to the edynamo output file.
!   
  character(len=10) :: lev_sequence
!
! lon_sequence is a string indicating ordering of the longitude
! coordinate lon, and of the field arrays along this dimension.
! lon_sequece can have 1 of 2 values:
!
!   '-180to180' means lon(1) is -180 deg west longitude, lon(nlon) is +180 east
!   'zeroto360' means lon(1) is 0 deg west longitude, lon(nlon) is 360 deg east
!
! Note that TIMEGCM convention is '-180to180' and WACCM convention is 'zeroto360'
! (this is treating similarly to lev_sequence above)
!
  character(len=9) :: lon_sequence

end module edyn_geogrid
