module dyn_grid
!-----------------------------------------------------------------------
!
! Define grid and decomposition for Eulerian spectral dynamics.
!
! Original code: John Drake and Patrick Worley
!
!-----------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use pmgrid,           only: plat, plev, plon, plevp
use physconst,        only: rair, rearth, ra
use spmd_utils,       only: masterproc, iam

use pio,              only: file_desc_t
use cam_initfiles,    only: initial_file_get_id

use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

#if (defined SPMD)
use spmd_dyn,         only: spmdinit_dyn
#endif

implicit none
private
save

public :: &
   dyn_grid_init,            &
   dyn_grid_find_gcols,      &! find nearest column for given lat/lon
   dyn_grid_get_colndx,      &! global lat and lon coordinate and MPI process indices
                              ! corresponding to a specified global column index
   dyn_grid_get_elem_coords, &! coordinates of a specified element (latitude)
                              ! of the dynamics grid (lat slice of the block)
   get_block_bounds_d,       &! first and last indices in global block ordering
   get_block_gcol_d,         &! global column indices for given block
   get_block_gcol_cnt_d,     &! number of columns in given block
   get_block_levels_d,       &! vertical levels in column
   get_block_lvl_cnt_d,      &! number of vertical levels in column
   get_block_owner_d,        &! process "owning" given block
   get_dyn_grid_parm,        &
   get_dyn_grid_parm_real1d, &
   get_gcol_block_d,         &! global block indices and local columns
                              ! index for given global column index
   get_gcol_block_cnt_d,     &! number of blocks containing data
                              ! from a given global column index
   get_horiz_grid_d,         &! horizontal grid coordinates
   get_horiz_grid_dim_d,     &! horizontal dimensions of dynamics grid
   physgrid_copy_attributes_d

! The Eulerian dynamics grids
integer, parameter, public :: dyn_decomp       = 101

integer, parameter, public :: ptimelevels = 3  ! number of time levels in the dycore

integer :: ngcols_d = 0     ! number of dynamics columns

!========================================================================================
contains
!========================================================================================

subroutine dyn_grid_init

   ! Initialize dynamics grid

   use pspect,          only: ptrm, ptrn, ptrk, pnmax, pmmax, pspt
   use comspe,          only: lpspt, numm, locm, lnstart, nstart, nlen, &
                              alp, dalp, lalp, ldalp
   use scanslt,         only: nlonex, platd, j1
   use gauaw_mod,       only: gauaw
   use commap,          only: sq, rsq, slat, w, cs, href, ecref, clat, clon, &
                              latdeg, londeg, xm
   use time_manager,    only: get_step_size
   use scamMod,         only: scmlat, scmlon, single_column
   use hycoef,          only: hycoef_init, hypi, hypm, hypd, nprlev
   use ref_pres,        only: ref_pres_init
   use eul_control_mod, only: ifax, trig, eul_nsplit

   ! Local variables
   type(file_desc_t), pointer :: fh_ini

   real(r8) zsi(plat)      ! sine of latitudes
   real(r8) zw(plat)       ! Gaussian weights
   real(r8) zra2           ! ra squared
   real(r8) zalp(2*pspt)   ! Legendre function array
   real(r8) zdalp(2*pspt)  ! Derivative array
   real(r8) zslat          ! sin of lat  and cosine of colatitude

   integer i           ! longitude index
   integer j           ! Latitude index
   integer k           ! Level index
   integer kk          ! Level index
   integer kkk         ! Level index
   integer m,lm,mr,lmr ! Indices for legendre array
   integer n           ! Index for legendre array
   integer nkk         ! Print control variables
   integer ik1         ! Print index temporary variable
   integer ik2         ! Print index temporary variable
   integer itmp        ! Dimension of polynomial arrays temporary.
   integer iter        ! Iteration index
   real(r8) :: zdt     ! Time step for settau

   integer :: irow        ! Latitude pair index
   integer :: lat         ! Latitude index

   real(r8) :: xlat    ! Latitude (radians)
   real(r8) :: pi      ! Mathematical pi (3.14...)
   real(r8) :: dtime   ! timestep size [seconds]

   character(len=*), parameter :: sub='dyn_grid_init'
   !-----------------------------------------------------------------------

   ! File handle for initial file.  Needed for vertical coordinate data.
   fh_ini => initial_file_get_id()

   ! Compute truncation parameters
   call trunc()

#if (defined SPMD)
   call spmdinit_dyn()
#endif

   ! Initialize hybrid coordinate arrays
   call hycoef_init(fh_ini)

   ! Initialize reference pressures
   call ref_pres_init(hypi, hypm, nprlev)


   dtime = get_step_size()
   zdt = dtime/eul_nsplit

   ! Initialize horizontal diffusion coefficients
   call hdinti(rearth, zdt)

   if (.not. single_column) then

      if (pmmax > plon/2) then
         call endrun (sub//': ERROR: mmax=ptrm+1 .gt. plon/2')
      end if
   end if

   ! NMAX dependent arrays
   zra2 = ra*ra
   do j = 2, pnmax
      sq(j)  = j*(j-1)*zra2
      rsq(j) = 1._r8/sq(j)
   end do
   sq(1)  = 0._r8
   rsq(1) = 0._r8

   ! MMAX dependent arrays
   do j = 1, pmmax
      xm(j) = j-1
   end do

   ! Integration matrices of hydrostatic equation(href) and conversion
   ! term(a).  href computed as in ccm0 but isothermal bottom ecref
   ! calculated to conserve energy

   do k = 1, plev
      do kk = 1, plev
         href(kk,k)  = 0._r8
         ecref(kk,k) = 0._r8
      end do
   end do

   ! Mean atmosphere energy conversion term is consistent with continiuty
   ! Eq.  In ecref, 1st index = column; 2nd index = row of matrix.
   ! Mean atmosphere energy conversion term is energy conserving

   do k = 1, plev
      ecref(k,k) = 0.5_r8/hypm(k) * hypd(k)
      do kk = 1, k-1
         ecref(kk,k) = 1._r8/hypm(k) * hypd(kk)
      end do
   end do

   ! Reference hydrostatic integration matrix consistent with conversion
   ! term for energy conservation.  In href, 1st index = column;
   ! 2nd index = row of matrix.

   do k = 1, plev
      do kk = k, plev
         href(kk,k) = ecref(k,kk)*hypd(kk)/hypd(k)
      end do
   end do

   href = href*rair

   if (single_column) then

      do j = 1, plat
         slat(j) = 1.0_r8 * sin(4.0_r8*atan(1.0_r8)*scmlat/180._r8)
         w(j)    = 2.0_r8/plat
         cs(j)   = 10._r8 - slat(j)*slat(j)
      end do

      xlat = asin(slat(1))
      clat(1) = xlat

      clat(1)=scmlat*atan(1._r8)/45._r8
      latdeg(1) = clat(1)*45._r8/atan(1._r8)
      clon(1,1)   = 4.0_r8*atan(1._r8)*mod((scmlon+360._r8),360._r8)/180._r8
      londeg(1,1) = mod((scmlon+360._r8),360._r8)

   else

      ! Gaussian latitude dependent arrays
      call gauaw(zsi, zw, plat)
      do irow = 1, plat/2
         slat(irow)        = zsi(irow)
         w(irow)           = zw(irow)
         w(plat-irow+1)    = zw(irow)
         cs(irow)          = 1._r8 - zsi(irow)*zsi(irow)
         xlat              = asin(slat(irow))
         clat(irow)        = -xlat
         clat(plat-irow+1) = xlat
      end do

      do lat = 1, plat
         latdeg(lat) = clat(lat)*45._r8/atan(1._r8)
      end do

      ! Compute constants related to Legendre transforms
      ! Compute and reorder ALP and DALP

      allocate(alp(pspt,plat/2))
      allocate(dalp(pspt,plat/2))

      do j = 1, plat/2
         zslat = slat(j)
         itmp  = 2*pspt - 1
         call phcs(zalp, zdalp, itmp, zslat)
         call reordp(j, itmp, zalp, zdalp)
      end do

      ! Copy and save local ALP and DALP

      allocate(lalp(lpspt,plat/2))
      allocate(ldalp(lpspt,plat/2))

      do j = 1, plat/2
         do lm = 1, numm(iam)
            m   = locm(lm,iam)
            mr  = nstart(m)
            lmr = lnstart(lm)
            do n = 1, nlen(m)
               lalp(lmr+n,j)  = alp(mr+n,j)
               ldalp(lmr+n,j) = dalp(mr+n,j)
            end do
         end do
      end do

      ! Mirror latitudes south of south pole

      lat = 1
      do j = j1-2, 1, -1
         nlonex(j) = plon
         lat = lat + 1
      end do
      nlonex(j1-1) = plon     ! south pole

      ! Real latitudes

      j = j1
      do lat = 1, plat
         nlonex(j) = plon
         j = j + 1
      end do
      nlonex(j1+plat) = plon  ! north pole

      ! Mirror latitudes north of north pole

      lat = plat
      do j = j1+plat+1, platd
         nlonex(j) = plon
         lat = lat - 1
      end do

      ! Longitude array

      pi = 4.0_r8*atan(1.0_r8)
      do lat = 1, plat
         do i = 1, plon
            londeg(i,lat) = (i-1)*360._r8/plon
            clon(i,lat)   = (i-1)*2.0_r8*pi/plon
         end do
      end do

      ! Set up trigonometric tables for fft

      do j = 1, plat
         call set99(trig(1,j), ifax(1,j), plon)
      end do
   end if

   ! Define the CAM grids (must be before addfld calls)
   call define_cam_grids()

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'EULERIAN dycore -- Done grid and decomposition initialization'
      write(iulog,*) '  Truncation Parameters: M =',ptrm,'  N =',ptrn,'  K =',ptrk
      write(iulog,*) '  zdt, dtime=', zdt, dtime
      write(iulog,*) ' '
   end if

end subroutine dyn_grid_init

!========================================================================================

   subroutine get_block_bounds_d(block_first,block_last)

!-----------------------------------------------------------------------
!
!
! Purpose: Return first and last indices used in global block ordering
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------
   use pmgrid, only: plat

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(out) :: block_first  ! first (global) index used for blocks
   integer, intent(out) :: block_last   ! last (global) index used for blocks

!-----------------------------------------------------------------------
!  latitude slice block
   block_first = 1
   block_last  = plat

   return
   end subroutine get_block_bounds_d

!
!========================================================================
!
   subroutine get_block_gcol_d(blockid,size,cdex)

!-----------------------------------------------------------------------
!
!
! Purpose: Return list of dynamics column indices in given block
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------
   use pmgrid,     only: plat, plon

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid      ! global block id
   integer, intent(in) :: size         ! array size

   integer, intent(out):: cdex(size)   ! global column indices
!---------------------------Local workspace-----------------------------
!
    integer i,j                            ! loop indices
    integer n                              ! column index
!-----------------------------------------------------------------------
! block == latitude slice
   if (size < plon) then
      write(iulog,*)'GET_BLOCK_GCOL_D: array not large enough (', &
                          size,' < ',plon,' ) '
      call endrun
   else
      n = (blockid-1)*plon
      do i = 1,plon
         n = n + 1
         cdex(i) = n
      end do
   end if
!
   return
   end subroutine get_block_gcol_d
!
!========================================================================
!
   integer function get_block_gcol_cnt_d(blockid)

!-----------------------------------------------------------------------
!
!
! Purpose: Return number of dynamics columns in indicated block
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------
   use pmgrid, only: plon

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
   get_block_gcol_cnt_d = plon

   return
   end function get_block_gcol_cnt_d

!
!========================================================================
!
   integer function get_block_lvl_cnt_d(blockid,bcid)

!-----------------------------------------------------------------------
!
!
! Purpose: Return number of levels in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block

!-----------------------------------------------------------------------
!  latitude slice block
   get_block_lvl_cnt_d = plev + 1

   return
   end function get_block_lvl_cnt_d
!
!========================================================================
!
   subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

!-----------------------------------------------------------------------
!
!
! Purpose: Return level indices in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block
   integer, intent(in) :: lvlsiz   ! dimension of levels array

   integer, intent(out) :: levels(lvlsiz) ! levels indices for block

!---------------------------Local workspace-----------------------------
!
    integer k                      ! loop index
!-----------------------------------------------------------------------
!  latitude slice block
   if (lvlsiz < plev + 1) then
      write(iulog,*)'GET_BLOCK_LEVELS_D: levels array not large enough (', &
                          lvlsiz,' < ',plev + 1,' ) '
      call endrun
   else
      do k=0,plev
         levels(k+1) = k
      end do
      do k=plev+2,lvlsiz
         levels(k) = -1
      end do
   end if

   return
   end subroutine get_block_levels_d

!
!========================================================================
!
   subroutine get_gcol_block_d(gcol,cnt,blockid,bcid,localblockid)

!-----------------------------------------------------------------------
!
!
! Purpose: Return global block index and local column index
!          for global column index
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------
   use pmgrid,     only: plat, plon

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: gcol     ! global column index
   integer, intent(in) :: cnt      ! size of blockid and bcid arrays

   integer, intent(out) :: blockid(cnt) ! block index
   integer, intent(out) :: bcid(cnt)    ! column index within block
   integer, intent(out), optional :: localblockid(cnt)
!---------------------------Local workspace-----------------------------
!
    integer jb                     ! loop index
!-----------------------------------------------------------------------
!  latitude slice block
   if (cnt < 1) then
      write(iulog,*)'GET_GCOL_BLOCK_D: arrays not large enough (', &
                          cnt,' < ',1,' ) '
      call endrun
   else
      blockid(1) = (gcol-1)/plon + 1
      bcid(1)    = gcol - (blockid(1)-1)*plon
      do jb=2,cnt
         blockid(jb) = -1
         bcid(jb)    = -1
      end do
   end if
!
   return
   end subroutine get_gcol_block_d
!
!========================================================================
!
   integer function get_gcol_block_cnt_d(gcol)

!-----------------------------------------------------------------------
!
!
! Purpose: Return number of blocks contain data for the vertical column
!          with the given global column index
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: gcol     ! global column index
!-----------------------------------------------------------------------
!  latitude slice block
   get_gcol_block_cnt_d = 1

   return
   end function get_gcol_block_cnt_d
!
!========================================================================
!
   integer function get_block_owner_d(blockid)

!-----------------------------------------------------------------------
!
!
! Purpose: Return id of processor that "owns" the indicated block
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------
#if ( defined SPMD )
   use spmd_dyn, only: proc
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
#if (defined SPMD)
   get_block_owner_d = proc(blockid)
#else
   get_block_owner_d = 0
#endif

   return
   end function get_block_owner_d
!
!========================================================================
!
   subroutine get_horiz_grid_dim_d(hdim1_d,hdim2_d)

!-----------------------------------------------------------------------
!
!
! Purpose: Returns declared horizontal dimensions of computational grid.
!          Note that global column ordering is assumed to be compatible
!          with the first dimension major ordering of the 2D array.
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------
   use pmgrid,     only: plat, plon

!------------------------------Arguments--------------------------------
   integer, intent(out) :: hdim1_d           ! first horizontal dimension
   integer, intent(out) :: hdim2_d           ! second horizontal dimension
!-----------------------------------------------------------------------
   if (ngcols_d == 0) then
      ngcols_d = plat*plon
   end if
   hdim1_d = plon
   hdim2_d = plat

   return
   end subroutine get_horiz_grid_dim_d
!
!========================================================================
!
   subroutine get_horiz_grid_d(size,clat_d_out,clon_d_out,area_d_out, &
                               wght_d_out,lat_d_out,lon_d_out)

!-----------------------------------------------------------------------
!
!
! Purpose: Return latitude and longitude (in radians), column surface
!          area (in radians squared) and surface integration weights
!          for global column indices that will be passed to/from physics
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------
   use pmgrid,        only: plat, plon
   use commap,        only: clat, clon, londeg, latdeg, w
   use physconst,     only: pi, spval
   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in)   :: size             ! array sizes

   real(r8), intent(out), optional :: clat_d_out(size) ! column latitudes
   real(r8), intent(out), optional :: clon_d_out(size) ! column longitudes
   real(r8), intent(out), optional :: area_d_out(size) ! column surface
                                                       !  area
   real(r8), intent(out), optional :: wght_d_out(size) ! column integration
                                                       !  weight
   real(r8), intent(out), optional :: lat_d_out(size)  ! column deg latitudes
   real(r8), intent(out), optional :: lon_d_out(size)  ! column deg longitudes
!---------------------------Local workspace-----------------------------
!
    integer i,j                      ! loop indices
    integer n                        ! column index
    real(r8) :: ns_vert(2,plon)      ! latitude grid vertices
    real(r8) :: ew_vert(2,plon)      ! longitude grid vertices
    real(r8) :: del_theta            ! difference in latitude at a grid point
    real(r8) :: del_phi              ! difference in longitude at a grid point
    real(r8), parameter :: degtorad=pi/180_r8
!-----------------------------------------------------------------------
    if(present(clon_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1, plon
                n = n + 1
                clon_d_out(n) = clon(i,j)
             end do
          end do
       else if(size == plon) then
          clon_d_out(:) = clon(:,1)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
    if(present(clat_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1, plon
                n = n + 1
                clat_d_out(n) = clat(j)
             end do
          end do
       else if(size == plat) then
          clat_d_out(:) = clat(:)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
    if ( ( present(wght_d_out) ) ) then

       if(size==plat) then
          wght_d_out(:) = (0.5_r8*w(:)/plon)* (4.0_r8*pi)
       else if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1, plon
                n = n + 1
                wght_d_out(n) = ( 0.5_r8*w(j)/plon ) * (4.0_r8*pi)
             end do
          end do
       end if
    end if
    if ( present(area_d_out) ) then
       if(size < ngcols_d) then
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
       n = 0
       do j = 1,plat

          ! First, determine vertices of each grid point.
          ! Verticies are ordered as follows:
          ! ns_vert: 1=lower left, 2 = upper left
          ! ew_vert: 1=lower left, 2 = lower right

          ! Latitude vertices
          ns_vert(:,:) = spval
          if (j .eq. 1) then
             ns_vert(1,:plon) = -90.0_r8
          else
             ns_vert(1,:plon) = (latdeg(j) + latdeg(j-1) )*0.5_r8
          end if

          if (j .eq. plat) then
             ns_vert(2,:plon) =  90.0_r8
          else
             ns_vert(2,:plon) = (latdeg(j) + latdeg(j+1) )*0.5_r8
          end if

          ! Longitude vertices
          ew_vert(:,:)       = spval
          ew_vert(1,1)       = (londeg(1,j) - 360.0_r8 + londeg(plon,j))*0.5_r8
          ew_vert(1,2:plon)  = (londeg(1:plon-1,j)+ londeg(2:plon,j))*0.5_r8
          ew_vert(2,:plon-1) = ew_vert(1,2:plon)
          ew_vert(2,plon)    = (londeg(plon,j) + (360.0_r8 + londeg(1,j)))*0.5_r8

          do i = 1, plon
             n = n + 1
             del_phi = sin( ns_vert(2,i)*degtorad ) - sin( ns_vert(1,i)*degtorad )
             del_theta = ( ew_vert(2,i) - ew_vert(1,i) )*degtorad
             area_d_out(n) = del_theta*del_phi
          end do

       end do
     end if
    if(present(lon_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1, plon
                n = n + 1
                lon_d_out(n) = londeg(i,j)
             end do
          end do
       else if(size == plon) then
          lon_d_out(:) = londeg(:,1)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
    if(present(lat_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1, plon
                n = n + 1
                lat_d_out(n) = latdeg(j)
             end do
          end do
       else if(size == plat) then
          lat_d_out(:) = latdeg(:)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
!
    return
  end subroutine get_horiz_grid_d


!#######################################################################
   function get_dyn_grid_parm_real2d(name) result(rval)
     use commap, only : londeg, clon
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:,:)

     if(name.eq.'clon') then
        rval => clon
     else if(name.eq.'londeg') then
        rval => londeg
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real2d

!#######################################################################
   function get_dyn_grid_parm_real1d(name) result(rval)
     use commap, only : latdeg, clat, w
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:)

     if(name.eq.'clat') then
        rval => clat
     else if(name.eq.'latdeg') then
        rval => latdeg
     else if(name.eq.'w') then
        rval => w
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real1d




   integer function get_dyn_grid_parm(name) result(ival)
     use pmgrid, only : beglat, endlat, plat, plon, plev, plevp
     character(len=*), intent(in) :: name

     if(name.eq.'beglat' .or. name .eq. 'beglatxy') then
        ival = beglat
     else if(name.eq.'endlat' .or. name .eq. 'endlatxy') then
        ival = endlat
     else if(name.eq.'plat') then
        ival = plat
     else if(name.eq.'plon' .or. name .eq. 'endlonxy') then
        ival = plon
     else if(name.eq.'plev') then
        ival = plev
     else if(name.eq.'plevp') then
        ival = plevp
     else if(name .eq. 'beglonxy') then
	ival = 1
     else
        ival = -1
     end if


   end function get_dyn_grid_parm

!#######################################################################

!-------------------------------------------------------------------------------
! This returns the lat/lon information (and corresponding MPI task numbers (owners))
! of the global model grid columns nearest to the input satellite coordinate (lat,lon)
!-------------------------------------------------------------------------------
subroutine dyn_grid_find_gcols( lat, lon, nclosest, owners, indx, jndx, rlat, rlon, idyn_dists )
  use spmd_utils,     only: iam
  use shr_const_mod,  only: SHR_CONST_PI, SHR_CONST_REARTH
  use pmgrid,         only: plon, plat

  real(r8), intent(in) :: lat
  real(r8), intent(in) :: lon
  integer, intent(in)  :: nclosest
  integer, intent(out) :: owners(nclosest)
  integer, intent(out) :: indx(nclosest)
  integer, intent(out) :: jndx(nclosest)

  real(r8),optional, intent(out) :: rlon(nclosest)
  real(r8),optional, intent(out) :: rlat(nclosest)
  real(r8),optional, intent(out) :: idyn_dists(nclosest)

  real(r8) :: dist            ! the distance (in radians**2 from lat, lon)
  real(r8) :: latr, lonr      ! lat, lon inputs converted to radians
  integer  :: ngcols
  integer  :: i, j

  integer :: blockid(1), bcid(1), lclblockid(1)

  real(r8), allocatable :: clat_d(:), clon_d(:), distmin(:)
  integer, allocatable :: igcol(:)
  real(r8), parameter :: rad2deg = 180._r8/SHR_CONST_PI

  latr = lat/rad2deg
  lonr = lon/rad2deg

  ngcols = plon*plat
  allocate( clat_d(1:ngcols) )
  allocate( clon_d(1:ngcols) )
  allocate( igcol(nclosest) )
  allocate( distmin(nclosest) )

  call get_horiz_grid_d(ngcols, clat_d_out=clat_d, clon_d_out=clon_d)

  igcol(:)    = -999
  distmin(:) = 1.e10_r8

  do i = 1,ngcols

     ! Use the Spherical Law of Cosines to find the great-circle distance.
     dist = acos(sin(latr) * sin(clat_d(i)) + cos(latr) * cos(clat_d(i)) * cos(clon_d(i) - lonr)) * SHR_CONST_REARTH
     do j = nclosest, 1, -1
        if (dist < distmin(j)) then

           if (j < nclosest) then
              distmin(j+1) = distmin(j)
              igcol(j+1)    = igcol(j)
           end if

           distmin(j) = dist
           igcol(j)    = i
        else
           exit
        end if
     end do

  end do

  do i = 1,nclosest

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))

     if ( iam==owners(i) ) then
        ! get global lat and lon coordinate indices from global column index
        ! -- plon is global number of longitude grid points
        jndx(i) = (igcol(i)-1)/plon + 1
        indx(i) = igcol(i) - (jndx(i)-1)*plon
     else
        jndx(i) = -1
        indx(i) = -1
     end if

     if ( present(rlat) ) rlat(i) = clat_d(igcol(i)) * rad2deg
     if ( present(rlon) ) rlon(i) = clon_d(igcol(i)) * rad2deg

     if (present(idyn_dists)) then
        idyn_dists(i) = distmin(i)
     end if

  end do

  deallocate( clat_d )
  deallocate( clon_d )
  deallocate( igcol )
  deallocate( distmin )

end subroutine dyn_grid_find_gcols

!#######################################################################
subroutine dyn_grid_get_colndx( igcol, nclosest, owners, indx, jndx )
  use spmd_utils, only: iam
  use pmgrid,     only: plon

  integer, intent(in)  :: nclosest
  integer, intent(in)  :: igcol(nclosest)
  integer, intent(out) :: owners(nclosest)
  integer, intent(out) :: indx(nclosest)
  integer, intent(out) :: jndx(nclosest)

  integer  :: i
  integer :: blockid(1), bcid(1), lclblockid(1)

  do i = 1,nclosest

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))

     if ( iam==owners(i) ) then
        ! get global lat and lon coordinate indices from global column index
        ! -- plon is global number of longitude grid points
        jndx(i) = (igcol(i)-1)/plon + 1
        indx(i) = igcol(i) - (jndx(i)-1)*plon
     else
        jndx(i) = -1
        indx(i) = -1
     endif

  end do

end subroutine dyn_grid_get_colndx
!#######################################################################

! this returns coordinates of a latitude slice of the block corresponding
! to latitude index latndx

subroutine dyn_grid_get_elem_coords( latndx, rlon, rlat, cdex )
  use commap, only : clat, clon
  use pmgrid, only : plon

  integer, intent(in) :: latndx ! lat  index

  real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the latndx slice
  real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the latndx slice
  integer, optional, intent(out) :: cdex(:) ! global column index

  integer :: i,ii,j

  if (present(cdex)) cdex(:) = -1
  if (present(rlat)) rlat(:) = -999._r8
  if (present(rlon)) rlon(:) = -999._r8

  j = latndx
  ii=0
  do i = 1,plon
     ii = ii+1
     if (present(cdex)) cdex(ii) = i + (j-1)*plon
     if (present(rlat)) rlat(ii) = clat(j)
     if (present(rlon)) rlon(ii) = clon(i,1)
  end do

end subroutine dyn_grid_get_elem_coords

!#######################################################################

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
  use cam_grid_support, only: max_hcoordname_len

  ! Dummy arguments
  character(len=max_hcoordname_len),          intent(out) :: gridname
  character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)

  gridname = 'gauss_grid'
  allocate(grid_attribute_names(4))
  grid_attribute_names(1) = 'gw'
  grid_attribute_names(2) = 'ntrm'
  grid_attribute_names(3) = 'ntrn'
  grid_attribute_names(4) = 'ntrk'

end subroutine physgrid_copy_attributes_d

!========================================================================================
! Private Methods
!========================================================================================


subroutine trunc()
!-----------------------------------------------------------------------
!
! Purpose:
! Check consistency of truncation parameters and evaluate pointers
! and displacements for spectral arrays
!
! Method:
!
! Author:
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!-----------------------------------------------------------------------

   use pspect,         only: ptrm, ptrn, ptrk, pmmax
   use comspe,         only: nstart, nlen, locm, lnstart

!---------------------------Local variables-----------------------------
!
   integer m              ! loop index
!
!-----------------------------------------------------------------------
!
! trunc first evaluates truncation parameters for a general pentagonal
! truncation for which the following parameter relationships are true
!
! 0 .le. |m| .le. ptrm
!
! |m| .le. n .le. |m|+ptrn for |m| .le. ptrk-ptrn
!
! |m| .le. n .le. ptrk     for (ptrk-ptrn) .le. |m| .le. ptrm
!
! Most commonly utilized truncations include:
!  1: triangular  truncation for which ptrk=ptrm=ptrn
!  2: rhomboidal  truncation for which ptrk=ptrm+ptrn
!  3: trapezoidal truncation for which ptrn=ptrk .gt. ptrm
!
! Simple sanity check
! It is necessary that ptrm .ge. ptrk-ptrn .ge. 0
!
   if (ptrm.lt.(ptrk-ptrn)) then
      call endrun ('TRUNC: Error in truncation parameters.  ntrm < (ptrk-ptrn)')
   end if
   if (ptrk.lt.ptrn) then
      call endrun ('TRUNC: Error in truncation parameters.  ptrk < ptrn')
   end if
!
! Evaluate pointers and displacement info based on truncation params
!
   nstart(1) = 0
   nlen(1) = ptrn + 1
   do m=2,pmmax
      nstart(m) = nstart(m-1) + nlen(m-1)
      nlen(m) = min0(ptrn+1,ptrk+2-m)
   end do
!
! Assign wavenumbers  and spectral offsets if not SPMD
!
#if ( ! defined SPMD )
   do m=1,pmmax
      locm(m,0) = m
      lnstart(m) = nstart(m)
   enddo
#endif

end subroutine trunc

!========================================================================================

subroutine define_cam_grids()
  use pspect,           only: ptrm, ptrn, ptrk
  use pmgrid,           only: beglat, endlat, plon, plat
  use commap,           only: londeg, latdeg, w
  use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
  use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register

  ! Local variables
  integer                      :: i, j, ind
  integer(iMap),       pointer :: grid_map(:,:)
  integer(iMap)                :: latmap(endlat - beglat + 1)
  type(horiz_coord_t), pointer :: lat_coord
  type(horiz_coord_t), pointer :: lon_coord
  real(r8),            pointer :: rattval(:)

  nullify(grid_map)
  nullify(lat_coord)
  nullify(lon_coord)
  nullify(rattval)

  ! Dynamics Grid
  ! Make grid and lat maps (need to do this because lat indices are distributed)
  ! Note that for this dycore, some pes may be inactive
  if(endlat >= beglat) then
    allocate(grid_map(4, (plon * (endlat - beglat + 1))))
    ind = 0
    do i = beglat, endlat
      do j = 1, plon
        ind = ind + 1
        grid_map(1, ind) = j
        grid_map(2, ind) = i
        grid_map(3, ind) = j
        grid_map(4, ind) = i
      end do
    end do
    ! Do we need a lat map?
    if ((beglat /= 1) .or. (endlat /= plat)) then
      do i = beglat, endlat
        latmap(i - beglat + 1) = i
      end do
    end if
  else
    allocate(grid_map(4, 0))
  end if

  ! Create the lat coordinate
  if ((beglat /= 1) .or. (endlat /= plat)) then
    lat_coord => horiz_coord_create('lat', '', plat, 'latitude',              &
         'degrees_north', beglat, endlat, latdeg(beglat:endlat), map=latmap)
  else
    lat_coord => horiz_coord_create('lat', '', plat, 'latitude',              &
         'degrees_north', beglat, endlat, latdeg(beglat:endlat))
  end if

  ! Create the lon coordinate
  lon_coord => horiz_coord_create('lon', '', plon, 'longitude',               &
       'degrees_east', 1, plon, londeg(1:plon, 1))

  call cam_grid_register('gauss_grid', dyn_decomp, lat_coord, lon_coord,      &
       grid_map, unstruct=.false.)

  allocate(rattval(size(w)))
  rattval = w
  call cam_grid_attribute_register('gauss_grid', 'gw', 'gauss weights', 'lat', rattval)
  nullify(rattval) ! belongs to attribute

  ! Scalar variable 'attributes'
  call cam_grid_attribute_register('gauss_grid', 'ntrm',                      &
      'spectral truncation parameter M', ptrm)
  call cam_grid_attribute_register('gauss_grid', 'ntrn',                      &
      'spectral truncation parameter N', ptrn)
  call cam_grid_attribute_register('gauss_grid', 'ntrk',                      &
      'spectral truncation parameter K', ptrk)
  ! These belong to the grid now
  nullify(grid_map)
  nullify(lat_coord)
  nullify(lon_coord)

end subroutine define_cam_grids

!========================================================================================

end module dyn_grid
