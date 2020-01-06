module dyn_grid
!-----------------------------------------------------------------------
!
! Define dynamics computational grid and decomposition.
!
! Original code: John Drake and Patrick Worley
!
!-----------------------------------------------------------------------

use shr_kind_mod,       only: r8=>shr_kind_r8
use pmgrid,             only: plon, plat, plev, plevp, splon, spmd_on
use commap,             only: w, w_staggered, clat, clat_staggered, clon, &
                              latdeg, londeg, latdeg_st, londeg_st
use constituents,       only: pcnst
use physconst,          only: pi, rearth, omega, spval
use spmd_utils,         only: iam
use spmd_dyn,           only: spmdinit_dyn, proc, lonrangexy, latrangexy
use time_manager,       only: get_step_size

use pio,                only: file_desc_t
use cam_initfiles,      only: initial_file_get_id

use dynamics_vars,      only: t_fvdycore_state, t_fvdycore_grid, grid_vars_init
use dyn_internal_state, only: get_dyn_state

use cam_abortutils,     only: endrun
use cam_logfile,        only: iulog

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

! The FV dynamics grids
integer, parameter, public :: dyn_decomp         = 101
integer, parameter, public :: dyn_stagger_decomp = 102 !Backward compatibility
integer, parameter, public :: dyn_ustag_decomp   = 102
integer, parameter, public :: dyn_vstag_decomp   = 103
integer, parameter, public :: dyn_zonal_decomp   = 104

integer, parameter, public :: ptimelevels = 2  ! number of time levels in the dycore

integer :: ngcols_d = 0     ! number of dynamics columns

type(t_fvdycore_grid), pointer :: grid

!========================================================================================
contains
!========================================================================================

subroutine dyn_grid_init()

   ! Initialize FV grid, decomposition, and PILGRIM communications

   use hycoef,         only: hycoef_init, hyai, hybi, hypi, hypm, nprlev
   use ref_pres,       only: ref_pres_init

   ! Local variables
   type(t_fvdycore_state), pointer :: state

   type(file_desc_t), pointer :: fh_ini

   integer :: i, k, lat

   real(r8) :: dt
   real(r8) :: dp
   real(r8) :: sum

   character(len=*), parameter :: sub='dyn_grid_init'
   !-----------------------------------------------------------------------

   ! Assign pointer to FV internal state object
   state => get_dyn_state()
   ! Assign pointer to grid object stored in state object
   grid => state%grid

   ! Get file handle for initial file and first consistency check
   fh_ini => initial_file_get_id()

   ! Set grid size parameters
   grid%im       = plon
   grid%jm       = plat
   grid%km       = plev
   grid%kmax     = plev + 1
   grid%nq       = pcnst
   grid%ntotq    = pcnst

   ! Initialize hybrid coordinate arrays
   call hycoef_init(fh_ini)

   ! Initialize reference pressures
   call ref_pres_init(hypi, hypm, nprlev)

   ! Hybrid coordinate info for FV grid object
   allocate(grid%ak(plev+1), grid%bk(plev+1))
   grid%ks = plev
   do k = 1, plev+1
      grid%ak(k) = hyai(k) * 1.e5_r8
      grid%bk(k) = hybi(k)
      if (grid%bk(k) == 0._r8) grid%ks = k-1
   end do
   grid%ptop = grid%ak(1)
   grid%pint = grid%ak(grid%ks+1)

   ! Initialize the grid decomposition and PILGRIM communications
   call spmdinit_dyn(state%jord, grid)

   ! Initialize FV specific grid object variables
   dt = get_step_size()
   call grid_vars_init(pi, rearth, omega, dt, state%fft_flt, &
                       state%am_geom_crrct, grid)

   ! initialize commap variables

   ! latitudes for cell centers
   dp = 180._r8/(plat-1)
   do lat = 1, plat
      latdeg(lat) = -90._r8 + (lat-1)*dp
      clat(lat) = latdeg(lat)*pi/180._r8
   end do

   ! latitudes for the staggered grid
   do lat = 1, plat-1
      clat_staggered(lat) = (clat(lat) + clat(lat+1)) / 2._r8
      latdeg_st     (lat) = clat_staggered(lat)*180._r8/pi
   end do

   ! Weights are defined as cos(phi)*(delta-phi)
   ! For a sanity check, the sum of w across all lats should be 2.
   do lat = 2, plat-1
      w(lat) = sin(clat_staggered(lat)) - sin(clat_staggered(lat-1))
   end do
   w(1) = sin(clat_staggered(1)) + 1._r8
   w(plat) = w(1)

   sum = 0._r8
   do lat=1,plat
      sum = sum + w(lat)
   end do
   if (abs(sum - 2._r8) > 1.e-8_r8) then
      write(iulog,*) sub//': ERROR: weights do not sum to 2. sum=', sum
      call endrun(sub//': ERROR: weights do not sum to 2.')
   end if

   dp = pi / real(plat-1,r8)
   do lat = 1, plat-1
      w_staggered(lat) = sin(clat(lat+1)) - sin(clat(lat))
   end do

   sum = 0._r8
   do lat = 1, plat-1
      sum = sum + w_staggered(lat)
   end do

   if (abs(sum - 2._r8) > 1.e-8_r8) then
      write(iulog,*) sub//': ERROR: staggered weights do not sum to 2. sum=', sum
      call endrun(sub//': ERROR: staggered weights do not sum to 2.')
   end if

   ! longitudes for cell centers
   do lat = 1, plat
      do i = 1, plon
         londeg(i,lat) = (i-1)*360._r8/plon
         clon(i,lat)   = (i-1)*2._r8*pi/plon
      end do
   end do

   ! longitudes for staggered grid
   do lat = 1, plat
      do i = 1, splon
         londeg_st(i,lat) = (i-1.5_r8)*360._r8/splon
      end do
   end do

   ! Define the CAM grids
   call define_cam_grids()

end subroutine dyn_grid_init

!========================================================================================

subroutine get_block_bounds_d(block_first, block_last)

   ! Return first and last indices used in global block ordering

   ! Arguments
   integer, intent(out) :: block_first  ! first (global) index used for blocks
   integer, intent(out) :: block_last   ! last (global) index used for blocks
   !---------------------------------------------------------------------------

   block_first = 1
   if (spmd_on .eq. 1) then
      ! Assume 1 block per subdomain
      block_last  = grid%nprxy_x*grid%nprxy_y
   else
      ! latitude slice block
      block_last  = plat
   end if

end subroutine get_block_bounds_d

!========================================================================================

subroutine get_block_gcol_d(blockid, size, cdex)

   ! Return list of dynamics column indices in given block

   ! Arguments
   integer, intent(in) :: blockid      ! global block id
   integer, intent(in) :: size         ! array size

   integer, intent(out):: cdex(size)   ! global column indices

   ! Local workspace
   integer :: i,j                        ! block coordinates
   integer :: blksiz                     ! block size
   integer :: k,l                        ! loop indices
   integer :: n                          ! column index
   character(len=*), parameter :: sub='get_block_gcol_d'
   !---------------------------------------------------------------------------

   if (spmd_on .eq. 1) then
      j = (blockid-1) / grid%nprxy_x + 1
      i = blockid - (j-1) * grid%nprxy_x
#if ( defined SPMD )
      blksiz = (lonrangexy(2,i)-lonrangexy(1,i)+1) *       &
               (latrangexy(2,j)-latrangexy(1,j)+1)
      if (size < blksiz) then
         write(iulog,*) sub//': ERROR: array not large enough (', &
                      size,' < ',blksiz,' ) '
         call endrun(sub//': ERROR: array not large enough')
      else
         n = 0
         do k=latrangexy(1,j),latrangexy(2,j)
            do l=lonrangexy(1,i),lonrangexy(2,i)
               n = n + 1
               cdex(n) = l + (k-1)*plon
            end do
         end do
      end if
#endif
   else
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
   end if

end subroutine get_block_gcol_d

!========================================================================================

integer function get_block_gcol_cnt_d(blockid)

   ! Return number of dynamics columns in indicated block

   ! Arguments
   integer, intent(in) :: blockid  ! global block id

   ! Local workspace
   integer :: i, j
   !---------------------------------------------------------------------------

   if (spmd_on .eq. 1) then
      j = (blockid-1) / grid%nprxy_x + 1
      i = blockid - (j-1) * grid%nprxy_x
#if ( defined SPMD )
      get_block_gcol_cnt_d = (lonrangexy(2,i)-lonrangexy(1,i)+1) *       &
           (latrangexy(2,j)-latrangexy(1,j)+1)
#endif
   else
      get_block_gcol_cnt_d = plon
   end if

end function get_block_gcol_cnt_d

!========================================================================================

integer function get_block_lvl_cnt_d(blockid, bcid)

   ! Return number of levels in indicated column. If column
   ! includes surface fields, then it is defined to also
   ! include level 0.

   ! Arguments
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid     ! column index within block
   !-----------------------------------------------------------------------

   ! latitude slice block
   get_block_lvl_cnt_d = plev + 1

end function get_block_lvl_cnt_d

!========================================================================================

subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

   ! Return level indices in indicated column. If column
   ! includes surface fields, then it is defined to also
   ! include level 0.

   ! Arguments
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block
   integer, intent(in) :: lvlsiz   ! dimension of levels array

   integer, intent(out) :: levels(lvlsiz) ! levels indices for block

   ! Local workspace
   integer :: k                      ! loop index
   character(len=*), parameter :: sub='get_block_levels_d'
   !-----------------------------------------------------------------------

   ! latitude slice block
   if (lvlsiz < plev + 1) then
      write(iulog,*) sub//': ERROR: levels array not large enough (', &
                          lvlsiz,' < ',plev + 1,' ) '
      call endrun(sub//': ERROR: levels array not large enough')
   else
      do k = 0, plev
         levels(k+1) = k
      end do
      do k = plev+2, lvlsiz
         levels(k) = -1
      end do
   end if

end subroutine get_block_levels_d

!========================================================================================

subroutine get_gcol_block_d(gcol, cnt, blockid, bcid, localblockid)

   ! Return global block index and local column index
   ! for global column index

   ! Arguments
   integer, intent(in) :: gcol     ! global column index
   integer, intent(in) :: cnt      ! size of blockid and bcid arrays

   integer, intent(out) :: blockid(cnt) ! block index
   integer, intent(out) :: bcid(cnt)    ! column index within block
   integer, intent(out), optional :: localblockid(cnt)

   ! Local workspace
   integer :: i, j, ii, jj                ! loop indices
   integer :: glon, glat                  ! global longitude and latitude indices
   integer :: ddlon                       ! number of longitudes in block
   character(len=*), parameter :: sub='get_gcol_block_d'
   !---------------------------------------------------------------------------

   ! lon/lat block
   if (cnt < 1) then
      write(iulog,*) sub//': ERROR: arrays not large enough (', cnt,' < ',1,' )'
      call endrun(sub//': ERROR: arrays not large enough')
   else
      if (spmd_on == 1) then
         ! Determine global latitude and longitude coordinate indices from
         ! global column index
         glat = (gcol-1)/plon + 1
         glon = gcol - ((glat-1)*plon)

         ! Determine block coordinates (ii,jj), where ii ranges from 1 to
         ! nprxy_x and jj ranges from 1 to nprxy_y.
#if ( defined SPMD )
         ii = 0
         do i = 1, grid%nprxy_x
            if (lonrangexy(1,i) <= glon .and. glon <= lonrangexy(2,i)) ii=i
         end do
         jj = 0
         do j = 1, grid%nprxy_y
            if (latrangexy(1,j) <= glat .and. glat <= latrangexy(2,j)) jj=j
         end do
         if (ii == 0 .or. jj == 0) then
            write(iulog,*) sub//': ERROR: could not find block indices for (', &
                           glon,',',glat,' ) '
            call endrun(sub//': ERROR: could not find block indices')
         end if

         ! Global block index
         blockid(1) = (jj-1)*grid%nprxy_x+ii

         ! Local coordinates in block
         j = glat-latrangexy(1,jj)+1
         i = glon-lonrangexy(1,ii)+1
         ddlon = lonrangexy(2,ii)-lonrangexy(1,ii)+1

         ! Local column index in block
         bcid(1) = (j-1)*ddlon+i
#endif
      else
         glat = (gcol-1)/plon + 1
         glon = gcol - ((glat-1)*plon)

         blockid(1) = glat
         bcid(1)    = glon
      end if

      do j=2,cnt
         blockid(j) = -1
         bcid(j)    = -1
      end do

   end if

end subroutine get_gcol_block_d

!========================================================================================

integer function get_gcol_block_cnt_d(gcol)

   ! Return number of blocks contain data for the vertical column
   ! with the given global column index

   ! Arguments
   integer, intent(in) :: gcol     ! global column index
   !---------------------------------------------------------------------------

   !  lon/lat block
   get_gcol_block_cnt_d = 1

end function get_gcol_block_cnt_d

!========================================================================================

integer function get_block_owner_d(blockid)

   ! Return id of processor that "owns" the indicated block

   ! Arguments
   integer, intent(in) :: blockid  ! global block id
   !---------------------------------------------------------------------------

   ! latitude slice block
#if (defined SPMD)
   if (spmd_on .eq. 1) then
      get_block_owner_d = blockid - 1
   else
      get_block_owner_d = proc(blockid)
   end if
#else
   get_block_owner_d = 0
#endif

end function get_block_owner_d

!========================================================================================

subroutine get_horiz_grid_dim_d(hdim1_d, hdim2_d)

   ! Returns declared horizontal dimensions of computational grid.
   ! Note that global column ordering is assumed to be compatible
   ! with the first dimension major ordering of the 2D array.

   ! Arguments
   integer, intent(out) :: hdim1_d           ! first horizontal dimension
   integer, intent(out) :: hdim2_d           ! second horizontal dimension
   !---------------------------------------------------------------------------

   if (ngcols_d == 0) then
      ngcols_d = plat*plon
   end if
   hdim1_d = plon
   hdim2_d = plat

end subroutine get_horiz_grid_dim_d

!========================================================================================

subroutine get_horiz_grid_d(size, clat_d_out, clon_d_out, area_d_out, wght_d_out, &
                            lat_d_out, lon_d_out)

   ! Return latitude and longitude (in radians), column surface
   ! area (in radians squared) and surface integration weights
   ! for global column indices that will be passed to/from physics

   ! Arguments
   integer, intent(in)   :: size             ! array sizes

   real(r8), intent(out), optional :: clat_d_out(size) ! column latitudes
   real(r8), intent(out), optional :: clon_d_out(size) ! column longitudes

   real(r8), intent(out), optional :: area_d_out(size) ! column surface
                                                       !  area
   real(r8), intent(out), optional :: wght_d_out(size) ! column integration
                                                       !  weight
   real(r8), intent(out), optional :: lat_d_out(size)  ! column deg latitudes
   real(r8), intent(out), optional :: lon_d_out(size)  ! column deg longitudes

   ! Local workspace
   integer  :: i, j            ! loop indices
   integer  :: n               ! column index
   real(r8) :: ns_vert(2,plon) ! latitude grid vertices
   real(r8) :: ew_vert(2,plon) ! longitude grid vertices
   real(r8) :: del_theta       ! difference in latitude at a grid point
   real(r8) :: del_phi         ! difference in longitude at a grid point
   real(r8), parameter :: degtorad=pi/180.0_r8 ! convert degrees to radians
   character(len=128)  :: errormsg
   character(len=*), parameter :: sub='get_horiz_grid_d'
   !----------------------------------------------------------------------------

   if (present(clon_d_out)) then
      if (size == ngcols_d) then
         n = 0
         do j = 1, plat
            do i = 1, plon
               n = n + 1
               clon_d_out(n) = clon(i,j)
            end do
         end do
      else if(size == plon) then
         clon_d_out(:) = clon(:,1)
      else
         write(errormsg, '(a,4(i0,a))')'clon_d_out array size incorrect (',    &
               size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plon,') '
         call endrun(sub//': ERROR: '//errormsg)
      end if
   end if

   if (present(clat_d_out)) then
      if (size == ngcols_d) then
         n = 0
         do j = 1, plat
            do i = 1, plon
               n = n + 1
               clat_d_out(n) = clat(j)
            end do
         end do
      else if (size == plat) then
         clat_d_out(:) = clat(:)
      else
         write(errormsg, '(a,4(i0,a))')'clat_d_out array size incorrect (',    &
               size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plat,') '
         call endrun(sub//': ERROR: '//errormsg)
      end if
   end if

   if (size==plat .and. present(wght_d_out)) then
      wght_d_out(:) = w(:)
      return
   end if

   if ( ( present(area_d_out) ) .or. ( present(wght_d_out) ) ) then
      if ((size < ngcols_d) .and. present(area_d_out)) then
         write(errormsg, '(a,2(i0,a))')'area_d_out array size incorrect (',  &
             size, ' /= ', ngcols_d, ') '
         call endrun(sub//': ERROR: '//errormsg)
      else if ((size < ngcols_d) .and. present(area_d_out)) then
         write(errormsg, '(a,2(i0,a))')'wght_d_out array size incorrect (',  &
             size, ' /= ', ngcols_d, ') '
         call endrun(sub//': ERROR: '//errormsg)
      end if

      n = 0
      do j = 1, plat

         ! First, determine vertices of each grid point.
         ! Verticies are ordered as follows:
         ! ns_vert: 1=lower left, 2 = upper left
         ! ew_vert: 1=lower left, 2 = lower right

         ! Latitude vertices
         ns_vert(:,:) = spval
         if (j .eq. 1) then
            ns_vert(1,:plon)    = -90._r8 + (latdeg(1) - latdeg(2))*0.5_r8
         else
            ns_vert(1,:plon)    = (latdeg(j) + latdeg(j-1) )*0.5_r8
         end if

         if (j .eq. plat) then
            ns_vert(2,:plon) =  90._r8 + (latdeg(plat) - latdeg(plat-1))*0.5_r8
         else
            ns_vert(2,:plon)    = (latdeg(j) + latdeg(j+1) )*0.5_r8
         end if

         ! Longitude vertices
         ew_vert(:,:) = spval
         ew_vert(1,1)          = (londeg(1,j) - 360._r8 + londeg(plon,j))*0.5_r8
         ew_vert(1,2:plon)  = (londeg(1:plon-1,j)+ londeg(2:plon,j))*0.5_r8
         ew_vert(2,:plon-1) = ew_vert(1,2:plon)
         ew_vert(2,plon)    = (londeg(plon,j) + (360._r8 + londeg(1,j)))*0.5_r8

         do i = 1,plon
            n = n + 1

            if (j .eq. 1) then
               del_phi = -sin( latdeg(j)*degtorad )    + sin( ns_vert(2,i)*degtorad )
            else if (j .eq. plat) then
               del_phi =  sin( latdeg(j)*degtorad )    - sin( ns_vert(1,i)*degtorad )
            else
               del_phi =  sin( ns_vert(2,i)*degtorad ) - sin( ns_vert(1,i)*degtorad )
            end if

            del_theta = ( ew_vert(2,i) - ew_vert(1,i) )*degtorad

            if (present(area_d_out)) area_d_out(n) = del_theta*del_phi
            if (present(wght_d_out)) wght_d_out(n) = del_theta*del_phi
         end do
      end do
   end if

   if (present(lon_d_out)) then
      if (size == ngcols_d) then
         n = 0
         do j = 1, plat
            do i = 1, plon
               n = n + 1
               lon_d_out(n) = londeg(i,j)
            end do
         end do
      else if(size == plon) then
         lon_d_out(:) = londeg(:,1)
      else
         write(errormsg, '(a,4(i0,a))')'lon_d_out array size incorrect (',    &
             size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plon,') '
         call endrun(sub//': ERROR: '//errormsg)
      end if
   end if

   if (present(lat_d_out)) then
      if (size == ngcols_d) then
         n = 0
         do j = 1, plat
            do i = 1, plon
               n = n +  1
               lat_d_out(n) = latdeg(j)
            end do
         end do
      else if (size == plat) then
         lat_d_out(:) = latdeg(:)
      else
         write(errormsg, '(a,4(i0,a))')'lat_d_out array size incorrect (',    &
             size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plat,') '
         call endrun(sub//': ERROR: '//errormsg)
      end if
   end if

end subroutine get_horiz_grid_d

!========================================================================================

function get_dyn_grid_parm_real2d(name) result(rval)

   character(len=*), intent(in) :: name
   real(r8), pointer :: rval(:,:)

   if(name.eq.'clon') then
      rval => clon
   else if(name.eq.'londeg') then
      rval => londeg
   else if(name.eq.'londeg_st') then
      rval => londeg_st
   else
      nullify(rval)
   end if
end function get_dyn_grid_parm_real2d

!========================================================================================

function get_dyn_grid_parm_real1d(name) result(rval)

   character(len=*), intent(in) :: name
   real(r8), pointer :: rval(:)

   if(name.eq.'clat') then
      rval => clat
   else if(name.eq.'latdeg') then
      rval => latdeg
   else if(name.eq.'latdeg_st') then
      rval => latdeg_st
   else if(name.eq.'clatdeg_staggered') then
      rval => latdeg_st
   else if(name.eq.'w') then
      rval => w
   else if(name.eq.'w_staggered') then
      rval => w_staggered
   else
      nullify(rval)
   end if
end function get_dyn_grid_parm_real1d

!========================================================================================

integer function get_dyn_grid_parm(name) result(ival)

   character(len=*), intent(in) :: name

   if (name.eq.'splon') then
      ival = splon
   else if (name.eq.'beglonxy') then
      ival = grid%ifirstxy
   else if (name.eq.'endlonxy') then
      ival = grid%ilastxy
   else if (name.eq.'beglatxy') then
      ival = grid%jfirstxy
   else if (name.eq.'endlatxy') then
      ival = grid%jlastxy
   else if (name.eq.'plat') then
      ival = plat
   else if (name.eq.'plon') then
      ival = plon
   else if (name.eq.'plev') then
      ival = plev
   else if (name.eq.'plevp') then
      ival = plevp
   else
      ival = -1
   end if

end function get_dyn_grid_parm

!========================================================================================

subroutine dyn_grid_find_gcols(lat, lon, nclosest, owners, indx, &
                               jndx, rlat, rlon, idyn_dists)

   ! Return the lat/lon information (and corresponding MPI task numbers (owners))
   ! of the global model grid columns nearest to the input coordinate (lat,lon)

   ! arguments
   real(r8), intent(in) :: lat
   real(r8), intent(in) :: lon
   integer, intent(in)  :: nclosest
   integer, intent(out) :: owners(nclosest)
   integer, intent(out) :: indx(nclosest)
   integer, intent(out) :: jndx(nclosest)

   real(r8),optional, intent(out) :: rlon(nclosest)
   real(r8),optional, intent(out) :: rlat(nclosest)
   real(r8),optional, intent(out) :: idyn_dists(nclosest)

   ! local variables
   real(r8) :: dist            ! the distance (in radians**2 from lat, lon)
   real(r8) :: latr, lonr      ! lat, lon inputs converted to radians
   integer  :: ngcols
   integer  :: i, j

   integer :: blockid(1), bcid(1), lclblockid(1)

   real(r8), allocatable :: clat_d(:), clon_d(:), distmin(:)
   integer, allocatable :: igcol(:)
   real(r8), parameter :: rad2deg = 180._r8/pi
   !----------------------------------------------------------------------------

   latr = lat/rad2deg
   lonr = lon/rad2deg

   ngcols = plon*plat
   allocate( clat_d(1:ngcols) )
   allocate( clon_d(1:ngcols) )
   allocate( igcol(nclosest) )
   allocate( distmin(nclosest) )

   call get_horiz_grid_d(ngcols, clat_d_out=clat_d, clon_d_out=clon_d)

   igcol(:)   = -999
   distmin(:) = 1.e10_r8

   do i = 1, ngcols

      ! Use the Spherical Law of Cosines to find the great-circle distance.
      dist = acos(sin(latr) * sin(clat_d(i)) + cos(latr) * cos(clat_d(i)) * &
             cos(clon_d(i) - lonr)) * rearth
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
      endif

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

!========================================================================================

subroutine dyn_grid_get_colndx(igcol, nclosest, owners, indx, jndx)

   ! arguments
   integer, intent(in)  :: nclosest
   integer, intent(in)  :: igcol(nclosest)
   integer, intent(out) :: owners(nclosest)
   integer, intent(out) :: indx(nclosest)
   integer, intent(out) :: jndx(nclosest)

   integer  :: i
   integer :: blockid(1), bcid(1), lclblockid(1)

   do i = 1,nclosest

      call  get_gcol_block_d(igcol(i), 1, blockid, bcid, lclblockid)
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

   end do

end subroutine dyn_grid_get_colndx

!========================================================================================

subroutine dyn_grid_get_elem_coords( latndx, rlon, rlat, cdex )

   ! return coordinates of a latitude slice of the block corresponding
   ! to latitude index latndx

   ! arguments
   integer, intent(in) :: latndx ! lat  index

   real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the latndx slice
   real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the latndx slice
   integer, optional, intent(out) :: cdex(:) ! global column index

   integer :: i, ii, j
   !----------------------------------------------------------------------------

   if (present(cdex)) cdex(:) = -1
   if (present(rlat)) rlat(:) = -999._r8
   if (present(rlon)) rlon(:) = -999._r8

   j  = latndx
   ii = 0
   do i = grid%ifirstxy, grid%ilastxy
      ii = ii+1
      if (present(cdex)) cdex(ii) = i + (j-1)*plon
      if (present(rlat)) rlat(ii) = clat(j)
      if (present(rlon)) rlon(ii) = clon(i,1)
   end do

end subroutine dyn_grid_get_elem_coords

!========================================================================================

subroutine define_cam_grids()

  use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
  use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
  use cam_grid_support, only: cam_grid_attribute_copy

  integer                      :: i, j, ind
  integer                      :: beglonxy, endlonxy
  integer                      :: beglatxy, endlatxy
  type(horiz_coord_t), pointer :: lat_coord
  type(horiz_coord_t), pointer :: lon_coord
  type(horiz_coord_t), pointer :: slat_coord
  type(horiz_coord_t), pointer :: slon_coord
  type(horiz_coord_t), pointer :: zlon_coord
  integer(iMap),   allocatable :: coord_map(:)
  integer(iMap),       pointer :: grid_map(:,:)
  real(r8)                     :: zlon_bnds(2,1)
  real(r8), allocatable        :: latvals(:)
  real(r8),            pointer :: rattval(:)
  logical                      :: is_lon_distributed
  !-----------------------------------------------------------------------------

  ! Note: not using get_horiz_grid_dim_d or get_horiz_grid_d since those
  !       are deprecated ('cause I said so' -- goldy)

  nullify(lat_coord)
  nullify(lon_coord)
  nullify(slat_coord)
  nullify(slon_coord)
  nullify(zlon_coord)
  nullify(grid_map)
  nullify(rattval)

  beglonxy = grid%ifirstxy
  endlonxy = grid%ilastxy
  beglatxy = grid%jfirstxy
  endlatxy = grid%jlastxy

  if (iam >= grid%npes_xy) then
    ! NB: On inactive PEs, beglonxy should be one and endlonxy should be zero
    if (beglonxy /= 1) then
      call endrun("DEFINE_CAM_GRIDS: ERROR: Bad beglonxy")
    end if
    if (endlonxy /= 0) then
      call endrun("DEFINE_CAM_GRIDS: ERROR: Bad endlonxy")
    end if
    ! NB: On inactive PEs, beglatxy should be one and endlatxy should be zero
    if (beglatxy /= 1) then
      call endrun("DEFINE_CAM_GRIDS: ERROR: Bad beglatxy")
    end if
    if (endlatxy /= 0) then
      call endrun("DEFINE_CAM_GRIDS: ERROR: Bad endlatxy")
    end if
  end if

  ! Figure out if lon and slon are distributed dimensions
  is_lon_distributed = (grid%nprxy_x > 1)

  ! Grid for cell centers
  ! Make a map
  allocate(grid_map(4, ((endlonxy - beglonxy + 1) * (endlatxy - beglatxy + 1))))
  ind = 0
  do i = beglatxy, endlatxy
    do j = beglonxy, endlonxy
      ind = ind + 1
      grid_map(1, ind) = j
      grid_map(2, ind) = i
      grid_map(3, ind) = j
      grid_map(4, ind) = i
    end do
  end do
  ! Cell-centered latitude coordinate
  allocate(coord_map(endlatxy - beglatxy + 1))
  if (endlatxy >= beglatxy) then
    if (beglonxy == 1) then
      coord_map = (/ (i, i = beglatxy, endlatxy) /)
    else
      coord_map = 0
    end if
  end if
  lat_coord => horiz_coord_create('lat', '', plat, 'latitude',                &
       'degrees_north', beglatxy, endlatxy, latdeg(beglatxy:endlatxy),        &
       map=coord_map)
  deallocate(coord_map)

  ! Cell-centered longitude coordinate
  if (is_lon_distributed) then
    allocate(coord_map(endlonxy - beglonxy + 1))
    if (endlonxy >= beglonxy) then
      if (beglatxy == 1) then
        coord_map = (/ (i, i = beglonxy, endlonxy) /)
      else
        coord_map = 0
      end if
    end if
    lon_coord => horiz_coord_create('lon', '', plon, 'longitude',             &
         'degrees_east', beglonxy, endlonxy, londeg(beglonxy:endlonxy,1),     &
         map=coord_map)
    deallocate(coord_map)
  else
    lon_coord => horiz_coord_create('lon', '', plon, 'longitude',             &
         'degrees_east', beglonxy, endlonxy, londeg(beglonxy:endlonxy,1))
  end if
  ! Cell-centered grid
  call cam_grid_register('fv_centers', dyn_decomp, lat_coord, lon_coord,      &
       grid_map, unstruct=.false.)
  allocate(rattval(size(w)))
  rattval = w
  call cam_grid_attribute_register('fv_centers', 'gw', 'latitude weights', 'lat', rattval)
  nullify(rattval)
  nullify(grid_map) ! Belongs to the grid

  ! Staggered grid for U_S
  ! Make a map
  allocate(grid_map(4, ((endlonxy - beglonxy + 1) * (endlatxy - beglatxy + 1))))
  ind = 0
  do i = beglatxy, endlatxy
    do j = beglonxy, endlonxy
      ind = ind + 1
      grid_map(1, ind) = j
      grid_map(2, ind) = i
      grid_map(3, ind) = j
      if ((i == beglatxy) .and. (beglatxy == 1)) then
        grid_map(4, ind) = 0
      else
        grid_map(4, ind) = i - 1
      end if
    end do
  end do

  ! Staggered latitudes 'skip' the first one so they are 'off by one'
  ! This means we always must have a coordinate map
  allocate(coord_map(endlatxy - beglatxy + 1))
  ! NB: coord_map(1) == 0 when beglat == 1, that element is not output
  do i = 1, size(coord_map)
    if (beglonxy == 1) then
      coord_map(i) = i + beglatxy - 2
    else
      coord_map(i) = 0
    end if
  end do
  if (iam .lt. grid%npes_xy) then
    allocate(latvals(beglatxy:endlatxy))
    if (beglatxy == 1) then
      latvals(1) = 0
      latvals(2:endlatxy) = latdeg_st(1:endlatxy-1)
    else
      i = beglatxy - 1 ! Stupid NAG 'error'
      latvals(beglatxy:endlatxy) = latdeg_st(i:endlatxy-1)
    end if
  else
    allocate(latvals(0))
  end if
  slat_coord => horiz_coord_create('slat', '', (plat - 1),                    &
       'staggered latitude', 'degrees_north', beglatxy, endlatxy, latvals,    &
       map=coord_map)
  deallocate(coord_map)
  deallocate(latvals)

  call cam_grid_register('fv_u_stagger', dyn_ustag_decomp, slat_coord,        &
       lon_coord, grid_map, unstruct=.false.)
  call cam_grid_attribute_register('fv_u_stagger', 'w_stag',                  &
       'staggered latitude weights', 'slat', w_staggered)
  nullify(grid_map) ! Belongs to the grid

  ! Staggered grid for V_S
  ! Make a map (need to do this because lat indices are distributed)
  allocate(grid_map(4, ((endlonxy - beglonxy + 1) * (endlatxy - beglatxy + 1))))
  ind = 0
  do i = beglatxy, endlatxy
    do j = beglonxy, endlonxy
      ind = ind + 1
      grid_map(1, ind) = j
      grid_map(2, ind) = i
      grid_map(3, ind) = j
      grid_map(4, ind) = i
    end do
  end do
  ! Staggered longitude coordinate
  if (is_lon_distributed) then
    allocate(coord_map(endlonxy - beglonxy + 1))
    if (endlonxy >= beglonxy) then
      if (beglatxy == 1) then
        coord_map = (/ (i, i = beglonxy, endlonxy) /)
      else
        coord_map = 0
      end if
    end if
    slon_coord => horiz_coord_create('slon', '', plon, 'staggered longitude', &
         'degrees_east', beglonxy, endlonxy, londeg_st(beglonxy:endlonxy,1),  &
         map=coord_map)
    deallocate(coord_map)
  else
    slon_coord => horiz_coord_create('slon', '', plon, 'staggered longitude', &
         'degrees_east', beglonxy, endlonxy, londeg_st(beglonxy:endlonxy,1))
  end if
  call cam_grid_register('fv_v_stagger', dyn_vstag_decomp, lat_coord,         &
       slon_coord, grid_map, unstruct=.false.)
  nullify(grid_map) ! Belongs to the grid

  ! Zonal mean grid
  ! Make a map
  allocate(grid_map(4, (endlatxy - beglatxy + 1)))
  ind = 0
  do i = beglatxy, endlatxy
    ind = ind + 1
    grid_map(1, ind) = 1
    grid_map(2, ind) = i
    grid_map(3, ind) = 1
    grid_map(4, ind) = i
  end do
  ! We need a special, size-one "longigude" coordinate
  ! NB: This is never a distributed coordinate so calc even on inactive PEs
  zlon_bnds(1,1) = minval(londeg)
  zlon_bnds(2,1) = maxval(londeg)
  allocate(latvals(1)) ! Really for a longitude
  latvals(1) = 0._r8
  zlon_coord => horiz_coord_create('zlon', '', 1, 'longitude',                &
       'degrees_east', 1, 1, latvals(1:1), bnds=zlon_bnds)
  deallocate(latvals)
  ! Zonal mean grid
  call cam_grid_register('fv_centers_zonal', dyn_zonal_decomp, lat_coord,     &
       zlon_coord, grid_map, unstruct=.false., zonal_grid=.true.)
  ! Make sure 'gw' attribute shows up even if all variables are zonal mean
  call cam_grid_attribute_copy('fv_centers', 'fv_centers_zonal', 'gw')
  nullify(grid_map) ! Belongs to the grid

end subroutine define_cam_grids

!========================================================================================

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
  use cam_grid_support, only: max_hcoordname_len

  ! Dummy arguments
  character(len=max_hcoordname_len),          intent(out) :: gridname
  character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)

  gridname = 'fv_centers'
  allocate(grid_attribute_names(1))
  grid_attribute_names(1) = 'gw'

end subroutine physgrid_copy_attributes_d

!========================================================================================

end module dyn_grid
