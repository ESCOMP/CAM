module dyn_grid
!-------------------------------------------------------------------------------
!
! Define SE computational grids on the dynamics decomposition.
!

! The grid used by the SE dynamics is called the GLL grid.  It is
! decomposed into elements which correspond to "blocks" in the
! physics/dynamics coupler terminology.  The columns in this grid are
! located at the Gauss-Lobatto-Legendre (GLL) quadrature points.  The GLL
! grid will also be used by the physics if the CSLAM advection is not used.
! If CSLAM is used for tracer advection then it uses an FVM grid and the
! physics will either use the same FVM grid or an FVM grid with a different
! number of equal area subcells.  The FVM grid used by the physics is
! referred to as the "physgrid".
!
! Module responsibilities:
!
! . Provide the physics/dynamics coupler (in module phys_grid) with data for the
!   physics grid on the dynamics decomposition.
!
! . Create CAM grid objects that are used by the I/O functionality to read
!   data from an unstructured grid format to the dynamics data structures, and
!   to write from the dynamics data structures to unstructured grid format.  The
!   global column ordering for the unstructured grid is determined by the SE dycore.
!
!-------------------------------------------------------------------------------

use shr_kind_mod,           only: r8 => shr_kind_r8, shr_kind_cl
use spmd_utils,             only: masterproc, iam, mpicom, mstrid=>masterprocid
use spmd_utils,             only: npes, mpi_integer, mpi_real8
use constituents,           only: pcnst
use physconst,              only: pi
use cam_initfiles,          only: initial_file_get_id
use cam_grid_support,       only: iMap
use physics_column_type,    only: physics_column_t

use cam_logfile,            only: iulog
use cam_abortutils,         only: endrun

use pio,                    only: file_desc_t

use dimensions_mod,         only: globaluniquecols, nelem, nelemd, nelemdmax
use dimensions_mod,         only: ne, np, npsq, fv_nphys, nlev, ntrac
use element_mod,            only: element_t
use fvm_control_volume_mod, only: fvm_struct
use hybvcoord_mod,          only: hvcoord_t
use prim_init,              only: prim_init1
use edge_mod,               only: initEdgeBuffer
use edgetype_mod,           only: EdgeBuffer_t
use time_mod,               only: TimeLevel_t
use dof_mod,                only: UniqueCoords, UniquePoints

implicit none
private
save

integer, parameter :: dyn_decomp = 101 ! The SE dynamics grid
integer, parameter :: fvm_decomp = 102 ! The FVM (CSLAM) grid
integer, parameter :: physgrid_d = 103 ! physics grid on dynamics decomp
integer, parameter :: ini_decomp = 104 ! alternate dynamics grid for reading initial file

character(len=3), protected :: ini_grid_name

! Name of horizontal grid dimension in initial file.
character(len=6), protected :: ini_grid_hdim_name = ''

integer, parameter :: ptimelevels = 2

type (TimeLevel_t)         :: TimeLevel     ! main time level struct (used by tracers)
type (hvcoord_t)           :: hvcoord
type(element_t),   pointer :: elem(:) => null()  ! local GLL elements for this task
type(fvm_struct),  pointer :: fvm(:) => null()   ! local FVM elements for this task

public :: dyn_decomp
public :: ini_grid_name
public :: ini_grid_hdim_name
public :: ptimelevels
public :: TimeLevel
public :: hvcoord
public :: elem
public :: fvm
public :: edgebuf

public :: dyn_grid_init              ! Initialize the dynamics grid
public :: get_dyn_grid_info          ! Return physics grid column information
public :: physgrid_copy_attributes_d ! Attributes to copy to physics grid

!!XXgoldyXX: v These need to be removed to complete the weak scaling transition.
public :: get_horiz_grid_d
public :: get_horiz_grid_dim_d
public :: get_dyn_grid_parm
public :: get_dyn_grid_parm_real1d
!!XXgoldyXX: ^ These need to be removed to complete the weak scaling transition.
public :: dyn_grid_get_elem_coords ! get coords of a specified block element

! Note: dyn_grid_get_colndx is not implemenented in SE
public :: dyn_grid_get_colndx ! get element block/column and MPI process indices

! Namelist variables controlling grid writing.
! Read in dyn_readnl from dyn_se_inparm group.
character(len=16),          public :: se_write_grid_file   = 'no'
character(len=shr_kind_cl), public :: se_grid_filename     = ''
logical,                    public :: se_write_gll_corners = .false.

type(physics_column_t), allocatable, target :: local_dyn_columns(:)

! number of global dynamics columns. Set by SE dycore init.
integer :: ngcols_d = 0
! number of global elements. Set by SE dycore init.
integer :: nelem_d = 0

real(r8), parameter :: rad2deg = 180.0_r8/pi

type(EdgeBuffer_t) :: edgebuf

!=============================================================================
contains
!=============================================================================

subroutine dyn_grid_init()

   ! Initialize SE grid, and decomposition.

   use hycoef,              only: hycoef_init, hypi, hypm, nprlev, &
                                  hyam, hybm, hyai, hybi, ps0
   use ref_pres,            only: ref_pres_init
   use spmd_utils,          only: MPI_MAX, MPI_INTEGER, mpicom
   use time_manager,        only: get_nstep, get_step_size
   use dp_mapping,          only: dp_init, dp_write
   use native_mapping,      only: do_native_mapping, create_native_mapping_files

   use parallel_mod,        only: par
   use hybrid_mod,          only: hybrid_t, init_loop_ranges, &
                                  get_loop_ranges, config_thread_region
   use control_mod,         only: qsplit, rsplit
   use time_mod,            only: tstep, nsplit
   use fvm_mod,             only: fvm_init2, fvm_init3, fvm_pg_init
   use dimensions_mod,      only: irecons_tracer
   use comp_gll_ctr_vol,    only: gll_grid_write
   use air_composition,     only: thermodynamic_active_species_num

   ! Local variables

   type(file_desc_t), pointer  :: fh_ini

   integer                     :: qsize_local
   integer                     :: k

   type(hybrid_t)              :: hybrid
   integer                     :: ierr
   integer                     :: dtime

   real(r8), allocatable       ::clat(:), clon(:), areaa(:)
   integer                     :: nets, nete

   character(len=*), parameter :: sub = 'dyn_grid_init'
   !----------------------------------------------------------------------------

   ! Get file handle for initial file and first consistency check
   fh_ini => initial_file_get_id()

   ! Initialize hybrid coordinate arrays
   call hycoef_init(fh_ini, psdry=.true.)

   hvcoord%hyam = hyam
   hvcoord%hyai = hyai
   hvcoord%hybm = hybm
   hvcoord%hybi = hybi
   hvcoord%ps0  = ps0
   do k = 1, nlev
      hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
   end do

   ! Initialize reference pressures
   call ref_pres_init(hypi, hypm, nprlev)

   if (iam < par%nprocs) then

      call prim_init1(elem, fvm, par, TimeLevel)
      if (fv_nphys > 0) then
         call dp_init(elem, fvm)
      end if

      if (fv_nphys > 0) then
         qsize_local = thermodynamic_active_species_num + 3
      else
         qsize_local = pcnst + 3
      end if

      call initEdgeBuffer(par, edgebuf, elem, qsize_local*nlev, nthreads=1)

   else  ! auxiliary processes

      globaluniquecols = 0
      nelem     = 0
      nelemd    = 0
      nelemdmax = 0
   endif

   ! nelemdmax is computed on the dycore comm, we need it globally.
   ngcols_d = nelemdmax
   call MPI_Allreduce(ngcols_d, nelemdmax, 1, MPI_INTEGER, MPI_MAX, mpicom, ierr)
   ! All pes might not have the correct global grid size
   call MPI_Allreduce(globaluniquecols, ngcols_d, 1, MPI_INTEGER, MPI_MAX, mpicom, ierr)
   ! All pes might not have the correct number of elements
   call MPI_Allreduce(nelem, nelem_d, 1, MPI_INTEGER, MPI_MAX, mpicom, ierr)

   ! nelemd (# of elements on this task) is set by prim_init1
   call init_loop_ranges(nelemd)

   ! Dynamics timestep
   !
   !  Note: dtime = timestep for physics/dynamics coupling
   !        tstep = the dynamics timestep:
   dtime = get_step_size()
   tstep = dtime / real(nsplit*qsplit*rsplit, r8)
   TimeLevel%nstep = get_nstep()*nsplit*qsplit*rsplit

   ! initial SE (subcycled) nstep
   TimeLevel%nstep0 = 0

   ! determine whether initial file uses 'ncol' or 'ncol_d'
   call get_hdim_name(fh_ini, ini_grid_hdim_name)

   ! Define the dynamics and physics grids on the dynamics decompostion.
   ! Physics grid on the physics decomposition is defined in phys_grid_init.
   call define_cam_grids()

   if (fv_nphys > 0) then

      ! ================================================
      ! finish fvm initialization
      ! ================================================

      if (iam < par%nprocs) then
         hybrid = config_thread_region(par,'serial')
         call get_loop_ranges(hybrid, ibeg=nets, iend=nete)

         ! initialize halo coordinate variables for cslam and physgrid
         call fvm_init2(elem, fvm, hybrid, nets, nete)
         call fvm_pg_init(elem, fvm, hybrid, nets, nete, irecons_tracer)
         call fvm_init3(elem, fvm, hybrid, nets, nete, irecons_tracer)
      end if

   end if

   ! write grid and mapping files
   if (se_write_gll_corners) then
      call write_grid_mapping(par, elem)
   end if

   if (trim(se_write_grid_file) /= "no") then
      if (fv_nphys > 0) then
         call dp_write(elem, fvm, trim(se_write_grid_file), trim(se_grid_filename))
      else
         call gll_grid_write(elem, trim(se_write_grid_file), trim(se_grid_filename))
      end if
   end if

   if (do_native_mapping) then

      allocate(areaA(ngcols_d))
      allocate(clat(ngcols_d),clon(ngcols_d))
      call get_horiz_grid_int(ngcols_d, clat_d_out=clat, clon_d_out=clon, area_d_out=areaA)

      ! Create mapping files using SE basis functions
      call create_native_mapping_files(par, elem, 'native', ngcols_d, clat, clon, areaa)
      call create_native_mapping_files(par, elem, 'bilin', ngcols_d, clat, clon, areaa)

      deallocate(areaa, clat, clon)
   end if

   call mpi_barrier(mpicom, ierr)

end subroutine dyn_grid_init

!==============================================================================

subroutine get_dyn_grid_info(hdim1_d, hdim2_d, num_lev,                       &
     index_model_top_layer, index_surface_layer, unstructured, dyn_columns)
   !------------------------------------------------------------
   !
   ! get_dyn_grid_info returns physics grid column information
   !
   !------------------------------------------------------------
   use shr_const_mod,          only: SHR_CONST_PI
   use cam_abortutils,         only: endrun
   use spmd_utils,             only: iam
   use coordinate_systems_mod, only: spherical_polar_t
   ! Dummy arguments
   integer,          intent(out)   :: hdim1_d ! # longitudes or grid size
   integer,          intent(out)   :: hdim2_d ! # latitudes or 1
   integer,          intent(out)   :: num_lev ! # levels
   integer,          intent(out)   :: index_model_top_layer
   integer,          intent(out)   :: index_surface_layer
   logical,          intent(out)   :: unstructured
   ! dyn_columns will contain a copy of the physics column info local to this
   ! dynamics task
   type(physics_column_t), allocatable, intent(out) :: dyn_columns(:)
   ! Local variables
   integer                         :: lindex
   integer                         :: gindex
   integer                         :: elem_ind, col_ind, ii, jj
   integer                         :: num_local_cols
   type(spherical_polar_t)         :: coord
   real(r8)                        :: dcoord
   real(r8),         parameter     :: radtodeg = 180.0_r8 / SHR_CONST_PI
   real(r8),         parameter     :: degtorad = SHR_CONST_PI / 180.0_r8
   character(len=*), parameter     :: subname = 'get_dyn_grid_info'

   unstructured = .true. ! SE is an unstructured dycore
   if (fv_nphys > 0) then ! physics uses an FVM grid
      num_local_cols = nelemd * fv_nphys * fv_nphys
   else
      num_local_cols = 0
      do elem_ind = 1, nelemd
         num_local_cols = num_local_cols + elem(elem_ind)%idxP%NumUniquePts
      end do
   end if
   if (allocated(local_dyn_columns)) then
      ! Check for correct number of columns
      if (size(local_dyn_columns) /= num_local_cols) then
         call endrun(subname//': called with inconsistent column numbers')
      end if
   else
      allocate(local_dyn_columns(num_local_cols))
      if (fv_nphys > 0) then ! physics uses an FVM grid
         hdim1_d = nelem * fv_nphys * fv_nphys
      else
         hdim1_d = ngcols_d
      end if
      hdim2_d = 1
      num_lev = nlev
      index_model_top_layer = 1
      index_surface_layer = nlev
      lindex = 0
      do elem_ind = 1, nelemd
         if (fv_nphys > 0) then ! physics uses an FVM grid
            do col_ind = 0, (fv_nphys * fv_nphys) - 1
               lindex = lindex + 1
               ii = MOD(col_ind, fv_nphys) + 1
               jj = (col_ind / fv_nphys) + 1
               coord = fvm(elem_ind)%center_cart_physgrid(ii, jj)
               local_dyn_columns(lindex)%lat_rad = coord%lat
               dcoord = local_dyn_columns(lindex)%lat_rad * radtodeg
               local_dyn_columns(lindex)%lat_deg = dcoord
               local_dyn_columns(lindex)%lon_rad = coord%lon
               dcoord = local_dyn_columns(lindex)%lon_rad * radtodeg
               local_dyn_columns(lindex)%lon_deg = dcoord
               local_dyn_columns(lindex)%area =                               &
                    fvm(elem_ind)%area_sphere_physgrid(ii,jj)
               local_dyn_columns(lindex)%weight =                             &
                    local_dyn_columns(lindex)%area
               ! File decomposition
               gindex = ((elem(elem_ind)%GlobalId-1) * fv_nphys * fv_nphys) + &
                    col_ind + 1
               local_dyn_columns(lindex)%global_col_num = gindex
               ! Note, coord_indices not used for unstructured dycores
               ! Dynamics decomposition
               local_dyn_columns(lindex)%dyn_task = iam
               local_dyn_columns(lindex)%local_dyn_block = elem_ind
               local_dyn_columns(lindex)%global_dyn_block =                   &
                    elem(elem_ind)%GlobalId
               allocate(local_dyn_columns(lindex)%dyn_block_index(1))
               local_dyn_columns(lindex)%dyn_block_index(1) = col_ind + 1
            end do
         else
            do col_ind = 1, elem(elem_ind)%idxP%NumUniquePts
               lindex = lindex + 1
               ii = elem(elem_ind)%idxP%ia(col_ind)
               jj = elem(elem_ind)%idxP%ja(col_ind)

               dcoord = elem(elem_ind)%spherep(ii,jj)%lat
               local_dyn_columns(lindex)%lat_rad = dcoord
               dcoord = local_dyn_columns(lindex)%lat_rad * radtodeg
               local_dyn_columns(lindex)%lat_deg = dcoord
               dcoord = elem(elem_ind)%spherep(ii,jj)%lon
               local_dyn_columns(lindex)%lon_rad = dcoord
               dcoord = local_dyn_columns(lindex)%lon_rad * radtodeg
               local_dyn_columns(lindex)%lon_deg = dcoord
               local_dyn_columns(lindex)%area =                               &
                    1.0_r8 / elem(elem_ind)%rspheremp(ii,jj)
               local_dyn_columns(lindex)%weight = local_dyn_columns(lindex)%area
               ! File decomposition
               gindex = elem(elem_ind)%idxP%UniquePtoffset + col_ind - 1
               local_dyn_columns(lindex)%global_col_num = gindex
               ! Note, coord_indices not used for unstructured dycores
               ! Dynamics decomposition
               local_dyn_columns(lindex)%dyn_task = iam
               local_dyn_columns(lindex)%local_dyn_block = elem_ind
               local_dyn_columns(lindex)%global_dyn_block =                   &
                    elem(elem_ind)%GlobalId
               allocate(local_dyn_columns(lindex)%dyn_block_index(1))
               local_dyn_columns(lindex)%dyn_block_index(1) = col_ind
            end do
         end if
      end do
   end if
   ! Copy the information to the output array
   if (allocated(dyn_columns)) then
      deallocate(dyn_columns)
   end if
   allocate(dyn_columns(lindex))
   do lindex = 1, num_local_cols
      dyn_columns(lindex) = local_dyn_columns(lindex)
   end do

   end subroutine get_dyn_grid_info

!==============================================================================

subroutine get_horiz_grid_dim_d(hdim1_d,hdim2_d)

   ! Returns declared horizontal dimensions of computational grid.
   ! For non-lon/lat grids, declare grid to be one-dimensional,
   ! i.e., (ngcols_d x 1)

   !------------------------------Arguments--------------------------------
   integer, intent(out)           :: hdim1_d ! first horizontal dimension
   integer, intent(out), optional :: hdim2_d ! second horizontal dimension
   !-----------------------------------------------------------------------

   if (fv_nphys > 0) then
      hdim1_d = fv_nphys*fv_nphys*nelem_d
   else
      hdim1_d = ngcols_d
   end if
   if (present(hdim2_d)) then
      hdim2_d = 1
   end if

end subroutine get_horiz_grid_dim_d

!=========================================================================================

subroutine get_horiz_grid_int(nxy, clat_d_out, clon_d_out, area_d_out, &
     wght_d_out, lat_d_out, lon_d_out)

   ! Return global arrays of latitude and longitude (in radians), column
   ! surface area (in radians squared) and surface integration weights for
   ! global column indices that will be passed to/from physics

   ! arguments
   integer, intent(in)   :: nxy                     ! array sizes

   real(r8), intent(out),         optional :: clat_d_out(:) ! column latitudes
   real(r8), intent(out),         optional :: clon_d_out(:) ! column longitudes
   real(r8), intent(out), target, optional :: area_d_out(:) ! column surface

   real(r8), intent(out), target, optional :: wght_d_out(:) ! column integration weight
   real(r8), intent(out),         optional :: lat_d_out(:)  ! column degree latitudes
   real(r8), intent(out),         optional :: lon_d_out(:)  ! column degree longitudes

   ! local variables
   real(r8),           pointer :: area_d(:)
   real(r8),           pointer :: temp(:)
   character(len=SHR_KIND_CL)  :: errormsg
   character(len=*), parameter :: sub = 'get_horiz_grid_int'
   !----------------------------------------------------------------------------

   ! check that nxy is set to correct size for global arrays
   if (fv_nphys > 0) then
      if (nxy < fv_nphys*fv_nphys*nelem_d) then
         write(errormsg, *) sub//': arrays too small; Passed',     &
            nxy, ', needs to be at least', fv_nphys*fv_nphys*nelem_d
         call endrun(errormsg)
      end if
   else
      if (nxy < ngcols_d) then
         write(errormsg,*) sub//': arrays not large enough; ',     &
            'Passed', nxy, ', needs to be at least', ngcols_d
         call endrun(errormsg)
      end if
   end if

   if ( present(area_d_out) ) then
      if (size(area_d_out) /= nxy) then
         call endrun(sub//': bad area_d_out array size')
      end if
      area_d => area_d_out
      call create_global_area(area_d)

   else if ( present(wght_d_out) ) then
      if (size(wght_d_out) /= nxy) then
         call endrun(sub//': bad wght_d_out array size')
      end if
      area_d => wght_d_out
      call create_global_area(area_d)

   end if

   ! If one of area_d_out  or wght_d_out was present, then it was computed
   ! above.  If they were *both* present, then do this:
   if ( present(area_d_out) .and. present(wght_d_out) ) then
      wght_d_out(:) = area_d_out(:)
   end if

   if (present(clon_d_out)) then
      if (size(clon_d_out) /= nxy) then
         call endrun(sub//': bad clon_d_out array size in dyn_grid')
      end if
   end if

   if (present(clat_d_out)) then

      if (size(clat_d_out) /= nxy) then
         call endrun('bad clat_d_out array size in dyn_grid')
      end if

      if (present(clon_d_out)) then
         call create_global_coords(clat_d_out, clon_d_out, lat_d_out, lon_d_out)
      else
         allocate(temp(nxy))
         call create_global_coords(clat_d_out, temp, lat_d_out, lon_d_out)
         deallocate(temp)
      end if

   else if (present(clon_d_out)) then

      allocate(temp(nxy))
      call create_global_coords(temp, clon_d_out, lat_d_out, lon_d_out)
      deallocate(temp)

   end if

end subroutine get_horiz_grid_int

!=========================================================================================

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)

   ! create list of attributes for the physics grid that should be copied
   ! from the corresponding grid object on the dynamics decomposition

   use cam_grid_support, only: max_hcoordname_len

   ! Dummy arguments
   character(len=max_hcoordname_len),          intent(out) :: gridname
   character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)

   if (fv_nphys > 0) then
      gridname = 'physgrid_d'
      allocate(grid_attribute_names(2))
      grid_attribute_names(1) = 'fv_nphys'
      grid_attribute_names(2) = 'ne'
   else
      gridname = 'GLL'
      allocate(grid_attribute_names(2))
      grid_attribute_names(1) = 'np'
      grid_attribute_names(2) = 'ne'
   end if

end subroutine physgrid_copy_attributes_d

!=========================================================================================

integer function get_dyn_grid_parm(name) result(ival)

   ! This function is in the process of being deprecated, but is still needed
   ! as a dummy interface to satisfy external references from some chemistry routines.

   use pmgrid,          only: plat, plev

   character(len=*), intent(in) :: name
   !----------------------------------------------------------------------------

   if (name.eq.'plat') then
      ival = plat
   else if(name.eq.'plon') then
      if (fv_nphys>0) then
         ival = fv_nphys*fv_nphys*nelem_d
      else
         ival = ngcols_d
      end if
   else if(name.eq.'plev') then
      ival = plev

   else
      ival = -1
   end if

end function get_dyn_grid_parm

!=========================================================================================

subroutine dyn_grid_get_colndx(igcol, ncols, owners, col, lbk)

   ! For each global column index return the owning task.  If the column is owned
   ! by this task, then also return the local block number and column index in that
   ! block.
   !
   ! NOTE: this routine needs to be updated for the physgrid

   integer, intent(in)  :: ncols
   integer, intent(in)  :: igcol(ncols)
   integer, intent(out) :: owners(ncols)
   integer, intent(out) :: col(ncols)
   integer, intent(out) :: lbk(ncols)

   !----------------------------------------------------------------------------

   owners = (igcol * 0) -1 ! Kill compiler warnings
   col    = -1             ! Kill compiler warnings
   lbk    = -1             ! Kill compiler warnings
   call endrun('dyn_grid_get_colndx: not implemented for unstructured grids')

end subroutine dyn_grid_get_colndx

!=========================================================================================

subroutine dyn_grid_get_elem_coords(ie, rlon, rlat, cdex)

   ! Returns coordinates of a specified block element of the dyn grid
   !
   ! NB: This routine only uses the GLL points (i.e, it ignores the physics
   !     grid). This is probably OK as current use is only for dyn_decomp
   !     variables in history.

   integer, intent(in) :: ie ! block element index

   real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the element
   real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the element
   integer, optional, intent(out) :: cdex(:) ! global column index

   integer :: sb,eb, ii, i,j, icol, igcol
   real(r8), allocatable :: clat(:), clon(:)
   !----------------------------------------------------------------------------

   if (fv_nphys > 0) then
      call endrun('dyn_grid_get_colndx: not implemented for the FVM physics grid')
   end if

   sb = elem(ie)%idxp%UniquePtOffset
   eb = sb + elem(ie)%idxp%NumUniquePts-1

   allocate( clat(sb:eb), clon(sb:eb) )
   call UniqueCoords( elem(ie)%idxP, elem(ie)%spherep, clat(sb:eb), clon(sb:eb) )

   if (present(cdex)) cdex(:) = -1
   if (present(rlat)) rlat(:) = -999._r8
   if (present(rlon)) rlon(:) = -999._r8

   do ii=1,elem(ie)%idxp%NumUniquePts
      i=elem(ie)%idxp%ia(ii)
      j=elem(ie)%idxp%ja(ii)
      icol = i+(j-1)*np
      igcol = elem(ie)%idxp%UniquePtoffset+ii-1
      if (present(cdex)) cdex(icol) = igcol
      if (present(rlat)) rlat(icol) = clat( igcol )
      if (present(rlon)) rlon(icol) = clon( igcol )
   end do

   deallocate( clat, clon )

end subroutine dyn_grid_get_elem_coords

!=========================================================================================
! Private routines.
!=========================================================================================

subroutine get_hdim_name(fh_ptr, grid_hdim_name)
   use pio, only: pio_inq_dimid, pio_seterrorhandling
   use pio, only: PIO_BCAST_ERROR, PIO_NOERR

   ! Determine whether the supplied file uses 'ncol' or 'ncol_d' horizontal
   ! dimension in the unstructured grid.  It is also possible when using
   ! analytic initial conditions that the file only contains
   ! vertical coordinates.
   ! Return 'none' if that is the case.

   ! Arguments
   type(file_desc_t),   pointer  :: fh_ptr
   character(len=6), intent(out) :: grid_hdim_name ! horizontal dimension name

   ! local variables
   integer  :: ierr, pio_errtype
   integer  :: ncol_did

   character(len=*), parameter :: sub = 'get_hdim_name'
   !----------------------------------------------------------------------------

   ! Set PIO to return error flags.
   call pio_seterrorhandling(fh_ptr, PIO_BCAST_ERROR, pio_errtype)

   ! Check for ncol_d first just in case the file also contains fields on
   ! the physics grid.
   ierr = pio_inq_dimid(fh_ptr, 'ncol_d', ncol_did)
   if (ierr == PIO_NOERR) then

      grid_hdim_name = 'ncol_d'

   else

      ! if 'ncol_d' not in file, check for 'ncol'
      ierr = pio_inq_dimid(fh_ptr, 'ncol', ncol_did)

      if (ierr == PIO_NOERR) then

         grid_hdim_name = 'ncol'

      else

         grid_hdim_name = 'none'

      end if
   end if

   ! Return PIO to previous error handling.
   call pio_seterrorhandling(fh_ptr, pio_errtype)

end subroutine get_hdim_name

subroutine define_cam_grids()

   ! Create grid objects on the dynamics decomposition for grids used by
   ! the dycore.  The decomposed grid object contains data for the elements
   ! in each task and information to map that data to the global grid.
   !
   ! Notes on dynamic memory management:
   !
   ! . Coordinate values and the map passed to the horiz_coord_create
   !   method are copied to the object.  The memory may be deallocated
   !   after the object is created.
   !
   ! . The area values passed to cam_grid_attribute_register are only pointed
   !   to by the attribute object, so that memory cannot be deallocated.  But the
   !   map is copied.
   !
   ! . The grid_map passed to cam_grid_register is just pointed to.
   !   Cannot be deallocated.

   use cam_grid_support, only: horiz_coord_t, horiz_coord_create
   use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
   use dimensions_mod,   only: nc

   ! Local variables
   integer                      :: i, ii, j, k, ie, mapind
   character(len=8)             :: latname, lonname, ncolname, areaname

   type(horiz_coord_t), pointer :: lat_coord
   type(horiz_coord_t), pointer :: lon_coord
   integer(iMap),       pointer :: grid_map(:,:)

   real(r8),        allocatable :: pelat_deg(:)  ! pe-local latitudes (degrees)
   real(r8),        allocatable :: pelon_deg(:)  ! pe-local longitudes (degrees)
   real(r8),        pointer     :: pearea(:) => null()  ! pe-local areas
   real(r8)                     :: areaw(np,np)
   integer(iMap)                :: fdofP_local(npsq,nelemd) ! pe-local map for dynamics decomp
   integer(iMap),   allocatable :: pemap(:)                 ! pe-local map for PIO decomp

   integer                      :: ncols_fvm, ngcols_fvm
   real(r8),        allocatable :: fvm_coord(:)
   real(r8),            pointer :: fvm_area(:)
   integer(iMap),       pointer :: fvm_map(:)

   integer                      :: ncols_physgrid, ngcols_physgrid
   real(r8),        allocatable :: physgrid_coord(:)
   real(r8),            pointer :: physgrid_area(:)
   integer(iMap),       pointer :: physgrid_map(:)
   !----------------------------------------------------------------------------

   !-----------------------
   ! Create GLL grid object
   !-----------------------

   ! Calculate the mapping between element GLL points and file order
   fdofp_local = 0_iMap
   do ie = 1, nelemd
      do ii = 1, elem(ie)%idxP%NumUniquePts
         i = elem(ie)%idxP%ia(ii)
         j = elem(ie)%idxP%ja(ii)
         fdofp_local((np*(j-1))+i,ie) = elem(ie)%idxP%UniquePtoffset + ii - 1
      end do
   end do

   allocate(pelat_deg(np*np*nelemd))
   allocate(pelon_deg(np*np*nelemd))
   allocate(pearea(np*np*nelemd))
   allocate(pemap(np*np*nelemd))

   pemap = 0_iMap
   ii = 1
   do ie = 1, nelemd
      areaw = 1.0_r8 / elem(ie)%rspheremp(:,:)
      pearea(ii:ii+npsq-1) = reshape(areaw, (/ np*np /))
      pemap(ii:ii+npsq-1) = fdofp_local(:,ie)
      do j = 1, np
         do i = 1, np
            pelat_deg(ii) = elem(ie)%spherep(i,j)%lat * rad2deg
            pelon_deg(ii) = elem(ie)%spherep(i,j)%lon * rad2deg
            ii = ii + 1
         end do
      end do
   end do

   ! If using the physics grid then the GLL grid will use the names with
   ! '_d' suffixes and the physics grid will use the unadorned names.
   ! This allows fields on both the GLL and physics grids to be written to history
   ! output files.
   if (trim(ini_grid_hdim_name) == 'ncol_d') then
      latname  = 'lat_d'
      lonname  = 'lon_d'
      ncolname = 'ncol_d'
      areaname = 'area_d'
   else
      latname  = 'lat'
      lonname  = 'lon'
      ncolname = 'ncol'
      areaname = 'area'
   end if
   lat_coord => horiz_coord_create('lat_d', 'ncol_d', ngcols_d,  &
         'latitude', 'degrees_north', 1, size(pelat_deg), pelat_deg, map=pemap)
   lon_coord => horiz_coord_create('lon_d', 'ncol_d', ngcols_d,  &
         'longitude', 'degrees_east', 1, size(pelon_deg), pelon_deg, map=pemap)

   ! Map for GLL grid
   allocate(grid_map(3,npsq*nelemd))
   grid_map = 0_iMap
   mapind = 1
   do j = 1, nelemd
      do i = 1, npsq
         grid_map(1, mapind) = i
         grid_map(2, mapind) = j
         grid_map(3, mapind) = pemap(mapind)
         mapind = mapind + 1
      end do
   end do

   ! The native SE GLL grid
   call cam_grid_register('GLL', dyn_decomp, lat_coord, lon_coord,           &
         grid_map, block_indexed=.false., unstruct=.true.)
   call cam_grid_attribute_register('GLL', 'area_d', 'gll grid areas', &
         'ncol_d', pearea, map=pemap)
   call cam_grid_attribute_register('GLL', 'np', '', np)
   call cam_grid_attribute_register('GLL', 'ne', '', ne)

   ! If dim name is 'ncol', create INI grid
   ! We will read from INI grid, but use GLL grid for all output
   if (trim(ini_grid_hdim_name) == 'ncol') then

      lat_coord => horiz_coord_create('lat', 'ncol', ngcols_d,  &
         'latitude', 'degrees_north', 1, size(pelat_deg), pelat_deg, map=pemap)
      lon_coord => horiz_coord_create('lon', 'ncol', ngcols_d,  &
         'longitude', 'degrees_east', 1, size(pelon_deg), pelon_deg, map=pemap)

      call cam_grid_register('INI', ini_decomp, lat_coord, lon_coord,         &
         grid_map, block_indexed=.false., unstruct=.true.)
      call cam_grid_attribute_register('INI', 'area', 'ini grid areas', &
               'ncol', pearea, map=pemap)

      ini_grid_name = 'INI'
   else
      ! The dyn_decomp grid can be used to read the initial file.
      ini_grid_name = 'GLL'
   end if

   ! Coordinate values and maps are copied into the coordinate and attribute objects.
   ! Locally allocated storage is no longer needed.
   deallocate(pelat_deg)
   deallocate(pelon_deg)
   deallocate(pemap)

   ! pearea cannot be deallocated as the attribute object is just pointing
   ! to that memory.  It can be nullified since the attribute object has
   ! the reference.
   nullify(pearea)

   ! grid_map cannot be deallocated as the cam_filemap_t object just points
   ! to it.  It can be nullified.
   nullify(grid_map)

   !---------------------------------
   ! Create FVM grid object for CSLAM
   !---------------------------------

   if (ntrac > 0) then

      ncols_fvm = nc * nc * nelemd
      ngcols_fvm = nc * nc * nelem_d
      allocate(fvm_coord(ncols_fvm))
      allocate(fvm_map(ncols_fvm))
      allocate(fvm_area(ncols_fvm))

      do ie = 1, nelemd
         k = 1
         do j = 1, nc
            do i = 1, nc
               mapind = k + ((ie - 1) * nc * nc)
               fvm_coord(mapind) = fvm(ie)%center_cart(i,j)%lon*rad2deg
               fvm_map(mapind) = k + ((elem(ie)%GlobalId-1) * nc * nc)
               fvm_area(mapind) = fvm(ie)%area_sphere(i,j)
               k = k + 1
            end do
         end do
      end do
      lon_coord => horiz_coord_create('lon_fvm', 'ncol_fvm', ngcols_fvm,      &
           'longitude', 'degrees_east', 1, size(fvm_coord), fvm_coord,        &
           map=fvm_map)

      do ie = 1, nelemd
         k = 1
         do j = 1, nc
            do i = 1, nc
               mapind = k + ((ie - 1) * nc * nc)
               fvm_coord(mapind) = fvm(ie)%center_cart(i,j)%lat*rad2deg
               k = k + 1
            end do
         end do
      end do
      lat_coord => horiz_coord_create('lat_fvm', 'ncol_fvm', ngcols_fvm,      &
           'latitude', 'degrees_north', 1, size(fvm_coord), fvm_coord,        &
           map=fvm_map)

      ! Map for FVM grid
      allocate(grid_map(3, ncols_fvm))
      grid_map = 0_iMap
      mapind = 1
      do j = 1, nelemd
         do i = 1, nc*nc
            grid_map(1, mapind) = i
            grid_map(2, mapind) = j
            grid_map(3, mapind) = fvm_map(mapind)
            mapind = mapind + 1
         end do
      end do

      ! create FVM (CSLAM) grid object
      call cam_grid_register('FVM', fvm_decomp, lat_coord, lon_coord,         &
           grid_map, block_indexed=.false., unstruct=.true.)
      call cam_grid_attribute_register('FVM', 'area_fvm', 'fvm grid areas',   &
           'ncol_fvm', fvm_area, map=fvm_map)
      call cam_grid_attribute_register('FVM', 'nc', '', nc)
      call cam_grid_attribute_register('FVM', 'ne', '', ne)

      deallocate(fvm_coord)
      deallocate(fvm_map)
      nullify(fvm_area)
      nullify(grid_map)

   end if

   !------------------------------------------------------------------
   ! Create grid object for physics grid on the dynamics decomposition
   !------------------------------------------------------------------

   if (fv_nphys > 0) then

      ncols_physgrid = fv_nphys * fv_nphys * nelemd
      ngcols_physgrid = fv_nphys * fv_nphys * nelem_d
      allocate(physgrid_coord(ncols_physgrid))
      allocate(physgrid_map(ncols_physgrid))
      allocate(physgrid_area(ncols_physgrid))

      do ie = 1, nelemd
         k = 1
         do j = 1, fv_nphys
            do i = 1, fv_nphys
               mapind = k + ((ie - 1) * fv_nphys * fv_nphys)
               physgrid_coord(mapind) = fvm(ie)%center_cart_physgrid(i,j)%lon*rad2deg
               physgrid_map(mapind) = k + ((elem(ie)%GlobalId-1) * fv_nphys * fv_nphys)
               physgrid_area(mapind) = fvm(ie)%area_sphere_physgrid(i,j)
               k = k + 1
            end do
         end do
      end do
      lon_coord => horiz_coord_create('lon', 'ncol', ngcols_physgrid,      &
           'longitude', 'degrees_east', 1, size(physgrid_coord), physgrid_coord,   &
           map=physgrid_map)

      do ie = 1, nelemd
         k = 1
         do j = 1, fv_nphys
            do i = 1, fv_nphys
               mapind = k + ((ie - 1) * fv_nphys * fv_nphys)
               physgrid_coord(mapind) = fvm(ie)%center_cart_physgrid(i,j)%lat*rad2deg
               k = k + 1
            end do
         end do
      end do
      lat_coord => horiz_coord_create('lat', 'ncol', ngcols_physgrid,      &
           'latitude', 'degrees_north', 1, size(physgrid_coord), physgrid_coord,   &
           map=physgrid_map)

      ! Map for physics grid
      allocate(grid_map(3, ncols_physgrid))
      grid_map = 0_iMap
      mapind = 1
      do j = 1, nelemd
         do i = 1, fv_nphys*fv_nphys
            grid_map(1, mapind) = i
            grid_map(2, mapind) = j
            grid_map(3, mapind) = physgrid_map(mapind)
            mapind = mapind + 1
         end do
      end do

      ! create physics grid object
      call cam_grid_register('physgrid_d', physgrid_d, lat_coord, lon_coord, &
           grid_map, block_indexed=.false., unstruct=.true.)
      call cam_grid_attribute_register('physgrid_d', 'area_physgrid', 'physics grid areas',   &
           'ncol', physgrid_area, map=physgrid_map)
      call cam_grid_attribute_register('physgrid_d', 'fv_nphys', '', fv_nphys)
      call cam_grid_attribute_register('physgrid_d', 'ne',       '', ne)

      deallocate(physgrid_coord)
      deallocate(physgrid_map)
      nullify(physgrid_area)
      nullify(grid_map)

   end if

   nullify(lat_coord)         ! Belongs to grid
   nullify(lon_coord)         ! Belongs to grid

end subroutine define_cam_grids

!========================================================================================

subroutine write_grid_mapping(par, elem)

   use parallel_mod,  only: parallel_t
   use cam_pio_utils, only: cam_pio_createfile, pio_subsystem
   use pio,           only: pio_def_dim, var_desc_t, pio_int, pio_def_var, &
                            pio_enddef, pio_closefile, pio_initdecomp, io_desc_t, &
                            pio_write_darray, pio_freedecomp
   use dof_mod,       only: createmetadata

   ! arguments
   type(parallel_t), intent(in) :: par
   type(element_t),  intent(in) :: elem(:)

   ! local variables
   integer, parameter :: npm12 = (np-1)*(np-1)

   type(file_desc_t) :: nc
   type(var_desc_t)  :: vid
   type(io_desc_t)   :: iodesc
   integer :: dim1, dim2, ierr, i, j, ie, cc, base, ii, jj
   integer :: subelement_corners(npm12*nelemd,4)
   integer :: dof(npm12*nelemd*4)
   !----------------------------------------------------------------------------

   ! Create a CS grid mapping file for postprocessing tools

   ! write meta data for physics on GLL nodes
   call cam_pio_createfile(nc, 'SEMapping.nc', 0)

   ierr = pio_def_dim(nc, 'ncenters', npm12*nelem_d, dim1)
   ierr = pio_def_dim(nc, 'ncorners', 4, dim2)
   ierr = pio_def_var(nc, 'element_corners', PIO_INT, (/dim1,dim2/), vid)

   ierr = pio_enddef(nc)
   call createmetadata(par, elem, subelement_corners)

   jj=0
   do cc = 0, 3
      do ie = 1, nelemd
         base = ((elem(ie)%globalid-1)+cc*nelem_d)*npm12
         ii=0
         do j = 1, np-1
            do i = 1, np-1
               ii=ii+1
               jj=jj+1
               dof(jj) = base+ii
            end do
         end do
      end do
   end do

   call pio_initdecomp(pio_subsystem, pio_int, (/nelem_d*npm12,4/), dof, iodesc)

   call pio_write_darray(nc, vid, iodesc, &
                         reshape(subelement_corners, (/nelemd*npm12*4/)), ierr)

   call pio_freedecomp(nc, iodesc)

   call pio_closefile(nc)

end subroutine write_grid_mapping

!=========================================================================================

subroutine create_global_area(area_d)
   use dp_mapping, only: dp_reorder, dp_allocate, dp_deallocate

   ! Gather global array of column areas for the physics grid,
   ! reorder to global column order, then broadcast it to all tasks.

   ! Input variables
   real(r8), pointer           :: area_d(:)

   ! Local variables
   real(r8)                    :: areaw(np,np)
   real(r8), allocatable       :: rbuf(:), dp_area(:,:)
   integer                     :: rdispls(npes), recvcounts(npes)
   integer                     :: ncol
   integer                     :: ie, sb, eb, i, j, k
   integer                     :: ierr
   integer                     :: ibuf
   character(len=*), parameter :: sub = 'create_global_area'
   !----------------------------------------------------------------------------

   if (masterproc) then
      write(iulog, *) sub//': INFO: Non-scalable action: gathering global area in SE dycore.'
   end if

   if (fv_nphys > 0) then ! physics uses an FVM grid

      ! first gather all data onto masterproc, in mpi task order (via
      ! mpi_gatherv) then redorder into globalID order (via dp_reorder)
      ncol = fv_nphys*fv_nphys*nelem_d
      allocate(rbuf(ncol))
      allocate(dp_area(fv_nphys*fv_nphys,nelem_d))

      do ie = 1, nelemd
         k = 1
         do j = 1, fv_nphys
            do i = 1, fv_nphys
               dp_area(k,ie) = fvm(ie)%area_sphere_physgrid(i,j)
               k = k + 1
            end do
         end do
      end do

      call mpi_gather(nelemd*fv_nphys*fv_nphys, 1, mpi_integer, recvcounts, 1, &
                      mpi_integer, mstrid, mpicom, ierr)
      ! Figure global displacements
      if (masterproc) then
         rdispls(1) = 0
         do ie = 2, npes
            rdispls(ie) = rdispls(ie-1) + recvcounts(ie-1)
         end do
         ! Check to make sure we counted correctly
         if (rdispls(npes) + recvcounts(npes) /= ncol) then
            call endrun(sub//': bad rdispls array size')
         end if
      end if

      ! Gather up the areas onto the masterproc
      call mpi_gatherv(dp_area, fv_nphys*fv_nphys*nelemd, mpi_real8, rbuf, &
                       recvcounts, rdispls, mpi_real8, mstrid, mpicom, ierr)

      ! Reorder to global order
      call dp_allocate(elem)
      if (masterproc) call dp_reorder(rbuf, area_d)
      call dp_deallocate()

      ! Send everyone else the data
      call mpi_bcast(area_d, ncol, mpi_real8, mstrid, mpicom, ierr)

      deallocate(dp_area)

   else ! physics is on the GLL grid

      allocate(rbuf(ngcols_d))
      do ie = 1, nelemdmax
         if (ie <= nelemd) then
            rdispls(iam+1)    = elem(ie)%idxp%UniquePtOffset - 1
            eb                = rdispls(iam+1) + elem(ie)%idxp%NumUniquePts
            recvcounts(iam+1) = elem(ie)%idxP%NumUniquePts
            areaw = 1.0_r8 / elem(ie)%rspheremp(:,:)
            call UniquePoints(elem(ie)%idxP, areaw, area_d(rdispls(iam+1)+1:eb))
         else
            rdispls(iam+1) = 0
            recvcounts(iam+1) = 0
         end if

         ibuf = rdispls(iam+1)
         call mpi_allgather(ibuf, 1, mpi_integer, rdispls, &
            1, mpi_integer, mpicom, ierr)

         ibuf = recvcounts(iam+1)
         call mpi_allgather(ibuf, 1, mpi_integer, recvcounts, &
            1, mpi_integer, mpicom, ierr)

         sb = rdispls(iam+1) + 1
         eb = rdispls(iam+1) + recvcounts(iam+1)

         rbuf(1:recvcounts(iam+1)) = area_d(sb:eb)
         call mpi_allgatherv(rbuf, recvcounts(iam+1), mpi_real8, area_d,       &
            recvcounts(:), rdispls(:), mpi_real8, mpicom, ierr)
      end do

   end if

   deallocate(rbuf)

end subroutine create_global_area

!=========================================================================================

subroutine create_global_coords(clat, clon, lat_out, lon_out)
   use dp_mapping, only: dp_reorder, dp_allocate, dp_deallocate

   ! Gather global arrays of column coordinates for the physics grid,
   ! reorder to global column order, then broadcast to all tasks.

   ! arguments
   real(r8),           intent(out) :: clat(:)
   real(r8),           intent(out) :: clon(:)
   real(r8), optional, intent(out) :: lat_out(:)
   real(r8), optional, intent(out) :: lon_out(:)

   ! Local variables
   real(r8), allocatable           :: rbuf(:), dp_lon(:,:), dp_lat(:,:)
   integer                         :: rdispls(npes), recvcounts(npes)
   integer                         :: ie, sb, eb, i, j, k
   integer                         :: ierr
   integer                         :: ibuf
   integer                         :: ncol
   character(len=*), parameter :: sub='create_global_coords'
   !----------------------------------------------------------------------------

   if (masterproc) then
      write(iulog, *) sub//': INFO: Non-scalable action: Creating global coords in SE dycore.'
   end if

   clat(:) = -iam
   clon(:) = -iam
   if (present(lon_out)) then
      lon_out(:) = -iam
   end if
   if (present(lat_out)) then
      lat_out(:) = -iam
   end if

   if (fv_nphys > 0) then  ! physics uses an FVM grid

      ! first gather all data onto masterproc, in mpi task order (via
      ! mpi_gatherv) then redorder into globalID order (via dp_reorder)

      ncol = fv_nphys*fv_nphys*nelem_d
      allocate(rbuf(ncol))
      allocate(dp_lon(fv_nphys*fv_nphys,nelem_d))
      allocate(dp_lat(fv_nphys*fv_nphys,nelem_d))

      do ie = 1, nelemd
         k = 1
         do j = 1, fv_nphys
            do i = 1, fv_nphys
               dp_lon(k,ie) = fvm(ie)%center_cart_physgrid(i,j)%lon   ! radians
               dp_lat(k,ie) = fvm(ie)%center_cart_physgrid(i,j)%lat
               k = k + 1
            end do
         end do
      end do

      call mpi_gather(nelemd*fv_nphys*fv_nphys, 1, mpi_integer, recvcounts, &
                      1, mpi_integer, mstrid, mpicom, ierr)

      ! Figure global displacements
      if (masterproc) then
         rdispls(1) = 0
         do ie = 2, npes
            rdispls(ie) = rdispls(ie-1) + recvcounts(ie-1)
         end do
         ! Check to make sure we counted correctly
         if (rdispls(npes) + recvcounts(npes) /= ncol) then
            call endrun(sub//': bad rdispls array size')
         end if
      end if

      ! Gather up global latitudes
      call mpi_gatherv(dp_lat, fv_nphys*fv_nphys*nelemd, mpi_real8, rbuf,     &
                       recvcounts, rdispls, mpi_real8, mstrid, mpicom, ierr)

      ! Reorder to global order
      call dp_allocate(elem)
      if (masterproc) call dp_reorder(rbuf, clat)

      ! Send everyone else the data
      call mpi_bcast(clat, ncol, mpi_real8, mstrid, mpicom, ierr)

      ! Gather up global longitudes
      call mpi_gatherv(dp_lon, fv_nphys*fv_nphys*nelemd, mpi_real8, rbuf,     &
                       recvcounts, rdispls, mpi_real8, mstrid, mpicom, ierr)

      ! Reorder to global order
      if (masterproc) call dp_reorder(rbuf, clon)
      call dp_deallocate()

      ! Send everyone else the data
      call mpi_bcast(clon, ncol, mpi_real8, mstrid, mpicom, ierr)

      ! Create degree versions if requested
      if (present(lat_out)) then
         lat_out(:) = clat(:) * rad2deg
      end if
      if (present(lon_out)) then
         lon_out(:) = clon(:) * rad2deg
      end if

      deallocate(dp_lon)
      deallocate(dp_lat)

   else ! physics uses the GLL grid

      allocate(rbuf(ngcols_d))

      do ie = 1, nelemdmax

         if(ie <= nelemd) then
            rdispls(iam+1)    = elem(ie)%idxp%UniquePtOffset - 1
            eb                = rdispls(iam+1) + elem(ie)%idxp%NumUniquePts
            recvcounts(iam+1) = elem(ie)%idxP%NumUniquePts

            call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep, &
                              clat(rdispls(iam+1)+1:eb), clon(rdispls(iam+1)+1:eb))

            if (present(lat_out)) then
               lat_out(rdispls(iam+1)+1:eb) = clat(rdispls(iam+1)+1:eb) * rad2deg
            end if

            if (present(lon_out)) then
               lon_out(rdispls(iam+1)+1:eb) = clon(rdispls(iam+1)+1:eb) * rad2deg
            end if

         else
            rdispls(iam+1) = 0
            recvcounts(iam+1) = 0
         end if

         ibuf = rdispls(iam+1)
         call mpi_allgather(ibuf, 1, mpi_integer, rdispls, &
                            1, mpi_integer, mpicom, ierr)

         ibuf = recvcounts(iam+1)
         call mpi_allgather(ibuf, 1, mpi_integer, recvcounts, &
                            1, mpi_integer, mpicom, ierr)

         sb = rdispls(iam+1) + 1
         eb = rdispls(iam+1) + recvcounts(iam+1)

         rbuf(1:recvcounts(iam+1)) = clat(sb:eb)  ! whats going to happen if end=0?
         call mpi_allgatherv(rbuf, recvcounts(iam+1), mpi_real8, clat,         &
                             recvcounts(:), rdispls(:), mpi_real8, mpicom, ierr)

         if (present(lat_out)) then
            rbuf(1:recvcounts(iam+1)) = lat_out(sb:eb)
            call mpi_allgatherv(rbuf, recvcounts(iam+1), mpi_real8, lat_out,    &
                                recvcounts(:), rdispls(:), mpi_real8, mpicom, ierr)
         end if

         rbuf(1:recvcounts(iam+1)) = clon(sb:eb)
         call mpi_allgatherv(rbuf, recvcounts(iam+1), mpi_real8, clon,         &
                             recvcounts(:), rdispls(:), mpi_real8, mpicom, ierr)

         if (present(lon_out)) then
            rbuf(1:recvcounts(iam+1)) = lon_out(sb:eb)
            call mpi_allgatherv(rbuf, recvcounts(iam+1), mpi_real8, lon_out,    &
                                recvcounts(:), rdispls(:), mpi_real8, mpicom, ierr)
         end if

      end do  ! ie = 1, nelemdmax

   end if  ! (fv_nphys > 0)

end subroutine create_global_coords

!=============================================================================
!==
!!!!!! DUMMY INTERFACE TO TEST WEAK SCALING FIX, THIS SHOULD GO AWAY
!==
!=============================================================================

subroutine get_horiz_grid_d(nxy, clat_d_out, clon_d_out, area_d_out, &
                            wght_d_out, lat_d_out, lon_d_out)

   ! Return global arrays of latitude and longitude (in radians), column
   ! surface area (in radians squared) and surface integration weights for
   ! global column indices that will be passed to/from physics

   ! arguments
   integer, intent(in)   :: nxy                     ! array sizes

   real(r8), intent(out),         optional :: clat_d_out(:) ! column latitudes
   real(r8), intent(out),         optional :: clon_d_out(:) ! column longitudes
   real(r8), intent(out), target, optional :: area_d_out(:) ! column surface

   real(r8), intent(out), target, optional :: wght_d_out(:) ! column integration weight
   real(r8), intent(out),         optional :: lat_d_out(:)  ! column degree latitudes
   real(r8), intent(out),         optional :: lon_d_out(:)  ! column degree longitudes
   character(len=*), parameter :: subname = 'get_horiz_grid_d'

   call endrun(subname//': NOT SUPPORTED WITH WEAK SCALING FIX')
end subroutine get_horiz_grid_d

!==============================================================================

function get_dyn_grid_parm_real1d(name) result(rval)

   ! This routine is not used for SE, but still needed as a dummy interface to satisfy
   ! references from mo_synoz.F90

   character(len=*), intent(in) :: name
   real(r8), pointer :: rval(:)

   if(name.eq.'w') then
      call endrun('get_dyn_grid_parm_real1d: w not defined')
   else if(name.eq.'clat') then
      call endrun('get_dyn_grid_parm_real1d: clat not supported')
   else if(name.eq.'latdeg') then
      call endrun('get_dyn_grid_parm_real1d: latdeg not defined')
   else
      nullify(rval)
   end if
end function get_dyn_grid_parm_real1d

!==============================================================================

end module dyn_grid
