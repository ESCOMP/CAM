!===============================================================================
!===============================================================================
module soil_erod_mod
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,      only: iulog
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun

  implicit none
  private

  public :: soil_erod_init
  public :: soil_erodibility
  public :: soil_erod_fact

  real(r8), allocatable ::  soil_erodibility(:,:)  ! soil erodibility factor
  real(r8) :: soil_erod_fact                       ! tuning parameter for dust emissions

contains

  !=============================================================================
  !=============================================================================
  subroutine soil_erod_init( dust_emis_fact, soil_erod_file )
    use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
    use ppgrid,           only: begchunk, endchunk, pcols
    use mo_constants,     only: pi, d2r
    use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use ioFileMod,        only: getfil
    use netcdf,           only: nf90_open, nf90_nowrite, nf90_close, nf90_noerr,  &
                                nf90_inq_dimid, nf90_inquire_dimension,           &
                                nf90_inq_varid, nf90_get_var
    use cam_shmem_mod,    only: cam_shmem_alloc_r8_2d, cam_shmem_fence,           &
                                cam_shmem_free, cam_shmem_is_leader,              &
                                cam_shmem_leader_comm, cam_shmem_npes_per_node
#ifdef SPMD
    use mpishorthand,     only: mpicom, mpiint, mpir8
#endif

    real(r8),         intent(in) :: dust_emis_fact
    character(len=*), intent(in) :: soil_erod_file

    real(r8), pointer     ::  soil_erodibility_in(:,:) => null()  ! node-shared input field
    real(r8), allocatable :: dst_lons(:)
    real(r8), allocatable :: dst_lats(:)
    character(len=cl)     :: infile
    integer :: did, vid, nlat, nlon
    integer :: ncid, iret
    integer :: win_serod = -1   ! per-node shared-memory window for soil_erodibility_in

    type(interp_type) :: lon_wgts, lat_wgts
    real(r8) :: to_lats(pcols), to_lons(pcols)
    integer :: c, ncols, ierr
    real(r8), parameter :: zero=0._r8, twopi=2._r8*pi

    soil_erod_fact = dust_emis_fact

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
       write(iulog,*) 'soil_erod_mod: soil_erod_fact = ', soil_erod_fact
    end if

    ! for soil erodibility in mobilization, apply inside CAM instead of lsm.
    ! read in soil erodibility factors, similar to Zender's boundary conditions

    ! Get file name.  
    call getfil(soil_erod_file, infile, 0)

    ! Read the source-grid soil erodibility on masterproc only (serial NetCDF) and
    ! hold it in per-node MPI-3 shared memory: one physical copy per node instead
    ! of one per MPI rank.  At high resolution (e.g. ne1024) this removes
    ! (ranks_per_node - 1) redundant nlon*nlat copies of the input field.
    if (masterproc) then
       iret = nf90_open(trim(infile), NF90_NOWRITE, ncid)
       if (iret /= nf90_noerr) call endrun('soil_erod_init: failed to open '//trim(infile))
       iret = nf90_inq_dimid( ncid, 'lon', did )
       iret = nf90_inquire_dimension( ncid, did, len=nlon )
       iret = nf90_inq_dimid( ncid, 'lat', did )
       iret = nf90_inquire_dimension( ncid, did, len=nlat )
    end if

#ifdef SPMD
    call mpibcast( nlon, 1, mpiint, 0, mpicom )
    call mpibcast( nlat, 1, mpiint, 0, mpicom )
#endif

    allocate(dst_lons(nlon))
    allocate(dst_lats(nlat))
    call cam_shmem_alloc_r8_2d( soil_erodibility_in, win_serod, nlon, nlat )
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: soil_erodibility_in held in per-node shared memory; ', &
            cam_shmem_npes_per_node()-1, ' redundant copies/node avoided'
    end if

    ! Open the window epoch, fill on the node leader (masterproc), then publish.
    call cam_shmem_fence( win_serod )

    if (masterproc) then
       iret = nf90_inq_varid( ncid, 'lon', vid )
       iret = nf90_get_var( ncid, vid, dst_lons )
       iret = nf90_inq_varid( ncid, 'lat', vid )
       iret = nf90_get_var( ncid, vid, dst_lats )
       iret = nf90_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
       iret = nf90_get_var( ncid, vid, soil_erodibility_in )
       iret = nf90_close( ncid )
    end if

#ifdef SPMD
    ! Source coordinates to every rank; the large field to every node leader.  The
    ! closing fence then publishes it to every rank on each node.
    call mpibcast( dst_lons, nlon, mpir8, 0, mpicom )
    call mpibcast( dst_lats, nlat, mpir8, 0, mpicom )
    if (cam_shmem_is_leader()) then
       call mpibcast( soil_erodibility_in, nlon*nlat, mpir8, 0, cam_shmem_leader_comm() )
    end if
#endif
    call cam_shmem_fence( win_serod )

    !-----------------------------------------------------------------------
    !     	... convert to radians and setup regridding
    !-----------------------------------------------------------------------
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

    allocate( soil_erodibility(pcols,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to allocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to allocate soil_erodibility_in')
    end if

    !-----------------------------------------------------------------------
    !     	... regrid ..
    !-----------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(soil_erodibility_in(:,:), nlon,nlat , soil_erodibility(:,c), ncols, lon_wgts,lat_wgts)

       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    ! Release the node-shared input field (collective over the node communicator).
    call cam_shmem_free( soil_erodibility_in, win_serod )

    deallocate( dst_lats )
    deallocate( dst_lons )

  end  subroutine soil_erod_init

end module soil_erod_mod
