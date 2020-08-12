module restart_dynamics

! Write and read dynamics fields from the restart file.  For exact restart
! it is necessary to write all element data, including duplicate columns,
! to the file.

    use cam_abortutils,   only: endrun
    use cam_grid_support, only: cam_grid_header_info_t, cam_grid_id, cam_grid_write_attr, &
                                cam_grid_write_var, cam_grid_get_decomp, cam_grid_dimensions, max_hcoordname_len
    use cam_logfile,      only: iulog
    use cam_pio_utils,    only: cam_pio_handle_error
    use dyn_comp,         only: dyn_import_t, dyn_export_t
    use dyn_grid,         only: mytile
    use fv_arrays_mod,    only: fv_atmos_type
    use pio,              only: file_desc_t, var_desc_t
    use shr_kind_mod,     only: r8 => shr_kind_r8, i8 => shr_kind_i8
    use spmd_utils,       only: masterproc

    implicit none
    private

    public :: init_restart_dynamics, write_restart_dynamics, read_restart_dynamics

    type(var_desc_t) :: udesc, vdesc, tdesc, psdesc, phisdesc, usdesc,vsdesc,delpdesc,omegadesc

    integer :: ncol_d_dimid, ncol_d_ew_dimid, ncol_d_ns_dimid,  nlev_dimid, nlevp_dimid
    type(var_desc_t), allocatable :: qdesc(:)
    integer :: is,ie,js,je


!=======================================================================
contains
!=======================================================================

subroutine init_restart_dynamics(File, dyn_out)

    use constituents,  only: cnst_name, pcnst
    use hycoef,        only: init_restart_hycoef
    use pio,           only: pio_unlimited, pio_double, pio_def_dim, &
                             pio_seterrorhandling, pio_bcast_error, &
                             pio_def_var, &
                             pio_inq_dimid

    ! arguments
    type(file_desc_t),  intent(inout) :: file
    type(dyn_export_t), intent(in)    :: dyn_out

    ! local variables
    integer :: vdimids(2)
    integer :: ierr, i, err_handling
    integer :: time_dimid
    integer :: is,ie,js,je
    type (fv_atmos_type),  pointer :: Atm(:)
    
    integer :: grid_id,grid_id_ns,grid_id_ew
    type(cam_grid_header_info_t) :: info,info_ew,info_ns

    !---------------------------------------------------------------------------

    Atm=>dyn_out%atm

    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    call init_restart_hycoef(File, vdimids)

    call pio_seterrorhandling(File, pio_bcast_error, err_handling)

    ierr = PIO_Def_Dim(File, 'time', PIO_UNLIMITED, time_dimid)

    grid_id = cam_grid_id('FFSL')
    call cam_grid_write_attr(File, grid_id, info)
    ncol_d_dimid = info%get_hdimid(1)

    grid_id_ew = cam_grid_id('FFSL_EW')
    call cam_grid_write_attr(File, grid_id_ew, info_ew)
    ncol_d_ew_dimid = info_ew%get_hdimid(1)

    grid_id_ns = cam_grid_id('FFSL_NS')
    call cam_grid_write_attr(File, grid_id_ns, info_ns)
    ncol_d_ns_dimid = info_ns%get_hdimid(1)

    nlev_dimid  = vdimids(1)

    ierr = PIO_Def_Var(File, 'U', pio_double, (/ncol_d_dimid, nlev_dimid/), Udesc)
    ierr = PIO_Def_Var(File, 'V', pio_double, (/ncol_d_dimid, nlev_dimid/), Vdesc)
    ierr = PIO_Def_Var(File, 'US', pio_double, (/ncol_d_ns_dimid, nlev_dimid/), USdesc)
    ierr = PIO_Def_Var(File, 'VS', pio_double, (/ncol_d_ew_dimid, nlev_dimid/), VSdesc)
    ierr = PIO_Def_Var(File, 'T', pio_double, (/ncol_d_dimid, nlev_dimid/), Tdesc)
    ierr = PIO_Def_Var(File, 'OMEGA', pio_double, (/ncol_d_dimid, nlev_dimid/), omegadesc)
    ierr = PIO_Def_Var(File, 'DELP', pio_double, (/ncol_d_dimid, nlev_dimid/), delpdesc)
    ierr = PIO_Def_Var(File, 'PS', pio_double, (/ncol_d_dimid/), PSdesc)
    ierr = PIO_Def_Var(File, 'PHIS', pio_double, (/ncol_d_dimid/), phisdesc)

    allocate(Qdesc(pcnst))

    do i = 1, pcnst
        ierr = PIO_Def_Var(File, cnst_name(i), pio_double, (/ncol_d_dimid, nlev_dimid/), Qdesc(i))
    end do

    call pio_seterrorhandling(File, err_handling)

end subroutine init_restart_dynamics

!=======================================================================

subroutine write_restart_dynamics(File, dyn_out)

    use hycoef,          only: write_restart_hycoef
    use constituents,    only: pcnst
    use dimensions_mod,  only: nlev
    use pio,             only: pio_offset_kind, io_desc_t, pio_double, pio_write_darray
    use time_manager,    only: get_curr_time, get_curr_date

    ! arguments
    type(file_desc_t), intent(inout) :: File
    type(dyn_export_t), intent(in)  :: dyn_out

    ! local variables
    integer(pio_offset_kind), parameter :: t_idx = 1
    type (fv_atmos_type),  pointer :: Atm(:)

    type(io_desc_t),pointer :: iodesc3d,iodesc3d_ns,iodesc3d_ew,iodesc
    integer :: m, ierr
    integer :: array_lens_3d(3), array_lens_2d(2)
    integer :: file_lens_2d(2), file_lens_1d(1)
    integer :: grid_id,grid_id_ns,grid_id_ew
    integer :: grid_dimlens(2),grid_dimlens_ew(2),grid_dimlens_ns(2)
    integer :: ilen,jlen

    !---------------------------------------------------------------------------

    call write_restart_hycoef(File)

    Atm=>dyn_out%atm
    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    grid_id = cam_grid_id('FFSL')
    grid_id_ew = cam_grid_id('FFSL_EW')
    grid_id_ns = cam_grid_id('FFSL_NS')
    
    ! write coordinate variables for unstructured FFSL, NS and EW restart grid 
    ! (restart grids have tile based global indicies with duplicate edge points
    ! being given uniq indicies. All duplicate point written out to restart file)
    ! - io overhead = 6 tile edges are duplicated and read from the file
    ! instead of mpi gathers to fill in duplicates.

    call cam_grid_write_var(File, grid_id)
    call cam_grid_write_var(File, grid_id_ew)
    call cam_grid_write_var(File, grid_id_ns)
    
    ! create map for distributed write
    call cam_grid_dimensions(grid_id, grid_dimlens)
    call cam_grid_dimensions(grid_id_ew, grid_dimlens_ew)
    call cam_grid_dimensions(grid_id_ns, grid_dimlens_ns)

    ilen=ie-is+1
    jlen=je-js+1

    ! create map for distributed write of 2D fields
    array_lens_2d = (/ilen,jlen/)
    file_lens_1d  = (/grid_dimlens(1)/)
    call cam_grid_get_decomp(grid_id, array_lens_2d, file_lens_1d, pio_double, iodesc)
    ! Write PHIS 
    call PIO_Write_Darray(File, phisdesc, iodesc, Atm(mytile)%phis(is:ie,js:je), ierr)
    ! Write PS 
    call PIO_Write_Darray(File, psdesc, iodesc, Atm(mytile)%ps(is:ie,js:je), ierr)

    array_lens_3d = (/ilen,jlen,nlev/)
    file_lens_2d  = (/grid_dimlens(1), nlev/)
    call cam_grid_get_decomp(grid_id, array_lens_3d, file_lens_2d, pio_double, iodesc3d)
    ! Write U a-grid
    call PIO_Write_Darray(File, Udesc, iodesc3d, Atm(mytile)%ua(is:ie,js:je,1:nlev), ierr)
    ! Write V a-grid
    call PIO_Write_Darray(File, Vdesc, iodesc3d, Atm(mytile)%va(is:ie,js:je,1:nlev) , ierr)
    ! Write OMEGA a-grid
    call PIO_Write_Darray(File, Omegadesc, iodesc3d, Atm(mytile)%omga(is:ie,js:je,1:nlev), ierr)
    ! Write DELP a-grid
    call PIO_Write_Darray(File, delpdesc, iodesc3d, Atm(mytile)%delp(is:ie,js:je,1:nlev), ierr)
    ! Write PT a-grid
    call PIO_Write_Darray(File, Tdesc, iodesc3d, Atm(mytile)%pt(is:ie,js:je,1:nlev), ierr)
    ! Write Tracers a-grid
    do m = 1, pcnst
       call PIO_Write_Darray(File, Qdesc(m), iodesc3d, Atm(mytile)%q(is:ie,js:je,1:nlev,m), ierr)
    end do

    deallocate(qdesc)
    
   ! create map for distributed write of 3D NS fields
    array_lens_3d = (/ilen ,(jlen+1), nlev/)
    file_lens_2d  = (/grid_dimlens_ns(1), nlev/)
    call cam_grid_get_decomp(grid_id_ns, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ns)

    !WRITE US
    call PIO_Write_Darray(File, USdesc, iodesc3d_ns, Atm(mytile)%u(is:ie,js:je+1,1:nlev), ierr)

   ! create map for distributed write of 3D EW fields
    array_lens_3d = (/(ilen+1), jlen, nlev /)
    file_lens_2d  = (/grid_dimlens_ew(1), nlev/)
    call cam_grid_get_decomp(grid_id_ew, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ew)

    !WRITE VS
    call PIO_Write_Darray(File, VSdesc, iodesc3d_ew, Atm(mytile)%v(is:ie+1,js:je,1:nlev), ierr)

end subroutine write_restart_dynamics

!=======================================================================

subroutine read_restart_dynamics(File, dyn_in, dyn_out)

    use cam_history_support,    only: max_fieldname_len
    use constituents,  only: cnst_name, pcnst
    use dimensions_mod,only: npy,npx,nlev
    use dyn_comp,      only: dyn_init
    use dyn_grid,      only: Atm
    use mpp_domains_mod, only: mpp_update_domains, DGRID_NE, mpp_get_boundary
    use pio,           only: file_desc_t, pio_double, &
                             pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
                             pio_read_darray, file_desc_t, io_desc_t, pio_double,pio_offset_kind,&
                             pio_seterrorhandling, pio_bcast_error

   ! arguments
   type(File_desc_t), intent(inout) :: File
   type(dyn_import_t), intent(out)  :: dyn_in
   type(dyn_export_t), intent(out)  :: dyn_out

   ! local variables
   integer(pio_offset_kind), parameter :: t_idx = 1

   integer :: tl
   integer :: i, k, m, j
   integer :: ierr, err_handling
   integer :: fnlev
   integer :: ncols_d_ns, ncols_d_ew, ncols_d

   integer :: ncol_d_dimid
   integer :: ncol_d_ns_dimid
   integer :: ncol_d_ew_dimid

   type(var_desc_t) :: omegadesc
   type(var_desc_t) :: delpdesc
   type(var_desc_t) :: udesc
   type(var_desc_t) :: vdesc
   type(var_desc_t) :: usdesc
   type(var_desc_t) :: vsdesc
   type(var_desc_t) :: tdesc
   type(var_desc_t) :: psdesc
   type(var_desc_t) :: phisdesc
   type(var_desc_t), allocatable :: qdesc(:)
   type(io_desc_t),pointer :: iodesc2d, iodesc3d,iodesc3d_ns,iodesc3d_ew
   integer :: array_lens_3d(3), array_lens_2d(2)
   integer :: file_lens_2d(2), file_lens_1d(1)
   integer :: grid_id,grid_id_ns,grid_id_ew,ilen,jlen
   integer :: grid_dimlens(2),grid_dimlens_ns(2),grid_dimlens_ew(2)

   real(r8),    allocatable :: ebuffer(:,:)
   real(r8),    allocatable :: nbuffer(:,:)

   character(len=*), parameter               :: sub = 'read_restart_dynamics'
   character(len=256) :: errormsg
   !----------------------------------------------------------------------------

   ! Note1: the hybrid coefficients are read from the same location as for an
   !        initial run (e.g., dyn_grid_init).

   ! Note2: the dyn_in and dyn_out objects are not associated with the Atm dynamics
   !        object until dyn_init is called.  Until the restart is better integrated
   !        into dyn_init we just access Atm directly from the dyn_grid
   !        module.  FV3 dyn_init calls an fv3 diagnostic init routine that tries to access
   !        surface pressure in the Atm structure and at the top of read_restart PS hasn't 
   !        been read in yet.

   tl = 1

   is = Atm(mytile)%bd%is
   ie = Atm(mytile)%bd%ie
   js = Atm(mytile)%bd%js
   je = Atm(mytile)%bd%je

   call pio_seterrorhandling(File, pio_bcast_error, err_handling)

   ierr = PIO_Inq_DimID(File, 'lev', nlev_dimid)
   ierr = PIO_Inq_dimlen(File, nlev_dimid, fnlev)
   if (nlev /= fnlev) then
      write(errormsg, *) ': Restart file nlev dimension does not match model levels:',&
           'file nlev=',fnlev,', model nlev=',nlev
      call endrun(sub//trim(errormsg))
   end if

   ! variable descriptors of required dynamics fields
   ierr = PIO_Inq_varid(File, 'DELP',     delpdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find DELP')
   ierr = PIO_Inq_varid(File, 'OMEGA',     omegadesc)
   call cam_pio_handle_error(ierr, sub//': cannot find OMEGA')
   ierr = PIO_Inq_varid(File, 'U',     udesc)
   call cam_pio_handle_error(ierr, sub//': cannot find UA')
   ierr = PIO_Inq_varid(File, 'V',     Vdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find VA')
   ierr = PIO_Inq_varid(File, 'US',     usdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find US')
   ierr = PIO_Inq_varid(File, 'VS',     Vsdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find VS')
   ierr = PIO_Inq_varid(File, 'T',     tdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find T')
   ierr = PIO_Inq_varid(File, 'PS', psdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find PS')
   ierr = PIO_Inq_varid(File, 'PHIS', phisdesc)
   call cam_pio_handle_error(ierr, sub//': cannot find PHIS')
   allocate(qdesc(pcnst))
   do m = 1, pcnst
      ierr = PIO_Inq_varid(File, trim(cnst_name(m)), Qdesc(m))
      call cam_pio_handle_error(ierr, sub//': cannot find '//trim(cnst_name(m)))
   end do
   
   ! check whether the restart fields on the GLL grid contain unique columns
   ! or the fv3 task structure (ncol_d_ns = (ie-is+1)*(je-js+2)+npes columns)
   ! or the fv3 task structure (ncol_d_ew = (ie-is+2)*(je-js+1)+npes columns)
   
   ierr = PIO_Inq_DimID(File, 'ncol_d', ncol_d_dimid)
   call cam_pio_handle_error(ierr, sub//': cannot find ncol_d')
   ierr = PIO_Inq_dimlen(File, ncol_d_dimid, ncols_d)
   
   ierr = PIO_Inq_DimID(File, 'ncol_d_ns', ncol_d_ns_dimid)
   call cam_pio_handle_error(ierr, sub//': cannot find ncol_d_ns')
   ierr = PIO_Inq_dimlen(File, ncol_d_ns_dimid, ncols_d_ns)
   
   ierr = PIO_Inq_DimID(File, 'ncol_d_ew', ncol_d_ew_dimid)
   call cam_pio_handle_error(ierr, sub//': cannot find ncol_d_ew')
   ierr = PIO_Inq_dimlen(File, ncol_d_ew_dimid, ncols_d_ew)
   
   grid_id = cam_grid_id('FFSL')
   grid_id_ns = cam_grid_id('FFSL_NS')
   grid_id_ew = cam_grid_id('FFSL_EW')
   call cam_grid_dimensions(grid_id, grid_dimlens)
   call cam_grid_dimensions(grid_id_ew, grid_dimlens_ew)
   call cam_grid_dimensions(grid_id_ns, grid_dimlens_ns)

   if (ncols_d /= grid_dimlens(1)) then
      write(errormsg, *) ':Restart file ncol_d dimension does not match number of model A-Grid columns',&
           'Restart ncols_d=',ncols_d,', A-Grid ncols=',grid_dimlens(1)
      call endrun(sub//trim(errormsg))
   end if

   if (ncols_d_ns /= grid_dimlens_ns(1)) then
      write(errormsg, *) ':Restart file ncol_d dimension does not match number of model D-Grid ns columns',&
           'Restart ncols_d_ns=',ncols_d_ns,', D-Grid ns ncols=',grid_dimlens_ns(1)
      call endrun(sub//trim(errormsg))
   end if

   if (ncols_d_ew /= grid_dimlens_ew(1)) then
      write(errormsg, *) ':Restart file ncol_d dimension does not match number of model D-Grid ew columns',&
           'Restart ncols_d_ew=',ncols_d_ew,', D-Grid ew ncols=',grid_dimlens_ew(1)
      call endrun(sub//trim(errormsg))
   end if

   ilen = ie-is+1
   jlen = je-js+1
   ! create map for distributed write of 2D fields
    array_lens_2d = (/ilen,jlen/)
    file_lens_1d  = (/grid_dimlens(1)/)
    call cam_grid_get_decomp(grid_id, array_lens_2d, file_lens_1d, pio_double, iodesc2d)

   ! create map for distributed write of 3D fields
    array_lens_3d = (/ilen, jlen,nlev/)
    file_lens_2d  = (/grid_dimlens(1), nlev/)
    call cam_grid_get_decomp(grid_id, array_lens_3d, file_lens_2d, pio_double, iodesc3d)
    
   ! create map for distributed write of 3D NS fields
    array_lens_3d = (/ilen, jlen+1, nlev/)
    file_lens_2d  = (/grid_dimlens_ns(1), nlev/)
    call cam_grid_get_decomp(grid_id_ns, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ns)
    
   ! create map for distributed write of 3D EW fields
    array_lens_3d = (/ilen+1, jlen, nlev/)
    file_lens_2d  = (/grid_dimlens_ew(1), nlev/)
    call cam_grid_get_decomp(grid_id_ew, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ew)
    
    ! PS
    call PIO_Read_Darray(File, psdesc, iodesc2d,atm(mytile)%ps(is:ie,js:je), ierr)
    ! PHIS
    call PIO_Read_Darray(File, phisdesc, iodesc2d, atm(mytile)%phis(is:ie,js:je), ierr)
    ! OMEGA
    call PIO_Read_Darray(File, omegadesc, iodesc3d,Atm(mytile)%omga(is:ie,js:je,1:nlev), ierr)
    ! DELP
    call PIO_Read_Darray(File, delpdesc, iodesc3d, atm(mytile)%delp(is:ie,js:je,1:nlev), ierr)
    ! T
    call PIO_Read_Darray(File, Tdesc, iodesc3d,atm(mytile)%pt(is:ie,js:je,1:nlev) , ierr)
    ! V
    call PIO_Read_Darray(File, Vdesc, iodesc3d, atm(mytile)%va(is:ie,js:je,1:nlev), ierr)
    ! U
    call PIO_Read_Darray(File, Udesc, iodesc3d, atm(mytile)%ua(is:ie,js:je,1:nlev), ierr)
    ! tracers
    do m = 1, pcnst
       call PIO_Read_Darray(File, Qdesc(m), iodesc3d, atm(mytile)%q(is:ie,js:je,1:nlev,m), ierr)
    end do

    deallocate(qdesc)

    ! US and VS  After reading unique points on D grid call get_boundary routine to fill
    ! missing points on the north and east block boundaries which are duplicated between
    ! adjacent blocks.

    allocate(ebuffer(npy+2,nlev))
    allocate(nbuffer(npx+2,nlev))
    nbuffer  = 0._r8
    ebuffer  = 0._r8
    ! US
    call PIO_Read_Darray(File, USdesc, iodesc3d_ns, atm(mytile)%u(is:ie,js:je+1,1:nlev), ierr)
    ! VS
    call PIO_Read_Darray(File, VSdesc, iodesc3d_ew, atm(mytile)%v(is:ie+1,js:je,1:nlev), ierr)
    ! US/VS duplicates
    call mpp_get_boundary(atm(mytile)%u, atm(mytile)%v, atm(mytile)%domain, ebuffery=ebuffer,  &
         nbufferx=nbuffer, gridtype=DGRID_NE )
    do k=1,nlev
       do i=is,ie
          atm(mytile)%u(i,je+1,k) = nbuffer(i-is+1,k)
       enddo
       do j=js,je
          atm(mytile)%v(ie+1,j,k) = ebuffer(j-js+1,k)
       enddo
    enddo
    deallocate(ebuffer)
    deallocate(nbuffer)

    ! Update halo points on each processor

    call mpp_update_domains( Atm(mytile)%phis, Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%ps,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%u,atm(mytile)%v,   Atm(mytile)%domain, gridtype=DGRID_NE, complete=.true. )
    call mpp_update_domains( atm(mytile)%pt,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%delp,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%omga,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%q,    Atm(mytile)%domain )

    call dyn_init(dyn_in, dyn_out)

    call pio_seterrorhandling(File, err_handling)
   
   
  end subroutine read_restart_dynamics

end module restart_dynamics
