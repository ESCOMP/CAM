module restart_dynamics

    use shr_kind_mod,  only: r8 => shr_kind_r8, i8 => shr_kind_i8
    use pio,           only: file_desc_t, var_desc_t
    use dyn_comp,      only: dyn_import_t, dyn_export_t
    use fv_arrays_mod, only: fv_atmos_type
    use dyn_grid,      only: get_horiz_grid_dim_d,mytile,mybindex,mybindex_ew,mybindex_ns, &
                             mygindex,mygindex_ew,mygindex_ns, &
                             mylindex,mylindex_ew,mylindex_ns, &
                             mygindexdups,mygindexdups_ew,mygindexdups_ns
    use cam_grid_support, only: cam_grid_header_info_t, cam_grid_id, cam_grid_write_attr, &
                            cam_grid_write_var, cam_grid_get_decomp, cam_grid_dimensions, max_hcoordname_len
    use cam_pio_utils,    only: cam_pio_handle_error
    use cam_abortutils,   only: endrun
    use mpp_mod,           only: mpp_pe, mpp_npes, mpp_chksum
    use spmd_utils,      only: masterproc,iam

    implicit none
    private

    public :: init_restart_dynamics, write_restart_dynamics, read_restart_dynamics

    type(var_desc_t) :: udesc, vdesc, tdesc, psdesc, phisdesc, usdesc,vsdesc,delpdesc,omegadesc

    integer :: ncol_d_dimid, ncol_d_ew_dimid, ncol_d_ns_dimid,  nlev_dimid, nlevp_dimid, npz
    type(var_desc_t), allocatable :: qdesc(:)
    integer(i8) :: checksum
    integer :: is,ie,js,je,isc,iec,jsc,jec,isd,ied,jsd,jed


!=======================================================================
contains
!=======================================================================

subroutine init_restart_dynamics(File, dyn_out)

    use pio,           only: pio_global, pio_unlimited, pio_double, pio_def_dim, &
                             pio_seterrorhandling, pio_bcast_error, pio_noerr, &
                             pio_put_att, pio_def_var, pio_initdecomp, &
                             pio_setdebuglevel, pio_inq_dimid
    use cam_pio_utils, only: pio_subsystem
    use constituents,  only: cnst_name, pcnst
    use hycoef,        only: init_restart_hycoef
    use spmd_utils,    only : npes
    use dyn_grid,      only : uniqpts_glob_ns,uniqpts_glob_ew,uniqpts_glob
    use dimensions_mod,     only: npx, npy, npz
    implicit none

    type(file_desc_t),  intent(inout) :: file
    type(dyn_export_t), intent(in)    :: dyn_out

    integer :: vdimids(2)
    integer :: ierr, i, err_handling
    integer :: time_dimid
    integer :: is,ie,js,je
    type (fv_atmos_type),  pointer :: Atm(:)
    
    integer :: grid_id,grid_id_ns,grid_id_ew,glob_pts_tiles
    type(cam_grid_header_info_t) :: info,info_ew,info_ns

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

    use pio,             only: pio_offset, pio_offset_kind, io_desc_t, pio_double, pio_write_darray, &
                               pio_setframe, pio_put_var, pio_initdecomp, pio_setframe, pio_def_dim, &
                               pio_freedecomp, pio_enddef
    use cam_pio_utils,   only: pio_subsystem
    use hycoef,          only: write_restart_hycoef
    use constituents,    only: pcnst
    use time_manager,    only: get_curr_time, get_curr_date
    use time_manager_mod,only: &
         date_to_string, set_date

    implicit none

    type(file_desc_t) :: File
    type(dyn_export_t), intent(in)  :: dyn_out


    ! local variables
    integer(pio_offset_kind), parameter :: t_idx = 1
    type (fv_atmos_type),  pointer :: Atm(:)

    type(io_desc_t),pointer :: iodesc2d, iodesc3d,iodesc3d_ns,iodesc3d_ew,iodesc
    real(r8), allocatable :: var3d(:,:,:), var3d_ew(:,:,:), var3d_ns(:,:,:), var2d(:,:)
    integer :: i,j,k,m, ierr
    integer :: ncols_d
    integer, pointer :: ldof(:),ldof_ew(:),ldof_ns(:)
    integer :: vsize2d, vsize3d
    integer :: ndcur, nscur
    real(kind=r8) :: time
    integer :: yy             ! CAM current year
    integer :: mm             ! CAM current month
    integer :: dd             ! CAM current day
    integer :: tt             ! CAM current time of day (sec)
    integer, dimension(6) :: date = (/ 0, 0, 0, 0, 0, 0 /)
    character(len=32) :: timestamp   
    integer :: array_lens_3d(3), array_lens_2d(2), array_lens_1d(1)
    integer :: file_lens_2d(2), file_lens_1d(1)
    integer :: grid_id,grid_id_ns,grid_id_ew
    integer :: grid_dimlens(2),grid_dimlens_ew(2),grid_dimlens_ns(2)
    integer :: ilen,jlen

    call write_restart_hycoef(File)

    Atm=>dyn_out%atm
    npz    = Atm(mytile)%npz
    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je
    isc = Atm(mytile)%bd%isc
    iec = Atm(mytile)%bd%iec
    jsc = Atm(mytile)%bd%jsc
    jec = Atm(mytile)%bd%jec

    isd = Atm(mytile)%bd%isd
    ied = Atm(mytile)%bd%ied
    jsd = Atm(mytile)%bd%jsd
    jed = Atm(mytile)%bd%jed


    grid_id = cam_grid_id('FFSL')
    grid_id_ew = cam_grid_id('FFSL_EW')
    grid_id_ns = cam_grid_id('FFSL_NS')
    
    ! write coordinate variables for unstructured FFSL, NS and EW restart grid (restart grids have tile based global indicies with duplicate edge points
    ! being given uniq indicies. All duplicate point written out to restart file) - io overhead = 6 tile edges are duplicated and read from the file
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
    checksum=mpp_chksum(atm(mytile)%phis(is:ie,js:je))
    if (masterproc) write(6,*)'writing PHIS is:ie,js:je CHKSUM=',checksum
    allocate(var2d(ilen,jlen))
    array_lens_2d = (/ilen,jlen/)
    file_lens_1d  = (/grid_dimlens(1)/)
    var2d=Atm(mytile)%phis(is:ie,js:je)
    call cam_grid_get_decomp(grid_id, array_lens_2d, file_lens_1d, pio_double, iodesc)
    call PIO_Write_Darray(File, phisdesc, iodesc, var2d, ierr)

    checksum=mpp_chksum(atm(mytile)%ps(is:ie,js:je))
    if (masterproc) write(6,*)'writing PS is:ie,js:je CHKSUM=',checksum
    var2d=Atm(mytile)%ps(is:ie,js:je)
    call PIO_Write_Darray(File, psdesc, iodesc, var2d, ierr)
    deallocate(var2d)

    checksum=mpp_chksum(atm(mytile)%ua(is:ie,js:je,1:npz))   
    if (masterproc) write(6,*)'writing UA is:ie,js:je CHKSUM=',checksum
    allocate(var3d(ilen,npz,jlen))
    array_lens_3d = (/ilen,npz,jlen/)
    file_lens_2d  = (/grid_dimlens(1), npz/)
    call cam_grid_get_decomp(grid_id, array_lens_3d, file_lens_2d, pio_double, iodesc3d)
    var3d=RESHAPE(Atm(mytile)%ua(is:ie,js:je,1:npz),(/ilen,npz,jlen/),ORDER=(/1,3,2/))
    call PIO_Write_Darray(File, Udesc, iodesc3d, var3d, ierr)

    checksum=mpp_chksum(atm(mytile)%va(is:ie,js:je,1:npz))
    if (masterproc) write(6,*)'writing VA is:ie,js:je CHKSUM=',checksum
    var3d=RESHAPE(Atm(mytile)%va(is:ie,js:je,1:npz),(/ilen,npz,jlen/),ORDER=(/1,3,2/))
    call PIO_Write_Darray(File, Vdesc, iodesc3d, var3d , ierr)

    checksum=mpp_chksum(atm(mytile)%omga(is:ie,js:je,1:npz))
    if (masterproc) write(6,*)'writing omga is:ie,js:je CHKSUM=',checksum
    var3d=RESHAPE(Atm(mytile)%omga(is:ie,js:je,1:npz),(/ilen,npz,jlen/),ORDER=(/1,3,2/))
    call PIO_Write_Darray(File, Omegadesc, iodesc3d, var3d, ierr)

    checksum=mpp_chksum(atm(mytile)%delp(is:ie,js:je,1:npz))
    if (masterproc) write(6,*)'writing delp is:ie,js:je CHKSUM=',checksum
    var3d=RESHAPE(Atm(mytile)%delp(is:ie,js:je,1:npz),(/ilen,npz,jlen/),ORDER=(/1,3,2/))
    call PIO_Write_Darray(File, delpdesc, iodesc3d, var3d, ierr)

    checksum=mpp_chksum(atm(mytile)%pt(is:ie,js:je,1:npz))
    if (masterproc) write(6,*)'writing T is:ie,js:je CHKSUM=',checksum
    var3d=RESHAPE(Atm(mytile)%pt(is:ie,js:je,1:npz),(/ilen,npz,jlen/),ORDER=(/1,3,2/))
    call PIO_Write_Darray(File, Tdesc, iodesc3d, var3d, ierr)

    do m = 1, pcnst
       checksum=mpp_chksum(atm(mytile)%q(is:ie,js:je,1:npz,m))
       if (masterproc) write(6,*)'writing Q is:ie,js:je ',m,' of ',pcnst,' CHKSUM=',checksum
       var3d=RESHAPE(Atm(mytile)%q(is:ie,js:je,1:npz,m),(/ilen,npz,jlen/),ORDER=(/1,3,2/))
       call PIO_Write_Darray(File, Qdesc(m), iodesc3d, var3d, ierr)
    end do

    deallocate(qdesc)
    deallocate(var3d)
    
   ! create map for distributed write of 3D NS fields
    checksum=mpp_chksum(atm(mytile)%u(is:ie,js:je+1,1:npz))
    if (masterproc) write(6,*)'writing US is:ie,js:je+1 CHKSUM=',checksum
    allocate(var3d_ns(ilen,npz,(jlen+1)))
    array_lens_3d = (/ilen , npz, (jlen+1)/)
    file_lens_2d  = (/grid_dimlens_ns(1), npz/)
    call cam_grid_get_decomp(grid_id_ns, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ns)
    var3d_ns=RESHAPE(Atm(mytile)%u(is:ie,js:je+1,1:npz),(/ilen,npz,jlen+1/),ORDER=(/1,3,2/))
    call PIO_Write_Darray(File, USdesc, iodesc3d_ns, var3d_ns, ierr)
    deallocate(var3d_ns)

   ! create map for distributed write of 3D EW fields

    checksum=mpp_chksum(atm(mytile)%v(is:ie+1,js:je,1:npz))
    if (masterproc) write(6,*)'writing VS is:ie+1,js:je CHKSUM=',checksum
    allocate(var3d_ew((ilen+1),npz,jlen))
    array_lens_3d = (/(ilen+1), npz, jlen /)
    file_lens_2d  = (/grid_dimlens_ew(1), npz/)
    call cam_grid_get_decomp(grid_id_ew, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ew)
    var3d_ew=RESHAPE(Atm(mytile)%v(is:ie+1,js:je,1:npz),(/(ilen+1),npz,jlen/),ORDER=(/1,3,2/))
    call PIO_Write_Darray(File, VSdesc, iodesc3d_ew, var3d_ew, ierr)
    deallocate(var3d_ew)

end subroutine write_restart_dynamics

!=======================================================================

subroutine read_restart_dynamics(File, dyn_in, dyn_out)
    use pio,           only: file_desc_t, pio_global, pio_double, pio_offset, &
                             pio_get_att, pio_inq_dimid, pio_inq_dimlen, pio_initdecomp, pio_inq_varid, &
                             pio_read_darray, pio_setframe, file_desc_t, io_desc_t, pio_double,pio_offset_kind,&
                             pio_seterrorhandling, pio_bcast_error, io_desc_t

    use dyn_comp,      only: dyn_init
    use dyn_grid,      only: Atm
    use constituents,  only: cnst_name, pcnst
    use cam_pio_utils, only: pio_subsystem
    use cam_logfile,   only: iulog
    use dimensions_mod,only: npz
    use cam_history_support,    only: max_fieldname_len
    use ncdio_atm,              only: infld
    use mpp_domains_mod, only: mpp_update_domains, domain2D, DGRID_NE

   ! arguments
   type(File_desc_t), intent(inout) :: File
   type(dyn_import_t), intent(out)  :: dyn_in
   type(dyn_export_t), intent(out)  :: dyn_out

   ! local variables
   integer(pio_offset_kind), parameter :: t_idx = 1

   integer :: tl
   integer :: i, k, m, j
   integer :: ierr, err_handling
   integer :: fnlev, fnc
   integer :: hdim_len, ncols_d_ns, ncols_d_ew, ncols_d

   integer :: npz_dimid
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
   type(io_desc_t),pointer :: iodesc2d, iodesc3d,iodesc3d_ns,iodesc3d_ew,iodesc,iodesc3d_ns_rst,iodesc3d_ew_rst
   integer :: array_lens_3d(3), array_lens_2d(2), array_lens_1d(1)
   integer :: file_lens_2d(2), file_lens_1d(1)
   integer :: grid_id,grid_id_ns,grid_id_ew,ilen,jlen,grid_id_ns_rst,grid_id_ew_rst
   integer :: grid_dimlens(2),grid_dimlens_ns(2),grid_dimlens_ew(2),grid_dimlens_ns_rst(2),grid_dimlens_ew_rst(2)
   character(len=max_hcoordname_len) :: dimname1, dimname2

   real(r8),    allocatable :: var2d(:,:)
   real(r8),    allocatable :: var3d(:,:,:)
   real(r8),    allocatable :: var3d_ns(:,:,:)
   real(r8),    allocatable :: var3d_ew(:,:,:)
   real(r8),    allocatable :: var3d_ns_tmp(:,:,:)
   real(r8),    allocatable :: var3d_ew_tmp(:,:,:)
   real(r8),    allocatable :: var3d_ns_rst(:,:,:)
   real(r8),    allocatable :: var3d_ew_rst(:,:,:)

   logical :: found

   character(len=*), parameter               :: sub = 'read_restart_dynamics'
   character(len=max_fieldname_len)          :: fieldname
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

   npz    = Atm(mytile)%npz
   is = Atm(mytile)%bd%is
   ie = Atm(mytile)%bd%ie
   js = Atm(mytile)%bd%js
   je = Atm(mytile)%bd%je
   isc = Atm(mytile)%bd%isc
   iec = Atm(mytile)%bd%iec
   jsc = Atm(mytile)%bd%jsc
   jec = Atm(mytile)%bd%jec
   isd = Atm(mytile)%bd%isd
   ied = Atm(mytile)%bd%ied
   jsd = Atm(mytile)%bd%jsd
   jed = Atm(mytile)%bd%jed

   call pio_seterrorhandling(File, pio_bcast_error, err_handling)

   ierr = PIO_Inq_DimID(File, 'lev', nlev_dimid)
   ierr = PIO_Inq_dimlen(File, nlev_dimid, fnlev)
   if (npz /= fnlev) then
      write(iulog,*) 'Restart file nlev does not match model. nlev (file, namelist):', &
                     fnlev, npz
      call endrun(sub//': Restart file nlev does not match model.')
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
   ! variable descriptors of required dynamics fields
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
   
   
   ! ierr = PIO_Inq_DimID(File, 'ncol_d_ns_rst', ncol_d_ns_rst_dimid)
   ! call cam_pio_handle_error(ierr, sub//': cannot find ncol_d_ns_rst')
   ! ierr = PIO_Inq_dimlen(File, ncol_d_ns_rst_dimid, ncols_d_ns_rst)
   
   ! ierr = PIO_Inq_DimID(File, 'ncol_d_ew_rst', ncol_d_ew_rst_dimid)
   ! call cam_pio_handle_error(ierr, sub//': cannot find ncol_d_ew_rst')
   ! ierr = PIO_Inq_dimlen(File, ncol_d_ew_rst_dimid, ncols_d_ew_rst)

   ierr = PIO_Inq_DimID(File, 'ncol_d_ns', ncol_d_ns_dimid)
   call cam_pio_handle_error(ierr, sub//': cannot find ncol_d_ns')
   ierr = PIO_Inq_dimlen(File, ncol_d_ns_dimid, ncols_d_ns)
   
   ierr = PIO_Inq_DimID(File, 'ncol_d_ew', ncol_d_ew_dimid)
   call cam_pio_handle_error(ierr, sub//': cannot find ncol_d_ew')
   ierr = PIO_Inq_dimlen(File, ncol_d_ew_dimid, ncols_d_ew)
   
   grid_id = cam_grid_id('FFSL')
   grid_id_ns = cam_grid_id('FFSL_NS')
   grid_id_ew = cam_grid_id('FFSL_EW')
   grid_id_ns_rst = cam_grid_id('FFSL_NS_RST')
   grid_id_ew_rst = cam_grid_id('FFSL_EW_RST')
   call cam_grid_dimensions(grid_id, grid_dimlens)
   call cam_grid_dimensions(grid_id_ew, grid_dimlens_ew)
   call cam_grid_dimensions(grid_id_ns, grid_dimlens_ns)
   call cam_grid_dimensions(grid_id_ew_rst, grid_dimlens_ew_rst)
   call cam_grid_dimensions(grid_id_ns_rst, grid_dimlens_ns_rst)
   
   if (ncols_d_ns /= grid_dimlens_ns(1)) then
      write(iulog,*) 'Restart file ncol_d_ns does not match model. ncols_d_ns (file, model):',&
           ncols_d_ns, grid_dimlens_ns(1)
      call endrun(sub//': Restart file ncols_d_ns does not match model.')
   end if
   
   ilen = ie-is+1
   jlen = je-js+1
   ! create map for distributed write of 2D fields
    array_lens_2d = (/ilen,jlen/)
    file_lens_1d  = (/grid_dimlens(1)/)
    call cam_grid_get_decomp(grid_id, array_lens_2d, file_lens_1d, pio_double, iodesc2d)

   ! create map for distributed write of 3D fields
    array_lens_3d = (/ilen,npz, jlen/)
    file_lens_2d  = (/grid_dimlens(1), npz/)
    call cam_grid_get_decomp(grid_id, array_lens_3d, file_lens_2d, pio_double, iodesc3d)
    
   ! create map for distributed write of 3D NS fields
    array_lens_3d = (/ilen, npz, jlen+1/)
    file_lens_2d  = (/grid_dimlens_ns(1), npz/)
    call cam_grid_get_decomp(grid_id_ns, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ns)
    
   ! create map for distributed write of 3D NS RST fields (reading dups - map has gindex for dups only all others 0)
    array_lens_3d = (/ilen, npz, jlen+1/)
    file_lens_2d  = (/grid_dimlens_ns_rst(1), npz/)
    call cam_grid_get_decomp(grid_id_ns_rst, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ns_rst)

   ! create map for distributed write of 3D EW fields
    array_lens_3d = (/ilen+1, npz, jlen/)
    file_lens_2d  = (/grid_dimlens_ew(1), npz/)
    call cam_grid_get_decomp(grid_id_ew, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ew)
    
   ! create map for distributed write of 3D EW RST fields  (reading dups - map has gindex for dups only all others 0)
    array_lens_3d = (/ilen+1, npz, jlen/)
    file_lens_2d  = (/grid_dimlens_ew_rst(1), npz/)
    call cam_grid_get_decomp(grid_id_ew_rst, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ew_rst)
    
    allocate(var2d(is:ie,js:je))
    var2d = 0._r8
    ! PS
    call PIO_Read_Darray(File, psdesc, iodesc2d, var2d, ierr)
    atm(mytile)%ps(is:ie,js:je) = var2d
    ! PHIS
    call PIO_Read_Darray(File, phisdesc, iodesc2d, var2d, ierr)
    atm(mytile)%phis(is:ie,js:je) = var2d
    deallocate(var2d)

    
    allocate(var3d(is:ie,npz,js:je))
    var3d = 0._r8

    ! OMEGA
    call PIO_Read_Darray(File, omegadesc, iodesc3d, var3d, ierr)
    Atm(mytile)%omga(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))

    ! DELP
    call PIO_Read_Darray(File, delpdesc, iodesc3d, var3d, ierr)
    atm(mytile)%delp(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))

    ! T
    call PIO_Read_Darray(File, Tdesc, iodesc3d, var3d, ierr)
    atm(mytile)%pt(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))

    ! V
    call PIO_Read_Darray(File, Vdesc, iodesc3d, var3d, ierr)
    atm(mytile)%va(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))

    ! U
    call PIO_Read_Darray(File, Udesc, iodesc3d, var3d, ierr)
    atm(mytile)%ua(is:ie,js:je,1:npz)   =RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
    ! tracers
    do m = 1, pcnst
       call PIO_Read_Darray(File, Qdesc(m), iodesc3d, var3d, ierr)
       atm(mytile)%q(is:ie,js:je,1:npz,m) = RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
    end do

    deallocate(var3d)
    deallocate(qdesc)

    ! US - For ew and ns grids PIO 1 can not read a single file point to
    !      multiple processors.  For now I am reading twice, once to get the
    !      uniq points distributed and the second time to read those points
    !      required on more than one processor. (ie those points at the 
    !      n/s e/w DGRID boundaries between processors

    allocate(var3d_ns(is:ie, npz,js:je+1))
    allocate(var3d_ns_tmp(is:ie, npz,js:je+1))
    allocate(var3d_ns_rst(is:ie, npz, js:je+1))
    var3d_ns = 0._r8
    var3d_ns_tmp = 0._r8
    var3d_ns_rst = 0._r8
    call PIO_Read_Darray(File, USdesc, iodesc3d_ns, var3d_ns_tmp, ierr)
    var3d_ns=var3d_ns_tmp
    ! US hack to read in duplicate points on adjacent processor
    ! should fill in zeros in the decomposed global array from previous read
    call PIO_Read_Darray(File, USdesc, iodesc3d_ns_rst, var3d_ns_rst, ierr)
    where(var3d_ns_tmp.eq.0)
       var3d_ns = var3d_ns_rst
    end where
    atm(mytile)%u(is:ie,js:je+1,1:npz) = RESHAPE(var3d_ns,(/ilen,jlen+1,npz/),ORDER=(/1,3,2/))

    deallocate(var3d_ns)
    deallocate(var3d_ns_rst)
    deallocate(var3d_ns_tmp)

    allocate(var3d_ew(is:ie+1, npz, js:je))
    allocate(var3d_ew_tmp(is:ie+1, npz, js:je))
    allocate(var3d_ew_rst(is:ie+1, npz, js:je))
    var3d_ew = 0._r8
    var3d_ew_tmp = 0._r8
    var3d_ew_rst = 0._r8

    ! VS
    call PIO_Read_Darray(File, VSdesc, iodesc3d_ew, var3d_ew_tmp, ierr)
    var3d_ew=var3d_ew_tmp


    ! VS hack to read in duplicate points on adjacent processor
    ! should fill in zeros in the decomposed global array from previous read
    call PIO_Read_Darray(File, VSdesc, iodesc3d_ew_rst, var3d_ew_rst, ierr)
    where(var3d_ew_tmp.eq.0)
       var3d_ew = var3d_ew_rst
    end where
    atm(mytile)%v(is:ie+1,js:je,1:npz) = RESHAPE(var3d_ew,(/ilen+1,jlen,npz/),ORDER=(/1,3,2/))

    deallocate(var3d_ew)
    deallocate(var3d_ew_tmp)
    deallocate(var3d_ew_rst)

    ! Update halo points on each processor and print out checksums

    call mpp_update_domains( Atm(mytile)%phis, Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%ps,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%u,atm(mytile)%v,   Atm(mytile)%domain, gridtype=DGRID_NE, complete=.true. )
    call mpp_update_domains( atm(mytile)%pt,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%delp,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%omga,   Atm(mytile)%domain )
    call mpp_update_domains( atm(mytile)%q,    Atm(mytile)%domain )

    checksum=mpp_chksum(atm(mytile)%phis(is:ie,js:je))
    if (masterproc) write(6,*)'reading PHIS is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%ps(is:ie,js:je))
    if (masterproc) write(6,*)'reading PS is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%ua(is:ie,js:je,:))
    if (masterproc) write(6,*)'reading UA is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%va(is:ie,js:je,:))
    if (masterproc) write(6,*)'reading VA is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%u(is:ie,js:je+1,:))
    if (masterproc) write(6,*)'reading US is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%v(is:ie+1,js:je,:))
    if (masterproc) write(6,*)'reading VS is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%pt(is:ie,js:je,:))
    if (masterproc) write(6,*)'reading T is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%delp(is:ie,js:je,:))
    if (masterproc) write(6,*)'reading delp is:ie,js:je CHKSUM=',checksum
    checksum=mpp_chksum(atm(mytile)%omga(is:ie,js:je,:))
    if (masterproc) write(6,*)'reading omga is:ie,js:je CHKSUM=',checksum
    do m = 1, pcnst
       checksum=mpp_chksum(atm(mytile)%q(is:ie,js:je,:,m))
       if (masterproc) write(6,*)'reading Q is:ie,js:je ',m,' of ',pcnst,' CHKSUM=',checksum
    end do

    call dyn_init(dyn_in, dyn_out)

    call pio_seterrorhandling(File, err_handling)
   
   
  end subroutine read_restart_dynamics

end module restart_dynamics
