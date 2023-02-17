module atm_stream_ndep

  !-----------------------------------------------------------------------
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic ndep file handling and
  ! interpolation.
  !-----------------------------------------------------------------------
  !
  use ESMF              , only : ESMF_Clock, ESMF_Mesh
  use ESMF              , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF              , only : ESMF_Finalize, ESMF_LogFoundError
  use nuopc_shr_methods , only : chkerr
  use dshr_strdata_mod  , only : shr_strdata_type
  use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use spmd_utils        , only : mpicom, masterproc, iam
  use spmd_utils        , only : mpi_character, mpi_integer
  use cam_logfile       , only : iulog
  use cam_abortutils    , only : endrun

  implicit none
  private

  public :: stream_ndep_init      ! position datasets for dynamic ndep
  public :: stream_ndep_interp    ! interpolates between two years of ndep file data

  private :: stream_ndep_check_units   ! Check the units and make sure they can be used

  type(shr_strdata_type) :: sdat_ndep                      ! input data stream
  logical, public        :: stream_ndep_is_initialized = .false.
  character(len=CS)      :: stream_varlist_ndep(2)
  type(ESMF_Clock)       :: model_clock

  character(len=*), parameter :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine stream_ndep_init(model_mesh, model_clock, rc)
    !
    ! Initialize data stream information.

    ! Uses:
    use shr_nl_mod       , only : shr_nl_find_group_name
    use dshr_strdata_mod , only : shr_strdata_init_from_inline

    ! input/output variables
    type(ESMF_CLock), intent(in)  :: model_clock
    type(ESMF_Mesh) , intent(in)  :: model_mesh
    integer         , intent(out) :: rc

    ! local variables
    integer                 :: nu_nml                 ! unit for namelist file
    integer                 :: nml_error              ! namelist i/o error flag
    character(len=CL)       :: stream_ndep_data_filename
    character(len=CL)       :: stream_ndep_mesh_filename
    integer                 :: stream_ndep_year_first ! first year in stream to use
    integer                 :: stream_ndep_year_last  ! last year in stream to use
    integer                 :: stream_ndep_year_align ! align stream_year_firstndep with
    integer                 :: ierr
    character(*), parameter :: subName = "('stream_ndep_init')"
    !-----------------------------------------------------------------------

    namelist /ndep_stream_nl/      &
         stream_ndep_data_filename, &
         stream_ndep_mesh_filename, &
         stream_ndep_year_first,    &
         stream_ndep_year_last,     &
         stream_ndep_year_align

    rc = ESMF_SUCCESS

    ! Default values for namelist
    stream_ndep_data_filename = ' '
    stream_ndep_mesh_filename = ' '
    stream_ndep_year_first    = 1 ! first year in stream to use
    stream_ndep_year_last     = 1 ! last  year in stream to use
    stream_ndep_year_align    = 1 ! align stream_ndep_year_first with this model year

    ! For now variable list in stream data file is hard-wired
    stream_varlist_ndep = (/'NDEP_NHx_month', 'NDEP_NOy_month'/)

    ! Read ndep_stream namelist
    if (masterproc) then
       open( newunit=nu_nml, file='atm_in', status='old', iostat=nml_error )
       call shr_nl_find_group_name(nu_nml, 'ndep_stream_nl', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=ndep_stream_nl, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(' ERROR reading ndep_stream_nl namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(' ERROR finding ndep_stream_nl namelist'//errMsg(sourcefile, __LINE__))
       end if
       close(nu_nml)
    endif
    call mpi_bcast(stream_ndep_mesh_filename, len(stream_ndep_mesh_filename), mpi_character, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_mesh_filename")
    call mpi_bcast(stream_ndep_data_filename, len(stream_ndep_data_filename), mpi_character, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_data_filename")
    call mpi_bcast(stream_ndep_year_first, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_year_first")
    call mpi_bcast(stream_ndep_year_last, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_year_last")
    call mpi_bcast(stream_ndep_year_align, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_year_align")

    if (masterproc) then
       write(iulog,'(a)'   ) ' '
       write(iulog,'(a,i8)')  'stream ndep settings:'
       write(iulog,'(a,a)' )  '  stream_ndep_data_filename = ',trim(stream_ndep_data_filename)
       write(iulog,'(a,a)' )  '  stream_ndep_mesh_filename = ',trim(stream_ndep_mesh_filename)
       write(iulog,'(a,a,a)') '  stream_varlist_ndep       = ',trim(stream_varlist_ndep(1)), trim(stream_varlist_ndep(2))
       write(iulog,'(a,i8)')  '  stream_ndep_year_first    = ',stream_ndep_year_first
       write(iulog,'(a,i8)')  '  stream_ndep_year_last     = ',stream_ndep_year_last
       write(iulog,'(a,i8)')  '  stream_ndep_year_align    = ',stream_ndep_year_align
       write(iulog,'(a)'   )  ' '
    endif

    ! Read in units
    call stream_ndep_check_units(stream_ndep_data_filename)

    ! Initialize the cdeps data type sdat_ndep
    call shr_strdata_init_from_inline(sdat_ndep,                    &
         my_task             = iam,                                 &
         logunit             = iulog,                               &
         compname            = 'ATM',                               &
         model_clock         = model_clock,                         &
         model_mesh          = model_mesh,                          &
         stream_meshfile     = trim(stream_ndep_mesh_filename),     &
         stream_filenames    = (/trim(stream_ndep_data_filename)/), &
         stream_yearFirst    = stream_ndep_year_first,              &
         stream_yearLast     = stream_ndep_year_last,               &
         stream_yearAlign    = stream_ndep_year_align,              &
         stream_fldlistFile  = stream_varlist_ndep,                 &
         stream_fldListModel = stream_varlist_ndep,                 &
         stream_lev_dimname  = 'null',                              &
         stream_mapalgo      = 'bilinear',                          &
         stream_offset       = 0,                                   &
         stream_taxmode      = 'cycle',                             &
         stream_dtlimit      = 1.0e30_r8,                           &
         stream_tintalgo     = 'linear',                            &
         stream_name         = 'Nitrogen deposition data ',         &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine stream_ndep_init

  !================================================================
  subroutine stream_ndep_check_units( stream_fldFileName_ndep)

    !--------------------------------------------------------
    ! Check that units are correct on the file and if need any conversion
    !--------------------------------------------------------

     use cam_pio_utils , only : cam_pio_createfile, cam_pio_openfile, cam_pio_closefile, pio_subsystem
     use pio           , only : file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim
     use pio           , only : pio_bcast_error, pio_seterrorhandling, pio_inq_varid, pio_get_att
     use pio           , only : PIO_NOERR, PIO_NOWRITE

    ! Arguments
    character(len=*), intent(in)  :: stream_fldFileName_ndep  ! ndep filename
    !
    ! Local variables
    type(file_desc_t) :: File     ! NetCDF filehandle for ndep file
    type(var_desc_t)  :: vardesc  ! variable descriptor
    integer           :: ierr     ! error status
    integer           :: err_handling ! temporary
    character(len=CS) :: ndepunits! ndep units
    !-----------------------------------------------------------------------

    call cam_pio_openfile( File, trim(stream_fldFileName_ndep), PIO_NOWRITE)
    call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
    ierr = pio_inq_varid(File, stream_varlist_ndep(1), vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(' ERROR finding variable: '//trim(stream_varlist_ndep(1))//" in file: "// &
            trim(stream_fldFileName_ndep)//errMsg(sourcefile, __LINE__))
    else
       ierr = PIO_get_att(File, vardesc, "units", ndepunits)
    end if
    call pio_seterrorhandling(File, err_handling)
    call cam_pio_closefile(File)

    ! Now check to make sure they are correct
    if (.not. trim(ndepunits) == "g(N)/m2/s"  )then
       call endrun(' ERROR in units for nitrogen deposition equal to: '//trim(ndepunits)//" not units expected"// &
            errMsg(sourcefile, __LINE__))
    end if

  end subroutine stream_ndep_check_units

  !================================================================
  subroutine stream_ndep_interp(cam_out, rc)

    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use camsrfexch       , only : cam_out_t
    use ppgrid           , only : begchunk, endchunk
    use time_manager     , only : get_curr_date
    use phys_grid        , only : get_ncols_p

    ! input/output variables
    type(cam_out_t) , intent(inout)  :: cam_out(begchunk:endchunk)
    integer         , intent(out)    :: rc

    ! local variables
    integer :: i,c,g
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    real(r8), pointer :: dataptr1d_nhx(:)
    real(r8), pointer :: dataptr1d_noy(:)
    !-----------------------------------------------------------------------

    ! Advance sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(sdat_ndep, ymd=mcdate, tod=sec, logunit=iulog, istr='ndepdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, stream_varlist_ndep(1), fldptr1=dataptr1d_nhx, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, stream_varlist_ndep(2), fldptr1=dataptr1d_noy, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    g = 1
    do c = begchunk,endchunk
       do i = 1,get_ncols_p(c)
          cam_out(c)%nhx_nitrogen_flx(i) = dataptr1d_nhx(g)
          cam_out(c)%noy_nitrogen_flx(i) = dataptr1d_noy(g)
          g = g + 1
       end do
    end do

  end subroutine stream_ndep_interp

end module atm_stream_ndep
