module atm_stream_ndep

  !-----------------------------------------------------------------------
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic ndep file handling and
  ! interpolation.
  !-----------------------------------------------------------------------
  !
  use ESMF
  use dshr_strdata_mod , only : shr_strdata_type
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_sys_mod      , only : shr_sys_abort
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdMod          , only : mpicom, masterproc, iam
  use cam_logfile      , only : iulog
  use cam_abortutils   , only : endrun

  implicit none
  private

  public :: stream_ndep_init      ! position datasets for dynamic ndep
  public :: stream_ndep_interp    ! interpolates between two years of ndep file data

  private :: stream_ndep_check_units   ! Check the units and make sure they can be used

  type(shr_strdata_type) :: sdat_ndep                      ! input data stream
  character(len=CS)      :: stream_varlist_ndep(2)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine stream_ndep_init(model_mesh, model_clock, rc)
    !
    ! Initialize data stream information.

    ! Uses:
    use shr_nl_mod       , only : shr_nl_find_group_name
    use shr_mpi_mod      , only : shr_mpi_bcast
    use dshr_strdata_mod , only : shr_strdata_init_from_inline

    ! input/output variables
    type(ESMF_Mesh)  , intent(in)  :: model_mesh
    type(ESMF_Clock) , intent(in)  :: model_clock
    integer          , intent(out) :: rc

    ! local variables
    integer                 :: nu_nml                 ! unit for namelist file
    integer                 :: nml_error              ! namelist i/o error flag
    character(len=CL)       :: stream_datafile_ndep
    character(len=CL)       :: stream_meshfile_ndep
    integer                 :: stream_year_first_ndep ! first year in stream to use
    integer                 :: stream_year_last_ndep  ! last year in stream to use
    integer                 :: model_year_align_ndep  ! align stream_year_firstndep with
    integer                 :: stream_nflds
    character(*), parameter :: subName = "('stream_ndep_init')"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    namelist /ndep_stream_nml/    &
         stream_datafile_ndep,    &
         stream_meshfile_ndep,    &
         stream_year_first_ndep,  &
         stream_year_last_ndep,   &
         model_year_align_ndep

    ! Default values for namelist
    stream_datafile_ndep   = ' '
    stream_meshfile_ndep   = ' '
    stream_year_first_ndep = 1                ! first year in stream to use
    stream_year_last_ndep  = 1                ! last  year in stream to use
    model_year_align_ndep  = 1                ! align stream_year_first_ndep with this model year

    ! For now variable list in stream data file is hard-wired
    stream_varlist_ndep = (/'NDEP_Nhx', 'NDEP_Noy'/)

    ! Read ndepdyn_nml namelist
    if (masterproc) then
       open( newunit=nu_nml, file='atm_in', status='old', iostat=nml_error )
       call shr_nl_find_group_name(nu_nml, 'ndepdyn_nml', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=ndep_stream_nml, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(' ERROR reading ndepdyn_nml namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(' ERROR finding ndepdyn_nml namelist'//errMsg(sourcefile, __LINE__))
       end if
       close(nu_nml)
    endif
    call shr_mpi_bcast(stream_meshfile_ndep   , mpicom)
    call shr_mpi_bcast(stream_datafile_ndep   , mpicom)
    call shr_mpi_bcast(stream_year_first_ndep , mpicom)
    call shr_mpi_bcast(stream_year_last_ndep  , mpicom)
    call shr_mpi_bcast(model_year_align_ndep  , mpicom)

    if (masterproc) then
       write(iulog,'(a)'   ) ' '
       write(iulog,'(a,i8)') 'stream ndep settings:'
       write(iulog,'(a,a)' ) '  stream_datafile_ndep    = ',trim(stream_datafile_ndep)
       write(iulog,'(a,a)' ) '  stream_meshfile_ndep    = ',trim(stream_meshfile_ndep)
       write(iulog,'(a,a)' ) '  stream_varlist_ndep     = ',trim(stream_varlist_ndep(1))
       write(iulog,'(a,i8)') '  stream_year_first_ndep  = ',stream_year_first_ndep
       write(iulog,'(a,i8)') '  stream_year_last_ndep   = ',stream_year_last_ndep
       write(iulog,'(a,i8)') '  model_year_align_ndep   = ',model_year_align_ndep
       write(iulog,'(a)'   ) ' '
    endif

    ! Read in units
    call stream_ndep_check_units(stream_datafile_ndep)

    ! Initialize the cdeps data type sdat_ndep
    call shr_strdata_init_from_inline(sdat_ndep,               &
         my_task             = iam,                            &
         logunit             = iulog,                          &
         compname            = 'ATM',                          &
         model_clock         = model_clock,                    &
         model_mesh          = model_mesh,                     &
         stream_meshfile     = trim(stream_meshfile_ndep),     &
         stream_filenames    = (/trim(stream_datafile_ndep)/), &
         stream_yearFirst    = stream_year_first_ndep,         &
         stream_yearLast     = stream_year_last_ndep,          &
         stream_yearAlign    = model_year_align_ndep,          &
         stream_fldlistFile  = stream_varlist_ndep,            &
         stream_fldListModel = stream_varlist_ndep,            &
         stream_lev_dimname  = 'null',                         &
         stream_mapalgo      = 'bilinear',                     &
         stream_offset       = 0,                              &
         stream_taxmode      = 'cycle',                        &
         stream_dtlimit      = 1.5_r8,                         &
         stream_tintalgo     = 'linear',                       &
         stream_name         = 'Nitrogen deposition data ',    &
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

    use ncdio_pio     , only : ncd_pio_openfile, ncd_inqvid, ncd_getatt, ncd_pio_closefile, ncd_nowrite
    use ncdio_pio     , only : file_desc_t, var_desc_t
    use shr_log_mod   , only : errMsg => shr_log_errMsg

    ! Arguments
    character(len=*), intent(in)  :: stream_fldFileName_ndep  ! ndep filename
    !
    ! Local variables
    type(file_desc_t) :: ncid     ! NetCDF filehandle for ndep file
    type(var_desc_t)  :: vardesc  ! variable descriptor
    integer           :: varid    ! variable index
    logical           :: readvar  ! If variable was read
    character(len=CS) :: ndepunits! ndep units
    !-----------------------------------------------------------------------

    call ncd_pio_openfile( ncid, trim(stream_fldFileName_ndep), ncd_nowrite )
    call ncd_inqvid(ncid, stream_varlist_ndep(1), varid, vardesc, readvar=readvar)
    if ( readvar ) then
       call ncd_getatt(ncid, varid, "units", ndepunits)
    else
       call endrun(' ERROR finding variable: '//trim(stream_varlist_ndep(1))//" in file: "// &
            trim(stream_fldFileName_ndep)//errMsg(sourcefile, __LINE__))
    end if
    call ncd_pio_closefile( ncid )

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

    ! Get pointer for stream data that is time and spatially interpolate to model time and grid
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
