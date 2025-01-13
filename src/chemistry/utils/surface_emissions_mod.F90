module surface_emissions_mod
  !---------------------------------------------------------------
  ! 	... surface emissions module
  !---------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8, shr_kind_cl
  use spmd_utils,    only : masterproc
  use cam_abortutils,only : endrun
  use ioFileMod,     only : getfil
  use cam_logfile,   only : iulog
  use tracer_data,   only : trfld,trfile
  use infnan,        only : nan, assignment(=)
  use cam_history,   only : addfld, outfld, add_default, horiz_only, fieldname_len

  implicit none

  type :: emission
     integer           :: bufndx
     real(r8)          :: scalefactor
     character(len=256):: filename
     character(len=16) :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type emission

  private

  public :: surface_emissions_readnl
  public :: surface_emissions_reg
  public :: surface_emissions_init
  public :: surface_emissions_adv
  public :: surface_emissions_set

  integer, parameter :: NMAX=50

  type(emission), allocatable :: emissions(:)
  integer                     :: n_emis_files = 0
  integer                     :: n_pbuf_flds = 0

  character(len=shr_kind_cl) :: emissions_specifier(NMAX) = ' '
  character(len=24) :: emissions_type
  integer :: emissions_cycle_yr
  integer :: emissions_fixed_ymd
  integer :: emissions_fixed_tod

  character(len=fieldname_len) :: names(NMAX) = ' '
  character(len=32) :: units(NMAX) = ' '
  integer :: indexes(NMAX) = -1
  integer :: n_diags = 0

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine surface_emissions_readnl(nlfile)

    use namelist_utils, only : find_group_name
    use spmd_utils,     only : mpicom, masterprocid, mpi_integer, mpi_character

    character(len=*), intent(in) :: nlfile ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'surface_emissions_readnl'

    namelist /surface_emissions_opts/ emissions_specifier, emissions_type, emissions_cycle_yr, &
                                       emissions_fixed_ymd, emissions_fixed_tod

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'surface_emissions_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, surface_emissions_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(emissions_specifier,len(emissions_specifier(1))*NMAX, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(emissions_type, len(emissions_type), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(emissions_cycle_yr, 1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(emissions_fixed_ymd, 1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(emissions_fixed_tod, 1, mpi_integer, masterprocid, mpicom, ierr)

  end subroutine surface_emissions_readnl

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine surface_emissions_reg( )
    use m_MergeSorts,   only : IndexSort
    use physics_buffer, only : pbuf_add_field, dtype_r8, pbuf_get_index
    use ppgrid,         only : pcols

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer  :: astat
    integer  :: j, l, m, n, i, nn, kk                   ! Indices
    character(len=16)  :: spc_name
    character(len=256) :: filename

    character(len=16)  :: emis_species(NMAX)
    character(len=256) :: emis_filenam(NMAX)
    integer  :: emis_indexes(NMAX)
    integer  :: indx(NMAX)
    real(r8) :: emis_scalefactor(NMAX)

    character(len=256) :: tmp_string = ' '
    character(len=32) :: xchr = ' '
    real(r8) :: xdbl


    integer :: err
    character(len=32) :: bname

    kk = 0
    nn = 0
    indx(:) = 0
    emis_species = ' '
    emis_indexes = -1
    emis_filenam = 'NONE'

    count_emis: do n=1,size(emissions_specifier)
       if ( len_trim(emissions_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(emissions_specifier(n),'->')
       spc_name = trim(adjustl(emissions_specifier(n)(:i-1)))

       ! need to parse out scalefactor ...
       tmp_string = adjustl(emissions_specifier(n)(i+2:))
       j = scan( tmp_string, '*' )
       if (j>0) then
          xchr = tmp_string(1:j-1) ! get the multipler (left of the '*')
          read( xchr, * ) xdbl   ! convert the string to a real
          tmp_string = adjustl(tmp_string(j+1:)) ! get the filepath name (right of the '*')
       else
          xdbl = 1._r8
       endif
       filename = trim(tmp_string)

       bname = trim(spc_name)//'_srfemis'

       m = pbuf_get_index(bname,errcode=err)
       if (m<1) then
          call pbuf_add_field(bname, 'physpkg', dtype_r8, (/pcols/), m)
          kk = kk+1
          names(kk) = bname
          indexes(kk) = m
       endif

       nn = nn+1
       emis_species(nn) = spc_name
       emis_filenam(nn) = filename
       emis_indexes(nn) = m
       emis_scalefactor(nn) = xdbl

       indx(n)=n
    enddo count_emis

    n_diags = kk
    n_emis_files = nn

    if (masterproc) write(iulog,*) 'srf_emis_inti: n_emis_files = ',n_emis_files

    allocate( emissions(n_emis_files), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'srf_emis_inti: failed to allocate emissions array; error = ',astat
       call endrun('srf_emis_inti: failed to allocate emissions array')
    end if

    !-----------------------------------------------------------------------
    ! Sort the input files so that the emissions sources are summed in the
    ! same order regardless of the order of the input files in the namelist
    !-----------------------------------------------------------------------
    if (n_emis_files > 0) then
       call IndexSort(n_emis_files, indx, emis_filenam)
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the emission type array
    !-----------------------------------------------------------------------
    do m=1,n_emis_files
       emissions(m)%bufndx           = emis_indexes(indx(m))
       emissions(m)%species          = emis_species(indx(m))
       emissions(m)%filename         = emis_filenam(indx(m))
       emissions(m)%scalefactor      = emis_scalefactor(indx(m))
    enddo
  end subroutine surface_emissions_reg

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine surface_emissions_init(pbuf2d)
    use tracer_data,   only : trcdata_init
    use cam_pio_utils, only : cam_pio_openfile
    use pio,           only : pio_inquire, pio_nowrite, pio_closefile, pio_inq_varndims
    use pio,           only : pio_inq_varname, pio_inq_vardimid, pio_inq_dimid
    use pio,           only : file_desc_t, pio_get_att, PIO_NOERR, PIO_GLOBAL
    use pio,           only : pio_seterrorhandling, PIO_BCAST_ERROR,PIO_INTERNAL_ERROR
    use string_utils,  only : GLC
    use physics_buffer,only : physics_buffer_desc, pbuf_set_field

    !--------------------------------------------------------
    !	... Dummy arguments
    !--------------------------------------------------------
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: ierr, astat, l, m, n
    logical :: unstructured
    integer :: vid, nvars, isec, num_dims_emis
    integer :: vndims
    logical, allocatable :: is_sector(:)
    type(file_desc_t) :: ncid
    character(len=32)  :: varname
    character(len=256) :: locfn
    character(len=80) :: file_interp_type = ' '
    integer, allocatable :: dimids(:)
    integer :: time_dimid, ncol_dimid

    character(len=32) :: emis_type = ' '
    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.
    real(r8) :: xnan

    xnan = nan
    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    files_loop: do m = 1, n_emis_files

       emissions(m)%nsectors = 0
       call getfil (emissions(m)%filename, locfn, 0)
       call cam_pio_openfile ( ncid, trim(locfn), PIO_NOWRITE)
       ierr = pio_inquire (ncid, nVariables=nvars)

       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ierr = pio_inq_dimid( ncid, 'ncol', ncol_dimid )
       unstructured = ierr==PIO_NOERR
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

       allocate(is_sector(nvars))
       is_sector(:) = .false.

       if (unstructured) then
          ierr = pio_inq_dimid( ncid, 'time', time_dimid )
       end if

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, vndims)

          if (unstructured) then
             num_dims_emis = 2
          else
             num_dims_emis = 3
          endif

          if( vndims < num_dims_emis ) then
             cycle
          elseif( vndims > num_dims_emis ) then
             ierr = pio_inq_varname (ncid, vid, varname)
             write(iulog,*) 'srf_emis_inti: Skipping variable ', trim(varname),', ndims = ',vndims, &
                  ' , species=',trim(emissions(m)%species)
             cycle
          end if

          if (unstructured) then
             allocate( dimids(vndims) )
             ierr = pio_inq_vardimid( ncid, vid, dimids )
             if ( any(dimids(:)==ncol_dimid) .and. any(dimids(:)==time_dimid) ) then
                emissions(m)%nsectors = emissions(m)%nsectors+1
                is_sector(vid)=.true.
             endif
             deallocate(dimids)
          else
             emissions(m)%nsectors = emissions(m)%nsectors+1
             is_sector(vid)=.true.
          end if

       enddo

       allocate( emissions(m)%sectors(emissions(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'srf_emis_inti: failed to allocate emissions(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1

       do vid = 1,nvars
          if( is_sector(vid) ) then
             ierr = pio_inq_varname(ncid, vid, emissions(m)%sectors(isec))
             isec = isec+1
          endif
       enddo
       deallocate(is_sector)

       ! Global attribute 'input_method' overrides the srf_emis_type namelist setting on
       ! a file-by-file basis.  If the emis file does not contain the 'input_method'
       ! attribute then the srf_emis_type namelist setting is used.
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ierr = pio_get_att(ncid, PIO_GLOBAL, 'input_method', file_interp_type)
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
       if ( ierr == PIO_NOERR) then
          l = GLC(file_interp_type)
          emis_type(1:l) = file_interp_type(1:l)
          emis_type(l+1:) = ' '
       else
          emis_type = trim(emissions_type)
       endif

       call pio_closefile (ncid)

       allocate(emissions(m)%file%in_pbuf(size(emissions(m)%sectors)))
       emissions(m)%file%in_pbuf(:) = .false.

       call trcdata_init( emissions(m)%sectors, &
                          emissions(m)%filename, filelist, datapath, &
                          emissions(m)%fields,  &
                          emissions(m)%file, &
                          rmv_file, emissions_cycle_yr, &
                          emissions_fixed_ymd, emissions_fixed_tod, trim(emis_type) )

       emissions(m)%units = emissions(m)%fields(1)%units

       call pbuf_set_field(pbuf2d, emissions(m)%bufndx, xnan)

       set_units: do n = 1,n_diags
          if (trim(emissions(m)%species)//'_srfemis'==names(n)) then
             units(n) = emissions(m)%fields(1)%units
             exit set_units
          end if
       end do set_units

    enddo files_loop

    do n = 1, n_diags
       call addfld(names(n), horiz_only, 'A', units(n), 'pbuf surf emis '//trim(names(n)))
    end do

  end subroutine surface_emissions_init

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine surface_emissions_adv( pbuf2d, state )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use tracer_data,  only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    !--------------------------------------------------------
    !	... Dummy arguments
    !--------------------------------------------------------
    type(physics_state), intent(in):: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    do m = 1,n_emis_files
       call advance_trcdata( emissions(m)%fields, emissions(m)%file, state, pbuf2d  )
       call pbuf_set_field(pbuf2d, emissions(m)%bufndx, 0._r8)
    end do

  end subroutine surface_emissions_adv

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine surface_emissions_set( lchnk, ncol, pbuf )
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    !--------------------------------------------------------
    !	... Dummy arguments
    !--------------------------------------------------------
    integer, intent(in) :: ncol
    integer, intent(in) :: lchnk
    type(physics_buffer_desc), pointer :: pbuf(:)

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  :: isec, m, n
    real(r8), pointer :: flux(:)

    !--------------------------------------------------------
    !	... set non-zero emissions
    !--------------------------------------------------------
    do m = 1,n_emis_files
       call pbuf_get_field(pbuf, emissions(m)%bufndx, flux)
       do isec = 1,emissions(m)%nsectors
          flux(:ncol) = flux(:ncol) + emissions(m)%scalefactor*emissions(m)%fields(isec)%data(:ncol,1,lchnk)
       enddo
    end do

    do n = 1, n_diags
       call pbuf_get_field(pbuf, indexes(n), flux)
       call outfld(names(n), flux(:ncol), ncol, lchnk)
    end do

  end subroutine surface_emissions_set

end module surface_emissions_mod
