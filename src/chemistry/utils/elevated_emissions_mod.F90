module elevated_emissions_mod
  !---------------------------------------------------------------
  ! 	... elevalted emissions module
  !---------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8, shr_kind_cl
  use spmd_utils,    only : masterproc
  use cam_abortutils,only : endrun
  use ioFileMod,     only : getfil
  use cam_logfile,   only : iulog
  use tracer_data,   only : trfld,trfile
  use infnan,        only : nan, assignment(=)
  use cam_history,   only : addfld, outfld, add_default, fieldname_len

  implicit none

  type :: emission
     integer           :: bufndx
     real(r8)          :: scalefactor
     character(len=256):: filename
     character(len=16) :: species
     character(len=32) :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type emission

  private

  public :: elevated_emissions_readnl
  public :: elevated_emissions_reg
  public :: elevated_emissions_init
  public :: elevated_emissions_adv
  public :: elevated_emissions_set

  integer, parameter :: NMAX=50

  type(emission), allocatable :: elev_emis(:)
  integer                     :: n_emis_files = 0
  integer                     :: n_pbuf_flds = 0

  character(len=shr_kind_cl) :: elev_emis_specifier(NMAX) = ' '
  character(len=24) :: elev_emis_type
  integer :: elev_emis_cycle_yr
  integer :: elev_emis_fixed_ymd
  integer :: elev_emis_fixed_tod

  character(len=fieldname_len) :: names(NMAX) = ' '
  character(len=32) :: units(NMAX) = ' '
  integer :: indexes(NMAX) = -1
  integer :: n_diags = 0

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine elevated_emissions_readnl(nlfile)

    use namelist_utils, only : find_group_name
    use spmd_utils,     only : mpicom, masterprocid, mpi_integer, mpi_character

    character(len=*), intent(in) :: nlfile ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr, i
    logical :: logmsg
    character(len=*), parameter :: subname = 'elevated_emissions_readnl'

    namelist /elevated_emissions_opts/ elev_emis_specifier, elev_emis_type, elev_emis_cycle_yr, &
                                       elev_emis_fixed_ymd, elev_emis_fixed_tod

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'elevated_emissions_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, elevated_emissions_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)

    end if

    ! Broadcast namelist variables
    call mpi_bcast(elev_emis_specifier,len(elev_emis_specifier(1))*NMAX, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(elev_emis_type, len(elev_emis_type), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(elev_emis_cycle_yr, 1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(elev_emis_fixed_ymd, 1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(elev_emis_fixed_tod, 1, mpi_integer, masterprocid, mpicom, ierr)

    logmsg = .false.
    if (masterproc) then
       do i = 1,NMAX
          if (len_trim(elev_emis_specifier(i))>0) then
             logmsg = .true.
             write(iulog,'(2a)') subname,': elev_emis_specifier: ',trim(elev_emis_specifier(i))
          endif
       enddo
       if (logmsg) then
          write(iulog,*) subname,': elev_emis_type: ',elev_emis_type
          write(iulog,*) subname,': elev_emis_cycle_yr: ',elev_emis_cycle_yr
          write(iulog,*) subname,': elev_emis_fixed_ymd: ',elev_emis_fixed_ymd
          write(iulog,*) subname,': elev_emis_fixed_tod: ',elev_emis_fixed_tod
       endif
    endif

  end subroutine elevated_emissions_readnl

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine elevated_emissions_reg( )
    use m_MergeSorts,   only : IndexSort
    use physics_buffer, only : pbuf_add_field, dtype_r8, pbuf_get_index
    use ppgrid,         only : pcols, pver

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

    character(len=*), parameter :: prefix = 'elevated_emissions_reg: '
    kk = 0
    nn = 0
    indx(:) = 0
    emis_species = ' '
    emis_indexes = -1
    emis_filenam = 'NONE'

    count_emis: do n=1,size(elev_emis_specifier)
       if ( len_trim(elev_emis_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(elev_emis_specifier(n),'->')
       spc_name = trim(adjustl(elev_emis_specifier(n)(:i-1)))

       ! need to parse out scalefactor ...
       tmp_string = adjustl(elev_emis_specifier(n)(i+2:))
       j = scan( tmp_string, '*' )
       if (j>0) then
          xchr = tmp_string(1:j-1) ! get the multipler (left of the '*')
          read( xchr, * ) xdbl   ! convert the string to a real
          tmp_string = adjustl(tmp_string(j+1:)) ! get the filepath name (right of the '*')
       else
          xdbl = 1._r8
       endif
       filename = trim(tmp_string)

       bname = trim(spc_name)//'_elevemis'

       m = pbuf_get_index(bname,errcode=err)
       if (m<1) then
          call pbuf_add_field(bname, 'physpkg', dtype_r8, (/pcols,pver/), m)
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

    if (masterproc) write(iulog,*) prefix,' n_emis_files = ',n_emis_files

    allocate( elev_emis(n_emis_files), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'elev_emis_inti: failed to allocate emissions array; error = ',astat
       call endrun('elev_emis_inti: failed to allocate emissions array')
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
       elev_emis(m)%bufndx           = emis_indexes(indx(m))
       elev_emis(m)%species          = emis_species(indx(m))
       elev_emis(m)%filename         = emis_filenam(indx(m))
       elev_emis(m)%scalefactor      = emis_scalefactor(indx(m))
    enddo
  end subroutine elevated_emissions_reg

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine elevated_emissions_init(pbuf2d)
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

       elev_emis(m)%nsectors = 0
       call getfil (elev_emis(m)%filename, locfn, 0)
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
             num_dims_emis = 3
          else
             num_dims_emis = 4
          endif

          if( vndims < num_dims_emis ) then
             cycle
          elseif( vndims > num_dims_emis ) then
             ierr = pio_inq_varname (ncid, vid, varname)
             write(iulog,*) 'elev_emis_inti: Skipping variable ', trim(varname),', ndims = ',vndims, &
                  ' , species=',trim(elev_emis(m)%species)
             cycle
          end if

          if (unstructured) then
             allocate( dimids(vndims) )
             ierr = pio_inq_vardimid( ncid, vid, dimids )
             if ( any(dimids(:)==ncol_dimid) .and. any(dimids(:)==time_dimid) ) then
                elev_emis(m)%nsectors = elev_emis(m)%nsectors+1
                is_sector(vid)=.true.
             endif
             deallocate(dimids)
          else
             elev_emis(m)%nsectors = elev_emis(m)%nsectors+1
             is_sector(vid)=.true.
          end if

       enddo

       allocate( elev_emis(m)%sectors(elev_emis(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'elev_emis_inti: failed to allocate elev_emis(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1

       do vid = 1,nvars
          if( is_sector(vid) ) then
             ierr = pio_inq_varname(ncid, vid, elev_emis(m)%sectors(isec))
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
          emis_type = trim(elev_emis_type)
       endif

       call pio_closefile (ncid)

       allocate(elev_emis(m)%file%in_pbuf(size(elev_emis(m)%sectors)))
       elev_emis(m)%file%in_pbuf(:) = .false.

       call trcdata_init( elev_emis(m)%sectors, &
                          elev_emis(m)%filename, filelist, datapath, &
                          elev_emis(m)%fields,  &
                          elev_emis(m)%file, &
                          rmv_file, elev_emis_cycle_yr, &
                          elev_emis_fixed_ymd, elev_emis_fixed_tod, trim(emis_type) )

       elev_emis(m)%units = elev_emis(m)%fields(1)%units

       call pbuf_set_field(pbuf2d, elev_emis(m)%bufndx, xnan)

       set_units: do n = 1,n_diags
          if (trim(elev_emis(m)%species)//'_elevemis'==names(n)) then
             units(n) = elev_emis(m)%fields(1)%units
             exit set_units
          end if
       end do set_units

    enddo files_loop

    do n = 1, n_diags
       call addfld(names(n), (/ 'lev' /), 'A', units(n), 'pbuf elev emis '//trim(names(n)))
    end do

  end subroutine elevated_emissions_init

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine elevated_emissions_adv( pbuf2d, state )
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
       call advance_trcdata( elev_emis(m)%fields, elev_emis(m)%file, state, pbuf2d  )
       call pbuf_set_field(pbuf2d, elev_emis(m)%bufndx, 0._r8)
    end do

  end subroutine elevated_emissions_adv

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine elevated_emissions_set( lchnk, ncol, pbuf )
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
    real(r8), pointer :: flux(:,:)

    !--------------------------------------------------------
    !	... set non-zero emissions
    !--------------------------------------------------------
    do m = 1,n_emis_files
       call pbuf_get_field(pbuf, elev_emis(m)%bufndx, flux)
       do isec = 1,elev_emis(m)%nsectors
          flux(:ncol,:) = flux(:ncol,:) + elev_emis(m)%scalefactor*elev_emis(m)%fields(isec)%data(:ncol,:,lchnk)
       enddo
    end do

    do n = 1, n_diags
       call pbuf_get_field(pbuf, indexes(n), flux)
       call outfld(names(n), flux(:ncol,:), ncol, lchnk)
    end do

  end subroutine elevated_emissions_set

end module elevated_emissions_mod
