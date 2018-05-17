module gcr_ionization

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use cam_abortutils, only : endrun
  use spmd_utils,     only : masterproc
  use tracer_data,    only : trfld,trfile
  use cam_logfile,    only : iulog
  use physics_buffer, only : physics_buffer_desc
  use physics_types,  only : physics_state
  use ppgrid,         only : begchunk, endchunk
  use ppgrid,         only : pcols, pver
  use tracer_data,    only : trcdata_init, advance_trcdata

  implicit none
  private 
  public :: gcr_ionization_readnl
  public :: gcr_ionization_init
  public :: gcr_ionization_adv
  public :: gcr_ionization_ionpairs

  type(trfld), pointer :: fields(:)
  type(trfile), save :: file

  character(len=32)  :: specifier(1) = 'prod'
  character(len=256) :: filename = 'NONE'
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  logical :: has_gcr_ionization = .false.

contains
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'gcr_ionization_readnl'

    character(len=16)  ::  gcr_ionization_fldname
    character(len=256) ::  gcr_ionization_filename
    character(len=256) ::  gcr_ionization_datapath
    character(len=256) ::  gcr_ionization_filelist
    character(len=32)  ::  gcr_ionization_datatype
    integer            ::  gcr_ionization_cycle_yr
    integer            ::  gcr_ionization_fixed_ymd
    integer            ::  gcr_ionization_fixed_tod

    namelist /gcr_ionization_nl/ &
         gcr_ionization_fldname, &
         gcr_ionization_filename, &
         gcr_ionization_datapath, &
         gcr_ionization_filelist, &
         gcr_ionization_datatype, &
         gcr_ionization_cycle_yr, &
         gcr_ionization_fixed_ymd, &
         gcr_ionization_fixed_tod

    gcr_ionization_fldname = specifier(1)
    gcr_ionization_filename = filename
    gcr_ionization_datapath = datapath
    gcr_ionization_filelist = filelist
    gcr_ionization_datatype = datatype
    gcr_ionization_cycle_yr = cycle_yr
    gcr_ionization_fixed_ymd = fixed_ymd
    gcr_ionization_fixed_tod = fixed_tod

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'gcr_ionization_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, gcr_ionization_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(gcr_ionization_fldname,  len(gcr_ionization_fldname),  mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_filename, len(gcr_ionization_filename), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_filelist, len(gcr_ionization_filelist), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_datapath, len(gcr_ionization_datapath), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_datatype, len(gcr_ionization_datatype), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_cycle_yr, 1, mpiint,  0, mpicom)
    call mpibcast(gcr_ionization_fixed_ymd,1, mpiint,  0, mpicom)
    call mpibcast(gcr_ionization_fixed_tod,1, mpiint,  0, mpicom)
#endif

    ! Update module variables with user settings.
    specifier(1) = gcr_ionization_fldname
    filename  = gcr_ionization_filename
    filelist  = gcr_ionization_filelist
    datapath  = gcr_ionization_datapath
    datatype  = gcr_ionization_datatype
    cycle_yr  = gcr_ionization_cycle_yr
    fixed_ymd = gcr_ionization_fixed_ymd
    fixed_tod = gcr_ionization_fixed_tod

    ! Turn on galactic cosmic rays if user has specified an input dataset.
    if (len_trim(filename) > 0 .and. filename.ne.'NONE') has_gcr_ionization = .true.

  end subroutine gcr_ionization_readnl

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_init()

    if (.not.has_gcr_ionization) return
    
    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype )

  end subroutine gcr_ionization_init

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_adv( pbuf2d, state )
    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (.not.has_gcr_ionization) return

    call advance_trcdata( fields, file, state, pbuf2d )

  end subroutine gcr_ionization_adv

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_ionpairs( ncol, lchnk, ionpairs )

    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    real(r8), intent(out) :: ionpairs(:,:)

    ionpairs(:,:) = 0._r8

    if (.not.has_gcr_ionization) return

    ionpairs(:ncol,:) = fields(1)%data(:ncol,:,lchnk)

  end subroutine gcr_ionization_ionpairs


end module gcr_ionization
