!-------------------------------------------------------------------------------
! Manages reading Upper Boundary Conditions (UBCs) from file
!-------------------------------------------------------------------------------
module upper_bc_file

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_kind_mod,  only: cx => shr_kind_cx ! 512
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc
  use cam_abortutils,only: endrun
  use cam_history,   only: addfld, horiz_only, outfld, fieldname_len

  use tracer_data,   only: trfld,trfile,MAXTRCRS

  implicit none
  private

  public :: upper_bc_file_readnl ! read namelist options
  public :: upper_bc_file_init   ! initialize
  public :: upper_bc_file_adv    ! advance data reader
  public :: upper_bc_file_get    ! returns UBC values
  public :: upper_bc_file_specified ! TRUE if UBC file is specified

  logical, protected :: upper_bc_file_specified = .false.

  ! private data members
  character(len=cx) :: ubc_file_path = 'NONE'
  character(len=32) :: ubc_file_input_type = 'NONE'
  integer :: ubc_file_cycle_yr  = -huge(1)
  integer :: ubc_file_fixed_ymd = -huge(1)
  integer :: ubc_file_fixed_tod = -huge(1)

  type(trfld), pointer :: fields(:) => null()
  type(trfile)         :: file

  integer :: num_ubc_flds = 0
  real(r8), allocatable :: ubc_fact(:)
  character(len=fieldname_len), allocatable :: hist_names(:)

contains

  !---------------------------------------------------------------------------
  ! read namelist options
  !---------------------------------------------------------------------------
  subroutine upper_bc_file_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils, only : mpicom, masterprocid, mpi_character, mpi_integer

    character(len=*), intent(in) :: nlfile

    integer :: unitn, ierr
    character(len=*), parameter :: prefix = 'upper_bc_file_readnl: '

    namelist /upper_bc_file_opts/ ubc_file_path, ubc_file_input_type
    namelist /upper_bc_file_opts/ ubc_file_cycle_yr, ubc_file_fixed_ymd, ubc_file_fixed_tod

    if (masterproc) then
       ! read namelist
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'upper_bc_file_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, upper_bc_file_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(prefix//'upper_bc_file_opts: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(ubc_file_path, len(ubc_file_path), mpi_character, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_file_path')
    call mpi_bcast(ubc_file_input_type, len(ubc_file_input_type), mpi_character, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_file_input_type')
    call mpi_bcast(ubc_file_fixed_ymd, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_file_fixed_ymd')
    call mpi_bcast(ubc_file_fixed_tod, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_file_fixed_tod')
    call mpi_bcast(ubc_file_cycle_yr,  1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_file_cycle_yr')

    upper_bc_file_specified = ubc_file_path /= 'NONE'

    if (masterproc) then
       write(iulog,*) prefix,'upper_bc_file_specified: ',upper_bc_file_specified
       write(iulog,*) prefix,'ubc_file_path = '//trim(ubc_file_path)
       write(iulog,*) prefix,'ubc_file_input_type = '//trim(ubc_file_input_type)
       write(iulog,*) prefix,'ubc_file_cycle_yr = ',ubc_file_cycle_yr
       write(iulog,*) prefix,'ubc_file_fixed_ymd = ',ubc_file_fixed_ymd
       write(iulog,*) prefix,'ubc_file_fixed_tod = ',ubc_file_fixed_tod
    end if

  end subroutine upper_bc_file_readnl

  !---------------------------------------------------------------------------
  ! initialize
  !---------------------------------------------------------------------------
  subroutine upper_bc_file_init( flds_list )
    use tracer_data, only: trcdata_init
    use constituents,only: cnst_get_ind, cnst_mw
    use physconst,   only: mwdry
    use string_utils,only: to_lower
    use ref_pres,    only: do_molec_diff

    character(len=*), intent(in) :: flds_list(:) ! flds specifier list

    integer :: m, ndx, ierr
    character(len=*), parameter :: prefix = 'upper_bc_file_init: '

    num_ubc_flds = size(flds_list)
    upper_bc_file_specified = upper_bc_file_specified .and. (num_ubc_flds>0)

    if (.not.upper_bc_file_specified) return

    allocate( ubc_fact(num_ubc_flds), stat=ierr )
    if (ierr /= 0) call endrun(prefix//'allocate error : ubc_fact')
    ubc_fact(:) = -huge(1._r8)

    allocate(file%in_pbuf(num_ubc_flds), stat=ierr)
    if (ierr /= 0) call endrun(prefix//'allocate error : file%in_pbuf')
    file%in_pbuf(:) = .false.

    call trcdata_init( flds_list, ubc_file_path, ' ', ' ', fields, file, .false., &
         ubc_file_cycle_yr, ubc_file_fixed_ymd, ubc_file_fixed_tod, ubc_file_input_type)

    if (do_molec_diff) then
       file%top_bndry = .true.
    else
       file%top_layer = .true.
    endif

    allocate(hist_names(num_ubc_flds), stat=ierr)
    if (ierr /= 0) call endrun(prefix//'allocate error : hist_names')
    hist_names = ' '

    do m = 1,num_ubc_flds

       call cnst_get_ind(trim(fields(m)%fldnam), ndx, abort=.true.)

       select case ( to_lower(trim(fields(m)%units)) )
       case ('k','kg/kg','kg kg-1','mmr')
          ubc_fact(m) = 1._r8
       case ('mol/mol','mole/mole','mol mol-1','vmr')
          ubc_fact(m) = cnst_mw(ndx)/mwdry
       case default
          call endrun('upper_bc_file_get: units are not recognized')
       end select

       hist_names(m) = trim(fields(m)%fldnam)//'_fubc'
       if ( to_lower(trim(fields(m)%units)) == 'k' ) then
          call addfld(hist_names(m), horiz_only, 'I', 'K', trim(fields(m)%fldnam)//' at upper boundary' )
       else
          call addfld(hist_names(m), horiz_only, 'I', 'kg/kg', trim(fields(m)%fldnam)//' at upper boundary' )
       end if

    end do

  end subroutine upper_bc_file_init

  !---------------------------------------------------------------------------
  ! advance data reader
  !---------------------------------------------------------------------------
  subroutine upper_bc_file_adv(pbuf2d, state)
    use tracer_data,    only : advance_trcdata
    use physics_types,  only : physics_state
    use physics_buffer, only : physics_buffer_desc

    ! args
    type(physics_state),    intent(in) :: state(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (.not.upper_bc_file_specified) return

    call advance_trcdata( fields, file, state, pbuf2d )

  end subroutine upper_bc_file_adv

  !---------------------------------------------------------------------------
  ! returns UBC values
  !---------------------------------------------------------------------------
  subroutine upper_bc_file_get(lchnk, ncol, val)

    integer, intent(in) :: ncol, lchnk
    real(r8), intent(out) :: val(:,:)

    integer :: m

    if (.not.upper_bc_file_specified) return

    do m = 1,num_ubc_flds
       val(:ncol,m) = ubc_fact(m)*fields(m)%data(:ncol,1,lchnk)
       call outfld( trim(hist_names(m)), val(:ncol,m), ncol, lchnk )
    enddo


  end subroutine upper_bc_file_get

end module upper_bc_file
