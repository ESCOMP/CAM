!--------------------------------------------------------------------------------
! Provids electron fluxes read from input data set
!--------------------------------------------------------------------------------
module mee_fluxes
  use shr_kind_mod,   only : r8 => shr_kind_r8, cl=> shr_kind_cl
  use spmd_utils,     only : masterproc
  use cam_logfile,    only : iulog
  use cam_abortutils, only : endrun
  use input_data_utils, only : time_coordinate
  use pio, only : file_desc_t, var_desc_t, pio_get_var, pio_inq_varid
  use pio, only : PIO_NOWRITE, pio_inq_dimid, pio_inq_dimlen
  use infnan, only: isnan

  implicit none

  private
  public :: mee_fluxes_readnl
  public :: mee_fluxes_init
  public :: mee_fluxes_final
  public :: mee_fluxes_adv     ! read and time interpolate fluxes
  public :: mee_fluxes_extract ! interpolate flux to column L-shell
  public :: mee_fluxes_active  ! true when input flux file is specified
  public :: mee_fluxes_denergy ! energy bin widths
  public :: mee_fluxes_energy  ! center of each energy bin
  public :: mee_fluxes_nenergy ! number of energy bins

  real(r8),protected, pointer :: mee_fluxes_denergy(:) => null()
  real(r8),protected, pointer :: mee_fluxes_energy(:) => null()
  integer, protected :: mee_fluxes_nenergy
  logical, protected :: mee_fluxes_active = .false.

  real(r8), allocatable :: lshell(:)
  real(r8), allocatable :: indata(:,:,:)
  real(r8), allocatable :: influx(:,:)
  logical , allocatable :: valflx(:,:)

  character(len=cl) :: mee_fluxes_filepath = 'NONE'
  logical :: mee_fluxes_fillin = .false.

  type(time_coordinate) :: time_coord
  integer :: nlshells

  type(file_desc_t) :: file_id
  type(var_desc_t) :: flux_var_id

contains

  !-----------------------------------------------------------------------------
  ! read namelist options
  !-----------------------------------------------------------------------------
  subroutine mee_fluxes_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use spmd_utils,     only: mpicom, mpi_character, mpi_logical, masterprocid

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'mee_fluxes_readnl'

    namelist /mee_fluxes_opts/ mee_fluxes_filepath, mee_fluxes_fillin

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'mee_fluxes_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, mee_fluxes_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(mee_fluxes_filepath, len(mee_fluxes_filepath), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(mee_fluxes_fillin, 1, mpi_logical, masterprocid, mpicom, ierr)

    mee_fluxes_active = mee_fluxes_filepath /= 'NONE'

    if ( masterproc ) then
       if ( mee_fluxes_active ) then
          write(iulog,*) subname//':: Input electron fluxes filepath: '//trim(mee_fluxes_filepath)
          write(iulog,*) subname//':: Fill in missing fluxes with vdk-derived fluxes: ', mee_fluxes_fillin
       else
          write(iulog,*) subname//':: Electron fluxes are not prescribed'
       end if
    end if

  end subroutine mee_fluxes_readnl

  !-----------------------------------------------------------------------------
  ! intialize -- allocate memory and read coordinate data
  !-----------------------------------------------------------------------------
  subroutine mee_fluxes_init()
    use cam_pio_utils,  only : cam_pio_openfile
    use ioFileMod,      only : getfil

    character(len=cl) :: filen
    integer :: ierr, dimid, varid
    real(r8), allocatable :: logdelta(:)
    character(len=*), parameter :: subname = 'mee_fluxes_init: '

    if (.not.mee_fluxes_active) return

    call time_coord%initialize( mee_fluxes_filepath, force_time_interp=.true. )

    call getfil( mee_fluxes_filepath, filen, 0 )
    call cam_pio_openfile( file_id, filen, PIO_NOWRITE )

    ierr = pio_inq_dimid(file_id, 'energy', dimid)
    ierr = pio_inq_dimlen(file_id, dimid, mee_fluxes_nenergy )

    ierr = pio_inq_dimid(file_id, 'lshell', dimid)
    ierr = pio_inq_dimlen(file_id, dimid, nlshells )

    ierr = pio_inq_varid(file_id, 'RBSP_flux_scaled', flux_var_id)

    allocate( indata( mee_fluxes_nenergy, nlshells, 2 ), stat=ierr)
    if (ierr/=0) call endrun(subname//'not able to allocate indata')
    allocate( influx( mee_fluxes_nenergy, nlshells ), stat=ierr )
    if (ierr/=0) call endrun(subname//'not able to allocate influx')
    allocate( valflx( mee_fluxes_nenergy, nlshells ), stat=ierr )
    if (ierr/=0) call endrun(subname//'not able to allocate valflx')
    allocate( mee_fluxes_energy( mee_fluxes_nenergy ), stat=ierr )
    if (ierr/=0) call endrun(subname//'not able to allocate mee_fluxes_energy')
    allocate( mee_fluxes_denergy( mee_fluxes_nenergy ), stat=ierr )
    if (ierr/=0) call endrun(subname//'not able to allocate mee_fluxes_denergy')
    allocate( logdelta( mee_fluxes_nenergy ), stat=ierr )
    if (ierr/=0) call endrun(subname//'not able to allocate logdelta')
    allocate( lshell( nlshells ), stat=ierr )
    if (ierr/=0) call endrun(subname//'not able to allocate lshell')


    ierr = pio_inq_varid(file_id, 'energy', varid)
    ierr = pio_get_var( file_id, varid, mee_fluxes_energy)

    ierr = pio_inq_varid(file_id, 'lshell', varid)
    ierr = pio_get_var( file_id, varid, lshell)

    logdelta(2:) =  log(mee_fluxes_energy(2:mee_fluxes_nenergy))-log(mee_fluxes_energy(1:mee_fluxes_nenergy-1))
    logdelta(1) = logdelta(2)
    mee_fluxes_denergy(:) = exp( log(mee_fluxes_energy(:)) + 0.5_r8*logdelta(:) ) &
                          - exp( log(mee_fluxes_energy(:)) - 0.5_r8*logdelta(:) )

    deallocate(logdelta)

    call read_fluxes()

  end subroutine mee_fluxes_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_fluxes_final()
    use cam_pio_utils, only : cam_pio_closefile

    if (.not.mee_fluxes_active) return

    call cam_pio_closefile(file_id)

    deallocate(indata)
    deallocate(influx)
    deallocate(valflx)
    deallocate(lshell)

    deallocate(mee_fluxes_energy)
    deallocate(mee_fluxes_denergy)

    nullify(mee_fluxes_energy)
    nullify(mee_fluxes_denergy)

  end subroutine mee_fluxes_final

  !-----------------------------------------------------------------------------
  ! time interpolate the input fluxes
  ! reads data as needed
  !-----------------------------------------------------------------------------
  subroutine mee_fluxes_adv

    if (.not.mee_fluxes_active) return

    call time_coord%advance()

    if ( time_coord%read_more() ) then
       call read_fluxes( )
    endif

    influx(:,:) = 0._r8

    valflx(:,:) = (.not.isnan(indata(:,:,1))) .and. (.not.isnan(indata(:,:,2)))
    where ( valflx(:,:) )
       influx(:,:) = time_coord%wghts(1)*indata(:,:,1) + time_coord%wghts(2)*indata(:,:,2)
    end where

    if (any(isnan(influx))) then
       call endrun('mee_fluxes_adv -- influx has NaNs')
    end if

  end subroutine mee_fluxes_adv

  !-----------------------------------------------------------------------------
  ! linear interpolate fluxes in L-shell where the fluxes are valid
  !-----------------------------------------------------------------------------
  subroutine mee_fluxes_extract( l_shell, fluxes, valid )

    real(r8), intent(in) :: l_shell
    real(r8), intent(out) :: fluxes(mee_fluxes_nenergy)
    logical, intent(out) :: valid(mee_fluxes_nenergy)

    integer :: i, ndx1, ndx2
    logical :: found
    real(r8) :: wght1,wght2

    valid(:) = .not. mee_fluxes_fillin
    fluxes(:) = 0._r8

    if (.not.mee_fluxes_active) return

    found = .false.

    findloop: do i = 1,nlshells-1
       if ( l_shell>=lshell(i) .and. l_shell<=lshell(i+1) ) then
          ndx1=i
          ndx2=i+1
          wght2 = (l_shell-lshell(ndx1))/(lshell(ndx2)-lshell(ndx1))
          wght1 = 1._r8 - wght2
          found = .true.
          exit findloop
       endif
    end do findloop

    if (found) then
       if (mee_fluxes_fillin) then
          valid(:) = valflx(:,ndx1) .and. valflx(:,ndx2)
       end if
       where( valid(:) )
          fluxes(:) = wght1*influx(:,ndx1) + wght2*influx(:,ndx2)
       end where
    end if

  end subroutine mee_fluxes_extract

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine read_fluxes()

    ! local vars
    integer :: ierr, cnt(4), start(4)

    cnt = (/1, mee_fluxes_nenergy, nlshells, 2/)

    ! use the 50 percentile level data ( index 3 )
    start = (/3, 1, 1, time_coord%indxs(1)/)

    ! float RBSP_flux_scaled(time, lshell, energy, percentiles) ;
    ierr = pio_get_var( file_id, flux_var_id, start, cnt, indata )

  end subroutine read_fluxes

end module mee_fluxes
