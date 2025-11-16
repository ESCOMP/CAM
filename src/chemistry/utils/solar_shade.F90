module solar_shade
  use shr_kind_mod, only: r8 => shr_kind_r8
  use solar_irrad_data, only: nbins, we
  use physics_types,only : physics_state
  use ppgrid, only : pcols, begchunk, endchunk
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc
  use cam_abortutils, only: endrun

  implicit none
  private
  public :: solar_shade_readnl
  public :: solar_shade_init

  real(r8), public, protected, allocatable :: sun_shade(:,:,:) ! chunk, col, wavelen dependent

  real(r8) :: soll0 = 1._r8  ! namelist variable
  real(r8) :: soll1 = 1._r8  ! namelist variable
  real(r8) :: soll2 = 1._r8  ! namelist variable

contains

  !============================================================================
  subroutine solar_shade_readnl(nlfile)
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils,      only: mpicom, masterprocid, mpi_real8, mpi_success

    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    character(len=*), parameter :: prefix = 'solar_shade_readnl: '

    namelist /solar_shade_opts/ soll0, soll1, soll2

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'solar_shade_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, solar_shade_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(prefix//'ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    call mpi_bcast(soll0, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(prefix//'ERROR bcast soll0')
    end if
    call mpi_bcast(soll1, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(prefix//'ERROR bcast soll1')
    end if
    call mpi_bcast(soll2, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(prefix//'ERROR bcast soll2')
    end if

    if (masterproc) then
       write(iulog,*) prefix,'factor soll0 = ', soll0
       write(iulog,*) prefix,'factor soll1 = ', soll1
       write(iulog,*) prefix,'factor soll2 = ', soll2
    end if

  end subroutine solar_shade_readnl

  !============================================================================
  subroutine solar_shade_init()
    use phys_grid,     only : get_ncols_p, get_rlat_all_p, get_rlon_all_p

    integer :: lchnk, ncol, icol
    real(r8) :: clat(pcols)     ! current latitudes(radians)
    real(r8) :: modval

    allocate(sun_shade(nbins, pcols, begchunk:endchunk))
    sun_shade = 1._r8

    ! copy of Ben's insolation function
    do lchnk=begchunk,endchunk
       ncol = get_ncols_p(lchnk)
       call get_rlat_all_p(lchnk, ncol, clat)
       do icol = 1,ncol
          modval =          (1._r8-soll0)
          modval = modval + (1._r8-soll1)*sin(clat(icol))
          modval = modval + (1._r8-soll2)*0.5_r8*(3._r8*(sin(clat(icol)))**2-1._r8)
          sun_shade(:nbins, icol, lchnk) = 1._r8-modval
       end do
    end do

  end subroutine solar_shade_init

end module solar_shade
