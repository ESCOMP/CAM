!===============================================================================
! Age of air test tracers
! provides dissipation rate and surface fluxes for diagnostic constituents
!===============================================================================

module aoa_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog
  use ref_pres,     only: pref_mid_norm
  use time_manager, only: get_curr_date, get_start_date
  use time_manager, only: is_leapyear, timemgr_get_calendar_cf, get_calday

  implicit none
  private

  ! Public interfaces
  public :: aoa_tracers_register         ! register constituents
  public :: aoa_tracers_implements_cnst  ! true if named constituent is implemented by this package
  public :: aoa_tracers_init_cnst        ! initialize constituent field
  public :: aoa_tracers_init             ! initialize history fields, datasets
  public :: aoa_tracers_timestep_init    ! place to perform per timestep initialization
  public :: aoa_tracers_timestep_tend    ! calculate tendencies
  public :: aoa_tracers_readnl           ! read namelist options

  ! Private module data

  integer, parameter :: ncnst=3  ! number of constituents implemented by this module

  ! constituent names
  character(len=6), parameter :: c_names(ncnst) = (/'AOAMF ', 'HORZ  ', 'VERT  '/)

  ! constituent source/sink names
  character(len=8), parameter :: src_names(ncnst) = (/'AOAMFSRC', 'HORZSRC ', 'VERTSRC '/)

  integer :: ifirst = -1 ! global index of first constituent
  integer :: ixaoa  = -1 ! global index for AOAMFSRC tracer
  integer :: ixht   = -1 ! global index for HORZ tracer
  integer :: ixvt   = -1 ! global index for VERT tracer

  ! Data from namelist variables
  logical :: aoa_tracers_flag  = .false.    ! true => turn on test tracer code, namelist variable
  logical :: aoa_read_from_ic_file = .true. ! true => tracers initialized from IC file

  real(r8),  parameter ::  treldays = 15._r8
  real(r8),  parameter ::  vert_offset = 10._r8

  ! 15-days used for diagnostic of transport circulation and K-tensors
  ! relaxation (in the original papers PM-1987 and YSGD-2000) => Zonal Mean
  ! to evaluate eddy-fluxes for 2D-diagnostics, here relaxation to the GLOBAL MEAN  IC
  ! it may help to keep gradients but will rule-out 2D-transport diagnostics
  ! in km  to avoid negative values of  vertical tracers
  ! VERT(k) = -7._r8*alog(hyam(k)+hybm(k)) + vert_offset

  ! PM-1987:
  ! Plumb, R. A., and J. D. Mahlman (1987), The zonally averaged transport
  ! characteristics of the GFDL general circulation/transport model,
  ! J. Atmos.Sci.,44, 298-327

  ! YSGD-2000:
  ! Yudin, Valery A., Sergey P. Smyshlyaev, Marvin A. Geller, Victor L. Dvortsov, 2000:
  ! Transport Diagnostics of GCMs and Implications for 2D Chemistry-Transport Model of
  ! Troposphere and Stratosphere. J. Atmos. Sci., 57, 673-699.
  ! doi: http://dx.doi.org/10.1175/1520-0469(2000)057<0673:TDOGAI>2.0.CO;2

  real(r8) :: qrel_vert(pver) = -huge(1._r8)  ! = -7._r8*log(pref_mid_norm(k)) + vert_offset

  integer :: yr0 = -huge(1)
  real(r8) :: calday0 = -huge(1._r8)
  real(r8) :: years = -huge(1._r8)

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine aoa_tracers_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use cam_abortutils, only: endrun
    use spmd_utils,     only: mpicom, masterprocid, mpi_logical, mpi_success

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aoa_tracers_readnl'

    namelist /aoa_tracers_nl/ aoa_tracers_flag, aoa_read_from_ic_file

    !-----------------------------------------------------------------------------

    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aoa_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aoa_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(aoa_tracers_flag, 1, mpi_logical, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: aoa_tracers_flag')
    end if
    call mpi_bcast(aoa_read_from_ic_file, 1, mpi_logical, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: aoa_read_from_ic_file')
    end if

  endsubroutine aoa_tracers_readnl

!================================================================================

  subroutine aoa_tracers_register
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents
    !
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    integer :: k

    if (.not. aoa_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixaoa, readiv=aoa_read_from_ic_file, &
                  longname='mixing ratio LB tracer')

    call cnst_add(c_names(2), mwdry, cpair, 1._r8, ixht,   readiv=aoa_read_from_ic_file, &
                  longname='horizontal tracer')
    call cnst_add(c_names(3), mwdry, cpair, 0._r8, ixvt,   readiv=aoa_read_from_ic_file, &
                  longname='vertical tracer')

    ifirst = ixaoa

    do k = 1,pver
       qrel_vert(k) = -7._r8*log(pref_mid_norm(k)) + vert_offset
    enddo

  end subroutine aoa_tracers_register

!===============================================================================

  function aoa_tracers_implements_cnst(name)
    !-----------------------------------------------------------------------
    !
    ! Purpose: return true if specified constituent is implemented by this package
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: aoa_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    aoa_tracers_implements_cnst = .false.

    if (.not. aoa_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          aoa_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function aoa_tracers_implements_cnst

!===============================================================================

  subroutine aoa_tracers_init_cnst(name, latvals, lonvals, mask, q)

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize test tracers mixing ratio fields
    !  This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, latvals, lonvals, mask, q)
       endif
    end do

  end subroutine aoa_tracers_init_cnst

!===============================================================================

  subroutine aoa_tracers_init

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default

    integer :: m, mm
    integer :: yr, mon, day, sec, ymd

    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    ! Set names of tendencies and declare them as history variables

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))
       call addfld(src_names(m),  (/ 'lev' /), 'A', 'kg/kg/s', trim(cnst_name(mm))//' source/sink')

       call add_default (cnst_name(mm), 1, ' ')
       call add_default (src_names(m),  1, ' ')
    end do

    call get_start_date(yr, mon, day, sec)

    ymd = yr*10000 + mon*100 + day

    yr0 = yr
    calday0 = get_calday(ymd, sec)

  end subroutine aoa_tracers_init

!===============================================================================

  subroutine aoa_tracers_timestep_init( phys_state )
    !-----------------------------------------------------------------------
    ! Provides a place to reinitialize diagnostic constituents HORZ and VERT
    !-----------------------------------------------------------------------

    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state

    type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state

    integer c, i, k, ncol
    integer yr, mon, day, tod,  ymd
    real(r8) :: calday, dpy
    !--------------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    call get_curr_date (yr,mon,day,tod)

    if ( day == 1 .and. tod == 0) then
       if (masterproc) then
         write(iulog,*) 'AGE_OF_AIR_CONSTITUENTS: RE-INITIALIZING HORZ/VERT CONSTITUENTS'
       endif

       do c = begchunk, endchunk
          ncol = phys_state(c)%ncol
          do k = 1, pver
             do i = 1, ncol
                phys_state(c)%q(i,k,ixht) = 2._r8 + sin(phys_state(c)%lat(i))
                phys_state(c)%q(i,k,ixvt) = qrel_vert(k)
             end do
          end do
       end do

    end if

    ymd = yr*10000 + mon*100 + day
    calday = get_calday(ymd, tod)

    dpy = 365._r8
    if (timemgr_get_calendar_cf() == 'gregorian' .and. is_leapyear(yr)) then
       dpy = 366._r8
    end if
    years = (yr-yr0) + (calday-calday0)/dpy

  end subroutine aoa_tracers_timestep_init

!===============================================================================

  subroutine aoa_tracers_timestep_tend(state, ptend, dt)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use cam_history,   only: outfld

    ! Arguments
    type(physics_state), intent(in)    :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(in)    :: dt                 ! timestep size (sec)

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk                          ! chunk identifier
    integer :: ncol                           ! no. of column in chunk
    real(r8) :: qrel                          ! value to be relaxed to
    real(r8) :: xhorz                         ! updated value of HORZ
    real(r8) :: xvert                         ! updated value of VERT
    logical  :: lq(pcnst)
    real(r8) :: teul                          ! relaxation in  1/sec*dt/2 = k*dt/2
    real(r8) :: wimp                          !     1./(1.+ k*dt/2)
    real(r8) :: wsrc                          !  teul*wimp

    real(r8) :: xmmr
    real(r8), parameter :: mmr0 = 1.0e-6_r8   ! initial lower boundary mmr
    real(r8), parameter :: per_yr = 0.02_r8   ! fractional increase per year

    !------------------------------------------------------------------

    teul = .5_r8*dt/(86400._r8 * treldays)   ! 1/2 for the semi-implicit scheme if dt=time step
    wimp = 1._r8/(1._r8 +teul)
    wsrc = teul*wimp

    if (.not. aoa_tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    lq(:)     = .FALSE.
    lq(ixaoa) = .TRUE.
    lq(ixht)  = .TRUE.
    lq(ixvt)  = .TRUE.

    call physics_ptend_init(ptend,state%psetcols, 'aoa_tracers', lq=lq)

    lchnk = state%lchnk
    ncol  = state%ncol

    ! AOAMF
    xmmr = mmr0*(1._r8 + per_yr*years)
    ptend%q(1:ncol,pver,ixaoa) = (xmmr - state%q(1:ncol,pver,ixaoa)) / dt

    do k = 1, pver
       do i = 1, ncol

          ! HORZ
          qrel              = 2._r8 + sin(state%lat(i))          ! qrel  should zonal mean
          xhorz             = state%q(i,k,ixht)*wimp + wsrc*qrel ! Xnew = weight*3D-tracer + (1.-weight)*1D-tracer
          ptend%q(i,k,ixht) = (xhorz - state%q(i,k,ixht)) / dt   ! Xnew = weight*3D-tracer + (1.-weight)*2D-tracer  zonal mean
                                                                 !  Can be still used .... to diagnose fluxes OT-tracers
          ! VERT
          qrel              = qrel_vert(k)                       ! qrel  should zonal mean
          xvert             = wimp*state%q(i,k,ixvt) + wsrc*qrel
          ptend%q(i,k,ixvt) = (xvert - state%q(i,k,ixvt)) / dt

       end do

    end do

    ! record tendencies on history files
    call outfld (src_names(1), ptend%q(:,:,ixaoa),  pcols, lchnk)
    call outfld (src_names(2), ptend%q(:,:,ixht),   pcols, lchnk)
    call outfld (src_names(3), ptend%q(:,:,ixvt),   pcols, lchnk)

  end subroutine aoa_tracers_timestep_tend

!===========================================================================

  subroutine init_cnst_3d(m, latvals, lonvals, mask, q)

    integer,  intent(in)  :: m          ! global constituent index
    real(r8), intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8), intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8), intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol,plev)

    integer :: j, k, gsize
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'AGE-OF-AIR CONSTITUENTS: INITIALIZING ',cnst_name(m),m
    end if

    if (m == ixaoa) then

       ! AOAMF
       q(:,:) = 0.0_r8

    else if (m == ixht) then

       ! HORZ
       gsize = size(q, 1)
       do j = 1, gsize
          q(j,:) = 2._r8 + sin(latvals(j))
       end do

    else if (m == ixvt) then

       ! VERT
       do k = 1, pver
          do j = 1, size(q,1)
             q(j,k) = qrel_vert(k)
          end do
       end do

    end if

  end subroutine init_cnst_3d

!=====================================================================

end module aoa_tracers
