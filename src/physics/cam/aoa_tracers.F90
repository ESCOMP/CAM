!===============================================================================
! Age of air test tracers
! provides dissipation rate and surface fluxes for diagnostic constituents
!===============================================================================

module aoa_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver, begchunk, endchunk
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog
  use ref_pres,     only: pref_mid_norm
  use time_manager, only: get_curr_date, get_start_date
  use time_manager, only: is_leapyear, timemgr_get_calendar_cf, get_calday
  use pio,          only: file_desc_t, var_desc_t, pio_double, pio_noerr
  use pio,          only: pio_put_var, pio_def_var, pio_inq_varid, pio_get_var
  use cam_abortutils, only: endrun

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

  public :: aoa_tracers_define_restart
  public :: aoa_tracers_write_restart
  public :: aoa_tracers_read_restart

  ! Private module data

  integer, parameter :: ncnst=3  ! number of constituents implemented by this module

  ! constituent names
  character(len=4), parameter :: c_names(ncnst) = (/'AOA1', 'HORZ', 'VERT'/)

  ! constituent source/sink names
  character(len=7), parameter :: src_names(ncnst) = (/'AOA1SRC', 'HORZSRC', 'VERTSRC'/)

  integer :: ifirst = -1 ! global index of first constituent
  integer :: ixaoa  = -1 ! global index for AOA1SRC tracer
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

  real(r8), parameter :: NOTSET = -huge(1._r8)
  real(r8) :: mmr0 = NOTSET   ! initial lower boundary mmr

  type(var_desc_t) :: mmr0_desc

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine aoa_tracers_readnl(nlfile)

    use namelist_utils, only: find_group_name
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

    if (masterproc) then
       write(iulog,*) subname//' : aoa_tracers_flag: ',aoa_tracers_flag
       write(iulog,*) subname//' : aoa_read_from_ic_file: ',aoa_read_from_ic_file
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
                  longname='mixing ratio LB tracer', cam_outfld=.false., mixtype='dry')
    call cnst_add(c_names(2), mwdry, cpair, 1._r8, ixht,  readiv=aoa_read_from_ic_file, &
                  longname='horizontal tracer', cam_outfld=.false., mixtype='dry')
    call cnst_add(c_names(3), mwdry, cpair, 0._r8, ixvt,  readiv=aoa_read_from_ic_file, &
                  longname='vertical tracer', cam_outfld=.false., mixtype='dry')

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
    character(len=cl) :: errstr
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    ! if aoa_read_from_ic_file is FALSE this routine is called
    ! if aoa_read_from_ic_file is TRUE and the AOA tracer is not found in the IC file this routine is called
    ! abort the run if the user expects the IC to include AOA tracers and not found in the IC file
    if (aoa_read_from_ic_file) then
       write(errstr,*) 'AGE_OF_AIR_CONSTITUENTS ERROR: '//trim(name)//' not found in IC file'
       if (masterproc) then
          write (iulog,*) trim(errstr)
          write (iulog,*) ' --> Need to provide IC which includes AOA tracers or set aoa_read_from_ic_file to .FALSE.'
       end if
       call endrun(trim(errstr))
    end if

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

    use physics_types,  only: physics_state
    use gmean_mod,      only: gmean
    use time_manager,   only: is_first_step, is_first_restart_step

    type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state

    integer c, i, k, ncol, lchnk
    integer yr, mon, day, tod,  ymd
    real(r8) :: calday, dpy
    real(r8) :: mmr_arr(pcols,begchunk:endchunk)
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

    ! if AOA1 is not found in the IC file mmr0 is set to small value
    ! if AOA1 is initialized from IC file mmr0 is not set
    ! --> set mmr0 to global mean in lowest layer
    if (mmr0==NOTSET) then
       do lchnk = begchunk, endchunk
          ncol = phys_state(lchnk)%ncol
          ! AOA1 mmr in lowest model layer
          mmr_arr(:ncol,lchnk) = phys_state(lchnk)%q(:ncol,pver,ixaoa)
       end do
       call gmean(mmr_arr,mmr0) ! global mean
    end if

    if (is_first_step().or.is_first_restart_step()) then
       if (masterproc) then
         write(iulog,'(a,e20.12)') 'AGE_OF_AIR_CONSTITUENTS: mmr0 set to: ',mmr0
       end if
    end if

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
    real(r8), parameter :: per_yr = 0.02_r8   ! fractional increase per year

    real(r8) :: mmr_out(pcols,pver,ncnst)

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

    ! AOA1
    xmmr = mmr0*(1._r8 + per_yr*years)

    ! Lower boundary
    ptend%q(1:ncol,pver,ixaoa) = (xmmr - state%q(1:ncol,pver,ixaoa)) / dt

    ! Set upper boundary AOA1 tendency so that upper boundary AOA1 state will end up 0.0 mmr
    ptend%q(1:ncol,1,ixaoa) = -state%q(1:ncol,1,ixaoa) / dt

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

    ! output mixing ratios to history
    mmr_out(:ncol,:,1) = state%q(:ncol,:,ixaoa) + dt*ptend%q(1:ncol,:,ixaoa)
    mmr_out(:ncol,:,2) = state%q(:ncol,:,ixht)  + dt*ptend%q(1:ncol,:,ixht)
    mmr_out(:ncol,:,3) = state%q(:ncol,:,ixvt)  + dt*ptend%q(1:ncol,:,ixvt)

    call outfld (c_names(1), mmr_out(:,:,1),  pcols, lchnk)
    call outfld (c_names(2), mmr_out(:,:,2),  pcols, lchnk)
    call outfld (c_names(3), mmr_out(:,:,3),  pcols, lchnk)

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
       write(iulog,*) 'AGE_OF_AIR CONSTITUENTS: INITIALIZING ',cnst_name(m),m
    end if

    if (m == ixaoa) then

       ! AOA1
       mmr0 = 1.e-6_r8 ! initial lower boundary mmr (small non-zero value)
       q(:,:) = 0.0_r8 ! initialize AOA1 mixing ratios to zero

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
 subroutine aoa_tracers_define_restart(file)

   ! define variables to be written to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr

   ierr = pio_def_var(File, 'aoa_mmr0', pio_double, mmr0_desc)
   if (ierr/=pio_noerr) then
      call endrun('aoa_tracers_define_restart: ERROR pio_def_var aoa_mmr0')
   end if

 end subroutine aoa_tracers_define_restart

!=====================================================================
 subroutine aoa_tracers_write_restart(file)

   ! write variables to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr

   ierr = pio_put_var(File, mmr0_desc, mmr0)
   if (ierr/=pio_noerr) then
      call endrun('aoa_tracers_write_restart: ERROR pio_put_var aoa_mmr0')
   end if

 end subroutine aoa_tracers_write_restart

 !===============================================================================
 subroutine aoa_tracers_read_restart(file)

   ! read variables from restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   type(var_desc_t) :: vardesc

   ierr = pio_inq_varid(File, 'aoa_mmr0', vardesc)
   if (ierr/=pio_noerr) then
      call endrun('aoa_tracers_read_restart: ERROR pio_inq_varid aoa_mmr0')
   end if

   ierr = pio_get_var(File, vardesc, mmr0)
   if (ierr/=pio_noerr) then
      call endrun('aoa_tracers_read_restart: ERROR pio_get_var aoa_mmr0')
   end if
 end subroutine aoa_tracers_read_restart

end module aoa_tracers
