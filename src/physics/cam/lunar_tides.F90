! Module to add lunar tide forcing for the M2
! gravitational lunar tide based on the empirical 
! formula of Champman and Lindzen (1970) and as
! implemented in Pedatella, Liu, and Richmond (2012, JGR)

module lunar_tides
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use phys_control,   only: use_simple_phys
  use spmd_utils,     only: masterproc
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog

  implicit none
  private

  public :: lunar_tides_readnl
  public :: lunar_tides_init
  public :: lunar_tides_tend
  
  ! lunar tides forcing option
  logical  :: apply_lunar_tides = .false.

contains

  !==========================================================================
  !==========================================================================
  subroutine lunar_tides_readnl(nlfile)
    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: mpicom, mstrid=>masterprocid, mpi_logical

    ! File containing namelist input.
    character(len=*), intent(in) :: nlfile

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: sub = 'lunar_tides_readnl'

    namelist /lunar_tides_opts/ apply_lunar_tides

    if (use_simple_phys) return

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'lunar_tides_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, lunar_tides_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(sub // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    call mpi_bcast(apply_lunar_tides, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: apply_lunar_tides")

    if (masterproc) then
       write(iulog,*) 'lunar_tides_readnl: apply_lunar_tides: ',apply_lunar_tides
    end if

  end subroutine lunar_tides_readnl

  !==========================================================================
  !==========================================================================
  subroutine lunar_tides_init
    use cam_history, only: addfld
    use time_manager,only: timemgr_get_calendar_cf

    if (apply_lunar_tides) then
       if (timemgr_get_calendar_cf().ne.'gregorian') then
          call endrun('lunar_tides_init: calendar must be gregorian')
       endif
       call addfld('UT_LUNAR', (/ 'lev' /), 'A','m/s2','Zonal wind tendency due to lunar tides')
       call addfld('VT_LUNAR', (/ 'lev' /), 'A','m/s2','Meridional wind tendency due to lunar tides')
    end if
 
  end subroutine lunar_tides_init
  
  !==========================================================================
  !==========================================================================
  subroutine lunar_tides_tend( state, ptend )
    use time_manager, only: get_curr_date, get_julday
    use physconst,    only: pi, rearth
    use ppgrid,       only: pver
    use cam_history,  only: outfld

    type(physics_state), intent(in) :: state
    type(physics_ptend), intent(out):: ptend

    integer  :: tod,yr,mm,dd
    real(r8) :: jd,nu,lt,lun_lt

    integer :: i, k

    real(r8), parameter :: deg2hrs = 1._r8/15._r8
    real(r8), parameter :: rad2deg = 180._r8/pi
    real(r8), parameter :: rad2hrs = rad2deg*deg2hrs
    real(r8), parameter :: tod2hrs = 24._r8/86400._r8
    real(r8), parameter :: hrs2rad = 1._r8/rad2hrs
    
    if (apply_lunar_tides) then

       call physics_ptend_init(ptend, state%psetcols, "Lunar Tides", lu=.true., lv=.true. )

       ! calculate the current date:
       call get_curr_date(yr,mm,dd,tod)
       ! convert date to Julian centuries
       jd = get_julday(yr,mm,dd,tod)
       ! calculation relies on time from noon on December 31, 1899, so
       ! subtract 2415020, which corresponds to the Julian date for Dec. 31 1899.
       jd = jd - 2415020._r8
       jd = jd / 36525._r8 ! convert to julian centuries

       ! Calculate the lunar local time (nu) based on the the time
       ! in Julian centuries using the formula given in Chapman and Lindzen (1970)
       nu = -9.26009_r8 + 445267.12165_r8*jd+0.00168_r8*jd*jd !nu in degrees

       do i=1,state%ncol
          ! solar local time (hours)
          lt = real(tod,kind=r8)*tod2hrs + state%lon(i)*rad2hrs

          ! lunar local time
          lun_lt = lt - nu*deg2hrs ! hours
          lun_lt = lun_lt*hrs2rad ! radians

          do k=1,pver
             ! Calculate the M2 lunar tide forcing in the zonal and meridional directions.
             ! The forcing is calculated based on the gradient of the M2 tidal
             ! potential, which is given in Chapman and Lindzen (1970).
             ! Additional details on the derivation of the forcing are in 
             ! Pedatella, Liu, and Richmond (2012) 
             ptend%u(i,k) = (-1._r8/((state%zm(i,k)+rearth)*cos(state%lat(i))))*2.456_r8*3._r8 *  &
                  ((state%zm(i,k)+rearth)/rearth)**2*cos(state%lat(i))*cos(state%lat(i))*2._r8*sin(2._r8*lun_lt)
             ptend%v(i,k) = (1._r8/(state%zm(i,k)+rearth))*2.456_r8*3._r8 * &
                  ((state%zm(i,k)+rearth)/rearth)**2*cos(2._r8*lun_lt)*2._r8*cos(state%lat(i))*sin(state%lat(i))
          end do
       end do
       
       call outfld('UT_LUNAR', ptend%u(:state%ncol,:), state%ncol, state%lchnk)
       call outfld('VT_LUNAR', ptend%v(:state%ncol,:), state%ncol, state%lchnk)

    end if

  end subroutine lunar_tides_tend

end module lunar_tides
