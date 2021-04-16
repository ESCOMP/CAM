module utils_mod
  use shr_kind_mod   ,only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use cam_logfile    ,only: iulog
  use cam_abortutils ,only: endrun
  use esmf           ,only: ESMF_FIELD
  use edyn_mpi       ,only: mlon0,mlon1,mlat0,mlat1, lon0,lon1,lat0,lat1
  use edyn_params    ,only: finit

  implicit none
  private

  public :: boxcar_ave
  public :: check_ncerr
  public :: check_alloc

contains

  !-----------------------------------------------------------------------
  subroutine boxcar_ave(x,y,lon,lat,mtime,itime,ibox)
    !
    ! perform boxcar average
    !
    ! Args:
    integer,  intent(in)  :: lon
    integer,  intent(in)  :: lat
    integer,  intent(in)  :: mtime
    integer,  intent(in)  :: itime
    integer,  intent(in)  :: ibox
    real(r8), intent(in)  :: x(lon,lat,mtime)
    real(r8), intent(out) :: y(lon,lat)

    ! Local:
    integer :: i, iset, iset1

    if (ibox > mtime) then
       call endrun('boxcar_ave: ibox > mtime')
    endif
    !
    iset = itime - ibox/2
    if (iset < 1) iset = 1
    iset1 = iset + ibox
    if (iset1 > mtime) then
       iset1 = mtime
       iset = iset1 - ibox
    end if
    y(:,:) = 0._r8
    do i=iset,iset1
       y(:,:) = y(:,:) + x(:,:,i)
    end do
    if (ibox > 0) y(:,:) = y(:,:)/ibox
    !
  end subroutine boxcar_ave

  !-----------------------------------------------------------------------
  subroutine check_alloc(ierror, subname, varname, lonp1, latp1, ntimes, lw)
    use spmd_utils, only: masterproc
    integer,           intent(in) :: ierror
    character(len=*),  intent(in) :: subname
    character(len=*),  intent(in) :: varname
    integer, optional, intent(in) :: lonp1
    integer, optional, intent(in) :: latp1
    integer, optional, intent(in) :: ntimes
    integer, optional, intent(in) :: lw
    ! Local variable
    character(len=cl) :: errmsg

    if (ierror /= 0) then
       write(errmsg, '(">>> ",a,": error allocating ",a)')                   &
            trim(subname), trim(varname)
       if (present(lonp1)) then
          write(errmsg(len_trim(errmsg)+1:), '(", lonp1 = ",i0)') lonp1
       end if
       if (present(latp1)) then
          write(errmsg(len_trim(errmsg)+1:), '(", latp1 = ",i0)') latp1
       end if
       if (present(ntimes)) then
          write(errmsg(len_trim(errmsg)+1:), '(", ntimes = ",i0)') ntimes
       end if
       if (present(lw)) then
          write(errmsg(len_trim(errmsg)+1:), '(", lw = ",i0)') lw
       end if
       if (masterproc) then
          write(iulog, *) trim(errmsg)
       end if
       call endrun(trim(errmsg))
    end if

  end subroutine check_alloc

  !-----------------------------------------------------------------------
  subroutine check_ncerr(istat, subname, msg)
    use pio, only: pio_noerr
    !
    ! Handle a netcdf lib error:
    !
    integer,          intent(in) :: istat
    character(len=*), intent(in) :: subname
    character(len=*), intent(in) :: msg
    !
    ! Local variable
    character(len=cl) :: errmsg
    !
    if (istat /= pio_noerr) then
       write(iulog,"(/72('-'))")
       write(iulog,"('>>> Error from netcdf library:')")
       write(iulog,"(a,': Error getting ',a)") trim(subname), trim(msg)

       write(iulog,"('istat=',i5)") istat
       write(iulog,"(72('-')/)")
       write(errmsg, '("NetCDF Error in ",a,": ",2a,", istat = ",i0)')        &
            trim(subname), 'Error getting ', trim(msg), istat
       call endrun(trim(errmsg))
    end if
  end subroutine check_ncerr

end module utils_mod
