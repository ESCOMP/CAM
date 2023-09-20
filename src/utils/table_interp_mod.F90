!----------------------------------------------------------------------------
! Utility module used for interpolation of aerosol optics table
!  NOTE: Results will be set to table edges for interpolations beyond
!        the edges -- no extropolations
!----------------------------------------------------------------------------
module table_interp_mod
  use shr_kind_mod, only: r8=>shr_kind_r8

  implicit none

  private
  public :: table_interp
  public :: table_interp_wghts
  public :: table_interp_calcwghts

  ! overload the interpolation routines
  interface table_interp
     module procedure interp1d
     module procedure interp2d
     module procedure interp4d
  end interface table_interp

  ! interpolation weights and indices
  type :: table_interp_wghts
     real(r8) :: wt1
     real(r8) :: wt2
     integer  :: ix1
     integer  :: ix2
  end type table_interp_wghts

contains

  !--------------------------------------------------------------------------
  ! 1-D interpolation
  !--------------------------------------------------------------------------
  pure function interp1d( ncol, nxs, xwghts, tbl ) result(res)

    integer, intent(in)  :: ncol                         ! number of model columns
    integer, intent(in)  :: nxs                          ! table size
    real(r8), intent(in) :: tbl(nxs)                     ! table values to be interpolated
    type(table_interp_wghts), intent(in) :: xwghts(ncol) ! interpolation weights and indices

    real(r8) :: res(ncol)

    integer :: i

    do i = 1,ncol

       res(i) = xwghts(i)%wt1*tbl(xwghts(i)%ix1) &
              + xwghts(i)%wt2*tbl(xwghts(i)%ix2)

    end do

  end function interp1d

  !--------------------------------------------------------------------------
  ! 2-D interpolation
  !--------------------------------------------------------------------------
  pure function interp2d( ncoef, ncol, nxs, nys, xwghts, ywghts, tbl ) result(res)

    integer, intent(in)  :: ncoef                        ! number chebyshev coefficients
    integer, intent(in)  :: ncol                         ! number of model columns
    integer, intent(in)  :: nxs                          ! table x-dimension size
    integer, intent(in)  :: nys                          ! table y-dimension size
    real(r8), intent(in) :: tbl(ncoef,nxs,nys)           ! table values to be interpolated
    type(table_interp_wghts), intent(in) :: xwghts(ncol) ! x interpolation weights and indices
    type(table_interp_wghts), intent(in) :: ywghts(ncol) ! y interpolation weights and indices

    real(r8) :: res(ncoef,ncol)

    real(r8) :: fx(ncoef,2)

    integer :: i

    do i = 1,ncol

       ! interp x dir
       fx(:,1) = xwghts(i)%wt1*tbl(:,xwghts(i)%ix1,ywghts(i)%ix1) & ! @ y1
               + xwghts(i)%wt2*tbl(:,xwghts(i)%ix2,ywghts(i)%ix1)
       fx(:,2) = xwghts(i)%wt1*tbl(:,xwghts(i)%ix1,ywghts(i)%ix2) & ! @ y2
               + xwghts(i)%wt2*tbl(:,xwghts(i)%ix2,ywghts(i)%ix2)

       ! interp y dir
       res(:,i) = ywghts(i)%wt1*fx(:,1) + ywghts(i)%wt2*fx(:,2)

    end do

  end function interp2d

  !--------------------------------------------------------------------------
  ! 4-D interpolation
  !--------------------------------------------------------------------------
  pure function interp4d( ncol, nxs, nys, nzs, nts, xwghts, ywghts, zwghts, twghts, tbl ) result(res)

    integer, intent(in)  :: ncol                         ! number of model columns
    integer, intent(in)  :: nxs                          ! table x-dimension size
    integer, intent(in)  :: nys                          ! table y-dimension size
    integer, intent(in)  :: nzs                          ! table z-dimension size
    integer, intent(in)  :: nts                          ! table t-dimension size
    real(r8), intent(in) :: tbl(nxs,nys,nzs,nts)         ! table values to be interpolated
    type(table_interp_wghts), intent(in) :: xwghts(ncol) ! x interpolation weights and indices
    type(table_interp_wghts), intent(in) :: ywghts(ncol) ! y interpolation weights and indices
    type(table_interp_wghts), intent(in) :: zwghts(ncol) ! z interpolation weights and indices
    type(table_interp_wghts), intent(in) :: twghts(ncol) ! t interpolation weights and indices

    real(r8) :: res(ncol)

    real(r8) :: fx(8)
    real(r8) :: fy(4)
    real(r8) :: fz(2)

    integer :: i

    do i = 1,ncol

       ! interp x dir
       fx(1) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix1,zwghts(i)%ix1,twghts(i)%ix1) & ! @ y1, z1, t1
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix1,zwghts(i)%ix1,twghts(i)%ix1)
       fx(2) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix2,zwghts(i)%ix1,twghts(i)%ix1) & ! @ y2, z1, t1
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix2,zwghts(i)%ix1,twghts(i)%ix1)

       fx(3) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix1,zwghts(i)%ix2,twghts(i)%ix1) & ! @ y1, z2, t1
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix1,zwghts(i)%ix2,twghts(i)%ix1)
       fx(4) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix2,zwghts(i)%ix2,twghts(i)%ix1) & ! @ y2, z2, t1
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix2,zwghts(i)%ix2,twghts(i)%ix1)

       fx(5) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix1,zwghts(i)%ix1,twghts(i)%ix2) & ! @ y1, z1, t2
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix1,zwghts(i)%ix1,twghts(i)%ix2)
       fx(6) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix2,zwghts(i)%ix1,twghts(i)%ix2) & ! @ y2, z1, t2
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix2,zwghts(i)%ix1,twghts(i)%ix2)

       fx(7) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix1,zwghts(i)%ix2,twghts(i)%ix2) & ! @ y1, z2, t2
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix1,zwghts(i)%ix2,twghts(i)%ix2)
       fx(8) = xwghts(i)%wt1*tbl(xwghts(i)%ix1,ywghts(i)%ix2,zwghts(i)%ix2,twghts(i)%ix2) & ! @ y2, z2, t2
             + xwghts(i)%wt2*tbl(xwghts(i)%ix2,ywghts(i)%ix2,zwghts(i)%ix2,twghts(i)%ix2)

       ! interp y dir
       fy(1) = ywghts(i)%wt1*fx(1) + ywghts(i)%wt2*fx(2) ! @ z1, t1
       fy(2) = ywghts(i)%wt1*fx(3) + ywghts(i)%wt2*fx(4) ! @ z2, t1
       fy(3) = ywghts(i)%wt1*fx(5) + ywghts(i)%wt2*fx(6) ! @ z1, t2
       fy(4) = ywghts(i)%wt1*fx(7) + ywghts(i)%wt2*fx(8) ! @ z2, t2

       ! interp z dir
       fz(1) = zwghts(i)%wt1*fy(1) + zwghts(i)%wt2*fy(2) ! @ t1
       fz(2) = zwghts(i)%wt1*fy(3) + zwghts(i)%wt2*fy(4) ! @ t2

       ! interp t dir
       res(i) = twghts(i)%wt1*fz(1) + twghts(i)%wt2*fz(2)

    end do

  end function interp4d

  !--------------------------------------------------------------------------
  ! determines interpolation weights and indices for given values at the model columns
  !--------------------------------------------------------------------------
  pure function table_interp_calcwghts( ngrid, xgrid, ncols, xcols ) result(wghts)

    integer,  intent(in) :: ngrid        ! number of grid point values
    real(r8), intent(in) :: xgrid(ngrid) ! grid point values
    integer,  intent(in) :: ncols        ! number of model columns
    real(r8), intent(in) :: xcols(ncols) ! values at the model columns

    type(table_interp_wghts) :: wghts(ncols) ! interpolations weights at the model columns

    integer :: i
    real(r8) :: xs(ncols)

    xs(:) = xcols(:)

    ! do not extrapolate beyond the edges of the table
    where(xs < xgrid(1))
       xs = xgrid(1)
    end where
    where(xs > xgrid(ngrid))
       xs = xgrid(ngrid)
    end where

    do i = 1,ncols
       wghts(i)%ix2 = find_index(ngrid,xgrid,xs(i))
       wghts(i)%ix1 = wghts(i)%ix2 - 1
       wghts(i)%wt1 = (xgrid(wghts(i)%ix2)-xs(i)) &
                     /(xgrid(wghts(i)%ix2)-xgrid(wghts(i)%ix1))
       wghts(i)%wt2 = 1._r8 - wghts(i)%wt1
    end do

  end function table_interp_calcwghts

  ! private methods
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  ! determines last index of grid vals of which is greater then or equal to
  ! value vx
  !--------------------------------------------------------------------------
  pure function find_index( nvals, vals, vx ) result(res)
    integer,  intent(in) :: nvals
    real(r8), intent(in) :: vals(nvals)
    real(r8), intent(in) :: vx
    integer :: res

    integer :: ndx

    res = -1

    find_ndx: do ndx = 2, nvals
       if (vals(ndx)>=vx) then
          res = ndx
          exit find_ndx
       end if
    end do find_ndx

  end function find_index

end module table_interp_mod
