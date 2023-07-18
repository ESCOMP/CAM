module table_interp_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: table_interp
  public :: table_interp_wghts
  public :: table_interp_updwghts

  interface table_interp
     module procedure interp2d
  end interface table_interp

  type :: table_interp_wghts
     real(r8) :: wt1
     real(r8) :: wt2
     integer  :: ix1
     integer  :: ix2
  end type table_interp_wghts

contains

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  pure function interp2d( ncoef,ncol,nxs,nys, xwghts,ywghts, tbl ) result(res)

    integer, intent(in)  :: ncoef,ncol,nxs,nys
    real(r8), intent(in) :: tbl(ncoef,nxs,nys)
    type(table_interp_wghts), intent(in) :: xwghts(ncol)
    type(table_interp_wghts), intent(in) :: ywghts(ncol)

    real(r8) :: res(ncoef,ncol)

    real(r8) :: fx(ncoef,2)

    integer :: i

    do i = 1,ncol

       fx(:,1) = xwghts(i)%wt1*tbl(:,xwghts(i)%ix1,ywghts(i)%ix1) &
               + xwghts(i)%wt2*tbl(:,xwghts(i)%ix2,ywghts(i)%ix1)
       fx(:,2) = xwghts(i)%wt1*tbl(:,xwghts(i)%ix1,ywghts(i)%ix2) &
               + xwghts(i)%wt2*tbl(:,xwghts(i)%ix2,ywghts(i)%ix2)

       res(:,i) = ywghts(i)%wt1*fx(:,1) + ywghts(i)%wt2*fx(:,2)

    end do


  end function interp2d

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  subroutine table_interp_updwghts( ngrid, xgrid, ncols, xcols, wghts )
    integer,  intent(in) :: ngrid
    real(r8), intent(in) :: xgrid(ngrid)
    integer,  intent(in) :: ncols
    real(r8), intent(in) :: xcols(ncols)
    type(table_interp_wghts), intent(inout) :: wghts(ncols)

    integer :: i

    do i = 1,ncols
       wghts(i)%ix2 = find_index(ngrid,xgrid,xcols(i))
       wghts(i)%ix1 = wghts(i)%ix2 - 1
       wghts(i)%wt1 = (xgrid(wghts(i)%ix2)-xcols(i)) &
                    /(xgrid(wghts(i)%ix2)-xgrid(wghts(i)%ix1))
       wghts(i)%wt2 = 1._8 - wghts(i)%wt1
    end do

  end subroutine table_interp_updwghts

  ! private methods
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  pure function find_index( nvals, vals, vx ) result(ndx)
    integer,  intent(in) :: nvals
    real(r8), intent(in) :: vals(nvals)
    real(r8), intent(in) :: vx

    integer :: ndx

    find_ndx: do ndx = 1, nvals-1
       if (vals(ndx)>vx) exit find_ndx
    end do find_ndx

  end function find_index

end module table_interp_mod
