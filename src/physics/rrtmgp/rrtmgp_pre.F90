module rrtmgp_pre
 use ccpp_kinds, only: kind_phys

 public :: rrtmgp_pre_run

CONTAINS

!> \section arg_table_rrtmgp_pre_run Argument Table
!! \htmlinclude rrtmgp_pre_run.html
!!
  subroutine rrtmgp_pre_run(coszrs, nstep, iradsw, iradlw, irad_always, &
                  ncol, idxday, nday, idxnite, nnite, dosw, dolw, errmsg, errflg)
     use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
     real(kind_phys), dimension(:), intent(in) :: coszrs
     integer,            intent(in) :: nstep
     integer,            intent(in) :: iradsw
     integer,            intent(in) :: iradlw
     integer,            intent(in) :: irad_always
     integer,            intent(in) :: ncol
     integer,           intent(out) :: nday
     integer,           intent(out) :: nnite
     integer, dimension(:), intent(out) :: idxday
     integer, dimension(:), intent(out) :: idxnite
     logical,               intent(out) :: dosw
     logical,               intent(out) :: dolw
     character(len=*),  intent(out) :: errmsg
     integer,           intent(out) :: errflg

     ! Local variables
     integer :: idx

     ! Set error variables
     errflg = 0
     errmsg = ''

     ! Gather night/day column indices.
     nday = 0
     nnite = 0
     do idx = 1, ncol
        if ( coszrs(idx) > 0.0_kind_phys ) then
           nday = nday + 1
           idxday(nday) = idx
        else
           nnite = nnite + 1
           idxnite(nnite) = idx
        end if
     end do

     ! Determine if we're going to do longwave and/or shortwave this timestep
     dosw = nstep == 0  .or.  iradsw == 1                     &
                       .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1)   &
                       .or. nstep <= irad_always 

     dolw = nstep == 0  .or.  iradlw == 1                     &
                       .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1)   &
                       .or. nstep <= irad_always
  end subroutine rrtmgp_pre_run

end module rrtmgp_pre
