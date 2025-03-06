module rrtmgp_pre
 use ccpp_kinds, only: kind_phys
 use cam_logfile, only: iulog

 public :: rrtmgp_pre_run
 public :: radiation_do_ccpp

CONTAINS

!> \section arg_table_rrtmgp_pre_run Argument Table
!! \htmlinclude rrtmgp_pre_run.html
!!
  subroutine rrtmgp_pre_run(coszrs, nstep, dtime, iradsw, iradlw, irad_always, ncol, &
                  nextsw_cday, idxday, nday, idxnite, nnite, dosw, dolw, errmsg, errflg)
     use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
     use time_manager,         only: get_curr_calday
     real(kind_phys), dimension(:), intent(in) :: coszrs
     integer,            intent(in) :: dtime
     integer,            intent(in) :: nstep
     integer,            intent(in) :: iradsw
     integer,            intent(in) :: iradlw
     integer,            intent(in) :: irad_always
     integer,            intent(in) :: ncol
     integer,           intent(out) :: nday
     integer,           intent(out) :: nnite
     real(kind_phys),   intent(out) :: nextsw_cday
     integer, dimension(:), intent(out) :: idxday
     integer, dimension(:), intent(out) :: idxnite
     logical,               intent(out) :: dosw
     logical,               intent(out) :: dolw
     character(len=*),  intent(out) :: errmsg
     integer,           intent(out) :: errflg

     ! Local variables
     integer :: idx
     integer :: offset
     integer :: nstep_next
     logical :: dosw_next
     real(kind_phys) :: caldayp1

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
     call radiation_do_ccpp('sw', nstep, iradsw, irad_always, dosw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if
     call radiation_do_ccpp('lw', nstep, iradlw, irad_always, dolw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     ! Get time of next radiation calculation - albedos will need to be
     ! calculated by each surface model at this time
     nextsw_cday = -1._kind_phys
     dosw_next = .false.
     offset = 0
     nstep_next = nstep
     do while (.not. dosw_next)
        nstep_next = nstep_next + 1
        offset = offset + dtime
        call radiation_do_ccpp('sw', nstep_next, iradsw, irad_always, dosw_next, errmsg, errflg)
        if (errflg /= 0) then
           return
        end if
        if (dosw_next) then
           nextsw_cday = get_curr_calday(offset=offset)
        end if
     end do
     if(nextsw_cday == -1._kind_phys) then
        errflg = 1
        errmsg = 'next calendar day with shortwave calculation not found'
        return
     end if

     ! determine if next radiation time-step not equal to next time-step
     if (nstep >= 1) then
        caldayp1 = get_curr_calday(offset=int(dtime))
        if (caldayp1 /= nextsw_cday) nextsw_cday = -1._kind_phys
     end if
  end subroutine rrtmgp_pre_run

!================================================================================================

subroutine radiation_do_ccpp(op, nstep, irad, irad_always, radiation_do, errmsg, errflg)

   ! Return true if the specified operation is done this timestep.

   character(len=*),  intent(in) :: op             ! name of operation
   integer,           intent(in) :: nstep
   integer,           intent(in) :: irad
   integer,           intent(in) :: irad_always
   integer,          intent(out) :: errflg
   character(len=*), intent(out) :: errmsg
   logical,          intent(out) :: radiation_do   ! return value

   !-----------------------------------------------------------------------

   ! Set error variables
   errflg = 0
   errmsg = ''

   select case (op)
      case ('sw') ! do a shortwave heating calc this timestep?
         radiation_do = nstep == 0  .or.  irad == 1                     &
                     .or. (mod(nstep-1,irad) == 0  .and.  nstep /= 1)   &
                     .or. nstep <= irad_always
      case ('lw') ! do a longwave heating calc this timestep?
         radiation_do = nstep == 0  .or.  irad == 1                     &
                     .or. (mod(nstep-1,irad) == 0  .and.  nstep /= 1)   &
                     .or. nstep <= irad_always
      case default
         errflg = 1
         errmsg = 'radiation_do_ccpp: unknown operation:'//op
   end select

end subroutine radiation_do_ccpp

end module rrtmgp_pre
