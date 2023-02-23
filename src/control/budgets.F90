module budgets

! Metadata manager for the budgets.

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog
use cam_thermo,       only: thermo_budget_vars, thermo_budget_vars_descriptor, &
                                 thermo_budget_vars_unit, thermo_budget_vars_massv, thermo_budget_num_vars
use cam_history,      only: addfld, add_default, horiz_only
use cam_history_support, only: max_fieldname_len,ptapes

implicit none
private
save

interface budget_add
  module procedure budget_stage_add
  module procedure budget_diff_add
end interface budget_add

! Public interfaces
public :: &
   budget_init,           &! initialize budget variables
   budget_add,            &! add a budget to the list of budgets
   budget_ind_byname,     &! return budget index given name
   budget_get_global,     &! return budget global
   budget_put_global,     &! put budget global
   budget_readnl,         &! budget_readnl: read cam thermo namelist
   is_budget               ! return logical if budget_defined


! Public data

integer, parameter, public           :: budget_array_max  = 500     ! number of budget diffs
integer,           public            :: budget_num     = 0 !
character(len=64), public, protected :: budget_name(budget_array_max)     ! budget names
character(len=128),public, protected :: budget_longname(budget_array_max) ! long name of budgets
character(len=128),public, protected :: budget_stagename(budget_array_max) ! long name of budgets
character(len=64), public, protected :: budget_stg1name(budget_array_max)
character(len=64), public, protected :: budget_stg2name(budget_array_max)

integer,           public            :: thermo_budget_averaging_n = 1
integer,           public            :: thermo_budget_histfile_num = 1
logical,           public            :: thermo_budget_history = .false.
character(len=8),  public            :: thermo_budget_averaging_option = 'NONE'
integer,           private           :: stepsize
!
! Constants for each budget

character*3, public :: budget_optype(budget_array_max)! stage or difference or sum
character*3, public :: budget_pkgtype(budget_array_max)! phy or dyn

!==============================================================================================
CONTAINS
!==============================================================================================
  
subroutine budget_stage_add (name, pkgtype, longname, cslam)
  use dycore,          only: dycore_is  

   character(len=*), intent(in)           :: &
      name      ! budget name used as variable name in history file output (8 char max)
   character(len=*), intent(in)           :: &
      pkgtype      ! budget type either phy or dyn
   character(len=*), intent(in), optional :: &
      longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)
   logical,          intent(in), optional :: &
      cslam  ! true => CSLAM used to transport mass tracers

   character (len=128)                    :: errmsg
   character (len=max_fieldname_len)      :: str1
   character (len=128)                    :: str2, str3
   logical                                :: thermo_budget_hist
   logical                                :: cslamtr        ! using cslam transport for mass tracers
   integer                                :: ivars
   character(len=*), parameter            :: sub='budget_stage_add'
  !-----------------------------------------------------------------------

   if (thermo_budget_history) then
   if (present(cslam)) then
      cslamtr=cslam
   else
      cslamtr = .false.
   end if
   do ivars=1, thermo_budget_num_vars

      write(str1,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(name))
      write(str2,*) TRIM(ADJUSTL(thermo_budget_vars_descriptor(ivars)))," ", &
           TRIM(ADJUSTL(longname))
      write(str3,*) TRIM(ADJUSTL(thermo_budget_vars_unit(ivars)))

      budget_num = budget_num+1
      ! set budget name and constants
      budget_name(budget_num) = trim(str1)
      if (present(longname)) then
         budget_longname(budget_num) = trim(str2)
      else
         budget_longname(budget_num) = trim(str1)
      end if
      
      budget_optype(budget_num)='stg'
      budget_pkgtype(budget_num)=pkgtype
      budget_stagename(budget_num)= trim(name)
      
      if (pkgtype=='phy') then
         call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)),gridname='physgrid')
      else
         if (dycore_is('SE')) then
            if (cslamtr .and. thermo_budget_vars_massv(ivars)) then
               call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)),gridname='FVM')
            else
               call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)),gridname='GLL')
            end if
         else if (dycore_is('MPAS')) then
            call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)),gridname='mpas_cell')
         else
            call endrun(sub//'budget_add is only supported for MPAS and SE dycores')
            call endrun(errmsg)
         end if
      end if
      call add_default(TRIM(ADJUSTL(str1)), thermo_budget_histfile_num, 'N') 
   end do
end if
 end subroutine budget_stage_add
 
!!$!==============================================================================

subroutine budget_diff_add (name, stg1name, stg2name, pkgtype, optype, longname, cslam)
  use dycore,          only: dycore_is  


   ! Register a budget.

   character(len=*), intent(in) :: &
      name,stg1name,stg2name   ! budget name used as variable name in history file output (8 char max)

   character(len=*), intent(in) :: &
      pkgtype    ! budget type either phy or dyn

   character(len=*), intent(in) :: &
      optype    !  dif (difference) or sum or stg (stage)

   character(len=*), intent(in), optional :: &
      longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)

   logical,          intent(in), optional :: &
      cslam       ! true => use cslam to transport mass variables

   character(len=*), parameter :: sub='budget_diff_add'
   character(len=128) :: errmsg
   character(len=1)   :: opchar
   character (len=256)         :: str1, str2, str3, strstg1, strstg2
   integer :: ivars
   logical :: cslamtr        ! using cslam transport for mass tracers
   !-----------------------------------------------------------------------

   if (thermo_budget_history) then
   if (present(cslam)) then
      cslamtr=cslam
   else
      cslamtr = .false.
   end if

! register history budget variables
   do ivars=1, thermo_budget_num_vars
      
      write(str1,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(name))
      write(strstg1,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(stg1name))
      write(strstg2,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(stg2name))
      write(str2,*) TRIM(ADJUSTL(thermo_budget_vars_descriptor(ivars)))," ", &
           TRIM(ADJUSTL(longname))
      write(str3,*) TRIM(ADJUSTL(thermo_budget_vars_unit(ivars)))
      
      
      budget_num = budget_num + 1
      budget_pkgtype(budget_num)=pkgtype
      
      ! set budget name and constants
      budget_name(budget_num) = trim(str1)
      if (present(longname)) then
         budget_longname(budget_num) = trim(str2)
      else
         budget_longname(budget_num) = trim(str1)
      end if
      if (optype=='dif') opchar='-'
      if (optype=='sum') opchar='+'
      if (optype=='stg') then
         write(errmsg,*) sub//': FATAL: bad value optype should be sum of dif:', optype
         call endrun(errmsg)
      end if
      budget_stg1name(budget_num) = trim(adjustl(strstg1))
      budget_stg2name(budget_num) = trim(adjustl(strstg2))
      budget_stagename(budget_num)= trim(adjustl(strstg1))//trim(opchar)//trim(adjustl(strstg2))
      budget_optype(budget_num)=optype
      
      
      
      if (pkgtype=='phy') then
         call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)), &
              gridname='physgrid',op=optype,op_f1name=TRIM(ADJUSTL(strstg1)),op_f2name=TRIM(ADJUSTL(strstg2)))
      else
         if (dycore_is('SE')) then
            if (cslamtr .and. thermo_budget_vars_massv(ivars)) then
               call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)), &
                    gridname='FVM',op=optype,op_f1name=TRIM(ADJUSTL(strstg1)),op_f2name=TRIM(ADJUSTL(strstg2)))
            else
               call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)), &
                    gridname='GLL',op=optype,op_f1name=TRIM(ADJUSTL(strstg1)),op_f2name=TRIM(ADJUSTL(strstg2)))
            end if
         else if (dycore_is('MPAS')) then
            call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'N', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)), &
                 gridname='mpas_cell',op=optype,op_f1name=TRIM(ADJUSTL(strstg1)),op_f2name=TRIM(ADJUSTL(strstg2)))
         else
            call endrun(sub//'budget_add is only supported for MPAS and SE dycores')
            call endrun(errmsg)
         end if
      end if
      call add_default(TRIM(ADJUSTL(str1)), thermo_budget_histfile_num, 'N') 
   end do
end if
 end subroutine budget_diff_add

!==============================================================================================

subroutine budget_init()
  use time_manager,         only:  get_step_size

  stepsize=get_step_size()

end subroutine budget_init
!==============================================================================

subroutine budget_get_global (name, me_idx, global)
  
  use cam_history,          only: get_field_properties
  use cam_history_support,  only: active_entry
  use cam_thermo,           only: thermo_budget_vars_massv
  
  ! Get the global integral of a budget.  Optional abort argument allows returning
  ! control to caller when budget name is not found.  Default behavior is
  ! to call endrun when name is not found.
  
  !-----------------------------Arguments---------------------------------
  character(len=*),  intent(in)  :: name    ! budget name
  integer,           intent(in)  :: me_idx  ! mass energy variable index
  real(r8),          intent(out) :: global  ! global budget index (in q array)
  
  !---------------------------Local workspace-----------------------------
  type (active_entry), pointer   :: tape(:) => null()          ! history tapes
  character (len=max_fieldname_len) :: str1
  character(len=128)             :: errmsg
  integer                        :: b_ind                      ! hentry index
  integer                        :: f(ptapes),ff               ! hentry index
  integer                        :: idx,pidx,midx              ! substring index for sum dif char
  integer                        :: m                          ! budget index
  logical                        :: found                      ! true if global integral found

  character(len=*), parameter    :: sub='budget_get_global'
  !-----------------------------------------------------------------------

  str1=''
  write(str1,*) TRIM(ADJUSTL(name))
  ! check for stagename short format (stg1//op/stg2) where stg1 is name without thermo string appended
  midx=index(str1, '-')
  pidx=index(str1, '+')
  idx=midx+pidx
  if (idx /= 0 .and. (midx==0 .or. pidx==0)) then
     write(str1,*) TRIM(ADJUSTL(thermo_budget_vars(me_idx)))//"_"//trim(adjustl(str1(1:idx)))// &
                   TRIM(ADJUSTL(thermo_budget_vars(me_idx)))//"_"//TRIM(ADJUSTL(str1(idx+1:)))
  end if
  b_ind=budget_ind_byname(trim(adjustl(str1)))

  if (idx>0 .and. budget_optype(b_ind) == 'stg') call endrun(sub//'FATAL not a difference budget but name contains + or - character')

  write(str1,*) TRIM(ADJUSTL(budget_name(b_ind)))

  ! Find budget name in list and return global value
  call get_field_properties(trim(adjustl(str1)), found, tape_out=tape, ff_out=ff, f_out=f)
  if (found.and.f(thermo_budget_histfile_num)>0) then
     call tape(thermo_budget_histfile_num)%hlist(f(thermo_budget_histfile_num))%get_global(global)
     if (.not. thermo_budget_vars_massv(me_idx)) global=global/stepsize
  else
     write(errmsg,*) sub//': FATAL: name not found: ', trim(name)
     call endrun(errmsg)
  end if
  
end subroutine budget_get_global
!==============================================================================
subroutine budget_put_global (name, me_idx, global)
  
  use cam_history,          only: get_field_properties
  use cam_history_support,  only: active_entry
  use cam_thermo,           only: thermo_budget_vars_massv

  ! Get the global integral of a budget.  Optional abort argument allows returning
  ! control to caller when budget name is not found.  Default behavior is
  ! to call endrun when name is not found.
  
  !-----------------------------Arguments---------------------------------
  character(len=*),  intent(in)  :: name    ! budget name
  integer,           intent(in)  :: me_idx  ! mass energy variable index
  real(r8),          intent(in)  :: global  ! global budget index (in q array)
  
  !---------------------------Local workspace-----------------------------
  type (active_entry), pointer :: tape(:) => null()          ! history tapes
  integer                      :: m                          ! budget index
  integer                      :: f(ptapes),ff               ! hentry index
  character(len=*), parameter  :: sub='budget_put_global'
  character(len=128)           :: errmsg
  character (len=128)          :: str1
  logical                      :: found                      ! true if global integral found
  real(r8)                     :: global_normalized
  !-----------------------------------------------------------------------
  
  ! append thermo field to stage name
  write(str1,*) TRIM(ADJUSTL(thermo_budget_vars(me_idx))),"_",TRIM(ADJUSTL(name))

  ! Find budget name in list and push global value to hentry
  call get_field_properties(trim(str1), found, tape_out=tape, ff_out=ff, f_out=f)
  if (found.and.f(thermo_budget_histfile_num)>0) then
     if (.not. thermo_budget_vars_massv(me_idx)) global_normalized=global/stepsize
     call tape(thermo_budget_histfile_num)%hlist(f(thermo_budget_histfile_num))%put_global(global_normalized)
  else
     write(errmsg,*) sub//': FATAL: name not found: ', trim(name)
     call endrun(errmsg)
  end if
  
end subroutine budget_put_global
!==============================================================================
function budget_ind_byname (name)

   ! Get the index of a budget.  Optional abort argument allows returning
   ! control to caller when budget name is not found.  Default behavior is
   ! to call endrun when name is not found.

   !-----------------------------Arguments---------------------------------
   character(len=*),  intent(in)  :: name  ! budget name

   !---------------------------Local workspace-----------------------------
   integer                        :: budget_ind_byname   ! function return
   integer                        :: m                   ! budget index
   character(len=*), parameter    :: sub='budget_ind_byname'
   !-----------------------------------------------------------------------
   ! Find budget name in list
   budget_ind_byname  = -1
   do m = 1, budget_array_max
      if (trim(adjustl(name)) == trim(adjustl(budget_name(m))).or.trim(adjustl(name)) == trim(adjustl(budget_stagename(m)))) then
         budget_ind_byname  = m
         return
      end if
   end do
   if (budget_ind_byname  == -1) then
      write(iulog,*)'ind_byname failed, name=',trim(name),'budget_name='
      call endrun()
   end if

!==============================================================================
 end function budget_ind_byname

 function is_budget(name)

   ! Get the index of a budget.  Optional abort argument allows returning
   ! control to caller when budget name is not found.  Default behavior is
   ! to call endrun when name is not found.

   !-----------------------------Arguments---------------------------------
   character(len=*),  intent(in)  :: name  ! budget name

   !---------------------------Local workspace-----------------------------
   logical                        :: is_budget           ! function return
   integer                        :: m                   ! budget index
   character(len=*), parameter :: sub='is_budget'
   !-----------------------------------------------------------------------

   ! Find budget name in list of defined budgets

   is_budget = .false.
   do m = 1, budget_array_max
      if (trim(name) == trim(budget_name(m)).or.trim(name) == trim(budget_stagename(m))) then
         is_budget = .true.
         return
      end if
   end do
 end function is_budget

 !===========================================================================
 ! Read namelist variables.
 subroutine budget_readnl(nlfile)
   use namelist_utils,  only: find_group_name
   use spmd_utils,      only: masterproc, mpicom, masterprocid
   use spmd_utils,      only: mpi_character, mpi_logical, mpi_integer
   use cam_logfile,     only: iulog
   use shr_string_mod, only: shr_string_toUpper

   ! Dummy argument: filepath for file containing namelist input
   character(len=*), intent(in) :: nlfile

   ! Local variables
   integer                     :: unitn, ierr
   integer,          parameter :: lsize = 76
   integer,          parameter :: fsize = 23
   character(len=*), parameter :: subname = 'budget_readnl :: '
   character(len=8)            :: period
   logical                     :: thermo_budgeting

   namelist /thermo_budget_nl/  thermo_budget_averaging_option, thermo_budget_averaging_n,    &
        thermo_budget_history, thermo_budget_histfile_num
   !-----------------------------------------------------------------------

   if (masterproc) then
      open(newunit=unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'thermo_budget_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, thermo_budget_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname//'ERROR reading namelist, thermo_budget_nl')
         end if
      end if
      close(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(thermo_budget_averaging_option, len(thermo_budget_averaging_option), mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: thermo_budget_averaging_option")
   call mpi_bcast(thermo_budget_history         , 1                                  , mpi_logical  , masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: thermo_budget_history")
   call mpi_bcast(thermo_budget_histfile_num    , 1                                  , mpi_integer  , masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: thermo_budget_histfile_num")
   call mpi_bcast(thermo_budget_averaging_n     , 1                                  , mpi_integer  , masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: thermo_budget_averaging_n")

   if (trim(shr_string_toUpper(thermo_budget_averaging_option)) == 'NONE') then
      thermo_budgeting=.false.
   else
      thermo_budgeting=.true.
   end if

   ! Write out thermo_budget options
   if (masterproc) then
      if (trim(thermo_budget_averaging_option) == 'NSTEP' ) then
         period='step'
      else if (trim(thermo_budget_averaging_option) == 'NHOUR' ) then
         period='hour'
      else if (trim(thermo_budget_averaging_option) == 'NDAY' ) then
         period='day'
      else if (trim(thermo_budget_averaging_option) == 'NMONTH' ) then
         period='month'
      else if (trim(thermo_budget_averaging_option) == 'NYEAR' ) then
         period='year'
      else
         period=''
      end if

      if (trim(thermo_budget_averaging_option) == 'ENDOFRUN' ) then
         write(iulog,*)'Thermo thermo_budgets will be written at the end of the run'
      else
         if (thermo_budget_averaging_n == 1) then
            write(iulog,*)'Thermo thermo_budgets will be written every ',period
         else
            write(iulog,*)'Thermo thermo_budgets will be written every ',thermo_budget_averaging_n,' ',trim(period)//'s'
         end if
      end if

      if(thermo_budget_history.and..not.thermo_budgeting) then
         write(iulog,*)subname//": FATAL: thermo_budget_averaging_option =",thermo_budget_averaging_option
         call endrun(subname//": FATAL: thermo_budget averaging option must not be set to NONE when requesting thermo_budget history output")
      end if

   end if
 end subroutine budget_readnl

end module budgets
