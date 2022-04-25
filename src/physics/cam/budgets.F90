
module budgets

! Metadata manager for the budgets.

use shr_kind_mod,     only: r8 => shr_kind_r8
use shr_const_mod,    only: shr_const_rgas
use spmd_utils,       only: masterproc
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

! Public interfaces
public :: &
   budget_stage_add,      &! add a budget to the list of budgets
   budget_diff_add,       &! add a budget to the list of budgets
   budget_num_avail,      &! returns the number of available slots in the budget array
   budget_get_ind,        &! get the index of a budget
   budget_chk_dim,        &! check that number of budgets added equals dimensions (budget_array_max)
   budget_name_byind,     &! return name of a budget
   budget_longname_byind, &! return longnamee of a budget
   budget_type_byind,     &! return stage or difference type of a budget
   budget_info_byind,     &! return stage or difference type of a budget
   budget_outfld           ! Returns true if default CAM output was specified in the budget_stage_add calls.

! Public data

integer, parameter, public :: budget_array_max  = 30     ! number of budget diffs

integer,           public            :: budget_cnt(budget_array_max)      ! outfld this stage
integer,           public            :: budget_num = 0  ! 
logical,           public, protected :: budget_out(budget_array_max)      ! outfld this stage
character(len=16), public, protected :: budget_name(budget_array_max)     ! budget names
character(len=128),public, protected :: budget_longname(budget_array_max) ! long name of budgets
integer,           public, protected :: budget_s1_ind(budget_array_max)
integer,           public, protected :: budget_s2_ind(budget_array_max)
character(len=16), public, protected :: budget_s1name(budget_array_max)
character(len=16), public, protected :: budget_s2name(budget_array_max)

!
! Constants for each budget

!character*3, public, protected :: budget_type(budget_array_max)! stage or difference
character*3, public :: budget_type(budget_array_max)! stage or difference

!++bee - temporary... These names should be declared in the module that makes the addfld and outfld calls.
! Lists of budget names and diagnostics
character(len=16), public :: physbudget    (budget_array_max)  ! budgets after physics  (FV core only)
character(len=16), public :: dynbudget    (budget_array_max)   ! budgets before physics (FV core only)

!==============================================================================================
CONTAINS
!==============================================================================================

subroutine budget_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical


   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: sub = 'budget_readnl'

   !-----------------------------------------------------------------------------

!!$   if (masterproc) then
!!$      unitn = getunit()
!!$      open( unitn, file=trim(nlfile), status='old' )
!!$      call find_group_name(unitn, 'budgets_nl', status=ierr)
!!$      if (ierr == 0) then
!!$         read(unitn, budgets_nl, iostat=ierr)
!!$         if (ierr /= 0) then
!!$            call endrun(sub//': FATAL: reading namelist')
!!$         end if
!!$      end if
!!$      close(unitn)
!!$      call freeunit(unitn)
!!$   end if

!!$   if (masterproc) then
!!$      write(iulog,*)'Summary of budget module options:'
!!$   end if

end subroutine budget_readnl


subroutine budget_stage_add (name, ind, longname, cam_outfld)

   ! Register a budget.

   character(len=*), intent(in) :: &
      name      ! budget name used as variable name in history file output (8 char max)
   integer, intent(out)   :: ind    ! global budget index (in q array)

   character(len=*), intent(in), optional :: &
      longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)
   logical,          intent(in), optional :: &
      cam_outfld  ! true => default CAM output of budget in kg/kg

   character(len=*), parameter :: sub='budget_stage_add'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   ! set budget index and check validity
   budget_num = budget_num+1
   ind  = budget_num
   if (budget_num > budget_array_max) then
      write(errmsg, *) sub//': FATAL: budget stage index greater than budget stage max=', budget_array_max
      call endrun(errmsg)
   end if

   ! set budget name and constants
   budget_name(ind) = name
   if (present(longname)) then
      budget_longname(ind) = longname
   else
      budget_longname(ind) = name
   end if

   ! set outfld type 
   ! (false: the module declaring the budget is responsible for outfld calls)
   if (present(cam_outfld)) then
      budget_out(ind) = cam_outfld
   else
      budget_out(ind) = .false.
   end if
   budget_type(ind)='stg'
end subroutine budget_stage_add

!==============================================================================
subroutine budget_diff_add (name, istage1, istage2, longname, cam_outfld)

   ! Register a budget.

   character(len=*), intent(in) :: &
      name      ! budget name used as variable name in history file output (8 char max)

   integer, intent(in)   :: istage1,istage2    ! global budget stage index (in te_budgets array)

   character(len=*), intent(in), optional :: &
      longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)

   logical,          intent(in), optional :: &
      cam_outfld  ! true => default CAM output of budget in kg/kg

   character(len=*), parameter :: sub='budget_diff_add'
   character(len=128) :: errmsg
   integer            :: ind
   !-----------------------------------------------------------------------

   ! set budget index and check validity
   budget_num = budget_num+1
   ind  = budget_num
   if (budget_num > budget_array_max) then
      write(errmsg, *) sub//': FATAL: budget diff index greater than budget_array_max=', budget_array_max
      call endrun(errmsg)
   end if

   ! set budget name and constants
   budget_name(ind) = name
   if (present(longname)) then
      budget_longname(ind) = longname
   else
      budget_longname(ind) = name
   end if

   budget_s1_ind(ind) = istage1
   budget_s2_ind(ind) = istage2
   budget_s1name(ind) = budget_name_byind(istage1)
   budget_s2name(ind) = budget_name_byind(istage2)

   ! set outfld type 
   ! (false: the module declaring the budget is responsible for outfld calls)
   if (present(cam_outfld)) then
      budget_out(ind) = cam_outfld
   else
      budget_out(ind) = .false.
   end if
   budget_type(ind)='dif'
 end subroutine budget_diff_add
!==============================================================================

function budget_num_avail()

   ! return number of available slots in the budget array

   integer budget_num_avail

   budget_num_avail = budget_array_max - budget_num

end function budget_num_avail

!==============================================================================================

character*3 function budget_type_byind(ind)

   ! Return the type of a budget stage or difference

   !-----------------------------Arguments---------------------------------
   integer, intent(in)   :: ind    ! global budget index (in te array)

   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='budget_type_byind'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------
   if (ind > 0 .and. ind <= budget_array_max) then
      budget_type_byind = budget_type(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for budget index=', ind
      call endrun(errmsg)
   end if

end function budget_type_byind

!==============================================================================================

subroutine budget_info_byind(ind, name, longname, stage1, istage1, stage2, istage2)

   ! Return the mixing ratio name of a budget 

   !-----------------------------Arguments---------------------------------
   integer, intent(in)   :: ind    ! global budget index (in te array)
   character(len=*), intent(out), optional :: &
      name,      &! budget name
      longname,  &! budget long_name
      stage1,    &! stage1 name value for difference budget
      stage2      ! stage2 name value for difference budget
   integer, intent(out), optional :: &
      istage1,   &! stage1 index for difference budget
      istage2     ! stage2 index for difference budget

   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='budget_name_byind'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------
   if (ind > 0 .and. ind <= budget_array_max) then
      if (present(name)) name=budget_name(ind)
      if (present(longname)) longname=budget_longname(ind)
      if (present(stage1))   stage1=budget_s1name(ind)
      if (present(stage2))   stage2=budget_s2name(ind)
      if (present(istage1))  istage1=budget_s1_ind(ind)
      if (present(istage2))  istage2=budget_s2_ind(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for budget index=', ind
      call endrun(errmsg)
   end if


 end subroutine budget_info_byind

!==============================================================================================

character*16 function budget_name_byind(ind)

   ! Return the mixing ratio name of a budget 

   !-----------------------------Arguments---------------------------------
   integer, intent(in)   :: ind    ! global budget index (in te array)
   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='budget_name_byind'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (ind > 0 .and. ind <= budget_array_max) then
      budget_name_byind = budget_name(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for budget index=', ind
      call endrun(errmsg)
   end if

end function budget_name_byind

!==============================================================================================

character*128 function budget_longname_byind(ind)

   ! Return the mixing ratio name of a budget 

   !-----------------------------Arguments---------------------------------
   integer, intent(in)   :: ind    ! global budget index (in te array)

   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='budget_name_byind'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (ind > 0 .and. ind <= budget_array_max) then
      budget_longname_byind = budget_longname(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for budget index=', ind
      call endrun(errmsg)
   end if

end function budget_longname_byind

!==============================================================================

subroutine budget_get_ind (name, ind, abort)

   ! Get the index of a budget.  Optional abort argument allows returning
   ! control to caller when budget name is not found.  Default behavior is
   ! to call endrun when name is not found.

   !-----------------------------Arguments---------------------------------
   character(len=*),  intent(in)  :: name  ! budget name
   integer,           intent(out) :: ind   ! global budget index (in q array)
   logical, optional, intent(in)  :: abort ! optional flag controlling abort

   !---------------------------Local workspace-----------------------------
   integer :: m                                   ! budget index
   logical :: abort_on_error
   character(len=*), parameter :: sub='budget_get_ind'
   !-----------------------------------------------------------------------

   ! Find budget name in list
   do m = 1, budget_array_max
      if (name == budget_name(m)) then
         ind  = m
         return
      end if
   end do

   ! Unrecognized name
   abort_on_error = .true.
   if (present(abort)) abort_on_error = abort

   if (abort_on_error) then
      write(iulog, *) sub//': FATAL: name:', name,  ' not found in list:', budget_name(:)
      call endrun(sub//': FATAL: name not found')
   end if

   ! error return
   ind = -1

end subroutine budget_get_ind

!==============================================================================


subroutine budget_chk_dim

   ! Check that the number of registered budgets is budget_array_max
   ! Write budget list to log file.

   integer :: i, m
   character(len=*), parameter :: sub='budget_chk_dim'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   ! if (budget_num /= budget_array_max) then
   !    write(errmsg, *) sub//': FATAL: number of added budgets (',budget_num, &
   !                     ') not equal to budget_array_max (', budget_array_max, ')'
   !    call endrun (errmsg)
   ! endif

   if (masterproc) then
      write(iulog,*) 'Budget stages list:'
      do i = 1, budget_num
         write(iulog,'(2x,i4,2x,a8,2x,a128)') i, budget_name(i), budget_longname(i)
      end do
      write(iulog,*) 'Budgets  list:'
      do i = 1, budget_num
         write(iulog,'(2x,i4,2x,a8,2x,a128)') i, budget_name(i), budget_longname(i)
      end do
   end if

   ! ! Set names of physics and dynamics budget
   ! do m=1,budget_array_max
   !    physbudget    (m)  = trim(budget_name(m))//'AP'
   !    dynbudget    (m)  = trim(budget_name(m))//'BP'
   ! end do

end subroutine budget_chk_dim

function budget_outfld(m)

   ! Query whether default CAM outfld calls should be made.

   !----------------------------------------------------------------------- 
   integer, intent(in) :: m                ! budget index

   logical             :: budget_outfld  ! true => use default CAM outfld calls

   character(len=*), parameter :: sub='budget_outfld'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (m > 0 .and. m <= budget_array_max) then
      budget_outfld = budget_out(m)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for budget diff index=', m
      call endrun(errmsg)
   end if

 end function budget_outfld

!==============================================================================

end module budgets
