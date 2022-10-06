
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

interface budget_add
  module procedure budget_stage_add
  module procedure budget_diff_add
end interface budget_add

interface budget_info
  module procedure budget_info_byind
  module procedure budget_info_byname
end interface budget_info

! Public interfaces
public :: &
   budget_init,           &! initialize budget variables
   budget_add,            &! add a budget to the list of budgets
   budget_num_avail,      &! returns the number of available slots in the budget array
   budget_chk_dim,        &! check that number of budgets added equals dimensions (budget_array_max)
   budget_name_byind,     &! return name of a budget
   budget_ind_byname,     &! return budget index given name
   budget_longname_byind, &! return longnamee of a budget
   budget_type_byind,     &! return stage or difference type of a budget
   budget_info,           &! return budget info by ind
   budget_cnt_adjust,     &! advance or reset budget count
   budget_count,          &! return budget count 
   is_budget,             &! return budget count 
   budget_get_global,     &! return budget count 
   budget_put_global,     &! return budget count 
   budget_outfld           ! Returns true if default CAM output was specified in the budget_stage_add calls.

! Public data

integer, parameter, public :: budget_array_max  = 60     ! number of budget diffs
integer, parameter, public :: budget_me_varnum  =  7     ! tot,se,ke,wv,wl,wi

integer,           public            :: budget_cnt(budget_array_max)      ! budget counts for normalization
integer,           public            :: budget_subcycle(budget_array_max) ! budget_subcycle counts
integer,           public            :: budget_num_dyn = 0 ! 
integer,           public            :: budget_num_phy = 0 ! 
integer,           public            :: budget_num     = 0 ! 
integer,           public            :: budget_state_ind(budget_array_max)      ! 
logical,           public, protected :: budget_out(budget_array_max)      ! outfld this stage
character(len=64), public, protected :: budget_name(budget_array_max)     ! budget names
character(len=128),public, protected :: budget_longname(budget_array_max) ! long name of budgets
character(len=128),public, protected :: budget_stagename(budget_array_max) ! long name of budgets
integer,           public, protected :: budget_stg1index(budget_array_max)
integer,           public, protected :: budget_stg2index(budget_array_max)
character(len=64), public, protected :: budget_stg1name(budget_array_max)
character(len=64), public, protected :: budget_stg2name(budget_array_max)
character(len=64), public, protected :: budget_me_names(budget_me_varnum)
integer,           public, protected :: budget_stg1stateidx(budget_array_max)
integer,           public, protected :: budget_stg2stateidx(budget_array_max)
real(r8),          public, protected :: budget_globals(budget_array_max,budget_me_varnum)

!
! Constants for each budget

!character*3, public, protected :: budget_type(budget_array_max)! stage or difference
character*3, public :: budget_optype(budget_array_max)! stage or difference or sum
character*3, public :: budget_pkgtype(budget_array_max)! phy or dyn

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


end subroutine budget_readnl


subroutine budget_stage_add (name, pkgtype, longname, outfld, subcycle)

   ! Register a budget.

   character(len=*), intent(in) :: &
      name      ! budget name used as variable name in history file output (8 char max)
   character(len=*), intent(in) :: &
      pkgtype      ! budget type either phy or dyn

   character(len=*), intent(in), optional :: &
      longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)
   logical,          intent(in), optional :: &
      outfld  ! true => default CAM output of budget in kg/kg
   logical,          intent(in), optional :: &
      subcycle  ! true => This budget is subcycled
   integer  :: state_idx    ! dyn/phy state budget index (in q array)
   character(len=*), parameter :: sub='budget_stage_add'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   ! set budget index and check validity
   if (pkgtype=='phy') then
      budget_num_phy = budget_num_phy+1
      state_idx  = budget_num_phy
   else if (pkgtype=='dyn') then
      budget_num_dyn = budget_num_dyn+1
      state_idx  = budget_num_dyn
   else
      call endrun('unknown budget pkgtype')
   end if
   budget_num = budget_num+1

   if (budget_num > budget_array_max) then
      write(errmsg, *) sub//': FATAL: budget stage index greater than budget stage max=', budget_array_max
      call endrun(errmsg)
   end if

   ! set budget name and constants
   budget_name(budget_num) = name
   if (present(longname)) then
      budget_longname(budget_num) = longname
   else
      budget_longname(budget_num) = name
   end if

   ! set outfld type 
   ! (false: the module declaring the budget is responsible for outfld calls)
   if (present(outfld)) then
      budget_out(budget_num) = outfld
   else
      budget_out(budget_num) = .false.
   end if
   budget_optype(budget_num)='stg'
   budget_pkgtype(budget_num)=pkgtype
   budget_state_ind(budget_num)=state_idx
   if (present(subcycle)) then
      budget_subcycle(budget_num)=subcycle
   else
      budget_subcycle(budget_num)=.false.
   end if
   budget_stagename(budget_num)= trim(name)

end subroutine budget_stage_add

!!$!==============================================================================

subroutine budget_diff_add (name, stg1name, stg2name, pkgtype, optype, longname, outfld, subcycle)

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
      outfld  ! true => default CAM output of budget in kg/kg

   logical,          intent(in), optional :: &
      subcycle  ! true => if this budget is subcycled

   character(len=*), parameter :: sub='budget_diff_add'
   character(len=128) :: errmsg
   character(len=1)   :: opchar
   integer            :: state_idx
   !-----------------------------------------------------------------------
   ! set budget index and check validity
   if (pkgtype=='phy') then
      budget_num_phy=budget_num_phy+1
      state_idx = budget_num_phy
   else if (pkgtype=='dyn') then
      budget_num_dyn=budget_num_dyn+1
      state_idx = budget_num_dyn
   else
      call endrun('bad budget pkgtype')
   end if
   budget_num= budget_num+1
   budget_pkgtype(budget_num)=pkgtype
   if (budget_num > budget_array_max) then
      write(errmsg, *) sub//': FATAL: budget diff index:',budget_num,' greater than budget_array_max=', budget_array_max
      call endrun(errmsg)
   end if

   ! set budget name and constants
   budget_name(budget_num) = name
   if (present(longname)) then
      budget_longname(budget_num) = longname
   else
      budget_longname(budget_num) = name
   end if
   if (optype=='dif') opchar='-'
   if (optype=='sum') opchar='+'
   if (optype=='stg') then
      write(errmsg,*) sub//': FATAL: bad value optype should be sum of dif:', optype
      call endrun(errmsg)
   end if
   budget_stg1name(budget_num) = trim(stg1name)
   budget_stg2name(budget_num) = trim(stg2name)
   budget_stagename(budget_num)= trim(stg1name)//opchar//trim(stg2name)
   budget_stg1index(budget_num) = budget_ind_byname(trim(stg1name))
   budget_stg2index(budget_num) = budget_ind_byname(trim(stg2name))
   budget_stg1stateidx(budget_num) = budget_state_ind(budget_stg1index(budget_num))
   budget_stg2stateidx(budget_num) = budget_state_ind(budget_stg2index(budget_num))
   ! set outfld type 
   ! (false: the module declaring the budget is responsible for outfld calls)
   if (present(outfld)) then
      budget_out(budget_num) = outfld
   else
      budget_out(budget_num) = .false.
   end if

   budget_optype(budget_num)=optype
   budget_state_ind(budget_num)=state_idx
   if (present(subcycle)) then
      budget_subcycle(budget_num)=subcycle
   else
      budget_subcycle(budget_num)=.false.
   end if

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
      budget_type_byind = budget_optype(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for budget index=', ind
      call endrun(errmsg)
   end if

end function budget_type_byind

!==============================================================================================

subroutine budget_info_byname(name, budget_ind, longname, stg1name, stg1stateidx, stg1index, stg2name, stg2stateidx, stg2index, optype, pkgtype,state_ind,subcycle,outfld)

   ! Return the mixing ratio name of a budget 

   !-----------------------------Arguments---------------------------------
   character(len=*), intent(in) :: name
   character(len=*), intent(out), optional :: &
      longname,  &! budget long_name
      stg1name,  &! stage1 name value for difference budget
      stg2name    ! stage2 name value for difference budget
   integer, intent(out), optional :: &
      budget_ind,   &! budget array index
      state_ind,    &! state budget array index 
      stg1stateidx, &! stage1 index for difference budget
      stg2stateidx, &! stage2 index for difference budget
      stg1index,    &! stage1 budget index
      stg2index      ! stage2 budget index
   character(len=3), intent(out), optional :: &
      optype,      &! budget type difference or stage
      pkgtype     ! physics or dynamics budget
   logical, intent(out), optional :: &
      subcycle,    &!
      outfld

   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='budget_info_byname'
   character(len=128) :: errmsg
   integer :: b_ind
   !-----------------------------------------------------------------------
   b_ind=budget_ind_byname(trim(name))
   if (b_ind > 0 .and. b_ind <= budget_array_max) then
      if (present(budget_ind)) budget_ind=b_ind
      if (present(longname)) longname=budget_longname(b_ind)
      if (present(optype)) optype=budget_optype(b_ind)
      if (present(pkgtype)) pkgtype=budget_pkgtype(b_ind)
      if (present(state_ind)) state_ind=budget_state_ind(b_ind)
      if (present(subcycle)) subcycle=budget_subcycle(b_ind)
      if (present(outfld)) outfld=budget_out(b_ind)
      if (budget_optype(b_ind)=='dif' .or. budget_optype(b_ind)=='sum') then
         if (present(stg1name))stg1name=budget_stg1name(b_ind)
         if (present(stg2name))stg2name=budget_stg2name(b_ind)
         if (present(stg1stateidx)) stg1stateidx=budget_stg1stateidx(b_ind)
         if (present(stg2stateidx)) stg2stateidx=budget_stg2stateidx(b_ind)
         if (present(stg1index)) stg1index=budget_stg1index(b_ind)
         if (present(stg2index)) stg2index=budget_stg2index(b_ind)
      else
         if (present(stg1name).or.present(stg2name).or.present(stg1stateidx).or.present(stg2stateidx) &
             .or.present(stg1index).or.present(stg2index)) &
         call endrun(sub//': stage1/2 info not applicable for a budget that is not a difference or sum')
      end if
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value name:',name,' budget index=', b_ind
      call endrun(errmsg)
   end if
 end subroutine budget_info_byname

 subroutine budget_info_byind(budget_ind, name, longname, stg1name, stg1stateidx, stg1index, stg2name, stg2stateidx, stg2index, optype, pkgtype,state_ind,subcycle,outfld)

   ! Return the mixing ratio name of a budget 

   !-----------------------------Arguments---------------------------------
   integer, intent(in)  :: budget_ind
   character(len=*), intent(out), optional :: &
      name,      &! budget long_name
      longname,  &! budget long_name
      stg1name,  &! stage1 name value for difference budget
      stg2name    ! stage2 name value for difference budget
   integer, intent(out), optional :: &
      state_ind,   &! state budget array index 
      stg1stateidx,&! stage1 index for difference budget
      stg2stateidx,&! stage2 index for difference budget
      stg1index,   &! stage1 budget index
      stg2index     ! stage2 budget index
   character(len=3), intent(out), optional :: &
      optype,      &! budget type difference or stage
      pkgtype       ! physics or dynamics budget
   logical, intent(out), optional :: &
      subcycle,    &!
      outfld

   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='budget_info_byind'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------
   if (budget_ind > 0 .and. budget_ind <= budget_array_max) then
      if (present(outfld)) outfld=budget_out(budget_ind)
      if (present(name)) name=budget_name(budget_ind)
      if (present(longname)) longname=budget_longname(budget_ind)
      if (present(optype)) optype=budget_optype(budget_ind)
      if (present(pkgtype)) pkgtype=budget_pkgtype(budget_ind)
      if (present(state_ind)) state_ind=budget_state_ind(budget_ind)
      if (present(subcycle)) subcycle=budget_subcycle(budget_ind)
      if (budget_optype(budget_ind)=='dif' .or. budget_optype(budget_ind)=='sum') then
         if (present(stg1name))stg1name=budget_stg1name(budget_ind)
         if (present(stg2name))stg2name=budget_stg2name(budget_ind)
         if (present(stg1stateidx)) stg1stateidx=budget_stg1stateidx(budget_ind)
         if (present(stg2stateidx)) stg2stateidx=budget_stg2stateidx(budget_ind)
         if (present(stg1index)) stg1index=budget_stg1index(budget_ind)
         if (present(stg2index)) stg2index=budget_stg2index(budget_ind)
      else
         if (present(stg1name).or.present(stg2name).or.present(stg1stateidx).or.present(stg2stateidx) &
             .or.present(stg1index).or.present(stg2index)) &
         call endrun(sub//': stage1/2 info not applicable for a budget that is not a difference or sum')
      end if
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value name:',name,' budget index=', budget_ind
      call endrun(errmsg)
   end if

 end subroutine budget_info_byind

!==============================================================================================

subroutine budget_cnt_adjust(ind,reset)

   ! Return the mixing ratio name of a budget 

   !-----------------------------Arguments---------------------------------
   integer, intent(in)   :: ind    ! global budget index (in te array)
   logical, intent(in),optional   :: reset    ! reset budget_cnt
   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='budget_cnt_adjust'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------
   if (ind > 0 .and. ind <= budget_array_max) then
      budget_cnt(ind)=budget_cnt(ind)+1
      if (present(reset)) then
         if (reset) then
            budget_cnt(ind)=0
         end if
      end if
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad index value for budget_cnt_adjust=', ind
      call endrun(errmsg)
   end if


 end subroutine budget_cnt_adjust
!==============================================================================================
subroutine budget_init()

   ! Initial budget module variables.
  
  budget_cnt(:) = 0._r8
  budget_subcycle(:) = 0._r8
  budget_num_dyn = 0
  budget_num_phy = 0
  budget_num = 0
  budget_state_ind(:) = 0  
  budget_out(:)  = .false. 
  budget_name(:)    = 'UNSET'
  budget_longname(:)= 'UNSET'
  budget_stg1index(:) = 0
  budget_stg2index(:) = 0
  budget_stg1name(:)= 'UNSET'
  budget_stg2name(:)= 'UNSET'
  budget_subcycle(:)= .false.
  
end subroutine budget_init
!==============================================================================================


character*64 function budget_name_byind(ind)

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

subroutine budget_get_global (name, me_idx, global, abort)

   ! Get the global integral of a budget.  Optional abort argument allows returning
   ! control to caller when budget name is not found.  Default behavior is
   ! to call endrun when name is not found.

   !-----------------------------Arguments---------------------------------
   character(len=*),  intent(in)  :: name    ! budget name
   integer,           intent(in)  :: me_idx  ! mass energy variable index
   real(r8),          intent(out) :: global  ! global budget index (in q array)
   logical, optional, intent(in)  :: abort   ! optional flag controlling abort

   !---------------------------Local workspace-----------------------------
   integer :: m                                   ! budget index
   logical :: abort_on_error
   character(len=*), parameter :: sub='budget_get_global'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   ! Find budget name in list
   do m = 1, budget_array_max
      if (trim(name) == trim(budget_stagename(m)).or.trim(name)==trim(budget_name(m))) then
         global  = budget_globals(m,me_idx)
         return
      end if
   end do

   ! Unrecognized name
   abort_on_error = .true.
   if (present(abort)) abort_on_error = abort
   
   if (abort_on_error) then
      write(errmsg,*) sub//': FATAL: name not found: ', trim(name)
      call endrun(errmsg)
   end if

 end subroutine budget_get_global
!==============================================================================
!==============================================================================
subroutine budget_put_global (name, me_idx, global, abort)

   ! store the global integral of a budget.  Optional abort argument allows returning
   ! control to caller when budget name is not found.  Default behavior is
   ! to call endrun when name is not found.

   !-----------------------------Arguments---------------------------------
   character(len=*),  intent(in)  :: name  ! budget name
   integer,           intent(in)  :: me_idx! mass energy variable index
   real(r8),          intent(out) :: global   ! global budget index (in q array)
   logical, optional, intent(in)  :: abort ! optional flag controlling abort

   !---------------------------Local workspace-----------------------------
   integer :: m                                   ! budget index
   logical :: abort_on_error
   character(len=*), parameter :: sub='budget_put_ind'
   !-----------------------------------------------------------------------

   ! Find budget name in list
   do m = 1, budget_array_max
      if (trim(name) == trim(budget_stagename(m)).or.trim(name)==trim(budget_name(m))) then
         budget_globals(m,me_idx) = global
         return
      end if
   end do

   ! Unrecognized name
   abort_on_error = .true.
   if (present(abort)) abort_on_error = abort
   
   if (abort_on_error) then
      call endrun(sub//': FATAL: name not found')
   end if

 end subroutine budget_put_global
!==============================================================================

subroutine budget_get_ind (name, budget_ind, abort)

   ! Get the index of a budget.  Optional abort argument allows returning
   ! control to caller when budget name is not found.  Default behavior is
   ! to call endrun when name is not found.

   !-----------------------------Arguments---------------------------------
   character(len=*),  intent(in)  :: name  ! budget name
   integer,           intent(out) :: budget_ind   ! global budget index (in q array)
   logical, optional, intent(in)  :: abort ! optional flag controlling abort

   !---------------------------Local workspace-----------------------------
   integer :: m                                   ! budget index
   logical :: abort_on_error
   character(len=*), parameter :: sub='budget_get_ind'
   !-----------------------------------------------------------------------

   ! Find budget name in list
   do m = 1, budget_array_max
      if (trim(name) == trim(budget_name(m)).or.trim(name)==trim(budget_stagename(m))) then
         budget_ind  = m
         return
      end if
   end do

   ! Unrecognized name
   abort_on_error = .true.
   if (present(abort)) abort_on_error = abort
   
   if (abort_on_error) then
      call endrun(sub//': FATAL: name not found')
   end if

end subroutine budget_get_ind
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
   character(len=*), parameter :: sub='budget_ind_byname'
   !-----------------------------------------------------------------------

   ! Find budget name in list

   budget_ind_byname  = -1
   do m = 1, budget_array_max
      if (trim(name) == trim(budget_name(m)).or.trim(name) == trim(budget_stagename(m))) then
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

subroutine budget_chk_dim

   ! Check that the number of registered budgets is budget_array_max
   ! Write budget list to log file.

   integer :: i, m
   character(len=*), parameter :: sub='budget_chk_dim'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'Budgets  list:'
      do i = 1, budget_num
         write(iulog,'(2x,i4,2x,a8,2x,a128)') i, trim(budget_name(i)), trim(budget_longname(i))
      end do
   end if

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
!==============================================================================
function budget_count(ind)

   ! Query whether default CAM outfld calls should be made.

   !----------------------------------------------------------------------- 
   integer, intent(in) :: ind                ! budget index

   integer             :: budget_count  ! true => use default CAM outfld calls

   character(len=*), parameter :: sub='budget_count'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (ind > 0 .and. ind <= budget_array_max) then
      budget_count = budget_cnt(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for budget diff index=', ind
      call endrun(errmsg)
   end if

 end function budget_count

!==============================================================================

end module budgets
