module budgets

! Metadata manager for the budgets.

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog
use cam_thermo,       only: thermo_budget_num_vars

implicit none
private
save

interface budget_add
  module procedure budget_stage_add
  module procedure budget_diff_add
end interface budget_add

!interface budget_info
! module procedure abudget_info_byind
!  module procedure budget_info_byname
!end interface budget_info

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
   budget_info_byname,    &! return budget info by name
   budget_cnt_adjust,     &! advance or reset budget count
   budget_count,          &! return budget count
   is_budget,             &! return budget count
   budget_get_global,     &! return budget count
   budget_put_global,     &! return budget count
   budget_write,          &! write_budget: time to write global budget

   budget_readnl,         &! budget_readnl: read cam thermo namelist
   budget_outfld           ! Returns true if default CAM output was specified in the budget_stage_add calls.


! Public data

integer, parameter, public :: budget_array_max  = 100     ! number of budget diffs
integer,           public            :: budget_cnt(budget_array_max)      ! budget counts for normalization
logical,           public            :: budget_subcycle(budget_array_max) ! budget_subcycle counts
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
integer,           public, protected :: budget_stg1stateidx(budget_array_max)
integer,           public, protected :: budget_stg2stateidx(budget_array_max)
real(r8),          public, protected :: budget_globals(budget_array_max,thermo_budget_num_vars)

integer,           public, protected :: thermo_budget_averaging_n = 1
integer,           public, protected :: thermo_budget_histfile_num = 1
logical,           public, protected :: thermo_budget_history = .false.
character(len=8),  public, protected :: thermo_budget_averaging_option = 'NONE'

!
! Constants for each budget

!character*3, public, protected :: budget_type(budget_array_max)! stage or difference
character*3, public :: budget_optype(budget_array_max)! stage or difference or sum
character*3, public :: budget_pkgtype(budget_array_max)! phy or dyn

!==============================================================================================
CONTAINS
!==============================================================================================

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

 subroutine budget_info(budget_ind, name, longname, stg1name, stg1stateidx, stg1index, stg2name, stg2stateidx, stg2index, optype, pkgtype,state_ind,subcycle,outfld)

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

 end subroutine budget_info

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
  budget_subcycle(:) = .false.
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

 logical function budget_write (step_offset)
   !
   !-----------------------------------------------------------------------
   !
   ! Purpose: Set flags that will initiate dump to IC file when OUTFLD and
   ! WSHIST are called
   !
   !-----------------------------------------------------------------------
   !
   use shr_kind_mod,   only: r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
   use shr_string_mod, only: shr_string_toUpper
   use time_manager,   only: timemgr_time_ge, timemgr_time_inc, get_curr_date, is_first_restart_step
   use time_manager,   only: get_step_size, get_nstep, is_last_step, is_first_step
   use time_manager,   only: get_start_date, get_stop_date
   use cam_abortutils, only: endrun
   use spmd_utils,     only: masterproc, mstrid=>masterprocid, mpicom
   use spmd_utils,     only: mpi_integer, mpi_real8, mpi_logical, mpi_character
   use cam_logfile,    only: iulog
   use shr_cal_mod,    only: shr_cal_ymd2date
   !
   ! Input/Output arguments
   !-----------------------
   integer, optional,    intent(in)  :: step_offset

   ! Local values
   !----------------
   character(len=*), parameter :: subname = 'budget_write :: '

   integer, save :: YMD_Next,Sec_Next, &
        YMD_Start,Sec_Start,YMD_Stop,Sec_Stop
   logical, save :: initialized=.false.
   integer :: YMD,Sec,YMD_Curr,Sec_Curr,YMD_Curr_woff,Sec_Curr_woff
   integer :: Year,Month,Day
   integer :: dtime             ! timestep size
   integer :: nstep             ! current timestep number
   integer :: offset            ! offset for writing thermo budget.
   logical :: Update_Budget

   !--------------------------------------------------------------

   budget_write  = .false.
   if (trim(shr_string_toUpper(thermo_budget_averaging_option)) == 'NONE') return

   offset=0
   if (present(step_offset)) offset=step_offset

   nstep = get_nstep()
   dtime = get_step_size()

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec_Curr)
   call shr_cal_ymd2date(Year,Month,Day,YMD_Curr)

   call get_curr_date(Year,Month,Day,Sec_Curr_woff,offset=offset)
   call shr_cal_ymd2date(Year,Month,Day,YMD_Curr_woff)

   if (masterproc) write(iulog,*)'budget_write YMD_Curr, Sec_Curr, offset',YMD_Curr,Sec_Curr,offset

   ! Initialize budget update times on first step
   if (.not. initialized) then
      ! Get Start time
      !--------------------
      call get_start_date(Year,Month,Day,Sec_Start)
      call shr_cal_ymd2date(Year,Month,Day,YMD_Start)

      ! Get End time
      !--------------------
      call get_stop_date(Year,Month,Day,Sec_Stop)
      call shr_cal_ymd2date(Year,Month,Day,YMD_Stop)

      ! Get Next Update time
      !--------------------
      if (thermo_budget_averaging_option == 'ENDOFRUN') then
         YMD_Next=YMD_Stop
         Sec_Next=Sec_Stop
      else
         YMD=YMD_Curr
         Sec=Sec_Curr
         if (thermo_budget_averaging_option == 'NSTEP') then
            call timemgr_time_inc(YMD,Sec,              &
              YMD_Next,Sec_Next,inc_s=dtime*thermo_budget_averaging_n)
         elseif (thermo_budget_averaging_option == 'NHOUR') then
            call timemgr_time_inc(YMD,Sec,              &
              YMD_Next,Sec_Next,inc_h=thermo_budget_averaging_n)
         elseif(thermo_budget_averaging_option == 'NDAY'   ) then
            call timemgr_time_inc(YMD,Sec,              &
              YMD_Next,Sec_Next,inc_d=thermo_budget_averaging_n)
         elseif(thermo_budget_averaging_option == 'NMONTH'   ) then
            call get_curr_date(Year,Month,Day,Sec_Curr)
            if (thermo_budget_averaging_n+Month.gt.12) then
               Year=Year+(thermo_budget_averaging_n+Month)/12
               Month=mod(thermo_budget_averaging_n+Month,12)
            else
               Month=thermo_budget_averaging_n+Month
            end if
            call shr_cal_ymd2date(Year,Month,Day,YMD_Next)
            Sec_Next=Sec_Curr
         elseif(thermo_budget_averaging_option == 'NYEAR'   ) then
            call get_curr_date(Year,Month,Day,Sec_Curr)
            call shr_cal_ymd2date(Year+thermo_budget_averaging_n,Month,Day,YMD_Next)
            Sec_Next=Sec_Curr
         end if

         if (masterproc) write(iulog,*)'init calc of next budget write ymdc/secc/ymdn/secn:',YMD_Curr,Sec_Curr,YMD_Next,Sec_Next
      end if

      initialized=.true.
   end if

   
   ! If an offset is present don't reset YMD_Next,Sec_Next just return budget_write using offset
   !--------------------------------------------------------------
   if (present(step_offset)) then

      call timemgr_time_ge(YMD_Next,Sec_Next,            &
           YMD_Curr_woff ,Sec_Curr_woff           ,update_budget)
!jt      budget_write  = ( nstep+(abs(offset/dtime))==1 .or. ((nstep /= 0).and.update_budget) )
      budget_write  = ((nstep /= 0).and.update_budget)
      if (masterproc) write(iulog,*)'checking for budget_write w/offset:',budget_write

   else
      ! When past the NEXT time, Update budget
      !--------------------------------------------------------------
      call timemgr_time_ge(YMD_Next,Sec_Next,            &
           YMD_Curr ,Sec_Curr           ,Update_Budget)
      if (masterproc) write(iulog,*)'checking for update_budget:',Update_Budget
      
      ! Reset YMD_Next and Sec_Next for next update
      !--------------------------------------------------------------
      if (Update_Budget) then
         if (thermo_budget_averaging_option == 'ENDOFRUN') then
            YMD_Next=YMD_Stop
            Sec_Next=Sec_Stop
         else
            YMD=YMD_Next
            Sec=Sec_Next
            if (thermo_budget_averaging_option == 'NSTEP') then
               call timemgr_time_inc(YMD,Sec,              &
                    YMD_Next,Sec_Next,inc_s=dtime*thermo_budget_averaging_n)
            elseif (thermo_budget_averaging_option == 'NHOUR') then
               call timemgr_time_inc(YMD,Sec,              &
                    YMD_Next,Sec_Next,inc_h=thermo_budget_averaging_n)
            elseif(thermo_budget_averaging_option == 'NDAY'   ) then
               call timemgr_time_inc(YMD,Sec,              &
                    YMD_Next,Sec_Next,inc_d=thermo_budget_averaging_n)
            elseif(thermo_budget_averaging_option == 'NMONTH'   ) then
               call get_curr_date(Year,Month,Day,Sec_Curr)
               if (thermo_budget_averaging_n+Month.gt.12) then
                  Year=Year+(thermo_budget_averaging_n+Month)/12
                  Month=mod(thermo_budget_averaging_n+Month,12)
               else
                  Month=thermo_budget_averaging_n+Month
               end if
               call shr_cal_ymd2date(Year,Month,Day,YMD_Next)
               Sec_Next=Sec_Curr
            elseif(thermo_budget_averaging_option == 'NYEAR'   ) then
               call get_curr_date(Year,Month,Day,Sec_Curr)
               call shr_cal_ymd2date(Year+thermo_budget_averaging_n,Month,Day,YMD_Next)
               Sec_Next=Sec_Curr
            end if
            if (masterproc) write(iulog,*)'curr gt next, reset next,new values ymdn/secn',YMD_Next,Sec_Next
         end if
      end if
!jt      budget_write  = ( nstep+(abs(offset/dtime))==1 .or. ((nstep /= 0).and.update_budget) )
      budget_write  = ((nstep /= 0).and.update_budget)
   end if

   return
 end function budget_write

 !===========================================================================
 ! Read namelist variables.
 subroutine budget_readnl(nlfile)
   use namelist_utils,  only: find_group_name
   use spmd_utils,      only: masterproc, mpicom, masterprocid
   use spmd_utils,      only: mpi_character, mpi_logical, mpi_real8, mpi_integer
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
