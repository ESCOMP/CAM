module budgets
  !----------------------------------------------------------------------------
  !
  ! Adds support for energy and mass snapshots and budgets using cam_history api.
  !
  ! Public functions/subroutines:
  !
  ! budget_init
  ! e_m_snapshot
  ! e_m_budget
  ! budget_ind_byname
  ! budget_get_global
  ! budget_readnl
  ! is_budget
  !-----------------------------------------------------------------------

  use cam_abortutils,      only: endrun
  use cam_history,         only: addfld, add_default, horiz_only
  use cam_history_support, only: max_fieldname_len
  use cam_logfile,         only: iulog
  use cam_thermo,          only: thermo_budget_vars, thermo_budget_vars_descriptor, &
       thermo_budget_vars_unit, thermo_budget_vars_massv, thermo_budget_num_vars,teidx,wvidx,wlidx,wiidx
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use shr_kind_mod,        only: cl => shr_kind_cl
  use spmd_utils,          only: masterproc, masterprocid, mpicom

  implicit none
  private
  save

  ! Public interfaces
  public :: &
       budget_init,           &! initialize budget variables
       e_m_snapshot,          &! define a snapshot and add to history buffer
       e_m_budget,            &! define a budget and add to history buffer
       budget_ind_byname,     &! return budget index given name
       budget_get_global,     &! get global budget from history buffer
       budget_readnl,         &! read budget namelist setting
       is_budget               ! return logical if budget_defined


  ! Public data

  integer, parameter, public           :: budget_array_max  = 500     ! number of budget diffs

  integer,           public            :: budget_num     = 0 !
  character(cl),     public, protected :: budget_name(budget_array_max)     ! budget names
  character(cl),     public, protected :: budget_longname(budget_array_max) ! long name of budgets
  character(cl),     public, protected :: budget_stagename(budget_array_max) ! long name of budgets
  character(cl),     public, protected :: budget_stg1name(budget_array_max)
  character(cl),     public, protected :: budget_stg2name(budget_array_max)

  integer,           public            :: thermo_budget_histfile_num = 1
  logical,           public            :: thermo_budget_history = .false.
  real(r8),          private           :: dstepsize
  !
  ! Constants for each budget

  character*3, public :: budget_optype(budget_array_max)! stage or difference or sum
  character*3, public :: budget_pkgtype(budget_array_max)! phy or dyn

  !==============================================================================================
CONTAINS
  !==============================================================================================
  !
  ! Read namelist variables.
  subroutine budget_readnl(nlfile)
    use dycore,          only: dycore_is
    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpi_character, mpi_logical, mpi_integer
    use shr_string_mod,  only: shr_string_toUpper
    use string_utils,    only: int2str
    
    ! Dummy argument: filepath for file containing namelist input
    character(len=*), intent(in) :: nlfile

    ! Local variables
    integer                     :: unitn, ierr
    character(len=*), parameter :: subname = 'budget_readnl :: '

    namelist /thermo_budget_nl/  thermo_budget_history, thermo_budget_histfile_num
    !-----------------------------------------------------------------------

    if (masterproc) then
       open(newunit=unitn, file=trim(nlfile), status='old')
       call find_group_name(unitn, 'thermo_budget_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, thermo_budget_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname//'ERROR reading namelist, thermo_budget_nl, errcode = '//int2str(ierr))
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(thermo_budget_history         , 1  , mpi_logical  , masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: thermo_budget_history")
    call mpi_bcast(thermo_budget_histfile_num    , 1  , mpi_integer  , masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: thermo_budget_histfile_num")

    ! Write out thermo_budget options
    if (masterproc) then
       if (thermo_budget_history) then
          if (dycore_is('EUL').or.dycore_is('FV').or.dycore_is('FV3')) then
             call endrun(subname//'ERROR thermodynamic budgets not implemented for this dycore')
          else
             write(iulog,*)'Thermo budgets will be written to the log file and diagnostics saved to history file:',&
               thermo_budget_histfile_num
          end if
       end if
    end if
  end subroutine budget_readnl

  !==============================================================================================

  subroutine budget_init()
    use time_manager,         only:  get_step_size

    dstepsize=get_step_size()

  end subroutine budget_init

  !==============================================================================================

  subroutine e_m_snapshot (name, pkgtype, longname, cslam)
    use dycore,          only: dycore_is  

    character(len=*), intent(in)           :: &
         name      ! budget name used as variable name in history file output (8 char max)
    character(len=*), intent(in)           :: &
         pkgtype      ! budget type either phy or dyn
    character(len=*), intent(in)           :: &
         longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)
    logical,          intent(in), optional :: &
         cslam  ! true => CSLAM used to transport mass tracers

    character (cl)                         :: errmsg
    character (len=max_fieldname_len)      :: name_str
    character (cl)                         :: desc_str, units_str
    character (cl)                         :: gridname
    logical                                :: cslamtr        ! using cslam transport for mass tracers
    integer                                :: ivars
    character(len=*), parameter            :: sub='e_m_snapshot'
    !-----------------------------------------------------------------------

    if (thermo_budget_history) then
       if (present(cslam)) then
          cslamtr=cslam
       else
          cslamtr = .false.
       end if
       do ivars=1, thermo_budget_num_vars
          write(name_str,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(name))
          write(desc_str,*) TRIM(ADJUSTL(thermo_budget_vars_descriptor(ivars)))," ", &
               TRIM(ADJUSTL(longname))
          write(units_str,*) TRIM(ADJUSTL(thermo_budget_vars_unit(ivars)))

          budget_num = budget_num+1
          ! set budget name and constants
          budget_name(budget_num) = trim(name_str)
          budget_longname(budget_num) = trim(desc_str)

          budget_optype(budget_num)='stg'
          budget_pkgtype(budget_num)=pkgtype
          budget_stagename(budget_num)= trim(name)

          if (pkgtype=='phy') then
             gridname='physgrid'
          else
             if (dycore_is('SE')) then
                if (cslamtr .and. thermo_budget_vars_massv(ivars)) then
                   gridname='FVM'
                else
                   gridname='GLL'
                end if
             else if (dycore_is('MPAS')) then
                 gridname='mpas_cell'
             else
                write(errmsg, *) sub, ': budget_add is only supported for MPAS and SE dycores'
                call endrun(errmsg)
             end if
          end if
          call addfld (TRIM(ADJUSTL(name_str)),   horiz_only, 'N', TRIM(ADJUSTL(units_str)),TRIM(ADJUSTL(desc_str)),gridname=gridname)
          call add_default(TRIM(ADJUSTL(name_str)), thermo_budget_histfile_num, 'N')
       end do
    end if
  end subroutine e_m_snapshot

  !==============================================================================

  subroutine e_m_budget (name, stg1name, stg2name, pkgtype, optype, longname, cslam)
    use dycore,          only: dycore_is  


    ! Register a budget.

    character(len=*), intent(in) :: &
         name,stg1name,stg2name   ! budget name used as variable name in history file output (8 char max)

    character(len=*), intent(in) :: &
         pkgtype    ! budget type either phy or dyn

    character(len=*), intent(in) :: &
         optype    !  dif (difference) or sum

    character(len=*), intent(in) :: &
         longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)

    logical,          intent(in), optional :: &
         cslam       ! true => use cslam to transport mass variables

    character(len=*), parameter            :: sub='e_m_budget'
    character(cl)                          :: errmsg
    character(len=1)                       :: opchar
    character (len=max_fieldname_len)      :: name_str
    character (cl)                         :: desc_str, units_str
    character (cl)                         :: gridname
    character (cl)                         :: strstg1, strstg2
    integer                                :: ivars
    logical                                :: cslamtr        ! using cslam transport for mass tracers
    !-----------------------------------------------------------------------

    if (thermo_budget_history) then
       if (present(cslam)) then
          cslamtr=cslam
       else
          cslamtr = .false.
       end if

       ! register history budget variables
       do ivars=1, thermo_budget_num_vars
          write(name_str,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(name))
          write(strstg1,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(stg1name))
          write(strstg2,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(stg2name))
          write(desc_str,*) TRIM(ADJUSTL(thermo_budget_vars_descriptor(ivars)))," ", &
               TRIM(ADJUSTL(longname))
          write(units_str,*) TRIM(ADJUSTL(thermo_budget_vars_unit(ivars)))

          budget_num = budget_num + 1
          budget_pkgtype(budget_num)=pkgtype

          ! set budget name and constants
          budget_name(budget_num) = trim(name_str)
          budget_longname(budget_num) = trim(desc_str)

          if (optype=='dif') opchar='-'
          if (optype=='sum') opchar='+'
          if (optype=='stg') then
             write(errmsg,*) sub, ': FATAL: bad value optype should be sum of dif:', optype
             call endrun(errmsg)
          end if
          budget_stg1name(budget_num) = trim(adjustl(strstg1))
          budget_stg2name(budget_num) = trim(adjustl(strstg2))
          budget_stagename(budget_num)= trim(adjustl(strstg1))//trim(opchar)//trim(adjustl(strstg2))
          budget_optype(budget_num)=optype

          if (pkgtype=='phy') then
            gridname='physgrid'
          else
             if (dycore_is('SE')) then
                if (cslamtr .and. thermo_budget_vars_massv(ivars)) then
                   gridname='FVM'
                else
                   gridname='GLL'
                end if
             else if (dycore_is('MPAS')) then
                gridname='mpas_cell'
             else
                write(errmsg, *) sub, ': budget_add is only supported for MPAS and SE dycores'
                call endrun(errmsg)
             end if
          end if
          call addfld (TRIM(ADJUSTL(name_str)),   horiz_only, 'N', TRIM(ADJUSTL(units_str)),TRIM(ADJUSTL(desc_str)), &     
               gridname=gridname,op=optype,op_f1name=TRIM(ADJUSTL(strstg1)),op_f2name=TRIM(ADJUSTL(strstg2)))
          call add_default(TRIM(ADJUSTL(name_str)), thermo_budget_histfile_num, 'N')
       end do
    end if
  end subroutine e_m_budget

  !==============================================================================

  subroutine budget_get_global (name, me_idx, global)

    use cam_history,          only: get_field_properties
    use cam_history_support,  only: active_entry,ptapes
    use cam_thermo,           only: thermo_budget_vars_massv

    ! Get the global integral of a budget. Endrun will be called
    ! when name is not found.
    !-----------------------------Arguments---------------------------------
    character(len=*),  intent(in)  :: name    ! budget name
    integer,           intent(in)  :: me_idx  ! mass energy variable index
    real(r8),          intent(out) :: global  ! global integral of the budget field

    !---------------------------Local workspace-----------------------------
    type (active_entry), pointer   :: tape(:)                    ! history tapes
    character (len=max_fieldname_len) :: name_str
    character(cl)                  :: errmsg
    integer                        :: b_ind                      ! budget index
    integer                        :: h_ind(ptapes)              ! hentry index
    integer                        :: m_ind                      ! masterlist index
    integer                        :: idx,pidx,midx,uidx         ! substring index for sum dif char
    integer                        :: m                          ! budget index
    logical                        :: found                      ! true if global integral found

    character(len=*), parameter    :: sub='budget_get_global'
    !-----------------------------------------------------------------------
    ! Initialize tape pointer here to avoid initialization only on first invocation
    nullify(tape)
    
    name_str=''
    write(name_str,*) TRIM(ADJUSTL(name))

    midx=index(name_str, '-')
    pidx=index(name_str, '+')
    idx=midx+pidx

    ! check for budget using stagename short format (stg1//op//stg2) where stg1 is name without thermo string appended
    if (idx /= 0 .and. (midx==0 .or. pidx==0)) then
       write(name_str,*) TRIM(ADJUSTL(thermo_budget_vars(me_idx)))//"_"//trim(adjustl(name_str(1:idx)))// &
            TRIM(ADJUSTL(thermo_budget_vars(me_idx)))//"_"//TRIM(ADJUSTL(name_str(idx+1:)))
    end if

    uidx=index(name_str, '_')
    if (uidx == 0) then
       !This is a stage name need to append the type of thermo variable using input index
       write(name_str,*) TRIM(ADJUSTL(thermo_budget_vars(me_idx)))//"_"//trim(adjustl(name_str(1:)))
    end if

    b_ind=budget_ind_byname(trim(adjustl(name_str)))

    if (b_ind < 0) call endrun(sub//': FATAL field name '//name//' not found'//' looked for '//trim(adjustl(name_str)))

    write(name_str,*) TRIM(ADJUSTL(budget_name(b_ind)))

    ! Find budget name in list and return global value
    call get_field_properties(trim(adjustl(name_str)), found, tape_out=tape, ff_out=m_ind, f_out=h_ind)

    if (found.and.h_ind(thermo_budget_histfile_num)>0) then
       call tape(thermo_budget_histfile_num)%hlist(h_ind(thermo_budget_histfile_num))%get_global(global)
       if (.not. thermo_budget_vars_massv(me_idx)) &
            global=global/dstepsize
    else
       write(errmsg,*) sub, ': FATAL: name not found: ', trim(name)
       call endrun(errmsg)
    end if

  end subroutine budget_get_global
  !==============================================================================
  function budget_ind_byname (name)
    !
    ! Get the index of a budget.  Ret -1 for not found
    !-----------------------------Arguments---------------------------------
    character(len=*),  intent(in)  :: name  ! budget name

    !---------------------------Local workspace-----------------------------
    integer                        :: budget_ind_byname   ! function return
    integer                        :: m                   ! budget index
    !-----------------------------------------------------------------------
    ! Find budget name in list
    budget_ind_byname  = -1
    do m = 1, budget_array_max
       if (trim(adjustl(name)) == trim(adjustl(budget_name(m))).or.trim(adjustl(name)) == trim(adjustl(budget_stagename(m)))) then
          budget_ind_byname  = m
          return
       end if
    end do
  end function budget_ind_byname

  !==============================================================================

  function is_budget(name)

    ! Get the index of a budget.  

    !-----------------------------Arguments---------------------------------
    character(len=*),  intent(in)  :: name  ! budget name

    !---------------------------Local workspace-----------------------------
    logical                        :: is_budget           ! function return
    integer                        :: m                   ! budget index
    !-----------------------------------------------------------------------

    ! Find budget name in list of defined budgets

    is_budget = .false.
    do m = 1, budget_array_max
       if (trim(adjustl(name)) == trim(adjustl(budget_name(m))).or.trim(adjustl(name)) == trim(adjustl(budget_stagename(m)))) then
          is_budget = .true.
          return
       end if
    end do
  end function is_budget

  !===========================================================================

end module budgets
