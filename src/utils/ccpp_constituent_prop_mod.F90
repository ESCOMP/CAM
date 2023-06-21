! This module is a placeholder for the CCPP generated module of the same name
module ccpp_constituent_prop_mod

   implicit none
   private

   !Define stub version of constituent properties mod
   type, public :: ccpp_constituent_prop_ptr_t
      logical, private :: thermo_active = .false.
   contains
      procedure :: standard_name     => ccpt_get_standard_name
      procedure :: is_thermo_active  => ccpt_is_thermo_active
      procedure :: set_thermo_active => ccpt_set_thermo_active
   end type ccpp_constituent_prop_ptr_t

   !CCPP properties init routine
   public  :: ccpp_const_props_init

   !Public properties DDT variable:
   type(ccpp_constituent_prop_ptr_t), allocatable, public :: ccpp_const_props(:)

contains

!+++++++++++++++++++++++++++++++++++++++
!CCPP constituent properties DDT methods
!+++++++++++++++++++++++++++++++++++++++

   subroutine ccpt_get_standard_name(this, std_name, errcode, errmsg)
      ! Return this constituent's standard name

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(in)  :: this
      character(len=*),                     intent(out) :: std_name
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      std_name = 'Not Used!'

      !Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = 'Still Not Used!'
      end if

   end subroutine ccpt_get_standard_name

   !------

   subroutine ccpt_is_thermo_active(this, val_out, errcode, errmsg)

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t), intent(in)  :: this
      logical,                              intent(out) :: val_out
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      !Pass back thermo active property:
      val_out = this%thermo_active

      !Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = 'Still Not Used!'
      end if

   end subroutine ccpt_is_thermo_active

   !------

   subroutine ccpt_set_thermo_active(this, thermo_flag, errcode, errmsg)
      ! Set whether this constituent is thermodynamically active, which
      ! means that certain physics schemes will use this constitutent
      ! when calculating thermodynamic quantities (e.g. enthalpy).

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(inout) :: this
      logical,                              intent(in)    :: thermo_flag
      integer,          optional,           intent(out)   :: errcode
      character(len=*), optional,           intent(out)   :: errmsg

      !Set thermodynamically active flag for this constituent:
      this%thermo_active = thermo_flag

      !Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = 'Still Not Used!'
      end if

   end subroutine ccpt_set_thermo_active

!+++++++++++++++++++++++++++++++++++++++++++++
!CCPP constituents stub initialization routine
!+++++++++++++++++++++++++++++++++++++++++++++

subroutine ccpp_const_prop_init()

    !Use statements:
    use constituents,    only: pcnst
    use cam_abortutils,  only: handle_allocate_error
    use air_composition, only: thermodynamic_active_species_idx

    !Local variables:
    integer :: ierr
    integer :: m

    character(len=*), parameter :: subname = 'ccpp_const_prop_init:'

    !Allocate constituents object:
    allocate(ccpp_const_props(pcnst), stat=ierr)

    !Check if allocation succeeded:
    call handle_allocate_error(ierr, subname, 'ccpp_const_props(pcnst)')

    !Set "thermo_active" property:
    do m = 1,pcnst
       if (any(thermodynamic_active_species_idx == m)) then
          call ccpp_const_props(m)%set_thermo_active(.true.)
       else
          call ccpp_const_props(m)%set_thermo_active(.false.)
       end if
    end do

end subroutine ccpp_const_prop_init

end module ccpp_constituent_prop_mod
