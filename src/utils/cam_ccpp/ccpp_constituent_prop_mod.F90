! This module is the CAM version of the CCPP generated module of the same name
module ccpp_constituent_prop_mod

   implicit none
   private

   ! Define CAM version of constituent properties mod
   type, public :: ccpp_constituent_prop_ptr_t
      logical, private :: thermo_active = .false.
      logical, private :: water_species = .false.
      logical, private :: species_is_dry
      logical, private :: species_is_wet
      character(len=256) :: std_name = ''

      contains
      procedure :: standard_name     => ccp_get_standard_name
      procedure :: set_standard_name => ccp_set_standard_name
      procedure :: is_thermo_active  => ccp_is_thermo_active
      procedure :: is_water_species  => ccp_is_water_species
      procedure :: set_thermo_active => ccp_set_thermo_active
      procedure :: set_water_species => ccp_set_water_species
      procedure :: is_dry            => ccp_is_dry
      procedure :: set_dry           => ccp_set_dry
      procedure :: is_wet            => ccp_is_wet
      procedure :: set_wet           => ccp_set_wet

   end type ccpp_constituent_prop_ptr_t

   type, public :: ccpp_constituent_properties_t
      contains
      procedure :: instantiate
   end type ccpp_constituent_properties_t

   ! CCPP properties init routine
   public :: ccpp_const_props_init

   ! Public properties DDT variable:
   type(ccpp_constituent_prop_ptr_t), allocatable, public :: ccpp_const_props(:)
   integer, public, parameter :: int_unassigned = HUGE(-1)

contains

!+++++++++++++++++++++++++++++++++++++++
!CCPP constituent properties DDT methods
!+++++++++++++++++++++++++++++++++++++++

   subroutine ccp_get_standard_name(this, std_name, errcode, errmsg)
      ! Return this constituent's standard name

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(in)  :: this
      character(len=*),                     intent(out) :: std_name
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      std_name = this%std_name

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_get_standard_name

   !------

   subroutine ccp_set_standard_name(this, std_name, errcode, errmsg)
      ! Set this constituent's standard name

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(inout)  :: this
      character(len=*),                     intent(in)  :: std_name
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      this%std_name = std_name

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_set_standard_name

   !------

   subroutine ccp_is_thermo_active(this, val_out, errcode, errmsg)

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(in)  :: this
      logical,                              intent(out) :: val_out
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      ! Pass back thermo active property:
      val_out = this%thermo_active

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_is_thermo_active

   !------

   subroutine ccp_is_water_species(this, val_out, errcode, errmsg)

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(in)  :: this
      logical,                              intent(out) :: val_out
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      ! Pass back water species property:
      val_out = this%water_species

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_is_water_species

   !------

   subroutine ccp_is_dry(this, val_out, errcode, errmsg)

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(in)  :: this
      logical,                              intent(out) :: val_out
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      ! Pass back water species property:
      val_out = this%species_is_dry

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_is_dry

   !------

   subroutine ccp_set_thermo_active(this, thermo_flag, errcode, errmsg)
      ! Set whether this constituent is thermodynamically active, which
      ! means that certain physics schemes will use this constitutent
      ! when calculating thermodynamic quantities (e.g. enthalpy).

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(inout) :: this
      logical,                              intent(in)    :: thermo_flag
      integer,          optional,           intent(out)   :: errcode
      character(len=*), optional,           intent(out)   :: errmsg

      ! Set thermodynamically active flag for this constituent:
      this%thermo_active = thermo_flag

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_set_thermo_active

   !------

   subroutine ccp_set_water_species(this, water_flag, errcode, errmsg)
      ! Set whether this constituent is a water species, which means
      ! that this constituent represents a particular phase or type
      ! of water in the atmosphere.

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(inout) :: this
      logical,                              intent(in)    :: water_flag
      integer,          optional,           intent(out)   :: errcode
      character(len=*), optional,           intent(out)   :: errmsg

      ! Set thermodynamically active flag for this constituent:
      this%water_species = water_flag

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_set_water_species

   subroutine ccp_set_dry(this, dry_flag, errcode, errmsg)
      ! Set whether this constituent is a dry species or not using the dry_flag which is passed in

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(inout) :: this
      logical,                              intent(in)    :: dry_flag
      integer,          optional,           intent(out)   :: errcode
      character(len=*), optional,           intent(out)   :: errmsg

      ! Set dry_flag for this constituent:
      this%species_is_dry = dry_flag

      ! Provide err values if requested:
      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_set_dry

   subroutine ccp_is_wet(this, val_out, errcode, errmsg)

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(in)  :: this
      logical,                              intent(out) :: val_out
      integer,          optional,           intent(out) :: errcode
      character(len=*), optional,           intent(out) :: errmsg

      val_out = this%species_is_wet

      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_is_wet

   !------

   subroutine ccp_set_wet(this, wet_flag, errcode, errmsg)

      ! Dummy arguments
      class(ccpp_constituent_prop_ptr_t),   intent(inout) :: this
      logical,                              intent(in)    :: wet_flag
      integer,          optional,           intent(out)   :: errcode
      character(len=*), optional,           intent(out)   :: errmsg

      this%species_is_wet = wet_flag

      if(present(errcode)) then
         errcode = 0
      end if
      if(present(errmsg)) then
         errmsg = ''
      end if

   end subroutine ccp_set_wet

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!CAM-equivalent CCPP constituents initialization routine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine ccpp_const_props_init(ix_qv)

    ! Use statements:
    use constituents,    only: pcnst, cnst_name, cnst_get_ind, cnst_get_type_byind
    use cam_abortutils,  only: handle_allocate_error
    use air_composition, only: dry_air_species_num
    use air_composition, only: thermodynamic_active_species_idx

    integer, intent(in) :: ix_qv
    ! Local variables:
    integer :: ierr
    integer :: m
    integer :: ix_cldliq, ix_cldice, ix_numliq, ix_numice

    character(len=*), parameter :: subname = 'ccpp_const_prop_init:'

    ! Allocate constituents object:
    allocate(ccpp_const_props(pcnst), stat=ierr)

    ! Check if allocation succeeded:
    call handle_allocate_error(ierr, subname, 'ccpp_const_props(pcnst)')

    ! Set "thermo_active" property:
    do m = 1,pcnst
       if(any(thermodynamic_active_species_idx == m)) then
          call ccpp_const_props(m)%set_thermo_active(.true.)
       end if
    end do

    ! Set "water_species" property:
    do m=1,pcnst
       if(any(thermodynamic_active_species_idx(dry_air_species_num+1:) == m)) then
          call ccpp_const_props(m)%set_water_species(.true.)
       end if
    end do

    ! Set dry/wet properties:
    do m=1,pcnst
       if (cnst_get_type_byind(m).eq.'dry') then
          call ccpp_const_props(m)%set_dry(.true.)
          call ccpp_const_props(m)%set_wet(.false.)
       else
          call ccpp_const_props(m)%set_dry(.false.)
          call ccpp_const_props(m)%set_wet(.true.)
       end if
    end do

    ! Set "std_name" for all constituents using cnst_name as baseline.
    ! This allows chemistry species (O3, NO, NO2, ...) to be looked up
    ! via ccpp_constituent_index using their constituent name directly.
    do m=1,pcnst
       call ccpp_const_props(m)%set_standard_name(trim(cnst_name(m)))
    end do

    ! Override with proper CCPP standard names for known species:
    call ccpp_const_props(ix_qv)%set_standard_name( &
         'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water')

    call cnst_get_ind('CLDLIQ', ix_cldliq, abort=.false.)
    if (ix_cldliq > 0) then
       call ccpp_const_props(ix_cldliq)%set_standard_name( &
            'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water')
    end if

    call cnst_get_ind('CLDICE', ix_cldice, abort=.false.)
    if (ix_cldice > 0) then
       call ccpp_const_props(ix_cldice)%set_standard_name( &
            'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water')
    end if

    call cnst_get_ind('NUMLIQ', ix_numliq, abort=.false.)
    if (ix_numliq > 0) then
       call ccpp_const_props(ix_numliq)%set_standard_name( &
            'mass_number_concentration_of_cloud_liquid_wrt_moist_air_and_condensed_water')
    end if

    call cnst_get_ind('NUMICE', ix_numice, abort=.false.)
    if (ix_numice > 0) then
       call ccpp_const_props(ix_numice)%set_standard_name( &
            'mass_number_concentration_of_ice_wrt_moist_air_and_condensed_water')
    end if

end subroutine ccpp_const_props_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!CAM-equivalent stub so dynamic constituents register routines can build
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine instantiate(this, std_name, long_name, diag_name, units, vertical_dim,  &
        advected, default_value, min_value, molar_mass, water_species,         &
        mixing_ratio_type, errcode, errmsg)
      use ccpp_kinds, only: kind_phys

      ! Dummy arguments
      class(ccpp_constituent_properties_t), intent(inout) :: this
      character(len=*),                     intent(in)    :: std_name
      character(len=*),                     intent(in)    :: long_name
      character(len=*),                     intent(in)    :: diag_name
      character(len=*),                     intent(in)    :: units
      character(len=*),                     intent(in)    :: vertical_dim
      logical, optional,                    intent(in)    :: advected
      real(kind_phys), optional,            intent(in)    :: default_value
      real(kind_phys), optional,            intent(in)    :: min_value
      real(kind_phys), optional,            intent(in)    :: molar_mass
      logical, optional,                    intent(in)    :: water_species
      character(len=*), optional,           intent(in)    :: mixing_ratio_type
      integer,                              intent(out)   :: errcode
      character(len=*),                     intent(out)   :: errmsg

      ! STUB DOES NOTHING

end subroutine instantiate

end module ccpp_constituent_prop_mod
