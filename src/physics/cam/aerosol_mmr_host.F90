module aerosol_mmr_host
! Host-model specific module for handling aerosol concentrations.
! This is the CAM implementation.
!
! Each host model provides a module of the same name implementing:
! 1) an aero_host_binding_t type containing host-model specific references
!    e.g., state and pbuf in CAM; CCPP constituents array in SIMA.
!    This allows the aerosol_state to store and pass this handle without
!    referencing host-model data structures directly, making them portable.
! 2) rad_cnst_get_aer_mmr, rad_cnst_get_mode_num, bin_num, bin_mmr
!    subroutines to retrieve data given the above host handle.

use shr_kind_mod,   only: r8 => shr_kind_r8
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc

implicit none
private

! This host-binding type contains host-model specific data references
! (e.g., physics_state and physics_buffer_desc for CAM)
! needed to retrieve aerosol fields.
!
! This type is referenced in the shared abstract aerosol interface code,
! which can then be passed around without referencing host-model types.
! Only this module needs to be aware of what is inside this type.
type :: aero_host_binding_t
   type(physics_state),       pointer :: state   => null()
   type(physics_buffer_desc), pointer :: pbuf(:) => null()
end type aero_host_binding_t

! define generic interface for MMR retrieval
interface rad_cnst_get_aer_mmr
   module procedure rad_cnst_get_aer_mmr_by_idx
   module procedure rad_cnst_get_mam_mmr_by_idx
   module procedure rad_cnst_get_aer_mmr_by_idx_host
   module procedure rad_cnst_get_mam_mmr_by_idx_host
end interface

! generic interfaces dispatching between the legacy (state, pbuf) variants
! and the host-binding handle variants
interface rad_cnst_get_mode_num
   module procedure rad_cnst_get_mode_num_cam
   module procedure rad_cnst_get_mode_num_host
end interface

interface rad_cnst_get_bin_num
   module procedure rad_cnst_get_bin_num_cam
   module procedure rad_cnst_get_bin_num_host
end interface

interface rad_cnst_get_bin_mmr_by_idx
   module procedure rad_cnst_get_bin_mmr_by_idx_cam
   module procedure rad_cnst_get_bin_mmr_by_idx_host
end interface

! values for constituents with requested value of zero
real(r8), allocatable, target :: zero_cols(:,:)

public :: aero_host_binding_t
public :: aero_host_binding   ! build a handle from host data structures
public :: aerosol_mmr_init    ! allocate zero_cols
public :: get_host_idx
public :: resolve_mode_idx
public :: resolve_bin_idx
public :: resolve_bulk_idx
public :: rad_cnst_get_aer_mmr
public :: rad_cnst_get_mam_mmr_idx
public :: rad_cnst_get_mode_num
public :: rad_cnst_get_mode_num_idx
public :: rad_cnst_get_bin_mmr_by_idx
public :: rad_cnst_get_bin_num
public :: rad_cnst_get_carma_mmr_idx
public :: rad_cnst_get_bin_mmr
public :: rad_aer_diag_init
public :: rad_aer_diag_out

!==============================================================================
contains
!==============================================================================

subroutine aerosol_mmr_init()
   use ppgrid, only: pcols, pver
   ! Allocate zero_cols array (must be called after ppgrid is set up)
   if (.not. allocated(zero_cols)) then
      allocate(zero_cols(pcols,pver))
      zero_cols = 0._r8
   end if
end subroutine aerosol_mmr_init

!================================================================================================

function aero_host_binding(state, pbuf) result(host)

   ! Build a host-binding handle from CAM host data structures.
   ! Called from host-side wiring only (aerosol_instances_mod); the
   ! resulting handle is stored opaquely by the aerosol_state objects.

   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   type(aero_host_binding_t) :: host

   host%state => state
   host%pbuf  => pbuf

end function aero_host_binding

!================================================================================================

integer function get_host_idx(source, name, routine)

   ! get index of name in internal CAM array; either the constituent array
   ! or the physics buffer

   use physics_buffer, only: pbuf_get_index
   use constituents,   only: cnst_get_ind
   use cam_abortutils, only: endrun

   character(len=*), intent(in) :: source
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: routine  ! name of calling routine

   integer :: idx
   integer :: errcode
   !-----------------------------------------------------------------------------

   if (source(1:1) == 'N') then

      idx = pbuf_get_index(trim(name),errcode)
      if (errcode < 0) then
         call endrun(routine//' ERROR: cannot find physics buffer field '//trim(name))
      end if

   else if (source(1:1) == 'A') then

      call cnst_get_ind(trim(name), idx, abort=.false.)
      if (idx < 0) then
         call endrun(routine//' ERROR: cannot find constituent field '//trim(name))
      end if

   else if (source(1:1) == 'Z') then

      idx = -1

   else

      call endrun(routine//' ERROR: invalid source for specie '//trim(name))

   end if

   get_host_idx = idx

end function get_host_idx

!===========================

subroutine resolve_mode_idx(modes)

   ! Initialize the mode definitions by looking up the relevent indices in the
   ! constituent and pbuf arrays, and getting the physprop IDs

   use phys_prop,      only: physprop_get_id
   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: modes_t

   ! Arguments
   type(modes_t), intent(inout) :: modes

   ! Local variables
   integer :: m, ispec, nspec

   character(len=*), parameter :: routine = 'resolve_mode_idx'
   !-----------------------------------------------------------------------------

   do m = 1, modes%nmodes

      ! indices for number mixing ratio components
      modes%comps(m)%idx_num_a = get_host_idx(modes%comps(m)%source_num_a, modes%comps(m)%camname_num_a, routine)
      modes%comps(m)%idx_num_c = get_host_idx(modes%comps(m)%source_num_c, modes%comps(m)%camname_num_c, routine)

      ! allocate memory for species
      nspec = modes%comps(m)%nspec
      allocate( &
         modes%comps(m)%idx_mmr_a(nspec), &
         modes%comps(m)%idx_mmr_c(nspec), &
         modes%comps(m)%idx_props(nspec)  )

      do ispec = 1, nspec

         ! indices for species mixing ratio components
         modes%comps(m)%idx_mmr_a(ispec) = get_host_idx(modes%comps(m)%source_mmr_a(ispec), &
                                                   modes%comps(m)%camname_mmr_a(ispec), routine)
         modes%comps(m)%idx_mmr_c(ispec) = get_host_idx(modes%comps(m)%source_mmr_c(ispec), &
                                                   modes%comps(m)%camname_mmr_c(ispec), routine)

         ! get physprop ID
         modes%comps(m)%idx_props(ispec) = physprop_get_id(modes%comps(m)%props(ispec))
         if (modes%comps(m)%idx_props(ispec) == -1) then
            call endrun(routine//' : ERROR idx not found for '//trim(modes%comps(m)%props(ispec)))
         end if

      end do

   end do

end subroutine resolve_mode_idx

!===========================

subroutine resolve_bin_idx(bins)

   ! Initialize the bin definitions by looking up the relevent indices in the
   ! constituent and pbuf arrays, and getting the physprop IDs

   use phys_prop,      only: physprop_get_id
   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: bins_t

   ! Arguments
   type(bins_t), intent(inout) :: bins

   ! Local variables
   integer :: m, ispec, nspec

   character(len=*), parameter :: routine = 'resolve_bin_idx'
   !-----------------------------------------------------------------------------

   do m = 1, bins%nbins

      ! indices for number mixing ratio components
      bins%comps(m)%idx_num_a = get_host_idx(bins%comps(m)%source_num_a, bins%comps(m)%camname_num_a, routine)
      bins%comps(m)%idx_num_c = get_host_idx(bins%comps(m)%source_num_c, bins%comps(m)%camname_num_c, routine)
      if ( bins%comps(m)%source_mass_a /= 'NOTSET' .and. bins%comps(m)%camname_mass_a /= 'NOTSET' ) then
         bins%comps(m)%idx_mass_a = get_host_idx(bins%comps(m)%source_mass_a, bins%comps(m)%camname_mass_a, routine)
      endif
      if ( bins%comps(m)%source_mass_c /= 'NOTSET' .and. bins%comps(m)%camname_mass_c /= 'NOTSET' ) then
         bins%comps(m)%idx_mass_c = get_host_idx(bins%comps(m)%source_mass_c, bins%comps(m)%camname_mass_c, routine)
      endif

      ! allocate memory for species
      nspec = bins%comps(m)%nspec
      allocate( &
         bins%comps(m)%idx_mmr_a(nspec), &
         bins%comps(m)%idx_mmr_c(nspec), &
         bins%comps(m)%idx_props(nspec)  )

      do ispec = 1, nspec

         ! indices for species mixing ratio components
         bins%comps(m)%idx_mmr_a(ispec) = get_host_idx(bins%comps(m)%source_mmr_a(ispec), &
                                                   bins%comps(m)%camname_mmr_a(ispec), routine)
         bins%comps(m)%idx_mmr_c(ispec) = get_host_idx(bins%comps(m)%source_mmr_c(ispec), &
                                                   bins%comps(m)%camname_mmr_c(ispec), routine)

         ! get physprop ID
         bins%comps(m)%idx_props(ispec) = physprop_get_id(bins%comps(m)%props(ispec))
         if (bins%comps(m)%idx_props(ispec) == -1) then
            call endrun(routine//' : ERROR idx not found for '//trim(bins%comps(m)%props(ispec)))
         end if

      end do

   end do

end subroutine resolve_bin_idx

!===========================

subroutine resolve_bulk_idx(aerlist)

   ! Resolve host-specific indices for bulk aerosols.
   ! Must be called before list_resolve_physprops (which resolves physprop IDs).

   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: aerlist_t

   type(aerlist_t), intent(inout) :: aerlist

   integer :: i
   character(len=*), parameter :: routine = 'resolve_bulk_idx'
   !-----------------------------------------------------------------------------

   do i = 1, aerlist%numaerosols
      aerlist%aer(i)%idx = get_host_idx(aerlist%aer(i)%source, aerlist%aer(i)%camname, routine)
   end do

end subroutine resolve_bulk_idx

!================================================================================================

subroutine rad_cnst_get_aer_mmr_by_idx(list_idx, aer_idx, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the aerosol from the specified
   ! climate or diagnostic list.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use radiative_aerosol_definitions, only: N_DIAG, aerlist_t, bulk_aerosol_list

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: aer_idx
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: idx
   character(len=1) :: source
   type(aerlist_t), pointer :: aerlist
   character(len=*), parameter :: subname = 'rad_cnst_get_aer_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => bulk_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid input aerosol index
   if (aer_idx < 1  .or.  aer_idx > aerlist%numaerosols) then
      write(iulog,*) subname//': aer_idx= ', aer_idx, '  numaerosols= ', aerlist%numaerosols
      call endrun(subname//': aerosol list index out of range')
   end if

   ! Get data source
   source = aerlist%aer(aer_idx)%source
   idx    = aerlist%aer(aer_idx)%idx
   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_aer_mmr_by_idx

!================================================================================================

subroutine rad_cnst_get_aer_mmr_by_idx_host(list_idx, aer_idx, host, mmr)

   ! Host-binding handle variant: unpack the handle and delegate.

   integer,                   intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                   intent(in) :: aer_idx
   type(aero_host_binding_t), intent(in) :: host
   real(r8),                  pointer    :: mmr(:,:)

   type(physics_buffer_desc), pointer :: pbuf(:)

   pbuf => host%pbuf
   call rad_cnst_get_aer_mmr_by_idx(list_idx, aer_idx, host%state, pbuf, mmr)

end subroutine rad_cnst_get_aer_mmr_by_idx_host

!================================================================================================

subroutine rad_cnst_get_mam_mmr_by_idx(list_idx, mode_idx, spec_idx, phase, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the modal aerosol specie from the specified
   ! climate or diagnostic list.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use radiative_aerosol_definitions, only: N_DIAG, modelist_t, modal_aerosol_list, modes

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: mode_idx    ! mode index
   integer,                     intent(in) :: spec_idx    ! index of specie in the mode
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => modal_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > modes%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', modes%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Get data source
   if (phase == 'a') then
      source = modes%comps(m_idx)%source_mmr_a(spec_idx)
      idx    = modes%comps(m_idx)%idx_mmr_a(spec_idx)
   else if (phase == 'c') then
      source = modes%comps(m_idx)%source_mmr_c(spec_idx)
      idx    = modes%comps(m_idx)%idx_mmr_c(spec_idx)
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_mam_mmr_by_idx

!================================================================================================

subroutine rad_cnst_get_mam_mmr_by_idx_host(list_idx, mode_idx, spec_idx, phase, host, mmr)

   ! Host-binding handle variant: unpack the handle and delegate.

   integer,                   intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                   intent(in) :: mode_idx    ! mode index
   integer,                   intent(in) :: spec_idx    ! index of specie in the mode
   character(len=1),          intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(aero_host_binding_t), intent(in) :: host
   real(r8),                  pointer    :: mmr(:,:)

   type(physics_buffer_desc), pointer :: pbuf(:)

   pbuf => host%pbuf
   call rad_cnst_get_mam_mmr_by_idx(list_idx, mode_idx, spec_idx, phase, host%state, pbuf, mmr)

end subroutine rad_cnst_get_mam_mmr_by_idx_host

!================================================================================================

subroutine rad_cnst_get_bin_mmr_by_idx_cam(list_idx, bin_idx, spec_idx, phase, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the modal aerosol specie from the specified
   ! climate or diagnostic list.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use radiative_aerosol_definitions, only: N_DIAG, binlist_t, sectional_aerosol_list, bins

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: bin_idx    ! mode index
   integer,                     intent(in) :: spec_idx    ! index of specie in the mode
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: s_idx
   integer :: idx
   integer :: lchnk
   character(len=1) :: source
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_mmr_by_idx'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sectional_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   s_idx = slist%idx(bin_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > bins%comps(s_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', bins%comps(s_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Get data source
   if (phase == 'a') then
      source = bins%comps(s_idx)%source_mmr_a(spec_idx)
      idx    = bins%comps(s_idx)%idx_mmr_a(spec_idx)
   else if (phase == 'c') then
      source = bins%comps(s_idx)%source_mmr_c(spec_idx)
      idx    = bins%comps(s_idx)%idx_mmr_c(spec_idx)
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   lchnk = state%lchnk

   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_bin_mmr_by_idx_cam

!================================================================================================

subroutine rad_cnst_get_bin_mmr_by_idx_host(list_idx, bin_idx, spec_idx, phase, host, mmr)

   ! Host-binding handle variant: unpack the handle and delegate.

   integer,                   intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                   intent(in) :: bin_idx     ! bin index
   integer,                   intent(in) :: spec_idx    ! index of specie in the bin
   character(len=1),          intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(aero_host_binding_t), intent(in) :: host
   real(r8),                  pointer    :: mmr(:,:)

   type(physics_buffer_desc), pointer :: pbuf(:)

   pbuf => host%pbuf
   call rad_cnst_get_bin_mmr_by_idx_cam(list_idx, bin_idx, spec_idx, phase, host%state, pbuf, mmr)

end subroutine rad_cnst_get_bin_mmr_by_idx_host

!================================================================================================

subroutine rad_cnst_get_mam_mmr_idx(mode_idx, spec_idx, idx)

   ! Return constituent index of mam specie mass mixing ratio for aerosol modes in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: modelist_t, modes, modal_aerosol_list

   ! Arguments
   integer, intent(in)  :: mode_idx    ! mode index
   integer, intent(in)  :: spec_idx    ! index of specie in the mode
   integer, intent(out) :: idx         ! index of specie in the constituent array

   ! Local variables
   integer :: m_idx
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mam_mmr_idx'
   !-----------------------------------------------------------------------------

   ! assume climate list (i.e., species are in the constituent array)
   mlist => modal_aerosol_list(0)

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > modes%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', modes%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Assume data source is interstitial since that's what's in the constituent array
   idx    = modes%comps(m_idx)%idx_mmr_a(spec_idx)

end subroutine rad_cnst_get_mam_mmr_idx

!================================================================================================

subroutine rad_cnst_get_carma_mmr_idx(bin_idx, spec_idx, idx)

   ! Return constituent index of camra species mass mixing ratio for aerosol bins in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: binlist_t, bins, sectional_aerosol_list

   ! Arguments
   integer, intent(in)  :: bin_idx     ! bin index
   integer, intent(in)  :: spec_idx    ! index of specie in the bin
   integer, intent(out) :: idx         ! index of specie in the constituent array

   ! Local variables
   integer :: b_idx
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_carma_mmr_idx'
   !-----------------------------------------------------------------------------

   ! assume climate list (i.e., species are in the constituent array)
   slist => sectional_aerosol_list(0)

   ! Check for valid bin index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding bin in the bin definition object
   b_idx = slist%idx(bin_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > bins%comps(b_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', bins%comps(b_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   ! Assume data source is interstitial since that's what's in the constituent array
   idx = bins%comps(b_idx)%idx_mmr_a(spec_idx)

end subroutine rad_cnst_get_carma_mmr_idx

!================================================================================================

subroutine rad_cnst_get_bin_mmr(list_idx, bin_idx, phase, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the aerosol bin from the specified
   ! climate or diagnostic list.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use radiative_aerosol_definitions, only: N_DIAG, binlist_t, sectional_aerosol_list, bins

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: bin_idx     ! bin index
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   character(len=1) :: source
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_mmr'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sectional_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid bin index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding bin in the bin definition object
   m_idx = slist%idx(bin_idx)

   ! Get data source
   if (phase == 'a') then
      source = bins%comps(m_idx)%source_mass_a
      idx    = bins%comps(m_idx)%idx_mass_a
   else if (phase == 'c') then
      source = bins%comps(m_idx)%source_mass_c
      idx    = bins%comps(m_idx)%idx_mass_c
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_bin_mmr

!================================================================================================

subroutine rad_cnst_get_mode_num_cam(list_idx, mode_idx, phase, state, pbuf, num)

   ! Return pointer to number mixing ratio for the aerosol mode from the specified
   ! climate or diagnostic list.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use radiative_aerosol_definitions, only: N_DIAG, modelist_t, modal_aerosol_list, modes

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: mode_idx    ! mode index
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: num(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_num'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => modal_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Get data source
   if (phase == 'a') then
      source = modes%comps(m_idx)%source_num_a
      idx    = modes%comps(m_idx)%idx_num_a
   else if (phase == 'c') then
      source = modes%comps(m_idx)%source_num_c
      idx    = modes%comps(m_idx)%idx_num_c
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   select case( source )
   case ('A')
      num => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, num)
   case ('Z')
      num => zero_cols
   end select

end subroutine rad_cnst_get_mode_num_cam

!================================================================================================

subroutine rad_cnst_get_mode_num_host(list_idx, mode_idx, phase, host, num)

   ! Host-binding handle variant: unpack the handle and delegate.

   integer,                   intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                   intent(in) :: mode_idx    ! mode index
   character(len=1),          intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(aero_host_binding_t), intent(in) :: host
   real(r8),                  pointer    :: num(:,:)

   type(physics_buffer_desc), pointer :: pbuf(:)

   pbuf => host%pbuf
   call rad_cnst_get_mode_num_cam(list_idx, mode_idx, phase, host%state, pbuf, num)

end subroutine rad_cnst_get_mode_num_host

!================================================================================================

subroutine rad_cnst_get_bin_num_cam(list_idx, bin_idx, phase, state, pbuf, num)

   ! Return pointer to number mixing ratio for the aerosol bin from the specified
   ! climate or diagnostic list.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use radiative_aerosol_definitions, only: N_DIAG, binlist_t, sectional_aerosol_list, bins

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in) :: bin_idx     ! bin index
   character(len=1),            intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: num(:,:)

   ! Local variables
   integer :: m_idx
   integer :: idx
   character(len=1) :: source
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_cnst_get_bin_num'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sectional_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Check for valid bin index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding bin in the bin definition object
   m_idx = slist%idx(bin_idx)

   ! Get data source
   if (phase == 'a') then
      source = bins%comps(m_idx)%source_num_a
      idx    = bins%comps(m_idx)%idx_num_a
   else if (phase == 'c') then
      source = bins%comps(m_idx)%source_num_c
      idx    = bins%comps(m_idx)%idx_num_c
   else
      write(iulog,*) subname//': phase= ', phase
      call endrun(subname//': unrecognized phase; must be "a" or "c"')
   end if

   select case( source )
   case ('A')
      num => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, num)
   case ('Z')
      num => zero_cols
   end select

end subroutine rad_cnst_get_bin_num_cam

!================================================================================================

subroutine rad_cnst_get_bin_num_host(list_idx, bin_idx, phase, host, num)

   ! Host-binding handle variant: unpack the handle and delegate.

   integer,                   intent(in) :: list_idx    ! index of the climate or a diagnostic list
   integer,                   intent(in) :: bin_idx     ! bin index
   character(len=1),          intent(in) :: phase       ! 'a' for interstitial, 'c' for cloud borne
   type(aero_host_binding_t), intent(in) :: host
   real(r8),                  pointer    :: num(:,:)

   type(physics_buffer_desc), pointer :: pbuf(:)

   pbuf => host%pbuf
   call rad_cnst_get_bin_num_cam(list_idx, bin_idx, phase, host%state, pbuf, num)

end subroutine rad_cnst_get_bin_num_host

!================================================================================================

subroutine rad_cnst_get_mode_num_idx(mode_idx, cnst_idx)

   ! Return constituent index of mode number mixing ratio for the aerosol mode in
   ! the climate list.

   ! This is a special routine to allow direct access to information in the
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.

   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: modelist_t, modes, modal_aerosol_list

   ! Arguments
   integer,  intent(in)  :: mode_idx    ! mode index
   integer,  intent(out) :: cnst_idx    ! constituent index

   ! Local variables
   integer :: m_idx
   character(len=1) :: source
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_cnst_get_mode_num'
   !-----------------------------------------------------------------------------

   ! assume climate list
   mlist => modal_aerosol_list(0)

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = mlist%idx(mode_idx)

   ! Check that source is 'A' which means the index is for the constituent array
   source = modes%comps(m_idx)%source_num_a
   if (source /= 'A') then
      write(iulog,*) subname//': source= ', source
      call endrun(subname//': requested mode number index not in constituent array')
   end if

   ! Return index in constituent array
   cnst_idx = modes%comps(m_idx)%idx_num_a

end subroutine rad_cnst_get_mode_num_idx
!================================================================================================

subroutine rad_aer_diag_init(alist)

! Add diagnostic fields to the master fieldlist.

   use cam_history,    only: addfld, fieldname_len, horiz_only
   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: aerlist_t

   type(aerlist_t), intent(inout) :: alist

   integer :: i, naer
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   naer = alist%numaerosols
   if (naer == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = alist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, naer

      ! construct names for mass per layer diagnostic fields
      name = 'm_' // trim(alist%aer(i)%camname) // trim(suffix)
      alist%aer(i)%mass_name = name
      long_name = trim(alist%aer(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), (/ 'lev' /), 'A', 'kg/m^2', trim(long_name))

      ! construct names for column burden diagnostic fields
      name = 'cb_' // trim(alist%aer(i)%camname) // trim(suffix)
      long_name = trim(alist%aer(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), horiz_only, 'A', 'kg/m^2', trim(long_name))

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_aer_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_aer_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_aer_diag_init

!================================================================================================

subroutine rad_aer_diag_out(list_idx, state, pbuf)

   ! Output the mass per layer, and total column burdens for aerosol
   ! constituents in either the climate or diagnostic lists.

   use ppgrid,         only: pcols, pver
   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use physconst,      only: rga
   use cam_history,    only: outfld
   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use radiative_aerosol_definitions, only: N_DIAG, aerlist_t, bulk_aerosol_list

   ! Arguments
   integer,                     intent(in) :: list_idx
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)

   ! Local variables
   integer :: i, naer, lchnk, ncol
   integer :: idx
   character(len=1)  :: source
   character(len=32) :: name, cbname
   real(r8)          :: mass(pcols,pver)
   real(r8)          :: cb(pcols)
   real(r8), pointer :: mmr(:,:)
   type(aerlist_t), pointer :: aerlist
   character(len=*), parameter :: subname = 'rad_aer_diag_out'
   !-----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointer with requested aerosol list
   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => bulk_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   naer = aerlist%numaerosols
   do i = 1, naer

      source = aerlist%aer(i)%source
      idx    = aerlist%aer(i)%idx
      name   = aerlist%aer(i)%mass_name
      ! construct name for column burden field by replacing the 'm_' prefix by 'cb_'
      cbname = 'cb_' // name(3:len_trim(name))

      select case( source )
      case ('A')
         mmr => state%q(:,:,idx)
      case ('N')
         call pbuf_get_field(pbuf, idx, mmr)
      end select

      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga
      call outfld(trim(name), mass, pcols, lchnk)

      cb(:ncol) = sum(mass(:ncol,:),2)
      call outfld(trim(cbname), cb, pcols, lchnk)

   end do

end subroutine rad_aer_diag_out

!================================================================================================

end module aerosol_mmr_host
