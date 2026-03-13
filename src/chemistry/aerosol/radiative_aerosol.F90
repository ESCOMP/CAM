module radiative_aerosol

!------------------------------------------------------------------------------------------------
!
! Facade module for aerosol definitions and queries.
!
! Provides query routines (rad_aer_get_info*, rad_aer_get_props*, etc.) and
! property-access routines that wrap phys_prop lookups.
! Init is via rad_aer_init (uses host-specific module aerosol_mmr_cam).
!
!------------------------------------------------------------------------------------------------

implicit none
private

! Generic interface for aerosol info queries.
interface rad_aer_get_info
   module procedure rad_aer_get_info
   module procedure rad_aer_get_info_by_mode
   module procedure rad_aer_get_info_by_mode_spec
   module procedure rad_aer_get_info_by_spectype
end interface

interface rad_aer_get_props
   module procedure rad_aer_get_props_by_idx
   module procedure rad_aer_get_mam_props_by_idx
end interface

! Public routines — aerosol queries (rad_aer_* naming)
public :: rad_aer_get_info
public :: rad_aer_get_info_by_mode, rad_aer_get_info_by_mode_spec
public :: rad_aer_get_info_by_spectype
public :: rad_aer_get_info_by_bin, rad_aer_get_info_by_bin_spec
public :: rad_aer_get_mode_idx, rad_aer_get_spec_idx
public :: rad_aer_num_name
public :: rad_aer_get_mode_props
public :: rad_aer_get_props
public :: rad_aer_get_bin_props_by_idx
public :: rad_aer_get_bin_props
public :: rad_aer_get_idx
public :: print_aerosol_lists
public :: rad_aer_readnl
public :: rad_aer_init

!==============================================================================
contains
!==============================================================================

function rad_aer_num_name(list_idx, spc_name_in, num_name_out, mode_out, spec_out ) result(found)
   use radiative_aerosol_definitions, only: modelist_t, modal_aerosol_list, modes

  ! for a given species name spc_name_in return (optionals):
  !   num_name_out -- corresponding number density species name
  !   mode_out -- corresponding mode number
  !   spec_out -- corresponding species number within the mode

  integer,           intent(in)  :: list_idx ! index of the climate or a diagnostic list
  character(len=*),  intent(in)  :: spc_name_in
  character(len=*),  intent(out) :: num_name_out
  integer, optional, intent(out) :: mode_out
  integer, optional, intent(out) :: spec_out

  logical :: found

  ! Local variables
  type(modelist_t), pointer :: m_list ! local pointer to mode list of interest
  integer :: n,m, mm
  integer :: nmodes
  integer :: nspecs
  character(len= 32) :: spec_name

  found = .false.

  m_list => modal_aerosol_list(list_idx)
  nmodes = m_list%nmodes

  do n = 1,nmodes
     mm = m_list%idx(n)
     nspecs = modes%comps(mm)%nspec
     do m = 1,nspecs
        spec_name = modes%comps(mm)%camname_mmr_a(m)
        if (spc_name_in == spec_name) then
           num_name_out = modes%comps(mm)%camname_num_a
           found = .true.
           if (present(mode_out)) then
              mode_out = n
           endif
           if (present(spec_out)) then
              spec_out = m
           endif
           return
        endif
     enddo
  enddo

end function rad_aer_num_name

!================================================================================================

subroutine rad_aer_get_info(list_idx, aernames, naero, nmodes, nbins)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: aerlist_t, modelist_t, binlist_t, &
      bulk_aerosol_list, modal_aerosol_list, sectional_aerosol_list

   ! Return info about aerosol lists (gas info handled in rad_constituents)

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=64), optional, intent(out) :: aernames(:)
   integer,           optional, intent(out) :: naero
   integer,           optional, intent(out) :: nmodes
   integer,           optional, intent(out) :: nbins

   ! Local variables
   type(aerlist_t),  pointer :: a_list ! local pointer to aerosol list of interest
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest
   type(binlist_t),  pointer :: s_list ! local pointer to bin list of interest

   integer          :: i
   integer          :: arrlen  ! length of assumed shape array

   character(len=*), parameter :: subname = 'rad_aer_get_info'
   !-----------------------------------------------------------------------------

   a_list => bulk_aerosol_list(list_idx)
   m_list => modal_aerosol_list(list_idx)
   s_list => sectional_aerosol_list(list_idx)

   ! number of bulk aerosols in list
   if (present(naero)) then
      naero = a_list%numaerosols
   endif

   ! number of aerosol modes in list
   if (present(nmodes)) then
      nmodes = m_list%nmodes
   endif

   ! number of aerosol bins in list
   if (present(nbins)) then
      nbins = s_list%nbins
   endif

   ! names of aerosols in list
   if (present(aernames)) then

      ! check that output array is long enough
      arrlen = size(aernames)
      if (arrlen < a_list%numaerosols) then
         write(iulog,*) subname//': ERROR: naero=', a_list%numaerosols, '  arrlen=', arrlen
         call endrun(subname//': ERROR: aernames too short')
      end if

      do i = 1, a_list%numaerosols
         aernames(i) = a_list%aer(i)%camname
      end do

   end if

end subroutine rad_aer_get_info

!================================================================================================

subroutine rad_aer_get_info_by_mode(list_idx, m_idx, &
   mode_type, num_name, num_name_cw, nspec)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: modelist_t, modal_aerosol_list, modes

   ! Return info about modal aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of mode in the specified list
   character(len=32), optional, intent(out) :: mode_type   ! type of mode (as used in MAM code)
   character(len=32), optional, intent(out) :: num_name    ! name of interstitial number mixing ratio
   character(len=32), optional, intent(out) :: num_name_cw ! name of cloud borne number mixing ratio
   integer,           optional, intent(out) :: nspec       ! number of species in the mode

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer          :: nmodes
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_aer_get_info_by_mode'
   !-----------------------------------------------------------------------------

   m_list => modal_aerosol_list(list_idx)

   ! check for valid mode index
   nmodes = m_list%nmodes
   if (m_idx < 1 .or. m_idx > nmodes) then
      write(iulog,*) subname//': ERROR - invalid mode index: ', m_idx
      call endrun(subname//': ERROR - invalid mode index')
   end if

   ! get index into the mode definition object
   mm = m_list%idx(m_idx)

   ! mode type
   if (present(mode_type)) then
      mode_type = modes%types(mm)
   endif

   ! number of species in the mode
   if (present(nspec)) then
      nspec = modes%comps(mm)%nspec
   endif

   ! name of interstitial number mixing ratio
   if (present(num_name)) then
      num_name = modes%comps(mm)%camname_num_a
   endif

   ! name of cloud borne number mixing ratio
   if (present(num_name_cw)) then
      num_name_cw = modes%comps(mm)%camname_num_c
   endif

end subroutine rad_aer_get_info_by_mode

!================================================================================================

subroutine rad_aer_get_info_by_bin(list_idx, m_idx, &
   bin_name, num_name, num_name_cw, mmr_name, mmr_name_cw, nspec)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: binlist_t, sectional_aerosol_list, bins

   ! Return info about CARMA aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of bin in the specified list
   character(len=*),  optional, intent(out) :: bin_name
   character(len=32), optional, intent(out) :: num_name    ! name of interstitial number mixing ratio
   character(len=32), optional, intent(out) :: num_name_cw ! name of cloud borne number mixing ratio
   character(len=32), optional, intent(out) :: mmr_name    ! name of interstitial mass mixing ratio
   character(len=32), optional, intent(out) :: mmr_name_cw ! name of cloud borne mass mixing ratio
   integer,           optional, intent(out) :: nspec       ! number of species in the mode

   ! Local variables
   type(binlist_t), pointer :: s_list ! local pointer to mode list of interest

   integer          :: nbins
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_aer_get_info_by_bin'
   !-----------------------------------------------------------------------------

   s_list => sectional_aerosol_list(list_idx)

   ! check for valid mode index
   nbins = s_list%nbins
   if (m_idx < 1 .or. m_idx > nbins) then
      write(iulog,*) subname//': ERROR - invalid bin index: ', m_idx
      call endrun(subname//': ERROR - invalid bin index')
   end if

   ! get index into the mode definition object
   mm = s_list%idx(m_idx)

   ! number of species in the mode
   if (present(nspec)) then
      nspec = bins%comps(mm)%nspec
   endif

   ! bin name
   if (present(bin_name)) then
      bin_name = bins%names(m_idx)
   end if

   ! name of interstitial number mixing ratio
   if (present(num_name)) then
      num_name = bins%comps(mm)%camname_num_a
   endif

   ! name of cloud borne number mixing ratio
   if (present(num_name_cw)) then
      num_name_cw = bins%comps(mm)%camname_num_c
   endif

   ! name of interstitial mass mixing ratio
   if (present(mmr_name)) then
      mmr_name = bins%comps(mm)%camname_mass_a
   endif

   ! name of cloud borne mass mixing ratio
   if (present(mmr_name_cw)) then
      mmr_name_cw = bins%comps(mm)%camname_mass_c
   endif

end subroutine rad_aer_get_info_by_bin

!================================================================================================
subroutine rad_aer_get_info_by_bin_spec(list_idx, m_idx, s_idx, &
   spec_type, spec_morph, spec_name, spec_name_cw)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: binlist_t, sectional_aerosol_list, bins

   ! Return info about CARMA aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of bin in the specified list
   integer,                     intent(in)  :: s_idx       ! index of species in the specified mode
   character(len=32), optional, intent(out) :: spec_type   ! type of species
   character(len=32), optional, intent(out) :: spec_morph  ! type of species
   character(len=32), optional, intent(out) :: spec_name   ! name of interstitial species
   character(len=32), optional, intent(out) :: spec_name_cw ! name of cloud borne species

   ! Local variables
   type(binlist_t), pointer :: s_list ! local pointer to mode list of interest
   integer          :: nbins,  nspec
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_aer_get_info_by_bin_spec'
   !-----------------------------------------------------------------------------

   s_list => sectional_aerosol_list(list_idx)

   ! check for valid mode index
   nbins = s_list%nbins
   if (m_idx < 1 .or. m_idx > nbins) then
      write(iulog,*) subname//': ERROR - invalid bin index: ', m_idx
      call endrun(subname//': ERROR - invalid bin index')
   end if

   ! get index into the mode definition object
   mm = s_list%idx(m_idx)

   ! check for valid species index
   nspec = bins%comps(mm)%nspec
   if (s_idx < 1 .or. s_idx > nspec) then
      write(iulog,*) subname//': ERROR - invalid specie index: ', s_idx
      call endrun(subname//': ERROR - invalid specie index')
   end if

   if (present(spec_type)) then
      spec_type = bins%comps(mm)%type(s_idx)
   endif
   if (present(spec_morph)) then
      spec_morph = bins%comps(mm)%morph(s_idx)
   endif
   if (present(spec_name)) then
      spec_name = bins%comps(mm)%camname_mmr_a(s_idx)
   endif
   if (present(spec_name_cw)) then
      spec_name_cw = bins%comps(mm)%camname_mmr_c(s_idx)
   endif

end subroutine rad_aer_get_info_by_bin_spec

!================================================================================================
subroutine rad_aer_get_info_by_mode_spec(list_idx, m_idx, s_idx, &
   spec_type, spec_name, spec_name_cw)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: modelist_t, modal_aerosol_list, modes

   ! Return info about modal aerosol lists

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: m_idx       ! index of mode in the specified list
   integer,                     intent(in)  :: s_idx       ! index of specie in the specified mode
   character(len=32), optional, intent(out) :: spec_type   ! type of specie
   character(len=32), optional, intent(out) :: spec_name   ! name of interstitial specie
   character(len=32), optional, intent(out) :: spec_name_cw ! name of cloud borne specie

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer          :: nmodes
   integer          :: nspec
   integer          :: mm

   character(len=*), parameter :: subname = 'rad_aer_get_info_by_mode_spec'
   !-----------------------------------------------------------------------------

   m_list => modal_aerosol_list(list_idx)

   ! check for valid mode index
   nmodes = m_list%nmodes
   if (m_idx < 1 .or. m_idx > nmodes) then
      write(iulog,*) subname//': ERROR - invalid mode index: ', m_idx
      call endrun(subname//': ERROR - invalid mode index')
   end if

   ! get index into the mode definition object
   mm = m_list%idx(m_idx)

   ! check for valid specie index
   nspec = modes%comps(mm)%nspec
   if (s_idx < 1 .or. s_idx > nspec) then
      write(iulog,*) subname//': ERROR - invalid specie index: ', s_idx
      call endrun(subname//': ERROR - invalid specie index')
   end if

   ! specie type
   if (present(spec_type)) then
      spec_type = modes%comps(mm)%type(s_idx)
   endif

   ! interstitial specie name
   if (present(spec_name)) then
      spec_name = modes%comps(mm)%camname_mmr_a(s_idx)
   endif

   ! cloud borne specie name
   if (present(spec_name_cw)) then
      spec_name_cw = modes%comps(mm)%camname_mmr_c(s_idx)
   endif

end subroutine rad_aer_get_info_by_mode_spec

!================================================================================================

subroutine rad_aer_get_info_by_spectype(list_idx, spectype, mode_idx, spec_idx)
   use radiative_aerosol_definitions, only: modelist_t, modal_aerosol_list, modes

   ! Return info about modes in the specified climate/diagnostics list

   ! Arguments
   integer,                     intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),            intent(in)  :: spectype    ! species type
   integer,           optional, intent(out) :: mode_idx    ! index of a mode that contains a specie of spectype
   integer,           optional, intent(out) :: spec_idx    ! index of the species of spectype

   ! Local variables
   type(modelist_t), pointer :: m_list ! local pointer to mode list of interest

   integer  :: i, nmodes, m_idx, nspec, ispec
   logical  :: found_spectype

   character(len=*), parameter :: subname = 'rad_aer_get_info_by_spectype'
   !-----------------------------------------------------------------------------

   m_list => modal_aerosol_list(list_idx)

   ! number of modes in specified list
   nmodes = m_list%nmodes

   ! loop through modes in specified climate/diagnostic list
   found_spectype = .false.
   do i = 1, nmodes

      ! get index of the mode in the definition object
      m_idx = m_list%idx(i)

      ! number of species in the mode
      nspec = modes%comps(m_idx)%nspec

      ! loop through species looking for spectype
      do ispec = 1, nspec

         if (trim(modes%comps(m_idx)%type(ispec)) == trim(spectype)) then
            if (present(mode_idx)) mode_idx = i
            if (present(spec_idx)) spec_idx = ispec
            found_spectype = .true.
            exit
         end if
      end do

      if (found_spectype) exit
   end do

   if (.not. found_spectype) then
      if (present(mode_idx)) mode_idx = -1
      if (present(spec_idx)) spec_idx = -1
   end if

end subroutine rad_aer_get_info_by_spectype

!================================================================================================

function rad_aer_get_mode_idx(list_idx, mode_type) result(mode_idx)
   use radiative_aerosol_definitions, only: modelist_t, modal_aerosol_list, modes

   ! Return mode index of the specified type in the specified climate/diagnostics list.
   ! Return -1 if not found.

   ! Arguments
   integer,           intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),  intent(in)  :: mode_type   ! mode type

   ! Return value
   integer                        :: mode_idx    ! mode index

   ! Local variables
   type(modelist_t), pointer :: m_list

   integer  :: i, nmodes, m_idx

   character(len=*), parameter :: subname = 'rad_aer_get_mode_idx'
   !-----------------------------------------------------------------------------

   ! if mode type not found return -1
   mode_idx = -1

   ! specified mode list
   m_list => modal_aerosol_list(list_idx)

   ! number of modes in specified list
   nmodes = m_list%nmodes

   ! loop through modes in specified climate/diagnostic list
   do i = 1, nmodes

      ! get index of the mode in the definition object
      m_idx = m_list%idx(i)

      ! look in mode definition object (modes) for the mode types
      if (trim(modes%types(m_idx)) == trim(mode_type)) then
         mode_idx = i
         exit
      end if
   end do

end function rad_aer_get_mode_idx

!================================================================================================

function rad_aer_get_spec_idx(list_idx, mode_idx, spec_type) result(spec_idx)
   use radiative_aerosol_definitions, only: modelist_t, mode_component_t, modal_aerosol_list, modes

   ! Return specie index of the specified type in the specified mode of the specified
   ! climate/diagnostics list.  Return -1 if not found.

   ! Arguments
   integer,           intent(in)  :: list_idx    ! index of the climate or a diagnostic list
   integer,           intent(in)  :: mode_idx    ! mode index
   character(len=*),  intent(in)  :: spec_type   ! specie type

   ! Return value
   integer                        :: spec_idx    ! specie index

   ! Local variables
   type(modelist_t),       pointer :: m_list
   type(mode_component_t), pointer :: mode_comps

   integer  :: i, m_idx, nspec

   character(len=*), parameter :: subname = 'rad_aer_get_spec_idx'
   !-----------------------------------------------------------------------------

   ! if specie type not found return -1
   spec_idx = -1

   ! modes in specified list
   m_list => modal_aerosol_list(list_idx)

   ! get index of the specified mode in the definition object
   m_idx = m_list%idx(mode_idx)

   ! object containing the components of the mode
   mode_comps => modes%comps(m_idx)

   ! number of species in specified mode
   nspec = mode_comps%nspec

   ! loop through species in specified mode
   do i = 1, nspec

      ! look in mode definition object (modes) for the mode types
      if (trim(mode_comps%type(i)) == trim(spec_type)) then
         spec_idx = i
         exit
      end if
   end do

end function rad_aer_get_spec_idx

!================================================================================================

integer function rad_aer_get_idx(list_idx, aer_name)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: N_DIAG, aerlist_t, bulk_aerosol_list

   ! Return the index of aerosol aer_name in the list specified by list_idx.

    ! Arguments
   integer,             intent(in) :: list_idx    ! 0 for climate list, 1-N_DIAG for diagnostic lists
   character(len=*),    intent(in) :: aer_name    ! aerosol name (in state or pbuf)

   ! Local variables
   integer :: i, aer_idx
   type(aerlist_t), pointer :: aerlist
   character(len=*), parameter :: subname = "rad_aer_get_idx"
   !-------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => bulk_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   ! Get index in aerosol list for requested name
   aer_idx = -1
   do i = 1, aerlist%numaerosols
      if (trim(aer_name) == trim(aerlist%aer(i)%camname)) then
         aer_idx = i
         exit
      end if
   end do

   if (aer_idx == -1) call endrun(subname//": ERROR - name not found")

   rad_aer_get_idx = aer_idx

end function rad_aer_get_idx

!================================================================================================

subroutine rad_aer_get_props_by_idx(list_idx, &
   aer_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, num_to_mass_aer)
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use phys_prop,      only: physprop_get, ot_length
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: N_DIAG, aerlist_t, bulk_aerosol_list

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in)  :: list_idx ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: aer_idx  ! index of the aerosol
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:)
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:)
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)
   character(len=20), optional, intent(out) :: aername
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer
   real(r8),          optional, intent(out) :: dispersion_aer
   real(r8),          optional, intent(out) :: num_to_mass_aer

   real(r8),          optional, pointer     :: r_sw_ext(:,:)
   real(r8),          optional, pointer     :: r_sw_scat(:,:)
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)
   real(r8),          optional, pointer     :: r_lw_abs(:,:)
   real(r8),          optional, pointer     :: mu(:)

   ! Local variables
   integer :: idx
   character(len=*), parameter :: subname = 'rad_aer_get_props_by_idx'
   type(aerlist_t), pointer :: aerlist
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      aerlist => bulk_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   if (aer_idx < 1 .or. aer_idx > aerlist%numaerosols) then
      write(iulog,*) subname//': aerosol list index out of range: ', aer_idx ,' list index: ',list_idx
      call endrun(subname//': aer_idx out of range')
   end if

   idx = aerlist%aer(aer_idx)%physprop_id

   if (present(opticstype))        call physprop_get(idx, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(idx, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(idx, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(idx, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(idx, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(idx, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(idx, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(idx, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(idx, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(idx, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(idx, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(idx, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(idx, refindex_aer_lw=refindex_aer_lw)

   if (present(aername))           call physprop_get(idx, aername=aername)
   if (present(density_aer))       call physprop_get(idx, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(idx, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(idx, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(idx, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(idx, num_to_mass_aer=num_to_mass_aer)

   if (present(r_lw_abs))          call physprop_get(idx, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(idx, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(idx, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(idx, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(idx, mu=mu)

end subroutine rad_aer_get_props_by_idx

!================================================================================================

subroutine rad_aer_get_mam_props_by_idx(list_idx, &
   mode_idx, spec_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, &
   num_to_mass_aer, spectype)
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use phys_prop,      only: physprop_get, ot_length
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: N_DIAG, modelist_t, modal_aerosol_list, modes

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: mode_idx  ! mode index
   integer,                     intent(in)  :: spec_idx  ! index of specie in the mode
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:)
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:)
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)

   real(r8),          optional, pointer     :: r_sw_ext(:,:)
   real(r8),          optional, pointer     :: r_sw_scat(:,:)
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)
   real(r8),          optional, pointer     :: r_lw_abs(:,:)
   real(r8),          optional, pointer     :: mu(:)

   character(len=20), optional, intent(out) :: aername
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer
   real(r8),          optional, intent(out) :: dispersion_aer
   real(r8),          optional, intent(out) :: num_to_mass_aer
   character(len=32), optional, intent(out) :: spectype

   ! Local variables
   integer :: m_idx, idx
   type(modelist_t), pointer   :: mlist
   character(len=*), parameter :: subname = 'rad_aer_get_mam_props_by_idx'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => modal_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
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

   idx = modes%comps(m_idx)%idx_props(spec_idx)

   if (present(opticstype))        call physprop_get(idx, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(idx, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(idx, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(idx, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(idx, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(idx, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(idx, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(idx, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(idx, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(idx, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(idx, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(idx, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(idx, refindex_aer_lw=refindex_aer_lw)

   if (present(r_lw_abs))          call physprop_get(idx, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(idx, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(idx, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(idx, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(idx, mu=mu)

   if (present(aername))           call physprop_get(idx, aername=aername)
   if (present(density_aer))       call physprop_get(idx, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(idx, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(idx, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(idx, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(idx, num_to_mass_aer=num_to_mass_aer)

   if (present(spectype)) spectype = modes%comps(m_idx)%type(spec_idx)

end subroutine rad_aer_get_mam_props_by_idx

!================================================================================================

subroutine rad_aer_get_bin_props_by_idx(list_idx, &
   bin_idx, spec_idx,  opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, &
   num_to_mass_aer, spectype, specmorph)
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use phys_prop,      only: physprop_get, ot_length
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: N_DIAG, binlist_t, sectional_aerosol_list, bins

   ! Return requested properties for the aerosol from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,                     intent(in)  :: bin_idx   ! mode index
   integer,                     intent(in)  :: spec_idx  ! index of specie in the mode
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:)
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:)
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)
   complex(r8),       optional, pointer     :: refindex_aer_sw(:)
   complex(r8),       optional, pointer     :: refindex_aer_lw(:)

   real(r8),          optional, pointer     :: r_sw_ext(:,:)
   real(r8),          optional, pointer     :: r_sw_scat(:,:)
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)
   real(r8),          optional, pointer     :: r_lw_abs(:,:)
   real(r8),          optional, pointer     :: mu(:)

   character(len=20), optional, intent(out) :: aername
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer
   real(r8),          optional, intent(out) :: dispersion_aer
   real(r8),          optional, intent(out) :: num_to_mass_aer
   character(len=32), optional, intent(out) :: spectype
   character(len=32), optional, intent(out) :: specmorph

   ! Local variables
   integer :: m_idx, idx
   type(binlist_t), pointer :: slist
   character(len=*), parameter :: subname = 'rad_aer_get_bin_props_by_idx'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sectional_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the index for the corresponding mode in the mode definition object
   m_idx = slist%idx(bin_idx)

   ! Check for valid specie index
   if (spec_idx < 1  .or.  spec_idx > bins%comps(m_idx)%nspec) then
      write(iulog,*) subname//': spec_idx= ', spec_idx, '  nspec= ', bins%comps(m_idx)%nspec
      call endrun(subname//': specie list index out of range')
   end if

   idx = bins%comps(m_idx)%idx_props(spec_idx)

   if (present(opticstype))        call physprop_get(idx, opticstype=opticstype)

   if (present(sw_hygro_ext))      call physprop_get(idx, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(idx, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(idx, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(idx, lw_hygro_abs=lw_hygro_ext)

   if (present(sw_nonhygro_ext))   call physprop_get(idx, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(idx, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(idx, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(idx, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(idx, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(idx, lw_abs=lw_ext)

   if (present(refindex_aer_sw))   call physprop_get(idx, refindex_aer_sw=refindex_aer_sw)
   if (present(refindex_aer_lw))   call physprop_get(idx, refindex_aer_lw=refindex_aer_lw)

   if (present(r_lw_abs))          call physprop_get(idx, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(idx, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(idx, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(idx, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(idx, mu=mu)

   if (present(aername))           call physprop_get(idx, aername=aername)
   if (present(density_aer))       call physprop_get(idx, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(idx, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(idx, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(idx, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(idx, num_to_mass_aer=num_to_mass_aer)

   if (present(spectype)) spectype = bins%comps(m_idx)%type(spec_idx)
   if (present(specmorph)) specmorph = bins%comps(m_idx)%morph(spec_idx)

end subroutine rad_aer_get_bin_props_by_idx

!================================================================================================

subroutine rad_aer_get_mode_props(list_idx, mode_idx, opticstype, &
   extpsw, abspsw, asmpsw, absplw, refrtabsw, &
   refitabsw, refrtablw, refitablw, ncoef, prefr, &
   prefi, sigmag, dgnum, dgnumlo, dgnumhi, &
   rhcrystal, rhdeliques)

   use shr_kind_mod,   only: r8 => shr_kind_r8
   use phys_prop,      only: physprop_get, ot_length
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: N_DIAG, modelist_t, modal_aerosol_list

   ! Return requested properties for the mode from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,             intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,             intent(in)  :: mode_idx  ! mode index

   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),  optional, pointer     :: extpsw(:,:,:,:)
   real(r8),  optional, pointer     :: abspsw(:,:,:,:)
   real(r8),  optional, pointer     :: asmpsw(:,:,:,:)
   real(r8),  optional, pointer     :: absplw(:,:,:,:)
   real(r8),  optional, pointer     :: refrtabsw(:,:)
   real(r8),  optional, pointer     :: refitabsw(:,:)
   real(r8),  optional, pointer     :: refrtablw(:,:)
   real(r8),  optional, pointer     :: refitablw(:,:)
   integer,   optional, intent(out) :: ncoef
   integer,   optional, intent(out) :: prefr
   integer,   optional, intent(out) :: prefi
   real(r8),  optional, intent(out) :: sigmag
   real(r8),  optional, intent(out) :: dgnum
   real(r8),  optional, intent(out) :: dgnumlo
   real(r8),  optional, intent(out) :: dgnumhi
   real(r8),  optional, intent(out) :: rhcrystal
   real(r8),  optional, intent(out) :: rhdeliques

   ! Local variables
   integer :: idx
   type(modelist_t), pointer :: mlist
   character(len=*), parameter :: subname = 'rad_aer_get_mode_props'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      mlist => modal_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (mode_idx < 1  .or.  mode_idx > mlist%nmodes) then
      write(iulog,*) subname//': mode_idx= ', mode_idx, '  nmodes= ', mlist%nmodes
      call endrun(subname//': mode list index out of range')
   end if

   ! Get the physprop index for the requested mode
   idx = mlist%idx_props(mode_idx)

   if (present(opticstype))  call physprop_get(idx, opticstype=opticstype)
   if (present(extpsw))      call physprop_get(idx, extpsw=extpsw)
   if (present(abspsw))      call physprop_get(idx, abspsw=abspsw)
   if (present(asmpsw))      call physprop_get(idx, asmpsw=asmpsw)
   if (present(absplw))      call physprop_get(idx, absplw=absplw)

   if (present(refrtabsw))   call physprop_get(idx, refrtabsw=refrtabsw)
   if (present(refitabsw))   call physprop_get(idx, refitabsw=refitabsw)
   if (present(refrtablw))   call physprop_get(idx, refrtablw=refrtablw)
   if (present(refitablw))   call physprop_get(idx, refitablw=refitablw)

   if (present(ncoef))       call physprop_get(idx, ncoef=ncoef)
   if (present(prefr))       call physprop_get(idx, prefr=prefr)
   if (present(prefi))       call physprop_get(idx, prefi=prefi)
   if (present(sigmag))      call physprop_get(idx, sigmag=sigmag)
   if (present(dgnum))       call physprop_get(idx, dgnum=dgnum)
   if (present(dgnumlo))     call physprop_get(idx, dgnumlo=dgnumlo)
   if (present(dgnumhi))     call physprop_get(idx, dgnumhi=dgnumhi)
   if (present(rhcrystal))   call physprop_get(idx, rhcrystal=rhcrystal)
   if (present(rhdeliques))  call physprop_get(idx, rhdeliques=rhdeliques)

end subroutine rad_aer_get_mode_props

!================================================================================================

subroutine rad_aer_get_bin_props(list_idx, bin_idx, opticstype, &
   extpsw, abspsw, asmpsw, absplw, corefrac, nfrac, &
   wgtpct, nwtp, bcdust, nbcdust, kap, nkap, relh, nrelh, &
   sw_hygro_ext_wtp, sw_hygro_ssa_wtp, sw_hygro_asm_wtp, lw_hygro_ext_wtp, &
   sw_hygro_coreshell_ext, sw_hygro_coreshell_ssa, sw_hygro_coreshell_asm, lw_hygro_coreshell_ext, dryrad )
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use phys_prop,      only: physprop_get, ot_length
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: N_DIAG, binlist_t, sectional_aerosol_list

   ! Return requested properties for the bin from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,             intent(in)  :: list_idx  ! index of the climate or a diagnostic list
   integer,             intent(in)  :: bin_idx   ! mode index

   character(len=ot_length), optional, intent(out) :: opticstype

   real(r8),  optional, pointer     :: extpsw(:,:)
   real(r8),  optional, pointer     :: abspsw(:,:)
   real(r8),  optional, pointer     :: asmpsw(:,:)
   real(r8),  optional, pointer     :: absplw(:,:)
   real(r8),  optional, pointer     :: corefrac(:)
   integer,   optional, intent(out) :: nfrac

   real(r8),  optional, pointer     :: sw_hygro_ext_wtp(:,:)
   real(r8),  optional, pointer     :: sw_hygro_ssa_wtp(:,:)
   real(r8),  optional, pointer     :: sw_hygro_asm_wtp(:,:)
   real(r8),  optional, pointer     :: lw_hygro_ext_wtp(:,:)
   real(r8),  optional, pointer     :: sw_hygro_coreshell_ext(:,:,:,:,:)
   real(r8),  optional, pointer     :: sw_hygro_coreshell_ssa(:,:,:,:,:)
   real(r8),  optional, pointer     :: sw_hygro_coreshell_asm(:,:,:,:,:)
   real(r8),  optional, pointer     :: lw_hygro_coreshell_ext(:,:,:,:,:)
   real(r8),  optional, pointer     :: wgtpct(:)
   real(r8),  optional, pointer     :: bcdust(:)
   real(r8),  optional, pointer     :: kap(:)
   real(r8),  optional, pointer     :: relh(:)
   integer,   optional, intent(out) :: nwtp
   integer,   optional, intent(out) :: nbcdust
   integer,   optional, intent(out) :: nkap
   integer,   optional, intent(out) :: nrelh
   real(r8),  optional, intent(out) :: dryrad

   ! Local variables
   integer :: idx
   type(binlist_t),  pointer   :: slist
   character(len=*), parameter :: subname = 'rad_aer_get_bin_props'
   !------------------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      slist => sectional_aerosol_list(list_idx)
   else
      write(iulog,*) subname//': list_idx = ', list_idx
      call endrun(subname//': list_idx out of range')
   endif

   ! Check for valid mode index
   if (bin_idx < 1  .or.  bin_idx > slist%nbins) then
      write(iulog,*) subname//': bin_idx= ', bin_idx, '  nbins= ', slist%nbins
      call endrun(subname//': bin list index out of range')
   end if

   ! Get the physprop index for the requested bin
   idx = slist%idx_props(bin_idx)

   if (present(opticstype))  call physprop_get(idx, opticstype=opticstype)
   if (present(extpsw))      call physprop_get(idx, extpsw2=extpsw)
   if (present(abspsw))      call physprop_get(idx, abspsw2=abspsw)
   if (present(asmpsw))      call physprop_get(idx, asmpsw2=asmpsw)
   if (present(absplw))      call physprop_get(idx, absplw2=absplw)
   if (present(corefrac))    call physprop_get(idx, corefrac=corefrac)
   if (present(nfrac))       call physprop_get(idx, nfrac=nfrac)

   if (present(sw_hygro_ext_wtp))       call physprop_get(idx, sw_hygro_ext_wtp=sw_hygro_ext_wtp)
   if (present(sw_hygro_ssa_wtp))       call physprop_get(idx, sw_hygro_ssa_wtp=sw_hygro_ssa_wtp)
   if (present(sw_hygro_asm_wtp))       call physprop_get(idx, sw_hygro_asm_wtp=sw_hygro_asm_wtp)
   if (present(lw_hygro_ext_wtp))       call physprop_get(idx, lw_hygro_abs_wtp=lw_hygro_ext_wtp)
   if (present(sw_hygro_coreshell_ext)) call physprop_get(idx, sw_hygro_coreshell_ext=sw_hygro_coreshell_ext)
   if (present(sw_hygro_coreshell_ssa)) call physprop_get(idx, sw_hygro_coreshell_ssa=sw_hygro_coreshell_ssa)
   if (present(sw_hygro_coreshell_asm)) call physprop_get(idx, sw_hygro_coreshell_asm=sw_hygro_coreshell_asm)
   if (present(lw_hygro_coreshell_ext)) call physprop_get(idx, lw_hygro_coreshell_abs=lw_hygro_coreshell_ext)
   if (present(wgtpct))                 call physprop_get(idx, wgtpct=wgtpct)
   if (present(bcdust))                 call physprop_get(idx, bcdust=bcdust)
   if (present(kap))                    call physprop_get(idx, kap=kap)
   if (present(relh))                   call physprop_get(idx, relh=relh)
   if (present(nwtp))                   call physprop_get(idx, nwtp=nwtp)
   if (present(nbcdust))                call physprop_get(idx, nbcdust=nbcdust)
   if (present(nkap))                   call physprop_get(idx, nkap=nkap)
   if (present(nrelh))                  call physprop_get(idx, nrelh=nrelh)
   if (present(dryrad))                 call physprop_get(idx, dryrad_aer=dryrad)

end subroutine rad_aer_get_bin_props

!================================================================================================

subroutine print_aerosol_lists(aer_list, m_list, s_list)
   use cam_logfile,    only: iulog
   use radiative_aerosol_definitions, only: nl, aerlist_t, modelist_t, binlist_t, modes, bins

   ! Print summary of bulk, modal, and bin aerosol lists.

   type(aerlist_t),  intent(in) :: aer_list
   type(modelist_t), intent(in) :: m_list
   type(binlist_t),  intent(in) :: s_list

   integer :: i, idx

   if (len_trim(aer_list%list_id) == 0) then
      write(iulog,*) nl//' bulk aerosol list for climate calculations'
   else
      write(iulog,*) nl//' bulk aerosol list for diag'//aer_list%list_id//' calculations'
   end if

   do i = 1, aer_list%numaerosols
      write(iulog,*) '  '//trim(aer_list%aer(i)%source)//':'//trim(aer_list%aer(i)%camname)//&
                     ' optics and phys props in :'//trim(aer_list%aer(i)%physprop_file)
   enddo

   if (len_trim(m_list%list_id) == 0) then
      write(iulog,*) nl//' modal aerosol list for climate calculations'
   else
      write(iulog,*) nl//' modal aerosol list for diag'//m_list%list_id//' calculations'
   end if

   do i = 1, m_list%nmodes
      idx = m_list%idx(i)
      write(iulog,*) '  '//trim(modes%names(idx))
   enddo

   if (len_trim(s_list%list_id) == 0) then
      write(iulog,*) nl//' bin aerosol list for climate calculations'
   else
      write(iulog,*) nl//' bin aerosol list for diag'//s_list%list_id//' calculations'
   end if

   do i = 1, s_list%nbins
      idx = s_list%idx(i)
      write(iulog,*) '  '//trim(bins%names(idx))
   enddo

end subroutine print_aerosol_lists

!================================================================================================

! Parse aerosol mode/bin definitions, accumulate physprop files,
! and initialize aerosol lists (phase 1).
!
! Called from rad_cnst_readnl after namelist I/O, broadcast, and
! parse_rad_specifier / active_calls have been set.
!
! In SIMA, this will read aerosol-specific namelists directly
! (rad_aerosol / rad_aer_diag_N instead of rad_climate / rad_diag_N).
subroutine rad_aer_readnl(mode_defs, bin_defs)
   use phys_prop,      only: physprop_accum_unique_files
   use spmd_utils,     only: masterproc
   use radiative_aerosol_definitions, only: &
      cs1, verbose, N_DIAG, modes, bins, &
      active_calls, bulk_aerosol_list, modal_aerosol_list, sectional_aerosol_list, &
      radcnst_namelist, parse_mode_defs, parse_bin_defs, &
      list_populate, print_modes, print_bins

   ! Arguments
   character(len=cs1), intent(inout) :: mode_defs(:)
   character(len=cs1), intent(inout) :: bin_defs(:)

   ! Local variables
   integer :: i
   character(len=2) :: suffix
   character(len=1), pointer   :: ctype(:)
   character(len=*), parameter :: subname = 'rad_aer_readnl'
   !-----------------------------------------------------------------------------

   ! Parse mode definition strings
   call parse_mode_defs(mode_defs, modes)

   ! Parse bin definition strings
   call parse_bin_defs(bin_defs, bins)

   ! Set the list_id fields for aerosol lists
   do i = 0, N_DIAG
      if (active_calls(i)) then
         if (i > 0) then
            write(suffix, fmt = '(i2.2)') i
         else
            suffix='  '
         end if
         bulk_aerosol_list(i)%list_id      = suffix
         modal_aerosol_list(i)%list_id     = suffix
         sectional_aerosol_list(i)%list_id = suffix
      end if
   end do

   ! Accumulate unique physprop files — bulk aerosol species
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call physprop_accum_unique_files(radcnst_namelist(i)%radname, radcnst_namelist(i)%type)
      endif
   enddo

   ! Accumulate physprop files for mode species
   do i = 1, modes%nmodes
      allocate(ctype(modes%comps(i)%nspec))
      ctype = 'A'
      call physprop_accum_unique_files(modes%comps(i)%props, ctype)
      deallocate(ctype)
   end do

   ! Accumulate physprop files for bin species
   do i = 1, bins%nbins
      allocate(ctype(bins%comps(i)%nspec))
      ctype = 'A'
      call physprop_accum_unique_files(bins%comps(i)%props, ctype)
      deallocate(ctype)
   end do

   ! Initialize aerosol lists (populate from namelist specifiers)
   do i = 0, N_DIAG
      if (active_calls(i)) then
         ! has to be done at readnl phase as information on structure of the lists will be needed
         ! in physics/chemistry initialization.
         call list_populate(radcnst_namelist(i), bulk_aerosol_list(i), modal_aerosol_list(i), sectional_aerosol_list(i))

         if (masterproc .and. verbose) then
            call print_aerosol_lists(bulk_aerosol_list(i), modal_aerosol_list(i), sectional_aerosol_list(i))
         end if
      end if
   end do

   if (masterproc .and. verbose) call print_modes(modes)
   if (masterproc .and. verbose) call print_bins(bins)

end subroutine rad_aer_readnl

!================================================================================================

! Complete aerosol initialization (phase 2).
! Reads physprop files, resolves constituent indices for modes/bins,
! finishes aerosol list init, and registers aerosol diagnostic fields.
!
! Called from physpkg before rad_cnst_init (gas init).
subroutine rad_aer_init()
   use phys_prop,      only: physprop_init
   use radiative_aerosol_definitions, only: &
      N_DIAG, modes, bins, active_calls, &
      bulk_aerosol_list, modal_aerosol_list, sectional_aerosol_list, list_resolve_physprops

   !REMOVECAM: aerosol_mmr_cam handles CAM-specific index resolution
   use aerosol_mmr_cam, only: aerosol_mmr_init, &
      resolve_mode_idx, resolve_bin_idx, resolve_bulk_idx, &
      rad_aer_diag_init
   !REMOVECAM_END

   integer :: i
   character(len=*), parameter :: subname = 'rad_aer_init'
   !-----------------------------------------------------------------------------

   ! Initialize a zero target for the 'Z' type of aerosol MMR.
   call aerosol_mmr_init()

   ! Read physical properties from data files
   call physprop_init()

   !REMOVECAM: resolve host-specific indices (CAM uses pbuf and state)
   call resolve_mode_idx(modes)
   call resolve_bin_idx(bins)
   !REMOVECAM_END

   ! Resolve physprop indices for aerosol lists
   do i = 0, N_DIAG
      if (active_calls(i)) then
         !REMOVECAM: resolve host-specific indices (CAM uses pbuf and state)
         call resolve_bulk_idx(bulk_aerosol_list(i))
         !REMOVECAM_END
         call list_resolve_physprops(bulk_aerosol_list(i), modal_aerosol_list(i), sectional_aerosol_list(i))
      end if
   end do

   !REMOVECAM: history add calls for radiative aerosol diagnostics.
   call rad_aer_diag_init(bulk_aerosol_list(0))
   !REMOVECAM_END

end subroutine rad_aer_init

!================================================================================================

end module radiative_aerosol
