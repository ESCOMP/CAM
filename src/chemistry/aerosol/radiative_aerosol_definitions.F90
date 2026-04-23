!-----------------------------------------------------------------------------
! Core aerosol definitions for radiative calculations: shared constants,
! types, data, parsing, and initialization routines for both modal and
! sectional (bin) aerosol representations.
!
! This module is the lowest-level shared module in the aerosol hierarchy.
! It will be shared with CAM-SIMA.
!-----------------------------------------------------------------------------
module radiative_aerosol_definitions

  use shr_kind_mod,   only: shr_kind_cl

  implicit none
  private
  save

  public :: parse_mode_defs, parse_bin_defs  ! parse mode and bin definitions for aerosol.
  public :: parse_rad_specifier              ! parse rad_climate and rad_diag_N specifiers into rad_cnst_namelist_t.
  public :: list_populate                    ! populate aerosol list structures from parsed namelist (run before register)
  public :: list_resolve_physprops           ! resolve physprop indices into aerosol list structures
  public :: print_modes, print_bins

  !===========================
  ! Named constants for mode/species/morph validation
  ! These categories and definitions are used throughout the aerosol models,
  ! not just in radiative_aerosol.
  !===========================
  integer, public, parameter :: num_mode_types = 9
  integer, public, parameter :: num_spec_types = 8
  character(len=14), public, parameter :: mode_type_names(num_mode_types) = (/ &
     'accum         ', 'aitken        ', 'primary_carbon', 'fine_seasalt  ', &
     'fine_dust     ', 'coarse        ', 'coarse_seasalt', 'coarse_dust   ', &
     'coarse_strat  '  /)
  character(len=9),  public, parameter :: spec_type_names(num_spec_types) = (/ &
     'sulfate  ', 'ammonium ', 'nitrate  ', 'p-organic', &
     's-organic', 'black-c  ', 'seasalt  ', 'dust     '/)

  integer, public, parameter :: num_bin_morphs  = 2
  character(len=8), public, parameter :: bin_morph_names(num_bin_morphs) = &
       (/ 'shell   ', 'core    ' /)

  !===========================
  ! Shared constants (shared with rad_constituents for gases) part 1.
  !===========================
  logical,          public, parameter :: verbose = .true.
  character(len=1), public, parameter :: newline = achar(10)

  !===========================
  ! Types
  !===========================
!! \section arg_table_rad_cnst_namelist_t
!! \htmlinclude rad_cnst_namelist_t.html
  ! type to provide access to the data parsed from the rad_climate and rad_diag_* strings
  type, public :: rad_cnst_namelist_t
     integer :: ncnst
     character(len=  1),         pointer :: source(:)  ! 'A' for state (advected), 'N' for pbuf (non-advected),
                                               ! 'M' for mode, 'Z' for zero
     character(len= 64),         pointer :: camname(:) ! name registered in pbuf or constituents
     character(len=shr_kind_cl), pointer :: radname(:) ! radname is the name as identfied in radiation,
                                               ! must be one of (rgaslist if a gas) or
                                               ! (/fullpath/filename.nc if an aerosol)
     character(len=  1),         pointer :: type(:)    ! 'A' if aerosol, 'G' if gas, 'M' if mode
  end type rad_cnst_namelist_t

!! \section arg_table_mode_component_t
!! \htmlinclude mode_component_t.html
  ! type to provide access to the components of a mode
  type, public :: mode_component_t
     integer :: nspec
     ! For "source" variables below, value is:
     ! 'N' if in pbuf (non-advected)
     ! 'A' if in state (advected)

     ! source of interstitial number conc field
     character(len=  1) :: source_num_a
     ! name registered in pbuf or constituents for number mixing ratio of interstitial species
     character(len= 32) :: camname_num_a
     ! source of cloud borne number conc field
     character(len=  1) :: source_num_c
     ! name registered in pbuf or constituents for number mixing ratio of cloud borne species
     character(len= 32) :: camname_num_c
     ! source of interstitial specie mmr fields
     character(len=  1),         pointer :: source_mmr_a(:)
     ! name registered in pbuf or constituents for mmr of interstitial components
     character(len= 32),         pointer :: camname_mmr_a(:)
     ! source of cloud borne specie mmr fields
     character(len=  1),         pointer :: source_mmr_c(:)
     ! name registered in pbuf or constituents for mmr of cloud borne components
     character(len= 32),         pointer :: camname_mmr_c(:)
     ! specie type (as used in MAM code)
     character(len= 32),         pointer :: type(:)
     ! file containing specie properties
     character(len=shr_kind_cl), pointer :: props(:)

     ! index in pbuf or constituents for number mixing ratio of interstitial species
     integer          :: idx_num_a
     ! index in pbuf for number mixing ratio of interstitial species
     integer          :: idx_num_c
     ! index in pbuf or constituents for mmr of interstitial species
     integer, pointer :: idx_mmr_a(:)
     ! index in pbuf for mmr of interstitial species
     integer, pointer :: idx_mmr_c(:)
     ! ID used to access physical properties of mode species from phys_prop module
     integer, pointer :: idx_props(:)
  end type mode_component_t

!! \section arg_table_modes_t
!! \htmlinclude modes_t.html
  ! type to provide access to all modes
  type, public :: modes_t
     integer :: nmodes
     character(len= 32),     pointer :: names(:) ! names used to identify a mode in the climate/diag lists
     character(len= 32),     pointer :: types(:) ! type of mode (as used in MAM code)
     type(mode_component_t), pointer :: comps(:) ! components which define the mode
  end type modes_t

!! \section arg_table_bin_component_t
!! \htmlinclude bin_component_t.html
  ! type to provide access to the components of a bin
  type, public :: bin_component_t
     integer :: nspec
     ! For "source" variables below, value is:
     ! 'N' if in pbuf (non-advected)
     ! 'A' if in state (advected)
     character(len=  1) :: source_num_a  ! source of interstitial number conc field
     character(len= 32) :: camname_num_a ! name registered in pbuf or constituents for number mixing ratio of interstitial species
     character(len=  1) :: source_num_c  ! source of cloud borne number conc field
     character(len= 32) :: camname_num_c ! name registered in pbuf or constituents for number mixing ratio of cloud borne species

     character(len=  1) :: source_mass_a  ! source of interstitial number conc field
     character(len= 32) :: camname_mass_a ! name registered in pbuf or constituents for number mixing ratio of interstitial species
     character(len=  1) :: source_mass_c  ! source of cloud borne number conc field
     character(len= 32) :: camname_mass_c ! name registered in pbuf or constituents for number mixing ratio of cloud borne species

     ! source of interstitial mmr field
     character(len=  1),         pointer :: source_mmr_a(:)
     ! name registered in pbuf or constituents for mmr species
     character(len= 32),         pointer :: camname_mmr_a(:)
     ! source of cloud borne specie mmr fields
     character(len=  1),         pointer :: source_mmr_c(:)
     ! name registered in pbuf or constituents for mmr of cloud borne components
     character(len= 32),         pointer :: camname_mmr_c(:)
     ! species type
     character(len= 32),         pointer :: type(:)
     ! species morphology
     character(len= 32),         pointer :: morph(:)
     ! file containing specie properties
     character(len=shr_kind_cl), pointer :: props(:)

     ! index in pbuf or constituents for number mixing ratio of interstitial species
     integer          :: idx_num_a
     ! index in pbuf for number mixing ratio of cloud-borne species
     integer          :: idx_num_c
     ! index in pbuf or constituents for mass mixing ratio of interstitial species
     integer          :: idx_mass_a
     ! index in pbuf for mass mixing ratio of cloud-borne species
     integer          :: idx_mass_c

     ! index in pbuf or constituents for mmr of interstitial species
     integer, pointer :: idx_mmr_a(:)
     ! index in pbuf or constituents for mmr of cloud-borne species
     integer, pointer :: idx_mmr_c(:)
     ! ID used to access physical properties of mode species from phys_prop module
     integer, pointer :: idx_props(:)
  end type bin_component_t

!! \section arg_table_bins_t
!! \htmlinclude bins_t.html
  ! type to provide access to all bins
  type, public :: bins_t
     integer :: nbins
     character(len= 32),    pointer :: names(:) ! names used to identify a mode in the climate/diag lists
     type(bin_component_t), pointer :: comps(:) ! components which define the mode
  end type bins_t

!! \section arg_table_aerosol_t
!! \htmlinclude aerosol_t.html
  ! Storage for bulk aerosol components in the climate/diagnostic lists
  type, public :: aerosol_t
     character(len=1)           :: source         ! A for state (advected), N for pbuf (non-advected), Z for zero
     character(len=64)          :: camname        ! name of constituent in physics state or buffer
     character(len=shr_kind_cl) :: physprop_file  ! physprop filename
     character(len=32)          :: mass_name      ! name for mass per layer field in history output
     integer                    :: idx            ! index of constituent in physics state or buffer
     integer                    :: physprop_id    ! ID used to access physical properties from phys_prop module
  end type aerosol_t

!! \section arg_table_aerlist_t
!! \htmlinclude aerlist_t.html
  type, public :: aerlist_t
     integer                  :: numaerosols  ! number of aerosols
     character(len=2)         :: list_id      ! set to "  " for climate list, or two character integer
                                              ! (include leading zero) to identify diagnostic list
     type(aerosol_t), pointer :: aer(:)       ! dimension(numaerosols)
  end type aerlist_t

!! \section arg_table_modelist_t
!! \htmlinclude modelist_t.html
  ! storage for modal aerosol components in the climate/diagnostic lists
  type, public :: modelist_t
     ! number of modes
     integer                             :: nmodes

     ! set to "  " for climate list, or two character integer
     ! (include leading zero) to identify diagnostic list
     ! used to construct history field names and descriptions
     character(len=2)                    :: list_id

     ! index of the mode in the mode definition object
     integer,                    pointer :: idx(:)
     ! physprop filename
     character(len=shr_kind_cl), pointer :: physprop_files(:)
     ! index of the mode properties in the physprop object
     integer,                    pointer :: idx_props(:)
  end type modelist_t

!! \section arg_table_binlist_t
!! \htmlinclude binlist_t.html
  ! storage for bin aerosol components in the climate/diagnostic lists
  type, public :: binlist_t
     ! number of bins
     integer            :: nbins

     ! set to "  " for climate list, or two character integer
     ! (include leading zero) to identify diagnostic list
     ! used to construct history field names and descriptions
     character(len=2)   :: list_id

     ! index of the bin in the bin definition object
     integer,   pointer :: idx(:)
     ! physprop filename
     character(len=shr_kind_cl), pointer :: physprop_files(:)
     ! index of the bin properties in the physprop object
     integer,   pointer :: idx_props(:)
  end type binlist_t

  ! max number of strings in mode definitions
  integer, public, parameter :: n_mode_str = 120

  ! max number of strings in bin definitions
  integer, public, parameter :: n_bin_str = 640

  !===========================
  ! Shared constants (shared with rad_constituents for gases)
  ! These have CCPP framework metadata attached to them as
  ! physics/chemistry CCPP schemes make use of these quantities.
  !===========================
!> \section arg_table_radiative_aerosol_definitions  Argument Table
!! \htmlinclude radiative_aerosol_definitions.html
  ! maximum number of diagnostic lists
  integer, public, parameter :: N_DIAG = 10

  ! max number of externally mixed entities in the climate/diag lists
  integer, public, parameter :: n_rad_cnst = 80

  ! climate list identifier (to keep CCPP framework happy)
  integer, public, parameter :: id_climate = 0

  !===========================
  ! Aerosol-specific module data.
  !===========================
  ! namelist data container per climate/diagnostic list.
  type(rad_cnst_namelist_t), public :: radcnst_namelist(id_climate:N_DIAG)

  ! flag for whether diagnostic lists are active
  logical, public :: active_calls(id_climate:N_DIAG) = .false.

  type(modes_t), public, target :: modes  ! mode definitions
  type(bins_t),  public, target :: bins   ! bin definitions

  ! list of bulk aerosols used in climate/diagnostic calculations
  type(aerlist_t),  public, target :: bulk_aerosol_list(id_climate:N_DIAG)

  ! list of aerosol modes used in climate/diagnostic calculations
  type(modelist_t), public, target :: modal_aerosol_list(id_climate:N_DIAG)

  ! list of aerosol bins used in climate/diagnostic calcs
  type(binlist_t),  public, target :: sectional_aerosol_list(id_climate:N_DIAG)

!==============================================================================
contains
!==============================================================================

subroutine list_populate(namelist, aerlist, modal_aerosol_list, sectional_aerosol_list)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use spmd_utils,     only: masterproc

   ! Populate aerosol list structures from parsed namelist specifiers.
   ! IMPORTANT: Must run at readnl time (before phys_register), because
   ! phys_register routines (e.g., modal_aero_data_reg) query
   ! modal_aerosol_list(0)%nmodes via rad_aer_get_info.
   ! Do NOT merge with list_resolve_physprops.
   !
   ! Gas initialization is handled in rad_constituents.
   type(rad_cnst_namelist_t), intent(in) :: namelist ! parsed namelist input for climate or diagnostic lists

   type(aerlist_t),        intent(inout) :: aerlist
   type(modelist_t),       intent(inout) :: modal_aerosol_list
   type(binlist_t),        intent(inout) :: sectional_aerosol_list

   ! Local variables
   integer :: ii, m, naero, nmodes, nbins
   integer :: ba_idx, ma_idx, sa_idx
   integer :: istat
   character(len=*), parameter :: subname = 'list_populate'
   !-----------------------------------------------------------------------------

   ! Determine the number of bulk aerosols and aerosol modes in the list
   naero = 0
   nmodes = 0
   nbins = 0
   do ii = 1, namelist%ncnst
      if (trim(namelist%type(ii)) == 'A') naero  = naero + 1
      if (trim(namelist%type(ii)) == 'M') nmodes = nmodes + 1
      if (trim(namelist%type(ii)) == 'B') nbins = nbins + 1
   end do
   aerlist%numaerosols = naero
   modal_aerosol_list%nmodes      = nmodes
   sectional_aerosol_list%nbins       = nbins

   ! allocate storage for the aerosol and mode lists
   allocate( &
      aerlist%aer(aerlist%numaerosols),      &
      modal_aerosol_list%idx(modal_aerosol_list%nmodes),           &
      modal_aerosol_list%physprop_files(modal_aerosol_list%nmodes), &
      modal_aerosol_list%idx_props(modal_aerosol_list%nmodes),     &
      sectional_aerosol_list%idx(sectional_aerosol_list%nbins),           &
      sectional_aerosol_list%physprop_files(sectional_aerosol_list%nbins), &
      sectional_aerosol_list%idx_props(sectional_aerosol_list%nbins),     &
      stat=istat)
   if (istat /= 0) call endrun(subname//': allocate ERROR; aero list components')

   if (masterproc .and. verbose) then
      if (len_trim(aerlist%list_id) == 0) then
         write(iulog,*) newline//' '//subname//': namelist input for climate list'
      else
         write(iulog,*) newline//' '//subname//': namelist input for diagnostic list:'//aerlist%list_id
      end if
   end if

   ! Loop over the radiatively active components specified in the namelist
   ba_idx = 0
   ma_idx = 0
   sa_idx = 0
   do ii = 1, namelist%ncnst

      ! Skip gas entries (handled in rad_constituents)
      if (namelist%type(ii) == 'G') cycle

      if (masterproc .and. verbose) &
         write(iulog,*) "  rad namelist spec: "// trim(namelist%source(ii)) &
         //":"//trim(namelist%camname(ii))//":"//trim(namelist%radname(ii))

      ! Check that the source specifier is legal.
      if (namelist%source(ii) /= 'A' .and. namelist%source(ii) /= 'M' .and. &
          namelist%source(ii) /= 'N' .and. namelist%source(ii) /= 'Z' .and. &
          namelist%source(ii) /= 'B' ) then
         call endrun(subname//": source must either be A, B, M, N or Z:"//&
                     " illegal specifier in namelist input: "//namelist%source(ii))
      end if

      ! Add component to appropriate list (modal or bulk aerosol)
      if (namelist%type(ii) == 'A') then

         ! Add to bulk aerosol list
         ba_idx = ba_idx + 1

         aerlist%aer(ba_idx)%source        = namelist%source(ii)
         aerlist%aer(ba_idx)%camname       = namelist%camname(ii)
         aerlist%aer(ba_idx)%physprop_file = namelist%radname(ii)

      else if (namelist%type(ii) == 'M') then

         ! Add to modal aerosol list
         ma_idx = ma_idx + 1

         ! Look through the mode definitions for the name of the specified mode.  The
         ! index into the modes object all the information relevent to the mode definition.
         modal_aerosol_list%idx(ma_idx) = -1
         do m = 1, modes%nmodes
            if (trim(namelist%camname(ii)) == trim(modes%names(m))) then
               modal_aerosol_list%idx(ma_idx) = m
               exit
            end if
         end do
         if (modal_aerosol_list%idx(ma_idx) == -1) &
            call endrun(subname//' ERROR cannot find mode name '//trim(namelist%camname(ii)))

         ! Also save the name of the physprop file
         modal_aerosol_list%physprop_files(ma_idx) = namelist%radname(ii)

      else if (namelist%type(ii) == 'B') then

         ! Add to bin aerosol list
         sa_idx = sa_idx + 1

         ! Look through the bin definitions for the name of the specified bin.  The
         ! index into the bins object all the information relevent to the bin definition.
         sectional_aerosol_list%idx(sa_idx) = -1
         do m = 1, bins%nbins
            if (trim(namelist%camname(ii)) == trim(bins%names(m))) then
               sectional_aerosol_list%idx(sa_idx) = m
               exit
            end if
         end do
         if (sectional_aerosol_list%idx(sa_idx) == -1) &
            call endrun(subname//' ERROR cannot find bin name '//trim(namelist%camname(ii)))

         ! Also save the name of the physprop file
         sectional_aerosol_list%physprop_files(sa_idx) = namelist%radname(ii)

      end if
   end do

end subroutine list_populate

!===========================

subroutine list_resolve_physprops(aerlist, modal_aerosol_list, sectional_aerosol_list)

   ! Resolve physprop indices for bulk aerosols, modes, and bins.
   ! IMPORTANT: Must run at init time (after physprop_init), because
   ! physprop_get_id requires physprop files to have been read.
   ! Do NOT merge with list_populate.
   !
   ! Host-specific index resolution (get_cam_idx) is handled
   ! separately by the host module (e.g. aerosol_mmr_cam).

   use phys_prop, only: physprop_get_id

   type(aerlist_t),        intent(inout) :: aerlist
   type(modelist_t),       intent(inout) :: modal_aerosol_list
   type(binlist_t),        intent(inout) :: sectional_aerosol_list

   ! Local variables
   integer :: i
   character(len=*), parameter :: subname = 'list_resolve_physprops'
   !-----------------------------------------------------------------------------

   ! Loop over bulk aerosols
   do i = 1, aerlist%numaerosols

      ! get the physprop_id from the phys_prop module
      aerlist%aer(i)%physprop_id = physprop_get_id(aerlist%aer(i)%physprop_file)

   end do

   ! Loop over modes
   do i = 1, modal_aerosol_list%nmodes

      ! get the physprop_id from the phys_prop module
      modal_aerosol_list%idx_props(i) = physprop_get_id(modal_aerosol_list%physprop_files(i))

   end do

   ! Loop over bins
   do i = 1, sectional_aerosol_list%nbins

      ! get the physprop_id from the phys_prop module
      sectional_aerosol_list%idx_props(i) = physprop_get_id(sectional_aerosol_list%physprop_files(i))

   end do

end subroutine list_resolve_physprops

!===========================

subroutine parse_mode_defs(nl_in, modes)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog

   ! Parse the mode definition specifiers.  The specifiers are of the form:
   !
   ! 'mode_name:mode_type:=',
   !  'source_num_a:camname_num_a:source_num_c:camname_num_c:num_mr:+',
   !  'source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file[:+]'[,]
   !  ['source_mmr_a:camname_mmr_a:source_mmr_c:camname_mmr_c:spec_type:prop_file][:+][']


   character(len=*), intent(inout) :: nl_in(:)    ! namelist input (blanks are removed on output)
   type(modes_t),    intent(inout) :: modes       ! structure containing parsed input

   ! Local variables
   integer :: m
   integer :: istat
   integer :: nmodes, nstr
   integer :: mbeg, mcur
   integer :: nspec, ispec
   integer :: strlen, iend, ipos
   logical :: num_mr_found
   character(len=*), parameter :: subname = 'parse_mode_defs'
   character(len=len(nl_in(1))) :: tmpstr
   character(len=1)  :: tmp_src_a
   character(len=32) :: tmp_name_a
   character(len=1)  :: tmp_src_c
   character(len=32) :: tmp_name_c
   character(len=32) :: tmp_type
   !-------------------------------------------------------------------------

   ! Determine number of modes defined by counting number of strings that are
   ! terminated by ':='
   ! (algorithm stops counting at first blank element).
   nmodes = 0
   nstr = 0
   do m = 1, n_mode_str

      if (len_trim(nl_in(m)) == 0) exit
      nstr = nstr + 1

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(nl_in(m))
      nl_in(m) = tmpstr
      do
         strlen = len_trim(nl_in(m))
         ipos = index(nl_in(m), ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = nl_in(m)(:ipos-1) // nl_in(m)(ipos+1:strlen)
         nl_in(m) = tmpstr
      end do
      ! count strings with ':=' terminator
      if (nl_in(m)(strlen-1:strlen) == ':=') nmodes = nmodes + 1

   end do
   modes%nmodes = nmodes

   ! return if no modes defined
   if (nmodes == 0) return

   ! allocate components that depend on nmodes
   allocate( &
      modes%names(nmodes),  &
      modes%types(nmodes),  &
      modes%comps(nmodes),  &
      stat=istat )
   if (istat > 0) then
      write(iulog,*) subname//': ERROR: cannot allocate storage for modes.  nmodes=', nmodes
      call endrun(subname//': ERROR allocating storage for modes')
   end if

   mcur = 1              ! index of current string being processed

   ! loop over modes
   do m = 1, nmodes

      mbeg = mcur  ! remember the first string of a mode

      ! check that first string in mode definition is ':=' terminated
      iend = len_trim(nl_in(mcur))
      if (nl_in(mcur)(iend-1:iend) /= ':=') call parse_error('= not found', nl_in(mcur))

      ! count species in mode definition.  definition will contain 1 string with
      ! with a ':+' terminator for each specie
      nspec = 0
      mcur = mcur + 1
      do
         iend = len_trim(nl_in(mcur))
         if (nl_in(mcur)(iend-1:iend) /= ':+') exit
         nspec = nspec + 1
         mcur = mcur + 1
      end do

      ! a mode must have at least one specie
      if (nspec == 0) call parse_error('mode must have at least one specie', nl_in(mbeg))

      ! allocate components that depend on number of species
      allocate( &
         modes%comps(m)%source_mmr_a(nspec),  &
         modes%comps(m)%camname_mmr_a(nspec), &
         modes%comps(m)%source_mmr_c(nspec),  &
         modes%comps(m)%camname_mmr_c(nspec), &
         modes%comps(m)%type(nspec),          &
         modes%comps(m)%props(nspec),         &
         stat=istat)

      if (istat > 0) then
         write(iulog,*) subname//': ERROR: cannot allocate storage for species.  nspec=', nspec
         call endrun(subname//': ERROR allocating storage for species')
      end if

      ! initialize components
      modes%comps(m)%nspec         = nspec
      modes%comps(m)%source_num_a  = ' '
      modes%comps(m)%camname_num_a = ' '
      modes%comps(m)%source_num_c  = ' '
      modes%comps(m)%camname_num_c = ' '
      do ispec = 1, nspec
         modes%comps(m)%source_mmr_a(ispec)  = ' '
         modes%comps(m)%camname_mmr_a(ispec) = ' '
         modes%comps(m)%source_mmr_c(ispec)  = ' '
         modes%comps(m)%camname_mmr_c(ispec) = ' '
         modes%comps(m)%type(ispec)          = ' '
         modes%comps(m)%props(ispec)         = ' '
      end do

      ! return to first string in mode definition
      mcur = mbeg
      tmpstr = nl_in(mcur)

      ! mode name
      ipos = index(tmpstr, ':')
      if (ipos < 2) call parse_error('mode name not found', tmpstr)
      modes%names(m) = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! mode type
      ipos = index(tmpstr, ':')
      if (ipos == 0) call parse_error('mode type not found', tmpstr)
      ! check for valid mode type
      call check_mode_type(tmpstr, 1, ipos-1)
      modes%types(m) = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! mode type must be followed by '='
      if (tmpstr(1:1) /= '=') call parse_error('= not found', tmpstr)

      ! move to next string
      mcur = mcur + 1
      tmpstr = nl_in(mcur)

      ! process mode component strings
      num_mr_found = .false.   ! keep track of whether number mixing ratio component is found
      ispec = 0                ! keep track of the number of species found
      do

         ! source of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find source field first', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_a = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_a = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! source of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find a source field', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_c = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_c = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! component type
         ipos = scan(tmpstr, ': ')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)

         if (tmpstr(:ipos-1) == 'num_mr') then

            ! there can only be one number mixing ratio component
            if (num_mr_found) call parse_error('more than 1 number component', nl_in(mcur))

            num_mr_found = .true.
            modes%comps(m)%source_num_a  = tmp_src_a
            modes%comps(m)%camname_num_a = tmp_name_a
            modes%comps(m)%source_num_c  = tmp_src_c
            modes%comps(m)%camname_num_c = tmp_name_c
            tmpstr                       = tmpstr(ipos+1:)

         else

            ! check for valid specie type
            call check_specie_type(tmpstr, 1, ipos-1)
            tmp_type = tmpstr(:ipos-1)
            tmpstr   = tmpstr(ipos+1:)

            ! get the properties file
            ipos = scan(tmpstr, ': ')
            if (ipos == 0) call parse_error('next separator not found', tmpstr)
            ! check for valid filename -- must have .nc extension
            if (tmpstr(ipos-3:ipos-1) /= '.nc') &
               call parse_error('filename not valid', tmpstr)

            ispec = ispec + 1
            modes%comps(m)%source_mmr_a(ispec)  = tmp_src_a
            modes%comps(m)%camname_mmr_a(ispec) = tmp_name_a
            modes%comps(m)%source_mmr_c(ispec)  = tmp_src_c
            modes%comps(m)%camname_mmr_c(ispec) = tmp_name_c
            modes%comps(m)%type(ispec)          = tmp_type
            modes%comps(m)%props(ispec)         = tmpstr(:ipos-1)
            tmpstr                              = tmpstr(ipos+1:)
         end if

         ! check if there are more components.  either the current character is
         ! a ' ' which means this string is the final mode component, or the character
         ! is a '+' which means there are more components
         if (tmpstr(1:1) == ' ') exit

         if (tmpstr(1:1) /= '+') &
               call parse_error('+ field not found', tmpstr)

         ! continue to next component...
         mcur = mcur + 1
         tmpstr = nl_in(mcur)
      end do

      ! check that a number component was found
      if (.not. num_mr_found) call parse_error('number component not found', nl_in(mbeg))

      ! check that the right number of species were found
      if (ispec /= nspec) call parse_error('component parsing got wrong number of species', nl_in(mbeg))

      ! continue to next mode...
      mcur = mcur + 1
      tmpstr = nl_in(mcur)
   end do

   !------------------------------------------------------------------------------------------------
   contains
   !------------------------------------------------------------------------------------------------

   subroutine parse_error(msg, str)

      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: str

      write(iulog,*) subname//': ERROR: '//msg
      write(iulog,*) ' input string: '//trim(str)
      call endrun(subname//': ERROR: '//msg)

   end subroutine parse_error

   !------------------------------------------------------------------------------------------------

   subroutine check_specie_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie

      integer :: i

      do i = 1, num_spec_types
         if (str(ib:ie) == trim(spec_type_names(i))) return
      end do

      call parse_error('specie type not valid', str(ib:ie))

   end subroutine check_specie_type

   !------------------------------------------------------------------------------------------------

   subroutine check_mode_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie  ! begin, end character of mode type substring

      integer :: i

      do i = 1, num_mode_types
         if (str(ib:ie) == trim(mode_type_names(i))) return
      end do

      call parse_error('mode type not valid', str(ib:ie))

   end subroutine check_mode_type

   !------------------------------------------------------------------------------------------------

end subroutine parse_mode_defs

!===========================

subroutine parse_bin_defs(nl_in, bins)
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog

   ! Parse the bin definition specifiers.

   character(len=*), intent(inout) :: nl_in(:)    ! namelist input (blanks are removed on output)
   type(bins_t),    intent(inout) :: bins       ! structure containing parsed input

   ! Local variables
   logical :: num_mr_found, mass_mr_found
   integer :: m
   integer :: istat
   integer :: nbins, nstr
   integer :: mbeg, mcur
   integer :: nspec, ispec
   integer :: strlen, iend, ipos
   character(len=*), parameter :: subname = 'parse_bin_defs'
   character(len=len(nl_in(1))) :: tmpstr
   character(len=1)  :: tmp_src_a
   character(len=32) :: tmp_name_a
   character(len=1)  :: tmp_src_c
   character(len=32) :: tmp_name_c
   character(len=32) :: tmp_type
   character(len=32) :: tmp_morph
   !-------------------------------------------------------------------------

   ! Determine number of bins defined by counting number of strings that are
   ! terminated by ':='
   ! (algorithm stops counting at first blank element).
   nbins = 0
   nstr = 0
   do m = 1, n_bin_str

      if (len_trim(nl_in(m)) == 0) exit
      nstr = nstr + 1

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(nl_in(m))
      nl_in(m) = tmpstr
      do
         strlen = len_trim(nl_in(m))
         ipos = index(nl_in(m), ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = nl_in(m)(:ipos-1) // nl_in(m)(ipos+1:strlen)
         nl_in(m) = tmpstr
      end do
      ! count strings with ':=' terminator
      if (nl_in(m)(strlen-1:strlen) == ':=') nbins = nbins + 1

   end do
   bins%nbins = nbins

   ! return if no bins defined
   if (nbins == 0) return

   ! allocate components that depend on nmodes
   allocate( &
      bins%names(nbins),  &
      bins%comps(nbins),  &
      stat=istat )
   if (istat > 0) then
      write(iulog,*) subname//': ERROR: cannot allocate storage for bins.  nbins=', nbins
      call endrun(subname//': ERROR allocating storage for bins')
   end if

   mcur = 1              ! index of current string being processed

   ! loop over bins
   bins_loop: do m = 1, nbins

      mbeg = mcur  ! remember the first string of a bin

      ! check that first string in bin definition is ':=' terminated
      iend = len_trim(nl_in(mcur))
      if (nl_in(mcur)(iend-1:iend) /= ':=') call parse_error('= not found', nl_in(mcur))

      ! count species in bin definition.  definition will contain 1 string with
      ! with a ':+' terminator for each specie
      nspec = 0
      mcur = mcur + 1
      do
         iend = len_trim(nl_in(mcur))
         if (nl_in(mcur)(iend-1:iend) /=    ':+') exit
         if (nl_in(mcur)(iend-4:iend) /= 'mmr:+') nspec = nspec + 1
         mcur = mcur + 1
      end do

      ! a bin must have at least one specie
      if (nspec == 0) call parse_error('bin must have at least one specie', nl_in(mbeg))

      ! allocate components that depend on number of species
      allocate( &
         bins%comps(m)%source_mmr_a(nspec),  &
         bins%comps(m)%camname_mmr_a(nspec), &
         bins%comps(m)%source_mmr_c(nspec),  &
         bins%comps(m)%camname_mmr_c(nspec), &
         bins%comps(m)%type(nspec),          &
         bins%comps(m)%morph(nspec),          &
         bins%comps(m)%props(nspec),         &
         stat=istat)

      if (istat > 0) then
         write(iulog,*) subname//': ERROR: cannot allocate storage for species.  nspec=', nspec
         call endrun(subname//': ERROR allocating storage for species')
      end if

      ! initialize components
      bins%comps(m)%nspec         = nspec
      bins%comps(m)%source_num_a  = ' '
      bins%comps(m)%camname_num_a = ' '
      bins%comps(m)%source_num_c  = ' '
      bins%comps(m)%camname_num_c = ' '
      bins%comps(m)%source_mass_a  = 'NOTSET'
      bins%comps(m)%camname_mass_a = 'NOTSET'
      bins%comps(m)%source_mass_c  = 'NOTSET'
      bins%comps(m)%camname_mass_c = 'NOTSET'
      do ispec = 1, nspec
         bins%comps(m)%source_mmr_a(ispec)  = ' '
         bins%comps(m)%camname_mmr_a(ispec) = ' '
         bins%comps(m)%source_mmr_c(ispec)  = ' '
         bins%comps(m)%camname_mmr_c(ispec) = ' '
         bins%comps(m)%type(ispec)          = ' '
         bins%comps(m)%props(ispec)         = ' '
      end do

      ! return to first string in mode definition
      mcur = mbeg
      tmpstr = nl_in(mcur)

      ! bin name
      ipos = index(tmpstr, ':')
      if (ipos < 2) call parse_error('bin name not found', tmpstr)
      bins%names(m)  = tmpstr(:ipos-1)
      tmpstr         = tmpstr(ipos+1:)

      ! bin name must be followed by '='
      if (tmpstr(1:1) /= '=') call parse_error('= not found', tmpstr)

      ! move to next string
      mcur = mcur + 1
      tmpstr = nl_in(mcur)

      ! process bin component strings
      num_mr_found = .false.        ! keep track of whether number mixing ratio component is found
      mass_mr_found = .false.       ! keep track of whether number mixing ratio component is found
      ispec = 0                     ! keep track of the number of species found
      comps_loop: do

         ! source of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find source field first', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_a = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of interstitial component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_a = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! source of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos < 2) call parse_error('expect to find a source field', tmpstr)
         ! check for valid source
         if (tmpstr(:ipos-1) /= 'A' .and. tmpstr(:ipos-1) /= 'N' .and. tmpstr(:ipos-1) /= 'Z') &
            call parse_error('source must be A, N or Z', tmpstr)
         tmp_src_c = tmpstr(:ipos-1)
         tmpstr    = tmpstr(ipos+1:)

         ! name of cloud borne component
         ipos = index(tmpstr, ':')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)
         tmp_name_c = tmpstr(:ipos-1)
         tmpstr     = tmpstr(ipos+1:)

         ! component type
         ipos = scan(tmpstr, ': ')
         if (ipos == 0) call parse_error('next separator not found', tmpstr)

         if (tmpstr(:ipos-1) == 'num') then

            ! there can only be one number mixing ratio component
            if (num_mr_found) call parse_error('more than 1 number component', nl_in(mcur))

            num_mr_found = .true.
            bins%comps(m)%source_num_a  = tmp_src_a
            bins%comps(m)%camname_num_a = tmp_name_a
            bins%comps(m)%source_num_c  = tmp_src_c
            bins%comps(m)%camname_num_c = tmp_name_c
            tmpstr                      = tmpstr(ipos+1:)

         else if (tmpstr(:ipos-1) == 'mmr') then

            ! there can only be one number mixing ratio component
            if (mass_mr_found) call parse_error('more than 1 mass mixing ratio component', nl_in(mcur))

            mass_mr_found = .true.
            bins%comps(m)%source_mass_a  = tmp_src_a
            bins%comps(m)%camname_mass_a = tmp_name_a
            bins%comps(m)%source_mass_c  = tmp_src_c
            bins%comps(m)%camname_mass_c = tmp_name_c
            tmpstr                       = tmpstr(ipos+1:)

         else

            ! check for valid species type
            call check_bin_type(tmpstr, 1, ipos-1)
            tmp_type = tmpstr(:ipos-1)
            tmpstr   = tmpstr(ipos+1:)

            ipos = index(tmpstr, ':')
            if (ipos == 0) call parse_error('next separator not found', tmpstr)

            ! check for valid species type
            call check_bin_morph(tmpstr, 1, ipos-1)
            tmp_morph = tmpstr(:ipos-1)
            tmpstr   = tmpstr(ipos+1:)

            ! get the properties file
            ipos = scan(tmpstr, ': ')
            if (ipos == 0) call parse_error('next separator not found', tmpstr)

             ! check for valid filename -- must have .nc extension
            if (tmpstr(ipos-3:ipos-1) /= '.nc') &
               call parse_error('filename not valid', tmpstr)

            ispec = ispec + 1

            bins%comps(m)%source_mmr_a(ispec)  = tmp_src_a
            bins%comps(m)%camname_mmr_a(ispec) = tmp_name_a
            bins%comps(m)%source_mmr_c(ispec)  = tmp_src_c
            bins%comps(m)%camname_mmr_c(ispec) = tmp_name_c
            bins%comps(m)%type(ispec)          = tmp_type
            bins%comps(m)%morph(ispec)         = tmp_morph

            bins%comps(m)%props(ispec)         = tmpstr(:ipos-1)
            tmpstr                             = tmpstr(ipos+1:)

         endif

         ! check if there are more components.  either the current character is
         ! a ' ' which means this string is the final mode component, or the character
         ! is a '+' which means there are more components
         if (tmpstr(1:1) == ' ') then
            exit comps_loop
         endif

         if (tmpstr(1:1) /= '+') &
               call parse_error('+ field not found', tmpstr)

         ! continue to next component...
         mcur = mcur + 1
         tmpstr = nl_in(mcur)
      end do comps_loop


      ! check that a number component was found
      if (.not. num_mr_found) call parse_error('number component not found', nl_in(mbeg))

      ! check that the right number of species were found
      if (ispec /= nspec) then
         write(*,*) 'ispec, nspec = ',ispec, nspec
         call parse_error('component parsing got wrong number of species', nl_in(mbeg))
      endif

      ! continue to next bin...
      mcur = mcur + 1
      tmpstr = nl_in(mcur)
   end do bins_loop

   !------------------------------------------------------------------------------------------------
   contains
   !------------------------------------------------------------------------------------------------

   subroutine parse_error(msg, str)

      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: str

      write(iulog,*) subname//': ERROR: '//msg
      write(iulog,*) ' input string: '//trim(str)
      call endrun(subname//': ERROR: '//msg)

   end subroutine parse_error

   !------------------------------------------------------------------------------------------------

   subroutine check_bin_morph(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie

      integer :: i

      do i = 1, num_bin_morphs
         if (str(ib:ie) == trim(bin_morph_names(i))) return
      end do

      call parse_error('bin morph not valid', str(ib:ie))

   end subroutine check_bin_morph

   !------------------------------------------------------------------------------------------------
   subroutine check_bin_type(str, ib, ie)

      character(len=*), intent(in) :: str
      integer,          intent(in) :: ib, ie  ! begin, end character of mode type substring

      integer :: i

      do i = 1, num_spec_types
         if (str(ib:ie) == trim(spec_type_names(i))) return
      end do

      call parse_error('bin species type not valid', str(ib:ie))

   end subroutine check_bin_type

   !------------------------------------------------------------------------------------------------

end subroutine parse_bin_defs

!===========================

subroutine parse_rad_specifier(specifier, namelist_data)
    use cam_abortutils, only: endrun

!-----------------------------------------------------------------------------
! Parse the radiation namelist specifiers.
!-----------------------------------------------------------------------------

    character(len=*), dimension(:), intent(in) :: specifier
    type(rad_cnst_namelist_t),   intent(inout) :: namelist_data

    ! Local variables
    integer                    :: number, i, j
    integer                    :: ipos, strlen
    integer                    :: astat
    character(len=shr_kind_cl) :: tmpstr
    character(len=1)           :: source(n_rad_cnst)
    character(len=64)          :: camname(n_rad_cnst)
    character(len=shr_kind_cl) :: radname(n_rad_cnst)
    character(len=1)           :: type(n_rad_cnst)
    !-------------------------------------------------------------------------

    number = 0

    parse_loop: do i = 1, n_rad_cnst
      if ( len_trim(specifier(i)) == 0 ) then
         exit parse_loop
      endif

      ! There are no fields in the input strings in which a blank character is allowed.
      ! To simplify the parsing go through the input strings and remove blanks.
      tmpstr = adjustl(specifier(i))
      do
         strlen = len_trim(tmpstr)
         ipos = index(tmpstr, ' ')
         if (ipos == 0 .or. ipos > strlen) exit
         tmpstr = tmpstr(:ipos-1) // tmpstr(ipos+1:strlen)
      end do

      ! Locate the ':' separating source from camname.
      j = index(tmpstr, ':')
      source(i) = tmpstr(:j-1)
      tmpstr = tmpstr(j+1:)

      ! locate the ':' separating camname from radname
      j = scan(tmpstr, ':')

      camname(i) = tmpstr(:j-1)
      radname(i) = tmpstr(j+1:)

      ! determine the type of constituent
      if (source(i) == 'M') then
         type(i) = 'M'
      else if (source(i) == 'B') then
         type(i) = 'B'
      else if(index(radname(i),".nc") .gt. 0) then
         type(i) = 'A'
      else
         type(i) = 'G'
      end if

      number = number+1
    end do parse_loop

    namelist_data%ncnst = number

    if (number == 0) return

    allocate(namelist_data%source (number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%source')
    allocate(namelist_data%camname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%camname')
    allocate(namelist_data%radname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%radname')
    allocate(namelist_data%type(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%type')

    namelist_data%source(:namelist_data%ncnst)  = source (:namelist_data%ncnst)
    namelist_data%camname(:namelist_data%ncnst) = camname(:namelist_data%ncnst)
    namelist_data%radname(:namelist_data%ncnst) = radname(:namelist_data%ncnst)
    namelist_data%type(:namelist_data%ncnst)    = type(:namelist_data%ncnst)

end subroutine parse_rad_specifier

!===========================

subroutine print_modes(modes)
   use cam_logfile,    only: iulog

   type(modes_t), intent(inout) :: modes

   integer :: i, m
   !---------------------------------------------------------------------------------------------

   write(iulog,*)' Mode Definitions'

   do m = 1, modes%nmodes

      write(iulog,*) newline//' name=',trim(modes%names(m)),'  type=',trim(modes%types(m))
      write(iulog,*) ' src_a=',trim(modes%comps(m)%source_num_a),'  num_a=',trim(modes%comps(m)%camname_num_a), &
                     ' src_c=',trim(modes%comps(m)%source_num_c),'  num_c=',trim(modes%comps(m)%camname_num_c)

      do i = 1, modes%comps(m)%nspec

         write(iulog,*) ' src_a=',trim(modes%comps(m)%source_mmr_a(i)), '  mmr_a=',trim(modes%comps(m)%camname_mmr_a(i)), &
                       '  src_c=',trim(modes%comps(m)%source_mmr_c(i)), '  mmr_c=',trim(modes%comps(m)%camname_mmr_c(i)), &
                       '  type=',trim(modes%comps(m)%type(i))
         write(iulog,*) '     prop file=', trim(modes%comps(m)%props(i))
      end do

   end do

end subroutine print_modes

!===========================

subroutine print_bins(bins)
   use cam_logfile,    only: iulog

   type(bins_t), intent(inout) :: bins

   integer :: i, m
   !---------------------------------------------------------------------------------------------

   write(iulog,*)' Bin Definitions'

   do m = 1, bins%nbins

      write(iulog,*) newline//' name=',trim(bins%names(m))

      do i = 1, bins%comps(m)%nspec

         write(iulog,*) ' src_a=',trim(bins%comps(m)%source_mmr_a(i)), '  mmr_a=',trim(bins%comps(m)%camname_mmr_a(i)), &
                       '  type=',trim(bins%comps(m)%type(i))
         write(iulog,*) '     prop file=', trim(bins%comps(m)%props(i))
      end do

   end do

end subroutine print_bins

!===========================

end module radiative_aerosol_definitions
