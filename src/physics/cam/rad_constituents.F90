module rad_constituents

!------------------------------------------------------------------------------------------------
!
! Gas-only radiative constituent handling and cloud optics settings.
!
! Provides: namelist I/O (shared gas+aerosol namelist), gas list init,
! gas MMR retrieval (state/pbuf), gas diagnostics output, and cloud optics
! public variables (iceopticsfile, liqopticsfile, etc.).
!
! Aerosol handling is in radiative_aerosol (facade) backed by
! radiative_aerosol_definitions (core definitions) and aerosol_mmr_cam
! (CAM-specific MMR retrieval).
!
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: rga
use physics_types,  only: physics_state
use phys_control,   only: use_simple_phys
use radconstants,   only: nradgas, rad_gas_index
use cam_history,    only: addfld, fieldname_len, outfld, horiz_only
use physics_buffer, only: physics_buffer_desc, pbuf_get_field

use cam_abortutils, only: endrun
use cam_logfile,    only: iulog

! Import from radiative_aerosol_definitions (core definitions)
use radiative_aerosol_definitions, only: cs1, N_DIAG, n_rad_cnst, verbose, nl, &
                            rad_cnst_namelist_t, radcnst_namelist, active_calls, &
                            n_mode_str, n_bin_str, parse_rad_specifier

!REMOVECAM
use aerosol_mmr_cam, only: get_cam_idx
!REMOVECAM_END

use radiative_aerosol, only: rad_aer_readnl

implicit none
private
save

! Storage for gas components in the climate/diagnostic lists

type :: gas_t
   character(len=1)  :: source       ! A for state (advected), N for pbuf (non-advected), Z for zero
   character(len=64) :: camname      ! name of constituent in physics state or buffer
   character(len=32) :: mass_name    ! name for mass per layer field in history output
   integer           :: idx          ! index from constituents or from pbuf
end type gas_t

type :: gaslist_t
   integer                :: ngas
   character(len=2)       :: list_id  ! set to "  " for climate list, or two character integer
                                      ! (include leading zero) to identify diagnostic list
   type(gas_t), pointer   :: gas(:)   ! dimension(ngas) where ngas = nradgas is from radconstants
end type gaslist_t

type(gaslist_t), target :: gaslist(0:N_DIAG)  ! gasses used in climate/diagnostic calculations

! values for constituents with requested value of zero
real(r8), allocatable, target :: zero_cols(:,:)

! Public interfaces — routines in this module
public :: &
   rad_cnst_readnl,             &! read namelist values and parse
   rad_cnst_init,               &! find optics files and all constituents
   rad_cnst_get_info,           &! gas+aerosol info wrapper
   rad_cnst_get_gas,            &! return pointer to mmr for gasses
   rad_cnst_out                  ! output constituent diagnostics (mass per layer and column burden)

character(len=cs1), public :: iceopticsfile, liqopticsfile
character(len=32),  public :: icecldoptics,liqcldoptics
logical,            public :: oldcldoptics = .false.

! Namelist variables
character(len=cs1), dimension(n_mode_str) :: mode_defs   = ' '
character(len=cs1), dimension(n_bin_str) :: bin_defs   = ' '
character(len=cs1) :: rad_climate(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_1(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_2(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_3(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_4(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_5(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_6(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_7(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_8(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_9(n_rad_cnst) = ' '
character(len=cs1) :: rad_diag_10(n_rad_cnst) = ' '

!==============================================================================
contains

subroutine rad_cnst_readnl(nlfile)

   ! Read rad_cnst_nl namelist group.  Parse input.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=2) :: suffix
   character(len=*), parameter :: subname = 'rad_cnst_readnl'

   namelist /rad_cnst_nl/ mode_defs,     &
                          bin_defs,      &
                          rad_climate,   &
                          rad_diag_1,    &
                          rad_diag_2,    &
                          rad_diag_3,    &
                          rad_diag_4,    &
                          rad_diag_5,    &
                          rad_diag_6,    &
                          rad_diag_7,    &
                          rad_diag_8,    &
                          rad_diag_9,    &
                          rad_diag_10,   &
                          iceopticsfile, &
                          liqopticsfile, &
                          icecldoptics,  &
                          liqcldoptics,  &
                          oldcldoptics

   !-----------------------------------------------------------------------------

   if (use_simple_phys) return

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rad_cnst_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rad_cnst_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast (mode_defs,     len(mode_defs(1))*n_mode_str,     mpichar, 0, mpicom)
   call mpibcast (bin_defs,      len(bin_defs(1))*n_bin_str,       mpichar, 0, mpicom)
   call mpibcast (rad_climate,   len(rad_climate(1))*n_rad_cnst,   mpichar, 0, mpicom)
   call mpibcast (rad_diag_1,    len(rad_diag_1(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_2,    len(rad_diag_2(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_3,    len(rad_diag_3(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_4,    len(rad_diag_4(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_5,    len(rad_diag_5(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_6,    len(rad_diag_6(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_7,    len(rad_diag_7(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_8,    len(rad_diag_8(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_9,    len(rad_diag_9(1))*n_rad_cnst,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_10,   len(rad_diag_10(1))*n_rad_cnst,   mpichar, 0, mpicom)
   call mpibcast (iceopticsfile, len(iceopticsfile),               mpichar, 0, mpicom)
   call mpibcast (liqopticsfile, len(liqopticsfile),               mpichar, 0, mpicom)
   call mpibcast (liqcldoptics,  len(liqcldoptics),                mpichar, 0, mpicom)
   call mpibcast (icecldoptics,  len(icecldoptics),                mpichar, 0, mpicom)
   call mpibcast (oldcldoptics,  1,                                mpilog , 0, mpicom)
#endif

   ! Parse the namelist input strings

   ! Lists of externally mixed entities for climate and diagnostic calculations
   do i = 0,N_DIAG
      select case (i)
      case(0)
         call parse_rad_specifier(rad_climate, radcnst_namelist(i))
      case (1)
         call parse_rad_specifier(rad_diag_1, radcnst_namelist(i))
      case (2)
         call parse_rad_specifier(rad_diag_2, radcnst_namelist(i))
      case (3)
         call parse_rad_specifier(rad_diag_3, radcnst_namelist(i))
      case (4)
         call parse_rad_specifier(rad_diag_4, radcnst_namelist(i))
      case (5)
         call parse_rad_specifier(rad_diag_5, radcnst_namelist(i))
      case (6)
         call parse_rad_specifier(rad_diag_6, radcnst_namelist(i))
      case (7)
         call parse_rad_specifier(rad_diag_7, radcnst_namelist(i))
      case (8)
         call parse_rad_specifier(rad_diag_8, radcnst_namelist(i))
      case (9)
         call parse_rad_specifier(rad_diag_9, radcnst_namelist(i))
      case (10)
         call parse_rad_specifier(rad_diag_10, radcnst_namelist(i))
      end select
   enddo

   ! were there any constituents specified for the nth diagnostic call?
   ! if so, radiation will make a call with those consituents
   active_calls(:) = (radcnst_namelist(:)%ncnst > 0)

   ! Aerosol init phase 1: parse mode/bin defs, accumulate physprop files,
   ! set aerosol list_id fields, initialize aerosol lists
   call rad_aer_readnl(mode_defs, bin_defs)

   ! Gas init phase 1: set gas list_id fields and populate gas lists
   if (masterproc) write(iulog,*) nl//subname//': Radiation gas constituent lists:'
   do i = 0, N_DIAG
      if (active_calls(i)) then
         if (i > 0) then
            write(suffix, fmt = '(i2.2)') i
         else
            suffix='  '
         end if
         gaslist(i)%list_id = suffix

         call gas_list_populate(radcnst_namelist(i), gaslist(i))

         if (masterproc .and. verbose) then
            call print_gas_list(gaslist(i))
         end if
      end if
   end do

end subroutine rad_cnst_readnl

subroutine rad_cnst_init()

   ! The initialization of the gas and aerosol lists is finished by
   ! 1) read the physprop files
   ! 2) find the index of each constituent in the constituent or physics buffer arrays
   ! 3) find the index of the aerosol constituents used to access its properties from the
   !    physprop module.

   integer :: i
   logical, parameter :: stricttest = .true.
   character(len=*), parameter :: subname = 'rad_cnst_init'
   !-----------------------------------------------------------------------------

   ! memory to point to if zero value requested
   allocate(zero_cols(pcols,pver))
   zero_cols = 0._r8

   ! Resolve constituent indices for gas lists
   if (masterproc) write(iulog,*) nl//subname//': checking for radiative gas constituents'
   do i = 0, N_DIAG
      if (active_calls(i)) then
         call gas_list_resolve_cnst_idx(gaslist(i))
      end if
   end do

   ! Check that all gases supported by the radiative transfer code have been specified.
   if (stricttest) then
      do i = 1, nradgas
         if (gaslist(0)%gas(i)%source .eq. 'Z' ) then
            call endrun(subname//': list of radiative gasses must include all radiation gasses for the climate specication')
         endif
      enddo
   endif

   ! Initialize gas history output for climate diagnostic quantities
   call rad_gas_diag_init(gaslist(0))

end subroutine rad_cnst_init

subroutine gas_list_populate(namelist, gaslist)

   ! Populate gas list from parsed namelist data.
   ! Must run at readnl time for consistency with aerosol list_populate.
   ! Do NOT merge with gas_list_resolve_cnst_idx.

   type(rad_cnst_namelist_t), intent(in)    :: namelist
   type(gaslist_t),           intent(inout) :: gaslist

   ! Local variables
   integer :: ii, igas, istat
   character(len=*), parameter :: routine = 'gas_list_populate'
   !-----------------------------------------------------------------------------

   ! nradgas is set by the radiative transfer code
   gaslist%ngas = nradgas

   allocate(gaslist%gas(gaslist%ngas), stat=istat)
   if (istat /= 0) call endrun(routine//': allocate ERROR; gas list components')

   ! Initialize sources to zero (default for unspecified gases)
   do igas = 1, gaslist%ngas
      gaslist%gas(igas)%source  = 'Z'
      gaslist%gas(igas)%camname = ' '
   end do

   ! Populate gas entries from 'G' type namelist entries
   do ii = 1, namelist%ncnst
      if (namelist%type(ii) /= 'G') cycle

      if (masterproc .and. verbose) &
         write(iulog,*) "  rad namelist spec: "// trim(namelist%source(ii)) &
         //":"//trim(namelist%camname(ii))//":"//trim(namelist%radname(ii))

      ! rad_gas_index will abort on illegal names
      igas = rad_gas_index(namelist%radname(ii))

      gaslist%gas(igas)%source  = namelist%source(ii)
      gaslist%gas(igas)%camname = namelist%camname(ii)
   end do

end subroutine gas_list_populate

!================================================================================================

subroutine gas_list_resolve_cnst_idx(gaslist)

   ! Resolve constituent indices for gas list entries.
   ! Must run at init time (after constituent registration).
   ! Do NOT merge with gas_list_populate.

   type(gaslist_t), intent(inout) :: gaslist

   ! Local variables
   integer :: i
   character(len=*), parameter :: routine = 'gas_list_resolve_cnst_idx'
   !-----------------------------------------------------------------------------

   do i = 1, gaslist%ngas
      gaslist%gas(i)%idx = get_cam_idx(gaslist%gas(i)%source, gaslist%gas(i)%camname, routine)
   end do

end subroutine gas_list_resolve_cnst_idx

!================================================================================================

subroutine rad_cnst_get_info(list_idx, gasnames, use_data_o3, ngas)

   ! Gas variant of rad_cnst_get_info;
   ! aerosol information moved to radiative_aerosol::rad_aer_get_info.

   ! Arguments
   integer,                     intent(in)  :: list_idx
   character(len=64), optional, intent(out) :: gasnames(:)
   logical,           optional, intent(out) :: use_data_o3
   integer,           optional, intent(out) :: ngas

   ! Local variables
   type(gaslist_t),  pointer :: g_list
   integer          :: i, igas, gaslen
   character(len=1) :: source
   character(len=*), parameter :: subname = 'rad_cnst_get_info'
   !-----------------------------------------------------------------------------

   ! Handle gas arguments locally
   g_list => gaslist(list_idx)

   if (present(ngas)) then
      ngas = g_list%ngas
   endif

   if (present(gasnames)) then
      gaslen = size(gasnames)
      if (gaslen < g_list%ngas) then
         write(iulog,*) subname//': ERROR: ngas=', g_list%ngas, '  gaslen=', gaslen
         call endrun(subname//': ERROR: gasnames too short')
      end if
      do i = 1, g_list%ngas
         gasnames(i) = g_list%gas(i)%camname
      end do
   end if

   if (present(use_data_o3)) then
      igas = rad_gas_index('O3')
      source = g_list%gas(igas)%source
      use_data_o3 = .false.
      if (source == 'N') use_data_o3 = .true.
   endif

end subroutine rad_cnst_get_info

!================================================================================================

subroutine rad_cnst_get_gas(list_idx, gasname, state, pbuf, mmr)

   ! Return pointer to mass mixing ratio for the gas from the specified
   ! climate or diagnostic list.

   ! Arguments
   integer,                     intent(in) :: list_idx    ! index of the climate or a diagnostic list
   character(len=*),            intent(in) :: gasname
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   real(r8),                    pointer    :: mmr(:,:)

   ! Local variables
   integer :: lchnk
   integer :: igas
   integer :: idx
   character(len=1) :: source
   type(gaslist_t), pointer :: list
   character(len=*), parameter :: subname = 'rad_cnst_get_gas'
   !-----------------------------------------------------------------------------

   if (list_idx >= 0 .and. list_idx <= N_DIAG) then
      list => gaslist(list_idx)
   else
      write(iulog,*) subname//': list_idx =', list_idx
      call endrun(subname//': list_idx out of bounds')
   endif

   lchnk = state%lchnk

   ! Get index of gas in internal arrays.  rad_gas_index will abort if the
   ! specified gasname is not recognized by the radiative transfer code.
   igas = rad_gas_index(trim(gasname))

   ! Get data source
   source = list%gas(igas)%source
   idx    = list%gas(igas)%idx
   select case( source )
   case ('A')
      mmr => state%q(:,:,idx)
   case ('N')
      call pbuf_get_field(pbuf, idx, mmr)
   case ('Z')
      mmr => zero_cols
   end select

end subroutine rad_cnst_get_gas

subroutine rad_cnst_out(list_idx, state, pbuf)

   ! Output the mass per layer, and total column burdens for gas
   ! constituents in either the climate or diagnostic lists.
   ! Aerosol output is now handled by rad_aer_diag_out in aerosol_mmr_cam.

   ! Arguments
   integer,                     intent(in) :: list_idx
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)

   ! Local variables
   integer :: i, ngas, lchnk, ncol
   integer :: idx
   character(len=1)  :: source
   character(len=32) :: name, cbname
   real(r8)          :: mass(pcols,pver)
   real(r8)          :: cb(pcols)
   real(r8), pointer :: mmr(:,:)
   type(gaslist_t), pointer :: g_list
   character(len=*), parameter :: subname = 'rad_cnst_out'
   !-----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointer with requested gas list
   g_list => gaslist(list_idx)

   ngas = g_list%ngas
   do i = 1, ngas

      source = g_list%gas(i)%source
      idx    = g_list%gas(i)%idx
      name   = g_list%gas(i)%mass_name
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

end subroutine rad_cnst_out

subroutine rad_gas_diag_init(glist)

! Add diagnostic fields to the master fieldlist.

   type(gaslist_t), intent(inout) :: glist

   integer :: i, ngas
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   ngas = glist%ngas
   if (ngas == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = glist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, ngas

      ! construct names for mass per layer diagnostics
      name = 'm_' // trim(glist%gas(i)%camname) // trim(suffix)
      glist%gas(i)%mass_name = name
      long_name = trim(glist%gas(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), (/ 'lev' /), 'A', 'kg/m^2', trim(long_name))

      ! construct names for column burden diagnostics
      name = 'cb_' // trim(glist%gas(i)%camname) // trim(suffix)
      long_name = trim(glist%gas(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), horiz_only, 'A', 'kg/m^2', trim(long_name))

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_gas_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_gas_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_gas_diag_init

!================================================================================================

subroutine print_gas_list(glist)

   ! Print summary of gas list.

   use radconstants, only: gascnst=>gaslist

   type(gaslist_t),  intent(in) :: glist

   integer :: i

   if (len_trim(glist%list_id) == 0) then
      write(iulog,*) nl//' gas list for climate calculations'
   else
      write(iulog,*) nl//' gas list for diag'//glist%list_id//' calculations'
   end if

   do i = 1, nradgas
      if (glist%gas(i)%source .eq. 'N') then
         write(iulog,*) '  '//glist%gas(i)%source//':'//gascnst(i)//' has pbuf name:'//&
                        trim(glist%gas(i)%camname)
      else if (glist%gas(i)%source .eq. 'A') then
         write(iulog,*) '  '//glist%gas(i)%source//':'//gascnst(i)//' has constituents name:'//&
                        trim(glist%gas(i)%camname)
      endif
   enddo

end subroutine print_gas_list

!================================================================================================


end module rad_constituents
