!--------------------------------------------------------------------------------
! For bulk aerosol representation.
! Here each aerosol is treated as a separate bin.
!--------------------------------------------------------------------------------
module bulk_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use string_utils, only : to_lower

  use aerosol_properties_mod, only: aerosol_properties

  use radiative_aerosol, only: rad_aer_get_info, rad_aer_get_props
  use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)

  implicit none

  private

  public :: bulk_aerosol_properties

  type, extends(aerosol_properties) :: bulk_aerosol_properties

     private

   contains

     procedure :: number_transported
     procedure :: get
     procedure :: amcube
     procedure :: actfracs
     procedure :: num_names
     procedure :: mmr_names
     procedure :: amb_num_name
     procedure :: amb_mmr_name
     procedure :: species_type
     procedure :: icenuc_updates_num
     procedure :: icenuc_updates_mmr
     procedure :: apply_number_limits
     procedure :: hetfrz_species
     procedure :: physprop_id
     procedure :: soluble
     procedure :: min_mass_mean_rad
     procedure :: bin_name
     procedure :: scav_diam
     procedure :: resuspension_resize
     procedure :: rebin_bulk_fluxes
     procedure :: hydrophilic
     procedure :: model_is

     final :: destructor

  end type bulk_aerosol_properties

  interface bulk_aerosol_properties
     procedure :: constructor
  end interface bulk_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(list_idx) result(newobj)

    integer, optional, intent(in) :: list_idx ! radiation list index (0=climate)
    type(bulk_aerosol_properties), pointer :: newobj

    integer,allocatable :: nspecies(:)
    real(r8),allocatable :: alogsig(:)
    real(r8),allocatable :: f1(:)
    integer :: ierr, naero, i
    integer :: list_idx_loc
    real(r8) :: dispersion_val

    list_idx_loc = 0
    if (present(list_idx)) list_idx_loc = list_idx

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    call rad_aer_get_info(list_idx_loc, naero=naero)

    ! Here treat each aerosol as a separate bin
    allocate( nspecies(naero),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( alogsig(naero),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( f1(naero),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    ! Bulk aerosols have 1 chemical species in each bin
    nspecies(:) = 1

    ! Read actual dispersion (sigma_logr) from physprop files
    do i = 1, naero
       call rad_aer_get_props(list_idx_loc, i, dispersion_aer=dispersion_val)
       alogsig(i) = log(dispersion_val)
    end do
    f1(:) = 1._r8

    ! For bulk aerosols, the number of bins and total number of constituents are
    ! the same (naero) -- one constituent (species and mass) per bin.
    call newobj%initialize(nbin=naero, ncnst=naero, nspec=nspecies, nmasses=nspecies, &
                           alogsig=alogsig, f1=f1, f2=f1, ierr=ierr, list_idx=list_idx_loc)

    deallocate(nspecies)
    deallocate(alogsig)
    deallocate(f1)

    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(bulk_aerosol_properties), intent(inout) :: self

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! returns number of transported aerosol constituents
  !------------------------------------------------------------------------------
  integer function number_transported(self)
    class(bulk_aerosol_properties), intent(in) :: self
    ! to be implemented later
    number_transported = -1
  end function number_transported


  !------------------------------------------------------------------------
  ! returns aerosol properties:
  !  density
  !  hygroscopicity
  !  species type
  !  species name
  !  short wave species refractive indices
  !  long wave species refractive indices
  !  species morphology
  !------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, density, hygro, &
                 spectype, specname, specmorph, refindex_sw, refindex_lw, num_to_mass_aer, &
                 dryrad)

    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    real(r8), optional, intent(out) :: density ! density (kg/m3)
    real(r8), optional, intent(out) :: hygro   ! hygroscopicity
    character(len=*), optional, intent(out) :: spectype  ! species type
    character(len=*), optional, intent(out) :: specname  ! species name
    character(len=*), optional, intent(out) :: specmorph ! species morphology
    complex(r8), pointer, optional, intent(out) :: refindex_sw(:) ! short wave species refractive indices
    complex(r8), pointer, optional, intent(out) :: refindex_lw(:) ! long wave species refractive indices
    real(r8), optional, intent(out) :: num_to_mass_aer ! ratio of number to mass concentration
    real(r8), optional, intent(out) :: dryrad  ! dry radius (m)

    character(len=20) :: aername

    if (present(density)) then
       call rad_aer_get_props(self%list_idx_, bin_ndx,  density_aer=density)
    end if

    if (present(hygro)) then
       call rad_aer_get_props(self%list_idx_, bin_ndx,  hygro_aer=hygro)
    end if
    if (present(spectype)) then

       call rad_aer_get_props(self%list_idx_, bin_ndx,  aername=aername)

       select case ( to_lower( aername(:4) ) )
       case('dust')
          spectype = 'dust'
       case('sulf','volc')
          spectype = 'sulfate'
       case('bcar','bcph')
          spectype = 'black-c'
       case('ocar','ocph')
          spectype = 'p-organic'
       case('sslt','seas','ssam','sscm')
          spectype = 'seasalt'
       case default
          spectype = 'UNKNOWN'
          call endrun('ERROR: bulk_aerosol_properties_mod%get aername not recognized : '//aername)
       end select

    end if
    if (present(specmorph)) then
      call endrun('ERROR: bulk_aerosol_properties_mod%get specmorph not yet implemented')
    end if
    if (present(specname)) then
       call rad_aer_get_props(self%list_idx_, bin_ndx,  aername=specname)
    end if
    if (present(refindex_sw)) then
       call rad_aer_get_props(self%list_idx_, bin_ndx,  refindex_aer_sw=refindex_sw)
    end if
    if (present(refindex_lw)) then
       call rad_aer_get_props(self%list_idx_, bin_ndx,  refindex_aer_lw=refindex_lw)
    end if
    if (present(num_to_mass_aer)) then
       call rad_aer_get_props(self%list_idx_, bin_ndx,  num_to_mass_aer=num_to_mass_aer)
    end if
    if (present(dryrad)) then
       call rad_aer_get_props(self%list_idx_, bin_ndx,  dryrad_aer=dryrad)
    end if

  end subroutine get

  !------------------------------------------------------------------------
  ! returns the physprop ID for a given bin (aerosol) index
  !------------------------------------------------------------------------
  integer function physprop_id(self, bin_ndx)
    use radiative_aerosol, only: rad_aer_bulk_physprop_id

    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx

    physprop_id = rad_aer_bulk_physprop_id(self%list_idx_, bin_ndx)

  end function physprop_id

  !------------------------------------------------------------------------------
  ! returns radius^3 (m3) of a given bin number
  !------------------------------------------------------------------------------
  pure elemental real(r8) function amcube(self, bin_ndx, volconc, numconc)

    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    real(r8), intent(in) :: volconc ! volume conc (m3/m3)
    real(r8), intent(in) :: numconc ! number conc (1/m3)

    amcube = nan  ! to be implemented later if needed

  end function amcube

  !------------------------------------------------------------------------------
  ! returns mass and number activation fractions
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )

    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin index
    real(r8),intent(in) :: smc       ! critical supersaturation for particles of bin radius
    real(r8),intent(in) :: smax      ! maximum supersaturation for multiple competing aerosols
    real(r8),intent(out) :: fn       ! activation fraction for aerosol number
    real(r8),intent(out) :: fm       ! activation fraction for aerosol mass

    ! to be implemented later if needed
    call endrun('ERROR: bulk_aerosol_properties_mod%actfracs not yet implemented')

  end subroutine actfracs

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol number dens
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens

    ! to be implemented later if needed
    call endrun('ERROR: bulk_aerosol_properties_mod%num_names not yet implemented')

  end subroutine num_names

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol MMR
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR

    ! to be implemented later if needed
    call endrun('ERROR: bulk_aerosol_properties_mod%mmr_names not yet implemented')

  end subroutine mmr_names

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol number dens

    ! to be implemented later if needed
    call endrun('ERROR: bulk_aerosol_properties_mod%amb_num_name not yet implemented')

  end subroutine amb_num_name

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol MMR

    ! to be implemented later if needed
    call endrun('ERROR: bulk_aerosol_properties_mod%amb_mmr_name not yet implemented')

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  ! returns species type
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx            ! bin number
    integer, intent(in) :: species_ndx        ! species number
    character(len=*), intent(out) :: spectype ! species type

    call self%get(bin_ndx, species_ndx, spectype=spectype)

  end subroutine species_type

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to given aerosol bin number
  !------------------------------------------------------------------------------
  function icenuc_updates_num(self, bin_ndx) result(res)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    logical :: res

    ! to be implemented later if needed
    res = .false.

  end function icenuc_updates_num

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to a given species within a bin
  !------------------------------------------------------------------------------
  function icenuc_updates_mmr(self, bin_ndx, species_ndx) result(res)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    logical :: res

    ! to be implemented later if needed
    res = .false.
  end function icenuc_updates_mmr

  !------------------------------------------------------------------------------
  ! apply max / min to number concentration
  !------------------------------------------------------------------------------
  subroutine apply_number_limits( self, naerosol, vaerosol, istart, istop, m )
    class(bulk_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(in)    :: vaerosol(:)  ! volume conc (m3/m3)
    integer,  intent(in) :: istart          ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop           ! stop column index
    integer,  intent(in) :: m               ! mode or bin index

    call endrun('ERROR: bulk_aerosol_properties_mod%apply_number_limits not yet implemented')

  end subroutine apply_number_limits

  !------------------------------------------------------------------------------
  ! returns TRUE if species `spc_ndx` in aerosol subset `bin_ndx` contributes to
  ! the particles' ability to act as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_species(self, bin_ndx, spc_ndx) result(res)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    integer, intent(in) :: spc_ndx  ! species number

    logical :: res

    ! to be implemented later if needed
    res = .false.
  end function hetfrz_species

  !------------------------------------------------------------------------------
  ! returns TRUE if soluble
  !------------------------------------------------------------------------------
  logical function soluble(self,bin_ndx)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    character(len=20) :: aername
    logical :: primary_carbon ! primary carbons (CB1 and OC1) are hydrophobic

    call rad_aer_get_props(self%list_idx_, bin_ndx, aername=aername)

    aername = to_lower(aername)

    primary_carbon = (aername=='bcpho') .or. (aername=='ocpho')
    soluble = .not. primary_carbon

  end function soluble

  !------------------------------------------------------------------------------
  ! returns minimum mass mean radius (meters)
  !------------------------------------------------------------------------------
  function min_mass_mean_rad(self,bin_ndx,species_ndx) result(minrad)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    real(r8) :: minrad  ! meters

    minrad = 0._r8

  end function min_mass_mean_rad

  !------------------------------------------------------------------------------
  ! returns name for a given aerosol bin
  !------------------------------------------------------------------------------
  function bin_name(self, bin_ndx) result(name)
    use aerosol_properties_mod, only: aero_name_len

    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number

    character(len=aero_name_len) :: name
    character(len=64), allocatable :: names(:)
    integer :: naer, astat

    call rad_aer_get_info(self%list_idx_, naero=naer)

    allocate( names(naer), stat=astat)
    if( astat/= 0 ) call endrun('bulk_aerosol_properties_mod%bin_name: names allocate error')

    call rad_aer_get_info(self%list_idx_, aernames=names)

    name = names(bin_ndx)

    deallocate(names)

  end function bin_name

  !------------------------------------------------------------------------------
  ! returns scavenging diameter (cm) for a given aerosol bin number
  !------------------------------------------------------------------------------
  function scav_diam(self, bin_ndx) result(diam)

    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: diam

    diam = nan ! to be implemented later if needed

  end function scav_diam

  !------------------------------------------------------------------------------
  ! adjust aerosol concentration tendencies to create larger sizes of aerosols
  ! during resuspension
  !------------------------------------------------------------------------------
  subroutine resuspension_resize(self, dcondt)

    class(bulk_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: dcondt(:)

    dcondt = nan ! to be implemented later if needed

  end subroutine resuspension_resize

  !------------------------------------------------------------------------------
  ! returns bulk deposition fluxes of the specified species type
  ! rebinned to specified diameter limits
  !------------------------------------------------------------------------------
  subroutine rebin_bulk_fluxes(self, bulk_type, dep_fluxes, diam_edges, bulk_fluxes, &
                               error_code, error_string)

    class(bulk_aerosol_properties), intent(in) :: self
    character(len=*),intent(in) :: bulk_type       ! aerosol type to rebin
    real(r8), intent(in) :: dep_fluxes(:)          ! kg/m2
    real(r8), intent(in) :: diam_edges(:)          ! meters
    real(r8), intent(out) :: bulk_fluxes(:)        ! kg/m2
    integer,  intent(out) :: error_code            ! error code (0 if no error)
    character(len=*), intent(out) :: error_string  ! error string

    ! to be implemented later if needed
    call endrun('ERROR: bulk_aerosol_properties_mod%rebin_bulk_fluxes not yet implemented')

  end subroutine rebin_bulk_fluxes

  !------------------------------------------------------------------------------
  ! Returns TRUE if bin is hydrophilic, otherwise FALSE
  !------------------------------------------------------------------------------
  logical function hydrophilic(self, bin_ndx)
    class(bulk_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx ! bin number

    hydrophilic = self%soluble(bin_ndx)

  end function hydrophilic

  !------------------------------------------------------------------------------
  ! returns TRUE if bulk aerosol representation
  !------------------------------------------------------------------------------
  pure logical function model_is(self, query)
    class(bulk_aerosol_properties), intent(in) :: self
    character(len=*),               intent(in) :: query

    if (trim(query) == 'BAM' .or. trim(query) == 'bam') then
       model_is = .true.
    else if (trim(query) == 'bulk_model') then
       model_is = .true.
    else
       model_is = .false.
    end if

  end function model_is

end module bulk_aerosol_properties_mod
