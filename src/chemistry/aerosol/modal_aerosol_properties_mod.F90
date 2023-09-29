module modal_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props, rad_cnst_get_aer_props

  implicit none

  private

  public :: modal_aerosol_properties

  type, extends(aerosol_properties) :: modal_aerosol_properties
     private
     real(r8), allocatable :: exp45logsig_(:)
     real(r8), allocatable :: voltonumblo_(:)
     real(r8), allocatable :: voltonumbhi_(:)
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
     procedure :: optics_params
     procedure :: nbins_rlist
     procedure :: nspecies_per_bin_rlist
     procedure :: alogsig_rlist
     procedure :: soluble
     procedure :: min_mass_mean_rad
     procedure :: bin_name

     final :: destructor
  end type modal_aerosol_properties

  interface modal_aerosol_properties
     procedure :: constructor
  end interface modal_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)

    type(modal_aerosol_properties), pointer :: newobj

    integer :: m, nmodes, ncnst_tot
    real(r8) :: dgnumlo
    real(r8) :: dgnumhi
    integer,allocatable :: nspecies(:)
    real(r8),allocatable :: sigmag(:)
    real(r8),allocatable :: alogsig(:)
    real(r8),allocatable :: f1(:)
    real(r8),allocatable :: f2(:)
    integer :: ierr

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    call rad_cnst_get_info(0, nmodes=nmodes)

    allocate(nspecies(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(alogsig(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( f1(nmodes),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( f2(nmodes),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    allocate(sigmag(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%exp45logsig_(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%voltonumblo_(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%voltonumbhi_(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    ncnst_tot = 0

    do m = 1, nmodes
       call rad_cnst_get_info(0, m, nspec=nspecies(m))

       ncnst_tot =  ncnst_tot + nspecies(m) + 1

       call rad_cnst_get_mode_props(0, m, sigmag=sigmag(m), &
                                    dgnumhi=dgnumhi, dgnumlo=dgnumlo )

       alogsig(m) = log(sigmag(m))

       newobj%exp45logsig_(m) = exp(4.5_r8*alogsig(m)*alogsig(m))

       f1(m) = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
       f2(m) = 1._r8 + 0.25_r8*alogsig(m)

       newobj%voltonumblo_(m) = 1._r8 / ( (pi/6._r8)* &
            (dgnumlo**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
       newobj%voltonumbhi_(m) = 1._r8 / ( (pi/6._r8)* &
            (dgnumhi**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )

    end do

    call newobj%initialize(nmodes,ncnst_tot,nspecies,nspecies,alogsig,f1,f2,ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    deallocate(nspecies)
    deallocate(alogsig)
    deallocate(sigmag)
    deallocate(f1)
    deallocate(f2)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_properties), intent(inout) :: self

    if (allocated(self%exp45logsig_)) then
       deallocate(self%exp45logsig_)
    end if
    if (allocated(self%voltonumblo_)) then
       deallocate(self%voltonumblo_)
    end if
    if (allocated(self%voltonumbhi_)) then
       deallocate(self%voltonumbhi_)
    end if

    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! returns number of transported aerosol constituents
  !------------------------------------------------------------------------------
  integer function number_transported(self)
    class(modal_aerosol_properties), intent(in) :: self
    ! to be implemented later
    number_transported = -1
  end function number_transported

  !------------------------------------------------------------------------
  ! returns aerosol properties:
  !  density
  !  hygroscopicity
  !  species type
  !  short wave species refractive indices
  !  long wave species refractive indices
  !  species morphology
  !------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, list_ndx, density, hygro, &
                 spectype, specmorph, refindex_sw, refindex_lw)

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    integer, optional, intent(in) :: list_ndx  ! climate or a diagnostic list number
    real(r8), optional, intent(out) :: density ! density (kg/m3)
    real(r8), optional, intent(out) :: hygro   ! hygroscopicity
    character(len=*), optional, intent(out) :: spectype  ! species type
    character(len=*), optional, intent(out) :: specmorph ! species morphology
    complex(r8), pointer, optional, intent(out) :: refindex_sw(:) ! short wave species refractive indices
    complex(r8), pointer, optional, intent(out) :: refindex_lw(:) ! long wave species refractive indices

    integer :: ilist

    if (present(list_ndx)) then
       ilist = list_ndx
    else
       ilist = 0
    end if

    call rad_cnst_get_aer_props(ilist, bin_ndx, species_ndx, &
                                density_aer=density, hygro_aer=hygro, spectype=spectype, &
                                refindex_aer_sw=refindex_sw, refindex_aer_lw=refindex_lw)

    if (present(specmorph)) then
       specmorph = 'UNKNOWN'
    end if

  end subroutine get

  !------------------------------------------------------------------------
  ! returns optics type and table parameters
  !------------------------------------------------------------------------
  subroutine optics_params(self, list_ndx, bin_ndx, opticstype, extpsw, abspsw, asmpsw, absplw, &
       refrtabsw, refitabsw, refrtablw, refitablw, ncoef, prefr, prefi, sw_hygro_ext_wtp, &
       sw_hygro_ssa_wtp, sw_hygro_asm_wtp, lw_hygro_ext_wtp, wgtpct, nwtp, &
       sw_hygro_coreshell_ext, sw_hygro_coreshell_ssa, sw_hygro_coreshell_asm, lw_hygro_coreshell_ext, &
       corefrac, bcdust, kap, relh, nfrac, nbcdust, nkap, nrelh )

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: list_ndx            ! rad climate/diags list

    character(len=*), optional, intent(out) :: opticstype

    ! refactive index table parameters
    real(r8),  optional, pointer     :: extpsw(:,:,:,:) ! short wave specific extinction
    real(r8),  optional, pointer     :: abspsw(:,:,:,:) ! short wave specific absorption
    real(r8),  optional, pointer     :: asmpsw(:,:,:,:) ! short wave asymmetry factor
    real(r8),  optional, pointer     :: absplw(:,:,:,:) ! long wave specific absorption
    real(r8),  optional, pointer     :: refrtabsw(:,:)  ! table of short wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitabsw(:,:)  ! table of short wave imaginary refractive indices for aerosols
    real(r8),  optional, pointer     :: refrtablw(:,:)  ! table of long wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitablw(:,:)  ! table of long wave imaginary refractive indices for aerosols
    integer,   optional, intent(out) :: ncoef  ! number of chebychev polynomials
    integer,   optional, intent(out) :: prefr  ! number of real refractive indices in table
    integer,   optional, intent(out) :: prefi  ! number of imaginary refractive indices in table

    ! hygrowghtpct table parameters
    real(r8),  optional, pointer     :: sw_hygro_ext_wtp(:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_ssa_wtp(:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_asm_wtp(:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_ext_wtp(:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: wgtpct(:)   ! weight precent of H2SO4/H2O solution
    integer,   optional, intent(out) :: nwtp        ! number of weight precent values

    ! hygrocoreshell table parameters
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ext(:,:,:,:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ssa(:,:,:,:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_asm(:,:,:,:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_coreshell_ext(:,:,:,:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: corefrac(:) ! core fraction dimension values
    real(r8),  optional, pointer     :: bcdust(:)   ! bc/(bc + dust) fraction dimension values
    real(r8),  optional, pointer     :: kap(:)      ! hygroscopicity dimension values
    real(r8),  optional, pointer     :: relh(:)     ! relative humidity dimension values
    integer,   optional, intent(out) :: nfrac       ! core fraction dimension size
    integer,   optional, intent(out) :: nbcdust     ! bc/(bc + dust) fraction dimension size
    integer,   optional, intent(out) :: nkap        ! hygroscopicity dimension size
    integer,   optional, intent(out) :: nrelh       ! relative humidity dimension size

    ! refactive index table parameters
    call rad_cnst_get_mode_props(list_ndx, bin_ndx, &
                                 opticstype=opticstype, &
                                 extpsw=extpsw, &
                                 abspsw=abspsw, &
                                 asmpsw=asmpsw, &
                                 absplw=absplw, &
                                 refrtabsw=refrtabsw, &
                                 refitabsw=refitabsw, &
                                 refrtablw=refrtablw, &
                                 refitablw=refitablw, &
                                 ncoef=ncoef, &
                                 prefr=prefr, &
                                 prefi=prefi)

    ! hygrowghtpct table parameters
    if (present(sw_hygro_ext_wtp)) then
       nullify(sw_hygro_ext_wtp)
    end if
    if (present(sw_hygro_ssa_wtp)) then
       nullify(sw_hygro_ssa_wtp)
    end if
    if (present(sw_hygro_asm_wtp)) then
       nullify(sw_hygro_asm_wtp)
    end if
    if (present(lw_hygro_ext_wtp)) then
       nullify(lw_hygro_ext_wtp)
    end if
    if (present(wgtpct)) then
       nullify(wgtpct)
    end if
    if (present(nwtp)) then
       nwtp = -1
    end if

    ! hygrocoreshell table parameters
    if (present(sw_hygro_coreshell_ext)) then
       nullify(sw_hygro_coreshell_ext)
    end if
    if (present(sw_hygro_coreshell_ssa)) then
       nullify(sw_hygro_coreshell_ssa)
    end if
    if (present(sw_hygro_coreshell_asm)) then
       nullify(sw_hygro_coreshell_asm)
    end if
    if (present(lw_hygro_coreshell_ext)) then
       nullify(lw_hygro_coreshell_ext)
    end if
    if (present(corefrac)) then
       nullify(corefrac)
    end if
    if (present(bcdust)) then
       nullify(bcdust)
    end if
    if (present(kap)) then
       nullify(kap)
    end if
    if (present(relh)) then
       nullify(relh)
    end if
    if (present(nfrac)) then
       nfrac = -1
    end if
    if (present(nbcdust)) then
       nbcdust = -1
    end if
    if (present(nkap)) then
       nkap = -1
    end if
    if (present(nrelh)) then
       nrelh = -1
    end if

  end subroutine optics_params

  !------------------------------------------------------------------------------
  ! returns radius^3 (m3) of a given bin number
  !------------------------------------------------------------------------------
  pure elemental real(r8) function amcube(self, bin_ndx, volconc, numconc)

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    real(r8), intent(in) :: volconc ! volume conc (m3/m3)
    real(r8), intent(in) :: numconc ! number conc (1/m3)

    amcube = (3._r8*volconc/(4._r8*pi*self%exp45logsig_(bin_ndx)*numconc))

  end function amcube

  !------------------------------------------------------------------------------
  ! returns mass and number activation fractions
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )
    use shr_spfn_mod, only: erf => shr_spfn_erf
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin index
    real(r8),intent(in) :: smc       ! critical supersaturation for particles of bin radius
    real(r8),intent(in) :: smax      ! maximum supersaturation for multiple competing aerosols
    real(r8),intent(out) :: fn       ! activation fraction for aerosol number
    real(r8),intent(out) :: fm       ! activation fraction for aerosol mass

    real(r8) :: x,y
    real(r8), parameter :: twothird = 2._r8/3._r8
    real(r8), parameter :: sq2      = sqrt(2._r8)

    x=twothird*(log(smc)-log(smax))/(sq2*self%alogsig(bin_ndx))
    y=x-1.5_r8*sq2*self%alogsig(bin_ndx)

    fn = 0.5_r8*(1._r8-erf(x))
    fm = 0.5_r8*(1._r8-erf(y))

  end subroutine actfracs

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol number dens
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens

    call rad_cnst_get_info(0,bin_ndx, num_name=name_a, num_name_cw=name_c)
  end subroutine num_names

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol MMR
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR

    call rad_cnst_get_info(0, bin_ndx, species_ndx, spec_name=name_a, spec_name_cw=name_c)
  end subroutine mmr_names

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol number dens

    call rad_cnst_get_info(0,bin_ndx, num_name=name)

  end subroutine amb_num_name

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol MMR

    call rad_cnst_get_info(0, bin_ndx, species_ndx, spec_name=name)

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  ! returns species type
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: spectype ! species type

    call rad_cnst_get_info(0, bin_ndx, species_ndx, spec_type=spectype)

  end subroutine species_type

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to given aerosol bin number
  !------------------------------------------------------------------------------
  function icenuc_updates_num(self, bin_ndx) result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    logical :: res

    character(len=aero_name_len) :: spectype
    character(len=aero_name_len) :: modetype
    integer :: spc_ndx

    res = .false.

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)
    if (.not.(modetype=='coarse' .or. modetype=='coarse_dust')) then
       return
    end if

    do spc_ndx = 1, self%nspecies(bin_ndx)
       call self%species_type( bin_ndx, spc_ndx, spectype)
       if (spectype=='dust') res = .true.
    end do

  end function icenuc_updates_num

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to a given species within a bin
  !------------------------------------------------------------------------------
  function icenuc_updates_mmr(self, bin_ndx, species_ndx) result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    logical :: res

    character(len=32) :: spectype
    character(len=32) :: modetype

    res = .false.

    if (species_ndx>0) then

       call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)
       if (.not.(modetype=='coarse' .or. modetype=='coarse_dust')) then
          return
       end if

       call self%species_type( bin_ndx, species_ndx, spectype)
       if (spectype=='dust') res = .true.
    end if

  end function icenuc_updates_mmr

  !------------------------------------------------------------------------------
  ! apply max / min to number concentration
  !------------------------------------------------------------------------------
  subroutine apply_number_limits( self, naerosol, vaerosol, istart, istop, m )
    class(modal_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(in)    :: vaerosol(:)  ! volume conc (m3/m3)
    integer,  intent(in) :: istart          ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop           ! stop column index
    integer,  intent(in) :: m               ! mode or bin index

    integer :: i

    ! adjust number so that dgnumlo < dgnum < dgnumhi
    ! -- the diameter falls within the lower and upper limits which are
    !    represented by voltonumhi and voltonumblo values, respectively
    do i = istart, istop
       naerosol(i) = max(naerosol(i), vaerosol(i)*self%voltonumbhi_(m))
       naerosol(i) = min(naerosol(i), vaerosol(i)*self%voltonumblo_(m))
    end do

  end subroutine apply_number_limits

  !------------------------------------------------------------------------------
  ! returns TRUE if species `spc_ndx` in aerosol subset `bin_ndx` contributes to
  ! the particles' ability to act as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_species(self, bin_ndx, spc_ndx) result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    integer, intent(in) :: spc_ndx  ! species number

    logical :: res

    character(len=aero_name_len) :: mode_name, species_type

    res = .false.

    call rad_cnst_get_info(0, bin_ndx, mode_type=mode_name)

    if ((trim(mode_name)/='aitken')) then

       call self%species_type(bin_ndx, spc_ndx, species_type)

       if ((trim(species_type)=='black-c').or.(trim(species_type)=='dust')) then

          res = .true.

       end if

    end if

  end function hetfrz_species

  !------------------------------------------------------------------------------
  ! returns TRUE if soluble
  !------------------------------------------------------------------------------
  logical function soluble(self,bin_ndx)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    character(len=aero_name_len) :: mode_name

    call rad_cnst_get_info(0, bin_ndx, mode_type=mode_name)

    soluble = trim(mode_name)/='primary_carbon'

  end function soluble

  !------------------------------------------------------------------------------
  ! returns minimum mass mean radius (meters)
  !------------------------------------------------------------------------------
  function min_mass_mean_rad(self,bin_ndx,species_ndx) result(minrad)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    real(r8) :: minrad  ! meters

    integer :: nmodes
    character(len=aero_name_len) :: species_type, mode_type

    call self%species_type(bin_ndx, species_ndx, spectype=species_type)
    select case ( trim(species_type) )
    case('dust')
       call rad_cnst_get_info(0, bin_ndx, mode_type=mode_type)
       select case ( trim(mode_type) )
       case ('accum','fine_dust')
          minrad = 0.258e-6_r8
       case ('coarse','coarse_dust')
          minrad = 1.576e-6_r8
       case default
          minrad = -huge(1._r8)
       end select
    case('black-c')
       call rad_cnst_get_info(0, nmodes=nmodes)
       if (nmodes==3) then
          minrad = 0.04e-6_r8
       else
          minrad = 0.067e-6_r8 ! from emission size
       endif
    case default
       minrad = -huge(1._r8)
    end select

  end function min_mass_mean_rad

  !------------------------------------------------------------------------------
  ! returns the total number of bins for a given radiation list index
  !------------------------------------------------------------------------------
  function nbins_rlist(self, list_ndx)  result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx  ! radiation list number

    integer :: res

    call rad_cnst_get_info(list_ndx, nmodes=res)

  end function nbins_rlist

  !------------------------------------------------------------------------------
  ! returns number of species in a bin for a given radiation list index
  !------------------------------------------------------------------------------
  function nspecies_per_bin_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    integer :: res

    call rad_cnst_get_info(list_ndx, bin_ndx, nspec=res)

  end function nspecies_per_bin_rlist

  !------------------------------------------------------------------------------
  ! returns the natural log of geometric standard deviation of the number
  ! distribution for radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function alogsig_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: res

    real(r8) :: sig

    call rad_cnst_get_mode_props(list_ndx, bin_ndx, sigmag=sig)
    res = log(sig)

  end function alogsig_rlist

  !------------------------------------------------------------------------------
  ! returns name for a given radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function bin_name(self, list_ndx,  bin_ndx) result(name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    character(len=32) name

    call rad_cnst_get_info(list_ndx, bin_ndx, mode_type=name)

  end function bin_name

end module modal_aerosol_properties_mod
