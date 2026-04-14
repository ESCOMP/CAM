module modal_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use radiative_aerosol, only: rad_aer_get_info, rad_aer_get_mode_props, rad_aer_get_props

  implicit none

  private

  public :: modal_aerosol_properties

  type, extends(aerosol_properties) :: modal_aerosol_properties
     private
     real(r8), allocatable :: exp45logsig_(:)
     real(r8), allocatable :: voltonumblo_(:)
     real(r8), allocatable :: voltonumbhi_(:)
     integer,  allocatable :: sulfate_mode_ndxs_(:)
     integer,  allocatable :: dust_mode_ndxs_(:)
     integer,  allocatable :: ssalt_mode_ndxs_(:)
     integer,  allocatable :: ammon_mode_ndxs_(:)
     integer,  allocatable :: nitrate_mode_ndxs_(:)
     integer,  allocatable :: msa_mode_ndxs_(:)
     integer,  allocatable :: bcarbon_mode_ndxs_(:,:)
     integer,  allocatable :: porganic_mode_ndxs_(:,:)
     integer,  allocatable :: sorganic_mode_ndxs_(:,:)
     integer,  allocatable :: mode_size_order_(:)
     integer :: num_soa_ = 0
     integer :: num_poa_ = 0
     integer :: num_bc_ = 0
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
  end type modal_aerosol_properties

  interface modal_aerosol_properties
     procedure :: constructor
  end interface modal_aerosol_properties

  logical, parameter :: debug = .false.

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(list_idx) result(newobj)

    integer, optional, intent(in) :: list_idx ! radiation list index (0=climate)
    type(modal_aerosol_properties), pointer :: newobj

    integer :: l, m, nmodes, ncnst_tot, mm, itmp
    integer :: list_idx_loc
    real(r8) :: dgnumlo_val
    real(r8) :: dgnumhi_val
    real(r8) :: dgnum_val
    real(r8) :: rhcrystal_val, rhdeliques_val
    integer,allocatable :: nspecies(:)
    real(r8),allocatable :: sigmag(:)
    real(r8),allocatable :: alogsig(:)
    real(r8),allocatable :: f1(:)
    real(r8),allocatable :: f2(:)
    real(r8),allocatable :: dgnum_arr(:)
    real(r8),allocatable :: dgnumhi_arr(:)
    real(r8),allocatable :: dgnumlo_arr(:)
    real(r8),allocatable :: rhcrystal_arr(:)
    real(r8),allocatable :: rhdeliques_arr(:)
    integer :: ierr

    character(len=aero_name_len) :: spectype

    integer :: npoa, nsoa, nbc

    list_idx_loc = 0
    if (present(list_idx)) list_idx_loc = list_idx

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    call rad_aer_get_info(list_idx_loc, nmodes=nmodes)

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
    allocate(dgnum_arr(nmodes),dgnumhi_arr(nmodes),dgnumlo_arr(nmodes), &
             rhcrystal_arr(nmodes),rhdeliques_arr(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%mode_size_order_(nmodes),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    ncnst_tot = 0

    do m = 1, nmodes
       call rad_aer_get_info(list_idx_loc, m, nspec=nspecies(m))

       ncnst_tot =  ncnst_tot + nspecies(m) + 1

       call rad_aer_get_mode_props(list_idx_loc, m, sigmag=sigmag(m), &
                                    dgnum=dgnum_val, dgnumhi=dgnumhi_val, dgnumlo=dgnumlo_val, &
                                    rhcrystal=rhcrystal_val, rhdeliques=rhdeliques_val)

       dgnum_arr(m) = dgnum_val
       dgnumhi_arr(m) = dgnumhi_val
       dgnumlo_arr(m) = dgnumlo_val
       rhcrystal_arr(m) = rhcrystal_val
       rhdeliques_arr(m) = rhdeliques_val

       alogsig(m) = log(sigmag(m))

       newobj%exp45logsig_(m) = exp(4.5_r8*alogsig(m)*alogsig(m))

       f1(m) = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
       f2(m) = 1._r8 + 0.25_r8*alogsig(m)

       newobj%voltonumblo_(m) = 1._r8 / ( (pi/6._r8)* &
            (dgnumlo_val**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
       newobj%voltonumbhi_(m) = 1._r8 / ( (pi/6._r8)* &
            (dgnumhi_val**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )

    end do

    ! compute mode_size_order_: indices sorted by dgnum_ descending (largest first)
    do m = 1, nmodes
       newobj%mode_size_order_(m) = m
    end do
    do m = 1, nmodes-1
       do l = m+1, nmodes
          if (dgnum_arr(newobj%mode_size_order_(l)) > dgnum_arr(newobj%mode_size_order_(m))) then
             itmp = newobj%mode_size_order_(m)
             newobj%mode_size_order_(m) = newobj%mode_size_order_(l)
             newobj%mode_size_order_(l) = itmp
          end if
       end do
    end do

    call newobj%initialize(nmodes,ncnst_tot,nspecies,nspecies,alogsig,f1,f2,ierr,list_idx_loc, &
                           dgnum=dgnum_arr,dgnumhi=dgnumhi_arr,dgnumlo=dgnumlo_arr, &
                           rhcrystal=rhcrystal_arr,rhdeliques=rhdeliques_arr)

    npoa = 0
    nsoa = 0
    nbc = 0

    m = 1
    do l = 1,newobj%nspecies(m)
       mm = newobj%indexer(m,l)
       call newobj%species_type(m, l, spectype)
       select case ( trim(spectype) )
       case('p-organic')
          npoa = npoa + 1
       case('s-organic')
          nsoa = nsoa + 1
       case('black-c')
          nbc = nbc + 1
       end select
    end do

    newobj%num_soa_ = nsoa
    newobj%num_poa_ = npoa
    newobj%num_bc_ = nbc

    allocate(newobj%sulfate_mode_ndxs_(newobj%nbins()),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%dust_mode_ndxs_(newobj%nbins()),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%ssalt_mode_ndxs_(newobj%nbins()),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%ammon_mode_ndxs_(newobj%nbins()),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%nitrate_mode_ndxs_(newobj%nbins()),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%msa_mode_ndxs_(newobj%nbins()),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    newobj%sulfate_mode_ndxs_ = 0
    newobj%dust_mode_ndxs_ = 0
    newobj%ssalt_mode_ndxs_ = 0
    newobj%ammon_mode_ndxs_ = 0
    newobj%nitrate_mode_ndxs_ = 0
    newobj%msa_mode_ndxs_ = 0

    allocate(newobj%porganic_mode_ndxs_(newobj%nbins(),npoa),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%sorganic_mode_ndxs_(newobj%nbins(),nsoa),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate(newobj%bcarbon_mode_ndxs_(newobj%nbins(),nbc),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    newobj%porganic_mode_ndxs_ = 0._r8
    newobj%sorganic_mode_ndxs_ = 0._r8
    newobj%bcarbon_mode_ndxs_ = 0._r8

    do m = 1,newobj%nbins()
       npoa = 0
       nsoa = 0
       nbc = 0

       do l = 1,newobj%nspecies(m)
          mm = newobj%indexer(m,l)
          call newobj%species_type(m, l, spectype)

          select case ( trim(spectype) )
          case('sulfate')
             newobj%sulfate_mode_ndxs_(m) = mm
          case('dust')
             newobj%dust_mode_ndxs_(m) = mm
          case('nitrate')
             newobj%nitrate_mode_ndxs_(m) = mm
          case('ammonium')
             newobj%ammon_mode_ndxs_(m) = mm
          case('seasalt')
             newobj%ssalt_mode_ndxs_(m) = mm
          case('msa')
             newobj%msa_mode_ndxs_(m) = mm
          case('p-organic')
             npoa = npoa + 1
             newobj%porganic_mode_ndxs_(m,npoa)  = mm
          case('s-organic')
             nsoa = nsoa + 1
             newobj%sorganic_mode_ndxs_(m,nsoa)  = mm
          case('black-c')
             nbc = nbc + 1
             newobj%bcarbon_mode_ndxs_(m,nbc)  = mm
          end select

       end do
    end do

    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    deallocate(nspecies)
    deallocate(alogsig)
    deallocate(sigmag)
    deallocate(f1)
    deallocate(f2)
    deallocate(dgnum_arr)
    deallocate(dgnumhi_arr)
    deallocate(dgnumlo_arr)
    deallocate(rhcrystal_arr)
    deallocate(rhdeliques_arr)

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
    if (allocated(self%mode_size_order_)) then
       deallocate(self%mode_size_order_)
    end if

    if (allocated(self%sulfate_mode_ndxs_)) then
       deallocate(self%sulfate_mode_ndxs_)
    end if
    if (allocated(self%dust_mode_ndxs_)) then
       deallocate(self%dust_mode_ndxs_)
    end if
    if (allocated(self%ssalt_mode_ndxs_)) then
       deallocate(self%ssalt_mode_ndxs_)
    end if
    if (allocated(self%ammon_mode_ndxs_)) then
       deallocate(self%ammon_mode_ndxs_)
    end if
    if (allocated(self%nitrate_mode_ndxs_)) then
       deallocate(self%nitrate_mode_ndxs_)
    end if
    if (allocated(self%msa_mode_ndxs_)) then
       deallocate(self%msa_mode_ndxs_)
    end if
    if (allocated(self%porganic_mode_ndxs_)) then
       deallocate(self%porganic_mode_ndxs_)
    end if
    if (allocated(self%sorganic_mode_ndxs_)) then
       deallocate(self%sorganic_mode_ndxs_)
    end if
    if (allocated(self%bcarbon_mode_ndxs_)) then
       deallocate(self%bcarbon_mode_ndxs_)
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
  !  species name
  !  short wave species refractive indices
  !  long wave species refractive indices
  !  species morphology
  !------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, density, hygro, &
                 spectype, specname, specmorph, refindex_sw, refindex_lw, num_to_mass_aer, &
                 dryrad)
    use cam_abortutils, only: endrun

    class(modal_aerosol_properties), intent(in) :: self
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

    call rad_aer_get_props(self%list_idx_, bin_ndx, species_ndx, &
                                density_aer=density, hygro_aer=hygro, spectype=spectype, &
                                refindex_aer_sw=refindex_sw, refindex_aer_lw=refindex_lw)

    if (present(specname)) then
       call rad_aer_get_info(self%list_idx_, bin_ndx, species_ndx, spec_name=specname)
    end if

    if (present(specmorph)) then
       specmorph = 'UNKNOWN'
    end if

    if (present(num_to_mass_aer)) then
       ! num_to_mass_aer for modal aerosols should not be read from file
       call endrun('modal_aerosol_properties_mod%get: num_to_mass_aer should not be read from file for modal aerosols')
    end if

    if (present(dryrad)) then
       ! dryrad for modal aerosols should not be read from file
       call endrun('modal_aerosol_properties_mod%get: dryrad should not be read from file for modal aerosols')
    end if

  end subroutine get

  !------------------------------------------------------------------------
  ! returns the physprop ID for a given bin (mode) index
  !------------------------------------------------------------------------
  integer function physprop_id(self, bin_ndx)
    use radiative_aerosol, only: rad_aer_mode_physprop_id

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx

    physprop_id = rad_aer_mode_physprop_id(self%list_idx_, bin_ndx)

  end function physprop_id

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

    call rad_aer_get_info(self%list_idx_,bin_ndx, num_name=name_a, num_name_cw=name_c)
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

    call rad_aer_get_info(self%list_idx_, bin_ndx, species_ndx, spec_name=name_a, spec_name_cw=name_c)
  end subroutine mmr_names

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol number dens

    call rad_aer_get_info(self%list_idx_,bin_ndx, num_name=name)

  end subroutine amb_num_name

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol MMR

    call rad_aer_get_info(self%list_idx_, bin_ndx, species_ndx, spec_name=name)

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  ! returns species type
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: spectype ! species type

    call rad_aer_get_info(self%list_idx_, bin_ndx, species_ndx, spec_type=spectype)

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

    call rad_aer_get_info(self%list_idx_, bin_ndx, mode_type=modetype)
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

       call rad_aer_get_info(self%list_idx_, bin_ndx, mode_type=modetype)
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

    call rad_aer_get_info(self%list_idx_, bin_ndx, mode_type=mode_name)

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

    call rad_aer_get_info(self%list_idx_, bin_ndx, mode_type=mode_name)

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
       call rad_aer_get_info(self%list_idx_, bin_ndx, mode_type=mode_type)
       select case ( trim(mode_type) )
       case ('accum','fine_dust')
          minrad = 0.258e-6_r8
       case ('coarse','coarse_dust')
          minrad = 1.576e-6_r8
       case default
          minrad = -huge(1._r8)
       end select
    case('black-c')
       call rad_aer_get_info(self%list_idx_, nmodes=nmodes)
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
  ! returns name for a given aerosol bin
  !------------------------------------------------------------------------------
  function bin_name(self, bin_ndx) result(name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number

    character(len=32) :: name

    call rad_aer_get_info(self%list_idx_, bin_ndx, mode_type=name)

  end function bin_name

  !------------------------------------------------------------------------------
  ! returns scavenging diameter (cm) for a given aerosol bin number
  !------------------------------------------------------------------------------
  function scav_diam(self, bin_ndx) result(diam)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: diam

    diam = self%dgnum(bin_ndx)

  end function scav_diam

  !------------------------------------------------------------------------------
  ! adjust aerosol concentration tendencies to create larger sizes of aerosols
  ! during resuspension
  !------------------------------------------------------------------------------
  subroutine resuspension_resize(self, dcondt)

    class(modal_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: dcondt(:)

    integer :: i
    character(len=4) :: spcstr

    call accumulate_to_larger_mode( 'SO4', self%sulfate_mode_ndxs_, dcondt )
    call accumulate_to_larger_mode( 'DUST',self%dust_mode_ndxs_,dcondt )
    call accumulate_to_larger_mode( 'NACL',self%ssalt_mode_ndxs_,dcondt )
    call accumulate_to_larger_mode( 'MSA', self%msa_mode_ndxs_, dcondt )
    call accumulate_to_larger_mode( 'NH4', self%ammon_mode_ndxs_, dcondt )
    call accumulate_to_larger_mode( 'NO3', self%nitrate_mode_ndxs_, dcondt )

    spcstr = '    '
    do i = 1,self%num_soa_
       write(spcstr,'(i4)') i
       call accumulate_to_larger_mode( 'SOA'//adjustl(spcstr), self%sorganic_mode_ndxs_(:,i), dcondt )
    enddo
    spcstr = '    '
    do i = 1,self%num_poa_
       write(spcstr,'(i4)') i
       call accumulate_to_larger_mode( 'POM'//adjustl(spcstr), self%porganic_mode_ndxs_(:,i), dcondt )
    enddo
    spcstr = '    '
    do i = 1,self%num_bc_
       write(spcstr,'(i4)') i
       call accumulate_to_larger_mode( 'BC'//adjustl(spcstr), self%bcarbon_mode_ndxs_(:,i), dcondt )
    enddo

  contains

    !------------------------------------------------------------------------------
    subroutine accumulate_to_larger_mode( spc_name, lptr, prevap )

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      character(len=*), intent(in) :: spc_name
      integer,  intent(in) :: lptr(:)
      real(r8), intent(inout) :: prevap(:)

      integer :: m,n, nl,ns

      logical, parameter :: debug = .false.

      ! find constituent index of the largest mode for the species
      loop1: do m = 1,self%nbins()-1
         nl = lptr(self%mode_size_order_(m))
         if (nl>0) exit loop1
      end do loop1

      if (.not. nl>0) return

      ! accumulate the smaller modes into the largest mode
      do n = m+1,self%nbins()
         ns = lptr(self%mode_size_order_(n))
         if (ns>0) then
            prevap(nl) = prevap(nl) + prevap(ns)
            prevap(ns) = 0._r8
            if (masterproc .and. debug) then
               write(iulog,'(a,i3,a,i3)') trim(spc_name)//' mode number accumulate ',ns,'->',nl
            endif
         endif
      end do

    end subroutine accumulate_to_larger_mode
    !------------------------------------------------------------------------------

  end subroutine resuspension_resize

  !------------------------------------------------------------------------------
  ! returns bulk deposition fluxes of the specified species type
  ! rebinned to specified diameter limits
  !------------------------------------------------------------------------------
  subroutine rebin_bulk_fluxes(self, bulk_type, dep_fluxes, diam_edges, bulk_fluxes, &
                               error_code, error_string)
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)

    class(modal_aerosol_properties), intent(in) :: self
    character(len=*),intent(in) :: bulk_type       ! aerosol type to rebin
    real(r8), intent(in) :: dep_fluxes(:)          ! kg/m2
    real(r8), intent(in) :: diam_edges(:)          ! meters
    real(r8), intent(out) :: bulk_fluxes(:)        ! kg/m2
    integer,  intent(out) :: error_code            ! error code (0 if no error)
    character(len=*), intent(out) :: error_string  ! error string

    real(r8) :: dns_dst ! kg/m3
    real(r8) :: sigma_g, vmd, tmp, massfrac_bin(size(bulk_fluxes))
    real(r8) :: Ntype, Mtype, Mtotal, Ntot
    integer :: k,l,m,mm, nbulk
    logical :: has_type, type_not_found

    character(len=aero_name_len) :: spectype
    character(len=aero_name_len) :: modetype

    real(r8), parameter :: sqrtwo = sqrt(2._r8)
    real(r8), parameter :: onethrd = 1._r8/3._r8

    error_code = 0
    error_string = ' '

    type_not_found = .true.

    nbulk = size(bulk_fluxes)

    bulk_fluxes(:) = 0.0_r8

    do m = 1,self%nbins()
       Mtype = 0._r8
       Mtotal = 0._r8
       mm = self%indexer(m,0)
       Ntot = dep_fluxes(mm) ! #/m2

       has_type = .false.

       do l = 1,self%nspecies(m)
          mm = self%indexer(m,l)
          call self%get(m,l, spectype=spectype, density=dns_dst) ! kg/m3
          if (spectype==bulk_type) then
             Mtype = dep_fluxes(mm) ! kg/m2
             has_type = .true.
             type_not_found = .false.
          end if
          Mtotal = Mtotal + dep_fluxes(mm) ! kg/m2
       end do
       mode_has_type: if (has_type) then
          call rad_aer_get_info(self%list_idx_, m, mode_type=modetype)
          if (Ntot>1.e-40_r8 .and. Mtype>1.e-40_r8 .and. Mtotal>1.e-40_r8) then

             call rad_aer_get_mode_props(self%list_idx_, m, sigmag=sigma_g)
             tmp = sqrtwo*log(sigma_g)

             ! type number concentration
             Ntype = Ntot * Mtype/Mtotal ! #/m2

             ! volume median diameter (meters)
             vmd = (6._r8*Mtype/(pi*Ntype*dns_dst))**onethrd * exp(1.5_r8*(log(sigma_g))**2)

             massfrac_bin = 0._r8

             do k = 1,nbulk
                massfrac_bin(k) = 0.5_r8*( erf((log(diam_edges(k+1)/vmd))/tmp) &
                                - erf((log(diam_edges(k  )/vmd))/tmp) )
                bulk_fluxes(k) = bulk_fluxes(k) + massfrac_bin(k) * Mtype
             end do

             if (debug) then
                if (abs(1._r8-sum(massfrac_bin)) > 1.e-6_r8) then
                   write(*,*) 'rebin_bulk_fluxes WARNING mode-num, massfrac_bin, sum(massfrac_bin) = ', &
                        m, massfrac_bin, sum(massfrac_bin)
                end if
             end if

          end if
       end if mode_has_type
    end do

    if (type_not_found) then
       bulk_fluxes(:) = nan
       error_code = 1
       write(error_string,*) 'aerosol_properties::rebin_bulk_fluxes ERROR : ',trim(bulk_type),' not found'
    end if

  end subroutine rebin_bulk_fluxes

  !------------------------------------------------------------------------------
  ! Returns TRUE if bin is hydrophilic, otherwise FALSE
  !------------------------------------------------------------------------------
  logical function hydrophilic(self, bin_ndx)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx ! bin number

    character(len=aero_name_len) :: modetype

    call rad_aer_get_info(self%list_idx_, bin_ndx, mode_type=modetype)

    hydrophilic = (trim(modetype) == 'accum')

  end function hydrophilic

  !------------------------------------------------------------------------------
  ! returns TRUE if modal aerosol representation
  !------------------------------------------------------------------------------
  pure logical function model_is(self, query)
    class(modal_aerosol_properties), intent(in) :: self
    character(len=*),               intent(in) :: query

    if (trim(query) == 'MAM' .or. trim(query) == 'mam') then
       model_is = .true.
    else if (trim(query) == 'modal') then
       model_is = .true.
    else
       model_is = .false.
    end if

  end function model_is

end module modal_aerosol_properties_mod
