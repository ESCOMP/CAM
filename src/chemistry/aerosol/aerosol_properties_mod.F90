module aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: aerosol_properties

  !> aerosol_properties defines the configuration of any aerosol package (using
  !! any aerosol representation) based on user specification. These values are
  !! set during initialization and do not vary during the simulation.
  !!
  !! Each aerosol package (e.g., MAM, CARMA, etc) must extend the abstract
  !! aerosol_properties class to define the details of their configuration. Any
  !! package must implement each of the deferred procedures of the abstract
  !! aerosol_properties class, may include additional private data members and
  !! type-bound procedures, and may override functions of the abstract class.
  !!
  !! Please see the modal_aerosol_properties module for an example of how the
  !! aerosol_properties class can be extended for a specific aerosol package.
  type, abstract :: aerosol_properties
     private
     integer :: nbins_ = 0  ! number of aerosol bins
     integer :: ncnst_tot_ = 0 ! total number of constituents
     integer, allocatable :: nmasses_(:) ! number of species masses
     integer, allocatable :: nspecies_(:) ! number of species
     integer, allocatable :: indexer_(:,:) ! unique indices of the aerosol elements
     real(r8), allocatable :: alogsig_(:) ! natural log of geometric deviation of the number distribution for aerosol bin
     real(r8), allocatable :: f1_(:) ! eq 28 Abdul-Razzak et al 1998
     real(r8), allocatable :: f2_(:) ! eq 29 Abdul-Razzak et al 1998
     ! Abdul-Razzak, H., S.J. Ghan, and C. Rivera-Carpio, A parameterization of aerosol activation,
     ! 1, Singleaerosoltype. J. Geophys. Res., 103, 6123-6132, 1998.
     real(r8) :: soa_equivso4_factor_ = -huge(1._r8)
     real(r8) :: pom_equivso4_factor_ = -huge(1._r8)
   contains
     procedure :: initialize => aero_props_init
     procedure,private :: nbins_0list
     procedure(aero_nbins_rlist), deferred :: nbins_rlist
     generic :: nbins => nbins_0list,nbins_rlist
     procedure :: ncnst_tot
     procedure,private :: nspecies_per_bin
     procedure(aero_nspecies_rlist), deferred :: nspecies_per_bin_rlist
     procedure,private :: nspecies_all_bins
     generic :: nspecies => nspecies_all_bins,nspecies_per_bin,nspecies_per_bin_rlist
     procedure,private :: n_masses_all_bins
     procedure,private :: n_masses_per_bin
     generic :: nmasses => n_masses_all_bins,n_masses_per_bin
     procedure :: indexer
     procedure :: maxsat
     procedure(aero_amcube), deferred :: amcube
     procedure :: alogsig_0list
     procedure(aero_alogsig_rlist), deferred :: alogsig_rlist
     generic :: alogsig =>  alogsig_0list,alogsig_rlist
     procedure(aero_number_transported), deferred :: number_transported
     procedure(aero_props_get), deferred :: get
     procedure(aero_actfracs), deferred :: actfracs
     procedure(aero_num_names), deferred :: num_names
     procedure(aero_mmr_names), deferred :: mmr_names
     procedure(aero_amb_num_name), deferred :: amb_num_name
     procedure(aero_amb_mmr_name), deferred :: amb_mmr_name
     procedure(aero_species_type), deferred :: species_type
     procedure(aero_icenuc_updates_num), deferred :: icenuc_updates_num
     procedure(aero_icenuc_updates_mmr), deferred :: icenuc_updates_mmr
     procedure(aero_apply_num_limits), deferred :: apply_number_limits
     procedure(aero_hetfrz_species), deferred :: hetfrz_species
     procedure :: soa_equivso4_factor ! SOA Hygroscopicity / Sulfate Hygroscopicity
     procedure :: pom_equivso4_factor ! POM Hygroscopicity / Sulfate Hygroscopicity
     procedure(aero_soluble), deferred :: soluble
     procedure(aero_min_mass_mean_rad), deferred :: min_mass_mean_rad
     procedure(aero_optics_params), deferred :: optics_params
     procedure(aero_bin_name), deferred :: bin_name

     procedure :: final=>aero_props_final
  end type aerosol_properties

  integer,public, parameter :: aero_name_len = 32 ! common length of aersols names, species, etc

  abstract interface

     !------------------------------------------------------------------------------
     ! returns number of transported aerosol constituents
     !------------------------------------------------------------------------------
     integer function aero_number_transported(self)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
     end function aero_number_transported

     !------------------------------------------------------------------------
     ! returns aerosol properties:
     !  density
     !  hygroscopicity
     !  species type
     !  short wave species refractive indices
     !  long wave species refractive indices
     !  species morphology
     !------------------------------------------------------------------------
     subroutine aero_props_get(self, bin_ndx, species_ndx, list_ndx, density, hygro, &
                               spectype, specmorph, refindex_sw, refindex_lw)
       import :: aerosol_properties, r8
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx             ! bin index
       integer, intent(in) :: species_ndx         ! species index
       integer, optional, intent(in) :: list_ndx  ! climate or a diagnostic list number
       real(r8), optional, intent(out) :: density ! density (kg/m3)
       real(r8), optional, intent(out) :: hygro   ! hygroscopicity
       character(len=*), optional, intent(out) :: spectype  ! species type
       character(len=*), optional, intent(out) :: specmorph ! species morphology
       complex(r8), pointer, optional, intent(out) :: refindex_sw(:) ! short wave species refractive indices
       complex(r8), pointer, optional, intent(out) :: refindex_lw(:) ! long wave species refractive indices

     end subroutine aero_props_get

     !------------------------------------------------------------------------
     ! returns optics type and table parameters
     !------------------------------------------------------------------------
     subroutine aero_optics_params(self, list_ndx, bin_ndx, opticstype, extpsw, abspsw, asmpsw, absplw, &
          refrtabsw, refitabsw, refrtablw, refitablw, ncoef, prefr, prefi, sw_hygro_ext_wtp, &
          sw_hygro_ssa_wtp, sw_hygro_asm_wtp, lw_hygro_ext_wtp, wgtpct, nwtp, &
          sw_hygro_coreshell_ext, sw_hygro_coreshell_ssa, sw_hygro_coreshell_asm, lw_hygro_coreshell_ext, &
          corefrac, bcdust, kap, relh, nfrac, nbcdust, nkap, nrelh )

       import :: aerosol_properties, r8

       class(aerosol_properties), intent(in) :: self
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

     end subroutine aero_optics_params

     !------------------------------------------------------------------------
     ! returns species type
     !------------------------------------------------------------------------
     subroutine aero_species_type(self, bin_ndx, species_ndx, spectype)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       integer, intent(in) :: species_ndx       ! species number
       character(len=*), intent(out) :: spectype ! species type

     end subroutine aero_species_type

     !------------------------------------------------------------------------
     ! returns mass and number activation fractions
     !------------------------------------------------------------------------
     subroutine aero_actfracs(self, bin_ndx, smc, smax, fn, fm )
       import :: aerosol_properties, r8
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx   ! bin index
       real(r8),intent(in) :: smc       ! critical supersaturation for particles of bin radius
       real(r8),intent(in) :: smax      ! maximum supersaturation for multiple competing aerosols
       real(r8),intent(out) :: fn       ! activation fraction for aerosol number
       real(r8),intent(out) :: fm       ! activation fraction for aerosol mass

     end subroutine aero_actfracs

     !------------------------------------------------------------------------
     ! returns constituents names of aerosol number mixing ratios
     !------------------------------------------------------------------------
     subroutine aero_num_names(self, bin_ndx, name_a, name_c)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol number dens
       character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens
     end subroutine aero_num_names

     !------------------------------------------------------------------------
     ! returns constituents names of aerosol mass mixing ratios
     !------------------------------------------------------------------------
     subroutine aero_mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       integer, intent(in) :: species_ndx       ! species number
       character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol MMR
       character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR
     end subroutine aero_mmr_names

     !------------------------------------------------------------------------
     ! returns constituent name of ambient aerosol number mixing ratios
     !------------------------------------------------------------------------
     subroutine aero_amb_num_name(self, bin_ndx, name)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx          ! bin number
       character(len=*), intent(out) :: name  ! constituent name of ambient aerosol number dens

     end subroutine aero_amb_num_name

     !------------------------------------------------------------------------
     ! returns constituent name of ambient aerosol mass mixing ratios
     !------------------------------------------------------------------------
     subroutine aero_amb_mmr_name(self, bin_ndx, species_ndx, name)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       integer, intent(in) :: species_ndx       ! species number
       character(len=*), intent(out) :: name   ! constituent name of ambient aerosol MMR

     end subroutine aero_amb_mmr_name

     !------------------------------------------------------------------------------
     ! returns radius^3 (m3) of a given bin number
     !------------------------------------------------------------------------------
     pure elemental real(r8) function aero_amcube(self, bin_ndx, volconc, numconc)
       import :: aerosol_properties, r8

       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx  ! bin number
       real(r8), intent(in) :: volconc ! volume conc (m3/m3)
       real(r8), intent(in) :: numconc ! number conc (1/m3)

     end function aero_amcube

     !------------------------------------------------------------------------------
     ! returns TRUE if Ice Nucleation tendencies are applied to given aerosol bin number
     !------------------------------------------------------------------------------
     function aero_icenuc_updates_num(self, bin_ndx) result(res)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx  ! bin number

       logical :: res

     end function aero_icenuc_updates_num

     !------------------------------------------------------------------------------
     ! returns TRUE if Ice Nucleation tendencies are applied to a given species within a bin
     !------------------------------------------------------------------------------
     function aero_icenuc_updates_mmr(self, bin_ndx, species_ndx) result(res)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx     ! bin number
       integer, intent(in) :: species_ndx ! species number

       logical :: res

     end function aero_icenuc_updates_mmr

     !------------------------------------------------------------------------------
     ! apply max / min to number concentration
     !------------------------------------------------------------------------------
     subroutine aero_apply_num_limits( self, naerosol, vaerosol, istart, istop, m )
       import :: aerosol_properties, r8
       class(aerosol_properties), intent(in) :: self
       real(r8), intent(inout) :: naerosol(:)  ! number conc (1/m3)
       real(r8), intent(in)    :: vaerosol(:)  ! volume conc (m3/m3)
       integer,  intent(in) :: istart          ! start column index (1 <= istart <= istop <= pcols)
       integer,  intent(in) :: istop           ! stop column index
       integer,  intent(in) :: m               ! mode or bin index

     end subroutine aero_apply_num_limits

     !------------------------------------------------------------------------------
     ! returns TRUE if species `spc_ndx` in aerosol subset `bin_ndx` contributes to
     ! the particles' ability to act as heterogeneous freezing nuclei
     !------------------------------------------------------------------------------
     function aero_hetfrz_species(self, bin_ndx, spc_ndx) result(res)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx  ! bin number
       integer, intent(in) :: spc_ndx  ! species number

       logical :: res

     end function aero_hetfrz_species

     !------------------------------------------------------------------------------
     ! returns minimum mass mean radius (meters)
     !------------------------------------------------------------------------------
     function aero_min_mass_mean_rad(self,bin_ndx,species_ndx) result(minrad)
       import :: aerosol_properties, r8
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       integer, intent(in) :: species_ndx       ! species number

       real(r8) :: minrad  ! meters

     end function aero_min_mass_mean_rad

     !------------------------------------------------------------------------------
     ! returns TRUE if soluble
     !------------------------------------------------------------------------------
     logical function aero_soluble(self,bin_ndx)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number

     end function aero_soluble

     !------------------------------------------------------------------------------
     ! returns the total number of bins for a given radiation list index
     !------------------------------------------------------------------------------
     function aero_nbins_rlist(self, list_ndx)  result(res)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: list_ndx  ! radiation list number

       integer :: res

     end function aero_nbins_rlist

     !------------------------------------------------------------------------------
     ! returns number of species in a bin for a given radiation list index
     !------------------------------------------------------------------------------
     function aero_nspecies_rlist(self, list_ndx,  bin_ndx)  result(res)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: list_ndx ! radiation list number
       integer, intent(in) :: bin_ndx  ! bin number

       integer :: res

     end function aero_nspecies_rlist

     !------------------------------------------------------------------------------
     ! returns the natural log of geometric standard deviation of the number
     ! distribution for radiation list number and aerosol bin
     !------------------------------------------------------------------------------
     function aero_alogsig_rlist(self, list_ndx,  bin_ndx)  result(res)
       import :: aerosol_properties, r8
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: list_ndx ! radiation list number
       integer, intent(in) :: bin_ndx  ! bin number

       real(r8) :: res

     end function aero_alogsig_rlist

     !------------------------------------------------------------------------------
     ! returns name for a given radiation list number and aerosol bin
     !------------------------------------------------------------------------------
     function aero_bin_name(self, list_ndx,  bin_ndx) result(name)
       import :: aerosol_properties, r8
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: list_ndx ! radiation list number
       integer, intent(in) :: bin_ndx  ! bin number

       character(len=32) name

     end function aero_bin_name

  end interface

contains

  !------------------------------------------------------------------------------
  ! object initializer
  !------------------------------------------------------------------------------
  subroutine aero_props_init(self, nbin, ncnst, nspec, nmasses, alogsig, f1,f2, ierr )
    class(aerosol_properties), intent(inout) :: self
    integer, intent(in) :: nbin               ! number of bins
    integer, intent(in) :: ncnst              ! total number of constituents
    integer, intent(in) :: nspec(nbin)        ! number of species in each bin
    integer, intent(in) :: nmasses(nbin)      ! number of masses in each bin
    real(r8),intent(in) :: alogsig(nbin)      ! natural log of the standard deviation (sigma) of the aerosol bins
    real(r8),intent(in) :: f1(nbin)           ! eq 28 Abdul-Razzak et al 1998
    real(r8),intent(in) :: f2(nbin)           ! eq 29 Abdul-Razzak et al 1998
    integer,intent(out) :: ierr

    integer :: imas,ibin,indx
    character(len=*),parameter :: prefix = 'aerosol_properties::aero_props_init: '

    real(r8), parameter :: spechygro_so4 = 0.507_r8          ! Sulfate hygroscopicity
    real(r8), parameter :: spechygro_soa = 0.14_r8           ! SOA hygroscopicity
    real(r8), parameter :: spechygro_pom = 0.1_r8            ! POM hygroscopicity

    ierr = 0

    allocate(self%nspecies_(nbin),stat=ierr)
    if( ierr /= 0 ) then
       return
    end if
    allocate(self%nmasses_(nbin),stat=ierr)
    if( ierr /= 0 ) then
       return
    end if
    allocate(self%alogsig_(nbin),stat=ierr)
    if( ierr /= 0 ) then
       return
    end if
    allocate(self%f1_(nbin),stat=ierr)
    if( ierr /= 0 ) then
       return
    end if
    allocate(self%f2_(nbin),stat=ierr)
    if( ierr /= 0 ) then
       return
    end if

    allocate( self%indexer_(nbin,0:maxval(nmasses)),stat=ierr )
    if( ierr /= 0 ) then
       return
    end if

    ! Local indexing compresses the mode and number/mass indices into one index.
    ! This indexing is used by the pointer arrays used to reference state and pbuf
    ! fields. We add number = 0, total mass = 1 (if available), and mass from each
    ! constituency into mm.

    self%indexer_ = -1
    indx = 0

    do ibin=1,nbin
       do imas = 0,nmasses(ibin)
          indx = indx+1
          self%indexer_(ibin,imas) = indx
       end do
    end do

    self%nbins_ = nbin
    self%ncnst_tot_ = ncnst
    self%nmasses_(:) = nmasses(:)
    self%nspecies_(:) = nspec(:)
    self%alogsig_(:) = alogsig(:)
    self%f1_(:) = f1(:)
    self%f2_(:) = f2(:)

    self%soa_equivso4_factor_ = spechygro_soa/spechygro_so4
    self%pom_equivso4_factor_ = spechygro_pom/spechygro_so4

  end subroutine aero_props_init

  !------------------------------------------------------------------------------
  ! Object clean
  !------------------------------------------------------------------------------
  subroutine aero_props_final(self)
    class(aerosol_properties), intent(inout) :: self

    if (allocated(self%nspecies_)) then
       deallocate(self%nspecies_)
    end if
    if (allocated(self%nmasses_)) then
       deallocate(self%nmasses_)
    end if
    if (allocated(self%indexer_)) then
       deallocate(self%indexer_)
    endif
    if (allocated(self%alogsig_)) then
       deallocate(self%alogsig_)
    endif
    if (allocated(self%f1_)) then
       deallocate(self%f1_)
    endif
    if (allocated(self%f2_)) then
       deallocate(self%f2_)
    endif

    self%nbins_ = 0
    self%ncnst_tot_ = 0

  end subroutine aero_props_final

  !------------------------------------------------------------------------------
  ! returns number of species in a bin
  !------------------------------------------------------------------------------
  pure function nspecies_per_bin(self,bin_ndx) result(val)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer :: val

    val = self%nspecies_(bin_ndx)
  end function nspecies_per_bin

  !------------------------------------------------------------------------------
  ! returns number of species for all bins
  !------------------------------------------------------------------------------
  pure function nspecies_all_bins(self) result(arr)
    class(aerosol_properties), intent(in) :: self
    integer :: arr(self%nbins_)

    arr(:) = self%nspecies_(:)

  end function nspecies_all_bins

  !------------------------------------------------------------------------------
  ! returns number of species masses in a given bin number
  !------------------------------------------------------------------------------
  pure function n_masses_per_bin(self,bin_ndx) result(val)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer :: val

    val = self%nmasses_(bin_ndx)
  end function n_masses_per_bin

  !------------------------------------------------------------------------------
  ! returns an array of number of species masses for all bins
  !------------------------------------------------------------------------------
  pure function n_masses_all_bins(self) result(arr)
    class(aerosol_properties), intent(in) :: self
    integer :: arr(self%nbins_)

    arr(:) = self%nmasses_(:)
  end function n_masses_all_bins

  !------------------------------------------------------------------------------
  ! returns a single index for given bin and species
  !------------------------------------------------------------------------------
  pure integer function indexer(self,bin_ndx,species_ndx)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    indexer = self%indexer_(bin_ndx,species_ndx)
  end function indexer

  !------------------------------------------------------------------------------
  ! returns the total number of bins
  !------------------------------------------------------------------------------
  pure function nbins_0list(self) result(nbins)
    class(aerosol_properties), intent(in) :: self
    integer :: nbins

    nbins = self%nbins_
  end function nbins_0list

  !------------------------------------------------------------------------------
  ! returns number of constituents (or elements) totaled across all bins
  !------------------------------------------------------------------------------
  pure integer function ncnst_tot(self)
    class(aerosol_properties), intent(in) :: self

    ncnst_tot = self%ncnst_tot_
  end function ncnst_tot

  !------------------------------------------------------------------------------
  ! returns the natural log of geometric standard deviation of the number distribution for aerosol bin
  !------------------------------------------------------------------------------
  pure real(r8) function alogsig_0list(self, bin_ndx)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    alogsig_0list = self%alogsig_(bin_ndx)
  end function alogsig_0list

  !------------------------------------------------------------------------------
  ! returns maximum supersaturation
  !------------------------------------------------------------------------------
  function maxsat(self, zeta,eta,smc) result(smax)

    !-------------------------------------------------------------------------
    ! Calculates maximum supersaturation for multiple competing aerosols.
    !
    ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844., 2000
    !-------------------------------------------------------------------------

    class(aerosol_properties), intent(in) :: self
    real(r8), intent(in)  :: zeta(self%nbins_) ! Abdul-Razzak and Ghan eq 10
    real(r8), intent(in)  :: eta(self%nbins_)  ! Abdul-Razzak and Ghan eq 11
    real(r8), intent(in)  :: smc(self%nbins_)  ! critical supersaturation

    real(r8) :: smax ! maximum supersaturation

    integer  :: m
    integer  :: nbins
    real(r8) :: sum, g1, g2, g1sqrt, g2sqrt

    real(r8), parameter :: small_maxsat = 1.e-20_r8 ! for weak forcing
    real(r8), parameter :: large_maxsat = 1.e20_r8  ! for small eta

    smax=0.0_r8
    nbins = self%nbins_

    check_loop: do m=1,nbins
       if((zeta(m) > 1.e5_r8*eta(m)) .or. (smc(m)*smc(m) > 1.e5_r8*eta(m))) then
          ! weak forcing -- essentially none activated
          smax=small_maxsat
       else
          ! significant activation of this mode -- calc activation all modes
          exit check_loop
       endif
       ! No significant activation in any mode.  Do nothing.
       if (m == nbins) return
    enddo check_loop

    sum=0.0_r8

    do m=1,nbins
       if(eta(m) > 1.e-20_r8)then
          ! from Abdul-Razzak and Ghan 2000
          g1=zeta(m)/eta(m)
          g1sqrt=sqrt(g1)
          g1=g1sqrt*g1
          g2=smc(m)/sqrt(eta(m)+3._r8*zeta(m))
          g2sqrt=sqrt(g2)
          g2=g2sqrt*g2
          sum=sum+(self%f1_(m)*g1+self%f2_(m)*g2)/(smc(m)*smc(m))
       else
          sum=large_maxsat
       endif
    enddo

    smax=1._r8/sqrt(sum)

  end function maxsat

  !------------------------------------------------------------------------------
  ! returns the ratio of SOA Hygroscopicity / Sulfate Hygroscopicity
  !------------------------------------------------------------------------------
  pure real(r8) function soa_equivso4_factor(self)
    class(aerosol_properties), intent(in) :: self

    soa_equivso4_factor = self%soa_equivso4_factor_

  end function soa_equivso4_factor

  !------------------------------------------------------------------------------
  ! returns the ratio of POM Hygroscopicity / Sulfate Hygroscopicity
  !------------------------------------------------------------------------------
  pure real(r8) function pom_equivso4_factor(self)
    class(aerosol_properties), intent(in) :: self

    pom_equivso4_factor = self%pom_equivso4_factor_

  end function pom_equivso4_factor

end module aerosol_properties_mod
