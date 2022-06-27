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
     integer :: nbins_ = 0
     integer :: ncnst_tot_ = 0
     integer, allocatable :: nmasses_(:)
     integer, allocatable :: nspecies_(:)
     integer, allocatable :: indexer_(:,:)
     real(r8), allocatable :: alogsig_(:)
     real(r8), allocatable :: f1_(:)
     real(r8), allocatable :: f2_(:)
   contains
     procedure :: initialize => aero_props_init
     procedure :: nbins
     procedure :: ncnst_tot
     procedure :: nspecies
     procedure,private :: n_masses_all_bins
     procedure,private :: n_masses_per_bin
     generic :: nmasses => n_masses_all_bins,n_masses_per_bin
     procedure :: indexer
     procedure :: maxsat
     procedure(aero_amcube), deferred :: amcube
     procedure :: alogsig
     procedure(aero_props_get), deferred :: get
     procedure(aero_actfracs), deferred :: actfracs
     procedure(aero_num_names), deferred :: num_names
     procedure(aero_mmr_names), deferred :: mmr_names
     procedure(aero_apply_num_limits), deferred :: apply_number_limits

     procedure :: final=>aero_props_final
  end type aerosol_properties

  interface

     !------------------------------------------------------------------------
     ! returns aerosol properties:
     !  density
     !  hygroscopicity
     !------------------------------------------------------------------------
     subroutine aero_props_get(self, bin_ndx, species_ndx, density,hygro)
       import :: aerosol_properties, r8
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx             ! bin index
       integer, intent(in) :: species_ndx         ! species index
       real(r8), optional, intent(out) :: density ! density (kg/m3)
       real(r8), optional, intent(out) :: hygro   ! hygroscopicity
     end subroutine aero_props_get

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
     ! returns constituents names of aersol number mixing ratios
     !------------------------------------------------------------------------
     subroutine aero_num_names(self, bin_ndx, name_a, name_c)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       character(len=32), intent(out) :: name_a ! constituent name of ambient aerosol number dens
       character(len=32), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens
     end subroutine aero_num_names

     !------------------------------------------------------------------------
     ! returns constituents names of aersol mass mixing ratios
     !------------------------------------------------------------------------
     subroutine aero_mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
       import :: aerosol_properties
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       integer, intent(in) :: species_ndx       ! species number
       character(len=32), intent(out) :: name_a ! constituent name of ambient aerosol MMR
       character(len=32), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR
     end subroutine aero_mmr_names

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

  end interface

contains

  !------------------------------------------------------------------------------
  ! object initializer
  !------------------------------------------------------------------------------
  subroutine aero_props_init(self, nbin, ncnst, nspec, nmasses, alogsig, f1,f2 )
    class(aerosol_properties), intent(inout) :: self
    integer, intent(in) :: nbin               ! number of bins
    integer, intent(in) :: ncnst              ! total number of constituents
    integer, intent(in) :: nspec(nbin)        ! number of species in each bin
    integer, intent(in) :: nmasses(nbin)      ! number of masses in each bin
    real(r8),intent(in) :: alogsig(nbin)      ! natural log of the standard deviation (sigma) of the aerosol bins
    real(r8),intent(in) :: f1(nbin)           ! abdul-razzak functions of width
    real(r8),intent(in) :: f2(nbin)           ! abdul-razzak functions of width

    integer :: l,m,mm

    allocate(self%nspecies_(nbin))
    allocate(self%nmasses_(nbin))
    allocate(self%alogsig_(nbin))
    allocate(self%f1_(nbin))
    allocate(self%f2_(nbin))

    allocate( self%indexer_(nbin,0:maxval(nmasses)) )

    ! Local indexing compresses the mode and number/mass indices into one index.
    ! This indexing is used by the pointer arrays used to reference state and pbuf
    ! fields. We add number = 0, total mass = 1 (if available), and mass from each
    ! constituency into mm.

    self%indexer_ = -1
    mm = 0

    do m=1,nbin
       do l = 0,nmasses(m)
          mm = mm+1
          self%indexer_(m,l) = mm
       end do
    end do

    self%nbins_ = nbin
    self%ncnst_tot_ = ncnst
    self%nmasses_(:) = nmasses(:)
    self%nspecies_(:) = nspec(:)
    self%alogsig_(:) = alogsig(:)
    self%f1_(:) = f1(:)
    self%f2_(:) = f2(:)

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
  pure integer function nspecies(self,bin_ndx)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    nspecies = self%nspecies_(bin_ndx)
  end function nspecies

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
  pure integer function nbins(self)
    class(aerosol_properties), intent(in) :: self

    nbins = self%nbins_
  end function nbins

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
  pure real(r8) function alogsig(self, bin_ndx)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    alogsig = self%alogsig_(bin_ndx)
  end function alogsig

  !------------------------------------------------------------------------------
  ! returns maximum supersaturation
  !------------------------------------------------------------------------------
  function maxsat(self, zeta,eta,smc) result(smax)

    !-------------------------------------------------------------------------
    ! Calculates maximum supersaturation for multiple competing aerosols.
    !
    ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
    !-------------------------------------------------------------------------

    class(aerosol_properties), intent(in) :: self
    real(r8), intent(in)  :: zeta(self%nbins_) ! Abdul-Razzak and Ghan eq 10
    real(r8), intent(in)  :: eta(self%nbins_)  ! Abdul-Razzak and Ghan eq 11
    real(r8), intent(in)  :: smc(self%nbins_)  ! critical supersaturation

    real(r8) :: smax ! maximum supersaturation

    integer  :: m
    integer  :: nbins
    real(r8) :: sum, g1, g2, g1sqrt, g2sqrt

    smax=0.0_r8
    nbins = self%nbins_

    check_loop: do m=1,nbins
       if((zeta(m) > 1.e5_r8*eta(m)) .or. (smc(m)*smc(m) > 1.e5_r8*eta(m))) then
          ! weak forcing -- essentially none activated
          smax=1.e-20_r8
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
          g1=zeta(m)/eta(m)
          g1sqrt=sqrt(g1)
          g1=g1sqrt*g1
          g2=smc(m)/sqrt(eta(m)+3._r8*zeta(m))
          g2sqrt=sqrt(g2)
          g2=g2sqrt*g2
          sum=sum+(self%f1_(m)*g1+self%f2_(m)*g2)/(smc(m)*smc(m))
       else
          sum=1.e20_r8
       endif
    enddo

    smax=1._r8/sqrt(sum)

  end function maxsat

end module aerosol_properties_mod
