module aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: aerosol_properties

  type, abstract :: aerosol_properties
     private
     integer :: nbins_ = 0
     integer :: ncnst_tot_ = 0
     integer, allocatable :: nmasses_(:)
     integer, allocatable :: nspecies_(:)
     integer, allocatable :: indexer_(:,:)
     real(r8), allocatable :: amcubecoefs_(:)
     real(r8), allocatable :: alogsig_(:)
     real(r8), allocatable :: f1_(:)
     real(r8), allocatable :: f2_(:)
   contains
     procedure :: initialize => aero_props_init
     procedure :: nbins
     procedure :: ncnst_tot
     procedure :: nspecies
     procedure,private :: nmassesa
     procedure,private :: nmassesv
     generic :: nmasses => nmassesa,nmassesv
     procedure :: indexer
     procedure :: maxsat
     procedure :: amcubecoef
     procedure :: amcube
     procedure :: alogsig
     procedure(aero_props_get), deferred :: get
     procedure(aero_actfracs), deferred :: actfracs
     procedure(aero_num_names), deferred :: get_num_names
     procedure(aero_mmr_names), deferred :: get_mmr_names

     procedure :: final=>aero_props_final
  end type aerosol_properties

  interface

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_props_get(self, m,l, density,hygro)
       import
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: m,l
       real(r8), optional, intent(out) :: density
       real(r8), optional, intent(out) :: hygro
     end subroutine aero_props_get

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_actfracs(self, m, smc, smax, fn, fm )
       import
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: m
       real(r8),intent(in) :: smc
       real(r8),intent(in) :: smax
       real(r8),intent(out) :: fn, fm

     end subroutine aero_actfracs

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_num_names(self, m, name_a, name_c)
       import
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: m
       character(len=32), intent(out) :: name_a, name_c
     end subroutine aero_num_names

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_mmr_names(self, m,l, name_a, name_c)
       import
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: m,l
       character(len=32), intent(out) :: name_a, name_c
     end subroutine aero_mmr_names

  end interface

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_props_init(self, nbin, ncnst, nspec, nmasses, amcubecoefs, alogsig, f1,f2 )
    class(aerosol_properties), intent(inout) :: self
    integer, intent(in) :: nbin
    integer, intent(in) :: ncnst
    integer, intent(in) :: nspec(nbin)
    integer, intent(in) :: nmasses(nbin)
    real(r8),intent(in) :: amcubecoefs(nbin)
    real(r8),intent(in) :: alogsig(nbin)
    real(r8),intent(in) :: f1(nbin)
    real(r8),intent(in) :: f2(nbin)

    integer :: l,m,mm

    allocate(self%nspecies_(nbin))
    allocate(self%nmasses_(nbin))
    allocate(self%amcubecoefs_(nbin))
    allocate(self%alogsig_(nbin))
    allocate(self%f1_(nbin))
    allocate(self%f2_(nbin))

    allocate( self%indexer_(nbin,0:maxval(nmasses)) )

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
    self%amcubecoefs_(:) = amcubecoefs(:)
    self%alogsig_(:) = alogsig(:)
    self%f1_(:) = f1(:)
    self%f2_(:) = f2(:)

  end subroutine aero_props_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_props_final(self)
    class(aerosol_properties), intent(inout) :: self

    if (allocated(self%nspecies_)) then
       deallocate(self%nspecies_)
    end if
    if (allocated(self%nmasses_)) then
       deallocate(self%nmasses_)
    end if
    if (allocated(self%amcubecoefs_)) then
       deallocate(self%amcubecoefs_)
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
  !------------------------------------------------------------------------------
  pure integer function nspecies(self,m)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: m

    nspecies = self%nspecies_(m)
  end function nspecies

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure function nmassesv(self,m) result(val)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    integer :: val

    val = self%nmasses_(m)
  end function nmassesv

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure function nmassesa(self) result(arr)
    class(aerosol_properties), intent(in) :: self
    integer :: arr(self%nbins_)

    arr(:) = self%nmasses_(:)
  end function nmassesa

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure integer function indexer(self,m,l)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: m,l

    indexer = self%indexer_(m,l)
  end function indexer

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure integer function nbins(self)
    class(aerosol_properties), intent(in) :: self

    nbins = self%nbins_
  end function nbins

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure integer function ncnst_tot(self)
    class(aerosol_properties), intent(in) :: self

    ncnst_tot = self%ncnst_tot_
  end function ncnst_tot

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function amcubecoef(self, m)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: m

    amcubecoef = self%amcubecoefs_(m)
  end function amcubecoef

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function alogsig(self, m)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: m

    alogsig = self%alogsig_(m)
  end function alogsig

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function amcube(self, m, volconc, numconc)

    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8), intent(in) :: volconc
    real(r8), intent(in) :: numconc

    amcube = self%amcubecoefs_(m)*volconc/numconc

  end function amcube

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function maxsat(self, zeta,eta,smc) result(smax)

    !-------------------------------------------------------------------------
    ! Calculates maximum supersaturation for multiple competing aerosol modes.
    !
    ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
    !-------------------------------------------------------------------------

    class(aerosol_properties), intent(in) :: self
    real(r8), intent(in)  :: zeta(self%nbins_)
    real(r8), intent(in)  :: eta(self%nbins_)
    real(r8), intent(in)  :: smc(self%nbins_) ! critical supersaturation for number mode radius

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
