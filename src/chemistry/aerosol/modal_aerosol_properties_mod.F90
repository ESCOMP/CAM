module modal_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props, rad_cnst_get_aer_props

  implicit none

  private

  public :: modal_aerosol_properties

  type, extends(aerosol_properties) :: modal_aerosol_properties
     private
     real(r8), allocatable :: sigmag_amode(:)
     real(r8), allocatable :: f1(:)
     real(r8), allocatable :: f2(:)
     real(r8), allocatable :: exp45logsig_(:)
     real(r8), allocatable :: voltonumblo_(:)
     real(r8), allocatable :: voltonumbhi_(:)
   contains
     procedure :: abdraz_f1
     procedure :: abdraz_f2
     procedure :: get
     procedure :: exp45logsig
     procedure :: voltonumblo
     procedure :: voltonumbhi
     procedure :: amcube
     procedure :: actfracs
     procedure :: get_num_names
     procedure :: get_mmr_names
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
    integer,allocatable :: nmasses(:)
    real(r8),allocatable :: amcubecoefs(:)
    real(r8),allocatable :: alogsig(:)

    allocate(newobj)

    call rad_cnst_get_info(0, nmodes=nmodes)

    allocate(nspecies(nmodes))
    allocate(nmasses(nmodes))
    allocate(amcubecoefs(nmodes))
    allocate(alogsig(nmodes))

    allocate(newobj%sigmag_amode(nmodes))
    allocate(newobj%f1(nmodes))
    allocate(newobj%f2(nmodes))
    allocate(newobj%exp45logsig_(nmodes))
    allocate(newobj%voltonumblo_(nmodes))
    allocate(newobj%voltonumbhi_(nmodes))

    ncnst_tot = 0

    do m = 1, nmodes
       call rad_cnst_get_info(0, m, nspec=nspecies(m))

       ncnst_tot =  ncnst_tot + nspecies(m) + 1
       nmasses(m) = nspecies(m)

       call rad_cnst_get_mode_props(0, m, sigmag=newobj%sigmag_amode(m), &
                                    dgnumhi=dgnumhi, dgnumlo=dgnumlo )

       alogsig(m) = log(newobj%sigmag_amode(m))

       newobj%exp45logsig_(m) = exp(4.5_r8*alogsig(m)*alogsig(m))

       amcubecoefs(m)=3._r8/(4._r8*pi*newobj%exp45logsig_(m))

       newobj%f1(m) = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
       newobj%f2(m) = 1._r8 + 0.25_r8*alogsig(m)

       newobj%voltonumblo_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumlo**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
       newobj%voltonumbhi_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumhi**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )

    end do

    call newobj%initialize(nmodes,ncnst_tot,nspecies,nmasses,amcubecoefs,alogsig)
    deallocate(nspecies)
    deallocate(nmasses)
    deallocate(amcubecoefs)
    deallocate(alogsig)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_properties), intent(inout) :: self

    deallocate(self%sigmag_amode)
    deallocate(self%f1)
    deallocate(self%f2)
    deallocate(self%exp45logsig_)
    deallocate(self%voltonumblo_)
    deallocate(self%voltonumbhi_)

    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function exp45logsig(self,m)
    class(modal_aerosol_properties), intent(in) :: self
    integer,intent(in) :: m
    exp45logsig = self%exp45logsig_(m)
  end function exp45logsig

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function voltonumblo(self,m)
    class(modal_aerosol_properties), intent(in) :: self
    integer,intent(in) :: m
    voltonumblo = self%voltonumblo_(m)
  end function voltonumblo

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function voltonumbhi(self,m)
    class(modal_aerosol_properties), intent(in) :: self
    integer,intent(in) :: m
    voltonumbhi = self%voltonumbhi_(m)
  end function voltonumbhi

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get(self, m,l, density,hygro)

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m,l
    real(r8), optional, intent(out) :: density
    real(r8), optional, intent(out) :: hygro

    call rad_cnst_get_aer_props(0, m, l, density_aer=density, hygro_aer=hygro)

  end subroutine get

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f1(self,m) result(f)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = self%f1(m)
  end function abdraz_f1

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f2(self,m) result(f)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = self%f2(m)
  end function abdraz_f2

  !------------------------------------------------------------------------------
  ! amcube is overridden to keep MAM b4b
  !------------------------------------------------------------------------------
  pure real(r8) function amcube(self, m, volconc, numconc)

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8), intent(in) :: volconc
    real(r8), intent(in) :: numconc

    amcube = (3._r8*volconc/(4._r8*pi*self%exp45logsig_(m)*numconc))

  end function amcube

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine actfracs(self, m, smc, smax, fn, fm )
    use shr_spfn_mod, only: erf => shr_spfn_erf
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8),intent(in) :: smc
    real(r8),intent(in) :: smax
    real(r8),intent(out) :: fn, fm

    real(r8) :: x,y
    real(r8), parameter :: twothird = 2._r8/3._r8
    real(r8), parameter :: sq2      = sqrt(2._r8)

    x=twothird*(log(smc)-log(smax))/(sq2*self%alogsig(m))
    y=x-1.5_r8*sq2*self%alogsig(m)

    fn = 0.5_r8*(1._r8-erf(x))
    fm = 0.5_r8*(1._r8-erf(y))

  end subroutine actfracs

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine get_num_names(self, m, name_a, name_c)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    character(len=32), intent(out) :: name_a, name_c

    call rad_cnst_get_info(0, m, num_name=name_a, num_name_cw=name_c)
  end subroutine get_num_names

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine get_mmr_names(self, m,l, name_a, name_c)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m,l
    character(len=32), intent(out) :: name_a, name_c

    call rad_cnst_get_info(0, m, l, spec_name=name_a, spec_name_cw=name_c)
  end subroutine get_mmr_names

end module modal_aerosol_properties_mod
