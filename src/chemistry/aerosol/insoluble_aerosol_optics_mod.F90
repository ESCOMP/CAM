!-------------------------------------------------------------------------------
! Insoluble (non-hygroscopic) aerosol optical properties
!-------------------------------------------------------------------------------
module insoluble_aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  use aerosol_optics_mod, only: aerosol_optics
  use aerosol_properties_mod, only: aerosol_properties
  use aerosol_state_mod, only: aerosol_state

  implicit none

  private

  public :: insoluble_aerosol_optics

  type, extends(aerosol_optics) :: insoluble_aerosol_optics
     real(r8), pointer :: lw_abs(:)
     real(r8), pointer :: sw_ext(:)
     real(r8), pointer :: sw_ssa(:)
     real(r8), pointer :: sw_asm(:)

     ! aerosol mass mixing ratio
     real(r8), pointer :: mmr(:,:)

   contains

     procedure :: sw_props
     procedure :: lw_props

     final :: destructor

  end type insoluble_aerosol_optics

  interface insoluble_aerosol_optics
     procedure :: constructor
  end interface insoluble_aerosol_optics

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(aero_props, aero_state, ibin) result(newobj)

    class(aerosol_properties),intent(in) :: aero_props   ! aerosol_properties object
    class(aerosol_state),     intent(in) :: aero_state      ! aerosol_state object
    integer, intent(in) :: ibin   ! bin number

    type(insoluble_aerosol_optics), pointer :: newobj

    integer :: ierr

    allocate(newobj, stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    ! get mode properties
    call aero_props%optics_params(ibin, &
         sw_insoluble_ext=newobj%sw_ext, &
         sw_insoluble_ssa=newobj%sw_ssa, &
         sw_insoluble_asm=newobj%sw_asm, &
         lw_insoluble_ext=newobj%lw_abs )

    call aero_state%get_ambient_mmr(species_ndx=1, bin_ndx=ibin, mmr=newobj%mmr)

  end function constructor


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)

    type(insoluble_aerosol_optics), intent(inout) :: self

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! returns short wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)

    class(insoluble_aerosol_optics), intent(in) :: self

    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pext(ncol) ! parameterized specific extinction (m2/kg)
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)
    real(r8),intent(out) :: palb(ncol) ! parameterized single scattering albedo
    real(r8),intent(out) :: pasm(ncol) ! parameterized asymmetry factor

    pext(:ncol) = self%sw_ext(iwav) * self%mmr(:ncol,ilev)
    palb(:ncol) = self%sw_ssa(iwav)
    pasm(:ncol) = self%sw_asm(iwav)

    pabs(:ncol) = pext(:ncol) * ( 1._r8 - palb(:ncol) )

  end subroutine sw_props

  !------------------------------------------------------------------------------
  ! returns long wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine lw_props(self, ncol, ilev, iwav, pabs)

    class(insoluble_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)

    pabs(:ncol) = self%lw_abs(iwav) * self%mmr(:ncol,ilev)

  end subroutine lw_props

end module insoluble_aerosol_optics_mod
