module insoluble_aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  use aerosol_optics_mod, only: aerosol_optics
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: insoluble_aerosol_optics

  type, extends(aerosol_optics) :: insoluble_aerosol_optics
     real(r8), pointer :: lw_abs(:)
     real(r8), pointer :: sw_ext(:)
     real(r8), pointer :: sw_ssa(:)
     real(r8), pointer :: sw_asm(:)
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
  function constructor(aero_props, ilist, ibin) result(newobj)

    class(aerosol_properties),intent(in), target :: aero_props   ! aerosol_properties object
    integer, intent(in) :: ilist  ! climate or a diagnostic list number
    integer, intent(in) :: ibin   ! bin number

    type(insoluble_aerosol_optics), pointer :: newobj

    integer :: ierr

    allocate(newobj, stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    ! get mode properties
    call aero_props%optics_params(ilist, ibin, &
         sw_nonhygro_ext=newobj%sw_ext, &
         sw_nonhygro_ssa=newobj%sw_ssa, &
         sw_nonhygro_asm=newobj%sw_asm, &
         lw_nonhygro_ext=newobj%lw_abs )

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

    pext(:ncol) = self%sw_ext(iwav)
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

    pabs(:ncol) = self%lw_abs(iwav)

  end subroutine lw_props

end module insoluble_aerosol_optics_mod
