module aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private
  public :: aerosol_optics

  !> aerosol_optics defines interfaces to optical properties of any aerosol package
  !!
  !! Each aerosol optics type must extend the abstract aerosol_optics class
  !! to define details of how aerosol optics properties are derived from
  !! aerosol states.
  type, abstract :: aerosol_optics

   contains

     procedure(aeropts_sw_props),deferred :: sw_props
     procedure(aeropts_lw_props),deferred :: lw_props

  end type aerosol_optics

  abstract interface

     !------------------------------------------------------------------------------
     ! returns short wave aerosol optics properties
     !------------------------------------------------------------------------------
     subroutine aeropts_sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)
       import :: aerosol_optics, r8

       class(aerosol_optics), intent(in) :: self
       integer, intent(in) :: ncol        ! number of columns
       integer, intent(in) :: ilev        ! vertical level index
       integer, intent(in) :: iwav        ! wave length index
       real(r8),intent(out) :: pext(ncol) ! parameterized specific extinction (m2/kg)
       real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)
       real(r8),intent(out) :: palb(ncol) ! parameterized asymmetry factor
       real(r8),intent(out) :: pasm(ncol) ! parameterized single scattering albedo

     end subroutine aeropts_sw_props

     !------------------------------------------------------------------------------
     ! returns long wave aerosol optics properties
     !------------------------------------------------------------------------------
     subroutine aeropts_lw_props(self, ncol, ilev, iwav, pabs)
       import :: aerosol_optics, r8

       class(aerosol_optics), intent(in) :: self
       integer, intent(in) :: ncol        ! number of columns
       integer, intent(in) :: ilev        ! vertical level index
       integer, intent(in) :: iwav        ! wave length index
       real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)

     end subroutine aeropts_lw_props

  end interface

end module aerosol_optics_mod
