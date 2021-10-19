module dyn_tests_utils
!----------------------------------------------------------------------- 
! 
! Utility data (and code) for dynamics testing
!
! The public items in this module are items used both by internal code
!   (e.g., analytic initial conditions) and by infrastructure which uses
!   the internal code (e.g., read_inidat). They cannot be members of the 
!   internal code because that is conditionally compiled.
!
!-----------------------------------------------------------------------


  implicit none
  private
  save

  integer, parameter :: vc_moist_pressure = 0 ! Moist pressure vertical coord
  integer, parameter :: vc_dry_pressure   = 1 ! Dry pressure vertical coord
  integer, parameter :: vc_height         = 2 ! Height vertical coord

  integer, parameter :: vc_str_lgth = 108     ! Length of string in 

  integer :: vc_dycore  !vertical coordinate of dynamical core - set in dyn_comp.F90
  integer :: vc_physics !vertical coordinate of physics - set in physconst.F90

  public :: vc_moist_pressure, vc_dry_pressure, vc_height, string_vc
  public :: vc_dycore, vc_physics, vc_str_lgth

contains
  subroutine string_vc(vc,str)
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog
    use string_utils,   only: int2str
    integer,                     intent(in)  :: vc    
    character (len=vc_str_lgth), intent(out) :: str
    
    select case (vc)
    case(vc_moist_pressure)
      str = 'Moist pressure/mass vertical coordinate'
    case(vc_dry_pressure)
      str = 'Dry pressure/mass vertical coordinate'
    case(vc_height)
      str = 'Height (z) vertical coordinate'
    case default
      write(iulog,*) 'string_vc: invalid vc= ',vc
      call endrun('string_vc: invalid vc ='//trim(int2str(vc)))
    end select
  end subroutine string_vc
end module dyn_tests_utils
