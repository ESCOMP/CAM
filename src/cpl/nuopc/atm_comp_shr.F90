module atm_comp_shr

  ! Model mesh info is here in order to be leveraged by CDEPS in line calls

  use ESMF
  use shr_kind_mod, only : r8 => shr_kind_r8, cl=>shr_kind_cl

  implicit none
  public

  type(ESMF_Clock)  :: model_clock    ! model clock
  type(ESMF_Mesh)   :: model_mesh     ! model_mesh
  character(len=cl) :: model_meshfile ! model mesh file 

end module atm_comp_shr
