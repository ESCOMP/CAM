module cam_instance

  implicit none
  public

  integer          , public :: atm_id
  integer          , public :: inst_index
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix

!===============================================================================
CONTAINS
!===============================================================================

  subroutine cam_instance_init(atm_id_in, inst_name_in, inst_index_in, inst_suffix_in)

    integer          , intent(in) :: atm_id_in
    character(len=*) , intent(in) :: inst_name_in
    integer          , intent(in) :: inst_index_in
    character(len=*) , intent(in) :: inst_suffix_in

    ! The following sets the module variables
    atm_id      = atm_id_in
    inst_name   = inst_name_in
    inst_index  = inst_index_in
    inst_suffix = inst_suffix_in

  end subroutine cam_instance_init

end module cam_instance
