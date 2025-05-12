module dycore_budget
implicit none

public :: print_budget

!=========================================================================================
contains
!=========================================================================================

subroutine print_budget(hstwr)

  use spmd_utils,             only: masterproc
  use cam_abortutils,         only: endrun  
  use cam_budget,             only: thermo_budget_histfile_num, thermo_budget_history

  ! arguments
  logical, intent(in) :: hstwr(:)
  character(len=*), parameter :: subname = 'dycore_budget:print_budgets:'

  !--------------------------------------------------------------------------------------

  if (masterproc .and. thermo_budget_history .and. hstwr(thermo_budget_histfile_num)) then
     call endrun(subname//' is not implemented for the FV dycore')
  end if
end subroutine print_budget

end module dycore_budget
