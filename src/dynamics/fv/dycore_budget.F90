module dycore_budget
use shr_kind_mod, only: r8=>shr_kind_r8
implicit none

public :: print_budget
real(r8), parameter :: eps      = 1.0E-9_r8
real(r8), parameter :: eps_mass = 1.0E-12_r8

!=========================================================================================
contains
!=========================================================================================

subroutine print_budget(hstwr)

  use spmd_utils,             only: masterproc
  use cam_abortutils,         only: endrun  
  use cam_logfile,            only: iulog
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use budgets,                only: budget_get_global, is_budget, thermo_budget_histfile_num, thermo_budget_history
  use cam_thermo,             only: thermo_budget_vars_descriptor, thermo_budget_num_vars, thermo_budget_vars_massv, &
                                    teidx, seidx, keidx, poidx
  use cam_thermo,             only: teidx, seidx, keidx, poidx
  use cam_thermo,             only: thermo_budget_vars_descriptor, thermo_budget_num_vars, thermo_budget_vars_massv

  ! arguments
  logical, intent(in) :: hstwr(:)

  ! Local variables
  character(len=*), parameter :: subname = 'check_energy:print_budgets'

  !--------------------------------------------------------------------------------------

  if (masterproc .and. thermo_budget_history .and. hstwr(thermo_budget_histfile_num)) then
     call endrun(subname//' is not implemented for the FV dycore')
  end if
end subroutine print_budget
!=========================================================================================
function abs_diff(a,b,pf)
  real(r8),                   intent(in) :: a,b
  character(LEN=5), optional, intent(out):: pf
  real(r8)                               :: abs_diff
  if (abs(b)>eps) then
    abs_diff = abs((b-a)/b)
  else
    abs_diff = abs(b-a)
  end if
  If (present(pf)) then
    if (abs_diff>eps) then
      pf = ' FAIL'
    else
      pf = ' PASS'
    end if
  end if
end function abs_diff
end module dycore_budget
