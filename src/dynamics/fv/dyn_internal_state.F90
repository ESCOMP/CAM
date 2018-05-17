module dyn_internal_state

! Container for the FV dynamics internal state.

use dynamics_vars, only : t_fvdycore_state, t_fvdycore_grid, &
                          t_fvdycore_constants, t_fvdycore_vars

implicit none
private
save

public :: &
   get_dyn_state, &
   get_dyn_state_grid, &
   get_dyn_state_vars, &
   get_dyn_state_constants

type (t_fvdycore_state), target :: dyn_state 

!========================================================================================
contains
!========================================================================================

function get_dyn_state() result(dynstate)
   type(t_fvdycore_state), pointer :: dynstate
   dynstate => dyn_state
end function get_dyn_state

!========================================================================================

function get_dyn_state_grid() result(grid)
   type(t_fvdycore_grid), pointer :: grid
   grid => dyn_state%grid
end function get_dyn_state_grid

!========================================================================================

function get_dyn_state_vars() result(vars)
   type(t_fvdycore_vars), pointer :: vars
   vars => dyn_state%vars
end function get_dyn_state_vars

!========================================================================================

function get_dyn_state_constants() result(constants)
   type(t_fvdycore_constants), pointer :: constants
   constants => dyn_state%constants
end function get_dyn_state_constants

!========================================================================================

end module dyn_internal_state
