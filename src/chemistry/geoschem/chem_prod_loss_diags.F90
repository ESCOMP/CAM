module chem_prod_loss_diags
  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : clscnt1, clscnt4, gas_pcnst, clsmap, permute
  use ppgrid, only : pver
  use chem_mods, only : rxntot
  use cam_history, only : addfld, outfld, add_default
  !use mo_tracname, only : solsym

  implicit none

  private
  public :: chem_prod_loss_diags_init
  public :: chem_prod_loss_diags_out

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine chem_prod_loss_diags_init

  end subroutine chem_prod_loss_diags_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine chem_prod_loss_diags_out( ncol, lchnk, base_sol, reaction_rates, prod_in, loss_in, xhnm )

    integer,  intent(in) :: ncol, lchnk
    real(r8), intent(in) :: base_sol(ncol,pver,gas_pcnst)
    real(r8), intent(in) :: reaction_rates(ncol,pver,max(1,rxntot))
    real(r8), intent(in) :: prod_in(ncol,pver,max(1,clscnt4))
    real(r8), intent(in) :: loss_in(ncol,pver,max(1,clscnt4))
    real(r8), intent(in) :: xhnm(ncol,pver)

  end subroutine chem_prod_loss_diags_out

end module chem_prod_loss_diags

