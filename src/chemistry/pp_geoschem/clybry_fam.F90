!-----------------------------------------------------------------------
!
! Manages the adjustment of ClOy and BrOy family components in response
! to conservation issues resulting from advection.
!
! Created by: Francis Vitt
! Date: 21 May 2008
! Modified by Stacy Walters
! Date: 13 August 2008
!-----------------------------------------------------------------------

module clybry_fam

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use ppgrid,        only : pcols, pver
  use chem_mods,     only : gas_pcnst, adv_mass
  use constituents,  only : pcnst
  use short_lived_species,only: set_short_lived_species,get_short_lived_species

  implicit none

  save

  private
  public :: clybry_fam_set
  public :: clybry_fam_adj
  public :: clybry_fam_init

  integer :: id_cly,id_bry

  integer :: id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2
  integer :: id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr

  logical :: has_clybry

contains

  !------------------------------------------
  !------------------------------------------
  subroutine clybry_fam_init

    !use mo_chem_utls, only : get_spc_ndx
    implicit none

    integer :: ids(16)

    !id_cly = get_spc_ndx('CLY')
    !id_bry = get_spc_ndx('BRY')

    !id_cl = get_spc_ndx('CL')
    !id_clo = get_spc_ndx('CLO')
    !id_hocl = get_spc_ndx('HOCL')
    !id_cl2 = get_spc_ndx('CL2')
    !id_cl2o2 = get_spc_ndx('CL2O2')
    !id_oclo = get_spc_ndx('OCLO')
    !id_hcl = get_spc_ndx('HCL')
    !id_clono2 = get_spc_ndx('CLONO2')

    !id_br = get_spc_ndx('BR')
    !id_bro = get_spc_ndx('BRO')
    !id_hbr = get_spc_ndx('HBR')
    !id_brono2 = get_spc_ndx('BRONO2')
    !id_brcl = get_spc_ndx('BRCL')
    !id_hobr = get_spc_ndx('HOBR')

    !ids = (/ id_cly,id_bry, &
    !         id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2, &
    !         id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr /)

    !has_clybry = all( ids(:) > 0 )

  endsubroutine clybry_fam_init

!--------------------------------------------------------------
! set the ClOy and BrOy mass mixing ratios
!  - this is call before advection
!--------------------------------------------------------------
  subroutine clybry_fam_set( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : get_nstep
    use physics_buffer, only : physics_buffer_desc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !real(r8) :: wrk(ncol,pver,2)
    !real(r8) :: mmr(pcols,pver,gas_pcnst)
    !integer  :: n, m

    if (.not. has_clybry) return

  end subroutine clybry_fam_set

!--------------------------------------------------------------
! adjust the ClOy and BrOy individual family members 
!  - this is call after advection
!--------------------------------------------------------------
  subroutine clybry_fam_adj( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : is_first_step
    use physics_buffer, only : physics_buffer_desc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

  end subroutine clybry_fam_adj

!--------------------------------------------------------------
! private methods
!--------------------------------------------------------------

!--------------------------------------------------------------
! compute the mass mixing retio of ClOy
!--------------------------------------------------------------
  function cloy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: cloy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    cloy = 0._r8
  
  end function cloy

!--------------------------------------------------------------
! compute the mass mixing retio of BrOy
!--------------------------------------------------------------
  function broy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: broy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    broy = 0._r8

  end function broy

end module clybry_fam
