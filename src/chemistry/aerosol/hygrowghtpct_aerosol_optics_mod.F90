module hygrowghtpct_aerosol_optics_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_optics_mod, only: aerosol_optics
  use aerosol_state_mod, only: aerosol_state
  use aerosol_properties_mod, only: aerosol_properties
  use table_interp_mod, only: table_interp, table_interp_wghts, table_interp_calcwghts

  implicit none

  private
  public :: hygrowghtpct_aerosol_optics

  !> hygrowghtpct_aerosol_optics
  !! Table look up implementation of aerosol_optics to parameterize aerosol
  !! radiative properties in terms of weight precent of H2SO4/H2O solution
  type, extends(aerosol_optics) :: hygrowghtpct_aerosol_optics

     real(r8), allocatable :: totalmmr(:,:) ! total mmr of the aerosol
     real(r8), allocatable :: wgtpct(:,:)   ! weight precent of H2SO4/H2O solution

     real(r8), pointer :: sw_hygro_ext_wtp(:,:) ! short wave extinction table
     real(r8), pointer :: sw_hygro_ssa_wtp(:,:) ! short wave single-scatter albedo table
     real(r8), pointer :: sw_hygro_asm_wtp(:,:) ! short wave asymmetry table
     real(r8), pointer :: lw_hygro_abs_wtp(:,:) ! long wave absorption table

     real(r8), pointer :: tbl_wgtpct(:) ! weight precent dimenstion values

     integer :: nwtp ! weight precent dimenstion size

   contains

     procedure :: sw_props
     procedure :: lw_props

     final :: destructor

  end type hygrowghtpct_aerosol_optics

  interface hygrowghtpct_aerosol_optics
     procedure :: constructor
  end interface hygrowghtpct_aerosol_optics

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(aero_props, aero_state, ilist, ibin, ncol, nlev, wgtpct_in) result(newobj)

    class(aerosol_properties),intent(in) :: aero_props ! aerosol_properties object
    class(aerosol_state),intent(in) :: aero_state      ! aerosol_state object
    integer, intent(in) :: ilist  ! climate or a diagnostic list number
    integer, intent(in) :: ibin   ! bin number
    integer, intent(in) :: ncol   ! number of columns
    integer, intent(in) :: nlev   ! number of levels
    real(r8),intent(in) :: wgtpct_in(ncol,nlev) ! sulfate weight percent

    type(hygrowghtpct_aerosol_optics), pointer :: newobj

    integer :: ierr, nspec
    integer :: ispec
    integer :: i,k

    real(r8), pointer :: specmmr(:,:)   ! species mass mixing ratio

    allocate(newobj, stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%totalmmr(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%wgtpct(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    ! weight precent of H2SO4/H2O solution
    newobj%wgtpct(:ncol,:nlev) = wgtpct_in(:ncol,:nlev)

    call aero_props%optics_params(ilist, ibin, wgtpct=newobj%tbl_wgtpct, nwtp=newobj%nwtp)

    nspec = aero_props%nspecies(ilist, ibin)

    newobj%totalmmr(:,:) = 0._r8

    do ispec = 1,nspec

       call aero_state%get_ambient_mmr(ilist,ispec,ibin,specmmr)
       newobj%totalmmr(:ncol,:nlev) = newobj%totalmmr(:ncol,:nlev) + specmmr(:ncol,:nlev)

    end do

    call aero_props%optics_params(ilist, ibin,  &
         sw_hygro_ext_wtp=newobj%sw_hygro_ext_wtp, &
         sw_hygro_ssa_wtp=newobj%sw_hygro_ssa_wtp, &
         sw_hygro_asm_wtp=newobj%sw_hygro_asm_wtp, &
         lw_hygro_ext_wtp=newobj%lw_hygro_abs_wtp)

  end function constructor

  !------------------------------------------------------------------------------
  ! returns short wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)

    class(hygrowghtpct_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pext(ncol) ! parameterized specific extinction (m2/kg)
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)
    real(r8),intent(out) :: palb(ncol) ! parameterized asymmetry factor
    real(r8),intent(out) :: pasm(ncol) ! parameterized single scattering albedo

    integer :: icol
    type(table_interp_wghts) :: wghts(ncol)

    wghts = table_interp_calcwghts( self%nwtp, self%tbl_wgtpct, ncol, self%wgtpct(:ncol,ilev) )
    pext = table_interp( ncol, self%nwtp, wghts, self%sw_hygro_ext_wtp(:,iwav) )
    pabs = (1._r8 - table_interp( ncol, self%nwtp, wghts, self%sw_hygro_ssa_wtp(:,iwav)))*pext
    pasm = table_interp( ncol, self%nwtp, wghts, self%sw_hygro_asm_wtp(:,iwav) )

    do icol = 1, ncol

       pext(icol) = pext(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = pabs(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = max(0._r8,pabs(icol))
       pabs(icol) = min(pext(icol),pabs(icol))

       palb(icol) = 1._r8-pabs(icol)/max(pext(icol),1.e-40_r8)
       palb(icol) = 1._r8-pabs(icol)/max(pext(icol),1.e-40_r8)

    end do

  end subroutine sw_props

  !------------------------------------------------------------------------------
  ! returns long wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine lw_props(self, ncol, ilev, iwav, pabs)

    class(hygrowghtpct_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)

    integer :: icol
    type(table_interp_wghts) :: wghts(ncol)

    wghts = table_interp_calcwghts( self%nwtp, self%tbl_wgtpct, ncol, self%wgtpct(:ncol,ilev) )

    pabs = table_interp( ncol, self%nwtp, wghts, self%lw_hygro_abs_wtp(:,iwav) )

    do icol = 1, ncol

       pabs(icol) = pabs(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = max(0._r8,pabs(icol))

    end do

  end subroutine lw_props

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)

    type(hygrowghtpct_aerosol_optics), intent(inout) :: self

    deallocate(self%totalmmr)
    deallocate(self%wgtpct)

    nullify(self%tbl_wgtpct)
    nullify(self%sw_hygro_ext_wtp)
    nullify(self%sw_hygro_ssa_wtp)
    nullify(self%sw_hygro_asm_wtp)
    nullify(self%lw_hygro_abs_wtp)

  end subroutine destructor

end module hygrowghtpct_aerosol_optics_mod
